#!/usr/bin/env nextflow

/* This workflow is a adapted from two previous pipeline written in Snakemake:
- https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms
- https://github.com/nanoporetech/pipeline-nanopore-denovo-isoforms
*/

import groovy.json.JsonBuilder;
import nextflow.util.BlankSeparatedList;
import java.util.ArrayList;
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'
include { reference_assembly } from './subworkflows/reference_assembly'
include { denovo_assembly } from './subworkflows/denovo_assembly'
include { gene_fusions } from './subworkflows/JAFFAL/gene_fusions'
include { differential_expression } from './subworkflows/differential_expression'



process summariseConcatReads {
    // concatenate fastq and fastq.gz in a dir write stats

    label "isoforms"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        tuple val(meta.sample_id), path("${meta.sample_id}.fastq"), emit: input_reads
        tuple val(meta.sample_id), path('*.stats'), emit: summary
    script:
    """

    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} >  ${meta.sample_id}.fastq
    """
}

process getVersions {
    label "isoforms"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import aplanat; print(f'aplanat,{aplanat.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    python -c "import sklearn; print(f'scikit-learn,{sklearn.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    python -c "import pychopper; print(f'pychopper,{pychopper.__version__}')" >> versions.txt
    gffread --version | sed 's/^/gffread,/' >> versions.txt
    seqkit version | head -n 1 | sed 's/ /,/' >> versions.txt
    stringtie --version | sed 's/^/stringtie,/' >> versions.txt
    gffcompare --version | head -n 1 | sed 's/ /,/' >> versions.txt
    spoa --version | sed 's/^/spoa,/' >> versions.txt
#     isONclust2 version | sed 's/ version: /,/' >> versions.txt
    """
}


process getParams {
    label "isoforms"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
        println('test')
        println(params.workDir)
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process preprocess_reads {
    /*
    Concatenate reads from a sample directory.
    Optionally classify, trim, and orient cDNA reads using pychopper
    */

    label "isoforms"
    cpus 4

    input:
        tuple val(sample_id), path(input_reads)
    output:
         tuple val(sample_id), path("${sample_id}_full_length_reads.fastq"), emit: full_len_reads
         path '*.tsv',  emit: report
    script:
        """
        pychopper -t ${params.threads} ${params.pychopper_opts} ${input_reads} ${sample_id}_full_length_reads.fastq
        mv pychopper.tsv ${sample_id}_pychopper.tsv
        workflow-glue generate_pychopper_stats --data ${sample_id}_pychopper.tsv --output .

        # Add sample id column
        sed "1s/\$/\tsample_id/; 1 ! s/\$/\t${sample_id}/" ${sample_id}_pychopper.tsv > tmp
        mv tmp ${sample_id}_pychopper.tsv
        """
}

process build_minimap_index{
    /*
    Build minimap index from reference genome
    */
    label "isoforms"
    cpus params.threads

    input:
        path reference
    output:
        path "genome_index.mmi", emit: index
    script:
    """
    minimap2 -t ${params.threads} ${params.minimap_index_opts} -I 1000G -d "genome_index.mmi" ${reference}
    """
}

process split_bam{
    /*
    Partition BAM file into loci or bundles with `params.bundle_min_reads` minimum size
    If no splitting required, just create single symbolic link to a single bundle.

    Output tuples containing `sample_id` so bundles can be combined later in th pipeline.
    */

    label 'isoforms'
    cpus params.threads

    input:
        tuple val(sample_id), path(bam)
    output:
        tuple val(sample_id), path('*.bam'), emit: bundles
    script:
    """
    n=`samtools view -c $bam`
    if [[ n -lt 1 ]]
    then
        echo 'There are no reads mapping for $sample_id. Exiting!'
        exit 1
    fi

    re='^[0-9]+\$'

    if [[ $params.bundle_min_reads =~ \$re ]]
    then
        echo "Bundling up the bams"
        seqkit bam -j ${params.threads} -N ${params.bundle_min_reads} ${bam} -o  bam_bundles/
        let i=1
        for b in bam_bundles/*.bam; do
            echo \$b
            newname="${sample_id}_batch_\${i}.bam"
            mv \$b \$newname
           ((i++))
        done
    else
        echo 'no bundling'
        ln -s ${bam} ${sample_id}_batch_1.bam
    fi
    """
}


process assemble_transcripts{
    /*
    Assemble transcripts using stringtie.
    Take aligned reads in bam format that may be a chunk of a larger alignment file.
    Optionally use reference annotation to guide assembly.

    Output gff annotation files in a tuple with `sample_id` for combining into samples later in the pipeline.
    */
    label 'isoforms'
    cpus params.threads

    input:
        tuple val(sample_id), path(bam)
        path ref_annotation
    output:
        tuple val(sample_id), path('*.gff'), emit: gff_bundles
    script:
        def G_FLAG = ref_annotation.name.startsWith('OPTIONAL_FILE') ? '' : "-G ${ref_annotation}"
        def prefix =  bam.name.split(/\./)[0]

    """
    stringtie --rf ${G_FLAG} -L -v -p ${task.cpus} ${params.stringtie_opts} \
    -o  ${prefix}.gff -l ${prefix} ${bam} 2>/dev/null
     """
}


process merge_gff_bundles{
    /*
    Merge gff bundles into a single gff file per sample.
    */
    label 'isoforms'

    input:
        tuple val(sample_id), path (gff_bundle)
    output:
        tuple val(sample_id), path('*.gff'), emit: gff
    script:
    def merged_gff = "transcripts_${sample_id}.gff"
    """
    echo '##gff-version 2' >> $merged_gff;
    echo '#pipeline-nanopore-isoforms: stringtie' >> $merged_gff;

    for fn in ${gff_bundle};
    do
        grep -v '#' \$fn >> $merged_gff

    done
    """
}

process run_gffcompare{
    /*
    Compare query and reference annotations.
    If ref_annotation is an optional file, just make an empty directory to satisfy
    the requirements of the downstream processes.
    */

    label 'isoforms'

    input:
       tuple val(sample_id), path(query_annotation)
       path ref_annotation
    output:
        tuple val(sample_id), path("${sample_id}_gffcompare"), emit: gffcmp_dir
        path ("${sample_id}_annotated.gtf"), emit: gtf, optional: true
    script:
    def out_dir = "${sample_id}_gffcompare"

    if ( ref_annotation.name.startsWith('OPTIONAL_FILE') ){
        """
        mkdir $out_dir
        """
    } else {
        """
        mkdir $out_dir
        echo "Doing comparison of reference annotation: ${ref_annotation} and the query annotation"

        gffcompare -o ${out_dir}/str_merged -r ${ref_annotation} \
            ${params.gffcompare_opts} ${query_annotation}

        workflow-glue generate_tracking_summary --tracking $out_dir/str_merged.tracking \
            --output_dir ${out_dir} --annotation ${ref_annotation}

        mv *.tmap $out_dir
        mv *.refmap $out_dir
        cp ${out_dir}/str_merged.annotated.gtf ${sample_id}_annotated.gtf
        """
    }
}


process get_transcriptome{
        /*
        Write out a transcriptome file based on the query gff annotations.
        */
        label 'isoforms'

        input:
            tuple val(sample_id), path(transcripts_gff), path(gffcmp_dir), path(reference_seq)
        output:
            tuple val(sample_id), path("*.fas"), emit: transcriptome

        script:
        def transcriptome = "${sample_id}_transcriptome.fas"
        def merged_transcriptome = "${sample_id}_merged_transcriptome.fas"
        """
        gffread -g ${reference_seq} -w ${transcriptome} ${transcripts_gff}
        if  [ "\$(ls -A $gffcmp_dir)" ];
            then
                gffread -F -g ${reference_seq} -w ${merged_transcriptome} $gffcmp_dir/str_merged.annotated.gtf
        fi
        """
}

process merge_transcriptomes {
    // Merge the transcriptomes from all samples
    label 'isoforms'
    input:
        path "query_annotations/*"
        path ref_annotation
        path ref_genome
    output:
        path "non_redundant.fasta", emit: fasta
        path "stringtie.gtf", emit: gtf
    """
    stringtie --merge -G $ref_annotation -p ${task.cpus} -o stringtie.gtf query_annotations/*
    seqkit subseq --feature "transcript" --gtf-tag "transcript_id" --gtf stringtie.gtf $ref_genome > temp_transcriptome.fasta
    seqkit rmdup -s < temp_transcriptome.fasta > temp_del_repeats.fasta
    cat temp_del_repeats.fasta | sed 's/>.* />/'  | sed -e 's/_[0-9]* \\[/ \\[/' > temp_rm_empty_seq.fasta
    awk 'BEGIN {RS = ">" ; FS = "\\n" ; ORS = ""} \$2 {print ">"\$0}' temp_rm_empty_seq.fasta > non_redundant.fasta
    rm temp_transcriptome.fasta
    rm temp_del_repeats.fasta
    rm temp_rm_empty_seq.fasta
    """
}


process makeReport {

    label "isoforms"

    input:
        path versions
        path "params.json"
        path "pychopper_report/*"
        path"jaffal_csv/*"
        val sample_ids
        path seq_summaries
        path "aln_stats/*"
        path gffcmp_dir
        path "gff_annotation/*"
        path "de_report/*"
        path "seqkit/*"
    output:
        path("wf-transcriptomes-*.html"), emit: report
    script:
        // Convert the sample_id arrayList.
        sids = new BlankSeparatedList(sample_ids)
        def report_name = "wf-transcriptomes-report.html"
        def OPT_DENOVO = params.transcriptome_source == "denovo" ? "--denovo" : ''
    """
    if [ -f "de_report/OPTIONAL_FILE" ]; then
        dereport=""
    else
        dereport="--de_report true --de_stats "seqkit/*""
        mv de_report/*.gtf de_report/stringtie_merged.gtf
    fi
    if [ -f "gff_annotation/OPTIONAL_FILE" ]; then
        OPT_GFF=""
    else
        OPT_GFF="--gffcompare_dir ${gffcmp_dir} --gff_annotation gff_annotation/*"
        
    fi
    if [ -f "jaffal_csv/OPTIONAL_FILE" ]; then
        OPT_JAFFAL_CSV=""
    else
        OPT_JAFFAL_CSV="--jaffal_csv jaffal_csv/*"
    fi
    if [ -f "aln_stats/OPTIONAL_FILE" ]; then
        OPT_ALN=""
    else
        OPT_ALN="--alignment_stats aln_stats/*"
    fi
    if [ -f "pychopper_report/OPTIONAL_FILE" ]; then
        OPT_PC_REPORT=""
    else
        OPT_PC_REPORT="--pychop_report pychopper_report/*"
    fi
    workflow-glue report --report $report_name \
    --versions $versions \
    --params params.json \
    \$OPT_ALN \
    \$OPT_PC_REPORT \
    --sample_ids $sids \
    --summaries $seq_summaries \
    \$OPT_GFF \
    --isoform_table_nrows $params.isoform_table_nrows \
    \$OPT_JAFFAL_CSV \
    $OPT_DENOVO \
    \$dereport 

    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    label "isoforms"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}

//meh
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters (missing protocol or profile will exit the run.)
if (params.input) { 
    ch_input = file(params.input) 
} else {
    exit 1, 'Input samplesheet not specified!'
}

// Function to check if running offline
def isOffline() {
    try {
        return NXF_OFFLINE as Boolean
    }
    catch( Exception e ) {
        return false
    }
}

def ch_guppy_model  = Channel.empty()
def ch_guppy_config = Channel.empty()

if (!params.skip_basecalling) {
    if (!params.guppy_config) {
        if (!params.flowcell) { exit 1, "Please specify a valid flowcell identifier for basecalling!" }
        if (!params.kit)      { exit 1, "Please specify a valid kit identifier for basecalling!"      }
    } else if (file(params.guppy_config).exists()) {
        ch_guppy_config = Channel.fromPath(params.guppy_config)
    }

    if (params.guppy_model) {
        if (file(params.guppy_model).exists()) {
            ch_guppy_model = Channel.fromPath(params.guppy_model)
        }
    }
} else {
    if (!params.skip_demultiplexing) {
        if (!params.barcode_kit) {
            params.barcode_kit = 'Auto'
        }

        def qcatBarcodeKitList = ['Auto', 'RBK001', 'RBK004', 'NBD103/NBD104',
                                'NBD114', 'NBD104/NBD114', 'PBC001', 'PBC096',
                                'RPB004/RLB001', 'PBK004/LWB001', 'RAB204', 'VMK001', 'DUAL']

        if (params.barcode_kit && qcatBarcodeKitList.contains(params.barcode_kit)) {
            if (params.input_path) {
                ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true)
            } else {
                exit 1, "Please specify a valid input fastq file to perform demultiplexing!"
            }
        } else {
            exit 1, "Please provide a barcode kit to demultiplex with qcat. Valid options: ${qcatBarcodeKitList}"
        }
    }
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$baseDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GET_NANOLYSE_FASTA    } from './modules/local/get_nanolyse_fasta'
include { GUPPY                 } from './modules/local/guppy'
include { DEMUX_FAST5           } from './modules/local/demux_fast5'
include { QCAT                  } from './modules/local/qcat'
include { MULTIQC               } from './modules/local/multiqc'

/*
* SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
*/

include { INPUT_CHECK                      } from './subworkflows/local/input_check'
include { PREPARE_GENOME                   } from './subworkflows/local/prepare_genome'
include { QCBASECALL_PYCOQC_NANOPLOT       } from './subworkflows/local/qcbasecall_pycoqc_nanoplot'
include { QCFASTQ_NANOPLOT_FASTQC          } from './subworkflows/local/qcfastq_nanoplot_fastqc'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
* MODULE: Installed directly from nf-core/modules
*/
include { NANOLYSE                    } from './modules/nf-core/modules/nanolyse/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
* SUBWORKFLOW: Consisting entirely of nf-core/modules
*/

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []
workflow nanoseq {
    if (params.input_path) {
        ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true)
    } else {
        ch_input_path = 'not_changed'
    }
    

    /*
    * Create empty software versions channel to mix
    */
    ch_software_versions = Channel.empty()

    /*
    * SUBWORKFLOW: Read in samplesheet, validate and stage input files
    */
    INPUT_CHECK ( ch_input, ch_input_path )
        .set { ch_sample }

    if (!params.skip_basecalling) {
        ch_sample
            .first()
            .map { it[0] }
            .set { ch_sample_name }

        /*
        * MODULE: Basecalling and demultipexing using Guppy
        */
        GUPPY ( ch_input_path, ch_sample_name, ch_guppy_config.ifEmpty([]), ch_guppy_model.ifEmpty([]) )
        ch_guppy_summary = GUPPY.out.summary
        ch_software_versions = ch_software_versions.mix(GUPPY.out.versions.ifEmpty(null))

        if (params.skip_demultiplexing) {
            ch_sample
                .map { it -> [ it[0], it[0].id, it[2], it[3], it[4], it[5] ] }
                .set { ch_sample }
        }

        GUPPY.out.fastq
            .flatten()
            .map { it -> [ it, it.baseName.substring(0,it.baseName.lastIndexOf('.')) ] }
            .join(ch_sample, by: 1) // join on barcode
            .map { it -> [ it[2], it[1], it[3], it[4], it[5], it[6] ] }
            .set { ch_fastq }
        if (params.output_demultiplex_fast5) {

            /*
            * MODULE: Demultiplex fast5 files using ont_fast5_api/demux_fast5
            */
            DEMUX_FAST5 ( ch_input_path, ch_guppy_summary )
            ch_software_versions = ch_software_versions.mix(DEMUX_FAST5.out.versions.ifEmpty(null))
        }
    } else {
        ch_guppy_summary = Channel.empty()

        if (!params.skip_demultiplexing) {

            /*
            * MODULE: Demultipexing using qcat
            */
            QCAT ( ch_input_path )
            ch_fastq = Channel.empty()
            QCAT.out.fastq
                .flatten()
                .map { it -> [ it, it.baseName.substring(0,it.baseName.lastIndexOf('.'))] }
                .join(ch_sample, by: 1) // join on barcode
                .map { it -> [ it[2], it[1], it[3], it[4], it[5], it[6] ] }
                .set { ch_fastq }
            ch_software_versions = ch_software_versions.mix(QCAT.out.versions.ifEmpty(null))
        } else {
            ch_fastq = Channel.empty()
        }
    }

    if (params.run_nanolyse) {
        ch_fastq
            .map { it -> [ it[0], it[1] ] }
            .set { ch_fastq_nanolyse }

        if (!params.nanolyse_fasta) {
            if (!isOffline()) {
                GET_NANOLYSE_FASTA ().set { ch_nanolyse_fasta }
            } else {
                exit 1, "NXF_OFFLINE=true or -offline has been set so cannot download lambda.fasta.gz file for running NanoLyse! Please explicitly specify --nanolyse_fasta."
            }
        } else {
            ch_nanolyse_fasta = file(params.nanolyse_fasta, checkIfExists: true)
        }
        /*
        * MODULE: DNA contaminant removal using NanoLyse
        */
        NANOLYSE ( ch_fastq_nanolyse, ch_nanolyse_fasta )
        NANOLYSE.out.fastq
            .join( ch_sample )
            .map { it -> [ it[0], it[1], it[3], it[4], it[5], it[6] ]}
            .set { ch_fastq }
        ch_software_versions = ch_software_versions.mix(NANOLYSE.out.versions.first().ifEmpty(null))
    }

    ch_pycoqc_multiqc = Channel.empty()
    ch_fastqc_multiqc = Channel.empty()
    if (!params.skip_qc) {
        if (!params.skip_basecalling) {

            /*
            * SUBWORKFLOW: Basecalling QC with PycoQC and Nanoplot
            */
            QCBASECALL_PYCOQC_NANOPLOT ( ch_guppy_summary , params.skip_pycoqc, params.skip_nanoplot )
            ch_software_versions = ch_software_versions.mix(QCBASECALL_PYCOQC_NANOPLOT.out.pycoqc_version.first().ifEmpty(null))
            ch_pycoqc_multiqc    = QCBASECALL_PYCOQC_NANOPLOT.out.pycoqc_multiqc.ifEmpty([])
        }

        /*
        * SUBWORKFLOW: Fastq QC with Nanoplot and fastqc
        */
        QCFASTQ_NANOPLOT_FASTQC ( ch_fastq, params.skip_nanoplot, params.skip_fastqc)
        if (params.skip_basecalling) {
            ch_software_versions = ch_software_versions.mix(QCFASTQ_NANOPLOT_FASTQC.out.nanoplot_version.first().ifEmpty(null))
        }
        ch_software_versions = ch_software_versions.mix(QCFASTQ_NANOPLOT_FASTQC.out.fastqc_version.first().ifEmpty(null))
        ch_fastqc_multiqc    = QCFASTQ_NANOPLOT_FASTQC.out.fastqc_multiqc.ifEmpty([])
    }

    print('~~~~~~~~~~~~~~~~~meh~~~~~~~~~~~~~~~~~~')
    ch_fastq.view()
    
}

// workflow module
workflow pipeline {
    take:
        reads
        ref_genome
        ref_annotation
        jaffal_refBase
        jaffal_genome
        jaffal_annotation
        condition_sheet
        ref_transcriptome
    main:
        map_sample_ids_cls = {it ->
        /* Harmonize tuples
        output:
            tuple val(sample_id), path('*.gff')
        When there are multiple paths, will emit:
            [sample_id, [path, path ..]]
        when there's a single path, this:
            [sample_id, path]
        This closure makes both cases:
            [[sample_id, path][sample_id, path]].
        */
            if (it[1].getClass() != java.util.ArrayList){
                // If only one path, `it` will be [sample_id, path]
                return [it]
            }
            l = [];
            for (x in it[1]){
                l.add(tuple(it[0], x))
            }
            return l
        }

        summariseConcatReads(reads)
        sample_ids = summariseConcatReads.out.summary.flatMap({it -> it[0]})

        software_versions = getVersions()
        workflow_params = getParams()

        print('~~~~~~~~~~~~~~~~~meh~~~~~~~~~~~~~~~~~~')
        summariseConcatReads.out.input_reads.view()
        sample_ids.view()

        if (!params.direct_rna){
            preprocess_reads(summariseConcatReads.out.input_reads)
            full_len_reads = preprocess_reads.out.full_len_reads
            pychopper_report = preprocess_reads.out.report.collectFile(keepHeader: true)
        }
        else{
            full_len_reads = summariseConcatReads.out.input_reads
            pychopper_report = file("$projectDir/data/OPTIONAL_FILE")
        }
        if (params.transcriptome_source != "precomputed"){
        
            if (params.transcriptome_source == "denovo"){
                log.info("Doing de novo assembly")
                assembly = denovo_assembly(full_len_reads, ref_genome)

            } else {
                build_minimap_index(ref_genome)
                log.info("Doing reference based transcript analysis")
                assembly = reference_assembly(build_minimap_index.out.index, ref_genome, full_len_reads)
            }
            assembly_stats = assembly.stats.map{ it -> it[1]}.collect()

            split_bam(assembly.bam)
        
            assemble_transcripts(split_bam.out.bundles.flatMap(map_sample_ids_cls), ref_annotation)

            merge_gff_bundles(assemble_transcripts.out.gff_bundles.groupTuple())

            use_ref_ann = !ref_annotation.name.startsWith('OPTIONAL_FILE')

            run_gffcompare(merge_gff_bundles.out.gff, ref_annotation)

            if (params.transcriptome_source == "denovo"){
                // Use the per-sample, de novo-assembled CDS
                seq_for_transcriptome_build = assembly.cds
            }else {
                // For reference based assembly, there is only one reference
                // So map this reference to all sample_ids
                seq_for_transcriptome_build = sample_ids.flatten().combine(Channel.fromPath(params.ref_genome))
            }

            get_transcriptome(
                merge_gff_bundles.out.gff
                .join(run_gffcompare.out.gffcmp_dir)
                .join(seq_for_transcriptome_build))

            gff_compare = run_gffcompare.out.gffcmp_dir.map{ it -> it[1]}.collect()
            merge_gff = merge_gff_bundles.out.gff.map{ it -> it[1]}.collect()
            results = Channel.empty()
        }else
        {
            gff_compare = file("$projectDir/data/OPTIONAL_FILE")
            merge_gff = file("$projectDir/data/OPTIONAL_FILE")
            assembly_stats = file("$projectDir/data/OPTIONAL_FILE")
            use_ref_ann = false
            results = Channel.empty()
        }
        if (jaffal_refBase){
                gene_fusions(full_len_reads, jaffal_refBase, jaffal_genome, jaffal_annotation)
                jaffal_out = gene_fusions.out.results_csv.collectFile(keepHeader: true, name: 'jaffal.csv')
            }else{
                jaffal_out = file("$projectDir/data/OPTIONAL_FILE")
        }


        if (params.de_analysis){

            if (!params.ref_transcriptome){
                merge_transcriptomes(run_gffcompare.output.gtf.collect(), ref_annotation, ref_genome)
                transcriptome = merge_transcriptomes.out.fasta
                gtf = merge_transcriptomes.out.gtf
            }
            else {
                transcriptome = ref_transcriptome
                gtf = Channel.fromPath(ref_annotation)
            }
            check_match = Channel.fromPath(params.condition_sheet)
            check_condition_sheet = check_match.splitCsv(header: true).map{ row -> tuple(
                row.sample_id)
            }
            check_condition_sheet.join(summariseConcatReads.out.input_reads, failOnMismatch: true)
            de = differential_expression(transcriptome, summariseConcatReads.out.input_reads, condition_sheet, gtf)
            de_report = de.all_de
            count_transcripts_file = de.count_transcripts
            dtu_plots = de.dtu_plots
            de_outputs = de.de_outputs
        } else{
            de_report = file("$projectDir/data/OPTIONAL_FILE")
            count_transcripts_file = file("$projectDir/data/OPTIONAL_FILE")
        }
        
        makeReport(
            software_versions,
            workflow_params,
            pychopper_report,
            jaffal_out,
            summariseConcatReads.out.summary.map{it->it[0]}.collect(),
            summariseConcatReads.out.summary.map{it->it[1]}.collect(),
            assembly_stats,
            gff_compare,
            merge_gff,
            de_report,
            count_transcripts_file)

        report = makeReport.out.report
        
        results = results.concat(makeReport.out.report)
       
       if (use_ref_ann){
            results = run_gffcompare.output.gffcmp_dir.concat(
                      assembly.stats,
                      get_transcriptome.out.transcriptome.flatMap(map_sample_ids_cls))
                      .map {it -> it[1]}
                      .concat(results)

       }
    
        if (!use_ref_ann && params.transcriptome_source == "reference-guided"){
            results =  assembly.stats.concat(
                       get_transcriptome.out.transcriptome.flatMap(map_sample_ids_cls))
                      .map {it -> it[1]}
                      .concat(results)

        }
        if (params.transcriptome_source == "denovo"){
            results = assembly.cds.concat(
                       assembly.stats,
                      seq_for_transcriptome_build,
                      get_transcriptome.out.transcriptome.flatMap(map_sample_ids_cls),
                      assembly.opt_qual_ch.flatMap {
                          it ->
                          l = []
                          for (x in it[1..-1]){
                              l.add(tuple(it[0], x))
                          }
                        return l
                      })
                      .map {it -> it[1]}
                      .concat(results)
        }
        if (params.jaffal_refBase){
            results = results
                .concat(gene_fusions.out.results
                .map {it -> it[1]})
        }

        if (params.de_analysis){
            results = results.concat(de.dtu_plots, de_outputs)
        }

    emit:
        results
        telemetry = workflow_params
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.disable_ping == false) {
            Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    fastq = file(params.fastq, type: "file")

    error = null

    if (!fastq.exists()) {
        error = "--fastq: File doesn't exist, check path."
    }

    if (params.transcriptome_source == "precomputed" && !params.ref_transcriptome){
        error = "As transcriptome source parameter is precomputed you must include a ref_transcriptome parameter"
    }
    if (params.transcriptome_source == "reference-guided" && !params.ref_genome){
        error = "As transcriptome source is reference guided you must include a ref_genome parameter"
    }
    if (params.ref_genome){
        ref_genome = file(params.ref_genome, type: "file")
        if (!ref_genome.exists()) {
            error = "--ref_genome: File doesn't exist, check path."
        }
    }else {
        ref_genome = file("$projectDir/data/OPTIONAL_FILE")
    }

    if (params.transcriptome_source == "denovo" && params.ref_annotation) {
        error = "Reference annotation with de denovo assembly is not supported"
    }

    if (params.ref_annotation){
        ref_annotation = file(params.ref_annotation, type: "file")

        if (!ref_annotation.exists()) {
            error = "--ref_annotation: File doesn't exist, check path."
        }
    }else{
        ref_annotation = file("$projectDir/data/OPTIONAL_FILE")
    }
    if (params.jaffal_refBase){
        jaffal_refBase = file(params.jaffal_refBase, type: "dir")
        if (!jaffal_refBase.exists()) {
            error = "--jaffa_refBase: Directory doesn't exist, check path."
        }
     }else{
        jaffal_refBase = null
     }
    ref_transcriptome = file("$projectDir/data/OPTIONAL_FILE")
    if (params.ref_transcriptome){
        log.info("Reference Transcriptome provided will be used for differential expression.")
        ref_transcriptome = file(params.ref_transcriptome, type:"file")
    }
    if (params.de_analysis){
        if (!params.ref_annotation){
            error = "You must provide a reference annotation."
        }
        if (!params.condition_sheet){

            error = "You must provide a condition_sheet or set de_analysis to false."
        }
        condition_sheet = file(params.condition_sheet, type:"file")
    } else{
        condition_sheet = file("$projectDir/data/OPTIONAL_FILE")
    }
    if (error){
        throw new Exception(error)
    }else{
        // reads = fastq_ingress([
        //     "input":params.fastq,
        //     "sample":params.sample,
        //     "sample_sheet":params.sample_sheet])

        //meh 
        nanoseq()

        // pipeline(reads, ref_genome, ref_annotation,
        //     jaffal_refBase, params.jaffal_genome, params.jaffal_annotation,
        //     condition_sheet, ref_transcriptome)

        // output(pipeline.out.results)
    }
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }
    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
