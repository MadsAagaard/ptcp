#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * KG Vejle — PTCP modules
 * =============================================================================
 * Commands follow PacificBiosciences/ptcp v3.2.1 (main.wdl, common.wdl,
 * trgt.wdl, paraphase.wdl, sawfish.wdl). Tools run from conda envs, one per
 * process, as in the germline pipeline (params.condaEnvs; ptcp and sawfish221
 * are added in conf/ptcp.config). samtools and bcftools come from PATH.
 * Tools taken from the ptcp env while another env is active are called by
 * full path (${params.condaEnvs.ptcp}bin/): pbindex, minimap2 (HBA step).
 * PacBio's PTCP scripts (f8_inversion.py, havanno, smn_homology,
 * ptcpqc_json_to_csv.py) come from the lrsCustomTools release resolved at
 * launch (params.ptcpScripts) and run as python3 <script> in the step's env.
 *
 * Read sets, as in the germline pipeline: HiFi and fail reads are aligned
 * separately; the fail BAM goes to TRGT only (--fail-reads), every other step
 * uses the HiFi BAM. ptcp-qc and the repeat-reads extract get a merged copy.
 * File prefixes follow lrsFunctions.nf: prefixHifi (HiFi-derived),
 * prefixFail (fail alignment), prefixStr (TRGT and anything using both sets).
 *
 * Output is aggregated per run (pureTarget layout): params.outBase(meta)
 * returns the same run folder for every sample (conf/ptcp.config), so each
 * publishDir below collects all samples of the run in one subfolder. Every
 * published name starts with the NPN, so samples cannot overwrite each other.
 * Only the TRGT plots, which are many SVGs per sample, get a per-sample
 * subfolder (Plots/<npn>/), as in the pureTarget pipeline.
 *
 * Local changes to the PTCP commands, all marked "KG:" in the scripts:
 *   - file names follow the germline prefix convention (<id>.<genomeVersion>.<tag>)
 *   - sex comes from meta.sex (fail-closed) instead of the PTCP default of F
 *   - pbmm2 uses params.genomeMmi (same reference as the germline pipeline)
 *   - trgt plot failures are logged per locus instead of being silently skipped
 *   - the ptcp-qc manifest is built from the joined per-sample channel,
 *     not by walking per-sample tarballs
 * =============================================================================
 */

include { karyotypeFlag;
          f8SexFlag;
          prefixHifi;
          prefixFail;
          prefixStr } from './lrsFunctions.nf'


////////////////////////////////////////////
/////// ------- RUN BOOKKEEPING ------- ////
////////////////////////////////////////////
// Same as the germline pipeline, but published to runInfo only — not to the
// shared lrsDocuments/labWareExport folders, which belong to the germline flow.

process write_input_summary {
    label "low"
    publishDir (
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/"},
        mode: 'copy',
        pattern: "*.txt")

    input:
    val(summary_ch)

    output:
    path("*.txt")

    script:
    """
    printf '%s\\n' '${summary_ch}' > ${params.ssBase}.PTCP.${params.inputReadSet}.input.allSamples.summary.txt
    """
}

process write_dropped_samples_summary {
    label "low"
    publishDir (
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/"},
        mode: 'copy',
        pattern: "*.txt")

    input:
    val(summary_ch)

    output:
    path("*.txt")

    script:
    """
    printf '%s\\n' '${summary_ch}' > ${params.ssBase}.PTCP.${params.inputReadSet}.dropped.samples.summary.txt
    """
}

process write_analyzed_samples_summary {
    label "low"
    publishDir (
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/"},
        mode: 'copy',
        pattern: "*.txt")

    input:
    val(summary_ch)

    output:
    path("*.txt")

    script:
    """
    printf '%s\\n' '${summary_ch}' > ${params.ssBase}.PTCP.${params.inputReadSet}.analyzed.samples.summary.txt
    """
}

// PTCP addition: samplesheet rows for which no uBAM was found (wrong NPN,
// testlist mismatch, data not yet in the archive). Without this they vanish
// from the run without a trace.
process write_missing_samples_summary {
    label "low"
    publishDir (
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/"},
        mode: 'copy',
        pattern: "*.txt")

    input:
    val(summary_ch)

    output:
    path("*.txt")

    script:
    """
    printf '%s\\n' '${summary_ch}' > ${params.ssBase}.PTCP.${params.inputReadSet}.missing.samples.summary.txt
    """
}

process inputFiles_symlinks_ubam {
    label "low"
    tag "$meta.id"
    publishDir (
        path: {"${params.outBase(meta)}/documents/inputSymlinks/"},
        mode: 'symlink',
        pattern: '*.{bam,pbi}',
        overwrite: true)

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path(data)

    script:
    """
    """
}

process symlinks_ubam_dropped {
    label "low"
    tag "$meta.id"

    publishDir (
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/dropped_samples_ubam_symlinks/"},
        mode: 'symlink',
        pattern: '*.{bam,pbi}',
        overwrite: true)

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path(data)

    script:
    """
    """
}


////////////////////////////////////////////
/////// ------- PREPROCESS + ALN ------- ///
////////////////////////////////////////////
// As in the germline pipeline: HiFi and fail reads are aligned separately.
// HiFi BAM -> everything; fail BAM -> TRGT (--fail-reads) only. A merged copy
// is made for ptcp-qc and the repeat-reads extract (see samtools_merge_allReads).
// PTCP instead merges all uBAMs (pbmerge) and aligns once; the alignments are
// the same, since pbmm2 aligns each read independently.

// germline create_fofn
process create_fofn {
    label "low"
    tag "$meta.id"

    publishDir (
        path: {"${params.outBase(meta)}/documents/"},
        mode: 'copy',
        pattern: "*.fofn",
        overwrite: true)

    input:
    tuple val(meta), path(data)   // unmapped BAMs

    output:
    tuple val(meta), path("${meta.id}.fofn"),                emit: allReads
    tuple val(meta), path("${meta.id}.hifi_reads.fofn"),     emit: hifiReads
    tuple val(meta), path("${meta.id}.fail_reads.fofn"),     emit: failReads

    script:
    """
    realpath ${data} > ${meta.id}.fofn
    grep 'fail_reads' ${meta.id}.fofn > ${meta.id}.fail_reads.fofn || true
    grep 'hifi_reads' ${meta.id}.fofn > ${meta.id}.hifi_reads.fofn || true
    """
}

// PTCP common.wdl pbmm2_align_hifi, on the HiFi uBAMs (germline names and fofn input)
process pbmm2_align_hifi {
    label "veryHigh"
    tag "$meta.id"
    conda "${params.condaEnvs.pbmm2}"

    publishDir (
        path: {"${params.outBase(meta)}/alignments/"},
        mode: 'copy',
        pattern: "*.pbmm2.merged.bam*")

    publishDir (
        path: {"${params.outBase(meta)}/QC/pbmm2/"},
        mode: 'copy',
        pattern: "*.pbmm2.{log,alarms.json,mapping_stats.report.json}")

    input:
    tuple val(meta), path(fofn)

    output:
    tuple val(meta),
          path("${prefixHifi(meta)}.pbmm2.merged.bam"),
          path("${prefixHifi(meta)}.pbmm2.merged.bam.bai"),                       emit: bam
    tuple val(meta), path("${prefixHifi(meta)}.pbmm2.merged.bam.pbi"),            emit: pbi
    tuple val(meta), path("${prefixHifi(meta)}.pbmm2.mapping_stats.report.json"), emit: report
    tuple val(meta), path("${prefixHifi(meta)}.pbmm2.log"),                       emit: log
    tuple val(meta), path("${prefixHifi(meta)}.pbmm2.alarms.json"), optional: true, emit: alarms

    script:
    def prefix = prefixHifi(meta)
    // KG: params.genomeMmi instead of the FASTA (same GRCh38_no_alt_analysis_set reference)
    """
    pbmm2 align \
    --log-level ${params.ptcpLogLevel} \
    --log-file ${prefix}.pbmm2.log \
    --alarms ${prefix}.pbmm2.alarms.json \
    --num-threads ${task.cpus} \
    --sort \
    --preset HiFi \
    --report-json ${prefix}.pbmm2.mapping_stats.report.json \
    --sample ${meta.id} \
    --unmapped \
    ${params.genomeMmi} \
    ${fofn} \
    ${prefix}.pbmm2.merged.bam

    ${params.condaEnvs.ptcp}bin/pbindex --num-threads ${task.cpus} ${prefix}.pbmm2.merged.bam
    """
}

// Same command on the fail uBAMs (germline pbmm2_align_failedOnly)
process pbmm2_align_failedOnly {
    label "intermediate"
    tag "$meta.id"
    conda "${params.condaEnvs.pbmm2}"

    publishDir (
        path: {"${params.outBase(meta)}/alignments/failedReads/"},
        mode: 'copy',
        pattern: "*.pbmm2.bam*")

    publishDir (
        path: {"${params.outBase(meta)}/QC/pbmm2/"},
        mode: 'copy',
        pattern: "*.pbmm2.{log,alarms.json,mapping_stats.report.json}")

    input:
    tuple val(meta), path(fofn)

    output:
    tuple val(meta),
          path("${prefixFail(meta)}.pbmm2.bam"),
          path("${prefixFail(meta)}.pbmm2.bam.bai"),                              emit: bam
    tuple val(meta), path("${prefixFail(meta)}.pbmm2.bam.pbi"),                   emit: pbi
    tuple val(meta), path("${prefixFail(meta)}.pbmm2.mapping_stats.report.json"), emit: report
    tuple val(meta), path("${prefixFail(meta)}.pbmm2.log"),                       emit: log
    tuple val(meta), path("${prefixFail(meta)}.pbmm2.alarms.json"), optional: true, emit: alarms

    script:
    def prefix = prefixFail(meta)
    """
    pbmm2 align \
    --log-level ${params.ptcpLogLevel} \
    --log-file ${prefix}.pbmm2.log \
    --alarms ${prefix}.pbmm2.alarms.json \
    --num-threads ${task.cpus} \
    --sort \
    --preset HiFi \
    --report-json ${prefix}.pbmm2.mapping_stats.report.json \
    --sample ${meta.id} \
    --unmapped \
    ${params.genomeMmi} \
    ${fofn} \
    ${prefix}.pbmm2.bam

    ${params.condaEnvs.ptcp}bin/pbindex --num-threads ${task.cpus} ${prefix}.pbmm2.bam
    """
}

// HiFi + fail alignments merged, for the two PTCP steps that need every read:
//  - ptcp-qc: matches TRGT spanning reads to the mapped BAM by read name; a
//    fail read missing from the mapped BAM is logged as an error and dropped
//    from the per-allele read support (ptcp-qc readers/utils.rs). It also
//    reports HiFi vs non-HiFi read counts (rq >= 0.99), as on PTCP's BAM.
//  - the repeat-reads extract (PTCP: reads overlapping repeats, all reads).
// Not published: it duplicates the two alignments above.
process samtools_merge_allReads {
    label "low"
    tag "$meta.id"

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta),
          path("${prefixStr(meta)}.pbmm2.merged.bam"),
          path("${prefixStr(meta)}.pbmm2.merged.bam.bai"), emit: bam

    script:
    def prefix = prefixStr(meta)
    """
    samtools merge \
    --threads ${task.cpus} \
    -c -p \
    -o ${prefix}.pbmm2.merged.bam \
    ${data.hifiBam} ${data.failBam}

    samtools index -@ ${task.cpus} ${prefix}.pbmm2.merged.bam
    """
}


///////////////////////////////////////////////////
/////// ------- REPEAT EXPANSIONS (TRGT) ------- //
///////////////////////////////////////////////////

// PTCP trgt.wdl get_reads_overlapping_repeats, on HiFi + fail reads (data.allBam)
process trgt5_ptcp_repeatReads {
    label "low"
    tag "$meta.id"
    conda "${params.condaEnvs.ptcp}"

    publishDir (
        path: {"${params.outBase(meta)}/repeatExpansions/TRGT/repeatReads/"},
        mode: 'copy',
        pattern: "*.repeats.bam*")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta),
          path("${prefixStr(meta)}.trgt5.ptcp.repeats.bam"),
          path("${prefixStr(meta)}.trgt5.ptcp.repeats.bam.bai"), emit: repeat_reads

    script:
    def prefix = "${prefixStr(meta)}.trgt5.ptcp"
    """
    awk -v sz=${params.ptcpRepeatPad} \
    '{if (\$2 < sz) { print \$1"\\t0\\t"\$3 + sz} else {print \$1"\\t"\$2 - sz"\\t"\$3 + sz} }' \
    ${params.ptcpTrgtBed} \
    > ${prefix}.expanded.bed

    samtools view \
    --write-index \
    --use-index \
    --target-file ${prefix}.expanded.bed \
    --output ${prefix}.repeats.bam##idx##${prefix}.repeats.bam.bai \
    ${data.allBam}
    """
}

// PTCP trgt.wdl trgt_genotype. KG: HiFi BAM + --fail-reads (germline setup)
// instead of one merged BAM. Needs TRGT >= 5.1.0. With --preset targeted TRGT
// already uses reads of any rq, so the reads used are the same as in PTCP.
process trgt5_ptcp {
    label "medium"
    tag "$meta.id"
    conda "${params.condaEnvs.trgt51}"

    publishDir (
        path: {"${params.outBase(meta)}/repeatExpansions/TRGT/"},
        mode: 'copy',
        pattern: "*.trgt5.ptcp.{vcf,vcf.gz,sorted.spanning.bam,sorted.spanning.bam.bai}")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta),
          path("${prefixStr(meta)}.trgt5.ptcp.vcf"),
          path("${prefixStr(meta)}.trgt5.ptcp.sorted.spanning.bam"),
          path("${prefixStr(meta)}.trgt5.ptcp.sorted.spanning.bam.bai"), emit: trgt_full
    tuple val(meta), path("${prefixStr(meta)}.trgt5.ptcp.vcf.gz"),         emit: vcf_gz

    script:
    def karyotype = karyotypeFlag(meta)      // KG: fail-closed; PTCP defaults missing sex to XX
    def failReads = data.failBam ? "--fail-reads ${data.failBam}" : ""
    def prefix    = "${prefixStr(meta)}.trgt5.ptcp"
    """
    trgt --verbose genotype \
    --threads ${task.cpus} \
    --preset targeted \
    --genome ${params.genomeFasta} \
    $karyotype \
    --reads ${data.hifiBam} \
    $failReads \
    --repeats ${params.ptcpTrgtBed} \
    --output-prefix ${prefix}

    gunzip -c ${prefix}.vcf.gz > ${prefix}.vcf
    samtools sort -o ${prefix}.sorted.spanning.bam ${prefix}.spanning.bam
    samtools index ${prefix}.sorted.spanning.bam
    """
}

// PTCP trgt.wdl trgt_plot — one task per (show, plot type), every TRID in the VCF.
// plotMode: [show, plotType], e.g. ['meth', 'waterfall']
// Published per sample (Plots/<npn>/): one run holds hundreds of SVGs per sample.
process trgt5_ptcp_plots {
    label "low"
    tag "$meta.id ${plotMode[0]}_${plotMode[1]}"
    conda "${params.condaEnvs.trgt51}"

    publishDir (
        path: {"${params.outBase(meta)}/repeatExpansions/TRGT/Plots/${meta.id}/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)
    each plotMode

    output:
    tuple val(meta),
          path("${prefixStr(meta)}.trgt5.ptcp.${plotMode[0]}_${plotMode[1]}"),
          path("${prefixStr(meta)}.trgt5.ptcp.${plotMode[0]}_${plotMode[1]}.trgt_plots.zip"), emit: plots

    script:
    def show      = plotMode[0]
    def plotType  = plotMode[1]
    def outPrefix = "${prefixStr(meta)}.trgt5.ptcp.${show}_${plotType}"
    // KG: SVGs go in a per-mode folder (published) and the zip as in PTCP.
    // KG: a failed locus is reported on stderr instead of aborting the loop.
    """
    mkdir -p ${outPrefix}

    bcftools query -f '%INFO/TRID\\n' ${data.vcf} | while read -r i; do
      trgt plot \
      --repeat-id "\$i" \
      --genome ${params.genomeFasta} \
      --repeats ${params.ptcpTrgtBed} \
      --spanning-reads ${data.spanBam} \
      --vcf ${data.vcf} \
      --show ${show} \
      --plot-type ${plotType} \
      --image "${outPrefix}/${outPrefix}.\$(printf '%s' "\$i").trgt_plot.svg" \
      || echo "WARNING: trgt plot failed for \$i (${show} ${plotType})" >&2
    done

    python3 -m zipfile -c ${outPrefix}.trgt_plots.zip ${outPrefix}/*.trgt_plot.svg
    """
}


///////////////////////////////////////////////////
/////// ------- PARAPHASE + F8 ------- ////////////
///////////////////////////////////////////////////

// PTCP paraphase.wdl run_paraphase (+ havanno when params.ptcpVariantList is set), HiFi BAM.
// Paraphase uses only rq >= 0.99 reads, so the HiFi BAM gives the PTCP result.
// Paraphase names its files from the read-group SM tag, i.e. meta.id
// (pbmm2 --sample), so files inside the folder are <id>.paraphase.*
process paraphase_ptcp {
    label "lowCPU"
    tag "$meta.id"
    conda "${params.condaEnvs.paraphase35}"

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/paraphase/"},
        mode: 'copy',
        pattern: "*.{paraphase,havanno.json}")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${prefixHifi(meta)}.paraphase"),                         emit: paraphase_dir
    tuple val(meta),
          path("${prefixHifi(meta)}.paraphase/${meta.id}.paraphase.bam"),
          path("${prefixHifi(meta)}.paraphase/${meta.id}.paraphase.bam.bai"),
          path("${prefixHifi(meta)}.paraphase/${meta.id}.paraphase.json"),         emit: paraphase
    tuple val(meta), path("${prefixHifi(meta)}.havanno.json"), optional: true,      emit: havanno

    script:
    def prefix  = prefixHifi(meta)
    // KG: Paraphase 3.5.0 (PTCP 3.2.1 ships 3.4.0). Same CLI and config format;
    // 3.5 changes the minimap2 realignment settings and long-deletion phasing,
    // and adds genes_in_region to the JSON. Calls can differ from PTCP 3.4.0.
    def havanno = params.ptcpVariantList
        ? "python3 ${params.ptcpScripts}/havanno --variant-vcf ${params.ptcpVariantList} --paraphase-dir ${prefix}.paraphase > ${prefix}.havanno.json"
        : ""
    """
    paraphase \
    --bam ${data.hifiBam} \
    --reference ${params.genomeFasta} \
    --out ${prefix}.paraphase \
    --genome ${params.ptcpGenomeVersion} \
    --threads ${task.cpus} \
    --config ${params.ptcpParaphaseConfig} \
    --write-nocalls-in-vcf \
    --targeted

    ${havanno}
    """
}

// PTCP paraphase.wdl call_f8
process f8inversion_ptcp {
    label "low"
    tag "$meta.id"
    conda "${params.condaEnvs.ptcp}"

    publishDir (
        path: {"${params.outBase(meta)}/specialAnalysis/paraphase/f8inversion/"},
        mode: 'copy')

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta),
          path("${prefixHifi(meta)}.f8inversion.vcf"),
          path("${prefixHifi(meta)}.f8inversion.json"), emit: f8

    script:
    def sex    = f8SexFlag(meta)      // KG: fail-closed; PTCP defaults missing sex to F
    def prefix = prefixHifi(meta)
    """
    mkdir -p out

    python3 ${params.ptcpScripts}/f8_inversion.py \
    --bam ${data.paraphaseBam} \
    --prefix ${prefix} \
    --genome ${params.ptcpGenomeVersion} \
    --json \
    $sex \
    --out out/

    mv out/${prefix}.f8inversion.vcf ${prefix}.f8inversion.vcf
    mv out/${prefix}.f8inversion.json ${prefix}.f8inversion.json
    """
}


///////////////////////////////////////////////////
/////// ------- HBA STRUCTURAL VARIANTS ------- ///
///////////////////////////////////////////////////

// PTCP sawfish.wdl run_sawfish — HBA only: extract HiFi-quality reads in the
// HBA region, remap with minimap2, sawfish discover (no CNV) + joint-call.
// HiFi BAM: the rq >= 0.99 filter would drop the fail reads anyway.
process sawFish2_ptcp_HBA {
    label "intermediate"
    tag "$meta.id"
    conda "${params.condaEnvs.sawfish221}"

    publishDir (
        path: {"${params.outBase(meta)}/structuralVariants/HBA/"},
        mode: 'copy',
        pattern: "*.{sawfishSV.HBA.vcf.gz,sawfishSV.HBA.vcf.gz.tbi,HBA.minimap2.bam,HBA.minimap2.bam.bai}")

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta),
          path("${prefixHifi(meta)}.sawfishSV.HBA.vcf.gz"),
          path("${prefixHifi(meta)}.sawfishSV.HBA.vcf.gz.tbi"),   emit: sv_vcf
    tuple val(meta),
          path("${prefixHifi(meta)}.HBA.minimap2.bam"),
          path("${prefixHifi(meta)}.HBA.minimap2.bam.bai"),       emit: minimap2_bam

    script:
    def prefix = prefixHifi(meta)
    """
    samtools view -h -e "[rq]>=0.99" -T ${params.genomeFasta} ${data.hifiBam} ${params.ptcpHbaRegion} | \
    samtools fastq -T MM,ML - | \
    ${params.condaEnvs.ptcp}bin/minimap2 -a -y -x map-hifi -Y -E 1,0 ${params.genomeFasta} - | \
    samtools view -b -h | \
    samtools sort > ${prefix}.HBA.minimap2.bam

    samtools index ${prefix}.HBA.minimap2.bam

    sawfish discover \
    --disable-cnv \
    --threads ${task.cpus} \
    --ref ${params.genomeFasta} \
    --bam ${prefix}.HBA.minimap2.bam \
    --output-dir ${meta.id}.sawfishDiscover.HBA

    sawfish joint-call \
    --threads ${task.cpus} \
    --sample ${meta.id}.sawfishDiscover.HBA \
    --output-dir ${meta.id}.sawfishSV.HBA

    mv ${meta.id}.sawfishSV.HBA/genotyped.sv.vcf.gz     ${prefix}.sawfishSV.HBA.vcf.gz
    mv ${meta.id}.sawfishSV.HBA/genotyped.sv.vcf.gz.tbi ${prefix}.sawfishSV.HBA.vcf.gz.tbi
    """
}


///////////////////////////////////////////////////
/////// ------- QUALITY CONTROL (ptcp-qc) ------- /
///////////////////////////////////////////////////

// PTCP common.wdl ptcp_qc — run-level, all samples together (per-sample
// qc.<id>.json + qc.aggregate.json + CSV summaries), published to the run
// folder's QC/ptcpQC/ next to the per-sample pbmm2 reports.
// records: one map per sample with the staged file names (see PTCP_QC.nf)
process ptcpQC {
    label "intermediate"
    conda "${params.condaEnvs.ptcp}"

    publishDir (
        path: {"${params.runOutBase}/QC/ptcpQC/"},
        mode: 'copy',
        saveAs: { fn -> file(fn).name })

    input:
    tuple val(records), path(files)

    output:
    path("ptcp_qc/qc.*.json"),                          emit: json
    path("ptcp_qc/qc.*.csv"),                           emit: csv
    path("${params.ssBase}.PTCP.ptcpQC.manifest.json"), emit: manifest

    script:
    def manifestName = "${params.ssBase}.PTCP.ptcpQC.manifest.json"
    // ptcp-qc takes the sample name from each file name by stripping a fixed
    // suffix (<sample>.mapped.bam, <sample>.trgt.vcf, ...; ptcp-qc
    // sample_context.rs). A required file (mapped BAM) with another name stops
    // it; an optional one (TRGT, Paraphase, F8, Sawfish) is skipped with only a
    // warning, so its QC would silently be missing. The staged files are
    // therefore linked under those names in qc_input/, as PTCP names them.
    def links = []
    def link  = { String src, String id, String suffix ->
        links << "ln -s ../${src} qc_input/${id}${suffix}"
        "qc_input/${id}${suffix}".toString()
    }
    def manifest = [ sample_names: [], mapped_bams: [], trgt_vcfs: [], trgt_spanning_bams: [],
                     paraphase_bams: [], paraphase_jsons: [], f8_jsons: [], sawfish_vcfs: [] ]
    records.each { r ->
        manifest.sample_names       << r.id
        manifest.mapped_bams        << link(r.mappedBam,     r.id, '.mapped.bam')
        link(r.mappedBai,     r.id, '.mapped.bam.bai')
        manifest.trgt_vcfs          << link(r.trgtVcf,       r.id, '.trgt.vcf')
        manifest.trgt_spanning_bams << link(r.trgtSpanBam,   r.id, '.trgt.sorted.spanning.bam')
        link(r.trgtSpanBai,   r.id, '.trgt.sorted.spanning.bam.bai')
        manifest.paraphase_bams     << link(r.paraphaseBam,  r.id, '.paraphase.bam')
        link(r.paraphaseBai,  r.id, '.paraphase.bam.bai')
        manifest.paraphase_jsons    << link(r.paraphaseJson, r.id, '.paraphase.json')
        manifest.f8_jsons           << link(r.f8Json,        r.id, '.f8inversion.json')
        manifest.sawfish_vcfs       << link(r.sawfishVcf,    r.id, '.sv.vcf.gz')
        link(r.sawfishTbi,    r.id, '.sv.vcf.gz.tbi')
    }
    // Same keys as the PTCP manifest. Paraphase VCFs are omitted, as in PTCP
    // (ptcp-qc does not consume them yet).
    def manifestJson = groovy.json.JsonOutput.toJson(manifest)
    def linkCmds     = links.join('\n    ')
    // KG: smn_homology exits non-zero when any file fails (e.g. a model with
    // another feature count). That must not fail ptcp-qc or the run: warn,
    // keep the unannotated JSONs; the reports then show "model supplied, no
    // result" for the affected samples. test-batch writes each file atomically.
    def smn = params.ptcpSmnRegression
        ? """if ! python3 ${params.ptcpScripts}/smn_homology test-batch --model ${params.ptcpSmnRegression} --input-dir ptcp_qc --in-place; then
        echo "WARNING: smn_homology failed for one or more samples (see above) - SMN homology adjustment missing for those" >&2
    fi"""
        : ""
    """
    mkdir -p ptcp_qc qc_input
    ${linkCmds}

    printf '%s\\n' '${manifestJson}' > ${manifestName}

    ptcp-qc -v analyze \
    --manifest ${manifestName} \
    --targets-bed ${params.ptcpQcBed} \
    --output-prefix ptcp_qc/qc \
    --threads ${task.cpus} \
    --json

    ${smn}

    python3 ${params.ptcpScripts}/ptcpqc_json_to_csv.py --inputs ptcp_qc/qc*.json \
    -o ptcp_qc/qc.csv
    """
}


///////////////////////////////////////////////////
/////// ------- REPORTS ------- ///////////////////
///////////////////////////////////////////////////

// One HTML report per sample (that sample's QC and results only, plus anonymous
// run medians) and one run QC report (all samples, QC and callability, no
// results), built by lrsCustomTools summaryScripts/ptcp/ptcp_report_collector.py
// from the ptcp-qc JSONs and the havanno JSONs. Published to <run>/reports/.
// sampleRows: [[id, sex M|K, testlist, cells, runs], ...]; runMeta: map (see main.nf).
// Target shares are compared with params.ptcpReportReference (lab reference,
// ptcp_build_reference.py) or, without it, with the sample's own SMRT Cell.
process ptcpReports {
    label "low"
    conda "${params.condaEnvs.ptcp}"

    publishDir (
        path: {"${params.runOutBase}/reports/"},
        mode: 'copy',
        saveAs: { fn -> file(fn).name })

    input:
    path(qcJsons, stageAs: 'ptcp_qc/*')     // qc.<npn>.json + qc.aggregate.json
    path(havanno, stageAs: 'havanno/*')     // may be empty (havanno skipped)
    val(sampleRows)
    val(runMeta)

    output:
    path("reports/*.html"), emit: html
    path("reports/*.json"), emit: json

    script:
    // sampleRows: [id, sex, testlist, cells, runs] (cells / runs comma-separated)
    def tsv      = (["id\tsex\ttestlist\tcells\truns"] + sampleRows.collect { it.join('\t') }).join('\n')
    def metaJson = groovy.json.JsonOutput.toJson(runMeta).replace("'", "'\\''")
    def thr      = params.ptcpReportThresholds ? "--thresholds ${params.ptcpReportThresholds}" : ""
    def ref      = params.ptcpReportReference  ? "--reference ${params.ptcpReportReference}" : ""
    def scripts  = params.ptcpReportScripts
    """
    printf '%s\\n' '${tsv}' > samples.tsv
    printf '%s\\n' '${metaJson}' > run_meta.json

    python3 ${scripts}/ptcp_report_collector.py \
    --qc-dir ptcp_qc \
    --havanno \$(ls havanno/*.json 2>/dev/null) \
    --samples samples.tsv \
    --run-meta run_meta.json \
    ${thr} ${ref} \
    --sample-template ${scripts}/ptcp_sample_report_template.html \
    --run-template ${scripts}/ptcp_run_qc_report_template.html \
    --outdir reports
    """
}
