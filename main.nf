#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

import java.util.Locale

/*
 * =============================================================================
 * KG Vejle — PacBio PureTarget Carrier Panel (PTCP), local Nextflow implementation
 * =============================================================================
 * Input: one sample sheet, either the SMRT Link CSV used to sequence the run
 * (as the pureTarget pipeline and pacbio.autoLaunch.sh) or LabWare metadata
 * (as the germline pipeline). Meta, channels, naming and read handling
 * follow the KG Vejle germline LRS pipeline. Tool commands follow
 * PacificBiosciences/ptcp v3.2.1; tools run from conda envs (params.condaEnvs).
 * Output is aggregated per run, as in the pureTarget pipeline (conf/ptcp.config).
 * =============================================================================
 */

def runDate      = new Date().format('yyMMdd')
def runDateTime  = new Date().format('yyMMdd HH:mm:ss')
def runUser      = System.getenv('USER') ?: 'unknown'
def runID        = "${runDate}.${runUser}"

def readModeText = params.hifiReads ? 'HiFi only (TRGT runs without --fail-reads)'
                 : 'HiFi + fail reads, aligned separately (fail reads: TRGT --fail-reads only)'

log.info """\
======================================================
Clinical Genetics Vejle: PacBio PureTarget Carrier Panel (PTCP ${params.ptcpVersion})
======================================================
Genome        : ${params.genome}
GenomeDir     : ${params.refFilesDir}
Read mode     : ${readModeText}
Name tag      : ${params.strTag}
RunID         : ${runID}
Script start  : ${runDateTime}
Genome FASTA  : ${params.genomeFasta}
Samplesheet   : ${params.samplesheet}
Env pbmm2     : ${params.condaEnvs.pbmm2}
Env TRGT      : ${params.condaEnvs.trgt51}
Env Paraphase : ${params.condaEnvs.paraphase35}
Env Sawfish   : ${params.condaEnvs.sawfish221}
Env ptcp      : ${params.condaEnvs.ptcp}
lrsShared     : ${params.lrsSharedDir}
lrsCustomTools: ${params.lrsCustomToolsBase}
PTCP scripts  : ${params.ptcpScripts}
TRGT panel    : ${params.ptcpTrgtBed}
Paraphase cfg : ${params.ptcpParaphaseConfig}
ptcp-qc BED   : ${params.ptcpQcBed}
Variant list  : ${params.ptcpVariantList ?: 'none (havanno skipped)'}
SMN model     : ${params.ptcpSmnRegression ?: 'none (SMN homology correction skipped)'}
uBAM search   : ${params.input ?: params.dataArchive}
OutputDirBase : ${params.outputDirBase}
Run output    : ${params.runOutBase}
workDir       : ${workflow.workDir}
error mode    : ${params.errorMode}
min input GB  : ${params.ptcpMinGB}
Reports       : ${params.ptcpReports ? params.ptcpReportScripts : 'off'}
testlist      : must start with '${params.ptcpTestlist}' (SL- prefix ignored); uBAMs matched on NPN + testlist
"""


///////////////////////////////////////////////////
/////// ------- INPUT VALIDATION ------- //////////
///////////////////////////////////////////////////

if (!params.samplesheet) {
    exit 1, """
    USER INPUT ERROR: --samplesheet is required (SMRT Link sample sheet or LabWare metadata).
    """.stripIndent()
}

if (params.genome != 'hg38') {
    exit 1, """
    USER INPUT ERROR: PTCP is implemented for hg38 only (--genome ${params.genome}).
    """.stripIndent()
}

if (params.failedReads) {
    exit 1, """
    USER INPUT ERROR: --failedReads is not supported for PTCP. Use the default
    (HiFi + fail) or --hifiReads.
    """.stripIndent()
}

if (params.aligned) {
    exit 1, """
    USER INPUT ERROR: --aligned is not supported. PTCP starts from uBAMs.
    """.stripIndent()
}

def unsupported = [
    customSS  : params.customSS,
    jointSS   : params.jointSS,
    familySS  : params.familySS,
    jointCall : params.jointCall,
    hpo       : params.hpo,
    singleOnly: params.singleOnly,
    intrefOnly: params.intrefOnly,
].findAll { k, v -> v }.keySet()

if (unsupported) {
    exit 1, """
    USER INPUT ERROR: not supported for PTCP: ${unsupported.collect { '--' + it }.join(', ')}.
    PTCP runs on one sample sheet (--samplesheet): SMRT Link CSV or LabWare metadata.
    """.stripIndent()
}

// PTCP scripts and reference files from the lrsCustomTools release: stop at
// launch if any is missing, instead of failing tasks hours into the run.
def ptcpRequired = [
    "${params.ptcpScripts}/f8_inversion.py",
    "${params.ptcpScripts}/havanno",
    "${params.ptcpScripts}/smn_homology",
    "${params.ptcpScripts}/ptcpqc_json_to_csv.py",
    params.ptcpTrgtBed,
    params.ptcpParaphaseConfig,
    params.ptcpQcBed,
    params.ptcpVariantList,
    params.ptcpSmnRegression,
].findAll { it }
def ptcpMissing = ptcpRequired.findAll { !file(it).exists() }

if (ptcpMissing) {
    exit 1, """
    SETUP ERROR: PTCP files not found (lrsCustomTools release: ${params.lrsCustomToolsBase}):
    ${ptcpMissing.join('\n    ')}
    """.stripIndent()
}

// Reports (lrsCustomTools summaryScripts/ptcp/): missing files disable the
// reports with a warning; they never stop the analysis.
def ptcpReportFiles = [
    "${params.ptcpReportScripts}/ptcp_report_collector.py",
    "${params.ptcpReportScripts}/ptcp_sample_report_template.html",
    "${params.ptcpReportScripts}/ptcp_run_qc_report_template.html",
    params.ptcpReportThresholds,
    params.ptcpReportReference,
].findAll { it }
def ptcpReportMissing = ptcpReportFiles.findAll { !file(it).exists() }
def runReports = params.ptcpReports && !ptcpReportMissing
if (params.ptcpReports && ptcpReportMissing) {
    log.warn "PTCP reports DISABLED — not found in ${params.ptcpReportScripts}:\n    ${ptcpReportMissing.join('\n    ')}"
}

if (params.allReads) {
    log.warn "--allReads is deprecated and ignored. Default is HiFi + fail reads."
}


///////////////////////////////////////////////////
/////// ------- SAMPLESHEET PARSING ------- ///////
///////////////////////////////////////////////////

// Two sample sheet formats, told apart by content:
//
//   SMRT Link CSV (the sheet used to sequence the run, as read by the
//   pureTarget pipeline): every line after the header
//       Bio Sample Name,Plate Well,Adapter,Adapter2
//   is a sample; column 1 is <npn>_<material>_<testlist>_<gender>.
//
//   LabWare metadata (*_metadata.txt, as for the germline pipeline): one sample
//   per line, first tab column
//       rekv_npn_material_testlist_gender_proband_intRef[_expectedCount]
//   Read as pacbio.ingestMetadata.sh reads it: anchored from the right on the
//   proband marker (T/F), so a testlist containing '_' is reassembled intact.
//
// Lines are trimmed (also removes Windows CR); empty lines and LabWare lines
// starting with '#' are skipped.
def SS_HEADER  = 'Bio Sample Name,Plate Well,Adapter,Adapter2'
def ssBaseName = params.ssBase

def sexFromGender = { gender, sampleId ->
    def g = (gender ?: '').toString().trim().toLowerCase()
    if (g in ['m', 'male'])           return 'male'
    if (g in ['k', 'f', 'female'])    return 'female'
    throw new IllegalArgumentException(
        "Sample ${sampleId}: gender field is '${gender}', expected M or K. " +
        "Fix the samplesheet — the pipeline will not guess a karyotype.")
}

// <npn>_<material>_<testlist>_<gender>, testlist may contain '_' (SMRT Link
// Bio Sample Name, and the start of every uBAM file name). An empty gender
// (name ending in '_') is kept as '': the sex always comes from the sample
// sheet, and an empty gender in a SMRT Link sheet stops the run at launch.
// null if fewer than 4 fields or no NPN/testlist.
def parseSampleName = { String name ->
    def t = name.split('_', -1) as List
    if (t.size() < 4) return null
    def r = [ npn: t[0], material: t[1], testlist: t[2..-2].join('_'), gender: t[-1] ]
    (r.npn && r.testlist) ? r : null
}

// Testlist as compared between sample sheet and uBAM names: '_' counts as '-'
// (LabWare testlists sometimes carry '_', see pacbio.ingestMetadata.sh) and a
// leading 'SL-' is dropped (LabWare has SL-PT-KENDT, the SMRT Link Bio Sample
// Name and so the uBAM name has PT-KENDT).
def normTestlist = { String tl -> (tl ?: '').replace('_', '-').replaceFirst(/^SL-/, '') }

// Join key between sample sheet rows and uBAMs: NPN + testlist, as in the
// pureTarget pipeline, so another PureTarget test of the same NPN (or its WGS
// uBAMs) never enters this sample's analysis.
def matchKey = { String npn, String tl -> "${npn}.${normTestlist(tl)}".toString() }

// LabWare metadata line -> map, or null if the proband marker is not found.
def parseLabwareLine = { String line ->
    def f = line.split('\t')[0].trim().tokenize('_')
    def n = f.size()
    def r = null
    if (n >= 7 && f[n - 2] in ['T', 'F']) {          // ... gender proband intRef
        r = [ gender: f[n - 3], proband: f[n - 2], intRef: f[n - 1], expectedCount: null, tlEnd: n - 4 ]
    }
    else if (n >= 8 && f[n - 3] in ['T', 'F']) {     // ... gender proband intRef expectedCount
        r = [ gender: f[n - 4], proband: f[n - 3], intRef: f[n - 2], expectedCount: f[n - 1], tlEnd: n - 5 ]
    }
    if (r == null || r.tlEnd < 3) return null
    [ rekv: f[0], npn: f[1], material: f[2], testlist: f[3..r.tlEnd].join('_'),
      gender: r.gender, proband: r.proband, intRef: r.intRef, expectedCount: r.expectedCount ]
}

def ssLines   = file(params.samplesheet, checkIfExists: true).readLines()*.trim().findAll { it }
def headerIdx = ssLines.findIndexOf { it == SS_HEADER }
def ssFormat  = headerIdx >= 0 ? 'smrtlink' : 'labware'
def badLines  = []
def ssRows

if (ssFormat == 'smrtlink') {
    ssRows = ssLines.drop(headerIdx + 1).collect { line ->
        def cols = line.split(/\s*,\s*/, -1)
        def sn   = parseSampleName(cols[0])
        if (sn == null) { badLines << line; return null }
        [ id       : sn.npn,
          material : sn.material,
          testlist : sn.testlist,
          gender   : sn.gender,
          well     : cols.size() > 1 ? cols[1] : null,
          adapter  : cols.size() > 2 ? cols[2] : null,
          adapter2 : cols.size() > 3 ? cols[3] : null ]
    }
}
else {
    ssRows = ssLines.findAll { !it.startsWith('#') }.collect { line ->
        def lw = parseLabwareLine(line)
        if (lw == null) { badLines << line; return null }
        [ id       : lw.npn,
          material : lw.material,
          testlist : lw.testlist,
          gender   : lw.gender,
          rekv     : lw.rekv,
          proband  : lw.proband,
          intRef   : lw.intRef ]
    }
}

if (badLines) {
    def expected = ssFormat == 'smrtlink'
        ? "Bio Sample Name must be <npn>_<material>_<testlist>_<gender>."
        : "No '${SS_HEADER}' line, so read as LabWare metadata;\n    expected rekv_npn_material_testlist_gender_proband_intRef[_expectedCount]."
    exit 1, """
    USER INPUT ERROR: ${badLines.size()} sample line(s) in ${params.samplesheet} could not be read.
    ${expected}
      ${badLines.take(5).join('\n      ')}
    """.stripIndent()
}

ssRows = ssRows
    .collect { row -> row + [ sex: sexFromGender(row.gender, row.id), ssBase: ssBaseName, ssFormat: ssFormat ] }
    .unique { it.id }      // a sample listed twice is analysed once (published names are keyed on NPN)

if (!ssRows) {
    exit 1, "USER INPUT ERROR: no sample lines in ${params.samplesheet}."
}

log.info "Samplesheet format: ${ssFormat == 'smrtlink' ? 'SMRT Link CSV' : 'LabWare metadata'}, ${ssRows.size()} sample(s)"


///////////////////////////////////////////////////
/////// ------- INPUT DISCOVERY ------- ///////////
///////////////////////////////////////////////////

// -------------------------------------------------------------------------
// uBAM discovery (germline globs, PureTarget archive; --input searches
// another folder). uBAM name: <npn>_<material>_<testlist>_<gender>.<movie>.<readset>.<barcode>.bam
// -------------------------------------------------------------------------
def root     = params.input ?: params.dataArchive
def inputBam = params.hifiReads
    ? (params.input ? "${root}/**/*.hifi_reads.*.bam" : "${root}/**/hifi_reads/*.hifi_reads.*.bam")
    : (params.input ? "${root}/**/*.bam"              : "${root}/**/*_reads/*.bam")

// summary lines, same columns as the germline pipeline
def summaryLine = { meta ->
    def gb = String.format(Locale.US, "%.2f", (meta.totalsizeGB as double))
    "${meta.id}\t${meta.nBams}\t${params.inputReadSet}\t${gb}\t${meta.testlist}"
}
def withHeader = { lines ->
    (["sample\tbamcount\treadSet\ttotal_gb\ttestlist"] + lines).join("\n")
}

channel.fromList(ssRows)
    .map { row -> tuple(matchKey(row.id, row.testlist), row) }
    .set { samplesheet_join }

// uBAMs: only PureTarget testlists (normalised testlist starts with
// params.ptcpTestlist, default 'PT': SL-PT*, PT*), keyed on NPN + testlist.
channel.fromPath(inputBam, followLinks: true)
    .map { bam -> tuple(parseSampleName(bam.baseName.tokenize('.')[0]), bam) }
    .filter { sn, bam -> sn && normTestlist(sn.testlist).startsWith(params.ptcpTestlist) }
    .map { sn, bam -> tuple(matchKey(sn.npn, sn.testlist), bam) }
    .set { ubam_by_id }

samplesheet_join
    .combine(ubam_by_id, by: 0)
    .map { key, metaSS, bam -> tuple(metaSS, bam) }
    .groupTuple(sort: true)
    .map { meta, bams ->
        long   totalBytes = (bams.sum { it.size() } as long)
        double totalGB    = totalBytes / (1024.0 * 1024 * 1024)
        tuple(meta + [nBams: bams.size(), totalsizeGB: totalGB], bams)
    }
    .set { ubam_ss_merged }

// ---- samplesheet rows with no matching uBAM ----------------------------------
samplesheet_join
    .join(ubam_ss_merged.map { meta, bams -> tuple(matchKey(meta.id, meta.testlist), true) }, remainder: true)
    .filter { key, metaSS, found -> metaSS != null && found == null }
    .map { key, metaSS, found -> "${metaSS.id}\t${metaSS.testlist}" }
    .collect()
    .map { lines -> (["sample\ttestlist"] + lines).join("\n") }
    .set { ubam_missing_ch }

// ---- summary lines -------------------------------------------------------------
ubam_ss_merged.map { meta, bams -> summaryLine(meta) }.collect()
    .map(withHeader).set { ubam_size_summary_ch }

// ---- size gate (PTCP-specific threshold, see conf/ptcp.config) ----------------
ubam_ss_merged
    .branch { meta, bams ->
        keep : (meta.totalsizeGB as double) >= (params.ptcpMinGB as double)
        drop : true
    }
    .set { ubam_ss_merged_size_split }

ubam_ss_merged_size_split.drop.map { meta, bams -> summaryLine(meta) }.collect()
    .map(withHeader).set { ubam_size_dropped_ch }

ubam_ss_merged_size_split.keep.map { meta, bams -> summaryLine(meta) }.collect()
    .map(withHeader).set { ubam_size_keep_ch }

// meta.sex was validated at parse time; nothing more to do here.
ubam_ss_merged_size_split.keep.set { finalUbamInput }

// ---- run metadata for the reports ---------------------------------------------
// Summary text (header + rows) -> list of maps; an empty channel -> [].
def summaryRows = { txt ->
    def lines = (txt ?: '').readLines().findAll { it.trim() }
    if (lines.size() < 2) return []
    def head = lines[0].split('\t') as List
    lines.drop(1).collect { l -> [head, (l.split('\t') as List)].transpose().collectEntries { [(it[0] == 'sample' ? 'id' : it[0]): it[1]] } }
}
// SMRT Cell (movie) and run folder per sample, from its uBAMs:
// <npn>_<material>_<testlist>_<gender>.<movie>.<readset>.<barcode>.bam in <run>/<readset>_reads/.
// Kept in its own channel (not in meta), so adding it does not change any task hash.
finalUbamInput
    .map { meta, bams -> [ meta.id,
                           bams.collect { it.name.tokenize('.')[1] }.unique().sort().join(','),
                           bams.collect { it.parent.parent.name }.unique().sort().join(',') ] }
    .set { sample_cells_ch }

// Run folders the analysed uBAMs came from (archive layout <run>/<readset>_reads/<uBAM>)
finalUbamInput
    .flatMap { meta, bams -> bams.collect { it.parent.parent.name } }
    .unique()
    .collect()
    .ifEmpty([])
    .set { run_folders_ch }

ubam_missing_ch.ifEmpty('')
    .combine(ubam_size_dropped_ch.ifEmpty(''))
    .combine(run_folders_ch.map { [it] })
    .map { missingTxt, droppedTxt, folders ->
        [ samplesheet           : params.ssBase,
          samplesheet_format    : ssFormat == 'smrtlink' ? 'SMRT Link CSV' : 'LabWare metadata',
          samples_in_sheet      : ssRows.size(),
          runfolder             : folders.sort().join(', '),
          run_output            : params.runOutBase,
          read_mode             : params.hifiReads ? 'HiFi only' : 'HiFi + fail',
          ptcp_version          : params.ptcpVersion,
          genome_version_tag    : params.genomeVersion,
          tag_hifi              : params.tagHifi,
          tag_str               : params.strTag,
          link_prefix           : '../',
          lrscustomtools_release: file(params.lrsCustomToolsBase).name,
          min_gb                : params.ptcpMinGB,
          tools                 : params.ptcpToolVersions,
          smn_model             : params.ptcpSmnRegression ? file(params.ptcpSmnRegression).name : null,
          missing               : summaryRows(missingTxt),
          dropped               : summaryRows(droppedTxt) ]
    }
    .set { run_meta_ch }


///////////////////////////////////////////////////
/////// ------- MODULES / SUBWORKFLOWS ------- ////
///////////////////////////////////////////////////

include {
        write_input_summary;
        write_dropped_samples_summary;
        write_analyzed_samples_summary;
        write_missing_samples_summary;
        symlinks_ubam_dropped;
        ptcpReports;
        } from './modules/ptcpModules.nf'

include { PREPROCESS } from './subworkflows/PREPROCESS.nf'
include { PTCP }       from './subworkflows/PTCP.nf'
include { PTCP_QC }    from './subworkflows/PTCP_QC.nf'


///////////////////////////////////////////////////
/////// ------- MAIN WORKFLOW ------- /////////////
///////////////////////////////////////////////////

workflow {

    write_input_summary(ubam_size_summary_ch)
    write_analyzed_samples_summary(ubam_size_keep_ch)
    write_dropped_samples_summary(ubam_size_dropped_ch)
    symlinks_ubam_dropped(ubam_ss_merged_size_split.drop)

    ubam_missing_ch
        .filter { txt -> txt.readLines().size() > 1 }
        .set { ubam_missing_nonEmpty_ch }

    ubam_missing_nonEmpty_ch.subscribe { txt ->
        log.warn "Samplesheet samples with no matching uBAM (see runInfo *.missing.samples.summary.txt):\n${txt}"
    }
    write_missing_samples_summary(ubam_missing_nonEmpty_ch)

    PREPROCESS(finalUbamInput)

    PTCP(PREPROCESS.out.alignedFinal)

    PTCP_QC(PTCP.out.qcInput)

    // ---- reports: per sample + run QC -------------------------------------------
    if (runReports) {
        PTCP.out.qcInput
            .map { meta, d -> [meta.id, meta.sex == 'male' ? 'M' : 'K', meta.testlist] }
            .join(sample_cells_ch)                       // + cells, runs
            .collect(flat: false)
            .set { report_samples_ch }

        ptcpReports(
            PTCP_QC.out.json,
            PTCP.out.havanno.map { meta, f -> f }.collect().ifEmpty([]),
            report_samples_ch,
            run_meta_ch)
    }
}


///////////////////////////////////////////////////
/////// ------- COMPLETION ------- ////////////////
///////////////////////////////////////////////////

workflow.onComplete {

    // ---- failure manifest (same as germline) --------------------------------
    // A sample that fails is missing from the run-level ptcp-qc report too.
    def failed = (workflow.stats.failedCount ?: 0) + (workflow.stats.ignoredCount ?: 0)

    if (failed > 0) {
        try {
            def f = file("${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/FAILED_TASKS.PTCP.txt")
            f.parent.mkdirs()
            f.text = """\
                Run          : ${workflow.runName}
                Completed    : ${workflow.complete}
                Succeeded    : ${workflow.stats.succeedCount}
                Failed       : ${workflow.stats.failedCount}
                Ignored      : ${workflow.stats.ignoredCount}
                Error mode   : ${params.errorMode}
                Run output   : ${params.runOutBase}

                Per-task detail is in the trace file under NextflowReports/.
                Output for the affected samples is INCOMPLETE and must not be
                reported without review. Failed samples are also absent from
                the run-level ptcp-qc report.
                """.stripIndent()
            log.warn "*** ${failed} task(s) failed or were ignored — see ${f} ***"
        }
        catch (Exception e) {
            log.warn "Could not write FAILED_TASKS.PTCP.txt: ${e.message}"
        }
    }
}
