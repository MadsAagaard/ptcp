#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * PTCP_QC — ptcp-qc across all samples in the run (PTCP common.wdl ptcp_qc).
 * Output goes to the run folder's QC/ptcpQC/: qc.<npn>.json per sample,
 * qc.aggregate.json and the CSV summaries, so no per-sample copy is needed
 * with the aggregated (pureTarget) output layout.
 *
 * PTCP builds the ptcp-qc manifest by tarring each sample's outputs and
 * walking the extracted folders. Here the manifest is built directly from
 * the joined per-sample channel, so sample order is identical in every list.
 */

include { ptcpQC } from "../modules/ptcpModules.nf"


workflow PTCP_QC {

    take:
    qcInput     // tuple(meta, [hifiBam, ..., allBam, allBai, ..., sawfishVcf, sawfishTbi]) from PTCP

    main:

    qcInput
        .map { meta, d ->
            // staged file names; ptcpQC links them under the names ptcp-qc expects
            def record = [
                id            : meta.id,
                mappedBam     : d.allBam.name,     // HiFi + fail (see samtools_merge_allReads)
                mappedBai     : d.allBai.name,
                trgtVcf       : d.vcf.name,
                trgtSpanBam   : d.spanBam.name,
                trgtSpanBai   : d.spanBai.name,
                paraphaseBam  : d.paraphaseBam.name,
                paraphaseBai  : d.paraphaseBai.name,
                paraphaseJson : d.paraphaseJson.name,
                f8Json        : d.f8Json.name,
                sawfishVcf    : d.sawfishVcf.name,
                sawfishTbi    : d.sawfishTbi.name,
            ]
            // indexes are staged next to their files so ptcp-qc can open them
            def files = [ d.allBam, d.allBai, d.vcf, d.spanBam, d.spanBai,
                          d.paraphaseBam, d.paraphaseBai, d.paraphaseJson,
                          d.f8Json, d.sawfishVcf, d.sawfishTbi ]
            tuple(record, files)
        }
        .collect(flat: false)                   // emits nothing if no sample survived
        .map { rows ->
            def sorted = rows.sort { it[0].id }
            tuple(sorted.collect { it[0] }, sorted.collect { it[1] }.flatten())
        }
        .set { ptcpQC_input_ch }

    ptcpQC(ptcpQC_input_ch)

    emit:
    json     = ptcpQC.out.json
    csv      = ptcpQC.out.csv
    manifest = ptcpQC.out.manifest
}
