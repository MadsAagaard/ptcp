#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * PREPROCESS — uBAM -> two separate alignments (as in the germline pipeline).
 *   HiFi BAM -> every PTCP step
 *   fail BAM -> TRGT only (--fail-reads), plus the merged copy for ptcp-qc
 * A sample without fail_reads uBAMs is run HiFi-only (warning in the log).
 */

include {
        create_fofn;
        inputFiles_symlinks_ubam;
        pbmm2_align_hifi;
        pbmm2_align_failedOnly;
        } from "../modules/ptcpModules.nf"


workflow PREPROCESS {

    take:
    finalUbamInput          // tuple(meta, [uBAMs])

    main:

    inputFiles_symlinks_ubam(finalUbamInput)
    create_fofn(finalUbamInput)

    pbmm2_align_hifi(create_fofn.out.hifiReads)

    if (params.hifiReads) {

        // HiFi-only mode: no fail alignment, and TRGT runs without --fail-reads.
        fail_ch = Channel.empty()
    }
    else {

        // An empty fail fofn would fail pbmm2; such a sample runs HiFi-only.
        create_fofn.out.failReads
            .branch { meta, fofn ->
                withFail : fofn.size() > 0
                noFail   : true
            }
            .set { failFofn }

        failFofn.noFail.subscribe { meta, fofn ->
            log.warn "Sample ${meta.id}: no fail_reads uBAMs found — TRGT runs on HiFi reads only for this sample."
        }

        pbmm2_align_failedOnly(failFofn.withFail)

        // fail BAM and index as one value, so the join below pads a sample
        // without fail reads with a single null
        fail_ch = pbmm2_align_failedOnly.out.bam
            .map { meta, failBam, failBai -> tuple(meta, [failBam, failBai]) }
    }

    pbmm2_align_hifi.out.bam
        .join(fail_ch, remainder: true)
        .map { meta, hifiBam, hifiBai, fail ->
            tuple(meta, [
                hifiBam : hifiBam,
                hifiBai : hifiBai,
                failBam : fail ? fail[0] : null,
                failBai : fail ? fail[1] : null
            ])
        }
        .set { alignedFinal_ch }

    emit:
    alignedFinal = alignedFinal_ch      // tuple(meta, [hifiBam, hifiBai, failBam, failBai])
}
