# KG Vejle PacBio PureTarget Carrier Panel (PTCP) pipeline

## General info

Local Nextflow implementation of PacBio's PureTarget Carrier Pipeline
([PacificBiosciences/ptcp](https://github.com/PacificBiosciences/ptcp) v3.2.1),
used at Clinical Genetics, Vejle.

- Input is one sample sheet: the SMRT Link CSV of the run (as for the
  pureTarget pipeline and `pacbio.autoLaunch.sh puretarget`) or LabWare
  metadata (as for the germline pipeline). uBAMs are searched in the
  PureTarget archive (`pureTarget_rawdata`).
- Meta map, channel structure, read handling and file naming follow the KG
  Vejle germline LRS pipeline (`MadsAagaard/pacbioGermline`).
- Output is aggregated per run, as in the pureTarget pipeline: one folder per
  sample sheet, with one subfolder per output type holding every sample.
- Tool commands follow the PTCP WDL tasks. Tools run from conda envs, one per
  process, as in the germline pipeline; samtools and bcftools come from PATH.
- PacBio's PTCP scripts and reference files come from an lrsCustomTools
  release (`ptcp/pacbio/`, `resources/ptcp/`), vendored there unmodified from
  PTCP 3.2.1. See "Scripts and reference files" below.
- Shared settings (reference, roots, labels, slurm profile, error strategy,
  condaEnvs) come from the same shared config as the germline pipeline:
  `/data/configs/nextflowLRS_refactor_v1.config`. PTCP-only settings are in
  `conf/ptcp.config`.

## Read handling

As in the germline pipeline, HiFi and fail reads are **aligned separately**.
The fail reads are used by TRGT only; every other step uses the HiFi BAM:

```
uBAMs ─ create_fofn ─┬─ pbmm2_align_hifi ───────┬─ trgt5_ptcp (--reads) ─ trgt5_ptcp_plots (x4)
                     │                          ├─ paraphase_ptcp ─ f8inversion_ptcp
                     │                          ├─ sawFish2_ptcp_HBA
                     │                          └─┐
                     └─ pbmm2_align_failedOnly ─┬─┼─ trgt5_ptcp (--fail-reads)
                                                └─┴─ samtools_merge_allReads ─┬─ trgt5_ptcp_repeatReads
                                                                              └─ ptcpQC (whole run) ─ ptcpReports
```

PTCP itself merges all uBAMs and aligns once. The results are the same:

- pbmm2 aligns each read independently, so two alignments hold the same
  records as one.
- TRGT with `--preset targeted` uses reads of any quality, so `--reads` HiFi
  plus `--fail-reads` gives it the same reads as PTCP's merged BAM. This needs
  TRGT 5.1.0, the first version with `--fail-reads`.
- Paraphase and the Sawfish HBA step use only `rq >= 0.99` reads, so the HiFi
  BAM is what they used in PTCP as well.
- ptcp-qc gets HiFi + fail merged (`samtools_merge_allReads`, not published).
  It matches TRGT's spanning reads to the mapped BAM by read name; with the
  HiFi BAM alone, every fail read TRGT used would be logged as an error and
  dropped from the per-allele read counts, and the HiFi/non-HiFi read counts
  would change. The PTCP "reads overlapping repeats" BAM is cut from the same
  merged copy.

A sample without fail_reads uBAMs runs HiFi-only (TRGT without `--fail-reads`,
ptcp-qc on the HiFi BAM), with a warning in the log. `--hifiReads` does the
same for every sample and skips the fail alignment.

## Tool versions

```
Step        Env                             Version     PTCP 3.2.1
pbmm2       condaEnvs.pbmm2                 1.16.99     SMRT Link 25.1   (26.1.99 with --latestVersions)
TRGT        condaEnvs.trgt51                5.1.0       5.0.0            (5.1.0 adds --fail-reads)
Paraphase   condaEnvs.paraphase35           3.5.0       3.4.0
Sawfish     condaEnvs.sawfish221 *          2.2.1       2.2.1
ptcp-qc     condaEnvs.ptcp *                1.0.0       1.0.0
PTCP scripts (f8_inversion.py, havanno, smn_homology, ptcpqc_json_to_csv.py)
            lrsCustomTools ptcp/pacbio/     PTCP tag 3.2.1   (python3 of the step's env)
pbindex, minimap2 (HBA step)
            condaEnvs.ptcp *                pbtk / minimap2 from bioconda
samtools, bcftools                          from PATH
```

TRGT 5.1.0 is needed for `--fail-reads`; its changelog lists only that and
BCF output support as changes from 5.0.0. Paraphase is the other version that
differs from PTCP. 3.5.0 has the same command
line and config format as 3.4.0, but changes the minimap2 settings used for
realignment, improves phasing across long deletions, and adds a
`genes_in_region` field to the JSON. Calls can therefore differ from PTCP 3.4.0
and should be covered by validation.

Env names follow the shared config's `params.condaEnvs`. The two marked * are
not in the shared config yet and are added by `conf/ptcp.config`; they can move
to the shared config's `condaEnvs` block unchanged. `sawfish221` is needed
because the shared `condaEnvs.sawfish2` is 2.1.1 unless `--latestVersions`.

## Scripts and reference files

The four PTCP scripts and the PTCP reference files (TRGT panel, Paraphase
config, ptcp-qc BED, variant list) are read from an lrsCustomTools release:

```
<release>/ptcp/pacbio/                 f8_inversion.py, havanno, smn_homology, ptcpqc_json_to_csv.py
<release>/resources/ptcp/v3.2.1/       trgt/, paraphase/, ptcp-qc/, variant_list/
```

Both paths hang off `params.lrsCustomToolsBase` from the shared config, the same
base the germline custom scripts use. The shared config resolves
`/data/mmaj/scripts_local/lrsCustomTools/releases/current` when it loads, as the
lrsCustomTools README prescribes:

- each run is pinned to the release that was current at launch, even if a new
  release is made while it runs;
- command lines name the real release folder, so after a new release `-resume`
  reruns every task that uses a PTCP script or reference file (TRGT, Paraphase,
  F8, ptcp-qc) instead of reusing old results; merging and alignment stay cached;
- a broken `current` link stops the run at launch, and `main.nf` also checks
  at launch that every script and reference file exists.

The launch log shows the release used (`lrsCustomTools:` line).
`--lrsCustomToolsBase <folder>` runs against another checkout instead, e.g. a
`lrsCustomToolsPull.sh test` folder when testing script changes.

## Analysis steps

| Step | Process | PTCP source |
|---|---|---|
| Split uBAMs into HiFi / fail | `create_fofn` | common.wdl `collect_inputs` |
| Alignment, HiFi and fail separately | `pbmm2_align_hifi`, `pbmm2_align_failedOnly` | common.wdl `pbmm2_align_hifi` |
| HiFi + fail merged (ptcp-qc, repeat reads) | `samtools_merge_allReads` | (local) |
| Reads overlapping repeats (±3 kb) | `trgt5_ptcp_repeatReads` | trgt.wdl `get_reads_overlapping_repeats` |
| Repeat genotyping | `trgt5_ptcp` | trgt.wdl `trgt_genotype` |
| Motif/methylation, allele/waterfall plots | `trgt5_ptcp_plots` | trgt.wdl `trgt_plot` |
| Paraphase (+ havanno) | `paraphase_ptcp` | paraphase.wdl `run_paraphase` |
| F8 intron 1/22 inversions | `f8inversion_ptcp` | paraphase.wdl `call_f8` |
| HBA structural variants | `sawFish2_ptcp_HBA` | sawfish.wdl `run_sawfish` |
| ptcp-qc (+ SMN homology, CSV) | `ptcpQC` | common.wdl `ptcp_qc` |
| Per-sample reports + run QC report (HTML) | `ptcpReports` | (local, lrsCustomTools `summaryScripts/ptcp/`) |

---

# Setup (once)

1. Build the `ptcp` conda env (the Python packages the PTCP scripts import,
   pbtk (pbindex), minimap2, and the ptcp-qc 1.0.0 binary with its MD5 checked):
   ```
   bash env/build_ptcp_env.sh
   ```
   Default location: `/lnx01_data3/shared/programmer/miniconda3/envs/ptcp`
   (override with `ENV=/path`). The script records the ptcp-qc version and the
   full package list in `<env>/ptcp_versions.txt`. Needs network access to the
   conda channels and github.com.
2. Deploy an lrsCustomTools release that contains `ptcp/` (tag `ptcp-v3.2.1`
   or later):
   ```
   lrsCustomToolsPull.sh release
   ```
3. Optional: SMN1/2 homology model → `--ptcpSmnRegression /path/model.json`.
   PTCP 3.2.1 ships neither a model nor its training JSONs (the smn_homology
   README says `ptcp-qc-outs/` and `puretarget-v8.json` are included; they are
   not). The model has to be trained with `smn_homology train` on ptcp-qc
   JSONs. Use the `.json` model format: it stores only the regression
   coefficients and does not depend on the scikit-learn version (`.pkl`
   does). The model file can live in lrsCustomTools `resources/ptcp/`; the
   training JSONs are patient data and must not. See "SMN1/2 homology model"
   under Reports for how the model is used and how to train one.

4. Reports: the collector and the two HTML templates must be in the
   lrsCustomTools release under `summaryScripts/ptcp/`
   (`ptcp_report_collector.py`, `ptcp_sample_report_template.html`,
   `ptcp_run_qc_report_template.html`). If they are missing the run still
   completes, without reports, and the launch log says so.
5. Optional: target-share reference → `--ptcpReportReference /path/ref.json`,
   built once from validation samples with `ptcp_build_reference.py` (see
   Reports). The file holds medians and ranges only, no sample names, so it
   can live in lrsCustomTools `resources/ptcp/`.

Nothing in the shared config needs to change (see Tool versions for two env
names that could move there).

---

# Usage

## Samplesheet

Either format; the pipeline tells them apart by content and logs which one it
read (`Samplesheet format:` line). Lines are trimmed (Windows line endings are
fine) and empty lines dropped.

**SMRT Link CSV** — the sheet used to sequence the run. Every line after the
header

```
Bio Sample Name,Plate Well,Adapter,Adapter2
```

is a sample; the Bio Sample Name is `<npn>_<material>_<testlist>_<gender>`.

**LabWare metadata** (`*_metadata.txt`) — used when the sheet has no SMRT Link
header. One sample per line (first tab column), lines starting with `#`
skipped:

```
rekv_npn_material_testlist_gender_proband_intRef[_expectedCount]
0000012345_123456789012_78_SL-PTCP_K_T_noInfo
```

Read as `pacbio.ingestMetadata.sh` reads it: anchored from the right on the
proband marker (`T`/`F`), so both the 7- and 8-field forms work and a testlist
containing `_` is kept whole. proband and intRef are read but not used (no
family analysis in PTCP).

In both formats gender must be M or K. A line that cannot be read, or any
other gender, stops the run at launch with the offending lines listed. A
sample listed twice is analysed once.

uBAMs are found in the PureTarget archive or in `--input` and matched to a
sheet row on NPN + testlist, so WGS uBAMs and other PureTarget tests of the
same NPN never enter the carrier analysis. Testlists are compared with `_` as
`-` and a leading `SL-` dropped: the LabWare testlist `SL-PT-KENDT` matches the
uBAM name testlist `PT-KENDT` (the SMRT Link Bio Sample Name has no `SL-`).
Only uBAMs whose testlist then starts with `PT` (`params.ptcpTestlist`) are
used. The gender in the uBAM name is not used and may be empty
(`<npn>_<material>_<testlist>_.m84...bam`); the sex always comes from the sheet.
Samplesheet rows with no matching uBAM are written to
`runInfo/<date>_<samplesheet>/*.missing.samples.summary.txt` and logged.

## Options

```
--samplesheet   [path]  SMRT Link CSV or LabWare metadata. Required
--input         [path]  Folder to search for uBAMs. Default: the PureTarget archive
--hifiReads     [bool]  HiFi only (default: HiFi + fail, aligned separately)
--ptcpTestlist  [str]   testlist must start with this, SL- prefix ignored (default PT)
--ptcpMinGB     [num]   Minimum combined uBAM size per sample (default 0 = keep all)
--ptcpSmnRegression [path]  SMN homology model; omitted = correction skipped
--ptcpVariantList   [path]  havanno variant list; set to null to skip havanno
--lrsCustomToolsBase [path] lrsCustomTools checkout to use (default: releases/current)
--ptcpReports   [bool]  HTML reports (default true; false = skip)
--ptcpReportThresholds [path]  approved QC/repeat thresholds (YAML/JSON); omitted = placeholders
--ptcpReportReference  [path]  target-share reference (ptcp_build_reference.py); omitted = same-cell comparison
--outdir        [path]  Output base (shared-config behaviour)
--errorMode     [str]   strict (default) | cohort
-profile slurm          Run on the KG Vejle SLURM cluster
```

Not supported (the run stops with an error): `--customSS`, `--jointSS`,
`--familySS`, `--jointCall`, `--hpo`, `--singleOnly`, `--intrefOnly`,
`--failedReads`, `--aligned`.

## Examples

```
nextflow run MadsAagaard/<ptcp-repo> -r main --samplesheet /path/to/run_samplesheet.csv -profile slurm
nextflow run MadsAagaard/<ptcp-repo> -r main --samplesheet /path/to/<name>_metadata.txt -profile slurm
nextflow run MadsAagaard/<ptcp-repo> -r main --samplesheet /path/to/run_samplesheet.csv --input /path/to/run/folder -profile slurm
```

---

# Output

Aggregated per run, as in the pureTarget pipeline: one folder per sample
sheet, and in it one subfolder per output type holding every sample of the run:

`<outputDirBase>/pureTargetCarrier/<samplesheet>_<readSet>/`

(`params.runOutBase` in `conf/ptcp.config`, the same for both sheet formats;
`params.outBase(meta)` returns it for every sample, so the publishDirs keep the
germline form). The launch log shows it on the `Run output` line. Sample
identity is in the file names, which follow the germline convention
`<npn>.hg38v4LRS.<tag>`: `HifiReads` for HiFi-derived output, `failedReads` for
the fail alignment, and `AllReadsNew` for TRGT, which uses both
(`HifiReadsOnly` with `--hifiReads`). Every name starts with the NPN and the
sheet is deduplicated on NPN at launch, so samples cannot overwrite each other.

```
alignments/                             <npn>.hg38v4LRS.HifiReads.pbmm2.merged.bam(.bai/.pbi)
alignments/failedReads/                 <npn>.hg38v4LRS.failedReads.pbmm2.bam(.bai/.pbi)
documents/                              <npn>.fofn, <npn>.hifi_reads.fofn, <npn>.fail_reads.fofn
documents/inputSymlinks/                uBAM symlinks
QC/pbmm2/                               pbmm2 logs and mapping reports, HiFi and fail
QC/ptcpQC/                              qc.<npn>.json, qc.aggregate.json, qc.{trgt,paraphase,coverage,aggregate}.csv,
                                        <samplesheet>.PTCP.ptcpQC.manifest.json
repeatExpansions/TRGT/                  ...AllReadsNew.trgt5.ptcp.vcf(.gz), .sorted.spanning.bam(.bai)
repeatExpansions/TRGT/Plots/<npn>/      ...AllReadsNew.trgt5.ptcp.<motifs|meth>_<allele|waterfall>/*.svg + .zip
repeatExpansions/TRGT/repeatReads/      ...AllReadsNew.trgt5.ptcp.repeats.bam(.bai)
specialAnalysis/paraphase/              ...HifiReads.paraphase/ (one folder per sample), ...HifiReads.havanno.json
specialAnalysis/paraphase/f8inversion/  ...HifiReads.f8inversion.vcf/.json
structuralVariants/HBA/                 ...HifiReads.sawfishSV.HBA.vcf.gz(.tbi), ...HifiReads.HBA.minimap2.bam(.bai)
reports/                                <npn>.PTCP.sampleReport.html (+ .json), <samplesheet>.PTCP.runQC.html (+ .json)
```

The TRGT plots are the one per-sample subfolder (`Plots/<npn>/`, as in the
pureTarget pipeline): four plot sets of one SVG per locus per sample would
otherwise put thousands of files in one folder.

Run bookkeeping, under `<outputDirBase>/runInfo/<date>_<samplesheet>/` (sheet
name without `.csv`/`.txt`), as in the germline pipeline:

```
*.PTCP.*.summary.txt             input / analyzed / dropped / missing samples
dropped_samples_ubam_symlinks/   uBAMs of samples dropped by the size gate
NextflowReports/                 timeline, report, trace, DAG
FAILED_TASKS.PTCP.txt            only if tasks failed or were ignored
```

---

# Reports

`ptcpReports` runs after ptcp-qc and writes to `reports/` in the run folder:

- `<npn>.PTCP.sampleReport.html`, one per sample, for the clinical
  interpreter. It holds that sample's QC and every result for that sample:
  SMN1/SMN2, GBA1, CYP21A2, HBB, RPGR (Paraphase + havanno), HBA (Paraphase +
  Sawfish), F8 intron 1/22 inversions, and FMR1, FXN, AFF2, ARX repeats. It
  shows nothing about any other sample: the only figures from other samples
  are anonymous comparison values (see Target share).
- `<samplesheet>.PTCP.runQC.html`, one per run, for the lab and bioinformatics.
  It names every sample, so it carries **no results**: reads, on-target,
  read quality, the sex check (sheet vs X/Y ratios), reads per carrier target
  and whether each caller produced a call, plus run integrity checks (sheet vs
  ptcp-qc, missing and dropped samples). Samples are grouped by the SMRT Cell
  their uBAMs came from; samples topped up on several cells are listed as
  "multiple cells".

The `.json` next to each HTML is the payload the HTML was built from. The
collector emits numbers only; every threshold, verdict and sentence is in the
templates. Until `--ptcpReportThresholds` points at a file with `approved_by`
and `approved_date`, every threshold is shown as a placeholder.

Links in the reports (BAMs, VCFs, TRGT plots) are relative to `reports/`, so
they work when the report is opened inside the run folder.

## Target share

Absolute reads per target depend on how deeply the sample was sequenced, so
they are not comparable between a sample alone on a cell, one in a full pool,
or one topped up on a second run. Each sample report therefore also shows the
**target share**: the target's reads as a fraction of the sample's on-target
reads (HiFi + fail, as ptcp-qc counts them). A share below 0.5× the comparison
(`share_ratio_warn`, placeholder) sets the sample's QC to review.

The comparison is, in order:

1. the lab reference given with `--ptcpReportReference` (median and range per
   target over validation samples);
2. without a reference: the median and **interquartile range** of the other
   samples of the same SMRT Cell, this sample left out, when there are at
   least 5 of them (collector `--batch-min`). Min/max are not shown, because
   an extreme is one identifiable sample's value;
3. otherwise no comparison; only the absolute read thresholds apply.

SMRT Cell and run folder come from the uBAM names and paths
(`<npn>_..._.<movie>.hifi_reads.bc*.bam` in `<runfolder>/hifi_reads/`); they
are passed to the reports in a separate channel, not in `meta`, so `-resume`
is not affected.

Building a reference (after ptcp-qc has run on the validation samples; any
number of runs):

```
python3 ptcp_build_reference.py --name kgvejle_ptcp_ref_v1 \
    --out ptcp_target_share_reference.v1.json \
    <run1>/QC/ptcpQC/qc.*.json <run2>/QC/ptcpQC/qc.*.json [--exclude <npn> ...]
```

`qc.aggregate.json` among the inputs is skipped. Samples below 2000 on-target
reads (`--min-on-target`), without reads at any carrier target, or listed with
`--exclude` are left out; the per-sample shares are printed to stderr for
review. At least 3 samples are needed. A reference is only valid for the
ptcp-qc version and targets BED it was built with; the sample report warns on
a version mismatch.

## SMN1/2 homology model

Paraphase counts SMN1 and SMN2 copies from the distinct haplotypes it sees.
When every SMN1 copy is identical and every SMN2 copy is identical
(`is_possible_homology`), a 1:1 sample cannot be told apart from 2:2 by
haplotypes, and 2:0 from 4:0. `smn_homology` predicts the expected reads per
copy from the other targets' depth and re-estimates the copy numbers.

The sample report:

- uses the Paraphase call by default;
- uses the depth model's copy numbers **only when `is_possible_homology` is
  true**, and says so;
- when Paraphase saw distinct haplotypes and the model still disagrees, keeps
  the Paraphase call and sets it to review;
- shows the model's result (or "not run" / "no result") next to the call.

In `ptcpQC` a failing `smn_homology` (for example a sample missing a target)
is a warning in the task log, not a failed run; the affected sample report
shows "no result".

Model details (`smn_homology` in PTCP 3.2.1): a linear regression without
intercept on 38 features: on-target HiFi reads, total HiFi reads, and per
locus the coverage, its square, and coverage divided by each read total, over
9 loci. AFF2 and ARX are named FRAXE_AFF2 / PRTS_ARX in the script but AFF2 /
ARX in ptcp-qc, so they never enter the model. The model JSON stores
coefficients by position only, so a model must be used with the same script
version it was trained with. Training uses only samples with ≥ 2000 on-target
HiFi reads, no possible homology, and every haplotype > 5 reads per copy.

38 coefficients from 10 samples are underdetermined. The recommended route:
ask PacBio for the training JSONs behind `puretarget-v8.json`, train on theirs
plus ours, validate on our samples against MLPA, and store the model in
lrsCustomTools `resources/ptcp/`.

# Local changes to the PTCP commands

Marked `KG:` in `modules/ptcpModules.nf`.

- **HiFi and fail reads aligned separately**, fail reads to TRGT through
  `--fail-reads` (TRGT 5.1.0), as in the germline pipeline. PTCP merges all
  uBAMs and aligns once. See "Read handling" for why the results are the same.
- **Sex is fail-closed.** PTCP treats a missing or unrecognised sex as female
  (TRGT XX, F8 `--sex F`). Here sex comes from `meta.sex`, validated at
  samplesheet parse time; anything other than M/K stops the run.
- **Tools from conda envs** instead of the PTCP image (see Tool versions).
  pbindex and minimap2 (HBA step) are called by full path from the `ptcp` env
  (`${params.condaEnvs.ptcp}bin/`), because another env is active in those steps.
- **PTCP scripts from lrsCustomTools**, run as `python3 <release>/ptcp/pacbio/<script>`
  with the Python of the step's env (the PTCP image has them on PATH).
- **pbmm2 uses `params.genomeMmi`** (the germline pipeline's index of the same
  `GRCh38_no_alt_analysis_set` reference) instead of indexing the FASTA.
- **trgt plot failures are logged per locus** (stderr) instead of being
  silently skipped by the WDL.
- **ptcp-qc manifest** is built from the joined per-sample channel instead of
  tarring and walking sample folders. Keys and order match PTCP; Paraphase
  VCFs are omitted as in PTCP.
- **File names** follow the germline prefix convention.

# Differences from the germline pipeline

- Same read handling (separate HiFi and fail alignments, fail reads to TRGT
  only), plus a merged copy for ptcp-qc; no hiPhase, and the pbmm2 BAMs are
  published instead of hiPhase BAMs.
- Output is aggregated per run (pureTarget layout), not one folder per sample.
- `minHifiGB` (35 GB, WGS) is not used; PTCP has its own `ptcpMinGB`.
- Run summaries are written to `runInfo/` only, not to the shared
  `lrsDocuments/summaryData` or `labWareExport` folders, and their names carry
  `.PTCP.`.
- No symlink maintenance in `onComplete`.
- PTCP processes run with `bash -l -eo pipefail` instead of the slurm
  profile's plain `bash -l`. A plain `bash -l` script carries on after a failed
  command, and only the last command's exit status decides whether the task
  failed. `-u` is left out because it breaks on unset variables in
  `/etc/profile.d` under a login shell.
