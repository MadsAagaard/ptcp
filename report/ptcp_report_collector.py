#!/usr/bin/env python3
"""
ptcp_report_collector.py
─────────────────────────────────────────────────────────────────────────────
Builds the PTCP (PureTarget Carrier Panel) report payloads from one run's
ptcp-qc and havanno output, and renders them into the two HTML templates:

  <npn>.PTCP.sampleReport.{json,html}   one per sample. Everything about THIS
                                        sample: its QC, Paraphase genes, F8,
                                        HBA, repeat expansions. The only run
                                        information in it is anonymous (target-
                                        share reference or same-cell median/IQR), so a
                                        reader of one sample's report never sees
                                        another sample's identity or results.
  <sheet>.PTCP.runQC.{json,html}        one per run. QC and callability for
                                        every sample, NO genotypes, copy numbers
                                        or repeat sizes.

Inputs (all from the run output folder of the KG Vejle PTCP pipeline):
  --qc-dir        QC/ptcpQC/ with qc.<npn>.json and qc.aggregate.json (ptcp-qc 1.0.0)
  --havanno       <npn>...havanno.json files (optional; matched on the leading NPN)
  --samples       TSV with header: id  sex  testlist  [cells  runs]  (sex M / K, from the sample
                  sheet; cells / runs = comma-separated SMRT Cell movies / run folders of its uBAMs)
  --reference     optional target-share reference (ptcp_build_reference.py); without it, target
                  shares are compared with the other samples of the same SMRT Cell (n >= --batch-min)
  --run-meta      JSON written by the pipeline (sheet, run folder, tools, missing/dropped, ...)
  --thresholds    optional YAML/JSON thresholds file, passed through verbatim

DESIGN RULES (as the germline / somatic / Kinnex collectors)
  1. The collector emits NUMBERS and names. Every threshold comparison, verdict,
     colour and sentence lives in the HTML templates.
  2. A missing input is a warning and an empty section, never a traceback.
  3. A sample report is built ONLY from that sample's files plus anonymous run
     medians; the run QC report never carries a call.

Field names follow ptcp-qc 1.0.0 as seen on real output (2026-10-15):
  - smn1 paraphase_results has no total_cn; the call is in smn_info
    (smn1_cn, smn2_cn, smn2_del78_cn — may be null — read_numbers, haplotypes).
  - Haplotype entries carry n_copy only when > 1; havanno repeats such copies
    as <hap>_cp2.
  - havanno wraps its loci in "annotations".
  - trgt_results are keyed by TRID, which can differ from the locus group
    (group AFF2 -> FRAXE_AFF2; group ARX -> EIEE1_ARX and PRTS_ARX).
  - sawfish_results is absent when Sawfish called nothing.
"""

import argparse
import json
import re
import statistics
import sys
from datetime import datetime
from pathlib import Path

__version__ = "261016"

PARAPHASE_GROUPS = ["smn1", "gba", "cyp21", "hba", "hbb", "rpgr", "f8inv1", "f8inv22"]


def warn(msg):
    print(f"WARNING: {msg}", file=sys.stderr)


def load_json(path, what):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception as e:  # noqa: BLE001 — any unreadable input becomes an empty section
        warn(f"{what} not readable ({path}): {e}")
        return None


def load_thresholds(path):
    if not path:
        return {}
    text = Path(path).read_text()
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        try:
            import yaml  # optional; only needed for a YAML thresholds file
        except ImportError:
            warn(f"thresholds file {path} is not JSON and PyYAML is not installed — ignored")
            return {}
        data = yaml.safe_load(text)
    if isinstance(data, dict):
        data["source_file"] = str(path)
        return data
    warn(f"thresholds file {path} is not a mapping — ignored")
    return {}


def stat_median(block):
    return block.get("median") if isinstance(block, dict) else None


def reads_summary(reads):
    """locus_results[*].reads -> counts + medians (no distributions)."""
    if not isinstance(reads, dict):
        return None
    st = reads.get("stats") or {}
    return {
        "n": reads.get("n"),
        "hifi": reads.get("hifi"),
        "fail": reads.get("fail"),
        "read_quality_median": stat_median(st.get("read_quality")),
        "read_passes_median": stat_median(st.get("read_passes")),
        "read_length_median": stat_median(st.get("read_length")),
    }


# ── Paraphase loci ────────────────────────────────────────────────────────────

def havanno_entry(a):
    v = a.get("variants") or {}
    return {
        "num_pathogenic": a.get("num_pathogenic_variants"),
        "num_pseudogene": a.get("num_pseudogene_variants"),
        "insertion_size": a.get("total_insertion_size"),
        "deletion_size": a.get("total_deletion_size"),
        "pathogenic": v.get("pathogenic") or [],
        "pseudogene": v.get("pseudogene") or [],
    }


def paraphase_locus(group, locus, havanno_locus):
    pr_all = (locus or {}).get("paraphase_results") or {}
    pr = pr_all.get(group)
    if not isinstance(pr, dict):
        return {"present": False}

    haps, unassigned = [], None
    for name, h in (pr.get("haplotypes") or {}).items():
        if name == "Unassigned":
            unassigned = h.get("total_count")
            continue
        haps.append({
            "name": name,
            "reads_total": h.get("total_count"),
            "reads_unique": h.get("unique_count"),
            "reads_fractional": h.get("fractional_count"),
            "n_copy": h.get("n_copy", 1),
            "havanno": None,
        })

    unmatched = []
    if havanno_locus is not None:
        names = {h["name"] for h in haps}
        for hname, a in havanno_locus.items():
            base = re.sub(r"_cp\d+$", "", hname)
            if hname in names:
                next(h for h in haps if h["name"] == hname)["havanno"] = havanno_entry(a)
            elif base not in names:
                unmatched.append(hname)          # never silently dropped
        for h in haps:
            if h["havanno"] is None:
                unmatched.append(f"{h['name']} (no havanno entry)")

    out = {
        "present": True,
        "total_cn": pr.get("total_cn"),
        "total_reads": pr.get("total_reads"),
        "unique_reads": pr.get("unique_reads"),
        "unassigned_reads": unassigned,
        "haplotypes": haps,
        "havanno_present": havanno_locus is not None,
        "havanno_unmatched": unmatched,
        "subgroups": {k: reads_summary(v.get("reads")) for k, v in ((locus or {}).get("subgroups") or {}).items()},
    }
    for key in ("genotype", "genotype_adjusted", "sv_calls", "f8_info", "smn_info"):
        if key in pr:
            out[key] = pr[key]
    if "homology_adjustment" in pr_all:
        out["homology_adjustment"] = pr_all["homology_adjustment"]
    if group == "hba":
        sw = ((locus or {}).get("sawfish_results") or {}).get("hba") or {}
        out["sawfish"] = sw.get("sv_breakpoints") or []
    return out


# ── TRGT loci ─────────────────────────────────────────────────────────────────

def trgt_loci(locus_results):
    out = {}
    for group, locus in (locus_results or {}).items():
        for trid, t in ((locus or {}).get("trgt_results") or {}).items():
            alleles = []
            for key in sorted((t.get("alleles") or {}).keys(), key=lambda k: int(k) if str(k).isdigit() else 99):
                a = t["alleles"][key]
                alleles.append({
                    "allele": key,
                    "coverage": a.get("coverage"),
                    "size": a.get("size"),
                    "ci": a.get("ci"),
                    "motif_count": a.get("motif_count"),
                    "motif_spans": a.get("motif_spans"),
                    "purity": a.get("purity"),
                    "methylation": a.get("methylation"),     # absent when TRGT reports none
                })
            out[trid] = {"group": group, "repeat_unit": t.get("repeat_unit"),
                         "total_reads": t.get("total_reads"), "alleles": alleles}
    return out


# ── Per-sample payload ───────────────────────────────────────────────────────

def sample_files(sid, meta):
    """Relative links from the reports/ folder to this sample's outputs."""
    gv, th, ts = meta.get("genome_version_tag", "hg38v4LRS"), meta.get("tag_hifi", "HifiReads"), meta.get("tag_str", "AllReadsNew")
    p = meta.get("link_prefix", "../")
    h, a = f"{sid}.{gv}.{th}", f"{sid}.{gv}.{ts}"
    return {
        "hifi_bam":         f"{p}alignments/{h}.pbmm2.merged.bam",
        "paraphase_json":   f"{p}specialAnalysis/paraphase/{h}.paraphase/{sid}.paraphase.json",
        "paraphase_bam":    f"{p}specialAnalysis/paraphase/{h}.paraphase/{sid}.paraphase.bam",
        "havanno_json":     f"{p}specialAnalysis/paraphase/{h}.havanno.json",
        "f8_json":          f"{p}specialAnalysis/paraphase/f8inversion/{h}.f8inversion.json",
        "f8_vcf":           f"{p}specialAnalysis/paraphase/f8inversion/{h}.f8inversion.vcf",
        "sawfish_vcf":      f"{p}structuralVariants/HBA/{h}.sawfishSV.HBA.vcf.gz",
        "hba_minimap2_bam": f"{p}structuralVariants/HBA/{h}.HBA.minimap2.bam",
        "trgt_vcf":         f"{p}repeatExpansions/TRGT/{a}.trgt5.ptcp.vcf.gz",
        "trgt_spanning_bam": f"{p}repeatExpansions/TRGT/{a}.trgt5.ptcp.sorted.spanning.bam",
        "trgt_plot":        f"{p}repeatExpansions/TRGT/Plots/{sid}/{a}.trgt5.ptcp.{{mode}}/{a}.trgt5.ptcp.{{mode}}.{{trid}}.trgt_plot.svg",
        "ptcp_qc_json":     f"{p}QC/ptcpQC/qc.{sid}.json",
    }


def build_sample(sid, sample_qc, agg_row, havanno, sheet_row, comp, run_meta, thresholds):
    stats = (sample_qc or {}).get("stats") or {}
    loci = (sample_qc or {}).get("locus_results") or {}
    ann = (havanno or {}).get("annotations") or {}
    qc = {
        "total_reads": stats.get("total_reads"),
        "on_target": stats.get("on_target"),
        "off_target": stats.get("off_target"),
        "unmapped": stats.get("unmapped"),
        "sex_ratios": stats.get("sex_ratios"),
        "on_target_fraction": (agg_row or {}).get("on_target_fraction"),
        "duplicate_fraction": (agg_row or {}).get("duplicate_fraction"),
        "read_length_median": stat_median((agg_row or {}).get("read_length_stats")),
        "read_quality_median": stat_median((agg_row or {}).get("read_quality_stats")),
    }
    return {
        "report_type": "ptcp_sample",
        "metadata": {
            "id": sid,
            "testlist": (sheet_row or {}).get("testlist"),
            "sex_reported": (sheet_row or {}).get("sex"),
            "samplesheet": run_meta.get("samplesheet"),
            "runfolder": run_meta.get("runfolder"),
            "read_mode": run_meta.get("read_mode"),
            "genome_version": (sample_qc or {}).get("genome_version"),
            "ptcp_version": run_meta.get("ptcp_version"),
            "ptcp_qc_version": (sample_qc or {}).get("ptcp-qc_version"),
            "tools": run_meta.get("tools", []),
            "lrscustomtools_release": run_meta.get("lrscustomtools_release"),
            "smn_model": run_meta.get("smn_model"),
            "cells": [c for c in ((sheet_row or {}).get("cells") or "").split(",") if c],
            "runs": [r for r in ((sheet_row or {}).get("runs") or "").split(",") if r],
            "havanno_variant_list": (havanno or {}).get("input_variant_list"),
            "collector_version": __version__,
            "report_date": datetime.now().strftime("%Y-%m-%d"),
            "run_qc_report": f"{run_meta.get('samplesheet', 'run')}.PTCP.runQC.html",
            "inputs_missing": [k for k, v in (("ptcp-qc sample JSON", sample_qc), ("havanno JSON", havanno),
                                              ("sample sheet row", sheet_row), ("aggregate sample_stats row", agg_row)) if v is None],
        },
        "qc": qc,
        "loci": {g: reads_summary(l.get("reads")) for g, l in loci.items()},
        "target_shares": target_shares(sample_qc),
        "comparison": comp,
        "paraphase": {g: paraphase_locus(g, loci.get(g), ann.get(g) if havanno else None) for g in PARAPHASE_GROUPS},
        "trgt": trgt_loci(loci),
        "files": sample_files(sid, run_meta),
        "thresholds": thresholds,
        "_demo": run_meta.get("demo_note"),
        "_sources": {
            "qc": f"QC/ptcpQC/qc.{sid}.json (stats, locus_results[*].reads) · qc.aggregate.json (sample_stats) · target-share comparison: "
                  + {"reference": f"reference {comp.get('label')}", "batch": f"SMRT Cell {comp.get('label')} (leave-one-out)"}.get(comp.get("source"), "none"),
            "paraphase": f"QC/ptcpQC/qc.{sid}.json paraphase_results · {Path(havanno.get('input_paraphase_dir', '')).name if havanno else 'no havanno JSON'}",
            "f8": f"QC/ptcpQC/qc.{sid}.json paraphase_results.f8inv1/f8inv22.f8_info",
            "hba": f"QC/ptcpQC/qc.{sid}.json locus_results.hba (paraphase_results, sawfish_results)",
            "trgt": f"QC/ptcpQC/qc.{sid}.json trgt_results",
        },
    }


# ── Run payloads ─────────────────────────────────────────────────────────────

CARRIER_TARGETS = ["smn1", "gba", "cyp21", "hba", "hbb", "rpgr", "f8inv1", "f8inv22",
                   "FXS_FMR1", "FRDA_FXN", "AFF2", "ARX"]


def quantile(sorted_vals, q):
    """linear interpolation (numpy default); same as ptcp_build_reference.py"""
    if not sorted_vals:
        return None
    pos = (len(sorted_vals) - 1) * q
    lo, hi = int(pos), min(int(pos) + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def target_shares(sample_qc):
    """Share of the sample's on-target reads per carrier target (HiFi + fail),
    the same definition as ptcp_build_reference.py. Depth-independent."""
    st = ((sample_qc or {}).get("stats") or {}).get("on_target") or {}
    tot = st.get("total")
    if not tot:
        return None
    loci = (sample_qc or {}).get("locus_results") or {}
    return {t: ((loci.get(t) or {}).get("reads") or {}).get("n", 0) / tot for t in CARRIER_TARGETS}


def on_target_fraction(sample_qc):
    st = (sample_qc or {}).get("stats") or {}
    return (st.get("on_target") or {}).get("total", 0) / st["total_reads"] if st.get("total_reads") else None


def single_cell(row):
    """The sample's SMRT Cell when all its uBAMs come from one cell, else None."""
    cells = [c for c in ((row or {}).get("cells") or "").split(",") if c]
    return cells[0] if len(cells) == 1 else None


def comparison(sid, sample_qcs, sheet, reference, batch_min):
    """What this sample's target shares are compared with, anonymously:
       1. the lab reference (validation samples), when supplied;
       2. otherwise the other samples of the same SMRT Cell in this run
          (leave-one-out median), when there are at least batch_min of them;
       3. otherwise nothing, with the reason.
    Never carries another sample's name or individual values."""
    if reference:
        sample_ver = (sample_qcs.get(sid) or {}).get("ptcp-qc_version")
        return {
            "source": "reference",
            "label": reference.get("name"),
            "n": reference.get("n_samples"),
            "created": reference.get("created"),
            "targets": reference.get("targets") or {},
            "on_target_fraction": reference.get("on_target_fraction"),
            "version_mismatch": bool(sample_ver and reference.get("ptcp_qc_versions")
                                     and sample_ver not in reference["ptcp_qc_versions"]),
            "ptcp_qc_versions": reference.get("ptcp_qc_versions"),
        }
    cell = single_cell(sheet.get(sid))
    if not cell:
        return {"source": "none", "reason": "no reference supplied, and the sample is not from a single SMRT Cell"}
    others = [o for o in sheet if o != sid and single_cell(sheet.get(o)) == cell and sample_qcs.get(o)]
    if len(others) < batch_min:
        return {"source": "none", "label": cell, "n": len(others),
                "reason": f"no reference supplied, and only {len(others)} other sample(s) of this SMRT Cell in the run (minimum {batch_min})"}
    # Spread as the interquartile range, never min / max: an extreme is one
    # identifiable other sample's value, a quartile of >= 5 samples is not.
    sh = [target_shares(sample_qcs[o]) for o in others]
    sh = [s for s in sh if s]
    targets = {}
    for t in CARRIER_TARGETS:
        v = sorted(s[t] for s in sh)
        targets[t] = {"median": statistics.median(v), "min": quantile(v, 0.25), "max": quantile(v, 0.75), "n": len(v), "spread": "IQR"}
    of = sorted(x for x in (on_target_fraction(sample_qcs[o]) for o in others) if x is not None)
    return {"source": "batch", "label": cell, "n": len(sh), "targets": targets,
            "on_target_fraction": {"median": statistics.median(of), "min": quantile(of, 0.25), "max": quantile(of, 0.75),
                                   "n": len(of), "spread": "IQR"} if of else None}


def callability(sample_qc):
    """Per target: did the caller produce a call? Status only, never the call."""
    loci = (sample_qc or {}).get("locus_results") or {}
    out = {}
    for g in PARAPHASE_GROUPS:
        pr = ((loci.get(g) or {}).get("paraphase_results") or {}).get(g)
        if g == "smn1":
            si = (pr or {}).get("smn_info") or {}
            out[g] = si.get("smn1_cn") is not None
        elif g.startswith("f8inv"):
            out[g] = ((pr or {}).get("f8_info") or {}).get("has_inversion") is not None
        else:
            out[g] = pr is not None and pr.get("total_cn") is not None
    for trid, t in trgt_loci(loci).items():
        out[trid] = any((a.get("coverage") or 0) > 0 for a in t["alleles"])
    return out


def build_run(aggregate, sample_qcs, sheet, run_meta, thresholds, reference):
    agg_rows = {r.get("sample_name"): r for r in (aggregate or {}).get("sample_stats") or []}
    ids = list(dict.fromkeys(list(sheet.keys()) + list(agg_rows.keys())))
    samples = []
    for sid in ids:
        sq = sample_qcs.get(sid)
        st = (sq or {}).get("stats") or {}
        ar = agg_rows.get(sid) or {}
        samples.append({
            "id": sid,
            "testlist": (sheet.get(sid) or {}).get("testlist"),
            "sex_reported": (sheet.get(sid) or {}).get("sex"),
            "in_sample_sheet": sid in sheet,
            "in_aggregate": sid in agg_rows,
            "has_sample_json": sq is not None,
            "total_reads": ar.get("total_reads", st.get("total_reads")),
            "on_target": st.get("on_target"),
            "on_target_reads": ar.get("on_target_reads"),
            "on_target_fraction": ar.get("on_target_fraction"),
            "duplicate_fraction": ar.get("duplicate_fraction"),
            "read_length_median": stat_median(ar.get("read_length_stats")),
            "read_quality_median": stat_median(ar.get("read_quality_stats")),
            "sex_ratios": st.get("sex_ratios"),
            "loci": {g: (l.get("reads") or {}).get("n") for g, l in ((sq or {}).get("locus_results") or {}).items()},
            "callable": callability(sq) if sq else None,
            "target_shares": target_shares(sq),
            "cells": [c for c in ((sheet.get(sid) or {}).get("cells") or "").split(",") if c],
            "runs": [r for r in ((sheet.get(sid) or {}).get("runs") or "").split(",") if r],
        })
    return {
        "report_type": "ptcp_run_qc",
        "metadata": dict(run_meta, collector_version=__version__, report_date=datetime.now().strftime("%Y-%m-%d"),
                         ptcp_qc_version=(aggregate or {}).get("ptcp-qc_version")),
        "run": {
            # sheet rows = analysed + no uBAM found + dropped by the size gate
            "samples_in_sheet": run_meta.get("samples_in_sheet",
                                             len(sheet) + len(run_meta.get("missing", [])) + len(run_meta.get("dropped", []))),
            "ptcp_qc_samples": (aggregate or {}).get("summary", {}).get("sample_count"),
            "missing": run_meta.get("missing", []),
            "dropped": run_meta.get("dropped", []),
            "total_reads": (aggregate or {}).get("summary", {}).get("total_reads"),
        },
        "locus_stats": (aggregate or {}).get("locus_stats") or {},
        "samples": samples,
        "reference": {k: reference.get(k) for k in ("name", "created", "n_samples", "ptcp_qc_versions", "targets", "on_target_fraction")} if reference else None,
        "thresholds": thresholds,
        "_demo": run_meta.get("demo_note"),
        "_sources": {
            "run": "launch parameters · runInfo/<date>_<sheet>/*.PTCP.*.summary.txt",
            "qc": "QC/ptcpQC/qc.aggregate.json (sample_stats, locus_stats) · qc.<npn>.json (stats, locus reads)",
        },
    }


# ── HTML rendering (same mechanism as the germline collector) ────────────────

def render_html(payload, template_path, out_path):
    if not template_path or not Path(template_path).exists():
        warn(f"HTML template not found: {template_path} — HTML not written")
        return
    template = Path(template_path).read_text(encoding="utf-8")
    data = json.dumps(payload, indent=1, ensure_ascii=False)
    data = data.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    sentinel = "/* DATA_INJECT */"
    updated = re.sub(r"/\* DATA_INJECT \*/\s*const DATA = [^\n]*;",
                     lambda _: f"{sentinel}\nconst DATA = {data};", template, count=1)
    if updated == template:
        warn(f"injection point not found in {template_path} — HTML not written")
        return
    Path(out_path).write_text(updated, encoding="utf-8")
    print(f"  {out_path}")


def write_json(payload, path):
    Path(path).write_text(json.dumps(payload, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"  {path}")


def read_sheet(path):
    rows = {}
    if not path:
        return rows
    lines = [l.rstrip("\n") for l in Path(path).read_text().splitlines() if l.strip()]
    header = lines[0].split("\t")
    for line in lines[1:]:
        r = dict(zip(header, line.split("\t")))
        rows[r["id"]] = r
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-dir", required=True, help="folder with qc.<npn>.json and qc.aggregate.json")
    ap.add_argument("--havanno", nargs="*", default=[], help="havanno JSONs, file names starting with the NPN")
    ap.add_argument("--samples", required=True, help="TSV: id, sex (M/K), testlist")
    ap.add_argument("--run-meta", required=True, help="run metadata JSON from the pipeline")
    ap.add_argument("--thresholds", default=None, help="thresholds YAML/JSON (optional)")
    ap.add_argument("--reference", default=None, help="target-share reference JSON from ptcp_build_reference.py (optional)")
    ap.add_argument("--batch-min", type=int, default=5,
                    help="without a reference: minimum number of OTHER samples of the same SMRT Cell for a batch comparison (default 5)")
    ap.add_argument("--sample-template", required=True)
    ap.add_argument("--run-template", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args()

    qc_dir, outdir = Path(args.qc_dir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_meta = load_json(args.run_meta, "run metadata") or {}
    thresholds = load_thresholds(args.thresholds)
    sheet = read_sheet(args.samples)
    aggregate = load_json(qc_dir / "qc.aggregate.json", "ptcp-qc aggregate")
    agg_rows = {r.get("sample_name"): r for r in (aggregate or {}).get("sample_stats") or []}

    sample_qcs = {}
    for sid in sheet:
        p = qc_dir / f"qc.{sid}.json"
        sample_qcs[sid] = load_json(p, f"ptcp-qc sample JSON for {sid}") if p.exists() else None
        if sample_qcs[sid] is None:
            warn(f"no ptcp-qc sample JSON for {sid}")

    havanno = {}
    for p in args.havanno:
        sid = Path(p).name.split(".")[0]
        havanno[sid] = load_json(p, f"havanno JSON {p}")

    reference = load_json(args.reference, "target-share reference") if args.reference else None
    if reference and reference.get("reference_type") != "ptcp_target_share":
        warn(f"{args.reference} is not a ptcp_target_share reference — ignored")
        reference = None
    sheet_name = run_meta.get("samplesheet", "run")
    print("Sample reports:")
    for sid in sheet:
        payload = build_sample(sid, sample_qcs.get(sid), agg_rows.get(sid), havanno.get(sid),
                               sheet.get(sid), comparison(sid, sample_qcs, sheet, reference, args.batch_min),
                               run_meta, thresholds)
        write_json(payload, outdir / f"{sid}.PTCP.sampleReport.json")
        render_html(payload, args.sample_template, outdir / f"{sid}.PTCP.sampleReport.html")

    print("Run QC report:")
    run_payload = build_run(aggregate, sample_qcs, sheet, run_meta, thresholds, reference)
    write_json(run_payload, outdir / f"{sheet_name}.PTCP.runQC.json")
    render_html(run_payload, args.run_template, outdir / f"{sheet_name}.PTCP.runQC.html")


if __name__ == "__main__":
    main()
