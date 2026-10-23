#!/usr/bin/env python3
"""
make_mock_run.py — synthetic PTCP run for testing the report templates.

Takes ONE real ptcp-qc sample JSON and its havanno JSON as the shape donor and
writes eight synthetic samples (DEMO_S01..S08) into <out>/, laid out like a
run folder (QC/ptcpQC/, specialAnalysis/paraphase/) with an aggregate JSON,
a samples.tsv and a run_meta.json. Every call is overwritten, so the donor's
results do not survive; its read statistics do (scaled), which keeps every
field ptcp-qc writes present in the mock.

    python3 make_mock_run.py --qc qc.<npn>.json --havanno <npn>...havanno.json --out mockrun/

States exercised:
  S01 K  normal                        S06 K  CYP21A2 variant on a gene hap + HBA SEA (Sawfish-adjusted)
  S02 M  SMA carrier, single X loci    S07 K  sheet female, X/Y ratios male (sex mismatch)   [cell B]
  S03 K  HBA -a3.7 (Paraphase+Sawfish) S08 M  low yield: SMN no call, targets < 50 reads     [cell B]
  S04 K  FMR1 premutation (1 AGG); no SMN depth-model result although a model is set
  S05 K  F8 intron 22 inversion carrier
  S09 K  identical SMN copies: Paraphase 1:1, depth model 2:2 (used); topped up on cells A + B
  S10 K  SMN depth model disagrees with distinct haplotypes (review); GBA share 0.3x
  Cell A: S01-S06, S09 (part), S10. S11 is in the sheet without uBAMs.
"""
import argparse
import copy
import json
import statistics
from pathlib import Path


def load(p):
    return json.loads(Path(p).read_text())


def scale_reads(qc, f):
    st = qc["stats"]
    st["total_reads"] = round(st["total_reads"] * f)
    for k in ("on_target", "off_target", "unmapped"):
        for kk in st[k]:
            st[k][kk] = round(st[k][kk] * f)
    for loc in qc["locus_results"].values():
        for r in [loc.get("reads")] + [s.get("reads") for s in (loc.get("subgroups") or {}).values()]:
            if r:
                for kk in ("n", "hifi", "fail"):
                    r[kk] = round(r[kk] * f)
        for pr in (loc.get("paraphase_results") or {}).values():
            if not isinstance(pr, dict):
                continue
            for kk in ("total_reads", "unique_reads"):
                if kk in pr:
                    pr[kk] = round(pr[kk] * f)
            for h in (pr.get("haplotypes") or {}).values():
                for kk in ("total_count", "unique_count"):
                    h[kk] = round(h[kk] * f)
                h["fractional_count"] = round(h["fractional_count"] * f, 2)
            fi = pr.get("f8_info")
            if fi:
                for kk in fi["read_counts"]:
                    fi["read_counts"][kk] = round(fi["read_counts"][kk] * f)
        for t in (loc.get("trgt_results") or {}).values():
            t["total_reads"] = round(t["total_reads"] * f)
            for a in t["alleles"].values():
                a["coverage"] = round(a["coverage"] * f)


def set_smn(qc, s1, s2, r1, r2, failed=False):
    pr = qc["locus_results"]["smn1"]["paraphase_results"]["smn1"]
    hs = list(k for k in pr["haplotypes"] if k != "Unassigned")
    if failed:
        pr["smn_info"].update(smn1_cn=None, smn2_cn=None, read_numbers={"smn1": None, "smn2": None, "smn2_del78": None})
        return
    # rename haplotypes so SMN1 / SMN2 copies exist as named haps
    new = {"Unassigned": pr["haplotypes"].get("Unassigned", {"total_count": 0, "unique_count": 0, "fractional_count": 0})}
    smn1_haps, smn2_haps = [], []
    for i in range(s1):
        n = f"smn1_smn1hap{i + 1}"; smn1_haps.append(n)
        new[n] = {"total_count": r1 // max(s1, 1), "unique_count": r1 // max(s1, 1), "fractional_count": r1 // max(s1, 1)}
    for i in range(s2):
        n = f"smn1_smn2hap{i + 1}"; smn2_haps.append(n)
        new[n] = {"total_count": r2 // max(s2, 1), "unique_count": r2 // max(s2, 1), "fractional_count": r2 // max(s2, 1)}
    pr["haplotypes"] = new
    pr["smn_info"] = {"smn1_cn": s1, "smn2_cn": s2, "smn2_del78_cn": 0,
                      "read_numbers": {"smn1": r1, "smn2": r2, "smn2_del78": 0},
                      "haplotypes": {"smn1": smn1_haps, "smn2": smn2_haps, "smn2_del78": []}}
    return smn1_haps + smn2_haps


def set_ha(qc, possible, adj1, adj2, pred, model_cn=None, status="has_all_loci_for_depth_adjustment"):
    """smn_homology's homology_adjustment block, as test-batch writes it into the JSON."""
    prr = qc["locus_results"]["smn1"]["paraphase_results"]
    if status != "has_all_loci_for_depth_adjustment":
        prr["homology_adjustment"] = {"status": status}
        return
    pr, si = prr["smn1"], prr["smn1"]["smn_info"]
    ha = {"status": status, "smn1_unique_haplotypes": len(si["haplotypes"]["smn1"]),
          "smn2_unique_haplotypes": len(si["haplotypes"]["smn2"]), "is_possible_homology": possible,
          "predicted_coverage_per_copy": pred}
    for h, v in pr["haplotypes"].items():
        if h == "Unassigned":
            continue
        cn = v.get("n_copy", 1)
        per = v["fractional_count"] / cn
        mcn = (model_cn or {}).get(h, cn)
        ha[h] = {"paraphase_copy_number": cn, "predicted_copy_number": mcn, "prediction_confidence": 0.97,
                 "observed_coverage": int(per * cn), "observed_coverage_per_copy": int(per), "obs_exp_ratio": round(per / pred, 4)}
    ha.update(success=True, adjusted_smn1_cn=adj1, adjusted_smn2_cn=adj2,
              matches_paraphase=(adj1 == si["smn1_cn"] and adj2 == si["smn2_cn"]))
    prr["homology_adjustment"] = ha


def male(qc):
    qc["stats"]["sex_ratios"].update(ratio_x_auto=0.51, ratio_y_auto=0.30)
    loci = qc["locus_results"]
    for trid_group, trid in (("FXS_FMR1", "FXS_FMR1"), ("AFF2", "FRAXE_AFF2")):
        t = loci[trid_group]["trgt_results"][trid]
        t["alleles"]["1"]["coverage"] = 0          # single X allele
    for g in ("f8inv1", "f8inv22"):
        loci[g]["paraphase_results"][g]["f8_info"]["inversion_genotype"] = "0"
    rp = loci["rpgr"]["paraphase_results"]["rpgr"]
    rp["total_cn"] = 1
    rp["haplotypes"].pop("rpgr_hap2", None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    ap.add_argument("--havanno", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    base_qc, base_hv = load(a.qc), load(a.havanno)
    out = Path(a.out)
    (out / "QC/ptcpQC").mkdir(parents=True, exist_ok=True)
    (out / "specialAnalysis/paraphase").mkdir(parents=True, exist_ok=True)

    samples, cells = [], {}
    CELL_A, CELL_B = ("m84328_261001_120000_s1", "r84328_20261001_DEMO_A"), ("m84328_261008_120000_s2", "r84328_20261008_DEMO_B")
    def new(sid, sex, f, cell=CELL_A):
        cells[sid] = cell if isinstance(cell, list) else [cell]
        qc, hv = copy.deepcopy(base_qc), copy.deepcopy(base_hv)
        qc["sample_name"] = sid
        hv["input_paraphase_dir"] = f"{sid}.hg38v4LRS.HifiReads.paraphase"
        scale_reads(qc, f)
        haps = set_smn(qc, 2, 2, round(170 * f), round(175 * f))
        hv["annotations"]["smn1"] = {h: {"num_pathogenic_variants": 0, "num_pseudogene_variants": 0, "total_insertion_size": 0,
                                         "total_deletion_size": 0, "variants": {"pathogenic": [], "pseudogene": []}} for h in haps}
        if sex == "M":
            male(qc)
            hv["annotations"]["rpgr"].pop("rpgr_hap2", None)
        samples.append((sid, sex, qc, hv))
        return qc, hv

    def smn_annot(hv, haps):
        hv["annotations"]["smn1"] = {h: {"num_pathogenic_variants": 0, "num_pseudogene_variants": 0, "total_insertion_size": 0,
                                         "total_deletion_size": 0, "variants": {"pathogenic": [], "pseudogene": []}} for h in haps}

    qc, hv = new("DEMO_S01", "K", 1.00)
    set_ha(qc, False, 2, 2, 86.0)                          # distinct haplotypes, model agrees
    qc, hv = new("DEMO_S02", "M", 0.93)
    smn_annot(hv, set_smn(qc, 1, 2, 81, 176))
    set_ha(qc, False, 1, 2, 84.0)                          # carrier, model agrees
    qc, hv = new("DEMO_S03", "K", 1.08)
    hba = qc["locus_results"]["hba"]
    hba["paraphase_results"]["hba"].update(total_cn=3, genotype="-a/aa", genotype_adjusted="-a/aa",
                                          sv_calls={"hba_3p7delhap1": "3p7del,chr16:172240-176044"})
    hba["sawfish_results"] = {"hba": {"sv_breakpoints": [{"chrom": "chr16", "start": 172240, "end": 176044, "length": 3804, "svtype": "DEL",
                                                          "genotype": "0/1", "allele_depth": [188, 92], "annotation": {"name": "3p7del", "score": 1}}]}}
    qc, hv = new("DEMO_S04", "K", 1.04)
    fm = qc["locus_results"]["FXS_FMR1"]["trgt_results"]["FXS_FMR1"]["alleles"]["1"]
    fm.update(motif_count="71_1", size=219, ci={"min": 210, "max": 228}, coverage=168, methylation=0.05)
    qc, hv = new("DEMO_S05", "K", 0.97)
    qc["locus_results"]["f8inv22"]["paraphase_results"]["f8inv22"]["f8_info"].update(
        has_inversion=True, inversion_genotype="0/1", read_counts={"00": 120, "01": 96, "10": 81, "11": 230, "x": 50})
    qc, hv = new("DEMO_S06", "K", 1.11)
    h3 = hv["annotations"]["cyp21"]["cyp21_hap3"]
    h3["variants"]["pathogenic"] = ["CYP21A2_p.Ile173Asn_p.Ile173Asn_c.518T>A"]; h3["num_pathogenic_variants"] = 1
    hba = qc["locus_results"]["hba"]
    hba["paraphase_results"]["hba"].update(genotype="aa/aa", genotype_adjusted="--/aa")
    hba["sawfish_results"] = {"hba": {"sv_breakpoints": [{"chrom": "chr16", "start": 165396, "end": 184700, "length": 19304, "svtype": "DEL",
                                                          "genotype": "0/1", "allele_depth": [82, 186], "annotation": {"name": "SEA", "score": 1}}]}}
    qc, hv = new("DEMO_S07", "K", 0.95, CELL_B)
    qc["stats"]["sex_ratios"].update(ratio_x_auto=0.49, ratio_y_auto=0.31)
    qc, hv = new("DEMO_S08", "M", 0.08, CELL_B)
    set_smn(qc, 0, 0, 0, 0, failed=True)
    set_ha(qc, False, None, None, None, status="missing_loci_for_depth_adjustment")
    # identical SMN copies: Paraphase 1:1, depth model 2:2 -> call from the model; topped up on a second cell
    qc, hv = new("DEMO_S09", "K", 1.02, [CELL_A, CELL_B])
    smn_annot(hv, set_smn(qc, 1, 1, 176, 172))
    set_ha(qc, True, 2, 2, 87.0, {"smn1_smn1hap1": 2, "smn1_smn2hap1": 2})
    # distinct haplotypes, model disagrees (Paraphase 2:2, one SMN2 haplotype at twice the depth -> model 2:3)
    # -> review, call stays Paraphase; and a weak GBA capture
    qc, hv = new("DEMO_S10", "K", 1.00)
    h = qc["locus_results"]["smn1"]["paraphase_results"]["smn1"]["haplotypes"]["smn1_smn2hap2"]
    for k in ("total_count", "unique_count", "fractional_count"):
        h[k] = h[k] * 2
    set_ha(qc, False, 2, 3, 86.0, {"smn1_smn2hap2": 2})
    g = qc["locus_results"]["gba"]["reads"]
    for k in ("n", "hifi", "fail"):
        g[k] = round(g[k] * 0.3)

    # aggregate
    rows, loci = [], {}
    for sid, sex, qc, hv in samples:
        st = qc["stats"]
        rows.append({"sample_name": sid, "total_reads": st["total_reads"], "on_target_reads": st["on_target"]["total"],
                     "off_target_reads": st["off_target"]["total"], "unmapped_reads": st["unmapped"]["total"],
                     "duplicate_reads": 0, "duplicate_fraction": 0,
                     "on_target_fraction": round(st["on_target"]["total"] / st["total_reads"], 4),
                     "read_length_stats": {"median": 5066}, "read_quality_stats": {"median": 99.9949 if sid != "DEMO_S08" else 99.71},
                     "coverage": {"total_targets": 48}})
        for g, l in qc["locus_results"].items():
            loci.setdefault(g, []).append(l["reads"]["n"])
    agg = {"report_type": "aggregate", "genome_version": "hg38", "ptcp-qc_version": base_qc["ptcp-qc_version"],
           "summary": {"sample_count": len(samples), "total_reads": {"total": sum(r["total_reads"] for r in rows), "hifi": None, "fail": None}},
           "locus_stats": {g: {"samples_with_coverage": sum(1 for v in vs if v), "total_reads": sum(vs),
                               "coverage_stats": {"min": min(vs), "median": statistics.median(vs), "max": max(vs)}} for g, vs in loci.items()},
           "sample_stats": rows}
    (out / "QC/ptcpQC/qc.aggregate.json").write_text(json.dumps(agg, indent=1))
    for sid, sex, qc, hv in samples:
        (out / f"QC/ptcpQC/qc.{sid}.json").write_text(json.dumps(qc, indent=1))
        (out / f"specialAnalysis/paraphase/{sid}.hg38v4LRS.HifiReads.havanno.json").write_text(json.dumps(hv, indent=1))
    # analysed samples only, as the pipeline writes it; DEMO_S09 is in the sheet but had no uBAM (run_meta "missing")
    (out / "samples.tsv").write_text("id\tsex\ttestlist\tcells\truns\n" + "".join(
        f"{s[0]}\t{s[1]}\tPT-KENDT\t{','.join(c[0] for c in cells[s[0]])}\t{','.join(c[1] for c in cells[s[0]])}\n" for s in samples))
    meta = {"samplesheet": "DEMO_sheet", "runfolder": "r84328_DEMO", "read_mode": "HiFi + fail", "ptcp_version": "3.2.1",
            "genome_version_tag": "hg38v4LRS", "tag_hifi": "HifiReads", "tag_str": "AllReadsNew", "link_prefix": "../",
            "run_output": "pureTargetCarrier/DEMO_sheet_AllAndHifi", "samples_in_sheet": 11, "lrscustomtools_release": "releases/<tag>", "min_gb": 0,
            "missing": [{"id": "DEMO_S11", "testlist": "PT-KENDT"}], "dropped": [], "smn_model": "DEMO_smn_model.json",
            "tools": [{"tool": "pbmm2", "version": "1.16.99", "ptcp": "SMRT Link 25.1"}, {"tool": "TRGT", "version": "5.1.0", "ptcp": "5.0.0"},
                      {"tool": "Paraphase", "version": "3.5.0", "ptcp": "3.4.0"}, {"tool": "Sawfish", "version": "2.2.1", "ptcp": "2.2.1"},
                      {"tool": "ptcp-qc", "version": "1.0.0", "ptcp": "1.0.0"}],
            "demo_note": "SYNTHETIC run (DEMO_S01–S10, two SMRT Cells) built by make_mock_run.py to exercise every report state. Not patient data."}
    (out / "run_meta.json").write_text(json.dumps(meta, indent=1))
    print(f"mock run written to {out}")


if __name__ == "__main__":
    main()
