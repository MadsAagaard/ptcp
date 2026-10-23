#!/usr/bin/env python3
"""
ptcp_build_reference.py
─────────────────────────────────────────────────────────────────────────────
Builds the PTCP target-share reference used by the reports: for every carrier
target, the share of the sample's on-target reads that falls on that target,
summarised over a set of reference (validation) samples.

A share is independent of sequencing depth, so a sample can be compared with
the reference whether it was sequenced alone, in a mixed batch or topped up
on a second run. Low share = that target captured worse than usual in this
sample, which absolute read thresholds alone do not show.

    share(target) = locus_results[target].reads.n / stats.on_target.total
    (both HiFi + fail reads, as ptcp-qc counts them on the merged BAM)

The output holds medians and ranges only. Sample names are printed to stderr
for review but are NEVER written to the reference file, so the file can be
versioned in lrsCustomTools resources/ptcp/.

Usage
  ptcp_build_reference.py --name kgvejle_ptcp_ref_v1 \
      --out ptcp_target_share_reference.v1.json \
      run1/QC/ptcpQC/qc.*.json run2/QC/ptcpQC/qc.*.json [--exclude NPN ...]

  qc.aggregate.json files among the inputs are skipped. Samples are excluded
  when listed with --exclude, below --min-on-target on-target reads, or
  missing any carrier target. A sample present twice (same name) is used once.
"""

import argparse
import json
import statistics
import sys
from datetime import datetime
from pathlib import Path

__version__ = "261016"

# ptcp-qc locus groups of the carrier panel (12 genes; F8 has two targets)
CARRIER_TARGETS = ["smn1", "gba", "cyp21", "hba", "hbb", "rpgr", "f8inv1", "f8inv22",
                   "FXS_FMR1", "FRDA_FXN", "AFF2", "ARX"]


def quantile(sorted_vals, q):
    """linear interpolation, as numpy's default"""
    if not sorted_vals:
        return None
    pos = (len(sorted_vals) - 1) * q
    lo, hi = int(pos), min(int(pos) + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jsons", nargs="+", help="ptcp-qc per-sample JSONs (qc.<npn>.json)")
    ap.add_argument("--name", required=True, help="reference name/version shown in the reports")
    ap.add_argument("--out", required=True)
    ap.add_argument("--exclude", nargs="*", default=[], help="sample names to leave out")
    ap.add_argument("--min-on-target", type=int, default=2000, help="minimum on-target reads (default 2000)")
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    shares = {t: [] for t in CARRIER_TARGETS}
    on_fracs = []
    used, excluded, seen, versions, beds = [], {}, set(), set(), set()
    for p in a.jsons:
        try:
            d = json.loads(Path(p).read_text())
        except Exception as e:  # noqa: BLE001
            print(f"WARNING: {p} not readable: {e}", file=sys.stderr)
            continue
        if d.get("report_type") != "sample":
            continue
        sid = d.get("sample_name")
        why = None
        st = (d.get("stats") or {}).get("on_target") or {}
        loci = d.get("locus_results") or {}
        if sid in seen:
            why = "duplicate"
        elif sid in a.exclude:
            why = "excluded on request"
        elif (st.get("total") or 0) < a.min_on_target:
            why = f"on-target reads {st.get('total')} < {a.min_on_target}"
        else:
            missing = [t for t in CARRIER_TARGETS if not ((loci.get(t) or {}).get("reads") or {}).get("n")]
            if missing:
                why = "no reads at " + ",".join(missing)
        if why:
            excluded[why] = excluded.get(why, 0) + 1
            print(f"  exclude  {sid}: {why}", file=sys.stderr)
            continue
        seen.add(sid)
        used.append(sid)
        versions.add(d.get("ptcp-qc_version"))
        beds.add(d.get("targets_bed"))
        if (d.get("stats") or {}).get("total_reads"):
            on_fracs.append(st["total"] / d["stats"]["total_reads"])
        row = []
        for t in CARRIER_TARGETS:
            s = loci[t]["reads"]["n"] / st["total"]
            shares[t].append(s)
            row.append(f"{t}={s:.4f}")
        print(f"  use      {sid}: on-target {st['total']}  " + " ".join(row), file=sys.stderr)

    if len(used) < 3:
        sys.exit(f"ERROR: only {len(used)} usable sample(s); a reference needs at least 3.")
    if len(versions) > 1 or len(beds) > 1:
        print(f"WARNING: mixed ptcp-qc versions {sorted(map(str, versions))} / BEDs {sorted(map(str, beds))}", file=sys.stderr)

    targets = {}
    for t, v in shares.items():
        s = sorted(v)
        targets[t] = {"median": round(statistics.median(s), 5), "q1": round(quantile(s, 0.25), 5),
                      "q3": round(quantile(s, 0.75), 5), "min": round(s[0], 5), "max": round(s[-1], 5), "n": len(s)}
    of = sorted(on_fracs)
    ref = {
        "reference_type": "ptcp_target_share",
        "name": a.name,
        "created": datetime.now().strftime("%Y-%m-%d"),
        "builder_version": __version__,
        "n_samples": len(used),
        "n_excluded": excluded,
        "ptcp_qc_versions": sorted(map(str, versions)),
        "targets_bed": sorted(map(str, beds)),
        "share_definition": "locus_results[target].reads.n / stats.on_target.total (HiFi + fail reads)",
        "targets": targets,
        "on_target_fraction": {"median": round(statistics.median(of), 4), "min": round(of[0], 4), "max": round(of[-1], 4),
                               "n": len(of)} if of else None,
    }
    Path(a.out).write_text(json.dumps(ref, indent=1))
    print(f"reference '{a.name}': {len(used)} samples -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
