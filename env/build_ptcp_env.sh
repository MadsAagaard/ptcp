#!/usr/bin/env bash
#
# Build the 'ptcp' conda env used by the KG Vejle PTCP pipeline.
#
#   1. conda env from env/ptcp_env.yml (Python packages for the PTCP scripts,
#      pbtk for pbindex, minimap2)
#   2. ptcp-qc release binary (MD5 checked against the release's .md5 file)
#   3. a version record in $ENV/ptcp_versions.txt, for documentation
#
# The PTCP scripts themselves (f8_inversion.py, havanno, smn_homology,
# ptcpqc_json_to_csv.py) and the PTCP reference files are NOT installed here:
# they come from the lrsCustomTools release (ptcp/pacbio/, resources/ptcp/).
#
# Usage:
#   bash env/build_ptcp_env.sh                       # default location below
#   ENV=/some/other/path bash env/build_ptcp_env.sh
#
# Needs network access to the conda channels and github.com.
# Re-running on an existing env is refused; remove it first if you mean to rebuild.

set -euo pipefail

ENV="${ENV:-/lnx01_data3/shared/programmer/miniconda3/envs/ptcp}"
PTCPQC_VERSION="${PTCPQC_VERSION:-1.0.0}"
CONDA="${CONDA:-$(command -v mamba || command -v conda)}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -e "$ENV" ]; then
    echo "ERROR: $ENV already exists. Remove it first to rebuild." >&2
    exit 1
fi

echo ">> Creating conda env at $ENV (using $CONDA)"
"$CONDA" env create -p "$ENV" -f "$HERE/ptcp_env.yml"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo ">> ptcp-qc $PTCPQC_VERSION"
QC="ptcp-qc-v${PTCPQC_VERSION}-x86_64-unknown-linux-gnu"
curl -sSfL -o "$TMP/ptcpqc.tar.gz" \
    "https://github.com/PacificBiosciences/ptcp-qc/releases/download/${PTCPQC_VERSION}/${QC}.tar.gz"
tar -xzf "$TMP/ptcpqc.tar.gz" -C "$TMP"
( cd "$TMP/$QC" && md5sum -c ptcp-qc.md5 )
install -m 755 "$TMP/$QC/ptcp-qc" "$ENV/bin/ptcp-qc"

echo ">> Smoke test"
# the packages the PTCP scripts import (f8_inversion.py: pysam;
# smn_homology: numpy, pandas, scipy, scikit-learn)
"$ENV/bin/python3" -c "import pysam, numpy, pandas, scipy, sklearn"
"$ENV/bin/ptcp-qc" --version
"$ENV/bin/pbindex" --version | head -1
"$ENV/bin/minimap2" --version

echo ">> Writing $ENV/ptcp_versions.txt"
{
    echo "built:            $(date -Iseconds) by ${USER:-unknown}"
    echo "ptcp-qc:          $("$ENV/bin/ptcp-qc" --version)"
    echo
    "$CONDA" list -p "$ENV" --explicit --md5
} > "$ENV/ptcp_versions.txt"

echo ">> Done: $ENV"
