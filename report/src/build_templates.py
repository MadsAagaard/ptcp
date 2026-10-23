#!/usr/bin/env python3
"""Assemble the two self-contained PTCP report templates from src/.

  src/style.css + src/helpers.js + src/scrollspy.js are shared; the bodies are
  src/sample_body.html and src/run_body.html. Run from anywhere:

      python3 report/src/build_templates.py

The assembled templates are what the collector and the pipeline use.
"""
from pathlib import Path

SRC = Path(__file__).resolve().parent
OUT = SRC.parent
parts = {
    "<!--STYLE-->": (SRC / "style.css").read_text(),
    "<!--HELPERS-->": (SRC / "helpers.js").read_text(),
    "<!--SCROLLSPY-->": (SRC / "scrollspy.js").read_text(),
}
for body, name in (("sample_body.html", "ptcp_sample_report_template.html"),
                   ("run_body.html", "ptcp_run_qc_report_template.html")):
    text = (SRC / body).read_text()
    for marker, content in parts.items():
        assert text.count(marker) == 1, f"{body}: {marker} must appear once"
        text = text.replace(marker, content.rstrip("\n"))
    (OUT / name).write_text(text)
    print(OUT / name)
