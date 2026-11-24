#!/usr/bin/env python3
"""Render NanoVI's run report (--report; F60, D42) from a results directory.

Reads only what a run publishes (sample_outputs/, consolidated_output/, qc/fastplong/, pipeline_info/; see
report_model.py) and writes into --outdir: index.html and one page per section (an offline site: fonts, styles and
the Vega libraries in assets/), figures/ (SVG, PDF, 600 dpi PNG and the plotted values as CSV) and
nanovi_report.pdf. A Typst failure gives a warning and no PDF. Two renders of the same results with the same
--date give identical files.
"""

import argparse
import json
import platform
import shutil
import sys
import time
from datetime import date, datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

ASSETS = Path(__file__).resolve().parent.parent / "assets" / "report"
PACKAGES = ["pandas", "numpy", "matplotlib", "jinja2", "altair", "vl-convert-python", "pyyaml"]


def report_versions():
    """The report's own software, in software_versions.yml's form (tool -> version)."""
    out = {"python": platform.python_version()}
    for p in PACKAGES:
        try:
            out[p] = version(p)
        except PackageNotFoundError:
            pass
    from report_render import typst_version
    t = typst_version()
    if t:
        out["typst"] = t
    return out


def report_date(text):
    """The date the report carries, and SOURCE_DATE_EPOCH for Typst: --date, else today (UTC)."""
    d = date.fromisoformat(text) if text else datetime.now(timezone.utc).date()
    return d.isoformat(), int(datetime(d.year, d.month, d.day, tzinfo=timezone.utc).timestamp())


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", help="NanoVI results directory (the run's --output_dir)")
    p.add_argument("--outdir", help="Where the report goes (default: <results>/report)")
    p.add_argument("--assets", default=str(ASSETS), help="Report templates and references (default: %(default)s)")
    p.add_argument("--params", help="JSON of the run's parameters, listed in the methods (optional)")
    p.add_argument("--metadata", help="Sample metadata, TSV or CSV: sample_id and optionally group (optional)")
    p.add_argument("--date", help="Date the report carries, YYYY-MM-DD (default: today, UTC)")
    p.add_argument("--no-pdf", action="store_true", help="Write the HTML site and figures only")
    p.add_argument("--versions", metavar="PROCESS", help="Print a versions.yml block for PROCESS and exit")
    args = p.parse_args(argv)
    if not args.versions and not args.results:
        p.error("--results is required")
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.versions:
        print(f'"{args.versions}":')
        for k, v in report_versions().items():
            print(f'    {k}: "{v}"')
        return 0
    import report_figures
    import report_model
    import report_render
    t0 = time.time()
    day, epoch = report_date(args.date)
    outdir = Path(args.outdir or Path(args.results) / "report")
    outdir.mkdir(parents=True, exist_ok=True)
    params = json.loads(Path(args.params).read_text(encoding="utf-8")) if args.params else None
    run = report_model.load_results(args.results, args.metadata)
    data = report_model.compute(run)
    fig_dir = outdir / "figures"
    if fig_dir.exists():
        shutil.rmtree(fig_dir)
    records = report_figures.draw_all(data, fig_dir)
    model = report_model.build(run, data, records, args.assets, params, {"REPORT": report_versions()}, day)
    pdf = None if args.no_pdf else report_render.write_pdf(outdir, model, args.assets, fig_dir, epoch)
    if pdf is None and (outdir / report_render.PDF_NAME).exists():
        (outdir / report_render.PDF_NAME).unlink()
    model["run"]["pdf"] = report_render.PDF_NAME if pdf else ""
    report_render.write_site(outdir, model, args.assets)
    print(f"Report for {len(run['samples'])} samples in {outdir} ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
