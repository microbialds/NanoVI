"""The run report's two documents, from the report model of report_model.py (adapted from BugBuster's
bin/summary_pdf.py and bin/summary_site.py, gene2dis/BugBuster dev @ de57f34):

- write_pdf: Typst templates (assets/report/typst/) reading the model as data.json, compiled with the static Inter
  TTFs only (--font-path, --ignore-system-fonts) and no Typst packages, so it needs no network. A missing or failing
  typst gives a warning and no PDF, never a failed run.
- write_site: a multi-page HTML site from the Jinja templates in assets/report/html/, one page per section; tables
  sortable and filterable (assets/report/static/tables.js, no library); read fate and composition interactive
  (report_interactive.py) with the Vega libraries bundled once in assets/; every figure links its static SVG, PDF,
  PNG and CSV in figures/. Nothing is loaded from the network, so the site works from file://.
"""

import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import figure_style as fs
import report_interactive
from report_model import die, warn

FONT_FACES = ("Regular", "Bold", "Italic", "BoldItalic", "Medium", "SemiBold")
WOFF2_FACES = ("Regular", "Bold", "Italic", "SemiBold")
PDF_NAME = "nanovi_report.pdf"


def pdf_view(model):
    """The model without what only the HTML reads: the figures' interactive data and the tables' HTML rows."""
    secs = []
    for sec in model["sections"]:
        figs = {k: {a: v for a, v in f.items() if a != "interactive"} for k, f in sec["figures"].items()}
        tabs = {k: ({a: v for a, v in t.items() if a != "rows"} if t else t) for k, t in sec["tables"].items()}
        secs.append({**sec, "figures": figs, "tables": tabs})
    return {**model, "sections": secs}


def typst_version():
    exe = shutil.which("typst")
    if exe is None:
        return None
    res = subprocess.run([exe, "--version"], capture_output=True, text=True)
    return res.stdout.split()[1] if res.returncode == 0 and len(res.stdout.split()) > 1 else None


def write_pdf(outdir, model, assets, fig_dir, epoch=None):
    """Compile <outdir>/nanovi_report.pdf; returns its path, or None (with a warning)."""
    typst = shutil.which("typst")
    if typst is None:
        warn(f"typst is not installed: no {PDF_NAME} (the pipeline's report image has it)")
        return None
    build = Path(outdir) / ".pdf_build"
    if build.exists():
        shutil.rmtree(build)
    shutil.copytree(Path(assets) / "typst", build)
    (build / "fonts").mkdir()
    fonts = fs.font_files(FONT_FACES)
    if not fonts:
        warn(f"no Inter TTF in {fs.font_dir()}: the PDF falls back to Typst's own fonts")
    for f in fonts:  # static TTFs only: with InterVariable in the folder Typst picks it
        shutil.copy(f, build / "fonts" / f.name)
    (build / "figures").mkdir()
    for sec in model["sections"]:
        for f in sec["figures"].values():
            for p in f["parts"]:
                shutil.copy(Path(fig_dir) / f"{p['key']}.svg", build / "figures" / f"{p['key']}.svg")
    methods = next((s for s in model["sections"] if s["id"] == "methods"), None)
    (build / "references.bib").write_text(methods["facts"].get("bibtex", "") if methods else "", encoding="utf-8")
    (build / "data.json").write_text(json.dumps(pdf_view(model), default=str, sort_keys=True), encoding="utf-8")
    out = Path(outdir).resolve() / PDF_NAME
    t0 = time.time()
    env = dict(os.environ)
    if epoch is not None:
        env["SOURCE_DATE_EPOCH"] = str(epoch)
    elif not env.get("SOURCE_DATE_EPOCH", "").strip().isdigit():
        env.pop("SOURCE_DATE_EPOCH", None)  # typst rejects an empty or malformed value
    try:
        res = subprocess.run([typst, "compile", "--ignore-system-fonts", "--font-path", "fonts", "report.typ",
                              str(out)], cwd=build, capture_output=True, text=True, env=env)
    finally:
        shutil.rmtree(build, ignore_errors=True)
    if res.returncode != 0:
        warn(f"typst failed, no {PDF_NAME}:\n{res.stderr.strip()}")
        if out.exists():
            out.unlink()
        return None
    for line in res.stderr.strip().splitlines():
        if line.strip():
            warn(f"typst: {line}")
    print(f"{PDF_NAME} written ({time.time() - t0:.1f} s)")
    return out


def write_site(outdir, model, assets):
    try:
        from jinja2 import Environment, FileSystemLoader, StrictUndefined, select_autoescape
    except ImportError:
        die("the HTML report needs jinja2 (conda-forge jinja2; the pipeline's report image has it)")
    t0 = time.time()
    outdir = Path(outdir)
    site_assets = outdir / "assets"
    if site_assets.exists():
        shutil.rmtree(site_assets)
    (site_assets / "fonts").mkdir(parents=True)
    for f in sorted((Path(assets) / "static").iterdir()):
        shutil.copy(f, site_assets / f.name)
    for f in fs.font_files(WOFF2_FACES, "woff2"):
        shutil.copy(f, site_assets / "fonts" / f.name)
    interactive = report_interactive.available()
    if not interactive:
        warn("altair / vl-convert not installed: the HTML report shows static figures only")
    specs = {}
    for sec in model["sections"]:
        for key, f in sec["figures"].items():
            specs[key] = report_interactive.spec_json(key, f.get("interactive")) if interactive else None
    bundle = any(specs.values())
    if bundle:
        (site_assets / "vega-bundle.js").write_text(report_interactive.javascript_bundle(), encoding="utf-8")
    env = Environment(loader=FileSystemLoader(str(Path(assets) / "html")), autoescape=select_autoescape(["html"]),
                      undefined=StrictUndefined, trim_blocks=True, lstrip_blocks=True, keep_trailing_newline=True)
    nav = [(s["page"], s["title"]) for s in model["sections"]]
    for sec in model["sections"]:
        html = env.get_template(f"{sec['id']}.html").render(run=model["run"], s=sec, nav=nav, specs=specs,
                                                            bundle=bundle, sections=model["sections"])
        (outdir / sec["page"]).write_text(html, encoding="utf-8")
    print(f"HTML report written: {len(model['sections'])} pages ({time.time() - t0:.1f} s)")
