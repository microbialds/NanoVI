"""Rendering the run report (bin/nanovi_report.py; F60, D42): the HTML site and figures everywhere the Python
packages are installed, the PDF where typst is (the report image, conda's environment). Two renders of the same
results give identical files, and no page makes a network request."""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip('matplotlib')
pytest.importorskip('jinja2')
pytest.importorskip('yaml')

from .report_fixture import make_results

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
TYPST = shutil.which('typst') is not None
# Attributes and CSS that make a browser fetch something; a Vega-Lite spec's "$schema" URL is data, not a request
REQUEST = re.compile(r'''(?:\b(?:src|href|action|poster|data)\s*=\s*["']?\s*(?:https?:)?//)|url\(\s*["']?\s*(?:https?:)?//|@import''',
                     re.I)


def render(results, outdir):
    res = subprocess.run([sys.executable, os.path.join(BIN, 'nanovi_report.py'), '--results', str(results),
                          '--outdir', str(outdir), '--date', '2026-11-06'] + ([] if TYPST else ['--no-pdf']),
                         capture_output=True, text=True, env={**os.environ, 'MPLCONFIGDIR': str(outdir) + '.mpl'})
    assert res.returncode == 0, res.stderr
    return res


def files(root):
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(Path(root).rglob('*')) if p.is_file()}


@pytest.fixture(scope='module')
def rendered(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('report')
    results = make_results(tmp / 'results', groups=True, cell=True, reads=True)
    render(results, tmp / 'a')
    render(results, tmp / 'b')
    return tmp


def test_two_renders_are_byte_identical(rendered):
    a, b = files(rendered / 'a'), files(rendered / 'b')
    assert a.keys() == b.keys()
    assert [k for k in a if a[k] != b[k]] == []


def test_site_has_a_page_per_section_and_every_figure_file(rendered):
    out = rendered / 'a'
    for page in ('index', 'composition', 'flagged', 'groups', 'cell', 'reads', 'methods'):
        assert (out / f'{page}.html').exists()
    for key in ('read_fate', 'composition_species', 'composition_genus', 'flagged_species', 'groups_heatmap',
                'cell_vs_read'):
        for ext in ('svg', 'pdf', 'png', 'csv'):
            assert (out / 'figures' / f'{key}.{ext}').stat().st_size > 0, (key, ext)
    assert not (out / 'figures' / 'heatmap_species.svg').exists()  # three samples: no heatmap


def test_no_page_requests_the_network(rendered):
    out = rendered / 'a'
    for p in list(out.glob('*.html')) + [out / 'assets' / 'site.css', out / 'assets' / 'tables.js']:
        text = p.read_text(encoding='utf-8')
        assert not REQUEST.search(text), (p.name, REQUEST.search(text).group(0))
    for p in out.glob('*.html'):
        for ref in re.findall(r'''\b(?:src|href)="([^"]+)"''', p.read_text(encoding='utf-8')):
            if not ref.startswith('#'):
                assert (out / ref).exists(), (p.name, ref)


def test_interactive_figures_are_inline_with_the_bundled_libraries(rendered):
    out = rendered / 'a'
    pytest.importorskip('altair')
    pytest.importorskip('vl_convert')
    assert (out / 'assets' / 'vega-bundle.js').stat().st_size > 100_000
    index = (out / 'index.html').read_text(encoding='utf-8')
    assert 'vegaEmbed("#vis-read_fate"' in index and 'src="assets/vega-bundle.js"' in index
    comp = (out / 'composition.html').read_text(encoding='utf-8')
    assert comp.count('vegaEmbed(') == 2  # species and genus bars; the flagged and cell figures stay static


def test_figure_csv_holds_the_plotted_values(rendered):
    import pandas as pd
    fate = pd.read_csv(rendered / 'a' / 'figures' / 'read_fate.csv')
    s1 = fate[(fate['sample'] == 'S1') & (fate['fate'] == 'reported species')]
    assert s1['pct_of_passed_reads'].iloc[0] == pytest.approx(100 * 100 / 110)
    comp = pd.read_csv(rendered / 'a' / 'figures' / 'composition_species.csv')
    assert comp.groupby('sample')['abundance'].sum().round(12).to_dict() == {'S1': 1.0, 'S10': 1.0, 'S2': 0.0}


@pytest.mark.skipif(not TYPST, reason='typst is not installed (the report image has it)')
def test_pdf_is_written_with_inter_embedded(rendered):
    pdf = rendered / 'a' / 'nanovi_report.pdf'
    data = pdf.read_bytes()
    assert data.startswith(b'%PDF')
    if shutil.which('pdffonts'):
        fonts = subprocess.run(['pdffonts', str(pdf)], capture_output=True, text=True).stdout
        lines = fonts.splitlines()[2:]
        assert 'Inter-Regular' in fonts and all(line.split()[-5] == 'yes' for line in lines)
        assert all('+Inter-' in line.split()[0] for line in lines), fonts  # no other font, e.g. Typst's raw font
    if shutil.which('pdfinfo'):
        pages = int(re.search(r'Pages:\s+(\d+)', subprocess.run(['pdfinfo', str(pdf)], capture_output=True,
                                                                  text=True).stdout).group(1))
        assert pages >= 7


def test_a_typst_failure_leaves_no_pdf_and_the_site(tmp_path):
    results = make_results(tmp_path / 'results')
    assets = tmp_path / 'assets'
    shutil.copytree(os.path.join(BIN, '..', 'assets', 'report'), assets)
    (assets / 'typst' / 'report.typ').write_text('#this-is-not-defined\n')
    res = subprocess.run([sys.executable, os.path.join(BIN, 'nanovi_report.py'), '--results', str(results),
                          '--outdir', str(tmp_path / 'out'), '--assets', str(assets), '--date', '2026-11-06'],
                         capture_output=True, text=True, env={**os.environ, 'MPLCONFIGDIR': str(tmp_path / 'mpl')})
    assert res.returncode == 0, res.stderr
    assert 'WARNING' in res.stderr and 'nanovi_report.pdf' in res.stderr
    assert not (tmp_path / 'out' / 'nanovi_report.pdf').exists()
    assert (tmp_path / 'out' / 'index.html').exists()
    assert 'nanovi_report.pdf' not in (tmp_path / 'out' / 'index.html').read_text(encoding='utf-8')
