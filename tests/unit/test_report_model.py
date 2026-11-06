"""The run report's model (bin/report_model.py; F60, D42): every number equals its source table, and a missing
optional table leaves its section out."""
import os
import sys

import pandas as pd
import pytest

pytest.importorskip('yaml')

from .report_fixture import FASTPLONG, PLACEHOLDER_COUNTS, SPECIES, make_results

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
ASSETS = os.path.abspath(os.path.join(BIN, '..', 'assets', 'report'))
sys.path.insert(0, BIN)

import report_model as rm  # noqa: E402


def model_of(tmp_path, params=None, **kw):
    run = rm.load_results(make_results(tmp_path / 'results', **kw))
    data = rm.compute(run)
    return run, data, rm.build(run, data, [], ASSETS, params, {'REPORT': {'python': '3.13.0', 'typst': '0.15.1'}},
                               '2026-11-06')


def section(model, sid):
    return next((s for s in model['sections'] if s['id'] == sid), None)


def test_samples_in_natural_order(tmp_path):
    run, _, _ = model_of(tmp_path)
    assert run['samples'] == ['S1', 'S2', 'S10']


def test_metadata_orders_samples_by_group_then_naturally(tmp_path):
    results = make_results(tmp_path / 'results')
    meta = tmp_path / 'meta.csv'
    meta.write_text('sample_id,group,site\nS10,case,a\nS1,control,b\nS2,case,c\nX9,case,d\n')
    run = rm.load_results(results, meta)
    assert run['samples'] == ['S2', 'S10', 'S1']
    assert run['group'] == {'S1': 'control', 'S2': 'case', 'S10': 'case'}
    assert list(run['metadata']['sample_id']) == ['S10', 'S1', 'S2']


def test_read_fate_equals_the_species_tables_and_fastplong(tmp_path):
    _, data, _ = model_of(tmp_path)
    fate = data['fate'].set_index('sample')
    for s, rows in SPECIES.items():
        below, unassigned = PLACEHOLDER_COUNTS[s]
        raw, passed, short, long_ = FASTPLONG[s]
        assert fate.loc[s, 'species_reads'] == sum(r[2] for r in rows)
        assert fate.loc[s, 'below_detection_reads'] == below
        assert fate.loc[s, 'unassigned_reads'] == unassigned
        assert (fate.loc[s, 'raw_reads'], fate.loc[s, 'passed_reads']) == (raw, passed)
        assert (fate.loc[s, 'too_short'], fate.loc[s, 'too_long']) == (short, long_)
        total = sum(r[2] for r in rows) + below + unassigned
        assert fate.loc[s, 'pct_species'] == pytest.approx(100 * sum(r[2] for r in rows) / total)
        assert fate.loc[s, 'n_species'] == len(rows)


def test_read_fate_without_counts_keeps_fastplong_only(tmp_path):
    _, data, model = model_of(tmp_path, counts=False)
    fate = data['fate']
    assert fate['species_reads'].isna().all() and fate['pct_species'].isna().all()
    assert list(fate['passed_reads']) == [FASTPLONG[s][1] for s in ['S1', 'S2', 'S10']]
    assert not section(model, 'overview')['facts']['has_counts']


def test_composition_equals_the_consolidated_table(tmp_path):
    run, data, _ = model_of(tmp_path)
    otu = pd.read_csv(tmp_path / 'results' / 'consolidated_output' / 'otu_table_abundance.tsv', sep='\t',
                      index_col=0)
    comp = data['composition_species']
    for s in run['samples']:
        c = comp[comp['sample'] == s].set_index('taxon')['abundance']
        for t in data['shown_species']:
            assert c[t] == otu.loc[t, s]
        rest = otu.drop(index=list(data['shown_species']) + ['Below detection', 'Unassigned'], errors='ignore')[s]
        assert c['Other'] == pytest.approx(rest.sum())
    genus = data['composition_genus']
    s1 = genus[genus['sample'] == 'S1'].set_index('taxon')['abundance']
    assert s1['Alpha'] == pytest.approx(0.9) and s1['Beta'] == pytest.approx(0.1)


def test_top_taxa_is_the_union_of_each_samples_top_capped_by_mean():
    wide = pd.DataFrame({'a': [0.5, 0.3, 0.2, 0.0], 'b': [0.0, 0.1, 0.1, 0.8]}, index=['t1', 't2', 't3', 't4'])
    assert rm.top_taxa(wide, per_sample=1) == ['t4', 't1']
    assert rm.top_taxa(wide, per_sample=2) == ['t4', 't1', 't2']
    assert rm.top_taxa(wide, per_sample=2, cap=2) == ['t4', 't1']


def test_split_parts_gives_equal_parts_of_at_most_60():
    assert [len(p) for p in rm.split_parts(list(range(117)))] == [59, 58]
    assert [len(p) for p in rm.split_parts(list(range(60)))] == [60]
    assert [len(p) for p in rm.split_parts(list(range(61)))] == [31, 30]
    assert sum(rm.split_parts(list(range(130))), []) == list(range(130))


def test_flags_follow_the_shared_fraction_and_identical_16s(tmp_path):
    _, data, model = model_of(tmp_path)
    fl = data['flagged'].set_index('species')
    assert set(fl.index) == {'Alpha two', 'Beta three'}
    assert fl.loc['Alpha two', 'samples_shared'] == 1 and fl.loc['Alpha two', 'max_shared_read_fraction'] == 0.7
    assert fl.loc['Beta three', 'samples_identical_16s'] == 1
    assert fl.loc['Beta three', 'identical_16s_species'] == 'Beta four'
    facts = section(model, 'flagged')['facts']
    assert (facts['n_flagged'], facts['n_shared'], facts['n_identical']) == (2, 1, 1)


def test_optional_sections_appear_only_with_their_tables(tmp_path):
    _, _, model = model_of(tmp_path / 'plain')
    assert [s['id'] for s in model['sections']] == ['overview', 'composition', 'flagged', 'methods']
    _, _, model = model_of(tmp_path / 'all', groups=True, cell=True, reads=True)
    assert [s['id'] for s in model['sections']] == ['overview', 'composition', 'flagged', 'groups', 'cell', 'reads',
                                                    'methods']


def test_groups_table_lists_multi_member_groups_and_cell_columns_only_with_them(tmp_path):
    _, data, model = model_of(tmp_path / 'g', groups=True)
    gl = data['groups_long']
    assert list(gl['members']) == ['Alpha one; Alpha two'] and list(gl['sample']) == ['S1']
    assert 'cell_abundance' not in section(model, 'groups')['tables']['groups']['header']
    _, _, model = model_of(tmp_path / 'gc', groups=True, cell=True)
    t = section(model, 'groups')['tables']['groups']
    assert t['header'][-3:] == ['cell_abundance', 'copy_number', 'copy_number_level']
    assert t['rows'][0][t['header'].index('cell_abundance')] == '0.8'


def test_cell_section_shows_the_comment_lines_and_the_values(tmp_path):
    run, data, model = model_of(tmp_path, cell=True)
    sec = section(model, 'cell')
    assert sec['facts']['explanation'].startswith('Cell-level abundance: read-level abundance / 16S copy number')
    left = sec['tables']['left_out']
    assert left['rows'][0] == ['S1', 'Not included: Below detection (5 reads) and Unassigned (5 reads).']
    cl = data['cell_long']
    s1 = cl[cl['sample'] == 'S1'].set_index('species')
    assert s1['abundance'].sum() == pytest.approx(1.0)
    assert s1.loc['Beta three', 'read_abundance'] == 0.1 and s1.loc['Beta three', 'copy_number'] == 2.0


def test_read_assignment_summary_counts_reads_and_pairs(tmp_path):
    _, data, _ = model_of(tmp_path, reads=True)
    per = data['reads'].set_index('sample')
    assert per.loc['S1', 'reads'] == 5
    assert per.loc['S1', 'p >= 0.99'] == 1 and per.loc['S1', '0.5 <= p < 0.99'] == 1 and per.loc['S1', 'p < 0.5'] == 1
    assert per.loc['S1', 'top: Below detection'] == 1 and per.loc['S1', 'top: Unassigned'] == 1
    assert per.loc['S1', 'reads on 2+ species'] == 2
    assert per.loc['S10', 'p >= 0.99'] == 2 and per.loc['S2', 'top: Unassigned'] == 1
    pairs = data['read_pairs']
    assert set(map(tuple, pairs[['species_a', 'species_b', 'reads']].values.tolist())) == {
        ('Alpha one', 'Alpha two', 1), ('Alpha one', 'Beta three', 1)}


def test_overview_totals_equal_the_read_fate(tmp_path):
    _, data, model = model_of(tmp_path)
    f = section(model, 'overview')['facts']
    assert f['raw_reads'] == str(sum(v[0] for v in FASTPLONG.values()))
    assert f['passed_reads'] == str(sum(v[1] for v in FASTPLONG.values()))
    assert f['species_reads'] == str(200) and f['n_species'] == 4
    assert section(model, 'overview')['lists']['attention'][0].startswith('1 sample with no reported species: S2')


def test_methods_list_tools_with_checked_references_and_credit_emu(tmp_path):
    _, _, model = model_of(tmp_path, params={'estimator': 'vb', 'N': 3})
    me = section(model, 'methods')
    tools = me['tables']['tools']
    names = [r[0] for r in tools['rows']]
    assert names[:2] == ['NanoVI', 'Nextflow']
    assert {'fastplong', 'minimap2', 'SAMtools', 'NumPy', 'SciPy', 'Typst'} <= set(names)
    assert 'Python' not in names
    emu = names.index('Emu (read likelihood and EM estimator)')
    assert tools['rows'][emu][2] == 'method' and tools['keys'][emu] == ['curry2022']
    cited = [r['key'] for r in me['lists']['references']]
    assert 'curry2022' in cited and 'curiqueo2026' in cited and 'li2018' in cited
    assert '@article{curry2022,' in me['facts']['bibtex']
    params = me['tables']['params']['rows']
    assert ['Estimator', '--estimator', 'vb'] in params and ['Secondary alignments (-N)', '--N', '3'] in params


def test_every_cited_key_is_in_references_and_has_a_doi():
    refs, tool_refs = rm.load_references(ASSETS)
    for tool, info in tool_refs.items():
        for k in info['keys']:
            assert k in refs, (tool, k)
    for k, r in refs.items():
        assert r.get('doi') and r.get('title') and r.get('year'), k


def test_a_missing_software_versions_file_still_gives_the_report_s_own(tmp_path):
    _, _, model = model_of(tmp_path, versions=False)
    names = [r[0] for r in section(model, 'methods')['tables']['tools']['rows']]
    assert 'Typst' in names and 'NanoVI' not in names


def test_pdf_tables_are_capped_with_a_note():
    nums = rm.Numbering()
    df = pd.DataFrame({'a': range(100), 'b': range(100)})
    t = rm.table_entry(nums, df, 'x', 'caption', tsv='x.tsv')
    assert len(t['rows']) == 100 and len(t['pdf']['rows']) == rm.PDF_MAX_ROWS
    assert 'first 60 of 100 rows' in t['pdf']['note'] and 'x.tsv' in t['pdf']['note']


def test_a_run_with_no_flagged_species_gives_an_empty_flag_table(tmp_path, monkeypatch):
    from . import report_fixture as fx
    monkeypatch.setitem(fx.SPECIES, 'S1', [('Alpha one', 0.6, 60, '', 0.1), ('Alpha two', 0.4, 40, '', 0.2)])
    monkeypatch.setitem(fx.SPECIES, 'S10', [('Gamma five', 1.0, 100, '', 0.0)])
    _, data, model = model_of(tmp_path)
    assert data['flagged'].empty and list(data['flagged'].columns) == rm.FLAGGED_COLUMNS
    sec = section(model, 'flagged')
    assert sec['facts']['n_flagged'] == 0 and sec['tables']['flagged']['rows'] == []


def test_overview_share_uses_the_tables_total_without_fastplong(tmp_path):
    results = make_results(tmp_path / 'results')
    for f in (results / 'qc' / 'fastplong').iterdir():
        f.unlink()
    run = rm.load_results(results)
    data = rm.compute(run)
    f = section(rm.build(run, data, [], ASSETS), 'overview')['facts']
    assert f['passed_reads'] == '' and f['counted_reads'] == '240'
    assert f['pct_species'] == f'{100 * 200 / 240:.1f}'


def test_counts_read_the_same_as_int_or_float():
    assert rm.fmt(5000) == rm.fmt(5000.0) == '5,000'
    assert rm.fmt(999) == rm.fmt(999.0) == '999'


def test_read_assignment_summary_streams_the_file(tmp_path):
    _, data, _ = model_of(tmp_path, reads=True)
    run = rm.load_results(tmp_path / 'results')
    lines = list(rm.read_groups(run['reads']['S1']))
    assert [r for r, _ in lines] == ['r1', 'r2', 'r3', 'r4', 'r5']
    assert lines[2][1] == [('Alpha two', 0.3), ('Alpha one', 0.3), ('Beta three', 0.4)]
