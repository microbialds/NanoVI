"""Interactive Vega-Lite (Altair) versions of the run report's read-fate and composition figures (D42), for the
HTML site: hover for values, a group menu when --report_metadata gives groups, zoom on the read-count axis. Built
from the same data as the static figure (its record's 'interactive' entry), with every sample in one chart, at a
screen scale (type 1.8x the print size). Adapted from BugBuster's bin/summary_interactive.py (gene2dis/BugBuster dev
@ de57f34). A spec whose JSON is larger than SPEC_BUDGET is not embedded: the page shows the static figure.
"""

import json

import numpy as np

import figure_style as fs
from report_model import warn

SPEC_BUDGET = 1_000_000  # bytes of spec JSON (data inline) per figure
SCALE = 1.8
STEP = 9 * SCALE         # px per sample row: the 6 pt label height at screen scale


def available():
    try:
        import altair  # noqa: F401 (availability check)
        import vl_convert  # noqa: F401
    except ImportError:
        return False
    return True


def config():
    small, base = fs.SMALL_PT * SCALE, fs.BASE_PT * SCALE
    return {
        "font": fs.FONT, "background": "white", "view": {"stroke": None},
        "axis": {"labelFontSize": small, "titleFontSize": base, "labelFont": fs.FONT, "titleFont": fs.FONT,
                 "domainWidth": 0.8, "grid": False, "tickWidth": 0.8, "tickSize": 3, "titleFontWeight": "normal"},
        "legend": {"labelFontSize": small, "titleFontSize": small, "titleFontWeight": "normal", "labelFont": fs.FONT,
                   "titleFont": fs.FONT, "gradientThickness": 10},
        "header": {"labelFontSize": small, "labelFont": fs.FONT},
        "range": {"category": fs.CATEGORY},
    }


def _group_select(alt, groups):
    names = list(dict.fromkeys(g for g in groups if g != ""))
    if len(names) < 2:
        return None
    return alt.selection_point(fields=["group"], name="grp", bind=alt.binding_select(
        options=[None] + names, labels=["all"] + names, name="Group "))


def _tip_group(frame):
    import altair as alt
    return [alt.Tooltip("group:N")] if (frame["group"] != "").any() else []


def read_fate(d):
    """Shares of the passed reads by fate, and raw and passed reads on a log axis (drag or scroll to zoom)."""
    import altair as alt
    f = d["frame"].copy()
    samples = d["samples"]
    y = alt.Y("sample:N", sort=samples, title=None, axis=alt.Axis(ticks=False, domain=False, labelOverlap=False))
    sel = _group_select(alt, f["group"])
    panels = []
    if d["has_pct"]:
        names = {"pct_species": "reported species", "pct_below_detection": "Below detection",
                 "pct_unassigned": "Unassigned"}
        long = f.melt(id_vars=["sample", "group"], value_vars=list(names), var_name="fate", value_name="pct")
        long["fate"] = long["fate"].map(names)
        long["order"] = long["fate"].map({v: i for i, v in enumerate(names.values())})
        left = alt.Chart(long).mark_bar().encode(
            x=alt.X("pct:Q", stack="zero", title="% of the reads that passed fastplong",
                    scale=alt.Scale(domain=[0, 100])),
            y=y, order="order:Q",
            color=alt.Color("fate:N", scale=alt.Scale(domain=list(names.values()),
                                                      range=[fs.CATEGORY[0], fs.CATEGORY[1], "#9A9A9A"]),
                            legend=alt.Legend(title=None, orient="bottom", columns=3, labelLimit=400)),
            tooltip=[alt.Tooltip("sample:N")] + _tip_group(f) + [alt.Tooltip("fate:N"),
                                                                 alt.Tooltip("pct:Q", format=".2f")])
        if sel is not None:
            left = left.add_params(sel).transform_filter(sel)
        panels.append(left.properties(width=300, height=alt.Step(STEP)))
    if d["has_raw"]:
        lo, hi = d["domain"]
        zoom = alt.selection_interval(bind="scales", encodings=["x"])
        base = alt.Chart(f)
        raw = base.mark_bar(color="#D0D0D0").encode(
            x=alt.X("raw_reads:Q", scale=alt.Scale(type="log", domain=[lo, hi]), title="reads (log scale)",
                    axis=alt.Axis(format="~s")),
            x2=alt.X2(datum=lo), y=y if not panels else alt.Y("sample:N", sort=samples, axis=None),
            tooltip=[alt.Tooltip("sample:N")] + _tip_group(f) + [alt.Tooltip("raw_reads:Q", format=","),
                                                                 alt.Tooltip("passed_reads:Q", format=",")])
        passed = base.mark_bar(color=fs.CATEGORY[2], size=STEP * 0.45).encode(
            x="passed_reads:Q", x2=alt.X2(datum=lo),
            y=y if not panels else alt.Y("sample:N", sort=samples, axis=None),
            tooltip=[alt.Tooltip("sample:N")] + _tip_group(f) + [alt.Tooltip("raw_reads:Q", format=","),
                                                                 alt.Tooltip("passed_reads:Q", format=",")])
        layers = [raw.add_params(zoom), passed]
        if sel is not None and not panels:
            layers[0] = layers[0].add_params(sel)
        if sel is not None:
            layers = [x.transform_filter(sel) for x in layers]
        panels.append(alt.layer(*layers).properties(width=260, height=alt.Step(STEP)))
    return alt.hconcat(*panels, spacing=10) if len(panels) > 1 else panels[0]


def bars(d):
    """Stacked composition bars, one per sample, taxa in the static figure's order and colours."""
    import altair as alt
    f = d["frame"].copy()
    f["pct"] = 100 * f["abundance"]
    order = {t: i for i, t in enumerate(d["taxa"])}
    f["order"] = f["taxon"].map(order)
    sel = _group_select(alt, f["group"])
    chart = alt.Chart(f).mark_bar().encode(
        x=alt.X("pct:Q", stack="zero", title="% of the reads on reported species", scale=alt.Scale(domain=[0, 100])),
        y=alt.Y("sample:N", sort=d["samples"], title=None, axis=alt.Axis(ticks=False, domain=False,
                                                                         labelOverlap=False)),
        order="order:Q",
        color=alt.Color("taxon:N", sort=d["taxa"], scale=alt.Scale(domain=d["taxa"], range=d["colors"]),
                        legend=alt.Legend(title=d["rank"], orient="right", labelLimit=320, symbolType="square")),
        tooltip=[alt.Tooltip("sample:N")] + _tip_group(f) + [alt.Tooltip("taxon:N", title=d["rank"]),
                                                             alt.Tooltip("pct:Q", format=".2f", title="%")])
    if sel is not None:
        chart = chart.add_params(sel).transform_filter(sel)
    return chart.properties(width=520, height=alt.Step(STEP))


def heatmap(d):
    """Species x samples heatmap of log10 abundance; grey cells not detected; samples as columns."""
    import altair as alt
    df = d["frame"].copy()
    feat = d["feature"]
    df["log10"] = np.where(df["abundance"] > 0, np.log10(df["abundance"].where(df["abundance"] > 0)), np.nan)
    s_enc = {"field": "sample", "type": "nominal", "sort": d["samples"], "title": None}
    f_enc = {"field": feat, "type": "nominal", "sort": d["rows"], "title": None}
    enc = {"x": alt.X(**s_enc, axis=alt.Axis(ticks=False, domain=False, labelAngle=-90, labelOverlap=False)),
           "y": alt.Y(**f_enc, axis=alt.Axis(ticks=False, domain=False, labelLimit=320, labelFontStyle="italic"))}
    tip = [alt.Tooltip("sample:N")] + _tip_group(df) + [alt.Tooltip(f"{feat}:N"),
                                                        alt.Tooltip("abundance:Q", format=".3g")]
    base = alt.Chart(df)
    sel = _group_select(alt, df["group"])
    bg = base.mark_rect(color=fs.NOT_DETECTED).encode(**enc, tooltip=tip)
    cells = base.transform_filter("isValid(datum.log10)").mark_rect().encode(
        **enc, tooltip=tip,
        color=alt.Color("log10:Q", scale=alt.Scale(scheme="viridis", domain=d["domain"], clamp=True),
                        legend=alt.Legend(title="log10 relative abundance", orient="bottom", direction="horizontal",
                                          gradientLength=160)))
    layers = [bg, cells]
    if sel is not None:
        layers = [layers[0].add_params(sel).transform_filter(sel), layers[1].transform_filter(sel)]
    sstep = max(STEP, min(36, 720 / max(len(d["samples"]), 1)))
    return alt.layer(*layers).properties(width=alt.Step(sstep), height=alt.Step(STEP))


BUILDERS = {"read_fate": read_fate, "bars": bars, "heatmap": heatmap}


def spec_json(key, data):
    """The figure's Vega-Lite spec as JSON text, or None (no interactive version, or above the size budget)."""
    if not data or data.get("kind") not in BUILDERS:
        return None
    import altair as alt
    alt.data_transformers.disable_max_rows()
    spec = BUILDERS[data["kind"]](data).to_dict()
    spec["config"] = config()
    spec["padding"] = {"left": 12, "top": 4, "right": 4, "bottom": 4}
    # inline in a <script>: no '<', '>' or '&' may appear raw (a label with '</script>' would end the block);
    # \u escapes stay valid JSON and JavaScript
    text = (json.dumps(spec, separators=(",", ":"), sort_keys=True, default=str)
            .replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026"))
    if len(text) > SPEC_BUDGET:
        warn(f"{key}: interactive spec {len(text) / 1e6:.1f} MB is above the {SPEC_BUDGET / 1e6:g} MB budget; the "
             "HTML shows the static figure")
        return None
    return text


def javascript_bundle():
    """vega + vega-lite + vega-embed in one script (vl-convert's bundle), for offline use."""
    import vl_convert as vlc
    return vlc.javascript_bundle()
