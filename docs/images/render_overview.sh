#!/bin/bash
#
# Refreshes docs/images/nanovi_overview.svg's embedded fonts and renders
# nanovi_overview.png from it, in two pinned containers (nothing installed on
# the host; needs docker). Run it after editing the SVG. Adapted from
# gene2dis/BugBuster's docs/images/render_overview.sh (dev @ de57f34).
#
#  1. The run report's image, ghcr.io/microbialds/nanovi-report under the tag
#     conf/base.config gives it (Inter 4.1 found by bin/figure_style.py's
#     font_dir(), fontTools 4.68): subsets Inter Regular and SemiBold to the
#     characters the SVG's text uses and writes them as @font-face data into its
#     <style id="fonts"> block, so the SVG shows Inter where Inter is not
#     installed (GitHub) and its text stays text. WOFF 1.0, not figure_style's
#     WOFF2: subsetting to WOFF2 needs brotli, which the image does not have.
#     SOURCE_DATE_EPOCH fixes the fonts' timestamps, so two runs write the same file.
#     Until a release has published that tag, build it locally first:
#     docker build -t ghcr.io/microbialds/nanovi-report:<tag> containers/report/
#  2. librsvg + font-inter 4.1 (Seqera Containers, conda-forge librsvg=2.62.3
#     font-inter=4.1 fontconfig; sha256:d5ea0d67...): rsvg-convert at 1.75x
#     (2,240 px wide). librsvg takes Inter from fontconfig, not from the
#     @font-face data.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
BIN_DIR="${REPO_DIR}/bin"
REPORT_IMAGE=$(grep -oE "ghcr.io/microbialds/nanovi-report:[0-9.]+" "${REPO_DIR}/conf/base.config" | head -n 1)
RSVG_IMAGE="community.wave.seqera.io/library/librsvg_font-inter_fontconfig@sha256:d5ea0d6797046ae4dab9a7a5fba9d1293fce013e98029a5672687aa6ef36b03b"
SVG="nanovi_overview.svg"
PNG="nanovi_overview.png"

if [ -z "${REPORT_IMAGE}" ]; then
    echo "no nanovi-report image tag in conf/base.config" >&2
    exit 1
fi

docker run --rm -i -u "$(id -u):$(id -g)" -v "${SCRIPT_DIR}:/work" -v "${BIN_DIR}:/nv_bin:ro" -w /work \
    -e MPLCONFIGDIR=/tmp -e SOURCE_DATE_EPOCH=1767225600 "${REPORT_IMAGE}" python -I - "${SVG}" <<'EOF'
import base64, html, io, re, sys
from pathlib import Path
from fontTools import subset
from fontTools.ttLib import TTFont
sys.path.insert(0, "/nv_bin")
from figure_style import font_dir

svg_path = Path(sys.argv[1])
svg = svg_path.read_text(encoding="utf-8")
body = re.sub(r'<style id="fonts">.*?</style>', '', svg, flags=re.S)
chars = set()
for text in re.findall(r"<text\b[^>]*>(.*?)</text>", body, flags=re.S):
    chars.update(html.unescape(re.sub(r"<[^>]+>", "", text)))
chars.add(" ")
fonts = font_dir()
rules = []
for face, weight in (("Regular", 400), ("SemiBold", 600)):
    font = TTFont(fonts / f"Inter-{face}.ttf")
    options = subset.Options()
    options.flavor = "woff"
    options.layout_features = ["kern", "liga", "calt"]
    options.name_IDs = ["*"]
    sub = subset.Subsetter(options)
    sub.populate(text="".join(sorted(chars)))
    sub.subset(font)
    out = io.BytesIO()
    font.flavor = "woff"
    font.save(out)
    data = base64.b64encode(out.getvalue()).decode()
    rules.append(f"@font-face{{font-family:Inter;font-weight:{weight};font-style:normal;"
                 f"src:url(data:font/woff;base64,{data}) format('woff');}}")
svg, n = re.subn(r'<style id="fonts">.*?</style>', lambda m: '<style id="fonts">' + "".join(rules) + "</style>",
                 svg, count=1, flags=re.S)
if n != 1:
    sys.exit(f'{svg_path}: no <style id="fonts"></style> element to fill; add it after the <desc> element')
svg_path.write_text(svg, encoding="utf-8")
print(f"{svg_path}: Inter 400 + 600 subset to {len(chars)} characters, {svg_path.stat().st_size // 1024} KB")
EOF

docker run --rm -u "$(id -u):$(id -g)" -v "${SCRIPT_DIR}:/work" -w /work "${RSVG_IMAGE}" \
    rsvg-convert --zoom 1.75 --format png --output "${PNG}" "${SVG}"
echo "${PNG}: $(( $(wc -c < "${SCRIPT_DIR}/${PNG}") / 1024 )) KB"
