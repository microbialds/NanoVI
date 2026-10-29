#!/usr/bin/env bash
# Fail when the pipeline's version strings disagree (F33). Run from the repository root.
set -euo pipefail

manifest=$(sed -n "s/^ *version *= *'\([^']*\)'.*/\1/p" nextflow.config)
changelog=$(sed -n 's/^## \[\([0-9][^]]*\)\].*/\1/p' CHANGELOG.md | head -n 1)
dockerfile=$(sed -n 's/^LABEL version="\([^"]*\)".*/\1/p' containers/Dockerfile)
singularity=$(sed -n 's/^ *Version \(.*\)$/\1/p' containers/Singularity.def)
python=$(sed -n 's/^__version__ = "\([^"]*\)".*/\1/p' bin/_version.py)
images=$(grep -o "ghcr.io/microbialds/nanovi-python:[A-Za-z0-9._-]*" conf/base.config | sed 's/.*://' | sort -u || true)

status=0
printf '%-48s %s\n' \
    "manifest.version (nextflow.config)" "${manifest:-MISSING}" \
    "top release (CHANGELOG.md)" "${changelog:-MISSING}" \
    "LABEL version (containers/Dockerfile)" "${dockerfile:-MISSING}" \
    "Version (containers/Singularity.def)" "${singularity:-MISSING}" \
    "__version__ (bin/_version.py)" "${python:-MISSING}" \
    "nanovi-python tag (conf/base.config)" "$(echo ${images:-MISSING})"
for v in "$changelog" "$dockerfile" "$singularity" "$python" "$images"; do
    if [ "$v" != "$manifest" ] || [ -z "$manifest" ]; then status=1; fi
done
if [ "$status" -ne 0 ]; then
    echo "Version strings disagree: set them all to the same version." >&2
else
    echo "All version strings are ${manifest}."
fi

# The cached minimap2 index is named after ext.minimap2_version (F26), which must
# match the minimap2 that the alignment label pins
mm2_ext=$(sed -n "s/^ *ext.minimap2_version *= *'\([^']*\)'.*/\1/p" conf/base.config)
mm2_conda=$(grep -o "bioconda::minimap2=[0-9.]*" conf/base.config | sed 's/.*=//' | sort -u || true)
printf '%-48s %s\n' \
    "ext.minimap2_version (conf/base.config)" "${mm2_ext:-MISSING}" \
    "minimap2 in the conda line (conf/base.config)" "$(echo ${mm2_conda:-MISSING})"
if [ -z "$mm2_ext" ] || [ "$mm2_ext" != "$mm2_conda" ]; then
    echo "ext.minimap2_version does not match the pinned minimap2." >&2
    status=1
fi
exit "$status"
