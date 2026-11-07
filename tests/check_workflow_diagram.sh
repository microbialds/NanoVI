#!/bin/bash
#
# Keeps docs/workflow.md in step with the code (F59): every process or subworkflow
# that main.nf or a local subworkflow includes (under its alias when it has one) must
# be named in one of the page's Mermaid diagrams or tables (prose does not count), and
# every process the page shows (an upper-case Mermaid node id, or a backticked
# upper-case name in any column of a table) must still be included somewhere.
# Includes inside // or /* */ comments do not count.
# Functions (validateParameters, samplesheetToList) are camelCase and are not checked;
# file nodes have mixed-case ids. Adapted from gene2dis/BugBuster's
# tests/bin/test_workflow_diagram.sh (dev @ de57f34). Run by CI's version-check job.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/.." && pwd )"
DOC="${DIAGRAM_DOC:-${REPO_DIR}/docs/workflow.md}"

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

PASS=0
FAIL=0

# Names included by the given Nextflow files (default: main.nf and the local
# subworkflows): the alias after "as" when there is one, else the name; include blocks
# may span lines and hold several names separated by ";". Comments are removed first
included_names() {
    if [ "$#" -eq 0 ]; then
        set -- "${REPO_DIR}/main.nf" "${REPO_DIR}"/subworkflows/local/*.nf
    fi
    cat "$@" \
        | sed -E 's#//.*$##' \
        | tr '\n' ' ' \
        | sed -E 's#/\*([^*]|\*+[^*/])*\*+/##g' \
        | grep -oE 'include[[:space:]]*\{[^}]*\}' \
        | sed -E 's/^include[[:space:]]*\{//; s/\}$//' \
        | tr ';' '\n' \
        | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//' \
        | awk 'NF { print ($2 == "as" ? $3 : $1) }' \
        | grep -E '^[A-Z][A-Z0-9_]*$' \
        | sort -u
}

# check_doc <doc>: prints one line per missing or stale name, returns 1 if any
check_doc() {
    local doc="$1" bad=0 name
    # the Mermaid blocks and the table rows
    awk '/^```mermaid/ { m = 1; next } m && /^```/ { m = 0 } m || /^\|/' "${doc}" > "${TMP_DIR}/named.txt"
    while read -r name; do
        if ! grep -qw -- "${name}" "${TMP_DIR}/named.txt"; then
            echo "missing: ${name} is included by the code but not in a diagram or table of ${doc##*/}"
            bad=1
        fi
    done < "${TMP_DIR}/included.txt"
    # the processes the page shows: Mermaid node ids (labels, edge texts and style
    # lines removed) and the backticked upper-case names in any column of a table
    {
        awk '/^```mermaid/ { m = 1; next } m && /^```/ { m = 0 } m' "${doc}" \
            | grep -vE '^[[:space:]]*(%%|classDef|class|flowchart|direction|subgraph|end)\b' \
            | sed -E 's/"[^"]*"//g; s/\|[^|]*\|//g; s/\[[^]]*\]//g; s/\([^)]*\)//g; s/\{[^}]*\}//g' \
            | grep -oE '\b[A-Z][A-Z0-9_]+\b' || true
        awk -F'|' '/^\|/ { for (i = 2; i < NF; i++) print $i }' "${doc}" | grep -oE '`[A-Z][A-Z0-9_]+`' | tr -d '`' || true
    } | sort -u > "${TMP_DIR}/shown.txt"
    while read -r name; do
        if ! grep -qx -- "${name}" "${TMP_DIR}/included.txt"; then
            echo "stale: ${doc##*/} shows ${name}, which the code no longer includes"
            bad=1
        fi
    done < "${TMP_DIR}/shown.txt"
    return "${bad}"
}

included_names > "${TMP_DIR}/included.txt"
n_names=$(wc -l < "${TMP_DIR}/included.txt")

# A commented-out include is not an include
printf '%s\n' "include { FILTER_READS } from './a'" "// include { OLD_LINE } from './b'" \
    "/* include { OLD_BLOCK } from './c'" "   include { OLD_BLOCK2 } from './d' */" > "${TMP_DIR}/comments.nf"
if [ "$(included_names "${TMP_DIR}/comments.nf" | tr '\n' ' ')" = "FILTER_READS " ]; then
    echo "✓ includes inside // and /* */ comments are not counted"
    PASS=$((PASS + 1))
else
    echo "✗ commented-out includes were counted: $(included_names "${TMP_DIR}/comments.nf" | tr '\n' ' ')"
    FAIL=$((FAIL + 1))
fi

# The parser must find the subworkflows and processes the pipeline is built from
if [ "${n_names}" -ge 16 ] && grep -qx ABUNDANCE "${TMP_DIR}/included.txt" \
        && grep -qx REPORT "${TMP_DIR}/included.txt" \
        && grep -qx GENERATE_PHYLOSEQ_TABLES "${TMP_DIR}/included.txt" \
        && grep -qx BUILD_DB "${TMP_DIR}/included.txt" \
        && ! grep -qx validateParameters "${TMP_DIR}/included.txt"; then
    echo "✓ ${n_names} process and subworkflow names read from main.nf and subworkflows/local"
    PASS=$((PASS + 1))
else
    echo "✗ include parsing: ${n_names} names, expected >= 16 with ABUNDANCE, REPORT, GENERATE_PHYLOSEQ_TABLES and BUILD_DB"
    FAIL=$((FAIL + 1))
fi

# The real page
if check_doc "${DOC}" > "${TMP_DIR}/doc.log"; then
    echo "✓ ${DOC#"${REPO_DIR}"/} names every included process and subworkflow, and no stale one"
    PASS=$((PASS + 1))
else
    echo "✗ ${DOC#"${REPO_DIR}"/} is out of step with the code:"
    sed 's/^/    /' "${TMP_DIR}/doc.log"
    FAIL=$((FAIL + 1))
fi

# The check itself fails on a missing and on a stale name
sed -E 's/\bIDENTICAL_16S\b/IDENTICAL_GONE/g' "${DOC}" > "${TMP_DIR}/missing.md"
if ! check_doc "${TMP_DIR}/missing.md" > "${TMP_DIR}/missing.log" && grep -q 'missing: IDENTICAL_16S' "${TMP_DIR}/missing.log"; then
    echo "✓ a process left out of the page is reported"
    PASS=$((PASS + 1))
else
    echo "✗ removing IDENTICAL_16S from the page was not reported"
    FAIL=$((FAIL + 1))
fi
# A box left in a diagram and a row left in a table after the code drops the process
# are reported; a file name or an acronym in backticks in prose is not a process
awk '{ print } /^```mermaid/ && !done { print "    LIKELIHOOD --> GET_CIGAR[\"GET_CIGAR\"]"; done = 1 }' "${DOC}" > "${TMP_DIR}/stale_box.md"
{ cat "${DOC}"; printf '\n| `COMPUTE_LOGP` | `vi_stats` | removed process | |\n'; } > "${TMP_DIR}/stale_row.md"
sed -E '/^\| `ABUNDANCE` \|/ s/`FILTER_READS`/`FILTER_READS`, `GET_CIGAR`/' "${DOC}" > "${TMP_DIR}/stale_cell.md"
for case in "stale_box.md GET_CIGAR" "stale_row.md COMPUTE_LOGP" "stale_cell.md GET_CIGAR"; do
    set -- ${case}
    if ! check_doc "${TMP_DIR}/$1" > "${TMP_DIR}/stale.log" && grep -q "stale: .*shows $2," "${TMP_DIR}/stale.log"; then
        echo "✓ $2, which the code no longer includes, is reported ($1)"
        PASS=$((PASS + 1))
    else
        echo "✗ a stale $2 was not reported ($1)"
        FAIL=$((FAIL + 1))
    fi
done
{ cat "${DOC}"; echo 'The table is a `TSV` with `CIGAR` statistics from the `BAM`.'; } > "${TMP_DIR}/acronym.md"
if check_doc "${TMP_DIR}/acronym.md" > "${TMP_DIR}/acronym.log"; then
    echo "✓ acronyms in backticks in prose are not taken for processes"
    PASS=$((PASS + 1))
else
    echo "✗ an acronym in prose was reported:"
    sed 's/^/    /' "${TMP_DIR}/acronym.log"
    FAIL=$((FAIL + 1))
fi
# A name left only in prose (LIKELIHOOD taken out of every diagram and table, then
# mentioned in a sentence) counts as missing
sed -E '/^\|/ s/\bLIKELIHOOD\b/LIKE/g; /^```mermaid/,/^```/ s/\bLIKELIHOOD\b/LIKE/g' "${DOC}" > "${TMP_DIR}/prose.md"
echo 'LIKELIHOOD is mentioned in prose only.' >> "${TMP_DIR}/prose.md"
if ! check_doc "${TMP_DIR}/prose.md" > "${TMP_DIR}/prose.log" && grep -q 'missing: LIKELIHOOD' "${TMP_DIR}/prose.log"; then
    echo "✓ a process named only in prose is reported"
    PASS=$((PASS + 1))
else
    echo "✗ LIKELIHOOD named only in prose was not reported"
    FAIL=$((FAIL + 1))
fi

echo
echo "${PASS} passed, ${FAIL} failed"
[ "${FAIL}" -eq 0 ]
