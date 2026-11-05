"""
Ambiguity groups: table rows that the sample cannot separate (F57, D39).

estimate_abundance.py --groups links two rows of the abundance table a and b
when the reads they share carry a large part of either one's estimate:

    f_a(b) = Σ_r φ_ra · [read r's share on b ≥ SHARED_READ_MIN_SPLIT] / E_a

where φ_ra is read r's responsibility on row a (summed over the row's tax IDs,
from the final fit: VB's φ or EM's posterior) and E_a = Σ_r φ_ra its expected
reads. a and b are linked when max(f_a(b), f_b(a)) ≥ GROUP_MIN_SHARED; the
maximum, so that a small row whose estimate rests on reads it shares with a
large one is caught. Groups are the connected components of the links, so every
row is in exactly one group (a row that links to nothing is a group of one).
As for shared_read_fraction, a read's share on b is judged among the reported
rows. write_abundance_table.py --groups-output writes them as
<sample>_groups.tsv, annotated with the species that share an identical 16S
with a member in the database and are not reported (identical_16s.tsv, D24).

f is per partner, so it is not shared_read_fraction, which counts a read as
shared when all the other rows together take SHARED_READ_MIN_SPLIT of it.
"""

import numpy as np
import pandas as pd
from scipy import sparse

from variational_inference_bayesian import SHARED_READ_MIN_SPLIT, row_groups

# Threshold on max(f_a(b), f_b(a)) (D39); estimate_abundance.py --group-min-shared is for tests
GROUP_MIN_SHARED = 0.25


def table_rows(reported_ids, row_keys):
    """The table row (0, 1, …, in order of first appearance) of each reported column, and the row count.

    The numbering of variational_inference_bayesian.row_groups, so the groups and
    the read summary agree on rows; without row_keys every species is a row of its own.
    """
    rows = row_groups(list(reported_ids), row_keys) or list(range(len(reported_ids)))
    return rows, max(rows) + 1 if rows else 0


def row_shares(idx, phi, row_of_column, n_rows):
    """Each read's responsibility on each row, as a sparse (reads x rows) matrix, and each read's total.

    idx, phi: the padded responsibilities, idx the entry's column among the
    reported species or -1 (padding, or a species not reported); entries with
    idx -1 are left out, so a read's total is its responsibility on reported rows.
    The tax IDs of a row add up.
    """
    reads, cols = np.nonzero(idx >= 0)
    rows = np.asarray(row_of_column, dtype=np.int64)[idx[reads, cols]]
    share = sparse.csr_matrix((phi[reads, cols], (reads, rows)), shape=(idx.shape[0], n_rows))
    share.sum_duplicates()
    return share, np.asarray(share.sum(axis=1)).ravel()


def _entry_fraction(share, total):
    """Each stored entry of share over its read's total (nan where the total is 0)."""
    read_of_entry = np.repeat(np.arange(share.shape[0]), np.diff(share.indptr))
    with np.errstate(divide="ignore", invalid="ignore"):
        return share.data / total[read_of_entry]


def pairwise_shared(idx, phi, row_of_column, n_rows, shares=None):
    """Expected reads per row and the reads each row shares with each other row.

    Returns (expected, shared): expected[a] = Σ_r φ_ra and
    shared[a, b] = Σ_r φ_ra over the reads whose share on b, among the reported
    rows, is at least SHARED_READ_MIN_SPLIT (shared[a, a] = 0). shares: the
    result of row_shares, when the caller has it.
    """
    share, total = shares if shares is not None else row_shares(idx, phi, row_of_column, n_rows)
    is_split = _entry_fraction(share, total) >= SHARED_READ_MIN_SPLIT
    split = sparse.csr_matrix((is_split.astype(np.float64), share.indices, share.indptr), shape=share.shape)
    expected = np.asarray(share.sum(axis=0)).ravel()
    shared = np.asarray((share.T @ split).todense())
    np.fill_diagonal(shared, 0.0)
    return expected, shared


def link_groups(expected, shared, min_shared=GROUP_MIN_SHARED):
    """Group number of each row and, per group, the largest f between two members.

    Rows a and b are linked when max(f_a(b), f_b(a)) ≥ min_shared; groups are
    numbered by their first row. The largest f of a group of one is None.
    """
    n = len(expected)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(expected[:, None] > 0, shared / expected[:, None], 0.0)
    strength = np.maximum(f, f.T)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(*np.nonzero(np.triu(strength >= min_shared, k=1))):
        parent[find(a)] = find(b)
    roots = [find(a) for a in range(n)]
    number = {}
    group_of_row = [number.setdefault(r, len(number)) for r in roots]
    members = [[] for _ in number]
    for a, g in enumerate(group_of_row):
        members[g].append(a)
    max_shared = [float(np.triu(strength[np.ix_(m, m)], k=1).max()) if len(m) > 1 else None for m in members]
    return group_of_row, max_shared


def ambiguity_groups(idx, phi, reported_ids, row_keys=None, min_shared=GROUP_MIN_SHARED):
    """The groups of the reported species, as written into the abundance JSON.

    A list with one entry per group: "tax_ids" (its species), "shared_reads"
    (their expected reads on reads split with rows outside the group: the
    group's share of the read, where the other groups take at least
    SHARED_READ_MIN_SPLIT of it, the D27 rule judged between groups) and
    "max_shared" (the largest f between two of its rows; None for a group of one
    row). Expected reads and α are in the JSON's read summary.
    """
    row_of_column, n_rows = table_rows(reported_ids, row_keys)
    share, total = row_shares(idx, phi, row_of_column, n_rows)
    expected, shared = pairwise_shared(idx, phi, row_of_column, n_rows, shares=(share, total))
    group_of_row, max_shared = link_groups(expected, shared, min_shared)
    n_groups = len(max_shared)
    to_group = sparse.csr_matrix((np.ones(n_rows), (np.arange(n_rows), group_of_row)), shape=(n_rows, n_groups))
    group_share = share @ to_group
    group_share.sum_duplicates()
    split = (group_share.data > 0) & (_entry_fraction(group_share, total) <= 1.0 - SHARED_READ_MIN_SPLIT)
    group_shared = np.bincount(group_share.indices[split], weights=group_share.data[split], minlength=n_groups)
    groups = [{"tax_ids": [], "shared_reads": float(group_shared[g]), "max_shared": m}
              for g, m in enumerate(max_shared)]
    for sid, r in zip(reported_ids, row_of_column):
        groups[group_of_row[r]]["tax_ids"].append(sid)
    return groups


GROUP_TABLE_COLUMNS = ["group", "members", "n_members", "abundance", "lowest_common_rank", "lowest_common_taxon"]
RELATIVES_COLUMN = "identical_16s_relatives"


def write_groups_table(groups, read_summary, species_df, taxonomy_df, path, binomial=False, intervals=False):
    """Write <sample>_groups.tsv from the JSON's groups and the sample's species table.

    species_df: the table freq_to_lineage_df wrote (its rows, lineages,
    estimated_counts and identical_16s_species columns, when present). Each
    group's rows are found there through the tax IDs' row keys
    (taxonomy.species_row_keys), so a group's abundance and estimated_counts are
    the sums of its rows' values. Groups go by decreasing abundance, then the
    Below detection and Unassigned rows of the species table. The group's
    shared_read_fraction is its shared reads over its expected reads, and with
    intervals ci_lower and ci_upper come from its summed α, as for a summed row;
    a group of one row takes them from its species row, so the two tables agree
    exactly. A tax ID missing from the read summary counts 0, as in the species table.
    """
    from taxonomy import (BELOW_DETECTION_LABEL, IDENTICAL_COLUMN, INTERVAL_COLUMNS, SHARED_COLUMN,
                          TAXONOMY_RANKS, UNASSIGNED_LABEL, interval_bounds, species_row_keys)

    counts = "estimated_counts" in species_df.columns
    relatives = IDENTICAL_COLUMN in species_df.columns
    is_placeholder = species_df["species"].isin([BELOW_DETECTION_LABEL, UNASSIGNED_LABEL])
    taxa = species_df[~is_placeholder]
    row_of_key = {tuple(r): i for i, r in zip(taxa.index, taxa[TAXONOMY_RANKS].itertuples(index=False, name=None))}
    reported_names = set(taxa["species"])
    summary = {str(k): v for k, v in read_summary.items()}
    posterior = lambda t, col: float(summary.get(t, {}).get(col, 0.0))
    total_alpha = sum(posterior(t, "alpha") for t in summary) if intervals else None
    # The row of every tax ID, in one pass over the taxonomy
    all_ids = [str(t) for group in groups for t in group["tax_ids"]]
    row_of_id = dict(zip(all_ids, (row_of_key[k] for k in species_row_keys(all_ids, taxonomy_df, binomial))))

    out = []
    for group in groups:
        tax_ids = [str(t) for t in group["tax_ids"]]
        rows = taxa.loc[sorted({row_of_id[t] for t in tax_ids})]
        rows = rows.iloc[np.lexsort((rows["species"].to_numpy(), -rows["abundance"].to_numpy(dtype=float)))]
        expected = sum(posterior(t, "expected_reads") for t in tax_ids)
        entry = {
            "members": "; ".join(rows["species"]),
            "n_members": len(rows),
            "abundance": float(rows["abundance"].astype(float).sum()),
            "lowest_common_rank": "",
            "lowest_common_taxon": "",
            "max_shared": group["max_shared"],
            SHARED_COLUMN: group["shared_reads"] / expected if expected > 0 else 0.0,
        }
        for rank in TAXONOMY_RANKS:
            values = set(rows[rank])
            if len(values) == 1 and "" not in values:
                entry["lowest_common_rank"], entry["lowest_common_taxon"] = rank, values.pop()
                break
        if counts:
            entry["estimated_counts"] = int(rows["estimated_counts"].astype(int).sum())
        if relatives:
            names = {n for v in rows[IDENTICAL_COLUMN] for n in str(v).split("; ") if n}
            entry[RELATIVES_COLUMN] = "; ".join(sorted(names - reported_names))
        if intervals:
            entry["alpha"] = sum(posterior(t, "alpha") for t in tax_ids)
        if len(rows) == 1:
            # The species row's own values (summed there in another order)
            entry[SHARED_COLUMN] = rows[SHARED_COLUMN].iloc[0] if SHARED_COLUMN in rows else entry[SHARED_COLUMN]
            if intervals:
                entry["species_ci"] = (rows["ci_lower"].iloc[0], rows["ci_upper"].iloc[0])
        out.append(entry)
    out.sort(key=lambda e: (-e["abundance"], e["members"]))
    for i, entry in enumerate(out, 1):
        entry["group"] = f"G{i}"

    df = pd.DataFrame(out) if out else pd.DataFrame(columns=GROUP_TABLE_COLUMNS)
    if intervals and out:
        lower, upper = interval_bounds(df["alpha"].to_numpy(dtype=float), total_alpha)
        for i, ci in df.get("species_ci", pd.Series(dtype=object)).items():
            if isinstance(ci, tuple):
                lower[i], upper[i] = ci
        df["ci_lower"], df["ci_upper"] = lower, upper
    placeholders = species_df[is_placeholder]
    df = pd.concat([df, pd.DataFrame({
        "group": placeholders["species"], "members": placeholders["species"], "abundance": 0.0,
        **({"estimated_counts": placeholders["estimated_counts"]} if counts else {}),
    })], ignore_index=True)

    columns = GROUP_TABLE_COLUMNS + (["estimated_counts"] if counts else []) + \
        ([RELATIVES_COLUMN] if relatives else []) + ["max_shared", SHARED_COLUMN] + \
        (INTERVAL_COLUMNS if intervals else [])
    for col in columns:
        if col not in df.columns:
            df[col] = np.nan
    df = df[columns]
    df["n_members"] = df["n_members"].astype("Int64")
    if counts:
        df["estimated_counts"] = df["estimated_counts"].astype("Int64")
    df.to_csv(path, sep="\t", index=False)
    return df
