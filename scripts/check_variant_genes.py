#!/usr/bin/env python3
"""Check variant_genes against variants in a built database, at full scale.

    python3 -I scripts/check_variant_genes.py dbs_verify/fv.2021.snps.sqlite [...]

Read-only, stdlib only. Streams both tables in (chr, pos) order and merges them, so
memory stays flat on 36M-row files. Exits non-zero if any invariant fails:

  1. every variant_genes row belongs to a variant (no orphans);
  2. each gene's csq_mask is within its variant's csq_mask;
  3. a variant's stored gene (variants.ensembl_gene) has a variant_genes row;
  4. that gene's most severe term (derived from its mask) is variants.consequence;
  5. meta n_variant_genes is the row count.

It also prints how many genes SNPs have, how many SNPs have nearby-only genes, and how
many gene IDs have no symbol in the genes table.
"""
from __future__ import annotations

import itertools
import sqlite3
import sys
from collections import Counter

FLANK = ("upstream_gene_variant", "downstream_gene_variant")


def check(path: str) -> bool:
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    terms = con.execute("SELECT bit, term, score FROM so_terms").fetchall()
    bit_of = {t: b for b, t, _ in terms}
    score = {b: s for b, _, s in terms}
    flank = sum(1 << bit_of[t] for t in FLANK if t in bit_of)

    def top(mask: int) -> int:
        return max((b for b in score if mask >> b & 1), key=score.get)

    vg = con.cursor().execute(
        "SELECT chr, pos, ensembl_gene, csq_mask FROM variant_genes ORDER BY chr, pos, ensembl_gene")
    groups = itertools.groupby(vg, key=lambda r: (r[0], r[1]))
    nxt = next(groups, None)

    n_var = n_rows = orphans = not_subset = missing = mismatch = 0
    per, near_only = Counter(), 0
    for chrom, pos, gene, csq, mask in con.cursor().execute(
            "SELECT chr, pos, ensembl_gene, consequence, csq_mask FROM variants ORDER BY chr, pos"):
        n_var += 1
        while nxt is not None and nxt[0] < (chrom, pos):   # gene rows before this SNP: orphans
            orphans += sum(1 for _ in nxt[1])
            nxt = next(groups, None)
        genes = {}
        if nxt is not None and nxt[0] == (chrom, pos):
            genes = {g: m for _, _, g, m in nxt[1]}
            nxt = next(groups, None)
        n_rows += len(genes)
        per[min(len(genes), 4)] += 1
        if genes and all(m & ~flank == 0 for m in genes.values()):
            near_only += 1
        not_subset += sum(1 for m in genes.values() if m & ~mask)
        if gene is not None:
            if gene not in genes:
                missing += 1
            elif bit_of.get(csq) != top(genes[gene]):
                mismatch += 1
    while nxt is not None:
        orphans += sum(1 for _ in nxt[1])
        nxt = next(groups, None)
    n_rows += orphans

    meta = dict(con.execute("SELECT meta_name, meta_value FROM meta"))
    no_symbol = con.execute(
        "SELECT count(*) FROM (SELECT DISTINCT ensembl_gene FROM variant_genes) "
        "WHERE ensembl_gene NOT IN (SELECT ensembl_id FROM genes)").fetchone()[0]
    n_genes = con.execute("SELECT count(DISTINCT ensembl_gene) FROM variant_genes").fetchone()[0]

    checks = [
        ("variant_genes rows with no variant", orphans),
        ("gene masks outside their variant's mask", not_subset),
        ("stored gene with no variant_genes row", missing),
        ("stored gene whose top term != consequence", mismatch),
        ("meta n_variant_genes != rows", int(meta.get("n_variant_genes", -1)) != n_rows),
    ]
    print(f"{path}")
    print(f"  variants {n_var:,}   variant_genes {n_rows:,}   distinct genes {n_genes:,}")
    print("  genes per SNP: " + "  ".join(f"{k}{'+' if k == 4 else ''}: {v:,}" for k, v in sorted(per.items())))
    multi = sum(v for k, v in per.items() if k >= 2)
    print(f"  SNPs with 2+ genes: {multi:,} ({100 * multi / n_var:.1f}%)   nearby-only: {near_only:,}")
    print(f"  gene IDs with no symbol in genes: {no_symbol:,}")
    for name, bad in checks:
        print(f"  {'ok  ' if not bad else 'FAIL'} {name}: {int(bad):,}")
    return not any(bad for _, bad in checks)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(0 if all([check(p) for p in sys.argv[1:]]) else 1)
