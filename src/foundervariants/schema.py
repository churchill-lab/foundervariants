"""DDL for the one artifact this project produces.

`fv.<release>.snps.sqlite` is simultaneously the served store and the R/qtl2 variant
database, so this is the single source of truth for both. There is no intermediate
format and no export step.

**The `variants` table is an R/qtl2 contract.** `qtl2::create_variant_query_func()`
requires `snp_id`, `chr`, `pos` (INTEGER basepairs), `sdp` (INTEGER) and — for
`top_snps()` — `alleles`, under exactly those names, indexed by `(chr, pos)`. Those
columns and that index are not ours to rename or reorder; everything after `sdp` is
pass-through payload that R/qtl2 carries along untouched.
"""
from __future__ import annotations

SCHEMA_VERSION = "compact-v1"

# One row per SNP. Columns through `sdp` are the R/qtl2 contract (see module docstring);
# `sdpn` + `confidence` are the canonical genotype from which ref, alleles, sdp and the
# per-founder bases/confidence are all derived on read (see encode.py).
DDL_VARIANTS = """
CREATE TABLE variants(
        snp_id TEXT, chr TEXT, pos INTEGER, ref TEXT, alleles TEXT, sdp INTEGER,
        sdpn INTEGER, confidence INTEGER, ensembl_gene TEXT, consequence TEXT,
        csq_mask INTEGER, type TEXT)
"""

# Self-contained gene resolution: no runtime dependency on the sibling ensimpl service.
DDL_GENES = """
CREATE TABLE genes(ensembl_id TEXT, ensembl_version TEXT, symbol TEXT,
        name TEXT, chromosome TEXT, start_position INTEGER, end_position INTEGER, strand INTEGER)
"""

# `bit` is the term's position in every variant's csq_mask, assigned by row order in
# so_terms.tsv. Stored here so a database always carries the mapping its masks were
# built with.
# Kept on one line so the stored DDL text stays byte-identical to databases built by the
# previous two-step pipeline (SQLite records CREATE statements verbatim).
DDL_SO_TERMS = (
    "CREATE TABLE so_terms(term TEXT, description TEXT, accession TEXT, display TEXT,"
    " impact TEXT, score INTEGER, bit INTEGER)"
)

DDL_META = "CREATE TABLE meta(meta_name TEXT, meta_value TEXT)"

# R/qtl2 reads this table for its own provenance display.
DDL_DESCRIPTION = "CREATE TABLE description(description TEXT, genome_build TEXT, date_created TEXT)"

DDL_TABLES = (DDL_VARIANTS, DDL_GENES, DDL_SO_TERMS, DDL_META, DDL_DESCRIPTION)

# Built after the bulk insert. chr_pos serves both R/qtl2's range query and ours.
DDL_INDEXES = (
    "CREATE INDEX chr_pos ON variants(chr,pos)",
    "CREATE INDEX ix_snp ON variants(snp_id)",
    "CREATE INDEX gi_sym ON genes(symbol COLLATE NOCASE)",
    "CREATE INDEX gi_loc ON genes(chromosome,start_position,end_position)",
)
