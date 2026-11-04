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

SCHEMA_VERSION = "compact-v2"   # v2 added variant_genes; `variants` is unchanged since v1
CSQ_MODEL = "per-variant bitmask + per-gene bitmask (variant_genes)"

# One row per SNP. Columns through `sdp` are the R/qtl2 contract (see module docstring);
# `sdpn` + `confidence` are the canonical genotype from which ref, alleles, sdp and the
# per-founder bases/confidence are all derived on read (see encode.py).
DDL_VARIANTS = """
CREATE TABLE variants(
        snp_id TEXT, chr TEXT, pos INTEGER, ref TEXT, alleles TEXT, sdp INTEGER,
        sdpn INTEGER, confidence INTEGER, ensembl_gene TEXT, consequence TEXT,
        csq_mask INTEGER, type TEXT)
"""

# Every gene VEP annotated at a SNP, one row per (SNP, gene), with that gene's own
# csq_mask (same bits as variants.csq_mask). `variants` keeps only the gene of the single
# most severe consequence, which hides genes that overlap: on GRCm39, 2.8% of SNPs sit in
# more than one gene body, 10.3% counting VEP's 5 kb up/downstream flanks.
#
# - Keyed by (chr, pos): unique in `variants` for every release built.
# - A gene's top consequence is derivable from its mask (so_terms scores are unique), so
#   it isn't stored. "Near only" (just up/downstream bits) is derivable too.
# - Gene-less CSQ entries (intergenic, regulatory) have no row; their bits are only in
#   variants.csq_mask.
# - A separate table, so `SELECT * FROM variants` (what R/qtl2 runs) is untouched.
DDL_VARIANT_GENES = """
CREATE TABLE variant_genes(chr TEXT, pos INTEGER, ensembl_gene TEXT, csq_mask INTEGER,
        PRIMARY KEY(chr, pos, ensembl_gene)) WITHOUT ROWID
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

DDL_TABLES = (DDL_VARIANTS, DDL_VARIANT_GENES, DDL_GENES, DDL_SO_TERMS, DDL_META,
              DDL_DESCRIPTION)

# Built after the bulk insert. chr_pos serves both R/qtl2's range query and ours.
DDL_INDEXES = (
    "CREATE INDEX chr_pos ON variants(chr,pos)",
    "CREATE INDEX ix_snp ON variants(snp_id)",
    "CREATE INDEX gi_sym ON genes(symbol COLLATE NOCASE)",
    "CREATE INDEX gi_loc ON genes(chromosome,start_position,end_position)",
    # v2: symbols for variant_genes IDs. Unique: one row per gene in every release.
    "CREATE UNIQUE INDEX gi_id ON genes(ensembl_id)",
)
