"""The R/qtl2 contract, asserted in Python so a schema change cannot break R silently.

`qtl2::create_variant_query_func(db)` is called with NO field overrides, which means the
column names below are load-bearing. R is not available in CI, so these tests encode what
R requires and replicate the query it actually issues:

    SELECT * FROM variants WHERE chr=='<chr>' AND pos>=<start_bp> AND pos<=<end_bp>

with start/end in Mbp multiplied by 1e6 by the query func, and pos divided by 1e6 on the
way out. Anything that breaks these assertions breaks every downstream scan1snps() call.
"""
from __future__ import annotations

import sqlite3

import pytest

# What create_variant_query_func() and top_snps() need, under exactly these names.
REQUIRED_COLUMNS = {"snp_id", "chr", "pos", "alleles", "sdp"}


@pytest.fixture
def raw(db_path):
    con = sqlite3.connect(db_path)
    yield con
    con.close()


def columns(con, table="variants"):
    return {r[1]: r[2].upper() for r in con.execute(f"PRAGMA table_info({table})")}


class TestSchemaContract:
    def test_the_table_is_called_variants(self, raw):
        names = {r[0] for r in raw.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        assert "variants" in names, "qtl2 defaults to table_name='variants'"

    def test_required_columns_exist_under_canonical_names(self, raw):
        assert set(columns(raw)) >= REQUIRED_COLUMNS

    def test_pos_is_an_integer_column(self, raw):
        # The query func multiplies Mbp by 1e6 and compares against pos, so pos must be
        # basepairs stored as an integer. Storing Mbp here would break every query.
        assert columns(raw)["pos"] == "INTEGER"

    def test_sdp_is_an_integer_column(self, raw):
        # index_snps() groups on sdp and snpprobs.cpp treats it as an integer bitfield.
        assert columns(raw)["sdp"] == "INTEGER"

    def test_chr_is_text(self, raw):
        # qtl2 compares chr to a character vector: chr=='19'.
        assert columns(raw)["chr"] == "TEXT"

    def test_chr_pos_index_exists(self, raw):
        idx = {r[0] for r in raw.execute(
            "SELECT name FROM sqlite_master WHERE type='index'")}
        assert "chr_pos" in idx

    def test_description_table_present(self, raw):
        row = raw.execute("SELECT description, genome_build FROM description").fetchone()
        assert row[0].startswith("FounderVariants")
        assert row[1]


class TestTheQueryRQtl2Issues:
    def _qv(self, con, chrom, start_mbp, end_mbp):
        """Replicate create_variant_query_func()'s SQL and its Mbp<->bp conversion."""
        rows = con.execute(
            "SELECT * FROM variants WHERE chr==? AND pos>=? AND pos<=?",
            (chrom, int(start_mbp * 1e6), int(end_mbp * 1e6))).fetchall()
        cols = [d[0] for d in con.execute("SELECT * FROM variants LIMIT 0").description]
        out = [dict(zip(cols, r, strict=True)) for r in rows]
        for r in out:
            r["pos"] = r["pos"] / 1e6          # what qtl2 returns to the caller
        return out

    def test_returns_rows_for_a_window(self, raw):
        res = self._qv(raw, "19", 0.001, 0.004)
        assert [r["snp_id"] for r in res] == ["rs001", "rs002", "rs003"]

    def test_pos_comes_back_in_mbp(self, raw):
        res = self._qv(raw, "19", 0.001, 0.004)
        assert res[0]["pos"] == pytest.approx(0.0015)

    def test_chromosome_is_respected(self, raw):
        assert self._qv(raw, "5", 0.0, 1.0)
        assert self._qv(raw, "19", 0.0, 0.0001) == []

    def test_empty_window_is_empty_not_an_error(self, raw):
        assert self._qv(raw, "19", 900.0, 901.0) == []

    def test_sdp_is_polymorphic_everywhere(self, raw):
        """qtl2 treats sdp 0 and 255 as monomorphic and index_snps() would drop them."""
        bad = raw.execute(
            "SELECT count(*) FROM variants WHERE sdp < 1 OR sdp > 254").fetchone()[0]
        assert bad == 0

    def test_alleles_is_populated_for_top_snps(self, raw):
        missing = raw.execute(
            "SELECT count(*) FROM variants WHERE alleles IS NULL OR alleles=''").fetchone()[0]
        assert missing == 0

    def test_snp_id_is_never_null(self, raw):
        n = raw.execute("SELECT count(*) FROM variants WHERE snp_id IS NULL").fetchone()[0]
        assert n == 0

    def test_index_snps_style_grouping_works(self, raw):
        """index_snps() partitions a window by distinct sdp; that must be meaningful."""
        res = self._qv(raw, "19", 0.0, 0.01)
        groups = {r["sdp"] for r in res}
        assert len(groups) > 1


class TestDdlIsTheSingleSourceOfTruth:
    def test_build_and_tests_use_the_same_ddl(self, raw):
        """The fixture is built from schema.DDL_TABLES, so this asserts the real schema."""
        from foundervariants import schema
        assert "snp_id TEXT" in schema.DDL_VARIANTS
        assert "chr_pos" in " ".join(schema.DDL_INDEXES)
        # Column order in the DDL is the order R sees from SELECT *.
        order = [r[1] for r in raw.execute("PRAGMA table_info(variants)")]
        assert order[:6] == ["snp_id", "chr", "pos", "ref", "alleles", "sdp"]
