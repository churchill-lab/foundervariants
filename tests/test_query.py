"""The read layer: selectors, filters, decoding, and the errors the API maps to HTTP."""
from __future__ import annotations

import pytest

from foundervariants import encode as enc
from foundervariants import query as Q

CAST, PWK, WSB = 5, 6, 7


def ids(res):
    return [v["variants_id"] for v in res["variants"]]


class TestConnect:
    def test_opens_read_only(self, db_path):
        import sqlite3
        con = Q.connect(db_path)
        with pytest.raises(sqlite3.OperationalError):
            con.execute("INSERT INTO variants VALUES (1,1,1,1,1,1,1,1,1,1,1,1)")
        con.close()

    def test_meta_is_a_flat_dict(self, con):
        m = Q.meta(con)
        assert m["release"] == "9999"
        assert m["grc"] == "GRCm39"
        assert m["founders"].split(",") == enc.FOUNDERS


class TestGenes:
    def test_prefix_search_is_case_insensitive(self, con):
        # lowercase query, symbol-ordered results
        assert [g["symbol"] for g in Q.genes_search(con, "pt")] == ["Pten", "Ptgs2"]

    def test_search_honours_limit(self, con):
        assert len(Q.genes_search(con, "Pt", limit=1)) == 1

    def test_exact_lookup(self, con):
        g = Q.gene(con, "Pten")
        assert (g["chr"], g["start"], g["end"]) == ("19", 1000, 5000)

    def test_exact_lookup_is_case_insensitive(self, con):
        assert Q.gene(con, "pTeN")["ensembl_id"] == "ENSMUSG00000013663"

    def test_missing_gene_is_none(self, con):
        assert Q.gene(con, "NoSuchGene") is None


class TestSelectors:
    def test_region(self, con):
        assert ids(Q.search_variants(con, region="19:1000-3000")) == ["rs001", "rs002"]

    def test_region_with_unit_suffix(self, con):
        # Guards the parse_region fix end-to-end, not just in the parser.
        assert ids(Q.search_variants(con, region="19:0-0.004Mb")) == \
               ids(Q.search_variants(con, region="19:0-4000"))

    def test_gene_resolves_to_its_span(self, con):
        # rs003 also sits in Pten but is low-confidence, so the default filter drops it.
        assert ids(Q.search_variants(con, gene_symbol="Pten")) == ["rs001", "rs002"]
        assert ids(Q.search_variants(con, gene_symbol="Pten", low_confidence=True)) == \
               ["rs001", "rs002", "rs003"]

    def test_gene_is_scoped_to_that_gene_only(self, con):
        # Actb is on chr5; nothing from chr19 may leak in.
        assert ids(Q.search_variants(con, gene_symbol="Actb")) == ["rs006", "5:600_G_A"]

    def test_rsid(self, con):
        assert ids(Q.search_variants(con, rsid="rs005")) == ["rs005"]

    def test_positional_id_works_as_an_rsid(self, con):
        assert ids(Q.search_variants(con, rsid="5:600_G_A")) == ["5:600_G_A"]

    def test_no_selector_raises_valueerror(self, con):
        # api.py -> HTTP 422
        with pytest.raises(ValueError):
            Q.search_variants(con)

    def test_unknown_gene_raises_lookuperror(self, con):
        # api.py -> HTTP 404
        with pytest.raises(LookupError):
            Q.search_variants(con, gene_symbol="Nope")

    def test_gene_wins_over_region_when_both_given(self, con):
        res = Q.search_variants(con, gene_symbol="Actb", region="19:1000-3000")
        assert ids(res) == ["rs006", "5:600_G_A"]   # chr5, i.e. the gene won

    def test_limit(self, con):
        assert len(Q.search_variants(con, region="19:0-99999", limit=2)["variants"]) == 2

    def test_results_are_position_ordered(self, con):
        pos = [v["pos"] for v in Q.search_variants(con, region="19:0-99999")["variants"]]
        assert pos == sorted(pos)


class TestConfidenceFilter:
    def test_high_confidence_only_by_default(self, con):
        # rs003 is the low-confidence one and must be absent.
        assert "rs003" not in ids(Q.search_variants(con, gene_symbol="Pten"))

    def test_low_confidence_can_be_included(self, con):
        assert "rs003" in ids(Q.search_variants(con, gene_symbol="Pten", low_confidence=True))

    def test_decoded_low_conf_founders_are_named(self, con):
        v = Q.search_variants(con, rsid="rs003", low_confidence=True)["variants"][0]
        assert v["high_confidence"] is False
        assert v["low_conf_founders"] == ["CAST_EiJ"]


class TestConsequenceFilter:
    def test_single_term(self, con):
        res = Q.search_variants(con, region="19:0-99999", consequence="stop_gained")
        assert ids(res) == ["rs002"]

    def test_multi_term_is_a_union(self, con):
        res = Q.search_variants(con, region="19:0-99999",
                                consequence="stop_gained,missense_variant")
        assert set(ids(res)) == {"rs001", "rs002", "rs005"}

    def test_matches_a_non_top_term(self, con):
        # rs001's top term is missense, but it also carries intron_variant. The mask
        # must match on any term, not only the one in the `consequence` column.
        res = Q.search_variants(con, rsid="rs001", consequence="intron_variant")
        assert ids(res) == ["rs001"]

    def test_unknown_term_raises_valueerror(self, con):
        with pytest.raises(ValueError, match="unknown consequence"):
            Q.search_variants(con, region="19:0-99999", consequence="not_a_term")

    def test_whitespace_and_empties_tolerated(self, con):
        res = Q.search_variants(con, region="19:0-99999", consequence=" stop_gained , ")
        assert ids(res) == ["rs002"]


class TestFounderPatterns:
    def test_private_to_wsb(self, con):
        res = Q.search_variants(con, region="19:0-99999", private_to="WSB_EiJ")
        assert ids(res) == ["rs005"]

    def test_private_to_finds_nothing_when_shared(self, con):
        res = Q.search_variants(con, region="19:2000-3000", private_to="WSB_EiJ")
        assert res["count"] == 0

    def test_differ(self, con):
        # rs004 bases, A/J..WSB, are C C T T G G C C:
        #   A_J=C  NOD=T   -> differ
        #   NZO=G  CAST=G  -> same
        assert ids(Q.search_variants(con, rsid="rs004", differ="A_J,NOD_ShiLtJ")) == ["rs004"]
        assert Q.search_variants(con, rsid="rs004", differ="NZO_HlLtJ,CAST_EiJ")["count"] == 0

    def test_differ_is_order_independent(self, con):
        a = Q.search_variants(con, rsid="rs004", differ="A_J,NOD_ShiLtJ")["count"]
        b = Q.search_variants(con, rsid="rs004", differ="NOD_ShiLtJ,A_J")["count"]
        assert a == b == 1

    def test_exact_sdpn(self, con):
        target = enc.bases_to_sdpn(list("TTTTTTTA"))
        assert ids(Q.search_variants(con, region="19:0-99999", sdpn=target)) == ["rs005"]

    def test_unknown_founder_raises_valueerror(self, con):
        with pytest.raises(ValueError, match="unknown founder"):
            Q.search_variants(con, region="19:0-99999", private_to="Nope_J")

    def test_differ_needs_exactly_two_founders(self, con):
        with pytest.raises(ValueError):
            Q.search_variants(con, region="19:0-99999", differ="A_J")

    def test_filters_compose(self, con):
        res = Q.search_variants(con, region="19:0-99999",
                                private_to="WSB_EiJ", consequence="missense_variant")
        assert ids(res) == ["rs005"]


class TestDecode:
    def test_founder_bases_are_decoded_from_sdpn(self, con):
        v = Q.search_variants(con, rsid="rs001")["variants"][0]
        assert v["founders"]["PWK_PhJ"] == "A"
        assert v["founders"]["C57BL_6J"] == "G"
        assert set(v["founders"]) == set(enc.FOUNDERS)

    def test_derived_sdp_matches_the_stored_column(self, con):
        # sdp is stored (R/qtl2 needs it) AND derivable. They must agree.
        for vid in ("rs001", "rs002", "rs004", "rs005"):
            v = Q.search_variants(con, rsid=vid)["variants"][0]
            assert v["sdp"] == enc.sdp_from_sdpn(v["sdpn"], v["ref"])

    def test_multiallelic_shows_three_groups(self, con):
        v = Q.search_variants(con, rsid="rs004")["variants"][0]
        assert len(v["sdpn_display"].split(":")) == 3

    def test_count_matches_the_rows_returned(self, con):
        res = Q.search_variants(con, region="19:0-99999")
        assert res["count"] == len(res["variants"])


class TestSoTerms:
    def test_ordered_by_severity(self, con):
        scores = [t["score"] for t in Q.so_terms(con)]
        assert scores == sorted(scores, reverse=True)

    def test_bits_are_distinct(self, con):
        bits = [t["bit"] for t in Q.so_terms(con)]
        assert len(bits) == len(set(bits))
