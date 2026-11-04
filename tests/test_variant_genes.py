"""Every gene per SNP (`variant_genes`, schema compact-v2).

`variants` keeps one gene per SNP: the gene of the most severe consequence. These tests
pin the side table that holds the rest, from the build's per-gene masks through the query
layer and API, and that adding it left `variants`, the R/qtl2 contract, untouched.
"""
from __future__ import annotations

import csv
import io
import json

import pytest
from conftest import NO_SYMBOL_GENE, VARIANTS, make_db
from fastapi.testclient import TestClient

from foundervariants import build, schema
from foundervariants import query as Q

PTEN, PTGS2, ACTB = "ENSMUSG00000013663", "ENSMUSG00000000001", "ENSMUSG00000000002"
GM42418, SNHG1, AAA1 = "ENSMUSG00000098178", "ENSMUSG00000088000", "ENSMUSG00000077777"


@pytest.fixture
def con(db_path):
    c = Q.connect(db_path)
    yield c
    c.close()


def by_id(res):
    return {v["variants_id"]: v for v in res["variants"]}


def all_variants(con):
    return by_id(Q.search_variants(con, region="19:0-99999", low_confidence=True)) | \
        by_id(Q.search_variants(con, region="5:0-99999", low_confidence=True))


class TestGenesPerVariant:
    def test_overlapping_genes_are_all_listed_most_severe_first(self, con):
        v = all_variants(con)["rs001"]
        assert [g["gene_id"] for g in v["genes"]] == [PTEN, GM42418]
        pten, gm = v["genes"]
        assert (pten["symbol"], pten["consequence"], pten["impact"]) == \
            ("Pten", "missense_variant", "MODERATE")
        assert pten["consequences"] == ["missense_variant", "intron_variant"]   # by severity
        assert (gm["symbol"], gm["consequence"], gm["near"]) == ("Gm42418", "intron_variant", False)

    def test_a_flank_only_gene_is_near_not_in(self, con):
        v = all_variants(con)["rs002"]
        assert [(g["symbol"], g["near"]) for g in v["genes"]] == [("Pten", False), ("Snhg1", True)]
        assert v["genes"][1]["consequence"] == "upstream_gene_variant"

    def test_a_tie_puts_the_stored_gene_first(self, con):
        # Aaa1 sorts before Ptgs2 and has the same top term; genes[0] must still be the
        # gene in variants.ensembl_gene, so it agrees with gene_id / max_csq.
        v = all_variants(con)["rs004"]
        assert [g["gene_id"] for g in v["genes"]] == [PTGS2, AAA1]

    def test_first_gene_always_agrees_with_the_stored_summary(self, con):
        for v in all_variants(con).values():
            if v["gene_id"]:
                assert v["genes"][0]["gene_id"] == v["gene_id"]
                assert v["genes"][0]["consequence"] == v["max_csq"]

    def test_gene_symbol_of_the_stored_gene(self, con):
        vs = all_variants(con)
        assert vs["rs001"]["gene_symbol"] == "Pten"
        assert vs["5:600_G_A"]["gene_symbol"] is None

    def test_a_gene_missing_from_the_genes_table_keeps_its_id(self, con):
        v = all_variants(con)["rs006"]
        assert {g["gene_id"]: g["symbol"] for g in v["genes"]} == {ACTB: "Actb", NO_SYMBOL_GENE: None}

    def test_every_consequence_at_the_snp_most_severe_first(self, con):
        vs = all_variants(con)
        assert vs["rs001"]["consequences"] == ["missense_variant", "intron_variant"]
        assert vs["rs002"]["consequences"] == ["stop_gained", "upstream_gene_variant"]
        assert vs["5:600_G_A"]["consequences"] == []
        for v in vs.values():
            assert (v["consequences"] or [None])[0] == v["max_csq"]
            assert "_csq_mask" not in v   # internal; never leaks into a response

    def test_a_snp_with_no_gene_has_an_empty_list(self, con):
        assert all_variants(con)["5:600_G_A"]["genes"] == []


class TestOlderFiles:
    """A compact-v1 file (no variant_genes) is still served, with the same response shape."""

    @pytest.fixture
    def v1(self, tmp_path):
        c = Q.connect(make_db(str(tmp_path / "fv.9999.snps.sqlite"), variant_genes=False))
        yield c
        c.close()

    def test_the_stored_gene_only(self, v1):
        v = all_variants(v1)["rs001"]
        assert [(g["gene_id"], g["symbol"], g["consequence"], g["near"]) for g in v["genes"]] == \
            [(PTEN, "Pten", "missense_variant", False)]

    def test_consequences_too(self, v1):
        assert all_variants(v1)["rs001"]["consequences"] == ["missense_variant", "intron_variant"]

    def test_no_gene_no_entry(self, v1):
        assert all_variants(v1)["5:600_G_A"]["genes"] == []


class TestGeneResolution:
    def test_by_ensembl_id_in_any_case(self, con):
        assert Q.gene(con, "ensmusg00000013663")["symbol"] == "Pten"

    def test_a_symbol_with_two_loci_is_ambiguous(self, con):
        with pytest.raises(Q.AmbiguousGeneError) as ei:
            Q.gene(con, "U6")
        assert ei.value.total == 3
        assert [(c["chr"], c["ensembl_id"]) for c in ei.value.candidates] == \
            [("5", "ENSMUSG00000064001"), ("19", "ENSMUSG00000064002"),
             ("X", "ENSMUSG00000064003")]   # 1..19 numerically, then X

    def test_an_ambiguous_gene_is_still_a_lookup_error(self, con):
        # callers that only know LookupError keep working
        with pytest.raises(LookupError):
            Q.gene(con, "U6")

    def test_an_exact_case_match_wins(self, con):
        assert Q.gene(con, "Mir1a")["ensembl_id"] == "ENSMUSG00000070001"
        assert Q.gene(con, "MIR1A")["ensembl_id"] == "ENSMUSG00000070002"

    def test_case_only_differences_are_ambiguous_without_an_exact_match(self, con):
        with pytest.raises(Q.AmbiguousGeneError):
            Q.gene(con, "mir1a")

    def test_unknown(self, con):
        assert Q.gene(con, "Nope") is None
        assert Q.gene(con, "ENSMUSG00000000000") is None

    def test_search_returns_the_resolved_gene(self, con):
        res = Q.search_variants(con, gene_symbol="Gm42418")
        assert res["gene"]["ensembl_id"] == GM42418
        assert [v["variants_id"] for v in res["variants"]] == ["rs001"]

    def test_search_by_ensembl_id(self, con):
        res = Q.search_variants(con, gene_symbol="ENSMUSG00000064002")
        assert res["gene"]["chr"] == "19"

    def test_a_region_search_has_no_gene(self, con):
        assert "gene" not in Q.search_variants(con, region="19:0-99999")


class TestApi:
    @pytest.fixture
    def client(self, db_path, monkeypatch):
        import os
        monkeypatch.setenv("FV_DB_DIR", os.path.dirname(db_path))
        from foundervariants.api import create_app
        with TestClient(create_app()) as c:
            yield c

    def test_ambiguous_gene_is_409_with_candidates(self, client):
        r = client.get("/9999/variants?gene=U6")
        assert r.status_code == 409
        d = r.json()["detail"]
        assert d["symbol"] == "U6" and d["total"] == 3
        assert {c["ensembl_id"] for c in d["candidates"]} == \
            {"ENSMUSG00000064001", "ENSMUSG00000064002", "ENSMUSG00000064003"}
        assert "U6" in d["message"]

    def test_one_of_its_loci_by_id(self, client):
        r = client.get("/9999/variants?gene=ENSMUSG00000064002")
        assert r.status_code == 200
        assert r.json()["data"]["gene"]["chr"] == "19"

    def test_json_has_genes(self, client):
        data = client.get("/9999/variants?region=19:1000-3000").json()["data"]
        v = {x["variants_id"]: x for x in data["variants"]}
        assert [g["symbol"] for g in v["rs001"]["genes"]] == ["Pten", "Gm42418"]

    def test_flat_export_lists_every_gene(self, client):
        r = client.get("/9999/variants?region=19:1000-3000&format=csv")
        rows = {x["variants_id"]: x for x in csv.DictReader(io.StringIO(r.text))}
        assert rows["rs001"]["gene_symbol"] == "Pten"
        assert rows["rs001"]["genes"] == \
            f"Pten:{PTEN}:missense_variant|Gm42418:{GM42418}:intron_variant"
        assert rows["rs002"]["genes"].endswith(f"|~Snhg1:{SNHG1}:upstream_gene_variant")
        assert rows["rs002"]["consequences"] == "stop_gained|upstream_gene_variant"

    def test_ndjson_export_is_one_full_variant_per_line(self, client):
        r = client.get("/9999/variants?region=19:1000-3000&format=ndjson")
        assert r.status_code == 200
        assert r.headers["content-type"].startswith("application/x-ndjson")
        assert 'filename="foundervariants_9999.ndjson"' in r.headers["content-disposition"]
        assert r.text.endswith("\n")
        lines = r.text.splitlines()
        rows = [json.loads(x) for x in lines]   # every line parses on its own
        assert [v["variants_id"] for v in rows] == ["rs001", "rs002"]
        assert [g["symbol"] for g in rows[0]["genes"]] == ["Pten", "Gm42418"]   # nesting intact
        assert rows == client.get("/9999/variants?region=19:1000-3000&format=json").json()

    def test_ndjson_of_nothing_is_empty(self, client):
        r = client.get("/9999/variants?region=19:99000-99999&format=ndjson")
        assert r.status_code == 200 and r.text == ""

    def test_json_body_is_compact_and_complete(self, client):
        r = client.get("/9999/variants?region=19:0-99999")
        assert r.headers["content-type"] == "application/json"
        assert ", " not in r.text[:200]          # compact separators
        body = r.json()
        assert body["meta"] == {"release": "9999", "count": body["data"]["count"]}


class TestRqtl2ContractUnchanged:
    def test_variants_columns_are_the_v1_columns(self, con):
        cols = [r[1] for r in con.execute("PRAGMA table_info(variants)")]
        assert cols == ["snp_id", "chr", "pos", "ref", "alleles", "sdp", "sdpn", "confidence",
                        "ensembl_gene", "consequence", "csq_mask", "type"]

    def test_select_star_from_variants_is_unaffected(self, con):
        # what R/qtl2's query function runs: one row per SNP, no join, no extra columns
        rows = con.execute("SELECT * FROM variants WHERE chr=='19' AND pos>=0 AND pos<=99999").fetchall()
        assert len(rows) == sum(1 for v in VARIANTS if v[1] == "19")
        assert all(len(r) == 12 for r in rows)

    def test_variant_genes_is_a_separate_table(self):
        assert "variant_genes" not in schema.DDL_VARIANTS
        assert "chr_pos" in " ".join(schema.DDL_INDEXES)


class _Rec:
    """Just enough of a cyvcf2 record for parse_chrom."""

    def __init__(self, pos, ref, alts, gts, csq, rid=None):
        self.POS, self.REF, self.ALT, self.ID = pos, ref, alts, rid
        self.genotypes = gts
        self.INFO = {"CSQ": csq}

    def format(self, name):
        return [[1]] * len(self.genotypes) if name == "FI" else None


class TestBuildPerGeneMasks:
    """parse_chrom on a stub VCF: the per-gene masks, from the same CSQ loop as csq_mask."""

    TERMBIT = {"missense_variant": 0, "intron_variant": 1, "upstream_gene_variant": 2,
               "intergenic_variant": 3, "synonymous_variant": 4}
    SCORE = {"missense_variant": 29, "intron_variant": 14, "upstream_gene_variant": 11,
             "intergenic_variant": 1, "synonymous_variant": 15}

    def run(self, rec):
        # 7 samples: every founder except C57BL/6J (the reference, never a sample)
        fidx = [0, None, 1, 2, 3, 4, 5, 6]
        stats: dict = {}
        out = list(build.parse_chrom(lambda chrom: [rec], "19", fidx, (0, 1, 2), self.SCORE,
                                     self.TERMBIT, True, False, stats))
        return out, stats

    def test_per_gene_masks(self):
        ref_gt, alt_gt = [0, 0, False], [1, 1, False]
        gts = [ref_gt] * 5 + [alt_gt] + [ref_gt]           # PWK (6th sample) carries T
        csq = ",".join([
            "T|missense_variant|GENE_A",
            "T|intron_variant|GENE_B",
            "T|intron_variant&upstream_gene_variant|GENE_A",
            "T|upstream_gene_variant|GENE_C",
            "C|synonymous_variant|GENE_D",                  # an allele no founder carries
            "T|intergenic_variant|",                        # no gene
        ])
        (row, genes), = self.run(_Rec(100, "G", ["T", "C"], gts, csq, "rs1"))[0]
        mask = sum(1 << self.TERMBIT[t] for t in
                   ("missense_variant", "intron_variant", "upstream_gene_variant", "intergenic_variant"))
        assert row[8:11] == ("GENE_A", "missense_variant", mask)
        assert genes == [("19", 100, "GENE_A", 0b111), ("19", 100, "GENE_B", 0b010),
                         ("19", 100, "GENE_C", 0b100)]   # sorted by gene; D ignored; no "" row

    def test_stats_count_gene_rows(self):
        gts = [[0, 0, False]] * 6 + [[1, 1, False]]
        out, stats = self.run(_Rec(5, "A", ["G"], gts, "G|intron_variant|X,G|intron_variant|Y"))
        assert stats["gene_rows"] == 2 and stats["kept"] == 1
        assert [g[2] for g in out[0][1]] == ["X", "Y"]

    def test_a_variants_row_is_unchanged_in_shape(self):
        gts = [[0, 0, False]] * 6 + [[1, 1, False]]
        (row, _), = self.run(_Rec(5, "A", ["G"], gts, "G|intron_variant|X"))[0]
        assert len(row) == 12 and row[0] == "19:5_A_G" and row[11] == "snp"
