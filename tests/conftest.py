"""Test fixtures.

Every fixture builds a tiny database from `schema.py` itself — the same DDL a real
build uses — so the tests exercise the production schema without a 24 GB VCF, without
`cyvcf2`, and without a network. A fixture database is a few kilobytes and builds in
milliseconds.

The variants below are hand-computed, not captured from a run: `sdpn` values are written
out in base 4 so the expected founder bases are readable, and the expected `sdp` follows
from them. That way a test failure means the code changed, not that a golden file drifted.
"""
from __future__ import annotations

import sqlite3

import pytest

from foundervariants import encode as enc
from foundervariants import schema

# The 40-term vocabulary is loaded from a TSV in production; here it is a 4-term stand-in.
# `bit` is the term's index, exactly as build._load_so assigns it from row order.
SO_TERMS = [
    # term, description, accession, display, impact, score, bit
    ("stop_gained", "premature stop", "SO:0001587", "Stop gained", "HIGH", 36, 0),
    ("missense_variant", "changes an aa", "SO:0001583", "Missense variant", "MODERATE", 22, 1),
    ("synonymous_variant", "no aa change", "SO:0001819", "Synonymous variant", "LOW", 15, 2),
    ("intron_variant", "in an intron", "SO:0001627", "Intron variant", "MODIFIER", 8, 3),
    # A flank term: a gene with only this is "near", not "in" (variant_genes)
    ("upstream_gene_variant", "5' of a gene", "SO:0001631", "Upstream gene variant", "MODIFIER", 5, 4),
]
BIT = {t[0]: t[6] for t in SO_TERMS}

GENES = [
    # ensembl_id, version, symbol, name, chrom, start, end, strand
    ("ENSMUSG00000013663", "1", "Pten", "phosphatase and tensin homolog", "19", 1000, 5000, 1),
    ("ENSMUSG00000000001", "2", "Ptgs2", "prostaglandin synthase 2", "19", 8000, 9000, -1),
    ("ENSMUSG00000000002", "1", "Actb", "actin beta", "5", 100, 900, 1),
    # Overlaps Pten: rs001 is in both genes
    ("ENSMUSG00000098178", "1", "Gm42418", "predicted gene", "19", 1200, 1800, -1),
    # Just past rs002: VEP calls rs002 upstream of it (near, not in)
    ("ENSMUSG00000088000", "1", "Snhg1", "small nucleolar RNA host gene 1", "19", 2600, 2900, 1),
    # Overlaps Ptgs2 at rs004 with the SAME top term: a tie, which Ptgs2 must win
    ("ENSMUSG00000077777", "1", "Aaa1", "tie partner", "19", 8400, 8700, 1),
    # One symbol, two loci (like U6): ambiguous
    ("ENSMUSG00000064001", "1", "U6", "U6 spliceosomal RNA", "5", 50, 60, 1),
    ("ENSMUSG00000064002", "1", "U6", "U6 spliceosomal RNA", "19", 7000, 7100, 1),
    ("ENSMUSG00000064003", "1", "U6", "U6 spliceosomal RNA", "X", 10, 20, 1),   # X sorts last
    # Symbols differing only in case: an exact-case match wins, otherwise ambiguous
    ("ENSMUSG00000070001", "1", "Mir1a", "microRNA", "5", 700, 800, 1),
    ("ENSMUSG00000070002", "1", "MIR1A", "microRNA", "5", 900, 950, 1),
]
NO_SYMBOL_GENE = "ENSMUSG00000199999"   # annotated by VEP, absent from the genes table


def _sdpn(bases: str) -> int:
    """8 founder bases (A/J..WSB) -> sdpn, via the production encoder."""
    return enc.bases_to_sdpn(list(bases))


# (snp_id, chr, pos, ref, founder bases A/J..WSB, confidence, gene, top term, terms)
VARIANTS = [
    # All founders REF except PWK -> sdp = 1<<6 = 64. This is the ARCHITECTURE.md example.
    ("rs001", "19", 1500, "G", "GGGGGGAG", 255, "ENSMUSG00000013663",
     "missense_variant", ["missense_variant", "intron_variant"]),
    # A/J and NOD differ -> sdp = 1<<0 | 1<<3 = 9
    ("rs002", "19", 2500, "G", "TGGTGGGG", 255, "ENSMUSG00000013663",
     "stop_gained", ["stop_gained", "upstream_gene_variant"]),
    # Low confidence: CAST's FI flag is 0 -> confidence != 255
    ("rs003", "19", 3500, "A", "ACAAAAAA", 223, "ENSMUSG00000013663",
     "synonymous_variant", ["synonymous_variant"]),
    # Multiallelic: three distinct bases among the founders. SDP cannot express this.
    ("rs004", "19", 8500, "C", "CCTTGGCC", 255, "ENSMUSG00000000001",
     "intron_variant", ["intron_variant"]),
    # WSB alone differs -> private to WSB. sdp = 1<<7 = 128
    ("rs005", "19", 8600, "T", "TTTTTTTA", 255, "ENSMUSG00000000001",
     "missense_variant", ["missense_variant"]),
    # A different chromosome, so region/gene scoping is actually tested.
    ("rs006", "5", 500, "G", "GGGGGGGC", 255, "ENSMUSG00000000002",
     "intron_variant", ["intron_variant"]),
    # No rsID: the positional id form build.py falls back to.
    ("5:600_G_A", "5", 600, "G", "AGGGGGGG", 255, None, None, []),
]

# Every gene per SNP, as build.py writes variant_genes: (chr, pos) -> {gene: terms}.
# Each SNP's stored gene (VARIANTS) is the gene of its top term, and every gene's terms
# are within the SNP's terms: the two invariants the build guarantees.
VARIANT_GENES = {
    ("19", 1500): {"ENSMUSG00000013663": ["missense_variant", "intron_variant"],
                   "ENSMUSG00000098178": ["intron_variant"]},
    ("19", 2500): {"ENSMUSG00000013663": ["stop_gained"],
                   "ENSMUSG00000088000": ["upstream_gene_variant"]},
    ("19", 3500): {"ENSMUSG00000013663": ["synonymous_variant"]},
    ("19", 8500): {"ENSMUSG00000000001": ["intron_variant"],
                   "ENSMUSG00000077777": ["intron_variant"]},
    ("19", 8600): {"ENSMUSG00000000001": ["missense_variant"]},
    ("5", 500): {"ENSMUSG00000000002": ["intron_variant"], NO_SYMBOL_GENE: ["intron_variant"]},
}


def _mask(terms) -> int:
    m = 0
    for t in terms:
        m |= 1 << BIT[t]
    return m


META = {
    "schema_version": schema.SCHEMA_VERSION,
    "release": "9999",
    "name": "REL-9999",
    "variant": "SNP",
    "grc": "GRCm39",
    "created": "2026-01-01T00:00:00+00:00",
    "founders": ",".join(enc.FOUNDERS),
    "gene_annotation_file": "test.genes.tsv",
    "gene_annotation_version": "Ensembl TEST",
    "source_vcf": "test.vcf.gz",
    "chroms": "19,5",
    "n_genes": str(len(GENES)),
    "n_variants": str(len(VARIANTS)),
    "n_consequences": "9",
    "n_variant_genes": str(sum(len(g) for g in VARIANT_GENES.values())),
    "built_with": "tests.conftest",
    "served_from": "single-sqlite",
    "csq_model": schema.CSQ_MODEL,
}


def make_db(path: str, release: str = "9999", variant_genes: bool = True) -> str:
    """Create a fixture database at `path` using the production DDL.

    `variant_genes=False` builds a schema compact-v1 file (before variant_genes existed),
    which the server must still serve.
    """
    con = sqlite3.connect(path)
    for ddl in schema.DDL_TABLES:
        if variant_genes or "variant_genes" not in ddl:
            con.execute(ddl)
    if variant_genes:
        con.executemany("INSERT INTO variant_genes VALUES (?,?,?,?)", [
            (c, p, g, _mask(terms)) for (c, p), genes in VARIANT_GENES.items()
            for g, terms in genes.items()])

    for snp_id, chrom, pos, ref, bases, conf, gene, top, terms in VARIANTS:
        sdpn = _sdpn(bases)
        mask = 0
        for t in terms:
            mask |= 1 << BIT[t]
        con.execute("INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", (
            snp_id, chrom, pos, ref,
            enc.alleles_string(sdpn, ref),     # derived, exactly as build.py derives it
            enc.sdp_from_sdpn(sdpn, ref),      # derived
            sdpn, conf, gene, top, mask, "snp"))

    con.executemany("INSERT INTO genes VALUES (?,?,?,?,?,?,?,?)", GENES)
    con.executemany("INSERT INTO so_terms VALUES (?,?,?,?,?,?,?)", SO_TERMS)
    meta = dict(META, release=release, name=f"REL-{release}")
    con.executemany("INSERT INTO meta VALUES (?,?)", list(meta.items()))
    con.execute("INSERT INTO description VALUES (?,?,?)",
                (f"FounderVariants REL-{release}", meta["grc"], meta["created"]))
    for ddl in schema.DDL_INDEXES:
        con.execute(ddl)
    con.commit()
    con.close()
    return path


@pytest.fixture
def db_factory():
    """The builder itself, for tests that need a custom directory layout."""
    return make_db


@pytest.fixture
def db_path(tmp_path):
    """One fixture database, named the way discovery expects."""
    return make_db(str(tmp_path / "fv.9999.snps.sqlite"))


@pytest.fixture
def con(db_path):
    """A read-only connection, opened the way the API opens them."""
    from foundervariants import query
    c = query.connect(db_path)
    yield c
    c.close()


@pytest.fixture
def db_dir(tmp_path):
    """A directory holding TWO releases, so multi-release discovery is exercised."""
    d = tmp_path / "dbs"
    d.mkdir()
    make_db(str(d / "fv.9999.snps.sqlite"), release="9999")
    make_db(str(d / "fv.8888.snps.sqlite"), release="8888")
    return str(d)


@pytest.fixture
def client(db_dir, monkeypatch):
    """A TestClient over a freshly created app bound to db_dir."""
    from fastapi.testclient import TestClient

    from foundervariants import api
    monkeypatch.setenv("FV_DB_DIR", db_dir)
    with TestClient(api.create_app()) as c:
        yield c
