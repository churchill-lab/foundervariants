"""The encodings — the part of the system everything else trusts.

`sdpn` is the only genotype actually stored; `sdp`, `alleles`, `ref` and the per-founder
bases are all derived from it on read. If these functions are wrong, every answer the
service gives is wrong, and no downstream test would notice because they all use the
same encoder. So these are checked against values computed by hand.
"""
from __future__ import annotations

import pytest

from foundervariants import encode as enc

A_J, C57, S129, NOD, NZO, CAST, PWK, WSB = range(8)


def test_founder_order_is_the_documented_one():
    # The A-H order is baked into both encodings and into every sdp/sdpn already stored.
    # Reordering this list silently reinterprets every database in existence.
    assert enc.FOUNDERS == [
        "A_J", "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ",
        "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ",
    ]
    assert len(enc.FOUNDERS) == 8


class TestSdpn:
    """Base-4 over A/C/G/T, A/J most-significant, 2 bits per founder."""

    def test_worked_example_from_the_docs(self):
        # ARCHITECTURE.md documents sdpn=43682 as "only PWK carries A, everyone else G".
        bases = enc.sdpn_to_bases(43682)
        assert bases[PWK] == "A"
        assert [b for i, b in enumerate(bases) if i != PWK] == ["G"] * 7

    def test_round_trip(self):
        for seq in ("AAAAAAAA", "TTTTTTTT", "ACGTACGT", "GGGGGGAG", "CCTTGGCC"):
            assert enc.sdpn_to_sequence(enc.bases_to_sdpn(list(seq))) == seq

    def test_bit_layout_matches_the_shift_formula(self):
        # founder i lives at bits 2*(7-i): the formula the SQL predicates use.
        sdpn = enc.bases_to_sdpn(list("ACGTACGT"))
        for i in range(8):
            assert (sdpn >> (2 * (7 - i))) & 3 == enc.founder_digit(sdpn, i)
            assert enc.founder_shift(i) == 2 * (7 - i)

    def test_range_is_16_bit(self):
        assert enc.bases_to_sdpn(list("AAAAAAAA")) == 0
        assert enc.bases_to_sdpn(list("TTTTTTTT")) == 65535

    def test_b6_position_carries_ref(self):
        # B6 is forced to REF at build time, so ref is recoverable from sdpn alone.
        sdpn = enc.bases_to_sdpn(list("TGGTGGGG"))
        assert enc.sdpn_to_bases(sdpn)[C57] == "G"


class TestSdp:
    """One bit per founder, A/J = bit 0, set when the founder differs from REF."""

    def test_only_pwk_differs(self):
        sdpn = enc.bases_to_sdpn(list("GGGGGGAG"))
        assert enc.sdp_from_sdpn(sdpn, "G") == 1 << PWK == 64

    def test_a_j_is_the_least_significant_bit(self):
        sdpn = enc.bases_to_sdpn(list("TGGGGGGG"))
        assert enc.sdp_from_sdpn(sdpn, "G") == 1

    def test_wsb_is_the_most_significant_bit(self):
        sdpn = enc.bases_to_sdpn(list("GGGGGGGT"))
        assert enc.sdp_from_sdpn(sdpn, "G") == 128

    def test_two_founders(self):
        sdpn = enc.bases_to_sdpn(list("TGGTGGGG"))
        assert enc.sdp_from_sdpn(sdpn, "G") == (1 << A_J) | (1 << NOD) == 9

    def test_b6_bit_is_always_clear(self):
        # B6 == REF by construction, so bit 1 can never be set.
        for seq in ("TGGTGGGG", "GGGGGGAG", "TTTTTTTT"):
            sdp = enc.sdp_from_sdpn(enc.bases_to_sdpn(list(seq)), seq[C57])
            assert not sdp & (1 << C57)

    def test_stays_in_the_rqtl2_polymorphic_range(self):
        # qtl2 treats 0 and 255 as monomorphic; a stored site should never be either.
        for seq in ("TGGTGGGG", "GGGGGGAG", "CCTTGGCC", "TTTTTTTA"):
            sdp = enc.sdp_from_sdpn(enc.bases_to_sdpn(list(seq)), seq[C57])
            assert 1 <= sdp <= 254


class TestConfidence:
    """Same bit trick as sdp, over the Sanger per-founder FI flags."""

    def test_all_high(self):
        assert enc.calc_sdp([1] * 8) == 255
        assert enc.conf_bits(255) == [1] * 8

    def test_one_low_founder_round_trips(self):
        bits = [1, 1, 1, 1, 1, 0, 1, 1]      # CAST low
        conf = enc.calc_sdp(bits)
        assert conf != 255
        assert enc.conf_bits(conf) == bits
        assert enc.conf_bits(conf)[CAST] == 0


class TestAllelesAndDisplay:
    def test_alleles_lists_ref_first(self):
        sdpn = enc.bases_to_sdpn(list("TGGTGGGG"))
        assert enc.alleles_string(sdpn, "G").startswith("G")
        assert "T" in enc.alleles_string(sdpn, "G")

    def test_multiallelic_alleles_carry_every_base(self):
        sdpn = enc.bases_to_sdpn(list("CCTTGGCC"))
        s = enc.alleles_string(sdpn, "C")
        assert set("CTG") <= set(s)

    def test_sdp_partition_splits_the_founders(self):
        # sdp=64 is PWK alone: the H-letter for PWK is G (7th letter, A-H).
        assert enc.sdp_to_partition(64) == "ABCDEFH:G"

    def test_sdpn_display_shows_more_groups_than_sdp_can(self):
        # THE point of sdpn: a 3-allele site is 3 groups, which one sdp bit cannot express.
        display = enc.sdpn_to_display(enc.bases_to_sdpn(list("CCTTGGCC")))
        assert len(display.split(":")) == 3

    def test_sdpn_display_biallelic_matches_the_sdp_partition(self):
        sdpn = enc.bases_to_sdpn(list("TGGTGGGG"))
        sdp = enc.sdp_from_sdpn(sdpn, "G")
        assert set(enc.sdpn_to_display(sdpn).split(":")) == set(enc.sdp_to_partition(sdp).split(":"))


class TestParseRegion:
    @pytest.mark.parametrize("text,expect", [
        ("19:5000000-5200000", ("19", 5_000_000, 5_200_000)),
        ("19:5.0-5.2Mb", ("19", 5_000_000, 5_200_000)),
        ("chr19:5.0-5.2Mb", ("19", 5_000_000, 5_200_000)),
        ("19:5.0-5.2M", ("19", 5_000_000, 5_200_000)),
        ("19:100-200kb", ("19", 100_000, 200_000)),
        ("19:500kb-2Mb", ("19", 500_000, 2_000_000)),   # each side names its own unit
        ("19:5.0Mb-5.2", ("19", 5_000_000, 5_200_000)),  # unit on the left instead
        ("X:1-1000", ("X", 1, 1000)),
        ("MT:1-16000", ("MT", 1, 16_000)),
    ])
    def test_accepted_forms(self, text, expect):
        r = enc.parse_region(text)
        assert (r.chrom, r.start, r.end) == expect

    def test_a_single_trailing_unit_applies_to_both_ends(self):
        """Regression: '19:5.0-5.2Mb' once parsed as 5bp..5.2Mb.

        The suffix was applied only to `end`, so the window silently started at the
        beginning of the chromosome. It returned a superset instead of raising, which
        is why it went unnoticed — and this is the exact form the docs advertise.
        """
        r = enc.parse_region("19:5.0-5.2Mb")
        assert r.start == 5_000_000, "trailing unit must scale the start too"
        assert r.end == 5_200_000
        assert enc.parse_region("19:100-200kb").start == 100_000

    @pytest.mark.parametrize("bad", ["", "19", "nonsense", "19:abc-def", "19:5.0"])
    def test_rejects_junk_with_valueerror(self, bad):
        # api.py maps ValueError to HTTP 422, so this is the contract behind that.
        with pytest.raises(ValueError):
            enc.parse_region(bad)


class TestFounderPredicates:
    """SQL fragments that do bit arithmetic on sdpn — the founder-pattern queries."""

    def test_differ_names_both_shifts(self):
        sql = enc.differ_predicate("CAST_EiJ", "PWK_PhJ")
        assert str(enc.founder_shift(CAST)) in sql
        assert str(enc.founder_shift(PWK)) in sql

    def test_private_mentions_the_founder_shift(self):
        assert str(enc.founder_shift(WSB)) in enc.private_predicate("WSB_EiJ")

    def test_predicates_are_valid_sqlite(self, con):
        # A predicate that does not parse would fail only at request time in production.
        for sql in (enc.differ_predicate("CAST_EiJ", "PWK_PhJ"),
                    enc.private_predicate("WSB_EiJ")):
            con.execute(f"SELECT count(*) FROM variants WHERE {sql}").fetchone()
