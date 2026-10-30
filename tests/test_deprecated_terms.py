"""The deprecated-consequence-term mapping.

Releases annotated by different VEP versions must land on ONE vocabulary, or the same
biological event becomes unfindable under the name the current docs use. VEP 75 emits
`initiator_codon_variant`; VEP 104 emits `start_lost`. See build.DEPRECATED_TERMS.
"""
from __future__ import annotations

import csv
import os

import pytest

from foundervariants.build import DEPRECATED_TERMS, canonical_term

# The real vocabulary lives one level up, outside the repo (a shared input). Tests that
# need it skip rather than fail when it is not there, so the suite stays self-contained.
SO_TSV = os.path.join(os.path.dirname(__file__), "..", "..", "data", "so_terms.tsv")
has_tsv = pytest.mark.skipif(not os.path.exists(SO_TSV),
                             reason="../data/so_terms.tsv not present")


def load_terms():
    with open(SO_TSV) as f:
        return {r["term"]: r for r in csv.DictReader(f, delimiter="\t")}


class TestMapping:
    def test_the_known_rename(self):
        assert canonical_term("initiator_codon_variant") == "start_lost"

    def test_current_terms_pass_through_unchanged(self):
        for t in ("start_lost", "missense_variant", "stop_gained", "intron_variant"):
            assert canonical_term(t) == t

    def test_unknown_terms_pass_through_unchanged(self):
        # A term we have never seen must not be silently rewritten.
        assert canonical_term("some_future_so_term") == "some_future_so_term"

    def test_mapping_is_not_chained(self):
        """No target may itself be deprecated — one pass must fully canonicalise."""
        for old, new in DEPRECATED_TERMS.items():
            assert new not in DEPRECATED_TERMS, f"{old} -> {new} -> ... is a chain"
            assert canonical_term(canonical_term(old)) == canonical_term(old)

    def test_nothing_maps_to_itself(self):
        for old, new in DEPRECATED_TERMS.items():
            assert old != new


@has_tsv
class TestAgainstTheRealVocabulary:
    def test_every_target_is_in_the_vocabulary(self):
        """A target missing from so_terms.tsv would score -1 and set no csq_mask bit."""
        terms = load_terms()
        for old, new in DEPRECATED_TERMS.items():
            assert new in terms, f"{old} maps to {new}, absent from so_terms.tsv"

    def test_no_deprecated_term_remains_in_the_vocabulary(self):
        """A deprecated term in the TSV is a dead filter entry: it can never match."""
        terms = load_terms()
        for old in DEPRECATED_TERMS:
            assert old not in terms, (
                f"{old} is both deprecated and still in so_terms.tsv — it would appear "
                "in the UI filter and always return zero rows")

    def test_start_lost_is_ranked_high(self):
        """Ensembl 75 ranked initiator_codon_variant 7th, right after stop_lost.

        That is the slot start_lost occupies now. If start_lost is ever demoted, the 143
        GRCm38 variants that map onto it would be under-ranked, which is the bug this
        whole mapping exists to fix.
        """
        terms = load_terms()
        sl, stop_lost = terms["start_lost"], terms["stop_lost"]
        assert sl["impact"] == "HIGH"
        assert int(sl["score"]) < int(stop_lost["score"])
        # and above the MODERATE block
        assert int(sl["score"]) > int(terms["missense_variant"]["score"])
        assert int(sl["score"]) > int(terms["start_retained_variant"]["score"])

    def test_bit_assignment_is_row_order_and_start_lost_is_stable(self):
        """csq_mask bits come from row order; start_lost must stay at bit 6.

        Every built database's masks depend on this. Removing initiator_codon_variant was
        only safe because it was the LAST row.
        """
        with open(SO_TSV) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        bits = {r["term"]: i for i, r in enumerate(rows)}
        assert bits["start_lost"] == 6

    def test_scores_are_unique_enough_to_rank(self):
        """Ties are tolerated by the builder, but a target must not tie its neighbours."""
        terms = load_terms()
        assert int(terms["start_lost"]["score"]) != int(terms["stop_lost"]["score"])
