"""
The UI's deep-link contract (the /linking page, web/static/LINKING.md): the page URL
other sites link to.

The parser lives in web/index.html between `// <deeplink>` markers and touches no DOM,
so these tests cut that block out and run it under node — the shipped code, not a
Python re-implementation of it. Skipped (not failed) where node isn't installed.
"""
from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess

import pytest

PAGE = os.path.join(os.path.dirname(__file__), "..", "src", "foundervariants", "web", "index.html")
RELS = [{"release": "1505", "grc": "GRCm38"}, {"release": "2021", "grc": "GRCm39"}]


def _node() -> str | None:
    if os.environ.get("FV_NODE"):
        return os.environ["FV_NODE"]
    # nvm installs node behind a shell function, invisible to a non-interactive PATH
    return shutil.which("node") or (sorted(glob.glob(os.path.expanduser("~/.nvm/versions/node/*/bin/node"))) or [None])[-1]


NODE = _node()
pytestmark = pytest.mark.skipif(NODE is None, reason="node not installed")


def _block() -> str:
    with open(PAGE, encoding="utf-8") as f:
        src = f.read()
    start, end = src.index("// <deeplink>"), src.index("// </deeplink>")
    return src[start:end]


def js(*exprs: str):
    """Evaluate JS expressions against the page's deep-link block; returns their JSON values."""
    harness = (_block() + f"\nconst RELS={json.dumps(RELS)};\n"
               "console.log(JSON.stringify(JSON.parse(process.argv[1]).map(e=>eval(e))));")
    out = subprocess.run([NODE, "-e", harness, json.dumps(list(exprs))],
                         capture_output=True, text=True, check=True)
    res = json.loads(out.stdout)
    return res[0] if len(res) == 1 else res


def parse(search: str) -> dict:
    return js(f"parseLink({json.dumps(search)}, RELS)")


class TestGuideExamples:
    """Every "Try it" link in the served guide must be a link the page actually accepts."""

    def links(self):
        import re
        with open(os.path.join(os.path.dirname(PAGE), "static", "LINKING.md"), encoding="utf-8") as f:
            md = f.read()
        return re.findall(r"\[([^\]]+)\]\((\?[^)]+)\)", md)

    def test_there_are_examples(self):
        assert len(self.links()) >= 5

    def test_each_example_parses(self):
        for text, href in self.links():
            p = parse(href)
            assert p["q"], (text, href, p)
            if "mistakes" in text:
                assert p["warnings"], href          # the deliberately broken one
            else:
                assert p["release"] and p["warnings"] == [], (href, p)


class TestNoLink:
    def test_a_bare_page_url_changes_nothing(self):
        p = parse("")
        assert p == {"release": None, "q": "", "consequence": [], "private_to": "", "sdpn": "",
                     "low_confidence": False, "mark": [], "limit": 2000, "warnings": []}


class TestRelease:
    @pytest.mark.parametrize("given,want", [
        ("1505", "1505"), ("2021", "2021"),
        ("GRCm38", "1505"), ("grcm39", "2021"),   # assembly names, any case
        ("mm10", "1505"), ("MM39", "2021"),       # UCSC names
    ])
    def test_release_by_id_assembly_or_ucsc_name(self, given, want):
        assert parse(f"?release={given}")["release"] == want

    def test_unknown_release_warns_and_falls_back(self):
        p = parse("?release=GRCh38&q=Pten")
        assert p["release"] is None and p["q"] == "Pten"
        assert p["warnings"] == ['unknown release "GRCh38"']


class TestDefaultRelease:
    """With no release= in the link, the page opens the newest assembly (first in byNewest)."""

    def order(self, rels):
        return js(f"byNewest({json.dumps(rels)}).map(r=>r.release)")

    def test_grcm39_is_the_default_over_grcm38(self):
        # the server lists 1505 (GRCm38) first; that order must not decide the default
        assert self.order(RELS)[0] == "2021"
        assert self.order(list(reversed(RELS)))[0] == "2021"

    def test_a_newer_assembly_takes_over_without_code_changes(self):
        assert self.order(RELS + [{"release": "2510", "grc": "GRCm40"}]) == ["2510", "2021", "1505"]

    def test_same_assembly_prefers_the_higher_release_numerically(self):
        rels = [{"release": "999", "grc": "GRCm39"}, {"release": "2021", "grc": "GRCm39"}]
        assert self.order(rels) == ["2021", "999"]

    def test_unknown_assembly_sorts_last(self):
        assert self.order([{"release": "x", "grc": ""}] + RELS)[-1] == "x"

    def test_explicit_release_still_wins(self):
        assert parse("?release=GRCm38&q=Pten")["release"] == "1505"


class TestWhatToSearch:
    def test_q_is_used_verbatim(self):
        assert parse("?q=19:5.0-5.2Mb")["q"] == "19:5.0-5.2Mb"

    @pytest.mark.parametrize("key,val", [("region", "11:100-102Mb"), ("gene", "Pten"), ("rsid", "rs30847235")])
    def test_api_style_selectors(self, key, val):
        assert parse(f"?{key}={val}")["q"] == val

    def test_q_wins_over_the_others(self):
        assert parse("?gene=Pten&q=Brca1&region=1:1-2&chr=2&pos=5")["q"] == "Brca1"

    def test_chr_pos_defaults_to_50kb_each_side_and_marks_pos(self):
        p = parse("?chr=10&pos=128350502")
        assert p["q"] == "10:128300502-128400502"
        assert p["mark"] == ["128350502"]
        assert p["warnings"] == []

    def test_chr_prefix_commas_and_window_units(self):
        p = parse("?chr=chr10&pos=128,350,502&window=1Mb")
        assert p["q"] == "10:127350502-129350502"

    def test_window_in_kb(self):
        assert parse("?chr=X&pos=70500000&window=25kb")["q"] == "X:70475000-70525000"

    def test_window_never_runs_below_position_1(self):
        assert parse("?chr=1&pos=3000&window=10kb")["q"] == "1:1-13000"

    def test_chr_start_end_with_units(self):
        p = parse("?chr=11&start=100.5Mb&end=102Mb")
        assert p["q"] == "11:100500000-102000000" and p["mark"] == []

    def test_start_after_end_is_rejected(self):
        p = parse("?chr=11&start=5Mb&end=4Mb")
        assert p["q"] == "" and len(p["warnings"]) == 1

    def test_chr_alone_is_not_enough(self):
        p = parse("?chr=11")
        assert p["q"] == "" and "chr needs" in p["warnings"][0]

    def test_bad_pos_is_rejected(self):
        p = parse("?chr=11&pos=abc")
        assert p["q"] == "" and p["mark"] == [] and len(p["warnings"]) == 1


class TestMark:
    def test_positions_and_rsids_normalised_and_deduped(self):
        p = parse("?q=10:1-2Mb&mark=RS123, 1500000 ,1.5Mb,rs123")
        assert p["mark"] == ["rs123", "1500000"]

    def test_mark_is_comma_separated_so_no_thousands_separators(self):
        # documented: "1,500,000" in mark is three positions, not one
        assert parse("?q=1:1-2Mb&mark=1,500,000")["mark"] == ["1", "500", "0"]

    def test_mark_merges_with_pos(self):
        assert parse("?chr=1&pos=500000&mark=500000,600000")["mark"] == ["500000", "600000"]

    def test_unusable_mark_warns(self):
        p = parse("?q=Pten&mark=peak")
        assert p["mark"] == [] and p["warnings"] == ['can\'t mark "peak"']


class TestFilters:
    def test_consequence_list(self):
        assert parse("?q=Pten&consequence=missense_variant, stop_gained,")["consequence"] == \
            ["missense_variant", "stop_gained"]

    @pytest.mark.parametrize("given,want", [
        ("CAST", "CAST_EiJ"), ("cast_eij", "CAST_EiJ"), ("A/J", "A_J"), ("a_j", "A_J"),
        ("B6", "C57BL_6J"), ("129", "129S1_SvImJ"), ("WSB_EiJ", "WSB_EiJ"),
    ])
    def test_private_to_accepts_short_and_full_names(self, given, want):
        assert parse(f"?q=Pten&private_to={given.replace('/', '%2F')}")["private_to"] == want

    def test_unknown_founder_warns(self):
        p = parse("?q=Pten&private_to=BALB")
        assert p["private_to"] == "" and p["warnings"] == ['unknown founder "BALB"']

    def test_sdpn(self):
        assert parse("?q=Pten&sdpn=43658")["sdpn"] == "43658"
        p = parse("?q=Pten&sdpn=4x")
        assert p["sdpn"] == "" and len(p["warnings"]) == 1

    @pytest.mark.parametrize("val,want", [("true", True), ("1", True), ("YES", True), ("on", True),
                                          ("false", False), ("0", False), ("", False)])
    def test_low_confidence(self, val, want):
        assert parse(f"?q=Pten&low_confidence={val}")["low_confidence"] is want


class TestLimit:
    def test_limit(self):
        assert parse("?q=Pten&limit=5000")["limit"] == 5000

    def test_limit_is_capped(self):
        p = parse("?q=Pten&limit=999999")
        assert p["limit"] == 25000 and p["warnings"] == ["limit capped at 25000"]

    @pytest.mark.parametrize("bad", ["0", "-5", "abc", "1.5"])
    def test_bad_limit_keeps_the_default(self, bad):
        p = parse(f"?q=Pten&limit={bad}")
        assert p["limit"] == 2000 and len(p["warnings"]) == 1


class TestCanonicalURL:
    """What the address bar shows after a search must parse back to the same state."""

    def test_round_trip(self):
        state, back = js(
            "parseLink('?release=GRCm38&chr=10&pos=128350502&window=1Mb&consequence=missense_variant"
            "&private_to=CAST&sdpn=7&low_confidence=1&mark=rs1&limit=9000', RELS)",
            "(s=>parseLink(linkParams(s).toString(), RELS))(parseLink('?release=GRCm38&chr=10&pos=128350502"
            "&window=1Mb&consequence=missense_variant&private_to=CAST&sdpn=7&low_confidence=1&mark=rs1"
            "&limit=9000', RELS))",
        )
        assert back == state

    def test_canonical_form_uses_q_and_omits_defaults(self):
        qs = js("linkParams(parseLink('?release=mm10&chr=10&pos=128350502', RELS)).toString()")
        assert qs == "release=1505&q=10%3A128300502-128400502&mark=128350502"
