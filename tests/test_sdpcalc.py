"""
The /sdp calculator's conversions (web/static/sdpcalc.js) must agree exactly with
encode.py, the code that built the databases. Runs the shipped JS under node and compares
it with the Python on random values and, where the databases are present, on real rows.
Skipped (not failed) where node isn't installed.
"""
from __future__ import annotations

import glob
import json
import os
import random
import shutil
import sqlite3
import subprocess

import pytest

from foundervariants import encode as enc

ROOT = os.path.join(os.path.dirname(__file__), "..")
CALC = os.path.abspath(os.path.join(ROOT, "src", "foundervariants", "web", "static", "sdpcalc.js"))
DBS = os.path.join(ROOT, "dbs")


def _node() -> str | None:
    if os.environ.get("FV_NODE"):
        return os.environ["FV_NODE"]
    # nvm installs node behind a shell function, invisible to a non-interactive PATH
    return shutil.which("node") or (sorted(glob.glob(os.path.expanduser("~/.nvm/versions/node/*/bin/node"))) or [None])[-1]


NODE = _node()
pytestmark = pytest.mark.skipif(NODE is None, reason="node not installed")


def js(expr: str, data=None):
    """Evaluate one JS expression with the calculator's functions in scope and `data` bound."""
    prog = (f"const m=require({json.dumps(CALC)});Object.assign(globalThis,m);"
            "const data=JSON.parse(process.argv[1]);"
            f"console.log(JSON.stringify(({expr})));")
    out = subprocess.run([NODE, "-e", prog, json.dumps(data)], capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


EDGES = [0, 1, 4, 5461, 21843, 21789, 43690, 65535, 16383, 49152]
VALUES = EDGES + random.Random(20261021).sample(range(65536), 3000)


class TestMatchesEncodePy:
    def test_decoding_every_field(self):
        got = js("""data.map(n => { const b = sdpnToBases(n), s = sdpFromBases(b);
                     return [b.join(''), sdpnToDisplay(n), s, sdpToPartition(s), allelesString(b),
                             basesToSdpn(b)]; })""", VALUES)
        for n, (seq, disp, sdp, part, alleles, back) in zip(VALUES, got, strict=True):
            ref = enc.sdpn_to_bases(n)[1]          # B6 is the reference
            assert seq == enc.sdpn_to_sequence(n), n
            assert disp == enc.sdpn_to_display(n), n
            assert sdp == enc.sdp_from_sdpn(n, ref), n
            assert part == enc.sdp_to_partition(sdp), n
            assert alleles == enc.alleles_string(n, ref), n
            assert back == n

    @pytest.mark.parametrize("release", ["1505", "2021"])
    def test_real_rows(self, release):
        path = os.path.join(DBS, f"fv.{release}.snps.sqlite")
        if not os.path.exists(path):
            pytest.skip("database not present")
        con = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
        rows = con.execute("SELECT sdpn, ref, sdp, alleles FROM variants "
                           "WHERE rowid % 9973 = 0 LIMIT 4000").fetchall()
        con.close()
        assert len(rows) > 1000
        got = js("data.map(n => { const b = sdpnToBases(n); return [b[1], sdpFromBases(b), allelesString(b)]; })",
                 [r[0] for r in rows])
        for (sdpn, ref, sdp, alleles), (b_base, js_sdp, js_alleles) in zip(rows, got, strict=True):
            assert b_base == ref, sdpn          # B6's base IS the reference base
            assert js_sdp == sdp, sdpn          # the stored SDP, recomputed in the browser
            assert js_alleles == alleles, sdpn


class TestEditingAnyFieldRoundTrips:
    """Typing back what the calculator shows must leave every founder's base unchanged."""

    def test_pattern_sdp_split_sequence_and_sdpn(self):
        got = js("""data.map(n => { const b = sdpnToBases(n), s = sdpFromBases(b);
                     return [parsePattern(sdpnToDisplay(n), b).bases.join(''),
                             parseSdp(String(s), b).bases.join(''),
                             parseSplit(sdpToPartition(s), b).bases.join(''),
                             parseSequence(b.join('').toLowerCase()).bases.join(''),
                             parseSdpn(String(n)).bases.join('')]; })""", VALUES)
        for n, row in zip(VALUES, got, strict=True):
            assert set(row) == {enc.sdpn_to_sequence(n)}, (n, row)


class TestEditingChangesThings:
    def test_pattern_assigns_bases_to_new_groups(self):
        b = list("ACCCCCCC")                       # rs29358506: A/J differs
        r = js("parsePattern('ABCDEF:G:H', data)", b)["bases"]
        assert r[1] == "C" and len(set(r[:6])) == 1 and len({r[0], r[6], r[7]}) == 3
        assert js("sdpnToDisplay(basesToSdpn(data))", r) == "ABCDEF:G:H"

    def test_sdp_with_the_reference_bit_is_complemented(self):
        r = js("parseSdp('2', data)", list("CCCCCCCC"))
        assert "always 0" in r["note"]
        assert js("sdpFromBases(data)", r["bases"]) == 253

    def test_split_with_b_in_the_second_group_is_swapped(self):
        r = js("parseSplit('ACDEFGH:B', data)", list("CCCCCCCC"))
        assert "swapped" in r["note"] and js("sdpToPartition(sdpFromBases(data))", r["bases"]) == "B:ACDEFGH"

    def test_sdp_keeps_existing_alt_bases(self):
        # G=A, H=T (SDP 192); re-entering 192 must not collapse them to one alt base
        r = js("parseSdp('192', data)", list("CCCCCCAT"))
        assert "".join(r["bases"]) == "CCCCCCAT"

    @pytest.mark.parametrize("call", [
        "parseSdpn('65536')", "parseSdpn('-1')", "parseSdpn('1.5')", "parseSdpn('')",
        "parseSequence('ACGT')", "parseSequence('ACGTACGN')",
        "parsePattern('AB:CD', data)", "parsePattern('AB:CDEFGHA', data)",
        "parsePattern('A:B:C:D:EFGH', data)", "parsePattern('AB:CDEFGX', data)",
        "parseSdp('256', data)", "parseSdp('x', data)",
        "parseSplit('ABCDEFGH', data)", "parseSplit('A:B:CDEFGH', data)",
    ])
    def test_bad_input_is_an_error_not_a_change(self, call):
        r = js(call, list("ACCCCCCC"))
        assert "error" in r and "bases" not in r
