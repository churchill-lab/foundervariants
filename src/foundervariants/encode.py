"""
Canonical encoding + query helpers for FounderVariants.

Design premise (validated on real data): for SNPs, (sdpn, ref) fully determines
each founder's base, so `alt`, `sdp`, per-founder alleles, and per-founder
confidence are DERIVED, not stored. This module is the single source of truth,
shared by the builder, the query layer, the API, and the R/qtl2 export.

Founder order (A..H) is fixed. Two encodings, opposite digit order (matches the
legacy code exactly, so output validates against fv.2021.snps.db3):

  SDPN : base-4 over A/C/G/T, A_J = MOST-significant digit  -> int 0..65535 (16-bit)
  SDP  : base-2 ref/alt,      A_J = LEAST-significant bit    -> int 0..255  (qtl2)
  conf : base-2 FI flags,     A_J = LEAST-significant bit    -> int 0..255
"""
from __future__ import annotations
import re
from dataclasses import dataclass

FOUNDERS: list[str] = [
    "A_J", "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ",
    "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ",
]
FOUNDER_LETTER = dict(zip(FOUNDERS, "ABCDEFGH"))
REFERENCE_STRAIN = "C57BL_6J"  # index 1: not in the VCF, forced to REF
BASES = ["A", "C", "G", "T"]
BASE_TO_DIGIT = {b: i for i, b in enumerate(BASES)}


# ---------------------------------------------------------------- SDPN (base-4)
def bases_to_sdpn(bases8: list[str]) -> int:
    """8 founder bases (A_J..WSB) -> sdpn. A_J is the most-significant digit."""
    return int("".join(str(BASE_TO_DIGIT[b]) for b in bases8), 4)


def sdpn_to_bases(sdpn: int) -> list[str]:
    """sdpn -> 8 founder bases (A_J..WSB)."""
    return [BASES[(sdpn >> (2 * (7 - i))) & 3] for i in range(8)]


def founder_digit(sdpn: int, i: int) -> int:
    """base-4 digit (0..3) for founder index i (0=A_J .. 7=WSB)."""
    return (sdpn >> (2 * (7 - i))) & 3


def founder_shift(i: int) -> int:
    """bit shift to extract founder i's 2-bit digit from sdpn (for SQL bit-ops)."""
    return 2 * (7 - i)


# ------------------------------------------------------------ SDP / confidence
def calc_sdp(bits8: list[int]) -> int:
    """base-2 with A_J as LSB (matches qtl2::calc_sdp and legacy calc_sdp)."""
    return int("".join(str(b) for b in reversed(bits8)), 2)


def sdp_from_sdpn(sdpn: int, ref: str) -> int:
    """Derive the qtl2/Broman SDP (ref-anchored: 1 = founder base != REF)."""
    return calc_sdp([0 if b == ref else 1 for b in sdpn_to_bases(sdpn)])


def conf_bits(confidence: int) -> list[int]:
    """confidence int -> 8 per-founder FI flags (A_J..WSB), A_J = LSB."""
    return [(confidence >> i) & 1 for i in range(8)]


# --------------------------------------------------------------------- display
def founder_alleles(sdpn: int) -> dict[str, str]:
    """{founder_name: base} decoded from sdpn."""
    return dict(zip(FOUNDERS, sdpn_to_bases(sdpn)))


def alleles_string(sdpn: int, ref: str) -> str:
    """R/qtl2-style 'alleles' = REF|alt1/alt2 among the 8 founders (ref-anchored)."""
    alts: list[str] = []
    for b in sdpn_to_bases(sdpn):
        if b != ref and b not in alts:
            alts.append(b)
    return f"{ref}|{'/'.join(alts)}" if alts else ref


def sdp_to_partition(sdp: int) -> str:
    """Human display of an SDP, e.g. 'ABEH:CDFG' (0-group : 1-group)."""
    zero = "".join(FOUNDER_LETTER[FOUNDERS[i]] for i in range(8) if not (sdp >> i) & 1)
    one = "".join(FOUNDER_LETTER[FOUNDERS[i]] for i in range(8) if (sdp >> i) & 1)
    return f"{zero}:{one}"


def sdpn_to_sequence(sdpn: int) -> str:
    """Raw base string across the 8 founders (A_J..WSB), e.g. 'AGCTTTGT'."""
    return "".join(sdpn_to_bases(sdpn))


def sdpn_to_display(sdpn: int) -> str:
    """Character/fill-string form of SDPN: founders grouped by shared base,
    largest group first, e.g. 'DEFH:BG:A:C' (analogous to sdp_to_partition)."""
    groups: dict[str, str] = {}
    for i, b in enumerate(sdpn_to_bases(sdpn)):
        groups[b] = groups.get(b, "") + FOUNDER_LETTER[FOUNDERS[i]]
    return ":".join(v for _, v in sorted(groups.items(), key=lambda kv: len(kv[1]), reverse=True))


# ---------------------------------------------------------------- region parse
_REGION = re.compile(
    r"^\s*(?:chr)?\s*(?P<chr>[0-9]{1,2}|X|Y|MT)\s*[:\-]\s*"
    r"(?P<start>[\d.]+)\s*(?P<su>mb|kb|m|k)?\s*[-:]\s*"
    r"(?P<end>[\d.]+)\s*(?P<eu>mb|kb|m|k)?\s*$",
    re.IGNORECASE,
)
_MULT = {"mb": 1_000_000, "m": 1_000_000, "kb": 1_000, "k": 1_000, None: 1, "": 1}


@dataclass
class Region:
    chrom: str
    start: int
    end: int


def parse_region(text: str) -> Region:
    m = _REGION.match(text)
    if not m:
        raise ValueError(
            f"cannot parse region {text!r} (try '19:5.0-5.2Mb' or '19:5000000-5200000')"
        )
    g = m.groupdict()
    # A unit named on one side applies to BOTH: "19:5.0-5.2Mb" means 5.0Mb-5.2Mb, not
    # 5bp-5.2Mb. Writing the suffix once, at the end, is the common form — and getting
    # this wrong silently widens the window to the start of the chromosome instead of
    # failing, so it returns a superset rather than an error.
    su, eu = (g["su"] or "").lower(), (g["eu"] or "").lower()
    su, eu = su or eu, eu or su
    s = float(g["start"]) * _MULT[su]
    e = float(g["end"]) * _MULT[eu]
    chrom = g["chr"].upper() if g["chr"].upper() in ("X", "Y", "MT") else g["chr"]
    return Region(chrom, int(round(s)), int(round(e)))


# --------------------------------------------- founder-pattern query -> SQL
def differ_predicate(founder_a: str, founder_b: str) -> str:
    """SQL predicate: two founders carry different bases."""
    ia, ib = FOUNDERS.index(founder_a), FOUNDERS.index(founder_b)
    return f"((sdpn>>{founder_shift(ia)})&3) <> ((sdpn>>{founder_shift(ib)})&3)"


def private_predicate(founder: str) -> str:
    """SQL predicate: this founder's base differs from all 7 others."""
    i = FOUNDERS.index(founder)
    parts = [
        f"((sdpn>>{founder_shift(i)})&3) <> ((sdpn>>{founder_shift(j)})&3)"
        for j in range(8) if j != i
    ]
    return "(" + " AND ".join(parts) + ")"
