/*
 * SDP / SDPN conversions for the /sdp page. Pure functions, no DOM.
 *
 * Mirrors src/foundervariants/encode.py exactly (tests/test_sdpcalc.py runs this file
 * under node and compares it with the Python on thousands of values). Founder order
 * A..H = A_J, C57BL_6J, 129S1_SvImJ, NOD_ShiLtJ, NZO_HlLtJ, CAST_EiJ, PWK_PhJ, WSB_EiJ.
 *
 *   SDPN: base-4 over A/C/G/T (A=0 C=1 G=2 T=3), A_J = MOST significant digit -> 0..65535
 *   SDP : base-2, 1 = base differs from the reference, A_J = LEAST significant bit -> 0..255
 *
 * The reference is C57BL/6J (B): the reference genome is B6, so in this database B's base
 * is always the reference base and B's SDP bit is always 0.
 */
"use strict";

const SDP_LETTERS = "ABCDEFGH";
const SDP_BASES = "ACGT";
const SDP_REF = 1;  // index of C57BL/6J

function basesToSdpn(bases) {
  let n = 0;
  for (const b of bases) n = n * 4 + SDP_BASES.indexOf(b);
  return n;
}

function sdpnToBases(n) {
  const out = [];
  for (let i = 0; i < 8; i++) out.push(SDP_BASES[(n >> (2 * (7 - i))) & 3]);
  return out;
}

// 1 where the founder's base differs from the reference (B's base); A_J = bit 0
function sdpFromBases(bases) {
  const ref = bases[SDP_REF];
  let sdp = 0;
  bases.forEach((b, i) => { if (b !== ref) sdp |= 1 << i; });
  return sdp;
}

// "ref group:alt group", e.g. "BCDEFGH:A" (encode.sdp_to_partition)
function sdpToPartition(sdp) {
  let zero = "", one = "";
  for (let i = 0; i < 8; i++) (((sdp >> i) & 1) ? (one += SDP_LETTERS[i]) : (zero += SDP_LETTERS[i]));
  return zero + ":" + one;
}

// founders grouped by shared base, largest group first; ties keep first-appearance order
// (encode.sdpn_to_display: a stable sort over an insertion-ordered dict)
function sdpnToDisplay(n) {
  const groups = new Map();
  sdpnToBases(n).forEach((b, i) => groups.set(b, (groups.get(b) || "") + SDP_LETTERS[i]));
  return [...groups.values()].sort((a, b) => b.length - a.length).join(":");
}

// R/qtl2-style "REF|alt1/alt2" (encode.alleles_string), alts in founder order
function allelesString(bases) {
  const ref = bases[SDP_REF], alts = [];
  for (const b of bases) if (b !== ref && !alts.includes(b)) alts.push(b);
  return alts.length ? `${ref}|${alts.join("/")}` : ref;
}

// ---- parsing what a person types; each returns {bases} or {error} ------------------

function parseSdpn(text) {
  const t = String(text).trim();
  if (!/^\d+$/.test(t) || +t > 65535) return {error: "SDPN is a whole number from 0 to 65535"};
  return {bases: sdpnToBases(+t)};
}

function parseSequence(text) {
  const t = String(text).trim().toUpperCase();
  if (!/^[ACGT]{8}$/.test(t)) return {error: "8 bases, one per founder A–H, each A, C, G or T"};
  return {bases: t.split("")};
}

// "A:BC:DEFGH": which founders share a base. It does not say WHICH base, so: B's group
// keeps the current reference base, and each other group keeps its first founder's
// current base when no earlier group has taken it, else the first unused base.
function parsePattern(text, current) {
  const groups = String(text).toUpperCase().replace(/\s+/g, "").split(":").filter(Boolean);
  const seen = groups.join("");
  if (!/^[A-H]+$/.test(seen) || seen.length !== 8 || new Set(seen).size !== 8)
    return {error: "use each founder letter A–H exactly once, groups separated by \":\""};
  if (groups.length > 4) return {error: "at most 4 groups: there are only 4 bases"};
  const out = current.slice(), used = new Set();
  const order = groups.slice().sort((a, b) => (b.includes("B") ? 1 : 0) - (a.includes("B") ? 1 : 0));
  for (const g of order) {
    const want = g.includes("B") ? current[SDP_REF] : current[SDP_LETTERS.indexOf(g[0])];
    const base = !used.has(want) ? want : [...SDP_BASES].find(b => !used.has(b));
    used.add(base);
    for (const l of g) out[SDP_LETTERS.indexOf(l)] = base;
  }
  return {bases: out};
}

// SDP number -> bases. Bit 0 founders take the reference base; bit 1 founders keep a
// non-reference base they already have, else take one (the first alt already in use, or
// the first base that isn't the reference). B is the reference, so an SDP with B's bit
// set is complemented: 255 - sdp splits the founders identically.
function parseSdp(text, current) {
  const t = String(text).trim();
  if (!/^\d+$/.test(t) || +t > 255) return {error: "SDP is a whole number from 0 to 255"};
  let sdp = +t, note = "";
  if (sdp & (1 << SDP_REF)) {
    note = `B (C57BL/6J) is the reference, so its bit is always 0: ${sdp} is the same split as ${255 - sdp}.`;
    sdp = 255 - sdp;
  }
  return {bases: basesFromSplit(sdp, current), note};
}

// "BCDEFGH:A": reference group, then the founders that differ
function parseSplit(text, current) {
  const parts = String(text).toUpperCase().replace(/\s+/g, "").split(":");
  const all = parts.join("");
  if (parts.length !== 2 || !/^[A-H]*$/.test(all) || all.length !== 8 || new Set(all).size !== 8)
    return {error: "two groups, reference:different, using each founder letter A–H once"};
  let sdp = 0;
  for (const l of parts[1]) sdp |= 1 << SDP_LETTERS.indexOf(l);
  let note = "";
  if (sdp & (1 << SDP_REF)) {
    note = "B (C57BL/6J) is the reference, so it belongs in the first group; the groups were swapped.";
    sdp = 255 - sdp;
  }
  return {bases: basesFromSplit(sdp, current), note};
}

function basesFromSplit(sdp, current) {
  const ref = current[SDP_REF];
  const alt = current.find((b, i) => ((sdp >> i) & 1) && b !== ref) || [...SDP_BASES].find(b => b !== ref);
  return current.map((b, i) => (((sdp >> i) & 1) ? (b !== ref ? b : alt) : ref));
}

// let node's require() load this file for the tests; browsers ignore it
if (typeof module !== "undefined") module.exports = {
  basesToSdpn, sdpnToBases, sdpFromBases, sdpToPartition, sdpnToDisplay, allelesString,
  parseSdpn, parseSequence, parsePattern, parseSdp, parseSplit,
};
