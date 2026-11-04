"""
Read side: pure functions over a read-only SQLite connection, returning plain
dicts/lists. Usable as a library with only the stdlib: `from foundervariants import query`.

The single compact SQLite is both the served store and the R/qtl2 file. Each
returned variant decodes per-founder bases, sdp, alleles, and per-founder
confidence ON READ from (sdpn, confidence). Consequence filtering uses the
per-variant `csq_mask` bitmask (no transcript-level table). Every gene at a SNP, each
with its own consequences, comes from `variant_genes` (see attach_genes).
"""
from __future__ import annotations
import sqlite3
from typing import Any

from . import encode as enc

# columns read from the variants table (R/qtl2 names + our extras)
_VAR_COLS = ["snp_id", "chr", "pos", "ref", "sdpn", "confidence", "consequence", "ensembl_gene",
             "csq_mask"]


def connect(path: str) -> sqlite3.Connection:
    """Open a read-only, immutable connection (no locking, safe to fan out)."""
    return sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)


def meta(con) -> dict[str, str]:
    return dict(con.execute("SELECT meta_name, meta_value FROM meta").fetchall())


def so_terms(con) -> list[dict]:
    cols = ["term", "description", "accession", "display", "impact", "score", "bit"]
    rows = con.execute(f"SELECT {','.join(cols)} FROM so_terms ORDER BY score DESC").fetchall()
    return [dict(zip(cols, r)) for r in rows]


def _term_bits(con) -> dict[str, int]:
    return dict(con.execute("SELECT term, bit FROM so_terms").fetchall())


_GENE_SQL = "SELECT ensembl_id,symbol,chromosome,start_position,end_position,strand FROM genes"
_GENE_KEYS = ["ensembl_id", "symbol", "chr", "start", "end", "strand"]
_GENE_ORDER = ("ORDER BY symbol, CASE chromosome WHEN 'X' THEN 20 WHEN 'Y' THEN 21 WHEN 'MT' THEN 22"
               " ELSE CAST(chromosome AS INTEGER) END, start_position")   # 1..19, X, Y, MT
MAX_CANDIDATES = 50   # loci listed in an ambiguous-gene error (U6 alone names 486)


class AmbiguousGeneError(LookupError):
    """A symbol names more than one locus. Carries the candidates so a caller can offer them."""

    def __init__(self, symbol: str, candidates: list[dict], total: int):
        super().__init__(f"gene {symbol!r} names {total} loci; search by Ensembl ID instead")
        self.symbol, self.candidates, self.total = symbol, candidates, total


def genes_search(con, q: str, limit: int = 25) -> list[dict]:
    rows = con.execute(f"{_GENE_SQL} WHERE symbol LIKE ? {_GENE_ORDER} LIMIT ?",
                       [q + "%", limit]).fetchall()
    return [dict(zip(_GENE_KEYS, r, strict=True)) for r in rows]


def gene(con, key: str) -> dict | None:
    """Resolve an Ensembl gene ID or a symbol to one gene.

    Symbols match case-insensitively, but an exact-case match wins over the others. Raises
    AmbiguousGeneError when a symbol still names more than one locus: picking one silently
    would show the wrong place.
    """
    key = key.strip()
    if key.upper().startswith("ENSMUSG"):
        r = con.execute(f"{_GENE_SQL} WHERE ensembl_id = ?", [key.upper()]).fetchone()
        return dict(zip(_GENE_KEYS, r, strict=True)) if r else None
    rows = con.execute(f"{_GENE_SQL} WHERE symbol = ? COLLATE NOCASE {_GENE_ORDER}",
                       [key]).fetchall()
    exact = [r for r in rows if r[1] == key]
    rows = exact or rows
    if len(rows) > 1:
        raise AmbiguousGeneError(key, [dict(zip(_GENE_KEYS, r, strict=True)) for r in rows[:MAX_CANDIDATES]],
                            len(rows))
    return dict(zip(_GENE_KEYS, rows[0], strict=True)) if rows else None


# Terms that only place a SNP NEAR a gene (VEP's 5 kb flanks), not in it.
_FLANK_TERMS = ("upstream_gene_variant", "downstream_gene_variant")


def _has_table(con, name: str) -> bool:
    return con.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
                       [name]).fetchone() is not None


def attach_genes(con, variants: list[dict]) -> None:
    """Add `consequences`, `genes` and `gene_symbol` to each variant.

    `consequences` is every SO term at the SNP (its csq_mask), most severe first, so
    `consequences[0]` is `max_csq`. `genes` is every gene VEP annotated at the SNP.

    Each entry is one gene with its own most severe consequence, all its terms, and
    `near` (true when the SNP is only up/downstream of it). Ordered most severe first; on a
    tie the gene stored in `variants` (`gene_id`) comes first, so `genes[0]` always agrees
    with `max_csq` / `gene_id`.

    Databases built before variant_genes existed (schema compact-v1) get the one stored
    gene, so the response shape doesn't depend on the file being served.
    """
    if not variants:
        return
    terms = con.execute("SELECT bit, term, score, impact FROM so_terms").fetchall()
    bit_term = {b: (t, s, imp) for b, t, s, imp in terms}
    flank = sum(1 << b for b, t, _, _ in terms if t in _FLANK_TERMS)
    decoded: dict[int, tuple] = {}   # mask -> (names, top term, impact, score); masks repeat a lot

    def decode_mask(mask: int) -> tuple:
        d = decoded.get(mask)
        if d is None:
            hits = sorted((bit_term[b] for b in bit_term if mask >> b & 1), key=lambda x: -x[1])
            d = decoded[mask] = ([t for t, _, _ in hits], *(hits[0] if hits else (None, -1, None)))
        return d

    for v in variants:
        mask = v.pop("_csq_mask", None)
        if mask is not None:
            v["consequences"] = decode_mask(mask)[0]
    symbols: dict[str, str | None] = {}
    by_pos: dict[tuple[str, int], list[tuple[str, int]]] = {}
    if _has_table(con, "variant_genes"):
        keys = sorted({(v["chr"], v["pos"]) for v in variants})
        for chrom in sorted({c for c, _ in keys}):
            pos = [p for c, p in keys if c == chrom]
            for i in range(0, len(pos), 5000):   # PK seeks; chunked under SQLite's variable limit
                chunk = pos[i:i + 5000]
                sql = ("SELECT vg.pos, vg.ensembl_gene, vg.csq_mask, g.symbol FROM variant_genes vg "
                       "LEFT JOIN genes g ON g.ensembl_id = vg.ensembl_gene "
                       f"WHERE vg.chr = ? AND vg.pos IN ({','.join('?' * len(chunk))})")
                for p, gid, mask, sym in con.execute(sql, [chrom, *chunk]):
                    by_pos.setdefault((chrom, p), []).append((gid, mask))
                    symbols[gid] = sym
    else:
        bits = {t: b for b, (t, _, _) in bit_term.items()}
        for v in variants:
            if v["gene_id"] and v["max_csq"] in bits:
                by_pos[(v["chr"], v["pos"])] = [(v["gene_id"], 1 << bits[v["max_csq"]])]
        ids = sorted({v["gene_id"] for v in variants if v["gene_id"]})
        for i in range(0, len(ids), 5000):
            chunk = ids[i:i + 5000]
            symbols.update(con.execute(
                f"SELECT ensembl_id, symbol FROM genes WHERE ensembl_id IN ({','.join('?' * len(chunk))})",
                chunk).fetchall())

    for v in variants:
        out = []
        for gid, mask in by_pos.get((v["chr"], v["pos"]), []):
            names, top, score, impact = decode_mask(mask)
            out.append({"gene_id": gid, "symbol": symbols.get(gid), "consequence": top,
                        "impact": impact, "consequences": names,
                        "near": mask & ~flank == 0, "_score": score})
        out.sort(key=lambda g: (-g["_score"], g["gene_id"] != v["gene_id"],
                                (g["symbol"] or "").lower(), g["gene_id"]))
        for g in out:
            del g["_score"]
        v["genes"] = out
        v["gene_symbol"] = symbols.get(v["gene_id"]) if v["gene_id"] else None


def decode(row: dict) -> dict:
    """Enrich a raw variant row: decode sdpn -> founders / sdp / alleles / confidence."""
    sdpn, ref, conf = row["sdpn"], row["ref"], row["confidence"]
    bases = enc.sdpn_to_bases(sdpn)
    fis = enc.conf_bits(conf)
    sdp = enc.sdp_from_sdpn(sdpn, ref)
    return {
        "variants_id": row["snp_id"],
        "chr": row["chr"], "pos": row["pos"], "ref": ref,
        "alleles": enc.alleles_string(sdpn, ref),
        "sdpn": sdpn, "sdpn_display": enc.sdpn_to_display(sdpn), "sdpn_seq": enc.sdpn_to_sequence(sdpn),
        "sdp": sdp, "sdp_partition": enc.sdp_to_partition(sdp),
        "confidence": conf, "high_confidence": conf == 255,
        "founders": dict(zip(enc.FOUNDERS, bases)),
        "low_conf_founders": [f for f, ok in zip(enc.FOUNDERS, fis) if not ok],
        "max_csq": row["consequence"], "gene_id": row["ensembl_gene"],
        "_csq_mask": row.get("csq_mask"),   # decoded to `consequences` by attach_genes
    }


def search_variants(
    con, *, region: str | None = None, gene_symbol: str | None = None,
    rsid: str | None = None, consequence: str | None = None,
    low_confidence: bool = False, differ: str | None = None,
    private_to: str | None = None, sdpn: int | None = None, limit: int = 1000,
) -> dict[str, Any]:
    """Resolve a variant query. Raises ValueError / LookupError on bad input.

    `gene_symbol` takes a symbol or an Ensembl gene ID; an ambiguous symbol raises
    AmbiguousGeneError (a LookupError) listing the candidate loci.
    """
    where: list[str] = []
    params: list[Any] = []

    if gene_symbol:
        g = gene(con, gene_symbol)
        if not g:
            raise LookupError(f"gene {gene_symbol!r} not found")
        where.append("chr=? AND pos BETWEEN ? AND ?")
        params += [g["chr"], g["start"], g["end"]]
    elif region:
        r = enc.parse_region(region)
        where.append("chr=? AND pos BETWEEN ? AND ?")
        params += [r.chrom, r.start, r.end]
    elif rsid:
        where.append("snp_id = ?")
        params.append(rsid)
    else:
        raise ValueError("supply region, gene, or rsid")

    if not low_confidence:
        where.append("confidence = 255")
    if sdpn is not None:
        where.append("sdpn = ?"); params.append(sdpn)
    if differ:
        a, b = _founders(differ, 2)
        where.append(enc.differ_predicate(a, b))
    if private_to:
        (f,) = _founders(private_to, 1)
        where.append(enc.private_predicate(f))
    if consequence:
        bits = _term_bits(con)
        mask = 0
        for t in (x.strip() for x in consequence.split(",") if x.strip()):
            if t not in bits:
                raise ValueError(f"unknown consequence {t!r}")
            mask |= (1 << bits[t])
        where.append("(csq_mask & ?) <> 0"); params.append(mask)

    sql = (f"SELECT {','.join(_VAR_COLS)} FROM variants "
           f"WHERE {' AND '.join(where)} ORDER BY pos LIMIT ?")
    params.append(limit)
    rows = [dict(zip(_VAR_COLS, r)) for r in con.execute(sql, params).fetchall()]
    variants = [decode(r) for r in rows]
    attach_genes(con, variants)
    out: dict[str, Any] = {"count": len(rows), "variants": variants}
    if gene_symbol:
        out["gene"] = g   # the resolved locus: lets a caller show each SNP's effect IN this gene
    return out


def _founders(csv_names: str, n: int) -> list[str]:
    names = [x.strip() for x in csv_names.split(",") if x.strip()]
    if len(names) != n:
        raise ValueError(f"expected {n} founder name(s), got {names}")
    for name in names:
        if name not in enc.FOUNDERS:
            raise ValueError(f"unknown founder {name!r}; valid: {enc.FOUNDERS}")
    return names
