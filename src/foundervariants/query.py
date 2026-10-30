"""
Read side: pure functions over a read-only SQLite connection, returning plain
dicts/lists. Usable as a library with only the stdlib: `from foundervariants import query`.

The single compact SQLite is both the served store and the R/qtl2 file. Each
returned variant decodes per-founder bases, sdp, alleles, and per-founder
confidence ON READ from (sdpn, confidence). Consequence filtering uses the
per-variant `csq_mask` bitmask (no transcript-level table).
"""
from __future__ import annotations
import sqlite3
from typing import Any

from . import encode as enc

# columns read from the variants table (R/qtl2 names + our extras)
_VAR_COLS = ["snp_id", "chr", "pos", "ref", "sdpn", "confidence", "consequence", "ensembl_gene"]


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


def genes_search(con, q: str, limit: int = 25) -> list[dict]:
    cols = ["ensembl_id", "symbol", "chromosome", "start_position", "end_position", "strand"]
    rows = con.execute(
        f"SELECT {','.join(cols)} FROM genes WHERE symbol LIKE ? ORDER BY symbol LIMIT ?",
        [q + "%", limit]).fetchall()
    out = ["ensembl_id", "symbol", "chr", "start", "end", "strand"]
    return [dict(zip(out, r)) for r in rows]


def gene(con, symbol: str) -> dict | None:
    r = con.execute(
        "SELECT ensembl_id,symbol,chromosome,start_position,end_position,strand "
        "FROM genes WHERE symbol = ? COLLATE NOCASE LIMIT 1", [symbol]).fetchone()
    if not r:
        return None
    return dict(zip(["ensembl_id", "symbol", "chr", "start", "end", "strand"], r))


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
    }


def search_variants(
    con, *, region: str | None = None, gene_symbol: str | None = None,
    rsid: str | None = None, consequence: str | None = None,
    low_confidence: bool = False, differ: str | None = None,
    private_to: str | None = None, sdpn: int | None = None, limit: int = 1000,
) -> dict[str, Any]:
    """Resolve a variant query. Raises ValueError / LookupError on bad input."""
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
    return {"count": len(rows), "variants": [decode(r) for r in rows]}


def _founders(csv_names: str, n: int) -> list[str]:
    names = [x.strip() for x in csv_names.split(",") if x.strip()]
    if len(names) != n:
        raise ValueError(f"expected {n} founder name(s), got {names}")
    for name in names:
        if name not in enc.FOUNDERS:
            raise ValueError(f"unknown founder {name!r}; valid: {enc.FOUNDERS}")
    return names
