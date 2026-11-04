"""
FastAPI app for FounderVariants (SQLite-only).

Discovers `fv.*.snps.sqlite` in FV_DB_DIR (reads each `meta`), serves release-first
read-only queries on worker threads (a fresh read-only/immutable connection per
request — safe to fan out), exports result sets as csv/tsv/json, and offers the
same SQLite as the R/qtl2 download. Mirrors the ensimpl service shape.
"""
from __future__ import annotations
import csv
import glob
import io
import json
import os
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, HTMLResponse, Response
from starlette.concurrency import run_in_threadpool

from . import __version__
from . import encode as enc
from . import log as _log
from . import query as Q
from .log import get_logger

LOG = get_logger()
WEB = os.path.join(os.path.dirname(__file__), "web")
# ndjson (JSON Lines): one complete variant per line, nested genes/consequences intact,
# readable line by line (jq, pandas read_json(lines=True), jsonlite::stream_in).
_EXPORT_MEDIA = {"csv": "text/csv", "tsv": "text/tab-separated-values", "json": "application/json",
                 "ndjson": "application/x-ndjson"}


def _discover(db_dir: str) -> dict[str, dict]:
    reg: dict[str, dict] = {}
    for path in sorted(glob.glob(os.path.join(db_dir, "fv.*.snps.sqlite"))):
        try:
            con = Q.connect(path)
            meta = Q.meta(con)
            con.close()
            rel = meta.get("release") or os.path.basename(path).split(".")[1]
            reg[str(rel)] = {"path": path, "meta": meta}
            LOG.info("loaded release %s from %s", rel, os.path.basename(path))
        except Exception as e:  # noqa: BLE001 — skip half-written DBs, keep serving
            LOG.warning("skipping %s: %s", path, e)
    return reg


@asynccontextmanager
async def lifespan(app: FastAPI):
    db_dir = os.environ.get("FV_DB_DIR", "dbs")
    app.state.registry = _discover(db_dir)
    LOG.info("FounderVariants %s — %d release(s): %s",
             __version__, len(app.state.registry), ",".join(app.state.registry) or "none")
    yield


def create_app() -> FastAPI:
    # Configure logging here too: under Docker/systemd uvicorn imports this module
    # directly, so the CLI (which calls log.configure) never runs. basicConfig is a
    # no-op when handlers already exist, so `serve -v` still wins.
    _log.configure(os.environ.get("FV_LOG_LEVEL", "INFO"))
    base = os.environ.get("FV_BASE_PATH", "")
    # Swagger (docs_url) is replaced by a vendored Scalar page at /docs (see below).
    app = FastAPI(title="FounderVariants", version=__version__,
                  root_path=base, docs_url=None, lifespan=lifespan)
    origins = os.environ.get("FV_CORS_ORIGINS", "").split(",") if os.environ.get("FV_CORS_ORIGINS") else ["*"]
    app.add_middleware(GZipMiddleware, minimum_size=1000)
    app.add_middleware(CORSMiddleware, allow_origins=origins, allow_methods=["GET"], allow_headers=["*"])

    def db(request: Request, release: str) -> dict:
        entry = request.app.state.registry.get(str(release))
        if not entry:
            raise HTTPException(404, f"release {release!r} not found; have {list(request.app.state.registry)}")
        return entry

    async def run(entry, fn):
        def work():
            con = Q.connect(entry["path"])
            try:
                return fn(con)
            finally:
                con.close()
        return await run_in_threadpool(work)

    # ---------------------------------------------------------------- routes
    @app.get("/health")
    def health(request: Request):
        return {"status": "ok", "version": __version__, "releases": list(request.app.state.registry)}

    @app.get("/genomes")
    def genomes(request: Request):
        return {"data": [e["meta"] for e in request.app.state.registry.values()]}

    @app.get("/{release}/meta")
    async def meta(request: Request, release: str):
        e = db(request, release)
        n = await run(e, lambda c: c.execute("SELECT count(*) FROM variants").fetchone()[0])
        return {"data": {"meta": e["meta"], "n_variants": n, "founders": enc.FOUNDERS}}

    @app.get("/{release}/so_terms")
    async def so_terms(request: Request, release: str):
        return {"data": await run(db(request, release), Q.so_terms)}

    @app.get("/{release}/genes")
    async def genes(request: Request, release: str, q: str, limit: int = 25):
        return {"data": await run(db(request, release), lambda c: Q.genes_search(c, q, limit))}

    @app.get("/{release}/variants")
    async def variants(
        request: Request, release: str,
        region: str | None = None, gene: str | None = None, rsid: str | None = None,
        consequence: str | None = None, low_confidence: bool = False,
        differ: str | None = None, private_to: str | None = None,
        sdpn: int | None = None, limit: int = 1000, format: str = "",
    ):
        e = db(request, release)
        try:
            res = await run(e, lambda c: Q.search_variants(
                c, region=region, gene_symbol=gene, rsid=rsid, consequence=consequence,
                low_confidence=low_confidence, differ=differ, private_to=private_to,
                sdpn=sdpn, limit=limit))
        except Q.AmbiguousGeneError as ex:
            # 409, not 404: the gene exists, at more than one place. `candidates` lets a
            # client offer the loci; searching by one's Ensembl ID resolves it.
            raise HTTPException(409, {"message": str(ex), "symbol": ex.symbol,
                                      "total": ex.total, "candidates": ex.candidates}) from ex
        except LookupError as ex:
            raise HTTPException(404, str(ex))
        except ValueError as ex:
            raise HTTPException(422, str(ex))
        exported = _export(res["variants"], format, f"foundervariants_{release}")
        if exported:
            return exported
        # Serialized here, not by FastAPI: its encoder walks every nested object in Python
        # (thousands of variants, each with a genes list), which cost more than the query.
        # The payload is already plain JSON types, so the stdlib C encoder is exact.
        body = await run_in_threadpool(
            json.dumps, {"data": res, "meta": {"release": release, "count": res["count"]}},
            separators=(",", ":"))
        return Response(body, media_type="application/json")

    @app.get("/{release}/download/rqtl2")
    def download_rqtl2(request: Request, release: str):
        e = db(request, release)
        return FileResponse(e["path"], filename=os.path.basename(e["path"]),
                            media_type="application/x-sqlite3")

    def page(name: str) -> HTMLResponse:
        path = os.path.join(WEB, name)
        if not os.path.exists(path):
            return HTMLResponse("<h1>FounderVariants</h1><p>UI not built yet.</p>")
        with open(path) as f:
            return HTMLResponse(f.read().replace("%%BASE%%", base))

    @app.get("/", response_class=HTMLResponse)
    def ui():
        return page("index.html")

    # How to link into the UI: renders static/LINKING.md in the browser
    @app.get("/linking", response_class=HTMLResponse, include_in_schema=False)
    def linking():
        return page("linking.html")

    # SDP/SDPN explainer + calculator: static/sdpcalc.js does the maths, static/SDP.md
    # is the explanation
    @app.get("/sdp", response_class=HTMLResponse, include_in_schema=False)
    def sdp():
        return page("sdp.html")

    # A plain route, not app.mount(StaticFiles): with root_path set (FV_BASE_PATH), a
    # Starlette Mount only matches the PREFIXED path, but Apache strips the prefix before
    # forwarding, so the docs page's scalar.js 404'd. Routes match the stripped path.
    @app.get("/static/{name}", include_in_schema=False)
    def static(name: str):
        path = os.path.join(WEB, "static", name)
        if os.path.basename(name) != name or name.startswith(".") or not os.path.isfile(path):
            raise HTTPException(404)
        return FileResponse(path)

    @app.get("/docs", include_in_schema=False)
    def api_docs():
        return HTMLResponse(_scalar_html(base))

    return app


def _scalar_html(base: str) -> str:
    """Interactive API reference via vendored Scalar (no CDN), bound to /openapi.json.
    `base` prefixes the script + schema URLs so it works at the root or behind a proxy."""
    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<title>FounderVariants · API reference</title>"
        "<style>body{margin:0}</style></head><body>"
        f"<script id='api-reference' data-url='{base}/openapi.json'></script>"
        f"<script src='{base}/static/scalar.js'></script>"
        "</body></html>"
    )


def _gene_cell(g: dict) -> str:
    """One gene for a flat export: `Nudt8:ENSMUSG00000110949:synonymous_variant`, `~` = near only."""
    return f"{'~' if g['near'] else ''}{g['symbol'] or ''}:{g['gene_id']}:{g['consequence'] or ''}"


def _flatten(v: dict) -> dict:
    row = {k: v[k] for k in ("variants_id", "chr", "pos", "ref", "alleles", "sdp", "sdpn",
                             "sdp_partition", "sdpn_display", "confidence", "max_csq", "gene_id",
                             "gene_symbol")}
    row["consequences"] = "|".join(v.get("consequences", []))   # every term, most severe first
    # every gene at the SNP, most severe first ("|"-separated; see _gene_cell)
    row["genes"] = "|".join(_gene_cell(g) for g in v.get("genes", []))
    row.update({f"founder_{k}": val for k, val in v["founders"].items()})
    row["low_conf_founders"] = "|".join(v["low_conf_founders"])
    return row


def _export(rows: list[dict], fmt: str, stem: str):
    if fmt not in _EXPORT_MEDIA:
        return None
    if fmt == "json":
        body = json.dumps(rows).encode()
    elif fmt == "ndjson":
        body = "".join(json.dumps(r, separators=(",", ":")) + "\n" for r in rows).encode()
    else:
        flat = [_flatten(r) for r in rows]
        delim = "," if fmt == "csv" else "\t"
        buf = io.StringIO()
        if flat:
            w = csv.DictWriter(buf, fieldnames=list(flat[0].keys()), delimiter=delim)
            w.writeheader()
            w.writerows(flat)
        body = buf.getvalue().encode()
    return Response(body, media_type=_EXPORT_MEDIA[fmt],
                    headers={"Content-Disposition": f'attachment; filename="{stem}.{fmt}"'})


app = create_app()
