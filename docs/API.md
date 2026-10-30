# FounderVariants — API Reference

FastAPI service over the per-release SQLite databases. Paths are **release-first**
(`/{release}/…`, e.g. `/2021/…`). All responses are JSON unless a `format` is
requested. Base URL below assumes the local dev server; behind a proxy set
`FV_BASE_PATH` and prefix accordingly — see
[DEPLOYMENT.md](DEPLOYMENT.md#reverse-proxy).

To link people to the **web UI** with a search filled in (e.g. from a GWAS peak), use
the page parameters in the linking guide, served at `/linking` (source:
[`web/static/LINKING.md`](../src/foundervariants/web/static/LINKING.md)); this document covers the data API.

Companion documents: [ARCHITECTURE.md](ARCHITECTURE.md) for the encodings and storage
model, [DATABASES.md](DATABASES.md) for what each release contains and where it came
from, [DEPLOYMENT.md](DEPLOYMENT.md) for running and deploying the service.

```
uv run foundervariants serve --dir dbs        # http://127.0.0.1:8000
```

Interactive docs: `GET /docs` — a vendored **Scalar** API reference (self-contained,
no CDN), bound to the auto-generated `GET /openapi.json`. (Swagger is disabled.)

---

## Service endpoints

### `GET /health`
Liveness; touches no database.
```json
{ "status": "ok", "version": "0.1.6", "releases": ["1505", "2021"] }
```

### `GET /genomes`
All available releases with full provenance (`meta`).
```json
{ "data": [ { "release": "2021", "grc": "GRCm39", "created": "2026-08-03T00:01:46+00:00",
  "gene_annotation_file": "ensimpl.104.mm.genes.tsv", "n_variants": "35805856", ... } ] }
```

### `GET /{release}/meta`
```json
{ "data": { "meta": { ...provenance... }, "n_variants": 35805856,
            "founders": ["A_J","C57BL_6J","129S1_SvImJ","NOD_ShiLtJ","NZO_HlLtJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"] } }
```

### `GET /{release}/so_terms`
Consequence vocabulary (drives the UI filter). Each has a `bit` (its position in `csq_mask`).
```json
{ "data": [ { "term": "missense_variant", "display": "Missense variant",
              "impact": "MODERATE", "score": 22, "bit": 17 }, ... ] }
```

### `GET /{release}/genes?q=<prefix>&limit=25`
Self-contained gene-symbol autocomplete.
```json
{ "data": [ { "ensembl_id": "ENSMUSG00000013663", "symbol": "Pten",
              "chr": "19", "start": 32734897, "end": 32803560, "strand": 1 } ] }
```

---

## `GET /{release}/variants` — the primary endpoint

Exactly one of **`region`**, **`gene`**, or **`rsid`** selects the variants; the rest
filter.

| param | type | meaning |
|---|---|---|
| `region` | str | `19:5.0-5.2Mb`, `19:5000000-5200000`, `chr19:5.0-5.2M`… (Mb/kb suffixes ok) |
| `gene` | str | gene **symbol** (resolved via the local `genes` table → its span) |
| `rsid` | str | exact `snp_id` (rsID, or `chr:pos_ref_alt` positional id) |
| `consequence` | str | comma-separated SO terms; variant matches if it has **any** (via `csq_mask`) |
| `low_confidence` | bool | default `false` → only `confidence=255`; `true` includes low-confidence |
| `private_to` | str | founder name; variants where that founder differs from **all** others |
| `differ` | str | `FounderA,FounderB`; variants where the two founders carry different bases |
| `sdpn` | int | exact SDPN match |
| `limit` | int | max rows (default 1000) |
| `format` | str | `csv` \| `tsv` \| `json` → download the result set instead of the JSON envelope |

Founder names (for `private_to` / `differ`): `A_J, C57BL_6J, 129S1_SvImJ,
NOD_ShiLtJ, NZO_HlLtJ, CAST_EiJ, PWK_PhJ, WSB_EiJ`.

### Response
```json
{
  "data": {
    "count": 1,
    "variants": [
      {
        "variants_id": "rs30847235", "chr": "19", "pos": 32735002, "ref": "G",
        "alleles": "G|T",
        "sdpn": 60330, "sdpn_display": "BCEFGH:AD", "sdpn_seq": "TGGTGGGG",
        "sdp": 9, "sdp_partition": "BCEFGH:AD",
        "confidence": 255, "high_confidence": true,
        "founders": { "A_J":"T","C57BL_6J":"G","129S1_SvImJ":"G","NOD_ShiLtJ":"T",
                      "NZO_HlLtJ":"G","CAST_EiJ":"G","PWK_PhJ":"G","WSB_EiJ":"G" },
        "low_conf_founders": [],
        "max_csq": "5_prime_UTR_variant", "gene_id": "ENSMUSG00000013663"
      }
    ]
  },
  "meta": { "release": "2021", "count": 1 }
}
```

Errors: `404` (unknown release / gene / rsid not found), `422` (bad region, unknown
founder or consequence, or no selector given). Note founder-pattern queries scoped
by `region`/`gene` are fast; genome-wide (`region=19:0-200Mb&private_to=…`) scans the
chromosome (~1 s).

### Examples

```bash
# region
curl 'http://127.0.0.1:8000/2021/variants?region=19:5.0-5.2Mb&limit=100'
# gene + multi-term consequence filter, high-confidence only
curl 'http://127.0.0.1:8000/2021/variants?gene=Pten&consequence=missense_variant,stop_gained'
# founder pattern: private to WSB, download as CSV
curl -OJ 'http://127.0.0.1:8000/2021/variants?region=19:0-61Mb&private_to=WSB_EiJ&format=csv'
# rsID
curl 'http://127.0.0.1:8000/2021/variants?rsid=rs30847235'
```

```python
import httpx
r = httpx.get("http://127.0.0.1:8000/2021/variants",
              params={"gene": "Pten", "consequence": "missense_variant"})
for v in r.json()["data"]["variants"]:
    print(v["variants_id"], v["founders"], v["max_csq"])
```

---

## `GET /{release}/download/rqtl2`
Streams the release's SQLite file — which **is** the served database. Use it locally
with R/qtl2:

```r
library(qtl2)
qv <- create_variant_query_func("fv.2021.snps.sqlite")   # no field overrides needed
res <- qv("19", 5.0, 5.2)     # start/end Mbp; pos returned in Mbp
# columns: snp, chr, pos, alleles, sdp, sdpn, ensembl_gene, consequence, ...
idx <- qtl2::index_snps(map, res)   # groups by sdp for scan1snps
```

---

## CSV/TSV export columns

`format=csv|tsv` flattens each variant: `variants_id, chr, pos, ref, alleles, sdp,
sdpn, sdp_partition, sdpn_display, confidence, max_csq, gene_id`, plus one
`founder_<STRAIN>` column per founder and `low_conf_founders` (|-joined).
`format=json` returns the bare decoded-variant array.

---

## Environment

| var | default | meaning |
|---|---|---|
| `FV_DB_DIR` | `dbs` | directory scanned for `fv.*.snps.sqlite` |
| `FV_BASE_PATH` | `""` | reverse-proxy prefix (e.g. `/foundervariants`) |
| `FV_CORS_ORIGINS` | `*` | comma-separated allowed origins |
| `FV_LOG_LEVEL` | `INFO` | logging level |
