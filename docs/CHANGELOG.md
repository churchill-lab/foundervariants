# Changelog

Notable changes, newest first. The topical documents carry the *why* in depth —
[ARCHITECTURE.md](ARCHITECTURE.md) for design decisions, [DATABASES.md](DATABASES.md) for
data provenance and validation, [DEPLOYMENT.md](DEPLOYMENT.md) for operations. This file is
the chronology, so those documents do not have to read like a diary.

## 2026-11-04 (0.2.0): every gene per SNP (schema compact-v2)

### Added: every gene at a SNP, each with its own consequence

`variants` stores one gene per SNP, the gene of the single most severe consequence. Genes
overlap, so that hid the rest: on GRCm39, 2.8% of SNPs are in more than one gene body and
10.3% have more than one gene counting VEP's 5 kb flanks. At rs37654964 the UI showed
"synonymous · ENSMUSG00000110949" (Nudt8). The SNP is also synonymous in Gm49405 (Nudt8
only won a tie, on CSQ order) and in an exon of the lncRNA Gm16312. Gene searches misled
too: searching Gm49405 listed consequences that were really about Nudt8.

- **Build:** the CSQ loop also ORs each term into its gene's own mask, written to a new
  `variant_genes(chr, pos, ensembl_gene, csq_mask)` table (24.8 M rows on GRCm39, 21.0 M
  on GRCm38). Plus a unique index on `genes(ensembl_id)` for symbols. Schema `compact-v2`.
- **R/qtl2: unchanged.** `variants`, its index and every row are identical to production
  in both releases, and `create_variant_query_func()` returns identical data frames. See
  [DATABASES.md](DATABASES.md#the-compact-v2-rebuild-every-gene-per-snp-2026-11-04).
- **API:** each variant gains `genes` (every gene, its own most severe consequence,
  impact, all its terms, and `near` for flank-only) and `gene_symbol`. `max_csq` and
  `gene_id` are unchanged, and `genes[0]` always agrees with them. A `gene=` search also
  returns the locus searched as `data.gene`. CSV/TSV exports gain `gene_symbol` and
  `genes` columns. Files built before compact-v2 are still served, with the one stored
  gene.
- **API: every consequence term too.** Each variant gains `consequences`, all the SO
  terms in its `csq_mask`, most severe first (exports: a `consequences` column). The UI
  shows the extra terms as "+N", e.g. "synonymous_variant +5", as MUSTER does.
- **API: faster.** `/variants` now serializes its JSON with the stdlib encoder instead
  of FastAPI's, which walked every nested object in Python and cost more than the query.
  A 2,000-row search is ~45% faster than before over HTTP (≈78 → 43 ms), despite
  carrying the gene lists. The rest of the UI's ~0.7 s for 2,000 rows is the browser
  drawing the table, unchanged.
- **UI:** the gene column lists symbols, not IDs: two, then "+N", with nearby-only genes
  faded and every gene's consequences on hover. A gene search shows the consequence *in
  that gene*, adding "· more severe in X" when the SNP's worst effect is elsewhere.
- **Size:** 4.27 → 5.28 GB (GRCm39), 4.27 → 5.13 GB (GRCm38).
- Tests: 236 → 269.

### Added: a hover card in the results table

The UI's details were native `title` tooltips: slow to appear, unstyled, one per cell.
They're replaced by one floating card, after MUSTER's genotype-grid tooltip, built from
the row's data when the pointer enters a cell:

- **a founder's base:** the strain, its call (reference/alternate), confidence, site,
  alleles, reference, variant ID, consequence and gene, with the founder's column
  highlighted;
- **consequence or genes:** every consequence at the SNP, then gene by gene with each
  gene's own terms, Ensembl ID and "nearby" flag (up to 8 genes; the rest are counted);
- **SDP/SDPN:** which founders carry which allele;
- **column headers:** what the column means.

Terms use their readable names ("Synonymous variant") with an impact dot. The table no
longer carries thousands of `title` attributes.

### Added: NDJSON export

`format=ndjson` (and an NDJSON button): JSON Lines, one complete variant per line, with
its `genes` and `consequences` nested, which CSV can only flatten. Reads line by line with
`jq`, `pandas.read_json(lines=True)` and R's `jsonlite::stream_in()` (each checked on a
real export). JSON-LD was considered and kept out of the API; see ARCHITECTURE.md's
decisions log.

### Changed: ambiguous gene symbols are refused, not guessed

`gene=U6` used to resolve to whichever of U6's 486 loci SQLite returned first, with no
sign that there were others. About 110 GRCm39 symbols name several loci, nearly all
multi-copy small RNAs. Now:

- `gene=` also takes an Ensembl gene ID (`ENSMUSG…`).
- An ambiguous symbol gets **409** with the candidate loci (first 50, in genome order). An
  exact-case match wins over case-insensitive ones (`Mir1a` vs `MIR1A`).
- The UI offers the loci as buttons; picking one searches its Ensembl ID, so the link is
  unambiguous. The autocomplete does the same when a symbol appears twice in its list.

### Tooling

- `scripts/check_variant_genes.py`: checks `variant_genes` against `variants` on a built
  file (stdlib only, streaming, about 3 min per release).
- `scripts/check_rqtl2_same.R`: R/qtl2's view of two files must be identical.
- `scripts/compare_builds.sh` now expects the compact-v2 additions and still requires
  every production object to be byte-identical.
- The new tests include `parse_chrom` on a stub VCF record, so the
  per-gene masks are tested without the 24 GB VCF.

## Unreleased (docs only)

### DATABASES.md: why VCF `FILTER` is ignored; multiallelic sites across releases

From an email review whose author checked both databases against the source VCFs and
found no bugs but two undocumented behaviours. Both re-verified here before
documenting. The `FILTER` breakdown was reproduced exactly by re-running the build's
keep rule on chr19 of both VCFs. The multiallelic counts (3 bases: 177,225 / 180,678;
4 bases: 1 / 321) and the rs52115084 and rs233533684 examples match exactly. The
per-founder pattern at REL1505's four-base sites matches in shape (CAST 100, PWK 26,
WSB 20 here, against 95/26/19; likely a definitional difference). One suggestion was
not taken as stated: REL2021's low-confidence calls are *not* generally "the
reference" (37% are), so the docs describe the pattern as specific to those sites.
No code or data changed.

## 2026-10-21 (0.1.6)

Published `churchilllab/foundervariants:0.1.6` (linux/amd64, 55.3 MB, digest
`sha256:246b2695…`). Checked on the image under `/foundervariants`: every page and asset
200; the pop-up (open, Use in query, row clicks, Esc/backdrop, Enter, phone layout); the
row order; `/sdp` (editing, links, explanation); and the earlier checks (linking, Copy
link, GRCm39 default, early typing, keys), with no browser errors.

### SDP/SDPN explained, and a calculator that fills the SDPN filter

**On the search page:** an **SDP calc** button next to the SDPN filter opens a pop-up
calculator. Set each founder's base (A/C/G/T buttons per founder), or type into any of
SDPN, the 8 bases, the SDPN groups (`A:BC:DEFGH`), SDP, or the SDP groups
(`reference:different`); every other field updates. **Use in query** (or Enter) puts the
SDPN in the filter and reruns the search. Clicking a row's SDP or SDPN number opens the
pop-up on that SNP's pattern, so "every SNP in this region with this pattern" is two
clicks; a modifier-click still opens `/sdp` in a new tab. The pop-up opens on the current
filter value, closes on Esc or a click outside, and keeps its buttons visible on phones.

**`/sdp`** explains both encodings with real rows (rs29358506: SDP 1, SDPN 5461;
rs237124876, where SDP `ABCDEF:GH` hides that PWK and WSB carry different bases and
SDPN `ABCDEF:G:H` does not), with the same calculator. The column headers and a
top-bar "SDP/SDPN" link go there. `/sdp?sdpn=` (or `?sdp=`, `?bases=`) opens a pattern.

The calculator's rows run Bases A–H, SDP, SDP groups, SDPN, SDPN groups, Alleles: the
bases everything is calculated from, then each encoding with its groups.

Inputs that don't determine the bases are resolved predictably. SDP and the groupings
say *which* founders share an allele, not *which base*, so B's group keeps the
reference base and other founders keep a base they already have where it fits. An SDP
with B's bit set is complemented (B6 is always the reference here; 255 - x is the same
split), with a note saying so.

Built as files: `static/sdpcalc.js` (pure conversions), `static/sdpcalc-ui.js` (the
component, used by both pages), `static/sdpcalc.css`, `static/SDP.md` (the
explanation), `web/sdp.html`; `static/site.css` now holds the styles shared by
`/linking` and `/sdp`. `tests/test_sdpcalc.py` runs the shipped JS under node against
`encode.py` on 3,010 values in every field and on 4,000 real rows per release, which
recompute each row's stored SDP and alleles exactly. It also checks that editing any
field round-trips, and that bad input is an error, never a change. Browser tests (Chrome)
caught one bug before release: the `/sdp` page's own address-bar update erased its
`?sdpn=` link before it was read, fixed by reading the link first.

## 2026-10-17 (0.1.5)

Published `churchilllab/foundervariants:0.1.5` (linux/amd64, 55.3 MB, digest
`sha256:181061b3…`). Checked on the image under `/foundervariants`: the linking guide,
its assets and Copy link (checks below), plus the earlier checks (GRCm39 default, keys,
early typing, deep links, `/docs`).

### A "Link to us" page, and a Copy link button

**`/linking`** now serves the linking guide on the site itself, for whoever builds links
into it (e.g. the MPD team). The main page links to it from the top bar ("Link to us")
and the footer. The text lives in one place: `docs/LINKING.md` moved to
`src/foundervariants/web/static/LINKING.md` (`docs/` is excluded from the image), and
`web/linking.html` renders it in the browser with a vendored copy of marked 12.0.2 (MIT,
35 KB, no CDN at runtime, same approach as the Scalar bundle). `docs/LINKING.md` is now a
pointer. Relative links in the guide resolve against the site's base path, so they work
behind the proxy. The guide gained a "Try it" list of clickable example links, and
its old links to repo files (`API.md`, `ARCHITECTURE.md`, `CHANGELOG.md`) were
replaced, since they would 404 on the web.

**Copy link** in the results header copies the current view's URL, which the address bar
already holds. It uses the Clipboard API (HTTPS or localhost), falls back to a hidden
textarea, and as a last resort shows the URL to copy by hand.

Tests: the page and its assets are served, with and without a base path; the main page
links to it; the guide has no repo-relative links; and every "Try it" link is run
through the page's own `parseLink`, so a broken example fails the suite. In Chrome, at
the root and under `/foundervariants`: the guide renders (7 sections, 4 tables), in-page
anchors and Try-it links work (896 rows, peak marked), Copy link puts exactly the
address-bar URL on the clipboard, the theme carries across both pages, there is no
sideways scroll at 390 px, and there are no console errors or failed requests.

## 2026-10-17 (0.1.4)

Published `churchilllab/foundervariants:0.1.4` (linux/amd64, 55.3 MB, digest
`sha256:41434ce7…`). Checked on the image with `FV_BASE_PATH=/foundervariants`: default
release GRCm39, explicit `release=GRCm38` honoured, and the 0.1.3 browser checks (keys,
early typing, deep links, `/docs`) all pass.

### The UI defaults to GRCm39

The page opened on whichever release the server listed first, and the server lists
databases in filename order, so `fv.1505` (GRCm38) won. It now opens on the **newest
assembly**, currently GRCm39 (release 2021), which is also listed first in the release
menu. The rule is "newest assembly, then highest release id" (`byNewest`), not a
hardcoded `2021`, so a future GRCm40 database becomes the default when it is added.
Explicit links are unaffected: `release=GRCm38` still opens GRCm38. **Links that omit
`release` now open GRCm39 instead of GRCm38**, which is why LINKING.md now says to always
pass it; the MPD DO GWAS files are GRCm38. The address bar always records the release,
so URLs copied from the page are unaffected. Covered by 5 new tests in
`tests/test_deeplink.py`, and checked in Chrome.

## 2026-10-17 (0.1.3)

### Keyboard selection in the gene suggestions

The suggestion list under the search box only took mouse clicks. Now **↓/↑** move a
highlight (↓ wraps from the last row to the first; ↑ from the first returns to the typed
text), **Enter** searches the highlighted gene or, with nothing highlighted, the typed
text, and **Esc** closes the list. Hover moves the same highlight, so mouse and keyboard
don't fight. The box and list carry ARIA combobox/listbox roles and
`aria-activedescendant`, so screen readers announce the highlighted gene. The stylesheet
already had a `.acitem.sel` style that nothing applied.

Tested with scripted key presses in real Chrome (Playwright, using the installed Chrome),
which caught two timing bugs in the first version, both fixed:
- **Stale highlight.** Typing more while an older list showed a highlighted row, then
  pressing Enter quickly, searched that stale row (typed `Pten`, got `Cacna1c`). Typing
  now clears the highlight immediately.
- **List reopening.** The debounced lookup queued by the last keystroke landed after Enter
  or Esc and reopened the list. The debounce is now cancellable, and Enter and Esc cancel
  it; Esc does so even before the list has appeared.

Mouse clicks on suggestions and the example chips still work.

### Fixed: typing during page load was erased

Found by the same browser tests, run against the 0.1.3 image (slower under amd64
emulation): typing `Cac` right after the page opened left `ac` in the box. At startup the
deep-link code (0.1.1) wrote the link's search text and filters into the controls even
when no link was given, so it reset them to empty after the startup fetches and erased
anything the user had already typed or changed. A real user on a slow connection would
lose keystrokes the same way. It now sets only what the link specifies. Verified in
Chrome against the image: typing `Pten` immediately on load keeps it and searches it
(554 rows), and deep links still apply.

### Released 0.1.3

Ships the keyboard selection and early-typing fix above and the API-link fix below.
Published `churchilllab/foundervariants:0.1.3` (linux/amd64, 55.3 MB, digest
`sha256:bb16999d…`). Checked on the image with `FV_BASE_PATH=/foundervariants`: health
0.1.3, `scalar.js` 200, `href="/foundervariants/docs"`, and the browser tests above.

### Fixed: the UI's "API" link 404'd from the slash-less address

The link was `href="docs"`, which is relative. From `…/foundervariants/` it resolved
correctly, but from `…/foundervariants`, the form people type and where the
`/foundersnps` redirect lands, it resolved to `/docs`, which 404s. It is now
`href="%%BASE%%/docs"`, substituted server-side like the page's other base-path uses,
so it is absolute either way. Tests cover it with and without a base path.

## 2026-10-17

### Deployment docs now match the live server; `/foundervariants` is the main address

The docs and config described a planned deployment that differed from the real one. Read
from the server's `/etc/httpd` and the running container, then corrected throughout
(`apache-foundervariants.conf`, both compose files, `.env.example`, `server-bundle/`,
README, DEPLOYMENT.md, LINKING.md):

| | was documented | was live | now |
|---|---|---|---|
| base path | `/foundervariants` | `/foundersnps` | **`/foundervariants`**, with `/foundersnps` redirecting to it |
| host port | 28223 | 28224 | **28224** (ensimpl3 holds 28223) |
| server directory | `/srv/foundervariants` | a site-specific path | site-specific; set by `HOST_DB_DIR` in `.env` |
| host | `PUBLIC_HOST` placeholder | churchilllab.jax.org | same (Apache: `conf.d/vhost/churchilllab.conf`) |

**Making `/foundervariants` the main address** took three changes together: the Apache
`<Location>`, its `X-Forwarded-Prefix`, and the container's `FV_BASE_PATH`. Apache was
switched first, while the container still had `FV_BASE_PATH=/foundersnps`. That half-state
was misleading: `/foundervariants/` returned 200, but the page fetched its data from
`/foundersnps/...`, which 404'd, so in a browser it showed no release and no results.
`FV_BASE_PATH` was then changed in the server `.env` (backup `.env.bak-20261017`) and the
container recreated. Environment variables only apply to a new container; a restart keeps
the old value.

**The redirect** from the old address is
`RedirectMatch 301 "^/foundersnps(/.*)?$" "/foundervariants$1"`. The first version written
on the server had `"^foundersnps..."`, without the leading slash, which can never match
because request paths start with `/`. It fails silently: no config error, just a 404. The
corrected line was tested in `httpd:2.4` in front of the app with the same `<Location>`
block, and needs a root edit and an Apache reload on the server.

The unproxied legacy container (`churchilllab/foundervariants:0.0.2-amd64`, port 15001,
up 19 months) was stopped. It had been started with auto-remove, so stopping it also
deleted it; the image is still on the server.

### Fixed: the API docs page was broken behind the proxy

`/docs` loads `static/scalar.js`, which 404'd on the live site. The static files were
served by a Starlette `Mount`, and with `root_path` set (`FV_BASE_PATH`) a Mount matches
only the prefixed path (`/foundervariants/static/...`). Apache strips the prefix before
forwarding, so the request never matched. Ordinary routes match either form, which is why
only static files failed and why it went unnoticed: no test ran with a base path. Static
files are now served by a plain route, which refuses anything that isn't a file directly
in `web/static/`. New tests in `tests/test_api.py` (`TestBehindAPathStrippingProxy`) run
the app the way the server does. They cover static files at the stripped path, the prefix
in generated pages, API routes, and traversal attempts.

### Fixed: the old-name redirect pointed at `http://`

With 0.1.2 live, `/foundersnps/...` redirected to `http://churchilllab.jax.org/foundervariants/...`,
and port 80 is closed, so old links failed after the redirect. TLS is terminated before
Apache, so Apache sees plain HTTP and builds absolute redirect URLs with that scheme.
The fix is an absolute target:
`RedirectMatch 301 "^/foundersnps(/.*)?$" "https://churchilllab.jax.org/foundervariants$1"`.
The local `httpd:2.4` test could not catch this because it served plain HTTP end to end.
A 301 is cached, so a browser that already followed the bad redirect may need its cache
cleared.

### Released 0.1.2

Published `churchilllab/foundervariants:0.1.2` (linux/amd64, 55.3 MB, digest
`sha256:617c86bb…`) with the static-file fix above. Before the push, the image ran with
`FV_BASE_PATH=/foundervariants` and prefix-stripped requests, as on the server:
`/static/scalar.js` returned 200 (`text/javascript`), an encoded traversal returned 404,
`/docs` and `/` carried the prefix, the peak region returned 896 variants, and `cyvcf2`
and `uv` were absent. Version bumped in the same places as 0.1.1. 0.1.0 and 0.1.1 remain
in the registry.

### Released 0.1.1

Published `churchilllab/foundervariants:0.1.1` (linux/amd64, 55.3 MB, digest
`sha256:2654bf9c…`), carrying the UI deep links below. Version bumped in `pyproject.toml`,
`__init__.py`, `uv.lock`, both compose files, `.env.example`, `server-bundle/` and the
deploy docs. `0.1.0` is untouched in the registry, so rollback is a one-line `.env` change.
Before the push, the image itself ran against the real `dbs/`: `/health` reported 0.1.1
with both releases, the deep-link page rendered the peak example in headless Chrome (896
rows, rs29358506 highlighted), `/static/scalar.js` served, and `cyvcf2` and `uv` were absent.
No database changes.

## 2026-10-16

### The web UI can be linked to

Requested so that Mouse Phenome Database GWAS pages for Diversity Outbred studies (e.g.
CSNA03) can link from a peak to the founders' genotypes there, as a DO-appropriate
alternative to MUSter. The page now reads its URL on load, fills the search box and
filters, and runs the search:

```
/?release=GRCm38&chr=10&pos=128350502            # ±50 kb around a peak, peak row highlighted
```

Parameters: `release` (release id, assembly name, or UCSC name), one selector (`q`,
`region`, `gene`, `rsid`, `chr`+`pos`+`window`, or `chr`+`start`+`end`), `mark`
(highlight positions or rsIDs), the existing filters (`consequence`, `private_to`, `sdpn`,
`low_confidence`), and `limit` (2,000 default, 25,000 max; was fixed at 2,000). Full
reference and the reasoning: [LINKING.md](LINKING.md). Every search also writes itself
back to the address bar, so any view can be shared.

**The default window is 50 kb, not the "±1 Mb" first suggested.** Measured at the
motivating peak (chr10:128.35 Mb, GRCm38), ±1 Mb holds 21,111 high-confidence founder
SNPs. With rows in position order and a 2,000-row cap, that would have shown only the
window's left edge, and the peak would never appear. ±50 kb holds 896. When a marked
position is cut off by the cap anyway, the page says so instead of silently showing the
wrong stretch.

Bad parameters are ignored and named in a note under the search box. They never block the
rest of the link. **No server change:** `/` already ignored its query string, and a new
test pins that a link returns the same page. The parser is pure JS between `// <deeplink>`
markers in `web/index.html`. `tests/test_deeplink.py` (50 tests) runs that exact block
under node, so the shipped code is what's tested (skipped where node is absent).
Checked end to end in headless Chrome against the real databases: the default window
renders 896 rows with rs29358506 highlighted; ±1 Mb at the default limit shows the
"past the first 2000 rows" note; ±1 Mb with `limit=25000` renders all 21,111 rows with
the peak highlighted; bad founder and consequence names are reported.

## 2026-09-09

### Corrected databases installed; container image built

The rebuilt databases verified on 2026-09-01 were swapped into `dbs/`. The superseded
2026-08-03 files were moved to **`dbs_backup_20260803/`** — deliberately *outside* `dbs/`,
both so the serving glob cannot see them and so `rsync dbs/ SERVER:…` ships 8.5 GB rather
than 17 GB. Verified after the swap: `PRAGMA quick_check` ok on both, 1505 now reporting
171 `start_lost` and 0 null consequences, `so_terms` at 39, and R/qtl2 drop-in re-confirmed
at four windows.

First build of the container image, **linux/amd64** (cross-built on an arm64 Mac for the
amd64 server). Confirmed inside the image: `fastapi` and `uvicorn` present, `cyvcf2`
absent, no databases baked in. Ran it under `compose.yaml` against the real `dbs/` mount —
healthy in 8 s, both releases discovered, `/`, `/docs` and `/static/scalar.js` all serving.

**Dockerfile: `uv` is now mounted per-`RUN` instead of `COPY`ed in.** The layer breakdown
showed `COPY /uv /uvx /bin/` accounting for 55 MB — a build-time-only tool sitting in the
runtime image. `RUN --mount=from=ghcr.io/astral-sh/uv:latest,…` keeps it out of every layer
and leaves no build tooling in the finished image. **80.5 MB → 58.2 MB.**

Published as `churchilllab/foundervariants:0.1.0` (digest `sha256:95e67f8f…`), verified
present in the registry as `linux/amd64` via `buildx imagetools inspect`.

Added **`server-bundle/`**: `compose.prod.yaml`, a pre-filled `.env` with the two
host-specific lines marked, the Apache block, and a `README-DEPLOY.txt`. The server needs
no source checkout — that bundle plus the rsynced databases is everything.

### Fixed: `apache-foundervariants.conf` had an invalid directive

The `<Location>` block carried a standalone `ProxyTimeout 1800`, which Apache rejects:

```
AH00526: Syntax error ... ProxyTimeout not allowed in <Location> context
```

`ProxyTimeout` is valid in server-config and `<VirtualHost>` context only. Replaced with
the per-route form, `ProxyPass http://localhost:28223 timeout=1800`, which is legal inside
`<Location>` and — better on a shared vhost — scopes the timeout to this route instead of
every proxy in the file. The config had never been exercised against a real Apache when it
was written. Every other directive in the block was audited against Apache's context rules
and is valid where it sits.

## 2026-09-01

Four changes, listed in dependency order. Together they remove a build stage, add a test
suite, and fix two correctness bugs — one of which the new tests found.

### Removed the DuckDB build stage

The pipeline was `VCF → DuckDB intermediate → pack → SQLite`. DuckDB had been the *served*
store in the pre-2026-08-03 design; after the pivot to a single SQLite per release its only
remaining job was to reconstruct each variant's `csq_mask` by aggregating 110M
transcript-level consequence rows (`bit_or(1 << bit) … GROUP BY variants_id`) — recomputing
something the VCF parse already had in hand and discarded.

Now `VCF → SQLite`, one streaming pass:

- `build.py` accumulates `csq_mask` while reading each record's CSQ entries and writes the
  served SQLite directly, in 100k-row batches, indexes last.
- `pack.py` **deleted**. `schema.py` rewritten as the DDL for the one artifact. The CLI
  surface is now just `build` and `serve`.
- The `build` extra dropped from `cyvcf2 + pyarrow + duckdb` to **`cyvcf2` alone**; the
  4 GB write-then-read of an intermediate is gone.

**Cost:** a correction to `so_terms.tsv` or the builder now needs a full re-parse (~48 min)
instead of a cheap re-convert (~3 min). Accepted deliberately — and as it turned out, the
re-convert could never have fixed the vocabulary bug below anyway.

**Verified:** chr19 rebuilt and diffed against the existing database — 881,574 rows, **0
differences** across all 12 columns, `csq_mask` identical, stored schema text
byte-identical, R/qtl2 drop-in confirmed.

### Added a test suite

`tests/` went from empty to **144 tests running in ~2 s**, with no network and no VCF.
Fixture databases are built from `schema.py` itself, so the production schema is what is
under test, and expectations are computed by hand rather than captured from a run. Per-file
coverage is in the [README](../README.md#development).

Also fixed here: `-v` is a **top-level** flag (`foundervariants -v build …`), which every
documentation example had wrong; a test now pins it.

### Fixed `parse_region`: a trailing unit applied only to the range end

```
before:  19:5.0-5.2Mb  ->  start=5        end=5200000
         19:100-200kb  ->  start=100      end=200000
after:   19:5.0-5.2Mb  ->  start=5000000  end=5200000
```

Found by the new tests. Every abbreviated region query — the exact syntax the README and
API reference advertise — was silently scanning from the start of the chromosome. It
returned a *superset*, never an error, which is why it went unnoticed. A unit named on
either side now applies to both; mixed forms (`19:500kb-2Mb`) are unaffected.

### One consequence vocabulary across all releases

Two related defects, both stemming from Ensembl having renamed a consequence term.

1. **172 variants in release 1505 had the wrong consequence.** A junk `so_terms.tsv` row
   (`score = -1`, `impact = NULL`, `accession = NULL` — still visible in the legacy
   `fv.1505.snps.db3`) made `initiator_codon_variant` unable to win the top-consequence
   comparison at all. 55 rows ended up with `consequence = NULL`; the rest got whatever
   lower-severity term won by default, including 23 reported as `missense_variant`. The old
   two-step pipeline could not repair this by re-converting, because the top consequence is
   chosen during the **parse**.
2. **Filtering could not span the two builds.** VEP 75 (GRCm38) emits
   `initiator_codon_variant`; VEP 104 (GRCm39) emits `start_lost`. So
   `consequence=start_lost` returned 167 rows on GRCm39 and a **silent zero** on GRCm38 —
   in a tool whose purpose is comparing founder genotypes across builds.

Fixed by canonicalising deprecated terms at a single point in the parse
(`build.DEPRECATED_TERMS`), so severity ranking, the top-consequence choice and the
`csq_mask` bit can never disagree about which term they mean. `initiator_codon_variant` was
removed from `so_terms.tsv` (39 terms now) — safe **only** because it was the last row, bit
39; `start_lost` stays at bit 6, so no existing mask changed meaning. The retired name now
returns `422 unknown consequence` rather than an empty result.

An interim patch had scored the term `21 / LOW`. That was invented rather than sourced, and
wrong: **Ensembl 75 — the VEP version that annotated REL-1505 — ranked it 7th most severe**,
right after `stop_lost` and right before `inframe_insertion`, exactly the slot modern
Ensembl gives `start_lost` (33 / HIGH). e75 had no IMPACT column at all; its ordering *was*
the severity, so no Ensembl release ever rated it LOW. At 21 it ranked "start codon
destroyed" *below* `start_retained_variant` (22).

That distinction was not cosmetic: at 21 the term still lost to `missense_variant` (29) and
`splice_region_variant` (27), so only 143 rows were repaired. At its correct rank it wins
those too — the 23 + 5 rows that take the total to 171 top-consequence corrections out of
172 affected (the last is correctly outranked by `stop_gained`). Full evidence, sources, and
the one case VEP 75 cannot distinguish:
[DATABASES.md](DATABASES.md#deprecated-vep-terms-are-canonicalised-on-read).

### Databases

Both releases were rebuilt and compared against production at full scale before anything
was swapped; production was never written to. Result: **2021 identical (0 rows differ),
1505 differs on 172 rows — all of them the start-codon fix above, and none of them
genotype data.** Method and full results:
[DATABASES.md § Validation](DATABASES.md#validation).

R/qtl2 was re-verified on both rebuilt files with `create_variant_query_func()` and no
field overrides, including all 603 rows of 1505's chrMT.

Tooling for this now lives in `scripts/`: `rebuild_verify.sh` builds into `dbs_verify/`,
`compare_builds.sh` diffs it against `dbs/` with an indexed merge join across all 12
columns plus a key-coverage check. (A two-directional whole-table `EXCEPT` is equivalent but
has to sort 36M twelve-column tuples — it ran over an hour before being abandoned.)

`build_cache/`, holding 6.8 GB of intermediates from the removed stage, was deleted.

### Deployment artifacts (earlier the same day)

The repo had no container or proxy configuration. Added `Dockerfile`, `.dockerignore`,
`compose.yaml`, `compose.prod.yaml`, `apache-foundervariants.conf`, and an annotated
`.env.example`, plus [DEPLOYMENT.md](DEPLOYMENT.md) and [DATABASES.md](DATABASES.md).
Conventions: host port `28223`, base path `/foundervariants`, image
`churchilllab/foundervariants:0.1.0`, server directory `/srv/foundervariants`.

Supporting fixes: `FV_LOG_LEVEL` is now actually honoured (under Docker/systemd uvicorn
imports the app directly, so the CLI never ran and startup discovery logged nothing), and
`pack`-time `meta` keys were merged rather than appended so `meta` holds one row per key.
