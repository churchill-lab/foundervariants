# Building, running, and deploying FounderVariants

End-to-end: set up locally, build databases from Sanger VCFs, verify, run the server, build
and publish the container image, copy databases to a server, deploy behind a reverse proxy,
and operate it.

Companion documents: [`DATABASES.md`](DATABASES.md) records how each existing database was
created and where its inputs came from, [`ARCHITECTURE.md`](ARCHITECTURE.md) explains the
encodings and storage decisions, and [`API.md`](API.md) is the endpoint reference.

## What actually runs

One process serving one directory of SQLite files.

```
             ┌─ Apache (public host, TLS) ───────────────────────┐
  browser ──▶│  <Location /foundervariants>  strips the prefix   │
             │  ProxyPass http://localhost:HOST_PORT             │
             └──────────────┬───────────────────────────────────-┘
                            ▼
             ┌─ Docker container ───────────────────────────────┐
             │  uvicorn foundervariants.api:app --workers 2     │
             │  reads FV_DB_DIR=/data/dbs                       │
             └──────────────┬───────────────────────────────────┘
                            ▼  bind mount, read-only
              /srv/foundervariants/dbs/fv.*.snps.sqlite  (on the host)
```

Three things follow from this shape and explain most operational surprises:

- **The image never contains databases.** They are mounted read-only at runtime, so one
  image serves any set of releases: the image is `python:3.12-slim` plus `fastapi`,
  `uvicorn` and the source, while each database is ~4.3 GB. It also means the registry
  never carries the data — databases are copied separately (step 6).
- **Database discovery happens once, at startup.** Availability is whatever
  `fv.*.snps.sqlite` files exist in the directory when the process starts. Adding a file to
  a running service does nothing until you restart it.
- **Databases are opened `immutable=1`.** Every request opens its own
  `file:…?mode=ro&immutable=1` connection, which is why requests fan out across workers
  with no locking — but it is also a promise to SQLite that the bytes will not change while
  open. Replace a database by writing a *new* file and restarting, never by rewriting the
  file a running server holds open (see step 6).

## Prerequisites

| Where | Needs |
|---|---|
| Development machine | Python 3.12+, [uv](https://docs.astral.sh/uv/), Docker (only to build images) |
| Development machine, to *build* a release | the `build` extra (`cyvcf2`), the source VCF (~24 GB), and ~4.5 GB free per release |
| Server | Docker + the Compose plugin, Apache with `mod_proxy`, `mod_proxy_http`, `mod_headers`, and ~4.3 GB disk per release |

No database engine to install and no application server to configure — SQLite is in the
standard library and uvicorn ships in the image. The server needs **none** of the build
dependencies: serving is pure stdlib `sqlite3` plus FastAPI.

## 1. Local setup

```bash
git clone <repo-url> foundervariants && cd foundervariants
uv sync                       # query + serve (fastapi + uvicorn; query core is stdlib sqlite3)
uv sync --extra build         # + cyvcf2, only needed to build from a VCF
```

The extra is separate on purpose: a machine that only serves databases never installs
`cyvcf2`, and the import graph is checked to keep it that way (see
[README](../README.md#development)).

## 2. Build a database

Full detail — the source VCFs, their provenance, what gets filtered out and why, and the
exact commands behind the two databases that exist — is in
[`DATABASES.md`](DATABASES.md). The short version, for **GRCm39 / release 2021**:

```bash
uv run foundervariants -v build \
  --vcf ~/work/RESOURCE/Sanger/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.rsID.vcf.gz \
  --genes-tsv ../data/ensimpl.104.mm.genes.tsv \
  --so-tsv    ../data/so_terms.tsv \
  --release 2021 --grc GRCm39 --gene-version "Ensembl 104 (ensimpl.104.Mm)" \
  --out dbs
```

Writes `dbs/fv.2021.snps.sqlite` — one file that is both the served store and the R/qtl2
variant database. Roughly **50 minutes** and ~4.3 GB for a whole genome, in a single
streaming pass with no intermediate file; the VCF parse is essentially all of it (index
creation is under a minute).

Note that `-v` is a top-level flag and goes **before** the subcommand.

> `--genes-tsv` and `--so-tsv` live in the **parent** project directory (`../data/`), not in
> this repo — they are shared inputs, and the gene TSVs are legacy ensimpl v1/v2 exports
> that current ensimpl no longer produces. See
> [`DATABASES.md`](DATABASES.md#gene-annotation-the-input-that-no-longer-regenerates).

Useful flags:

- `--chroms 19` — build one chromosome (or `1,2,X`). A single-chromosome build takes about a
  minute and is the right way to smoke-test a pipeline change before spending an hour.
- `--keepno` / `--no-keeplow` — change which founder calls are kept. The defaults reproduce
  the legacy database exactly; changing them makes the output no longer comparable to it.

There is no separate export, conversion or post-processing step: a correction to
`so_terms.tsv`, the gene TSV, or the builder itself is applied by re-running `build`. That
costs a full re-parse, which is the deliberate trade for a pipeline with one input, one
output and one dependency.

### Verify the build

Every database describes itself; there is no manifest to consult:

```bash
sqlite3 dbs/fv.2021.snps.sqlite "SELECT meta_name || ' = ' || meta_value FROM meta ORDER BY meta_name"
sqlite3 dbs/fv.2021.snps.sqlite "SELECT count(*) FROM variants"     # must equal meta n_variants
```

And confirm the same file still works as an R/qtl2 variant database, which is the half of
its job that the API cannot exercise:

```r
library(qtl2)
qv  <- create_variant_query_func("dbs/fv.2021.snps.sqlite")   # no field overrides
res <- qv("19", 5.0, 5.2)                                     # start/end in Mbp
stopifnot(nrow(res) > 0, all(res$sdp >= 1 & res$sdp <= 254))
idx <- index_snps(map, res)                                   # groups by sdp for scan1snps
```

## 3. Test

```bash
uv run pytest            # unit + API tests, no network, ~2 s
uv run ruff check .      # lint
```

The suite builds tiny fixture databases from `schema.py` and exercises the encodings, the
query layer, the whole HTTP API, database discovery, and the R/qtl2 column contract — see
the [README](../README.md#development) for the per-file breakdown. It does **not** cover the
VCF parse, which is I/O against a 24 GB file; that is exercised by actually building a
chromosome and diffing it (below).

`ruff` currently reports 26 pre-existing style findings and no bugs, so treat a non-zero
exit from it as noise rather than a gate.

A smoke test against real data is still worth running before any deploy, because it
exercises the two things that actually break in deployment (discovery finding the files,
and a query returning rows against a 4 GB database):

```bash
uv run foundervariants serve --dir dbs --port 8099 &
curl -s localhost:8099/health                                     # {"status":"ok",...,"releases":[...]}
curl -s 'localhost:8099/2021/variants?gene=Pten&limit=2' | head -c 400
curl -s 'localhost:8099/2021/genes?q=Pt&limit=3'
kill %1
```

Check that `releases` lists every release you expect — the list is a real integrity check,
because a file that cannot be opened or has no readable `meta` is **skipped with a warning
rather than failing startup**, so a broken database looks like a missing one.

The unit suite cannot catch this class of problem: it proves the code is right, while the
smoke test proves *this host's* databases are readable and discoverable.

## 4. Run the server locally

**Foreground** (Ctrl-C to stop):

```bash
uv run foundervariants -v serve --dir dbs --host 127.0.0.1 --port 8000
```

Then open <http://127.0.0.1:8000/> for the UI, `/docs` for the API reference.

**Auto-reload while editing code:**

```bash
uv run foundervariants serve --dir dbs --reload
```

**Background (daemon-ish) for a quick local session:**

```bash
nohup uv run foundervariants serve --dir dbs --port 8000 > fv.log 2>&1 &
echo $! > fv.pid               # remember the PID so you can stop it cleanly

tail -f fv.log                 # watch it
kill "$(cat fv.pid)"           # stop it
```

`pkill -f "foundervariants.api:app"` also works, but it will match every instance — prefer
the PID file if you run more than one.

For anything long-lived, use Docker (below) or a systemd unit rather than `nohup`; neither
survives a reboot on its own.

<details>
<summary>systemd unit, if you want it running without Docker</summary>

```ini
# /etc/systemd/system/foundervariants.service
[Unit]
Description=FounderVariants API
After=network.target

[Service]
Type=exec
User=foundervariants
WorkingDirectory=/srv/foundervariants/app
Environment=FV_DB_DIR=/srv/foundervariants/dbs
Environment=FV_BASE_PATH=/foundervariants
Environment=FV_LOG_LEVEL=INFO
ExecStart=/usr/local/bin/uv run uvicorn foundervariants.api:app --host 127.0.0.1 --port 28224 --workers 2
Restart=on-failure

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now foundervariants
sudo systemctl status foundervariants
sudo journalctl -u foundervariants -f
```

This route needs a source checkout and `uv sync` on the server, which is the reason Docker
is the recommended path.

</details>

## 5. Build the container image

The image holds only the application. Databases are mounted at runtime.

**Local build and run** — `compose.yaml` builds from source:

```bash
cp .env.example .env             # adjust HOST_PORT, HOST_DB_DIR, FV_BASE_PATH
docker compose up --build -d     # build + start, detached
docker compose logs -f           # follow the log
docker compose ps                # health status
docker compose restart           # pick up newly added .sqlite files
docker compose down              # stop and remove
```

It serves on `http://localhost:${HOST_PORT}` (default `28224`).

`.dockerignore` keeps `dbs/` out of the build context. Without it Docker would try to
upload ~8.5 GB of databases before starting the build — if a build appears to hang at
"transferring context", that file is the first thing to check.

**Publish a multi-architecture image** for a server that is not your machine's
architecture. The default builder cannot emit multi-arch manifests, so create a container
builder once:

```bash
docker buildx create --name fv-builder --driver docker-container --use
```

Then build and push in one step:

```bash
docker login                     # needs push rights to the target registry namespace
docker buildx build --platform linux/amd64,linux/arm64 \
  -t churchilllab/foundervariants:0.1.6 --push .
docker buildx imagetools inspect churchilllab/foundervariants:0.1.6   # confirm both platforms
```

Bump the tag per release rather than moving an existing one — a moved tag makes "what is
running in production" unanswerable. Keep the tag in step with `version` in
`pyproject.toml`.

## 6. Copy databases to the server

The registry carries the image, never the data. Databases go separately:

```bash
rsync -avP --partial dbs/ SERVER:/srv/foundervariants/dbs/
```

`--partial` matters — these are ~4.3 GB files and a dropped connection should resume rather
than restart.

**Do not pass `--inplace`.** Default rsync writes a temporary file and renames it over the
target, which leaves the running container holding the old inode: it keeps serving the
previous version, intact, until you restart it. `--inplace` rewrites the bytes underneath
that open `immutable=1` connection instead, and SQLite has been told those bytes will not
change — the result is wrong answers or a corruption error rather than a clean failure.

Copying is otherwise safe against a live service: the mount is read-only and discovery is
startup-only. That also means **you must restart to pick up new or replaced files**
(step 8).

## 7. Deploy

The server needs no source checkout — two files and the databases.

```bash
ssh SERVER
mkdir -p /srv/foundervariants && cd /srv/foundervariants
# copy compose.prod.yaml and a .env here (scp from your machine)
```

A minimal server `.env`:

```bash
IMAGE_NAME=churchilllab/foundervariants:0.1.6
HOST_PORT=28224                        # must be free; check first
HOST_DB_DIR=/srv/foundervariants/dbs   # absolute path to the copied databases
FV_DB_DIR=/data/dbs                    # mount point inside the container
FV_BASE_PATH=/foundervariants              # must equal the Apache <Location>
FV_LOG_LEVEL=INFO
```

Check the port is actually free before starting — a collision with an existing service is
the most common first-deploy failure. On churchilllab, ensimpl2 holds 28222 and ensimpl3
28223, which is why FounderVariants is on 28224.

```bash
ss -lntp | grep 28224 || echo "free"
```

Start it:

```bash
docker compose -f compose.prod.yaml pull
docker compose -f compose.prod.yaml up -d
docker compose -f compose.prod.yaml ps          # wait for health: healthy
curl -s localhost:28224/health                  # {"status":"ok","version":"…","releases":["1505","2021"]}
```

`compose.prod.yaml` pulls a published image and never builds, so it works on a host with no
source. It sets `pull_policy: always`, `restart: unless-stopped`, and a healthcheck that
hits `/health` using the interpreter already in the image.

### Reverse proxy

Copy the `<Location>` block from
[`../apache-foundervariants.conf`](../apache-foundervariants.conf) into the site's Apache
config (on churchilllab it lives in `/etc/httpd/conf.d/vhost/churchilllab.conf`), adjusting
the hostname, path and port if you deploy elsewhere. Then:

```bash
sudo apachectl configtest && sudo systemctl reload httpd
curl -s https://churchilllab.jax.org/foundervariants/health
curl -sI https://churchilllab.jax.org/foundervariants/static/scalar.js | head -1   # expect 200
```

Three things must agree, and nothing checks them for you:

| Apache | `.env` |
|---|---|
| `<Location /foundervariants>` | `FV_BASE_PATH=/foundervariants` |
| `ProxyPass http://localhost:28224` | `HOST_PORT=28224` |

**Trailing slashes are load-bearing.** Inside `<Location>`, the location path and the
`ProxyPass` target must both have a trailing slash or both lack one. Mixing them makes
Apache append the matched remainder to a target that already ends in `/`, so
`/foundervariants/genomes` reaches the app as `//genomes` and 404s — while
`/foundervariants` itself still returns the UI, so the page loads and only its API calls
fail.

### A second name: `/foundersnps`

The site was first served at `/foundersnps`, and links to that address exist. The main
address is now `/foundervariants`; keep the old one working by adding the `RedirectMatch`
line from the same config file at `<VirtualHost>` level, next to the `<Location>` block:

```apache
RedirectMatch 301 "^/foundersnps(/.*)?$" "https://churchilllab.jax.org/foundervariants$1"
```

Everything under `/foundersnps` gets a 301 to the same path under `/foundervariants`, and
the **query string carries over**, so deep links (the `/linking` guide), API calls
and CSV downloads all work under either name. **The pattern must start with `^/`.**
Request paths always begin with a slash, so `"^foundersnps(/.*)?$"` never matches: no
error, the old address just 404s. **The target must be an absolute `https://` URL.** TLS
is terminated before this Apache, so it sees plain HTTP and turns a relative target such
as `/foundervariants$1` into `http://churchilllab.jax.org/...`. Port 80 is closed, so the
browser follows the redirect and then cannot connect. It is a redirect, not a second `<Location>` proxy, on
purpose:
- **The app has one base path.** It writes `FV_BASE_PATH` into its pages, so a page served
  under a second proxied prefix would still call back to `/foundervariants`.
- **One canonical address** keeps bookmarks, shared links and access logs in one place.

**Renaming the main address takes three changes together:** the `<Location>` path and
`X-Forwarded-Prefix` in Apache, and `FV_BASE_PATH` in `.env`, which only takes effect when
the container is recreated (`docker compose -f compose.prod.yaml up -d`). If Apache and
`FV_BASE_PATH` disagree, the page itself still loads (200), but it fetches its data from
the other path and shows nothing.

Browsers cache a 301, so use 302 while the decision is still open. Verify after a reload:

```bash
curl -sI 'https://churchilllab.jax.org/foundersnps/?release=GRCm38&chr=10&pos=128350502' | grep -i -E '^(HTTP|location)'
# HTTP/1.1 301 ...   location: https://churchilllab.jax.org/foundervariants/?release=GRCm38&chr=10&pos=128350502
# The location MUST start with https:// (see above).
```

### The R/qtl2 download is a 4 GB proxied response

`GET /{release}/download/rqtl2` streams the entire database file back through Apache and
the app. It works, and `curl -OJ` gets a normal response with a `Content-Length`. But a
**browser** sends `Accept-Encoding: gzip`, and the app's `GZipMiddleware` then compresses
the stream on the fly, which drops `Content-Length` — so the browser shows no progress bar,
cannot resume, and burns CPU compressing a file that is already mostly incompressible
integers.

If this endpoint sees real use, serve the files off disk instead: the commented `Alias`
block at the bottom of `apache-foundervariants.conf` publishes `HOST_DB_DIR` directly, so
Apache handles range requests, resume and `Content-Length` natively. The API endpoint stays
useful for scripted single-shot fetches.

## 8. Operate

```bash
cd /srv/foundervariants
docker compose -f compose.prod.yaml ps          # status + health
docker compose -f compose.prod.yaml logs -f     # follow the log
docker compose -f compose.prod.yaml restart     # restart (also: pick up new .sqlite files)
docker compose -f compose.prod.yaml stop        # stop, keep the container
docker compose -f compose.prod.yaml start       # start it again
docker compose -f compose.prod.yaml down        # stop and remove
```

**Deploy a new version:**

```bash
# on the dev machine: build and push a new tag, then on the server
sed -i 's/foundervariants:0.1.5/foundervariants:0.1.6/' .env
docker compose -f compose.prod.yaml pull
docker compose -f compose.prod.yaml up -d       # recreates with the new image
```

**Roll back** by putting the previous tag in `.env` and repeating — old images stay in the
registry, which is the argument for immutable tags.

**Add or replace databases:** rsync them in (step 6), then `restart`. Nothing else. A
replaced file needs the restart before anyone sees the new data; until then the old inode
is still being served, which is a feature, not a race.

**Uninstall:** `docker compose -f compose.prod.yaml down`, remove the Apache block, reload
Apache. The databases are untouched — the mount is read-only.

## Configuration

Every setting is an environment variable read at startup, so nothing requires a code
change. `.env.example` is the annotated template; the same table is in the
[README](../README.md#configuration-environment-variables).

| Variable | Default | Purpose |
|---|---|---|
| `FV_DB_DIR` | `dbs` | Directory scanned for `fv.*.snps.sqlite` (availability = files present) |
| `FV_BASE_PATH` | *(empty)* | Reverse-proxy prefix (FastAPI `root_path`); empty at the domain root |
| `FV_CORS_ORIGINS` | *(empty = all)* | Comma-separated allowed browser origins |
| `FV_LOG_LEVEL` | `INFO` | `DEBUG` / `INFO` / `WARNING` / …; discovery logs at INFO, skipped databases at WARNING |

`foundervariants serve` sets `FV_DB_DIR` from `--dir` and `FV_LOG_LEVEL` from `-v`; in
Docker they come from `.env` / the compose file.

Compose-only settings (`IMAGE_NAME`, `PLATFORM`, `HOST_PORT`, `HOST_DB_DIR`) configure the
container rather than the app.

## Troubleshooting

**Intermittent "connection reset by peer", fine on retry (every app on the host).**
Investigated 2026-10-17. Not FounderVariants: the failed requests never reach Apache's
backends, and the static site root fails the same way. The front load balancer
keeps idle connections to Apache on port 80 and reuses them.
Apache closes idle connections after its default `KeepAliveTimeout` of 5 s. When the
load balancer reuses one Apache has already closed, it resets the client. Requests 3 s apart
failed 0/30; 7 s apart, 10/30. The fix is to make Apache outlast the load balancer's idle
timeout (often several minutes; check yours). As root, in a new
`/etc/httpd/conf.d/00-keepalive.conf`:
```apache
KeepAlive On
KeepAliveTimeout 400
MaxKeepAliveRequests 0
```
then `apachectl configtest && apachectl graceful`. **Applied on the production host on
2026-10-17; the 7 s-gap test went from 10/30 failures to 0/30.** Prefer `graceful` over a full restart
for config changes generally: a restart drops every connection the load balancer holds and
causes a burst of resets.

**`/health` works locally but the public URL 404s.** uvicorn logs the path it actually
received, which settles this immediately:

```bash
docker compose -f compose.prod.yaml logs --tail=30
```

`GET //genomes … 404` means the proxy is doubling the slash — see the trailing-slash rule
above. `GET /foundervariants/health … 404` means the prefix is not being stripped. No line
at all means the request never reached the container: wrong port, or the `<Location>` is not
active.

**`releases: []` or fewer than expected.** The directory is empty, the mount path is wrong,
the files are unreadable, or a file's name does not match the `fv.*.snps.sqlite` glob that
discovery uses. Compare `ls "$HOST_DB_DIR"/fv.*.snps.sqlite | wc -l` with what `/health`
reports, then run with `FV_LOG_LEVEL=DEBUG` — every skipped file is logged with the reason.
A half-copied database is the usual cause, and it is skipped rather than fatal by design.

**A newly copied database is not served.** Discovery runs at startup only. Restart.

**A database was replaced and the old data is still coming back.** Same cause, same fix —
the running process is holding the previous inode open. Restart.

**Port already in use.** Something else holds `HOST_PORT` — `ss -lntp | grep <port>`. Pick
another and update `.env` *and* the Apache `ProxyPass` together.

**The API reference at `/docs` is blank.** The vendored Scalar bundle is not being served.
Confirm with `curl -sI https://churchilllab.jax.org/foundervariants/static/scalar.js` — a 404 there
with a working `/docs` page means the proxy prefix is not reaching static assets.

**A founder-pattern query takes about a second.** That is the documented cost, not a hang:
genome-wide `private_to` / `differ` / `sdpn` filters scan the chromosome because they test
bit arithmetic on `sdpn` rather than an index. Region- or gene-scoped versions of the same
query are milliseconds. See [`ARCHITECTURE.md`](ARCHITECTURE.md#6-performance-full-358m-row-grcm39-sqlite).

**The R/qtl2 download is slow or shows no progress in a browser.** Expected — see the 4 GB
download note in step 7, and prefer the Apache `Alias` for bulk transfer.

**Image runs locally but not on the server.** Architecture mismatch — build with
`--platform linux/amd64` (or multi-arch) rather than relying on your machine's default.

**`docker build` hangs on "transferring context".** `.dockerignore` is missing or was
edited — the build context is trying to include `dbs/`.
