FounderVariants — server bundle
================================
Everything the server needs. It does NOT need a source checkout.

FILES HERE
  compose.prod.yaml            pulls/runs the image; never builds
  .env                         host settings — EDIT the two marked lines
  apache-foundervariants.conf  reverse-proxy block, as live on churchilllab.jax.org

SEPARATELY, the databases (~8.5 GB total) must be rsynced to HOST_DB_DIR.
They are NOT in this bundle and are never baked into the image.

FULL WALKTHROUGH: docs/DEPLOYMENT.md in the repo (steps 5-8).

THE IMAGE IS ON DOCKER HUB
  churchilllab/foundervariants:0.2.0
  linux/amd64, 55.4 MB, digest sha256:4aea2aae4f1abafaf52199cda954e2657f5d424796f904264502833d7a8214ca
  Needs schema compact-v2 databases for the gene features (it still serves older files,
  with one gene per SNP).
  compose.prod.yaml sets pull_policy: always, so `up -d` fetches it.
  If the repo is private, run `docker login` on the server first.

QUICK SEQUENCE
  1. ss -lntp | grep 28224 || echo free
  1b. docker compose -f compose.prod.yaml pull
  2. docker compose -f compose.prod.yaml up -d
  3. docker compose -f compose.prod.yaml ps        # wait for "healthy"
  4. curl -s localhost:28224/health                # expect 2 releases
  5. add the <Location> block to Apache, then:
     sudo apachectl configtest && sudo systemctl reload httpd
  6. curl -s https://churchilllab.jax.org/foundervariants/health

THREE THINGS MUST AGREE (nothing checks them for you)
  Apache <Location /foundervariants>       ==   FV_BASE_PATH in .env
  Apache ProxyPass ...:28224           ==   HOST_PORT in .env
  Trailing slashes: the <Location> path and the ProxyPass target must BOTH
  have one or BOTH lack one. Mixing them makes /foundervariants/genomes
  arrive as //genomes and 404 -- while /foundervariants itself still loads,
  so the page appears fine and only its API calls fail.

AFTER ADDING OR REPLACING A DATABASE
  Discovery runs once at startup. You must restart:
     docker compose -f compose.prod.yaml restart
  Do NOT rsync with --inplace: the app opens files with SQLite immutable=1.
  Default rsync writes a temp file and renames, which is safe.
