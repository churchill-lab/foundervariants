# syntax=docker/dockerfile:1

# Image for the FounderVariants HTTP API. It contains only the app and its two runtime
# dependencies (fastapi, uvicorn) — the databases are NOT baked in. They are mounted
# read-only at runtime (see compose.yaml), so one image serves whatever releases are
# present on the host. This matters more here than for most services: each database is
# ~4 GB, so baking them in would mean a ~9 GB image per release set.
#
# The build extra (cyvcf2) is deliberately NOT installed. Building a release parses a
# 24 GB VCF for ~50 minutes and is a development-machine job; the served path is pure
# stdlib sqlite3.
FROM python:3.12-slim

# uv drives the install so the container's dependency versions match uv.lock exactly.
# It is MOUNTED per-RUN rather than COPYed in: the uv binary is ~55 MB and is needed only
# at build time, so copying it would leave it in a layer and account for two thirds of the
# finished image. With the mount, the runtime image contains no build tooling at all.
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    # Install into the system prefix so `uvicorn` is on PATH (no `uv run` needed at runtime).
    UV_PROJECT_ENVIRONMENT=/usr/local

WORKDIR /app

# Install dependencies first, in their own layer, so they stay cached until the lockfile
# changes — editing source below does not re-resolve or re-download the dependency tree.
COPY pyproject.toml uv.lock README.md ./
RUN --mount=from=ghcr.io/astral-sh/uv:latest,source=/uv,target=/usr/local/bin/uv \
    uv sync --frozen --no-install-project --no-dev

# Then add the source and install the project itself. src/ carries the vendored Scalar
# bundle (web/static/scalar.js, ~3.6 MB) that renders /docs without a CDN.
COPY src/ ./src/
RUN --mount=from=ghcr.io/astral-sh/uv:latest,source=/uv,target=/usr/local/bin/uv \
    uv sync --frozen --no-dev

# Databases are mounted here read-only at runtime; the app reads FV_DB_DIR.
ENV FV_DB_DIR=/data/dbs

EXPOSE 8000

# Two worker processes by default: every request opens its own read-only + immutable
# SQLite connection, which parallelizes cleanly across processes. Tune --workers for the
# host (or override the command in compose.yaml).
#
# uvicorn's access log is left ENABLED on purpose. The app has no request-logging
# middleware of its own, and the one line uvicorn prints per request — including the path
# as the container actually received it — is what diagnoses reverse-proxy prefix problems.
CMD ["uvicorn", "foundervariants.api:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "2"]
