"""foundervariants CLI: build the databases, serve the API."""
from __future__ import annotations
import argparse
import os
import sys

from . import log


def _build(a) -> int:
    from .build import build  # cyvcf2 only needed here
    build(
        a.vcf, a.chroms.split(",") if a.chroms != "all" else ["all"],
        a.out, a.genes_tsv, a.so_tsv,
        name=a.name, release=a.release, grc=a.grc, gene_version=a.gene_version,
        keeplow=a.keeplow, keepno=a.keepno,
    )
    return 0


def _serve(a) -> int:
    import uvicorn
    os.environ["FV_DB_DIR"] = a.dir
    # The app re-reads this at import time, which is what --reload subprocesses see.
    os.environ.setdefault("FV_LOG_LEVEL", "DEBUG" if a.verbose else "INFO")
    uvicorn.run("foundervariants.api:app", host=a.host, port=a.port, reload=a.reload)
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="foundervariants", description=__doc__)
    p.add_argument("-v", "--verbose", action="count", default=0)
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="build a release's SQLite (served + R/qtl2) from a Sanger VCF")
    b.add_argument("--vcf", required=True)
    b.add_argument("--genes-tsv", required=True, help="ensimpl-exported gene annotation TSV")
    b.add_argument("--so-tsv", required=True)
    b.add_argument("--out", default="dbs", help="output directory (default: dbs)")
    b.add_argument("--chroms", default="all", help="'all' or comma list, e.g. 1,2,X")
    b.add_argument("-r", "--release", required=True, help="release id, e.g. 2021")
    b.add_argument("--grc", required=True, help="assembly, e.g. GRCm39")
    b.add_argument("--name", default=None, help="human name; default REL-<release>")
    b.add_argument("--gene-version", default="", help="gene annotation source/version for meta")
    b.add_argument("--keepno", action="store_true", help="keep no-call/het records")
    b.add_argument("--no-keeplow", dest="keeplow", action="store_false",
                   help="drop low-confidence records (default keeps them)")
    b.set_defaults(func=_build)

    s = sub.add_parser("serve", help="serve the API + web UI")
    s.add_argument("--dir", default=os.environ.get("FV_DB_DIR", "dbs"))
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--port", type=int, default=8000)
    s.add_argument("--reload", action="store_true")
    s.set_defaults(func=_serve)

    a = p.parse_args(argv)
    if getattr(a, "name", None) is None and getattr(a, "release", None) is not None:
        a.name = f"REL-{a.release}"
    log.configure("DEBUG" if a.verbose else "INFO")
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
