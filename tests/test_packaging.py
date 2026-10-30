"""Properties of the package that deployment depends on.

These are cheap and they guard decisions that are otherwise only visible when a
container is already built and running.
"""
from __future__ import annotations

import subprocess
import sys


def test_serving_does_not_import_the_vcf_reader():
    """The image installs no build extra, so importing the app must not need cyvcf2.

    Run in a subprocess: an in-process check would pass simply because another test
    already imported it.
    """
    code = (
        "import sys, foundervariants.api;"
        "bad=[m for m in ('cyvcf2',) if m in sys.modules];"
        "print('LEAKED:'+','.join(bad) if bad else 'CLEAN')"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "CLEAN", out.stdout


def test_query_layer_is_stdlib_only():
    """`from foundervariants import query` must work with nothing installed but stdlib."""
    code = (
        "import sys, foundervariants.query;"
        "third=[m for m in sys.modules if m.split('.')[0] in "
        "('fastapi','starlette','uvicorn','cyvcf2')];"
        "print('LEAKED:'+','.join(sorted(set(third))) if third else 'CLEAN')"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "CLEAN", out.stdout


def test_cli_exposes_exactly_build_and_serve():
    out = subprocess.run([sys.executable, "-m", "foundervariants.cli", "--help"],
                         capture_output=True, text=True)
    assert "{build,serve}" in out.stdout, out.stdout


def test_verbose_is_a_top_level_flag_only():
    """`foundervariants -v build ...` is the valid form, not `build ... -v`.

    Worth pinning because the natural guess is the other way round — the docs had it
    wrong once. Only `--help` is invoked here, so a regression can never start a server.
    """
    top = subprocess.run([sys.executable, "-m", "foundervariants.cli", "--help"],
                         capture_output=True, text=True)
    assert "--verbose" in top.stdout, top.stdout

    for sub in ("build", "serve"):
        out = subprocess.run([sys.executable, "-m", "foundervariants.cli", sub, "--help"],
                             capture_output=True, text=True)
        assert "--verbose" not in out.stdout, (
            f"-v became a {sub} flag; the README and DEPLOYMENT.md say it is top-level")


def test_version_is_exposed():
    import foundervariants
    assert foundervariants.__version__
