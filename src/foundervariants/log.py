"""Timing-first logging (mirrors the ensimpl pattern)."""
from __future__ import annotations
import logging
import time
from contextlib import contextmanager

_LOG = logging.getLogger("foundervariants")


def configure(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
    )


def get_logger() -> logging.Logger:
    return _LOG


@contextmanager
def step(msg: str):
    """Log a step with its wall-clock duration."""
    t0 = time.perf_counter()
    _LOG.info("%s ...", msg)
    try:
        yield
    finally:
        _LOG.info("%s (%.1fs)", msg, time.perf_counter() - t0)
