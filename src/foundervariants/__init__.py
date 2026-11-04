"""FounderVariants — DO/CC founder SNP genotypes as a fast local service.

The read side (`foundervariants.query`) is usable as a library against a built
`fv.<release>.snps.sqlite` with nothing installed at all — it is pure stdlib
`sqlite3`. The VCF build pipeline (the `build` extra) is a separate concern.
"""

__version__ = "0.2.0"
