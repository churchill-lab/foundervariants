<!-- Rendered below the calculator on the site's /sdp page (web/sdp.html). Relative links
     resolve against the site's base path. The numbers here are real database rows; the
     conversions are web/static/sdpcalc.js, tested against encode.py. -->
## How SDP and SDPN are calculated

Both numbers describe **which of the 8 founders carry which allele** at a SNP. They differ
in how much they keep: SDP records only *same as the reference or not*, SDPN records each
founder's actual base.

### The founders and their letters

| letter | A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|---|
| strain | A/J | C57BL/6J | 129S1/SvImJ | NOD/ShiLtJ | NZO/HlLtJ | CAST/EiJ | PWK/PhJ | WSB/EiJ |

**B (C57BL/6J) is the reference.** The mouse reference genome is B6, so B's base is the
reference base at every SNP in this database.

### SDP: strain distribution pattern

The encoding R/qtl2 uses. Each founder gets **one bit: 0 if its base matches the reference,
1 if it differs**. The bits are weighted by founder, **A/J as the lowest**:

| founder | A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|---|
| weight | 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 |

**SDP = the sum of the weights of the founders that differ**, a number from 0 to 255.
Because B is the reference, B's bit is always 0 here.

The page shows it as two groups, **reference : different**. Example, rs29358506 (the peak
of a DO GWAS on chr10): A/J carries A and the other seven carry C, so only A differs.

> SDP = **1** &nbsp;·&nbsp; shown as `BCDEFGH:A`

SDP can't tell *which* non-reference base a founder has. At **rs237124876** (chr11, GRCm39)
PWK carries A and WSB carries T, both differing from the reference C:

> SDP = 64 + 128 = **192** &nbsp;·&nbsp; shown as `ABCDEF:GH`, as if G and H shared an allele

R/qtl2's own founder SNP files may set 1 for the *minor* allele instead of the
non-reference one. A site's SDP there can then be 255 minus the one here. Both split the
founders the same way.

### SDPN: SDP with nucleotides

Each founder's actual base becomes **one base-4 digit: A = 0, C = 1, G = 2, T = 3**. The 8
digits make one number from 0 to 65535, with **A/J as the most significant digit**:

| founder | A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|---|
| weight | 16384 | 4096 | 1024 | 256 | 64 | 16 | 4 | 1 |

**SDPN = the sum of each founder's digit times its weight.** Note the founder order is the
**opposite of SDP's**: A/J is SDP's *lowest* bit but SDPN's *highest* digit.

rs29358506, bases `A C C C C C C C` (digits 0 1 1 1 1 1 1 1):

> SDPN = 0×16384 + 1×(4096 + 1024 + 256 + 64 + 16 + 4 + 1) = **5461** &nbsp;·&nbsp; shown as `BCDEFGH:A`

rs237124876, bases `C C C C C C A T` (digits 1 1 1 1 1 1 0 3):

> SDPN = 16384 + 4096 + 1024 + 256 + 64 + 16 + 0×4 + 3×1 = **21843** &nbsp;·&nbsp; shown as `ABCDEF:G:H`

The SDPN display groups the founders that share a base, **largest group first**. Unlike
SDP it keeps G and H apart. Hover over it in the results to see the 8 bases.

### Converting between them

- **SDPN → everything** (the bases, SDP, the groupings): always possible, because SDPN holds
  every founder's base and B's base is the reference.
- **SDP → SDPN: not possible on its own.** SDP doesn't record which bases are involved. The
  calculator keeps the bases it already has where they fit, and picks the rest.
- **A grouping like `A:BC:DEFGH` → SDPN:** the same; it says which founders share a base,
  not which base. B's group keeps the reference base.

### What SDPN makes possible

The founder-pattern filters on the [search page](./) are simple arithmetic on SDPN, which
is why they are fast:

- **private to** a founder: that founder's digit differs from all 7 others
- **SDPN =** a number: exactly that pattern, nucleotides included
- **differ** (API only): two founders' digits differ; see the [API reference](docs)
