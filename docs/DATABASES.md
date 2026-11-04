# How every database was created

Each `fv.<release>.snps.sqlite` in `dbs/` is one Sanger Mouse Genomes release, reduced to
the 8 Diversity Outbred / Collaborative Cross founders. This file records where each one
came from — the source VCF, the gene annotation, the SO-term vocabulary, the exact command,
and what was dropped along the way.

For the encodings and the storage rationale see [`ARCHITECTURE.md`](ARCHITECTURE.md); for
building and shipping them see [`DEPLOYMENT.md`](DEPLOYMENT.md).

**Every database describes itself.** Its `meta` table is written when it is created and is
the only source of truth; no outside file is consulted to answer what a database is:

```bash
sqlite3 dbs/fv.2021.snps.sqlite "SELECT meta_name || ' = ' || meta_value FROM meta ORDER BY meta_name"
```

```
built_with = foundervariants.build streaming cyvcf2->sqlite
chroms = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X
created = 2026-08-03T00:01:46+00:00
csq_model = per-variant bitmask
founders = A_J,C57BL_6J,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ
gene_annotation_file = ensimpl.104.mm.genes.tsv
gene_annotation_version = Ensembl 104 (ensimpl.104.Mm)
grc = GRCm39
n_consequences = 109712490
n_genes = 53921
n_variants = 35805856
name = REL-2021
release = 2021
schema_version = compact-v1
served_from = single-sqlite
source_vcf = mgp_REL2021_snps.rsID.vcf.gz
variant = SNP
```

The API reads exactly this at startup and exposes it at `/genomes` and `/{release}/meta`, so
what a client sees about provenance is what the file says about itself, not a hand-kept
list.

## The `meta` contract

All of it is written by `foundervariants.build._write_meta`, from the VCF and the CLI
arguments. One row per key.

| key | meaning |
|---|---|
| `schema_version` | internal schema version — currently `compact-v2` (`compact-v1` had no `variant_genes`; its `variants` table is identical) |
| `release` / `name` | release id (`2021`) and human name (`REL-2021`) |
| `variant` | `SNP` — the type of variant in this file |
| `grc` | genome build (`GRCm39`) — **coordinates are only comparable within one** |
| `created` | UTC ISO-8601, when the VCF was parsed |
| `founders` | the 8 founder names in canonical A–H order, as the encodings assume them |
| `source_vcf` | basename of the Sanger VCF parsed |
| `gene_annotation_file` | the gene TSV baked into the `genes` table |
| `gene_annotation_version` | which Ensembl/ensimpl release that TSV came from |
| `chroms` | chromosomes actually present in the output, in build order |
| `n_genes` | rows in `genes` |
| `n_variants` | rows in `variants` |
| `n_consequences` | transcript-level consequence rows seen — **not stored**, see below |
| `n_variant_genes` | rows in `variant_genes`: one per (SNP, gene) |
| `built_with` | the pipeline shape that produced it |
| `served_from` | `single-sqlite` — this file is both the served store and the R/qtl2 file |
| `csq_model` | `per-variant bitmask + per-gene bitmask (variant_genes)` — how consequences are represented (`per-variant bitmask` before `compact-v2`) |

`n_consequences` is deliberately a record of something the file does **not** contain: the
110 M transcript-level consequence rows that were folded into each variant's `csq_mask` as
they were read. Keeping the count makes the reduction auditable — you can see what was
folded away without storing any of it.

> **Two known differences in the currently shipped files**, both from the older two-step
> builder that produced them on 2026-08-03:
>
> - They carry `n_variants` **twice** (same value). `query.meta()` builds a dict, so this
>   was never visible through the API — only in a raw `SELECT * FROM meta`.
> - Their `built_with` names a staging step the pipeline no longer has.
>
> Both clear on the next rebuild. Neither affects the data: the current builder was
> verified to reproduce these files' `variants` rows exactly (see
> [Validation](#validation)).

## Source data: the Sanger Mouse Genomes VCFs

Both databases come from the Mouse Genomes Project, downloaded to a local
mirror (`RESOURCE/Sanger/`). The URLs below are the ones recorded in
`../data/<release>.snps.meta.tsv` when the files were fetched (2023); the facts in the table
are read from the VCF headers themselves and are current.

| | REL-2021 | REL-1505 |
|---|---|---|
| directory | `REL-2112-v8-SNPs_Indels/` | `REL-1505-SNPs_Indels/` |
| SNP VCF | `mgp_REL2021_snps.rsID.vcf.gz` | `mgp.v5.merged.snps_all.dbSNP142.vcf.gz` |
| size | 24.0 GB (+ `.csi`) | 23.0 GB (+ `.tbi`) |
| download URL (as recorded) | `https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/` | `http://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/` |
| reference | `Mus_musculus.GRCm39.dna.toplevel.fa.gz` | `ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa` |
| genome build | GRCm39 | GRCm38 |
| VEP | v104, cache `104_GRCm39`, dbSNP 150, GENCODE M27 | Ensembl 75 |
| strains in the VCF | 52 | 36 |
| contigs | 20 — `1`–`19`, `X` | 22 — `1`–`19`, `X`, `Y`, `MT` |

Two consequences of that last row, both real and neither a bug:

- **GRCm39 has no Y and no MT.** The REL-2021 SNP VCF simply carries no such contigs, so
  `fv.2021.snps.sqlite` has none either. This matches the legacy `fv.2021.snps.db3`. Only
  the GRCm38 build covers Y and MT — and thinly: 3,309 Y and 603 MT variants survive
  filtering, against ~1.4–2.8 M on an autosome.
- **The `chroms` meta key is measured, not declared.** `build` intersects the requested
  chromosome list with `vcf.seqnames`, so a release that lacks a contig produces a database
  that says so rather than one with an empty chromosome.

### The founders, and the one that isn't in the file

The pipeline needs 8 founders. **Both VCFs contain only 7 of them.**

| founder | in REL-2021 | in REL-1505 |
|---|:--:|:--:|
| A/J, 129S1/SvImJ, NOD/ShiLtJ, NZO/HlLtJ, CAST/EiJ, PWK/PhJ, WSB/EiJ | ✓ | ✓ |
| **C57BL/6J** | — | — |

C57BL/6J *is* the reference genome, so it is never a sample: `build.parse_chrom` skips
founder index 1 outright and assigns it REF, with confidence 1. Everything downstream
depends on this — it is why REF is recoverable from `sdpn` alone, and why `sdp` bit 1 is
always 0 (see [`ARCHITECTURE.md`](ARCHITECTURE.md#2-encodings)).

**Do not mistake `C57BL_6NJ` for it.** Both VCFs carry a `C57BL_6NJ` column; 6NJ is a
different substrain and is *not* used. The lookup is by exact name against
`encode.FOUNDERS`, so `C57BL_6J` simply misses — which is the correct outcome here, but
means a founder renamed in a future Sanger release would silently become a no-call rather
than an error.

### The CSQ field order differs between the two releases

VEP writes its field layout into the `CSQ` INFO header, and the two releases do not agree:

```
REL-2021 (VEP 104):  Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|…
REL-1505 (VEP 75):   Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|…
```

`build._csq_indices` therefore parses the header and looks fields up **by name**, never by
position. This is the single reason one pipeline handles both releases, and it is the thing
to preserve when adding a third.

### What gets dropped, and why the counts differ from the VCF

`parse_chrom` keeps a record only if it is a clean, informative, biallelic-per-founder SNP:

1. REF and every ALT must be a single base (multi-base records are skipped before counting).
2. Every founder call must be **homozygous** — any het, half-call or missing genotype drops
   the record, unless `--keepno`.
3. Low-confidence calls (Sanger's per-sample `FI` flag) are **kept** by default; pass
   `--no-keeplow` to require all 8 founders high-confidence at build time instead. (The API
   filters on this at query time via `confidence`, which is why keeping them is the default.)
4. At least one founder must differ from REF — a site where all 8 match the reference
   carries no founder information.

**The VCF `FILTER` column is deliberately ignored.** `FILTER` is one flag for the whole
site, set to fail if *any* strain fails, and REL1505 has 36 strains. Founder quality is
already captured per founder by `FI`, stored as `confidence`. A site-level filter would
only throw away good founder calls. Checked on chr19, re-running the keep rule against the
source VCFs (the same 886,707 and 881,574 rows as the databases):

| REL1505 chr19, kept rows | rows |
|---|--:|
| `FILTER = PASS` | 758,436 |
| non-PASS, all 8 founders pass `FI` (the failure is in other strains) | 83,486 |
| non-PASS, some founder fails `FI` (already in `confidence`) | 44,785 |
| non-PASS with no strain failing `FI` | 0 |

For example, chr17:34,265,003 is `MinDP` in REL1505 only because of LEWES_EiJ; all the
founders pass. REL2021 doesn't use `FILTER` at all: every kept chr19 row is `PASS`, and its
quality information is entirely in `FI`. (Raised in an email review, 2026-10;
re-verified here.)

The result, from the build logs:

| release | SNP records seen | kept | kept | consequence rows seen |
|---|--:|--:|--:|--:|
| 2021 | 83,212,162 | 35,805,856 | 43.0% | 109,712,490 |
| 1505 | 78,772,544 | 36,218,870 | 46.0% | 86,055,840 |

The consequence rows are counted **one row per SO term**, where the source VCF `&`-joins
multiple terms into one CSQ entry. That roughly doubles the count against the legacy
database's string-joined form; only ~0.2% are true duplicates, which is not worth
de-duplicating given that none of these rows are stored (they are aggregated into
`csq_mask`).

### Multiallelic sites, and how the two releases differ

Every kept row has a clean homozygous call for all 8 founders, so these are exact counts of
sites where the founders carry 3 or 4 different bases:

| bases among the 8 founders | REL2021 (GRCm39) | REL1505 (GRCm38) |
|---|--:|--:|
| 2 | 35,628,630 | 36,037,871 |
| 3 | 177,225 | 180,678 |
| 4 | 1 | 321 |

SDP can't represent these (it records only *reference or not*); SDPN can. See `/sdp` on
the site.

**REL2021's only four-base site** is rs52115084 (chr1:192,017,092; GRCm38
chr1:192,334,784), in an intron of *Kcnh1*: `TCTTAGTG` in both releases, all 8 calls
high-confidence.

**Why REL1505 has 321 and REL2021 has 1.** 296 of the 321 have an rsID, and 197 of those
are in REL2021. At 187 of the 197, one founder carried the allele REL2021 no longer has.
At 152 of those, that founder's REL2021 call is the reference base with low confidence
(CAST 100, PWK 26, WSB 20, NZO 4, A/J 2). For example, rs233533684 (chr4:21,921,554) is
CAST = C in REL1505 and CAST = G, the reference, flagged low-confidence in REL2021. Which
release is right at these sites hasn't been investigated.

**What this means for users comparing releases:** a founder call can differ between
releases where one of them is low-confidence. The API and search page return
**all-high-confidence rows only, by default**, so rows like rs233533684's REL2021 row are
hidden unless `low_confidence=true`. Don't read REL2021's low-confidence calls as "no
call, reported as the reference" in general: across all REL2021 rows, 37% of
low-confidence founder calls are the reference base (64% in REL1505). The pattern above is
specific to these sites.

### Indels exist in both releases and are not built

Both source directories carry an indel VCF (`mgp_REL2021_indels.rsID.vcf.gz`, 5.7 GB;
`mgp.v5.merged.indels.dbSNP142.normed.vcf.gz`, 1.8 GB). Nothing consumes them yet. The
schema is indel-ready — `variants.type` exists and is `'snp'` for every current row — but
`sdpn` is intrinsically a SNP encoding (base-4 over A/C/G/T), so an indel dataset needs a
parallel genotype representation, not just a new `type` value. Note also that the recorded
metadata for the 2021 indel file says `VEP Ensembl 107` where the SNP file says 104, and its
`VARIANT` field says `SNP` — a typo in that 2023 TSV, not in the VCF.

## Gene annotation: the input that no longer regenerates

`genes` is baked in at build time from a TSV, which is what makes gene resolution
self-contained — the service has **no runtime dependency** on the sibling ensimpl project.

| database | gene TSV | recorded as | genes |
|---|---|---|--:|
| `fv.2021` | `../data/ensimpl.104.mm.genes.tsv` | `Ensembl 104 (ensimpl.104.Mm)` | 53,921 |
| `fv.1505` | `../data/ensimpl.75.mm.genes.tsv` | `Ensembl 75 (ensimpl.75.Mm)` | 38,924 |

**These TSVs are legacy exports and cannot currently be regenerated.** They were produced by
ensimpl v1/v2, whose `.db3` files still sit beside them in `../data/`. Current ensimpl
(v3.0.0) does not export a `genes.tsv` at all. Refreshing the gene annotation therefore
means either dumping the columns from an ensimpl v3 database directly
(`ensembl_id, ensembl_version, symbol, name, chromosome, start_position, end_position,
strand`) or adding an export path to ensimpl — not re-running an existing command.

Both `gene_annotation_file` and `gene_annotation_version` are recorded in `meta` precisely
so a refresh is traceable: swapping the TSV and rebuilding changes what the file says about
itself.

Matching the assembly matters more than matching the release number. `fv.1505` is GRCm38 and
uses the Ensembl 75 export; `fv.2021` is GRCm39 and uses Ensembl 104. Pairing a GRCm38
variant set with GRCm39 gene coordinates would produce a database that looks fine and
answers `?gene=` wrongly.

## The SO-term vocabulary, and why its row order is load-bearing

`../data/so_terms.tsv` is 39 terms with `term, description, accession, display, impact,
score`. It does two jobs:

- **`score` ranks severity.** The highest-scoring term at a variant becomes its
  `consequence` column (`transcript_ablation` = 39 down to the modifiers).
- **Row order assigns `bit`.** `build._load_so` numbers the terms `0..38` *in file order*,
  and that number is the bit position in every variant's `csq_mask`.

> **Reordering, inserting, or removing a row in `so_terms.tsv` silently invalidates the
> `csq_mask` of every existing database.** The mapping is stored in each file's `so_terms`
> table, so an old database keeps answering with its own bits — but a mask computed from a
> new TSV will not mean the same thing. Append new terms at the end, or rebuild every
> database.

Because the bits and the top-consequence choice are both decided during the VCF parse, any
correction to this file costs a full rebuild — the trade named in
[`DEPLOYMENT.md`](DEPLOYMENT.md#2-build-a-database).

### Deprecated VEP terms are canonicalised on read

We build releases annotated by very different VEP versions from one vocabulary, and
Ensembl renames consequence terms. `build.DEPRECATED_TERMS` maps the old name to the
current one during the parse:

| old term (as the VCF emits it) | canonicalised to | renamed in |
|---|---|---|
| `initiator_codon_variant` (SO:0001582) | `start_lost` (SO:0002012) | Ensembl ~87 |

VEP 75 (GRCm38 / REL-1505) emits `initiator_codon_variant`; VEP 104 (GRCm39 / REL-2021)
emits `start_lost` for the same event. Canonicalising means **one vocabulary spans every
release**, so `consequence=start_lost` matches on both builds. Without it the same
biological event was unfindable under the name the current Ensembl docs use — it returned
zero rows on GRCm38, silently. The mapping is applied at a single point in
`build.parse_chrom`, so severity ranking, the top-consequence choice and the `csq_mask` bit
can never disagree about which term they mean.

`initiator_codon_variant` has been **removed** from `so_terms.tsv` accordingly: after
canonicalisation no variant can carry it, and leaving it in would put a filter entry in the
UI that always returns nothing. Removing it was only safe because it was the **last** row
(bit 39) — deleting any other row would have shifted every subsequent term's `csq_mask`
bit. `start_lost` remains at bit 6, so no existing mask changed meaning.

> **Caveat, and it is a real one.** Modern Ensembl split the old term into `start_lost`
> (HIGH) and `start_retained_variant` (LOW). VEP 75 could not distinguish them, so a
> GRCm38 variant whose start codon is rescued by an alternative in-frame ATG is called
> `start_lost` here. That is a limit of the source annotation, not of the mapping — VEP 75
> simply did not record the distinction.

#### Why `start_lost`, and not a score of its own

The junk row this replaced had `score = -1`, `impact = NULL`, `accession = NULL` — visible
in the legacy `fv.1505.snps.db3` and the cause of the 172 wrong rows described under
[Validation](#validation). A first patch gave it `21 / LOW`, which was invented rather than
sourced, and was wrong in a way worth recording:

- **Ensembl 75 — the very VEP version that annotated REL-1505 — ranked
  `initiator_codon_variant` 7th most severe**, immediately after `stop_lost` and
  immediately before `inframe_insertion`. That is exactly the slot modern Ensembl gives
  `start_lost` (score 33, HIGH). e75 had no IMPACT column at all; its ordering *was* the
  severity, so no Ensembl release has ever rated this term LOW.
- At `21` it sat **below `start_retained_variant` (22)** — ranking "start codon destroyed"
  as less severe than "start codon changed but retained" — and tied with
  `stop_retained_variant`, a synonymous-class event.

Sources: [Ensembl 75 predicted data](http://feb2014.archive.ensembl.org/info/genome/variation/predicted_data.html)
· [Ensembl current consequences](https://jun2026.archive.ensembl.org/info/genome/variation/prediction/predicted_data.html)
(the old term no longer appears at all).

## The two databases that exist

| | `fv.2021.snps.sqlite` | `fv.1505.snps.sqlite` |
|---|---|---|
| release / build | 2021 / **GRCm39** | 1505 / **GRCm38** |
| variants | 35,805,856 | 36,218,870 |
| chromosomes | 1–19, X | 1–19, X, Y, MT |
| genes | 53,921 (Ensembl 104) | 38,924 (Ensembl 75) |
| SO terms | 40 | 40 |
| file size | 4.27 GB | 4.27 GB |
| `created` | 2026-08-03T00:01:46Z | 2026-08-03T00:55:07Z |
| parse time | 48.1 min | 45.9 min |

To reproduce either one with the current builder:

```bash
uv run foundervariants -v build \
  --vcf ~/work/RESOURCE/Sanger/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.rsID.vcf.gz \
  --genes-tsv ../data/ensimpl.104.mm.genes.tsv --so-tsv ../data/so_terms.tsv \
  --release 2021 --grc GRCm39 --gene-version "Ensembl 104 (ensimpl.104.Mm)" \
  --out dbs

uv run foundervariants -v build \
  --vcf ~/work/RESOURCE/Sanger/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
  --genes-tsv ../data/ensimpl.75.mm.genes.tsv --so-tsv ../data/so_terms.tsv \
  --release 1505 --grc GRCm38 --gene-version "Ensembl 75 (ensimpl.75.Mm)" \
  --out dbs
```

The two files currently in `dbs/` were produced on 2026-08-03 by an earlier two-step
builder that staged the parse before writing the SQLite. That staging step was removed on
2026-09-01 after being shown to add nothing (see [Validation](#validation)); the commands
above are the current, single-step equivalent and produce the same `variants` rows.

A stale `build_cache/` directory may still hold that old builder's staging files (~7.3 GB).
Nothing reads them any more — delete the directory.

### Validation

`fv.2021.snps.sqlite` was checked row-for-row against the legacy `fv.2021.snps.db3`
(11 GB, the database this project replaces): **35,805,856 variants, with `sdpn`, derived
`sdp`, and `confidence` matching on every row — 0 mismatches.** The `csq_mask` filter was
separately verified to return the same variant sets as the old transcript-level
`consequences` table.

Both files were then confirmed drop-in in real R at full scale —
`create_variant_query_func()` with **no field overrides**, `chr11:100–102 Mb` returning
~29 k rows with `pos` in Mbp and `sdp` in 1–253, feeding `index_snps()`.

`fv.1505` has no legacy counterpart to diff against at this level; it was validated
structurally (same schema, same encodings, R/qtl2 drop-in) and by the same build code that
reproduced 2021 exactly.

**The single-step builder was validated against those files in turn** (2026-09-01), which
is how the removed staging step was shown to be redundant. Rebuilding chr19 of release 2021
and diffing it against `dbs/fv.2021.snps.sqlite`:

| check | result |
|---|---|
| rows | 881,574 both sides |
| `seen` / `kept` / consequence counts | 1,929,721 / 881,574 / 3,089,283 — identical to the original build log |
| all 12 columns, both directions (`EXCEPT`) | **0 differences** |
| `csq_mask` mismatches | 0 |
| `sdp` mismatches | 0 |
| `so_terms` (term, bit, score) | 0 differences |
| `genes` | 0 differences |
| stored schema text (`.schema`) | byte-identical |
| `create_variant_query_func()`, no overrides | 2,034 rows for `19:5.0–5.2 Mb`, `pos` in Mbp, `sdp` 4–233, `index_snps()` → 17 distinct SDPs |

### The full-scale rebuild (2026-09-01)

Both releases were then rebuilt end-to-end with the single-step builder into `dbs_verify/`
and compared against the production files, every row on all 12 columns:

| | release 2021 | release 1505 |
|---|--:|--:|
| rows, production / rebuilt | 35,805,856 / 35,805,856 | 36,218,870 / 36,218,870 |
| key pairs matched 1:1, no orphans either side | ✓ | ✓ |
| stored schema text | identical | identical |
| **rows differing on any column** | **0** | **172** (0.0005%) |
| of which genotype columns (`ref`/`alleles`/`sdp`/`sdpn`/`confidence`) | 0 | **0** |

**No genotype data changed anywhere.** Every difference is in `consequence` and `csq_mask`;
the founder alleles and the R/qtl2 columns are identical in both releases.

The counting method was cross-checked against a two-directional whole-table `EXCEPT` on an
earlier rebuild: it reported the differing-row count in *each* direction (the expected
signature of rows whose key matches but a column differs) and `0` for 2021, agreeing with
the indexed merge join the comparison script now uses. Plus, on the rebuilt files: all 42
chromosome×release
aggregate fingerprints matching (row count, `count(DISTINCT snp_id)`, and sums of
`pos`/`sdp`/`sdpn`/`confidence`/`csq_mask`), and R/qtl2 `create_variant_query_func()`
verified with no field overrides at four windows including `chr11:100–102 Mb` (29,151
rows) and all 603 rows of 1505's chrMT.

**The 172 differing rows are a fix, not a regression**, and they are all the same bug:
every one of them is a variant that destroys a start codon.

| production reported | rebuilt reports | rows |
|---|---|--:|
| `(NULL)` | `start_lost` | 55 |
| `intron_variant` | `start_lost` | 34 |
| `missense_variant` | `start_lost` | 23 |
| `upstream_gene_variant` | `start_lost` | 16 |
| `non_coding_transcript_exon_variant` | `start_lost` | 14 |
| `5_prime_UTR_variant` | `start_lost` | 12 |
| `downstream_gene_variant` | `start_lost` | 8 |
| `splice_region_variant` | `start_lost` | 5 |
| `3_prime_UTR_variant` | `start_lost` | 2 |
| `NMD_transcript_variant` | `start_lost` | 2 |
| `stop_gained` (unchanged; only `csq_mask` moved) | `stop_gained` | 1 |

All 172 carried the deprecated term, and the mask migration is exact: production had 172
variants with bit 39 set and none with bit 6; the rebuild has 172 with bit 6 and none with
bit 39. 171 of them now report `start_lost` as their *top* consequence — the one exception
is correctly outranked by `stop_gained` (36 > 33).

The cause was the junk `so_terms.tsv` row. Its `score = -1` made the term unable to win the
top-consequence comparison at all (`sc > max_score`, with `max_score` initialised to `-1`).
Correcting the TSV and re-converting could not repair it: **the top consequence is chosen
during the parse, so the old two-step pipeline baked the broken ranking into its
intermediate and only the `so_terms` table received the correction.** The shipped 1505 file
therefore had a `so_terms` table that disagreed with its own `consequence` column.

Release 2021 is untouched because the term had been **renamed** — VEP 75 (GRCm38) emits
`initiator_codon_variant`, VEP 104 (GRCm39) emits `start_lost` — which is also the second,
sharper defect this uncovered. The two builds did not share a vocabulary, so no query could
span them:

| query | release 2021 (GRCm39) | release 1505 (GRCm38) |
|---|---|---|
| `consequence=start_lost`, before | matches | **silent zero** |
| `consequence=initiator_codon_variant`, before | silent zero | matches |
| either query, after | `start_lost` matches; the retired name returns `422 unknown consequence` | same |

Asking for start-codon-destroying variants returned a **silent zero** on GRCm38 — the same
question, unanswerable on one of the two builds. Now the query works on both, and the
retired name fails loudly rather than quietly.

(The `consequence=` filter tests `csq_mask`, so it matches a variant carrying the term on
*any* transcript — 172 variants in 1505. That is a broader set than the 171 whose *top*
consequence is `start_lost`, which is the `consequence` column shown in results.)

**Getting the severity right is what made the fix complete.** An interim patch scored the
term `21 / LOW`, which repaired only 143 rows: at 21 the term still lost to
`missense_variant` (29) and `splice_region_variant` (27). At its correct rank — `start_lost`
at 33 — it wins those too, which is the 23 + 5 rows that take the count from 143 to 171.
A variant that abolishes the start codon was being reported as a missense variant.

Reproduce a single chromosome with `--chroms 19 --out <scratch>` and compare against
`dbs/`, or do it at full scale with the two scripts kept for exactly this purpose:

```bash
./scripts/rebuild_verify.sh          # rebuild BOTH releases into dbs_verify/ (~100 min)
./scripts/compare_builds.sh          # diff dbs_verify/ against dbs/, read-only on both
./scripts/compare_builds.sh 2021     # or just one release
```

`rebuild_verify.sh` never writes to `dbs/`, so production keeps serving throughout.
`compare_builds.sh` does an indexed merge join on `(chr, pos, snp_id)` across all 12
columns plus a key-coverage check: when `matched 1:1` equals both row counts and
`DIFFERING` is 0, the two `variants` tables are identical as sets of rows. (A
two-directional whole-table `EXCEPT` is equivalent but has to sort 36 M twelve-column
tuples — it ran over an hour before being abandoned in favour of the join.)

**Before swapping: `matched 1:1` must equal both row counts, and every `DIFFERING` row
must be one you have explained and accepted.** Physical row order is deliberately not
compared — see the note below. Expected `meta` differences are only `created` (a new
timestamp) and, for files predating 2026-09-01, `built_with` and the duplicate
`n_variants` row.

**The files are not byte-identical, by design.** The old builder's convert step ordered
rows by `chr` as TEXT, so production stores them `1, 10, 11, …, 19, 2, …`; the single-step
builder writes in canonical order `1, 2, 3, …, 19, X`. Same rows, different physical
order, which is why 2021's rebuild differs by 8 KB of page packing (1505's happens to land
on exactly the same size). Nothing depends on physical order: every query and R/qtl2's own
range query go through the `chr_pos` index, and nothing in the codebase does
`ORDER BY rowid`.

### The compact-v2 rebuild: every gene per SNP (2026-11-04)

Both releases were rebuilt to add `variant_genes`, the per-gene consequences that
`variants` folds into one gene (see
[ARCHITECTURE.md](ARCHITECTURE.md#variant_genes--every-gene-at-a-snp-schema-compact-v2)).
The requirement was that R/qtl2 see **exactly** the same data, so the rebuilt files were
compared against production with `scripts/compare_builds.sh`, checked with
`scripts/check_variant_genes.py`, and read through R/qtl2 itself with
`scripts/check_rqtl2_same.R`:

| | release 2021 (GRCm39) | release 1505 (GRCm38) |
|---|--:|--:|
| `variants` rows, production / rebuilt | 35,805,856 / 35,805,856 | 36,218,870 / 36,218,870 |
| key pairs matched 1:1, no orphans | ✓ | ✓ |
| **`variants` rows differing on any of the 12 columns** | **0** | **0** |
| stored schema of every production object (`variants`, `chr_pos`, …) | identical | identical |
| objects added | `variant_genes`, `gi_id` | `variant_genes`, `gi_id` |
| `genes`, `so_terms` | 0 diffs | 0 diffs |
| `create_variant_query_func()`, no overrides: data frames, old vs new | identical (chr1, chr2:50–60, chr11:100–102, chr19, chrX; 5.1 M rows) | identical (same windows plus chrY, chrMT; 5.2 M rows) |
| `index_snps()` on chr11:100–102 Mb | 29,151 SNPs → 204 index SNPs | 28,907 → 226 |
| `variant_genes` rows | 24,812,793 | 21,041,348 |
| distinct genes in it | 53,650 | 41,773 |
| SNPs in 2+ gene bodies | 994,882 (2.8%) | 681,452 (1.9%) |
| SNPs with 2+ genes, counting 5 kb flanks | 3,685,106 (10.3%) | 2,563,444 (7.1%) |
| per-gene invariants (`check_variant_genes.py`, all 5) | ✓ | ✓ |
| file size, compact-v1 → v2 | 4.27 → **5.28 GB** | 4.27 → **5.13 GB** |
| parse time | 57.5 min | 56.2 min |

`meta` differs only where it should: `created`, `schema_version` (`compact-v2`),
`csq_model`, and the new `n_variant_genes`.

Query speed on the full GRCm39 file, in-process, is unchanged for normal searches (10–35 ms
for a gene or a 2,000-row window); the 25,000-row maximum takes ~0.6 s, of which ~0.15 s
is attaching genes. Over HTTP, the larger responses first made `/variants` ~40% slower
(FastAPI's JSON encoder); serializing with the stdlib encoder instead made it ~45% faster
than before the change (2,000 rows: ≈78 → 43 ms). That depends on the `gi_id` index: without it each symbol lookup scans `genes`, and
the same query took 71 s on a chr19 prototype.

**A gap this made visible, not one it caused:** VEP annotated some genes that the baked-in
ensimpl gene list lacks, so they have no symbol and show as Ensembl IDs. In `variant_genes`
that's 1,428 gene IDs on GRCm39 and 3,747 on GRCm38 (Ensembl 75's list is the thinner one).
Production already had it in `variants.ensembl_gene`: 885 genes on 63,198 SNPs (GRCm39) and
3,655 genes on 667,716 SNPs (GRCm38) show only an ID today.

### Why 4.27 GB, and what it is not

(Written for schema `compact-v1`. `compact-v2` adds `variant_genes`, about 1 GB per release;
the reasoning below is unchanged.)

SQLite's row store is the floor here, and the honest history is worth recording because the
size was mis-estimated twice during design — ~1.3 GB was predicted, ~4 GB is real. A row
store does not compress the repetitive TEXT columns (`ensembl_gene`, `consequence`, `ref`,
`alleles`, `type`) the way a column store would; measured on the real data, the same
variants in a columnar format come to well under 1 GB.

That cost was accepted deliberately. R/qtl2 can only read SQLite, so a SQLite file has to
exist regardless — and once it exists, making it *also* the served store means **one
artifact instead of two**. An earlier design that served a columnar copy alongside the
R/qtl2 SQLite came to 7.7 GB per release for the pair. Against the 11 GB legacy database,
what ships is both smaller and a single file.

A further reduction to ~2–2.5 GB is available — int-encode `ensembl_gene` and `consequence`,
drop the derivable `ref`, `alleles` and `type` — and has **not** been applied, because it
would make the file materially less browsable from R, which is half its purpose.

## Adding a release

A new build records its own metadata; there is nothing else to update:

```bash
uv run foundervariants -v build --vcf <new-snps.vcf.gz> \
  --genes-tsv ../data/<matching-assembly-genes.tsv> --so-tsv ../data/so_terms.tsv \
  --release <id> --grc <GRCm3x> --gene-version "<Ensembl NN (source)>" --out dbs
```

Before you start, check the four things that have actually broken between releases:

1. **Founder names in the VCF header** must match `encode.FOUNDERS` exactly — a rename
   becomes a silent no-call, not an error.
2. **The `CSQ` header** must name its fields; parsing is by name, so a layout change is
   fine but a missing `Consequence`/`Gene`/`Feature`/`Allele` field is not.
3. **The gene TSV must match the new assembly**, not merely be newer.
4. **`so_terms.tsv` must only ever grow at the end** — see the warning above.

Then verify against `meta`, and — since half the file's job is outside the API — open it in
R once. Both checks are in
[`DEPLOYMENT.md`](DEPLOYMENT.md#verify-the-build).

Drop the finished file into the serving directory and restart; discovery is a directory scan
for `fv.*.snps.sqlite` with no manifest to edit.
