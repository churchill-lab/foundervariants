# FounderVariants — Architecture & Data Model

This document explains **how** the system is built and **why** it's shaped this
way. For the HTTP surface see [API.md](API.md); for the provenance of the
databases that exist see [DATABASES.md](DATABASES.md); for building, shipping and
operating the service see [DEPLOYMENT.md](DEPLOYMENT.md); for a quickstart see the
[README](../README.md).

---

## 1. What it is

Given a genomic region, gene, or rsID, return — fast — the allele carried by each
of the **8 Diversity Outbred / Collaborative Cross founder strains** at every SNP,
together with the strain-distribution pattern (SDP / SDPN), the variant's
functional consequence, and its gene. The same artifact is a drop-in **R/qtl2**
variant database.

Scope decisions (deliberate):

- **8 founders only.** Not the ~50 Sanger strains, not MUSTER's 600+. This is a
  founder-genotype tool for DO/CC QTL work.
- **SNPs first, indel-ready.** SDPN is intrinsically a SNP concept (base-4 over
  A/C/G/T). The schema/pipeline can take an indel dataset as a parallel type later.
- **No LLM/chat layer.** Query-by-location-and-filter.

The 8 founders, in fixed order (letters A–H used in displays):

| idx | strain | letter | notes |
|--:|---|:--:|---|
| 0 | A/J | A | |
| 1 | C57BL/6J | B | **reference** — not in the VCF; forced to REF |
| 2 | 129S1/SvImJ | C | |
| 3 | NOD/ShiLtJ | D | |
| 4 | NZO/HlLtJ | E | |
| 5 | CAST/EiJ | F | |
| 6 | PWK/PhJ | G | |
| 7 | WSB/EiJ | H | |

---

## 2. Encodings

Only **`sdpn`** and **`confidence`** are stored per variant; `sdp`, `alt`,
per-founder bases, per-founder confidence, and `alleles` are **derived on read**
(single source of truth: `encode.py`). This is the core compactness decision.

### SDPN — the canonical genotype (Vincent)

Each founder's *actual base* A/C/G/T maps to a base-4 digit (`A=0, C=1, G=2, T=3`).
The 8 digits, **A/J most-significant → WSB least-significant**, form a 16-bit
integer (0–65535). Two bits per founder:

```
founder i base = (sdpn >> (2*(7-i))) & 3        # 0=A_J .. 7=WSB
```

Worked example — `sdpn = 43682`:

```
43682 = 0b 10 10 10 10 10 10 00 10
          A_J C57 129 NOD NZO CAST PWK WSB
           G   G   G   G   G   G    A   G      (digit 2=G, 0=A)
```
So only **PWK** carries `A`; every other founder carries `G` (= REF, since B6=G).

Because B6 is forced to REF, **REF is always founder 1's base** (`sdpn_to_bases(sdpn)[1]`),
so `ref` itself is derivable too.

### SDP — the R/qtl2 pattern (Broman)

One bit per founder: `0` if the founder's base == REF, `1` if it differs.
**A/J = bit 0 (LSB) … WSB = bit 7 (MSB)**, integer 0–255. Ref-anchored (B6 is
always 0). Derived from `sdpn` + `ref`:

```python
sdp = calc_sdp([0 if base == ref else 1 for base in sdpn_to_bases(sdpn)])
```

For `sdpn=43682` (only PWK differs, PWK = bit 6): `sdp = 2^6 = 64`.

> **Compatibility note.** `qtl2::calc_sdp` sets the bit for the *minor* allele
> among the 8 founders; we set it for *non-REF*. These agree except at the rare
> sites where REF is the minor allele (there our `sdp` = `255 − qtl2's`), which is
> an equivalent 2-group partition for mapping. This matches the original
> `fv.2021` behavior and is verified drop-in in R.

### Display strings (the "fill strings")

- `sdp_to_partition(64)` → `ABCDEFH:G` — founders with bit 0 : founders with bit 1.
- `sdpn_to_display(sdpn)` → founders grouped by shared base, largest group first,
  e.g. `DEFH:BG:A:C`. At biallelic sites this matches the SDP partition; at
  **multiallelic** sites it shows 3–4 groups where SDP (2 groups) cannot — this is
  SDPN's whole point.
- `sdpn_to_sequence(sdpn)` → raw 8-base string `AGCTTTGT` (A/J..WSB).

### Confidence

The per-founder Sanger **FI** flag (1 = high, 0 = low) encoded exactly like SDP
(A/J = bit 0), integer 0–255. `confidence == 255` ⇒ all 8 founders high-confidence.
`conf_bits(confidence)` recovers the per-founder flags for shading in the UI.

### `csq_mask` — consequences without a giant table

A variant can hit many transcripts, so VEP emits **one consequence row per
(allele × transcript)** — ~3 per SNP, **110M rows genome-wide**. That transcript-level
table is what made the legacy DB huge, and the UI only ever shows the *top*
consequence and *filters by type*.

So instead of storing 110M rows, each variant carries a **64-bit `csq_mask`**: bit
`k` set ⇔ the variant has SO term with `bit=k` (mapping stored in `so_terms.bit`,
≤40 terms). Multi-term filtering is one bitwise test:

```sql
WHERE (csq_mask & :mask) <> 0        -- :mask = OR of the selected terms' bits
```

Verified to return results identical to the old transcript-level table. The
per-transcript breakdown is intentionally dropped (link out to Ensembl/VEP for
that); the top consequence per variant is kept in the `consequence` column.

---

## 3. Storage: one SQLite per release

### The decision (and the honest history)

We benchmarked row-store against column-store options on the real 35.8M-row data, and
researched the field (TileDB/GenomicsDB/Hail solve a many-samples problem we don't have —
we have 8 fixed founders collapsed into one 16-bit int). Findings:

- Region / gene / rsID lookups: **≤16 ms on every option tried** — SQLite's B-tree is
  genuinely excellent here.
- Genome-wide SDPN/founder-pattern scans: **~12–40 ms columnar vs ~1.4 s SQLite** (~100×).
- On disk (variants only): **under 1 GB columnar vs ~4 GB SQLite** — a row store doesn't
  compress the repetitive TEXT columns the way dictionary encoding does.

**Chosen: SQLite-only, one file per release.** Rationale: R/qtl2 can *only* read
SQLite, so a SQLite file is mandatory anyway; making it *also* the served store
means **one artifact, not two**, which is simplest to operate and ship. The cost is
that genome-wide founder-pattern scans take ~1.4 s instead of ~12 ms — acceptable
for an occasional operation; region/gene/rsID (the daily drivers) stay fast.

**SQLite is the only store in the system** — at build time, at serve time, and in R.
Nothing stages, mirrors or converts the data on its way to that file; `build.py` writes it
directly and `query.py` reads it with stdlib `sqlite3`. The one runtime dependency beyond
the standard library is FastAPI, and the one build dependency is a VCF reader.

Size reality: **~4 GB per file** (vs the 11 GB legacy DB, and vs the 7.7 GB of an
earlier two-file design that served a columnar copy alongside the R/qtl2 SQLite).
It is *smaller and one file*, but not tiny —
SQLite's row store is the floor here. An optional ~2–2.5 GB variant (int-encode
`ensembl_gene`/`consequence`, drop derivable `ref`/`alleles`/`type`) is possible but
reduces direct R-browsability; not currently applied.

### Schema (the served + R/qtl2 SQLite)

```sql
-- one row per SNP; doubles as the R/qtl2 `variants` table
CREATE TABLE variants(
  snp_id TEXT,       -- rsID when present, else chr:pos_ref_alt  (R/qtl2 id)
  chr TEXT,          -- '1'..'19','X','Y','MT'
  pos INTEGER,       -- basepairs  (R/qtl2 query func divides by 1e6)
  ref TEXT,
  alleles TEXT,      -- 'REF|alt1/alt2' among founders (R/qtl2 display)
  sdp INTEGER,       -- 0..255, A/J=bit0 (R/qtl2)
  sdpn INTEGER,      -- 0..65535, base-4 A/C/G/T, A/J=MSB (canonical)
  confidence INTEGER,-- FI bits, 255 = all high-confidence
  ensembl_gene TEXT, -- gene of the top consequence
  consequence TEXT,  -- top SO term at the variant
  csq_mask INTEGER,  -- 64-bit bitmask over SO terms (multi-term filtering)
  type TEXT          -- 'snp'
);
CREATE INDEX chr_pos ON variants(chr,pos);   -- R/qtl2 + region queries
CREATE INDEX ix_snp  ON variants(snp_id);    -- rsID lookup

CREATE TABLE genes(ensembl_id, ensembl_version, symbol, name,
                   chromosome, start_position, end_position, strand);  -- self-contained gene resolution
CREATE TABLE so_terms(term, description, accession, display, impact, score, bit); -- bit = csq_mask position
CREATE TABLE meta(meta_name, meta_value);        -- provenance (see below)
CREATE TABLE description(description, genome_build, date_created);  -- R/qtl2 metadata table
```

`create_variant_query_func()` needs `snp_id, chr, pos, alleles, sdp` — all present
with canonical names, so **no field overrides**. Everything else is passed through.

### Provenance (`meta`)

Recorded at build time so a database is self-describing and reproducible:
`release`, `grc`, `created` (UTC ISO-8601), `founders`, `source_vcf`,
`gene_annotation_file` + `gene_annotation_version` (which ensimpl/Ensembl export the
genes came from — so annotations can be refreshed), `n_variants`, `chroms`,
`schema_version`, `csq_model`.

---

## 4. Build pipeline

One pass, one output file, one dependency.

```
Sanger VCF ──cyvcf2 stream, per chromosome──▶ fv.<release>.snps.sqlite
 24 GB, dev machine only                      4.27 GB, the only artifact
```

`build.py` (streaming, per chromosome, bounded memory — nothing is staged):
1. Locate the 8 founders in the VCF header (B6 is absent → forced to REF).
2. Parse CSQ field indices **by name** (the field order differs between Ensembl 75
   and 104 — parsing by name handles both).
3. Per record, keep **homozygous** founder calls only; drop any-het/no-call (unless
   `--keepno`); keep low-confidence by default (`--keeplow`); drop records where no
   founder differs from REF. Compute `sdpn` and `confidence`.
4. While reading the record's CSQ entries — the moment every SO term for the variant is
   in hand — canonicalise each term through `DEPRECATED_TERMS` (VEP 75's
   `initiator_codon_variant` → VEP 104's `start_lost`), then accumulate `csq_mask`
   (`mask |= 1 << bit`) and track the top-scoring term and its gene. The transcript-level
   rows are counted for `meta.n_consequences` and then dropped on the floor; nothing
   materializes 110M rows.
5. Derive `alleles` and `sdp` from `sdpn` + `ref` and insert the finished row, in
   100k-row batches.
6. After the last chromosome: load self-contained `genes` (from the ensimpl TSV) and
   `so_terms` (with its bit assignment), write `meta` + `description`, then build the
   indexes.

Indexes are created last, after the bulk insert, rather than maintained during it.

Validated twice over. The pipeline reproduces the legacy `fv.2021.snps.db3` **exactly** —
35,805,856 variants with `sdpn`, derived `sdp`, and `confidence` matching on every row
(0 mismatches). And when the staging step this pipeline used to have was removed
(2026-09-01), the single-step result was diffed column-for-column against the databases
the two-step builder had produced: **0 differences** across all 12 columns, identical
`csq_mask`, byte-identical stored schema. See
[DATABASES.md](DATABASES.md#validation).

---

## 5. Serving

`api.py` is a FastAPI app-factory. At startup it scans `FV_DB_DIR` for
`fv.*.snps.sqlite`, reads each `meta`, and indexes by release. Each request opens a
**fresh read-only, immutable** SQLite connection on a worker thread
(`file:…?mode=ro&immutable=1`) — no locking, safe to fan out across workers. The
query layer (`query.py`) is pure stdlib `sqlite3` and is usable as a library
(`from foundervariants import query`). Middleware: gzip, permissive CORS (read-only).
Reverse-proxy prefix via `FV_BASE_PATH`.

Gene resolution is **self-contained** (own `genes` table) — no runtime dependency on
the sibling ensimpl service; ensimpl is only the *source* of the gene annotation TSV
at build time.

---

## 6. Performance (full 35.8M-row GRCm39, SQLite)

| query | latency |
|---|--:|
| region (2 Mb) | ~2–16 ms |
| gene / rsID | ~ms |
| multi-term consequence filter (region-scoped) | ~ms |
| genome-wide founder-pattern (`private_to`, `differ`) | ~1–1.5 s |

---

## 7. Design decisions log

- **8 founders only** — founder-genotype tool, not a strain browser.
- **SQLite-only, one file** — R/qtl2 needs SQLite anyway; one artifact beats two.
- **One store, one pass, no staging** (2026-09-01) — the build used to write an
  intermediate analytical database and convert it, which made sense while the
  transcript-level consequences were a real queried table. Once they became a bitmask,
  that step existed only to re-derive, via a 110M-row `GROUP BY`, something the VCF parse
  already had in hand and threw away. Removing it deleted a module, two dependencies and
  a 4 GB write-then-read per build, and changed the output not at all (verified
  column-for-column). The cost: a `so_terms.tsv` correction now needs a re-parse instead
  of a cheap re-convert — worth it for a vocabulary that has changed once.
- **`sdpn` canonical, everything else derived** — smaller, and eliminates a class of
  column-mismatch bugs the legacy code had (`sample_AJ` vs `sample_A_J`).
- **`csq_mask` instead of 110M transcript rows** — keeps the file to one manageable
  file and makes multi-term filtering a single bitwise test; drop-in Ensembl links
  cover the rare need for per-transcript detail.
- **Canonical R/qtl2 column names** (`snp_id`,`sdp`,…) — drop-in with no overrides,
  cleaner than the old DB which needed `id_field="variants_id"`.
- **Gene annotations baked in with provenance** — no runtime ensimpl dependency;
  refreshable via a newer TSV, tracked in `meta`.
- **Parse VCF CSQ by field name** — robust across Ensembl VEP versions.
- **One consequence vocabulary across all releases** (2026-09-01) — Ensembl renames terms,
  so a release annotated by VEP 75 and one by VEP 104 described the same event with
  different words (`initiator_codon_variant` vs `start_lost`). Filtering could not span the
  two builds: `consequence=start_lost` silently returned zero rows on GRCm38. Deprecated
  terms are now mapped to current ones at a single point in the parse, so severity, the top
  consequence and the `csq_mask` bit cannot disagree. See
  [DATABASES.md](DATABASES.md#deprecated-vep-terms-are-canonicalised-on-read) for the
  evidence behind the ranking and the one case the source annotation cannot distinguish.
