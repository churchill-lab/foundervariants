"""
Streaming Sanger VCF -> the single compact SQLite (served store AND R/qtl2 database).

One pass, one output file, one dependency (cyvcf2). Whole-genome capable: each
chromosome streams through in batches, so memory stays bounded regardless of release
size.

`csq_mask` is accumulated per variant *while the VCF's CSQ field is being read* — the
point at which every SO term for that variant is already in hand. The transcript-level
consequence rows are therefore never materialized; only their count is kept, as
`n_consequences` in `meta`, so the reduction stays auditable.

Legacy semantics replicated exactly (validated vs fv.2021.snps.db3):
  B6 forced to REF; homozygous-only; drop any-het/no-call (unless keep_no_call);
  keep low-confidence by default; drop records where no founder differs from REF.
"""
from __future__ import annotations

import csv
import datetime as _dt
import itertools
import os
import sqlite3

from . import encode as enc
from . import schema
from .log import get_logger, step

LOG = get_logger()
CANONICAL_CHROMS = [str(i) for i in range(1, 20)] + ["X", "Y", "MT"]
BATCH = 100_000

# Ensembl has renamed consequence terms over the years, and we build releases annotated
# by very different VEP versions from the same vocabulary. VEP 75 (GRCm38/REL-1505) emits
# `initiator_codon_variant`; VEP 104 (GRCm39/REL-2021) emits `start_lost` for the same
# event, and Ensembl has since dropped the old term from its table entirely.
#
# Canonicalising on read means ONE vocabulary spans every release: `consequence=start_lost`
# matches on GRCm38 and GRCm39 alike. Without it the same biological event is unfindable
# under the name the current docs use, silently returning zero rows on the older build.
#
# Caveat: modern Ensembl split the old term into `start_lost` (HIGH) and
# `start_retained_variant` (LOW). VEP 75 could not distinguish them, so a variant whose
# start codon is rescued by an alternative in-frame ATG is called `start_lost` here. That
# is a limit of the source annotation, not of the mapping.
DEPRECATED_TERMS = {
    "initiator_codon_variant": "start_lost",   # renamed in Ensembl ~87
}


def canonical_term(term: str) -> str:
    """Map a deprecated VEP consequence term to the current one (identity if current)."""
    return DEPRECATED_TERMS.get(term, term)

_INSERT_VARIANT = f"INSERT INTO variants VALUES ({','.join('?' * 12)})"


def _load_so(path):
    """SO vocabulary -> (rows, term->severity score, term->csq_mask bit).

    Bits are assigned by ROW ORDER in the TSV, and are stored in each database's
    so_terms table. Reordering the file changes what every existing csq_mask means:
    append new terms at the end.
    """
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return (rows,
            {r["term"]: int(r["score"]) for r in rows},
            {r["term"]: i for i, r in enumerate(rows)})


def _load_genes(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _csq_indices(vcf):
    """Locate the CSQ subfields we need BY NAME.

    VEP writes its field layout into the header and the layout differs between
    releases (Ensembl 75 puts Consequence 5th, Ensembl 104 puts it 2nd), so position
    is never assumed.
    """
    desc = vcf.get_header_type("CSQ")["Description"]
    fmt = desc.split("Format:")[-1].strip().strip('"').split("|")
    idx = {name: i for i, name in enumerate(fmt)}
    return idx["Allele"], idx["Consequence"], idx["Gene"]


def parse_chrom(vcf, chrom, fidx, ci, so2score, termbit, keeplow, keepno, stats):
    """Yield one ready-to-insert `variants` row per kept SNP on `chrom`, in position order.

    `stats` is updated in place with 'seen', 'kept' and 'csq' once the generator is
    exhausted.
    """
    i_allele, i_csq, i_gene = ci
    need = max(i_allele, i_csq, i_gene)
    seen = kept = ncsq = 0
    for rec in vcf(chrom):
        ref, alts = rec.REF, rec.ALT
        if len(ref) != 1 or any(len(a) != 1 for a in alts):
            continue
        seen += 1
        gts = rec.genotypes
        fi_arr = rec.format("FI")
        bases, ints, fis, no_call = [], [], [], False
        for k, fi_i in enumerate(fidx):
            if k == 1:  # C57BL_6J = the reference; never a sample in the VCF
                bases.append(ref); ints.append(0); fis.append(1); continue
            a0, a1, _ = gts[fi_i]
            fis.append(1 if (fi_arr is not None and int(fi_arr[fi_i][0]) == 1) else 0)
            if a0 < 0 or a1 < 0 or a0 != a1:  # missing / half / het
                bases.append("."); ints.append(0); no_call = True
            else:
                bases.append(ref if a0 == 0 else alts[a0 - 1]); ints.append(1 if a0 > 0 else 0)
        if (no_call and not keepno) or (sum(fis) < 8 and not keeplow) or sum(ints) == 0:
            continue

        sdpn = enc.bases_to_sdpn(bases)
        founder_alt = {b for b in bases if b != ref and b != "."}
        mask, max_csq, max_score, gene_of_max = 0, None, -1, None
        for entry in (rec.INFO.get("CSQ") or "").split(","):
            p = entry.split("|")
            if len(p) <= need:
                continue
            allele, csq, gene = p[i_allele], p[i_csq], p[i_gene]
            if not csq or allele not in founder_alt:
                continue
            for raw_term in csq.split("&"):
                ncsq += 1
                # One canonicalisation point, so severity ranking, the top-consequence
                # choice and the csq_mask bit all agree on the same term.
                term = canonical_term(raw_term)
                sc = so2score.get(term, -1)
                if sc > max_score:
                    max_score, max_csq, gene_of_max = sc, term, (gene or None)
                bit = termbit.get(term)
                if bit is not None:
                    mask |= 1 << bit

        yield (rec.ID or f"{chrom}:{rec.POS}_{ref}_{'|'.join(alts)}",
               chrom, rec.POS, ref, enc.alleles_string(sdpn, ref),
               enc.sdp_from_sdpn(sdpn, ref), sdpn, enc.calc_sdp(fis),
               gene_of_max, max_csq, mask, "snp")
        kept += 1
    stats.update(seen=seen, kept=kept, csq=ncsq)


def build(vcf_path, chroms, out_dir, genes_tsv, so_tsv, *, name, release, grc,
          gene_version, keeplow=True, keepno=False):
    from cyvcf2 import VCF  # imported here so the query/serve path needs no cyvcf2

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"fv.{release}.snps.sqlite")
    if os.path.exists(out_path):
        os.remove(out_path)

    so_rows, so2score, termbit = _load_so(so_tsv)
    genes = _load_genes(genes_tsv)
    vcf = VCF(vcf_path)
    smp = {s: i for i, s in enumerate(vcf.samples)}
    fidx = [smp.get(f) for f in enc.FOUNDERS]
    ci = _csq_indices(vcf)
    seqs = set(vcf.seqnames)
    chroms = [c for c in (CANONICAL_CHROMS if chroms == ["all"] else chroms) if c in seqs]

    con = sqlite3.connect(out_path)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    for ddl in schema.DDL_TABLES:
        con.execute(ddl)

    total_v = total_c = 0
    with step(f"build {release} ({grc}) variants"):
        for chrom in chroms:
            stats: dict[str, int] = {}
            rows = parse_chrom(vcf, chrom, fidx, ci, so2score, termbit,
                               keeplow, keepno, stats)
            while batch := list(itertools.islice(rows, BATCH)):
                con.executemany(_INSERT_VARIANT, batch)
            con.commit()
            total_v += stats["kept"]
            total_c += stats["csq"]
            LOG.info(f"chr{chrom:<3} seen {stats['seen']:,} kept {stats['kept']:,} "
                     f"csq {stats['csq']:,}")

    _load_reference_tables(con, genes, so_rows, termbit)
    _write_meta(con, genes_tsv, gene_version, vcf_path, chroms,
                name=name, release=release, grc=grc, n_genes=len(genes),
                n_variants=total_v, n_consequences=total_c)
    with step("index"):
        for ddl in schema.DDL_INDEXES:
            con.execute(ddl)
    con.commit()
    con.close()
    LOG.info("SQLite %.2fGB -> %s", os.path.getsize(out_path) / 1e9, out_path)
    return out_path


def _load_reference_tables(con, genes, so_rows, termbit):
    con.executemany("INSERT INTO genes VALUES (?,?,?,?,?,?,?,?)", [
        (g["ensembl_id"], g.get("ensembl_version", ""), g["symbol"], g.get("name", ""),
         g["chromosome"], int(g["start_position"]), int(g["end_position"]), int(g["strand"]))
        for g in genes])
    con.executemany("INSERT INTO so_terms VALUES (?,?,?,?,?,?,?)", [
        (r["term"], r["description"], r["accession"], r["display"], r["impact"],
         int(r["score"]), termbit[r["term"]]) for r in so_rows])


def _write_meta(con, genes_tsv, gene_version, vcf_path, chroms, *,
                name, release, grc, n_genes, n_variants, n_consequences):
    created = _dt.datetime.now(_dt.UTC).isoformat(timespec="seconds")
    con.executemany("INSERT INTO meta VALUES (?,?)", [
        ("schema_version", schema.SCHEMA_VERSION),
        ("release", release), ("name", name), ("variant", "SNP"),
        ("grc", grc), ("created", created),
        ("founders", ",".join(enc.FOUNDERS)),
        ("gene_annotation_file", os.path.basename(genes_tsv)),
        ("gene_annotation_version", gene_version),
        ("source_vcf", os.path.basename(vcf_path)),
        ("chroms", ",".join(chroms)),
        ("n_genes", str(n_genes)), ("n_variants", str(n_variants)),
        # Counted, not stored: the transcript-level consequence rows folded into csq_mask.
        ("n_consequences", str(n_consequences)),
        ("built_with", "foundervariants.build streaming cyvcf2->sqlite"),
        ("served_from", "single-sqlite"), ("csq_model", "per-variant bitmask"),
    ])
    con.execute("INSERT INTO description VALUES (?,?,?)",
                [f"FounderVariants {name}", grc, created])
