#!/bin/bash
# Full rebuild of both releases into dbs_verify/, for comparison against the production
# databases in dbs/ before anything is swapped. Production is never written to.
set -euo pipefail
cd "$(dirname "$0")/.."

OUT=dbs_verify
SANGER="$HOME/work/RESOURCE/Sanger"
mkdir -p "$OUT"

echo "$(date '+%H:%M:%S') START full rebuild -> $OUT (production dbs/ untouched)"

echo "$(date '+%H:%M:%S') building 2021 (GRCm39)..."
uv run --extra build foundervariants -v build \
  --vcf "$SANGER/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.rsID.vcf.gz" \
  --genes-tsv ../data/ensimpl.104.mm.genes.tsv \
  --so-tsv    ../data/so_terms.tsv \
  --release 2021 --grc GRCm39 --gene-version "Ensembl 104 (ensimpl.104.Mm)" \
  --out "$OUT"

echo "$(date '+%H:%M:%S') building 1505 (GRCm38)..."
uv run --extra build foundervariants -v build \
  --vcf "$SANGER/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz" \
  --genes-tsv ../data/ensimpl.75.mm.genes.tsv \
  --so-tsv    ../data/so_terms.tsv \
  --release 1505 --grc GRCm38 --gene-version "Ensembl 75 (ensimpl.75.Mm)" \
  --out "$OUT"

echo "$(date '+%H:%M:%S') DONE"
ls -la "$OUT"
