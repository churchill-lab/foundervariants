#!/bin/bash
# Compare freshly rebuilt databases in dbs_verify/ against the production ones in dbs/,
# at full scale, before anything is swapped. Read-only on both sides.
#
#   ./scripts/compare_builds.sh            # both releases
#   ./scripts/compare_builds.sh 2021       # one release
#
# Method: an indexed merge join on (chr, pos, snp_id) comparing all 12 columns, plus a
# key-coverage check proving the join is 1:1 with no unmatched rows on either side. That
# pair of checks is equivalent to a two-directional EXCEPT but runs in minutes instead of
# hours — a whole-table EXCEPT has to sort 36M twelve-column tuples and took >1h before
# being abandoned. If "matched" equals both row counts AND "differing" is 0, the two
# variants tables are identical as sets of rows.
#
# Physical row order is NOT compared, and does not matter: every query and R/qtl2 itself
# goes through the chr_pos index, and nothing does ORDER BY rowid.
set -uo pipefail
cd "$(dirname "$0")/.."

RELEASES=${@:-"2021 1505"}
FAIL=0

for rel in $RELEASES; do
  OLD="dbs/fv.${rel}.snps.sqlite"
  NEW="dbs_verify/fv.${rel}.snps.sqlite"
  echo "═══════════════════════════════════════════════════════════════"
  echo "  RELEASE $rel     production: $OLD"
  echo "                   rebuilt   : $NEW"
  echo "═══════════════════════════════════════════════════════════════"
  if [[ ! -f "$NEW" ]]; then echo "  !! rebuilt file missing"; FAIL=1; continue; fi
  stat -f "  %z bytes  %N" "$OLD" "$NEW"

  # Schema compact-v2 ADDS variant_genes and the gi_id index; everything production has
  # (above all the R/qtl2 `variants` table and chr_pos) must be byte-identical. Statements
  # are compared whole (sqlite_master.sql), so a multi-line CREATE can't half-match.
  echo "  -- stored schema text (production objects must be identical) --"
  V2_NEW="'variant_genes','gi_id'"
  SCHEMA_SQL="SELECT type||' '||name||': '||replace(sql,char(10),' ') FROM sqlite_master WHERE sql IS NOT NULL"
  if diff <(sqlite3 "$OLD" "$SCHEMA_SQL ORDER BY name") \
          <(sqlite3 "$NEW" "$SCHEMA_SQL AND name NOT IN ($V2_NEW) ORDER BY name") > /dev/null; then
    echo "     IDENTICAL"
  else
    echo "     DIFFERS:"
    diff <(sqlite3 "$OLD" "$SCHEMA_SQL ORDER BY name") \
         <(sqlite3 "$NEW" "$SCHEMA_SQL AND name NOT IN ($V2_NEW) ORDER BY name") | sed 's/^/       /'
    FAIL=1
  fi
  echo "  -- added in compact-v2 --"
  sqlite3 "$NEW" "SELECT name FROM sqlite_master WHERE name IN ($V2_NEW) ORDER BY name" | sed 's/^/     /'

  echo "  -- rows: coverage, then all 12 columns --"
  sqlite3 <<SQL | sed 's/^/     /'
ATTACH '$OLD' AS o;
ATTACH '$NEW' AS n;
.mode list
SELECT 'variants  old = ' || (SELECT count(*) FROM o.variants)
     || '   new = '       || (SELECT count(*) FROM n.variants);
SELECT 'matched 1:1   = ' || count(*) FROM o.variants o
  JOIN n.variants n ON o.chr=n.chr AND o.pos=n.pos AND o.snp_id=n.snp_id;
SELECT 'DIFFERING     = ' || count(*) FROM o.variants o
  JOIN n.variants n ON o.chr=n.chr AND o.pos=n.pos AND o.snp_id=n.snp_id
 WHERE o.ref IS NOT n.ref OR o.alleles IS NOT n.alleles OR o.sdp IS NOT n.sdp
    OR o.sdpn IS NOT n.sdpn OR o.confidence IS NOT n.confidence
    OR o.ensembl_gene IS NOT n.ensembl_gene OR o.consequence IS NOT n.consequence
    OR o.csq_mask IS NOT n.csq_mask OR o.type IS NOT n.type;
SELECT 'genes     old = ' || (SELECT count(*) FROM o.genes)
     || '   new = '       || (SELECT count(*) FROM n.genes)
     || '   diffs = '     || (SELECT count(*) FROM (SELECT * FROM o.genes EXCEPT SELECT * FROM n.genes));
SELECT 'so_terms  old = ' || (SELECT count(*) FROM o.so_terms)
     || '   new = '       || (SELECT count(*) FROM n.so_terms)
     || '   diffs = '     || (SELECT count(*) FROM
        (SELECT term,bit,score FROM o.so_terms EXCEPT SELECT term,bit,score FROM n.so_terms));
SELECT 'sdp range new = ' || (SELECT min(sdp)||'..'||max(sdp) FROM n.variants)
     || '   out-of-range(1..254) = ' || (SELECT count(*) FROM n.variants WHERE sdp<1 OR sdp>254);
SQL

  echo "  -- if DIFFERING is non-zero, what changed --"
  sqlite3 <<SQL | sed 's/^/     /'
ATTACH '$OLD' AS o;
ATTACH '$NEW' AS n;
.mode column
.headers on
SELECT coalesce(o.consequence,'(NULL)') AS production_had,
       coalesce(n.consequence,'(NULL)') AS rebuilt_has, count(*) AS rows
FROM o.variants o JOIN n.variants n ON o.chr=n.chr AND o.pos=n.pos AND o.snp_id=n.snp_id
WHERE o.consequence IS NOT n.consequence OR o.ensembl_gene IS NOT n.ensembl_gene
GROUP BY 1,2 ORDER BY rows DESC LIMIT 20;
SQL

  echo "  -- meta (expect: created; compact-v2 adds n_variant_genes, schema_version, csq_model) --"
  diff <(sqlite3 "$OLD" "SELECT meta_name||' = '||meta_value FROM meta ORDER BY 1") \
       <(sqlite3 "$NEW" "SELECT meta_name||' = '||meta_value FROM meta ORDER BY 1") \
       | sed 's/^/     /' || true
  echo
done

echo "  Next: scripts/check_variant_genes.py on each rebuilt file (per-gene invariants)."
echo "═══════════════════════════════════════════════════════════════"
if [[ $FAIL -eq 0 ]]; then
  echo "  To swap: 'matched 1:1' must equal BOTH row counts, and every"
  echo "  DIFFERING row must be one you have explained and accepted."
else
  echo "  SCHEMA MISMATCH — do not swap."
fi
echo "═══════════════════════════════════════════════════════════════"
