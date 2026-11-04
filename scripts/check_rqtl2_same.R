#!/usr/bin/env Rscript
# R/qtl2 must see exactly the same data in a rebuilt file as in production.
#
#   Rscript scripts/check_rqtl2_same.R dbs/fv.2021.snps.sqlite dbs_verify/fv.2021.snps.sqlite
#
# Runs qtl2::create_variant_query_func() with NO field overrides on both files and compares
# the returned data frames, across windows on several chromosomes (whole chromosomes where
# that's cheap). Then feeds the rebuilt file's rows to index_snps(). Exits non-zero on any
# difference. Read-only on both files.
suppressMessages(library(qtl2))
a <- commandArgs(TRUE)
if (length(a) != 2) stop("usage: check_rqtl2_same.R OLD.sqlite NEW.sqlite")
qo <- create_variant_query_func(a[1])
qn <- create_variant_query_func(a[2])

windows <- list(c("1", 0, 200), c("2", 50, 60), c("11", 100, 102), c("19", 0, 70),
                c("X", 0, 200))
bad <- 0
for (w in windows) {
  chr <- w[1]; lo <- as.numeric(w[2]); hi <- as.numeric(w[3])
  o <- qo(chr, lo, hi); n <- qn(chr, lo, hi)
  o <- o[order(o$pos, o$snp), , drop = FALSE]; n <- n[order(n$pos, n$snp), , drop = FALSE]
  rownames(o) <- NULL; rownames(n) <- NULL
  same <- identical(names(o), names(n)) && isTRUE(all.equal(o, n, check.attributes = FALSE))
  if (!same) bad <- bad + 1
  cat(sprintf("chr%-2s %6.1f-%6.1f Mbp  rows old=%9d new=%9d  identical=%s\n",
              chr, lo, hi, nrow(o), nrow(n), same))
}
if (bad) { cat(bad, "window(s) differ\n"); quit(status = 1) }
n <- qn("11", 100, 102)
map <- list("11" = setNames(seq(100, 102, by = 0.5), paste0("m", 1:5)))
s <- index_snps(map, n)
cat(sprintf("index_snps(chr11:100-102): %d SNPs -> %d index SNPs; sdp %d-%d; columns: %s\n",
            nrow(s), length(unique(s$index)), min(n$sdp), max(n$sdp),
            paste(names(n), collapse = ",")))
cat("R/qtl2 sees identical data in both files\n")
