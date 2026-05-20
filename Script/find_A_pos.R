#!/usr/bin/env Rscript
## Find the A-site position per fragment length from tri-repeat pileup data.
##
## Robust to per-frame noise: instead of taking an independent argmax in each
## codon frame (which can pick a different tract when one frame is noisy at
## that length), we find the position p in the search window that maximises
##     density(p-1) + density(p) + density(p+1)
## and assign the three frame peaks to (p-1, p, p+1) - guaranteed consecutive,
## one peak per frame.
##
## Usage (compatible with the existing snakemake rule):
##   Rscript find_A_pos.R <start_pos.tsv> <L_1> <L_2> <5p|3p> <out.tsv> <out.pdf>
##
## Output TSV format (no header):
##   <length>\t<|peak_F2|>\t<|peak_F3|>\t<|peak_F1|>
## - same column order as the legacy script (drop-in for DT fitting)

library(data.table)

args         <- commandArgs(trailingOnly = TRUE)
start_pos    <- args[1]
L_1          <- as.integer(args[2])
L_2          <- as.integer(args[3])
A_site_end   <- args[4]
out_file_tsv <- args[5]
out_file_pdf <- args[6]

## --------------------------------------------------------------------------
## read + row-normalise pileup, then average per length
## --------------------------------------------------------------------------

table_pos  <- fread(start_pos, sep = "\t", stringsAsFactors = FALSE)[, -203]
length_pos <- as.numeric(gsub("L:", "", table_pos$V1))
table_pos  <- table_pos[, -1]
colnames(table_pos) <- as.character(c(-100:0, 1:100))

tp         <- rowSums(table_pos) > 10
table_pos  <- subset(table_pos, tp)
length_pos <- length_pos[tp]

sum_pos    <- sweep(table_pos, MARGIN = 1, rowSums(table_pos), FUN = "/")
ss.pos     <- split(1:nrow(sum_pos), length_pos)
sum_pos.l  <- sapply(ss.pos, function(x) colMeans(sum_pos[x, ]))
rownames(sum_pos.l) <- as.character(c(-100:0, 1:100))

l <- as.character(L_1:L_2)
l <- l[l %in% colnames(sum_pos.l)]

## --------------------------------------------------------------------------
## search window for the triplet (tighter than the legacy pos_am-derived
## window). Biological rationale: typical ribosome footprint ~30 nt with
## the A-site ~15 nt from the 5' end, so the 5'-anchored peak sits at
## ~-15 nt from the codon centre. We allow +/- 5 nt around that.
## In 3p-anchored mode the same logic gives a window around +15.
## --------------------------------------------------------------------------

window <- if (A_site_end == "5p") c(-20L, -10L) else c(10L, 20L)

## triplet centres p s.t. (p-1, p, p+1) all lie inside the window.
## As p slides by 1, the mod-3 frame assignment of (p-1, p, p+1) cycles
## through (F1,F2,F3) -> (F2,F3,F1) -> (F3,F1,F2), so all three cyclic
## configurations are tested by this loop.
search_p <- (window[1] + 1L):(window[2] - 1L)

frame_of_pos <- function(p) {
  idx <- p + 101L
  c("F1", "F2", "F3")[((idx - 1L) %% 3L) + 1L]
}

## --------------------------------------------------------------------------
## consensus peak per length
## --------------------------------------------------------------------------

consensus_per_length <- lapply(l, function(L) {
  dens <- sum_pos.l[, L]
  names(dens) <- rownames(sum_pos.l)
  triplet_sums <- vapply(search_p, function(p) {
    sum(dens[as.character(c(p - 1L, p, p + 1L))], na.rm = TRUE)
  }, numeric(1))
  p_best  <- search_p[which.max(triplet_sums)]
  triplet <- c(p_best - 1L, p_best, p_best + 1L)
  frames  <- frame_of_pos(triplet)
  list(F1 = triplet[frames == "F1"],
       F2 = triplet[frames == "F2"],
       F3 = triplet[frames == "F3"])
})
names(consensus_per_length) <- l

## --------------------------------------------------------------------------
## per-length diagnostic PDF (matches legacy layout, adds star markers
## on the chosen consensus peaks)
## --------------------------------------------------------------------------

pdf(out_file_pdf)
par(mfrow = c(2, 2), pty = "s")
for (L in l) {
  pk <- consensus_per_length[[L]]
  plot(1:nrow(sum_pos.l), sum_pos.l[, L],
       xlim = c(70, 130), main = paste0("L = ", L, " nt"),
       xaxt = "n", lty = "blank", cex = 0,
       xlab = "Position", ylab = "Mean footprint density")
  abline(v = 101, col = "grey")
  axis(at = 1:nrow(sum_pos.l), labels = rownames(sum_pos.l), side = 1,
       cex.axis = 0.4)
  lines(seq(1, nrow(sum_pos.l), 3),
        sum_pos.l[seq(1, nrow(sum_pos.l), 3), L],
        col = "darkred",   type = "h")
  lines(seq(2, nrow(sum_pos.l), 3),
        sum_pos.l[seq(2, nrow(sum_pos.l), 3), L],
        col = "darkblue",  type = "h")
  lines(seq(3, nrow(sum_pos.l), 3),
        sum_pos.l[seq(3, nrow(sum_pos.l), 3), L],
        col = "darkgreen", type = "h")
  for (fr in c("F1", "F2", "F3")) {
    p <- pk[[fr]]
    if (length(p) == 1L) {
      points(p + 101L, sum_pos.l[as.character(p), L],
             pch = 8, cex = 1.0, lwd = 1.2, col = "black")
      text(p + 101L, sum_pos.l[as.character(p), L],
           labels = sprintf("%+d", p), pos = 3, cex = 0.55)
    }
  }
}
dev.off()

## --------------------------------------------------------------------------
## write inferred.tsv (legacy column order: |F2|, |F3|, |F1|)
## --------------------------------------------------------------------------

A_pos <- t(sapply(consensus_per_length, function(x) {
  c(abs(x$F2), abs(x$F3), abs(x$F1))
}))
rownames(A_pos) <- l
write.table(A_pos, file = out_file_tsv, sep = "\t", quote = FALSE,
            col.names = FALSE)
