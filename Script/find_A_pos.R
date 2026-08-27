#!/usr/bin/env Rscript
## Find the A-site offset per fragment length from the start-codon pileup.
##
## COORDINATE CONVENTION (this is where it is easy to go wrong)
## -----------------------------------------------------------
## compute_profile_all.pl centres the metagene on field 5 of the GTF
## start_codon record, i.e. on the LAST nucleotide of the AUG. So metagene
## position q corresponds to CDS offset q + 2 for a 5'-anchored read, and to
## CDS offset (q + 2 - L) for a 3'-anchored one.
##
## Initiating ribosomes hold the AUG in the P site, so their A-site codon
## begins at CDS offset 3. CountingFullSeq_Apos.pl needs an offset A with
##     (read CDS offset) + A = 3            [5p]
##     (read CDS offset) + L - A = 3        [3p]
## Both reduce to
##     A = |q - 1|
## so the reported offset is the modal 5'/3'-end position shifted by ONE
## nucleotide, not the modal position itself.
##
## Ragged ends put reads at q in {peak-1, peak, peak+1}; subtracting one gives
## the three reported offsets {peak-2, peak-1, peak}. They are consecutive by
## construction, one per codon frame, and all three derive from a single peak
## estimate - so a noisy minor frame cannot pull one of them onto a different
## tract. CountingFullSeq_Apos.pl looks them up by residue mod 3, which for a
## read at q equals (q - 1) mod 3, i.e. the residue of its own offset member.
##
## Usage:
##   Rscript find_A_pos.R <start_pos.tsv> <L_1> <L_2> <5p|3p> <out.tsv> <out.pdf>
##
## Output TSV (no header):  <length>\t<A|res 0>\t<A|res 1>\t<A|res 2>

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
## Search window for the initiation peak. A ~30 nt footprint with the A site
## ~15 nt from the 5' end puts the 5'-anchored peak near -14; +/- 5 nt around
## that. In 3p-anchored mode the mirrored logic gives a window around +15.
## --------------------------------------------------------------------------

window <- if (A_site_end == "5p") c(-20L, -10L) else c(10L, 20L)

## --------------------------------------------------------------------------
## one peak per length -> three consecutive offsets, one per frame
## --------------------------------------------------------------------------

consensus_per_length <- lapply(l, function(L) {
  dens <- sum_pos.l[, L]
  names(dens) <- rownames(sum_pos.l)
  peak <- as.integer(names(which.max(dens[as.character(window[1]:window[2])])))
  if (peak %in% window) {
    message(sprintf(paste0("WARNING: L = %s: peak sits on the edge of the search ",
                           "window (%+d in [%+d, %+d]); the true peak may lie outside"),
                    L, peak, window[1], window[2]))
  }
  members <- c(peak - 2L, peak - 1L, peak)     # modal positions minus one
  list(peak = peak,
       off  = setNames(members, as.character(members %% 3L)))
})
names(consensus_per_length) <- l

## --------------------------------------------------------------------------
## per-length diagnostic PDF: bars coloured by codon frame, the peak circled
## and the three reported offsets starred
## --------------------------------------------------------------------------

pdf(out_file_pdf)
par(mfrow = c(2, 2), pty = "s")
for (L in l) {
  pk <- consensus_per_length[[L]]
  plot(1:nrow(sum_pos.l), sum_pos.l[, L],
       xlim = c(70, 130), main = paste0("L = ", L, " nt"),
       xaxt = "n", lty = "blank", cex = 0,
       xlab = "Position (0 = last nt of AUG)", ylab = "Mean footprint density")
  abline(v = 101, col = "grey")
  axis(at = 1:nrow(sum_pos.l), labels = rownames(sum_pos.l), side = 1,
       cex.axis = 0.4)
  for (k in 1:3) {
    lines(seq(k, nrow(sum_pos.l), 3), sum_pos.l[seq(k, nrow(sum_pos.l), 3), L],
          col = c("darkred", "darkblue", "darkgreen")[k], type = "h")
  }
  points(pk$peak + 101L, sum_pos.l[as.character(pk$peak), L],
         pch = 1, cex = 1.8, lwd = 1.2, col = "black")
  for (m in pk$off) {
    points(m + 101L, sum_pos.l[as.character(m), L],
           pch = 8, cex = 1.0, lwd = 1.2, col = "black")
    text(m + 101L, sum_pos.l[as.character(m), L],
         labels = sprintf("A=%d", abs(m)), pos = 3, cex = 0.55)
  }
}
dev.off()

## --------------------------------------------------------------------------
## write inferred.tsv, columns ordered by residue mod 3 (the key used by
## CountingFullSeq_Apos.pl)
## --------------------------------------------------------------------------

A_pos <- t(sapply(consensus_per_length, function(x)
  abs(c(x$off[["0"]], x$off[["1"]], x$off[["2"]]))))
rownames(A_pos) <- l
write.table(A_pos, file = out_file_tsv, sep = "\t", quote = FALSE,
            col.names = FALSE)
