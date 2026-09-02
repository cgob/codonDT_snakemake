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
##   Rscript find_A_pos.R <start_pos.tsv> <L_1> <L_2> <5p|3p> <out.tsv> <out.pdf> \
##                        [<win_lo> <win_hi>]
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
win_args     <- if (length(args) >= 8) args[7:8] else NULL

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

## Taken from config.yaml (A_site_window) when supplied; the old hard-coded
## monosome defaults are kept as a fallback. The window is in METAGENE
## coordinates (0 = last nt of the AUG), so its sign has to match A_site_end:
## a 5'-anchored read sits upstream of the AUG (negative), a 3'-anchored one
## downstream (positive). Getting that pairing wrong used to fail silently, so
## it is now an error.
if (!is.null(win_args)) {
  window <- as.integer(win_args)
  if (any(is.na(window)))
    stop("A_site_window must be two integers, got: ", paste(win_args, collapse = ", "))
} else {
  window <- if (A_site_end == "5p") c(-20L, -10L) else c(10L, 20L)
}
if (window[1] >= window[2])
  stop(sprintf("A_site_window must be increasing, got [%d, %d]", window[1], window[2]))
if (A_site_end == "5p" && window[2] > 0)
  stop(sprintf(paste0("A_site_end is 5p (read 5' end, upstream of the AUG) but ",
                      "A_site_window = [%d, %d] is not negative. See the ",
                      "A_site_window comment in config.yaml."), window[1], window[2]))
if (A_site_end == "3p" && window[1] < 0)
  stop(sprintf(paste0("A_site_end is 3p (read 3' end, downstream of the AUG) but ",
                      "A_site_window = [%d, %d] is not positive. See the ",
                      "A_site_window comment in config.yaml."), window[1], window[2]))
message(sprintf("A_site_end = %s, search window = [%+d, %+d]",
                A_site_end, window[1], window[2]))

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

## Plotted in metagene coordinates (not row indices), over a range wide enough
## to contain the peak wherever it actually is. The old fixed xlim = c(70, 130)
## was in index units, i.e. -31..+29, which cut off the initiation peak of any
## fragment long enough to reach further upstream - a 64 nt disome peaks near
## -42 and was simply off the page, so a window sitting on noise looked fine.
## The search window is now shaded and the global maximum of the plotted range
## dashed, so a peak outside the window is visible at a glance.
pdf(out_file_pdf)
par(mfrow = c(2, 2), pty = "s")
pp <- as.integer(rownames(sum_pos.l))
for (L in l) {
  pk <- consensus_per_length[[L]]
  Ln <- as.integer(L)
  if (A_site_end == "5p") {
    xlo <- min(window[1] - 15L, -(Ln + 15L)); xhi <- max(window[2] + 15L, 40L)
  } else {
    xlo <- min(window[1] - 15L, -40L);        xhi <- max(window[2] + 15L, Ln + 15L)
  }
  xlo <- max(xlo, min(pp)); xhi <- min(xhi, max(pp))

  plot(pp, sum_pos.l[, L],
       xlim = c(xlo, xhi), main = paste0("L = ", L, " nt"),
       xaxt = "n", lty = "blank", cex = 0,
       xlab = "Position (0 = last nt of AUG)", ylab = "Mean footprint density")
  usr <- par("usr")
  rect(window[1], usr[3], window[2], usr[4], col = rgb(1, 0, 0, 0.10), border = NA)
  abline(v = 0, col = "grey")
  axis(side = 1, at = seq(-100L, 100L, by = 5L), cex.axis = 0.5)
  for (k in 1:3) {
    idx <- seq(k, length(pp), 3)
    lines(pp[idx], sum_pos.l[idx, L],
          col = c("darkred", "darkblue", "darkgreen")[k], type = "h")
  }
  sel <- pp >= xlo & pp <= xhi
  gi  <- which(sel)[which.max(sum_pos.l[sel, L])]
  abline(v = pp[gi], lty = 2, col = "grey40")
  if (pp[gi] < window[1] || pp[gi] > window[2])
    mtext(sprintf("global max %+d OUTSIDE window", pp[gi]),
          side = 3, line = -1, cex = 0.55, col = "red")
  points(pk$peak, sum_pos.l[as.character(pk$peak), L],
         pch = 1, cex = 1.8, lwd = 1.2, col = "black")
  for (m in pk$off) {
    points(m, sum_pos.l[as.character(m), L],
           pch = 8, cex = 1.0, lwd = 1.2, col = "black")
    text(m, sum_pos.l[as.character(m), L],
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
