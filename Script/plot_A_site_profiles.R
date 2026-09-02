#!/usr/bin/env Rscript
## Diagnostic: start-codon metagene profile for EVERY fragment length, plotted
## over a wide window so a peak anywhere in the pile-up is visible.
##
## Unlike find_A_pos.R this does no peak selection and applies no search
## window, and it plots every length present rather than only L1:L2. That is
## what is needed to CHOOSE the size range and the search window: a population
## outside the current thresholds is invisible in the find_A_pos.R output.
##
## compute_profile_all.pl records the read end named by A_site_end: the 5' end
## ($b[3]) for 5p, the 3' end ($b[3] + $length) for 3p. So the pile-up is in
## whichever anchored coordinate the run used, and the OTHER end is at
## peak +/- L. Both are shown, because the one that comes out flat against
## length identifies which end the population is really anchored on.
##
## Usage:
##   Rscript plot_A_site_profiles.R <A_site_pos.tsv> <out.pdf> \
##                                  [xlo xhi] [win_lo win_hi] [5p|3p]

suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
inf  <- args[1]; outf <- args[2]
xlo  <- if (length(args) >= 3) as.integer(args[3]) else -100L
xhi  <- if (length(args) >= 4) as.integer(args[4]) else   60L
win  <- if (length(args) >= 6) as.integer(args[5:6]) else NULL
A_site_end <- if (length(args) >= 7) args[7] else "5p"

anchor <- if (A_site_end == "3p") "3'" else "5'"
other  <- if (A_site_end == "3p") "5'" else "3'"
sgn    <- if (A_site_end == "3p") -1L else 1L

tp <- fread(inf, sep = "\t", stringsAsFactors = FALSE)[, -203]
lp <- as.numeric(gsub("L:", "", tp$V1)); tp <- tp[, -1]
pos <- c(-100:0, 1:100); colnames(tp) <- as.character(pos)
nreads <- rowSums(tp)
keep <- nreads > 10; tp <- subset(tp, keep); lp <- lp[keep]; nreads <- nreads[keep]

sp <- sweep(tp, 1, rowSums(tp), "/")
ss <- split(seq_len(nrow(sp)), lp)
m  <- sapply(ss, function(x) colMeans(sp[x, ]))
Ls <- as.integer(colnames(m))
ntx <- sapply(ss, length)
rds <- sapply(ss, function(x) sum(nreads[x]))

## peak at the anchored end, and the position that implies for the other end
pkA <- apply(m, 2, function(d) pos[which.max(d)])
pkO <- pkA + sgn * Ls
if (!is.null(win)) {
  wsel <- pos >= win[1] & pos <= win[2]
  pwA <- apply(m, 2, function(d) pos[wsel][which.max(d[wsel])])
  pwO <- pwA + sgn * Ls
} else { pwA <- NULL; pwO <- NULL }

## consensus of the in-window peaks (falling back to the global ones): a
## length deviating from it is flagged, without assuming any particular value
ref  <- if (!is.null(pwA)) pwA else pkA
refO <- if (!is.null(pwO)) pwO else pkO
medA <- median(ref); medO <- median(refO)

pdf(outf, width = 11, height = 8)

## ---- page 1: summary ------------------------------------------------------
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

plot(Ls, rds, type = "h", lwd = 4, col = "grey35",
     xlab = "Fragment length (nt)", ylab = "Reads in metagene",
     main = "Reads per length")

plot(Ls, ref, pch = 19, col = "darkblue", xlab = "Fragment length (nt)",
     ylab = sprintf("%s-end peak position", anchor),
     main = sprintf("%s end (the anchored coordinate)", anchor))
abline(h = medA, col = "darkgreen", lty = 2)
legend("bottomleft", bty = "n", cex = 0.7,
       legend = sprintf("median %+d; flat here = %s-end anchored", medA, anchor))

plot(Ls, refO, pch = 19, col = "darkred", xlab = "Fragment length (nt)",
     ylab = sprintf("implied %s-end position", other),
     main = sprintf("%s end (peak %s L)", other, if (sgn > 0) "+" else "-"))
abline(h = medO, col = "darkgreen", lty = 2)
legend("topleft", bty = "n", cex = 0.7,
       legend = sprintf("median %+d; flat here = %s-end anchored", medO, other))

plot(Ls, ntx, type = "h", lwd = 4, col = "grey55",
     xlab = "Fragment length (nt)", ylab = "Transcripts contributing",
     main = "Support per length (n transcripts)")

## ---- per-length panels ----------------------------------------------------
par(mfrow = c(3, 2), mar = c(4, 4, 2.5, 1))
for (j in seq_along(Ls)) {
  L <- Ls[j]; d <- m[, j]
  plot(pos, d, type = "n", xlim = c(xlo, xhi),
       xlab = sprintf("Read %s end, position (0 = last nt of AUG)", anchor),
       ylab = "Mean density",
       main = sprintf("L = %d nt   (%d reads, %d transcripts)", L, rds[j], ntx[j]))
  usr <- par("usr")
  if (!is.null(win)) rect(win[1], usr[3], win[2], usr[4], col = rgb(1, 0, 0, 0.09), border = NA)
  abline(v = 0, col = "grey")
  for (k in 1:3) {
    idx <- seq(k, length(pos), 3)
    lines(pos[idx], d[idx], col = c("darkred", "darkblue", "darkgreen")[k], type = "h")
  }
  abline(v = pkA[j], lty = 2, lwd = 2)
  mtext(sprintf("global max %+d  ->  %s end %+d", pkA[j], other, pkO[j]),
        side = 3, line = -1.2, cex = 0.6,
        col = if (abs(pkA[j] - medA) <= 2) "darkgreen" else "red")
  if (!is.null(pwA)) {
    points(pwA[j], m[as.character(pwA[j]), j], pch = 1, cex = 2, lwd = 2, col = "blue")
    mtext(sprintf("in-window %+d  ->  %s end %+d   (median %+d)", pwA[j], other, pwO[j], medA),
          side = 3, line = -2.1, cex = 0.6,
          col = if (abs(pwA[j] - medA) <= 2) "darkgreen" else "red")
  }
  axis(side = 1, at = seq(-100, 100, by = 10), cex.axis = 0.6)
}
invisible(dev.off())
cat(sprintf("wrote %s  (%d lengths, %s-anchored, xlim %d..%d)\n",
            outf, length(Ls), anchor, xlo, xhi))
