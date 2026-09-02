#!/usr/bin/env Rscript
## Diagnostic: start-codon metagene profile for EVERY fragment length, plotted
## over a wide window so a peak anywhere in the pileup is visible.
##
## Unlike find_A_pos.R this does no peak selection and applies no search
## window - it just shows where the density actually is, so the window can be
## chosen from the data instead of assumed.
##
## Usage:
##   Rscript plot_A_site_profiles.R <A_site_pos.tsv> <out.pdf> [xlo] [xhi] [win_lo win_hi]
##
## If a window is given, the peak WITHIN it is circled alongside the global
## maximum. That matters for low-count lengths, where a lone spurious spike can
## win the unconstrained argmax while the real peak sits plainly in the window.

suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
inf  <- args[1]; outf <- args[2]
xlo  <- if (length(args) >= 3) as.integer(args[3]) else -100L
xhi  <- if (length(args) >= 4) as.integer(args[4]) else  60L
win  <- if (length(args) >= 6) as.integer(args[5:6]) else NULL

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

## per-length global peak, and where that puts the 3' end
pk5 <- apply(m, 2, function(d) pos[which.max(d)])
pk3 <- pk5 + Ls
## peak restricted to the window, if one was supplied
if (!is.null(win)) {
  wsel <- pos >= win[1] & pos <= win[2]
  pw5 <- apply(m, 2, function(d) pos[wsel][which.max(d[wsel])])
  pw3 <- pw5 + Ls
} else { pw5 <- NULL; pw3 <- NULL }

pdf(outf, width = 11, height = 8)

## ---- page 1: summary ------------------------------------------------------
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

plot(Ls, rds, type = "h", lwd = 4, col = "grey35",
     xlab = "Fragment length (nt)", ylab = "Reads in metagene",
     main = "Reads per length")

plot(Ls, pk5, pch = 19, col = "darkblue", xlab = "Fragment length (nt)",
     ylab = "5' peak position", main = "5'-anchored peak vs length",
     ylim = range(c(pk5, -50, 0)))
abline(h = -15, col = "red", lty = 2)
lines(Ls, 22 - Ls, col = "darkgreen", lty = 3)
legend("bottomleft", bty = "n", cex = 0.7,
       legend = c("flat at -15 = 5'-anchored", "diagonal = 3'-anchored at +22"),
       col = c("red", "darkgreen"), lty = c(2, 3))

plot(Ls, pk3, pch = 19, col = "darkred", xlab = "Fragment length (nt)",
     ylab = "3' end of the peak (peak + L)", main = "3'-anchored peak vs length",
     ylim = range(c(pk3, 0, 40)))
abline(h = 22, col = "darkgreen", lty = 2)
legend("topleft", bty = "n", cex = 0.7,
       legend = "flat at +22 = 3'-anchored", col = "darkgreen", lty = 2)

plot(Ls, ntx, type = "h", lwd = 4, col = "grey55",
     xlab = "Fragment length (nt)", ylab = "Transcripts contributing",
     main = "Support per length (n transcripts)")

## ---- per-length panels ----------------------------------------------------
par(mfrow = c(3, 2), mar = c(4, 4, 2.5, 1))
for (j in seq_along(Ls)) {
  L <- Ls[j]; d <- m[, j]
  plot(pos, d, type = "n", xlim = c(xlo, xhi),
       xlab = "Position (0 = last nt of AUG)", ylab = "Mean density",
       main = sprintf("L = %d nt   (%d reads, %d transcripts)", L, rds[j], ntx[j]))
  ## candidate windows, for reference only
  usr <- par("usr")
  rect(-20, usr[3], -10, usr[4], col = rgb(1, 0, 0, 0.08), border = NA)
  rect(-50, usr[3], -25, usr[4], col = rgb(0, 0, 1, 0.07), border = NA)
  abline(v = 0, col = "grey")
  for (k in 1:3) {
    idx <- seq(k, length(pos), 3)
    lines(pos[idx], d[idx], col = c("darkred", "darkblue", "darkgreen")[k], type = "h")
  }
  abline(v = pk5[j], lty = 2, lwd = 2)
  mtext(sprintf("global max %+d  ->  3' end %+d", pk5[j], pk3[j]),
        side = 3, line = -1.2, cex = 0.6,
        col = if (pk3[j] >= 18 && pk3[j] <= 26) "darkgreen" else "red")
  if (!is.null(pw5)) {
    points(pw5[j], m[as.character(pw5[j]), j], pch = 1, cex = 2, lwd = 2, col = "blue")
    mtext(sprintf("in-window %+d  ->  3' end %+d", pw5[j], pw3[j]),
          side = 3, line = -2.1, cex = 0.6,
          col = if (pw3[j] >= 18 && pw3[j] <= 26) "darkgreen" else "red")
  }
  axis(side = 1, at = seq(-100, 100, by = 10), cex.axis = 0.6)
}
invisible(dev.off())
cat(sprintf("wrote %s  (%d lengths, xlim %d..%d)\n", outf, length(Ls), xlo, xhi))
