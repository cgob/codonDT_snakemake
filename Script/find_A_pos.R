#!/usr/bin/env Rscript
## A-site offset per fragment length from the start-codon pileup.
## Metagene 0 is the last nt of the AUG; for a modal read-end position q the
## offset is |q - 1|, reported for the three frames {peak-2, peak-1, peak}.
## Usage: find_A_pos.R <start_pos.tsv> <L1> <L2> <5p|3p> <out.tsv> <out.pdf>
##        [<win_lo> <win_hi>] [<offset_shift>] [<fixed_offset>]
## Output (no header): <length> <A res 0> <A res 1> <A res 2>

library(data.table)

args         <- commandArgs(trailingOnly = TRUE)
start_pos    <- args[1]
L_1          <- as.integer(args[2])
L_2          <- as.integer(args[3])
A_site_end   <- args[4]
out_file_tsv <- args[5]
out_file_pdf <- args[6]
win_args     <- if (length(args) >= 8) args[7:8] else NULL
off_shift    <- if (length(args) >= 9) as.integer(args[9]) else 0L
## fixed offset: uniform at every length, start-codon peak ignored
fixed_off    <- if (length(args) >= 10) as.integer(args[10]) else NA_integer_
if (length(args) >= 10 && is.na(fixed_off))
  stop("fixed_offset must be an integer, got: ", args[10])
if (!is.na(fixed_off) && off_shift != 0L)
  stop("fixed_offset and offset_shift are mutually exclusive: the forced offset ",
       "is already absolute, shifting it as well is almost certainly a mistake")
if (is.na(off_shift))
  stop("offset_shift must be an integer, got: ", args[9])
if (off_shift %% 3L != 0L)
  stop(sprintf(paste0("offset_shift must be a whole number of codons (a multiple ",
                      "of 3), got %+d"),
               off_shift))

## pileup, row-normalised and averaged per length

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

## search window for the initiation peak (metagene coordinates; sign must match A_site_end)
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
if (is.na(fixed_off)) {
  message(sprintf("A_site_end = %s, search window = [%+d, %+d], offset shift = %+d nt",
                  A_site_end, window[1], window[2], off_shift))
} else {
  message(sprintf(paste0("A_site_end = %s, FORCED uniform offset = %d nt at every ",
                         "length (start-codon peak ignored, search window unused)"),
                  A_site_end, fixed_off))
}

## one peak per length -> three offsets, one per frame

consensus_per_length <- lapply(l, function(L) {
  dens <- sum_pos.l[, L]
  names(dens) <- rownames(sum_pos.l)
  peak <- if (!is.na(fixed_off)) fixed_off + 1L else
            as.integer(names(which.max(dens[as.character(window[1]:window[2])])))
  if (is.na(fixed_off) && peak %in% window) {
    message(sprintf(paste0("WARNING: L = %s: peak sits on the edge of the search ",
                           "window (%+d in [%+d, %+d]); the true peak may lie outside"),
                    L, peak, window[1], window[2]))
  }
  members <- c(peak - 2L, peak - 1L, peak) + off_shift   # modal positions minus one
  list(peak = peak,
       off  = setNames(members, as.character(members %% 3L)))
})
names(consensus_per_length) <- l

## diagnostic PDF: pileup per length coloured by frame, search window shaded
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

## inferred offsets, columns ordered by residue mod 3

A_pos <- t(sapply(consensus_per_length, function(x)
  abs(c(x$off[["0"]], x$off[["1"]], x$off[["2"]]))))
rownames(A_pos) <- l
write.table(A_pos, file = out_file_tsv, sep = "\t", quote = FALSE,
            col.names = FALSE)
