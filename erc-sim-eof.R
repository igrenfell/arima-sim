library(MASS)
library(stats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(lmtest)
library(stringr)
library(lubridate)
library(ROCR)
library(pROC)
library(maps)

###--- Configuration ---###
yearseq      <- 2006:2020
ndays.all    <- 365 * length(yearseq)
fod_dir      <- "I:\\workspace\\fod-2025" ##This is in the input folder for you
input_dir    <- "I:\\workspace\\fsim-climate\\station-data-complete\\PsuedoStationsPlusERC\\BASE_2020"
output_dir   <- "I:\\workspace\\ercm\\output-sim\\eof-output"

###--- Load GridMet NFDRS data ---###
setwd(input_dir)
flist     <- Sys.glob("*.csv")
flist     <- flist[grep("NFDRS_output", flist)]
nstations <- length(flist)

ercmat.gm   <- matrix(NA, nrow = ndays.all, ncol = nstations)
namevec     <- rep(NA, nstations)
dayseq.all  <- rep(1:365, length(yearseq))
yearseq.all <- rep(yearseq, each = 365)

for (curstation in 1:nstations) {
  curdat <- read.csv(flist[curstation])
  datev  <- as.Date(curdat$Date)
  yearn  <- as.numeric(format(datev, "%Y"))
  curdat <- curdat[yearn %in% yearseq, ]
  doy    <- strftime(curdat$Date, format = "%j")
  curdat.sub <- curdat[doy != "366", ]
  ercmat.gm[, curstation] <- curdat.sub$ERC
  namevec[curstation]     <- gsub(".csv", "", flist[curstation])
  print(curstation / nstations)
}

###--- Pyrome-station mapping ---###
subnamevec         <- gsub("NFDRS_output_", "", namevec)
pyrome.station.mat <- read.csv("Pyrome-station-selection - Sheet1.csv")
py.sub.name.mat    <- gsub("[^A-Za-z0-9]", "", pyrome.station.mat$Pyrome)
py.sub.name.mat2   <- gsub("NFDRSoutput", "", py.sub.name.mat)
station.num.vec    <- pyrome.station.mat[, 3]

py.test      <- 1:128
py.test.cols <- rep(NA, length(py.test))
stat.test    <- rep(NA, length(py.test))

for (i in 1:length(py.test)) {
  idx             <- which(py.sub.name.mat2 == subnamevec[i])
  py.test.cols[i] <- idx
  stat.test[i]    <- station.num.vec[idx]
}

##For testing, can ignore
#stat.test[38] <- "353228"
#stat.test     <- as.character(stat.test)

stat.pad     <- str_pad(stat.test, 6, pad = "0")
file.testvec <- paste0("NFDRS_", stat.pad, ".csv")

###--- Initialise per-pyrome containers ---###
summary.gm.list     <- vector("list", length(py.test))
sd.gm.vec           <- rep(0, length(py.test))
logistic.summary.gm <- vector("list", length(py.test))
dailymean.gm        <- matrix(NA, nrow = 365, ncol = length(py.test))
dailysd.gm          <- matrix(NA, nrow = 365, ncol = length(py.test))
auc.gm              <- rep(NA, length(py.test))

lft.thresh <- read.csv("large-fire-threshold-lorenz.csv")
thresh.vec <- lft.thresh$LFT
setwd("I:\\workspace\\fod-2025")
fod.mat <- read.csv("FOD_Pyromes_merge.csv")

###--- Per-pyrome: logistic models and daily climatology ---###
pdf(file.path(output_dir, "roc-curves-eof.pdf"), width = 7, height = 7)
for (curtest in 1:length(py.test)) {
  setwd(input_dir)
  print(curtest)
  if (!file.exists(file.testvec[curtest])) next

  ind.test <- py.test.cols[curtest]
  gm.erc   <- ercmat.gm[, ind.test]

  setwd(output_dir)

  ind.vals   <- as.numeric(gsub("py", "", subnamevec[curtest]))
  thresh.ind <- which(ind.vals == lft.thresh$Pyrome)
  cur.thresh <- thresh.vec[thresh.ind]

  sd.gm.vec[curtest]         <- sd(gm.erc)
  summary.gm.list[[curtest]] <- summary(gm.erc)

  pyrome.vec <- fod.mat$PYROME
  pyrome.unq <- unique(pyrome.vec)
  py.zero    <- str_pad(py.test[curtest], 3, pad = "0")
  fod.sub    <- fod.mat[pyrome.vec == py.test[curtest], ]
  fod.sub    <- fod.sub[fod.sub$FIRE_SIZE >= cur.thresh, ]

  if (nrow(fod.sub) > 0) {
    disc.date.split <- strsplit(fod.sub$DISCOVERY_DATE, " ")
    dmy.vec  <- rep(NA, nrow(fod.sub))
    year.fod <- rep(NA, nrow(fod.sub))
    for (i in seq_along(dmy.vec)) {
      temp.dmy    <- as.Date(disc.date.split[[i]][1], format = "%m/%d/%Y")
      dmy.vec[i]  <- as.character(temp.dmy)
      year.fod[i] <- as.numeric(format(temp.dmy, "%Y"))
    }
    dmy.sub      <- dmy.vec[year.fod %in% yearseq]
    nfires.table <- table(dmy.sub)

    startday <- 1
    dateseq  <- character(ndays.all)
    for (curyear in yearseq) {
      tempseq <- as.character(seq(as.Date(paste0("1/1/", curyear), format = "%m/%d/%Y"),
                                  by = "day", length.out = 365))
      dateseq[startday:(startday + 364)] <- tempseq
      startday <- startday + 365
    }
    all.dmy <- as.Date(dateseq)
    nfbyday <- data.frame(date = all.dmy, nfires = 0L)
    for (curday in seq_along(nfires.table)) {
      tempdate <- names(nfires.table[curday])
      if (yday(tempdate) != 366)
        nfbyday$nfires[nfbyday$date == tempdate] <- nfires.table[curday]
    }

    isfday  <- as.numeric(nfbyday$nfires > 0)
    gridmod <- glm(isfday ~ gm.erc, family = "binomial")
    basemod <- glm(isfday ~ 1,      family = "binomial")

    A      <- logLik(gridmod)
    Base   <- logLik(basemod)
    p.valA <- pchisq(-2 * (as.numeric(A) - as.numeric(Base)), df = 1, lower.tail = FALSE)

    logistic.summary.gm[[curtest]] <- summary(gridmod)

    predictgm <- predict(gridmod, type = "response")
    pgm <- prediction(predictgm, isfday)

    rocgm <- roc(as.numeric(isfday), predictgm,
      percent = FALSE, boot.n = 1000, ci.alpha = 0.9, stratified = FALSE,
      plot = TRUE, grid = TRUE, show.thres = TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
      print.auc = TRUE, print.thres.col = "blue", ci = TRUE, ci.type = "bars",
      print.thres.cex = 0.7)

    auc.gm[curtest] <- rocgm$auc
  }

  if (!is.na(gm.erc[1])) {
    ercdf    <- data.frame(erc = gm.erc, year = as.character(yearseq.all), day = dayseq.all)
    ercdaily <- na.exclude(ercdf) %>%
      group_by(day) %>%
      summarize(mean = mean(erc), sd = sd(erc), .groups = "drop")
    dailymean.gm[, curtest] <- ercdaily$mean
    dailysd.gm[, curtest]   <- ercdaily$sd
  }
}

dev.off()  # close roc-curves-eof.pdf

###--- Compute daily means, SDs and standardised AR residuals ---###
dailymeanmat.gm   <- matrix(NA, nrow = 365,       ncol = nstations)
dailysdmat.gm     <- matrix(NA, nrow = 365,       ncol = nstations)
residsmat.gm      <- matrix(NA, nrow = ndays.all, ncol = nstations)
residsmat.filt.gm <- matrix(NA, nrow = ndays.all, ncol = nstations)
armod.list.gm     <- vector("list", nstations)
sdvec.resids.gm   <- rep(NA, nstations)
sdvec.filt.gm     <- rep(NA, nstations)

for (curstation in 1:nstations) {
  tempvec <- na.exclude(ercmat.gm[, curstation])
  if (sd(tempvec) <= 1e-16) { print(curstation / nstations); next }

  ercdf <- data.frame(erc = ercmat.gm[, curstation],
                      year = as.character(yearseq.all), day = dayseq.all)
  ercsummary <- na.exclude(ercdf) %>%
    group_by(day) %>%
    summarize(mean = mean(erc), sd = sd(erc), .groups = "drop")

  meanvec <- rep(ercsummary$mean, length(yearseq))
  sdvec   <- rep(pmax(ercsummary$sd, 1e-6), length(yearseq))

  dailymeanmat.gm[, curstation] <- ercsummary$mean
  dailysdmat.gm[, curstation]   <- pmax(ercsummary$sd, 1e-6)

  std_residvec <- (ercmat.gm[, curstation] - meanvec) / sdvec
  armod.gm     <- ar(std_residvec, order.max = 30, method = "yule-walker")

  armod.list.gm[[curstation]]     <- armod.gm
  residsmat.gm[, curstation]      <- std_residvec
  residsmat.filt.gm[, curstation] <- armod.gm$resid
  sdvec.resids.gm[curstation]     <- sd(std_residvec)
  sdvec.filt.gm[curstation]       <- sd(na.exclude(armod.gm$resid))
  print(curstation / nstations)
}

###--- Identify valid stations ---### -- they should all be valid but just making sure
validmat.filt.gm <- data.frame(residsmat.filt.gm)
validmat.filt.gm[is.na(validmat.filt.gm)] <- 0

colsd.gm <- apply(validmat.filt.gm, 2, sd); colsd.gm[is.na(colsd.gm)] <- 0
valid.gm <- colsd.gm > 1e-16

nonzero.gm         <- validmat.filt.gm[, valid.gm]
nonzero.mean.gm    <- dailymeanmat.gm[, valid.gm]
nonzero.sd.gm      <- dailysdmat.gm[, valid.gm]
sdresids.valid.gm  <- sdvec.resids.gm[valid.gm]
armod.valid.gm     <- armod.list.gm[valid.gm]
valid.station.inds <- which(valid.gm)
n.valid            <- length(valid.station.inds)

###--- EOF decomposition of AR filter residuals ---###
eof <- prcomp(nonzero.gm, center = TRUE, scale. = FALSE)

var.explained <- eof$sdev^2 / sum(eof$sdev^2)
cum.var       <- cumsum(var.explained)
n.pcs.80      <- min(which(cum.var >= 0.80))
cat(sprintf("%d PCs explain 80%% of variance\n", n.pcs.80))

setwd(output_dir)
png("eof-variance-explained.png", width = 700, height = 400)
plot(cum.var, type = "l",
     xlab = "Number of PCs", ylab = "Cumulative variance explained",
     main = "EOF decomposition of AR filter residuals")
abline(h = 0.80, col = "red", lty = 2)
abline(v = n.pcs.80, col = "blue", lty = 2)
dev.off()

###--- Fit t-distribution to each PC score ---###
n.pcs    <- ncol(eof$x)
df.pc    <- rep(NA, n.pcs)
scale.pc <- rep(NA, n.pcs)
loc.pc   <- rep(NA, n.pcs)

for (pc in 1:n.pcs) {
  scores <- eof$x[, pc]
  tryfit <- try(fitdistr(scores, "t"), silent = TRUE)
  if (!inherits(tryfit, "try-error")) {
    df.pc[pc]    <- tryfit$estimate["df"]
    scale.pc[pc] <- tryfit$estimate["s"]
    loc.pc[pc]   <- tryfit$estimate["m"]
  } else {
    df.pc[pc]    <- 1000
    scale.pc[pc] <- sd(scores)
    loc.pc[pc]   <- mean(scores)
  }
  print(pc / n.pcs)
}

cat("Summary of fitted df across PCs:\n")
print(summary(df.pc))

setwd(output_dir)
png("eof-df-by-pc.png", width = 700, height = 400)
plot(1:n.pcs, pmin(df.pc, 50), type = "o", pch = 16, cex = 0.6,
     xlab = "PC", ylab = "Fitted df (capped at 50)",
     main = "t-distribution df by PC order")
abline(h = 30, col = "red", lty = 2)
abline(v = n.pcs.80, col = "blue", lty = 2)
legend("bottomright", legend = c("df = 30 threshold", "80% variance"),
       col = c("red", "blue"), lty = 2, bty = "n")
dev.off()

###--- Diagnostics (uncomment to inspect specific PCs) ---###
# pcs.diag <- c(1, 2, 3, n.pcs.80, n.pcs.80 + 5)
#
# # QQ plots and ACF for selected PCs
# setwd(output_dir)
# png("eof-pc-diagnostics.png", width = 1400, height = 800)
# par(mfrow = c(3, length(pcs.diag)))
# for (pc in pcs.diag) {
#   scores <- eof$x[, pc]
#   hist(scores, breaks = 50, freq = FALSE,
#        main = sprintf("PC%d (df=%.1f)", pc, df.pc[pc]), xlab = "")
#   curve(dt((x - loc.pc[pc]) / scale.pc[pc], df = df.pc[pc]) / scale.pc[pc],
#         add = TRUE, col = "red", lwd = 2)
#   qqplot(qt(ppoints(length(scores)), df = df.pc[pc]) * scale.pc[pc] + loc.pc[pc],
#          scores, main = paste0("PC", pc, " QQ vs t(", round(df.pc[pc], 1), ")"),
#          xlab = "Theoretical", ylab = "Sample", pch = 16, cex = 0.3)
#   abline(0, 1, col = "red")
#   acf(scores, lag.max = 60, main = paste0("PC", pc, " ACF"))
# }
# dev.off()
#
# # Spatial loading map for leading PCs
# coords <- read.csv(file.path(input_dir, "pyrome-station-coords.csv"))
# pyrome.df.idx  <- match(py.test.cols, valid.station.inds)
# coords$load.pc1 <- eof$rotation[pyrome.df.idx, 1]
# coords$load.pc2 <- eof$rotation[pyrome.df.idx, 2]
# coords <- coords[!is.na(coords$LATITUDE) & !is.na(pyrome.df.idx), ]
# usa <- map_data("state")
# for (pc in 1:2) {
#   load.col <- paste0("load.pc", pc)
#   p <- ggplot() +
#     geom_polygon(data = usa, aes(x = long, y = lat, group = group),
#                  fill = "grey92", colour = "grey60", linewidth = 0.3) +
#     geom_point(data = coords, aes(x = LONGITUDE, y = LATITUDE,
#                colour = coords[[load.col]]), size = 3) +
#     scale_colour_gradient2(name = "Loading", low = "blue", mid = "white", high = "red") +
#     coord_fixed(1.3, xlim = c(-125, -65), ylim = c(24, 50)) +
#     theme_bw() +
#     labs(title = sprintf("EOF%d spatial loading (%.1f%% var)", pc,
#                          var.explained[pc] * 100))
#   ggsave(sprintf("eof-loading-pc%d.png", pc), p, width = 10, height = 6, dpi = 150)
# }

###--- Simulate: draw iid from per-PC t-distributions, reconstruct ---###
# Each PC score is treated as iid in time (the AR filter has already removed
# temporal autocorrelation). Spatial correlation is recovered automatically
# through the EOF rotation.
nyear.sim <- 1001
n.sim     <- nyear.sim * 365

t1 <- Sys.time()
pc.sim <- matrix(NA, nrow = n.sim, ncol = n.pcs)
for (pc in 1:n.pcs) {
  u            <- runif(n.sim)
  pc.sim[, pc] <- qt(u, df = df.pc[pc]) * scale.pc[pc] + loc.pc[pc]
}

# Reconstruct residuals: sim PC scores x rotation^T (+ center, which ≈ 0)
residsim.eof <- pc.sim %*% t(eof$rotation) +
                matrix(eof$center, nrow = n.sim, ncol = n.valid, byrow = TRUE)
t2 <- Sys.time(); print(t2 - t1)

###--- AR-filter simulated innovations and reconstruct ERC ---###
output.ercmean.gm <- matrix(NA, nrow = n.sim, ncol = n.valid)

for (curstation in 1:n.valid) {
  armod <- armod.valid.gm[[curstation]]
  sr    <- sdresids.valid.gm[curstation]

  curfilt <- stats::filter(residsim.eof[, curstation], armod$ar,
                            sides = 1, method = "recursive")
  curfilt <- curfilt * sr / sd(na.exclude(curfilt))

  daily_mean <- rep(nonzero.mean.gm[, curstation], nyear.sim)
  daily_sd   <- rep(nonzero.sd.gm[, curstation],   nyear.sim)

  outvec.erc <- daily_mean + daily_sd * curfilt
  outvec.erc[outvec.erc < 0] <- 0
  output.ercmean.gm[, curstation] <- outvec.erc
  print(curstation / n.valid)
}

###--- Write simulation output ---###
setwd(output_dir)

colnames(output.ercmean.gm) <- namevec[valid.gm]

outmat.gm.sub <- output.ercmean.gm[-(1:365), ]

#write.table(outmat.gm.sub,        "erc-sim-eof-gridmet.csv",     sep = ",", row.names = FALSE)
write.table(round(outmat.gm.sub), "erc-sim-eof-gridmet-int.csv", sep = ",", row.names = FALSE)

var.output.gm    <- var(outmat.gm.sub)
eigen.var.output <- eigen(var.output.gm)
eigen.var.gm     <- eigen(cov(nonzero.gm))

png(file.path(output_dir, "eigenvalue-comparison-eof.png"), width = 600, height = 500)
plot(eigen.var.gm$values, eigen.var.output$values,
     xlab = "Observed eigenvalues", ylab = "Simulated eigenvalues",
     main = "Eigenvalue comparison: observed vs simulated")
abline(0, 1, col = "red")
dev.off()
cat(sprintf("Eigenvalue correlation: %.4f\n", cor(eigen.var.gm$values, eigen.var.output$values)))

###--- Write daily climatology ---###
py.test.names          <- paste0("py", py.test)
colnames(dailymean.gm) <- py.test.names
colnames(dailysd.gm)   <- py.test.names

#write.table(dailymean.gm, "dailymean-gm.csv", sep = ",", row.names = FALSE)
#write.table(dailysd.gm,   "dailysd-gm.csv",   sep = ",", row.names = FALSE)

###--- Simulation diagnostics ---###
med.gm <- apply(outmat.gm.sub, 2, median)
png(file.path(output_dir, "median-erc-eof.png"), width = 700, height = 400)
plot(med.gm, main = "Median simulated ERC by station (EOF)", ylab = "Median ERC")
dev.off()

###--- Spatially coherent seasonal block resampling ---###
# Alternative anomaly model. The EOF/AR/t pipeline above reproduces the observed
# ACF to within 0.006 at lags 1-15, yet still splits above-P80 runs into ~1.7x
# too many episodes at a matched exceedance rate. The reason is not persistence
# and not tail weight, but asymmetry: a linear AR driven by SYMMETRIC innovations
# has exactly zero third-order asymmetry (every third moment vanishes when
# E[eps^3] = 0, so a heavy-tailed but symmetric t buys none of it), while ERC's
# rises and falls are strongly asymmetric - above P80 the mean daily fall is
# ~3.5x the mean daily rise. The ACF cannot see this: autocovariance is symmetric
# in lag, gamma(k) = gamma(-k), so no second-order statistic distinguishes a slow
# ramp then a crash from a crash then a slow ramp. The observed AR innovations
# are skewed (-0.73), which fitting a symmetric t discards by construction.
# (A linear process is strictly time-reversible iff its innovations are Gaussian;
# with symmetric non-Gaussian innovations it is reversible to third order, which
# is the order that run-length geometry depends on.)
# Restoring innovation asymmetry helps but does not rescue this - drawing
# from the *exact* observed innovations closes only ~17% of the run-length gap,
# and the observed innovations are heavy-tailed (kurtosis 4.4), not light-tailed;
# the light-tailed marginal (kurtosis 2.55) is emergent from the saturating
# dynamics, not something an innovation law can install.
#
# Block resampling drops the parametric temporal model. One historical source
# year is drawn per (simulated year, seasonal block) and applied to EVERY pyrome
# at once, so observed cross-pyrome correlation and within-block run-length
# structure both carry through unchanged. blocklen sets the trade-off: shorter
# blocks give more distinct sequences but break more runs at the joins
# (91 days ~ 15^4 distinct combinations per simulated year).
#
# Note this cannot exceed historical ERC *values*. Because run lengths above a
# quantile threshold are invariant to any monotone transform of the marginal,
# extra severity can be added afterwards via a per-day-of-year monotone
# transform at zero cost to the run-length structure recovered here.

block_resample_rows <- function(nyr.out, nyear.hist, blocklen = 91, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # fold the remainder days into the final block: 365 = 4*91 + 1, so a plain
  # integer division would leave a spurious 1-day block at the end of each year
  grp <- pmin(((0:364) %/% blocklen) + 1, 365 %/% blocklen)
  ug  <- unique(grp)
  idx <- integer(nyr.out * 365)
  r   <- 0L
  for (y in seq_len(nyr.out)) {
    for (g in ug) {
      ii  <- which(grp == g)
      idx[r + seq_along(ii)] <- (sample.int(nyear.hist, 1) - 1L) * 365L + ii
      r   <- r + length(ii)
    }
  }
  idx
}

blocklen.sim <- 91
nyear.blk    <- nyear.sim - 1   # no burn-in needed: no recursive filter transient
row.idx.blk  <- block_resample_rows(nyear.blk, length(yearseq),
                                    blocklen = blocklen.sim, seed = 20260814)

outmat.blk <- ercmat.gm[row.idx.blk, valid.station.inds, drop = FALSE]
colnames(outmat.blk)  <- namevec[valid.gm]
storage.mode(outmat.blk) <- "integer"   # observed ERC is integer-valued already

setwd(output_dir)
write.table(outmat.blk, "erc-sim-block-gridmet-int.csv", sep = ",", row.names = FALSE)

# Spatial structure: block resampling reproduces it by construction, whereas
# applying station-specific AR filters *after* imposing spatial correlation on
# the innovations decoheres the field (differing transfer functions per station).
obsmat.valid   <- ercmat.gm[, valid.station.inds, drop = FALSE]
cor.obs        <- cor(obsmat.valid)
eigen.var.blk  <- eigen(var(outmat.blk))
eigen.var.obs  <- eigen(var(obsmat.valid))
cat(sprintf("Correlation-matrix RMSE vs observed:  EOF=%.4f   block=%.4f\n",
            sqrt(mean((cor(outmat.gm.sub) - cor.obs)^2)),
            sqrt(mean((cor(outmat.blk)    - cor.obs)^2))))
cat(sprintf("Leading eigenvalue: observed=%.0f  EOF=%.0f  block=%.0f\n",
            eigen.var.obs$values[1], eigen.var.output$values[1],
            eigen.var.blk$values[1]))

###--- Duration-dependent drift model ---###
# Keeps everything the EOF pipeline gets right - the same spatially correlated
# innovations (residsim.eof) and the same per-pyrome AR coefficients - but
# replaces the purely linear recursion with
#   x_t = sum_j phi_j x_{t-j} + g(state_{t-1}, dur_{t-1}) + eps_t
# where state = 1{ERC >= P80} and dur = consecutive days held in that state.
#
# g is estimable precisely because Yule-Walker only forces E[e_t x_{t-j}] = 0.
# State and duration are NONLINEAR functions of the past, so E[e_t | state, dur]
# is free to be non-zero, and that conditional mean is exactly the asymmetry the
# linear model discards. The feedback (dur depends on past states, which depend
# on past x) makes the recursion time-irreversible - the missing ingredient the
# hazard diagnostic points to. Unlike a two-stage probit this needs no coupling
# of 128 binary fields and leaves the continuous marginal free to settle where
# the dynamics put it, so the ERC values the downstream logistic fire model
# consumes are still generated rather than imposed.

dur.breaks <- c(1, 2, 3, 4, 7, 14, Inf)
n.buck     <- length(dur.breaks) - 1L
DMAX       <- 400L
buck.lut   <- as.integer(cut(1:DMAX, dur.breaks, right = FALSE, include.lowest = TRUE))

# P80 in ERC space per pyrome (reused by the acceptance test below)
thr.p80 <- vapply(valid.station.inds, function(s) {
  v  <- na.exclude(ercmat.gm[, s]); su <- sort(unique(v))
  su[min(which(vapply(su, function(x) mean(v <= x), 0) >= 0.80))]
}, numeric(1))

# on the standardised scale the state threshold varies with day-of-year
q.obs <- matrix(NA_real_, 365, n.valid)
for (s in 1:n.valid)
  q.obs[, s] <- (thr.p80[s] - nonzero.mean.gm[, s]) / nonzero.sd.gm[, s]

resids.valid <- residsmat.gm[, valid.gm, drop = FALSE]

## --- estimate g from the observed innovations ------------------------------
g.arr <- array(0, dim = c(n.valid, 2L, n.buck))
for (s in 1:n.valid) {
  x   <- resids.valid[, s]
  phi <- armod.valid.gm[[s]]$ar
  p   <- length(phi)
  X   <- embed(x, p + 1L)
  e   <- as.numeric(X[, 1] - X[, -1, drop = FALSE] %*% phi)
  st  <- as.integer(x >= rep(q.obs[, s], length(yearseq)))
  rr  <- rle(st)
  dv  <- unlist(lapply(rr$lengths, seq_len))   # days so far in the current state
  ii  <- (p + 1L):length(x)
  f1  <- factor(st[ii - 1L] + 1L, levels = 1:2)
  f2  <- factor(buck.lut[pmin(dv[ii - 1L], DMAX)], levels = 1:n.buck)
  ag  <- tapply(e, list(f1, f2), mean);   ag[is.na(ag)]   <- 0
  cnt <- tapply(e, list(f1, f2), length); cnt[is.na(cnt)] <- 0
  g.arr[s, , ] <- ag - sum(ag * cnt) / sum(cnt)   # centre so mean drift is zero
}
dimnames(g.arr) <- list(NULL, c("below", "above"),
                        c("1", "2", "3", "4-6", "7-13", "14+"))
cat("\nFitted drift g by state x duration (mean over pyromes, sd units):\n")
print(round(apply(g.arr, c(2, 3), mean), 4))

## --- simulate (vectorised across pyromes; t must stay sequential) ----------
pmax.ar <- max(vapply(armod.valid.gm, function(a) length(a$ar), 1L))
PHI     <- matrix(0, pmax.ar, n.valid)
for (s in 1:n.valid) { ph <- armod.valid.gm[[s]]$ar; PHI[seq_along(ph), s] <- ph }

q.sim    <- q.obs[rep(1:365, nyear.sim), , drop = FALSE]
xs       <- matrix(0, n.sim, n.valid)
guard    <- 8 * sdresids.valid.gm      # stability clamp on the feedback loop
st.prev  <- integer(n.valid)
dur.prev <- rep(1L, n.valid)

t1 <- Sys.time()
for (t in (pmax.ar + 1L):n.sim) {
  lin      <- colSums(PHI * xs[(t - 1L):(t - pmax.ar), , drop = FALSE])
  gk       <- g.arr[cbind(1:n.valid, st.prev + 1L, buck.lut[pmin(dur.prev, DMAX)])]
  v        <- pmax(pmin(lin + gk + residsim.eof[t, ], guard), -guard)
  xs[t, ]  <- v
  st.now   <- as.integer(v >= q.sim[t, ])
  dur.prev <- ifelse(st.now == st.prev, dur.prev + 1L, 1L)
  st.prev  <- st.now
}
cat("duration-drift recursion: "); print(Sys.time() - t1)

# Deliberately NOT rescaled to sdresids.valid.gm the way the EOF path is: the
# state threshold is evaluated during the recursion, so a post-hoc global rescale
# would silently change which days were above P80. Report the ratio instead - if
# it drifts far from 1 the drift term is inflating variance and the exceedance
# rate below will show it.
sd.ratio <- vapply(1:n.valid, function(s) sd(xs[, s]), 0) / sdresids.valid.gm
cat(sprintf("sd(simulated anomaly)/sd(observed anomaly): mean=%.3f  range=[%.3f, %.3f]\n",
            mean(sd.ratio), min(sd.ratio), max(sd.ratio)))

output.dd <- matrix(NA_real_, n.sim, n.valid)
for (s in 1:n.valid) {
  e <- rep(nonzero.mean.gm[, s], nyear.sim) + rep(nonzero.sd.gm[, s], nyear.sim) * xs[, s]
  e[e < 0] <- 0
  output.dd[, s] <- e
}
outmat.dd <- round(output.dd[-(1:365), , drop = FALSE])
colnames(outmat.dd)      <- namevec[valid.gm]
storage.mode(outmat.dd)  <- "integer"
write.table(outmat.dd, "erc-sim-durdrift-gridmet-int.csv", sep = ",", row.names = FALSE)

p.dd <- vapply(1:n.valid, function(s) mean(outmat.dd[, s]    >= thr.p80[s]), 0)
p.ob <- vapply(1:n.valid, function(s) mean(obsmat.valid[, s] >= thr.p80[s]), 0)
p.ef <- vapply(1:n.valid, function(s) mean(outmat.gm.sub[, s] >= thr.p80[s]), 0)
cat(sprintf("exceedance rate at P80: observed=%.4f  EOF=%.4f  durdrift=%.4f\n",
            mean(p.ob), mean(p.ef), mean(p.dd)))
cat(sprintf("Correlation-matrix RMSE vs observed: durdrift=%.4f\n",
            sqrt(mean((cor(outmat.dd) - cor.obs)^2))))
rm(output.dd); gc(verbose = FALSE)

###--- State-dependent innovation skewness (conditional anamorphosis) ---###
# The rise/fall asymmetry is NOT in the innovation variance. Conditional on the
# lagged level, sd is U-shaped and actually HIGHEST at the top decile (0.826 vs
# 0.642 at the median, ratio 1.29), so damping variance when high would be
# anti-realistic. The asymmetry is in the innovation SKEWNESS, which swings from
# +0.27 at low levels to -1.25 at high levels: once ERC is high, sharp crashes
# are available but large further rises are not, because the fuel is near
# saturation. A symmetric t has zero skewness at every level, which is why the
# EOF pipeline reproduces none of this.
#
# Construction, a copula-style conditional anamorphosis:
#   1. take the spatially correlated EOF innovation for this pyrome-day
#   2. convert it to a uniform via its own marginal (rank transform), which
#      leaves the cross-pyrome rank dependence untouched
#   3. map that uniform through the empirical quantile function of the OBSERVED
#      innovations for this pyrome's current (level, duration) cell
# Step 3 installs the state-dependent skewness. Steps 2-3 are monotone per
# pyrome so the spatial rank structure survives; realised Pearson correlation is
# measured below rather than assumed.
#
# Memory: this holds another n.sim x n.valid matrix. At nyear.sim = 1001 that is
# ~370 MB on top of the EOF, block and drift outputs - free what you don't need
# or lower nyear.sim if the machine is tight.

nlev.cs <- 8L; ndur.cs <- 3L
ncell   <- nlev.cs * ndur.cs

BRK   <- matrix(NA_real_, nlev.cs - 1L, n.valid)
POOLl <- vector("list", n.valid * ncell)
LENv  <- integer(n.valid * ncell)
kk    <- 0L
for (s in 1:n.valid) {
  x   <- resids.valid[, s]
  phi <- armod.valid.gm[[s]]$ar
  p   <- length(phi)
  X   <- embed(x, p + 1L)
  e   <- as.numeric(X[, 1] - X[, -1, drop = FALSE] %*% phi)
  l   <- X[, 2]                                   # lagged level
  st  <- as.integer(x >= rep(q.obs[, s], length(yearseq)))
  rr  <- rle(st)
  dv  <- unlist(lapply(rr$lengths, seq_len))
  d   <- dv[(p + 1L):length(x) - 1L]
  br  <- as.numeric(quantile(l, seq(0, 1, length.out = nlev.cs + 1L))[-c(1, nlev.cs + 1L)])
  BRK[, s] <- br
  lb  <- findInterval(l, br) + 1L                 # duplicated breaks -> empty cells
  db  <- ifelse(d < 2L, 1L, ifelse(d < 5L, 2L, 3L))
  for (a in 1:nlev.cs) for (b in 1:ndur.cs) {
    ii <- which(lb == a & db == b)
    if (length(ii) < 25) ii <- which(lb == a)     # fall back to level-only
    if (length(ii) < 5)  ii <- seq_along(e)       # then to unconditional
    kk <- kk + 1L
    POOLl[[kk]] <- sort(e[ii]); LENv[kk] <- length(ii)
  }
}
POOL <- unlist(POOLl); rm(POOLl)
OFFv <- c(0L, cumsum(LENv)[-length(LENv)])        # flat index order: s, then a, then b

# consume residsim.eof in place: convert to uniforms via per-pyrome ranks
for (s in 1:n.valid)
  residsim.eof[, s] <- rank(residsim.eof[, s], ties.method = "average") / (n.sim + 1)

xs.cs    <- matrix(0, n.sim, n.valid)
st.prev  <- integer(n.valid)
dur.prev <- rep(1L, n.valid)
sidx     <- 1:n.valid
soff     <- (sidx - 1L) * ncell
t1 <- Sys.time()
for (t in (pmax.ar + 1L):n.sim) {
  lev  <- colSums(BRK < rep(xs.cs[t - 1L, ], each = nlev.cs - 1L)) + 1L
  dbk  <- ifelse(dur.prev < 2L, 1L, ifelse(dur.prev < 5L, 2L, 3L))
  fidx <- dbk + (lev - 1L) * ndur.cs + soff
  eps  <- POOL[OFFv[fidx] + pmax(1L, ceiling(residsim.eof[t, ] * LENv[fidx]))]
  lin  <- colSums(PHI * xs.cs[(t - 1L):(t - pmax.ar), , drop = FALSE])
  v    <- pmax(pmin(lin + eps, guard), -guard)
  xs.cs[t, ] <- v
  st.now   <- as.integer(v >= q.sim[t, ])
  dur.prev <- ifelse(st.now == st.prev, dur.prev + 1L, 1L)
  st.prev  <- st.now
}
cat("conditional-anamorphosis recursion: "); print(Sys.time() - t1)

output.cs <- matrix(NA_real_, n.sim, n.valid)
for (s in 1:n.valid) {
  e <- rep(nonzero.mean.gm[, s], nyear.sim) + rep(nonzero.sd.gm[, s], nyear.sim) * xs.cs[, s]
  e[e < 0] <- 0
  output.cs[, s] <- e
}
outmat.cs <- round(output.cs[-(1:365), , drop = FALSE])
colnames(outmat.cs)     <- namevec[valid.gm]
storage.mode(outmat.cs) <- "integer"
rm(output.cs); gc(verbose = FALSE)
write.table(outmat.cs, "erc-sim-condskew-gridmet-int.csv", sep = ",", row.names = FALSE)

p.cs  <- vapply(1:n.valid, function(s) mean(outmat.cs[, s] >= thr.p80[s]), 0)
cs.sd <- vapply(1:n.valid, function(s) sd(xs.cs[, s]), 0) / sdresids.valid.gm
cat(sprintf("sd ratio: mean=%.3f  exceedance rate=%.4f (observed %.4f)\n",
            mean(cs.sd), mean(p.cs), mean(p.ob)))
cat(sprintf("Correlation-matrix RMSE vs observed: condskew=%.4f\n",
            sqrt(mean((cor(outmat.cs) - cor.obs)^2))))

# did the asymmetry actually make it into the output?
rf_ratio <- function(mat, thr) {
  r <- vapply(seq_along(thr), function(s) {
    x <- as.numeric(mat[, s]); dd <- diff(x); hi <- which(x[-length(x)] >= thr[s])
    if (length(hi) < 10) return(c(NA_real_, NA_real_))
    c(mean(dd[hi][dd[hi] > 0]), mean(dd[hi][dd[hi] < 0]))
  }, numeric(2))
  mean(abs(r[2, ]), na.rm = TRUE) / mean(r[1, ], na.rm = TRUE)
}
cat(sprintf("above-P80 |fall|/rise: observed=%.3f  EOF=%.3f  condskew=%.3f  block=%.3f\n",
            rf_ratio(obsmat.valid,  thr.p80), rf_ratio(outmat.gm.sub, thr.p80),
            rf_ratio(outmat.cs,     thr.p80), rf_ratio(outmat.blk,    thr.p80)))

###--- Burn block analysis ---###
burn_blocks <- function(erc_vec, ref_vec = NULL, p_thresh = 0.80) {
  ref   <- na.exclude(if (is.null(ref_vec)) erc_vec else ref_vec)
  if (length(ref) == 0) return(list(above = integer(0), below = integer(0)))
  above <- ecdf(ref)(erc_vec) >= p_thresh
  runs  <- rle(as.vector(na.exclude(above)))
  list(
    above = runs$lengths[ runs$values],
    below = runs$lengths[!runs$values]
  )
}

obs_blocks <- lapply(valid.station.inds, function(s) {
  burn_blocks(ercmat.gm[, s])
})

sim_blocks <- lapply(seq_along(valid.station.inds), function(s) {
  burn_blocks(outmat.gm.sub[, s], ref_vec = ercmat.gm[, valid.station.inds[s]])
})

blk_blocks <- lapply(seq_along(valid.station.inds), function(s) {
  burn_blocks(outmat.blk[, s], ref_vec = ercmat.gm[, valid.station.inds[s]])
})

dd_blocks <- lapply(seq_along(valid.station.inds), function(s) {
  burn_blocks(outmat.dd[, s], ref_vec = ercmat.gm[, valid.station.inds[s]])
})

cs_blocks <- lapply(seq_along(valid.station.inds), function(s) {
  burn_blocks(outmat.cs[, s], ref_vec = ercmat.gm[, valid.station.inds[s]])
})

obs_above <- unlist(lapply(obs_blocks, `[[`, "above"))
obs_below <- unlist(lapply(obs_blocks, `[[`, "below"))
sim_above <- unlist(lapply(sim_blocks, `[[`, "above"))
sim_below <- unlist(lapply(sim_blocks, `[[`, "below"))
blk_above <- unlist(lapply(blk_blocks, `[[`, "above"))
blk_below <- unlist(lapply(blk_blocks, `[[`, "below"))
dd_above  <- unlist(lapply(dd_blocks,  `[[`, "above"))
dd_below  <- unlist(lapply(dd_blocks,  `[[`, "below"))
cs_above  <- unlist(lapply(cs_blocks,  `[[`, "above"))
cs_below  <- unlist(lapply(cs_blocks,  `[[`, "below"))

cat("Observed above-P80 block lengths:\n");   print(summary(obs_above))
cat("EOF above-P80 block lengths:\n");        print(summary(sim_above))
cat("Dur-drift above-P80 block lengths:\n");  print(summary(dd_above))
cat("Cond-skew above-P80 block lengths:\n");  print(summary(cs_above))
cat("Block-resampled above-P80 lengths:\n");  print(summary(blk_above))
cat("Observed below-P80 block lengths:\n");   print(summary(obs_below))
cat("EOF below-P80 block lengths:\n");        print(summary(sim_below))
cat("Dur-drift below-P80 block lengths:\n");  print(summary(dd_below))
cat("Cond-skew below-P80 block lengths:\n");  print(summary(cs_below))
cat("Block-resampled below-P80 lengths:\n");  print(summary(blk_below))

ab.list <- list(observed = obs_above, EOF = sim_above, durdrift = dd_above,
                condskew = cs_above, block = blk_above)
cat(sprintf("\n%-22s %9s %9s %9s %9s %9s\n", "above-P80 statistic",
            names(ab.list)[1], names(ab.list)[2], names(ab.list)[3],
            names(ab.list)[4], names(ab.list)[5]))
fmt_row <- function(label, f, digits = 3) {
  v <- vapply(ab.list, f, 0)
  cat(sprintf(paste0("%-22s", paste(rep(paste0(" %9.", digits, "f"), 5), collapse = "")), label,
              v[1], v[2], v[3], v[4], v[5]), "\n", sep = "")
}
fmt_row("mean log run length", function(v) mean(log(v)))
fmt_row("mean run length",     function(v) mean(v))
fmt_row("median run length",   function(v) median(v), 1)
fmt_row("% episodes >=14d",    function(v) mean(v >= 14) * 100, 2)
fmt_row("% episodes >=30d",    function(v) mean(v >= 30) * 100, 2)

setwd(output_dir)

png("burn-blocks-above-p80-eof.png", width = 1300, height = 760)
par(mfrow = c(2, 3))
hist(obs_above, breaks = 30, main = "Observed: runs above P80",
     xlab = "Block length (days)", col = "steelblue", border = "white", xlim = c(0, 150))
hist(sim_above, breaks = 30, main = "EOF/AR/t: runs above P80",
     xlab = "Block length (days)", col = "tomato",    border = "white", xlim = c(0, 150))
hist(dd_above,  breaks = 30, main = "Duration-drift: runs above P80",
     xlab = "Block length (days)", col = "mediumpurple", border = "white", xlim = c(0, 150))
hist(cs_above,  breaks = 30, main = "Cond-skew: runs above P80",
     xlab = "Block length (days)", col = "goldenrod", border = "white", xlim = c(0, 150))
hist(blk_above, breaks = 30, main = "Block resample: runs above P80",
     xlab = "Block length (days)", col = "darkseagreen", border = "white", xlim = c(0, 150))
dev.off()

png("burn-blocks-below-p80-eof.png", width = 1300, height = 760)
par(mfrow = c(2, 3))
hist(obs_below, breaks = 30, main = "Observed: runs below P80",
     xlab = "Block length (days)", col = "steelblue", border = "white", xlim = c(0, 500))
hist(sim_below, breaks = 30, main = "EOF/AR/t: runs below P80",
     xlab = "Block length (days)", col = "tomato",    border = "white", xlim = c(0, 500))
hist(dd_below,  breaks = 30, main = "Duration-drift: runs below P80",
     xlab = "Block length (days)", col = "mediumpurple", border = "white", xlim = c(0, 500))
hist(cs_below,  breaks = 30, main = "Cond-skew: runs below P80",
     xlab = "Block length (days)", col = "goldenrod", border = "white", xlim = c(0, 500))
hist(blk_below, breaks = 30, main = "Block resample: runs below P80",
     xlab = "Block length (days)", col = "darkseagreen", border = "white", xlim = c(0, 500))
dev.off()

grp.cols  <- c("steelblue", "tomato", "mediumpurple", "goldenrod", "darkseagreen")
grp.names <- c("observed", "EOF", "durdrift", "condskew", "block")

above_vec <- c(log(obs_above), log(sim_above), log(dd_above),
               log(cs_above), log(blk_above))
os_vec    <- factor(rep(grp.names, c(length(obs_above), length(sim_above),
                                     length(dd_above), length(cs_above),
                                     length(blk_above))),
                    levels = grp.names)
png("burn-block-above-eof.png", width = 600, height = 450)
par(mfrow = c(1, 1))
boxplot(above_vec ~ os_vec, main = "Above-P80 run lengths", xlab = "",
        ylab = "log run length (days)", col = grp.cols)
dev.off()

below_vec <- c(log(obs_below), log(sim_below), log(dd_below),
               log(cs_below), log(blk_below))
os_vec.b  <- factor(rep(grp.names, c(length(obs_below), length(sim_below),
                                     length(dd_below), length(cs_below),
                                     length(blk_below))),
                    levels = grp.names)
png("burn-block-below-eof.png", width = 600, height = 450)
par(mfrow = c(1, 1))
boxplot(below_vec ~ os_vec.b, main = "Below-P80 run lengths", xlab = "",
        ylab = "log run length (days)", col = grp.cols)
dev.off()

###--- ECDF comparison ---###
png("burn-block-ecdf-eof.png", width = 700, height = 450)
ec <- list(ecdf(obs_above), ecdf(sim_above), ecdf(dd_above),
           ecdf(cs_above),  ecdf(blk_above))
plot(ec[[1]], main = "ECDF of above-P80 block lengths", xlab = "Block length (days)",
     ylab = "Cumulative probability", xlim = c(0, 60), col = grp.cols[1])
for (i in 2:5) plot(ec[[i]], add = TRUE, col = grp.cols[i])
legend("bottomright", c("observed", "EOF/AR/t", "duration-drift",
                        "cond-skew", "block resample"),
       col = grp.cols, lty = 1, bty = "n")
dev.off()

###--- Hazard of run termination: where the fragmentation happens ---###
# The EOF model ends nascent runs at ~1.8x the observed rate on day 1, then the
# hazards converge by day 10. It is a startup-stickiness failure, not a tail
# failure, which is why the longest simulated runs look fine.
hazard  <- function(v, K = 30) vapply(1:K, function(k) sum(v == k) / sum(v >= k), 0)
haz.mat <- cbind(observed = hazard(obs_above),
                 EOF      = hazard(sim_above),
                 durdrift = hazard(dd_above),
                 condskew = hazard(cs_above),
                 block    = hazard(blk_above))
png("burn-block-hazard.png", width = 700, height = 450)
matplot(1:30, haz.mat, type = "o", pch = 16, cex = 0.6, lty = 1, col = grp.cols,
        xlab = "Run length k (days)", ylab = "P(run ends at k | run >= k)",
        main = "Hazard of above-P80 run termination")
legend("topright", colnames(haz.mat), col = grp.cols, lty = 1, pch = 16, bty = "n")
dev.off()
cat("\nHazard of above-P80 run termination:\n")
print(round(haz.mat[c(1, 2, 3, 5, 7, 10, 14, 20, 30), ], 4))

###--- Acceptance test: is the observed record a plausible 15-yr draw? ---###
# The t-tests below are kept for continuity but cannot answer that question.
# They compare ~32k observed episodes against millions of simulated ones and
# treat every episode as independent, when episodes are correlated both serially
# within a pyrome and spatially across pyromes. With n that lopsided, any
# trivial difference returns p < 2e-16 (it does: t ~ 93).
#
# The correct reference is the model's OWN 15-year sampling distribution, which
# preserves the spatial dependence, the serial dependence, and the run
# truncation that a short record imposes. If the observed statistic falls
# outside that band, the discrepancy is structural rather than sampling noise.

# thr.p80 (integer P80 per pyrome) is defined above, with the duration-drift model
window_stat <- function(mat, thr) {
  a <- unlist(lapply(seq_along(thr), function(s) {
    r <- rle(mat[, s] >= thr[s]); r$lengths[r$values]
  }))
  c(meanlog = mean(log(a)), mean = mean(a), median = median(a),
    p90 = unname(quantile(a, 0.90)), frac30 = mean(a >= 30))
}

acceptance_null <- function(mat, nyr.win = length(yearseq), nwin = 200,
                            thr = thr.p80, seed = 1) {
  set.seed(seed)
  nyr.tot <- nrow(mat) %/% 365
  starts  <- sample.int(nyr.tot - nyr.win + 1, nwin, replace = TRUE)
  t(vapply(starts, function(y) {
    window_stat(mat[((y - 1) * 365 + 1):((y - 1 + nyr.win) * 365), , drop = FALSE], thr)
  }, numeric(5)))
}

report_acceptance <- function(nullmat, obsvec, label) {
  cat("\n", label, "\n", sep = "")
  for (v in names(obsvec)) {
    d <- nullmat[, v]; o <- obsvec[[v]]
    cat(sprintf("  %-7s observed=%9.3f | model 15-yr window: mean=%9.3f  95%%CI=[%9.3f,%9.3f]  observed pctile=%5.1f%%\n",
                v, o, mean(d), quantile(d, 0.025), quantile(d, 0.975),
                mean(d <= o) * 100))
  }
}

obs.win <- window_stat(obsmat.valid, thr.p80)
report_acceptance(acceptance_null(outmat.gm.sub), obs.win, "EOF/AR/t model:")
report_acceptance(acceptance_null(outmat.dd),     obs.win, "Duration-drift model:")
report_acceptance(acceptance_null(outmat.cs),     obs.win, "Cond-skew model:")
report_acceptance(acceptance_null(outmat.blk),    obs.win, "Block resample model:")

t.test(log(obs_above), log(sim_above))
t.test(log(obs_below), log(sim_below))

