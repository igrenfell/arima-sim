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

# seed so runs are comparable: this single draw feeds residsim.eof, and
# therefore the EOF, duration-drift and regime-switching outputs alike
set.seed(20260820)

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


###--- Write daily climatology ---###
py.test.names          <- paste0("py", py.test)
colnames(dailymean.gm) <- py.test.names
colnames(dailysd.gm)   <- py.test.names

#write.table(dailymean.gm, "dailymean-gm.csv", sep = ",", row.names = FALSE)
#write.table(dailysd.gm,   "dailysd-gm.csv",   sep = ",", row.names = FALSE)


###--- Shared setup for the regime-switching model ---###
# These were previously defined inside the block-resampling and duration-drift
# sections; both models are gone, so the definitions they contributed live here.

obsmat.valid <- ercmat.gm[, valid.station.inds, drop = FALSE]
cor.obs      <- cor(obsmat.valid)
resids.valid <- residsmat.gm[, valid.gm, drop = FALSE]

# P80 in ERC space per pyrome (reused by the acceptance test below)
thr.p80 <- vapply(valid.station.inds, function(s) {
  v  <- na.exclude(ercmat.gm[, s]); su <- sort(unique(v))
  su[min(which(vapply(su, function(x) mean(v <= x), 0) >= 0.80))]
}, numeric(1))

# on the standardised scale the state threshold varies with day-of-year
q.obs <- matrix(NA_real_, 365, n.valid)
for (s in 1:n.valid)
  q.obs[, s] <- (thr.p80[s] - nonzero.mean.gm[, s]) / nonzero.sd.gm[, s]

q.sim <- q.obs[rep(1:365, nyear.sim), , drop = FALSE]
guard <- 8 * sdresids.valid.gm      # stability clamp on the feedback loop
p.ob  <- vapply(1:n.valid, function(s) mean(obsmat.valid[, s] >= thr.p80[s]), 0)


###--- Regime-switching AR + state-dependent innovation skewness ---###
# The rise/fall asymmetry is NOT in the innovation variance. Conditional on the
# lagged level, sd is U-shaped and actually HIGHEST at the top decile (0.826 vs
# 0.642 at the median, ratio 1.29), so damping variance when high would be
# anti-realistic. The asymmetry is in the innovation SKEWNESS, which swings from
# +0.27 at low levels to -1.25 at high levels: once ERC is high, sharp crashes
# are available but large further rises are not, because the fuel is near
# saturation. A symmetric t has zero skewness at every level, which is why the
# EOF pipeline reproduces none of this.
#
# Two mechanisms together (each alone plateaus: varying phi only = 29%, varying
# eps only = 49%; both = 54%):
#  (i) regime-dependent AR coefficients, fit separately above/below threshold.
#      A FIXED phi was the binding constraint - at day 1 of a run the linear
#      predictor averages the previous p days, which were mostly BELOW
#      threshold, so it drags the process back down and no innovation draw can
#      compensate. The fitted above-regime AR is genuinely more persistent.
# (ii) a copula-style conditional anamorphosis for the innovation:
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

# Regime AR order. There is a real trade-off here, measured over 300 sim years:
#   p    %run-length gap closed    acf10   acf20   acf40      (observed: .312 .235 .163)
#   1           66.4               0.114   0.050   0.016
#   3           62.4               0.190   0.078   0.019
#   10          57.6               0.324   0.176   0.061
#   20          54.3               0.320   0.241   0.114
# Short orders post the best run-length number by DESTROYING the long memory
# that the original pipeline already got right - that is relocating the error,
# not fixing it. p = 20 keeps the ACF close through lag 20 and still beats the
# fixed-phi version (49.5% -> 54.3%). Raise only with the ACF columns in view.
p.reg   <- 20L
ndur.cs <- 6L

# level bins are on the MARGIN from threshold (x_{t-1} - q), not the raw level:
# q varies with day-of-year, so a given level is above threshold in summer and
# below it in winter. Quantile bins are also threshold-agnostic - with P80 the
# threshold falls inside bin 7 of 8, so 61% of above-threshold days shared a
# pool with below-threshold days. Anchoring a breakpoint at 0 fixes that.
mbrk.q  <- c(0.20, 0.45, 0.65)
mbrk.f  <- c(-0.5, -0.2, 0, 0.2, 0.5, 1.0)
durbin  <- function(d) ifelse(d < 2L, 1L, ifelse(d < 3L, 2L, ifelse(d < 4L, 3L,
                       ifelse(d < 7L, 4L, ifelse(d < 14L, 5L, 6L)))))

B0 <- matrix(0, p.reg + 1L, n.valid)   # below-threshold regime: intercept + phi
B1 <- matrix(0, p.reg + 1L, n.valid)   # above-threshold regime
Elist <- Mlist <- Dlist <- vector("list", n.valid)
for (s in 1:n.valid) {
  x  <- resids.valid[, s]; n <- length(x)
  X  <- embed(x, p.reg + 1L); y <- X[, 1]; Lg <- X[, -1, drop = FALSE]
  lag.idx <- seq_len(n - p.reg) + p.reg - 1L          # index of t-1
  qq <- q.obs[((lag.idx - 1L) %% 365L) + 1L, s]
  st <- as.integer(x >= rep(q.obs[, s], length(yearseq)))
  rr <- rle(st); dv <- unlist(lapply(rr$lengths, seq_len))
  reg <- st[lag.idx]                                  # regime from state at t-1
  Z  <- cbind(1, Lg)
  fitreg <- function(k) {
    i <- which(reg == k)
    if (length(i) < 5 * (p.reg + 1L)) return(NULL)
    tryCatch(qr.solve(Z[i, , drop = FALSE], y[i]), error = function(e) NULL)
  }
  b0 <- fitreg(0L); b1 <- fitreg(1L)
  if (is.null(b0)) b0 <- qr.solve(Z, y)
  if (is.null(b1)) b1 <- b0
  B0[, s] <- b0; B1[, s] <- b1
  pred <- ifelse(reg == 1L, as.numeric(Z %*% b1), as.numeric(Z %*% b0))
  Elist[[s]] <- y - pred
  Mlist[[s]] <- Lg[, 1] - qq
  Dlist[[s]] <- dv[lag.idx]
}
cat(sprintf("regime AR persistence: below sum(phi)=%.4f  above sum(phi)=%.4f\n",
            mean(colSums(B0[-1, , drop = FALSE])), mean(colSums(B1[-1, , drop = FALSE]))))

BRK.cs  <- sort(unique(c(as.numeric(quantile(unlist(Mlist), mbrk.q)), mbrk.f)))
nlev.cs <- length(BRK.cs) + 1L
ncell   <- nlev.cs * ndur.cs
POOLl <- vector("list", n.valid * ncell)
LENv  <- integer(n.valid * ncell)
kk    <- 0L
for (s in 1:n.valid) {
  lb <- findInterval(Mlist[[s]], BRK.cs) + 1L
  db <- durbin(Dlist[[s]]); e <- Elist[[s]]
  for (a in 1:nlev.cs) for (b in 1:ndur.cs) {
    ii <- which(lb == a & db == b)
    if (length(ii) < 25) ii <- which(lb == a)     # fall back to margin-only
    if (length(ii) < 5)  ii <- seq_along(e)       # then to unconditional
    kk <- kk + 1L
    POOLl[[kk]] <- sort(e[ii]); LENv[kk] <- length(ii)
  }
}
POOL <- unlist(POOLl); rm(POOLl, Elist, Mlist, Dlist)
OFFv <- c(0L, cumsum(LENv)[-length(LENv)])        # flat index order: s, then a, then b

# consume residsim.eof in place: convert to uniforms via per-pyrome ranks
for (s in 1:n.valid)
  residsim.eof[, s] <- rank(residsim.eof[, s], ties.method = "average") / (n.sim + 1)

xs.cs    <- matrix(0, n.sim, n.valid)
st.prev  <- integer(n.valid)
dur.prev <- rep(1L, n.valid)
sidx     <- 1:n.valid
soff     <- (sidx - 1L) * ncell
a0 <- B0[1, ]; a1 <- B1[1, ]
C0 <- B0[-1, , drop = FALSE]; C1 <- B1[-1, , drop = FALSE]
t1 <- Sys.time()
for (t in (p.reg + 1L):n.sim) {
  xl   <- xs.cs[(t - 1L):(t - p.reg), , drop = FALSE]
  marg <- xs.cs[t - 1L, ] - q.sim[t - 1L, ]
  # regime-dependent AR: both predictions are formed and selected, which is
  # cheaper per step than rebuilding a coefficient matrix 365k times
  lin  <- ifelse(marg >= 0, a1 + colSums(C1 * xl), a0 + colSums(C0 * xl))
  lev  <- findInterval(marg, BRK.cs) + 1L
  fidx <- durbin(dur.prev) + (lev - 1L) * ndur.cs + soff
  eps  <- POOL[OFFv[fidx] + pmax(1L, ceiling(residsim.eof[t, ] * LENv[fidx]))]
  v    <- pmax(pmin(lin + eps, guard), -guard)
  xs.cs[t, ] <- v
  st.now   <- as.integer(v >= q.sim[t, ])
  dur.prev <- ifelse(st.now == st.prev, dur.prev + 1L, 1L)
  st.prev  <- st.now
}
cat("regime-switching recursion: "); print(Sys.time() - t1)

# Physical ceiling on ERC. ERC takes no wind or slope input, so it is a function
# of the fuel-moisture classes against fixed fuel-model parameters and attains a
# maximum when every class sits at its dry floor. That value was computed by
# running NFDRS4_cli on this project's own NFDRSInit config against 4 years of
# synthetic extreme heat / near-zero RH / zero precipitation: all four dead
# classes equilibrate at 2.35%, herb at 30%, woody at 60%, GSI at 0, and ERC
# plateaus at 121.08 (daily value at obsHour; 125.85 as an hourly maximum).
#
# Independent corroboration: the observed maximum in the BASE_2020 GridMet
# dataset is exactly 121, i.e. that record touched the ceiling.
#
# NB this is the ceiling for fuel model Y, which is what all 149 NFDRSInit files
# actually set (fuelModel = "Y"). The customFuelModel block in those files is
# parameterised from NFDRS88 fuel model G but is inert unless fuelModel = "C".
# Fuel model Y confirmed 2026-08-26, so 121 is the operative ceiling. The G-
# parameterised block is a red herring for a future reader - ignore it.
#
# Truncation fixes physical realism only. It cannot change any above-P80 run
# length: P80 sits far below the ceiling, so a clipped day is still above P80.
erc.max <- 121

output.cs <- matrix(NA_real_, n.sim, n.valid)
n.clip <- 0L
for (s in 1:n.valid) {
  e <- rep(nonzero.mean.gm[, s], nyear.sim) + rep(nonzero.sd.gm[, s], nyear.sim) * xs.cs[, s]
  e[e < 0] <- 0
  n.clip <- n.clip + sum(e > erc.max)
  e[e > erc.max] <- erc.max
  output.cs[, s] <- e
}
cat(sprintf("ERC truncated at %d: %d of %d values clipped (%.5f%%)\n",
            erc.max, n.clip, n.sim * n.valid, 100 * n.clip / (n.sim * n.valid)))
outmat.cs <- round(output.cs[-(1:365), , drop = FALSE])
colnames(outmat.cs)     <- namevec[valid.gm]
storage.mode(outmat.cs) <- "integer"
rm(output.cs); gc(verbose = FALSE)
write.table(outmat.cs, "erc-sim-regime-gridmet-int.csv", sep = ",", row.names = FALSE)

p.cs  <- vapply(1:n.valid, function(s) mean(outmat.cs[, s] >= thr.p80[s]), 0)
cs.sd <- vapply(1:n.valid, function(s) sd(xs.cs[, s]), 0) / sdresids.valid.gm
cat(sprintf("sd ratio: mean=%.3f  exceedance rate=%.4f (observed %.4f)\n",
            mean(cs.sd), mean(p.cs), mean(p.ob)))
cat(sprintf("Correlation-matrix RMSE vs observed: regime=%.4f\n",
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
cat(sprintf("above-P80 |fall|/rise: observed=%.3f  regime=%.3f\n",
            rf_ratio(obsmat.valid, thr.p80), rf_ratio(outmat.cs, thr.p80)))


###--- Burn block analysis: observed vs regime-switching model ---###
burn_blocks <- function(erc_vec, ref_vec = NULL, p_thresh = 0.80) {
  ref   <- na.exclude(if (is.null(ref_vec)) erc_vec else ref_vec)
  if (length(ref) == 0) return(list(above = integer(0), below = integer(0)))
  above <- ecdf(ref)(erc_vec) >= p_thresh
  runs  <- rle(as.vector(na.exclude(above)))
  list(above = runs$lengths[ runs$values],
       below = runs$lengths[!runs$values])
}

obs_blocks <- lapply(valid.station.inds, function(s) burn_blocks(ercmat.gm[, s]))
cs_blocks  <- lapply(seq_along(valid.station.inds), function(s) {
  burn_blocks(outmat.cs[, s], ref_vec = ercmat.gm[, valid.station.inds[s]])
})

obs_above <- unlist(lapply(obs_blocks, `[[`, "above"))
obs_below <- unlist(lapply(obs_blocks, `[[`, "below"))
cs_above  <- unlist(lapply(cs_blocks,  `[[`, "above"))
cs_below  <- unlist(lapply(cs_blocks,  `[[`, "below"))

cat("\nObserved above-P80 block lengths:\n"); print(summary(obs_above))
cat("Regime   above-P80 block lengths:\n");   print(summary(cs_above))
cat("Observed below-P80 block lengths:\n");   print(summary(obs_below))
cat("Regime   below-P80 block lengths:\n");   print(summary(cs_below))

cat(sprintf("\n%-22s %10s %10s\n", "above-P80 statistic", "observed", "regime"))
for (r in list(c("mean log run length", "meanlog"), c("mean run length", "mean"),
               c("median run length", "median"),    c("% episodes >=14d", "f14"),
               c("% episodes >=30d", "f30"))) {
  f <- switch(r[2], meanlog = function(v) mean(log(v)), mean = mean,
              median = function(v) median(v),
              f14 = function(v) mean(v >= 14) * 100, f30 = function(v) mean(v >= 30) * 100)
  cat(sprintf("%-22s %10.3f %10.3f\n", r[1], f(obs_above), f(cs_above)))
}

setwd(output_dir)
grp.cols <- c("steelblue", "tomato")

png("burn-blocks-p80-regime.png", width = 900, height = 700)
par(mfrow = c(2, 2))
hist(obs_above, breaks = 30, main = "Observed: runs above P80", xlab = "Block length (days)",
     col = grp.cols[1], border = "white", xlim = c(0, 150))
hist(cs_above,  breaks = 30, main = "Regime: runs above P80",   xlab = "Block length (days)",
     col = grp.cols[2], border = "white", xlim = c(0, 150))
hist(obs_below, breaks = 30, main = "Observed: runs below P80", xlab = "Block length (days)",
     col = grp.cols[1], border = "white", xlim = c(0, 500))
hist(cs_below,  breaks = 30, main = "Regime: runs below P80",   xlab = "Block length (days)",
     col = grp.cols[2], border = "white", xlim = c(0, 500))
dev.off()

png("burn-block-boxplot-regime.png", width = 700, height = 450)
par(mfrow = c(1, 2))
boxplot(c(log(obs_above), log(cs_above)) ~
        factor(rep(c("observed","regime"), c(length(obs_above), length(cs_above))),
               levels = c("observed","regime")),
        main = "Above-P80", xlab = "", ylab = "log run length (days)", col = grp.cols)
boxplot(c(log(obs_below), log(cs_below)) ~
        factor(rep(c("observed","regime"), c(length(obs_below), length(cs_below))),
               levels = c("observed","regime")),
        main = "Below-P80", xlab = "", ylab = "log run length (days)", col = grp.cols)
dev.off()

png("burn-block-ecdf-regime.png", width = 700, height = 450)
plot(ecdf(obs_above), main = "ECDF of above-P80 block lengths",
     xlab = "Block length (days)", ylab = "Cumulative probability",
     xlim = c(0, 60), col = grp.cols[1])
plot(ecdf(cs_above), add = TRUE, col = grp.cols[2])
legend("bottomright", c("observed", "regime"), col = grp.cols, lty = 1, bty = "n")
dev.off()

###--- Hazard of run termination ---###
# Localises the remaining error: the model still ends runs too readily on day 1,
# converging to the observed hazard by about day 6. NB the pooled curve overstates
# real duration dependence - about half its decline is a pooling artefact across
# heterogeneous pyromes, so validate against a stratified hazard if it matters.
hazard  <- function(v, K = 30) vapply(1:K, function(k) sum(v == k) / sum(v >= k), 0)
haz.mat <- cbind(observed = hazard(obs_above), regime = hazard(cs_above))
png("burn-block-hazard-regime.png", width = 700, height = 450)
matplot(1:30, haz.mat, type = "o", pch = 16, cex = 0.6, lty = 1, col = grp.cols,
        xlab = "Run length k (days)", ylab = "P(run ends at k | run >= k)",
        main = "Hazard of above-P80 run termination")
legend("topright", colnames(haz.mat), col = grp.cols, lty = 1, pch = 16, bty = "n")
dev.off()
cat("\nHazard of above-P80 run termination:\n")
print(round(haz.mat[c(1, 2, 3, 5, 7, 10, 14, 20, 30), ], 4))

###--- Acceptance test: is the observed record a plausible 15-yr draw? ---###
# A t-test on pooled episodes cannot answer this: it compares ~32k observed
# episodes against millions of simulated ones and treats every episode as
# independent, when they are correlated both serially within a pyrome and
# spatially across pyromes. Any trivial difference then returns p < 2e-16.
# The correct reference is the model's OWN 15-year sampling distribution, which
# preserves the spatial dependence, the serial dependence, and the run
# truncation a short record imposes.
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
report_acceptance(acceptance_null(outmat.cs), obs.win, "Regime-switching model:")
