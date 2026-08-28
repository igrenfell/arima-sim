library(MASS)
library(stats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(lmtest)
library(stringr)
library(lubridate)
library(ggplot2)
library(forecast)
library(sf)
library(pROC)

###read nfdrs data

#setwd("G:\\Workspace\\ercm")
#stationcat <- read.csv("virtualStnCatalog_update.csv")
#statnum <- stationcat[,2]
#tatnum.0 <- str_pad(statnum, 6, pad = "0")
#stationcat[,2] <- statnum.0
#write.table(stationcat, "virtualStnCatalog_update_zero.csv", row.names = FALSE, sep = ",", quote= F)
#setwd("G:\\Workspace\\ercm\\NFDRS_output")
###--- Configuration: 2026 pyrome scheme ---###
# 2026 pyromes are numbered globally across three geodatabases:
#   CONUS_v2.gdb  CONUS_Pyromes_Expanded  141 features  PYROME   1-141
#   AK_v2.gdb     AK_Pyromes_Expand         7 features  PYROME 142-148
#   HI_v2.gdb     HI_Pyrome                 1 feature   PYROME     149
# = 149, matching the 149 NFDRS files. The filename carries that same number:
# PYROME_0001NFDRSoutputMinMaxMed.csv -> PYROME 1. Column order below is set by
# pyrome number, not by glob order, and every column is labelled with it.
#
# NB the numbering CHANGED from the previous 128-pyrome scheme. Results are not
# comparable pyrome-for-pyrome; pyinfo$ORIG_PYROME_NO carries the crosswalk back
# to the old numbers, but it is many-to-one (14 of the new CONUS pyromes are
# splits of old ones) so old results must be aggregated, not renamed.

nfdrs_dir  <- "I:/workspace/ercm/ERC_2026/From_Matt_20260812/NFDRS4HourlyAndDailyOut"
pyrome_dir <- "I:/workspace/ercm/Pyromes-2026/2.1_Pyromes"
output_dir <- "I:/workspace/ercm/output-sim/fsim-generator"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
conus.only <- TRUE     # TRUE restricts the run to PYROME 1-141

rd <- function(gdb, lyr, region, orig) {
  d <- st_drop_geometry(st_read(file.path(pyrome_dir, gdb), layer = lyr, quiet = TRUE))
  data.frame(PYROME = as.integer(d$PYROME), REGION = region,
             NAME   = if ("NAME" %in% names(d)) d$NAME else NA_character_,
             ORIG_PYROME_NO = as.integer(d[[orig]]))
}
pyinfo <- rbind(rd("CONUS_v2.gdb", "CONUS_Pyromes_Expanded", "CONUS", "ORIG_PYROME_NO"),
                rd("AK_v2.gdb",    "AK_Pyromes_Expand",      "AK",    "ORIG_PYROME"),
                rd("HI_v2.gdb",    "HI_Pyrome",              "HI",    "ORIG_PYROME"))
pyinfo <- pyinfo[order(pyinfo$PYROME), ]
cat("pyrome lookup:", nrow(pyinfo), "total -",
    sum(pyinfo$REGION == "CONUS"), "CONUS,", sum(pyinfo$REGION == "AK"), "AK,",
    sum(pyinfo$REGION == "HI"), "HI\n")

setwd(nfdrs_dir)
flist <- Sys.glob("*.csv")
flist <- flist[grep("NFDRS", flist)]
flist <- flist[grep("MinMax", flist)]

pynum <- as.integer(sub("^PYROME_0*([0-9]+)NFDRS.*$", "\\1", flist))
if (anyNA(pynum))
  stop("could not parse a pyrome number from: ",
       paste(flist[is.na(pynum)], collapse = ", "))
ord   <- order(pynum)
flist <- flist[ord]; pynum <- pynum[ord]

if (conus.only) {
  keep  <- pynum %in% pyinfo$PYROME[pyinfo$REGION == "CONUS"]
  flist <- flist[keep]; pynum <- pynum[keep]
}
want <- if (conus.only) pyinfo$PYROME[pyinfo$REGION == "CONUS"] else pyinfo$PYROME

if (any(duplicated(pynum)))
  stop("duplicate pyrome numbers in file list: ",
       paste(unique(pynum[duplicated(pynum)]), collapse = ", "))
if (length(setdiff(pynum, pyinfo$PYROME)))
  stop("NFDRS file(s) with no pyrome polygon: ",
       paste(setdiff(pynum, pyinfo$PYROME), collapse = ", "))
if (length(setdiff(want, pynum)))
  stop("pyrome(s) with no NFDRS file: ", paste(setdiff(want, pynum), collapse = ", "))

nstations <- length(flist)
pyinfo    <- pyinfo[match(pynum, pyinfo$PYROME), ]   # aligned to column order
cat("stations:", nstations, " pyromes", min(pynum), "-", max(pynum), "\n")

# Two windows. The NFDRS files hold yearseq.data; every model below is fitted on
# yearseq alone. Keeping the whole record loaded costs nothing extra on the read
# and lets the percentile table contrast the fitting window against the full
# weather history.
yearseq.data <- 2003:2025          # everything the NFDRS files contain
yearseq      <- 2011:2025          # model-fitting window
stopifnot(all(yearseq %in% yearseq.data))

ndays.data <- 365*length(yearseq.data)
ndays.all  <- 365*length(yearseq)

ercmat.full <- matrix(NA, nrow = ndays.data, ncol = nstations)
namevec <- paste0("py", pynum)      # column identity = new PYROME number

dayseq.data      <- rep(1:365, length(yearseq.data))
yearseq.data.all <- rep(yearseq.data, each = 365)
dayseq.all       <- rep(1:365, length(yearseq))
yearseq.all      <- rep(yearseq, each = 365)
cat("data window:", min(yearseq.data), "-", max(yearseq.data),
    " | model window:", min(yearseq), "-", max(yearseq), "\n")

for(curstation in 1:nstations)
{
  fname  <- flist[curstation]
  curdat <- read.csv(fname)

  ##--- fail fast rather than silently poisoning the climatology and AR fits ---##
  # NB: `$` partial-matches on data frames, so curdat$Date would silently return
  # the DateTime column. Reference the real column name.
  if (!all(c("DateTime", "ERC") %in% names(curdat)))
    stop("[", fname, "] missing required column(s): ",
         paste(setdiff(c("DateTime", "ERC"), names(curdat)), collapse = ", "))

  # Timestamps are ISO basic UTC, e.g. 20030101T230000Z. Take the date portion by
  # position instead of building a datetime: the stamp is 23:00Z, so any timezone
  # conversion can roll the calendar day.
  datevals <- curdat$DateTime
  datev    <- as.Date(substr(datevals, 1, 8), format = "%Y%m%d")
  if (anyNA(datev))
    stop("[", fname, "] ", sum(is.na(datev)), " DateTime value(s) failed to parse; first: '",
         datevals[which(is.na(datev))[1]], "'")

  yearn    <- as.numeric(format(datev, "%Y"))
  subyears <- yearn %in% yearseq.data
  curdat   <- curdat[subyears, ]
  doy      <- format(datev[subyears], "%j")

  #ditch leap days (drops Dec 31 of leap years, leaving 365 days/year)
  curdat.sub <- curdat[doy != "366", ]

  if (nrow(curdat.sub) != ndays.data)
    stop("[", fname, "] has ", nrow(curdat.sub), " rows after restricting to ",
         min(yearseq.data), "-", max(yearseq.data),
         " and dropping day 366; expected ", ndays.data)

  ercvec <- curdat.sub$ERC
  if (!is.numeric(ercvec))
    stop("[", fname, "] ERC did not read as numeric (class '", class(ercvec),
         "'); check for text or malformed values in the column")
  if (anyNA(ercvec))
    stop("[", fname, "] ERC contains ", sum(is.na(ercvec)),
         " NA value(s); first at row ", which(is.na(ercvec))[1])
  if (any(ercvec < 0))
    stop("[", fname, "] ERC contains ", sum(ercvec < 0),
         " negative value(s) (missing-data sentinel or corruption); first is ",
         ercvec[which(ercvec < 0)[1]], " at row ", which(ercvec < 0)[1])

  # the file's own StationID must agree with the pyrome number in its filename,
  # otherwise a column would be silently mislabelled
  if ("StationID" %in% names(curdat.sub)) {
    sid <- as.integer(sub("^PYROME_0*([0-9]+)$", "\\1", curdat.sub$StationID[1]))
    if (!is.na(sid) && sid != pynum[curstation])
      stop("[", fname, "] StationID '", curdat.sub$StationID[1],
           "' disagrees with pyrome ", pynum[curstation], " parsed from the filename")
  }

  ercmat.full[, curstation] <- ercvec

  print(curstation/nstations)
  
  
  
}

# everything downstream models the fitting window; ercmat.full is kept only for
# the percentile comparison at the end
ercmat      <- ercmat.full[yearseq.data.all %in% yearseq, , drop = FALSE]
ercmat.base <- ercmat
stopifnot(nrow(ercmat) == ndays.all)
cat("ercmat.full:", paste(dim(ercmat.full), collapse = " x "),
    " -> model window ercmat:", paste(dim(ercmat), collapse = " x "), "\n")

###--- FOD: per-pyrome large-fire occurrence vs ERC ---###
# Fits, for every pyrome, a logistic model of "was there a large fire today"
# against that day's ERC, and records the AUC. This is the validation that ERC
# actually discriminates large-fire days for the new pyrome geography.
#
# Two things differ from the old 128-pyrome version of this code:
#  1. The 2026 FOD carries DISCOVERY_DATE as a proper ISO date, so the old
#     strsplit / "%m/%d/%Y" parsing is gone, and PYROME is already attributed to
#     the new 149-pyrome scheme (see FOD_Pyromes_2026_merge.csv).
#  2. The large-fire thresholds are still keyed to the OLD numbering (1-137), so
#     they are mapped across with pyrome-crosswalk-2026-to-old.csv. That mapping
#     is many-to-one: the 22 new pyromes that came from a split inherit their
#     parent's threshold, which is an approximation, not a per-pyrome refit.

fod_file  <- "I:/workspace/fod-2026/FOD_Pyromes_2026_merge.csv"
lft_file  <- file.path(nfdrs_dir, "..", "large-fire-threshold-lorenz.csv")
xw_file   <- "I:/workspace/ercm/Pyromes-2026/pyrome-crosswalk-2026-to-old.csv"
roc_pdf   <- TRUE       # write per-pyrome ROC curves to a PDF

## --- large-fire threshold per NEW pyrome, via the old-number crosswalk ------
if (!file.exists(lft_file))
  lft_file <- "I:/workspace/fsim-climate/station-data-complete/PsuedoStationsPlusERC/BASE_2020/large-fire-threshold-lorenz.csv"
lft.raw <- read.csv(lft_file)
xw      <- read.csv(xw_file)
lft.map <- merge(xw[, c("PYROME", "OLD_PYROME")], lft.raw,
                 by.x = "OLD_PYROME", by.y = "Pyrome", all.x = TRUE)
lft.vec <- lft.map$LFT[match(pynum, lft.map$PYROME)]
if (anyNA(lft.vec))
  stop("no large-fire threshold for pyrome(s): ",
       paste(pynum[is.na(lft.vec)], collapse = ", "))
cat("large-fire thresholds mapped for", length(lft.vec), "pyromes  (range ",
    min(lft.vec), "-", max(lft.vec), " acres)\n", sep = " ")

## --- the ERC calendar: 365 days per year, dropping Dec 31 of leap years ----
# seq(Jan 1, length.out = 365) reproduces exactly that convention, so the fire
# day-series lines up row-for-row with ercmat.
all.dates <- as.Date(unlist(lapply(yearseq, function(y)
  as.character(seq(as.Date(paste0(y, "-01-01")), by = "day", length.out = 365)))))
stopifnot(length(all.dates) == ndays.all)

## --- FOD, restricted to the modelled years --------------------------------
fod <- read.csv(fod_file)
fod$DISC <- as.Date(fod$DISCOVERY_DATE)
if (anyNA(fod$DISC))
  stop(sum(is.na(fod$DISC)), " DISCOVERY_DATE value(s) failed to parse; first: '",
       fod$DISCOVERY_DATE[which(is.na(fod$DISC))[1]], "'")
fod <- fod[!is.na(fod$PYROME) & fod$FIRE_YEAR %in% yearseq, ]
cat("FOD rows in", min(yearseq), "-", max(yearseq), "with a pyrome:", nrow(fod), "\n")

# The FOD ends before the weather record does. Any modelled year with no FOD
# coverage would otherwise contribute 365 spurious "no large fire" days and bias
# the logistic fits, so those days are excluded from the fit rather than counted
# as zeros.
fod.years <- sort(unique(fod$FIRE_YEAR))
fit.keep  <- yearseq.all %in% fod.years
if (!all(yearseq %in% fod.years))
  cat("NOTE: no FOD coverage for", paste(setdiff(yearseq, fod.years), collapse = ", "),
      "- excluding", sum(!fit.keep), "of", ndays.all, "days from the logistic fits\n")

## --- per-pyrome logistic model --------------------------------------------
auc.vec        <- rep(NA_real_, nstations)
lrt.p          <- rep(NA_real_, nstations)
erc.coef       <- rep(NA_real_, nstations)
erc.intercept  <- rep(NA_real_, nstations)

nfire.days     <- rep(NA_integer_, nstations)
nfires.used    <- rep(NA_integer_, nstations)
logistic.list  <- vector("list", nstations)

if (roc_pdf) pdf(file.path(output_dir, "roc-curves-pyromes-2026.pdf"), width = 7, height = 7)
for (i in seq_len(nstations)) {
  py  <- pynum[i]
  erc <- ercmat.base[, i]

  fs  <- fod[fod$PYROME == py & fod$FIRE_SIZE >= lft.vec[i], ]
  nfires.used[i] <- nrow(fs)

  # map fire dates onto the ERC calendar; NA = a date the calendar omits
  # (Dec 31 of a leap year), which is dropped rather than misassigned
  idx <- match(fs$DISC, all.dates)
  isfday <- integer(ndays.all)
  isfday[idx[!is.na(idx)]] <- 1L
  nfire.days[i] <- sum(isfday)

  # need both classes present, and enough events for the fit to mean anything
  if (nfire.days[i] < 10 || nfire.days[i] == ndays.all) next

  yy <- isfday[fit.keep]; xx <- erc[fit.keep]
  fitmod  <- glm(yy ~ xx, family = "binomial")
  nullmod <- glm(yy ~ 1,  family = "binomial")
  logistic.list[[i]] <- summary(fitmod)
  erc.coef[i] <- coef(fitmod)[["xx"]]
  erc.intercept[i] <- coef(fitmod)[["(Intercept)"]]

  lrt.p[i] <- pchisq(-2 * (as.numeric(logLik(nullmod)) - as.numeric(logLik(fitmod))),
                     df = 1, lower.tail = FALSE)

  pr <- predict(fitmod, type = "response")   # aligned to fit.keep days
  # NB no bootstrap CI here: the old code passed boot.n = 1000 / ci = TRUE but
  # only ever kept rocgm$auc, so the bootstrap cost 149x for nothing.
  rc <- pROC::roc(yy, pr, quiet = TRUE, plot = roc_pdf, legacy.axes = TRUE,
                  print.auc = TRUE, grid = TRUE,
                  main = sprintf("pyrome %d (%s)  n large fires = %d",
                                 py, pyinfo$REGION[i], nfires.used[i]))
  auc.vec[i] <- as.numeric(rc$auc)
}
if (roc_pdf) dev.off()

## --- summary --------------------------------------------------------------
fit.ok <- !is.na(auc.vec)
cat("\nlogistic ERC -> large-fire-day models fitted:", sum(fit.ok), "of", nstations, "\n")
cat("skipped (fewer than 10 large-fire days):",
    paste(pynum[!fit.ok], collapse = ", "), "\n")
cat("\nAUC across pyromes:\n"); print(summary(auc.vec[fit.ok]))
cat("AUC < 0.6 (weak discrimination):", sum(auc.vec[fit.ok] < 0.6), "pyromes\n")
cat("positive ERC coefficient:", sum(erc.coef[fit.ok] > 0), "of", sum(fit.ok), "\n")
cat("likelihood-ratio p < 0.05:", sum(lrt.p[fit.ok] < 0.05), "of", sum(fit.ok), "\n")

logistic.summary <- data.frame(
  PYROME      = pynum,
  REGION      = pyinfo$REGION,
  NAME        = pyinfo$NAME,
  OLD_PYROME  = pyinfo$ORIG_PYROME_NO,
  LFT_acres   = lft.vec,
  n_fires     = nfires.used,
  n_fire_days = nfire.days,
  erc_coef    = erc.coef,
  erc_intercept  =erc.intercept,
  lrt_p       = lrt.p,
  AUC         = auc.vec)
write.csv(logistic.summary, file.path(output_dir, "logistic-erc-firedays-2026.csv"),
          row.names = FALSE)
cat("\nwrote logistic-erc-firedays-2026.csv\n")
print(utils::head(logistic.summary[order(-logistic.summary$AUC),
                                   c("PYROME","REGION","LFT_acres","n_fire_days","AUC")], 5))

###--- Preliminary: fire counts by pyrome ---###
# Counts are over the modelling window (yearseq) and over the pyromes actually
# being run, so they respond to conus.only. "large" means FIRE_SIZE >= that
# pyrome's Lorenz large-fire threshold.
fires.by.pyrome <- do.call(rbind, lapply(seq_len(nstations), function(i) {
  py <- pynum[i]
  f  <- fod[fod$PYROME == py, ]
  fl <- f[f$FIRE_SIZE >= lft.vec[i], ]
  sz <- f$FIRE_SIZE[!is.na(f$FIRE_SIZE)]
  data.frame(
    PYROME        = py,
    REGION        = pyinfo$REGION[i],
    NAME          = pyinfo$NAME[i],
    LFT_acres     = lft.vec[i],
    n_fires_all   = nrow(f),
    n_fires_large = nrow(fl),
    fires_per_yr  = round(nrow(f)  / length(yearseq), 1),
    large_per_yr  = round(nrow(fl) / length(yearseq), 2),
    acres_total   = if (length(sz)) round(sum(sz))    else NA_real_,
    acres_median  = if (length(sz)) median(sz)        else NA_real_,
    acres_max     = if (length(sz)) round(max(sz))    else NA_real_,
    n_fire_days   = nfire.days[i],
    stringsAsFactors = FALSE)
}))

write.csv(fires.by.pyrome, file.path(output_dir, "fires-by-pyrome-2026.csv"),
          row.names = FALSE)

cat("\n=== fires by pyrome,", min(yearseq), "-", max(yearseq), "===\n")
cat("pyromes:", nrow(fires.by.pyrome),
    "  total fires:", sum(fires.by.pyrome$n_fires_all),
    "  large fires:", sum(fires.by.pyrome$n_fires_large), "\n")
cat("\nall fires per pyrome:\n");   print(summary(fires.by.pyrome$n_fires_all))
cat("large fires per pyrome:\n");   print(summary(fires.by.pyrome$n_fires_large))
cat("\npyromes with zero fires:",
    paste(fires.by.pyrome$PYROME[fires.by.pyrome$n_fires_all == 0], collapse = ", "), "\n")
cat("pyromes with < 10 large fires:",
    paste(fires.by.pyrome$PYROME[fires.by.pyrome$n_fires_large < 10], collapse = ", "), "\n")

cat("\n-- 10 busiest pyromes by total fire count --\n")
print(fires.by.pyrome[order(-fires.by.pyrome$n_fires_all),
      c("PYROME","NAME","n_fires_all","fires_per_yr","n_fires_large","acres_total")][1:10, ],
      row.names = FALSE)
cat("\n-- 10 sparsest pyromes by total fire count --\n")
print(fires.by.pyrome[order(fires.by.pyrome$n_fires_all),
      c("PYROME","NAME","n_fires_all","fires_per_yr","n_fires_large","acres_total")][1:10, ],
      row.names = FALSE)
cat("\nwrote fires-by-pyrome-2026.csv\n")

###--- Bridge: hand the new-pyrome data to the regime pipeline ------------###
# Everything below is the regime-switching simulator, lifted unchanged from
# erc-sim-regime.R. It expects the ERC matrix under the name ercmat.gm and
# builds its own standardised residuals, AR fits, EOF and t-distributed PCs.
#
# NOTE this REPLACES the old mvrnorm-on-covariance simulator that
# erc-generator-fsim.r uses. Two substantive differences:
#   * residuals are standardised by the daily SD, not just centred on the daily
#     mean, and the reconstruction multiplies dailysd back in. The old generator
#     omitted that and ran its P80 about 2 ERC low, worst in the high-variance
#     southwestern pyromes.
#   * output is truncated at the fuel-model-Y ceiling (erc.max = 121).
ercmat.gm <- ercmat            # model window only (yearseq), not the full record

# the lifted climatology block writes these; keep the names it expects
py.test      <- pynum
dailymean.gm <- matrix(NA, nrow = 365, ncol = nstations)
dailysd.gm   <- matrix(NA, nrow = 365, ncol = nstations)

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
nyear.sim <- 10000
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

obsmat.valid <- ercmat.gm[,   valid.station.inds, drop = FALSE]   # fitting window
obsmat.full  <- ercmat.full[, valid.station.inds, drop = FALSE]   # whole record
cor.obs      <- cor(obsmat.valid)
resids.valid <- residsmat.gm[, valid.gm, drop = FALSE]

# P80 in ERC space per pyrome, taken over the FULL historical record
# (yearseq.data), not the fitting window. This single value is the threshold
# everywhere downstream: it sets the above/below regime state that drives the
# simulation, and it is the burn-block threshold for both observed and simulated
# streams. Using the full record means the threshold does not move when the
# fitting window changes, and it is estimated from 8395 days rather than 5475.
thr.p80 <- vapply(valid.station.inds, function(s) {
  v  <- na.exclude(ercmat.full[, s]); su <- sort(unique(v))
  su[min(which(vapply(su, function(x) mean(v <= x), 0) >= 0.80))]
}, numeric(1))

# on the standardised scale the state threshold varies with day-of-year
q.obs <- matrix(NA_real_, 365, n.valid)
for (s in 1:n.valid)
  q.obs[, s] <- (thr.p80[s] - nonzero.mean.gm[, s]) / nonzero.sd.gm[, s]

q.sim <- q.obs[rep(1:365, nyear.sim), , drop = FALSE]
guard <- 8 * sdresids.valid.gm      # stability clamp on the feedback loop
p.ob   <- vapply(1:n.valid, function(s) mean(obsmat.valid[, s] >= thr.p80[s]), 0)
p.ob.f <- vapply(1:n.valid, function(s) mean(obsmat.full[,  s] >= thr.p80[s]), 0)
cat(sprintf("P80 threshold from %d-%d: mean %.2f ERC\n",
            min(yearseq.data), max(yearseq.data), mean(thr.p80)))
cat(sprintf("exceedance at that threshold: full record %.4f | fitting window %.4f\n",
            mean(p.ob.f), mean(p.ob)))


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
            rf_ratio(obsmat.full, thr.p80), rf_ratio(outmat.cs, thr.p80)))


###--- Burn block analysis: observed vs regime-switching model ---###
# Explicit threshold rather than a per-call ecdf, so there is exactly one
# definition of "above P80" in this script: thr.p80, from the full record.
burn_blocks <- function(erc_vec, thr) {
  above <- na.exclude(as.vector(erc_vec >= thr))
  if (!length(above)) return(list(above = integer(0), below = integer(0)))
  runs <- rle(above)
  list(above = runs$lengths[ runs$values],
       below = runs$lengths[!runs$values])
}

# observed baseline uses the FULL record, matching the threshold's basis
obs_blocks <- lapply(seq_along(valid.station.inds), function(s)
  burn_blocks(obsmat.full[, s], thr.p80[s]))
cs_blocks  <- lapply(seq_along(valid.station.inds), function(s)
  burn_blocks(outmat.cs[, s],   thr.p80[s]))
obs_blocks.fit <- lapply(seq_along(valid.station.inds), function(s)
  burn_blocks(obsmat.valid[, s], thr.p80[s]))

obs_above <- unlist(lapply(obs_blocks, `[[`, "above"))
obs_below <- unlist(lapply(obs_blocks, `[[`, "below"))
cs_above  <- unlist(lapply(cs_blocks,  `[[`, "above"))
cs_below  <- unlist(lapply(cs_blocks,  `[[`, "below"))

cat("\nObserved above-P80 block lengths:\n"); print(summary(obs_above))
cat("Regime   above-P80 block lengths:\n");   print(summary(cs_above))
cat("Observed below-P80 block lengths:\n");   print(summary(obs_below))
cat("Regime   below-P80 block lengths:\n");   print(summary(cs_below))

obs_above.fit <- unlist(lapply(obs_blocks.fit, `[[`, "above"))
cat(sprintf("\n%-22s %10s %10s %10s\n", "above-P80 statistic",
            "obs(full)", "obs(fit)", "regime"))
for (r in list(c("mean log run length", "meanlog"), c("mean run length", "mean"),
               c("median run length", "median"),    c("% episodes >=14d", "f14"),
               c("% episodes >=30d", "f30"))) {
  f <- switch(r[2], meanlog = function(v) mean(log(v)), mean = mean,
              median = function(v) median(v),
              f14 = function(v) mean(v >= 14) * 100, f30 = function(v) mean(v >= 30) * 100)
  cat(sprintf("%-22s %10.3f %10.3f %10.3f\n", r[1],
              f(obs_above), f(obs_above.fit), f(cs_above)))
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

###--- Acceptance test: is the observed record a plausible draw? ---###
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

# window length follows the observed record being compared (the full record),
# so the run-truncation a finite record imposes is matched
acceptance_null <- function(mat, nyr.win = length(yearseq.data), nwin = 200,
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
    cat(sprintf("  %-7s observed=%9.3f | model %d-yr window: mean=%9.3f  95%%CI=[%9.3f,%9.3f]  observed pctile=%5.1f%%\n",
                v, o, length(yearseq.data), mean(d), quantile(d, 0.025),
                quantile(d, 0.975), mean(d <= o) * 100))
  }
}

obs.win <- window_stat(obsmat.full, thr.p80)
report_acceptance(acceptance_null(outmat.cs), obs.win, "Regime-switching model:")

###--- ERC percentiles: full record vs fitting window vs regime simulation ---###
# full = observed ERC over the whole NFDRS record (yearseq.data)
# fit  = observed ERC over the model-fitting window (yearseq)
# sim  = regime-switching simulation
# The sim-vs-fit gap is the validation that matters: the simulator is built from
# the fitting window and should reproduce its distribution. The old
# mvrnorm generator ran P80 about 2 ERC low here because it never scaled by the
# daily SD; this reconstruction does, so that bias should be gone.
pctl <- c(0.50, 0.80, 0.90, 0.95, 0.99)
stopifnot(sum(valid.gm) == ncol(outmat.cs))

qcols <- function(m, tag) {
  z <- t(apply(m, 2, quantile, probs = pctl, na.rm = TRUE))
  colnames(z) <- paste0(tag, "_p", sub("^0\\.", "", format(pctl)))
  as.data.frame(z, row.names = NULL)
}

erc.pctl <- cbind(
  data.frame(PYROME = pynum[valid.gm],
             REGION = pyinfo$REGION[valid.gm],
             NAME   = pyinfo$NAME[valid.gm],
             row.names = NULL),
  qcols(ercmat.full[, valid.gm, drop = FALSE], "full"),
  qcols(ercmat.gm[,   valid.gm, drop = FALSE], "fit"),
  qcols(outmat.cs,                             "sim"))

erc.pctl$d_fit_full <- round(erc.pctl$fit_p80 - erc.pctl$full_p80, 2)
erc.pctl$d_sim_fit  <- round(erc.pctl$sim_p80 - erc.pctl$fit_p80,  2)

write.csv(erc.pctl, file.path(output_dir, "erc-percentiles-by-pyrome-regime.csv"),
          row.names = FALSE)

cat("\n=== ERC 80th percentile by pyrome (regime model) ===\n")
cat("full record:", min(yearseq.data), "-", max(yearseq.data),
    " | fitting window:", min(yearseq), "-", max(yearseq),
    " | simulation:", nyear.sim, "years\n\n")
print(rbind(full = summary(erc.pctl$full_p80),
            fit  = summary(erc.pctl$fit_p80),
            sim  = summary(erc.pctl$sim_p80)))
cat(sprintf("\nmean P80: full=%.2f  fit=%.2f  sim=%.2f\n",
            mean(erc.pctl$full_p80), mean(erc.pctl$fit_p80), mean(erc.pctl$sim_p80)))
cat(sprintf("fit - full: mean %+.2f\n", mean(erc.pctl$d_fit_full)))
cat(sprintf("sim - fit : mean %+.2f  median %+.2f  |  within +/-2 ERC: %d of %d\n",
            mean(erc.pctl$d_sim_fit), median(erc.pctl$d_sim_fit),
            sum(abs(erc.pctl$d_sim_fit) <= 2), nrow(erc.pctl)))
cat("correlation of P80 (sim vs fit):",
    round(cor(erc.pctl$sim_p80, erc.pctl$fit_p80), 4), "\n")

cat("\n-- 10 pyromes where the simulated P80 departs most from the fit window --\n")
print(erc.pctl[order(-abs(erc.pctl$d_sim_fit)),
      c("PYROME","NAME","full_p80","fit_p80","sim_p80","d_sim_fit")][1:10, ],
      row.names = FALSE)

png(file.path(output_dir, "erc-p80-fit-vs-sim-regime.png"), width = 900, height = 450)
par(mfrow = c(1, 2))
plot(erc.pctl$full_p80, erc.pctl$fit_p80, pch = 16, cex = 0.7, col = "steelblue",
     xlab = "P80, full record", ylab = "P80, fitting window",
     main = "Fitting window vs full record"); abline(0, 1, col = "red")
plot(erc.pctl$fit_p80, erc.pctl$sim_p80, pch = 16, cex = 0.7, col = "tomato",
     xlab = "P80, fitting window", ylab = "P80, simulated",
     main = "Regime simulation vs fitting window"); abline(0, 1, col = "red")
dev.off()
cat("\nwrote erc-percentiles-by-pyrome-regime.csv and erc-p80-fit-vs-sim-regime.png\n")

###--- Visual check: simulated seasons vs the input record, 5 pyromes ---###
# One representative pyrome per region, chosen by polygon centroid rather than by
# name. Deliberately spans very different ERC regimes: a wet coastal PNW pyrome,
# Mediterranean southern California, the arid Great Basin, the humid Southeast,
# and the continental upper Midwest.
#
# One image PER PYROME rather than a stacked 5-panel figure, so each can be
# viewed at full size. Files are named seasons-<type>-py<N>-<region>.png.
show.py <- c("Pacific NW"    = 1,     # NW Washington                 (-123, 48)
             "Southern CA"   = 33,    # S. California Mtns / N. Baja  (-118, 34)
             "Great Basin"   = 44,    # Carbonate Woodland/Sagebrush  (-115, 39)
             "Southeast"     = 140,   # Florida Flatwoods/Okefenokee  ( -82, 29)
             "Upper Midwest" = 91)    # N. Minnesota Lakes & Forests  ( -93, 47)

# column lookups: j indexes the valid subset (outmat.cs, thr.p80),
#                 i indexes the full station set (ercmat.gm, ercmat.full)
jj <- match(show.py, pynum[valid.gm])
if (anyNA(jj))
  stop("pyrome(s) not in the valid simulated set: ",
       paste(show.py[is.na(jj)], collapse = ", "))
ii <- valid.station.inds[jj]

nyr.obs <- length(yearseq)
nyr.sm  <- nrow(outmat.cs) %/% 365
doy.obs <- rep(1:365, nyr.obs)
doy.sim <- rep(1:365, nyr.sm)
yr.obs  <- rep(seq_along(yearseq), each = 365)
yr.sim  <- rep(seq_len(nyr.sm),    each = 365)

n.show <- 2                      # years of each drawn on the overlay
set.seed(4)
pick.obs <- sort(sample(nyr.obs, n.show))
pick.sim <- sort(sample(nyr.sm,  n.show))
obs.thin <- adjustcolor("steelblue", alpha.f = 0.75)
sim.thin <- adjustcolor("tomato",    alpha.f = 0.75)
shade    <- adjustcolor("tomato",    alpha.f = 0.18)
mo.start <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
mo.lab   <- c("J","F","M","A","M","J","J","A","S","O","N","D")

# slug for filenames: "Pacific NW" -> "pacific-nw"
slug <- tolower(gsub("[^A-Za-z0-9]+", "-", names(show.py)))
fig.name <- function(type, k)
  file.path(output_dir, sprintf("seasons-%s-py%d-%s.png", type, show.py[k], slug[k]))

# Checkpoint the plot inputs so the figures can be re-drawn without repeating
# the ~25 min simulation. Re-plot from:
#   d <- readRDS(file.path(output_dir, "season-plot-inputs.rds"))
saveRDS(list(obs = ercmat.gm[, ii, drop = FALSE],
             sim = outmat.cs[, jj, drop = FALSE],
             thr = thr.p80[jj], show.py = show.py,
             yearseq = yearseq, nyr.sm = nyr.sm),
        file.path(output_dir, "season-plot-inputs.rds"))

## ---- Figure set 1: observed (solid) vs simulated (dashed) seasons -------
# Thin lines are individual years. The bold line is the OBSERVED seasonal
# average over all record years (not just the n.show drawn), so it is the
# climatology rather than a mean of the small sample shown.
for (k in seq_along(show.py)) {
  o <- ercmat.gm[, ii[k]]; s <- as.numeric(outmat.cs[, jj[k]])
  yl <- c(0, max(o, s) * 1.02); thr <- thr.p80[jj[k]]
  clim.o <- tapply(o, doy.obs, mean)
  png(fig.name("overlay", k), width = 1100, height = 480, res = 110)
  par(mar = c(3.2, 3.8, 2.8, 0.8), mgp = c(2.3, 0.7, 0))
  plot(NA, xlim = c(1, 365), ylim = yl, xaxt = "n", xlab = "", ylab = "ERC",
       main = sprintf("py%d  %s      P80 = %.0f", show.py[k],
                      names(show.py)[k], thr), cex.main = 1.05)
  axis(1, at = mo.start, labels = mo.lab)
  for (y in pick.obs) lines(1:365, o[yr.obs == y], col = obs.thin, lty = 1, lwd = 1.2)
  for (y in pick.sim) lines(1:365, s[yr.sim == y], col = sim.thin, lty = 2, lwd = 1.2)
  lines(1:365, clim.o, col = "steelblue4", lty = 1, lwd = 2.8)
  abline(h = thr, col = "grey25", lty = 3)
  legend("topleft", bty = "n", cex = 0.78, seg.len = 2.6,
         legend = c(sprintf("observed, %d yr", n.show),
                    sprintf("simulated, %d yr", n.show),
                    "observed seasonal mean"),
         col = c(obs.thin, sim.thin, "steelblue4"),
         lty = c(1, 2, 1), lwd = c(1.2, 1.2, 2.8))
  dev.off()
}

## ---- Figure set 2: one season each, above-threshold spells shaded ------
set.seed(11)
par(mfrow = c(1,1))
for (k in seq_along(show.py)) {
  o <- ercmat.gm[, ii[k]]; s <- as.numeric(outmat.cs[, jj[k]])
  yl <- c(0, max(o, s) * 1.02); thr <- thr.p80[jj[k]]
  oy <- sample(nyr.obs, 2); sy <- sample(pick.sim, 2)
  png(fig.name("example", k), width = 1100, height = 620, res = 110)
  par(mfrow = c(1, 1), mar = c(3.0, 3.8, 2.4, 0.8), mgp = c(2.1, 0.7, 0))
  srcvec <- c(1,1,2,2)
  
    src <- srcvec[i]
    v1 <- o[yr.obs == oy[1]]
    v2 <- o[yr.obs == oy[2]]
    v3 <-  s[yr.sim == sy[1]]
    v4 <-  s[yr.sim == sy[2]]

    #v <- if (src == 1) o[yr.obs == oy] else s[yr.sim == sy]
    plot(NA, xlim = c(1, 365), ylim = yl, xaxt = "n", xlab = "", ylab = "ERC",
         main = sprintf("py%d %s - %s", show.py[k], names(show.py)[k],
                        if (src == 1) paste("observed", yearseq[oy])
                        else paste("simulated yr", sy)), cex.main = 0.98)
    axis(1, at = mo.start, labels = mo.lab)
    r <- rle(v3 >= thr); e <- cumsum(r$lengths); b <- e - r$lengths + 1
    for (m in which(r$values))
      rect(b[m] - 0.5, yl[1], e[m] + 0.5, yl[2], col = shade, border = NA)
    #lines(1:365, v, col = if (src == 1) "steelblue" else "tomato",
    #      lty = if (src == 1) 1 else 2, lwd = 1.4)
      lines(1:365, v1, col = "blue",
          lty = if (src == 1) 1 else 2, lwd = 1.4)
      lines(1:365, v2, col = "steelblue",
          lty = if (src == 1) 1 else 2, lwd = 1.4)
      lines(1:365, v3, col = "tomato",
          lty = 2, lwd = 1.4)
        lines(1:365, v4, col = "red",
          lty = 2, lwd = 1.4)
  
  
    abline(h = thr, col = "grey25", lty = 3)
    #legend("topleft", bty = "n", cex = 0.75,
    #       legend = sprintf("%d d above P80, longest %d",
    #                        sum(v >= thr),
    #                        if (any(r$values)) max(r$lengths[r$values]) else 0))
    legend("topleft", legend = c("Obs 1", "Obs 2", "Sim 1", "Sim 2"), 
    lty = c(1, 1, 2, 2), col = c("blue", "steelblue", "tomato", "red")
  )
    dev.off()

}


## ---- Figure set 3: ERC distribution, observed vs simulated ------------
for (k in seq_along(show.py)) {
  o <- ercmat.gm[, ii[k]]; s <- as.numeric(outmat.cs[, jj[k]])
  png(fig.name("ecdf", k), width = 560, height = 500, res = 110)
  par(mar = c(3.6, 3.8, 2.6, 0.8), mgp = c(2.2, 0.7, 0))
  plot(ecdf(o), main = sprintf("py%d %s", show.py[k], names(show.py)[k]),
       xlab = "ERC", ylab = "cumulative prob", col = "steelblue",
       xlim = c(0, max(o, s)), cex.main = 0.98, pch = NA, verticals = TRUE)
  plot(ecdf(s), add = TRUE, col = "tomato", pch = NA, verticals = TRUE)
  abline(v = thr.p80[jj[k]], col = "grey25", lty = 3)
  legend("bottomright", c("observed", "simulated"),
         col = c("steelblue", "tomato"), lty = 1, bty = "n", cex = 0.78)
  dev.off()
}

## ---- numeric companion so the figures can be checked, not just eyeballed
mo.of.doy <- as.integer(format(as.Date(0:364, origin = "2001-01-01"), "%m"))
season.chk <- do.call(rbind, lapply(seq_along(show.py), function(k) {
  o <- ercmat.gm[, ii[k]]; s <- as.numeric(outmat.cs[, jj[k]])
  co <- tapply(o, doy.obs, mean); csm <- tapply(s, doy.sim, mean)
  ro <- rle(o >= thr.p80[jj[k]]); rs <- rle(s >= thr.p80[jj[k]])
  data.frame(region = names(show.py)[k], PYROME = show.py[k],
             P80 = round(thr.p80[jj[k]]),
             #peak_mo_obs = mo.of.doy[which.max(co)],
             #peak_mo_sim = mo.of.doy[which.max(csm)],
             peak_erc_obs = round(max(co), 1), peak_erc_sim = round(max(csm), 1),
             mean_obs = round(mean(o), 1),     mean_sim = round(mean(s), 1),
             days_above_yr_obs = round(sum(o >= thr.p80[jj[k]]) / nyr.obs, 1),
             days_above_yr_sim = round(sum(s >= thr.p80[jj[k]]) / nyr.sm, 1),
             longest_obs = if (any(ro$values)) max(ro$lengths[ro$values]) else 0,
             longest_sim = if (any(rs$values)) max(rs$lengths[rs$values]) else 0,
             stringsAsFactors = FALSE)
}))
write.csv(season.chk, file.path(output_dir, "seasons-check-5pyromes.csv"),
          row.names = FALSE)
cat("\n=== seasonal comparison, 5 representative pyromes ===\n")
print(season.chk, row.names = FALSE)
cat(sprintf("\nwrote %d season figures (overlay/example/ecdf x %d pyromes),\n",
            3L * length(show.py), length(show.py)))
cat("      seasons-check-5pyromes.csv and season-plot-inputs.rds\n")
