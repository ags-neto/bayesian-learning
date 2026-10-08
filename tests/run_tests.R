#!/usr/bin/env Rscript
# run_tests.R -- test suite for the bayesian-learning repository.
#
# What it does: loads `forestfire.txt` exactly as the report does, runs the
# three JAGS models of "Bayesian Learning and Monte Carlo Simulation - Final
# Project.Rmd" (in the same order and with the same priors), and checks a set of
# facts that are known independently of the run (data shape, summary statistics
# of the response, posterior location/scale, convergence, DIC ordering).
#
# It performs one negative check: with a corrupted copy of the data the
# analysis must abort with an error instead of producing numbers.
#
# Exit code: 0 when every check passes, 1 otherwise.
#
# The model definitions are copied verbatim from the .Rmd, except that the
# categorical variable `month` is coded as integer 1..12 because JAGS cannot
# take a character vector. That coding does not change the fitted values: the
# `month` regression column is a single numeric vector either way, and only its
# linear coefficient (itself uninterpreted in the report) depends on the coding.
#
# Environment variables (defaults are the values used in the report):
#   BL_BURNIN (2000)  BL_SAMPLE (5000)  BL_LOG (tests/run_tests.log)

suppressWarnings(suppressMessages({
  ok_pkgs <- sapply(c("runjags", "dplyr"), requireNamespace, quietly = TRUE)
}))
if (!all(ok_pkgs)) {
  cat("MISSING R PACKAGES:", paste(names(ok_pkgs)[!ok_pkgs], collapse = ", "), "\n")
  quit(status = 1)
}
suppressWarnings(suppressMessages({
  library(runjags)
  library(dplyr)
}))

# ---------------------------------------------------------------- assertions
checks <- new.env(parent = emptyenv())
checks$rows <- list()
ok_pkgs <- NULL

record <- function(what, expected, obtained, pass) {
  checks$rows[[length(checks$rows) + 1L]] <- data.frame(
    check = what, expected = expected, obtained = obtained,
    pass = if (isTRUE(pass)) "yes" else "NO", stringsAsFactors = FALSE)
}

expect_true <- function(what, condition, expected, obtained = if (isTRUE(condition)) "TRUE" else "FALSE") {
  record(what, expected, obtained, isTRUE(condition))
}

expect_at_least <- function(what, got, floor_value) {
  record(what, paste0("> ", floor_value), as.character(got), isTRUE(got > floor_value))
}

fnum <- function(x, d = 3) formatC(as.numeric(x), format = "f", digits = d)

# ------------------------------------------------------------- repository map
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1L) args[[1]] else ".", mustWork = TRUE)
data_file <- file.path(root, "forestfire.txt")
rmd_file <- file.path(root, "Bayesian Learning and Monte Carlo Simulation - Final Project.Rmd")
pdf_file <- file.path(root, "Bayesian Learning and Monte Carlo Simulation - Final Project.pdf")

log_file <- Sys.getenv("BL_LOG", file.path(root, "tests", "run_tests.log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
log_con <- file(log_file, open = "wt")

emit <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_con)
}

burnin <- as.integer(Sys.getenv("BL_BURNIN", "5000"))
sample_n <- as.integer(Sys.getenv("BL_SAMPLE", "5000"))
set.seed(as.integer(Sys.getenv("BL_SEED", "20260926")))

emit("== bayesian-learning | run_tests.R ==")
emit("R ", R.version.string, " | ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | root=", root)

# 1. static facts about the checked-in artefacts -----------------------------
emit("\n-- artefacts --")
for (f in c(data_file, rmd_file, pdf_file)) {
  expect_true(paste0("exists: ", basename(f)), file.exists(f), "file",
              if (file.exists(f)) "present" else "missing")
}
rmd_lines <- length(readLines(rmd_file))
pdf_size <- file.info(pdf_file)$size
emit("Rmd lines: ", rmd_lines, " | PDF bytes: ", pdf_size)
expect_at_least("Rmd line count", rmd_lines, 300)

# 2. the data ----------------------------------------------------------------
emit("\n-- data (forestfire.txt) --")
run_analysis <- function(dat) {
  forest <- dat %>% group_by(month)
  y <- forest$area
  N <- length(y)
  s0 <- 0.01
  s1 <- 1 / var(forest$FFMC); s2 <- 1 / var(forest$DMC); s3 <- 1 / var(forest$DC)
  s4 <- 1 / var(forest$ISI); s5 <- 1 / var(forest$temp); s6 <- 1 / var(forest$RH)
  s7 <- 1 / var(forest$wind); s8 <- 1 / var(forest$rain)
  # JAGS cannot take a character vector; code the categorical variable as 1..12
  month4jags <- as.integer(factor(forest$month, levels = sort(unique(forest$month))))

  model_string <- "model {
for (i in 1:N){
   y[i] ~ dnorm(beta0 + beta1*FFMC[i] +
              beta2*DMC[i]+beta3*DC[i]+beta4*ISI[i]+beta5*temp[i]+beta6*RH[i]+beta7*wind[i]+beta8*rain[i]+beta9*month[i], invsigma2)
}
beta0 ~ dnorm(mu0, sigma_0)
beta1 ~ dnorm(mu1, sigma_1)
beta2 ~ dnorm(mu2, sigma_2)
beta3 ~ dnorm(mu3, sigma_3)
beta4 ~ dnorm(mu4, sigma_4)
beta5 ~ dnorm(mu5, sigma_5)
beta6 ~ dnorm(mu6, sigma_6)
beta7 ~ dnorm(mu7, sigma_7)
beta8 ~ dnorm(mu8, sigma_8)
beta9 ~ dnorm(mu9, sigma_9)
invsigma2 ~ dgamma(a, b)
sigma <- sqrt(pow(invsigma2, -1))
}"
  the_data <- list("y" = y, "FFMC" = forest$FFMC, "DMC" = forest$DMC, "DC" = forest$DC,
                   "ISI" = forest$ISI, "temp" = forest$temp, "RH" = forest$RH,
                   "wind" = forest$wind, "rain" = forest$rain, "month" = month4jags, "N" = N,
                   "mu0" = 0, "sigma_0" = s0,
                   "mu1" = 0, "sigma_1" = s1, "mu2" = 0, "sigma_2" = s2,
                   "mu3" = 0, "sigma_3" = s3, "mu4" = 0, "sigma_4" = s4,
                   "mu5" = 0, "sigma_5" = s5, "mu6" = 0, "sigma_6" = s6,
                   "mu7" = 0, "sigma_7" = s7, "mu8" = 0, "sigma_8" = s8,
                   "mu9" = 0, "sigma_9" = 0.001, "a" = 0.001, "b" = 0.001)
  # least-squares fit, used only to spread the initial values of the chains
  x1 <- data.frame(y = y, FFMC = forest$FFMC, DMC = forest$DMC, DC = forest$DC,
                   ISI = forest$ISI, temp = forest$temp, RH = forest$RH,
                   wind = forest$wind, rain = forest$rain,
                   month = as.factor(month4jags))
  coef1 <- stats::coef(stats::lm(y ~ ., data = x1))
  init1 <- function(chain) {
    set.seed(2000L + chain)
    list(beta0 = coef1[[1]] * (1 + stats::rnorm(1, 0, 0.1)),
         beta1 = coef1[[2]], beta2 = coef1[[3]], beta3 = coef1[[4]],
         beta4 = coef1[[5]], beta5 = coef1[[6]], beta6 = coef1[[7]],
         beta7 = coef1[[8]], beta8 = coef1[[9]], beta9 = coef1[[10]],
         invsigma2 = 1 / stats::var(y))
  }
  posterior <- run.jags(model_string, n.chains = 2, data = the_data,
                        monitor = c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5",
                                    "beta6", "beta7", "beta8", "beta9", "sigma"),
                        inits = init1,
                        burnin = burnin, sample = sample_n)

  model_string1 <- "model {
for (i in 1:N){
   y[i] ~ dnorm(beta0 +beta5*temp[i]+beta6*RH[i]+beta7*wind[i]+beta8*rain[i], invsigma2)
}
beta0 ~ dnorm(mu0, sigma_0)
beta5 ~ dnorm(mu5, sigma_5)
beta6 ~ dnorm(mu6, sigma_6)
beta7 ~ dnorm(mu7, sigma_7)
beta8 ~ dnorm(mu8, sigma_8)
invsigma2 ~ dgamma(a, b)
sigma <- sqrt(pow(invsigma2, -1))
}"
  the_data1 <- list("y" = y, "temp" = forest$temp, "RH" = forest$RH,
                    "wind" = forest$wind, "rain" = forest$rain, "N" = N,
                    "mu0" = 0, "sigma_0" = s0,
                    "mu5" = 0, "sigma_5" = s5, "mu6" = 0, "sigma_6" = s6,
                    "mu7" = 0, "sigma_7" = s7, "mu8" = 0, "sigma_8" = s8,
                    "a" = 0.001, "b" = 0.001)
  x2 <- data.frame(y = y, temp = forest$temp, RH = forest$RH,
                   wind = forest$wind, rain = forest$rain)
  coef2 <- stats::coef(stats::lm(y ~ ., data = x2))
  init2 <- function(chain) {
    set.seed(3000L + chain)
    list(beta0 = coef2[[1]] * (1 + stats::rnorm(1, 0, 0.1)),
         beta5 = coef2[[2]], beta6 = coef2[[3]], beta7 = coef2[[4]],
         beta8 = coef2[[5]], invsigma2 = 1 / stats::var(y))
  }
  posterior1 <- run.jags(model_string1, n.chains = 2, data = the_data1,
                         monitor = c("beta0", "beta5", "beta6", "beta7", "beta8", "sigma"),
                         inits = init2,
                         burnin = burnin, sample = sample_n)

  model_string2 <- "model {
for (i in 1:N){
   y[i] ~ dnorm(beta0 +beta5, invsigma2)
}
beta0 ~ dnorm(mu0, sigma_0)
beta5 ~ dnorm(mu5, sigma_5)
invsigma2 ~ dgamma(a, b)
sigma <- sqrt(pow(invsigma2, -1))
}"
  the_data2 <- list("y" = y, "temp" = forest$temp, "N" = N,
                    "mu0" = 0, "sigma_0" = s0,
                    "mu5" = 0, "sigma_5" = s5,
                    "a" = 0.001, "b" = 0.001)
  coef3 <- stats::coef(stats::lm(y ~ 1, data = data.frame(y = y)))
  init3 <- function(chain) {
    set.seed(4000L + chain)
    list(beta0 = coef3[[1]] * (1 + stats::rnorm(1, 0, 0.1)),
         beta5 = coef3[[1]] / mean(forest$temp),
         invsigma2 = 1 / stats::var(y))
  }
  posterior2 <- run.jags(model_string2, n.chains = 2, data = the_data2,
                         monitor = c("beta0", "beta5", "sigma"),
                         inits = init3,
                         burnin = burnin, sample = sample_n)

  # DIC = Dbar + pD. runjags returns, per observation, the MCMC samples of the
  # deviance and the penalty (pD/n), so Dbar is the sum over observations of the
  # posterior mean of the per-observation deviance. It also prints the rounded
  # "Penalized deviance" line, which the run was checked against.
  dic <- function(p) {
    invisible(capture.output(d <- extract.runjags(p, "dic")))
    dev <- unlist(d$deviance)
    Dbar <- sum(colMeans(matrix(dev, ncol = length(dev))))
    as.numeric(Dbar + d$penalty[[1]])
  }
  summarise <- function(p, vars) {
    s <- summary(p)
    if (is.list(s)) s <- s$statistics
    out <- list(mean = s[vars, "Mean"], lower = s[vars, "Lower95"], upper = s[vars, "Upper95"])
    names(out$mean) <- vars; names(out$lower) <- vars; names(out$upper) <- vars
    out
  }
  psrf_of <- function(p, vars) {
    s <- summary(p)
    if (is.list(s)) {
      r <- s$psrf[, "Point est."]
    } else {
      r <- s[vars, "psrf"]
    }
    r
  }
  list(
    n_obs = N,
    dim = dim(dat),
    area_zero = sum(dat$area == 0),
    area_mean = mean(dat$area),
    area_sd = sd(dat$area),
    area_max = max(dat$area),
    m1 = summarise(posterior, c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5",
                                "beta6", "beta7", "beta8", "beta9", "sigma")),
    m2 = summarise(posterior1, c("beta0", "beta5", "beta6", "beta7", "beta8", "sigma")),
    m3 = summarise(posterior2, c("beta0", "beta5", "sigma")),
    psrf = list(m1 = psrf_of(posterior, c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5",
                                          "beta6", "beta7", "beta8", "beta9", "sigma")),
                m2 = psrf_of(posterior1, c("beta0", "beta5", "beta6", "beta7", "beta8", "sigma")),
                m3 = psrf_of(posterior2, c("beta0", "beta5", "sigma"))),
    dic = c(model1_all = dic(posterior), model2_weather = dic(posterior1),
            model3_temp = dic(posterior2))
  )
}

forestfire <- read.table(data_file)
emit("dim: ", nrow(forestfire), " x ", ncol(forestfire), " | columns: ",
     paste(names(forestfire), collapse = ","))
expect_true("data rows", nrow(forestfire) == 517, "517", as.character(nrow(forestfire)))
expect_true("data columns", ncol(forestfire) == 13, "13", as.character(ncol(forestfire)))
expect_true("required columns present",
            all(c("X", "Y", "month", "day", "FFMC", "DMC", "DC", "ISI", "temp", "RH",
                  "wind", "rain", "area") %in% names(forestfire)),
            "13 named columns", paste(names(forestfire), collapse = ","))
expect_true("no missing values", !any(is.na(forestfire)), "0 NA", as.character(sum(is.na(forestfire))))
expect_true("area is right-skewed and zero-inflated",
            max(forestfire$area) > 1000 && sum(forestfire$area == 0) > 200,
            "max>1000 & >200 zeros",
            paste0("max=", max(forestfire$area), " zeros=", sum(forestfire$area == 0)))

# 3. the analysis -------------------------------------------------------------
emit("\n-- models (JAGS) -- burnin=", burnin, " sample=", sample_n, " chains=2")
t0 <- Sys.time()
res <- run_analysis(forestfire)
emit("elapsed: ", fnum(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
emit("posterior means  model1(all):  ", paste(names(res$m1$mean), fnum(res$m1$mean, 3), sep = "=", collapse = " "))
emit("posterior means  model2(wx):   ", paste(names(res$m2$mean), fnum(res$m2$mean, 3), sep = "=", collapse = " "))
emit("posterior means  model3(temp): ", paste(names(res$m3$mean), fnum(res$m3$mean, 3), sep = "=", collapse = " "))

expect_true("posterior sigma (model 3) within 5 of the raw sd of area (63.66)",
            abs(res$m3$mean[["sigma"]] - 63.66) <= 5, "|mean-63.66|<=5",
            fnum(res$m3$mean[["sigma"]], 2))
expect_true("posterior interval for beta5 (model 2) covers the report's mean (0.918)",
            res$m2$lower[["beta5"]] < 0.918 && res$m2$upper[["beta5"]] > 0.918,
            "0.918 inside interval",
            paste0(fnum(res$m2$lower[["beta5"]]), " .. ", fnum(res$m2$upper[["beta5"]])))
expect_true("posterior mean of beta5 (model 3) is near the report's 3.2",
            abs(res$m3$mean[["beta5"]] - 3.2) < 1.0, "3.2 +/- 1.0",
            fnum(res$m3$mean[["beta5"]], 3))
# Recorded, not corrected: the report text calls temperature "the only significant
# predictor" and cites the interval (0.199, 0.916) for it, but the table it prints for
# that model (beta5 lower95 -19.4, upper95 16.2) has an interval that contains zero, and
# the paragraph below it says so in the same terms. Both statements are in the .Rmd and
# both are left as they are. The same report concludes from a DIC spread of five units
# (5768 / 5763 / 5765) that model 2 is the most plausible; the run below reproduces that
# ordering but not that level of confidence.
expect_true("model 1 interval for beta1 (FFMC) covers zero, as the report states",
            res$m1$lower[["beta1"]] < 0 && res$m1$upper[["beta1"]] > 0,
            "lower<0<upper",
            paste0(fnum(res$m1$lower[["beta1"]]), " .. ", fnum(res$m1$upper[["beta1"]])))
max_psrf <- max(unlist(list(res$psrf$m1, res$psrf$m2, res$psrf$m3)))
expect_true("all three models converge (psrf <= 1.05)", max_psrf <= 1.05, "<=1.05",
            fnum(max_psrf, 4))
expect_true("model 1 exposes all 11 monitored parameters",
            length(res$m1$mean) == 11, "11", as.character(length(res$m1$mean)))
emit("DIC: ", paste(names(res$dic), fnum(res$dic, 1), sep = "=", collapse = " "))
expect_true("DIC of the three models is finite and near the report's 5763-5768",
            all(is.finite(res$dic)) && all(res$dic > 5700 & res$dic < 5800),
            "within 5700-5800", paste(fnum(res$dic, 1), collapse = "/"))
dic_ok <- all(is.finite(res$dic)) && all(res$dic > 5000) && which.min(res$dic) == 2
expect_true("DIC finite and model 2 (weather) has the lowest DIC, as the report argues",
            dic_ok, "min at model2_weather",
            paste(names(res$dic)[which.min(res$dic)], fnum(min(res$dic, na.rm = TRUE), 1)))
emit("DIC ordering obtained: ", paste(names(res$dic)[order(res$dic)], collapse = " < "))

# 4. negative check: corrupted data must fail, not produce numbers -----------
emit("\n-- negative check: corrupted data --")
corrupt <- forestfire[, setdiff(names(forestfire), "temp")]
neg <- tryCatch({
  run_analysis(corrupt)
  "NO ERROR - the analysis returned numbers from corrupted data"
}, error = function(e) paste("aborted as expected:", conditionMessage(e)))
neg_pass <- !grepl("NO ERROR", neg)
emit(neg)
expect_true("corrupted data aborts the analysis", neg_pass,
            "error raised", if (neg_pass) "error raised" else neg)

missing <- tryCatch({
  read.table(file.path(tempdir(), "does-not-exist.txt"))
  "NO ERROR - missing file was accepted"
}, error = function(e) paste("aborted as expected:", conditionMessage(e)))
expect_true("missing data file aborts the load", grepl("aborted", missing),
            "error raised", if (grepl("aborted", missing)) "error raised" else missing)

# 5. the same thing end to end: the suite run against a corrupted repository must
#    itself exit non-zero, otherwise a broken checkout would still report green
suite_failed <- tryCatch({
  bad_root <- file.path(tempdir(), "corrupted-repo")
  dir.create(bad_root, showWarnings = FALSE)
  file.copy(data_file, bad_root, overwrite = TRUE)
  bad <- read.table(data_file)
  bad$temp <- NULL
  utils::write.table(bad, file.path(bad_root, "forestfire.txt"), row.names = FALSE)
  for (nm in c("Bayesian Learning and Monte Carlo Simulation - Final Project.Rmd",
               "Bayesian Learning and Monte Carlo Simulation - Final Project.pdf")) {
    file.copy(file.path(root, nm), bad_root, overwrite = TRUE)
  }
  dir.create(file.path(bad_root, "tests"), showWarnings = FALSE)
  self <- file.path(root, "tests", "run_tests.R")
  child_log <- file.path(bad_root, "tests", "run_tests.log")
  child_stdout <- file.path(bad_root, "child-output.txt")
  code <- suppressWarnings(system2("Rscript", c(shQuote(self), shQuote(bad_root)),
                                   env = c("LANG=C", "LC_ALL=C",
                                           paste0("BL_LOG=", shQuote(child_log)),
                                           "BL_BURNIN=200", "BL_SAMPLE=500"),
                                   stdout = child_stdout, stderr = child_stdout))
  child_text <- paste(readLines(child_stdout), collapse = " ")
  if (code == 0) "NO ERROR - the suite passed on corrupted data"
  else if (!grepl("Column `month` is not found|'x' is NULL|Error", child_text))
    paste("exit", code, "without a visible error")
  else paste0("suite exited ", code, " and aborted on the corrupted column, as expected")
}, error = function(e) paste("aborted as expected:", conditionMessage(e)))
expect_true("the test suite fails on a corrupted repository (exit != 0)",
            grepl("^suite exited ", suite_failed), "non-zero exit", suite_failed)

# --------------------------------------------------------------------- report
emit("\n== verification table ==")
tbl <- do.call(rbind, checks$rows)
widths <- sapply(tbl, function(x) max(nchar(as.character(x)), nchar("obtained")))
pad <- function(x, w) formatC(as.character(x), width = w, flag = "-")
for (i in seq_len(nrow(tbl))) {
  emit(paste0("| ", pad(tbl$check[i], widths[["check"]]), " | ",
              pad(tbl$expected[i], widths[["expected"]]), " | ",
              pad(tbl$obtained[i], widths[["obtained"]]), " | ",
              pad(tbl$pass[i], widths[["pass"]]), " |"))
}
failed <- sum(tbl$pass != "yes")
emit("checks: ", nrow(tbl), " | failed: ", failed)
emit(if (failed == 0) "RESULT: PASS" else "RESULT: FAIL")
close(log_con)
cat("log written to ", log_file, "\n", sep = "")
quit(save = "no", status = if (failed == 0) 0 else 1)
