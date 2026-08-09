## ============================================================================
## beran_ebvk.R
##
## Shared core functions for the Escobar-Bach & Van Keilegom (2019/2023)
## conditional cure / heavy-tail extrapolation estimator, implemented from
## their eqs. (2.6)-(2.8). Used by every Phase 2 script.
##
## Two versions of the base conditional-KM layer are provided:
##   - make_beran():    single continuous covariate, Epanechnikov kernel
##   - make_beran_md(): multivariate product kernel (Epanechnikov on
##                      continuous covariates, Aitchison-Aitken discount
##                      kernel on categorical/discrete covariates)
##
## On top of either base layer, gamma_hat_fn()/p_hat_fn()/fit_at_x0() apply
## the EBVK tail-index and extrapolated-cure-probability correction.
##
## Note on condSURV::Beran: the write-up (Sec. 5.1) describes the base layer
## as "via condSURV::Beran". condSURV pulls in a heavy, CRAN-only dependency
## chain (np, doParallel, doRNG, foreach) that isn't reachable from this
## sandboxed/offline environment. Beran's (1981) kernel-weighted KM estimator
## is a short, closed-form product-limit formula, so it is reimplemented
## directly here; it is algebraically identical to what condSURV::Beran
## computes for a single continuous covariate with a fixed bandwidth.
## ============================================================================

## ---- Epanechnikov kernel ----
epan <- function(u) ifelse(abs(u) < 1, 0.75 * (1 - u^2), 0)

## ---- Aitchison-Aitken discount kernel for an unordered categorical covariate ----
## weight = 1 for a match, lambda/(m-1)-scaled discount otherwise (m = #categories)
aitchison_aitken <- function(xk, x0k, lambda, m) {
  ifelse(xk == x0k, 1 - lambda, lambda / (m - 1))
}

## ---- single-covariate Beran (1981) conditional KM ----
## T: observed times, delta: event indicator, X: covariate, x0: eval point,
## h: bandwidth. Returns a function Fn(t) = P(T <= t | X = x0).
make_beran <- function(T, delta, X, x0, h) {
  w_raw <- epan((x0 - X) / h)
  if (sum(w_raw) <= 0) stop("bandwidth too small: no observations in window")
  W <- w_raw / sum(w_raw)
  ord <- order(T)
  Ts <- T[ord]; ds <- delta[ord]; Ws <- W[ord]
  cumWprev <- c(0, cumsum(Ws)[-length(Ws)])
  denom <- 1 - cumWprev
  denom[denom <= 1e-12] <- 1e-12
  factor <- 1 - Ws / denom
  factor[factor < 0] <- 0
  logfac <- ifelse(ds == 1, log(factor), 0)
  cumlog <- cumsum(logfac)
  function(t) {
    idx <- which(Ts <= t)
    if (length(idx) == 0) return(0)
    1 - exp(cumlog[max(idx)])
  }
}

## ---- multivariate product-kernel Beran conditional KM ----
## X: data.frame of covariates (numeric columns = continuous -> Epanechnikov;
##    columns flagged in `categorical` = discrete -> Aitchison-Aitken).
## x0: named list/row with the same column names, h: named vector of
## per-column bandwidths (continuous) / lambda values (categorical).
make_beran_md <- function(T, delta, X, x0, h, categorical = character(0)) {
  stopifnot(all(names(x0) %in% names(X)), all(names(h) %in% names(X)))
  vars <- names(x0)
  w_raw <- rep(1, nrow(X))
  for (v in vars) {
    if (v %in% categorical) {
      m <- length(unique(X[[v]]))
      w_raw <- w_raw * aitchison_aitken(X[[v]], x0[[v]], h[[v]], m)
    } else {
      w_raw <- w_raw * epan((x0[[v]] - X[[v]]) / h[[v]])
    }
  }
  if (sum(w_raw) <= 0) stop("bandwidth too small: no observations in window")
  W <- w_raw / sum(w_raw)
  ess <- 1 / sum(W^2)   # Kish's effective sample size
  ord <- order(T)
  Ts <- T[ord]; ds <- delta[ord]; Ws <- W[ord]
  cumWprev <- c(0, cumsum(Ws)[-length(Ws)])
  denom <- 1 - cumWprev
  denom[denom <= 1e-12] <- 1e-12
  factor <- 1 - Ws / denom
  factor[factor < 0] <- 0
  logfac <- ifelse(ds == 1, log(factor), 0)
  cumlog <- cumsum(logfac)
  Fn <- function(t) {
    idx <- which(Ts <= t)
    if (length(idx) == 0) return(0)
    1 - exp(cumlog[max(idx)])
  }
  attr(Fn, "ess") <- ess
  Fn
}

## ---- EBVK eq (2.6): tail-index (Hill-type) ratio estimator ----
## Truncated at 0.1 from below (validated fix, see Phase 2 synthetic check):
## the paper's y1^(-1/gamma) term blows up / goes out of [0,1] range without it.
GAMMA_FLOOR <- 0.1

gamma_hat_fn <- function(Fn, tau_n, y2) {
  num <- Fn(y2^2 * tau_n) - Fn(y2 * tau_n)
  den <- Fn(y2 * tau_n) - Fn(tau_n)
  ratio <- num / den
  if (!is.finite(ratio) || ratio <= 0) return(NA_real_)
  -1 / (log(ratio) / log(y2))
}

## ---- EBVK eq (2.7): extrapolated cure-probability estimator ----
p_hat_fn <- function(Fn, tau_n, y1, gam) {
  if (is.na(gam)) return(NA_real_)
  Fn_tau <- Fn(tau_n)
  Fn_tau + (Fn_tau - Fn(y1 * tau_n)) / (y1^(-1 / gam) - 1)
}

## ---- tuning-constant grid, as in the paper's search over (y1, y2) ----
G_GRID <- seq(0.25, 0.89, by = 0.02)

## ---- eq (2.8) selection rule, algebraically simplified: ----
## argmin_k sum_j (Pk - Pj)^2  <=>  argmin_k |Pk - mean(P)|
## (the "quadratic errors criterion"), restricted to grid points whose
## implied p_hat falls in the theoretically valid range [Fn(tau_n), 1]
## (validated fix, see Phase 2 real-data application).
fit_at_x0 <- function(Fn, tau_n, grid = G_GRID) {
  beran_p <- Fn(tau_n)
  gam_grid <- pmax(sapply(grid, function(y2) gamma_hat_fn(Fn, tau_n, y2)), GAMMA_FLOOR)
  pmat <- matrix(NA_real_, length(grid), length(grid))
  for (j in seq_along(grid)) for (i in seq_along(grid))
    pmat[i, j] <- p_hat_fn(Fn, tau_n, grid[i], gam_grid[j])
  valid <- is.finite(pmat) & pmat >= beran_p - 1e-8 & pmat <= 1
  if (!any(valid)) return(list(beran = beran_p, corrected = NA_real_, gamma = NA_real_))
  mean_p <- mean(pmat[valid])
  d <- abs(pmat - mean_p); d[!valid] <- Inf
  best <- which(d == min(d), arr.ind = TRUE)[1, ]
  list(beran = beran_p, corrected = pmat[best[1], best[2]], gamma = gam_grid[best[2]])
}

## ---- convenience wrapper: single covariate, one evaluation point ----
fit_one <- function(Tobs, delta, X, x0, h) {
  Fn <- make_beran(Tobs, delta, X, x0, h)
  tau_n <- max(Tobs)
  fit_at_x0(Fn, tau_n)
}

## ---- convenience wrapper: multivariate, one evaluation point ----
fit_one_md <- function(Tobs, delta, Xdf, x0, h, categorical = character(0)) {
  Fn <- make_beran_md(Tobs, delta, Xdf, x0, h, categorical)
  tau_n <- max(Tobs)
  res <- fit_at_x0(Fn, tau_n)
  res$ess <- attr(Fn, "ess")
  res
}
