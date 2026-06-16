############################################################################################
# Multivariate rescaling via Gaussian optimal transport of dependence structure (linear
# Monge–Kantorovich transform) with different options for handling non-stationarity
# of correlations.

MRSmk <- function(o.c, m.c, m.p, mhat.c, mhat.p,
  cor_change = c("logdelta", "fisherz"), ties.method = "random") {

  cor_change <- match.arg(cor_change)

  stopifnot(is.matrix(o.c), is.matrix(m.c), is.matrix(m.p), is.matrix(mhat.c),
    is.matrix(mhat.p))
  stopifnot(ncol(o.c) == ncol(m.c), ncol(m.c) == ncol(m.p),
    ncol(mhat.c) == ncol(mhat.p))
  stopifnot(ncol(mhat.c) == ncol(o.c), ncol(mhat.p) == ncol(o.c))

  if (anyNA(o.c) || anyNA(m.c) || anyNA(m.p) || anyNA(mhat.c) || 
    anyNA(mhat.p)) {
    stop("Inputs contain NA; handle missingness before calling.")
  }
  if (ncol(o.c) <= 1 || ncol(m.c) <= 1 || ncol(mhat.c) <= 1 ||
    ncol(mhat.p) <= 1) {
    stop("Inputs must have more than one column.")
  }

  symm <- function(A) (A + t(A)) / 2

  # Symmetric positive definite (SPD) helper
  ensure_spd <- function(S, tol = 1e-12) {
    S <- (S + t(S)) / 2
    S <- as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
    S <- (S + t(S)) / 2
    eg <- eigen(S, symmetric = TRUE)
    lam <- pmax(eg$values, tol)
    S <- eg$vectors %*% diag(lam, nrow(S)) %*% t(eg$vectors)
    (S + t(S)) / 2
  }

  # Projects to nearest PD correlation matrix (Matrix::nearPD, corr=TRUE) and
  # then renormalizes to enforce unit diagonal via SPD-preserving scaling
  ensure_cor <- function(R, tol = 1e-12) {
    R <- (R + t(R)) / 2
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    R <- (R + t(R)) / 2
    d <- sqrt(pmax(diag(R), tol))
    Dinv <- diag(1 / d, nrow(R))
    R <- (Dinv %*% R %*% Dinv + t(Dinv %*% R %*% Dinv)) / 2
    diag(R) <- 1
    R
  }

  # Clip correlations away from +/-1 for numerical stability.
  clip_cor <- function(R, eps = 1e-7) {
    R <- pmin(pmax(R, -1 + eps), 1 - eps)
    symm(R)
  }

  eig_clip <- function(S, tol = 1e-12) {
    S <- ensure_spd(S, tol = tol)
    eg <- eigen(S, symmetric = TRUE)
    list(U = eg$vectors, lam = pmax(eg$values, tol))
  }

  sqrtm_spd <- function(S, tol = 1e-12) {
    ec <- eig_clip(S, tol = tol)
    (ec$U %*% diag(sqrt(ec$lam), nrow(S)) %*% t(ec$U) + 
      t(ec$U %*% diag(sqrt(ec$lam), nrow(S)) %*% t(ec$U))) / 2
  }

  invsqrtm_spd <- function(S, tol = 1e-12) {
    ec <- eig_clip(S, tol = tol)
    (ec$U %*% diag(1 / sqrt(ec$lam), nrow(S)) %*% t(ec$U) + 
      t(ec$U %*% diag(1 / sqrt(ec$lam), nrow(S)) %*% t(ec$U))) / 2
  }

  logm_spd <- function(S, tol = 1e-12) {
    ec <- eig_clip(S, tol = tol)
    (ec$U %*% diag(log(ec$lam), nrow(S)) %*% t(ec$U) + 
      t(ec$U %*% diag(log(ec$lam), nrow(S)) %*% t(ec$U))) / 2
  }

  expm_sym <- function(A) {
    A <- symm(A)
    eg <- eigen(A, symmetric = TRUE)
    symm(eg$vectors %*% diag(exp(eg$values), nrow(A)) %*% t(eg$vectors))
  }

  # Linear Monge–Kantorovich (Gaussian OT) factor on correlation matrices
  mk_factor <- function(Rm, Ro) {
    Rm <- ensure_cor(Rm)
    Ro <- ensure_cor(Ro)
    Rm_s  <- sqrtm_spd(Rm)
    Rm_is <- invsqrtm_spd(Rm)
    A <- Rm_s %*% Ro %*% Rm_s
    T <- Rm_is %*% sqrtm_spd(A) %*% Rm_is
    symm(T)
  }

  # Fisher-Z delta:
  # z = atanh(r), add delta in z-space, invert with tanh.
  # Avoid atanh(1) on diagonal by forcing diag to (1-epsZ) before atanh.
  fisherz_delta <- function(Roc, Rmc, Rmp, epsZ = 1e-7) {

    Roc <- ensure_cor(clip_cor(Roc, epsZ))
    Rmc <- ensure_cor(clip_cor(Rmc, epsZ))
    Rmp <- ensure_cor(clip_cor(Rmp, epsZ))

    # Make diagonals strictly inside (-1,1) before atanh
    diag(Roc) <- 1 - epsZ
    diag(Rmc) <- 1 - epsZ
    diag(Rmp) <- 1 - epsZ

    Zoc <- atanh(Roc)
    Zmc <- atanh(Rmc)
    Zmp <- atanh(Rmp)

    Rt  <- tanh(Zoc + (Zmp - Zmc))
    Rt  <- ensure_cor(clip_cor(Rt, epsZ))
    diag(Rt) <- 1
    Rt
  }

  # Van der Waerden normal score transform
  vdw_normal_scores <- function(X, ties.method = "random") {
    if (!is.matrix(X)) X <- as.matrix(X)
    n <- nrow(X)
    if (n < 2) stop("Need at least 2 rows to compute ranks/normal scores.")

    Z <- apply(X, 2, function(col) {
      r <- rank(col, ties.method = ties.method, na.last = "keep")
      u <- r / (n + 1)
      qnorm(u)
    })

    dimnames(Z) <- dimnames(X)
    Z
  }

  # Rank reordering
  rank_reorder <- function(source, target) {
    if (!is.matrix(source) || !is.numeric(source)) {
      stop("`source` must be a numeric matrix.")
    }
    if (!is.matrix(target) || !is.numeric(target)) {
      stop("`target` must be a numeric matrix.")
    }
    if (nrow(source) != nrow(target)) {
      stop("`source` and `target` must have the same number of rows (cases).")
    }

    out <- source
    p <- ncol(source)

    for (c in seq_len(p)) {
      x <- source[, c]
      y <- target[, c]

      idx <- which(!is.na(x) & !is.na(y))
      if (length(idx) <= 1L) next

      # Order rows by target ranks, then assign sorted source values to those rows.
      ord_rows <- idx[order(y[idx], method = "radix")]
      x_sorted <- sort(x[idx], na.last = TRUE)

      out[ord_rows, c] <- x_sorted
    }

    out
  }

  # Normalize datasets
  oZ <- vdw_normal_scores(o.c, ties.method = ties.method)
  cZ <- vdw_normal_scores(m.c, ties.method = ties.method)
  pZ <- vdw_normal_scores(m.p, ties.method = ties.method)
  mhatcZ <- vdw_normal_scores(mhat.c, ties.method = ties.method)
  mhatpZ <- vdw_normal_scores(mhat.p, ties.method = ties.method)

  # Correlation matrices
  Ro.c <- cor(oZ)
  Rm.c <- cor(cZ)
  Rm.p <- cor(pZ)
  Rmhat.c <- cor(mhatcZ)
  Rmhat.p <- cor(mhatpZ)

  # Target projected correlation Ro.p*
  Ro.p <- switch(
    cor_change,

    logdelta = {
      # Ro.p = exp( log(Ro.c) + (log(Rm.p) - log(Rm.c)) )
      # (log/exp are applied via eigenvalues with clipping to avoid log(0))
      Lo.c <- logm_spd(ensure_cor(Ro.c))
      Lm.c <- logm_spd(ensure_cor(Rm.c))
      Lm.p <- logm_spd(ensure_cor(Rm.p))
      ensure_cor(expm_sym(Lo.c + (Lm.p - Lm.c)))
    },

    fisherz = {
      fisherz_delta(Ro.c, Rm.c, Rm.p)
    }

  )

  # MK maps in correlation space
  T.c <- mk_factor(Rmhat.c, Ro.c)
  T.p <- mk_factor(Rmhat.p, Ro.p)

  # Apply to normalized model anomalies
  cZ_hat <- mhatcZ %*% T.c
  pZ_hat <- mhatpZ %*% T.p

  # Reorder
  mhat.c <- rank_reorder(mhat.c, cZ_hat)
  mhat.p <- rank_reorder(mhat.p, pZ_hat)

  out <- list(mhat.c = mhat.c, mhat.p = mhat.p)
  out
}


############################################################################################
