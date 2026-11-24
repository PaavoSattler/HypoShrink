
#' Check a Companion Representation
#'
#' Checks whether two matrix-vector pairs preserve the same quadratic
#' function, up to a specified numerical tolerance.
#'
#' @param H Original numeric hypothesis matrix.
#' @param L Numeric matrix of the proposed companion representation.
#' @param y Optional right-hand side for H. Defaults to zero.
#' @param ytilde Optional right-hand side for L. Defaults to zero.
#' @param tol Nonnegative relative tolerance.
#'
#' @return A list containing:
#' \describe{
#'   \item{valid}{Whether all three companion conditions are satisfied.}
#'   \item{gram_error}{Relative Gram-matrix error.}
#'   \item{cross_error}{Relative cross-product error.}
#'   \item{norm_error}{Relative squared-norm error.}
#'   \item{rank}{Numerical rank of H, determined by QR.}
#'   \item{minimal}{Whether L has exactly rank(H) rows.}
#' }
#'
#' @details
#' A valid companion representation satisfies, within numerical tolerance,
#' \deqn{L^T L = H^T H,}
#' \deqn{L^T \tilde{y} = H^T y,}
#' and
#' \deqn{\|\tilde{y}\|^2 = \|y\|^2.}
#'
#' These conditions preserve the quadratic function for every parameter
#' vector. The check is restricted to equality without a scale factor.
#' Minimality is assessed separately using the numerical QR rank of H.
#'
#' @examples
#' H <- rbind(c(1, -1), c(2, -2))
#' result <- CompanionHypothesis(H)
#' checkCompanion(H, result$L)
#'
#' @export
checkCompanion <- function(
    H, L, y = NULL, ytilde = NULL,
    tol = sqrt(.Machine$double.eps)) {

  if (!is.matrix(H) || !is.numeric(H) ||
      nrow(H) == 0L || ncol(H) == 0L ||
      any(!is.finite(H))) {
    stop("H must be a non-empty finite numeric matrix.")
  }

  if (!is.matrix(L) || !is.numeric(L) ||
      ncol(L) != ncol(H) || any(!is.finite(L))) {
    stop("L must be a finite numeric matrix with the same number of columns as H.")
  }

  if (is.null(y)) {
    y <- rep(0, nrow(H))
  }
  if (is.null(ytilde)) {
    ytilde <- rep(0, nrow(L))
  }

  if (!is.numeric(y) || !is.null(dim(y)) ||
      length(y) != nrow(H) || any(!is.finite(y))) {
    stop("y must be a finite numeric vector matching the rows of H.")
  }

  if (!is.numeric(ytilde) || !is.null(dim(ytilde)) ||
      length(ytilde) != nrow(L) || any(!is.finite(ytilde))) {
    stop("ytilde must be a finite numeric vector matching the rows of L.")
  }

  if (!is.numeric(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol < 0) {
    stop("tol must be a single finite nonnegative number.")
  }

  # Scale inputs to avoid overflow and underflow in quadratic products.
  scale_H <- max(abs(H), abs(L))
  scale_y <- max(abs(y), abs(ytilde))

  Hs <- if (scale_H > 0) H / scale_H else H
  Ls <- if (scale_H > 0) L / scale_H else L
  ys <- if (scale_y > 0) y / scale_y else y
  yts <- if (scale_y > 0) ytilde / scale_y else ytilde

  vecnorm <- function(x) {
    sqrt(sum(x^2))
  }

  relative_error <- function(x1, x2, scale) {
    if (scale == 0) {
      return(0)
    }
    vecnorm(x1 - x2) / scale
  }

  G1 <- crossprod(Hs)
  G2 <- crossprod(Ls)

  gram_error <- relative_error(
    G1, G2,
    max(vecnorm(G1), vecnorm(G2))
  )

  C1 <- drop(crossprod(Hs, ys))
  C2 <- drop(crossprod(Ls, yts))

  cross_error <- relative_error(
    C1, C2,
    max(
      vecnorm(Hs) * vecnorm(ys),
      vecnorm(Ls) * vecnorm(yts)
    )
  )

  N1 <- sum(ys^2)
  N2 <- sum(yts^2)

  norm_error <- if (max(N1, N2) == 0) {
    0
  } else {
    abs(N1 - N2) / max(N1, N2)
  }

  r <- qr(H)$rank

  list(
    valid = max(gram_error, cross_error, norm_error) <= tol,
    gram_error = gram_error,
    cross_error = cross_error,
    norm_error = norm_error,
    rank = r,
    minimal = nrow(L) == r
  )
}
