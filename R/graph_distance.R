# Distances between directed graphs.

#' Normalised structural Hamming distance between two directed graphs.
#' Additions and deletions count 1 each; a strict reversal (A has i -> j only,
#' B has j -> i only) counts once instead of twice. The total is divided by the
#' number of edges in A plus B; two empty graphs have distance 0.
nSHD <- function(A, B) {
  stopifnot(
    is.matrix(A), is.matrix(B), identical(dim(A), dim(B)),
    !anyNA(A), !anyNA(B),
    all(A %in% c(0, 1)), all(B %in% c(0, 1)),
    all(diag(A) == 0), all(diag(B) == 0)
  )
  reversals <- sum(A == 1 & t(A) == 0 & B == 0 & t(B) == 1)
  denominator <- sum(A) + sum(B)
  if (denominator == 0) return(0)
  (sum(A != B) - reversals) / denominator
}
