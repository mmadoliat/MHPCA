#' @title  Penalty Function
#'
#' @description
#' Calculate the penalty matrix for `mvmfd` objects.
#'
#' @param data an object of class `mvmfd`.
#' @param devorder The order of the derivative.
#' @param type The type of penalty. The types "coefpen" and "basispen" is supported.
#' @return The penalty matrix.
#' @importFrom fda eval.penalty
pen_fun <- function(data, devorder = 2, type) {
  init_pen_check(data, devorder, type)
  blocks <- lapply(seq_len(data$nvar), function(i) {
    b <- data$basis$basis[[i]]
    # penalty for each argument of variable i
    P <- lapply(seq_len(b$dimSupp), function(k) {
      if (type == "basispen") {
        eval.penalty(b$basis[[k]], Lfdobj = devorder)
      } else {
        L <- diff(diag(b$nbasis[k]), differences = devorder)
        t(L) %*% L
      }
    })
    if (length(P) == 1) {
      P[[1]]
    } else if (length(P) == 2) {
      P[[2]] %x% diag(nrow(P[[1]])) + diag(nrow(P[[2]])) %x% P[[1]]
    } else {
      stop("Penalties are implemented for 1- and 2-dimensional domains only.", call. = FALSE)
    }
  })
  block_diag(blocks)
}

# a function to check the validity
init_pen_check <- function(data, devorder, type) {
  stopifnot(
    is.mvmfd(data),
    devorder >= 1,
    type == "basispen" | type == "coefpen"
  )
}
