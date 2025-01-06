I_alpha <- function(data, lambda) {
  p <- length(data$coefs)
  I <- c()
  for (i in 1:p) {
    # I is the diagonal matrix with penalty
    I_sub <- lambda[i][[1]] * diag(prod(data$basis$nbasis[[i]]))
    # generate block diagonal matrix
    if (i == 1) {
      I <- I_sub
    } else {
      I <- rbind(cbind(I, matrix(0, nrow = nrow(I), ncol = ncol(I_sub))), cbind(matrix(0, nrow = nrow(I_sub), ncol = ncol(I)), I_sub))
    }
  }
  return(I)
}


# GCVs_rc <- function(data, G, G_half, B_tilde, penalty, p, v, lambda) {
#   s_alpha_tilde <- G_half %*% solve(G + I_alpha(data, lambda) %*% penalty) %*% G_half
#   gcv_score <- (sum(((diag(dim(s_alpha_tilde)[1]) - s_alpha_tilde) %*% (t(B_tilde) %*% v))^2) / ((1 - sum(diag(s_alpha_tilde)) / dim(G)[1])^2)) / dim(G)[1]
#   return(gcv_score)
# }

# a function to check the validity of initializer
init_mfd_check <- function(argval, X, basis, method) {
  stopifnot(
    (is.numeric(argval) | is.list(argval) | is.null(argval)),
    (is.matrix(X) | is.array(X) | is.vector(X)),
    (is.basismfd(basis) | is.basis(basis)),
    (method == "coefs" | method == "data")
  )
}


# A function to check the validity of initializer
init_mfd_list_check <- function(mfd_list) {
  if (!all(sapply(mfd_list, is.mfd))) {
    stop("All the elements of the inputs list must have the class of `mfd`")
  }
  n <- mfd_list[[1]]$nobs
  for (y in mfd_list) {
    if (n != y$nobs) stop("The number of observations in all variables should be equal.")
  }
}

convert_to_nfd <- function(x) {
  if (is.nfd(x)) {
    return(x)  # Already an 'nfd' object
  } else if (is.matrix(x) || is.array(x)) {
    return(nfd(x))  # Convert to 'nfd'
  } else {
    stop("All elements must be either 'nfd', 'matrix', or 'array'.")
  }
}

# A function to check the validity of initializer
init_nfd_list_check <- function(nfd_list) {
  nfd_list <- lapply(nfd_list, convert_to_nfd)
  if (!all(sapply(nfd_list, is.nfd))) {
    stop("All the elements of the inputs list must have the class of `nfd`")
  }
  n <- nrow(nfd_list[[1]])
  for (y in nfd_list) {
    if (n != nrow(y)) stop("The number of observations in all variables should be equal.")
  }
}

# A function to check the validity of initializer
init_hd_list_check <- function(hd_list) {
  mfd_list <- list()
  mvmfd_list <- list()
  nfd_list <- list()
  mvnfd_list <- list()
  for (i in 1:length(hd_list)) {
    if (is.mfd(hd_list[[i]])){
      mfd_list[[length(mfd_list)+1]] <- hd_list[[i]]
    } else if (is.mvmfd(hd_list[[i]])) {
      mvmfd_list[[length(mvmfd_list)+1]] <- hd_list[[i]]
    } else if (is.nfd(hd_list[[i]])){
      nfd_list[[length(nfd_list)+1]] <- hd_list[[i]]
    } else if (is.mvnfd(hd_list[[i]])) {
      mvnfd_list[[length(mvnfd_list)+1]] <- hd_list[[i]]
    }
    else {
      stop("All the elements of the inputs list must have the class of `mfd` , `mvmfd`, `nfd` or `mvnfd`")
    }
  }
  if (length(mfd_list) > 0 && length(mvmfd_list) > 0) {
    stop("You can only provide either 'mfd' or 'mvmfd' objects, not both.")
  }
  if (length(mvmfd_list) > 1) {
    stop("You can only provide one 'mvmfd' objects.")
  }
  n <- hd_list[[1]]$nobs
  for (y in hd_list) {
    if (n != y$nobs) stop("The number of observations in all variables should be equal.")
  }
}
# Function to check the validity of evaluation
eval_mvbasismf_validity_check <- function(evalarg, nvar) {
  if (!is.list(evalarg) & !is.numeric(evalarg)) {
    stop("evalarg must be a list or numeric vector")
  }
  if (is.numeric(evalarg)) {
    if (nvar != 1) {
      stop("evalarg is allowed to be a numeric if nvar = 1.")
    } else {
      evalarg <- list(list(evalarg))
    }
  }
  if (!all(sapply(evalarg, function(x) is.numeric(x) | is.list(x)))) {
    stop("evalarg list elements must be a list or numeric vector")
  }
  if (length(evalarg) != nvar) {
    stop("length of evalarg is not equal to nvar.")
  }
}


# Function to check the validity of initializer
init_mvbasismfd_check <- function(basis) {
  if (is.list(basis)) {
    if (!all(sapply(basis, function(x) {
      return(is.basis(x) | is.basismfd(x))
    }))) {
      stop("All the elements of basis list must be basisfd or basismfd object.")
    }
  }
}


init_basismfd_check <- function(basis) {
  if (!is.basis(basis) & is.list(basis)) {
    if (!all(sapply(basis, is.basis))) {
      stop("All elements of the basis list must be `basisfd` objects.")
    }
  }
}


eval_basismf_validity_check <- function(evalarg, dimSupp) {
  if (!is.list(evalarg) & !is.numeric(evalarg)) {
    stop("evalarg must be a list or numeric vector.")
  }
  if (!all(sapply(evalarg, is.numeric))) {
    stop("evalarg must be a list of numeric vectors.")
  }
  if (is.numeric(evalarg)) {
    evalarg <- list(evalarg)
  }
  if (length(evalarg) != dimSupp) {
    stop("Length of evalarg list must be equal to dimSupp.")
  }
}
