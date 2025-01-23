mean_mvnfd <- function(mvnfd_obj) {
  p <- mvnfd_obj$nvar
  mvlist <- lapply(1:p, function(j) mean(mvnfd_obj[, j]))
  return(Mvnfd(mvlist))
}

sd_mvnfd <- function(mvnfd_obj) {
  p <- mvnfd_obj$nvar
  mvlist <- lapply(1:p, function(j) sd(mvnfd_obj[, j]))
  return(Mvnfd(mvlist))
}


#' @title  Addition of two `mvnfd` objects
#'
#' @param obj1 An `mvnfd` object
#' @param obj2 An `mvnfd` object
#' @return An `mvnfd` object
#' @seealso \code{\link{mvnfd}},\code{\link{mvbasismfd}}
#' @export
"+.mvnfd" <- function(obj1, obj2) {
  if (!is.mvnfd(obj1) | !is.mvnfd(obj2)) {
    stop("Both objects must be object of type `mvnfd`")
  }
  if (obj1$nvar != obj2$nvar){
    stop("Both objects must have same number of variables")
  }
  for (i in seq_len(obj1$nvar)) {
    if (all(dim(obj1$data[[i]]) != dim(obj2$data[[i]]))){
      stop("Both objects must have same dimension")
    }
  }
  nfd_list <- lapply(seq_len(obj1$nvar),function(i){
    return(nfd(obj1$data[[i]] + obj2$data[[2]]))
  })
  return(Mvnfd(nfd_list))
}


#' @title  Subtraction of two `mvnfd` objects
#'
#' @param obj1 An `mvnfd` object
#' @param obj2 An optional `mvnfd` object
#' @return An `mvnfd` object
#' @seealso \code{\link{mvnfd}},\code{\link{mvbasismfd}}
#' @export
"-.mvnfd" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
      nfd_obj <- lapply(obj1$data,function(x) {
        data <- (-1)*x
        nfd(data)
        }
      )
    return(Mvnfd(nfd_obj))
  }
  return(obj1 + (-1)*obj2)
}

#'  Multiplication of an `mvnfd` object with a scalar
#'
#'
#' @param obj1 An `mvnfd` object or a scalar
#' @param obj2 An `mvnfd` object or a scalar
#' @return An `mvnfd` object
#' @seealso \code{\link{mvnfd}},\code{\link{nfd}}
#' @export
"*.mvnfd" <- function(obj1, obj2) {
  if (xor(is.mvnfd(obj1), is.mvnfd(obj2))) {
    if (xor(is.double(obj1), is.double(obj2))) {
      if (is.double(obj1)) {
        temp <- obj1
        obj1 <- obj2
        obj2 <- temp
      }
    }
    p <- obj1$nvar
    mvlist <- list()
    for (j in 1:p) {
      mvlist[[j]] <- obj1[,j] * obj2
    }
  } else {
    stop("One object must be an mvnfd, and the other one a scalar")
  }
  return(Mvnfd(mvlist))
}

#'  Extract subsets of an `mvnfd` object
#'
#' @param mvnfd_obj An `mvnfd` object
#' @param i An index or indices specifying the subsets to extract for the first dimension
#' @param j An index or indices specifying the subsets to extract for the second dimension

#' @return An `mvnfd` object containing the specified subsets
#' @seealso \code{\link{mvnfd}},\code{\link{nfd}}
#' @export
"[.mvnfd" <- function(mvnfd_obj, i = NULL, j = NULL) {
  if (is.null(i)) i <- seq_len(mvnfd_obj$nobs)
  if (is.null(j)) j <- seq_len(mvnfd_obj$nvar)
  
  if (any(i > mvnfd_obj$nobs) || any(i < 1)) {
    stop("Subscript 'i' (observations) out of bounds", call. = FALSE)
  }
  
  if (any(j > mvnfd_obj$nvar) || any(j < 1)) {
    stop("Subscript 'j' (variables) out of bounds", call. = FALSE)
  }
  
  selected_data <- mvnfd_obj$data[j]
  
  nfd_objs_list <- lapply(selected_data, function(mat) {
    subset_mat <- mat[i, , drop = FALSE]
    nfd(subset_mat)
  })
  
  if (length(nfd_objs_list) > 1){
    return(Mvnfd(nfd_objs_list))
  } else {
    return(nfd_objs_list[[1]])
  }
  
}

#' @title Compute the inner product between two objects of class `mvnfd`
#'
#' @param mvnfd_obj1 An `mvnfd` object
#' @param mvnfd_obj2 An `mvnfd` object
#' @return The inner products matrix between the two `mvnfd` objects
#' @seealso \code{\link{mvnfd}}
#' @export
inprod_mvnfd <- function(mvnfd_obj1, mvnfd_obj2) {
  if (!is.mvnfd(mvnfd_obj1) || !is.mvnfd(mvnfd_obj2)) stop("Both objects must be of `mvnfd` object type")
  if (mvnfd_obj1$nvar != mvnfd_obj2$nvar) stop("Both objects must have same number of variables")
  return(do.call("+",lapply(seq_along(mvnfd_obj1$data),function(i) inprod_nfd(mvnfd_obj1[,i],mvnfd_obj2[,i]))))
}

#' @title Compute the inner product between two objects of class `mvnfd`
#'
#' @param mvnfd_obj An `mvnfd` object
#' @return The norm of `mvnfd` objects
#' @seealso \code{\link{mvnfd}}
#' @export
norm_mvnfd <- function(mvnfd_obj1) {
  return(sqrt(diag(inprod_mvnfd(mvnfd_obj1,mvnfd_obj1))))
}

center_mvnfd <- function(mvnfd_obj) {
  mmvnfd <- mean(mvnfd_obj)
  nfd_list <- lapply(seq_len(mvnfd_obj$nvar),function(i){
    if (is.matrix(mvnfd_obj$data[[i]])){
      centered_data <- sweep(mvnfd_obj$data[[i]],2,mmvnfd$data[[i]],"-")
      nfd(centered_data)
      } else if (is.array(mvnfd_obj$data[[i]])){
        # I'm not sure about this code check it later for array 
        flattened_current <- apply(mvnfd_obj$data[[i]], 3, as.vector)
        flattened_mean <- apply(mmvnfd$data[[i]], 3, as.vector)
        row_means <- rowMeans(flattened_mean)
        centered_flattened <- sweep(flattened_current, 1, row_means, "-")
        centered_data <- centered_flattened
        nfd(centered_data)
      }
    }
    )
  return(Mvnfd(nfd_list))
}


#' Scale an `mvnfd` Object
#'
#' This function scales an `mvnfd` object by calculating a scaling factor based on the variance of its evaluations 
#' or using a provided weight. It returns a new scaled `mvnfd` object.
#'
#' @param mvnfd_obj An object of class `mvnfd`.
#'@param weight An optional numeric vector of scaling factors for each variable. If NULL, scaling factors are calculated automatically.
#'
#' @return A scaled mvnfd object.
scale_mvnfd <- function(mvnfd_obj,weight = NULL){
  if (!is.null(weight) && length(weight) != mvnfd_obj$nvar) stop("The length of weight vector must be the number of nfd variables")
  nfd_list <- list()
  for (i in 1:mvnfd_obj$nvar){
    nfd_list[[i]] <- scale_nfd(mvnfd_obj[,i],weight = weight[i])
  }
  mvnfd_scaled <- Mvnfd(nfd_list)
  return(mvnfd_scaled)
}
