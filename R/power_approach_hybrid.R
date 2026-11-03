#' @importFrom expm sqrtm
#' @importFrom utils  txtProgressBar setTxtProgressBar
#' @importFrom parallel clusterExport parLapplyLB makeCluster clusterEvalQ stopCluster
#' @importFrom foreach getDoParWorkers

sparse_pen_fun <- function(y, tuning_parameter, type, alpha = 3.7, group_sizes = NULL, gamma = 2.5) {
  if (length(tuning_parameter) > 1) {
    n_groups <- length(tuning_parameter)
    # If group_sizes is not provided, assume equal-sized groups.
    if (is.null(group_sizes)) {
      if (length(y) %% n_groups != 0) {
        stop("Length of y is not divisible by number of groups. Please provide group_sizes.")
      }
      group_sizes <- rep(length(y) / n_groups, n_groups)
    }
    
    if (sum(group_sizes) != length(y)) {
      stop("The sum of group_sizes must equal the length of y.")
    }
    
    # Initialize output and compute indices for splitting.
    result <- vector("list", n_groups)
    start <- 1
    for (i in seq_len(n_groups)) {
      end <- start + group_sizes[i] - 1
      y_sub <- y[start:end]
      # Recursive call: note that tuning_parameter[i] is a single number.
      result[[i]] <- sparse_pen_fun(y_sub, tuning_parameter[i], type, alpha)
      start <- end + 1
    }
    return(unlist(result))
  }
  
  y_sorted <- sort(abs(y))
  lambda <- y_sorted[tuning_parameter]
  
  if (tuning_parameter == 0) {
    return(y)
  }
  
  if (type == "soft") {
    return(sign(y) * pmax(abs(y) - lambda, 0))
  } else if (type == "hard") {
    return(ifelse(abs(y) > lambda, y, 0))
  } else if (type == "SCAD") {
    res <- ifelse(abs(y) <= 2 * lambda,
                  sign(y) * pmax(abs(y) - lambda, 0),
                  ifelse(abs(y) <= alpha * lambda,
                         ((alpha - 1) * y - sign(y) * alpha * lambda) / (alpha - 2),
                         y))
    return(res)
  } else if (type == "adaptive") {
    weights <- 1 / (abs(y)^gamma + 10^(-6))
    thr <- lambda * weights
    return(sign(y) * pmax(abs(y) - thr, 0))
    
  } else {
    stop("Unknown penalty type.")
  }
}


# One rank-one sequential power iteration (u = scores, fv / nfv = functional /
# non-functional loadings), used by the final fit and by the CV and eBIC
# tuning. Three modes:
#   fit  (cv_flag = FALSE): loadings are smoothed (S_smooth) and penalized;
#        returns fv (back-transformed and normalized in the smoothing norm),
#        nfv and u.
#   cv_u (cv_flag, pen_u): only u is penalized; returns u.
#   cv_v (cv_flag, pen_nfd / pen_fd): only the loadings are penalized; returns
#        fv (back-transformed) and nfv.
init_sequential_hybrid <- function(fdata,
                                   nfdata,
                                   sparse_tuning_result_u,
                                   sparse_tuning_result_nfd,
                                   sparse_tuning_result_fd,
                                   sparse_tuning_type_u,
                                   sparse_tuning_type_nfd,
                                   sparse_tuning_type_fd,
                                   S_smooth = NULL,
                                   S_2_inverse = NULL,
                                   G_half_inverse = NULL,
                                   G_half = NULL,
                                   cv_flag = FALSE,
                                   pen_u = FALSE,
                                   pen_nfd = FALSE,
                                   pen_fd = FALSE,
                                   tol, max_iter) {
  mode <- if (cv_flag == TRUE && pen_u == TRUE) {
    "cv_u"
  } else if (cv_flag == TRUE && (pen_nfd == TRUE || pen_fd == TRUE)) {
    "cv_v"
  } else if (cv_flag == FALSE) {
    "fit"
  } else {
    return(NULL)  # CV without any penalty: nothing to compute
  }
  has_f <- !is.null(fdata)
  has_nf <- !is.null(nfdata)

  # scores from loadings
  scores <- function(fv, nfv) {
    u <- 0
    if (has_f) u <- u + fdata %*% fv
    if (has_nf) u <- u + nfdata %*% nfv
    u
  }
  # the score update: optional sparsity on u, then unit length
  update_u <- function(u, penalize) {
    if (penalize) {
      u <- sparse_pen_fun(y = u, tuning_parameter = sparse_tuning_result_u, sparse_tuning_type_u)
    }
    u / norm(u, type = "2")
  }
  # loadings from scores; in the fit they are smoothed, and outside cv_u they
  # are penalized (functional loadings in the original, not G_half, scale)
  loadings <- function(u) {
    penalize <- mode != "cv_u"
    fv <- NULL
    nfv <- NULL
    if (has_f) {
      fv <- if (mode == "fit") S_smooth %*% t(fdata) %*% u else t(fdata) %*% u
      if (penalize && pen_fd == TRUE) {
        fv <- G_half %*% sparse_pen_fun(y = G_half_inverse %*% fv, tuning_parameter = sparse_tuning_result_fd,
                                        type = sparse_tuning_type_fd)
      }
    }
    if (has_nf) {
      nfv <- t(nfdata) %*% u
      if (penalize && pen_nfd == TRUE) {
        nfv <- sparse_pen_fun(y = nfv, tuning_parameter = sparse_tuning_result_nfd, type = sparse_tuning_type_nfd)
      }
    }
    # joint unit length
    coef_norm <- sqrt((if (has_f) norm(fv, "2")^2 else 0) + (if (has_nf) norm(nfv, "2")^2 else 0))
    list(fv = if (has_f) fv / coef_norm, nfv = if (has_nf) nfv / coef_norm)
  }

  fv_old <- if (has_f) svd(fdata)$v[, 1] else NULL
  nfv_old <- if (has_nf) svd(nfdata)$v[, 1] else NULL
  errors <- 10^60
  iter <- 0
  while (errors > tol && iter < max_iter) {
    iter <- iter + 1
    u_old <- update_u(scores(fv_old, nfv_old), pen_u == TRUE)
    v <- loadings(u_old)
    u_new <- update_u(scores(v$fv, v$nfv), mode != "cv_v" && pen_u == TRUE)
    errors <- sum((u_new - u_old)^2)
    fv_old <- v$fv
    nfv_old <- v$nfv
  }

  if (mode == "cv_u") return(u_new)

  # Final normalization: functional loadings back to the original scale, with
  # the smoothing norm (fit) or the Euclidean norm (cv_v).
  fv <- v$fv
  nfv <- v$nfv
  coef_norm <- if (has_nf) sqrt(0 + norm(nfv, "2")^2) else 0
  if (has_f) {
    fv <- G_half_inverse %*% fv
    f_norm <- if (mode == "fit") as.numeric(sqrt(t(fv) %*% S_2_inverse %*% fv)) else norm(fv, "2")
    coef_norm <- sqrt(coef_norm^2 + f_norm^2)
    fv <- fv / coef_norm
  }
  if (has_nf) nfv <- nfv / coef_norm
  if (mode == "cv_v") return(list(fv = fv, nfv = nfv))
  list(fv = fv, nfv = nfv, u = u_new)
}

#joint power for smoothing
init_joint_hybrid = function(fdata, nfdata, S_smooth = NULL, S_2_inverse = NULL, G_half_inverse = NULL, G_half = NULL, n = n, tol,max_iter){
  fv_old <- if (!is.null(fdata)) svd(fdata)$v[, 1:n] else NULL
  nfv_old <- if (!is.null(nfdata)) svd(nfdata)$v[, 1:n] else NULL
  
  errors = 10^60
  iter = 0
  while (errors > tol && iter < max_iter) {
    iter <- iter + 1
    u_old <- 0
    if (!is.null(fdata)) {
      u_old <- u_old + fdata %*% fv_old
    }
    if (!is.null(nfdata)) {
      u_old <- u_old + nfdata %*% nfv_old
    }
    u_old = sweep(u_old,2,sqrt(diag(t(u_old)%*%u_old)),"/")
    M1 <- if(!is.null(fdata)) as.matrix(t(fdata)%*%u_old) else NULL
    M2 <- if(!is.null(nfdata)) as.matrix(t(nfdata)%*%u_old) else NULL
    M <- rbind(M1,M2)
    M <- qr.Q(qr(rbind(M)))
    # double check normalizaron 
    if (!is.null(fdata) && !is.null(nfdata)){
      fv_new <- S_smooth%*%M[1:(nrow(M) - ncol(nfdata)),]
      nfv_new <- M[(nrow(M) - (ncol(nfdata)-1)):nrow(M),]
      fv_new_back = G_half_inverse %*% fv_new
      fv_weight <- apply(t(fv_new_back) %*% S_2_inverse %*% fv_new_back,2,norm,"2")^2
      nfv_weight <- apply(nfv_new,2,norm,"2")^2
      coef_norm <- sqrt(fv_weight + nfv_weight)
      fv_new_back <- sweep(fv_new_back,2,coef_norm,"/") 
      fv_new <- G_half %*% fv_new_back
      nfv_new <- sweep(nfv_new,2,coef_norm,"/") 
    } else if (!is.null(fdata) && is.null(nfdata)){
      fv_new <- S_smooth%*%M
      nfv_new <- NULL
    } else if (is.null(fdata) && !is.null(nfdata)){
      nfv_new <- M
      fv_new <- NULL
    }
    u_new <- 0
    if (!is.null(fdata)) {
      u_new <- u_new + fdata %*% fv_new
    }
    if (!is.null(nfdata)) {
      u_new <- u_new + nfdata %*% nfv_new
    }
    u_new = sweep(u_new,2,sqrt(diag(t(u_new)%*%u_new)),"/")

    errors = sum((u_new - u_old)^2)
    fv_old = fv_new
    nfv_old = nfv_new
  }
  # normalizing fv and nfv
  ########################################## 
  if (!is.null(fdata)) {
    fv_new_back = G_half_inverse %*% fv_new
    fv_weight <- apply(t(fv_new_back) %*% S_2_inverse %*% fv_new_back,2,norm,"2")^2
  } else {
    fv_weight <- 0
  }
  if (!is.null(nfdata)){
    nfv_weight <- apply(nfv_new,2,norm,"2")^2
  } else {
    nfv_weight <- 0
  }
  coef_norm <- sqrt(fv_weight + nfv_weight)
  fv_new_back <- if(!is.null(fdata)) sweep(fv_new_back,2,coef_norm,"/") else NULL
  nfv_new <- if(!is.null(nfdata)) sweep(nfv_new,2,coef_norm,"/") else NULL
  ##########################################
  
  return(list(fv_new_back, nfv_new, u_new))
}



cv_local_hybrid <- function(
    fdata,
    nfdata,
    G_half,
    G_half_inverse,
    sparse_tuning_single_u,
    sparse_tuning_single_nfd,
    sparse_tuning_single_fd,
    shuffled_row_u,
    shuffled_row_nfd,
    shuffled_row_fd,
    group_size_u,
    group_size_nfd,
    group_size_fd,
    settings
) {
  nfold_u <- settings$nfold_u
  nfold_nfd <- settings$nfold_nfd
  nfold_fd <- settings$nfold_fd
  sparse_tuning_type_u <- settings$sparse_type_u
  sparse_tuning_type_nfd <- settings$sparse_type_nfd
  sparse_tuning_type_fd <- settings$sparse_type_fd
  pen_u <- settings$pen_u
  pen_nfd <- settings$pen_nfd
  pen_fd <- settings$pen_fd
  tol <- settings$tol
  max_iter <- settings$max_iter

  # Precompute fdata in transformed form (if fdata is not NULL)
  data_double_tilde <- if (!is.null(fdata)) t(fdata %*% G_half) else NULL

  # Extract row indices from shuffled_row_u (for u-penalty CV)
  shuffled_row_f  <- shuffled_row_u$f
  shuffled_row_nf <- shuffled_row_u$nf

  # Extract group sizes from group_size_u (for u-penalty CV)
  group_size_f  <- group_size_u$f
  group_size_nf <- group_size_u$nf

  # Initialize accumulators for CV errors
  error_score_sparse_u   <- 0
  error_score_sparse_nfd <- 0
  error_score_sparse_fd <- 0
  #-------------------------------
  # 1) Cross-validation over 'u'
  #-------------------------------
  err_u <- NULL
  if (pen_u) {
    fold_errors_u <- numeric(nfold_u)
    for (k in seq_len(nfold_u)) {
      rows_to_remove_f <- if (!is.null(fdata)) {
        shuffled_row_f[((k-1)*group_size_f + 1) : (k * group_size_f)]
      } else NULL

      rows_to_remove_nf <- if (!is.null(nfdata)) {
        shuffled_row_nf[((k-1)*group_size_nf + 1) : (k * group_size_nf)]
      } else NULL

      # Training/testing splits
      if (!is.null(fdata)) {
        fdata_train <- data_double_tilde[-rows_to_remove_f, , drop = FALSE]
        fdata_test  <- data_double_tilde[ rows_to_remove_f, , drop = FALSE]

      } else {
        fdata_train <- NULL
        fdata_test  <- NULL

      }

      if (!is.null(nfdata)) {
        nfdata_train <- nfdata[, -rows_to_remove_nf, drop = FALSE]
        nfdata_test  <- nfdata[,  rows_to_remove_nf, drop = FALSE]
      } else {
        nfdata_train <- NULL
        nfdata_test  <- NULL
      }

      # Fit model with current penalty
      
      u_test <- init_sequential_hybrid(
        fdata                 = if (!is.null(fdata_train)) t(fdata_train) else NULL,
        nfdata                = nfdata_train,
        sparse_tuning_result_u= sparse_tuning_single_u,
        sparse_tuning_type_u  = sparse_tuning_type_u,
        cv_flag               = TRUE,
        pen_nfd          = FALSE,
        pen_u            = pen_u,
        tol = tol, max_iter = max_iter
      )

      # Compute error on the held-out fold
      if (!is.null(fdata)) {
        fv_test_smooth_back    <- (data_double_tilde %*% u_test)[rows_to_remove_f, , drop = FALSE]
        fdata_test_smooth_back <- t(data_double_tilde)[, rows_to_remove_f, drop = FALSE]
        fv_error <- sum((t(fdata_test_smooth_back) - fv_test_smooth_back %*% t(u_test))^2)
      } else {
        fv_error <- 0
      }

      if (!is.null(nfdata)) {
        nfv_test  <- t(nfdata_test) %*% u_test
        nfv_error <- sum((t(nfdata_test) - nfv_test %*% t(u_test))^2)
      } else {
        nfv_error <- 0
      }

      error_score_sparse_u <- error_score_sparse_u + (fv_error + nfv_error)
      fold_errors_u[k] <- fv_error + nfv_error
    }

    # Normalization
    normalization_factor <- if (!is.null(fdata)) {
      ncol(fdata)
    } else if (!is.null(nfdata)) {
      nrow(nfdata)
    } else {
      1
    }
    err_u <- error_score_sparse_u / normalization_factor
    errs_u <- fold_errors_u/normalization_factor
    mean_u <- mean(errs_u)
    sd_u <- sd(errs_u)
    se_u <- sd_u/sqrt(nfold_u)
  } else {
    mean_u <- sd_u <- se_u <- NULL
  }

  #----------------------------------------------
  # 2) Cross-validation over 'nfd' terms
  #----------------------------------------------

  err_nfd <- NULL

  if (pen_nfd) {
    fold_errors_nfd <- numeric(nfold_nfd)
    
    for (k in seq_len(nfold_nfd)) {
      rows_to_remove_nfd <- shuffled_row_nfd[((k-1)*group_size_nfd + 1) : (k * group_size_nfd)]
      nfdata_train <- nfdata[-rows_to_remove_nfd, , drop = FALSE]
      nfdata_test  <- nfdata[ rows_to_remove_nfd, , drop = FALSE]
      

      # Fit model with current penalty
      v_test <- init_sequential_hybrid(
        fdata                 = NULL, 
        nfdata                = nfdata_train,
        sparse_tuning_result_nfd = sparse_tuning_single_nfd,
        sparse_tuning_type_nfd   = sparse_tuning_type_nfd,
        G_half_inverse = G_half_inverse,
        G_half = G_half,
        cv_flag               = TRUE,
        pen_nfd          = pen_nfd,
        pen_fd           = FALSE,
        pen_u            = FALSE,
        tol = tol, max_iter = max_iter
      )


      nfv_test <- v_test$nfv

      # Combine contributions into u_test
      u_test <- 0
      
      u_test <- u_test + nfdata_test %*% nfv_test


      nfv_error <- sum((t(nfdata_test) - nfv_test %*% t(u_test))^2)


      error_score_sparse_nfd <- error_score_sparse_nfd + nfv_error # (fv_error + nfv_error)
      fold_errors_nfd[k] <-  nfv_error

    }

    # Normalization
    normalization_factor <- nrow(nfdata)
    err_nfd <- error_score_sparse_nfd / normalization_factor
    errs_nfd <- fold_errors_nfd/normalization_factor
    mean_nfd <- mean(errs_nfd)
    sd_nfd <- sd(errs_nfd)
    se_nfd <- sd_nfd/sqrt(nfold_nfd)

  } else {
    mean_nfd <- sd_nfd <- se_nfd <- NULL
  }

  #----------------------------------------------
  # 3) Cross-validation over 'fd' terms
  #----------------------------------------------

  err_fd <- NULL

  if (pen_fd) {
    fold_errors_fd <- numeric(nfold_fd)

    for (k in seq_len(nfold_fd)) {
      rows_to_remove_fd <- shuffled_row_fd[((k-1)*group_size_fd + 1) : (k * group_size_fd)]
      # Training/testing splits
      fdata_train <- data_double_tilde[, -rows_to_remove_fd, drop = FALSE]
      fdata_test  <- data_double_tilde[,  rows_to_remove_fd, drop = FALSE]

      # Fit model with current penalty
      v_test <- init_sequential_hybrid(
        fdata                 = t(fdata_train), 
        nfdata                = NULL, #nfdata_train,
        sparse_tuning_result_fd  = sparse_tuning_single_fd,
        sparse_tuning_type_fd    = sparse_tuning_type_fd,
        G_half_inverse = G_half_inverse,
        G_half = G_half,
        cv_flag               = TRUE,
        pen_nfd          = FALSE,
        pen_fd           = pen_fd,
        pen_u            = FALSE,
        tol = tol, max_iter = max_iter
      )

      fv_test_smooth_back <- v_test$fv

      # Combine contributions into u_test
      u_test <- 0
      u_test <- u_test + t(fdata_test) %*% fv_test_smooth_back

      # Compute errors for this fold

      fdata_test_smooth_back <- t(data_double_tilde)[rows_to_remove_fd, , drop = FALSE]
      fv_error <- sum((t(fdata_test_smooth_back) - fv_test_smooth_back %*% t(u_test))^2)

      error_score_sparse_fd <- error_score_sparse_fd + fv_error 
      fold_errors_fd[k] <-  fv_error

    }

    # Normalization
    normalization_factor <- ncol(fdata)
    err_fd <- error_score_sparse_fd / normalization_factor
    errs_fd <- fold_errors_fd/normalization_factor
    mean_fd <- mean(errs_fd)
    sd_fd <- sd(errs_fd)
    se_fd <- sd_fd/sqrt(nfold_fd)

  } else {
    mean_fd <- sd_fd <- se_fd <- NULL
  }

  #--------------------------------
  # Return the same output structure
  #--------------------------------
  list(err_u = err_u, err_nfd = err_nfd,err_fd = err_fd, se_u = se_u, se_nfd = se_nfd, se_fd = se_fd)
}


#computing gcv score for smoothing tuning
gcv_local_hybrid = function(fdata, nfdata, hd_obj, G, G_half, S_smooth, u, smooth_tuning) {
  mvmfd_obj <-  hd_obj$mf 
  p = mvmfd_obj$nvar
  indices <- sapply(1:p, function(i) prod(mvmfd_obj$basis$nbasis[[i]]))
  C_subtilde = fdata %*% G_half
  # print(dim(u))
  if (all(smooth_tuning == 0)) {
    error_smooth_score <- 0
  } else {
    error_smooth_score <- 0
    start_index <- 1
    
    for (i in 1:p) {
      end_index <- start_index + indices[i] - 1
      
      if (smooth_tuning[i] == 0) {
        error_smooth_score_i <- 0
      } else {
        s_alpha_tilde_i <- S_smooth[start_index:end_index, start_index:end_index]
        C_subtilde_i <- C_subtilde[, start_index:end_index]
        error_smooth_score_i <- sum(((diag(indices[i]) - s_alpha_tilde_i) %*% (t(C_subtilde_i) %*% u))^2) /  
          (1 - sum(diag(s_alpha_tilde_i)) / indices[i])^2
      }
      
      error_smooth_score <- error_smooth_score + error_smooth_score_i
      start_index <- end_index + 1
    }
  }
  
  return(error_smooth_score)
}

# Function to handle smooth tuning selection with progress bar and GCV score calculation
handle_smooth_tuning_hybrid <- function(fdata, 
                                        nfdata, 
                                        G_half, 
                                        G, 
                                        S_smooth, 
                                        S_2_inverse, 
                                        G_half_inverse, 
                                        hd_obj, 
                                        sparse_tuning_selection_u = NULL, 
                                        sparse_tuning_selection_nfd = NULL,
                                        sparse_tuning_selection_fd = NULL,
                                        smooth_tuning,
                                        CV_score_smooth,
                                        power_type = "sequential",
                                        n = NULL,
                                        pb,
                                        count,
                                        pen_nfd = F,
                                        pen_fd = F,
                                        pen_u = F,
                                        settings) {
  # pen_* are not taken from settings: the smoothing GCV is run without the
  # sparsity penalties unless a caller asks for them.
  sparse_tuning_type_u <- settings$sparse_type_u
  sparse_tuning_type_nfd <- settings$sparse_type_nfd
  sparse_tuning_type_fd <- settings$sparse_type_fd
  tol <- settings$tol
  max_iter <- settings$max_iter
  if (!is.null(hd_obj$mf)){
    mvmfd_obj <- hd_obj$mf
    gcv_scores <- NULL
    smooth_tuning_selection <- NULL
    index_selection <- NULL
    if (is.null(smooth_tuning)) {
      count <- count + 1
      setTxtProgressBar(pb, count)
      smooth_tuning_selection <- expand.grid(lapply(rep(0, mvmfd_obj$nvar), function(x) x[1]))
      index_selection <- 1
    } else {
      for (smooth_index in 1:dim(smooth_tuning)[1]) {
        count <- count + 1
        setTxtProgressBar(pb, count)
        if (all(smooth_tuning == 0)) {
          S_smooth[[smooth_index]] <- diag(dim(G)[1])
        }
        if (power_type == "sequential") {
          test_temp <- init_sequential_hybrid(fdata = fdata %*% G_half, 
                                              nfdata = nfdata, 
                                              sparse_tuning_result_u = sparse_tuning_selection_u,
                                              sparse_tuning_result_nfd = sparse_tuning_selection_nfd, 
                                              sparse_tuning_result_fd = sparse_tuning_selection_fd,
                                              sparse_tuning_type_u = sparse_tuning_type_u,
                                              sparse_tuning_type_nfd = sparse_tuning_type_nfd,
                                              sparse_tuning_type_fd = sparse_tuning_type_fd,
                                              S_smooth = S_smooth[[smooth_index]], 
                                              S_2_inverse = S_2_inverse[[smooth_index]],
                                              G_half_inverse =  G_half_inverse,
                                              G_half = G_half,
                                              pen_nfd = pen_nfd,
                                              pen_fd = pen_fd,
                                              pen_u = pen_u,
                                              tol = tol, max_iter = max_iter)
        } else {
          test_temp <- init_joint_hybrid(fdata %*% G_half, nfdata ,S_smooth[[smooth_index]], S_2_inverse[[smooth_index]], G_half_inverse, G_half, n = n,tol = tol, max_iter = max_iter)
        }
        u_temp <- test_temp[[3]]
        smooth_score <- gcv_local_hybrid(fdata, nfdata, hd_obj, G, G_half, S_smooth[[smooth_index]], u_temp, smooth_tuning = smooth_tuning[smooth_index, ])
        gcv_scores <- c(gcv_scores, smooth_score)
        if (smooth_score <= CV_score_smooth) {
          CV_score_smooth <- smooth_score
          smooth_tuning_selection <- smooth_tuning[smooth_index, ]
          index_selection <- smooth_index
        }
      }
    }
    return(list(smooth_tuning_selection = smooth_tuning_selection, index_selection = index_selection, gcv_scores = gcv_scores))
  } else {
    return(NULL)
  }
  
  
}


handle_sparse_tuning_hybrid <- function(
    fdata,
    nfdata,
    G_half,
    G_half_inverse,
    sparse_tuning_u,
    sparse_tuning_nfd,
    sparse_tuning_fd,
    shuffled_row_u,
    shuffled_row_nfd,
    shuffled_row_fd,
    group_size_u,
    group_size_nfd,
    group_size_fd,
    CV_score_sparse_u,
    CV_score_sparse_nfd,
    CV_score_sparse_fd,
    pb,
    cl,   # cluster passed in from outside (NULL for sequential)
    settings
) {
  sparse_tuning_type_u <- settings$sparse_type_u
  sparse_tuning_type_nfd <- settings$sparse_type_nfd
  sparse_tuning_type_fd <- settings$sparse_type_fd
  nfold_u <- settings$nfold_u
  nfold_nfd <- settings$nfold_nfd
  nfold_fd <- settings$nfold_fd
  pen_u <- settings$pen_u
  pen_nfd <- settings$pen_nfd
  pen_fd <- settings$pen_fd
  sparse_iter <- settings$sparse_iter
  tol <- settings$tol
  max_iter <- settings$max_iter
  cv.pick <- settings$cv.pick

  applyFun <- if (is.null(cl)) lapply else function(x, FUN) parallel::parLapplyLB(cl, x, FUN)
  # Collect the names of the objects to export to the workers:
  data_vars <- c(
    "fdata", "nfdata",
    "G_half", "G_half_inverse",
    "nfold_u", "nfold_nfd", "nfold_fd",
    "sparse_tuning_type_u", "sparse_tuning_type_nfd", "sparse_tuning_type_fd",
    "shuffled_row_u", "shuffled_row_nfd", "shuffled_row_fd",
    "group_size_u", "group_size_nfd", "group_size_fd","tol","max_iter",
    "pen_u","pen_fd","pen_nfd","settings"
  )
  
  # Quick exit if nothing to tune:
  if (is.null(sparse_tuning_u) &&
      is.null(sparse_tuning_nfd) &&
      is.null(sparse_tuning_fd)) {
    setTxtProgressBar(pb, 4)
    return(list(
      sparse_tuning_selection_u   = 0,
      sparse_tuning_selection_nfd = NULL,
      sparse_tuning_selection_fd  = NULL,
      cv_scores_u                 = NULL,
      cv_scores_nfd               = NULL,
      cv_scores_fd                = NULL,
      CV_score_sparse_u           = CV_score_sparse_u,
      CV_score_sparse_nfd         = CV_score_sparse_nfd,
      CV_score_sparse_fd          = CV_score_sparse_fd
    ))
  }
  # Cross-validate the candidates `cands` of block `block` (variable j), with
  # the other blocks held at their current selections.
  evaluate_cv <- function(block, j, cands, current) {
    if (!is.null(cl)) {
      parallel::clusterExport(cl, varlist = c(data_vars, "current", "j"), envir = environment())
    }
    rows <- applyFun(cands, function(candidate) {
      trial <- with_candidate(current, block, j, candidate)
      res <- cv_local_hybrid(
        fdata   = fdata,
        nfdata  = nfdata,
        G_half  = G_half,
        G_half_inverse = G_half_inverse,
        sparse_tuning_single_u   = trial$u,
        sparse_tuning_single_nfd = trial$nfd,
        sparse_tuning_single_fd  = trial$fd,
        shuffled_row_u   = shuffled_row_u,
        shuffled_row_nfd = shuffled_row_nfd,
        shuffled_row_fd  = shuffled_row_fd,
        group_size_u   = group_size_u,
        group_size_nfd = group_size_nfd,
        group_size_fd  = group_size_fd,
        settings = settings
      )
      c(err = res[[paste0("err_", block)]], se = res[[paste0("se_", block)]])
    })
    m <- do.call(rbind, rows)
    list(score = m[, "err"], se = m[, "se"])
  }

  tuned <- tune_sparse_blocks(
    candidates = list(u = sparse_tuning_u, nfd = sparse_tuning_nfd, fd = sparse_tuning_fd),
    pen = list(u = pen_u, nfd = pen_nfd, fd = pen_fd),
    method = "cv",
    evaluate = evaluate_cv,
    pick = function(score, se) pick_cv_index(score, se, cv.pick),
    pb = pb, count = 0, passes = sparse_iter
  )
  # CV_score_sparse_*: score of the last selected candidate (the last
  # variable for nfd / fd); the value passed in if the block was not tuned.
  last_or <- function(block, default) {
    if (is.null(tuned$last_score[[block]])) default else tuned$last_score[[block]]
  }
  list(
    sparse_tuning_selection_u   = tuned$selection$u,
    sparse_tuning_selection_nfd = tuned$selection$nfd,
    sparse_tuning_selection_fd  = tuned$selection$fd,
    cv_scores_u                 = tuned$scores$u,
    cv_scores_nfd               = tuned$scores$nfd,
    cv_scores_fd                = tuned$scores$fd,
    CV_score_sparse_u           = last_or("u", CV_score_sparse_u),
    CV_score_sparse_nfd         = last_or("nfd", CV_score_sparse_nfd),
    CV_score_sparse_fd          = last_or("fd", CV_score_sparse_fd),
    cv_se_u                     = tuned$ses$u,
    cv_se_nfd                   = tuned$ses$nfd,
    cv_se_fd                    = tuned$ses$fd
  )
}



cv_gcv_sequential_hybrid <- function(
    fdata,
    nfdata,
    hd_obj,
    smooth_tuning,
    sparse_tuning_u,
    sparse_tuning_nfd,
    sparse_tuning_fd,
    G,
    G_half,
    G_half_inverse,
    S_smooth,
    S_2_inverse,
    cl,
    settings
) {
  # Validate sparse_selection_method
  sparse_selection_method <- match.arg(settings$sparse_selection_method, c("cv", "ebic"))
  #-----------------------------------------
  # Branch based on selection method
  #-----------------------------------------
  n_iter <- (if (is.null(smooth_tuning)) 1 else nrow(smooth_tuning)) +
    (if (is.null(sparse_tuning_u)) 1 else length(sparse_tuning_u)) +
    (if (is.null(sparse_tuning_nfd)) 1 else sum(sapply(sparse_tuning_nfd, length)))+
    (if (is.null(sparse_tuning_fd)) 1 else  sum(sapply(sparse_tuning_fd, length)))
  pb <- txtProgressBar(min = 0, max =  n_iter, style = 3, width = 50, char = "=")
  if (sparse_selection_method == "ebic") {
    
    # Use eBIC for sparse parameter selection
    message("Using eBIC for sparse tuning parameter selection...")

    # Smoothing used while tuning sparsity. With a single smoothing value its
    # matrix is known. With a GCV grid, smoothing is selected only after
    # sparsity, so sparsity is tuned without smoothing (lambda = 0), as the CV
    # path does, rather than with an arbitrary grid point.
    if (!is.null(fdata) && length(S_smooth) > 1) {
      S_2_none <- solve(G)
      S_smooth_ebic <- G_half %*% S_2_none %*% G_half
      S_2_inverse_ebic <- solve(S_2_none)
    } else {
      S_smooth_ebic <- S_smooth
      S_2_inverse_ebic <- S_2_inverse
    }

    sparse_res <- eBIC_sequential_hybrid(
      # eBIC fits in the same transformed coordinates as the final fit
      fdata = if (!is.null(fdata)) fdata %*% G_half else NULL,
      nfdata = nfdata,
      hd_obj = hd_obj,
      sparse_tuning_u = sparse_tuning_u,
      sparse_tuning_nfd = sparse_tuning_nfd,
      sparse_tuning_fd = sparse_tuning_fd,
      sparse_tuning_type_u = settings$sparse_type_u,
      sparse_tuning_type_nfd = settings$sparse_type_nfd,
      sparse_tuning_type_fd = settings$sparse_type_fd,
      G_half = G_half,
      G_half_inverse = G_half_inverse,
      S_smooth = S_smooth_ebic,
      S_2_inverse = S_2_inverse_ebic,
      pen_u = settings$pen_u,
      pen_nfd = settings$pen_nfd,
      pen_fd = settings$pen_fd,
      xi = settings$ebic_xi,
      N_type = settings$ebic_N_type,
      RSS_type = settings$ebic_RSS_type,
      tol = settings$tol,
      max_iter = settings$max_iter,
      verbose = TRUE
    )
    
    sparse_tuning_selection_u <- sparse_res$sparse_tuning_selection_u
    sparse_tuning_selection_nfd <- sparse_res$sparse_tuning_selection_nfd
    sparse_tuning_selection_fd <- sparse_res$sparse_tuning_selection_fd
    
    # For eBIC, we return eBIC scores instead of CV scores
    cv_scores_u <- NULL
    cv_scores_nfd <- NULL
    cv_scores_fd <- NULL
    cv_se_u <- NULL
    cv_se_nfd <- NULL
    cv_se_fd <- NULL
    
    eBIC_scores_u <- sparse_res$eBIC_scores_u
    eBIC_scores_nfd <- sparse_res$eBIC_scores_nfd
    eBIC_scores_fd <- sparse_res$eBIC_scores_fd
    
  } else {
  # Prepare fold assignments
  ncf   <- if (!is.null(fdata))  ncol(fdata)  else NULL
  ncnf  <- if (!is.null(nfdata)) ncol(nfdata) else NULL
  shuffled_row_u   <- list(f = sample(ncf), nf = sample(ncnf))
  group_size_u     <- list(f = length(shuffled_row_u$f)/settings$nfold_u,
                           nf= length(shuffled_row_u$nf)/settings$nfold_u)
  nrf   <- if (!is.null(fdata))  nrow(fdata)  else nrow(nfdata)
  shuffled_row_nfd <- shuffled_row_fd <- sample(nrf)
  group_size_nfd   <- length(shuffled_row_nfd)/settings$nfold_nfd
  group_size_fd    <- length(shuffled_row_nfd)/settings$nfold_fd
  # Determine progress-bar length
  # n_iter <- (if (is.null(smooth_tuning)) 1 else nrow(smooth_tuning)) +
  #   (if (is.null(sparse_tuning_u)) 1 else length(sparse_tuning_u)) +
  #   (if (is.null(sparse_tuning_nfd)) 1 else sum(sapply(sparse_tuning_nfd, length)))+
  #   (if (is.null(sparse_tuning_fd)) 1 else  sum(sapply(sparse_tuning_fd, length)))
  # pb <- txtProgressBar(min = 0, max =  n_iter, style = 3, width = 50, char = "=")

  # 1) Sparse tuning (parallel)
  sparse_res <- handle_sparse_tuning_hybrid(
    fdata   = fdata,      nfdata = nfdata,
    G_half  = G_half,     G_half_inverse = G_half_inverse,
    sparse_tuning_u   = sparse_tuning_u,
    sparse_tuning_nfd = sparse_tuning_nfd,
    sparse_tuning_fd  = sparse_tuning_fd,
    shuffled_row_u   = shuffled_row_u,
    shuffled_row_nfd = shuffled_row_nfd,
    shuffled_row_fd  = shuffled_row_fd,
    group_size_u   = group_size_u,
    group_size_nfd = group_size_nfd,
    group_size_fd  = group_size_fd,
    CV_score_sparse_u   = Inf,
    CV_score_sparse_nfd = Inf,
    CV_score_sparse_fd  = Inf,
    pb          = pb,
    cl           = cl,
    settings = settings
  )
  sparse_tuning_selection_u   <- sparse_res$sparse_tuning_selection_u
  sparse_tuning_selection_nfd <- sparse_res$sparse_tuning_selection_nfd
  sparse_tuning_selection_fd  <- sparse_res$sparse_tuning_selection_fd
  cv_scores_u                 <- sparse_res$cv_scores_u
  cv_scores_nfd               <- sparse_res$cv_scores_nfd
  cv_scores_fd                <- sparse_res$cv_scores_fd
  cv_se_u <- sparse_res$cv_se_u
  cv_se_nfd <- sparse_res$cv_se_nfd
  cv_se_fd <- sparse_res$cv_se_fd
  
  eBIC_scores_u <- NULL
  eBIC_scores_nfd <- NULL
  eBIC_scores_fd <- NULL
  
  # close(pb)
  
  }
  
  #pb_smooth <- txtProgressBar(min = 0, max = n_smooth, style = 3, width = 50, char = "=")
  
  # 2) Smooth tuning (GCV)
  count0 <- (if (is.null(sparse_tuning_u)) 1 else length(sparse_tuning_u)) +
    (if (is.null(sparse_tuning_nfd)) 1 else sum(sapply(sparse_tuning_nfd, length))) +
    (if (is.null(sparse_tuning_fd)) 1 else sum(sapply(sparse_tuning_fd, length)))
  smooth_res <- handle_smooth_tuning_hybrid(
    fdata                     = fdata,
    nfdata                    = nfdata,
    G_half                    = G_half,
    G                         = G,
    S_smooth                  = S_smooth,
    S_2_inverse               = S_2_inverse,
    G_half_inverse            = G_half_inverse,
    hd_obj                    = hd_obj,
    sparse_tuning_selection_u   = sparse_tuning_selection_u,
    sparse_tuning_selection_nfd = sparse_tuning_selection_nfd,
    sparse_tuning_selection_fd  = sparse_tuning_selection_fd,
    smooth_tuning             = smooth_tuning,
    CV_score_smooth           = Inf,
    power_type                = "sequential",
    pb                         = pb,
    count                      = count0,
    settings = settings
  )
  #close(pb_smooth)
  close(pb)
  smooth_tuning_selection <- smooth_res$smooth_tuning_selection
  index_selection         <- smooth_res$index_selection
  gcv_scores              <- smooth_res$gcv_scores
  
  list(
   sparse_tuning_selection_u = sparse_tuning_selection_u,
   sparse_tuning_selection_nfd = sparse_tuning_selection_nfd,
   sparse_tuning_selection_fd = sparse_tuning_selection_fd,
   eBIC_scores_u = eBIC_scores_u,
   eBIC_scores_nfd = eBIC_scores_nfd,
   eBIC_scores_fd = eBIC_scores_fd,
   cv_scores_u = cv_scores_u,
   cv_scores_nfd = cv_scores_nfd,
   cv_scores_fd = cv_scores_fd,
   gcv_scores = gcv_scores,
   index_selection = index_selection,
   cv_se_u = cv_se_u,
   cv_se_nfd = cv_se_nfd,
   cv_se_fd = cv_se_fd,
   selection_method = sparse_selection_method
 )
}

# Function for gcv_joint
gcv_joint_hybrid <- function(fdata, nfdata, G_half, G, S_smooth, S_2_inverse, G_half_inverse, hd_obj, smooth_tuning, n, settings) {
  mvmfd_obj <- hd_obj$mf
  CV_score_smooth <- Inf
  pb <- txtProgressBar(min = 0, max = if (is.null(smooth_tuning)) 1 else dim(smooth_tuning)[1], style = 3, width = 50, char = "=")
  smooth_tuning_result <- handle_smooth_tuning_hybrid(fdata, nfdata, G_half, G, S_smooth, S_2_inverse, G_half_inverse, hd_obj, smooth_tuning = smooth_tuning, CV_score_smooth = CV_score_smooth, power_type = "joint", n = n, pb = pb, count = 0, settings = settings)
  smooth_tuning_selection <- smooth_tuning_result$smooth_tuning_selection
  index_selection <- smooth_tuning_result$index_selection
  gcv_scores <- smooth_tuning_result$gcv_scores
  
  close(pb)
  return(list(smooth_tuning_selection = smooth_tuning_selection, index_selection = index_selection, gcv_scores = gcv_scores))
}

# Function to handle variance calculation and update
handle_variance_update_hybrid <- function(i, n, C, nf_data, G, fv_total, nfv_total,hd_obj, all_equal_check, 
                                   sparse_tuning, pc, lsv, fv, nfv, u, G_half, test_result, temp_count, fp) {
  
  mvmfd_obj <- hd_obj$mf
  if (!is.null(mvmfd_obj)){
    for (j in 1:fp) {
      index_start <- (temp_count + 1)
      index_end <- (temp_count + prod(mvmfd_obj$basis$nbasis[[j]]))
      if (i == 1) {
        pc[[j]] <- fv[index_start:index_end, ]
      } else {
        pc[[j]] <- cbind(pc[[j]], fv[index_start:index_end, ])
      }
      temp_count <- index_end
    }
  }
  
  lsv = cbind(lsv, u)
  fv_total = cbind(fv_total, fv)
  nfv_total = cbind(nfv_total, nfv)
  ### correct it later
  if (i == 1 || all(all_equal_check)  & (is.null(sparse_tuning) || all(unique(unlist(sparse_tuning)) == 0))) {
    CGfv <- if (!is.null(mvmfd_obj)) C %*% G %*% fv else NULL
    Xnfv <- if (!is.null(hd_obj$nf)) nf_data %*% nfv else NULL
    if (!is.null(CGfv) && !is.null(Xnfv)) variance <- (t(Xnfv)%*%Xnfv+2*t(Xnfv)%*% CGfv+t(CGfv)%*% CGfv)/(mvmfd_obj$nobs - 1)
    if (!is.null(CGfv) && is.null(Xnfv)) variance <- (t(CGfv)%*% CGfv)/(mvmfd_obj$nobs - 1)
    if (is.null(CGfv) && !is.null(Xnfv)) variance <- (t(Xnfv)%*%Xnfv)/(hd_obj$nf$nobs - 1)
  } else {
    if (!is.null(mvmfd_obj)){
      
      CGfv <- C %*% G %*% fv_total
      fvTGfv <- t(fv_total) %*% G %*% fv_total
      CGfvp <- C %*% G %*% fv_total[,-i]
      fvTGfvp <- t(fv_total[,-i]) %*% G %*% fv_total[,-i]
    } else {
      CGfv <- fvTGfv <- CGfvp <-  fvTGfvp <- 0
    }
    if (!is.null(hd_obj$nf)){
      Xnfv <- nf_data %*% nfv_total 
      nfvTnfv <- t(nfv_total) %*% nfv_total
      Xnfvp <- nf_data %*% nfv_total[,-i] 
      nfvTnfvp <- t(nfv_total[,-i]) %*% nfv_total[,-i]
    } else {
      Xnfv <- nfvTnfv <- Xnfvp <- nfvTnfvp <- 0
    }
    
    G_pc_inv <- solve(nfvTnfv + fvTGfv)
    G_pc_invp <- solve(nfvTnfvp + fvTGfvp)
    total_variance <- sum(diag((Xnfv + CGfv)%*%G_pc_inv%*%t(Xnfv + CGfv))) 
    total_variance_p <- sum(diag((Xnfvp + CGfvp)%*%G_pc_invp%*%t(Xnfvp + CGfvp))) 
    nn <- if(!is.null(mvmfd_obj)) (mvmfd_obj$nobs - 1) else (hd_obj$nf$nobs - 1)
    variance = (total_variance - total_variance_p) / nn
  }
  return(list(pc = pc, lsv = lsv, fv_total = fv_total, nfv_total = nfv_total, variance = variance))
}

# sequential power algorithm
sequential_power_hybrid <- function(hd_obj,
                                    n,
                                    smooth_tuning,
                                    sparse_tuning_u,
                                    sparse_tuning_nfd,
                                    sparse_tuning_fd,
                                    settings) {
  smooth_tuning_type <- settings$smoothing_type
  centerfns <- settings$centerfns
  sparse_CV <- settings$sparse_CV
  smooth_GCV <- settings$smooth_GCV
  sparse_tuning_type_u <- settings$sparse_type_u
  sparse_tuning_type_nfd <- settings$sparse_type_nfd
  sparse_tuning_type_fd <- settings$sparse_type_fd
  pen_u <- settings$pen_u
  pen_nfd <- settings$pen_nfd
  pen_fd <- settings$pen_fd
  tol <- settings$tol
  max_iter <- settings$max_iter
  sparse_selection_method <- settings$sparse_selection_method
  ncor <- settings$ncor

  use_cluster <- ncor > 1L
  if (use_cluster){
    cl <- parallel::makeCluster(
      ncor,
      outfile = if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
    )
    on.exit({
      parallel::stopCluster(cl)
      #closeAllConnections()
    }, add = TRUE)
    
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(
        library(MHPCA, quietly = TRUE, warn.conflicts = FALSE)
      )
    })
  } else {
    cl <- NULL
  }
 

  
  #######centralize########
  if (centerfns) hd_obj <- center_hd(hd_obj)
  mvmfd_obj <- hd_obj$mf
  nf_obj <- hd_obj$nf
  if (!is.null(nf_obj)) {
    nf_data <- do.call("cbind",nf_obj$data)
  } else {
    nf_data <- NULL
  }
  if (!is.null(mvmfd_obj)){
    C <- do.call("rbind", lapply(mvmfd_obj$coefs, coef_matrix)) 
    p <- mvmfd_obj$nvar 
    smooth_penalty <-  pen_fun(mvmfd_obj, type = smooth_tuning_type) 
    C <- t(C)
    lsv <- c()
    pc <- list()
    variance <- vector()
    smooth_tuning_result  <- list()
    sparse_tuning_result_u  <- list()
    sparse_tuning_result_nfd  <- list()
    sparse_tuning_result_fd  <- list()
    G <- as.matrix(mvmfd_obj$basis$gram)
    G_half <- sqrtm(G)
    G_half_inverse = solve(G_half)
    all_equal_check <- sapply(smooth_tuning, function(x) length(unique(x)) == 1)
  } else {
    C <- NULL 
    p <- NULL
    smooth_penalty <-  NULL
    lsv <- c()
    pc <- list()
    variance <- vector()
    smooth_tuning_result  <- list()
    sparse_tuning_result_u <- list()
    sparse_tuning_result_nfd <- list()
    sparse_tuning_result_fd <- list()
    G <- NULL
    G_half <- NULL
    G_half_inverse = NULL
    all_equal_check <- NULL
  }
  ######### smoothing parameters and the loop over components #########
  # smooth_GCV = TRUE ("grid"): `smooth_tuning` is a grid searched by GCV. The
  #   smoothing matrices for the whole grid are computed once, up front, and
  #   the one selected by GCV is used for the fit.
  # smooth_GCV = FALSE: `smooth_tuning` gives the value to use for each PC; the
  #   smoothing matrix is computed for that PC inside the loop.
  smooth_grid <- !(smooth_GCV == FALSE)
  zero_smoothing <- function() {
    if (!is.null(mvmfd_obj)) expand.grid(lapply(rep(0,mvmfd_obj$nvar), function(x) x[1])) else expand.grid(lapply(rep(0,nf_obj$nvar), function(x) x[1]))
  }

  S_smooth <- S_2_inverse <- list()
  if (smooth_grid) {
    smooth_tuning_temp <- if (is.null(smooth_tuning)) zero_smoothing() else expand.grid(smooth_tuning)
    message("Preprocessing...")
    pb <- txtProgressBar(min = 0, max = dim(smooth_tuning_temp)[1], style = 3, width = 50, char = "=")
    #####Do the following computation in advance to save computational cost#####
    if (!is.null(mvmfd_obj)){
      for (smooth_index in 1:dim(smooth_tuning_temp)[1]) {
        setTxtProgressBar(pb, smooth_index)
        D <- I_alpha(mvmfd_obj, smooth_tuning_temp[smooth_index, ]) %*% smooth_penalty
        S_2 <- solve(G + D)
        S_2_inverse[[smooth_index]] = solve(S_2)
        S_smooth[[smooth_index]] <- G_half %*% (S_2) %*% G_half
      }
    } else {
      setTxtProgressBar(pb, 1)
      S_2_inverse = NULL
      S_smooth <- NULL
    }
    close(pb)
  }

  fv_total = c()
  nfv_total = c()
  # Scores are stored per PC in lists; with fixed smoothing and fixed sparsity
  # (smooth_GCV = FALSE, sparse_CV = FALSE) there is one score per PC, so they
  # are kept as a vector, and GCV scores likewise whenever smoothing is fixed.
  concat_scores <- !smooth_grid && sparse_CV == FALSE
  GCV_score <- if (smooth_grid) list() else c()
  if (concat_scores) {
    CV_score_u <- CV_score_nfd <- CV_score_fd <- cv_se_u <- cv_se_nfd <- cv_se_fd <- c()
    eBIC_score_u <- eBIC_score_nfd <- eBIC_score_fd <- c()
  } else {
    CV_score_u <- CV_score_nfd <- CV_score_fd <- cv_se_u <- cv_se_nfd <- cv_se_fd <- list()
    eBIC_score_u <- eBIC_score_nfd <- eBIC_score_fd <- list()
  }

  for (i in 1:n) {
    message(sprintf("Computing the %s PC...", ordinal(i)))
    if (!smooth_grid) {
      smooth_tuning_temp <- if (is.null(smooth_tuning)) zero_smoothing() else expand.grid(lapply(smooth_tuning, function(x) x[i]))
    }
    if (i == 1) {
      C_temp = C
      nf_data_temp <- nf_data
    } else {
      if (!is.null(mvmfd_obj)) {
        b_original = t(C_temp)%*%u
        C_temp = C_temp - u%*%t(b_original)
      } else {
        b_original = NULL
        C_temp = NULL
      }
      if (!is.null(nf_data)) nf_data_temp <- nf_data_temp - u %*% t(t(nf_data_temp)%*%u)
    }
    if (!smooth_grid) {
      if (!is.null(mvmfd_obj)){
        D <- I_alpha(mvmfd_obj, smooth_tuning_temp) %*% smooth_penalty
        S_2 <- solve(G + D)
        S_2_inverse[[1]] = solve(S_2)
        S_smooth[[1]] <- G_half %*% (S_2) %*% G_half
      } else {
        S_2_inverse[[1]] = NULL
        S_smooth[[1]] <- NULL
      }
    }

    if (!is.null(sparse_tuning_u)) {
      sparse_tuning_temp_u <- if (sparse_CV == FALSE) sparse_tuning_u[i] else sparse_tuning_u
    }
    if (!is.null(sparse_tuning_nfd)) {
      if (sparse_CV == FALSE) {
        sparse_tuning_temp_nfd <- lapply(sparse_tuning_nfd,function(x) x[i])
      } else {
        sparse_tuning_temp_nfd <- sparse_tuning_nfd
      }
    }
    if (!is.null(sparse_tuning_fd)) {
      if (sparse_CV == FALSE) {
        sparse_tuning_temp_fd <- lapply(sparse_tuning_fd,function(x) x[i])
      } else {
        sparse_tuning_temp_fd <- sparse_tuning_fd
      }
    }

    cv_result = cv_gcv_sequential_hybrid(
      # tune on the deflated data of this component, as used by the fit below
      fdata = C_temp,
      nfdata = nf_data_temp,
      hd_obj = hd_obj,
      smooth_tuning = if (is.null(smooth_tuning)) smooth_tuning else smooth_tuning_temp,
      sparse_tuning_u = if (is.null(sparse_tuning_u)) sparse_tuning_u else sparse_tuning_temp_u,
      sparse_tuning_nfd = if (is.null(sparse_tuning_nfd)) sparse_tuning_nfd else sparse_tuning_temp_nfd,
      sparse_tuning_fd = if (is.null(sparse_tuning_fd)) sparse_tuning_fd else sparse_tuning_temp_fd,
      G = G,
      G_half = G_half,
      G_half_inverse = G_half_inverse,
      S_smooth = S_smooth,
      S_2_inverse = S_2_inverse,
      cl = cl,
      settings = settings
    )
    sparse_result_u = cv_result$sparse_tuning_selection_u
    sparse_result_nfd = cv_result$sparse_tuning_selection_nfd
    sparse_result_fd = cv_result$sparse_tuning_selection_fd
    smooth_result_index = cv_result$index_selection

    # Store CV or eBIC scores based on selection method
    if (concat_scores) {
      if (sparse_selection_method == "cv") {
        CV_score_u <- c(CV_score_u, cv_result$cv_scores_u)
        CV_score_nfd <- c(CV_score_nfd, cv_result$cv_scores_nfd)
        CV_score_fd <- c(CV_score_fd, cv_result$cv_scores_fd)
        cv_se_u <- c(cv_se_u, cv_result$cv_se_u)
        cv_se_fd <- c(cv_se_fd, cv_result$cv_se_fd)
        cv_se_nfd <- c(cv_se_nfd, cv_result$cv_se_nfd)
      } else {
        eBIC_score_u <- c(eBIC_score_u, list(cv_result$eBIC_scores_u))
        eBIC_score_nfd <- c(eBIC_score_nfd, list(cv_result$eBIC_scores_nfd))
        eBIC_score_fd <- c(eBIC_score_fd, list(cv_result$eBIC_scores_fd))
      }
    } else {
      if (sparse_selection_method == "cv") {
        CV_score_u[[i]] <- cv_result$cv_scores_u
        CV_score_nfd[[i]] <- cv_result$cv_scores_nfd
        CV_score_fd[[i]] <- cv_result$cv_scores_fd
        cv_se_u[[i]] <- cv_result$cv_se_u
        cv_se_fd[[i]] <- cv_result$cv_se_fd
        cv_se_nfd[[i]] <- cv_result$cv_se_nfd
      } else {
        eBIC_score_u[[i]] <- cv_result$eBIC_scores_u
        eBIC_score_nfd[[i]] <- cv_result$eBIC_scores_nfd
        eBIC_score_fd[[i]] <- cv_result$eBIC_scores_fd
      }
    }
    if (smooth_grid) {
      GCV_score[[i]] = cv_result$gcv_scores
    } else {
      GCV_score = c(GCV_score, cv_result$gcv_scores)
    }

    # Smoothing matrix for the fit: the one selected by GCV, or the PC's own.
    smooth_k <- if (smooth_grid) smooth_result_index else 1
    CG_temp <- if(!is.null(mvmfd_obj)) C_temp%*%G_half else NULL
    test_result = init_sequential_hybrid(fdata = CG_temp,
                                         nfdata = nf_data_temp,
                                         sparse_tuning_result_u = sparse_result_u,
                                         sparse_tuning_result_nfd = sparse_result_nfd,
                                         sparse_tuning_result_fd = sparse_result_fd,
                                         sparse_tuning_type_u = sparse_tuning_type_u,
                                         sparse_tuning_type_nfd = sparse_tuning_type_nfd,
                                         sparse_tuning_type_fd = sparse_tuning_type_fd,
                                         S_smooth = S_smooth[[smooth_k]],
                                         S_2_inverse = S_2_inverse[[smooth_k]],
                                         G_half_inverse = G_half_inverse,
                                         G_half = G_half,
                                         pen_nfd = pen_nfd,
                                         pen_fd = pen_fd,
                                         pen_u = pen_u,
                                         tol = tol, max_iter = max_iter)
    u = test_result[[3]]
    fv = test_result[[1]]
    nfv = test_result[[2]]
    smooth_tuning_result[[i]] = smooth_tuning_temp[smooth_result_index, ]
    sparse_tuning_result_u[[i]] = sparse_result_u
    sparse_tuning_result_nfd[[i]] = sparse_result_nfd
    sparse_tuning_result_fd[[i]] = sparse_result_fd
    # Name the selections by variable.
    names( sparse_tuning_result_u[[i]]) <- paste0("Var",seq_len(length( sparse_tuning_result_u[[i]])))
    if (!is.null(sparse_result_nfd)){
      names(sparse_tuning_result_nfd[[i]]) <- paste0("Var",seq_len(length(sparse_tuning_result_nfd[[i]])))
    }
    if (!is.null(sparse_result_fd)){
      names(sparse_tuning_result_fd[[i]]) <- paste0("Var",seq_len(length(sparse_tuning_result_fd[[i]])))
    }
    temp_count <- 0
    variance_result <- handle_variance_update_hybrid(i, n, C, nf_data,G, fv_total, nfv_total,hd_obj, all_equal_check, c(sparse_tuning_u,sparse_tuning_nfd,sparse_tuning_fd), pc, lsv, fv, nfv,u, G_half, test_result, temp_count, p)
    pc <- variance_result$pc
    lsv <- variance_result$lsv
    fv_total <- variance_result$fv_total
    nfv_total <- variance_result$nfv_total
    variance[i] <- variance_result$variance
  }
  # Name the selections by component.
  names(sparse_tuning_result_u) <- paste0("PC",seq_len(length(sparse_tuning_result_u)))
  if (!length(sparse_tuning_result_nfd) == 0){
    names(sparse_tuning_result_nfd) <- paste0("PC",seq_len(length(sparse_tuning_result_nfd)))
  }
  if (!length(sparse_tuning_result_fd) == 0){
    names(sparse_tuning_result_fd) <- paste0("PC",seq_len(length(sparse_tuning_result_fd)))
  }
  if (is.null(smooth_tuning)) {
    GCV_score = NULL
  }
  if (isTRUE(sparse_CV) && sparse_selection_method == "cv"){
    cv_tables <- function(block, grid, scores, ses, enabled) {
      if (!is.null(grid) && isTRUE(enabled)) format_cv_scores(block, grid, scores, ses) else NULL
    }
    combined_fd <- cv_tables("fd", sparse_tuning_fd, CV_score_fd, cv_se_fd, pen_fd)
    combined_nfd <- cv_tables("nfd", sparse_tuning_nfd, CV_score_nfd, cv_se_nfd, pen_nfd)
    combined_u <- cv_tables("u", sparse_tuning_u, CV_score_u, cv_se_u, pen_u)
  } else {
    combined_fd <- combined_nfd <- combined_u <- NULL
  }
  
  # Format eBIC scores for output
  if (isTRUE(sparse_CV) && sparse_selection_method == "ebic") {
    combined_eBIC_u <- eBIC_score_u
    combined_eBIC_nfd <- eBIC_score_nfd
    combined_eBIC_fd <- eBIC_score_fd
    
    if (length(combined_eBIC_u) > 0) {
      names(combined_eBIC_u) <- paste0("PC", seq_len(length(combined_eBIC_u)))
    }
    if (length(combined_eBIC_nfd) > 0) {
      names(combined_eBIC_nfd) <- paste0("PC", seq_len(length(combined_eBIC_nfd)))
    }
    if (length(combined_eBIC_fd) > 0) {
      names(combined_eBIC_fd) <- paste0("PC", seq_len(length(combined_eBIC_fd)))
    }
  } else {
    combined_eBIC_u <- combined_eBIC_nfd <- combined_eBIC_fd <- NULL
  }
  
   
  
  return(list(
    pc_fd = pc, 
    pc_nfd = nfv_total,
    lsv = lsv, 
    variance = variance, 
    smooth_tuning_result = smooth_tuning_result, 
    sparse_tuning_result_u = sparse_tuning_result_u,
    sparse_tuning_result_nfd = sparse_tuning_result_nfd, 
    sparse_tuning_result_fd = sparse_tuning_result_fd,
    CV_score_u = combined_u,
    CV_score_nfd = combined_nfd, 
    CV_score_fd = combined_fd, 
    eBIC_score_u = combined_eBIC_u,
    eBIC_score_nfd = combined_eBIC_nfd,
    eBIC_score_fd = combined_eBIC_fd,
    GCV_score = GCV_score,
    selection_method = sparse_selection_method
    ))
} 



# joint smooth and sparse power algorithm
joint_power_hybrid <- function(hd_obj, n, smooth_tuning, settings) {
  smooth_tuning_type <- settings$smoothing_type
  centerfns <- settings$centerfns
  tol <- settings$tol
  max_iter <- settings$max_iter

  #######centralize########
  if (centerfns) hd_obj <- center_hd(hd_obj)
  mvmfd_obj <- hd_obj$mf
  nf_obj <- hd_obj$nf
  nf_data <- if (!is.null(nf_obj)) do.call("cbind",nf_obj$data) else NULL
  if (!is.null(mvmfd_obj)){
    p <- mvmfd_obj$nvar
    smooth_penalty <- pen_fun(mvmfd_obj, type = smooth_tuning_type)
    C <- do.call("rbind", lapply(mvmfd_obj$coefs, coef_matrix))  
    ########some initial setting#######
    C <- t(C)
    lsv <- c()
    pc <- list()
    variance <- vector()
    smooth_tuning_result  <- list()
    sparse_tuning_result  <- list()
    G <- as.matrix(mvmfd_obj$basis$gram)
    G_half <- sqrtm(G)
    G_half_inverse = solve(G_half)
    ###################################
  } else {
    p <- NULL
    smooth_penalty <- NULL
    C <- NULL  
    ########some initial setting#######
    lsv <- c()
    pc <- list()
    variance <- vector()
    smooth_tuning_result  <- list()
    sparse_tuning_result  <- list()
    G <- NULL
    G_half <- NULL
    G_half_inverse = NULL
    ###################################
  }
  
  
  #####smoothing parameter#######
  if (is.null(smooth_tuning)) {
    smooth_tuning_temp = if (!is.null(mvmfd_obj)) expand.grid(lapply(rep(0,mvmfd_obj$nvar), function(x) x[1])) else NULL
  } else{
    smooth_tuning_temp <- expand.grid(smooth_tuning)
  }
  ####################################
  
  ####joint power####
  # I_a <- D <- S_2 <- S_smooth <- S_2_inverse <- list()
  S_smooth <- S_2_inverse <- list()
  n_iter1 <- dim(smooth_tuning_temp)[1]
  if (!is.null(mvmfd_obj)){
    message("Preprocessing...")
    
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = n_iter1, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
  }
  
  #####Do the following computation in advance to save computational cost#####
  if (!is.null(mvmfd_obj)){
    for (smooth_index in 1:dim(smooth_tuning_temp)[1]) {
      setTxtProgressBar(pb, smooth_index)
      D <- I_alpha(mvmfd_obj, smooth_tuning_temp[smooth_index, ]) %*% smooth_penalty
      S_2 <- solve(G + D)
      S_2_inverse[[smooth_index]] = solve(S_2)
      S_smooth[[smooth_index]] <- G_half %*% (S_2) %*% G_half
    }
    close(pb)
  } else {
      D <- NULL
      S_2 <- NULL
      S_2_inverse = NULL
      S_smooth <- NULL
  }
  
  if (!is.null(mvmfd_obj)){
    message("Computing PCs...")
  }
  

  cv_result = gcv_joint_hybrid(
    fdata = C, 
    nfdata = nf_data,
    hd_obj = hd_obj, 
    smooth_tuning = if (is.null(smooth_tuning)) smooth_tuning else smooth_tuning_temp, 
    G = G, 
    G_half = G_half, 
    G_half_inverse = G_half_inverse, 
    S_smooth = S_smooth, 
    S_2_inverse = S_2_inverse, n = n, settings = settings
  )
  smooth_result_index = cv_result[[2]]
  GCV_score = cv_result[[3]]
  CG_temp <- if (!is.null(mvmfd_obj)) C%*%G_half else NULL
  test_result = init_joint_hybrid(CG_temp, nf_data,S_smooth[[smooth_result_index]], S_2_inverse[[smooth_result_index]], G_half_inverse, G_half, n = n,tol = tol, max_iter = max_iter)
  u = test_result[[3]]
  fv = test_result[[1]]
  nfv <- test_result[[2]]
  smooth_tuning_result = smooth_tuning_temp[smooth_result_index, ]
  temp_count <- 0
  
  temp_count <- 0
  if (!is.null(mvmfd_obj)){
    for (j in 1:p) {
      index_start <- (temp_count + 1)
      index_end <- (temp_count + prod(mvmfd_obj$basis$nbasis[[j]]))
      pc[[j]] <- fv[index_start:index_end, ]
      temp_count <- index_end
    }
  } else {
    pc <- NULL
    }
  
  
  lsv = cbind(lsv, u)
  fv_total = fv
  nfv_total = nfv
  for (k in 1:n) {
    CGfv <- if (!is.null(mvmfd_obj)) C %*% G %*% fv[,k] else NULL
    Xnfv <- if (!is.null(nf_data)) nf_data %*% nfv[,k] else NULL
    if (!is.null(mvmfd_obj) && !is.null(nf_data) ) {
      variance[k] <- (t(Xnfv)%*%Xnfv+2*t(Xnfv)%*% CGfv+t(CGfv)%*% CGfv)/(mvmfd_obj$nobs - 1)
    } else if (!is.null(mvmfd_obj) && is.null(nf_data) ){
      variance[k] <- (t(CGfv)%*% CGfv)/(mvmfd_obj$nobs - 1)
    } else if (is.null(mvmfd_obj) && !is.null(nf_data) ){
      variance[k] <- (t(Xnfv)%*%Xnfv)/(nf_obj$nobs - 1)
    } 
  }
  
  return(list(pc_fd = pc, pc_nfd = nfv_total,lsv = lsv, variance = variance, smooth_tuning_result = smooth_tuning_result, GCV_score = GCV_score))
} 

# Coefficients of one functional variable as a (basis functions x observations)
# matrix; 2-D variables store them as an array (nbasis1 x nbasis2 x nobs).
coef_matrix <- function(coefs) {
  if (is.matrix(coefs)) coefs else apply(coefs, length(dim(coefs)), as.vector)
}
