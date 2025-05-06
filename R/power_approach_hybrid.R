#' @importFrom expm sqrtm
#' @importFrom utils  txtProgressBar setTxtProgressBar
#' @importFrom parallel clusterExport parLapplyLB makeCluster clusterEvalQ stopCluster

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
  
  # Otherwise, tuning_parameter is a single number.
  if (length(y) < tuning_parameter) {
    stop("tuning parameter should be less than number of instances.")
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


# sequential power algorithm
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
                                   penalize_u = FALSE,
                                   penalize_nfd = FALSE,
                                   penalize_fd = FALSE) {
  
  fv_old <- if (!is.null(fdata)) svd(fdata)$v[, 1] else NULL
  nfv_old <- if (!is.null(nfdata)) svd(nfdata)$v[, 1] else NULL
  errors <- 10^60
  
  while (errors > 0.00001) {
    
    u_old <- 0
    if (!is.null(fdata)) {
      u_old <- u_old + fdata %*% fv_old
    }
    if (!is.null(nfdata)) {
      u_old <- u_old + nfdata %*% nfv_old
    }

    if (penalize_u==TRUE) {
       u_old <- sparse_pen_fun(y = u_old, tuning_parameter = sparse_tuning_result_u, sparse_tuning_type_u)
    } 
    u_old <- u_old / norm(u_old, type = "2")
    if (cv_flag == TRUE && penalize_u == TRUE) { # incorporate power algorithm in CV sparse tuning selection for u
      

      fv_new <- if (!is.null(fdata)) t(fdata) %*% u_old else NULL

      if (!is.null(nfdata)) {
        nfv_new <- t(nfdata) %*% u_old
      }

      if (!is.null(fdata) && !is.null(nfdata)) {
        coef_norm <- sqrt(norm(fv_new, "2")^2 + norm(nfv_new, "2")^2)
        fv_new <- fv_new / coef_norm
        nfv_new <- nfv_new / coef_norm
      } else if (!is.null(fdata) && is.null(nfdata)) {
        coef_norm <- norm(fv_new, "2")
        fv_new <- fv_new / coef_norm
      } else if (is.null(fdata) && !is.null(nfdata)) {
        coef_norm <- norm(nfv_new, "2")
        nfv_new <- nfv_new / coef_norm
      }

      u_new <- 0
      if (!is.null(fdata)) {
        u_new <- u_new + fdata %*% fv_new
      }
      if (!is.null(nfdata)) {
        u_new <- u_new + nfdata %*% nfv_new
      }

      if (penalize_u==TRUE) {
        u_new <-  sparse_pen_fun(y = u_new, tuning_parameter = sparse_tuning_result_u, sparse_tuning_type_u)
      } 
      u_new <- u_new / norm(u_new, type = "2")
      errors <- sum((u_new - u_old)^2)
      if (!is.null(fdata)) fv_old <- fv_new
      if (!is.null(nfdata)) nfv_old <- nfv_new
      u_old <- u_new
    } else if(cv_flag == TRUE && (penalize_nfd == TRUE || penalize_fd == TRUE)){ # incorporate power algorithm in CV sparse tuning selection for nfd or fd
      if (!is.null(fdata)){
        fv_new <- t(fdata) %*% u_old # v_tilde
        
        # v_tilde
        fv_new <- if(penalize_fd==TRUE) G_half%*%sparse_pen_fun(y = G_half_inverse%*%fv_new, tuning_parameter = sparse_tuning_result_fd, type = sparse_tuning_type_fd) else fv_new 
        #fv_new_back <- G_half_inverse %*% fv_new # real v 
        fv_weight <- norm(fv_new,"2")^2 #fv_weight <- norm(fv_new_back,"2") 
      } else {
        fv_weight <- 0
        fv_new <- NULL
      }
      
      
      if (!is.null(nfdata)) {
        nfv_new <- t(nfdata) %*% u_old
        
        nfv_new <- if(penalize_nfd==TRUE) sparse_pen_fun(y = nfv_new, tuning_parameter = sparse_tuning_result_nfd, type = sparse_tuning_type_nfd) else nfv_new
        nfv_weight <- norm(nfv_new,"2")^2
      } else {
        nfv_weight <- 0
        nfv_new <- NULL
      }
      coef_norm <- sqrt(fv_weight + nfv_weight)
      if (!is.null(fdata) && !is.null(nfdata)) {
        
        #fv_new <- G_half %*%(fv_new_back / coef_norm)
        fv_new <- (fv_new / coef_norm)
        nfv_new <- nfv_new / coef_norm
      } else if (!is.null(fdata) && is.null(nfdata)) {
        
        #fv_new <- G_half %*%(fv_new_back / coef_norm)
        fv_new <- (fv_new / coef_norm)
      } else if (is.null(fdata) && !is.null(nfdata)) {
        
        nfv_new <- nfv_new / coef_norm
      }
      
      u_new <- 0
      if (!is.null(fdata)) {
        u_new <- u_new + fdata %*% fv_new
      }
      if (!is.null(nfdata)) {
        u_new <- u_new + nfdata %*% nfv_new
      }
      
      u_new <- u_new / norm(u_new, type = "2")
      errors <- sum((u_new - u_old)^2)
      if (!is.null(fdata)) fv_old <- fv_new
      if (!is.null(nfdata)) nfv_old <- nfv_new
      u_old <- u_new
      
      
    } else if(cv_flag == FALSE) { # incorporate power algorithm smoothing
      
      
      if (!is.null(fdata)) {
        fv_new <- S_smooth %*% t(fdata) %*% u_old
        fv_new <- if (penalize_fd==TRUE) G_half%*%sparse_pen_fun(y = G_half_inverse%*%fv_new, tuning_parameter = sparse_tuning_result_fd, type = sparse_tuning_type_fd) else fv_new
        fv_weight <- norm(fv_new,"2")^2
        } else {
        fv_new <- NULL
        fv_weight <- 0
      }
      
      if (!is.null(nfdata)) {
        nfv_new <- t(nfdata) %*% u_old
        nfv_new <- if (penalize_nfd==TRUE) sparse_pen_fun(y = nfv_new, tuning_parameter = sparse_tuning_result_nfd, type = sparse_tuning_type_nfd) else nfv_new
        nfv_weight <-  norm(nfv_new, "2")^2
        } else {
        nfv_new <- NULL
        nfv_weight <- 0
        }
      

      # normalizing fv and nfv
      ##########################################
     
      coef_norm <- sqrt(fv_weight + nfv_weight)
      if (!is.null(fdata)) {
        #fv_new_back <- fv_new_back / coef_norm
        #fv_new <- G_half %*% fv_new_back
        fv_new <- fv_new/coef_norm
      }
      if (!is.null(nfdata)) {
        nfv_new <- nfv_new / coef_norm
      }
      ##########################################
      u_new <- 0
      if (!is.null(fdata)) {
        u_new <- u_new + fdata %*% fv_new
      }
      if (!is.null(nfdata)) {
        u_new <- u_new + nfdata %*% nfv_new
      }

      if (penalize_u ==TRUE) {
        u_new <-  sparse_pen_fun(y = u_new, tuning_parameter = sparse_tuning_result_u, sparse_tuning_type_u)
      } 
      u_new <- u_new / norm(u_new, type = "2")

      errors <- sum((u_new - u_old)^2)
      if (!is.null(fdata)) fv_old <- fv_new
      if (!is.null(nfdata)) nfv_old <- nfv_new
    }
  }
  #
  if (cv_flag == TRUE && penalize_u==TRUE) {
    return(u_new)
  } else if (cv_flag == TRUE && (penalize_nfd==TRUE || penalize_fd==TRUE)){
    coef_norm <- 0
    if (is.null(nfdata)){
      nfv_new <- NULL
      nf_weight <- 0
    } else {
      coef_norm <- sqrt(coef_norm + norm(nfv_new,"2")^2)
    }
    if (is.null(fdata)){
      fv_new <- NULL 
      f_weight <- 0
    }  else {
      fv_new <- G_half_inverse%*%fv_new # output is read fv not tilde
      coef_norm <- sqrt(coef_norm^2 + norm(fv_new,"2")^2)
      fv_new <- fv_new/coef_norm
    }
    if (!is.null(nfdata)){
      nfv_new <- nfv_new/coef_norm
    }
    return(list(fv = fv_new,nfv = nfv_new))
  } else if (cv_flag == FALSE){
    coef_norm <- 0
    
    if (is.null(nfdata)){
      nfv_new <- NULL
    } else {
      coef_norm <- sqrt(coef_norm^2 + norm(nfv_new,"2")^2)
    }
    if (!is.null(fdata)){
      fv_new <- G_half_inverse %*% fv_new
      coef_norm <- sqrt(coef_norm^2 + as.numeric(sqrt(t(fv_new) %*% S_2_inverse %*% fv_new))^2)
     
      fv_new <- fv_new / coef_norm
    }  else {
      fv_new <- NULL
    } 
    if (!is.null(nfdata)){
      nfv_new <- nfv_new/coef_norm
    }
    return(list(fv = fv_new, nfv = nfv_new, u = u_new))
  }
}

#joint power for smoothing
init_joint_hybrid = function(fdata, nfdata, S_smooth = NULL, S_2_inverse = NULL, G_half_inverse = NULL, G_half = NULL, n = n){
  fv_old <- if (!is.null(fdata)) svd(fdata)$v[, 1:n] else NULL
  nfv_old <- if (!is.null(nfdata)) svd(nfdata)$v[, 1:n] else NULL
  
  errors = 10^60
  while (errors > 10^-10) {
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
    K_fold_u,
    K_fold_nfd,
    K_fold_fd,
    sparse_tuning_single_u,
    sparse_tuning_single_nfd,
    sparse_tuning_single_fd,
    sparse_tuning_type_u,
    sparse_tuning_type_nfd,
    sparse_tuning_type_fd,
    shuffled_row_u,
    shuffled_row_nfd,
    shuffled_row_fd,
    group_size_u,
    group_size_nfd,
    group_size_fd,
    penalize_nfd = FALSE,
    penalize_fd  = FALSE,
    penalize_u   = FALSE
) {

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
  if (penalize_u) {
    for (k in seq_len(K_fold_u)) {
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
        penalize_nfd          = FALSE,
        penalize_u            = penalize_u
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
  }

  #----------------------------------------------
  # 2) Cross-validation over 'nfd' terms
  #----------------------------------------------

  err_nfd <- NULL

  if (penalize_nfd) {

    for (k in seq_len(K_fold_nfd)) {
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
        penalize_nfd          = penalize_nfd,
        penalize_fd           = FALSE,
        penalize_u            = FALSE
      )


      nfv_test <- v_test$nfv

      # Combine contributions into u_test
      u_test <- 0
      
      u_test <- u_test + nfdata_test %*% nfv_test


      nfv_error <- sum((t(nfdata_test) - nfv_test %*% t(u_test))^2)


      error_score_sparse_nfd <- error_score_sparse_nfd + nfv_error # (fv_error + nfv_error)

    }

    # Normalization
    normalization_factor <- nrow(nfdata)
    err_nfd <- error_score_sparse_nfd / normalization_factor

  }

  #----------------------------------------------
  # 3) Cross-validation over 'fd' terms
  #----------------------------------------------

  err_fd <- NULL

  if (penalize_fd) {

    for (k in seq_len(K_fold_fd)) {
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
        penalize_nfd          = FALSE,
        penalize_fd           = penalize_fd,
        penalize_u            = FALSE
      )

      fv_test_smooth_back <- v_test$fv

      # Combine contributions into u_test
      u_test <- 0
      u_test <- u_test + t(fdata_test) %*% fv_test_smooth_back

      # Compute errors for this fold

      fdata_test_smooth_back <- t(data_double_tilde)[rows_to_remove_fd, , drop = FALSE]
      fv_error <- sum((t(fdata_test_smooth_back) - fv_test_smooth_back %*% t(u_test))^2)

      error_score_sparse_fd <- error_score_sparse_fd + fv_error 

    }

    # Normalization
    normalization_factor <- ncol(fdata)
    err_fd <- error_score_sparse_fd / normalization_factor

  }

  #--------------------------------
  # Return the same output structure
  #--------------------------------
  list(err_u = err_u, err_nfd = err_nfd,err_fd = err_fd)
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
                                        sparse_tuning_type_u = NULL,
                                        sparse_tuning_type_nfd = NULL,
                                        sparse_tuning_type_fd = NULL,
                                        smooth_tuning, 
                                        CV_score_smooth, 
                                        power_type = "sequential", 
                                        n = NULL, 
                                        pb, 
                                        count, 
                                        penalize_nfd = F,
                                        penalize_fd = F,
                                        penalize_u = F) {
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
                                              penalize_nfd = penalize_nfd,
                                              penalize_fd = penalize_fd,
                                              penalize_u = penalize_u)
        } else {
          test_temp <- init_joint_hybrid(fdata %*% G_half, nfdata ,S_smooth[[smooth_index]], S_2_inverse[[smooth_index]], G_half_inverse, G_half, n = n)
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
    sparse_tuning_type_u,
    sparse_tuning_type_nfd,
    sparse_tuning_type_fd,
    K_fold_u,
    K_fold_nfd,
    K_fold_fd,
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
    penalize_nfd = FALSE,
    penalize_fd  = FALSE,
    penalize_u   = FALSE,
    cl   # <<— cluster passed in from outside
) {
  # Collect the names of the objects to export to the workers:
  data_vars <- c(
    "fdata", "nfdata",
    "G_half", "G_half_inverse",
    "K_fold_u", "K_fold_nfd", "K_fold_fd",
    "sparse_tuning_type_u", "sparse_tuning_type_nfd", "sparse_tuning_type_fd",
    "shuffled_row_u", "shuffled_row_nfd", "shuffled_row_fd",
    "group_size_u", "group_size_nfd", "group_size_fd"
  )
  
  # Quick exit if nothing to tune:
  if (is.null(sparse_tuning_u) &&
      is.null(sparse_tuning_nfd) &&
      is.null(sparse_tuning_fd)) {
    setTxtProgressBar(pb, 4)
    return(list(
      sparse_tuning_selection_u   = 0,
      sparse_tuning_selection_nfd = 0,
      sparse_tuning_selection_fd  = NULL,
      cv_scores_u                 = NULL,
      cv_scores_nfd               = NULL,
      cv_scores_fd                = NULL,
      CV_score_sparse_u           = CV_score_sparse_u,
      CV_score_sparse_nfd         = CV_score_sparse_nfd,
      CV_score_sparse_fd          = CV_score_sparse_fd
    ))
  }
  
  count <- 0
  
  # --- 1) tune 'u' ---
  if (!is.null(sparse_tuning_u)) {
    parallel::clusterExport(cl, varlist = c(data_vars, "sparse_tuning_u"), envir = environment())
    cv_scores_u_list <- parallel::parLapplyLB(cl, sparse_tuning_u, function(candidate_u) {
      res <- cv_local_hybrid(
        fdata   = fdata,
        nfdata  = nfdata,
        G_half  = G_half,
        G_half_inverse = G_half_inverse,
        K_fold_u       = K_fold_u,
        K_fold_nfd     = K_fold_nfd,
        K_fold_fd      = K_fold_fd,
        sparse_tuning_single_u   = candidate_u,
        sparse_tuning_single_nfd = NULL,
        sparse_tuning_single_fd  = NULL,
        sparse_tuning_type_u     = sparse_tuning_type_u,
        sparse_tuning_type_nfd   = sparse_tuning_type_nfd,
        sparse_tuning_type_fd    = sparse_tuning_type_fd,
        shuffled_row_u   = shuffled_row_u,
        shuffled_row_nfd = shuffled_row_nfd,
        shuffled_row_fd  = shuffled_row_fd,
        group_size_u   = group_size_u,
        group_size_nfd = group_size_nfd,
        group_size_fd  = group_size_fd,
        penalize_nfd = FALSE,
        penalize_fd  = FALSE,
        penalize_u   = TRUE
      )
      res$err_u
    })
    cv_scores_u    <- unlist(cv_scores_u_list)
    best_index_u   <- which.min(cv_scores_u)
    best_u         <- sparse_tuning_u[best_index_u]
    CV_score_sparse_u <- cv_scores_u[best_index_u]
    count <- count + length(sparse_tuning_u)
    setTxtProgressBar(pb, count)
  } else {
    best_u      <- 0
    cv_scores_u <- NULL
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
  
  # --- 2) tune 'nfd' ---
  if (!is.null(sparse_tuning_nfd)) {
    parallel::clusterExport(cl, varlist = c(data_vars, "sparse_tuning_nfd", "best_u"), envir = environment())
    cv_scores_nfd_list <- parallel::parLapplyLB(cl, sparse_tuning_nfd, function(candidate_nfd) {
      res <- cv_local_hybrid(
        fdata   = fdata,
        nfdata  = nfdata,
        G_half  = G_half,
        G_half_inverse = G_half_inverse,
        K_fold_u       = K_fold_u,
        K_fold_nfd     = K_fold_nfd,
        K_fold_fd      = K_fold_fd,
        sparse_tuning_single_u   = best_u,
        sparse_tuning_single_nfd = candidate_nfd,
        sparse_tuning_single_fd  = NULL,
        sparse_tuning_type_u     = sparse_tuning_type_u,
        sparse_tuning_type_nfd   = sparse_tuning_type_nfd,
        sparse_tuning_type_fd    = sparse_tuning_type_fd,
        shuffled_row_u   = shuffled_row_u,
        shuffled_row_nfd = shuffled_row_nfd,
        shuffled_row_fd  = shuffled_row_fd,
        group_size_u   = group_size_u,
        group_size_nfd = group_size_nfd,
        group_size_fd  = group_size_fd,
        penalize_nfd = TRUE,
        penalize_fd  = FALSE,
        penalize_u   = penalize_u
      )
      res$err_nfd
    })
    cv_scores_nfd      <- unlist(cv_scores_nfd_list)
    best_index_nfd     <- which.min(cv_scores_nfd)
    best_nfd           <- sparse_tuning_nfd[best_index_nfd]
    CV_score_sparse_nfd <- cv_scores_nfd[best_index_nfd]
    count <- count + length(sparse_tuning_nfd)
    setTxtProgressBar(pb, count)
  } else {
    best_nfd       <- 0
    cv_scores_nfd  <- NULL
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
  
  # --- 3) tune 'fd' ---
  if (!is.null(sparse_tuning_fd)) {
    parallel::clusterExport(cl, varlist = c(data_vars, "sparse_tuning_fd", "best_u", "best_nfd"), envir = environment())
    cv_scores_fd_list <- parallel::parLapplyLB(cl, seq_len(nrow(sparse_tuning_fd)), function(i) {
      candidate_fd <- sparse_tuning_fd[i, , drop = TRUE]
      res <- cv_local_hybrid(
        fdata   = fdata,
        nfdata  = nfdata,
        G_half  = G_half,
        G_half_inverse = G_half_inverse,
        K_fold_u       = K_fold_u,
        K_fold_nfd     = K_fold_nfd,
        K_fold_fd      = K_fold_fd,
        sparse_tuning_single_u   = best_u,
        sparse_tuning_single_nfd = best_nfd,
        sparse_tuning_single_fd  = candidate_fd,
        sparse_tuning_type_u     = sparse_tuning_type_u,
        sparse_tuning_type_nfd   = sparse_tuning_type_nfd,
        sparse_tuning_type_fd    = sparse_tuning_type_fd,
        shuffled_row_u   = shuffled_row_u,
        shuffled_row_nfd = shuffled_row_nfd,
        shuffled_row_fd  = shuffled_row_fd,
        group_size_u   = group_size_u,
        group_size_nfd = group_size_nfd,
        group_size_fd  = group_size_fd,
        penalize_nfd = penalize_nfd,
        penalize_fd  = TRUE,
        penalize_u   = penalize_u
      )
      res$err_fd
    })
    cv_scores_fd    <- unlist(cv_scores_fd_list)
    best_index_fd   <- which.min(cv_scores_fd)
    best_fd         <- sparse_tuning_fd[best_index_fd, , drop = TRUE]
    CV_score_sparse_fd <- cv_scores_fd[best_index_fd]
    count <- count + nrow(sparse_tuning_fd)
    setTxtProgressBar(pb, count)
  } else {
    best_fd      <- NULL
    cv_scores_fd <- NULL
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
  
  list(
    sparse_tuning_selection_u   = best_u,
    sparse_tuning_selection_nfd = best_nfd,
    sparse_tuning_selection_fd  = best_fd,
    cv_scores_u                 = cv_scores_u,
    cv_scores_nfd               = cv_scores_nfd,
    cv_scores_fd                = cv_scores_fd,
    CV_score_sparse_u           = CV_score_sparse_u,
    CV_score_sparse_nfd         = CV_score_sparse_nfd,
    CV_score_sparse_fd          = CV_score_sparse_fd
  )
}

# Function for cv_gcv_sequential  !!!!!!!!!!!!!!!!!!!!! not complete
cv_gcv_sequential_hybrid_old <- function(fdata, 
                                     nfdata, 
                                     hd_obj, 
                                     smooth_tuning, 
                                     sparse_tuning_u,
                                     sparse_tuning_nfd,
                                     sparse_tuning_fd,
                                     sparse_tuning_type_u, 
                                     sparse_tuning_type_nfd, 
                                     sparse_tuning_type_fd, 
                                     K_fold_u, 
                                     K_fold_nfd,
                                     K_fold_fd,
                                     G, 
                                     G_half, 
                                     G_half_inverse, 
                                     S_smooth, 
                                     S_2_inverse, 
                                     penalize_nfd = FALSE,
                                     penalize_fd = FALSE,
                                     penalize_u = FALSE,
                                     n_cores) {
  
  mvmfd_obj <- hd_obj$mf
  CV_score_sparse_u <- CV_score_sparse_nfd <- CV_score_sparse_fd <- CV_score_smooth <- Inf
  result <- c()
  count <- 0
  nc <- 0 
  ncf <- if (!is.null(fdata))  ncol(fdata) else NULL
  ncnf <- if (!is.null(nfdata))  ncol(nfdata) else NULL
  shuffled_row_f <- sample(ncf)
  len_f <- if (!is.null(fdata)) length(shuffled_row_f) else 0
  shuffled_row_nf <- sample(ncnf)
  len_nf <- if (!is.null(nfdata)) length(shuffled_row_nf) else 0
  group_size_f <- len_f / K_fold_u 
  group_size_nf <- len_nf / K_fold_u
  group_size_u <- list(f = group_size_f,nf = group_size_nf)
  shuffled_row_u <- list(f = shuffled_row_f, nf = shuffled_row_nf)
  
  nrf <- if (!is.null(fdata))  nrow(fdata) else nrow(nfdata)
  shuffled_row_nfd <- shuffled_row_fd <- sample(nrf)
  group_size_nfd <- group_size_fd <- length(shuffled_row_nfd)/K_fold_nfd
  group_size_fd <- length(shuffled_row_nfd)/K_fold_fd
  
  n_iter <- (if (is.null(smooth_tuning)) 1 else dim(smooth_tuning)[1]) +( 
    (if (is.null(sparse_tuning_u)) 1 else length(sparse_tuning_u)) + 
    (if (is.null(sparse_tuning_nfd)) 1 else length(sparse_tuning_nfd)) + 
    (if (is.null(sparse_tuning_fd)) 1 else length(sparse_tuning_fd)))
  
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  # Handle sparse tuning
  
  sparse_tuning_result <- handle_sparse_tuning_hybrid(fdata = fdata, 
                                                      nfdata = nfdata, 
                                                      G_half = G_half, 
                                                      G_half_inverse = G_half_inverse,
                                                      sparse_tuning_u = sparse_tuning_u, 
                                                      sparse_tuning_nfd = sparse_tuning_nfd,
                                                      sparse_tuning_fd = sparse_tuning_fd,
                                                      sparse_tuning_type_u = sparse_tuning_type_u, 
                                                      sparse_tuning_type_nfd = sparse_tuning_type_nfd,
                                                      sparse_tuning_type_fd = sparse_tuning_type_fd,
                                                      K_fold_u = K_fold_u,
                                                      K_fold_nfd = K_fold_nfd,
                                                      K_fold_fd = K_fold_fd,
                                                      shuffled_row_u = shuffled_row_u, 
                                                      shuffled_row_nfd = shuffled_row_nfd,
                                                      shuffled_row_fd = shuffled_row_fd,
                                                      group_size_u = group_size_u, 
                                                      group_size_nfd = group_size_nfd,
                                                      group_size_fd = group_size_fd,
                                                      CV_score_sparse_u = CV_score_sparse_u,
                                                      CV_score_sparse_nfd = CV_score_sparse_nfd,
                                                      CV_score_sparse_fd = CV_score_sparse_fd,
                                                      pb = pb, 
                                                      penalize_nfd = penalize_nfd,
                                                      penalize_fd = penalize_fd, 
                                                      penalize_u = penalize_u,
                                                      n_cores = n_cores)
  
  sparse_tuning_selection_u <- sparse_tuning_result$sparse_tuning_selection_u
  sparse_tuning_selection_nfd <- sparse_tuning_result$sparse_tuning_selection_nfd
  sparse_tuning_selection_fd <- sparse_tuning_result$sparse_tuning_selection_fd
  cv_scores_u <- sparse_tuning_result$cv_scores_u
  cv_scores_nfd <- sparse_tuning_result$cv_scores_nfd
  cv_scores_fd <- sparse_tuning_result$cv_scores_fd
  CV_score_sparse_u <- sparse_tuning_result$CV_score_sparse_u
  CV_score_sparse_nfd <- sparse_tuning_result$CV_score_sparse_nfd
  CV_score_sparse_fd <- sparse_tuning_result$CV_score_sparse_fd
  
  # Handle smooth tuning
  countt <- (if (is.null(sparse_tuning_u)) 1 else length(sparse_tuning_u)) + 
    (if (is.null(sparse_tuning_nfd)) 1 else length(sparse_tuning_nfd)) + 
    (if (is.null(sparse_tuning_fd)) 1 else length(sparse_tuning_fd))
    
  smooth_tuning_result <- handle_smooth_tuning_hybrid(fdata = fdata, 
                                                      nfdata = nfdata,
                                                      G_half = G_half,
                                                      G = G, 
                                                      S_smooth = S_smooth,
                                                      S_2_inverse =  S_2_inverse, 
                                                      G_half_inverse = G_half_inverse,
                                                      hd_obj =  hd_obj, 
                                                      sparse_tuning_selection_u = sparse_tuning_selection_u, 
                                                      sparse_tuning_selection_nfd = sparse_tuning_selection_nfd,
                                                      sparse_tuning_selection_fd = sparse_tuning_selection_fd,
                                                      sparse_tuning_type_u = sparse_tuning_type_u, 
                                                      sparse_tuning_type_nfd = sparse_tuning_type_nfd,
                                                      sparse_tuning_type_fd = sparse_tuning_type_fd,
                                                      smooth_tuning = smooth_tuning, 
                                                      CV_score_smooth = CV_score_smooth, 
                                                      power_type = "sequential", 
                                                      pb = pb, 
                                                      count = countt
                                                      )
  smooth_tuning_selection <- smooth_tuning_result$smooth_tuning_selection
  index_selection <- smooth_tuning_result$index_selection
  gcv_scores <- smooth_tuning_result$gcv_scores
  
  close(pb)
  result <- list(sparse_tuning_selection_u = sparse_tuning_selection_u, 
                 smooth_tuning_selection = smooth_tuning_selection, 
                 index_selection = index_selection, cv_scores_u = cv_scores_u,
                 cv_scores_nfd = cv_scores_nfd, cv_scores_fd = cv_scores_fd,
                 gcv_scores = gcv_scores, sparse_tuning_selection_nfd = sparse_tuning_selection_nfd,
                 sparse_tuning_selection_fd = sparse_tuning_selection_fd)
  
  return(result)
}

cv_gcv_sequential_hybrid <- function(
    fdata,
    nfdata,
    hd_obj,
    smooth_tuning,
    sparse_tuning_u,
    sparse_tuning_nfd,
    sparse_tuning_fd,
    sparse_tuning_type_u,
    sparse_tuning_type_nfd,
    sparse_tuning_type_fd,
    K_fold_u,
    K_fold_nfd,
    K_fold_fd,
    G,
    G_half,
    G_half_inverse,
    S_smooth,
    S_2_inverse,
    penalize_nfd = FALSE,
    penalize_fd  = FALSE,
    penalize_u   = FALSE,
    cl
) {
  # Prepare fold assignments
  ncf   <- if (!is.null(fdata))  ncol(fdata)  else NULL
  ncnf  <- if (!is.null(nfdata)) ncol(nfdata) else NULL
  shuffled_row_u   <- list(f = sample(ncf), nf = sample(ncnf))
  group_size_u     <- list(f = length(shuffled_row_u$f)/K_fold_u,
                           nf= length(shuffled_row_u$nf)/K_fold_u)
  nrf   <- if (!is.null(fdata))  nrow(fdata)  else nrow(nfdata)
  shuffled_row_nfd <- shuffled_row_fd <- sample(nrf)
  group_size_nfd   <- length(shuffled_row_nfd)/K_fold_nfd
  group_size_fd    <- length(shuffled_row_nfd)/K_fold_fd
  
  # Determine progress‐bar length
  n_iter <- (if (is.null(smooth_tuning)) 1 else nrow(smooth_tuning)) +
    (if (is.null(sparse_tuning_u)) 1 else length(sparse_tuning_u)) +
    (if (is.null(sparse_tuning_nfd)) 1 else length(sparse_tuning_nfd)) +
    (if (is.null(sparse_tuning_fd)) 1 else nrow(sparse_tuning_fd))
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  
  # 1) Sparse tuning (parallel)
  sparse_res <- handle_sparse_tuning_hybrid(
    fdata   = fdata,      nfdata = nfdata,
    G_half  = G_half,     G_half_inverse = G_half_inverse,
    sparse_tuning_u   = sparse_tuning_u,
    sparse_tuning_nfd = sparse_tuning_nfd,
    sparse_tuning_fd  = sparse_tuning_fd,
    sparse_tuning_type_u   = sparse_tuning_type_u,
    sparse_tuning_type_nfd = sparse_tuning_type_nfd,
    sparse_tuning_type_fd  = sparse_tuning_type_fd,
    K_fold_u   = K_fold_u,
    K_fold_nfd = K_fold_nfd,
    K_fold_fd  = K_fold_fd,
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
    penalize_nfd = penalize_nfd,
    penalize_fd  = penalize_fd,
    penalize_u   = penalize_u,
    cl           = cl
  )
  sparse_tuning_selection_u   <- sparse_res$sparse_tuning_selection_u
  sparse_tuning_selection_nfd <- sparse_res$sparse_tuning_selection_nfd
  sparse_tuning_selection_fd  <- sparse_res$sparse_tuning_selection_fd
  cv_scores_u                 <- sparse_res$cv_scores_u
  cv_scores_nfd               <- sparse_res$cv_scores_nfd
  cv_scores_fd                <- sparse_res$cv_scores_fd
  
  # 2) Smooth tuning (GCV)
  count0 <- (if (is.null(sparse_tuning_u)) 1 else length(sparse_tuning_u)) +
    (if (is.null(sparse_tuning_nfd)) 1 else length(sparse_tuning_nfd)) +
    (if (is.null(sparse_tuning_fd)) 1 else nrow(sparse_tuning_fd))
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
    sparse_tuning_type_u        = sparse_tuning_type_u,
    sparse_tuning_type_nfd      = sparse_tuning_type_nfd,
    sparse_tuning_type_fd       = sparse_tuning_type_fd,
    smooth_tuning             = smooth_tuning,
    CV_score_smooth           = Inf,
    power_type                = "sequential",
    pb                         = pb,
    count                      = count0
  )
  smooth_tuning_selection <- smooth_res$smooth_tuning_selection
  index_selection         <- smooth_res$index_selection
  gcv_scores              <- smooth_res$gcv_scores
  
  close(pb)
  list(
    sparse_tuning_selection_u   = sparse_tuning_selection_u,
    sparse_tuning_selection_nfd = sparse_tuning_selection_nfd,
    sparse_tuning_selection_fd  = sparse_tuning_selection_fd,
    cv_scores_u                 = cv_scores_u,
    cv_scores_nfd               = cv_scores_nfd,
    cv_scores_fd                = cv_scores_fd,
    gcv_scores                  = gcv_scores,
    index_selection             = index_selection
  )
}

# Function for gcv_joint
gcv_joint_hybrid <- function(fdata, nfdata, G_half, G, S_smooth, S_2_inverse, G_half_inverse, hd_obj, smooth_tuning, n) {
  mvmfd_obj <- hd_obj$mf
  CV_score_smooth <- Inf
  pb <- txtProgressBar(min = 0, max = dim(smooth_tuning)[1], style = 3, width = 50, char = "=")
  smooth_tuning_result <- handle_smooth_tuning_hybrid(fdata, nfdata, G_half, G, S_smooth, S_2_inverse, G_half_inverse, hd_obj, smooth_tuning = smooth_tuning, CV_score_smooth = CV_score_smooth, power_type = "joint", n = n, pb = pb, count = 0)
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
                                    smooth_tuning_type, 
                                    sparse_tuning_u, 
                                    sparse_tuning_nfd,
                                    sparse_tuning_fd,
                                    sparse_tuning_type_u,
                                    sparse_tuning_type_nfd, 
                                    sparse_tuning_type_fd,
                                    centerfns, 
                                    alpha_orth, 
                                    K_fold_u,
                                    K_fold_nfd,
                                    K_fold_fd,
                                    sparse_CV, 
                                    smooth_GCV, 
                                    penalize_nfd = FALSE,
                                    penalize_fd = FALSE,
                                    penalize_u = FALSE,
                                    n_cores) {
  cl <- parallel::makeCluster(
    n_cores,
    outfile = if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
  )
  on.exit({
    parallel::stopCluster(cl)
    closeAllConnections()
  }, add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(
      library(MHPCA, quietly = TRUE, warn.conflicts = FALSE)
    )
  })

  
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
    C <- do.call("rbind",mvmfd_obj$coefs) 
    p <- mvmfd_obj$nvar 
    smooth_penalty <-  MHPCA:::pen_fun(mvmfd_obj, type = smooth_tuning_type) 
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
  #########matrix input of smoothing parameters###########
  if (smooth_GCV == FALSE) {
    fv_total = c()
    nfv_total <- c()
    # I_a <- D <- S_2 <- S_smooth <- S_2_inverse <- list()
    S_smooth <- S_2_inverse <- list()
    GCV_score = c()
    if(sparse_CV == FALSE){
      CV_score_u <- CV_score_nfd <- CV_score_fd <- c()
    } else{
      CV_score_u <- CV_score_nfd <- CV_score_fd <- list()
    }
    for (i in 1:n) {
      cat(sprintf("Computing the %s PC...\n", ordinal_msg(i)))
      if (is.null(smooth_tuning)) {
        smooth_tuning_temp = if (!is.null(mvmfd_obj)) expand.grid(lapply(rep(0,mvmfd_obj$nvar), function(x) x[1])) else expand.grid(lapply(rep(0,nf_obj$nvar), function(x) x[1]))
      } else{
        smooth_tuning_temp = expand.grid(lapply(smooth_tuning, function(x) x[i]))
      }
      if (i == 1) {
        if (!is.null(mvmfd_obj)) C_temp = C
        #if (!is.null(nf_obj)) nf_data_temp <- nf_data
        nf_data_temp <- nf_data
      } 
      else{
        if (!is.null(mvmfd_obj)) {
          b_original = t(C_temp)%*%u
          C_temp = C_temp - u%*%t(b_original)
        } else {
          b_original = NULL
          C_temp = NULL
        }
        
        if (!is.null(nf_data)) nf_data_temp <- nf_data_temp - u %*% t(t(nf_data_temp)%*%u) else nf_data_temp <- NULL
      }
      if (!is.null(mvmfd_obj)){
        D <- MHPCA:::I_alpha(mvmfd_obj, smooth_tuning_temp) %*% smooth_penalty
        S_2 <- solve(G + D)
        S_2_inverse[[1]] = solve(S_2)
        S_smooth[[1]] <- G_half %*% (S_2) %*% G_half
      } else {
        D <- NULL
        S_2 <- NULL
        S_2_inverse[[1]] = NULL
        S_smooth[[1]] <- NULL
      }
      
      
      if (!is.null(sparse_tuning_u)) {
        sparse_tuning_temp_u <- if (sparse_CV == FALSE) sparse_tuning_u[i] else sparse_tuning_u
      }
      if (!is.null(sparse_tuning_nfd)) {
        sparse_tuning_temp_nfd <- if (sparse_CV == FALSE) sparse_tuning_nfd[i] else sparse_tuning_nfd
      }
      # if (!is.null(sparse_tuning_fd)) {
      #   sparse_tuning_temp_fd <- if (sparse_CV == FALSE) sparse_tuning_fd[i] else sparse_tuning_fd
      # }
      if (!is.null(sparse_tuning_fd)) {
        if (sparse_CV == FALSE) {
          # Convert the list/vector into a matrix and extract the i-th row as a 1-row matrix.
          #smooth_tuning_temp = expand.grid(lapply(smooth_tuning, function(x) x[i]))
          sparse_tuning_temp_fd <- as.matrix(expand.grid(lapply(sparse_tuning_fd,function(x) x[i])))
          #sparse_tuning_temp_fd <- as.matrix(expand.grid(sparse_tuning_fd))[i, , drop = FALSE]
        } else {
          sparse_tuning_temp_fd <- as.matrix(expand.grid(sparse_tuning_fd))
        }
      }
      
      
      cv_result = cv_gcv_sequential_hybrid(
        fdata = C_temp, 
        nfdata = nf_data,
        hd_obj = hd_obj, 
        smooth_tuning = if (is.null(smooth_tuning)) smooth_tuning else smooth_tuning_temp, 
        sparse_tuning_u = if (is.null(sparse_tuning_u)) sparse_tuning_u else sparse_tuning_temp_u, 
        sparse_tuning_nfd = if (is.null(sparse_tuning_nfd)) sparse_tuning_nfd else sparse_tuning_temp_nfd, 
        sparse_tuning_fd = if (is.null(sparse_tuning_fd)) sparse_tuning_fd else sparse_tuning_temp_fd,
        sparse_tuning_type_u = sparse_tuning_type_u, 
        sparse_tuning_type_nfd = sparse_tuning_type_nfd,
        sparse_tuning_type_fd = sparse_tuning_type_fd,
        K_fold_u = K_fold_u, 
        K_fold_nfd = K_fold_nfd,
        K_fold_fd = K_fold_fd,
        G = G, 
        G_half = G_half, 
        G_half_inverse = G_half_inverse, 
        S_smooth = S_smooth, 
        S_2_inverse = S_2_inverse, 
        penalize_nfd = penalize_nfd,
        penalize_fd = penalize_fd,
        penalize_u = penalize_u, 
        cl = cl
      )
      sparse_result_u = cv_result$sparse_tuning_selection_u
      sparse_result_nfd = cv_result$sparse_tuning_selection_nfd
      sparse_result_fd = cv_result$sparse_tuning_selection_fd
      smooth_result_index = cv_result$index_selection
      if (sparse_CV == FALSE) {
        CV_score_u = c(CV_score_u, cv_result$cv_scores_u)
        CV_score_nfd = c(CV_score_nfd, cv_result$cv_scores_nfd)
        CV_score_fd = c(CV_score_fd, cv_result$cv_scores_fd)
      } else{
        CV_score_u[[i]] = cv_result$cv_scores_u
        CV_score_nfd[[i]] = cv_result$cv_scores_nfd
        CV_score_fd[[i]] = cv_result$cv_scores_fd
      }
      GCV_score = c(GCV_score, cv_result$gcv_scores)
      CG_temp <- if(!is.null(mvmfd_obj)) C_temp%*%G_half else NULL
      
      test_result = init_sequential_hybrid(fdata = CG_temp, 
                                           nfdata = nf_data_temp,
                                           sparse_tuning_result_u = sparse_result_u, 
                                           sparse_tuning_result_nfd = sparse_result_nfd,
                                           sparse_tuning_result_fd = sparse_result_fd,
                                           sparse_tuning_type_u = sparse_tuning_type_u,
                                           sparse_tuning_type_nfd = sparse_tuning_type_nfd,
                                           sparse_tuning_type_fd = sparse_tuning_type_fd,
                                           S_smooth = S_smooth[[1]], 
                                           S_2_inverse = S_2_inverse[[1]], 
                                           G_half_inverse = G_half_inverse, 
                                           G_half = G_half, 
                                           penalize_nfd = penalize_nfd,
                                           penalize_fd = penalize_fd,
                                           penalize_u =  penalize_u)
      u = test_result[[3]]
      fv = test_result[[1]]
      nfv <- test_result[[2]]
      smooth_tuning_result[[i]] = smooth_tuning_temp[smooth_result_index, ]
      sparse_tuning_result_u[[i]] = sparse_result_u
      sparse_tuning_result_nfd[[i]] = sparse_result_nfd
      sparse_tuning_result_fd[[i]] = sparse_result_fd
      temp_count <- 0
      variance_result <- handle_variance_update_hybrid(i, n, C, nf_data,G, fv_total, nfv_total,hd_obj, all_equal_check, c(sparse_tuning_u,sparse_tuning_nfd,sparse_tuning_fd), pc, lsv, fv, nfv,u, G_half, test_result, temp_count, p)
      pc <- variance_result$pc
      lsv <- variance_result$lsv
      fv_total <- variance_result$fv_total
      nfv_total <- variance_result$nfv_total
      variance[i] <- variance_result$variance
    }
  } 
  #########sequential inputs of smoothing parameters###########
  else{
    
    if (is.null(smooth_tuning)) {
      smooth_tuning_temp = if (!is.null(mvmfd_obj)) expand.grid(lapply(rep(0,mvmfd_obj$nvar), function(x) x[1])) else expand.grid(lapply(rep(0,nf_obj$nvar), function(x) x[1]))
    } else{
      smooth_tuning_temp <- expand.grid(smooth_tuning)
    }
    S_smooth <- S_2_inverse <- list()
    cat("Preprocessing...\n")
    
    n_iter1 <- dim(smooth_tuning_temp)[1]
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = n_iter1, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
    #####Do the following computation in advance to save computational cost#####
    
    if (!is.null(mvmfd_obj)){
      for (smooth_index in 1:dim(smooth_tuning_temp)[1]) {
        setTxtProgressBar(pb, smooth_index)
        D <- MHPCA:::I_alpha(mvmfd_obj, smooth_tuning_temp[smooth_index, ]) %*% smooth_penalty
        S_2 <- solve(G + D)
        S_2_inverse[[smooth_index]] = solve(S_2)
        S_smooth[[smooth_index]] <- G_half %*% (S_2) %*% G_half
        
      }
    } else {
      setTxtProgressBar(pb, 1)
      D <- NULL
      S_2 <- NULL
      S_2_inverse = NULL
      S_smooth <- NULL
    }
    close(pb)
    
    fv_total = c()
    nfv_total = c()
    GCV_score = list()
    if(sparse_CV == FALSE){
      CV_score_u <- CV_score_nfd <- CV_score_fd <-  c()
    } else{
      CV_score_u <- CV_score_nfd <- CV_score_fd <- list()
    }
    
    for (i in 1:n) {
      cat(sprintf("Computing the %s PC...\n", ordinal_msg(i)))
      if (i == 1) {
        C_temp = C
        nf_data_temp <- nf_data
      } 
      else{
        if (!is.null(mvmfd_obj)){
          b_original = t(C_temp)%*%u
          C_temp = C_temp - u%*%t(b_original)
          
        } else {
          b_original = NULL
          C_temp = NULL
          
        }
        if (!is.null(nf_data)) nf_data_temp <- nf_data_temp - u %*% t(t(nf_data_temp)%*%u)
      }
      
      if (!is.null(sparse_tuning_u)) {
        sparse_tuning_temp_u <- if (sparse_CV == FALSE) sparse_tuning_u[i] else sparse_tuning_u
      }
      if (!is.null(sparse_tuning_nfd)) {
        sparse_tuning_temp_nfd <- if (sparse_CV == FALSE) sparse_tuning_nfd[i] else sparse_tuning_nfd
      }
      if (!is.null(sparse_tuning_fd)) {
        if (sparse_CV == FALSE) {
          # Convert the list/vector into a matrix and extract the i-th row as a 1-row matrix.
          sparse_tuning_temp_fd <- as.matrix(expand.grid(lapply(sparse_tuning_fd,function(x) x[i])))
        } else {
          sparse_tuning_temp_fd <- as.matrix(expand.grid(sparse_tuning_fd))
        }
      }
      
      cv_result = cv_gcv_sequential_hybrid(
        fdata = C_temp,
        nfdata = nf_data,
        hd_obj = hd_obj, 
        smooth_tuning = if (is.null(smooth_tuning)) smooth_tuning else smooth_tuning_temp, 
        sparse_tuning_u = if (is.null(sparse_tuning_u)) sparse_tuning_u else sparse_tuning_temp_u, 
        sparse_tuning_nfd = if (is.null(sparse_tuning_nfd)) sparse_tuning_nfd else sparse_tuning_temp_nfd,
        sparse_tuning_fd = if (is.null(sparse_tuning_fd)) sparse_tuning_fd else sparse_tuning_temp_fd,
        sparse_tuning_type_u = sparse_tuning_type_u, 
        sparse_tuning_type_nfd = sparse_tuning_type_nfd,
        sparse_tuning_type_fd = sparse_tuning_type_fd,
        K_fold_u = K_fold_u, 
        K_fold_nfd = K_fold_nfd,
        K_fold_fd = K_fold_fd,
        G = G, 
        G_half = G_half, 
        G_half_inverse = G_half_inverse, 
        S_smooth = S_smooth, 
        S_2_inverse = S_2_inverse, 
        penalize_nfd = penalize_nfd,
        penalize_fd = penalize_fd,
        penalize_u = penalize_u,
        cl = cl
      )
      
      sparse_result_u = cv_result$sparse_tuning_selection_u
      sparse_result_nfd = cv_result$sparse_tuning_selection_nfd
      sparse_result_fd = cv_result$sparse_tuning_selection_fd
      smooth_result_index = cv_result$index_selection
      if (sparse_CV == FALSE) {
        CV_score_u = c(CV_score_u, cv_result$cv_scores_u)
        CV_score_nfd = c(CV_score_nfd, cv_result$cv_scores_nfd)
        CV_score_fd = c(CV_score_fd, cv_result$cv_scores_fd)
      } else{
        CV_score_u[[i]] = cv_result$cv_scores_u
        CV_score_nfd[[i]] = cv_result$cv_scores_nfd
        CV_score_fd[[i]] = cv_result$cv_scores_fd
      }
      GCV_score[[i]] = cv_result$gcv_scores
      CG_temp <- if(!is.null(mvmfd_obj)) C_temp%*%G_half else NULL
      
      test_result = init_sequential_hybrid(fdata = CG_temp, 
                                           nfdata =nf_data_temp,
                                           sparse_tuning_result_u = sparse_result_u, 
                                           sparse_tuning_result_nfd = sparse_result_nfd,
                                           sparse_tuning_result_fd = sparse_result_fd,
                                           sparse_tuning_type_u = sparse_tuning_type_u, 
                                           sparse_tuning_type_nfd = sparse_tuning_type_nfd,
                                           sparse_tuning_type_fd = sparse_tuning_type_fd,
                                           S_smooth = S_smooth[[smooth_result_index]], 
                                           S_2_inverse = S_2_inverse[[smooth_result_index]], 
                                           G_half_inverse = G_half_inverse, 
                                           G_half = G_half, 
                                           penalize_nfd = penalize_nfd,
                                           penalize_fd = penalize_fd, 
                                           penalize_u = penalize_u)

      u = test_result[[3]]
      fv = test_result[[1]]
      nfv = test_result[[2]]
      smooth_tuning_result[[i]] = smooth_tuning_temp[smooth_result_index, ]
      sparse_tuning_result_u[[i]] = sparse_result_u
      sparse_tuning_result_nfd[[i]] = sparse_result_nfd
      sparse_tuning_result_fd[[i]] = sparse_result_fd
      temp_count <- 0
      variance_result <- handle_variance_update_hybrid(i, n, C, nf_data,G, fv_total, nfv_total,hd_obj, all_equal_check, c(sparse_tuning_u,sparse_tuning_nfd,sparse_tuning_fd), pc, lsv, fv, nfv,u, G_half, test_result, temp_count, p)
      pc <- variance_result$pc
      lsv <- variance_result$lsv
      fv_total <- variance_result$fv_total
      nfv_total <- variance_result$nfv_total
      variance[i] <- variance_result$variance
    }
    if (is.null(smooth_tuning)) {
      GCV_score = NULL
    }
    if (is.null(sparse_tuning_u) && is.null(sparse_tuning_nfd)) {
      CV_score = NULL
    }
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
    CV_score_u = CV_score_u,
    CV_score_nfd = CV_score_nfd, 
    CV_score_fd = CV_score_fd, 
    GCV_score = GCV_score
    ))
} 



# joint smooth and sparse power algorithm
joint_power_hybrid <- function(hd_obj, n, smooth_tuning, smooth_tuning_type, centerfns, alpha_orth) {
  
  print("*****")
  #######centralize########
  if (centerfns) hd_obj <- center_hd(hd_obj)
  mvmfd_obj <- hd_obj$mf
  nf_obj <- hd_obj$nf
  nf_data <- if (!is.null(nf_obj)) do.call("cbind",nf_obj$data) else NULL
  if (!is.null(mvmfd_obj)){
    p <- mvmfd_obj$nvar
    smooth_penalty <- MHPCA:::pen_fun(mvmfd_obj, type = smooth_tuning_type)
    C <- do.call("rbind",mvmfd_obj$coefs)  
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
    cat("Preprocessing...\n")
    
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
      D <- MHPCA:::I_alpha(mvmfd_obj, smooth_tuning_temp[smooth_index, ]) %*% smooth_penalty
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
    cat(sprintf("Computing PCs...\n"))
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
    S_2_inverse = S_2_inverse, n = n
  )
  smooth_result_index = cv_result[[2]]
  GCV_score = cv_result[[3]]
  CG_temp <- if (!is.null(mvmfd_obj)) C%*%G_half else NULL
  test_result = init_joint_hybrid(CG_temp, nf_data,S_smooth[[smooth_result_index]], S_2_inverse[[smooth_result_index]], G_half_inverse, G_half, n = n)
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
