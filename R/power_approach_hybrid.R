#' @importFrom expm sqrtm
#' @importFrom utils  txtProgressBar setTxtProgressBar


# sequential power algorithm
init_sequential_hybrid <- function(fdata,
                                   nfdata,
                                   sparse_tuning_result,
                                   sparse_tuning_type,
                                   S_smooth = NULL,
                                   S_2_inverse = NULL,
                                   G_half_inverse = NULL,
                                   G_half = NULL,
                                   cv_flag = FALSE,
                                   penalize_u = FALSE,
                                   penalize_nfd = FALSE) {
  fv_old <- if (!is.null(fdata)) svd(fdata)$v[, 1] else NULL
  nfv_old <- if (!is.null(nfdata)) svd(nfdata)$v[, 1] else NULL
  errors <- 10^60
  while (errors > 0.00001) {
    y <- 0
    if (!is.null(fdata)) {
      y <- y + fdata %*% fv_old
    }
    if (!is.null(nfdata)) {
      y <- y + nfdata %*% nfv_old
    }

    u_old <- if (penalize_u) {
      sparse_pen_fun(y = y, tuning_parameter = sparse_tuning_result, sparse_tuning_type)
    } else {
      y
    }
    u_old <- u_old / norm(u_old, type = "2")
    if (cv_flag == TRUE) { # incorporate power algorithm in CV sparse tuning selection

      fv_new <- if (!is.null(fdata)) t(fdata) %*% u_old else NULL

      if (!is.null(nfdata)) {
        nfv_new <- t(nfdata) %*% u_old
        nfv_new <- if (penalize_nfd) sparse_pen_fun(y = nfv_new, tuning_parameter = sparse_tuning_result, type = sparse_tuning_type) else nfv_new
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

      y <- 0
      if (!is.null(fdata)) {
        y <- y + fdata %*% fv_new
      }
      if (!is.null(nfdata)) {
        y <- y + nfdata %*% nfv_new
      }

      u_new <- if (penalize_u) {
        sparse_pen_fun(y = y, tuning_parameter = sparse_tuning_result, sparse_tuning_type)
      } else {
        y
      }
      u_new <- u_new / norm(u_new, type = "2")
      errors <- sum((u_new - u_old)^2)
      if (!is.null(fdata)) fv_old <- fv_new
      if (!is.null(nfdata)) nfv_old <- nfv_new
      u_old <- u_new
    } else { # incorporate power algorithm smoothing
      fv_new <- if (!is.null(fdata)) S_smooth %*% t(fdata) %*% u_old else NULL
      if (!is.null(nfdata)) {
        nfv_new <- t(nfdata) %*% u_old
        nfv_new <- if (penalize_nfd) sparse_pen_fun(y = nfv_new, tuning_parameter = sparse_tuning_result, type = sparse_tuning_type) else nfv_new
      } else {
        nfv_new <- NULL
      }

      # normalizing fv and nfv
      ##########################################
      if (!is.null(fdata)) {
        fv_new_back <- G_half_inverse %*% fv_new
        fv_weight <- as.numeric(sqrt(t(fv_new_back) %*% S_2_inverse %*% fv_new_back))^2
      } else {
        fv_weight <- 0
      }
      nfv_weight <- if (!is.null(nfdata)) norm(nfv_new, "2")^2 else 0
      coef_norm <- sqrt(fv_weight + nfv_weight)
      if (!is.null(fdata)) {
        fv_new_back <- fv_new_back / coef_norm
        fv_new <- G_half %*% fv_new_back
      }
      if (!is.null(nfdata)) {
        nfv_new <- nfv_new / coef_norm
      }
      ##########################################
      y <- 0
      if (!is.null(fdata)) {
        y <- y + fdata %*% fv_new
      }
      if (!is.null(nfdata)) {
        y <- y + nfdata %*% nfv_new
      }

      u_new <- if (penalize_u) {
        sparse_pen_fun(y = y, tuning_parameter = sparse_tuning_result, sparse_tuning_type)
      } else {
        y
      }
      u_new <- u_new / norm(u_new, type = "2")

      errors <- sum((u_new - u_old)^2)
      if (!is.null(fdata)) fv_old <- fv_new
      if (!is.null(nfdata)) nfv_old <- nfv_new
    }
  }
  if (cv_flag == TRUE) {
    return(u_new)
  } else {
    if (is.null(nfdata)) nfv_new <- NULL
    if (!is.null(fdata)) fv_new <- G_half_inverse %*% fv_new else fv_new <- NULL
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

#computing cv score for sparse tuning
cv_local_hybrid = function(fdata, nfdata, G_half, K_fold, sparse_tuning_single, sparse_tuning_type, shuffled_row, group_size, penalize_nfd = FALSE, penalize_u = FALSE){
  
  shuffled_row_f <- shuffled_row$f
  shuffled_row_nf <- shuffled_row$nf
  group_size_f <- group_size$f
  group_size_nf <- group_size$nf
  data_double_tilde = if (!is.null(fdata)) t(fdata%*%G_half) else NULL
  error_score_sparse = 0
  for (k in 1:K_fold) {
    rows_to_remove_f <- if (!is.null(fdata)) shuffled_row_f[((k-1)*group_size_f+1):((k)*group_size_f)]
    rows_to_remove_nf <- if (!is.null(nfdata)) shuffled_row_nf[((k-1)*group_size_nf+1):((k)*group_size_nf)]
    fdata_train = if (!is.null(fdata)) data_double_tilde[-rows_to_remove_f, ,drop = FALSE] else NULL
    fdata_test = if (!is.null(fdata)) data_double_tilde[rows_to_remove_f, ,drop = FALSE] else NULL
    nfdata_train = if (!is.null(nfdata)) nfdata[,-rows_to_remove_nf , drop = FALSE] else NULL
    nfdata_test = if (!is.null(nfdata)) nfdata[,rows_to_remove_nf , drop = FALSE] else NULL
    u_test = init_sequential_hybrid(t(fdata_train), nfdata_train ,sparse_tuning_single, sparse_tuning_type, cv_flag = TRUE, penalize_nfd = penalize_nfd, penalize_u = penalize_u)
    fv_test = if (!is.null(fdata)) fdata_test%*%u_test else NULL
    nfv_test = if (!is.null(nfdata)) t(nfdata_test)%*%u_test else NULL
    fv_test_smooth_back = if (!is.null(fdata)) (data_double_tilde%*%u_test)[rows_to_remove_f, , drop = FALSE] else NULL
    fdata_test_smooth_back = if (!is.null(fdata)) t(data_double_tilde)[, rows_to_remove_f, drop = FALSE] else NULL
    fv_error <- if (!is.null(fdata)) sum((t(fdata_test_smooth_back)-fv_test_smooth_back%*%t(u_test))^2) else 0
    nfv_error <- if (!is.null(nfdata)) sum((t(nfdata_test)-nfv_test%*%t(u_test))^2) else 0
    error_score_sparse = error_score_sparse + fv_error + nfv_error
  }
  normalization_factor <- if (!is.null(fdata)) {
    ncol(fdata)
  } else if (!is.null(nfdata)) {
    nrow(nfdata)
  } else {
    1  
  }
  return(error_score_sparse / normalization_factor)

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
handle_smooth_tuning_hybrid <- function(fdata, nfdata, G_half, G, S_smooth, S_2_inverse, G_half_inverse, 
                                        hd_obj, sparse_tuning_selection = NULL, sparse_tuning_type = NULL, smooth_tuning, 
                                        CV_score_smooth, power_type = "sequential", n = NULL, pb, 
                                        count, penalize_nfd = F, penalize_u = FALSE) {
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
          test_temp <- init_sequential_hybrid(fdata %*% G_half, nfdata, sparse_tuning_selection, sparse_tuning_type, S_smooth[[smooth_index]], 
                                              S_2_inverse[[smooth_index]], G_half_inverse, G_half,penalize_nfd = penalize_nfd, penalize_u = penalize_u)
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

# Function to handle sparse tuning selection
handle_sparse_tuning_hybrid <- function(fdata, nfdata,G_half, sparse_tuning, sparse_tuning_type, K_fold, shuffled_row, group_size, CV_score_sparse, pb, penalize_nfd = FALSE, penalize_u = FALSE) {
  count <- 0
  cv_scores <- c()
  sparse_tuning_selection <- NULL
  #browser()
  if (is.null(sparse_tuning)) {
    count <- count + 1
    setTxtProgressBar(pb, count)
    cv_scores <- NULL
    sparse_tuning_selection <- 0
  } else {
    for (sparse_tuning_single in sparse_tuning) {
      count <- count + 1
      setTxtProgressBar(pb, count)
      sparse_score <- cv_local_hybrid(fdata, nfdata, G_half, K_fold, sparse_tuning_single, sparse_tuning_type, shuffled_row, group_size, penalize_nfd = penalize_nfd, penalize_u = penalize_u)
      cv_scores <- c(cv_scores, sparse_score)
      if (sparse_score < CV_score_sparse) {
        CV_score_sparse <- sparse_score
        sparse_tuning_selection <- sparse_tuning_single
      }
    }
  }
  return(list(sparse_tuning_selection = sparse_tuning_selection, cv_scores = cv_scores, CV_score_sparse = CV_score_sparse))
}

# Function for cv_gcv_sequential
cv_gcv_sequential_hybrid <- function(fdata, nfdata, hd_obj, smooth_tuning, sparse_tuning, sparse_tuning_type, K_fold, G, G_half, G_half_inverse, S_smooth, S_2_inverse, penalize_nfd = FALSE, penalize_u = FALSE) {
  mvmfd_obj <- hd_obj$mf
  CV_score_sparse <- CV_score_smooth <- Inf
  result <- c()
  count <- 0
  nc <- 0 
  if (!is.null(fdata)) ncf <- ncol(fdata)
  if (!is.null(nfdata)) ncnf <- ncol(nfdata)
  
  shuffled_row_f <- sample(ncf)
  len_f <- if (!is.null(fdata)) length(shuffled_row_f) else 0
  shuffled_row_nf <- sample(ncnf)
  len_nf <- if (!is.null(nfdata)) length(shuffled_row_nf) else 0
  group_size_f <- len_f / K_fold 
  group_size_nf <- len_nf / K_fold
  group_size <- list(f = group_size_f,nf = group_size_nf)
  shuffled_row <- list(f = shuffled_row_f, nf = shuffled_row_nf)
  n_iter <- (if (is.null(smooth_tuning)) 1 else dim(smooth_tuning)[1]) + (if (is.null(sparse_tuning)) 1 else length(sparse_tuning))
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  
  # Handle sparse tuning
  sparse_tuning_result <- handle_sparse_tuning_hybrid(fdata, nfdata, G_half, sparse_tuning, sparse_tuning_type, K_fold, shuffled_row, group_size, CV_score_sparse, pb = pb, penalize_nfd = penalize_nfd, penalize_u = penalize_u)
  sparse_tuning_selection <- sparse_tuning_result$sparse_tuning_selection
  cv_scores <- sparse_tuning_result$cv_scores
  CV_score_sparse <- sparse_tuning_result$CV_score_sparse
  
  # Handle smooth tuning
  smooth_tuning_result <- handle_smooth_tuning_hybrid(fdata, nfdata, G_half, G, S_smooth, S_2_inverse, G_half_inverse, hd_obj, sparse_tuning_selection, sparse_tuning_type, smooth_tuning, CV_score_smooth, power_type = "sequential", pb = pb, count = (if (is.null(sparse_tuning)) 1 else length(sparse_tuning)))
  smooth_tuning_selection <- smooth_tuning_result$smooth_tuning_selection
  index_selection <- smooth_tuning_result$index_selection
  gcv_scores <- smooth_tuning_result$gcv_scores
  
  close(pb)
  result <- list(sparse_tuning_selection, smooth_tuning_selection, index_selection, cv_scores, gcv_scores)
  return(result)
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
  if (i == 1 || all(all_equal_check) || T & (is.null(sparse_tuning) || all(unique(sparse_tuning) == 0))) {
    CGfv <- if (!is.null(mvmfd_obj)) C %*% G %*% fv else NULL
    Xnfv <- if (!is.null(hd_obj$nf)) nf_data %*% nfv else NULL
    if (!is.null(CGfv) && !is.null(Xnfv)) variance <- (t(Xnfv)%*%Xnfv+2*t(Xnfv)%*% CGfv+t(CGfv)%*% CGfv)/(mvmfd_obj$nobs - 1)
    if (!is.null(CGfv) && is.null(Xnfv)) variance <- (t(CGfv)%*% CGfv)/(mvmfd_obj$nobs - 1)
    if (is.null(CGfv) && !is.null(Xnfv)) variance <- (t(Xnfv)%*%Xnfv)/(hd_obj$nf$nobs - 1)
    #variance <- (t(Xnfv)%*%Xnfv+2*t(Xnfv)%*% CGfv+t(CGfv)%*% CGfv)/(mvmfd_obj$nobs - 1)
    # CGv <- C %*% G %*% v
    # variance <- t(CGv)%*%CGv / (mvmfd_obj$nobs - 1)
  } else {
    # it is completely wrong for now, fix it later by writting its math
    CGvfv <- if (!is.null(mvmfd_obj)) C %*% G %*% fv_total
    G_pc = if (!is.null(mvmfd_obj)) t(fv_total)%*%G%*%fv_total 
    coef_pc = CGv %*% solve(G_pc)
    total_variance = sum(diag((coef_pc%*%t(fv_total)) %*% G %*% t(coef_pc%*%t(fv_total))))
    G_pc_pre = t(v_total[, -i])%*%G%*%v_total[, -i]
    coef_pc_pre = CGv[, -i] %*% solve(G_pc_pre)
    total_variance_previous = sum(diag((coef_pc_pre%*%t(v_total[, -i])) %*% G %*% t(coef_pc_pre%*%t(v_total[, -i]))))
    variance = (total_variance - total_variance_previous) / (mvmfd_obj$nobs - 1)
  }
  return(list(pc = pc, lsv = lsv, fv_total = fv_total, nfv_total = nfv_total, variance = variance))
}

# sequential power algorithm
sequential_power_hybrid <- function(hd_obj, n, smooth_tuning, smooth_tuning_type, sparse_tuning, sparse_tuning_type, centerfns, alpha_orth, K_fold, sparse_CV, smooth_GCV, penalize_nfd = FALSE, penalize_u = FALSE) {
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
    sparse_tuning_result  <- list()
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
    sparse_tuning_result  <- list()
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
      CV_score = c()
    } else{
      CV_score = list()
    }
    for (i in 1:n) {
      cat(sprintf("Computing the %s PC...\n", ordinal_msg(i)))
      if (is.null(smooth_tuning)) {
        smooth_tuning_temp = if (!is.null(mvmfd_obj)) expand.grid(lapply(rep(0,mvmfd_obj$nvar), function(x) x[1])) else expand.grid(lapply(rep(0,mvnfd_obj$nvar), function(x) x[1]))
      } else{
        smooth_tuning_temp = expand.grid(lapply(smooth_tuning, function(x) x[i]))
      }
      if (i == 1) {
        if (!is.null(mvmfd_obj)) C_temp = C
        if (!is.null(nf_obj)) nf_data_temp <- nf_data
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
      
      
      if (!is.null(sparse_tuning)) {
        sparse_tuning_temp <- if (sparse_CV == FALSE) sparse_tuning[i] else sparse_tuning
      }
      
      cv_result = cv_gcv_sequential_hybrid(
        fdata = C_temp, 
        nfdata = nf_data,
        hd_obj = hd_obj, 
        smooth_tuning = if (is.null(smooth_tuning)) smooth_tuning else smooth_tuning_temp, 
        sparse_tuning = if (is.null(sparse_tuning)) sparse_tuning else sparse_tuning_temp, 
        sparse_tuning_type = sparse_tuning_type, 
        K_fold = K_fold, 
        G = G, 
        G_half = G_half, 
        G_half_inverse = G_half_inverse, 
        S_smooth = S_smooth, 
        S_2_inverse = S_2_inverse, 
        penalize_nfd = penalize_nfd, penalize_u = penalize_u
      )
      sparse_result = cv_result[[1]]
      smooth_result_index = cv_result[[3]]
      if (sparse_CV == FALSE) {
        CV_score = c(CV_score, cv_result[[4]])
      } else{
        CV_score[[i]] = cv_result[[4]]
      }
      GCV_score = c(GCV_score, cv_result[[5]])
      CG_temp <- if(!is.null(mvmfd_obj)) C_temp%*%G_half else NULL
      test_result = init_sequential_hybrid(CG_temp, nf_data_temp,sparse_result, sparse_tuning_type, S_smooth[[1]], S_2_inverse[[1]], G_half_inverse, G_half, penalize_nfd = penalize_nfd, penalize_u =  penalize_u)
      u = test_result[[3]]
      fv = test_result[[1]]
      nfv <- test_result[[2]]
      smooth_tuning_result[[i]] = smooth_tuning_temp[smooth_result_index, ]
      sparse_tuning_result[[i]] = sparse_result
      temp_count <- 0
      variance_result <- handle_variance_update_hybrid(i, n, C, nf_data,G, fv_total, nfv_total,hd_obj, all_equal_check, sparse_tuning, pc, lsv, fv, nfv,u, G_half, test_result, temp_count, p)
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
      smooth_tuning_temp = if (!is.null(mvmfd_obj)) expand.grid(lapply(rep(0,mvmfd_obj$nvar), function(x) x[1])) else expand.grid(lapply(rep(0,mvnfd_obj$nvar), function(x) x[1]))
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
        close(pb)
      }
    } else {
      D <- NULL
      S_2 <- NULL
      S_2_inverse = NULL
      S_smooth <- NULL
    }
   
    
    fv_total = c()
    nfv_total = c()
    GCV_score = list()
    if(sparse_CV == FALSE){
      CV_score = c()
    } else{
      CV_score = list()
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
      
      if (!is.null(sparse_tuning)) {
        sparse_tuning_temp <- if (sparse_CV == FALSE) sparse_tuning[i] else sparse_tuning
      }
      
      cv_result = cv_gcv_sequential_hybrid(
        fdata = C_temp,
        nfdata = nf_data,
        hd_obj = hd_obj, 
        smooth_tuning = if (is.null(smooth_tuning)) smooth_tuning else smooth_tuning_temp, 
        sparse_tuning = if (is.null(sparse_tuning)) sparse_tuning else sparse_tuning_temp, 
        sparse_tuning_type = sparse_tuning_type, 
        K_fold = K_fold, 
        G = G, 
        G_half = G_half, 
        G_half_inverse = G_half_inverse, 
        S_smooth = S_smooth, 
        S_2_inverse = S_2_inverse, 
        penalize_nfd = penalize_nfd, 
        penalize_u = penalize_u
      )
      
      sparse_result = cv_result[[1]]
      smooth_result_index = cv_result[[3]]
      if (sparse_CV == FALSE) {
        CV_score = c(CV_score, cv_result[[4]])
      } else{
        CV_score[[i]] = cv_result[[4]]
      }
      GCV_score[[i]] = cv_result[[5]]
      CG_temp <- if(!is.null(mvmfd_obj)) C_temp%*%G_half else NULL
      test_result = init_sequential_hybrid(CG_temp, nf_data_temp,sparse_result, sparse_tuning_type, S_smooth[[smooth_result_index]], S_2_inverse[[smooth_result_index]], G_half_inverse, G_half, penalize_nfd = penalize_nfd, penalize_u = penalize_u)

      u = test_result[[3]]
      fv = test_result[[1]]
      nfv = test_result[[2]]
      smooth_tuning_result[[i]] = smooth_tuning_temp[smooth_result_index, ]
      sparse_tuning_result[[i]] = sparse_result
      temp_count <- 0
      variance_result <- handle_variance_update_hybrid(i, n, C, nf_data,G, fv_total, nfv_total,hd_obj, all_equal_check, sparse_tuning, pc, lsv, fv, nfv,u, G_half, test_result, temp_count, p)
      pc <- variance_result$pc
      lsv <- variance_result$lsv
      fv_total <- variance_result$fv_total
      nfv_total <- variance_result$nfv_total
      variance[i] <- variance_result$variance
    }
    if (is.null(smooth_tuning)) {
      GCV_score = NULL
    }
    if (is.null(sparse_tuning)) {
      CV_score = NULL
    }
  }
  return(list(pc_fd = pc, pc_nfd = nfv_total,lsv = lsv, variance = variance, smooth_tuning_result = smooth_tuning_result, sparse_tuning_result = sparse_tuning_result, CV_score = CV_score, GCV_score = GCV_score))
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
  

  # cv_result = gcv_joint(data = C, mvmfd_obj = mvmfd_obj, smooth_tuning = smooth_tuning, G = G, G_half = G_half, G_half_inverse = G_half_inverse, S_smooth = S_smooth, S_2_inverse = S_2_inverse, n = n)
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
