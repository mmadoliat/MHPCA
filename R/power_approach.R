#' @importFrom expm sqrtm
#' @importFrom utils  txtProgressBar setTxtProgressBar
#' 
#sparse penalty function
sparse_pen_fun <- function(y, tuning_parameter, type, alpha = 3.7) {
  y_sorted <- sort(abs(y))
  lambda = y_sorted[tuning_parameter]
  if (tuning_parameter == 0) {
    return(y)
  }
  if (type == "soft") {
    return(sign(y) * pmax(abs(y) - lambda, 0))
  } 
  else if (type == "hard") {
    return(ifelse(abs(y) > lambda, y, 0))
  } 
  else if (type == "SCAD") {
    res <- ifelse(abs(y) <= 2 * lambda,
                  sign(y) * pmax(abs(y) - lambda, 0),
                  ifelse(abs(y) <= alpha * lambda,
                         ((alpha - 1) * y - sign(y) * alpha * lambda) / (alpha - 2),
                         y))
    return(res)
  }
}

#sequential power algorithm 
init_sequential = function(data, sparse_tuning_result, sparse_tuning_type, S_smooth = NULL, S_2_inverse = NULL, G_half_inverse = NULL, G_half = NULL, cv_flag = FALSE){
  b_old = svd(data)$v[, 1]
  errors = 10^60
  while (errors > 0.00001) {
    v_new = sparse_pen_fun(y = data%*%b_old, tuning_parameter = sparse_tuning_result, sparse_tuning_type)
    if (cv_flag == TRUE) { #incorporate power algorithm in CV sparse tuning selection
      b_new = t(data)%*%v_new
      b_new = b_new / norm(b_new, type = "2")
      errors = sum((b_new - b_old)^2)
      b_old = b_new
    } else{ #implement power algorithm to obtain functional pc
      b_new = S_smooth%*%t(data)%*%v_new
      b_new_back = G_half_inverse %*% b_new
      b_new_back = b_new_back/as.numeric(sqrt(t(b_new_back) %*% S_2_inverse %*% b_new_back))
      b_new = G_half%*%b_new_back
      errors = sum((b_new - b_old)^2)
      b_old = b_new
    }
  }
  v_new_rescale = v_new/norm(v_new, type = "2")
  if (cv_flag == TRUE) {
    return(v_new_rescale)
  }
  else{
    b_new = G_half_inverse%*%b_new
    b_new = b_new %*% solve(sqrt(t(b_new) %*% S_2_inverse %*% b_new))
    return(list(b_new, v_new_rescale, v_new))
  }
}

#joint power for smoothing
init_joint = function(data, S_smooth = NULL, S_2_inverse = NULL, G_half_inverse = NULL, G_half = NULL, n = n){
  b_old = svd(data)$v[, 1:n]
  errors = 10^60
  while (errors > 0.00001) {
    v_new = data%*%b_old
    b_new = S_smooth%*%qr.Q(qr(as.matrix(t(data)%*%v_new)))
    b_new_back = G_half_inverse %*% b_new
    b_new_back = b_new_back%*%diag(1/sqrt(diag(t(b_new_back) %*% S_2_inverse %*% b_new_back)))
    b_new = G_half%*%b_new_back
    errors = sum((b_new - b_old)^2)
    b_old = b_new
  }
  v_new_rescale = v_new%*%diag(1/sqrt(diag(t(v_new)%*%v_new)))
  b_new = G_half_inverse%*%b_new
  b_new = b_new %*% diag(1/sqrt(diag(t(b_new) %*% S_2_inverse %*% b_new)))
  return(list(b_new, v_new_rescale, v_new))
}


#computing cv score for sparse tuning
cv_local = function(data, G_half, K_fold, sparse_tuning_single, sparse_tuning_type, shuffled_row, group_size){
  data_double_tilde = t(data%*%G_half)
  error_score_sparse = 0
  for (k in 1:K_fold) {
    rows_to_remove = shuffled_row[((k-1)*group_size+1):((k)*group_size)]
    data_train = data_double_tilde[-rows_to_remove, ]
    data_test = data_double_tilde[rows_to_remove, ]
    v_test = init_sequential(t(data_train), sparse_tuning_single, sparse_tuning_type, cv_flag = TRUE)
    b_test = data_test%*%v_test
    b_test_smooth_back = (data_double_tilde%*%v_test)[rows_to_remove, ]
    data_test_smooth_back = t(data_double_tilde)[, rows_to_remove]
    error_score_sparse = error_score_sparse + sum((t(data_test_smooth_back)-b_test_smooth_back%*%t(v_test))^2)
  }
  return(error_score_sparse/ncol(data))
}

#computing gcv score for smoothing tuning
gcv_local = function(data, hd_obj, G, G_half, S_smooth, v, smooth_tuning) {
  p = hd_obj$nvar
  indices <- sapply(1:p, function(i) prod(hd_obj$basis$nbasis[[i]]))
  B_subtilde = data %*% G_half
  if (all(smooth_tuning == 0)) {
    error_smooth_score <- 0
  } else {
    if (p == 1) {
      error_smooth_score <- (sum(((diag(dim(S_smooth)[1]) - S_smooth) %*% (t(B_subtilde) %*% v))^2) / 
                               ((1 - sum(diag(S_smooth)) / dim(G)[1])^2)) / dim(G)[1]
    } else {
      error_smooth_score <- 0
      start_index <- 1
      for (i in 1:p) {
        end_index <- start_index + indices[i] - 1
        s_alpha_tilde_i <- S_smooth[start_index:end_index, start_index:end_index]
        B_subtilde_i <- B_subtilde[, start_index:end_index]
        error_smooth_score_i <- sum(((diag(indices[i]) - s_alpha_tilde_i) %*% (t(B_subtilde_i) %*% v))^2) / 
          (1 - sum(diag(s_alpha_tilde_i)) / indices[i])^2
        error_smooth_score <- error_smooth_score + error_smooth_score_i
        start_index <- end_index + 1
      }
    }
  }
  return(error_smooth_score)
}



# cv and gcv tuning selection process in sequential power
cv_gcv_sequential = function(data, hd_obj, smooth_tuning, sparse_tuning, sparse_tuning_type, K_fold, G, G_half, G_half_inverse, S_smooth, S_2_inverse){
  CV_score_sparse = CV_score_smooth = 10^60
  result = c()
  count = 0
  shuffled_row = sample(ncol(data))
  group_size = length(shuffled_row) / K_fold
  n_iter <- dim(smooth_tuning)[1] + length(sparse_tuning)
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar

  cv_scores = gcv_scores = c()
  for (sparse_tuning_single in sparse_tuning) {
    count = count +1
    setTxtProgressBar(pb, count)
    sparse_score = cv_local(data, G_half, K_fold, sparse_tuning_single, sparse_tuning_type, shuffled_row, group_size)
    cv_scores = c(cv_scores, sparse_score) 
    if (sparse_score <= CV_score_sparse) {
      CV_score_sparse = sparse_score
      sparse_tuning_selection = sparse_tuning_single
    }
  }
  for (smooth_index in 1:dim(smooth_tuning)[1]) {
    count = count + 1
    setTxtProgressBar(pb, count)
    if (all(smooth_tuning == 0)) {
      S_smooth[[smooth_index]] = diag(dim(G)[1])
    }
    test_temp = init_sequential(data%*%G_half, sparse_tuning_selection, sparse_tuning_type, S_smooth[[smooth_index]], S_2_inverse[[smooth_index]], G_half_inverse, G_half)
    v_temp = test_temp[[2]]
    smooth_score = gcv_local(data, hd_obj, G, G_half, S_smooth[[smooth_index]], v_temp, smooth_tuning=smooth_tuning[smooth_index, ])
    gcv_scores = c(gcv_scores,smooth_score)
    if (smooth_score <= CV_score_smooth) {
      CV_score_smooth = smooth_score
      smooth_tuning_selection = smooth_tuning[smooth_index, ]
      index_selection = smooth_index
    }
  }
  # }
  close(pb) # Close the connection
  result = list(sparse_tuning_selection, smooth_tuning_selection, index_selection, cv_scores, gcv_scores)
  return(result)
}

# gcv tuning selection process in joint power
gcv_joint = function(data, hd_obj, smooth_tuning, G, G_half, G_half_inverse, S_smooth, S_2_inverse, n){
  CV_score_smooth = 10^60
  result = c()
  count = 0
  n_iter <- dim(smooth_tuning)[1]
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  gcv_scores = c()
  for (smooth_index in 1:dim(smooth_tuning)[1]) {
  count = count + 1
  setTxtProgressBar(pb, count)
  if (all(smooth_tuning == 0)) {
    S_smooth[[smooth_index]] = diag(dim(G)[1])
  }
  test_temp = init_joint(data%*%G_half, S_smooth[[smooth_index]], S_2_inverse[[smooth_index]], G_half_inverse, G_half, n = n)
  v_temp = test_temp[[2]]
  smooth_score = gcv_local(data, hd_obj, G, G_half, S_smooth[[smooth_index]], v_temp, smooth_tuning=smooth_tuning[smooth_index, ])
  gcv_scores = c(gcv_scores, smooth_score)
  if (smooth_score <= CV_score_smooth) {
    CV_score_smooth = smooth_score
    smooth_tuning_selection = smooth_tuning[smooth_index, ]
    index_selection = smooth_index
  }
  }
  close(pb) # Close the connection
  result = list(smooth_tuning_selection, index_selection, gcv_scores)
  return(result)
}


ordinal_msg <- function(i) {
  if (i == 1) {
    return(paste0(i, "st"))
  } else if (i == 2) {
    return(paste0(i, "nd"))
  } else if (i == 3) {
    return(paste0(i, "rd"))
  }
  else {
    return(paste0(i, "th"))
  }
}

# Define a function to process hd_obj
centralized_hd <- function(hd_obj, centerfns = TRUE) {
  p <- hd_obj$nvar
  B <- c()
  
  for (i in 1:p) {
    coefs_i <- hd_obj$coefs[[i]]
    
    if (is.matrix(coefs_i)) {
      c <- if (centerfns) coefs_i - rowMeans(coefs_i) else coefs_i
    } else {
      cc <- apply(coefs_i, 3, as.vector)
      c <- if (centerfns) cc - rowMeans(cc) else cc
    }
    
    B <- rbind(B, c)
  }
  
  return(B)
}

# sequential power algorithm
sequential_power <- function(hd_obj, n, smooth_tuning, smooth_tuning_type, sparse_tuning, sparse_tuning_type, centerfns, alpha_orth, K_fold, sparse_CV, smooth_GCV) {
  p <- hd_obj$nvar
  smooth_penalty <- MHPCA:::pen_fun(hd_obj, type = smooth_tuning_type)

  #######centralize########
  B <- centralized_hd(hd_obj, centerfns)
  ########some initial setting#######
  B <- t(B)
  lsv <- c()
  pc <- list()
  variance <- vector()
  smooth_tuning_result  <- list()
  sparse_tuning_result  <- list()
  G <- as.matrix(hd_obj$basis$gram)
  G_half <- sqrtm(G)
  G_half_inverse = solve(G_half)
  if (is.null(smooth_tuning)) {
    for (i in 1:p) {
      smooth_tuning <- c(smooth_tuning, list(2^seq(-20, 20, length.out = 10)))
    }
  }
  
  #########matrix input of smoothing parameters###########
  if (smooth_GCV == FALSE) {
    b_total = v_total = c()
    I_a <- D <- S_2 <- S_smooth <- S_2_inverse <- list()
    GCV_score = c()
    if(sparse_CV == FALSE){
      CV_score = c()
    } else{
      CV_score = list()
    }
    for (i in 1:n) {
      cat(sprintf("Computing the %s PC...\n", ordinal_msg(i)))
      smooth_tuning_temp = expand.grid(lapply(smooth_tuning, function(x) x[i]))
      if (i == 1) {
        B_temp = B
      } 
      else{
        b_original = t(B_temp)%*%v
        B_temp = B_temp - v%*%t(b_original)
      }
      I_a[[1]] <- MHPCA:::I_alpha(hd_obj, smooth_tuning_temp)
      D[[1]] <- I_a[[1]] %*% smooth_penalty
      S_2[[1]] <- solve(G + D[[1]])
      S_2_inverse[[1]] = solve(S_2[[1]])
      S_smooth[[1]] <- G_half %*% (S_2[[1]]) %*% G_half
      
      if(sparse_CV == FALSE){
        sparse_tuning_temp = sparse_tuning[i]
      } else{
        sparse_tuning_temp = sparse_tuning
      }
      cv_result = cv_gcv_sequential(data = B_temp, hd_obj = hd_obj, smooth_tuning = smooth_tuning_temp, sparse_tuning = sparse_tuning_temp, sparse_tuning_type = sparse_tuning_type, K_fold = K_fold, G = G, G_half = G_half, G_half_inverse = G_half_inverse, S_smooth = S_smooth, S_2_inverse = S_2_inverse)
      sparse_result = cv_result[[1]]
      smooth_result_index = cv_result[[3]]
      if (sparse_CV == FALSE) {
        CV_score = c(CV_score, cv_result[[4]])
      } else{
        CV_score[[i]] = cv_result[[4]]
      }
      GCV_score = c(GCV_score, cv_result[[5]])
      test_result = init_sequential(B_temp%*%G_half, sparse_result, sparse_tuning_type, S_smooth[[1]], S_2_inverse[[1]], G_half_inverse, G_half)
      
      v = test_result[[2]]
      b = test_result[[1]]
      smooth_tuning_result[[i]] = smooth_tuning_temp[smooth_result_index, ]
      sparse_tuning_result[[i]] = sparse_result
      temp_count <- 0
      for (j in 1:p) {
        index_start <- (temp_count + 1)
        index_end <- (temp_count + prod(hd_obj$basis$nbasis[[j]]))
        if (i == 1) {
          pc[[j]] <- b[index_start:index_end, ]
        } else {
          pc[[j]] <- cbind(pc[[j]], b[index_start:index_end, ])
        }
        temp_count <- index_end
      }
      lsv = cbind(lsv, v)
      b_total = cbind(b_total, b)
      v_total = cbind(v_total, test_result[[3]])
      if (i == 1) {
        variance[i] <- (norm(B %*% G %*% b, type = "2") / sqrt(hd_obj$nobs - 1))^2
      }
      else{
        G_pc = t(b_total)%*%G%*%b_total
        coef_pc = B %*% G %*% b_total %*% solve(G_pc)
        total_variance = sum(diag((coef_pc%*%t(b_total)) %*% G %*% t(coef_pc%*%t(b_total))))
        G_pc_pre = t(b_total[, -i])%*%G%*%b_total[, -i]
        coef_pc_pre = B %*% G %*% b_total[, -i] %*% solve(G_pc_pre)
        total_variance_previous = sum(diag((coef_pc_pre%*%t(b_total[, -i])) %*% G %*% t(coef_pc_pre%*%t(b_total[, -i]))))
        variance[i] = (total_variance - total_variance_previous) / (hd_obj$nobs - 1)
      }
    }
  } 
  #########sequential inputs of smoothing parameters###########
  else{
    smooth_tuning <- expand.grid(smooth_tuning)
    
    I_a <- D <- S_2 <- S_smooth <- S_2_inverse <- list()
    cat("Preprocessing...\n")
    n_iter1 <- dim(smooth_tuning)[1]
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = n_iter1, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
    #####Do the following computation in advance to save computational cost#####
    for (smooth_index in 1:dim(smooth_tuning)[1]) {
      setTxtProgressBar(pb, smooth_index)
      I_a[[smooth_index]] <- MHPCA:::I_alpha(hd_obj, smooth_tuning[smooth_index, ])
      D[[smooth_index]] <- I_a[[smooth_index]] %*% smooth_penalty
      S_2[[smooth_index]] <- solve(G + D[[smooth_index]])
      S_2_inverse[[smooth_index]] = solve(S_2[[smooth_index]])
      S_smooth[[smooth_index]] <- G_half %*% (S_2[[smooth_index]]) %*% G_half
    }
    close(pb)
    b_total = v_total = c()
    GCV_score = list()
    if(sparse_CV == FALSE){
      CV_score = c()
    } else{
      CV_score = list()
    }
    for (i in 1:n) {
      cat(sprintf("Computing the %s PC...\n", ordinal_msg(i)))
      if (i == 1) {
        B_temp = B
      } 
      else{
        b_original = t(B_temp)%*%v
        B_temp = B_temp - v%*%t(b_original)
      }
      
      if(sparse_CV == FALSE){
        sparse_tuning_temp = sparse_tuning[i]
      } else{
        sparse_tuning_temp = sparse_tuning
      }
      
      cv_result = cv_gcv_sequential(data = B_temp, hd_obj = hd_obj, smooth_tuning = smooth_tuning, sparse_tuning = sparse_tuning_temp, sparse_tuning_type = sparse_tuning_type, K_fold = K_fold, G = G, G_half = G_half, G_half_inverse = G_half_inverse, S_smooth = S_smooth, S_2_inverse = S_2_inverse)
      sparse_result = cv_result[[1]]
      smooth_result_index = cv_result[[3]]
      if (sparse_CV == FALSE) {
        CV_score = c(CV_score, cv_result[[4]])
      } else{
        CV_score[[i]] = cv_result[[4]]
      }
      GCV_score[[i]] = cv_result[[5]]
      test_result = init_sequential(B_temp%*%G_half, sparse_result, sparse_tuning_type, S_smooth[[smooth_result_index]], S_2_inverse[[smooth_result_index]], G_half_inverse, G_half)
      
      v = test_result[[2]]
      b = test_result[[1]]
      smooth_tuning_result[[i]] = smooth_tuning[smooth_result_index, ]
      sparse_tuning_result[[i]] = sparse_result
      temp_count <- 0
      for (j in 1:p) {
        index_start <- (temp_count + 1)
        index_end <- (temp_count + prod(hd_obj$basis$nbasis[[j]]))
        if (i == 1) {
          pc[[j]] <- b[index_start:index_end, ]
        } else {
          pc[[j]] <- cbind(pc[[j]], b[index_start:index_end, ])
        }
        temp_count <- index_end
      }
      lsv = cbind(lsv, v)
      b_total = cbind(b_total, b)
      v_total = cbind(v_total, test_result[[3]])
      if (i == 1) {
        variance[i] <- (norm(B %*% G %*% b, type = "2") / sqrt(hd_obj$nobs - 1))^2
      }
      else{
        G_pc = t(b_total)%*%G%*%b_total
        coef_pc = B %*% G %*% b_total %*% solve(G_pc)
        total_variance = sum(diag((coef_pc%*%t(b_total)) %*% G %*% t(coef_pc%*%t(b_total))))
        G_pc_pre = t(b_total[, -i])%*%G%*%b_total[, -i]
        coef_pc_pre = B %*% G %*% b_total[, -i] %*% solve(G_pc_pre)
        total_variance_previous = sum(diag((coef_pc_pre%*%t(b_total[, -i])) %*% G %*% t(coef_pc_pre%*%t(b_total[, -i]))))
        variance[i] = (total_variance - total_variance_previous) / (hd_obj$nobs - 1)
      }
      
    }
  }
  return(list(pc, lsv, variance, smooth_tuning_result, sparse_tuning_result, CV_score, GCV_score))
} 


# joint smooth and sparse power algorithm
joint_power <- function(hd_obj, n, smooth_tuning, smooth_tuning_type, centerfns, alpha_orth) {
  p <- hd_obj$nvar
  smooth_penalty <- MHPCA:::pen_fun(hd_obj, type = smooth_tuning_type)
  
  #######centralize########
  B <- centralized_hd(hd_obj, centerfns)
  #########################
  
  ########some initial setting#######
  B <- t(B)
  lsv <- c()
  pc <- list()
  variance <- vector()
  smooth_tuning_result  <- list()
  sparse_tuning_result  <- list()
  G <- as.matrix(hd_obj$basis$gram)
  G_half <- sqrtm(G)
  G_half_inverse = solve(G_half)
  ###################################
  
  #####smoothing parameter#######
  if (is.null(smooth_tuning)) {
    for (i in 1:p) {
      smooth_tuning <- c(smooth_tuning, list(2^seq(-20, 20, length.out = 10)))
    }
  }
  smooth_tuning <- expand.grid(smooth_tuning)
  ####################################
  
  ####joint power####
  I_a <- D <- S_2 <- S_smooth <- S_2_inverse <- list()
  cat("Preprocessing...\n")
  n_iter1 <- dim(smooth_tuning)[1]
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter1, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  #####Do the following computation in advance to save computational cost#####
  for (smooth_index in 1:dim(smooth_tuning)[1]) {
    setTxtProgressBar(pb, smooth_index)
    I_a[[smooth_index]] <- MHPCA:::I_alpha(hd_obj, smooth_tuning[smooth_index, ])
    D[[smooth_index]] <- I_a[[smooth_index]] %*% smooth_penalty
    S_2[[smooth_index]] <- solve(G + D[[smooth_index]])
    S_2_inverse[[smooth_index]] = solve(S_2[[smooth_index]])
    S_smooth[[smooth_index]] <- G_half %*% (S_2[[smooth_index]]) %*% G_half
  }
  close(pb)
  cat(sprintf("Computing PCs...\n"))
  cv_result = gcv_joint(data = B, hd_obj = hd_obj, smooth_tuning = smooth_tuning, G = G, G_half = G_half, G_half_inverse = G_half_inverse, S_smooth = S_smooth, S_2_inverse = S_2_inverse, n = n)
  smooth_result_index = cv_result[[2]]
  GCV_score = cv_result[[3]]
  test_result = init_joint(B%*%G_half, S_smooth[[smooth_result_index]], S_2_inverse[[smooth_result_index]], G_half_inverse, G_half, n = n)
  v = test_result[[2]]
  b = test_result[[1]]
  smooth_tuning_result = smooth_tuning[smooth_result_index, ]
  temp_count <- 0
  
  temp_count <- 0
  for (j in 1:p) {
    index_start <- (temp_count + 1)
    index_end <- (temp_count + prod(hd_obj$basis$nbasis[[j]]))
    pc[[j]] <- b[index_start:index_end, ]
    temp_count <- index_end
  }
  
  lsv = cbind(lsv, v)
  b_total = b
  v_total = test_result[[3]]
  for (k in 1:n) {
    if (k == 1) {
      variance[k] <- (norm(B %*% G %*% b[, k], type = "2") / sqrt(hd_obj$nobs - 1))^2
    } else{
      G_pc = t(b_total)%*%G%*%b_total
      coef_pc = B %*% G %*% b_total %*% solve(G_pc)
      total_variance = sum(diag((coef_pc%*%t(b_total)) %*% G %*% t(coef_pc%*%t(b_total))))
      G_pc_pre = t(b_total[, 1:(k-1)])%*%G%*%b_total[, 1:(k-1)]
      coef_pc_pre = B %*% G %*% b_total[, 1:(k-1)] %*% solve(G_pc_pre)
      total_variance_previous = sum(diag((coef_pc_pre%*%t(b_total[, 1:(k-1)])) %*% G %*% t(coef_pc_pre%*%t(b_total[, 1:(k-1)]))))
      variance[k] = (total_variance - total_variance_previous) / (hd_obj$nobs - 1)
    }
  }
  return(list(pc, lsv, variance, smooth_tuning_result, GCV_score))
} 