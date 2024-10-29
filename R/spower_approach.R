#' @importFrom expm sqrtm
#' @importFrom utils  txtProgressBar setTxtProgressBar
csparse_pen_fun <- function(y,tuning_parameter = 0, type,alpha = 3.7) {
  #sparsity penalty 
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

vsparse_pen_fun <- function(y,tuning_parameter = 0, type,alpha = 3.7) {
  #sparsity penalty 
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

power_algo = function(data,
                      sparse_tuning_result,
                      sparse_tuning_type,
                      S_2 = NULL,
                      G_half_inverse = NULL,
                      G_half = NULL,
                      type = "real",
                      vdata = NULL,
                      vsparse_tuning_result,
                      vsparse_tuning_type){
  # data = t(data)
  if (is.null(vdata)) {
  b_old = svd(data)$v[,1]
  errors = 10^60
  while (errors > 10^(-7)) {
    v_new = csparse_pen_fun(y = data%*%b_old,tuning_parameter = sparse_tuning_result,sparse_tuning_type) 
    #v_new <- sign(v_new)*v_new ### identifieability

    if (type == "CV") {
      b_new = t(data)%*%v_new
    } else{
      b_new = G_half%*%S_2%*%G_half%*%t(data)%*%v_new
    }
    b_new = b_new / norm(b_new,type = "2")
    errors = sum((b_new - b_old)^2)
    b_old = b_new
    #error <- sqrt(sum(csparse_pen_fun(y = data%*%b_old,tuning_parameter = sparse_tuning_result,sparse_tuning_type) - v_new)^2)
    
  }
  v_new = v_new/norm(v_new,type = "2")
  if (type == "CV") {
    return(v_new)
  }
  else{
    b_new = G_half_inverse%*%b_new
    b_new = b_new %*% solve(sqrt(t(b_new) %*% solve(S_2) %*% b_new))
    return(list(b_new,v_new))
  }
  } else {
    bv_old = svd(vdata)$v[,1]
    b_old = svd(data)$v[,1]
    errors = 10^60
    while (errors > 10^(-7)) {
      v_new = csparse_pen_fun(y = data%*%b_old,tuning_parameter = sparse_tuning_result,sparse_tuning_type) + vdata%*%bv_old
      #v_new <- sign(v_new)*v_new
      if (type == "CV") {
        b_new = t(data)%*%v_new
        bv_new = t(vdata)%*%v_new 
        # add sparsity to vdata v's
        #bv_new <- vsparse_pen_fun(y = t(vdata)%*%v_new,tuning_parameter = vsparse_tuning_result,vsparse_tuning_type)
      } else{
        b_new = G_half%*%S_2%*%G_half%*%t(data)%*%v_new
        bv_new = t(vdata)%*%v_new
        # add sparsity to vdata v's
        #bv_new <- vsparse_pen_fun(y = t(vdata)%*%v_new,tuning_parameter = vsparse_tuning_result,vsparse_tuning_type)
      }
      b_new = b_new / norm(b_new,type = "2")
      bv_new <- bv_new / norm(bv_new,type = "2")
      errors = sqrt(sum((b_new - b_old)^2) + sum((bv_new - bv_old)^2 ))
      b_old = b_new
      bv_old <- bv_new
      #error <- sqrt(sum(csparse_pen_fun(y = data%*%b_old,tuning_parameter = sparse_tuning_result,sparse_tuning_type) + vdata%*%bv_old - v_new)^2)
    }
    v_new = v_new/norm(v_new,type = "2")
    if (type == "CV") {
      return(v_new)
    }
    else{
      b_new = G_half_inverse%*%b_new
      b_new = b_new %*% solve(sqrt(t(b_new) %*% solve(S_2) %*% b_new))
      return(list(b_new,v_new,bv_new))
    }
    
  }
}

cv_score_sparse = function(data,G,S,K_fold,sparse_tuning_single,sparse_tuning_type,shuffled_row,group_size,S_back,S_back_inverse,vdata = NULL){
  data_double_tilde = t(data%*%G%*%S)
  error_score_sparse = 0
  for (k in 1:K_fold) {
    rows_to_remove = shuffled_row[((k-1)*group_size+1):((k)*group_size)]
    data_train = data_double_tilde[-rows_to_remove,]
    data_test = data_double_tilde[rows_to_remove,]
    v_test = power_algo(t(data_train),sparse_tuning_single,sparse_tuning_type,type = "CV",vdata = vdata)
    b_test = data_test%*%v_test
    b_test_smooth_back = S_back[rows_to_remove,]%*%(data_double_tilde%*%v_test)
    data_test_smooth_back = t(data_double_tilde)%*%S_back_inverse[,rows_to_remove]
    error_score_sparse = error_score_sparse + sum((t(data_test_smooth_back)-b_test_smooth_back%*%t(v_test))^2)
  }
  return(error_score_sparse/ncol(data))
}


cv_score_smooth = function(data,mvmfd_obj,G,G_half,S_smooth,v,smooth_tuning){
  p = mvmfd_obj$nvar
  B_subtilde = data%*%G_half
  if (all(smooth_tuning == 0)) {
    error_smooth_score <- 0
  } else {
    if (p == 1) {
      error_smooth_score <- (sum(((diag(dim(S_smooth)[1]) - S_smooth) %*% (t(B_subtilde) %*% v))^2) / ((1 - sum(diag(S_smooth)) / dim(G)[1])^2)) / dim(G)[1]
    }
    else{
      index1 = prod(mvmfd_obj$basis$nbasis[[1]])
      index2 = prod(mvmfd_obj$basis$nbasis[[2]])
      s_alpha_tilde_1 = S_smooth[1:index1,1:index1]
      s_alpha_tilde_2 = S_smooth[(1+index1):(index1+index2),(1+index1):(index1+index2)]
      B_subtilde_1 = B_subtilde[,1:index1]
      B_subtilde_2 = B_subtilde[,(1+index1):(index1+index2)]
      error_smooth_score_1 <- sum(((diag(index1) - s_alpha_tilde_1) %*% (t(B_subtilde_1) %*% v))^2)/(1-sum(diag(s_alpha_tilde_1))/index1)^2
      error_smooth_score_2 <- sum(((diag(index2) - s_alpha_tilde_2) %*% (t(B_subtilde_2) %*% v))^2)/(1-sum(diag(s_alpha_tilde_2))/index2)^2
      error_smooth_score = error_smooth_score_1 + error_smooth_score_2
    }
  }
  return(error_smooth_score)
}

# ignore this
cv_selection_joint = function(data,mvmfd_obj,smooth_tuning,sparse_tuning,sparse_tuning_type,K_fold,G,G_half,G_half_inverse,S,S_2,S_inverse,S_smooth,S_back,S_back_inverse,vdata = NULL){
  CV_score = 10^60
  result = c()
  count = 0
  shuffled_row = sample(ncol(data))
  group_size = length(shuffled_row) / K_fold
  n_iter <- dim(smooth_tuning)[1] * length(sparse_tuning)
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for (smooth_index in 1:dim(smooth_tuning)[1]) {
    if (all(smooth_tuning == 0)) {
      S[[smooth_index]] = solve(G_half)
      S_2[[smooth_index]] = solve(G)
      S_smooth[[smooth_index]] = diag(dim(G)[1])
      S_inverse[[smooth_index]] = G_half_inverse
      S_back_inverse[[smooth_index]] = diag(ncol(S_2[[smooth_index]]))
      S_back[[smooth_index]] = diag(ncol(S_2[[smooth_index]]))
    }
    for (sparse_tuning_single in sparse_tuning) {
      count = count + 1
      setTxtProgressBar(pb, count)
      if (sparse_tuning_single == 0) {
        sparse_score = 0
      } else{
        sparse_score = cv_score_sparse(data,G,S[[smooth_index]],K_fold,sparse_tuning_single,sparse_tuning_type,shuffled_row,group_size,S_back[[smooth_index]],S_back_inverse[[smooth_index]],vdata = vdata)
      }
      test_temp = power_algo(data%*%G_half,sparse_tuning_single,sparse_tuning_type,S_2[[smooth_index]],G_half_inverse,G_half,vdata)
      v_temp = test_temp[[2]]
      smooth_score = cv_score_smooth(data,mvmfd_obj,G,G_half,S_smooth[[smooth_index]],v_temp,smooth_tuning=smooth_tuning[smooth_index,])
      smooth_index_result = smooth_index
      sparse_index_result = sparse_tuning_single
      result = rbind(result, c(sparse_score,smooth_score,sparse_index_result,smooth_index_result))
      # print(c(sparse_score,smooth_score,sparse_index_result,smooth_index_result))
    }
  }
  close(pb) # Close the connection
  return(result)
}


cv_selection_marginal = function(data,mvmfd_obj,smooth_tuning,sparse_tuning,sparse_tuning_type,K_fold,G,G_half,G_half_inverse,S_2,S_smooth,vdata = NULL){
  # tuning parameter selection for sequential power
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
  
  for (sparse_tuning_single in sparse_tuning) {
    count = count +1
    setTxtProgressBar(pb, count)
    if (sparse_tuning_single == 0) {
      sparse_score = 0
    } else{
      sparse_score = cv_score_sparse(data,G,G_half_inverse,K_fold,sparse_tuning_single,sparse_tuning_type,shuffled_row,group_size,diag(ncol(data)),diag(ncol(data)),vdata = vdata)
    }
    if (sparse_score <= CV_score_sparse) {
      CV_score_sparse = sparse_score
      sparse_tuning_selection = sparse_tuning_single
    }
  }
  for (smooth_index in 1:dim(smooth_tuning)[1]) {
    count = count + 1
    setTxtProgressBar(pb, count)
    if (all(smooth_tuning == 0)) {
      S_2[[smooth_index]] = solve(G)
      S_smooth[[smooth_index]] = diag(dim(G)[1])
    }
    test_temp = power_algo(data%*%G_half,sparse_tuning_selection,sparse_tuning_type,S_2[[smooth_index]],G_half_inverse,G_half,vdata = vdata)
    v_temp = test_temp[[2]]
    smooth_score = cv_score_smooth(data,mvmfd_obj,G,G_half,S_smooth[[smooth_index]],v_temp,smooth_tuning=smooth_tuning[smooth_index,])
    if (smooth_score <= CV_score_smooth) {
      CV_score_smooth = smooth_score
      smooth_tuning_selection = smooth_tuning[smooth_index,]
      index_selection = smooth_index
    }
  }
  close(pb) # Close the connection
  result = list(sparse_tuning_selection,smooth_tuning_selection,index_selection)
  return(result)
}

ordinal <- function(i) {
  # printing 
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

# main function 
s_power_algorithm <- function(mvmfd_obj, 
                              n, 
                              smooth_tuning, 
                              smooth_tuning_type,
                              sparse_tuning, 
                              sparse_tuning_type,
                              centerfns, 
                              alpha_orth,
                              K_fold = 10,
                              cv_type = "marginal",
                              vdata = NULL) {
  #browser()
  p <- mvmfd_obj$nvar
  smooth_penalty <- ReMFPCA:::pen_fun(mvmfd_obj, type = smooth_tuning_type)
  B <- c()
  
  #######smoothing parameter##########
  if (is.null(smooth_tuning)) {
    for (i in 1:p) {
      smooth_tuning <- c(smooth_tuning, list(2^seq(-20, 20, length.out = 10)))
    }
  }
  smooth_tuning <- expand.grid(smooth_tuning)
  ####################################
  
  #######centralize########
  if (centerfns == TRUE) {
    for (i in 1:p) {
      if (is.matrix(mvmfd_obj$coefs[[i]])) {
        c <- mvmfd_obj$coefs[[i]] - rowMeans(mvmfd_obj$coefs[[i]])
        B <- rbind(B, c)
      } else {
        cc <- apply(mvmfd_obj$coefs[[i]], 3, as.vector)
        c <- cc - rowMeans(cc)
        B <- rbind(B, c)
      }
    }
  } else {
    for (i in 1:p) {
      if (is.matrix(mvmfd_obj$coefs[[i]])) {
        c <- mvmfd_obj$coefs[[i]]
        B <- rbind(B, c)
      } else {
        cc <- apply(mvmfd_obj$coefs[[i]], 3, as.vector)
        B <- rbind(B, cc)
      }
    }
  }
  #########################
  
  
  ########some initial setting#######
  B <- t(B)
  lsv <- c()
  pc <- list()
  variance <- vector()
  smooth_tuning_result  <- list()
  sparse_tuning_result  <- list()
  G <- as.matrix(mvmfd_obj$basis$gram)
  G_half <- sqrtm(G)
  G_half_inverse = solve(G_half)
  
  ###################################
  ####sequential power####
  if (alpha_orth == FALSE) {
    I_a <- D <- S_2 <- S <- S_smooth <- S_inverse <- S_back_inverse <- S_back <- list()
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
      I_a[[smooth_index]] <- ReMFPCA:::I_alpha(mvmfd_obj, smooth_tuning[smooth_index, ])
      D[[smooth_index]] <- I_a[[smooth_index]] %*% smooth_penalty
      S_2[[smooth_index]] <- solve(G + D[[smooth_index]])
      S[[smooth_index]] = sqrtm(S_2[[smooth_index]])
      S_inverse[[smooth_index]] = solve(S_2[[smooth_index]])
      S_smooth[[smooth_index]] <- G_half %*% (S_2[[smooth_index]]) %*% G_half
      S_back[[smooth_index]] = G_half%*%S[[smooth_index]]
      S_back_inverse[[smooth_index]] = solve(S_back[[smooth_index]])
    }
    close(pb)
    if (!is.null(vdata)){
      bv <- NULL
    }
    for (i in 1:n) {
      cat(sprintf("Computing the %s PC...\n", ordinal(i)))
      if (i == 1) {
        B_temp = B
        if (!is.null(vdata)) {
          vdata_temp <- vdata
        } else {
          vdata_temp <- NULL
        }
      } 
      else{
        b_original = t(B_temp)%*%v
        B_temp = B_temp - v%*%t(b_original)
        if (!is.null(vdata)) {
          bv_original = t(vdata)%*%v
          vdata_temp <- vdata_temp - v%*%t(bv_original)
        } else {
          vdata_temp <- NULL
        }
        
      }
      if (cv_type == "joint") {
        cv_result = cv_selection_joint(data = B_temp, mvmfd_obj = mvmfd_obj, smooth_tuning = smooth_tuning, sparse_tuning = sparse_tuning, sparse_tuning_type = sparse_tuning_type, K_fold = K_fold, S = S, S_inverse = S_inverse,G = G,G_half = G_half,S_2 = S_2,S_smooth = S_smooth,G_half_inverse = G_half_inverse,S_back = S_back, S_back_inverse = S_back_inverse,vdata = vdata_temp)
        sparse_result = cv_result[which.min(cv_result[,1]+cv_result[,2]),3]
        smooth_result_index = cv_result[which.min(cv_result[,1]+cv_result[,2]),4]
        test_result = power_algo(B_temp%*%G_half,sparse_result,sparse_tuning_type,S_2[[smooth_result_index]],G_half_inverse,G_half,vdata = vdata_temp)
      }
      else if (cv_type == "marginal") {
        cv_result = cv_selection_marginal(data = B_temp, mvmfd_obj = mvmfd_obj, smooth_tuning = smooth_tuning, sparse_tuning = sparse_tuning, sparse_tuning_type = sparse_tuning_type, K_fold = K_fold,G = G,G_half = G_half,G_half_inverse = G_half_inverse,S_2 = S_2,S_smooth = S_smooth,vdata = vdata_temp )
        sparse_result = cv_result[[1]]
        smooth_result_index = cv_result[[3]]
        test_result = power_algo(B_temp%*%G_half,sparse_result,sparse_tuning_type,S_2[[smooth_result_index]],G_half_inverse,G_half,vdata = vdata_temp)
      }
      
      # sparse_result = cv_result[which.min(cv_result[,1]+cv_result[,2]),3]
      # smooth_result_index = cv_result[which.min(cv_result[,1]+cv_result[,2]),4]
      # print(cbind(scale_sparse,scale_smooth,cv_result[,3:4]))
      # print(sparse_result)
      # print(smooth_result_index)
      # test_result = power_algo_v2(B_temp%*%G_half,sparse_result,sparse_tuning_type,S[[smooth_result_index]],G,S_2[[smooth_result_index]],S_inverse[[smooth_result_index]],G_half_inverse,G_half)
      v = test_result[[2]]
      b = test_result[[1]]
      if (!is.null(vdata)){
        bv <- cbind(bv,test_result[[3]])
      } else {
        bv <- NULL
      }
      smooth_tuning_result[[i]] = smooth_tuning[smooth_result_index,]
      sparse_tuning_result[[i]] = sparse_result
      temp_count <- 0
      for (j in 1:p) {
        index_start <- (temp_count + 1)
        index_end <- (temp_count + prod(mvmfd_obj$basis$nbasis[[j]]))
        if (i == 1) {
          pc[[j]] <- b[index_start:index_end, ]
        } else {
          pc[[j]] <- cbind(pc[[j]], b[index_start:index_end, ])
        }
        temp_count <- index_end
      }
      lsv = cbind(lsv,v)
      BGb <- B %*% G %*% b
      if (!is.null(vdata)){
        Xv <- vdata %*% bv[,i]
        variance[i] <- (t(Xv)%*%Xv+2*t(Xv)%*%BGb+t(BGb)%*%BGb)/(mvmfd_obj$nobs - 1)
      } else {
        variance[i] <- (t(BGb)%*%BGb)/(mvmfd_obj$nobs - 1)
      }
    }
  }
  return(list(pc, lsv, variance,smooth_tuning_result,sparse_tuning_result,bv))
} 
