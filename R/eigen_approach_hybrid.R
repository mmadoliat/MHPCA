#' @importFrom expm sqrtm
#' @importFrom utils  txtProgressBar setTxtProgressBar

eigen_approach_hybrid <- function(hd_obj, n, alpha, centerfns, penalty_type) {
  mvmfd_obj <- hd_obj$mf
  mvnfd_obj <- hd_obj$nf
  if (!is.null(mvmfd_obj)){
    m.rep <- mvmfd_obj$nobs
    p <- mvmfd_obj$nvar
    indices <- sapply(1:p, function(i) prod(mvmfd_obj$basis$nbasis[[i]]))
    if (is.null(alpha)) {
      for (i in 1:p) {
        alpha <- c(alpha, list(2^seq(-20, 20, length.out = 5)))
      }
    }
    if (p == 2) {
      gcv_row <- length(alpha[[1]])
      gcv_column <- length(alpha[[2]])
      index1 = mvmfd_obj$basis$nbasis[[1]]
      index2 = mvmfd_obj$basis$nbasis[[2]]
    }
    alpha <- expand.grid(alpha) 
    penalty <- pen_fun(mvmfd_obj, type = penalty_type)
    G <- as.matrix(mvmfd_obj$basis$gram)
    G_half <- expm::sqrtm(G)
    
    B <- c()
    B_c <- c()
    
    if (centerfns) {
      for (i in 1:p) {
        if (is.matrix(mvmfd_obj$coefs[[i]])) {
          B <- rbind(B, mvmfd_obj$coefs[[i]])
          B_c <- rbind(B_c, mvmfd_obj$coefs[[i]] - rowMeans(mvmfd_obj$coefs[[i]]))
        } else {
          cc <- apply(mvmfd_obj$coefs[[i]], 3, as.vector)
          B <- rbind(B, cc)
          B_c <- rbind(B_c, cc - rowMeans(cc))
        }
      }
     if (!is.null(mvnfd_obj)) { 
       mvnfd_obj <- center_mvnfd(mvnfd_obj)
       mvnfd_data <- do.call("cbind",mvnfd_obj$data)
     } else {
         mvnfd_obj <- NULL
         mvnfd_data <- NULL
       }
    } else {
      for (i in 1:p) {
        if (is.matrix(mvmfd_obj$coefs[[i]])) {
          B <- rbind(B, mvmfd_obj$coefs[[i]])
        } else {
          cc <- apply(mvmfd_obj$coefs[[i]], 3, as.vector)
          B <- rbind(B, cc)
        }
      }
      B_c <- B
      mvnfd_data <- if(!is.null(mvnfd_obj)) do.call("cbind",mvnfd_obj$data) else NULL
    }
    B <- t(B)
    B_c <- t(B_c)
    
  } else {
    m.rep <- mvnfd_obj$nobs
    p <- indices <-  penalty <-  G_half <- B <- B_c <- NULL
    alpha <- expand.grid(1)
    G <- matrix(nrow=0,ncol=0)
    mvnfd_obj <- if (centerfns) center_mvnfd(mvnfd_obj) 
    mvnfd_data <- do.call("cbind",mvnfd_obj$data)
    }
  # B_subtilde <- B_c %*% G_half
  # I_matrix <- diag(1, m.rep)
  # J <- matrix(1, m.rep, m.rep)
  # if (centerfns) {
  #   V <- (1 / (m.rep - 1)) * (t(B) %*% (I_matrix - (1 / m.rep) * J) %*% B)
  # } else {
  #   V <- (1 / (m.rep - 1)) * (t(B) %*% B)
  # }
  BG <- if (!is.null(mvmfd_obj)) B_c%*%G else NULL
  Z <- cbind(BG,mvnfd_data); tmp <- t(Z)%*%Z
  ZtZ <- t(Z)%*%Z
  V <- (1 / (m.rep - 1)) * ZtZ
  GCV_score <- 10^60
  GCVs <- c()
  # Initializes the progress bar
  n_iter <- if (!is.null(alpha)) dim(alpha)[1] else 1
  pb <- txtProgressBar(
    min = 0, # Minimum value of the progress bar
    max = n_iter, # Maximum value of the progress bar
    style = 3, # Progress bar style (also available style = 1 and style = 2)
    width = 50, # Progress bar width. Defaults to getOption("width")
    char = "="
  ) # Character used to create the bar
  
  for (j in 1:dim(alpha)[1]) {
    setTxtProgressBar(pb, j)
    I <- if (!is.null(mvmfd_obj)) I_alpha(mvmfd_obj, alpha[j, ]) else NULL
    D <- if (!is.null(mvmfd_obj)) I %*% penalty else NULL
    if (!is.null(mvmfd_obj) && !is.null(mvnfd_obj)){
      L <- as.matrix(Matrix::bdiag(t(chol(G + D)),diag(ncol(mvnfd_data))))
    } else if (!is.null(mvmfd_obj) && is.null(mvnfd_obj)) {
      L <- as.matrix(t(chol(G + D)))
    } else if (is.null(mvmfd_obj) && !is.null(mvnfd_obj)){
      L <- diag(ncol(mvnfd_data))
    }
   
    #L <- t(chol(G + D))
    S <-  as.matrix(solve(L))
    E <- eigen(S%*%V%*%t(S))
    #E <- eigen(S %*% t(G) %*% V %*% G %*% t(S))
    u <- E$vectors
    s_alpha <- if (!is.null(mvmfd_obj)) sqrtm(solve(G + D)) else NULL
    s_alpha_tilde <- if (!is.null(mvmfd_obj)) G_half %*% (solve(G + D)) %*% G_half else NULL
    b_temp <- c()
    for (k in 1:n) {
      b_temp <- cbind(b_temp, ((t(S) %*% u[, k]) %*% (t(u[, k]) %*% S %*% Matrix::bdiag(G,diag(ncol(mvnfd_data))) %*% t(S) %*% u[, k])^(-0.5)))
    }
    if (!is.null(mvmfd_obj) && !is.null(mvnfd_obj)){
      bv_temp <- as.matrix(b_temp[(nrow(b_temp) - (ncol(mvnfd_data)- 1)):nrow(b_temp), ])
      b_temp <- as.matrix(b_temp[1:(nrow(b_temp) - ncol(mvnfd_data)),])
      v_temp <- B_c %*% G %*% b_temp + mvnfd_data%*%bv_temp
    } else if (!is.null(mvmfd_obj) && is.null(mvnfd_obj)){
      b_temp <- as.matrix(b_temp)
      bv_temp <- NULL
      v_temp <- B_c %*% G %*% b_temp 
    } else if (is.null(mvmfd_obj) && !is.null(mvnfd_obj)){
      bv_temp <- as.matrix(b_temp)
      b_temp <- NULL
      v_temp <- mvnfd_data%*%bv_temp
    }

    v_temp <- as.matrix(v_temp)
    v_temp <- sweep(v_temp,2,sqrt(diag(t(v_temp)%*%v_temp)),"/")
    GCV_score_temp = if (!is.null(mvmfd_obj)) gcv_local(data = B_c, mvmfd_obj = mvmfd_obj, G = G, G_half = G_half, S_smooth = s_alpha_tilde, u = v_temp, smooth_tuning = alpha[j, ]) else NULL
    GCVs <- c(GCVs, GCV_score_temp)
    if (!is.null(mvmfd_obj)){
      if (GCV_score_temp < GCV_score) {
        b <- b_temp
        v <- v_temp
        GCV_score <- GCV_score_temp
        GCV_result <- alpha[j, ]
      }
    } else {
      b <- GCV_score <- GCV_result <- NULL
      v <- v_temp
    }
    
  }
  close(pb) # Close the connection
  if (!is.null(mvmfd_obj)){
    if (p == 2) {
      GCVs <- matrix(GCVs, nrow = gcv_row, ncol = gcv_column)
    }
  }
  
  temp_count <- 0
  pc <- list()
  if (!is.null(mvmfd_obj)){
    for (i in 1:p) {
      index_start <- (temp_count + 1)
      index_end <- (temp_count + prod(mvmfd_obj$basis$nbasis[[i]]))
      pc[[i]] <- b[index_start:index_end, ]
      temp_count <- temp_count + prod(mvmfd_obj$basis$nbasis[[i]])
    }
  } else {
    pc <- NULL
  }
  
  bv <- bv_temp
  #variance <- diag(t(b) %*% G %*% V %*% G %*% b)
  variance <- (diag(t(as.matrix(rbind(b,bv))) %*% V %*% as.matrix(rbind(b,bv))))
  lsv <- scale(v, center = FALSE, scale = sqrt(colSums(as.matrix(v)^2)))

  return(list(pc_fd = pc, lsv = lsv, variance = variance, GCV_result = GCV_result, GCVs = GCVs,pc_nfd = bv))
}
