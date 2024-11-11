#' @importFrom expm sqrtm
#' @importFrom utils  txtProgressBar setTxtProgressBar
eigen_approach <- function(hd_obj, n, alpha, centerfns, penalty_type) {
  m.rep <- hd_obj$nobs
  p <- hd_obj$nvar
  if (is.null(alpha)) {
    for (i in 1:p) {
      alpha <- c(alpha, list(2^seq(-20, 20, length.out = 5)))
    }
  }
  if (p == 2) {
    gcv_row <- length(alpha[[1]])
    gcv_column <- length(alpha[[2]])
    index1 = hd_obj$basis$nbasis[[1]]
    index2 = hd_obj$basis$nbasis[[2]]
  }
  alpha <- expand.grid(alpha) 
  penalty <- pen_fun(hd_obj, type = penalty_type)
  G <- as.matrix(hd_obj$basis$gram)
  G_half <- expm::sqrtm(G)

  B <- c()
  B_c <- c()
  if (centerfns) {
    for (i in 1:p) {
      if (is.matrix(hd_obj$coefs[[i]])) {
        B <- rbind(B, hd_obj$coefs[[i]])
        B_c <- rbind(B_c, hd_obj$coefs[[i]] - rowMeans(hd_obj$coefs[[i]]))
      } else {
        cc <- apply(hd_obj$coefs[[i]], 3, as.vector)
        B <- rbind(B, cc)
        B_c <- rbind(B_c, cc - rowMeans(cc))
      }
    }
  } else {
    for (i in 1:p) {
      if (is.matrix(hd_obj$coefs[[i]])) {
        B <- rbind(B, hd_obj$coefs[[i]])
      } else {
        cc <- apply(hd_obj$coefs[[i]], 3, as.vector)
        B <- rbind(B, cc)
      }
    }
    B_c <- B
  }
  B <- t(B)
  B_c <- t(B_c)
  B_subtilde <- B_c %*% G_half
  I_matrix <- diag(1, m.rep)
  J <- matrix(1, m.rep, m.rep)
  if (centerfns) {
    V <- (1 / (m.rep - 1)) * (t(B) %*% (I_matrix - (1 / m.rep) * J) %*% B)
  } else {
    V <- (1 / (m.rep - 1)) * (t(B) %*% B)
  }
  GCV_score <- 10^60
  GCVs <- c()
  # Initializes the progress bar
  n_iter <- dim(alpha)[1]
  pb <- txtProgressBar(
    min = 0, # Minimum value of the progress bar
    max = n_iter, # Maximum value of the progress bar
    style = 3, # Progress bar style (also available style = 1 and style = 2)
    width = 50, # Progress bar width. Defaults to getOption("width")
    char = "="
  ) # Character used to create the bar

  for (j in 1:dim(alpha)[1]) {
    setTxtProgressBar(pb, j)
    I <- I_alpha(hd_obj, alpha[j, ])
    D <- I %*% penalty
    L <- t(chol(G + D))
    S <- as.matrix(solve(L))
    E <- eigen(S %*% t(G) %*% V %*% G %*% t(S))
    u <- E$vectors
    s_alpha <- sqrtm(solve(G + D))
    s_alpha_tilde <- G_half %*% (s_alpha %*% s_alpha) %*% G_half
    b_temp <- c()
    for (k in 1:n) {
      b_temp <- cbind(b_temp, ((t(S) %*% u[, k]) %*% (t(u[, k]) %*% S %*% G %*% t(S) %*% u[, k])^(-0.5)))
    }
    v_temp <- B_c %*% G %*% b_temp

    if (all(alpha[j, ] == 0)) {
      GCV_score_temp <- 0
    } else {
      # GCV_score_temp <- (sum(((diag(dim(s_alpha_tilde)[1]) - s_alpha_tilde) %*% (t(B_subtilde) %*% v_temp))^2) / ((1 - sum(diag(s_alpha_tilde)) / dim(G)[1])^2)) / dim(G)[1]
      if (p == 1) {
        GCV_score_temp <- (sum(((diag(dim(s_alpha_tilde)[1]) - s_alpha_tilde) %*% (t(B_subtilde) %*% v_temp))^2) / ((1 - sum(diag(s_alpha_tilde)) / dim(G)[1])^2)) / dim(G)[1]
      }
      else{
        s_alpha_tilde_1 = s_alpha_tilde[1:index1,1:index1]
        s_alpha_tilde_2 = s_alpha_tilde[(1+index1):(index1+index2),(1+index1):(index1+index2)]
        B_subtilde_1 = B_subtilde[,1:index1]
        B_subtilde_2 = B_subtilde[,(1+index1):(index1+index2)]
        GCV_score_temp_1 <- sum(((diag(index1) - s_alpha_tilde_1) %*% (t(B_subtilde_1) %*% v_temp))^2)/(1-sum(diag(s_alpha_tilde_1))/index1)^2
        GCV_score_temp_2 <- sum(((diag(index2) - s_alpha_tilde_2) %*% (t(B_subtilde_2) %*% v_temp))^2)/(1-sum(diag(s_alpha_tilde_2))/index2)^2
        GCV_score_temp = GCV_score_temp_1 + GCV_score_temp_2
      }
    }
    GCVs <- c(GCVs, GCV_score_temp)

    if (GCV_score_temp < GCV_score) {
      b <- b_temp
      v <- v_temp
      GCV_score <- GCV_score_temp
      GCV_result <- alpha[j, ]
    }
  }
  close(pb) # Close the connection
  if (p == 2) {
    GCVs <- matrix(GCVs, nrow = gcv_row, ncol = gcv_column)
  }
  temp_count <- 0
  pc <- list()
  for (i in 1:p) {
    index_start <- (temp_count + 1)
    index_end <- (temp_count + prod(hd_obj$basis$nbasis[[i]]))
    pc[[i]] <- b[index_start:index_end, ]
    temp_count <- temp_count + prod(hd_obj$basis$nbasis[[i]])
  }
  variance <- diag(t(b) %*% G %*% V %*% G %*% b)
  sigma <- sqrt(variance)
  lsv <- (B_c %*% G %*% b) %*% solve(diag(sqrt(diag(t(B_c %*% G %*% b) %*% (B_c %*% G %*% b)))))

  bbbb <- c()
  for (k in 1:p) {
    bbbb <- rbind(bbbb, pc[[k]])
  }
  return(list(pc, lsv, variance, GCV_result, 0))
}
