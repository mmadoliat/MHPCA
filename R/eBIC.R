#' @title Extended Bayesian Information Criterion for Hybrid PCA
#'
#' @description
#' Computes the extended BIC (eBIC) for sparse hybrid functional-vector PCA.
#' Uses simplified eBIC formula: eBIC = -2log(L) + q*log(n) + 2*xi*q*log(p)
#'
#' @param fdata Functional data coefficient matrix (n x J), transformed (C %*% G_half)
#' @param nfdata Non-functional data matrix (n x m)
#' @param u_hat Estimated left singular vector (n x 1)
#' @param fv_hat Estimated functional right singular vector (J x 1), or NULL
#' @param nfv_hat Estimated non-functional right singular vector (m x 1), or NULL
#' @param hd_obj The hybrid data object
#' @param G_half Square root of Gram matrix for functional data
#' @param xi eBIC tuning parameter in [0, 1]. Default is 0.5.
#' @param N_type How to compute effective sample size
#' @param RSS_type Character: "direct" for direct RSS (Scenario 1) or 
#'        "refitted" for least squares refitted RSS (Scenario 2). Default is "direct".
#' @param tuning_stage Character: which tuning stage we are in.
#'        "u" for tuning gamma_u, "nfd" for tuning gamma_nfd, "fd" for tuning gamma_fd,
#'        "all" for using all components (default, backward compatible).
#'
#' @return A list containing eBIC, BIC, RSS, degrees of freedom, etc.
#'
#' @export
compute_eBIC_hybrid <- function(fdata,
                                nfdata,
                                u_hat,
                                fv_hat = NULL,
                                nfv_hat = NULL,
                                hd_obj = NULL,
                                G_half = NULL,
                                xi = 0.5,
                                N_type = "product",
                                RSS_type = "direct",
                                tuning_stage = "all") {
  
  # Validate inputs
  RSS_type <- match.arg(RSS_type, c("direct", "refitted"))
  tuning_stage <- match.arg(tuning_stage, c("all", "u", "nfd", "fd"))
  
  n <- length(u_hat)
  u_hat_vec <- as.vector(u_hat)
  
  # Dimensions
  if (!is.null(fdata)) {
    J <- if (is.matrix(fdata)) ncol(fdata) else length(fdata)
  } else {
    J <- 0
  }
  
  if (!is.null(nfdata)) {
    m <- ncol(nfdata)
  } else {
    m <- 0
  }
  
  # Count non-zero elements for each component
  tol_zero <- .Machine$double.eps^0.5
  k_u <- sum(abs(u_hat) > tol_zero)
  k_fv <- if (!is.null(fv_hat)) sum(abs(fv_hat) > tol_zero) else 0
  k_nfv <- if (!is.null(nfv_hat)) sum(abs(nfv_hat) > tol_zero) else 0
  
  # Determine component-specific df and dimension p based on tuning_stage
  if (tuning_stage == "u") {
    df <- k_u
    p_dim <- n
  } else if (tuning_stage == "nfd") {
    df <- k_nfv
    p_dim <- m
  } else if (tuning_stage == "fd") {
    df <- k_fv
    p_dim <- J
  } else {
    # "all" - backward compatible
    df <- k_u + k_fv + k_nfv
    p_dim <- n + J + m
  }
  
  # Effective sample size for likelihood computation
  total_features <- J + m
  if (is.numeric(N_type) && length(N_type) == 1) {
    N <- N_type
  } else if (N_type == "product") {
    N <- n * total_features
  } else if (N_type == "sum") {
    N <- n + total_features
  } else {
    N <- n * total_features
  }
  
  # Compute RSS based on selected method and tuning_stage
  RSS <- 0
  
  if (RSS_type == "direct") {
    #-----------------------------------------
    # Scenario 1: Direct RSS
    # RSS = ||Z - u_hat * v_hat^T||_F^2
    #-----------------------------------------
    
    # For tuning_stage specific RSS
    if (tuning_stage == "fd" || tuning_stage == "all" || tuning_stage == "u") {
      if (!is.null(fdata) && !is.null(fv_hat)) {
        fv_hat_vec <- as.vector(fv_hat)
        fdata_reconstructed <- u_hat_vec %*% t(fv_hat_vec)
        if (tuning_stage == "fd") {
          RSS <- sum((fdata - fdata_reconstructed)^2)
        } else {
          RSS <- RSS + sum((fdata - fdata_reconstructed)^2)
        }
      }
    }
    
    if (tuning_stage == "nfd" || tuning_stage == "all" || tuning_stage == "u") {
      if (!is.null(nfdata) && !is.null(nfv_hat)) {
        nfv_hat_vec <- as.vector(nfv_hat)
        nfdata_reconstructed <- u_hat_vec %*% t(nfv_hat_vec)
        if (tuning_stage == "nfd") {
          RSS <- sum((nfdata - nfdata_reconstructed)^2)
        } else {
          RSS <- RSS + sum((nfdata - nfdata_reconstructed)^2)
        }
      }
    }
    
  } else {
    #-----------------------------------------
    # Scenario 2: Refitted with Scale Estimation
    # Fix u and v, estimate scale s that minimizes ||X - s * u * v^T||^2
    # s = (u^T X v) / (||u||^2 * ||v||^2)
    # RSS = ||X - s * u * v^T||^2
    #-----------------------------------------
    
    u_norm_sq <- sum(u_hat_vec^2)
    
    # Functional component
    if (tuning_stage == "fd" || tuning_stage == "all" || tuning_stage == "u") {
      if (!is.null(fdata) && !is.null(fv_hat)) {
        fv_hat_vec <- as.vector(fv_hat)
        fv_norm_sq <- sum(fv_hat_vec^2)
        
        if (u_norm_sq > .Machine$double.eps && fv_norm_sq > .Machine$double.eps) {
          # s_f = (u^T fdata v) / (||u||^2 * ||v||^2)
          # u^T fdata v = t(u) %*% fdata %*% v
          numerator_f <- as.numeric(t(u_hat_vec) %*% fdata %*% fv_hat_vec)
          denominator_f <- u_norm_sq * fv_norm_sq
          s_f <- numerator_f / denominator_f
          
          # Reconstruction with scale
          fdata_reconstructed <- s_f * (u_hat_vec %*% t(fv_hat_vec))
          
          if (tuning_stage == "fd") {
            RSS <- sum((fdata - fdata_reconstructed)^2)
          } else {
            RSS <- RSS + sum((fdata - fdata_reconstructed)^2)
          }
        } else {
          if (tuning_stage == "fd") {
            RSS <- sum(fdata^2)
          } else {
            RSS <- RSS + sum(fdata^2)
          }
        }
      }
    }
    
    # Non-functional component
    if (tuning_stage == "nfd" || tuning_stage == "all" || tuning_stage == "u") {
      if (!is.null(nfdata) && !is.null(nfv_hat)) {
        nfv_hat_vec <- as.vector(nfv_hat)
        nfv_norm_sq <- sum(nfv_hat_vec^2)
        
        if (u_norm_sq > .Machine$double.eps && nfv_norm_sq > .Machine$double.eps) {
          # s_nf = (u^T nfdata v) / (||u||^2 * ||v||^2)
          numerator_nf <- as.numeric(t(u_hat_vec) %*% nfdata %*% nfv_hat_vec)
          denominator_nf <- u_norm_sq * nfv_norm_sq
          s_nf <- numerator_nf / denominator_nf
          
          # Reconstruction with scale
          nfdata_reconstructed <- s_nf * (u_hat_vec %*% t(nfv_hat_vec))
          
          if (tuning_stage == "nfd") {
            RSS <- sum((nfdata - nfdata_reconstructed)^2)
          } else {
            RSS <- RSS + sum((nfdata - nfdata_reconstructed)^2)
          }
        } else {
          if (tuning_stage == "nfd") {
            RSS <- sum(nfdata^2)
          } else {
            RSS <- RSS + sum(nfdata^2)
          }
        }
      }
    }
  }
  
  # Avoid log(0) issues
  if (RSS < .Machine$double.eps) {
    RSS <- .Machine$double.eps
  }
  
  # Adjust N for component-specific tuning stages
  if (tuning_stage == "fd") {
    N_eff <- n * J
  } else if (tuning_stage == "nfd") {
    N_eff <- n * m
  } else {
    N_eff <- N
  }
  
  # Log-likelihood (Gaussian model)
  sigma2_hat <- RSS / N_eff
  log_lik <- -N_eff / 2 * log(2 * pi * sigma2_hat) - N_eff / 2
  
  # Standard BIC: -2*log(L) + df * log(N)
  BIC <- -2 * log_lik + df * log(N_eff)
  
  # Simplified eBIC: BIC + 2 * xi * df * log(p_dim)
  # When xi = 0, this reduces to regular BIC
  if (p_dim > 0 && df > 0) {
    eBIC <- BIC + 2 * xi * df * log(p_dim)
  } else {
    eBIC <- BIC
  }
  
  return(list(
    eBIC = eBIC,
    BIC = BIC,
    RSS = RSS,
    RSS_type = RSS_type,
    df = df,
    k_u = k_u,
    k_fv = k_fv,
    k_nfv = k_nfv,
    log_likelihood = log_lik,
    sigma2_hat = sigma2_hat,
    N = N_eff,
    n = n,
    J = J,
    m = m
  ))
}


#' @title Sequential eBIC Tuning for Hybrid PCA
#'
#' @description
#' Performs sequential tuning of sparsity parameters using eBIC:
#' 1. First tune gamma_u (sparsity in u)
#' 2. Then tune gamma_nfd | gamma_u* (sparsity in non-functional v)
#' 3. Finally tune gamma_fd | gamma_u*, gamma_nfd* (sparsity in functional v)
#'
#' @param fdata Functional data coefficient matrix (n x J), transformed (C %*% G_half)
#' @param nfdata Non-functional data matrix (n x m)
#' @param hd_obj Hybrid data object
#' @param sparse_tuning_u Vector of candidate gamma values for u
#' @param sparse_tuning_nfd List of candidate gamma values for nfd (one vector per variable)
#' @param sparse_tuning_fd List of candidate gamma values for fd (one vector per variable)
#' @param sparse_tuning_type_u Type of sparsity penalty for u ("soft", "hard", "SCAD")
#' @param sparse_tuning_type_nfd Type of sparsity penalty for nfd
#' @param sparse_tuning_type_fd Type of sparsity penalty for fd
#' @param G_half Square root of Gram matrix
#' @param G_half_inverse Inverse of G_half
#' @param S_smooth Smoothing matrix (or list for GCV case)
#' @param S_2_inverse Inverse of S_2 matrix (or list for GCV case)
#' @param pen_u Logical: penalize u?
#' @param pen_nfd Logical: penalize nfd?
#' @param pen_fd Logical: penalize fd?
#' @param xi eBIC tuning parameter (0 to 1)
#' @param N_type How to compute effective sample size
#' @param RSS_type Character: "direct" or "refitted"
#' @param tol Convergence tolerance
#' @param max_iter Maximum iterations
#' @param pb Progress bar object
#' @param count Current progress count
#' @param verbose Logical: print progress?
#'
#' @return A list containing optimal tuning parameters and eBIC scores
#'
#' @export
eBIC_sequential_hybrid <- function(fdata,
                                   nfdata,
                                   hd_obj = NULL,
                                   sparse_tuning_u = NULL,
                                   sparse_tuning_nfd = NULL,
                                   sparse_tuning_fd = NULL,
                                   sparse_tuning_type_u = "soft",
                                   sparse_tuning_type_nfd = "soft",
                                   sparse_tuning_type_fd = "soft",
                                   G_half = NULL,
                                   G_half_inverse = NULL,
                                   S_smooth = NULL,
                                   S_2_inverse = NULL,
                                   pen_u = FALSE,
                                   pen_nfd = FALSE,
                                   pen_fd = FALSE,
                                   xi = 0.5,
                                   N_type = "product",
                                   RSS_type = "direct",
                                   tol = 1e-4,
                                   max_iter = 1000,
                                   pb = NULL,
                                   count = 0,
                                   verbose = TRUE) {
  
  # Handle S_smooth and S_2_inverse (may be list from GCV preprocessing)
  S_smooth_use <- if (is.list(S_smooth) && !is.matrix(S_smooth)) S_smooth[[1]] else S_smooth
  S_2_inverse_use <- if (is.list(S_2_inverse) && !is.matrix(S_2_inverse)) S_2_inverse[[1]] else S_2_inverse
  
  # Quick exit if nothing to tune
  if (is.null(sparse_tuning_u) &&
      is.null(sparse_tuning_nfd) &&
      is.null(sparse_tuning_fd)) {
    if (!is.null(pb)) setTxtProgressBar(pb, count + 4)
    return(list(
      sparse_tuning_selection_u   = 0,
      sparse_tuning_selection_nfd = NULL,
      sparse_tuning_selection_fd  = NULL,
      eBIC_scores_u               = NULL,
      eBIC_scores_nfd             = NULL,
      eBIC_scores_fd              = NULL
    ))
  }
  
  # Initialize best values
  best_u <- 0
  
  if (!is.null(sparse_tuning_nfd)) {
    d <- length(sparse_tuning_nfd)
    best_nfd <- numeric(d)
    names(best_nfd) <- paste0("Var", seq_len(d))
  } else {
    best_nfd <- NULL
    d <- 0
  }
  
  if (!is.null(sparse_tuning_fd)) {
    dd <- length(sparse_tuning_fd)
    best_fd <- numeric(dd)
    names(best_fd) <- paste0("Var", seq_len(dd))
  } else {
    best_fd <- NULL
    dd <- 0
  }
  
  eBIC_scores_u <- NULL
  eBIC_scores_nfd <- NULL
  eBIC_scores_fd <- NULL
  
  #-----------------------------------------
  # Stage 1: Tune gamma_u
  # tuning_stage = "u": df = k_u, RSS from both components
  #-----------------------------------------
  if (!is.null(sparse_tuning_u) && isTRUE(pen_u)) {
    n_gamma_u <- length(sparse_tuning_u)
    eBIC_scores_u <- numeric(n_gamma_u)
    
    for (idx in seq_len(n_gamma_u)) {
      candidate_u <- sparse_tuning_u[idx]
      
      fit <- init_sequential_hybrid(
        fdata = fdata,
        nfdata = nfdata,
        sparse_tuning_result_u = candidate_u,
        sparse_tuning_result_nfd = best_nfd,
        sparse_tuning_result_fd = best_fd,
        sparse_tuning_type_u = sparse_tuning_type_u,
        sparse_tuning_type_nfd = sparse_tuning_type_nfd,
        sparse_tuning_type_fd = sparse_tuning_type_fd,
        S_smooth = S_smooth_use,
        S_2_inverse = S_2_inverse_use,
        G_half_inverse = G_half_inverse,
        G_half = G_half,
        pen_u = TRUE,
        pen_nfd = pen_nfd,
        pen_fd = pen_fd,
        tol = tol,
        max_iter = max_iter
      )
      
      fv_hat <- fit[[1]]
      nfv_hat <- fit[[2]]
      u_hat <- fit[[3]]
      
      ebic_result <- compute_eBIC_hybrid(
        fdata = fdata,
        nfdata = nfdata,
        u_hat = u_hat,
        fv_hat = fv_hat,
        nfv_hat = nfv_hat,
        hd_obj = hd_obj,
        G_half = G_half,
        xi = xi,
        N_type = N_type,
        RSS_type = RSS_type,
        tuning_stage = "u"
      )
      
      eBIC_scores_u[idx] <- ebic_result$eBIC
    }
    
    best_index_u <- which.min(eBIC_scores_u)
    best_u <- sparse_tuning_u[best_index_u]
    names(eBIC_scores_u) <- as.character(sparse_tuning_u)
    
    count <- count + n_gamma_u
    if (!is.null(pb)) setTxtProgressBar(pb, count)
  } else {
    best_u <- 0
    eBIC_scores_u <- NULL
    count <- count + 1
    if (!is.null(pb)) setTxtProgressBar(pb, count)
  }
  
  #-----------------------------------------
  # Stage 2: Tune gamma_nfd (loop over variables)
  # tuning_stage = "nfd": df = k_nfv, RSS from non-functional only
  #-----------------------------------------
  if (!is.null(sparse_tuning_nfd) && isTRUE(pen_nfd)) {
    eBIC_scores_nfd <- vector("list", length = d)
    
    for (j in seq_len(d)) {
      cand_nfd <- sparse_tuning_nfd[[j]]
      n_gamma_nfd <- length(cand_nfd)
      eBIC_var <- numeric(n_gamma_nfd)
      
      for (idx in seq_len(n_gamma_nfd)) {
        candidate_nfd <- best_nfd
        candidate_nfd[j] <- cand_nfd[idx]
        
        fit <- init_sequential_hybrid(
          fdata = fdata,
          nfdata = nfdata,
          sparse_tuning_result_u = best_u,
          sparse_tuning_result_nfd = candidate_nfd,
          sparse_tuning_result_fd = best_fd,
          sparse_tuning_type_u = sparse_tuning_type_u,
          sparse_tuning_type_nfd = sparse_tuning_type_nfd,
          sparse_tuning_type_fd = sparse_tuning_type_fd,
          S_smooth = S_smooth_use,
          S_2_inverse = S_2_inverse_use,
          G_half_inverse = G_half_inverse,
          G_half = G_half,
          pen_u = pen_u,
          pen_nfd = TRUE,
          pen_fd = pen_fd,
          tol = tol,
          max_iter = max_iter
        )
        
        fv_hat <- fit[[1]]
        nfv_hat <- fit[[2]]
        u_hat <- fit[[3]]
        
        ebic_result <- compute_eBIC_hybrid(
          fdata = fdata,
          nfdata = nfdata,
          u_hat = u_hat,
          fv_hat = fv_hat,
          nfv_hat = nfv_hat,
          hd_obj = hd_obj,
          G_half = G_half,
          xi = xi,
          N_type = N_type,
          RSS_type = RSS_type,
          tuning_stage = "nfd"
        )
        
        eBIC_var[idx] <- ebic_result$eBIC
      }
      
      best_index_nfd <- which.min(eBIC_var)
      best_nfd[j] <- cand_nfd[best_index_nfd]
      eBIC_scores_nfd[[j]] <- eBIC_var
      names(eBIC_scores_nfd[[j]]) <- as.character(cand_nfd)
      
      count <- count + n_gamma_nfd
    }
    
    names(eBIC_scores_nfd) <- paste0("Var", seq_len(d))
    if (!is.null(pb)) setTxtProgressBar(pb, count)
  } else {
    best_nfd <- NULL
    eBIC_scores_nfd <- NULL
    count <- count + 1
    if (!is.null(pb)) setTxtProgressBar(pb, count)
  }
  
  #-----------------------------------------
  # Stage 3: Tune gamma_fd (loop over variables)
  # tuning_stage = "fd": df = k_fv, RSS from functional only
  #-----------------------------------------
  if (!is.null(sparse_tuning_fd) && isTRUE(pen_fd)) {
    eBIC_scores_fd <- vector("list", length = dd)
    
    for (j in seq_len(dd)) {
      cand_fd <- sparse_tuning_fd[[j]]
      n_gamma_fd <- length(cand_fd)
      eBIC_var <- numeric(n_gamma_fd)
      
      for (idx in seq_len(n_gamma_fd)) {
        candidate_fd <- best_fd
        candidate_fd[j] <- cand_fd[idx]
        
        fit <- init_sequential_hybrid(
          fdata = fdata,
          nfdata = nfdata,
          sparse_tuning_result_u = best_u,
          sparse_tuning_result_nfd = best_nfd,
          sparse_tuning_result_fd = candidate_fd,
          sparse_tuning_type_u = sparse_tuning_type_u,
          sparse_tuning_type_nfd = sparse_tuning_type_nfd,
          sparse_tuning_type_fd = sparse_tuning_type_fd,
          S_smooth = S_smooth_use,
          S_2_inverse = S_2_inverse_use,
          G_half_inverse = G_half_inverse,
          G_half = G_half,
          pen_u = pen_u,
          pen_nfd = pen_nfd,
          pen_fd = TRUE,
          tol = tol,
          max_iter = max_iter
        )
        
        fv_hat <- fit[[1]]
        nfv_hat <- fit[[2]]
        u_hat <- fit[[3]]
        
        ebic_result <- compute_eBIC_hybrid(
          fdata = fdata,
          nfdata = nfdata,
          u_hat = u_hat,
          fv_hat = fv_hat,
          nfv_hat = nfv_hat,
          hd_obj = hd_obj,
          G_half = G_half,
          xi = xi,
          N_type = N_type,
          RSS_type = RSS_type,
          tuning_stage = "fd"
        )
        
        eBIC_var[idx] <- ebic_result$eBIC
      }
      
      best_index_fd <- which.min(eBIC_var)
      best_fd[j] <- cand_fd[best_index_fd]
      eBIC_scores_fd[[j]] <- eBIC_var
      names(eBIC_scores_fd[[j]]) <- as.character(cand_fd)
      
      count <- count + n_gamma_fd
    }
    
    names(eBIC_scores_fd) <- paste0("Var", seq_len(dd))
    if (!is.null(pb)) setTxtProgressBar(pb, count)
  } else {
    best_fd <- NULL
    eBIC_scores_fd <- NULL
    count <- count + 1
    if (!is.null(pb)) setTxtProgressBar(pb, count)
  }
  
  return(list(
    sparse_tuning_selection_u = best_u,
    sparse_tuning_selection_nfd = best_nfd,
    sparse_tuning_selection_fd = best_fd,
    eBIC_scores_u = eBIC_scores_u,
    eBIC_scores_nfd = eBIC_scores_nfd,
    eBIC_scores_fd = eBIC_scores_fd
  ))
}