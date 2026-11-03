# Fit-wide settings of an Mhpca() fit. Built once in mhpca$initialize() and
# passed down to the algorithm functions as `settings`, instead of threading
# each option through every layer.
#
# To add a new fit-wide option: add the argument to Mhpca() and
# mhpca$initialize(), add it here, and read settings$<name> where it is used.
#
# Every argument defaults to NULL so each path only has to supply (and force)
# the options it uses: the joint/eigen path, for example, never evaluates
# sparse_CV, which mhpca$new() has no default for.
mhpca_settings <- function(smoothing_type = NULL, centerfns = NULL,
                           smooth_GCV = NULL, sparse_CV = NULL,
                           sparse_type_u = NULL, sparse_type_nfd = NULL, sparse_type_fd = NULL,
                           nfold_u = NULL, nfold_nfd = NULL, nfold_fd = NULL,
                           pen_u = NULL, pen_nfd = NULL, pen_fd = NULL,
                           sparse_iter = NULL, cv.pick = NULL,
                           sparse_selection_method = NULL,
                           ebic_xi = NULL, ebic_N_type = NULL, ebic_RSS_type = NULL,
                           tol = NULL, max_iter = NULL, ncor = NULL) {
  list(
    smoothing_type = smoothing_type,
    centerfns = centerfns,
    smooth_GCV = smooth_GCV,
    sparse_CV = sparse_CV,
    sparse_type_u = sparse_type_u,
    sparse_type_nfd = sparse_type_nfd,
    sparse_type_fd = sparse_type_fd,
    nfold_u = nfold_u,
    nfold_nfd = nfold_nfd,
    nfold_fd = nfold_fd,
    pen_u = pen_u,
    pen_nfd = pen_nfd,
    pen_fd = pen_fd,
    sparse_iter = sparse_iter,
    cv.pick = cv.pick,
    sparse_selection_method = sparse_selection_method,
    ebic_xi = ebic_xi,
    ebic_N_type = ebic_N_type,
    ebic_RSS_type = ebic_RSS_type,
    tol = tol,
    max_iter = max_iter,
    ncor = ncor
  )
}
