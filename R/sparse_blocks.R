# The three sparsity blocks that are tuned in turn: the scores (u), the
# non-functional loadings (nfd) and the functional loadings (fd). Everything
# that differs between blocks is described here, so the tuning loops (CV and
# eBIC) and the formatting of their scores can handle the blocks uniformly.
# To add a block, add an entry here and give its candidates to the tuners.
#
# per_variable  candidates are given per variable (a list of vectors), and a
#               value is selected for each variable; otherwise one vector of
#               candidates and one selected value (u).
# when_skipped  selection returned for a block that is not tuned (no grid or
#               its pen_* flag is FALSE): "zero" -> 0 (no sparsity),
#               "null" -> NULL.
# progress      when the progress bar is updated, per method: once per
#               "block" or after each "variable".
sparse_blocks <- list(
  u = list(
    per_variable = FALSE,
    when_skipped = "zero",
    progress = c(cv = "block", ebic = "block")
  ),
  nfd = list(
    per_variable = TRUE,
    when_skipped = "null",
    progress = c(cv = "block", ebic = "block")
  ),
  fd = list(
    per_variable = TRUE,
    when_skipped = "null",
    progress = c(cv = "variable", ebic = "block")
  )
)

# Starting selection of a block before tuning: 0 for u; for per-variable
# blocks a zero per variable (named Var1, Var2, ...) or NULL without a grid.
initial_selection <- function(block, candidates) {
  if (!sparse_blocks[[block]]$per_variable) return(0)
  if (is.null(candidates)) return(NULL)
  d <- length(candidates)
  stats::setNames(numeric(d), paste0("Var", seq_len(d)))
}

skipped_selection <- function(block) {
  switch(sparse_blocks[[block]]$when_skipped,
    zero = 0,
    null = NULL
  )
}

# `current` (a list with elements u, nfd, fd) with `value` substituted for
# block `block` (variable `j` for per-variable blocks).
with_candidate <- function(current, block, j, value) {
  if (sparse_blocks[[block]]$per_variable) {
    current[[block]][j] <- value
  } else {
    current[[block]] <- value
  }
  current
}

# Index of the selected candidate under the CV rule: "1se" takes the last
# candidate whose CV error is within one SE of the minimum (assumes candidates
# are ordered by increasing sparsity); otherwise the minimum.
pick_cv_index <- function(score, se, cv.pick) {
  if (cv.pick == "1se") {
    min_idx <- which.min(score)
    cutoff <- score[min_idx] + se[min_idx]
    max(which(score <= cutoff))
  } else {
    which.min(score)
  }
}

# Tune the blocks in turn (u, nfd, fd), `passes` times.
#
# candidates  list(u = , nfd = , fd = ) of candidate grids (NULL = no grid).
# pen         list(u = , nfd = , fd = ) of the pen_* flags.
# method      "cv" or "ebic" (selects the per-method behaviour above).
# evaluate    function(block, j, cands, current) returning
#             list(score = <one per candidate>, se = <one per candidate or NULL>);
#             j is NULL for blocks that are not per variable.
# pick        function(score, se) returning the index of the chosen candidate.
# pb, count   progress bar (or NULL) and its running count.
# name_scores name the per-variable score lists Var1, Var2, ... .
#
# Returns the selections, scores and SEs per block, the score of the last
# selected candidate per block (NULL if the block was never tuned), and count.
tune_sparse_blocks <- function(candidates, pen, method, evaluate, pick,
                               pb = NULL, count = 0, passes = 1,
                               name_scores = FALSE) {
  blocks <- names(sparse_blocks)
  empty <- stats::setNames(vector("list", length(blocks)), blocks)
  current <- empty
  for (b in blocks) current[b] <- list(initial_selection(b, candidates[[b]]))
  scores <- ses <- last_score <- empty
  update_progress <- function() if (!is.null(pb)) setTxtProgressBar(pb, count)

  for (pass in seq_len(passes)) {
    for (b in blocks) {
      spec <- sparse_blocks[[b]]
      if (!is.null(candidates[[b]]) && isTRUE(pen[[b]])) {
        if (spec$per_variable) {
          d <- length(candidates[[b]])
          block_scores <- block_ses <- vector("list", length = d)
          for (j in seq_len(d)) {
            cands <- candidates[[b]][[j]]
            ev <- evaluate(b, j, cands, current)
            block_scores[j] <- list(ev$score)
            block_ses[j] <- list(ev$se)
            idx <- pick(ev$score, ev$se)
            current[[b]][j] <- cands[idx]
            last_score[b] <- list(ev$score[idx])
            count <- count + length(cands)
            if (spec$progress[[method]] == "variable") update_progress()
          }
          if (spec$progress[[method]] == "block") update_progress()
          if (name_scores) names(block_scores) <- paste0("Var", seq_len(d))
        } else {
          cands <- candidates[[b]]
          ev <- evaluate(b, NULL, cands, current)
          block_scores <- ev$score
          block_ses <- ev$se
          idx <- pick(ev$score, ev$se)
          current[b] <- list(cands[idx])
          last_score[b] <- list(ev$score[idx])
          count <- count + length(cands)
          update_progress()
        }
        scores[b] <- list(block_scores)
        ses[b] <- list(block_ses)
      } else {
        current[b] <- list(skipped_selection(b))
        scores[b] <- list(NULL)
        ses[b] <- list(NULL)
        count <- count + 1
        update_progress()
      }
    }
  }
  list(selection = current, scores = scores, ses = ses,
       last_score = last_score, count = count)
}

# CV results of one block as tables with columns gamma / CV / SE: a list per
# PC (named PC1, PC2, ...), and within it a table per variable (Var1, ...) for
# per-variable blocks. `grid` is the candidate grid given by the user;
# `scores` / `ses` hold the CV errors / SEs per PC.
format_cv_scores <- function(block, grid, scores, ses) {
  sp <- list(grid)
  if (sparse_blocks[[block]]$per_variable) {
    out <- Map(function(cvs_list, ses_list, sps_list) {
      temp <- Map(function(cvs_i, ses_i, sps_i) {
        cbind(gamma = sps_i, CV = cvs_i, SE = ses_i)
      }, cvs_list, ses_list, sps_list)
      names(temp) <- paste0("Var", seq_len(length(temp)))
      temp
    }, scores, ses, sp)
  } else {
    out <- Map(function(cvs_i, ses_i, sp_i) {
      cbind(gamma = sp_i, CV = cvs_i, SE = ses_i)
    }, scores, ses, sp)
  }
  names(out) <- paste0("PC", seq_len(length(out)))
  out
}
