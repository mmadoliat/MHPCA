# Numeric snapshots with a tolerance, so tiny floating-point noise does not fail
# a test but any real change in results does.
#   small objects -> "json2": readable in _snaps/*.md, keeps attributes
#   large objects -> "serialize": compact, not human-readable
# ("deparse" is not used: it is re-parsed in an empty environment, so values
# such as 4:5 or character(0) cannot be read back.)
snap_tol <- 1e-6

expect_snapshot_small <- function(x) {
  # A few objects (e.g. zero-length named vectors inside lists) do not survive
  # a json2 round trip; those are stored with "serialize" instead.
  style <- if (json2_roundtrips(x)) "json2" else "serialize"
  testthat::expect_snapshot_value(x, style = style, tolerance = snap_tol)
}

json2_roundtrips <- function(x) {
  back <- tryCatch(
    jsonlite::unserializeJSON(jsonlite::serializeJSON(x, pretty = TRUE)),
    error = function(e) structure(list(), class = "roundtrip_failed")
  )
  length(waldo::compare(x, back, tolerance = snap_tol)) == 0
}

expect_snapshot_large <- function(x) {
  testthat::expect_snapshot_value(x, style = "serialize", tolerance = snap_tol)
}
