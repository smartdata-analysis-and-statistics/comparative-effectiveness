# Chapter 10: compatibility helper for mixgb 2.2.3.
# The installed mixgb_cv() uses N.class before defining it in the
# multiclass branch. Create a local copy and initialize N.class first.
# Do not modify the installed package or the cross-validation parameters.

mixgb_cv_ch10 <- local({
  f <- mixgb::mixgb_cv
  original_body <- paste(deparse(body(f), width.cutoff = 500L), collapse = "\n")

  target <- "xgb.params$num_class <- N.class"
  initialization <- "N.class <- length(levels(cc.data[[response]]))"

  occurrences <- lengths(regmatches(original_body, gregexpr(target, original_body, fixed = TRUE)))
  if (occurrences != 1L) {
    stop("Unexpected mixgb_cv() implementation: check the installed mixgb version.")
  }
  if (!grepl(initialization, original_body, fixed = TRUE)) {
    stop("Unexpected mixgb_cv() implementation: N.class initialization was not found.")
  }

  patched_body <- sub(
    target,
    paste0(initialization, "\n        ", target),
    original_body,
    fixed = TRUE
  )
  body(f) <- parse(text = patched_body)[[1L]]
  f
})
