test_that("resultsTable reports successful convergence for each participant", {
  set.seed(1)

  dat <- data.frame(
    subject = rep(1:4, each = 5),
    subscale = rep(c("A", "B"), each = 10),
    direction = rep(c("pos", "neg"), each = 10),
    response = c(2, 3, 4, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 2, 3, 4, 5, 4, 3)
  )

  results <- lasarsmodel(
    data = dat,
    response = response,
    subscale = subscale,
    subject = subject,
    direction = direction,
    rev_label = "neg",
    estDirectPref = TRUE,
    resp_opts = 5,
    subj_est_table = TRUE
  )

  results_df <- results$resultsTable$asDF
  success_count <- sum(results_df$convergence == "Success", na.rm = TRUE)

  expect_true("convergence" %in% names(results_df))
  expect_equal(success_count, nrow(results_df))
  expect_equal(success_count, 4L)
})
