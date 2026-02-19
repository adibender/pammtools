context("Transition probability robustness")

test_that("transition probabilities remain bounded and confidence intervals ordered", {
  fx <- get_multistate_fixture()
  newdata_sorted <- make_transition_newdata(fx$ped, shuffle = FALSE)

  out <- add_trans_prob(newdata_sorted, fx$pam, ci = TRUE)

  expect_true(all(out$trans_prob >= 0 & out$trans_prob <= 1))
  expect_true(all(out$trans_lower >= 0 & out$trans_lower <= 1))
  expect_true(all(out$trans_upper >= 0 & out$trans_upper <= 1))
  expect_true(all(out$trans_lower <= out$trans_prob))
  expect_true(all(out$trans_prob <= out$trans_upper))
  expect_false(anyNA(out$trans_prob))
  expect_false(anyNA(out$trans_lower))
  expect_false(anyNA(out$trans_upper))
})

test_that("transition probabilities are invariant to row order within groups", {
  fx <- get_multistate_fixture()

  sorted_input <- make_transition_newdata(fx$ped, shuffle = FALSE)
  shuffled_input <- make_transition_newdata(fx$ped, shuffle = TRUE, seed = 101)

  out_sorted <- add_trans_prob(sorted_input, fx$pam, ci = FALSE) |>
    dplyr::arrange(.data$transition, .data$tend)
  out_shuffled <- add_trans_prob(shuffled_input, fx$pam, ci = FALSE) |>
    dplyr::arrange(.data$transition, .data$tend)

  expect_equal(
    out_sorted$trans_prob,
    out_shuffled$trans_prob,
    tolerance = 1e-8
  )
})
