# Tests for create.stable.pop()
#
# Input-validation tests intentionally supply invalid arguments so that the
# function errors out during validation, before running the (expensive)
# calibration loop -- this keeps them fast regardless of population size.

base_args <- function() {
  list(
    max_age       = 4L,
    survival      = c(0.5, rep(0.8, 4)),
    pop_size      = 300,
    maturity_age  = 1L,
    litter_size   = 4,
    psi_nurse     = 0.3,
    psi_rest      = 0.8,
    male_behavior = NULL
  )
}

# ── Input validation ─────────────────────────────────────────────────────

test_that("errors on invalid max_age", {
  a <- base_args(); a$max_age <- -1
  expect_error(do.call(create.stable.pop, a), "max_age")
})

test_that("errors when survival has wrong length", {
  a <- base_args(); a$survival <- c(0.5, 0.8, 0.8)
  expect_error(do.call(create.stable.pop, a), "survival")
})

test_that("errors when survival values are out of [0, 1]", {
  a <- base_args(); a$survival[2] <- 1.5
  expect_error(do.call(create.stable.pop, a), "survival")
})

test_that("errors on non-positive pop_size", {
  a <- base_args(); a$pop_size <- 0
  expect_error(do.call(create.stable.pop, a), "pop_size")
})

test_that("errors on litter_size < 1", {
  a <- base_args(); a$litter_size <- 0.5
  expect_error(do.call(create.stable.pop, a), "litter_size")
})

test_that("errors when maturity ogive is not non-decreasing", {
  a <- base_args()
  bad_ogive <- c(0, 0.5, 0.3, 0.8, 1)
  a$maturity_age <- bad_ogive
  expect_error(do.call(create.stable.pop, a), "non-decreasing")
})

test_that("errors when maturity_age list is missing an element", {
  a <- base_args()
  a$maturity_age <- list(female = 1L)
  expect_error(do.call(create.stable.pop, a), "female.*male|male.*female")
})

test_that("errors on psi_nurse out of [0, 1]", {
  a <- base_args(); a$psi_nurse <- 1.2
  expect_error(do.call(create.stable.pop, a), "psi_nurse")
})

test_that("errors on female_fraction out of (0, 1)", {
  a <- base_args(); a$female_fraction <- 1
  expect_error(do.call(create.stable.pop, a), "female_fraction")
})

test_that("errors on invalid infertility values", {
  a <- base_args(); a$infertility <- 1
  expect_error(do.call(create.stable.pop, a), "infertility")
})

test_that("errors when pod_size is set without superpod_size", {
  a <- base_args(); a$pod_size <- 10
  expect_error(do.call(create.stable.pop, a), "superpod_size")
})

test_that("errors when stickiness_year is set without pod_size", {
  a <- base_args(); a$stickiness_year <- 0.9
  expect_error(do.call(create.stable.pop, a), "pod_size")
})

test_that("errors on invalid male_behavior", {
  a <- base_args()
  a$pod_size <- 10; a$superpod_size <- 3
  a$male_behavior <- "bogus"
  expect_error(do.call(create.stable.pop, a), "male_behavior")
})

test_that("errors when weaning_age exceeds max_age", {
  a <- base_args(); a$weaning_age <- 10
  expect_error(do.call(create.stable.pop, a), "weaning_age")
})

test_that("density_dependence = TRUE requires dd_max", {
  a <- base_args(); a$density_dependence <- TRUE
  expect_error(do.call(create.stable.pop, a), "dd_max")
})

test_that("density_dependence = TRUE rejects psi values of exactly 0 or 1", {
  a <- base_args()
  a$density_dependence <- TRUE
  a$dd_max <- 3
  a$psi_nurse <- 0
  expect_error(do.call(create.stable.pop, a), "psi_nurse.*psi_rest|density_dependence")
})

test_that("errors on non-positive check_interval", {
  a <- base_args(); a$check_interval <- 0
  expect_error(do.call(create.stable.pop, a), "check_interval")
})

test_that("warns (but does not error) when psi_nurse > psi_rest", {
  a <- base_args()
  a$psi_nurse   <- 0.9
  a$psi_rest    <- 0.1
  a$max_windows <- 1L  # keep it fast; convergence irrelevant to this test

  # A non-convergence warning is also expected at max_windows = 1; muffle it
  # so only the psi_nurse/psi_rest warning under test reaches expect_warning().
  run <- function() {
    withCallingHandlers(
      do.call(create.stable.pop, a),
      warning = function(w) {
        if (grepl("Did not converge", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  expect_warning(suppressMessages(run()), "psi_nurse.*psi_rest")
})

# ── Output structure: s0 calibration (density_dependence = FALSE) ─────────

test_that("s0 calibration returns a well-formed config", {
  cfg <- quick_config()

  expect_type(cfg, "list")
  expect_false(cfg$density_dependence)
  expect_true(is.numeric(cfg$s0) && cfg$s0 > 0 && cfg$s0 < 1)
  expect_length(cfg$survival, cfg$max_age + 1L)
  expect_equal(cfg$survival[1], cfg$s0)
  expect_true(cfg$final_N > 0)
  expect_true(cfg$years_simulated > 0)

  # Pass-through parameters are preserved verbatim
  expect_equal(cfg$max_age, 4L)
  expect_equal(cfg$litter_size, 4)
  expect_equal(cfg$psi_nurse, 0.3)
  expect_equal(cfg$psi_rest, 0.8)
})

test_that("s0 calibration is reproducible with a fixed seed", {
  cfg1 <- quick_config(seed = 917)
  cfg2 <- quick_config(seed = 917)
  expect_equal(cfg1$s0, cfg2$s0)
  expect_equal(cfg1$final_N, cfg2$final_N)
})

test_that("pod/superpod parameters pass through when supplied", {
  cfg <- quick_config_pods()
  expect_equal(cfg$pod_size, 10)
  expect_equal(cfg$superpod_size, 3)
  expect_equal(cfg$stickiness_year, 0.9)
  expect_equal(cfg$male_behavior, "random")
  expect_equal(cfg$weaning_age, 2L)
})

test_that("sex-specific maturity ogives (list input) are preserved", {
  ogive_f <- plogis(0:4, location = 1, scale = 0.4)
  cfg <- quick_config(maturity_age = list(female = ogive_f, male = 1L))
  expect_type(cfg$maturity_age, "list")
  expect_equal(cfg$maturity_age$female, ogive_f)
  expect_equal(cfg$maturity_age$male, 1L)
})

test_that("infertility is echoed back unchanged", {
  cfg <- quick_config(infertility = c(0.1, 0.05))
  expect_equal(cfg$infertility, c(0.1, 0.05))
})

# ── Output structure: density dependence (density_dependence = TRUE) ──────

test_that("DD calibration returns a well-formed config", {
  cfg <- quick_config_dd()

  expect_true(cfg$density_dependence)
  expect_true(is.na(cfg$s0_leslie))
  expect_true(is.numeric(cfg$theta_shift))
  expect_true(cfg$psi_nurse_K > 0 && cfg$psi_nurse_K < 1)
  expect_true(cfg$psi_rest_K  > 0 && cfg$psi_rest_K  < 1)
  expect_true(cfg$K_1plus > 0)
  # Survival is used as-supplied in DD mode (not recalibrated)
  expect_equal(cfg$survival, c(0.55, rep(0.8, 4)))
  expect_equal(cfg$s0, 0.55)
})
