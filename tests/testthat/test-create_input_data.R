test_that("create.pop.data basic functionality works", {
  # Simple, feasible life-history
  max_age <- 5
  survival <- rep(0.9, max_age)
  pop_number <- 2
  pop_size <- c(100, 200)
  mating_periodicity <- 1
  maturity_age <- 2
  litter_size <- 4
  female_fraction <- 0.5

  out <- create.pop.data(
    max_age = max_age,
    survival = survival,
    pop_number = pop_number,
    pop_size = pop_size,
    mating_periodicity = mating_periodicity,
    maturity_age = maturity_age,
    litter_size = litter_size,
    female_fraction = female_fraction,
    stable = TRUE
  )

  # Basic return structure
  expect_true(is.list(out))
  expect_named(out, c("numbers_at_age", "survival", "fecundity", "s0"))

  # numbers_at_age checks
  expect_true(is.data.frame(out$numbers_at_age))
  expect_equal(nrow(out$numbers_at_age), (max_age + 1) * pop_number)

  # Ages should be 0:max_age
  expect_equal(sort(unique(out$numbers_at_age$age)), 0:max_age)

  # Populations numbered 1..pop_number
  expect_equal(sort(unique(out$numbers_at_age$population)), 1:pop_number)

  # Population totals preserved
  pop_totals <- tapply(out$numbers_at_age$N, out$numbers_at_age$population, sum)
  expect_equal(as.numeric(pop_totals), pop_size)

  # Survival vector length is max_age + 1
  expect_equal(length(out$survival), max_age + 1)

  # max_age survival must be 0
  expect_equal(out$survival[max_age + 1], 0)

  # fecundity length correct
  expect_equal(length(out$fecundity), max_age + 1)

  # Fecundity before maturity_age is 0
  expect_true(all(out$fecundity[1:maturity_age] == 0))

  # Fecundity at and after maturity_age > 0
  expect_true(all(out$fecundity[(maturity_age+1):(max_age+1)] > 0))

  # s0 must be between 0 and 1
  expect_true(out$s0 > 0 && out$s0 < 1)
})


test_that("create.pop.data errors correctly when survival/fecundity cannot sustain lambda = 1", {
  max_age <- 5
  survival <- rep(0.2, max_age)  # too low to ever be stable
  pop_number <- 1
  pop_size <- 100
  mating_periodicity <- 2
  maturity_age <- 4
  litter_size <- 1  # too low
  female_fraction <- 0.5

  expect_error(
    create.pop.data(
      max_age = max_age,
      survival = survival,
      pop_number = pop_number,
      pop_size = pop_size,
      mating_periodicity = mating_periodicity,
      maturity_age = maturity_age,
      litter_size = litter_size,
      female_fraction = female_fraction,
      stable = TRUE
    ),
    regexp = "Unable to compute YOY survival"
  )
})


test_that("create.pop.data works when stable = FALSE and YOY survival supplied", {
  max_age <- 6
  survival <- rep(0.9, max_age + 1)  # includes s0 here
  pop_number <- 1
  pop_size <- 500
  mating_periodicity <- 1
  maturity_age <- 3
  litter_size <- 6
  YOY_survival <- 0.4

  out <- create.pop.data(
    max_age = max_age,
    survival = survival,
    pop_number = pop_number,
    pop_size = pop_size,
    mating_periodicity = mating_periodicity,
    maturity_age = maturity_age,
    litter_size = litter_size,
    YOY_survival = YOY_survival,
    stable = FALSE
  )

  # Ensure the overridden YOY survival was used
  expect_equal(out$survival[1], YOY_survival)

  # max_age survival must still be 0
  expect_equal(out$survival[max_age + 1], 0)

  # Population total preserved
  expect_equal(sum(out$numbers_at_age$N), pop_size)
})
