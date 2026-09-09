# Tests for simulate.pop()

# ── Input validation ─────────────────────────────────────────────────────

test_that("errors when sim_config is not a list", {
  expect_error(simulate.pop("not a list", num_years = 5), "sim_config.*list")
})

test_that("errors when sim_config is missing required fields", {
  expect_error(
    simulate.pop(list(max_age = 4), num_years = 5),
    "missing required fields"
  )
})

test_that("errors on invalid num_years", {
  cfg <- quick_config()
  expect_error(simulate.pop(cfg, num_years = 0), "num_years")
  expect_error(simulate.pop(cfg, num_years = -1), "num_years")
  expect_error(simulate.pop(cfg, num_years = 1.5), "num_years")
})

test_that("errors on non-numeric sample_years", {
  cfg <- quick_config()
  expect_error(simulate.pop(cfg, num_years = 5, sample_years = "last"), "sample_years")
})

test_that("warns when a snapshot year falls inside the burn-in window", {
  cfg <- quick_config()
  # burn_in = 2 * max_age = 8. A scalar sample_years means "last N years",
  # so an explicit vector is needed to request an early (in-burn-in) year.
  expect_warning(
    suppressMessages(simulate.pop(cfg, num_years = 5, sample_years = c(1L, 5L))),
    "burn-in"
  )
})

# ── Output structure (no pods) ────────────────────────────────────────────

test_that("simulation length equals burn-in + num_years", {
  cfg <- quick_config()
  sim <- quick_sim(cfg, num_years = 5)

  burn_in     <- 2L * cfg$max_age
  total_years <- burn_in + 5L

  expect_equal(range(sim$pop_summary$year), c(1L, total_years))
  expect_null(sim$pod_to_sp)
  expect_identical(sim$sim_config, cfg)
})

test_that("pop_summary has the expected columns", {
  cfg <- quick_config()
  sim <- quick_sim(cfg, num_years = 5)
  expect_named(sim$pop_summary, c("sex", "age", "N", "year"))
})

test_that("snapshots are stored at the requested years with expected columns", {
  cfg <- quick_config()
  sim <- quick_sim(cfg, num_years = 5, sample_years = 3)

  burn_in     <- 2L * cfg$max_age
  total_years <- burn_in + 5L
  expected_years <- seq.int(total_years - 2L, total_years)

  expect_equal(as.integer(names(sim$snapshots)), expected_years)
  expect_length(sim$snapshots, 3L)

  snap <- sim$snapshots[[1]]
  expect_true(all(c("id", "birth_year", "age", "sex", "mat_age", "mother_id",
                     "father_id", "breed_state", "fertile", "population",
                     "calf_id") %in% names(snap)))
  # s2_year is internal and must be dropped from snapshots
  expect_false("s2_year" %in% names(snap))
})

test_that("sample_years accepts an explicit vector of year indices", {
  cfg <- quick_config()
  sim <- suppressWarnings(quick_sim(cfg, num_years = 5, sample_years = c(9L, 10L)))
  expect_equal(as.integer(names(sim$snapshots)), c(9L, 10L))
})

test_that("sample_years = NULL produces no snapshots", {
  cfg <- quick_config()
  sim <- quick_sim(cfg, num_years = 5, sample_years = NULL)
  expect_length(sim$snapshots, 0L)
})

test_that("the founder generation is fully flushed by the burn-in", {
  cfg <- quick_config()
  sim <- quick_sim(cfg, num_years = 5, sample_years = 3)
  for (snap in sim$snapshots) {
    expect_equal(sum(snap$mother_id == 0L & snap$father_id == 0L), 0L)
  }
})

test_that("all individual ages fall within [0, max_age]", {
  cfg <- quick_config()
  sim <- quick_sim(cfg, num_years = 5, sample_years = 3)
  for (snap in sim$snapshots) {
    expect_true(all(snap$age >= 0L & snap$age <= cfg$max_age))
  }
})

# ── Output structure (with pods) ──────────────────────────────────────────

test_that("pod_to_sp is returned and snapshots carry pod/superpod columns", {
  cfg <- quick_config_pods()
  sim <- quick_sim(cfg, num_years = 5, sample_years = 3)

  expect_true(is.integer(sim$pod_to_sp) || is.numeric(sim$pod_to_sp))
  snap <- sim$snapshots[[1]]
  expect_true(all(c("pod", "superpod") %in% names(snap)))
})

test_that("dependent calves share their mother's pod and superpod", {
  cfg <- quick_config_pods(weaning_age = 2L)
  sim <- quick_sim(cfg, num_years = 5, sample_years = 3)

  for (snap in sim$snapshots) {
    dep <- snap[snap$age < 2L & snap$mother_id != 0L, ]
    if (nrow(dep) == 0L) next
    mother_rows <- match(dep$mother_id, snap$id)
    has_mother  <- !is.na(mother_rows)
    if (!any(has_mother)) next
    expect_equal(dep$pod[has_mother], snap$pod[mother_rows[has_mother]])
    expect_equal(dep$superpod[has_mother], snap$superpod[mother_rows[has_mother]])
  }
})

# ── Density dependence ────────────────────────────────────────────────────

test_that("DD simulation returns a depletion vector of the correct length", {
  cfg <- quick_config_dd()
  sim <- quick_sim(cfg, num_years = 5, sample_years = 2)

  burn_in     <- 2L * cfg$max_age
  total_years <- burn_in + 5L

  expect_true("depletion" %in% names(sim))
  expect_length(sim$depletion, total_years)
  expect_true(all(sim$depletion > 0))
})

test_that("non-DD simulation does not return a depletion vector", {
  cfg <- quick_config()
  sim <- quick_sim(cfg, num_years = 5)
  expect_false("depletion" %in% names(sim))
})

# ── Fishing mortality ──────────────────────────────────────────────────

test_that("F_t = NULL produces no fishing output (backward compatibility)", {
  cfg <- quick_config_dd()
  sim <- quick_sim(cfg, num_years = 5)
  expect_false("F_t" %in% names(sim))
  expect_false("selectivity" %in% names(sim))
})

test_that("scalar F_t is accepted and returned in output", {
  cfg <- quick_config_dd()
  sel <- c(0, rep(1, cfg$max_age))
  sim <- quick_sim(cfg, num_years = 5, F_t = 0.01, selectivity = sel)
  expect_true("F_t" %in% names(sim))
  expect_true("selectivity" %in% names(sim))
  expect_equal(sim$F_t, 0.01)
})

test_that("fishing reduces population size vs no fishing", {
  cfg <- quick_config_dd()
  sel <- c(0, rep(1, cfg$max_age))
  sim_nof <- quick_sim(cfg, num_years = 10, seed = 1111)
  sim_f   <- quick_sim(cfg, num_years = 10, seed = 1111,
                        F_t = 0.05, selectivity = sel)
  # Fished population should be smaller
  final_nof <- sim_nof$pop_summary[year == max(year), sum(N)]
  final_f   <- sim_f$pop_summary[year == max(year), sum(N)]
  expect_lt(final_f, final_nof)
})

test_that("errors on wrong-length F_t vector", {
  cfg <- quick_config_dd()
  sel <- c(0, rep(1, cfg$max_age))
  expect_error(
    simulate.pop(cfg, num_years = 5, F_t = c(0.01, 0.02), selectivity = sel),
    "F_t.*length"
  )
})

test_that("errors when F_t supplied without selectivity", {
  cfg <- quick_config_dd()
  expect_error(
    simulate.pop(cfg, num_years = 5, F_t = 0.01),
    "selectivity.*required"
  )
})

test_that("errors on negative F_t", {
  cfg <- quick_config_dd()
  sel <- c(0, rep(1, cfg$max_age))
  expect_error(
    simulate.pop(cfg, num_years = 5, F_t = -0.1, selectivity = sel),
    "F_t.*>= 0"
  )
})

test_that("errors on selectivity out of [0, 1]", {
  cfg <- quick_config_dd()
  expect_error(
    simulate.pop(cfg, num_years = 5, F_t = 0.01, selectivity = rep(1.5, cfg$max_age + 1)),
    "selectivity.*0 and 1"
  )
})

# ── Depleted initialization ───────────────────────────────────────────

test_that("init_depletion reduces initial population size", {
  cfg <- quick_config_dd()
  sel <- c(0, rep(1, cfg$max_age))
  sim <- quick_sim(cfg, num_years = 5, F_t = 0.01, selectivity = sel,
                    init_depletion = 0.5)
  first_N <- sim$pop_summary[year == 1, sum(N)]
  expect_lt(first_N, cfg$pop_size)
})

test_that("init_depletion requires DD = TRUE", {
  cfg <- quick_config()
  sel <- c(0, rep(1, cfg$max_age))
  expect_error(
    simulate.pop(cfg, num_years = 5, F_t = 0.01, selectivity = sel,
                 init_depletion = 0.5),
    "density_dependence.*TRUE"
  )
})

test_that("init_depletion requires F_t", {
  cfg <- quick_config_dd()
  expect_error(
    simulate.pop(cfg, num_years = 5, init_depletion = 0.5),
    "init_depletion.*F_t"
  )
})

# ── Orphan mortality ──────────────────────────────────────────────────

test_that("no orphaned dependent calves survive in snapshots", {
  cfg <- quick_config_pods(weaning_age = 2L)
  sim <- quick_sim(cfg, num_years = 10, sample_years = 3)
  for (snap in sim$snapshots) {
    wa <- cfg$weaning_age
    dep <- snap[age <= wa & mother_id != 0L]
    if (nrow(dep) > 0L) {
      orphans <- dep[!mother_id %in% snap$id]
      expect_equal(nrow(orphans), 0L)
    }
  }
})

# ── solve_F_sustainable ───────────────────────────────────────────────

test_that("solve_F_sustainable returns positive F for DD config", {
  cfg <- quick_config_dd()
  result <- suppressMessages(suppressWarnings(solve_F_sustainable(cfg)))
  expect_true(result$F_sustainable > 0)
  expect_true(result$F_leslie > 0)
  expect_true(result$interval_at_equilibrium > 0)
  expect_equal(length(result$stable_age_fished), cfg$max_age + 1)
})

test_that("solve_F_sustainable errors on DD = FALSE config", {
  cfg <- quick_config()
  expect_error(solve_F_sustainable(cfg), "density_dependence.*TRUE")
})
