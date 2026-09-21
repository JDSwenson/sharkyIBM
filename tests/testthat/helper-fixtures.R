# Shared fixtures for fast, deterministic tests of the create.stable.pop() ->
# simulate.pop() -> sample.pop() pipeline.
#
# Life-history parameters here are chosen purely for test speed (small
# max_age, high fecundity, loose convergence tolerances) and are NOT
# biologically realistic for dolphins. Do not reuse them for analysis.

# -- No pod structure, s0 calibration (density_dependence = FALSE) --
# pop_size = 900 (bumped up from 300): at pop_size = 300, this fixture's
# combination of a tiny population, litter_size = 4, and maturity_age = 1
# is demographically volatile enough that some seeds/RNG draw sequences run
# it to extinction well within the test horizon -- unrelated to correctness,
# just an artifact of small-N stochasticity. See quick_config_pods() below
# for why pod_size/superpod_size were also bumped, not just pop_size.
quick_config <- function(seed = 4821, max_age = 4L, ...) {
  set.seed(seed)
  survival <- c(0.5, rep(0.8, max_age))
  args <- list(
    max_age         = max_age,
    survival        = survival,
    pop_size        = 900,
    maturity_age    = 1L,
    litter_size     = 4,
    psi_nurse       = 0.3,
    psi_rest        = 0.8,
    weaning_age     = 1L,
    male_behavior   = NULL,
    check_interval  = 3L,
    growth_tol      = 0.05,
    stable_required = 2L,
    max_windows     = 8L
  )
  extra <- list(...)
  args[names(extra)] <- extra
  suppressMessages(suppressWarnings(do.call(create.stable.pop, args)))
}

# -- Pod / superpod structure, s0 calibration --
# pod_size = 20, superpod_size = 25 (bumped up from 10, 3): superpod
# occupancy drifts unevenly over time under stickiness_year-driven
# emigration -- with many small superpods (e.g. the old pod_size = 10,
# superpod_size = 3, giving ~20 superpods of nominal size 30), it's fairly
# likely that at least one ends up nearly empty by the time it's sampled,
# which makes `sample.pop(..., sampling = "superpod")` tests that expect an
# exact sample_size flaky (whichever superpod happens to be drawn for a
# trip might not have enough individuals left). This combination keeps
# exactly 2 superpods (so tests that need >=2, e.g. superpod_pool
# restriction, still run rather than skip) while keeping each comfortably
# larger than any sample_size used in the test suite.
quick_config_pods <- function(seed = 4821, max_age = 4L, ...) {
  quick_config(
    seed = seed, max_age = max_age,
    pod_size = 20, superpod_size = 25,
    stickiness_year = 0.9,
    male_behavior = "random",
    weaning_age = 2L,
    ...
  )
}

# -- Density dependence (density_dependence = TRUE) --
quick_config_dd <- function(seed = 4821, max_age = 4L, ...) {
  set.seed(seed)
  survival <- c(0.55, rep(0.8, max_age))
  args <- list(
    max_age            = max_age,
    survival           = survival,
    pop_size           = 300,
    maturity_age       = 1L,
    litter_size        = 4,
    rho                = -1.5,
    weaning_age        = 1L,
    male_behavior      = NULL,
    density_dependence = TRUE,
    dd_max             = 3,
    check_interval     = 3L,
    growth_tol         = 0.05,
    stable_required    = 2L,
    max_windows        = 8L
  )
  extra <- list(...)
  args[names(extra)] <- extra
  suppressMessages(suppressWarnings(do.call(create.stable.pop, args)))
}

quick_sim <- function(cfg, num_years = 5L, sample_years = 3L, seed = 4821, ...) {
  set.seed(seed)
  suppressMessages(suppressWarnings(
    simulate.pop(cfg, num_years = num_years, sample_years = sample_years, ...)
  ))
}
