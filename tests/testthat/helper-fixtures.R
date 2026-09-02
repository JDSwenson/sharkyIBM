# Shared fixtures for fast, deterministic tests of the create.stable.pop() ->
# simulate.pop() -> sample.pop() pipeline.
#
# Life-history parameters here are chosen purely for test speed (small
# max_age, high fecundity, loose convergence tolerances) and are NOT
# biologically realistic for dolphins. Do not reuse them for analysis.

# -- No pod structure, s0 calibration (density_dependence = FALSE) --
quick_config <- function(seed = 4821, max_age = 4L, ...) {
  set.seed(seed)
  survival <- c(0.5, rep(0.8, max_age))
  args <- list(
    max_age         = max_age,
    survival        = survival,
    pop_size        = 300,
    maturity_age    = 1L,
    litter_size     = 4,
    psi_nurse       = 0.3,
    psi_rest        = 0.8,
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
quick_config_pods <- function(seed = 4821, max_age = 4L, ...) {
  quick_config(
    seed = seed, max_age = max_age,
    pod_size = 10, superpod_size = 3,
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
    psi_nurse          = 0.3,
    psi_rest           = 0.8,
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
