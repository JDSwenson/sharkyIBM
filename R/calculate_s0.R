#' Find age-0 survival (s0) that stabilizes an age-structured population
#'
#' This function uses a two-phase approach to identify a value for survival of age-0 individuals (s0) that will result in a stable population over the course of a long individual-based simulation run via \code{simulate.pop}. The two phases include:
#' 1. **Leslie matrix** — builds a deterministic Leslie matrix from the supplied
#'    life-history parameters and solves analytically for the s0 that gives
#'    \eqn{\lambda = 1}. This is used as the starting point for the second phase.
#' 2. **Iterative simulation** — runs the stochastic individual-based simulation
#'    forward in time, adjusting s0 via bisection until convergence. Convergence
#'    occurs when \eqn{\lambda} = 1 +/- \code{growth_tol} for the number of
#'    consecutive cycles specified by \code{stable_required}.
#' To streamline the full simulation workflow, this function accepts as input all life history/biology parameters that are relevant to the population simulation - whether needed for this step or not - and outputs them in a list that is used as input to \code{simulate.pop}.
#'
#' Calculating s0 from a Leslie matrix with a dominant eigenvalue of 1 does not
#' account for the stochasticity of the individual-based simulation and is
#' unlikely to produce a stable population over the course of the simulation.
#' However, it is an excellent starting point to refine empirically via
#' iterative simulation. Overall, the two-phase approach implemented in this function reduces computation time while returning a value for s0 that will maintain approximately stable population growth for the full simulation.
#'
#' @param max_age Integer. Maximum age in the population. Any individuals that
#'   remain alive at max_age will be removed from the population before their
#'   next opportunity to mate.
#' @param survival Numeric vector of length \code{max_age + 1}. Annual survival
#'   probabilities for ages 0 through \code{max_age}. The age-0 value is a
#'   placeholder; it will be replaced by the optimized estimate.
#' @param pop_size Integer. Total starting population size. This number will be
#'   parsed into age-classes based on the stable age distribution for the
#'   population.
#' @param maturity_age Maturity specification. Can be:
#'   \itemize{
#'     \item An integer: knife-edged maturity at that age for both sexes.
#'     \item A numeric vector of length \code{max_age + 1}: cumulative
#'       probability of being mature at each age (0 through \code{max_age}),
#'       applied to both sexes.
#'     \item A list with \code{female} and \code{male} elements, each of which can be an integer or a numeric vector as above.
#'   }
#' @param litter_size Numeric. Mean litter size per breeding female.
#' @param psi_2 Numeric in \[0, 1\]. Probability that a female with a dependent
#'   calf (breeding state S2) becomes pregnant again the following year.
#'   Low values (~0.05-0.15) mean females rarely conceive while nursing.
#'   In combination with \code{psi_3}, this parameter is used to model
#'   reproduction as a Markov process similar to Jacobson et al. (2026). Default 0.1.
#' @param psi_3 Numeric in \[0, 1\]. Probability that a resting female
#'   (breeding state S3) becomes pregnant the following year. Controls how
#'   quickly females re-enter the breeding cycle after weaning. Default 0.7, which together with the default value for psi_2, results in a realized breeding interval of roughly 3-4 years, consistent with empirical observations of small-medium sized delphinids.
#' @param num_mates Integer vector or scalar. Number of mates per female.
#'   Not used in \code{calculate.s0} but captured for pass-through to
#'   \code{simulate.pop}.
#' @param female_fraction Numeric in (0, 1). Proportion female offspring i.e., the sex ratio of newborns
#' @param infertility Numeric (scalar or length-2 vector \code{c(female,
#'   male)}). Proportion of individuals that are permanently infertile.
#'   Default 0.
#' @param pod_size Integer or NULL. Initial number of individuals per pod. Over the course of the simulation, pod sizes will vary according to population dynamics.
#' @param superpod_size Integer or NULL. Number of pods per superpod.
#' @param stickiness_year Numeric (scalar or length-2 vector) or NULL.
#'   Between-year superpod fidelity. Can be made sex-specific via \code{c(female, male)}. Lower values mean that the individual is more likely to leave its superpod to join another superpod between simulation years.
#' @param male_behavior Character or NULL. Options are \code{"random"} (default) or \code{"strong_bull"}. \itemize{
#'     \item \code{random} specifies that all males have the chance to reproduce with any female within their superpod.
#'     \item \code{strong_bull} means that only one male can reproduce with females within its superpod each year. When dynamics are governed by age-based processes (as opposed to length-based processes), the bull will be chosen randomly from the mature males within each superpod and will remain the sole mating male within that superpod until his demise, at which point another mature male will be chosen for that superpod.
#'     \item \code{NULL} lumps all mature males across all superpods into a single pool of potential fathers with which all females - regardless of superpod - may reproduce.
#'     }
#' @param max_females Integer or NULL. Maximum number of females a single male
#'   can mate with per year. NULL (default) means no cap. Only used in
#'   \code{simulate.pop}; captured here for pass-through.
#' @param weaning_age Integer or NULL. Age at independence from mother. Prior to
#'   this age, each surviving calf will move with its mother anytime the mother
#'   changes superpods.
#' @param check_interval Integer. Years between population growth assessments to determine stability. Default 10.
#' @param growth_tol Numeric. Stability tolerance i.e., how far from 0 can the annual population growth rate vary and still be considered "stable".  Default 0.001.
#' @param stable_required Integer. Consecutive simulation years of stability required to determine that the population is stable. Default 5.
#' @param max_windows Integer. Safety cap on windows in case stable population
#'   growth is never achieved. Prevents the code from hanging indefinitely if it proves impossible to achieve stability via adjusting s0. Default 100.
#'
#' @return A list (returned invisibly) with calibration results and all shared
#'   parameters needed by \code{simulate.pop}.
#'
#' @importFrom data.table data.table set rbindlist
#' @importFrom stats runif rpois uniroot
#' @export
calculate.s0 <- function(max_age,
                         survival,
                         pop_size,
                         maturity_age,
                         litter_size,
                         psi_2 = 0.1,
                         psi_3 = 0.7,
                         num_mates       = 1L,
                         female_fraction = 0.5,
                         infertility     = 0,
                         pod_size        = NULL,
                         superpod_size   = NULL,
                         stickiness_year = NULL,
                         male_behavior   = "random",
                         max_females     = NULL,
                         weaning_age     = NULL,
                         check_interval  = 10L,
                         growth_tol      = 0.001,
                         stable_required = 5L,
                         max_windows     = 100L) {

  # ═══════════════════════════════════════════════════════════════════════════
  # INPUT VALIDATION
  # ═══════════════════════════════════════════════════════════════════════════
  # Check the most common user errors upfront so they fail with clear messages
  # rather than cryptic indexing errors deep in the simulation loop.

  # -- Core demographic parameters --
  if (!is.numeric(max_age) || length(max_age) != 1L || max_age < 1 ||
      max_age != round(max_age))
    stop("`max_age` must be a positive integer.")

  if (!is.numeric(survival) || length(survival) != max_age + 1L)
    stop("`survival` must be a numeric vector of length max_age + 1 (",
         max_age + 1L, "). Got length ", length(survival), ".")
  if (any(survival < 0 | survival > 1))
    stop("All `survival` values must be between 0 and 1.")

  if (!is.numeric(pop_size) || length(pop_size) != 1L || pop_size < 1)
    stop("`pop_size` must be a positive number.")

  if (!is.numeric(litter_size) || length(litter_size) != 1L || litter_size < 1)
    stop("`litter_size` must be >= 1.")

  # -- Maturity specification --
  # Validate structure: must be a scalar integer, a vector of length max_age+1,
  # or a list with $female and $male elements (each of which is one of the above).
  validate_maturity_input <- function(x, label) {
    if (length(x) == 1L && is.numeric(x) && x == round(x)) {
      if (x < 0 || x > max_age)
        stop(label, ": knife-edged maturity age must be between 0 and max_age (",
             max_age, "). Got ", x, ".")
    } else if (is.numeric(x)) {
      if (length(x) != max_age + 1L)
        stop(label, ": ogive vector must have length max_age + 1 (",
             max_age + 1L, "). Got length ", length(x), ".")
      if (any(x < 0 | x > 1))
        stop(label, ": ogive values must be between 0 and 1.")
      if (is.unsorted(x))
        stop(label, ": ogive must be non-decreasing (cumulative probability).")
    } else {
      stop(label, ": must be an integer or a numeric vector of length ",
           max_age + 1L, ".")
    }
  }

  if (is.list(maturity_age)) {
    if (is.null(maturity_age$female) || is.null(maturity_age$male))
      stop("`maturity_age` list must have both `female` and `male` elements.")
    validate_maturity_input(maturity_age$female, "maturity_age$female")
    validate_maturity_input(maturity_age$male,   "maturity_age$male")
  } else {
    validate_maturity_input(maturity_age, "maturity_age")
  }

  # -- Breeding cycle probabilities --
  if (!is.numeric(psi_2) || length(psi_2) != 1L || psi_2 < 0 || psi_2 > 1)
    stop("`psi_2` must be a single number between 0 and 1.")
  if (!is.numeric(psi_3) || length(psi_3) != 1L || psi_3 < 0 || psi_3 > 1)
    stop("`psi_3` must be a single number between 0 and 1.")

  # -- Sex ratio --
  if (!is.numeric(female_fraction) || length(female_fraction) != 1L ||
      female_fraction <= 0 || female_fraction >= 1)
    stop("`female_fraction` must be between 0 and 1 (exclusive).")

  # -- Infertility --
  if (!is.numeric(infertility) || !length(infertility) %in% c(1L, 2L))
    stop("`infertility` must be a numeric scalar or length-2 vector.")
  if (any(infertility < 0 | infertility >= 1))
    stop("`infertility` values must be >= 0 and < 1.")

  # -- Pod / superpod consistency --
  if (!is.null(pod_size) && is.null(superpod_size))
    stop("If `pod_size` is specified, `superpod_size` must also be specified.")
  if (is.null(pod_size) && !is.null(superpod_size))
    stop("If `superpod_size` is specified, `pod_size` must also be specified.")
  if (!is.null(pod_size) && (pod_size < 1 || pod_size != round(pod_size)))
    stop("`pod_size` must be a positive integer.")
  if (!is.null(superpod_size) && (superpod_size < 1 ||
                                   superpod_size != round(superpod_size)))
    stop("`superpod_size` must be a positive integer.")

  # -- Stickiness --
  if (!is.null(stickiness_year)) {
    if (!is.numeric(stickiness_year) || !length(stickiness_year) %in% c(1L, 2L))
      stop("`stickiness_year` must be a numeric scalar or length-2 vector.")
    if (any(stickiness_year < 0 | stickiness_year > 1))
      stop("`stickiness_year` values must be between 0 and 1.")
    if (is.null(pod_size))
      stop("`stickiness_year` requires `pod_size` and `superpod_size`.")
  }

  # -- Male behavior --
  if (!is.null(male_behavior) &&
      !male_behavior %in% c("random", "strong_bull"))
    stop('`male_behavior` must be NULL, "random", or "strong_bull".')
  if (!is.null(male_behavior) && is.null(pod_size))
    stop("`male_behavior` requires `pod_size` and `superpod_size`.")

  # -- Weaning age --
  if (!is.null(weaning_age)) {
    if (!is.numeric(weaning_age) || length(weaning_age) != 1L ||
        weaning_age < 0 || weaning_age != round(weaning_age))
      stop("`weaning_age` must be a non-negative integer.")
    if (weaning_age > max_age)
      stop("`weaning_age` (", weaning_age, ") must be <= max_age (", max_age, ").")
  }

  # -- Convergence parameters --
  if (check_interval < 1L) stop("`check_interval` must be >= 1.")
  if (growth_tol <= 0)     stop("`growth_tol` must be > 0.")
  if (stable_required < 1) stop("`stable_required` must be >= 1.")
  if (max_windows < 1)     stop("`max_windows` must be >= 1.")

  # ═══════════════════════════════════════════════════════════════════════════
  # INTERNAL HELPERS
  # ═══════════════════════════════════════════════════════════════════════════

  # Dominant eigenvalue of matrix A (= asymptotic population growth rate λ).
  mat_lambda <- function(A) Mod(eigen(A, only.values = TRUE)$values[1])

  # Stable age distribution from matrix A (right eigenvector of λ, normalised).
  mat_stable <- function(A) {
    w <- Mod(eigen(A)$vectors[, 1])
    w / sum(w)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PARSE MATURITY SPECIFICATION
  # ═══════════════════════════════════════════════════════════════════════════
  # Convert user input into two ogive vectors (one per sex), each of length

  # max_age + 1, giving P(mature | age) for ages 0:max_age.  A knife-edged
  # integer is expanded into a step-function ogive (0 below, 1 at/above).

  make_ogive <- function(x, max_a) {
    n <- max_a + 1L
    if (length(x) == 1L && x == round(x)) {
      # Knife-edged: everyone matures at exactly this age
      ogive <- rep(0, n)
      if (x <= max_a) ogive[seq(x + 1L, n)] <- 1
      return(ogive)
    }
    # Otherwise x is already a cumulative probability vector
    if (length(x) == n) return(x)
    stop("maturity_age: ogive vector must have length max_age + 1 (", n, ").")
  }

  if (is.list(maturity_age)) {
    ogive_f <- make_ogive(maturity_age$female, max_age)
    ogive_m <- make_ogive(maturity_age$male, max_age)
  } else {
    ogive_f <- make_ogive(maturity_age, max_age)
    ogive_m <- ogive_f
  }

  # Derive a probability mass function from the ogive CDF so we can sample
  # individual maturity ages at birth.  The PMF gives P(mature at exactly age a).
  # If the ogive doesn't reach 1.0, the leftover probability means "never mature"
  # — those individuals get mat_age = max_age + 1 (unreachable, so they never breed).
  sample_mat_ages <- function(ogive, n) {
    pmf     <- diff(c(0, ogive))
    p_never <- 1 - sum(pmf)
    if (p_never > 0.001) {
      # Some fraction will never mature — add a "never" category
      pmf  <- c(pmf, p_never)
      ages <- c(0:max_age, max_age + 1L)
    } else {
      # Absorb rounding error into the last age class
      ages <- 0:max_age
      pmf[length(pmf)] <- pmf[length(pmf)] + p_never
    }
    sample(ages, n, replace = TRUE, prob = pmf)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PARSE INFERTILITY
  # ═══════════════════════════════════════════════════════════════════════════
  # Split into sex-specific rates.  Infertile individuals are permanently
  # unable to reproduce — the flag is assigned at birth and never changes.

  if (length(infertility) == 1L) {
    infertility_f <- infertility_m <- infertility
  } else {
    infertility_f <- infertility[1]
    infertility_m <- infertility[2]
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # MARKOV BREEDING CYCLE — STATIONARY DISTRIBUTION
  # ═══════════════════════════════════════════════════════════════════════════
  # The breeding cycle is modelled as a three-state Markov chain:
  #   S1 = pregnant, S2 = with dependent calf, S3 = resting (no calf)
  #
  # Transition matrix Ψ (rows = current state, cols = next state):
  #        S1(preg)  S2(calf)  S3(rest)
  #   S1 [   0         1         0     ]   ← pregnant females always give birth
  #   S2 [  ψ₂         0       1-ψ₂   ]   ← nursing females may conceive or wean
  #   S3 [  ψ₃         0       1-ψ₃   ]   ← resting females may conceive or stay
  #
  # At equilibrium the stationary proportions are:
  #   π₁ = π₂ = ψ₃ / (2·ψ₃ + 1 - ψ₂)      [fraction pregnant = fraction with calf]
  #   π₃ = (1-ψ₂)·π₁ / ψ₃                   [fraction resting]
  #
  # π₂ is the fraction of mature females that breed in any given year, which
  # determines the fecundity row of the Leslie matrix below.

  pi_denom <- 2 * psi_3 + 1 - psi_2
  pi_1 <- psi_3 / pi_denom          # proportion pregnant
  pi_2 <- pi_1                       # proportion with calf (= breeding)
  pi_3 <- (1 - psi_2) * pi_1 / psi_3 # proportion resting

  message(sprintf(
    "Breeding cycle: psi_2=%.2f, psi_3=%.2f -> %.1f%% breeding/yr, avg interval=%.1f yrs",
    psi_2, psi_3, pi_2 * 100, 1 / pi_2
  ))

  # ═══════════════════════════════════════════════════════════════════════════
  # PHASE 1: LESLIE-MATRIX STARTING ESTIMATE FOR s0
  # ═══════════════════════════════════════════════════════════════════════════
  # Build a standard Leslie matrix where:
  #   - Row 1 = age-specific fecundity (expected female offspring per female)
  #   - Sub-diagonal = age-specific survival probabilities
  # Then solve for the s0 (age-0 survival) that makes the dominant
  # eigenvalue λ = 1 (i.e., zero population growth).
  #
  # Fecundity at each age = ogive(age) × litter_size × female_fraction × π₂
  #                         × (1 − infertility_f)
  # The ogive gives the fraction of females that are mature at each age.
  # π₂ is the fraction of mature females that breed each year (Markov stationary).
  # (1 − infertility_f) adjusts for the fraction that can never breed.

  n_classes <- max_age + 1L
  A <- matrix(0, nrow = n_classes, ncol = n_classes)

  ff <- litter_size * female_fraction * pi_2 * (1 - infertility_f)
  # Shift the ogive right by one year: newly mature females enter the breeding
  # cycle at S3 (resting) and cannot breed in their first year of maturity.
  # Only females who were already mature (ogive value at the PREVIOUS age)
  # contribute fecundity.  This aligns the Leslie matrix with the IBM's
  # one-year gestation delay for newly mature females.
  ogive_f_leslie <- c(0, ogive_f[seq_len(max_age)])
  f_vec <- ogive_f_leslie * ff   # age-specific fecundity (length max_age + 1)

  # Row 1: fecundity for each age class
  A[1, ] <- f_vec
  # Sub-diagonal: survival from age (i-1) to age i.
  # survival[i] is the probability of surviving at age (i-1), 1-indexed.
  for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i]

  # Check that a valid s0 exists in [0.01, 0.99].  If λ > 1 at both extremes
  # (or < 1 at both), no root exists and the life-history parameters are
  # infeasible for a stable population.
  A_lo <- A; A_lo[2, 1] <- 0.01
  A_hi <- A; A_hi[2, 1] <- 0.99

  if ((mat_lambda(A_lo) - 1) * (mat_lambda(A_hi) - 1) > 0) {
    stop(
      "No s0 in [0.01, 0.99] can stabilise this population. ",
      sprintf("Lambda = %.4f at s0=0.01 and %.4f at s0=0.99.",
              mat_lambda(A_lo), mat_lambda(A_hi))
    )
  }

  # Find the root: the s0 where λ(A) = 1.
  s0_leslie <- uniroot(
    function(s0) { A[2, 1] <- s0; mat_lambda(A) - 1 },
    interval = c(0.01, 0.99)
  )$root

  # Plug the analytical s0 back in and extract the stable age distribution.
  # This distribution determines how many individuals of each age to create
  # at initialisation, so the population starts near equilibrium.
  A[2, 1]  <- s0_leslie
  stable_A <- mat_stable(A)

  message(sprintf(
    "Phase 1 — Leslie estimate: s0 = %.4f  (lambda = %.6f)",
    s0_leslie, mat_lambda(A)
  ))

  # ═══════════════════════════════════════════════════════════════════════════
  # PHASE 2: INITIALISE THE INDIVIDUAL-BASED POPULATION
  # ═══════════════════════════════════════════════════════════════════════════
  # Create an initial population of individuals using the stable age
  # distribution from Phase 1.  Each individual gets:
  #   - An age (from the stable distribution)
  #   - A sex (randomly assigned based on female_fraction)
  #   - A personal maturity age (sampled from the ogive PMF)
  #   - A fertility flag (TRUE/FALSE based on infertility rate)
  #   - A breeding state (mature fertile females get a state from the Markov
  #     stationary distribution; everyone else gets NA)

  use_pods <- !is.null(pod_size)
  wa       <- if (is.null(weaning_age)) 0L else weaning_age

  # Allocate individuals to age classes proportional to the stable distribution.
  # pmax(..., 0) prevents negative counts from rounding.
  init_N    <- pmax(round(stable_A * pop_size), 0L)
  init_ages <- rep(0:max_age, times = init_N)
  n_init    <- sum(init_N)

  # Randomly assign sex to each individual
  init_sex <- sample(c("F", "M"), n_init,
                     prob = c(female_fraction, 1 - female_fraction),
                     replace = TRUE)

  # Draw a personal maturity age for each individual from the ogive PMF.
  # Males and females may have different ogives.
  init_mat_age <- integer(n_init)
  is_f <- init_sex == "F"
  init_mat_age[is_f]  <- sample_mat_ages(ogive_f, sum(is_f))
  init_mat_age[!is_f] <- sample_mat_ages(ogive_m, sum(!is_f))

  # Assign permanent fertility status.  Each individual independently draws
  # a Bernoulli trial at birth; those who "fail" are permanently infertile.
  init_fertile <- rep(TRUE, n_init)
  if (infertility_f > 0) init_fertile[is_f]  <- runif(sum(is_f))  >= infertility_f
  if (infertility_m > 0) init_fertile[!is_f] <- runif(sum(!is_f)) >= infertility_m

  # Assign breeding state to mature, fertile females from the Markov stationary
  # distribution.  Immature females, infertile females, and all males get NA.
  init_breed_state <- rep(NA_integer_, n_init)
  mature_f <- which(is_f & init_ages >= init_mat_age & init_fertile)
  if (length(mature_f) > 0L) {
    init_breed_state[mature_f] <- sample(
      1:3, length(mature_f), replace = TRUE,
      prob = c(pi_1, pi_2, pi_3)
    )
  }

  # Assemble the population data.table.  Note: calculate.s0() uses a
  # simplified population structure (no id, birth_year, or parentage) because
  # father tracking is unnecessary for calibrating s0.
  pop <- data.table(
    age         = init_ages,
    sex         = init_sex,
    mat_age     = init_mat_age,
    fertile     = init_fertile,
    breed_state = init_breed_state
  )

  # ── Pod / superpod setup ──
  # Pods are fixed-size family groups; superpods are communities of pods.
  # Individuals are evenly distributed across pods, and pods are evenly
  # distributed across superpods.
  pod_to_sp <- NULL

  if (use_pods) {
    n_pods <- max(1L, round(n_init / pod_size))
    n_sp   <- max(1L, ceiling(n_pods / superpod_size))

    # Assign individuals to pods using round-robin, then look up each
    # pod's superpod from the pod_to_sp mapping.
    pod_vec   <- rep(seq_len(n_pods), length.out = n_init)
    pod_to_sp <- rep(seq_len(n_sp), each = superpod_size, length.out = n_pods)

    set(pop, j = "pod",      value = pod_vec)
    set(pop, j = "superpod", value = pod_to_sp[pod_vec])
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PHASE 3: ITERATIVE BISECTION SEARCH FOR s0
  # ═══════════════════════════════════════════════════════════════════════════
  # The Leslie-matrix s0 is deterministic — the stochastic IBM will typically
  # grow or shrink slightly at that s0.  We bisect around it: if the population
  # grows, s0 is too high (lower it); if it shrinks, s0 is too low (raise it).
  #
  # The simulation runs in "windows" of check_interval years.  After each
  # window, we measure the mean log-growth rate.  If |growth| < growth_tol
  # for stable_required consecutive windows, we declare convergence.

  s0_low     <- max(0.01, s0_leslie - 0.10)
  s0_high    <- min(0.99, s0_leslie + 0.10)
  s0_current <- s0_leslie

  year_counter       <- 0L
  consecutive_stable <- 0L

  # Parse sex-specific between-year stickiness (probability of staying in
  # the same superpod from one year to the next).
  if (!is.null(stickiness_year)) {
    if (length(stickiness_year) == 1L) {
      stick_yr_F <- stick_yr_M <- stickiness_year
    } else {
      stick_yr_F <- stickiness_year[1]
      stick_yr_M <- stickiness_year[2]
    }
  }

  message(sprintf(
    "Phase 2 — searching (check every %d yrs, tol = %.4f, need %d stable)...",
    check_interval, growth_tol, stable_required
  ))

  for (win in seq_len(max_windows)) {

    # Apply the current s0 candidate to the survival vector
    surv_vec    <- survival
    surv_vec[1] <- s0_current
    # Track population size at the end of each year within this window
    window_N    <- integer(check_interval)

    for (w in seq_len(check_interval)) {
      year_counter <- year_counter + 1L

      # ── Survival ──
      # Each individual survives with probability surv_vec[their age + 1].
      # Age is 0-indexed but R vectors are 1-indexed, hence the +1.
      rates <- surv_vec[pop$age + 1L]
      alive <- runif(nrow(pop)) <= rates
      pop   <- pop[alive]

      # ── Aging ──
      # Survivors age by one year.  Individuals that exceed max_age are removed
      # (they have lived their full lifespan).
      set(pop, j = "age", value = pop$age + 1L)
      pop <- pop[pop$age <= max_age]

      if (nrow(pop) == 0L) stop("Population went extinct during s0 search.")

      # ── Between-year superpod reshuffling ──
      # Each eligible individual (age >= weaning_age) independently decides
      # whether to stay or leave its current superpod.  Movers join a random
      # pod in a DIFFERENT superpod (emigration, not within-community shuffle).
      if (use_pods && !is.null(stickiness_year)) {
        elig <- which(pop$age >= wa)
        if (length(elig) > 0L) {
          elig_sex  <- pop$sex[elig]
          stay_prob <- ifelse(elig_sex == "F", stick_yr_F, stick_yr_M)
          movers    <- elig[runif(length(elig)) > stay_prob]

          if (length(movers) > 0L) {
            current_sp <- pop$superpod[movers]
            all_pods   <- unique(pop$pod)
            # Build a lookup: for each superpod, which pods are NOT in it?
            other_pool <- lapply(
              split(all_pods, pod_to_sp[all_pods]),
              function(x) all_pods[!all_pods %in% x]
            )
            new_pods <- integer(length(movers))
            for (sp in unique(current_sp)) {
              mask <- which(current_sp == sp)
              pool <- other_pool[[as.character(sp)]]
              # Fallback: if only one superpod exists, allow within-superpod moves
              if (is.null(pool) || length(pool) == 0L) pool <- all_pods
              new_pods[mask] <- sample(pool, length(mask), replace = TRUE)
            }
            set(pop, i = movers, j = "pod",      value = new_pods)
            set(pop, i = movers, j = "superpod", value = pod_to_sp[new_pods])
          }
        }
      }

      # ── Markov breeding state transitions ──
      # IMPORTANT: We identify all females in each state BEFORE applying any
      # transitions.  This prevents double-transitions (e.g., S2 → S3 → S1)
      # within a single year.

      # Newly mature females enter the cycle at S3 (resting = ready to conceive).
      # Only fertile females can enter the breeding cycle.
      new_mature <- which(pop$sex == "F" & pop$age == pop$mat_age &
                            is.na(pop$breed_state) & pop$fertile)
      if (length(new_mature) > 0L) {
        set(pop, i = new_mature, j = "breed_state", value = 3L)
      }

      # Snapshot current states before any transitions
      mother_idx <- which(pop$breed_state == 1L)  # S1: pregnant → will give birth
      s2_idx     <- which(pop$breed_state == 2L)  # S2: with calf → may conceive or rest
      s3_idx     <- which(pop$breed_state == 3L)  # S3: resting → may conceive or stay

      # S2 females: conceive again (→ S1) with prob ψ₂, or wean calf and rest (→ S3)
      if (length(s2_idx) > 0L) {
        new_state_s2 <- ifelse(runif(length(s2_idx)) < psi_2, 1L, 3L)
        set(pop, i = s2_idx, j = "breed_state", value = new_state_s2)
      }
      # S3 females: conceive (→ S1) with prob ψ₃, or remain resting (→ S3)
      if (length(s3_idx) > 0L) {
        new_state_s3 <- ifelse(runif(length(s3_idx)) < psi_3, 1L, 3L)
        set(pop, i = s3_idx, j = "breed_state", value = new_state_s3)
      }
      # S1 females: give birth (→ S2, now with dependent calf)
      if (length(mother_idx) > 0L) {
        set(pop, i = mother_idx, j = "breed_state", value = 2L)
      }

      # ── Breeding: create offspring ──
      # Only proceed if there are both pregnant females (S1) and at least one
      # mature, fertile male in the population.
      has_males <- any(pop$sex == "M" & pop$age >= pop$mat_age & pop$fertile)

      if (length(mother_idx) > 0L && has_males) {
        n_mothers <- length(mother_idx)
        # Litter size is Poisson-distributed with minimum 1 (1 + Pois(λ-1))
        litter_sizes <- 1L + rpois(n_mothers, lambda = litter_size - 1)
        n_yoy        <- sum(litter_sizes)

        # Assign sex to each newborn
        yoy_sex <- sample(c("F", "M"), n_yoy,
                          prob = c(female_fraction, 1 - female_fraction),
                          replace = TRUE)

        # Assign a personal maturity age from the sex-appropriate ogive
        yoy_mat_age <- integer(n_yoy)
        yoy_is_f    <- yoy_sex == "F"
        if (any(yoy_is_f))  yoy_mat_age[yoy_is_f]  <- sample_mat_ages(ogive_f, sum(yoy_is_f))
        if (any(!yoy_is_f)) yoy_mat_age[!yoy_is_f] <- sample_mat_ages(ogive_m, sum(!yoy_is_f))

        # Assign permanent fertility status
        yoy_fertile <- rep(TRUE, n_yoy)
        if (infertility_f > 0 && any(yoy_is_f))
          yoy_fertile[yoy_is_f]  <- runif(sum(yoy_is_f))  >= infertility_f
        if (infertility_m > 0 && any(!yoy_is_f))
          yoy_fertile[!yoy_is_f] <- runif(sum(!yoy_is_f)) >= infertility_m

        # Build newborn data.table (no parentage tracking in calculate.s0)
        yoy <- data.table(
          age         = rep(0L, n_yoy),
          sex         = yoy_sex,
          mat_age     = yoy_mat_age,
          fertile     = yoy_fertile,
          breed_state = NA_integer_
        )

        # Offspring inherit their mother's pod and superpod
        if (use_pods) {
          mother_pods <- pop$pod[mother_idx]
          yoy_pods    <- rep(mother_pods, times = litter_sizes)
          set(yoy, j = "pod",      value = yoy_pods)
          set(yoy, j = "superpod", value = pod_to_sp[yoy_pods])
        }

        pop <- rbindlist(list(pop, yoy), use.names = TRUE)
      }

      # Record population size at end of this year
      window_N[w] <- nrow(pop)

    } # end inner loop (one assessment window)

    # ── Assess growth rate for this window ──
    # Mean year-over-year log-growth across the window.
    growth  <- mean(diff(log(window_N)))
    s0_used <- s0_current

    if (abs(growth) < growth_tol) {
      # Growth is within tolerance — count as a stable window
      consecutive_stable <- consecutive_stable + 1L
      status <- sprintf("STABLE (%d/%d)", consecutive_stable, stable_required)
    } else {
      # Growth deviates — reset the stability counter and bisect s0.
      # If population grew, s0 was too high; if it shrank, s0 was too low.
      consecutive_stable <- 0L
      if (growth > 0) s0_high <- s0_current else s0_low <- s0_current
      s0_current <- (s0_low + s0_high) / 2
      status <- sprintf("-> new s0 = %.4f", s0_current)
    }

    message(sprintf(
      "  window %3d  (yr %d-%d)  s0=%.4f  growth=%+.5f  N=%s  %s",
      win,
      year_counter - check_interval + 1L, year_counter,
      s0_used, growth,
      format(window_N[check_interval], big.mark = ","),
      status
    ))

    # If we've seen enough consecutive stable windows, declare convergence
    if (consecutive_stable >= stable_required) break

  } # end bisection loop

  # ═══════════════════════════════════════════════════════════════════════════
  # RETURN RESULTS
  # ═══════════════════════════════════════════════════════════════════════════

  if (consecutive_stable < stable_required) {
    warning(sprintf(
      "Did not converge after %d windows (%d years). Returning best estimate.",
      max_windows, year_counter
    ))
  } else {
    message(sprintf(
      "\nConverged: s0 = %.4f  (%d years simulated, final N = %s)",
      s0_current, year_counter, format(nrow(pop), big.mark = ",")
    ))
  }

  # Build the final survival vector with the calibrated s0 in position 1
  survival_out    <- survival
  survival_out[1] <- s0_current

  # Return a bundled config list containing both the calibration results
  # and all the user-specified parameters.  This list is designed to be
  # passed directly to simulate.pop() as sim_config, so the user only
  # specifies life-history parameters once.
  invisible(list(
    # ── Calibration results ──
    s0              = s0_current,      # calibrated age-0 survival
    survival        = survival_out,    # full survival vector with s0 applied
    s0_leslie       = s0_leslie,       # analytical starting estimate
    final_N         = nrow(pop),       # population size at end of calibration
    years_simulated = year_counter,    # total years run during calibration
    # ── Shared life-history parameters (passed through to simulate.pop) ──
    max_age         = max_age,
    pop_size        = pop_size,
    maturity_age    = maturity_age,    # original user input (integer, vector, or list)
    litter_size     = litter_size,
    psi_2           = psi_2,
    psi_3           = psi_3,
    num_mates       = num_mates,
    female_fraction = female_fraction,
    infertility     = infertility,
    # ── Social structure ──
    pod_size        = pod_size,
    superpod_size   = superpod_size,
    stickiness_year = stickiness_year,
    male_behavior   = male_behavior,
    max_females     = max_females,
    weaning_age     = weaning_age
  ))
}
