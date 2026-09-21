#' Create a stable population configuration for simulation
#'
#' Determines the parameter adjustments needed to stabilize an age-structured
#' population, then returns a bundled configuration object for use by
#' \code{simulate.pop()}.  Two stabilization modes are available:
#' \itemize{
#'   \item \strong{s0 calibration} (\code{density_dependence = FALSE}, default):
#'     Bisects age-0 survival (s0) until population growth is near zero.  The
#'     user supplies a survival vector with a placeholder for s0, and this
#'     function finds the value that balances births and deaths.
#'   \item \strong{Density dependence} (\code{density_dependence = TRUE}):
#'     The user supplies a complete, realistic survival curve (including s0),
#'     and the function solves for the conception level (\code{theta}, the
#'     logit of \code{psi_rest}) that makes the population stable at carrying
#'     capacity K.  The conception gap between nursing and resting states is
#'     set by \code{rho} (log-odds ratio), which is invariant to the density-
#'     dependent shift.  In \code{simulate.pop()}, a Pella--Tomlinson
#'     compensation mechanism then adjusts conception rates dynamically based
#'     on depletion.
#' }
#'
#' Both modes use a two-phase approach:
#' \enumerate{
#'   \item \strong{Leslie matrix} -- analytically solves for the starting
#'     estimate (s0 or theta).
#'   \item \strong{Iterative simulation} -- refines via bisection until
#'     convergence.
#' }
#'
#' @param max_age Integer. Maximum age in the population.  Default 40
#'   (eastern spinner dolphin longevity).
#' @param survival Numeric vector of length \code{max_age + 1}. Annual survival
#'   probabilities for ages 0 through \code{max_age}.  When
#'   \code{density_dependence = FALSE}, the age-0 value is a placeholder that
#'   will be replaced by the calibrated estimate.  When
#'   \code{density_dependence = TRUE}, all values including age-0 are used as
#'   supplied.
#' @param pop_size Integer. Total starting population size. Equal to carrying
#'   capacity K when density dependence is active.
#' @param maturity_age Maturity specification. Default 9 (eastern spinner
#'   dolphin).  Can be:
#'   \itemize{
#'     \item An integer: knife-edged maturity at that age for both sexes.
#'     \item A numeric vector of length \code{max_age + 1}: cumulative
#'       probability of being mature at each age (ogive CDF), applied to both
#'       sexes.
#'     \item A list with \code{female} and \code{male} elements, each of which
#'       can be an integer or numeric vector as above.
#'   }
#' @param litter_size Numeric. Mean litter size per breeding female. Default 1.
#'
#' @param rho Numeric.  Log-odds gap between nursing and resting conception:
#'   \code{logit(psi_nurse) - logit(psi_rest)}.  \strong{Only used when
#'   \code{density_dependence = TRUE}}, where it is the sole user-supplied
#'   quantity governing the breeding cycle shape.  Must be negative (nursing
#'   suppresses conception).  Default \code{-3.255}, derived from the ~2.9\%
#'   of eastern spinner females found simultaneously pregnant and lactating
#'   (Chivers 1992).
#'
#'   \strong{Why rho and not a probability ratio:} Density-dependent
#'   compensation shifts both \code{psi_nurse} and \code{psi_rest} by the
#'   same amount on the logit scale.  Under that transform, the log-odds gap
#'   (\code{rho}) is exactly invariant -- it means the same thing at every
#'   depletion level, including K.  A probability ratio (\code{psi_nurse /
#'   psi_rest}) is NOT invariant: it drifts as both probabilities shift
#'   toward the boundaries.
#'
#'   \strong{Not used when \code{density_dependence = FALSE}.}  In that mode,
#'   use \code{calving_interval} + \code{suppression_ratio} or direct
#'   \code{psi_nurse}/\code{psi_rest} instead.
#'
#' @param calving_interval Numeric or NULL. Target calving interval in years
#'   (average time between successive births).  \strong{Only used when
#'   \code{density_dependence = FALSE}.}  When supplied, the function solves
#'   for the underlying conception probabilities (\code{psi_nurse},
#'   \code{psi_rest}) using \code{suppression_ratio} via \code{uniroot()} on
#'   the Markov breeding stationary distribution.  Must be >= 2 (structural
#'   floor: 1 year pregnancy + 1 year dependency).
#'
#'   \strong{Cannot be used with \code{density_dependence = TRUE}.}  In DD
#'   mode, the calving interval at K is emergent from the equilibrium solve
#'   (the replacement constraint R0 = 1 pins the absolute conception level,
#'   and \code{rho} pins the shape).  To target an empirically observed
#'   calving interval, use \code{target_interval} + \code{target_depletion},
#'   which anchor compensation strength (\code{dd_max}) at the depletion
#'   where the interval was measured.
#' @param suppression_ratio Numeric in (0, 1\]. Ratio of conception probability
#'   while nursing to conception probability while resting
#'   (\code{psi_nurse / psi_rest}).  Controls the strength of lactational
#'   suppression.  Default 0.15.  \strong{Only used when
#'   \code{density_dependence = FALSE} and \code{calving_interval} is
#'   supplied.}
#' @param psi_nurse Numeric in \[0, 1\] or NULL.  Conception probability while
#'   nursing a dependent calf (lactational suppression).  Supply together with
#'   \code{psi_rest} as an alternative to \code{calving_interval}.
#'   \strong{Only used when \code{density_dependence = FALSE}.}  Default NULL.
#' @param psi_rest Numeric in \[0, 1\] or NULL.  Conception probability while
#'   resting (no dependent calf), or after calf death.  Should be >=
#'   \code{psi_nurse}.  Must be supplied together with \code{psi_nurse}.
#'   \strong{Only used when \code{density_dependence = FALSE}.}  Default NULL.
#'
#' @param num_mates Integer. Number of mates per female per breeding cycle.
#'   Passed through to \code{simulate.pop()}.
#' @param female_fraction Numeric in (0, 1). Fraction of offspring that are
#'   female.  Default 0.5.
#' @param infertility Numeric (scalar or length-2 vector \code{c(female,
#'   male)}). Proportion of permanently infertile individuals.  Default 0.
#' @param pod_size Integer or NULL. Number of individuals per pod.
#' @param superpod_size Integer or NULL. Number of pods per superpod.
#' @param stickiness_year Numeric (scalar or length-2 vector
#'   \code{c(female, male)}) or NULL.  Between-year superpod fidelity, i.e.,
#'   probability that individuals remain in the same superpod between
#'   simulation years.
#' @param male_behavior Character or NULL. \code{"random"} or
#'   \code{"strong_bull"}.
#' @param max_females Integer or NULL. Per-year cap on matings per male.
#' @param weaning_age Integer or NULL. Age at independence from mother.  Also
#'   determines the maximum number of years a mother stays in the "with
#'   dependent calf" breeding state (S2).  Default 2 (calves dependent at
#'   ages 0 and 1, matching eastern spinner dolphin lactation of ~19 months).
#'   If NULL, defaults to 1-year dependency for the breeding cycle (no
#'   cow-calf social following).
#' @param density_dependence Logical. If \code{TRUE}, use Pella--Tomlinson
#'   density dependence instead of s0 calibration.  Default \code{FALSE}.
#' @param z_pt Numeric. Pella--Tomlinson shape parameter.  Default 2.39
#'   (IWC convention; MNPL at ~0.6 K).
#' @param dd_max Numeric or NULL. Maximum logit-scale shift on conception
#'   probabilities at zero depletion.  Controls compensation strength.
#'   \strong{In most cases, leave NULL and use \code{target_interval} /
#'   \code{target_depletion} instead.}  Supplying \code{dd_max} directly
#'   bypasses the target-interval anchor and is intended for sensitivity
#'   analyses or theoretical runs.
#' @param target_interval Numeric. Target calving interval in years at
#'   the reference depletion \code{target_depletion}.  Default 2.84 (eastern
#'   spinner pregnancy rate reciprocal, Chivers 1992).  Used to solve
#'   \code{dd_max} when \code{density_dependence = TRUE} and \code{dd_max}
#'   is NULL.  Must be >= 2 and shorter than the calving interval at K.
#' @param target_depletion Numeric in (0, 1). Reference depletion level at
#'   which \code{target_interval} was observed.  Default 0.3 (30\% of K).
#' @param check_interval Integer. Years per assessment window. Default 10.
#' @param growth_tol Numeric. Stability tolerance. Default 0.001.
#' @param stable_required Integer. Consecutive stable windows needed. Default 5.
#' @param max_windows Integer. Safety cap on windows. Default 100.
#'
#' @details
#'
#' ## Calibration Modes
#'
#' **s0 calibration** (\code{density_dependence = FALSE}, default): Searches for
#' the age-0 survival rate (s0) that results in a stable population.  The user
#' provides a survival vector with any placeholder for age 0; this function
#' replaces it with the calibrated estimate.
#'
#' **Density dependence mode** (\code{density_dependence = TRUE}): The user
#' supplies a complete survival curve (including s0) and \code{rho} (the
#' log-odds gap between nursing and resting conception).  The function solves
#' for \code{theta} (logit of psi_rest at K) such that R0 = 1 (replacement),
#' then derives \code{psi_nurse_K} and \code{psi_rest_K}.  The calving interval
#' at K is \strong{emergent} -- it falls out from the equilibrium solve and
#' \code{rho}, and is NOT directly settable.
#'
#' To target an empirically measured calving interval (e.g., 2.84 yr from
#' Chivers 1992), use \code{target_interval} + \code{target_depletion}.
#' These anchor \code{dd_max} (compensation strength) so the model reproduces
#' the target interval at the stated depletion -- the depletion at which the
#' data was actually collected, NOT at K.
#'
#' ## Breeding Cycle Parameterization
#'
#' \strong{DD = TRUE:}
#' \itemize{
#'   \item Supply \code{rho} (default \code{-3.255}).  This is the sole
#'     user-supplied breeding-cycle quantity.  \code{theta} (the absolute
#'     level) is solved for replacement at K.
#'   \item \code{psi_nurse = inv_logit(theta + rho)},
#'     \code{psi_rest = inv_logit(theta)}.
#'   \item Calving interval at K is emergent.
#'   \item To target an observed interval, use \code{target_interval}.
#' }
#'
#' \strong{DD = FALSE:}
#' \itemize{
#'   \item Path A: supply \code{calving_interval} + \code{suppression_ratio}.
#'     Solves for \code{psi_nurse}/\code{psi_rest} via \code{uniroot()}.
#'   \item Path B: supply \code{psi_nurse} + \code{psi_rest} directly.
#' }
#'
#' ## Social Structure Initialization
#'
#' Pods and superpods are created deterministically during initialization. If
#' \code{pod_size = 20} and \code{superpod_size = 10}, then pods 1--10 belong
#' to superpod 1, pods 11--20 to superpod 2, etc.
#'
#' ## Output Structure
#'
#' The returned list bundles all parameters needed by \code{simulate.pop()},
#' including:
#' \itemize{
#'   \item \strong{Calibration results}: \code{s0} (or \code{theta}),
#'     \code{final_N}, \code{years_simulated}.
#'   \item \strong{Life history}: \code{max_age}, \code{survival},
#'     \code{litter_size}.
#'   \item \strong{Breeding}: \code{psi_nurse}, \code{psi_rest},
#'     \code{num_mates}, \code{weaning_age}.
#'   \item \strong{DD (if active)}: \code{theta}, \code{rho},
#'     \code{psi_nurse_K}, \code{psi_rest_K}, \code{z_pt}, \code{dd_max},
#'     \code{K_1plus}, \code{target_interval}, \code{target_depletion}.
#' }
#'
#' @return A list (returned invisibly) with calibration results and all shared
#'   parameters needed by \code{simulate.pop()}.
#'
#' @references
#' Caswell, H. (2001). Matrix Population Models (2nd ed.). Sinauer Associates.
#'
#' Pella, J. J., & Tomlinson, P. K. (1973). A generalized stock production
#' model. IATTC Bulletin, 13, 422-458.
#'
#' Barlow, J., & Boveng, P. (1991). Modeling age-specific mortality for marine
#' mammal populations. Marine Mammal Science, 7(1), 50-65.
#'
#' @importFrom data.table data.table set rbindlist
#' @importFrom stats runif rpois uniroot qlogis plogis
#' @export
create.stable.pop <- function(max_age         = 40L,
                              survival,
                              pop_size,
                              maturity_age    = 9L,
                              litter_size     = 1,
                              # -- Breeding cycle: DD=TRUE uses rho --
                              rho             = -3.255,
                              # -- Breeding cycle: DD=FALSE uses these --
                              calving_interval  = NULL,
                              suppression_ratio = 0.15,
                              psi_nurse       = NULL,
                              psi_rest        = NULL,
                              # -- Demographics --
                              num_mates       = 1L,
                              female_fraction = 0.5,
                              infertility     = 0,
                              # -- Social structure --
                              pod_size        = NULL,
                              superpod_size   = NULL,
                              stickiness_year = NULL,
                              male_behavior   = "random",
                              max_females     = NULL,
                              weaning_age     = 2L,
                              # -- Density dependence --
                              density_dependence = FALSE,
                              z_pt            = 2.39,
                              dd_max          = NULL,
                              target_interval  = 2.84,
                              target_depletion = 0.3,
                              # -- Convergence --
                              check_interval  = 10L,
                              growth_tol      = 0.001,
                              stable_required = 5L,
                              max_windows     = 100L) {

  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================

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

  # -- Breeding cycle parameterization --
  # DD=TRUE: rho (log-odds gap) is the sole breeding input; theta is solved.
  # DD=FALSE: calving_interval + suppression_ratio, or psi_nurse + psi_rest.
  has_ci  <- !is.null(calving_interval)
  has_psi <- !is.null(psi_nurse) || !is.null(psi_rest)

  if (density_dependence) {
    # DD mode: rho only. Error on calving_interval or psi_nurse/psi_rest.
    if (has_ci)
      stop("Cannot use `calving_interval` with `density_dependence = TRUE`. ",
           "In DD mode, the calving interval at K is emergent from the ",
           "equilibrium solve and `rho`. To target an observed calving ",
           "interval, use `target_interval` + `target_depletion`.")
    if (has_psi)
      stop("Cannot use `psi_nurse`/`psi_rest` with `density_dependence = TRUE`. ",
           "In DD mode, supply `rho` (log-odds gap) instead; psi values are ",
           "solved for the replacement condition R0 = 1.")
    if (!is.numeric(rho) || length(rho) != 1L)
      stop("`rho` must be a single numeric value.")
    if (rho >= 0)
      warning("`rho` is >= 0, meaning conception while nursing is at least as ",
              "likely as while resting. Lactational suppression requires rho < 0.")
  } else {
    # DD=FALSE: calving_interval or psi_nurse/psi_rest
    if (has_ci && has_psi)
      stop("Supply `calving_interval` or `psi_nurse`/`psi_rest`, not both.")
    if (has_psi && (is.null(psi_nurse) || is.null(psi_rest)))
      stop("Both `psi_nurse` and `psi_rest` must be supplied together.")
    if (!has_ci && !has_psi)
      stop("When `density_dependence = FALSE`, supply either ",
           "`calving_interval` (recommended) or both `psi_nurse` and ",
           "`psi_rest` to define the breeding cycle.")

    if (has_ci) {
      if (!is.numeric(calving_interval) || length(calving_interval) != 1L ||
          calving_interval <= 0)
        stop("`calving_interval` must be a positive number.")
      if (calving_interval < 2)
        stop("`calving_interval` must be >= 2 (structural floor: 1 yr ",
             "pregnancy + 1 yr dependency).")
      if (!is.numeric(suppression_ratio) || length(suppression_ratio) != 1L ||
          suppression_ratio <= 0 || suppression_ratio > 1)
        stop("`suppression_ratio` must be in (0, 1].")
    }
    if (has_psi) {
      if (!is.numeric(psi_nurse) || length(psi_nurse) != 1L ||
          psi_nurse < 0 || psi_nurse > 1)
        stop("`psi_nurse` must be a single number between 0 and 1.")
      if (!is.numeric(psi_rest) || length(psi_rest) != 1L ||
          psi_rest < 0 || psi_rest > 1)
        stop("`psi_rest` must be a single number between 0 and 1.")
      if (psi_nurse > psi_rest)
        warning("`psi_nurse` > `psi_rest`: lactational suppression is reversed.")
    }
  }

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
        weaning_age < 1 || weaning_age != round(weaning_age))
      stop("`weaning_age` must be a positive integer (>= 1).")
    if (weaning_age > max_age)
      stop("`weaning_age` (", weaning_age, ") must be <= max_age (", max_age, ").")
  }

  # -- Density dependence parameters --
  if (!is.logical(density_dependence) || length(density_dependence) != 1L)
    stop("`density_dependence` must be TRUE or FALSE.")

  if (density_dependence) {
    # dd_max takes precedence if supplied; otherwise solve from target_interval
    has_dd <- !is.null(dd_max)

    if (has_dd) {
      if (!is.numeric(dd_max) || length(dd_max) != 1L || dd_max <= 0)
        stop("`dd_max` must be a positive number.")
      message(
        "Note: `dd_max` supplied directly. In most cases, leave `dd_max` ",
        "NULL and use `target_interval` / `target_depletion` instead, ",
        "which anchors compensation strength to an observed calving interval."
      )
    } else {
      # Solve dd_max from target_interval (has a default of 2.84)
      if (!is.numeric(target_interval) || length(target_interval) != 1L ||
          target_interval <= 0)
        stop("`target_interval` must be a positive number.")
      if (target_interval < 2)
        stop("`target_interval` must be >= 2 (structural floor: 1 yr ",
             "pregnancy + 1 yr dependency).")
      if (!is.numeric(target_depletion) || length(target_depletion) != 1L ||
          target_depletion <= 0 || target_depletion >= 1)
        stop("`target_depletion` must be in (0, 1), exclusive.")
    }

    if (!is.numeric(z_pt) || length(z_pt) != 1L || z_pt <= 0)
      stop("`z_pt` must be a positive number.")
  }

  # -- Convergence parameters --
  if (check_interval < 1L) stop("`check_interval` must be >= 1.")
  if (growth_tol <= 0)     stop("`growth_tol` must be > 0.")
  if (stable_required < 1) stop("`stable_required` must be >= 1.")
  if (max_windows < 1)     stop("`max_windows` must be >= 1.")

  # ===========================================================================
  # INTERNAL HELPERS
  # ===========================================================================

  mat_lambda <- function(A) Mod(eigen(A, only.values = TRUE)$values[1])

  mat_stable <- function(A) {
    w <- Mod(eigen(A)$vectors[, 1])
    w / sum(w)
  }

  # -- Breeding cycle stationary distribution --
  # Builds the Markov transition matrix for the generalized breeding cycle
  # and solves for the stationary distribution.
  #
  # States: S1 (pregnant), S2(1)..S2(w) (with dependent calf, year 1..w), S3 (resting)
  # Dimension: w + 2
  #
  # Calf survival couples the mother's transition to her calf's fate:
  #   - Calf alive -> conception prob = psi_n (suppressed)
  #   - Calf dead  -> conception prob = psi_r (released)
  # The population-level effective rate at S2 year k is the mixture:
  #   psi_k = ell_k * psi_n + (1 - ell_k) * psi_r
  # where ell_k = sv[k] is the calf's survival probability for that year.
  breeding_stationary <- function(psi_n, psi_r, sv, w) {
    n_st <- w + 2L
    P <- matrix(0, nrow = n_st, ncol = n_st)

    # S1 -> S2(1): give birth (deterministic)
    P[1, 2] <- 1

    # S2(k) transitions
    for (k in seq_len(w)) {
      row <- k + 1L
      ell_k <- sv[k]
      psi_k <- ell_k * psi_n + (1 - ell_k) * psi_r

      P[row, 1] <- psi_k  # -> S1 (conceive)

      if (k < w) {
        # Calf alive AND no conception -> advance to S2(k+1)
        P[row, row + 1] <- ell_k * (1 - psi_n)
        # Calf dead AND no conception -> S3
        P[row, n_st] <- (1 - ell_k) * (1 - psi_r)
      } else {
        # Last dependent year: calf weaned regardless -> S3 if no conception
        P[row, n_st] <- 1 - psi_k
      }
    }

    # S3 -> S1 or stay S3
    P[n_st, 1]    <- psi_r
    P[n_st, n_st] <- 1 - psi_r

    # Stationary distribution: left eigenvector of P
    ev  <- eigen(t(P))
    idx <- which.min(abs(Mod(ev$values) - 1))
    pi  <- Mod(ev$vectors[, idx])
    pi / sum(pi)
  }

  # Solve for dd_max that reproduces a target calving interval at a reference

  # depletion level.  Called after theta_shift calibration (which determines
  # psi_nurse_K and psi_rest_K).  The calving interval = 1 / pi_1, where pi_1
  # is the stationary probability of S1 (pregnant) in the breeding chain with
  # conception rates shifted by dd_max * (1 - D^z) on the logit scale.
  solve_dd_max <- function(psi_n_K, psi_r_K, target_int, target_dep, z, sv, w) {
    logit_n    <- qlogis(psi_n_K)
    logit_r    <- qlogis(psi_r_K)
    shift_frac <- 1 - target_dep^z  # constant for this D and z

    # Calving interval at K (dd_max contribution = 0)
    pi_1_K      <- breeding_stationary(psi_n_K, psi_r_K, sv, w)[1]
    interval_K  <- 1 / pi_1_K

    if (target_int >= interval_K)
      stop(sprintf(
        paste0("`target_interval` (%.2f yr) must be shorter than the calving ",
               "interval at K (%.2f yr). Compensation can only shorten the ",
               "interval below K, not lengthen it."),
        target_int, interval_K
      ))

    obj <- function(dd) {
      delta <- dd * shift_frac
      pn_d  <- plogis(logit_n + delta)
      pr_d  <- plogis(logit_r + delta)
      pi_1  <- breeding_stationary(pn_d, pr_d, sv, w)[1]
      1 / pi_1 - target_int
    }

    # Lower bound: dd_max just above 0 (interval ≈ at-K value)
    # Upper bound: dd_max = 50 (extremely strong compensation)
    result <- tryCatch(
      uniroot(obj, interval = c(1e-4, 50), tol = 1e-6),
      error = function(e) {
        stop(sprintf(
          paste0("Cannot find dd_max that produces a %.2f-yr calving interval ",
                 "at D=%.2f. The target may be below the structural floor ",
                 "(2 yr). Error: %s"),
          target_int, target_dep, conditionMessage(e)
        ))
      }
    )
    result$root
  }

  # -- Solve for psi_nurse/psi_rest from a target calving interval --
  # Given a target calving interval (years) and a suppression ratio

  # (psi_nurse / psi_rest), solves for psi_rest such that the stationary
  # calving interval from the Markov breeding chain equals the target.
  # psi_nurse is then derived as suppression_ratio * psi_rest.
  #
  # The interval is monotonically decreasing in psi_rest (higher conception
  # rate -> shorter interval), so uniroot is well-posed.
  solve_psi_from_interval <- function(calving_int, sup_ratio, sv, w) {
    # Upper bound: both psi_rest and psi_nurse = sup_ratio * psi_rest must be < 1
    upper <- min(1 - 1e-6, (1 - 1e-6) / sup_ratio)

    obj <- function(psi_r) {
      psi_n <- sup_ratio * psi_r
      pi_1  <- breeding_stationary(psi_n, psi_r, sv, w)[1]
      1 / pi_1 - calving_int
    }

    # Check that the target interval is achievable
    # At psi_rest -> upper: shortest achievable interval
    shortest <- 1 / breeding_stationary(sup_ratio * upper, upper, sv, w)[1]
    if (calving_int <= shortest)
      stop(sprintf(
        paste0("`calving_interval` (%.2f yr) is shorter than the minimum ",
               "achievable interval (%.2f yr) for this survival curve and ",
               "suppression_ratio (%.2f). Try a smaller suppression_ratio or ",
               "a longer calving_interval."),
        calving_int, shortest, sup_ratio
      ))

    result <- tryCatch(
      uniroot(obj, interval = c(1e-6, upper), tol = 1e-6),
      error = function(e) {
        stop(sprintf(
          paste0("Cannot solve for psi_rest that produces a %.2f-yr calving ",
                 "interval with suppression_ratio=%.2f. Error: %s"),
          calving_int, sup_ratio, conditionMessage(e)
        ))
      }
    )
    psi_r <- result$root
    list(psi_nurse = sup_ratio * psi_r, psi_rest = psi_r)
  }

  # ===========================================================================
  # PARSE MATURITY SPECIFICATION
  # ===========================================================================

  make_ogive <- function(x, max_a) {
    n <- max_a + 1L
    if (length(x) == 1L && x == round(x)) {
      ogive <- rep(0, n)
      if (x <= max_a) ogive[seq(x + 1L, n)] <- 1
      return(ogive)
    }
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

  sample_mat_ages <- function(ogive, n) {
    pmf     <- diff(c(0, ogive))
    p_never <- 1 - sum(pmf)
    if (p_never > 0.001) {
      pmf  <- c(pmf, p_never)
      ages <- c(0:max_age, max_age + 1L)
    } else {
      ages <- 0:max_age
      pmf[length(pmf)] <- pmf[length(pmf)] + p_never
    }
    sample(ages, n, replace = TRUE, prob = pmf)
  }

  # ===========================================================================
  # PARSE INFERTILITY
  # ===========================================================================

  if (length(infertility) == 1L) {
    infertility_f <- infertility_m <- infertility
  } else {
    infertility_f <- infertility[1]
    infertility_m <- infertility[2]
  }

  # ===========================================================================
  # PARSE WEANING AGE FOR BREEDING CYCLE
  # ===========================================================================
  # weaning_age serves two roles:
  #   1. Breeding cycle: how many years a mother stays in S2 (with calf)
  #   2. Social structure: calves below weaning_age follow mother's pod
  # If NULL, the breeding cycle defaults to 1-year dependency, and there
  # is no cow-calf social following.

  wa_breed <- if (is.null(weaning_age)) 1L else as.integer(weaning_age)

  # ===========================================================================
  # RESOLVE BREEDING CYCLE PARAMETERIZATION
  # ===========================================================================
  # DD=TRUE:  rho is the sole breeding input; theta (logit psi_rest) is solved
  #           for R0 = 1 at K. psi_nurse/psi_rest are derived, not supplied.
  # DD=FALSE: calving_interval + suppression_ratio -> solve psi values
  #           OR psi_nurse + psi_rest directly.

  psi_from_interval <- FALSE

  if (!density_dependence) {
    # DD=FALSE: resolve psi values from calving_interval or direct supply
    if (has_ci) {
      psi_solved <- solve_psi_from_interval(calving_interval, suppression_ratio,
                                            survival, wa_breed)
      psi_nurse <- psi_solved$psi_nurse
      psi_rest  <- psi_solved$psi_rest
      psi_from_interval <- TRUE
      message(sprintf(
        "Solved from calving_interval=%.2f yr, suppression_ratio=%.2f: psi_nurse=%.4f, psi_rest=%.4f",
        calving_interval, suppression_ratio, psi_nurse, psi_rest
      ))
    }
    # If has_psi, psi_nurse/psi_rest are already set from user input.
  }
  # DD=TRUE: psi values are derived in Phase 1 below (from rho + theta solve).
  # No psi resolution needed here.

  # ===========================================================================
  # PHASE 1: LESLIE MATRIX -- STARTING ESTIMATE
  # ===========================================================================
  # Shift the female ogive right by one year: newly mature females enter at S3
  # (resting) and cannot breed in their first year of maturity.
  ogive_f_leslie <- c(0, ogive_f[seq_len(max_age)])
  n_classes <- max_age + 1L

  if (!density_dependence) {
    # -- Mode A: solve for s0 --
    # Build Leslie matrix with current psi_nurse/psi_rest and find the s0 that
    # makes the dominant eigenvalue lambda = 1. The fecundity row depends on pi_1
    # (breeding fraction), which itself depends on s0 through calf survival.

    A <- matrix(0, nrow = n_classes, ncol = n_classes)
    for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i]

    leslie_fn <- function(s0) {
      sv    <- survival
      sv[1] <- s0
      A[2, 1] <- s0
      pi_1 <- breeding_stationary(psi_nurse, psi_rest, sv, wa_breed)[1]
      ff   <- litter_size * female_fraction * pi_1 * (1 - infertility_f)
      A[1, ] <- ogive_f_leslie * ff
      mat_lambda(A) - 1
    }

    # Check that a valid s0 exists
    if (leslie_fn(0.01) * leslie_fn(0.99) > 0) {
      stop(
        "No s0 in [0.01, 0.99] can stabilise this population. ",
        sprintf("Lambda-1 = %.4f at s0=0.01 and %.4f at s0=0.99.",
                leslie_fn(0.01), leslie_fn(0.99))
      )
    }

    s0_leslie <- uniroot(leslie_fn, interval = c(0.01, 0.99))$root

    # Rebuild Leslie at the solved s0 for stable age distribution
    sv_leslie      <- survival
    sv_leslie[1]   <- s0_leslie
    A[2, 1]        <- s0_leslie
    pi_stat        <- breeding_stationary(psi_nurse, psi_rest, sv_leslie, wa_breed)
    pi_1           <- pi_stat[1]
    ff             <- litter_size * female_fraction * pi_1 * (1 - infertility_f)
    A[1, ]         <- ogive_f_leslie * ff
    stable_A       <- mat_stable(A)

    # Effective psi values for initialisation (unchanged from user input)
    psi_nurse_eff <- psi_nurse
    psi_rest_eff  <- psi_rest

    psi_label <- if (psi_from_interval) {
      sprintf("from calving_interval=%.2f", calving_interval)
    } else {
      "direct"
    }
    message(sprintf(
      "Breeding cycle (%s): psi_nurse=%.3f, psi_rest=%.3f -> %.1f%% breeding/yr, avg interval=%.1f yrs",
      psi_label, psi_nurse, psi_rest, pi_1 * 100, 1 / pi_1
    ))
    message(sprintf(
      "Phase 1 -- Leslie estimate: s0 = %.4f  (lambda = %.6f)",
      s0_leslie, mat_lambda(A)
    ))

  } else {
    # -- Mode B: solve for theta (density dependence) --
    # Survival is fixed (user-supplied, including s0). Find theta
    # (logit of psi_rest at K) that makes lambda(K) = 1.
    # rho (log-odds gap) is the sole user-supplied breeding quantity.

    A <- matrix(0, nrow = n_classes, ncol = n_classes)
    for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i]

    dd_leslie_fn <- function(theta) {
      pn <- plogis(theta + rho)
      pr <- plogis(theta)
      pi_1 <- breeding_stationary(pn, pr, survival, wa_breed)[1]
      ff   <- litter_size * female_fraction * pi_1 * (1 - infertility_f)
      A[1, ] <- ogive_f_leslie * ff
      mat_lambda(A) - 1
    }

    # Check feasibility
    if (dd_leslie_fn(-10) * dd_leslie_fn(10) > 0) {
      stop(
        "No theta in [-10, 10] can stabilise this population. ",
        sprintf("Lambda-1 = %.4f at theta=-10 and %.4f at theta=+10. ",
                dd_leslie_fn(-10), dd_leslie_fn(10)),
        "Check that the survival curve and `rho` are plausible."
      )
    }

    theta_leslie <- uniroot(dd_leslie_fn, interval = c(-10, 10))$root

    # Compute at-K conception rates and stationary distribution
    psi_nurse_eff <- plogis(theta_leslie + rho)
    psi_rest_eff  <- plogis(theta_leslie)

    # Set psi_nurse/psi_rest for downstream use (init, output)
    psi_nurse <- psi_nurse_eff
    psi_rest  <- psi_rest_eff

    pi_stat <- breeding_stationary(psi_nurse_eff, psi_rest_eff, survival, wa_breed)
    pi_1    <- pi_stat[1]

    # Rebuild Leslie with solved fecundity for stable age distribution
    ff     <- litter_size * female_fraction * pi_1 * (1 - infertility_f)
    A[1, ] <- ogive_f_leslie * ff
    stable_A <- mat_stable(A)

    # K_1plus: 1+ component of carrying capacity (excludes age-0)
    K_1plus <- round(pop_size * (1 - stable_A[1]))

    message(sprintf(
      "Breeding cycle (rho=%.3f): psi_nurse_K=%.4f, psi_rest_K=%.4f",
      rho, psi_nurse_eff, psi_rest_eff
    ))
    message(sprintf(
      "Phase 1 -- Leslie estimate: theta = %.4f  (lambda=%.6f)",
      theta_leslie, mat_lambda(A)
    ))
    message(sprintf(
      "  K_1plus = %s  (%.1f%% breeding/yr at K, avg interval=%.1f yrs)",
      format(K_1plus, big.mark = ","), pi_1 * 100, 1 / pi_1
    ))
  }

  # ===========================================================================
  # INITIALISE THE INDIVIDUAL-BASED POPULATION
  # ===========================================================================

  use_pods <- !is.null(pod_size)
  wa       <- if (is.null(weaning_age)) 0L else weaning_age

  init_N    <- pmax(round(stable_A * pop_size), 0L)
  init_ages <- rep(0:max_age, times = init_N)
  n_init    <- sum(init_N)

  init_sex <- sample(c("F", "M"), n_init,
                     prob = c(female_fraction, 1 - female_fraction),
                     replace = TRUE)

  init_mat_age <- integer(n_init)
  is_f <- init_sex == "F"
  init_mat_age[is_f]  <- sample_mat_ages(ogive_f, sum(is_f))
  init_mat_age[!is_f] <- sample_mat_ages(ogive_m, sum(!is_f))

  init_fertile <- rep(TRUE, n_init)
  if (infertility_f > 0) init_fertile[is_f]  <- runif(sum(is_f))  >= infertility_f
  if (infertility_m > 0) init_fertile[!is_f] <- runif(sum(!is_f)) >= infertility_m

  # Assign breeding states from the full stationary distribution.
  # pi_stat has wa_breed + 2 elements: S1, S2(1), ..., S2(wa_breed), S3
  # We use the effective psi values (adjusted if DD mode).
  surv_init <- if (!density_dependence) {
    sv <- survival; sv[1] <- s0_leslie; sv
  } else {
    survival
  }
  pi_stat_init <- breeding_stationary(psi_nurse_eff, psi_rest_eff,
                                       surv_init, wa_breed)
  n_breed_states <- length(pi_stat_init)

  init_breed_state <- rep(NA_integer_, n_init)
  init_s2_year     <- rep(NA_integer_, n_init)

  mature_f <- which(is_f & init_ages >= init_mat_age & init_fertile)
  if (length(mature_f) > 0L) {
    state_idx <- sample(seq_len(n_breed_states), length(mature_f),
                        replace = TRUE, prob = pi_stat_init)
    # Map to breed_state: 1=S1, 2=S2, 3=S3
    bs <- ifelse(state_idx == 1L, 1L,
                 ifelse(state_idx == n_breed_states, 3L, 2L))
    # s2_year = sub-state index for S2 mothers
    s2y <- ifelse(bs == 2L, state_idx - 1L, NA_integer_)
    init_breed_state[mature_f] <- bs
    init_s2_year[mature_f]     <- s2y
  }

  pop <- data.table(
    age         = init_ages,
    sex         = init_sex,
    mat_age     = init_mat_age,
    fertile     = init_fertile,
    breed_state = init_breed_state,
    s2_year     = init_s2_year
  )

  # -- Pod / superpod setup --
  pod_to_sp <- NULL
  if (use_pods) {
    n_pods <- max(1L, round(n_init / pod_size))
    n_sp   <- max(1L, ceiling(n_pods / superpod_size))

    pod_vec   <- rep(seq_len(n_pods), length.out = n_init)
    pod_to_sp <- rep(seq_len(n_sp), each = superpod_size, length.out = n_pods)

    set(pop, j = "pod",      value = pod_vec)
    set(pop, j = "superpod", value = pod_to_sp[pod_vec])
  }

  # ===========================================================================
  # PHASE 2: ITERATIVE BISECTION
  # ===========================================================================
  # Mode A (DD=FALSE): bisect s0; psi values are constant.
  # Mode B (DD=TRUE):  bisect theta; survival is constant.

  surv_vec <- survival

  if (!density_dependence) {
    s0_low     <- max(0.01, s0_leslie - 0.10)
    s0_high    <- min(0.99, s0_leslie + 0.10)
    s0_current <- s0_leslie
  } else {
    theta_low     <- theta_leslie - 2
    theta_high    <- theta_leslie + 2
    theta_current <- theta_leslie
  }

  year_counter       <- 0L
  consecutive_stable <- 0L

  if (!is.null(stickiness_year)) {
    if (length(stickiness_year) == 1L) {
      stick_yr_F <- stick_yr_M <- stickiness_year
    } else {
      stick_yr_F <- stickiness_year[1]
      stick_yr_M <- stickiness_year[2]
    }
  }

  message(sprintf(
    "Phase 2 -- searching (check every %d yrs, tol = %.4f, need %d stable)...",
    check_interval, growth_tol, stable_required
  ))

  for (win in seq_len(max_windows)) {

    # Apply the current candidate to the relevant parameter
    if (!density_dependence) {
      surv_vec[1]  <- s0_current
      psi_nurse_w  <- psi_nurse
      psi_rest_w   <- psi_rest
    } else {
      psi_nurse_w <- plogis(theta_current + rho)
      psi_rest_w  <- plogis(theta_current)
    }

    window_N <- integer(check_interval)

    for (w in seq_len(check_interval)) {
      year_counter <- year_counter + 1L

      # -- Survival --
      rates <- surv_vec[pop$age + 1L]
      alive <- runif(nrow(pop)) <= rates
      pop   <- pop[alive]

      # -- Aging --
      set(pop, j = "age", value = pop$age + 1L)
      pop <- pop[pop$age <= max_age]

      if (nrow(pop) == 0L) stop("Population went extinct during calibration.")

      # -- Between-year superpod reshuffling --
      if (use_pods && !is.null(stickiness_year)) {
        elig <- which(pop$age >= wa)
        if (length(elig) > 0L) {
          elig_sex  <- pop$sex[elig]
          stay_prob <- ifelse(elig_sex == "F", stick_yr_F, stick_yr_M)
          movers    <- elig[runif(length(elig)) > stay_prob]

          if (length(movers) > 0L) {
            current_sp <- pop$superpod[movers]
            all_pods   <- unique(pop$pod)
            other_pool <- lapply(
              split(all_pods, pod_to_sp[all_pods]),
              function(x) all_pods[!all_pods %in% x]
            )
            new_pods <- integer(length(movers))
            for (sp in unique(current_sp)) {
              mask <- which(current_sp == sp)
              pool <- other_pool[[as.character(sp)]]
              if (is.null(pool) || length(pool) == 0L) pool <- all_pods
              # Note: sample(pool, n) would misbehave if pool has length 1 (R
              # reinterprets a length-1 numeric x as the range 1:x). Indexing
              # by position avoids this.
              new_pods[mask] <- pool[sample.int(length(pool), length(mask), replace = TRUE)]
            }
            set(pop, i = movers, j = "pod",      value = new_pods)
            set(pop, i = movers, j = "superpod", value = pod_to_sp[new_pods])
          }
        }
      }

      # -- Markov breeding state transitions --
      # Newly mature, fertile females enter at S3
      new_mature <- which(pop$sex == "F" & pop$age == pop$mat_age &
                            is.na(pop$breed_state) & pop$fertile)
      if (length(new_mature) > 0L) {
        set(pop, i = new_mature, j = "breed_state", value = 3L)
      }

      # Snapshot states before transitions (prevents double transitions)
      mother_idx <- which(pop$breed_state == 1L)   # S1: will give birth
      s2_idx     <- which(pop$breed_state == 2L)   # S2: with calf
      s3_idx     <- which(pop$breed_state == 3L)   # S3: resting

      # -- S2 transitions: calf-survival-dependent --
      # Each S2 mother's transition depends on whether her calf survived.
      # In the calibration we don't track individual calves, so we draw
      # calf survival stochastically from surv_vec based on s2_year.
      if (length(s2_idx) > 0L) {
        k <- pop$s2_year[s2_idx]
        calf_surv_prob <- surv_vec[k]
        calf_alive <- runif(length(s2_idx)) < calf_surv_prob

        # Conception probability depends on calf fate
        psi_eff <- ifelse(calf_alive, psi_nurse_w, psi_rest_w)
        conceive <- runif(length(s2_idx)) < psi_eff

        # Determine new state
        new_state <- rep(3L, length(s2_idx))      # default: S3 (resting)
        new_state[conceive] <- 1L                  # conceive -> S1
        stay_s2 <- !conceive & calf_alive & k < wa_breed
        new_state[stay_s2] <- 2L                   # stay S2

        set(pop, i = s2_idx, j = "breed_state", value = new_state)

        # Update s2_year: increment for stayers, clear for leavers
        new_s2y <- rep(NA_integer_, length(s2_idx))
        new_s2y[stay_s2] <- k[stay_s2] + 1L
        set(pop, i = s2_idx, j = "s2_year", value = new_s2y)
      }

      # -- S3 transitions --
      if (length(s3_idx) > 0L) {
        new_state_s3 <- ifelse(runif(length(s3_idx)) < psi_rest_w, 1L, 3L)
        set(pop, i = s3_idx, j = "breed_state", value = new_state_s3)
      }

      # -- S1 -> S2: give birth --
      if (length(mother_idx) > 0L) {
        set(pop, i = mother_idx, j = "breed_state", value = 2L)
        set(pop, i = mother_idx, j = "s2_year",     value = 1L)
      }

      # -- Breeding: create offspring --
      has_males <- any(pop$sex == "M" & pop$age >= pop$mat_age & pop$fertile)

      if (length(mother_idx) > 0L && has_males) {
        n_mothers <- length(mother_idx)
        litter_sizes <- 1L + rpois(n_mothers, lambda = litter_size - 1)
        n_yoy        <- sum(litter_sizes)

        yoy_sex <- sample(c("F", "M"), n_yoy,
                          prob = c(female_fraction, 1 - female_fraction),
                          replace = TRUE)

        yoy_mat_age <- integer(n_yoy)
        yoy_is_f    <- yoy_sex == "F"
        if (any(yoy_is_f))  yoy_mat_age[yoy_is_f]  <- sample_mat_ages(ogive_f, sum(yoy_is_f))
        if (any(!yoy_is_f)) yoy_mat_age[!yoy_is_f] <- sample_mat_ages(ogive_m, sum(!yoy_is_f))

        yoy_fertile <- rep(TRUE, n_yoy)
        if (infertility_f > 0 && any(yoy_is_f))
          yoy_fertile[yoy_is_f]  <- runif(sum(yoy_is_f))  >= infertility_f
        if (infertility_m > 0 && any(!yoy_is_f))
          yoy_fertile[!yoy_is_f] <- runif(sum(!yoy_is_f)) >= infertility_m

        yoy <- data.table(
          age         = rep(0L, n_yoy),
          sex         = yoy_sex,
          mat_age     = yoy_mat_age,
          fertile     = yoy_fertile,
          breed_state = NA_integer_,
          s2_year     = NA_integer_
        )

        if (use_pods) {
          mother_pods <- pop$pod[mother_idx]
          yoy_pods    <- rep(mother_pods, times = litter_sizes)
          set(yoy, j = "pod",      value = yoy_pods)
          set(yoy, j = "superpod", value = pod_to_sp[yoy_pods])
        }

        pop <- rbindlist(list(pop, yoy), use.names = TRUE)
      }

      window_N[w] <- nrow(pop)
    } # end inner loop

    # -- Assess growth rate for this window --
    growth <- mean(diff(log(window_N)))

    if (abs(growth) < growth_tol) {
      consecutive_stable <- consecutive_stable + 1L
      status <- sprintf("STABLE (%d/%d)", consecutive_stable, stable_required)
    } else {
      consecutive_stable <- 0L

      if (!density_dependence) {
        if (growth > 0) s0_high <- s0_current else s0_low <- s0_current
        s0_current <- (s0_low + s0_high) / 2
        status <- sprintf("-> new s0 = %.4f", s0_current)
      } else {
        if (growth > 0) theta_high <- theta_current else theta_low <- theta_current
        theta_current <- (theta_low + theta_high) / 2
        status <- sprintf("-> new theta = %.4f", theta_current)
      }
    }

    if (!density_dependence) {
      message(sprintf(
        "  window %3d  (yr %d-%d)  s0=%.4f  growth=%+.5f  N=%s  %s",
        win,
        year_counter - check_interval + 1L, year_counter,
        s0_current, growth,
        format(window_N[check_interval], big.mark = ","),
        status
      ))
    } else {
      message(sprintf(
        "  window %3d  (yr %d-%d)  theta=%+.4f  growth=%+.5f  N=%s  %s",
        win,
        year_counter - check_interval + 1L, year_counter,
        theta_current, growth,
        format(window_N[check_interval], big.mark = ","),
        status
      ))
    }

    if (consecutive_stable >= stable_required) break
  } # end bisection loop

  # ===========================================================================
  # RETURN RESULTS
  # ===========================================================================

  if (consecutive_stable < stable_required) {
    warning(sprintf(
      "Did not converge after %d windows (%d years). Returning best estimate.",
      max_windows, year_counter
    ))
  }

  if (!density_dependence) {
    # -- s0 calibration output --
    survival_out    <- survival
    survival_out[1] <- s0_current

    message(sprintf(
      "\nConverged: s0 = %.4f  (%d years simulated, final N = %s)",
      s0_current, year_counter, format(nrow(pop), big.mark = ",")
    ))

    invisible(list(
      # Calibration results
      s0              = s0_current,
      survival        = survival_out,
      s0_leslie       = s0_leslie,
      final_N         = nrow(pop),
      years_simulated = year_counter,
      # DD flag
      density_dependence = FALSE,
      # Shared life-history parameters
      max_age         = max_age,
      pop_size        = pop_size,
      maturity_age    = maturity_age,
      litter_size     = litter_size,
      psi_nurse       = psi_nurse,
      psi_rest        = psi_rest,
      num_mates       = num_mates,
      female_fraction = female_fraction,
      infertility     = infertility,
      # Social structure
      pod_size        = pod_size,
      superpod_size   = superpod_size,
      stickiness_year = stickiness_year,
      male_behavior   = male_behavior,
      max_females     = max_females,
      weaning_age     = weaning_age
    ))

  } else {
    # -- Density dependence output --
    psi_nurse_K <- plogis(theta_current + rho)
    psi_rest_K  <- plogis(theta_current)

    # Solve for dd_max from target calving interval (unless dd_max supplied)
    if (is.null(dd_max)) {
      dd_max <- solve_dd_max(psi_nurse_K, psi_rest_K, target_interval,
                             target_depletion, z_pt, survival, wa_breed)
      message(sprintf(
        "Anchored dd_max = %.4f  (target: %.2f-yr interval at D=%.2f)",
        dd_max, target_interval, target_depletion
      ))
    }

    # Report emergent calving interval at K
    pi_1_K <- breeding_stationary(psi_nurse_K, psi_rest_K, survival, wa_breed)[1]
    interval_K <- 1 / pi_1_K

    message(sprintf(paste0(
      "\nConverged: theta = %.4f, rho = %.3f",
      "\n  psi_nurse_K=%.4f, psi_rest_K=%.4f",
      "\n  Calving interval at K = %.2f yr (emergent)",
      "\n  %d years simulated, final N = %s"),
      theta_current, rho,
      psi_nurse_K, psi_rest_K,
      interval_K,
      year_counter, format(nrow(pop), big.mark = ",")
    ))

    invisible(list(
      # Calibration results
      s0              = survival[1],
      survival        = survival,
      theta           = theta_current,
      theta_leslie    = theta_leslie,
      rho             = rho,
      psi_nurse_K     = psi_nurse_K,
      psi_rest_K      = psi_rest_K,
      interval_K      = interval_K,
      final_N         = nrow(pop),
      years_simulated = year_counter,
      # DD parameters
      density_dependence = TRUE,
      z_pt            = z_pt,
      dd_max          = dd_max,
      target_interval  = target_interval,
      target_depletion = target_depletion,
      K_1plus         = K_1plus,
      # Shared life-history parameters (psi set to K values for simulate.pop)
      max_age         = max_age,
      pop_size        = pop_size,
      maturity_age    = maturity_age,
      litter_size     = litter_size,
      psi_nurse       = psi_nurse_K,
      psi_rest        = psi_rest_K,
      num_mates       = num_mates,
      female_fraction = female_fraction,
      infertility     = infertility,
      # Social structure
      pod_size        = pod_size,
      superpod_size   = superpod_size,
      stickiness_year = stickiness_year,
      male_behavior   = male_behavior,
      max_females     = max_females,
      weaning_age     = weaning_age
    ))
  }
}


# ═══════════════════════════════════════════════════════════════════════════════
# ibm_growth_rate — lightweight IBM for measuring growth at a given depletion
# ═══════════════════════════════════════════════════════════════════════════════
#
# Runs a stripped-down individual-based simulation at target_depletion with
# density-dependent conception, fishing mortality, and orphan mortality.
# Returns the mean exponential growth rate over the measurement window.
#
# Used by solve_F_sustainable() to empirically bisect F.
# Not exported; internal use only.

ibm_growth_rate <- function(sim_config, F_val, selectivity,
                            n_years = 60L, pop_size = NULL) {

  # ── Extract parameters from sim_config ──
  max_age     <- sim_config$max_age
  survival    <- sim_config$survival
  litter_size <- sim_config$litter_size
  female_fraction <- sim_config$female_fraction
  infertility <- sim_config$infertility
  dd_max      <- sim_config$dd_max
  z_pt        <- sim_config$z_pt
  rho         <- sim_config$rho
  theta       <- sim_config$theta
  target_dep  <- sim_config$target_depletion
  wa_breed    <- if (is.null(sim_config$weaning_age)) 1L else
                   as.integer(sim_config$weaning_age)

  if (is.null(pop_size)) pop_size <- sim_config$pop_size

  if (length(infertility) == 1L) {
    infertility_f <- infertility_m <- infertility
  } else {
    infertility_f <- infertility[1]
    infertility_m <- infertility[2]
  }

  # ── Parse maturity ogive ──
  make_ogive <- function(x, max_a) {
    n <- max_a + 1L
    if (length(x) == 1L && x == round(x)) {
      ogive <- rep(0, n)
      if (x <= max_a) ogive[seq(x + 1L, n)] <- 1
      return(ogive)
    }
    if (length(x) == n) return(x)
    stop("maturity_age: ogive must have length max_age + 1.")
  }
  maturity_age <- sim_config$maturity_age
  if (is.list(maturity_age)) {
    ogive_f <- make_ogive(maturity_age$female, max_age)
    ogive_m <- make_ogive(maturity_age$male,   max_age)
  } else {
    ogive_f <- ogive_m <- make_ogive(maturity_age, max_age)
  }

  sample_mat_ages <- function(ogive, n) {
    pmf     <- diff(c(0, ogive))
    p_never <- 1 - sum(pmf)
    if (p_never > 0.001) {
      pmf  <- c(pmf, p_never)
      ages <- c(0:max_age, max_age + 1L)
    } else {
      ages <- 0:max_age
      pmf[length(pmf)] <- pmf[length(pmf)] + p_never
    }
    sample(ages, n, replace = TRUE, prob = pmf)
  }

  breeding_stationary <- function(psi_n, psi_r, sv, w) {
    n_st <- w + 2L
    P <- matrix(0, nrow = n_st, ncol = n_st)
    P[1, 2] <- 1
    for (k in seq_len(w)) {
      row   <- k + 1L
      ell_k <- sv[k]
      psi_k <- ell_k * psi_n + (1 - ell_k) * psi_r
      P[row, 1] <- psi_k
      if (k < w) {
        P[row, row + 1] <- ell_k * (1 - psi_n)
        P[row, n_st]    <- (1 - ell_k) * (1 - psi_r)
      } else {
        P[row, n_st] <- 1 - psi_k
      }
    }
    P[n_st, 1]    <- psi_r
    P[n_st, n_st] <- 1 - psi_r
    ev  <- eigen(t(P))
    idx <- which.min(abs(Mod(ev$values) - 1))
    pi  <- Mod(ev$vectors[, idx])
    pi / sum(pi)
  }

  # ── DD-shifted conception at target_depletion (for initial age structure) ──
  logit_psi_nurse_K <- theta + rho
  logit_psi_rest_K  <- theta
  shift_init        <- dd_max * (1 - target_dep^z_pt)
  psi_nurse_dep     <- plogis(logit_psi_nurse_K + shift_init)
  psi_rest_dep      <- plogis(logit_psi_rest_K  + shift_init)

  # ── Build fished Leslie matrix for stable age distribution ──
  surv_fished <- survival * exp(-F_val * selectivity)
  ogive_f_leslie <- c(0, ogive_f[seq_len(max_age)])

  pi_stat <- breeding_stationary(psi_nurse_dep, psi_rest_dep,
                                  surv_fished, wa_breed)
  pi_1_dep <- pi_stat[1]
  n_breed_states <- length(pi_stat)

  ff <- litter_size * female_fraction * pi_1_dep * (1 - infertility_f)
  n_classes <- max_age + 1L

  A <- matrix(0, n_classes, n_classes)
  A[1, ] <- ogive_f_leslie * ff
  for (i in seq_len(max_age)) A[i + 1L, i] <- surv_fished[i]

  w_eig    <- Mod(eigen(A)$vectors[, 1])
  stable_A <- w_eig / sum(w_eig)

  # ── K_1plus scaled to this pop_size ──
  # At K, the unfished stable age distribution determines K_1plus
  psi_nurse_K <- plogis(logit_psi_nurse_K)
  psi_rest_K  <- plogis(logit_psi_rest_K)
  pi_stat_K   <- breeding_stationary(psi_nurse_K, psi_rest_K, survival, wa_breed)
  pi_1_K      <- pi_stat_K[1]
  ff_K        <- litter_size * female_fraction * pi_1_K * (1 - infertility_f)

  A_K <- matrix(0, n_classes, n_classes)
  A_K[1, ] <- ogive_f_leslie * ff_K
  for (i in seq_len(max_age)) A_K[i + 1L, i] <- survival[i]
  w_K      <- Mod(eigen(A_K)$vectors[, 1])
  stable_K <- w_K / sum(w_K)
  K_1plus  <- pop_size * (1 - stable_K[1])

  # ── Initialise population at target_depletion ──
  init_pop_size <- round(pop_size * target_dep)
  init_N    <- pmax(round(stable_A * init_pop_size), 0L)
  init_ages <- rep(0:max_age, times = init_N)
  n_init    <- sum(init_N)

  if (n_init == 0L) return(NA_real_)

  init_sex <- sample(c("F", "M"), n_init,
                     prob = c(female_fraction, 1 - female_fraction),
                     replace = TRUE)

  init_mat_age <- integer(n_init)
  is_f <- init_sex == "F"
  init_mat_age[is_f]  <- sample_mat_ages(ogive_f, sum(is_f))
  init_mat_age[!is_f] <- sample_mat_ages(ogive_m, sum(!is_f))

  init_fertile <- rep(TRUE, n_init)
  if (infertility_f > 0) init_fertile[is_f]  <- runif(sum(is_f))  >= infertility_f
  if (infertility_m > 0) init_fertile[!is_f] <- runif(sum(!is_f)) >= infertility_m

  # Breeding state assignment from stationary distribution
  init_breed_state <- rep(NA_integer_, n_init)
  init_s2_year     <- rep(NA_integer_, n_init)
  init_calf_id     <- rep(0L, n_init)

  mature_f <- which(is_f & init_ages >= init_mat_age & init_fertile)
  if (length(mature_f) > 0L) {
    state_idx <- sample(seq_len(n_breed_states), length(mature_f),
                        replace = TRUE, prob = pi_stat)
    bs <- ifelse(state_idx == 1L, 1L,
                 ifelse(state_idx == n_breed_states, 3L, 2L))
    s2y <- ifelse(bs == 2L, state_idx - 1L, NA_integer_)
    init_breed_state[mature_f] <- bs
    init_s2_year[mature_f]     <- s2y
  }

  pop <- data.table::data.table(
    id          = seq_len(n_init),
    age         = init_ages,
    sex         = init_sex,
    mat_age     = init_mat_age,
    mother_id   = 0L,
    breed_state = init_breed_state,
    fertile     = init_fertile,
    calf_id     = init_calf_id,
    s2_year     = init_s2_year
  )

  next_id    <- n_init + 1L
  surv_vec   <- survival
  sel_vec    <- selectivity
  use_weaning <- wa_breed > 0L
  weaning_age <- wa_breed

  # ── Main simulation loop ──
  N_vec    <- integer(n_years)
  burn_in  <- 20L  # let age structure settle before measuring

  for (yr in seq_len(n_years)) {

    # 1. Survival (natural + fishing)
    if (F_val > 0) {
      rates <- surv_vec[pop$age + 1L] *
        exp(-F_val * sel_vec[pop$age + 1L])
    } else {
      rates <- surv_vec[pop$age + 1L]
    }
    alive <- runif(nrow(pop)) <= rates
    pop   <- pop[alive]

    # 2. Aging + max_age removal
    data.table::set(pop, j = "age", value = pop$age + 1L)
    pop <- pop[pop$age <= max_age]

    if (nrow(pop) == 0L) return(NA_real_)

    # 2b. Orphan mortality
    if (use_weaning) {
      dep_idx <- which(pop$age <= weaning_age & pop$mother_id != 0L)
      if (length(dep_idx) > 0L) {
        orphans <- dep_idx[!pop$mother_id[dep_idx] %in% pop$id]
        if (length(orphans) > 0L) pop <- pop[-orphans]
      }
    }

    # 3. DD conception adjustment
    N_1plus <- sum(pop$age >= 1L)
    D_t     <- N_1plus / K_1plus
    delta_t <- dd_max * (1 - D_t^z_pt)
    psi_nurse_yr <- plogis(logit_psi_nurse_K + delta_t)
    psi_rest_yr  <- plogis(logit_psi_rest_K  + delta_t)

    # 4. Markov breeding state transitions
    # Newly mature females enter at S3
    new_mature <- which(pop$sex == "F" & pop$age == pop$mat_age &
                          is.na(pop$breed_state) & pop$fertile)
    if (length(new_mature) > 0L) {
      data.table::set(pop, i = new_mature, j = "breed_state", value = 3L)
    }

    # Snapshot states before transitions
    mother_idx <- which(pop$breed_state == 1L)
    s2_idx     <- which(pop$breed_state == 2L)
    s3_idx     <- which(pop$breed_state == 3L)

    # S2 transitions: calf-survival-dependent
    if (length(s2_idx) > 0L) {
      calf_ids <- pop$calf_id[s2_idx]
      calf_rows  <- match(calf_ids, pop$id)
      calf_alive <- !is.na(calf_rows)

      # Founder S2 mothers (calf_id = 0): stochastic calf survival
      founder_s2 <- calf_ids == 0L
      if (any(founder_s2)) {
        k_founder <- pop$s2_year[s2_idx[founder_s2]]
        calf_alive[founder_s2] <- runif(sum(founder_s2)) < surv_vec[k_founder]
      }

      k <- pop$s2_year[s2_idx]
      psi_eff <- ifelse(calf_alive, psi_nurse_yr, psi_rest_yr)
      conceive <- runif(length(s2_idx)) < psi_eff

      new_state <- rep(3L, length(s2_idx))
      new_state[conceive] <- 1L
      stay_s2 <- !conceive & calf_alive & k < wa_breed
      new_state[stay_s2] <- 2L

      data.table::set(pop, i = s2_idx, j = "breed_state", value = new_state)

      new_s2y <- rep(NA_integer_, length(s2_idx))
      new_s2y[stay_s2] <- k[stay_s2] + 1L
      data.table::set(pop, i = s2_idx, j = "s2_year", value = new_s2y)

      leaving_s2 <- which(new_state != 2L)
      if (length(leaving_s2) > 0L) {
        data.table::set(pop, i = s2_idx[leaving_s2], j = "calf_id", value = 0L)
      }
    }

    # S3 transitions
    if (length(s3_idx) > 0L) {
      new_state_s3 <- ifelse(runif(length(s3_idx)) < psi_rest_yr, 1L, 3L)
      data.table::set(pop, i = s3_idx, j = "breed_state", value = new_state_s3)
    }

    # S1 -> S2: give birth
    if (length(mother_idx) > 0L) {
      data.table::set(pop, i = mother_idx, j = "breed_state", value = 2L)
      data.table::set(pop, i = mother_idx, j = "s2_year",     value = 1L)
    }

    # 5. Create offspring (simplified: no father tracking, no pods)
    has_males <- any(pop$sex == "M" & pop$age >= pop$mat_age & pop$fertile)

    if (length(mother_idx) > 0L && has_males) {
      n_mothers    <- length(mother_idx)
      litter_sizes <- 1L + rpois(n_mothers, lambda = litter_size - 1)
      n_yoy        <- sum(litter_sizes)

      yoy_sex <- sample(c("F", "M"), n_yoy,
                        prob = c(female_fraction, 1 - female_fraction),
                        replace = TRUE)

      yoy_mat_age <- integer(n_yoy)
      yoy_is_f    <- yoy_sex == "F"
      if (any(yoy_is_f))  yoy_mat_age[yoy_is_f]  <- sample_mat_ages(ogive_f, sum(yoy_is_f))
      if (any(!yoy_is_f)) yoy_mat_age[!yoy_is_f] <- sample_mat_ages(ogive_m, sum(!yoy_is_f))

      yoy_fertile <- rep(TRUE, n_yoy)
      if (infertility_f > 0 && any(yoy_is_f))
        yoy_fertile[yoy_is_f]  <- runif(sum(yoy_is_f))  >= infertility_f
      if (infertility_m > 0 && any(!yoy_is_f))
        yoy_fertile[!yoy_is_f] <- runif(sum(!yoy_is_f)) >= infertility_m

      yoy_ids <- seq.int(next_id, length.out = n_yoy)

      # mother_id for each offspring
      off_mother_idx <- rep(seq_len(n_mothers), times = litter_sizes)
      off_mother_id  <- pop$id[mother_idx[off_mother_idx]]

      yoy <- data.table::data.table(
        id          = yoy_ids,
        age         = 0L,
        sex         = yoy_sex,
        mat_age     = yoy_mat_age,
        mother_id   = off_mother_id,
        breed_state = NA_integer_,
        fertile     = yoy_fertile,
        calf_id     = 0L,
        s2_year     = NA_integer_
      )

      # Set calf_id on mothers (first offspring is the dependent calf)
      first_yoy_per_mother <- yoy_ids[!duplicated(off_mother_idx)]
      data.table::set(pop, i = mother_idx, j = "calf_id",
                      value = first_yoy_per_mother)

      next_id <- next_id + n_yoy
      pop <- data.table::rbindlist(list(pop, yoy), use.names = TRUE)
    }

    N_vec[yr] <- nrow(pop)
  }

  # ── Compute growth rate from measurement window (skip burn-in) ──
  meas_start <- min(burn_in + 1L, n_years)
  meas_N     <- N_vec[meas_start:n_years]

  if (length(meas_N) < 2L || any(meas_N == 0L)) return(NA_real_)

  mean(diff(log(meas_N)))
}


#' Solve for the sustained fishing mortality that holds a DD population at its
#' target depletion
#'
#' Given a calibrated density-dependent configuration (from
#' \code{create.stable.pop(density_dependence = TRUE)}), this function finds the
#' constant instantaneous fishing mortality \eqn{F} such that the population's
#' growth rate is exactly zero at \code{sim_config$target_depletion}.  At that
#' depletion, density-dependent compensation boosts conception rates above the
#' carrying-capacity baseline; \eqn{F} removes exactly the surplus, so the
#' population sits at a stable, self-correcting equilibrium.
#'
#' The solution uses a two-phase approach: (1) a fast Leslie-matrix solve for
#' an analytical starting estimate, then (2) an IBM-based bisection that
#' empirically finds the F where a lightweight stochastic simulation produces
#' near-zero growth at \code{target_depletion}.  The IBM phase captures growth
#' sinks (orphan mortality, demographic stochasticity) that the Leslie matrix
#' misses.  If the IBM simulation is unstable (e.g., extreme test parameters),
#' the function falls back to the Leslie estimate with a warning.
#'
#' The function also returns the stable age distribution of the
#' \strong{fished} population, which \code{simulate.pop()} uses when
#' \code{init_depletion} is set to skip the fishing-down transient and start
#' directly at the depleted equilibrium.
#'
#' @param sim_config List returned by \code{create.stable.pop()} with
#'   \code{density_dependence = TRUE}.
#' @param selectivity Numeric vector of length \code{max_age + 1}, values in
#'   \code{[0, 1]}.  Age-specific vulnerability to the fishery.  Default:
#'   \code{c(0, rep(1, max_age))} -- flat across ages 1+, age-0 excluded
#'   (matches the purse-seine assumption that calves are not independently
#'   vulnerable).
#'
#' @return A named list:
#' \describe{
#'   \item{F_sustainable}{The constant instantaneous fishing mortality (apical
#'     rate) that holds the population stationary at \code{target_depletion}.
#'     Calibrated via IBM bisection when possible; Leslie fallback otherwise.}
#'   \item{F_leslie}{The Leslie-matrix estimate of sustainable F (analytical,
#'     typically higher than \code{F_sustainable} because the Leslie matrix
#'     does not capture IBM-specific growth sinks).}
#'   \item{selectivity}{The selectivity vector used (returned for downstream
#'     use in \code{simulate.pop()}).}
#'   \item{target_depletion}{Copied from \code{sim_config}.}
#'   \item{interval_at_equilibrium}{Theoretical mean calving interval at the
#'     fished equilibrium.}
#'   \item{stable_age_fished}{Stable age distribution of the fished population
#'     (eigenvector of the fished Leslie matrix).  Used by
#'     \code{simulate.pop()} when \code{init_depletion} is set.}
#' }
#'
#' @details
#' \strong{How it works.}
#' \enumerate{
#'   \item Compute DD-shifted conception rates at \code{target_depletion}:
#'     \code{psi_rest_D = inv_logit(theta + shift)},
#'     \code{psi_nurse_D = inv_logit(theta + rho + shift)}, where
#'     \code{shift = dd_max * (1 - D^z)}.
#'   \item \strong{Phase 1 (Leslie):} For a candidate \eqn{F}, build a Leslie
#'     matrix with survival \code{S(a) = survival[a] * exp(-F * selectivity[a])}
#'     and fecundity derived from the DD-shifted breeding stationary
#'     distribution.  \code{uniroot()} finds the \eqn{F} where the dominant
#'     eigenvalue equals 1.  This gives \code{F_leslie}.
#'   \item \strong{Phase 2 (IBM):} A lightweight IBM simulation
#'     (\code{ibm_growth_rate()}) runs at \code{target_depletion} with DD,
#'     fishing, and orphan mortality to measure the actual growth rate.
#'     Binary search over \eqn{F} finds the value where IBM growth is
#'     approximately zero.  This gives \code{F_sustainable}, which is
#'     typically lower than \code{F_leslie}.
#' }
#'
#' \strong{Self-correcting equilibrium.}  The solved \eqn{F} produces a stable
#' equilibrium: if stochastic fluctuation pushes depletion below the target,
#' compensation strengthens (growth exceeds removals, population recovers); if
#' depletion rises above the target, compensation weakens (removals exceed
#' growth, population declines back).
#'
#' @importFrom stats uniroot plogis
#' @export
solve_F_sustainable <- function(sim_config, selectivity = NULL) {

  # ── Input validation ──
  if (!is.list(sim_config))
    stop("`sim_config` must be a list (output of create.stable.pop()).")
  if (!isTRUE(sim_config$density_dependence))
    stop("`solve_F_sustainable()` requires `density_dependence = TRUE`.")
  if (is.null(sim_config$target_depletion))
    stop("`sim_config$target_depletion` is missing. Was `dd_max` supplied ",
         "directly instead of `target_interval`?")

  max_age  <- sim_config$max_age
  survival <- sim_config$survival

  if (is.null(selectivity)) {
    selectivity <- c(0, rep(1, max_age))
  }
  if (!is.numeric(selectivity) || length(selectivity) != max_age + 1L)
    stop("`selectivity` must be a numeric vector of length max_age + 1 (",
         max_age + 1L, ").")
  if (any(selectivity < 0 | selectivity > 1))
    stop("`selectivity` values must be between 0 and 1.")

  # ── Extract parameters ──
  theta_leslie <- sim_config$theta_leslie
  theta_ibm    <- sim_config$theta
  rho          <- sim_config$rho
  dd_max     <- sim_config$dd_max
  z_pt       <- sim_config$z_pt
  target_dep <- sim_config$target_depletion
  wa_breed   <- if (is.null(sim_config$weaning_age)) 1L else
                  as.integer(sim_config$weaning_age)
  litter_size     <- sim_config$litter_size
  female_fraction <- sim_config$female_fraction

  infertility <- sim_config$infertility
  infertility_f <- if (length(infertility) == 1L) infertility else infertility[1]

  # ── Parse maturity ogive ──
  make_ogive <- function(x, max_a) {
    n <- max_a + 1L
    if (length(x) == 1L && x == round(x)) {
      ogive <- rep(0, n)
      if (x <= max_a) ogive[seq(x + 1L, n)] <- 1
      return(ogive)
    }
    if (length(x) == n) return(x)
    stop("maturity_age: ogive must have length max_age + 1 (", n, ").")
  }
  maturity_age <- sim_config$maturity_age
  if (is.list(maturity_age)) {
    ogive_f <- make_ogive(maturity_age$female, max_age)
  } else {
    ogive_f <- make_ogive(maturity_age, max_age)
  }
  ogive_f_leslie <- c(0, ogive_f[seq_len(max_age)])

  # ── Local helpers ──
  mat_lambda <- function(A) Mod(eigen(A, only.values = TRUE)$values[1])
  mat_stable <- function(A) {
    w <- Mod(eigen(A)$vectors[, 1])
    w / sum(w)
  }

  breeding_stationary <- function(psi_n, psi_r, sv, w) {
    n_st <- w + 2L
    P <- matrix(0, nrow = n_st, ncol = n_st)
    P[1, 2] <- 1
    for (k in seq_len(w)) {
      row   <- k + 1L
      ell_k <- sv[k]
      psi_k <- ell_k * psi_n + (1 - ell_k) * psi_r
      P[row, 1] <- psi_k
      if (k < w) {
        P[row, row + 1] <- ell_k * (1 - psi_n)
        P[row, n_st]    <- (1 - ell_k) * (1 - psi_r)
      } else {
        P[row, n_st] <- 1 - psi_k
      }
    }
    P[n_st, 1]    <- psi_r
    P[n_st, n_st] <- 1 - psi_r
    ev  <- eigen(t(P))
    idx <- which.min(abs(Mod(ev$values) - 1))
    pi  <- Mod(ev$vectors[, idx])
    pi / sum(pi)
  }

  # ── DD-shifted conception at target depletion ──
  # Use theta_leslie for the Leslie matrix: it makes lambda = 1 at K, so the
  # DD boost at D < 1 correctly produces lambda > 1 (surplus growth).
  shift       <- dd_max * (1 - target_dep^z_pt)
  psi_nurse_D <- plogis(theta_leslie + rho + shift)
  psi_rest_D  <- plogis(theta_leslie + shift)

  n_classes <- max_age + 1L

  # ── Build fished Leslie matrix (with orphan mortality correction) ──
  # Dependent calves (age <= wa_breed post-aging) die if their mother dies.
  # Natural orphan mortality is already absorbed by the calibrated theta

  # (Phase 2 bisection accounts for it).  Here we add only the INCREMENTAL
  # orphan mortality from fishing: P(mother killed by fishing), averaged over
  # mother ages weighted by survivorship × maturity ogive.

  # Survivorship l(a): probability of reaching age a under natural mortality
  lx <- cumprod(c(1, survival[seq_len(max_age)]))
  mom_wt <- lx * ogive_f_leslie
  mom_wt_sum <- sum(mom_wt)

  build_leslie_F <- function(FF) {
    sv_F <- survival * exp(-FF * selectivity)

    # Fishing-only maternal survival: exp(-F * selectivity) per mother age
    fish_surv <- exp(-FF * selectivity)
    avg_mom_fish_surv <- if (mom_wt_sum > 0) {
      sum(mom_wt * fish_surv) / mom_wt_sum
    } else 1

    A <- matrix(0, n_classes, n_classes)
    for (i in seq_len(max_age)) {
      s_i <- sv_F[i]
      if (i <= wa_breed) s_i <- s_i * avg_mom_fish_surv
      A[i + 1L, i] <- s_i
    }
    pi_1 <- breeding_stationary(psi_nurse_D, psi_rest_D, sv_F, wa_breed)[1]
    ff   <- litter_size * female_fraction * pi_1 * (1 - infertility_f)
    A[1, ] <- ogive_f_leslie * ff
    A
  }

  # ── Leslie-based F solve (analytical starting estimate) ──
  lambda_at_F <- function(FF) mat_lambda(build_leslie_F(FF))

  lam_0 <- lambda_at_F(0)
  if (lam_0 <= 1)
    stop(sprintf(
      paste0("Population cannot sustain fishing at D=%.2f ",
             "(lambda=%.4f without fishing). ",
             "Check dd_max and target_depletion."),
      target_dep, lam_0
    ))

  F_leslie <- uniroot(function(FF) lambda_at_F(FF) - 1,
                       c(0, 2), tol = 1e-8)$root

  # ── IBM-based F bisection ──
  # The Leslie solve overestimates sustainable F because the IBM has growth

  # sinks (orphan mortality, demographic stochasticity) that the Leslie
  # matrix does not capture.  Use a lightweight IBM simulation to empirically
  # find F where growth ≈ 0 at target_depletion.
  ibm_pop_size <- 50000L

  message(sprintf(
    "Leslie F = %.4f. Starting IBM bisection (pop_size=%s)...",
    F_leslie, format(ibm_pop_size, big.mark = ",")
  ))

  # First check: can the IBM sustain any fishing at this depletion?
  growth_at_zero <- ibm_growth_rate(sim_config, F_val = 0,
                                     selectivity = selectivity,
                                     pop_size = ibm_pop_size)

  # If IBM simulation crashes (NA), fall back to Leslie estimate with warning
  ibm_ok <- !is.na(growth_at_zero)

  if (ibm_ok && growth_at_zero <= 0) {
    stop(sprintf(
      paste0("IBM growth at D=%.2f with F=0 is %.5f (non-positive). ",
             "dd_max is too weak for the IBM to produce surplus growth. ",
             "Consider IBM-based dd_max calibration or removing orphan ",
             "mortality."),
      target_dep, growth_at_zero
    ))
  }

  if (!ibm_ok) {

    warning("IBM growth check returned NA (population crashed). ",
            "Falling back to Leslie estimate.")
    F_solved <- F_leslie

  } else {

    message(sprintf("  IBM growth at F=0: %+.5f (positive — bisecting F...)",
                    growth_at_zero))

    F_low  <- 0
    F_high <- F_leslie * 2.5
    F_tol  <- 1e-4
    growth_tol <- 0.001
    max_iter   <- 30L
    F_mid      <- F_leslie  # start at Leslie estimate

    for (iter in seq_len(max_iter)) {
      growth <- ibm_growth_rate(sim_config, F_val = F_mid,
                                 selectivity = selectivity,
                                 pop_size = ibm_pop_size)

      # If IBM crashes mid-bisection, fall back to best estimate so far
      if (is.na(growth)) {
        warning(sprintf(
          "IBM sim returned NA at F=%.4f. Using current midpoint.", F_mid
        ))
        break
      }

      message(sprintf("  iter %2d  F=%.4f  growth=%+.5f", iter, F_mid, growth))

      if (abs(growth) < growth_tol) break

      if (growth > 0) {
        F_low <- F_mid
      } else {
        F_high <- F_mid
      }
      F_mid <- (F_low + F_high) / 2

      if ((F_high - F_low) < F_tol) break
    }

    F_solved <- F_mid
  }

  # ── Fished stable age distribution (from Leslie, for init_depletion) ──
  A_fished      <- build_leslie_F(F_solved)
  stable_fished <- mat_stable(A_fished)
  sv_F_eq       <- survival * exp(-F_solved * selectivity)
  pi_1          <- breeding_stationary(psi_nurse_D, psi_rest_D, sv_F_eq, wa_breed)[1]
  interval_eq   <- 1 / pi_1

  message(sprintf(
    paste0("Sustainable F = %.4f at D = %.2f  (calving interval = %.2f yr)\n",
           "  Leslie estimate was F = %.4f"),
    F_solved, target_dep, interval_eq, F_leslie
  ))

  list(
    F_sustainable          = F_solved,
    F_leslie               = F_leslie,
    selectivity            = selectivity,
    target_depletion       = target_dep,
    interval_at_equilibrium = interval_eq,
    stable_age_fished      = stable_fished
  )
}
