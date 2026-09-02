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
#'     and the function bisects a logit-scale offset on the conception
#'     probabilities (\code{psi_nurse}, \code{psi_rest}) until the population is
#'     stable at carrying capacity K.  In \code{simulate.pop()}, a
#'     Pella--Tomlinson compensation mechanism then adjusts conception rates
#'     dynamically based on depletion.
#' }
#'
#' Both modes use a two-phase approach:
#' \enumerate{
#'   \item \strong{Leslie matrix} -- analytically solves for the starting
#'     estimate (s0 or theta_shift).
#'   \item \strong{Iterative simulation} -- refines via bisection until
#'     convergence.
#' }
#'
#' @param max_age Integer. Maximum age in the population. Individuals will breed at this age, but will be removed from the population before breeding again.
#' @param survival Numeric vector of length \code{max_age + 1}. Annual survival
#'   probabilities for ages 0 through \code{max_age}.  When
#'   \code{density_dependence = FALSE}, the age-0 value is a placeholder that
#'   will be replaced by the calibrated estimate.  When
#'   \code{density_dependence = TRUE}, all values including age-0 are used as
#'   supplied.
#' @param pop_size Integer. Total starting population size. Equal to carrying capacity K when density dependence is active.
#' @param maturity_age Maturity specification. Can be:
#'   \itemize{
#'     \item An integer: knife-edged maturity at that age for both sexes.
#'     \item A numeric vector of length \code{max_age + 1}: cumulative
#'       probability of being mature at each age (ogive CDF), applied to both
#'       sexes.
#'     \item A list with \code{female} and \code{male} elements, each of which
#'       can be an integer or numeric vector as above.
#'   }
#' @param litter_size Numeric. Mean litter size per breeding female.
#' @param psi_nurse Numeric in \[0, 1\]. Conception probability while nursing
#'   a dependent calf (lactational suppression).  This is the suppressed rate;
#'   mothers whose calves die are released to \code{psi_rest}.
#' @param psi_rest Numeric in \[0, 1\]. Conception probability while resting
#'   (no dependent calf), or after calf death.  Should be >= \code{psi_nurse}.
#' @param num_mates Integer. Number of mates per female per breeding cycle.  Passed through to
#'   \code{simulate.pop()}.
#' @param female_fraction Numeric in (0, 1). Fraction of offspring that are female.
#' @param infertility Numeric (scalar or length-2 vector \code{c(female, male)}). Proportion of permanently infertile individuals.  Default 0.
#' @param pod_size Integer or NULL. Number of individuals per pod.
#' @param superpod_size Integer or NULL. Number of pods per superpod.
#' @param stickiness_year Numeric (scalar or length-2 vector \code{c(female, male)}) or NULL. Between-year superpod fidelity i.e., probability that individuals will switch superpods between simulation years.
#' @param male_behavior Character or NULL. \code{"random"} or
#'   \code{"strong_bull"}. "random" means that all mature males are eligible to mate with mature females from their superpod; "strong_bull" means that just one randomly selected male will be eligible to mate with mature females from their superpod. The "strong_bull" status will remain with the male who receives that designation until he dies, at which point another randomly selected male will take over as the strong bull for that superpod.
#' @param max_females Integer or NULL. Per-year cap on matings per male.
#' @param weaning_age Integer or NULL. Age at independence from mother.  Also
#'   determines the maximum number of years a mother stays in the "with
#'   dependent calf" breeding state (S2).  If NULL, defaults to 1-year
#'   dependency for the breeding cycle (no cow-calf social following).
#' @param density_dependence Logical. If \code{TRUE}, use Pella--Tomlinson
#'   density dependence instead of s0 calibration.  Default \code{FALSE}.
#' @param z_pt Numeric. Pella--Tomlinson shape parameter.  Default 2.39
#'   (IWC convention; MNPL at ~0.6 K).  Only used when
#'   \code{density_dependence = TRUE}. The default value for this parameter was modeled after cetaceans.
#' @param dd_max Numeric or NULL. Maximum logit-scale shift on conception
#'   probabilities at zero depletion (D -> 0).  Controls compensation strength.
#'   Larger values allow faster recovery from depletion.  When
#'   \code{density_dependence = TRUE}, supply either \code{dd_max} (direct) or
#'   \code{target_interval} (anchored), not both.
#'
#'   \strong{In most cases, leave this \code{NULL} and use \code{target_interval}
#'   / \code{target_depletion} instead.}  \code{dd_max} has no biologically
#'   meaningful default; a logit-scale shift is not something most users can
#'   reason about directly, and a plausible-looking value (e.g. 2) can be far
#'   too weak for a given life history (see \code{target_interval} below and
#'   the anchoring discussion in \code{@details}).  Supplying \code{dd_max}
#'   directly is intended for sensitivity analyses (sweeping compensation
#'   strength) or purely theoretical runs with no real calving-interval data
#'   to anchor to.  Supplying it will print a message recommending
#'   \code{target_interval} unless you are clearly running a sweep.
#' @param target_interval Numeric or NULL. \strong{Recommended over
#'   \code{dd_max} for real species work.} Target calving interval in years at
#'   the reference depletion \code{target_depletion}.  When provided (and
#'   \code{density_dependence = TRUE}), the function solves for the
#'   \code{dd_max} that reproduces this interval, rather than requiring
#'   \code{dd_max} directly.  Must be >= 2 (structural floor: 1 year pregnancy
#'   + 1 year dependency) and shorter than the calving interval at carrying
#'   capacity K.  Anchoring to an empirically observed calving interval is far
#'   more defensible than guessing a logit-scale shift.
#' @param target_depletion Numeric in (0, 1). Reference depletion level at
#'   which \code{target_interval} was observed.  Default 0.3 (30\% of K).
#'   Only used when \code{target_interval} is provided.
#' @param check_interval Integer. Years per assessment window. Default 10.
#' @param growth_tol Numeric. Stability tolerance. Default 0.001.
#' @param stable_required Integer. Consecutive stable windows needed. Default 5.
#' @param max_windows Integer. Safety cap on windows. Default 100.
#'
#' @details
#'
#' ## Calibration Modes
#'
#' **s0 calibration** (\code{density_dependence = FALSE}, default): This function
#' searches for the age-0 survival rate (s0) that results in a stable population.
#' The user provides a complete survival vector with any placeholder value for age 0;
#' this function replaces it with the calibrated estimate.  The algorithm uses a
#' two-phase approach:
#' \enumerate{
#'   \item Leslie matrix analytically estimates the starting s0 (or gives a bounds
#'     check for the bisection).
#'   \item Bisection iteratively refines s0 over assessment windows of
#'     \code{check_interval} years until population growth is near zero
#'     (\code{|r| < growth_tol}) for \code{stable_required} consecutive windows.
#' }
#' The Leslie matrix fecundity row uses a corrected maturity ogive (shifted right
#' by one year) so that newly mature females contribute zero fecundity at their
#' first mature age, matching the IBM behavior.
#'
#' **Density dependence mode** (\code{density_dependence = TRUE}): Instead of
#' calibrating s0, the user supplies a complete, empirically-informed survival
#' curve (including age-0) and this function calibrates a logit-scale offset
#' (\code{theta_shift}) on conception rates (\code{psi_nurse}, \code{psi_rest}).
#' The result is a population that is stable at carrying capacity K, with
#' Pella-Tomlinson compensation active in \code{simulate.pop()}.  Outputs include
#' \code{theta_shift}, \code{theta_shift_leslie} (analytical), \code{psi_nurse_K},
#' \code{psi_rest_K} (reference conception rates at K).  The calibration solves
#' for the offset that makes lambda(K) = 1 (growth rate at carrying capacity = 1).
#'
#' ## Two-Phase Approach
#'
#' **Phase 1: Leslie matrix** (analytical estimate)
#' Constructs an age-structured Leslie matrix based on the stable age distribution
#' and fecundity. Computes the dominant eigenvalue (growth rate lambda) and
#' extracts the stable age distribution. Uses this to estimate the starting s0
#' (or theta_shift) for bisection.  The Leslie matrix is built using:
#' \itemize{
#'   \item Survival probabilities from the \code{survival} input (or a trial value
#'     for s0/theta_shift in DD mode).
#'   \item Fecundity = (maturity ogive) * (breeding stationary proportion) * (litter_size).
#'   \item For s0 calibration: solves analytically for s0 that makes lambda = 1.
#'   \item For DD: uses uniroot to find theta_shift that makes lambda = 1.
#' }
#'
#' **Phase 2: Iterative bisection** (refinement via simulation)
#' The Leslie estimate (or bounds) seed a bisection search. The algorithm:
#' \enumerate{
#'   \item Runs \code{check_interval} years of simulation at the current candidate
#'     s0 (or theta_shift).
#'   \item Assesses population growth rate r over the window.
#'   \item If r > \code{growth_tol}: s0 is too high (or theta_shift too high);
#'     bisect down.
#'   \item If r < \code{-growth_tol}: s0 is too low (or theta_shift too low);
#'     bisect up.
#'   \item If \code{|r| < growth_tol}: window is stable; increment stable counter.
#'   \item Continues until \code{stable_required} consecutive stable windows occur,
#'     or \code{max_windows} is reached (safety limit).
#' }
#' Total calibration time is typically 100--220 years. Default tolerances
#' (\code{check_interval = 10}, \code{growth_tol = 0.001}, \code{stable_required = 5})
#' balance precision and computational cost.
#'
#' ## Initial Population Structure
#'
#' The initial population is sampled from the Leslie-stable age distribution.
#' Every individual is assigned:
#' \itemize{
#'   \item A unique ID (1, 2, ..., pop_size).
#'   \item Sex: Bernoulli(female_fraction).
#'   \item Age: Drawn from stable age distribution.
#'   \item Maturity age (\code{mat_age}): Drawn from maturity ogive PMF at birth.
#'     Later, the individual matures when calendar age >= mat_age.
#'   \item Fertility: Bernoulli(1 - infertility); if infertile, the individual
#'     never breeds.
#'   \item Pod and superpod: Assigned deterministically from fixed mapping,
#'     starting with pod 1 to superpod 1, pod 2 to superpod 1, etc.
#'   \item Breeding state (for females): S1, S2, or S3, drawn from the breeding
#'     stationary distribution. S2 females are assigned a placeholder calf (age 0,
#'     sex random, generated on the fly).
#'   \item mother_id = 0, father_id = 0 (founders; marked for later flush).
#' }
#'
#' ## Markov Breeding Cycle
#'
#' The function uses a numerically-derived Markov breeding stationary distribution
#' to initialize and reason about breeding state proportions. The breeding cycle
#' is parameterized by:
#' \itemize{
#'   \item \code{psi_nurse}, \code{psi_rest}: Conception rates while nursing vs.
#'     resting.
#'   \item \code{weaning_age}: Maximum years an S2 mother stays with her calf
#'     (default 1 if NULL).
#' }
#' The stationary distribution is computed from a Markov transition matrix built
#' by \code{breeding_stationary()} helper, accounting for calf survival rates. In
#' \code{simulate.pop()}, this ratio between conception states drives the actual
#' age-specific birth rates.
#'
#' ## Social Structure Initialization
#'
#' Pods and superpods are created deterministically during initialization. If
#' \code{pod_size = 20} and \code{superpod_size = 10}, then:
#' \itemize{
#'   \item Pods 1--10 belong to superpod 1.
#'   \item Pods 11--20 belong to superpod 2.
#'   \item Etc.
#' }
#' This mapping is stored in the \code{pod_to_sp} vector and passed to
#' \code{simulate.pop()} and \code{sample.pop()}.
#'
#' ## Mating and Fertility
#'
#' The function sets up the framework for \code{simulate.pop()} to handle mating:
#' \itemize{
#'   \item \code{num_mates}: Number of males each female mates with.
#'   \item \code{male_behavior}: "random" (polyandry from superpod pool) or
#'     "strong_bull" (one persistent male per superpod).
#'   \item \code{max_females}: Per-year cap on a male's mating partners (prevents
#'     monopolization).
#'   \item \code{infertility}: Proportion of infertile individuals (same or
#'     sex-specific). Infertile females never breed; infertile males are excluded
#'     from mating pools.
#' }
#'
#' ## Density Dependence Parameters
#'
#' When \code{density_dependence = TRUE}, the function calibrates:
#' \itemize{
#'   \item \code{theta_shift}: Logit-scale offset on conception rates that achieves
#'     stability at K.
#'   \item \code{psi_nurse_K}, \code{psi_rest_K}: Reference conception rates at K
#'     (shifted by theta_shift from user inputs).
#'   \item \code{z_pt}: Pella-Tomlinson shape (default 2.39, IWC convention places
#'     MNPL at 0.6K).
#'   \item \code{dd_max}: Maximum logit shift at zero depletion (controls
#'     compensation strength). Either supplied directly or solved from
#'     \code{target_interval}.
#' }
#'
#' **Anchoring dd_max to a target calving interval (recommended):** Rather
#' than choosing \code{dd_max} directly, supply \code{target_interval}
#' (observed calving interval in years) and \code{target_depletion} (the
#' depletion level at which that interval was observed).  After the
#' theta_shift calibration finds the at-K conception rates
#' (\code{psi_nurse_K}, \code{psi_rest_K}), the function solves for the
#' \code{dd_max} that reproduces the target interval at the reference
#' depletion using \code{uniroot()}.  This is far more principled than
#' choosing \code{dd_max} directly, because:
#' \itemize{
#'   \item Users can reason about calving intervals more intuitively than
#'     logit-scale shift magnitudes -- \code{dd_max} has no biologically
#'     meaningful default value.
#'   \item The anchored value is species-specific and grounded in empirical
#'     demographic data, rather than guessed.
#'   \item Plausible-looking placeholder values (e.g., dd_max = 2) can be far
#'     too weak for a given species' life history, producing a population
#'     that cannot recover from depletion.
#' }
#' \strong{In most cases, leave \code{dd_max = NULL} and use
#' \code{target_interval} / \code{target_depletion} instead.}  Supplying
#' \code{dd_max} directly is intended for sensitivity analyses (e.g., sweeping
#' compensation strength across several values) or theoretical runs with no
#' empirical calving-interval data to anchor to; a message is printed
#' whenever \code{dd_max} is supplied directly, as a reminder.
#'
#' In \code{simulate.pop()}, depletion D(t) = N_1+(t) / K_1+ is computed each year,
#' and conception rates are shifted as psi(t) = psi_K + delta_max * \code{[1 - D(t)^z]}.
#'
#' ## Output Structure
#'
#' The returned list bundles all parameters needed by \code{simulate.pop()},
#' including:
#' \itemize{
#'   \item \strong{Calibration results}: \code{s0} (or \code{theta_shift}),
#'     \code{s0_leslie}, \code{final_N}, \code{years_simulated}.
#'   \item \strong{Life history}: \code{max_age}, \code{survival}, \code{litter_size}.
#'   \item \strong{Demographics}: \code{pop_size}, \code{female_fraction},
#'     \code{maturity_age}, \code{infertility}.
#'   \item \strong{Breeding}: \code{psi_nurse}, \code{psi_rest},
#'     \code{num_mates}, \code{male_behavior}, \code{max_females},
#'     \code{weaning_age}.
#'   \item \strong{Social}: \code{pod_size}, \code{superpod_size},
#'     \code{stickiness_year}, \code{pod_to_sp}.
#'   \item \strong{DD (if active)}: \code{density_dependence = TRUE},
#'     \code{theta_shift}, \code{theta_shift_leslie}, \code{psi_nurse_K},
#'     \code{psi_rest_K}, \code{z_pt}, \code{dd_max}, \code{K_1plus},
#'     \code{target_interval}, \code{target_depletion} (if anchored).
#' }
#'
#' @return A list (returned invisibly) with calibration results and all shared
#'   parameters needed by \code{simulate.pop()}.
#'
#' @references
#' Caswell, H. (2001). Matrix Population Models: Construction, Analysis, and
#' Interpretation (2nd ed.). Sinauer Associates.
#'
#' Pella, J. J., & Tomlinson, P. K. (1973). A generalized stock production
#' model. Inter-American Tropical Tuna Commission Bulletin, 13, 422-458.
#'
#' Barlow, J., & Boveng, P. (1991). Modeling age-specific mortality for marine
#' mammal populations. Marine Mammal Science, 7(1), 50-65.
#'
#' @importFrom data.table data.table set rbindlist
#' @importFrom stats runif rpois uniroot qlogis plogis
#' @export
create.stable.pop <- function(max_age,
                              survival,
                              pop_size,
                              maturity_age,
                              litter_size,
                              psi_nurse       = 0.1,
                              psi_rest        = 0.7,
                              num_mates       = 1L,
                              female_fraction = 0.5,
                              infertility     = 0,
                              pod_size        = NULL,
                              superpod_size   = NULL,
                              stickiness_year = NULL,
                              male_behavior   = "random",
                              max_females     = NULL,
                              weaning_age     = NULL,
                              density_dependence = FALSE,
                              z_pt            = 2.39,
                              dd_max          = NULL,
                              target_interval  = NULL,
                              target_depletion = 0.3,
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

  # -- Breeding cycle probabilities --
  if (!is.numeric(psi_nurse) || length(psi_nurse) != 1L ||
      psi_nurse < 0 || psi_nurse > 1)
    stop("`psi_nurse` must be a single number between 0 and 1.")
  if (!is.numeric(psi_rest) || length(psi_rest) != 1L ||
      psi_rest < 0 || psi_rest > 1)
    stop("`psi_rest` must be a single number between 0 and 1.")
  if (psi_nurse > psi_rest)
    warning("`psi_nurse` > `psi_rest`: lactational suppression is reversed. ",
            "Usually psi_nurse <= psi_rest.")

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
    # Exactly one of dd_max or target_interval must be provided
    has_dd  <- !is.null(dd_max)
    has_ti  <- !is.null(target_interval)
    if (!has_dd && !has_ti)
      stop("When `density_dependence = TRUE`, supply either `dd_max` ",
           "(direct) or `target_interval` (anchored).")
    if (has_dd && has_ti)
      stop("Supply `dd_max` or `target_interval`, not both.")

    if (has_dd) {
      if (!is.numeric(dd_max) || length(dd_max) != 1L || dd_max <= 0)
        stop("`dd_max` must be a positive number.")
      message(
        "Note: `dd_max` was supplied directly. In most cases, `dd_max` ",
        "should be left NULL and `target_interval` (with `target_depletion`) ",
        "used instead, which anchors compensation strength to an observed ",
        "calving interval rather than an arbitrary logit-scale shift. ",
        "Supplying `dd_max` directly is intended for sensitivity analyses ",
        "or theoretical runs without empirical calving-interval data."
      )
    }
    if (has_ti) {
      if (!is.numeric(target_interval) || length(target_interval) != 1L ||
          target_interval <= 0)
        stop("`target_interval` must be a positive number (calving interval in years).")
      # Structural floor: pregnancy (1 yr) + at least 1 yr dependency
      if (target_interval < 2)
        stop("`target_interval` must be >= 2 (minimum: 1 year pregnancy + 1 year dependency).")
      if (!is.numeric(target_depletion) || length(target_depletion) != 1L ||
          target_depletion <= 0 || target_depletion >= 1)
        stop("`target_depletion` must be in (0, 1), exclusive.")
    }

    if (!is.numeric(z_pt) || length(z_pt) != 1L || z_pt <= 0)
      stop("`z_pt` must be a positive number.")
    if (psi_nurse == 0 || psi_rest == 0)
      stop("When `density_dependence = TRUE`, `psi_nurse` and `psi_rest` must be > 0 ",
           "(zero cannot be shifted on the logit scale).")
    if (psi_nurse == 1 || psi_rest == 1)
      stop("When `density_dependence = TRUE`, `psi_nurse` and `psi_rest` must be < 1 ",
           "(1.0 cannot be shifted on the logit scale).")
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

    message(sprintf(
      "Breeding cycle: psi_nurse=%.3f, psi_rest=%.3f -> %.1f%% breeding/yr, avg interval=%.1f yrs",
      psi_nurse, psi_rest, pi_1 * 100, 1 / pi_1
    ))
    message(sprintf(
      "Phase 1 -- Leslie estimate: s0 = %.4f  (lambda = %.6f)",
      s0_leslie, mat_lambda(A)
    ))

  } else {
    # -- Mode B: solve for theta_shift (density dependence) --
    # Survival is fixed (user-supplied, including s0). Find the logit-scale
    # offset on psi_nurse/psi_rest that makes lambda(K) = 1.
    logit_psi_nurse_base <- qlogis(psi_nurse)
    logit_psi_rest_base  <- qlogis(psi_rest)

    A <- matrix(0, nrow = n_classes, ncol = n_classes)
    for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i]

    dd_leslie_fn <- function(theta_shift) {
      pn <- plogis(logit_psi_nurse_base + theta_shift)
      pr <- plogis(logit_psi_rest_base  + theta_shift)
      pi_1 <- breeding_stationary(pn, pr, survival, wa_breed)[1]
      ff   <- litter_size * female_fraction * pi_1 * (1 - infertility_f)
      A[1, ] <- ogive_f_leslie * ff
      mat_lambda(A) - 1
    }

    # Check feasibility
    if (dd_leslie_fn(-10) * dd_leslie_fn(10) > 0) {
      stop(
        "No theta_shift in [-10, 10] can stabilise this population. ",
        sprintf("Lambda-1 = %.4f at theta=-10 and %.4f at theta=+10. ",
                dd_leslie_fn(-10), dd_leslie_fn(10)),
        "Check that the survival curve and breeding parameters are plausible."
      )
    }

    theta_shift_leslie <- uniroot(dd_leslie_fn, interval = c(-10, 10))$root

    # Compute at-K conception rates and stationary distribution
    psi_nurse_eff <- plogis(logit_psi_nurse_base + theta_shift_leslie)
    psi_rest_eff  <- plogis(logit_psi_rest_base  + theta_shift_leslie)

    pi_stat <- breeding_stationary(psi_nurse_eff, psi_rest_eff, survival, wa_breed)
    pi_1    <- pi_stat[1]

    # Rebuild Leslie with solved fecundity for stable age distribution
    ff     <- litter_size * female_fraction * pi_1 * (1 - infertility_f)
    A[1, ] <- ogive_f_leslie * ff
    stable_A <- mat_stable(A)

    # K_1plus: 1+ component of carrying capacity (excludes age-0)
    K_1plus <- round(pop_size * (1 - stable_A[1]))

    message(sprintf(
      "Breeding cycle (user-supplied): psi_nurse=%.3f, psi_rest=%.3f", psi_nurse, psi_rest
    ))
    message(sprintf(
      "Phase 1 -- Leslie estimate: theta_shift = %.4f  (psi_nurse_K=%.4f, psi_rest_K=%.4f, lambda=%.6f)",
      theta_shift_leslie, psi_nurse_eff, psi_rest_eff, mat_lambda(A)
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
  # Mode B (DD=TRUE):  bisect theta_shift; survival is constant.

  surv_vec <- survival

  if (!density_dependence) {
    s0_low     <- max(0.01, s0_leslie - 0.10)
    s0_high    <- min(0.99, s0_leslie + 0.10)
    s0_current <- s0_leslie
  } else {
    logit_psi_nurse_base <- qlogis(psi_nurse)
    logit_psi_rest_base  <- qlogis(psi_rest)
    theta_low     <- theta_shift_leslie - 2
    theta_high    <- theta_shift_leslie + 2
    theta_current <- theta_shift_leslie
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
      psi_nurse_w <- plogis(logit_psi_nurse_base + theta_current)
      psi_rest_w  <- plogis(logit_psi_rest_base  + theta_current)
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
    psi_nurse_K <- plogis(logit_psi_nurse_base + theta_current)
    psi_rest_K  <- plogis(logit_psi_rest_base  + theta_current)

    # Solve for dd_max from target calving interval if requested
    if (!is.null(target_interval)) {
      dd_max <- solve_dd_max(psi_nurse_K, psi_rest_K, target_interval,
                             target_depletion, z_pt, survival, wa_breed)
      message(sprintf(
        "Anchored dd_max = %.4f  (target: %.2f-yr interval at D=%.2f)",
        dd_max, target_interval, target_depletion
      ))
    }

    message(sprintf(paste0(
      "\nConverged: theta_shift = %.4f",
      "  (psi_nurse_K=%.4f, psi_rest_K=%.4f, %d years simulated, final N = %s)"),
      theta_current, psi_nurse_K, psi_rest_K,
      year_counter, format(nrow(pop), big.mark = ",")
    ))

    invisible(list(
      # Calibration results
      s0              = survival[1],
      survival        = survival,
      s0_leslie       = NA_real_,
      theta_shift     = theta_current,
      theta_shift_leslie = theta_shift_leslie,
      psi_nurse_K     = psi_nurse_K,
      psi_rest_K      = psi_rest_K,
      final_N         = nrow(pop),
      years_simulated = year_counter,
      # DD parameters
      density_dependence = TRUE,
      z_pt            = z_pt,
      dd_max          = dd_max,
      target_interval  = target_interval,
      target_depletion = if (!is.null(target_interval)) target_depletion else NULL,
      K_1plus         = K_1plus,
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
  }
}
