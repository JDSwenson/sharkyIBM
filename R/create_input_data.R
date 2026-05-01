#' Define parameters for population simulation
#'
#' Format and prepare input data for individual-based population simulation from user-defined input values.
#'
#' @param max_age Integer. Maximum age class (inclusive), so ages are 0:max_age.
#' @param survival Numeric vector. If stable = TRUE, length must be max_age
#'        (s1..s_max_age), and s0 is calculated. If stable = FALSE, length must
#'        be max_age + 1, including s0.
#' @param pop_size Integer vector (length = pop_number). Total size per population.
#' @param pop_number Integer. Number of populations.
#' @param mating_periodicity Integer. 1 = annual, 2 = biennial, etc.
#' @param YOY_survival Optional numeric. Overrides s0 if stable = FALSE.
#' @param stable Logical. If TRUE, solve for s0 using Euler-Lotka.
#' @param F_by_age Optional numeric vector (length max_age + 1). Direct fecundity.
#' @param maturity_age Optional integer. Used to derive fecundity if F_by_age is NULL. We do not need maturity age of males for this.
#' @param litter_size Optional numeric > 0. Mean pups per litter.
#' @param female_fraction Numeric in (0,1]. Default 0.5.
#' @param growth_params Either a single dataframe or a list of length pop_number with each list element corresponding to a dataframe with growth parameters for the population. If a single dataframe is used, then the growth parameters will be applied to each population; if a list of length pop_number is defined, then population-specific growth parameters can be defined. Each dataframe must have two rows, which define the sex-specific vonBertalanffy growth parameters for the species/population. If defined, these parameters will be used to assign a growth rate to each individual in the simulation at birth based on random draws from a multivariate normal distribution. Thereafter, growth proceeds deterministically, governed by the assigned growth rate. This dataframe, which will be formatted and saved for the simulation, must have one row for males and one for females and contain the following named columns:
#'   - L_inf: point estimate of L_inf (asymptotic maximum length)
#'   - L_inf_sd: standard deviation of L_inf
#'   - K: point estimate of K (growth coefficient)
#'   - K_sd: standard deviation surrounding K
#'   - t0: theoretical age at zero length
#'   - rho: correlation coefficient between asymptotic length (L_inf) and growth rate (K), bounded between -1 and 1. A negative value for rho implies that individuals that grow fast (large K) tend to have smaller asymptotic maximum lengths (small L_inf), while individuals that grow slower (small K) have larger asymptotic maximum lengths (large L_inf). A positive value for rho implies that individuals that grow fast (large K) also tend to reach larger maximum lengths (large L_inf). A value of 0 implies no correlation.
#'   - min_L0: minimum length-at-birth, in the same units as L_inf
#'   - max_L0: maximum length-at-birth, in the same units as L_inf
#'   - sex: to which sex do the parameters in this row apply? ("M" or "F")
#' @returns A list containing:
#'   - numbers_at_age: dataframe
#'   - survival: full survival vector including age 0 and max age (which is set to 0)
#'   - fecundity: vector of fecundity-by-age
#'   - s0: age 0 survival
#'
#'   Currently, create.pop.data takes a single vector of survival and a single vector of fecundity, which are applied to all populations. Future updates will allow for population-specific survival and fecundity.
#' @export
#'
#' @examples
#' input_pop <- create.pop.data(max_age = 10, survival = c(rep(0.8, times = 10)), pop_number = 3, pop_size = c(1000, 20000, 50000), mating_periodicity = 2, maturity_age = 5, litter_size = 7)
#' input_pop

create.pop.data <- function(
    max_age,
    survival,
    pop_number,
    pop_size,
    mating_periodicity,
    maturity_age = NULL,
    litter_size = NULL,
    YOY_survival = NULL,
    stable = TRUE,
    F_by_age = NULL,
    female_fraction = 0.5,
    infertility = 0,
    growth_params = NULL
) {

  ## --------------------------- Input checks --------------------------- ##
  stopifnot(is.numeric(max_age), max_age >= 1)
  stopifnot(length(pop_size) == pop_number)
  stopifnot(all(pop_size >= 0))
  stopifnot(infertility >= 0, infertility < 1)

  ## --------------------------- Survival input --------------------------- ##
  if (stable) {
    stopifnot(length(survival) == max_age | length(survival) == 1)
    if (length(survival) == 1)
      survival <- rep(survival, max_age)

    s1_to_smax <- as.numeric(survival)
    stopifnot(all(s1_to_smax > 0 & s1_to_smax <= 1))

  } else {
    stopifnot(length(survival) == max_age + 1 | length(survival) == 1)
    if (length(survival) == 1)
      survival <- c(YOY_survival, rep(survival, max_age))
    if (!is.null(YOY_survival)) survival[1] <- YOY_survival
  }

  ## --------------------------- Fecundity --------------------------- ##
  if (is.null(F_by_age)) {
    stopifnot(!is.null(maturity_age), !is.null(litter_size))
    stopifnot(maturity_age >= 0, maturity_age <= max_age)
    stopifnot(litter_size > 0)

    F_by_age <- numeric(max_age + 1)

    # Annualized fecundity in breeding ages
    F_by_age[(maturity_age + 1):(max_age + 1)] <-
      (litter_size * female_fraction) / mating_periodicity
  }

  ## --------------------------- Apply infertility --------------------------- ##
  # Infertility permanently removes a fraction of females from reproduction
  F_eff <- F_by_age * (1 - infertility)

  ## --------------------------- Solve for s0 --------------------------- ##
  if (stable) {

    if (max_age >= 2) {
      lx_no_s0 <- c(1, cumprod(s1_to_smax[1:(max_age - 1)]))
    } else {
      lx_no_s0 <- 1
    }

    Fa <- F_eff[2:(max_age + 1)]

    denom <- sum(lx_no_s0 * Fa)

    if (!is.finite(denom) || denom <= 0)
      stop("Euler–Lotka denominator non-positive")

    s0 <- 1 / denom

    if (s0 <= 0 || s0 > 1)
      stop(sprintf("Implied s0 = %.3f not in [0,1]", s0))

    survival_full <- c(s0, s1_to_smax)

  } else {
    survival_full <- as.numeric(survival)
  }

  ## Enforce terminal mortality
  survival_full[max_age + 1] <- 0

  ## --------------------------- Stable age distribution --------------------------- ##
  if (max_age >= 1) {
    lx <- cumprod(c(1, survival_full[1:max_age]))
  } else {
    lx <- 1
  }

  w <- lx / sum(lx)

  alloc_counts <- function(total, props) {
    raw <- total * props
    flo <- floor(raw)
    remainder <- total - sum(flo)
    if (remainder > 0) {
      idx <- order(raw - flo, decreasing = TRUE)
      flo[idx[seq_len(remainder)]] <- flo[idx[seq_len(remainder)]] + 1L
    }
    as.integer(flo)
  }

  ages <- 0:max_age
  out <- vector("list", pop_number)

  for (p in seq_len(pop_number)) {
    counts <- alloc_counts(pop_size[p], w)
    out[[p]] <- data.frame(
      population = p,
      age = ages,
      N = counts
    )
  }

  numbers_at_age <- do.call(rbind, out)

  ## --------------------------- Age-length formatting --------------------------- ##
  # Check that input requirements are met
  if(!is.null(growth_params)){
  if(class(growth_params) == "list"){
    if(length(growth_params) != pop_number){

      stop("growth_params must be a single dataframe or a list of length pop_number.")

      } else if(length(growth_params) == pop_number){

      growth_params_out <- growth_params

      }
  } else if(class(growth_params != "list")){
    if(length(growth_params) != 1){

      stop("growth_params must be a single dataframe or a list of length pop_number.")

      } else if(length(growth_params) == 1){

        # If a single dataframe is provided, copy the dataframe into a list of length pop_number
        growth_params_out <- rep(list(
          growth_params
        ), times = length(pop_number))
      }
    }
  }
  ## --------------------------- Return --------------------------- ##
  if(!is.null(growth_params)){
  list(
    numbers_at_age = numbers_at_age,
    survival = survival_full,
    fecundity = F_by_age,
    s0 = survival_full[1],
    litter_size = litter_size,
    growth_params = growth_params_out,
    infertility = infertility,
    mating_periodicity = mating_periodicity,
    female_fraction = female_fraction,
    maturity_age = maturity_age
  )
  } else{

    list(
      numbers_at_age = numbers_at_age,
      survival = survival_full,
      fecundity = F_by_age,
      s0 = survival_full[1],
      litter_size = litter_size,
      infertility = infertility,
      mating_periodicity = mating_periodicity,
      female_fraction = female_fraction,
      maturity_age = maturity_age
    )

  }
}
