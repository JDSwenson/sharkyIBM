
#' Define parameters for population simulation
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
#' @param maturity_age Optional integer. Used to derive fecundity if F_by_age is NULL.
#' @param litter_size Optional numeric > 0. Mean pups per litter.
#' @param female_fraction Numeric in (0,1]. Default 0.5.
#'
#' @returns A list containing:
#'   - numbers_at_age: dataframe
#'   - survival: full survival vector including age 0 and setting survival of max age to 0
#'   - fecundity: fecundity vector used
#'   - s0: age 0 survival
#' @export
#'
#' @examples
#' create.pop.data(max_age = 10, survival = c(rep(0.8, times = 10)), pop_number = 3, pop_size = c(1000, 20000, 50000), mating_periodicity = 2, maturity_age = 5, litter_size = 7)
create.pop.data <- function(max_age,
                            survival,
                            pop_number,
                            pop_size,
                            mating_periodicity,
                            maturity_age = NULL,
                            litter_size = NULL,
                            YOY_survival = NULL,
                            stable = TRUE,
                            F_by_age = NULL,
                            female_fraction = 0.5) {

  ## --------------------------- Input checks --------------------------- ##
  stopifnot(is.numeric(max_age), max_age >= 1)
  stopifnot(length(pop_size) == pop_number)
  stopifnot(all(pop_size >= 0))

  ## --------------------------- Survival input --------------------------- ##
  if (stable == TRUE) {
    stopifnot(length(survival) == max_age | length(survival) == 1)
    if(length(survival) == 1){
      survival <- c(rep(survival, times = max_age))
    }
    s1_to_smax <- as.numeric(survival)
    stopifnot(all(s1_to_smax > 0 & s1_to_smax <= 1))
  } else {
    stopifnot(length(survival) == max_age + 1 | length(survival) == 1)
    if(length(survival) == 1){
      survival <- c(YOY_survival, rep(survival, times = max_age))
    }
    if (!is.null(YOY_survival)) survival[1] <- YOY_survival
    stopifnot(all(survival > 0 & survival <= 1))
  }

  ## --------------------------- Fecundity --------------------------- ##
  if (is.null(F_by_age) == TRUE) {
    stopifnot(!is.null(maturity_age), !is.null(litter_size))
    stopifnot(maturity_age >= 0, maturity_age <= max_age)
    stopifnot(litter_size > 0)

    F_by_age <- numeric(max_age + 1)

    # Ages are 0:max_age; indices 1:(max_age+1)
    F_by_age[seq.int(from = maturity_age + 1L, to = max_age + 1L)] <-
      (litter_size * female_fraction) / mating_periodicity

  } else {
    stopifnot(length(F_by_age) == max_age + 1)
    stopifnot(all(F_by_age >= 0))
  }

  ## --------------------------- Solve for s0 if stable --------------------------- ##
  if (stable == TRUE) {

    if (max_age >= 2) {
      prod_pre_a <- c(1, cumprod(s1_to_smax[1:(max_age - 1)]))
    } else {
      prod_pre_a <- 1
    }

    # F_a for a = 1..max_age; indices 2:(max_age+1)
    Fa_age1plus <- F_by_age[2:(max_age + 1)]

    denom <- sum(prod_pre_a * Fa_age1plus)

    if (!is.finite(denom) || denom <= 0) {
      stop("Degenerate fecundity/survival: Euler-Lotka denominator non-positive.")
    }

    s0 <- 1 / denom

    if (s0 <= 0 || s0 > 1) {
      stop(sprintf(
        "Unable to compute YOY survival within [0,1]. Implied s0 = %.4f.
         Increase survival or fecundity, or set stable=FALSE.",
        s0
      ))
    }

    survival_full <- c(s0, s1_to_smax)

  } else {

    survival_full <- as.numeric(survival)
  }

  ## --------------------------- Enforce s[max_age] = 0 --------------------------- ##
  survival_full[max_age + 1L] <- 0

  ## --------------------------- Cumulative survivorship (l_x) --------------------------- ##
  if (max_age >= 2) {
    prod_pre_a_full <- c(1, cumprod(survival_full[2:max_age]))
  } else {
    prod_pre_a_full <- 1
  }

  l_vec <- c(1, survival_full[1] * prod_pre_a_full)

  ## --------------------------- Stable age distribution --------------------------- ##
  w <- l_vec / sum(l_vec)

  ## Integer allocation helper
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

  ## --------------------------- Build N-at-age --------------------------- ##
  ages <- 0:max_age
  out_list <- vector("list", pop_number)

  for (p in seq_len(pop_number)) {
    counts <- alloc_counts(pop_size[p], w)
    out_list[[p]] <- data.frame(
      population = p,
      age = ages,
      N = counts,
      stringsAsFactors = FALSE
    )
  }

  numbers_at_age <- do.call(rbind, out_list)

  ## --------------------------- Return --------------------------- ##
  list(
    numbers_at_age = numbers_at_age,
    survival = survival_full,     # includes computed s0 AND enforced s[max_age]=0
    fecundity = F_by_age,
    s0 = survival_full[1]
    )
}
