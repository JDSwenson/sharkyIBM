#' Find age-0 survival (s0) that stabilises an age-structured population
#'
#' Uses a two-phase approach:
#'
#' 1. **Leslie matrix** — build a deterministic Leslie matrix from the supplied
#'    life-history parameters and solve analytically for the s0 that gives
#'    \eqn{\lambda = 1}.
#' 2. **Iterative simulation** — run the stochastic individual-based simulation
#'    forward in time, adjusting s0 via bisection until convergence.
#'
#' @param max_age Integer. Maximum age in the population.
#' @param survival Numeric vector of length \code{max_age + 1}. Annual survival
#'   probabilities for ages 0 through \code{max_age}. The age-0 value is a
#'   placeholder; it will be replaced by the optimised estimate.
#' @param pop_size Integer. Total starting population size.
#' @param maturity_age Integer. Age at first reproduction (knife-edged).
#' @param litter_size Numeric. Mean litter size per breeding female.
#' @param psi_2 Numeric in \[0, 1\]. Probability that a female with a dependent
#'   calf (breeding state S2) becomes pregnant again the following year.
#'   Low values (~0.05–0.15) mean females rarely conceive while nursing.
#' @param psi_3 Numeric in \[0, 1\]. Probability that a resting female
#'   (breeding state S3) becomes pregnant the following year.  Controls how
#'   quickly females re-enter the breeding cycle after weaning.
#' @param num_mates Integer vector or scalar. Number of mates per female.
#'   Not used in \code{calculate.s0} but captured for pass-through.
#' @param female_fraction Numeric in (0, 1). Proportion female offspring.
#' @param infertility Reserved for future use.
#' @param pod_size Integer or NULL. Target individuals per pod.
#' @param superpod_size Integer or NULL. Pods per superpod.
#' @param stickiness_year Numeric (scalar or length-2) or NULL.
#'   Between-year superpod fidelity. Sex-specific via \code{c(female, male)}.
#' @param male_behavior Character or NULL. \code{"random"} or
#'   \code{"strong_bull"}.
#' @param weaning_age Integer or NULL. Age at independence from mother.
#' @param check_interval Integer. Years between growth assessments. Default 5.
#' @param growth_tol Numeric. Stability tolerance. Default 0.005.
#' @param stable_required Integer. Consecutive stable windows needed. Default 3.
#' @param max_windows Integer. Safety cap on windows. Default 100.
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
                         psi_2,
                         psi_3,
                         num_mates       = 1L,
                         female_fraction = 0.5,
                         infertility     = 0,
                         pod_size        = NULL,
                         superpod_size   = NULL,
                         stickiness_year = NULL,
                         male_behavior   = NULL,
                         weaning_age     = NULL,
                         check_interval  = 5L,
                         growth_tol      = 0.005,
                         stable_required = 3L,
                         max_windows     = 100L) {

  # ─── Internal helpers ────────────────────────────────────────────────────

  mat_lambda <- function(A) Mod(eigen(A, only.values = TRUE)$values[1])

  mat_stable <- function(A) {
    w <- Mod(eigen(A)$vectors[, 1])
    w / sum(w)
  }

  # ─── Markov breeding cycle stationary distribution ──────────────────────
  # Transition matrix Ψ (rows = from, cols = to):
  #        S1(preg)  S2(calf)  S3(rest)
  #   S1 [   0         1         0     ]
  #   S2 [  ψ2         0       1-ψ2   ]
  #   S3 [  ψ3         0       1-ψ3   ]
  #
  # Stationary distribution:
  #   π1 = π2 = ψ3 / (2·ψ3 + 1 - ψ2)
  #   π3 = (1-ψ2)·π1 / ψ3

  pi_denom <- 2 * psi_3 + 1 - psi_2
  pi_1 <- psi_3 / pi_denom          # proportion pregnant
  pi_2 <- pi_1                       # proportion with calf (= breeding)
  pi_3 <- (1 - psi_2) * pi_1 / psi_3 # proportion resting

  message(sprintf(
    "Breeding cycle: psi_2=%.2f, psi_3=%.2f -> %.1f%% breeding/yr, avg interval=%.1f yrs",
    psi_2, psi_3, pi_2 * 100, 1 / pi_2
  ))

  # ─── Phase 1: Leslie-matrix starting estimate ──────────────────────────

  n_classes <- max_age + 1L
  A <- matrix(0, nrow = n_classes, ncol = n_classes)

  # Expected female offspring per mature female per year:
  # pi_2 is the fraction of mature females that are breeding (in S2) each year
  ff    <- litter_size * female_fraction * pi_2
  f_vec <- c(rep(0, maturity_age), rep(ff, n_classes - maturity_age))

  A[1, ] <- f_vec
  for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i + 1L]

  # Feasibility check
  A_lo <- A; A_lo[2, 1] <- 0.01
  A_hi <- A; A_hi[2, 1] <- 0.99

  if ((mat_lambda(A_lo) - 1) * (mat_lambda(A_hi) - 1) > 0) {
    stop(
      "No s0 in [0.01, 0.99] can stabilise this population. ",
      sprintf("Lambda = %.4f at s0=0.01 and %.4f at s0=0.99.",
              mat_lambda(A_lo), mat_lambda(A_hi))
    )
  }

  s0_leslie <- uniroot(
    function(s0) { A[2, 1] <- s0; mat_lambda(A) - 1 },
    interval = c(0.01, 0.99)
  )$root

  A[2, 1]  <- s0_leslie
  stable_A <- mat_stable(A)

  message(sprintf(
    "Phase 1 — Leslie estimate: s0 = %.4f  (lambda = %.6f)",
    s0_leslie, mat_lambda(A)
  ))

  # ─── Phase 2: Initialise the individual-based population ────────────────

  use_pods <- !is.null(pod_size)
  wa       <- if (is.null(weaning_age)) 0L else weaning_age

  init_N    <- pmax(round(stable_A * pop_size), 0L)
  init_ages <- rep(0:max_age, times = init_N)
  n_init    <- sum(init_N)

  init_sex <- sample(c("F", "M"), n_init,
                     prob = c(female_fraction, 1 - female_fraction),
                     replace = TRUE)

  # Assign breeding state: mature females get state from stationary dist;
  # immature females and all males get NA.
  init_breed_state <- rep(NA_integer_, n_init)
  mature_f <- which(init_sex == "F" & init_ages >= maturity_age)
  if (length(mature_f) > 0L) {
    init_breed_state[mature_f] <- sample(
      1:3, length(mature_f), replace = TRUE,
      prob = c(pi_1, pi_2, pi_3)
    )
  }

  pop <- data.table(
    age         = init_ages,
    sex         = init_sex,
    breed_state = init_breed_state
  )

  # ── Pod setup ──
  pod_to_sp <- NULL

  if (use_pods) {
    n_pods <- max(1L, round(n_init / pod_size))
    n_sp   <- max(1L, ceiling(n_pods / superpod_size))

    pod_vec   <- rep(seq_len(n_pods), length.out = n_init)
    pod_to_sp <- rep(seq_len(n_sp), each = superpod_size, length.out = n_pods)

    set(pop, j = "pod",      value = pod_vec)
    set(pop, j = "superpod", value = pod_to_sp[pod_vec])
  }

  # ─── Phase 3: Iterative bisection search ────────────────────────────────

  s0_low     <- max(0.01, s0_leslie - 0.10)
  s0_high    <- min(0.99, s0_leslie + 0.10)
  s0_current <- s0_leslie

  year_counter       <- 0L
  consecutive_stable <- 0L

  # Parse sex-specific stickiness_year
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

    surv_vec    <- survival
    surv_vec[1] <- s0_current
    window_N    <- integer(check_interval)

    for (w in seq_len(check_interval)) {
      year_counter <- year_counter + 1L

      # ---- Survival ----
      rates <- surv_vec[pop$age + 1L]
      alive <- runif(nrow(pop)) <= rates
      pop   <- pop[alive]

      set(pop, j = "age", value = pop$age + 1L)
      pop <- pop[pop$age <= max_age]

      if (nrow(pop) == 0L) stop("Population went extinct during s0 search.")

      # ---- Between-year superpod reshuffling ----
      if (use_pods && !is.null(stickiness_year)) {
        elig <- which(pop$age >= wa)
        if (length(elig) > 0L) {
          elig_sex  <- pop$sex[elig]
          stay_prob <- ifelse(elig_sex == "F", stick_yr_F, stick_yr_M)
          movers    <- elig[runif(length(elig)) > stay_prob]

          if (length(movers) > 0L) {
            current_sp <- pop$superpod[movers]
            all_pods   <- unique(pop$pod)
            # Pre-compute pools of pods in OTHER superpods (once per superpod)
            other_pool <- lapply(
              split(all_pods, pod_to_sp[all_pods]),
              function(x) all_pods[!all_pods %in% x]
            )
            new_pods <- integer(length(movers))
            for (sp in unique(current_sp)) {
              mask <- which(current_sp == sp)
              pool <- other_pool[[as.character(sp)]]
              if (is.null(pool) || length(pool) == 0L) pool <- all_pods
              new_pods[mask] <- sample(pool, length(mask), replace = TRUE)
            }
            set(pop, i = movers, j = "pod",      value = new_pods)
            set(pop, i = movers, j = "superpod", value = pod_to_sp[new_pods])
          }
        }
      }

      # ---- Markov breeding state transitions ----
      # Newly mature females enter at S3 (resting)
      new_mature <- which(pop$sex == "F" & pop$age == maturity_age &
                            is.na(pop$breed_state))
      if (length(new_mature) > 0L) {
        set(pop, i = new_mature, j = "breed_state", value = 3L)
      }

      # Identify all breeding state groups BEFORE applying any transitions
      # (prevents double transitions like S2→S3→S1 in one year)
      mother_idx <- which(pop$breed_state == 1L)
      s2_idx     <- which(pop$breed_state == 2L)
      s3_idx     <- which(pop$breed_state == 3L)

      # Apply all transitions simultaneously
      if (length(s2_idx) > 0L) {
        new_state_s2 <- ifelse(runif(length(s2_idx)) < psi_2, 1L, 3L)
        set(pop, i = s2_idx, j = "breed_state", value = new_state_s2)
      }
      if (length(s3_idx) > 0L) {
        new_state_s3 <- ifelse(runif(length(s3_idx)) < psi_3, 1L, 3L)
        set(pop, i = s3_idx, j = "breed_state", value = new_state_s3)
      }
      # S1 → S2 (give birth)
      if (length(mother_idx) > 0L) {
        set(pop, i = mother_idx, j = "breed_state", value = 2L)
      }

      # ---- Breeding: create offspring ----
      has_males <- any(pop$sex == "M" & pop$age >= maturity_age)

      if (length(mother_idx) > 0L && has_males) {
        n_mothers    <- length(mother_idx)
        litter_sizes <- 1L + rpois(n_mothers, lambda = litter_size - 1)
        n_yoy        <- sum(litter_sizes)

        yoy_sex <- sample(c("F", "M"), n_yoy,
                          prob = c(female_fraction, 1 - female_fraction),
                          replace = TRUE)

        yoy <- data.table(
          age         = rep(0L, n_yoy),
          sex         = yoy_sex,
          breed_state = NA_integer_
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

    # ── Assess growth ──
    growth  <- mean(diff(log(window_N)))
    s0_used <- s0_current

    if (abs(growth) < growth_tol) {
      consecutive_stable <- consecutive_stable + 1L
      status <- sprintf("STABLE (%d/%d)", consecutive_stable, stable_required)
    } else {
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

    if (consecutive_stable >= stable_required) break

  } # end bisection loop

  # ── Final status ──
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

  survival_out    <- survival
  survival_out[1] <- s0_current

  invisible(list(
    # ── Calibration results ──
    s0              = s0_current,
    survival        = survival_out,
    s0_leslie       = s0_leslie,
    final_N         = nrow(pop),
    years_simulated = year_counter,
    # ── Shared life-history parameters ──
    max_age         = max_age,
    pop_size        = pop_size,
    maturity_age    = maturity_age,
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
    weaning_age     = weaning_age
  ))
}
