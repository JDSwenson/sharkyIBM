#' Find age-0 survival (s0) that stabilises an age-structured population
#'
#' Uses a two-phase approach:
#'
#' 1. **Leslie matrix** — build a deterministic Leslie matrix from the supplied
#'    life-history parameters and solve analytically for the s0 that gives
#'    \eqn{\lambda = 1}. This analytical value is a useful starting estimate
#'    but often differs from the value that stabilises the stochastic IBM
#'    because the IBM includes demographic stochasticity, a finite-population
#'    age structure, and a discrete breeding cycle.
#' 2. **Iterative simulation** — run the stochastic individual-based simulation
#'    forward in time, pausing every \code{check_interval} years to measure
#'    realised population growth and adjusting s0 via bisection. Convergence
#'    is declared after \code{stable_required} consecutive assessment windows
#'    each show near-zero growth.
#'
#' @param max_age Integer. Maximum age in the population.
#' @param survival Numeric vector of length \code{max_age + 1}. Annual survival
#'   probabilities for ages 0 through \code{max_age}. The age-0 value is a
#'   placeholder; it will be replaced by the optimised estimate.
#' @param pop_size Integer. Total starting population size (all ages combined).
#' @param mating_periodicity Integer. Breeding interval in years (1 = annual,
#'   2 = biennial, etc.).
#' @param maturity_age Integer. Age at first reproduction (knife-edged).
#' @param litter_size Numeric. Mean litter size per reproductive female.
#'   Realised sizes are drawn as \code{1 + rpois(n, litter_size - 1)},
#'   which guarantees a minimum of one offspring per breeding female.
#' @param num_mates Integer vector or scalar. Number of mates per female.
#'   Retained for signature consistency with \code{simulate.pop} but not used
#'   here: because father identity does not affect the total number of offspring
#'   per mother, mate assignment is unnecessary when only population growth rate
#'   is needed.
#' @param female_fraction Numeric in (0, 1). Proportion of offspring that are
#'   female. Default 0.5.
#' @param infertility Numeric in \[0, 1). Probability that an individual is
#'   permanently infertile. Retained for forward-compatibility; currently
#'   ignored (all individuals are treated as fertile).
#' @param stickiness Numeric in \[0, 1\] or NULL. Pod-fidelity probability.
#'   1 = individuals never leave their pod; 0 = individuals always move at each
#'   shuffle opportunity. NULL disables pod dynamics entirely.
#' @param sticky_age Integer or NULL. Minimum age at which pod shuffling applies.
#' @param sticky_interval Integer or NULL. Years between pod-shuffle events.
#' @param superpod_size Integer or NULL. Number of pods grouped into each
#'   superpod. Mating occurs within superpods.
#' @param male_behavior Character or NULL. `"random"` (any mature male in the
#'   superpod may sire offspring) or `"strong_bull"` (only the oldest male per
#'   superpod sires offspring). Not used in \code{calculate.s0} because father
#'   identity does not affect offspring count, but retained for
#'   forward-compatibility.
#' @param check_interval Integer. Years simulated between growth-rate
#'   assessments. Larger values smooth stochastic noise but slow convergence.
#'   Default 5.
#' @param growth_tol Numeric. A window is deemed stable when
#'   \eqn{|\bar{r}| < \code{growth\_tol}}, where \eqn{\bar{r}} is the mean
#'   annual log growth rate over the window. Default 0.005.
#' @param stable_required Integer. Number of consecutive stable windows required
#'   before declaring convergence. Default 3.
#' @param max_windows Integer. Safety cap on total assessment windows to prevent
#'   infinite loops. A warning is issued if this limit is reached. Default 100.
#'
#' @return A list (returned invisibly) with elements:
#' \describe{
#'   \item{s0}{Empirically optimised age-0 survival probability.}
#'   \item{survival}{Full survival vector with \code{s0} inserted at position 1.
#'     Pass this directly to \code{simulate.pop}.}
#'   \item{s0_leslie}{Leslie-matrix analytical starting estimate, for comparison.}
#'   \item{final_N}{Population size at the end of the search.}
#'   \item{years_simulated}{Total number of years simulated during the search.}
#' }
#'
#' @importFrom data.table data.table set rbindlist
#' @importFrom stats runif rpois uniroot
#' @export
calculate.s0 <- function(max_age,
                         survival,
                         pop_size,
                         mating_periodicity,
                         maturity_age,
                         litter_size,
                         num_mates       = 1L,
                         female_fraction = 0.5,
                         infertility     = 0,
                         stickiness      = NULL,
                         sticky_age      = NULL,
                         sticky_interval = NULL,
                         superpod_size   = NULL,
                         male_behavior   = NULL,
                         check_interval  = 5L,
                         growth_tol      = 0.005,
                         stable_required = 3L,
                         max_windows     = 100L) {

  # ─── Internal helper functions ─────────────────────────────────────────────

  # Returns the dominant eigenvalue (lambda) of a matrix. For a Leslie matrix,
  # this is the population growth rate: lambda > 1 means growing, < 1 declining.
  mat_lambda <- function(A) Mod(eigen(A, only.values = TRUE)$values[1])

  # Returns the stable age distribution of a Leslie matrix as a proportion vector
  # (sums to 1). This is the right eigenvector corresponding to lambda, normalised.
  # We use it to set the initial age structure of the IBM to a theoretically
  # stable configuration, which reduces the transient burn-in time needed.
  mat_stable <- function(A) {
    w <- Mod(eigen(A)$vectors[, 1])
    w / sum(w)
  }

  # ─── Phase 1: Leslie-matrix starting estimate ─────────────────────────────
  # We build a classic age-structured Leslie matrix to get an analytical estimate
  # of s0. This serves two purposes:
  #   (a) provides the starting value for the bisection search, dramatically
  #       reducing the number of simulation years needed to converge, and
  #   (b) provides the stable age distribution used to initialise the IBM.

  n_classes <- max_age + 1L
  A <- matrix(0, nrow = n_classes, ncol = n_classes)

  # Expected female offspring produced per mature female per year.
  # Dividing litter_size by mating_periodicity converts from offspring-per-litter
  # to offspring-per-year, and multiplying by female_fraction keeps only daughters.
  ff <- litter_size * female_fraction / mating_periodicity

  # Fecundity vector: zero before maturity, ff thereafter (knife-edged maturity).
  # This vector fills the first row of the Leslie matrix.
  f_vec <- c(rep(0, maturity_age), rep(ff, n_classes - maturity_age))

  # Fill the Leslie matrix:
  #   Row 1 (fecundity):  A[1, j] = expected female offspring from age-j females
  #   Sub-diagonal (survival): A[i+1, i] = probability of surviving from age i-1 to i
  # Note: survival[i+1] gives the rate for age i because survival is 1-indexed in R
  # and our age classes run from 0 (stored in position 1) to max_age (position max_age+1).
  A[1, ] <- f_vec
  for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i + 1L]

  # ── Feasibility check ──
  # Before searching for s0, confirm that at least one value in [0.01, 0.99]
  # can produce a stable population. We do this by evaluating lambda at extreme
  # s0 values and checking that one gives lambda > 1 and the other lambda < 1
  # (i.e., the root is bracketed). If both extremes produce the same sign, no
  # feasible s0 exists and we stop with an informative error.
  A_lo <- A; A_lo[2, 1] <- 0.01   # A[2,1] is the s0 position in the Leslie matrix
  A_hi <- A; A_hi[2, 1] <- 0.99

  if ((mat_lambda(A_lo) - 1) * (mat_lambda(A_hi) - 1) > 0) {
    stop(
      "No s0 in [0.01, 0.99] can stabilise this population. ",
      sprintf("Lambda = %.4f at s0=0.01 and %.4f at s0=0.99.",
              mat_lambda(A_lo), mat_lambda(A_hi))
    )
  }

  # ── Solve for s0 analytically via bisection (uniroot) ──
  # We find the s0 that makes lambda exactly equal to 1, i.e., the root of
  # f(s0) = lambda(A(s0)) - 1. uniroot() uses Brent's method, which is fast
  # and guaranteed to converge given that the root is bracketed.
  s0_leslie <- uniroot(
    function(s0) { A[2, 1] <- s0; mat_lambda(A) - 1 },
    interval = c(0.01, 0.99)
  )$root

  # Insert the analytical s0 into the Leslie matrix and extract the corresponding
  # stable age distribution. We use this distribution to set the initial age
  # structure of the IBM, which is more informative than an arbitrary flat
  # or uniform distribution.
  A[2, 1]  <- s0_leslie
  stable_A <- mat_stable(A)

  message(sprintf(
    "Phase 1 — Leslie estimate: s0 = %.4f  (lambda = %.6f)",
    s0_leslie, mat_lambda(A)
  ))

  # ─── Phase 2: Initialise the individual-based population ───────────────────
  # Flag for whether pod structure is active throughout the simulation.
  use_pods <- !is.null(stickiness)

  # Translate the stable age distribution into integer individual counts at each
  # age class. We scale by pop_size and round to the nearest integer.
  # pmax(..., 0) guards against any rounding artefacts that could produce negatives.
  init_N    <- pmax(round(stable_A * pop_size), 0L)
  init_ages <- rep(0:max_age, times = init_N)  # expand into one age per individual
  n_init    <- sum(init_N)

  # Randomly assign sex to all founding individuals.
  init_sex <- sample(c("F", "M"), n_init,
                     prob = c(female_fraction, 1 - female_fraction),
                     replace = TRUE)

  # Build the population as a data.table. We use data.table throughout the
  # simulation loop because it modifies columns in place (avoiding copies),
  # supports very fast row subsetting, and handles the large row counts expected
  # for million-individual populations.
  #
  # repro_cycle: determines which year within a breeding cycle each female breeds.
  # A female with repro_cycle == k breeds in years where (year %% mating_periodicity) + 1 == k.
  # This distributes breeding females roughly evenly across years, which is
  # important for biennial or longer breeding cycles. Males are assigned NA.
  pop <- data.table(
    age         = init_ages,
    sex         = init_sex,
    repro_cycle = ifelse(
      init_sex == "F",
      sample.int(mating_periodicity, n_init, replace = TRUE),
      NA_integer_
    )
  )

  # ── Pod setup (only if stickiness is defined) ──
  if (use_pods) {
    # Identify currently mature females — they each form the nucleus of one pod.
    # Using one pod per mature female at initialisation ensures a sensible number
    # of pods relative to population size.
    mat_f_idx <- which(pop$sex == "F" & pop$age >= maturity_age)
    n_pods    <- length(mat_f_idx)
    if (n_pods == 0L) stop("No mature females in the initial population.")

    # Assign each mature female a unique pod ID (1:n_pods), then randomly assign
    # all other individuals (juvenile males, juvenile females, and mature males)
    # to one of the existing pods.
    pod_vec                <- rep(NA_integer_, n_init)
    pod_vec[mat_f_idx]     <- seq_len(n_pods)
    pod_vec[is.na(pod_vec)] <- sample.int(n_pods, sum(is.na(pod_vec)), replace = TRUE)

    # Group pods into superpods. Superpods are the mating and sampling units.
    # We pre-compute a lookup vector (pod_to_sp) indexed by pod number so that
    # pod → superpod lookups during the loop are O(1) vector indexing rather than joins.
    n_sp      <- ceiling(n_pods / superpod_size)
    pod_to_sp <- rep(seq_len(n_sp), each = superpod_size, length.out = n_pods)

    # Add pod columns to the population table using data.table's set(), which
    # modifies the table in place without creating a copy.
    set(pop, j = "pod",      value = pod_vec)
    set(pop, j = "superpod", value = pod_to_sp[pod_vec])
    set(pop, j = "pod_year", value = 0L)  # year in which this individual last shuffled pods
  }

  # ─── Phase 3: Iterative bisection search ──────────────────────────────────
  # Strategy: run the IBM forward, assess growth every check_interval years,
  # and use bisection to narrow in on the s0 that produces zero growth.
  #
  # Bisection works by maintaining a bracket [s0_low, s0_high] that is guaranteed
  # to contain the true s0:
  #   - If the population grows (r > 0), s0 is too high → lower the ceiling.
  #   - If the population declines (r < 0), s0 is too low → raise the floor.
  #   - The next candidate is the midpoint of the updated bracket.
  #
  # We initialise the bracket around the Leslie estimate, which is usually close
  # enough that convergence is fast. The ±0.10 window covers typical discrepancies
  # between the analytical and realised s0; uniroot will extend it if needed.
  s0_low     <- max(0.01, s0_leslie - 0.10)
  s0_high    <- min(0.99, s0_leslie + 0.10)
  s0_current <- s0_leslie

  year_counter       <- 0L
  consecutive_stable <- 0L

  message(sprintf(
    "Phase 2 — searching (check every %d yrs, tol = %.4f, need %d stable)...",
    check_interval, growth_tol, stable_required
  ))

  for (win in seq_len(max_windows)) {

    # Apply the current candidate s0 by substituting it into position 1 of the
    # survival vector. Position 1 corresponds to age 0 because survival is indexed
    # 1:(max_age+1) in R for age classes 0:max_age.
    surv_vec    <- survival
    surv_vec[1] <- s0_current

    # Record total population size at the end of each year in this window.
    # We use this vector to compute the mean log growth rate after the window.
    window_N <- integer(check_interval)

    # ── Inner loop: simulate check_interval years ──
    for (w in seq_len(check_interval)) {
      year_counter <- year_counter + 1L

      # ---- Survival ----
      # Look up each individual's annual survival probability by direct vector
      # indexing: age 0 → surv_vec[1], age 1 → surv_vec[2], etc.
      # This avoids a table join and is substantially faster at large N.
      # Each individual then independently survives (Bernoulli draw via runif).
      rates <- surv_vec[pop$age + 1L]
      alive <- runif(nrow(pop)) <= rates
      pop   <- pop[alive]              # retain survivors only

      # Age everyone up by one year, then remove individuals past max_age.
      # Individuals at max_age survive and breed in the current year (via the
      # fecundity first-row of the Leslie matrix), but are then removed so that
      # the age structure matches the Leslie framework.
      set(pop, j = "age", value = pop$age + 1L)
      pop <- pop[pop$age <= max_age]

      if (nrow(pop) == 0L) stop("Population went extinct during s0 search.")

      # ---- Pod shuffling (only active if stickiness is defined) ----
      # At every sticky_interval years, individuals at or above sticky_age have
      # a (1 - stickiness) probability of switching to a new pod. Individuals
      # that switch are assigned a random pod from the current pool, and their
      # superpod membership is updated via the pre-computed pod_to_sp lookup.
      # pod_year tracks the last year each individual shuffled, so we can check
      # whether sticky_interval years have elapsed.
      if (use_pods && !is.null(sticky_interval) && !is.null(sticky_age)) {
        elig <- which(pop$age >= sticky_age &
                        (year_counter - pop$pod_year) >= sticky_interval)
        if (length(elig) > 0L) {
          # Probabilistic move: each eligible individual stays with probability
          # stickiness, moves with probability 1 - stickiness.
          movers <- elig[runif(length(elig)) > stickiness]
          if (length(movers) > 0L) {
            new_pods <- sample(unique(pop$pod), length(movers), replace = TRUE)
            set(pop, i = movers, j = "pod",      value = new_pods)
            set(pop, i = movers, j = "superpod", value = pod_to_sp[new_pods])
          }
          # Reset pod_year for all eligible individuals (movers and stayers alike),
          # so the next shuffle opportunity is correctly timed from this year.
          set(pop, i = elig, j = "pod_year", value = year_counter)
        }
      }

      # ---- Breeding ----
      # Identify which breeding-cycle value applies this year. The modular
      # arithmetic cycles through 1:mating_periodicity each year, so roughly
      # 1/mating_periodicity of mature females breed in any given year.
      cycle <- (year_counter %% mating_periodicity) + 1L

      # Select eligible mothers: mature, fertile (all for now), and in the
      # correct breeding-cycle year.
      mother_idx <- which(
        pop$sex == "F" &
          pop$age >= maturity_age &
          pop$repro_cycle == cycle
      )

      # Require at least one mature male to exist somewhere in the population
      # before allowing any breeding. For calculate.s0 we do not track individual
      # fathers, so the only question is whether males are present at all.
      has_males <- any(pop$sex == "M" & pop$age >= maturity_age)

      if (length(mother_idx) > 0L && has_males) {
        n_mothers <- length(mother_idx)

        # Draw a litter size for each mother. Using 1 + Poisson(lambda - 1)
        # guarantees a minimum litter of 1 offspring per breeding female, while
        # preserving the intended mean. (A plain Poisson(lambda) would allow
        # zero-offspring litters, which is biologically inappropriate here.)
        litter_sizes <- 1L + rpois(n_mothers, lambda = litter_size - 1)
        n_yoy        <- sum(litter_sizes)

        # Assign sex to offspring with the user-specified female fraction.
        yoy_sex <- sample(c("F", "M"), n_yoy,
                          prob = c(female_fraction, 1 - female_fraction),
                          replace = TRUE)

        # Create the young-of-year (YOY) table. We omit parentage columns here
        # because they are not needed for computing population growth rate,
        # and recording them for every individual at large population sizes
        # would slow this function substantially.
        yoy <- data.table(
          age         = rep(0L, n_yoy),
          sex         = yoy_sex,
          repro_cycle = ifelse(
            yoy_sex == "F",
            # Each female offspring is randomly assigned to one cycle phase,
            # distributing future breeding effort evenly across years.
            sample.int(mating_periodicity, n_yoy, replace = TRUE),
            NA_integer_
          )
        )

        # If pods are active, offspring inherit their mother's pod and superpod.
        # rep(..., times = litter_sizes) replicates each mother's pod ID once
        # per offspring she produced, creating the correct pod assignment for
        # the full YOY cohort.
        if (use_pods) {
          mother_pods <- pop$pod[mother_idx]
          yoy_pods    <- rep(mother_pods, times = litter_sizes)
          set(yoy, j = "pod",      value = yoy_pods)
          set(yoy, j = "superpod", value = pod_to_sp[yoy_pods])
          set(yoy, j = "pod_year", value = year_counter)
        }

        # Append YOY to the living population. rbindlist() is data.table's
        # equivalent of bind_rows() but faster; use.names = TRUE aligns columns
        # by name rather than position, which is safer when column sets differ
        # (e.g., pod columns present only when use_pods is TRUE).
        pop <- rbindlist(list(pop, yoy), use.names = TRUE)
      }

      # Record total population size at the end of this year.
      window_N[w] <- nrow(pop)

    } # end inner (annual) loop

    # ── Assess growth over the completed window ──
    # Mean log growth rate: mean of log(N[t+1]/N[t]) across consecutive year pairs.
    # This is equivalent to (log(N_end) - log(N_start)) / (check_interval - 1),
    # but the mean-of-differences is more robust to mid-window fluctuations.
    # Negative values = declining, positive = growing, near zero = stable.
    growth  <- mean(diff(log(window_N)))
    s0_used <- s0_current   # save for display before potentially updating

    if (abs(growth) < growth_tol) {
      # This window qualifies as stable. Accumulate toward the required streak.
      # We do NOT adjust s0 during stable windows; the current value is kept
      # for the next window so the streak can continue building.
      consecutive_stable <- consecutive_stable + 1L
      status <- sprintf("STABLE (%d/%d)", consecutive_stable, stable_required)

    } else {
      # Growth is outside tolerance — bisect.
      # Reset the streak counter because we are changing s0; the population
      # will need time to adjust to the new value before stability can be assessed.
      consecutive_stable <- 0L

      # Update the bisection bracket: the side that produced this candidate s0
      # is moved inward to the current value, and the new candidate is the midpoint.
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

    # ── Convergence check ──
    if (consecutive_stable >= stable_required) break

  } # end bisection (outer) loop

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

  # Build and return the complete, updated survival vector.
  # The caller can pass result$survival directly to simulate.pop().
  survival_out    <- survival
  survival_out[1] <- s0_current

  invisible(list(
    s0              = s0_current,
    survival        = survival_out,
    s0_leslie       = s0_leslie,
    final_N         = nrow(pop),
    years_simulated = year_counter
  ))
}
