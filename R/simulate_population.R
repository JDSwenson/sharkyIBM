#' Simulate an age-structured population and sample individuals for CKMR
#'
#' Runs a forward-time, individual-based simulation with full parentage tracking.
#' The simulation initialises from the stable age distribution implied by a
#' Leslie matrix, then loops through annual cycles of survival, optional pod
#' reshuffling, and breeding.  Individuals can be sampled in user-specified years,
#' producing a dataframe suitable for close-kin mark–recapture analysis.
#'
#' @param max_age Integer. Maximum age in the population.
#' @param survival Numeric vector of length \code{max_age + 1}. Annual survival
#'   probabilities for ages 0 through \code{max_age}.  Pass the calibrated
#'   vector returned by \code{calculate.s0()$survival}.
#' @param pop_size Integer. Total starting population size (all ages combined).
#' @param mating_periodicity Integer. Breeding cycle length in years (1 = annual,
#'   2 = biennial, etc.).
#' @param maturity_age Integer. Age at first reproduction (knife-edged maturity).
#' @param litter_size Numeric. Mean litter size per reproductive female.  Realised
#'   sizes are drawn as \code{1 + rpois(n, litter_size - 1)}, guaranteeing at
#'   least one offspring per breeding female.
#' @param num_mates Integer vector or scalar. Number of mates per female per
#'   breeding event, sampled with replacement from this vector each year.
#'   Controls the rate of half-siblings in the sampled pedigree: \code{num_mates = 1}
#'   produces only full siblings within a litter; larger values allow half-siblings.
#' @param num_years Integer. Total number of years to simulate.  Should be large
#'   enough that sampling years occur well after founders have died (see Details).
#' @param female_fraction Numeric in (0, 1). Proportion of offspring that are
#'   female.  Default 0.5.
#' @param infertility Numeric in \[0, 1). Probability that an individual is
#'   permanently infertile.  Reserved for future use; currently all individuals
#'   are treated as fertile.
#' @param stickiness Numeric in \[0, 1\] or NULL. Pod-fidelity probability.
#'   1 = individuals never leave their pod; 0 = individuals always move at each
#'   shuffle opportunity.  NULL disables pod dynamics entirely.
#' @param sticky_age Integer or NULL. Minimum age at which pod shuffling begins.
#' @param sticky_interval Integer or NULL. Years between pod-shuffle events.
#' @param superpod_size Integer or NULL. Number of pods per superpod.  Mating
#'   and (optionally) sampling occur within superpods.
#' @param male_behavior Character or NULL. \code{"random"} (all mature males
#'   in the superpod may sire offspring) or \code{"strong_bull"} (only the
#'   oldest male per superpod sires offspring).  NULL behaves like \code{"random"}
#'   without pod restriction.
#' @param sample_size Integer (scalar or vector).
#'   \itemize{
#'     \item If \code{sampling = "random"}: scalar giving total individuals to
#'           sample per year.
#'     \item If \code{sampling = "superpod"}: vector of per-superpod sample
#'           sizes; \code{length(sample_size)} superpods are randomly chosen.
#'   }
#' @param sample_years Integer (scalar or vector).
#'   \itemize{
#'     \item Scalar: sample the last \code{sample_years} years of the simulation.
#'     \item Vector: sample in these specific simulation years.
#'     \item NULL: no sampling (only population metrics are returned).
#'   }
#' @param sampling Character. \code{"random"} (uniform random sample from the
#'   full population) or \code{"superpod"} (sample within randomly chosen
#'   superpods).  Default \code{"random"}.
#' @param process_by Character. \code{"age"} (default).  \code{"length"} is
#'   reserved for future length-based dynamics.
#' @param growth_params Reserved for future length-based dynamics.
#' @param popstructure Character. \code{"panmictic"} (default) or
#'   \code{"structured"}.  Reserved for future multi-population support.
#' @param movement_array Reserved for future movement/dispersal support.
#'
#' @details
#' **Founders and burn-in:**
#' Individuals in the starting population are "founders" whose parents are
#' unknown (coded as \code{mother_id = 0, father_id = 0}).  Founder-parented
#' individuals persist for up to \code{max_age + maturity_age} years.  To
#' ensure that all sampled individuals have known (non-founder) parents,
#' sampling should not begin before year \code{max_age + maturity_age}.  A
#' warning is issued if \code{sample_years} violates this.
#'
#' **Pod dynamics:**
#' When \code{stickiness} is not NULL, the population is organised into pods
#' (family groups) that are aggregated into superpods (mating/sampling units).
#' Offspring inherit their mother's pod.  At intervals of
#' \code{sticky_interval} years, individuals at or above \code{sticky_age}
#' may move to a new pod with probability \code{1 - stickiness}.  If a
#' superpod loses all mature males, females in that superpod mate with males
#' from other superpods (cross-superpod fallback) to prevent reproductive
#' failure.
#'
#' @return A named list (returned invisibly) with two elements:
#' \describe{
#'   \item{samples}{A \code{data.table} of sampled individuals with columns:
#'     \code{id}, \code{birth_year}, \code{age}, \code{sex}, \code{mother_id},
#'     \code{father_id}, \code{capture_year}, \code{population}, and (if pods
#'     are active) \code{pod}, \code{superpod}.  Empty if no sampling occurs.}
#'   \item{pop_summary}{A \code{data.table} of population counts by year, age,
#'     and sex: columns \code{year}, \code{sex}, \code{age}, \code{N}.}
#' }
#'
#' @importFrom data.table data.table set rbindlist
#' @importFrom stats runif rpois
#' @export
simulate.pop <- function(max_age,
                         survival,
                         pop_size,
                         mating_periodicity,
                         maturity_age,
                         litter_size,
                         num_mates,
                         num_years,
                         female_fraction = 0.5,
                         infertility     = 0,
                         stickiness      = NULL,
                         sticky_age      = NULL,
                         sticky_interval = NULL,
                         superpod_size   = NULL,
                         male_behavior   = NULL,
                         sample_size     = NULL,
                         sample_years    = NULL,
                         sampling        = "random",
                         process_by      = "age",
                         growth_params   = NULL,
                         popstructure    = "panmictic",
                         movement_array  = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════

  # INPUT PARSING AND SETUP
  # ═══════════════════════════════════════════════════════════════════════════

  use_pods <- !is.null(stickiness)

  # ── Parse sample_years into an integer vector of simulation years ──
  # Scalar → sample the last N years; vector → specific years; NULL → none.
  if (is.null(sample_years) || is.null(sample_size)) {
    s_years <- integer(0)
  } else if (length(sample_years) == 1L) {
    s_years <- seq.int(num_years - sample_years + 1L, num_years)
  } else {
    s_years <- as.integer(sample_years)
  }
  do_sampling <- length(s_years) > 0L

  # Warn if any sampling year is early enough that founder-parented individuals
  # may still be alive.  For CKMR, every sampled individual's parents should
  # themselves be trackable (non-founder) individuals.
  min_safe_year <- max_age + maturity_age
  if (do_sampling && min(s_years) < min_safe_year) {
    warning(sprintf(
      paste0(
        "Sampling begins in year %d, but founder-parented individuals may ",
        "persist until year %d (max_age + maturity_age). Consider increasing ",
        "num_years or shifting sample_years."
      ),
      min(s_years), min_safe_year
    ))
  }

  # Build the survival lookup vector.  Entry i corresponds to age i-1, so
  # surv_vec[age + 1] gives the correct rate for any individual.
  surv_vec <- survival

  # ═══════════════════════════════════════════════════════════════════════════
  # STABLE AGE DISTRIBUTION (Leslie matrix)
  # ═══════════════════════════════════════════════════════════════════════════
  # We use the (calibrated) survival vector to construct a Leslie matrix and
  # extract its stable age distribution.  Initialising the IBM from this
  # distribution minimises transient dynamics at the start of the simulation.

  n_classes <- max_age + 1L
  A <- matrix(0, nrow = n_classes, ncol = n_classes)

  # Expected female offspring per mature female per year
  ff    <- litter_size * female_fraction / mating_periodicity
  f_vec <- c(rep(0, maturity_age), rep(ff, n_classes - maturity_age))

  A[1, ] <- f_vec
  for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i + 1L]
  A[2, 1] <- survival[1]

  # Right eigenvector of the dominant eigenvalue, normalised to sum to 1
  w        <- Mod(eigen(A)$vectors[, 1])
  stable_A <- w / sum(w)

  # ═══════════════════════════════════════════════════════════════════════════
  # INITIALISE POPULATION
  # ═══════════════════════════════════════════════════════════════════════════

  # Convert the stable proportions into integer counts at each age class.
  init_N    <- pmax(round(stable_A * pop_size), 0L)
  init_ages <- rep(0:max_age, times = init_N)
  n_init    <- sum(init_N)

  init_sex <- sample(c("F", "M"), n_init,
                     prob = c(female_fraction, 1 - female_fraction),
                     replace = TRUE)

  # The population is stored as a data.table for in-place modification and
  # fast subsetting, which is critical at million-individual scales.
  #
  # Key columns:
  #   id          – unique integer identifier, never recycled
  #   birth_year  – simulation year of birth (0 for founders)
  #   age         – current age in years
  #   sex         – "F" or "M"
  #   mother_id   – id of mother (0 = founder / unknown)
  #   father_id   – id of father (0 = founder / unknown)
  #   repro_cycle – which phase of the breeding cycle this female breeds in
  #                 (1:mating_periodicity for females, NA for males)
  #   fertile     – TRUE/FALSE (always TRUE for now; reserved for infertility)
  #   population  – population label (always 1 for now; reserved for structure)
  pop <- data.table(
    id          = seq_len(n_init),
    birth_year  = 0L,
    age         = init_ages,
    sex         = init_sex,
    mother_id   = 0L,
    father_id   = 0L,
    repro_cycle = ifelse(
      init_sex == "F",
      sample.int(mating_periodicity, n_init, replace = TRUE),
      NA_integer_
    ),
    fertile    = TRUE,
    population = 1L
  )

  # ── Pod / superpod initialisation ──────────────────────────────────────────
  # Pods are small family groups; superpods aggregate pods into mating and
  # sampling units.  The mapping from pod → superpod is fixed at initialisation
  # and stored as a pre-computed lookup vector (pod_to_sp) for O(1) access.
  pod_to_sp <- NULL
  n_sp      <- 0L

  if (use_pods) {
    # Each mature female seeds one pod; all other individuals are assigned to
    # an existing pod at random.
    mat_f_idx <- which(pop$sex == "F" & pop$age >= maturity_age)
    n_pods    <- length(mat_f_idx)
    if (n_pods == 0L) stop("No mature females in the initial population.")

    pod_vec            <- rep(NA_integer_, n_init)
    pod_vec[mat_f_idx] <- seq_len(n_pods)
    pod_vec[is.na(pod_vec)] <- sample.int(n_pods, sum(is.na(pod_vec)),
                                          replace = TRUE)

    # Pre-compute the persistent pod → superpod mapping.
    n_sp      <- ceiling(n_pods / superpod_size)
    pod_to_sp <- rep(seq_len(n_sp), each = superpod_size, length.out = n_pods)

    set(pop, j = "pod",      value = pod_vec)
    set(pop, j = "superpod", value = pod_to_sp[pod_vec])
    set(pop, j = "pod_year", value = 0L)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PRE-ALLOCATE OUTPUT STORAGE
  # ═══════════════════════════════════════════════════════════════════════════

  # Population metrics: one summary table per year, combined at the end.
  pop_counts <- vector("list", num_years)

  # Sampled individuals: pre-allocate a list with one slot per sample year.
  sample_list <- if (do_sampling) vector("list", length(s_years)) else list()
  sample_idx  <- 0L

  # Running ID counter — ensures every individual across all years gets a
  # unique id, which is essential for unambiguous parentage tracking.
  next_id <- n_init + 1L

  # Bull registry for strong_bull mode: integer vector indexed by superpod
  # number.  bull_registry[sp] holds the id of the current dominant male in
  # superpod sp (0 = vacant, needs election).
  bull_registry <- NULL
  if (use_pods && !is.null(male_behavior) && male_behavior == "strong_bull") {
    bull_registry <- rep(0L, n_sp)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # MAIN SIMULATION LOOP
  # ═══════════════════════════════════════════════════════════════════════════
  # Each iteration represents one year:
  #   1. Survival — stochastic mortality based on age-specific rates
  #   2. Aging    — all survivors age by one year; those past max_age die
  #   3. Pod shuffling (if active) — probabilistic pod reassignment
  #   4. Breeding — mate assignment, offspring creation with full parentage
  #   5. Sampling (if this is a sample year) — draw individuals, store metadata
  #   6. Metrics  — record numbers-at-age-and-sex

  for (yr in seq_len(num_years)) {

    # ─── 1 & 2. Survival + aging ─────────────────────────────────────────────
    # Look up each individual's survival probability via direct vector indexing
    # (no table join), then draw independent Bernoulli outcomes.
    rates <- surv_vec[pop$age + 1L]
    alive <- runif(nrow(pop)) <= rates
    pop   <- pop[alive]

    # Age all survivors by one year, then remove any that exceed max_age.
    set(pop, j = "age", value = pop$age + 1L)
    pop <- pop[pop$age <= max_age]

    if (nrow(pop) == 0L) stop("Population went extinct in year ", yr, ".")

    # ─── 3. Pod shuffling ────────────────────────────────────────────────────
    # At every sticky_interval years, individuals at or above sticky_age have
    # a (1 - stickiness) probability of moving to a randomly chosen pod.
    # pod_year tracks the simulation year of each individual's last shuffle
    # event, so the interval is measured individually.
    if (use_pods && !is.null(sticky_interval) && !is.null(sticky_age)) {
      elig <- which(pop$age >= sticky_age &
                      (yr - pop$pod_year) >= sticky_interval)
      if (length(elig) > 0L) {
        movers <- elig[runif(length(elig)) > stickiness]
        if (length(movers) > 0L) {
          new_pods <- sample(unique(pop$pod), length(movers), replace = TRUE)
          set(pop, i = movers, j = "pod",      value = new_pods)
          set(pop, i = movers, j = "superpod", value = pod_to_sp[new_pods])
        }
        # Reset the shuffle clock for all eligible individuals (movers and
        # stayers alike) so the next opportunity is correctly timed.
        set(pop, i = elig, j = "pod_year", value = as.integer(yr))
      }
    }

    # ─── 4. Breeding ─────────────────────────────────────────────────────────
    # Identify the breeding-cycle phase for this year.  The modular arithmetic
    # cycles through 1:mating_periodicity, so roughly 1/mating_periodicity of
    # mature females breed each year.
    cycle <- (yr %% mating_periodicity) + 1L

    mother_rows <- which(
      pop$sex == "F" &
        pop$age >= maturity_age &
        pop$repro_cycle == cycle
        # & pop$fertile  # uncomment when infertility is implemented
    )

    mature_male_mask <- pop$sex == "M" & pop$age >= maturity_age
    has_males <- any(mature_male_mask)

    if (length(mother_rows) > 0L && has_males) {
      n_mothers <- length(mother_rows)

      # Each mother independently draws her number of mates and litter size.
      n_mates_vec  <- sample(num_mates, n_mothers, replace = TRUE)
      litter_sizes <- 1L + rpois(n_mothers, lambda = litter_size - 1)
      max_nm       <- max(n_mates_vec)

      # ── Build father matrix ──
      # father_mat is an (n_mothers × max_nm) integer matrix where row i holds
      # the IDs of the males that mother i mates with.  Only the first
      # n_mates_vec[i] columns are used for mother i; the rest are ignored
      # during offspring assignment.
      #
      # We fill this matrix differently depending on pod structure and male
      # behavior, looping over superpods (not individual mothers) to keep the
      # operation fast even at large population sizes.
      father_mat <- matrix(NA_integer_, nrow = n_mothers, ncol = max_nm)

      if (use_pods) {
        mother_superpods <- pop$superpod[mother_rows]
        all_mature_ids   <- pop$id[mature_male_mask]
        all_mature_sps   <- pop$superpod[mature_male_mask]

        if (!is.null(male_behavior) && male_behavior == "strong_bull") {
          # ── Strong-bull mode ──
          # One dominant male per superpod sires all offspring in that superpod.
          # Check which existing bulls are still alive; elect new bulls for any
          # superpod whose bull has died.

          # Mark dead bulls as vacant
          alive_ids <- pop$id
          bull_registry[!bull_registry %in% c(0L, alive_ids)] <- 0L

          # Elect new bulls for vacant superpods (oldest mature male wins)
          vacant <- which(bull_registry == 0L)
          if (length(vacant) > 0L) {
            mature_males_dt <- data.table(
              id       = all_mature_ids,
              age      = pop$age[mature_male_mask],
              superpod = all_mature_sps
            )
            candidates <- mature_males_dt[
              superpod %in% vacant,
              .(bull_id = id[which.max(age)]),
              by = superpod
            ]
            if (nrow(candidates) > 0L) {
              bull_registry[candidates$superpod] <- candidates$bull_id
            }
          }

          # Fill father_mat: every column gets the bull for that mother's superpod
          for (sp in unique(mother_superpods)) {
            sp_mask <- which(mother_superpods == sp)
            bull_id <- bull_registry[sp]
            # If this superpod still has no bull (no mature males at all),
            # fall back to any living bull from another superpod.
            if (bull_id == 0L) bull_id <- sample(all_mature_ids, 1L)
            father_mat[sp_mask, ] <- bull_id
          }

        } else {
          # ── Random mating within superpod ──
          # Build a list of available father IDs per superpod for fast lookup,
          # then fill father_mat by superpod (one vectorised sample call per
          # superpod rather than per mother).
          father_by_sp <- split(all_mature_ids, all_mature_sps)

          for (sp in unique(mother_superpods)) {
            sp_mask   <- which(mother_superpods == sp)
            n_sp_moms <- length(sp_mask)
            pool      <- father_by_sp[[as.character(sp)]]
            # Cross-superpod fallback: if no mature males in this superpod,
            # allow mating with males from any superpod.
            if (is.null(pool) || length(pool) == 0L) pool <- all_mature_ids
            father_mat[sp_mask, ] <- sample(pool, n_sp_moms * max_nm,
                                            replace = TRUE)
          }
        }

      } else {
        # ── No pods: panmictic mating ──
        # All mature males form a single mating pool.
        all_father_ids   <- pop$id[mature_male_mask]
        father_mat[]     <- sample(all_father_ids, n_mothers * max_nm,
                                   replace = TRUE)
      }

      # ── Create offspring ──
      n_yoy <- sum(litter_sizes)

      # Expand mother index: one entry per offspring, pointing back to the
      # row in mother_rows that produced it.
      off_mother_idx <- rep(seq_len(n_mothers), times = litter_sizes)

      # For each offspring, pick one of its mother's mates uniformly.
      # ceiling(runif(n) * k) produces uniform integers in 1:k without the
      # overhead of per-offspring sample() calls.
      off_n_mates <- n_mates_vec[off_mother_idx]
      mate_col    <- as.integer(ceiling(runif(n_yoy) * off_n_mates))

      # Look up mother and father IDs
      off_mother_id <- pop$id[mother_rows[off_mother_idx]]
      off_father_id <- father_mat[cbind(off_mother_idx, mate_col)]

      yoy_sex <- sample(c("F", "M"), n_yoy,
                        prob = c(female_fraction, 1 - female_fraction),
                        replace = TRUE)

      yoy <- data.table(
        id          = seq.int(next_id, length.out = n_yoy),
        birth_year  = as.integer(yr),
        age         = 0L,
        sex         = yoy_sex,
        mother_id   = off_mother_id,
        father_id   = off_father_id,
        repro_cycle = ifelse(
          yoy_sex == "F",
          sample.int(mating_periodicity, n_yoy, replace = TRUE),
          NA_integer_
        ),
        fertile    = TRUE,
        population = 1L
      )

      # Offspring inherit their mother's pod and superpod membership.
      if (use_pods) {
        off_pods <- pop$pod[mother_rows[off_mother_idx]]
        set(yoy, j = "pod",      value = off_pods)
        set(yoy, j = "superpod", value = pod_to_sp[off_pods])
        set(yoy, j = "pod_year", value = as.integer(yr))
      }

      next_id <- next_id + n_yoy
      pop <- rbindlist(list(pop, yoy), use.names = TRUE)

    } # end breeding

    # ─── 5. Sampling ─────────────────────────────────────────────────────────
    if (do_sampling && yr %in% s_years) {
      sample_idx <- sample_idx + 1L

      if (sampling == "random") {
        # Uniform random sample from the entire living population.
        n_to_draw <- min(sample_size, nrow(pop))
        draw_rows <- sample.int(nrow(pop), n_to_draw)
        sampled   <- pop[draw_rows]

      } else if (sampling == "superpod") {
        # Choose random superpods, then sample within each.
        # sample_size is a vector: one element per superpod to sample.
        n_pods_to_sample <- length(sample_size)
        chosen_sps <- sample(unique(pop$superpod), n_pods_to_sample)

        sampled_list_inner <- vector("list", n_pods_to_sample)
        for (k in seq_len(n_pods_to_sample)) {
          sp_pop  <- pop[pop$superpod == chosen_sps[k]]
          n_draw  <- min(sample_size[k], nrow(sp_pop))
          sampled_list_inner[[k]] <- sp_pop[sample.int(nrow(sp_pop), n_draw)]
        }
        sampled <- rbindlist(sampled_list_inner)
      }

      set(sampled, j = "capture_year", value = as.integer(yr))
      sample_list[[sample_idx]] <- sampled
    }

    # ─── 6. Record population metrics ────────────────────────────────────────
    # Aggregate counts by sex and age for ground-truthing CKMR estimates.
    yr_counts <- pop[, .N, by = .(sex, age)]
    set(yr_counts, j = "year", value = as.integer(yr))
    pop_counts[[yr]] <- yr_counts

    # Progress message every 10 years (or first/last year)
    if (yr == 1L || yr == num_years || yr %% 10L == 0L) {
      message(sprintf(
        "  year %4d  |  N = %s", yr, format(nrow(pop), big.mark = ",")
      ))
    }

  } # end main loop

  # ═══════════════════════════════════════════════════════════════════════════
  # COMPILE AND RETURN OUTPUTS
  # ═══════════════════════════════════════════════════════════════════════════

  pop_summary <- rbindlist(pop_counts)

  if (do_sampling) {
    samples_df <- rbindlist(sample_list)
  } else {
    # Return an empty table with the correct schema so downstream code doesn't break
    samples_df <- data.table(
      id = integer(), birth_year = integer(), age = integer(),
      sex = character(), mother_id = integer(), father_id = integer(),
      capture_year = integer(), population = integer()
    )
  }

  message(sprintf(
    "Simulation complete: %d years, final N = %s, %d individuals sampled.",
    num_years, format(nrow(pop), big.mark = ","), nrow(samples_df)
  ))

  invisible(list(
    samples     = samples_df,
    pop_summary = pop_summary
  ))
}
