#' Simulate an age-structured population and store snapshots for sampling
#'
#' Runs a forward-time, individual-based simulation with full parentage tracking
#' and a Markov breeding cycle.  The simulation runs a burn-in of \code{max_age}
#' years (to flush founders), then \code{num_years} additional years.  At
#' user-specified years, a snapshot of the full living population is stored and
#' returned for use by \code{sample.pop()}.
#'
#' @param sim_config List returned by \code{calculate.s0()}.
#' @param num_years Integer. Simulation years after the burn-in.
#' @param sample_years Integer (scalar, vector, or NULL).
#'
#' @details
#' **Markov breeding cycle:** Each mature female carries a breeding state:
#' S1 (pregnant), S2 (with dependent calf), or S3 (resting).  Each year,
#' states transition probabilistically:
#' \itemize{
#'   \item S1 → S2 (probability 1): pregnant female gives birth.
#'   \item S2 → S1 (probability \code{psi_2}): conceive while nursing.
#'   \item S2 → S3 (probability \code{1 - psi_2}): wean calf, rest.
#'   \item S3 → S1 (probability \code{psi_3}): conceive from resting.
#'   \item S3 → S3 (probability \code{1 - psi_3}): remain resting.
#' }
#' Offspring are created when a female transitions S1 → S2.  Newly mature
#' females enter the cycle at S3.
#'
#' @return A named list (returned invisibly):
#' \describe{
#'   \item{pop_summary}{data.table: \code{year, sex, age, N} for all years.}
#'   \item{snapshots}{Named list of data.tables at snapshot years.}
#'   \item{pod_to_sp}{Integer vector mapping pod → superpod.}
#'   \item{sim_config}{Passed through for \code{sample.pop()}.}
#' }
#'
#' @importFrom data.table data.table set rbindlist copy
#' @importFrom stats runif rpois
#' @export
simulate.pop <- function(sim_config,
                         num_years,
                         sample_years = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # EXTRACT PARAMETERS
  # ═══════════════════════════════════════════════════════════════════════════

  max_age         <- sim_config$max_age
  survival        <- sim_config$survival
  pop_size        <- sim_config$pop_size
  maturity_age    <- sim_config$maturity_age
  litter_size     <- sim_config$litter_size
  psi_2           <- sim_config$psi_2
  psi_3           <- sim_config$psi_3
  num_mates       <- sim_config$num_mates
  female_fraction <- sim_config$female_fraction
  pod_size_target <- sim_config$pod_size
  superpod_size   <- sim_config$superpod_size
  stickiness_year <- sim_config$stickiness_year
  male_behavior   <- sim_config$male_behavior
  weaning_age     <- sim_config$weaning_age

  # ═══════════════════════════════════════════════════════════════════════════
  # SETUP
  # ═══════════════════════════════════════════════════════════════════════════

  use_pods    <- !is.null(pod_size_target)
  use_weaning <- !is.null(weaning_age)
  wa          <- if (use_weaning) weaning_age else 0L

  total_years <- max_age + num_years

  # ── Parse sample_years ──
  if (is.null(sample_years)) {
    s_years <- integer(0)
  } else if (length(sample_years) == 1L) {
    s_years <- seq.int(total_years - sample_years + 1L, total_years)
  } else {
    s_years <- as.integer(sample_years)
  }
  do_snapshots <- length(s_years) > 0L

  if (do_snapshots && min(s_years) <= max_age) {
    warning(sprintf(
      "Snapshot year %d is within the %d-year burn-in. Founders may be present.",
      min(s_years), max_age
    ))
  }

  surv_vec <- survival

  # Parse sex-specific stickiness_year
  if (!is.null(stickiness_year)) {
    if (length(stickiness_year) == 1L) {
      stick_yr_F <- stick_yr_M <- stickiness_year
    } else {
      stick_yr_F <- stickiness_year[1]
      stick_yr_M <- stickiness_year[2]
    }
  }

  # ── Markov stationary distribution for initialisation ──
  pi_denom <- 2 * psi_3 + 1 - psi_2
  pi_1 <- psi_3 / pi_denom
  pi_2 <- pi_1
  pi_3 <- (1 - psi_2) * pi_1 / psi_3

  # ═══════════════════════════════════════════════════════════════════════════
  # STABLE AGE DISTRIBUTION (Leslie matrix)
  # ═══════════════════════════════════════════════════════════════════════════

  n_classes <- max_age + 1L
  A <- matrix(0, nrow = n_classes, ncol = n_classes)

  ff    <- litter_size * female_fraction * pi_2
  f_vec <- c(rep(0, maturity_age), rep(ff, n_classes - maturity_age))

  A[1, ] <- f_vec
  for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i + 1L]
  A[2, 1] <- survival[1]

  w        <- Mod(eigen(A)$vectors[, 1])
  stable_A <- w / sum(w)

  # ═══════════════════════════════════════════════════════════════════════════
  # INITIALISE POPULATION
  # ═══════════════════════════════════════════════════════════════════════════

  init_N    <- pmax(round(stable_A * pop_size), 0L)
  init_ages <- rep(0:max_age, times = init_N)
  n_init    <- sum(init_N)

  init_sex <- sample(c("F", "M"), n_init,
                     prob = c(female_fraction, 1 - female_fraction),
                     replace = TRUE)

  # Breeding state: mature females from stationary dist; others NA
  init_breed_state <- rep(NA_integer_, n_init)
  mature_f <- which(init_sex == "F" & init_ages >= maturity_age)
  if (length(mature_f) > 0L) {
    init_breed_state[mature_f] <- sample(
      1:3, length(mature_f), replace = TRUE,
      prob = c(pi_1, pi_2, pi_3)
    )
  }

  pop <- data.table(
    id          = seq_len(n_init),
    birth_year  = 0L,
    age         = init_ages,
    sex         = init_sex,
    mother_id   = 0L,
    father_id   = 0L,
    breed_state = init_breed_state,
    fertile     = TRUE,
    population  = 1L
  )

  # ── Pod / superpod initialisation ──
  pod_to_sp <- NULL
  n_sp      <- 0L

  if (use_pods) {
    n_pods <- max(1L, round(n_init / pod_size_target))
    n_sp   <- max(1L, ceiling(n_pods / superpod_size))

    pod_vec   <- rep(seq_len(n_pods), length.out = n_init)
    pod_to_sp <- rep(seq_len(n_sp), each = superpod_size, length.out = n_pods)

    set(pop, j = "pod",      value = pod_vec)
    set(pop, j = "superpod", value = pod_to_sp[pod_vec])
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PRE-ALLOCATE OUTPUT STORAGE
  # ═══════════════════════════════════════════════════════════════════════════

  pop_counts <- vector("list", total_years)
  snapshots  <- if (do_snapshots) vector("list", length(s_years)) else list()
  snap_names <- character(0)
  snap_idx   <- 0L
  next_id    <- n_init + 1L

  # Bull registry for strong_bull mode
  bull_registry <- NULL
  if (use_pods && !is.null(male_behavior) && male_behavior == "strong_bull") {
    bull_registry <- rep(0L, n_sp)
  }

  message(sprintf(
    "Starting simulation: %d years burn-in + %d years = %d total  (N0 = %s)",
    max_age, num_years, total_years, format(n_init, big.mark = ",")
  ))

  # ═══════════════════════════════════════════════════════════════════════════
  # MAIN SIMULATION LOOP
  # ═══════════════════════════════════════════════════════════════════════════

  for (yr in seq_len(total_years)) {

    # ─── 1 & 2. Survival + aging ─────────────────────────────────────────
    rates <- surv_vec[pop$age + 1L]
    alive <- runif(nrow(pop)) <= rates
    pop   <- pop[alive]

    set(pop, j = "age", value = pop$age + 1L)
    pop <- pop[pop$age <= max_age]

    if (nrow(pop) == 0L) stop("Population went extinct in year ", yr, ".")

    # ─── 3. Between-year superpod reshuffling + cow-calf following ───────
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
            new_pods[mask] <- sample(pool, length(mask), replace = TRUE)
          }
          set(pop, i = movers, j = "pod",      value = new_pods)
          set(pop, i = movers, j = "superpod", value = pod_to_sp[new_pods])
        }
      }

      # Cow-calf following
      if (use_weaning) {
        dep_idx <- which(pop$age < weaning_age & pop$mother_id != 0L)
        if (length(dep_idx) > 0L) {
          dep_mother_ids     <- pop$mother_id[dep_idx]
          mother_rows_in_pop <- match(dep_mother_ids, pop$id)
          has_living_mother  <- !is.na(mother_rows_in_pop)
          if (any(has_living_mother)) {
            update_idx   <- dep_idx[has_living_mother]
            mother_rows_ <- mother_rows_in_pop[has_living_mother]
            set(pop, i = update_idx, j = "pod",      value = pop$pod[mother_rows_])
            set(pop, i = update_idx, j = "superpod", value = pop$superpod[mother_rows_])
          }
        }
      }
    }

    # ─── 4. Markov breeding state transitions ────────────────────────────
    # Newly mature females enter at S3
    new_mature <- which(pop$sex == "F" & pop$age == maturity_age &
                          is.na(pop$breed_state))
    if (length(new_mature) > 0L) {
      set(pop, i = new_mature, j = "breed_state", value = 3L)
    }

    # Identify all breeding state groups BEFORE applying any transitions
    # (prevents double transitions like S2→S3→S1 in one year)
    mother_rows <- which(pop$breed_state == 1L)
    s2_idx      <- which(pop$breed_state == 2L)
    s3_idx      <- which(pop$breed_state == 3L)

    # Apply all transitions simultaneously
    if (length(s2_idx) > 0L) {
      new_state_s2 <- ifelse(runif(length(s2_idx)) < psi_2, 1L, 3L)
      set(pop, i = s2_idx, j = "breed_state", value = new_state_s2)
    }
    if (length(s3_idx) > 0L) {
      new_state_s3 <- ifelse(runif(length(s3_idx)) < psi_3, 1L, 3L)
      set(pop, i = s3_idx, j = "breed_state", value = new_state_s3)
    }
    # S1 → S2: give birth
    if (length(mother_rows) > 0L) {
      set(pop, i = mother_rows, j = "breed_state", value = 2L)
    }

    # ─── 5. Create offspring ─────────────────────────────────────────────
    mature_male_mask <- pop$sex == "M" & pop$age >= maturity_age
    has_males <- any(mature_male_mask)

    if (length(mother_rows) > 0L && has_males) {
      n_mothers <- length(mother_rows)

      n_mates_vec  <- sample(num_mates, n_mothers, replace = TRUE)
      litter_sizes <- 1L + rpois(n_mothers, lambda = litter_size - 1)
      max_nm       <- max(n_mates_vec)

      # ── Build father matrix ──
      father_mat <- matrix(NA_integer_, nrow = n_mothers, ncol = max_nm)

      if (use_pods) {
        mother_superpods <- pop$superpod[mother_rows]
        all_mature_ids   <- pop$id[mature_male_mask]
        all_mature_sps   <- pop$superpod[mature_male_mask]

        if (!is.null(male_behavior) && male_behavior == "strong_bull") {
          alive_ids <- pop$id
          bull_registry[!bull_registry %in% c(0L, alive_ids)] <- 0L

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

          for (sp in unique(mother_superpods)) {
            sp_mask <- which(mother_superpods == sp)
            bull_id <- bull_registry[sp]
            if (bull_id == 0L) bull_id <- sample(all_mature_ids, 1L)
            father_mat[sp_mask, ] <- bull_id
          }

        } else {
          father_by_sp <- split(all_mature_ids, all_mature_sps)
          for (sp in unique(mother_superpods)) {
            sp_mask   <- which(mother_superpods == sp)
            n_sp_moms <- length(sp_mask)
            pool      <- father_by_sp[[as.character(sp)]]
            if (is.null(pool) || length(pool) == 0L) pool <- all_mature_ids
            father_mat[sp_mask, ] <- sample(pool, n_sp_moms * max_nm,
                                            replace = TRUE)
          }
        }

      } else {
        all_father_ids <- pop$id[mature_male_mask]
        father_mat[]   <- sample(all_father_ids, n_mothers * max_nm,
                                 replace = TRUE)
      }

      # ── Create offspring ──
      n_yoy <- sum(litter_sizes)

      off_mother_idx <- rep(seq_len(n_mothers), times = litter_sizes)
      off_n_mates    <- n_mates_vec[off_mother_idx]
      mate_col       <- as.integer(ceiling(runif(n_yoy) * off_n_mates))

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
        breed_state = NA_integer_,
        fertile     = TRUE,
        population  = 1L
      )

      if (use_pods) {
        off_pods <- pop$pod[mother_rows[off_mother_idx]]
        set(yoy, j = "pod",      value = off_pods)
        set(yoy, j = "superpod", value = pod_to_sp[off_pods])
      }

      next_id <- next_id + n_yoy
      pop <- rbindlist(list(pop, yoy), use.names = TRUE)

    } # end breeding

    # ─── 6. Snapshot ─────────────────────────────────────────────────────
    if (do_snapshots && yr %in% s_years) {
      snap_idx <- snap_idx + 1L
      snapshots[[snap_idx]] <- copy(pop)
      snap_names <- c(snap_names, as.character(yr))
    }

    # ─── 7. Record population metrics ────────────────────────────────────
    yr_counts <- pop[, .N, by = .(sex, age)]
    set(yr_counts, j = "year", value = as.integer(yr))
    pop_counts[[yr]] <- yr_counts

    if (yr == 1L || yr == max_age || yr == total_years || yr %% 10L == 0L) {
      phase <- if (yr <= max_age) "burn-in" else "post-burn-in"
      message(sprintf(
        "  year %4d (%s)  |  N = %s", yr, phase,
        format(nrow(pop), big.mark = ",")
      ))
    }

  } # end main loop

  # ═══════════════════════════════════════════════════════════════════════════
  # COMPILE AND RETURN
  # ═══════════════════════════════════════════════════════════════════════════

  pop_summary <- rbindlist(pop_counts)
  if (do_snapshots) names(snapshots) <- snap_names

  message(sprintf(
    "Simulation complete: %d total years (%d burn-in + %d post), final N = %s, %d snapshots stored.",
    total_years, max_age, num_years,
    format(nrow(pop), big.mark = ","), length(snapshots)
  ))

  invisible(list(
    pop_summary = pop_summary,
    snapshots   = snapshots,
    pod_to_sp   = pod_to_sp,
    sim_config  = sim_config
  ))
}
