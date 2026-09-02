#' Simulate an age-structured population and store snapshots for sampling
#'
#' Runs a forward-time, individual-based simulation with full parentage tracking
#' and a Markov breeding cycle coupled to calf survival.  The simulation runs a
#' burn-in of \code{2 * max_age} years (to flush founders and let the age
#' structure equilibrate), after which \code{num_years} additional years are run.
#' A snapshot of the full living population for each year specified by
#' \code{sample_years} will be stored and returned for use by
#' \code{sample.pop()}.  This function is meant to be run after
#' \code{create.stable.pop()} and prior to \code{sample.pop()}.
#'
#' When density dependence is active (\code{sim_config$density_dependence =
#' TRUE}), conception probabilities are adjusted each year based on depletion
#' relative to carrying capacity, using a Pella--Tomlinson compensation
#' mechanism on the logit scale.
#'
#' @param sim_config List returned by \code{create.stable.pop()}.  Contains all
#'   life-history and social structure parameters plus calibration results.
#' @param num_years Integer. Simulation years to run after the burn-in.
#' @param sample_years Integer (scalar, vector, or NULL).
#'   \itemize{
#'     \item Scalar: store snapshots for the last \code{sample_years} years of
#'       the simulation.
#'     \item Vector: store snapshots at exactly these year indices (1-indexed,
#'       counting from the start of burn-in).
#'     \item NULL: no snapshots (but then why are you running this?).
#'   }
#'
#' @details
#' **Markov breeding cycle with calf-survival coupling:** Each mature female
#' carries a breeding state: S1 (pregnant), S2 (with dependent calf), or
#' S3 (resting).  The S2 state can last up to \code{weaning_age} years (default
#' 1 if NULL).  Each year in S2, the mother's conception probability depends on
#' whether her individual calf survived:
#' \itemize{
#'   \item Calf alive: conception probability = \code{psi_nurse} (suppressed).
#'   \item Calf dead: conception probability = \code{psi_rest} (released).
#' }
#' This creates a natural compensatory feedback: when calf mortality is high,
#' mothers breed sooner, partially offsetting the loss.
#'
#' @return A named list (returned invisibly):
#' \describe{
#'   \item{pop_summary}{data.table: sex-specific numbers-at-age for all simulation years, including burn-in.}
#'   \item{snapshots}{Named list of data.tables containing metadata for every simulated individual in each year specified by \code{sample_years}.}
#'   \item{pod_to_sp}{Integer vector mapping pod -> superpod. This is used in \code{sample.pop} to shuffle animals around between sets and trips.}
#'   \item{sim_config}{Passed through for \code{sample.pop()}.}
#' }
#'
#' @importFrom data.table data.table set rbindlist copy
#' @importFrom stats runif rpois
#' @exportS3Method NULL
#' @rawNamespace export(simulate.pop)
simulate.pop <- function(sim_config,
                         num_years,
                         sample_years = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # INPUT VALIDATION
  # ═══════════════════════════════════════════════════════════════════════════

  if (!is.list(sim_config))
    stop("`sim_config` must be a list (typically the output of create.stable.pop()).")

  required_fields <- c("max_age", "survival", "pop_size", "maturity_age",
                        "litter_size", "psi_nurse", "psi_rest", "num_mates",
                        "female_fraction", "infertility", "density_dependence")
  missing_fields <- setdiff(required_fields, names(sim_config))
  if (length(missing_fields) > 0L)
    stop("`sim_config` is missing required fields: ",
         paste(missing_fields, collapse = ", "),
         ". Did you pass the output of create.stable.pop()?")

  if (!is.numeric(num_years) || length(num_years) != 1L ||
      num_years < 1 || num_years != round(num_years))
    stop("`num_years` must be a positive integer.")

  if (!is.null(sample_years) && !is.numeric(sample_years))
    stop("`sample_years` must be NULL, a positive integer, or an integer vector.")

  # ═══════════════════════════════════════════════════════════════════════════
  # EXTRACT PARAMETERS FROM sim_config
  # ═══════════════════════════════════════════════════════════════════════════

  max_age         <- sim_config$max_age
  survival        <- sim_config$survival
  pop_size        <- sim_config$pop_size
  maturity_age    <- sim_config$maturity_age
  litter_size     <- sim_config$litter_size
  psi_nurse       <- sim_config$psi_nurse
  psi_rest        <- sim_config$psi_rest
  num_mates       <- sim_config$num_mates
  female_fraction <- sim_config$female_fraction
  pod_size_target <- sim_config$pod_size
  superpod_size   <- sim_config$superpod_size
  stickiness_year <- sim_config$stickiness_year
  male_behavior   <- sim_config$male_behavior
  max_females     <- sim_config$max_females
  weaning_age     <- sim_config$weaning_age
  infertility     <- sim_config$infertility

  # ── Density dependence parameters ──
  use_dd <- isTRUE(sim_config$density_dependence)
  if (use_dd) {
    psi_nurse_K      <- sim_config$psi_nurse_K
    psi_rest_K       <- sim_config$psi_rest_K
    logit_psi_nurse_K <- qlogis(psi_nurse_K)
    logit_psi_rest_K  <- qlogis(psi_rest_K)
    z_pt             <- sim_config$z_pt
    dd_max           <- sim_config$dd_max
    K_1plus          <- sim_config$K_1plus
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SETUP
  # ═══════════════════════════════════════════════════════════════════════════

  use_pods    <- !is.null(pod_size_target)
  use_weaning <- !is.null(weaning_age)
  wa          <- if (use_weaning) weaning_age else 0L

  # Effective weaning age for the breeding cycle (S2 duration)
  wa_breed <- if (is.null(weaning_age)) 1L else as.integer(weaning_age)

  # When DD is active, psi_nurse_K/psi_rest_K are the at-K rates.
  # When DD is off, use the raw psi_nurse/psi_rest values.
  psi_nurse_init <- if (use_dd) psi_nurse_K else psi_nurse
  psi_rest_init  <- if (use_dd) psi_rest_K  else psi_rest

  # Total simulation = burn-in (2 × max_age) + post-burn-in (num_years)
  burn_in     <- 2L * max_age
  total_years <- burn_in + num_years

  # ── Parse sample_years ──
  if (is.null(sample_years)) {
    s_years <- integer(0)
  } else if (length(sample_years) == 1L) {
    s_years <- seq.int(total_years - sample_years + 1L, total_years)
  } else {
    s_years <- as.integer(sample_years)
  }
  do_snapshots <- length(s_years) > 0L

  if (do_snapshots && min(s_years) <= burn_in) {
    warning(sprintf(
      "Snapshot year %d is within the %d-year burn-in. Population may not be at equilibrium.",
      min(s_years), burn_in
    ))
  }

  surv_vec <- survival

  # ── Parse sex-specific stickiness_year ──
  if (!is.null(stickiness_year)) {
    if (length(stickiness_year) == 1L) {
      stick_yr_F <- stick_yr_M <- stickiness_year
    } else {
      stick_yr_F <- stickiness_year[1]
      stick_yr_M <- stickiness_year[2]
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PARSE MATURITY SPECIFICATION
  # ═══════════════════════════════════════════════════════════════════════════

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

  # ═══════════════════════════════════════════════════════════════════════════
  # PARSE INFERTILITY
  # ═══════════════════════════════════════════════════════════════════════════

  if (length(infertility) == 1L) {
    infertility_f <- infertility_m <- infertility
  } else {
    infertility_f <- infertility[1]
    infertility_m <- infertility[2]
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # BREEDING CYCLE — STATIONARY DISTRIBUTION
  # ═══════════════════════════════════════════════════════════════════════════
  # Same breeding_stationary helper as in create.stable.pop(): builds the
  # (wa_breed + 2)-state Markov chain and solves for the stationary dist.

  breeding_stationary <- function(psi_n, psi_r, sv, w_b) {
    n_st <- w_b + 2L
    P <- matrix(0, nrow = n_st, ncol = n_st)
    P[1, 2] <- 1
    for (k in seq_len(w_b)) {
      row   <- k + 1L
      ell_k <- sv[k]
      psi_k <- ell_k * psi_n + (1 - ell_k) * psi_r
      P[row, 1] <- psi_k
      if (k < w_b) {
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

  pi_stat <- breeding_stationary(psi_nurse_init, psi_rest_init, surv_vec, wa_breed)
  pi_1    <- pi_stat[1]
  n_breed_states <- length(pi_stat)

  # ═══════════════════════════════════════════════════════════════════════════
  # STABLE AGE DISTRIBUTION (Leslie matrix)
  # ═══════════════════════════════════════════════════════════════════════════

  n_classes <- max_age + 1L
  A <- matrix(0, nrow = n_classes, ncol = n_classes)

  ff <- litter_size * female_fraction * pi_1 * (1 - infertility_f)
  ogive_f_leslie <- c(0, ogive_f[seq_len(max_age)])
  f_vec <- ogive_f_leslie * ff

  A[1, ] <- f_vec
  for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i]

  w_eig    <- Mod(eigen(A)$vectors[, 1])
  stable_A <- w_eig / sum(w_eig)

  # ═══════════════════════════════════════════════════════════════════════════
  # INITIALISE POPULATION
  # ═══════════════════════════════════════════════════════════════════════════

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

  # Assign breeding states from the full stationary distribution
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
    # Founder S2 mothers: calf_id = 0 (no tracked calf)
  }

  pop <- data.table(
    id          = seq_len(n_init),
    birth_year  = 0L,
    age         = init_ages,
    sex         = init_sex,
    mat_age     = init_mat_age,
    mother_id   = 0L,
    father_id   = 0L,
    breed_state = init_breed_state,
    fertile     = init_fertile,
    population  = 1L,
    calf_id     = init_calf_id,
    s2_year     = init_s2_year
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

  depletion_vec <- if (use_dd) numeric(total_years) else NULL

  bull_registry <- NULL
  if (use_pods && !is.null(male_behavior) && male_behavior == "strong_bull") {
    bull_registry <- rep(0L, n_sp)
  }

  # Columns to drop from snapshots (internal tracking only)
  internal_cols <- c("s2_year")

  message(sprintf(
    "Starting simulation: %d years burn-in (2 x max_age) + %d years = %d total  (N0 = %s)",
    burn_in, num_years, total_years, format(n_init, big.mark = ",")
  ))

  # ═══════════════════════════════════════════════════════════════════════════
  # MAIN SIMULATION LOOP
  # ═══════════════════════════════════════════════════════════════════════════
  # Each iteration = one year.  The order of operations:
  #   1. Survival (stochastic; each individual survives independently)
  #   2. Aging (deterministic; survivors age by one year)
  #   3. Between-year superpod reshuffling + cow-calf following
  #   4. Density-dependent conception adjustment (if DD active)
  #   5. Markov breeding state transitions (calf-survival-dependent)
  #   6. Create offspring (with full parentage tracking)
  #   7. Store snapshot (if this year is a snapshot year)
  #   8. Record population summary statistics

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
            # Note: sample(pool, n) would misbehave if pool has length 1 (R
            # reinterprets a length-1 numeric x as the range 1:x). Indexing
            # by position avoids this.
            new_pods[mask] <- pool[sample.int(length(pool), length(mask), replace = TRUE)]
          }
          set(pop, i = movers, j = "pod",      value = new_pods)
          set(pop, i = movers, j = "superpod", value = pod_to_sp[new_pods])
        }
      }

      # Cow-calf following: dependent calves follow mother's pod/superpod
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

    # ─── 4. Density-dependent conception adjustment ─────────────────────
    if (use_dd) {
      N_1plus <- sum(pop$age >= 1L)
      D_t     <- N_1plus / K_1plus
      delta_t <- dd_max * (1 - D_t^z_pt)
      psi_nurse_yr <- plogis(logit_psi_nurse_K + delta_t)
      psi_rest_yr  <- plogis(logit_psi_rest_K  + delta_t)
      depletion_vec[yr] <- D_t
    } else {
      psi_nurse_yr <- psi_nurse
      psi_rest_yr  <- psi_rest
    }

    # ─── 5. Markov breeding state transitions ────────────────────────────
    # IMPORTANT: All state groups are identified BEFORE any transitions are
    # applied.  This prevents double transitions within a single year.

    # Newly mature, fertile females enter at S3 (resting)
    new_mature <- which(pop$sex == "F" & pop$age == pop$mat_age &
                          is.na(pop$breed_state) & pop$fertile)
    if (length(new_mature) > 0L) {
      set(pop, i = new_mature, j = "breed_state", value = 3L)
    }

    # Snapshot current states before any transitions
    mother_rows <- which(pop$breed_state == 1L)   # S1: pregnant → will give birth
    s2_idx      <- which(pop$breed_state == 2L)   # S2: with calf
    s3_idx      <- which(pop$breed_state == 3L)   # S3: resting

    # ── S2 transitions: calf-survival-dependent ──
    # For each S2 mother, check if her calf is still alive. The calf's fate
    # determines whether the mother uses psi_nurse (suppressed) or psi_rest
    # (released from lactational suppression).
    if (length(s2_idx) > 0L) {
      calf_ids <- pop$calf_id[s2_idx]

      # Check calf survival by looking up calf_id in the current population
      calf_rows  <- match(calf_ids, pop$id)
      calf_alive <- !is.na(calf_rows)

      # Founder S2 mothers (calf_id = 0): no tracked calf, use stochastic
      # calf survival based on s2_year as proxy for calf age
      founder_s2 <- calf_ids == 0L
      if (any(founder_s2)) {
        k_founder <- pop$s2_year[s2_idx[founder_s2]]
        calf_alive[founder_s2] <- runif(sum(founder_s2)) < surv_vec[k_founder]
      }

      k <- pop$s2_year[s2_idx]

      # Conception probability depends on calf fate
      psi_eff <- ifelse(calf_alive, psi_nurse_yr, psi_rest_yr)
      conceive <- runif(length(s2_idx)) < psi_eff

      # Determine new state:
      #   Conceive → S1 (pregnant)
      #   Not conceive, calf alive, k < wa_breed → stay S2 (calf still dependent)
      #   Not conceive, calf alive, k >= wa_breed → S3 (calf weaned)
      #   Not conceive, calf dead → S3 (released)
      new_state <- rep(3L, length(s2_idx))
      new_state[conceive] <- 1L
      stay_s2 <- !conceive & calf_alive & k < wa_breed
      new_state[stay_s2] <- 2L

      set(pop, i = s2_idx, j = "breed_state", value = new_state)

      # Update s2_year: increment for stayers, clear for leavers
      new_s2y <- rep(NA_integer_, length(s2_idx))
      new_s2y[stay_s2] <- k[stay_s2] + 1L
      set(pop, i = s2_idx, j = "s2_year", value = new_s2y)

      # Clear calf_id for mothers leaving S2
      leaving_s2 <- which(new_state != 2L)
      if (length(leaving_s2) > 0L) {
        set(pop, i = s2_idx[leaving_s2], j = "calf_id", value = 0L)
      }
    }

    # ── S3 transitions ──
    if (length(s3_idx) > 0L) {
      new_state_s3 <- ifelse(runif(length(s3_idx)) < psi_rest_yr, 1L, 3L)
      set(pop, i = s3_idx, j = "breed_state", value = new_state_s3)
    }

    # ── S1 → S2: pregnant females give birth (deterministic) ──
    if (length(mother_rows) > 0L) {
      set(pop, i = mother_rows, j = "breed_state", value = 2L)
      set(pop, i = mother_rows, j = "s2_year",     value = 1L)
      # calf_id will be set below after offspring are created
    }

    # ─── 6. Create offspring ─────────────────────────────────────────────
    mature_male_mask <- pop$sex == "M" & pop$age >= pop$mat_age & pop$fertile
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
          # ── Strong bull mode ──
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
              # Indexing by position avoids sample()'s length-1 numeric footgun.
              .(bull_id = id[sample.int(.N, 1L)]),
              by = superpod
            ]
            if (nrow(candidates) > 0L) {
              bull_registry[candidates$superpod] <- candidates$bull_id
            }
          }

          for (sp in unique(mother_superpods)) {
            sp_mask <- which(mother_superpods == sp)
            bull_id <- bull_registry[sp]
            if (bull_id == 0L)
              bull_id <- all_mature_ids[sample.int(length(all_mature_ids), 1L)]
            father_mat[sp_mask, ] <- bull_id
          }

        } else {
          # ── Random mating mode ──
          father_by_sp <- split(all_mature_ids, all_mature_sps)
          for (sp in unique(mother_superpods)) {
            sp_mask   <- which(mother_superpods == sp)
            n_sp_moms <- length(sp_mask)
            pool      <- father_by_sp[[as.character(sp)]]
            if (is.null(pool) || length(pool) == 0L) pool <- all_mature_ids
            father_mat[sp_mask, ] <-
              pool[sample.int(length(pool), n_sp_moms * max_nm, replace = TRUE)]
          }
        }

      } else {
        # ── No pod structure: global random mating ──
        all_father_ids <- pop$id[mature_male_mask]
        father_mat[]   <-
          all_father_ids[sample.int(length(all_father_ids), n_mothers * max_nm,
                                     replace = TRUE)]
      }

      # ── Assign offspring to parents ──
      n_yoy <- sum(litter_sizes)
      off_mother_idx <- rep(seq_len(n_mothers), times = litter_sizes)
      off_n_mates    <- n_mates_vec[off_mother_idx]
      mate_col       <- as.integer(ceiling(runif(n_yoy) * off_n_mates))

      off_mother_id <- pop$id[mother_rows[off_mother_idx]]
      off_father_id <- father_mat[cbind(off_mother_idx, mate_col)]

      # ── Enforce max_females cap ──
      if (!is.null(max_females)) {
        fid_tab  <- table(off_father_id)
        over_ids <- as.integer(names(fid_tab[fid_tab > max_females]))

        if (length(over_ids) > 0L) {
          if (use_pods) {
            avail_by_sp <- split(all_mature_ids, all_mature_sps)
          }
          for (fid in over_ids) {
            idx     <- which(off_father_id == fid)
            keep    <- sample(idx, max_females)
            redo    <- setdiff(idx, keep)
            for (ri in redo) {
              if (use_pods) {
                mom_sp <- pop$superpod[mother_rows[off_mother_idx[ri]]]
                pool   <- avail_by_sp[[as.character(mom_sp)]]
                pool   <- pool[pool != fid]
                if (is.null(pool) || length(pool) == 0L)
                  pool <- all_mature_ids[all_mature_ids != fid]
              } else {
                pool <- all_father_ids[all_father_ids != fid]
              }
              if (length(pool) > 0L)
                off_father_id[ri] <- pool[sample.int(length(pool), 1L)]
            }
          }
        }
      }

      # ── Assign sex, maturity age, and fertility to newborns ──
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

      yoy <- data.table(
        id          = yoy_ids,
        birth_year  = as.integer(yr),
        age         = 0L,
        sex         = yoy_sex,
        mat_age     = yoy_mat_age,
        mother_id   = off_mother_id,
        father_id   = off_father_id,
        breed_state = NA_integer_,
        fertile     = yoy_fertile,
        population  = 1L,
        calf_id     = 0L,
        s2_year     = NA_integer_
      )

      # Offspring inherit their mother's pod and superpod
      if (use_pods) {
        off_pods <- pop$pod[mother_rows[off_mother_idx]]
        set(yoy, j = "pod",      value = off_pods)
        set(yoy, j = "superpod", value = pod_to_sp[off_pods])
      }

      # Set calf_id on mothers to their first offspring (the one that drives
      # the breeding cycle). For litter_size > 1, the first offspring is the
      # "dependent" calf.
      first_yoy_per_mother <- yoy_ids[!duplicated(off_mother_idx)]
      set(pop, i = mother_rows, j = "calf_id", value = first_yoy_per_mother)

      next_id <- next_id + n_yoy
      pop <- rbindlist(list(pop, yoy), use.names = TRUE)

    } # end breeding

    # ─── 7. Snapshot ─────────────────────────────────────────────────────
    # Store a deep copy, dropping internal columns
    if (do_snapshots && yr %in% s_years) {
      snap_idx <- snap_idx + 1L
      snap_pop <- copy(pop)
      # Drop internal tracking columns
      for (col in internal_cols) {
        if (col %in% names(snap_pop)) set(snap_pop, j = col, value = NULL)
      }
      snapshots[[snap_idx]] <- snap_pop
      snap_names <- c(snap_names, as.character(yr))
    }

    # ─── 8. Record population metrics ────────────────────────────────────
    yr_counts <- pop[, .N, by = .(sex, age)]
    set(yr_counts, j = "year", value = as.integer(yr))
    pop_counts[[yr]] <- yr_counts

    if (yr == 1L || yr == burn_in || yr == total_years || yr %% 10L == 0L) {
      phase <- if (yr <= burn_in) "burn-in" else "post-burn-in"
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
    total_years, burn_in, num_years,
    format(nrow(pop), big.mark = ","), length(snapshots)
  ))

  out <- list(
    pop_summary = pop_summary,
    snapshots   = snapshots,
    pod_to_sp   = pod_to_sp,
    sim_config  = sim_config
  )
  if (use_dd) out$depletion <- depletion_vec
  invisible(out)
}
