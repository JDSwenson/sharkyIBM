#' Simulate an age-structured population and store snapshots for sampling
#'
#' Runs a forward-time, individual-based simulation with full parentage tracking
#' and a Markov breeding cycle.  The simulation runs a burn-in of
#' \code{2 * max_age} years (to flush founders and let the age structure
#' equilibrate), after which \code{num_years} additional years are run. A snapshot of the full living population for each year specified by \code{sample_years} will be stored and returned for use by \code{sample.pop()}. This function is meant to be run after \code{calculate.s0} and prior to \code{sample.pop}.
#'
#' @param sim_config List returned by \code{calculate.s0()}.  Contains all
#'   life-history and social structure parameters plus the calibrated survival
#'   vector.
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
#' **Markov breeding cycle:** Each mature female carries a breeding state:
#' S1 (pregnant), S2 (with dependent calf), or S3 (resting).  Each year,
#' states transition probabilistically:
#' \itemize{
#'   \item S1 -> S2 (probability 1): pregnant female gives birth.
#'   \item S2 -> S1 (probability \code{psi_2}): conceive while nursing.
#'   \item S2 -> S3 (probability \code{1 - psi_2}): wean calf, rest.
#'   \item S3 -> S1 (probability \code{psi_3}): conceive from resting.
#'   \item S3 -> S3 (probability \code{1 - psi_3}): remain resting.
#' }
#' Offspring are created when a female transitions S1 -> S2.  Newly mature
#' females enter the cycle at S3.
#'
#' @return A named list (returned invisibly):
#' \describe{
#'   \item{pop_summary}{data.table: sex-specific numbers-at-age for all simulation years, including burn-in.}
#'   \item{snapshots}{Named list of data.tables containing metadata for every simulated individual in each year specified by \code{sample.years.}
#'   \item{pod_to_sp}{Integer vector mapping pod -> superpod.} This is used in \code{sample.pop} to shuffle animals around between sets and trips.
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
  # INPUT VALIDATION
  # ═══════════════════════════════════════════════════════════════════════════

  if (!is.list(sim_config))
    stop("`sim_config` must be a list (typically the output of calculate.s0()).")

  required_fields <- c("max_age", "survival", "pop_size", "maturity_age",
                        "litter_size", "psi_2", "psi_3", "num_mates",
                        "female_fraction", "infertility")
  missing_fields <- setdiff(required_fields, names(sim_config))
  if (length(missing_fields) > 0L)
    stop("`sim_config` is missing required fields: ",
         paste(missing_fields, collapse = ", "),
         ". Did you pass the output of calculate.s0()?")

  if (!is.numeric(num_years) || length(num_years) != 1L ||
      num_years < 1 || num_years != round(num_years))
    stop("`num_years` must be a positive integer.")

  if (!is.null(sample_years) && !is.numeric(sample_years))
    stop("`sample_years` must be NULL, a positive integer, or an integer vector.")

  # ═══════════════════════════════════════════════════════════════════════════
  # EXTRACT PARAMETERS FROM sim_config
  # ═══════════════════════════════════════════════════════════════════════════
  # All life-history and social structure parameters were specified once in
  # calculate.s0() and bundled into sim_config.  We unpack them here.

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
  max_females     <- sim_config$max_females
  weaning_age     <- sim_config$weaning_age
  infertility     <- sim_config$infertility

  # ═══════════════════════════════════════════════════════════════════════════
  # SETUP
  # ═══════════════════════════════════════════════════════════════════════════

  use_pods    <- !is.null(pod_size_target)
  use_weaning <- !is.null(weaning_age)
  wa          <- if (use_weaning) weaning_age else 0L

  # Total simulation length = burn-in (2 × max_age years) + post-burn-in (num_years).
  # The burn-in serves two purposes:
  #   1. Flush all founders (mother_id = 0, father_id = 0) — requires max_age years.
  #   2. Let the IBM's age structure settle from the Leslie-derived initial
  #      distribution to its true stochastic equilibrium — requires ~2× max_age.
  # Using 2× max_age absorbs both the founder flush and the transient adjustment,
  # so snapshots are taken from a population at demographic equilibrium.
  burn_in     <- 2L * max_age
  total_years <- burn_in + num_years

  # ── Parse sample_years ──
  # If sample_years is a single integer, store the last N years.
  # If it's a vector, store snapshots at those exact year indices.
  if (is.null(sample_years)) {
    s_years <- integer(0)
  } else if (length(sample_years) == 1L) {
    s_years <- seq.int(total_years - sample_years + 1L, total_years)
  } else {
    s_years <- as.integer(sample_years)
  }
  do_snapshots <- length(s_years) > 0L

  # Warn if any snapshots fall within the burn-in period
  if (do_snapshots && min(s_years) <= burn_in) {
    warning(sprintf(
      "Snapshot year %d is within the %d-year burn-in. Population may not be at equilibrium.",
      min(s_years), burn_in
    ))
  }

  surv_vec <- survival

  # ── Parse sex-specific stickiness_year ──
  # Stickiness controls the probability that an individual stays in its current
  # superpod between years.  Can be a single value (same for both sexes) or
  # c(female, male) for sex-specific rates.
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
  # Same parsing logic as calculate.s0(): convert user input into two ogive
  # vectors (female and male), each of length max_age + 1, giving P(mature|age).

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

  # Derive PMF from ogive for sampling individual maturity ages at birth.
  # See calculate.s0() for detailed explanation of the PMF derivation.
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
  # MARKOV BREEDING CYCLE — STATIONARY DISTRIBUTION
  # ═══════════════════════════════════════════════════════════════════════════
  # Used to initialise breeding states for the starting population.
  # See calculate.s0() for the full derivation of the stationary distribution.

  pi_denom <- 2 * psi_3 + 1 - psi_2
  pi_1 <- psi_3 / pi_denom          # proportion pregnant
  pi_2 <- pi_1                       # proportion with calf (= breeding)
  pi_3 <- (1 - psi_2) * pi_1 / psi_3 # proportion resting

  # ═══════════════════════════════════════════════════════════════════════════
  # STABLE AGE DISTRIBUTION (Leslie matrix)
  # ═══════════════════════════════════════════════════════════════════════════
  # Build the same Leslie matrix as calculate.s0() to derive the stable age
  # distribution for initialising the population.  This ensures the population
  # starts near equilibrium rather than experiencing artificial transient
  # dynamics during the burn-in.

  n_classes <- max_age + 1L
  A <- matrix(0, nrow = n_classes, ncol = n_classes)

  # Age-specific fecundity = P(mature at age) × litter_size × P(female offspring)
  #                          × P(breeding this year) × P(fertile)
  ff <- litter_size * female_fraction * pi_2 * (1 - infertility_f)
  # Shift the ogive right by one year to match IBM timing: newly mature females
  # cannot breed in their first year of maturity (they enter at S3).
  ogive_f_leslie <- c(0, ogive_f[seq_len(max_age)])
  f_vec <- ogive_f_leslie * ff

  A[1, ] <- f_vec                                                # fecundity row
  for (i in seq_len(max_age)) A[i + 1L, i] <- survival[i]       # survival sub-diagonal
  # Note: survival[1] is already s0, so A[2,1] is set correctly by the loop

  # Extract the stable age distribution (right eigenvector of λ₁, normalised)
  w        <- Mod(eigen(A)$vectors[, 1])
  stable_A <- w / sum(w)

  # ═══════════════════════════════════════════════════════════════════════════
  # INITIALISE POPULATION
  # ═══════════════════════════════════════════════════════════════════════════
  # Create the initial population with:
  #   - Ages drawn from the stable age distribution
  #   - Unique IDs (sequential integers starting at 1)
  #   - Founder parentage markers (mother_id = 0, father_id = 0)
  #   - Individual maturity ages sampled from the sex-appropriate ogive
  #   - Permanent fertility flags (based on infertility rate)
  #   - Breeding states for mature fertile females (from Markov stationary dist)

  init_N    <- pmax(round(stable_A * pop_size), 0L)
  init_ages <- rep(0:max_age, times = init_N)
  n_init    <- sum(init_N)

  init_sex <- sample(c("F", "M"), n_init,
                     prob = c(female_fraction, 1 - female_fraction),
                     replace = TRUE)

  # Assign individual maturity ages from the sex-appropriate ogive
  init_mat_age <- integer(n_init)
  is_f <- init_sex == "F"
  init_mat_age[is_f]  <- sample_mat_ages(ogive_f, sum(is_f))
  init_mat_age[!is_f] <- sample_mat_ages(ogive_m, sum(!is_f))

  # Assign permanent fertility status
  init_fertile <- rep(TRUE, n_init)
  if (infertility_f > 0) init_fertile[is_f]  <- runif(sum(is_f))  >= infertility_f
  if (infertility_m > 0) init_fertile[!is_f] <- runif(sum(!is_f)) >= infertility_m

  # Assign Markov breeding state to mature, fertile females only.
  # Immature females, infertile females, and all males get NA (not in cycle).
  init_breed_state <- rep(NA_integer_, n_init)
  mature_f <- which(is_f & init_ages >= init_mat_age & init_fertile)
  if (length(mature_f) > 0L) {
    init_breed_state[mature_f] <- sample(
      1:3, length(mature_f), replace = TRUE,
      prob = c(pi_1, pi_2, pi_3)
    )
  }

  pop <- data.table(
    id          = seq_len(n_init),
    birth_year  = 0L,                      # founders are born in "year 0"
    age         = init_ages,
    sex         = init_sex,
    mat_age     = init_mat_age,            # age at which this individual matures
    mother_id   = 0L,                      # 0 = founder (no tracked mother)
    father_id   = 0L,                      # 0 = founder (no tracked father)
    breed_state = init_breed_state,        # 1=S1/pregnant, 2=S2/calf, 3=S3/resting, NA=not breeding
    fertile     = init_fertile,            # FALSE = permanently infertile
    population  = 1L                       # population ID (for future multi-pop support)
  )

  # ── Pod / superpod initialisation ──
  # Pods are family groups of ~pod_size individuals.
  # Superpods are communities of superpod_size pods — they serve as both
  # mating units (males mate within their superpod) and sampling units.
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

  pop_counts <- vector("list", total_years)  # one summary per year
  snapshots  <- if (do_snapshots) vector("list", length(s_years)) else list()
  snap_names <- character(0)
  snap_idx   <- 0L
  next_id    <- n_init + 1L     # running counter for unique individual IDs

  # Bull registry: tracks the dominant male (strong_bull mode) in each superpod.
  # Value of 0 means the position is vacant and needs to be filled.
  bull_registry <- NULL
  if (use_pods && !is.null(male_behavior) && male_behavior == "strong_bull") {
    bull_registry <- rep(0L, n_sp)
  }

  message(sprintf(
    "Starting simulation: %d years burn-in (2 x max_age) + %d years = %d total  (N0 = %s)",
    burn_in, num_years, total_years, format(n_init, big.mark = ",")
  ))

  # ═══════════════════════════════════════════════════════════════════════════
  # MAIN SIMULATION LOOP
  # ═══════════════════════════════════════════════════════════════════════════
  # Each iteration = one year.  The order of operations within a year is:
  #   1. Survival (stochastic; each individual survives independently)
  #   2. Aging (deterministic; survivors age by one year)
  #   3. Between-year superpod reshuffling + cow-calf following
  #   4. Markov breeding state transitions
  #   5. Create offspring (with full parentage tracking)
  #   6. Store snapshot (if this year is a snapshot year)
  #   7. Record population summary statistics

  for (yr in seq_len(total_years)) {

    # ─── 1 & 2. Survival + aging ─────────────────────────────────────────
    # Each individual survives with probability surv_vec[age + 1].
    # Survivors age by one year; those exceeding max_age are removed.
    rates <- surv_vec[pop$age + 1L]
    alive <- runif(nrow(pop)) <= rates
    pop   <- pop[alive]

    set(pop, j = "age", value = pop$age + 1L)
    pop <- pop[pop$age <= max_age]

    if (nrow(pop) == 0L) stop("Population went extinct in year ", yr, ".")

    # ─── 3. Between-year superpod reshuffling + cow-calf following ───────
    # Individuals above weaning age independently decide whether to emigrate.
    # Movers join a random pod in a DIFFERENT superpod.
    # Dependent calves (below weaning_age) follow their mother's movements.
    if (use_pods && !is.null(stickiness_year)) {
      elig <- which(pop$age >= wa)
      if (length(elig) > 0L) {
        elig_sex  <- pop$sex[elig]
        stay_prob <- ifelse(elig_sex == "F", stick_yr_F, stick_yr_M)
        movers    <- elig[runif(length(elig)) > stay_prob]

        if (length(movers) > 0L) {
          current_sp <- pop$superpod[movers]
          all_pods   <- unique(pop$pod)
          # For each superpod, find pods in OTHER superpods (emigration targets)
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

      # Cow-calf following: dependent calves are reassigned to their mother's
      # current pod/superpod.  If the mother is dead, the calf stays put.
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
    # IMPORTANT: All state groups are identified BEFORE any transitions are
    # applied.  This prevents double transitions (e.g., S2 → S3 → S1) within
    # a single year.

    # Newly mature, fertile females enter the breeding cycle at S3 (resting)
    new_mature <- which(pop$sex == "F" & pop$age == pop$mat_age &
                          is.na(pop$breed_state) & pop$fertile)
    if (length(new_mature) > 0L) {
      set(pop, i = new_mature, j = "breed_state", value = 3L)
    }

    # Snapshot current states before any transitions
    mother_rows <- which(pop$breed_state == 1L)   # S1: pregnant → will give birth
    s2_idx      <- which(pop$breed_state == 2L)   # S2: with calf
    s3_idx      <- which(pop$breed_state == 3L)   # S3: resting

    # S2 → S1 (prob ψ₂) or S2 → S3 (prob 1-ψ₂)
    if (length(s2_idx) > 0L) {
      new_state_s2 <- ifelse(runif(length(s2_idx)) < psi_2, 1L, 3L)
      set(pop, i = s2_idx, j = "breed_state", value = new_state_s2)
    }
    # S3 → S1 (prob ψ₃) or S3 → S3 (prob 1-ψ₃)
    if (length(s3_idx) > 0L) {
      new_state_s3 <- ifelse(runif(length(s3_idx)) < psi_3, 1L, 3L)
      set(pop, i = s3_idx, j = "breed_state", value = new_state_s3)
    }
    # S1 → S2: pregnant females give birth (deterministic)
    if (length(mother_rows) > 0L) {
      set(pop, i = mother_rows, j = "breed_state", value = 2L)
    }

    # ─── 5. Create offspring ─────────────────────────────────────────────
    # Only mature, fertile males are eligible to sire offspring.
    mature_male_mask <- pop$sex == "M" & pop$age >= pop$mat_age & pop$fertile
    has_males <- any(mature_male_mask)

    if (length(mother_rows) > 0L && has_males) {
      n_mothers <- length(mother_rows)

      # Each mother draws her number of mates from the num_mates vector.
      # num_mates can be a scalar (all mothers get the same) or a vector
      # (each mother randomly samples from the available values).
      n_mates_vec  <- sample(num_mates, n_mothers, replace = TRUE)
      # Litter size: minimum 1, drawn from 1 + Poisson(litter_size - 1)
      litter_sizes <- 1L + rpois(n_mothers, lambda = litter_size - 1)
      max_nm       <- max(n_mates_vec)

      # ── Build father matrix ──
      # Each row = one mother; each column = one of her potential mates.
      # Later, each offspring randomly selects one column as its actual father.
      father_mat <- matrix(NA_integer_, nrow = n_mothers, ncol = max_nm)

      if (use_pods) {
        # Superpod-based mating: males can only mate with females in their superpod
        mother_superpods <- pop$superpod[mother_rows]
        all_mature_ids   <- pop$id[mature_male_mask]
        all_mature_sps   <- pop$superpod[mature_male_mask]

        if (!is.null(male_behavior) && male_behavior == "strong_bull") {
          # ── Strong bull mode ──
          # One dominant male per superpod sires all offspring.
          # Bulls hold their position until death; vacancies are filled by
          # randomly selecting a mature male from that superpod.

          # Remove dead bulls from the registry
          alive_ids <- pop$id
          bull_registry[!bull_registry %in% c(0L, alive_ids)] <- 0L

          # Fill vacant positions by randomly electing from mature males
          # in the superpod (random selection gives realistic multi-year tenure)
          vacant <- which(bull_registry == 0L)
          if (length(vacant) > 0L) {
            mature_males_dt <- data.table(
              id       = all_mature_ids,
              age      = pop$age[mature_male_mask],
              superpod = all_mature_sps
            )
            candidates <- mature_males_dt[
              superpod %in% vacant,
              .(bull_id = sample(id, 1L)),
              by = superpod
            ]
            if (nrow(candidates) > 0L) {
              bull_registry[candidates$superpod] <- candidates$bull_id
            }
          }

          # Assign the bull as father for all mothers in each superpod
          for (sp in unique(mother_superpods)) {
            sp_mask <- which(mother_superpods == sp)
            bull_id <- bull_registry[sp]
            # Fallback: if no bull registered (e.g., no mature males in that
            # superpod), randomly pick any mature male from the population
            if (bull_id == 0L) bull_id <- sample(all_mature_ids, 1L)
            father_mat[sp_mask, ] <- bull_id
          }

        } else {
          # ── Random mating mode ──
          # Each mother's mates are sampled with replacement from mature males
          # in her superpod.  Any male can sire with multiple females.
          father_by_sp <- split(all_mature_ids, all_mature_sps)
          for (sp in unique(mother_superpods)) {
            sp_mask   <- which(mother_superpods == sp)
            n_sp_moms <- length(sp_mask)
            pool      <- father_by_sp[[as.character(sp)]]
            # Fallback: if no males in this superpod, use all males
            if (is.null(pool) || length(pool) == 0L) pool <- all_mature_ids
            father_mat[sp_mask, ] <- sample(pool, n_sp_moms * max_nm,
                                            replace = TRUE)
          }
        }

      } else {
        # ── No pod structure: global random mating ──
        # All mature males are in one global pool
        all_father_ids <- pop$id[mature_male_mask]
        father_mat[]   <- sample(all_father_ids, n_mothers * max_nm,
                                 replace = TRUE)
      }

      # ── Assign offspring to parents ──
      n_yoy <- sum(litter_sizes)

      # Map each offspring back to its mother's index in mother_rows
      off_mother_idx <- rep(seq_len(n_mothers), times = litter_sizes)
      # Each offspring randomly picks one of its mother's mates as its father
      off_n_mates    <- n_mates_vec[off_mother_idx]
      mate_col       <- as.integer(ceiling(runif(n_yoy) * off_n_mates))

      off_mother_id <- pop$id[mother_rows[off_mother_idx]]
      off_father_id <- father_mat[cbind(off_mother_idx, mate_col)]

      # ── Enforce max_females cap ──
      # If max_females is set, no single male should sire more than max_females
      # offspring per year.  Excess offspring are reassigned to other available
      # males from the same superpod (or global pool if no pods).
      if (!is.null(max_females)) {
        fid_tab  <- table(off_father_id)
        over_ids <- as.integer(names(fid_tab[fid_tab > max_females]))

        if (length(over_ids) > 0L) {
          # Build superpod-specific male pools for reassignment
          if (use_pods) {
            avail_by_sp <- split(all_mature_ids, all_mature_sps)
          }
          for (fid in over_ids) {
            idx     <- which(off_father_id == fid)
            keep    <- sample(idx, max_females)  # randomly keep max_females
            redo    <- setdiff(idx, keep)         # reassign the rest
            for (ri in redo) {
              if (use_pods) {
                # Reassign to another male from the same superpod
                mom_sp <- pop$superpod[mother_rows[off_mother_idx[ri]]]
                pool   <- avail_by_sp[[as.character(mom_sp)]]
                pool   <- pool[pool != fid]
                if (is.null(pool) || length(pool) == 0L)
                  pool <- all_mature_ids[all_mature_ids != fid]
              } else {
                pool <- all_father_ids[all_father_ids != fid]
              }
              if (length(pool) > 0L)
                off_father_id[ri] <- sample(pool, 1L)
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

      # Build the newborn data.table with full parentage
      yoy <- data.table(
        id          = seq.int(next_id, length.out = n_yoy),
        birth_year  = as.integer(yr),
        age         = 0L,
        sex         = yoy_sex,
        mat_age     = yoy_mat_age,
        mother_id   = off_mother_id,
        father_id   = off_father_id,
        breed_state = NA_integer_,       # newborns are not in the breeding cycle
        fertile     = yoy_fertile,
        population  = 1L
      )

      # Offspring inherit their mother's pod and superpod
      if (use_pods) {
        off_pods <- pop$pod[mother_rows[off_mother_idx]]
        set(yoy, j = "pod",      value = off_pods)
        set(yoy, j = "superpod", value = pod_to_sp[off_pods])
      }

      next_id <- next_id + n_yoy
      pop <- rbindlist(list(pop, yoy), use.names = TRUE)

    } # end breeding

    # ─── 6. Snapshot ─────────────────────────────────────────────────────
    # Store a deep copy of the current population for later sampling.
    # copy() is essential — without it, later in-place modifications to pop
    # would corrupt the stored snapshot.
    if (do_snapshots && yr %in% s_years) {
      snap_idx <- snap_idx + 1L
      snapshots[[snap_idx]] <- copy(pop)
      snap_names <- c(snap_names, as.character(yr))
    }

    # ─── 7. Record population metrics ────────────────────────────────────
    # Aggregate counts by sex and age for each year (used for pop_summary)
    yr_counts <- pop[, .N, by = .(sex, age)]
    set(yr_counts, j = "year", value = as.integer(yr))
    pop_counts[[yr]] <- yr_counts

    # Print progress at key milestones
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

  # Combine yearly summary tables into one data.table
  pop_summary <- rbindlist(pop_counts)
  if (do_snapshots) names(snapshots) <- snap_names

  message(sprintf(
    "Simulation complete: %d total years (%d burn-in + %d post), final N = %s, %d snapshots stored.",
    total_years, burn_in, num_years,
    format(nrow(pop), big.mark = ","), length(snapshots)
  ))

  # Return simulation output.  This list is designed to be passed to
  # sample.pop().  sim_config is passed through so sampling can access
  # parameters like weaning_age.
  invisible(list(
    pop_summary = pop_summary,    # data.table: year × sex × age counts
    snapshots   = snapshots,      # named list of full-population data.tables
    pod_to_sp   = pod_to_sp,      # integer vector: pod index → superpod index
    sim_config  = sim_config      # original config (for sample.pop)
  ))
}
