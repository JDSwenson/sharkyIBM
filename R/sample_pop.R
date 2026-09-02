#' Sample individuals from population snapshots using a hierarchical scheme
#'
#' Takes population snapshots produced by \code{simulate.pop()} and applies a
#' sampling design based on purse seine operations: multiple trips per year,
#' multiple sets per trip, with social reshuffling between sets and trips.
#' This function operates on stored snapshots and is cheap to run — it can be
#' called repeatedly with different sampling parameters without re-running the full population simulation.
#'
#' @param sim_output List returned by \code{simulate.pop()}.  Must contain
#'   \code{snapshots}, \code{pod_to_sp}, and \code{sim_config}.
#' @param n_trips Integer. Number of sampling trips per year.
#' @param n_sets Integer. Number of purse seine sets per trip.
#' @param sample_size Integer. Number of individuals sampled.  Interpretation
#'   depends on \code{sample_per}.
#' @param sample_per Character. \code{"set"} (sample \code{sample_size}
#'   individuals per set) or \code{"trip"} (sample \code{sample_size}
#'   individuals total per trip, distributed across sets).
#' @param sampling Character. Determines how sampling occurs during each "trip". Options are \code{"random"} (uniform random from the full
#'   population) or \code{"superpod"} (sample from a single randomly chosen
#'   superpod per trip). The latter is meant to emulate a purse seine vessel that encounters a single superpod and follows that same superpod throughout the trip.
#' @param stickiness_set Numeric (scalar or length-2 vector), bounded (0,1). Probability of staying in the same superpod between consecutive sets within a trip.  Scalar = same for both sexes; \code{c(female, male)} for sex-specific rates. Default 1 (no reshuffling between sets).
#' @param stickiness_trip Numeric (scalar or length-2 vector), bounded (0,1). Probability of staying in the same superpod between consecutive trips within a year. Default 1 (no reshuffling between trips).
#' @param superpod_pool NULL or a list of integer vectors.  Controls which
#'   superpods each trip can encounter.
#'   \itemize{
#'     \item NULL: all superpods are available for all trips.
#'     \item List of length \code{n_trips}: each element is a vector of
#'           superpod IDs available to that trip (e.g., for spatial structure
#'           where northern trips encounter different communities than
#'           southern trips).
#'   }
#'
#' @details
#'
#' ## Hierarchical Sampling Structure
#'
#' The function implements a three-level sampling hierarchy mimicking purse seine
#' operations:
#' \itemize{
#'   \item **Year**: All snapshots in the requested year(s) are processed.
#'   \item **Trip** (n_trips per year): Each trip represents a deployment of a
#'     sampling vessel.
#'   \item **Set** (n_sets per trip): Each set represents a single purse seine
#'     operation during a trip.
#' }
#'
#' Within each trip, a single superpod is encountered (or a random subset if
#' \code{stickiness_trip} < 1.0 reshuffles between trips). Within each set, the
#' same superpod is encountered; however, \code{stickiness_set} controls whether
#' individuals remain in that superpod between sets.
#'
#' ## Sampling Modes
#'
#' **sampling = "random"**: At each set, a sample is drawn uniformly at random
#' from the entire living population (across all superpods). This is fast and
#' useful for baseline comparisons but unrealistic for marine mammal sampling.
#'
#' **sampling = "superpod"**: At the start of each trip, a superpod is chosen
#' (uniformly at random from \code{superpod_pool} if provided, or all available
#' superpods otherwise). Sampling within that trip occurs only from that superpod
#' and its nearby communities. At the start of each set within the trip,
#' \code{stickiness_set} determines whether the catch remains in the same superpod
#' or moves to a different one. This realistically models social structure and
#' catch patterns.
#'
#' ## Stickiness Parameters
#'
#' **stickiness_trip**: Between-trip fidelity within a year (0--1). Each trip
#' starts with a randomly selected superpod from \code{superpod_pool}.  Between
#' trips, each individual in the current superpod remains with probability
#' \code{stickiness_trip}; those that move (probability 1 - stickiness_trip)
#' emigrate to a different superpod chosen uniformly at random. Can be scalar
#' (both sexes) or sex-specific \code{c(female, male)}.
#'
#' **stickiness_set**: Between-set fidelity within a trip (0--1). Within each
#' trip, the same superpod is targeted for all sets. Before each set, individuals
#' remain in their current superpod with probability \code{stickiness_set}; those
#' that move emigrate to a different superpod.  Can be scalar or sex-specific.
#'
#' These parameters control **duplication rates**: when \code{stickiness_set} is
#' high (close to 1.0), the same individuals are likely to be caught in multiple
#' sets within a trip, increasing duplication.  Conversely, low stickiness means
#' more turnover and more unique individuals per trip.
#'
#' ## Sampling Modes: sample_per
#'
#' **sample_per = "set"**: Each set independently samples \code{sample_size}
#' individuals (without replacement within a set). Total per trip = sample_size *
#' n_sets. Allows different catch sizes per set.
#'
#' **sample_per = "trip"**: The total per trip is \code{sample_size}, distributed
#' evenly across n_sets (sample_size / n_sets per set). If sample_size is not
#' evenly divisible, the first few sets get an extra individual to reach
#' sample_size total. Enforces a per-trip quota.
#'
#' ## Spatial Structure: superpod_pool
#'
#' Controls which superpods are accessible to each trip:
#' \itemize{
#'   \item NULL (default): All superpods are equally available to all trips.
#'   \item List of length n_trips: \code{superpod_pool[[i]]} is a vector of
#'     superpod IDs available to trip i. Useful for enforcing spatial structure
#'     (e.g., northern trips encounter communities 1--30, southern trips encounter
#'     communities 31--60).
#' }
#'
#' ## Cow-Calf Dynamics
#'
#' Dependent calves (age < \code{weaning_age}) follow their mother's pod and
#' superpod and are sampled accordingly. Once weaned, calves are treated as
#' independent individuals and subject to the same reshuffling as adults.
#'
#' ## Output Structure
#'
#' Returns a \code{data.table} with one row per catch event (individual × set).
#' Columns include:
#' \itemize{
#'   \item \strong{Individual data}: \code{id}, \code{birth_year}, \code{age},
#'     \code{sex}, \code{mother_id}, \code{father_id}.
#'   \item \strong{Population}: \code{population} (always 1 for single-population
#'     simulations).
#'   \item \strong{Social structure}: \code{pod}, \code{superpod} (if used).
#'   \item \strong{Sampling metadata}: \code{year}, \code{trip}, \code{set}.
#' }
#'
#' **Duplicate detection**: The same individual may appear multiple times (one row
#' per set where caught). To count unique samples per trip (or year), deduplicate
#' by \code{id}:
#' \code{unique(samples, by = "id")} or \code{samples[!duplicated(id)]}.
#'
#' **Exclusions**: Internal columns (\code{s2_year}, \code{breed_state},
#' \code{fertile}, \code{mat_age}) are dropped before return to keep the output
#' focused on sampling-relevant fields.
#'
#' @return A \code{data.table} of sampled individuals with columns:
#'   \code{id}, \code{birth_year}, \code{age}, \code{sex}, \code{mother_id},
#'   \code{father_id}, \code{population}, \code{year}, \code{trip}, \code{set},
#'   and (if pods are used) \code{pod}, \code{superpod}.
#'
#'   The same individual may appear multiple times if caught in different sets.
#'   Deduplicate by \code{id} (within or across years) to count unique samples.
#'
#' @references
#' Barlow, J., & Boveng, P. (1991). Modeling age-specific mortality for marine
#' mammal populations. Marine Mammal Science, 7(1), 50-65.
#'
#' Caswell, H. (2001). Matrix Population Models: Construction, Analysis, and
#' Interpretation (2nd ed.). Sinauer Associates.
#'
#' @importFrom data.table data.table set rbindlist copy uniqueN
#' @importFrom stats runif
#' @export
#' @rawNamespace export(sample.pop)
sample.pop <- function(sim_output,
                       n_trips,
                       n_sets,
                       sample_size,
                       sample_per      = "set",
                       sampling        = "superpod",
                       stickiness_set  = 1,
                       stickiness_trip = 1,
                       superpod_pool   = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # INPUT VALIDATION
  # ═══════════════════════════════════════════════════════════════════════════

  if (!is.list(sim_output))
    stop("`sim_output` must be a list (the output of simulate.pop()).")

  required_fields <- c("snapshots", "sim_config")
  missing_fields <- setdiff(required_fields, names(sim_output))
  if (length(missing_fields) > 0L)
    stop("`sim_output` is missing required fields: ",
         paste(missing_fields, collapse = ", "),
         ". Did you pass the output of simulate.pop()?")

  if (!is.numeric(n_trips) || length(n_trips) != 1L || n_trips < 1 ||
      n_trips != round(n_trips))
    stop("`n_trips` must be a positive integer.")

  if (!is.numeric(n_sets) || length(n_sets) != 1L || n_sets < 1 ||
      n_sets != round(n_sets))
    stop("`n_sets` must be a positive integer.")

  if (!is.numeric(sample_size) || length(sample_size) != 1L || sample_size < 1 ||
      sample_size != round(sample_size))
    stop("`sample_size` must be a positive integer.")

  if (!sample_per %in% c("set", "trip"))
    stop('`sample_per` must be "set" or "trip".')

  if (!sampling %in% c("random", "superpod"))
    stop('`sampling` must be "random" or "superpod".')

  # Stickiness values must be probabilities
  if (!is.numeric(stickiness_set) || !length(stickiness_set) %in% c(1L, 2L))
    stop("`stickiness_set` must be a numeric scalar or length-2 vector.")
  if (any(stickiness_set < 0 | stickiness_set > 1))
    stop("`stickiness_set` values must be between 0 and 1.")

  if (!is.numeric(stickiness_trip) || !length(stickiness_trip) %in% c(1L, 2L))
    stop("`stickiness_trip` must be a numeric scalar or length-2 vector.")
  if (any(stickiness_trip < 0 | stickiness_trip > 1))
    stop("`stickiness_trip` values must be between 0 and 1.")

  # ═══════════════════════════════════════════════════════════════════════════
  # EXTRACT SIMULATION COMPONENTS
  # ═══════════════════════════════════════════════════════════════════════════

  snapshots   <- sim_output$snapshots
  pod_to_sp   <- sim_output$pod_to_sp
  sim_config  <- sim_output$sim_config
  weaning_age <- sim_config$weaning_age

  use_pods    <- !is.null(pod_to_sp)
  use_weaning <- !is.null(weaning_age)

  # Must have at least one snapshot to sample from
  if (length(snapshots) == 0L) {
    stop("No snapshots available. Run simulate.pop() with sample_years specified.")
  }

  # ── Parse sex-specific stickiness ──
  # Between-set stickiness controls reshuffling between consecutive purse
  # seine sets within the same trip.  Higher values = more duplication.
  if (length(stickiness_set) == 1L) {
    stick_set_F <- stick_set_M <- stickiness_set
  } else {
    stick_set_F <- stickiness_set[1]
    stick_set_M <- stickiness_set[2]
  }

  # Between-trip stickiness controls reshuffling between trips within a year.
  if (length(stickiness_trip) == 1L) {
    stick_trip_F <- stick_trip_M <- stickiness_trip
  } else {
    stick_trip_F <- stickiness_trip[1]
    stick_trip_M <- stickiness_trip[2]
  }

  # ── Validate superpod_pool ──
  # If provided, must be a list with one element per trip, each element being
  # a vector of superpod IDs that trip can encounter.
  if (!is.null(superpod_pool)) {
    if (!is.list(superpod_pool))
      stop("`superpod_pool` must be NULL or a list.")
    if (length(superpod_pool) != n_trips)
      stop("`superpod_pool` must be NULL or a list of length n_trips (",
           n_trips, "). Got length ", length(superpod_pool), ".")
  }

  # ── Per-set sample size ──
  # In "trip" mode, distribute sample_size across sets (ceiling to ensure
  # we get at least sample_size total, then trim excess at end of trip).
  if (sample_per == "trip") {
    per_set <- ceiling(sample_size / n_sets)
  } else {
    per_set <- sample_size
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # INTERNAL HELPER: RESHUFFLE SUPERPOD MEMBERSHIP
  # ═══════════════════════════════════════════════════════════════════════════
  # Reshuffles the working population between sets or trips.  Each eligible
  # individual (age >= weaning_age) independently decides whether to stay in
  # its current superpod or move to a random pod in a different superpod.
  # Dependent calves follow their mothers after the reshuffle.
  #
  # This is the same mechanism used in simulate.pop() for between-year
  # reshuffling, but applied here with set-level or trip-level stickiness.

  reshuffle <- function(pop, stick_F, stick_M) {
    if (!use_pods) return(pop)

    wa <- if (use_weaning) weaning_age else 0L
    elig <- which(pop$age >= wa)
    if (length(elig) == 0L) return(pop)

    # Each eligible individual stays with probability stick_F/M
    elig_sex  <- pop$sex[elig]
    stay_prob <- ifelse(elig_sex == "F", stick_F, stick_M)
    movers    <- elig[runif(length(elig)) > stay_prob]

    if (length(movers) > 0L) {
      current_sp <- pop$superpod[movers]
      all_pods   <- unique(pop$pod)
      # Build lookup: for each superpod, which pods are in OTHER superpods?
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
        # reinterprets a length-1 numeric x as the range 1:x). Indexing by
        # position avoids this.
        new_pods[mask] <- pool[sample.int(length(pool), length(mask), replace = TRUE)]
      }
      set(pop, i = movers, j = "pod",      value = new_pods)
      set(pop, i = movers, j = "superpod", value = pod_to_sp[new_pods])
    }

    # Cow-calf following: calves below weaning_age are reassigned to their
    # mother's current pod/superpod (if mother is alive)
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

    pop
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # MAIN SAMPLING LOOP
  # ═══════════════════════════════════════════════════════════════════════════
  # Nested loop: snapshots (years) → trips → sets
  #
  # Key design decisions:
  #   - Within a trip, the vessel encounters the SAME superpod across all sets.
  #     The superpod is chosen once per trip, not per set.
  #   - Between sets, individuals reshuffle with probability (1 - stickiness_set).
  #     High stickiness = individuals stay in the same superpod between sets,
  #     increasing the chance of recapturing the same individuals.
  #   - Between trips, a separate reshuffle occurs with stickiness_trip.
  #   - The working_pop is a mutable copy — reshuffling modifies pod/superpod
  #     assignments in place without affecting the stored snapshot.

  all_samples <- vector("list", length(snapshots) * n_trips * n_sets)
  sample_counter <- 0L

  for (snap_i in seq_along(snapshots)) {
    snap_year   <- as.integer(names(snapshots)[snap_i])
    # Deep copy so reshuffling doesn't modify the stored snapshot
    working_pop <- copy(snapshots[[snap_i]])

    for (trip in seq_len(n_trips)) {

      # Determine which superpods are reachable for this trip.
      # If superpod_pool is set, each trip has a restricted pool (spatial structure).
      if (!is.null(superpod_pool)) {
        avail_sps <- superpod_pool[[trip]]
      } else {
        avail_sps <- if (use_pods) unique(working_pop$superpod) else NULL
      }

      # Within a trip, the vessel encounters ONE superpod across all sets.
      # This represents operating in the same geographic area during a trip.
      chosen_sp <- NULL
      if (sampling == "superpod" && use_pods) {
        # Indexing by position avoids sample()'s length-1 numeric footgun.
        chosen_sp <- avail_sps[sample.int(length(avail_sps), 1L)]
      }

      # Track total samples drawn this trip (for "trip" mode bookkeeping)
      trip_sample_count <- 0L

      for (s in seq_len(n_sets)) {

        # ── Determine draw size for this set ──
        if (sample_per == "trip") {
          # In "trip" mode, distribute sample_size across sets, capping at
          # the remaining quota for this trip
          remaining   <- sample_size - trip_sample_count
          this_draw   <- min(per_set, remaining)
          if (this_draw <= 0L) next
        } else {
          this_draw <- per_set
        }

        # ── Draw individuals ──
        if (!is.null(chosen_sp)) {
          # Superpod sampling: only draw from the trip's chosen superpod
          sp_rows <- which(working_pop$superpod == chosen_sp)
          n_avail <- length(sp_rows)
          n_draw  <- min(this_draw, n_avail)
          if (n_draw > 0L) {
            draw_rows <- sp_rows[sample.int(n_avail, n_draw)]
          } else {
            next  # skip if superpod is empty
          }

        } else {
          # Random sampling: draw from the entire population
          n_draw    <- min(this_draw, nrow(working_pop))
          draw_rows <- sample.int(nrow(working_pop), n_draw)
        }

        # Extract the sampled individuals and tag with metadata
        sampled <- working_pop[draw_rows]
        set(sampled, j = "year", value = snap_year)
        set(sampled, j = "trip", value = as.integer(trip))
        set(sampled, j = "set",  value = as.integer(s))

        sample_counter <- sample_counter + 1L
        all_samples[[sample_counter]] <- sampled

        trip_sample_count <- trip_sample_count + n_draw

        # ── Between-set reshuffling (except after the last set) ──
        # Individuals may leave the superpod between consecutive sets.
        # This controls within-trip duplication: high stickiness_set means
        # individuals tend to stay, increasing recapture probability.
        if (s < n_sets) {
          working_pop <- reshuffle(working_pop, stick_set_F, stick_set_M)
        }

      } # end sets

      # ── Between-trip reshuffling (except after the last trip) ──
      # Larger-scale reshuffling between separate trips to the study area.
      if (trip < n_trips) {
        working_pop <- reshuffle(working_pop, stick_trip_F, stick_trip_M)
      }

    } # end trips
  } # end snapshots

  # ═══════════════════════════════════════════════════════════════════════════
  # COMPILE AND RETURN
  # ═══════════════════════════════════════════════════════════════════════════

  # Trim unused pre-allocated list slots and combine all samples
  all_samples <- all_samples[seq_len(sample_counter)]
  samples_df  <- rbindlist(all_samples, use.names = TRUE)

  # Select only the columns relevant for downstream CKMR analysis.
  # Internal columns like breed_state, mat_age, fertile, and calf_id are dropped.
  keep_cols <- c("id", "birth_year", "age", "sex", "mother_id", "father_id",
                 "population", "year", "trip", "set")
  if (use_pods) keep_cols <- c(keep_cols, "pod", "superpod")
  samples_df <- samples_df[, .SD, .SDcols = intersect(keep_cols, names(samples_df))]

  n_unique <- uniqueN(samples_df$id)
  message(sprintf(
    "Sampling complete: %d total captures, %d unique individuals across %d year(s).",
    nrow(samples_df), n_unique, length(snapshots)
  ))

  samples_df
}
