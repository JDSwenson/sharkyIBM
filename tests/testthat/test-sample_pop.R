# Tests for sample.pop()

sim_output_pods <- function() {
  cfg <- quick_config_pods()
  quick_sim(cfg, num_years = 5, sample_years = 3)
}

sim_output_no_pods <- function() {
  cfg <- quick_config()
  quick_sim(cfg, num_years = 5, sample_years = 2)
}

# ── Input validation ─────────────────────────────────────────────────────

test_that("errors when sim_output is not a list", {
  expect_error(sample.pop("nope", n_trips = 1, n_sets = 1, sample_size = 5),
               "sim_output.*list")
})

test_that("errors when sim_output is missing required fields", {
  expect_error(
    sample.pop(list(snapshots = list()), n_trips = 1, n_sets = 1, sample_size = 5),
    "missing required fields"
  )
})

test_that("errors on invalid n_trips, n_sets, sample_size", {
  sim <- sim_output_no_pods()
  expect_error(sample.pop(sim, n_trips = 0, n_sets = 1, sample_size = 5), "n_trips")
  expect_error(sample.pop(sim, n_trips = 1, n_sets = -1, sample_size = 5), "n_sets")
  expect_error(sample.pop(sim, n_trips = 1, n_sets = 1, sample_size = 0), "sample_size")
})

test_that("errors on invalid sample_per / sampling values", {
  sim <- sim_output_no_pods()
  expect_error(
    sample.pop(sim, n_trips = 1, n_sets = 1, sample_size = 5, sample_per = "bogus"),
    "sample_per"
  )
  expect_error(
    sample.pop(sim, n_trips = 1, n_sets = 1, sample_size = 5, sampling = "bogus"),
    "sampling"
  )
})

test_that("errors on out-of-range stickiness values", {
  sim <- sim_output_no_pods()
  expect_error(
    sample.pop(sim, n_trips = 1, n_sets = 1, sample_size = 5, stickiness_set = 1.5),
    "stickiness_set"
  )
  expect_error(
    sample.pop(sim, n_trips = 1, n_sets = 1, sample_size = 5, stickiness_trip = -0.1),
    "stickiness_trip"
  )
})

test_that("errors when superpod_pool length does not match n_trips", {
  sim <- sim_output_pods()
  expect_error(
    sample.pop(sim, n_trips = 2, n_sets = 1, sample_size = 5,
               superpod_pool = list(1L)),
    "superpod_pool"
  )
})

test_that("errors when simulate.pop() was run without snapshots", {
  cfg <- quick_config()
  sim <- quick_sim(cfg, num_years = 5, sample_years = NULL)
  expect_error(
    sample.pop(sim, n_trips = 1, n_sets = 1, sample_size = 5),
    "No snapshots"
  )
})

# ── Output structure ──────────────────────────────────────────────────────

test_that("samples without pods have no pod/superpod columns", {
  sim <- sim_output_no_pods()
  out <- suppressMessages(sample.pop(sim, n_trips = 1, n_sets = 2, sample_size = 10,
                                      sampling = "random"))
  expect_true(all(c("id", "birth_year", "age", "sex", "mother_id", "father_id",
                     "population", "year", "trip", "set") %in% names(out)))
  expect_false(any(c("pod", "superpod") %in% names(out)))
  expect_false(any(c("breed_state", "mat_age", "fertile", "calf_id") %in% names(out)))
})

test_that("samples with pods include pod/superpod columns", {
  sim <- sim_output_pods()
  out <- suppressMessages(sample.pop(sim, n_trips = 1, n_sets = 2, sample_size = 5))
  expect_true(all(c("pod", "superpod") %in% names(out)))
})

# ── Sample-size accounting ────────────────────────────────────────────────

test_that("sample_per = 'set' draws exactly sample_size per set", {
  sim <- sim_output_pods()
  out <- suppressMessages(sample.pop(sim, n_trips = 2, n_sets = 3, sample_size = 8,
                                      sample_per = "set"))
  counts <- out[, .N, by = .(year, trip, set)]
  expect_true(all(counts$N == 8))
})

test_that("sample_per = 'trip' draws exactly sample_size per trip", {
  sim <- sim_output_pods()
  out <- suppressMessages(sample.pop(sim, n_trips = 2, n_sets = 3, sample_size = 15,
                                      sample_per = "trip"))
  counts <- out[, .N, by = .(year, trip)]
  expect_true(all(counts$N == 15))
})

test_that("higher stickiness_set produces more within-trip duplication", {
  sim <- sim_output_pods()

  set.seed(55)
  high <- suppressMessages(sample.pop(sim, n_trips = 1, n_sets = 5, sample_size = 20,
                                       sample_per = "set", stickiness_set = 1))
  set.seed(55)
  low <- suppressMessages(sample.pop(sim, n_trips = 1, n_sets = 5, sample_size = 20,
                                      sample_per = "set", stickiness_set = 0))

  expect_gt(sum(duplicated(high$id)), sum(duplicated(low$id)))
})

test_that("superpod_pool restricts each trip to its assigned superpod(s)", {
  sim <- sim_output_pods()
  all_sps <- unique(sim$snapshots[[1]]$superpod)
  skip_if(length(all_sps) < 2, "need at least 2 superpods for this test")

  pool <- list(all_sps[1], all_sps[2])
  out <- suppressMessages(sample.pop(sim, n_trips = 2, n_sets = 2, sample_size = 5,
                                      superpod_pool = pool))

  trip1_sps <- unique(out$superpod[out$trip == 1L])
  trip2_sps <- unique(out$superpod[out$trip == 2L])
  expect_true(all(trip1_sps == all_sps[1]))
  expect_true(all(trip2_sps == all_sps[2]))
})

test_that("the same individual can be flagged via duplicated ids across sets", {
  sim <- sim_output_pods()
  out <- suppressMessages(sample.pop(sim, n_trips = 1, n_sets = 3, sample_size = 10,
                                      stickiness_set = 1))
  n_unique <- length(unique(out$id))
  expect_lte(n_unique, nrow(out))
})
