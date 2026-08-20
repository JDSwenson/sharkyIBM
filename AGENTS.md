# sharkyIBM Project — Conversation Summary

**Last updated:** Thursday, August 20, 2026
**Project:** Individual-Based Population Simulation for Close-Kin Mark–Recapture (CKMR) Analysis
**Species:** Pantropical spotted dolphins, eastern spinner dolphins

---

## Overview

This is a dolphin population simulation package (`sharkyIBM` — placeholder name) designed to:
- Simulate age-structured populations forward in time with full parentage tracking
- Test sampling schemes (modeled on purse seine operations) and life-history attributes with CKMR models
- Ground-truth CKMR estimates with known population demographics
- Estimate the number of usable (unique) samples obtainable per sampling trip given dolphin social cohesion

### Three-function pipeline

1. **`R/calculate_s0.R`** — Calibrates age-0 survival (s0) via Leslie-matrix + IBM bisection. Returns a bundled config object with all shared parameters.
2. **`R/simulate_population.R`** — Runs the forward-time IBM with automatic `max_age`-year burn-in, Markov breeding cycle, between-year superpod reshuffling, and parentage tracking. Returns population snapshots at user-specified years.
3. **`R/sample_pop.R`** — Takes snapshots from `simulate.pop()` and applies a hierarchical sampling scheme: trips → sets → samples, with stickiness-based reshuffling between sets and trips. Can be re-run cheaply with different sampling parameters without re-simulating.

### Parked/Legacy (not wired in)

- `R/helper_functions.R` — Older dplyr-based breeding functions with length-based growth (von Bertalanffy, MVN-correlated L_inf/K). Not called by current pipeline.
- `R/create_input_data.R` — Parked. Uses Euler-Lotka for s0, superseded by `calculate_s0.R`.

---

## Key Decisions & Design Principles

### 1. Performance First
- All code uses `data.table` for in-place column updates and fast subsetting
- Survival lookups use direct vector indexing (`surv_vec[age + 1]`)
- No father-identity tracking in `calculate.s0()` — only offspring counts matter for growth rate
- Simulation runs once; sampling experiments are cheap (operate on stored snapshots)

### 2. Leslie Matrix as Foundation
- Both `calculate.s0()` and `simulate.pop()` build a Leslie matrix to extract the stable age distribution for initialisation
- `calculate.s0()` also uses the Leslie matrix for an analytical s0 starting estimate
- Fecundity in the Leslie matrix uses the Markov stationary proportion π₂ (fraction of mature females breeding per year)

### 3. Empirical s0 Calibration
- `calculate.s0()` runs one continuous simulation with ~5-year assessment windows
- Bisects s0 when growth deviates from zero; converges after 3 consecutive stable windows (|r| < 0.005)

### 4. Bundled Config Object
- `calculate.s0()` returns a list with all shared parameters + calibration results
- `simulate.pop()` accepts this list as `sim_config`, extracting parameters internally
- Users specify life-history params once; simulation-length and sampling params are separate

### 5. Automatic Burn-in
- `simulate.pop()` runs `max_age` years of burn-in before `num_years` of post-burn-in
- Total simulation = `max_age + num_years`
- Guarantees all founders (mother_id = 0, father_id = 0) are dead before snapshots
- `pop_summary` covers ALL years (burn-in + post)

### 6. Markov Breeding Cycle
- Replaces the former deterministic `mating_periodicity`/`repro_cycle` system
- Each mature female carries a `breed_state`: S1 (pregnant), S2 (with calf), S3 (resting)
- **Transition matrix Ψ:**
  - S1 → S2 (prob 1): pregnant female gives birth
  - S2 → S1 (prob ψ₂): conceive while nursing
  - S2 → S3 (prob 1−ψ₂): wean calf, rest
  - S3 → S1 (prob ψ₃): conceive from resting
  - S3 → S3 (prob 1−ψ₃): remain resting
- `psi_2` and `psi_3` are user-defined; calving interval emerges from the Markov chain
- Newly mature females enter at S3 (resting)
- **Important implementation detail:** All state indices (S1, S2, S3) must be identified BEFORE applying transitions to prevent double transitions (S2→S3→S1) in a single year
- Stationary distribution: π₁ = π₂ = ψ₃/(2ψ₃ + 1 − ψ₂), π₃ = (1−ψ₂)·π₁/ψ₃

### 7. Pod / Superpod Social Structure
- **Pods** = family groups of `pod_size` individuals (initialisation unit)
- **Superpods** = communities of `superpod_size` pods (mating + sampling units)
- Pod-to-superpod mapping is fixed at initialisation
- Offspring inherit mother's pod and superpod

### 8. Three Levels of Stickiness
- **`stickiness_year`** (in `simulate.pop`): between-year superpod fidelity. Affects mating pools. Sex-specific (scalar or `c(female, male)`).
- **`stickiness_trip`** (in `sample.pop`): between-trip reshuffling within a year. Sex-specific.
- **`stickiness_set`** (in `sample.pop`): between-set reshuffling within a trip. Controls duplication rate. Sex-specific.
- Movers go to a random pod in a DIFFERENT superpod (emigration, not within-community reshuffling)
- Within a trip, the vessel encounters the SAME superpod across all sets. Stickiness controls whether individuals leave that community between sets.

### 9. Cow-Calf Dynamics
- **`weaning_age`** (integer or NULL): age at which a calf becomes independent
- Below `weaning_age`: calves excluded from reshuffling, follow mother's pod/superpod
- Also serves as the threshold for stickiness-based reshuffling eligibility
- Orphaned calves (mother dead) stay in their current pod
- Active in all three functions (calculate_s0, simulate_pop, sample_pop)

### 10. Mating Modes
- **`"random"`**: any mature male in the superpod may sire offspring
- **`"strong_bull"`**: persistent bull registry — oldest mature male per superpod sires all offspring; new bull elected when incumbent dies
- Cross-superpod fallback: if a superpod has no mature males, females mate with males from other superpods

---

## Function Specifications

### `calculate.s0()`

```r
calculate.s0(
  max_age, survival, pop_size, maturity_age, litter_size,
  psi_2, psi_3,
  num_mates = 1L, female_fraction = 0.5, infertility = 0,
  pod_size = NULL, superpod_size = NULL, stickiness_year = NULL,
  male_behavior = NULL, weaning_age = NULL,
  check_interval = 5L, growth_tol = 0.005, stable_required = 3L,
  max_windows = 100L
)
```

**Output:** List with 19 elements — calibration results (`s0`, `survival`, `s0_leslie`, `final_N`, `years_simulated`) plus all shared parameters.

### `simulate.pop()`

```r
simulate.pop(sim_config, num_years, sample_years = NULL)
```

**Output:** List with:
- `pop_summary` — data.table: `year, sex, age, N` for all years
- `snapshots` — named list of population data.tables at snapshot years
- `pod_to_sp` — integer vector mapping pod → superpod
- `sim_config` — passed through for `sample.pop()`

### `sample.pop()`

```r
sample.pop(
  sim_output, n_trips, n_sets, sample_size,
  sample_per = "set", sampling = "superpod",
  stickiness_set = 1, stickiness_trip = 1,
  superpod_pool = NULL
)
```

- `sample_per`: `"set"` (sample_size per set) or `"trip"` (sample_size per trip, distributed across sets)
- `superpod_pool`: NULL (all available) or list of int vectors per trip (spatial structure)
- Within a trip, the same superpod is encountered across all sets

**Output:** data.table with `id, birth_year, age, sex, mother_id, father_id, population, pod, superpod, year, trip, set`. May contain duplicates (same individual in multiple sets).

---

## Example Usage

```r
library(data.table)
source("R/calculate_s0.R")
source("R/simulate_population.R")
source("R/sample_pop.R")

result <- calculate.s0(
  max_age = 30, survival = c(0.65, rep(0.96, 30)),
  pop_size = 50000, maturity_age = 8, litter_size = 1,
  psi_2 = 0.1, psi_3 = 0.7,
  num_mates = 1L,
  pod_size = 20, superpod_size = 10,
  stickiness_year = c(0.95, 0.85),
  male_behavior = "random", weaning_age = 3L
)

sim <- simulate.pop(result, num_years = 50, sample_years = 5)

samples <- sample.pop(
  sim, n_trips = 3, n_sets = 5, sample_size = 40,
  sample_per = "set", sampling = "superpod",
  stickiness_set = 0.8, stickiness_trip = 0.5,
  superpod_pool = list(1:50, 51:100, 101:125)
)

# Usable (unique) samples
unique_samples <- samples[!duplicated(id)]
```

---

## Test Results (August 20, 2026)

✓ **Bundled config:** 19-element list with all shared params
✓ **Burn-in:** 30-year burn-in + 50-year post = 80 total; zero founders in snapshots
✓ **Snapshots:** Stored at correct years (76–80); full population with all columns
✓ **Pop summary:** Covers all 80 years
✓ **Population stability:** ~50K population stable across 80 years (growth ≈ −0.0007)
✓ **Markov breeding states:** Observed S1=32.5%, S2=30.0%, S3=37.5% (expected: 30.4/30.4/39.1)
✓ **Cow-calf dynamics:** 8492/8492 calves match mother's superpod (0 mismatches)
✓ **Stickiness effect:** 24.3% duplication at stickiness_set=1.0 vs 3.4% at stickiness_set=0.0
✓ **Per-trip sampling:** Exactly 100 samples per trip distributed across sets
✓ **Superpod pools:** Trip 1 sampled only from superpods 1–10; Trip 2 only from 11–20
✓ **Sex-specific stickiness:** Accepts `c(female, male)` vector
✓ **Strong bull mating:** 6,116 unique fathers vs 25,582 in random mode (76% fewer)
✓ **Markov transition bug caught and fixed:** indices must be computed before applying state changes

---

## Next Steps

1. **Scaling test** — run at 500K–1M individuals to benchmark performance and memory
2. **CKMR kin-pair identification** — build a function to find parent-offspring, full-sibling, and half-sibling pairs in samples and compare to ground-truth pedigree
3. **Sensitivity analysis** — vary stickiness_set across a range and plot unique-sample yield per trip
4. **End-to-end pipeline** — connect `calculate.s0()` → `simulate.pop()` → `sample.pop()` → CKMR validation

### Deferred Features
- Length-based maturity/survival (partial implementation in legacy `helper_functions.R`)
- Multi-population structure with movement
- Infertility dynamics
- Rename package from "sharkyIBM"

---

## File Structure

```
sharkyIBM/
├── R/
│   ├── calculate_s0.R          [active — s0 calibration + bundled config]
│   ├── simulate_population.R   [active — IBM with burn-in + Markov breeding + snapshots]
│   ├── sample_pop.R            [active — hierarchical sampling (trips/sets/stickiness)]
│   ├── helper_functions.R      [legacy — not wired in; has length-based growth]
│   ├── create_input_data.R     [parked — superseded]
├── AGENTS.md                   [this file]
```
