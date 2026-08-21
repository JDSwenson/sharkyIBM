# sharkyIBM Project — Conversation Summary

**Last updated:** Friday, August 21, 2026
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
2. **`R/simulate_population.R`** — Runs the forward-time IBM with automatic `2 × max_age`-year burn-in, Markov breeding cycle, between-year superpod reshuffling, and parentage tracking. Returns population snapshots at user-specified years.
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
- `calculate.s0()` runs one continuous simulation with assessment windows of `check_interval` years
- Bisects s0 when growth deviates from zero; converges after `stable_required` consecutive stable windows (|r| < `growth_tol`)
- Default tolerances: check_interval=10, growth_tol=0.001, stable_required=5 (typically ~100–220 years to converge)
- **Fecundity correction (Aug 21):** The Leslie matrix shifts the maturity ogive right by one year (`ogive_f_leslie = c(0, ogive_f[1:max_age])`) so that newly mature females contribute zero fecundity at their first mature age — matching the IBM, where they enter at S3 and cannot breed until the following year. This gives a more accurate analytical starting estimate for the bisection.

### 4. Bundled Config Object
- `calculate.s0()` returns a list with all shared parameters + calibration results
- `simulate.pop()` accepts this list as `sim_config`, extracting parameters internally
- Users specify life-history params once; simulation-length and sampling params are separate

### 5. Automatic Burn-in
- `simulate.pop()` runs `2 × max_age` years of burn-in before `num_years` of post-burn-in
- Total simulation = `2 × max_age + num_years`
- The burn-in serves two purposes: (1) flush all founders (mother_id = 0, father_id = 0), which requires max_age years, and (2) let the IBM's age structure settle from the Leslie-derived initial distribution to its true stochastic equilibrium, which requires ~2× max_age
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
- **Important implementation detail:** mother_rows (S1 females who give birth) is identified BEFORE new-mature females enter the cycle, so newly mature females NEVER breed in their first year of maturity. This is biologically appropriate (one-year gestation delay) but creates a mismatch with the Leslie matrix which applies π₂ fecundity at all mature ages.
- Stationary distribution: π₁ = π₂ = ψ₃/(2ψ₃ + 1 − ψ₂), π₃ = (1−ψ₂)·π₁/ψ₃

### 7. Maturity Specification
- **`maturity_age`** can be:
  - An integer: knife-edged maturity at that age for both sexes (backward compatible)
  - A numeric vector of length `max_age + 1`: cumulative probability of being mature at each age (ogive CDF)
  - A list with `female` and `male` elements, each of which can be an integer or ogive vector
- At birth, each individual is assigned a personal `mat_age` sampled from the ogive PMF
- This allows sex-specific maturity (e.g., logistic ogive for females, knife-edged for males)
- The ogive is also used in the Leslie matrix fecundity row for the analytical s0 estimate

### 8. Infertility
- **`infertility`** = scalar (both sexes) or `c(female, male)`
- At birth, each individual draws a permanent `fertile` flag (Bernoulli trial)
- Infertile females never enter the breeding cycle (breed_state stays NA)
- Infertile males are excluded from mating pools
- Leslie matrix scales fecundity by `(1 − infertility_f)`

### 9. Pod / Superpod Social Structure
- **Pods** = family groups of `pod_size` individuals (initialisation unit)
- **Superpods** = communities of `superpod_size` pods (mating + sampling units)
- Pod-to-superpod mapping is fixed at initialisation
- Offspring inherit mother's pod and superpod

### 10. Three Levels of Stickiness
- **`stickiness_year`** (in `simulate.pop`): between-year superpod fidelity. Affects mating pools. Sex-specific (scalar or `c(female, male)`).
- **`stickiness_trip`** (in `sample.pop`): between-trip reshuffling within a year. Sex-specific.
- **`stickiness_set`** (in `sample.pop`): between-set reshuffling within a trip. Controls duplication rate. Sex-specific.
- Movers go to a random pod in a DIFFERENT superpod (emigration, not within-community reshuffling)
- Within a trip, the vessel encounters the SAME superpod across all sets. Stickiness controls whether individuals leave that community between sets.

### 11. Cow-Calf Dynamics
- **`weaning_age`** (integer or NULL): age at which a calf becomes independent
- Below `weaning_age`: calves excluded from reshuffling, follow mother's pod/superpod
- Also serves as the threshold for stickiness-based reshuffling eligibility
- Orphaned calves (mother dead) stay in their current pod
- Active in all three functions (calculate_s0, simulate_pop, sample_pop)

### 12. Mating Modes
- **`"random"`**: any mature, fertile male in the superpod may sire offspring. `num_mates` controls polyandry (how many males each female mates with).
- **`"strong_bull"`**: persistent bull registry — one mature male per superpod sires all offspring. New bulls are elected **randomly** from mature males (not oldest), giving realistic multi-year tenure. Bulls serve until death.
- **`max_females`** (NULL or integer): per-year cap on how many females a single male can mate with. Enforced by reassigning excess offspring to other males in the same superpod.
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
  male_behavior = NULL, max_females = NULL, weaning_age = NULL,
  check_interval = 10L, growth_tol = 0.001, stable_required = 5L,
  max_windows = 100L
)
```

**Output:** List with 20 elements — calibration results (`s0`, `survival`, `s0_leslie`, `final_N`, `years_simulated`) plus all shared parameters.

**Input validation:** Checks types, lengths, ranges, and consistency (e.g., `pod_size` requires `superpod_size`, ogive must be non-decreasing CDF, etc.).

### `simulate.pop()`

```r
simulate.pop(sim_config, num_years, sample_years = NULL)
```

**Output:** List with:
- `pop_summary` — data.table: `year, sex, age, N` for all years
- `snapshots` — named list of population data.tables at snapshot years (columns: `id, birth_year, age, sex, mat_age, mother_id, father_id, breed_state, fertile, population, pod, superpod`)
- `pod_to_sp` — integer vector mapping pod → superpod
- `sim_config` — passed through for `sample.pop()`

**Input validation:** Checks `sim_config` has required fields, `num_years` is positive, etc.

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

**Output:** data.table with `id, birth_year, age, sex, mother_id, father_id, population, pod, superpod, year, trip, set`. May contain duplicates (same individual in multiple sets). Internal columns (`mat_age`, `fertile`, `breed_state`) are dropped.

**Input validation:** Checks types, ranges, and `superpod_pool` length matches `n_trips`.

---

## Example Usage

```r
library(data.table)
source("R/calculate_s0.R")
source("R/simulate_population.R")
source("R/sample_pop.R")

# Sex-specific maturity ogive example
ogive_f <- plogis(0:30, location = 9, scale = 1.5)

result <- calculate.s0(
  max_age = 30, survival = c(0.65, rep(0.96, 30)),
  pop_size = 50000,
  maturity_age = list(female = ogive_f, male = 10L),
  litter_size = 1,
  psi_2 = 0.1, psi_3 = 0.7,
  infertility = c(0.10, 0.05),
  num_mates = 1L,
  pod_size = 20, superpod_size = 10,
  stickiness_year = c(0.95, 0.85),
  male_behavior = "random", max_females = 50, weaning_age = 3L
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

## Test Results (August 21, 2026)

### Original tests (Aug 20)
✓ **Bundled config:** All shared params passed through correctly
✓ **Burn-in:** max_age-year burn-in flushes all founders; zero founders in snapshots
✓ **Snapshots:** Stored at correct years; full population with all columns
✓ **Pop summary:** Covers all years (burn-in + post)
✓ **Markov breeding states:** Observed proportions match expected stationary distribution
✓ **Cow-calf dynamics:** All calves match mother's superpod (0 mismatches)
✓ **Stickiness effect:** Higher stickiness_set → more duplication
✓ **Per-trip sampling:** Exactly sample_size per trip distributed across sets
✓ **Superpod pools:** Trips correctly restricted to specified superpods
✓ **Sex-specific stickiness:** Accepts `c(female, male)` vector
✓ **Strong bull mating:** Fewer unique fathers vs random mode
✓ **Markov transition bug caught and fixed:** indices must be computed before applying state changes

### New features (Aug 21)
✓ **Maturity ogive:** Knife-edged (integer), CDF vector, and list(female, male) all work
✓ **Sex-specific maturity:** Female logistic ogive + male knife-edged produces expected mat_age distributions in snapshots
✓ **Infertility:** 10% female / 5% male rates confirmed in snapshots
✓ **Infertility effect on s0:** Higher s0 needed to compensate for fewer breeders
✓ **max_females cap:** Enforced per-year; no father exceeds cap within a single breeding year
✓ **Strong bull random election:** Mean tenure ~7.7 years, max 23 years (vs ~1 year with old oldest-male approach)
✓ **Input validation:** 9/9 common error types caught with clear messages (wrong lengths, out-of-range values, missing config fields, inconsistent parameters)
✓ **Backward compatibility:** All existing call patterns work unchanged
✓ **Annotations:** Thorough code annotations added to all three scripts

### Population decline investigation & fix (Aug 21)
⚠ **Systematic decline confirmed:** 5 replicate runs with old defaults showed consistent negative growth (mean r ≈ -0.0006/yr, ~4% decline over 80 years)
⚠ **Not demographic stochasticity:** Decline persisted (and worsened) at N=200K
✓ **Root cause identified & fixed:** Leslie matrix assigned π₂ fecundity at first mature age, but IBM newly mature females enter at S3 and cannot breed until the following year (0% S2 at age 8 vs 30.4% expected). Fixed by shifting ogive right by one year in the Leslie matrix fecundity row.
✓ **Survival indexing bug fixed:** Sub-diagonal used `survival[i+1]` instead of `survival[i]` — off by one. Didn't affect flat survival vectors but would have produced wrong results with age-varying survival.
✓ **Default tolerances tightened:** check_interval=5→10, growth_tol=0.005→0.001, stable_required=3→5. With the corrected Leslie estimate, bisection now converges to s0≈0.590 (vs old 0.574) in ~100–220 years of calibration. Population growth is near-zero.
✓ **Burn-in doubled:** Changed from `max_age` to `2 × max_age` years. The first ~60 years show a transient as the age structure adjusts from the Leslie equilibrium to the IBM's true stochastic equilibrium. With `2 × max_age` burn-in, this transient is fully absorbed before snapshots.

---

## Known Issues / Open Questions

1. **Residual s0 drift** — After fixes, the population shows a very slight drift (typically < ±0.03%/yr). This is within the convergence tolerance (0.001) and much improved from the old -0.06%/yr.

---

## Next Steps

1. **CKMR kin-pair identification** — Build a function to find parent-offspring, full-sibling, and half-sibling pairs in samples and compare to ground-truth pedigree
4. **Sensitivity analysis** — Vary stickiness_set across a range and plot unique-sample yield per trip
5. **Scaling test** — Run at 500K–1M individuals to benchmark performance and memory
6. **End-to-end pipeline** — Connect `calculate.s0()` → `simulate.pop()` → `sample.pop()` → CKMR validation

### Deferred Features
- Length-based maturity/survival (partial implementation in legacy `helper_functions.R`)
- Multi-population structure with movement
- `prime_bull_age` parameter for controlling strong bull election eligibility age
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
