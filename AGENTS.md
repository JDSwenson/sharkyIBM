# sharkyIBM Project — Conversation Summary

**Last updated:** Tuesday, September 1, 2026
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

1. **`R/create_stable_pop.R`** (`create.stable.pop()`) — Stabilizes an age-structured population via one of two modes: (a) s0 calibration (`density_dependence = FALSE`, default) or (b) density-dependent conception adjustment (`density_dependence = TRUE`). Returns a bundled config object with all shared parameters.
2. **`R/simulate_population.R`** (`simulate.pop()`) — Runs the forward-time IBM with automatic `2 × max_age`-year burn-in, Markov breeding cycle, optional Pella-Tomlinson density dependence, between-year superpod reshuffling, and parentage tracking. Returns population snapshots at user-specified years.
3. **`R/sample_pop.R`** (`sample.pop()`) — Takes snapshots from `simulate.pop()` and applies a hierarchical sampling scheme: trips → sets → samples, with stickiness-based reshuffling between sets and trips. Can be re-run cheaply with different sampling parameters without re-simulating.

### Parked/Legacy (not wired in)

- `R/helper_functions.R` — Older dplyr-based breeding functions with length-based growth (von Bertalanffy, MVN-correlated L_inf/K). Not called by current pipeline.
- `R/create_input_data.R` — Parked. Uses Euler-Lotka for s0, superseded by `create_stable_pop.R`.

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

### 6. Markov Breeding Cycle with Calf-Survival Coupling
- Each mature female carries a `breed_state`: S1 (pregnant), S2 (with dependent calf), S3 (resting)
- **Parameterization:** `psi_nurse` (conception while nursing, suppressed) and `psi_rest` (conception while resting/after calf death, unsuppressed). Replaces the former `psi_2`/`psi_3`.
- **Calf-survival coupling (matching dolphin_population_model.qmd):** In `simulate.pop()`, each S2 mother's conception probability depends on whether *her individual calf* survived:
  - Calf alive → conception prob = `psi_nurse` (lactational suppression)
  - Calf dead → conception prob = `psi_rest` (released from suppression)
  - This creates compensatory feedback: higher calf mortality → mothers breed sooner
- **Population-level equivalence:** The effective conception rate at each S2 year k is the mixture: ψ_k = ℓ_k·ψ_nurse + (1−ℓ_k)·ψ_rest, where ℓ_k = surv_vec[k] (calf survival). In `create.stable.pop()` calibration (no individual IDs), calf survival is drawn stochastically from surv_vec.
- **Variable-length S2:** S2 lasts up to `weaning_age` years (default 1 if NULL). Each year, if the calf is alive and below weaning_age, the mother stays in S2. At weaning_age or if calf dies, the mother transitions to S1 (conceive) or S3 (rest).
  - weaning_age=1: equivalent to old 3-state model (S2 lasts exactly 1 year)
  - weaning_age=2: matches QMD's 4-state cycle (WITH NEWBORN + WITH YEARLING)
- **Individual tracking in simulate.pop():** Each S2 mother has a `calf_id` column linking to her dependent calf's ID. Founder S2 mothers (calf_id=0) use stochastic calf survival via `s2_year`.
- **Stationary distribution:** Computed numerically by `breeding_stationary()` helper, which builds the (wa_breed+2)-state Markov transition matrix and solves for the left eigenvector. Used for Leslie matrix fecundity and breeding state initialization.
- Newly mature females enter at S3 (resting)
- **Important implementation detail:** All state indices (S1, S2, S3) must be identified BEFORE applying transitions to prevent double transitions within a single year
- **Important implementation detail:** mother_rows (S1 females who give birth) is identified BEFORE new-mature females enter the cycle, so newly mature females NEVER breed in their first year of maturity.

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

### 9. Density Dependence (Pella-Tomlinson)
- **Two stabilization modes** in `create.stable.pop()`:
  - `density_dependence = FALSE` (default): bisect s0 to stabilize (original behavior)
  - `density_dependence = TRUE`: user supplies complete survival curve (incl. realistic s0); function bisects a logit-scale offset (`theta_shift`) on psi_2/psi_3 to find the conception level that makes λ(K) = 1
- **DD mechanism in `simulate.pop()`:** each year, computes depletion D(t) = N_1+(t) / K_1+, then shifts conception log-odds by Δ(t) = dd_max × (1 - D(t)^z_pt)
  - At K: Δ = 0, calibrated psi_2_K/psi_3_K apply directly
  - Below K: Δ > 0, conception rates increase (compensation)
  - Above K: Δ < 0, conception rates decrease
- **Logit-scale shift** preserves the ratio between psi_2 and psi_3 (odds ratio invariant to depletion), same approach as the dolphin_population_model.qmd assessment model
- **K_1plus** = 1+ component of carrying capacity (age ≥ 1), derived from stable age distribution. Excludes age-0 to avoid feedback (newborn count is a direct function of breeding rate)
- **User-supplied psi_2/psi_3** are reference values; the calibration adjusts them via theta_shift. Adjusted values stored as `psi_2_K`, `psi_3_K` in the config
- **DD parameters:** `z_pt` (Pella-Tomlinson shape, default 2.39), `dd_max` (max logit shift, required when DD=TRUE)
- **DD NOT active during calibration** — calibration simulates at K where Δ = 0; DD only activates in `simulate.pop()`

### 10. Pod / Superpod Social Structure
- **Pods** = family groups of `pod_size` individuals (initialisation unit)
- **Superpods** = communities of `superpod_size` pods (mating + sampling units)
- Pod-to-superpod mapping is fixed at initialisation
- Offspring inherit mother's pod and superpod

### 11. Three Levels of Stickiness
- **`stickiness_year`** (in `simulate.pop`): between-year superpod fidelity. Affects mating pools. Sex-specific (scalar or `c(female, male)`).
- **`stickiness_trip`** (in `sample.pop`): between-trip reshuffling within a year. Sex-specific.
- **`stickiness_set`** (in `sample.pop`): between-set reshuffling within a trip. Controls duplication rate. Sex-specific.
- Movers go to a random pod in a DIFFERENT superpod (emigration, not within-community reshuffling)
- Within a trip, the vessel encounters the SAME superpod across all sets. Stickiness controls whether individuals leave that community between sets.

### 12. Cow-Calf Dynamics
- **`weaning_age`** (integer or NULL): age at which a calf becomes independent
- Below `weaning_age`: calves excluded from reshuffling, follow mother's pod/superpod
- Also serves as the threshold for stickiness-based reshuffling eligibility
- Orphaned calves (mother dead) stay in their current pod
- Active in all three functions (calculate_s0, simulate_pop, sample_pop)

### 13. Mating Modes
- **`"random"`**: any mature, fertile male in the superpod may sire offspring. `num_mates` controls polyandry (how many males each female mates with).
- **`"strong_bull"`**: persistent bull registry — one mature male per superpod sires all offspring. New bulls are elected **randomly** from mature males (not oldest), giving realistic multi-year tenure. Bulls serve until death.
- **`max_females`** (NULL or integer): per-year cap on how many females a single male can mate with. Enforced by reassigning excess offspring to other males in the same superpod.
- Cross-superpod fallback: if a superpod has no mature males, females mate with males from other superpods

---

## Function Specifications

### `create.stable.pop()`

```r
create.stable.pop(
  max_age, survival, pop_size, maturity_age, litter_size,
  psi_nurse = 0.1, psi_rest = 0.7,
  num_mates = 1L, female_fraction = 0.5, infertility = 0,
  pod_size = NULL, superpod_size = NULL, stickiness_year = NULL,
  male_behavior = NULL, max_females = NULL, weaning_age = NULL,
  density_dependence = FALSE, z_pt = 2.39, dd_max = NULL,
  check_interval = 10L, growth_tol = 0.001, stable_required = 5L,
  max_windows = 100L
)
```

**Output (DD=FALSE):** List with calibration results (`s0`, `survival`, `s0_leslie`, `final_N`, `years_simulated`, `density_dependence = FALSE`) plus all shared parameters (including `psi_nurse`, `psi_rest`).

**Output (DD=TRUE):** List with calibration results (`s0` [user-supplied], `survival` [unchanged], `s0_leslie = NA`, `theta_shift`, `theta_shift_leslie`, `psi_nurse_K`, `psi_rest_K`, `final_N`, `years_simulated`, `density_dependence = TRUE`, `z_pt`, `dd_max`, `K_1plus`) plus all shared parameters.

**Input validation:** Checks types, lengths, ranges, and consistency (e.g., `pod_size` requires `superpod_size`, ogive must be non-decreasing CDF, DD requires dd_max and non-zero psi values, etc.).

### `simulate.pop()`

```r
simulate.pop(sim_config, num_years, sample_years = NULL)
```

**Output:** List with:
- `pop_summary` — data.table: `year, sex, age, N` for all years
- `snapshots` — named list of population data.tables at snapshot years (columns: `id, birth_year, age, sex, mat_age, mother_id, father_id, breed_state, fertile, population, calf_id, pod, superpod`; `s2_year` is internal and excluded from snapshots)
- `pod_to_sp` — integer vector mapping pod → superpod
- `sim_config` — passed through for `sample.pop()`
- `depletion` — numeric vector of D(t) for each year (only when DD=TRUE)

**Input validation:** Checks `sim_config` has required fields (including `density_dependence`), `num_years` is positive, etc.

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
source("R/create_stable_pop.R")
source("R/simulate_population.R")
source("R/sample_pop.R")

# ── Example 1: s0 calibration (default, DD=FALSE) ──
ogive_f <- plogis(0:30, location = 9, scale = 1.5)

result <- create.stable.pop(
  max_age = 30, survival = c(0.65, rep(0.96, 30)),
  pop_size = 50000,
  maturity_age = list(female = ogive_f, male = 10L),
  litter_size = 1,
  psi_nurse = 0.1, psi_rest = 0.7,
  infertility = c(0.10, 0.05),
  num_mates = 1L,
  pod_size = 20, superpod_size = 10,
  stickiness_year = c(0.95, 0.85),
  male_behavior = "random", max_females = 50, weaning_age = 3L
)

# ── Example 2: density dependence (DD=TRUE) with Siler curve ──
# Siler U-shaped survival: M(a) = M_calf*exp(-b*a) + M_adult + senescent
ages <- 0:40
M <- 0.15 * exp(-0.5 * ages) + 0.04
M[ages >= 28] <- M[ages >= 28] +
  0.10 * (exp(0.3 * (ages[ages >= 28] - 28)) - 1) /
         (exp(0.3 * (40 - 28)) - 1)
surv_siler <- exp(-M)

result_dd <- create.stable.pop(
  max_age = 40, survival = surv_siler,
  pop_size = 50000,
  maturity_age = 9L,
  litter_size = 1,
  psi_nurse = 0.1, psi_rest = 0.7,
  pod_size = 20, superpod_size = 10,
  stickiness_year = 0.9,
  male_behavior = "random", weaning_age = 2L,
  density_dependence = TRUE, z_pt = 2.39, dd_max = 3.0
)

sim <- simulate.pop(result_dd, num_years = 50, sample_years = 5)

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

### Calf-survival-dependent breeding (Sep 1)
✓ **New parameterization:** `psi_nurse`/`psi_rest` replace `psi_2`/`psi_3`; calf survival couples to mother's breeding transition
✓ **breeding_stationary() helper:** Builds (w+2)-state Markov chain, solves for stationary distribution; reproduces old closed-form for wa_breed=1
✓ **Compensatory feedback verified:** Low calf survival (s0=0.4) → shorter calving interval (3.07 yr) vs high calf survival (s0=0.9) → 3.87 yr
✓ **Individual calf tracking:** S2 mothers correctly link to their calf via `calf_id`; all S2 mothers in snapshots have living calves; calf ages 0–1 for wa_breed=2
✓ **Variable-length S2:** weaning_age=2 gives 2-year dependency; weaning_age=NULL defaults to 1-year
✓ **Breeding state proportions match stationary distribution:** S1=0.246 vs 0.232, S2=0.400 vs 0.414, S3=0.354 vs 0.355
✓ **DD mode works with new parameterization:** theta_shift calibration converges; at-K psi_nurse_K/psi_rest_K stored
✓ **DD recovery from 50% K:** Population recovers to K within ~60 years (same as before)
✓ **sample.pop() compatibility:** Internal columns (s2_year) excluded from samples; calf_id excluded via keep_cols whitelist
✓ **Backward compatibility note:** `psi_2`/`psi_3` parameters are REMOVED; existing code must switch to `psi_nurse`/`psi_rest`

### Density dependence (Sep 1)
✓ **Renamed function:** `calculate.s0()` → `create.stable.pop()`, file `calculate_s0.R` → `create_stable_pop.R`
✓ **Dual-mode design:** `density_dependence = FALSE` bisects s0 (original); `density_dependence = TRUE` bisects theta_shift on psi_2/psi_3 logits
✓ **Leslie Phase 1:** Solves for theta_shift analytically via uniroot on Leslie matrix (DD=TRUE) or s0 (DD=FALSE)
✓ **IBM Phase 2:** Bisects theta_shift (DD=TRUE) or s0 (DD=FALSE) with same window/tolerance machinery
✓ **simulate.pop() DD support:** Yearly depletion computation, Pella-Tomlinson logit shift on psi_2/psi_3, depletion vector in output
✓ **Siler survival curve test:** U-shaped mortality (M_adult=0.04, M_calf_excess=0.15, b_juv=0.5, a_sen_onset=28, max_age=40) — calibration converges, population stable at K ≈ 50K with DD active
✓ **DD recovery test:** Population started at 50% K recovers smoothly to K within ~60 years
✓ **Backward compatibility:** DD=FALSE mode produces identical behavior to old calculate.s0()
✓ **sample.pop() compatibility:** Works with DD simulation output unchanged

---

## Known Issues / Open Questions

1. **Residual s0 drift** — After fixes, the population shows a very slight drift (typically < ±0.03%/yr). This is within the convergence tolerance (0.001) and much improved from the old -0.06%/yr.

---

## Next Steps

1. **CKMR kin-pair identification** — Build a function to find parent-offspring, full-sibling, and half-sibling pairs in samples and compare to ground-truth pedigree
4. **Sensitivity analysis** — Vary stickiness_set across a range and plot unique-sample yield per trip
5. **Scaling test** — Run at 500K–1M individuals to benchmark performance and memory
6. **End-to-end pipeline** — Connect `create.stable.pop()` → `simulate.pop()` → `sample.pop()` → CKMR validation

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
│   ├── create_stable_pop.R     [active — s0 calibration OR DD conception calibration + bundled config]
│   ├── simulate_population.R   [active — IBM with burn-in + Markov breeding + DD + snapshots]
│   ├── sample_pop.R            [active — hierarchical sampling (trips/sets/stickiness)]
│   ├── helper_functions.R      [legacy — not wired in; has length-based growth]
│   ├── create_input_data.R     [parked — superseded]
├── AGENTS.md                   [this file]
```
