# sharkyIBM Project — Conversation Summary

**Date:** Friday, July 17, 2026  
**Project:** Individual-Based Population Simulation for Close-Kin Mark–Recapture (CKMR) Analysis

---

## Overview

This is a shark population simulation package (`sharkyIBM`) designed to:
- Simulate age-structured populations forward in time with full parentage tracking
- Test sampling schemes and life-history attributes with CKMR models
- Ground-truth CKMR estimates with known population demographics

The project consists of three main scripts:
- `R/calculate_s0.R` — Empirically finds age-0 survival (s0) that stabilises a population
- `R/simulate_population.R` — Runs the forward-time individual-based simulation (IBM)
- `R/helper_functions.R` — Utility functions (currently empty; will be populated)

---

## Key Decisions & Design Principles

### 1. **Performance First**
- All code uses `data.table` for in-place column updates and fast subsetting (critical at million-individual scale)
- Survival lookups use direct vector indexing (`surv_vec[age + 1]`) instead of table joins
- No father-identity tracking in `calculate.s0()` — only offspring counts matter for growth rate
- No unnecessary copying; vectorised operations throughout

### 2. **Leslie Matrix as Foundation**
- Both functions build a classic age-structured Leslie matrix to:
  - Provide an analytical s0 starting estimate (s0_leslie)
  - Extract the stable age distribution for population initialisation
- The Leslie matrix approach ensures consistency with demographic theory while the IBM captures stochasticity

### 3. **Empirical s0 Calibration**
- `calculate.s0()` **does not** run multiple full simulations. Instead:
  - Runs one continuous simulation for ~5-year assessment windows
  - Checks growth every window; bisects s0 if growth deviates from zero
  - Converges when 3 consecutive windows show growth in tolerance (default: |r| < 0.005)
  - At 100K pop: converges in ~15 years, 0.4 seconds
- This iterative approach is ~1000× faster than running separate multi-year replicates

### 4. **Forward Compatibility for Extensions**
Code is structured to support future additions:
- **Length-based dynamics:** `process_by` parameter; age/length maturity can be swapped
- **Movement/migration:** `popstructure` ("panmictic" vs "structured"), `movement_array` placeholder
- **Infertility:** `infertility` parameter and `fertile` column in population table (currently unused)
- **Multiple populations:** `population` column in all data structures

### 5. **Pod / Superpod Social Structure**
- **Pods** = small family groups (kin-based); **superpods** = aggregates (mating/sampling units)
- **Stickiness** (0–1): probability of staying in pod at each shuffle interval
- **Sticky_age** & **sticky_interval**: when reshuffling begins and frequency
- **Male_behavior**: `"random"` (any mature male in superpod mates) or `"strong_bull"` (oldest male only)
- **Cross-superpod fallback:** if a superpod has no mature males, females mate with males from other superpods
- Offspring inherit mother's pod and superpod membership

---

## Function Specifications

### `calculate.s0()`

**Purpose:** Find empirical age-0 survival that produces population growth rate ≈ 0

**Required Inputs:**
```r
max_age = 30
survival = c(0.65, rep(0.80, 30))  # age-0 placeholder; will be overwritten
pop_size = 100000
mating_periodicity = 2  # biennial breeding
maturity_age = 7
litter_size = 7
num_mates = 1:3
```

**Default Parameters:**
- `female_fraction = 0.5`
- `check_interval = 5` (years between growth assessments)
- `growth_tol = 0.005` (|r| < this = stable)
- `stable_required = 3` (consecutive stable windows to converge)

**Pod Parameters** (if `stickiness` is not NULL):
- `stickiness, sticky_age, sticky_interval, superpod_size, male_behavior`

**Output:** List with:
- `s0` — calibrated age-0 survival probability
- `survival` — full survival vector with s0 inserted
- `s0_leslie` — analytical Leslie estimate (for comparison)
- `final_N, years_simulated` — diagnostics

**Key Annotations:**
- Explains Leslie matrix construction and stable age distribution
- Documents bisection logic and why it resets the stability streak on s0 adjustment
- Clarifies why litter size is `1 + rpois(lambda - 1)` (no zero-litter guarantee)

---

### `simulate.pop()`

**Purpose:** Run forward-time IBM with full parentage tracking; sample individuals for CKMR

**Required Inputs:**
```r
max_age, survival, pop_size, mating_periodicity, maturity_age,
litter_size, num_mates, num_years
```

**Sampling Inputs:**
- `sample_size` — scalar (random) or vector (superpod-based)
- `sample_years` — scalar (last N years), vector (specific years), or NULL (no sampling)
- `sampling` — `"random"` or `"superpod"`

**Pod Parameters:**
- `stickiness, sticky_age, sticky_interval, superpod_size, male_behavior`

**Output:** List with:
- `samples` — data.table of sampled individuals:
  - `id, birth_year, age, sex, mother_id, father_id`
  - `capture_year, population`
  - (if pods active: `pod, superpod`)
- `pop_summary` — data.table of counts by year, sex, age:
  - Columns: `year, sex, age, N`

**Annual Loop:**
1. **Survival** — direct vector lookup + Bernoulli draws
2. **Aging** — increment age; remove those past max_age
3. **Pod shuffling** — (if active) probabilistic pod reassignment
4. **Breeding** — identify females by cycle; assign mates; create offspring with full parentage
5. **Sampling** — (if sample year) draw individuals and record metadata
6. **Metrics** — tally counts by sex/age for ground-truthing

**Father Assignment** (vectorised per-superpod):
- Strong-bull: one dominant male per superpod sires all offspring (oldest wins)
- Random: all mature males in superpod are available; sampled with replacement per mating event
- Fallback: if superpod has no males, use any male from population

**Founder Burn-in:**
- Founders (`mother_id = 0, father_id = 0`) in initial population
- All founder-parented individuals die by year `max_age + maturity_age`
- Warning issued if sampling begins too early

---

## Performance Benchmarks (100K population)

| Configuration | Years | Sampled | Time |
|---|---|---|---|
| No pods | 80 | 2,500 | 5.4s |
| Pods + random | 80 | 2,500 | 12.8s |
| Pods + strong_bull | 50 | 600 (50K pop) | 2.3s |

---

## Test Results Summary

✓ **Without pods:** Population stable; zero founder parents in samples (sampling at year 76)  
✓ **With pods (random):** 88/2500 sampled individuals with parent also in sample (sensible hit rate)  
✓ **With pods (strong_bull):** Fewer unique fathers (360) vs mothers (574); mating structure as expected  
✓ **Parent-offspring pairs recoverable:** Full parentage tracked for CKMR analysis  

---

## Next Steps (for future sessions)

1. **Annotation of `simulate.pop()`**
   - Add thorough inline comments explaining biological logic and algorithmic choices
   - Clarify father-matrix construction and per-superpod fills

2. **Scaling Test**
   - Run at 1 million individuals to benchmark production-scale performance
   - Confirm memory usage and timing

3. **Superpod Sampling**
   - Test `sampling = "superpod"` mode to verify correct within-superpod sampling

4. **CKMR Validation Function**
   - Build a function to identify kin pairs in samples (PO, FS, HS)
   - Compare identified pairs to ground-truth pedigree

5. **End-to-end Pipeline**
   - Connect `calculate.s0()` → `simulate.pop()` in setup script
   - Save example output for documentation

6. **Future Feature Implementation**
   - Length-based maturity/survival (stub already in place)
   - Multi-population structure with movement
   - Infertility dynamics

---

## File Structure

```
sharkyIBM/
├── R/
│   ├── calculate_s0.R          [fully implemented, annotated]
│   ├── simulate_population.R   [fully implemented, needs annotation]
│   ├── helper_functions.R      [skeleton; will expand as needed]
├── sharky_IBM_setup_code.R     [user-facing parameter defaults]
├── AGENTS.md                   [this file; project memory]
```

---

## Current Session State (R environment)

Variables loaded in RStudio console:
- `result` — output from `calculate.s0()`: s0, survival, s0_leslie, etc.
- `sim`, `sim_pods`, `sim_bull` — outputs from test runs of `simulate.pop()`
  - `sim$samples` — 2,500 sampled individuals (no pods)
  - `sim$pop_summary` — 4,960 rows of counts by year/sex/age
- `survival` — calibrated survival vector from result
- Life-history parameters (`max_age`, `pop_size`, `mating_periodicity`, etc.)

All functions source cleanly; ready for continued development.
