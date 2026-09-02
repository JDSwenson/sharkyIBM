# sharkyIBM Documentation Roadmap

This document maps parameters and concepts across the three main functions and their documentation.

---

## Function Documentation Structure

All three functions (`create.stable.pop()`, `simulate.pop()`, `sample.pop()`) now have comprehensive `@details` sections in their roxygen2 docstrings. You can access them with `?function_name`.

---

## create.stable.pop()

**Purpose**: Calibrate population parameters to stability; bundle configuration for simulation.

**@details Sections**:
1. **Calibration Modes** — Explains s0 calibration vs. density dependence modes
   - When to use each mode
   - Leslie matrix phase (analytical)
   - Bisection phase (simulation-based refinement)
   - Default tolerances and typical calibration time

2. **Two-Phase Approach** — Deep dive into Leslie matrix and bisection
   - Phase 1: Leslie matrix construction and eigenvalue solution
   - Phase 2: Iterative bisection logic (growth rates, tolerance windows)

3. **Initial Population Structure** — How the population is initialized
   - Individual attributes (ID, sex, age, maturity, fertility)
   - Breeding state initialization
   - Founder marking (mother_id=0, father_id=0)

4. **Markov Breeding Cycle** — Breeding state parameterization
   - psi_nurse vs. psi_rest conception rates
   - weaning_age effects on S2 duration
   - Stationary distribution calculation

5. **Social Structure Initialization** — Pods, superpods, deterministic mapping
   - pod_size, superpod_size
   - pod_to_sp vector creation

6. **Mating and Fertility** — Framework for num_mates, male_behavior, max_females, infertility

7. **Density Dependence Parameters** — theta_shift, psi_nurse_K, psi_rest_K, z_pt, dd_max

8. **Output Structure** — Complete list of returned parameters (bundled config)

**Key Parameters Documented**:
- max_age, survival, pop_size, litter_size
- maturity_age (knife-edged, ogive, sex-specific)
- psi_nurse, psi_rest, weaning_age
- num_mates, male_behavior, max_females, infertility, female_fraction
- pod_size, superpod_size, stickiness_year
- density_dependence, z_pt, dd_max
- check_interval, growth_tol, stable_required, max_windows

---

## simulate.pop()

**Purpose**: Run forward-time individual-based simulation with full parentage tracking.

**@details Sections**:
1. **Life History Parameters** — max_age, survival, pop_size, litter_size
   - How survival indexing works (survival[age + 1])
   - Role of survival[0] (s0) in each calibration mode
   - Litter size mechanics (Poisson with minimum 1)
   - Only first offspring tracked as dependent calf

2. **Maturity and Breeding Parameters** — maturity_age, psi_nurse, psi_rest, num_mates, female_fraction, infertility
   - Three maturity specifications (knife-edged, ogive, sex-specific)
   - Personal mat_age drawn at birth
   - How conception probabilities change with density dependence
   - infertility effects on Leslie matrix fecundity

3. **Markov Breeding Cycle with Calf-Survival Coupling** — Detailed state transitions
   - S1 (pregnant) → S2 (with calf)
   - S2 aging and conception (depends on calf survival)
   - S2 → S1 or S3 transitions
   - S3 → S1 with high probability
   - Newly mature females start at S3 (never breed first year)
   - calf_id tracking for dependent calf
   - s2_year internal column (calf survival by mother's age at conception)

4. **Social Structure** — pod_size, superpod_size, stickiness_year, weaning_age
   - Fixed hierarchical structure (pods nested in superpods)
   - Dependent calves follow mother's pod/superpod
   - stickiness_year mechanics (ages >= weaning_age eligible, emigrate to different superpod)
   - Sex-specific stickiness (c(female, male))

5. **Mating Systems** — male_behavior, max_females, num_mates
   - "random": polyandry from superpod pool
   - "strong_bull": persistent bull per superpod (randomly elected)
   - max_females cap enforcement
   - Cross-superpod fallback

6. **Density Dependence (Pella-Tomlinson)** — Full mathematical explanation
   - Conception rate adjustment equation: psi(t) = psi_K + delta_max * [1 - D(t)^z]
   - Logit-scale shift (preserves odds ratios)
   - Depletion D(t) = N_1+(t) / K_1+
   - MNPL at D ≈ 0.6 with z = 2.39 (IWC convention)
   - When density_dependence = FALSE, constant conception rates

7. **Automatic Burn-in and Snapshots** — 2*max_age burn-in rationale, snapshot mechanics
   - Founders flushed (~max_age years)
   - Age structure equilibration (~max_age years)
   - pop_summary covers all years
   - Snapshots exclude internal columns (s2_year)

**Key Parameters Documented**: All life history, demographic, breeding, social, mating, and DD parameters from sim_config

---

## sample.pop()

**Purpose**: Apply hierarchical sampling (trips/sets) to population snapshots with social reshuffling.

**@details Sections**:
1. **Hierarchical Sampling Structure** — Three-level hierarchy (year → trip → set)
   - Within each trip: single superpod (or reshuffled between sets/trips)
   - Realistic purse seine operation emulation

2. **Sampling Modes** — sampling = "random" vs. "superpod"
   - Random: uniform sampling from entire population
   - Superpod: sample from targeted community (more realistic)

3. **Stickiness Parameters** — stickiness_trip, stickiness_set
   - Between-trip fidelity within a year
   - Between-set fidelity within a trip
   - Effects on duplication rates
   - Sex-specific options

4. **Sampling Modes: sample_per** — "set" vs. "trip"
   - "set": sample_size per set (total = sample_size * n_sets)
   - "trip": sample_size total per trip (distributed across sets)

5. **Spatial Structure: superpod_pool** — NULL vs. list of pools per trip
   - Enforce spatial heterogeneity
   - Useful for regional sampling designs

6. **Cow-Calf Dynamics** — Dependent calves follow mother until weaning_age

7. **Output Structure** — Columns, duplicate detection, exclusions
   - One row per catch event (individual × set)
   - Deduplication by id for unique counts
   - Internal columns excluded

**Key Parameters Documented**: All sampling design parameters + social structure effects

---

## Parameter Cross-Reference

| Parameter | Defined in | Used in | Documentation |
|-----------|-----------|---------|---|
| max_age | create.stable.pop() param | simulate.pop(), sample.pop() | create_stable_pop.R §8; simulate_population.R §1, §8 |
| survival | create.stable.pop() param | simulate.pop() | create_stable_pop.R §2; simulate_population.R §1 |
| pop_size | create.stable.pop() param | simulate.pop() | create_stable_pop.R §4; simulate_population.R §2 |
| litter_size | create.stable.pop() param | simulate.pop() | create_stable_pop.R §3; simulate_population.R §1 |
| maturity_age | create.stable.pop() param | simulate.pop() | create_stable_pop.R §3; simulate_population.R §2 |
| psi_nurse | create.stable.pop() param | simulate.pop() | create_stable_pop.R §5; simulate_population.R §3, §6 |
| psi_rest | create.stable.pop() param | simulate.pop() | create_stable_pop.R §5; simulate_population.R §3, §6 |
| num_mates | create.stable.pop() param | simulate.pop() | create_stable_pop.R §6; simulate_population.R §5 |
| female_fraction | create.stable.pop() param | simulate.pop() | create_stable_pop.R §2; simulate_population.R §2 |
| infertility | create.stable.pop() param | simulate.pop() | create_stable_pop.R §6; simulate_population.R §2, §6 |
| pod_size | create.stable.pop() param | simulate.pop(), sample.pop() | create_stable_pop.R §7; simulate_population.R §4; sample_pop.R §1 |
| superpod_size | create.stable.pop() param | simulate.pop(), sample.pop() | create_stable_pop.R §7; simulate_population.R §4; sample_pop.R §1 |
| stickiness_year | create.stable.pop() param | simulate.pop() | create_stable_pop.R §8; simulate_population.R §4 |
| male_behavior | create.stable.pop() param | simulate.pop() | create_stable_pop.R §6; simulate_population.R §5 |
| max_females | create.stable.pop() param | simulate.pop() | create_stable_pop.R §6; simulate_population.R §5 |
| weaning_age | create.stable.pop() param | simulate.pop(), sample.pop() | create_stable_pop.R §5; simulate_population.R §3, §4, §6; sample_pop.R §6 |
| density_dependence | create.stable.pop() param | simulate.pop() | create_stable_pop.R §1, §7; simulate_population.R §6 |
| z_pt | create.stable.pop() param | simulate.pop() | create_stable_pop.R §7; simulate_population.R §6 |
| dd_max | create.stable.pop() param | simulate.pop() | create_stable_pop.R §7; simulate_population.R §6 |
| stickiness_set | sample.pop() param | sample.pop() | sample_pop.R §3, §7 |
| stickiness_trip | sample.pop() param | sample.pop() | sample_pop.R §3, §7 |
| sample_per | sample.pop() param | sample.pop() | sample_pop.R §4 |
| sampling | sample.pop() param | sample.pop() | sample_pop.R §2 |
| superpod_pool | sample.pop() param | sample.pop() | sample_pop.R §5 |

---

## How to Use This Documentation

1. **For a parameter you've set**: Look it up in the table above, then open the function's documentation (`?function_name`) and navigate to the @details section.

2. **For a conceptual question** (e.g., "How does the Markov breeding cycle work?"): 
   - Check `?simulate.pop` and look for "Markov Breeding Cycle with Calf-Survival Coupling"
   - Or check `?create.stable.pop` and look for "Markov Breeding Cycle"

3. **For sampling design questions** (e.g., "What does stickiness_set do?"):
   - Check `?sample.pop` and look for "Stickiness Parameters"

4. **For calibration/stability questions** (e.g., "What's the difference between s0 calibration and density dependence?"):
   - Check `?create.stable.pop` and look for "Calibration Modes"

---

## Files Modified (Sep 2, 2026)

- `R/create_stable_pop.R` — Added comprehensive 8-section @details
- `R/simulate_population.R` — Added comprehensive 8-section @details
- `R/sample_pop.R` — Added comprehensive 7-section @details
- `DOCUMENTATION_ROADMAP.md` (this file) — Created for easy navigation

---

## Next Steps

- Generate `.Rd` files via `roxygen2::roxygenise()`
- Verify all cross-references in roxygen2 metadata
- Consider adding a vignette with worked examples linking all three functions
