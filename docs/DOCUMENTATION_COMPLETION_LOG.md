# sharkyIBM Documentation Completion Log

**Date**: September 2, 2026  
**Status**: ✓ Complete and passing `devtools::check()`

---

## Summary

Added comprehensive `@details` sections to all three main functions in the sharkyIBM package. Every parameter flowing through from `create.stable.pop()` to `simulate.pop()` and `sample.pop()` is now documented with clear explanations of its role, interpretation, and effects.

**Final package check**: 0 errors, 0 warnings, 1 note (about top-level TODO files, which is expected and harmless).

---

## Documentation Additions by Function

### 1. create.stable.pop.R

**File**: `R/create_stable_pop.R`

**Added**: Comprehensive 8-section `@details` documentation covering:

1. **Calibration Modes** (s0 vs. density dependence)
   - When to use each mode
   - User responsibilities (placeholder vs. complete survival vector)
   - Leslie matrix + bisection pipeline

2. **Two-Phase Approach** (Leslie matrix + iterative refinement)
   - Leslie matrix construction and analytical s0/theta_shift estimate
   - Bisection window mechanics (check_interval, growth_tol, stable_required)
   - Typical calibration time (100–220 years)

3. **Initial Population Structure**
   - Individual attribute assignment (ID, sex, age, mat_age, fertility)
   - Breeding state initialization from stationary distribution
   - Founder marking for flush detection

4. **Markov Breeding Cycle**
   - psi_nurse vs. psi_rest conception rates
   - weaning_age effects on S2 duration
   - Stationary distribution computation

5. **Social Structure Initialization**
   - pod_size, superpod_size deterministic assignment
   - pod_to_sp vector construction

6. **Mating and Fertility**
   - num_mates, male_behavior ("random" vs. "strong_bull")
   - max_females cap enforcement
   - infertility effects on fecundity

7. **Density Dependence Parameters**
   - theta_shift, psi_nurse_K, psi_rest_K calibration
   - z_pt (Pella-Tomlinson shape, IWC convention)
   - dd_max (compensation strength)
   - K_1+ calculation

8. **Output Structure**
   - Complete list of returned config parameters
   - Organization into calibration, demographics, breeding, social, DD categories

**Parameters Documented**: max_age, survival, pop_size, litter_size, maturity_age, psi_nurse, psi_rest, num_mates, female_fraction, infertility, pod_size, superpod_size, stickiness_year, male_behavior, max_females, weaning_age, density_dependence, z_pt, dd_max, check_interval, growth_tol, stable_required

---

### 2. simulate_population.R

**File**: `R/simulate_population.R`

**Added**: Comprehensive 8-section `@details` documentation covering:

1. **Life History Parameters** (from sim_config)
   - max_age removal mechanics (breed at max_age, then removed)
   - survival indexing (survival[age + 1])
   - pop_size as carrying capacity
   - litter_size Poisson mechanics with minimum 1
   - Dependent calf tracking (first offspring only)

2. **Maturity and Breeding Parameters**
   - maturity_age specifications (knife-edged, ogive, sex-specific list)
   - Personal mat_age draw at birth
   - psi_nurse vs. psi_rest conception probabilities
   - Density dependence adjustment (psi_nurse_K, psi_rest_K)
   - num_mates polyandry
   - female_fraction sex ratio
   - infertility effects on breeding pool and Leslie fecundity

3. **Markov Breeding Cycle with Calf-Survival Coupling**
   - Five detailed state transition steps (S1→S2, S2 aging, S2→S1/S3, S3→S1, newborns→S3)
   - calf_id dependent calf tracking
   - s2_year internal column (calf survival by mother's conception age)
   - Newly mature females start at S3 (never breed first year)
   - weaning_age effects on S2 duration and social following

4. **Social Structure**
   - pod_size, superpod_size hierarchical nesting
   - pod_to_sp fixed mapping
   - stickiness_year between-year superpod fidelity
   - Eligible ages (>= weaning_age) and emigration mechanics
   - Dependent calves follow mother's pod/superpod

5. **Mating Systems**
   - male_behavior: "random" (polyandry) vs. "strong_bull" (persistent bull)
   - max_females per-male cap enforcement
   - Cross-superpod fallback

6. **Density Dependence (Pella-Tomlinson)**
   - Full equation: psi(t) = psi_K + delta_max * [1 - D(t)^z]
   - Logit-scale shift mechanics (preserves odds ratios)
   - Depletion D(t) = N_1+(t) / K_1+
   - z_pt convention (2.39 → MNPL at 0.6K)
   - When density_dependence = FALSE, constant conception rates

7. **Automatic Burn-in and Snapshots**
   - 2*max_age burn-in rationale (founders + equilibration)
   - pop_summary coverage (all years, burn-in + post)
   - Snapshot column exclusions (s2_year excluded)

8. **Output Structure** (already in @return, now cross-referenced in @details)

**Parameters Documented**: All parameters from sim_config + num_years, sample_years

---

### 3. sample_pop.R

**File**: `R/sample_pop.R`

**Added**: Comprehensive 7-section `@details` documentation covering:

1. **Hierarchical Sampling Structure**
   - Three-level hierarchy: year → trip → set
   - Single superpod per trip (with optional reshuffling)
   - Purse seine operation emulation

2. **Sampling Modes** (sampling = "random" vs. "superpod")
   - Random: uniform from entire population (unrealistic)
   - Superpod: targeted community sampling (realistic for marine mammals)

3. **Stickiness Parameters** (stickiness_trip, stickiness_set)
   - Between-trip fidelity within a year
   - Between-set fidelity within a trip
   - Duplication rate effects (high stickiness → more duplicates)
   - Sex-specific options via c(female, male)

4. **Sampling Modes: sample_per** ("set" vs. "trip")
   - "set": sample_size per set (total = sample_size * n_sets)
   - "trip": sample_size total per trip (distributed across sets)

5. **Spatial Structure: superpod_pool** (NULL vs. list)
   - NULL: all superpods available to all trips
   - List: enforces spatial heterogeneity (e.g., north vs. south trips)

6. **Cow-Calf Dynamics**
   - Dependent calves follow mother's pod/superpod
   - Weaning-age threshold for independence

7. **Output Structure and Deduplication**
   - Column list and internal exclusions
   - Duplicate detection (same individual in multiple sets)
   - Unique count via \code{unique(samples, by = "id")}

**Parameters Documented**: n_trips, n_sets, sample_size, sample_per, sampling, stickiness_set, stickiness_trip, superpod_pool

---

## Cross-Reference Documentation

Created **DOCUMENTATION_ROADMAP.md** with:
- Function-by-function section listing, @details subsections, and parameters documented in each
- Comprehensive parameter cross-reference table showing:
  - Parameter name
  - Where it's defined (which function)
  - Where it's used (which downstream functions)
  - Where it's documented in @details (file and section)
- Usage guidance for finding documentation

---

## Quality Assurance

### Roxygen2 Processing
- ✓ All three `.R` files successfully processed by roxygen2
- ✓ `.Rd` files generated correctly in `man/` directory
- ✓ Fixed LaTeX/Greek letter rendering issues by using text notation and `\code{}` formatting

### Package Checks
- ✓ `devtools::check()` runs successfully
- ✓ 0 errors, 0 warnings, 1 note (about non-standard top-level files, expected)
- ✓ All tests pass
- ✓ No missing documentation entries
- ✓ Rd cross-references valid

### Documentation Accessibility
- All functions callable via `?function_name` in R
- `@details` sections display cleanly with subsection headers (## Markdown style)
- Code examples and parameter usage clear and consistent

---

## Files Modified

| File | Change | Lines Added |
|------|--------|-------------|
| R/create_stable_pop.R | Added @details (8 subsections) | ~190 |
| R/simulate_population.R | Added @details (8 subsections) | ~170 |
| R/sample_pop.R | Added @details (7 subsections) | ~130 |
| DOCUMENTATION_ROADMAP.md | Created (new file) | ~150 |
| DOCUMENTATION_COMPLETION_LOG.md | Created (this file) | ~200 |

**Total documentation additions**: ~840 lines

---

## Key Features of Documentation

### Clarity and Completeness
- Every parameter flowing through the package is explained
- Technical concepts (Markov states, Leslie matrix, Pella-Tomlinson, density dependence) covered at appropriate depth
- Connections between functions clearly explained

### Consistency
- Same parameter used in multiple functions is explained consistently
- Code formatting consistent across all three functions (`\code{}` for parameter names and R objects)
- References to other functions use standard roxygen2 link format

### Accessibility
- Hierarchical organization with clear subsection headers
- Mix of technical depth and plain English (no jargon-heavy prose)
- Equations presented in readable `\preformatted{}` and `\itemize{}` formats (no LaTeX rendering errors)

### Maintainability
- DOCUMENTATION_ROADMAP.md serves as a guide for future maintainers
- Cross-reference table makes it easy to find where parameters are documented
- Clear file/section notation in roadmap table

---

## Next Steps (Parked)

1. Generate package website via pkgdown (includes searchable documentation)
2. Create vignettes with worked examples linking all three functions
3. Add example code blocks to @examples sections (currently empty)
4. Consider expanding DOCUMENTATION_ROADMAP.md into a formal Getting Started guide

---

## Conclusion

All three main functions now have comprehensive, checked, and passing documentation. Every parameter that flows through the package is explained at appropriate length and technical depth. The DOCUMENTATION_ROADMAP serves as a user-friendly navigation guide.

**Status**: ✅ **READY FOR PUBLICATION**
