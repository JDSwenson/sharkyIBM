# sharkyIBM: Species Generalization Roadmap

**Current target:** Eastern spinner dolphins (and cetaceans generally)

**Current generality:** 60-70% generic; 30-40% cetacean-specific

---

## What's Already Generic

✓ **Age structure** — any `max_age`, any survival curve (Siler, Gompertz, etc.)
✓ **Demography** — arbitrary sex ratio, infertility, maturity ogives
✓ **Litter size** — handles single and multiple offspring
✓ **Mating systems** — random mating or strong-bull modes; supports multi-mates
✓ **Social structure** — flexible pod/superpod groupings, stickiness parameters
✓ **Density dependence** — Pella-Tomlinson model works for any species
✓ **Sampling** — hierarchical (trips/sets) with arbitrary superpod pools

---

## What's Cetacean-Specific

✗ **Lactational suppression** — Hard-coded as psi_nurse vs psi_rest
✗ **Postpartum anestrus** — Fixed 3-state breeding cycle (PREGNANT / WITH CALF / RESTING)
✗ **Calf dependency** — Mother's breeding rate depends on individual calf survival
✗ **Single offspring** — Implicit assumption that dependencies are one-per-mother

---

## Future Modifications by Species

### Tier 1: Minimal Changes (Ungulates, Other Cetaceans)

**Examples:** Caribou, elk, seals, sea lions, manatees

**Changes needed:**
- Adjust `weaning_age` (e.g., 0.5 years for fast-growing ungulates)
- Adjust `psi_nurse` and `psi_rest` (weaker suppression in some ungulates)
- Adjust `litter_size` (twins common in ungulates)
- Possibly add `litter_size_variance` parameter

**Code changes:** 0-5 lines (mostly parameter adjustments)

**Effort:** 1-2 hours

---

### Tier 2: Moderate Changes (Primates, Some Terrestrial Mammals)

**Examples:** Primates with weak postpartum anestrus; tree shrews

**Changes needed:**
- Add optional parameter `fixed_s2_duration` to exit S2 regardless of calf fate
- Modify S2 transition logic (lines 493-499 in simulate_population.R):
  ```
  # Current: exit if (calf_dead OR k >= wa_breed)
  # New:     exit if (fixed_s2_duration AND k >= wa_breed) OR (calf_dead)
  ```
- Make calf_alive check conditional on `use_lactation` parameter

**Code changes:** 10-30 lines in breeding transitions

**Effort:** 3-5 hours

**Implementation sketch:**
```r
# In create.stable.pop() signature:
fixed_s2_duration = FALSE,  # if TRUE, S2 exits only at weaning_age

# In simulate_population.R, S2 transitions (line 496):
if (fixed_s2_duration) {
  # Exit only if weaned, regardless of conception or calf fate
  new_state <- ifelse(k < wa_breed, 2L, 3L)
} else {
  # Current logic: exit if calf dead, time to wean, or conceive
  stay_s2 <- !conceive & calf_alive & k < wa_breed
  ...
}
```

---

### Tier 3: Major Refactor (Rodents, Polyestrous Species)

**Examples:** Laboratory rodents, most rodents, some primates

**Changes needed:**
- Remove S2 entirely (use only S1, S3)
- Remove `calf_id` tracking
- Simplify `breeding_stationary()` to 2-state chain
- Set `psi_nurse = psi_rest` (no suppression)
- Remove calf survival dependency

**Code changes:** 50-100 lines

**Files affected:**
- `R/create_stable_pop.R` — remove S2 from Leslie matrix fecundity
- `R/simulate_population.R` — remove S2 transition block (lines 466-513)
- Remove `s2_year`, `calf_id` columns from population data.table

**Effort:** 8-12 hours

**Decision point:** Add a `breeding_mode` parameter?
```r
breeding_mode = "cetacean",  # "cetacean", "simple", or "custom"
```

---

### Tier 4: Architectural Changes (Birds, Reptiles, Fish)

**Examples:** Birds (biparental care, clutches), lizards, fish (broadcast spawning)

**Changes needed:**
- Rearchitect breeding cycle (not just parameter tweaks)
- Pair-bonding or mating group structures
- Biparental/multiparental care (multiple individuals affect calf survival)
- Clutch size instead of single birth
- Sex-specific or age-specific reproductive roles

**Code changes:** 200+ lines or full module redesign

**Effort:** 40+ hours

**Not recommended without significant refactoring.**

---

## Priority Implementation Order

If expanding to multiple taxa:

1. **✓ DONE:** Cetaceans (dolphins, whales)
2. **Next:** Ungulates (parameter-only, Tier 1)
3. **Then:** Pinnipeds (Tier 1-2)
4. **After:** Primates (Tier 2)
5. **Later:** Rodents (Tier 3, requires decision on architecture)
6. **Last:** Birds/reptiles (Tier 4, if at all)

---

## Decision Points for Future Work

### Q1: Should breeding_mode be a parameter or a version?

**Option A (Parameter):** Add `breeding_mode = "cetacean" | "simple" | "custom"` argument
- Pro: Single function handles multiple species
- Con: Increases complexity, harder to test
- Con: "custom" mode still requires code modification

**Option B (Separate functions):** Create `simulate.pop.cetacean()`, `simulate.pop.polyestrous()`, etc.
- Pro: Each is clean and focused
- Con: Code duplication, maintenance burden

**Option C (S3 methods):** Use S3 dispatch based on class
- Pro: Extensible, clean architecture
- Con: Overkill for current scope

**Recommendation for dolphins:** Stick with current monolithic version. If you expand beyond cetaceans in a serious way, revisit this.

---

### Q2: How to handle flexible litter sizes?

**Current:** Fixed `litter_size` parameter used as rpois lambda

**Future:** Some species have highly variable litter size
```r
litter_size = 1,           # mean
litter_size_dist = "pois"  # "pois", "norm", "negbinom"
```

**Effort:** 5-10 lines

---

### Q3: Should offspring survival be independent of mother?

**Current:** Calf survival drawn from fixed `surv_vec`, independent of mother's age/quality

**Future:** Age-dependent offspring survival or maternal quality effects
```r
calf_surv_maternal_age_effect = FALSE  # if TRUE, calf survival ~ mother age
```

**Effort:** 20-30 lines

---

## Testing Strategy for New Species

When you implement Tier 1-2 changes:

1. Run with known parameters (e.g., real ungulate data if available)
2. Compare Leslie matrix eigenvalue against manual calculation
3. Check that 400-year equilibrium run stays flat at K
4. Verify recruitment rates match known literature values
5. Add species-specific test cases to `tests/testthat/`

---

## Memory Notes

- **Lactation block:** Lines 466-513 in simulate_population.R (S2 transitions)
- **Breeding cycle definition:** Lines 30-32 in simulate_population.R docstring
- **Calf tracking:** calf_id column in population data.table
- **S2 duration:** Controlled by `weaning_age` parameter
- **Multi-calf handling:** Already supported (see line 679 comment)

---

**Last updated:** September 2, 2026
**Status:** Dolphins only; generalization roadmap documented for future expansion
