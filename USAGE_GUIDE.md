# sharkyIBM Usage Guide

## Quick Start

The sharkyIBM package provides a complete workflow to simulate dolphin populations and sample from them. Here's the three-step process:

### 1. Calibrate a Stable Population

```r
library(sharkyIBM)

# Define your life-history parameters
survival <- exp(-c(0.15, rep(0.04, 40)))  # Siler mortality curve
ogive_f <- plogis(0:40, location = 9, scale = 1/1.5)  # Maturity at age 9

# Calibrate the population (estimates density-dependence parameters)
config <- create.stable.pop(
  max_age = 40,
  survival = survival,
  maturity_age = ogive_f,
  litter_size = 1,
  pop_size = 20000,
  rho = -3.255,  # Lactational suppression
  weaning_age = 2L,
  density_dependence = TRUE,
  target_interval = 2.84,
  target_depletion = 0.3
)
```

### 2. Run a Simulation

```r
# Optionally solve for sustainable fishing
f_result <- solve_F_sustainable(config)

# Simulate for 100 years
sim <- simulate.pop(
  sim_config = config,
  num_years = 100,
  F_t = f_result$F_sustainable,
  init_depletion = 0.3
)
```

### 3. Sample the Population

```r
# Simulate a fishing operation with social structure
samples <- sample.pop(
  sim_output = sim,
  n_trips = 5,
  n_sets = 4,
  sample_size = 30,
  stickiness_set = 0.8,
  stickiness_trip = 0.6
)

# Unique individuals captured
unique_samples <- samples[!duplicated(samples$id), ]
```

## Key Concepts

### Markov Breeding Cycle
Females cycle through three reproductive states:
- **PREGNANT (S1):** Carrying a fetus
- **WITH DEPENDENT CALF (S2):** Nursing, with lactational suppression
- **RESTING (S3):** Ready to breed again

When a calf dies, the mother is released from suppression early (compensatory breeding).

### Density Dependence
The population self-corrects through conception rate feedback:
- When population is small → higher conception rates → faster recovery
- When population is large → lower conception rates → slower growth

Calibrated to match observed calving intervals at a reference depletion level.

### Individual-Based Simulation
Each dolphin is tracked with:
- Unique ID, age, sex
- Reproductive state and history
- Mother and father IDs (complete pedigree)
- Social group membership (pod, superpod)

This allows realistic sampling that mimics fishing operations and captures social structure effects.

## Parameters for Eastern Spinner Dolphins

All recommended parameters are derived from the `dolphin_population_model.qmd` specification:

```r
# Life history
max_age <- 40
M_adult <- 0.04
M_calf_excess <- 0.15
b_juv <- 0.5
a_sen_onset <- 28

# Maturity
a50_maturity <- 9
maturity_slope <- 1.5

# Breeding
rho <- -3.255  # log-odds of conception while nursing
target_interval <- 2.84  # calving interval at D=0.3
target_depletion <- 0.3
z_pt <- 2.39  # Pella-Tomlinson shape

# Social structure
pod_size <- 20
superpod_size <- 10
stickiness_year <- 0.9
```

## Output Interpretation

### Calibration (`config`)
- `theta`: Log-odds of conception (internal parameter)
- `interval_K`: Emergent calving interval at carrying capacity
- `dd_max`: Density-dependence compensation strength
- `K_1plus`: Age-1+ animals at carrying capacity

### Simulation (`sim`)
- `pop_summary`: Population counts by year, age, sex
- `snapshots`: Individual-level data at specified years
- `depletion`: Depletion trajectory (N_1+ / K_1+)

### Samples (`samples`)
- One row per captured individual
- Columns: id, age, sex, birth_year, mother_id, father_id, trip, set, year
- Duplicates indicate same individual caught in multiple sets (due to social stickiness)

## Advanced Features

### Fishing Mortality
```r
# Define fishing parameters
F_t <- rep(0.02, 100)  # Fishing rate per year
selectivity <- c(0, rep(1, 40))  # Age 0 not vulnerable

sim <- simulate.pop(config, num_years = 100, F_t = F_t, selectivity = selectivity)
```

### Depleted Initialization
Start population directly at depletion without long transient:
```r
sim <- simulate.pop(config, num_years = 100, F_t = F_result$F_sustainable, 
                    init_depletion = 0.3)
```

### Custom Sampling Scheme
```r
# Different sampling modes
samples <- sample.pop(
  sim_output = sim,
  sampling = "random",  # vs "superpod" (default)
  sample_per = "trip",  # vs "set" (default)
  stickiness_set = 0.5,
  stickiness_trip = 0.3
)
```

## For More Details

See the full vignette: `vignette("using_sharkyibm")`

Or read the function documentation:
- `?create.stable.pop()` — Calibration and density-dependence
- `?simulate.pop()` — Running the IBM
- `?sample.pop()` — Hierarchical sampling
- `?solve_F_sustainable()` — Fishing mortality solutions
