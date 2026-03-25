# Add a function to create age_length_df
# Ultimately, probably want to add this to the create_input_data.R function
# Want the output from create_input_data to include a list element that corresponds to simulated lengths/ages


# Function to sample fathers
# Function to randomly sample fathers from the same population as each mother when mating occurs
sample_fathers <- function(population, fathers_by_pop, infertility = NULL) {
  sample(fathers_by_pop[[population]], size = 1)
}

# Function to create offspring from each mating event where individuals only mate within their population
# CONFIRMED that when popstructure == "structured", mating only occurs within each population - March 25, 2026
# Optimized run time - March 25, 2026
create.YOY <- function(mothers2, fathers_by_population, ff, year){

  # UPDATED MARCH 17, 2026 to use a zero-truncated Poisson distribution (assist from Copilot) so that no mating event is assigned 0 offspring.
  # CONFIRMED that males are not being assigned a mate and then failing to breed - March 25, 2026.

  if(popstructure == "structured"){
  # Create a list where each set of potential fathers are in different list elements corresponding to their population
  fathers_by_population <- fathers %>%
    group_by(population) %>%
    summarise(indv_name = list(indv_name)) %>%
    deframe()

  # Each row of mothers2 corresponds to a mating event with a different father. This function first assigns a father to each row of mothers2; then, each mating event is assigned a number of offspring based on the object ff; finally, the dataframe is expanded with uncount() so that each row corresponds to a single offspring between the mother and father.
  mating_events <- mothers2 %>%
    mutate(
      father = map_chr(population, ~ sample_fathers(.x, fathers_by_population)),
      # zero-truncated Poisson:
      num_off = {
        p0 <- exp(-ff)
        u  <- runif(n())
        qpois(p0 + u * (1 - p0), ff)
      }
    ) %>%
    group_by(mother, father) %>%
    tidyr::uncount(num_off) %>%   # repeats rows 'num_off' times
    ungroup()

} else if(popstructure == "panmictic"){
  mating_events <- mothers2 %>%
    mutate(
      father = map_chr(1:n(), ~sample(fathers$indv_name, size = 1, replace = TRUE)),
      # zero-truncated Poisson:
      num_off = {
        p0 <- exp(-ff)
        u  <- runif(n())
        qpois(p0 + u * (1 - p0), ff)
      }
    ) %>%
    group_by(mother, father) %>%
    tidyr::uncount(num_off) %>%   # repeats rows 'num_off' times
    ungroup()
}

  # How many instances of mating?
  n <- nrow(mating_events)

  # Add indv_name, age, sex, repro_cycle, etc. to dataframe of mating instances
  YOY_df <- mating_events
  YOY_df$indv_name <- replicate(
    n,
    paste(sample(letters, 20, replace = TRUE), collapse = "")
  )

  YOY_df$age <- 0L
  YOY_df$birth_year <- year

  YOY_df$sex <- sample(
    c("F", "M"),
    size = n,
    replace = TRUE,
    prob = c(female_fraction, 1 - female_fraction)
  )

  YOY_df <- YOY_df %>% mutate(
    repro_cycle = case_when(
      sex == "F" ~ sample(
        c(1:mating_periodicity),
        size = n(),
        replace = T
      ),
      TRUE ~ NA),
    fertile = rbinom(n(), size = 1, prob = 1 - infertility) == 1
  ) %>%
    dplyr::select(indv_name, sex, age, birth_year, population, repro_cycle, mother, father)

  return(YOY_df)
}
