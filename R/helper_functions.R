# Add a function to create age_length_df
# Ultimately, probably want to add this to the create_input_data.R function
# Want the output from create_input_data to include a list element that corresponds to simulated lengths/ages


# Function to sample fathers
# Function to randomly sample fathers from the same population as each mother when mating occurs
# sample_fathers <- function(population, fathers_by_pop) {
#   sample(fathers_by_pop[[population]], size = 1, replace = T)
# }

# Function to create offspring from each mating event where individuals only mate within their population
# CONFIRMED that when popstructure == "structured", mating only occurs within each population - March 25, 2026
# Optimized run time - March 25, 2026
create.YOY <- function(mothers2, fathers, litters, year){

  # The code below splits the litter size - stored in litters - among separate "mating events" - the number of rows in mothers2 - based on draws from a multinomial distribution. THEN a father is assigned to all mating events that produce at least one offspring. This rectifies the foundational issue we had previously where fathers could be assigned to a mating event that produces zero offspring, which would have been fine, except fathers are also randomly assigned to females. Collectively, this meant that some males were not assigned to females and others were assigned to females but produced zero offspring, resulting in two moments to lose successful male breeders. The first filter is expected, but the second was unintentional. We want to be able to control this better in the future. The process below fixes this.

  if(popstructure == "structured"){

  # Create a list where each set of potential fathers are in different list elements corresponding to their population
  fathers_by_population <- fathers %>%
    group_by(population) %>%
    summarise(indv_name = list(indv_name)) %>%
    deframe()

  # The dataframe mating_events has one row per YOY, with mother, father, and population in the correct(ly labeled) columns, so we can use this dataframe as the foundation of YOY_df.
  mating_events <- mothers2 %>%
    left_join(litters, by = c("mother", "population")) %>%
    group_by(mother) %>%
    mutate(
      litter_per_mating = rmultinom(
        n = 1,
        size = as.integer(litter_size[1]),
        prob = rep(1/n(), n())
      )[,1]
    ) %>%
    ungroup() %>%
    filter(litter_per_mating > 0) %>%
    mutate(
      father = map_chr(
        population,
        ~ sample(fathers_by_population[[.x]], size = 1, replace = TRUE)
      )
    ) %>%
    tidyr::uncount(litter_per_mating)

} else if(popstructure == "panmictic"){

  # The dataframe mating_events has one row per YOY, with mother, father, and population in the correct(ly labeled) columns, so we can use this dataframe as the foundation of YOY_df.
  # TO DO: Can I adjust the below function so that no mating events produce zero offspring?
  mating_events <- mothers2 %>%
    left_join(litters, by = c("mother", "population")) %>%
    group_by(mother) %>%
    mutate(
      litter_per_mating = rmultinom(
        n = 1,
        size = as.integer(litter_size[1]),
        prob = rep(1/n(), n())
      )[,1]
    ) %>%
    ungroup() %>%
    filter(litter_per_mating > 0) %>%
    mutate(father = sample(fathers$indv_name, size = n(), replace = TRUE)) %>%
    tidyr::uncount(litter_per_mating)

}

  # How many instances of mating?
  n <- nrow(mating_events)

  # Add indv_name, age, sex, repro_cycle, etc. to dataframe of mating instances
  YOY_df <- mating_events %>% dplyr::select(population, mother, father)

  YOY_df$indv_name <-
    indv_name <- sprintf(
      "%03d_%020d",
      year,
      seq_len(n)
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
    dplyr::select(indv_name, sex, age, birth_year, population, repro_cycle, fertile, mother, father)

  return(YOY_df)
}
