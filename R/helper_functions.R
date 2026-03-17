# Add a function to create age_length_df
# Ultimately, probably want to add this to the create_input_data.R function
# Want the output from create_input_data to include a list element that corresponds to simulated lengths/ages


# Function to sample fathers
# Function to randomly sample fathers from the same population as each mother when mating occurs
sample_fathers <- function(population, fathers_by_pop, infertility = NULL) {
  sample(fathers_by_pop[[population]], size = 1)
}

# Function to create offspring from each mating event where individuals only mate within their population
createYOY.byPop <- function(mothers2, fathers_by_population, ff, year){

  # Creates tibble of mating events (kinda) and then expands so that each row corresponds to a single offspring between the mother and father.
  # UPDATED MARCH 17, 2026 to use a zero-truncated Poisson distribution (assist from Copilot) so that no mating event is assigned 0 offspring.
  mating_events <- mothers2 %>%
    lazy_dt() %>%
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
    ungroup() %>%
    as_tibble()

  #### PICK UP HERE AFTER MARCH 17

  YOY.df <- mating.events %>%
    lazy_dt() %>%
    select(-num.off) %>%
    mutate(
      indv.name = map_chr(
        1:n(),
        ~paste(sample(letters, size = 20, replace = T),
               collapse="")
      ),
      age.x = 0,
      birth.year = year,

      repro.cycle = map_dbl(
        1:n(),
        ~sample(1:mating.periodicity, size = 1, replace = TRUE)
      ),

      sex = map_chr(
        1:n(),
        ~sample(c('F','M'),
                size = 1,
                prob = c(init.prop.female, 1-init.prop.female))
      ),

      indv_length = map_dbl(
        1:n(),
        ~runif(n = 1, min = length.at.birth[1], max = length.at.birth[2])
      ),

      beta_0 = case_when(
        population == "MX" ~ MX.beta.0,
        population == "ES" ~ ES.beta.0,
        population == "EC" ~ EC.beta.0,
        TRUE ~ NA),

      beta_1 = case_when(
        population == "MX" ~ MX.beta.1,
        population == "ES" ~ ES.beta.1,
        population == "EC" ~ EC.beta.1,
        TRUE ~ NA)
    ) %>%
    left_join(age.length.df, by = "age.x") %>%
    as_tibble()

  return(YOY.df)
}

# Create offspring from each mating event where individuals only mate within their population
createYOY.panmictic <- function(mothers, fathers, ff, year){

  mating.events <- mothers %>%
    lazy_dt() %>%
    mutate(father.x = map_chr(1:n(), ~sample(fathers$indv.name, size = 1, replace = TRUE))) %>%
    ungroup() %>%
    mutate(num.off = rpois(n(), ff)) %>%
    group_by(mother.x, father.x) %>%
    slice(rep(1:n(), num.off)) %>%
    ungroup() %>%
    as_tibble()

  YOY.df <- mating.events %>%
    lazy_dt() %>%
    select(-num.off) %>%
    mutate(
      indv.name = map_chr(
        1:n(),
        ~paste(sample(letters, size = 20, replace = T),
               collapse="")
      ),
      age.x = 0,
      birth.year = year,

      repro.cycle = map_dbl(
        1:n(),
        ~sample(1:mating.periodicity, size = 1, replace = TRUE)
      ),

      sex = map_chr(
        1:n(),
        ~sample(c('F','M'),
                size = 1,
                prob = c(init.prop.female, 1-init.prop.female))
      ),

      indv_length = map_dbl(
        1:n(),
        ~runif(n = 1, min = length.at.birth[1], max = length.at.birth[2])
      ),

      beta_0 = case_when(
        population == "MX" ~ MX.beta.0,
        population == "ES" ~ ES.beta.0,
        population == "EC" ~ EC.beta.0,
        TRUE ~ NA),

      beta_1 = case_when(
        population == "MX" ~ MX.beta.1,
        population == "ES" ~ ES.beta.1,
        population == "EC" ~ EC.beta.1,
        TRUE ~ NA)
    ) %>%
    left_join(age.length.df, by = "age.x") %>%
    as_tibble()

  return(YOY.df)
}

