# This script for now contains all helper functions that are called in simulate.pop
# After all helper functions have been added here, then I will decide how to split them into different scripts as necessary.

#---------------------------------- Create.YOY ----------------------------------
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

  if(!is.null(stickiness)){

    # Add indv_name, age, sex, repro_cycle, etc. to dataframe of mating instances
    YOY_df <- mating_events %>% dplyr::select(population, pod, mother, father)

    } else{

      YOY_df <- mating_events %>% dplyr::select(population, mother, father)

      }

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

  YOY_df <- assign.growth(
    df = YOY_df,
    growth_params = growth_params
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
  )

  # Account for pod structure
  if(!is.null(stickiness)){

    YOY_df <- YOY_df %>%
      dplyr::select(indv_name, sex, age, birth_year, population, pod, repro_cycle, fertile, mother, father, L_inf, K, length)

  } else {

    YOY_df <- YOY_df %>%
      dplyr::select(indv_name, sex, age, birth_year, population, repro_cycle, fertile, mother, father, L_inf, K, length)

    }

  return(YOY_df)
}

#---------------------------------- Assign growth at birth ----------------------------------
assign.growth <- function(df, growth_params) {

  df <- df %>%
    dplyr::left_join(
      growth_params,
      by = c("population", "sex")
    )

  # MVN draw per individual
  draws <- purrr::pmap_dfr(
    list(df$L_inf, df$L_inf_sd,
         df$K, df$K_sd,
         df$rho),
    function(mu_L, sd_L, mu_K, sd_K, rho) {

      Sigma <- matrix(
        c(sd_L^2,
          rho * sd_L * sd_K,
          rho * sd_L * sd_K,
          sd_K^2),
        nrow = 2, byrow = TRUE
      )

      x <- MASS::mvrnorm(1, mu = c(mu_L, mu_K), Sigma = Sigma)

      tibble::tibble(
        L_inf = x[1],
        K     = x[2]
      )
    }
  )


  df <- df %>%
#    dplyr::select(indv_name, population, mother, father, age, birth_year, sex, min_L0, max_L0) %>%
    dplyr::select(-c(L_inf, L_inf_sd, K, K_sd, t0, rho)) %>%
    dplyr::bind_cols(draws) %>%
    dplyr::mutate(
      length = runif(n(), min_L0, max_L0)
      ) %>%
    dplyr::select(
      -c(min_L0, max_L0)
    )

  return(df)
}


#---------------------------------- Update length ----------------------------------
update.length.vb <- function(length, L_inf, K) {
  L_inf - (L_inf - length) * exp(-K)
}
