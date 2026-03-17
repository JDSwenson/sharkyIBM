#' Simulate age-structured shark populations and sampling
#'
#' Runs a forward-time, age-structured simulation across one or more populations,
#' including reproduction, growth, and survival, with optional sampling of
#' individuals in specified years. Returns yearly population metrics and sampled
#' individuals suitable for downstream analysis.
#'
#' @param init_pop_size Integer vector. Initial population size per population (order should match columns of `Nages`).
#' @param init_prop_female Numeric scalar in \[0, 1\]. Initial proportion female in the population.
#' @param Nages Integer matrix. Rows are ages, columns are populations; entries give counts used to build initial ages.
#' @param mating_periodicity Integer (>= 1). Female reproductive cycle periodicity (e.g., 1 = annual, 2 = biennial).
#' @param repro_age Integer (>= 0). Minimum female age for reproduction (knife-edged maturity).
#' @param YOY_survival Numeric in \[0, 1\]. Survival parameter for young-of-year (if not overridden by a survival table).
#' @param juvenile_survival Numeric in \[0, 1\]. Survival parameter for juveniles (if not overridden).
#' @param adult_survival Numeric in \[0, 1\]. Survival parameter for adults (if not overridden).
#' @param max_age Integer (>= 0). Individuals at or above this age are removed (die of senescence).
#' @param num_mates Integer vector or scalar. Number of mates per mother (sampled with replacement).
#' @param ff List or numeric. Fecundity-related parameters passed to offspring-generation helpers.
#' @param burn_in Integer (>= 0). Number of initial years simulated before analysis years.
#' @param num_years Integer (>= 1). Number of analysis years after burn-in.
#' @param age_length_df Data frame/tibble. Must contain at least `age.x`, `mean_length`, and `age_length_sd`.
#' @param movement_array Three-dimensional array of age-based movement probabilities to locations where individual animals will be sampled. Dimensions should be \[age, sampling_location, population\]. In other words, each population should have its own matrix, ordered the same as the populations. The array should include age 0 individuals, so the number of rows should be max_age+1, and the movement probability for age a should be in the row a+1. For example, the probability to find age 6 individuals from population 1 at sampling location 4 would be specified in the movement_array \[7, 4, 1\].
#' @param infertility Numeric in \[0, 1\] specifying the proportion of the population that is infertile throughout their lives. This characteristic is assigned at birth and stays with an individual throughout its life.
#'
#' @details
#' This function orchestrates initialization, breeding, growth, survival, and optional
#' sampling of individuals. It relies on several helper functions (e.g., for offspring
#' generation, reproduction probabilities, dispersal, and sampling) that are defined
#' elsewhere in the package.
#'
#' The return value is a list with population metrics per year and any sampled
#' individuals.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{pop.size}: data frame/tibble of population metrics by year and population.
#'   \item \code{samples.df}: data frame/tibble of sampled individuals (may be empty if no sampling occurs).
#' }
#' The list is returned invisibly.
#'
#' @examples
#' \dontrun{
#' # Minimal scaffold showing inputs (uses tiny sizes and fake tables)
#' init_pop_size <- c(MX = 10, ES = 8, EC = 12)
#' Nages <- matrix(c(3,4,3,
#'                   4,2,4,
#'                   3,2,5), nrow = 3, byrow = TRUE)
#' colnames(Nages) <- names(init_pop_size)
#'
#' age_length_df <- tibble::tibble(
#'   age = 0:10,
#'   mean_length = seq(80, 130, length.out = 11),
#'   age_length_sd = rep(8, 11)
#' )
#'
#' out <- simulate.pop(
#'   init_pop_size = init_pop_size,
#'   init_prop_female = 0.5,
#'   Nages = Nages,
#'   mating_periodicity = 1,
#'   repro_age = 5,
#'   YOY_survival = 0.7,
#'   juvenile_survival = 0.85,
#'   adult_survival = 0.92,
#'   max_age = 20,
#'   num_mates = 1:2,
#'   ff = list(),
#'   burn_in = 0,
#'   num_years = 1,
#'   age_length_df = age_length_df
#' )
#' }
#'
#' @seealso
#' Helper functions you will define and document separately, e.g.,
#' offspring generation and sampling utilities.
#'
#' @export

simulate.pop <- function(input_data,
                         mating_periodicity,
                         f_maturity_age,
                         m_maturity_age,
                         num_mates,
                         num_years,
                         female_fraction = 0.5,
                         age_length_df = NULL,
                         movement_array = NULL,
                         infertility = NULL
                         ) {

  # Save initial values as distinct R objects
  init_pop_size <- input_data$numbers_at_age
  max_age <- max(input_data$numbers_at_age$age)
  f <- max(input_data$fecundity)
  YOY_survival <- input_data$s0
  juvenile_survival <- input_data$survival[2:(maturity_age-1)]
  adult_survival <- input_data$survival[maturity_age:(max_age-1)]
  ff <- f/female_fraction * mating_periodicity/mean(num_mates) # Female fecundity per breeding event at equilibrium


  # Make initial
  init_ages <- rep(init_pop_size$age, times = init_pop_size$N)
  init_pops <- rep(init_pop_size$population, times = init_pop_size$N)
  init_sex <- sample(
    c("F", "M"),
    size = sum(init_pop_size$N),
    prob = c(female_fraction, 1 - female_fraction),
    replace = T)
  init_repro_cycle <- sample(
    c(1:mating_periodicity),
    size = sum(init_pop_size$N),
    replace = T)

  # Summarize population numbers
  total_pop_sizes_df <- init_pop_size %>% group_by(population) %>%
    reframe(total = sum(N, na.rm = T)) %>%
    arrange(population)

  total_pop_sizes_vec <- rep(total_pop_sizes_df$population, times = total_pop_sizes_df$total)

  ###############################################`
  ####---------Set up initial population-----####
  ###############################################`
  # Initial population
  init_pop <- tibble(
    # Assign a random 20 character string for each individual's name
    indv_name = map_chr(
      1:sum(init_pop_size$N), # sum because we have multiple populations
      ~paste(sample(letters, size = 20, replace = T),
             collapse="")
    ),
    birth_year = -1, # This is a place holder for individuals born within the simulation
    age = init_ages, # Assign ages based on the stable age distribution
    mother = "xxxxx", # The individuals in the initial population do not have known mothers
    father = "xxxxx", # The individuals in the initial population do not have known fathers
    sex = init_sex,
    population = total_pop_sizes_vec) %>%
    mutate(repro_cycle = case_when(
      sex == "F" ~ init_repro_cycle[row_number()],
      TRUE ~ NA),
      fertile = rbinom(n(), size = 1, prob = 1 - infertility) == 1
      )

  # TO DO: Add code to quickly and easily simulate initial starting lengths for all individuals from age_length_df (if supplied). Might want to just add to the tibble above, if can simulate in create_input_data.R script.
# if(!is.null(age_length_df)){
#   # Join with age.length table, assign age, and repro probability
#   init_pop2 <- init_pop %>%
#     lazy_dt() %>%
#     left_join(age_length_df, by = "age") %>%
#     mutate(indv_length = rtruncnorm(n(), mean = mean_length, sd = age_length_sd, a = 0.2)) %>% #Assign individual length -- make sure nobody grows backwards, so set lower limit of 0.2
#     mutate(beta_0 = case_when(
#       population == "MX" ~ MX.beta.0,
#       population == "ES" ~ ES.beta.0,
#       population == "EC" ~ EC.beta.0,
#       TRUE ~ NA),
#       beta_1 = case_when(
#         population == "MX" ~ MX.beta.1,
#         population == "ES" ~ ES.beta.1,
#         population == "EC" ~ EC.beta.1,
#         TRUE ~ NA)) %>% # Save values for growth curve so we can vectorize with case_when
#     as_tibble() %>%
#     mutate(repro_prob = case_when( # Store probability of reproduction
#       age < 5 ~ 0, # No individuals younger than age 5 will reproduce (5 is an arbitrary number)
#       age >= 5 ~ repro.prob(beta.0 = beta_0, beta.1 = beta_1, TLflex = indv_length),
#       TRUE ~ NA))
# }

  ####----------Breeding----------####
  repro_cycle_vec <- rep(1:mating_periodicity, times = num_years) # Generate a vector which will be used to determine if it is an even or odd breeding year (or a 1/3 breeding year)

  ####--------- For year 0 breeding
  #------------Mothers------------#
  # Mothers with knife-edged maturity to allow the population to stay stable (at least for now)
  if(!is.null(f_maturity_age)){
  mothers <- init_pop %>% filter(sex == 'F',
                                 age >= f_maturity_age,
                                 fertile, # filters to only keep indvs with fertile == T, but is faster without the conditional statement
                                 repro_cycle == repro_cycle_vec[1]) # Determine which females are available to breed in this year
  }

  # TO DO: create vector of mothers from age-specific fecundity vector
  mothers <- mothers %>% mutate(num_mates = sample(num_mates, size = n(), replace = TRUE)) # Assign random number of mates to each mother

  # END HERE March 9, 2026

  # Make a new dataframe where each row corresponds to an instance of mating
  mothers2 <- mothers %>%
    lazy_dt() %>%
    group_by(indv_name) %>%
    slice(rep(1:n(), num_mates)) %>%
    ungroup() %>%
    select(indv_name, population) %>%
    rename(mother = indv_name) %>%
    as_tibble()

  #------------Fathers------------#
  fathers <- init_pop %>% filter(sex=='M',
                                 init_pop$age >= m_maturity_age
                                 ) %>% # Uncomment for age-based maturity
    select(indv_name, population)

  if(popstructure == "structured"){

    # Create a list where each set of potential fathers are in different list elements corresponding to their population
    # Confirmed that this works
    fathers_by_population <- fathers %>%
      group_by(population) %>%
      summarise(indv_name = list(indv_name)) %>%
      deframe()

    # PICK UP WITH createYOY.byPop after March 17, 2026

    # Create dataframe of mating events and generate initial offspring from each mating event
    YOY_df <- createYOY.byPop(mothers2, fathers_by_population, ff, year = 0)

  } else if(popstructure == "panmictic"){

    # Create dataframe of mating events and generate initial offspring from each mating event
    YOY.df <- createYOY.init.panmictic(mothers2, fathers, ff)

  }

  # This dataframe holds the population at the end of the first year of the simulation
  year.end.pop.0 <- bind_rows(init.pop2, YOY.df)

  # And finally, we assign age-specific mortality rates to each individual and then determine whether they survive into the next year or not
  loopy.pop <- year.end.pop.0 %>% left_join(survival.df, by = "age.x") %>%
    rename(survival_rate = survival.rate) %>%
    mutate(
      survival = case_when(
        runif(nrow(year.end.pop.0)) <= survival_rate ~ "S",
        .default = "M"
      )
    )


  # At the end of year 0 ...
  print(paste("year 0 ", names(table(loopy.pop$population)), ": N_mothers=", table(mothers2$population), ", N_pups=", table(YOY.df$population), ", N_deaths=", table(loopy.pop$population[loopy.pop$survival=="M"]), ", Total N=", table(loopy.pop$population[loopy.pop$survival=="S"]) , sep=""))


  #############################################################`
  ####---------Loop through remaining simulation years-----####
  #############################################################`

  ####---------------####
  pop.size <- data.frame() # Initialize dataframe for storing population size

  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year -- to save space, we won't populate this for now
  plot.list <- list() # If we want to plot number of mature/breeding females/males each year for each iteration

  samples.df <- NULL

  parents.tibble <- tibble() # This will store info on offspring distribution per parent
  moms.temp = dads.temp <- NULL

  for(v in 1:(burn.in + Num.years)){ # Loop through all of the years in the simulation - the burn in and the years that matter

    # Bring in the data from the previous iteration, but only include those that survive (and leave out columns that need to be updated)
    data1 <- loopy.pop %>% dplyr::filter(survival == "S") %>%
      dplyr::select(indv.name,
                    birth.year,
                    age.x,
                    mother.x,
                    father.x,
                    sex,
                    repro.cycle,
                    population,
                    indv_length,
                    mean_growth_rate,
                    growth_rate_sd,
                    beta_0,
                    beta_1) %>%
      mutate(age.x = age.x+1, # Increase age by one year - happy birthday survivors!
             indv_length = indv_length + rtruncnorm(n = n(), mean = mean_growth_rate, sd = growth_rate_sd, a = 0.2) # Increase length based on age-specific growth rate
      ) %>%
      dplyr::select(-c(mean_growth_rate, growth_rate_sd)) %>% # Now remove growth rate so we can assign the proper age/length-specific growth rate (after advancing the age and length)
      left_join(age.length.df, by = "age.x") %>% # Assign appropriate growth rate
      dplyr::select(-c(mean_length, age_length_sd)) %>%
      mutate(repro_prob = case_when( # Male probability of breeding is based on length
        age.x < 5 ~ 0, # No individuals younger than age 5 will reproduce (5 is an arbitrary number)
        age.x >= 5 ~ repro.prob(beta.0 = beta_0, beta.1 = beta_1, TLflex = indv_length),
        .default = NA)) # Store probability of reproduction

    #If individuals are older than max.age, they will be killed

    ####----------Breeding----------####
    #------------Mothers------------#
    mothers <- data1 %>% filter(sex=='F',
                                age.x>=repro.age,
                                repro.cycle == repro.cycle.vec[v+1]) # Determine which females are available to breed in this year

    # Add column that contains the number of mates each mother will mate with this year
    mothers <- mothers %>% mutate(num.mates = sample(num.mates, size = nrow(mothers), replace = TRUE))

    # Make a new row where each row corresponds to an instance of mating
    mothers2 <- mothers %>%
      lazy_dt() %>%
      group_by(indv.name) %>%
      slice(rep(1:n(), num.mates)) %>%
      ungroup() %>%
      select(indv.name, population) %>%
      rename(mother.x = indv.name) %>%
      as_tibble()

    #------------Fathers------------#
    fathers <- data1 %>% dplyr::filter(sex=='M',
                                       #init.pop$age.x>=repro.age #Uncomment for age-based maturity
                                       runif(n()) <= repro_prob
    ) %>% # Determine which fathers are available to breed in this year
      dplyr::select(indv.name, population)

    if(popstructure == "structured"){

      # Create a list where each set of potential fathers are in different list elements corresponding to their population
      # Confirmed that this works
      fathers_by_population <- fathers %>%
        group_by(population) %>%
        summarise(indv.name = list(indv.name)) %>%
        deframe()

      # Create dataframe of mating events, and then assign offspring to each mating event
      YOY.df <- createYOY.byPop(mothers2, fathers_by_population, ff, v)

    } else if(popstructure == "panmictic"){

      # Create dataframe of mating events, and then assign offspring to each mating event
      YOY.df <- createYOY.panmictic(mothers2, fathers, ff, v)

    }
    #Only bother assigning a sampling location for the years we're taking samples; otherwise just slows down code.
    if(v >= min(sample.years)){

      # YOY from the same mother are assigned to the same sampling location
      YOY.df <- YOY.df %>% group_by(mother.x, population) %>%
        mutate(
          sampling_location = sample(
            sampling.locations,
            1,
            prob = c(dispersal_kernel(age = 0, birth_population = population[1])),
            replace = TRUE)
        ) %>%
        ungroup()

      #Pull out mothers and assign them the same sampling location as their offspring from this year
      mother.sample.df <- YOY.df %>% select(indv.name = mother.x, sampling_location) %>% distinct()

      mothers.df <- mother.sample.df %>% lazy_dt() %>%
        left_join(mothers, by = "indv.name") %>%
        select(-num.mates) %>%
        as_tibble()

      #Assign all other individuals assigned randomly
      loopy.pop <- data1 %>%
        lazy_dt() %>%
        filter(!indv.name %chin% YOY.df$mother.x) %>%
        mutate(
          sampling_location = map2_chr(
            age.x,
            population,
            ~sample(
              sampling.locations,
              1,
              prob = c(dispersal_kernel(.x, .y)),
              replace = TRUE
            ))
        ) %>%
        as_tibble() %>%
        bind_rows(YOY.df, mothers.df)


    } else {

      # No need to assign sampling location if we're not sampling this year
      loopy.pop <- bind_rows(data1, YOY.df)

    }

    # Assign survival or mortality based on age-specific survival probabilities
    loopy.pop <- loopy.pop %>% left_join(survival.df, by = "age.x") %>%
      select(-c(stable.age)) %>%
      rename(survival_rate = survival.rate) %>%
      mutate(
        survival = case_when(
          age.x >= max.age ~ "M",
          runif(nrow(loopy.pop)) <= survival_rate ~ "S",
          .default = "M"
        )
      )

    ###############################################`
    ####---------------Sampling----------------####
    ###############################################`
    if(v %in% sample.years){

      samples.df.temp <- loopy.pop %>%
        group_by(sampling_location) %>%
        group_map(~sample_fixed(.x, samples.vec[.y$sampling_location[1]]), .keep = TRUE) %>%
        bind_rows() %>%
        mutate(capture.year = v,
               iteration = iter)

      samples.df <- bind_rows(samples.df, samples.df.temp)

    }


    ###############################################`
    ####-------------Save metrics--------------####
    ###############################################`

    # Calculate number of produced offspring per mother and father this year
    moms.temp <- YOY.df %>% group_by(mother.x, population) %>%
      summarize(num.off = n()) %>%
      rename(parent = mother.x) %>%
      mutate(year = v, parent.sex = "mother")

    dads.temp <- YOY.df %>% group_by(father.x, population) %>%
      summarize(num.off = n()) %>%
      rename(parent = father.x) %>%
      mutate(year = v, parent.sex = "father")

    mothers3 <- mothers2 %>% group_by(mother.x, population) %>%
      summarize(num.off = n()) %>%
      rename(parent = mother.x)

    # Add to the tibble of offspring distribution - can use to check if/whether some indvs are reproducing much more
    #   parents.tibble <- rbind(parents.tibble, moms.temp, dads.temp)

    # Print info about the population to the console
    cat(paste("\nyear", v, " ", names(table(loopy.pop$population)),
              "N_mothers=", table(moms.temp$population),
              "N_fathers=", table(dads.temp$population),
              "\nN_pups=", table(YOY.df$population),
              "\nN_deaths=", table(loopy.pop$population[loopy.pop$survival=="M"]),
              "\nTotal N= ", table(loopy.pop$population[loopy.pop$survival=="S"]) , sep=" "))


    # Save the population size
    pop.size.vec.MX <- cbind.data.frame(year=v,
                                        population = "MX",
                                        population_size=nrow(data1[data1$population == "MX",]),
                                        Male.adult.pop = nrow(data1[data1$sex == "M" & data1$age.x >= repro.age & data1$population == "MX",]),
                                        Female.adult.pop = nrow(data1[data1$sex == "F" & data1$age.x >= repro.age & data1$population == "MX",]),
                                        Num.mothers = nrow(moms.temp[moms.temp$population == "MX",]),
                                        Num.fathers = nrow(dads.temp[dads.temp$population == "MX",]),
                                        iteration = iter)

    pop.size.vec.ES <- cbind.data.frame(year=v,
                                        population = "ES",
                                        population_size=nrow(data1[data1$population == "ES",]),
                                        Male.adult.pop = nrow(data1[data1$sex == "M" & data1$age.x >= repro.age & data1$population == "ES",]),
                                        Female.adult.pop = nrow(data1[data1$sex == "F" & data1$age.x >= repro.age & data1$population == "ES",]),
                                        Num.mothers = nrow(moms.temp[moms.temp$population == "ES",]),
                                        Num.fathers = nrow(dads.temp[dads.temp$population == "ES",]),
                                        iteration = iter)

    pop.size.vec.EC <- cbind.data.frame(year=v,
                                        population = "EC",
                                        population_size=nrow(data1[data1$population == "EC",]),
                                        Male.adult.pop = nrow(data1[data1$sex == "M" & data1$age.x >= repro.age & data1$population == "EC",]),
                                        Female.adult.pop = nrow(data1[data1$sex == "F" & data1$age.x >= repro.age & data1$population == "EC",]),
                                        Num.mothers = nrow(moms.temp[moms.temp$population == "EC",]),
                                        Num.fathers = nrow(dads.temp[dads.temp$population == "EC",]),
                                        iteration = iter)

    pop.size <- bind_rows(pop.size, pop.size.vec.MX, pop.size.vec.ES, pop.size.vec.EC)

    # For checking that male maturity isn't changing over the simulation (this was a bug earlier) ...
    # plot.list[[v]] <- data1 %>% dplyr::filter(sex=='M') %>%
    #   mutate(unif = runif(n()),
    #          repro = case_when(
    #            unif < repro_prob ~ "yes",
    #            .default = "no"
    #          )) %>%
    #   gghistogram(x = c("repro_prob")) +
    #   ggtitle(paste0("Year ", v))

    # When I figure out how to efficiently store this info for a large population, then I will uncomment this.
    #loopy.list[[v]] <- loopy.pop



  } # End loop over sim years

  # Label the list elements with the year
  # names(loopy.list) <- paste0("year.end.pop.", seq(1:(burn.in + Num.years)), "_iteration_", iter)

  return(invisible(list(pop.size, samples.df)))
}
