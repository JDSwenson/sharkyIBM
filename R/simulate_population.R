simulate.pop <- function(init_pop_size,
                         init_prop_female,
                         Nages,
                         mating_periodicity,
                         repro_age,
                         YOY_survival,
                         juvenile_survival,
                         adult_survival,
                         max_age,
                         num_mates,
                         ff,
                         burn_in,
                         num_years,
                         age_length_df
                         ) {

  ###############################################`
  ####---------Set up initial population-----####
  ###############################################`
  init_ages <- NULL

  # init_ages is a vector with the different ages of the individuals for all populations in a single vector.
  for(z in 1:ncol(Nages)){
    for(y in 1:nrow(Nages)){
      init.ages <- c(init.ages, rep(y, Nages[y, z]))
    }
  }

  # Initial population
  init.pop <- tibble(

    # Assign a random 20 character string for each individual's name
    indv.name = map_chr(
      1:sum(init.pop.size), # sum because we have multiple populations
      ~paste(sample(letters, size = 20, replace = T),
             collapse="")
    ),

    birth.year = -1, # This is a place holder for individuals born within the simulation
    age.x = init.ages, # Assign ages based on the stable age distribution
    mother.x = "xxxxx", # The individuals in the initial population do not have known mothers
    father.x = "xxxxx", # The individuals in the initial population do not have known fathers

    sex = map_chr( # Randomly draw sex based on the proportions set in the parameter section
      1:sum(init.pop.size),
      ~sample(c('F','M'),
              size = 1,
              prob = c(init.prop.female, 1-init.prop.female))
    ),

    repro.cycle = map_dbl(
      1:sum(init.pop.size),
      ~sample(1:mating.periodicity, size = 1) # Randomly assign whether this individual mother will breed in even or odd years (relevant for multiennial breeding only)
    )
  )


  # Assign populations based on row numbers
  total.pop.sizes <- cumsum(init.pop.size)
  row.numbers <- seq(1:sum(init.pop.size))

  init.population <- sapply(row.numbers, function(row_num) {
    for (i in seq_along(total.pop.sizes)) {
      if (row_num <= total.pop.sizes[i]) {
        return(populations[i])
      }
    }
  })

  init.pop2 <- init.pop %>%
    mutate(population = init.population)


  # Join with age.length table, assign age, and repro probability
  init.pop2 <- init.pop2 %>%
    lazy_dt() %>%
    left_join(age.length.df, by = "age.x") %>%
    mutate(indv_length = rtruncnorm(n(), mean = mean_length, sd = age_length_sd, a = 0.2)) %>% #Assign individual length -- make sure nobody grows backwards, so set lower limit of 0.2
    mutate(beta_0 = case_when(
      population == "MX" ~ MX.beta.0,
      population == "ES" ~ ES.beta.0,
      population == "EC" ~ EC.beta.0,
      TRUE ~ NA),
      beta_1 = case_when(
        population == "MX" ~ MX.beta.1,
        population == "ES" ~ ES.beta.1,
        population == "EC" ~ EC.beta.1,
        TRUE ~ NA)) %>% # Save values for growth curve so we can vectorize with case_when
    as_tibble() %>%
    mutate(repro_prob = case_when( # Store probability of reproduction
      age.x < 5 ~ 0, # No individuals younger than age 5 will reproduce (5 is an arbitrary number)
      age.x >= 5 ~ repro.prob(beta.0 = beta_0, beta.1 = beta_1, TLflex = indv_length),
      TRUE ~ NA))

  #head(init.pop2)

  # Confirmed - June 2024
  # Check that counts by ages and sex align with expectations
  # init.pop2 %>% dplyr::count(age.x, population) %>%
  #   arrange(population, age.x) %>%
  #   print(n = 100)
  #
  # init.pop2 %>% count(sex, population) %>%
  #   arrange(population, sex)

  ####----------Breeding----------####
  repro.cycle.vec <- rep(1:mating.periodicity, times = 100) # Generate a vector which will be used to determine if it is an even or odd breeding year (or a 1/3 breeding year)

  ####--------- For year 0 breeding
  #------------Mothers------------#
  # Mothers will have knife-edged maturity to allow the population to stay stable (at least for now)
  mothers <- init.pop2 %>% filter(sex=='F',
                                  age.x>=repro.age,
                                  repro.cycle == repro.cycle.vec[1])# Determine which females are available to breed in this year

  mothers <- mothers %>% mutate(num.mates = sample(num.mates, size = n(), replace = TRUE)) # Assign random number of mates to each mother

  # Make a new dataframe where each row corresponds to an instance of mating
  mothers2 <- mothers %>%
    lazy_dt() %>%
    group_by(indv.name) %>%
    slice(rep(1:n(), num.mates)) %>%
    ungroup() %>%
    select(indv.name, population) %>%
    rename(mother.x = indv.name) %>%
    as_tibble()

  #------------Fathers------------#
  fathers <- init.pop2 %>% filter(sex=='M',
                                  #init.pop$age.x>=repro.age # Uncomment for age-based maturity
                                  runif(n()) <= repro_prob
  ) %>% # Determine which fathers are available to breed in this year
    select(indv.name, population)

  if(popstructure == "structured"){

    # Create a list where each set of potential fathers are in different list elements corresponding to their population
    # Confirmed that this works
    fathers_by_population <- fathers %>%
      group_by(population) %>%
      summarise(indv.name = list(indv.name)) %>%
      deframe()

    # Create dataframe of mating events and generate initial offspring from each mating event
    YOY.df <- createYOY.init.byPop(mothers2, fathers_by_population, ff)

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
