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
                         num_mates,
                         num_years,
                         movement_array = NULL,
                         popstructure = "panmictic" #panmictic or "structured"
                         ){

  # Save initial values as distinct R objects
  init_pop_size <- input_data$numbers_at_age
  max_age <- max(input_data$numbers_at_age$age)
  f <- max(input_data$fecundity)
  litter_size <- input_data$litter_size

  # Save growth parameters as a tibble and specify population as the corresponding list index
  growth_params <- purrr::imap(
    input_data$growth_params,
    ~ dplyr::mutate(.x, population = .y)
  ) %>%
    purrr::list_rbind()

  infertility <- input_data$infertility
  mating_periodicity <- input_data$mating_periodicity
  female_fraction <- input_data$female_fraction
  maturity_age <- input_data$maturity_age

#  YOY_survival <- input_data$s0
#  juvenile_survival <- input_data$survival[2:(maturity_age)]
#  adult_survival <- input_data$survival[maturity_age:(max_age-1)]
#  ff <- f / female_fraction # Female fecundity per breeding event at equilibrium
  survival_df <- tibble(age = c(0:max_age), survival_rate = input_data$survival)

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
  init_fertile_vec <- runif(n = length(init_repro_cycle)) > infertility

  # The simulation should run long enough to get rid of founders, so we just need placeholders for length.
  # Here, we'll define length somewhat arbitrarily as L_inf/2
  init_growth_params <- growth_params %>% dplyr::select(sex, population, L_inf, K) %>%
    mutate(length = L_inf/2)

  # Summarize population numbers
  total_pop_sizes_df <- init_pop_size %>% group_by(population) %>%
    reframe(total = sum(N, na.rm = T)) %>%
    arrange(population)

  total_pop_sizes_vec <- rep(total_pop_sizes_df$population, times = total_pop_sizes_df$total)

  ###############################################`
  ####---------Set up initial population-----####
  ###############################################`
  # Initial population
# init_pop <- tibble(
#   indv_name   = sprintf("%020d", seq_len(sum(init_pop_size$N))),
#   birth_year  = -1L,
#   age         = init_ages,
#   mother      = "xxxxx",
#   father      = "xxxxx",
#   sex         = init_sex,
#   population  = total_pop_sizes_vec
# ) %>%
#   mutate(
#     repro_cycle = if_else(
#       sex == "F",
#       init_repro_cycle[row_number()],
#       NA_integer_
#     ),
#     fertile = init_fertile_vec
#   )

  init_pop <- tibble(
    indv_name   = sprintf("%020d", seq_len(sum(init_pop_size$N))),
    birth_year  = -1L,
    age         = init_ages,
    mother      = "xxxxx",
    father      = "xxxxx",
    sex         = init_sex,
    population  = total_pop_sizes_vec
  ) %>%
    mutate(
      repro_cycle = if_else(
        sex == "F",
        init_repro_cycle[row_number()],
        NA_integer_
      ),
      fertile = init_fertile_vec
    ) %>%
    left_join(init_growth_params, by = c("sex", "population"))

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
  repro_cycle_vec <- rep(1:mating_periodicity, times = num_years+1) # Generate a vector which will be used to determine if it is an even or odd breeding year (or a 1/3 breeding year)

  ####--------- For year 0 breeding
  #------------Mothers------------#
  # Mothers with knife-edged maturity to allow the population to stay stable (at least for now)
  if(!is.null(maturity_age)){
  mothers <- init_pop %>% filter(sex == 'F',
                                 age >= maturity_age,
                                 fertile, # filters to only keep indvs with fertile == T, but is faster without the conditional statement
                                 repro_cycle == repro_cycle_vec[1]) # Determine which females are available to breed in this year
  }

  # TO DO: create vector of mothers from age-specific fecundity vector
  mothers <- mothers %>% mutate(n_mates = sample(num_mates, size = n(), replace = TRUE)) # Assign random number of mates to each mother

  # Make a new dataframe where each row corresponds to an instance of mating
  # mothers2 <- mothers %>%
  #   lazy_dt() %>%
  #   group_by(indv_name) %>%
  #   slice(rep(1:n(), n_mates)) %>%
  #   ungroup() %>%
  #   select(indv_name, population) %>%
  #   rename(mother = indv_name) %>%
  #   as_tibble()

  mothers2 <- mothers %>%
    tidyr::uncount(n_mates, .remove = FALSE) %>%
    select(indv_name, population) %>%
    rename(mother = indv_name)

  # total litter size per reproductive female (guaranteed >= 1)
  litters <- mothers2 %>%
    distinct(mother, population) %>%
    mutate(
      litter_size = 1 + rpois(n(), lambda = litter_size - 1)
    )

  #------------Fathers------------#
  fathers <- init_pop %>% filter(sex=='M',
                                 age >= maturity_age
                                 ) %>% # Uncomment for age-based maturity
    select(indv_name, population)

    # Create dataframe of mating events and generate initial offspring from each mating event
    YOY_df <- create.YOY(mothers2, fathers, litters, year = 0)

  # This dataframe holds the population at the end of the first year of the simulation
  year_end_pop_0 <- bind_rows(init_pop, YOY_df)

  # And finally, we assign age-specific mortality rates to each individual and then determine whether they survive into the next year or not
  loopy_pop <- year_end_pop_0

  # At the end of year 0 ...
  message(paste("year 0, Population ", names(table(loopy_pop$population)), ": N_mothers=", table(mothers2$population), ", N_pups=", table(YOY_df$population), ", Total N=", table(loopy_pop$population) , sep=""))


  #############################################################`
  ####---------Loop through remaining simulation years-----####
  #############################################################`

  pop_size <- data.frame() # Initialize dataframe for storing population size

  loopy_list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year -- to save space, we won't populate this for now
  plot_list <- list() # If we want to plot number of mature/breeding females/males each year for each iteration

  samples_df <- NULL

  parents_tibble <- tibble() # This will store info on offspring distribution per parent
  moms_temp = dads_temp <- NULL

  for(v in 1:num_years){ # Loop through all of the years in the simulation - the burn in and the years that matter

    # Assign survival or mortality based on age-specific survival probabilities
    data1 <- loopy_pop %>% left_join(survival_df, by = "age") %>%
      mutate(survival = runif(n()) <= survival_rate) %>%
      dplyr::filter(survival) %>%
      mutate(age = age + 1,
             length = update.length.vb(length, L_inf, K)) %>%
      select(-c(survival_rate, survival))
    #If individuals are older than max_age, they will be killed after they reproduce

    ####----------Breeding----------####
    #------------Mothers------------#
    mothers <- data1 %>% filter(sex == 'F',
                                age >= maturity_age,
                                fertile, # filters to only keep indvs with fertile == T, but is faster without the conditional statement
                                repro_cycle == repro_cycle_vec[v + 1]) # Determine which females are available to breed in this year

    # Add column that contains the number of mates each mother will mate with this year
    mothers <- mothers %>% mutate(n_mates = sample(num_mates, size = nrow(mothers), replace = TRUE))

    # Make a new row where each row corresponds to an instance of mating
    # mothers2 <- mothers %>%
    #   lazy_dt() %>%
    #   group_by(indv_name) %>%
    #   slice(rep(1:n(), n_mates)) %>%
    #   ungroup() %>%
    #   select(indv_name, population) %>%
    #   rename(mother = indv_name) %>%
    #   as_tibble()

    mothers2 <- mothers %>%
      tidyr::uncount(n_mates, .remove = FALSE) %>%
      select(indv_name, population) %>%
      rename(mother = indv_name)

    # total litter size per reproductive female (guaranteed >= 1)
    litters <- mothers2 %>%
      distinct(mother, population) %>%
      mutate(
        litter_size = 1 + rpois(n(), lambda = litter_size - 1)
      )

    #------------Fathers------------#
    fathers <- data1 %>% dplyr::filter(sex=='M',
                                       age >= maturity_age
    ) %>% # Uncomment for age-based maturity
      select(indv_name, population)

    # Create dataframe of mating events and generate initial offspring from each mating event
    YOY_df <- create.YOY(mothers2, fathers, litters, year = v)

    #Only bother assigning a sampling location for the years we're taking samples; otherwise just slows down code.
    ## TO DO: ADD SAMPLE YEAR AND SITE INFORMATION THEN UNCOMMENT AND UPDATE BELOW

    # if(v >= min(sample_years)){
    #
    #   # YOY from the same mother are assigned to the same sampling location
    #   YOY_df <- YOY_df %>% group_by(mother, population) %>%
    #     mutate(
    #       sampling_location = sample(
    #         sampling_locations,
    #         1,
    #         prob = c(dispersal_kernel(age = 0, birth_population = population[1])),
    #         replace = TRUE)
    #     ) %>%
    #     ungroup()
    #
    #   #Pull out mothers and assign them the same sampling location as their offspring from this year
    #   mother_sample_df <- YOY_df %>% select(indv_name = mother, sampling_location) %>% distinct()
    #
    #   mothers_df <- mother_sample_df %>% lazy_dt() %>%
    #     left_join(mothers, by = "indv_name") %>%
    #     select(-num_mates) %>%
    #     as_tibble()
    #
    #   #Assign all other individuals assigned randomly
    #   loopy_pop <- data1 %>%
    #     lazy_dt() %>%
    #     filter(!indv_name %chin% YOY_df$mother) %>%
    #     mutate(
    #       sampling_location = map2_chr(
    #         age,
    #         population,
    #         ~sample(
    #           sampling_locations,
    #           1,
    #           prob = c(dispersal_kernel(.x, .y)),
    #           replace = TRUE
    #         ))
    #     ) %>%
    #     as_tibble() %>%
    #     bind_rows(YOY_df, mothers_df)
    #
    #
    # } else {
    #
    #   # No need to assign sampling location if we're not sampling this year
       loopy_pop <- bind_rows(data1, YOY_df)
    #
    # }

    ###############################################`
    ####---------------Sampling----------------####
    ###############################################`
    # if(v %in% sample_years){
    #
    #   samples_df_temp <- loopy_pop %>%
    #     group_by(sampling_location) %>%
    #     group_map(~sample_fixed(.x, samples_vec[.y$sampling_location[1]]), .keep = TRUE) %>%
    #     bind_rows() %>%
    #     mutate(capture_year = v)
    #
    #   samples_df <- bind_rows(samples_df, samples_df_temp)
    #
    # }


    ###############################################`
    ####-------------Save metrics--------------####
    ###############################################`

    # Calculate number of produced offspring per mother and father this year
    moms_temp <- YOY_df %>%
      count(mother, population, name = "num_off") %>%
      mutate(
        indv_name = mother,
        population,
        num_off,
        year = v,
        which_parent = "mother",
        .keep = "none"
      )

    dads_temp <- YOY_df %>%
      count(father, population, name = "num_off") %>%
      mutate(
        indv_name = father,
        population,
        num_off,
        year = v,
        which_parent = "father",
        .keep = "none"
      )

    # Add to the tibble of offspring distribution - can use to check if/whether some indvs are reproducing much more
    #   parents.tibble <- rbind(parents.tibble, moms.temp, dads.temp)

    # Print info about the population to the console
    message(paste("\nyear", v, " ", names(table(loopy_pop$population)),
              "N_mothers=", table(moms_temp$population),
              "N_fathers=", table(dads_temp$population),
              "\nN_pups=", table(YOY_df$population),
              "\nTotal N= ", table(loopy_pop$population), sep=" "))

    # Save the population size by age and sex
    pop_size <- loopy_pop %>% dplyr::count(population, sex, age)

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

  return(invisible(list(pop_size, samples_df)))
}
