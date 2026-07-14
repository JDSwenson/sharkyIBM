#' Return to this script after finalizing the simulate_population.R script. Want to emulate simualte.pop exactly, but calculate s0, so I'll just manually figure out s0 values that work until I get a reasonable version of the script running.

# Search for "for simulate_pop" to find sections that belong in that script and not this one

calculate.s0 <- function(max_age,
                         survival,
                         pop_size,
                         mating_periodicity,
                         maturity_age,
                         litter_size,
                         num_mates,
                         num_years,
                         osr = 0.5,
                         infertility = 0,
                         F_by_age = F,
                         movement_array = NULL,
                         popstructure = "panmictic", #panmictic or "structured"
                         stickiness = NULL,
                         sticky_age = NULL,
                         sticky_interval = NULL,
                         superpod_size = NULL,
                         male_behavior = NULL, # "mischievous" or "strong_bull"
                         sampling = "random", # "random" or "superpod"
                         sample_size = NULL,
                         sample_years = NULL,
                         growth_params = NULL,
                         process_by = "age" # age or length
){


  #--------------------- Find initial starting value for s0 ---------------------
  #--------------------- Construct exploratory Leslie matrix
  # Initialize matrix
  A <- matrix(0, nrow = max_age + 1, ncol = max_age + 1)

  # How many female offspring per year?
  ff <- litter_size*female_fraction/mating_periodicity

  # Create fecundity vector
  f_vec <- c(rep(0, times = maturity_age), rep(ff, times = (max_age + 1 - maturity_age)))

  # Turn above values into a Leslie matrix
  for(i in 1:(max_age + 1)){

    A[1, i] <- f_vec[i]

    if(i <= max_age){

      A[i + 1, i] <- survival$survival_rate[i]

    }
  }

  # Calculate lambda
  lam_A <- lambda(A)
  if(lam_A > 1) trajectory <- "growing" else trajectory <- "declining"

  # Check to see if it's possible to achieve stable population growth by altering YOY survival
  A_test <- A

  # If the population is growing, then set s0 to 0.01 and see if the population declines; if the population is declining, then set s0 to 0.99 and see if the population grows. If not, then that means there is unlikely to be a reasonable value for s0 that will produce a stable population ("unlikely" only bc the realized values in the simulation will be different).

  if(trajectory == "growing") A_test[2, 1] <- 0.01 else A_test[2, 1] <- 0.99

  lam_A_test <- lambda(A_test)

  if(lam_A_test > 1) trajectory_test <- "growing" else trajectory_test <- "declining"

  # If the population doesn't switch from growing to declining after giving an extreme value to s0, then stop and request new parameters.
  if (trajectory == trajectory_test) {
    stop("No value of s0 can produce stable population growth given current fecundity and survival. Even extreme values for s0 result in a population that is ", trajectory, ". Please adjust survival and/or fecundity values.")
  } # End error message


  #--------------------- Find value of s0 that produces lambda = 1 in the Leslie matrix
  find_s0 <- function(s0, A, lambda_fun) {
    A[2, 1] <- s0
    lambda_fun(A) - 1
  }

  s0_est <- uniroot(
    find_s0,
    interval = c(0.001, 1),
    A = A,
    lambda_fun = lambda
  )$root

  # Confirm that adding the estimated value of s0 produces lambda ~ 1
  A[2, 1] <- s0_est
  lam_updated <- lambda(A)

  # Save stable age distribution of updated Leslie Matrix to use as starting value
  stable_A <- stable.stage(A)

  message("Found starting value for s0. Beginning optimization.")

  #------------------------ Data setup ------------------------
  # Theoretical stable age distribution
  init_pop_size <- tibble(
    age = c(0:max_age),
    N = round(stable_A*pop_size, 0)
  )
  f <- max(f_vec)

  # For simulate_pop
  # s_years <- c(((num_years - sample_years)+1): num_years) # which years will be sampled?

  # Prepare to make initial individual-based tibble
  init_ages <- rep(init_pop_size$age, times = init_pop_size$N)
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
  if(process_by == "length"){
    init_growth_params <- growth_params %>% dplyr::select(sex, population, L_inf, K) %>%
      mutate(length = L_inf/2)
  }

  # Make initial individual-based tibble
  init_pop <- tibble(
    indv_name   = sprintf("%020d", seq_len(sum(init_pop_size$N))),
    birth_year  = -1L,
    age         = init_ages,
    mother      = "xxxxx",
    father      = "xxxxx",
    sex         = init_sex
  ) %>%
    mutate(
      repro_cycle = if_else(
        sex == "F",
        init_repro_cycle[row_number()],
        NA_integer_
      ),
      fertile = init_fertile_vec
    )

  if(process_by == "length"){

    init_pop <- init_pop %>%
      left_join(init_growth_params, by = c("sex"))

  }

  # Add pods if including pod structure
  if(!is.null(stickiness)){

    # Assign pod as a random number between 1 and the number of mature females, then assign juveniles and mature males after
    init_pop <- init_pop %>%
      group_by(sex) %>%
      mutate(
        pod = if_else(
          sex == "F" & age >= maturity_age,
          sample.int(n()),
          NA_integer_
        )
      ) %>%
      ungroup()

    # Save pod numbers as a vector
    init_pods <- init_pop %>% filter(is.na(pod) == F) %>% pull(pod)

    # Assign pods to juveniles and males -- grouping by age keeps the pod distribution relatively even.
    init_pop <- init_pop %>% group_by(age) %>%
      mutate(
      pod = ifelse(
        sex == "M" | age < maturity_age,
        sample(init_pods, size = n(), replace = T),
        pod
      )
    ) %>% ungroup()

  }


  # Comment out all of below bc I don't think we want breeding the first year.

  ####----------Breeding----------####
  # repro_cycle_vec <- rep(1:mating_periodicity, times = num_years+1) # Generate a vector which will be used to determine if it is an even or odd breeding year (or a 1/3 breeding year)
  #
  # ####--------- For year 0 breeding
  # #------------Mothers------------#
  # # Mothers with knife-edged maturity to allow the population to stay stable (at least for now)
  # if(!is.null(maturity_age)){
  #   mothers <- init_pop %>% filter(sex == 'F',
  #                                  age >= maturity_age,
  #                                  fertile, # filters to only keep indvs with fertile == T, but is faster without the conditional statement
  #                                  repro_cycle == repro_cycle_vec[1]) # Determine which females are available to breed in this year
  # }
  #
  # # TO DO: create vector of mothers from age-specific fecundity vector
  # mothers <- mothers %>% mutate(n_mates = sample(num_mates, size = n(), replace = TRUE)) # Assign random number of mates to each mother
  #
  # # Make a new dataframe where each row corresponds to an instance of mating
  # # mothers2 <- mothers %>%
  # #   lazy_dt() %>%
  # #   group_by(indv_name) %>%
  # #   slice(rep(1:n(), n_mates)) %>%
  # #   ungroup() %>%
  # #   select(indv_name, population) %>%
  # #   rename(mother = indv_name) %>%
  # #   as_tibble()
  #
  # # Keep pods in mothers dataframe if including pod structure
  # if(!is.null(stickiness)){
  #
  #   mothers2 <- mothers %>%
  #     tidyr::uncount(n_mates, .remove = FALSE) %>%
  #     select(indv_name, population, pod) %>%
  #     rename(mother = indv_name)
  # } else {
  #
  #   mothers2 <- mothers %>%
  #     tidyr::uncount(n_mates, .remove = FALSE) %>%
  #     select(indv_name, population) %>%
  #     rename(mother = indv_name)
  # }
  #
  # # total litter size per reproductive female (guaranteed >= 1)
  # litters <- mothers2 %>%
  #   distinct(mother, population) %>%
  #   mutate(
  #     litter_size = 1 + rpois(n(), lambda = litter_size - 1)
  #   )
  #
  # #------------Fathers------------#
  # fathers <- init_pop %>% filter(sex=='M',
  #                                age >= maturity_age,
  #                                fertile
  # ) %>% # Uncomment for age-based maturity
  #   select(indv_name, population)
  #
  # # Create dataframe of mating events and generate initial offspring from each mating event
  # YOY_df <- create.YOY.init(mothers2, fathers, litters, process_by = process_by)
  #
  # # Add pods for fathers if including pod structure
  # if(!is.null(stickiness)){
  #
  #   pods_vec <- YOY_df %>% pull(pod)
  #
  #   # Assign each father to one of the pods of the females he mated with
  #   pod_fathers <- YOY_df %>% distinct(father, .keep_all = T) %>%
  #     select(indv_name = father, father_pod = pod) # could change to slice_sample if we want to make the pod totally random, but I don't think it matters
  #
  #   # Join dataframe of fathers with assigned pod to init_pop
  #   init_pop <- init_pop %>% left_join(pod_fathers, by = "indv_name") %>%
  #     mutate(pod = ifelse(
  #       is.na(father_pod) != T,
  #       father_pod,
  #       pod
  #     )) %>%
  #     select(-father_pod) %>%
  #     mutate(pod = ifelse( # Juvenile founding males do not have a pod yet, so we'll assign one here.
  #       is.na(pod) == T,
  #       sample(pods_vec, size = n(), replace = T),
  #       pod
  #     ))
  #
  # }

  # This dataframe holds the population at the end of the first year of the simulation
  #year_end_pop_0 <- bind_rows(init_pop, YOY_df)
  year_end_pop_0 <- init_pop

  ### END HERE June 4, 2026

  # Create superpods -- number of pods in each superpod is specified via superpod_size
  if(!is.null(stickiness)){

    pods <- year_end_pop_0 %>%
      distinct(pod) %>%
      mutate(
        superpod = rep(
          seq_len(ceiling(n() / superpod_size)),
          each = superpod_size
        )[sample(n())]
      )

    loopy_pop <- year_end_pop_0 %>%
      left_join(pods, by = "pod") %>%
      mutate(pod_year = 0)

    # If male_behavior %in% c("mischievous", "family_oriented"), then we do not need to do anything special to assign superpods for now
    if(male_behavior == "strong_bull"){

      if(process_by == "age"){

        # Randomly pick one father for each superpod
        bulls <- data1 %>% dplyr::filter(sex=='M',
                                         age >= maturity_age) %>%
          group_by(superpod) %>%
          slice_sample(n = 1) %>%
          ungroup() %>%
          select(indv_name, population, pod, superpod)

      }else if(process_by == "length"){

        bulls <- data1 %>% dplyr::filter(sex=='M',
                                         age >= maturity_age) %>%
          group_by(superpod) %>%
          slice_max(length, n = 1, with_ties = FALSE) %>%
          ungroup() %>%
          select(indv_name, population, pod, superpod)

      }
    }
  } else {

    loopy_pop <- year_end_pop_0

  }

  # At the end of year 0 ...
  message(paste("year 0, Population ", names(table(loopy_pop$population)), ": N_mothers=", table(mothers2$population), ", N_pups=", table(YOY_df$population), ", Total N=", table(loopy_pop$population) , sep=""))

  #############################################################`
  ####---------Loop through remaining simulation years-----####
  #############################################################`

  pop_size <- data.frame() # Initialize dataframe for storing population size

  loopy_list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year -- to save space, we won't populate this for now
  #plot_list <- list() # If we want to plot number of mature/breeding females/males each year for each iteration

  samples_df <- NULL

  parents_tibble <- tibble() # This will store info on offspring distribution per parent
  moms_temp = dads_temp <- NULL

  for(v in 1:num_years){ # Loop through all of the years in the simulation - the burn in and the years that matter

    # Assign survival or mortality based on age-specific survival probabilities
    data1 <- loopy_pop %>% left_join(survival_df, by = "age") %>%
      mutate(survival = runif(n()) <= survival_rate) %>%
      dplyr::filter(survival) %>%
      select(-c(survival_rate, survival)) %>%
      mutate(age = age + 1)
    #If individuals are older than max_age, they will be killed after they reproduce

    if(process_by == "length"){

      data1 <- data1 %>% mutate(length = update.length.vb(length, L_inf, K))
    }

    if(!is.null(stickiness)){

      # Index superpods so that we can ensure that all pods stay in the same superpod
      superpod_index <- data1 %>% distinct(pod, superpod)

      # Make separate dataframe of individuals that need to (potentially) change pods
      temp_pod_df <- data1 %>% filter(age >= sticky_age,
                                      v - pod_year == sticky_interval) %>%
        mutate(
          change_pod = runif(n()) > stickiness
        ) %>%
        select(-superpod) %>%
        mutate(pod = ifelse(
          change_pod == T,
          sample(superpod_index$pod, replace = T),
          pod
        )) %>%
        left_join(superpod_index, by = "pod")


      data1 <- data1 %>% filter(age < sticky_age |
                                  v - pod_year != sticky_interval) %>%
        bind_rows(temp_pod_df)

    }

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

    # Keep pods if including pod structure
    if(!is.null(stickiness)){

      mothers2 <- mothers %>%
        tidyr::uncount(n_mates, .remove = FALSE) %>%
        select(indv_name, population, pod, superpod) %>%
        rename(mother = indv_name)

    } else{

      mothers2 <- mothers %>%
        tidyr::uncount(n_mates, .remove = FALSE) %>%
        select(indv_name, population) %>%
        rename(mother = indv_name)

    }

    # total litter size per reproductive female (guaranteed >= 1)
    litters <- mothers2 %>%
      distinct(mother, population) %>%
      mutate(
        litter_size = 1 + rpois(n(), lambda = litter_size - 1)
      )

    if(!is.null(male_behavior)){

      if(male_behavior == "mischievous"){

        # Same as default behavior
        fathers <- data1 %>% dplyr::filter(sex=='M',
                                           age >= maturity_age) %>%
          select(indv_name, population, pod, superpod)

      } else if(male_behavior == "family_oriented"){
        # Same as default behavior
        fathers <- data1 %>% dplyr::filter(sex=='M',
                                           age >= maturity_age) %>%
          select(indv_name, population, pod, superpod)

      } else if(male_behavior == "strong_bull"){

        # which bulls are still alive?
        living_bulls <- bulls[bulls$indv_name %in% data1$indv_name,]

        if(process_by == "age"){

          # Randomly pick one bull per superpod IF there isn't already a living bull for that superpod
          bulls <- data1 %>% dplyr::filter(sex=='M',
                                           age >= maturity_age,
                                           !superpod %in% living_bulls$superpod) %>%
            group_by(superpod) %>%
            slice_sample(n = 1) %>%
            ungroup() %>%
            #        distinct(superpod, .keep_all = T) %>%
            select(indv_name, population, pod, superpod) %>%
            bind_rows(living_bulls)

          # Just one father for each superpod
          fathers <- bulls

        } else if(process_by == "length"){

          bulls <- data1 %>% dplyr::filter(sex=='M',
                                           age >= maturity_age,
                                           !superpod %in% living_bulls$superpod) %>%
            group_by(superpod) %>%
            slice_max(length, n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            select(indv_name, population, pod, superpod)

          # Just one father for each superpod
          fathers <- bulls

        }
      }
    } else {
      fathers <- data1 %>% dplyr::filter(sex=='M',
                                         age >= maturity_age) %>%
        select(indv_name, population)

    }

    # Confirm that each superpod has at least one mature male and one mature female
    if(!is.null(stickiness)){

      ### PICK UP HERE AFTER MAY 1:

      # Identify superpods that lack a reproductive male or female
      missing_superpods <- setdiff(
        union(fathers$superpod, mothers$superpod),
        intersect(fathers$superpod, mothers$superpod)
      )

      # Identify viable superpods i.e., superpods with both a reproductive male and female
      viable_superpods <- unique(data1$superpod[!data1$superpod %in% missing_superpods])

      # Move all indvs in the missing superpods to a viable superpod
      data1 <- data1 %>% mutate(
        superpod = ifelse(
          superpod %in% missing_superpods,
          sample(viable_superpods, size = n(), replace = T),
          superpod
        )
      )

      # Double-check
      # length(unique(data1$superpod[!data1$superpod %in% missing_superpods]))
      # length(unique(data1$superpod))
      # length(missing_superpods)

      # Save vector of updated superpods
      new_superpods <- data1 %>% dplyr::select(indv_name, superpod)

      # Update superpods in mothers and fathers list
      mothers2 <- mothers2 %>% rename(indv_name = mother, old_superpod = superpod) %>%
        left_join(new_superpods, by = "indv_name") %>%
        select(-old_superpod) %>%
        rename(mother = indv_name)

      fathers <- fathers %>% rename(old_superpod = superpod) %>%
        left_join(new_superpods, by = "indv_name") %>%
        select(-old_superpod)

      # Double-check
      # mothers2 %>% rename(indv_name = mother, old_superpod = superpod) %>%
      #   left_join(new_superpods, by = "indv_name") %>%
      #   filter(old_superpod %in% missing_superpods)

    }
    # Create dataframe of mating events and generate initial offspring from each mating event
    YOY_df <- create.YOY(mothers2, fathers, litters, year = v, process_by = process_by)

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

    loopy_pop %>% count(age)

    ###############################################`
    ####---------------Sampling----------------####
    ###############################################`
    if(v %in% s_years){

      if(sampling == "random"){

        samples_df_temp <- loopy_pop %>%
          slice_sample(n = sample_size) %>%
          mutate(capture_year = v)


      } else if(sampling == "superpod"){

        # Step 1: choose superpods
        chosen_pods <- loopy_pop %>%
          distinct(superpod) %>%
          slice_sample(n = length(sample_size)) %>%
          mutate(n_to_sample = sample_size)

        # Step 2: sample within each superpod
        samples_df_pod <- loopy_pop %>%
          inner_join(chosen_pods, by = "superpod")


        # Split data by superpod (only chosen ones are present)
        pod_list <- split(samples_df_pod, samples_df_pod$superpod)

        # Named vector of sample sizes (superpod → n)
        n_vec <- chosen_pods$n_to_sample
        names(n_vec) <- chosen_pods$superpod

        # Sample within each superpod
        sampled_list <- Map(
          function(df, n) df[sample.int(nrow(df), n), ],
          pod_list,
          n_vec
        )

        # Recombine into a single data frame
        samples_df_temp <- dplyr::bind_rows(sampled_list)

      }

      samples_df <- bind_rows(samples_df, samples_df_temp)

    } # End sampling


    ###############################################`
    ####-------------Save metrics--------------####
    ###############################################`

    # Calculate number of produced offspring per mother and father this year
    # moms_temp <- YOY_df %>%
    #   count(mother, population, name = "num_off") %>%
    #   mutate(
    #     indv_name = mother,
    #     population,
    #     num_off,
    #     year = v,
    #     which_parent = "mother",
    #     .keep = "none"
    #   )
    #
    # dads_temp <- YOY_df %>%
    #   count(father, population, name = "num_off") %>%
    #   mutate(
    #     indv_name = father,
    #     population,
    #     num_off,
    #     year = v,
    #     which_parent = "father",
    #     .keep = "none"
    #   )

    # Add to the tibble of offspring distribution - can use to check if/whether some indvs are reproducing much more
    #   parents.tibble <- rbind(parents.tibble, moms.temp, dads.temp)

    # Print info about the population to the console
    message(paste("\nyear", v, " ", names(table(loopy_pop$population)),
                  "N_mothers=", table(moms_temp$population),
                  "N_fathers=", table(dads_temp$population),
                  "\nN_pups=", table(YOY_df$population),
                  "\nTotal N= ", table(loopy_pop$population), sep=" "))

    # Save the population size by age and sex
    pop_size_temp <- loopy_pop %>% dplyr::count(population, sex, age) %>%
      mutate(year = v)

    pop_size <- bind_rows(pop_size, pop_size_temp)

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


} # End calculate.s0 function
