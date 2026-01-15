# Create list containing input data for population simulation
create.pop.data <- function(max_age, survival, pop_size, pop_number, stable = T){
  # max_age should be a single number

  # IF wanting to decide own survival rates and age-specific population sizes, then
  ## stable = F
  ## survival should be a vector of length max_age + 1 (bc age 0 is included)
  ## pop_size should be a tibble with ncols=pop_number and nrows=max_age

  # IF wanting a stable population, then
  ## stable = T
  ## survival should be a vector of length max_age; the function will adjust survival of age 0 to stabilize population size
  ## pop_size will be created

  # Errors:
  ## if stable = F, then nrows(popsize) must equal max_age and ncols(popsize) must equal pop_number

}
