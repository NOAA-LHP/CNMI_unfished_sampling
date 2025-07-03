### -----------------------------------------------------------------------------------------
#	June 2025
#
#  this script reads in the assumed population and life history parameters for each BMUS
#	from BMUS_LH_sim_pop_parameters.xlsx and the population simulation functions modified from
#	LH_Sampling to create the simulated populations used for the parameter estimation relative error
#	analysis.
### -----------------------------------------------------------------------------------------


#----- load packages, adjust options, and define working directory

# rm(list=ls())
   if (!require("pacman")) install.packages("pacman")
    pacman::p_load(reshape, dplyr, ggplot2, magrittr, assertthat, gridExtra, readxl, this.path)

# turn off scientific notation in console output
  options(scipen=999)
  this_dir <- this.path::here(.. = 0)

#----- read in the modified LH Sampling functions
  source(paste0(this_dir, "/modify_LH_sample_function.R"))

#----- if exists, read in previous simulated pops
  #   load(paste0(this_dir, "/simulated_pops.RData"))

#----- read in the excel sheet containing true population parameters
  all_BMUS_pop_parms <- as.data.frame(read_excel(paste0(this_dir, "/BMUS_LH_sim_pop_parameters.xlsx"), sheet = "parms"))	# str(all_BMUS_pop_parms)


#----- for each BMUS, simulate the population and save in an R space

###----- APRU
    APRU_pop_parms <- as.list(all_BMUS_pop_parms$APRU)
    names(APRU_pop_parms) <- all_BMUS_pop_parms$Parameter

    APRU <- do.call(simulate_population_harvest, APRU_pop_parms)

###----- ETCA
    ETCA_pop_parms <- as.list(all_BMUS_pop_parms$ETCA)
    names(ETCA_pop_parms) <- all_BMUS_pop_parms$Parameter

    ETCA <- do.call(simulate_population_harvest, ETCA_pop_parms)

###----- ETCO
    ETCO_pop_parms <- as.list(all_BMUS_pop_parms$ETCO)
    names(ETCO_pop_parms) <- all_BMUS_pop_parms$Parameter

    ETCO <- do.call(simulate_population_harvest, ETCO_pop_parms)

###----- PRAU
    PRAU_pop_parms <- as.list(all_BMUS_pop_parms$PRAU)
    names(PRAU_pop_parms) <- all_BMUS_pop_parms$Parameter

    PRAU <- do.call(simulate_population_harvest, PRAU_pop_parms)

###----- PRFI
    PRFI_pop_parms <- as.list(all_BMUS_pop_parms$PRFI)
    names(PRFI_pop_parms) <- all_BMUS_pop_parms$Parameter

    PRFI <- do.call(simulate_population_harvest, PRFI_pop_parms)

###----- PRFL
    PRFL_pop_parms <- as.list(all_BMUS_pop_parms$PRFL)
    names(PRFL_pop_parms) <- all_BMUS_pop_parms$Parameter

    PRFL <- do.call(simulate_population_harvest, PRFL_pop_parms)

###----- PRZO
    PRZO_pop_parms <- as.list(all_BMUS_pop_parms$PRZO)
    names(PRZO_pop_parms) <- all_BMUS_pop_parms$Parameter

    PRZO <- do.call(simulate_population_harvest, PRZO_pop_parms)


  save(APRU, ETCA, ETCO, PRAU, PRFI, PRFL, PRZO, all_BMUS_pop_parms, file = paste0(this_dir, "/simulated_pops_30Jun.RData"))

# 




#----- for each BMUS, perform the bootstrap sampling from the simulated populations
#		based on the samples we have available originally and save in an R space
#

## read in the samples we have
  have_lengths <- read.csv(paste0(this_dir,"\\cnmi_unfished_samples_update.csv"),header=T)
  have_lengths <- have_lengths[,c(2,5,10)]
  names(have_lengths)[] <- c('species','length','type')



###----- APRU

  full_name <- "Aphareus rutilans"
  bin_width <- 5
  have <- subset(have_lengths, species ==  full_name)		

  sample_plan <- APRU_original <- get_plan(have, bin_width, original=TRUE)
  sim_output <- APRU
  n_boots <- 1000
  age_max <- sim_output$parameters$age_max
  save_bootstraps <- FALSE

  APRU_sim_fit_org <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 


###----- ETCA

  full_name <- "Etelis carbunculus"
  bin_width <- 5
  have <- subset(have_lengths, species ==  full_name)		

  sample_plan <- ETCA_original <- get_plan(have, bin_width, original=TRUE)
  sim_output <- ETCA
  n_boots <- 1000
  age_max <- sim_output$parameters$age_max
  save_bootstraps <- FALSE

  ETCA_sim_fit_org <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 



###----- ETCO

  full_name <- "Etelis coruscans"
  bin_width <- 5
  have <- subset(have_lengths, species ==  full_name)		

  sample_plan <- ETCO_original <- get_plan(have, bin_width, original=TRUE)
  sim_output <- ETCO
  n_boots <- 1000
  age_max <- sim_output$parameters$age_max
  save_bootstraps <- FALSE

  ETCO_sim_fit_org <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 



###----- PRAU

  full_name <- "Pristipomoides auricilla"
  bin_width <- 2
  have <- subset(have_lengths, species ==  full_name)		

  sample_plan <- PRAU_original <- get_plan(have, bin_width, original=TRUE)
  sim_output <- PRAU
  n_boots <- 1000
  age_max <- sim_output$parameters$age_max
  save_bootstraps <- FALSE

  PRAU_sim_fit_org <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 


###----- PRFI

  full_name <- "Pristipomoides filamentosus"
  bin_width <- 2
  have <- subset(have_lengths, species ==  full_name)		

  sample_plan <- PRFI_original <- get_plan(have, bin_width, original=TRUE)
  sim_output <- PRFI
  n_boots <- 1000
  age_max <- sim_output$parameters$age_max
  save_bootstraps <- FALSE

  PRFI_sim_fit_org <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 


###----- PRFL

  full_name <- "Pristipomoides flavipinnis"
  bin_width <- 2
  have <- subset(have_lengths, species ==  full_name)		

  sample_plan <- PRFL_original <- get_plan(have, bin_width, original=TRUE)
  sim_output <- PRFL
  n_boots <- 1000
  age_max <- sim_output$parameters$age_max
  save_bootstraps <- FALSE

  PRFL_sim_fit_org <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 


###----- PRZO

  full_name <- "Pristipomoides zonatus"
  bin_width <- 2
  have <- subset(have_lengths, species ==  full_name)		

  sample_plan <- PRZO_original <- get_plan(have, bin_width, original=TRUE)
  sim_output <- PRZO
  n_boots <- 1000
  age_max <- sim_output$parameters$age_max
  save_bootstraps <- FALSE

  PRZO_sim_fit_org <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 


#  save

 save(APRU_sim_fit_org, ETCA_sim_fit_org, ETCO_sim_fit_org, 
	PRAU_sim_fit_org, PRFI_sim_fit_org, PRFL_sim_fit_org, 
	PRZO_sim_fit_org , file = paste0(this_dir, "/original_sample_fits_2July.RData"))


































####