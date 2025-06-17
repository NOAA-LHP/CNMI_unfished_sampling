### -----------------------------------------------------------------------------------------
#	June 16 2025, write some codes to explore which length bins we need samples from
#		for the unfished areas of CNMI
#
#
### -----------------------------------------------------------------------------------------


# rm(list=ls())
 if (!require("pacman")) install.packages("pacman")
  pacman::p_load(reshape, dplyr, ggplot2, magrittr, assertthat)


this_dir <- "C:\\Users\\Erin.Bohaboy\\Documents\\CNMI\\NMarianas_Cruise_Jul2025\\CNMI_unfished_sampling"


## read in the existing LH_sampling functions
source(paste0(this_dir, "\\LHSample_functions.R"))


# ------- first, read in the .csv of samples by species and length that we already have, simplify

have_lengths <- read.csv(paste0(this_dir,"\\cnmi_unfished_samples.csv"),header=T)
BMUS_names <- c("Aphareus rutilans", "Caranx ignobilis", "Caranx lugubris", "Etelis carbunculus", "Etelis coruscans", 
				"Lethrinus rubrioperculatus", "Lutjanus kasmira", "Pristipomoides auricilla", 
				"Pristipomoides filamentosus", "Pristipomoides flavipinnis", "Pristipomoides sieboldii", 
				"Pristipomoides zonatus", "Variola louti")

have_lengths <- subset(have_lengths, ScientificName %in% BMUS_names)
summary(as.factor(have_lengths$ScientificName))			# there are no ignobilis, rubrio, kasmira, or sieboldii

have_lengths <- have_lengths[,c(2,5)]
names(have_lengths)[] <- c('species','length')

# ------- 

# ---  Using onaga as an example, what do we have, what would POS look like for an unfished population, and what do we need?

# -- 1. 
  have <- subset(have_lengths, species == 'Etelis coruscans')
  hist_have <- hist(have$length, breaks = seq(0,120,5),include.lowest=TRUE, right=FALSE,plot=TRUE)

  we_have <- data.frame('binL' = hist_have$breaks[1:(length(hist_have$breaks)-1)], 'have' = hist_have$counts)


# -- 2. simulate population
  onaga <- simulate_population_harvest(Linf = 100, Linf_sd = 2.5, 
                                     M = 0.125, F = 0.01, Lorenzen = TRUE, 
                                     mincat = 10, catsd = 2.5, maxcat = 200, maxcatsd = 0, 
                                     L0 = 10, L0_sd = 2.5, k = 0.14, k_sd = 0, 
                                     Amax = 55, age_max = 50, N = 100000)

  # take a look
  # onaga_hist_pop <- hist(onaga$population$length, breaks = seq(0,120,5),include.lowest=TRUE, right=FALSE,plot=TRUE)
  #  this is mildly annoying because you can see ages 0 to 3 in the length comps, then a lot of fish build up around Linf

# -- step 3. what would a POS approach look like, if we could sample directly from the population?
#  	use the find_POS function

  plan <- find_POS(sim_output = onaga, samp_size = 300, Lbin_width = 5)		#  str(plan)


# -- 4. put together do figures and tables.

  summary <- merge(x=plan, y=we_have, by= 'binL', all.x = TRUE)
  summary$need <- pmax(0,(summary$plan-summary$have))


### Erin stop here 16June















































####