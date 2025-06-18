### -----------------------------------------------------------------------------------------
#	June 16 2025, write some codes to explore which length bins we need samples from
#		for the unfished areas of CNMI
#
#
### -----------------------------------------------------------------------------------------


# rm(list=ls())
 if (!require("pacman")) install.packages("pacman")
  pacman::p_load(reshape, dplyr, ggplot2, magrittr, assertthat, gridExtra)


this_dir <- "C:\\Users\\Erin.Bohaboy\\Documents\\CNMI\\NMarianas_Cruise_Jul2025\\CNMI_unfished_sampling"


## read in LH_sampling functions
# source(paste0(this_dir, "\\LHSample_functions.R"))


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



# ===== QUESTION 1: ===================
#  Using onaga as an example, what do we have, what would POS look like for an unfished population, and what do we need?

# -- 1. 
  have <- subset(have_lengths, species == 'Etelis coruscans')
  hist_have <- hist(have$length, breaks = seq(0,120,5),include.lowest=TRUE, right=FALSE,plot=TRUE)

  we_have <- data.frame('binL' = hist_have$breaks[1:(length(hist_have$breaks)-1)], 'have' = hist_have$counts)


# -- 2. simulate population
  onaga <- simulate_population_harvest(Linf = 100, Linf_sd = 2.5, 
                                     M = 0.125, F = 0.01, Lorenzen = TRUE, 
                                     mincat = 10, catsd = 2.5, maxcat = 200, maxcatsd = 0, 
                                     L0 = 10, L0_sd = 2.5, k = 0.14, k_sd = 0.008, 
                                     Amax = 55, age_max = 20, N = 100000, Linf_k_cor_TF = TRUE)

  # take a look
  # onaga_hist_pop <- hist(onaga$population$length, breaks = seq(0,120,5),include.lowest=TRUE, right=FALSE,plot=TRUE)
  #  this is mildly annoying because you can see ages 0 to 3 in the length comps, then a lot of fish build up around Linf
  #  save this for later to save time if still doing dev 
  # save(onaga, file = paste0(this_dir, "\\onaga_sim_dev.RData"))


# -- step 3. what would a POS approach look like, if we could sample directly from the population?
#  	use the find_POS function

  plan <- find_POS(sim_output = onaga, samp_size = 300, Lbin_width = 5)		#  str(plan)


# -- 4. put together do figures and tables.

  summary <- merge(x=plan, y=we_have, by= 'binL', all.x = TRUE)
  summary$need <- pmax(0,(summary$plan-summary$have))
  summary$best <- pmax(summary$plan, summary$have)

# so....  that's the answer. If using a strict POS, if we want N = whatever, here we are.


# ===== QUESTION 2: ===================
#  A more interesting question is, if we fit von B using a given sampling plan, what would the error be?

# -- 1. repeats step 1 and 2 from above
# -- 2. 

# -- 3. define the sample plan (lower bin edge and number of samples per bin).
#		this might be what we have or what we want.

 sample_plan <- data.frame(binL = summary$binL, nsamps = summary$have)		# str(sample_plan)

 # use a new function to bootstrap sample the population according to this plan and fit the von Bert

  sim_output <- onaga
  n_boots <- 100
  age_max <- sim_output$parameters$age_max		# define age max which is for CV_L_old
  save_bootstraps <- FALSE

  sample_plan_1_results <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 

  names(sample_plan_1_results)

sample_plan_1_results$parameter_summary_all_boots
sample_plan_1_results$pop_parms


# -- 4. make some figures

   plan <- sample_plan_1_results$params_input_output$sample_plan 
   parm_summary <- sample_plan_1_results$parameter_summary_all_boots			#str(parm_summary)
   parm_pop <- sample_plan_1_results$pop_parms
    
   ticks <- plan$binL
   ymax <- max(plan$nsamps)
   xmax <- max(plan$binL)
   tick_height <- ymax/50

   Tot_N <- paste0("N = ",sum(plan$nsamps))


    p_plan <- ggplot(data=plan, aes(x=binL, y=nsamps )) +
		    geom_bar(stat="identity",fill='#008998') +
		    scale_x_continuous(n.breaks=10) +
		    theme_LH_bar() +
		    annotate("text", x=(0.1*xmax), y=(0.9*ymax), label = Tot_N) +
 		    labs(title="", subtitle="", y="Number samples", x="Lower edge of length bin", caption="")
    
    p_Linf <- ggplot(data=subset(parm_summary, parm_name=='Linf'), aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = parm_pop$Linf, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="cm", x="", caption="")

    p_k <- ggplot(data=subset(parm_summary, parm_name=='K1'), aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = parm_pop$K1, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="", x="", caption="")

    p_L0 <- ggplot(data=subset(parm_summary, parm_name=='L0'), aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = parm_pop$L0, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="cm", x="", caption="")

    ggsave(paste0(this_dir,"/figures/onaga_1A.png"),grid.arrange(grobs=list(p_plan), ncol=1),width=7.5, height=3, units="in")	
    ggsave(paste0(this_dir,"/figures/onaga_1B.png"),grid.arrange(grobs=list(p_Linf, p_k, p_L0), ncol=3),width=7.5, height=5, units="in")
	 
    


# -- repeat steps 3 and 4 for a different sampling approach

 sample_plan <- data.frame(binL = summary$binL, nsamps = summary$best)		# str(sample_plan)

 # use a new function to bootstrap sample the population according to this plan and fit the von Bert

  sim_output <- onaga
  n_boots <- 100
  age_max <- sim_output$parameters$age_max		# define age max which is for CV_L_old
  save_bootstraps <- FALSE

  sample_plan_2_results <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 

# -- 4. make some figures

   plan <- sample_plan_2_results$params_input_output$sample_plan 
   parm_summary <- sample_plan_2_results$parameter_summary_all_boots			#str(parm_summary)
   parm_pop <- sample_plan_2_results$pop_parms
    
   ticks <- plan$binL
   ymax <- max(plan$nsamps)
   xmax <- max(plan$binL)
   tick_height <- ymax/50

   Tot_N <- paste0("N = ",sum(plan$nsamps))


    p_plan <- ggplot(data=plan, aes(x=binL, y=nsamps )) +
		    geom_bar(stat="identity",fill='#008998') +
		    scale_x_continuous(n.breaks=10) +
		    theme_LH_bar() +
		    annotate("text", x=(0.1*xmax), y=(0.9*ymax), label = Tot_N) +
 		    labs(title="", subtitle="", y="Number samples", x="Lower edge of length bin", caption="")
    
    p_Linf <- ggplot(data=subset(parm_summary, parm_name=='Linf'), aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = parm_pop$Linf, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="cm", x="", caption="")

    p_k <- ggplot(data=subset(parm_summary, parm_name=='K1'), aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = parm_pop$K1, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="", x="", caption="")

    p_L0 <- ggplot(data=subset(parm_summary, parm_name=='L0'), aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = parm_pop$L0, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="cm", x="", caption="")

    ggsave(paste0(this_dir,"/figures/onaga_2A.png"),grid.arrange(grobs=list(p_plan), ncol=1),width=7.5, height=3, units="in")	
    ggsave(paste0(this_dir,"/figures/onaga_2B.png"),grid.arrange(grobs=list(p_Linf, p_k, p_L0), ncol=3),width=7.5, height=5, units="in")
	 
    




































####