### -----------------------------------------------------------------------------------------
#	June 16 2025, write some codes to explore which length bins we need samples from
#		for the unfished areas of CNMI
#
#
### -----------------------------------------------------------------------------------------


# rm(list=ls())
 if (!require("pacman")) install.packages("pacman")
  pacman::p_load(reshape, dplyr, ggplot2, magrittr, assertthat, gridExtra, readxl)

# turn off scientific notation in console output
options(scipen=999)


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


# ----- read in our assumptions about the true population parameters

  all_BMUS_pop_parms <- as.data.frame(read_excel(paste0(this_dir, "\\BMUS_LH_sim_pop_parameters.xlsx"), sheet = "parms"))

  # pick the BMUS we are working with, put parameters into a named list
  these_pop_parms <- as.list(all_BMUS_pop_parms$ETCO)
  names(these_pop_parms) <- all_BMUS_pop_parms$Parameter


# ===== QUESTION 1: ===================
#  Using onaga as an example, what do we have, what would POS look like for an unfished population, and what do we need?

# -- 1. 
  have <- subset(have_lengths, species == 'Etelis coruscans')
  # hist_have <- hist(have$length, breaks = seq(0,120,5),include.lowest=TRUE, right=FALSE,plot=TRUE)

  we_have <- data.frame('binL' = hist_have$breaks[1:(length(hist_have$breaks)-1)], 'have' = hist_have$counts)


# -- 2. simulate population
  onaga <- do.call(simulate_population_harvest, these_pop_parms)		#onaga$parameters

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
    



    p_plan <- ggplot(data=plan, aes(x=binL, y=nsamps )) +
		    geom_bar(stat="identity",fill='#008998') +
		    scale_x_continuous(n.breaks=10) +
		    theme_LH_bar() +
		    annotate("text", x=(0.1*max(plan$binL)), y=(0.9*max(plan$nsamps)), label = paste0("N = ",sum(plan$nsamps))) +
 		    labs(title="", subtitle="", y="Number samples", x="Lower edge of length bin", caption="")
 
#  error in absolute
   
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
    ggsave(paste0(this_dir,"/figures/onaga_1B.png"),grid.arrange(grobs=list(p_k, p_L0, p_Linf), ncol=3),width=5, height=2, units="in")
	 
  
#  relative error  (obs-true/true)

   parm_summary_raw <- sample_plan_1_results$parameter_summary_all_boots[c(1,2,6),]			#str(parm_summary)
   parm_pop <- as.numeric(sample_plan_1_results$pop_parms[c(1,2,4)])
    
   rel_summary <- data.frame(parm_name = parm_summary_raw$parm_name, 
					lower95 = (parm_summary_raw[,"lower95"] - parm_pop)/parm_pop,
					lower50 = (parm_summary_raw[,"lower50"] - parm_pop)/parm_pop,
					avg = (parm_summary_raw[,"avg"] - parm_pop)/parm_pop,
					upper50 = (parm_summary_raw[,"upper50"] - parm_pop)/parm_pop,   
					upper95 = (parm_summary_raw[,"upper95"] - parm_pop)/parm_pop)

   p_rel <- ggplot(data=rel_summary, aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="relative error", x="", caption="")

   ggsave(paste0(this_dir,"/figures/onaga_1B_rel.png"),p_rel,width=5, height=2, units="in")
	

# -- repeat steps 3 and 4 for an ideal sampling approach

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
    


    p_plan <- ggplot(data=plan, aes(x=binL, y=nsamps )) +
		    geom_bar(stat="identity",fill='#008998') +
		    scale_x_continuous(n.breaks=10) +
		    theme_LH_bar() +
		    annotate("text", x=(0.1*max(plan$binL)), y=(0.9*max(plan$nsamps)), label = paste0("N = ",sum(plan$nsamps))) +
 		    labs(title="", subtitle="", y="Number samples", x="Lower edge of length bin", caption="")
 
#  error in absolute
   
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
    ggsave(paste0(this_dir,"/figures/onaga_2B.png"),grid.arrange(grobs=list(p_k, p_L0, p_Linf), ncol=3),width=5, height=2, units="in")
	 
  
#  relative error  (obs-true/true)

   parm_summary_raw <- sample_plan_2_results$parameter_summary_all_boots[c(1,2,6),]			#str(parm_summary)
   parm_pop <- as.numeric(sample_plan_2_results$pop_parms[c(1,2,4)])
    
   rel_summary <- data.frame(parm_name = parm_summary_raw$parm_name, 
					lower95 = (parm_summary_raw[,"lower95"] - parm_pop)/parm_pop,
					lower50 = (parm_summary_raw[,"lower50"] - parm_pop)/parm_pop,
					avg = (parm_summary_raw[,"avg"] - parm_pop)/parm_pop,
					upper50 = (parm_summary_raw[,"upper50"] - parm_pop)/parm_pop,   
					upper95 = (parm_summary_raw[,"upper95"] - parm_pop)/parm_pop)

   p_rel <- ggplot(data=rel_summary, aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="relative error", x="", caption="")

   ggsave(paste0(this_dir,"/figures/onaga_2B_rel.png"),p_rel,width=5, height=2, units="in")
	



# -- repeat steps 3 and 4 for an imaginary sample addition

 sample_plan <- data.frame(binL = summary$binL, nsamps = summary$have)	# add in some samples
 sample_plan$nsamps[6] <- 3

 # use a new function to bootstrap sample the population according to this plan and fit the von Bert

  sim_output <- onaga
  n_boots <- 100
  age_max <- sim_output$parameters$age_max		# define age max which is for CV_L_old
  save_bootstraps <- FALSE

  sample_plan_3_results <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 

# -- 4. make some figures

   plan <- sample_plan_3_results$params_input_output$sample_plan 
   parm_summary <- sample_plan_3_results$parameter_summary_all_boots			#str(parm_summary)
   parm_pop <- sample_plan_3_results$pop_parms
    

    p_plan <- ggplot(data=plan, aes(x=binL, y=nsamps )) +
		    geom_bar(stat="identity",fill='#008998') +
		    scale_x_continuous(n.breaks=10) +
		    theme_LH_bar() +
		    annotate("text", x=(0.1*max(plan$binL)), y=(0.9*max(plan$nsamps)), label = paste0("N = ",sum(plan$nsamps))) +
 		    labs(title="", subtitle="", y="Number samples", x="Lower edge of length bin", caption="")
 
#  error in absolute
   
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

    ggsave(paste0(this_dir,"/figures/onaga_3A.png"),grid.arrange(grobs=list(p_plan), ncol=1),width=7.5, height=3, units="in")	
    ggsave(paste0(this_dir,"/figures/onaga_3B.png"),grid.arrange(grobs=list(p_k, p_L0, p_Linf), ncol=3),width=5, height=2, units="in")
	 
  
#  relative error  (obs-true/true)

   parm_summary_raw <- sample_plan_3_results$parameter_summary_all_boots[c(1,2,6),]			#str(parm_summary)
   parm_pop <- as.numeric(sample_plan_3_results$pop_parms[c(1,2,4)])
    
   rel_summary <- data.frame(parm_name = parm_summary_raw$parm_name, 
					lower95 = (parm_summary_raw[,"lower95"] - parm_pop)/parm_pop,
					lower50 = (parm_summary_raw[,"lower50"] - parm_pop)/parm_pop,
					avg = (parm_summary_raw[,"avg"] - parm_pop)/parm_pop,
					upper50 = (parm_summary_raw[,"upper50"] - parm_pop)/parm_pop,   
					upper95 = (parm_summary_raw[,"upper95"] - parm_pop)/parm_pop)

   p_rel <- ggplot(data=rel_summary, aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="relative error", x="", caption="")

   ggsave(paste0(this_dir,"/figures/onaga_3B_rel.png"),p_rel,width=5, height=2, units="in")
	




# -- repeat steps 3 and 4 for perfect POS only

 sample_plan <- data.frame(binL = summary$binL, nsamps = summary$plan)

 # use a new function to bootstrap sample the population according to this plan and fit the von Bert

  sim_output <- onaga
  n_boots <- 100
  age_max <- sim_output$parameters$age_max		# define age max which is for CV_L_old
  save_bootstraps <- FALSE

  sample_plan_4_results <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 

# -- 4. make some figures

   plan <- sample_plan_4_results$params_input_output$sample_plan 
   parm_summary <- sample_plan_4_results$parameter_summary_all_boots			#str(parm_summary)
   parm_pop <- sample_plan_4_results$pop_parms
    

    p_plan <- ggplot(data=plan, aes(x=binL, y=nsamps )) +
		    geom_bar(stat="identity",fill='#008998') +
		    scale_x_continuous(n.breaks=10) +
		    theme_LH_bar() +
		    annotate("text", x=(0.1*max(plan$binL)), y=(0.9*max(plan$nsamps)), label = paste0("N = ",sum(plan$nsamps))) +
 		    labs(title="", subtitle="", y="Number samples", x="Lower edge of length bin", caption="")
 
#  error in absolute
   
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

    ggsave(paste0(this_dir,"/figures/onaga_4A.png"),grid.arrange(grobs=list(p_plan), ncol=1),width=7.5, height=3, units="in")	
    ggsave(paste0(this_dir,"/figures/onaga_4B.png"),grid.arrange(grobs=list(p_k, p_L0, p_Linf), ncol=3),width=5, height=2, units="in")
	 
  
#  relative error  (obs-true/true)

   parm_summary_raw <- sample_plan_4_results$parameter_summary_all_boots[c(1,2,6),]			#str(parm_summary)
   parm_pop <- as.numeric(sample_plan_4_results$pop_parms[c(1,2,4)])
    
   rel_summary <- data.frame(parm_name = parm_summary_raw$parm_name, 
					lower95 = (parm_summary_raw[,"lower95"] - parm_pop)/parm_pop,
					lower50 = (parm_summary_raw[,"lower50"] - parm_pop)/parm_pop,
					avg = (parm_summary_raw[,"avg"] - parm_pop)/parm_pop,
					upper50 = (parm_summary_raw[,"upper50"] - parm_pop)/parm_pop,   
					upper95 = (parm_summary_raw[,"upper95"] - parm_pop)/parm_pop)

   p_rel <- ggplot(data=rel_summary, aes(x=parm_name)) +
		geom_boxplot(aes(ymin = lower95, lower = lower50, middle = avg, 
			upper = upper50, ymax = upper95), stat = "identity", fill='#E8E8E8') +
 		theme_LH_bar() +
		geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
		labs(title="", subtitle="", y="relative error", x="", caption="")

   ggsave(paste0(this_dir,"/figures/onaga_4B_rel.png"),p_rel,width=5, height=2, units="in")
	



































####