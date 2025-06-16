#CNMI unfished sampling needs

## S1: Auricilla, catch small, F = 0.5M

onaga <- simulate_population_harvest(Linf = 100, Linf_sd = 2.5, 
                                     M = 0.125, F = 0.0625, Lorenzen = TRUE, 
                                     mincat = 10, catsd = 2.5, maxcat = 200, maxcatsd = 0, 
                                     L0 = 10, L0_sd = 2.5, k = 0.14, k_sd = 0, 
                                     Amax = 55, age_max = 20, N = 100000)
summary(onaga$harvest$length)
str(onaga$harvest$length)
onaga_hist_length_harv <- hist(onaga$harvest$length, breaks = seq(0,120,5),include.lowest=TRUE, right=FALSE,plot=TRUE)

onaga_sample <- LH_sample(sim_output = onaga, n_boots = 1000, samp_size = 400, 
                          sample_type = 'POS', supp_large = FALSE, supp_small = FALSE, constrained = FALSE, 
                          save_bootstraps = TRUE, Amax = 55, age_max = 20, Lbin_width = 5) 




#POS from length data (all CNMI)


