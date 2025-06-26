### -----------------------------------------------------------------------------------------
#	June 2025
#	  modify LH_sample functions and add some new functions
#	  these will be used in R markdown to create summary reports
### -----------------------------------------------------------------------------------------





# --------------------------------------------------------------------------------------------------------------------------------------
# ------- READ IN POPULATION SIMULATION FUNCTION  -----------------------
#  updated Sep24: put a floor on L0 draws to avoid negative length fish, should reduce/eliminate warnings()
#	add a true population CV at A0, Amax, Age_max
#	age_max <- 20				# choose age_max (the plus group for pop dynamics) a priori
#					# this is only used for calculating CV_L
#  updated Jun17, 2025: add correlation between L0 and k. new argument Linf_k_cor_TF, new output is Linf,k cor
#				also put floor on k (10%).

  if (1==2) {

	Linf = 100
	Linf_sd = 2.5
	M = 0.125
	F = 0.01
	Lorenzen = 1
	mincat = 10
	catsd = 2.5
	maxcat = 200
	maxcatsd = 0
	L0 = 10
	L0_sd = 2.5
	k = 0.14						
	k_sd = 0.008					# k_sd/k		# try for CV about 5%
	Amax = 55
	age_max = 20
	N = 1000
	Linf_k_cor_TF = TRUE

	test_sph <- simulate_population_harvest(Linf,Linf_sd, M, Lorenzen, F, mincat, catsd, maxcat, maxcatsd, 
                                        L0, L0_sd, k, k_sd, Amax, age_max, N, Linf_k_cor_TF)

	test_sph$parameters


}


simulate_population_harvest <- function(Linf,Linf_sd, M, Lorenzen, F, mincat, catsd, maxcat, maxcatsd, 
                                        L0, L0_sd, k, k_sd, Amax, age_max, N, Linf_k_cor_TF){
  
  # create Amax + 1 cohorts of N Age zero fish, put them all into a list object
  all_cohorts <- list()
  all_harvest <- list()
  #  hr = 1-exp(-F)
  #  mr = 1-exp(-M)
  
  # --------  IF USING LORENZEN LENGTH-DEPENDENT NATURAL MORTALITY 
  # Solve for M1 now
  if (Lorenzen == TRUE | Lorenzen == 1) {
    survship <- data.frame(Ages = seq(0,Amax,1), L_age=9999, M_L_hat=9999, 
                           Mr_L_hat = 9999, N_age=9999)
    # calc A0
    A0 = log(1-L0/Linf)/k
    
    # grow the fish
    for (i in 1:(Amax+1)) {
      survship[i,2] = Linf*(1-exp(-k*(survship[i,1]-A0)))
    }
    
    M1 = 10		# reset M1 to something reasonable
    # set up a function to do the Lorenzen distribution
    sq_resid_M1 = function(M1) {
      for (i in 1:(Amax+1)) {
        survship[i,3] = M1/survship[i,2]
      }
      survship$Mr_L_hat = 1-exp(-1*survship$M_L_hat)
      survship[1,5] = 10000
      for (i in 2:(Amax+1)) {
        survship[i,5] = survship[(i-1),5]*(1-survship[(i-1),4])
      }
      M_overall_hat = -1*(log(survship[(Amax+1),5]/10000)/Amax)
      sq_resid = (M - M_overall_hat)^2
      return(sq_resid)
    }
    
    # run optimize for 1 parmater to fit M1
    fit = optimize(sq_resid_M1, interval=c(0.01,100))
    M1 <- fit$minimum 
  }
  # -----------  END solve Lorenzen M1
  
  # ---- Create Amax+1 Cohorts
  
  for (i in 1:(Amax+1)) {
    cohort <- matrix(nrow=N,ncol=Amax+4)
    
    # assign L0, Linf, k to each fish, calculate A0
    cohort[,3] = pmax(rnorm(N, mean = k, sd = k_sd), rep(0.1*k,N))						# summary(cohort[,3])
    cohort[,4] = pmax(rnorm(N, mean = L0, sd = L0_sd),rep(0.1*L0,N))

    ## ** JUNE 2025 UPDATE
    #  Jensen 1997 Eq. 4 theoretical expected relationship between Linf and K
    #	    Linf = C1*k^(-1/LW_beta)	where LW_beta is the theoretical allometric W-L scalar, assume = 3 here		
    #	 so, if we want them to correlate; mean Linf to sample from is calculated from k and solved value of Jensen_C1

      if(Linf_k_cor_TF == TRUE | Linf_k_cor_TF == 1) {
	  Jensen_C1 <- Linf/(k^(-1/3))
	  cohort[,2] = rnorm(N, mean = (Jensen_C1*(cohort[,3]^(-1/3))), sd = Linf_sd)
 	  }
	if(Linf_k_cor_TF == FALSE) {
    	  cohort[,2] = rnorm(N, mean = Linf, sd = Linf_sd)
	  }
     # check Linf=f(k) and correlation (save for later, will be last i only in parameters list)
     # plot(cohort[,3], cohort[,2]) 	# cor(cohort[,3],cohort[,2])
        k_Linf_cor_value <- cor(cohort[,3],cohort[,2])
    
    cohort[,1]= (log(1-(cohort[,4]/cohort[,2])))/cohort[,3]

    if (i == 1) {
	sample_cohort <- as.data.frame(cohort[,1:4])
	names(sample_cohort)[] <- c("a0","Linf","k","L0")
	}

    #  end update

    # ----- step A. grow the fish, all survive
    for (A in 1:Amax) {
      cohort[,A+4] = cohort[,2]*(1-exp(-(cohort[,3])*(A-cohort[,1])))
    }
    # ----  step B. apply natural mortality - Baranov - depends on fixed or Lorenzen
    
    selex_cohort = pmin((1-pnorm(mincat,cohort[,3:(Amax+3)],catsd)),pnorm(maxcat ,cohort[,3:(Amax+3)],maxcatsd))   
    
    if (Lorenzen == FALSE) {
      mp_cohort <- (M/(F*selex_cohort + M))*(1-exp(-F*selex_cohort - M))
    } else {
      M_L_cohort <- M1 / cohort[,4:(Amax+4)]
      mp_cohort <- (M_L_cohort/(F*selex_cohort + M_L_cohort))*(1-exp(-F*selex_cohort - M_L_cohort))
    }
    # --- Nat mort survival probabilities
    surv_mr = matrix(rbinom(N*(Amax+1), size = 1, prob = 1-(mp_cohort)),nrow=N, ncol=Amax+1)
    
    # all fish used to survive, now we kill them based on length. So N needs to be big.
    
    # modify mr into a cumulative cross-product
    for (A in 1:(Amax)) {
      surv_mr[,A+1]=surv_mr[,A]*surv_mr[,A+1]
    }
    # lengths for fish that survive natural mortality, 0s for those that don't
    mr_surv_lengths = cohort[,4:(Amax+4)]*surv_mr
    
    # ----- step C. apply fishing mortality and collect catch
    # will depend on fishery selectivity at length as well as Lorenzen vs. not
    
    if (Lorenzen == FALSE) {
      hp_cohort <- ((F*selex_cohort)/(F*selex_cohort+M))*(1-exp(-F*selex_cohort-M))
    } else {
      hp_cohort <- ((F*selex_cohort)/(F*selex_cohort+M_L_cohort))*(1-exp(-F*selex_cohort-M_L_cohort))
    }
    
    catch_cohort = matrix(rbinom(N*(Amax+1), size = 1, prob = hp_cohort),nrow=N, ncol=Amax+1)
    
    #  each fish only gets caught once
    # modify catch_cohort into additive vector over ages so once a fish gets "caught" the second time, catch_cohort > 1
    for (A in 1:(Amax)) {
      catch_cohort[,A+1]=catch_cohort[,A]+catch_cohort[,A+1]
    }
    # this is goofy- take another additive product so that the age a fish gets caught first = 1, subsequent ages > 1
    for (A in 1:(Amax)) {
      catch_cohort[,A+1]=catch_cohort[,A]+catch_cohort[,A+1]
    }
    # the age a fish gets caught first = 1, all others = 0
    catch_once_cohort <- replace(catch_cohort, catch_cohort > 1, 0)
    
    # lengths for fish that survive natural mortality AND GET CAUGHT, 0s for those that don't
    catch_lengths = mr_surv_lengths*catch_once_cohort
    #  save it in the all_harvest list
    all_harvest[[i]] <- catch_lengths
    
    # apply fishing mortality by taking the fish that got caught (identified above) out of the population of fish that survived natural mortality
    hr_surv <- 1-catch_once_cohort
    # turn this into a cumulative cross-product
    for (A in 1:(Amax)) {
      hr_surv[,A+1]=hr_surv[,A]*hr_surv[,A+1]
    }
    # lengths for fish that survive, 0s for those that don't
    cohort_lengths = mr_surv_lengths*hr_surv
    #  save it
    all_cohorts[[i]] <- cohort_lengths
  }
  
  # now pull one age from each cohort to build the overall population
  # create a population object using Age0
  population <- data.frame(age=0, length = all_cohorts[[1]][,1])
  #drop dead fish
  population <- subset(population , length != 0)
  
  for (i in 1:(Amax)) {
    each_age <- data.frame(age=i, length = all_cohorts[[(i+1)]][,(i+1)])
    #drop dead fish
    each_age <- subset(each_age , length != 0)
    #save
    population <- rbind(population, each_age)
  }
  #  we have a population
  
  #  repeat for harvest 
  harvest <- data.frame(age=0, length = all_harvest[[1]][,1])
  #drop dead fish
  harvest <- subset(harvest , length != 0)
  for (i in 1:(Amax)) {
    each_h_age <- data.frame(age=i, length = all_harvest[[(i+1)]][,(i+1)])
    #drop dead fish
    each_h_age <- subset(each_h_age , length != 0)
    #save
    harvest <- rbind(harvest, each_h_age)
  }
  #  we have a harvest
  
  # add a warning message if L0_sd was large compared to L0 such that we have negative length fish
  if (min(population$length)<0) { print("Warning: negative length fish, L0_sd too large compared to L0") }
  
  # make a quick average M at age object
  A0 = log(1-L0/Linf)/k
  Avg_age <- data.frame(Ages = seq(0,Amax,1), L_age=9999, M_age = 9999, Selex = 9999)
  # grow the fish
  for (i in 1:(Amax+1)) {
    Avg_age[i,2] = Linf*(1-exp(-k*(Avg_age[i,1]-A0)))
  }
  
  if (Lorenzen == FALSE) {
    Avg_age$M_age = M
  } else {
    Avg_age$M_age = M1 / Avg_age$L_age
  }
  
  Avg_age$Selex <- pmin((1-pnorm(mincat,Avg_age$L_age,catsd)),pnorm(maxcat,Avg_age$L_age,maxcatsd))
  
  if (Lorenzen == TRUE | Lorenzen == 1) {
    M1 = M1
  } else {
    M1 = 'NA'
  }
  
  # calculate CV length at age
  
  all_big_fish_pop <- subset(population, age >= age_max)
  all_big_fish_harv <- subset(harvest, age >= age_max)
  all_big_fish <- rbind(all_big_fish_pop, all_big_fish_harv)
  CV_L_age_max <- sd(all_big_fish$length)/mean(all_big_fish$length)
  
  all_0_fish_pop <- subset(population, age == 0)
  all_0_fish_harv <- subset(harvest, age == 0)
  all_0_fish <- rbind(all_0_fish_pop, all_0_fish_harv)
  CV_L_0 <- sd(all_0_fish$length)/mean(all_0_fish$length)		# in theory should be very close to L0_sd / L0
  
  params <- list(Linf = Linf,
                 Linf_sd = Linf_sd,
                 Lorenzen = Lorenzen,
                 M = M,
                 M1 = M1,
                 F = F, 
                 mincat = mincat , 
                 catsd = catsd , 
                 maxcat = maxcat , 
                 maxcatsd = maxcatsd , 
                 L0 = L0 , 
                 L0_sd= L0_sd, 
                 k = k,
                 k_sd = k_sd, 
                 Amax = Amax,
                 age_max = age_max,
                 CV_L_age_max =  CV_L_age_max,
                 CV_L_0 =  CV_L_0,
                 N = N,
		     Linf_k_cor_TF = Linf_k_cor_TF,
		     k_Linf_cor_value = k_Linf_cor_value)
  
  return(list(population = population, harvest = harvest, Avg_age = Avg_age, parameters = params, sample_cohort = sample_cohort))
  
}	#  --------------------------------  end function


# --------------------------------------------------------------------------------------------------------------------------------------



## function 2.
##	given what we think we know about the population, what is the prescribed POS assuming a given total N?

if (1==2) {
#  -----------------  arguments, typical values
	sim_output <- onaga		# output of the simulate_population_harvest function
   	samp_size <- 300			# sample size for each bootstrap
   	Lbin_width <- 5			# length bin width
}



find_POS <- function(sim_output, samp_size, Lbin_width) {
  
   # for an unfished population and assume all lengths are fully selected, so "harvest" (e.g., sampling) will be 
   #		proportional to population
   population_harv <- population_true <- sim_output$population
  
  #  bin harvest and population lengths		#head(population_harv)
  brks <- seq(0,1.1*max(population_true$length),Lbin_width) 
  age <- seq(0,sim_output$parameters$Amax,1)
  population_harv$binL <- cut(population_harv[,"length"], breaks=brks, labels=brks[1:(length(brks)-1)],right=F)
  
  population_harv <- population_harv %>% 		
    group_by(binL) %>%
    mutate(prop=n()) %>%
    group_by() %>%
    mutate(n=sum(n())) %>%
    mutate(propL=prop/n)					#str(population_harv)	head(population_harv) 
  
  population_harv <- mutate(population_harv, binL_num = as.character(binL))
  population_harv$binL_num <- as.numeric(population_harv$binL_num)		#population_harv[1:500,]
  
  population_harv_df <- as.data.frame(population_harv)
  
    basic_POS_props <- unique(population_harv_df[c("binL","propL")])		#sum(basic_POS_props$prop1)
    basic_POS_props <- basic_POS_props[order(basic_POS_props$binL),]		#str(basic_POS_samps)
    basic_POS_samps <- mutate(basic_POS_props, basic_POS_n = round((samp_size*propL),0), 
                              binL_num = as.numeric(as.character(binL)))	#basic_POS_samps[9,3] <- 40
    # 
    if (sum(basic_POS_samps$basic_POS_n) > samp_size) {
      while(sum(basic_POS_samps$basic_POS_n) > samp_size) {
        most_index <- which(basic_POS_samps$basic_POS_n == max(basic_POS_samps$basic_POS_n))[1]
        basic_POS_samps[most_index,3] <- basic_POS_samps[most_index,3] - 1
      }
    }
    if (sum(basic_POS_samps$basic_POS_n) < samp_size) {
      while(sum(basic_POS_samps$basic_POS_n) < samp_size) {
        most_index <- which(basic_POS_samps$basic_POS_n == max(basic_POS_samps$basic_POS_n))[1]
        basic_POS_samps[most_index,3] <- basic_POS_samps[most_index,3] + 1
      }
    }
 
  
    POS_plan <- basic_POS_samps[,c("binL","basic_POS_n")]				# str(POS_plan)
	names(POS_plan)[] <- c("binL","plan")
    POS_plan$binL <- as.numeric(as.character(POS_plan$binL))
    rownames(POS_plan) <- NULL

    return(POS_plan)

  }


# --------------------------------------------------------------------------------------------------------------------------------------



## function 3.
##	Fit vonB with a given total sample size and length comp- what is the error?

# args
if (1==2) {
#  define the sample plan

  sample_plan <- data.frame(binL = summary$binL, nsamps = summary$have)		# str(sample_plan)
  sim_output <- onaga
  n_boots <- 7
  age_max <- sim_output$parameters$age_max		# define age max which is for CV_L_old
  save_bootstraps <- TRUE

#  try it

  try_it <- fit_plan(sample_plan, sim_output, n_boots, age_max, save_bootstraps) 

# names(try_it)			#try_it$parameter_outputs

}


#  -----------------  BEGIN READ IN FUNCTION

fit_plan <- function(sample_plan, sim_output, n_boots, age_max, save_bootstraps) {
  
  #  Preliminaries:
  #  Get system time
  time_start <- Sys.time()
  print(time_start)
  count_error <- 0
  
  #  define population_true and population_harvest
  population_true <- sim_output$population				# str(population_true)	
  population_harv <- sim_output$population			
  
  #  fit the von Bertalanffy parameters to the entire population now
  pop_fit <- nls(length ~ Linf* (1-exp(-(K1*(age-a0)))),
                                   data=population_true, start = list (Linf = 100, K1 = 0.14, a0=0.1))

  pop_parms <- data.frame(Linf = as.numeric(coef(pop_fit)[1]),
				  K1 = as.numeric(coef(pop_fit)[2]),          
				  a0 = as.numeric(coef(pop_fit)[3]), 
				  L0 =  mean(subset(population_true, age == 0)$length),
				  sd_L0 = sd(subset(population_true, age == 0)$length),
				  L_age_max = mean(subset(population_true, age == age_max)$length),
				  sd_L_age_max = sd(subset(population_true, age == 0)$length))

  #  make empty lists to hold each of the things we are interested in:
  list_boot_mods <- list()
  list_boot_preds <- list()
  list_boot_CVs <- list()
  list_boot_SDs <- list()
  list_boot_samps <- list()
  list_some_boot_samps <- list()
  
  #  bin population lengths (get bins from our plan)
  brks <- plan$binL
  age <- seq(0,sim_output$parameters$Amax,1)
  population_harv$binL <- cut(population_harv[,"length"], breaks=brks, labels=brks[1:(length(brks)-1)],right=F)
  
  # add back in a numeric binL column, return to df
  population_harv <- mutate(population_harv, binL_num = as.character(binL))
  population_harv$binL_num <- as.numeric(population_harv$binL_num)		#population_harv[1:500,]
  
  population_harv_df <- as.data.frame(population_harv)			# head(population_harv_df)	# str(population_harv_df)
  
  #  sample according to our plan
  #  BEGIN BOOTSTRAP   --------------------------
  
  for(i in 1:n_boots) {
    
	# take row 1 of the df to get the column names and formats
	sample_build <- population_harv_df[1,]
      for (pp in 1:nrow(sample_plan)) {
          sample_me_binL <- sample_plan$binL[pp]
          sample_me_n	<- sample_plan$nsamps[pp]
          sample_me_pop_harv <- subset(population_harv_df, binL_num == sample_me_binL)		#View(sample_build)	
          sample_me <- sample_n(sample_me_pop_harv, sample_me_n, replace=FALSE)
          sample_build <- rbind(sample_build, sample_me)
          }
      
      #drop row 1 of sample_build
      sample_boot <- sample_build[-1,]		#str(sample_boot)
    
      if (save_bootstraps == TRUE) { 
	  list_boot_samps[[i]] <- sample_boot
	  }

    # ------------------   now we have the sample. Solve growth equation.		# use_sample <- list_boot_samps[[i]]
    
    # check to see if the model will minimize
      trynls <- try(nls(length ~ Linf* (1-exp(-(K1*(age-a0)))),
                        data=sample_boot, start = list (Linf = 100, K1 = 0.14, a0=0.1)), silent=T)
      
      # if the model will not minimize (use is.error from assertthat), increase the counter, restart loop 
      #  ---- ERIN MODIFY this step LATER so the loop will automatically do additional i's
      #		to make up for each sample that didn't minimize.		 
      # if the model does minimize, refit the model and proceed with rest of loop
      
      if (is.error(trynls)) { 	
        count_error <- count_error + 1 
      } else
      {
        list_boot_mods[[i]] <- nls(length ~ Linf* (1-exp(-(K1*(age-a0)))),
                                   data=sample_boot, start = list (Linf = 100, K1 = 0.14, a0=0.1))
       
        pred.length_use_sample <- coef(list_boot_mods[[i]])[1]*(1 - exp(-coef(list_boot_mods[[i]])[2]*(age - coef(list_boot_mods[[i]])[3])))
        pred_lengths <- data.frame(age=age, length=pred.length_use_sample)
        list_boot_preds[[i]] <- pred_lengths
        CV_boot <- data.frame(age=age, CV_age = as.numeric(NA))
        for (a in 1 :(length(age)+1)) {
          select_fish <- subset(sample_boot, age == (a-1))
          CV_a <- sd(select_fish$length)/mean(select_fish$length)
          CV_boot[a,2] <- CV_a
        }
        list_boot_CVs[[i]] <- CV_boot
        
        SD_boot <- data.frame(age=age, SD_age = as.numeric(NA))
        for (a in 1 :(length(age)+1)) {
          select_fish <- subset(sample_boot, age == (a-1))
          SD_a <- sd(select_fish$length)
          SD_boot[a,2] <- SD_a
        }
        list_boot_SDs[[i]] <- SD_boot
        
      }  
          
  }	# ------ END boot i
  
  
  #  make a warning message
  if (count_error > 0) {
    warn_message <- paste('warning: bootstraps convergence fail ', count_error, sep="")
    print(warn_message)
  }
  
  parameter_outputs <- data.frame(Linf = numeric(), K1 = numeric(), a0 = numeric(), CV_L_a0 = numeric(), CV_L_age_max = numeric())
  
  for (i in 1: length(list_boot_mods)) {
    if (is.null(list_boot_mods[[i]])) {
      # print('von bert nls error')
    } else {
      parameter_outputs[i,1] <- as.numeric(coef(list_boot_mods[[i]])[1])
      parameter_outputs[i,2] <- as.numeric(coef(list_boot_mods[[i]])[2])
      parameter_outputs[i,3] <- as.numeric(coef(list_boot_mods[[i]])[3]) 
      CV_age_boot <- list_boot_CVs[[i]]
      SD_age_boot <- list_boot_SDs[[i]]
      
      #  	 calculate CV_L_age0 and CV_L_ageMax
      #		UPDATE Nov 2 for the manuscript, we assumed SD Length was constant (linear) over ages, so assuming CV is linear here 
      #			is setting us up for poor CV_L0 estimates. Erin modify so SD is assumed linear with length.
      #	 if the slope is 0, nls will throw an error, deal with that, use intercept (avg).
      
      
        trynls <- try(nls(CV_age ~ b1*age + b0 ,
                          data=CV_age_boot, start = list (b1=-0.001, b0=0.05)), silent=T)
        if (is.error(trynls)) { 
          parameter_outputs[i,4] <- mean(CV_age_boot$CV_age, na.rm = TRUE)
          parameter_outputs[i,5] <- mean(CV_age_boot$CV_age, na.rm = TRUE)
        } else {
          fit_line_cv <- nls(CV_age ~ b1*age + b0 ,
                             data=CV_age_boot, start = list (b1=-0.001, b0=0.05))
          parameter_outputs[i,4] <- max(0.001,as.numeric(coef(fit_line_cv)[2]))
          parameter_outputs[i,5] <- max(0.001,as.numeric(coef(fit_line_cv)[2]) + age_max*as.numeric(coef(fit_line_cv)[1]))
        } 
     }   
  }	#  --- END create parameter_outputs

  # back-calculate L0 for easy comparison
  parameter_outputs <- mutate(parameter_outputs, L0 = Linf*(1-exp(a0*K1)))
  
  #    Extract out summaries of all bootstrap runs, put in a df
  
  extract_summary <- data.frame(parm_name = character(), lower95 = numeric(), lower50 = numeric(), 
                                avg = numeric(), upper50 = numeric(), upper95 = numeric(), stringsAsFactors=FALSE)
  
  for (j in 1: ncol(parameter_outputs)) {
    extract_summary[j,1] <- colnames(parameter_outputs)[j]
    extract_summary[j,2] <- as.numeric(quantile(parameter_outputs[,j],probs = 0.025, na.rm = TRUE))
    extract_summary[j,3] <- as.numeric(quantile(parameter_outputs[,j],probs = 0.25, na.rm = TRUE))
    extract_summary[j,4] <- mean(parameter_outputs[,j],na.rm=TRUE)
    extract_summary[j,5] <- as.numeric(quantile(parameter_outputs[,j],probs = 0.75, na.rm = TRUE))
    extract_summary[j,6] <- as.numeric(quantile(parameter_outputs[,j],probs = 0.975, na.rm = TRUE))
  }
  
  # sample_plan, sim_output, n_boots, age_max
  params_input_output <- list(sim_output_name = substitute(sim_output), 
                              n_boots = n_boots, 
                              sample_plan = sample_plan, 
                              age_max = age_max,
                              boots_nls_fail = count_error)
  
  if (save_bootstraps == TRUE) {
    
    return(list(list_boot_samps = list_boot_samps,
                list_boot_preds = list_boot_preds,
                list_boot_mods = list_boot_mods, 
                list_boot_CVs = list_boot_CVs,
                list_boot_SDs = list_boot_SDs,
                parameter_outputs = parameter_outputs,
                params_input_output = params_input_output,
                simulation_params = sim_output$parameters,
                parameter_summary_all_boots = extract_summary,
                list_some_boot_samps = list_some_boot_samps,
		    pop_parms = pop_parms))
  }
  
  if (save_bootstraps == FALSE) {
    
    return(list(list_boot_samps = 'NA',
                list_boot_preds = list_boot_preds,
                list_boot_mods = list_boot_mods, 
                list_boot_CVs = list_boot_CVs,
                list_boot_SDs = list_boot_SDs,
                parameter_outputs = parameter_outputs,
                params_input_output = params_input_output,
                simulation_params = sim_output$parameters,
                parameter_summary_all_boots = extract_summary,
                list_some_boot_samps=list_some_boot_samps,
		    pop_parms = pop_parms))
  }
  
  #  print out runtime
  run_time_calc <- Sys.time() - time_start
  print(run_time_calc)
  
} 		# ---------- END FUNCTION





#  ggplot theme

theme_LH_bar <- function(){
  font <- 'sans'   #assign font family up front
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      panel.grid.major.y = element_line(colour = "#D8D8D8"), #major y gridlines in gray
      panel.grid.major.x = element_blank(),  
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(0, 12, 0, 0), # margin top, right, bottom, left in pts
      
      
      #text elements
      plot.title = element_text(size=10),             # do not save space for the title
	strip.text = element_text(size=10),
	strip.text.x = element_text(margin = margin(10, 0, 1, 0)),      

      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 14),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 10),               #font size
      
      axis.title.y = element_text(angle = 90, margin = margin(l = 4), vjust = 0),
	axis.title.x = element_text(angle = 0, margin = margin(t=5, r=0, b=0, l=0), vjust = 0),
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 10),                #font size
      
      # move labels (numbers) closer to axis
      axis.text.y = element_text(margin=margin(t = 0, r = 0, b = 0, l= 10, unit = "pt")),
      axis.text.x = element_text(margin=margin(t = -18, r = 0, b = 0, l= 0, unit = "pt")),
     
	# axis.ticks.x = element_line(color = "black"),  # show x-axis ticks
	# axis.ticks.length = unit(1, "pt"),             # control length of ticks

     legend.position = "none",
#      legend.text = element_blank()
 
#      legend.title = element_blank(), # put legend right, top
#      legend.justification = "top",
      legend.text = element_text(size = 8)
    )
}










#### --------------


