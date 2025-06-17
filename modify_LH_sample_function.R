#### ---- modify the LH_sample function, make some new ones.




## function 1.
##	given what we think we know about the population, what is the prescribed POS assuming a given total N?

#  -----------------  arguments, typical values
	sim_output <- onaga		# output of the simulate_population_harvest function
   	samp_size <- 300			# sample size for each bootstrap
   	Lbin_width <- 5			# length bin width

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






 
  #  add size bins for the population
  population_true$binL <- cut(population_true[,"length"], breaks=brks, labels=brks[1:(length(brks)-1)],right=F)
  population_true$binL <- as.character(population_true$binL)
  population_true$binL <- as.numeric(population_true$binL)		# head(population_true)	#str(population_true)
  
  #  reset error counter
  count_error <- 0 
  
  # ---- prepare some FOS error warnings outside of the bootstrap loop
  #		this sample won't actually get used in the bootstrap sampling, it is just here so errors don't come out
  #		nsamp times
  if(sample_type == 'FOS') {
    # determine how many fish in each bin in our harvest
    count_harv_samps_per_bin <- summary(population_harv$binL)
    count_harv_samps_per_bin_df <- data.frame(bin_lower = as.numeric(names(count_harv_samps_per_bin)),
                                              bin_counts = as.numeric(count_harv_samps_per_bin))
    bins_w_samps_harv <- nrow(subset(count_harv_samps_per_bin_df, bin_counts != 0))
    min_samps_per_bin <- ceiling(samp_size/bins_w_samps_harv)
    
    # try to take that many samples per bin
    samps_per_bin <- min_samps_per_bin
    FOS <- population_harv %>% dplyr::group_by(binL) %>% sample_n(pmin(n(),samps_per_bin ))
    while(nrow(FOS) < samp_size) {
      samps_per_bin <- samps_per_bin + 1
      FOS <- population_harv %>% dplyr::group_by(binL) %>% sample_n(pmin(n(),samps_per_bin ))
    }
    FOS <- FOS %>% 
      group_by(binL) %>%
      mutate(prop=n()) %>%
      select(age, length, binL, prop, binL_num)
    # coerce to dataframe
    FOS <- as.data.frame(FOS)
    
    # return warnings to alert user if they are asking for a large number of supplemental samples relative to the
    #		number of samples per bin in the unsupplemented design
    if(supp_small == TRUE) {
      if(supp_small_n_per_bin > max(FOS$prop)) {
        warn_message <- paste('warning: max unsupplemented FOS samps per Lbin', max(FOS$prop),
                              'user requested', supp_small_n_per_bin, 'additional samps per Lbin', sep=" ")
        print(warn_message)
      }
    }
    if(supp_large == TRUE) {
      if(supp_large_n_per_bin > max(FOS$prop)) {
        warn_message <- paste('warning: max unsupplemented FOS samps per Lbin', max(FOS$prop),
                              'user requested', supp_large_n_per_bin, 'additional samps per Lbin', sep=" ")
        print(warn_message)
      }
    }
  }
  
  
  for(i in 1:n_boots) {
    
    # print(paste("begin boot",i,sep=" "))
    
    if(sample_type == 'POS') {	# ---------------------  BEGIN POS SAMPLING
      
      
      #head(population_harv_df)		#str(population_harv_df)
      #	head(population_harv)		# str(basic_POS_samps)
      
      POS_build <- population_harv_df[1,]
      for (pp in 1:nrow(basic_POS_samps)) {
        if(basic_POS_samps$basic_POS_n[pp]>0){
          sample_me_binL <- basic_POS_samps$binL_num[pp]
          sample_me_n	<- basic_POS_samps$basic_POS_n[pp]
          sample_me_pop_harv <- subset(population_harv_df, binL_num == sample_me_binL)		#str(sample_me_pop_harv)	
          sample_me <- sample_n(sample_me_pop_harv, sample_me_n, replace=FALSE)
          POS_build <- rbind(POS_build, sample_me)
        }
      }
      #drop row 1 of POS_build
      POS_build <- POS_build[-1,]		#str(POS_build)
      
      #rearrange to make compatible with remainder of sampling...		#str(POS)
      POS <- POS_build %>% select(age, length, binL_num, n)
      colnames(POS) <- c('age','length','binL','n' )
      
      # re-count samps per bin
      POS <- POS %>% 		
        group_by(binL) %>%
        mutate(nbin=n()) %>%
        group_by()
      # coerce back to df
      POS <- as.data.frame(POS)
      # sort by length ascending
      POS <- POS[order(POS$length),]	#nrow(POS)

	need <- unique(POS[,c("binL","nbin")])


            
      #trim POS	head(POS)
      POS <- POS[,c(1:3)]
      
      
            
       # ----  save this sample if randomly chosen
    if(i %in% save_9boots) {
      save_index <- which(i == save_9boots)
      list_some_boot_samps[[save_index]] <- hist(use_sample$length, breaks = seq(0,(max(population_true$length)+Lbin_width),Lbin_width),
                                                 include.lowest=TRUE, right=FALSE,plot=FALSE)
    }
    


    # ------------------   now we have the sample. Solve growth equation.		# use_sample <- list_boot_samps[[i]]
    
    # if t0 constrained to 0
        
    # --- if we estimate L0 freely
    if (constrained == FALSE) {
      
      # check to see if the model will minimize
      trynls <- try(nls(length ~ Linf* (1-exp(-(K1*(age-a0)))),
                        data=use_sample, start = list (Linf = 100, K1 = 0.14, a0=0.1)), silent=T)
      
      # if the model will not minimize (use is.error from assertthat), increase the counter, restart loop 
      #  ---- ERIN MODIFY this step LATER so the loop will automatically do additional i's
      #		to make up for each sample that didn't minimize.		 
      # if the model does minimize, refit the model and proceed with rest of loop
      
      if (is.error(trynls)) { 	
        count_error <- count_error + 1 
      } else
      {
        list_boot_mods[[i]] <- nls(length ~ Linf* (1-exp(-(K1*(age-a0)))),
                                   data=use_sample, start = list (Linf = 100, K1 = 0.14, a0=0.1))
        # t0 was zero here, include estimate t0
        pred.length_use_sample <- coef(list_boot_mods[[i]])[1]*(1 - exp(-coef(list_boot_mods[[i]])[2]*(age - coef(list_boot_mods[[i]])[3])))
        pred_lengths <- data.frame(age=age, length=pred.length_use_sample)
        list_boot_preds[[i]] <- pred_lengths
        CV_boot <- data.frame(age=age, CV_age = as.numeric(NA))
        for (a in 1 :Amax+1) {
          select_fish <- subset(use_sample, age == (a-1))
          CV_a <- sd(select_fish$length)/mean(select_fish$length)
          CV_boot[a,2] <- CV_a
        }
        list_boot_CVs[[i]] <- CV_boot
        
        SD_boot <- data.frame(age=age, SD_age = as.numeric(NA))
        for (a in 1 :Amax+1) {
          select_fish <- subset(use_sample, age == (a-1))
          SD_a <- sd(select_fish$length)
          SD_boot[a,2] <- SD_a
        }
        list_boot_SDs[[i]] <- SD_boot
        
      }  # END if nls error check
      
      # save each bootstrap sample, watch out, this could get big.
      if (save_bootstraps == TRUE) {
        list_boot_samps[[i]] <- use_sample
      }
    }  # --- END if not constrained
    
    # print(i)
    
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
      
      if (constrained == FALSE) {
        parameter_outputs[i,3] <- as.numeric(coef(list_boot_mods[[i]])[3])
      }
      
      if (constrained == TRUE) {
        parameter_outputs[i,3] <- -1
      }
      
      CV_age_boot <- list_boot_CVs[[i]]
      SD_age_boot <- list_boot_SDs[[i]]
      
      #  	 calculate CV_L_age0 and CV_L_ageMax
      #		UPDATE Nov 2 for the manuscript, we assumed SD Length was constant (linear) over ages, so assuming CV is linear here 
      #			is setting us up for poor CV_L0 estimates. Erin modify so SD is assumed linear with length.
      #	 if the slope is 0, nls will throw an error, deal with that, use intercept (avg).
      
      if (SD_L_const == FALSE) {
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
        
      } else {
        
        trynls <- try(nls(SD_age ~ b1*age + b0 ,
                          data=SD_age_boot, start = list (b1=-0.001, b0=0.05)), silent=T)
        if (is.error(trynls)) { 
          mean_sd <- mean(SD_age_boot$SD_age, na.rm = TRUE)
          est_L0 <- parameter_outputs[i,1]*(1-exp(parameter_outputs[i,3]*parameter_outputs[i,2]))
          parameter_outputs[i,4] <- mean_sd/est_L0
          parameter_outputs[i,5] <- mean_sd/parameter_outputs[i,1]
          
        } else {
          
          fit_line_sd <- nls(SD_age ~ b1*age + b0 ,
                             data=SD_age_boot, start = list (b1=-0.001, b0=0.05))
          pred_sd_L0 <- max(0.001,as.numeric(coef(fit_line_sd)[2]))
          est_L0 <- parameter_outputs[i,1]*(1-exp(parameter_outputs[i,3]*parameter_outputs[i,2]))
          parameter_outputs[i,4] <- pred_sd_L0/est_L0
          
          pred_sd_Lmax <- max(0.001,as.numeric(coef(fit_line_sd)[2]) + age_max*as.numeric(coef(fit_line_sd)[1]))
          parameter_outputs[i,5] <- pred_sd_Lmax/parameter_outputs[i,1]
          
        } 
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
  
  
  params_input_output <- list(sim_output_name = substitute(sim_output), 
                              n_boots = n_boots, 
                              samp_size = samp_size, 
                              sample_type = sample_type, 
                              supp_large = supp_large, 
                              supp_large_n_per_bin = supp_large_n_per_bin, 
                              supp_small = supp_small, 
                              supp_small_n_per_bin = supp_small_n_per_bin,
                              supp_min_length = supp_min_length, 
                              constrained = constrained, 
                              t0 = t0,
                              SD_L_const = SD_L_const,
                              save_bootstraps = save_bootstraps, 
                              Amax = Amax , 
                              age_max = age_max, 
                              Lbin_width = Lbin_width,
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
                list_some_boot_samps = list_some_boot_samps))
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
                list_some_boot_samps=list_some_boot_samps))
  }
  
  #  print out runtime
  run_time_calc <- Sys.time() - time_start
  print(run_time_calc)
  
} 		# ---------- END FUNCTION



