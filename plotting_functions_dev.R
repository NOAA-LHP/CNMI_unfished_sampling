### -----------------------------------------------------------------------------------------
#	June 2025
#
#  Functions for use by RMarkdown to make data reports during the LH / BFish survey.
### -----------------------------------------------------------------------------------------


#---- DEV:

#----- load packages, adjust options, and define working directory

# rm(list=ls())
   if (!require("pacman")) install.packages("pacman")
    pacman::p_load(reshape, dplyr, ggplot2, magrittr, assertthat, gridExtra, readxl, this.path)

# turn off scientific notation in console output
  options(scipen=999)

  this_dir <- this.path::here(.. = 0)


load(paste0(this_dir, "/simulated_pops.RData"))


#  sim_pop_harv_output <- APRU			# names(APRU$parameters$M)




plot_pop_vb_parms <- function(sim_pop_harv_output) {
 
  # prelim calcs
  Jensen_C1 <- sim_pop_harv_output$parameters$Linf/(sim_pop_harv_output$parameters$k^(-1/3))
  plot_jensens <- data.frame(k_range = seq(0.01,0.9,0.005))
  plot_jensens$calc_Linf <- Jensen_C1*(plot_jensens$k_range^(-1/3))
  k_Linf_cor <- round(cor(sim_pop_harv_output$sample_cohort$k, sim_pop_harv_output$sample_cohort$Linf),2)

  plotme <- sample_n(sim_pop_harv_output$sample_cohort, 10000)

    # 4-panel plot in base R
    par(mfrow = c(2,2))
    par(omi=c(0.1,0.1,0.1,0.1)) #set outer margins
    par(mai=c(0.9,0.9, 0.3, 0.1)) #set inner margins

      hist(plotme$L0, main="L0", xlab="cm")
	  abline(v=median(sim_pop_harv_output$sample_cohort$L0), lwd=2, lty=2, col="red")
      hist(plotme$k, main="k", xlab="")
        abline(v=median(sim_pop_harv_output$sample_cohort$k), lwd=2, lty=2, col="red")
      hist(plotme$Linf, main="Linf", xlab="cm")
	  abline(v=median(sim_pop_harv_output$sample_cohort$Linf), lwd=2, lty=2, col="red")
      plot(plotme$k, plotme$Linf, main="cor(k,Linf)", xlab= "k", ylab="Linf")
        lines(plot_jensens$k_range, plot_jensens$calc_Linf, lwd=2, lty=3, col="red")
        mtext(k_Linf_cor,line=-2, side = 3, adj=0.9)

	plot_pop_vb_parms <- recordPlot()
      
  }



plot_fit_pop <- function(sim_pop_harv_output) {


  #  fit the von Bertalanffy parameters to the entire population now
  plot_fit <- nls(length ~ Linf* (1-exp(-(K1*(age-a0)))),
                                   data=sim_pop_harv_output$population, start = list (Linf = 100, K1 = 0.14, a0=0.1))

  plot_parms <- data.frame(Linf = as.numeric(coef(plot_fit)[1]),
				  K1 = as.numeric(coef(plot_fit)[2]),          
				  a0 = as.numeric(coef(plot_fit)[3]), 
				  L0 =  mean(subset(sim_pop_harv_output$population, age == 0)$length))

  plot_vb <- data.frame(age = seq(0, sim_pop_harv_output$parameters$Amax,0.1))
  plot_vb$length <-plot_parms$Linf*(1-exp(-(plot_parms$K1)*(plot_vb$age-plot_parms$a0)))
 

  plot_pop <- sample_n(sim_pop_harv_output$population, 10000)
  plot(plot_pop$age, plot_pop$length, main = "Population length at age", xlab = "years", ylab = "cm")
  lines(plot_vb$age, plot_vb$length, lty=1, lwd=2, col="red")

	plot_fit_pop <- recordPlot()

  }





plot_pop_n_at_age <- function(sim_pop_harv_output) {

  hist(sim_pop_harv_output$population$length,include.lowest=TRUE, right=FALSE,plot=TRUE,
		main = "Population N at length", xlab= "cm", ylab = "count" )

  plot_pop_n_at_age <- recordPlot()

  }



names(sim_pop_harv_output)
str(sim_pop_harv_output$sample_cohort)

try_it <- subset(sim_pop_harv_output$sample_cohort, k < 0.141 & k > 0.139)
      hist(try_it$Linf, xlab="cm")
median(try_it$Linf)

	  abline(v=median(sim_pop_harv_output$sample_cohort$L0), lwd=2, lty=2, col="red")


    # take a look
    plotme <- onaga 
    hist_pop <- hist(plotme$population$length,include.lowest=TRUE, right=FALSE,plot=TRUE)

  Jensen_C1 <- plotme$parameters$Linf/(plotme$parameters$k^(-1/3))
  plot_jensens <- data.frame(k_range = seq(0.01,0.9,0.005))
  plot_jensens$calc_Linf <- Jensen_C1*plot_jensens$k_range^(-1/3)
  k_Linf_cor <- round(cor(plotme$sample_cohort$k, plotme$sample_cohort$Linf),2)


  par(mfrow = c(2,2))

  hist(plotme$sample_cohort$L0, main="L0", xlab="cm")
  hist(plotme$sample_cohort$k, main="k", xlab="")
  hist(plotme$sample_cohort$Linf, main="Linf", xlab="cm")
  plot(plotme$sample_cohort$k, plotme$sample_cohort$Linf, main="cor(k,Linf)", xlab= "k", ylab="Linf")
  lines(plot_jensens$k_range, plot_jensens$calc_Linf, lwd=2, lty=3, col="red")
  mtext(k_Linf_cor,line=-2, side = 3, adj=0.9)



  #  fit the von Bertalanffy parameters to the entire population now
  plot_fit <- nls(length ~ Linf* (1-exp(-(K1*(age-a0)))),
                                   data=plotme$population, start = list (Linf = 100, K1 = 0.14, a0=0.1))



  plot_parms <- data.frame(Linf = as.numeric(coef(plot_fit)[1]),
				  K1 = as.numeric(coef(plot_fit)[2]),          
				  a0 = as.numeric(coef(plot_fit)[3]), 
				  L0 =  mean(subset(plotme$population, age == 0)$length))

  plot_vb <- data.frame(age = seq(0, plotme$parameters$Amax,0.1))
  plot_vb$length <-plot_parms$Linf*(1-exp(-(plot_parms$K1)*(plot_vb$age-plot_parms$a0)))
 

  plot_pop <- sample_n(plotme$population, 10000)
  plot(plot_pop$age, plot_pop$length, main = "Population length at age", xlab = "years", ylab = "cm")
  lines(plot_vb$age, plot_vb$length, lty=1, lwd=2, col="red")





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








###