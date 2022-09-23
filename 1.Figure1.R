require(pheatmap)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(ggplot2)
library(reshape2)
library(heatmap3)
library(forcats)
library(grid)
library(gridExtra)

library(ggscatter)

wd <- '/Users/marcofondi/Dropbox/PhTAC125/AnTReN/metabolomics/'
setwd(wd)

#import the data
  growth_data_zero <- read.table("growth_curve_0_intra.txt", sep = '\t', header = T, dec = ',')

growth_data_zero$mean <- (growth_data_zero$OD1+growth_data_zero$OD2)/2
growth_data_zero$sd <- apply(growth_data_zero[c(2, 3)],1, sd)


pgrowthzero <- ggplot(data=growth_data_zero, aes(x=T, y=mean)) +
  geom_point(alpha = 1/5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color="#6699CC",se = TRUE)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),  width=.1, position=position_dodge(0.05))+
#  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.3, fill = "grey70") +
  ylab('O.D. (600nm)')+
  xlab('Time (hours)')+
  theme_classic()
pgrowthzero


growth_data_fifteen <- read.table("growth_curve_15_extra.txt", sep = '\t', header = T, dec = ',')

growth_data_fifteen$mean <- (growth_data_fifteen$OD1+growth_data_fifteen$OD2)/2
growth_data_fifteen$sd <- apply(growth_data_fifteen[c(2, 3)],1, sd)

pgrowthfifteen <- ggplot(data=growth_data_fifteen, aes(x=T, y=mean)) +
  geom_point(alpha = 1/5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, color = "#CC6677")+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),  width=.1, position=position_dodge(0.1))+
  #  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.3, fill = "grey70") +
  ylab('O.D. (600nm)')+
  xlab('Time (hours)')+
  theme_classic()
pgrowthfifteen

grid.newpage()
grid.arrange(pgrowthzero,pgrowthfifteen, ncol=2, nrow=1)


#import the GG uptake rate data
GG_uptake_0 <- read.table("GG_uptake_0.txt", sep = '\t', header = T, dec = ',')

GG_uptake_0_plot <- ggplot(data=GG_uptake_0, aes(x=Time, y=mM, factor=Source, color=Source)) +
  geom_point(alpha = 1/5) +
  scale_color_manual(values=c('#332288','#999933'))+
  geom_smooth(method = "lm", formula = y ~ poly(x, 2),se = TRUE)+
  #  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),  width=.1, position=position_dodge(0.05))+
  #  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.3, fill = "grey70") +
  ylab('mM')+
  xlab('Time (hours)')+
  theme_classic()+
  theme(legend.position = "none") 
GG_uptake_0_plot

#import the data
GG_uptake_15 <- read.table("GG_uptake_15.txt", sep = '\t', header = T, dec = ',')

GG_uptake_15_plot <- ggplot(data=GG_uptake_15, aes(x=Time, y=mM, factor=Source, color=Source)) +
  geom_point(alpha = 1/5) +
  #geom_line(aes(y = nM, colour = "Glutamate"))+
  scale_color_manual(values=c('#332288','#999933'))+
  geom_smooth(method = "lm", formula = y ~ poly(x, 2),se = TRUE)+
  #  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),  width=.1, position=position_dodge(0.05))+
  #  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.3, fill = "grey70") +
  ylab('mM')+
  xlab('Time (hours)')+
  theme_classic() +
  theme(legend.position = "none") 

GG_uptake_15_plot

grid.newpage()
grid.arrange(GG_uptake_0_plot,GG_uptake_15_plot, ncol=2, nrow=1)


grid.newpage()
#grid.arrange(GG_uptake_0_plot,GG_uptake_15_plot, ncol=2, nrow=1)
grid.arrange( pgrowthzero,pgrowthfifteen,GG_uptake_0_plot,GG_uptake_15_plot, ncol=2, nrow=2)


# COMPUTE uptake rates for G and G at 0°

OD_to_gL_scaling_factor <- .74;
hours_zero <-  141;
volume_zero <- 1.6 
GG_uptake_0

# %consumed nmol for glutamate (MW = 147.13 g/mol - initial concentration = 5g/l - volume = 1.6 l)
 Glutamate_consumed_moles_zero <- GG_uptake_0$mM[1] -  GG_uptake_0$mM[5]
# %consumed nmol for gluconate (MW = 196.16 g/mol - initial concentration = 5g/l - volume = 1.6 l)
 Gluconate_consumed_moles_zero <- GG_uptake_0$mM[13] -  GG_uptake_0$mM[18]
# %overall biomasss (final OD * scaling factor * total volume (1.6 L)
biomass_zero =growth_data_zero$mean[3]*OD_to_gL_scaling_factor*volume_zero
#                    %Growth rate (ln(OD_final) - ln(OD_init)/hours)
                    mu_zero = (log(growth_data_zero$mean[3]) - log(growth_data_zero$mean[1]))/hours_zero;
                    mu_zero_average = (log(growth_data_zero$mean[6]) - log(growth_data_zero$mean[1]))/244;
                    
                    #                    %Glutamate yield
                    Glutamate_yield_zero <- biomass_zero/Glutamate_consumed_moles_zero
#                    %Gluconate yield
                    Gluconate_yield_zero <- biomass_zero/Gluconate_consumed_moles_zero
#                    %Glutamate uptake rate (Growht rate / Yield)
                    Glutamate_UR_zero <- mu_zero/Glutamate_yield_zero
#                    %Gluconate uptake rate
                    Gluconate_UR_zero <- mu_zero/Gluconate_yield_zero

                                      
                    
# COMPUTE uptake rates for G and G at 15°
                    
OD_to_gL_scaling_factor <- .74
hours_fifteen <- 14
volume_fifteen <- 1.9
GG_uptake_15
                    
                    # %consumed nmol for glutamate (MW = 147.13 g/mol - initial concentration = 5g/l - volume = 1.6 l)
                    Glutamate_consumed_moles_fifteen <- GG_uptake_15$mM[1] -  GG_uptake_15$mM[5]
                    # %consumed nmol for gluconate (MW = 196.16 g/mol - initial concentration = 5g/l - volume = 1.6 l)
                    Gluconate_consumed_moles_fifteen <- GG_uptake_15$mM[13] -  GG_uptake_15$mM[18]
                    # %overall biomasss (final OD * scaling factor * total volume (1.6 L)
                    biomass_fifteen =growth_data_fifteen$mean[3]*OD_to_gL_scaling_factor*volume_fifteen
                    #                    %Growth rate (ln(OD_final) - ln(OD_init)/hours)
                    mu_fifteen = (log(growth_data_fifteen$mean[3]) - log(growth_data_fifteen$mean[1]))/hours_fifteen;
                    mu_fifteen_average = (log(growth_data_fifteen$mean[6]) - log(growth_data_fifteen$mean[1]))/40;
                    
                    #                    %Glutamate yield
                    Glutamate_yield_fifteen <- biomass_fifteen/Glutamate_consumed_moles_fifteen
                    #                    %Gluconate yield
                    Gluconate_yield_fifteen <- biomass_fifteen/Gluconate_consumed_moles_fifteen
                    #                    %Glutamate uptake rate (Growht rate / Yield)
                    Glutamate_UR_fifteen <- mu_fifteen/Glutamate_yield_fifteen
                    #                    %Gluconate uptake rate
                    Gluconate_UR_fifteen <- mu_fifteen/Gluconate_yield_fifteen
                    
                    
                    
                    #average growth rate over the entire gorwth
                    
                    
                    mu_fifteen_1 = (log(growth_data_fifteen$OD1[6]) - log(growth_data_fifteen$OD1[1]))/39;
                    mu_fifteen_2 = (log(growth_data_fifteen$OD2[6]) - log(growth_data_fifteen$OD2[1]))/39;
                    mu15 <- c(mu_fifteen_1, mu_fifteen_2)
                    sd(mu15)
                    mean(mu15)
                    
                    mu_zero_1 = (log(growth_data_zero$OD1[6]) - log(growth_data_zero$OD1[1]))/240;
                    mu_zero_2 = (log(growth_data_zero$OD2[6]) - log(growth_data_zero$OD2[1]))/240;
                    mu0 <- c(mu_zero_1, mu_zero_2)
                    sd(mu0)
                    
                    
                    
                    #exponential phase growth rate over the entire gorwth
                    mu_fifteen_1_exponential = (log(growth_data_fifteen$OD1[3]) - log(growth_data_fifteen$OD1[1]))/hours_fifteen;
                    mu_fifteen_2_exponential = (log(growth_data_fifteen$OD2[3]) - log(growth_data_fifteen$OD2[1]))/hours_fifteen;
                    mu15_exponential <- c(mu_fifteen_1_exponential, mu_fifteen_2_exponential)
                    sd(mu15_exponential)
                    mean(mu15_exponential)
                    
                    
                    mu_zero_1_exponential = (log(growth_data_zero$OD1[3]) - log(growth_data_zero$OD1[1]))/hours_zero;
                    mu_zero_2_exponential = (log(growth_data_zero$OD2[3]) - log(growth_data_zero$OD2[1]))/hours_zero;
                    mu0_exponential <- c(mu_zero_1_exponential, mu_zero_2_exponential)
                    sd(mu0_exponential)
                    mean(mu0_exponential)
                    
                    

                    #plot bar 
                    
                    df2 <- data.frame(temp=rep(c("0", "15"), each=1),
                                      test=rep(c("exp", "model"),each=2),
                                      len=c(mean(mu0_exponential), mean(mu15_exponential), 0.0219,  0.1355),
                                      deviation= c(sd(mu15_exponential), sd(mu0_exponential), 0, 0))
                    
                    
                   p <-  ggplot(data=df2, aes(x=temp, y=len, fill=test)) +
                      geom_bar(stat="identity", color="black", position=position_dodge())+
                      geom_errorbar(aes(ymin=len-deviation, ymax=len+deviation), width=.2,
                                    position=position_dodge(.9))
                    p + scale_fill_brewer(palette="Blues")+theme_minimal()+
                      ylab("Growth rate (h^-1")
                    xlab("")
                    
                    