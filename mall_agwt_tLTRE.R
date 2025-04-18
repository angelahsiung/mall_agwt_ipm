library(matrixStats)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)

# Load model output
mall.ipm <-readRDS("mall_ipm_2005-2020_output.rda")

agwt.ipm <- readRDS("agwt_ipm_1992_2020_output.rda")


###################################################
####### contribution to variance of lambda  #######
###################################################
# Code adapted from Koons et al. (2017). Citation below.
## Koons, D. N., T. W. Arnold, and M. Schaub. 2017. “Understanding the Demographic Drivers of Realized Population Growth Rates.”  Ecological Applications 27: 2102–2115.

# Turn tLTRE into a function
tltre.fun <- function(ipm){
  
# # Calculate proportional population sizes
n.years <- ncol(ipm$sims.list$SH.am)
samples <- nrow(ipm$sims.list$SH.am)
draws <- ipm$sims.list # Dig out MCMC samples

SH.am <- ipm$sims.list$SH.am
SN.am <- ipm$sims.list$SN.am
SH.af <- ipm$sims.list$SH.af
SN.af <- ipm$sims.list$SN.af
SH.jm <- ipm$sims.list$SH.jm
SN.jm <- ipm$sims.list$SN.jm
SH.jf <- ipm$sims.list$SH.jf
SN.jf <- ipm$sims.list$SN.jf
R <- ipm$sims.list$R
R <- R[, 1:28]
N_HY_mal <- ipm$sims.list$N_HY_mal
N_AHY_mal <- ipm$sims.list$N_AHY_mal
N_HY_fem <- ipm$sims.list$N_HY_fem
N_AHY_fem <- ipm$sims.list$N_AHY_fem

lam <- ipm$sims.list$lambda

tempvar_lam <- rowVars(lam)
mean(tempvar_lam)
quantile(tempvar_lam,0.05)
quantile(tempvar_lam,0.95)

lambda <- expression(((n.hy.fem + n.ahy.fem)*R*SH.jm*(SN.jm^0.25) + (n.hy.fem + n.ahy.fem)*R*SH.jf*(SN.jf^0.25) + (n.hy.mal + n.ahy.mal)*SH.am*SN.am + (n.hy.fem + n.ahy.fem)*SH.af*SN.af)/(n.hy.fem + n.ahy.fem + n.hy.mal + n.ahy.mal))

D(lambda, "R")

# Calculate age-sex class-specific proportions of abundances at each time step and for each of the saved MCMC samples. 
n.hy.mal <- draws$N_HY_mal[, 1:n.years] / draws$N_tot[, 1:n.years]
n.ahy.mal <- draws$N_AHY_mal[, 1:n.years] / draws$N_tot[, 1:n.years]
n.hy.fem <- draws$N_HY_fem[, 1:n.years] / draws$N_tot[, 1:n.years]
n.ahy.fem <- draws$N_AHY_fem[, 1:n.years] / draws$N_tot[, 1:n.years]

# Calculate the transient sensitivities for each demographic parameter, evaluated at temporal means of each parameter. 
sens_SH_am <- sens_SN_am <- sens_SH_af <-sens_SN_af <- sens_SH_jm <- sens_SN_jm <- sens_SH_jf <- sens_SN_jf <- sens_R <-sens_n_hy_mal <-sens_n_ahy_mal <-sens_n_hy_fem <-sens_n_ahy_fem <- matrix(0,samples,1)

mu <- list(SH.am = rowMeans(SH.am),
           SN.am = rowMeans(SN.am),
           SH.af = rowMeans(SH.af),
           SN.af = rowMeans(SN.af),
           SH.jm = rowMeans(SH.jm),
           SN.jm = rowMeans(SN.jm),
           SH.jf = rowMeans(SH.jf),
           SN.jf = rowMeans(SN.jf),
           R = rowMeans(R),
           n.hy.mal = rowMeans(n.hy.mal),
           n.ahy.mal = rowMeans(n.ahy.mal),
           n.hy.fem = rowMeans(n.hy.fem),
           n.ahy.fem = rowMeans(n.ahy.fem))

sens_SH_am <- eval(D(lambda, "SH.am"), envir=mu)
sens_SN_am <- eval(D(lambda, "SN.am"), envir=mu)
sens_SH_af <- eval(D(lambda, "SH.af"), envir=mu)
sens_SN_af <- eval(D(lambda, "SN.af"), envir=mu)
sens_SH_jm <- eval(D(lambda, "SH.jm"), envir=mu)
sens_SN_jm <- eval(D(lambda, "SN.jm"), envir=mu)
sens_SH_jf <- eval(D(lambda, "SH.jf"), envir=mu)
sens_SN_jf <- eval(D(lambda, "SN.jf"), envir=mu)
sens_R <- eval(D(lambda, "R"), envir=mu)
sens_N_HY_fem <- eval(D(lambda, "n.hy.fem"), envir=mu)
sens_N_AHY_fem <- eval(D(lambda, "n.ahy.fem"), envir=mu)
sens_N_HY_mal <- eval(D(lambda, "n.hy.mal"), envir=mu)
sens_N_AHY_mal<- eval(D(lambda, "n.ahy.mal"), envir=mu)


# Calculate the contributions of temporal process variation and covariations in the demographic parameters to temporal variation in the realized population growth rates.

cont_SH_am <- cont_SN_am <- cont_SH_af <-cont_SN_af <- cont_SH_jm <- cont_SN_jm <- cont_SH_jf <- cont_SN_jf <- cont_R <-cont_N_HY_mal <-cont_N_AHY_mal <-cont_N_HY_fem <-cont_N_AHY_fem <- matrix(0,samples,1)

for (j in 1:samples){
  dp_stoch <- cbind(SH.am[j,],SN.am[j,],SH.af[j,],SN.af[j,], SH.jm[j,], SN.jm[j,], SH.jf[j,], SN.jf[j,], R[j,],n.hy.fem[j,],n.ahy.fem[j,],n.hy.mal[j,], n.ahy.mal[j,])

  dp_varcov <- var(dp_stoch)
  sensvec <- c(sens_SH_am[j], sens_SN_am[j], sens_SH_af[j], sens_SN_af[j], sens_SH_jm[j], sens_SN_jm [j], sens_SH_jf[j], sens_SN_jf[j], sens_R[j], sens_N_HY_fem[j], sens_N_AHY_fem[j], sens_N_HY_mal[j], sens_N_AHY_mal[j])
  contmatrix <- matrix(0,13,13)
   for (k in 1:13){
     for (m in 1:13){
      contmatrix <- dp_varcov*sensvec[k]*sensvec[m]
     }
   }
  contributions <- rowSums(contmatrix)
  cont_SH_am[j] <- contributions[1]
  cont_SN_am[j] <- contributions[2]
  cont_SH_af[j] <- contributions[3]
  cont_SN_af[j] <- contributions[4]
  cont_SH_jm[j] <- contributions[5]
  cont_SN_jm[j] <- contributions[6]
  cont_SH_jf[j] <- contributions[7]
  cont_SN_jf[j] <- contributions[8]
  cont_R[j] <- contributions[9]
  cont_N_HY_fem[j] <- contributions[10]
  cont_N_AHY_fem[j] <- contributions[11]
  cont_N_HY_mal[j] <- contributions[12]
  cont_N_AHY_mal[j] <- contributions[13]
}

# Calculate the mean and Bayesian credible intervals from the  posterior distributions of the LTRE contributions.
totalcont <- cont_SH_am + cont_SN_am + cont_SH_af + cont_SN_af + cont_SH_jm + cont_SN_jm + cont_SH_jf + cont_SN_jf + cont_R + cont_N_HY_mal + cont_N_AHY_mal + cont_N_HY_fem + cont_N_AHY_fem
mean(totalcont)
quantile(totalcont,0.05)
quantile(totalcont,0.95)

cont.means <- c(mean(cont_SH_am), mean(cont_SN_am), mean(cont_SH_af), mean(cont_SN_af), mean(cont_SH_jm), mean(cont_SN_jm), mean(cont_SH_jf), mean(cont_SN_jf), mean(cont_R), mean(cont_N_HY_fem), mean(cont_N_AHY_fem), mean(cont_N_HY_mal), mean(cont_N_AHY_mal))

cont.lower<- c(quantile(cont_SH_am, 0.05), quantile(cont_SN_am, 0.05), quantile(cont_SH_af, 0.05), quantile(cont_SN_af, 0.05), quantile(cont_SH_jm, 0.05), quantile(cont_SN_jm, 0.05), quantile(cont_SH_jf, 0.05), quantile(cont_SN_jf, 0.05), quantile(cont_R, 0.05), quantile(cont_N_HY_fem, 0.05), quantile(cont_N_AHY_fem, 0.05), quantile(cont_N_HY_mal, 0.05), quantile(cont_N_AHY_mal, 0.05))

cont.upper <- c(quantile(cont_SH_am, 0.95), quantile(cont_SN_am, 0.95), quantile(cont_SH_af, 0.95), quantile(cont_SN_af, 0.95), quantile(cont_SH_jm, 0.95), quantile(cont_SN_jm, 0.95), quantile(cont_SH_jf, 0.95), quantile(cont_SN_jf, 0.95), quantile(cont_R, 0.95), quantile(cont_N_HY_fem, 0.95), quantile(cont_N_AHY_fem, 0.95), quantile(cont_N_HY_mal, 0.95), quantile(cont_N_AHY_mal, 0.95))

cont.table <- as.data.frame(cbind(Parameter = c("S_FW_AM", "S_SS_AM", "S_FW_AF", "S_SS_AF", "S_FW_JM", "S_SS_JM", "S_FW_JF", "S_SS_JF", "Productivity", "Proportion_JF", "Proportion_AF", "Proportion_JM", "Proportion_AM"), Contribution = as.numeric(cont.means), Lower = as.numeric(cont.lower), Upper = as.numeric(cont.upper)))

cont.table <- cont.table %>% 
                mutate_at(c('Contribution', "Lower", "Upper"), as.numeric) 

return(cont.table)
}

#Run tltre function for both species

mall.cont.table <- tltre.fun(mall.ipm)
agwt.cont.table <- tltre.fun(agwt.ipm)

# Custom plotting theme
theme_tltre <- function(){ theme(
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 7),
  axis.text.x = element_text(angle = 90, vjust = 0.5),
  plot.title = element_text(size=10),
  panel.background = element_rect(fill = "white", colour = "white",
                                  linewidth = 1, linetype = "solid"),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(colour = "black")
)
}

# Plotting
mall.tltre.plot <- ggplot(mall.cont.table, aes(x = Parameter, y = Contribution)) +
  theme_tltre() +
  geom_col() +
  geom_linerange( aes(x=Parameter, ymin=Lower, ymax=Upper), colour="orange", alpha=0.9, linewidth=1.1)


agwt.tltre.plot <- ggplot(agwt.cont.table, aes(x = Parameter, y = Contribution)) +
  theme_tltre() +
  geom_col() +
  geom_linerange( aes(x=Parameter, ymin=Lower, ymax=Upper), colour="orange", alpha=0.9, linewidth=1.1)

# ggsave("figures/agwt_contribution_plot_2005_2020.jpg", units="cm", width=8, height=7, dpi=600)

mall.tltre.plot | agwt.tltre.plot


# Combine tables then plot using facet_wrap
cont.table <- data.frame(cbind(Species = rep(c("Mallard", "Green-winged Teal"), each = 13), rbind(mall.cont.table, agwt.cont.table)))

ggplot(cont.table, aes(x = Parameter, y = Contribution)) +
theme(strip.text = element_text(size = 9)) +
  theme_tltre() +
  geom_col() +
  geom_linerange( aes(x=Parameter, ymin=Lower, ymax=Upper), colour="orange", alpha=0.9, linewidth=0.5)+
  facet_wrap(~Species, scales = "free_y")


ggsave("figures/contribution_plot_new_both.jpg", units="cm", width=9, height=7, dpi=600)


#########################################################################
####### contribution to differences in lambda in successive years #######
#########################################################################
# AGWT
ipm <- agwt.ipm

lambda <- expression(((N_HY_fem + N_AHY_fem)*R*SH.jm*(SN.jm^0.25) + (N_HY_fem + N_AHY_fem)*R*SH.jf*(SN.jf^0.25) + (N_HY_mal + N_AHY_mal)*SH.am*SN.am + (N_HY_fem + N_AHY_fem)*SH.af*SN.af)/(N_HY_fem+N_AHY_fem+N_HY_mal+N_AHY_mal))

D(lambda, "R")
# Compute differences and means of demographic rates and population structure between successive years
n.years <- length(ipm$mean$SH.am)
n.draws <- ipm$mcmc.info$n.samples # Determine number of MCMC draws
draws <- ipm$sims.list # Dig out MCMC samples
diff <- array(NA, dim=c(n.draws, n.years-1, 13))
dimnames(diff)[[3]] <- c("SH.am", "SN.am", "SH.af", "SN.af", "SH.jm", "SN.jm", "SH.jf", "SN.jf", "R", "N_HY_fem", "N_AHY_fem", "N_HY_mal", "N_AHY_mal")

# Function to compute differences over successive time steps
getDiff <- function(x) x[,2:n.years] - x[,1:(n.years-1)]
diff[,,"SH.am"] <- getDiff(draws$SH.am)
diff[,,"SN.am"] <- getDiff(draws$SN.am)
diff[,,"SH.af"] <- getDiff(draws$SH.af)
diff[,,"SN.af"] <- getDiff(draws$SN.af)
diff[,,"SH.jm"] <- getDiff(draws$SH.jm)
diff[,,"SN.jm"] <- getDiff(draws$SN.jm)
diff[,,"SH.jf"] <- getDiff(draws$SH.jf)
diff[,,"SN.jf"] <- getDiff(draws$SN.jf)
diff[,,"R"] <- getDiff(draws$R)
diff[,,"N_HY_fem"] <- getDiff(draws$N_HY_fem)
diff[,,"N_AHY_fem"] <- getDiff(draws$N_AHY_fem)
diff[,,"N_HY_mal"] <- getDiff(draws$N_HY_mal)
diff[,,"N_AHY_mal"] <- getDiff(draws$N_AHY_mal)
diff.lambda <- getDiff(draws$lambda)
mean.diff.lambda <- apply(diff.lambda, 2, mean)
mean.diff.lambda <- mean.diff.lambda %>% 
  as.data.frame() %>% 
  mutate(Year = c(1993:2019)) %>% 
  pivot_longer(!Year, values_to = "lambda")

ggplot(mean.diff.lambda, aes(x = Year, y = lambda)) +
  geom_bar(stat = "identity")
# Function to compute means over successive time steps, store them in a list
getMn <- function(x) (x[,2:n.years] + x[,1:(n.years-1)]) / 2
means <- list(SH.am=getMn(draws$SH.am), SN.am=getMn(draws$SN.am), SH.af=getMn(draws$SH.af), SN.af = getMn(draws$SN.af), SH.jm=getMn(draws$SH.jm), SN.jm=getMn(draws$SN.jm), SH.jf=getMn(draws$SH.jf), SN.jf=getMn(draws$SN.jf), R=getMn(draws$R), N_HY_fem = getMn(draws$N_HY_fem), N_AHY_fem = getMn(draws$N_AHY_fem), N_HY_mal = getMn(draws$N_HY_mal), N_AHY_mal = getMn(draws$N_AHY_mal))


# Compute sensitivities
senss <- array(NA, dim=c(n.draws, n.years-1, 13))
dimnames(senss)[[3]] <- c("SH.am", "SN.am", "SH.af", "SN.af", "SH.jm", "SN.jm", "SH.jf", "SN.jf", "R", "N_HY_fem", "N_AHY_fem", "N_HY_mal", "N_AHY_mal")
senss[,,"SH.am"] <- eval(D(lambda, "SH.am"), envir=means)
senss[,,"SN.am"] <- eval(D(lambda, "SN.am"), envir=means)
senss[,,"SH.af"] <- eval(D(lambda, "SH.af"), envir=means)
senss[,,"SN.af"] <- eval(D(lambda, "SN.af"), envir=means)
senss[,,"SH.jm"] <- eval(D(lambda, "SH.jm"), envir=means)
senss[,,"SN.jm"] <- eval(D(lambda, "SN.jm"), envir=means)
senss[,,"SH.jf"] <- eval(D(lambda, "SH.jf"), envir=means)
senss[,,"SN.jf"] <- eval(D(lambda, "SN.jf"), envir=means)
senss[,,"R"] <- eval(D(lambda, "R"), envir=means)
senss[,,"N_HY_fem"] <- eval(D(lambda, "N_HY_fem"), envir=means)
senss[,,"N_AHY_fem"] <- eval(D(lambda, "N_AHY_fem"), envir=means)
senss[,,"N_HY_mal"] <- eval(D(lambda, "N_HY_mal"), envir=means)
senss[,,"N_AHY_mal"] <- eval(D(lambda, "N_AHY_mal"), envir=means)

conts <- diff*senss
ci.conts <- apply(conts, c(2,3), cri)
mean.conts <- apply(conts, c(2,3), mean)
mean.conts <- mean.conts %>% 
  as.data.frame() %>% 
  # select(SH.am:R) %>% 
  mutate(Year = c(1993:2019)) %>% 
  pivot_longer(!Year, names_to = "Param", values_to = "Contribution")

ggplot(mean.conts, aes(x = Year, y = Contribution, fill = Param))+
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired")

ci.conts.new <- as.data.frame.table(ci.conts)
ci.conts.new <- data.frame(ci.conts.new) %>% 
                  mutate(quants = rep(c('lower', 'upper'), 351),
                         Year = rep(rep(1993:2019, each = 2), 13)) %>% 
                  rename(Param = Var3) %>% 
                  select(Param, Year, quants, Freq) %>% 
                  # relocate(Parameter, quants) %>% 
                  # pivot_longer(!c(Parameter, quants), names_to = 'year', values_to = 'CI') %>% 
                  pivot_wider(names_from = quants, values_from = Freq)

mean.conts <- full_join(mean.conts, ci.conts.new, by = c("Year" = "Year", "Param" = "Param"))


## Adding plots to illustrate relationship between env covs and contribution 

# Load environmental covariates
envi_cov <- readRDS("data/envi_cov.rda")
# standardize precipitation in southern states
prcp <- (envi_cov$avg_prcp-mean(envi_cov$avg_prcp))/sd(envi_cov$avg_prcp)
# standardize average # of days where max temp below freezing from mid-latitude cities
dx32 <- (envi_cov$avg_dx32_ml-mean(envi_cov$avg_dx32_ml))/sd(envi_cov$avg_dx32_ml)

# May pond count estimates in thousands (1990-2019)
ponds <- c(3508.5, 3200, 3608.9, 3611.7, 5984.8, 6335.4, 7482.2, 7458.2, 4586.9, 6704.3, 3946.9, 4640.4, 2720.0, 5190.1, 3919.6, 5381.2, 6093.9, 7002.7, 4431.4, 6434.0, 6665.0, 8132.2, 5544.0, 6891.7, 7181.2, 6307.7, 5012.5, 6096.0, 5227.4, 4990.3)

# subset ponds to 1992-2019
ponds <- ponds[3:length(ponds)] 
ponds.std <- (ponds-mean(ponds))/sd(ponds) # standardize

# subset precip to 1992-2019
prcp.new <- prcp[2:(length(prcp)-1)] 
# prcp.new <- prcp[15:(length(prcp)-1)] 

prcp.diff <- prcp.new[2:length(prcp.new)]-prcp.new[1:(length(prcp.new)-1)]

mean.conts.R <- mean.conts[mean.conts$Param == "R",]
# mean.conts.R <- cbind(mean.conts.R, lower = ci.conts[1,,'R'], upper = ci.conts[2,,'R'])
nrow(mean.conts.R)
length(prcp.diff)

mean.conts.R <- cbind(mean.conts.R, Precipitation = prcp.diff)

mean.conts.R <- cbind(mean.conts.R, SN.jf = mean.conts[mean.conts$Param == 'SN.jf',])

pond.diff <- ponds.std[2:length(ponds.std)]-ponds.std[1:(length(ponds.std)-1)]
#pond.diff <- ponds.std[15:length(ponds.std)]-ponds.std[14:(length(ponds.std)-1)]

mean.conts.R <- cbind(mean.conts.R, Ponds = pond.diff)

Rcont.precip <- ggplot(mean.conts.R, aes(x = Precipitation, y = Contribution)) +
    # geom_point(aes(x = Precipitation,  y = SN.jf.Contribution)) +
    theme_tltre() +
    geom_pointrange(aes(x = Precipitation, ymin = lower, ymax = upper)) +
  theme(axis.text.x = element_text(angle = 0), axis.text = element_text(size = 16), axis.title = element_text(size = 20)) +
  xlab(paste0('Interannual difference in winter ', '\n', 'precipitation in the south')) +
  ylab('Contribution of reproduction to population growth') +
  annotate("text", x=-3.5, y=1, label= "(b)", size = 13) 

  
Rcont.ponds <- ggplot(mean.conts.R, aes(x = Ponds, y = Contribution)) +
  # geom_point(aes(x = Precipitation,  y = SN.jf.Contribution)) +
  theme_tltre() +
  geom_pointrange(aes(x = Ponds, ymin = lower, ymax = upper)) +
  theme(axis.text.x = element_text(angle = 0), axis.text = element_text(size = 16), axis.title = element_text(size = 20)) +
  xlab("Interannual difference in breeding habitat") +
  ylab('Contribution of reproduction to population growth') +
  annotate("text", x=-2, y=1, label= "(b)", size = 13) 
  


Rcont.precip | Rcont.ponds

ggsave('figures/agwt_Rcont_cov_relationship.jpg', width = 12, height = 8)

 
