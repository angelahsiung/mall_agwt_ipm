##################################################
### The purpose of this script is to load the data
### needed to run the IPM for mallard and 
### green-winged teal, run the model and plot the output
### Last update date: 3/26/2025
###################################################

# Prepare packages
list.of.packages <- c("tidyverse", "rjags","jagsUI", "coda","parallel","doParallel", "foreach", "here", "MCMCvis", "ggridges", "jagshelper", "GGally", "patchwork", "gridExtra", "ggpubr", "rphylopic", "ggdist", "todor")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
lapply(list.of.packages, require, character.only = TRUE)

#------ Load and prepare data -----------
## Band-recovery marrays

# # Mallard
MALL.marray <- readRDS('data/MALL_marray.rda')
MALL.release <- readRDS('data/MALL_release.rda')
# dim(MALL.marray)
# dim(MALL.release)

# # # subset to years 2000-2020 for mallards
MALL.marray <- MALL.marray[10:nrow(MALL.marray), 10:ncol(MALL.marray),,]
MALL.release <- MALL.release[10:nrow(MALL.release),,]

# Green-winged Teal
# AGWT.release <- readRDS(file = "data/AGWT_release.rda")
# AGWT.marray <- readRDS(file = "data/AGWT_marray.rda")
# dim(AGWT.marray)

# Subset to 2000-2020 for green-winged teal
# AGWT.release <- AGWT.release[10:nrow(AGWT.release),,]
# AGWT.marray <- AGWT.marray[10:nrow(AGWT.marray), 10:ncol(AGWT.marray),,]

# Bpop estimates (divided by 10,000)
y_esa <- readRDS(file = "data/MALL_ESA_Bpop.rda") # Eastern survey area
y_tsa <- readRDS(file = "data/MALL_TSA_Bpop.rda") # Traditional survey area

# y_esa <- readRDS(file = "data/AGWT_ESA_Bpop.rda")
# y_tsa <- readRDS(file = "data/AGWT_TSA_Bpop.rda")

# replace NAs with 0 for y_esa for summing
y_esa[is.na(y_esa)] <- 0

# Add estimates from both survey areas
y_tot <- y_esa + y_tsa

# Subset for year 2000-2020 
y_tot <- y_tot[10:(length(y_tot)-1)]

# load wing data (females)
jf.wing <- readRDS(file = "data/MALL_juv_female_wing.rda")
female.wing <- readRDS(file = "data/MALL_all_female_wing.rda")

# jf.wing <- readRDS(file = "data/AGWT_juv_female_wing.rda")
# female.wing <- readRDS(file = "data/AGWT_all_female_wing.rda")

# Subset to 2000-2020
jf.wing <- jf.wing[10:nrow(jf.wing),]
female.wing <- female.wing[10:nrow(female.wing),]


# Load environmental covariates
envi_cov <- readRDS("data/envi_cov.rda")

# standardize precipitation in southern states
prcp <- (envi_cov$avg_prcp-mean(envi_cov$avg_prcp))/sd(envi_cov$avg_prcp)

# standardize average # of days where max temp below freezing from mid-latitude cities
dx32 <- (envi_cov$avg_dx32_ml-mean(envi_cov$avg_dx32_ml))/sd(envi_cov$avg_dx32_ml)

# May pond count estimates in thousands (1992-2019)
ponds <- c(3608.9, 3611.7, 5984.8, 6335.4, 7482.2, 7458.2, 4586.9, 6704.3, 3946.9, 4640.4, 2720.0, 5190.1, 3919.6, 5381.2, 6093.9, 7002.7, 4431.4, 6434.0, 6665.0, 8132.2, 5544.0, 6891.7, 7181.2, 6307.7, 5012.5, 6096.0, 5227.4, 4990.3)

# subset ponds to 2000-2019
ponds <- ponds[9:length(ponds)]
ponds.std <- (ponds-mean(ponds))/sd(ponds) # standardize

#------ Running IPM ----------
# Bundle data
nyrs <- dim(MALL.marray)[1]
nClass <- dim(MALL.marray)[4]
# nyrs <- dim(AGWT.marray)[1]
# nClass <- dim(AGWT.marray)[4]

bugs.data <- list(nyrs=dim(MALL.marray)[1],
                  #pre-hunting band-recoveries (calendar year 2000-2020)
                  recovmat.am = MALL.marray[,,2,3], recovmat.af = MALL.marray[,,2,4], recovmat.jm = MALL.marray[,,2,1], recovmat.jf = MALL.marray[,,2,2]
                  #post-hunting band-recoveries (calendar year 2001-2020)
                  , recovmatP.am = MALL.marray[2:nyrs,2:(nyrs+1),1,3], recovmatP.af = MALL.marray[2:nyrs,2:(nyrs+1),1,4], recovmatP.jm = MALL.marray[2:nyrs,2:(nyrs+1),1,1], recovmatP.jf = MALL.marray[2:nyrs,2:(nyrs+1),1,2]
                  #pre-hunting total bandings
                  , relmat.am = MALL.release[,2,3], relmat.af = MALL.release[,2,4], relmat.jm = MALL.release[,2,1], relmat.jf = MALL.release[,2,2]
                  #post-hunting total bandings
                  , relmatP.am = MALL.release[2:nyrs,1,3], relmatP.af = MALL.release[2:nyrs,1,4], relmatP.jm = MALL.release[2:nyrs,1,1], relmatP.jf = MALL.release[2:nyrs,1,2]
                  # wing data
                  , W.jv = jf.wing$Number, W.tot = female.wing$sum
                  #Bpop estimates
                  , y_t = y_tot
                  #, begin.esa.year = min(which(y_esa!=0))
                  # envi cov
                  , prcp = prcp[10:length(prcp)], dx32 = dx32[10:length(dx32)]
                  , ponds = c(ponds.std[1:length(ponds.std)],NA)
                  # envi cov from previous year
                  , prev.prcp = prcp[9:(length(prcp)-1)], prev.dx32 = dx32[9:(length(dx32)-1)]
)

# bugs.data <- list(nyrs=dim(AGWT.marray)[1],
#                   #pre-hunting band-recoveries (calendar year 2000-2020)
#                   recovmat.am = AGWT.marray[,,2,3], recovmat.af = AGWT.marray[,,2,4], recovmat.jm = AGWT.marray[,,2,1], recovmat.jf = AGWT.marray[,,2,2]
#                   #post-hunting band-recoveries (calendar year 2001-2020)
#                   , recovmatP.am = AGWT.marray[2:nyrs,2:(nyrs+1),1,3], recovmatP.af = AGWT.marray[2:nyrs,2:(nyrs+1),1,4], recovmatP.jm = AGWT.marray[2:nyrs,2:(nyrs+1),1,1], recovmatP.jf = AGWT.marray[2:nyrs,2:(nyrs+1),1,2]
#                   #pre-hunting total bandings
#                   , relmat.am = AGWT.release[,2,3], relmat.af = AGWT.release[,2,4], relmat.jm = AGWT.release[,2,1], relmat.jf = AGWT.release[,2,2]
#                   #post-hunting total bandings
#                   , relmatP.am = AGWT.release[2:nyrs,1,3], relmatP.af = AGWT.release[2:nyrs,1,4], relmatP.jm = AGWT.release[2:nyrs,1,1], relmatP.jf = AGWT.release[2:nyrs,1,2]
#                   # wing data
#                   , W.jv = jf.wing$Number, W.tot = female.wing$sum
#                   #Bpop estimates
#                   , y_t = y_tot
#                   # envi cov
#                   , prcp = prcp[10:length(prcp)], dx32 = dx32[10:length(dx32)]
#                   , ponds = c(ponds.std, NA)
#                   # envi cov from previous year
#                   , prev.prcp = prcp[9:length(prcp)-1], prev.dx32 = dx32[9:length(dx32)-1]
#                   )

# Initial values
inits <- function(){list(f.am=runif(nyrs,0,0.5),f.af=runif(nyrs,0,0.5),f.jm=runif(nyrs,0,0.5),f.jf=runif(nyrs,0,0.5)
                         , beta = rnorm(1, 0, 1)
                         # initial pop for mallard
                         , n_HY_mal = runif(1, 0, 50), n_AHY_mal = runif(1, 0, 50), n_HY_fem = runif(1, 0, 25), n_AHY_fem = runif(1, 0, 25)
                        # initial pop for agwt
                           # , n_HY_mal = runif(1, 0, 200), n_AHY_mal = runif(1, 0, 200), n_HY_fem = runif(1, 0, 200), n_AHY_fem = runif(1, 0, 200)
                          , alpha_prcp = rnorm(4, 0, 1), alpha_dx32 = rnorm(4, 0, 1)
                         , beta_pond = rnorm(1, 0, 1), beta_prcp = rnorm(1, 0, 1), beta_dx32 = rnorm(1, 0, 1)
                          ,gamma_prcp = rnorm(4, 0, 1), gamma_dx32 = rnorm(4, 0, 1)

   )}

# Parameters
parameters <- c("SH.am", "SN.am", "SH.af", "SN.af", "SH.jm", "SN.jm", "SH.jf", "SN.jf"
                , "S.am", "S.af", "S.jm", "S.jf", "f.am", "f.af", "f.jm", "f.jf"
                ,"R", "sd.R", "v", "q"
                , "N_tot", "N_HY_fem", "N_AHY_fem", "N_HY_mal", "N_AHY_mal", "sig.t", "beta", "delta", "lambda"
                , "muH.am", "muH.af", "muH.jm", "muH.jf", "muN.am", "muN.af", "muN.jm", "muN.jf"
                , "alpha_prcp", "gamma_prcp"
                , "alpha_dx32", "gamma_dx32"
                ,"beta_pond", "beta_prcp", "beta_dx32"
)

# MCMC settings
ni <- 300
nt <- 1
nb <- 100
nc <- 3
na <- 1000

# Run model
## MALL
mall.ipm <- jagsUI(bugs.data, inits=inits, parameters, "mall_agwt_ipm_clean.jags", n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE, store.data=TRUE)

# # Save model output
# saveRDS(mall.ipm, "mall_ipm_2000-2020_output.rda")

## AGWT
agwt.ipm <- jagsUI(bugs.data, inits=inits, parameters, "mall_agwt_ipm_clean.jags", n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE, store.data=TRUE)

# Save output
# saveRDS(agwt.ipm, 'agwt_ipm_2000_2020_output.rda')


#------Plotting model results--------
## Load model output (if needed)
# mall.ipm <- readRDS("mall_ipm_2000_2020_output.rda")
# agwt.ipm <- readRDS("agwt_ipm_2000_2020_output.rda")

# Set color palette
cbbPalette <- c("#E69F00","#009E73", "#CC79A7","#F0E442", "#0072B2", "#D55E00")

# custom plot theme
theme_ipm <- function(){ theme(
  axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
  axis.text.y = element_text(size = 5),
  #axis.title= element_blank(),
  plot.title = element_text(size=6, hjust = 0.5),
  #legend.position = "none",
  #legend.position = "bottom",
  #legend.title = element_text( size = 18),
  #legend.text = element_text(size = 16),
  text = element_text(size = 7),
  panel.background = element_rect(fill = "white", colour = "white",
                                  linewidth = 0.5, linetype = "solid"),
  panel.grid.major = element_line(linewidth = 0.2, linetype = 'solid',
                                  colour = "lightgrey"),
  panel.border = element_blank(),
  axis.line = element_line(colour = "black", linewidth = 0.2),
  legend.position = "bottom"
  )
}

# Grab phylopic
# mall.id <- get_uuid(name = "Anas platyrhynchos", n = 1)
# mall.img <- pick_phylopic(name = "Anas platyrhynchos", n = 1)
# agwt.id <- get_uuid(name = "Anas crecca", n = 4)
# agwt.img <- pick_phylopic(name = "Anas crecca", uuid = agwt.id, n = 4)

# ----- Plotting seasonal survival
# Data frame with both hunting and post-hunting survival
mall.h.ph <- data.frame(cbind(Year = rep(c(2000:2020), 8), 
                              Class = rep(rep(c("Juvenile Male", "Juvenile Female", "Adult Male", "Adult Female"), each = length(mall.ipm$mean$R)), 2), 
                              Season = rep(c("Fall/Winter", "Spring/Summer"), each = 84 #116
                                           ), 
                              Mean = c(mall.ipm$mean$SH.jm, NA, mall.ipm$mean$SH.jf, NA, mall.ipm$mean$SH.am, NA, mall.ipm$mean$SH.af, NA, NA, mall.ipm$mean$SN.jm, NA, mall.ipm$mean$SN.jf, NA, mall.ipm$mean$SN.am, NA, mall.ipm$mean$SN.af), 
                              Lower = c(mall.ipm$q2.5$SH.jm, NA, mall.ipm$q2.5$SH.jf, NA, mall.ipm$q2.5$SH.am, NA, mall.ipm$q2.5$SH.af, NA, NA, mall.ipm$q2.5$SN.jm, NA, mall.ipm$q2.5$SN.jf, NA, mall.ipm$q2.5$SN.am, NA, mall.ipm$q2.5$SN.af), 
                              Upper = c(mall.ipm$q97.5$SH.jm, NA, mall.ipm$q97.5$SH.jf, NA, mall.ipm$q97.5$SH.am, NA, mall.ipm$q97.5$SH.af, NA, NA, mall.ipm$q97.5$SN.jm, NA, mall.ipm$q97.5$SN.jf, NA, mall.ipm$q97.5$SN.am, NA, mall.ipm$q97.5$SN.af)))

mall.h.ph <- mall.h.ph %>% 
                mutate_at(c("Mean", "Lower", "Upper"), as.numeric)

#plot mallard seasonal survivals
ggplot(mall.h.ph, aes(x = Year, y = Mean, group = Season)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper, color = Season), size = 0.2, alpha = 0.75, position=position_dodge(width=0.35)) +
  theme_ipm() +
  scale_color_manual(values=c("royalblue3", "firebrick3"))+ 
  ylab("Survival probability") +
  scale_x_discrete(breaks = seq(1992, 2020, by = 2)) +
  facet_wrap(~Class)

# save plot
ggsave("figures/mall_seasonal_survivals_2000_2020.jpg", units="cm", width=10, height=10, dpi=600)

# green-winged teal seasonal survival data frame
agwt.h.ph <- data.frame(cbind(Year = rep(c(2000:2020), 8), 
                              Class = rep(rep(c("Juvenile Male", "Juvenile Female", "Adult Male", "Adult Female"), each = length(agwt.ipm$mean$R)), 2), 
                              Season = rep(c("Fall/Winter", "Spring/Summer"), each = 84), 
                              Mean = c(agwt.ipm$mean$SH.jm, NA, agwt.ipm$mean$SH.jf, NA, agwt.ipm$mean$SH.am, NA, agwt.ipm$mean$SH.af, NA, NA, agwt.ipm$mean$SN.jm, NA, agwt.ipm$mean$SN.jf, NA, agwt.ipm$mean$SN.am, NA, agwt.ipm$mean$SN.af), 
                              Lower = c(agwt.ipm$q2.5$SH.jm, NA, agwt.ipm$q2.5$SH.jf, NA, agwt.ipm$q2.5$SH.am, NA, agwt.ipm$q2.5$SH.af, NA, NA, agwt.ipm$q2.5$SN.jm, NA, agwt.ipm$q2.5$SN.jf, NA, agwt.ipm$q2.5$SN.am, NA, agwt.ipm$q2.5$SN.af), 
                              Upper = c(agwt.ipm$q97.5$SH.jm, NA, agwt.ipm$q97.5$SH.jf, NA, agwt.ipm$q97.5$SH.am, NA, agwt.ipm$q97.5$SH.af, NA, NA, agwt.ipm$q97.5$SN.jm, NA, agwt.ipm$q97.5$SN.jf, NA, agwt.ipm$q97.5$SN.am, NA, agwt.ipm$q97.5$SN.af)))

agwt.h.ph <- agwt.h.ph %>% 
  mutate_at(c("Mean", "Lower", "Upper"), as.numeric)

# plot green-winged teal seasonal survivals
ggplot(agwt.h.ph, aes(x = Year, y = Mean, group = Season)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper, color = Season), size = 0.2, alpha = 0.75, position=position_dodge(width=0.35)) +
  theme_ipm() +
  scale_color_manual(values=c("royalblue3", "firebrick3"))+ 
  ylab("Survival probability") +
  scale_x_discrete(breaks = seq(1992, 2020, by = 2)) +
  facet_wrap(~Class)

# save plot
ggsave("figures/agwt_seasonal_survivals_2000_2020.jpg", units="cm", width=10, height=10, dpi=600)

# ------ Plotting annual survival
# adult and juvenile annual survival
# MALL
S.a.mall.out <- data.frame(cbind(Year = rep(c(2000:2019),4), 
                 Class = rep(c("Juvenile Male", "Juvenile Female", "Adult Male", "Adult Female"), each = length(2000:2019)), 
                 Mean = c(mall.ipm$mean$S.jm, mall.ipm$mean$S.jf, mall.ipm$mean$S.am, mall.ipm$mean$S.af), 
                 Lower = c(mall.ipm$q2.5$S.jm, mall.ipm$q2.5$S.jf,mall.ipm$q2.5$S.am, mall.ipm$q2.5$S.af), 
                 Upper = c(mall.ipm$q97.5$S.jm, mall.ipm$q97.5$S.jf, mall.ipm$q97.5$S.am, mall.ipm$q97.5$S.af)))

S.a.mall.out <- S.a.mall.out %>% 
  mutate_at(c('Mean', 'Lower', 'Upper'), as.numeric)

mall.a.plot <- ggplot(S.a.mall.out, aes(x = Year, y = Mean, group = Class)) +
  #geom_point(aes(color = Class)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper, color = Class),size = 0.01, alpha = 0.75) +
  theme_ipm() +
  theme(legend.position = "none") +
  ylab("Annual survival") +
  ylim(0,1) +
  scale_color_manual(values=cbbPalette) +
  scale_x_discrete(breaks = seq(1992, 2019, by = 2)) +
  facet_wrap(~Class)

ggsave("figures/mall_annual_survival_2000_2020.jpg", units="cm", width=10, height=10, dpi=600)

# AGWT
S.a.agwt.out <- data.frame(cbind(Class = rep(c("Juvenile Male", "Juvenile Female", "Adult Male", "Adult Female"), each= length(2000:2019)), Year = rep(2000:2019, 4), S_a.mean = c(agwt.ipm$mean$S.jm, agwt.ipm$mean$S.jf, agwt.ipm$mean$S.am, agwt.ipm$mean$S.af), S_a.lower = c(agwt.ipm$q2.5$S.jm, agwt.ipm$q2.5$S.jf,agwt.ipm$q2.5$S.am, agwt.ipm$q2.5$S.af), S_a.upper = c(agwt.ipm$q97.5$S.jm, agwt.ipm$q97.5$S.jf, agwt.ipm$q97.5$S.am, agwt.ipm$q97.5$S.af)))

S.a.agwt.out <- S.a.agwt.out %>% 
  mutate_at(c('S_a.mean', 'S_a.lower', 'S_a.upper'), as.numeric)

agwt.a.plot <- ggplot(S.a.agwt.out, aes(x = Year, y = S_a.mean, group = Class)) +
  #geom_point(aes(color = Class)) +
  geom_pointrange(aes(ymin = S_a.lower, ymax = S_a.upper, color = Class),size = 0.01,alpha = 0.75,position=position_dodge(width=0.35)) +
  theme_ipm() +
  theme(legend.position = "none") +
  ylab("") +
  ylim(0,1) +
  scale_x_discrete(breaks = seq(1992, 2019, by = 2)) +
  scale_color_manual(values=cbbPalette) + 
  facet_wrap(~Class)

ggsave("figures/agwt_annual_survival_2000_2019.jpg", units="cm", width=10, height=10, dpi=600)

ggarrange(mall.a.plot, agwt.a.plot, ncol=2, legend="none") + theme(plot.margin = margin(0.1,0.2,0.2,0.2, "cm")) 

# ------ Plotting recovery probabilities
# MALL
mall.f.out <- data.frame(cbind(Class = rep(c("Juvenile Male", "Juvenile Female", "Adult Male", "Adult Female"), each = length(2000:2020)), Year = rep(2000:2020, 4), f.mean = c(mall.ipm$mean$f.jm, mall.ipm$mean$f.jf, mall.ipm$mean$f.am, mall.ipm$mean$f.af), f.lower = c(mall.ipm$q2.5$f.jm, mall.ipm$q2.5$f.jf, mall.ipm$q2.5$f.am, mall.ipm$q2.5$f.af), f.upper = c(mall.ipm$q97.5$f.jm, mall.ipm$q97.5$f.jf, mall.ipm$q97.5$f.am, mall.ipm$q97.5$f.af)))

mall.f.out <- mall.f.out %>% 
                mutate_at(c('f.mean', 'f.lower', 'f.upper'), as.numeric)


mall.recov.plot <- ggplot(mall.f.out, aes(x = Year, y = f.mean, group = Class)) +
                    geom_pointrange(aes(ymin = f.lower, ymax = f.upper, color = Class), size = 0.01, alpha = 0.7) +
  ggtitle("Mallard") +
                    ylab("Recovery probability") +
  # add_phylopic(img = mall.img, color = "black", x = 17, y = 0.22, ysize = 0.05) +
  scale_color_manual(values=cbbPalette)+ 
  ylim(0, 0.25) +
  xlab("") +
  guides(color=guide_legend(nrow=1)) +
  theme_ipm()

# AGWT
agwt.f.out <- data.frame(cbind(Class = rep(c("Juvenile Male", "Juvenile Female", "Adult Male", "Adult Female"), each = length(2000:2020)), Year = rep(2000:2020, 4), f.mean = c(agwt.ipm$mean$f.jm, agwt.ipm$mean$f.jf, agwt.ipm$mean$f.am, agwt.ipm$mean$f.af), f.lower = c(agwt.ipm$q2.5$f.jm, agwt.ipm$q2.5$f.jf, agwt.ipm$q2.5$f.am, agwt.ipm$q2.5$f.af), f.upper = c(agwt.ipm$q97.5$f.jm, agwt.ipm$q97.5$f.jf, agwt.ipm$q97.5$f.am, agwt.ipm$q97.5$f.af)))

agwt.f.out <- agwt.f.out %>% 
  mutate_at(c('f.mean', 'f.lower', 'f.upper'), as.numeric)


agwt.recov.plot <- ggplot(agwt.f.out, aes(x = Year, y = f.mean, group = Class)) +
  geom_pointrange(aes(ymin = f.lower, ymax = f.upper, color = Class), size = 0.01, alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # add_phylopic(img = agwt.img, color = "black", x = 17, y = 0.01, ysize = 0.03) +
  ggtitle("Green-winged Teal") + 
  ylab("Recovery probability") +
  xlab("") +
  scale_color_manual(values=cbbPalette)+ 
  ylim(0, 0.15) +
  guides(color=guide_legend(nrow=1)) +
  theme_ipm()

# Putting them together
mall.recov.plot / agwt.recov.plot + plot_layout(guides = "collect", widths = 1.5, heights = unit(c(3, 3), 'cm')) & theme(legend.position = 'bottom')

ggsave("figures/recovery_prob_plot_both.jpg", units="cm", width=11, height=11, dpi=600)
# ggsave("figures/agwt_recovery_plot_new_1992_2020.jpg", units="cm", width=10, height=10, dpi=600)

# ------ Plotting age ratio (productivity)

# custom theme
theme_recruit <- function(){ theme(
  axis.text.x = element_text(size = 5, angle = 60, vjust = 0.5, hjust = 0.5),
  axis.text.y = element_text(size = 5),
  #plot.title = element_text(size=16, hjust = 0.5),
  text = element_text(size = 6),
  panel.background = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(colour = "black"), 
  legend.position = "bottom",
  legend.margin = margin(0,0,0,0),
  legend.box.margin = margin(0,0,0,0)
)
}
sp.col <- c("Dark green", "Chocolate4")

# Model output for both species
R.out <- data.frame(cbind(Year = 2000:2020, Species = c(rep(c("Mallard", "Green-winged Teal"), each = length(agwt.ipm$mean$R))), R.mean = c(mall.ipm$mean$R, agwt.ipm$mean$R), R.lower = c(mall.ipm$q2.5$R, agwt.ipm$q2.5$R), R.upper = c(mall.ipm$q97.5$R, agwt.ipm$q97.5$R)))

R.out <- R.out %>% 
  mutate_at(c('Year', 'R.mean', 'R.lower', 'R.upper'), as.numeric)

# Plot
R.plot <- ggplot(R.out, aes(x = Year, y = R.mean, group = Species, fill = Species)) +
  # add_phylopic(img = agwt.img, color = "black", x = 2019, y = 3, ysize = 0.8) +
  # add_phylopic(img = mall.img, color = "black", x = 2019, y = 0.35, ysize = 0.7) +
  geom_line(aes(color = Species), size = 0.5) +
  geom_ribbon(aes(ymin = R.lower, ymax = R.upper), alpha = 0.3) +
  theme_recruit() +
  scale_color_manual(values=sp.col) +
  scale_fill_manual(values = sp.col) +
  scale_x_continuous(labels = 2000:2020, breaks = seq(2000, 2020, 1)) +
  ylim(0, 5) +
  ylab("Age ratio")


ggsave("figures/age_ratio_new_both.jpg", units="cm", width = 9, height=7, dpi=600)

# vulnerability
v.out <- data.frame(cbind(Year = 2000:2020, Species = c(rep(c("Mallard", "Green-winged Teal"), each = length(agwt.ipm$mean$v))), v.mean = c(mall.ipm$mean$v, agwt.ipm$mean$v), v.lower = c(mall.ipm$q2.5$v, agwt.ipm$q2.5$v), v.upper = c(mall.ipm$q97.5$v, agwt.ipm$q97.5$v)))

v.out <- v.out %>% 
  mutate_at(c('Year','v.mean', 'v.lower', 'v.upper'), as.numeric)

ggplot(v.out, aes(x = Year, y = v.mean, group = Species, fill = Species)) +
  # add_phylopic(img = agwt.img, color = "black", x = 2019, y = 0.4, ysize = 1) +
  # add_phylopic(img = mall.img, color = "black", x = 2019, y = 3.5, ysize = 1) +
  geom_line(aes(color = Species), linewidth = 0.5) +
  geom_ribbon(aes(ymin = v.lower, ymax = v.upper), alpha = 0.3) +
  theme_recruit() +
  scale_color_manual(values=sp.col) +
  scale_fill_manual(values = sp.col) +
  ylim(0,6) +
  ylab("Vulnerability")

ggsave("figures/vulnerability_new_both.jpg", units="cm", width=8, height=7, dpi=600)

# overall population (bpop)
Bpop.both <- data.frame(cbind(Year = 2000:2020, Species = c(rep(c("Mallard", "Green-winged Teal"), each = length(agwt.ipm$mean$N_tot))), mean = c(mall.ipm$mean$N_tot, agwt.ipm$mean$N_tot), lower = c(mall.ipm$q2.5$N_tot, agwt.ipm$q2.5$N_tot), upper = c(mall.ipm$q97.5$N_tot, agwt.ipm$q97.5$N_tot)))

Bpop.both <- Bpop.both %>% 
  mutate_at(c('mean', 'lower', 'upper'), as.numeric)

ggplot(Bpop.both, aes(x = Year, y = mean, group = Species)) +
  geom_line(aes(color = Species), linewidth = 0.5) +
  # add_phylopic(img = mall.img, color = "black", x = 22, y = 1200, ysize = 150) +
  # add_phylopic(img = agwt.img, color = "black", x =26, y = 140, ysize = 220) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Species), alpha = 0.5) +
  theme_recruit()+
  scale_color_manual(values=sp.col) +
  scale_fill_manual(values = sp.col) +
  ylab("Population index (10,000s)")

ggsave("figures/total_pop_size_both.jpg", units="cm", width=8, height=7, dpi=600)

# Plotting population by age-sex class
mall.Bpop.out <- data.frame(cbind(Class = rep(c("Juvenile Male", "Juvenile Female", "Adult Male", "Adult Female"), each = length(mall.ipm$mean$N_AHY_fem)), Year = rep(2005:2020, 4), Bpop.mean = c(mall.ipm$mean$N_HY_mal, mall.ipm$mean$N_HY_fem, mall.ipm$mean$N_AHY_mal, mall.ipm$mean$N_AHY_fem), Bpop.lower = c(mall.ipm$q2.5$N_HY_mal, mall.ipm$q2.5$N_HY_fem, mall.ipm$q2.5$N_AHY_mal, mall.ipm$q2.5$N_AHY_fem), Bpop.upper = c(mall.ipm$q97.5$N_HY_mal, mall.ipm$q97.5$N_HY_fem, mall.ipm$q97.5$N_AHY_mal, mall.ipm$q97.5$N_AHY_fem)))

mall.Bpop.out <- mall.Bpop.out %>% 
                    mutate_at(c("Bpop.mean", "Bpop.lower", "Bpop.upper"), as.numeric)

mall.Bpop.plot <- ggplot(mall.Bpop.out, aes(x = Year, y = Bpop.mean, group = Class)) +
  #add_phylopic(img = agwt.img, color = "black", x = 5.5, y = 500, ysize = 100) +
  geom_line(aes(color = Class)) +
  geom_ribbon(aes(ymin = Bpop.lower, ymax = Bpop.upper, fill = Class ), alpha = 0.5)


agwt.Bpop.out <- data.frame(cbind(Class = rep(c("Juvenile Male", "Juvenile Female", "Adult Male", "Adult Female"), each = length(agwt.ipm$mean$N_AHY_fem)), Year = rep(1992:2020, 4), Bpop.mean = c(agwt.ipm$mean$N_HY_mal, agwt.ipm$mean$N_HY_fem, agwt.ipm$mean$N_AHY_mal, agwt.ipm$mean$N_AHY_fem), Bpop.lower = c(agwt.ipm$q2.5$N_HY_mal, agwt.ipm$q2.5$N_HY_fem, agwt.ipm$q2.5$N_AHY_mal, agwt.ipm$q2.5$N_AHY_fem), Bpop.upper = c(agwt.ipm$q97.5$N_HY_mal, agwt.ipm$q97.5$N_HY_fem, agwt.ipm$q97.5$N_AHY_mal, agwt.ipm$q97.5$N_AHY_fem)))

agwt.Bpop.out <- agwt.Bpop.out %>% 
  mutate_at(c("Bpop.mean", "Bpop.lower", "Bpop.upper"), as.numeric)

agwt.Bpop.plot <- ggplot(agwt.Bpop.out, aes(x = Year, y = Bpop.mean, group = Class)) +
  #add_phylopic(img = agwt.img, color = "black", x = 5.5, y = 500, ysize = 100) +
  geom_line(aes(color = Class)) +
  geom_ribbon(aes(ymin = Bpop.lower, ymax = Bpop.upper, fill = Class ), alpha = 0.5) +
  theme_recruit()


## Plotting environmental covariate effects

theme_cov <- function(){ theme(
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.title= element_text(size = 6),
  plot.title = element_text(size=8, hjust = 0.5),
  plot.margin = margin(0,0,0,0),
  text = element_text(size = 6),
  panel.background = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = "bottom",
  legend.key.size = unit(0.3, "cm")
)
}
# mall.ipm <- mall.ipm2
# covariate effects on fall/winter survival plot

alpha.prcp <- mall.ipm$sims.list$alpha_prcp %>%
  as.data.frame() %>% 
  rename("Adult Male" = V1, "Adult Female" = V2, "Juvenile Male" = V3, "Juvenile Female" = V4) %>% 
  pivot_longer(everything(), names_to = "Class", values_to = "Covariate_effect")

alpha.dx32 <- mall.ipm$sims.list$alpha_dx32 %>%
  as.data.frame() %>% 
  rename("Adult Male" = V1, "Adult Female" = V2, "Juvenile Male" = V3, "Juvenile Female" = V4) %>% 
  pivot_longer(everything(), names_to = "Class", values_to = "Covariate_effect")

alpha.all <- bind_rows(alpha.prcp, alpha.dx32) %>% 
  mutate(Covariate = rep(c("Precipitation", "dx32"), each = nrow(alpha.prcp)))

alpha.plot <- ggplot(alpha.all, aes(x = Covariate_effect, y = Covariate, fill = Class)) +
  stat_halfeye(interval_size_range = c(0.25, 0.45), position = "dodge")+
  theme_cov() +
  geom_vline(xintercept = 0, linetype="dashed", color = "black", linewidth=0.3) +
  scale_fill_manual(values=cbbPalette) +
  annotate("text", x=-0.9, y=2.5, label= "(a)") +
  ylab("Environmental Covariate") +
  xlab("Covariate effect estimate") +
  ggtitle("Fall-Winter Survival") +
  xlim(-1.5, 1.5) 

ggsave("mall_cov_effects_fall_winter_surv.jpg", units="cm", width=10, height=7, dpi=600)

# covariate effects on spring/summer survival plot

gamma.prcp <- mall.ipm$sims.list$gamma_prcp %>%
  as.data.frame() %>% 
  rename("Adult Male" = V1, "Adult Female" = V2, "Juvenile Male" = V3, "Juvenile Female" = V4) %>% 
  pivot_longer(everything(), names_to = "Class", values_to = "Covariate_effect")

gamma.dx32 <- mall.ipm$sims.list$gamma_dx32 %>%
  as.data.frame() %>% 
  rename("Adult Male" = V1, "Adult Female" = V2, "Juvenile Male" = V3, "Juvenile Female" = V4) %>% 
  pivot_longer(everything(), names_to = "Class", values_to = "Covariate_effect")

gamma.all <- bind_rows(gamma.prcp, gamma.dx32) %>% 
  mutate(Covariate = rep(c("Precipitation", "dx32"), each = nrow(gamma.prcp)))
  
gamma.plot <- ggplot(gamma.all, aes(x = Covariate_effect, y = Covariate, fill = Class)) +
  stat_halfeye(interval_size_range = c(0.25, 0.45), position = "dodge")+
  theme_cov() +
  geom_vline(xintercept = 0, linetype="dashed", color = "black", linewidth=0.3) +
  scale_fill_manual(values=cbbPalette) +
  annotate("text", x=-1.1, y=2.5, label= "(b)") +
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.5,
  #          alpha = .2,fill = "gray") +
  # annotate("text", x=c(1.3, 0.3), y=c(1.375, 2.375), label= "*") +
  ggtitle("Spring-Summer Survival") +
  xlim(-1.5, 1.5) +
  ylab("Environmental Covariate") +
  xlab("Covariate effect estimate")

alpha.plot + gamma.plot + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("mall_cov_effects_surv_2000_2020.jpg", units="cm", width=10, height=7, dpi=600)

# covariate effects on productivity plot
mall.beta.prcp <- mall.ipm$sims.list$beta_prcp %>%
  as.data.frame()

mall.beta.dx32 <-  mall.ipm$sims.list$beta_dx32 %>%
  as.data.frame()

# mall.beta.nao <-  mall.ipm$sims.list$beta_nao %>%
#   as.data.frame()

mall.beta.pond <-  mall.ipm$sims.list$beta_pond %>%
  as.data.frame()

mall.beta.all <- bind_rows(mall.beta.prcp, mall.beta.dx32, mall.beta.pond) %>% 
            rename("Covariate_effect" = ".") %>% 
            mutate(Covariate = rep(c("Precipitation", "dx32", "Pond"), each = nrow(mall.beta.prcp)))
  
mall.beta.plot <- ggplot(mall.beta.all, aes(x = Covariate_effect, y = Covariate)) +
  stat_halfeye(fill = "dodgerblue2", interval_size_range = c(0.25, 0.45), position = "dodge")+
  theme_cov() +
  geom_vline(xintercept = 0, linetype="dashed", color = "black", linewidth=0.3) +
  annotate("text", x = 0.4, y = 4, label = "*") +
  xlim(-0.5, 0.5) +
  annotate("text", x=-0.45, y=3.75, label= "(a)") +
  annotate("text", x=0.4, y=3, label= "*") +
  labs(x = "Covariate effect estimate", y = "Environmetnal covariate") +
  ggtitle("Mallard")

# combine all plots
# (alpha.plot + gamma.plot)/mall.beta.plot +  plot_layout(guides = "collect", widths = 1.5, heights = unit(c(6, 3), 'cm')) & theme(legend.position = 'bottom')

ggsave("figures/cov_effect_prod_mall_new_2005_2020.jpg",  units="cm", width=12, height=12, dpi=600)

# ------- Calculating proportion of posterior distribution greater than 0
mall.alpha.prcp <- mall.ipm$sims.list$alpha_prcp
mall.alpha.dx32 <- mall.ipm$sims.list$alpha_dx32
mall.alpha.nao <- mall.ipm$sims.list$alpha_nao
mall.gamma.prcp <- mall.ipm$sims.list$gamma_prcp
mall.gamma.dx32 <- mall.ipm$sims.list$gamma_dx32
mall.beta.prcp <- mall.ipm$sims.list$beta_prcp
mall.beta.dx32 <- mall.ipm$sims.list$beta_dx32
mall.beta.pond <- mall.ipm$sims.list$beta_pond

1-apply(mall.alpha.prcp, 2, function(c) ecdf(c)(0))
1-apply(mall.alpha.dx32, 2, function(c) ecdf(c)(0))
1-apply(mall.gamma.prcp, 2, function(c) ecdf(c)(0))
1-apply(mall.gamma.dx32, 2, function(c) ecdf(c)(0))

prob.less.0 <- ecdf(mall.beta.pond)
(prob.greater.0 <- 1-prob.less.0(0))
prob.less.0 <- ecdf(mall.beta.prcp)
(prob.greater.0 <- 1-prob.less.0(0))
prob.less.0 <- ecdf(mall.beta.dx32)
(prob.greater.0 <- 1-prob.less.0(0))

## AGWT cov plot

agwt.alpha.prcp <- agwt.ipm$sims.list$alpha_prcp %>%
  as.data.frame() %>% 
  rename("Adult Male" = V1, "Adult Female" = V2, "Juvenile Male" = V3, "Juvenile Female" = V4) %>% 
  pivot_longer(everything(), names_to = "Class", values_to = "Covariate_effect")

agwt.alpha.dx32 <- agwt.ipm$sims.list$alpha_dx32 %>%
  as.data.frame() %>% 
  rename("Adult Male" = V1, "Adult Female" = V2, "Juvenile Male" = V3, "Juvenile Female" = V4) %>% 
  pivot_longer(everything(), names_to = "Class", values_to = "Covariate_effect")

agwt.alpha.all <- bind_rows(agwt.alpha.prcp, agwt.alpha.dx32) %>% 
  mutate(Covariate = rep(c("Precipitation", "dx32"), each = nrow(agwt.alpha.prcp)))

agwt.alpha.plot <- ggplot(agwt.alpha.all, aes(x = Covariate_effect, y = Covariate, fill = Class)) +
  stat_halfeye(interval_size_range = c(0.25, 0.45), position = "dodge")+
  theme_cov() +
  geom_vline(xintercept = 0, linetype="dashed", color = "black", linewidth=0.3) +
  scale_fill_manual(values=cbbPalette) +
  annotate("text", x=-1.1, y=2.5, label= "(a)") +
  ggtitle("Fall-Winter Survival") +
  labs(x = "Covariate effect estimate", y = "Environmetnal covariate") +
  xlim(-1.5, 1.5) 

# covariate effects on spring/summer survival plot

agwt.gamma.prcp <- agwt.ipm$sims.list$gamma_prcp %>%
  as.data.frame() %>% 
  rename("Adult Male" = V1, "Adult Female" = V2, "Juvenile Male" = V3, "Juvenile Female" = V4) %>% 
  pivot_longer(everything(), names_to = "Class", values_to = "Covariate_effect")

agwt.gamma.dx32 <- agwt.ipm$sims.list$gamma_dx32 %>%
  as.data.frame() %>% 
  rename("Adult Male" = V1, "Adult Female" = V2, "Juvenile Male" = V3, "Juvenile Female" = V4) %>% 
  pivot_longer(everything(), names_to = "Class", values_to = "Covariate_effect")

agwt.gamma.all <- bind_rows(agwt.gamma.prcp, agwt.gamma.dx32) %>% 
  mutate(Covariate = rep(c("Precipitation", "dx32"), each = nrow(agwt.gamma.prcp)))

agwt.gamma.plot <- ggplot(agwt.gamma.all, aes(x = Covariate_effect, y = Covariate, fill = Class)) +
  stat_halfeye(interval_size_range = c(0.25, 0.45), position = "dodge")+
  theme_cov() +
  geom_vline(xintercept = 0, linetype="dashed", color = "black", linewidth=0.3) +
  scale_fill_manual(values=cbbPalette) +
 annotate("text", x = -1.1, y = 2.5, label = "(b)") +
  ggtitle("Spring-Summer Survival") +
  labs(x = "Covariate effect estimate", y = "Environmetnal covariate") +
  xlim(-1.5, 1.5) 


agwt.alpha.plot + agwt.gamma.plot + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("agwt_cov_effects_surv_2000_2020.jpg", units = "cm", width = 10, height = 7, dpi = 600)

# covariate effects on productivity plot
agwt.beta.prcp <- agwt.ipm$sims.list$beta_prcp %>%
  as.data.frame()

agwt.beta.dx32 <-  agwt.ipm$sims.list$beta_dx32 %>%
  as.data.frame()

agwt.beta.pond <-  agwt.ipm$sims.list$beta_pond %>%
  as.data.frame()

agwt.beta.all <- bind_rows(agwt.beta.prcp, agwt.beta.dx32, agwt.beta.pond) %>% 
  rename("Covariate_effect" = ".") %>% 
  mutate(Covariate = rep(c("Precipitation", "dx32", "Pond"), each = nrow(agwt.beta.prcp)))

agwt.beta.plot <- ggplot(agwt.beta.all, aes(x = Covariate_effect, y = Covariate)) +
  stat_halfeye(fill = "dodgerblue2", interval_size_range = c(0.25, 0.45)) +
  theme_cov() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  xlim(-0.5, 0.5) +
  annotate("text", x = -0.45, y = 3.75, label= "(b)") +
  labs(x = "Covariate effect estimate", y = "Environmetnal covariate") +
  ggtitle("Green-winged teal") 

mall.beta.plot + agwt.beta.plot + plot_layout(guides = "collect")

ggsave("figures/cov_effects_prod_both.jpg",  units = "cm", width = 10, height = 5, dpi = 600)

# Calculating proportion of posterior distribution greater than 0
agwt.alpha.prcp <- agwt.ipm$sims.list$alpha_prcp
agwt.alpha.dx32 <- agwt.ipm$sims.list$alpha_dx32
agwt.gamma.prcp <- agwt.ipm$sims.list$gamma_prcp
agwt.gamma.dx32 <- agwt.ipm$sims.list$gamma_dx32
agwt.beta.prcp <- agwt.ipm$sims.list$beta_prcp
agwt.beta.dx32 <- agwt.ipm$sims.list$beta_dx32
agwt.beta.pond <- agwt.ipm$sims.list$beta_pond

1-apply(agwt.alpha.prcp, 2, function(c) ecdf(c)(0))
1-apply(agwt.alpha.dx32, 2, function(c) ecdf(c)(0))
1-apply(agwt.gamma.prcp, 2, function(c) ecdf(c)(0))
1-apply(agwt.gamma.dx32, 2, function(c) ecdf(c)(0))


prob.less.0 <- ecdf(agwt.beta.pond)
(prob.greater.0 <- 1-prob.less.0(0))
prob.less.0 <- ecdf(agwt.beta.prcp)
(prob.greater.0 <- 1-prob.less.0(0))
prob.less.0 <- ecdf(agwt.beta.dx32)
(prob.greater.0 <- 1-prob.less.0(0))


# Make band summary tables for appendix
dims <- dim(MALL.marray)
totRecov <- array(NA, c(dims[1], dims[3], dims[4]))
for(i in 1:dims[3]){
  for(j in 1:dims[4]){
    totRecov[,i,j] <- rowSums(AGWT.marray[, 1:dims[1], i, j])
  }
}

totRecov <- apply(totRecov, 2, identity)
totRecov <- totRecov %>% 
  data.frame() %>% 
  mutate(class = rep(c('JM', 'JF', 'AM', 'AF'), each = dims[1])) %>% 
  rename('Pre-Season' = 'X2', 'Post-Season' = 'X1') %>% 
  write.csv('AGWT_recoveries.csv')

totRelease <- apply(AGWT.release, 2, identity)
totRelease <- totRelease %>% 
  data.frame() %>% 
  mutate(class = rep(c('JM', 'JF', 'AM', 'AF'), each = dims[1])) %>% 
  write.csv('AGWT_releases.csv')
