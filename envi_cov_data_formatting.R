library(tidyverse)
library(GGally)
library(rnoaa)

# Load climate data by state
ga.dat <- read.csv("raw_dat/environmental_covariates/GA_ClimDat.csv")
fl.dat <- read.csv("raw_dat/environmental_covariates/FL_ClimDat.csv")
sc.dat <- read.csv("raw_dat/environmental_covariates/SC_ClimDat.csv")
nc.dat <- read.csv("raw_dat/environmental_covariates/NC_ClimDat.csv")
tn.dat <- read.csv("raw_dat/environmental_covariates/TN_ClimDat.csv")
al.dat <- read.csv("raw_dat/environmental_covariates/AL_ClimDat.csv")
ms.dat <- read.csv("raw_dat/environmental_covariates/MS_ClimDat.csv")
ar.dat <- read.csv("raw_dat/environmental_covariates/AR_ClimDat.csv")
la.dat <- read.csv("raw_dat/environmental_covariates/LA_ClimDat.csv")
tx.dat <- read.csv("raw_dat/environmental_covariates/TX_ClimDat.csv")
ky.dat <- read.csv("raw_dat/environmental_covariates/KY_ClimDat.csv")
mo.dat <- read.csv("raw_dat/environmental_covariates/MO_ClimDat.csv")
nm.dat <- read.csv("raw_dat/environmental_covariates/NM_ClimDat.csv")
az.dat <- read.csv("raw_dat/environmental_covariates/AZ_ClimDat.csv")
ok.dat <- read.csv("raw_dat/environmental_covariates/OK_ClimDat.csv")
ks.dat <- read.csv("raw_dat/environmental_covariates/KS_ClimDat.csv")


# Climate data from mid-latitute cities
chicago <- read.csv("raw_dat/environmental_covariates/Chicago.csv")
nplatte <- read.csv("raw_dat/environmental_covariates/NorthPlatte.csv")
desmoines <- read.csv("raw_dat/environmental_covariates/DesMoines.csv")
cleveland <- read.csv("raw_dat/environmental_covariates/Cleveland.csv")

# Combine all
dat <- rbind(ga.dat, fl.dat, sc.dat, nc.dat, tn.dat, al.dat, ms.dat, ar.dat, la.dat, tx.dat, ky.dat, mo.dat, nm.dat, az.dat, ok.dat, ks.dat)

# Weather stations in Mexico
mx.dat <- dat[str_sub(dat$NAME, -2, -1)=="MX",]

# Remove data from Mexico
dat <- dat[!(str_sub(dat$NAME, -2, -1)=="MX"),]
dat <- dat %>% 
  distinct(.keep_all = TRUE)

dat$YEAR <- str_sub(dat$DATE, 1, 4)
dat$MONTH <- str_sub(dat$DATE, 6, 7)
dat$STATE <- str_sub(dat$NAME, -5, -4)

# Identify southern states
states <- c("GA", "FL", "SC", "NC", "TN", "AL", "MS", "AR", "LA", "TX", "KY", "MO", "NM", "AZ", "OK", "KS")

# Identify winter months
w.months <- c("01", "09", "10", "11", "12")

#subset only southern states
dat <- dat[dat$STATE%in%states, ]

# Subset winter covariates
# DT32: number of days per month where minimum temp was below freezing
# DX32: number of days per month where maximum temp was below freezing
# PRCP: monthly amount of precipitation (in)
winter <- dat %>%
              select(STATION, NAME, DT32, DX32, PRCP, YEAR, MONTH, STATE) %>%
              filter(MONTH%in%c("01", "09", "10", "11", "12"))

winter$YEAR <- as.numeric(winter$YEAR)

# January of following year are considered part of the previous winter
year.new <- rep(NA, nrow(winter))
for(i in 1:nrow(winter)){
  year.new[i] <- ifelse(winter$MONTH[i]%in%c("01"), winter[i,6]-1, winter[i,6])
}

# Incorporate newly assigned year to dataframe
winter <- winter %>% 
        mutate(YEAR2 = year.new) %>% 
        drop_na() %>% 
        group_by(STATE, NAME, YEAR2) %>% 
        summarise(n_months = n()) %>% 
  pivot_wider(names_from = YEAR2, values_from = n_months) %>% 
  select(-c('1990', '2021')) %>% 
  mutate_if(is.integer, as.numeric)

winter$Sum_mo <- rowSums(winter[,3:ncol(winter)])
winter <- winter %>% 
          filter(Sum_mo==150)

# Winter precipitation table
# winter.prcp <- winter %>% 
#                select(NAME, PRCP, YEAR2, MONTH, STATE) %>% 
#                group_by(STATE, NAME, YEAR2) %>% 
#                summarise(n_months = n()) %>% 
#               pivot_wider(names_from = YEAR2, values_from = n_months) %>% 
#               drop_na() %>% 
#               select(-c('1990', '2021')) %>% 
#               mutate_if(is.integer, as.numeric) %>% 
#               mutate(Sum_mo = rowSums(winter[,3:ncol(winter)]))
# 
# winter.prcp$Sum_mo <- rowSums(winter.prcp[,3:ncol(winter.prcp)])
# winter.prcp <- winter.prcp %>% 
#                 filter(Sum_mo==180)
# 
# write.csv(winter.prcp, "Winter_precipitation.csv")
# 
# # Number of days minimum temp below freezing
# winter.dt32 <- winter %>% 
#   select(NAME, DT32, YEAR2, MONTH, STATE) %>% 
#   group_by(STATE, NAME, YEAR2) %>% 
#   summarise(n_months = n()) %>% 
#   pivot_wider(names_from = YEAR2, values_from = n_months) %>% 
#   drop_na() %>% 
#   select(-c('1990', '2021')) %>% 
#   mutate_if(is.integer, as.numeric)
# 
# winter.dt32$Sum_mo <- rowSums(winter.dt32[,3:ncol(winter.dt32)])
# winter.dt32 <- winter.dt32 %>% 
#   filter(Sum_mo==180)
# 
# # Number of days maximum temp below freezing
# winter.dx32 <- winter %>% 
#   select(NAME, DX32, YEAR2, MONTH, STATE) %>% 
#   group_by(STATE, NAME, YEAR2) %>% 
#   summarise(n_months = n()) %>% 
#   pivot_wider(names_from = YEAR2, values_from = n_months) %>% 
#   drop_na() %>% 
#   select(-c('1990', '2021')) %>% 
#   mutate_if(is.integer, as.numeric)
# 
# winter.dx32$Sum_mo <- rowSums(winter.dx32[,3:ncol(winter.dx32)])
# winter.dx32 <- winter.dx32 %>% 
#   filter(Sum_mo==180)
# winter.prcp$NAME==winter.dt32$NAME
# 
# # Join all dataframes and keep stations have complete sets of data
# all <- winter.prcp[,c(1:2,33)] %>% 
#   full_join(winter.dt32[,c(1:2,33)], by = "NAME") %>% 
#   full_join(winter.dx32[,c(1:2,33)], by = "NAME") %>% 
#   filter(Sum_mo==180 & Sum_mo.x==180 & Sum_mo.y==180)

stations <- c("MACON MIDDLE GA REGIONAL AIRPORT, GA US", "ORLANDO INTERNATIONAL AIRPORT, FL US", "COLUMBIA METROPOLITAN AIRPORT, SC US", "SALISBURY 9 WNW, NC US", "MURFREESBORO 5 N, TN US", "BIRMINGHAM AIRPORT, AL US", "JACKSON INTERNATIONAL AIRPORT, MS US", "LITTLE ROCK AIRPORT ADAMS FIELD, AR US", "ALEXANDRIA, LA US","ABILENE REGIONAL AIRPORT, TX US","LEXINGTON BLUEGRASS AIRPORT, KY US", "ROLLA MISSOURI S AND T, MO US", "ALBUQUERQUE INTERNATIONAL AIRPORT, NM US", "PHOENIX AIRPORT, AZ US","OKLAHOMA CITY WILL ROGERS WORLD AIRPORT, OK US", "HUTCHINSON 10 SW, KS US")



winter.new <- dat %>% 
              filter(NAME%in%stations & MONTH%in%w.months) %>% 
              mutate_at("YEAR", as.numeric) %>% 
              mutate(YEAR2 = case_when(MONTH%in%c("01") ~ YEAR-1,
                                       TRUE ~ YEAR)) %>% 
              filter(!(YEAR2%in%c(1990, 2021)))

avg <- winter.new %>% 
          group_by(YEAR2) %>% 
          summarise(avg_prcp = mean(PRCP), avg_dt32 = mean(DT32), avg_dx32 = mean(DX32))

# Palmer drought severity index

pdsi <- read.csv("raw_dat/environmental_covariates/PDSI.csv")
pdsi$Year <- str_sub(pdsi$ID, -4, -1)
pdsi$Year <- as.numeric(pdsi$Year)

years <- c(1991:2021)
pdsi.new <- pdsi %>% 
  select(State, Year, January, September:December) %>% 
  filter(Year%in%years & State %in%states) %>% 
  group_by(State) %>% 
  mutate(January = lead(January, order_by = State))

mean.pdsi <- pdsi.new %>% 
                group_by(Year) %>% 
                summarise(across(where(is.numeric), mean))

mean.pdsi <- cbind(mean.pdsi, Mean.pdsi = rowMeans(mean.pdsi[,2:6]))

mean.pdsi <- mean.pdsi[complete.cases(mean.pdsi),]
#mean.pdsi <- mean.pdsi[1:nrow(mean.pdsi)-1,]

# NAO
nao <- read.csv("raw_dat/environmental_covariates/NAO.csv")

# subset from 1991-2021
nao <- nao %>% 
        filter(Year%in%years)

nao.jan <- nao %>% 
  select(January)

nao.sepdec <- nao %>% 
  select(Year,September:December)

nao.new <- cbind(nao.sepdec, lead(nao.jan))

mean.nao <- cbind(nao.new, Mean.nao = rowMeans(nao.new[,2:6]))

mean.nao <- mean.nao[complete.cases(mean.nao),]
#mean.nao <- mean.nao[1:nrow(mean.nao)-1,]


# Midlatitude 
midlat <- rbind(chicago, nplatte, desmoines, cleveland)

midlat$YEAR <- str_sub(midlat$DATE, 1, 4)
midlat$MONTH <- str_sub(midlat$DATE, 6, 7)
midlat$STATE <- str_sub(midlat$NAME, -5, -4)
midlat$YEAR <- as.numeric(midlat$YEAR)
midlat <- midlat %>%
  select(STATION, NAME, DT32, DX32, PRCP, SNOW, YEAR, MONTH, STATE) %>%
  filter(MONTH%in%c("01", "09", "10", "11", "12"))


midlat.new <- midlat %>% 
  filter(MONTH%in%w.months) %>% 
  mutate_at("YEAR", as.numeric) %>% 
  mutate(YEAR2 = case_when(MONTH%in%c("01") ~ YEAR-1,
                           TRUE ~ YEAR)) %>% 
  filter(!(YEAR2%in%c(1989, 1990, 2021, 2022)))

midlat %>% filter(if_any(everything(), is.na)) # Only one row contains NA: Chicago Midway Airport snow record September 2012

midlat.avg <- midlat.new %>% 
  group_by(YEAR2) %>% 
  summarise(avg_prcp_ml = mean(PRCP), avg_dt32_ml = mean(DT32), avg_dx32_ml = mean(DX32), avg_snow_ml = mean(SNOW, na.rm = TRUE))

# Putting everything in one table
new.dat <- as.data.frame(cbind(avg, Mean.nao = mean.nao$Mean.nao, Mean.pdsi = mean.pdsi$Mean.pdsi, midlat.avg[2:5]))

saveRDS(new.dat, "data/envi_cov.rda")

# May pond counts (estimates in thousands 1990-2019)
ponds <- c(3508.5, 3200, 3608.9, 3611.7, 5984.8, 6335.4, 7482.2, 7458.2, 4586.9, 6704.3, 3946.9, 4640.4, 2720.0, 5190.1, 3919.6, 5381.2, 6093.9, 7002.7, 4431.4, 6434.0, 6665.0, 8132.2, 5544.0, 6891.7, 7181.2, 6307.7, 5012.5, 6096.0, 5227.4, 4990.3)

# Standardize ponds
ponds.std <- (ponds-mean(ponds))/sd(ponds)

# correlation between envi covs (without ponds)
ggpairs(new.dat[,2:ncol(new.dat)])

# correlation with ponds (off set by one year because want to see whether winter conditions from previous year are correlated with pond counts of following year)
new.dat2 <- cbind(new.dat, ponds = c(ponds.std[2:length(ponds.std)], NA))
ggpairs(new.dat2)
