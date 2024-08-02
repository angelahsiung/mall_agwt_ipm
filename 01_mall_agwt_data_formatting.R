##############################################################
### The purpose of this script is to format the band-recovery, 
### harvested wings, and Bpop data for the IPM. The code in this
### script is for formatting both mallard and green-winged teal data. 
### To run code for each species, comment out the code for the other.
### Date last updated: 07/26/2024
##############################################################

# Prepare packages
list.of.packages <- c("tidyverse", "here", "todor")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) {
  install.packages(new.packages)
}
lapply(list.of.packages, require, character.only = TRUE)

# Run code for review
todor(
  todo_types = NULL,
  search_path = getwd(),
  # file = "mall_agwt_data_formatting.R",
  output = "markers"
)

# Read in banding (and recovery) data
# Mallard
# mall.bandings <- read.csv("raw_dat/GameBirds_MALL_Bandings_1991_2021.csv")
# mall.recov <- read.csv("raw_dat/GameBirds_MALL_B1991_R1991_2021.csv") ## NOTE: May need to extract this file first before loading it into R. Had to compress it before uploading it on Github because it was too big

# Green-winged teal
agwt.recov <- read.csv("raw_dat/GameBirds_AGWT_B1991_R1991_2021.csv")
agwt.bandings <- read.csv("raw_dat/GameBirds_AGWT_1991_2021_bandings.csv")

#######################################
####### Recovery data summaries #######
#######################################

# Only using status 3 birds
# recov <- mall.recov[mall.recov$Status == 3, ]
recov <- agwt.recov[agwt.recov$Status==3,]

# Subset birds recovered before 2021
recov <- recov[recov$R.Year < 2021, ]

# Get rid of events where month is not 1-12
recov <- recov[recov$R.Month < 13, ]
recov <- recov[recov$B.Month < 13, ]

# Only retain birds that were "shot"
recov <- recov[recov$How.Obt == 1, ]

# Get rid of unknown hunting seasons survived
recov <- recov[recov$Hunt..Season.Surv. != 99, ]
# unknown.hunt <- recov[recov$Hunt..Season.Surv.==99,]

# Look at age-cutoff
unique(recov$MIN_AGE_AT_ENC)

age.table <- as.data.frame(with(recov, table(MIN_AGE_AT_ENC)))

ggplot(recov, aes(x = MIN_AGE_AT_ENC)) +
  geom_histogram(binwidth = 1)

## It looks like, based on reported minimum age at encounter, number of individuals beyond age 15 have count of <100.

# Subset only individuals with minimum age of 0-15 at encounter
recov <- recov[recov$MIN_AGE_AT_ENC < 16, ]

# Excluding birds recovered outside of Flyways 1-6
recov <- recov[recov$R.Flyway %in% c(1:6), ]

# Subset only birds banded in February, March, July, August, and September
recov <- recov[recov$B.Month %in% c(2:3, 7:9), ]

# Subset only recoveries during hunting season (Sep-Jan)
recov <- recov[recov$R.Month %in% c(1, 9:12), ]

# Subset data from Atlantic, Mississippi and Central flyways, and corresponding provinces in CA
## REVIEW: Does this subsetting make sense? Would this be excluding birds banded in Canada and recovered in the US and vise versa?
recov.US <- recov[recov$R.Flyway %in% c(1:3) & recov$B.Flyway %in% c(1:3), ]
recov.CA <- recov[recov$R.Flyway == 6 & recov$RRegion..STA %in% c("NS", "NB", "PE", "LAB", "PQ", "ONT", "MAN", "SK", "AB", "NU", "NT") &
  recov$B.Flyway == 6 & recov$BRegion..STA %in% c("NS", "NB", "PE", "LAB", "PQ", "ONT", "MAN", "SK", "AB", "NU", "NT"), ] # Nova Scotia, New Brunswick, Prince Edward Island, Newfoundland and Labrador, Quebec, Ontario, Manitoba, Alberta, Saskatchewan, Nunavut, Northwest territories

recov <- rbind(recov.US, recov.CA)

# Get rid of unknown age at banding
recov <- recov[recov$Age..VAGE != "Unknown", ]


#####################################################
######## Starting new data frame for analysis #######
#####################################################

# Code adapted from Saunders et al. (2018). Citation below.
## Saunders, S. P., M. T. Farr, A. D. Wright, C. A. Bahlai, J. W. Ribeiro, S. Rossman, A. L. Sussman, T. W. Arnold, and E. F. Zipkin. 2019. “Disentangling Data Discrepancies with Integrated Population Models.” Ecology 100: 1–14.

# Bring in B.Month, convert to season
clean <- as.data.frame(matrix(NA, nrow = length(recov$B.Month), ncol = 1))

clean[recov$B.Month %in% c(7:9), ] <- 2 # pre-season banding (Jul-Sep)
clean[recov$B.Month %in% c(2:3), ] <- 1 # post-season banding (Feb-Mar)

# Bring in B.Year
clean[, 2] <- recov$B.Year

# Bring in recovery year and account for recoveries occurring in Jan
clean[, 3] <- NA
clean[recov$R.Month >= 2, 3] <- recov[recov$R.Month >= 2, "R.Year"]

# recoveries in January belong in the hunting season of previous year
clean[recov$R.Month < 2, 3] <- recov[recov$R.Month < 2, "R.Year"] - 1

# Bring in age
## REVIEW: Does age assignment makes sense based on the age in teh data and the month the bird was banded
clean[, 4] <- NA
clean[recov$Age..VAGE == "Local" | recov$Age..VAGE == "Hatch Year", 4] <- 1 # 1 indicates juvenile
clean[recov$Age..VAGE == "After Hatch Year" & recov$B.Month %in% c(2:3), 4] <- 1 # AHY banded post-season is considered juvenile
clean[recov$Age..VAGE == "Second Year", 4] <- 1
clean[recov$Age..VAGE == "Juvenile (obsolete)", 4] <- 1
clean[recov$Age..VAGE == "After Hatch Year" & recov$B.Month %in% c(7:9), 4] <- 2 # 2 indicates adult
clean[recov$Age..VAGE == "After Second Year" | recov$Age..VAGE == "After Third Year", 4] <- 2

# remove unknowns
clean <- clean[!is.na(clean$V4), ]

# Locals that are banded in Jul-Sep
recov[recov$Age..VAGE == "Local" & recov$B.Month %in% c(7:9), ]

# Locals banded in Feb-Mar (there shouldn't be any)
recov[recov$Age..VAGE == "Local" & recov$B.Month %in% c(2:3), ] # none

# Bring in sex and convert to sex-age class (1 - Juvenile male, 2 - Juvenile female, 3 - Adult male, 4 - adult female)
clean[, 5] <- NA
clean[recov$Sex..VSEX == "Male" & clean[, 4] == 1, 5] <- 1
clean[recov$Sex..VSEX == "Female" & clean[, 4] == 1, 5] <- 2
clean[recov$Sex..VSEX == "Male" & clean[, 4] == 2, 5] <- 3
clean[recov$Sex..VSEX == "Female" & clean[, 4] == 2, 5] <- 4

# Remove unknowns
clean <- clean[!is.na(clean$V5), ]

# Plot recoveries by class
clean_age_class <- clean %>%
  group_by(V3, V5) %>%
  summarize(Number_recov = n())

clean_age_class$V5 <- as.factor(clean_age_class$V5)

ggplot(clean_age_class, aes(x = V3, y = Number_recov, fill = V5)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5))

# dummy column for marray for summing later
clean[, ncol(clean) + 1] <- 1

colnames(clean) <- c("BSeason", "BYear", "RYear", "Age", "Class", "Dummy")


##################################
######## Creating m-array ########
##################################
## REVIEW: Make sure indexing is correct in the marray

Year <- sort(unique(clean$BYear))
nYear <- length(Year)
Season <- sort(unique(clean$BSeason))
nSeason <- length(Season)
Class <- unique(clean$Class)
Class <- sort(Class)
nClass <- length(Class)

marray <- array(NA,
  dim = c(nYear, nYear, nSeason, nClass),
  dimnames = list(
    Year, Year, c("post-season", "pre-season"),
    c("Juvenile_Male", "Juvenile_Female", "Adult_Male", "Adult_Female")
  )
)


for (s in 1:nSeason) {
  for (c in 1:nClass) {
    for (b in 1:nYear) {
      for (r in 1:nYear) {
        marray[b, r, s, c] <- sum(clean[clean$BYear == Year[b] & clean$RYear == Year[r] & clean$BSeason == Season[s] & clean$Class == Class[c], ncol(clean)])
      }
    }
  }
}

# juvenile m-arrays
View(marray[1:nYear, 1:nYear, 1, 1]) # juvenile male banded in Feb-Mar and recovered in subsequent hunting season
View(marray[1:nYear, 1:nYear, 2, 1]) # juvenile male banded in Jul-Sep
View(marray[1:nYear, 1:nYear, 1, 2]) # juvenile female banded in Feb-Mar
View(marray[1:nYear, 1:nYear, 2, 2]) # juvenile female banded in Jul-Sep

# adult m-arrays
View(marray[1:nYear, 1:nYear, 1, 3]) # adult males banded in Feb-Mar
View(marray[1:nYear, 1:nYear, 2, 3]) # adult males banded in Jul-Sep
View(marray[1:nYear, 1:nYear, 1, 4]) # adult female banded in Feb-Mar
View(marray[1:nYear, 1:nYear, 2, 4]) # adult female banded in Jul-Sep

########################################
####### Banding data summaries #########
########################################

# Only using status 3 birds
# band <- mall.bandings[mall.bandings$Status == 3, ]
band <- agwt.bandings[agwt.bandings$Status==3,]

# Subset birds banded before 2021
band <- band[band$B.Year < 2021, ]

# Subset B.Month for both banding seasons

band <- band[band$B.Month %in% c(2:3, 7:9), ]

# Subset data from Atlantic and Mississippi flyways, and corresponding provinces in CA

band.US <- band[band$B.Flyway %in% c(1:3), ]
band.CA <- band[band$B.Flyway == 6 & band$Region..State %in% c("Quebec", "Nova Scotia", "New Brunswick", "Newfoundland and Labrador and St. Pierre et Miquelon", "Ontario", "Manitoba", "Prince Edward Island", "Saskatchewan", "Alberta", "Nunavut", "Northwest Territories"), ]

band <- rbind(band.US, band.CA)

# Get rid of unknown age at banding
band <- band[band$Age..VAGE != "Unknown", ]

## clean data frame, bring in B.Month
clean.band <- as.data.frame(matrix(NA, nrow = length(band$B.Month), ncol = 1))

clean.band[band$B.Month %in% c(7:9), ] <- 2 # pre-season banding (Jul-Sep)
clean.band[band$B.Month %in% c(2:3), ] <- 1 # post-season banding (Feb-Mar)

# Bring in B.Year
clean.band[, 2] <- band$B.Year

# Bring in age
clean.band[, 3] <- NA
clean.band[band$Age..VAGE == "Local" | band$Age..VAGE == "Hatch Year", 3] <- 1 # 1 indicates juvenile
clean.band[band$Age..VAGE == "After Hatch Year" & band$B.Month %in% c(2:3), 3] <- 1 # AHY banded in Feb and Mar are considered juveniles
clean.band[band$Age..VAGE == "Second Year", 3] <- 1
clean.band[band$Age..VAGE == "After Hatch Year" & band$B.Month %in% c(7:9), 3] <- 2 # 2 indicates adult # AHY banded in Jul-Sep are considered adults
clean.band[band$Age..VAGE == "After Second Year" | band$Age..VAGE == "After Third Year", 3] <- 2

# Bring in sex
clean.band[band$Sex..VSEX == "Male" & clean.band[, 3] == 1, 4] <- 1 # juvenile males
clean.band[band$Sex..VSEX == "Female" & clean.band[, 3] == 1, 4] <- 2 # juvenile females
clean.band[band$Sex..VSEX == "Male" & clean.band[, 3] == 2, 4] <- 3 # adult males
clean.band[band$Sex..VSEX == "Female" & clean.band[, 3] == 2, 4] <- 4 # adult females

# Remove unknown sex
band <- band[!is.na(clean.band$V4), ]
clean.band <- clean.band[!is.na(clean.band$V4), ]

# Add dummy column

clean.band[, 5] <- as.numeric(band$Count.of.Birds)
colnames(clean.band) <- c("BSeason", "BYear", "Age", "Class", "Dummy")

# Sorting
Year.bands <- sort(unique(clean.band$BYear))
nYear.bands <- length(Year.bands)
Season.bands <- sort(unique(clean.band$BSeason))
nSeason.bands <- length(Season.bands)
Class.bands <- sort(unique(clean.band$Class))
nClass.bands <- length(Class.bands)

# Total number of bandings each year

bands <- array(NA,
  dim = c(nYear, 1, nSeason, nClass),
  dimnames = list(
    Year, "Banded", c("post-season", "pre-season"),
    c("Juvenile_Male", "Juvenile_Female", "Adult_Male", "Adult_Female")
  )
)

for (s in 1:nSeason.bands) {
  for (c in 1:nClass.bands) {
    for (b in 1:nYear.bands) {
      bands[b, , s, c] <- sum(clean.band[clean.band$BYear == Year.bands[b] & clean.band$BSeason == Season.bands[s] & clean.band$Class == Class.bands[c], ncol(clean.band)])
    }
  }
}

# Total bands not recovered each year
nonrecov <- array(NA,
  dim = c(nYear, 1, nSeason, nClass),
  dimnames = list(Year, "Not Recovered", c("post-season", "pre-season"), c("Juvenile_Male", "Juvenile_Female", "Adult_Male", "Adult_Female"))
)

for (s in 1:nSeason.bands) {
  for (c in 1:nClass.bands) {
    for (b in 1:nYear.bands) {
      nonrecov[b, , s, c] <- bands[b, , s, c] - sum(marray[b, , s, c])
    }
  }
}

# Incorporate nonrecovery in marray
MALL.marray <- array(NA,
  dim = c(nYear, nYear + 1, nSeason, nClass),
  dimnames = list(
    Year, c(Year, "Nonrecov"), c("post-season", "pre-season"),
    c("Juvenile_Male", "Juvenile_Female", "Adult_Male", "Adult_Female")
  )
)
MALL.release <- array(NA,
  dim = c(nYear, nSeason, nClass),
  dimnames = list(
    Year, c("post-season", "pre-season"),
    c("Juvenile_Male", "Juvenile_Female", "Adult_Male", "Adult_Female")
  )
)

for (s in 1:nSeason) {
  for (c in 1:nClass) {
    for (b in 1:nYear) {
      for (r in 1:(nYear + 1)) {
        if (r <= nYear) {
          MALL.marray[b, r, s, c] <- marray[b, r, s, c]
        } else {
          MALL.marray[b, r, s, c] <- nonrecov[b, , s, c]
        }
        MALL.release[b, s, c] <- sum(MALL.marray[b, , s, c])
      }
    }
  }
}

saveRDS(MALL.marray, file = "data/MALL_marray_2.rda")
saveRDS(MALL.release, file = "data/MALL_release_2.rda")

## Green-winged teal
# AGWT.marray<-array(NA,dim=c(nYear,nYear+1,nSeason,nClass),
#                    dimnames =list(Year, c(Year,"Nonrecov"), c("post-season", "pre-season"),
#                                   c("Juvenile_Male", "Juvenile_Female", "Adult_Male", "Adult_Female")))
# AGWT.release<-array(NA,dim=c(nYear,nSeason,nClass),
#                     dimnames =list(Year, c("post-season", "pre-season"),
#                                    c("Juvenile_Male", "Juvenile_Female", "Adult_Male", "Adult_Female")))
#
# for (s in 1:nSeason){
#   for (c in 1:nClass){
#     for (b in 1:nYear){
#       for (r in 1:(nYear+1)){
#         if(r <= nYear){
#           AGWT.marray[b,r,s,c] <- marray[b,r,s,c]
#         }else{
#           AGWT.marray[b,r,s,c] <- nonrecov[b,,s,c]
#         }
#         AGWT.release[b,s,c] <- sum(AGWT.marray[b,,s,c])
#       }}}}
#
#  saveRDS(AGWT.marray, file = "AGWT_marray_new.rda")
#  saveRDS(AGWT.release, file = "AGWT_release_new.rda")



## Additional data summaries
# # Bandings by month
#
# band_by_mo <- band %>%
#   filter(!is.na(Count.of.Birds)) %>%
#   group_by(B.Month, Region..Flyway) %>%
#   summarize(Number_band = sum(Count.of.Birds))
#
#
# ## Plotting number of birds banded by month and by flyway
# ggplot(band_by_mo, aes(x = B.Month, y = Number_band, fill = Region..Flyway)) +
#   geom_col() +
#   scale_x_discrete(name ="Month banded",
#                    limits= factor(c(1:12))) +
#   ggtitle("Number of total Mallards banded by month in 1960-2021") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_manual(values=cbbPalette)+
#   guides(fill=guide_legend(title="Flyway"))
#
#
#
# # Bandings by year
#
# band_by_yr <- band %>%
#   filter(!is.na(Count.of.Birds)) %>%
#   group_by(B.Year, Region..Flyway) %>%
#   summarize(Number_band = sum(Count.of.Birds))
#
# ggplot(band_by_yr, aes(x = B.Year, y = Number_band, fill = Region..Flyway)) +
#   geom_col() +
#   scale_x_discrete(name ="Year banded",
#                    limits= c(1991:2021)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle("Number of total green-winged teals banded by year in 1991-2021") +
#   theme(plot.title = element_text(size = 16, hjust = 0.5),
#         legend.position = "bottom") +
#   scale_fill_manual(values=cbbPalette)+
#   guides(fill=guide_legend(title="Flyway"))
#
# band_by_sex <- band %>%
#   filter(!is.na(Count.of.Birds)) %>%
#   group_by(Sex..VSEX, B.Year) %>%
#   summarize(Number_band = sum(Count.of.Birds))
#
# ggplot(band_by_sex, aes(x = B.Year, y = Number_band, fill = Sex..VSEX)) +
#   geom_col() +
#   scale_x_discrete(name ="Year banded",
#                    limits= c(1991:2021)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle("Number of total green-winged teals banded grouped by sex in 1991-2021") +
#   theme(plot.title = element_text(size = 16, hjust = 0.5)) +
#   scale_fill_manual(values=cbbPalette)+
#   guides(fill=guide_legend(title="Sex"))
#
# band_by_age <- band %>%
#   filter(!is.na(Count.of.Birds)) %>%
#   group_by(Age..VAGE, B.Year) %>%
#   summarize(Number_band = sum(Count.of.Birds))
#
# ggplot(band_by_age, aes(x = B.Year, y = Number_band, fill = Age..VAGE)) +
#   geom_col() +
#   scale_x_discrete(name ="Year banded",
#                    limits= c(1991:2021)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle("Number of total green-winged teals banded grouped by age in 1991-2021") +
#   theme(plot.title = element_text(size = 16, hjust = 0.5),
#         legend.position = "bottom") +
#   scale_fill_manual(values=cbbPalette)+
#   guides(fill=guide_legend(title="Sex"))
# # Birds banded in January per year
#
# jan.band <- band[band$B.Month==1, ]
#
# jan.band_by_yr <- jan.band %>%
#   filter(!is.na(Count.of.Birds)) %>%
#   group_by(B.Year, Region..Flyway) %>%
#   summarize(Number_band = sum(Count.of.Birds))
#
# ggplot(jan.band_by_yr, aes(x = B.Year, y = Number_band, fill = Region..Flyway)) +
#   geom_col() +
#   scale_x_discrete(name ="Year banded",
#                    limits= c(1960:2021)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle("Number of Mallards banded in January by year in 1960-2021") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_manual(values=cbbPalette)+
#   guides(fill=guide_legend(title="Flyway"))


#################################
####### Format wing data ########
#################################
# import raw wing data
# mall.wings <- read.csv("raw_dat/PCS_MALL_1991_2020.csv")
agwt.wings <- read.csv("raw_dat/PCS_AGWT_1991_2020.csv")

# Remove unknown age or sex
# mall.wings <- mall.wings[!(mall.wings$age_code==0 | mall.wings$sex_code==0),] # 889384 records remaining
agwt.wings <- agwt.wings[!(agwt.wings$age_code == 0 | agwt.wings$sex_code == 0), ] # 889384 records remaining


# Exclude states in the Pacific Flyway and Alaska
# mall.wings$state <- trimws(mall.wings$state) # get rid of extra space at the end of states
# mall.wings <- mall.wings[!(mall.wings$state%in%c("AK", "CA", "OR", "WA", "UT", "ID", "NV", "AZ")),]
agwt.wings$state <- trimws(agwt.wings$state) # get rid of extra space at the end of states
agwt.wings <- agwt.wings[!(agwt.wings$state %in% c("AK", "CA", "OR", "WA", "UT", "ID", "NV", "AZ")), ]

# Summarize wings by cohort
# cohort_table <- mall.wings %>%
#   group_by(cohort, Season) %>%
#   summarize(Number = n())

cohort_table <- agwt.wings %>%
  group_by(cohort, Season) %>%
  summarize(Number = n())

jm.wing <- cohort_table[cohort_table$cohort == "Immature Male", ]
am.wing <- cohort_table[cohort_table$cohort == "Adult Male", ]
jf.wing <- cohort_table[cohort_table$cohort == "Immature Female", ]
af.wing <- cohort_table[cohort_table$cohort == "Adult Female", ]

male.wing <- cohort_table %>%
  filter(cohort %in% c("Immature Male", "Adult Male")) %>%
  group_by(Season) %>%
  summarise(sum = sum(Number))
female.wing <- cohort_table %>%
  filter(cohort %in% c("Immature Female", "Adult Female")) %>%
  group_by(Season) %>%
  summarise(sum = sum(Number))

saveRDS(jf.wing, file = "data/AGWT_juv_female_wing_new.rda")
saveRDS(female.wing, file = "data/AGWT_all_female_wing_new.rda")

# plot cohorts

ggplot(cohort_table, aes(x = cohort, y = Number)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
  ggtitle("Number wings by cohort") +
  theme(plot.title = element_text(hjust = 0.5))

################################
###### Bpop survey #############
###############################

# Breeding population estimates and standard errors (in thousands) for Mallards and Green-winged teal from the traditional survey area (strata 1-18, 20-50, 75-77) and eastern survey area.

# Read in stratum-specific estimates
## mallard
dat.tsa <- read.csv("wbphs_traditionalarea_estimates_forDistribution.csv")
dat.esa <- read.csv("easternsurvey_mall_agwt.csv")
dat.esa <- dat.esa[dat.esa$MASAlpha == "mall" & dat.esa$ReportingScale == "EasternCA", ]
dat.tsa <- dat.tsa[dat.tsa$survey_species == "MALL" & dat.tsa$stratum %in% c(14:max(dat.tsa$stratum)), ]

## agwt
# dat.tsa <- read.csv("wbphs_traditionalarea_estimates_forDistribution.csv")
# dat.esa <- read.csv("easternsurvey_mall_agwt.csv")
# dat.esa <- dat.esa[dat.esa$MASAlpha=="agwt" & dat.esa$ReportingScale == "EasternCA",]
# dat.tsa <- dat.tsa[dat.tsa$survey_species=="AGWT" & dat.tsa$stratum%in%c(14:max(dat.tsa$stratum)),]


dat.tsa <- dat.tsa %>%
  group_by(survey_year) %>%
  summarise(Sum = sum(estimate))

# remove 2022 estimate
dat.tsa <- dat.tsa[-nrow(dat.tsa), ]
dat.esa <- dat.esa[-nrow(dat.esa), ]

# divide by 10,000
dat.tsa$Sum <- dat.tsa$Sum / 10000
y_esa <- y_esa / 10000


y_esa <- c(rep(NA, 1998 - 1955), dat.esa$Est)

y_tsa <- dat.tsa$Sum
year.trad <- c(1955:2019)
T <- length(year.trad)

# subsetting data from 1991-2021
y_esa <- c(y_esa[37:length(y_esa) - 1], NA, NA) # adding NAs to years with no survey
y_tsa <- c(y_tsa[37:length(y_tsa) - 1], NA, NA)
saveRDS(y_tsa, "MALL_TSA_Bpop.rda")
saveRDS(y_esa, "MALL_ESA_Bpop.rda")
# saveRDS(y_tsa, "AGWT_TSA_Bpop.rda")
# saveRDS(y_esa, "AGWT_ESA_Bpop.rda")


####################
###### Misc. #######
####################

# Creating table with flyway and states within each flyway
atl.states <- unique(mall.recov.all[mall.recov.all$RFly..Flyway == "Atlantic", ]$RRegion..State)
miss.states <- unique(mall.recov.all[mall.recov.all$RFly..Flyway == "Mississippi", ]$RRegion..State)
cen.states <- unique(mall.recov.all[mall.recov.all$RFly..Flyway == "Central", ]$RRegion..State)
pac.states <- unique(mall.recov.all[mall.recov.all$RFly..Flyway == "Pacific", ]$RRegion..State)
ca.prov <- unique(mall.recov.all[mall.recov.all$RFly..Flyway == "Canada", ]$RRegion..State)

flyway.states <- as.data.frame(cbind(Flyway = c(rep("Atlantic", length(atl.states)), rep("Mississippi", length(miss.states)), rep("Central", length(cen.states)), rep("Pacific", length(pac.states)), rep("Canada", length(ca.prov))), States = c(atl.states, miss.states, cen.states, pac.states, ca.prov)))

write.csv(flyway.states, "Flyway_States.csv")
