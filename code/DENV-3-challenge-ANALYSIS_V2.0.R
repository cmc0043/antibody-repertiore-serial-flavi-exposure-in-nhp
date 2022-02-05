#### Code Overview ####

## Viral Loads
# Authors: Chelsea Crooks & Hunter Ries & Luis Haddock

## Overview
# - csv file downloaded from ZEST open research portal
# - clean csv data 
# - Save cleaned dataframes as CSV files
# - Plot the dataframes as line plots 
# - X: DPI
# - Y: viral load (copies DENV-3 vRNA/mL plasma)  

## Reference 
# Animal IDs: 
# - Cohort 1 (flavi-naive): 043-101, 043-102, 043-103, 043-104, 043-105, 043-106, 043-107, 043-108
# - Cohort 2 (ZIKV-exposed): 044-109, 034-102, 046-101, 044-110, 044-114, 036-102, 044-116, 044-118, 046-103, 044-103, 044-104, 046-104, 035-110, 035-109
# - Cohort 3 (DENV-2/ZIKV-exposed): 042-101, 042-102, 042-104, 042-106, 042-107
# Study on ZEST Open Research portal: ZIKV-043 (https://openresearch.labkey.com/project/ZEST/OConnor/ZIKV-043/begin.view?)


## Challenge data
# - DENV-3 Isolate: Dengue virus/H. sapiens-tc/NIC/DENV3-6629
# - DENV-3 Sequence: SRA link to come
# - Challenge dose: 10^4 PFU
# - Challenge route: Subcutaneous
# - Days viral loads performed: 0-10, 15 days post-infection


## Input: 
# 1. Raw viral load data from ZIKV portal (see "README-viral-load-raw-data-search-parameters.txt"): 
#   - https://openresearch.labkey.com/ > ZEST > Private > Viral loads in EHR > select Animal IDs above + DENV-3 assay
# - `C*-VL-raw.csv`

## Output: 
# 1. CSV files with all viral loads:       
#   - `C*-[virus]-VL-clean.csv`

# 2. Figures showing viral loads over time:   
#   - ``
  
#### Session prep ####
## clear Global Environment
rm(list = ls())

## install packages and load libraries as required
if(!require(tidyverse)){
  install.packages("tidyverse",dependencies = T)
  library(tidyverse)
}
if(!require(ggplot2)){
  install.packages("ggplot2",dependencies = T)
  library(ggplot2)
}
if(!require(MESS)){
  install.packages("MESS",dependencies = T)
  library(MESS)
}
if(!require(gridExtra)){
  install.packages("gridExtra",dependencies = T)
  library(gridExtra)
}
if(!require(grid)){
  install.packages("grid",dependencies = T)
  library(grid)
}
if(!require(lubridate)){
  install.packages("lubridate",dependencies = T)
  library(lubridate)
}
if(!require(ggsignif)){
  install.packages("ggsignif",dependencies = T)
  library(ggsignif)
}
if(!require(data.table)){
  install.packages("data.table",dependencies = T)
  library(data.table)
}

#### Import + set variables ####
## set WD
setwd("~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_raw/viral-loads")

## import data

# list number of cohorts in analysis, cohort names, and cohort datafiles
cohort_number <- (1:3)
number_of_cohorts <- length(cohort_number)

nameC1 <- "Flavivirus-naive"
nameC2 <- "ZIKV-exposed"
nameC3 <- "DENV-2/ZIKV-exposed"
  # nameCn <- "xxx"

raw_data_C1 <- read_csv("C1-VL-raw.csv")
raw_data_C2 <- read_csv("C2-VL-raw.csv")
raw_data_C3 <- read_csv("C3-VL-raw.csv")
  # raw_data_Cn <- read_csv("xxx.csv")

# range of dpi values to include in plots and stats analysis
dpimin <- as.numeric(0)
dpimax <- as.numeric(15)

# sample source
samplesource = "Plasma"

# specify virus and viral load assay
virus = "DENV-3"
VLassay = "DENV3" # options: "DENV3", "Lanciotti_ZIKV_universal", "DENV2"

## import animal data
# CSV should be at least a three column CSV with "animal", "hexcode", "challenge-date-a", and as column headers 
# if animal has received multiple challenges, additional columns can be added using "challenge-date-b", "-c", "-d", etc.
# CSV should have one row for each animal on the study  
# leave cells blank if animal wasn't challenged

animal_study_data <- read_csv("043-VL-study-data.csv")

## write out data file names

# Cleaned df C1 (file name should be cohort-assay-VL-clean)
cleanC1csv = "~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned/C1-DENV-3-VL-clean.csv"

# Cleaned df C2
cleanC2csv = "~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned/C2-DENV-3-VL-clean.csv"

# Cleaned df C3
cleanC3csv = "~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned/C3-DENV-3-VL-clean.csv"

# Path for storing plots
plotpng = "~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/figures/VL/"

#### Clean animal data ####

## clean challenge information  

# make sure challenge dates are dates
challenge_dates <- as.data.frame(animal_study_data)
challenge_dates$dateA <- mdy(challenge_dates$dateA)
challenge_dates$dateB <- mdy(challenge_dates$dateB)
challenge_dates$dateC <- mdy(challenge_dates$dateC)

# A dates
datesA <- distinct(challenge_dates, animal, dateA)
chaldateA <- challenge_dates$dateA
names(chaldateA) <- challenge_dates$animal

# B dates
datesB <- distinct(challenge_dates, animal, dateB)
chaldateB <- challenge_dates$dateB
names(chaldateB) <- challenge_dates$animal

# C dates
datesC <- distinct(challenge_dates, animal, dateC)
chaldateC <- challenge_dates$dateC
names(chaldateC) <- challenge_dates$animal


## clean color scheme information
color_scheme <- as.data.frame(animal_study_data)
colors <- distinct(color_scheme, animal, hexcode)
anihex <- colors$hexcode
names(anihex) <- colors$animal


#### Clean Cohort 1 (Flavi-naive) ####
## Remove remove unneeded columns
raw_data_C1$`Key` <- NULL
raw_data_C1$`Nucleic Acid` <- NULL
raw_data_C1$`Viral Load Replicates` <- NULL
raw_data_C1$`Comment` <- NULL
raw_data_C1$`Experiment Number` <- NULL
raw_data_C1$`RNA Isolation Method` <- NULL

## rename columns
raw_data_C1 <- raw_data_C1 %>%
  rename(animal = `Participant ID`) %>%
  rename(date = `Date`) %>%
  rename(amount = `Viral Load`) %>%
  rename(assay = `Assay`) %>%
  rename(equivocal = `Equivocal`) %>%
  rename(source = `Sample Source`)

## convert to data frame
raw_data_C1 <- as.data.frame(raw_data_C1)

## factor by animal and by assay
raw_data_C1$animal <- as.factor(raw_data_C1$animal)
raw_data_C1$assay <- as.factor(raw_data_C1$assay)

## viral titers are numbers, ya know
raw_data_C1$amount <- as.numeric(raw_data_C1$amount)

## dates are dates
raw_data_C1$date <- ymd(raw_data_C1$date)

## make sure all VLs come from plasma and then remove source column
raw_data_C1 <- filter(raw_data_C1, raw_data_C1$source == samplesource)
raw_data_C1$source <- NULL

## remove all rows that are 'equivocal'

# replace empty "NA" values with 0 
is.na_replace_0_C1 <- raw_data_C1$equivocal                               
is.na_replace_0_C1[is.na(is.na_replace_0_C1)] <- 0 
raw_data_C1$equivocal <- is.na_replace_0_C1

# convert to factor 
raw_data_C1$equivocal <- as.factor(raw_data_C1$equivocal)

# select only rows that are NOT equivocal & then remove equivocal column
raw_data_C1 <- filter(raw_data_C1, raw_data_C1$equivocal == "0")
raw_data_C1$equivocal <- NULL

## select only rows that use the DENV3 assay
raw_data_C1 <- filter(raw_data_C1, raw_data_C1$assay == VLassay)

## Create df with C1 animals 
C1 <- unique(raw_data_C1[c("animal")])

## Create list with each animal per group
list_C1 <- split(raw_data_C1,raw_data_C1$animal)

## calculate dpi -- there must be a VL for day 0 in order for function to work

# create a list with an integer for each indexes in list_C1
length_list_C1 <- 1:length(list_C1)

# for loop to calculate DPI based on the challenge date and to create the DPI column
if (length(C1$animal) > 0){
for (i in length_list_C1){
  n = i
  length_C1 <- length(list_C1[[n]]$date)
  for (i in 1:length_C1) {
    list_C1[[n]]$dpi[i] <- 
      as.numeric(difftime(
        list_C1[[n]]$date[i],
        chaldateA[[as.character(unique(list_C1[[n]]$animal))]],
        units = "days"))
  }
} 

# for loop to remove any DPIs outside the specified min max range
for (i in length_list_C1){
  n = i
  #list_C1[[n]]$dpi <- as.numeric(list_C1[[n]]$dpi)
  list_C1[[n]] <- filter(list_C1[[n]], list_C1[[n]]$dpi >= dpimin)
  list_C1[[n]] <- filter(list_C1[[n]], list_C1[[n]]$dpi <= dpimax)
}

## one big, happy df & export to clean folder
df_C1 <- bind_rows(list_C1)
write.csv(df_C1, cleanC1csv, row.names = FALSE)
} else {
  no_VL <- "No VL data for this cohort"
  write.csv(no_VL, cleanC1csv, row.names = FALSE)
}

#### Clean Cohort 2 (ZIKV-exposed) ####
## Remove remove unneeded columns
raw_data_C2$`Key` <- NULL
raw_data_C2$`Nucleic Acid` <- NULL
raw_data_C2$`Viral Load Replicates` <- NULL
raw_data_C2$`Experiment Number` <- NULL
raw_data_C2$`RNA Isolation Method` <- NULL

## rename columns
raw_data_C2 <- raw_data_C2 %>%
  rename(animal = `Participant ID`) %>%
  rename(date = `Date`) %>%
  rename(amount = `Viral Load`) %>%
  rename(assay = `Assay`) %>%
  rename(equivocal = `Equivocal`) %>%
  rename(source = `Sample Source`) %>%
  rename(comment = `Comment`)

## convert to data frame
raw_data_C2 <- as.data.frame(raw_data_C2)

## factor by animal and by assay
raw_data_C2$animal <- as.factor(raw_data_C2$animal)
raw_data_C2$assay <- as.factor(raw_data_C2$assay)

## viral titers are numbers, ya know
raw_data_C2$amount <- as.numeric(raw_data_C2$amount)

## dates are dates
raw_data_C2$date <- ymd(raw_data_C2$date)

## make sure all VLs come from plasma and then remove source column
raw_data_C2 <- filter(raw_data_C2, raw_data_C2$source == "Plasma")
raw_data_C2$source <- NULL

## remove all rows that are 'equivocal'

# replace empty "NA" values with 0 
is.na_replace_0_C2 <- raw_data_C2$equivocal                               
is.na_replace_0_C2[is.na(is.na_replace_0_C2)] <- 0 
raw_data_C2$equivocal <- is.na_replace_0_C2

# convert to factor 
raw_data_C2$equivocal <- as.factor(raw_data_C2$equivocal)

# select only rows that are NOT equivocal & then remove equivocal column
raw_data_C2 <- filter(raw_data_C2, raw_data_C2$equivocal == "0")
raw_data_C2$equivocal <- NULL

## select only rows that use the DENV3 assay
raw_data_C2 <- filter(raw_data_C2, raw_data_C2$assay == "DENV3")

## remove all rows that are 'equivocal'

# replace empty "NA" values with 0 
is.na_replace_0_C2 <- raw_data_C2$comment                               
is.na_replace_0_C2[is.na(is.na_replace_0_C2)] <- 0 
raw_data_C2$comment <- is.na_replace_0_C2

# select only rows that are NOT Sodium citrate and then remove the comment column
raw_data_C2 <- filter(raw_data_C2, raw_data_C2$comment != "Sodium citrate")
raw_data_C2$comment <- NULL

## Create df with C2 animals 
C2 <- unique(raw_data_C2[c("animal")])

## Create list with each animal per group
list_C2 <- split(raw_data_C2,raw_data_C2$animal)

## calculate dpi -- there must be a VL for day 0 in order for function to work

# create a list with an integer for each indexes in list_C2
length_list_C2 <- 1:length(list_C2)

# for loop to calculate DPI and to create the DPI column
for (i in length_list_C2){
  n = i
  length_C2 <- length(list_C2[[n]]$date)
  for (i in 1:length_C2) {
    list_C2[[n]]$dpi[1] <- 
      as.numeric(difftime(
        list_C2[[n]]$date[1],
        max(list_C2[[n]]$date),
        units = "days"))
    list_C2[[n]]$dpi[i] <- 
      sum(as.numeric(difftime(
        list_C2[[n]]$date[i],
        max(list_C2[[n]]$date),
        units = "days")) +
          (-1*list_C2[[n]]$dpi[1]))
  }
  list_C2[[n]]$dpi[1]=0
}

## one big, happy df & export to clean folder
df_C2 <- bind_rows(list_C2)
write.csv(df_C2,"~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned/C2-DENV-3-VL-clean.csv", row.names = FALSE)

#### Clean Cohort 3 (DENV-2/ZIKV-exposed) ####

## Remove remove unneeded columns
raw_data_C3$`Key` <- NULL
raw_data_C3$`Nucleic Acid` <- NULL
raw_data_C3$`Viral Load Replicates` <- NULL
raw_data_C3$`Comment` <- NULL
raw_data_C3$`Experiment Number` <- NULL
raw_data_C3$`RNA Isolation Method` <- NULL

## rename columns
raw_data_C3 <- raw_data_C3 %>%
  rename(animal = `Participant ID`) %>%
  rename(date = `Date`) %>%
  rename(amount = `Viral Load`) %>%
  rename(assay = `Assay`) %>%
  rename(equivocal = `Equivocal`) %>%
  rename(source = `Sample Source`)

## convert to data frame
raw_data_C3 <- as.data.frame(raw_data_C3)

## factor by animal and by assay
raw_data_C3$animal <- as.factor(raw_data_C3$animal)
raw_data_C3$assay <- as.factor(raw_data_C3$assay)

## viral titers are numbers, ya know
raw_data_C3$amount <- as.numeric(raw_data_C3$amount)

## dates are dates
raw_data_C3$date <- ymd(raw_data_C3$date)

## make sure all VLs come from plasma and then remove source column
raw_data_C3 <- filter(raw_data_C3, raw_data_C3$source == "Plasma")
raw_data_C3$source <- NULL

## remove all rows that are 'equivocal'

# replace empty "NA" values with 0 
is.na_replace_0_C3 <- raw_data_C3$equivocal                               
is.na_replace_0_C3[is.na(is.na_replace_0_C3)] <- 0 
raw_data_C3$equivocal <- is.na_replace_0_C3

# convert to factor 
raw_data_C3$equivocal <- as.factor(raw_data_C3$equivocal)

# select only rows that are NOT equivocal & then remove equivocal column
raw_data_C3 <- filter(raw_data_C3, raw_data_C3$equivocal == "0")
raw_data_C3$equivocal <- NULL

## select only rows that use the DENV3 assay
raw_data_C3 <- filter(raw_data_C3, raw_data_C3$assay == "DENV3")

## Create df with C3 animals 
C3 <- unique(raw_data_C3[c("animal")])

## Create list with each animal per group
list_C3 <- split(raw_data_C3,raw_data_C3$animal)

## calculate dpi -- there must be a VL for day 0 in order for function to work

# create a list with an integer for each indexes in list_C3
length_list_C3 <- 1:length(list_C3)

# for loop to calculate DPI and to create the DPI column
for (i in length_list_C3){
  n = i
  length_C3 <- length(list_C3[[n]]$date)
  for (i in 1:length_C3) {
    list_C3[[n]]$dpi[1] <- 
      as.numeric(difftime(
        list_C3[[n]]$date[1],
        max(list_C3[[n]]$date),
        units = "days"))
    list_C3[[n]]$dpi[i] <- 
      sum(as.numeric(difftime(
        list_C3[[n]]$date[i],
        max(list_C3[[n]]$date),
        units = "days")) +
          (-1*list_C3[[n]]$dpi[1]))
  }
  list_C3[[n]]$dpi[1]=0
}

## one big, happy df & export to clean folder
df_C3 <- bind_rows(list_C3)
write.csv(df_C3,"~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned/C3-DENV-3-VL-clean.csv", row.names = FALSE)


#### Plot C1 ####
if (length(C1$animal) > 0){
## change 0s to 0.001 and log transform
df_C1$amount <- as.character(df_C1$amount)
df_C1$amount[df_C1$amount == "0"] <- "0.001"
df_C1$amount <- as.numeric(df_C1$amount)
df_C1$amount <- log10(df_C1$amount)
df_C1$amount <- as.character(df_C1$amount)
df_C1$amount[df_C1$amount == "-3"] <- "0"
df_C1$amount <- as.numeric(df_C1$amount)
df_C1$dpi <- as.integer(df_C1$dpi)
df_C1$amount[df_C1$amount == -5] <- NA
color_break_C1 <- as.character(C1$animal)

## Plot
Plot_C1 <- ggplot(na.omit(df_C1), aes(dpi, amount, group = animal, color = animal)) + 
  geom_line(size = .5) + 
  scale_y_continuous(
    breaks = c(2, 3, 4, 5, 6, 7, 8),
    labels = c(2, 3, 4, 5, 6, 7, 8),
    expand = c(0, 0)) + 
  scale_x_continuous(
    breaks = dpimin:dpimax, 
    labels = dpimin:dpimax,
    expand = c(0, 0)) + 
  coord_cartesian(
    ylim = c(2, 8),
    xlim = c((dpimin - 1),(dpimax + 1))) + 
  scale_color_manual(
    values = anihex, 
    breaks = color_break_C1) +
  labs(
    x = "Days post-infection",
    y = paste("log10", virus, "vRNA copies/mL", sep = " "),
    title = paste(nameC1, "cohort")) + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),     #removes ver. grid lines
    panel.grid.minor = element_blank(),     #removes hor. grid lines
    legend.position = "right",              #place legend on right
    legend.title = element_blank(),         #remove legend title
    plot.title = element_text(hjust = 0),   #move title left
    plot.margin = margin(10, 10, 10, 20),   #give us some room!
    strip.text.x = element_text(size = 10), #change text size
    strip.text.y = element_text(size = 10)) #change text size

# visualize plot in R studio
Plot_C1

# write out figure to png
ggsave(
  paste("C1", virus, "VL.png", sep="-"),
  plot = Plot_C1,
  device = NULL,
  path = plotpng,
  scale = 1,
 # width = NA,
 # height = NA,
 # units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

} else {print("No viral load data to plot")}
# Plot rendering & export found in `All Plots` below

#### Plot C2 ####

## change 0s to 0.001 and log transform
df_C2$amount <- as.character(df_C2$amount)
df_C2$amount[df_C2$amount == "0"] <- "0.001"
df_C2$amount <- as.numeric(df_C2$amount)
df_C2$amount <- log10(df_C2$amount)
df_C2$amount <- as.character(df_C2$amount)
df_C2$amount[df_C2$amount == "-3"] <- "0"
df_C2$amount <- as.numeric(df_C2$amount)
df_C2$dpi <- as.integer(df_C2$dpi)
df_C2$amount[df_C2$amount == -5] <- NA
color_break_C2 <- as.character(C2$animal)

## Plot
Plot_C2 <- ggplot(na.omit(df_C2), aes(dpi, amount, group = animal, color = animal)) + 
  geom_line(size = .5) + 
  scale_y_continuous(
    breaks = c(2, 3, 4, 5, 6, 7, 8),
    labels = c(2, 3, 4, 5, 6, 7, 8),
    expand = c(0, 0)) + 
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    expand = c(0, 0)) + 
  coord_cartesian(
    ylim = c(2, 8),
    xlim = c(-1, 16)) + 
  scale_color_manual(
    values = anihex, 
    breaks = color_break_C2) +
  labs(
    x = "Days post-infection",
    y = expression(paste("log"[10]," DENV-3 vRNA copies/mL")),
    title = "ZIKV-exposed cohort") + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),     #removes ver. grid lines
    panel.grid.minor = element_blank(),     #removes hor. grid lines
    legend.position = "right",              #place legend on right
    legend.title = element_blank(),         #remove legend title
    plot.title = element_text(hjust = 0),   #move title left
    plot.margin = margin(10, 10, 10, 20),   #give us some room!
    strip.text.x = element_text(size = 10), #change text size
    strip.text.y = element_text(size = 10)) #change text size

# Plot_C2
# Plot rendering & export found in `All Plots` below


#### Plot C3 ####

## change 0s to 0.001 and log transform
df_C3$amount <- as.character(df_C3$amount)
df_C3$amount[df_C3$amount == "0"] <- "0.001"
df_C3$amount <- as.numeric(df_C3$amount)
df_C3$amount <- log10(df_C3$amount)
df_C3$amount <- as.character(df_C3$amount)
df_C3$amount[df_C3$amount == "-3"] <- "0"
df_C3$amount <- as.numeric(df_C3$amount)
df_C3$dpi <- as.integer(df_C3$dpi)
df_C3$amount[df_C3$amount == -5] <- NA
color_break_C3 <- as.character(C3$animal)

## Plot
Plot_C3 <- ggplot(na.omit(df_C3), aes(dpi, amount, group = animal, color = animal)) + 
  geom_line(size = .5) + 
  scale_y_continuous(
    breaks = c(2, 3, 4, 5, 6, 7, 8),
    labels = c(2, 3, 4, 5, 6, 7, 8),
    expand = c(0, 0)) + 
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    expand = c(0, 0)) + 
  coord_cartesian(
    ylim = c(2, 8),
    xlim = c(-1, 16)) + 
  scale_color_manual(
    values = anihex,
    breaks = color_break_C3) + 
  labs(
    x = "Days post-infection",
    y = expression(paste("log"[10]," DENV-3 vRNA copies/mL")),
    title = "DENV-2/ZIKV-exposed cohort") + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),     #removes ver. grid lines
    panel.grid.minor = element_blank(),     #removes hor. grid lines
    legend.position = "right",              #place legend on right
    legend.title = element_blank(),         #remove legend title
    plot.title = element_text(hjust = 0),   #move title left
    plot.margin = margin(10, 10, 10, 20),   #give us some room!
    strip.text.x = element_text(size = 10), #change text size
    strip.text.y = element_text(size = 10)) #change text size

# Plot_C3
# Plot rendering & export found in `All Plots` below


#### AUC ####
## groups are already split

### calculate AUC

## truncate y-axis to 100

# write function to eliminate rows that have a VL < 100
lod_cutoff <- function(x){
  x <- x%>%
    filter(amount > 100)
    return (x)
}

# list_C1 contains df for each animal in the naive cohort    
C1_truncate_100 <- map(list_C1,
                      ~lod_cutoff(.))

# list_C2 contains df for each animal in the Zika-exposed cohort    
C2_truncate_100 <- map(list_C2,
                       ~lod_cutoff(.))

# list_C3 contains df for each animal in the DENV-2/Zika-exposed cohort    
C3_truncate_100 <- map(list_C3,
                       ~lod_cutoff(.))


## write function to calculate the AUC for each animal

## C1 - Flavi-naive
# create empty data frame for AUC values
aunderc_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "AUC")
colnames(aunderc_C1) <- x

#calculate AUC for each animal and output to empty dataframe aunderc_C1
for (i in seq_along(C1_truncate_100)){
    aunderc_C1[nrow(aunderc_C1)+1,] = c(i, auc(C1_truncate_100[[i]]$dpi, C1_truncate_100[[i]]$amount))
}

# replace the index with the animal names from data frame C1, add column with log transformed data
aunderc_C1 <- aunderc_C1 %>%
  add_column(C1,
             .after = "index")%>%
  add_column(log_AUC = log10(aunderc_C1$AUC))
  
aunderc_C1 <- aunderc_C1[-c(1)]

## C2 - ZIKV-exposed
# create empty data frame for AUC values
aunderc_C2 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "AUC")
colnames(aunderc_C2) <- x

#calculate AUC for each animal and output to empty dataframe aunderc_C2
for (i in seq_along(C2_truncate_100)){
  aunderc_C2[nrow(aunderc_C2)+1,] = c(i, auc(C2_truncate_100[[i]]$dpi, C2_truncate_100[[i]]$amount))
}

# replace the index with the animal names from data frame C2, add column with log transformed data
aunderc_C2 <- aunderc_C2 %>%
  add_column(C2,
             .after = "index")%>%
  add_column(log_AUC = log10(aunderc_C2$AUC))

aunderc_C2 <- aunderc_C2[-c(1)]

## C3 - DENV-2/ZIKV-exposed
# create empty data frame for AUC values
aunderc_C3 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "AUC")
colnames(aunderc_C3) <- x

#calculate AUC for each animal and output to empty dataframe aunderc_C3
for (i in seq_along(C3_truncate_100)){
  aunderc_C3[nrow(aunderc_C3)+1,] = c(i, auc(C3_truncate_100[[i]]$dpi, C3_truncate_100[[i]]$amount))
}

# replace the index with the animal names from data frame C3, add column with log transformed data
aunderc_C3 <- aunderc_C3 %>%
  add_column(C3,
             .after = "index")%>%
  add_column(log_AUC = log10(aunderc_C3$AUC))

aunderc_C3 <- aunderc_C3[-c(1)]

#### Peak ####

## C1 - Flavi-naive
# create empty data frame for AUC values
peak_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "peak")
colnames(peak_C1) <- x

#calculate AUC for each animal and output to empty dataframe aunderc_C1
for (i in seq_along(list_C1)){
  peak_C1[nrow(peak_C1)+1,] = c(i, max(list_C1[[i]]$amount))
}

# replace the index with the animal names from data frame C3, add column with log transformed data
peak_C1 <- peak_C1 %>%
  add_column(C1,
             .after = "index")%>%
  add_column(log_peak = log10(peak_C1$peak))

peak_C1 <- peak_C1[-c(1)]

## C2 - ZIKV-exposed
# create empty data frame for AUC values
peak_C2 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "peak")
colnames(peak_C2) <- x

#calculate AUC for each animal and output to empty dataframe aunderc_C2
for (i in seq_along(list_C2)){
  peak_C2[nrow(peak_C2)+1,] = c(i, max(list_C2[[i]]$amount))
}

# replace the index with the animal names from data frame C3, add column with log transformed data
peak_C2 <- peak_C2 %>%
  add_column(C2,
             .after = "index")%>%
  add_column(log_peak = log10(peak_C2$peak))

peak_C2 <- peak_C2[-c(1)]

## C3 - DENV-2/ZIKV-exposed
# create empty data frame for AUC values
peak_C3 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "peak")
colnames(peak_C3) <- x

#calculate AUC for each animal and output to empty dataframe aunderc_C3
for (i in seq_along(list_C3)){
  peak_C3[nrow(peak_C3)+1,] = c(i, max(list_C3[[i]]$amount))
}

# replace the index with the animal names from data frame C3, add column with log transformed data
peak_C3 <- peak_C3 %>%
  add_column(C3,
             .after = "index")%>%
  add_column(log_peak = log10(peak_C3$peak))

peak_C3 <- peak_C3[-c(1)]

#### Time to peak ####

## C1 - Flavi-naive
# create empty data frame for AUC values
ttpeak_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "ttpeak")
colnames(ttpeak_C1) <- x

#calculate AUC for each animal and output to empty dataframe aunderc_C1
for (i in seq_along(list_C1)){
  ttpeak_C1[nrow(ttpeak_C1)+1,] = c(i, list_C1[[i]]$dpi[which.max(list_C1[[i]]$amount)])
}

# replace the index with the animal names from data frame C3, add column with log transformed data
ttpeak_C1 <- ttpeak_C1 %>%
  add_column(C1,
             .after = "index")

ttpeak_C1 <- ttpeak_C1[-c(1)]

## C2 - ZIKV-exposed
# create empty data frame for AUC values
ttpeak_C2 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "ttpeak")
colnames(ttpeak_C2) <- x

#calculate AUC for each animal and output to empty dataframe aunderc_C2
for (i in seq_along(list_C2)){
  ttpeak_C2[nrow(ttpeak_C2)+1,] = c(i, list_C2[[i]]$dpi[which.max(list_C2[[i]]$amount)])
}

# replace the index with the animal names from data frame C3, add column with log transformed data
ttpeak_C2 <- ttpeak_C2 %>%
  add_column(C2,
             .after = "index")

ttpeak_C2 <- ttpeak_C2[-c(1)]

## C3 - DENV-2/ZIKV-exposed
# create empty data frame for AUC values
ttpeak_C3 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "ttpeak")
colnames(ttpeak_C3) <- x

#calculate AUC for each animal and output to empty dataframe aunderc_C3
for (i in seq_along(C1_truncate_100)){
  ttpeak_C3[nrow(ttpeak_C3)+1,] = c(i, list_C3[[i]]$dpi[which.max(list_C3[[i]]$amount)])
}

# replace the index with the animal names from data frame C3, add column with log transformed data
ttpeak_C3 <- ttpeak_C3 %>%
  add_column(C3,
             .after = "index")

ttpeak_C3 <- ttpeak_C3[-c(1)]

#### Duration ####

## C1 - Flavi-naive
# create empty data frame for duration values
duration_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "duration")
colnames(duration_C1) <- x

#calculate duration for each animal and output to empty dataframe duration_C1
for (i in seq_along(C1_truncate_100)){
  duration_C1[nrow(duration_C1)+1,] = c(i, length(C1_truncate_100[[i]]$dpi[which(C1_truncate_100[[i]]$amount > 100)]))
}

# replace the index with the animal names from data frame C1, add column with log transformed data
duration_C1 <- duration_C1 %>%
  add_column(C1,
             .after = "index")

duration_C1 <- duration_C1[-c(1)]

## C2 - ZIKV-exposed
# create empty data frame for duration values
duration_C2 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "duration")
colnames(duration_C2) <- x

#calculate duration for each animal and output to empty dataframe duration_C2
for (i in seq_along(C2_truncate_100)){
  duration_C2[nrow(duration_C2)+1,] = c(i, length(C2_truncate_100[[i]]$dpi[which(C2_truncate_100[[i]]$amount > 100)]))
}

# replace the index with the animal names from data frame C2, add column with log transformed data
duration_C2 <- duration_C2 %>%
  add_column(C2,
             .after = "index")

duration_C2 <- duration_C2[-c(1)]

## C3 - DENV-2/ZIKV-exposed
# create empty data frame for duration values
duration_C3 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("index", "duration")
colnames(duration_C3) <- x

#calculate duration for each animal and output to empty dataframe duration_C3
for (i in seq_along(C3_truncate_100)){
  duration_C3[nrow(duration_C3)+1,] = c(i, length(C3_truncate_100[[i]]$dpi[which(C3_truncate_100[[i]]$amount > 100)]))
}

# replace the index with the animal names from data frame C3, add column with log transformed data
duration_C3 <- duration_C3 %>%
  add_column(C3,
             .after = "index")

duration_C3 <- duration_C3[-c(1)]

#### ALL ####
## C1 - Flavi-naive - stats
C1_stats_df <- aunderc_C1 %>%
  add_column(peak_C1[2:3], ttpeak_C1[2], duration_C1[2])%>%
  add_column(group = "C1",
             .after = "animal")

C1_stats_summary <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("statistic", "C1")
colnames(C1_stats_summary) <- x
  
C1_stats_summary[nrow(C1_stats_summary)+1,] = c("mean_AUC", mean(C1_stats_df$log_AUC))
C1_stats_summary[nrow(C1_stats_summary)+1,] = c("mean_peak", mean(C1_stats_df$log_peak))
C1_stats_summary[nrow(C1_stats_summary)+1,] = c("mean_ttpeak", mean(C1_stats_df$ttpeak))
C1_stats_summary[nrow(C1_stats_summary)+1,] = c("mean_duration", mean(C1_stats_df$duration))

## C2 - Zika-exposed - stats
C2_stats_df <- aunderc_C2 %>%
  add_column(peak_C2[2:3], ttpeak_C2[2], duration_C2[2])%>%
  add_column(group = "C2",
             .after = "animal")

C2_stats_summary <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("statistic", "C2")
colnames(C2_stats_summary) <- x

C2_stats_summary[nrow(C2_stats_summary)+1,] = c("mean_AUC", mean(C2_stats_df$log_AUC))
C2_stats_summary[nrow(C2_stats_summary)+1,] = c("mean_peak", mean(C2_stats_df$log_peak))
C2_stats_summary[nrow(C2_stats_summary)+1,] = c("mean_ttpeak", mean(C2_stats_df$ttpeak))
C2_stats_summary[nrow(C2_stats_summary)+1,] = c("mean_duration", mean(C2_stats_df$duration))

## C3 - DENV-2/ZIKV-exposed- stats
C3_stats_df <- aunderc_C3 %>%
  add_column(peak_C3[2:3], ttpeak_C3[2], duration_C3[2])%>%
  add_column(group = "C3",
             .after = "animal")

C3_stats_summary <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("statistic", "C3")
colnames(C3_stats_summary) <- x

C3_stats_summary[nrow(C3_stats_summary)+1,] = c("mean_AUC", mean(C3_stats_df$log_AUC))
C3_stats_summary[nrow(C3_stats_summary)+1,] = c("mean_peak", mean(C3_stats_df$log_peak))
C3_stats_summary[nrow(C3_stats_summary)+1,] = c("mean_ttpeak", mean(C3_stats_df$ttpeak))
C3_stats_summary[nrow(C3_stats_summary)+1,] = c("mean_duration", mean(C3_stats_df$duration))


## one big, happy df (summary and all)

# means for each statistic for each group
ALL_stats_summary <- C1_stats_summary %>%
  add_column(C2_stats_summary[2], C3_stats_summary[2])
ALL_stats_summary <- setDT(ALL_stats_summary)
ALL_stats_summary <- dcast(melt(ALL_stats_summary, id.vars = "statistic"), variable ~ statistic)
write.csv(ALL_stats_summary,"~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned/043-denv-3-viremia-stats-summary.csv", row.names = FALSE)

# values for each animal
ALL_stats_df <- rbind(C1_stats_df, C2_stats_df, C3_stats_df)
write.csv(ALL_stats_df,"~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned/043-denv-3-viremia-stats-all.csv", row.names = FALSE)

# create CSV for copy paste into the master spreadsheet
master_stats_df <- select(ALL_stats_df, animal, AUC, peak, ttpeak, duration)
master_stats_df$animal <- as.character(master_stats_df$animal)
master_stats_df <- master_stats_df[order(master_stats_df$animal),]
write.csv(master_stats_df,"~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned/043-denv-3-master-viremia-stats-all.csv", row.names = FALSE)

#### AUC Plot ####
AUC_plot <- ggplot(ALL_stats_df, aes(group, log_AUC, group = group, color = animal)) + 
  geom_point(position = position_jitter(w = 0.05, h = 0)) + 
  #geom_errorbar(data = ALL_stats_summary, aes(ALL_stats_summary$variable, ymax = mean_AUC, ymin = mean_AUC),
   #            size = 0.5, inherit.aes = F, width = 0.2) + 
  geom_signif(comparisons =  list(c("Naive", "ZIKV-history"), 
                                  c("Naive", "DENV-2/ZIKV-history"), 
                                  c("ZIKV-history", "DENV-2/ZIKV-history")), 
              test='t.test', map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "NS"=2),
              y_position = c(7.87, 8.67, 9.47), tip_length = 0) + 
  coord_cartesian(ylim = c(2, 10)) +
  scale_color_manual(
    values = anihex,
    breaks = c(color_break_C1, color_break_C2, color_break_C3)) + #Orange
  labs(
    y = expression(paste("Area Under the Curve (log"[10],")"))) + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),     #removes ver. grid lines
    panel.grid.minor = element_blank(),     #removes hor. grid lines
    legend.position = "none",               #remove legend
    legend.title = element_blank(),         #remove legend title
    plot.title = element_text(hjust = 0),   #move title left
    plot.margin = margin(10, 10, 10, 20),   #room for two-line x-axis title
    strip.text.x = element_text(size = 10), #change text size
    strip.text.y = element_text(size = 10), #change text size
    axis.title.x = element_blank(),         #remove x-axis title
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

AUC_plot

## T-test
ALL_ttest.df <- data.frame(ALL_df$group, ALL_df$log_AUC)
ALL_ttest.df <- ALL_ttest.df %>%
  rename(group = ALL_df.group) %>%
  rename(log_AUC = ALL_df.log_AUC)

AUC_ttest.dfN <- filter(ALL_ttest.df, group=="Naive")
AUC_ttest.df1 <- filter(ALL_ttest.df, group=="ZIKV-history")
AUC_ttest.df2 <- filter(ALL_ttest.df, group=="DENV-2/ZIKV-history")

# the all-mighty p-value
AUC_ttestN1 <- t.test(AUC_ttest.dfN$log_AUC, AUC_ttest.df1$log_AUC)
AUC_ttestN2 <- t.test(AUC_ttest.dfN$log_AUC, AUC_ttest.df2$log_AUC)
AUC_ttest12 <- t.test(AUC_ttest.df1$log_AUC, AUC_ttest.df2$log_AUC)
cat("Naive/C1 p-value:",AUC_ttestN1$p.value,"\n");cat("Naive/C2 p-value:",AUC_ttestN2$p.value,"\n");cat("   C1/C2 p-value:",AUC_ttest12$p.value,"\n")

#### Peak Plot ####
Peak_plot <- ggplot(ALL_df, aes(group, log_Peak, group = group, color = organism)) + 
  geom_point(position = position_jitter(w = 0.05, h = 0)) + 
  geom_errorbar(data = ALL_df, aes(group, ymax = ymean_Peak, ymin = ymean_Peak),
                size = 0.5, inherit.aes = F, width = 0.2) + 
  geom_signif(comparisons =  list(c("Naive", "ZIKV-history"), 
                                  c("Naive", "DENV-2/ZIKV-history"), 
                                  c("ZIKV-history", "DENV-2/ZIKV-history")), 
              test='t.test', map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "NS"=2),
              y_position = c(7.87, 8.67, 9.47), tip_length = 0) + 
  coord_cartesian(ylim = c(2, 10)) +
  scale_color_manual(
    values = c(r14059 = "#339933", rh3008 = "#33CC66", rh3009 = "#006633", rh3012 = "#33FF99", rh3018 = "#339966", rh3019 = "#00CC00", rh3035 = "#00CC66", rh3039 = "#00FF00", #Green
               r03041 = "#3366FF", r05029 = "#003366", rh2005 = "#0033FF", rh2799 = "#6699FF", rh2809 = "#3366CC", rhbc48 = "#3399CC", rhbd23 = "#00CCFF", rhbe47 = "#0033CC", #Blue
               r04118 = "#FF6600", r04159 = "#993300", r05055 = "#FF9933", r07035 = "#FFCC99", r11026 = "#FF9900"), #Orange
    breaks = c("r14059", "rh3008", "rh3009", "rh3012", "rh3018", "rh3019", "rh3035", "rh3039", #Green
               "r03041", "r05029", "rh2005", "rh2799", "rh2809", "rhbc48", "rhbd23", "rhbe47",  #Blue
               "r04118", "r04159", "r05055", "r07035", "r11026")) + #Orange
  labs(
    y = expression(paste("Peak viremia (log"[10],")"))) + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),     #removes ver. grid lines
    panel.grid.minor = element_blank(),     #removes hor. grid lines
    legend.position = "none",               #remove legend
    legend.title = element_blank(),         #remove legend title
    plot.title = element_text(hjust = 0),   #move title left
    plot.margin = margin(10, 10, 10, 20),   #room for two-line x-axis title
    strip.text.x = element_text(size = 10), #change text size
    strip.text.y = element_text(size = 10), #change text size
    axis.title.x = element_blank(),         #remove x-axis title
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## T-test
ALL_ttest.df_Peak <- data.frame(ALL_df$group, ALL_df$log_Peak)
ALL_ttest.df_Peak <- ALL_ttest.df_Peak %>%
  rename(group = ALL_df.group) %>%
  rename(log_Peak = ALL_df.log_Peak)

ALL_ttest.df_Peak_N <- filter(ALL_ttest.df_Peak, group=="Naive")
ALL_ttest.df_Peak_1 <- filter(ALL_ttest.df_Peak, group=="ZIKV-history")
ALL_ttest.df_Peak_2 <- filter(ALL_ttest.df_Peak, group=="DENV-2/ZIKV-history")

# the all-mighty p-value
ALL_ttest.df_Peak_N1 <- t.test(ALL_ttest.df_Peak_N$log_Peak, ALL_ttest.df_Peak_1$log_Peak)
ALL_ttest.df_Peak_N2 <- t.test(ALL_ttest.df_Peak_N$log_Peak, ALL_ttest.df_Peak_2$log_Peak)
ALL_ttest.df_Peak_12 <- t.test(ALL_ttest.df_Peak_1$log_Peak, ALL_ttest.df_Peak_2$log_Peak)

#### ttPeak Plot ####
ttPeak_plot <- ggplot(ALL_df, aes(group, ttPeak, group = group, color = organism)) + 
  geom_point(position = position_jitter(w = 0.05, h = 0)) + 
  geom_errorbar(data = ALL_df, aes(group, ymax = ymean_ttPeak, ymin = ymean_ttPeak),
                size = 0.5, inherit.aes = F, width = 0.2) + 
  geom_signif(comparisons =  list(c("Naive", "ZIKV-history"), 
                                  c("Naive", "DENV-2/ZIKV-history"), 
                                  c("ZIKV-history", "DENV-2/ZIKV-history")), 
              test='t.test', map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "NS"=2),
              y_position = c(11, 12.5, 14), tip_length = 0) + 
  coord_cartesian(ylim = c(0, 15)) +
  scale_color_manual(
    values = c(r14059 = "#339933", rh3008 = "#33CC66", rh3009 = "#006633", rh3012 = "#33FF99", rh3018 = "#339966", rh3019 = "#00CC00", rh3035 = "#00CC66", rh3039 = "#00FF00", #Green
               r03041 = "#3366FF", r05029 = "#003366", rh2005 = "#0033FF", rh2799 = "#6699FF", rh2809 = "#3366CC", rhbc48 = "#3399CC", rhbd23 = "#00CCFF", rhbe47 = "#0033CC", #Blue
               r04118 = "#FF6600", r04159 = "#993300", r05055 = "#FF9933", r07035 = "#FFCC99", r11026 = "#FF9900"), #Orange
    breaks = c("r14059", "rh3008", "rh3009", "rh3012", "rh3018", "rh3019", "rh3035", "rh3039", #Green
               "r03041", "r05029", "rh2005", "rh2799", "rh2809", "rhbc48", "rhbd23", "rhbe47",  #Blue
               "r04118", "r04159", "r05055", "r07035", "r11026")) + #Orange
  labs(
    y = expression(paste("Time to peak viremia (days)"))) + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),     #removes ver. grid lines
    panel.grid.minor = element_blank(),     #removes hor. grid lines
    legend.position = "none",               #remove legend
    legend.title = element_blank(),         #remove legend title
    plot.title = element_text(hjust = 0),   #move title left
    plot.margin = margin(10, 10, 10, 20),   #room for two-line x-axis title
    strip.text.x = element_text(size = 10), #change text size
    strip.text.y = element_text(size = 10), #change text size
    axis.title.x = element_blank(),         #remove x-axis title
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## T-test
ALL_ttest.df_ttPeak <- data.frame(ALL_df$group, ALL_df$ttPeak)
ALL_ttest.df_ttPeak <- ALL_ttest.df_ttPeak %>%
  rename(group = ALL_df.group) %>%
  rename(ttPeak = ALL_df.ttPeak)

ttest_ttPeak.dfN <- filter(ALL_ttest.df_ttPeak, group=="Naive")
ttest_ttPeak.df1 <- filter(ALL_ttest.df_ttPeak, group=="ZIKV-history")
ttest_ttPeak.df2 <- filter(ALL_ttest.df_ttPeak, group=="DENV-2/ZIKV-history")

# the all-mighty p-value
ttPeak_ttestN1 <- t.test(ttest_ttPeak.dfN$ttPeak, ttest_ttPeak.df1$ttPeak)
ttPeak_ttestN2 <- t.test(ttest_ttPeak.dfN$ttPeak, ttest_ttPeak.df2$ttPeak)
ttPeak_ttest12 <- t.test(ttest_ttPeak.df1$ttPeak, ttest_ttPeak.df2$ttPeak)
cat("     Naive/ZIKV-history p-value:",ttPeak_ttestN1$p.value,"\n");cat("Naive/DENV-2/ZIKV-history p-value:",ttPeak_ttestN2$p.value,"\n");cat(" ZIKV/DENV-2/ZIKV-history p-value:",ttPeak_ttest12$p.value,"\n")
cat("                  Naive mean value:",ALL_df$ymean_ttPeak[which(ALL_df$organism=="rh3019")],"\n");cat("Naive/DENV-2/ZIKV-history mean value:",ALL_df$ymean_ttPeak[which(ALL_df$organism=="rh2799")],"\n");cat(" ZIKV/DENV-2/ZIKV-history mean value:",ALL_df$ymean_ttPeak[which(ALL_df$organism=="r11026")],"\n")

#### Dur Plot ####
dur_plot <- ggplot(ALL_df, aes(group, dur, group = group, color = organism)) + 
  geom_point(position = position_jitter(w = 0.05, h = 0)) + 
  geom_errorbar(data = ALL_df, aes(group, ymax = ymean_dur, ymin = ymean_dur),
                size = 0.5, inherit.aes = F, width = 0.2) + 
  geom_signif(comparisons =  list(c("Naive", "ZIKV-history"), 
                                  c("Naive", "DENV-2/ZIKV-history")), 
              test='t.test', map_signif_level = F,
              y_position = c(11, 12.5), tip_length = 0) + 
  geom_signif(comparisons =  list(c("ZIKV-history", "DENV-2/ZIKV-history")), 
              test='t.test', map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "NS"=2),
              y_position = c(14), tip_length = 0) + 
  coord_cartesian(ylim = c(0, 15)) +
  scale_color_manual(
    values = c(r14059 = "#339933", rh3008 = "#33CC66", rh3009 = "#006633", rh3012 = "#33FF99", rh3018 = "#339966", rh3019 = "#00CC00", rh3035 = "#00CC66", rh3039 = "#00FF00", #Green
               r03041 = "#3366FF", r05029 = "#003366", rh2005 = "#0033FF", rh2799 = "#6699FF", rh2809 = "#3366CC", rhbc48 = "#3399CC", rhbd23 = "#00CCFF", rhbe47 = "#0033CC", #Blue
               r04118 = "#FF6600", r04159 = "#993300", r05055 = "#FF9933", r07035 = "#FFCC99", r11026 = "#FF9900"), #Orange
    breaks = c("r14059", "rh3008", "rh3009", "rh3012", "rh3018", "rh3019", "rh3035", "rh3039", #Green
               "r03041", "r05029", "rh2005", "rh2799", "rh2809", "rhbc48", "rhbd23", "rhbe47",  #Blue
               "r04118", "r04159", "r05055", "r07035", "r11026")) + #Orange
  labs(
    y = expression(paste("Duration of detectable viremia (days)"))) + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),     #removes ver. grid lines
    panel.grid.minor = element_blank(),     #removes hor. grid lines
    legend.position = "none",               #remove legend
    legend.title = element_blank(),         #remove legend title
    plot.title = element_text(hjust = 0),   #move title left
    plot.margin = margin(10, 10, 10, 20),   #room for two-line x-axis title
    strip.text.x = element_text(size = 10), #change text size
    strip.text.y = element_text(size = 10), #change text size
    axis.title.x = element_blank(),         #remove x-axis title
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## T-test
ALL_ttest.df_dur <- data.frame(ALL_df$group, ALL_df$dur)
ALL_ttest.df_dur <- ALL_ttest.df_dur %>%
  rename(group = ALL_df.group) %>%
  rename(dur = ALL_df.dur)

ttest_dur.dfN <- filter(ALL_ttest.df_dur, group=="Naive")
ttest_dur.df1 <- filter(ALL_ttest.df_dur, group=="ZIKV-history")
ttest_dur.df2 <- filter(ALL_ttest.df_dur, group=="DENV-2/ZIKV-history")

# the all-mighty p-value
dur_ttestN1 <- t.test(ttest_dur.dfN$dur, ttest_dur.df1$dur)
dur_ttestN2 <- t.test(ttest_dur.dfN$dur, ttest_dur.df2$dur)
dur_ttest12 <- t.test(ttest_dur.df1$dur, ttest_dur.df2$dur)
cat("     Naive/ZIKV-history p-value:",dur_ttestN1$p.value,"\n");cat("Naive/DENV-2/ZIKV-history p-value:",dur_ttestN2$p.value,"\n");cat(" ZIKV/DENV-2/ZIKV-history p-value:",dur_ttest12$p.value,"\n")
cat("                  Naive mean value:",ALL_df$ymean_dur[which(ALL_df$organism=="rh3019")],"\n");cat("Naive/DENV-2/ZIKV-history mean value:",ALL_df$ymean_dur[which(ALL_df$organism=="rh2799")],"\n");cat(" ZIKV/DENV-2/ZIKV-history mean value:",ALL_df$ymean_dur[which(ALL_df$organism=="r11026")],"\n")

#### All plots ####
Plot_NAIVE # 600 x 300
Plot_Blue
Plot_Orange
grid.arrange(Plot_NAIVE,Plot_Blue,Plot_Orange, nrow = 3) # 900 x 600 vert
grid.arrange(Plot_NAIVE,Plot_Blue,Plot_Orange, nrow = 1) # 1500 x 300 horiz

AUC_plot # 300 x 350
cat("AUC:\n");cat("       Naive vs. ZIKV-history p-value:",AUC_ttestN1$p.value,"\n");cat("Naive vs. DENV-2/ZIKV-history p-value:",AUC_ttestN2$p.value,"\n");cat(" ZIKV vs. DENV-2/ZIKV-history p-value:",AUC_ttest12$p.value,"\n");cat("                     Naive mean value:",ALL_df$ymean_AUC[which(ALL_df$organism=="rh3019")],"\n");cat("              ZIKV-history mean value:",ALL_df$ymean_AUC[which(ALL_df$organism=="rh2799")],"\n");cat("       DENV-2/ZIKV-history mean value:",ALL_df$ymean_AUC[which(ALL_df$organism=="r11026")],"\n")

Peak_plot # 300 x 350
cat("Peak:\n");cat("       Naive vs. ZIKV-history p-value:",ALL_ttest.df_Peak_N1$p.value,"\n");cat("Naive vs. DENV-2/ZIKV-history p-value:",ALL_ttest.df_Peak_N2$p.value,"\n");cat(" ZIKV vs. DENV-2/ZIKV-history p-value:",ALL_ttest.df_Peak_12$p.value,"\n");cat("                     Naive mean value:",ALL_df$ymean_Peak[which(ALL_df$organism=="rh3019")],"\n");cat("              ZIKV-history mean value:",ALL_df$ymean_Peak[which(ALL_df$organism=="rh2799")],"\n");cat("       DENV-2/ZIKV-history mean value:",ALL_df$ymean_Peak[which(ALL_df$organism=="r11026")],"\n")

ttPeak_plot # 300 x 350
cat("ttPeak:\n");cat("       Naive vs. ZIKV-history p-value:",ttPeak_ttestN1$p.value,"\n");cat("Naive vs. DENV-2/ZIKV-history p-value:",ttPeak_ttestN2$p.value,"\n");cat(" ZIKV vs. DENV-2/ZIKV-history p-value:",ttPeak_ttest12$p.value,"\n");cat("                     Naive mean value:",ALL_df$ymean_ttPeak[which(ALL_df$organism=="rh3019")],"\n");cat("              ZIKV-history mean value:",ALL_df$ymean_ttPeak[which(ALL_df$organism=="rh2799")],"\n");cat("       DENV-2/ZIKV-history mean value:",ALL_df$ymean_ttPeak[which(ALL_df$organism=="r11026")],"\n")

dur_plot # 300 x 350
cat("Dur:\n");cat("       Naive vs. ZIKV-history p-value:",dur_ttestN1$p.value,"\n");cat("Naive vs. DENV-2/ZIKV-history p-value:",dur_ttestN2$p.value,"\n");cat(" ZIKV vs. DENV-2/ZIKV-history p-value:",dur_ttest12$p.value,"\n");cat("                     Naive mean value:",ALL_df$ymean_dur[which(ALL_df$organism=="rh3019")],"\n");cat("              ZIKV-history mean value:",ALL_df$ymean_dur[which(ALL_df$organism=="rh2799")],"\n");cat("       DENV-2/ZIKV-history mean value:",ALL_df$ymean_dur[which(ALL_df$organism=="r11026")],"\n")

# all
grid.arrange(AUC_plot, dur_plot, Peak_plot, ttPeak_plot, ncol = 4, nrow = 1)  # 1200 x 350
grid.arrange(AUC_plot, dur_plot, Peak_plot, ttPeak_plot, ncol = 2, nrow = 2)  # 700 x 700
cat("AUC:\n");cat("       Naive vs. ZIKV-history p-value:",AUC_ttestN1$p.value,"\n");cat("Naive vs. DENV-2/ZIKV-history p-value:",AUC_ttestN2$p.value,"\n");cat(" ZIKV vs. DENV-2/ZIKV-history p-value:",AUC_ttest12$p.value,"\n");cat("                     Naive mean value:",ALL_df$ymean_AUC[which(ALL_df$organism=="rh3019")],"\n");cat("              ZIKV-history mean value:",ALL_df$ymean_AUC[which(ALL_df$organism=="rh2799")],"\n");cat("       DENV-2/ZIKV-history mean value:",ALL_df$ymean_AUC[which(ALL_df$organism=="r11026")],"\n");cat("\nPeak:\n");cat("       Naive vs. ZIKV-history p-value:",ALL_ttest.df_Peak_N1$p.value,"\n");cat("Naive vs. DENV-2/ZIKV-history p-value:",ALL_ttest.df_Peak_N2$p.value,"\n");cat(" ZIKV vs. DENV-2/ZIKV-history p-value:",ALL_ttest.df_Peak_12$p.value,"\n");cat("                     Naive mean value:",ALL_df$ymean_Peak[which(ALL_df$organism=="rh3019")],"\n");cat("              ZIKV-history mean value:",ALL_df$ymean_Peak[which(ALL_df$organism=="rh2799")],"\n");cat("       DENV-2/ZIKV-history mean value:",ALL_df$ymean_Peak[which(ALL_df$organism=="r11026")],"\n");cat("\nttPeak:\n");cat("       Naive vs. ZIKV-history p-value:",ttPeak_ttestN1$p.value,"\n");cat("Naive vs. DENV-2/ZIKV-history p-value:",ttPeak_ttestN2$p.value,"\n");cat(" ZIKV vs. DENV-2/ZIKV-history p-value:",ttPeak_ttest12$p.value,"\n");cat("                     Naive mean value:",ALL_df$ymean_ttPeak[which(ALL_df$organism=="rh3019")],"\n");cat("              ZIKV-history mean value:",ALL_df$ymean_ttPeak[which(ALL_df$organism=="rh2799")],"\n");cat("       DENV-2/ZIKV-history mean value:",ALL_df$ymean_ttPeak[which(ALL_df$organism=="r11026")],"\n");cat("\nDur:\n");cat("       Naive vs. ZIKV-history p-value:",dur_ttestN1$p.value,"\n");cat("Naive vs. DENV-2/ZIKV-history p-value:",dur_ttestN2$p.value,"\n");cat(" ZIKV vs. DENV-2/ZIKV-history p-value:",dur_ttest12$p.value,"\n");cat("                     Naive mean value:",ALL_df$ymean_dur[which(ALL_df$organism=="rh3019")],"\n");cat("              ZIKV-history mean value:",ALL_df$ymean_dur[which(ALL_df$organism=="rh2799")],"\n");cat("       DENV-2/ZIKV-history mean value:",ALL_df$ymean_dur[which(ALL_df$organism=="r11026")],"\n")


