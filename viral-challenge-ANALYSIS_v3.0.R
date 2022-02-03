#### Code Overview ####

## Viremia plotting and stats analysis
# Authors: Chelsea Crooks & Hunter Ries & Luis Haddock

## Overview
# - csv file downloaded from ZEST open research portal
# - Clean data and save cleaned dataframes as CSV files
# - Plot the dataframes as line plots 
# - Calculate viremia statistics (AUC, peak, duration, time to peak) and save data 
# - Plot viremia statistics and save plots
# - Save p-values from viremia statistics

## Input: 
# 1. Raw viral load data from ZIKV portal (see "README-viral-load-raw-data-search-parameters.txt"): 
  # - https://openresearch.labkey.com/ > ZEST > Private > Viral loads in EHR > select Animal IDs
  # - `C*-VL-raw.csv`
# 2. Study-specific data 
  # - csv file with the animal, hexcode color, and challenge date
  
## Output: 
# 1. CSV files with all viral loads:       
  # - `C*-[virus]-VL-clean.csv`
# 2. CSV files with the viremia parameters (AUC, peak, time to peak, duration)
  # - `[virus]-viremia-stats-all.csv` [parameter values for each animal]
  # - `[virus]-viremia-stats-sumary.csv` [mean values for each cohort]
  # - `[virus]-master-viremia-stats-all.csv` [specific format that I wannt data in to copy and paste into a shared spreadsheet]
  # - `[virus]-VL-stat-significance.csv`
# 3. Figures showing viral loads over time:   
  # - `C*-[virus]-VL.png`
# 4. Figures with viremia parameters AUC, peak, time to peak, duration)
  # - `[virus]-AUC-stats.png`
  # - `[virus]-peak-stats.png`
  # - `[virus]-ttpeak-stats.png`
  # - `[virus]-duration-stats.png`
  
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

#### Set parameters & Import data ####

## GOAL: THIS IS THE ONLY SECTION WHERE PARAMETERS NEED TO BE CHANGED IN ORDER FOR CODE TO RUN PROPERLY
## REALITY: If 3 cohorts are imported, but only 1 or 2 have VLs for a given virus, then things in the stats section need to be changed
# Items that need to be changed later in the code are marked with "changethis" so that they can be IDed with command+F

# set WD to folder where raw VL data is stored
setwd("~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_raw/viral-loads")

## Define variables for desired analysis

# list number of cohorts in analysis, cohort names, and cohort datafiles
cohort_number <- (1:3)
number_of_cohorts <- length(cohort_number)
cohort_comparisons <- list(c(1,2), c(1,3), c(2,3)) # include all unique combinations of cohort numbers

# identify cohort names; must be as many list items as there are cohorts in cohort_number
cohortnames <- c("Flavivirus-naive", "ZIKV-exposed", "DENV-2/ZIKV-exposed")

# range of dpi values to include in plots and stats analysis
dpimin <- as.numeric(0)
dpimax <- as.numeric(15)

# if different axes are preferred for plots as compared to stats analysis they can be modified here
# default is to have the same values as the parameters above
plotdpimin <- dpimin # as.numeric() or dpimin
plotdpimax <- dpimax # as.numeric() or dpimax

# sample source for viral loads
samplesource = "Plasma"

# specify virus and viral load assay
virus = "DENV-3" # how the virus will be refered to in plots and filenames; should be no spaces in name
VLassay = "DENV3" # options: "DENV3", "Lanciotti_ZIKV_universal", "DENV2"

# path to folder for saving the clean viral load data
cleanVLcsv = "~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned/"

# path to folder for saving the data plots
plotpng = "~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/figures/VL/"

## import animal data

# CSV should be at least a three column CSV with "animal", "hexcode", "challenge-date-a", and as column headers 
# if animal has received multiple challenges, additional columns can be added using "challenge-date-b", "-c", "-d", etc.
# CSV should have one row for each animal on the study  
# leave cells blank if animal wasn't challenged

animal_study_data <- read_csv("043-VL-study-data.csv")

# import data as data frame
challenge_dates <- as.data.frame(animal_study_data)

# specify which challenge date you want to use from the animal study csv and make sure its a date
challenge_dates$dateA <- mdy(challenge_dates$dateA) # change two instances here
chaldate <- challenge_dates$dateA # change one instance here
names(chaldate) <- challenge_dates$animal

## import viral load data

# import CSV files 
# CSV files must be in the working directory
# CSV file names must have syntax "*-VL-raw.csv"
rawVLcsvs = list.files(pattern="*-VL-raw.csv")
rawVLdata = lapply(rawVLcsvs, read.delim, sep=",")

#### Defining functions and color scheme ####

## function to clean up the raw data from the portal
# takes data frame of raw viral loads as input (a)
VLdataclean <- function(a){
  
  a <- raw_data_C1
  ## remove unnecessary columns
  a$`Key` <- NULL
  a$`Nucleic.Acid` <- NULL
  a$`Viral.Load.Replicates` <- NULL
  a$`Comment` <- NULL
  a$`Experiment.Number` <- NULL
  a$`RNA.Isolation.Method` <- NULL
  
  ## rename remaining columns
  a <- a %>%
    rename(animal = `Participant.ID`) %>%
    rename(date = `Date`) %>%
    rename(amount = `Viral.Load`) %>%
    rename(assay = `Assay`) %>%
    rename(equivocal = `Equivocal`) %>%
    rename(source = `Sample.Source`)
  
  # make sure its a dataframe
  a <- as.data.frame(a)
  
  ## factor by animal and by assay
  a$animal <- as.factor(a$animal)
  a$assay <- as.factor(a$assay)
  
  ## viral titers are numbers
  a$amount <- as.numeric(a$amount)
  
  ## dates are dates
  a$date <- ymd(a$date)
  
  ## make sure all VLs come from the specified sample source and then remove source column
  a <- filter(a, a$source == samplesource)
  a$source <- NULL
  
  ## remove all rows that are 'equivocal'
  
  # select only rows that are NOT equivocal & then remove equivocal column
  a <- filter(a, a$equivocal == "")
  a$equivocal <- NULL
  
  ## select only rows that use the DENV3 assay
  a <- filter(a, a$assay == VLassay)
}

## function to add a column to the raw VL data from the portal and add a column with the DPI of the data point
# input is a list of sequential integers, with one integer for each animal in the group
addDPI <- function(b){
  for (y in 1:b) { # for loop to cycle through each of the indicies in the list
    list_C1[[n]]$dpi[y] <- # create a dpi column
      as.numeric(difftime( # dpi column is numeric and is the difference in time between
        list_C1[[n]]$date[y], # the date in the date column in the specified list
        chaldate[[as.character(unique(list_C1[[n]]$animal))]], # and the challenge date of the animal
        units = "days"))
  }
  return(list_C1[[n]]) 
  }

## modify the data so that it is able to be plotted on a log scale
# takes a dataframe of all of the cleaned viral load data for the cohort
logtransform <- function(d){
  d$amount <- as.character(d$amount)
  d$amount[d$amount == "0"] <- ("0.001") # change all 0 values to 0.001 so that they are plotted on log scale
  d$amount <- as.numeric(d$amount)
  d$amount <- log10(d$amount) # log transform VLs
  d$amount <- as.character(d$amount) 
  d$amount[d$amount == "-3"] <- ("0") # set values that are -3 on log scale -- which were previously 0 -- to 0
  d$amount <- as.numeric(d$amount) 
  #d$amount[d$amount == -5] <- (NA) 
  d$dpi <- as.integer(d$dpi) 
  return(d)
}

## ggplot code for plotting viremia data
# takes log transformed data frame as input
viremiaplot <- function(e){
  Plot_C1 <- ggplot(na.omit(e), aes(dpi, amount, group = animal, color = animal)) + 
    geom_line(size = .5) + 
    scale_y_continuous(
      breaks = c(2, 3, 4, 5, 6, 7, 8),
      labels = c(2, 3, 4, 5, 6, 7, 8),
      expand = c(0, 0)) + 
    #scale_x_continuous(
      #breaks = plotdpimin:plotdpimax, 
      #labels = plotdpimin:plotdpimax,
      #expand = c(0, 0)) + 
    coord_cartesian(
      ylim = c(2, 8),
      xlim = c((plotdpimin - 1),(plotdpimax + 1))) + 
    scale_color_manual(
      values = anihex, 
      breaks = color_break_C1) +
    labs(
      x = "Days post-infection",
      y = paste("log10", virus, "vRNA copies/mL", sep = " "),
      title = paste(cohortnames[[x]], "cohort")) + 
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
}

## ggsave code to save the plot output 
# takes a ggplot list a input
png_VLplot <- function(f){
  ggsave(
    paste(paste("C", cohort_number[[x]], sep=''), virus, "VL.png", sep="-"), # file name
    plot = f,
    device = NULL,
    path = plotpng,
    scale = 1,
    width = 6,
    height = 4,
    units = c("in"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
  )
}

## define color scheme information using data from the input animmal study csv
color_scheme <- as.data.frame(animal_study_data)
colors <- distinct(color_scheme, animal, hexcode)
anihex <- colors$hexcode
names(anihex) <- colors$animal

#### Clean + Plot VL data ####

## for loop to perform these functions for each of the cohorts defined above
for (i in cohort_number){
  x=i
  raw_data_C1 <- as.data.frame(rawVLdata[[x]])

  ## call function to clean data
  raw_data_C1 <- VLdataclean(raw_data_C1)

  ## Create df with cohort animals
  raw_data_C1$animal <- as.character(raw_data_C1$animal)
  raw_data_C1$animal <- as.factor(raw_data_C1$animal)
  C1 <- unique(raw_data_C1[c("animal")])

  ## Create list with each animal in the group
  list_C1 <- split(raw_data_C1,raw_data_C1$animal)

  ## calculate dpi -- there must be a VL for day 0 in order for function to work

  # create a list with an integer for each indexes in list_C1
  length_list_C1 <- 1:length(list_C1)

  # for loop to calculate DPI based on the challenge date and to create the DPI column
  # if/else 1 - make sure that there is some VL data for the cohort, else return "no VL data"
  # if/else 2 - if there is some VL, make sure that there are more than 2 data points, else return "insufficient VL data"
  if (length(C1$animal) > 0){
  for (i in length_list_C1){ 
    n = i
    length_C1_dates <- length(list_C1[[n]]$date)
    if (length_C1_dates > 2){
    list_C1[[n]] <- addDPI(length_C1_dates) 
    list_C1[[n]] <- filter(list_C1[[n]], list_C1[[n]]$dpi >= dpimin) # remove data points outside the defined min (eliminates pre-infection VLs and VLs from previous studies)
    list_C1[[n]] <- filter(list_C1[[n]], list_C1[[n]]$dpi <= dpimax) # remove data points outside the defined min (eliminates VLs from future studies)
    } else {
      list_C1[[n]]$date <- NULL
      list_C1[[n]]$assay <- NULL
      list_C1[[n]]$amount <- NULL
      #insufficient_VL <- "Insufficient VL data for this cohort"
      #write.csv(insufficient_VL, paste(cleanVLcsv, paste(paste("C", cohort_number[[x]], sep=''), virus, "VL-clean.csv", sep = "-"), sep = ''), row.names = FALSE)
      }
    }
  }
    ## one big, happy df & export to clean folder
    df_C1 <- bind_rows(list_C1)
    if (ncol(df_C1) > 1){
    write.csv(df_C1, paste(cleanVLcsv, paste(paste("C", cohort_number[[x]], sep=''), virus, "VL-clean.csv", sep = "-"), sep = ''), row.names = FALSE)
    } else {
      no_VL <- "No or insufficent VL data for this cohort"
      write.csv(no_VL, paste(cleanVLcsv, paste(paste("C", cohort_number[[x]], sep=''), "no", virus, "VL-clean.csv", sep = "-"), sep = ''), row.names = FALSE)
      }
  
  #### Plot viremia data ####
  
 if (ncol(df_C1) > 1){ # check to make sure there is 
    ## change 0s to 0.001 and log transform
    df_C1 <- logtransform(df_C1)
    
    ## ID color scheme
    color_break_C1 <- as.character(C1$animal)
    
    ## Plot
    Plot_C1 <- viremiaplot(df_C1)
    
    # write out figure to png
    png_VLplot(Plot_C1)
    
    } else {print("No viral load data to plot")}
  }
#### Stats Set-up ####

## set working directory to the location of the cleaned data files
setwd(cleanVLcsv)

## import cleaned data files for each cohort to be used

# import CSV files 
# CSV files must be in the working directory
# CSV file names must have syntax "*-VL-raw.csv"
virusfiletype <- paste("*", virus, "VL-clean.csv", sep="-")
cleanVLcsvs = list.files(pattern=virusfiletype)
cleanVLdata = lapply(cleanVLcsvs, read.delim, sep=",")

cleanVLdata <- Filter(function(x) {nrow(x) >= 2}, cleanVLdata) # import only data files that have more than two rows
lencleanVL <- length(cleanVLdata)
stats_cohorts <- (1:lencleanVL)

## write function to eliminate rows that have a VL < 100
lod_cutoff <- function(x){
  x <- x%>%
    filter(amount > 100)
  return (x)
}

## write function to plot stats (note - max of 6 cohorts can be plotted, otherwise need to change the y coordinates of signif bars)
# g - name of stats comparison (log_AUC, log_peak, ttpeak, duration)
# j - name of mean stats comparison (mean_auc, mean_peak, mean_ttpeak, mean_duration)
# k - y-axis label
# l - y position of the significance bars, options defined below
y_position1 <- c(7.87, 8.37, 8.87, 9.37, 9.87, 10.37)
y_position2 <- c(11, 12, 13, 14, 15, 16)
# m - y axis limits, options defined below
ylim1 <- c(2, 10)
ylim2 <- c(0, 15)

statplot <- function(g, j, k, l, m){
  stat_plot <- ggplot(ALL_stats, aes(group, g, group = group, color = animal)) + 
    geom_point(position = position_jitter(w = 0.05, h = 0)) + 
    geom_errorbar(data = ALL_summary_stats, aes(x = group, ymax = j, ymin = j),
                  size = 0.5, inherit.aes = F, width = 0.2) + 
    geom_signif(comparisons =  cohort_comparisons, # changethis if there is only one cohort
                test='t.test', map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "NS"=2),
                y_position = l, tip_length = 0) + 
    coord_cartesian(ylim = m) +
    scale_x_discrete(
      breaks = unique(ALL_stats$group), 
      labels = cohortnames) + # changethis for variable number of cohorts (can use cohortnames)
    scale_color_manual(
      values = anihex,
      breaks = c(animal_study_data$animal)) + 
    labs(
      y = k) + 
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
}
  
## write function to save the statistics plots
# o - name of the statistic as a string
# p - ggplot list
png_statplot <- function(o, p){
  ggsave(
    paste(virus, o, "stats.png", sep="-"),
    plot = p,
    device = NULL,
    path = plotpng,
    scale = 1,
    width = 4,
    height = 4,
    units = c("in"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
  )
}

# create empty data frame for p-values values
stat_pvalues <- data.frame(matrix(ncol = 0, nrow = lencleanVL))
stat_pvalues <- stat_pvalues%>%
  add_column(comparisons = c("C1/C2","C1/C3","C2/C3"))#"C1/C2", "C1/C3", "C2/C3")) # changethis for variable number of cohorts

## T-test
# q - column in the ALL_stats dataframe (ALL_stats$log_auc, ALL_stats$log_peak, ALL_stats$ttpeak, ALL_stats$duration)
# aa - dataframe (stats_pvalue) to add data to
# bb - name of the column in the p-value datafram (auc, peak, ttpeak, duration)

ttest <- function(q, aa, bb){
  ALL_ttest <- data.frame(ALL_stats$group, q)
  
  ttest_C1 <- filter(ALL_ttest, ALL_stats.group=="C1")
  ttest_C2 <- filter(ALL_ttest, ALL_stats.group=="C2")
  ttest_C3 <- filter(ALL_ttest, ALL_stats.group=="C3") # changethis for variable number of cohorts
  
  # the all-mighty p-value
  ttest12 <- t.test(ttest_C1[2], ttest_C2[2])
  ttest13 <- t.test(ttest_C1[2], ttest_C3[2]) # changethis for variable number of cohorts
  ttest23 <- t.test(ttest_C2[2], ttest_C3[2]) # changethis for variable number of cohorts
  
  # add p-value to dataframe
  #assign("bb", auc)
  
  temp <- c(ttest12$p.value, ttest13$p.value, ttest23$p.value) # changethis for variable number of cohorts
  
  aa <- aa%>%
    add_column(temp)%>%
    rename(!!bb := temp)
  }

#### AUC ####

aunderc <- vector("list",number_of_cohorts)

# df for summary stats
auc_summary <- data.frame(matrix(ncol = 2, nrow = 0))
u <- c("group", "mean_auc")
colnames(auc_summary) <- u

for (i in stats_cohorts){
  x=i
  clean_data_C1 <- as.data.frame(cleanVLdata[[x]])
  list_C1 <- split(clean_data_C1,clean_data_C1$animal)
  C1 <- unique(clean_data_C1[c("animal")])
  
  # list_C1 contains df for each animal in the cohort 
  # C1_truncate_100 contains all VLs > 100 copies/ml (the LOD of the assay)
  C1_truncate_100 <- map(list_C1,
                      ~lod_cutoff(.))

  ## write function to calculate the AUC for each animal

  # create empty data frame for AUC values
  aunderc_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
  newdf <- c("index", "auc")
  colnames(aunderc_C1) <- newdf
  
  #calculate AUC for each animal and output to empty dataframe aunderc_C1
  for (y in seq_along(C1_truncate_100)){
      aunderc_C1[nrow(aunderc_C1)+1,] = c(y, auc(C1_truncate_100[[y]]$dpi, C1_truncate_100[[y]]$amount))
  }
  
  # replace the index with the animal names from data frame C1, add column with log transformed data
  aunderc_C1 <- aunderc_C1 %>%
    add_column(C1,
               .after = "index")%>%
    add_column(log_auc = log10(aunderc_C1$auc))%>%
    add_column(group = paste("C", x, sep=""),
               .after = "animal")
    
  aunderc_C1 <- aunderc_C1[-c(1)]
  
  # save the output to the auc list
  aunderc[[x]] <- aunderc_C1
  
  # add mean of the cohort's AUC to summary data frame
  auc_summary[nrow(auc_summary)+1,] = c(unique(aunderc[[x]]$group), mean(aunderc[[x]]$log_auc))
  }

# convert mean values into numeric
auc_summary$mean_auc <- as.numeric(auc_summary$mean_auc)

#### Peak ####

peak <- vector("list",number_of_cohorts)

# df for summary stats
peak_summary <- data.frame(matrix(ncol = 2, nrow = 0))
t <- c("group", "mean_peak")
colnames(peak_summary) <- t

for (i in stats_cohorts){
  x=i
  clean_data_C1 <- as.data.frame(cleanVLdata[[x]])
  list_C1 <- split(clean_data_C1,clean_data_C1$animal)
  C1 <- unique(clean_data_C1[c("animal")])

  # create empty data frame for peak values
  peak_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
  newdf <- c("index", "peak")
  colnames(peak_C1) <- newdf
  
  #calculate peak for each animal and output to empty dataframe peak_C1
  for (y in seq_along(list_C1)){
    peak_C1[nrow(peak_C1)+1,] = c(y, max(list_C1[[y]]$amount))
  }
  
  # replace the index with the animal names from data frame C1, add column with log transformed data
  peak_C1 <- peak_C1 %>%
    add_column(C1,
               .after = "index")%>%
    add_column(log_peak = log10(peak_C1$peak))%>%
    add_column(group = paste("C", x, sep=""),
               .after = "animal")
  
  peak_C1 <- peak_C1[-c(1)]
  
  # save the output to the peak list
  peak[[x]] <- peak_C1
  
  # add mean of the cohort's peak to summary data frame
  peak_summary[nrow(peak_summary)+1,] = c(unique(peak[[x]]$group), mean(peak[[x]]$log_peak))
}

# convert mean values into numeric
peak_summary$mean_peak <- as.numeric(peak_summary$mean_peak)

#### Time to peak ####

ttpeak <- vector("list",number_of_cohorts)

# df for summary stats
ttpeak_summary <- data.frame(matrix(ncol = 2, nrow = 0))
s <- c("group", "mean_ttpeak")
colnames(ttpeak_summary) <- s

for (i in stats_cohorts){
  x=i
  clean_data_C1 <- as.data.frame(cleanVLdata[[x]])
  list_C1 <- split(clean_data_C1,clean_data_C1$animal)
  C1 <- unique(clean_data_C1[c("animal")])

  # create empty data frame for ttpeak values
  ttpeak_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
  newdf <- c("index", "ttpeak")
  colnames(ttpeak_C1) <- newdf
  
  #calculate AUC for each animal and output to empty dataframe ttpeak_C1
  for (y in seq_along(list_C1)){
    ttpeak_C1[nrow(ttpeak_C1)+1,] = c(y, list_C1[[y]]$dpi[which.max(list_C1[[y]]$amount)])
  }
  
  # replace the index with the animal names from data frame C1, add column with log transformed data
  ttpeak_C1 <- ttpeak_C1 %>%
    add_column(C1,
               .after = "index")%>%
    add_column(group = paste("C", x, sep=""),
               .after = "animal")
  
  ttpeak_C1 <- ttpeak_C1[-c(1)] 
  
  # save the output to the ttpeak list
  ttpeak[[x]] <- ttpeak_C1
  
  # add mean of the cohort's ttpeak to summary data frame
  ttpeak_summary[nrow(ttpeak_summary)+1,] = c(unique(ttpeak[[x]]$group), mean(ttpeak[[x]]$ttpeak))
}

# convert mean values into numeric
ttpeak_summary$mean_ttpeak <- as.numeric(ttpeak_summary$mean_ttpeak)

#### Duration ####

duration <- vector("list",number_of_cohorts)

# df for summary stats
duration_summary <- data.frame(matrix(ncol = 2, nrow = 0))
r <- c("group", "mean_duration")
colnames(duration_summary) <- r

for (i in stats_cohorts){
  x=i
  clean_data_C1 <- as.data.frame(cleanVLdata[[x]])
  list_C1 <- split(clean_data_C1,clean_data_C1$animal)
  C1 <- unique(clean_data_C1[c("animal")])
  
  # list_C1 contains df for each animal in the cohort 
  # C1_truncate_100 contains all VLs > 100 copies/ml (the LOD of the assay)    
  C1_truncate_100 <- map(list_C1,
                         ~lod_cutoff(.))
  
  # create empty data frame for duration values
  duration_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
  newdf <- c("index", "duration")
  colnames(duration_C1) <- newdf
  
  #calculate duration for each animal and output to empty dataframe duration_C1
  for (y in seq_along(C1_truncate_100)){
    duration_C1[nrow(duration_C1)+1,] = c(y, ((max(C1_truncate_100[[y]]$dpi)+1)-(min(C1_truncate_100[[y]]$dpi))))
  }
  
  # replace the index with the animal names from data frame C1, add column with log transformed data
  duration_C1 <- duration_C1 %>%
    add_column(C1,
               .after = "index")%>%
    add_column(group = paste("C", x, sep=""),
               .after = "animal")
  
  duration_C1 <- duration_C1[-c(1)]
  
  # save the output to the duration list
  duration[[x]] <- duration_C1
  
  # add mean of the cohort's duration to summary data frame
  duration_summary[nrow(duration_summary)+1,] = c(unique(duration[[x]]$group), mean(duration[[x]]$duration))
}

# convert mean values into numeric
duration_summary$mean_duration <- as.numeric(duration_summary$mean_duration)

#### Summary Stats ####

# combine all of the summmary stats into one data frame
ALL_summary_stats <- bind_cols(auc_summary, 
                               peak_summary[2], 
                               ttpeak_summary[2], 
                               duration_summary[2])

# combine all of the stats dataframes together
aunderc_all <- bind_rows(aunderc[1:lencleanVL])
peak_all <- bind_rows(peak[1:lencleanVL])
ttpeak_all <- bind_rows(ttpeak[1:lencleanVL])
duration_all <- bind_rows(duration[1:lencleanVL])

ALL_stats <- bind_cols(aunderc_all, peak_all[3:4], ttpeak_all[3], duration_all[3])

## write out data

# means for each statistic for each group
write.csv(ALL_summary_stats, paste(virus, "viremia-stats-summary.csv", sep="-"), row.names = FALSE)

# values for each animal
write.csv(ALL_stats, paste(virus, "viremia-stats-all.csv", sep="-"), row.names = FALSE)

# create CSV for direct copy paste into the master google sheet
master_stats <- ALL_stats[-c(2,4,6)]
write.csv(master_stats, paste(virus, "master-viremia-stats-all.csv", sep="-"), row.names = FALSE)

#### Stat Plots ####

# execute function to plot stats
# g - name of stats comparison (ALL_stats$log_AUC, ALL_stats$log_peak, ALL_stats$ttpeak, ALL_stats$duration)
# j - name of mean stats comparison (ALL_summary_stats$mean_auc, ALL_summary_stats$mean_peak, ALL_summary_stats$mean_ttpeak, ALL_summary_stats$mean_duration)
# k - y-axis label 
  # expression(paste("Area Under the Curve (log"[10],")"))
  # expression(paste("Peak viremia (log"[10],")"))
  # expression(paste("Time to peak viremia (days)"))
  # expression(paste("Duration of detectable viremia (days)"))
# l - y position of the significance bars, options defined above (y_position1, y_position2)
# m - y axis limits, options defined above (ylim1, ylim2)

auc_plot <- statplot(ALL_stats$log_auc, ALL_summary_stats$mean_auc, expression(paste("Area Under the Curve (log"[10],")")), y_position1, ylim1) #changethis to edit plot axes
save_auc_plot <- png_statplot("AUC", auc_plot)

peak_plot <- statplot(ALL_stats$log_peak, ALL_summary_stats$mean_peak, expression(paste("Peak viremia (log"[10],")")), y_position1, ylim1) #changethis to edit plot axes
save_peak_plot <- png_statplot("peak", peak_plot)

ttpeak_plot <- statplot(ALL_stats$ttpeak, ALL_summary_stats$mean_ttpeak, expression(paste("Time to peak viremia (days)")), y_position2, ylim2) #changethis to edit plot axes
save_ttpeak_plot <- png_statplot("ttpeak", ttpeak_plot)

duration_plot <- statplot(ALL_stats$duration, ALL_summary_stats$mean_duration, expression(paste("Duration of detectable viremia (days)")), y_position2, ylim2) #changethis to edit plot axes
save_duration_plot <- png_statplot("duration", duration_plot)

#### Save t-test p-values ####

# append each of the stats to the empty p-value data frame
stat_pvalues <- ttest(ALL_stats$log_auc, stat_pvalues, "auc")
stat_pvalues <- ttest(ALL_stats$log_peak, stat_pvalues, "peak")
stat_pvalues <- ttest(ALL_stats$ttpeak, stat_pvalues, "ttpeak")
stat_pvalues <- ttest(ALL_stats$duration, stat_pvalues, "duration")

# write out csv with the p-value data

write.csv(stat_pvalues, paste(cleanVLcsv, paste(virus,"VL-stat-significance.csv", sep = "-"), sep = ""), row.names = FALSE)
