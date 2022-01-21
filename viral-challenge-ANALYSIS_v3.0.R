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

#### Set parameters & Import data ####

# set WD to folder where raw VL data is stored
setwd("~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_raw/viral-loads")

## Define variables for desired analysis

# list number of cohorts in analysis, cohort names, and cohort datafiles
cohort_number <- (1:3)
number_of_cohorts <- length(cohort_number)

# identify cohort names; must be as many list items as there are cohorts in cohort_number
cohortnames <- c("Flavivirus-naive", "ZIKV-exposed", "DENV-2/ZIKV-exposed")

# range of dpi values to include in plots and stats analysis
dpimin <- as.numeric(0)
dpimax <- as.numeric(15)

# sample source
samplesource = "Plasma"

# specify virus and viral load assay
virus = "DENV-3"
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

# specify which challenge date you want to use and make sure its a date
challenge_dates$dateA <- mdy(challenge_dates$dateA) # change two instances here
chaldate <- challenge_dates$dateA # change one instance here
names(chaldate) <- challenge_dates$animal

## import viral load data

# import CSV files 
# CSV files must be in the working directory
# CSV file names must have syntax "*-VL-raw.csv"
rawVLcsvs = list.files(pattern="*-VL-raw.csv")
rawVLdata = lapply(rawVLcsvs, read.delim, sep=",")

# if you want to import all of the VLs from all cohorts into a single file

#tbl <-
  #list.files(pattern = "*-VL-raw.csv") %>% 
  #map_df(~read_csv(.))


#### Defining functions and color scheme ####

VLdataclean <- function(a){
  a$`Key` <- NULL
  a$`Nucleic.Acid` <- NULL
  a$`Viral.Load.Replicates` <- NULL
  a$`Comment` <- NULL
  a$`Experiment.Number` <- NULL
  a$`RNA.Isolation.Method` <- NULL
  
  ## rename columns
  a <- a %>%
    rename(animal = `Participant.ID`) %>%
    rename(date = `Date`) %>%
    rename(amount = `Viral.Load`) %>%
    rename(assay = `Assay`) %>%
    rename(equivocal = `Equivocal`) %>%
    rename(source = `Sample.Source`)
  
  a <- as.data.frame(a)
  
  ## factor by animal and by assay
  a$animal <- as.factor(a$animal)
  a$assay <- as.factor(a$assay)
  
  ## viral titers are numbers, ya know
  a$amount <- as.numeric(a$amount)
  
  ## dates are dates
  a$date <- ymd(a$date)
  
  ## make sure all VLs come from plasma and then remove source column
  a <- filter(a, a$source == samplesource)
  a$source <- NULL
  
  ## remove all rows that are 'equivocal'
  
  # replace empty "NA" values with 0 
  #is.na_replace_0_C1 <- a$equivocal                               
  #is.na_replace_0_C1[is.na(is.na_replace_0_C1)] <- 0 
  #a$equivocal <- is.na_replace_0_C1
  
  # convert to factor 
  #a$equivocal <- as.factor(a$equivocal)
  
  # select only rows that are NOT equivocal & then remove equivocal column
  #a <- filter(a, a$equivocal == "0")
  a$equivocal <- NULL
  
  ## select only rows that use the DENV3 assay
  a <- filter(a, a$assay == VLassay)
}

addDPI <- function(b){
  for (y in 1:b) {
    list_C1[[n]]$dpi[y] <- 
      as.numeric(difftime(
        list_C1[[n]]$date[y],
        chaldate[[as.character(unique(list_C1[[n]]$animal))]],
        units = "days"))
  }
  return(list_C1[[n]])
  }

DPIcutoffs <- function(c){
  for (i in c){
    n = i
    #list_C1[[n]]$dpi <- as.numeric(list_C1[[n]]$dpi)
    list_C1[[n]] <- filter(list_C1[[n]], list_C1[[n]]$dpi >= dpimin)
    list_C1[[n]] <- filter(list_C1[[n]], list_C1[[n]]$dpi <= dpimax)
  }
  return(list_C1[[n]])
}

logtransform <- function(dd){
  dd$amount <- as.character(dd$amount)
  dd$amount[dd$amount == "0"] <- ("0.001")
  dd$amount <- as.numeric(dd$amount)
  dd$amount <- log10(dd$amount)
  dd$amount <- as.character(dd$amount) 
  dd$amount[dd$amount == "-3"] <- ("0")
  dd$amount <- as.numeric(dd$amount) 
  dd$amount[dd$amount == -5] <- (NA)
  dd$dpi <- as.integer(dd$dpi) 
  return(dd)
}

viremiaplot <- function(e){
  Plot_C1 <- ggplot(na.omit(e), aes(dpi, amount, group = animal, color = animal)) + 
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

pngplot <- function(f){
  ggsave(
    paste(paste("C", cohort_number[[x]], sep=''), virus, "VL.png", sep="-"),
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

## define color scheme information
color_scheme <- as.data.frame(animal_study_data)
colors <- distinct(color_scheme, animal, hexcode)
anihex <- colors$hexcode
names(anihex) <- colors$animal

#### Clean VL data ####

## for loop to perform these functions for each of the cohorts defined above
for (i in cohort_number){
  x=i
  raw_data_C1 <- as.data.frame(rawVLdata[[x]])

  ## call function to clean data
  raw_data_C1 <- VLdataclean(raw_data_C1)

  ## Create df with cohort animals
  C1 <- unique(raw_data_C1[c("animal")])

  ## Create list with each animal per group
  list_C1 <- split(raw_data_C1,raw_data_C1$animal)

  ## calculate dpi -- there must be a VL for day 0 in order for function to work

  # create a list with an integer for each indexes in list_C1
  length_list_C1 <- 1:length(list_C1)

  # for loop to calculate DPI based on the challenge date and to create the DPI column
  # if/else 1 - 
  # if/else 2 - 
  if (length(C1$animal) > 0){
  for (i in length_list_C1){
    n = i
    length_C1_dates <- length(list_C1[[n]]$date)
    if (length_C1_dates > 2){
    list_C1[[n]] <- addDPI(length_C1_dates) 
    } else {
      insufficient_VL <- "Insufficient VL data for this cohort"
      write.csv(insufficient_VL, paste(cleanVLcsv, paste(paste("C", cohort_number[[x]], sep=''), virus, "VL-clean.csv", sep = "-"), sep = ''), row.names = FALSE)
      }
    }
    # for loop to remove any DPIs outside the specified min max range
    list_C1[[n]] <- DPIcutoffs(length_list_C1)
  
    ## one big, happy df & export to clean folder
    df_C1 <- bind_rows(list_C1)
    write.csv(df_C1, paste(cleanVLcsv, paste(paste("C", cohort_number[[x]], sep=''), virus, "VL-clean.csv", sep = "-"), sep = ''), row.names = FALSE)
    } else {
      no_VL <- "No VL data for this cohort"
      write.csv(no_VL, paste(cleanVLcsv, paste(paste("C", cohort_number[[x]], sep=''), virus, "VL-clean.csv", sep = "-"), sep = ''), row.names = FALSE)
      }
  
  #### Plot viremia data ####
  
 if (length(C1$animal) > 0){
    ## change 0s to 0.001 and log transform
    df_C1 <- logtransform(df_C1)
    
    ## ID color scheme
    color_break_C1 <- as.character(C1$animal)
    
    ## Plot
    Plot_C1 <- viremiaplot(df_C1)
    
    # visualize plot in R studio
    Plot_C1
    
    # write out figure to png
    pngplot(Plot_C1)
    
    } else {print("No viral load data to plot")}
  }
#### Stats Set-up ####

## set working directory to the location of the cleaned data files
setwd("~/research/chelsea_crooks/projects/zikv_denv_interactions/antibody-repertiore-serial-flavi-exposure-in-nhp/data_analysis/VL-cleaned")

## import cleaned data files for each cohort to be used

# import CSV files 
# CSV files must be in the working directory
# CSV file names must have syntax "*-VL-raw.csv"
virusfiletype <- paste("*", virus, "VL-clean.csv", sep="-")
cleanVLcsvs = list.files(pattern=virusfiletype)
cleanVLdata = lapply(cleanVLcsvs, read.delim, sep=",")

# write function to eliminate rows that have a VL < 100
lod_cutoff <- function(x){
  x <- x%>%
    filter(amount > 100)
  return (x)
}
  
#### AUC ####

aunderc <- vector("list",number_of_cohorts)
for (i in cohort_number){
  x=i
  clean_data_C1 <- as.data.frame(cleanVLdata[[x]])
  list_C1 <- split(clean_data_C1,clean_data_C1$animal)
  C1 <- unique(clean_data_C1[c("animal")])
  
  # list_C1 contains df for each animal in the naive cohort    
  C1_truncate_100 <- map(list_C1,
                      ~lod_cutoff(.))

  ## write function to calculate the AUC for each animal

  # create empty data frame for AUC values
  aunderc_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
  newdf <- c("index", "AUC")
  colnames(aunderc_C1) <- newdf
  
  #calculate AUC for each animal and output to empty dataframe aunderc_C1
  for (y in seq_along(C1_truncate_100)){
      aunderc_C1[nrow(aunderc_C1)+1,] = c(y, auc(C1_truncate_100[[y]]$dpi, C1_truncate_100[[y]]$amount))
  }
  
  # replace the index with the animal names from data frame C1, add column with log transformed data
  aunderc_C1 <- aunderc_C1 %>%
    add_column(C1,
               .after = "index")%>%
    add_column(log_AUC = log10(aunderc_C1$AUC))
    
  aunderc_C1 <- aunderc_C1[-c(1)]
  
  # save the output to the auc list
  aunderc[[x]] <- aunderc_C1
  }

#### Peak ####

peak <- vector("list",number_of_cohorts)
for (i in cohort_number){
  x=i
  clean_data_C1 <- as.data.frame(cleanVLdata[[x]])
  list_C1 <- split(clean_data_C1,clean_data_C1$animal)
  C1 <- unique(clean_data_C1[c("animal")])
  # list_C1 contains df for each animal in the naive cohort    
  
  # create empty data frame for AUC values
  peak_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
  newdf <- c("index", "peak")
  colnames(peak_C1) <- newdf
  
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
  
  # save the output to the auc list
  peak[[x]] <- peak_C1
}

#### Time to peak ####

ttpeak <- vector("list",number_of_cohorts)
for (i in cohort_number){
  x=i
  clean_data_C1 <- as.data.frame(cleanVLdata[[x]])
  list_C1 <- split(clean_data_C1,clean_data_C1$animal)
  C1 <- unique(clean_data_C1[c("animal")])

  # create empty data frame for AUC values
  ttpeak_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
  newdf <- c("index", "ttpeak")
  colnames(ttpeak_C1) <- newdf
  
  #calculate AUC for each animal and output to empty dataframe aunderc_C1
  for (i in seq_along(list_C1)){
    ttpeak_C1[nrow(ttpeak_C1)+1,] = c(i, list_C1[[i]]$dpi[which.max(list_C1[[i]]$amount)])
  }
  
  # replace the index with the animal names from data frame C3, add column with log transformed data
  ttpeak_C1 <- ttpeak_C1 %>%
    add_column(C1,
               .after = "index")
  
  ttpeak_C1 <- ttpeak_C1[-c(1)] 
  
  # save the output to the auc list
  ttpeak[[x]] <- ttpeak_C1
}

#### Duration ####

duration <- vector("list",number_of_cohorts)
for (i in cohort_number){
  x=i
  clean_data_C1 <- as.data.frame(cleanVLdata[[x]])
  list_C1 <- split(clean_data_C1,clean_data_C1$animal)
  C1 <- unique(clean_data_C1[c("animal")])
  
  # list_C1 contains df for each animal in the naive cohort    
  C1_truncate_100 <- map(list_C1,
                         ~lod_cutoff(.))
  
  # create empty data frame for duration values
  duration_C1 <- data.frame(matrix(ncol = 2, nrow = 0))
  newdf <- c("index", "duration")
  colnames(duration_C1) <- newdf
  
  #calculate duration for each animal and output to empty dataframe duration_C1
  for (w in seq_along(C1_truncate_100)){
    duration_C1[nrow(duration_C1)+1,] = c(w, length(C1_truncate_100[[w]]$dpi[which(C1_truncate_100[[w]]$amount > 100)]))
  }
  
  # replace the index with the animal names from data frame C1, add column with log transformed data
  duration_C1 <- duration_C1 %>%
    add_column(C1,
               .after = "index")
  
  duration_C1 <- duration_C1[-c(1)]
  
  # save the output to the auc list
  duration[[x]] <- duration_C1
}

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


