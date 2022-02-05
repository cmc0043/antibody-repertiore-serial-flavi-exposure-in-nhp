# antibody-repertiore-serial-flavi-exposure-in-nhp
code and data regarding our project to assessing serial flavivirus infection in non-human primates

## code: 
- `viral-challenge-ANALYSIS_v3.0` -- code used to clean up viral load data; plot viral loads; calculate peak, time to peak, duration, and AUC of viremia; and run t test on viremia parameters

- `DENV-3-challenge-ANALYSIS_v2.0` -- previous, archived version of the code

## data_raw/viremia
- `C*-VL-raw.csv` -- raw data used to create all viral load plots
- `README-viral-load-raw-data-search-parameters` -- specific information about how data was acquired

## data_derived/viremia
outputs from `viral-challenge-ANALYSIS_v3.0`
- `C*-[virus]-VL-cleaned` cleaned viral load data 
- CSV files with the viremia parameters (AUC, peak, time to peak, duration)
    - `[virus]-viremia-stats-all.csv` [parameter values for each animal]
    - `[virus]-viremia-stats-sumary.csv` [mean values for each cohort]
    - `[virus]-master-viremia-stats-all.csv` [specific format that I want data in to copy and paste into a shared spreadsheet]
    - `[virus]-VL-stat-significance.csv` [p-values for each parameter comparison]

## figures 
Figures showing viral loads over time:   
  - `C*-[virus]-VL.png`

Figures with viremia parameters AUC, peak, time to peak, duration)
  - `[virus]-AUC-stats.png`
  - `[virus]-peak-stats.png`
  - `[virus]-ttpeak-stats.png`
  - `[virus]-duration-stats.png`