# ======================================================================
# NHANES 2011-2014 MIMS PERCENTILES
# Youth aged 6-17 years
#
# IMPORTANT:
# This script DOES NOT refit a GAMLSS model.
#
# Participant MIMS values are converted to age- and sex-specific
# percentiles using the published BCT distribution parameters from:
#
# Belcher BR et al.
# U.S. Population-referenced Percentiles for Wrist-Worn
# Accelerometer-derived Activity.
#
# ======================================================================


# ======================================================================
# 0. PACKAGES
# ======================================================================

library(haven)
library(dplyr)
library(gamlss.dist)
library(openxlsx)

options(scipen = 999)


# ======================================================================
# 1. FILE PATHS
# ======================================================================

data_dir <- "C:/NHANES"

pax_g_file  <- file.path(data_dir, "PAXMIN_G.xpt")
pax_h_file  <- file.path(data_dir, "PAXMIN_H.xpt")

demo_g_file <- file.path(data_dir, "DEMO_G.xpt")
demo_h_file <- file.path(data_dir, "DEMO_H.xpt")

output_dir <- file.path(
  data_dir,
  "NHANES_MIMS_percentiles_age6_17"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ======================================================================
# 2. PREPARE DEMOGRAPHIC DATA
# ======================================================================

prepare_demo <- function(file, cycle_name) {
  
  read_xpt(file) %>%
    
    transmute(
      
      SEQN,
      
      cycle = cycle_name,
      
      ageyr = RIDAGEYR,
      
      sex = case_when(
        RIAGENDR == 1 ~ "Male",
        RIAGENDR == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      
      WTMEC2YR
    ) %>%
    
    filter(
      ageyr >= 6,
      ageyr <= 17,
      !is.na(sex)
    )
}


demo_g <- prepare_demo(
  demo_g_file,
  "2011-2012"
)

demo_h <- prepare_demo(
  demo_h_file,
  "2013-2014"
)


# ======================================================================
# 3. PROCESS ONE PAXMIN CYCLE
# ======================================================================

process_pax_cycle <- function(
    pax_file,
    demo_data) {
  
  message(
    "Processing ",
    basename(pax_file)
  )
  
  
  pax <- read_xpt(pax_file) %>%
    
    # Restrict immediately to eligible participants
    semi_join(
      demo_data %>% select(SEQN),
      by = "SEQN"
    ) %>%
    
    select(
      SEQN,
      PAXDAYM,
      PAXMTSM,
      PAXPREDM,
      PAXFLGSM
    )
  
  
  # --------------------------------------------------------------------
  # Minute-level validity criteria
  # --------------------------------------------------------------------
  
  valid_minutes <- pax %>%
    
    filter(
      
      # Keep assessment days 2 through 8
      PAXDAYM %in% 2:8,
      
      # MIMS must be available
      !is.na(PAXMTSM),
      
      # Remove uncomputable MIMS
      PAXMTSM != -0.01,
      
      # Keep wake wear and sleep wear
      PAXPREDM %in% c(1, 2)
    ) %>%
    
    # Remove minutes with any PAX quality flag
    filter(
      is.na(PAXFLGSM) |
        trimws(as.character(PAXFLGSM)) == ""
    )
  
  
  # --------------------------------------------------------------------
  # Aggregate minute-level MIMS to day
  # --------------------------------------------------------------------
  
  day_level <- valid_minutes %>%
    
    group_by(
      SEQN,
      PAXDAYM
    ) %>%
    
    summarise(
      
      valid_wear_minutes = n(),
      
      daily_MIMS =
        sum(
          PAXMTSM,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) %>%
    
    # Valid day criterion
    filter(
      valid_wear_minutes >= 600
    )
  
  
  # --------------------------------------------------------------------
  # Aggregate valid days to participant
  # --------------------------------------------------------------------
  
  person_level <- day_level %>%
    
    group_by(SEQN) %>%
    
    summarise(
      
      n_valid_days = n(),
      
      avg_daily_MIMS =
        mean(
          daily_MIMS,
          na.rm = TRUE
        ),
      
      mean_valid_wear_minutes =
        mean(
          valid_wear_minutes,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) %>%
    
    # Participant must have at least 3 valid days
    filter(
      n_valid_days >= 3
    ) %>%
    
    inner_join(
      demo_data,
      by = "SEQN"
    )
  
  
  list(
    day_level = day_level,
    person_level = person_level
  )
}


# ======================================================================
# 4. PROCESS 2011-2012
# ======================================================================

g <- process_pax_cycle(
  pax_g_file,
  demo_g
)

saveRDS(
  g$day_level,
  file.path(
    output_dir,
    "NHANES_2011_2012_day_level.rds"
  )
)

saveRDS(
  g$person_level,
  file.path(
    output_dir,
    "NHANES_2011_2012_person_level.rds"
  )
)

rm(g)
gc()


# ======================================================================
# 5. PROCESS 2013-2014
# ======================================================================

h <- process_pax_cycle(
  pax_h_file,
  demo_h
)

saveRDS(
  h$day_level,
  file.path(
    output_dir,
    "NHANES_2013_2014_day_level.rds"
  )
)

saveRDS(
  h$person_level,
  file.path(
    output_dir,
    "NHANES_2013_2014_person_level.rds"
  )
)

rm(h)
gc()


# ======================================================================
# 6. COMBINE CYCLES
# ======================================================================

person_g <- readRDS(
  file.path(
    output_dir,
    "NHANES_2011_2012_person_level.rds"
  )
)

person_h <- readRDS(
  file.path(
    output_dir,
    "NHANES_2013_2014_person_level.rds"
  )
)


analytic <- bind_rows(
  person_g,
  person_h
) %>%
  
  filter(
    ageyr >= 6,
    ageyr <= 17,
    n_valid_days >= 3,
    !is.na(avg_daily_MIMS),
    avg_daily_MIMS > 0
  )


# ======================================================================
# 7. PUBLISHED BCT PARAMETERS
#
# Belcher et al.
# Youth ages 3-19.
#
# Columns:
#   mu     = location / median
#   sigma  = scale
#   nu     = Box-Cox power / skewness
#   tau    = degrees of freedom / kurtosis
#
# Only ages 6-17 are entered because that is the current analytic sample.
# ======================================================================


# ----------------------------------------------------------------------
# MALES
# ----------------------------------------------------------------------

male_reference <- data.frame(
  
  sex = "Male",
  
  ageyr = 6:17,
  
  mu = c(
    20613,
    20365,
    19828,
    19110,
    18168,
    16955,
    15686,
    14774,
    14179,
    13737,
    13424,
    13250
  ),
  
  sigma = c(
    0.15,
    0.16,
    0.17,
    0.17,
    0.17,
    0.19,
    0.21,
    0.22,
    0.22,
    0.22,
    0.22,
    0.22
  ),
  
  nu = c(
    0.87,
    0.87,
    0.86,
    0.84,
    0.81,
    0.78,
    0.74,
    0.69,
    0.61,
    0.51,
    0.37,
    0.19
  ),
  
  tau = c(
    2.69,
    2.66,
    2.64,
    2.61,
    2.59,
    2.57,
    2.54,
    2.52,
    2.49,
    2.47,
    2.44,
    2.42
  )
)


# ----------------------------------------------------------------------
# FEMALES
# ----------------------------------------------------------------------

female_reference <- data.frame(
  
  sex = "Female",
  
  ageyr = 6:17,
  
  mu = c(
    20706,
    20658,
    20144,
    19365,
    18425,
    17348,
    16280,
    15461,
    14923,
    14612,
    14408,
    14209
  ),
  
  sigma = c(
    0.13,
    0.13,
    0.14,
    0.14,
    0.15,
    0.16,
    0.17,
    0.18,
    0.19,
    0.21,
    0.21,
    0.22
  ),
  
  nu = c(
    1.11,
    1.06,
    1.02,
    0.97,
    0.93,
    0.88,
    0.84,
    0.79,
    0.75,
    0.70,
    0.66,
    0.61
  ),
  
  tau = c(
    2.26,
    2.34,
    2.42,
    2.50,
    2.58,
    2.66,
    2.74,
    2.82,
    2.90,
    2.98,
    3.06,
    3.14
  )
)


published_reference <- bind_rows(
  male_reference,
  female_reference
)


# ======================================================================
# 8. MERGE PUBLISHED PARAMETERS WITH PARTICIPANTS
# ======================================================================

analytic <- analytic %>%
  
  left_join(
    published_reference,
    by = c(
      "sex",
      "ageyr"
    )
  )


# Check for any unmatched ages/sex
if (any(is.na(analytic$mu))) {
  
  stop(
    "At least one participant did not match a published age/sex parameter."
  )
}


# ======================================================================
# 9. CALCULATE CONTINUOUS SEX- AND AGE-SPECIFIC MIMS PERCENTILE
#
# Belcher et al. selected the Box-Cox t (BCT) distribution.
#
# percentile = 100 x CDF(MIMS | age, sex)
# ======================================================================

analytic <- analytic %>%
  
  mutate(
    
    MIMS_percentile_raw =
      
      100 *
      
      pBCT(
        
        q = avg_daily_MIMS,
        
        mu = mu,
        
        sigma = sigma,
        
        nu = nu,
        
        tau = tau
      )
  )


# ======================================================================
# 10. OPTIONAL: CONSTRAIN TO A 1-100 REPORTING SCALE
#
# This is only if you explicitly want no value <1.
#
# I recommend retaining MIMS_percentile_raw as the primary value.
# ======================================================================

analytic <- analytic %>%
  
  mutate(
    
    MIMS_percentile_1_100 =
      
      pmin(
        100,
        
        pmax(
          1,
          MIMS_percentile_raw
        )
      )
  )


# ======================================================================
# 11. OPTIONAL INTEGER PERCENTILE
#
# Example:
#   47.3 -> 47
#
# Use only if you need an integer category.
# Keep continuous percentile for analysis.
# ======================================================================

analytic <- analytic %>%
  
  mutate(
    
    MIMS_percentile_integer =
      pmin(
        100,
        
        pmax(
          1,
          round(
            MIMS_percentile_raw
          )
        )
      )
  )


# ======================================================================
# 12. QC: COMPARE PUBLISHED MEDIAN WITH CDF = 50
#
# At MIMS = mu, BCT percentile should be approximately 50.
# ======================================================================

reference_check <- published_reference %>%
  
  mutate(
    
    percentile_at_mu =
      100 *
      
      pBCT(
        q = mu,
        mu = mu,
        sigma = sigma,
        nu = nu,
        tau = tau
      )
  )


print(reference_check)


# ======================================================================
# 13. QC: RECREATE SELECTED PUBLISHED CENTILES
#
# These should approximately reproduce Supplemental Tables S2 and S3.
# Small differences are possible because published parameters are rounded.
# ======================================================================

published_reference_check <- published_reference %>%
  
  rowwise() %>%
  
  mutate(
    
    P3 =
      qBCT(
        0.03,
        mu,
        sigma,
        nu,
        tau
      ),
    
    P5 =
      qBCT(
        0.05,
        mu,
        sigma,
        nu,
        tau
      ),
    
    P10 =
      qBCT(
        0.10,
        mu,
        sigma,
        nu,
        tau
      ),
    
    P25 =
      qBCT(
        0.25,
        mu,
        sigma,
        nu,
        tau
      ),
    
    P50 =
      qBCT(
        0.50,
        mu,
        sigma,
        nu,
        tau
      ),
    
    P75 =
      qBCT(
        0.75,
        mu,
        sigma,
        nu,
        tau
      ),
    
    P90 =
      qBCT(
        0.90,
        mu,
        sigma,
        nu,
        tau
      ),
    
    P95 =
      qBCT(
        0.95,
        mu,
        sigma,
        nu,
        tau
      ),
    
    P97 =
      qBCT(
        0.97,
        mu,
        sigma,
        nu,
        tau
      )
  ) %>%
  
  ungroup()


# ======================================================================
# 14. OPTIONAL: CREATE FULL 1st-99th REFERENCE TABLE
#
# This produces MIMS values corresponding to each integer percentile
# for each age and sex.
# ======================================================================

reference_1_99 <- expand.grid(
  
  sex = c(
    "Male",
    "Female"
  ),
  
  ageyr = 6:17,
  
  percentile = 1:99,
  
  stringsAsFactors = FALSE
) %>%
  
  left_join(
    published_reference,
    by = c(
      "sex",
      "ageyr"
    )
  ) %>%
  
  mutate(
    
    MIMS_reference_value =
      
      qBCT(
        
        p = percentile / 100,
        
        mu = mu,
        
        sigma = sigma,
        
        nu = nu,
        
        tau = tau
      )
  )


# ======================================================================
# 15. ANALYTIC SAMPLE SUMMARY
# ======================================================================

sample_summary <- analytic %>%
  
  group_by(
    sex,
    ageyr
  ) %>%
  
  summarise(
    
    N = n(),
    
    mean_valid_days =
      mean(n_valid_days),
    
    mean_MIMS =
      mean(avg_daily_MIMS),
    
    SD_MIMS =
      sd(avg_daily_MIMS),
    
    mean_percentile =
      mean(MIMS_percentile_raw),
    
    SD_percentile =
      sd(MIMS_percentile_raw),
    
    median_percentile =
      median(MIMS_percentile_raw),
    
    .groups = "drop"
  )


print(
  sample_summary,
  n = Inf
)


# ======================================================================
# 16. SAVE PARTICIPANT-LEVEL DATA
# ======================================================================

write.csv(
  
  analytic,
  
  file.path(
    output_dir,
    "NHANES_2011_2014_MIMS_percentiles_age6_17.csv"
  ),
  
  row.names = FALSE
)


# ======================================================================
# 17. SAVE REFERENCE TABLES
# ======================================================================

write.csv(
  
  published_reference,
  
  file.path(
    output_dir,
    "Belcher_published_BCT_parameters_age6_17.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  reference_1_99,
  
  file.path(
    output_dir,
    "Belcher_MIMS_reference_percentiles_1_99_age6_17.csv"
  ),
  
  row.names = FALSE
)


# ======================================================================
# 18. EXCEL OUTPUT
# ======================================================================

wb <- createWorkbook()


addWorksheet(
  wb,
  "Participant percentiles"
)

writeData(
  wb,
  "Participant percentiles",
  analytic
)


addWorksheet(
  wb,
  "Published parameters"
)

writeData(
  wb,
  "Published parameters",
  published_reference
)


addWorksheet(
  wb,
  "Reference percentiles 1-99"
)

writeData(
  wb,
  "Reference percentiles 1-99",
  reference_1_99
)


addWorksheet(
  wb,
  "Published centile check"
)

writeData(
  wb,
  "Published centile check",
  published_reference_check
)


addWorksheet(
  wb,
  "Sample summary"
)

writeData(
  wb,
  "Sample summary",
  sample_summary
)


saveWorkbook(
  
  wb,
  
  file.path(
    output_dir,
    "NHANES_MIMS_percentile_results.xlsx"
  ),
  
  overwrite = TRUE
)


# ======================================================================
# 19. SESSION INFORMATION
# ======================================================================

sink(
  file.path(
    output_dir,
    "R_sessionInfo.txt"
  )
)

sessionInfo()

sink()