# nhanes-youth-mims-percentiles
R code to process NHANES 2011–2014 wrist accelerometer data and calculate age- and sex-specific MIMS percentiles in youth aged 6–17 years

## Overview

This code:
1. Processes NHANES 2011–2014 minute-level accelerometer data.
2. Restricts participants to ages 6–17 years.
3. Applies prespecified accelerometer validity criteria.
4. Calculates average daily MIMS.
5. Assigns sex- and age-specific MIMS percentile ranks using published
   BCT reference parameters.

## Accelerometer validity criteria

- PAXDAYM = 2–8
- PAXMTSM != -0.01
- PAXFLGSM blank/missing
- PAXPREDM = 1 or 2
- Valid day: ≥600 valid wear minutes
- Valid participant: ≥3 valid days

## Required data

- PAXMIN_G.xpt
- PAXMIN_H.xpt
- DEMO_G.xpt
- DEMO_H.xpt

NHANES public-use data are available from the National Center for Health Statistics.

## Percentile calculation

The original GAMLSS reference model is not refit.
Published sex- and age-specific Box-Cox t parameters are used to calculate
continuous MIMS percentile ranks.
