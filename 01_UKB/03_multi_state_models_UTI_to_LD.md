# Project Documentation & Pipeline: Liver Disease & UTI in UK Biobank

This document describes the R-based processing pipeline for the UK Biobank (UKB) Liver Disease & UTI Analysis project.

Objective: Use a 5-state multi-state (illness–death–type) model where everyone starts in State A (Healthy) at baseline (no prevalent diabetes, no prevalent liver disease, no viral hepatitis B15–B19).

States:
- A: Healthy (baseline state)
- C: Diabetes (intermediate state)
- B: Liver diseases (K70–K77) — *final/absorbing*
- D: Liver diseases after diabetes (K70–K77 occurring after diabetes) — *final/absorbing*
- E: Death — *final/absorbing*

Allowed transitions:
- 1. A -> B: incident liver disease (direct from healthy)
- 2. A -> C: incident diabetes
- 3. A -> E: death before diabetes/liver disease
- 4. C -> D: incident liver disease after diabetes
- 5. C -> E: death after diabetes (before liver disease)

The hazards of primary interest are the two liver-disease transitions: A->B (1) and C->D (4), while death (E) acts as a competing absorbing outcome.

Key input is the file: `ukb_analysis_liver_uti_final_cohort_479546_participants_27Jan2026.rds` from `02_cox_models_UTI_to_LD.md`. Content: UKB-derived dataset prepared by me with relevant variables for the analysis, it was filtered already to exclude participants with prevalent liver disease and prevalent viral hepatitis cases.

It contains n=479,546 participants with data on UTIs, liver disease, demographics, lifestyle, comorbidities, and biomarkers. 

The dataset still need some processing before the analysis:
- Remove participants with baseline any diabetes
- Remove participants with date of incident liver disease = date incident any diabetes
- If date of death = date of incident liver disease, add a small number (0.01) to the days of death counted after baseline
- If date of death = date of any indicent diabetes, add a small number (0.01) to the days of death counted after baseline

# R Script: Final Data Preparation and Multi-State Models 

```R
# ==============================================================================
# 1. GLOBAL CONFIGURATION & LIBRARIES 
# ==============================================================================
library(tidyverse) # data manipulation and visualization
library(survival) # to run Cox models
packageVersion("survival") # 3.8.6
library(mstate) # multi-state models
packageVersion("mstate") # 0.3.3
options(width = 200, scipen = 999) # improve console readability

# ==============================================================================
# SECTION 2: DATA LOADING 
# ==============================================================================

# Load the file produced in 02_cox_models_UTI_to_LD.md
df <- readRDS("ukb_analysis_liver_uti_final_cohort_479546_participants_27Jan2026.rds")
dim(df) # 479546    119

# ==============================================================================
# SECTION 3: DEFINING FOLLOW-UP TIME VARIABLES FOR DIABETES & COHORT FILTRATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 Define Date for Incident Diabetes (ICD-10 Codes E10-E14)
# ------------------------------------------------------------------------------

# Identify the Incident Event Date (earliest occurrence of any diabetes-related ICD-10 code E10-E14)
df$ANY_DIABETES_INCIDENT_DATE <- do.call(pmin, c(df[paste0("E", 10:14, "_date")], na.rm = TRUE))

# Verification: Does anyone have an incident date AFTER their follow-up ended?
df %>% filter(ANY_DIABETES_INCIDENT == 1 & ANY_DIABETES_INCIDENT_DATE > date_end_of_followup) %>%  nrow()

# One participant had an incident diabetes date recorded after end of follow-up and was excluded.
# remove this participant since incident date is after follow-up end
df <- df %>% filter(!(ANY_DIABETES_INCIDENT == 1 & ANY_DIABETES_INCIDENT_DATE > date_end_of_followup))

# ------------------------------------------------------------------------------
# 3.2 Remove Participants According to Criteria (Prevalent Diabetes, Incident Dates)
# ------------------------------------------------------------------------------

# Remove participants with baseline any diabetes
df %>% filter(ANY_DIABETES_PREVAL == 1) %>% nrow() # 24545 participants with prevalent diabetes at baseline
df <- df %>% filter(ANY_DIABETES_PREVAL == 0) 
dim(df) # 455000    120

# Remove participants with date of incident liver disease = date incident any diabetes
df %>% filter(ANY_LIVER_INCIDENT_DATE == ANY_DIABETES_INCIDENT_DATE & !is.na(ANY_LIVER_INCIDENT_DATE) & !is.na(ANY_DIABETES_INCIDENT_DATE)) %>% nrow() # 436 participants with same date
df <- df %>% filter(!(ANY_LIVER_INCIDENT_DATE == ANY_DIABETES_INCIDENT_DATE & !is.na(ANY_LIVER_INCIDENT_DATE) & !is.na(ANY_DIABETES_INCIDENT_DATE)))
dim(df) # 454564    120

# ------------------------------------------------------------------------------
# 3.3 Create Time-to-Event Variables for Incident Liver Disease, Diabetes, Death
# Note: if no event, time will be date_end_of_followup
# ------------------------------------------------------------------------------

# Time (days) to incident liver disease (K70-K77)
df <- df %>% mutate(LIVER_DAYS = as.numeric(dplyr::if_else(ANY_LIVER_INCIDENT == 1, ANY_LIVER_INCIDENT_DATE - date_attendance, date_end_of_followup - date_attendance)))

# Time (days) to incident any diabetes (E10-E14)
df <- df %>% mutate(DIABETES_DAYS = as.numeric(dplyr::if_else(ANY_DIABETES_INCIDENT == 1, ANY_DIABETES_INCIDENT_DATE - date_attendance, date_end_of_followup - date_attendance)))

# Time (days) to death
df <- df %>% mutate(DEATH_DAYS = as.numeric(dplyr::if_else(death_status == 1, death_date - date_attendance, date_end_of_followup - date_attendance)))

# Any person with liver disease AND diabetes AND death on the same day?
df %>% filter(ANY_LIVER_INCIDENT == 1 & ANY_DIABETES_INCIDENT == 1 & death_status == 1 & LIVER_DAYS == DIABETES_DAYS & DIABETES_DAYS == DEATH_DAYS) %>% nrow() # 0 participants

# ------------------------------------------------------------------------------
# 3.4 Adjust Death Dates if Coinciding with Incident Events for Liver Disease 
# ------------------------------------------------------------------------------

# If date of death = date of incident liver disease, add epsilon of 0.5 (12 hs) to the days of death counted after baseline (in case both events happen)
df %>% filter(ANY_LIVER_INCIDENT == 1 & death_status == 1 & LIVER_DAYS == DEATH_DAYS) %>% nrow() # 207 participants with same date

df <- df %>%
  mutate(DEATH_DAYS = if_else(ANY_LIVER_INCIDENT == 1 & death_status == 1 & LIVER_DAYS == DEATH_DAYS,
                              DEATH_DAYS + 0.5,
                              DEATH_DAYS))

# PROBLEM: we also need to adjust DIABETES_DAYS for cases with no events (ANY_DIABETES_INCIDENT == 0) since it should match DEATH_DAYS now, example:
df %>% filter(row_id == "6294") %>% select(row_id, ANY_DIABETES_INCIDENT, DIABETES_DAYS, death_status, DEATH_DAYS, ANY_LIVER_INCIDENT, LIVER_DAYS)

df <- df %>% mutate(DIABETES_DAYS = if_else(ANY_DIABETES_INCIDENT == 0 & death_status == 1 & DIABETES_DAYS == DEATH_DAYS - 0.5, DEATH_DAYS, DIABETES_DAYS))

df %>% filter(row_id == "6294") %>% select(row_id, ANY_DIABETES_INCIDENT, DIABETES_DAYS, death_status, DEATH_DAYS, ANY_LIVER_INCIDENT, LIVER_DAYS)

# ------------------------------------------------------------------------------
# 3.4 Adjust Death Dates if Coinciding with Incident Events for Diabetes
# ------------------------------------------------------------------------------

# 2. If date of death = date of any incident diabetes, add epsilon of 0.5 (12 hs) to the days of death counted after baseline (in case both events happen)
df %>% filter(ANY_DIABETES_INCIDENT == 1 & death_status == 1 & DIABETES_DAYS == DEATH_DAYS) %>% nrow() # 130 participants with same date

df <- df %>%
  mutate(DEATH_DAYS = if_else(ANY_DIABETES_INCIDENT == 1 & death_status == 1 & DIABETES_DAYS == DEATH_DAYS,
                              DEATH_DAYS + 0.5,
                              DEATH_DAYS))

# PROBLEM: we also need to adjust LIVER_DAYS for cases with no events (ANY_LIVER_INCIDENT == 0) since it should match DEATH_DAYS now, example:
df %>% filter(row_id == "584") %>% select(row_id, ANY_LIVER_INCIDENT, LIVER_DAYS, death_status, DEATH_DAYS, ANY_DIABETES_INCIDENT, DIABETES_DAYS)

df <- df %>% mutate(LIVER_DAYS = if_else(ANY_LIVER_INCIDENT == 0 & death_status == 1 & LIVER_DAYS == DEATH_DAYS - 0.5, DEATH_DAYS, LIVER_DAYS))

df %>% filter(row_id == "584") %>% select(row_id, ANY_LIVER_INCIDENT, LIVER_DAYS, death_status, DEATH_DAYS, ANY_DIABETES_INCIDENT, DIABETES_DAYS)

# ==============================================================================
# SECTION 4: CREATING FINAL DATASET FOR MULTI-STATE MODELING
# ==============================================================================

# Recode variables as numeric
df$UTI_PREVAL <- as.integer(as.character(df$UTI_PREVAL))
df$ANY_LIVER_INCIDENT <- as.integer(as.character(df$ANY_LIVER_INCIDENT))
df$ANY_DIABETES_INCIDENT <- as.integer(as.character(df$ANY_DIABETES_INCIDENT))
df$death_status <- as.integer(as.character(df$death_status))

# Extract only columns of interest for the multi-state model
df <- df %>% 
  select(row_id, 
  ANY_LIVER_INCIDENT, LIVER_DAYS, # LD event and time
  death_status, DEATH_DAYS, # death event and time
  ANY_DIABETES_INCIDENT, DIABETES_DAYS, # diabetes event and time
  UTI_PREVAL, age_recruitment, sex, ethnicity, household_income, BMI, smoking_status, I10_PREVAL, N18_PREVAL)

# Missing-data handling (we need the real N for the multi-state model)
df <- df %>% drop_na()

str(df)
# 'data.frame':   450586 obs. of  16 variables:

# ==============================================================================
# SECTION 5: BASIC MODEL WITH HEALTHY, LIVER DISEASE AND DEATH TRANSITIONS
# ==============================================================================

# List the Covariates to keep
covs <- c("age_recruitment", "sex", "ethnicity", "household_income", "BMI", "smoking_status", "I10_PREVAL", "N18_PREVAL", "UTI_PREVAL")

# Create the Transition Matrix
tmat <- trans.illdeath(names=c("Healthy","Liver","Dead"))

# View Matrix
print(tmat)
#          to
# from      Healthy Liver Dead
#   Healthy      NA     1    2
#   Liver        NA    NA    3
#   Dead         NA    NA   NA

# Getting Long Format Data for Multi-State Modelling
mssimdat <- msprep(
    data = df,
    trans= tmat,
    time = c(NA, "LIVER_DAYS", "DEATH_DAYS"),  # Times to transitions
    status = c(NA, "ANY_LIVER_INCIDENT", "death_status"),  # Status indicating the transition
    keep = covs,
    id = "row_id" # column name containing these subject ids
)

# Check First 4 Rows
head(mssimdat, 4)

# Overview of the Transitions
events(mssimdat)
# $Frequencies
#          to
# from      Healthy  Liver   Dead no event total entering
#   Healthy       0  13922  32201   404463         450586
#   Liver         0      0   3506    10416          13922
#   Dead          0      0      0    35707          35707

# $Proportions
#          to
# from         Healthy      Liver       Dead   no event
#   Healthy 0.00000000 0.03089754 0.07146471 0.89763774
#   Liver   0.00000000 0.00000000 0.25183163 0.74816837
#   Dead    0.00000000 0.00000000 0.00000000 1.00000000

# Add Transition-Specific Covariates
mssimdat <- expand.covs(mssimdat, covs, append = TRUE, longnames = FALSE)

colnames(mssimdat) # extension .i refers to transition number i
#  [1] "row_id"              "from"                "to"                  "trans"               "Tstart"              "Tstop"               "time"                "status"             
#  [9] "age_recruitment"     "sex"                 "ethnicity"           "household_income"    "BMI"                 "smoking_status"      "I10_PREVAL"          "N18_PREVAL"         
# [17] "UTI_PREVAL"          "age_recruitment.1"   "age_recruitment.2"   "age_recruitment.3"   "sex.1"               "sex.2"               "sex.3"               "ethnicity1.1"       
# [25] "ethnicity1.2"        "ethnicity1.3"        "ethnicity2.1"        "ethnicity2.2"        "ethnicity2.3"        "ethnicity3.1"        "ethnicity3.2"        "ethnicity3.3"       
# [33] "ethnicity4.1"        "ethnicity4.2"        "ethnicity4.3"        "household_income1.1" "household_income1.2" "household_income1.3" "household_income2.1" "household_income2.2"
# [41] "household_income2.3" "household_income3.1" "household_income3.2" "household_income3.3" "household_income4.1" "household_income4.2" "household_income4.3" "household_income5.1"
# [49] "household_income5.2" "household_income5.3" "BMI.1"               "BMI.2"               "BMI.3"               "smoking_status.1"    "smoking_status.2"    "smoking_status.3"   
# [57] "I10_PREVAL.1"        "I10_PREVAL.2"        "I10_PREVAL.3"        "N18_PREVAL.1"        "N18_PREVAL.2"        "N18_PREVAL.3"        "UTI_PREVAL.1"        "UTI_PREVAL.2"       
# [65] "UTI_PREVAL.3"       

# Calculate Cox with "stratification" on the different possible transitions 
cox_mstat <- coxph(Surv(Tstart,Tstop,status) ~ 
          UTI_PREVAL.1 + UTI_PREVAL.2 + UTI_PREVAL.3 +
          age_recruitment.1 + age_recruitment.2 + age_recruitment.3 + 
          sex.1 + sex.2 + sex.3 +
          ethnicity1.1 + ethnicity1.2 + ethnicity1.3 + ethnicity2.1 + ethnicity2.2 + ethnicity2.3 + ethnicity3.1 + ethnicity3.2 + ethnicity3.3 + ethnicity4.1 + ethnicity4.2 + ethnicity4.3 +
          household_income1.1 + household_income1.2 + household_income1.3 + household_income2.1 + household_income2.2 + household_income2.3 + household_income3.1 + household_income3.2 + household_income3.3 + household_income4.1 + household_income4.2 + household_income4.3 + household_income5.1 + household_income5.2 + household_income5.3 +
          BMI.1 + BMI.2 + BMI.3 +
          smoking_status.1 + smoking_status.2 + smoking_status.3 +
          I10_PREVAL.1 + I10_PREVAL.2 + I10_PREVAL.3 +
          N18_PREVAL.1 + N18_PREVAL.2 + N18_PREVAL.3 +
          strata(trans) + cluster(row_id), data=mssimdat, method="breslow")

# Check Final Rows for Primary Predictor UTI_PREVAL
summary(cox_mstat)$coefficients[1:3, ] %>% as.data.frame() %>% mutate(`Pr(>|z|)` = round(`Pr(>|z|)`, 5))
#                    coef exp(coef)   se(coef)  robust se          z Pr(>|z|)
# UTI_PREVAL.1 0.51634479  1.675891 0.05832626 0.05855077  8.8187525  0.00000
# UTI_PREVAL.2 0.40696585  1.502253 0.03934614 0.03970344 10.2501414  0.00000
# UTI_PREVAL.3 0.09382291  1.098365 0.11571433 0.12706706  0.7383732  0.46029

# Check All Results and Export to .txt
summary(cox_mstat)

sink("Multi_State_Model_Liver_Disease_Death_Healthy_Results.txt")
print(summary(cox_mstat))
sink()

# Export Summary of Multi-State Cox Model
broom::tidy(cox_mstat) %>% write_csv("Multi_State_Model_Liver_Disease_Death_Healthy_Results.csv")

# ==============================================================================
# SECTION 6: FINAL MODEL WITH HEALTHY, LIVER DISEASE, DIABETES AND DEATH TRANSITIONS
# ==============================================================================

# Create the Transition Matrix
tmat <- transMat(x = list(c(2, 3, 4), c(), c(2, 4), c()), 
                 names = c("Healthy", "Liver", "Diabetes", "Dead"))

# View Matrix
print(tmat)
#           to
# from       Healthy Liver Diabetes Dead
#   Healthy       NA     1        2    3
#   Liver         NA    NA       NA   NA
#   Diabetes      NA     4       NA    5
#   Dead          NA    NA       NA   NA

# Getting Long Format Data for Multi-State Modelling
mssimdat <- msprep(
    data = df,
    trans= tmat,
    time = c(NA, "LIVER_DAYS", "DIABETES_DAYS", "DEATH_DAYS"),  # Times to transitions
    status = c(NA, "ANY_LIVER_INCIDENT", "ANY_DIABETES_INCIDENT", "death_status"),  # Status indicating the transition
    keep = covs,
    id = "row_id" # column name containing these subject ids
)

# Check First 4 Rows
head(mssimdat, 4)

# Overview of the Transitions
events(mssimdat)
# $Frequencies
#           to
# from       Healthy  Liver Diabetes   Dead no event total entering
#   Healthy        0  12652    22454  29007   386473         450586
#   Liver          0      0        0      0    13922          13922
#   Diabetes       0   1270        0   3194    17990          22454
#   Dead           0      0        0      0    32201          32201

# $Proportions
#           to
# from          Healthy      Liver   Diabetes       Dead   no event
#   Healthy  0.00000000 0.02807899 0.04983288 0.06437617 0.85771196
#   Liver    0.00000000 0.00000000 0.00000000 0.00000000 1.00000000
#   Diabetes 0.00000000 0.05656008 0.00000000 0.14224637 0.80119355
#   Dead     0.00000000 0.00000000 0.00000000 0.00000000 1.00000000

# Add Transition-Specific Covariates
mssimdat <- expand.covs(mssimdat, covs, append = TRUE, longnames = FALSE)

colnames(mssimdat) # extension .i refers to transition number i
#  [1] "row_id"              "from"                "to"                  "trans"               "Tstart"              "Tstop"               "time"                "status"             
#  [9] "age_recruitment"     "sex"                 "ethnicity"           "household_income"    "BMI"                 "smoking_status"      "I10_PREVAL"          "N18_PREVAL"         
# [17] "UTI_PREVAL"          "age_recruitment.1"   "age_recruitment.2"   "age_recruitment.3"   "age_recruitment.4"   "age_recruitment.5"   "sex.1"               "sex.2"              
# [25] "sex.3"               "sex.4"               "sex.5"               "ethnicity1.1"        "ethnicity1.2"        "ethnicity1.3"        "ethnicity1.4"        "ethnicity1.5"       
# [33] "ethnicity2.1"        "ethnicity2.2"        "ethnicity2.3"        "ethnicity2.4"        "ethnicity2.5"        "ethnicity3.1"        "ethnicity3.2"        "ethnicity3.3"       
# [41] "ethnicity3.4"        "ethnicity3.5"        "ethnicity4.1"        "ethnicity4.2"        "ethnicity4.3"        "ethnicity4.4"        "ethnicity4.5"        "household_income1.1"
# [49] "household_income1.2" "household_income1.3" "household_income1.4" "household_income1.5" "household_income2.1" "household_income2.2" "household_income2.3" "household_income2.4"
# [57] "household_income2.5" "household_income3.1" "household_income3.2" "household_income3.3" "household_income3.4" "household_income3.5" "household_income4.1" "household_income4.2"
# [65] "household_income4.3" "household_income4.4" "household_income4.5" "household_income5.1" "household_income5.2" "household_income5.3" "household_income5.4" "household_income5.5"
# [73] "BMI.1"               "BMI.2"               "BMI.3"               "BMI.4"               "BMI.5"               "smoking_status.1"    "smoking_status.2"    "smoking_status.3"   
# [81] "smoking_status.4"    "smoking_status.5"    "I10_PREVAL.1"        "I10_PREVAL.2"        "I10_PREVAL.3"        "I10_PREVAL.4"        "I10_PREVAL.5"        "N18_PREVAL.1"       
# [89] "N18_PREVAL.2"        "N18_PREVAL.3"        "N18_PREVAL.4"        "N18_PREVAL.5"        "UTI_PREVAL.1"        "UTI_PREVAL.2"        "UTI_PREVAL.3"        "UTI_PREVAL.4"       
# [97] "UTI_PREVAL.5"       

# Calculate Cox with "stratification" on the different possible transitions 
cox_mstat <- coxph(Surv(Tstart,Tstop,status) ~ 
          UTI_PREVAL.1 + UTI_PREVAL.2 + UTI_PREVAL.3 + UTI_PREVAL.4 + UTI_PREVAL.5 +
          age_recruitment.1 + age_recruitment.2 + age_recruitment.3 + age_recruitment.4 + age_recruitment.5 +
          sex.1 + sex.2 + sex.3 + sex.4 + sex.5 +
          ethnicity1.1 + ethnicity1.2 + ethnicity1.3 + ethnicity1.4 + ethnicity1.5 +
          ethnicity2.1 + ethnicity2.2 + ethnicity2.3 + ethnicity2.4 + ethnicity2.5 + 
          ethnicity3.1 + ethnicity3.2 + ethnicity3.3 + ethnicity3.4 + ethnicity3.5 + 
          ethnicity4.1 + ethnicity4.2 + ethnicity4.3 + ethnicity4.4 + ethnicity4.5 +
          household_income1.1 + household_income1.2 + household_income1.3 + household_income1.4 + household_income1.5 +
          household_income2.1 + household_income2.2 + household_income2.3 + household_income2.4 + household_income2.5 +
          household_income3.1 + household_income3.2 + household_income3.3 + household_income3.4 + household_income3.5 +
          household_income4.1 + household_income4.2 + household_income4.3 + household_income4.4 + household_income4.5 +
          household_income5.1 + household_income5.2 + household_income5.3 + household_income5.4 + household_income5.5 +
          BMI.1 + BMI.2 + BMI.3 + BMI.4 + BMI.5 +
          smoking_status.1 + smoking_status.2 + smoking_status.3 + smoking_status.4 + smoking_status.5 +
          I10_PREVAL.1 + I10_PREVAL.2 + I10_PREVAL.3 + I10_PREVAL.4 + I10_PREVAL.5 +
          N18_PREVAL.1 + N18_PREVAL.2 + N18_PREVAL.3 + N18_PREVAL.4 + N18_PREVAL.5 +
          strata(trans) + cluster(row_id), data=mssimdat, method="breslow")

# Check Final Rows for Primary Predictor UTI_PREVAL
summary(cox_mstat)$coefficients[1:5, ] %>% as.data.frame() %>% mutate(`Pr(>|z|)` = round(`Pr(>|z|)`, 5))
#                   coef exp(coef)   se(coef)  robust se        z Pr(>|z|)
# UTI_PREVAL.1 0.5120426  1.668696 0.06212364 0.06239661 8.206257  0.00000
# UTI_PREVAL.2 0.4049177  1.499179 0.04712138 0.05053693 8.012313  0.00000
# UTI_PREVAL.3 0.3844434  1.468797 0.04252514 0.04280635 8.980990  0.00000
# UTI_PREVAL.4 0.3697273  1.447340 0.16969100 0.17136072 2.157596  0.03096
# UTI_PREVAL.5 0.3461459  1.413609 0.10392444 0.10953919 3.160019  0.00158

# Check All Results and Export to .txt
summary(cox_mstat)

sink("Multi_State_Model_Liver_Disease_Death_Diabetes_Healthy_Results.txt")
print(summary(cox_mstat))
sink()

# Export Summary of Multi-State Cox Model
broom::tidy(cox_mstat) %>% write_csv("Multi_State_Model_Liver_Disease_Death_Diabetes_Healthy_Results.csv")
```

# R Script: Multi-State Models with GP registered participants (n=218599)

```R
# ==============================================================================
# 1. GLOBAL CONFIGURATION & LIBRARIES 
# ==============================================================================
library(tidyverse) # data manipulation and visualization
library(survival) # to run Cox models
packageVersion("survival") # 3.8.6
library(mstate) # multi-state models
packageVersion("mstate") # 0.3.3
options(width = 200, scipen = 999) # improve console readability

# ==============================================================================
# SECTION 2: DATA LOADING 
# ==============================================================================

# Load the file produced in 02_cox_models_UTI_to_LD.md
df <- readRDS("ukb_analysis_liver_uti_final_cohort_479546_participants_27Jan2026.rds")
dim(df) # 479546    120

# ==============================================================================
# SECTION 2.1: REMOVE PARTICIPANTS WITH NO GP REGISTRATION 
# ==============================================================================

df <- df %>% filter(gp_reg == "Yes")
dim(df) # 218599    120

# ==============================================================================
# SECTION 3: DEFINING FOLLOW-UP TIME VARIABLES FOR DIABETES & COHORT FILTRATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 Define Date for Incident Diabetes (ICD-10 Codes E10-E14)
# ------------------------------------------------------------------------------

# Identify the Incident Event Date (earliest occurrence of any diabetes-related ICD-10 code E10-E14)
df$ANY_DIABETES_INCIDENT_DATE <- do.call(pmin, c(df[paste0("E", 10:14, "_date")], na.rm = TRUE))

# Verification: Does anyone have an incident date AFTER their follow-up ended?
df %>% filter(ANY_DIABETES_INCIDENT == 1 & ANY_DIABETES_INCIDENT_DATE > date_end_of_followup) %>%  nrow()

# 1 participant died at 2016-12-31 (date_end_of_followup and death_date cols) but had incident diabetes date recorded as 2017-01-08, as well as I10_date (hypertension) 
# remove this participant since incident date is after follow-up end
df <- df %>% filter(!(ANY_DIABETES_INCIDENT == 1 & ANY_DIABETES_INCIDENT_DATE > date_end_of_followup))

# ------------------------------------------------------------------------------
# 3.2 Remove Participants According to Criteria (Prevalent Diabetes, Incident Dates)
# ------------------------------------------------------------------------------

# Remove participants with baseline any diabetes
df %>% filter(ANY_DIABETES_PREVAL == 1) %>% nrow() # 11028 participants with prevalent diabetes at baseline
df <- df %>% filter(ANY_DIABETES_PREVAL == 0) 
dim(df) # 207570    120

# Remove participants with date of incident liver disease = date incident any diabetes
df %>% filter(ANY_LIVER_INCIDENT_DATE == ANY_DIABETES_INCIDENT_DATE & !is.na(ANY_LIVER_INCIDENT_DATE) & !is.na(ANY_DIABETES_INCIDENT_DATE)) %>% nrow() # 146 participants with same date
df <- df %>% filter(!(ANY_LIVER_INCIDENT_DATE == ANY_DIABETES_INCIDENT_DATE & !is.na(ANY_LIVER_INCIDENT_DATE) & !is.na(ANY_DIABETES_INCIDENT_DATE)))
dim(df) # 207424    120

# ------------------------------------------------------------------------------
# 3.3 Create Time-to-Event Variables for Incident Liver Disease, Diabetes, Death
# Note: if no event, time will be date_end_of_followup
# ------------------------------------------------------------------------------

# Time (days) to incident liver disease (K70-K77)
df <- df %>% mutate(LIVER_DAYS = as.numeric(dplyr::if_else(ANY_LIVER_INCIDENT == 1, ANY_LIVER_INCIDENT_DATE - date_attendance, date_end_of_followup - date_attendance)))

# Time (days) to incident any diabetes (E10-E14)
df <- df %>% mutate(DIABETES_DAYS = as.numeric(dplyr::if_else(ANY_DIABETES_INCIDENT == 1, ANY_DIABETES_INCIDENT_DATE - date_attendance, date_end_of_followup - date_attendance)))

# Time (days) to death
df <- df %>% mutate(DEATH_DAYS = as.numeric(dplyr::if_else(death_status == 1, death_date - date_attendance, date_end_of_followup - date_attendance)))

# Any person with liver disease AND diabetes AND death on the same day?
df %>% filter(ANY_LIVER_INCIDENT == 1 & ANY_DIABETES_INCIDENT == 1 & death_status == 1 & LIVER_DAYS == DIABETES_DAYS & DIABETES_DAYS == DEATH_DAYS) %>% nrow() # 0 participants

# ------------------------------------------------------------------------------
# 3.4 Adjust Death Dates if Coinciding with Incident Events for Liver Disease 
# ------------------------------------------------------------------------------

# If date of death = date of incident liver disease, add epsilon of 0.5 (12 hs) to the days of death counted after baseline (in case both events happen)
df %>% filter(ANY_LIVER_INCIDENT == 1 & death_status == 1 & LIVER_DAYS == DEATH_DAYS) %>% nrow() # 74 participants with same date

df <- df %>%
  mutate(DEATH_DAYS = if_else(ANY_LIVER_INCIDENT == 1 & death_status == 1 & LIVER_DAYS == DEATH_DAYS,
                              DEATH_DAYS + 0.5,
                              DEATH_DAYS))

# PROBLEM: we also need to adjust DIABETES_DAYS for cases with no events (ANY_DIABETES_INCIDENT == 0) since it should match DEATH_DAYS now, example:
df %>% filter(row_id == "6294") %>% select(row_id, ANY_DIABETES_INCIDENT, DIABETES_DAYS, death_status, DEATH_DAYS, ANY_LIVER_INCIDENT, LIVER_DAYS)

df <- df %>% mutate(DIABETES_DAYS = if_else(ANY_DIABETES_INCIDENT == 0 & death_status == 1 & DIABETES_DAYS == DEATH_DAYS - 0.5, DEATH_DAYS, DIABETES_DAYS))

df %>% filter(row_id == "6294") %>% select(row_id, ANY_DIABETES_INCIDENT, DIABETES_DAYS, death_status, DEATH_DAYS, ANY_LIVER_INCIDENT, LIVER_DAYS)

# ------------------------------------------------------------------------------
# 3.4 Adjust Death Dates if Coinciding with Incident Events for Diabetes
# ------------------------------------------------------------------------------

# 2. If date of death = date of any incident diabetes, add epsilon of 0.5 (12 hs) to the days of death counted after baseline (in case both events happen)
df %>% filter(ANY_DIABETES_INCIDENT == 1 & death_status == 1 & DIABETES_DAYS == DEATH_DAYS) %>% nrow() # 55 participants with same date

df <- df %>%
  mutate(DEATH_DAYS = if_else(ANY_DIABETES_INCIDENT == 1 & death_status == 1 & DIABETES_DAYS == DEATH_DAYS,
                              DEATH_DAYS + 0.5,
                              DEATH_DAYS))

# PROBLEM: we also need to adjust LIVER_DAYS for cases with no events (ANY_LIVER_INCIDENT == 0) since it should match DEATH_DAYS now, example:
df %>% filter(row_id == "584") %>% select(row_id, ANY_LIVER_INCIDENT, LIVER_DAYS, death_status, DEATH_DAYS, ANY_DIABETES_INCIDENT, DIABETES_DAYS)

df <- df %>% mutate(LIVER_DAYS = if_else(ANY_LIVER_INCIDENT == 0 & death_status == 1 & LIVER_DAYS == DEATH_DAYS - 0.5, DEATH_DAYS, LIVER_DAYS))

df %>% filter(row_id == "584") %>% select(row_id, ANY_LIVER_INCIDENT, LIVER_DAYS, death_status, DEATH_DAYS, ANY_DIABETES_INCIDENT, DIABETES_DAYS)

# ==============================================================================
# SECTION 4: CREATING FINAL DATASET FOR MULTI-STATE MODELING
# ==============================================================================

# Recode variables as numeric
df$UTI_PREVAL <- as.integer(as.character(df$UTI_PREVAL))
df$ANY_LIVER_INCIDENT <- as.integer(as.character(df$ANY_LIVER_INCIDENT))
df$ANY_DIABETES_INCIDENT <- as.integer(as.character(df$ANY_DIABETES_INCIDENT))
df$death_status <- as.integer(as.character(df$death_status))

# Extract only columns of interest for the multi-state model
df <- df %>% 
  select(row_id, 
  ANY_LIVER_INCIDENT, LIVER_DAYS, # LD event and time
  death_status, DEATH_DAYS, # death event and time
  ANY_DIABETES_INCIDENT, DIABETES_DAYS, # diabetes event and time
  UTI_PREVAL, age_recruitment, sex, ethnicity, household_income, BMI, smoking_status, I10_PREVAL, N18_PREVAL)

# Missing-data handling (we need the real N for the multi-state model)
df <- df %>% drop_na()

str(df)

# ==============================================================================
# SECTION 5: BASIC MODEL WITH HEALTHY, LIVER DISEASE AND DEATH TRANSITIONS
# ==============================================================================

# List the Covariates to keep
covs <- c("age_recruitment", "sex", "ethnicity", "household_income", "BMI", "smoking_status", "I10_PREVAL", "N18_PREVAL", "UTI_PREVAL")

# Create the Transition Matrix
tmat <- trans.illdeath(names=c("Healthy","Liver","Dead"))

# View Matrix
print(tmat)
#          to
# from      Healthy Liver Dead
#   Healthy      NA     1    2
#   Liver        NA    NA    3
#   Dead         NA    NA   NA

# Getting Long Format Data for Multi-State Modelling
mssimdat <- msprep(
    data = df,
    trans= tmat,
    time = c(NA, "LIVER_DAYS", "DEATH_DAYS"),  # Times to transitions
    status = c(NA, "ANY_LIVER_INCIDENT", "death_status"),  # Status indicating the transition
    keep = covs,
    id = "row_id" # column name containing these subject ids
)

# Check First 4 Rows
head(mssimdat, 4)

# Overview of the Transitions
events(mssimdat)
# $Frequencies
#          to
# from      Healthy  Liver   Dead no event total entering
#   Healthy       0   6874  14148   184711         205733
#   Liver         0      0   1554     5320           6874
#   Dead          0      0      0    15702          15702

# $Proportions
#          to
# from         Healthy      Liver       Dead   no event
#   Healthy 0.00000000 0.03341224 0.06876874 0.89781902
#   Liver   0.00000000 0.00000000 0.22606925 0.77393075
#   Dead    0.00000000 0.00000000 0.00000000 1.00000000

# Add Transition-Specific Covariates
mssimdat <- expand.covs(mssimdat, covs, append = TRUE, longnames = FALSE)

colnames(mssimdat) # extension .i refers to transition number i
#  [1] "row_id"              "from"                "to"                  "trans"               "Tstart"              "Tstop"               "time"                "status"             
#  [9] "age_recruitment"     "sex"                 "ethnicity"           "household_income"    "BMI"                 "smoking_status"      "I10_PREVAL"          "N18_PREVAL"         
# [17] "UTI_PREVAL"          "age_recruitment.1"   "age_recruitment.2"   "age_recruitment.3"   "sex.1"               "sex.2"               "sex.3"               "ethnicity1.1"       
# [25] "ethnicity1.2"        "ethnicity1.3"        "ethnicity2.1"        "ethnicity2.2"        "ethnicity2.3"        "ethnicity3.1"        "ethnicity3.2"        "ethnicity3.3"       
# [33] "ethnicity4.1"        "ethnicity4.2"        "ethnicity4.3"        "household_income1.1" "household_income1.2" "household_income1.3" "household_income2.1" "household_income2.2"
# [41] "household_income2.3" "household_income3.1" "household_income3.2" "household_income3.3" "household_income4.1" "household_income4.2" "household_income4.3" "household_income5.1"
# [49] "household_income5.2" "household_income5.3" "BMI.1"               "BMI.2"               "BMI.3"               "smoking_status.1"    "smoking_status.2"    "smoking_status.3"   
# [57] "I10_PREVAL.1"        "I10_PREVAL.2"        "I10_PREVAL.3"        "N18_PREVAL.1"        "N18_PREVAL.2"        "N18_PREVAL.3"        "UTI_PREVAL.1"        "UTI_PREVAL.2"       
# [65] "UTI_PREVAL.3"       

# Calculate Cox with "stratification" on the different possible transitions 
cox_mstat <- coxph(Surv(Tstart,Tstop,status) ~ 
          UTI_PREVAL.1 + UTI_PREVAL.2 + UTI_PREVAL.3 +
          age_recruitment.1 + age_recruitment.2 + age_recruitment.3 + 
          sex.1 + sex.2 + sex.3 +
          ethnicity1.1 + ethnicity1.2 + ethnicity1.3 + ethnicity2.1 + ethnicity2.2 + ethnicity2.3 + ethnicity3.1 + ethnicity3.2 + ethnicity3.3 + ethnicity4.1 + ethnicity4.2 + ethnicity4.3 +
          household_income1.1 + household_income1.2 + household_income1.3 + household_income2.1 + household_income2.2 + household_income2.3 + household_income3.1 + household_income3.2 + household_income3.3 + household_income4.1 + household_income4.2 + household_income4.3 + household_income5.1 + household_income5.2 + household_income5.3 +
          BMI.1 + BMI.2 + BMI.3 +
          smoking_status.1 + smoking_status.2 + smoking_status.3 +
          I10_PREVAL.1 + I10_PREVAL.2 + I10_PREVAL.3 +
          N18_PREVAL.1 + N18_PREVAL.2 + N18_PREVAL.3 +
          strata(trans) + cluster(row_id), data=mssimdat, method="breslow")

# Check Final Rows for Primary Predictor UTI_PREVAL
summary(cox_mstat)$coefficients[1:3, ] %>% as.data.frame() %>% mutate(`Pr(>|z|)` = round(`Pr(>|z|)`, 5))
#                   coef exp(coef)   se(coef)  robust se        z Pr(>|z|)
# UTI_PREVAL.1 0.4720566  1.603288 0.08169161 0.08182864 5.768844  0.00000
# UTI_PREVAL.2 0.3951835  1.484657 0.05718321 0.05779848 6.837264  0.00000
# UTI_PREVAL.3 0.3431128  1.409328 0.15774745 0.17223150 1.992161  0.04635

# Check All Results and Export to .txt
summary(cox_mstat)

sink("Multi_State_Model_Liver_Disease_Death_Healthy_Results_GP_Registered_Participants.txt")
print(summary(cox_mstat))
sink()

# Export Summary of Multi-State Cox Model
broom::tidy(cox_mstat) %>% write_csv("Multi_State_Model_Liver_Disease_Death_Healthy_Results_GP_Registered_Participants.csv")

# ==============================================================================
# SECTION 6: FINAL MODEL WITH HEALTHY, LIVER DISEASE, DIABETES AND DEATH TRANSITIONS
# ==============================================================================

# Create the Transition Matrix
tmat <- transMat(x = list(c(2, 3, 4), c(), c(2, 4), c()), 
                 names = c("Healthy", "Liver", "Diabetes", "Dead"))

# View Matrix
print(tmat)
#           to
# from       Healthy Liver Diabetes Dead
#   Healthy       NA     1        2    3
#   Liver         NA    NA       NA   NA
#   Diabetes      NA     4       NA    5
#   Dead          NA    NA       NA   NA

# Getting Long Format Data for Multi-State Modelling
mssimdat <- msprep(
    data = df,
    trans= tmat,
    time = c(NA, "LIVER_DAYS", "DIABETES_DAYS", "DEATH_DAYS"),  # Times to transitions
    status = c(NA, "ANY_LIVER_INCIDENT", "ANY_DIABETES_INCIDENT", "death_status"),  # Status indicating the transition
    keep = covs,
    id = "row_id" # column name containing these subject ids
)

# Check First 4 Rows
head(mssimdat, 4)

# Overview of the Transitions
events(mssimdat)
# $Frequencies
#           to
# from       Healthy  Liver Diabetes   Dead no event total entering
#   Healthy        0   6314    10491  12750   176178         205733
#   Liver          0      0        0      0     6874           6874
#   Diabetes       0    560        0   1398     8533          10491
#   Dead           0      0        0      0    14148          14148

# $Proportions
#           to
# from          Healthy      Liver   Diabetes       Dead   no event
#   Healthy  0.00000000 0.03069026 0.05099328 0.06197353 0.85634293
#   Liver    0.00000000 0.00000000 0.00000000 0.00000000 1.00000000
#   Diabetes 0.00000000 0.05337909 0.00000000 0.13325708 0.81336384
#   Dead     0.00000000 0.00000000 0.00000000 0.00000000 1.00000000

# Add Transition-Specific Covariates
mssimdat <- expand.covs(mssimdat, covs, append = TRUE, longnames = FALSE)

colnames(mssimdat) # extension .i refers to transition number i
#  [1] "row_id"              "from"                "to"                  "trans"               "Tstart"              "Tstop"               "time"                "status"             
#  [9] "age_recruitment"     "sex"                 "ethnicity"           "household_income"    "BMI"                 "smoking_status"      "I10_PREVAL"          "N18_PREVAL"         
# [17] "UTI_PREVAL"          "age_recruitment.1"   "age_recruitment.2"   "age_recruitment.3"   "age_recruitment.4"   "age_recruitment.5"   "sex.1"               "sex.2"              
# [25] "sex.3"               "sex.4"               "sex.5"               "ethnicity1.1"        "ethnicity1.2"        "ethnicity1.3"        "ethnicity1.4"        "ethnicity1.5"       
# [33] "ethnicity2.1"        "ethnicity2.2"        "ethnicity2.3"        "ethnicity2.4"        "ethnicity2.5"        "ethnicity3.1"        "ethnicity3.2"        "ethnicity3.3"       
# [41] "ethnicity3.4"        "ethnicity3.5"        "ethnicity4.1"        "ethnicity4.2"        "ethnicity4.3"        "ethnicity4.4"        "ethnicity4.5"        "household_income1.1"
# [49] "household_income1.2" "household_income1.3" "household_income1.4" "household_income1.5" "household_income2.1" "household_income2.2" "household_income2.3" "household_income2.4"
# [57] "household_income2.5" "household_income3.1" "household_income3.2" "household_income3.3" "household_income3.4" "household_income3.5" "household_income4.1" "household_income4.2"
# [65] "household_income4.3" "household_income4.4" "household_income4.5" "household_income5.1" "household_income5.2" "household_income5.3" "household_income5.4" "household_income5.5"
# [73] "BMI.1"               "BMI.2"               "BMI.3"               "BMI.4"               "BMI.5"               "smoking_status.1"    "smoking_status.2"    "smoking_status.3"   
# [81] "smoking_status.4"    "smoking_status.5"    "I10_PREVAL.1"        "I10_PREVAL.2"        "I10_PREVAL.3"        "I10_PREVAL.4"        "I10_PREVAL.5"        "N18_PREVAL.1"       
# [89] "N18_PREVAL.2"        "N18_PREVAL.3"        "N18_PREVAL.4"        "N18_PREVAL.5"        "UTI_PREVAL.1"        "UTI_PREVAL.2"        "UTI_PREVAL.3"        "UTI_PREVAL.4"       
# [97] "UTI_PREVAL.5"       

# Calculate Cox with "stratification" on the different possible transitions 
cox_mstat <- coxph(Surv(Tstart,Tstop,status) ~ 
          UTI_PREVAL.1 + UTI_PREVAL.2 + UTI_PREVAL.3 + UTI_PREVAL.4 + UTI_PREVAL.5 +
          age_recruitment.1 + age_recruitment.2 + age_recruitment.3 + age_recruitment.4 + age_recruitment.5 +
          sex.1 + sex.2 + sex.3 + sex.4 + sex.5 +
          ethnicity1.1 + ethnicity1.2 + ethnicity1.3 + ethnicity1.4 + ethnicity1.5 +
          ethnicity2.1 + ethnicity2.2 + ethnicity2.3 + ethnicity2.4 + ethnicity2.5 + 
          ethnicity3.1 + ethnicity3.2 + ethnicity3.3 + ethnicity3.4 + ethnicity3.5 + 
          ethnicity4.1 + ethnicity4.2 + ethnicity4.3 + ethnicity4.4 + ethnicity4.5 +
          household_income1.1 + household_income1.2 + household_income1.3 + household_income1.4 + household_income1.5 +
          household_income2.1 + household_income2.2 + household_income2.3 + household_income2.4 + household_income2.5 +
          household_income3.1 + household_income3.2 + household_income3.3 + household_income3.4 + household_income3.5 +
          household_income4.1 + household_income4.2 + household_income4.3 + household_income4.4 + household_income4.5 +
          household_income5.1 + household_income5.2 + household_income5.3 + household_income5.4 + household_income5.5 +
          BMI.1 + BMI.2 + BMI.3 + BMI.4 + BMI.5 +
          smoking_status.1 + smoking_status.2 + smoking_status.3 + smoking_status.4 + smoking_status.5 +
          I10_PREVAL.1 + I10_PREVAL.2 + I10_PREVAL.3 + I10_PREVAL.4 + I10_PREVAL.5 +
          N18_PREVAL.1 + N18_PREVAL.2 + N18_PREVAL.3 + N18_PREVAL.4 + N18_PREVAL.5 +
          strata(trans) + cluster(row_id), data=mssimdat, method="breslow")

# Check Final Rows for Primary Predictor UTI_PREVAL
summary(cox_mstat)$coefficients[1:5, ] %>% as.data.frame() %>% mutate(`Pr(>|z|)` = round(`Pr(>|z|)`, 5))
#                   coef exp(coef)   se(coef)  robust se         z Pr(>|z|)
# UTI_PREVAL.1 0.4827189  1.620474 0.08567197 0.08577517 5.6277237  0.00000
# UTI_PREVAL.2 0.3628978  1.437489 0.06815061 0.07145826 5.0784578  0.00000
# UTI_PREVAL.3 0.3602716  1.433719 0.06208254 0.06234198 5.7789557  0.00000
# UTI_PREVAL.4 0.2476175  1.280970 0.27173197 0.27292412 0.9072761  0.36426
# UTI_PREVAL.5 0.5005726  1.649666 0.14730815 0.15602078 3.2083717  0.00133

# Check All Results and Export to .txt
summary(cox_mstat)

sink("Multi_State_Model_Liver_Disease_Death_Diabetes_Healthy_Results_GP_Registered_Participants.txt")
print(summary(cox_mstat))
sink()

# Export Summary of Multi-State Cox Model
broom::tidy(cox_mstat) %>% write_csv("Multi_State_Model_Liver_Disease_Death_Diabetes_Healthy_Results_GP_Registered_Participants.csv")
```