# Project Documentation & Pipeline: Liver Disease & UTI in UK Biobank

This document describes the R-based processing pipeline for the UK Biobank (UKB) Liver Disease & UTI Analysis project.

Objective: Use cox models for time-to-event (survival) models connecting urinary tract infections (UTIs) to incident liver disease (LD). In addition, check various covariates distributions.

Key input is the file: `ukb_analysis_liver_uti_inc_prev_23Jan2026.rds` from `01_extract_data_UKB.md`. Content: UKB-derived dataset prepared by me with relevant variables for the analysis. 

It contains n=486,748 participants with data on UTIs, liver disease, demographics, lifestyle, comorbidities, and biomarkers. The variables still need some processing before the analysis.

# R Script: Final Data Preparation and Cox Models

```R
# ==============================================================================
# 1. GLOBAL CONFIGURATION & LIBRARIES 
# ==============================================================================
library(tidyverse) # data manipulation and visualization
library(report)  # automated reporting of descriptive statistics
library(survival) # to run Cox models
packageVersion("survival") # 3.8.6
options(width = 200, scipen = 999) # improve console readability

# ==============================================================================
# SECTION 2: DATA LOADING & COHORT FILTRATION
# ==============================================================================

# Load the file produced in 01_extract_data_UKB.md
df <- readRDS("ukb_analysis_liver_uti_inc_prev_23Jan2026.rds")
dim(df) # 486748    117

# Check how many prevalent liver disease and viral hepatitis cases at baseline
df %>% summarise(n_liver_preval = sum(ANY_LIVER_PREVAL == 1, na.rm = TRUE), n_viral_hep_preval = sum(ANY_VIRAL_HEP_PREVAL == 1, na.rm = TRUE))
#   n_liver_preval n_viral_hep_preval
# 1           4738               2837

df %>% filter(ANY_LIVER_PREVAL == 1 | ANY_VIRAL_HEP_PREVAL == 1) %>% nrow() # 7202 participants with either condition

# Exclude prevalent liver disease and viral hepatitis cases at baseline
df <- df %>% filter(ANY_LIVER_PREVAL == 0 & ANY_VIRAL_HEP_PREVAL == 0)
dim(df) # 479546    117

# ==============================================================================
# EXTRA SECTION : save these outputs in a separate folder
# ==============================================================================

dir.create("Outputs_UTI_to_LD_analysis", showWarnings = FALSE)
setwd("Outputs_UTI_to_LD_analysis")

# ==============================================================================
# SECTION 3: DEFINING FOLLOW-UP TIME VARIABLES
# ==============================================================================

# Define Study Census Date 
admin_end <- as.Date("2023-09-01")

# Establish the End of Follow-up Date 
# Rationale: The date the participant is no longer "under observation" for non-event reasons
# Includes: Death or the administrative end of the study
df$date_end_of_followup <- pmin(df$death_date, admin_end, na.rm = TRUE)

# Verification: Check distribution
summary(as.numeric(df$date_end_of_followup - df$date_attendance) / 365.25)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#  0.01095 13.72758 14.47228 14.07830 15.18686 17.47023 

# ==============================================================================
# SECTION 4: RECODING & CLEANING ALL VARIABLES
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Define Date for Incident Liver (ICD-10 Codes K70-K77)
# ------------------------------------------------------------------------------

# Identify the Incident Event Date (earliest occurrence of any liver-related ICD-10 code K70-K77)
df$ANY_LIVER_INCIDENT_DATE <- do.call(pmin, c(df[paste0("K", 70:77, "_date")], na.rm = TRUE))

# Verification: Does anyone have an incident date AFTER their follow-up ended?
df %>% filter(ANY_LIVER_INCIDENT == 1 & ANY_LIVER_INCIDENT_DATE > date_end_of_followup) %>%  nrow()
# O, no participants with incident dates after follow-up end

# Summary table of incident cases for each K-code
df %>% 
  filter(ANY_LIVER_INCIDENT == 1) %>% 
  summarise(Total_Incident_Cases = n(), across(paste0("K", 70:77, "_INCIDENT"), ~sum(as.numeric(as.character(.x)), na.rm = TRUE)))
#   Total_Incident_Cases K70_INCIDENT K71_INCIDENT K72_INCIDENT K73_INCIDENT K74_INCIDENT K75_INCIDENT K76_INCIDENT K77_INCIDENT
# 1                16561         1509          154         1056          118         2001         1534        13752           24

# ------------------------------------------------------------------------------
# 4.2 Categorical Variable Recoding & Factoring
# ------------------------------------------------------------------------------

df <- df %>%
  mutate(
    # 1. Sex
    sex = factor(recode(sex, "Female" = 0, "Male" = 1)),

    # 2. Ethnicity (Aggregated into 5 major clinical groups)
    ethnicity = case_when(
      ethnicity %in% c("British", "Irish", "White", "Any other white background") ~ "White",
      ethnicity %in% c("White and Black Caribbean", "White and Black African", "White and Asian", "Any other mixed background", "Mixed") ~ "Mixed",
      ethnicity %in% c("Chinese", "Indian", "Pakistani", "Bangladeshi", "Any other Asian background", "Asian or Asian British") ~ "Asian",
      ethnicity %in% c("Caribbean", "African", "Any other Black background", "Black or Black British") ~ "Black",
      ethnicity %in% c("Other ethnic group") ~ "Other",
      TRUE ~ NA_character_ # Set "Do not know" and "Prefer not to answer" as NA
    ),
    ethnicity = factor(ethnicity, levels = c("White", "Asian", "Black", "Mixed", "Other")),

    # 3. Smoking Status
    smoking_status = case_when(
      smoking_status %in% c("Yes, on most or all days", "Only occasionally") ~ "Current",
      smoking_status == "No" ~ "No",
      TRUE ~ NA_character_
    ),
    smoking_status = factor(smoking_status, levels = c("No", "Current")),
    
    # 4. Alcohol Frequency (Use "One to three times a month" as reference)
    alcohol_frequency = na_if(as.character(alcohol_frequency), "Prefer not to answer"),
    alcohol_frequency = factor(alcohol_frequency, 
                               levels = c("One to three times a month", "Never", "Special occasions only", "Once or twice a week", "Three or four times a week", "Daily or almost daily")),
    
    # 5. Household Income (Grouping missingness into 'Unknown')
    household_income = case_when(
      household_income %in% c("Prefer not to answer", "Do not know") | is.na(household_income) ~ "Unknown",
      TRUE ~ as.character(household_income)
    ),
    household_income = factor(household_income, 
                              levels = c("Less than 18,000", "18,000 to 30,999", "31,000 to 51,999", "52,000 to 100,000", "Greater than 100,000", "Unknown")),
    
    # 6. Physical Activity (IPAQ Scoring Protocol)
    physical_activity = case_when(
      !is.finite(physical_activity) ~ NA_character_,
      physical_activity < 600 ~ "<600",
      physical_activity < 3000 ~ "600–<3000",
      physical_activity >= 3000 ~ "≥3000"
    ),
    physical_activity = factor(physical_activity, levels = c("<600", "600–<3000", "≥3000")),
    
    # 7. BMI Categories (WHO standard)
    BMI_category = case_when(
      is.na(BMI) ~ NA_character_,
      BMI < 18.5 ~ "<18.5",
      BMI < 25   ~ "18.5–<25",
      BMI < 30   ~ "25–<30",
      BMI >= 30  ~ "≥30"
    ),
    BMI_category = factor(BMI_category, levels = c("18.5–<25", "<18.5", "25–<30", "≥30")) # Reference set to healthy
  )

  # Recode gp_reg (Na=No, 1 or >1=Yes) for sensitivity analysis 5 (GP registered subjects only)
df$gp_reg <- case_when(is.na(df$gp_reg) ~ "No", df$gp_reg >= 1 ~ "Yes")
df$gp_reg <- factor(df$gp_reg, levels = c("No", "Yes"))
table(df$gp_reg, useNA = "ifany")

# Check the new variables distributions with report package
report::report(df %>% select(sex, ethnicity, smoking_status, alcohol_frequency, household_income, physical_activity, BMI_category, gp_reg))
# The data contains 479546 observations of the following 7 variables:
#   - sex: 2 levels, namely 0 (n = 260434, 54.31%) and 1 (n = 219112, 45.69%)
#   - ethnicity: 5 levels, namely White (n = 451990, 94.25%), Asian (n = 10713, 2.23%), Black (n = 7532, 1.57%), Mixed (n = 2799, 0.58%), Other (n = 4244, 0.89%) and missing (n = 2268, 0.47%)
#   - smoking_status: 2 levels, namely No (n = 428599, 89.38%), Current (n = 50040, 10.43%) and missing (n = 907, 0.19%)
#   - alcohol_frequency: 6 levels, namely One to three times a month (n = 53318, 11.12%), Never (n = 37945, 7.91%), Special occasions only (n = 54955, 11.46%), Once or twice a week (n = 123779, 25.81%),
# Three or four times a week (n = 110944, 23.14%), Daily or almost daily (n = 97555, 20.34%) and missing (n = 1050, 0.22%)
#   - household_income: 6 levels, namely Less than 18,000 (n = 92713, 19.33%), 18,000 to 30,999 (n = 104180, 21.72%), 31,000 to 51,999 (n = 106882, 22.29%), 52,000 to 100,000 (n = 83350, 17.38%), Greater
# than 100,000 (n = 22137, 4.62%) and Unknown (n = 70284, 14.66%)
#   - physical_activity: 3 levels, namely <600 (n = 67955, 14.17%), 600–<3000 (n = 187807, 39.16%), ≥3000 (n = 114940, 23.97%) and missing (n = 108844, 22.70%)
#   - BMI_category: 4 levels, namely 18.5–<25 (n = 155911, 32.51%), <18.5 (n = 2461, 0.51%), 25–<30 (n = 203307, 42.40%), ≥30 (n = 115935, 24.18%) and missing (n = 1932, 0.40%)
#   - gp_reg: 2 levels, namely No (n = 260947, 54.42%) and Yes (n = 218599, 45.58%)

# Save the analytic-ready dataset
saveRDS(df, file = "ukb_analysis_liver_uti_final_cohort_479546_participants_27Jan2026.rds")

# ==============================================================================
# SECTION 5: ASSESSING DATA ATTRITION ACROSS MODELS
# ==============================================================================

library(dtrackr) # for tracking data attrition

#-----------------------------
# 5.1 Stepwise attrition table for Models 0 to 4 
#-----------------------------

# Note: the model saved in the end will be Model 4 
df_m4 <- df %>% 
  track() %>%
  comment("UTI prevalent status") %>% filter(!is.na(UTI_PREVAL)) %>%
  comment("Age at recruitment") %>% filter(!is.na(age_recruitment)) %>%
  comment("Sex") %>% filter(!is.na(sex)) %>%
  comment("Ethnicity") %>% filter(!is.na(ethnicity)) %>%
  comment("Household income") %>% filter(!is.na(household_income)) %>%
  comment("BMI") %>% filter(!is.na(BMI)) %>%
  comment("Smoking status") %>% filter(!is.na(smoking_status)) %>%
  comment("Alcohol frequency") %>% filter(!is.na(alcohol_frequency)) %>%
  comment("Prevalent diabetes\n(ICD-10 E10-E14)") %>% filter(!is.na(ANY_DIABETES_PREVAL)) %>%
  comment("Prevalent hypertension\n(ICD-10 I10)") %>% filter(!is.na(I10_PREVAL)) %>%
  comment("Prevalent CKD\n(ICD-10 N18)") %>% filter(!is.na(N18_PREVAL)) 

dim(df_m4) # 474719  (final dataset for model 4)

df_m4 %>% flowchart(filename = "data_attrition_flowchart_model_1_to_4.html", format = "pdf", fontsize = 8)

#-----------------------------
# 5.2 Stepwise attrition table for models Sensitivity 1 to 4 
# Note: these are additional models beyond Model 4 or modifications of it
#-----------------------------

# Sensitivity 1 (physical activity)
df_m4 %>%
  track() %>%
  comment("Physical activity") %>% filter(!is.na(physical_activity)) %>% 
  flowchart(filename = "data_attrition_flowchart_sensitivity_1.html", format = "pdf", fontsize = 8)

# Sensitivity 2 (biomarkers)
df_m4 %>%
  track() %>%
  comment("Creatinine") %>% filter(!is.na(creatinine)) %>% 
  comment("C-reactive protein") %>% filter(!is.na(c_reactive_protein)) %>%
  comment("HDL cholesterol") %>% filter(!is.na(hdl_cholesterol)) %>%
  comment("Triglycerides") %>% filter(!is.na(triglycerides)) %>%
  flowchart(filename = "data_attrition_flowchart_sensitivity_2.html", format = "pdf", fontsize = 8)

# Sensitivity 3 (alcohol numeric)
 df %>% 
  track() %>%
  comment("UTI prevalent status") %>% filter(!is.na(UTI_PREVAL)) %>%
  comment("Age at recruitment") %>% filter(!is.na(age_recruitment)) %>%
  comment("Sex") %>% filter(!is.na(sex)) %>%
  comment("Ethnicity") %>% filter(!is.na(ethnicity)) %>%
  comment("Household income") %>% filter(!is.na(household_income)) %>%
  comment("BMI") %>% filter(!is.na(BMI)) %>%
  comment("Smoking status") %>% filter(!is.na(smoking_status)) %>%
  comment("Alcohol g (continuous)") %>% filter(!is.na(alcohol_g)) %>%
  comment("Prevalent diabetes\n(ICD-10 E10-E14)") %>% filter(!is.na(ANY_DIABETES_PREVAL)) %>%
  comment("Prevalent hypertension\n(ICD-10 I10)") %>% filter(!is.na(I10_PREVAL)) %>%
  comment("Prevalent CKD\n(ICD-10 N18)") %>% filter(!is.na(N18_PREVAL)) %>%
  flowchart(filename = "data_attrition_flowchart_sensitivity_3.html", format = "pdf", fontsize = 8)

# Sensitivity 4 (current drinkers)
df_m4 %>%
  track() %>%
  comment("Current drinkers only") %>% filter(!alcohol_frequency %in% c("Never", "Special occasions only")) %>%
  flowchart(filename = "data_attrition_flowchart_sensitivity_4.html", format = "pdf", fontsize = 8)

# Sensitivity 5 (GP registered subjects only)
df_m4 %>%
  track() %>%
  comment("GP registered subjects only") %>% filter(gp_reg == "Yes") %>%
  flowchart(filename = "data_attrition_flowchart_sensitivity_5.html", format = "pdf", fontsize = 8)

# ==============================================================================
# SECTION 6: CREATE TABLE WITH VARIABLES DISTRIBUTION (N=479,546) BY UTI STATUS
# ==============================================================================

library(gtsummary)
library(readr)
library(forcats)

as01 <- function(x) parse_integer(as.character(x))

#-----------------------------
# 6.1 Prep/derive variables for the exact rows we want
#-----------------------------

df_tbl1 <- df %>%
  mutate(
    UTI_PREVAL = factor(as01(UTI_PREVAL),
                        levels = c(0, 1),
                        labels = c("No prevalent UTI", "Prevalent UTI")),

    # Men as a single dichotomous row
    men = factor(if_else(as01(sex) == 1, "Yes", "No"), levels = c("No", "Yes")),

    # Current smoking as a single dichotomous row
    current_smoking = factor(if_else(smoking_status == "Current", "Yes", "No"),
                             levels = c("No", "Yes")),

    # Physical activity with the exact level labels & order
    physical_activity = factor(
      physical_activity,
      levels = c("<600", "600–<3000", "≥3000"),
      labels = c("Low (<600)", "Moderate (600–<3000)", "High (≥3000)")
    ),

    # Make all 0/1 disease variables consistent Yes/No
    across(
      c(ANY_LIVER_INCIDENT, K70_INCIDENT, K71_INCIDENT, K72_INCIDENT, K73_INCIDENT,
        K74_INCIDENT, K75_INCIDENT, K76_INCIDENT, K77_INCIDENT,
        ANY_DIABETES_PREVAL, E11_PREVAL, I10_PREVAL, N18_PREVAL, death_status),
      ~ factor(as01(.x), levels = c(0, 1), labels = c("No", "Yes"))
    )
  )

#-----------------------------
# 6.2 Force the variable order 
#-----------------------------

# Variable order
var_order <- c("age_recruitment", "men", "BMI", "ethnicity", "household_income", "alcohol_frequency", "current_smoking",
                "ANY_LIVER_INCIDENT", "K70_INCIDENT", "K71_INCIDENT", "K72_INCIDENT", "K73_INCIDENT", "K74_INCIDENT",
                "K75_INCIDENT", "K76_INCIDENT", "K77_INCIDENT", "ANY_DIABETES_PREVAL", "E11_PREVAL", "I10_PREVAL",
                "N18_PREVAL", "physical_activity", "creatinine", "c_reactive_protein", "hdl_cholesterol", "triglycerides", "death_status")

# Dichotomous rows (single-line Yes counts)
dicho_rows <- c("men", "current_smoking", "ANY_LIVER_INCIDENT", "K70_INCIDENT", "K71_INCIDENT", "K72_INCIDENT", "K73_INCIDENT", "K74_INCIDENT", 
                "K75_INCIDENT", "K76_INCIDENT", "K77_INCIDENT", "ANY_DIABETES_PREVAL", "E11_PREVAL", "I10_PREVAL", "N18_PREVAL", "death_status")

#-----------------------------
# 6.3 Build Table 1 with exact labels 
#-----------------------------

tbl1 <- df_tbl1 %>%
  select(UTI_PREVAL, all_of(var_order)) %>%
  tbl_summary(
    by = UTI_PREVAL,
    type  = list(all_of(dicho_rows) ~ "dichotomous"),
    value = list(all_of(dicho_rows) ~ "Yes"),
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)",
      all_dichotomous() ~ "{n} ({p}%)"
    ),
    label = list(
      age_recruitment     ~ "Age at enrollment, median (IQR), years",
      men                ~ "Men, n (%)",
      BMI                ~ "BMI, median (IQR), kg/m2",
      ethnicity          ~ "Ethnicity, n (%)",
      household_income   ~ "Household income, n (%)",
      alcohol_frequency  ~ "Alcohol intake, n (%)",
      current_smoking    ~ "Current smoking, n (%)",
      ANY_LIVER_INCIDENT ~ "Any liver disease (K70-77)",
      K70_INCIDENT       ~ "Alcoholic liver disease (K70)",
      K71_INCIDENT       ~ "Toxic liver disease (K71)",
      K72_INCIDENT       ~ "Hepatic failure (K72) incidence",
      K73_INCIDENT       ~ "Chronic hepatitis, not elsewhere classified (K73)",
      K74_INCIDENT       ~ "Fibrosis/cirrhosis of liver (K74)",
      K75_INCIDENT       ~ "Other inflammatory liver diseases (K75)",
      K76_INCIDENT       ~ "Other diseases of liver (K76)",
      K77_INCIDENT       ~ "Liver disorders in diseases classified elsewhere (K77)",
      ANY_DIABETES_PREVAL~ "Any diabetes (E10-E14)",
      E11_PREVAL         ~ "Type 2 diabetes mellitus (E11)",
      I10_PREVAL         ~ "Prevalence of essential (primary) hypertension (I10)",
      N18_PREVAL         ~ "Prevalence of chronic kidney disease (N18)",
      physical_activity  ~ "Physical activity (total metabolic equivalent of task, min/week)",
      creatinine         ~ "Serum creatinine, median (IQR), mmol/L",
      c_reactive_protein ~ "C-reactive protein, median (IQR), mg/L",
      hdl_cholesterol    ~ "HDL cholesterol, median (IQR), mmol/L",
      triglycerides      ~ "Triglycerides, median (IQR), mmol/L",
      death_status       ~ "Deaths during follow-up"
    ),
    missing = "ifany"
  ) %>%
  bold_labels()

# Export table as CSV
write.csv(as_tibble(tbl1, col_labels = TRUE), "table1_479546_participants_by_uti.csv", row.names = FALSE)

# ==============================================================================
# SECTION 7: FUNCTIONS TO RUN COX MODELS
# ==============================================================================

# Main function, runs Cox models with different follow-up horizons and covariates
run_cox_model <- function(df, outcome, outcome_date, primary_predictor, covariates = NULL) {  
  # 1. Global Data Setup 
  cat("==============================\n")
  cat("Running Cox model for outcome:", outcome, "and primary predictor:", primary_predictor, "\n")
  cat("--------- Covariates:", ifelse(is.null(covariates), "None (Crude model)", paste(covariates, collapse = ", ")), "\n")
  
  # Create working copy to avoid side effects
  df_clean <- df 

  # Remove missing data from covariates and primary predictor
  if(!is.null(covariates)) {
    df_clean <- df_clean %>% drop_na(all_of(c(covariates, primary_predictor)))
  } else {
    df_clean <- df_clean %>% drop_na(all_of(primary_predictor))
  }
  
  # Standardize primary predictor (the main exposure)
  df_clean$primary_predictor <- as.integer(as.character(df_clean[[primary_predictor]]))
  
  # Create Event/Time columns once (Full follow-up)
  df_clean$event_status <- as.integer(as.character(df_clean[[outcome]]))
  
  # Calculate full follow-up end date
  # Logic: if event occurred, use event date, else use end_of_followup
  df_clean$final_date <- dplyr::if_else(df_clean$event_status == 1, df_clean[[outcome_date]], df_clean$date_end_of_followup)
  
  # Calculate baseline raw days
  df_clean$time_days_full <- as.numeric(df_clean$final_date - df_clean$date_attendance)
  
  # 2. Helper Function 
  fit_horizon <- function(data_model, cutoff_years, label) {
    
    # Efficiently truncate time and status without deep copying dates
    # If cutoff is infinite (Full model), use original data
    if (is.infinite(cutoff_years)) {
      data_model$time <- data_model$time_days_full
      data_model$status <- data_model$event_status
      cutoff_days <- Inf
    } else {
      cutoff_days <- cutoff_years * 365.25
      # Status is 1 ONLY if event happened AND it happened before cutoff
      data_model$status <- ifelse(data_model$event_status == 1 & data_model$time_days_full <= cutoff_days, 1, 0)
      # Time is the minimum of actual time or cutoff
      data_model$time <- pmin(data_model$time_days_full, cutoff_days)
    }

    # Define Formula
    fml <- stats::reformulate(termlabels = c("primary_predictor", covariates), response = "Surv(time, status)")
    
    # Fit Model, use 'data_model' here which is the local subset/modified version
    model <- coxph(fml, data = data_model)
    
    # Tidy results
    model_tidy <- broom::tidy(model, exp = TRUE, conf.int = TRUE)
    
    # Add metadata
    model_tidy$n_primary_predictor <- sum(data_model[[primary_predictor]] == 1, na.rm = TRUE)
    model_tidy$n_outcome <- nobs(model)
    model_tidy$n_participants <- model$n
    
    # Check PH Assumption
    zph <- cox.zph(model)
    violated <- zph$table %>% 
      as.data.frame() %>% 
      rownames_to_column("var") %>% 
      filter(p < 0.05) %>% 
      pull(var)
    
    if(length(violated) == 0) {
      ph_msg <- "PH assumption met"
    } else {
      ph_msg <- paste("PH assumption violated for:", paste(violated, collapse = ", "))
    }
    model_tidy$PH_assumption <- ph_msg
    
    cat("--------- ", label, ": ", ph_msg, "\n", sep = "")
    
    return(list(model = model, model_summary = summary(model), model_tidy = model_tidy, ph_assumption = ph_msg))
  }
  
  # 3. Execution  
  # Define horizons: Name = Years
  horizons <- list("full" = Inf, "horizon_5y" = 5, "horizon_10y" = 10)
  
  # Run all models using lapply
  results <- lapply(names(horizons), function(h_name) {
    fit_horizon(df_clean, horizons[[h_name]], h_name)
  })
  names(results) <- names(horizons)
  
  return(results)
}

# Function to extract tidy results for main predictor from model results
extract_main_predictor_results <- function(res, model_name) {
  bind_rows(
    res$full$model_tidy %>% mutate(follow_up = "Full"),
    res$horizon_5y$model_tidy  %>% mutate(follow_up = "5-year"),
    res$horizon_10y$model_tidy %>% mutate(follow_up = "10-year")
  ) %>%
    filter(term == "primary_predictor") %>% 
    mutate(model = model_name,
          PH_meeting = ifelse(PH_assumption == "PH assumption met", "Yes", "No"),
          estimate = round(estimate, 2),
          conf.low = round(conf.low, 2),
          conf.high = round(conf.high, 2),
          p.value = formatC(p.value, format = "e", digits = 2),
          HR_CI = paste0(estimate, " (", conf.low, "–", conf.high, ")")) %>%
    select(model, follow_up, HR_CI, p.value, n_participants, n_primary_predictor, n_outcome, PH_meeting, PH_assumption)
    }

# Function to create basic plot of Cumulative Incidence (will work with all horizons)
cumulative_incidence_plot <- function(df, outcome, outcome_date, primary_predictor, max_years = 17.5, 
                                      legend_title = "UTI status", legend_labs = c("No prevalent UTI", "Prevalent UTI"),
                                      outcome_name_y = "incident liver disease (K70–K77)") {  

  # 1. Setup Data once (similar to run_cox_model function)
  # Create working copy to avoid side effects
  df_clean <- df 
  
  # Standardize primary predictor (the main exposure)
  df_clean$primary_predictor <- as.integer(as.character(df_clean[[primary_predictor]]))
  
  # Create Event/Time columns once (Full follow-up)
  df_clean$event_status <- as.integer(as.character(df_clean[[outcome]]))
  
  # Calculate full follow-up end date
  # Logic: if event occurred, use event date, else use end_of_followup
  df_clean$final_date <- dplyr::if_else(df_clean$event_status == 1, df_clean[[outcome_date]], df_clean$date_end_of_followup)
  
  # Calculate baseline raw days
  df_clean$time_days_full <- as.numeric(df_clean$final_date - df_clean$date_attendance)

  # Fit survival curve once (for full follow-up)
  fit <- survfit(Surv(time_days_full, event_status) ~ primary_predictor, data = df_clean)

  # 2. Create plot to be modified for each horizon later
  # this stopped working due to updates with ggplot2 and survminer
  # I am leaving it here commented because it is the best way to create these plots
  # library(ggsurvfit)
  # p <- survfit2(Surv(time_days_full, event_status) ~ primary_predictor, data = df_clean) %>%
  #       ggsurvfit(type = "risk") +
  #       add_confidence_interval() +
  #       add_risktable(size = 4) +
  #       scale_ggsurvfit() 
  
  # alternative using survminer
  library(survminer) 
  p <- ggsurvplot(fit, data = df_clean,
    fun = "event",
    palette = c("#609DEC", "#E599E9"),  
    conf.int = TRUE,
    surv.scale = "percent",
    xscale = 365.25,                          
    xlim = c(-150, max_years * 365.25),         
    break.time.by = 5 * 365.25,    
    xlab = "Years since baseline",
    ylab = paste("Cumulative incidence of \n", outcome_name_y),
    risk.table = TRUE,
    legend.title = legend_title,
    legend.labs = legend_labs,
    risk.table.col = "strata",
    risk.table.y.text = FALSE,
    risk.table.height = 0.25
  )

  # Return plot
  return(p)
}

# Function to save forest plot as SVG
save_forest_svg <- function(fit, file, data = NULL, width = 12, height = 8) {
  grDevices::svg(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  p <- survminer::ggforest(fit, data = data)
  print(p)  # <- crucial inside loops/devices
  invisible(p)
}

# ==============================================================================
# SECTION 8: ANALYSES WITH ANY_LIVER_INCIDENT (ICD-10 K70-K77) AS OUTCOME
# ==============================================================================

#-----------------------------
# 8.0 Preparations of a folder to save results
#-----------------------------

# Create new results directory and set as working directory
dir.create("Results_UTI_PREVAL_to_ANY_LIVER_INCIDENT", showWarnings = FALSE)
setwd("Results_UTI_PREVAL_to_ANY_LIVER_INCIDENT")

#-----------------------------
# 8.1 Check crude model and plot cumulative incidence 
#-----------------------------

# Run function with no covariates (crude model)
res <- run_cox_model(df, outcome = "ANY_LIVER_INCIDENT", outcome_date = "ANY_LIVER_INCIDENT_DATE", primary_predictor = "UTI_PREVAL")

rbind(res$full$model_tidy, res$horizon_5y$model_tidy, res$horizon_10y$model_tidy)
#   term              estimate std.error statistic  p.value conf.low conf.high n_primary_predictor n_outcome n_participants PH_assumption    
# 1 primary_predictor     2.12    0.0495     15.2  2.19e-52     1.93      2.34                6161     16561         479546 PH assumption met
# 2 primary_predictor     2.25    0.0958      8.48 2.19e-17     1.87      2.72                6161      3990         479546 PH assumption met
# 3 primary_predictor     2.29    0.0599     13.9  9.77e-44     2.04      2.58                6161     10275         479546 PH assumption met

# Cumulative Incidence by UTI status
svg("CumInc_UTI_PREVAL_to_ANY_LIVER_INCIDENT.svg", width = 6, height = 6)
cumulative_incidence_plot(df, outcome = "ANY_LIVER_INCIDENT", outcome_date = "ANY_LIVER_INCIDENT_DATE", primary_predictor = "UTI_PREVAL", 
                          max_years = 17.5, 
                          legend_title = "UTI status", legend_labs = c("No prevalent UTI", "Prevalent UTI"),
                          outcome_name_y = "incident liver disease (K70–K77)")

dev.off()

#-----------------------------
# 8.2 Run all models 0-4 and sensitivities 1-4, and the new sensitivity 5 (GP registered subjects only)
#-----------------------------

# Define covariate sets hierarchically
covs_demo   <- c("age_recruitment", "sex")
covs_socio  <- c(covs_demo, "ethnicity", "household_income")
covs_life   <- c(covs_socio, "BMI", "smoking_status", "alcohol_frequency")
covs_comorb <- c(covs_life, "ANY_DIABETES_PREVAL", "I10_PREVAL", "N18_PREVAL")

# Define Sensitivity sets
covs_sens1 <- c(covs_comorb, "physical_activity")
covs_sens2 <- c(covs_comorb, "creatinine", "c_reactive_protein", "hdl_cholesterol", "triglycerides")
covs_sens3 <- c(setdiff(covs_comorb, "alcohol_frequency"), "alcohol_g")

# Prepare specialized dataframe for Sensitivity 4
df_drinkers <- df %>%  
  filter(!alcohol_frequency %in% c("Never", "Special occasions only")) %>% 
  mutate(across(where(is.factor), droplevels))

# Prepare specialized dataframe for Sensitivity 5 (GP registered subjects only)
df_gp_registered <- df %>%
  filter(gp_reg == "Yes") %>%
  mutate(across(where(is.factor), droplevels))

# List of models to run
# Format: list(Name, Covariates, Dataframe(optional))
model_runs <- list(
  list(name = "Model_0_Crude",              covs = NULL,        data = df),
  list(name = "Model_1_Demographics",       covs = covs_demo,   data = df),
  list(name = "Model_2_Sociodemographics",  covs = covs_socio,  data = df),
  list(name = "Model_3_Lifestyle",          covs = covs_life,   data = df),
  list(name = "Model_4_Comorbidity",        covs = covs_comorb, data = df),
  list(name = "Sens_1_Physical_activity",   covs = covs_sens1,  data = df),
  list(name = "Sens_2_Biomarkers",          covs = covs_sens2,  data = df),
  list(name = "Sens_3_Alcohol_numeric",     covs = covs_sens3,  data = df),
  list(name = "Sens_4_Current_drinkers",    covs = covs_comorb, data = df_drinkers),
  list(name = "Sens_5_GP_registered",       covs = covs_comorb, data = df_gp_registered)
)

# Run all models and store results
results_all_models <- data.frame()

for (model_info in model_runs) {

  # Run Cox for this model
  res <- run_cox_model(
    df = model_info$data,
    outcome = "ANY_LIVER_INCIDENT",
    outcome_date = "ANY_LIVER_INCIDENT_DATE",
    primary_predictor = "UTI_PREVAL",
    covariates = model_info$covs
  )

  # Save individual model result
  res$model <- NULL # drop model object before saving (keeps only summaries)
  saveRDS(res, paste0("Cox_", model_info$name, "_UTI_PREVAL_to_ANY_LIVER_INCIDENT.rds"))

  # Create forest plots
  save_forest_svg(
    res$full$model,
    paste0("Cox_", model_info$name, "_UTI_PREVAL_to_ANY_LIVER_INCIDENT_full_forest_plot.svg"),
    data = model.frame(res$full$model)
  )

  save_forest_svg(
    res$horizon_5y$model,
    paste0("Cox_", model_info$name, "_UTI_PREVAL_to_ANY_LIVER_INCIDENT_5y_forest_plot.svg"),
    data = model.frame(res$horizon_5y$model)
  )

  save_forest_svg(
    res$horizon_10y$model,
    paste0("Cox_", model_info$name, "_UTI_PREVAL_to_ANY_LIVER_INCIDENT_10y_forest_plot.svg"),
    data = model.frame(res$horizon_10y$model)
  )

  # Extract result for main predictor and store
  results_all_models <- bind_rows(results_all_models, extract_main_predictor_results(res, model_info$name))
}

# Any model with PH violation in all follow-up horizons?
results_all_models %>%
  group_by(model) %>%
  summarise(all_violated = all(PH_meeting == "No")) %>%
  filter(all_violated == TRUE)
#   model                   all_violated
# 1 Sens_4_Current_drinkers TRUE          ---> We will run a stratified analysis with BMI of Sens_4_Current_drinkers to address this

# Export results for primary predictor across all models
results_all_models %>% write_csv("Cox_UTI_PREVAL_to_ANY_LIVER_INCIDENT.csv")

#-----------------------------
# 8.3 Create a BMI-stratified analysis for Sens 4 (Current drinkers) due to PH violation
#-----------------------------

# Modify df_drinkers to remove missing data on covariates
df_drinkers <- df_drinkers %>% drop_na(all_of(covs_comorb))
dim(df_drinkers) # 383088    119

# Check the distribution of BMI categories 
df_drinkers %>%
  drop_na(BMI_category) %>%
  count(UTI_PREVAL, BMI_category, ANY_LIVER_INCIDENT) %>%
  pivot_wider(names_from = ANY_LIVER_INCIDENT,  values_from = n, values_fill = 0) %>%
  rename(`No incident LD` = `0`, `Incident LD` = `1`) %>%
  mutate(Total = `No incident LD` + `Incident LD`) %>%
  mutate(BMI_category = fct_relevel(BMI_category, "<18.5", "18.5–<25", "25–<30", "≥30")) %>%
  arrange(UTI_PREVAL, BMI_category)
#   UTI_PREVAL BMI_category `No incident LD` `Incident LD`  Total
# 1 0          <18.5                    1625            45   1670
# 2 0          18.5–<25               125067          2346 127413
# 3 0          25–<30                 160509          4951 165460
# 4 0          ≥30                     79773          4664  84437
# 5 1          <18.5                      23             0     23 -> category dropped due to sparse data
# 6 1          18.5–<25                 1166            54   1220
# 7 1          25–<30                   1711            79   1790
# 8 1          ≥30                       972           103   1075

# Drop level <18.5, sparse data there
df_drinkers <- df_drinkers %>% filter(BMI_category != "<18.5") %>%  mutate(BMI_category = factor(BMI_category, levels = c("18.5–<25", "25–<30", "≥30")))
dim(df_drinkers) # 381395    119

# Run Cox with Sensitivity 4 scheme within each BMI category
results_BMI_strat <- data.frame()

for (category in levels(df_drinkers$BMI_category)) {
  res <- run_cox_model(
    df = df_drinkers %>% filter(BMI_category == category),
    outcome = "ANY_LIVER_INCIDENT",
    outcome_date = "ANY_LIVER_INCIDENT_DATE",
    primary_predictor = "UTI_PREVAL",
    covariates = setdiff(covs_comorb, "BMI") # remove BMI from covariates since we are stratifying by it
  )
  
  # Save individual model result
  res$model <- NULL # drop model object before saving (keeps only summaries)
  saveRDS(res, paste0("Cox_BMI_stratified_Sens_4_Current_drinkers_UTI_PREVAL_to_ANY_LIVER_INCIDENT_", make.names(category), ".rds"))

  # Extract result for main predictor and store
  results_BMI_strat <- bind_rows(results_BMI_strat, extract_main_predictor_results(res, paste0("Sens_4_Current_drinkers_BMI_", category)))

  # Save forest plot ONLY for 5-year horizon
  svg(paste0("Cox_BMI_stratified_Sens_4_Current_drinkers_UTI_PREVAL_to_ANY_LIVER_INCIDENT_", make.names(category), "_5y_forest_plot.svg"), width = 12, height = 8)
  survminer::ggforest(res$horizon_5y$model, data = model.frame(res$horizon_5y$model))
  dev.off()
}

# Export results for primary predictor across all models
results_BMI_strat %>% write_csv("Cox_BMI_stratified_Sens_4_Current_drinkers_UTI_PREVAL_to_ANY_LIVER_INCIDENT.csv")
```
