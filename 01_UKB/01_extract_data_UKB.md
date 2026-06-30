# Project Documentation & Pipeline: Liver Disease & UTI in UK Biobank

This document describes the R-based processing pipeline for the UK Biobank (UKB) Liver Disease & UTI Analysis project.

Objective: Extract ICD-10 diagnoses, derive prevalent/incident case definitions, and assemble clinical covariates for downstream time-to-event (survival) models.

Key inputs:
- Cleaned UKB dataset. File: `240112_ukb_data_clean.rda`. Content: UKB-derived dataset prepared in 2024 (prepared by William)
- Consent withdrawals (exclude from analyses). File: `w55469_20250818.csv`. Content: list of participants who have withdrawn consent (prepared by Fumi)

To identify and verify relevant UKB fields/column names:
- Local reference: `ukb_column_names_data_487044_28360.txt`. Use this to map field IDs and locate candidate variables in the workspace.
- Online reference: `https://biobank.ndph.ox.ac.uk/ukb/search.cgi`. Useful for searching UKB field descriptions, encodings, and data-coding schemes.

# Extract Data 

```R
# ==============================================================================
# SECTION 1. GLOBAL CONFIGURATION & LIBRARIES
# ==============================================================================
library(tidyverse) # data manipulation and visualization
library(report)  # automated reporting of descriptive statistics
options(width = 200, scipen = 999) # improve console readability

# Define file paths 
PATH_DATA_RDA   <- "240112_ukb_data_clean.rda"
PATH_WITHDRAWAL <- "w55469_20250818.csv"

# ==============================================================================
# SECTION 2. DATA LOADING & PRIMARY FILTERING
# ==============================================================================
# Load the pre-cleaned UKB dataset
load(PATH_DATA_RDA)
dim(ukb_data_clean) # 487044  28360

# Handle Withdrawn Consent: 
# UKB requires immediate removal of participants who withdraw consent
withdrawn_ids <- read.csv(PATH_WITHDRAWAL, header = FALSE) %>% pull(V1)
length(withdrawn_ids) # 566 withdrawn IDs

ukb_data <- ukb_data_clean %>%
  filter(!eid %in% withdrawn_ids) %>%
  mutate(row_id = row_number()) # Create internal ID for efficient joining

# Clean up memory, remove the original large object
rm(ukb_data_clean)

# Check dataset dimensions after filtering
dim(ukb_data) # 486748  28361 (296 participants removed)

# ==============================================================================
# SECTION 3. CENSORING & DATE VALIDATION
# ==============================================================================
# Rationale: We need to define the end of the study follow-up
# We identify the most recent date in the dataset, excluding obvious errors 

date_cols <- grep("^date_", names(ukb_data), value = TRUE) %>%
  setdiff(grep("date_of_attending_assessment_centre", ., value = TRUE))
length(date_cols) # 1830 date columns to check

# Identify the census date (excluding the known error date 2037-07-07)
ukb_data %>%
  summarise(across(all_of(date_cols), ~ max(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "field", values_to = "max_date") %>%
  filter(max_date < as.Date("2037-01-01")) %>% # Exclude erroneous dates
  arrange(desc(max_date)) %>%
  head(3) # View top 3 most recent valid dates

# STUDY_CENSUS_DATE is set to 2023-09-01 based on data availability
STUDY_CENSUS_DATE <- as.Date("2023-09-01")

# ==============================================================================
# SECTION 4: DATA EXTRACTION FUNCTIONS
# ==============================================================================

# Extract ICD10 Binary Indicators and First-Reported Dates
# data: The UKB dataframe
# codes: A vector of ICD10 prefixes (e.g., "K70")
# details: Searches for 'date_<code>' and 'source_of_report_of_<code>'
extract_icd10_summary <- function(data, codes) {
  output <- data %>% select(row_id)
  
  for (code in codes) {
    message("Status: Processing ICD10 variable [", code, "]")
    
    # Identify variables using regex patterns
    d_var <- grep(paste0("^date_", code), names(data), value = TRUE, ignore.case = TRUE)
    s_var <- grep(paste0("^source_of_report_of_", code), names(data), value = TRUE, ignore.case = TRUE)
    
    # Build logic: 1 if any source column has data, 0 otherwise
    tmp <- data %>%
      select(row_id, all_of(d_var), all_of(s_var)) %>%
      mutate(
        !!paste0(code, "_date") := if (length(d_var) > 0) .[[d_var[1]]] else as.Date(NA),
        !!code := as.integer(if_any(all_of(s_var), ~ !is.na(.x) & .x != "" & .x != "."))
      ) %>%
      select(row_id, all_of(paste0(code, "_date")), all_of(code)) %>%
      mutate(!!code := factor(!!sym(code), levels = c(0, 1)))
    
    output <- left_join(output, tmp, by = "row_id")
  }
  return(output)
}

# ==============================================================================
# SECTION 5: FEATURE EXTRACTION
# ==============================================================================

# --- 5.1 Primary ICD10 Variable Groups ---
codes_liver    <- paste0("K", 70:77)
codes_diabetes <- paste0("E", 10:14)
codes_other    <- c("N39", "I10", "N18", "B15", "B16", "B17", "B18", "B19")

icd10_extracted <- extract_icd10_summary(ukb_data, c(codes_liver, codes_diabetes, codes_other))

# Quick verification of extraction, K70 should have 2349 cases
icd10_extracted$K70 %>% table(useNA = "ifany") # 0: 484399, 1: 2349

# --- 5.2 Specific UTI (N39.0) Extraction ---
# Rationale: N390 (UTI, site not specified) requires searching through the 'Diagnoses - ICD10' (f41270) and 'Date of first diagnosis' (f41280) arrays
diag_slots <- grep("diagnoses_icd10_f41270", names(ukb_data), value = TRUE)

uti_long <- ukb_data %>%
  select(row_id, all_of(diag_slots)) %>%
  pivot_longer(-row_id, names_to = "slot", values_to = "icd_code") %>%
  filter(icd_code == "N390") %>%
  mutate(idx = str_extract(slot, "\\d+_\\d+$"),
         date_field = paste0("date_of_first_inpatient_diagnosis_icd10_f41280_", idx))

# Verify UTI cases count, should be 30611
uti_long %>% distinct(row_id) %>% nrow() # 30611, matches expectation

# Match diagnosis dates back from original columns
uti_dates <- uti_long %>%
  rowwise() %>%
  mutate(UTI_date = ukb_data[[date_field]][ukb_data$row_id == row_id]) %>%
  group_by(row_id) %>%
  summarise(UTI_date = min(UTI_date, na.rm = TRUE), .groups = "drop") %>%
  mutate(UTI = factor(1, levels = c(0, 1)))

# Safely joining and handling missing factor levels
uti_final <- data.frame(row_id = ukb_data$row_id) %>%
  left_join(uti_dates, by = "row_id") %>%
  mutate(UTI = fct_expand(UTI, "0")) %>% # add the '0' level first
  mutate(UTI = replace_na(UTI, "0")) %>% # fill NAs with the character '0'
  mutate(UTI = fct_relevel(UTI, "1", "0")) # ensure 1 is the reference level

# Check row_ids 
head(uti_final, 5)

# --- 5.3 Demographic & Laboratory Covariates ---
# Rationale: Mapping cryptic UKB field IDs to human-readable names for analysis
field_map <- c(
  gp_reg                = "gp_registration_records_f42038_0_0",
  age_recruitment       = "age_at_recruitment_f21022_0_0",
  date_attendance       = "date_of_attending_assessment_centre_f53_0_0",
  sex                   = "sex_f31_0_0",
  BMI                   = "body_mass_index_bmi_f21001_0_0",
  ethnicity             = "ethnic_background_f21000_0_0",
  household_income      = "average_total_household_income_before_tax_f738_0_0",
  assessment_centre     = "uk_biobank_assessment_centre_f54_0_0",
  alcohol_frequency     = "alcohol_intake_frequency_f1558_0_0",
  alcohol_g             = "alcohol_f100022_0_0",
  smoking_status        = "current_tobacco_smoking_f1239_0_0",
  physical_activity     = "summed_met_minutes_per_week_for_all_activity_f22040_0_0",
  creatinine            = "creatinine_f23478_0_0",
  c_reactive_protein    = "creactive_protein_f30710_0_0",
  hdl_cholesterol       = "hdl_cholesterol_f30760_0_0",
  triglycerides         = "triglycerides_f30870_0_0",
  date_lost_to_followup = "date_lost_to_followup_f191_0_0",
  date_of_death         = "date_of_death_f40000_0_0",
  age_at_death          = "age_at_death_f40007_0_0"
)

covariates <- ukb_data %>% select(row_id, all_of(field_map))

# --- 5.4 Death Records ---
death_info <- ukb_data %>%
  select(row_id, death_date = date_of_death_f40000_0_0, death_age = age_at_death_f40007_0_0) %>%
  mutate(death_status = factor(as.integer(!is.na(death_date)), levels = c(0, 1)))

death_info %>% summarise(total_deaths = sum(as.integer(as.character(death_status)))) # 42841 deaths

# ==============================================================================
# SECTION 6: TEMPORAL CLASSIFICATION (INCIDENT VS PREVALENT)
# ==============================================================================

# Merge all extracted data
df_final <- list(icd10_extracted, covariates, uti_final, death_info) %>%
  reduce(left_join, by = "row_id")

# Check structure
str(df_final)

# Rationale: 
# Prevalent = Diagnosis occurred ON or BEFORE the assessment centre visit
# Incident  = Diagnosis occurred AFTER the assessment centre visit
classify_temporal <- function(df, codes) {
  for (code in codes) {
    d_col <- paste0(code, "_date")
    if (!d_col %in% names(df)) next
    
    df <- df %>%
      mutate(
        !!paste0(code, "_PREVAL")  := factor(as.integer(!is.na(!!sym(d_col)) & !!sym(d_col) <= date_attendance), levels = c(0, 1)),
        !!paste0(code, "_INCIDENT") := factor(as.integer(!is.na(!!sym(d_col)) & !!sym(d_col) > date_attendance), levels = c(0, 1))
      )
  }
  return(df)
}

analysis_codes <- c(codes_liver, codes_diabetes, codes_other, "UTI")
df_final <- classify_temporal(df_final, analysis_codes)

# --- Aggregate Clinical Blocks ---
# Create group-level variables (e.g., "Any Liver Disease" or "Any Diabetes")
create_clinical_block <- function(df, name, code_range) {
  p_cols <- paste0(code_range, "_PREVAL")
  i_cols <- paste0(code_range, "_INCIDENT")
  
  df %>% mutate(
    !!paste0(name, "_PREVAL")  := factor(as.integer(if_any(all_of(p_cols), ~ .x == "1")), levels = c(0, 1)),
    !!paste0(name, "_INCIDENT") := factor(as.integer(if_any(all_of(i_cols), ~ .x == "1")), levels = c(0, 1))
  )
}

df_final <- df_final %>%
  create_clinical_block("ANY_LIVER", codes_liver) %>%
  create_clinical_block("ANY_DIABETES", codes_diabetes) %>%
  create_clinical_block("ANY_VIRAL_HEP", paste0("B", 15:19))

# ==============================================================================
# SECTION 7: FINAL VERIFICATIONS
# ==============================================================================

# Verify PREVAL distributions
p_cols <- grep("_PREVAL$", names(df_final), value = TRUE)

report::report(df_final %>% select(all_of(p_cols))) 
# The data contains 486748 observations of the following 25 variables:
#   - K70_PREVAL: 2 levels, namely 0 (n = 486030, 99.85%) and 1 (n = 718, 0.15%)
#   - K71_PREVAL: 2 levels, namely 0 (n = 486684, 99.99%) and 1 (n = 64, 0.01%)
#   - K72_PREVAL: 2 levels, namely 0 (n = 486597, 99.97%) and 1 (n = 151, 0.03%)
#   - K73_PREVAL: 2 levels, namely 0 (n = 486562, 99.96%) and 1 (n = 186, 0.04%)
#   - K74_PREVAL: 2 levels, namely 0 (n = 486047, 99.86%) and 1 (n = 701, 0.14%)
#   - K75_PREVAL: 2 levels, namely 0 (n = 485117, 99.66%) and 1 (n = 1631, 0.34%)
#   - K76_PREVAL: 2 levels, namely 0 (n = 484653, 99.57%) and 1 (n = 2095, 0.43%)
#   - K77_PREVAL: 2 levels, namely 0 (n = 486731, 100.00%) and 1 (n = 17, 3.49e-03%)
#   - E10_PREVAL: 2 levels, namely 0 (n = 484213, 99.48%) and 1 (n = 2535, 0.52%)
#   - E11_PREVAL: 2 levels, namely 0 (n = 473565, 97.29%) and 1 (n = 13183, 2.71%)
#   - E12_PREVAL: 2 levels, namely 0 (n = 486747, 100.00%) and 1 (n = 1, 2.05e-04%)
#   - E13_PREVAL: 2 levels, namely 0 (n = 486663, 99.98%) and 1 (n = 85, 0.02%)
#   - E14_PREVAL: 2 levels, namely 0 (n = 465586, 95.65%) and 1 (n = 21162, 4.35%)
#   - N39_PREVAL: 2 levels, namely 0 (n = 459283, 94.36%) and 1 (n = 27465, 5.64%)
#   - I10_PREVAL: 2 levels, namely 0 (n = 356250, 73.19%) and 1 (n = 130498, 26.81%)
#   - N18_PREVAL: 2 levels, namely 0 (n = 480898, 98.80%) and 1 (n = 5850, 1.20%)
#   - B15_PREVAL: 2 levels, namely 0 (n = 485692, 99.78%) and 1 (n = 1056, 0.22%)
#   - B16_PREVAL: 2 levels, namely 0 (n = 486426, 99.93%) and 1 (n = 322, 0.07%)
#   - B17_PREVAL: 2 levels, namely 0 (n = 486467, 99.94%) and 1 (n = 281, 0.06%)
#   - B18_PREVAL: 2 levels, namely 0 (n = 486420, 99.93%) and 1 (n = 328, 0.07%)
#   - B19_PREVAL: 2 levels, namely 0 (n = 485464, 99.74%) and 1 (n = 1284, 0.26%)
#   - UTI_PREVAL: 2 levels, namely 0 (n = 480342, 98.68%) and 1 (n = 6406, 1.32%)
#   - ANY_LIVER_PREVAL: 2 levels, namely 0 (n = 482010, 99.03%) and 1 (n = 4738, 0.97%)
#   - ANY_DIABETES_PREVAL: 2 levels, namely 0 (n = 461379, 94.79%) and 1 (n = 25369, 5.21%)
#   - ANY_VIRAL_HEP_PREVAL: 2 levels, namely 0 (n = 483911, 99.42%) and 1 (n = 2837, 0.58%)

# Verify INCIDENT distributions
i_cols <- grep("_INCIDENT", names(df_final), value = TRUE)

report::report(df_final %>% select(all_of(i_cols))) 
# The data contains 486748 observations of the following 25 variables:
#   - K70_INCIDENT: 2 levels, namely 0 (n = 485117, 99.66%) and 1 (n = 1631, 0.34%)
#   - K71_INCIDENT: 2 levels, namely 0 (n = 486588, 99.97%) and 1 (n = 160, 0.03%)
#   - K72_INCIDENT: 2 levels, namely 0 (n = 485549, 99.75%) and 1 (n = 1199, 0.25%)
#   - K73_INCIDENT: 2 levels, namely 0 (n = 486600, 99.97%) and 1 (n = 148, 0.03%)
#   - K74_INCIDENT: 2 levels, namely 0 (n = 484465, 99.53%) and 1 (n = 2283, 0.47%)
#   - K75_INCIDENT: 2 levels, namely 0 (n = 485071, 99.66%) and 1 (n = 1677, 0.34%)
#   - K76_INCIDENT: 2 levels, namely 0 (n = 472533, 97.08%) and 1 (n = 14215, 2.92%)
#   - K77_INCIDENT: 2 levels, namely 0 (n = 486723, 99.99%) and 1 (n = 25, 5.14e-03%)
#   - E10_INCIDENT: 2 levels, namely 0 (n = 484098, 99.46%) and 1 (n = 2650, 0.54%)
#   - E11_INCIDENT: 2 levels, namely 0 (n = 454601, 93.40%) and 1 (n = 32147, 6.60%)
#   - E12_INCIDENT: 2 levels, namely 0 (n = 486741, 100.00%) and 1 (n = 7, 1.44e-03%)
#   - E13_INCIDENT: 2 levels, namely 0 (n = 486241, 99.90%) and 1 (n = 507, 0.10%)
#   - E14_INCIDENT: 2 levels, namely 0 (n = 482497, 99.13%) and 1 (n = 4251, 0.87%)
#   - N39_INCIDENT: 2 levels, namely 0 (n = 455846, 93.65%) and 1 (n = 30902, 6.35%)
#   - I10_INCIDENT: 2 levels, namely 0 (n = 420292, 86.35%) and 1 (n = 66456, 13.65%)
#   - N18_INCIDENT: 2 levels, namely 0 (n = 462336, 94.98%) and 1 (n = 24412, 5.02%)
#   - B15_INCIDENT: 2 levels, namely 0 (n = 486660, 99.98%) and 1 (n = 88, 0.02%)
#   - B16_INCIDENT: 2 levels, namely 0 (n = 486659, 99.98%) and 1 (n = 89, 0.02%)
#   - B17_INCIDENT: 2 levels, namely 0 (n = 486566, 99.96%) and 1 (n = 182, 0.04%)
#   - B18_INCIDENT: 2 levels, namely 0 (n = 486218, 99.89%) and 1 (n = 530, 0.11%)
#   - B19_INCIDENT: 2 levels, namely 0 (n = 486696, 99.99%) and 1 (n = 52, 0.01%)
#   - UTI_INCIDENT: 2 levels, namely 0 (n = 462543, 95.03%) and 1 (n = 24205, 4.97%)
#   - ANY_LIVER_INCIDENT: 2 levels, namely 0 (n = 469315, 96.42%) and 1 (n = 17433, 3.58%)
#   - ANY_DIABETES_INCIDENT: 2 levels, namely 0 (n = 451154, 92.69%) and 1 (n = 35594, 7.31%)
#   - ANY_VIRAL_HEP_INCIDENT: 2 levels, namely 0 (n = 485880, 99.82%) and 1 (n = 868, 0.18%)

# ==============================================================================
# SECTION 8: EXPORT
# ==============================================================================

# Save the final dataset for downstream analysis
saveRDS(df_final, file = "ukb_analysis_liver_uti_inc_prev_23Jan2026.rds")

# ------------------------------------------------------------------------------
# NOTE ON OVERLAPPING STATUS:
# It is possible for an individual to be flagged for both PREVALENT and INCIDENT status within the same clinical block 
#
#   A. AS COVARIATES (Adjusting for baseline health):
#      Use the "_PREVAL" flags. This accounts for the participant's medical 
#      history at the time of recruitment
#
#   B. AS OUTCOMES (Incident Analysis):
#      When analyzing a specific incident outcome (e.g., LIVER_ANY_INCIDENT), 
#      individuals flagged as LIVER_ANY_PREVAL should be EXCLUDED from the 
#      risk set to ensure we are only modeling true incident events
# ------------------------------------------------------------------------------
```