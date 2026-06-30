# Plot for the table with all models and follow-up periods
library(dplyr)
library(forcats)
library(ggplot2)
library(stringr)
library(readr)

# Read the data
df <- read.csv("Cox_UTI_PREVAL_to_ANY_LIVER_INCIDENT.csv",stringsAsFactors = FALSE)

# Define model and follow-up period mappings
model_map <- c(
  "Model_0_Crude"             = "Model 0 - Crude",
  "Model_1_Demographics"      = "Model 1 - Demographics",
  "Model_2_Sociodemographics" = "Model 2 - Socio-demographics",
  "Model_3_Lifestyle"         = "Model 3 - Lifestyle",
  "Model_4_Comorbidity"       = "Model 4 - Comorbidity",
  "Sens_1_Physical_activity"  = "Sensitivity 1 - Physical activity",
  "Sens_2_Biomarkers"         = "Sensitivity 2 - Biomarkers",
  "Sens_3_Alcohol_numeric"    = "Sensitivity 3 - Alcohol numeric",
  "Sens_4_Current_drinkers"   = "Sensitivity 4 - Current drinkers",
  "Sens_5_GP_registered"      = "Sensitivity 5 - GP registered"
)

followup_map <- c("Full" = "17.5-year (full)", "5-year" = "5-year", "10-year" = "10-year")

# Process the data for plotting
df_plot <- df %>%
  mutate(
    HR_CI = str_replace_all(HR_CI, "–", "-"),
    HR    = as.numeric(str_extract(HR_CI, "^[0-9.]+")),
    lower = as.numeric(str_match(HR_CI, "\\(([^-]+)-([^)]+)\\)")[, 2]),
    upper = as.numeric(str_match(HR_CI, "\\(([^-]+)-([^)]+)\\)")[, 3]),
    model_clean = recode(model, !!!model_map),
    followup = recode(follow_up, !!!followup_map)
  ) %>%
  group_by(model) %>%
  mutate(
    model_label = paste0(first(model_clean), " (n=", format(first(n_participants), big.mark = ","), ")")
  ) %>%
  ungroup()

model_order <- c("Model 0 - Crude", "Model 1 - Demographics", "Model 2 - Socio-demographics", "Model 3 - Lifestyle", "Model 4 - Comorbidity",
                "Sensitivity 1 - Physical activity", "Sensitivity 2 - Biomarkers", "Sensitivity 3 - Alcohol numeric", "Sensitivity 4 - Current drinkers", "Sensitivity 5 - GP registered")

df_plot <- df_plot %>%
  mutate(
    model_clean = factor(model_clean, levels = model_order),
    model_label = factor(
      model_label,
      levels = df_plot %>%
        distinct(model_clean, model_label) %>%
        arrange(model_clean) %>%
        pull(model_label)
    ),
    model_label = fct_rev(model_label),
    followup = factor(followup, levels = c("5-year", "10-year", "17.5-year (full)"))
  )

pd <- position_dodge(width = 0.65)

p <- ggplot(df_plot, aes(x = HR, y = model_label, colour = followup)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_errorbar( aes(xmin = lower, xmax = upper), position = pd, width = 0.15, linewidth = 0.7) +
  geom_point(position = pd, size = 2.8) +
  scale_color_manual( values = c( "17.5-year (full)" = "#1b9e77", "5-year" = "#d95f02", "10-year" = "#7570b3")) +
  scale_x_log10(limits = c(1, 4), breaks = c(1, 1.25, 1.5, 2, 2.5, 3, 4), labels = c("1", "1.25", "1.5", "2", "2.5", "3", "4")) +
  labs(x = "Hazard ratio for incident liver disease (95% CI)", y = "Prevalent UTI association models", colour = "Follow-up period") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = "top", axis.text.y = element_text(size = 10))

svg("Cox_UTI_PREVAL_to_ANY_LIVER_INCIDENT_plot_13Mar.svg", width = 8, height = 6)
p
dev.off()
