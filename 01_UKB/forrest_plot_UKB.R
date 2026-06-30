# =============================================================================
# Forest plot — prior UTI and incident liver disease (UK Biobank)
# 10 models x 3 follow-up windows
# =============================================================================
library(ggplot2)

# ---- 1. Model structure (all 10 models) ------------------------------------
models <- tibble::tribble(
  ~order, ~model,                                  ~detail,                                           ~n,
  1L, "Model 0: Crude",                      "Unadjusted",                                       479546,
  2L, "Model 1: Demographics",               "Age at recruitment, sex",                          479546,
  3L, "Model 2: Socioeconomic",              "Model 1 + ethnicity, income",                      477278,
  4L, "Model 3: Lifestyle",                  "Model 2 + BMI, smoking, alcohol frequency",        474719,
  5L, "Model 4: Comorbidities",              "Model 3 + diabetes, hypertension, CKD",            474719,
  6L, "Sensitivity 1: Physical activity",    "Model 4 + physical activity",                      368072,
  7L, "Sensitivity 2: Biomarkers",           "Model 4 + creatinine, CRP, HDL-C, triglycerides",  226332,
  8L, "Sensitivity 3: Alcohol (continuous)", "Model 4 + alcohol as g/day",                       67497,
  9L, "Sensitivity 4: Current drinkers only","Model 4 excluding abstainers",                     383088,
  10L,"Sensitivity 5: GP-registered subset", "Model 4 restricted to GP records",                 216547
)
models$yc <- max(models$order) - models$order + 1   # Model 0 at the top

# ---- 2. Estimates (data from the table) -------------------------------
est <- tibble::tribble(
  ~order, ~horizon,     ~events,  ~hr,   ~lo,   ~hi,   ~p,
  1L,  "5 years",     3990,  2.25, 1.87, 2.72, "<0.0001",
  1L,  "10 years",   10275,  2.29, 2.04, 2.58, "<0.0001",
  1L,  "17.5 years", 16561,  2.12, 1.93, 2.34, "<0.0001",
  2L,  "5 years",     3990,  2.21, 1.83, 2.67, "<0.0001",
  2L,  "10 years",   10275,  2.24, 1.99, 2.52, "<0.0001",
  2L,  "17.5 years", 16561,  2.07, 1.88, 2.28, "<0.0001",
  3L,  "5 years",     3964,  2.08, 1.73, 2.52, "<0.0001",
  3L,  "10 years",   10222,  2.11, 1.88, 2.38, "<0.0001",
  3L,  "17.5 years", 16475,  1.96, 1.78, 2.16, "<0.0001",
  4L,  "5 years",     3904,  1.85, 1.52, 2.24, "<0.0001",
  4L,  "10 years",   10120,  1.90, 1.69, 2.15, "<0.0001",
  4L,  "17.5 years", 16326,  1.78, 1.61, 1.96, "<0.0001",
  5L,  "5 years",     3904,  1.66, 1.37, 2.01, "<0.0001",
  5L,  "10 years",   10120,  1.75, 1.55, 1.98, "<0.0001",
  5L,  "17.5 years", 16326,  1.66, 1.50, 1.84, "<0.0001",
  6L,  "5 years",     2888,  1.68, 1.34, 2.11, "<0.0001",
  6L,  "10 years",    7408,  1.74, 1.51, 2.01, "<0.0001",
  6L,  "17.5 years", 11949,  1.59, 1.41, 1.80, "<0.0001",
  7L,  "5 years",     1919,  2.05, 1.59, 2.63, "<0.0001",
  7L,  "10 years",    4824,  1.71, 1.44, 2.04, "<0.0001",
  7L,  "17.5 years",  7837,  1.61, 1.39, 1.86, "<0.0001",
  8L,  "5 years",      586,  2.19, 1.38, 3.47, "0.0009",
  8L,  "10 years",    1328,  1.82, 1.30, 2.57, "0.0006",
  8L,  "17.5 years",  1989,  1.64, 1.21, 2.21, "0.001",
  9L,  "5 years",     2864,  1.77, 1.39, 2.26, "<0.0001",
  9L,  "10 years",    7513,  1.73, 1.48, 2.02, "<0.0001",
  9L,  "17.5 years", 12242,  1.61, 1.42, 1.84, "<0.0001",
  10L, "5 years",     2376,  1.42, 1.09, 1.84, "0.009",
  10L, "10 years",    5520,  1.56, 1.32, 1.85, "<0.0001",
  10L, "17.5 years",  7986,  1.56, 1.36, 1.80, "<0.0001"
)

# ---- 3. Styling for the three follow-up windows ----------------------------
hlev  <- c("5 years", "10 years", "17.5 years")
cols  <- c("5 years" = "#9ecae1", "10 years" = "#4292c6", "17.5 years" = "#08519c")
dodge <- c("5 years" = 0.24, "10 years" = 0.00, "17.5 years" = -0.24)

est$horizon <- factor(est$horizon, levels = hlev)
est$yc <- models$yc[match(est$order, models$order)]
est$y  <- est$yc + dodge[as.character(est$horizon)]

bands <- models[models$order %% 2 == 1, ]          # alternating row shading

ylabs <- sprintf("%s\n%s\nn = %s",                 # model name / adjustment / N
                 models$model, models$detail,
                 formatC(models$n, big.mark = ",", format = "d"))

# ---- 4. Plot ---------------------------------------------------------------
p <- ggplot(est, aes(x = hr, y = y)) +
  geom_rect(data = bands, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = yc - 0.5, ymax = yc + 0.5),
            fill = "#f4f4f4") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#9a9a9a", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = horizon),
                 height = 0.16, linewidth = 0.9, show.legend = FALSE) +
  geom_point(aes(fill = horizon), shape = 21, size = 3.8,
             colour = "#333333", stroke = 0.4) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_fill_manual(values = cols, name = "Follow-up") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3), labels = c("1", "1.5", "2", "3")) +
  scale_y_continuous(breaks = models$yc, labels = ylabs,
                     expand = expansion(add = 0.7)) +
  coord_cartesian(xlim = c(0.95, 3.7), clip = "off") +
  labs(x = "Hazard ratio (95% CI)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid          = element_blank(),
    panel.grid.major.x  = element_line(colour = "#e8e8e8", linewidth = 0.4),
    axis.text.y         = element_text(size = 12,  colour = "#333333",
                                       lineheight = 0.95, hjust = 0),
    axis.ticks.y        = element_blank(),
    axis.title.y        = element_blank(),
    axis.text.x         = element_text(size = 15,  colour = "#222222"),
    axis.title.x        = element_text(size = 17,  colour = "#222222"),
    axis.line.x         = element_line(colour = "#555555", linewidth = 0.5),
    axis.ticks.x        = element_line(colour = "#555555", linewidth = 0.4),
    legend.position     = c(0.92, 0.97),
    legend.background   = element_rect(fill = adjustcolor("white", alpha.f = 0.95),
                                       colour = "#dddddd", linewidth = 0.3),
    legend.title        = element_text(face = "bold", size = 15),
    legend.text         = element_text(size = 14),
    legend.key          = element_blank(),
    legend.key.size     = unit(7, "mm"),
    plot.margin         = margin(8, 14, 8, 8)
  )

print(p)

# ---- 5. Save ---------------------------------------------------------------
ggsave("forest_plot_UTI_liver.svg", p, width = 10, height = 9, units = "in", device = function(filename, ...) svglite::svglite(filename, ...))  