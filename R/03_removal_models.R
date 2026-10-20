# =============================================================================
# 03_removal_models.R -- how fast are walnuts removed, and does dose matter?
#
# 1. Descriptive removal curves (by station and by dose, mean +/- SE).
# 2. For each station, fit a linear model and an exponential decay model of
#    walnuts remaining over time; compare them with R^2 and AIC.
# 3. One-way ANOVA (and Tukey HSD) of the station-level removal rates by dose,
#    using the better-fitting functional form (linear, as in the paper).
# =============================================================================

# ---- 1. Descriptive curves ------------------------------------------------------
p_stations <- ggplot(station_counts,
                     aes(hours, walnuts, colour = dose, group = station)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = DOSE_COLORS) +
  labs(title = "Walnut removal at each station",
       x = "Hours since first check", y = "Walnuts remaining", colour = "Urine dose")
save_fig(p_stations, "03_removal_by_station")

dose_means <- station_counts %>%
  group_by(dose, hours) %>%
  summarise(mean = mean(walnuts), se = sd(walnuts) / sqrt(n()), .groups = "drop")

p_means <- ggplot(dose_means, aes(hours, mean, colour = dose)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 1.5, alpha = 0.6) +
  scale_colour_manual(values = DOSE_COLORS) +
  labs(title = "Mean walnut removal by urine dose (+/- 1 SE)",
       x = "Hours since first check", y = "Mean walnuts remaining", colour = "Urine dose")
save_fig(p_means, "03_removal_by_dose")

# ---- 2. Linear vs exponential decay, per station ------------------------------
fit_station <- function(d) {
  lin <- lm(walnuts ~ hours, data = d)
  exp_fit <- tryCatch(
    nls(walnuts ~ a * exp(-k * hours), data = d,
        start = list(a = max(d$walnuts), k = 0.01)),
    error = function(e) NULL
  )
  pseudo_r2 <- function(fit) {
    1 - sum((d$walnuts - predict(fit))^2) / sum((d$walnuts - mean(d$walnuts))^2)
  }
  tibble(
    linear_slope   = unname(coef(lin)["hours"]),
    linear_r2      = summary(lin)$r.squared,
    linear_aic     = AIC(lin),
    exp_k          = if (is.null(exp_fit)) NA_real_ else unname(coef(exp_fit)["k"]),
    exp_pseudo_r2  = if (is.null(exp_fit)) NA_real_ else pseudo_r2(exp_fit),
    exp_aic        = if (is.null(exp_fit)) NA_real_ else AIC(exp_fit)
  )
}

station_fits <- station_counts %>%
  group_by(station, block, dose) %>%
  group_modify(~ fit_station(.x)) %>%
  ungroup()
save_table(station_fits, "03_station_model_fits")

if (any(is.na(station_fits$exp_aic))) {
  message("Exponential model failed to converge for: ",
          paste(station_fits$station[is.na(station_fits$exp_aic)], collapse = ", "))
}

# Table 1 of the paper: mean R^2 and AIC by dose
model_comparison <- station_fits %>%
  group_by(dose) %>%
  summarise(
    mean_linear_r2 = mean(linear_r2),
    mean_exp_r2    = mean(exp_pseudo_r2, na.rm = TRUE),
    mean_linear_aic = mean(linear_aic),
    mean_exp_aic    = mean(exp_aic, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(better_model = ifelse(mean_linear_aic < mean_exp_aic, "Linear", "Exponential"))
print(model_comparison)
save_table(model_comparison, "03_model_comparison_table1")

aic_long <- station_fits %>%
  select(station, dose, Linear = linear_aic, Exponential = exp_aic) %>%
  pivot_longer(c(Linear, Exponential), names_to = "model", values_to = "aic")

p_aic <- ggplot(aic_long, aes(model, aic, fill = model)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.7) +
  facet_wrap(~ dose, nrow = 1) +
  labs(title = "AIC by station: linear vs exponential decay (lower is better)",
       x = NULL, y = "AIC") +
  theme(legend.position = "none")
save_fig(p_aic, "03_aic_comparison", width = 10)

# ---- 3. Do removal rates differ by dose? --------------------------------------
# Unit of analysis = station (n = 16), so walnuts within a station are not
# treated as independent here. Block is ignored, as in the paper; the
# sensitivity analysis in 05_sensitivity.R adds it.
slope_anova <- aov(linear_slope ~ dose, data = station_fits)
print(summary(slope_anova))
slope_tukey <- TukeyHSD(slope_anova)
print(slope_tukey)

save_table(data.frame(term = trimws(rownames(summary(slope_anova)[[1]])),
                      summary(slope_anova)[[1]], row.names = NULL,
                      check.names = FALSE),
           "03_slope_anova")
save_table(data.frame(comparison = rownames(slope_tukey$dose), slope_tukey$dose,
                      row.names = NULL), "03_slope_tukey")

p_slopes <- ggplot(station_fits, aes(dose, linear_slope, colour = dose)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.1, size = 3) +
  scale_colour_manual(values = DOSE_COLORS) +
  labs(title = "Linear removal rate by urine dose (one point per station)",
       x = "Urine dose", y = "Slope (walnuts per hour)") +
  theme(legend.position = "none")
save_fig(p_slopes, "03_removal_rates")
