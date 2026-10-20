# =============================================================================
# 05_sensitivity.R -- post-hoc checks added after the course (NOT in the paper)
#
# The paper's Kaplan-Meier / log-rank analysis treats each walnut as independent.
# Walnuts at the same station share the same squirrels and the same jar, so they
# are clustered; ignoring that can make p-values too small. These models ask
# whether the dose effect survives when the clustering and the block design are
# taken into account.
#
#   A. Cox model with block as a stratum and station-clustered (robust) SEs
#   B. Cox model with a station-level frailty (random effect)
#   C. Linear removal rates: ANOVA with block as a fixed effect
#
# Remaining limitation: removal times are interval-censored (known only to lie
# between two checks). An interval-censored model (e.g. icenReg) would be the
# next step; see README.
# =============================================================================

# Dose as a numeric score (0, 1, 2, 3) for a trend test, and as a factor.
walnut_surv <- walnut_surv %>%
  mutate(dose_score = as.integer(dose) - 1L)

# ---- A. Stratified Cox, robust SEs clustered by station -------------------------
cox_cluster <- coxph(Surv(hours, event) ~ dose + strata(block),
                     cluster = station, data = walnut_surv)
print(summary(cox_cluster))

cox_trend <- coxph(Surv(hours, event) ~ dose_score + strata(block),
                   cluster = station, data = walnut_surv)
print(summary(cox_trend))

hr_table <- function(fit, label) {
  s <- summary(fit)
  data.frame(
    model    = label,
    term     = rownames(s$conf.int),
    hazard_ratio = s$conf.int[, "exp(coef)"],
    ci_lower = s$conf.int[, "lower .95"],
    ci_upper = s$conf.int[, "upper .95"],
    p_robust = s$coefficients[, ncol(s$coefficients)],
    row.names = NULL
  )
}

# ---- B. Shared frailty by station ------------------------------------------------
cox_frailty <- coxph(Surv(hours, event) ~ dose + strata(block) +
                       frailty(station, distribution = "gamma"),
                     data = walnut_surv)
print(cox_frailty)

# Each station is a single block x dose cell, so after strata(block) and dose
# there is very little information left to estimate the station variance; it
# is reported as estimated (it can collapse to ~0).
frailty_summary <- summary(cox_frailty)
frailty_terms <- grep("^dose", rownames(frailty_summary$conf.int), value = TRUE)
save_table(data.frame(
  model            = "Cox, strata(block), gamma frailty(station)",
  term             = frailty_terms,
  hazard_ratio     = frailty_summary$conf.int[frailty_terms, "exp(coef)"],
  ci_lower         = frailty_summary$conf.int[frailty_terms, "lower .95"],
  ci_upper         = frailty_summary$conf.int[frailty_terms, "upper .95"],
  p_value          = frailty_summary$coefficients[frailty_terms, "p"],
  frailty_variance = cox_frailty$history[[1]]$theta,
  row.names = NULL
), "05_frailty")

save_table(bind_rows(
  hr_table(cox_cluster, "Cox, strata(block), cluster(station)"),
  hr_table(cox_trend,   "Cox trend, strata(block), cluster(station)")
), "05_cox_hazard_ratios")

# Proportional-hazards check for model A
ph_check <- cox.zph(cox_cluster)
print(ph_check)
save_table(data.frame(term = rownames(ph_check$table), ph_check$table,
                      row.names = NULL), "05_cox_ph_check")

# ---- C. Removal rates with block in the model -----------------------------------
slope_anova_block <- aov(linear_slope ~ dose + block, data = station_fits)
print(summary(slope_anova_block))
save_table(data.frame(term = trimws(rownames(summary(slope_anova_block)[[1]])),
                      summary(slope_anova_block)[[1]], row.names = NULL,
                      check.names = FALSE),
           "05_slope_anova_with_block")
