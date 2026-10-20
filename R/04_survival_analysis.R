# =============================================================================
# 04_survival_analysis.R -- time to walnut removal (Kaplan-Meier, log-rank)
#
# Each walnut is an observation; the "event" is removal. Walnuts still present
# at the last check are right-censored (all of them were in the 100% dose).
# Pairwise log-rank p-values are adjusted with Benjamini-Hochberg (the default
# of survminer::pairwise_survdiff).
# =============================================================================

km_fit <- survfit(Surv(hours, event) ~ dose, data = walnut_surv)
print(km_fit)

# Kaplan-Meier curves with a risk table and the global log-rank p-value
km_plot <- ggsurvplot(
  km_fit, data = walnut_surv,
  risk.table = TRUE, pval = TRUE, conf.int = FALSE,
  break.time.by = 24,                      # one tick per day, shared with the risk table
  xlab = "Hours since first check", ylab = "Proportion of walnuts remaining",
  legend.title = "Urine dose", legend.labs = DOSE_LABELS,
  palette = unname(DOSE_COLORS),
  ggtheme = theme_minimal(base_size = 13)
)
png(file.path(FIG_DIR, "04_kaplan_meier.png"), width = 8, height = 6.5,
    units = "in", res = 300)
print(km_plot)
dev.off()

# Global log-rank test
global_logrank <- survdiff(Surv(hours, event) ~ dose, data = walnut_surv)
print(global_logrank)
save_table(data.frame(
  test = "Global log-rank",
  chisq = global_logrank$chisq,
  df = length(global_logrank$n) - 1,
  p_value = pchisq(global_logrank$chisq, length(global_logrank$n) - 1, lower.tail = FALSE)
), "04_global_logrank")

# Pairwise log-rank tests (BH-adjusted)
pairwise <- pairwise_survdiff(Surv(hours, event) ~ dose, data = walnut_surv,
                              p.adjust.method = "BH")
print(pairwise)

pairwise_long <- as.data.frame(pairwise$p.value) %>%
  mutate(dose_1 = rownames(.)) %>%
  pivot_longer(-dose_1, names_to = "dose_2", values_to = "p_adjusted") %>%
  filter(!is.na(p_adjusted)) %>%
  mutate(dose_1 = sub("dose=", "", dose_1), dose_2 = sub("dose=", "", dose_2))
save_table(pairwise_long, "04_pairwise_logrank")

p_pairwise <- ggplot(pairwise_long, aes(dose_1, dose_2, fill = p_adjusted)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(p_adjusted < 0.001, "<0.001",
                               sprintf("%.3f", p_adjusted))), size = 5) +
  scale_fill_gradient(low = "#0072B2", high = "#F0E442", name = "Adjusted p") +
  labs(title = "Pairwise log-rank tests (BH-adjusted)", x = NULL, y = NULL)
save_fig(p_pairwise, "04_pairwise_logrank", width = 6.5, height = 5)
