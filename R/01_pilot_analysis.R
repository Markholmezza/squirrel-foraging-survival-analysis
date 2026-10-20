# =============================================================================
# 01_pilot_analysis.R -- preliminary assay (protocol validation)
#
# Design: 12 stations, 5 walnuts each; 6 controls (no urine) and 6 with a jar of
# 100% coyote urine. Walnuts counted on Day 0, Day 2 and Day 7.
# Question: is an urine effect detectable in the field at all?
# Test: Welch two-sample t test on Day-2 counts (paper: p = 0.01197).
# =============================================================================

pilot <- read.csv(PILOT_FILE)

pilot_long <- pilot %>%
  pivot_longer(c(day0, day2, day7), names_to = "day", values_to = "walnuts") %>%
  mutate(day = factor(day, levels = c("day0", "day2", "day7"),
                      labels = c("Day 0", "Day 2", "Day 7")))

# Figure: group means over time
pilot_means <- pilot_long %>%
  group_by(group, day) %>%
  summarise(mean_walnuts = mean(walnuts), .groups = "drop")

p_pilot <- ggplot(pilot_means, aes(day, mean_walnuts, colour = group, group = group)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = c("No urine" = "#0072B2", "Urine" = "#D55E00")) +
  labs(title = "Pilot assay: mean walnuts remaining",
       x = NULL, y = "Mean walnuts remaining (of 5)", colour = NULL)
save_fig(p_pilot, "01_pilot_means")

# Welch t test on Day 2
pilot_t <- t.test(day2 ~ group, data = pilot)
print(pilot_t)

save_table(
  data.frame(test = "Welch two-sample t test, Day 2 counts",
             t = unname(pilot_t$statistic),
             df = unname(pilot_t$parameter),
             p_value = pilot_t$p.value),
  "01_pilot_t_test"
)
