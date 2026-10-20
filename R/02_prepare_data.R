# =============================================================================
# 02_prepare_data.R -- read the main experiment and build analysis datasets
#
# Input:  data/BIOL247_Data.xlsx
#         one row per station ("Rep1_0", "Rep1_33", ... "Rep4_100"),
#         one column per check, column names are the check date-times
#         (e.g. "2025_12_03_16:00"; a column without a time means 12:00).
# Output: station_counts -- long format, one row per station x check
#         walnut_surv    -- one row per walnut (time-to-removal, event/censored)
# =============================================================================

raw <- read_excel(DATA_FILE)
names(raw) <- gsub("_", " ", names(raw))

# ---- Parse the check times from the column headers ---------------------------
time_cols <- names(raw)[-1]
time_cols <- ifelse(grepl(":", time_cols), time_cols, paste0(time_cols, " 12:00"))
check_times <- parse_date_time(time_cols, orders = c("Y m d H M", "Y m d I:M p"))
stopifnot("Some check times could not be parsed" = !anyNA(check_times))
names(raw)[-1] <- format(check_times, "%Y-%m-%d %H:%M:%S")

# ---- Long format: one row per station x check --------------------------------
station_counts <- raw %>%
  pivot_longer(-`Experimental Group`, names_to = "datetime", values_to = "walnuts") %>%
  mutate(
    datetime  = ymd_hms(datetime),
    block     = sub("_.*", "", `Experimental Group`),            # Rep1 ... Rep4
    dose      = factor(sub(".*_", "", `Experimental Group`),
                       levels = DOSE_LEVELS, labels = DOSE_LABELS),
    station   = paste(block, dose, sep = " / "),
    hours     = as.numeric(difftime(datetime, min(datetime), units = "hours"))
  ) %>%
  select(station, block, dose, datetime, hours, walnuts) %>%
  arrange(station, hours)

# ---- Data checks (fail loudly instead of analysing bad data) -----------------
design <- station_counts %>% distinct(block, dose)
stopifnot(
  "Expected 4 blocks x 4 doses = 16 stations" = nrow(design) == 16,
  "Dose labels did not parse"                 = !anyNA(station_counts$dose),
  "Missing walnut counts"                     = !anyNA(station_counts$walnuts)
)
count_increase <- station_counts %>%
  group_by(station) %>%
  summarise(any_increase = any(diff(walnuts) > 0), .groups = "drop")
stopifnot("A walnut count increased between checks (data entry error?)" =
            !any(count_increase$any_increase))

# ---- One row per walnut, for time-to-event analysis ---------------------------
# A walnut "fails" at the first check where it is no longer present.
# Walnuts still present at the last check are right-censored at that time.
# Note: the true removal time lies between two checks (interval censoring);
# like the paper, we use the check time. See README "Limitations".
removed <- station_counts %>%
  group_by(station) %>%
  arrange(hours, .by_group = TRUE) %>%
  mutate(n_removed = lag(walnuts, default = first(walnuts)) - walnuts) %>%
  ungroup() %>%
  filter(n_removed > 0) %>%
  uncount(n_removed) %>%
  mutate(event = 1L)

last_check <- max(station_counts$hours)
censored <- station_counts %>%
  filter(hours == last_check, walnuts > 0) %>%
  uncount(walnuts, .remove = FALSE) %>%
  mutate(event = 0L)

walnut_surv <- bind_rows(removed, censored) %>%
  select(station, block, dose, hours, event)

# ---- Descriptive summary (matches the Results section of the paper) ----------
study_summary <- tibble(
  stations             = n_distinct(station_counts$station),
  checks_per_station   = n_distinct(station_counts$hours),
  station_observations = nrow(station_counts),
  walnuts_deployed     = station_counts %>% filter(hours == 0) %>% pull(walnuts) %>% sum(),
  walnuts_removed      = sum(walnut_surv$event == 1),
  walnuts_censored     = sum(walnut_surv$event == 0),
  monitoring_hours     = last_check
)
print(study_summary)
save_table(study_summary, "02_study_summary")

stopifnot("Every deployed walnut should appear exactly once" =
            study_summary$walnuts_deployed ==
            study_summary$walnuts_removed + study_summary$walnuts_censored)

save_table(walnut_surv %>% count(dose, event, name = "walnuts"), "02_events_by_dose")
