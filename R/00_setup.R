# =============================================================================
# 00_setup.R -- packages, shared constants and helper functions
#
# Sourced by every other script (via run_all.R). Nothing here touches the data.
# =============================================================================

required_packages <- c(
  "readxl", "dplyr", "tidyr", "lubridate", "purrr",
  "ggplot2", "survival", "survminer"
)
missing <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing) > 0) {
  stop("Install missing packages first: install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(purrr)
  library(ggplot2)
  library(survival)
  library(survminer)
})

# ---- Paths -------------------------------------------------------------------
DATA_FILE   <- file.path("data", "BIOL247_Data.xlsx")   # main time-series experiment
PILOT_FILE  <- file.path("data", "pilot_counts.csv")    # preliminary assay
FIG_DIR     <- "figures"
OUT_DIR     <- "output"
dir.create(FIG_DIR, showWarnings = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE)

# ---- Design constants ----------------------------------------------------------
DOSE_LEVELS <- c("0", "33", "66", "100")                 # % coyote urine
DOSE_LABELS <- c("0%", "33%", "66%", "100%")

# One colour per dose (colour-blind-safe Okabe-Ito palette), one per block.
DOSE_COLORS  <- c("0%" = "#E69F00", "33%" = "#56B4E9",
                  "66%" = "#009E73", "100%" = "#D55E00")
BLOCK_COLORS <- c("Rep1" = "#0072B2", "Rep2" = "#F0E442",
                  "Rep3" = "#009E73", "Rep4" = "#CC79A7")

theme_set(theme_minimal(base_size = 14))

# ---- Helpers -----------------------------------------------------------------

# Save a ggplot to figures/ with consistent size.
save_fig <- function(plot, name, width = 8, height = 5) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
  invisible(plot)
}

# Write a data frame to output/ as CSV.
save_table <- function(df, name) {
  write.csv(df, file.path(OUT_DIR, paste0(name, ".csv")), row.names = FALSE)
  invisible(df)
}
