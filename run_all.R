# Run the whole analysis from the repository root:
#   Rscript run_all.R
# Figures go to figures/, tables to output/.

set.seed(2025)
source("R/00_setup.R")
source("R/01_pilot_analysis.R")
source("R/02_prepare_data.R")
source("R/03_removal_models.R")
source("R/04_survival_analysis.R")
source("R/05_sensitivity.R")

cat("\nDone. Figures in figures/, tables in output/.\n")
writeLines(capture.output(sessionInfo()), file.path("output", "sessionInfo.txt"))
