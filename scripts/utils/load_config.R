################################################################################
# load_config.R — Source this at the top of any R analysis script:
#   source("scripts/utils/load_config.R")
#
# Provides:
#   cfg              — full config list from config/pipeline_config.yaml
#   PROJECT_DIR      — project root path
#   DATA_DIR         — data/ directory
#   RESULTS_DIR      — results_final/ directory
#   CONDITION_COLORS — named vector of hex colors for the 4 conditions
#   CONDITION_ORDER  — factor-level order for condition axis
################################################################################

if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml", repos = "https://cloud.r-project.org", quiet = TRUE)
}

# Locate config relative to project root (works when run via Rscript from project root
# or when called via source() from a script anywhere in the project tree)
.config_candidates <- c(
  "config/pipeline_config.yaml",
  file.path(dirname(dirname(dirname(sys.frame(1)$ofile))), "config/pipeline_config.yaml")
)
.config_path <- .config_candidates[file.exists(.config_candidates)][1]

if (is.na(.config_path)) {
  stop("Cannot locate config/pipeline_config.yaml. Run scripts from the project root directory.")
}

cfg <- yaml::read_yaml(.config_path)

PROJECT_DIR <- cfg$project$dir
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results_final")

# Condition color palette (Control/DLB/PDD/PD)
CONDITION_COLORS <- unlist(cfg$visualization$condition_colors)
CONDITION_ORDER  <- cfg$visualization$condition_order

# Seed
set.seed(cfg$project$seed)

cat(sprintf("[load_config] Config loaded from: %s\n", .config_path))
cat(sprintf("[load_config] Project dir: %s | Seed: %d\n",
            PROJECT_DIR, cfg$project$seed))
