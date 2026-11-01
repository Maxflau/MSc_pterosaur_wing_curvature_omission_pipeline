#######################################################################
# The block is preparing the data used to compare depositional settings
#######################################################################
library(ggplot2)
for (dd in c("output/plots", "output/plots_PDF")) {
dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("METRICS")) stop("Source Msc_stats_0_setup.R first.")

LABELS <- c(aspect_ratio = "Aspect ratio", r2_hat = "Second moment of area",
stress_root = "Root bending stress index",
trailing_edge = "Trailing edge", shape_complexity = "Shape complexity",
pareto_rank_ratio = "Pareto rank ratio")

ENV_VAR <- intersect(c("Depositional.settings.paleoenvironment", "Environment"),
colnames(performance_data_clean))[1]
if (is.na(ENV_VAR)) stop("No environmental column found.")
cat("Environmental factor:", ENV_VAR, "\n")

d <- performance_data_clean
d[[ENV_VAR]] <- trimws(as.character(d[[ENV_VAR]]))
d <- d[!is.na(d[[ENV_VAR]]) & d[[ENV_VAR]] != "", ]
keep_lv <- names(which(table(d[[ENV_VAR]]) >= MIN_N))
d <- d[d[[ENV_VAR]] %in% keep_lv, ]
d$setting <- d[[ENV_VAR]]
cat(sprintf("Settings with n >= %d: %d | specimens: %d\n\n", MIN_N,
length(keep_lv), nrow(d)))

############################################################
# The block is comparing the metrics by depositional setting
############################################################

