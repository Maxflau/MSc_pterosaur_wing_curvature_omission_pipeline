################################################################################################
# The block is checking dataset and helper function are available for setting up the figure
################################################################################################
if (!exists("performance_data_clean")) {
stop("performance_data_clean not found - source Msc_impossible_region_v3.R first.")
}
if (!exists("inside_envelope")) {
stop("inside_envelope() not found - source Msc_impossible_region_v3.R first.")
}
########################################################################################################################
# The block is preparing and ensure the presence of the Von Mises Stress column for the Biomechanical morphospace figure
########################################################################################################################
if (!"von_mises_stress" %in% colnames(performance_data_clean)) {
if (exists("calculate_von_mises_stress") && exists("outlines_list")) {
performance_data_clean$von_mises_stress <- vapply(seq_len(nrow(performance_data_clean)), function(i) {
s <- performance_data_clean$species[i]
o <- outlines_list[[s]]
if (is.null(o) || nrow(o) < 4) return(NA_real_)
mk <- if ("Mass_kg" %in% colnames(performance_data_clean))
performance_data_clean$Mass_kg[i] else NA_real_
ws <- if ("Wingspan_cm" %in% colnames(performance_data_clean))
performance_data_clean$Wingspan_cm[i] else NA_real_
tryCatch(calculate_von_mises_stress(o, mk, ws), error = function(e) NA_real_)
}, numeric(1))
cat(sprintf("computed from outlines: %d of %d specimens\n",  sum(is.finite(performance_data_clean$von_mises_stress)), nrow(performance_data_clean)))
}
if (!"von_mises_stress" %in% colnames(performance_data_clean) ||
all(is.na(performance_data_clean$von_mises_stress))) {
for (p in c("output/results/pareto_all_specimens.csv","output/results/performance_metrics_complete_CLADE.csv",
"output/results/Supplemental_performance/wing_classification.csv"))
{
if (!file.exists(p)) next
src <- read.csv(p, stringsAsFactors = FALSE)
if (!all(c("species", "von_mises_stress") %in% names(src))) next
m <- match(trimws(performance_data_clean$species), trimws(src$species))
performance_data_clean$von_mises_stress <- src$von_mises_stress[m]
cat(sprintf("  read from %s: %d matched\n", basename(p), sum(!is.na(m))))
break
}
}
if (!"von_mises_stress" %in% colnames(performance_data_clean) ||
sum(is.finite(performance_data_clean$von_mises_stress)) < 20) {
stop("von_mises_stress could not be built.\n","Source in this order, then re-run this file:\n",
" Msc_biomechanical_functions_v3.R  (calculate_von_mises_stress)\n"," Msc_outline_loading.R (outlines_list)\n"," Msc_performance_metrics_v3.R")
}
}
#############################################################################
# The lines are choosing which measure will colour the background surface
#############################################################################
if (!exists("BACKGROUND_VAR")) BACKGROUND_VAR <- "von_mises_stress"
stopifnot(is.character(BACKGROUND_VAR), length(BACKGROUND_VAR) == 1)
###########################################################################
# The block is preparing the variable which will be consider for the figure
###########################################################################
if (!BACKGROUND_VAR %in% colnames(performance_data_clean)) {
stop("BACKGROUND_VAR '", BACKGROUND_VAR, "' not in performance_data_clean. ",
"Available: ", paste(intersect(c("von_mises_stress","r2_hat","aspect_ratio", "trailing_edge","shape_complexity","pareto_rank_ratio"),
colnames(performance_data_clean)), collapse = ", "))
}
BACKGROUND_LABELS <- c( von_mises_stress  = "Von Mises stress\n(scale-invariant)",
r2_hat = "Second moment\nof area", aspect_ratio = "Aspect ratio",
trailing_edge = "Trailing edge", shape_complexity = "Shape complexity",
pareto_rank_ratio = "Pareto\nrank ratio")
surf_label <- if (BACKGROUND_VAR %in% names(BACKGROUND_LABELS)) {
BACKGROUND_LABELS[[BACKGROUND_VAR]]
} else BACKGROUND_VAR
cat(sprintf("\nSurface variable: %s | colour bar: %s\n", BACKGROUND_VAR, gsub("\n", " ", surf_label)))
