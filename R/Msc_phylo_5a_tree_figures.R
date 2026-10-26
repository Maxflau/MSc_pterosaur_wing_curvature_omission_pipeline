#########################################################################
# The block is making sure the plot folders and input table are ready
#########################################################################
for (dd in c("output/plots", "output/plots_PDF")) {
dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("performance_data_clean")) stop("performance_data_clean not found.")
if (!exists("PTERODACT")) {
PTERODACT <- c("Azhdarchoidea", "Euctenochasmatia", "Dsungaripteroidea",
"Ornithocheiromorpha", "Pteranodontia", "Basal pterodactyliform")
}

###################################################################################
# The block is following the same convention as Msc_phylo_v3_1.R's order_group 
# keeping Darwinoptera among the non-pterodactyliform 
###################################################################################
perf_pareto_df <- performance_data_clean
perf_pareto_df$panel_group <- ifelse(
perf_pareto_df$clade %in% PTERODACT & perf_pareto_df$clade != "Darwinoptera",
"Pterodactyliform", "non-pterodactyliform")

#####################################################
# The block is testing the real Pareto objective
#####################################################
METRICS_D <- c("aspect_ratio", "r2_hat", "von_mises_stress")
missing_d <- setdiff(METRICS_D, names(perf_pareto_df))
if (length(missing_d) > 0) {
stop("Missing column(s): ", paste(missing_d, collapse = ", "),
"\n  Run Msc_performance_metrics_v3.R first.")
}

d <- perf_pareto_df[is.finite(perf_pareto_df$aspect_ratio) & perf_pareto_df$aspect_ratio > 0 &
is.finite(perf_pareto_df$r2_hat) & perf_pareto_df$r2_hat > 0 &
is.finite(perf_pareto_df$von_mises_stress) & perf_pareto_df$von_mises_stress > 0 &
!is.na(perf_pareto_df$panel_group), ]
cat(sprintf("Specimens with all three metrics finite and positive: %d / %d\n",
nrow(d), nrow(perf_pareto_df)))
print(table(d$panel_group))
cat("\n")

GROUP_COLS <- c("non-pterodactyliform" = "#9575CD", "Pterodactyliform" = "#FFE413")

has_theoretical <- exists("theoretical_data") &&
all(c("self_intersecting") %in% names(theoretical_data))

if (!has_theoretical) {
cat("theoretical_data not found - envelope will show real specimens' own hull only.\n")
cat("  Source Msc_theoretical_spaces_v3.R first to include the theoretical shapes too.\n\n")
}
