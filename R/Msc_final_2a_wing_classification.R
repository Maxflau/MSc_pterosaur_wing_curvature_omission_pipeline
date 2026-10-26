#################################################################################################################################
# The lines are coding classifying the studied wing to assess which functional trade off hold back the optimality of per specimen
#################################################################################################################################

#####################################################################################################################
# The block is setting up the metric to cassify each taxon by which performance trade-off appears to hold it back
#####################################################################################################################
OBJECTIVES <- c(aspect_ratio = "max", aspect_ratio = "min",   r2_hat = "min", r2_hat = "max", von_mises_stress = "min", von_mises_stress = "max")
OBJ_MEANING <- c( aspect_ratio= "aerodynamic efficiency (lift-to-drag)", 
r2_hat = "distal area distribution (average lift in flapping)",
von_mises_stress = "root bending load (scale-invariant)")

############################################################################################
# The block is verifying  which performance objectives can be used for the ranking
############################################################################################
missing_obj <- setdiff(names(OBJECTIVES), colnames(performance_data_clean))
if (length(missing_obj) > 0) {
cat("Missing objective column(s):", paste(missing_obj, collapse = ", "), "\n")
cat("Present metrics:", paste(intersect( c("aspect_ratio", "r2_hat", "von_mises_stress", "trailing_edge", "shape_complexity"),
colnames(performance_data_clean)), collapse = ", "), "\n")
OBJECTIVES <- OBJECTIVES[!names(OBJECTIVES) %in% missing_obj]
cat("Dropped:", paste(missing_obj, collapse = ", "), "- running with", length(OBJECTIVES), "objectives.\n")
if (length(OBJECTIVES) < 2) {
stop("Fewer than two objectives available - a Pareto front needs at least ",
"two. Source the biomechanical functions script (calculate_von_mises_stress) ","and Msc_performance_metrics_v3.R first.")
}
cat("\n")
}
cat("Objectives in use:", paste(names(OBJECTIVES), OBJECTIVES, sep = " (", collapse = "), "), ")\n\n")

############################################################################
# The block is sorting out the theoretical shapes into broad outcome classes
############################################################################
if (exists("theoretical_data") &&
all(names(OBJECTIVES) %in% colnames(theoretical_data))) {
theoretical_data$geometry_class <- ifelse( theoretical_data$self_intersecting, "Geometrically impossible", "Possible")
obj_ranges <- lapply(names(OBJECTIVES), function(v) range(performance_data_clean[[v]], na.rm = TRUE))
names(obj_ranges) <- names(OBJECTIVES)
within_observed <- Reduce(`&`, lapply(names(OBJECTIVES), function(v) {
x <- theoretical_data[[v]]
is.finite(x) & x >= obj_ranges[[v]][1] & x <= obj_ranges[[v]][2]
}))
theoretical_data$realisation_class <- ifelse( theoretical_data$self_intersecting, "Geometrically impossible",
ifelse(within_observed, "Realised region", "Possible but unrealised"))
cat("THEORETICAL SHAPES:\n")
print(table(theoretical_data$realisation_class))
cat("\n'Possible but unrealised' is the non-circular result: shapes that could\n")
cat("exist and that no pterosaur occupied. The cause - mechanical, phylogenetic\n")
cat("or ecological - is NOT determined by this analysis.\n\n")
} else {
cat("theoretical_data unavailable or missing objectives - levels 1-2 skipped.\n\n")
}
