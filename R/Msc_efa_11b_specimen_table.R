##########################################################################
# The block is preparing the Pareto and classification columns for later csv files
##########################################################################
if (!exists("shape_perf")) {
} else {
PARETO_COLS <- c("goldberg_rank","strategy","limiting_objective", "distance_to_front","on_front","pareto_rank_ratio")
for (p in c("output/results/pareto_all_specimens.csv", "output/results/wing_classification.csv",
"output/results/Supplemental_performance/wing_classification.csv","output/results/reference_specimens_for_curvature.csv")) {
shape_perf <- join_missing(shape_perf, p, PARETO_COLS)
}
##########################################################################
# The block is creating of the diet combination for the EFA table
##########################################################################
if (!"Diet_combo" %in% names(shape_perf) &&
all(c("Diet.1","Diet.2") %in% names(shape_perf))) {
d1 <- trimws(as.character(shape_perf$Diet.1))
d2 <- trimws(as.character(shape_perf$Diet.2))
shape_perf$Diet_combo <- mapply(function(a, b) {
if (is.na(a) || a == "") return(NA_character_)
if (is.na(b) || b == "" || b == a) a else paste(sort(c(a, b)), collapse = "+")
}, d1, d2, USE.NAMES = FALSE)
cat(sprintf("    derived: Diet_combo (%d levels)\n",length(unique(na.omit(shape_perf$Diet_combo)))))
}
still_missing <- setdiff(c("family","goldberg_rank","strategy",  "limiting_objective","distance_to_front","Diet_combo"),
names(shape_perf))
if (length(still_missing) > 0) {
cat("  Still absent:", paste(still_missing, collapse = ", "), "\n")
if (any(c("strategy","limiting_objective","distance_to_front") %in% still_missing)) {
cat("    -> run Msc_wing_classification_v1.R\n")
}
if ("goldberg_rank" %in% still_missing) cat("    -> run Msc_pareto_front_v3.R\n")
}
##########################################################################
# The block is lining buidling the main EFA csv table
##########################################################################
ref_cols <- c("species","aspect_ratio","r2_hat","von_mises_stress","trailing_edge","shape_complexity",
"Order","clade","family","Environment","Diet.1","Diet.2","Midpoint","Period.Name",
"Depositional.settings.paleoenvironment","pareto_rank_ratio","goldberg_rank","on_front","strategy",
"limiting_objective","distance_to_front","Diet_combo","Time_Bin")
base_cols <- ref_cols[ref_cols %in% names(shape_perf)]
dropped <- setdiff(ref_cols, base_cols)
if (length(dropped) > 0) {
cat("  Columns not available:", paste(dropped, collapse = ", "), "\n")
}
sp_table <- shape_perf[, base_cols, drop = FALSE]
sp_table$shapePC1 <- round(shape_perf$shapePC1, 4)
sp_table$shapePC2 <- round(shape_perf$shapePC2, 4)
sp_table$PC1 <- pc1_vals
sp_table$PC2 <- pc2_vals
first_cols <- c("species","aspect_ratio","r2_hat", "von_mises_stress","trailing_edge","shape_complexity",
"Order","clade","family","Environment","Diet.1","Diet.2", "Midpoint","Period.Name",
"Depositional.settings.paleoenvironment","shapePC1","shapePC2","PC1","PC2")
first_cols <- first_cols[first_cols %in% names(sp_table)]
sp_table <- sp_table[, c(first_cols, setdiff(names(sp_table), first_cols))]

#################################################################################
# The block is lining exporting the csv tables, writing the loading along side
#################################################################################
out_path <- file.path(res_dir, "specimen_biomechanical_pca_table.csv")
write.csv(sp_table, out_path, row.names = FALSE)
cat(sprintf("  %d specimens x %d columns -> %s\n",
nrow(sp_table), ncol(sp_table), out_path))
if (exists("pca_performance")) {
write.csv(data.frame(variable = rownames(pca_performance$rotation), pca_performance$rotation[, 1:2]),"output/pca_loadings.csv", row.names = FALSE)
cat("  pca_loadings.csv refreshed from the current PCA.\n")
}
}
