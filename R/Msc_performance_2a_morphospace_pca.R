####################################################################################################
# The lines are coding building the performance morphospace and mark the observed envelope of points
####################################################################################################

######################################################################
# The block is picking the measurements used to build the morphospace
######################################################################
CORE_VARS <- c("aspect_ratio", "r2_hat","von_mises_stress", "trailing_edge", "shape_complexity")
missing_vars <- setdiff(CORE_VARS, colnames(performance_data))
if (length(missing_vars) > 0) {
stop("Missing performance columns: ", paste(missing_vars, collapse = ", "),
"\n  Present: ", paste(intersect(
c("aspect_ratio","r2_hat","von_mises_stress",
"trailing_edge","shape_complexity"),
colnames(performance_data)), collapse = ", "),
"\n  Run Msc_performance_metrics_v3.R first.")
}
PERFORMANCE_VARS <- CORE_VARS

######################################################################
# The block is building build the measurement matrix used for the PCA
######################################################################
performance_matrix <- performance_data %>% dplyr::select(all_of(PERFORMANCE_VARS)) %>% as.matrix()
complete_rows <- complete.cases(performance_matrix)
performance_matrix_clean <- performance_matrix[complete_rows, ]
performance_data_clean   <- performance_data[complete_rows, ]

cat(paste0("Specimens with complete performance data: ",
nrow(performance_matrix_clean), " / ", nrow(performance_data), "\n\n"))

if (any(performance_matrix_clean <= 0, na.rm = TRUE)) {
bad <- colnames(performance_matrix_clean)[
apply(performance_matrix_clean <= 0, 2, any)]
stop("Non-positive values in: ", paste(bad, collapse = ", "),
"\n  A zero usually means a degenerate outline.")
}
performance_matrix_log <- log(performance_matrix_clean)

######################################################################
# The block is checking whether any measurements are almost duplicates
######################################################################
cm <- cor(performance_matrix_log, method = "spearman")
high <- which(abs(cm) > 0.9 & upper.tri(cm), arr.ind = TRUE)
if (nrow(high) > 0) {
cat("WARNING - pairs correlating above 0.9:\n")
for (i in seq_len(nrow(high))) {
cat(sprintf("  %s / %s : rho = %+.3f\n", rownames(cm)[high[i,1]],
colnames(cm)[high[i,2]], cm[high[i,1], high[i,2]]))
}
cat("Consider dropping one of each pair before interpreting the axes.\n\n")
}

#######################################################################
# The block is linessing runing the PCA and store the first two axes
#######################################################################
pca_performance <- prcomp(performance_matrix_log, scale. = TRUE, center = TRUE)
performance_data_clean$PC1 <- pca_performance$x[, 1]
performance_data_clean$PC2 <- pca_performance$x[, 2]

#########################################################################################################
# The block is organizing the table explaining how much variance each PC axis explains the running total
#########################################################################################################
pca_summary <- summary(pca_performance)
variance_table <- data.frame(
Axis           = paste0("PC", seq_len(ncol(pca_performance$x))),
Eigenvalue     = pca_performance$sdev^2,
Variance_pct   = pca_summary$importance["Proportion of Variance", ] * 100,
Cumulative_pct = pca_summary$importance["Cumulative Proportion", ] * 100
)

cat("Variance explained by each PC axis:\n")
variance_table_display <- variance_table
variance_table_display$Eigenvalue     <- format_fixed(variance_table$Eigenvalue, 4)
variance_table_display$Variance_pct   <- format_fixed(variance_table$Variance_pct, 2)
variance_table_display$Cumulative_pct <- format_fixed(variance_table$Cumulative_pct, 2)
print(variance_table_display, row.names = FALSE)

cat(sprintf("\nFirst three axes collectively account for %s%% of total variance.\n",
format_fixed(variance_table$Cumulative_pct[min(3, nrow(variance_table))], 1)))
########################################################################
# The block is linessing saving the table into the supplemental folder
########################################################################
dir.create("output/results/Supplemental_performance", showWarnings = FALSE, recursive = TRUE)
write.csv(variance_table, "output/results/Supplemental_performance/performance_pca_variance_table.csv",
row.names = FALSE)
