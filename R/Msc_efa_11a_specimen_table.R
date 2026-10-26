####################################################################
# The block is creating of EFA specimen grid biomechanical table
####################################################################
if (!exists("shape_perf")) {
cat("  shape_perf not found - skipping.\n")
} else {
if (!exists("res_dir")) res_dir <- "output/results/Supplemental_EFA"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("  Specimens: %d\n", nrow(shape_perf)))
####################################################################
# The block is calculating PC1/PC2
####################################################################
PCA_VARS <- c("aspect_ratio","r2_hat","von_mises_stress","trailing_edge","shape_complexity")
missing_pca <- setdiff(PCA_VARS, names(shape_perf))
if (length(missing_pca) > 0) {
stop("Missing metrics: ", paste(missing_pca, collapse = ", "),"\n  Re-run Msc_performance_metrics_v3.R (with calculate_von_mises_stress sourced).")
}
#######################################
# The block is creating EFA PCA loading
#######################################
use_loadings <- FALSE
if (file.exists("output/pca_loadings.csv")) {
loadings_df <- read.csv("output/pca_loadings.csv", stringsAsFactors = FALSE); rownames(loadings_df) <- loadings_df$variable
if (all(PCA_VARS %in% rownames(loadings_df))) {
use_loadings <- TRUE
} else {
cat("  pca_loadings.csv predates the metric rewrite - ignored.\n")
}
}
if (use_loadings) {
loadings_mat <- as.matrix(loadings_df[PCA_VARS, c("PC1","PC2")])
ref_means <- sapply(PCA_VARS, function(v) mean(shape_perf[[v]], na.rm = TRUE))
ref_sds   <- sapply(PCA_VARS, function(v) sd(shape_perf[[v]],   na.rm = TRUE))
z_mat <- sweep(sweep(as.matrix(shape_perf[, PCA_VARS]), 2, ref_means, "-"), 2, ref_sds, "/")
pc_proj  <- z_mat %*% loadings_mat
pc1_vals <- round(pc_proj[, 1], 4)
pc2_vals <- round(pc_proj[, 2], 4)
cat("  PC1/PC2 from the rotation matrix.\n")
} else if (all(c("PC1","PC2") %in% names(shape_perf))) {
pc1_vals <- round(shape_perf$PC1, 4)
pc2_vals <- round(shape_perf$PC2, 4)
cat("  PC1/PC2 taken from shape_perf.\n")
} else {
stop("No usable loadings and no PC1/PC2 in shape_perf.")
}
join_missing <- function(target, path, cols, key = "species") {
if (!file.exists(path)) return(target)
src <- read.csv(path, stringsAsFactors = FALSE)
if (!key %in% names(src)) return(target)
src[[key]] <- trimws(as.character(src[[key]]))
want <- intersect(setdiff(cols, names(target)), names(src))
if (length(want) == 0) return(target)
m <- match(trimws(as.character(target[[key]])), src[[key]])
for (v in want) target[[v]] <- src[[v]][m]
cat(sprintf("    from %s: %s (%d matched)\n", basename(path),paste(want, collapse = ", "), sum(!is.na(m))))
target
}
cat("  Recovering missing columns:\n")
##############################################################################################
# The block is linking the EFA table with the abiotic and taxonomic data from the main CSV
##############################################################################################
if (file.exists("data/pteros_main_data.csv")) {
raw_db <- read.csv("data/pteros_main_data.csv", sep = ";", stringsAsFactors = FALSE)
sp_col <- grep("^SPECIES", names(raw_db), value = TRUE)[1]
if (!is.na(sp_col)) {
raw_db$species <- trimws(raw_db[[sp_col]])
m <- match(trimws(shape_perf$species), raw_db$species)
for (v in c("family")) {
tgt <- v
if (v %in% names(raw_db) && !tgt %in% names(shape_perf)) {
shape_perf[[tgt]] <- trimws(as.character(raw_db[[v]][m]))
cat(sprintf("    from pteros_main_data.csv: %s (%d matched)\n",
tgt, sum(!is.na(m))))
}
}
}
}
} 
