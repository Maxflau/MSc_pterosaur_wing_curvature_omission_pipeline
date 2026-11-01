##############################################################################
# The block is tabling 2: loadings + structure correlations (r, r^2) per PC
##############################################################################
n_axes_report <- 3   # how many PCs to tabulate; change as needed

loadings_raw <- pca_performance$rotation[, seq_len(n_axes_report), drop = FALSE]
eigenvalues  <- pca_performance$sdev[seq_len(n_axes_report)]^2

structure_r  <- sweep(loadings_raw, 2, sqrt(eigenvalues), `*`)   # correlation (r)
structure_r2 <- structure_r^2                                    # r^2

loadings_table <- data.frame(Variable = rownames(loadings_raw))
for (i in seq_len(n_axes_report)) {
loadings_table[[paste0("Loading_PC", i)]] <- format_fixed(loadings_raw[, i], 3)
loadings_table[[paste0("r_PC", i)]]        <- format_fixed(structure_r[, i], 3)
loadings_table[[paste0("r2_PC", i)]]       <- format_fixed(structure_r2[, i], 3)
}

cat("\nLoadings and structure correlations (r, r^2) per PC:\n")
print(loadings_table, row.names = FALSE)
cat("\n")

write.csv(loadings_table, "output/results/Supplemental_performance/performance_pca_loadings_table.csv",
row.names = FALSE)

cat("Correlation among log-transformed performance variables:\n")
print(round(cor(performance_matrix_log), 2))
cat("\n")
###################################################################################
# The block is building the envelope that contains most observed specimens
###################################################################################
ENVELOPE_LEVEL <- 0.95

pc_plane <- cbind(performance_data_clean$PC1, performance_data_clean$PC2)
envelope_center <- colMeans(pc_plane)
envelope_cov    <- cov(pc_plane)
envelope_cut    <- qchisq(ENVELOPE_LEVEL, df = 2)

performance_data_clean$mahalanobis_d2 <-
mahalanobis(pc_plane, center = envelope_center, cov = envelope_cov)

performance_data_clean$morphometric_outlier <-
performance_data_clean$mahalanobis_d2 > envelope_cut

cat(sprintf("Empirical envelope: %.0f%% chi-square ellipse (cutoff %.2f)\n",
ENVELOPE_LEVEL * 100, envelope_cut))
cat(sprintf("Morphometric outliers: %d / %d (%.1f%%)\n\n",
sum(performance_data_clean$morphometric_outlier),
nrow(performance_data_clean),
100 * mean(performance_data_clean$morphometric_outlier)))

if (any(performance_data_clean$morphometric_outlier)) {
cat("Specimens outside the envelope (inspect silhouettes and CSV rows):\n")
print(performance_data_clean %>%
dplyr::filter(morphometric_outlier) %>%
dplyr::arrange(desc(mahalanobis_d2)) %>%
dplyr::select(species, dplyr::any_of("clade"),
mahalanobis_d2, PC1, PC2) %>%
head(15))
cat("\n")
}

#########################################################################
# The block is saving the ellipse shape and the inside/outside test
#########################################################################
envelope_ellipse <- local({
theta <- seq(0, 2 * pi, length.out = 200)
circle <- cbind(cos(theta), sin(theta)) * sqrt(envelope_cut)
ev <- eigen(envelope_cov)
pts <- circle %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
data.frame(PC1 = pts[, 1] + envelope_center[1],
PC2 = pts[, 2] + envelope_center[2])
})

inside_envelope <- function(pc1, pc2) {
d2 <- mahalanobis(cbind(pc1, pc2), center = envelope_center, cov = envelope_cov)
d2 <= envelope_cut
}

cat("Empirical envelope defined. inside_envelope() and envelope_ellipse available.\n\n")
