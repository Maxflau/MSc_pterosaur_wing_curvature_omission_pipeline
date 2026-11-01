##################################################################################################
# The block is checking that all the objects this script needs already exist before continuing
##################################################################################################
required_objects <- c("df", "pca_dist", "summarise_pcs", "test_dispersion",
"perm_clade", "perm_dep", "perm_env",
"perm_diet1", "perm_diet2", "perm_dietcombo", "perm_full",
"D_FULL", "MV_META")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
stop("Missing: ", paste(missing_objects, collapse = ", "),
"\n  Source Msc_stat_permanova_2.R AND Msc_stat_permanova.R",
"\n  (via Msc_stats_0_setup.R) before this script.")
}
######################################################################
# The block is defining the groups compared in every test below
######################################################################
GROUPS <- c("clade", "Depositional", "Palaeoenvironment",
"Diet_primary", "Diet_secondary", "Diet_combined")
perm_results_pca <- list(
clade              = perm_clade,
Depositional       = perm_dep,
Palaeoenvironment  = perm_env,
Diet_primary       = perm_diet1,
Diet_secondary     = perm_diet2,
Diet_combined      = perm_dietcombo)
##########################################################################################
# The block is turning one PERMANOVA (adonis2) result into a tidy, labelled table row set
##########################################################################################
tidy_adonis2 <- function(res, factor_label) {
tbl           <- as.data.frame(res)
tbl$Term      <- rownames(tbl)
tbl$Factor    <- factor_label
rownames(tbl) <- NULL
tbl <- tbl[, c("Factor", "Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
names(tbl)[names(tbl) == "Pr(>F)"] <- "p_value"
tbl
}
######################################################################################################
# The block is turning one betadisper (variance homogeneity) test into a tidy, labelled table row set
######################################################################################################
tidy_betadisper <- function(dist_mat, groups, factor_label) {
bd            <- betadisper(dist_mat, groups)
pv            <- permutest(bd, permutations = 9999)
tbl           <- as.data.frame(pv$tab)
tbl$Term      <- rownames(tbl)
tbl$Factor    <- factor_label
rownames(tbl) <- NULL
tbl <- tbl[, c("Factor", "Term", "Df", "Sum Sq", "Mean Sq", "F", "N.Perm", "Pr(>F)")]
names(tbl)[names(tbl) == "Pr(>F)"]  <- "p_value"
names(tbl)[names(tbl) == "Sum Sq"]  <- "SumSq"
names(tbl)[names(tbl) == "Mean Sq"] <- "MeanSq"
tbl
}
OUT <- "output/results/Supplemental_performance/"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("\n\n=== EXPORTING CSV FILES ===\n")
######################################################################
# The block is summarising the PCs by group (01)
######################################################################
for (g in GROUPS) {
pc_tbl <- summarise_pcs(df, g)
write.csv(pc_tbl, file = paste0(OUT, "01_summary_pcs_", g, ".csv"), row.names = FALSE)
}
cat("01 - PC summary CSVs written.\n")
######################################################################
# The block is reporting the individual PERMANOVA results in PCA space (03)
######################################################################
pca_perm_list <- lapply(names(perm_results_pca), function(g) {
tidy_adonis2(perm_results_pca[[g]], g)
})
pca_perm_df <- do.call(rbind, pca_perm_list)
######################################################################
# The block is showing a fixed number of decimals so every row lines up
######################################################################
pca_perm_display <- pca_perm_df
pca_perm_display$SumOfSqs <- format_fixed(pca_perm_df$SumOfSqs, 4)
pca_perm_display$R2       <- format_fixed(pca_perm_df$R2, 4)
pca_perm_display$F        <- format_fixed(pca_perm_df$F, 3)
pca_perm_display$p_value  <- format_fixed(pca_perm_df$p_value, 4)
names(pca_perm_display)[names(pca_perm_display) == "R2"] <- "R²"
write.csv(pca_perm_display,
file      = paste0(OUT, "03_permanova_pca.csv"),
row.names = FALSE)
cat("03 - PCA PERMANOVA CSV written.\n")
######################################################################
# The block is fitting the full marginal model in PCA space (04)
######################################################################
full_model_tbl <- tidy_adonis2(perm_full, "Full model (PCA)")
full_model_tbl$SumOfSqs <- format_fixed(full_model_tbl$SumOfSqs, 4)
full_model_tbl$R2       <- format_fixed(full_model_tbl$R2, 4)
full_model_tbl$F        <- format_fixed(full_model_tbl$F, 3)
full_model_tbl$p_value  <- format_fixed(full_model_tbl$p_value, 4)
names(full_model_tbl)[names(full_model_tbl) == "R2"] <- "R²"
write.csv(full_model_tbl,
file      = paste0(OUT, "04_permanova_pca_full_model.csv"),
row.names = FALSE)
cat("04 - PCA full model CSV written.\n")
######################################################################
# The block is running betadisper in PCA space (07)
######################################################################
set.seed(42)
pca_disp_list <- lapply(GROUPS, function(g) tidy_betadisper(pca_dist, df[[g]], g))
pca_disp_df <- do.call(rbind, pca_disp_list)
##########################################################################
# The block is showing a fixed number of decimals so every row lines up
##########################################################################
pca_disp_display <- pca_disp_df
pca_disp_display$SumSq   <- format_fixed(pca_disp_df$SumSq, 4)
pca_disp_display$MeanSq  <- format_fixed(pca_disp_df$MeanSq, 4)
pca_disp_display$F       <- format_fixed(pca_disp_df$F, 3)
pca_disp_display$p_value <- format_fixed(pca_disp_df$p_value, 4)
write.csv(pca_disp_display,
file      = paste0(OUT, "07_betadisper_pca.csv"),
row.names = FALSE)
cat("07 - PCA betadisper CSV written.\n")
