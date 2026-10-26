#######################################################################################################
# The block is using the dplyr versions of common helper to create and manage the data with the tables
#######################################################################################################
select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
############################################################################################
# The lines are evaluating the load the PCA table and clean the names of the abiotic factors
############################################################################################
df <- read.csv("output/results/Supplemental_performance/performance_metrics_complete_CLADE.csv", stringsAsFactors = FALSE) %>%
rename(Depositional = Depositional.settings.paleoenvironment,Palaeoenvironment = Environment, Diet_primary = Diet.1, Diet_secondary = Diet.2)

###############################################################################
# The block is building the combined diet label if it is missing from the table
###############################################################################
if (!"Diet_combined" %in% names(df) && all(c("Diet_primary", "Diet_secondary") %in% names(df))) {
d1 <- trimws(as.character(df$Diet_primary)); d2 <- trimws(as.character(df$Diet_secondary))
df$Diet_combined <- ifelse(is.na(d2) | d2 == "" | d2 == d1, d1, paste(d1, d2, sep = " + "))
}
########################################################################
# The lines are choosing the biomechanical and group columns used below
########################################################################
BIO_VARS <- c("aspect_ratio","r2_hat","von_mises_stress","trailing_edge","shape_complexity")
GROUPS   <- c("clade","Depositional","Palaeoenvironment","Diet_primary","Diet_secondary","Diet_combined")

########################################################################
# The block is calculating the variance explained by the performance PCA
########################################################################
if (exists("pca_performance")) {
imp <- summary(pca_performance)$importance
var_exp <- c(PC1 = round(imp[2, 1] * 100, 2), PC2 = round(imp[2, 2] * 100, 2))
cum_var <- c(PC1 = var_exp[["PC1"]], PC2 = round(sum(imp[2, 1:2]) * 100, 2))
} else {
stop("pca_performance not found - source Msc_impossible_regions_v3.R.")
}
var_tbl <- data.frame(PC = names(var_exp), Var_pct = var_exp, Cumul_pct = cum_var)
cat("=== VARIANCE EXPLAINED (PC1 + PC2) ===\n"); print(var_tbl)
write.csv(var_tbl, "output/results/Supplemental_performance/01_pca_variance_explained.csv", row.names = FALSE)
##################################################################################################################
# The block is calculating how strongly each variable lines up with the PCs axes and Scale the selected measurements
##################################################################################################################
bio_scaled_mat <- scale(df[, BIO_VARS])

loadings_computed <- data.frame(
variable = BIO_VARS, PC1      = apply(bio_scaled_mat, 2, function(x) cor(x, df$PC1, use = "complete.obs")),
PC2      = apply(bio_scaled_mat, 2, function(x) cor(x, df$PC2, use = "complete.obs")))
cat("\n=== RECALCULATED LOADINGS (correlation with PC scores) ===\n")
print(loadings_computed)
write.csv(loadings_computed, "output/results/Supplemental_performance/00_pca_loadings_recalculated.csv", row.names = FALSE)
load_df <- loadings_computed
##########################################################
# The block is creating of the loading biplot(PC1 vs PC2)
##########################################################
arrow_scale <- max(abs(df[, c("PC1","PC2")])) * 0.6
p_biplot <- ggplot() +geom_point(data = df, aes(PC1, PC2, colour = clade), size = 1.5, alpha = 0.7) +
geom_segment(data = load_df, aes(x=0, y=0, xend=PC1*arrow_scale, yend=PC2*arrow_scale),
arrow = arrow(length = unit(0.2,"cm")), colour="black", linewidth=0.6) + geom_label_repel(data = load_df,
aes(x=PC1*arrow_scale, y=PC2*arrow_scale, label=variable), size=3, max.overlaps=20) +
labs(title  = "Biomechanical morphospace - PC1 vs PC2", x = paste0("PC1 (",var_exp[["PC1"]],"%)"), y = paste0("PC2 (",var_exp[["PC2"]],"%)"),
colour = "Clade") + theme_bw()
#####################################################################################################
# The lines are creating the supplementals folder here, since it is not created earlier in the driver
#####################################################################################################
dir.create("output/plots_PDF/supplementals", showWarnings = FALSE, recursive = TRUE)
ggsave("output/plots_PDF/supplementals/02_biplot_PC1_PC2.pdf", p_biplot, width=12, height=8)

###################################################################
# The block is summarising PC1 and PC2 within each group into a csv
###################################################################
summarise_pcs <- function(data, g) {
data %>% group_by(across(all_of(g))) %>%
summarise(N=n(), PC1_mean=round(mean(PC1),3), PC1_sd=round(sd(PC1),3),
PC2_mean=round(mean(PC2),3), PC2_sd=round(sd(PC2),3), .groups="drop")
}
for (g in GROUPS) {
tbl <- summarise_pcs(df, g)
cat("\n--- PC summary:", g, "---\n")
####################################################################################################################
# The block is showing the PCA summary output into the Rstudio console while being registered into the supplemental
####################################################################################################################
tbl_display <- tbl
tbl_display$PC1_mean <- format_fixed(tbl$PC1_mean, 3)
tbl_display$PC1_sd   <- format_fixed(tbl$PC1_sd, 3)
tbl_display$PC2_mean <- format_fixed(tbl$PC2_mean, 3)
tbl_display$PC2_sd   <- format_fixed(tbl$PC2_sd, 3)
print(tbl_display)
write.csv(tbl_display, paste0("output/results/Supplemental_performance/03_summary_pcs_", g, ".csv"), row.names = FALSE)
}
