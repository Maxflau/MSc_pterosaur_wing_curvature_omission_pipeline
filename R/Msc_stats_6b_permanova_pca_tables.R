######################################################################
# The block is betadisperring - biomechanical (multivariate) space
######################################################################
bio_dist <- D_FULL
bio_groups <- intersect(GROUPS, colnames(MV_META))
skipped_groups <- setdiff(GROUPS, bio_groups)
if (length(skipped_groups) > 0) {
cat("NOTE: skipping betadisper (biomech) for groups not in MV_META:",
paste(skipped_groups, collapse = ", "), "\n")
}
bio_disp_list <- lapply(bio_groups, function(g) {
grp <- as.character(MV_META[[g]])
ok <- !is.na(grp) & grp != ""
if (sum(ok) < 10 || length(unique(grp[ok])) < 2) {
cat("SKIPPED betadisper (biomech) for", g, "- too few usable specimens or levels\n")
return(NULL)
}
Ds <- as.dist(as.matrix(bio_dist)[ok, ok])
tidy_betadisper(Ds, factor(grp[ok]), g)
})
bio_disp_df <- do.call(rbind, bio_disp_list)
######################################################################
# The block is showing a fixed number of decimals so every row lines up
######################################################################
bio_disp_display <- bio_disp_df
bio_disp_display$SumSq   <- format_fixed(bio_disp_df$SumSq, 4)
bio_disp_display$MeanSq  <- format_fixed(bio_disp_df$MeanSq, 4)
bio_disp_display$F       <- format_fixed(bio_disp_df$F, 3)
bio_disp_display$p_value <- format_fixed(bio_disp_df$p_value, 4)
write.csv(bio_disp_display, file = paste0(OUT, "08_betadisper_biomech.csv"), row.names = FALSE)
cat("08 - Biomechanical betadisper CSV written.\n")
##########################################################################
# The block is kruskal-wallissing tests, univariate, per metric x group
##########################################################################
KW_METRICS <- intersect(c("aspect_ratio", "r2_hat", "von_mises_stress",
"trailing_edge", "shape_complexity",
"reynolds", "PC1", "PC2"), colnames(df))
cat("\nKruskal-Wallis metrics:", paste(KW_METRICS, collapse = ", "), "\n\n")
kw_results <- data.frame()
for (g in GROUPS) {
for (v in KW_METRICS) {
ok <- !is.na(df[[v]]) & !is.na(df[[g]]) & df[[g]] != ""
if (sum(ok) < 10 || length(unique(df[[g]][ok])) < 2) next
kw <- kruskal.test(df[[v]][ok], as.factor(df[[g]][ok]))
kw_results <- rbind(kw_results, data.frame(
Group = g, Variable = v,
H = round(kw$statistic, 3), df_kw = kw$parameter,
p_value = round(kw$p.value, 4),
sig = ifelse(kw$p.value < 0.001, "***",
ifelse(kw$p.value < 0.01,  "**",
ifelse(kw$p.value < 0.05,  "*",
ifelse(kw$p.value < 0.1,   ".", ""))))))
}
}
cat("\nKruskal-Wallis omnibus results:\n")
##########################################################################
# The block is showing a fixed number of decimals so every row lines up
##########################################################################
kw_results_display <- kw_results
if (nrow(kw_results_display) > 0) {
kw_results_display$H       <- format_fixed(kw_results$H, 3)
kw_results_display$p_value <- format_fixed(kw_results$p_value, 4)
}
print(kw_results_display, row.names = FALSE)
write.csv(kw_results_display, file = paste0(OUT, "09_kruskal_wallis.csv"), row.names = FALSE)
cat("09 - Kruskal-Wallis CSV written.\n\n")
