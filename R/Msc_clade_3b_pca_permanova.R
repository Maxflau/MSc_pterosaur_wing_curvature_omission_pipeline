###############################################################
# The block is conducting PERMANOVA tests along PC1 and PC2 axes
###############################################################
pc_dist <- dist(df[, c("PC1","PC2")], method = "euclidean")
set.seed(42)
perm_rows <- lapply(GROUPS, function(g) {
res <- adonis2(as.formula(paste("pc_dist ~", g)),data=df, permutations=9999, method="euclidean")
tbl <- as.data.frame(res); tbl$Term <- rownames(tbl); tbl$Factor <- g
tbl[!is.na(tbl$F), ]
})
perm_df <- do.call(rbind, perm_rows) %>%
select(Factor, Term, Df, SumOfSqs, R2, F, `Pr(>F)`) %>%
rename(p_value=`Pr(>F)`) %>%
mutate(sig = case_when(p_value<0.001~"***", p_value<0.01~"**", p_value<0.05~"*",   p_value<0.1~".", TRUE~""))
cat("\n=== PERMANOVA (PC1+PC2) ===\n")
#######################################################################
# The block is showing a fixed number of decimals so every row lines up
#######################################################################
perm_df_display <- perm_df
perm_df_display$SumOfSqs <- format_fixed(perm_df$SumOfSqs, 4)
perm_df_display$R2 <- format_fixed(perm_df$R2, 4)
perm_df_display$F <- format_fixed(perm_df$F, 3)
perm_df_display$p_value  <- format_fixed(perm_df$p_value, 4)
names(perm_df_display)[names(perm_df_display) == "R2"] <- "R²"
print(perm_df_display)
write.csv(perm_df_display, "output/results/Supplemental_performance/04_permanova_pc1pc2.csv", row.names = FALSE)

###############################################################################################
# The block is running pairwise PERMANOVA with Bonferroni correction on each grouping factor
#############################################################################################
set.seed(42)
for (g in c("clade","Diet_combined")) {
pw <- pairwise.adonis2(as.formula(paste("pc_dist ~", g)),
data=df, permutations=9999, p.adjust.m="bonferroni")
pw_rows <- lapply(names(pw)[names(pw) != "parent_call"], function(nm) {
tbl <- as.data.frame(pw[[nm]]); tbl$Comparison <- nm; tbl$Term <- rownames(tbl); tbl
})
pw_df <- do.call(rbind, pw_rows) %>%
filter(!Term %in% c("Residual","Total")) %>%
select(Comparison, Df, SumOfSqs, R2, F, `Pr(>F)`) %>%
rename(p_adj=`Pr(>F)`) %>%
mutate(sig = case_when(p_adj<0.001~"***", p_adj<0.01~"**", p_adj<0.05~"*",   p_adj<0.1~".", TRUE~""))
cat("\n--- Pairwise PERMANOVA:", g, "---\n")
##################################################################################
# The block is showing the exact statistical values in the supplemental csv table
##################################################################################
pw_df_display <- as.data.frame(pw_df)
pw_df_display$SumOfSqs <- format_fixed(pw_df$SumOfSqs, 4)
pw_df_display$R2       <- format_fixed(pw_df$R2, 4)
pw_df_display$F <- format_fixed(pw_df$F, 3)
pw_df_display$p_adj<- format_fixed(pw_df$p_adj, 4)
names(pw_df_display)[names(pw_df_display) == "R2"] <- "R²"
print(pw_df_display)
write.csv(pw_df_display, paste0("output/results/Supplemental_performance/05_pairwise_permanova_", g, ".csv"), row.names=FALSE)
}
##############################################################################################
# The block is evaluating the homogenity and dispersion across the clades and (a)biotic groups
##############################################################################################
disp_rows <- lapply(GROUPS, function(g) {
bd  <- betadisper(pc_dist, df[[g]])
pv  <- permutest(bd, permutations=9999)
tbl <- as.data.frame(pv$tab); tbl$Term <- rownames(tbl); tbl$Factor <- g
tbl[!is.na(tbl$F), ]
})
disp_df <- do.call(rbind, disp_rows) %>%
rename(p_value=`Pr(>F)`) %>%
select(Factor, Term, Df, `Sum Sq`, `Mean Sq`, F, p_value)
##################################################################################
# The block is showing the exact statistical values in the supplemental csv table
##################################################################################
disp_df_display <- disp_df
disp_df_display$`Sum Sq`  <- format_fixed(disp_df$`Sum Sq`, 4)
disp_df_display$`Mean Sq` <- format_fixed(disp_df$`Mean Sq`, 4)
disp_df_display$F         <- format_fixed(disp_df$F, 3)
disp_df_display$p_value   <- format_fixed(disp_df$p_value, 4)
write.csv(disp_df_display, "output/results/Supplemental_performance/06_betadisper_pc1pc2.csv", row.names=FALSE)

cat("\n=== All CSV files saved. Analysis complete. ===\n")
