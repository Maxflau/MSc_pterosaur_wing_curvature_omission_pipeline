############################################################################
# The block is building time-bin summary (requires Time_Bin and Midpoint)
############################################################################
if (!"Time_Bin" %in% names(sp) || !"Midpoint" %in% names(sp)) {
cat("  Time_Bin / Midpoint not in shape_perf - some K analyses skipped.\n")
has_time <- FALSE
} else {
has_time <- TRUE
}
##################################################################################################################
# The block is computing Pareto rank_ratio per Time Bin  and summarize each time bin and test simple through-time
##################################################################################################################
if (has_time && "pareto_rank_ratio" %in% names(sp)) {
cat("\n  K1: Temporal richness ~ pareto_rank_ratio\n")
time_summary <- sp %>%
filter(!is.na(Time_Bin), !is.na(pareto_rank_ratio)) %>% group_by(Time_Bin) %>%
summarise( n_species = n(), n_clades = n_distinct(clade), Midpoint= mean(Midpoint, na.rm=TRUE),
mean_opt = mean(pareto_rank_ratio,na.rm=TRUE), sd_opt = sd(pareto_rank_ratio,na.rm=TRUE),
mean_war = mean(aspect_ratio,na.rm=TRUE), mean_svm = mean(von_mises_stress,  na.rm=TRUE),     
SoV_PC1 = var(shapePC1, na.rm=TRUE), SoV_PC2 = var(shapePC2, na.rm=TRUE),
disparity = SoV_PC1 + SoV_PC2, mean_shapePC1 = mean(shapePC1, na.rm=TRUE),
mean_shapePC2 = mean(shapePC2,na.rm=TRUE), .groups= "drop") %>%
arrange(Midpoint) %>%
#####################################################################################################
# The block is calculating of the centroid distance from overall mean (morphospace expansion metric)
#####################################################################################################
mutate( centroid_dist = sqrt((mean_shapePC1 - mean(sp$shapePC1, na.rm=TRUE))^2 +(mean_shapePC2 - mean(sp$shapePC2, na.rm=TRUE))^2
))
save_csv(time_summary, "EFA_K_time_bin_summary")
cat(sprintf("  %d time bins analysed.\n", nrow(time_summary)))
################################################################################################
# The block is conducting Spearman test between number of species and the mean Pareto optimality
################################################################################################
if (nrow(time_summary) >= 4) {
cor_rich_opt  <- cor.test(time_summary$n_species, time_summary$mean_opt,method="spearman", exact=FALSE)
cor_rich_war  <- cor.test(time_summary$n_species, time_summary$mean_war,method="spearman", exact=FALSE)
cor_rich_disp <- cor.test(time_summary$n_species, time_summary$disparity,method="spearman", exact=FALSE)
cor_clade_opt <- cor.test(time_summary$n_clades, time_summary$mean_opt,method="spearman", exact=FALSE)
k1_results <- data.frame(
Test  = c("n_species ~ pareto_rank_ratio", "n_species ~ mean_WAR", "n_species ~ disparity (Stress)", "n_clades  ~ pareto_rank_ratio"),
rho = round(c(cor_rich_opt$estimate, cor_rich_war$estimate,
cor_rich_disp$estimate, cor_clade_opt$estimate), 4), p_value = round(c(cor_rich_opt$p.value, cor_rich_war$p.value,
cor_rich_disp$p.value, cor_clade_opt$p.value), 5), method = "Spearman", n_bins = nrow(time_summary)
) %>%
mutate(sig=case_when(p_value<0.001~"***",p_value<0.01~"**", p_value<0.05~"*",  p_value<0.1~".", TRUE~""))
save_csv(k1_results, "EFA_K1_richness_optimality_temporal_spearman")
cat("  K1 results:\n"); print(k1_results)
}
################################################################################################
# The block is lining undertaking regression analyses between optimality and aerodynamic metrics
################################################################################################
if (nrow(time_summary) >= 4) {
lm_disp_opt <- lm(disparity ~ mean_opt, data=time_summary)
lm_disp_war <- lm(disparity ~ mean_war, data=time_summary)
lm_disp_all <- lm(disparity ~ mean_opt + mean_war + mean_svm, data=time_summary)
k2_results <- rbind( data.frame(Model="disparity ~ mean_opt", R2=round(summary(lm_disp_opt)$r.squared,4),
Adj_R2=round(summary(lm_disp_opt)$adj.r.squared,4), F=round(summary(lm_disp_opt)$fstatistic[1],3),
p_value=round(pf(summary(lm_disp_opt)$fstatistic[1],
summary(lm_disp_opt)$fstatistic[2],
summary(lm_disp_opt)$fstatistic[3],  lower.tail=FALSE),5)),
data.frame(Model="disparity ~ mean_WAR",
R2=round(summary(lm_disp_war)$r.squared,4),
Adj_R2=round(summary(lm_disp_war)$adj.r.squared,4),F=round(summary(lm_disp_war)$fstatistic[1],3),
p_value=round(pf(summary(lm_disp_war)$fstatistic[1],summary(lm_disp_war)$fstatistic[2],summary(lm_disp_war)$fstatistic[3],lower.tail=FALSE),5)),
data.frame(Model="disparity ~ opt + WAR + SVM",
R2=round(summary(lm_disp_all)$r.squared,4), Adj_R2=round(summary(lm_disp_all)$adj.r.squared,4),
F=round(summary(lm_disp_all)$fstatistic[1],3), p_value=round(pf(summary(lm_disp_all)$fstatistic[1],
summary(lm_disp_all)$fstatistic[2], summary(lm_disp_all)$fstatistic[3],lower.tail=FALSE),5))
) %>% mutate(sig=case_when(p_value<0.001~"***",p_value<0.01~"**", p_value<0.05~"*",  p_value<0.1~".",TRUE~""))
names(k2_results)[names(k2_results) == "R2"]     <- "R\u00b2"
names(k2_results)[names(k2_results) == "Adj_R2"] <- "Adj_R\u00b2"
save_csv(k2_results, "EFA_K2_disparity_regression_models")
cat("  K2 results:\n"); print(k2_results)
}
}
