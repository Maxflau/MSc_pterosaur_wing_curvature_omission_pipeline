#################################################################################################################################
# The block is coding conducting taphonomic sensitivity analysis (via Krustal-Wallis) across the Biomechanical base morphospace
################################################################################################################################
run_taph_test <- function(pc_col, pc_label, tag) {
dat$.pc <- dat[[pc_col]]
group_means <- dat %>%
group_by(Pres_Quality) %>%
summarise(mean_pc = mean(.pc, na.rm = TRUE), n = n(), .groups = "drop")
#################################################
# The block is setting up the Krustal Wallis test
#################################################
kw <- kruskal.test(.pc ~ Pres_Quality, data = dat)
cat(sprintf("\n--- Taphonomic sensitivity: %s ---\n", pc_label))
print(kw)
y_lab_pos <- max(dat$.pc, na.rm = TRUE) * 1.05

#################################################################################################
# The block is constructing the structure of the taphonomic bar plot figures and Spearman test
################################################################################################
p1 <- ggplot(dat, aes(x = Pres_Quality, y = .pc, fill = Pres_Quality)) +
geom_boxplot(outlier.shape = NA, alpha = 0.2) + geom_jitter(width = 0.18, alpha = 0.6, size = 2, aes(color = Pres_Quality)) +
scale_fill_viridis(discrete = TRUE, end = 0.8) + scale_color_viridis(discrete = TRUE, end = 0.8, option = "D") +
ylab(sprintf("%s score", pc_label)) + xlab("Preservation quality (depositional environment)") +
ggtitle(sprintf("Taphonomic sensitivity: %s by preservation quality", pc_label)) +geom_text(data = group_means,
aes(x = Pres_Quality, y = y_lab_pos, label = paste0("mean=", round(mean_pc, 2), "\n(n=", n, ")")),
size = 3.5, inherit.aes = FALSE) + theme_minimal(base_size = 14) + guides(color = "none", fill = "none")
corsimple <- cor.test(dat$.pc, dat$Pres_Score, method = "spearman", exact = FALSE)
cat(sprintf("Raw correlation (%s): rho=%.3f, p=%.3g\n", pc_label, corsimple$estimate, corsimple$p.value))
p2 <- ggplot(dat, aes(x = Pres_Quality, y = .pc)) +
geom_jitter(width = 0.18, size = 2, alpha = 0.7, aes(color = Pres_Quality)) +
stat_smooth(method = "lm", aes(group = 1), size = 1.1, se = TRUE, color = "black") +
scale_color_viridis(discrete = TRUE, end = 0.8, option = "D") +
ylab(sprintf("%s score", pc_label)) + xlab("Preservation quality") +
annotate("text", x = 1, y = max(dat$.pc, na.rm = TRUE),
label = sprintf("rho=%.2f, p=%.2g", corsimple$estimate, corsimple$p.value),  size = 5, hjust = 0) +
theme_minimal(base_size = 13) + guides(color = "none")

#######################################################################################
# The block is building the linear model controlling for clade (partial correlation)
######################################################################################
mod <- lm(.pc ~ Pres_Score + clade, data = dat)
cat(sprintf("\n--- lm(%s ~ Pres_Score + clade) ---\n", pc_label))
print(summary(mod))

pc_clade_resid <- residuals(lm(.pc ~ clade, data = dat))
presscore_clade_resid <- residuals(lm(Pres_Score ~ clade, data = dat))
corresid <- cor.test(pc_clade_resid, presscore_clade_resid, method = "spearman", exact = FALSE)
cat(sprintf("Partial correlation after clade (%s): rho=%.3f, p=%.3g\n", pc_label, corresid$estimate, corresid$p.value))

p3 <- ggplot(dat, aes(x = presscore_clade_resid, y = pc_clade_resid, color = Pres_Quality)) + geom_point(size = 2, alpha = 0.7) +
stat_smooth(method = "lm", se = TRUE, color = "black") + scale_color_viridis(discrete = TRUE, end = 0.8, option = "D") +
xlab("Preservation quality (residuals after clade)") + ylab(sprintf("%s (residuals after clade)", pc_label)) +
annotate("text",
x = min(presscore_clade_resid, na.rm = TRUE),
y = max(pc_clade_resid, na.rm = TRUE),
label = sprintf("rho=%.2f, p=%.2g", corresid$estimate, corresid$p.value),
hjust = 0, size = 5) + theme_minimal(base_size = 13) + guides(color = "none")

#############################################################################################################
# The block is preparing the heatmap showing the mean PC score by exclusive clade and preservation quality
##############################################################################################################
hm <- dat %>%
group_by(clade, Pres_Quality) %>%
summarise(mean_pc = mean(.pc, na.rm = TRUE), n = n(), .groups = "drop")
hm$label <- with(hm, sprintf("%.1f\n(n=%d)", mean_pc, n))
hm$Pres_Quality <- factor(hm$Pres_Quality, levels = c("Low", "Medium", "Med-High", "High"))
hm$clade <- factor(hm$clade, levels = rev(sort(unique(dat$clade))))

p4 <- ggplot(hm, aes(x = Pres_Quality, y = clade, fill = mean_pc)) +
geom_tile(color = "white") + geom_text(aes(label = label), color = "black", size = 3.5) +
scale_fill_gradient2(low = "#1a9850", mid = "white", high = "#d73027", midpoint = 0, name = paste("Mean", pc_label)) +
scale_x_discrete(drop = FALSE) + theme_minimal(base_size = 12) +
labs(x = "Preservation quality", y = "Clade") + theme(axis.text.x = element_text(angle = 30, hjust = 1),panel.grid = element_blank())

#######################################################
# The block is coding lines saving the heatmap figure
######################################################
fig_box <- ggarrange(p1, p2, p3, ncol = 3, labels = c("A", "B", "C"))
ggsave(sprintf("output/plots_PDF/taphonomic_sensitivity_%s_box_corrs.pdf", tag),fig_box, width = 14, height = 5, dpi = 300)
ggsave(sprintf("output/plots_PDF/taphonomic_sensitivity_%s_heatmap.pdf", tag),p4, width = 8, height = 5, dpi = 300)

############################################################################################################################################
# The block is coding lines saving the statistical results into a csvof the taphonomic sensitivity across the biomechanical morphospace
############################################################################################################################################
hm_display <- hm
hm_display$mean_pc <- format_fixed(hm$mean_pc, 3)
names(hm_display)[names(hm_display) == "mean_pc"] <- paste0("mean_", pc_col)
write.csv(hm_display,
sprintf("output/results/Supplemental_performance/taphonomic_%s_by_clade_prescat.csv", tag), row.names = FALSE)
list(pc_col = pc_col, pc_label = pc_label, n = nrow(dat),
rho_raw = unname(corsimple$estimate), p_raw = corsimple$p.value,
rho_partial = unname(corresid$estimate), p_partial = corresid$p.value,
pc_clade_resid = pc_clade_resid, presscore_clade_resid = presscore_clade_resid)
}
res_pc1 <- run_taph_test("PC1", "PC1", "pc1")
res_pc2 <- run_taph_test("PC2", "PC2", "pc2")

###########################################################################
# The block is undertaking the comparison of the signals from PC1 and PC2
############################################################################
r23_raw <- cor(dat$PC1, dat$PC2, method = "spearman")
hw_raw <- hotelling_williams_test(res_pc1$rho_raw, res_pc2$rho_raw, r23_raw, nrow(dat))
r23_partial <- cor(res_pc1$pc_clade_resid, res_pc2$pc_clade_resid, method = "spearman")
hw_partial  <- hotelling_williams_test(res_pc1$rho_partial, res_pc2$rho_partial, r23_partial, nrow(dat))

################################################################################################
# The block is conducting Hotelling-Williams test for dependent correlations between PC1 and PC2
#################################################################################################
cat(sprintf("Raw signal: rho(PC1,Pres)=%.3f, rho(PC2,Pres)=%.3f, r(PC1,PC2)=%.3f -> t=%.3f, df=%d, p=%.4g -> %s\n",
res_pc1$rho_raw, res_pc2$rho_raw, r23_raw, hw_raw$t_stat, hw_raw$df, hw_raw$p_value,
ifelse(hw_raw$p_value < 0.05, "signal DIFFERS between PC1 and PC2", "no evidence of a difference")))
cat(sprintf("Partial signal: rho(PC1,Pres)=%.3f, rho(PC2,Pres)=%.3f, r(PC1,PC2)=%.3f -> t=%.3f, df=%d, p=%.4g -> %s\n",
res_pc1$rho_partial, res_pc2$rho_partial, r23_partial,hw_partial$t_stat, hw_partial$df, hw_partial$p_value,
ifelse(hw_partial$p_value < 0.05, "signal DIFFERS between PC1 and PC2", "no evidence of a difference")))
cat("NOTE: the Hotelling-Williams test formally compares Pearson correlations; it is\n")
cat("applied here to Spearman rho as a standard approximation and should be read\n")
cat("as indicative rather than exact.\n\n")

################################################################################################
# The block is preparing the csv containing the results of correlations test between PC1 and PC2
#################################################################################################
taph_pc_comparison <- data.frame( comparison = c("raw", "partial_after_clade"),
rho_PC1 = c(res_pc1$rho_raw, res_pc1$rho_partial), rho_PC2 = c(res_pc2$rho_raw, res_pc2$rho_partial),
r_PC1_PC2 = c(r23_raw, r23_partial), t_stat = c(hw_raw$t_stat, hw_partial$t_stat),
df = c(hw_raw$df, hw_partial$df), p_value = c(hw_raw$p_value, hw_partial$p_value))
write.csv(taph_pc_comparison,"output/results/Supplemental_performance/taphonomic_PC1_vs_PC2_signal_comparison.csv", row.names = FALSE)
