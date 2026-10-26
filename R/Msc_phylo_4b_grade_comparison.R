#############################################################################################
# The block is comparing the Pareto scores of the two broad grades in one summary violin plot
#############################################################################################
vdf <- as.data.frame(pd)
vdf$grade <- ifelse(vdf$clade %in% PTERODACT,
"Pterodactyliform", "non-pterodactyliform")

##################################################################################
# The lines are naing in grade or opt does not merely drop rows:
# a violin whose group is NA is drawn with the na.value colour (grey)
##################################################################################
vdf <- vdf[!is.na(vdf$grade) & is.finite(vdf$opt), ]

###############################################################################################################
# The block is explicitting factor levels, so ggplot cannot re-sort them out of step with the colour vector
###############################################################################################################
vdf$grade <- factor(vdf$grade,
levels = c("non-pterodactyliform", "Pterodactyliform"))

cat("\nGrade levels:\n"); print(table(vdf$grade, useNA = "always"))

if (nlevels(droplevels(vdf$grade)) == 2) {
GRADE_COLS <- c("non-pterodactyliform" = "#9575CD",
"Pterodactyliform"     = "#FFE413")

wx <- wilcox.test(opt ~ grade, data = vdf)
cat(sprintf("Wilcoxon: W = %.0f, p = %.3e\n", wx$statistic, wx$p.value))

yr_v <- range(vdf$opt, na.rm = TRUE)

p_violin <- ggplot(vdf, aes(x = grade, y = opt)) +
###########################################################################
# The block is filling inside the geom's own aes(), which sets fill="white"
###########################################################################
geom_violin(aes(fill = grade), alpha = 0.6, trim = FALSE,
colour = "grey25", linewidth = 0.4, scale = "width") +
geom_boxplot(width = 0.12, fill = "white", colour = "grey25",
alpha = 0.9, outlier.shape = NA, linewidth = 0.4) +
geom_jitter(width = 0.08, alpha = 0.35, size = 1.4, colour = "grey30") +
scale_fill_manual(values = GRADE_COLS, drop = FALSE, guide = "none") +
annotate("text", x = 1.5, y = yr_v[2] - 0.03 * diff(yr_v),
label = sprintf("Wilcoxon\nW = %.0f\np = %.2e",
wx$statistic, wx$p.value),
size = 3, hjust = 0.5, family = "mono") +
labs(title = "(c) Pareto optimality by grade",
subtitle = sprintf("Performance morphospace | n = %d vs %d",
sum(vdf$grade == levels(vdf$grade)[1]),
sum(vdf$grade == levels(vdf$grade)[2])),
x = NULL, y = "Pareto rank ratio") +
theme_classic(base_size = 12) +
theme(plot.title = element_text(face = "bold"),
plot.subtitle = element_text(size = 9, colour = "grey35"))

ggsave("output/plots_PDF/PERF_15C_optimality_violin.pdf",
p_violin, width = 7, height = 8)
cat("Written: PERF_15C_optimality_violin\n")
} else {
cat("Only", nlevels(droplevels(vdf$grade)), "grade(s) - violin skipped.\n")
print(unique(vdf$clade))
}
cat("\n")
