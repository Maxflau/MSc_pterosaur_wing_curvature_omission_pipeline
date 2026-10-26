######################################################################################################################
# The block is ccomparing wing aspect ratio/r2_hat/von_mises_stress from the uncurved  wing outline
# With those predicted  from skeleton alone to evaluate the gap from wing reconsruction and morphological information 
######################################################################################################################
if (!exists("POSTCRANIAL")) stop("Source Msc_stats_8_allometry_data_prep.R first.")
for (dd in c("output/plots_PDF")) dir.create(dd, showWarnings = FALSE, recursive = TRUE)

AERO_TARGETS <- intersect(c("aspect_ratio", "r2_hat", "von_mises_stress"),
colnames(POSTCRANIAL))
MIN_FIT <- 10

panels_pva <- list(); panels_gap <- list(); gap_rows <- list(); imp_rows <- list()

for (tg in AERO_TARGETS) {
d <- POSTCRANIAL[, c(tg, "clade", GM_ELEMENTS)] %>% drop_na()
if (nrow(d) < MIN_FIT) { cat(sprintf("Skipping %s (n < %d)\n", tg, MIN_FIT)); next }

rf <- tryCatch(
randomForest(as.formula(paste(tg, "~", paste(GM_ELEMENTS, collapse = "+"))),
data = d, ntree = 500, importance = TRUE),
error = function(e) NULL)
if (is.null(rf)) next

###########################################################################
# The block is predictin skeleton "reconstructed" values of each specimens
###########################################################################
d$reconstructed <- predict(rf)
d$gap <- d$reconstructed - d[[tg]]
r2 <- round(cor(d[[tg]], d$reconstructed, use = "complete.obs")^2, 3)

imp_df <- as.data.frame(importance(rf)); imp_df$element <- rownames(imp_df)
imp_df$target <- tg; imp_rows[[tg]] <- imp_df

gap_summary <- d %>% group_by(clade) %>%
summarise(n = n(), mean_gap = round(mean(gap, na.rm = TRUE), 5),
sd_gap = round(sd(gap, na.rm = TRUE), 5), .groups = "drop")
gap_summary$target <- tg; gap_rows[[tg]] <- gap_summary

panels_pva[[tg]] <- ggplot(d, aes(.data[[tg]], reconstructed)) +
geom_point(aes(colour = clade), alpha = 0.7, size = 2) +
geom_abline(aes(slope = 1, intercept = 0, linetype = "1:1 (skeletal based)"),
colour = "grey40") +
geom_smooth(aes(linetype = "Fitted trend"), method = "lm", se = FALSE, colour = "black") +
scale_linetype_manual(name = NULL,
values = c("1:1 (skeletal based)" = "dashed",
"Fitted trend" = "dotted")) +
labs(title = tg,
subtitle = bquote(R^2 == .(r2) ~ "(skeleton-only OOB reconstruction vs. outline-measured value)"),
x = paste("Actual", tg, "(measured from the digitized wing outline)"),
y = paste("Skeleton-reconstructed", tg)) +
theme_bw(base_size = 10) + theme(legend.position = "bottom")

panels_gap[[tg]] <- ggplot(d, aes(fct_infreq(as.factor(clade)), gap)) +
geom_boxplot(aes(fill = clade), alpha = 0.55, outlier.shape = NA) +
geom_jitter(size = 1, alpha = 0.45, width = 0.2) +
geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
labs(title = paste(tg, "gap by clade"),
subtitle = "Positive = skeleton over-predicts relative to the real wing outline",
x = "Clade", y = "Gap (skeleton-reconstructed - actual)") +
theme_bw(base_size = 9) +
theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

#####################################################################
# The block is abuilding a combined figure with reconstructed-vs-actual 
# on the top row and gap by clade on the bottom row
#####################################################################
if (length(panels_pva) > 0) {
top_row <- wrap_plots(panels_pva, nrow = 1, guides = "collect") &
theme(legend.position = "bottom")
bottom_row <- wrap_plots(panels_gap, nrow = 1)

p_combined <- (top_row / bottom_row) +
plot_layout(heights = c(1.1, 1)) +
plot_annotation(
title = "Skeleton-only reconstruction vs. actual (curved-wing) aerodynamic properties",
caption = paste("Top: each aerodynamic metric predicted from skeletal elements alone",
"(out-of-bag random-forest reconstruction, no outline/curvature information)",
"against its actual value, measured directly from the real, digitized wing outline;",
"dashed = 1:1 perfect reconstruction, dotted = fitted trend.",
"Bottom: reconstructed-minus-actual gap per clade; a positive gap means",
"the skeleton-only reconstruction overestimates performance relative to",
"the real, digitized wing outline."))
ggsave("output/plots_PDF/BIOMECH_GAP_01_skeleton_reconstruction_vs_actual.pdf",
p_combined, width = 5 * length(panels_pva), height = 10, limitsize = FALSE)
cat(sprintf("Written: BIOMECH_GAP_01_skeleton_reconstruction_vs_actual (%d targets)\n",
length(panels_pva)))
}

if (length(imp_rows) > 0) write_supp(do.call(rbind, imp_rows), "S15b_skeleton_aero_rf_importance")
if (length(gap_rows) > 0) write_supp(do.call(rbind, gap_rows), "S15b_skeleton_aero_gap_by_clade")
cat("\n")
