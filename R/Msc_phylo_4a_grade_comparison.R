
FIG_W <- 9; FIG_H <- 9
###############################################################################################
# The lines are rebuilding the subgroup labels and drawing one phylomorphospace panel at a time
###############################################################################################
ax1 <- sprintf("PC1 (%.1f%%)", summary(pca_performance)$importance[2, 1] * 100)
ax2 <- sprintf("PC2 (%.1f%%)", summary(pca_performance)$importance[2, 2] * 100)

################################################################################
# The block is rebuilding the panel groups from the clade names each time 
# avoiding empty panels if saved labels differ by even a small spelling change
################################################################################
if (!"clade" %in% names(pd)) stop("pd has no clade column.")
pd <- as.data.frame(pd)   # a tbl_df with list columns breaks geom_point
pd$panel_group <- ifelse(pd$clade %in% PTERODACT,
"Pterodactyliform", "Non-pterodactyliform")

cat("\nClades found in pd:\n"); print(table(pd$clade, useNA = "ifany"))
cat("\nPanel groups:\n");       print(table(pd$panel_group))

unknown <- setdiff(unique(pd$clade),
c(PTERODACT, "Darwinoptera", "Anurognathidae",
"Rhamphorhynchidae", "Non-breviquartossan", "Most Basal pterosaur"))
if (length(unknown) > 0) {
cat("\nClades not in PTERODACT and assumed basal:\n"); print(unknown)
}

make_panel <- function(group_name, letter, subtitle) {
grp <- pd[pd$panel_group == group_name, ]
cat(sprintf("\n  %s: %d tips\n", group_name, nrow(grp)))
if (nrow(grp) == 0) {
cat("  WARNING: no tips in this group - landscape and grey tips only.\n")
}
######################################################################################################
# The block is grading colour the colour of the phylomorphospace and taxonomic grouping of interests
######################################################################################################
grp_col <- if (group_name == "Pterodactyliform") "#FFE413" else "#FF00F6"

p <- ggplot() +
geom_raster(data = grid[!is.na(grid$opt), ],
aes(x = PC1, y = PC2, fill = opt), interpolate = TRUE) +
scale_fill_gradientn(
colours = c("#ffffff","#e8f4ea","#bfe0c6","#8fc79c","#54a468","#1b7837"),
name = "Pareto\noptimality", limits = c(0, 1), na.value = "white")

if (!is.null(wing_perf)) {
p <- p + geom_polygon(data = wing_perf, aes(x = PC1, y = PC2, group = grid_id),
fill = "grey55", colour = "grey45",
alpha = 0.15, linewidth = 0.18)
}

p <- p +
geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
colour = "grey35", alpha = 0.45, linewidth = 0.28) +
############################################################################
# The block is everying tip in pale grey first, so no branch ends in mid-air
# Edges covers the whole tree while grp is only the subgroup
############################################################################
geom_point(data = pd, aes(x = PC1, y = PC2),
shape = 21, fill = "grey82", colour = "grey50",
size = 1.6, stroke = 0.3)

if (nrow(grp) > 0) {
p <- p + geom_point(data = grp, aes(x = PC1, y = PC2),
shape = 21, fill = grp_col, colour = "grey15",
size = 1.7 , stroke = 0.4)
}

######################################################################
# The block is geoming_raster sets the panel limits from the grid alone; 
# tips falling outside it would be clipped and vanish without warning
######################################################################
p +
coord_fixed(ratio = 0.75, expand = FALSE,
xlim = range(c(grid$PC1, pd$PC1), na.rm = TRUE),
ylim = range(c(grid$PC2, pd$PC2), na.rm = TRUE)) +
labs(title = sprintf("(%s) %s Phylomorphospace", letter, group_name),
subtitle = sprintf("%s | %d of %d tips highlighted",
subtitle, nrow(grp), nrow(pd)),
x = ax1, y = ax2) +
theme_bw(base_size = 10) +
theme(plot.title = element_text(face = "bold", size = 16),
plot.subtitle = element_text(size = 9.5, colour = "grey35"),
panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.6),
panel.grid = element_blank())
}

############################################################################
# The block is writing the two tree panels after the function is defined,
# keeping this loop outside the function so it does not call itself forever
############################################################################
for (pn in list(
list(g = "Non-pterodactyliform", l = "a",
s = "Basal pterosaurs on the performance optimality landscape",
f = "PERF_15A_phylo_nonpterodact"),
list(g = "Pterodactyliform", l = "b",
s = "Derived pterosaurs on the performance optimality landscape",
f = "PERF_15B_phylo_pterodact"))) {
p <- make_panel(pn$g, pn$l, pn$s)
ggsave(sprintf("output/plots_PDF/%s.pdf", pn$f), p,
width = FIG_W, height = FIG_H)
cat("Written:", pn$f, "\n")
}
