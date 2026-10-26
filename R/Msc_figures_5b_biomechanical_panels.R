#################################################################################
# The block is making the structure of every  biomechanical PCA panel by metrics
#################################################################################
make_biomech_panel <- function(colour_var, panel_title, palette = "plasma",log_colour = FALSE, use_stress_palette = FALSE, hollow_points = FALSE)
{
if (!colour_var %in% colnames(performance_data_clean)) {
cat("SKIPPED", panel_title, "- column", colour_var, "does not exist.\n")
cat("  Present:", paste(intersect(names(METRIC_LABELS), colnames(performance_data_clean)), collapse = ", "), "\n")
return(NULL)
}
d <- performance_data_clean
v <- d[[colour_var]]
d$.clade <- if ("clade" %in% colnames(d)) as.character(d$clade) else NA_character_
if (log_colour && any(v <= 0, na.rm = TRUE)) {
cat(sprintf("  %s: %d non-positive values - linear scale used.\n",colour_var, sum(v <= 0, na.rm = TRUE)))
log_colour <- FALSE
}
d$.col <- if (log_colour) log(v) else v
d <- d[is.finite(d$.col), ]
if (nrow(d) < 5) {
cat("SKIPPED", panel_title, "- only", nrow(d), "finite values.\n"); return(NULL)
}
pt_label_raw <- if (colour_var %in% names(METRIC_LABELS)) {
METRIC_LABELS[[colour_var]]
} else colour_var
##################################################################
# The block is building the fill of the data point and metric grid
#################################################################
if (log_colour) pt_label_raw <- paste0("log(", gsub("\n", " ", pt_label_raw), ")")
pt_label <- paste0("Point fill - ", gsub("\n", " ", pt_label_raw))
bg_label <- paste0("Background - ", gsub("\n", " ", surf_label))
cat(sprintf("  %s: %d specimens | legend: %s | points: %s\n", colour_var, nrow(d), gsub("\n", " ", pt_label_raw),
if (hollow_points) "hollow (clade border only)" else "filled"))
p <- ggplot() +geom_tile(data = stress_grid, aes(x = PC1, y = PC2, fill = stress),alpha = 0.7) +
scale_fill_stress(name = bg_label, palette = "dark")
if (!is.null(wing_bg)) {
p <- p + geom_polygon(data = wing_bg, aes(x = PC1, y = PC2, group = grid_id), fill = "#2a2a2a", colour = "#1a1a1a",alpha = 0.18, linewidth = 0.25)
}
p <- p + ggnewscale::new_scale_colour()
if (has_hulls) {
p <- p + geom_path(data = hulls_clade, aes(x = PC1, y = PC2, group = grp, colour = grp), linewidth = 1.0, alpha = 0.9, show.legend = FALSE)
}
if (hollow_points) {
p <- p + geom_point(data = d, aes(x = PC1, y = PC2, colour = .clade), shape = 21, fill = NA, size = 3.2, stroke = 1.0, alpha = 0.9) +
scale_colour_manual(values = CLADE_COLS_MATCHED, na.value = "grey50",name = "Point border - Clade")
} else {
p <- p + ggnewscale::new_scale_fill()
fill_scale <- if (use_stress_palette) {
ggplot2::scale_fill_gradientn(colors = STRESS_PALETTES$dark, name = pt_label, na.value = "white")
} else { scale_fill_viridis_c(option = palette, name = pt_label)
}
p <- p + geom_point(data = d, aes(x = PC1, y = PC2, fill = .col, colour = .clade),shape = 21, size = 3.2, stroke = 0.9, alpha = 0.9) +
fill_scale + scale_colour_manual(values = CLADE_COLS_MATCHED, na.value = "grey50",name = "Point border - Clade")
}
p + labs(title = panel_title, subtitle = paste0("Surface: ", gsub("\n", " ", surf_label),". Dark outlines: wing profiles across the morphospace.",
if (hollow_points) " Points: clade border only, no fill." else ""),
x = pc_lab(1), y = pc_lab(2)) + theme_minimal(base_size = 14) +
theme(plot.title = element_text(hjust = 0, face = "bold", size = 16),
plot.subtitle = element_text(hjust = 0, size = 9, colour = "grey35"), legend.position = "right",
panel.grid.major = element_line(colour = "grey90"), panel.background = element_rect(fill = "white", colour = NA),
plot.background  = element_rect(fill = "white", colour = NA),panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.6))
}
####################################################
# The block is making the three PCA panel by metrics
####################################################
if (!"von_mises_stress" %in% colnames(performance_data_clean)) {
stop("von_mises_stress not found in performance_data_clean. ","Source Msc_performance_metrics_v3.R first.")
}
stress_col <- "von_mises_stress"
PANELS <- list(list(var = "aspect_ratio", title = "(a) Aspect ratio", pal = "plasma", log = TRUE, stress = FALSE,
hollow = FALSE, file = "CLADE_01A_aspect_ratio"),
list(var = "r2_hat", title = "(b) Second moment of area", pal = "viridis", log = FALSE, stress = FALSE,
hollow = FALSE, file = "CLADE_01B_second_moment"),
list(var = stress_col, title = "(c) Von Mises stress", pal = "mako", log = TRUE, stress = TRUE, hollow = TRUE, file = "CLADE_01C_von_mises_stress"))
plots <- list()
for (pn in PANELS) {
if (is.na(pn$var)) next
p <- make_biomech_panel(pn$var, pn$title, pn$pal, pn$log, pn$stress, pn$hollow)
if (is.null(p)) next
plots[[pn$file]] <- p
ggsave(sprintf("output/plots_PDF/biomechanic_plot/%s.pdf", pn$file), p, width = 18, height = 11)
cat("Biomechanic panel written:", pn$file, "\n")
}
