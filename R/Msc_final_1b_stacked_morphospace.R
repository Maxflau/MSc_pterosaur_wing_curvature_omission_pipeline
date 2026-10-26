#################################################################
# The block is saving the function to make the stack morphospaces
#################################################################
save_stack <- function(expr, name, w = 9, h = 11.5) {
for (fmt in c("png", "pdf")) {
path <- sprintf("output/plots%s/%s.%s",
if (fmt == "pdf") "_PDF" else "", name, fmt)
if (fmt == "png") png(path, width = w, height = h, units = "in", res = 300) else pdf(path, width = w, height = h)
layout(matrix(1:2, ncol = 1), heights = c(5.2, 1))
par(mar = c(1.5, 3.5, 2.5, 1), xpd = TRUE); eval(expr)
par(mar = c(0, 1, 0, 1)); plot.new()
legend("center", legend = names(PAL), pt.bg = PAL, pch = 21, col = "grey15",
ncol = 3, bty = "o", cex = 0.82, pt.cex = 1.3, title = "Clade")
dev.off() }
cat("  Written:", name, "\n")
}
pct <- function(p, k) if (is.null(p)) "" else
sprintf(" (%.1f%%)", summary(p)$importance[2, k] * 100)

#####################################################################################
# The block is applying the stratigraphic stacked morphospace for the EFA morphospace
#####################################################################################
if (exists("shape_scores") && all(c("shapePC1","shapePC2") %in% names(shape_scores))) {
sp <- if ("species" %in% names(shape_scores)) "species" else names(shape_scores)[1]
sd_ <- merge(temporal_data, shape_scores[, c(sp, "shapePC1", "shapePC2")], by.x = "species", by.y = sp)
save_stack(quote(stacked_panels(sd_, "shapePC1", "shapePC2",
paste0("Shape PC1", pct(if (exists("pca_shape")) pca_shape else NULL, 1)),
paste0("Shape PC2", pct(if (exists("pca_shape")) pca_shape else NULL, 2)),
"EFA shape morphospace through time")), "STACK_01_shape_morphospace_time")
} else cat("shape_scores unavailable - shape stack skipped.\n")

#############################################################################################
# The block is applying the stratigraphic stacked morphospace for the performance morphospace
#############################################################################################
if (all(c("PC1","PC2") %in% names(temporal_data))) {
save_stack(quote(stacked_panels(temporal_data, "PC1", "PC2",
paste0("PC1", pct(if (exists("pca_performance")) pca_performance else NULL, 1)),
paste0("PC2", pct(if (exists("pca_performance")) pca_performance else NULL, 2)),
"Performance morphospace through time")),"STACK_02_performance_morphospace_time")
} else cat("PC1/PC2 absent - performance stack skipped.\n")
