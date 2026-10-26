########################################################################################################
# The block is coding  setting up the exclusive clade colours and build a smooth background VMS surface
########################################################################################################
CLADE_COLS <- c("Anurognathidae" = "#46852F", "Azhdarchoidea"= "#DB8E00", "Basal pterodactyliform" = "#F4A7B9",
"Most Basal pterosaur" = "#0CE27F","Euctenochasmatia"= "#9F1212","Darwinoptera"= "#791444","Dsungaripteroidea" = "#60A293",
"Eudimorphodontoidea"= "#FF00F0", "Ornithocheiromorpha"= "#FFE413","Pteranodontia"= "#DA8AAA","Rhamphorhynchidae"= "#03240B")

###############################################################################################
# The block is running a safeguard to make sure the colour match the assigned exclusive clades
###############################################################################################
match_clade_colour <- function(clade_names) {
norm <- function(x) tolower(trimws(x))
lut <- CLADE_COLS
names(lut) <- norm(names(CLADE_COLS))
out <- unname(lut[norm(clade_names)])
names(out) <- clade_names
out
}

####################################################################
# The block is buidling the axis labels from the current PCA output
####################################################################
pc_lab <- function(k) {
sprintf("PC%d (%.1f%%)", k, summary(pca_performance)$importance[2, k] * 100)
}
if (!requireNamespace("mgcv", quietly = TRUE)) {
stop("Package 'mgcv' is required for the smoothed stress surface: install.packages(\"mgcv\")")
}
################################################################################################
# The block is smoothing the stress values across the morphospace grid for the background layer
################################################################################################
build_smooth_stress_grid <- function(value_col = "von_mises_stress") {
df <- data.frame(x = performance_data_clean$PC1, y = performance_data_clean$PC2,
z = performance_data_clean[[value_col]])
df <- df[complete.cases(df) & is.finite(df$z) & df$z > 0, ]
if (nrow(df) < 15) { cat("  Too few points for a smoothed stress surface.\n"); return(NULL) }
m <- mgcv::gam(log(z) ~ s(x, y, k = min(60, nrow(df) - 1)), data = df)
g <- data.frame(PC1 = stress_grid$PC1, PC2 = stress_grid$PC2)
g$stress <- exp(predict(m, newdata = data.frame(x = g$PC1, y = g$PC2)))
g
}
stress_grid_smooth <- build_smooth_stress_grid("von_mises_stress")
if (!is.null(stress_grid_smooth)) {
cat(sprintf("Smoothed stress surface (GAM) built: %d cells, range %.4g to %.4g\n\n",
nrow(stress_grid_smooth), min(stress_grid_smooth$stress), max(stress_grid_smooth$stress)))
}
