###################################################################################################################################
# The block is building the VMS background grid in PCA space. Preparation of a hared weighted-median grid builder
# This is used for both the performance (PC1/PC2) and EFA shape (shapePC1/2) morphospaces so stress surfaces are built the same way
###################################################################################################################################
weighted_median <- function(x, w) {
ok <- is.finite(x) & is.finite(w) & w > 0
if (!any(ok)) return(NA_real_)
x <- x[ok]; w <- w[ok]
o <- order(x); x <- x[o]; w <- w[o]
x[which(cumsum(w) >= sum(w) / 2)[1]]
}
build_weighted_stress_grid <- function(data, x_col, y_col, z_col,
pad_frac = 0.15, bw_frac = 0.08, n = 100) {
x <- data[[x_col]]; y <- data[[y_col]]; z <- data[[z_col]]
x_rng <- range(x, na.rm = TRUE); y_rng <- range(y, na.rm = TRUE)
padx <- diff(x_rng) * pad_frac; pady <- diff(y_rng) * pad_frac
g <- expand.grid(x = seq(x_rng[1] - padx, x_rng[2] + padx, length.out = n),
y = seq(y_rng[1] - pady, y_rng[2] + pady, length.out = n))
#####################################################################################################
# The lines are setting up  the smoothing width function from the size of the morphospace itself
####################################################################################################
bw <- bw_frac * max(diff(x_rng), diff(y_rng))
#######################################################################################################################
# The block is setting the cells far past the padded envelopefall back to the nearest specimen's value instead of NA
######################################################################################################################
g$z <- vapply(seq_len(nrow(g)), function(i) {
d2 <- (x - g$x[i])^2 + (y - g$y[i])^2
w <- exp(-d2 / (2 * bw^2))
if (sum(w) > 1e-6) weighted_median(z, w) else z[which.min(d2)]
}, numeric(1))
cat(sprintf("Weighted stress grid on (%s, %s): %d cells | bandwidth %.4g\n", x_col, y_col, nrow(g), bw))
list(grid = g, bandwidth = bw)
}
pc1_rng <- range(performance_data_clean$PC1, na.rm = TRUE)
pc2_rng <- range(performance_data_clean$PC2, na.rm = TRUE)
perf_grid <- build_weighted_stress_grid(performance_data_clean, "PC1", "PC2", BACKGROUND_VAR)
stress_grid <- perf_grid$grid
names(stress_grid)[names(stress_grid) == "x"] <- "PC1"
names(stress_grid)[names(stress_grid) == "y"] <- "PC2"
names(stress_grid)[names(stress_grid) == "z"] <- "stress"
bw <- perf_grid$bandwidth
stress_grid$viable <- inside_envelope(stress_grid$PC1, stress_grid$PC2)
cat(sprintf("Grid: %d cells | %d inside the envelope\n",nrow(stress_grid), sum(stress_grid$viable)))
cat(sprintf("Kernel bandwidth: %.4g | weighted median, nearest-neighbour fallback\n", bw))

############################################################
# The block is building the scale for the background surface
############################################################
sv <- stress_grid$stress[is.finite(stress_grid$stress)]
if (length(sv) == 0) stop("Background grid is entirely NA.")
surf_limits <- unname(quantile(sv, c(0.01, 0.99)))
if (surf_limits[1] == surf_limits[2]) surf_limits <- range(sv)
use_log <- all(sv > 0) && (surf_limits[2] / surf_limits[1] > 10)
###################################################################
# The block is setting up the colour palette for the VMS background
###################################################################
STRESS_PALETTES <- list( dark  = c("#5184a3","#327097","#19568a","#053d77","#052451","#020f23"),light = c("#64a2c7",
"#3d88b8","#1f69a8","#074b91","#072d64","#03142b"))
scale_fill_stress <- function(palette = "dark", name = surf_label) {
palette <- match.arg(palette, names(STRESS_PALETTES))
ggplot2::scale_fill_gradientn( colors   = STRESS_PALETTES[[palette]],
name = name, na.value = "white", trans = if (use_log) "log10" else "identity",
limits = surf_limits,  oob = scales::squish,
guide = ggplot2::guide_colourbar( barheight = grid::unit(60, "mm"), barwidth = grid::unit(5,  "mm"),
ticks.colour = "grey30",frame.colour = "grey30", title.position = "top",
label.theme    = ggplot2::element_text(size = 9, hjust = 0)))
}
cat(sprintf("Scale: %s | limits %.4g to %.4g\n\n",
if (use_log) "log10" else "linear", surf_limits[1], surf_limits[2]))
ENVELOPE_CAPTION <- paste0("Surface: ", gsub("\n", " ", surf_label), ".")
