######################################################################
# The block is summarizing how much trailing edge variation is present
######################################################################
if (!exists("write_supp")) stop("Source Msc_stats_0_setup.R first.")
if (!exists("outlines_list")) stop("outlines_list not found.")

for (d in c("output/plots", "output/plots_PDF")) {
dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

#############################################################
# The block is summarizing the measured trailing edge values
#############################################################
d <- performance_data_clean[is.finite(performance_data_clean$trailing_edge), ]
cat(sprintf("Specimens with a curvature value: %d of %d\n",
nrow(d), nrow(performance_data_clean)))
print(summary(d$trailing_edge))

CURV_LEVELS <- intersect(c("clade", "Order", "Diet_combined"), colnames(d))

dist_tab <- do.call(rbind, lapply(CURV_LEVELS, function(g) {
s <- d[!is.na(d[[g]]) & d[[g]] != "", ]
do.call(rbind, lapply(unique(s[[g]]), function(lv) {
x <- s$trailing_edge[s[[g]] == lv]
if (length(x) < MIN_N) return(NULL)
data.frame(level = g, group = lv, n = length(x),
mean = round(mean(x), 5), sd = round(sd(x), 5),
median = round(median(x), 5),
q05 = round(quantile(x, 0.05), 5),
q95 = round(quantile(x, 0.95), 5), stringsAsFactors = FALSE)
}))
}))
write_supp(dist_tab, "S23_curvature_distribution")

##################################################################
# The block is comparing each outline with its curvature-free hull
##################################################################
hull_metrics <- function(o) {
h <- chull(o$x, o$y)
oh <- data.frame(x = o$x[h], y = o$y[h])
oh <- rbind(oh, oh[1, ])
c(curv = calculate_trailing_edge(oh),
r2   = calculate_r2_hat(oh),
cplx = calculate_shape_complexity(oh),
ar   = calculate_aspect_ratio(oh))
}

hm <- t(vapply(performance_data_clean$species, function(s) {
o <- outlines_list[[s]]
if (is.null(o) || nrow(o) < 4) return(c(curv = NA, r2 = NA, cplx = NA, ar = NA))
tryCatch(hull_metrics(o), error = function(e) c(curv = NA, r2 = NA, cplx = NA, ar = NA))
}, numeric(4)))

cmp <- data.frame(
species = performance_data_clean$species,
clade = performance_data_clean$clade,
curvature_measured = performance_data_clean$trailing_edge,
curvature_hull = hm[, "curv"],
r2_measured = performance_data_clean$r2_hat, r2_hull = hm[, "r2"],
complexity_measured = performance_data_clean$shape_complexity,
complexity_hull = hm[, "cplx"], stringsAsFactors = FALSE)

cmp$curvature_lost <- cmp$curvature_measured - cmp$curvature_hull
cmp$r2_shift_pct <- 100 * (cmp$r2_hull - cmp$r2_measured) / cmp$r2_measured

ok <- is.finite(cmp$curvature_lost)

#################################################################################################################
# The block is checking that the difference is zero everywhere when outlines already went through chull() upstream
###################################################################################################################
if (all(abs(cmp$curvature_lost[ok]) < 1e-8)) {
# The outlines are already convex hulls (Msc_outline_extract_1.R applies chull()),
# so trailing_edge measures the bulge of a hull rather than membrane camber.
# This is stated as a methodological limitation; see methods/caveats text.
}

write_supp(cmp, "S24_curvature_free_comparison")

##########################################
# The block is saving the summary figures
##########################################
p1 <- ggplot(d, aes(x = trailing_edge)) +
geom_histogram(bins = 30, fill = "#2E8B57", colour = "grey30", linewidth = 0.25) +
geom_vline(xintercept = median(d$trailing_edge), linetype = "dashed",
colour = "grey20") +
labs(title = "Variation of the uncruved trailing edge across Pterosauria",
subtitle = sprintf("n = %d; dashed line: median %.4f",
nrow(d), median(d$trailing_edge)),
x = "Trailing edge (camber ratio)", y = "Specimens") + theme_bw(base_size = 12) + theme(panel.border = element_rect(colour = "grey30",
fill = NA, linewidth = 0.5))

ggsave("output/plots_PDF/CURVATURE_01_distribution.pdf", p1,
width = 8, height = 5)
