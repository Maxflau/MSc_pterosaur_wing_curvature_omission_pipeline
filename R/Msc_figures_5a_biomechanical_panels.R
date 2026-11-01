################################################################################################################
# The block is preparing the biomechanical morphospace panels which are stored in the biomechanical plot folder
################################################################################################################
dir.create("output/plots_PDF/biomechanic_plot", showWarnings = FALSE, recursive = TRUE)
METRIC_LABELS <- c( von_mises_stress = "Von Mises stress\n(scale-invariant)", r2_hat= "Second moment\nof area",
aspect_ratio = "Aspect ratio", trailing_edge = "Trailing edge", shape_complexity  = "Shape complexity",
pareto_rank_ratio = "Pareto\nrank ratio")
surf_label <- if (BACKGROUND_VAR %in% names(METRIC_LABELS)) {
METRIC_LABELS[[BACKGROUND_VAR]]
} else BACKGROUND_VAR
cat("Surface variable:", BACKGROUND_VAR, "| colour bar reads:", gsub("\n", " ", surf_label), "\n")

#################################################################################
# The block is building the clade outlines from the points shown within the plot
#################################################################################
build_hulls <- function(group_col = "clade") {
if (!group_col %in% colnames(performance_data_clean)) return(NULL)
d <- performance_data_clean[!is.na(performance_data_clean[[group_col]]) &
is.finite(performance_data_clean$PC1) &is.finite(performance_data_clean$PC2), ]
if (nrow(d) < 3) return(NULL)
do.call(rbind, lapply(unique(d[[group_col]]), function(g) {
dg <- d[d[[group_col]] == g, ]
if (nrow(dg) < 3) return(NULL)
h <- chull(dg$PC1, dg$PC2); h <- c(h, h[1])
data.frame(PC1 = dg$PC1[h], PC2 = dg$PC2[h], grp = g)
}))
}
###########################################################################
# The block is building the convex hulls of the exclusive clades
###########################################################################
hulls_clade <- build_hulls("clade")
has_hulls <- !is.null(hulls_clade) && nrow(hulls_clade) > 0
cat(sprintf("Clade hulls: %d groups\n\n",
if (has_hulls) length(unique(hulls_clade$grp)) else 0))
ALL_CLADES <- if ("clade" %in% colnames(performance_data_clean)) {
unique(as.character(performance_data_clean$clade))
} else character(0)
CLADE_COLS_MATCHED <- match_clade_colour(ALL_CLADES)
unmapped_cl <- ALL_CLADES[is.na(CLADE_COLS_MATCHED)]
if (length(unmapped_cl) > 0) {
cat("NOTE: no CLADE_COLS entry for:", paste(unmapped_cl, collapse = ", "),"- these get a grey border.\n")
CLADE_COLS_MATCHED[is.na(CLADE_COLS_MATCHED)] <- "grey50"
}
cat("Clade colour mapping in use:\n")
print(CLADE_COLS_MATCHED)
cat("\n")