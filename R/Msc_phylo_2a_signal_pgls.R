
library(nlme)

######################################################################
# The block is picking the variables that can be tested on the tree
######################################################################
SIGNAL_VARS <- intersect(
c("aspect_ratio", "r2_hat", "von_mises_stress",
"trailing_edge", "shape_complexity", "PC1", "PC2", "pareto_rank_ratio"),
colnames(phylo_data))

######################################################################################################################################
# The block is computing phylogenetic signal 
# (lambda: 0 means no tree pattern, 1 means Brownian motion; K below 1 means relatives are less alike than Brownian motion predicts)
######################################################################################################################################
signal_results <- data.frame()

for (v in SIGNAL_VARS) {
x <- phylo_data[[v]]
##############################################################################################
# The block is logging the metrics; PC scores and the rank ratio are already on a usable scale
##############################################################################################
if (!v %in% c("PC1", "PC2", "pareto_rank_ratio") && all(x > 0, na.rm = TRUE)) {
x <- log(x)
}
names(x) <- phy_pruned$tip.label
x <- x[is.finite(x)]
if (length(x) < 10) next

tree_v <- drop.tip(phy_pruned, setdiff(phy_pruned$tip.label, names(x)))

lam <- tryCatch(phylosig(tree_v, x, method = "lambda", test = TRUE),
error = function(e) NULL)
kk  <- tryCatch(phylosig(tree_v, x, method = "K", test = TRUE, nsim = 1000),
error = function(e) NULL)

signal_results <- rbind(signal_results, data.frame(
variable = v, n = length(x),
lambda   = if (!is.null(lam)) round(lam$lambda, 3) else NA,
lambda_p = if (!is.null(lam)) signif(lam$P, 3)      else NA,
K        = if (!is.null(kk))  round(kk$K, 3)        else NA,
K_p      = if (!is.null(kk))  signif(kk$P, 3)       else NA,
stringsAsFactors = FALSE))
}

cat("PHYLOGENETIC SIGNAL:\n")
######################################################################
# The block is showing a fixed number of decimals so every row lines up
######################################################################
signal_results_display <- signal_results
if (nrow(signal_results_display) > 0) {
signal_results_display$lambda   <- format_fixed(signal_results$lambda, 3)
signal_results_display$lambda_p <- format_fixed(signal_results$lambda_p, 4)
signal_results_display$K        <- format_fixed(signal_results$K, 3)
signal_results_display$K_p      <- format_fixed(signal_results$K_p, 4)
}
print(signal_results_display, row.names = FALSE)
cat("\nlambda near 1 with a significant p means related taxa resemble each other.\n")
cat("K below 1 means evolution is more labile than Brownian motion, which is\n")
cat("compatible with ecological convergence on a phylogenetically structured\n")
cat("background. Report both, they answer different questions.\n\n")
