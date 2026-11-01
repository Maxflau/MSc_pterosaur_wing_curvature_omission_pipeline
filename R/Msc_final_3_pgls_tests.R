###################################################################
# The lines are conducting the PGLS test onto the metrics, PCA axes
###################################################################

dir.create("output/results/Supplemental_performance", showWarnings = FALSE, recursive = TRUE)

######################################################################################
# The block is making species names easier to match between the table and the tree
######################################################################################
normalise_name <- function(x) {
tolower(gsub("[^A-Za-z0-9_]", "", gsub("[ ]+", "_", trimws(as.character(x)))))
}
if (!exists("performance_data_clean")) stop("performance_data_clean not found.")

######################################################################
# The block is setting up the tree for the PGLS statistical analyses
######################################################################
if (!exists("phy_pruned")) {
cand <- ls(envir = .GlobalEnv)[vapply(ls(envir = .GlobalEnv), function(x)
inherits(get(x, envir = .GlobalEnv), "phylo"), logical(1))]
if (length(cand) > 0) {
phy_pruned <- get(cand[1]); cat("Using phylo object:", cand[1], "\n")
} else if (file.exists("output/phylogenetics/Thomas_updated.nex")) {
tr0 <- read.nexus("output/phylogenetics/Thomas_updated.nex")
phy_pruned <- if (inherits(tr0, "multiPhylo")) tr0[[1]] else tr0
cat(sprintf("Tree read: %d tips\n", Ntip(phy_pruned)))
} else stop("No tree available in ", getwd())
}

###########################################################################################
# The block is preparing a single metric set: von_mises_stress replaces 
# both the old stress_root (non-dimensional) and the older pixel-dependent von_mises_stress
###########################################################################################
CORE <- c("aspect_ratio", "r2_hat", "von_mises_stress")
if (!all(CORE %in% names(performance_data_clean))) {
stop("Metric set incomplete. Needs: ", paste(CORE, collapse = ", "),
"\n  Present: ", paste(intersect(CORE, names(performance_data_clean)), collapse = ", "))
}
SET_NAME <- "von_mises_stress"

VARS <- intersect(c(CORE, "trailing_edge","shape_complexity", "pareto_rank_ratio",
"optimality_mean", "optimality_geom", "PC1", "PC2"),
colnames(performance_data_clean))
cat("Metric set:", SET_NAME, "\nVariables:", paste(VARS, collapse = ", "), "\n")

###############################################
# The block is aligning the data with the tree
##############################################
d   <- performance_data_clean
key <- normalise_name(d$species); tip <- normalise_name(phy_pruned$tip.label)
common <- intersect(key, tip)
cat(sprintf("Matched: %d of %d\n\n", length(common), nrow(d)))
if (length(common) < 20) stop("Too few matches.")

d  <- d[key %in% common, ]
tr <- drop.tip(phy_pruned, phy_pruned$tip.label[!tip %in% common])
tr$tip.label <- normalise_name(tr$tip.label)
d  <- d[match(tr$tip.label, normalise_name(d$species)), ]
######################################################################################
# The lines are replacing zero-length branches with a tiny value so the model can run
######################################################################################
if (any(tr$edge.length <= 0)) tr$edge.length[tr$edge.length <= 0] <- 1e-6 * max(tr$edge.length)

###########################################################
# The block is trying a few starting lambdas since corPagel 
# Chokes when the true value is far off
# Falls back to Brownian if nothing fits
##########################################################
fit_pgls <- function(df, trs, vname) {
for (start in c(0.5, 0.9, 0.1, 0.99)) {
m <- tryCatch(gls(y ~ Midpoint, data = df,
correlation = corPagel(start, phy = trs, form = ~species),
method = "ML"), error = function(e) NULL)
if (!is.null(m)) {
lam <- as.numeric(coef(m$modelStruct$corStruct, unconstrained = FALSE))
return(list(model = m, lambda = lam, method = sprintf("Pagel (start %.2f)", start)))
}
}
######################################################################################################
# The block is checking Brownian results lambda at 1 and overcorrects, 
# but it converges where the optimiser will not: reported as such rather than silently substituted
######################################################################################################
m <- tryCatch(gls(y ~ Midpoint, data = df,
correlation = corBrownian(phy = trs, form = ~species),
method = "ML"), error = function(e) NULL)
if (!is.null(m)) return(list(model = m, lambda = 1, method = "Brownian (fallback)"))

cat("    all PGLS attempts failed for", vname, "\n")
list(model = NULL, lambda = NA_real_, method = "failed")
}
