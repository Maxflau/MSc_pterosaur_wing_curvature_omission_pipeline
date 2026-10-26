##############################################################################################
# The block is pglsing - test whether time trends remain after the tree is taken into account
##############################################################################################
run_pgls <- function(response) {
if (!all(c(response, "Midpoint") %in% colnames(phylo_data))) return(invisible(NULL))

df <- phylo_data[, c(response, "Midpoint")]
names(df)[1] <- "y"
ok <- complete.cases(df)
df <- df[ok, ]
tr <- drop.tip(phy_pruned, phy_pruned$tip.label[!ok])
df$species <- tr$tip.label

cat(sprintf("--- %s ~ Midpoint ---\n", response))
cat("Uncorrected OLS:\n")
print(round(summary(lm(y ~ Midpoint, data = df))$coefficients, 5))

m_lam <- tryCatch(
gls(y ~ Midpoint, data = df,
correlation = corPagel(1, phy = tr, form = ~species), method = "ML"),
error = function(e) { cat("corPagel failed:", e$message, "\n"); NULL })

if (!is.null(m_lam)) {
cat("\nPGLS, Pagel's lambda estimated:\n")
print(round(summary(m_lam)$tTable, 5))
cat(sprintf("Estimated lambda: %.3f\n\n",
as.numeric(coef(m_lam$modelStruct$corStruct, unconstrained = FALSE))))
}
invisible(m_lam)
}

for (resp in intersect(c("pareto_rank_ratio", "PC1", "aspect_ratio","PC2", "r2_hat","von_mises_stress"),
colnames(phylo_data))) {
run_pgls(resp)
}

######################################################################
# The lines are testing clade differences against the phylogenetic null
######################################################################

##############################################################################################
# The block is cladessing are defined by the tree itself, so use a phylogenetic test
# That compares the observed clade differences with the right tree-based null model
##############################################################################################
if ("clade" %in% colnames(phylo_data)) {
y   <- setNames(phylo_data$PC1, phy_pruned$tip.label)
grp <- setNames(as.factor(phylo_data$clade), phy_pruned$tip.label)
ok  <- is.finite(y) & !is.na(grp)

pa <- tryCatch(
phylANOVA(drop.tip(phy_pruned, phy_pruned$tip.label[!ok]),
grp[ok], y[ok], nsim = 1000, posthoc = FALSE),
error = function(e) { cat("phylANOVA failed:", e$message, "\n"); NULL })

if (!is.null(pa)) {
cat(sprintf("PHYLOGENETIC ANOVA (PC1 ~ clade): F = %.3f, p = %.4f\n\n",
pa$F, pa$Pf))
}
}

######################################################################
# The block is testing robustness across grafted trees
######################################################################
if (!is.null(all_trees) && length(all_trees) > 1 &&
all(c("pareto_rank_ratio", "Midpoint") %in% colnames(phylo_data))) {
n_trees <- min(length(all_trees), 100)
cat(sprintf("Repeating the PGLS across %d trees...\n", n_trees))

slopes <- numeric(0); pvals <- numeric(0)
for (i in seq_len(n_trees)) {
tr <- all_trees[[i]]
kt <- normalise_name(tr$tip.label) %in% phylo_data$tip_key
tr <- drop.tip(tr, tr$tip.label[!kt])
dd <- phylo_data[normalise_name(tr$tip.label), c("pareto_rank_ratio", "Midpoint")]
names(dd)[1] <- "y"
if (any(!complete.cases(dd))) next
dd$species <- tr$tip.label
m <- tryCatch(gls(y ~ Midpoint, data = dd,
correlation = corPagel(1, phy = tr, form = ~species), method = "ML"),
error = function(e) NULL)
if (!is.null(m)) {
slopes <- c(slopes, summary(m)$tTable["Midpoint", "Value"])
pvals  <- c(pvals,  summary(m)$tTable["Midpoint", "p-value"])
}
}

if (length(slopes) > 0) {
cat(sprintf("Successful fits: %d\n", length(slopes)))
print(summary(slopes))
cat(sprintf("Trees with p < 0.05: %d / %d (%.1f%%)\n\n",
sum(pvals < 0.05), length(pvals), 100 * mean(pvals < 0.05)))
cat("Report this proportion: it is the honest measure of how much the result\n")
cat("depends on the placement of the 25 grafted taxa.\n\n")
}
}

phylo_results <- list(signal = signal_results, data = phylo_data, tree = phy_pruned)
cat("Phylogenetic analysis complete.\n\n")
