##############################################################################################
# The block is checking whether results depend on tree choice by comparing Thomas (updated) 
# against three other topologies (Smyth, Pegas, Andres). 
# Thomas has some known polyphyletic placements; the other three are more conservative/cladistic
#################################################################################################

TREE_SPECS <- c(Thomas = "output/phylogenetics/Thomas_updated.nex",
Smyth = "data/phylogenetics/Smyth",
Pegas = "data/phylogenetics/Pegas",
Andres = "data/phylogenetics/Andres_tree.tre")

norm_nm <- function(x) {
x <- gsub("[^A-Za-z0-9_]", "", gsub("[ ]+", "_", trimws(as.character(x))))
tolower(x)
}

prep_tree <- function(path) {
tr <- tryCatch(ape::read.nexus(path), error = function(e) ape::read.tree(path))
if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
tr <- ape::reorder.phylo(tr)
if (!ape::is.binary(tr)) tr <- ape::multi2di(tr)
if (any(tr$edge.length <= 0))
tr$edge.length[tr$edge.length <= 0] <- 1e-6 * max(tr$edge.length)
tr
}

pd_all <- performance_data_clean[is.finite(performance_data_clean$PC1) &
is.finite(performance_data_clean$PC2), ]
pd_all$key <- norm_nm(pd_all$species)

sens_rows <- list()
panels <- list()

for (nm in names(TREE_SPECS)) {
tr <- prep_tree(TREE_SPECS[[nm]])
tip_key <- norm_nm(tr$tip.label)
pd <- pd_all[pd_all$key %in% tip_key, ]
pd <- pd[!duplicated(pd$key), ]
keep <- tr$tip.label[tip_key %in% pd$key]
tr <- ape::drop.tip(tr, setdiff(tr$tip.label, keep))
pd <- pd[match(norm_nm(tr$tip.label), pd$key), ]

s1 <- setNames(pd$PC1, tr$tip.label); s2 <- setNames(pd$PC2, tr$tip.label)
a1 <- phytools::fastAnc(tr, s1); a2 <- phytools::fastAnc(tr, s2)
all1 <- c(s1, a1); all2 <- c(s2, a2)
nt <- ape::Ntip(tr)
child_nm <- ifelse(tr$edge[, 2] <= nt, tr$tip.label[tr$edge[, 2]], as.character(tr$edge[, 2]))
edges_t <- data.frame(x = all1[as.character(tr$edge[, 1])], y = all2[as.character(tr$edge[, 1])],
xend = all1[child_nm], yend = all2[child_nm])
edges_t <- edges_t[stats::complete.cases(edges_t), ]

k1 <- phytools::phylosig(tr, s1, method = "K", test = TRUE, nsim = 500)
k2 <- phytools::phylosig(tr, s2, method = "K", test = TRUE, nsim = 500)
sens_rows[[nm]] <- data.frame(tree = nm, n_tips = nt,
K_PC1 = k1$K, p_PC1 = k1$P, K_PC2 = k2$K, p_PC2 = k2$P)

panels[[nm]] <- ggplot() +
geom_raster(data = grid[!is.na(grid$opt), ], aes(x = PC1, y = PC2, fill = opt), interpolate = TRUE) +
scale_fill_gradientn(colours = c("#ffffff", "#e8f4ea", "#bfe0c6", "#8fc79c", "#54a468", "#1b7837"),
name = "Pareto\noptimality", limits = c(0, 1), na.value = "white") +
{ if (!is.null(wing_perf)) geom_polygon(data = wing_perf, aes(x = PC1, y = PC2, group = grid_id),
fill = "grey55", colour = "grey45", alpha = 0.15, linewidth = 0.18) } +
geom_segment(data = edges_t, aes(x = x, y = y, xend = xend, yend = yend),
colour = "grey35", alpha = 0.45, linewidth = 0.28) +
geom_point(data = pd, aes(x = PC1, y = PC2), shape = 21, fill = "black", colour = "black",
size = 1.7, stroke = 0.35) +
coord_fixed(ratio = 0.75, expand = FALSE,
xlim = range(c(grid$PC1, pd$PC1), na.rm = TRUE),
ylim = range(c(grid$PC2, pd$PC2), na.rm = TRUE)) +
labs(title = nm, subtitle = sprintf("K(PC1)=%.3f p=%.3f | K(PC2)=%.3f p=%.3f", k1$K, k1$P, k2$K, k2$P),
x = ax1, y = ax2) +
theme_bw(base_size = 9) +
theme(plot.title = element_text(face = "bold", size = 13),
plot.subtitle = element_text(size = 7.5, colour = "grey35"),
panel.grid = element_blank(),           legend.position = if (nm == "Andres") "right" else "none")
}

sens_tbl <- do.call(rbind, sens_rows)
dir.create("output/results", showWarnings = FALSE, recursive = TRUE)
write.csv(sens_tbl, "output/results/tree_topology_sensitivity_biomech.csv", row.names = FALSE)
cat("\nBlomberg's K across topologies (biomechanical space):\n"); print(sens_tbl)

combo <- (panels[["Thomas"]] | panels[["Smyth"]]) / (panels[["Pegas"]] | panels[["Andres"]]) +
patchwork::plot_annotation(
title = "Performance phylomorphospace across tree topologies",
subtitle = "Thomas (updated, default) vs Smyth vs Pegas vs Andres - same optimality landscape, different ancestral reconstructions")

ggsave("output/plots_PDF/PERF_15E_tree_topology_comparison.pdf", combo, width = 11, height = 11)
cat("Written: PERF_15E_tree_topology_comparison\n\n")
