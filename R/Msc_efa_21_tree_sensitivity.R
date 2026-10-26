#####################################################################################################################
# The block is assessing the sensitivity of the EFA optimality data across the different phylogenetic tree topologies
#####################################################################################################################
TREE_SPECS_EFA <- c(Thomas = "output/phylogenetics/Thomas_updated.nex", Smyth = "data/phylogenetics/Smyth",
Pegas = "data/phylogenetics/Pegas", Andres = "data/phylogenetics/Andres_tree.tre")
norm_nm_efa <- function(x) {
x <- gsub("[^A-Za-z0-9_]", "", gsub("[ ]+", "_", trimws(as.character(x))))
tolower(x)
}
prep_tree_efa <- function(path) {
tr <- tryCatch(ape::read.nexus(path), error = function(e) ape::read.tree(path))
if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
tr <- ape::reorder.phylo(tr)
if (!ape::is.binary(tr)) tr <- ape::multi2di(tr)
if (any(tr$edge.length <= 0))
tr$edge.length[tr$edge.length <= 0] <- 1e-6 * max(tr$edge.length)
tr
}

################################################################################################################################
# The block is setting up the statistical analyses for the EFA optimality data across the different phylogenetic tree topologies
################################################################################################################################
sp_all <- shape_perf[is.finite(shape_perf$shapePC1) & is.finite(shape_perf$shapePC2), ]
sp_all$key <- norm_nm_efa(sp_all$species)
sens_rows_efa <- list()
panels_efa <- list()
for (nm in names(TREE_SPECS_EFA)) {
tr <- prep_tree_efa(TREE_SPECS_EFA[[nm]])
tip_key <- norm_nm_efa(tr$tip.label)
sd <- sp_all[sp_all$key %in% tip_key, ]
sd <- sd[!duplicated(sd$key), ]
keep <- tr$tip.label[tip_key %in% sd$key]
tr <- ape::drop.tip(tr, setdiff(tr$tip.label, keep))
sd <- sd[match(norm_nm_efa(tr$tip.label), sd$key), ]
s1 <- setNames(sd$shapePC1, tr$tip.label); s2 <- setNames(sd$shapePC2, tr$tip.label)
a1 <- phytools::fastAnc(tr, s1); a2 <- phytools::fastAnc(tr, s2)
all1 <- c(s1, a1); all2 <- c(s2, a2)
nt <- ape::Ntip(tr)
child_nm <- ifelse(tr$edge[, 2] <= nt, tr$tip.label[tr$edge[, 2]], as.character(tr$edge[, 2]))
edges_efa <- data.frame(x = all1[as.character(tr$edge[, 1])], y = all2[as.character(tr$edge[, 1])],
xend = all1[child_nm], yend = all2[child_nm])
edges_efa <- edges_efa[stats::complete.cases(edges_efa), ]

###########################################################################################################################
# The block is preparing the figure structure for the EFA optimality data across the different phylogenetic tree topologies
###########################################################################################################################
k1 <- phytools::phylosig(tr, s1, method = "K", test = TRUE, nsim = 500)
k2 <- phytools::phylosig(tr, s2, method = "K", test = TRUE, nsim = 500)
sens_rows_efa[[nm]] <- data.frame(tree = nm, n_tips = nt, K_shapePC1 = k1$K, p_shapePC1 = k1$P, K_shapePC2 = k2$K, p_shapePC2 = k2$P)
panels_efa[[nm]] <- ggplot() + geom_raster(data = pareto_bg, aes(x = x, y = y, fill = z), interpolate = TRUE) +
scale_fill_gradientn(colours = c("#ffffff", "#e8f4ea", "#bfe0c6", "#8fc79c", "#54a468", "#1b7837"), name = "Pareto\noptimality", na.value = "white") +
{ if (nrow(impossible_region) > 0) geom_tile(data = impossible_region, aes(x = x, y = y),
fill = "grey80", alpha = 0.75, inherit.aes = FALSE) } + geom_polygon(data = grid_dense, aes(x = gx, y = gy, group = grid_id),
fill = "grey88", colour = "grey55", alpha = 0.38, linewidth = 0.06) + geom_segment(data = edges_efa, aes(x = x, y = y, xend = xend, yend = yend),
colour = "grey35", alpha = 0.45, linewidth = 0.28) + geom_point(data = sd, aes(x = shapePC1, y = shapePC2),
shape = 21, fill = "black",colour = "black", size = 1.7, stroke = 0.35) +
coord_fixed(expand = FALSE) + labs(title = nm, subtitle = sprintf("K(PC1)=%.3f p=%.3f | K(PC2)=%.3f p=%.3f", k1$K, k1$P, k2$K, k2$P),
x = lab1, y = lab2) + theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 13),
plot.subtitle = element_text(size = 7.5, colour = "grey35"), panel.grid = element_blank(), legend.position = if (nm == "Andres") "right" else "none")
}

###############################################################################################################################
# The block is saving the EFA optimality statistical signals across the phylomorphospace axes and showing them into the console
###############################################################################################################################
sens_tbl_efa <- do.call(rbind, sens_rows_efa)
dir.create("output/results", showWarnings = FALSE, recursive = TRUE)
write.csv(sens_tbl_efa, "output/results/tree_topology_sensitivity_efa.csv", row.names = FALSE)
cat("\nBlomberg's K across topologies (EFA shape space):\n"); print(sens_tbl_efa)

################################################################
# The block is saving the EFA optimality phylomorphospace figure
################################################################
combo_efa <- (panels_efa[["Thomas"]] | panels_efa[["Smyth"]]) / (panels_efa[["Pegas"]] | panels_efa[["Andres"]]) +
patchwork::plot_annotation( title = "EFA shape phylomorphospace across tree topologies",
subtitle = "Thomas (updated) vs Smyth vs Pegas vs Andres - same optimality landscape, different ancestral reconstructions")
ggsave("output/plots_PDF/CLADE_15E_tree_topology_comparison_EFA.pdf", combo_efa, width = 11, height = 11)
cat("Written: CLADE_15E_tree_topology_comparison_EFA\n\n")
