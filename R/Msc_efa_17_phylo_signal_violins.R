#############################################################################################
# The block is measuring phylogenetic signal, savingviolin plots for the taxonomic inclusive group comparison
#############################################################################################
if (exists("pruned_tree")) {
sp$species_clean <- gsub(" ","_",trimws(sp$species))
traits_phylo <- c("shapePC1","shapePC2",BIO_VARS[BIO_VARS %in% names(sp)])
k_results <- do.call(rbind, lapply(traits_phylo, function(tr) {
trait_df <- sp %>% select(species_clean, value=all_of(tr)) %>% filter(!is.na(value))
common   <- intersect(trait_df$species_clean, pruned_tree$tip.label)
if (length(common) < 3)
return(data.frame(trait=tr, n=length(common), K=NA, p_value=NA, sig=""))
vals <- setNames(trait_df$value[match(common,trait_df$species_clean)], common)
sub_tree <- ape::drop.tip(pruned_tree, pruned_tree$tip.label[!pruned_tree$tip.label %in% common])
vals <- vals[sub_tree$tip.label]
res <- tryCatch(phytools::phylosig(sub_tree, vals, method="K", test=TRUE, nsim=999), error=function(e) NULL)
if (is.null(res)) return(data.frame(trait=tr, n=length(common), K=NA, p_value=NA, sig=""))
pv <- round(res$P,6)
data.frame(trait=tr, n=length(common), K=round(res$K,6), p_value=pv,  sig=case_when(pv<0.001~"***",pv<0.01~"**",pv<0.05~"*",TRUE~""))
}))
save_csv(k_results, "EFA_blomberg_K_all_traits")
cat("  Blomberg's K results:\n"); print(k_results)
} else {
cat("  pruned_tree not found - phylogenetic signal skipped.\n")
cat("  Run the phylo section of Msc_perf_surfaces_v4.R first.\n")
}
