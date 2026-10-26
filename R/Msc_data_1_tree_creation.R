############################################################################################################
# The block is building the full tree from the Thomas reference adding missing species via spelling fix, 
#same-genus, reference-tree, or fallback-clade placement
############################################################################################################
library(phytools)
tr <- read.nexus("data/phylogenetics/Thomas")
if (inherits(tr, "multiPhylo")) {
tr <- tr[[1]]                       
attr(tr, "TipLabel") <- NULL       
class(tr) <- "phylo"
}
tr <- reorder(tr)                     
if (any(duplicated(tr$tip.label)))   
stop("Thomas tree itself has duplicated tip labels: ", paste(unique(tr$tip.label[duplicated(tr$tip.label)]), collapse=", "))
if (!is.binary(tr)) tr <- multi2di(tr)

##################################################################
# The block is renaming taxa with different orthograph and genus name
##################################################################
rename_map <- c(
"Gallodactylus_canjuersensis"="Cycnorhamphus_canjuersensis",     
"Hatzegopteryx_type_only"="Hatzegopteryx_thambema",              
"Shenzhoupterus_sanyainus"="Meilifeilong_sanyainus",
"Altmuehlopterus_rhamphastinus"="Germanodactylus_rhamphastinus",
"Jidapterus_edentus_"="Jidapterus_edentus",
"Haopterus_gracilis_"="Haopterus_gracilis",
"Hamipterus_tianshanensis_"="Hamipterus_tianshanensis",
"Plataleorhynchus_streptophorodo"="Plataleorhynchus_streptophorodon",
"Gladocephaloideus_jingangshanen"="Gladocephaloideus_jingangshanensis")   
for (old in names(rename_map)) {
i <- match(old, tr$tip.label); if (!is.na(i)) tr$tip.label[i] <- rename_map[[old]]
}
###############################################################################################
# The block is identifying congeneric and place sister taxa based on (Pegas / Andres) topologies
###############################################################################################
congeneric <- character(0)  
ref_sister <- c( "Simurghia_lamegoi" = "Simurghia_robusta", "Tupuxuara_leonardii" = "Tupuxuara_longicristatus",          
"Thalassodromeus_oberlii"  = "Thalassodromeus_sethi", "Haliskia_peterseni"  = "Ferrodraco_lentoni",                
"Anhanguera_santanae"= "Anhanguera_blittersdorffi", "Anhanguera_robustus"= "Anhanguera_piscator",              
"Liaoxipterus_brachyognathus" = "Istiodactylus_latidens","Eosipterus_yangi" = "Beipiaopterus_chenianus",           
"Kepodactylus_insperatus"= "Ardeadactylus_longicollum", "Tsogtopteryx_mongoliensis" = "Cryodrakon_boreas",            
"Gobiazhdarcho_tsogtbaatari"  = "Nipponopterus_mifunensis","Garudapterus_buffetauti"  = "Lusognathus_almadrava",          
"Longchengopterus_zhaoi" = "Nurhachius_ignaciobritoi", "Eotephradactylus_mcintirae"  = "Seazzadactylus_venieri",        
"Pachagnathus_benitoi" = "Raeticodactylus_filisurensis","Makrodactylus_oligodontus" = "Skiphosoura_bavarica",
"Yixianopterus_jingangshanensis" = "Mimodactylus_libanensis")    

ref_mrca <- list( Gladocephaloideus_jingangshanensis = c("Feilongus_youngi","Moganopterus_zhuiana"),      
Kryptodrakon_progenitor = c("Changchengopterus_pani","Propterodactylus_frankerlae"),   
Balaenognathus_maeuseri = c("Gnathosaurus_subulatus","Plataleorhynchus_streptophorodon"),
Orientognathus_chaoyangensis = c("Angustinaripterus_longicephalus","Dorygnathus_banthensis"),
Puntanipterus_globosus = c("Noripterus_parvus","Noripterus_complicidens"),
Gladocephaloideus_jingangshanensis = c("Feilongus_youngi","Moganopterus_zhuiana")) # andres
ref_mrca[["Gladocephaloideus_jingangshanensis"]] <- NULL

###############################################################################################
# The block is identifying congeneric and place sister taxa based on (Pegas / Andres) topologies
###############################################################################################
clade_anchor <- list( Rhamphinion_jenkinsi  = c("Dimorphodon_macronyx","Caelestiventus_hanseni"),
Cascocauda_rong  = c("Anurognathus_ammoni","Vesperopterylus_lamadongensis"),
Ceoptera_evansae = c("Darwinopterus_modularis","Kunpengopterus_sinensis"),
Archaeoistiodactylus_linglongtaensis = c("Darwinopterus_modularis","Kunpengopterus_sinensis"),
Laueropterus_vitriolus  = c("Skiphosoura_bavarica","Makrodactylus_oligodontus"),
Jianchangopterus_zhaoianus  = c("Darwinopterus_modularis","Kunpengopterus_sinensis"))

################################################################################
# The block is checking if the tree branch were grafted properly
################################################################################
graft_sister <- function(tree, new_tip, ref_tip) {
i <- match(ref_tip, tree$tip.label)
if (is.na(i)) { warning("missing ref tip: ", ref_tip); return(tree) }
el <- tree$edge.length[which(tree$edge[,2]==i)]
bind.tip(tree, new_tip, where=i, position=el/2)
}
graft_mrca <- function(tree, new_tip, anchors) {
anchors <- anchors[anchors %in% tree$tip.label]
if (length(anchors) < 2) { return(tree) }
bind.tip(tree, new_tip, where=getMRCA(tree, anchors))
}
for (nt in names(congeneric)) tr <- graft_sister(tr, nt, congeneric[[nt]])
for (nt in names(ref_sister)) tr <- graft_sister(tr, nt, ref_sister[[nt]])
for (nt in names(ref_mrca))   tr <- graft_mrca(tr, nt, ref_mrca[[nt]])
for (nt in names(clade_anchor)) tr <- graft_mrca(tr, nt, clade_anchor[[nt]])

#########################################################################
# The block is resolving a polyphetic issues affecting Dsungaripteroidea
########################################################################
if (all(c("Germanodactylus_cristatus","Germanodactylus_rhamphastinus") %in% tr$tip.label)) {
tr <- drop.tip(tr, c("Germanodactylus_cristatus","Germanodactylus_rhamphastinus"))
tr <- graft_mrca(tr, "Germanodactylus_cristatus", c("Dsungaripterus_weii","Noripterus_complicidens"))
tr <- graft_sister(tr, "Germanodactylus_rhamphastinus", "Germanodactylus_cristatus")
}

##############################################################
# The block is exporting the updated Thomas phylogenetic tree
#############################################################
added <- c(names(congeneric), names(ref_sister), names(ref_mrca), names(clade_anchor))
cat("Tips after grafting :", Ntip(tr), "\n")
cat("Successfully added  :", sum(added %in% tr$tip.label), "/", length(added), "\n")
fail <- added[!(added %in% tr$tip.label)]
if (length(fail)) cat("FAILED:\n", paste(" -", fail, collapse="\n"), "\n")
class(tr) <- "phylo"                
attr(tr, "TipLabel") <- NULL
write.tree(tr, file="output/phylogenetics/Thomas_updated.tre") 
write.nexus(tr, file="output/phylogenetics/Thomas_updated.nex")
cat("Saved -> Thomas_updated.nex (and .tre)\n")