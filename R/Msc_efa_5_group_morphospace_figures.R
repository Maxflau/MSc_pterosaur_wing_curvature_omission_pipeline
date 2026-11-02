#######################################################################################
# The block is saving the remaining group of morphospace figures for the abiotic plots
#######################################################################################
save_plot(make_group_morph("Environment", c("Archipelago"="#A300E9","Coastal environment"="#4F4F4F",
"Desert"="#B51963","Floodplain"="#FFC107", "Fluvial environment"="#FF5722","Forested wetland environment"="#FF9800",
"Marine environment"="#FAAF90","Semi-arid floodplain"="#4CAF50"), "Morphospace Occupation by Environment"),
"CLADE_11_morphospace_environment_convex_hulls")

save_plot(make_group_morph("Diet.1",c("Carnivore"="#E74C3C","Durophageous"="#4CAF50","Filter Feeder"="#7B1FA2",
"Generalist"="#144B00","Herbivore"="#80DEEA","Insectivore"="#FFA726","Piscivore"="#F5F5B8"),
"Morphospace Occupation by Primary Diet"),"CLADE_12_morphospace_diet_convex_hulls")

save_plot(make_group_morph("Diet_combined",
# The pair order is now alphabetised when Diet_combined is built (see Msc_efa_1b_data_join.R),
# so the same guild is never split across two reversed-order keys here (e.g. "Durophageous+Herbivore"
# replaces the old "Herbivore+Durophageous", and "Durophageous+Insectivore" replaces "Insectivore+Durophageous")
c("Carnivore+Durophageous"="#9575CD","Carnivore+Generalist"="#B0BEC5","Carnivore+Insectivore"="#42A5F5",
"Carnivore+Piscivore" ="#7B1FA2","Durophageous+Herbivore"="#F48FB1","Durophageous+Insectivore"="#CE93D8",
"Durophageous+Piscivore"="#FF8A65","Filter Feeder+Piscivore"="#FF7043","Generalist+Piscivore"="#EC407A",
"Insectivore+Piscivore"="#EF5350"),"Morphospace by Diet Combination", "Convex hulls by combined diet") +

################################################################
# The block is creating the labels for the dietary combinations
################################################################
guides(fill   = guide_none(),colour = guide_legend(title="Diet (combined)", ncol=1,override.aes=list(shape=21, size=4.5, fill=NA, stroke=1.0))) +
theme(legend.text  = element_text(size=8.5),legend.title = element_text(size=11, face="bold"),
plot.margin  = unit(c(8, 190, 8, 8), "pt")),"CLADE_13_morphospace_dietcombo_convex_hulls", w=14)
save_plot( make_group_morph("Depositional.settings.paleoenvironment", c("Aeolian"="#482072","Alluvial plain"="#4CAF50",
"Coastal/Shallow marine"="#D55E00","Fluviodeltaic"="#144B00","Lacustrine: large lake"="#FFC107",
"Lacustrine: small lake"="#F0CD5D","Lagoonal deposit"="#A300E9","playa"="#6A9DA3"),
"Morphospace Occupation by Depositional Setting"),"CLADE_14_morphospace_depositional_convex_hulls")

cat("\nAll morphospace figures saved.\n")
