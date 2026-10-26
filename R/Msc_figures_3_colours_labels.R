#############################################################################################################
# The lines are coding the colour of the exclusive clade data point and labels, and folders for later figures
#############################################################################################################

########################################################################
# The block assigning one fixed border colour for each exlcusive clades
########################################################################
CLADE_COLS <- c("Anurognathidae"= "#46852F", "Azhdarchoidea" = "#DB8E00","Basal pterodactyliform" = "#F4A7B9",
"Most Basal pterosaur"= "#0CE27F","Euctenochasmatia" = "#9F1212", "Darwinoptera"= "#791444",
"Dsungaripteroidea" = "#60A293","Eudimorphodontoidea" = "#FF00F0","Ornithocheiromorpha"= "#FFE413", 
"Pteranodontia"= "#DA8AAA","Rhamphorhynchidae"= "#03240B")

########################################################################
# The block is assigning colours for the  Palaeoenvironment categories
########################################################################
ENV_COLS <- c( "Archipelago" = "#A300E9", "Coastal environment"= "#4F4F4F", "Desert" = "#B51963",
"Floodplain" = "#FFC107","Fluvial environment"= "#FF5722", "Forested wetland environment"  = "#FF9800",
"Marine environment" = "#FAAF90", "Semi-arid floodplain"= "#4CAF50")

#################################################################
# The lines are assigning colours for the primary diet categories
#################################################################
DIET_COLS <- c( "Carnivore"= "#E74C3C", "Durophageous" = "#4CAF50", "Filter Feeder"= "#7B1FA2", "Generalist"= "#144B00",
"Herbivore"= "#80DEEA", "Insectivore"= "#FFA726", "Piscivore"= "#F5F5B8")

###############################################################
# The block is assigning for the dietary combination categories
###############################################################
DIETCOMBO_COLS <- c( "Carnivore+Durophageous"= "#9575CD", "Carnivore+Generalist"= "#B0BEC5",
"Carnivore+Insectivore"= "#42A5F5", "Carnivore+Piscivore"= "#7B1FA2","Durophageous+Piscivore"= "#FF8A65", 
"Filter Feeder+Piscivore"  = "#FF7043","Generalist+Piscivore"= "#EC407A", "Herbivore+Durophageous"= "#F48FB1",
"Insectivore+Carnivore"= "#26A69A", "Insectivore+Durophageous" = "#CE93D8","Insectivore+Piscivore"= "#EF5350",
"Piscivore+Carnivore"= "#FFA726", "Piscivore+Durophageous"= "#66BB6A", "Piscivore+Generalist"= "#FFEE58",
"Piscivore+Insectivore"= "#FFC107", "Piscivore+Piscivore"= "#FDD835")

########################################################################
# The block is assigning colours for the depositional setting categories
########################################################################
DEPOSITIONAL_COLS <- c("Aeolian"= "#482072","Alluvial plain"= "#4CAF50","Coastal/Shallow marine"  = "#D55E00",
"Fluviodeltaic"= "#144B00", "Lacustrine: large lake"  ="#FFC107","Lacustrine: small lake"= "#F0CD5D",
"Lagoonal deposit"= "#A300E9","playa"= "#6A9DA3")

########################################################################
# The block is making sure the colour match with the assigned categories
########################################################################
match_fixed_colour <- function(names_in_data, colour_lut) {
norm <- function(x) gsub("\\s+", "", tolower(trimws(x)))
lut <- colour_lut
names(lut) <- norm(names(colour_lut))
out <- unname(lut[norm(names_in_data)])
names(out) <- names_in_data
out
}

#########################################################################################
# The block is generating readable labels for the measurement names used in legends
#########################################################################################
METRIC_LABELS <- c(von_mises_stress = "Von Mises stress\n(scale-invariant)",r2_hat = "Second moment\nof area",
aspect_ratio = "Aspect ratio", trailing_edge = "Trailing edge",shape_complexity = "Shape complexity",
pareto_rank_ratio = "Pareto\nrank ratio")
surf_label <- if (BACKGROUND_VAR %in% names(METRIC_LABELS)) {
METRIC_LABELS[[BACKGROUND_VAR]]
} else BACKGROUND_VAR
cat("Surface variable:", BACKGROUND_VAR, "| colour bar reads:",gsub("\n", " ", surf_label), "\n\n")
