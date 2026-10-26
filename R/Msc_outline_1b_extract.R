
#################################################################################################################
# The block is onlying re-run extraction if the raw images are sitting here locally (too big to track in git).
# Otherwise we just use the pre-extracted outlines already in the .nts file
#################################################################################################################
image_folder <- "./data/pteros_out_finale/"
specimen_names <- complete_data$SPECIES..Pegas.

if (dir.exists(image_folder)) {
extraction_results <- extract_simple_outlines_fuzzy(
image_folder = image_folder,
specimen_names = specimen_names,
output_file = "data/Pterosaur_Outlines_Clean",
threshold = 50,
num_points = 100
)
} else {
cat("Image folder '", image_folder, "' not found - skipping outline extraction ",
"and using the pre-extracted data/Pterosaur_Outlines_Clean.nts instead.\n", sep = "")
}
