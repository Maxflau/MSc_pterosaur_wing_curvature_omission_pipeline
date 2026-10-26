###########################################################################
# The lines are loading the data table and auto-detect the key columns
###########################################################################
pterosaur_data <- read.csv("data/pteros_main_data.csv", sep=";", stringsAsFactors = FALSE)
############################################################################
# The block is looking for columns containing numerical data for the analyses
###########################################################################
wing_area_cols <- grep("Wing.area...cm2..EWM..", colnames(pterosaur_data), ignore.case = TRUE, value = TRUE)
war_cols <- grep("war|aspect", colnames(pterosaur_data), ignore.case = TRUE, value = TRUE)
wingspan_cols <- grep("wingspan", colnames(pterosaur_data), ignore.case = TRUE, value = TRUE)
mass_cols <- grep("mass", colnames(pterosaur_data), ignore.case = TRUE, value = TRUE)
cat("Found columns:\n")
cat("Wing area columns:", paste(wing_area_cols, collapse = ", "), "\n")
cat("WAR/Aspect ratio columns:", paste(war_cols, collapse = ", "), "\n")
cat("Wingspan columns:", paste(wingspan_cols, collapse = ", "), "\n")
cat("Mass columns:", paste(mass_cols, collapse = ", "), "\n")
########################################################################
# The block is keeping the first matching column found for each variable
########################################################################
wing_area_col <- if (length(wing_area_cols) > 0) wing_area_cols[1] else NULL
war_col <- if (length(war_cols) > 0) war_cols[length(war_cols)] else NULL
wingspan_col <- if (length(wingspan_cols) > 0) wingspan_cols[1] else NULL
mass_col <- if (length(mass_cols) > 0) mass_cols[1] else NULL
cat("Using columns:\n")
cat("Wing area:", wing_area_col, "\n")
cat("WAR:", war_col, "\n")
cat("Wingspan:", wingspan_col, "\n")
cat("Mass:", mass_col, "\n")
########################################################################
# The block is listing of anatomical measurements used for the analysis
########################################################################
key_vars <- c(
"Wing.length...cm.",
"measured.hind.limbs..cm.",
"Humerus",
"Ulna",
"Mc.IV",
"Wg.ph.1",
"Wg.ph.2",
"Wg.ph.3",
"Wg.ph.4"
)
##################################################################
# The block is adding the columns found automatically above
##################################################################
if (!is.null(wingspan_col)) key_vars <- c(key_vars, wingspan_col)
if (!is.null(wing_area_col)) key_vars <- c(key_vars, wing_area_col)
if (!is.null(war_col)) key_vars <- c(key_vars, war_col)
if (!is.null(mass_col)) key_vars <- c(key_vars, mass_col)
##################################################################
# The lines are keeping only the columns existing within the table
##################################################################
existing_vars <- key_vars[key_vars %in% colnames(pterosaur_data)]
##################################################################
# The block is identifying and context columns to keep
##################################################################
select_cols <- c("SPECIES..Pegas.", "Image_files_names", "Order", "clade", "family", "Environment", 
"Depositional.settings.paleoenvironment", "Diet.1", "Diet.2", "Midpoint", "Period.Name")
morph_data <- pterosaur_data %>% dplyr::select(all_of(select_cols), everything())
########################################################################
# The block is creating simple, fixed column names, easy to reuse later
########################################################################
if (!is.null(wingspan_col)) morph_data$Wingspan_cm <- morph_data[[wingspan_col]]
if (!is.null(wing_area_col)) morph_data$Wing_area_cm2 <- morph_data[[wing_area_col]]
if (!is.null(war_col)) morph_data$WAR <- morph_data[[war_col]]
if (!is.null(mass_col)) morph_data$Mass_kg <- morph_data[[mass_col]]
###########################################################################
# The block is removing extra spaces in columns for better readability in R
########################################################################
morph_data$Environment <- trimws(as.character(morph_data$Environment))
morph_data$Depositional.settings.paleoenvironment <- trimws(as.character(morph_data$Depositional.settings.paleoenvironment))
morph_data$Diet.1 <- trimws(as.character(morph_data$Diet.1))
morph_data$Diet.2 <- trimws(as.character(morph_data$Diet.2))
morph_data$Midpoint <- as.numeric(morph_data$Midpoint)
########################################################################
# The block is keeping only specimens with enough measured data
########################################################################
threshold_na <- length(existing_vars) * 0.7
complete_data <- morph_data %>%
rowwise() %>% mutate(na_count = sum(is.na(c_across(all_of(existing_vars))))) %>%
filter(na_count < threshold_na) %>% dplyr::select(-na_count)
cat(paste("Working with", nrow(complete_data), "specimens after filtering\n\n"))