results_list <- vector("list", length(outlines_list))
n_skipped <- 0
#################################################################################
# The block is going through each outline and calculate the main flight measures
#################################################################################
for (i in seq_along(outlines_list)) {
specimen_name <- names(outlines_list)[i]
outline <- outlines_list[[i]]
######################################################################
# The block is matching the outline to its row in the main data table
######################################################################
specimen_row <- complete_data[which(complete_data$SPECIES..Pegas. == specimen_name), ]
if (nrow(specimen_row) == 0) { n_skipped <- n_skipped + 1; next }
specimen_row <- specimen_row[1, ]
AR <- if ("WAR" %in% colnames(specimen_row)) as.numeric(specimen_row$WAR) else NA_real_
MASS_KG   <- if ("Mass_kg" %in% colnames(specimen_row)) as.numeric(specimen_row$Mass_kg) else NA_real_
WSPAN_CM  <- if ("Wingspan_cm" %in% colnames(specimen_row)) as.numeric(specimen_row$Wingspan_cm) else NA_real_
r2_hat_val <- calculate_r2_hat(outline)
metrics <- tryCatch({
data.frame(
species = specimen_name, aspect_ratio = AR, r2_hat = r2_hat_val,
#################################################################################################
# The block is vonning_mises_stress is our only stress metric now - stress_index 
# and the old stress_root are gone for good, no fallback or side-by-side comparison
#################################################################################################
von_mises_stress = calculate_von_mises_stress(outline, MASS_KG, WSPAN_CM),
trailing_edge = calculate_trailing_edge(outline),
shape_complexity = calculate_shape_complexity(outline),
aspect_ratio_outline = calculate_aspect_ratio(outline),
second_moment_DIAGNOSTIC = calculate_second_moment_DIAGNOSTIC(outline),
stringsAsFactors = FALSE
)
}, error = function(e) {
cat(paste("Error calculating metrics for", specimen_name, ":", e$message, "\n"))
NULL
})
results_list[[i]] <- metrics
if (i %% 20 == 0) cat(paste("Processed", i, "/", length(outlines_list), "specimens\n"))
}
performance_data <- dplyr::bind_rows(results_list)
cat(paste("\nCalculated performance metrics for", nrow(performance_data), "specimens\n"))
cat(paste("Outlines with no matching database row:", n_skipped, "\n\n"))
if (nrow(performance_data) == 0 || !"species" %in% colnames(performance_data)) {
stop("performance_data is empty or has no species column. Every metric call ",
"failed - check the error lines printed above.")
}
#################################################################################################
# The block is checking that the new measurements look sensible before later scripts use them
#################################################################################################
cat("r2_hat (expected roughly 0.30-0.70; a spread over orders of magnitude\nmeans outlines were not span-aligned):\n")
print(summary(performance_data$r2_hat))
cat("\nvon_mises_stress:\n")
print(summary(performance_data$von_mises_stress))
cat("\ntip_angle (degrees, expect roughly 50-120 per Walters et al.):\n")
print(summary(performance_data$tip_angle))
cat("\npitch_agility (r2_hat-based):\n")
print(summary(performance_data$pitch_agility))
wa <- complete_data$Wing_area_cm2[match(performance_data$species, complete_data$SPECIES..Pegas.)]
ok <- is.finite(performance_data$r2_hat) & is.finite(wa) & wa > 0
if (sum(ok) > 10) {
cat(sprintf("\nCorrelation r2_hat vs real wing area (log-log): %.3f\n",
cor(log(performance_data$r2_hat[ok]), log(wa[ok]))))
cat("Near zero = scale-invariant, as intended.\n")
}
ok2 <- is.finite(performance_data$aspect_ratio) & is.finite(performance_data$aspect_ratio_outline)
if (sum(ok2) > 10) {
cat(sprintf("\nCorrelation database AR vs outline AR: %.3f\n",
cor(performance_data$aspect_ratio[ok2],
performance_data$aspect_ratio_outline[ok2])))
}
cat("\n")
#################################################################################################
# The block is adding the taxonomic and environmental columns needed later in the pipeline
#################################################################################################
join_cols <- c("SPECIES..Pegas.", "Order", "clade", "family", "Environment","Diet.1", "Diet.2", "Midpoint",
"Period.Name", "Depositional.settings.paleoenvironment","Wingspan_cm", "Mass_kg", "Wing_area_cm2")
join_cols <- join_cols[join_cols %in% colnames(complete_data)]
performance_data <- performance_data %>%
left_join(complete_data %>% dplyr::select(all_of(join_cols)) %>% distinct(), by = c("species" = "SPECIES..Pegas."))
