#################################################################################
# The block is catchessing masses entered in g instead of kg
# checked against an external wingspan-based estimate only corrected
# when /1000 actually lines up with that estimate
#################################################################################
fix_mass_unit_errors <- function(data, mass_col = "Mass_kg", span_col = "Wingspan_cm", ratio_flag = 20, ratio_fix_tol = 5) {
ok <- is.finite(data[[mass_col]]) & data[[mass_col]] > 0 &
is.finite(data[[span_col]]) & data[[span_col]] > 0
wingspan_m <- data[[span_col]] / 100
mass_est   <- 0.2889 * wingspan_m^2.6562
ratio <- data[[mass_col]] / mass_est
flagged <- ok & is.finite(ratio) & ratio > ratio_flag
corrected_mass <- data[[mass_col]]
fixed <- rep(FALSE, nrow(data))
for (i in which(flagged)) {
candidate <- data[[mass_col]][i] / 1000
candidate_ratio <- candidate / mass_est[i]
if (candidate_ratio >= 1 / ratio_fix_tol && candidate_ratio <= ratio_fix_tol) {
corrected_mass[i] <- candidate
fixed[i] <- TRUE
}
}
data[[paste0(mass_col, "_corrected")]] <- corrected_mass
data$mass_est_kg         <- mass_est
data$mass_unit_flagged   <- flagged
data$mass_unit_corrected <- fixed
cat(sprintf("\nMass unit check (external wingspan-based estimator): %d specimens flagged (ratio > %gx)\n",
sum(flagged), ratio_flag))
cat(sprintf("  -> %d corrected (/1000 brought them within %gx of the estimate)\n",
sum(fixed), ratio_fix_tol))
cat(sprintf("  -> %d still flagged, NOT corrected - inspect manually\n",
sum(flagged) - sum(fixed)))
if (sum(flagged) > 0) {
cat("\nFlagged specimens:\n")
print(data.frame(
species = data$species[flagged], Mass_kg_original = round(data[[mass_col]][flagged], 3),
Mass_kg_estimated = round(mass_est[flagged], 3), Mass_kg_corrected = round(corrected_mass[flagged], 3),
Wingspan_cm = data[[span_col]][flagged], corrected = fixed[flagged]
))
}
cat("\n")
data
}
if (all(c("Mass_kg", "Wingspan_cm") %in% colnames(performance_data))) {
performance_data <- fix_mass_unit_errors(performance_data)
} else {
cat("Mass_kg/Wingspan_cm absent - mass unit check skipped.\n")
performance_data$Mass_kg_corrected <- performance_data$Mass_kg
}
cat(paste("\nFinal dataset:", nrow(performance_data), "rows x",
ncol(performance_data), "columns\n\n"))
