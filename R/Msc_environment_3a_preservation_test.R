#########################################################################
# The lines are loading the results table used for the preservation test
#########################################################################
if (exists("performance_data_clean")) {
dat <- performance_data_clean
} else {
dat <- read.csv("output/results/performance_metrics_complete_CLADE.csv", stringsAsFactors = FALSE)
}
###################################################################
# The block is turning depositional labels into preservation groups
###################################################################
env_raw <- trimws(dat$Depositional.settings.paleoenvironment)
dat$Pres_Quality <- dplyr::case_when(
grepl("Lagoon", env_raw, ignore.case = TRUE) ~ "High", grepl("Lacustrine:\\s*small", env_raw, ignore.case = TRUE) ~ "High",
grepl("Lacustrine:\\s*large",  env_raw, ignore.case = TRUE) ~ "Med-High",
grepl("Fluviodeltaic", env_raw, ignore.case = TRUE) ~ "Medium", grepl("Coastal|Shallow marine",  env_raw, ignore.case = TRUE) ~ "Medium",
grepl("Alluvial",   env_raw, ignore.case = TRUE) ~ "Low", grepl("playa",   env_raw, ignore.case = TRUE) ~ "Low",
grepl("Aeolian", env_raw, ignore.case = TRUE) ~ "Low",
TRUE ~ NA_character_   )
############################################################################################################
# The block is undertaking a sanity check reporting any unmatched depositional settings before dropping them
############################################################################################################
unmatched <- unique(env_raw[is.na(dat$Pres_Quality) & env_raw != ""])
if (length(unmatched) > 0) {
warning("Unmatched depositional settings (set to NA): ", paste(unmatched, collapse = " | "))
}

###################################################################################################
# The lines are keeping the depositional groups in a fixed order before plotting and testing them
###################################################################################################
dat$Pres_Quality <- factor(dat$Pres_Quality, levels = c("Low", "Medium", "Med-High", "High"))

############################################################################
# The block is testing conducted on both PC1 and PC2 a on the same specimens
############################################################################
dat <- dat %>% filter(!is.na(Pres_Quality) & !is.na(PC1) & !is.na(PC2))
dat$Pres_Score <- as.numeric(dat$Pres_Quality)   # Low..High coded 1..4
cat("\nSpecimens per preservation-quality category:\n")
print(table(dat$Pres_Quality))
cat("\nDepositional setting x preservation quality:\n")
print(table(trimws(dat$Depositional.settings.paleoenvironment), dat$Pres_Quality))

#####################################################################################################
# The block is applying of Hotelling-Williams test to compares two correlations sharing one variable
#####################################################################################################
hotelling_williams_test <- function(r12, r13, r23, n) {
det_r <- 1 - r12^2 - r13^2 - r23^2 + 2 * r12 * r13 * r23
num <- (r12 - r13) * sqrt((n - 1) * (1 + r23))
den <- sqrt(2 * ((n - 1) / (n - 3)) * det_r + ((r12 + r13)^2 / 4) * (1 - r23)^3)
t_stat <- num / den
df <- n - 3
p_value <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
data.frame(t_stat = round(t_stat, 4), df = df, p_value = signif(p_value, 4))
}
