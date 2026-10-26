select <- dplyr::select; filter <- dplyr::filter
CLADE_COLS <- c("Anurognathidae" = "#46852F","Azhdarchoidea" = "#DB8E00",
"Basal pterodactyliform" = "#F4A7B9","Most Basal pterosaur" = "#0CE27F",
"Euctenochasmatia" = "#9F1212","Darwinoptera"= "#791444",
"Dsungaripteroidea" = "#60A293","Eudimorphodontoidea"= "#FF00F0",
"Ornithocheiromorpha"  = "#FFE413","Pteranodontia" = "#DA8AAA",
"Rhamphorhynchidae"= "#03240B")

####################################################################
# The block is loading the performance metrics for the EFA figures
####################################################################
PERF_CSV <- "output/results/performance_metrics_complete_CLADE.csv"
if (exists("performance_data_clean")) {
perf <- performance_data_clean
cat(sprintf("perf: using performance_data_clean in memory (%d rows)\n", nrow(perf)))
} else if (file.exists(PERF_CSV)) {
perf <- read.csv(PERF_CSV, stringsAsFactors = FALSE)
cat(sprintf("perf: read from CSV (%d rows) - may predate the metric rewrite\n", nrow(perf)))
} else {
stop("Neither performance_data_clean nor ", PERF_CSV, " is available.")
}

####################################################################################
# The block is checking if the Von Mises Stress are computed neatly for the figures
####################################################################################
if (!"von_mises_stress" %in% names(perf)) {
if (exists("calculate_von_mises_stress") && exists("outlines_list")) {
perf$von_mises_stress <- vapply(seq_len(nrow(perf)), function(i) {
s <- perf$species[i] 
o <- outlines_list[[s]]
if (is.null(o) || nrow(o) < 4) return(NA_real_)
mk <- if ("Mass_kg" %in% colnames(perf)) perf$Mass_kg[i] else NA_real_
ws <- if ("Wingspan_cm" %in% colnames(perf)) perf$Wingspan_cm[i] else NA_real_
tryCatch(calculate_von_mises_stress(o, mk, ws), error = function(e) NA_real_)
}, numeric(1))
cat(sprintf("von_mises_stress recomputed for %d specimens\n",sum(is.finite(perf$von_mises_stress))))
} else {
stop("von_mises_stress absent and cannot be recomputed. Source in the following code in this order in case:\n",
"  <biomechanical functions script> (calculate_von_mises_stress) ->\n",
"  Msc_outline_loading.R -> Msc_performance_metrics_v3.R")
}
}
for (v in c("r2_hat", "pareto_rank_ratio")) {
if (!v %in% names(perf)) {
stop(v, " absent from perf. Re-run ", ifelse(v == "r2_hat", "Msc_performance_metrics_v3.R","Msc_pareto_front_v3.R"), " before this script.")
}
}
perm_csv <- read.csv("output/results/Supplemental_performance/04_permanova_pc1pc2.csv", stringsAsFactors = FALSE, check.names = FALSE)
get_perm <- function(factor_name) {
r <- perm_csv[perm_csv$Factor == factor_name, ]
if (nrow(r) == 0) return(list(r2 = NA, f = NA, p = NA, sig = ""))
##################################################################################
# The block is reading and prepare the data for designing the EFA based figures
##################################################################################
list(r2 = r[["R\u00b2"]][1], f = r$F[1],p = as.numeric(r$p_value[1]), sig = r$sig[1])
}

PERF_JOIN_COLS <- c("species","aspect_ratio","von_mises_stress","trailing_edge","shape_complexity","r2_hat",
"pareto_rank_ratio","on_front","strategy", "PC1","PC2","clade","Order","Diet.1","Diet.2",
"Environment","Depositional.settings.paleoenvironment","Midpoint","Period.Name")

build_shape_perf <- function() {
cols_to_add <- PERF_JOIN_COLS[PERF_JOIN_COLS %in% names(perf)]
if (exists("shape_scores")) {
cols_to_add <- c("species", cols_to_add[!cols_to_add %in% names(shape_scores)])
}
sp <- shape_scores %>% left_join(perf %>% select(any_of(cols_to_add)),
by = "species") %>% select(-ends_with(".y")) %>%
rename_with(~ str_remove(.x, "\\.x$")) %>% filter(!is.na(pareto_rank_ratio))
cat(sprintf("shape_perf: %d specimens | on Pareto front: %d | von_mises_stress: %d\n", nrow(sp), sum(sp$on_front, na.rm = TRUE),
sum(is.finite(sp$von_mises_stress))))
sp
}

if (exists("shape_scores")) {
shape_perf <- build_shape_perf()
} else {
message("shape_scores not found - run build_shape_perf() after MSc_rep.R")
}