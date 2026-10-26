########################################################################################################################
# The block is classifying the taxe based on the functional trade off limiting them to reach the pareto front
########################################################################################################################
obj_matrix <- as.matrix(performance_data_clean[, names(OBJECTIVES), drop = FALSE])
keep <- complete.cases(obj_matrix) & apply(is.finite(obj_matrix), 1, all)
if (sum(keep) < 10) {
stop("Only ", sum(keep), " specimens have complete objective data.")
}
cls_data <- performance_data_clean[keep, ]
M <- obj_matrix[keep, , drop = FALSE]
for (v in names(OBJECTIVES)) if (OBJECTIVES[[v]] == "min") M[, v] <- -M[, v]

##################################################################################################
# The block is determining how many  specimens can be ranked with the assigned objectives
##################################################################################################
cat(sprintf("Real taxa classified: %d\n\n", nrow(M)))
cls_data$dominance_count <- vapply(seq_len(nrow(M)), function(i) {
ge <- sweep(M, 2, M[i, ], ">="); gt <- sweep(M, 2, M[i, ], ">")
sum(apply(ge, 1, all) & apply(gt, 1, any))
}, numeric(1))
cls_data$on_front <- cls_data$dominance_count == 0
pct <- apply(M, 2, function(x) rank(x, ties.method = "average") / length(x))
colnames(pct) <- names(OBJECTIVES)
cls_data$limiting_objective  <- colnames(pct)[apply(pct, 1, which.min)]
cls_data$limiting_percentile <- apply(pct, 1, min)
cls_data$balance <- apply(pct, 1, function(p) 1 - (max(p) - min(p)))
Ms <- scale(M)
front_s <- Ms[cls_data$on_front, , drop = FALSE]
cls_data$distance_to_front <- vapply(seq_len(nrow(Ms)), function(i) {
min(sqrt(rowSums((sweep(front_s, 2, Ms[i, ], "-"))^2)))
}, numeric(1))
cls_data$strategy <- with(cls_data, ifelse(
on_front, "Best trade-off", ifelse(balance > 0.65, "Middling on all axes", paste0("Held back by ", limiting_objective))))
cat("Strategy classification:\n")
print(table(cls_data$strategy))
cat("\nLimiting objective across all taxa:\n")
print(table(cls_data$limiting_objective))
cat("\nMeaning of each objective:\n")
for (v in names(OBJECTIVES)) {
if (v %in% names(OBJ_MEANING)) cat(sprintf("  %-14s %s\n", v, OBJ_MEANING[[v]]))
}
cat("\n")

##################################################################################
# The block is grouping the Diet and clade data into the wing classification
##################################################################################
for (grp in intersect(c("clade", "Diet.1"), colnames(cls_data))) {
cat(sprintf("Limiting objective by %s:\n", grp))
tb <- table(cls_data[[grp]], cls_data$limiting_objective)
print(tb)
if (nrow(tb) > 1 && ncol(tb) > 1) {
ft <- fisher.test(tb, simulate.p.value = TRUE, B = 10000)
cat(sprintf("Fisher exact (simulated): p = %s\n\n", formatC(ft$p.value, format = "e", digits = 3)))
} else cat("\n")
}

if ("Wingspan_cm" %in% colnames(cls_data)) {
cat("Wingspan by limiting objective (cm):\n")
print(tapply(cls_data$Wingspan_cm, cls_data$limiting_objective, function(x) round(summary(x[is.finite(x)]), 1)))
cat("\n")
}
cat(" These tests treat specimens as independent. With Pagel's lambda near\n")
new_cols <- c("dominance_count", "on_front", "limiting_objective", "limiting_percentile", "balance", "distance_to_front", "strategy")
for (nm in new_cols) {
performance_data_clean[[nm]] <- NA
performance_data_clean[[nm]][keep] <- cls_data[[nm]]
}
##################################################################################################
# The block is saving the wing classification into a csv table in the supplemental folder
##################################################################################################
dir.create("output/results", showWarnings = FALSE, recursive = TRUE)
write.csv(cls_data[, intersect(c("species", "clade", "Diet.1", "Diet.2", "DIETCOMBO_COLS", "Wingspan_cm",
names(OBJECTIVES), new_cols), colnames(cls_data))],
"output/results/wing_classification.csv", row.names = FALSE)
cat("Added to performance_data_clean:", paste(new_cols, collapse = ", "), "\n")
cat("Written: output/results/wing_classification.csv\n\n")
