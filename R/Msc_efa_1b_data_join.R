##################################################################
# The block is checking the diet data has been well read
##################################################################
if (all(c("Diet.1", "Diet.2") %in% names(shape_perf)) && !"Diet_combined" %in% names(shape_perf)) {
d1 <- trimws(as.character(shape_perf$Diet.1)); d2 <- trimws(as.character(shape_perf$Diet.2))
shape_perf$Diet_combined <- ifelse(is.na(d2) | d2 == "" | d2 == d1, d1, paste(d1, d2, sep = "+"))
cat(sprintf("Diet_combined: %d distinct levels, %d with n >= 3\n",length(unique(shape_perf$Diet_combined)),
sum(table(shape_perf$Diet_combined) >= 3)))
print(sort(table(shape_perf$Diet_combined), decreasing = TRUE))
}
###############################################################################
# The block is making the time groups used by the later time-based plots
###############################################################################
if (!"Time_Bin" %in% names(shape_perf)) {
pcol <- intersect(c("Period.Name", "Period_Name"), names(shape_perf))[1]
if (!is.na(pcol)) {
pc <- trimws(as.character(shape_perf[[pcol]]))
pc[pc %in% c("Early Jurassic", "Middle Jurassic")] <- "Early+Middle Jurassic"
shape_perf$Time_Bin <- factor(pc, levels = c("Late Triassic", "Early+Middle Jurassic", "Late Jurassic","Early Cretaceous", "Late Cretaceous"))
cat(sprintf("Time_Bin built: %d specimens across %d bins\n",sum(!is.na(shape_perf$Time_Bin)),nlevels(droplevels(shape_perf$Time_Bin))))
} else {
cat("No period column - Time_Bin not built, Section 10 will skip.\n")
}
}
if (!"Midpoint" %in% names(shape_perf)) {
cat("WARNING: Midpoint absent - Section 10 needs it as well as Time_Bin.\n")
}
##################################################################
# The block is creating a pairwise distance table for later comparison tests If those lines are conducted betadisper phase failed
##################################################################
if (!exists("shape_dist")) {
shape_dist <- dist(as.matrix(shape_perf[, c("shapePC1", "shapePC2")]),method = "euclidean")
}
if (!exists("bio_dist")) {
BIO_FOR_DIST <- intersect(c("aspect_ratio", "r2_hat","von_mises_stress", "trailing_edge", "shape_complexity"),names(shape_perf))
M <- as.matrix(shape_perf[, BIO_FOR_DIST])
##################################################################
# The block is standardising the data
##################################################################
ok <- complete.cases(M) & apply(M > 0, 1, all)
if (sum(ok) < nrow(M)) {
cat(sprintf("bio_dist: %d of %d specimens have all metrics finite and positive.\n",sum(ok), nrow(M)))
cat("Rows with missing metrics are excluded; betadisper needs a matching subset.\n")
}
bio_mat  <- scale(log(M[ok, , drop = FALSE]))
bio_dist <- dist(bio_mat, method = "euclidean")
bio_keep <- ok         
cat(sprintf("bio_dist: %d specimens x %d metrics\n", sum(ok), length(BIO_FOR_DIST)))
}
