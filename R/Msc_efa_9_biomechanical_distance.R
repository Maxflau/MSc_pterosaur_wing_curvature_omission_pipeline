##########################################################################################################
# The block is generating a distance table from the biomechanical measures through log standardisation
##########################################################################################################
BIO_FOR_DIST <- setdiff(BIO_VARS, "pareto_rank_ratio")
M <- as.matrix(stat_df[, BIO_FOR_DIST])
bio_keep <- complete.cases(M) & apply(M > 0, 1, all)
bio_dist <- dist(scale(log(M[bio_keep, , drop=FALSE])), method="euclidean")
cat(sprintf("bio_dist: %d of %d specimens\n", sum(bio_keep), nrow(M)))

perm_bio_df <- do.call(rbind, Filter(Negate(is.null), lapply(GROUPS, function(g) {
gv <- stat_df[[g]][bio_keep]
ok <- !is.na(gv) & gv != ""
if (length(unique(gv[ok])) < 2) return(NULL)
d <- as.dist(as.matrix(bio_dist)[ok, ok])
tidy_perm(adonis2(d ~ factor(gv[ok]), permutations=9999), g)
})))
write.csv(perm_bio_df, file.path(res_dir,"permanova_biomechanical.csv"), row.names=FALSE)
cat("Written: permanova_biomechanical.csv\n")

###################################################################################################
# The block is assessing whether the spread of the studied differs not only because of dispersion 
# It also allow to check if PERMANOVA are sensitive to dispersion
# It allow to identify if the signal is produced by something else than the placement of the data
####################################################################################################

tidy_betadisp <- function(dist_mat, groups, label) {
ok <- !is.na(groups) & groups != ""
if (length(unique(groups[ok])) < 2) return(NULL)
d <- as.dist(as.matrix(dist_mat)[ok, ok])
bd <- betadisper(d, factor(groups[ok]))
pv <- permutest(bd, permutations=9999)
tbl <- as.data.frame(pv$tab); tbl$Term <- rownames(tbl); tbl$Factor <- label
tbl[!is.na(tbl$F), c("Factor","Term","Df","Sum Sq","Mean Sq","F","Pr(>F)")]
}

disp_shape <- bind_rows(lapply(GROUPS, function(g)
tidy_betadisp(shape_dist, stat_df[[g]], g)))
write.csv(disp_shape, file.path(res_dir,"betadisper_shape_pcs.csv"), row.names=FALSE)

###################################################################################################
# The block is bioing_dist covers only the bio_keep rows, so the grouping vector is subset to match
# passing the full column would fail on dimensions
###################################################################################################
disp_bio <- bind_rows(lapply(GROUPS, function(g)
tidy_betadisp(bio_dist, stat_df[[g]][bio_keep], g)))
write.csv(disp_bio, file.path(res_dir,"betadisper_biomechanical.csv"), row.names=FALSE)
cat("Written: betadisper CSVs\n")

####################################################################################
# The block is applying of Krustal Wallis and Dunn-post Hoc to test each measurement 
# across the main groups and save the results
####################################################################################
bio_available <- BIO_VARS[BIO_VARS %in% names(stat_df)]
bio_available <- bio_available[vapply(bio_available, function(v) sum(is.finite(stat_df[[v]])) >= 10, logical(1))]
if (length(bio_available) == 0) {
stop("No testable metric in stat_df. Check that Msc_performance_metrics_v3.R ran.")
}
cat("Metrics tested:", paste(bio_available, collapse = ", "), "\n")
cat("Factors tested:", paste(GROUPS, collapse = ", "), "\n")

kw_safe <- function(vals, grps) {
ok   <- is.finite(vals) & !is.na(grps) & grps != ""
vals <- vals[ok]; grps <- droplevels(as.factor(as.character(grps[ok])))
keep <- names(which(table(grps) >= 2))
ok2  <- grps %in% keep
vals <- vals[ok2]; grps <- droplevels(grps[ok2])
if (nlevels(grps) < 2 || length(vals) < 4)
return(data.frame(n = length(vals), k = nlevels(grps), H = NA_real_, df_kw = NA_real_, p_value = NA_real_, sig = ""))
res <- tryCatch(kruskal.test(vals, grps), error = function(e) NULL)
if (is.null(res))
return(data.frame(n = length(vals), k = nlevels(grps), H = NA_real_, df_kw = NA_real_, p_value = NA_real_, sig = ""))
pv <- res$p.value
data.frame(n = length(vals), k = nlevels(grps),
H = round(unname(res$statistic), 3),
df_kw = unname(res$parameter),
p_value = signif(pv, 6),
sig = case_when(pv < 0.001 ~ "***", pv < 0.01 ~ "**", pv < 0.05 ~ "*", pv < 0.1 ~ ".", TRUE ~ ""))
}

kw_results <- do.call(rbind, lapply(GROUPS, function(g) {
do.call(rbind, lapply(bio_available, function(v) {
cbind(Group = g, Variable = v, kw_safe(stat_df[[v]], stat_df[[g]]))
}))
}))
###################################################################################################
# The block is testing the five metrics across six factors consist of 42 tests
###################################################################################################
kw_results$p_adjusted_BH <- signif(p.adjust(kw_results$p_value, "BH"), 6)
kw_results$sig_BH <- ifelse(is.na(kw_results$p_adjusted_BH), "", ifelse(kw_results$p_adjusted_BH < 0.05, "*", ""))
write.csv(kw_results, file.path(res_dir, "kruskal_wallis_by_group.csv"), row.names = FALSE)
cat("Written: kruskal_wallis_by_group.csv\n")
sig_rows <- kw_results[!is.na(kw_results$p_adjusted_BH) & kw_results$p_adjusted_BH < 0.05, ]
cat(sprintf("\nSignificant after BH correction: %d of %d tests\n", nrow(sig_rows), nrow(kw_results)))
if (nrow(sig_rows) > 0) print(sig_rows[, c("Group","Variable","n","k","H", "p_value","p_adjusted_BH")],row.names = FALSE)

###################################################################################################
# The block is dunn-posthoccing were applied when the Krustal Wallis test weren't applicable
###################################################################################################
if (requireNamespace("FSA", quietly = TRUE) && nrow(sig_rows) > 0) {
dunn_results <- do.call(rbind, Filter(Negate(is.null), lapply(seq_len(nrow(sig_rows)), function(i) {
g <- sig_rows$Group[i]; v <- sig_rows$Variable[i]
sub <- stat_df[is.finite(stat_df[[v]]) & !is.na(stat_df[[g]]) & stat_df[[g]] != "", ]
keep <- names(which(table(sub[[g]]) >= 3))
sub <- sub[sub[[g]] %in% keep, ]
if (length(keep) < 2) return(NULL)
r <- tryCatch(FSA::dunnTest(reformulate(g, response = v), data = sub, method = "bh")$res,error = function(e) NULL)
if (is.null(r)) return(NULL)
data.frame(Group = g, Variable = v, Comparison = r$Comparison, Z = round(r$Z, 3), p_adj = signif(r$P.adj, 5),
sig = ifelse(r$P.adj < 0.05, "*", ""), stringsAsFactors = FALSE)
})))
if (!is.null(dunn_results)) {
write.csv(dunn_results, file.path(res_dir, "dunn_posthoc.csv"), row.names = FALSE)
cat(sprintf("Written: dunn_posthoc.csv (%d comparisons)\n", nrow(dunn_results)))
}
} else if (nrow(sig_rows) == 0) {
cat("No Kruskal-Wallis survived BH correction - no post-hoc run.\n")
} else {
cat("Package 'FSA' absent - Dunn post-hoc skipped.\n")
}
