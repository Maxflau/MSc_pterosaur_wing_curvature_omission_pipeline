######################################################################################
# The lines are measuring signal and compare ordinary and phylogenetic trends
######################################################################################
res <- do.call(rbind, lapply(VARS, function(v) {
  y <- d[[v]]
########################################################################################################
# The block is bounding rank ratios and composites on 0-1 scale
# Therefore  their residuals cannot be gaussian near the bounds: a logit puts them on an unbounded scale
########################################################################################################
if (v %in% c("pareto_rank_ratio", "optimality_mean", "optimality_geom")) {
  y <- qlogis(pmin(pmax(y, 0.001), 0.999))
} else if (!v %in% c("PC1", "PC2") && all(y > 0, na.rm = TRUE)) {
  y <- log(y)
}

ok <- is.finite(y) & is.finite(d$Midpoint)
if (sum(ok) < 20) { cat(sprintf("%-20s skipped (n=%d)\n", v, sum(ok))); return(NULL) }

trs <- drop.tip(tr, tr$tip.label[!ok])
lam <- tryCatch(phylosig(trs, setNames(y[ok], trs$tip.label),
                         method = "lambda", test = TRUE), error = function(e) NULL)

df <- data.frame(y = y[ok], Midpoint = d$Midpoint[ok], species = trs$tip.label)
o  <- summary(lm(y ~ Midpoint, data = df))$coefficients
f  <- fit_pgls(df, trs, v)

if (is.null(f$model)) {
  ps <- pse <- pp <- NA_real_
} else {
  tt <- summary(f$model)$tTable
  ps <- tt["Midpoint","Value"]; pse <- tt["Midpoint","Std.Error"]
  pp <- tt["Midpoint","p-value"]
}

cat(sprintf("%-20s n=%3d | OLS p=%.4g | PGLS p=%.4g | lambda=%.3f | %s\n",
            v, sum(ok), o[2,4], pp, f$lambda, f$method))

data.frame(variable = v, n = sum(ok),
           transform = if (v %in% c("pareto_rank_ratio","optimality_mean",
                                    "optimality_geom")) "logit"
           else if (v %in% c("PC1","PC2")) "none" else "log",
           signal_lambda = if (!is.null(lam)) round(lam$lambda, 3) else NA,
           signal_p      = if (!is.null(lam)) signif(lam$P, 4) else NA,
           ols_slope = round(o[2,1], 6), ols_p = signif(o[2,4], 4),
           pgls_slope = round(ps, 6), pgls_se = round(pse, 6),
           pgls_p = signif(pp, 4), pgls_lambda = round(f$lambda, 3),
           pgls_method = f$method, stringsAsFactors = FALSE)
}))

if (is.null(res)) stop("No variable produced a usable model.")
#######################################################################################################
# The block is evaluating the ratio of slopes measures how much of the uncorrected signal was ancestry
######################################################################################################
res$slope_ratio <- round(res$ols_slope / res$pgls_slope, 2)
res$verdict <- ifelse(is.na(res$pgls_p), "PGLS failed",
                      ifelse(res$ols_p < 0.05 & res$pgls_p < 0.05, "trend holds",
                             ifelse(res$ols_p < 0.05, "LOST after correction", "no trend either way")))

cat("\n")
######################################################################################
# The block is showing a fixed number of decimals so every row lines up
######################################################################################
res_display <- res
res_display$signal_lambda <- format_fixed(res$signal_lambda, 3)
res_display$ols_p         <- format_fixed(res$ols_p, 4)
res_display$pgls_p        <- format_fixed(res$pgls_p, 4)
res_display$pgls_lambda   <- format_fixed(res$pgls_lambda, 3)
res_display$slope_ratio   <- format_fixed(res$slope_ratio, 2)
print(res_display[, c("variable","n","transform","signal_lambda","ols_p","pgls_p",
                      "pgls_lambda","slope_ratio","verdict")], row.names = FALSE)


######################################################################################
# The block is saving the result into a csv table in the supplemental folder
######################################################################################
res$metric_set <- SET_NAME
write.csv(res, sprintf("output/results/Supplemental_performance/pgls_%s_metrics.csv", SET_NAME),row.names = FALSE)
cat(sprintf("\nWritten: pgls_%s_metrics.csv\n", SET_NAME))
cat("Report OLS and PGLS together: the gap between them is the result.\n\n")
