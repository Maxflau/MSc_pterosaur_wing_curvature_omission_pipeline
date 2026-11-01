#######################################################################################
# The block is checking whether preservation quality could be driving the main patterns
#######################################################################################
if (!exists("POSTCRANIAL")) stop("Source Msc_stats_06a_postcranial.R first.")
if (!"pareto_rank_ratio" %in% colnames(performance_data_clean)) {
stop("pareto_rank_ratio not found - run Msc_pareto_front_v3.R first.")
}

d <- POSTCRANIAL
d$completeness <- d$n_elements / length(GM_ELEMENTS)
d <- d[is.finite(d$completeness), ]

#########################################################################################################
# The block is stating the arbitrary completeness threshold and reporting its effect on sample size
# It keeps specimens with at least four fifths of elements measured at 0.8
#########################################################################################################
COMPLETENESS_CUT <- 0.8
d$well_preserved <- d$completeness >= COMPLETENESS_CUT

cat(sprintf("Completeness threshold: %.2f\n", COMPLETENESS_CUT))
cat(sprintf("Well preserved: %d | Poorly preserved: %d\n\n",
sum(d$well_preserved), sum(!d$well_preserved)))

if (sum(d$well_preserved) < 30) {
cat("WARNING: fewer than 30 well-preserved specimens. Lower COMPLETENESS_CUT\n")
cat("or treat the subset comparisons below as indicative only.\n\n")
}

################################################################################
# The block is testing whether preservation predicts front membership,
# The results would indicate a taphonomic bias in the Pareto result
################################################################################
front_rows <- NULL
if ("on_front" %in% colnames(d)) {
ok <- is.finite(d$completeness) & !is.na(d$on_front)
tb <- table(d$well_preserved[ok], d$on_front[ok])

if (all(dim(tb) == c(2, 2))) {
ft <- fisher.test(tb)
glmfit <- glm(on_front ~ completeness, data = d[ok, ], family = binomial)
cf <- summary(glmfit)$coefficients

front_rows <- data.frame(
test = c("Fisher: well-preserved vs on-front",
"Logistic: on_front ~ completeness"),
statistic = c(round(ft$estimate, 4), round(cf[2, 1], 4)),
p_value = signif(c(ft$p.value, cf[2, 4]), 4),
stringsAsFactors = FALSE)

cat("Preservation and Pareto front membership:\n")
print(tb)
cat(sprintf("Fisher odds ratio = %.3f, p = %.4g\n", ft$estimate, ft$p.value))
cat(sprintf("Logistic slope on completeness = %+.4f, p = %.4g\n\n",
cf[2, 1], cf[2, 4]))

#################################################################
# The block is showing a fixed number of decimals for readability
#################################################################
front_rows_display <- front_rows
front_rows_display$statistic <- format_fixed(front_rows$statistic, 4)
front_rows_display$p_value   <- format_fixed(front_rows$p_value, 4)
write_supp(front_rows_display, "S39_preservation_predicts_front")
}
}
