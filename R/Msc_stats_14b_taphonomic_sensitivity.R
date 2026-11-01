##############################################################
# The block is running every Kruskal-Wallis test both ways 
# To determine whether the subset agrees with the full sample
##############################################################
TEST_METRICS <- intersect(METRICS, colnames(d))
TEST_GROUPS  <- intersect(c("clade", "Order", "Diet.1", "Diet_combined"),
colnames(d))

kw_p <- function(sub, m, g) {
s <- sub[is.finite(sub[[m]]) & !is.na(sub[[g]]) & sub[[g]] != "", ]
keep_lv <- names(which(table(s[[g]]) >= MIN_N))
s <- s[s[[g]] %in% keep_lv, ]
if (length(keep_lv) < 2 || nrow(s) < 15) return(NA_real_)
kruskal.test(as.formula(paste0("`", m, "` ~ `", g, "`")), data = s)$p.value
}

sens <- do.call(rbind, lapply(TEST_GROUPS, function(g) {
do.call(rbind, lapply(TEST_METRICS, function(m) {
p_all  <- kw_p(d, m, g)
p_well <- kw_p(d[d$well_preserved, ], m, g)
if (is.na(p_all) && is.na(p_well)) return(NULL)
data.frame(group = g, metric = m,
n_all = sum(is.finite(d[[m]])),
n_well = sum(is.finite(d[[m]]) & d$well_preserved),
p_full_sample = signif(p_all, 4),
p_well_preserved = signif(p_well, 4),
stringsAsFactors = FALSE)
}))
}))

if (!is.null(sens)) {
sens$verdict <- with(sens, ifelse(
is.na(p_full_sample) | is.na(p_well_preserved), "not testable",
ifelse(p_full_sample < 0.05 & p_well_preserved < 0.05, "robust",
ifelse(p_full_sample >= 0.05 & p_well_preserved >= 0.05, "null in both",
ifelse(p_full_sample < 0.05, "LOST in well-preserved subset",
"appears only in well-preserved subset")))))

cat("Sensitivity of the group tests to preservation:\n")
print(table(sens$verdict))
fragile <- sens[grepl("LOST|appears only", sens$verdict), ]
if (nrow(fragile) > 0) {
cat("\nResults that change with the subset:\n")
print(fragile[, c("group", "metric", "p_full_sample", "p_well_preserved",
"verdict")], row.names = FALSE)
}
##############################################################################
# The block is showing a fixed number of decimals for interpretation clarity
##############################################################################
sens_display <- sens
sens_display$p_full_sample    <- format_fixed(sens$p_full_sample, 4)
sens_display$p_well_preserved <- format_fixed(sens$p_well_preserved, 4)
write_supp(sens_display, "S40_subset_sensitivity")
}

#########################################################################
# The block is running partial correlations holding completeness constant
# Where a vanishing link suggests preservation bias
###########################################################################
partial_cor <- function(x, y, z) {
ok <- is.finite(x) & is.finite(y) & is.finite(z)
if (sum(ok) < 20) return(c(NA, NA, NA))
rx <- residuals(lm(x[ok] ~ z[ok]))
ry <- residuals(lm(y[ok] ~ z[ok]))
raw <- cor(x[ok], y[ok], method = "spearman")
ct <- cor.test(rx, ry, method = "spearman", exact = FALSE)
c(raw, unname(ct$estimate), ct$p.value)
}

SIZE_REF <- intersect(c("Wingspan_cm", "SkeletalSize"), colnames(d))[1]

pc_tab <- NULL
if (!is.na(SIZE_REF)) {
pc_tab <- do.call(rbind, lapply(TEST_METRICS, function(m) {
r <- partial_cor(d[[m]], d[[SIZE_REF]], d$completeness)
if (all(is.na(r))) return(NULL)
data.frame(metric = m, reference = SIZE_REF,
rho_raw = round(r[1], 4),
rho_partial = round(r[2], 4),
p_partial = signif(r[3], 4),
change = round(r[2] - r[1], 4), stringsAsFactors = FALSE)
}))

if (!is.null(pc_tab)) {
pc_tab$artefact_suspected <- abs(pc_tab$change) > 0.1
cat("\nCorrelation with", SIZE_REF, "before and after controlling for completeness:\n")
#############################################################################
# The block is showing a fixed number of decimals for interpretation clarity
#############################################################################
pc_tab_display <- pc_tab
pc_tab_display$rho_raw     <- format_fixed(pc_tab$rho_raw, 4)
pc_tab_display$rho_partial <- format_fixed(pc_tab$rho_partial, 4)
pc_tab_display$p_partial   <- format_fixed(pc_tab$p_partial, 4)
pc_tab_display$change      <- format_fixed(pc_tab$change, 4)
print(pc_tab_display, row.names = FALSE)
write_supp(pc_tab_display, "S41_partial_correlations")
}
}
