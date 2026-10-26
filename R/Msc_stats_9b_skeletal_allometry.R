######################################################################################
# The block is assesing the allometric relationship between forelimb against hindlimb 
######################################################################################
FORE <- intersect(c("Humerus", "Ulna", "McIV", "WingPh1", "WingPh2",
"WingPh3", "WingPh4"), GM_ELEMENTS)
HIND <- intersect(c("Femur", "Tibia", "MtIII"), GM_ELEMENTS)

if (length(FORE) >= 3 && length(HIND) >= 2) {
POSTCRANIAL$ForelimbGM <- exp(rowMeans(log(POSTCRANIAL[, FORE]), na.rm = TRUE))
POSTCRANIAL$HindlimbGM <- exp(rowMeans(log(POSTCRANIAL[, HIND]), na.rm = TRUE))

ok <- is.finite(POSTCRANIAL$ForelimbGM) & is.finite(POSTCRANIAL$HindlimbGM) &
POSTCRANIAL$ForelimbGM > 0 & POSTCRANIAL$HindlimbGM > 0

m <- lm(log(POSTCRANIAL$ForelimbGM[ok]) ~ log(POSTCRANIAL$HindlimbGM[ok]))
cf <- summary(m)$coefficients
p_iso <- 2 * pt(abs((cf[2, 1] - 1) / cf[2, 2]), df = sum(ok) - 2, lower.tail = FALSE)

cat(sprintf("Forelimb on hindlimb: slope = %.3f (SE %.3f), R² = %.3f, n = %d\n",
cf[2, 1], cf[2, 2], summary(m)$r.squared, sum(ok)))
cat(sprintf("Against isometry (slope 1): p = %.4g -> %s\n\n", p_iso,
ifelse(p_iso >= 0.05, "isometric",
ifelse(cf[2, 1] > 1, "forelimb grows faster",
"hindlimb grows faster"))))

POSTCRANIAL$ForeHindRatio <- POSTCRANIAL$ForelimbGM / POSTCRANIAL$HindlimbGM
}

#############################################################################################
# The block is allometrying within each clade 
#A global slope can just reflect differences between clades, not real scaling within them 
#############################################################################################
if ("clade" %in% colnames(POSTCRANIAL)) {
d <- POSTCRANIAL[!is.na(POSTCRANIAL$clade), ]

by_clade <- do.call(rbind, lapply(unique(d$clade), function(cl) {
s <- d[d$clade == cl, ]
if (sum(is.finite(s$SkeletalSize)) < MIN_FIT) return(NULL)
res <- do.call(rbind, lapply(GM_ELEMENTS, fit_element, d = s))
if (is.null(res)) return(NULL)
cbind(clade = cl, res)
}))

if (!is.null(by_clade)) {
by_clade$p_adjusted_BH <- signif(p.adjust(by_clade$p_vs_isometry, "BH"), 4)
cat(sprintf("Within-clade allometry: %d clades with n >= %d\n",
length(unique(by_clade$clade)), MIN_FIT))
###########################################################################
# The block is showing a fixed number of decimals so every row lines up
###########################################################################
by_clade$slope         <- format_fixed(by_clade$slope, 4)
by_clade$ci_lower      <- format_fixed(by_clade$ci_lower, 4)
by_clade$ci_upper      <- format_fixed(by_clade$ci_upper, 4)
by_clade$r_squared     <- format_fixed(by_clade$r_squared, 4)
by_clade$p_adjusted_BH <- format_fixed(by_clade$p_adjusted_BH, 4)
write_supp(by_clade, "S13_skeletal_allometry_by_clade")
} else {
cat("No clade reached n >=", MIN_FIT, "for within-clade allometry.\n")
}
}

cat("\nNOTE: these regressions treat specimens as independent. With Pagel's\n")
cat("lambda near 0.7, a PGLS is the correct test for any slope reported.\n\n")
