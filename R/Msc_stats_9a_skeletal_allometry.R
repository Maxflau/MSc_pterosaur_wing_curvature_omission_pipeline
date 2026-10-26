#############################################################################################
# The lines are testing how each bone measurement changes as overall skeletal size increases
# MIN_FIT with fewer points than this gives an unstable slope
#############################################################################################
if (!exists("POSTCRANIAL")) stop("Source Msc_stats_06a_postcranial.R first.")

MIN_FIT <- 10   # fewer points than this gives an unstable slope

###################################################################################################
# The block is running a log-log regression of one element on skeletal size tested against isometry, 
# Where both variables are lengths so the isometric slope is 1
###################################################################################################
fit_element <- function(el, d, size_var = "SkeletalSize", expected = 1) {
x <- d[[size_var]]; y <- d[[el]]
ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
if (sum(ok) < MIN_FIT) return(NULL)

m  <- lm(log(y[ok]) ~ log(x[ok]))
cf <- summary(m)$coefficients
slope <- cf[2, 1]; se <- cf[2, 2]

#########################################################
# The block is assessing confidence interval on the slope
# In order to directly determine isometry 
#########################################################
ci <- slope + c(-1, 1) * qt(0.975, sum(ok) - 2) * se
t_iso <- (slope - expected) / se
p_iso <- 2 * pt(abs(t_iso), df = sum(ok) - 2, lower.tail = FALSE)

data.frame(element = el, n = sum(ok),
slope = round(slope, 4), se = round(se, 4),
ci_lower = round(ci[1], 4), ci_upper = round(ci[2], 4),
r_squared = round(summary(m)$r.squared, 4),
expected_slope = expected,
p_vs_isometry = signif(p_iso, 4),
stringsAsFactors = FALSE)
}

######################################################################
# The block is grouping allometry of each element
# In order to fit a size relationship per bone across the whole sample
######################################################################
allo <- do.call(rbind, lapply(GM_ELEMENTS, fit_element, d = POSTCRANIAL))
if (is.null(allo)) stop("No element fit succeeded - check SkeletalSize.")

allo$p_adjusted_BH <- signif(p.adjust(allo$p_vs_isometry, method = "BH"), 4)
allo$allometry <- ifelse(allo$p_adjusted_BH >= 0.05, "isometric",
ifelse(allo$slope > allo$expected_slope,
"positive allometry", "negative allometry"))
allo <- allo[order(-allo$slope), ]

cat("Element allometry against the geometric-mean skeletal size:\n")
#########################################################################
# The block is showing a fixed number of decimals so every row lines up
#########################################################################
allo_display <- allo
allo_display$slope         <- format_fixed(allo$slope, 4)
allo_display$ci_lower      <- format_fixed(allo$ci_lower, 4)
allo_display$ci_upper      <- format_fixed(allo$ci_upper, 4)
allo_display$r_squared     <- format_fixed(allo$r_squared, 4)
allo_display$p_adjusted_BH <- format_fixed(allo$p_adjusted_BH, 4)
print(allo_display[, c("element", "n", "slope", "ci_lower", "ci_upper", "r_squared",
"p_adjusted_BH", "allometry")], row.names = FALSE)
write_supp(allo_display, "S12_skeletal_allometry")

cat("\nSlope 1 = the element keeps its proportion as the animal grows.\n")
cat("Above 1 = it grows faster than the skeleton as a whole; below 1, slower.\n")
cat("This is where a disproportionately elongated wing finger would show up.\n\n")

#######################################################################################################
# The block is drawing all element log-log allometry regressions together in one combined panel set 
# Instead of separate files, so slopes are easy to compare
########################################################################################################
for (dd in c("output/plots_PDF")) dir.create(dd, showWarnings = FALSE, recursive = TRUE)

allo_points <- do.call(rbind, lapply(GM_ELEMENTS, function(el) {
x <- POSTCRANIAL$SkeletalSize; y <- POSTCRANIAL[[el]]
ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
if (sum(ok) < MIN_FIT) return(NULL)
data.frame(element = el, log_size = log(x[ok]), log_element = log(y[ok]))
}))

if (!is.null(allo_points) && requireNamespace("ggplot2", quietly = TRUE)) {
allo_lab <- setNames(
sprintf("%s (b=%s, R\u00b2=%s, %s)", allo$element,
format_fixed(allo$slope, 2), format_fixed(allo$r_squared, 2), allo$allometry),
allo$element)
allo_points$panel <- factor(allo_lab[allo_points$element],
levels = allo_lab[allo$element])

p_allo <- ggplot(allo_points, aes(log_size, log_element)) +
geom_point(size = 1.1, alpha = 0.5, colour = "grey30") +
geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
colour = "#b2182b", fill = "#b2182b", alpha = 0.15, linewidth = 0.7) +
geom_abline(slope = 1, intercept = NA, linetype = "dashed", colour = "grey45") +
facet_wrap(~ panel, scales = "free", ncol = 3) +
labs(title = "Skeletal allometry: each element against geometric-mean skeletal size",
subtitle = paste("Points = specimens (log-log). Red line = fitted slope (shaded = 95% CI).",
"Dashed grey = isometric expectation (slope = 1).",
"Facet label: fitted slope (b), R\u00b2, and isometry verdict (S12)."),
x = "log(Skeletal size)", y = "log(Element length)") +
theme_bw(base_size = 10) +
theme(strip.text = element_text(face = "bold", size = 8),
strip.background = element_rect(fill = "grey95", colour = "grey40"),
panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
plot.subtitle = element_text(size = 8, colour = "grey35"))

n_el <- length(unique(allo_points$panel))
ggsave("output/plots_PDF/ALLOMETRY_01_skeletal_elements_combined.pdf", p_allo,
width = 11, height = 3.6 * ceiling(n_el / 3), limitsize = FALSE)
cat("Written: ALLOMETRY_01_skeletal_elements_combined (", n_el, "elements, 1 combined figure)\n\n")
}
