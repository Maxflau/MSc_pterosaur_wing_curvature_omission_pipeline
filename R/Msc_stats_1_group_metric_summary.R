###########################################################################
# The lines are summarising each performance metric within each main group
###########################################################################
if (!exists("write_supp")) stop("Source Msc_stats_0_setup.R first.")

###########################################################################
# The block is calculating simple summary numbers for each group
###########################################################################
summarise_group <- function(d, group_col) {
d <- d[!is.na(d[[group_col]]) & d[[group_col]] != "", ]
levels_present <- unique(d[[group_col]])

out <- do.call(rbind, lapply(levels_present, function(lv) {
s <- d[d[[group_col]] == lv, ]
if (nrow(s) < MIN_N) return(NULL)

row <- data.frame(group = lv, n = nrow(s), stringsAsFactors = FALSE)
for (m in METRICS) {
x <- s[[m]][is.finite(s[[m]])]
row[[paste0("mean_", m)]]   <- if (length(x))     round(mean(x), 5)   else NA
row[[paste0("sd_", m)]]     <- if (length(x) > 1) round(sd(x), 5)     else NA
row[[paste0("median_", m)]] <- if (length(x))     round(median(x), 5) else NA
}
row
}))

if (is.null(out)) return(NULL)
out <- out[order(-out$n), ]
names(out)[1] <- group_col
out
}

####################################################
# The block is adding the combined diet categories
####################################################
SUMMARY_LEVELS <- intersect(c(GROUPS, "Diet_combined"),
colnames(performance_data_clean))

for (g in SUMMARY_LEVELS) {
summ <- summarise_group(performance_data_clean, g)
n_excluded <- length(unique(performance_data_clean[[g]])) -
if (is.null(summ)) 0 else nrow(summ)

cat(sprintf("Summary by %s: %d groups retained, %d excluded (n < %d)\n",
g, if (is.null(summ)) 0 else nrow(summ), n_excluded, MIN_N))
write_supp(summ, paste0("S1_summary_statistics_by_", tolower(g)))
}

###################################################################################################
# The block is sampling sizes showing how many specimens support each later comparison G
# Cosnidering the clade sizes are very unequal the table makes it more readable
###################################################################################################
size_tab <- do.call(rbind, lapply(GROUPS, function(g) {
tb <- table(performance_data_clean[[g]][!is.na(performance_data_clean[[g]])])
if (length(tb) == 0) return(NULL)  # column is all-NA, nothing to tabulate
data.frame(level = g, group = names(tb), n = as.integer(tb),
retained = as.integer(tb) >= MIN_N, stringsAsFactors = FALSE)
}))
size_tab <- size_tab[order(size_tab$level, -size_tab$n), ]

cat("\nGroup sizes:\n")
print(size_tab, row.names = FALSE)
write_supp(size_tab, "S1b_group_sample_sizes")

cat("\n")
