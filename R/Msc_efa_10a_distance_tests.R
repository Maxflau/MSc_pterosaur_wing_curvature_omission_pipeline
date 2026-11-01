################################################################
# The block is ensuring the directory exist for the analysis
################################################################
if (!exists("res_dir")) {
res_dir <- "output/results/Supplemental_EFA"; dir.create(res_dir, showWarnings=FALSE, recursive=TRUE)
}
required_time_cols <- c("Time_Bin", "Midpoint")
missing_time <- required_time_cols[!required_time_cols %in% names(shape_perf)]

if (length(missing_time) > 0) {
cat(sprintf("  Missing: %s - skipping Section 10.\n",paste(missing_time, collapse=", ")))
} else {
########################################################################################################################
# The block is preparing the data for the temporal disparity analysis by keeping only the rows needed for the figures
########################################################################################################################
bin_order <- c("Late Triassic", "Early+Middle Jurassic", "Late Jurassic","Early Cretaceous", "Late Cretaceous")
bin_order <- bin_order[bin_order %in% unique(shape_perf$Time_Bin)]
shape_temporal <- shape_perf %>% filter(!is.na(Time_Bin), !is.na(shapePC1), !is.na(shapePC2),
!is.na(clade)) %>% mutate(Time_Bin = factor(Time_Bin, levels=bin_order))
n_bins_efa <- nlevels(droplevels(shape_temporal$Time_Bin))
cat(sprintf("  %d specimens × %d bins: %s\n",nrow(shape_temporal), n_bins_efa, paste(levels(droplevels(shape_temporal$Time_Bin)), collapse=" | ")))

#####################################################################################
# The block is drawing a hull for each exclusive clades inside each time bin
###################################################################################
efa_hulls_time <- shape_temporal %>% group_by(Time_Bin, clade) %>% filter(n() >= 3) %>%
do({
idx <- chull(.$shapePC1, .$shapePC2); idx <- c(idx, idx[1])
data.frame(shapePC1=.$shapePC1[idx], shapePC2=.$shapePC2[idx])
}) %>% ungroup()

##########################################################################################
# The block is assigning a colour for all the exclusive clade presented in the analysis
##########################################################################################
clades_in_time <- sort(unique(shape_temporal$clade))
efa_time_cols  <- CLADE_COLS[clades_in_time]
missing_idx    <- is.na(efa_time_cols)
if (any(missing_idx)) {
extra <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(sum(missing_idx))
efa_time_cols[missing_idx] <- extra
}

#################################################################################################
# The block is building clear stratigraphic panel labels indicating the age and sample counts
#################################################################################################
bin_labels <- shape_temporal %>% group_by(Time_Bin) %>% summarise(mid  = round(mean(Midpoint, na.rm=TRUE), 0),
n    = n(), ncl  = n_distinct(clade), .groups="drop") %>%
mutate(label = paste0(as.character(Time_Bin),"\n(~", mid, " Ma, n=", n, ", clades=", ncl, ")"))
shape_temporal <- shape_temporal %>% left_join(bin_labels %>% select(Time_Bin, label), 
by="Time_Bin") %>% mutate(Time_label = factor(label, levels=bin_labels$label))
efa_hulls_time <- efa_hulls_time %>% left_join(bin_labels %>% select(Time_Bin, label), 
by="Time_Bin") %>%mutate(Time_label = factor(label, levels=bin_labels$label))
} 
