if (length(missing_time) > 0) {
################################################################################
# The lines are skipping further action here if columns are missing
################################################################################
} else {
###############################################################################
# The block is preparing of the plot
################################################################################
efa_time_plot <- ggplot() +
geom_point(data = shape_perf %>% filter(!is.na(shapePC1)), aes(shapePC1, shapePC2), colour="grey82", size=0.6, alpha=0.55, shape=16) +
geom_polygon(data = efa_hulls_time, aes(shapePC1, shapePC2, group=interaction(Time_label, clade), fill=clade,
colour=clade), alpha=0.28, linewidth=0.65) +
geom_point(data = shape_temporal, aes(shapePC1, shapePC2, fill=clade),  shape=21, colour="black", size=2.2, stroke=0.45) +
################################################################################
# The lines are correcting the time scale for the EFA time disparity analyse
################################################################################
scale_fill_manual(values=efa_time_cols, name="Clade",guide=guide_legend(ncol=2,
override.aes=list(shape=21, colour="black", size=3.5, stroke=0.6))) + scale_colour_manual(values=efa_time_cols, guide="none") +
##############################################################################################
# The block is following lines make sure the column match the references from the main dataset
##############################################################################################
facet_wrap(~ Time_label, ncol=3) +
labs(title   = "EFA Morphospace Through Geological Time", subtitle = sprintf("%d specimens across %d time bins",
nrow(shape_temporal), n_bins_efa), x = lab1, y = lab2) +coord_equal() +
theme_classic(base_size=10) + theme( strip.text = element_text(face="bold", size=8.5), 
strip.background = element_rect(fill="grey95", colour="grey60", linewidth=0.4),
legend.position = "right",legend.text = element_text(size=7.5),
legend.title = element_text(size=9, face="bold"), legend.key.size = unit(0.4, "cm"),
panel.border= element_rect(colour="grey40", fill=NA, linewidth=0.4),
panel.grid.major = element_line(colour="grey90", linewidth=0.2),
panel.spacing = unit(0.5, "cm"), plot.background  = element_rect(fill="white", colour=NA))

################################################################################
# The block is assessing of the shape disparity through time
################################################################################
efa_disparity <- shape_temporal %>% group_by(Time_Bin) %>% summarise(n= n(),n_clades = n_distinct(clade),
SoV_PC1  = var(shapePC1, na.rm=TRUE), SoV_PC2  = var(shapePC2, na.rm=TRUE),
disparity= SoV_PC1 + SoV_PC2, mean_opt = mean(pareto_rank_ratio, na.rm=TRUE),
Midpoint = mean(Midpoint,  na.rm=TRUE), .groups  = "drop") %>% arrange(Midpoint)
disp_time_plot <- ggplot(efa_disparity, aes(x=Midpoint, y=disparity)) +
geom_line(colour="#2980B9", linewidth=0.9) + geom_point(aes(size=n, fill=mean_opt), 
shape=21, colour="grey20", stroke=0.5) +geom_text(aes(label=paste0(as.character(Time_Bin),"\n(n=",n,")")),
vjust=-0.7, size=2.4, colour="grey30") +scale_fill_gradient(low="white", high="#1A5C2A",
name="Mean Pareto\nOptimality", limits=c(0,1)) + scale_size_continuous(name="N specimens", range=c(3,10)) +scale_x_reverse(name="Age (Ma)") +
labs(title = "EFA Morphospace Disparity Through Time", subtitle = "Sum of variances (shapePC1 + shapePC2)",
y = "Disparity (SoV)") + theme_bw(base_size=11)

#######################################################################################################
# The block is saving the summary table of the CLADE_EFA_morphospace_time and  CLADE_EFA_disparity_time
#######################################################################################################
write.csv(efa_disparity,file.path(res_dir, "efa_disparity_through_time.csv"), row.names=FALSE)
cat("efa_disparity_through_time.csv completed\n")
} 
