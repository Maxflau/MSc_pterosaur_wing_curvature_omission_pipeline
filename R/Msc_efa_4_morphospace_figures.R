###############################################################################
# The block is listing the metric to plot, then save the main morphospace figures
###############################################################################
metrics <- list(list(col="aspect_ratio",label="Aspect Ratio\n(WAR)",log=FALSE, pal="plasma", 
rev=FALSE, brk=c(6,7.5,10,12.5,15),clade_col=TRUE,  transp=FALSE, dark_bg=FALSE),
list(col="r2_hat",label="Second moment\nof area", log=FALSE, pal="viridis", rev=FALSE, brk=NULL,
clade_col=TRUE,  transp=FALSE, dark_bg=FALSE),
list(col="von_mises_stress",  label="Von Mises stress\n(log10)", log=TRUE,  pal="magma",   rev=FALSE, brk=NULL,
clade_col=TRUE,  transp=TRUE,  dark_bg=TRUE))

###############################################################################
# The block is saving each finished figure in the directory
###############################################################################
dir.create("output/plots_PDF", showWarnings=FALSE, recursive=TRUE)
save_plot <- function(p, name, w=11, h=8.5) {
base <- file.path("output/plots_PDF", name)
ggsave(paste0(base,".pdf"), p, width=w, height=h, dpi=300)
cat(sprintf("✓ Plot saved: %s (PDF)\n", name))
}

###################################################################################
# The block is plotting block of the EFA morphospace showing the metric distribution
###################################################################################
metric_labels <- c("1A","1B","1C")
for (i in seq_along(metrics)) {
m    <- metrics[[i]]
lbl  <- metric_labels[i]
p <- make_metric_plot( metric_col = m$col, metric_label = m$label, metric_log = m$log,
metric_pal = m$pal, metric_rev = m$rev, metric_breaks = m$brk, group_var = "clade",
group_cols  = CLADE_COLS, hull_fill  = FALSE, clade_colour = m$clade_col,
transparent_pts = m$transp, dark_bg  = m$dark_bg, legend_bottom = FALSE, perm_factor = "clade",
title  = sprintf("Plot %s: %s with smaller filled wing shapes", lbl, m$label) )
if (is.null(p)) next
save_plot(p, sprintf("CLADE_0%s_%s", lbl, m$col))
}

####################################################################################################
# The block is plotting block of the EFA morphospace showing the distribution of (a)biotic factors
####################################################################################################
make_group_morph <- function(group_var, group_cols, title, subtitle="", hull_alpha=0.18, legend_bottom=FALSE) {
hulls <- make_hulls(shape_perf, group_var)
pts   <- shape_perf %>% mutate(grp = .data[[group_var]]) %>% filter(!is.na(grp))
legend_pos <- if (legend_bottom) "bottom" else "right"
p <- ggplot()

###################################################################################
# The block is adding of the Von Mises stress grid
###################################################################################
if (!is.null(stress_bg)) {
p <- p + geom_raster(data=stress_bg, aes(x,y,fill=z), interpolate=TRUE) + 
scale_fill_stress(palette = "light", name = "Von Mises stress\n(log10)")
}
###################################################################################
# The block is setting up the Impossible region
###################################################################################
if (nrow(impossible_region) > 0)
p <- p + geom_tile(data=impossible_region, aes(x,y), fill="grey80", alpha=0.75)
p <- p + new_scale_fill() + geom_polygon(data=grid_dense, aes(gx,gy,group=grid_id), 
fill="white", colour="grey45", linewidth=0.06, alpha=0.50)
###################################################################################
# The block is creating convex hulls: semi-transparent coloured fill + matching border
###################################################################################
p <- p + geom_polygon(data=hulls,aes(shapePC1, shapePC2, group=grp, fill=grp, colour=grp),
alpha=hull_alpha, linewidth=0.60, show.legend=FALSE) + scale_fill_manual(values=group_cols, guide="none") +
scale_colour_manual(values=group_cols, guide="none")

###################################################################################
# The block is setting up the data points
###################################################################################
p <- p + new_scale_colour() +geom_point(data=pts,aes(shapePC1, shapePC2, colour=grp),
shape=21, fill=NA,size=2.5, stroke=0.80) + scale_colour_manual(values=group_cols, name=gsub("_", " ", group_var),
guide=guide_legend(order=1,ncol=if(legend_bottom) 3L else 1L, override.aes=list(shape=21,fill=NA, size=4.0, stroke=1.0)))
p + labs(title=title, subtitle=subtitle, x=lab1, y=lab2) + coord_fixed(ratio = 0.95, expand = FALSE) + base_theme  
}
cat("\n── Group morphospace plots (CLADE_10–14)\n")
save_plot(make_group_morph("clade", CLADE_COLS,"Morphospace Occupation by Clade","Convex hulls by clade"), 
"CLADE_10_morphospace_clade_convex_hulls")
