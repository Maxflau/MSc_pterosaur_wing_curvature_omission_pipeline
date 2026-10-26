#####################################################################################################################
# The block is evaluating the trailing edge gap across the exclusive clade (adapted vs predicted )
#####################################################################################################################
if ("trailing_edge" %in% names(sp) && !is.null(rf_results$shapePC2)) {
curv_preds <- setdiff(allo_preds, "trailing_edge")
d_curv <- sp[, c("trailing_edge","clade",curv_preds)] %>% drop_na()
rf_c <- tryCatch( randomForest(trailing_edge ~ ., data=d_curv %>% select(-clade),  ntree=500, importance=TRUE),
error=function(e) NULL)
if (!is.null(rf_c)) {
d_curv$curv_pred <- predict(rf_c)
d_curv$curv_gap  <- d_curv$curv_pred - d_curv$trailing_edge
##################################################################################################
# The block is assessing wing curvature actual vs predicted (same style as the shapePC panels above)
##################################################################################################
r2_curv <- round(cor(d_curv$trailing_edge, d_curv$curv_pred, use="complete.obs")^2, 3)
p_curv_pva <- ggplot(d_curv, aes(trailing_edge, curv_pred)) + geom_point(aes(colour="Specimen"), alpha=0.6, size=2) +
geom_abline(aes(slope=1, intercept=0, linetype="1:1 (perfect prediction)"), colour="grey40") +
geom_smooth(aes(linetype="Fitted trend"), method="lm", se=FALSE, colour="black") + scale_colour_manual(name=NULL, values=c(Specimen="#2980B9")) +
scale_linetype_manual(name=NULL, values=c("1:1 (perfect prediction)"="dashed","Fitted trend"="dotted")) +
labs(title="Wing curvature: predicted vs. actual",
subtitle=bquote(R^2 == .(r2_curv) ~ "(skeleton-based random-forest prediction vs. observed trailing-edge curvature)"),
x="Actual trailing-edge curvature", y="Predicted trailing-edge curvature") + theme_bw(base_size=12) +  theme(legend.position="bottom")
p_gap <- ggplot(d_curv, aes(fct_infreq(as.factor(clade)), curv_gap)) +  geom_boxplot(aes(fill=clade), alpha=0.55, outlier.shape=NA) +
geom_jitter(size=1.2, alpha=0.45, width=0.2) + geom_hline(yintercept=0, linetype="dashed", colour="red") +
labs(title="Curvature gap by clade (EFA space)", subtitle="Positive = skeleton predicts more curvature than observed",
x="Clade", y="Curvature gap (predicted - actual)") + theme_bw(base_size=11) + theme(axis.text.x=element_text(angle=45,
hjust=1), legend.position="none")
##################################################################################################
# The block is combining figure: actual-vs-predicted scatter (left) + gap-by-clade (right)
##################################################################################################
p_curv_combined <- (p_curv_pva + p_gap) + plot_annotation( title = "Wing curvature: predicted-vs-actual fit and gap by clade",
caption = paste("Left: predicted vs. actual trailing-edge curvature, with 1:1 and fitted-trend reference lines.",
"Right: predicted-minus-actual gap per exclusive clade; a positive gap means the skeleton", "predicts more curvature than observed."))
save_plot(p_curv_combined, "EFA_curvature_predicted_vs_actual_combined", w=16, h=7)
gap_summary <- d_curv %>% group_by(clade) %>% summarise(n=n(), mean_gap=round(mean(curv_gap,na.rm=TRUE),5),
sd_gap=round(sd(curv_gap,na.rm=TRUE),5), .groups="drop")
gap_summary$overall_r2 <- r2_curv
save_csv(gap_summary, "EFA_curvature_gap_summary")
}
}

################################################################################################
# The block is building one combined scatter figure with each biomechanical predictor faceted
# showing Pearson r, against shapePC1 and shapePC2
################################################################################################
build_efa_scatter_block <- function(resp, resp_lab) {
rows <- lapply(allo_preds, function(pred) {
d <- sp %>% filter(!is.na(.data[[pred]]), !is.na(.data[[resp]]))
if (nrow(d) < 5) return(NULL)
r_val <- round(cor(d[[pred]], d[[resp]], use="complete.obs"), 3)
data.frame(x=d[[pred]], y=d[[resp]], clade=d$clade, panel=sprintf("%s (r=%s)", pred, r_val))
})
dd <- do.call(rbind, rows)
if (is.null(dd) || nrow(dd) == 0) return(NULL)
ggplot(dd, aes(x, y, colour=clade)) + geom_point(size=1.6, alpha=0.7) + geom_smooth(method="lm", se=TRUE, colour="black", linewidth=0.6, alpha=0.15) +
facet_wrap(~panel, scales="free_x", ncol=3) + labs(title=paste("Biomechanical predictors vs.", resp_lab),
x="Predictor value", y=resp_lab, colour="Exclusive clade") + theme_bw(base_size=9) + theme(strip.text=element_text(face="bold",
size=7.5), legend.position="bottom")
}
p_pc1_block <- build_efa_scatter_block("shapePC1", efa_lab1)
p_pc2_block <- build_efa_scatter_block("shapePC2", efa_lab2)
if (!is.null(p_pc1_block) && !is.null(p_pc2_block)) {
p_scatter_combined <- (p_pc1_block / p_pc2_block) +  plot_layout(guides="collect") &
theme(legend.position="bottom")
p_scatter_combined <- p_scatter_combined + plot_annotation(
title="EFA scatter plots: biomechanical predictors vs. shape axes", caption=paste("Top: each biomechanical predictor against shapePC1.", 
"Bottom: each biomechanical predictor against shapePC2.","Facet labels give the predictor and its Pearson r;",
"line = linear fit with 95% CI, coloured by clade."))
n_pred <- length(allo_preds)
save_plot(p_scatter_combined, "EFA_scatter_biomech_predictors_combined", w=11, h=7*ceiling(n_pred/3))
}
