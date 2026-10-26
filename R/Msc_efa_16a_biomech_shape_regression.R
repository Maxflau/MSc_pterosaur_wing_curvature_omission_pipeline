#########################################################################################################################
# The block is assessing potential relationship between biomechanical measures to EFA shape positions and save the plots
#########################################################################################################################
if (!exists("sp") || !is.data.frame(sp)) {
if (!exists("shape_perf")) stop("sp/shape_perf not found - source Msc_efa_12_save_helpers.R first.")
sp <- as.data.frame(shape_perf)  # rebuild if missing/stale from a section-only run
}
allo_targets <- c("shapePC1","shapePC2","trailing_edge")
allo_preds   <- BIO_VARS[BIO_VARS %in% names(sp)]

######################################################################################################################
# The block is correlating biomechanical metrics with shape PC positions to evaluate which measures track each axis
######################################################################################################################
bio_mat   <- sp[, c(allo_preds,"shapePC1","shapePC2")]
corr_mat  <- cor(bio_mat, use="pairwise.complete.obs")
save_csv(as.data.frame(corr_mat), "EFA_biomech_shapePC_correlations")
cat("  Correlation (shapePC1):\n")
print(round(corr_mat["shapePC1", allo_preds], 3))
cat("  Correlation (shapePC2):\n")
print(round(corr_mat["shapePC2", allo_preds], 3))

###########################################################################################################################
# The block is fitting linear models of each biomechanical metric against shapePC1 and shapePC2, one predictor at a time
###########################################################################################################################
lm_results <- do.call(rbind, lapply(allo_preds, function(pred) {
do.call(rbind, lapply(c("shapePC1","shapePC2"), function(resp) {
d <- sp[,c(resp,pred)]; d <- d[complete.cases(d),]
if (nrow(d) < 5) return(NULL)
m <- lm(as.formula(paste(resp,"~",pred)), data=d)
s <- summary(m)
data.frame(Response=resp, Predictor=pred,
`R²`=round(s$r.squared,4),`Adj_R²`=round(s$adj.r.squared,4),
F_stat=round(s$fstatistic[1],3), p_value=round(pf(s$fstatistic[1],s$fstatistic[2],
s$fstatistic[3],lower.tail=FALSE),6), check.names = FALSE)
}))
}))
save_csv(lm_results, "EFA_lm_biomech_predictors_shapePC")

###################################################################################################################
# The block is running multiple regression with all biomechanical metrics together against shapePC1 and shapePC2
###################################################################################################################
for (resp in c("shapePC1","shapePC2")) {
d_multi <- sp[, c(resp, allo_preds)] %>% drop_na()
if (nrow(d_multi) < 10) next
m_multi <- lm(as.formula(paste(resp, "~", paste(allo_preds, collapse="+"))), data=d_multi)
s <- summary(m_multi)
coef_df <- as.data.frame(coef(s))
coef_df$term <- rownames(coef_df)
coef_df$Response <- resp
save_csv(coef_df, paste0("EFA_multiple_regression_", resp))
cat(sprintf("  Multiple regression %s: R²=%.3f, adj-R²=%.3f, F-p=%.4f\n", resp, s$r.squared, s$adj.r.squared,
pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)))
}

############################################################################################################################################
# The block is applying a random forest predicting shapePC1/shapePC2 from biomechanical metrics to see which measures matter most per clade
############################################################################################################################################
rf_results   <- list()
pva_panels   <- list()   # collects the two predicted-vs-actual panels to combine below
for (resp in c("shapePC1","shapePC2")) {
d_rf <- sp[, c(resp, allo_preds)] %>% drop_na()
if (nrow(d_rf) < 20) next
rf <- tryCatch(
randomForest(as.formula(paste(resp,"~",paste(allo_preds,collapse="+"))), data=d_rf, ntree=500, importance=TRUE), error=function(e) NULL)
if (!is.null(rf)) {
rf_results[[resp]] <- rf
imp_df <- as.data.frame(importance(rf))
imp_df$Variable <- rownames(imp_df)
imp_df$Response <- resp
save_csv(imp_df, paste0("EFA_rf_importance_", resp))
###############################################################
# The block is creating of the variable importance plot
###############################################################
p_imp <- ggplot(imp_df %>% slice_max(`%IncMSE`, n=10),  aes(reorder(Variable,`%IncMSE`), `%IncMSE`)) +
geom_col(fill="#2980B9", alpha=0.85) + coord_flip() + labs(title=paste("RF variable importance for", resp),
subtitle="Top 10 biomechanical predictors of EFA shape position", x="Biomechanical variable", y="% Increase in MSE") + theme_bw(base_size=12)
save_plot(p_imp, paste0("EFA_rf_importance_", resp), w=9, h=6)

###############################################################
# The block is keeping the predicted-vs-actual panel in memory so the two responses are combined into a single figure instead of two files
############################################################################
d_rf$predicted <- predict(rf)
r2 <- round(cor(d_rf[[resp]], d_rf$predicted, use="complete.obs")^2, 3)
axis_lab <- if (resp == "shapePC1") efa_lab1 else if (resp == "shapePC2") efa_lab2 else resp
pva_panels[[resp]] <- ggplot(d_rf, aes(.data[[resp]], predicted)) + geom_point(aes(colour="Specimen"), alpha=0.6, size=2) +
geom_abline(aes(slope=1, intercept=0, linetype="1:1 (Predicted tendency)"), colour="grey40") +
geom_smooth(aes(linetype="Current fitted trend"), method="lm", se=FALSE, colour="black") + scale_colour_manual(name=NULL,
values=c(Specimen="#2980B9")) +
scale_linetype_manual(name=NULL, values=c("1:1 (Predicted tendency)"="dashed","Current fitted trend"="dotted")) +
labs(title=paste("Random-forest predicted vs. actual:", resp), subtitle=bquote(R^2 == .(r2) ~ "(cross-validated predictions vs. observed values)"),
x=paste("Actual", axis_lab), y=paste("Predicted", axis_lab)) + theme_bw(base_size=12) + theme(legend.position="bottom")
}
}
##########################################################################################################################################
# The block is combining the shapePC1/shapePC2 "predicted vs actual" panels into one figure with a single shared, clearly labelled legend
##########################################################################################################################################
if (length(pva_panels) == 2) {
p_pva_combined <- (pva_panels[["shapePC1"]] + pva_panels[["shapePC2"]]) +  plot_layout(guides = "collect") & theme(legend.position = "bottom")
p_pva_combined <- p_pva_combined + plot_annotation(
title = "EFA shape axes: random-forest predicted vs. actual values",
caption = paste("Points = specimens; dashed line = 1:1 perfect-prediction reference;",
"dotted line = linear trend of predicted on actual.", "Left: shapePC1. Right: shapePC2."))
save_plot(p_pva_combined, "EFA_rf_predicted_vs_actual_combined", w=14, h=7)
} else if (length(pva_panels) == 1) {
save_plot(pva_panels[[1]], paste0("EFA_rf_predicted_vs_actual_", names(pva_panels)[1]), w=8, h=7)
}
