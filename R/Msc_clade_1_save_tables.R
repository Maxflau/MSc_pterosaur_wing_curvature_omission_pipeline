###########################################################################################################
# The lines are saving the main result tables into the Supplemental folder for later statistical analyses
##########################################################################################################
write.csv(performance_data_clean, "output/results/Supplemental_performance/performance_metrics_complete_CLADE.csv", row.names = FALSE)
cat("Well exported: output/results/performance_metrics_complete_CLADE.csv\n")

###############################################################
# The lines are saving the main performance results into a table
###############################################################
write.csv(performance_data_clean, "output/results/performance_metrics_complete_CLADE.csv", row.names = FALSE)
cat("Well exported: output/results/performance_metrics_complete_CLADE.csv\n")

################################################################################
# The block is saving the disparities of the performance metrics into a table
################################################################################
if (exists("disparity_metrics")) {
write.csv(disparity_metrics,"output/results/Supplemental_performance/disparity_metrics.csv", row.names = FALSE)
cat("Well exported: output/results/disparity_metrics.csv\n")
} else {
cat("SKIPPED disparity_metrics - not found (source Msc_temp_3_1.R first).\n")
}
###############################################################
# The block is saving a summary of the time-bin into a table
###############################################################
if (exists("temporal_metrics")) {
write.csv(temporal_metrics, "output/results/Supplemental_performance/temporal_metrics_CLADE.csv", row.names = FALSE)
cat("Well exported: output/results/Supplemental_performance/temporal_metrics_CLADE.csv\n")
} else {
cat("SKIPPED temporal_metrics - not found (source Msc_temp_3_1.R first).\n")
}
#######################################################
# The block is saving a theoretical shape into table
#######################################################
if (exists("theoretical_data")) {
write.csv(theoretical_data, "output/results/Supplemental_performance/theoretical_shape_space.csv", row.names = FALSE)
cat("Well exported: output/results/theoretical_shape_space.csv\n\n")
} else {
cat("SKIPPED theoretical_data - not found (source Msc_theoretical_spaces_v3.R first).\n\n")
}
