####################################################################
# Main driver for conducting the main analyses
####################################################################
rm(list=ls())
##############################################################################
# Line move automatically into the repository folder (if opened in RStudio)
##############################################################################
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
###########################################################
# 1) Preparation of the libraries and biomechanical functions 
##########################################################
source("R/Msc_setup_1_libraries.R")
source("R/Msc_setup_2_biomechanics_functions.R")
#######################################################
# 2) Phylogenetic tree creation and database sorting
#######################################################
source("R/Msc_data_1_tree_creation.R")
source("R/Msc_data_2_specimen_table.R")
###########################################
# 3) Outline extraction and preparation
###########################################
source("R/Msc_outline_1a_extract.R"); source("R/Msc_outline_1b_extract.R")
source("R/Msc_outline_2_load.R");source("R/Msc_outline_3_orientation.R")
###########################################################################################
# 4) Calculation of the flight performance metrics and preparing the theoretical spaces
###########################################################################################
source("R/Msc_performance_1a_metrics.R"); source("R/Msc_performance_1b_metrics.R")
source("R/Msc_performance_2a_morphospace_pca.R"); source("R/Msc_performance_2b_morphospace_pca.R")
source("R/Msc_performance_3_pareto_front.R");source("R/Msc_performance_4_outline_efourier.R")
source("R/Msc_performance_5_theoretical_morphospace.R"); source("R/Msc_performance_6_theoretical_pareto_map.R")
#######################################
# 5) Designing the layout of the figures
#######################################
source("R/Msc_figures_1a_setup.R"); source("R/Msc_figures_1b_setup.R")
source("R/Msc_figures_2_wing_background.R"); source("R/Msc_figures_3_colours_labels.R")
source("R/Msc_figures_4_clade_colours_stress.R");source("R/Msc_figures_5a_biomechanical_panels.R")
source("R/Msc_figures_5b_biomechanical_panels.R");source("R/Msc_figures_6_group_morphospace.R")
####################################################################
# 6) the phylogenetic analyses
####################################################################
source("R/Msc_phylo_1_prep.R"); source("R/Msc_phylo_2a_signal_pgls.R"); 
source("R/Msc_phylo_2b_signal_pgls.R"); source("R/Msc_phylo_3_morphospace_figures.R")
source("R/Msc_phylo_4a_grade_comparison.R"); source("R/Msc_phylo_4b_grade_comparison.R")
source("R/Msc_phylo_5a_tree_figures.R"); source("R/Msc_phylo_5b_tree_figures.R")
source("R/Msc_phylo_6_tree_sensitivity.R")
####################################################################
# 7) the temporal analyses and visualising the figures
####################################################################
source("R/Msc_temporal_1_time_bins.R"); source("R/Msc_temporal_2_disparity_hulls.R")
source("R/Msc_temporal_3_summary_plots.R");source("R/Msc_temporal_4_multipanel_plots.R")
source("R/Msc_temporal_5_plots.R"); source("R/Msc_temporal_6_size_disparity_panels.R")
source("R/Msc_final_1a_stacked_morphospace.R"); source("R/Msc_final_1b_stacked_morphospace.R")
####################################################################
# 8) taphonomic and palaeoenvironmental analyses
####################################################################
source("R/Msc_environment_1_summary.R"); source("R/Msc_environment_2_comparison.R"); 
source("R/Msc_environment_3a_preservation_test.R"); source("R/Msc_environment_3b_preservation_test.R")
######################################################################
# 9) statistical analyses of the performance PCA and the sensitivity analyses
#######################################################################
source("R/Msc_stats_0_setup.R");source("R/Msc_clade_1_save_tables.R");
source("R/Msc_clade_2_summary.R");source("R/Msc_clade_3a_pca_permanova.R"); 
source("R/Msc_clade_3b_pca_permanova.R"); source("R/Msc_stats_1_group_metric_summary.R")
source("R/Msc_stats_2_kruskal_dunn.R"); source("R/Msc_stats_3_disparity_by_group.R");
source("R/Msc_stats_4_permanova.R");source("R/Msc_stats_5_permanova_pc_scores.R")
source("R/Msc_stats_6a_permanova_pca_tables.R"); source("R/Msc_stats_6b_permanova_pca_tables.R")
source("R/Msc_stats_7_pairwise_variance_loadings.R");source("R/Msc_stats_8_allometry_data_prep.R"); 
source("R/Msc_stats_9a_skeletal_allometry.R"); source("R/Msc_stats_9b_skeletal_allometry.R");
source("R/Msc_stats_10_wing_vs_skeleton.R"); source("R/Msc_stats_10b_skeleton_aero_gap.R");
source("R/Msc_stats_11_phylo_vs_nonphylo.R");source("R/Msc_stats_12_curvature_distribution.R"); 
source("R/Msc_stats_13_depositional_data_prep.R");source("R/Msc_stats_14a_taphonomic_sensitivity.R"); 
source("R/Msc_stats_14b_taphonomic_sensitivity.R"); source("R/Msc_stats_15_optimality_through_time.R");
source("R/Msc_stats_16_optimality_bootstrap.R");source("R/Msc_stats_17_categorical_associations.R")
####################################################################
# 10) EFA-based figures
####################################################################
source("R/Msc_efa_0_parse_nts.R");source("R/Msc_efa_1a_data_join.R"); source("R/Msc_efa_1b_data_join.R"); 
source("R/Msc_efa_2_plot_helpers.R");source("R/Msc_efa_3_base_theme.R"); source("R/Msc_efa_4_morphospace_figures.R")
source("R/Msc_efa_5_group_morphospace_figures.R"); source("R/Msc_efa_6_phylo_morphospace.R")
####################################################################
# 11) EFA-based statistical analyses
####################################################################
source("R/Msc_efa_7_grade_comparison_figures.R")
source("R/Msc_efa_8_permanova_shape.R"); source("R/Msc_efa_9_biomechanical_distance.R")
source("R/Msc_efa_10a_distance_tests.R"); source("R/Msc_efa_10b_distance_tests.R") 
source("R/Msc_efa_11a_specimen_table.R"); source("R/Msc_efa_11b_specimen_table.R")
source("R/Msc_efa_12_save_helpers.R"); source("R/Msc_efa_13_loadings_variance_tables.R")
source("R/Msc_efa_14_permanova_all_groups.R"); source("R/Msc_efa_15_environment_curvature_tests.R")
source("R/Msc_efa_16a_biomech_shape_regression.R"); source("R/Msc_efa_16b_biomech_shape_regression.R"); 
source("R/Msc_efa_17_phylo_signal_violins.R");source("R/Msc_efa_18_time_bin_regression.R"); 
source("R/Msc_efa_19_morphospace_expansion.R");source("R/Msc_efa_20_pca_variance_loadings.R"); 
source("R/Msc_efa_21_tree_sensitivity.R")
####################################################################
# 12) wing classification and PGLS tests
####################################################################
source("R/Msc_final_2a_wing_classification.R"); source("R/Msc_final_2b_wing_classification.R")
source("R/Msc_final_3_pgls_tests.R")
