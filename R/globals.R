# Declare global variables to prevent R CMD check from complaining about
# non-standard evaluation used by dplyr, ggplot2, and tidyr.

utils::globalVariables(c(
  "Abundance", "Axis", "Axis_numeric", "cum_abundance", "Group", 
  "mean_Abundance", "Mean_Position", "N_samples", "poly_group", 
  "poly_x", "poly_y", "Q25_Abundance", "Q25_Position", "Q75_Abundance", 
  "Q75_Position", "Sample", "sample_idx", "SE_Abundance", "SE_Max", 
  "SE_Min", "SE_Position", "sum_Abundance", "Taxa_Group", "Taxon", 
  "TotalAbundance", "Value", "width", "x_pos"
))
