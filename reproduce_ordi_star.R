# Reproduce PCoA Star Plot

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(phyloseq)
library(microViz)

# Source all R files (since we are not building/installing the package yet)
source("R/utils-radar.R")
source("R/utils-colors.R")
source("R/taxa_star_plot.R")
source("R/plot_ordi_star.R")

# 1. Create Mock Data (Phyloseq Object)
set.seed(42)
n_samples <- 20
n_taxa <- 50

# OTU Table
otu_mat <- matrix(sample(0:100, n_samples * n_taxa, replace = TRUE), nrow = n_taxa)
rownames(otu_mat) <- paste0("Taxa", 1:n_taxa)
colnames(otu_mat) <- paste0("Sample", 1:n_samples)

# Sample Data
meta_df <- data.frame(
  SampleID = paste0("Sample", 1:n_samples),
  Group = sample(c("GroupA", "GroupB", "GroupC"), n_samples, replace = TRUE),
  row.names = paste0("Sample", 1:n_samples)
)

physeq <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  sample_data(meta_df)
)

print("Phyloseq object created.")

# 2. Test plot_ordi_star
tryCatch({
  # Test 1: Together (Default Colors)
  p1 <- plot_ordi_star(
    physeq = physeq,
    sample_var = "Group",
    view_type = "together"
  )
  
  print("PCoA Plot (Together, Default Colors) generated successfully!")
  
  # Test 2: Separate with SE and Low Alpha
  p2_se_alpha <- plot_ordi_star(
    physeq = physeq,
    sample_var = "Group",
    view_type = "separate",
    error_bar = "SE",
    fill_alpha = 0.1
  )
  
  print("PCoA Plot (Separate, SE, Alpha 0.1) generated successfully!")
  
  # Test 3: Separate (Original)
  p2_original <- plot_ordi_star(
    physeq = physeq,
    sample_var = "Group",
    colors_all = c("GroupA" = "red", "GroupB" = "blue", "GroupC" = "green"),
    view_type = "separate"
  )
  
  print("PCoA Plot (Separate) generated successfully!")
  
  if (inherits(p1, "ggplot") && inherits(p2_se_alpha, "ggplot") && inherits(p2_original, "ggplot")) {
    print("Verification PASSED: All results are ggplot objects.")
    
    # Check if p2_original has facets
    if (!is.null(p2_original$facet)) {
       print("Verification PASSED: Separate plot has facets.")
    } else {
       print("Verification FAILED: Separate plot does not have facets.")
    }
    
  } else {
    print("Verification FAILED: Results are not ggplot objects.")
  }
  
}, error = function(e) {
  print(paste("Verification FAILED with error:", e$message))
  traceback()
})

