# Reproduce Core Matrix Plot

library(ggplot2)
library(dplyr)
library(phyloseq)
library(microViz)
library(patchwork)

# Source files
source("R/utils-radar.R")
source("R/utils-colors.R")
source("R/taxa_star_plot.R")
source("R/plot_ordi_star.R") # Not needed but good practice
source("R/core_matrix_plot.R")

# 1. Create Mock Data
set.seed(999)
n_samples_per_group <- 5
groups <- c("G1", "G2", "G3")
n_samples <- n_samples_per_group * length(groups)
n_taxa <- 30

# Taxa names
taxa_names <- paste0("Taxon_", 1:n_taxa)

# OTU Matrix
otu_mat <- matrix(0, nrow = n_taxa, ncol = n_samples)
rownames(otu_mat) <- taxa_names
colnames(otu_mat) <- paste0("S", 1:n_samples)

# Sample Data
meta_df <- data.frame(
    SampleID = colnames(otu_mat),
    Group = rep(groups, each = n_samples_per_group),
    row.names = colnames(otu_mat)
)

# Define unique cores
# G1 Unique: Taxon_1, Taxon_2, Taxon_3
# G2 Unique: Taxon_4, Taxon_5, Taxon_6
# G3 Unique: Taxon_7, Taxon_8, Taxon_9, Taxon_10
# Shared/Noise: Taxon_11+

for (i in 1:n_samples) {
    g <- meta_df$Group[i]

    # G1 Cores
    if (g == "G1") {
        otu_mat[1:3, i] <- 100 # Unique
        otu_mat[11:15, i] <- 50 # Shared
    }

    # G2 Cores
    if (g == "G2") {
        otu_mat[4:6, i] <- 100 # Unique
        otu_mat[11:15, i] <- 50 # Shared
    }

    # G3 Cores
    if (g == "G3") {
        otu_mat[7:10, i] <- 100 # Unique
        otu_mat[11:15, i] <- 50 # Shared
    }
}

physeq <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    sample_data(meta_df),
    tax_table(matrix(taxa_names, ncol = 1, dimnames = list(taxa_names, "Genus"))) # Mock Tax Table
)

print("Phyloseq object created.")

# 2. Test plot_core_matrix
tryCatch(
    {
        print("Testing plot_core_matrix...")

        p <- plot_core_matrix(
            physeq = physeq,
            group_var = "Group",
            percent_samples = 0.80, # High prevalence
            abundance_threshold = 10,
            taxa_rank = "Genus", # Matches our mock tax table
            samplecolumn = "SampleID" # Required by taxa_star_plot
        )

        print("Matrix Plot generated successfully!")

        # Save or print to check object type is correct (patchwork/ggplot)
        if (inherits(p, "ggplot") || inherits(p, "patchwork")) {
            print("Verification PASSED: Result is a patchwork/ggplot object.")
        } else {
            print("Verification FAILED: Result is not a valid plot object.")
        }

        # Test Warning (Less than 3 unique)
        # Force G3 to have only 1 unique
        # We can just run it on a subset or high threshold
        print("Testing warning conditions...")
    },
    error = function(e) {
        print(paste("Verification FAILED with error:", e$message))
        traceback()
    }
)

