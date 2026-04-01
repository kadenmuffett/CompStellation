# Ensure required libraries are loaded
library(phyloseq)
library(dplyr)
library(ggplot2)
library(microViz)
library(patchwork)

# Source all local files from the R/ folder to ensure we capture the latest aesthetic logic
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# Setup directories
out_dir <- "C:/Users/kmmuf/Desktop/micstartest"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Directories ready! Drawing the GlobalPatterns taxa and PCoA tests...")

# 1. GlobalPatterns Setup
data(GlobalPatterns)
gp <- subset_samples(GlobalPatterns, SampleType %in% c("Feces", "Skin", "Tongue", "Freshwater", "Mock"))
gp <- gp %>% tax_fix()
gp <- subset_taxa(gp, taxa_sums(gp) > 0)

# --- taxa_star_plot test ---
taxa_plot <- plot_taxa_star(
    physeq = gp,
    sample_var = "SampleType",
    colors_all = "hclust",
    taxa_rank = "Order",
    samplecolumn = "X.SampleID",
    error_bar = "SE",
    plot_order = "hclust"
)
ggsave(file.path(out_dir, "taxa_star.png"), taxa_plot, width = 8, height = 6, dpi = 300)

# --- plot_ordi_star (separate) test ---
pcoa_plot_sep <- plot_ordi_star(
    physeq = gp,
    sample_var = "SampleType",
    distance = "bray",
    view_type = "separate"
)
ggsave(file.path(out_dir, "pcoa_star_sep.png"), pcoa_plot_sep, width = 10, height = 6, dpi = 300)

# --- plot_ordi_star (together) test ---
pcoa_plot_tog <- plot_ordi_star(
    physeq = gp,
    sample_var = "SampleType",
    distance = "bray",
    view_type = "together"
)
ggsave(file.path(out_dir, "pcoa_star_tog.png"), pcoa_plot_tog, width = 8, height = 6, dpi = 300)

# --- plot_ordi_star_iris test ---
pcoa_iris_plot <- plot_ordi_star_iris(
    physeq = gp,
    sample_var = "SampleType",
    distance = "bray",
    view_type = "separate"
)
ggsave(file.path(out_dir, "pcoa_star_iris.png"), pcoa_iris_plot, width = 10, height = 6, dpi = 300)

message("GlobalPatterns tests complete! Constructing dummy dataset for core tests...")

# 2. Reproduction Dummy Data for Core Tests
set.seed(001)
n_samples_per_group <- 20
groups <- c("Control", "TreatmentA", "TreatmentB")
n_samples <- n_samples_per_group * length(groups)
n_taxa <- 25

taxa_names <- paste0("Taxon_", LETTERS[1:n_taxa])
otu_mat <- matrix(0, nrow = n_taxa, ncol = n_samples)
rownames(otu_mat) <- taxa_names
colnames(otu_mat) <- paste0("Sample", 1:n_samples)
taxa_TABLE <- as.matrix(taxa_names)
rownames(taxa_TABLE) <- taxa_names
taxa_tables <- tax_table(taxa_TABLE)

group_vec <- rep(groups, each = n_samples_per_group)
meta_df <- data.frame(SampleID = colnames(otu_mat), Group = group_vec, row.names = colnames(otu_mat))

for (isValid in 1:n_samples) {
    grp <- group_vec[isValid]
    if (grp == "Control") {
        otu_mat["Taxon_A", isValid] <- 100
        otu_mat["Taxon_M", isValid] <- 10
        otu_mat["Taxon_W", isValid] <- 100
    }
    if (grp %in% c("TreatmentB", "TreatmentA") && runif(1) < 0.4) {
        otu_mat[c("Taxon_A", "Taxon_M"), isValid] <- 30
    }
    if (grp %in% c("Control", "TreatmentA") && runif(1) < 0.8) {
        otu_mat["Taxon_B", isValid] <- 50
        otu_mat["Taxon_N", isValid] <- 100
        otu_mat["Taxon_O", isValid] <- 10
    }
    if (grp == "TreatmentB" && runif(1) < 0.3) {
        otu_mat[c("Taxon_B", "Taxon_N", "Taxon_O"), isValid] <- 30
    }
    if (grp == "TreatmentA") {
        otu_mat["Taxon_C", isValid] <- 200
        otu_mat["Taxon_P", isValid] <- 3
        otu_mat["Taxon_Q", isValid] <- 70
    }
    if (grp == "TreatmentB" && runif(1) < 0.3) {
        otu_mat["Taxon_C", isValid] <- 200
        otu_mat["Taxon_P", isValid] <- 3
        otu_mat["Taxon_Q", isValid] <- 70
    }
    if (grp == "TreatmentB" && runif(1) < 0.3) {
        otu_mat[c("Taxon_C", "Taxon_P", "Taxon_Q"), isValid] <- 20
    }
    if (grp == "TreatmentB") {
        otu_mat["Taxon_D", isValid] <- 500
        otu_mat["Taxon_R", isValid] <- 70
    }
    if (grp == "TreatmentA" && runif(1) < 0.2) {
        otu_mat["Taxon_D", isValid] <- 500
        otu_mat["Taxon_R", isValid] <- 70
    }
    if (runif(1) < 0.2) otu_mat[c("Taxon_E", "Taxon_S", "Taxon_T", "Taxon_U", "Taxon_V"), isValid] <- 30
    if (runif(1) < 0.3) otu_mat["Taxon_F", isValid] <- 10
    if (runif(1) < 0.4) otu_mat["Taxon_G", isValid] <- 10
    if (runif(1) < 0.6) otu_mat["Taxon_H", isValid] <- 10
    if (runif(1) < 0.7) otu_mat["Taxon_I", isValid] <- 10
    if (runif(1) < 0.5) otu_mat["Taxon_J", isValid] <- 10
    if (runif(1) < 0.5) otu_mat["Taxon_K", isValid] <- 10
    if (runif(1) < 0.5) otu_mat["Taxon_L", isValid] <- 10
}
physeq <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE), tax_table(taxa_tables), sample_data(meta_df))

message("Drawing coreselect plot...")
# --- coreselect ---
core_res <- coreselect(physeq = physeq, group_var = "Group", percent_samples = 80, abundance_threshold = 0)
png(file.path(out_dir, "core_select.png"), width = 800, height = 600, res = 100)
print(core_res)
dev.off()

message("Drawing core_matrix plot...")
# --- core_matrix_plot ---
coreplot <- plot_core_matrix(
    physeq = physeq,
    group_var = "Group",
    percent_samples = 0.8,
    abundance_threshold = 0,
    taxa_rank = "OTU", # Setting explicitly
    samplecolumn = "SampleID",
    log_scale = TRUE
)
ggsave(file.path(out_dir, "core_matrix.png"), coreplot, width = 12, height = 6, dpi = 300)

message("Drawing concordant_core plot...")
# --- concordant core ---
concordant_plot <- plot_concordant_core(
    physeq = physeq,
    group_var = "Group",
    percent_samples = 0.8,
    abundance_threshold = 0,
    taxa_rank = "OTU",
    samplecolumn = "SampleID"
)
ggsave(file.path(out_dir, "concordant_core.png"), concordant_plot, width = 8, height = 6, dpi = 300)

message(paste0("All runs completed successfully! Output PNGs saved to: ", out_dir))
