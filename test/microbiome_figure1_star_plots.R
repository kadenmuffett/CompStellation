# ============================================================
# Figure 1: Microbiome Enterosignatures & Global Microbiome Analysis
#
# Panel order:
#   a) plot_taxa_star  — Human gut enterosignatures by geographic location
#   b) Barplot         — Mean enterosignature composition by geographic location
#   c) plot_ordi_star  — Ordination star plot of global host-associated microbiomes
#   d) hclust tree     — Hierarchical clustering of global host-associated microbiomes
#   e) PCoA scatter    — PCoA ordination of global host-associated microbiomes
#
# Input data:
#   Human panels (a, b): "metadata 1.tsv" + "aligned_model_samples.tsv"
#   Global panels (c, d, e): "x.tsv" (genus-level abundances across biomes)
#
# All figures are saved to the same directory as this script.
# Update `data_dir` below if your input files live elsewhere.
# ============================================================

# --- One-time package installation (run once, then leave commented out) ---
BiocManager::install("kadenmuffett/CompStellation")
# devtools::document("C:/Users/kmmuf/Desktop/microbiomeplots")
# devtools::install("C:/Users/kmmuf/Desktop/microbiomeplots")

# --- Libraries ---
library(CompStellation)
library(phyloseq)
library(microViz)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(vegan)
library(stringr)

# --- Directories ---
# out_dir: where all figures are saved (same folder as this script)
out_dir  <- dirname(rstudioapi::getSourceEditorContext()$path)

# data_dir: where your TSV input files live — update if needed
data_dir <- "."


# ============================================================
# PANELS a & b — Human Gut Enterosignatures
#   Data: "metadata 1.tsv" (geographic_location, Biome, etc.)
#         "aligned_model_samples.tsv" (per-sample enterosignature scores S1–S5)
# ============================================================

# Load metadata and enterosignature scores
metadata_antony <- read.delim(file.path(data_dir, "metadata 1.tsv"))
enter_antony    <- read.delim(file.path(data_dir, "aligned_model_samples.tsv"))

# Merge on internal_id
tenter_antony             <- as.data.frame(enter_antony)
tenter_antony$internal_id <- tenter_antony$X
anthonys_h_df             <- right_join(metadata_antony, tenter_antony, join_by("internal_id"))

# Select: internal_id, gender, geographic_location, and S1–S5 enterosignature columns
anthonys_human_minimal <- anthonys_h_df[, c(1, 2, 19, 43:47)]

# Subset to focal geographic locations (globally representative set)
focal_locations <- c(
  "Fiji", "USA", "China", "Cameroon", "Madagascar",
  "Ethiopia", "Ireland", "United Kingdom", "Indonesia",
  "Singapore", "Tanzania", "Peru"
)
anthonys_human <- subset(anthonys_human_minimal, geographic_location %in% focal_locations)

# Rename enterosignature columns with interpretable labels
anthonys_human$`Esc/Bif`  <- as.numeric(anthonys_human$S1..Esch.Bifi.)
anthonys_human$`Prev`     <- as.numeric(anthonys_human$S2..Prev.)
anthonys_human$`Bact`     <- as.numeric(anthonys_human$S3..Bact.)
anthonys_human$`Firm #1`  <- as.numeric(anthonys_human$S4..Firm.)
anthonys_human$`Firm #2`  <- as.numeric(anthonys_human$S5..Firm2.)

# Trim to working columns
anthonys_human_2 <- anthonys_human[, c(2, 3, 9:13)]
colnames(anthonys_human_2) <- c(
  "internal_id", "geographic_location",
  "Esc/Bif", "Prev", "Bact", "Firm #1", "Firm #2"
)


# ---- Panel a: plot_taxa_star — Human Gut Enterosignatures by Geography ----

fig1a <- plot_taxa_star(
  physeq       = anthonys_human_2,
  taxa_names   = c("Prev", "Bact", "Firm #1", "Esc/Bif", "Firm #2"),
  colors_all   = c("#D55E00","#009E73","#CC79A7", "#0072B2"),
  sample_var   = "geographic_location",
  samplecolumn = "internal_id",
  error_bar    = "SE",
  fill_alpha = 0.8,
  taxa_rank    = "OTU",
  plot_order   = "hclust",
  view_type    = "separate"
) + theme(
  plot.margin  = margin(0, 0, 0, 0, "cm"),
  axis.text.x  = element_text(size = 10, angle = 45),
  strip.text.x = element_text(family = "sans", size = 12, face = "bold")
)

ggsave(
  file.path(out_dir, "fig1a_human_enterosignature_star.png"),
  fig1a, width = 8, height = 6, dpi = 300, units = "in"
)


# ---- Panel b: Barplot — Mean Enterosignature Composition by Geography ----

# Pivot to long format and compute per-location means
anthonys_human_2long <- pivot_longer(
  anthonys_human_2,
  cols      = !c(geographic_location, internal_id),
  names_to  = "Enterosignature",
  values_to = "Percent"
)
anthonys_human_2long_sum <- anthonys_human_2long %>%
  group_by(geographic_location, Enterosignature) %>%
  summarise(Percent = mean(Percent), .groups = "drop")

# Order locations along a diet/industrialisation gradient (least → most Western)
anthonys_human_2long_sum$geographic_location <- factor(
  anthonys_human_2long_sum$geographic_location,
  levels = c(
    "Singapore", "Tanzania", "Cameroon", "Fiji", "Ethiopia",
    "Madagascar", "China", "USA", "United Kingdom", "Ireland",
    "Indonesia", "Peru"
  )
)
anthonys_human_2long_sum <- anthonys_human_2long_sum[
  order(anthonys_human_2long_sum$geographic_location),
]

fig1b <- ggplot(
  anthonys_human_2long_sum,
  aes(Percent, geographic_location, fill = Enterosignature)
) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.margin   = margin(l = -100),
    text            = element_text(face = "bold"),
    axis.text.y     = element_text(size = 12)
  )

ggsave(
  file.path(out_dir, "fig1b_human_enterosignature_barplot.png"),
  fig1b, height = 8, width = 5.5, units = "in", dpi = 300
)


# ============================================================
# PANELS c, d, e — Global Host-Associated Microbiomes
#   Data: "x.tsv" (genus-level abundance table across global biomes)
#         metadata_antony loaded above (includes Biome column)
# ============================================================

# Load genus-level abundance table
anthony_genera <- read.delim(file.path(data_dir, "x.tsv"))

# Build phyloseq sample data from metadata loaded above
metadata_antony_sam           <- sample_data(metadata_antony)
metadata_antony_sam@row.names <- metadata_antony_sam$internal_id

# Format OTU and taxonomy tables
rownames(anthony_genera) <- anthony_genera[, 1]
anthony_genera           <- anthony_genera[, -1]
taxa_TABLE               <- as.matrix(rownames(anthony_genera))
rownames(taxa_TABLE)     <- rownames(anthony_genera)
taxa_tables              <- tax_table(taxa_TABLE)

# Build full phyloseq object
physeq_global <- phyloseq(
  otu_table(anthony_genera, taxa_are_rows = TRUE),
  tax_table(taxa_tables),
  sample_data(metadata_antony_sam)
)

# Subset to host-associated biomes only
host_associated <- c(
  "Birds", "Phylloplane", "Fish", "Root",
  "Human_Nose", "Human_Oral", "Human_Skin", "Human_Vagina",
  "Bovine_Rumen", "Cat_Gut", "Dog_Gut",
  "Human_Gut", "Mouse_Gut", "Pig_Gut", "Porifera"
)
antgen <- subset_samples(physeq_global, Biome %in% host_associated)
antgen <- subset_taxa(antgen, taxa_sums(antgen) > 0)

# Cache ordination and phyloseq object for re-use
antord <- ordinate(antgen, method = "PCoA", distance = "bray")
saveRDS(antgen, file.path(out_dir, "antgen.rds"))

# Further cleaning: remove low-information / ambiguous categories
ant_gen_2 <- subset_samples(antgen,
  Biome != "Pig_Gut" & Biome != "Birds" & Biome != "Phylloplane"
)
plot_ordi_star_tree(
  physeq = ant_gen_3,
  sample_var = "Biome",
  ord = ant_gen_3_ord,  # Pass it here
  n_axes = 5,
  colors_all = c("maroon","purple3","cadetblue","orange3")
) <- subset_taxa(ant_gen_2, taxa_sums(ant_gen_2) > 0)

# Drop Eukaryota-assigned and taxonomically ambiguous entries
keep_taxa <- taxa_names(ant_gen_2)[
  !str_detect(tax_table(ant_gen_2)[, "ta1"], "Eukaryota")
]
ant_gen_2 <- prune_taxa(keep_taxa, ant_gen_2)
ant_gen_2 <- subset_taxa(ant_gen_2, ta1 != ";?;?")
ant_gen_3 <- prune_samples(sample_sums(ant_gen_2) > 0, ant_gen_2)

# Final ordination on cleaned object
ant_gen_3<-microbiome::transform(ant_gen_3, "compositional")
ant_gen_3_ord <- ordinate(ant_gen_3, method = "PCoA", distance = "bray")
saveRDS(ant_gen_3_ord, file.path(out_dir, "antord.rds"))


# ---- Panel c: plot_ordi_star — Global Microbiome Ordination Stars ----

fig1c <- plot_ordi_star(
  physeq     = ant_gen_3,
  ord        = ant_gen_3_ord,
  method     = "PCoA",
  colors_all = c("maroon","purple3","cadetblue","orange3"),
  sample_var = "Biome",
  view_type  = "separate",
  plot_order = "hclust",
  fill_alpha = 0.7
) +  theme(
  plot.margin  = margin(0, 0, 0, 0, "cm"),legend.position = "none",
  axis.text.x  = element_text(size = 10, angle = 45, color = "#595959"),
  strip.text.x = element_text(family = "sans", size = 12, face = "bold")
) + ylab("Mean Axis Position")
fig1c
ggsave(
  file.path(out_dir, "fig1c_global_microbiome_ordi_star.png"),
  fig1c, width = 10, height = 8, dpi = 300, units = "in"
)


# ---- Panel d: hclust Tree — Global Host-Associated Microbiomes ----
# NOTE: vegan::vegdist expects a sample x species matrix.
# If the line below throws an error, replace ant_gen_3_ord with:
#   t(otu_table(ant_gen_3))

d  <- vegan::vegdist(t(otu_table(ant_gen_3)), method = "bray")
hc <- hclust(d)
source("utils-colors.R")
source("ordi_star_hclust_tree.R")

# Example usage with an existing ordination:
braytree<-plot_ordi_star_tree(
  physeq = ant_gen_3,
  sample_var = "Biome",
  ord = ant_gen_3_ord,  # Pass it here
  n_axes = 5,
  colors_all = c("maroon","purple3","cadetblue","orange3")
)

png(
  file.path(out_dir, "fig1d_global_microbiome_hclust.png"),
  width = 6, height = 4, units = "in", res = 300
)
plot_ordi_star_tree(
  physeq = ant_gen_3,
  sample_var = "Biome",
  ord = ant_gen_3_ord,  # Pass it here
  n_axes = 5,
  colors_all = c("maroon","purple3","cadetblue","orange3")
)
dev.off()


# ---- Panel e: PCoA Scatter — Global Host-Associated Microbiomes ----
named_cov<-c("Bovine_Rumen" = "#F2DDE5","Human_Gut"="#E5BACA",
             "Mouse_Gut" = "#D898B0","Cat_Gut" = "#CA7595", "Dog_Gut" = "#B03060",
             "Human_Vagina" = "orange3", "Porifera" = "#A86EDE",
             "Fish" = "#D4B7EE", "Root"="purple3","Human_Oral"="#CADFDF",
             "Human_Nose" = "#94BEC0","Human_Skin"="cadetblue")
fig1e <- plot_ordination(ant_gen_3, ant_gen_3_ord, color = "Biome") +
  scale_color_manual(values = named_cov)+
  geom_point(size = 2.5) +
  theme_minimal() +
  theme(text = element_text(family = "serif")) +
  ggtitle("PCoA Axes 1 & 2")

ggsave(
  file.path(out_dir, "fig1e_global_microbiome_pcoa_scatter.png"),
  fig1e, width = 10, height = 8, dpi = 300, units = "in"
)
