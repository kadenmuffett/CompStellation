# CompStellation

<!-- badges: start -->
[![R-CMD-check](https://github.com/kadenmuffett/CompStellation/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kadenmuffett/CompStellation/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `CompStellation` package provides unique visualization techniques for microbiome data using `phyloseq` and `microViz`. Its flagship features involve replacing standard points or bars with "star" (radar) plots that show the relative abundance of core taxa or any taxa of interest. These multi-dimensional shapes can be arranged in core matrices or used to summarise several ordination axes at once, making it easier to compare community structures across samples or treatments (up to 4 groups/treatments for comparative core matrices).

<img width="1315" height="747" alt="image" src="https://github.com/user-attachments/assets/89dd71a7-bb12-43b5-979c-83a2a9a43128" />
*Mean position on the first five PCoA axes for each sample type in the GlobalPatterns dataset (`plot_ordi_star()`)*

Transparency note! I used Gemini and Claude to help identify and fix inconsistencies in the readme and between subcommands.

## Installation

`CompStellation` depends on `phyloseq` (Bioconductor) and `microViz` (GitHub), so those need to be in place first:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

BiocManager::install("phyloseq")
remotes::install_github("david-barnett/microViz")
remotes::install_github("kadenmuffett/CompStellation")
```

## Features

### `plot_taxa_star()`
Creates radar/star plots to visualize the relative abundance of specified taxa (or the top 4 most abundant taxa automatically, up to a maximum of 8) across different sample groups, complete with standard error or IQR shading overlays.

### `plot_ordi_star()`
Runs an ordination (PCoA by default) on a beta-diversity distance matrix such as Bray-Curtis, then summarises each group's position by plotting its **mean score on the first *n* ordination axes as the spokes of a star plot**. This lets you compare groups across five axes at once rather than the two a conventional ordination scatterplot shows. Individual samples and taxonomic composition are not drawn — for taxonomy alongside ordination, see `plot_ordi_star_iris()`.

### `plot_ordi_star_iris()` *(experimental)*
Places the `plot_ordi_star()` radar at the centre of a circular stacked bar chart ("iris"), where each segment of the outer ring is one sample's taxonomic composition. Useful for seeing ordination position and community composition in a single figure.

### `plot_core_matrix()`
Generates a facet matrix displaying the "unique core" ASVs across groups. Shows the abundance of taxa that are considered "core" to one treatment/group but completely missing from the core of the other groups.

### `plot_concordant_core()`
The complement of `plot_core_matrix()`: identifies the taxa that are core to *every* evaluated group and plots their relative abundance as a single star plot.

### `coreselect()`
Provides a presence/absence dot plot for assessing core microbiome thresholds across various groupings within your phyloseq object. This tool is exploratory and can be used for selection of appropriate taxa for `plot_taxa_star()` or to assess the appropriateness of `plot_core_matrix()` for your data.

> **Note:** the core functions (`coreselect()`, `plot_core_matrix()`, `plot_concordant_core()`) need groups with substantial taxonomic overlap and a meaningful exclusive core. `GlobalPatterns` does not have that structure — the vignette builds a small simulated dataset that does. Run `vignette("CompStellation")` for those examples.

## Examples

```r
library(CompStellation)
library(phyloseq)
library(microViz)

data("GlobalPatterns")
# Subset to a smaller physeq for demonstration
gp <- subset_samples(GlobalPatterns, SampleType %in% c("Feces", "Ocean", "Soil"))
# tax_fix() fills unassigned ranks so taxa can be aggregated; then drop empty taxa
gp <- tax_fix(gp)
gp <- subset_taxa(gp, taxa_sums(gp) > 0)

# Plot Taxa Star Plot with Standard Error shading
plot_taxa_star(
  physeq = gp,
  sample_var = "SampleType",
  taxa_rank = "Phylum",
  samplecolumn = "X.SampleID",
  error_bar = "SE"
)

# Summarise each group's mean position on the first 5 PCoA axes
# (Bray-Curtis distance), with the groups overlaid on one radar
plot_ordi_star(
  physeq = gp,
  sample_var = "SampleType",
  distance = "bray",
  view_type = "together"
)

# The same thing faceted, one panel per group
plot_ordi_star(
  physeq = gp,
  sample_var = "SampleType",
  distance = "bray",
  view_type = "separate"
)
```

For `plot_core_matrix()`, `plot_concordant_core()` and `coreselect()`, see `vignette("CompStellation")` — those functions need the simulated dataset built there rather than `GlobalPatterns`.

## Interpreting star plots

A few things are worth keeping in mind when reading these figures, and worth stating in a caption if you publish one:

- **Area is not a quantity.** The area enclosed by a star polygon depends on the order of the spokes, so it changes if the axes are reordered even though no value has changed. Compare spoke by spoke, not by overall size.
- **The radial axis in `plot_ordi_star()` is offset.** Ordination scores are negative as well as positive, so all axes are shifted by a constant to make them plottable on a radius. Radii are therefore not proportional to the scores, and are not comparable between axes. The tick labels show the original (unshifted) values.
- **`error_bar = "IQR"` and `error_bar = "SE"` answer different questions.** The IQR describes the spread of samples within a group (and is drawn around the group median); the SE describes the uncertainty of the group mean (and is drawn around the mean). The plot title records which one was used.
- **The abundances are compositional.** Relative abundances within a sample sum to 1, so the spokes are not independent — a rise in one taxon mechanically depresses the others, and the "Other" wedge is usually the largest single value.
- **These plots are descriptive.** They do not test whether groups differ. Pair them with an appropriate test (for example `vegan::adonis2()` on the same distance matrix) before making a claim about group separation.

## Getting help

If you encounter a clear bug, please file an issue with a minimal reproducible example at <https://github.com/kadenmuffett/CompStellation/issues>.
