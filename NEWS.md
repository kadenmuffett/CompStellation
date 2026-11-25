# CompStellation 0.1.0

* Initial CRAN release.
* Provides compositional and multidimensional star (radar) plots for ecological and microbiome datasets:
  * `plot_taxa_star()` for group-level relative abundance star plots of focal or top-ranked taxa with IQR or SE intervals.
  * `plot_ordi_star()` and experimental `plot_ordi_star_iris()` for displaying group centers across multiple ordination axes simultaneously.
  * `coreselect()`, `prev_check()`, `plot_core_matrix()`, and `plot_concordant_core()` for evaluating prevalence thresholds and visualizing group-specific, unique, and shared core microbiomes.
  * `vegan_to_phyloseq()` for converting `vegan`-style sample-by-taxon community tables and metadata into `phyloseq` objects.
* Includes the `cassiopea_ps` 16S rRNA microbiome dataset of *Cassiopea xamachana* and surrounding environmental compartments.
