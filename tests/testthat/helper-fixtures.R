# A small simulated community with a known core structure, used across the
# tests. GlobalPatterns is not suitable for the core functions (see README), and
# the package does not yet ship an example dataset.
make_sim <- function(seed = 1, n_per = 12) {
  set.seed(seed)
  grps <- c("Control", "TreatmentA", "TreatmentB")
  n_s <- n_per * length(grps)
  taxa <- paste0("Taxon_", LETTERS[1:12])

  om <- matrix(0, nrow = length(taxa), ncol = n_s,
               dimnames = list(taxa, paste0("Sample", seq_len(n_s))))
  gv <- rep(grps, each = n_per)

  for (i in seq_len(n_s)) {
    g <- gv[i]
    # Core to every group, so it shows up as the concordant core
    om["Taxon_A", i] <- 100
    om["Taxon_B", i] <- 40
    # Core to one group only, so it shows up as that group's unique core
    if (g == "Control")    om["Taxon_C", i] <- 200
    if (g == "TreatmentA") om["Taxon_D", i] <- 150
    if (g == "TreatmentB") om["Taxon_E", i] <- 300
    # Right-skewed: one sample per group carries a large outlier, which is what
    # makes the mean and the median differ
    if (i %% n_per == 1)   om["Taxon_F", i] <- 5000 else om["Taxon_F", i] <- 5
    # Noise
    for (k in c("Taxon_G", "Taxon_H", "Taxon_I")) {
      if (stats::runif(1) < 0.5) om[k, i] <- 10
    }
  }

  tt <- as.matrix(taxa)
  rownames(tt) <- taxa
  colnames(tt) <- "Taxon"

  md <- data.frame(
    SampleID = colnames(om),
    Group = gv,
    row.names = colnames(om)
  )

  phyloseq::phyloseq(
    phyloseq::otu_table(om, taxa_are_rows = TRUE),
    phyloseq::tax_table(tt),
    phyloseq::sample_data(md)
  )
}
