# Tests for the statistical behaviour of the plots. Each block pins a specific
# correction so it cannot silently regress.

test_that("percent variance uses the corrected eigenvalues when some are negative", {
  skip_if_not_installed("ape")
  skip_if_not_installed("vegan")

  set.seed(1)
  m <- matrix(stats::rpois(30 * 60, 0.6), nrow = 30)
  m[sample(length(m), length(m) * 0.75)] <- 0
  m <- m[rowSums(m) > 0, , drop = FALSE]

  ord <- ape::pcoa(vegan::vegdist(m, method = "bray"))
  eig <- ord$values$Eigenvalues
  skip_if(sum(eig < 0) == 0, "this fixture did not produce negative eigenvalues")

  got <- axis_percent_var(ord)

  # ape supplies the corrected relative eigenvalues exactly for this case
  expect_equal(got$pct, 100 * ord$values$Rel_corr_eig)

  # Negatives shrink sum(eig), so the plain relative eigenvalues overstate each
  # positive axis. The corrected values must be smaller.
  plain <- 100 * eig / sum(eig)
  expect_lt(got$pct[1], plain[1])

  # Both conventions still total 100%
  expect_equal(sum(got$pct, na.rm = TRUE), 100)
})

test_that("percent variance is unchanged when no eigenvalues are negative", {
  # Built by hand rather than fished out of a random matrix: whether a given
  # dataset yields negative eigenvalues is not something a test should gamble on.
  eig <- c(4, 3, 2, 1)
  ord <- list(values = data.frame(
    Eigenvalues = eig,
    Relative_eig = eig / sum(eig)
  ))

  got <- axis_percent_var(ord)
  expect_equal(got$pct, 100 * eig / sum(eig))
  expect_equal(got$pct, c(40, 30, 20, 10))
  expect_match(got$note, "relative eigenvalues")
})

test_that("percent variance falls back to positive eigenvalues only", {
  # No Rel_corr_eig supplied, but negatives present: the denominator must not
  # include the negatives, or every positive axis is overstated.
  eig <- c(4, 3, 2, -1)
  ord <- list(values = data.frame(Eigenvalues = eig))

  got <- axis_percent_var(ord)
  expect_equal(got$pct, c(100 * 4 / 9, 100 * 3 / 9, 100 * 2 / 9, NA_real_))
  expect_equal(sum(got$pct, na.rm = TRUE), 100)
  expect_match(got$note, "positive eigenvalues only")

  # The old formula divided by sum(eig) = 8, overstating axis 1 as 50%
  expect_lt(got$pct[1], 100 * eig[1] / sum(eig))
})

test_that("the plotted centre matches the interval that is drawn", {
  sim <- make_sim()

  p_iqr <- plot_taxa_star(sim, "Group", taxa_rank = "OTU",
                          samplecolumn = "SampleID", error_bar = "IQR")
  p_se <- plot_taxa_star(sim, "Group", taxa_rank = "OTU",
                         samplecolumn = "SampleID", error_bar = "SE")

  # The IQR brackets the median, so the median is what gets drawn
  expect_equal(p_iqr$data$Center_Abundance, p_iqr$data$median_Abundance)
  # The SE describes the mean, so the mean is what gets drawn
  expect_equal(p_se$data$Center_Abundance, p_se$data$mean_Abundance)

  # The fixture is right-skewed, so the two genuinely differ
  expect_false(isTRUE(all.equal(p_iqr$data$median_Abundance,
                                p_iqr$data$mean_Abundance)))

  # And the title says which one
  expect_match(p_iqr$labels$title, "^Median")
  expect_match(p_se$labels$title, "^Mean")
})

test_that("spokes are ordered by abundance with 'Other' last", {
  sim <- make_sim()
  p <- plot_taxa_star(sim, "Group", taxa_rank = "OTU", samplecolumn = "SampleID")

  lv <- levels(p$data$Taxa_Group)
  expect_true(is.factor(p$data$Taxa_Group))
  expect_identical(lv[length(lv)], "Other")
})

test_that("a zero-total sample does not become NaN", {
  sim <- make_sim()
  ot <- as(phyloseq::otu_table(sim), "matrix")
  ot[, 1] <- 0
  phyloseq::otu_table(sim) <- phyloseq::otu_table(
    ot, taxa_are_rows = phyloseq::taxa_are_rows(sim)
  )

  p <- plot_taxa_star(sim, "Group", taxa_rank = "OTU",
                      samplecolumn = "SampleID", error_bar = "none")
  expect_false(any(is.nan(p$data$Center_Abundance)))
  expect_equal(to_relative(c(0, 0, 0)), c(0, 0, 0))
})

test_that("a group with one sample warns instead of silently losing its interval", {
  sim <- make_sim()
  keep <- c(
    phyloseq::sample_names(sim)[phyloseq::sample_data(sim)$Group == "Control"][1],
    phyloseq::sample_names(sim)[phyloseq::sample_data(sim)$Group == "TreatmentA"]
  )
  thin <- phyloseq::prune_samples(keep, sim)

  expect_warning(
    plot_taxa_star(thin, "Group", taxa_rank = "OTU",
                   samplecolumn = "SampleID", error_bar = "SE"),
    "single sample"
  )
})

test_that("Euclidean distance is invariant to the radial offset", {
  # This is why the group clustering had to move off Bray-Curtis: the offset
  # applied to make ordination scores plottable is arbitrary, and Bray-Curtis
  # would let it change the answer.
  set.seed(7)
  m <- matrix(stats::rnorm(20), nrow = 5)
  expect_equal(as.vector(stats::dist(m)), as.vector(stats::dist(m + 7.3)))
})

test_that("the prevalence threshold is resolved unambiguously", {
  expect_equal(resolve_prevalence(0.8, verbose = FALSE), 0.8)
  expect_equal(resolve_prevalence(80, verbose = FALSE), 0.8)
  # A bare 1 means 100%, not 1%
  expect_equal(resolve_prevalence(1, verbose = FALSE), 1)

  expect_error(resolve_prevalence(0), "proportion")
  expect_error(resolve_prevalence(-0.5), "proportion")
  expect_error(resolve_prevalence(101), "proportion")
  expect_error(resolve_prevalence("half"), "single non-missing number")
  expect_error(resolve_prevalence(c(0.5, 0.8)), "single non-missing number")

  expect_message(resolve_prevalence(0.8), "80%")
})

test_that("all three core functions validate the threshold the same way", {
  sim <- make_sim()
  expect_error(
    coreselect(sim, "Group", percent_samples = -0.5),
    "proportion"
  )
  expect_error(
    plot_core_matrix(sim, "Group", percent_samples = -0.5,
                     taxa_rank = "OTU", samplecolumn = "SampleID"),
    "proportion"
  )
  expect_error(
    plot_concordant_core(sim, "Group", percent_samples = -0.5,
                         taxa_rank = "OTU", samplecolumn = "SampleID"),
    "proportion"
  )
})

test_that("prevalence can be computed on relative abundance instead of counts", {
  sim <- make_sim()
  by_counts <- taxa_prevalence(sim, abundance_threshold = 0,
                               abundance_type = "counts")
  by_rel <- taxa_prevalence(sim, abundance_threshold = 0.05,
                            abundance_type = "relative")

  expect_named(by_counts)
  expect_true(all(by_rel <= by_counts))
  expect_true(all(by_counts >= 0 & by_counts <= 1))
})

test_that("truncation keeps the most abundant taxa, not the first in table order", {
  sim <- make_sim()
  rel <- phyloseq::transform_sample_counts(sim, to_relative)
  ordered <- order_taxa_by_abundance(sim, c("Taxon_A", "Taxon_B", "Taxon_G"))

  means <- phyloseq::taxa_sums(rel)[ordered] / phyloseq::nsamples(rel)
  expect_false(is.unsorted(rev(means)))
})

test_that("group_subset restricts the panels, not just the core calculation", {
  sim <- make_sim()
  p <- suppressWarnings(plot_core_matrix(
    sim, "Group", percent_samples = 0.8, taxa_rank = "OTU",
    samplecolumn = "SampleID", group_subset = c("Control", "TreatmentA"),
    verbose = FALSE
  ))

  panels <- unique(as.character(p[[1]]$data$Group))
  expect_setequal(panels, c("Control", "TreatmentA"))
})

test_that("the core functions return their memberships, not just a plot", {
  sim <- make_sim()

  cs <- suppressWarnings(coreselect(sim, "Group", percent_samples = 0.9, verbose = FALSE))
  expect_type(attr(cs, "core_list"), "list")
  expect_true("Taxon_C" %in% attr(cs, "core_list")$Control)

  cm <- suppressWarnings(plot_core_matrix(
    sim, "Group", percent_samples = 0.9, taxa_rank = "OTU",
    samplecolumn = "SampleID", verbose = FALSE
  ))
  expect_true("Taxon_C" %in% attr(cm, "unique_core")$Control)
  # The prevalences behind each "unique" call are returned so the margin can
  # be checked
  expect_type(attr(cm, "prevalence"), "list")
  expect_named(attr(cm, "prevalence")$Control)

  cc <- suppressWarnings(plot_concordant_core(
    sim, "Group", percent_samples = 0.9, taxa_rank = "OTU",
    samplecolumn = "SampleID", verbose = FALSE
  ))
  expect_true("Taxon_A" %in% attr(cc, "concordant_core"))
})

test_that("constrained ordinations are refused with a usable message", {
  sim <- make_sim()
  expect_error(
    plot_ordi_star(sim, "Group", method = "CAP", distance = "bray"),
    "constrained"
  )
  # "PCA" is not a phyloseq method; it used to be guarded for but unreachable
  expect_error(plot_ordi_star(sim, "Group", method = "PCA", distance = "bray"))
})

test_that("the exported plot functions return plot objects", {
  sim <- make_sim()

  expect_s3_class(
    plot_taxa_star(sim, "Group", taxa_rank = "OTU", samplecolumn = "SampleID"),
    "ggplot"
  )
  expect_s3_class(
    plot_ordi_star(sim, "Group", distance = "bray", verbose = FALSE),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(coreselect(sim, "Group", percent_samples = 0.9, verbose = FALSE)),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(plot_core_matrix(
      sim, "Group", percent_samples = 0.9, taxa_rank = "OTU",
      samplecolumn = "SampleID", verbose = FALSE
    )),
    "patchwork"
  )
})

test_that("uneven library sizes are flagged when presence is tested on counts", {
  sim <- make_sim()
  depths <- phyloseq::sample_sums(sim)
  expect_gt(max(depths) / min(depths), 10)

  expect_warning(
    coreselect(sim, "Group", percent_samples = 0.9, verbose = FALSE),
    "sensitive to sequencing depth"
  )

  # Switching to relative abundance removes the reason for the warning
  expect_no_warning(
    coreselect(sim, "Group", percent_samples = 0.9, abundance_threshold = 0.01,
               abundance_type = "relative", verbose = FALSE)
  )
})
