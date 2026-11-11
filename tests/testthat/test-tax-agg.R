# Unit tests for internal taxonomic aggregation helper agg_taxa()

make_tax_sim <- function() {
  otu <- matrix(
    c(
      10, 20, 30, 40,
      5,  15, 25, 35,
      12, 18, 22, 28,
      8,  12, 16, 20,
      3,  7,  11, 15,
      2,  4,  6,  8
    ),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("ASV", 1:6), paste0("S", 1:4))
  )

  tax <- matrix(
    c(
      "Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Bacteroides",
      "Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Bacteroides",
      "Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "unclassified",
      "Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales", "g__",            NA,
      "Bacteria", "Firmicutes",   "Clostridia",  "Lachnospirales","Lachnospiraceae","Incertae Sedis",
      "Bacteria", "Firmicutes",   "Clostridia",  "Oscillospirales","Ruminococcaceae","Incertae Sedis"
    ),
    nrow = 6, byrow = TRUE,
    dimnames = list(
      paste0("ASV", 1:6),
      c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    )
  )

  sd <- data.frame(
    SampleID = paste0("S", 1:4),
    Group = c("A", "A", "B", "B"),
    row.names = paste0("S", 1:4)
  )

  phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax),
    phyloseq::sample_data(sd)
  )
}

test_that("agg_taxa sums counts, fills unclassified placeholders, and emits message", {
  ps <- make_tax_sim()

  expect_message(
    agg <- agg_taxa(ps, "Genus"),
    "Taxa aggregated to rank: Genus using an internal heuristic to remove unclassifieds. Some unclassifieds may remain."
  )

  # ASV1 and ASV2 share Bacteroides -> counts summed (10+5=15, 20+15=35, etc.)
  otu_agg <- methods::as(phyloseq::otu_table(agg), "matrix")
  expect_equal(as.numeric(otu_agg["Bacteroides", ]), c(15, 35, 55, 75))

  # ASV3 has Genus = "unclassified", filled from Family -> "Bacteroidaceae Family"
  expect_true("Bacteroidaceae Family" %in% phyloseq::taxa_names(agg))
  expect_equal(as.numeric(otu_agg["Bacteroidaceae Family", ]), c(12, 18, 22, 28))

  # ASV4 has Family = "g__" and Genus = NA, filled from Order -> "Bacteroidales Order"
  expect_true("Bacteroidales Order" %in% phyloseq::taxa_names(agg))
  expect_equal(as.numeric(otu_agg["Bacteroidales Order", ]), c(8, 12, 16, 20))

  # ASV5 and ASV6 both have "Incertae Sedis" (a placeholder), so they get filled
  # from their respective families: "Lachnospiraceae Family" and "Ruminococcaceae Family"
  expect_true("Lachnospiraceae Family" %in% phyloseq::taxa_names(agg))
  expect_true("Ruminococcaceae Family" %in% phyloseq::taxa_names(agg))

  # Total library sizes are preserved
  expect_equal(phyloseq::sample_sums(agg), phyloseq::sample_sums(ps))
})

test_that("agg_taxa disambiguates identical taxon names under different parent lineages", {
  ps <- make_tax_sim()
  tt <- methods::as(phyloseq::tax_table(ps), "matrix")
  # Give ASV5 and ASV6 the same non-placeholder genus name under different families
  tt["ASV5", "Genus"] <- "SharedGenus"
  tt["ASV6", "Genus"] <- "SharedGenus"
  phyloseq::tax_table(ps) <- phyloseq::tax_table(tt)

  agg <- suppressMessages(agg_taxa(ps, "Genus"))
  expect_true("SharedGenus (Lachnospiraceae)" %in% phyloseq::taxa_names(agg))
  expect_true("SharedGenus (Ruminococcaceae)" %in% phyloseq::taxa_names(agg))
})

test_that("plot_taxa_star, plot_core_matrix, and plot_concordant_core work with agg_taxa", {
  ps <- make_tax_sim()

  expect_message(
    p_taxa <- plot_taxa_star(ps, "Group", taxa_rank = "Phylum", samplecolumn = "SampleID"),
    "Taxa aggregated to rank: Phylum using an internal heuristic to remove unclassifieds. Some unclassifieds may remain."
  )
  expect_s3_class(p_taxa, "ggplot")

  expect_message(
    p_cc <- suppressWarnings(plot_concordant_core(
      ps, "Group", percent_samples = 0.5, taxa_rank = "Family",
      samplecolumn = "SampleID", verbose = FALSE
    )),
    "Taxa aggregated to rank: Family using an internal heuristic to remove unclassifieds. Some unclassifieds may remain."
  )
  expect_s3_class(p_cc, "ggplot")
})
