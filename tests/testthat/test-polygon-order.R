# geom_polygon() joins vertices in row order, not in x-axis order. If the
# factor levels on x disagree with the row order, the polygon crosses over
# itself and the star shows a shape the data does not have. This is easy to
# reintroduce whenever the spoke order is changed, so pin it.

vertices_are_in_axis_order <- function(p, x_var, group_var) {
  d <- p$data
  groups <- unique(as.character(d[[group_var]]))
  all(vapply(groups, function(g) {
    lv <- as.integer(d[[x_var]][as.character(d[[group_var]]) == g])
    identical(lv, sort(lv))
  }, logical(1)))
}

test_that("plot_taxa_star traces its spokes in axis order", {
  sim <- make_sim()

  for (eb in c("IQR", "SE", "none")) {
    p <- plot_taxa_star(sim, "Group", taxa_rank = "OTU",
                        samplecolumn = "SampleID", error_bar = eb)
    expect_true(
      vertices_are_in_axis_order(p, "Taxa_Group", "Group"),
      info = paste("error_bar =", eb)
    )
  }
})

test_that("plot_taxa_star vertex order survives a custom spoke order", {
  sim <- make_sim()
  p <- plot_taxa_star(sim, "Group", taxa_rank = "OTU", samplecolumn = "SampleID",
                      taxa_names = c("Taxon_F", "Taxon_A", "Taxon_C"))

  # "Other" must still be the last level, and the rows must follow the levels
  lv <- levels(p$data$Taxa_Group)
  expect_identical(lv[length(lv)], "Other")
  expect_true(vertices_are_in_axis_order(p, "Taxa_Group", "Group"))
})

test_that("what ggplot actually builds is a simple, non-crossing polygon", {
  sim <- make_sim()
  p <- plot_taxa_star(sim, "Group", taxa_rank = "OTU",
                      samplecolumn = "SampleID", error_bar = "IQR")

  built <- ggplot2::ggplot_build(p)$data[[1]]
  for (g in unique(built$group)) {
    xs <- built$x[built$group == g]
    expect_false(is.unsorted(xs),
                 label = paste("built polygon x-order for group", g))
  }
})

test_that("plot_ordi_star traces its axes in order", {
  sim <- make_sim()
  p <- plot_ordi_star(sim, "Group", distance = "bray", verbose = FALSE)
  expect_true(vertices_are_in_axis_order(p, "Axis", "Group"))
})

test_that("axis order holds beyond nine axes, where labels stop sorting numerically", {
  # "Axis 10" sorts before "Axis 2" alphabetically, so relying on dplyr's row
  # order silently breaks here.
  set.seed(3)
  n <- 24
  m <- matrix(stats::rpois(n * 80, 5), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), paste0("T", seq_len(80))))
  md <- data.frame(Group = rep(c("A", "B", "C"), length.out = n),
                   row.names = rownames(m))
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(m, taxa_are_rows = FALSE),
    phyloseq::sample_data(md)
  )

  p <- plot_ordi_star(ps, "Group", distance = "bray", n_axes = 12, verbose = FALSE)
  expect_gt(nlevels(p$data$Axis), 9)
  expect_true(vertices_are_in_axis_order(p, "Axis", "Group"))
})
