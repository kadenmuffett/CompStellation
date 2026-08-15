# Shared helpers for the core-microbiome functions (coreselect,
# plot_core_matrix, plot_concordant_core). These three functions previously
# carried near-identical copies of this logic, which had drifted apart.

#' Resolve a prevalence threshold to a proportion
#'
#' Accepts a threshold supplied either as a proportion in (0, 1] or as a
#' percentage in (1, 100], and returns a proportion. The chosen interpretation
#' is reported, because a bare `1` is ambiguous: it is read as 100%, not 1%.
#'
#' @param percent_samples The user-supplied threshold.
#' @param verbose Logical. If TRUE (default), report the resolved threshold.
#'
#' @return A single numeric value in (0, 1].
#' @keywords internal
resolve_prevalence <- function(percent_samples, verbose = TRUE) {
  if (!is.numeric(percent_samples) || length(percent_samples) != 1L ||
      is.na(percent_samples)) {
    stop("'percent_samples' must be a single non-missing number.", call. = FALSE)
  }

  if (percent_samples <= 0 || percent_samples > 100) {
    stop(
      "'percent_samples' must be a proportion in (0, 1] or a percentage in (1, 100]. Got ",
      percent_samples, ".",
      call. = FALSE
    )
  }

  prop <- if (percent_samples > 1) percent_samples / 100 else percent_samples

  if (verbose) {
    message(
      "Prevalence threshold: taxa must be present in at least ",
      format(prop * 100, trim = TRUE), "% of the samples in a group.",
      if (identical(percent_samples, 1)) {
        " (A bare 1 is read as the proportion 1, i.e. 100%; use 0.01 for 1%.)"
      } else {
        ""
      }
    )
  }

  prop
}

#' Per-taxon prevalence within a set of samples
#'
#' @param physeq_sub A phyloseq object already subset to a single group.
#' @param abundance_threshold The abundance a taxon must exceed in a sample to
#'   count as present. Strictly greater than, so the default of 0 means
#'   "any non-zero value".
#' @param abundance_type Either "counts" (the threshold is applied to the
#'   values as supplied) or "relative" (each sample is converted to proportions
#'   first, so the threshold is a relative abundance).
#'
#' @return A named numeric vector of prevalences in \[0, 1\], one per taxon.
#' @keywords internal
taxa_prevalence <- function(physeq_sub, abundance_threshold = 0,
                            abundance_type = c("counts", "relative")) {
  abundance_type <- match.arg(abundance_type)

  if (abundance_type == "relative") {
    physeq_sub <- phyloseq::transform_sample_counts(physeq_sub, to_relative)
  }

  otu_tab <- phyloseq::otu_table(physeq_sub)
  if (!phyloseq::taxa_are_rows(otu_tab)) otu_tab <- phyloseq::t(otu_tab)
  mat <- as.matrix(otu_tab)

  rowSums(mat > abundance_threshold) / ncol(mat)
}

#' Convert one sample's abundances to proportions
#'
#' Guards against division by zero, which would otherwise turn an empty sample
#' into a row of `NaN` that propagates silently into every downstream mean.
#'
#' @param x A numeric vector of abundances for one sample.
#'
#' @return A numeric vector of the same length.
#' @keywords internal
to_relative <- function(x) {
  s <- sum(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) x else x / s
}

#' Warn when library sizes are uneven enough to distort prevalence
#'
#' Presence/absence on untransformed counts is sensitive to sequencing depth:
#' a taxon registers as present more often in deeply sequenced samples, so a
#' prevalence-based core partly measures sequencing effort.
#'
#' @param physeq A phyloseq object.
#' @param ratio The fold-difference between the largest and smallest library
#'   at which to warn.
#'
#' @return Invisibly `NULL`, called for the side effect.
#' @keywords internal
warn_library_sizes <- function(physeq, ratio = 10) {
  depths <- phyloseq::sample_sums(physeq)
  depths <- depths[is.finite(depths) & depths > 0]
  if (length(depths) < 2) {
    return(invisible(NULL))
  }

  if (max(depths) / min(depths) >= ratio) {
    warning(
      "Library sizes vary more than ", ratio, "-fold (",
      format(round(min(depths))), " to ", format(round(max(depths))),
      "). Presence/absence, and therefore the core, is sensitive to sequencing ",
      "depth. Consider rarefying, or using abundance_type = \"relative\" with a ",
      "non-zero abundance_threshold.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

#' Order taxa by mean relative abundance
#'
#' Used when a core has to be truncated for legibility, so that the taxa kept
#' are the most abundant ones rather than whichever happened to come first in
#' the taxonomy table.
#'
#' @param physeq A phyloseq object.
#' @param taxa A character vector of taxa names to order.
#' @param samples Optional character vector of sample names to restrict to.
#'
#' @return `taxa`, reordered by decreasing mean relative abundance.
#' @keywords internal
order_taxa_by_abundance <- function(physeq, taxa, samples = NULL) {
  taxa <- intersect(taxa, phyloseq::taxa_names(physeq))
  if (length(taxa) < 2) {
    return(taxa)
  }

  ps <- physeq
  if (!is.null(samples)) ps <- phyloseq::prune_samples(samples, ps)
  ps <- phyloseq::transform_sample_counts(ps, to_relative)
  ps <- phyloseq::prune_taxa(taxa, ps)

  means <- phyloseq::taxa_sums(ps) / phyloseq::nsamples(ps)
  names(sort(means, decreasing = TRUE))
}

#' Put "Other" last in a taxon grouping
#'
#' The enclosed area of a star plot depends on the order of its spokes, so the
#' order should be deliberate rather than whatever `dplyr` happened to produce.
#' Named taxa are ordered as supplied and the pooled "Other" wedge is pinned to
#' the end.
#'
#' @param x A character vector of taxon-group labels.
#' @param levels_order Optional character vector giving the desired order of the
#'   named taxa.
#'
#' @return A factor.
#' @keywords internal
order_taxa_group <- function(x, levels_order = NULL) {
  x <- as.character(x)
  named <- setdiff(unique(x), "Other")

  if (!is.null(levels_order)) {
    named <- c(intersect(levels_order, named), setdiff(named, levels_order))
  } else {
    named <- sort(named)
  }

  lv <- c(named, if ("Other" %in% x) "Other")
  factor(x, levels = lv)
}
