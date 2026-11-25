#' Check Core Microbiome Counts Across Prevalence Thresholds
#'
#' Plots a faceted bar graph showing the number of core taxa in each sample
#' group across prevalence thresholds (`percent_samples`) from `0.25` to `1`,
#' incremented by `step_size`.
#'
#' @param physeq A phyloseq object containing your data.
#' @param group_var A character string. The name of the column in your sample
#'   data to use for grouping/faceting samples (e.g., `"Treatment"`, `"Body_cav"`).
#' @param abundance_threshold A numeric value (default `0`). The abundance a
#'   taxon must **exceed** in a sample for that sample to count as present.
#' @param abundance_type Either `"counts"` (default) or `"relative"`; see
#'   [coreselect()].
#' @param taxa_rank A character string specifying the taxonomic rank to
#'   aggregate to before evaluating prevalence (e.g., `"Family"`, `"Genus"`),
#'   or `"OTU"` (default) for unaggregated taxa.
#' @param step_size A single numeric value in `(0, 1]` (default `0.05`) giving
#'   the step increment between prevalence thresholds from `0.25` to `1`.
#' @param verbose Logical. If `TRUE` (default), report aggregation and attribute
#'   access instructions.
#' @param font Optional character string specifying the font family for plot
#'   text (e.g., `"sans"`, `"serif"`). If `NULL` (default), uses the active
#'   ggplot2 theme default.
#' @param label_color Character string specifying the color for all plot text
#'   and labels (default `"black"`).
#'
#' @return A `ggplot` object faceted by `group_var`, with the underlying core
#'   counts table attached as `attr(plot, "prev_check")`.
#'
#' @import phyloseq
#' @import ggplot2
#' @export
#'
#' @examples
#' data("cassiopea_ps")
#' prev_check(
#'   cassiopea_ps,
#'   group_var = "Body_cav",
#'   taxa_rank = "Family",
#'   step_size = 0.05
#' )
prev_check <- function(physeq, group_var, abundance_threshold = 0,
                       abundance_type = c("counts", "relative"),
                       taxa_rank = "OTU", step_size = 0.05,
                       verbose = TRUE, font = NULL,
                       label_color = "black") {
  if (!inherits(physeq, "phyloseq")) {
    stop("'physeq' must be a VALID phyloseq object.", call. = FALSE)
  }
  if (!group_var %in% phyloseq::sample_variables(physeq)) {
    stop(
      "'", group_var, "' is not a sample variable in the phyloseq object. Available: ",
      paste(phyloseq::sample_variables(physeq), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.character(taxa_rank) || length(taxa_rank) != 1L || is.na(taxa_rank) || !nzchar(taxa_rank)) {
    stop("'taxa_rank' must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.numeric(step_size) || length(step_size) != 1L || is.na(step_size) ||
      step_size <= 0 || step_size > 1) {
    stop("'step_size' must be a single number in (0, 1].", call. = FALSE)
  }

  abundance_type <- match.arg(abundance_type)
  if (abundance_type == "counts") warn_library_sizes(physeq)

  physeq_glom <- if (taxa_rank != "OTU") {
    agg_taxa(physeq, rank = taxa_rank, verbose = verbose)
  } else {
    if (verbose) message("Taxa not aggregated. Using OTUs.")
    physeq
  }

  sample_df <- data.frame(phyloseq::sample_data(physeq_glom), check.names = FALSE)
  groups <- unique(sample_df[[group_var]])
  thresholds <- round(seq(0.25, 1, by = step_size), 10)

  rows <- lapply(groups, function(g) {
    samples_g <- rownames(sample_df)[sample_df[[group_var]] == g]
    if (length(samples_g) == 0) return(NULL)
    sub_ps <- phyloseq::prune_samples(samples_g, physeq_glom)
    prev <- taxa_prevalence(sub_ps, abundance_threshold = abundance_threshold,
                            abundance_type = abundance_type)
    data.frame(
      Group = as.character(g),
      Threshold = thresholds,
      Core_Taxa = vapply(thresholds, function(t) sum(prev >= t), integer(1)),
      stringsAsFactors = FALSE
    )
  })
  plot_data <- do.call(rbind, rows)
  names(plot_data)[1] <- group_var

  p <- ggplot(plot_data, aes(x = Threshold, y = Core_Taxa)) +
    geom_col(width = step_size * 0.85, fill = "black") +
    facet_wrap(stats::as.formula(paste("~", group_var))) +
    theme_minimal() +
    labs(
      title = "Core Microbiome by Prevalence Thresholds",
      x = "Prevalence Threshold (percent_samples)",
      y = "Number of Core Taxa"
    ) +
    theme(
      text = element_text(color = label_color),
      axis.text = element_text(color = label_color),
      axis.title = element_text(color = label_color),
      strip.text = element_text(color = label_color),
      plot.title = element_text(color = label_color)
    )

  if (!is.null(font)) {
    p <- p + theme(text = element_text(family = font, color = label_color))
  }

  attr(p, "prev_check") <- plot_data
  if (verbose) {
    message("Attached attributes: access core counts table via attr(plot, \"prev_check\").")
  }
  p
}
