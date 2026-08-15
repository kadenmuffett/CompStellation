#' Plot Concordant Core Taxa
#'
#' This function identifies the concordant (overlapping) core taxa across all evaluated groups
#' and plots their relative abundance using `plot_taxa_star`.
#'
#' @param physeq A phyloseq object.
#' @param group_var A character string. The name of the column to use for grouping.
#' @param percent_samples The proportion of samples within a group that a taxon
#'   must be present in to be considered "core". Given either as a proportion in
#'   (0, 1] or as a percentage in (1, 100]; the interpretation chosen is
#'   reported, because a bare `1` means 100%, not 1%.
#' @param abundance_threshold A numeric value (default 0). The abundance a taxon
#'   must **exceed** in a sample for that sample to count as present.
#' @param abundance_type Either "counts" (default) or "relative"; see [coreselect()].
#' @param taxa_rank A character string. The taxonomic rank to use (e.g., "Genus").
#' @param samplecolumn A character string. The name of the sample ID column.
#' @param log_scale A logical value. If TRUE, applies a pseudo-log transformation to the y-axis.
#' @param group_subset A character vector specifying which groups to include (if NULL, all are used).
#' @param verbose Logical. If TRUE (default), report the resolved thresholds and
#'   the size of each group's core.
#' @param ... Additional arguments passed to `plot_taxa_star`, e.g. "colors_all".
#'
#' @return A ggplot object representing the star plot of concordant core taxa,
#'   with the shared core attached as the "concordant_core" attribute and the
#'   per-group core sizes as the "core_sizes" attribute.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
#' @examples
#' # plot_concordant_core(physeq, "Treatment", percent_samples = 0.5,
#' #                      taxa_rank = "Genus", samplecolumn = "SampleID")
plot_concordant_core <- function(physeq, group_var, percent_samples, abundance_threshold = 0,
                                 abundance_type = c("counts", "relative"), taxa_rank, samplecolumn,
                                 log_scale = FALSE, group_subset = NULL, verbose = TRUE, ...) {
    # --- 1. Input Validation ---

    if (!inherits(physeq, "phyloseq")) {
        stop("'physeq' must be a VALID phyloseq object.", call. = FALSE)
    }

    abundance_type <- match.arg(abundance_type)

    # Resolve and report the prevalence threshold (shared with coreselect() and
    # plot_core_matrix(), which previously validated it differently)
    percent_samples <- resolve_prevalence(percent_samples, verbose = verbose)

    if (abundance_type == "counts") warn_library_sizes(physeq)

    # Check groups
    sample_data_df <- data.frame(phyloseq::sample_data(physeq), check.names = FALSE)
    if (!group_var %in% names(sample_data_df)) {
        stop(
            "'", group_var, "' is not a sample variable in the phyloseq object. Available: ",
            paste(names(sample_data_df), collapse = ", "), ".",
            call. = FALSE
        )
    }

    groups <- unique(sample_data_df[[group_var]])
    if (!is.null(group_subset)) {
        groups <- intersect(groups, group_subset)
        if (length(groups) < 2) {
            stop("'group_subset' must contain at least 2 valid groups present in the data.", call. = FALSE)
        }
    }

    # --- 2. Identify Core Taxa per Group ---

    if (taxa_rank != "OTU") {
        physeq_glom <- microViz::tax_agg(physeq, rank = taxa_rank)
    } else {
        physeq_glom <- physeq
    }

    core_list <- list()

    for (g in groups) {
        samples_in_group <- rownames(sample_data_df)[sample_data_df[[group_var]] == g]
        if (length(samples_in_group) == 0) next

        physeq_sub <- phyloseq::prune_samples(samples_in_group, physeq_glom)
        physeq_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq_sub) > 0, physeq_sub)

        if (phyloseq::ntaxa(physeq_sub) == 0) next

        prevalence_prop <- taxa_prevalence(
            physeq_sub,
            abundance_threshold = abundance_threshold,
            abundance_type = abundance_type
        )

        core_taxa <- names(prevalence_prop)[prevalence_prop >= percent_samples]
        core_list[[as.character(g)]] <- core_taxa
        if (verbose) {
            message("Group ", g, ": ", length(core_taxa), " core taxa of ",
                    length(prevalence_prop), " present.")
        }
    }

    # --- 3. Intersect Core Taxa (Concordant Core) ---

    if (length(core_list) < 2) {
        stop("Could not identify valid core taxa lists for at least two groups.", call. = FALSE)
    }

    core_sizes <- vapply(core_list, length, integer(1))

    concordant_core <- Reduce(intersect, core_list)
    n_concordant <- length(concordant_core)

    if (verbose) {
        message("Concordant core: ", n_concordant, " taxa shared by all ",
                length(core_list), " evaluated groups.")
    }

    if (n_concordant == 0) {
        # Name the groups with the smallest cores: an empty intersection is
        # usually driven by one restrictive group, and the bare plot title
        # gives no way to tell which.
        empty <- names(core_sizes)[core_sizes == 0]
        message(
            "No concordant core. Core sizes per group: ",
            paste(names(core_sizes), core_sizes, sep = "=", collapse = ", "),
            if (length(empty) > 0) {
                paste0(". Group(s) with no core at all: ", paste(empty, collapse = ", "), ".")
            } else {
                "."
            }
        )
        p <- ggplot() +
             theme_void() +
             ggtitle("No Concordant Core Found")
        attr(p, "concordant_core") <- character(0)
        attr(p, "core_sizes") <- core_sizes
        return(p)
    }

    if (n_concordant < 3) {
        warning("There are fewer than 3 concordant core taxa (", n_concordant,
                "). Plot might look sparse.", call. = FALSE)
    } else if (n_concordant > 8) {
        # Keep the 8 most abundant rather than whichever 8 happened to come
        # first in the taxonomy table
        samples_in_evaluated_groups <- rownames(sample_data_df)[sample_data_df[[group_var]] %in% groups]
        concordant_core <- order_taxa_by_abundance(
            physeq_glom, concordant_core, samples = samples_in_evaluated_groups
        )[1:8]
        warning("There are more than 8 concordant core taxa (", n_concordant,
                "). Keeping the 8 with the highest mean relative abundance across the evaluated groups.",
                call. = FALSE)
        n_concordant <- length(concordant_core)
    }

    # --- 4. Validate Abundance ---

    physeq_rel <- phyloseq::transform_sample_counts(physeq_glom, to_relative)

    # Subsetting to ALL evaluated groups for abundance check
    samples_in_evaluated_groups <- rownames(sample_data_df)[sample_data_df[[group_var]] %in% groups]
    physeq_rel_sub <- phyloseq::prune_samples(samples_in_evaluated_groups, physeq_rel)
    
    # Pruning evaluated core
    core_abundances <- phyloseq::taxa_sums(phyloseq::prune_taxa(concordant_core, physeq_rel_sub)) / length(samples_in_evaluated_groups)
    combined_abundance <- sum(core_abundances, na.rm = TRUE)
    
    if (combined_abundance < 0.05) {
        warning(sprintf("The concordant core taxa make up less than 5%% of the average relative abundance across evaluated groups (%.2f%%). Plot shape may be minimal.", combined_abundance * 100), call. = FALSE)
    }

    # --- 5. Generate Plot ---
    
    # Subset the main object to the evaluated groups
    physeq_plot <- phyloseq::prune_samples(samples_in_evaluated_groups, physeq)

    p <- plot_taxa_star(
        physeq = physeq_plot,
        sample_var = group_var,
        taxa_rank = taxa_rank,
        taxa_names = concordant_core,
        samplecolumn = samplecolumn,
        log_scale = log_scale,
        ...
    )

    attr(p, "concordant_core") <- concordant_core
    attr(p, "core_sizes") <- core_sizes

    return(p)
}
