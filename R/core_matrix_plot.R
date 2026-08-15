#' Plot Matrix of Unique Core Taxa
#'
#' This function identifies the "unique core" taxa for each group (present in that group
#' but not in the core of any others) and generates a matrix of star plots.
#' Each column of the matrix corresponds to the unique core of a specific group.
#' Each row corresponds to a sample type (from the group variable).
#'
#' @section A caution about thresholds:
#' The unique core is a hard set difference, so a taxon at 0.849 prevalence in
#' one group and 0.851 in another is declared unique to the second. The
#' prevalence of every unique-core taxon in every group is returned as the
#' "prevalence" attribute of the result so the margin can be checked. See
#' [coreselect()] for the caution about sequencing depth, which applies here too.
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
#' @param ... Additional arguments, e.g. "colors_all" passed to `plot_taxa_star`.
#'
#' @return A patchwork object containing the matrix of plots, with the unique
#'   core memberships attached as the "unique_core" attribute and the per-group
#'   prevalence of those taxa as the "prevalence" attribute.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#'
#' @export
#'
#' @examples
#' # plot_core_matrix(physeq, "Treatment", percent_samples = 0.5,
#' #                  taxa_rank = "Genus", samplecolumn = "SampleID")
plot_core_matrix <- function(physeq, group_var, percent_samples, abundance_threshold = 0,
                             abundance_type = c("counts", "relative"), taxa_rank, samplecolumn,
                             log_scale = FALSE, group_subset = NULL, verbose = TRUE, ...) {
    # --- 1. Input Validation ---

    if (!inherits(physeq, "phyloseq")) {
        stop("'physeq' must be a VALID phyloseq object.", call. = FALSE)
    }

    abundance_type <- match.arg(abundance_type)

    # Resolve and report the prevalence threshold (shared with coreselect() and
    # plot_concordant_core(), which previously validated it differently)
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
    n_groups <- length(groups)

    if (n_groups > 4) {
        stop("This function supports a maximum of 4 groups/treatments.", call. = FALSE)
    }

    # Samples belonging to the groups actually being evaluated
    samples_evaluated <- rownames(sample_data_df)[sample_data_df[[group_var]] %in% groups]

    # --- 2. Identify Core Taxa per Group ---

    # Aggregating to rank FIRST is important if we want core at Genus level
    if (taxa_rank != "OTU") {
        physeq_glom <- microViz::tax_agg(physeq, rank = taxa_rank)
    } else {
        physeq_glom <- physeq
    }

    core_list <- list()
    prevalence_list <- list()

    for (g in groups) {
        # Subset
        samples_in_group <- rownames(sample_data_df)[sample_data_df[[group_var]] == g]
        if (length(samples_in_group) == 0) next

        physeq_sub <- phyloseq::prune_samples(samples_in_group, physeq_glom)
        # Prune empty taxa
        physeq_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq_sub) > 0, physeq_sub)

        if (phyloseq::ntaxa(physeq_sub) == 0) next

        prevalence_prop <- taxa_prevalence(
            physeq_sub,
            abundance_threshold = abundance_threshold,
            abundance_type = abundance_type
        )
        prevalence_list[[as.character(g)]] <- prevalence_prop

        core_taxa <- names(prevalence_prop)[prevalence_prop >= percent_samples]
        core_list[[as.character(g)]] <- core_taxa
        if (verbose) {
            message("Group ", g, ": ", length(core_taxa), " core taxa of ",
                    length(prevalence_prop), " present.")
        }
    }

    # --- 3. Identify Unique Core Taxa ---

    unique_core_list <- list()
    all_groups <- names(core_list)

    for (g in all_groups) {
        my_core <- core_list[[g]]
        other_groups <- setdiff(all_groups, g)

        # Union of all other cores
        other_cores <- unique(unlist(core_list[other_groups]))

        # Unique core
        unique_core <- setdiff(my_core, other_cores)
        unique_core_list[[g]] <- unique_core
    }

    # --- 4. Validation of Unique Core Size & Abundance ---

    physeq_rel <- phyloseq::transform_sample_counts(physeq_glom, to_relative)

    plots <- list()

    for (g in all_groups) {
        u_core <- unique_core_list[[g]]
        n_unique <- length(u_core)

        if (n_unique < 3) {
            warning("Group ", g, " has fewer than 3 unique core taxa (", n_unique,
                    "). Plot might look sparse.", call. = FALSE)
        } else if (n_unique > 8) {
            # Keep the 8 most abundant rather than whichever 8 happened to come
            # first in the taxonomy table
            samples_in_g <- rownames(sample_data_df)[sample_data_df[[group_var]] == g]
            u_core <- order_taxa_by_abundance(physeq_glom, u_core, samples = samples_in_g)[1:8]
            warning("Group ", g, " has more than 8 unique core taxa (", n_unique,
                    "). Keeping the 8 with the highest mean relative abundance in this group.",
                    call. = FALSE)
            n_unique <- length(u_core)
        }

        # --- 5. Generate Plots ---

        if (n_unique > 0) {
            # Check combined abundance
            samples_in_g <- rownames(sample_data_df)[sample_data_df[[group_var]] == g]
            physeq_rel_sub <- phyloseq::prune_samples(samples_in_g, physeq_rel)
            core_abundances <- phyloseq::taxa_sums(phyloseq::prune_taxa(u_core, physeq_rel_sub)) / length(samples_in_g)
            combined_abundance <- sum(core_abundances, na.rm = TRUE)
            
            if (combined_abundance < 0.05) {
                warning(sprintf("Group %s unique core taxa make up less than 5%% of the average relative abundance within the group (%.2f%%). Plot shape may be minimal.", g, combined_abundance * 100), call. = FALSE)
            }

            # Plot these taxa across the evaluated groups only. Passing the full
            # object here would facet in groups whose cores were never compared,
            # so the "unique" claim would not hold for those panels.
            p <- plot_taxa_star(
                physeq = phyloseq::prune_samples(samples_evaluated, physeq),
                sample_var = group_var,
                taxa_rank = taxa_rank,
                taxa_names = u_core,
                samplecolumn = samplecolumn,
                log_scale = log_scale,
                ...
            )

            # Customize Plot
            p <- p +
                facet_wrap(as.formula(paste("~", group_var)), ncol = 1, strip.position = "right") +
                ggtitle(paste("Unique Core:", g)) +
                theme(legend.position = "none")

            plots[[g]] <- p
        } else {
            plots[[g]] <- ggplot() +
                theme_void() +
                ggtitle(paste("No Unique Core:", g))
        }
    }
    # --- 6. Global Scaling ---

    # Calculate global maximum value for radial axis across all generated star plots
    global_max_y <- 0
    for (p in plots) {
        if (inherits(p, "ggplot") && !is.null(p$data) && nrow(p$data) > 0) {
            pb <- ggplot2::ggplot_build(p)
            for (layer_data in pb$data) {
                if ("ymax" %in% names(layer_data)) {
                    global_max_y <- max(global_max_y, max(layer_data$ymax, na.rm = TRUE), na.rm = TRUE)
                }
                if ("y" %in% names(layer_data)) {
                    global_max_y <- max(global_max_y, max(layer_data$y, na.rm = TRUE), na.rm = TRUE)
                }
            }
        }
    }

    # Apply global maximum limit to all valid plots to ensure scales are consistent
    if (global_max_y > 0) {
        for (g in names(plots)) {
            if (inherits(plots[[g]], "ggplot") && !is.null(plots[[g]]$data) && nrow(plots[[g]]$data) > 0) {
                if (log_scale) {
                    plots[[g]] <- plots[[g]] + ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.01), limits = c(0, global_max_y))
                } else {
                    plots[[g]] <- plots[[g]] + ggplot2::scale_y_continuous(limits = c(0, global_max_y))
                }
            }
        }
    }

    # --- 7. Combine Plots ---

    final_plot <- patchwork::wrap_plots(plots, nrow = 1)

    # Return the memberships and the prevalences behind them, so the margin on
    # each "unique" call can be checked rather than taken on trust.
    attr(final_plot, "unique_core") <- unique_core_list
    attr(final_plot, "prevalence") <- prevalence_list

    return(final_plot)
}
