#' Plot Matrix of Unique Core Taxa
#'
#' This function identifies the "unique core" taxa for each group (present in that group
#' but not in the core of any others) and generates a matrix of star plots.
#' Each column of the matrix corresponds to the unique core of a specific group.
#' Each row corresponds to a sample type (from the group variable).
#'
#' @param physeq A phyloseq object.
#' @param group_var A character string. The name of the column to use for grouping.
#' @param percent_samples A numeric value for prevalence threshold (0-1).
#' @param abundance_threshold A numeric value for abundance threshold.
#' @param taxa_rank A character string. The taxonomic rank to use (e.g., "Genus").
#' @param samplecolumn A character string. The name of the sample ID column.
#' @param log_scale A logical value. If TRUE, applies a pseudo-log transformation to the y-axis.
#' @param group_subset A character vector specifying which groups to include (if NULL, all are used).
#' @param ... Additional arguments, e.g. "colors_all" passed to `plot_taxa_star`.
#'
#' @return A patchwork object containing the matrix of plots.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#'
#' @export
#'
#' @examples
#' # plot_core_matrix(physeq, "Treatment", percent_samples = 0.5, taxa_rank = "Genus", samplecolumn = "SampleID", log_scale = FALSE)
plot_core_matrix <- function(physeq, group_var, percent_samples, abundance_threshold = 0, taxa_rank, samplecolumn, log_scale = FALSE, group_subset = NULL, ...) {
    # --- 1. Input Validation ---

    if (!inherits(physeq, "phyloseq")) {
        stop("Error: 'physeq' must be a VALID phyloseq object.")
    }

    # Normalize percent_samples
    if (percent_samples > 1) percent_samples <- percent_samples / 100

    # Check groups
    sample_data_df <- data.frame(phyloseq::sample_data(physeq))
    if (!group_var %in% names(sample_data_df)) {
        stop("Error: group_var not found in sample data.")
    }

    groups <- unique(sample_data_df[[group_var]])
    if (!is.null(group_subset)) {
        groups <- intersect(groups, group_subset)
        if (length(groups) < 2) {
            stop("Error: group_subset must contain at least 2 valid groups present in the data.")
        }
    }
    n_groups <- length(groups)

    if (n_groups > 4) {
        stop("Error: This function supports a maximum of 4 groups/treatments.")
    }

    # --- 2. Identify Core Taxa per Group ---

    # Aggregating to rank FIRST is important if we want core at Genus level
    if (taxa_rank != "OTU") {
        physeq_glom <- microViz::tax_agg(physeq, rank = taxa_rank)
    } else {
        physeq_glom <- physeq
    }

    core_list <- list()

    for (g in groups) {
        # Subset
        samples_in_group <- rownames(sample_data_df)[sample_data_df[[group_var]] == g]
        if (length(samples_in_group) == 0) next

        physeq_sub <- phyloseq::prune_samples(samples_in_group, physeq_glom)
        # Prune empty taxa
        physeq_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq_sub) > 0, physeq_sub)

        if (phyloseq::ntaxa(physeq_sub) == 0) next

        # Prevalence
        otu_tab <- phyloseq::otu_table(physeq_sub)
        if (!phyloseq::taxa_are_rows(otu_tab)) otu_tab <- phyloseq::t(otu_tab)

        mat <- as.matrix(otu_tab)
        present_mat <- mat > abundance_threshold
        prevalence_prop <- rowSums(present_mat) / ncol(mat)

        core_taxa <- names(prevalence_prop)[prevalence_prop >= percent_samples]
        core_list[[as.character(g)]] <- core_taxa
        message(paste("Core taxa for group", g, ":", paste(core_taxa, collapse = ", ")))
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

    physeq_rel <- phyloseq::transform_sample_counts(physeq_glom, function(x) {
        s <- sum(x, na.rm = TRUE)
        if (s == 0) return(x) else return(x / s)
    })

    plots <- list()

    for (g in all_groups) {
        u_core <- unique_core_list[[g]]
        n_unique <- length(u_core)

        if (n_unique < 3) {
            warning(paste("Group", g, "has fewer than 3 unique core taxa (", n_unique, "). Plot might look sparse."))
        } else if (n_unique > 8) {
            warning(paste("Group", g, "has more than 8 unique core taxa (", n_unique, "). Plot might be cluttered."))
            u_core <- u_core[1:8]
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
                warning(sprintf("Group %s unique core taxa make up less than 5%% of the average relative abundance within the group (%.2f%%). Plot shape may be minimal.", g, combined_abundance * 100))
            }

            # We need to plot these taxa across ALL sample types (groups).

            p <- plot_taxa_star(
                physeq = physeq,
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

    return(final_plot)
}
