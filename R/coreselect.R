#' Select and Plot Core Microbiome
#'
#' This function identifies "core" taxa within treatment groups—defined as taxa
#' present above an abundance threshold in a specified proportion of the samples
#' in a group—and visualizes their presence/absence across all groups.
#'
#' @section A caution about sequencing depth:
#' With `abundance_type = "counts"` (the default), presence is tested on the
#' abundances as supplied. If library sizes are uneven, a taxon registers as
#' present more often in deeply sequenced samples, so the core partly measures
#' sequencing effort rather than biology. The function warns when library sizes
#' vary more than tenfold. Rarefying beforehand, or using
#' `abundance_type = "relative"` with a non-zero `abundance_threshold`, avoids
#' this.
#'
#' @param physeq A phyloseq object containing your microbiome data.
#' @param group_var A character string. The name of the column in your sample data
#'   to use for grouping samples (e.g., "Treatment", "Location").
#' @param percent_samples The proportion of samples within a group that a taxon
#'   must be present in to be considered "core". Given either as a proportion in
#'   (0, 1] or as a percentage in (1, 100]; the interpretation chosen is
#'   reported, because a bare `1` means 100%, not 1%. Values outside those
#'   ranges are an error.
#' @param abundance_threshold A numeric value (default 0). The abundance a taxon
#'   must **exceed** in a sample for that sample to count as present. The
#'   default of 0 therefore means "any non-zero abundance".
#' @param abundance_type Either "counts" (default; `abundance_threshold` is
#'   applied to the abundances as supplied) or "relative" (each sample is
#'   converted to proportions first, so `abundance_threshold` is a relative
#'   abundance).
#' @param verbose Logical. If TRUE (default), report the resolved thresholds.
#'
#' @return A ggplot object representing the core taxa presence/absence plot.
#'   The core taxa themselves are attached as the "core_list" attribute.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
#' @examples
#' # data(GlobalPatterns)
#' # coreselect(GlobalPatterns, "SampleType", percent_samples = 50, abundance_threshold = 0)
coreselect <- function(physeq, group_var, percent_samples = 0.95, abundance_threshold = 0,
                       abundance_type = c("counts", "relative"), verbose = TRUE) {
    # --- 1. Input Validation ---

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

    abundance_type <- match.arg(abundance_type)

    # Resolve and report the prevalence threshold, so that an ambiguous value
    # such as 1 cannot be silently read as 1% when it means 100%
    percent_samples <- resolve_prevalence(percent_samples, verbose = verbose)

    if (abundance_type == "counts") warn_library_sizes(physeq)

    # Check number of groups
    sample_data_df <- data.frame(phyloseq::sample_data(physeq))
    groups <- unique(sample_data_df[[group_var]])
    n_groups <- length(groups)

    if (n_groups < 2) {
        warning("Factor '", group_var, "' has fewer than 2 groups. Comparison may not be meaningful.", call. = FALSE)
    } else if (n_groups > 7) {
        warning("Factor '", group_var, "' has more than 7 groups. Plot may be cluttered.", call. = FALSE)
    }

    # --- 2. Identify Core Taxa per Group ---

    core_list <- list()

    # Get the OTU table. We need samples as columns for easy subsetting logic or standard approaches.
    # phyloseq::otu_table always allows access.
    # Let's use standard phyloseq tools.

    for (g in groups) {
        # Subset to current group
        # We must quote the group value if it's a string, which it likely is.
        # Prune samples strictly.
        # We can use formatted strings for subset_samples, but logic is tricky with variable names.
        # Easier to indices.

        # Identify samples in this group
        samples_in_group <- rownames(sample_data_df)[sample_data_df[[group_var]] == g]

        if (length(samples_in_group) == 0) next

        # Prune physeq object to just these samples
        # We use prune_samples.
        physeq_sub <- phyloseq::prune_samples(samples_in_group, physeq)

        # Remove taxa not present in any sample of this subset (optional, but speeds up calc)
        physeq_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq_sub) > 0, physeq_sub)

        if (phyloseq::ntaxa(physeq_sub) == 0) {
            next
        }

        # Prevalence = fraction of samples in the group where the abundance
        # exceeds the threshold. Shared with plot_core_matrix() and
        # plot_concordant_core() so the three cannot drift apart.
        prevalence_prop <- taxa_prevalence(
            physeq_sub,
            abundance_threshold = abundance_threshold,
            abundance_type = abundance_type
        )

        # Filter core taxa
        core_taxa <- names(prevalence_prop)[prevalence_prop >= percent_samples]

        if (length(core_taxa) > 0) {
            core_list[[as.character(g)]] <- core_taxa
        }
    }

    # --- 3. Combine and Format for Plotting ---

    if (length(core_list) == 0) {
        stop("No core taxa found in any group with the designated thresholds.", call. = FALSE)
    }

    # Get all unique core taxa across all groups
    all_core_taxa <- unique(unlist(core_list))

    # Create a data frame for plotting
    # We want a grid of (Taxon, Group) -> Present/Absent

    plot_data <- expand.grid(
        Taxon = all_core_taxa,
        Group = as.character(groups),
        stringsAsFactors = FALSE
    )

    # Determine presence
    plot_data$IsCore <- mapply(function(t, g) {
        if (g %in% names(core_list)) {
            t %in% core_list[[g]]
        } else {
            FALSE
        }
    }, plot_data$Taxon, plot_data$Group)

    # Filter to only show points where IsCore lies?
    # Request: "points representing presense absence".
    # Usually this means a dot if present, nothing if absent.
    # So we keep the TRUE ones for geom_point.
    # But we might want the grid to be visible? usually just dots.

    plot_data_present <- plot_data[plot_data$IsCore, ]

    # --- 4. Plotting ---

    p <- ggplot(plot_data_present, aes(x = Group, y = Taxon)) +
        geom_point(size = 3) +
        theme_bw() + # Clean theme
        labs(
            title = paste0("Core Microbiome (Prevalence >= ", percent_samples * 100, "%)"),
            subtitle = paste0(
                "Present = ", abundance_type, " abundance > ", abundance_threshold,
                "; n = ", paste(vapply(
                    as.character(groups),
                    function(g) sum(sample_data_df[[group_var]] == g, na.rm = TRUE),
                    numeric(1)
                ), collapse = "/"), " samples per group"
            ),
            x = "Sample Group",
            y = "Taxon"
        ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8)
        )

    # Return the memberships alongside the plot: previously the only way to see
    # which taxa were core was to read them off the y-axis.
    attr(p, "core_list") <- core_list

    return(p)
}
