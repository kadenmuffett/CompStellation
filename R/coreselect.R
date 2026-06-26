#' Select and Plot Core Microbiome
#'
#' This function identifies "core" taxa within treatment groups—defined as taxa present
#' in a specified percentage of samples at a specified relative abundance—and visualizes
#' their presence/absence across all groups.
#'
#' @param physeq A phyloseq object containing your microbiome data.
#' @param group_var A character string. The name of the column in your sample data
#'   to use for grouping samples (e.g., "Treatment", "Location").
#' @param percent_samples A numeric value between 0 and 1 (or 0 and 100).
#'   The percentage of samples within a group that a taxon must be present in to be considered "core".
#'   If > 1, it is treated as a percentage (e.g., 50 means 50%). If <= 1, it is treated as a proportion (e.g., 0.5 means 50%).
#' @param abundance_threshold A numeric value (default 0).
#'   in a sample to be considered "present".
#'
#' @return A ggplot object representing the core taxa presence/absence plot.
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
coreselect <- function(physeq, group_var, percent_samples = 0.95, abundance_threshold = 0) {
    # --- 1. Input Validation ---

    if (!inherits(physeq, "phyloseq")) {
        stop("Error: 'physeq' must be a VALID phyloseq object.")
    }

    if (!group_var %in% phyloseq::sample_variables(physeq)) {
        stop("Error: '", group_var, "' is not a valid sample variable in the phyloseq object.")
    }

    # Normalize percent_samples to proportion 0-1
    if (percent_samples > 1) {
        percent_samples <- percent_samples / 100
    }

    if (percent_samples < 0 || percent_samples > 1) {
        stop("Error: 'percent_samples' must be between 0 and 1 (or 0 and 100).")
    }

    # Check number of groups
    sample_data_df <- data.frame(phyloseq::sample_data(physeq))
    groups <- unique(sample_data_df[[group_var]])
    n_groups <- length(groups)

    if (n_groups < 2) {
        warning("Warning: Factor '", group_var, "' has fewer than 2 groups. Comparison may not be meaningful.")
    } else if (n_groups > 7) {
        warning("Warning: Factor '", group_var, "' has more than 7 groups. Plot may be cluttered.")
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

        # Calculate prevalence for each taxon
        # Prevalence = fraction of samples where abundance > threshold

        # Extract OTU table
        otu_tab <- phyloseq::otu_table(physeq_sub)

        # Ensure taxa are rows for rowSums/rowMeans logic consistency
        if (!phyloseq::taxa_are_rows(otu_tab)) {
            otu_tab <- phyloseq::t(otu_tab)
        }

        # Convert to matrix
        mat <- as.matrix(otu_tab)

        # Binary matrix: 1 if > threshold, 0 otherwise
        present_mat <- mat > abundance_threshold

        # Calculate row sums (number of samples with presence)
        prevalence_counts <- rowSums(present_mat)

        # Calculate proportion
        n_samples_group <- ncol(mat)
        prevalence_prop <- prevalence_counts / n_samples_group

        # Filter core taxa
        core_taxa <- names(prevalence_prop)[prevalence_prop >= percent_samples]

        if (length(core_taxa) > 0) {
            core_list[[as.character(g)]] <- core_taxa
        }
    }

    # --- 3. Combine and Format for Plotting ---

    if (length(core_list) == 0) {
        stop("No core taxa found in any group with the designated thresholds.")
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
            subtitle = paste0("Abundance Threshold > ", abundance_threshold),
            x = "Sample Group",
            y = "Taxon"
        ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8)
        )

    return(p)
}
