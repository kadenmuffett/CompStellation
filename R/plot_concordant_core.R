#' Plot Concordant Core Taxa
#'
#' This function identifies the concordant (overlapping) core taxa across all evaluated groups
#' and plots their relative abundance using `plot_taxa_star`.
#'
#' @param physeq A phyloseq object.
#' @param group_var A character string. The name of the column to use for grouping.
#' @param percent_samples A numeric value for prevalence threshold (0-1).
#' @param abundance_threshold A numeric value for abundance threshold.
#' @param taxa_rank A character string. The taxonomic rank to use (e.g., "Genus").
#' @param samplecolumn A character string. The name of the sample ID column.
#' @param log_scale A logical value. If TRUE, applies a pseudo-log transformation to the y-axis.
#' @param group_subset A character vector specifying which groups to include (if NULL, all are used).
#' @param ... Additional arguments passed to `plot_taxa_star`, e.g. "colors_all".
#'
#' @return A ggplot object representing the star plot of concordant core taxa.
#'
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
#' @examples
#' # plot_concordant_core(physeq, "Treatment", percent_samples = 0.5, taxa_rank = "Genus", samplecolumn = "SampleID", log_scale = FALSE)
plot_concordant_core <- function(physeq, group_var, percent_samples, abundance_threshold = 0, taxa_rank, samplecolumn, log_scale = FALSE, group_subset = NULL, ...) {
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

        otu_tab <- phyloseq::otu_table(physeq_sub)
        if (!phyloseq::taxa_are_rows(otu_tab)) otu_tab <- phyloseq::t(otu_tab)

        mat <- as.matrix(otu_tab)
        present_mat <- mat > abundance_threshold
        prevalence_prop <- rowSums(present_mat) / ncol(mat)

        core_taxa <- names(prevalence_prop)[prevalence_prop >= percent_samples]
        core_list[[as.character(g)]] <- core_taxa
        message(paste("Core taxa for group", g, ":", paste(core_taxa, collapse = ", ")))
    }

    # --- 3. Intersect Core Taxa (Concordant Core) ---

    if (length(core_list) < 2) {
        stop("Error: Could not identify valid core taxa lists for at least two groups.")
    }

    concordant_core <- Reduce(intersect, core_list)
    n_concordant <- length(concordant_core)

    message("Concordant core taxa across evaluated groups: ", paste(concordant_core, collapse = ", "))

    if (n_concordant == 0) {
        p <- ggplot() +
             theme_void() +
             ggtitle("No Concordant Core Found")
        return(p)
    }

    if (n_concordant < 3) {
        warning(paste("There are fewer than 3 concordant core taxa (", n_concordant, "). Plot might look sparse."))
    } else if (n_concordant > 8) {
        warning(paste("There are more than 8 concordant core taxa (", n_concordant, "). Plot might be cluttered, limiting to top 8."))
        concordant_core <- concordant_core[1:8]
        n_concordant <- length(concordant_core)
    }

    # --- 4. Validate Abundance ---

    physeq_rel <- phyloseq::transform_sample_counts(physeq_glom, function(x) {
        s <- sum(x, na.rm = TRUE)
        if (s == 0) return(x) else return(x / s)
    })

    # Subsetting to ALL evaluated groups for abundance check
    samples_in_evaluated_groups <- rownames(sample_data_df)[sample_data_df[[group_var]] %in% groups]
    physeq_rel_sub <- phyloseq::prune_samples(samples_in_evaluated_groups, physeq_rel)
    
    # Pruning evaluated core
    core_abundances <- phyloseq::taxa_sums(phyloseq::prune_taxa(concordant_core, physeq_rel_sub)) / length(samples_in_evaluated_groups)
    combined_abundance <- sum(core_abundances, na.rm = TRUE)
    
    if (combined_abundance < 0.05) {
        warning(sprintf("The concordant core taxa make up less than 5%% of the average relative abundance across evaluated groups (%.2f%%). Plot shape may be minimal.", combined_abundance * 100))
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
    
    return(p)
}
