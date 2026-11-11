# Placeholder labels that mean "no real classification" (matched case-insensitively,
# after trimming whitespace).THIS MAY NOT BE EVERYTHING. CHECK YOUR TAXA.
unclassified_pattern <- paste0(
  "^(",
  "na|n/a|nan|null|none|nil|",
  "(unclassified|unknown|unidentified|unassigned|uncultured|uncharacterized)([ _.-].*)?|",
  "ambiguous[ _]taxa|(.*[ _])?incertae[ _]sedis([ _].*)?|.*metagenome|",
  "([a-z]|d_[0-9]+)__|",   # bare rank prefixes such as "g__" or "D_5__"
  "[-_?.]*",               # "", "-", "?", "__"
  ")$"
)

#' Aggregate a phyloseq object to a taxonomic rank
#'
#' Sums counts of all taxa that share a lineage down to `rank`. Placeholder
#' labels (NA, "NULL", "unclassified", "g__", ...) are treated as missing and
#' replaced with the nearest classified higher rank plus RANK NAME, e.g.
#' "Bacteroidaceae Family".
#' Taxa names in the result are the labels at `rank`.
#'
#' @param physeq A phyloseq object with a tax_table.
#' @param rank Name of the rank to aggregate to.
#' @return A phyloseq object with one taxon per label at `rank`. The phylogenetic
#'   tree and reference sequences are dropped, since they no longer apply.
#' @keywords internal
#' @noRd
agg_taxa <- function(physeq, rank) {
  tt <- methods::as(phyloseq::tax_table(physeq), "matrix")
  idx <- match(rank, colnames(tt))
  if (is.na(idx)) {
    stop("'", rank, "' is not a rank in the tax_table.", call. = FALSE)
  }
  tt <- tt[, seq_len(idx), drop = FALSE]

  # Normalise placeholders to NA
  tt[] <- trimws(tt)
  tt[is.na(tt) | grepl(unclassified_pattern, tt, ignore.case = TRUE)] <- NA

  # Fill gaps from the nearest classified higher rank: "Bacteroidaceae Family"
  last_known <- rep("Unclassified", nrow(tt))
  for (j in seq_len(ncol(tt))) {
    miss <- is.na(tt[, j])
    last_known[!miss] <- paste(tt[!miss, j], colnames(tt)[j])
    tt[miss, j] <- last_known[miss]
  }

  # Labels must be unique per lineage: a name shared by different lineages
  # (e.g. a genus listed under two families) gets its parent appended.
  lineage <- apply(tt, 1, paste, collapse = "|")
  first <- !duplicated(lineage)
  lin_label <- stats::setNames(tt[first, idx], lineage[first])
  clash <- duplicated(lin_label) | duplicated(lin_label, fromLast = TRUE)
  if (any(clash) && idx > 1) {
    lin_label[clash] <- paste0(lin_label[clash], " (", tt[first, idx - 1][clash], ")")
  }
  lin_label[] <- make.unique(lin_label, sep = " ")
  label <- unname(lin_label[lineage])
  tt[, idx] <- label  # keep the rank column in step with taxa_names for psmelt()

  # Sum counts per label
  otu <- methods::as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) otu <- t(otu)
  otu_agg <- rowsum(otu, group = label, reorder = FALSE)

  tt_agg <- tt[!duplicated(label), , drop = FALSE]
  rownames(tt_agg) <- label[!duplicated(label)]
  tt_agg <- tt_agg[rownames(otu_agg), , drop = FALSE]

  parts <- list(
    phyloseq::otu_table(otu_agg, taxa_are_rows = TRUE),
    phyloseq::tax_table(tt_agg)
  )
  sd <- phyloseq::sample_data(physeq, errorIfNULL = FALSE)
  if (!is.null(sd)) parts <- c(parts, list(sd))
  message("Taxa aggregated to rank: ", rank, " using an internal heuristic to remove unclassifieds. Some unclassifieds may remain.")
  do.call(phyloseq::phyloseq, parts)
}
