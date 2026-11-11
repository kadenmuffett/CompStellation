#' Convert a vegan-Style Data Frame or Matrix to a phyloseq Object
#'
#' Converts a sample-by-taxon community table (`vegan` style, where rows are
#' samples and columns are taxa) into a [`phyloseq::phyloseq`] object ready for
#' all `CompStellation` functions. Sample metadata can be supplied either as a
#' separate `metadata` data frame or as columns inside `df`.
#'
#' @param df A `data.frame` or `matrix` with samples as rows and taxa as columns.
#'   May also include sample metadata columns if `metadata` is `NULL`.
#' @param metadata Optional `data.frame` of sample metadata (rows = samples).
#'   If `NULL`, metadata columns are extracted from `df` using `metadata_cols`
#'   (or all non-numeric columns plus `sample_col`).
#' @param metadata_cols Optional character vector of column names in `df` to
#'   treat as sample metadata when `metadata` is `NULL`. Useful when some
#'   metadata columns are numeric (e.g., `"pH"`, `"Depth"`).
#' @param sample_col Optional character string giving a column name (in
#'   `metadata` or `df`) that holds unique sample IDs. If supplied, it is used
#'   as the row names of the resulting `phyloseq` object and retained in
#'   `sample_data` so it can be passed as `samplecolumn` to plotting functions.
#'   If `NULL` and no sample ID column exists, a `"SampleID"` column is added
#'   from the row names.
#' @param tax_table Optional `data.frame` or `matrix` of taxonomic assignments
#'   (rows = taxa matching the taxon columns of `df`, columns = taxonomic ranks
#'   such as `"Phylum"`, `"Family"`, `"Genus"`). If `NULL` (default), a
#'   single-column taxonomy table named `taxa_rank` is created from the taxon
#'   column names.
#' @param taxa_rank Character string naming the taxonomic rank column when
#'   `tax_table` is `NULL` (default `"OTU"`).
#'
#' @return A [`phyloseq::phyloseq`] object containing `otu_table`,
#'   `sample_data`, and `tax_table`.
#'
#' @import phyloseq
#' @export
#'
#' @examples
#' # 1. Single combined wide data.frame (samples as rows)
#' wide_df <- data.frame(
#'   SampleID = paste0("S", 1:6),
#'   Habitat  = rep(c("Reef", "Mangrove"), each = 3),
#'   TaxonA   = c(120, 95, 140, 10, 5, 12),
#'   TaxonB   = c(15, 20, 10, 110, 130, 90),
#'   TaxonC   = c(45, 50, 40, 45, 60, 55),
#'   TaxonD   = c(5, 8, 12, 30, 25, 40)
#' )
#' ps <- vegan_to_phyloseq(wide_df, sample_col = "SampleID", taxa_rank = "Genus")
#'
#' # 2. Separate vegan-style community matrix and metadata data.frame
#' comm_mat <- wide_df[, c("TaxonA", "TaxonB", "TaxonC", "TaxonD")]
#' meta_df  <- wide_df[, c("SampleID", "Habitat")]
#' ps2 <- vegan_to_phyloseq(comm_mat, metadata = meta_df, sample_col = "SampleID")
vegan_to_phyloseq <- function(df,
                              metadata = NULL,
                              metadata_cols = NULL,
                              sample_col = NULL,
                              tax_table = NULL,
                              taxa_rank = "OTU") {
  if (!is.data.frame(df) && !is.matrix(df)) {
    stop("'df' must be a data.frame or matrix with samples as rows.", call. = FALSE)
  }
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)

  if (!is.null(metadata)) {
    if (!is.data.frame(metadata) && !is.matrix(metadata)) {
      stop("'metadata' must be a data.frame or matrix.", call. = FALSE)
    }
    sample_dat <- as.data.frame(metadata, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(sample_dat) != nrow(df)) {
      stop("'df' and 'metadata' must have the same number of rows (samples).", call. = FALSE)
    }
    taxa_cols <- setdiff(names(df), c(sample_col, metadata_cols))
    taxa_df <- df[, taxa_cols, drop = FALSE]
  } else {
    if (!is.null(metadata_cols)) {
      missing_meta <- setdiff(metadata_cols, names(df))
      if (length(missing_meta) > 0) {
        stop("Columns in 'metadata_cols' not found in 'df': ",
             paste(missing_meta, collapse = ", "), ".", call. = FALSE)
      }
      meta_names <- unique(c(sample_col, metadata_cols))
    } else {
      is_num <- vapply(df, is.numeric, logical(1))
      meta_names <- unique(c(sample_col, names(df)[!is_num]))
    }
    taxa_cols <- setdiff(names(df), meta_names)
    if (length(taxa_cols) == 0) {
      stop("No numeric taxon columns found in 'df'.", call. = FALSE)
    }
    taxa_df <- df[, taxa_cols, drop = FALSE]
    sample_dat <- df[, meta_names, drop = FALSE]
  }

  if (!all(vapply(taxa_df, is.numeric, logical(1)))) {
    stop("All taxon columns must be numeric. Specify non-taxon columns via 'metadata' or 'metadata_cols'.",
         call. = FALSE)
  }

  # Resolve sample IDs and row names
  if (!is.null(sample_col)) {
    if (!sample_col %in% names(sample_dat)) {
      if (sample_col %in% names(df)) {
        sample_dat[[sample_col]] <- df[[sample_col]]
      } else {
        stop("'sample_col' (\"", sample_col, "\") not found in 'metadata' or 'df'.", call. = FALSE)
      }
    }
    sample_ids <- as.character(sample_dat[[sample_col]])
    if (anyDuplicated(sample_ids)) {
      stop("Sample IDs in '", sample_col, "' must be unique.", call. = FALSE)
    }
    rownames(taxa_df) <- sample_ids
    rownames(sample_dat) <- sample_ids
  } else {
    if (is.null(rownames(taxa_df)) || all(rownames(taxa_df) == as.character(seq_len(nrow(taxa_df))))) {
      rownames(taxa_df) <- paste0("Sample_", seq_len(nrow(taxa_df)))
    }
    rownames(sample_dat) <- rownames(taxa_df)
    if (ncol(sample_dat) == 0 || !"SampleID" %in% names(sample_dat)) {
      sample_dat$SampleID <- rownames(taxa_df)
    }
  }

  taxa_mat <- as.matrix(taxa_df)
  otu <- phyloseq::otu_table(taxa_mat, taxa_are_rows = FALSE)
  sdata <- phyloseq::sample_data(sample_dat)

  if (!is.null(tax_table)) {
    tax_mat <- as.matrix(tax_table)
    if (is.null(rownames(tax_mat))) {
      if (nrow(tax_mat) == ncol(taxa_mat)) {
        rownames(tax_mat) <- colnames(taxa_mat)
      } else {
        stop("'tax_table' must have rownames matching the taxon column names of 'df'.", call. = FALSE)
      }
    }
  } else {
    tax_mat <- matrix(colnames(taxa_mat), ncol = 1,
                      dimnames = list(colnames(taxa_mat), taxa_rank))
  }
  tax <- phyloseq::tax_table(tax_mat)

  phyloseq::phyloseq(otu, sdata, tax)
}
