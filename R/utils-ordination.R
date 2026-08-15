# Shared helpers for the ordination-based star plots (plot_ordi_star,
# plot_ordi_star_iris). Both functions previously carried their own copy of the
# score-extraction and eigenvalue logic.

#' Ordination methods accepted by the star plot functions
#'
#' These are the methods `phyloseq::ordinate()` understands. Note that "PCA" is
#' not among them: an unconstrained RDA is the phyloseq spelling of a PCA.
#'
#' @keywords internal
ordi_methods <- function() {
  c("PCoA", "MDS", "NMDS", "DPCoA", "DCA", "CCA", "RDA", "CAP")
}

#' Percent of variation explained by each ordination axis
#'
#' Bray-Curtis and other semimetric dissimilarities are not always Euclidean, so
#' `ape::pcoa()` can return negative eigenvalues (whether it does depends on the
#' data). When it does, dividing by the sum of *all* eigenvalues understates the
#' denominator, because the negatives subtract from it, and so overstates every
#' positive axis: on a test case with 5 of 29 eigenvalues negative, axis 1 reads
#' 11.14% that way against 9.30% by the corrected convention. Both conventions
#' total 100%; the difference is how the total is divided up.
#'
#' This prefers the corrected relative eigenvalues that `ape` supplies exactly
#' when a correction is warranted (`Rel_corr_eig`), and otherwise restricts the
#' denominator to the positive eigenvalues. Where there are no negative
#' eigenvalues it returns the same numbers as the plain relative eigenvalues.
#'
#' @param ord An ordination object.
#'
#' @return `NULL`, or a list with `pct` (numeric, percentages, `NA` for
#'   negative eigenvalues) and `note` (a character description of the
#'   convention used).
#' @keywords internal
axis_percent_var <- function(ord) {
  vals <- ord$values
  if (is.null(vals)) {
    return(NULL)
  }

  # ape::pcoa() supplies this exactly when a Lingoes/Cailliez correction applies
  if (!is.null(vals$Rel_corr_eig)) {
    return(list(
      pct = 100 * vals$Rel_corr_eig,
      note = "corrected relative eigenvalues"
    ))
  }

  eig <- vals$Eigenvalues
  if (is.null(eig)) {
    if (!is.null(vals$Relative_eig)) {
      return(list(pct = 100 * vals$Relative_eig, note = "relative eigenvalues"))
    }
    return(NULL)
  }

  if (any(eig < 0, na.rm = TRUE)) {
    pos <- !is.na(eig) & eig > 0
    pct <- rep(NA_real_, length(eig))
    pct[pos] <- 100 * eig[pos] / sum(eig[pos])
    return(list(
      pct = pct,
      note = paste0(
        "positive eigenvalues only (", sum(eig < 0, na.rm = TRUE),
        " eigenvalues are negative because the distance is not Euclidean)"
      )
    ))
  }

  list(pct = 100 * eig / sum(eig), note = "relative eigenvalues")
}

#' Extract site scores and axis labels from an ordination object
#'
#' @param ord An ordination object.
#' @param n_axes The number of axes wanted.
#' @param verbose Logical. If TRUE, report the eigenvalue convention used.
#'
#' @return A list with `scores` (a matrix of site scores, samples in rows),
#'   `labels` (axis labels including percent variation where available) and
#'   `n_axes` (the number of axes actually available).
#' @keywords internal
extract_ordination_axes <- function(ord, n_axes, verbose = TRUE) {
  percent_var <- NULL

  if (inherits(ord, c("pcoa", "dpcoa", "list")) && !is.null(ord$vectors)) {
    vecs <- ord$vectors
    pv <- axis_percent_var(ord)
    if (!is.null(pv)) {
      percent_var <- pv$pct
      if (verbose) message("Axis percentages use ", pv$note, ".")
    }
    axis_raw_names <- paste0("Axis ", seq_len(ncol(vecs)))
  } else if (inherits(ord, c("metaMDS", "monoMDS"))) {
    vecs <- ord$points
    axis_raw_names <- paste0("NMDS", seq_len(ncol(vecs)))
  } else if (inherits(ord, c("rda", "cca"))) {
    # scores() defaults to two axes, which silently truncated n_axes before
    vecs <- vegan::scores(ord, display = "sites", choices = seq_len(n_axes))
    vecs <- as.matrix(vecs)
    # Constrained axes live in $CCA; $CA holds the unconstrained residual axes
    eig <- if (!is.null(ord$CCA) && length(ord$CCA$eig) > 0) {
      ord$CCA$eig
    } else {
      ord$CA$eig
    }
    if (!is.null(eig)) {
      denom <- if (!is.null(ord$tot.chi)) ord$tot.chi else sum(eig)
      percent_var <- 100 * eig / denom
      if (verbose) message("Axis percentages are given as a share of total inertia.")
    }
    axis_raw_names <- colnames(vecs)
    if (is.null(axis_raw_names)) axis_raw_names <- paste0("Axis ", seq_len(ncol(vecs)))
  } else if (inherits(ord, "prcomp")) {
    vecs <- ord$x
    if (!is.null(ord$sdev)) percent_var <- 100 * ord$sdev^2 / sum(ord$sdev^2)
    axis_raw_names <- paste0("PC", seq_len(ncol(vecs)))
  } else {
    vecs <- tryCatch(
      as.matrix(vegan::scores(ord, display = "sites", choices = seq_len(n_axes))),
      error = function(e) ord$vectors
    )
    if (is.null(vecs)) {
      stop("Could not extract site scores from the supplied 'ord' object.", call. = FALSE)
    }
    axis_raw_names <- colnames(vecs)
    if (is.null(axis_raw_names)) axis_raw_names <- paste0("Axis ", seq_len(ncol(vecs)))
  }

  available <- min(n_axes, ncol(vecs))
  idx <- seq_len(available)

  if (!is.null(percent_var) && length(percent_var) >= available) {
    pct <- percent_var[idx]
    labels <- ifelse(
      is.na(pct),
      axis_raw_names[idx],
      paste0(axis_raw_names[idx], "\n(", sprintf("%.1f", pct), "%)")
    )
  } else {
    labels <- axis_raw_names[idx]
  }

  list(
    scores = vecs[, idx, drop = FALSE],
    labels = labels,
    n_axes = available
  )
}

#' Align sample data to the row order of an ordination's site scores
#'
#' A positional `cbind()` silently pairs each sample's coordinates with another
#' sample's metadata if the ordination dropped or reordered any samples, which a
#' user-supplied `ord` object may well have done. This matches on sample names
#' where they are available and errors where they cannot be reconciled.
#'
#' @param scores A matrix of site scores, ideally with sample names as rownames.
#' @param physeq The phyloseq object the ordination came from.
#'
#' @return The sample data as a data frame, reordered to match the score rows.
#' @keywords internal
align_sample_data <- function(scores, physeq) {
  sample_dat <- data.frame(phyloseq::sample_data(physeq), check.names = FALSE)
  score_names <- rownames(scores)

  if (is.null(score_names)) {
    if (nrow(scores) != nrow(sample_dat)) {
      stop(
        "The ordination returned ", nrow(scores), " rows for ", nrow(sample_dat),
        " samples, and has no sample names to match on.",
        call. = FALSE
      )
    }
    return(sample_dat)
  }

  absent <- setdiff(score_names, rownames(sample_dat))
  if (length(absent) > 0) {
    stop(
      "The ordination contains samples that are not in the phyloseq object: ",
      paste(utils::head(absent, 5), collapse = ", "),
      if (length(absent) > 5) ", ..." else "", ".",
      call. = FALSE
    )
  }

  sample_dat[score_names, , drop = FALSE]
}

#' Shift ordination scores onto a plottable radius
#'
#' A radar plot's radius cannot represent negative values, so all axes are
#' shifted by a single constant. The shift is reported so it can be undone in
#' the axis labels.
#'
#' @param x A data frame or matrix of scores.
#'
#' @return The offset to add (0 if all values are already non-negative).
#' @keywords internal
radial_offset <- function(x) {
  min_val <- min(as.matrix(x), na.rm = TRUE)
  if (!is.finite(min_val) || min_val >= 0) 0 else abs(min_val) * 1.1
}
