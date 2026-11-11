#' Upside-Down Jellyfish (*Cassiopea xamachana*) and Environmental Microbiome Dataset
#'
#' A `phyloseq` object containing 16S rRNA gene amplicon sequencing counts for
#' the top 3,000 most abundant amplicon sequence variants (ASVs) across 76
#' host-associated and environmental samples collected from 8 sites in the
#' Florida Keys.
#'
#' @format A `phyloseq` object with 3,000 taxa and 76 samples:
#' \describe{
#'   \item{`otu_table`}{3,000 ASVs (rows) by 76 samples (columns) of sequencing read counts.}
#'   \item{`tax_table`}{3,000 ASVs by 7 taxonomic ranks (`Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species`).}
#'   \item{`sample_data`}{76 samples by 21 sample metadata variables, including:
#'     \itemize{
#'       \item `Sample_NAME`: Unique sample identifier (character).
#'       \item `Location`: Sampling location code in the Florida Keys (`BPK`, `CK`, `GB`, `HHP`, `KL`, `LW`, `MK`, `VP`).
#'       \item `Body_cav`: Sample compartment / habitat (`Bell` = jellyfish bell, `GVC` = gastrovascular cavity, `Sub` = benthic substrate, `WTR` = surrounding water column).
#'       \item `Identity`, `Specimen`: Specimen/replicate identifiers.
#'       \item `Buffer`: Preservation buffer (`DESS`).
#'       \item `Species_HOST`, `Species_O`: Host species code (`CX` = *Cassiopea xamachana*, `NA` for environmental substrate/water controls).
#'       \item `Temp`, `Temperature_.C.`, `Salinity`, `pH`, `Depth..m.`, `Lat`, `Long`, `Date`, `Time`: Environmental and geographic measurements.
#'       \item `Diameter..cm.`, `OAL..cm.`, `Concentration..indv.m2.`: Host morphometric and density measurements.
#'       \item `Neg`: Logical indicator for negative control samples (`FALSE`).
#'     }
#'   }
#' }
#'
#' @usage data("cassiopea_ps")
#'
#' @examples
#' data("cassiopea_ps")
#' cassiopea_ps
"cassiopea_ps"
