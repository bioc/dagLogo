#' @rdname dagLogo-package
#' @aliases dagLogo-package
#' @docType package
#' @title Visualize significant conserved amino acid sequence pattern in groups
#' based on probability theory
#' @description dagLogo provides differential analysis of grouped/ungrouped amino acid usage
#' between an input set of aligned peptide sequences and a background set of 
#' aligned peptide sequences which can be generated in different ways. Results 
#' of Fisher's exact test and/or Z-test are visualized using a heatmap or DAG Logo.
#' @details 
#'   DAG: Differential Amino acid Group
#'   
#'   There are several differences between dagLogo from iceLogo:
#'   
#'   1. The sequence patterns can be grouped by charge, chemistry, hydrophobicity and etc.
#'   
#'   2. dagLogo accepts different length of aligned amino acid sequences.
#'   
#'   3. Except Random, regional (called restricted in dagLogo) and terminal 
#'   (called anchored) background model, the background sequence could be set to 
#'   other regions of the genes in inputs and complementary set of the proteome.
#' @author Jianhong Ou, Haibo Liu, Julie Lihua Zhu
#' 
#' Maintainer: Jianhong Ou <jianhong.ou@duke.edu>
#' 
#' @keywords internal
#' @examples
#'   data("seq.example")
#'   data("proteome.example")
#'   bg <- buildBackgroundModel(seq.example, proteome=proteome.example, numSubsamples=10L)
#'   t <- testDAU(seq.example, bg)
#'   dagLogo(t)
"_PACKAGE"
