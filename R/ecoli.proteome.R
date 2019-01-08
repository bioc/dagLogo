#' @docType data
#' @title the subset proteome of Escherichia coli
#' @description the subset proteome of Escherichia coli
#' @format
#' An object of Proteome for Escherichia coli proteome. 
#' The format is: A list with one data frame and an character.
#' 
#' *`proteome`: 'data.frame':     obs. of  4 variables
#' *`type`: 'character':   "UniProt"
#' 
#' The format of proteome is
#' *`ENTREZ_GENE`: a character vector, records entrez gene id
#' *`SEQUENCE`: a character vector, peptide sequences
#' *`ID`: a character vector, Uniprot ID
#' *`LEN`: a character vector, length of peptides
#' 
#' @details
#'used in the examples
#' Annotation data obtained by:
#'   library(UniProt.ws)
#'   taxId(UniProt.ws) <- 562
#'   proteome <- prepareProteome(UniProt.ws, species="Escherichia coli")
#' @examples
#'  data(ecoli.proteome)
#'  head(ecoli.proteome@proteome)
#'  ecoli.proteome@type
#' @keywords datasets
"ecoli.proteome"
