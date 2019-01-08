#' @docType data
#' @title the subset proteome of fruit fly
#' @description the subset proteome of fruit fly
#' @format An object of Proteome for fly subset proteome. 
#' The format is: A list with one data frame and an character.
#' 
#' *`proteome`: 'data.frame':    1406 obs. of  4 variables
#' *`type`: 'character':   "UniProt"
#' 
#'  The format of proteome is
#'  
#' *`ENTREZ_GENE`: a character vector, records entrez gene id
#' *`SEQUENCE`: a character vector, peptide sequences
#' *`ID`: a character vector, Uniprot ID
#' *`LEN`: a character vector, length of peptides
#' 
#' @details
#' used in the examples
#' Annotation data obtained by:
#'   library(UniProt.ws)
#'   taxId(UniProt.ws) <- 7227
#'   proteome <- prepareProteome(UniProt.ws)
#'   proteome@proteome <- proteome@proteome[sample(1:19902, 1406), ]
#' @examples
#' data(proteome.example)
#' head(proteome.example@proteome)
#' proteome.example@type
#' @keywords datasets
"proteome.example"