#' Prepare a Proteome object for background building
#' 
#' @description Create an object of \code{\link{Proteome}} Class by query
#' the UniProt database of an organism of a given species' scientific name,
#' or by using peptide sequences in a fasta file or in an AAStringSet object.
#' 
#' @param UniProt.ws An object of \code{\link{UniProt.ws}}.
#' @param fasta A fasta file name or an object of \code{\link{AAStringSet}}.
#' @param species An character vector of length (1) to designate the species
#' of the proteome
#' @importFrom Biostrings readAAStringSet
#' @export
#' @return An object of Proteome which contain protein sequence information.
#' @author Jianhong Ou
#' @seealso \code{\link{formatSequence}}, \code{\link{buildBackgroundModel}}
#' @examples
#' if(interactive()){
#'    library(UniProt.ws)
#'    availableUniprotSpecies("Drosophila melanogaster")
#'    UniProt.ws <- UniProt.ws(taxId=7227)
#'    proteome <- prepareProteomeByUniProtWS(UniProt.ws, species="Drosophila melanogaster")
#'  }
#' @keywords misc

prepareProteomeByUniProtWS <- function(UniProt.ws, fasta, species="unknown"){
    if(!missing(UniProt.ws) && class(UniProt.ws)=="UniProt.ws"){
        egs <- keys(UniProt.ws, "ENTREZ_GENE")
        cols <- c("SEQUENCE", "ID")
        proteome <- select(UniProt.ws, egs, cols, "ENTREZ_GENE")
        proteome$SEQUENCE <- gsub(" ", "", proteome$SEQUENCE, fixed=TRUE)
        proteome$LEN <- nchar(proteome$SEQUENCE)
        return(new("Proteome",
                   proteome=proteome,
                   type="UniProt",
                   species=species))
    } else {
        if(!missing(fasta)){
            if(length(fasta)==1 && class(fasta)=="character"){
                fasta <- readAAStringSet(fasta)
            }
            if(class(fasta)!="AAStringSet"){
                stop("fasta should be character or an object of AAStringSet",
                     call.=FALSE)
            }
            proteome <- data.frame(SEQUENCE=as.character(fasta),
                                   ID=names(fasta),
                                   stringsAsFactors=FALSE)
            return(new("Proteome", 
                       proteome=proteome,
                       type="fasta",
                       species=species))
        }
    }
    stop("Please check you inputs.", call.=FALSE)
}
