#' prepare proteome for background building
#' @description prepare proteome from UniProt webserver or a fasta file
#' @param UniProt.ws an object of UniProt.ws
#' @param fasta fasta file name or an object of AAStringSet
#' @param species an character to assign the species of the proteome
#' @import Biostrings
#' @export
#' @return an object of Proteome which contain protein sequence information.
#' @author Jianhong Ou
#' @seealso \code{\link{formatSequence}}, \code{\link{buildBackgroundModel}}
#' @examples
#' if(interactive()){
#'    library(UniProt.ws)
#'    availableUniprotSpecies("Drosophila melanogaster")
#'    UniProt.ws <- UniProt.ws(taxId=7227)
#'    proteome <- prepareProteome(UniProt.ws, species="Drosophila melanogaster")
#'  }
#' @keywords misc

prepareProteome <- function(UniProt.ws, fasta, species="unknown"){
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
                stop("fasta should be character or an object AAStringSet", call.=FALSE)
            }
            proteome <- data.frame(SEQUENCE=as.character(fasta),
                                   ID=names(fasta),
                                   stringsAsFactors=FALSE)
        }
        return(new("Proteome", 
                   proteome=proteome,
                   type="fasta",
                   species=species))
    }
    stop("Please check you inputs.", call.=FALSE)
}
