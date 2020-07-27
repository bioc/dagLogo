#' prepare proteome for background building
#' @description prepare proteome from UniProt webserver or a fasta file
#' @param source An object of \code{\link{UniProt.ws}} or A character "UniProt".
#' @param fasta fasta file name or an object of AAStringSet
#' @param species an character to assign the species of the proteome
#' @param ... parameters could be passed to \link{prepareProteomeByFTP}.
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

prepareProteome <- function(source, fasta, species="unknown", ...){
  if(!missing(source)){
    if(is(source, "UniProt.ws")){
      return(prepareProteomeByUniProtWS(UniProt.ws=source, species=species))
    }else{
      res <- lapply(species, function(.ele){
        prepareProteomeByFTP(source = "UniProt", species = .ele, ...)
      })
      prot <- res[[1]]
      if(length(res)>1){
        prot@proteome <- 
          do.call(rbind, lapply(res, function(.ele) .ele@proteome))
        prot@species <- species
      }
      return(prot)
    }
  }else{
    if(!missing(fasta)){
      return(prepareProteomeByUniProtWS(fasta=fasta, species=species))
    }
  }
  stop("Please check you inputs.", call.=FALSE)
}
