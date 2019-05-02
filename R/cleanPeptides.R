#' clean up peptides
#' @description clean up the input peptide subsequences. The function removes 
#' peptides which do NOT contain any anchoring amino acid. Adds peptide for
#' each additional anchor in each peptide, and allows multiple anchoring amino acids.
#' 
#' @param dat input data. The input dat contains two columns `symbol`, protein ID,
#' and `peptides`, peptide sequence.The anchoring amino acid must be in lower case.
#' @param anchors A vector of character, anchoring amino acid must be in lower case.
#' @return A data.frame with columns: `symbol`, `peptides` and `anchor`
#' @export
#' @author Jianhong Ou, Julie Zhu
#' @examples 
#' dat <- read.csv(system.file("extdata", "peptides2filter.csv", package="dagLogo"))
#' dat
#' dat.new <- cleanPeptides(dat, anchors = c("s", "t"))
#' @keywords misc

cleanPeptides <- function(dat, anchors){
  stopifnot(all(c("symbol", "peptides") %in% colnames(dat)))
  stopifnot(is.character(anchors) && length(anchors) > 0)
  
  ## specify 20 amino acids in one-letter symbol
  stopifnot(all(grepl("[galmfwkqespvicyhrndt]", anchors)))
  
  if(!is.data.frame(dat)){
    dat <- as.data.frame(dat, stringsAsFactors=FALSE)
  }
  ## find all the anchors in lower case
  dat <- dat[grepl(paste0("[", paste(anchors, collapse = ""), "]"), 
                   as.character(dat$peptides)), ]
  ## find the positions
  anchorPos <- gregexpr(paste0("[", paste(anchors, collapse = ""), "]"), 
                        as.character(dat$peptides))
  ## keep one lower case per row and change all the other lower case
  ## into upper case.
  len <- lengths(anchorPos)
  dat <- dat[rep(seq.int(nrow(dat)), len), ]
  anchorPos <- unlist(anchorPos)
  dat$anchor <- substr(as.character(dat$peptides), anchorPos, anchorPos)
  dat$peptides <- toupper(dat$peptides)
  substr(dat$peptides, anchorPos, anchorPos) <- dat$anchor
  dat
}
