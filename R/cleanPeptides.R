#' clean up peptides
#' @description clean the input file. The function removes peptides without any anchors, 
#' adds peptide for each additional anchor in each peptide, and allows multiple anchors as input.
#' 
#' @param dat input data. The input dat contains two columns symbol and peptides.
#' the anchor AA must be in lower case.
#' @param anchors a vector of character, anchor Amino Acid, must be lower case.
#' @return an data.frame with columns: `symbol`, `peptides` and `anchor`
#' @export
#' @author Jianhong Ou, Julie Zhu
#' @examples 
#' dat <- read.csv(system.file("extdata", "peptides2filter.csv", package="dagLogo"))
#' dat
#' dat.new <- cleanPeptides(dat, anchors = c("s", "t"))
#' @keywords misc

cleanPeptides <- function(dat, anchors){
  stopifnot(all(c("symbol", "peptides") %in% colnames(dat)))
  stopifnot(is.character(anchors))
  stopifnot(all(!grepl("[^a-z]", anchors)))
  stopifnot(length(anchors)>0)
  if(!is.data.frame(dat)){
    dat <- as.data.frame(dat, stringsAsFactors=FALSE)
  }
  ## find all the anchors in lower case
  dat <- dat[grepl(paste0("[", paste(anchors, collapse = ""), "]"), 
                   as.character(dat$peptides)), ]
  ## find the positions
  anchorPos <- gregexpr(paste0("[", paste(anchors, collapse = ""), "]"), 
                        as.character(dat$peptides))
  ## keep one lower case per row and change all the other lower case into capital.
  len <- lengths(anchorPos)
  dat <- dat[rep(seq.int(nrow(dat)), len), ]
  anchorPos <- unlist(anchorPos)
  dat$anchor <- substr(as.character(dat$peptides), anchorPos, anchorPos)
  dat$peptides <- toupper(dat$peptides)
  substr(dat$peptides, anchorPos, anchorPos) <- dat$anchor
  dat
}
