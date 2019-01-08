#' @docType class
#' @title Class \code{"dagPeptides"}
#' @description An object of class \code{"dagPeptides"} represents the information of peptides.
#' @section Objects from the Class: 
#'  Objects can be created by calls of the form 
#'  \code{new("dagPeptides", data, peptides, upstreamOffset, downstreamOffset, type)}.
#' @slot data Object of class \code{"data.frame"} The details of the input sequences. 
#' It includes the columns: IDs, anchorAA (anchor Amino Acid), anchorPos (anchor Position),
#' peptide (protein peptide), anchor, upstream, downstream (peptides in given upstream and 
#' downstream offset from anchor)
#' @slot peptides code{"matrix"} The input peptides. Each column contains one peptide in that position
#' @slot upstreamOffset \code{"numeric"} The upstream offset from anchor
#' @slot downstreamOffset \code{"numeric"} The downstream offset from anchor
#' @slot type \code{"charactger"} ID type of inputs
#' @keywords classes
#' @export
setClass("dagPeptides",
         representation(data="data.frame",
                        peptides="matrix",
                        upstreamOffset="numeric",
                        downstreamOffset="numeric",
                        type="character"),
         validity=function(object){
             re<-TRUE
             if(object@upstreamOffset < 0 || object@downstreamOffset < 0) 
                 re <- "upstreamOffset and downstreamOffset should be a integer greater than 0"
             peptides <- as.character(object@peptides)
             peptides <- peptides[(!is.na(peptides)) & (peptides!="NA")]
             if(!all(1==nchar(peptides)))
                 re <- "peptides must be a matrix with one amino acide in each position"
             re
         })



#' @docType class
#' @title Class \code{"Proteome"}
#' @description An object of class \code{"Proteome"} represents proteome of a given species.
#' @section Objects from the Class:
#'  Objects can be created by calls of the form 
#'  \code{new("Proteome", proteome, type, species)}.
#' @slot proteome Object of class \code{"data.frame"} the proteome of a given species, 
#' should include ids and peptide sequences.
#' @slot type \code{"character"} indicates how the object is prepared, 
#' could be "fasta" or "UniProt"
#' @slot species \code{"character"} the species
#' @keywords classes
#' @export

setClass("Proteome",
         representation(proteome="data.frame", type="character", species="character"),
         validity=function(object){
             re <- TRUE
             if(!object@type %in% c("fasta", "UniProt"))
                 re <- "type must be fasta or UniProt"
             if(object@type=="UniProt" && is.null(object@proteome$ENTREZ_GENE))
                 re <- "when type equals to UniProt, ENTREZ_GENE column is required for proteome"
             if(is.null(object@proteome$SEQUENCE))
                 re <- "proteome sequence is required"
             re
         })

#' @docType class
#' @title Class \code{"dagBackground"}
#' @description An object of class \code{"dagBackground"} represents background model.
#' @section Objects from the Class:
#'   Objects can be created by calls of the form 
#'   \code{new("dagBackground", background, permutationSize)}.
#' @slot background Object of class \code{"list"} records the background model
#' @slot permutationSize code{"integer"} permutation size of background
#' @keywords classes
#' @export

setClass("dagBackground",
         representation(background="list", 
                        permutationSize="integer"))

#' @docType class
#' @title Class \code{"testDAUresults"}
#' @description
#' An object of class \code{"testDAUresults"} represents background model.
#' @section Objects from the Class:
#'   Objects can be created by calls of the form 
#'  \code{new("dagBackground", group="character",
#'            difference="matrix",
#'            zscore="matrix",
#'            pvalue="matrix",
#'            background="matrix",
#'            motif="matrix",
#'            upstream="numeric",
#'            downstream="numeric")}.
#' @slot group Object of class \code{"character"} could be "null", "classic", 
#' "charge", "chemistry", "hydrophobicity"
#' @slot difference code{"matrix"} the difference of inputs from background for 
#' each amino acid in each position
#' @slot zscore code{"matrix"} z score for each amino acid in each position
#' @slot pvalue code{"matrix"} pvalue for each amino acid in each position
#' @slot background code{"matrix"} background frequencies for each amino acid in each position
#' @slot motif code{"matrix"} inputs frequencies for each amino acid in each position
#' @slot upstream \code{"numeric"} The upstream offset from anchor
#' @slot downstream \code{"numeric"} The downstream offset from anchor
#' @keywords classes
#' @export

setClass("testDAUresults",
         representation(group="character",
                        difference="matrix",
                        zscore="matrix",
                        pvalue="matrix",
                        background="matrix",
                        motif="matrix",
                        upstream="numeric",
                        downstream="numeric"),
         validity=function(object){
             re<-TRUE
             if(object@upstream < 0 || object@downstream < 0) 
                 re <- "upstream and downstream should be a integer greater than 0"
             
             if(ncol(object@zscore)==0 || ncol(object@difference)==0 || ncol(object@pvalue)==0)
                 re <- "slots zscore, difference and pvalue could not be empty"
             if(any(dim(object@zscore)!=dim(object@difference)) || any(dim(object@zscore)!=dim(object@pvalue)))
                 re <- "dim of slots zscore, difference and pvalue should be identical"
             re
         })
