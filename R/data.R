#' @docType data
#' @title An object of \code{\link{Proteome-class}} representing the 
#' \emph{Escherichia coli} proteome.
#' @description A dataset containing the \emph{E. coli} proteome.
#' @format
#' An object of  \code{\link{Proteome-class}} for Escherichia coli proteome. 
#' The format is: A list with one data frame and an character.
#' 
#' *`proteome`: 'data.frame':    13780 obs. of  4 variables
#' *`type`: 'character':   "UniProt"
#' *`species`: 'character': "Escherichia coli"
#' 
#' The format of proteome is
#' *`ENTREZ_GENE`: a character vector, records entrez gene id
#' *`SEQUENCE`: a character vector, peptide sequences
#' *`ID`: a character vector, Uniprot ID
#' *`LEN`: a character vector, length of peptides
#' 
#' @source \url{http://www.uniprot.org/}
#' 
#' @details
#' used as an example dataset
#'
#' Annotation data obtained by:
#' 
#'   library(UniProt.ws)
#'   
#'   taxId(UniProt.ws) <- 562
#'   
#'   proteome <- prepareProteome(UniProt.ws, species="Escherichia coli")
#'   
#' @examples
#'  data(ecoli.proteome)
#'  head(ecoli.proteome@proteome)
#'  ecoli.proteome@type
#' @keywords datasets
"ecoli.proteome"

#' @docType data
#' @title An object of \code{\link{Proteome-class}} representing the subset of 
#' \emph{Drosophila melanogaster} proteome.
#' @description The subset \code{\link{Proteome-class}} of fruit fly.
#' @format An object of \code{\link{Proteome-class}} for fly subset proteome. 
#' The format is: A list with one data frame and an character.
#' 
#' *`proteome`: 'data.frame':    1406 obs. of  4 variables
#' *`type`: 'character':   "UniProt"
#' *`species`: 'character': "Drosophila melanogaster"
#' 
#'  The format of proteome is
#'  
#' *`ENTREZ_GENE`: a character vector, records entrez gene id
#' *`SEQUENCE`: a character vector, peptide sequences
#' *`ID`: a character vector, Uniprot ID
#' *`LEN`: a character vector, length of peptides
#' 
#' @source \url{http://www.uniprot.org/}
#' 
#' @details
#' used as an example dataset
#' 
#' Annotation data obtained by:
#' 
#'   library(UniProt.ws)
#'   
#'   taxId(UniProt.ws) <- 7227
#'   
#'   proteome <- prepareProteome(UniProt.ws)
#'   
#'   proteome@proteome <- proteome@proteome[sample(1:19902, 1406), ]
#'   
#' @examples
#' data(proteome.example)
#' head(proteome.example@proteome)
#' proteome.example@type
#' @keywords datasets
"proteome.example"


#' @docType data
#' @title An object of \code{\link{dagPeptides-class}} representing acetylated lysine-containing peptides.
#' @description A dataset containing the acetylated lysine-containing peptides from \emph{Drosophila melanogaster}.
#' @format
#'   An object of \code{\link{dagPeptides-class}} Class
#'   The format is: A list.
#'   
#'     *`data`: 'data.frame':    732 obs. of  7 variables
#'     *`peptides`: 'matrix':   amnio acid in each position
#'     *`upstreamOffset`: an integer, upstream offset position
#'     *`downstreamOffset`: an integer, downstream offset position
#'     *`type`: "character", type of identifiers
#'   
#'   The format of data is
#'
#'     *`IDs`: a character vector, input identifiers
#'     *`anchorAA`: a character vector, anchor amino acid provided in inputs
#'     *`anchorPos`: a numeric vector, anchor position in the protein
#'     *`peptide`: a character vector, peptide sequences
#'     *`anchor`: a character vector, anchor amino acid in the protein
#'     *`upstream`: a character vector, upstream peptides
#'     *`downstream`: a character vector, downstream peptides
#'
#' @details
#'   used as an example dataset
#'   
#'   seq obtained by:
#'   
#'     mart <- useMart("ensembl", "dmelanogaster_gene_ensembl")
#'     
#'     dat <- read.csv(system.file("extdata", "dagLogoTestData.csv", package="dagLogo"))
#'     
#'     seq <- fetchSequence(as.character(dat$entrez_geneid), 
#'     
#'                        anchorPos=as.character(dat$NCBI_site), 
#'                        
#'                        mart=mart, 
#'                        
#'                        upstreamOffset=7, 
#'                        
#'                        downstreamOffset=7)
#' @examples
#'   data(seq.example)
#'   head(seq.example@peptides)
#'   seq.example@upstreamOffset
#'   seq.example@downstreamOffset
#' @keywords datasets
#' 
"seq.example"
