#' Create an object of \code{\link{Proteome}} Class.
#'
#' Create an object of \code{\link{Proteome}} Class by downloading a whole 
#' proteome data from UniProt for a given organism of an NCBI taxonomy ID or 
#' species' scientific name, or by using peptide sequences in a fasta file.
#'
#' @param source A character vector of length 1 or NULL. A database source from 
#' which the proteome sequences are to be downloaded. By default, currently it is
#' "UniProt". If it is NULL, then \code{fastaFile} has to be specified. The 
#' priority of \code{source} is higher than \code{fastaFile}.
#' @param taxonID An integer(1), specifying Taxonomy ID for a species of interest.
#' Check the NCBI taxonomy database: \url{https://www.ncbi.nlm.nih.gov/taxonomy} 
#' or the UniProt database \url{http://www.uniprot.org/taxonomy/}. At least one 
#' of the two parameters, \code{taxonID} and \code{species}, should be specified. 
#' If both are specified, \code{taxonID} will be used preferentially.
#' @param species A character vector of length 1. The Latin name of a species 
#' confirming to the Linnaean taxonomy nomenclature system. CAUTION: for species 
#' with different strains, attention should be paid. You can interactively choose
#' the right \code{taxonID} from an output list.
#' @param destDir A character vector of length 1. A destination directory with 
#' writing permission for saving downloaded sequences. Default is a temporary 
#' directory in the system's temporary directory.
#' @param fastaFile A character vector of length 1. A fasta file name from which
#' protein sequences are read in.
#' @param ... other parameters passing to the function \code{\link{download.file}}.
#'
#' @importFrom UniProt.ws availableUniprotSpecies lookupUniprotSpeciesFromTaxId
#' @import methods
#' @importFrom Biostrings readAAStringSet
#' @importFrom utils download.file read.table
#' @importFrom httr GET content status_code
#' @return An object of Proteome
#' @export
#' @author Haibo Liu
#' @examples 
#' \dontrun{
#' ## Prepare an objecto of Proteome Class for a proteome from the UniProt database
#' #' proteome <- prepareProteomeByFTP(source = "UniProt", species = "Homo sapiens")
#' }
#' ## Prepare an objecto of Proteome Class from a fasta file
#' fasta <- system.file("extdata", "HUMAN.fasta", package="dagLogo")
#' proteome <- prepareProteomeByFTP(source = NULL, species = "Homo sapiens", 
#' fastaFile=fasta)

prepareProteomeByFTP <- function(source = "UniProt", taxonID = NULL, species = NULL, 
                                 destDir=tempdir(check = TRUE), fastaFile, ...) 
{
  if (is.null(taxonID) && is.null(species))
  {
    stop("Parameter species or taxonID must be specified!", call. = FALSE)
  }
  if (!is.null(source))
  {
    if (!tolower(source) %in% "uniprot")
    {
      stop("Source must be UniProt.", call. = FALSE)
    } else
    {
      source <- "UniProt"
    }
    tempFile <- tempfile(pattern = gsub(" ", ".", paste0(source, species, sep =".")), 
                         tmpdir = destDir, fileext = ".fasta")
    if (!is.null(species) && is.null(taxonID))
    {
      speciesInfo <- availableUniprotSpecies(pattern = paste0("^", species, "$"))
      numberOfEntries <- nrow(speciesInfo)
      if ( numberOfEntries != 1) 
      {
        if(numberOfEntries == 0)
        {
          stop("No matching species in the UniProt database. ",
               "Parameter species might be wrong. ", 
               "Please specify a correct Latin name for the species!",
               call. = FALSE)
        } else 
        {
          if (interactive())
          {
            print(speciesInfo) 
            message("There is multiple matching entries.")
            taxonID <- readline(prompt=paste0("Please select the ",
                                              "right taxon ID based on the printout above: taxon ID: "))
          } else 
          {
            stop("There are multiple entries in the Uniprot database", 
                 " matching your input species. ", 
                 "Please specify the exact one by referring to the ", 
                 "NCBI taxonomy database: ",
                 "https://www.ncbi.nlm.nih.gov/taxonomy ", 
                 "or the UniProt database: http://www.uniprot.org/taxonomy.")
          }
        }
      } else 
      {
        taxonID <- speciesInfo[1, 1]
      }
    } 
    ## check validity of taxon ID 
    species <- lookupUniprotSpeciesFromTaxId(taxonID)
    message(paste0("Downloading data for species: ", species))
    
    proteome_id <- GET(paste0("https://rest.uniprot.org/uniprotkb/search?",
                              "fields=xref_proteomes&format=tsv&query=organism_id:",
                              taxonID))
    if(status_code(proteome_id)!=200){
      stop("Cannot get the organism id for given taxonID: ", taxonID)
    }
    proteome_id <- content(proteome_id, encoding = 'UTF-8')
    proteome_id <- read.table(text = proteome_id,
                              header = FALSE,
                              sep = ":",
                              skip = 1)
    if(nrow(proteome_id)==0){
      stop("Cannot get the organism id for given taxonID: ", taxonID)
    }
    proteome_id <- unique(proteome_id[, 1])
    if (length(proteome_id) == 1) {
      url <- paste0("https://rest.uniprot.org/uniprotkb/stream?compressed=false",
                    "&format=fasta&query=%28proteome%3A", proteome_id, "%29")
      download.file(url = url, destfile = tempFile, ...)
      
      fasta <- readAAStringSet(tempFile)
      isFromUniProt <- TRUE
    } else {
      stop("Can't find the proteome ID for the provided taxonID: ", 
           taxonID)
    }
  } else if (!missing(fastaFile)) ## prepare Proteome from a fasta file
  {
    if (length(fastaFile) == 1 && class(fastaFile) == "character") 
    {
      fasta <- readAAStringSet(fastaFile)
      isFromUniProt <- FALSE
    }
    if (class(fasta) != "AAStringSet") 
    {
      stop("fasta should be an object of AAStringSet",
           call. = FALSE)
    }
    
  } else
  {
    stop("At least one of arguments, source or fastaFile, should be provided!",
         call. = FALSE)
  }
  
  ## ID may be simplified by using gsub
  if (all(grepl("^.+?\\|(.+?)\\|.+", names(fasta))))
  {
    fasta_names <- gsub("^.+?\\|(.+?)\\|.+", "\\1", names(fasta), perl = TRUE)
  }else {
    fasta_names <- gsub("^([^ ]+).+", "\\1", names(fasta), perl = TRUE)
  }
  
  proteome <- data.frame(
    SEQUENCE = unname(as.character(fasta)),
    ID = fasta_names,
    stringsAsFactors = FALSE)
  
  new("Proteome",
      proteome = proteome,
      type = ifelse(isFromUniProt, "UniProt", "fasta"),
      species = species)
}
