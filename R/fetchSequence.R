#' Fetch protein/peptide sequences and create a \code{\link{dagPeptides-class}} object.
#'
#' This function fetches protein/peptide sequences from a Biomart database or 
#' from a \code{\link{Proteome-class}} object based on protein/peptide IDs and create 
#' a \code{\link{dagPeptides-class}} object following restriction as specified by 
#' parameters: anchorAA or anchorPos, upstreamOffset and downstreamOffset.
#'
#' @param IDs  A character vector containing protein/peptide IDs used to fetch 
#' sequences from a Biomart database or a \code{\link{Proteome-class}} object.
#' @param type A character vector of length 1. The available options are 
#' "entrezgene" and "uniprotswissprot" if parameter \code{mart} is missing;
#' otherwise it can be any type of IDs available in Biomart databases.
#' @param anchorAA A character vector of length 1 or the same length as that of
#' anchorPos, each element of which is a single letter symbol of amino acids, 
#' for example, "K" for lysine.
#' @param anchorPos A character or numeric vector. Each element of which is (1) a 
#' single-letter symbol of amino acid followed by the position of the anchoring 
#' amino acid in the target peptide/protein sequence, for example, "K123" for lysine 
#' at position 123 or the position of the anchoring amino acid in the target 
#' peptide/protein sequence, for example, "123" for an amino acid at position 123; 
#' or (2) a vector of subsequences containing the anchoring AAs.
#' @param mart A Biomart database name you want to connect to. Either of parameters 
#' \code{mart} or \code{proteome} should be provided.
#' @param proteome An object of \code{\link{Proteome-class}}. Either of parameters 
#' \code{mart} or \code{\link{Proteome-class}} should be provided.
#' @param upstreamOffset An integer, the upstream offset relative to
#' the anchoring position.
#' @param downstreamOffset An integer, the downstream offset relative
#' to the anchoring position.
#' @import biomaRt
#' @importFrom BiocGenerics start
#' @importFrom utils adist
#' @importFrom Biostrings AAString matchPattern
#' @import methods
#' @return An object of class \code{\link{dagPeptides-class}} 
#' @export
#' 
#' @examples
#' ## Case 1: You have both positions of the anchoring AAs and the identifiers 
#' ## of their enclosing peptide/protein sequences for fetching sequences using 
#' ## the fetchSequence function via the Biomart.
#' 
#' if (interactive())
#' {
#'     try({
#'     mart <- useMart("ensembl")
#'     fly_mart <-
#'        useDataset(mart = mart, dataset = "dmelanogaster_gene_ensembl")
#'     dat <- read.csv(system.file("extdata", "dagLogoTestData.csv",
#'                            package = "dagLogo"))
#'     seq <- fetchSequence(
#'        IDs = as.character(dat$entrez_geneid),
#'        anchorPos = as.character(dat$NCBI_site),
#'        mart = fly_mart,
#'        upstreamOffset = 7,
#'        downstreamOffset = 7)
#'    head(seq@peptides)
#'    })
#' }
#' 
#' 
#' ## Case 2: You don't have the exactly postion information, but You have the 
#' ## interesting peptide subsequences and the identifiers of their enclosing 
#' ## peptide/protein sequences for fetching sequences using the fetchSequence
#' ## function via the Biomart. In the following examples, the anchoring AAs 
#' ## are marked by asterisks. 
#' if (interactive())
#' {
#'     try({
#'         mart <- useMart("ensembl")
#'         fly_mart <-
#'             useDataset(mart = mart, dataset = "dmelanogaster_gene_ensembl")
#'         dat <- read.csv(system.file("extdata", "dagLogoTestData.csv",
#'                                     package = "dagLogo"))
#'         seq <- fetchSequence(
#'             IDs = as.character(dat$entrez_geneid),
#'             anchorAA = "*",
#'             anchorPos = as.character(dat$peptide),
#'             mart = fly_mart,
#'             upstreamOffset = 7,
#'             downstreamOffset = 7
#'         )
#'         head(seq@peptides)
#'     })
#' }
#' 
#' ## In following example, the anchoring AAs are lower-case "s" for amino acid 
#' ## serine.
#' if(interactive())
#' {
#'    try({
#'        dat <- read.csv(system.file("extdata", "peptides4dagLogo.csv",
#'                                    package = "dagLogo"))
#'         mart <- useMart("ensembl")
#'         human_mart <-
#'             useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
#'         seq <- fetchSequence(IDs = toupper(as.character(dat$symbol)),
#'                              type = "hgnc_symbol",
#'                              anchorAA = "s",
#'                              anchorPos = as.character(dat$peptides),
#'                              mart = human_mart,
#'                              upstreamOffset = 7,
#'                              downstreamOffset = 7)
#'         head(seq@peptides)
#'     })
#' }
#' 
#' 

fetchSequence <-function(IDs,
                        type = "entrezgene",
                        anchorAA = NULL,
                        anchorPos,
                        mart,
                        proteome,
                        upstreamOffset,
                        downstreamOffset) 
{
    if (missing(mart) && missing(proteome)) 
    {
        stop("Need mart or proteome for fetching sequences.", call. = FALSE)
    }
    if (!missing(mart) && class(mart) != "Mart") 
    {
        stop("mart should be an object of the Mart class", call. = FALSE)
    }
    if (!missing(proteome) && class(proteome) != "Proteome") 
    {
        stop("proteome should be an object of Proteome class.",
            "Try ?prepareProteome to get help.", call. = FALSE)
    }
    if (missing(upstreamOffset) || missing(downstreamOffset)) 
    {
        stop("Please provide the upstreamOffset and downstreamOffset positions 
            relative to the anchoring amino acid", call. = FALSE)
    }
    if (upstreamOffset < 0 || downstreamOffset < 0) 
    {
        stop("The upstreamOffset and downstreamOffset should be positive integers.",
             call. = FALSE)
    }
    # if (downstreamOffset > 20 || upstreamOffset > 20) 
    # {
    #     stop("The upstreamOffset and downstreamOffset should be a positive integer less than 20",
    #          call. = FALSE)
    # }
    if (missing(IDs) || missing(anchorPos)) 
    {
        stop("Missing required parameter IDs or anchorPos", call. = FALSE)
    }
    if (class(IDs) != "character" ||
        !inherits(anchorPos, c("character", "numeric", "integer"))) 
    {
        stop("IDs must be characters and anchorPos should be a character or number",
             call. = FALSE)
    }
    if (any(is.na(IDs)) || any(is.na(anchorPos))) 
    {
        stop("IDs or anchorPos contains NA", call. = FALSE)
    }
    if (length(IDs) != length(anchorPos)) 
    {
        stop("The length of IDs and anchorPos are not identical.", call. = FALSE)
    }
    if (length(anchorAA) > 0 && any(nchar(anchorAA) != 1)) 
    {
        stop("anchorAA must be a single amino acid or *", call. = FALSE)
    }
   
    searchAnchor <- FALSE
    anchor <- anchorPos
    if (class(anchorPos) == "character") 
    {
        ## removing leading and trailing "-"
        anchorPos <- gsub("^\\-+", "", anchorPos)
        anchorPos <- gsub("\\-+$", "", anchorPos)
        
        ## anchorPos eg. K135
        if (any(!grepl("^[A-Z]\\d+$", toupper(anchorPos)))) 
        {
             if (length(anchorAA) < 1 || any(grepl("[^*A-Z]", toupper(anchorPos)))) 
            {
                stop("anchorPos should be the amino acid followed by the position, eg. K123. ",
                    "Otherwise, anchorPos should be the strings of amino acid and 
                    anchorAA is the anchoring amino acid", call. = FALSE)
            }
            searchAnchor <- TRUE
            anchorPos <- strsplit(anchorPos, "")
            
            ## find relative index of the anchoring amino acids
            anchor <- mapply(function(x, y) {which(x == y)}, anchorPos, 
                             anchorAA, SIMPLIFY = FALSE)
            
            ## amino acid highlighted using "*": the AA immediately before "*"
            ## is the anchoring AA.
            if (any(anchorAA == "*")) 
            {
                if (length(anchorAA) == length(anchorPos)) 
                {
                    ## the index of the amino acid immediately before the "*"
                    ## the actual index should be index -1
                    anchor[anchorAA == "*"] <-
                        lapply(anchor[anchorAA == "*"], function(.ele){.ele - 1:length(.ele)})
                } else
                {
                    if (length(anchorAA) == 1) 
                    {
                        anchor <- sapply(anchor, function(.ele){.ele - 1:length(.ele)})
                    } else
                    {
                        stop("length of anchorAA must be 1 or equal to that of anchorPos.")
                    }
                }
                ## remove "*" 
                anchorPos <- lapply(anchorPos, function(.ele){.ele[.ele != "*"]})
            }
            ## and collapse into a single string and change to uppper cases
            anchorPos <- sapply(anchorPos, paste, collapse = "")
            anchorPos <- toupper(anchorPos)
        } else
        {
            anchorPos <- toupper(anchorPos)
            ## single character symbol of amino acid
            anchorAA <- substr(anchorPos, 1, 1)
            ## now position is only number
            anchorPos <- as.numeric(substring(anchorPos, 2))
            anchor <- anchorPos
        }
    }
    inputs <- data.frame(IDs, anchorAA, anchorPos, oid = seq_along(anchorPos))
    if(is.list(anchor)){
      ids <- lengths(anchor)
      inputs <- inputs[rep(seq_along(anchorPos), ids), , drop=FALSE]
      inputs$anchor <- unlist(anchor)
    }else{
      inputs$anchor <- anchor
    }
    
    if (!missing(mart))  ## retreive sequence from biomart
    {
        possibleTypeIds <- listFilters(mart = mart)
        if(!type[1] %in% possibleTypeIds[, 1]){
          if(type[1]=="entrezgene" && "entrezgene_id" %in% possibleTypeIds[, 1]){
            type <- "entrezgene_id"
          }else{
            ad <- adist(as.character(possibleTypeIds[, 1]), type[1])
            possibleType <- as.character(possibleTypeIds[which(ad==min(ad, na.rm = TRUE)), 1])
            stop("Invalid type argument. Use the listFilters function to select a valid type argument.",
                 paste("Do you mean", paste(possibleType, collapse = ", ")))
          }
        }
        protein <- getSequence(id = unique(as.character(IDs)),
                               type = type,
                               seqType = "peptide",
                               mart = mart)
        
        protein <- protein[!protein$peptide =="Sequence unavailable", ]
        
        if (nrow(protein) < 1)
        {
            stop("Too few retrieved protein sequences from Ensembl Biomart.",
                 "Make sure the IDs and their type are correct!", call. = FALSE)
        }
    } else  ## retreive sequence from Proteome
    {
        
        if (!(type %in% c("entrezgene", "uniprotswissprot"))) 
        {
            stop( "Only accept 'entrezgene' or 'uniprotswissprot' for type when 
                  using proteome for fetching sequences.", call. = FALSE)
        }
        if (type == "entrezgene") 
        {
            protein <-
                proteome@proteome[proteome@proteome[, "ID"] %in% 
                                      unique(as.character(IDs)), c(2, 1)]
            colnames(protein) <- c("peptide", "entrezgene")
        } else
        {
            protein <-
                proteome@proteome[proteome@proteome[, "ID"] %in% 
                                      unique(as.character(IDs)),
                                  c("SEQUENCE", "ID")]
            colnames(protein) <- c("peptide", "uniprotswissprot")
        }
        protein <- protein[!is.na(protein[, 2]), ]
    }
    
    dat <- merge(inputs, protein, by.x = 1, by.y = 2)
    dat$peptide <- toupper(dat$peptide)
    
    if (searchAnchor) 
    {
        ## get the absolute index for the first occurence of query peptides
        anchorPos <- mapply(function(.ele, .pep) {
            BiocGenerics::start(matchPattern(AAString(toupper(.ele)), .pep))
        }, dat$anchorPos, dat$peptide, SIMPLIFY = FALSE)
        ## get the absolute index for the anchoring amino acids in the first 
        ## occurence of query peptide
        anchorPos <- mapply(function(.pos, .anchor) {
            if (length(.pos) == 0) 
            {
                return(integer())
            }
            .pos[1] + .anchor - 1
        }, anchorPos, dat$anchor, SIMPLIFY = FALSE)
        
        ## remove peptide without queryed "anchoring peptide"
        dat <- dat[rep(seq.int(nrow(dat)), lengths(anchorPos)), , drop=FALSE]
        
        ## replacing query peptide with the absolute index of anchoring amino acid
        dat$anchorPos <- unlist(anchorPos)
        dat$anchorAA <- unlist(mapply(function(pep, pos) {
            substr(pep, pos, pos)},dat$peptide, dat$anchorPos))
    } 
    ## replace indice with symbols of anchoring amino acid
    dat$anchor <-unlist(mapply(function(pep, pos) {
                substr(pep, pos, pos)}, dat$peptide, dat$anchorPos))

    ## check sequence of NCBIsites
    if (!is.null(anchorAA)[1])
    {
        dat <- dat[toupper(dat$anchorAA) == dat$anchor, , drop=FALSE]
    }
    
    ## extract sequences for logo
    upstreamGuard <-
        paste0(rep.int("?", upstreamOffset), collapse = '')
    downstreamGuard <-
        paste0(rep.int("?", downstreamOffset), collapse = '')
    peptide.guarded <-
        paste0(upstreamGuard, dat$peptide, downstreamGuard)
    if(length(dat$anchorPos)==0){
      stop("Cannot find any sequence by given anchors.")
    }
    dat$upstream <-
        substr(peptide.guarded,
               dat$anchorPos,
               dat$anchorPos + upstreamOffset - 1)
    dat$downstream <-
        substr(peptide.guarded,
            dat$anchorPos + upstreamOffset + 1,
            dat$anchorPos + upstreamOffset + downstreamOffset)
    
    # unique dat by oid and upstream/anchor/downstream sequence
    dat <- dat[!duplicated(paste(dat$oid, dat$upstream, dat$downstream)), ]
    dat$oid <- NULL
    rownames(dat) <- NULL
    
    # convert logo sequences into character matrix
    seqchar.upstream <-
        do.call(rbind, strsplit(dat$upstream, "", fixed = TRUE))
    seqchar.upstream[seqchar.upstream == '?'] <- NA
    seqchar.downstream <-
        do.call(rbind, strsplit(dat$downstream, "", fixed = TRUE))
    seqchar.downstream[seqchar.downstream == '?'] <- NA
    seqchar <- cbind(seqchar.upstream, dat$anchor, seqchar.downstream)
    
    new("dagPeptides",
        data = dat,
        peptides = seqchar,
        upstreamOffset = upstreamOffset,
        downstreamOffset = downstreamOffset,
        type = type)
}