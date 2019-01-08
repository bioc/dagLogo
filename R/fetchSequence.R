#' @title fetch sequence by id
#' @description  fetch amino acid sequence by given identifiers via biomaRt or proteome prepared by \code{\link{prepareProteome}}
#' @param IDs A vector of Identifiers to retrieve peptides
#' @param type type of identifiers
#' @param anchorAA a vector of character, anchor Amino Acid
#' @param anchorPos a vector of character or numeric, anchor position, for example, K121. 
#' Or a vector of character with amino acid sequences. 
#' If AA sequences is used, the anchorAA must be the a vector of character with single AA for each.
#' @param mart an object of Mart
#' @param proteome an object of Proteome, output of \code{\link{prepareProteome}}
#' @param upstreamOffset an integer, upstream offset position
#' @param downstreamOffset an integer, downstream offset position
#' @import biomaRt
#' @export
#' @return an object of \code{\linkS4class{dagPeptides}}
#' @author Jianhong Ou, Alexey Stukalov, Julie Zhu
#' @seealso \code{\link{formatSequence}}
#' @examples
#' if(interactive()){
#'   mart <- useMart("ensembl", "dmelanogaster_gene_ensembl")
#'   dat <- read.csv(system.file("extdata", "dagLogoTestData.csv", package="dagLogo"))
#'   seq <- fetchSequence(as.character(dat$entrez_geneid[1:5]), 
#'                        anchorPos=as.character(dat$NCBI_site[1:5]), 
#'                        mart=mart, 
#'                        upstreamOffset=7, 
#'                        downstreamOffset=7)
#'   ## sample: use sequence as anchorPos
#'   sequences <- seq@peptides
#'   sequences[, 8] <- "k"
#'   sequences <- apply(sequences, 1, paste, collapse="")
#'   seq <- fetchSequence(as.character(seq@data$IDs),
#'                        anchorAA="k",
#'                        anchorPos=sequences,
#'                        mart=mart,
#'                        upstreamOffset=7,
#'                        downstreamOffset=7)
#'   ## sample: use sequence as anchorPos 2
#'   sequences <- cbind(seq@peptides[, 1:8], "*", seq@peptides[, 9:15])
#'   sequences <- apply(sequences, 1, paste, collapse="")
#'   seq <- fetchSequence(as.character(seq@data$IDs),
#'                        anchorAA="*",
#'                        anchorPos=sequences,
#'                        mart=mart,
#'                        upstreamOffset=7,
#'                        downstreamOffset=7)
#' }
#' @keywords misc

fetchSequence <- function(IDs, type="entrezgene", anchorAA=NULL, anchorPos,
                            mart, proteome, upstreamOffset, downstreamOffset){
    if(missing(mart) && missing(proteome)){
        stop("Need mart or proteome for fetching sequences.", call.=FALSE)
    }
    if(!missing(mart) && class(mart)!="Mart"){
        stop("mart should be an object of Mart", call.=FALSE)
    }
    if(!missing(proteome) && class(proteome)!="Proteome"){
        stop("proteome should be an object of Proteome. \n
                 Try ?prepareProteome to get help", call.=FALSE)
    }
    if(missing(upstreamOffset) || missing(downstreamOffset)){
        stop("Please indicate the upstreamOffset and downstreamOffset position offset anchor amino acid", call.=FALSE)
    }
    if(upstreamOffset < 0 || downstreamOffset < 0){
        stop("upstreamOffset and downstreamOffset should be a integer greater than 0", call.=FALSE)
    }
    if(downstreamOffset > 20 || upstreamOffset > 20){
        stop("upstreamOffset and downstreamOffset should be the offset of anchor amino acid, no greater than 20", call.=FALSE)
    }
    if(missing(IDs) || missing(anchorPos)){
        stop("Missing required arguments IDs or anchorPos", call.=FALSE)
    }
    if(class(IDs)!="character" || 
           !inherits(anchorPos, c("character", "numeric", "integer"))){
        stop("IDs must be characters and anchorPos should be character or number",
             call.=FALSE)
    }
    if(any(is.na(IDs)) || any(is.na(anchorPos))){
        stop("IDs or anchorPos contains NA", call.=FALSE)
    }
    if(length(IDs)!=length(anchorPos)){
        stop("length of IDs and anchorPos are not identical.", call.=FALSE)
    }
    if(length(anchorAA)>0){
        if(any(nchar(anchorAA)!=1)){
            stop("anchorAA must be a single amino acid or *", call.=FALSE)
        }
    }
    searchAnchor <- FALSE
    if(class(anchorPos)=="character"){
        anchorPos <- gsub("^\\-+", "", anchorPos)
        anchorPos <- gsub("\\-+$", "", anchorPos)
        if(any(!grepl("^[A-Z]\\d+$", toupper(anchorPos)))){
            if(length(anchorAA)<1 || any(grepl("[^*A-Z]", toupper(anchorPos)))){
                stop("anchorPos should be the amino acid followed by the position, eg. K123", 
                     "Otherwise, anchorPos should be the strings of amino acid and anchorAA is the anchor amino acid", call.=FALSE)
            }
            searchAnchor <- TRUE
            anchorPos <- strsplit(anchorPos, "")
            anchor <- mapply(function(x, y){
                which(x==y)
            }, anchorPos, anchorAA, SIMPLIFY = FALSE)
            if(any(anchorAA=="*")){
                ## asterisk
                if(length(anchorAA)==length(anchorPos)){
                    anchor[anchorAA=="*"] <- 
                        lapply(anchor[anchorAA=="*"], function(.ele) 
                            .ele - 1:length(.ele))
                }else{
                    if(length(anchorAA)==1){
                        anchor <- lapply(anchor, function(.ele) 
                            .ele - 1:length(.ele))
                    }else{
                        stop("length of anchorAA must be 1 or equal to anchorPos")
                    }
                }
                anchorPos <- lapply(anchorPos, function(.ele) .ele[.ele!="*"])
            }
            anchorPos <- sapply(anchorPos, paste, collapse="")
            anchorPos <- toupper(anchorPos)
        }else{
            anchorPos <- toupper(anchorPos)
            anchorAA <- substr(anchorPos, 1, 1)
            anchorPos <- as.numeric(substring(anchorPos, 2))
        }
    }
    inputs <- data.frame(IDs, anchorAA, anchorPos, oid=1:length(anchorPos))
    ## retreive sequence
    if(!missing(mart)){
        protein <- getSequence(id=unique(as.character(IDs)), 
                               type=type,
                               seqType="peptide",
                               mart=mart)
    }else{
        if(!(type %in% c("entrezgene", "UniProtKB_ID"))){
            stop("Only accept 'entrezgene' or 'UniProtKB_ID' for type when using proteome.",
                 call.=FALSE)
        }
        if(type=="entrezgene"){
            protein <- proteome@proteome[proteome@proteome[,"ENTREZ_GENE"] %in% unique(as.character(IDs)), c(2, 1)]
            colnames(protein) <- c("peptide", "entrezgene")
        }else{
            protein <- proteome@proteome[proteome@proteome[,"ID"] %in% unique(as.character(IDs)), c("SEQUENCE", "ID")]
            colnames(protein) <- c("peptide", "UniProtKB_ID")
        }
        protein <- protein[!is.na(protein[, 2]),]
    }
    
    dat <- merge(inputs, protein, by.x=1, by.y=2)
    dat$peptide <- toupper(dat$peptide)
    
    if(searchAnchor){
        anchorPos <- mapply(function(.ele, .pep){
            start(matchPattern(AAString(toupper(.ele)), .pep))
        }, dat$anchorPos, dat$peptide, SIMPLIFY = FALSE)
        anchorPos <- mapply(function(.pos, .anchor){
            if(length(.pos)==0){
                return(integer())
            }
            .pos[1] + .anchor - 1
        }, anchorPos, anchor[dat$oid], SIMPLIFY=FALSE)
        dat <- dat[rep(1:nrow(dat), sapply(anchorPos, length)), ]
        dat$anchorPos <- unlist(anchorPos)
        dat$anchorAA <- 
            unlist(mapply(function(pep, pos){substr(pep, pos, pos)},
                          dat$peptide, dat$anchorPos))
    }
    
    dat$anchor <- unlist(mapply(function(pep, pos){substr(pep, pos, pos)},
                              dat$peptide, dat$anchorPos))
    ##colnames(dat)==c("IDs", "anchorAA", "anchorPos", "peptide", "anchor")
    ## check sequence of NCBIsites
    if(!is.null(anchorAA)[1]){
        dat <- dat[toupper(dat$anchorAA)==dat$anchor, ]
    }
    ## extract sequences for logo
    upstreamGuard <- paste0( rep.int( "?", upstreamOffset ), collapse = '' )
    downstreamGuard <- paste0( rep.int( "?", downstreamOffset ), collapse = '' )
    peptide.guarded <- paste0( upstreamGuard, dat$peptide, downstreamGuard )
    dat$upstream <- substr(peptide.guarded, dat$anchorPos, dat$anchorPos+upstreamOffset-1)
    dat$downstream <- substr(peptide.guarded, dat$anchorPos+upstreamOffset+1, dat$anchorPos+upstreamOffset+downstreamOffset)
    # unique dat by oid and upstream/anchor/downstream sequence
    dat <- dat[!duplicated(paste(dat$oid, dat$upstream, dat$downstream)), ]
    dat$oid <- NULL
    rownames(dat) <- NULL
    # convert logo sequences into character matrix
    seqchar.upstream <- do.call(rbind, strsplit(dat$upstream, "", fixed=TRUE))
    seqchar.upstream[ seqchar.upstream == '?' ] <- NA
    seqchar.downstream <- do.call(rbind, strsplit(dat$downstream,"", fixed=TRUE))
    seqchar.downstream[ seqchar.downstream == '?' ] <- NA
    
    seqchar <- cbind(seqchar.upstream, dat$anchor, seqchar.downstream)
    new("dagPeptides", data=dat, peptides=seqchar, 
                   upstreamOffset=upstreamOffset, 
                   downstreamOffset=downstreamOffset, 
                   type=type)
}