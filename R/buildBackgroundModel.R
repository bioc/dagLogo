#' @title  build background model
#' @description build background model for dag test
#' @param dagPeptides an object of dagPeptides, output of \code{\link{fetchSequence}} or 
#' \code{\link{formatSequence}}
#' @param bg could be "wholeGenome", "inputSet" or "nonInputSet"
#' @param model could be "any" or "anchored"
#' @param targetPosition could be "any", "Nterminus" or "Cterminus"
#' @param uniqueSeq should the background sequence be unique?
#' @param permutationSize how many times should it samples
#' @param rand.seed random seed
#' @param replacement Should sampling be with replacement?
#' @param proteome an object of Proteome, output of \code{\link{prepareProteome}}
#' @export
#' @details The background could be generated from wholeGenome, inputSet or nonInputSet.
#' whole genome: randomly select subsequences from the whole genome with each subsequence 
#'  containing amino acids with same width of input sequences.
#'  anchored whole genome: randomly select subsequences from the whole genome with 
#'  each subsequence containing amino acids with same width of input sequences
#'  where the middle amino acids must contain anchor amino acid, e.g., K, which is
#'  specified by user.
#'  input set: same to whole genome, but only use protein sequence from input id and not 
#'  including the site specified in input sequences
#'  anchored input set: same to anchored whole genome, but only use protein sequences from input id, 
#'  and not including the site specified in input sequences.
#'  non-input set: whole genome - input set.
#'  anchored non-input set: whole genome - input set and the middle amino acids must contain anchor amino acid.
#' @return  an object of dagBackground which contains background and permutationSize.
#' @author Jianhong Ou, Alexey Stukalov, Julie Zhu
#' @seealso \code{\link{prepareProteome}}
#' @examples
#'   data("seq.example")
#'   data("proteome.example")
#'   bg <- buildBackgroundModel(seq.example, proteome=proteome.example)
#' @keywords misc

buildBackgroundModel <- function(dagPeptides, 
                                 bg=c("wholeGenome", "inputSet", "nonInputSet"),
                                 model=c("any", "anchored"),
                                 targetPosition=c("any", "Nterminus", "Cterminus"),
                                 uniqueSeq=TRUE,
                                 permutationSize=30L,
                                 rand.seed=1,
                                 replacement=FALSE,
                                 proteome){
    if(missing(dagPeptides) || class(dagPeptides)!="dagPeptides"){
        stop("dagPeptides should be an object of dagPeptides.", call.=FALSE)
    }
    bg <- match.arg(bg)
    targetPosition <- match.arg(targetPosition)
    model <- match.arg(model)
    permutationSize <- as.integer(permutationSize)
    if(permutationSize<2) stop("permutationSize should be greater than 1")
    
    length <- dagPeptides@upstreamOffset + dagPeptides@downstreamOffset + 1
    ###### generate random sequences
    ## TODO, howto remove fetchSequence from background
    if(bg!="inputSet"){
        if(missing(proteome) || class(proteome)!="Proteome"){
            stop("proteome should be an object of Proteome. \n
                 Try ?prepareProteome to get help", call.=FALSE)
        }
        if(bg=="wholeGenome"){
            SequenceStr <- proteome@proteome$SEQUENCE
        }else{
            proteome.s <- proteome@proteome[!proteome@proteome$SEQUENCE %in% dagPeptides@data$peptide,]
            SequenceStr <- proteome.s$SEQUENCE
        }
    }else{
        SequenceStr <- dagPeptides@data$peptide
    }
    if(model=="anchored"){
        anchorAA <- table(dagPeptides@data$anchorAA)
        anchorAA <- anchorAA[order(anchorAA, decreasing=TRUE)]
        if(length(anchorAA)>1){
            model <- "any"
            warning("anchor amino acid is not unique. model is set to any")
            anchorAA <- paste("[",paste(names(anchorAA), collapse=""),"]", sep="")
        }else{
            anchorAA <- names(anchorAA)[1]
            pattern <- paste("([A-Z]{", dagPeptides@upstreamOffset, "}", 
                             anchorAA, "[A-Z]{", dagPeptides@downstreamOffset, "})", sep="")
        }
    }
    if(model=="any"){
        pattern <- paste("([A-Z]{", length ,"})", sep="")
    }
    if(targetPosition=="Cterminus"){
        pattern <- paste(pattern, "$", sep="")
    }else{
        if(targetPosition=="Nterminus"){
            pattern <- paste("^", pattern, sep="")
        }
    }
    matches <- gregexpr(pattern, SequenceStr)
    matches <- unlist(regmatches(SequenceStr, matches))
    set.seed(rand.seed)
    n <- nrow(dagPeptides@data)
    if(length(matches)<n) 
        stop("too less matches in background. Please try different parameters.", 
             call.=FALSE)
    background <- lapply(seq_len(permutationSize), function(p){
        s <- sample(matches, n, replace=replacement, prob=NULL)
        if(uniqueSeq){
            s <- unique(s)
        }
        do.call( rbind, strsplit( s, "", fixed = TRUE ) )
    })
    new("dagBackground", 
        background=background,
        permutationSize=permutationSize)
}