#' Add a custom coloring or grouping scheme.
#' 
#' Add a custom coloring or grouping scheme for ungrouped or grouped amino acids
#' as desired.
#' 
#' @param color A named vector of character. This vector specifies
#' different colors for visualizing the different amino acids or amino acid groups.
#' @param symbol A named vector of character. This vector specifies the 
#' different symbols for visualizing the different amino acids or amino acid groups.
#' @param group A list or NULL. If only coloring amino acids of similar property is
#' desired, set \code{group} to NULL; otherwise \code{group} should be a list with
#' same names as those of \code{color} and \code{symbol}.
#'
#' @return Add the custom coloring or grouping scheme to the environment 
#' \code{cacheEnv}.
#' @export
#' 
#' @examples
#' ## Add a grouping scheme based on the BLOSUM50 level 3 
#' color = c(LVIMC = "#33FF00", AGSTP = "#CCFF00",
#'          FYW = '#00FF66', EDNQKRH = "#FF0066")
#' symbol = c(LVIMC = "L", AGSTP = "A", FYW = "F", EDNQKRH = "E")
#' group = list(
#'    LVIMC = c("L", "V", "I", "M", "C"), 
#'    AGSTP = c("A", "G", "S", "T", "P"),
#'    FYW = c("F", "Y", "W"),
#'    EDNQKRH = c("E", "D", "N", "Q", "K", "R", "H"))
#' addScheme(color = color, symbol = symbol, group = group) 
#' 
 
addScheme <- function(color = vector("character"), 
                     symbol = vector("character"),
                     group = NULL)
{
    if(!is.null(group))
    {
        if (length(color) < 1 || length(symbol) < 1 || length(group) < 1)
        {
            stop("Too few groups!")
        }
        if (length(color) != length(symbol) || length(color) != length(group) ||
            length(symbol) != length(group))
        {
            stop("Wrong grouping specification:", 
                 "The length of color, symbol and group should be the same!")
        }
        if(any(names(color) != names(symbol)) || any(names(color) != names(group)) ||
           any(names(symbol) != names(group)))
        {
            stop("The names of color, symbol and group should be the same!") 
        }
        custom_group <- list(color = color, symbol = symbol, group = group)
        assign("custom_group", custom_group, envir = cachedEnv)
    } else {
        if (length(color) != 20 || length(symbol) != 20 )
        {
            stop("The length of colors and symbols should be 20.")
        }
        custom <- list(color = color, symbol = symbol, group = NULL)
        assign("custom", custom, envir = cachedEnv)
    }
}


#' Differential usage test of amino acids or amino acid groups.
#'
#' Test differential usage of amino acids with or without grouping between
#' experimental sets and background sets.
#'
#' @param dagPeptides An object of Class \code{\link{dagPeptides-class}}.
#' @param dagBackground An object of Class \code{\link{dagBackground-class}}.
#' @param groupingScheme A character vector of length 1. Available choices are 
#' "no","bulkiness_Zimmerman","hydrophobicity_KD", "hydrophobicity_HW", 
#' "isoelectric_point_Zimmerman", "contact_potential_Maiorov",
#' "chemistry_property_Mahler", "consensus_similarity_SF", 
#' "volume_Bigelow", "structure_alignments_Mirny", "polarity_Grantham", 
#' "sequence_alignment_Dayhoff", "bulkiness_Zimmerman_group", "hydrophobicity_KD_group",
#' "hydrophobicity_HW_group", "charge_group", "contact_potential_Maiorov_group",
#' "chemistry_property_Mahler_group", "consensus_similarity_SF_group", 
#' "volume_Bigelow_group", "structure_alignments_Mirny_group", "polarity_Grantham_group",  
#' "sequence_alignment_Dayhoff_group", "custom" and "custom_group". If "custom" or
#' "custom_group" are used, users must define a grouping scheme using a list 
#' containing sublist named as "color", and "symbol" using the function
#' addScheme, with group set as "NULL" or a list with same names as those of color 
#' and symbol. No grouping was applied for the first 12 schemes. It is used to 
#' color AAs based on similarities or group amino acids into groups of similarities.
#' 
#' @param bgNoise A numeric vector of length 1 if not NA. It should be in 
#' the interval of (0, 1) when not NA.
#'
#' @return An object of Class \code{\link{testDAUresults-class}}.
#' @importFrom stats pnorm rgamma sd fisher.test
#' @import methods
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' dat <- unlist(read.delim(system.file(
#'                                    "extdata", "grB.txt", package = "dagLogo"),
#'                          header = FALSE, as.is = TRUE))
#'                          
#' ##prepare an object of Proteome Class from a fasta file
#' proteome <- prepareProteome(fasta = system.file("extdata",
#'                                                 "HUMAN.fasta",
#'                                                 package = "dagLogo"), 
#'                             species = "Homo sapiens")
#' ##prepare an object of dagPeptides Class
#' seq <- formatSequence(seq = dat, proteome = proteome, upstreamOffset = 14,
#'                      downstreamOffset = 15)
#' bg_fisher <- buildBackgroundModel(seq, background = "wholeProteome",
#'                                   proteome = proteome, testType = "fisher")
#' bg_ztest <- buildBackgroundModel(seq, background = "wholeProteome",
#'                                    proteome = proteome, testType = "ztest")
#' 
#' ## no grouping and distinct coloring scheme
#' t0 <- testDAU(seq, dagBackground = bg_ztest)
#' 
#' ## grouped by polarity index (Granthm, 1974)
#' t1 <- testDAU(dagPeptides = seq, dagBackground = bg_ztest, 
#'              groupingScheme = "polarity_Grantham_group")
#'              
#' ## grouped by charge.
#' t2 <- testDAU(dagPeptides = seq, dagBackground = bg_ztest, 
#'               groupingScheme = "charge_group")
#'               
#' ## grouped on the basis of the chemical property of side chains.
#' t3 <- testDAU(dagPeptides = seq, dagBackground = bg_ztest, 
#'               groupingScheme = "chemistry_property_Mahler_group")
#'               
#' ## grouped on the basis of hydrophobicity (Kyte and Doolittle, 1982)
#' t4 <- testDAU(dagPeptides = seq, dagBackground = bg_ztest, 
#'               groupingScheme = "hydrophobicity_KD_group")                                   


testDAU <- function(dagPeptides,
                    dagBackground,
                    groupingScheme = ls(envir = cachedEnv),
                    bgNoise = NA) 
{
    if (missing(dagPeptides) || class(dagPeptides) != "dagPeptides") 
    {
        stop("dagPeptides should be an object of dagPeptides class.\n
            Please try ?fetchSequence to get help.",
            call. = FALSE)
    }
    if (missing(dagBackground) || class(dagBackground) != "dagBackground") 
    {
        stop("dagBackground should be an object of dagBackground class .\n
            Please try ?buildBackgroundModel to get help.",
            call. = FALSE)
    }
    if (!is.na(bgNoise)) 
    {
        if (bgNoise < 0 || bgNoise > 1)
        {
            stop("The background noise is a number in the range of 0 to 1", 
                 call. = FALSE)
        }
    }

    exp <- dagPeptides@peptides
    bg <- dagBackground@background
    
    groupingScheme <- match.arg(groupingScheme)
    if (!groupingScheme %in% ls(envir = cachedEnv))
    {
        stop("Unknown grouping scheme used!")
    }
    
    group <- get(groupingScheme, envir = cachedEnv)$group
    if (is.null(group))
    {
        AA <- get(groupingScheme, envir = cachedEnv)$symbol
        coln <- as.character(AA)
    } else
    {
        coln <- names(get(groupingScheme, envir = cachedEnv)$group)
    }
    ## helper function used to group AAs based on groupingScheme
    convert <- function(x, gtype) 
    {
        d <- dim(x)
        x <- as.character(x)
        for (i in 1:length(gtype)) 
        {
            id <- x %in% gtype[[i]]
            name <- names(gtype)[i]
            x[id] <- name
        }
        matrix(x, nrow = d[1], ncol = d[2])
    }
    groupAA <- function(dat, group) 
    {
        if (grepl("group", group))
        {
           dat <- convert(dat, get(group, envir = cachedEnv)$group)
        }
        dat
    }
    
    bg <- lapply(bg, function(.bg) groupAA(.bg, groupingScheme))
    
    exp <- groupAA(exp, groupingScheme)
    if (ncol(exp) != ncol(bg[[1]]))
    {
        stop("The length of background is different from inputs", call. = FALSE)
    }
    
    ## get counts for each amino acid at each position if Fisher's exact test
    ## otherwise get proportions
    testType <- dagBackground@testType
    counts <- function(mat, coln, testType) 
    {
        num <- apply(mat, 2, function(.ele) {
            cnt <- table(.ele)[coln]
            ## just in case not all coln in the dataset
            names(cnt) <- coln 
            cnt[is.na(cnt)] <- 0
            total <- sum(cnt)
            list(freq = cnt, percent = cnt / total)})
    }
    ## a list of list
    bg <- lapply(bg, counts, coln, testType = testType)
 
    ## a list
    exp <- counts(exp, coln, testType = testType)
    ## get frequency of each amino acid at each position
    exp_freq <- do.call(cbind, lapply(exp, function(.ele){
       .ele$freq
    }))
    exp_freq[is.na(exp_freq)] <- 0
    ## get percentage of each amino acid at each position
    exp_percent <- do.call(cbind, lapply(exp, function(.ele){
        .ele$percent
    }))
    exp_percent[is.na(exp_percent)] <- 0
    rownames(exp_freq) <- rownames(exp_percent ) <- coln
    
    ## Add background noise to the background model
    if(testType == "z-test")
    {
        per_sample_bg_percent <- lapply(bg, function(.bg) {
                do.call(cbind, lapply(.bg, function(.ele){
                    .ele$percent
                }))
            })
        per_position_bg_percent <- lapply(1:ncol(exp_percent),function(i) {
            do.call(cbind, lapply(per_sample_bg_percent, function(.ele){
            .ele[, i]
            }))
        })
        
        if (!is.na(bgNoise)) 
        {
            rdirichlet <- function (n, alpha) {
                l <- length(alpha)
                x <-matrix(rgamma(l * n, alpha),
                           ncol = n,
                           byrow = TRUE)
                sm <- x %*% rep(1, n)
                t(x / as.vector(sm))
            }
            bg_percent <- lapply(per_position_bg_percent, function(col_pct) {
                (1 - bgNoise) * col_pct + bgNoise * 
                    rdirichlet(nrow(col_pct), rep(1, ncol(col_pct)))
            })
        }
        
        ## get mean proportion of each amino acid at each position
        mu_percent <- do.call(cbind, lapply(bg_percent, function(.bg) {
            apply(.bg, 1, mean, na.rm = TRUE)
        }))
        
        diff_percent <- exp_percent - mu_percent
        diff_percent[is.na(diff_percent)] <- 0
        n <-  dagBackground@numSubsamples
        std_percent <- sqrt(mu_percent*(1-mu_percent)/n)
        statistics <- diff_percent / std_percent
        pvalue <- 2 * pnorm(-abs(statistics))
    } else 
    {
        bg <- bg[[1]]
        bg_freq <- do.call(cbind, lapply(bg, function(.ele){
            .ele$freq
        }))
        mu_percent <- do.call(cbind, lapply(bg, function(.ele){
            .ele$percent
        }))

        diff_percent <- exp_percent - mu_percent
        
        ## Fisher exact test
        non_bg_freq <- -sweep(x= bg_freq, MARGIN= 2, STATS = colSums(bg_freq), FUN = "-") 
        non_exp_freq <- -sweep(x= exp_freq, MARGIN= 2, STATS = colSums(exp_freq), FUN = "-")
        
        testResults <- mapply(function(x, X, y, Y){
            fisher <- fisher.test(matrix(c(x, X, y, Y), nrow =2, byrow =TRUE),
                                  alternative = "two.sided")
            c(pvalues=fisher$p.value, statistics = fisher$estimate)
        }, exp_freq, non_exp_freq, bg_freq, non_bg_freq, SIMPLIFY = FALSE)
        
        testOut <- as.data.frame(do.call(rbind, testResults))
        colnames(testOut) <- c("pvalues", "statistics")
        pvalue <- matrix(testOut$pvalues, nrow = nrow(exp_percent), byrow =FALSE)
        statistics <- matrix(testOut$statistics, nrow = nrow(exp_percent), byrow =FALSE)
    }
    rownames(diff_percent) <- rownames(pvalue) <- rownames(statistics) <- coln
   
    coln <- c()
    if (dagPeptides@upstreamOffset > 0) 
    {
        coln <- paste("AA", -1 * (dagPeptides@upstreamOffset:1), sep = "")
    }
    coln <- c(coln, "AA0")
    if (dagPeptides@downstreamOffset > 0) 
    {
        coln <- c(coln,
                  paste("AA", 1:dagPeptides@downstreamOffset, sep = ""))
    }
    if (length(coln) == ncol(diff_percent)) 
    {
        colnames(diff_percent) <- colnames(statistics) <- colnames(pvalue) <- coln
    } else
    {
        coln <- paste("AA", 1:ncol(diff_percent), sep = "")
        colnames(diff_percent) <- colnames(statistics) <- colnames(pvalue) <- coln
    }

    ## return the test results as an object of testDAUresults Class
    new("testDAUresults",
        group = groupingScheme,
        difference = diff_percent,
        statistics = statistics,
        pvalue = pvalue,
        background = mu_percent,
        motif = exp_percent,
        testType = dagBackground@testType,
        upstreamOffset = dagPeptides@upstreamOffset,
        downstreamOffset = dagPeptides@downstreamOffset)
}