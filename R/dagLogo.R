#' Get the data for visualization.
#'
#' A helper function to get the data and the label for visualization.
#' 
#' @param type A character vector of length 1, the type of metrics to display 
#' on y-axis. The available options are "diff" and "statistics", which are 
#' differences in amino acid usage at each position between the inputSet and 
#' the backgroundSet, and the Z-scores or odds ratios when Z-test or Fisher's 
#' exact test is performed to test the differential usage of amino acid at each 
#' position between the two sets.
#' @param testDAUresults An object of \code{\link{testDAUresults-class}} Class (results
#' of testing differential amino acid usage).
#'
#' @return A list containing the following components:
#' label A character vector of length 1. The type of data for visualization.
#' dat   A matrix of numeric data for visualization.
#' @keywords internal
#' 
#' @author Haibo Liu
#' 
getData <- function(type, testDAUresults)
{
  testType <- testDAUresults@testType
  if (type == "diff")
  {
    dat <- testDAUresults@difference
    label <- "DAU %"
  } else
  {
    dat <- testDAUresults@statistics
    if (testType == "fisher")
    {
      label <- "OddsRatio"
    } else
    {
      label <- "Z-score"
    }
  }
  list(label = label, dat = dat)
}



#' Visualize daglogo using a heatmap.
#' 
#' Using a heatmap to visualize results of testing differential amino acid usage.
#' 
#' @param testDAUresults An object of \code{\link{testDAUresults-class}} (results of 
#' testing differential amino acid usage).
#' @param type A character vector of length 1, the type of metrics to display 
#' on y-axis. The available options are "diff" and "statistics", which are 
#' differences in amino acid usage at each position between the inputSet and 
#' the backgroundSet, and the Z-scores or odds ratios when Z-test or Fisher's 
#' exact test is performed to test the differential usage of amino acid at each 
#' position between the two sets.
#' @param ... other parameters passed to the\code{\link[pheatmap]{pheatmap}} function.
#' @import pheatmap
#' @return The output from the \code{\link[pheatmap]{pheatmap}} function.
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' data("seq.example")
#' data("proteome.example")
#' bg <- buildBackgroundModel(seq.example, proteome=proteome.example, 
#'                            numSubsamples=10)
#' t0 <- testDAU(seq.example, bg)
#' dagHeatmap(testDAUresults = t0, type = "diff")

dagHeatmap <-function(testDAUresults, type = c("diff", "statistics"), ...) 
{
  if (missing(testDAUresults) || class(testDAUresults) != "testDAUresults") 
  {
    stop("testDAUresults should be an object of testDAUresults\n
         Please try ?testDAU to get help.", call. = FALSE)
  }
  type <- match.arg(type)
  data <- getData(type, testDAUresults)
  dat <- data$dat[order(rowSums(data$dat), decreasing = TRUE), ]
  pheatmap(dat,
           ...,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           scale = "column",
           main = paste0("Heatmap showing ", data$label))
  }



#' @title retrieve color setting for logo
#' @description retrieve prepared color setting for logo
#' @param colorScheme could be 'null', 'charge', 'chemistry', 'classic' or 'hydrophobicity'
#' @return A character vector of color scheme
#' @author Jianhong Ou
#' @keywords figure

colorsets2 <- function(colorScheme=c("null", "classic", "charge", "chemistry", "hydrophobicity")){
    colorScheme <- match.arg(colorScheme)
    auto<-c('A'='#CCFF00', 'C'='#FFFF00', 'D'='#FF0000', 'E'='#FF0066', 
            'F'='#00FF66', 'G'='#FF9900', 'H'='#0066FF', 'I'='#66FF00', 
            'K'='#6600FF', 'L'='#33FF00', 'M'='#00FF00', 'N'='#CC00FF', 
            'P'='#FFCC00', 'Q'='#FF00CC', 'R'='#0000FF', 'S'='#FF3300', 
            'T'='#FF6600', 'V'='#99FF00', 'W'='#00CCFF', 'Y'='#00FFCC')
    classic <- c("nonpolar_aliphatic"="#000000",
                    "polar_uncharged"="#00811B",
                    "aromatic"="#2000C7",
                    "positively_charged"="#800080",
                    "negatively_charged"="#FFB32C")
    chemistry <- c("hydrophobic"="#000000",
                      "polar"="#00811B",
                      "basic"="#2000C7",
                      "neutral"="#800080",
                      "acidic"='#D00001')
    hydrophobicity <- c("hydrophilic"='#000000', 
                           "neutral"='#00811B', 
                           "hydrophobic"='#2000C7')
    charge<-c("positive"="#FFB32C", "neutral"="#2000C7", "negative"="#CCCCCC")
    switch(colorScheme,
           null=auto,
           classic=classic,
           charge=charge,
           chemistry=chemistry,
           hydrophobicity=hydrophobicity,
           auto)
}

#' Color sets for visualization.
#' 
#' Create color encoding for visualization of a peptide sequence logo.
#'
#' @param colorScheme A character vector of length 1, determining the color scheme
#' based on amino acid classification methods. The available \code{colorScheme} 
#' are ""no","bulkiness_Zimmerman","hydrophobicity_KD", "hydrophobicity_HW", 
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
#' addScheme, with group set as "NULL" or a list with same names as those of \code{color} 
#' and \code{symbol}. No grouping was applied for the first 12 schemes. It is used to 
#' color AAs based on similarities or group amino acids into groups of similarities.. 
#'
#' @seealso {\code{\link{addScheme}}}
#' @return A named character vector of colors
#' @export
#' @keywords internal
#' 
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' colorsets("polarity_Grantham_group")

colorsets <-function(colorScheme = ls(envir = cachedEnv)) 
{
  colorScheme <- match.arg(colorScheme)
  if (!colorScheme %in% ls(envir = cachedEnv))
  {
    stop("Unknown color scheme used!")
  }
  get(colorScheme, envir = cachedEnv)$color
}



#' @title convert group name to a single character
#' @description convert group name to a single character to shown in a logo
#' @param nameScheme could be "classic", "charge", "chemistry", "hydrophobicity"
#' @return A character vector of name scheme
#' @author Jianhong Ou
#' @keywords figure

nameHash <- function(nameScheme=c("classic", "charge", "chemistry", "hydrophobicity")){
    nameScheme <- match.arg(nameScheme)
    classic <- c("nonpolar_aliphatic"="M",
                 "polar_uncharged"="U",
                 "aromatic"="A",
                 "positively_charged"="P",
                 "negatively_charged"="N")
    chemistry <- c("hydrophobic"="H",
                   "polar"="P",
                   "basic"="B",
                   "neutral"="U",
                   "acidic"='A')
    hydrophobicity <- c("hydrophilic"='I', 
                        "neutral"='U', 
                        "hydrophobic"='O')
    charge<-c("positive"="P", "neutral"="U", "negative"="N")
    switch(nameScheme,
           classic=classic,
           charge=charge,
           chemistry=chemistry,
           hydrophobicity=hydrophobicity)
}


#' Get character symbols for grouped amino acids
#'
#' @param groupingScheme A character vector of length 1, determining the character
#' symbols used to represent amino acids grouped by their physical and chemical
#' properties. The available \code{groupingScheme} are "no","bulkiness_Zimmerman",
#' "hydrophobicity_KD", "hydrophobicity_HW", "isoelectric_point_Zimmerman",
#' "contact_potential_Maiorov", "chemistry_property_Mahler", "consensus_similarity_SF", 
#' "volume_Bigelow", "structure_alignments_Mirny", "polarity_Grantham", 
#' "sequence_alignment_Dayhoff", "bulkiness_Zimmerman_group", "hydrophobicity_KD_group",
#' "hydrophobicity_HW_group", "charge_group", "contact_potential_Maiorov_group",
#' "chemistry_property_Mahler_group", "consensus_similarity_SF_group", 
#' "volume_Bigelow_group", "structure_alignments_Mirny_group", "polarity_Grantham_group",  
#' "sequence_alignment_Dayhoff_group", and "custom". If "custom" is used, users 
#' must define a grouping scheme using a list containing sublist named as "color",
#' "symbol" and "group" using the function \code{\link{addScheme}}. No grouping 
#' was applied for the first 12 schemes. 
#'
#' @seealso {\code{\link{addScheme}}}
#'
#' @return A named character vector of character symbols if grouping is applied;
#' otherwise NULL.
#' @export
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu
#'
#' @examples
#' getGroupingSymbol("polarity_Grantham_group")

getGroupingSymbol <-function(groupingScheme = ls(envir = cachedEnv))
{
  groupingScheme <- match.arg(groupingScheme)
  if (!groupingScheme %in% ls(envir = cachedEnv))
  {
    stop("Unknown color scheme used!")
  }
  if (grepl("group", groupingScheme))
  {
    get(groupingScheme, envir = cachedEnv)$symbol
  } else 
  {
    NULL
  }
}


#' Create sequence logos.
#' 
#' Create sequence logos for visualizing results of testing differential usage 
#' of amino acids.
#'
#' @param testDAUresults An object of \code{\link{testDAUresults-class}} (results of
#' testing differential amino acid usage).
#' @param type A character vector of length 1. Type of statistics to display 
#' on y-axis, available choices are "diff" or "zscore".
#' @param pvalueCutoff A numeric vector of length 1. A cutoff of p-values.
#' @param groupingSymbol A named character vector.
#' @param font A character vector of length 1. Font type for displaying sequence
#' Logo.
#' @param legend A logical vector of length 1, indicating whether to show the 
#' legend.
#' @param title A character vector of length 1, main title for a plot.
#' @param labelRelativeToAnchor A logical vector of length 1, indicating whether
#' x-axis label should be adjusted relative to the anchoring position.
#' @param labels A character vector, x-axis labels.
#' @param fontface An integer, fontface of text for axis annotation and legends.
#' @param fontsize An integer, fontsize of text for axis annotation and legends.
#' @param alpha Alpha channel for transparency of low affinity letters.
#' @importFrom grDevices dev.size
#' @importFrom graphics plot.new
#' @import motifStack
#' @import grid
#' @importFrom grImport2 pictureGrob 
#'
#' @return A sequence Logo is plotted without returned values.
#' @export
#' @author Jianhong Ou, Haibo Liu
#'
#' @examples 
#' data('seq.example') 
#' data('proteome.example')
#' bg <- buildBackgroundModel(seq.example, proteome=proteome.example, 
#'                            numSubsamples=10, testType = "ztest")
#' t0 <- testDAU(seq.example, bg)
#' t1 <- testDAU(dagPeptides = seq.example, dagBackground = bg, 
#'               groupingScheme = "hydrophobicity_KD")
#' t2 <- testDAU(dagPeptides = seq.example, dagBackground = bg, 
#'              groupingScheme = "charge_group")
#' t3 <- testDAU(dagPeptides = seq.example, dagBackground = bg, 
#'              groupingScheme = "chemistry_property_Mahler")
#' t4 <- testDAU(dagPeptides = seq.example, dagBackground = bg, 
#'               groupingScheme = "hydrophobicity_KD_group")
#' dagLogo(t0)
#' dagLogo(t1, groupingSymbol = getGroupingSymbol(t1@group))
#' dagLogo(t2, groupingSymbol = getGroupingSymbol(t2@group))
#' dagLogo(t3, groupingSymbol = getGroupingSymbol(t3@group))
#' dagLogo(t4, groupingSymbol = getGroupingSymbol(t4@group))


dagLogo <- function(testDAUresults,
                        type = c("diff", "zscore"),
                        pvalueCutoff = 0.05,
                        groupingSymbol = getGroupingSymbol(testDAUresults@group),
                        font = "Helvetica-Bold",
                        fontface = "bold",
                        fontsize = 5,
                        title = NULL,
                        legend = FALSE,
                        labelRelativeToAnchor = FALSE,
                        labels = NULL, alpha=0.5) 
{
  if (missing(testDAUresults) || class(testDAUresults) != "testDAUresults") 
  {
    stop("testDAUresults should be an object of testDAUresults\n
         Please try ?testDAU to get help.", call. = FALSE)
  }
  type <- match.arg(type)
  data <- getData(type, testDAUresults)
  dat <- data$dat
  npos <- ncol(dat)
  ncha <- nrow(dat)
  if (!is.null(labels)) 
  {
    if (length(labels) < npos) 
    {
      stop("The length of labels is too short!", call. = FALSE)
    }
  }
  colset <- colorsets(testDAUresults@group)
  colset <- colset[rownames(dat)]
  if (any(is.na(colset)))
    stop("Not every symbol has its color setting.", call. = FALSE)
  if (!is.null(groupingSymbol)) 
  {
    rownames(dat) <- groupingSymbol[rownames(dat)]
    names(colset) <- groupingSymbol[names(colset)]
  }
  rname <- rownames(dat)
  if (max(nchar(rname)) > 1) 
  {
    stop("Please using the groupingSymbol to convert",
         " the symbols into single letters.",
         call. = FALSE)
  }
  key <- paste("x",
               ncha,
               font,
               paste(colset[rname], collapse = ""),
               paste(rname, collapse = ""),
               sep = "_")
  
  ## can't use ifelse here
  if(exists("tmp_motifStack_symbolsCache", envir = .globalEnv))
  {
    symbolsCache = get("tmp_motifStack_symbolsCache", envir = .globalEnv)   
  } else
  {
    symbolsCache = list()
  }
  if (!is.null(symbolsCache[[key]])) 
  {
    symbols <- symbolsCache[[key]]
  } else 
  {
    symbols <- motifStack:::coloredSymbols(ncha, font, colset[rname], rname, alpha = alpha)
    symbolsCache[[key]] <- symbols
    
    ## save symbolsCache to the environment variable for future use
    assign("tmp_motifStack_symbolsCache", symbolsCache, envir = .globalEnv)
  }
  
  dw <- ifelse(legend, 1 / (npos + 6), 1 / (npos + 2))
  
  ##check font height to fix the number size
  ## plot.new()  ## This doesn't refresh the setting for the new page of 
  ##                the grid graphics system
  grid.newpage()  ## update the plotting frame for the grid graphics system
  line1 <-
    as.numeric(convertUnit(unit(1, "strwidth", "W"), "npc"))
  cex <- dw / line1
  lwd <- cex / 3
  x0 <- dw
  x1 <- 4 / 5 * dw
  x2 <- 6 / 5 * dw
  x3 <- 1 / 5 * dw
  ##set ylim
  ## below x-axis
  datN <- apply(dat, 2, function(.col)
    sum(.col[.col < 0]))
  
  ## above x-axis
  datP <- apply(dat, 2, function(.col)
    sum(.col[.col > 0]))
  
  ylim <-
    c((as.integer(min(datN) / 0.05) - 1) * 0.05, 
      (as.integer(max(datP) / 0.05) + 1) * 0.05)
  
  remap <- function(x) {
    (ylim[1] - x) / (ylim[1] - ylim[2]) / (1 + dw)
  }
  reheight <- function(x) {
    abs(x / (ylim[1] - ylim[2]) / (1 + dw))
  }
  
  ##draw axis
  x.pos <- 0
  grid.text(0,
            x0,
            remap(0) + dw / 2,
            just = c(.5, .5),
            gp = gpar(fontsize=fontsize, fontface = fontface)
  )
  grid.lines(
    x = x0,
    y = c(remap(0) + dw, 1),
    arrow = arrow(length = unit(0.01, "npc")),
    gp = gpar(lwd = lwd/2)
  )
  grid.lines(
    x = x0,
    y = c(remap(0), 0),
    arrow = arrow(length = unit(0.01, "npc")),
    gp = gpar(lwd = lwd/2)
  )
  
  ##draw tick
  tick <- ifelse(type == "diff", 0.1, 10 ^ as.integer(log10(max(
    abs(as.numeric(dat))))))
  
  times <- ifelse(type == "diff", 100, 1)
  for (i in c(as.integer(min(datN) / tick):(-1), 1:as.integer(max(datP) / tick))) 
  {
    grid.lines(
      x = c(x2, x0),
      y = remap(i * tick) + dw / 2,
      gp = gpar(lwd = lwd/2)
    )
    grid.text(
      label = times * tick * i,
      x = x1,
      y = remap(i * tick) + dw / 2,
      just = c(1, .5),
      gp = gpar(fontsize=fontsize, fontface = fontface)
    )
  }
  
  ## add ylab
  grid.text(label= data$label, 
            x = x3, y = remap(0) + dw / 2,
            just = "centre",  rot = 90, 
            gp = gpar(fontsize=fontsize, fontface = fontface))
  
  x.pos <- x0 + dw / 2
  for (j in 1:npos) {
    heights <- dat[, j]
    id <- order(heights)
    heights <- heights[testDAUresults@pvalue[, j] < pvalueCutoff]
    id <-
      id[id %in% which(testDAUresults@pvalue[, j] < pvalueCutoff)]
    id1 <- order(heights)
    y.pos <- remap(sum(heights[heights < 0]))
    flag <- 0
    if (length(heights) > 0) {
      for (i in 1:length(heights)) {
        h <- reheight(heights[id1[i]])
        if (heights[id1[i]] > 0)
          flag <- flag + 1
        if (flag == 1) 
        {
          ## plot x-axis tick labels
          if (labelRelativeToAnchor) 
          {
            grid.text(
              j - 1 - testDAUresults@upstream,
              x.pos + dw / 2,
              y.pos + dw / 2,
              just = c(.5, .5),
              gp = gpar(fontsize=fontsize, fontface = fontface)
            )
          } else
          {
            if (!is.null(labels)) 
            {
              grid.text(
                labels[j],
                x.pos + dw / 2,
                y.pos + dw / 2,
                just = c(.5, .5),
                gp = gpar(fontsize=fontsize, fontface = fontface)
              )
            } else
            {
              grid.text(
                j,
                x.pos + dw / 2,
                y.pos + dw / 2,
                just = c(.5, .5),
                gp = gpar(fontsize=fontsize, fontface = fontface)
              )
            }
          }
          y.pos <- y.pos + dw
          flag <- flag + 1
        }
        ## plot symbols for amino acids 
        if (h > 0) {
          symid <- ifelse(heights[id1[i]] > 0, id[i], paste0(id[i], "_", alpha))
          grid.draw(
            grImport2::pictureGrob(
              symbols[[symid]],
              x.pos,
              y.pos,
              dw,
              h,
              just = c(0, 0),
              distort = TRUE
            )
          )
          y.pos <- y.pos + h
        }
      }
      
    }
    ## plot x-axis tick labels
    if (flag == 0)
    {
      grid.text(j,
                x.pos + dw / 2,
                y.pos + dw / 2,
                just = c(.5, .5),
                gp = gpar(fontsize=fontsize, fontface = fontface))
    }    
    x.pos <- x.pos + dw
  }
  if (!is.null(title))
  {
    grid.text(
      as.character(title),
      x = 0.5,
      y = 0.98,
      just = c(.5, .5),
      gp = gpar(col = "black", fontsize=fontsize*2, fontface = fontface))
  }
  
  ## plot legend
  if (legend) 
  {
    if (is.null(groupingSymbol)) 
    {
      groupingSymbol <- get(testDAUresults@group, envir = cachedEnv)$symbol
    }
    for (i in 1:length(groupingSymbol)) 
    {
      grid.text(
        groupingSymbol[i],
        x = (npos+2)*dw,
        y = .95 - i * 0.1 * lwd/2,
        just = c(.5, .5),
        gp = gpar(col = colset[groupingSymbol[i]], fontsize=fontsize, 
                  fontface = fontface))
      grid.text(
        names(groupingSymbol)[i],
        x = (npos+2)*dw + 0.1 * lwd/2,
        y = .95 - i * 0.1 * lwd/2,
        just = c(0, .5),
        gp = gpar(col = colset[groupingSymbol[i]], fontsize=fontsize, 
                  fontface = fontface))
    }
  }
  }
