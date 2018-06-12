
#################################################
## Functions for ClEvaR package
#################################################



#################################################
## Mutual Information
#################################################

# code based on https://stackoverflow.com/questions/21831953/r-package-available-for-adjusted-mutual-information

# *** clusterings should be vectors of POSITIVE (>0!!!) integers ***


#' Mutual Information
#' 
#' @param subject Vector of reference cluster assignments.
#' @param query Vector of cluster assignments for comparison.
#' @return 
MI <- function(subject, query) {
  # infotheo license: GPL (>= 3)
  infotheo::mutinformation(X = subject,
                           Y = query,
                           method = "emp")
}


#' Normalized Mutual Information
#' 
#' @param subject Vector of reference cluster assignments.
#' @param query Vector of cluster assignments for comparison.
#' @return 
NMI <- function(subject, query) {
  mi <- MI(subject, query)
  s1 <- tabulate(subject)
  s2 <- tabulate(query)
  N <- length(subject)
  h1 <- -sum(s1*log(s1/N))/N
  h2 <- -sum(s2*log(s2/N))/N
  nmi <- mi/max(h1, h2)
  return(nmi)
}

#' Expected Mutual Information
#' 
#' @param subject Vector of reference cluster assignments.
#' @param query Vector of cluster assignments for comparison.
#' @return 
EMI <- function(subject, query) {
  mi <- MI(subject, query)
  s1 <- tabulate(subject)
  s2 <- tabulate(query)
  N <- l1 <- length(subject)
  l2 <- length(query)
  s_emi <- 0
  for(i in 1:l1){
    for (j in 1:l2){
      min_nij <- max(1,s1[i]+s2[j]-N)
      max_nij <- min(s1[i],s2[j])
      n.ij <- seq(min_nij, max_nij)   #sequence of consecutive numbers
      t1<- (n.ij / N) * log((n.ij * N) / (s1[i]*s2[j]))
      t2 <- exp(lfactorial(s1[i]) + 
                  lfactorial(s2[j]) + 
                  lfactorial(N - s1[i]) +
                  lfactorial(N - s2[j]) - 
                  lfactorial(N) - 
                  lfactorial(n.ij) - 
                  lfactorial(s1[i] - n.ij) - 
                  lfactorial(s2[j] - n.ij) - 
                  lfactorial(N - s1[i] - s2[j] + n.ij))
      emi <- sum(t1*t2)
      s_emi <- s_emi + emi
    }
  }
  return(s_emi)
}

#' Adjusted Mutual Information
#' 
#' @param subject Vector of reference cluster assignments.
#' @param query Vector of cluster assignments for comparison.
#' @return 
AMI <- function(subject, query) {
  mi <- MI(subject, query)
  emi <- EMI(subject, query)
  s1 <- tabulate(subject)
  s2 <- tabulate(query)
  N <- length(subject)
  h1 <- -sum(s1*log(s1/N))/N
  h2 <- -sum(s2*log(s2/N))/N
  ami <- (mi-emi)/max(h1,h2)
  return(ami)
}



#################################################
## Visualization
#################################################

#' Make Donut Plot.
#' 
#' @param subject \code{Vector} of reference cluster assignments.
#' @param query \code{Vector} of cluster assignments for comparison with reference.
#' @param subquery \code{Vector} of lower level cluster assignments for comparison with reference. Defaults to \code{NULL}.
#' @param colorList=NULL A named \code{list} of colors for clusters. The \code{list} needs elements "query" and optionally "subquery".
#' @param pieLim \code{vector} of length 2 giving the \code{xlim} values for pies; defaults to \code{c(0, 4)}.
#' @param pieCut \code{integer} giving the radius of the query pie if subquery is defined; should be within the \code{pieLim} range; defaults to \code{2.5}.
#' @param piesPerRow Number of pies per row. 
#' @return A \code{list} with elements \code{donuts} (a ggplot object) and \code{data} (the underlying data.frame).
#' @examples
#' add(1, 1)
#' add(10, 1)
makeDonuts <- function(subject, 
                       query,
                       subquery = NULL,
                       colorList = NULL,
                       pieLim = c(0, 4),
                       pieCut = 2.5,
                       piesPerRow = ceiling(sqrt(length(unique(subject)))),
                       savePDF = FALSE) {
  donutList <- list()
  subject <- factor(subject)
  query <- factor(query)
  if (is.null(subquery)) {
    if (is.null(colorList)) {
      # define default colorList
      colorList <- list()
      queryColors <- rainbow(length(levels(query)))
      names(queryColors) <- levels(query)
      colorList$query <- queryColors
    }
    colorList$allColors <- colorList$query
    contingencyTable <- t(as.matrix(ftable(subject~query)))
    plotTable <- expand.grid(subject=levels(subject),
                             query=levels(query))
    plotTable$count <- as.vector(contingencyTable)
    plotTable$queryColor <- colorList$query[plotTable$query]
  } else {
    subquery <- factor(subquery)
    if (is.null(colorList)) {
      # define default colorList
      colorList <- list()
      queryColors <- grey.colors(length(levels(query)))
      names(queryColors) <- levels(query)
      colorList$query <- queryColors
      subqueryColors <- rainbow(length(levels(subquery)))
      names(subqueryColors) <- levels(subquery)
      colorList$subquery <- subqueryColors
    }
    colorList$allColors <- c(colorList$query, colorList$subquery)
    contingencyTable <- t(as.matrix(ftable(subject~query+subquery)))
    plotTable <- expand.grid(subject=levels(subject),
                             subquery=levels(subquery),
                             query=levels(query))
    plotTable$count <- as.vector(contingencyTable)
    plotTable$queryColor <- colorList$query[plotTable$query]
    plotTable$subqueryColor <- colorList$subquery[plotTable$subquery]
  }
  # finish table with information required for donut plotting
  plotTable <- subset(plotTable, count!=0)
  plotTable <- plotTable[order(plotTable$subject),]
  plotTable$prop <- plotTable$from <- plotTable$to <- 0
  
  for (i in levels(subject)) {
    index <- plotTable$subject==i
    plotTable$prop[index] <- plotTable$count[index] / sum(plotTable$count[index])
    plotTable$from[index] <- cumsum(c(0, plotTable$prop[index]))[1:length(plotTable$prop[index])]
    plotTable$to[index] <- cumsum(plotTable$prop[index])
  }
  # make ggplot object with donuts
  donutList$donuts <- ggplot2::ggplot(plotTable)
  donutList$donuts <- if (is.null(subquery)) {
    donutList$donuts + ggplot2::geom_rect(ggplot2::aes(fill = query, 
                                                    ymin = from, 
                                                    ymax = to, 
                                                    xmin = min(pieLim), 
                                                    xmax = max(pieLim)))
  } else {
    donutList$donuts + ggplot2::geom_rect(ggplot2::aes(fill = subquery, 
                                                    ymin = from, 
                                                    ymax = to, 
                                                    xmin = pieCut, 
                                                    xmax = max(pieLim))) + 
      ggplot2::geom_rect(ggplot2::aes(fill = query, 
                                      ymin = from, 
                                      ymax = to, 
                                      xmin = min(pieLim), 
                                      xmax = pieCut))
  }
  donutList$donuts <- donutList$donuts +
    ggplot2::coord_polar(theta = "y") + 
    ggplot2::xlim(pieLim) + 
    ggplot2::scale_fill_manual(values = colorList$allColors) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(axis.text = ggplot2::element_blank()) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
    ggplot2::theme(axis.line = ggplot2::element_blank()) + 
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = "") + 
    ggplot2::facet_wrap(~subject, ncol = piesPerRow)
  donutList$data <- plotTable
  return(donutList)
}


#' Plot the legend(s) for donut plots.
#' 
#' @param data \code{data.frame} returned from \code{makeDonuts} in list element \code{data}.
#' @examples
#' add(1, 1)
#' add(10, 1)
plotLegend <- function(data, subquery = FALSE) {
  if (subquery) {
    clusterColumn <- "subquery"
    colorColumn <- "subqueryColor"
  } else {
    clusterColumn <- "query"
    colorColumn <- "queryColor"
  }
  data <- subset(data, !duplicated(data[, clusterColumn]))[, c(clusterColumn, colorColumn)]
  data <- data[order(data[, 1]),]
  legendDf <- data.frame(x = as.numeric(data[,1]),
                         row.names = levels(data[,1]),
                         stringsAsFactors = FALSE)
  colnames(legendDf) <- clusterColumn
  legendColors <- list()
  legendColors[[clusterColumn]] <- data[,2]
  names(legendColors[[clusterColumn]]) <- levels(data[,1])

  pheatmap::pheatmap(legendDf,
                     cellwidth = 10,
                     border_color = NA,
                     color = legendColors[[clusterColumn]],
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     legend = FALSE,
                     labels_col = "")
}


#' Plot donuts.
#' 
#' @param subject \code{Vector} of reference cluster assignments.
#' @param query \code{Vector} of cluster assignments for comparison with reference.
#' @param subquery \code{Vector} of lower level cluster assignments for comparison with reference. Defaults to \code{NULL}.
#' @param savePDF Should plots be saved in PDF? Defaults to \code{FALSE}.
#' @param ... Further parameters used for function \code{makeDonuts}.
#' @examples
plotDonuts <- function(subject, query, subquery = NULL, savePDF = FALSE, ...) {
  require(ggplot2)
  donuts <- makeDonuts(subject, query, subquery, ...)
  if (savePDF) pdf(paste(gsub(":", "_", Sys.time()), "donut plot.pdf"))
  print(donuts$donuts) # the donut plot
  if (!savePDF) x11()
  plotLegend(data=donuts$data, subquery = FALSE)
  if (!is.null(subquery)) {
    if (!savePDF) x11()
    plotLegend(data=donuts$data, subquery = TRUE)
  }
  if (savePDF) dev.off()
}


# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
    layout <- t(layout)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




#' Plot a confusion heatmap.
#' 
#' @param subject \code{Vector} of reference cluster assignments.
#' @param query \code{Vector} of cluster assignments for comparison with reference.
#' @param logCounts Should counts in confusion matrix be \code{log10(counts+1)} transformed? Defaults to \code{FALSE}.
#' @param rankCounts Should counts in confusion matrix be ranked? Defaults to \code{FALSE}.
#' @param ... Further parameters for \code{pheatmap} function.
#' @return A \code{matrix} with possibly transformed counts.
#' @examples
#' a = c(rep("A", 1000), rep("B", 100), rep("C", 10))
#' b = c(rep("A", 500), rep("B", 595), rep("C", 15))
#' confusionHeatmap(a, b)
#' confusionHeatmap(a, b, logCounts = T)
#' confusionHeatmap(a, b, rankCounts = T)
confusionHeatmap <- function(subject,
                             query,
                             logCounts = FALSE,
                             rankCounts = FALSE,
                             ...) {
  confusionMatrix <- as.matrix(ftable(subject~query))
  if (logCounts) confusionMatrix <- log10(confusionMatrix+1)
  if (rankCounts) {
    ranks <- rank(confusionMatrix)
    confusionMatrix <- matrix(ranks,
                              dimnames = dimnames(confusionMatrix),
                              nrow = nrow(confusionMatrix),
                              ncol = ncol(confusionMatrix))
  }
  index <- order(apply(confusionMatrix, 2, max), decreasing = TRUE)
  confusionMatrix <- confusionMatrix[, index]
  index <- order(apply(confusionMatrix, 1, max), decreasing = TRUE)
  confusionMatrix <- confusionMatrix[index,]
  pheatmap(confusionMatrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           ...)
  return(confusionMatrix)
}
