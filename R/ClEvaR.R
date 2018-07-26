
#################################################
## Functions for ClEvaR package
#################################################



#################################################
## Mutual Information
#################################################

# code based on https://stackoverflow.com/questions/21831953/r-package-available-for-adjusted-mutual-information

# *** clusterings should be vectors of POSITIVE (>0!!!) integers ***

#' Calculate Shannon Entropy
#'
#' @param subject Vector of reference cluster assignments.
#' @return Mutual Information, a numeric vector of length 1.
#' @export
entropy <- function(subject) {
  s1 <- tabulate(subject)
  N <- length(subject)
  return(-sum(s1*log(s1/N))/N)
}

#' Mutual Information
#'
#' @param subject Vector of reference cluster assignments.
#' @param query Vector of cluster assignments for comparison.
#' @return Mutual Information, a numeric vector of length 1.
#' @export
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
#' @return Normalized Mutual Information, a numeric vector of length 1.
#' @export
NMI <- function(subject, query) {
  mi <- MI(subject, query)
  es <- entropy(subject)
  eq <- entropy(query)
  nmi <- mi/max(es, eq)
  return(nmi)
}

#' Expected Mutual Information
#'
#' @param subject Vector of reference cluster assignments.
#' @param query Vector of cluster assignments for comparison.
#' @return Expected Mutual Information, a numeric vector of length 1.
EMI <- function(subject, query) {
  mi <- MI(subject, query)
  s1 <- tabulate(subject)
  s2 <- tabulate(query)
  N <- length(subject)
  l1 <- length(s1)
  l2 <- length(s2)
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
#' @return Adjusted Mutual Information, a numeric vector of length 1.
#' @export
AMI <- function(subject, query) {
  mi <- MI(subject, query)
  emi <- EMI(subject, query)
  es <- entropy(subject)
  eq <- entropy(query)
  ami <- (mi-emi)/max(es, eq)
  return(ami)
}


#' Cluster Granularity
#'
#' Computes the proportion of query clusters assigned to each subject cluster.
#' @param subject Vector of reference cluster assignments.
#' @param query Vector of cluster assignments for comparison.
#' @param plot Display the query x subject matrix; defaults to \code{FALSE}.
#' @return Vector of cluster granularities for each subject cluster.
#' @export
clusterGranularity <- function(subject, query, plot = FALSE) {
  counts <- as.matrix(ftable(subject~query))
  rowProps <- counts/rowSums(counts)
  if (plot) pheatmap::pheatmap(rowProps)
  return(colSums(rowProps))
}


#' Cluster Integrity
#'
#' Computes if query clusters are exclusively assigned to subject cluster.
#' @param subject Vector of reference cluster assignments.
#' @param query Vector of cluster assignments for comparison.
#' @param maxIgnore Maximum count that is ignored in contingency matrix; defaults to 0.
#' @return Vector of cluster integrity for each subject cluster.
#' @export
clusterIntegrity <- function(subject, query, maxIgnore = 0) {
  counts <- as.matrix(ftable(subject~query))
  rowProps <- counts/rowSums(counts)
  propSums <- colSums(rowProps)
  clusterCount <- colSums(counts>maxIgnore)
  return(propSums/clusterCount)
}


#################################################
## FASTGenomics I/O
#################################################

#' Read expression data from FASTGenomics dense HDF5 file.
#'
#' @param fileName Name of FASTGenomics HDF5 file.
#' @param cellIDs Optional vector of cell IDs (i.e. row names) to be read from FASTGenomics HDF5 file.
#' @param geneIDs Optional vector of gene IDs (i.e. column names) to be read from FASTGenomics HDF5 file.
#' @param matrixName Dataset name of expression matrix in HDF5 file; defaults to "/matrix".
#' @param observationsName Dataset name of observation (i.e. cell) name vector in HDF5 file; defaults to "/obs_names".
#' @param variablesName Dataset name of variable (i.e. gene) name vector in HDF5 file; defaults to "/var_names".
#' @return A sparse Matrix of expression values in cells x genes format.
#' @export
readFGH5dense <- function(fileName,
                          cellIDs = NULL,
                          geneIDs = NULL,
                          matrixName = "/matrix",
                          observationsName = "/obs_names",
                          variablesName = "/var_names") {
  tmpH5 <- h5::h5file(fileName)
  cellIndex <- if (is.null(cellIDs)) {
    rep(TRUE, length(tmpH5[observationsName][]))
  } else {
    tmpH5[observationsName][] %in% cellIDs
  }
  geneIndex <- if (is.null(geneIDs)) {
    rep(TRUE, length(tmpH5[variablesName][]))
  } else {
    tmpH5[variablesName][] %in% geneIDs
  }
  cellIndex <- which(cellIndex)
  geneIndex <- which(geneIndex)
  exprs <- Matrix::Matrix(tmpH5[matrixName][cellIndex, geneIndex],
                          sparse = TRUE)
  rownames(exprs) <- tmpH5[observationsName][cellIndex]
  colnames(exprs) <- tmpH5[variablesName][geneIndex]
  h5::h5close(tmpH5)
  return(exprs)
}



#################################################
## Differentially Expressed Genes
#################################################

#' Select Marker Genes
#'
#' Select potential marker genes from differentially expressed genes on effect size and presence in other clusters
#'
#' @param FGDEGtab The raw FASTGenomics Differentially Expressed Gene table, e.g. the output from FASTGenomics calc_de_genes_nonparametric.
#' @param minES Minimal effect size for gene to be considered; defaults to 0.5.
#' @param maxClustersPerGene Maximum number of clusters that contain a potential marker gene; defaults to 1.
#' @param clusterColumn Column name for cluster assignment; defaults to "cluster_id".
#' @param geneColumn Column name for gene ID; defaults to "entrez_id".
#' @param effectColumn Column name for effect size; defaults to "effect.size".
#' @return The filtered FASTGenomics Differentially Expressed Gene table.
#' @export
selectMarkerGenes <- function(FGDEGtab,
                              minES = 0.5,
                              maxClustersPerGene = 1,
                              clusterColumn = "cluster_id",
                              geneColumn = "entrez_id",
                              effectColumn = "effect.size") {
  ESfiltered <- FGDEGtab[FGDEGtab[, effectColumn]>=minES, ]
  allGenes <- unique(ESfiltered[, geneColumn])
  clusterCounts <- sapply(allGenes, function(x, y) sum(x==y), y=ESfiltered[, geneColumn])
  index <- match(ESfiltered[, geneColumn], allGenes)
  ESfiltered$clusterHits <- clusterCounts[index]
  ESHitsFiltered <- subset(ESfiltered, clusterHits<=maxClustersPerGene)
  return(ESHitsFiltered)
}


#' Characterize Marker Thresholds
#'
#' Construct a matrix of clusters with marker genes captured with
#' ranges of effect size and cluster hit threshold values.
#' The resulting matrix should help at selecting optimal threshold values.
#'
#' @param FGDEGtab The raw FASTGenomics Differentially Expressed Gene table, e.g. the output from FASTGenomics calc_de_genes_nonparametric.
#' @param esThresh A vector of minimal effect sizes; defaults to \code{seq(0.5, 1, .01)}.
#' @param nThresh An optional vector of maximum cluster numbers that contain a potential marker gene; defaults to NULL, which takes the range between 1 and the maximum number of clusters found in the FGDGEtab.
#' @param clusterColumn Column name for cluster assignment; defaults to "cluster_id".
#' @param geneColumn Column name for gene ID; defaults to "entrez_id".
#' @param effectColumn Column name for effect size; defaults to "effect.size".
#' @return The marker threshold matrix.
#' @export
characterizeMarkerThresholds <- function(FGDEGtab,
                                         esThresh = seq(0.5, 1, .01),
                                         nThresh = NULL,
                                         clusterColumn = "cluster_id",
                                         geneColumn = "entrez_id",
                                         effectColumn = "effect.size") {
  if (is.null(nThresh)) {
    nClust <- length(unique(FGDEGtab[, clusterColumn]))
    propClust <- round(0.1*nClust)
    maxClust <- max(2, propClust)
    nThresh <- 1:maxClust
  }
  thresholdMatrix <- matrix(data = 0,
                            nrow = length(esThresh),
                            ncol = length(nThresh))
  for (i in 1:length(esThresh)) {
    for (j in 1:length(nThresh)) {
      tmp <- selectMarkerGenes(FGDEGtab = FGDEGtab,
                               minES = esThresh[i],
                               maxClustersPerGene = nThresh[j],
                               clusterColumn = clusterColumn,
                               geneColumn = geneColumn,
                               effectColumn = effectColumn)
      thresholdMatrix[i, j] <- length(unique(tmp[, clusterColumn]))
    }
  }
  rownames(thresholdMatrix) <- paste("ES=", esThresh, sep ="")
  colnames(thresholdMatrix) <- paste(nThresh, "cluster")
  return(thresholdMatrix)
}


#' Find Marker Thresholds
#'
#' Find optimized values for minES and maxClustersPerGene in marker threshold matrix.
#' Assumes that most stringent parameter combination is in the bottom left corner
#' and the least stringent combination is in the upper right corner of the
#' marker threshold matrix.
#'
#' @param MTM The marker threshold matrix.
#' @param targetValue The minimum target value that will be searched the the marker threshold matrix; defaults to \code{max(MTM)}.
#' @return A list with elements minES and maxClustersPerGene.
findMarkerThresholds <- function(MTM,
                                 targetValue = max(MTM)) {
  # flip the matrix so that high stringency has low index positions
  MTM <- MTM[seq(nrow(MTM), 1),]
  x <- apply(MTM, 2, function(x, t) any(x>=t), t = targetValue)
  y <- suppressWarnings(apply(MTM, 2, function(x, t) max(which(x>=t)), t = targetValue))
  hits_x <- which(x)
  hits_y <- y[hits_x]
  optimal_xy_pos <- which.min(hits_x+hits_y)
  tradeoff <- list(minES = as.numeric(sub("ES=", "", rownames(MTM)[hits_y[optimal_xy_pos]])),
                   maxClustersPerGene = as.numeric(sub(" cluster", "", colnames(MTM)[hits_x[optimal_xy_pos]])))
  return(tradeoff)
}


#' Construct a matrix of expression summary values per cluster.
#'
#' @param FGDEGtab The raw FASTGenomics Differentially Expressed Gene table, e.g. the output from FASTGenomics calc_de_genes_nonparametric.
#' @param FGexprs A sparse Matrix of expression values in cells x genes format.
#' @param FGassign data.frame with cluster assignments for cells.
#' @param Q Quantile used for expression value summary; defaults to 0.5, i.e. the median.
#' @param cellColumn Column name for cell name; defaults to "cell_id".
#' @param geneColumn Column name for gene ID; defaults to "entrez_id".
#' @param clusterColumn Column name for cluster assignment; defaults to "cluster_id".
#' @return A gene x cluster matrix of summarized cluster expression values.
#' @export
clusterExprsSummaryMatrix <- function(FGDEGtab,
                                      FGexprs,
                                      FGassign,
                                      Q = 0.5,
                                      cellColumn = "cell_id",
                                      geneColumn = "entrez_id",
                                      clusterColumn = "cluster_id") {
  FGassign <- FGassign[FGassign[, cellColumn] %in% rownames(FGexprs),]
  clusterAssignments <- split(x = as.character(FGassign[, cellColumn]),
                              f = FGassign[, clusterColumn])
  deGenes <- as.character(unique(FGDEGtab[, geneColumn]))
  deGenes <- deGenes[deGenes %in% colnames(FGexprs)]
  FGexprs <- FGexprs[, deGenes]
  FGexprs <- sapply(clusterAssignments, function(ca, M) M[ca,], M=FGexprs)
  qMat <- sapply(FGexprs, function(m, q) apply(m, 2, quantile, probs=q), q=Q)
  return(qMat)
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
#' @export
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


#' Multiple Plot Function
#'
#' Taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @param ... ggplot obejects to plot
#' @param plotlist An optional list of ggplot objects.
#' @param cols Number of columns in layout; defaults to 1.
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
#' @export
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
#' @export
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
  pheatmap::pheatmap(confusionMatrix,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     ...)
  return(confusionMatrix)
}
