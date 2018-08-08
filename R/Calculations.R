
#################################################
## Calculations for ClEvaR package
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


#' Compute the HBITM similarity between two clusterings.
#'
#' @param subject Vector of reference cluster assignments.
#' @param query Vector of cluster assignments for comparison.
#' @param weighted Use weights? Defaults to \code{TRUE}.
#' @return The (weighted) HBITM ("How blue is the matrix?") measure, a numeric vector of length 1.
#' @examples
#' a = c(rep("A", 1000), rep("B", 100), rep("C", 10))
#' b = c(rep("A", 500), rep("B", 595), rep("C", 15))
#' HBITM(subject = a, query = b, weighted = TRUE)
#' HBITM(subject = a, query = b, weighted = FALSE)
#'
#' data(zeisel2018)
#' HBITM(subject = zeisel2018$cell_metadata$class,
#'       query = zeisel2018$cell_metadata$cluster_name,
#'       weighted = TRUE)
#' HBITM(subject = zeisel2018$cell_metadata$class,
#'       query = zeisel2018$cell_metadata$cluster_name,
#'       weighted = FALSE)
#' @export
HBITM <- function(subject,
                  query,
                  weighted = TRUE) {
  m <- as.matrix(ftable(subject~query))
  if (ncol(m)>nrow(m)) m <- t(m)
  contrast_weights <- if (weighted) {
    rowSums(m)
  } else {
    rep(1, nrow(m))
  }
  max_index <- apply(m , 1, which.max)
  max_values <- sapply(1:nrow(m), function(x, m, i) m[x, i[x]],
                       m = m,
                       i = max_index)
  not_max_sum <- sapply(1:nrow(m), function(x, m, i) mean(m[x, -i[x]]),
                        m = m,
                        i = max_index)
  row_contrast <- (max_values-not_max_sum)/(max_values+not_max_sum)
  contrast_sum <- sum(row_contrast * contrast_weights)
  best_sum <- sum(contrast_weights)
  return(contrast_sum/best_sum)
}


#################################################
## Cluster Integrity
#################################################

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
## Cell Type Assignment
#################################################

#' Classify cells using marker genes.
#'
#' @param dataset A FASTGenomics Matrix of expression values.
#' @param markerGenes An optional data frame with scores of marker genes for cell types; the immune cell markers compiled from Newman et al., Nature Methods (2015) is provided as default.
#' @param geneColumn Gene column name in markerGenes data frame; defaults to \code{entrez_id}.
#' @param scoreColumn Score column name in markerGenes data frame; defaults to \code{score}.
#' @param classColumn Class column name in markerGenes data frame; defaults to \code{cell_type}.
#' @return A data frame with class assignment and scores for classes per cell.
#' @export
scoreCells <- function(dataset,
                       markerGenes = NULL,
                       geneColumn = "entrez_id",
                       scoreColumn = "score",
                       classColumn = "cell_type") {
  if (is.null(markerGenes)) {
    data("newman2015")
    markerGenes <- newman2015
  }
  signatures <- markerGenes[markerGenes[, geneColumn] %in% rownames(dataset),]
  signatures <- split(signatures, f = signatures[, classColumn])
  assignStats <- data.frame(cell_id = colnames(dataset),
                            stringsAsFactors = FALSE)
  for (i in names(signatures)) {
    assignStats[, i] <- apply(dataset[signatures[[i]][, geneColumn],], 2,
                              function(x, s) sum(x*s), s = signatures[[i]][, scoreColumn])
  }
  tmp <- apply(assignStats[, -1], 1, function(x) which(x==max(x))[1])
  assignStats[, classColumn] <- colnames(assignStats)[tmp+1]
  return(assignStats)
}
