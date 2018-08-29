#################################################
## FASTGenomics I/O for ClEvaR package
#################################################


#' Read FASTGenomics Matrix in HDF5 format.
#'
#' @param filename Name of FASTGenomics HDF5 file.
#' @param chunk Optional 1-based number of chunk to read.
#' @param chunkSize Number of cells per chunk, defaults to 10000.
#' @param fromCell Optional 1-based index of first cell to read.
#' @param toCell Optional 1-based index of last cell to read.
#' @return A sparse \code{Matrix} of expression values.
#' @export
readFGH5 <- function(filename,
                     chunk = NULL,
                     chunkSize = 10000,
                     fromCell = NULL,
                     toCell = NULL) {
  file_to_open <- h5::h5file(filename)
  dataset_names <- h5::list.datasets(file_to_open)
  h5::h5close(file_to_open)
  sparse_names <- c("/matrix/data",
                    "/matrix/indices",
                    "/matrix/indptr",
                    "/obs_names",
                    "/var_names")
  dense_names <- c("/matrix",
                   "/obs_names",
                   "/var_names")
  return_object <- if (all(sparse_names %in% dataset_names)) {
    readFGH5sparse(filename,
                   chunk = chunk,
                   chunkSize = chunkSize,
                   fromCell = fromCell,
                   toCell = toCell)
  } else if (all(dense_names %in% dataset_names)) {
    readFGH5dense(filename,
                  chunk = chunk,
                  chunkSize = chunkSize,
                  fromCell = fromCell,
                  toCell = toCell)
  }
  return(return_object)
}


#' Read FASTGenomics Matrix in sparse HDF5 format.
#'
#' @param filename Name of FASTGenomics HDF5 file.
#' @param chunk Optional 1-based number of chunk to read.
#' @param chunkSize Number of cells per chunk, defaults to 10000.
#' @param fromCell Optional 1-based index of first cell to read.
#' @param toCell Optional 1-based index of last cell to read.
#' @return A sparse \code{Matrix} of expression values.
readFGH5sparse <- function(filename,
                           chunk = NULL,
                           chunkSize = 10000,
                           fromCell = NULL,
                           toCell = NULL) {
  file_to_open <- h5::h5file(filename)
  observations <- file_to_open["/obs_names"][]
  variables <- file_to_open["/var_names"][]
  if (is.null(chunk) & is.null(fromCell) & is.null(toCell)) {
    return_object <- Matrix::sparseMatrix(x = file_to_open["/matrix/data"][],
                                          i = file_to_open["/matrix/indices"][],
                                          p = file_to_open["/matrix/indptr"][],
                                          index1 = FALSE,
                                          dimnames = list(variables, observations))
  } else {
    if (!is.null(fromCell) & !is.null(toCell)) {
      start_index <- fromCell
      stop_index <- toCell
      # what if indexes outside index range?
      cell_index <- seq(start_index, stop_index)
    }
    if (!is.null(chunk)) {
      start_index <- (chunk-1)*chunkSize+1
      stop_index <- min(chunk*chunkSize, length(observations))
      # what if chunk outside index range?
      cell_index <- seq(start_index, stop_index)
    }
    index_pointer <- file_to_open["/matrix/indptr"][c(cell_index[1], cell_index+1)]
    file_index <- seq(min(index_pointer)+1, max(index_pointer))
    expression_index <- file_to_open["/matrix/indices"][file_index]
    return_object <- Matrix::sparseMatrix(x = file_to_open["/matrix/data"][file_index],
                                          i = expression_index,
                                          p = index_pointer-min(index_pointer),
                                          index1 = FALSE)
    detected_index <- unique(expression_index) + 1
    colnames(return_object) <- observations[cell_index]
    rownames(return_object) <- variables[detected_index]
  }
  h5::h5close(file_to_open)
  return(return_object)
}


#' Read FASTGenomics Matrix in dense HDF5 format.
#'
#' @param filename Name of FASTGenomics HDF5 file.
#' @param chunk Optional 1-based number of chunk to read.
#' @param chunkSize Number of cells per chunk, defaults to 10000.
#' @param fromCell Optional 1-based index of first cell to read.
#' @param toCell Optional 1-based index of last cell to read.
#' @return A sparse \code{Matrix} of expression values.
readFGH5dense <- function(filename,
                          chunk = NULL,
                          chunkSize = 10000,
                          fromCell = NULL,
                          toCell = NULL) {
  file_to_open <- h5::h5file(filename)
  observations <- file_to_open["/obs_names"][]
  variables <- file_to_open["/var_names"][]
  cell_index <- seq(1, length(observations))
  if (!is.null(fromCell) & !is.null(toCell)) {
    start_index <- fromCell
    stop_index <- toCell
    # what if indexes outside index range?
    cell_index <- seq(start_index, stop_index)
  }
  if (!is.null(chunk)) {
    start_index <- (chunk-1)*chunkSize+1
    stop_index <- min(chunk*chunkSize, length(observations))
    # what if chunk outside index range?
    cell_index <- seq(start_index, stop_index)
  }
  dim_names <- list(variables, observations[cell_index])
  return_object <- Matrix::Matrix(file_to_open["/matrix"][cell_index, ],
                                  sparse = TRUE,
                                  dimnames = dim_names)
  h5::h5close(file_to_open)
  return(return_object)
}


#' Retrieve observations from FASTGenomics Matrix HDF5 file.
#'
#' Observations typically correspond to cells.
#' @param filename Name of FASTGenomics HDF5 file.
#' @return A vector of observation IDs.
#' @export
observationsFGH5 <- function(filename) {
  file_to_open <- h5::h5file(filename)
  observations <- file_to_open["/obs_names"][]
  h5::h5close(file_to_open)
  return(observations)
}


#' Retrieve variables from FASTGenomics Matrix HDF5 file.
#'
#' Variables typically correspond to genes.
#' @param filename Name of FASTGenomics HDF5 file.
#' @return A vector of variable IDs.
#' @export
variablesFGH5 <- function(filename) {
  file_to_open <- h5::h5file(filename)
  variables <- file_to_open["/var_names"][]
  h5::h5close(file_to_open)
  return(variables)
}


#' Determine the number of chunks of a FASTGenomics Matrix HDF5 file.
#'
#' @param filename Name of FASTGenomics HDF5 file.
#' @param chunkSize Number of cells per chunk, defaults to 10000.
#' @return The number of chunks.
#' @export
maxChunkFGH5 <- function(filename,
                         chunkSize = 10000) {
  observations <- observationsFGH5(filename)
  n_chunks <- ceiling(length(observations)/chunkSize)
  return(n_chunks)
}


#' Read FASTGenomics Autoencoder Output
#' @param filename Name of Autoencoder output file, defaults to \code{encoded_data.h5}.
readAEoutput <- function(filename = "encoded_data.h5") {
  require(h5)
  ftr <- h5file(filename)
  data <- ftr["/matrix"][]
  rownames(data) <- ftr["/obs_names"][]
  colnames(data) <- ftr["/var_names"][]
  h5close(ftr)
  return(data)
}
