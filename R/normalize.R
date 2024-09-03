#'
#' Normalization
#'
#' @description
#' Apply a specified type of normalization to the data.
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param method A string (default = "quant") specifying the method of normalization to
#' apply:
#' \enumerate{
#' \item Row-wise (sample-wise) normalization:
#' \itemize{
#' \item "mean": Samples are normalized to the mean value of all variables for a given
#' sample.
#' \item "median": samples are normalized to the median value of all variables for a given
#' sample.
#' \item "quantile": Quantile normalization.
#' }
#' \item Column-wise (metabolite-wise) normalization:
#' \itemize{
#' \item "meanCtr": Mean centering.
#' \item "medianCtr": Median centering.
#' \item "auto": Auto scaling.
#' \item "level": Level scaling.
#' \item "pareto": Pareto scaling.
#' \item "range": Range scaling.
#' \item "vast": Vast scaling.
#' }
#' }
#'
#' @details
#' Quantile normalization is generally recommended. Mean and median normalization are
#' going to be included as popular previous methods. No normalization is not recommended.
#' Boxplots are also generated for before and after the normalization to give a visual
#' indicator of the changes.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom limma normalizeQuantiles
#'
#' @return The normalized data.
#'
#' @autoglobal
#'
#' @export

normalize <- function(dataSet, method = "quantile") {

  ## create a boxplot for pre-normalization
  plot <- normPlot(dataSet = dataSet) +
    ggtitle("Pre-Normalization Boxplot")
  print(plot)

  attrnames <- attributes(dataSet)$attrnames

  ## select the numerical data
  dataPoints <- select(dataSet, -all_of(attrnames))

  if (method == "meanCtr") {

    ## mean centering
    ## subtract the average intensity of each compound from every value for that compound
    normDataPoints <- base::scale(dataPoints, center = TRUE, scale = FALSE)

  } else if (method == "medianCtr") {

    ## median centering
    ## subtract the median intensity of each compound from every value for that compound
    normDataPoints <- base::scale(dataPoints,
                                  center = apply(dataPoints, 2, median, na.rm = TRUE),
                                  scale = FALSE)

  } else if (method == "auto") {

    ## auto scaling
    normDataPoints <- base::scale(dataPoints, center = TRUE, scale = TRUE)

  } else if (method == "level") {

    ## level scaling
    normDataPoints <- base::scale(dataPoints, center = TRUE,
                                  scale = apply(dataPoints, 2, mean, na.rm = TRUE))

  } else if (method == "pareto") {

    ## Pareto scaling
    sigma <- apply(base::scale(dataPoints, center = TRUE, scale = FALSE), 2, function(col) {
      col <- col[!is.na(col)]
      sqrt(sum(col^2) / max(1L, length(col)-1L))
    })
    normDataPoints <- base::scale(dataPoints, center = TRUE, scale = sqrt(sigma))

  } else if (method == "range") {

    ## range scaling
    range_diff <- apply(dataPoints, 2, function(col) {
      res <- max(col, na.rm = TRUE) - min(col, na.rm = TRUE)
      res <- ifelse(res == 0, 1, res)
      return(res)
    })
    normDataPoints <- base::scale(dataPoints, center = TRUE, scale = range_diff)

  } else if (method == "vast") {

    ## vast scaling
    mu <- colMeans(dataPoints, na.rm = TRUE)
    sigma <- apply(base::scale(dataPoints, center = TRUE, scale = FALSE), 2, function(col) {
      col <- col[!is.na(col)]
      sqrt(sum(col^2) / max(1L, length(col)-1L))
    })
    normDataPoints <- base::scale(dataPoints, center = TRUE, scale = sigma^2/mu)

  } else if (method == "mean") {

    ## normalization by mean
    normDataPoints <- base::scale(t(dataPoints), center = TRUE, scale = FALSE)

  } else if (method == "median") {

    ## normalization by median
    normDataPoints <- base::scale(t(dataPoints),
                                  center = apply(dataPoints, 1, median, na.rm = TRUE),
                                  scale = FALSE)

  } else if (method == "quantile") {

    ## quantile normalization
    normDataPoints <- t(limma::normalizeQuantiles(t(dataPoints)))

  }

  ## recombine the labels and transformed data into a single data frame
  normDataSet <- cbind(dataSet[,attrnames], normDataPoints)
  attributes(normDataSet)$attrnames <- attrnames

  ## create a boxplot for post-normalization
  plot <- normPlot(normDataSet) +
    ggtitle("Post-Normalization Boxplot")
  print(plot)

  ## return pre-processed data
  return(normDataSet)
}
