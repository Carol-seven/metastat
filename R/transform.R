#'
#' Plotting a graph of mean versus variance
#'
#' @description
#' Take a set of metabolics data organized by column, calculate the mean and variance of
#' each column, and then plot those statistics.
#'
#' @param datMV A data frame containing the data signals.
#'
#' @param title A string with the desired title for the mean-variance plot.
#'
#' @import ggplot2
#' @importFrom stats var median
#'
#' @returns An object of class \code{plot}.
#'
#' @autoglobal
#'
#' @noRd

meanVarPlot <- function(datMV, title = "") {

  ## calculate the mean and variance for each protein individually
  plotData <- data.frame(t(sapply(datMV, function(x) {
    c(Mean = mean(x, na.rm = TRUE), Variance = var(x, na.rm = TRUE))
  })))

  ## plot the mean-variance relationship
  ggplot(plotData, aes(Mean, Variance)) +
    geom_point() +
    labs(title = title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}


##----------------------------------------------------------------------------------------
#'
#' Log-based transformation
#'
#' @description
#' Apply a logarithmic transformation to the data to stabilize the variance.
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param names A vector of strings (default = c("gender", "treatment", "replicate"))
#' specifying the names of the attribute columns.
#'
#' @param logFold An integer specifying the base for the log transformation.
#'
#' @details
#' The function executes the following:
#' \enumerate{
#' \item Plots the mean-variance relationship using \code{meanVariancePlot()}.
#' \item Log-transforms the data, using the specified base.
#' \item Plots the mean-variance relationship again for comparison.
#' }
#'
#' @returns The transformed data.
#'
#' @autoglobal
#'
#' @export

transform <- function(dataSet,
                      names = c("gender", "treatment", "replicate"),
                      logFold = 2) {

  ## organize the data for transformation
  dataPoints <- dataSet %>%
    select(-names)

  ## calculate and plot a mean-variance plot
  plotPre <- meanVarPlot(dataPoints, title = "Pre-Transformation")
  print(plotPre)

  ## take the log of the numerical data
  transDataPoints <- log(dataPoints, logFold)

  ## calculate and plot a mean-variance plot
  plotPost <- meanVarPlot(transDataPoints,  title = "Post-Transformation")
  print(plotPost)

  ## recombine the labels and transformed data into a single data frame
  transDataSet <- cbind(dataSet[,names], transDataPoints)

  ## return the transformed data
  return(transDataSet)
}
