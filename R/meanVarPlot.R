#'
#' Plotting a graph of mean versus variance
#'
#' @description
#' Take a set of metabolomics data organized by column, calculate the mean and variance of
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
