#'
#' Centering
#'
#' @description
#' Apply centering to the data.
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param names A vector of strings (default = c("gender", "treatment", "replicate"))
#' specifying the names of the attribute columns.
#'
#' @details
#' The function executes the following:
#' \enumerate{
#' \item Plots the mean-variance relationship.
#' \item Centers the data.
#' \item Plots the mean-variance relationship again for comparison.
#' }
#'
#' @returns The centered data.
#'
#' @autoglobal
#'
#' @export

center <- function(dataSet,
                   names = c("gender", "treatment", "replicate")) {

  ## organize the data for centering
  dataPoints <- dataSet %>%
    select(-any_of(names))

  ## calculate and plot a mean-variance plot
  plotPre <- meanVarPlot(dataPoints, title = "Pre-Centering")
  print(plotPre)

  ## centering
  centerData <- colMeans(dataPoints, na.rm=TRUE)
  centeredDataPoints <- sweep(dataPoints, 2L, centerData, check.margin = FALSE)

  ## calculate and plot a mean-variance plot
  plotPost <- meanVarPlot(centeredDataPoints,  title = "Post-Centering")
  print(plotPost)

  ## recombine the labels and centered data into a single data frame
  centeredDataSet <- cbind(dataSet[,names], centeredDataPoints)

  ## return the centered data
  return(centeredDataSet)
}
