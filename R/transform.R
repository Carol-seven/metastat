#'
#' Transformation
#'
#' @description
#' Apply a transformation to the data.
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param names A vector of strings (default = c("gender", "treatment", "replicate"))
#' specifying the names of the attribute columns.
#'
#' @param method A string (default = "log") specifying the method to be used for the
#' transformation:
#' \enumerate{
#' \item "log": Logarithm transformation.
#' \item "root": Root transformation.
#' }
#'
#' @param logFold An integer (default = 2) specifying the base for the log transformation
#' when \code{method = "log"}.
#'
#' @param root An integer (default = 2) specifying the degree of the root for the root
#' transformation when \code{method = "root"}. For example, set it to 2 for square root or
#' 3 for cube root.
#'
#' @details
#' The function executes the following:
#' \enumerate{
#' \item Plots the mean-variance relationship.
#' \item Transforms the data.
#' \item Plots the mean-variance relationship again for comparison.
#' }
#'
#' @returns The transformed data.
#'
#' @export

transform <- function(dataSet,
                      names = c("gender", "treatment", "replicate"),
                      method = "log", logFold = 2, root = 2) {

  ## organize the data for transformation
  dataPoints <- select(dataSet, -any_of(names))

  ## calculate and plot a mean-variance plot
  plotPre <- meanVarPlot(dataPoints, title = "Pre-Transformation")
  print(plotPre)

  ## take the log of the numerical data
  if (method == "log") {
    transDataPoints <- log(dataPoints, logFold)
  } else if (method == "root") {
    transDataPoints <- dataPoints^(1/root)
  }

  ## calculate and plot a mean-variance plot
  plotPost <- meanVarPlot(transDataPoints,  title = "Post-Transformation")
  print(plotPost)

  ## recombine the labels and transformed data into a single data frame
  transDataSet <- cbind(dataSet[,names], transDataPoints)

  ## return the transformed data
  return(transDataSet)
}
