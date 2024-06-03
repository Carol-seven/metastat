#'
#' Centering and/or scaling
#'
#' @description
#' Apply centering and/or scaling to the data.
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param names A vector of strings (default = c("gender", "treatment", "replicate"))
#' specifying the names of the attribute columns.
#'
#' @param center A boolean (default = TRUE) specifying whether to apply centering.
#'
#' @param scale A string (default = "auto") specifying which method of scaling to apply:
#' \enumerate{
#' \item "none" No scaling will be applied.
#' \item "auto": Auto scaling.
#' \item "level": Level scaling.
#' \item "pareto": Pareto scaling.
#' \item "range": Range scaling.
#' \item "vast": Vast scaling.
#' }
#' @param mu A numeric vector specifying the centering reference. The default centering
#' reference is set to the mean value.
#'
#' @param sigma A numeric vector specifying the scaling reference based on the corresponding
#' definition. The default scaling reference is set to the standard deviation when
#' \code{center = TRUE} and the adjusted root mean square when \code{center = FALSE}.
#'
#' @details
#' The function executes the following:
#' \enumerate{
#' \item Plots the mean-variance relationship.
#' \item Centers and/or scales the data.
#' \item Plots the mean-variance relationship again for comparison.
#' }
#'
#' @returns The centered and/or scaled data.
#'
#' @autoglobal
#'
#' @export

scale <- function(dataSet,
                  names = c("gender", "treatment", "replicate"),
                  center = TRUE, scale = "auto",
                  mu = NULL, sigma = NULL) {

  ## organize the data for centering
  dataPoints <- dataSet %>%
    select(-any_of(names))

  ## calculate and plot a mean-variance plot
  plotPre <- meanVarPlot(dataPoints, title = "Pre-Pretreatment")
  print(plotPre)

  if (is.null(mu)) {
    mu <- colMeans(dataPoints, na.rm = TRUE)
  }

  ## centering
  if (center) {
    dataPoints <- sweep(dataPoints, 2L, mu, check.margin = FALSE)
  }

  if (is.null(sigma)) {
    sigma <- apply(dataPoints, 2, function(col) {
      col <- col[!is.na(col)]
      sqrt(sum(col^2) / max(1L, length(col)-1L))
    })
  }

  ## scaling
  if (scale == "none") {
    result <- dataPoints
  } else if (scale == "auto") {
    result <- sweep(dataPoints, 2L, sigma, `/`, check.margin = FALSE)
  } else if (scale == "level") {
    result <- sweep(dataPoints, 2L, mu, `/`, check.margin = FALSE)
  } else if (scale == "pareto") {
    result <- sweep(dataPoints, 2L, sqrt(sigma), `/`, check.margin = FALSE)
  } else if (scale == "range") {
    range_diff <- apply(dataPoints, 2, function(col) {
      res <- max(col, na.rm = TRUE) - min(col, na.rm = TRUE)
      res <- ifelse(res == 0, 1, res)
      return(res)
    })
    result <- sweep(dataPoints, 2L, range_diff, `/`, check.margin = FALSE)
  } else if (scale == "vast") {
    result <- sweep(dataPoints, 2L, mu/sigma^2, `*`, check.margin = FALSE)
  }

  ## calculate and plot a mean-variance plot
  plotPost <- meanVarPlot(result,  title = "Post-Pretreatment")
  print(plotPost)

  ## recombine the labels and centered and/or scaled data into a single data frame
  result <- cbind(dataSet[,names], result)

  ## return the centered and/or scaled data
  return(result)
}
