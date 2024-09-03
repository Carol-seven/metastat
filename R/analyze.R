#'
#' Student's t-test
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param condCol A string (default = "merged_condition") specifying the column name in \code{dataSet} that contains
#' the conditions to compared.
#'
#' @param cond A string specifying which two conditions to compare. The order is
#' important, as the second condition serves as the reference for comparison. When there
#' are two conditions in \code{condCol} of \code{dataSet} and this argument is not
#' specified, the \code{cond} will automatically be selected by sorting the unique
#' values alphabetically and in ascending order.
#'
#' @import dplyr
#' @importFrom stats t.test
#'
#' @returns A data frame containing the differences in means and p-values for each
#' compound between the two conditions.
#'
#' @export

t_test <- function(dataSet, condCol = "merged_condition", cond) {

  ## check for exactly two conditions
  if (missing(cond)) {
    cond <- sort(unique(dataSet[[condCol]]))
    if (length(cond) != 2) {
      stop("Please provide exactly two conditions for comparison.")
    }
  } else if (length(cond) != 2) {
    stop("This analysis can only be performed on two conditions at a time.
         Please select exactly two conditions to compare.
         Alternatively, please consider using ANOVA.")
  }

  ## filter data set by the conditions
  filteredData <- dataSet %>%
    filter(.data[[condCol]] %in% cond)

  ## index of the two conditions
  indexA <- which(filteredData[[condCol]] == cond[1])
  indexB <- which(filteredData[[condCol]] == cond[2])

  ## the difference in means and the p-value of t-test
  result <- as.data.frame(apply(
    select(filteredData, -any_of(attributes(dataSet)$attrnames)), 2,
    function(x) {
      tryCatch(
        c("Difference" = mean(x[indexA])-mean(x[indexB]),
          "P-value" = t.test(x[indexA], x[indexB])$p.value),

        ## if an error is thrown, return the fold change and set the p-value to 'NA'.
        error = function(e) {
          message("Data are essentially constant.")
          c("Difference" = mean(x[indexA])-mean(x[indexB]), "P-value" = NA)
        }
      )
    }
  ))

  return(result)

}


##----------------------------------------------------------------------------------------
#'
#' Moderated t-test
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param condCol A string (default = "merged_condition") specifying the column name in \code{dataSet} that contains
#' the conditions to compared.
#'
#' @param cond A string specifying which two conditions to compare. The order is
#' important, as the second condition serves as the reference for comparison. When there
#' are two conditions in \code{condCol} of \code{dataSet} and this argument is not
#' specified, the \code{cond} will automatically be selected by sorting the unique
#' values alphabetically and in ascending order.
#'
#' @import dplyr
#' @importFrom stats model.matrix
#'
#' @returns A data frame containing the differences in means and p-values for each
#' compound between the two conditions.
#'
#' @export

mod_t_test <- function(dataSet, condCol = "merged_condition", cond) {

  ## check for exactly two conditions
  if (missing(cond)) {
    cond <- sort(unique(dataSet[[condCol]]))
    if (length(cond) != 2) {
      stop("Please provide exactly two conditions for comparison.")
    }
  } else if (length(cond) != 2) {
    stop("This analysis can only be performed on two conditions at a time.
         Please select exactly two conditions to compare.
         Alternatively, please consider using ANOVA.")
  }

  ## filter data set by the conditions
  filteredData <- dataSet %>%
    filter(.data[[condCol]] %in% cond)

  cond <- factor(filteredData[[condCol]], levels = cond, labels = LETTERS[1:2])
  design <- model.matrix(~ 0 + cond)

  ## fit linear model for each protein
  fit1 <- limma::lmFit(
    t(select(filteredData, -any_of(attributes(dataSet)$attrnames))), design)

  ## construct the contrast matrix
  cont.matrix <- limma::makeContrasts("condA - condB", levels = design)

  ## compute contrasts from linear model 'fit1'
  fit2 <- limma::contrasts.fit(fit1, cont.matrix)

  ## empirical Bayes statistics
  fit3 <- limma::eBayes(fit2)

  result <- as.data.frame(t(cbind(fit3$coefficients, fit3$p.value)))
  rownames(result) <- c("Difference", "P-value")

  return(result)

}


##----------------------------------------------------------------------------------------
#'
#' Wilcoxon signed-rank test
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param condCol A string (default = "merged_condition") specifying the column name in \code{dataSet} that contains
#' the conditions to compared.
#'
#' @param cond A string specifying which two conditions to compare. The order is
#' important, as the second condition serves as the reference for comparison. When there
#' are two conditions in \code{condCol} of \code{dataSet} and this argument is not
#' specified, the \code{cond} will automatically be selected by sorting the unique
#' values alphabetically and in ascending order.
#'
#' @import dplyr
#' @importFrom stats wilcox.test
#'
#' @returns A data frame containing the differences in means and p-values for each
#' compound between the two conditions.
#'
#' @export

wilcox_test <- function(dataSet, condCol = "merged_condition", cond) {

  ## check for exactly two conditions
  if (missing(cond)) {
    cond <- sort(unique(dataSet[[condCol]]))
    if (length(cond) != 2) {
      stop("Please provide exactly two conditions for comparison.")
    }
  } else if (length(cond) != 2) {
    stop("This analysis can only be performed on two conditions at a time.
         Please select exactly two conditions to compare.
         Alternatively, please consider using ANOVA.")
  }

  ## filter data set by the conditions
  filteredData <- dataSet %>%
    filter(.data[[condCol]] %in% cond)

  ## index of the two conditions
  indexA <- which(filteredData[[condCol]] == cond[1])
  indexB <- which(filteredData[[condCol]] == cond[2])

  ## the difference in means and the p-value of t-test
  result <- as.data.frame(apply(
    select(filteredData, -any_of(attributes(dataSet)$attrnames)), 2,
    function(x) {
      tryCatch(
        c("Difference" = mean(x[indexA])-mean(x[indexB]),
          "P-value" = wilcox.test(x[indexA], x[indexB])$p.value),

        ## if an error is thrown, return the fold change and set the p-value to 'NA'.
        error = function(e) {
          message("Data are essentially constant.")
          c("Difference" = mean(x[indexA])-mean(x[indexB]), "P-value" = NA)
        }
      )
    }
  ))

  return(result)

}

