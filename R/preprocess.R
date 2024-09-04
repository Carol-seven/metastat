#'
#' Loading and reformatting of metabolomics data
#'
#' @description
#' Read a data file, select columns necessary for analysis, and return the reformatted data.
#'
#' @param fileName The name of the .csv file containing metabolomics data (including the
#' path to the file, if needed).
#'
#' @param dataSet The raw data set, if already loaded in R.
#'
#' @param attrnames A vector of strings (default = c("gender", "treatment", "replicate"))
#' specifying the names of the attribute columns.
#'
#' @param zeroNA A boolean (default = TRUE) specifying whether 0's should be converted to
#' NA's.
#'
#' @details
#' The function executes the following:
#' \enumerate{
#' \item Reads the file.
#' \item Provides summary statistics and a histogram of all values reported in the data set.
#' \item Re-formats the data to present individual compounds as columns.
#' \item Stores the data as a \code{data.frame} and prints the levels of attributes to the
#' user.
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @import tibble
#' @import tidyr
#' @importFrom utils read.csv
#'
#' @returns A 2d dataframe.
#'
#' @autoglobal
#'
#' @export

preprocess <- function(fileName, dataSet = NULL,
                       attrnames = c("gender", "treatment", "replicate"),
                       zeroNA = TRUE) {
  if (missing(fileName)) {
    if (is.null(dataSet)) {
      stop("Either 'fileName' or 'dataSet' must be provided.")
    }
  } else {
    ## read in the metabolomics quantitative csv file
    dataSet <- read.csv(fileName)
  }

  dataPoints <- t(select(dataSet, -Compound))

  if (zeroNA) {
    dataPoints[dataPoints == 0] <- NA
  }

  colnames(dataPoints) <- dataSet$Compound

  ## reformatted data
  result <- as.data.frame(dataPoints) %>%
    rownames_to_column("name") %>%
    separate(name, into = attrnames, sep = "_") %>%
    unite("merged_condition", all_of(setdiff(attrnames, "replicate")), sep = "_", remove = FALSE) %>%
    mutate(replicate = factor(replicate, levels = sort(as.numeric(unique(replicate)))))

  ## generate a histogram of the log2-transformed values for full raw data set
  value <- data.frame(value = as.vector(log2(dataPoints)))
  plot <- ggplot(value) +
    geom_histogram(aes(x = value),
                   binwidth = 2,
                   # breaks = seq(floor(min(value, na.rm = TRUE)),
                   #              ceiling(max(value, na.rm = TRUE)), 1),
                   color = "black", fill = "gray") +
    # scale_x_continuous(breaks = seq(floor(min(value, na.rm = TRUE)),
    #                                 ceiling(max(value, na.rm = TRUE)), 2)) +
    labs(title = "Histogram of Full Raw Data Set",
         x = expression("log"[2]*"(Data)"), y = "Frequency") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  print(plot)

  ## print summary statistics for full raw data set
  cat("\nSummary of Full Data Signals (Raw):\n")
  print(summary(value$value))
  cat("\n")

  ## print levels of attributes
  for (k in 1:length(attrnames)) {
    cat(paste0("Levels of ", attrnames[k], ":"),
        eval(parse(text = paste0("unique(result$", attrnames[k], ")"))), "\n")
  }
  cat("\n")

  attr(result, "attrnames") <- c("merged_condition", attrnames)

  return(result)
}
