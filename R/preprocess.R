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
#' @param names A vector of strings (default = c("gender", "treatment", "replicate"))
#' specifying the names of the attribute columns.
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
#' @import tidyr
#' @importFrom utils read.csv
#'
#' @returns A 2d dataframe.
#'
#' @autoglobal
#'
#' @export

preprocess <- function(fileName, dataSet = NULL,
                       names = c("gender", "treatment", "replicate")) {
  if (missing(fileName)) {
    if (is.null(dataSet)) {
      stop("Either 'fileName' or 'dataSet' must be provided.")
    }
  } else {
    ## read in the metabolomics quantitative csv file
    dataSet <- read.csv(fileName)
  }

  ## long dataSet table
  dataSet.l <- dataSet %>%
    pivot_longer(cols = -Compound) %>%
    mutate(value = ifelse(value == 0, NA, value))

  ## reformatted data
  result <- dataSet.l %>%
    pivot_wider(id_cols = name, names_from = Compound, values_from = value) %>%
    separate(name, into = names, sep = "_") %>%
    mutate(replicate = factor(replicate, levels = sort(as.numeric(unique(replicate))))) %>%
    as.data.frame()

  ## generate a histogram of the log2-transformed values for full raw data set
  value <- data.frame(value = log2(dataSet.l$value))
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
  print(summary(dataSet.l$value))
  cat("\n")

  ## print levels of attributes
  for (k in 1:length(names)) {
    cat(paste0("Levels of ", names[k], ":"),
        eval(parse(text = paste0("unique(result$", names[k], ")"))), "\n")
  }
  cat("\n")

  return(result)
}
