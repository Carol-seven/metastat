#'
#' Plotting a graph of normalized data
#'
#' @description
#' Plot the boxplot
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param names A vector of strings (default = c("gender", "treatment", "replicate"))
#' specifying the names of the attribute columns.
#'
#' @import ggplot2
#'
#' @returns An object of class \code{plot}.
#'
#' @autoglobal
#'
#' @noRd

normPlot <- function(dataSet,
                     names = c("gender", "treatment", "replicate")) {

  plotData <- dataSet %>%
    pivot_longer(-names) %>%
    unite("attribute", names[-length(names)], sep = "-", remove = FALSE)

  ggplot(plotData, aes(x = attribute, y = value,
                       fill = if (length(unique(replicate)) == 1) attribute else replicate)) +
    geom_boxplot(varwidth = TRUE) +
    guides(fill = guide_legend(
      title = ifelse(length(unique(plotData$replicate)) == 1, "Attribute", "Replicate"))) +
    labs(title = "Normalization Boxplot") +
    xlab("Attribute") +
    ylab("Signal value") +
    scale_fill_brewer(palette = "RdYlBu") +
    theme_bw() +
    theme(legend.position = ifelse(length(unique(plotData$replicate)) == 1, "none", "bottom"),
          plot.title = element_text(hjust = .5))
}
