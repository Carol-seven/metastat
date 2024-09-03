#'
#' Plotting a graph of normalized data
#'
#' @description
#' Plot the boxplot
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @import ggplot2
#'
#' @returns An object of class \code{plot}.
#'
#' @autoglobal
#'
#' @noRd

normPlot <- function(dataSet) {

  attrnames <- attributes(dataSet)$attrnames

  plotData <- dataSet %>%
    pivot_longer(-all_of(attrnames))
    # unite("attribute", attrnames[attrnames != "replicate"], sep = "-", remove = FALSE)

  ggplot(plotData, aes(x = merged_condition, y = value,
                       fill = if (length(unique(replicate)) == 1) attribute else replicate)) +
    geom_boxplot(varwidth = TRUE) +
    guides(fill = guide_legend(
      title = ifelse(length(unique(plotData$replicate)) == 1, "Condition", "Replicate"))) +
    labs(title = "Normalization Boxplot") +
    xlab("Condition") +
    ylab("Signal value") +
    scale_fill_brewer(palette = "RdYlBu") +
    theme_bw() +
    theme(legend.position = ifelse(length(unique(plotData$replicate)) == 1, "none", "bottom"),
          plot.title = element_text(hjust = .5))
}
