#'
#' Filtering compounds
#'
#' @description
#' Apply a series of filtering steps to the data set.
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param listName A character vector of compounds to select or remove.
#'
#' @param regexName A character vector specifying compounds for regular expression pattern
#' matching to select or remove.
#'
#' @param removeList A boolean (default = TRUE) specifying whether the list of compounds
#' should be removed or selected.
#' \itemize{
#' \item TRUE: Remove the list of compounds from the data.
#' \item FALSE: Remove all compounds not in the list from the data.
#' }
#'
#' @param saveRm A boolean (default = TRUE) specifying whether to save removed data to
#' current working directory. This option only works when \code{removeList = TRUE}.
#'
#' @details
#' If both \code{listName} and \code{regexName} are provided, the compound names selected
#' or removed will be the union of those specified in \code{listName} and those matching
#' the regex pattern in \code{regexName}.
#'
#' @import dplyr
#' @importFrom utils write.csv
#'
#' @return The filtered data.
#'
#' @export

filterOutIn <- function(dataSet,
                        listName = c(),
                        regexName = c(),
                        removeList = TRUE,
                        saveRm = TRUE) {

  attrnames <- attributes(dataSet)$attrnames

  ## relabel the data frame
  filteredData <- select(dataSet, -any_of(attrnames))

  ## only list filter if listName is present
  if (length(listName) != 0) {
    listIndex <- which(colnames(filteredData) %in% listName)
  } else {
    listIndex <- NULL
  }

  ## only regex filter if regexName is present
  if (length(regexName) != 0) {
    regexIndex <- grep(paste(regexName, collapse = "|"), colnames(filteredData))
  } else {
    regexIndex <- NULL
  }

  ## combine compound names from list and regex filters
  unionName <- colnames(filteredData)[sort(union(listIndex, regexIndex))]

  ## create a dataframe of the data of compounds
  unionData <- select(dataSet, any_of(c(attrnames, unionName)))

  ## if contaminants are being removed
  if (removeList == TRUE) {

    if (saveRm) {

      ## save removed data to current working directory
      write.csv(unionData, file = "filtered_out_data.csv", row.names = FALSE)
    }

    ## remove all of the contaminants if they are present
    filteredData <- select(dataSet, -any_of(unionName))

    ## if certain proteins are being selected
  } else if (removeList == FALSE) {

    ## select only compounds of interest
    filteredData <- unionData
  }

  attributes(filteredData)$attrnames <- attrnames

  ## return the filtered data
  return(filteredData)
}


##----------------------------------------------------------------------------------------
#'
#' Filtering NA's post-imputation
#'
#' @description
#' Remove compounds with NA values.
#'
#' @param dataSet A data frame containing the data signals.
#'
#' @param saveRm A boolean (default = TRUE) specifying whether to save removed data to
#' current working directory.
#'
#' @details
#' If compounds that do not meet the imputation requirement are removed, a .csv file is
#' created with the removed data.
#'
#' @import dplyr
#' @importFrom utils write.csv
#'
#' @return The filtered data.
#'
#' @export

filterNA <- function(dataSet, saveRm = TRUE) {

  attrnames <- attributes(dataSet)$attrnames

  if (saveRm) {

    ## create a dataframe of the removed data
    removedData <- bind_cols(select(dataSet, attrnames),
                             select_if(dataSet, ~any(is.na(.))))

    ## save removed data to current working directory
    write.csv(removedData, file = "filtered_NA_data.csv", row.names = FALSE)
  }

  ## remove all of the contaminants if they are present
  filteredData <- dataSet %>% select_if(~!any(is.na(.)))
  attributes(filteredData)$attrnames <- attrnames

  ## return the filtered data
  return(filteredData)
}
