
#' Calculate pooled RMAD for normalization evaluation.
#'
#' #' @description
#' Calculate pooled relative median absolute deviation for each metabolic feature.
#'
#' @param FeatureTable Feature intensity table with features in row and samples in column (default).
#' @param GroupNames A character vector indicating the names of each group.
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FLASE} if samples are in row.
#' @param output \code{TRUE} will output the result table in current working directory
#'
#' @details
#' \code{FeatureTable} contains measured signal intensities of metabolic features,
#' with features in row and samples in column (default). The column names should
#' be sample names, and the first row should be sample group names (e.g. control, case).\cr
#' The first column should be unique feature identifiers.
#' An example of \code{FeatureTable} is provided as \code{TestingData} in this package.
#'
#' @return
#' This function will return a vector that contains the calculated PRMADs for all features.
#'
#' @export
#'
#' @references Yu, Huaxu, and Tao Huan. "MAFFIN: Metabolomics Sample Normalization
#' Using Maximal Density Fold Change with High-Quality Metabolic Features and Corrected
#' Signal Intensities." \emph{bioRxiv} (2021).
#'
#' @examples
#' prmad = EvaPRMAD(TestingData, GroupNames=c("HY", "SX", "SW", "YC"))

EvaPRMAD = function(FeatureTable, GroupNames, SampleInCol=TRUE, output=FALSE){

  # Transpose FeatureTable if samples are in row
  if (!SampleInCol) {
    FeatureTable = t(FeatureTable)
  }

  # Find names of sample groups
  group_seq = tolower(as.character(FeatureTable[1,-1]))
  GroupNames = tolower(GroupNames)
  temp = rep(TRUE, length(group_seq))
  for (i in 1:length(GroupNames)) {
    temp = temp | group_seq==GroupNames[i]
  }
  group_vector = group_seq[temp]
  # Convert feature intensities to numeric values
  # Remove the first row and column for downstream processing
  IntTable = FeatureTable[-1,-1]
  IntTable = IntTable[,temp]

  # Test if all cells in IntTable are numeric
  IntTable = tryCatch(sapply(IntTable, as.numeric),warning=function(w) w)
  if(is(IntTable,"warning")){
    print("Non-numeric value is found in feature intensities. Return NA.")
    return(NA)
  }

  pRMAD_each = c()
  for (i in 1:(nrow(FeatureTable)-1)) {
    d = as.numeric(IntTable[i,])
    pRMAD_each[i] = pooled_rMAD(d, group_vector)
  }
  return(pRMAD_each)
}


#' Calculate pooled RSD for normalization evaluation.
#'
#' @description
#' Calculate pooled relative standard deviation for each metabolic feature.
#'
#' @param FeatureTable Feature intensity table with features in row and samples in column (default).
#' @param GroupNames A character vector indicating the names of each group.
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FALSE} if samples are in row.
#' @param output \code{TRUE} will output the result table in current working directory
#'
#' @details
#' \code{FeatureTable} contains measured signal intensities of metabolic features,
#' with features in row and samples in column (default). The column names should
#' be sample names, and the first row should be sample group names (e.g. control, case).\cr
#' The first column should be unique feature identifiers.
#' An example of \code{FeatureTable} is provided as \code{TestingData} in this package.
#'
#' @return
#' This function will return a vector that contains the calculated PRSDs for all features.
#'
#' @export
#'
#' @references Yu, Huaxu, and Tao Huan. "MAFFIN: Metabolomics Sample Normalization
#' Using Maximal Density Fold Change with High-Quality Metabolic Features and Corrected
#' Signal Intensities." \emph{bioRxiv} (2021).
#'
#' @examples
#' prsd = EvaPRMAD(TestingData, GroupNames=c("HY", "SX", "SW", "YC"))


EvaPRSD = function(FeatureTable, GroupNames, SampleInCol=TRUE, output=FALSE){

  # Transpose FeatureTable if samples are in row
  if (!SampleInCol) {
    FeatureTable = t(FeatureTable)
  }

  # Find names of sample groups
  group_seq = tolower(as.character(FeatureTable[1,-1]))
  temp = rep(TRUE, length(group_seq))
  for (i in 1:length(GroupNames)) {
    temp = temp | group_seq==GroupNames[i]
  }
  group_vector = group_seq[temp]
  # Convert feature intensities to numeric values
  # Remove the first row and column for downstream processing
  IntTable = FeatureTable[-1,-1]
  IntTable = IntTable[,temp]

  # Test if all cells in IntTable are numeric
  IntTable = tryCatch(sapply(IntTable, as.numeric),warning=function(w) w)
  if(is(IntTable,"warning")){
    print("Non-numeric value is found in feature intensities. Return NA.")
    return(NA)
  }

  pRSD_each = c()
  for (i in 1:(nrow(FeatureTable)-1)) {
    d = as.numeric(IntTable[i,])
    pRSD_each[i] = pooled_rsd(d, group_vector)
  }
  return(pRSD_each)
}
