
#' Quantile normalization.
#'
#' @param FeatureTable Feature intensity table with samples in column and features in row (default).
#' @param IntThreshold Feature intensity threshold. Feature is detected when its intensity larger than this value.
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FLASE} if samples are in row.
#' @param output \code{TRUE} will output the result table in current working directory
#' @param OutputNormFactors \code{TRUE} will show the normalization factors after normalization
#' @param RunEvaluation \code{TRUE} will evaluate the normalization results using intragroup variation.
#'
#' @return
#' This function will return a list that contains two items: the normalized feature table,
#' and a set of normalization factors.
#' @export
#'
#' @examples
#' QuantileNormedTable = QuantileNorm(TestingData)


QuantileNorm = function(FeatureTable, IntThreshold=0, SampleInCol=TRUE, output=FALSE,
                   OutputNormFactors=FALSE, RunEvaluation=TRUE){
  message("Normalization is running...")

  # Transpose FeatureTable if samples are in row
  if (!SampleInCol) {
    FeatureTable = t(FeatureTable)
  }

  # Find names of sample groups
  group_seq = tolower(as.character(FeatureTable[1,-1]))
  temp = !(group_seq=="featurequality" | group_seq=="model" | group_seq=="sqc_points")
  group_seq = group_seq[temp]
  group_unique = unique(group_seq[-1])

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

  # Find sample files
  temp = !(grepl("SQC", group_seq,ignore.case = T) | group_seq=="blank" | group_seq=="rt" | group_seq=="qc")
  sample_table = data.matrix(IntTable[, temp])

  group_seq1 = tolower(as.character(FeatureTable[1,-1]))
  temp = !(grepl("SQC", group_seq1,ignore.case = T) | group_seq1=="blank" | group_seq1=="rt" | group_seq1=="qc" |
             group_seq1=="featurequality" | group_seq1=="model" | group_seq1=="sqc_points")
  FeatureTable_index = (2:ncol(FeatureTable))[c(temp)]
  group_vector = as.character(FeatureTable[1,FeatureTable_index])


  # Normalization Main
  # Data matrix for normalized data

  d = as.data.frame(preprocessCore::normalize.quantiles(as.matrix(sample_table)))

  for (i in 1:ncol(sample_table)) {
    FeatureTable[-1,FeatureTable_index[i]] = d[,i]
  }

  quality_exist = !is.null(FeatureTable$Quality)
  if (RunEvaluation) {
    pRMAD_each1 = c()
    pRMAD_each2 = c()
    for (i in 1:(nrow(FeatureTable)-1)) {
      if (quality_exist) {
        if(FeatureTable$Quality[i+1] == "low"){
          pRMAD_each1[i] = NaN
          pRMAD_each2[i] = NaN
          next
        }
      }

      d1 = as.numeric(sample_table[i,])
      d2 = as.numeric(FeatureTable[i+1, FeatureTable_index])
      pRMAD_each1[i] = pooled_rMAD(d1, group_vector)
      pRMAD_each2[i] = pooled_rMAD(d2, group_vector)
    }
    pRMAD1 = round(median(pRMAD_each1[!is.nan(pRMAD_each1)]), digits = 4)
    pRMAD2 = round(median(pRMAD_each2[!is.nan(pRMAD_each2)]), digits = 4)
    message(paste0("Median of PRMAD changed from ",
                   pRMAD1, " to ", pRMAD2, " after normalization."))
  }

  if (output) {
    write.csv(FeatureTable, "normalized data table.csv", row.names = F)
  }

  if (OutputNormFactors) {
    message("Normalization factors:")
    cat("Normalization factor is not available for this method.")
  }
  results = list(FeatureTable, NA)
  if (RunEvaluation) {
    results = list(FeatureTable, NA, pRMAD_each1, pRMAD_each2)
  }
  return(results)
  message("Normalization is done.")
}

