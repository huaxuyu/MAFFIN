
#' Maximal density fold change normalization
#'
#' @param FeatureTable Feature intensity table with samples in column and features in row (default).
#' @param IntThreshold Feature intensity threshold. Feature is detected when its intensity larger than this value.
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FLASE} if samples are in row.
#' @param output \code{TRUE} will output the result table in current working directory
#' @param OutputNormFactors \code{TRUE} will show the normalization factors after normalization
#'
#' @return
#' This function will return a list contains two items: the normalized feature table,
#' and a set of normalization factors.
#'
#' @export
#'
#' @examples
#' Please see GitHub for demo.


MDFCNorm = function(FeatureTable, IntThreshold=0, SampleInCol=TRUE, output=FALSE,
                    OutputNormFactors=TRUE){
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

  # Find the reference file with the most detected features
  f_number = c()
  for (i in 1:ncol(sample_table)) {
    f_number[i] = sum(sample_table[,i] > IntThreshold)
  }

  # Normalization Main
  # Data matrix for normalized data

  ref_file = as.numeric(sample_table[, which.max(f_number)])

  # Only use high-quality features if labeled
  quality_exist = !is.null(FeatureTable$Quality)
  if (quality_exist) {
    hq_table = sample_table[FeatureTable$Quality[-1] == "high", ]
    ref_file = ref_file[FeatureTable$Quality[-1] == "high"]
  }

  best_bw = bw_opt(hq_table, group_vector)

  f = c()
  for (i in 1:ncol(sample_table)) {
    if(i == which.max(f_number)){
      f[i] = 1
      next
    }
    # Use the commonly detected features for normalization
    fc_seq = as.numeric(hq_table[,i]) / ref_file
    fc_seq = fc_seq[as.numeric(hq_table[,i]) > IntThreshold & ref_file > IntThreshold]
    d = density(log2(fc_seq), bw = best_bw)
    norm_factor = 2^d$x[which.max(d$y)]
    f[i] = norm_factor
    FeatureTable[-1,FeatureTable_index[i]] = as.numeric(FeatureTable[-1,FeatureTable_index[i]]) / norm_factor
  }

  if (output) {
    write.csv(FeatureTable, "normalized data table.csv", row.names = F)
  }

  if (OutputNormFactors) {
    message("Normalization factors:")
    print(f)
  }
  results = list(FeatureTable, f)
  return(results)
  message("Normalization is done.")
}
