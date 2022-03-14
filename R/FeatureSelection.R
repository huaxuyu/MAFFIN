
#' High-quality feature selection
#' @description
#' Select high-quality features for quantitative analysis.
#'
#' @param FeatureTable Feature intensity table with samples in column and features in row (default).
#' @param BlankFilter High-quality when mean(sample intensities) > mean(blank intensities) * \code{BlankFilter}
#' @param RtRange Range of the defined retention time window, in minute.
#' @param QCRSD Relative standard deviation threshold for QC samples.
#' @param SQCcor Pearson's correlation threshold for serial QC samples (recommend: 0.8-0.9).
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FLASE} if samples are in row.
#' @param output \code{TRUE} will output the result table in current working directory
#' @param IntThreshold Feature intensity threshold. Feature is detected when its intensity larger than this value.
#'
#' @details  The first row should be sample names, and the second row should be group names.\cr
#' For group names, please use: \cr
#' "RT" for retention time column; \cr
#' "QC" for quality control samples between real samples (normal QC samples); \cr
#' "blank" for blank
#'
#' samples; \cr
#' "SQC_###" for serial QC samples with a certain loading amount.
#' For example, SQC_1.0 means a serial QC sample with injection volume as 1.0 uL.
#' @export
#' @return
#' This function will return the original feature table with an extra column named "Quality" to indicate the feature quality.
#'
#' @examples
#' SelectionTable = FeatureSelection(TestingData)


FeatureSelection = function(FeatureTable, BlankFilter=2, RtRange=c(0,100),
                            QCRSD=0.25, SQCcor=0.9, IntThreshold=0,
                            SampleInCol=TRUE, output=FALSE){
  message("Selecting high-quality features...")

  # Transpose FeatureTable if samples are in row
  if (!SampleInCol) {
    FeatureTable = t(FeatureTable)
  }
  filter.blank = filter.SQC = filter.RT = TRUE

  # Find names of sample groups
  group_seq = tolower(as.character(FeatureTable[1,-1]))
  group_unique = unique(group_seq[-1])

  # Convert feature intensities to numeric values
  # Remove the first row and column for downstream processing
  IntTable = FeatureTable[-1,-1]

  # Test if all cells in IntTable are numeric
  IntTable = tryCatch(sapply(IntTable, as.numeric),warning=function(w) w)
  if(is(IntTable,"warning")){
    print("Non-numeric value is found in feature intensities. Return NA.")
    return(NA)
  }

  # Generate a data frame to store the filtering results. 1 for pass, 0 for fail.
  # column 1: feature identifier
  # column 2: method blank filter result
  # column 3: retention time filter result
  # column 4: RSD in QC sample filter result
  # column 5: QC linearity filter result
  derep_table = data.frame(matrix(data = 1, nrow = nrow(IntTable), ncol = 5))
  colnames(derep_table) = c("Identifier", "MB", "RT", "QC_RSD", "Pearson_cor")
  derep_table[,1] = FeatureTable[-1,1]

  # Generate a table to store intensities from all blank samples
  blank_table = c()
  temp = group_seq == "blank"
  if(sum(temp) == 0){
    filter.blank = FALSE
    message("Blank data are not detected.")
  } else {
    blank_table = data.matrix(IntTable[, temp])
    if (length(sum(temp)) > 1) {
      blank_table = apply(blank_table, 1, mean)
    }
  }

  SQC_index = grep("SQC", group_unique, ignore.case = T)
  if(length(SQC_index) == 0){
    filter.SQC = FALSE
    message("Serial QC data are not detected.")
  } else if(length(SQC_index) < 5){
    filter.SQC = FALSE
    message("Serial QC data points are too less for evaluation (5 is required).")
  } else{
    SQC_table = data.frame(matrix(nrow = nrow(IntTable), ncol = length(SQC_index)))
    colnames(SQC_table) = group_unique[SQC_index]
    for (i in 1:length(SQC_index)) {
      SQC_table[,i] = apply(data.matrix(IntTable[,group_unique[SQC_index[i]] == group_seq]), 1, mean)
    }
    QC_conc = c()
    SQC.list = stringr::str_split(colnames(SQC_table), pattern = "_")
    for (i in 1:ncol(SQC_table)) {
      QC_conc[i] = as.numeric(SQC.list[[i]][2])
    }
  }

  temp = group_seq=="qc" | grepl("SQC", group_seq,ignore.case = T) | group_seq=="blank" | group_seq=="rt"
  sample_table = data.matrix(IntTable[, !temp])

  group_vector = group_seq[!temp]

  RT_v = data.matrix(IntTable[, group_seq=="rt"])
  if (length(RT_v) == 0) {
    filter.RT = FALSE
    message("Retention time data are not detected.")
  }

  temp = group_seq=="qc"
  if(sum(temp) != 0){
    MQC_table = data.matrix(IntTable[, temp])
    if(filter.blank){
      for (i in 1:nrow(IntTable)) {
        if(mean(sample_table[i,]) < blank_table[i]*BlankFilter) { derep_table[i,2] = 0 }
      }
    }
    if(sum(temp) >= 3){
      for (i in 1:nrow(IntTable)) {
        if(mean(MQC_table[i,] == 0)){
          derep_table[i,4] = 0
        } else {
          if(sd(MQC_table[i,])/mean(MQC_table[i,]) > QCRSD){
            derep_table[i,4] = 0}
        }
      }
    } else {
      message("QC data are not enough to calculate RSD (3 is required).")
    }
  } else {
    message("QC data (normal QC) are not detected.")
  }


  # Remove the features out of the defined retention time range
  for (i in 1:nrow(IntTable)) {
    # Check retention time
    if (filter.RT) {
      if(RT_v[i] < RtRange[1] | RT_v[i] > RtRange[2]) { derep_table[i,3] = 0 }
    }
  }

  all_filter = rep("low", nrow(IntTable))

  # Remove the features with low correlation in SQC samples
  if (filter.SQC) {
    for (i in 1:nrow(IntTable)) {
      # Check serial diluted QC linearity
      QC_int = as.numeric(SQC_table[i,])
      valid_int = which(QC_int > IntThreshold)
      QC_int = QC_int[valid_int]
      selected_QC_conc = QC_conc[valid_int]
      selected_int_points = length(selected_QC_conc)

      if(selected_int_points <= 3){
        derep_table[i,5] = 0
      } else if(cor(QC_int, selected_QC_conc) < SQCcor){
        derep_table[i,5] = 0
      }
      if (sum(derep_table[i,-1]) == 4) {
        all_filter[i] = "high"
      }
    }
  }

  FeatureTable$Quality = ""
  FeatureTable$Quality[1] = "FeatureQuality"
  FeatureTable$Quality[-1] = all_filter

  if (output) {
    write.csv(FeatureTable, "Feature_selection.csv", row.names = FALSE)
  }

  message("High-quality feature selection is done.")
  message(paste0(sum(all_filter=="high"), " features are selected from ", nrow(IntTable), "."))
  return(FeatureTable)
}
