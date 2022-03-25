
#' High-quality feature selection
#' @description
#' Select high-quality features for quantitative analysis.
#'
#' @param FeatureTable Data frame with features in row and samples in column (default).
#' @param BlankFilter A numeric value. High-quality when mean(sample intensities) > mean(blank intensities) * \code{BlankFilter}
#' @param RtRange A numeric vector indicating the range of the defined retention time window, in minute.
#' @param QCRSD A numeric value indicating the relative standard deviation threshold for QC samples.
#' @param SQCcor A numeric value indicating the Pearson's correlation threshold for serial QC samples (recommend: 0.8-0.9).
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FALSE} if samples are in row.
#' @param output \code{TRUE} will output the result table in current working directory
#' @param IntThreshold A numeric value indicating the feature intensity threshold. Feature is detected when its intensity larger than this value.
#'
#' @details
#' \code{FeatureTable} contains measured signal intensities of metabolic features,
#' with features in row and samples in column (default). The column names should
#' be sample names, and the first row should be sample group names (e.g. control, case).\cr
#' The first column should be unique feature identifiers.
#' For group names, please use: \cr
#' "RT" for retention time column; \cr
#' "QC" for quality control samples between real samples (normal QC samples); \cr
#' "blank" for blank samples; \cr
#' "SQC_###" for serial QC samples with a certain loading amount.
#' For example, SQC_1.0 means a serial QC sample with injection volume of 1.0 uL. \cr
#' An example of \code{FeatureTable} is provided as \code{TestingData} in this package.
#' @export
#' @return
#' This function will return the original data frame with an extra column
#' named "Quality" to indicate the feature quality.
#' @references Yu, Huaxu, and Tao Huan. "MAFFIN: Metabolomics Sample Normalization
#' Using Maximal Density Fold Change with High-Quality Metabolic Features and Corrected
#' Signal Intensities." \emph{bioRxiv} (2021).
#' @examples
#' selectedTable = FeatureSelection(TestingData)


FeatureSelection = function(FeatureTable, BlankFilter=2, RtRange=c(0,100),
                            QCRSD=0.25, SQCcor=0.9, IntThreshold=0,
                            SampleInCol=TRUE, output=FALSE){
  message("Selecting high-quality features...")

  # Transpose FeatureTable if samples are in row
  if (!SampleInCol) {
    FeatureTable = t(FeatureTable)
  }
  filter.blank = filter.RT = filter.QC = filter.SQC = TRUE

  # Find names of sample groups
  group_seq = tolower(as.character(FeatureTable[1,-1]))
  group_unique = unique(group_seq)

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
  derep_table = data.frame(matrix(data = TRUE, nrow = nrow(IntTable), ncol = 5))
  colnames(derep_table) = c("Identifier", "MB", "RT", "QC_RSD", "Pearson_cor")
  derep_table[,1] = FeatureTable[-1,1]

  # Calculate the average intensity from all blank samples
  blank_v = c()
  temp = group_seq=="blank"
  if(sum(temp) == 0){
    filter.blank = FALSE
    message("Blank data are not detected.")
  } else {
    blank_v = rowMeans(data.matrix(IntTable[, temp]))
  }

  # Find the column of retention times
  RT_v = c()
  temp = group_seq=="rt"
  if (sum(temp) == 0) {
    filter.RT = FALSE
    message("Retention time data are not detected.")
  } else {
    RT_v = data.matrix(IntTable[, temp])
  }

  # Calculate the relative standard deviation from normal QC samples
  QC_RSD_v = c()
  temp = group_seq=="qc"
  if(sum(temp) == 0){
    filter.QC = FALSE
    message("QC data are not detected.")
  } else if (sum(temp) < 3) {
    filter.QC = FALSE
    message("QC data are not enough to calculate RSD (3 is required).")
  } else {
    for (i in 1:nrow(IntTable)) {
      QC_RSD_v[i] = sd(IntTable[i, temp]) / mean(IntTable[i, temp])
      if (is.na(QC_RSD_v[i])) {
        QC_RSD_v[i] = 0
      }
    }
  }

  # Find intensities for serial QC samples
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
    QC_conc = c()
    SQC.list = stringr::str_split(colnames(SQC_table), pattern = "_")
    for (i in 1:ncol(SQC_table)) {
      QC_conc[i] = as.numeric(SQC.list[[i]][2])
    }
    for (i in 1:length(SQC_index)) {
      SQC_table[,i] = rowMeans(data.frame(IntTable[,group_seq == group_unique[SQC_index[i]]]))

    }
  }

  # Calculate the average intensity from real samples
  temp = group_seq=="qc" | grepl("SQC", group_seq,ignore.case = T) | group_seq=="blank" | group_seq=="rt"
  sample_v = rowMeans(data.matrix(IntTable[, !temp]))

  # Remove the features in blank
  if (filter.blank) {
    derep_table[,2] = sample_v > blank_v*BlankFilter
  }

  # Remove the features out of the defined retention time range
  if (filter.RT) {
    derep_table[,3] = RT_v > RtRange[1] & RT_v < RtRange[2]
  }

  # Remove the features with low RSD in QC samples
  if (filter.QC) {
    derep_table[,4] = QC_RSD_v < QCRSD
  }

  # Remove the features with low correlation in SQC samples
  if (filter.SQC) {
    for (i in 1:nrow(IntTable)) {
      # Check serial diluted QC linearity
      QC_int = as.numeric(SQC_table[i,])
      valid_int = which(QC_int > IntThreshold)
      QC_int = QC_int[valid_int]
      selected_QC_conc = QC_conc[valid_int]
      selected_int_points = length(selected_QC_conc)

      if(selected_int_points > 3){
        if (cor(QC_int, selected_QC_conc) > SQCcor) {
          derep_table[i,5] = TRUE
        } else {
          derep_table[i,5] = FALSE
        }
      } else {
        derep_table[i,5] = FALSE
      }
    }
  }

  FeatureTable$Quality = ""
  FeatureTable$Quality[1] = "FeatureQuality"
  reasons = c("Blank", "RT", "LowRSD", "LowCor")

  for (i in 1:nrow(derep_table)) {
    if (sum(derep_table[i,-1]) == 4) {
      FeatureTable$Quality[i+1] = "high-quality"
    } else {
      temp = paste(reasons[!derep_table[i,-1]], collapse=", ")
      FeatureTable$Quality[i+1] = paste("low-quality | ", temp)
    }
  }

  if (output) {
    write.csv(FeatureTable, "Feature_selection.csv", row.names = FALSE)
  }

  message("High-quality feature selection is done.")
  message(paste0(sum(FeatureTable$Quality[-1]=="high-quality"), " features are selected from ", nrow(IntTable), "."))
  return(FeatureTable)
}
