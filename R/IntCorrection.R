
#' Intensity correction using serial QC samples
#'
#' @description
#' Correct MS signal intensities using serial QC samples
#'
#' @param FeatureTable Data frame with features in row and samples in column (default).
#' @param IntThreshold A numeric value indicating the feature intensity threshold. Feature is detected when its intensity larger than this value.
#' @param LR_QC_points Minimum serial QC data points for quadratic regression.
#' @param QR_QC_points Minimum serial QC data points for cubic regression.
#' @param SQCcor Pearson's correlation threshold for serial QC samples (recommend: 0.8-0.9).
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FALSE} if samples are in row.
#' @param output \code{TRUE} will output the result table in current working directory.
#'
#' @details
#' \code{FeatureTable} contains measured signal intensities of metabolic features,
#' with features in row and samples in column (default). The column names should
#' be sample names, and the first row should be sample group names (e.g. control, case).\cr
#' The first column should be unique feature identifiers.
#' For group names, please use: \cr
#' "QC" for quality control samples between real samples (normal QC samples); \cr
#' "SQC_###" for serial QC samples with a certain loading amount.
#' For example, SQC_1.0 means a serial QC sample with injection volume of 1.0 uL. \cr
#' Please note, group names of real biological samples cannot be "RT" and "blank". \cr
#' An example of \code{FeatureTable} is provided as \code{TestingData} in this package.
#'
#' @return This function will return the original feature table with corrected intensities.
#' @export
#'
#' @references Yu, Huaxu, and Tao Huan. "MAFFIN: Metabolomics Sample Normalization
#' Using Maximal Density Fold Change with High-Quality Metabolic Features and Corrected
#' Signal Intensities." \emph{bioRxiv} (2021).
#'
#' @examples
#' intCorrectedTable = IntCorrection(TestingData)

IntCorrection = function(FeatureTable, IntThreshold=0, LR_QC_points=5, QR_QC_points=7,
                         SQCcor=0.9, SampleInCol=TRUE, output=FALSE){

  message("MS intensity correction is running...")

  model_pool = c("Uncali.","Linear", "Quadratic", "Cubic")

  # Separate SQC table, sample table from FeatureTable
  # Transpose FeatureTable if samples are in row
  if (!SampleInCol) {
    FeatureTable = t(FeatureTable)
  }

  # Find names of sample groups
  group_seq = tolower(as.character(FeatureTable[1,-1]))
  temp = !(group_seq=="featurequality")
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

  # Check SQC availability and separate SQC table
  SQC_index = grep("SQC", group_unique, ignore.case = T)
  if(length(SQC_index) == 0){
    message("Serial QC data are not detected. Intensity correction failed.")
    return(NA)
  } else if(length(SQC_index) < 5){
    message("Serial QC data are not enough (5 is required). Intensity correction failed.")
    return(NA)
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

  quality_not_exist = is.null(FeatureTable$Quality)
  # If the feature quality is not determined yet
  if(quality_not_exist){
    FeatureTable$Quality = "high"
    FeatureTable$Quality[1] = "Quality"
    # Check the Pearson's correlation of all the features
    for (i in 1:nrow(IntTable)) {
      QC_int = as.numeric(SQC_table[i,])
      valid_int = which(QC_int > IntThreshold)
      QC_int = QC_int[valid_int]
      selected_QC_conc = QC_conc[valid_int]
      selected_int_points = length(selected_QC_conc)

      if (cor(QC_int, selected_QC_conc) < SQCcor | length(QC_int) < 5){
        FeatureTable$Quality[i+1] = "low"
      }
    }
  }

  # Calibrate sample and QC intensities
  temp = !(grepl("SQC", group_seq,ignore.case = T) | group_seq=="blank" | group_seq=="rt" |
             group_seq=="featurequality")
  sample_table = IntTable[,temp]
  temp_quality = FeatureTable$Quality[-1]

  model = rep(NA, nrow(sample_table))
  SQCpoint = rep(NA, nrow(sample_table))

  pb <- txtProgressBar(min = 0, max = nrow(sample_table), style = 3)

  for (i in 1:nrow(sample_table)) {
    if (temp_quality[i] == "high") {
      QC_int = as.numeric(SQC_table[i,])
      valid_int = which(QC_int > IntThreshold)
      QC_int = QC_int[valid_int]
      selected_QC_conc = QC_conc[valid_int]
      selected_int_points = length(selected_QC_conc)

      if(selected_int_points < LR_QC_points){
        best_model = cross_validation(as.numeric(QC_int),selected_QC_conc,1)
      } else if(selected_int_points>=LR_QC_points & selected_int_points<QR_QC_points){
        best_model = cross_validation(as.numeric(QC_int),selected_QC_conc,2)
      } else {
        best_model = cross_validation(as.numeric(QC_int),selected_QC_conc,3)
      }

      calibrated_int = calibrate_intensity(QC_int,selected_QC_conc,best_model,as.numeric(sample_table[i,]))[[1]]

      for (j in 1:(ncol(sample_table))) {
        if(sample_table[i,j] != 0 & as.numeric(calibrated_int[j]) == 0){
          best_model = cross_validation(as.numeric(QC_int),selected_QC_conc,1)
          calibrated_int = calibrate_intensity(QC_int,selected_QC_conc,best_model,as.numeric(sample_table[i,]))[[1]]
          break
        }
      }
      FeatureTable[i+1, c(FALSE,temp, FALSE)] = calibrated_int
      model[i] = model_pool[best_model+1]
      SQCpoint[i] = selected_int_points
      setTxtProgressBar(pb, i)
    }

  }
  FeatureTable$Model = c("Model", model)
  FeatureTable$SQC_points = c("SQC_points", SQCpoint)
  close(pb)

  if (output) {
    write.csv(FeatureTable, "Intensity_correction.csv", row.names = FALSE)
  }

  return(FeatureTable)
}
