
#' MAFFFIN normalization
#'
#' @description
#' Perform sample normalization using MAFFIN algorithm. MAFFIN algorithm consists
#' of three modules: high-quality feature selection, MS signal intensity correction,
#' and maximal density fold change normalization.
#'
#' @param FeatureTable Data frame with features in row and samples in column (default).
#' @param BlankFilter A numeric value. High-quality when mean(sample intensities) > mean(blank intensities) * \code{BlankFilter}.
#' @param RtRange A numeric vector indicating the range of the defined retention time window, in minute.
#' @param QCRSD A numeric value indicating the relative standard deviation threshold for QC samples.
#' @param SQCcor A numeric value indicating the Pearson's correlation threshold for serial QC samples (recommend: 0.8-0.9).
#' @param IntThreshold A numeric value indicating the feature intensity threshold. Feature is detected when its intensity larger than this value.
#' @param LR_QC_points Minimum serial QC data points for quadratic regression.
#' @param QR_QC_points Minimum serial QC data points for cubic regression.
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FALSE} if samples are in row.
#' @param output \code{TRUE} will output the result table in current working directory.
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
#'
#' @return
#' This function will return a list that contains four items if \code{RunEvaluation = TRUE}:
#' the normalized feature table, normalization factors, PRMAD of original data,
#' and PRMAD of normalized data. The last two items will not be generated if
#' \code{RunEvaluation = FALSE}
#'
#' @export
#'
#' @references Yu, Huaxu, and Tao Huan. "MAFFIN: Metabolomics Sample Normalization
#' Using Maximal Density Fold Change with High-Quality Metabolic Features and Corrected
#' Signal Intensities." \emph{bioRxiv} (2021).
#'
#' @examples
#' MAFFINTable = MAFFINNorm(TestingData)


MAFFINNorm = function(FeatureTable, BlankFilter=2, RtRange=c(0,100),
                      QCRSD=0.25, SQCcor=0.9, IntThreshold=0, LR_QC_points=5,
                      QR_QC_points=7, SampleInCol=TRUE, output=FALSE){

  f1 = FeatureSelection(FeatureTable=FeatureTable, BlankFilter=BlankFilter, RtRange=RtRange,
                        QCRSD=QCRSD, SQCcor=SQCcor, IntThreshold=IntThreshold,
                        SampleInCol=SampleInCol)

  f2 = IntCorrection(FeatureTable=f1, IntThreshold=IntThreshold, LR_QC_points=LR_QC_points,
                     QR_QC_points=QR_QC_points, SQCcor=SQCcor, SampleInCol=SampleInCol)

  f3 = MDFCNorm(FeatureTable=f2, IntThreshold=IntThreshold, SampleInCol=TRUE)

  if (output) {
    write.csv(f3[[1]], "MAFFIN_result.csv", row.names = FALSE)
  }

  return(f3)
}
