
#' MAFFFIN normalization
#'
#' @param FeatureTable Feature intensity table with samples in column and features in row (default).
#' @param BlankFilter High-quality when mean(sample intensities) > mean(blank intensities) * \code{BlankFilter}.
#' @param RtRange Range of the defined retention time window, in minute.
#' @param QCRSD Relative standard deviation threshold for QC samples.
#' @param SQCcor Pearson's correlation threshold for serial QC samples (recommend: 0.8-0.9).
#' @param IntThreshold Feature intensity threshold. Feature is detected when its intensity larger than this value.
#' @param LR_QC_points Required data points for quadratic regression (>= this value)
#' @param QR_QC_points Required data points for cubic regression (>= this value)
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FLASE} if samples are in row.
#' @param output \code{TRUE} will output the result table in current working directory.
#'
#' @return
#' This function will return a list contains two items: the MAFFIN result table,
#' and a set of normalization factors.
#' @export
#'
#' @examples
#' Please see GitHub for demo.


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
    write.csv(f3[[1]], "MAFFIN_result.csv")
  }

  return(f3)
}
