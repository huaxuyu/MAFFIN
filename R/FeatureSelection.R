#' High-quality feature selection
#' @description
#' Select high-quality features for quantitative analysis.
#'
#' @param FeatureTable Feature intensity table with samples in column and features in row (default).
#' @param BlankFilter High-quality when mean(sample intensities) > mean(blank intensities) * \code{BlankFilter}
#' @param RtRange Range of the defined retention time window, in minute.
#' @param QCRSD Relative standard deviation threshold for QC samples
#' @param SQCcor Pearson's correlation threshold for serial QC samples (recommend: 0.8-0.9)
#' @param SampleInCol \code{TRUE} if samples are in column.
#' @param output \code{TRUE} will output the result table in current working directory
#'
#' @details  The first row should be sample names, and the second row should be group names.
#' For group names, please use: "RT" for retention time column,
#' "QC" for quality control samples between real samples (normal QC samples);
#' "blank" for blank samples;
#' "SQC_#amount#" for serial QC samples with a certain loading amount. For example, SQC_1.0 means a serial QC sample with injection volume as 1.0 uL;
#'
#' @return
#' This function will return the original feature table with an extra column named "Quality" to indicate the feature quality
#'
#' @examples
#' Please see GitHub for demo.


FeatureSelection = function(FeatureTable, BlankFilter=2, RtRange=c(0,100),
                            QCRSD=0.25, SQCcor=0.9, SampleInCol=TRUE, output=FALSE){

  f_table = read.csv(inputfile_name)
  group_seq = as.character(f_table[1,])
  group_unique = unique(group_seq[-1])
  f_table = as.matrix(sapply(f_table[-1,], as.numeric))

  derep_table = data.frame(matrix(data = 1, nrow = nrow(f_table), ncol = 5))
  colnames(derep_table) = c("Identifier", "MB", "RT", "QC_RSD", "Pearson_cor")
  derep_table[,1] = f_table[,1]

  message("Selecting high-quality features...")

  blank_table = c()

  if(length(grep("blank", group_seq, ignore.case = T)) != 0 & !is.na(blank_filter)){
    blank_table = data.matrix(f_table[, grep("blank", group_seq, ignore.case = T)])
    blank_table = apply(blank_table, 1, mean)
  } else if(length(grep("blank", group_seq, ignore.case = T)) == 0 & !is.na(blank_filter)){
    warning("Missed blank data. Filtering failed.")
  } else if(length(grep("blank", group_seq, ignore.case = T)) != 0 & is.na(blank_filter)){
    warning("Missed blank_filter parameter. Filtering failed.")
  }


  if(length(grep("SQC", group_seq, ignore.case = T)) != 0 & MRC){
    SQC_index = grep("SQC", group_unique, ignore.case = T)
    SQC_table = data.frame(matrix(nrow = nrow(f_table), ncol = length(SQC_index)))
    for (i in 1:length(SQC_index)) {
      SQC_table[,i] = apply(data.matrix(f_table[,group_unique[SQC_index[i]] == group_seq]), 1, mean)
    }
    colnames(SQC_table) = group_unique[SQC_index]
  } else if(length(grep("SQC", group_seq, ignore.case = T)) == 0 & MRC){
    warning("Serial QC data are missed or not correctly labeled.")
  } else if(length(grep("SQC", group_seq, ignore.case = T)) != 0 & !MRC){
    warning("MRC function is not turned on.")
  }

  QC_conc = c()
  SQC.list = stringr::str_split(colnames(SQC_table), pattern = "_")
  for (i in 1:ncol(SQC_table)) {
    QC_conc[i] = as.numeric(SQC.list[[i]][2])
  }

  temp = grepl("QC_mid", group_seq ,ignore.case = T) | grepl("SQC", group_seq,
                                                             ignore.case = T) | grepl("blank", group_seq, ignore.case = T) | grepl("rt", group_seq, ignore.case = T)
  sample_table = data.matrix(f_table[, !temp])
  sample_table2 = as.data.frame(sample_table)

  group_vector = group_seq[!temp][-1]

  RT_v = data.matrix(f_table[, grep("RT", group_seq, ignore.case = T)])

  if(length(grep("QC_mid", group_seq, ignore.case = T)) != 0){
    MQC_table = data.matrix(f_table[, grep("QC_mid", group_seq, ignore.case = T)])
    if(!is.na(blank_filter)){
      for (i in 1:nrow(f_table)) {
        if(mean(MQC_table[i,]) < blank_table[i]*blank_filter) { derep_table[i,2] = 0 }
      }
    }
    if(length(grep("QC_mid", group_seq, ignore.case = T)) >= 3){
      for (i in 1:nrow(f_table)) {
        if(mean(MQC_table[i,] == 0)){
          derep_table[i,4] = 0
        } else {
          if(sd(MQC_table[i,])/mean(MQC_table[i,]) > RSD){
            derep_table[i,4] = 0}
        }
      }
    } else {
      warning("QC data between samples are not enough to calculate RSD.")
    }

  } else {
    if(!is.na(blank_filter)){
      for (i in 1:nrow(f_table)) {
        if(mean(sample_table[i,]) < blank_table[i]*blank_filter) { derep_table[i,2] = 0 }
      }
    }
    warning("QC data between samples are missed or not correctly labeled.")
  }


  # Format of clean-up table. 1 for pass, 0 for fail.
  # column 1: feature identifier
  # column 2: method blank filter result
  # column 3: retention time filter result
  # column 4: RSD in QC sample filter result
  # column 5: QC linearity filter result

  all_filter = c()
  if(length(RT) == 2){
    for (i in 1:nrow(f_table)) {
      # Check retention time
      if(RT_v[i] < RT[1] | RT_v[i] > RT[2]) { derep_table[i,3] = 0 }
    }
  }

  for (i in 1:nrow(f_table)) {
    # Check serial diluted QC linearity
    QC_int = as.numeric(SQC_table[i,])
    valid_int = which(QC_int > int_threshold)
    QC_int = QC_int[valid_int]
    selected_QC_conc = QC_conc[valid_int]
    selected_int_points = length(selected_QC_conc)

    if(selected_int_points <= 3){
      derep_table[i,5] = 0
    } else if(cor(QC_int, selected_QC_conc) < QC_cor){
      derep_table[i,5] = 0
    }
    all_filter[i] = sum(derep_table[i,-1])
  }

  sample_table = sample_table[all_filter == 4,]
  SQC_table = SQC_table[all_filter == 4,]

  message("High-quality feature selection is done.")
  message(paste0(nrow(sample_table), " features are selected from ", nrow(f_table), "."))



}

