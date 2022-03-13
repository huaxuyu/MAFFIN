# This is to store the internal functions in MAFFIN package

# Calculate the pooled relative median absolute deviation
pooled_rMAD = function(data, group_vector) {
  unique_g = unique(group_vector)
  pMAD = 0
  for (i in 1:length(unique_g)) {
    d = data[group_vector == unique_g[i]]
    pMAD = pMAD + (length(d)-1)*mad(d)^2
  }
  pMAD = sqrt(pMAD / (length(group_vector) - length(unique_g)))
  prMAD = pMAD / mean(data)
  return(prMAD)
}


# Calculate the pooled relative standard deviation
pooled_rsd = function(data, group_vector) {
  unique_g = unique(group_vector)
  psd = 0
  for (i in 1:length(unique_g)) {
    d = data[group_vector == unique_g[i]]
    psd = psd + (length(d)-1)*var(d)
  }
  psd = sqrt(psd / (length(group_vector) - length(unique_g)))
  prsd = psd / mean(data)
  return(prsd)
}


# Optimize bandwidth
bw_opt = function(df, group_vector, bw_seq = seq(0.1,10,0.1)) {
  message("Optimizing bandwidth...")
  pb = txtProgressBar(min = 1, max = length(bw_seq), style = 3)
  pRMAD = c()

  for (b in 1:length(bw_seq)) {
    opt_df = df
    for (i in 2:ncol(df)) {
      fc_seq = as.numeric(df[,i]) / as.numeric(df[,1])
      fc_seq = fc_seq[as.numeric(df[,i]) != 0 & as.numeric(df[,1]) != 0]
      d = density(log2(fc_seq), bw = bw_seq[b])
      norm_factor = 2^d$x[which.max(d$y)]
      opt_df[,i] = as.numeric(df[,i]) / norm_factor
    }
    pRMAD_each = c()
    for (i in 1:nrow(opt_df)) {
      sample_data = as.numeric(opt_df[i,])
      pRMAD_each[i] = pooled_rMAD(sample_data, group_vector)
    }
    pRMAD[b] = median(pRMAD_each[!is.nan(pRMAD_each)])
    setTxtProgressBar(pb, b)
  }
  close(pb)
  return(bw_seq[which.min(pRMAD)])
}

#Create function for cross-validation
cross_validation = function(intensity_seq,conc_seq,order_number) {
  comparison_result = c()
  for (p in 2:(length(intensity_seq)-3)) {
    for (q in (p+2):(length(intensity_seq)-1)) {
      real_FC = conc_seq[q]/conc_seq[p]

      valid_data1 = as.numeric(intensity_seq[p])
      valid_data2 = as.numeric(intensity_seq[q])

      training_intensity_seq = intensity_seq[-c(p,q)]
      training_conc_seq = conc_seq[-c(p,q)]

      #Uncalibrated Ratio
      Uncali_FC = as.numeric(valid_data2/valid_data1)

      #Linear regression
      Linear_valid_data = calibrate_intensity(training_intensity_seq,training_conc_seq,1,c(valid_data1,valid_data2))
      Linear_calibrated_FC = Linear_valid_data[[1]][2]/Linear_valid_data[[1]][1]

      #Quadratic regression
      if(order_number >= 2){
        Quadratic_valid_data = calibrate_intensity(training_intensity_seq,training_conc_seq,2,c(valid_data1,valid_data2))
        Quadratic_Calibrated_FC = Quadratic_valid_data[[1]][2]/Quadratic_valid_data[[1]][1]
      } else{Quadratic_Calibrated_FC = 10000}


      #Cubic regression if order number is 3
      if(order_number >= 3){
        Cubic_valid_data = calibrate_intensity(training_intensity_seq,training_conc_seq,3,c(valid_data1,valid_data2))
        Cubic_Calibrated_FC = Cubic_valid_data[[1]][2]/Cubic_valid_data[[1]][1]
      } else{Cubic_Calibrated_FC = 10000}

      FC_diff = abs(c(Uncali_FC,Linear_calibrated_FC,Quadratic_Calibrated_FC,Cubic_Calibrated_FC)-real_FC)
      comparison_result = c(comparison_result, match(min(FC_diff),FC_diff))
    }
  }
  #Use lower order of regression if two models show same performance in cross-cvalidation
  return(as.numeric(names(sort(table(comparison_result),decreasing=TRUE)[1]))-1)
}


# Create function for MS signal calibration
# Need to input selected QC intensities, selected QC concentrations, regression model, and intensities for calibration
calibrate_intensity = function(s_QC_int, s_QC_conc, model_level, real_intensities){

  if(model_level != 0){
    calibrated_intensity = c()

    Re_coeff = lm(as.numeric(s_QC_int) ~ poly(s_QC_conc, model_level, raw = T))$coefficients

    for (int in 1:length(real_intensities)) {
      cali_int = 0
      if(real_intensities[int] == 0) {
        calibrated_intensity = c(calibrated_intensity,cali_int)
        next
      }
      Re_equation = polynom::polynomial(c(Re_coeff[1]-real_intensities[int],Re_coeff[2:length(Re_coeff)]))
      All_solutions = solve(Re_equation)
      All_solutions = Re(All_solutions[which(Im(All_solutions) == 0)])

      pre_cali_int = (All_solutions[All_solutions < tail(s_QC_conc,1) & All_solutions > 0])

      if(length(pre_cali_int) == 0){pre_cali_int = (All_solutions[All_solutions < 1.5*tail(s_QC_conc,1)])}
      if(length(pre_cali_int) == 0){pre_cali_int = (All_solutions[All_solutions < 2*tail(s_QC_conc,1)])}
      if(length(pre_cali_int) != 0){cali_int = max(pre_cali_int)}

      if(model_level == 1){cali_int = All_solutions}

      if(cali_int<0){cali_int = real_intensities[int]/s_QC_int[1]*s_QC_conc[1]}
      calibrated_intensity = c(calibrated_intensity,cali_int)
    }
    calibrated_intensity = list(calibrated_intensity)} else {calibrated_intensity = list(real_intensities)}

  return(calibrated_intensity)
}
