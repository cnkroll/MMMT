# Priyanka Rajashekar, Chuck Kroll, and Mohamed Mowafy
# SUNY ESF, Syracuse, NY, USA
# Contact: prajashe@syr.edu and cnkroll@esf.edu
# Modified Mood's Median Test
# Last Updated: August 2025

###
# This script contains 5 functions:
# prepro = preprocessing function which partitions annual data between the start and end month and calculate percentiles and percentage of zeros.
# MMMT = loop to perform Modified Mood's Median Test (MMMT) across a range of percentiles and percentage of zeros for datasets with zeros
# MMMT_once = calculates MMMT test statistic and p-values and flags for p-value ranges for a single sample
# FDR = calculates the field/global significance using a multiple False Discovery Rate (FDR) Test
# Plot_Func = routine to plot MMMT output
# Plot_Func2 = routine to plot MMMT output using the quantile ratio method
#
# These functions perform the MMMT and FDR test in:
# Mowafy, M, Kroll, C.N., and Vogel, R.M. (2023) 
# A Graphical Non-Parameteric Hydrologic Alteration Test Using Flow Duration Curves, 
# https://doi.org/10.1080/02626667.2023.2248109
# and in the recently submitted article:
# Rajashekar, P, Kroll, C.N., and Vogel, R.M. (2025) 
# Extending the MMMT temporal change detection tool to assess significant shifts 
# in diverse environmental variables, submitted to the Journal of the American Water
# Resources Association (JAWRA) August 2025.
###

###
# Function: prepro
###
# Function to take data frame of input data and partition data from
# Start_Month to End_Month into columns of a matrix. Note that Start_Month
# can be less than End_Month (within year) or greater than End_Month (over year).
# Then using the data in each column, the percentiles at different values of p are 
# calculated using Harrell-Davis Distribution-Free Quantile Estimator and the percentage of zeros are calculated for each year for datasets with zeros.
###

prepro<-function(data,Start_Month,End_Month,p){
  
  # Variable Description:
  # data = input data frame
  # Start_Month = starting month in which to begin analysis (1 = Jan, 2 = Feb, . . . 12 = Dec)
  # End_Month = ending month in which to end analysis (1 = Jan, 2 = Feb, . . . 12 = Dec)
  # p = cdfs at which to calculate percentiles
  # data = input data frame (dates formatted as.date in 1st column, numerical values in 2nd column)
  # data[,3] = month of date
  # data[,4] = year of date
  # datanew1 = data frame of rows of data between Start_Month and End_Month (Start_Month < EndMonth)
  # datanew1a = data frame of rows at Start_Month or before (Start_Month > End_Month)
  # datanew1b = data frame of rows at End_Month or after (Start_Month > End_Month)
  # datanew1c = merged datanew1a and datanew1b with a new column of the year of the period (Year_P)
  # datanew2 = list of datanew1, with each item in list a different year/time period of datanew1
  # output1 = matrix of data in order of occurence in each year (column) in data set
  # output2 = matrix of FDC percentiles estimated from output1 for each year
  # output 3 = matrix of percentage of zeros estimated from output1 for each year
  
  # Determine the month and year of each date
  
  data[,3]<-as.numeric(format(data[,1],format="%m"))
  data[,4]<-as.numeric(format(data[,1],format="%Y"))
  colnames(data)<-c("Date","Data","Month","Year")
  
  # Handle situations where Start_Month < End_Month and Start_Month > End_Month
  
  if (Start_Month < End_Month){
    
    datanew1<-data[(data[,3] >= Start_Month) & (data[,3] <= End_Month),]
    datanew2<-split(datanew1,f=datanew1$Year) 
    
  } else {
    
    datanew1a<-data[(data[,3] >= Start_Month),  ]
    datanew1a[,5]<-datanew1a$Year+1
    datanew1b<-data[(data[,3] <= End_Month), ]
    datanew1b[,5]<-datanew1b$Year
    colnames(datanew1a)[5]<-c("Year_P")
    colnames(datanew1b)[5]<-c("Year_P")
    datanew1c<-rbind(datanew1a,datanew1b)
    datanew2<-split(datanew1c,f=datanew1c$Year_P)
    
  }
  
  # output1 will have the raw data in each year (columns)
  # Limit amount of data in year to 10000 values in any one year
  
  n<-length(datanew2)
  output1<-matrix(NA,10000,n)
  
  for (i in seq(1,length(datanew2))){
    for (j in seq(1,dim(datanew2[[i]])[1])){
      output1[j,i]<-datanew2[[i]]$Data[j]
    }
  }
  
  # output2 will have the estimated percentile in each year (columns)
  
  output2<-matrix(NA,length(p),dim(output1)[2])
  
  for (i in seq(1,dim(output1)[2])){
    output2[,i]<-hdquantile(na.omit(output1[,i]), probs = p)
  }
  
  # Set columm names for output1 and output2 and row names for output2
  
  colnames(output1)<-names(datanew2)
  colnames(output2)<-names(datanew2)
  rownames(output2)<-p
  
  # Remove all rows of output1 which are all NAs
  
  NAall<-apply(output1,1,function(x) all(is.na(x)))
  output1<-output1[!NAall,]
 
  ## For datasets with zeros
  
  # output3 will have the percentage of zeroes in each year (columns)
  
  output3 <- matrix(0,1,dim(output1)[2])
  
  for (i in seq(1,dim(output1)[2])){
    
    # Check if any value is 0, excluding NA values
    if (any(output1 == 0, na.rm = TRUE)) {
      
      num_zeroes <- sum(output1[,i] == 0, na.rm = TRUE)
      
      total_obs <- sum(!is.na(output1[,i]))
      
      output3[,i]<- round((num_zeroes/total_obs)*100, digits = 2)
      
      colnames(output3)<-names(datanew2)
      
      # Setup output list
      
      output<-list("Percentiles"=output2,"Data"=output1,"Zero Data"=output3)
      
      
    } else {

      # Output for datasets with no zeros
            
      output<-list("Percentiles"=output2,"Data"=output1)
      
    }
    
  }
  
  return(output) 
  
}

###
# Function: MMMT
###
# Function to read in preprocessed data from 2 samples and to perform
# Modified Moods Median Test for each percentile or percentage of zeros (using the function MMMT_once), 
# returning p-value, test statistic and whether p-values are within flagged 
# ranges (< 0.01,0.01-0.05, and 0.05-0.10)
###

# Loop through each percentile and perform MMMT on each percentile
# Perform MMMT on percentage of zeros if dataset has zero values

MMMT<-function(pre,post,logsp_rlsp){
  
  #Variable Description:
  # logsp_rlsp = T/F flag where when T = performs MMMT on log transformed data 
  # F = performs MMMT on data in real space
  
  # Check if there are 0's in the either of the records, excluding NA values
  
  if (any(c(pre$Data, post$Data) == 0, na.rm = TRUE)) {
    
  Np<-dim(pre$Percentiles)[1]
  MMMT_out<-matrix(0,Np,3)
  
  Nz<-dim(pre$`Zero Data`)[1]
  MMMT_out1 <- matrix(0,Nz,3)
  
  for (i in seq(1,Np)){
    MMMT_out[i,]<-MMMT_once(pre$Percentiles[i,],post$Percentiles[i,],logsp_rlsp)
  }
  
  for (j in seq(1,Nz)){
    MMMT_out1[j,]<-MMMT_once(pre$`Zero Data`[j,],post$`Zero Data`[j,],logsp_rlsp)
  }
  
  # Set column and row names of MMMT_out
  
  colnames(MMMT_out)<-c("p-value","Test Statistic","p-value flag")
  rownames(MMMT_out)<-rownames(pre)
  colnames(MMMT_out1)<-c("p-value","Test Statistic","p-value flag")
  
  Pre_zeros <- round(median(pre$`Zero Data`), digits = 2)
  Post_zeros <- round(median(post$`Zero Data`), digits = 2)
  
  # Formatting the pre_zeros to have consistent decimal places in Pre_zeros and Post_zeros outputs
  
  Pre_zeros <- sprintf("%.1f", Pre_zeros)
  Post_zeros <- sprintf("%.1f", Post_zeros)
  
  zerodata_flag <- MMMT_out1[,3]
  
  MMMT_output <- list("Percentiles_output" = MMMT_out,"ZeroData_output" = MMMT_out1,
                      "Median_pre_zeros" = Pre_zeros, "Median_post_zeros" = Post_zeros,
                      "Zerodata_pvalue_flag" = zerodata_flag)
  
  } else{
  
    Np<-dim(pre$Percentiles)[1]
    MMMT_out<-matrix(0,Np,3)
    
    for (i in seq(1,Np)){
      MMMT_out[i,]<-MMMT_once(pre$Percentiles[i,],post$Percentiles[i,],logsp_rlsp)
    }
    
    # Set column and row names of MMMT_out
  
    colnames(MMMT_out)<-c("p-value","Test Statistic","p-value flag")
    rownames(MMMT_out)<-rownames(pre)
    
  MMMT_output <- list("Percentiles_output" = MMMT_out)  
    
  }
  
  return(MMMT_output)
  
}

###
# Function: MMMT_once
###
# Function to perform the Modified Mood's Median Test (MMMT; Fligner and Rust, 1982)
# on 2 samples. This function is called by MMMT for annual values of each percentile.
###

# Perform MMMT on a single sample

MMMT_once<-function(F1, F2, logsp_rlsp){
  
  # Variable Description:
  # F1 = input array of length m 
  # F2 = input array of length n
  # logsp_rlsp = T/F flag where when T = performs MMMT on log transformed data 
  # F = performs MMMT on data in real space
  # p = m/n = ratio of sample sizes
  # N = m + n = total sample size 
  # combined_record = a combined record of F1 and F2
  # M = median of combined sample
  # bN = Modified Mood's Median Test spread indicator
  # Z_combined = sorted combined sample
  # UN = upper value of combined sample to determine distribution spread
  # LN = lower value of combined sample to determine distribution spread
  # F1_cdf_UN/LN = cdf of F1 evaluated at UN/LN
  # F2_cdf_UN/LN = cdf of F2 evaluated at UN/LN
  # R_hat= (F1_cdf_UN - F1_cdf_LN)/(F2_cdf_UN - F2_cdf_LN) = Ratio of spreads
  # Var_hat = variance of test statistic
  # class_1/0.5/0 = count of how many values of F2 are greater than/equal to/less than M
  # Initial_T = Mood's Median Test test stastic (called "T" in Mowafy et al. (2022))
  # T_hat = Modified Mood's Median Test test statistic
  # d.f = n + m = degrees of freedom of test
  # p_value = p-value of test 
  # output_array = p_value, T_hat, and flags of whether p-value is < 0.01 (3), 0.01-0.05 (2)
  #                0.05 - 0.10 (1), and > 0.10 (0)
  
  # Determine size of sample of F1 and F2 
  
  m<-length(F1)
  n<-length(F2)
  
  # Calculate the ratio of sample sizes p
  
  p<-m/n
  
  # Calculate length of combined sample N 
  
  N<-m + n
  
  # Combine the pre and post sample 
  
  combined_record<-c(as.numeric(F1),as.numeric(F2))
  
  # Calculate the median for the combined record
  
  M<-median(combined_record)
  
  # Calculate rank of UN, bN (and rank of LN, N - bN)
  
  bN<-floor((0.5 * (N + 1)) + sqrt(N))
  
  # Sort the combined sample
  
  Z_Combined<-sort(combined_record)
  
  # Determine the values of UN and LN by ranked order of combined sample
  
  LN<-Z_Combined[N - bN]
  UN<-Z_Combined[bN]
  
  # If LN is less than a min value in a sample, make it equal to the min value in the sample
  # If UN is more than a max value in a sample, make it equal to the max value in the sample 
  
  # Check to make sure LN < UN; if not set T_hat = 3 and p-value = 0.005
  
  if (LN == UN | median(F1) == median(F2)){
    
    T_hat <- 0
    p_value <- 1.0
    
   } else if (max(F1) < min(F2) | max(F2) < min(F1)){
    T_hat <- 4
    p_value <- 0.001
    
    } else {
    
    if (logsp_rlsp == TRUE){  
      
    # Determine the cdf of UN and LN in each sample by doing a linear interpolation across
    # the log of the flows using a weibull plotting position (i/(N+1))
      
    if (UN > max(F1)){F1_cdf_UN <- 1} else {
    F1_cdf_UN<-approx(log(F1),rank(F1)/(length(F1)+1), log(UN),ties = mean, na.rm = TRUE)$y}

    if (UN > max(F2)){F2_cdf_UN <- 1} else {
    F2_cdf_UN<-approx(log(F2),rank(F2)/(length(F2)+1), log(UN),ties = mean, na.rm = TRUE)$y}
      
    if ((LN < min(F1)) | (LN == 0)){F1_cdf_LN <- 0} else {   
    F1_cdf_LN<-approx(log(F1),rank(F1)/(length(F1)+1), log(LN),ties = mean, na.rm = TRUE)$y}
    
    if(LN > max(F1)){F1_cdf_LN <- 1}
      
    if ((LN < min(F2)) | (LN == 0)){F2_cdf_LN <- 0} else { 
    F2_cdf_LN<-approx(log(F2),rank(F2)/(length(F2)+1), log(LN),ties = mean, na.rm = TRUE)$y}
    
    if(LN > max(F2)){F2_cdf_LN <- 1}
      
    } else {
      
    if (UN > max(F1)){F1_cdf_UN <- 1} else {
    F1_cdf_UN<-approx(F1,rank(F1)/(length(F1)+1), UN,ties = mean, na.rm = TRUE)$y}
      
    if (UN > max(F2)){F2_cdf_UN <- 1} else {
    F2_cdf_UN<-approx(F2,rank(F2)/(length(F2)+1),UN,ties = mean, na.rm = TRUE)$y}
      
      
    if (LN < min(F1)){F1_cdf_LN <- 0} else if (LN == 0 && min(F1)==0){
      F1_cdf_LN<-length(which(F1==0))/(length(F1)+1)
      F1_woz<-F1[which(F1>0)]
      F1_cdf_UN<-F1_cdf_LN+(1-F1_cdf_LN)*approx(F1_woz,rank(F1_woz)/(length(F1_woz)+1),UN,ties=mean,na.rm=TRUE)$y
      } else {
      F1_cdf_LN<-approx(F1,rank(F1)/(length(F1)+1), LN,ties = mean, na.rm = TRUE)$y}
    
    if(LN > max(F1)){F1_cdf_LN <- 1}
      
    if (LN < min(F2)){F2_cdf_LN <- 0} else if (LN==0 && min(F2)==0){ 
      F2_cdf_LN<-length(which(F2==0))/(length(F2)+1)
      F2_woz<-F2[which(F2>0)]
      F2_cdf_UN<-F2_cdf_LN+(1-F2_cdf_LN)*approx(F2_woz,rank(F2_woz)/(length(F2_woz)+1),UN,ties=mean,na.rm=TRUE)$y
      } else {
    F2_cdf_LN<-approx(F2,rank(F2)/(length(F2)+1), LN,ties = mean, na.rm = TRUE)$y}  
    
    if(LN > max(F2)){F2_cdf_LN <- 1}
      
    } 
    # Calculate R_hat and Var_hat 
    # If G_cdf_UN = G_cdf_LN (or F_cdf_UN = F_cdf_LN), Var_hat = 0.25*p (see Mowafy et al. (2022) for discussion)
    
    if (isTRUE(F1_cdf_UN == F1_cdf_LN) || isTRUE(F2_cdf_UN == F2_cdf_LN)){Var_hat = 0.25*p} else {
      R_hat<-(F1_cdf_UN - F1_cdf_LN)/(F2_cdf_UN - F2_cdf_LN) 
    Var_hat<-(0.25 * p * (1 + (p * R_hat^2)))/((1 + (p * R_hat))^2)
    }
    
    # Set up counters to determine original Mood's median test T
    
    class_1<-0
    class_0.5<-0
    class_0<-0
    
    # Check how many values in the post sample that are (greater than, equal, or less than) the median of combined record
    
    for (i in seq(1, length(F2))){
      if (isTRUE(F2[i] < M)){ class_1<-class_1 + 1 }
      else if (isTRUE(F2[i] == M)) { class_0.5<-class_0.5 + 1 }
      else { class_0<-class_0 + 1 }
    }
    
    # Calculate the initial test statistic T
    
    Initial_T<-((class_1 * 1) + (class_0.5 * 0.5) + (class_0 * 0))/n
    
    # Calculate improved distribution free test statistiv T_hat
    
    T_hat<- (sqrt(n) * (Initial_T - 0.5))/ sqrt(Var_hat)
    
    # Assign the degrees of freedom to m+n
    
    d.f<- n + m
    
    # Calculate the p-value for two-sided test
    
    p_value<- 2 *(1 -  (pt(abs(T_hat), d.f, lower.tail = TRUE)))
    
  } 
  
  # Create array to store the results
  
  output_array<-array(0, 3)
  
  output_array[1]<-p_value
  output_array[2]<-T_hat
  
  # Calculate the significance level flag (1 if sig. @ 1%,   2 if sig. @ 5%,  3 if sig. @ 10%,  and zero if not  sign.)
  
  if (isTRUE(p_value < 0.01)) {output_array[3]<-3} else if (isTRUE(p_value > 0.01 & p_value <0.05)){
    output_array[3]<-2} else if (isTRUE(p_value > 0.05 & p_value < 0.1)){
      output_array[3]<-1} else {
        output_array[3]<-0}
  
  # Return the output of the function
  
  return(output_array) 
  
}

###
# Function: FDR
###
# Function to calculate the field/global significance using a multiple False Discovery Rate (FDR) Test
###

FDR<-function(MMMT_out){
  
  # Variable Description:
  # MMMT_out = output from performing the MMMT test at each percentile
  # p_val_sort = p-values from MMMT across all percentiles sorted from smallest to largest
  # alpha_glob = array of global significance type 1 errors (alpha) = 0.1, 0.05, and 0.01
  # N = number of p-values
  # FDR_flag = flag of results of MMMT test where [3 = p-value < 0.01 ; 2 = 0.01 < p-value < 0.05 ; 1 = 0.05 < p-value < 0.10 ; 0 = p-value > 0.10]
  
  # Sort p-values from high to low
  
  p_val_sort <- sort(MMMT_out[,1],decreasing = F)
  
  # Set global type 1 (alpha) values
  
  alpha_glob <- c(0.10,0.05,0.01)
  
  # Determine the number of p-values
  
  N <- length(p_val_sort)
  
  # Set FDR_flag to 0 and override this if field/global test p-values are < 0.1
  
  FDR_flag <- 0
  for (j in seq(1,length(alpha_glob))){
    for (i in seq(1,N)){
      if (p_val_sort[i]< i*alpha_glob[j]/N) FDR_flag <- j
    }          
  }
  
  return(FDR_flag)
  
}

###
# Function: Plot_Func
###
# Function to make plot of MMMT output
###

Plot_Func<-function(pre,post,mood_modified,p,label,logy,ymin,ymax,ex,FDR_In,FDR_flag, pvalueflag_In, yaxis_custom, y_step, pdfname){
  
  # Variable Description:
  # pre = matrix of annual percentiles for pre-disturbance record calculated from prepro
  # post = matrix of annual percentiles for post-disturbance record calculated from prepro
  # mood_modified_mat = matrix of output from MMMT
  # p = array of cdfs 
  # label = array of labels for plots (main, ylab, xlab)
  # p_ex = 1 - p = exceedance probabilities
  # mycol_lightblue = colors for pre annual
  # mycol_pink = colors for post annual 
  # Np = length of p_ex
  # FDR_In = T/F flag where when T = FDR test results are plotted on figure and F = FDR test results are not plotted on figure
  # pvalueflag_In = T/F flag where when T = Median Annual Percent zero are plotted on figure and F = Median Annual Percent zero are not plotted on figure
  # FDR_flag = result of the field/global significance using a multiple False Discovery Rate (FDR) Test
  # Final_df = data frame with data from p_ex, median of pre for each percentile (Median_pre), 
  #            median of post for each percentile (Median_post), and p-value, Test Statistic, 
  #            and p-value flag for MMMT test
  # Final_df_split = list made by splitting Final_df based on p-value flag
  
  
  # Turn cdfs into exceedance probabilities  
  
  pdf(pdfname, width = 8, height = 6.5)
  
  if (ex==TRUE){
    p_ex<-1-p} else {p_ex<-p}
  
  # Set plot colors
  
  mycol_pink <-rgb(1, 0.8, 0.8, 0.35)
  mycol_lightblue <-rgb(0.4, 0.8, 1, 0.25)
  
  # Determine how many percentiles
  
  Np<-length(p_ex)
  
  # Create matrix of p, medians at each percentile for pre and post, p-value, test stat, and flag 
  
  pre<-as.matrix(pre)
  post<-as.matrix(post)
  Final_df<-data.frame(p_ex,apply(pre,1,median),apply(post,1,median),mood_modified)
  colnames(Final_df)<-c("p_ex","Median_pre","Median_post","p-value","Test Statistic","Flag")
  
  # Split Final_df by sig. Level flag for the sake of plotting in different color
  
  Final_df_split<-split(Final_df, Final_df$Flag)
  
  # Set the plotting space for one graph
  
  par(mfrow=c(1,1), mai=c(1.7, 0.65, 0.35, 0.35)) 
  
  # Remove scientific notation on y-axis labels
  
  options(scipen=999)
  
  # Plot the first year of the pre record individual annual percentiles (FDCs)  
  # If logy is TRUE, plot y-axis on a log-base10-scale; if not, plot y-axis on a linear-scale
  
  if (isTRUE(logy)){plot(p_ex, pre[,1], col = mycol_lightblue, ylim = c(ymin,ymax), log = "y", type = 'l', lwd =1, xaxt = "n", yaxt = "n", ann = FALSE)} else 
  {plot(p_ex, pre[,1], col = mycol_lightblue, ylim = c(ymin,ymax), type = 'l', lwd =1, xaxt = "n", yaxt = "n", ann = FALSE)}
  
  # Loop from the 2nd year up to the last year of the pre record, plotting percentiles for individual years
  
  for (i in seq(2, ncol(pre))){
    lines(p_ex, pre[,i], col = mycol_lightblue,  type = 'l', lwd =1)  }
  
  # Loop from the 1st year up to the last year of the post record, plotting percentiles for individual years
  
  for (i in seq(1, ncol(post))){
    lines(p_ex, post[,i], col = "snow2",  type = 'l', lwd =1)  }
  
  lines(Final_df_split$`0`$p_ex, Final_df_split$`0`$Median_pre, col = "blue",  type = 'p',pch = 15, lwd = 0.9, cex=0.9)  
  
  lines(Final_df_split$`0`$p_ex, Final_df_split$`0`$Median_post, col = "black",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  # Add the median annual FDCs for both pre and post records for percentiles that are significant at type 1 error of 10%,
  # but not significant at 5% (p-value flag = 1; 0.05 < p-value < 0.1)
  
  lines(Final_df_split$`1`$p_ex, Final_df_split$`1`$Median_pre, col = "yellow",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  lines(Final_df_split$`1`$p_ex, Final_df_split$`1`$Median_post, col = "yellow",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  # Add the median annual FDCs for both pre and post records for percentiles that are significant at type 1 error of 5%,
  # but not significant at 1% (p-value flag = 2; 0.01 < p-value < 0.05)
  
  lines(Final_df_split$`2`$p_ex, Final_df_split$`2`$Median_pre, col = "orange",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  lines(Final_df_split$`2`$p_ex, Final_df_split$`2`$Median_post, col = "orange",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  # Add the median annual FDCs for both pre and post records for percentiles that are significant at type 1 error of 1%
  # (p-value flag = 3; p-value < 0.01)
  
  lines(Final_df_split$`3`$p_ex, Final_df_split$`3`$Median_pre, col = "red",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  lines(Final_df_split$`3`$p_ex, Final_df_split$`3`$Median_post, col = "red",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  # Add x-axis and y-axis
  
  axis(1, seq(0, 1, 0.1), tck = 0.01)
  
  # Add a flag for drawing the y-axis
  # if F = keeps the original axis logic exactly as is, T= skips this line entirely so the user can add a custom axis in the main script using axis()

  if (isFALSE(yaxis_custom)) {
    if (isTRUE(logy)) {
      axis(2, 10^seq(log10(ymin),log10(ymax)), tck = 0.01)
    } else {
      axis(2,seq(ymin,ymax,(ymax-ymin)/4))
      }
  } else if (!is.null(y_step)) {
    axis(2, seq(ymin, ymax, by = y_step))   # use step size
  }
  
  # Add plot title and axis labels

  title(main = label[1],cex.main=0.8)
  title(ylab = label[2],line = 2)
  title(xlab = label[3],line = 2) 
  
  # Add legends
  
  legend("bottomleft", inset=c(0,-0.35), xpd=TRUE , 
         c(paste("Median ",label[4],sep=""),paste("Median ",label[5],sep=""), 
           paste("Individual Annual ",label[4]," (",dim(pre)[2]," Years)",sep=""), 
           paste("Individual Annual ",label[5]," (",dim(post)[2]," Years)",sep="")), 
         lty = c(NA,NA,1,1), lwd = c(1.2,1.2,1.2,2,2), pch = c(0,0,NA,NA) , 
         col = c("blue", "black", mycol_lightblue, "snow2" ),cex=0.8,bty = "n")
  
  legend("bottomleft", inset=c(0.75,-0.35), xpd=TRUE , c("               p-value <= 0.01", "    0.01 < p-value <= 0.05", "    0.05 < p-value <= 0.1","               p-value > 0.1"), bty="n", col = c("red", "orange", "yellow","blue"), pch = c(15,15,15,15), y.intersp = 0.8, x.intersp = 0.6,cex=0.8)
  legend("bottomleft", inset=c(0.767,-0.35), xpd=TRUE , "" , bty="n", col = "black", pch = 15, y.intersp = 0.8, x.intersp = 0.6,cex=0.8)
  
  # Add FDR results if FDR_In = TRUE
  
  if (isTRUE(FDR_In)) {
    if (FDR_flag == 3) {col_FDR = "red"} else 
      if (FDR_flag == 2) {col_FDR = "orange"} else
        if (FDR_flag == 1) {col_FDR = "yellow"} else
        {col_FDR = "blue"}
    
    text(p_ex[floor(0.2*length(p_ex))],(ymin+0.6*(ymax-ymin)),"FDR p-value",col = col_FDR)    
    
  }
  
  if (isTRUE(pvalueflag_In)) {
    if (MMMT_out[[5]] == 3) {col_percentzero = "red"} else 
      if (MMMT_out[[5]] == 2) {col_percentzero = "orange"} else
        if (MMMT_out[[5]] == 1) {col_percentzero = "yellow"} else
        {col_percentzero = "blue"}
    
    text(p_ex[floor(0.2*length(p_ex))],(ymin+0.5*(ymax-ymin)),"Median Annual Percent Zero",col = col_percentzero)    
    text(p_ex[floor(0.2*length(p_ex))],(ymin+0.45*(ymax-ymin)),
         paste0("Early = ",MMMT_out[[3]], "%, ", "Late = ",MMMT_out[[4]],"%"),col = col_percentzero)
  }

  dev.off()
  
}

###
# Function: Plot_Func2
###
# Function to make plot of MMMT output using the Quantile ratio method
###

Plot_Func2<-function(pre,post,mood_modified,p,label,ymin2,ymax2,ex,FDR_In,FDR_flag, yaxis_custom, y_step, pdfname2){
  
  # Variable description:
  # pre = matrix of annual percentiles for pre-disturbance record calculated from prepro
  # post = matrix of annual percentiles for post-disturbance record calculated from prepro
  # mood_modified_mat = matrix of output from MMMT
  # p = array of cdfs 
  # label = array of labels for plots (main, ylab, xlab)
  # p_ex = 1 - p = exceedance probabilities
  # mycol_lightblue = colors for pre annual
  # mycol_pink = colors for post annual 
  # Np = length of p_ex
  # Final_df = data frame with data from p_ex, median of pre for each percentile (Median_pre), 
  #            median of post for each percentile (Median_post), and p-value, Test Statistic, 
  #            and p-value flag for MMMT test
  # Final_df_split = list made by splitting Final_df based on p-value flag
  
  # Turn cdfs into exceedance probabilities  
  
  pdf(pdfname2, width = 8, height = 6.5)
  
  if (ex==TRUE){
    p_ex<-1-p} else {p_ex<-p}
  
  # Set plot colors
  
  mycol_pink <-rgb(1, 0.8, 0.8, 0.35)
  mycol_lightblue <-rgb(0.4, 0.8, 1, 0.25)
  
  # Determine how many percentiles
  
  Np<-length(p_ex)
  
  # Create matrix of p, medians at each percentile for pre and post, p-value, test stat, and flag 
  
  pre <- pre_data$Percentiles
  pre<-as.matrix(pre)
  
  post <- post_data$Percentiles
  post<-as.matrix(post)
  
  ## Final medians for each quantile for pre and post
  ## Then divide post median by pre median 
  ## [All series will be divided by pre median]
  
  pre_med<-apply(pre,1,median)
  post_med<-apply(post,1,median)
  rat_med<-post_med/pre_med
  
  ## Create data frame to hold results
  
  mood_modified <- MMMT_out
  Final_df<-data.frame(p_ex,rat_med,mood_modified)
  
  colnames(Final_df)<-c("p_ex","Median_ratio","p-value","Test Statistic","Flag")
  
  # Split Final_df by sig. Level flag for the sake of plotting in different color
  
  Final_df_split<-split(Final_df, Final_df$Flag)
  
  # Set the plotting space for one graph
  
  par(mfrow=c(1,1), mai=c(1.7, 0.65, 0.35, 0.35)) 
  
  # Remove scientific notation on y-axis labels
  
  options(scipen=999)
  
  ## Divide individual year FDCs (pre and post) by median pre FDC
  
  d_pre<-dim(pre)
  pre_rat<-matrix(NA,d_pre[1],d_pre[2])
  for (i in seq(1,d_pre[2])){
    pre_rat[,i]<-pre[,i]/pre_med  
  }
  d_post<-dim(post)
  post_rat<-matrix(NA,d_post[1],d_post[2])
  for (i in seq(1,d_pre[2])){
    post_rat[,i]<-post[,i]/pre_med  
  }
  
  # Plot the first year of the pre record individual annual percentiles (FDCs)  
  # If logy is TRUE, plot y-axis on a log-base10-scale; if not, plot y-axis on a linear-scale
  
  if (isTRUE(logy)){plot(p_ex, pre_rat[,1], col = mycol_lightblue, ylim = c(ymin2,ymax2), type = 'l', lwd =1, xaxt = "n", yaxt = "n", ann = FALSE)} else 
  {plot(p_ex, pre_rat[,1], col = mycol_lightblue, ylim = c(ymin2,ymax2), type = 'l', lwd =1, xaxt = "n", yaxt = "n", ann = FALSE)}
  
  # Loop from the 2nd year up to the last year of the pre record, plotting percentiles for individual years
  
  for (i in seq(2, ncol(pre))){
    lines(p_ex, pre_rat[,i], col = mycol_lightblue,  type = 'l', lwd =1)  }
  
  ## Add a line across at 1
  
  abline(h=1,col="black")
  
  # Loop from the 1st year up to the last year of the post record, plotting percentiles for individual years
  
  for (i in seq(1, ncol(post))){
    lines(p_ex, post_rat[,i], col = "snow2",  type = 'l', lwd =1)  }
  
  # Add the median annual FDCs for both pre and post records for percentiles that are 
  # not significant at a type 1 error of 10% (p-value flag = 0; p-value > 0.1)
  
  lines(Final_df_split$`0`$p_ex, Final_df_split$`0`$Median_ratio, col = "black",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  # Add the median annual FDCs for both pre and post records for percentiles that are significant at type 1 error of 10%,
  # but not significant at 5% (p-value flag = 1; 0.05 < p-value < 0.1)
  
  lines(Final_df_split$`1`$p_ex, Final_df_split$`1`$Median_ratio, col = "yellow",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  # Add the median annual FDCs for both pre and post records for percentiles that are significant at type 1 error of 5%,
  # but not significant at 1% (p-value flag = 2; 0.01 < p-value < 0.05)
  
  lines(Final_df_split$`2`$p_ex, Final_df_split$`2`$Median_ratio, col = "orange",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  # Add the median annual FDCs for both pre and post records for percentiles that are significant at type 1 error of 1%
  # (p-value flag = 3; p-value < 0.01)
  
  lines(Final_df_split$`3`$p_ex, Final_df_split$`3`$Median_ratio, col = "red",  type = 'p', pch = 15, lwd = 0.9, cex=0.9)  
  
  # Add x-axis and y-axis
  
  axis(1, seq(0, 1, 0.1), tck = 0.01)
  
  # Add a flag for drawing the y-axis
  # if F = keeps the original axis logic exactly as is, T= skips this line entirely so the user can add a custom axis in the main script using axis()
  
  if (isFALSE(yaxis_custom)) {
    if (isTRUE(logy)) {
      axis(2, 10^seq(log10(ymin2),log10(ymax2)), tck = 0.01)
    } else {
      axis(2,seq(ymin2,ymax2,(ymax2-ymin2)/4))
    }
  } else if (!is.null(y_step)) {
    axis(2, seq(ymin2, ymax2, by = y_step))   # use step size
  }
  
  # Add legends
  
  legend("bottomleft", inset=c(0,-0.35), xpd=TRUE , 
         c(paste("Median ",label[4],sep=""),paste("Median ",label[5],sep=""), 
           paste("Individual Annual ",label[4]," (",dim(pre)[2]," Years)",sep=""), 
           paste("Individual Annual ",label[5]," (",dim(post)[2]," Years)",sep="")), 
         lty = c(NA,NA,1,1), lwd = c(1.2,1.2,1.2,2,2), pch = c(0,0,NA,NA) , 
         col = c("blue", "black", mycol_lightblue, "snow2" ),cex=0.8,bty = "n")
  
  legend("bottomleft", inset=c(0.75,-0.35), xpd=TRUE , c("               p-value <= 0.01", "    0.01 < p-value <= 0.05", "    0.05 < p-value <= 0.1","               p-value > 0.1"), bty="n", col = c("red", "orange", "yellow","blue"), pch = c(15,15,15,15), y.intersp = 0.8, x.intersp = 0.6,cex=0.8)
  legend("bottomleft", inset=c(0.767,-0.35), xpd=TRUE , "" , bty="n", col = "black", pch = 15, y.intersp = 0.8, x.intersp = 0.6,cex=0.8)
  
  # Add FDR results if FDR_In = TRUE
  
  if (isTRUE(FDR_In)) {
    if (FDR_flag == 3) {col_FDR = "red"} else 
      if (FDR_flag == 2) {col_FDR = "orange"} else
        if (FDR_flag == 1) {col_FDR = "yellow"} else
        {col_FDR = "blue"}
    
    text(p_ex[floor(0.2*length(p_ex))],(ymin2+0.6*(ymax2-ymin2)),"FDR p-value",col = col_FDR)    
    
  }
  
  # Add plot title and axis labels
  
  title(main = label[1],cex.main=0.8)
  title(ylab = label[2],line = 2)
  title(xlab = label[3],line = 2) 
  
  dev.off()
  
}
