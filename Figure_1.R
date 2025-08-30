# Priyanka Rajashekar and Chuck Kroll
# SUNY College of Environmental Science and Forestry (ESF), Syracuse, NY, USA
# Contact: prajashe@syr.edu and cnkroll@esf.edu
# Modified Mood's Median Test on FDC Percentiles
# Figure 1: Annual Flow Duration Curves (FDCs) from daily streamflow records
# Last Updated: August 2025

# This example produces Figure 1 from the recently submitted article: 
# Rajashekar, P, Kroll, C.N., and Vogel, R.M. (2025) 
# Extending the MMMT temporal change detection tool to assess significant shifts 
# in diverse environmental variables, submitted to JAWRA August 2025.

# For this 1st case study, daily streamflow records from the U.S. Geological Survey 
# (USGS) Hydro-Climatic Data Network (HCDN) (https://waterdata.usgs.gov/nwis/sw) at 
# Sycamore Creek near Fort McDowell, AZ (09510200) are obtained. The dataset spans from 1980 to 2023. 

# To perform this analysis at another USGS streamflow site, the user would need to change
# information in the "Data Import" (line 69) and "Plotting" (line 133) sections of this script.

# Variable Description:
# site = USGS #, as a character string
# data = Input data for the variable of interest
# data1a = daily streamflow downloaded from USGS site for pre-disturbance/early record
# data2a = daily streamflow downloaded from USGS site for post-disturbance/late record
# data1 = data frame of pre-disturbance/early data (dates formatted as.Date in 1st column, numerical values in 2nd column)
# data2 = data frame of post-disturbance/late data (dates formatted as.Date in 1st column, numerical values in 2nd column)
# Start_Month = first month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# End_Month = last month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# min_quantile = minimum value for the FDC quantile range
# max_quantile = maximum value for the FDC quantile range
# increment = step size for the range of quantiles
# p = cdfs at which to estimate percentiles
# pre_data = list of pre-disturbance/early data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column)
# post_data = list of post-disturbance/late data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column)
# MMMT_out = matrix of output from Modified Mood's Median Test at each percentile (row) with columns of p-value, test statistic, and p-value flag
# MMMT_out1 = matrix of output from Modified Mood's Median Test for percentage of zeros in each year with columns of p-value, test statistic, and p-value flag
# FDR_flag = result of the field/global significance using a multiple False Discovery Rate (FDR) Test
# Pre_zeros = median annual percentage of zeros in the pre-disturbance/early record
# Post_zeros = median annual percentage of zeros in the post-disturbance/late record
# zerodata_flag = result of the significance on annual percentage of zeros using the MMMT Test
# label = array of labels for plot where [1] = plot table, [2] = x-axis title, [3] = y-axis title, [4] = data set 1 (pre-disturbance/early), and [5] = data set 2 (post-disturbance/late)
# yr1_1/yr1_2 = year of the first day of data1/data2
# yrn_1/yrn_2 = year of the last day of data1/data2
# ex = T/F flag where T = x-axis plots exceedance probabilities (1-cdf) and F = axis plots cdfs
# logy = T/F flag where when T = y-axis is plotted on a log scale and F = y-axis is not plotted on a log scale
# FDR = T/F flag where when T = FDR test results are plotted on figure and F = FDR test results are not plotted on figure
# ymin = minimum value on y axis
# ymax = maximum value on y axis
# y_steps = Step size to customize the y-axis
# yaxis_custom = T/F flag where when T = user can customize the step size of the y-axis, F= keeps the original axis logic exactly as is in the plot function
# pdfname = name of pdf file that the plot should be saved as

###
# Load packages and external functions
###

  # Load required packages 
  
  library(Hmisc)
  library(dataRetrieval)
  library(tidyr)
  library(dplyr)
  
  # Source functions in MMMT.R
  
  source("MMMT.R")

###
# Data Import
###

  # Enter the site number for Sycamore Creek near Fort McDowell, AZ
  # [For different USGS site in the NWIS system, the site number will change]
  
  site<-c("09510200") 
  
  # Download streamflow data (cfs) through the USGS portal using the DataRetrieval Package 
  # for specific time periods
  # [Can change date ranges here] 
  
  data1a<-readNWISdv(site[1],"00060","1980-10-01","1998-9-30")
  data1<-data1a[,3:4]
  data2a<-readNWISdv(site[1],"00060","1999-10-01","2023-9-30")
  data2<-data2a[,3:4]
  info<-readNWISsite(siteNumber = site[1])

###
# Preprocess pre- and post-disturbance/late data into a matrix with estimates of quantiles (percentiles; p)
# for each year or specific time period (e.g., water year) in each column
###

  # Enter the Start month and end month
  
  Start_Month<-10
  End_Month<-9
  
  # User must specify the range of FDC percentiles to estimate based on 1/(nmin+1) to nmin/(nmin+1)
  # where nmin is the minimum number of samples in any one year
  # for instance for a dataset with nmin = 365 points in an year, min can be equal to 0.005, 
  # max can be equal to 0.995 and increment = 0.005
  # Note: User can specify the quantiles as well
  
  min_quantile <- 0.005
  max_quantile <- 0.995
  increment <- 0.005
  p <- seq(min_quantile, max_quantile, increment)
  
  # Call function to preprocess pre-disturbance/early and post-disturbance/late data
  # (data between Start_Month and End_Month in each year in columns) and 
  # determine FDC percentiles for each year
  
  pre_data<-prepro(data1,Start_Month,End_Month,p)
  post_data<-prepro(data2,Start_Month,End_Month,p)

###
# MMMT Routine
###

  # Flag where when T = performs MMMT on log transformed data and F = performs MMMT on data in real space

  logsp_rlsp <- FALSE
  
  # Run MMMT on each percentile of pre_data and post_data and also on the annual percentage of zeros for datasets with zeros
  
  MMMT_out <- MMMT(pre_data,post_data,logsp_rlsp)
  
  # Calculate field/global significance using a multiple False Discovery Rate (FDR) Test
  # [3 = p-value < 0.01 ; 2 = 0.01 < p-value < 0.05 ; 1 = 0.05 < p-value < 0.10 ; 0 = p-value > 0.10] ]
  
  FDR_flag <- FDR(MMMT_out[[1]])

###
# Plotting
###
# Figure 1  
###

  # Set the name of the pdf to be produced
  
  pdfname <-paste0("Figure_1.pdf")   
  
  # Set labels for plot

  label<-as.character(array(0,5))
  label[1]<-paste(info[3],": USGS Site #",site[1],sep="")  # for streamflows
  label[2]<- paste0("Streamflow (cfs)")
  label[3] <- paste0("Exceedance Probability")
  
  # Enter the start and end years for the pre-disturbance/early and post-disturbance/late records
  
  yr1_1<-as.POSIXlt(data1[1,1])$year+1900
  yrn_1<-as.POSIXlt(data1[dim(data1)[1],1])$year+1900
  yr1_2<-as.POSIXlt(data2[1,1])$year+1900
  yrn_2<-as.POSIXlt(data2[dim(data2)[1],1])$year+1900
  
  label[4]<-paste("Early-record ",yr1_1,"-",yrn_1,sep="")
  label[5]<-paste("Late-record ",yr1_2,"-",yrn_2,sep="")
  
  # T/F flag where T = exceedance probabilities and F = non-exceedance probabilities
  
  ex<-TRUE
  
  # T/F flag where T = y-axis on a log scale and F = y-axis in real-space
  
  logy<-FALSE
  
  # T/F flag where T = FDR test results plotted on figure and F = FDR test results not plotted on figure
  
  FDR_In<-TRUE
  
  # T/F flag where T = MMMT test results for median annual percentage of zeros plotted on 
  # figure and F = not plotted on figure
  
  pvalueflag_In<-TRUE
  
  # T/F flag where T = customize the stepsize of the y-axis 
  # figure and F = keep the y-axis logic as is in the function
  
  yaxis_custom <- TRUE
  
  # Set the ymin and ymax
  
  ymin <- 0
  ymax <- 50
  
  # Customize the y-axis step size
  
  y_steps <- c(10)
  
  
  # Create plot and export as a pdf (name of pdf defined as pdfname) 
  
  Plot_Func(pre_data$Percentiles,post_data$Percentiles,MMMT_out[[1]],p,label,logy,ymin,ymax,ex,FDR_In,FDR_flag,pvalueflag_In, yaxis_custom, y_steps, pdfname)

  