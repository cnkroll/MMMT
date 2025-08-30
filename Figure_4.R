# Priyanka Rajashekar and Chuck Kroll
# SUNY College of Environmental Science and Forestry (ESF), Syracuse, NY, USA
# Contact: prajashe@syr.edu and cnkroll@esf.edu
# Modified Mood's Median Test on DC Percentiles
# Figure 4: Annual Duration Curves (DCs) from monthly water quality records
# Last Updated: August 2025

# This example produces Figure 4 from the recently submitted article: 
# Rajashekar, P, Kroll, C.N., and Vogel, R.M. (2025) 
# Extending the MMMT temporal change detection tool to assess significant shifts 
# in diverse environmental variables, submitted to JAWRA August 2025.

# For this 4th case study, water quality monitoring data (Total nitrogen (TN) concentrations) from the Chesapeake Bay Nontidal 
# Network (NTN; https://datahub.chesapeakebay.net/WaterQuality) station at Conestoga 
# River (WQN0273) in Lancaster County, PA, spanning from 2006 to 2022 are obtained.

# To perform this analysis at another site, the user would need to change information in 
# the "Data Import" (line 71) and "Plotting" (line 151) sections of this script. In this 
# example, data from the Chesapeake Bay NTN was downloaded and saved as a csv file, which was 
# read into R. See "Waterquality_data.csv" for file structure and column headings.

# Variable Description:
# data = data for the site selected for the entire period of record
# Start_date = first date in the record
# End_date = last date in the record
# Start_year = first year in the record
# End_year = last year in the record
# Difference_in_years = difference between start and end dates to split the record in half
# data1 = data frame of pre-disturbance data (formatted dates in 1st column, numerical values in 2nd column)
# data2 = data frame of post-disturbance data (formatted dates in 1st column, numerical values in 2nd column)
# Start_Month = first month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# End_Month = last month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# min_quantile = minimum value for the FDC quantile range
# max_quantile = maximum value for the FDC quantile range
# increment = step size for the range of quantiles
# p = cdfs at which to estimate percentiles
# pre_data = list of pre-disturbance data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column)
# post_data = list of post-disturbance data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column)
# logsp_rlsp: T/F flag where when T = performs MMMT on log transformed data and F = performs MMMT on data in real space
# MMMT_out = matrix of output from Modified Mood's Median Test at each percentile (row) with columns of p-value, test statistic, and p-value flag
# label = array of labels for plot using the original plot method where [1] = plot title, [2] = x-axis title, [3] = y-axis title, [4] = data set 1 (pre-disturbance), and [5] = data set 2 (post-disturbance)
# yr1_1/yr1_2 = year of the first day of data1/data2
# yrn_1/yrn_2 = year of the last day of data1/data2
# ex = T/F flag where T = x-axis plots exceedance probabilities (1-cdf) and F = axis plots cdfs
# logy = T/F flag where when T = y-axis is plotted on a log scale and F = y-axis is plotted in real space
# FDR_flag = result of the field/global significance using a multiple False Discovery Rate (FDR) Test
# FDR_In = T/F flag where when T = FDR test results are plotted on figure and F = FDR test results are not plotted on figure
# pvalueflag_In = T/F flag where when T = Median Annual Percent zero are plotted on figure and F = Median Annual Percent zero are not plotted on figure
# yaxis_custom = T/F flag where when T = user can customize the step size of the y-axis, F= keeps the original axis logic exactly as is in the plot function
# ymin = minimum value on y axis for the plot using the original plot method
# ymax = maximum value on y axis for the plot using the original plot method
# pdfname = name of pdf file that the plot using the original plot method should be saved as

###
# Load packages and external functions
###

  # Load required packages 

  library(Hmisc)
  library(dataRetrieval)
  library(tidyr)
  library(dplyr)
  library(lubridate)

  # Source functions in MMMT.R

  source("MMMT.R")

###
# Data Import
###

  # Read in the water quality data (TN) .csv file downloaded from Chesapeake Bay Non-tidal Network website. 
  # While additional information is contained in the example input file, 
  # only the date and total nitrogen concentration (mg/l) for each sample is used in this analysis.
  
  data <- read.csv("Waterquality_data.csv")

  # Date column is in character format. Convert this to a date format
  
  data$SampleDate <- as.Date(data$SampleDate, format = "%Y-%m-%d")
 
  # Enter the start and end dates of the record
  
  start_date <- data[1,9]
  end_date <- data[dim(data)[1],9]
  start_year <- as.numeric(format(start_date, "%Y"))
  end_year <- as.numeric(format(end_date, "%Y"))

  # Find the difference between the start and end period to split the record into halves

  Difference_in_years <- ceiling(((end_year - start_year)/2))

  # Extract the dates and total nitrogen concentration measurements   

  data <- as.data.frame(data[,c(9,20)])
  data$Year <- format(data[,1], "%Y")
  colnames(data) <- c("Date", "Total Nitrogen Concentration", "Year")

  # Split  the primary dataset into two distinct dataframes representing the 
  # pre-disturbance/early and post-disturbance/late periods by specifying the start and 
  # end dates for each period 

  data1 <- data[data$Year >= start_year & data$Year <= (start_year+(Difference_in_years-1)), ]
  data2 <- data[data$Year >= (start_year+Difference_in_years) & data$Year <= end_year, ]
   
###
# Preprocess pre- and post-disturbance data into a matrix with all
# data in a year or specific time period in each column
###
   
  # Enter the Start month and end month
  
  Start_Month<-1
  End_Month<-12
   
  # User must specify the range of FDC percentiles to estimate based on 1/(nmin+1) to nmin/(nmin+1)
  # where nmin is the minimum number of samples in any one year

  min_quantile <- 0.05
  max_quantile <- 0.95
  increment <- 0.05
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

  # Run MMMT on each row of pre_data and post_data

  MMMT_out<-MMMT(pre_data,post_data, logsp_rlsp)
  MMMT_out <- as.data.frame(MMMT_out)
   
  # Calculate field/global significance using a multiple False Discovery Rate (FDR) Test
  # [3 = p-value < 0.01 ; 2 = 0.01 < p-value < 0.05 ; 1 = 0.05 < p-value < 0.10 ; 0 = p-value > 0.10] ]

  FDR_flag <- FDR(MMMT_out)
   
###
# Plotting
###
# Figure 4  
###
  
  # Set the name of the pdf to be produced
  
  pdfname <-paste0("Figure_4.pdf") 

  # Set labels for plot  
  
  label<-as.character(array(0,5))
  label[1]<- paste0("Conestoga River at Conestoga, Pa")
  label[2]<- paste0("Total Nitrogen (mg/l)")
  label[3]<-"Non-Exceedance Probability (cdf)"

  # Enter the start and end years for the pre-disturbance/early and post-disturbance/late records
  
  yr1_1<-as.POSIXlt(data1[1,1])$year+1900
  yrn_1<-as.POSIXlt(data1[dim(data1)[1],1])$year+1900
  yr1_2<-as.POSIXlt(data2[1,1])$year+1900
  yrn_2<-as.POSIXlt(data2[dim(data2)[1],1])$year+1900
  
  label[4]<-paste("Early Record ",yr1_1,"-",yrn_1,sep="")
  label[5]<-paste("Late Record ",yr1_2,"-",yrn_2,sep="")

  # T/F flag where T = exceedance probabilities and F = non-exceedance probabilities

  ex<-FALSE

  # T/F flag where T = y-axis on a log scale and F = y-axis in real-space

  logy<-FALSE

  # T/F flag where T = FDR test results plotted on figure and F = FDR test results not plotted on figure

  FDR_In<-TRUE

  # T/F flag where T = MMMT test results for median annual percentage of zeros plotted on 
  # figure and F = not plotted on figure

  pvalueflag_In<-FALSE

  # T/F flag where T = customize the stepsize of the y-axis 
  # figure and F = keep the y-axis logic as is in the function

  yaxis_custom <- TRUE
   
  # Set ymin and ymax

  ymin<-floor(min(min(data1[,2], na.rm = TRUE), min(data2[,2], na.rm = TRUE)))
  ymax<-ceiling(max(max(data1[,2], na.rm = T), max(data2[,2], na.rm = T)))
   
  # Customize the y-axis step size
  
  y_steps <- c(2) 
    
  # Create plot and export as a pdf (name of pdf defined as pdfname) 
  
  Plot_Func(pre_data$Percentiles,post_data$Percentiles,MMMT_out,p,label,logy,ymin,ymax,ex,FDR_In,FDR_flag,pvalueflag_In, yaxis_custom,y_steps, pdfname)
  
  