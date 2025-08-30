# Priyanka Rajashekar and Chuck Kroll
# SUNY College of Environmental Science and Forestry (ESF), Syracuse, NY, USA
# Contact: prajashe@syr.edu and cnkroll@esf.edu
# Modified Mood's Median Test on DC Percentiles
# Figure 5: Annual Duration Curves (DCs) from daily maximum temperature records
# Last Updated: August 2025

# This example produces Figures 5a and 5b from the recently submitted article: 
# Rajashekar, P, Kroll, C.N., and Vogel, R.M. (2025) 
# Extending the MMMT temporal change detection tool to assess significant shifts 
# in diverse environmental variables, submitted to JAWRA August 2025.

# For this 5th case study, maximum daily temperature data are obtained from the NOAA National Centers for Environmental 
# Information (https://www.ncei.noaa.gov/cdo-web/search?datasetid=GHCND). The pre-disturbance record was chosen as the 1990-2006 and
# the post-disturbance record was chosen as the 2007-2023.

# To perform this analysis at another site, the user would need to change information in 
# the "Data Import" (line 70) and "Plotting" (line 137) sections of this script. In this 
# example, data from the NOAA NCEI was downloaded and saved as a csv file, which was 
# read into R. See "Temperature_Data.csv" for file structure and column headings.

# Variable Description:
# data = data for the site selected for the entire period of record
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
# ymin2 = minimum value on y axis for the plot using the quantile ratio method  
# ymax2 = maximum value on y axis for the plot using the quantile ratio method
# label2 = array of labels for plot using the quantile ratio method where [1] = plot title, [2] = x-axis title, [3] = y-axis title, [4] = data set 1 (pre-disturbance), and [5] = data set 2 (post-disturbance)
# pdfname2 = name of pdf file that the plot using the quantile ratio method should be saved as

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

  # Read in the temperature data downloaded from NOAA NCEI for Phoenix Airport, AZ.
  # The example input file (Temperature_Data.csv) contains the date (m/d/Y) and columns for daily average, maximum 
  # and minimum temperature (°F), and only daily maximum is used in this analysis.
  
  data <- read.csv(file="Temperature_Data.csv", header = T, sep = ",")

  # Format the dates using as.Date and extract only the columns for dates 
  # and the values of the variable of interest.

  data[,2] <- as.Date(data[,2], format = "%m/%d/%Y")
  data <- data[, c("DATE", "TMAX")]

  # Split  the primary dataset into two distinct dataframes representing the 
  # pre-disturbance/early and post-disturbance/late periods by specifying the start and 
  # end dates for each period 

  data1 <- data[data$DATE >= "1990-01-01" & data$DATE <= "2006-12-31", ]
  data2 <- data[data$DATE >= "2007-01-01" & data$DATE <= "2023-12-31", ]

###
# Preprocess pre-disturbance/early and post-disturbance/late data into a matrix with all
# data in a year or specific time period in each column
###

  # Enter the Start month and end month

  Start_Month<-1
  End_Month<-12

  # User must specify the range of FDC percentiles to estimate based on 1/(nmin+1) to nmin/(nmin+1)
  # where nmin is the minimum number of samples in any one year
  # for instance for a dataset with nmin = 365 points in an year, min can be equal to 0.005, 
  # max can be equal to 0.995 and increment = 0.005
  # Note: User can specify the quantiles as well

  min_quantile <- 0.005
  max_quantile <- 0.995
  increment <- 0.005
  p <- seq(min_quantile, max_quantile, increment)

  # Call function to preprocess pre-disturbance and post-disturbance data
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

  MMMT_out<-MMMT(pre_data,post_data,logsp_rlsp)
  MMMT_out <- as.data.frame(MMMT_out)

  # Calculate field/global significance using a multiple False Discovery Rate (FDR) Test
  # [3 = p-value < 0.01 ; 2 = 0.01 < p-value < 0.05 ; 1 = 0.05 < p-value < 0.10 ; 0 = p-value > 0.10] ]

  FDR_flag <- FDR(MMMT_out)

###
# Plotting
###
# Figure 5a  
###
  
  # Set the name of the pdf to be produced
  
  pdfname <-paste0("Figure_5a.pdf") 
  
  # Set labels for plot, determining the beginning and ending year for both data sets  
  
  label<-as.character(array(0,5))
  label[1]<-"Phoenix Airport, AZ, US"
  label[2]<-"Daily Maximum Temperature (°F)"
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

  ymin<-floor(min(min(data1$TMAX, na.rm = TRUE), min(data2$TMAX, na.rm = TRUE)))
  ymax<-ceiling(max(max(data1$TMAX, na.rm = T), max(data2$TMAX, na.rm = T)))

  # Rescale the ymin and ymax values

  ymin <- 40 
  ymax <- 130

  # Customize the y-axis step size
  
  y_steps <- c(20)
  
  # Create plot using the original plot function and export as a pdf (name of pdf defined as pdfname) 

  Plot_Func(pre_data$Percentiles,post_data$Percentiles,MMMT_out,p,label,logy,ymin,ymax,ex,FDR_In,FDR_flag,pvalueflag_In, yaxis_custom,y_steps, pdfname)

###
# Figure 5b  
###
  
  # Create plot using the new quantile ratio method. 
  # Note that commands for Figure 5a starting on line 148 should be run prior to running the following commands.
  
  # Set the name of the pdf to be produced
  
  pdfname2 <-paste0("Figure_5b.pdf") 
  
  # Calculate field/global significance using a multiple False Discovery Rate (FDR) Test
  # [3 = p-value < 0.01 ; 2 = 0.01 < p-value < 0.05 ; 1 = 0.05 < p-value < 0.10 ; 0 = p-value > 0.10] ]
  
  FDR_flag <- FDR(MMMT_out)
  
  # T/F flag where T = FDR test results plotted on figure and F = FDR test results not plotted on figure
  
  FDR_In<-TRUE
  
  # Set ymin and ymax

  ymin2<-0.9
  ymax2<-1.2  

  # Set labels for plot
  
  label2<-label
  label2[2]<-"Quantile Ratios of Daily Maximum Temperature"

  # Customize the y-axis step size
  
  y_steps <- c(0.10)
  
  # Create plot and export as a pdf (name of pdf defined as pdfname) 
  
  Plot_Func2(pre_data$Percentiles,post_data$Percentiles,MMMT_out[[1]],p,label2,ymin2,ymax2,ex,FDR_In,FDR_flag,yaxis_custom,y_steps,pdfname2)










