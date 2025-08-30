# Priyanka Rajashekar and Chuck Kroll
# SUNY College of Environmental Science and Forestry (ESF), Syracuse, NY, USA
# Contact: prajashe@syr.edu and cnkroll@esf.edu
# Modified Mood's Median Test on FDC Percentiles
# FIgure 3: Annual Duration Curves (DCs) from daily precipitation records
# Last Updated: August 2025

# This example produces Figure 2a from the recently submitted article: 
# Rajashekar, P, Kroll, C.N., and Vogel, R.M. (2025) 
# Extending the MMMT temporal change detection tool to assess significant shifts 
# in diverse environmental variables, submitted to JAWRA August 2025.

# For this 2nd case study, daily precipitation data from the Global Historical Climatology Network Daily 
# (GHCN-D; https://www.ncei.noaa.gov/cdo-web/search?datasetid=GHCND), compiled by 
# NOAA’s National Centers for Environmental Information for Riverside, CA (USC00047470) are obtained for the period 1901 to 2023.

# To perform this analysis at another site, the user would need to change information in 
# the "Data Import" (line 68) and "Plotting" (line 128) sections of this script. In this 
# example, data from NOAA NCEI was downloaded and saved as a csv file, which was 
# read into R. See "Precipitation_Riverside.csv" for file structure and column headings.

# Variable Description
# data = data for the site selected for the entire period of record
# data1 = data frame of pre-disturbance data 
# data2 = data frame of post-disturbance data 
# Start_Month = first month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# End_Month = last month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# min_quantile = minimum value for the FDC quantile range
# max_quantile = maximum value for the FDC quantile range
# increment = step size for the range of quantiles
# p = cdfs at which to estimate percentiles
# pre_data = list of pre-disturbance data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column) [[3]] percentage of zeroes in each year
# post_data = list of post-disturbance data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column) [[3]] percentage of zeroes in each year
# logsp_rlsp = flag where when T = Performs MMMT on log transformed dataset and F = Performs MMMT on dataset in realspace 
# MMMT_out = list of [[1]] output from Modified Mood's Median Test at each percentile (row) with columns of p-value, test statistic, and p-value flag, 
# [[2]] MMMT output with zero data [[3]] median percentage of zeros in the pre-disturbance/early period
# [[4]] median percentage of zeros in the post-disturbance/late period [[5]] pvalue flag for the zero data
# label = array of labels for plot where [1] = plot table, [2] = x-axis title, [3] = y-axis title, [4] = data set 1 (pre-disturbance), and [5] = data set 2 (post-disturbance)
# yr1_1/yr1_2 = year of the first day of data1/data2
# yrn_1/yrn_2 = year of the last day of data1/data2
# ex = T/F flag where T = x-axis plots exceedance probabilities (1-cdf) and F = axis plots cdfs
# logy = T/F flag where when T = y-axis is plotted on a log scale and F = y-axis is not plotted on a log scale
# FDR_flag = result of the field/global significance using a multiple False Discovery Rate (FDR) Test
# FDR_In = T/F flag where when T = FDR test results are plotted on figure and F = FDR test results are not plotted on figure
# pvalueflag_In = T/F flag where when T = Median Annual Percent zero are plotted on figure and F = Median Annual Percent zero are not plotted on figure
# yaxis_custom = T/F flag where when T = user can customize the step size of the y-axis, F= keeps the original axis logic exactly as is in the plot function
# ymin = minimum value on y axis
# ymax = maximum value on y axis
# y_steps = Step size to customize the y-axis
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
  
  # Read in the .csv file containing the precipitation (in) data for the study site Riverside, CA (USC00047470)

  data <- read.csv(file="Precipitation_Riverside.csv", header = T, sep = ",")   

  # format the date column
  
  data$DATE <- as.Date(data$DATE, format = "%m/%d/%Y")
  
  #Split  the primary dataset into two distinct dataframes representing the 
  #pre-disturbance/early and post-disturbance/late periods by specifying the start and 
  # end dates for each period   
  
  data1 <- data[data$DATE >= "1901-10-01" & data$DATE <= "1960-09-30", ]
  data2 <- data[data$DATE >= "1961-10-01" & data$DATE <= "2023-09-30", ]

###
# Preprocess pre-disturbance/early and post-disturbance/late data into a matrix with all
# data in a year or specific time period in each column
###

  # Start month and day

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

  # flag where when T = performs MMMT on log transformed data and F = performs MMMT on data in real space

  logsp_rlsp <- FALSE

  # Run MMMT on each row of pre_data and post_data

  MMMT_out <- MMMT(pre_data,post_data,logsp_rlsp)

  # Calculate field/global significance using a multiple False Discovery Rate (FDR) Test
  # [3 = p-value < 0.01 ; 2 = 0.01 < p-value < 0.05 ; 1 = 0.05 < p-value < 0.10 ; 0 = p-value > 0.10] ]

  FDR_flag <- FDR(MMMT_out[[1]])

###
# Plotting
###
# Figure 2a  
###
  
  # Set the name of the pdf to be produced
  
  pdfname <-paste0("Figure_2a.pdf") 
  
  # Set labels for plot  
  
  label<-as.character(array(0,5))
  label[1]<- "Precipitation in Riverside, CA"
  label[2]<-"Precipitation (inches)"
  label[3]<-"Non-Exceedance Probability (Cdf)"
  
  # Enter the start and end years for the pre-disturbance/early and post-disturbance/late records
  
  yr1_1<-as.POSIXlt(data1[1,1])$year+1900
  yrn_1<-as.POSIXlt(data1[dim(data1)[1],1])$year+1900
  yr1_2<-as.POSIXlt(data2[1,1])$year+1900
  yrn_2<-as.POSIXlt(data2[dim(data2)[1],1])$year+1900

  label[4]<-paste("Early-record ",yr1_1,"-",yrn_1,sep="")
  label[5]<-paste("Late-record ",yr1_2,"-",yrn_2,sep="")

  # T/F flag where T = exceedance probabilities and F = non-exceedance probabilities
  
  ex<-FALSE

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

  # set ymin and ymax

  ymin<-floor(min(data$PRCP, na.rm = T))
  ymax<-1.0

  # Customize the y-axis step size

  y_steps <- c(0.25)  

  
  # Create plot and export as a pdf (name of pdf defined as pdfname) 
  
  Plot_Func(pre_data$Percentiles,post_data$Percentiles,MMMT_out[[1]],p,label,logy,ymin,ymax,ex,FDR_In,FDR_flag,pvalueflag_In, yaxis_custom, y_steps, pdfname)




