# Priyanka Rajashekar and Chuck Kroll
# SUNY College of Environmental Science and Forestry (ESF), Syracuse, NY, USA
# Contact: prajashe@syr.edu and cnkroll@esf.edu
# Modified Mood's Median Test on DC Percentiles
# Figure 3a: Annual Duration Curves (DCs) from depth to groundwater records
# Last Updated: August 2025

# This example produces Figure 3a from the recently submitted article: 
# Rajashekar, P, Kroll, C.N., and Vogel, R.M. (2025) 
# Extending the MMMT temporal change detection tool to assess significant shifts 
# in diverse environmental variables, submitted to JAWRA August 2025.

# For this 3rd case study, depth to groundwater data are obtained from the USGS NWIS database at two regions namely
# Lexington, NE :USGS#404717099460501. The pre-disturbance record was chosen as the 1986-2005 and
# the post-disturbance record was chosen as the 2006-2024.

# To perform this analysis at another USGS groundwater site, the user would need to change
# information in the "Data Import" (line 75) and "Plotting" (line 175) sections of this script.

# Variable Description
# SiteNumber = Station site identifiers
# state_code = State where the station is located
# whats_available: data types (parm_cd) are available for the site
# Site_Info: Site information for daily values of groundwater measurements
# Start_date = first date in the record
# End_date = last date in the record
# Start_year = first year in the record
# End_year = last year in the record
# Difference_in_years = difference between start and end dates to split the record in half
# Dataset_list = Dataset for each site stored as a list
# Groundwater_Data = Groundwater data extracted for the site
# letters_suffix = Sequence of letters for naming figure PDF files
# y_steps = Step size to customize the y-axis 
# data1 = data frame of pre-disturbance data 
# data2 = data frame of post-disturbance data 
# Start_Month = first month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# End_Month = last month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# min_quantile = minimum value for the FDC quantile range
# max_quantile = maximum value for the FDC quantile range
# increment = step size for the range of quantiles
# p = cdfs at which to estimate percentiles
# pre_data = list of pre-disturbance data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column)
# post_data = list of post-disturbance data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column)
# logsp_rlsp = flag where when T = Performs MMMT on log transformed dataset and F = Performs MMMT on dataset in realspace 
# MMMT_out = matrix of output from Modified Mood's Median Test at each percentile (row) with columns of p-value, test statistic, and p-value flag
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

  # Enter the site number for LExington, NE
  # [For different USGS sites, the site number will change]

  SiteNumber <- c("404717099460501")

  # Enter the state in which the site is located

  state_code <- "NE"

  # Check the types of data (parm_cd) available at the site

  whats_available <- whatNWISdata(siteNumber = SiteNumber)

  # Extract the site info for daily values of depth to groundwater (ft; parm_cd = 72019)
  
  Site_Info <- whats_available %>%
  filter(data_type_cd == "dv")%>%
  filter(parm_cd %in% c("72019"))

  # Enter the start and end dates of the record
  # [Can change date ranges here]

  start_date <- as.Date("1986-01-01")
  end_date <- as.Date("2024-12-31")
  start_year <- as.numeric(format(start_date, "%Y"))
  end_year <- as.numeric(format(end_date, "%Y"))

  # Find the difference between the start and end period to split the record into halves
  
  Difference_in_years <- ceiling(((end_year - start_year)/2))

  # Download the data through the USGS portal using DataRetrieval Package for specific time periods
  
  Groundwater_Data <- readNWISdv(Site_Info[,2], parameterCd = Site_Info[,14], 
                                 statCd = Site_Info[,15], 
                                 startDate = start_date, endDate = end_date)

  # Extract the depth to groundwater and dates for each site  
  
  Groundwater_Data <- as.data.frame(Groundwater_Data[,3:4])
  Groundwater_Data$Year <- format(Groundwater_Data[,1], "%Y")
  colnames(Groundwater_Data) <- c("Date", "GW_Level", "Year")
  
  #Split  the primary dataset into two distinct dataframes representing the 
  #pre-disturbance/early and post-disturbance/late periods by specifying the start and 
  # end dates for each period 

  data1 <- Groundwater_Data[Groundwater_Data$Year >= start_year & Groundwater_Data$Year <= (start_year+Difference_in_years), ]
  data2 <- Groundwater_Data[Groundwater_Data$Year >= (start_year+Difference_in_years+1) & Groundwater_Data$Year <= end_year, ]
  
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
  
  MMMT_out<-MMMT(pre_data,post_data, logsp_rlsp)
  MMMT_out <- as.data.frame(MMMT_out)
  
  # Calculate field/global significance using a multiple False Discovery Rate (FDR) Test
  # [3 = p-value < 0.01 ; 2 = 0.01 < p-value < 0.05 ; 1 = 0.05 < p-value < 0.10 ; 0 = p-value > 0.10] ]

  FDR_flag <- FDR(MMMT_out)
  
###
# Plotting
###
# Figure 3a  
###
  
  # Set the name of the pdf to be produced
  
  pdfname <-paste0("Figure_3a.pdf") 

  # Set labels for plot, determining the beginning and ending year for both data sets  
  
  label<-as.character(array(0,5))
  label[1]<- paste0("Lexington, NE")
  label[2] <- paste0("Depth to Groundwater (ft)")
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

  # set ymin and ymax

  ymin<-floor(min(Groundwater_Data$GW_Level, na.rm = T))
  ymax<-ceiling(max(Groundwater_Data$GW_Level, na.rm = T))

  # Customize the y-axis step size

  y_steps <- c(2) 
  
  
  # Create plot and export as a pdf (name of pdf defined as pdfname) 
  
  Plot_Func(pre_data$Percentiles,post_data$Percentiles,MMMT_out,p,label,logy,ymin,ymax,ex,FDR_In,FDR_flag,pvalueflag_In,yaxis_custom,y_steps, pdfname)


