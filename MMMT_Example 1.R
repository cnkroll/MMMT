# Chuck Kroll
# February 2022
# Modified Mood's Median Test on FDC Percentiles
# Example #1: Annual Flow Duration Curves (FDCs) from daily streamflow records

# In this example, daily streamflow are obtained at the West Fork of the Tualatin River near Dilley, OR, USA
# USGS#14203500. A dam was installed on this river in 1971. The pre-disturbance record was chosen as the 1940-1960 water years
# and the post-disturbance record was chosen as the 1975-2020 water years. For methodology, see draft paper at
# www.esf.edu/ere/kroll/publications/, Mowafy, Kroll and Vogel, A Graphical Non-Parameteric Hydrologic Alteration Test Using Flow Duration Curves

# Variable Description
# site = USGS #, as a character string
# data1a = daily streamflow downloaded from USGS site for pre-disturbance record
# data2a = daily streamflow downloaded from USGS site for post-disturbance record
# data1 = data frame of pre-disturbance data (dates formatted as.Date in 1st column, numerical values in 2nd column)
# data2 = data frame of post-disturbance data (dates formatted as.Date in 1st column, numerical values in 2nd column)
# Start_Month = first month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# End_Month = last month of analysis (used to partition data sets) where 1 = January, 2 = February, . . . , 12 = December
# p = cdfs at which to estimate percentiles
# pre_data = list of pre-disturbance data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column)
# post_data = list of post-disturbance data output which contains: [[1]] annual percentile and [[2]] annual raw data for each year (column)
# MMMT_out = matrix of output from Modified Mood's Median Test at each percentile (row) with columns of p-value, test statistic, and p-value flag
# label = array of labels for plot where [1] = plot table, [2] = x-axis title, [3] = y-axis title, [4] = data set 1 (pre-disturbance), and [5] = data set 2 (post-disturbance)
# yr1_1/yr1_2 = year of the first day of data1/data2
# yrn_1/yrn_2 = year of the last day of data1/data2
# ex = T/F flag where T = x-axis plots exceedance probabilities (1-cdf) and F = axis plots cdfs
# logy = T/F flag where when T = y-axis is plotted on a log scale and F = y-axis is not plotted on a log scale
# ymin = minimum value on y axis
# ymax = maximum value on y axis
# pdfname = name of pdf file that the plot should be saved as

# Source functions in MMMT.R

source("MMMT.R")

# Load required package Hmisc

library(Hmisc)

# Load USGS dataRetrieval package to obtain streamflow data

library("dataRetrieval")

# Select site: USGS#14203500 Tualatin River near Dilley, OR, USA
# A dam was installed on this river in 1971

site<-c("14203500") 

data1a<-readNWISdv(site[1],"00060","1939-10-01","1965-9-30")
data1<-data1a[,3:4]
data2a<-readNWISdv(site[1],"00060","1975-10-01","2020-9-30")
data2<-data2a[,3:4]
info<-readNWISsite(siteNumber = site[1])

###
# Preprocess pre- and post-disturbance data into a matrix with all
# data in a year or specific time period (e.g., water year) in each column
###

# Start month and day

Start_Month<-10
End_Month<-9

# Set FDC percentiles to estimate 

p<-seq(0.005,.995,.005)

# Call function to preprocess pre-disturbance and post-disturbance data
# (data between Start_Month and End_Month in each year in columns) and 
# determine FDC percentiles for each year

pre_data<-prepro(data1,Start_Month,End_Month,p)
post_data<-prepro(data2,Start_Month,End_Month,p)

###
# MMMT Routine
###

# Run MMMT on each row of pre_data and post_data

MMMT_out<-MMMT(pre_data$Percentiles,post_data$Percentiles)

###
# Plotting Routine
###

# Set labels for plot, determining the beginning and ending year for both data sets

label<-as.character(array(0,5))
label[1]<-paste(info[3],": USGS Site #",site[1],sep="")
label[2]<-"Streamflow (cfs)"
label[3]<-"Exceedance Probability"
yr1_1<-as.POSIXlt(data1[1,1])$year+1900
yrn_1<-as.POSIXlt(data1[dim(data1)[1],1])$year+1900
yr1_2<-as.POSIXlt(data2[1,1])$year+1900
yrn_2<-as.POSIXlt(data2[dim(data2)[1],1])$year+1900
label[4]<-paste("Pre-disturbance ",yr1_1,"-",yrn_1,sep="")
label[5]<-paste("Post-disturbance ",yr1_2,"-",yrn_2,sep="")
ex<-TRUE
logy<-TRUE
ymin=0.1
ymax=10000
pdfname<-"Fig Ex 2.pdf"


# Create plot and export as a pdf (name of pdf defined as pdfname) 

Plot_Func(pre_data$Percentiles,post_data$Percentiles,MMMT_out,p,label,logy,ymin,ymax,ex,pdfname)


