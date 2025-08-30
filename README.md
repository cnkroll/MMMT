Name: Modified Mood’s Median Test (MMMT) Authors: Priyanka Rajashekar (prajashe@syr.edu), Charles N. Kroll (cnkroll@esf.edu) and Mohamed Mowafy (mowafymd@mail.uc.edu) 
Date: August 2025

Introduction:

This group of six R functions is designed to perform a Modification of Mood’s Median Test (MMMT) across a series of percentiles. It’s ideally suited for situations where one has a time series over many years of record, both before (pre-disturbance/early) and after (post-disturbance/late) some potential disturbance to the system. Initially developed for hydrologic Flow Duration Curves based on U.S. water years, the tool has been successfully applied to a variety of other contexts for quantifying system changes over time (e.g., intermittent streamflows, precipitation, groundwater depth, water quality, air temperature).
The 6 functions are:
  1. prepro = a preprocessing function that partitions raw data between user-defined start and end months on an annual basis, calculates percentiles for each year, and determines the percentage of zero values per year when zeros are present in the dataset.
  2. MMMT = a loop to perform Modified Mood's Median Test (MMMT) across a range of percentiles and median annual percentage of zeros when zeros are present in the dataset.
  3. MMMT_once = calculates MMMT test statistic, p-values and flags for p-value ranges for a single sample.
  4. FDR = calculates field/global significance using a multiple False Discovery Rate (FDR) Test
  5. Plot_Func = a routine to plot MMMT output across a range of percentiles. 
  6. Plot_Func2 =  a routine to plot quantile ratios of pre and post medians across a range of percentiles. 

Each function is described below.

Description of Functions:

  1. prepro(data, Start_Month, End_Month, p)

data = a data frame of 2 columns: (1) the date of each data point formatted as a DATE in R and (2) raw data for each date formatted as a numeric value. 

Start_Month = Month as an integer (1 = Jan, 2 = Feb, . . . 12 = Dec) in which to start analysis determined by the user. The annual analysis begins on the first day of the Start_Month. For a hydrologic water year in the US, which begins in October and ends in September, Start_Month = 10. 

End_Month = Month as an integer (1 = Jan, 2 = Feb, . . . 12 = Dec) in which to end analysis determined by the user. The annual analysis ends on the last day of the End_Month. For a hydrologic water year in the US, which begins in October and ends in September, End_Month = 9. 

p = array of probabilities (cdfs) as numeric values determined by the user at which percentiles are calculated using the Harrell and Davis (1982) distribution-free quantile estimator for data in each year between the Start_Month and End_Month. Note that the smallest p should be at least 1/(nmin+1) and the largest p should be no bigger than 1 – 1/(nmin+1), where nmin is the minimum number of data points in any one year.

The prepro function creates a list with two objects for non-zero datasets (Percentiles and Data) and three objects for datasets with zeros (Percentiles, Data, and Zero data):

Percentiles: a matrix of the calculated percentiles at each value in array “p” for each year of the record (between Start_Month and End_Month). Each column is a separate year.

Data: a matrix of the raw data used to calculate the Percentiles for each year of the record (between Start_Month and End_Month). Each year’s data is in the order that it appears in original data set (from the first data point of Start_Month to the last data point of End_Month). Each column is a separate year.

Zero Data: percentage of zero values in each year

The prepro function should be run on both the pre- and post-disturbance record as the outputted Percentiles matrices are employed in the function MMMT (and MMMT_once). The values of p should be the same when processing the pre- and post-disturbance records. There is a limit in any one year of 10000 data points for either record (this could be modified on line 86 of R script). NOTE: The user should be careful with the start and end date of their raw data and Start_Month and End_Month of the analysis so not to have a partial year of data in the first or last year of the record.

  2. MMMT(pre, post, logsp_rlsp)

pre = preprocessed data from running the prepro on the pre-disturbance/early record.

post = preprocessed data from running the prepro on the post-disturbance/late record.

logsp_rlsp = T/F flag where when T = performs MMMT on log transformed data, F = performs MMMT on data in real space

For non-zero datasets:
The MMMT function applies the Modified Mood’s Median Test (Fligner and Rust, 1982) to compare pre-disturbance and post-disturbance percentiles for each year. Its output is a list which contains one item, a matrix with three columns: the first contains the test statistic, the second gives the corresponding p-value, and the third provides a significance flag based on the p-value (3 = p < 0.01; 2 = 0.01 < p < 0.05; 1 = 0.05 < p < 0.10; 0 = p > 0.10). For each percentile, the function MMMT_once is executed. 

For datasets with zero values:
The MMMT function performs the Modified Mood’s Median Test (Fligner and Rust, 1982) to compare pre-disturbance and post-disturbance percentiles for each year, as well as the percentage of zero values between the two periods. The function returns a list of five objects:
i. A matrix containing the MMMT results for the percentiles, where the first column is the test statistic, the second column is the corresponding p-value, and the third column is a categorical significance flag (3 = p < 0.01; 2 = 0.01 < p < 0.05; 1 = 0.05 < p < 0.10; 0 = p > 0.10). For each percentile, the function MMMT_once is executed.
ii. The MMMT results for the annual percentage of zero values, including the test statistic, its p-value, and the corresponding significance flag.
iii. The median percentage of zero values during the pre-disturbance (early record) period.
iv. The median percentage of zero values during the post-disturbance (late record) period.
v. The significance flag associated with the p-value for the zero-value comparison.

  3. MMMT_once(F1, F2,  logsp_rlsp)

F1 = array of pre-disturbance values 

F2 = array of post-disturbance values 

logsp_rlsp = T/F flag where when T = performs MMMT on log transformed data, F = performs MMMT on data in real space
The MMMT_once function conducts the Modified Mood’s Median Test (MMMT) on two datasets, F1 and F2, and returns the corresponding test statistic and p-value. Its notation follows that of Mowafy et al. (2023). Although MMMT_once is called internally by the MMMT function, it can also be used independently to compute the MMMT test statistic and p-value for any pair of samples. The cumulative distribution functions (CDFs) within MMMT_once can be calculated either in logarithmic or real space, depending on the user-specified flag.

  4. FDR(MMMT_out)

MMMT_out = matrix of output from MMMT, including the first column of p-values for MMMT test at each percentile.

The FDR function performs test of field/global significance using results of the MMMT test across all percentiles using a multiple False Discovery Rate (FDR) Test (see Mowafy et al., 2023). The result is a flag which equals 3 if the p-value of the FDR test < 0.01, 2 if 0.01 < p-value of FDR test < 0.05, 1 if 0.05 < p-value of FDR test < 0.10, and 0 if p-value of FDR test > 0.10. Including this result is an option in the output figure.

  5. Plot_func(pre,post,mood_modified,p,label,logy,ymin,ymax,ex,FDR_In,FDR_flag, pvalueflag_In, yaxis_custom, y_step, pdfname)

pre = output matrix of annual percentiles from running prepro on pre-disturbance record. 

post = output matrix of annual percentiles from running prepro on post-disturbance record. 

mood_modified = output matrix from running MMMT using pre and post. 

p = array of cumulative distribution function values at which percentiles are calculated using the Harrell and Davis (1982) Distribution-Free Quantile Estimator for data in each year between the Start_Month and End_Month. 

label = array of 5 characters which contains: label[1] = Title of plot (placed on top of figure) label[2] = y-axis label label[3] = x-axis label label[4] = Name/description of pre-disturbance data set (for legend) label[5] = Name/description of post-disturbance data set (for legend) Note that the user may have to adjust label lengths to properly fit on figure. 

logy = logical character that is either TRUE or FALSE. When TRUE, the y-axis is plot on a log base-10 scale; when FALSE, the y-axis is plot on a linear scale. 

ymin = minimum value on y-axis

ymax = maximum value on y-axis 

ex = logical character that is either TRUE or FALSE. When TRUE, the exceedance probabilities for each percentile will be plotted on the x-axis; when FALSE, the cdfs for each percentile (1 – the exceedance probabilities) will be plotted on the x-axis.  

FDR_In = logical character that is either TRUE or FALSE. When TRUE, the result of the FDR test is included on the output figure; when FALSE, the result of the FDR test is not included. 

FDR_flag = output from the FDR function 

pvalueflag_In = logical character that is either TRUE or FALSE. When TRUE, the median annual percentages of zeros for both periods are included on the output figure; when FALSE, the values are not displayed on the figure. 

yaxis_custom = logical character that is either TRUE or FALSE. When TRUE, skips this line entirely so the user can add a custom axis in the main script using axis(); When FALSE, keeps the original axis logic exactly as is in the function 

y_steps = step size for the y-axis values

pdfname = character string with name of the pdf file to be created by Plot_func.

The Plot_func function creates a PDF showing a plot of the MMMT function output. Percentiles for each year are plotted from the pre- and post-disturbance matrices as light lines, and the medians at each percentile are shown as darker symbols. The MMMT p-values are used to color-code the medians for each percentile: red for p < 0.01, orange for 0.01 ≤ p < 0.05, yellow for 0.05 ≤ p < 0.10, and blue for p > 0.10. The same color scheme is applied to indicate the results of the FDR test (via the “FDR p-value” text) when FDR = TRUE, and to represent the MMMT test results on the median annual percentage of zeros when pvalueflag_In = TRUE.

  6. Plot_func2(pre,post,mood_modified,p,label,ymin2,ymax2,ex,FDR_In,FDR_flag,yaxis_custom,y_step,pdfname)

pre = output matrix of annual percentiles from running prepro on pre-disturbance record. 

post = output matrix of annual percentiles from running prepro on post-disturbance record. 

mood_modified = output matrix from running MMMT using pre and post. 

p = array of cumulative distribution function values at which percentiles are calculated using the Harrell and Davis (1982) Distribution-Free Quantile Estimator for data in each year between the Start_Month and End_Month. 

label = array of 5 characters which contains: label[1] = Title of plot (placed on top of figure) label[2] = y-axis label label[3] = x-axis label label[4] = Name/description of pre-disturbance data set (for legend) label[5] = Name/description of post-disturbance data set (for legend) Note that the user may have to adjust label lengths to properly fit on figure. 

ymin = minimum value on y-axis

ymax = maximum value on y-axis 

ex = logical character that is either TRUE or FALSE. When TRUE, the exceedance probabilities for each percentile will be plotted on the x-axis; when FALSE, the cdfs for each percentile (1 – the exceedance probabilities) will be plotted on the x-axis. 

FDR_In = logical character that is either TRUE or FALSE. When TRUE, the result of the FDR test is included on the output figure; when FALSE, the result of the FDR test is not included. 

FDR_flag = output from the FDR function  
yaxis_custom = logical character that is either TRUE or FALSE. When TRUE, skips this line entirely so the user can add a custom axis in the main script using axis(); When FALSE, keeps the original axis logic exactly as is in the function 

y_steps = step size for the y-axis values

pdfname = character string with name of the pdf file to be created by Plot_func.

The Plot_func2 function generates a PDF of a plot showing quantile ratios versus exceedance/non-exceedance probabilities. Quantile ratios for each year are calculated by dividing the pre- and post-disturbance matrices by the median of the pre-disturbance (early) period and are shown as light lines on the plot. The post-disturbance medians at each percentile are represented as the ratio of post median to pre median (darker symbols), while the pre-disturbance medians are represented as 1, since they are divided by themselves. The MMMT p-values are used to color-code the median quantile ratios for each percentile: red for p < 0.01, orange for 0.01 ≤ p < 0.05, yellow for 0.05 ≤ p < 0.10, and blue for p > 0.10.
References:
Fligner, M.A., Rust, S.W. 1982. A modification of Mood’s median test for the generalized Behrens-Fisher problem. Biometrika, 69, 221–226, doi.org/10.1093/biomet/69.1.221.

Harrell, F.E., and Davis, C.E. 1982. A new distribution-free quantile estimator. Biometrika, 69, 635-640, doi.org/10.1093/biomet/69.3.635.

Mowafy, M.H.M.A., Kroll, C.N., and Vogel, R.M. 2023. A graphical non-parametric hydrologic alteration test using flow duration curves. Hydrological Science Journal, doi.org/10.1080/02626667.2023.2248109.

Rajashekar, P., Kroll, C.N., and Vogel, R.M. (2025). Extending the MMMT temporal change detection tool to assess significant shifts in diverse environmental variables, submitted to JAWRA August 2025.
