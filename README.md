Name: Modified Mood’s Median Test (MMMT)
Authors: Charles N. Kroll (cnkroll@esf.edu) and Mohamed Mowafy (mowafymd@mail.uc.edu)
Date: June 2023

This group of five R functions is designed to perform a Modification of Mood’s Median Test (MMMT) across a series of percentiles. It’s ideally suited for situations where one has a time series over many years of record, both before (pre-disturbance) and after (post-disturbance) some potential disturbance to the system. While its initial development was for hydrologic Flow Duration Curves using US water years, we believe these functions may be useful for many other applications to assess change in a system over time (e.g., biotic, climatic, hydrologic, etc.).

The 5 functions are:
1)	prepro = a preprocessing function which partitions raw data between the start and end months at an annual basis and calculates percentiles for each year of data.
2)	MMMT = a loop to perform Modified Mood's Median Test (MMMT) across a range of percentiles.
3)	MMMT_once = calculates MMMT test statistic and p-values and flags for p-value ranges for a single sample.
4)	FDR = calculates field/global significance using a multiple False Discovery Rate (FDR) Test
5)	Plot_Func = a routine to plot MMMT output across a range of percentiles.
Each function is described below.

1) prepro(data, Start_Month, End_Month, p)

data = a data frame of 2 columns: (1) the date of each data point formatted as a DATE in R and (2) raw data for each date formatted as a numeric value. 
Start_Month = Month as an integer (1 = Jan, 2 = Feb, . . . 12 = Dec) in which to start analysis determined by the user. The annual analysis begins on the first day of the Start_Month. For a hydrologic water year in the US, which begins in October and ends in September, Start_Month = 10.
End_Month = Month as an integer (1 = Jan, 2 = Feb, . . . 12 = Dec) in which to end analysis determined by the user. The annual analysis ends on the last day of the End_Month. For a hydrologic water year in the US, which begins in October and ends in September, End_Month = 9.
p = array of probabilities (cdfs) as numeric values determined by the user at which percentiles are calculated using the Harrell and Davis (1982) distribution-free quantile estimator for data in each year between the Start_Month and End_Month. Note that the smallest p should be at least 1/(nmin+1) and the largest p should be no bigger than 1 – 1/(nmin+1), where nmin is the minimum number of data points in any one year. 

The prepro function creates a list with two objects: 
1)	Percentiles: a matrix of the calculated percentiles at each value in array “p” for each year of the record (between Start_Month and End_Month). Each column is a separate year.
2)	Data: a matrix of the raw data used to calculate the Percentiles for each year of the record (between Start_Month and End_Month). Each year’s data is in the order that it appears in original data set (from the first data point of Start_Month to the last data point of End_Month). Each column is a separate year.
	
The prepro function should be run on both the pre- and post-disturbance record as the outputted Percentiles matrices are employed in the function MMMT (and MMMT_once). The values of p should be the same when processing the pre- and post-disturbance records. There is a limit in any one year of 10000 data points for either record (this could be modified on line 72 of R script).
NOTE: The user should be careful with the start and end date of their raw data and Start_Month and End_Month of the analysis so not to have a partial year of data in the first or last year of the record.

2) MMMT(pre, post) 

pre = output matrix of Percentiles from running prepro (variable “Percentile” of prepro output list) on pre-disturbance record.
post = output matrix of Percentiles from running prepro (variable “Percentile” of prepro output list) on post-disturbance record.

The MMMT function takes a matrix of pre-disturbance and post-disturbance percentiles for each year, and for each percentile (a specific row for each pre and post matrix) it performs the Modified Mood’s Median Test (Fligner and Rust 1982). The output from MMMT is a matrix where the first column is the test statistic from MMMT, the second column is the p-value for the test statistic from MMMT, and the third column is a flag for p-value (3 = p-value < 0.01; 2 = 0.01 < p-value < 0.05; 1 = 0.05 < p-value < 0.10, and 0 = p-value > 0.10). For each percentile (row of input matrices), MMMT_once is called.

3) MMMT_once(F1, F2)

F1 = array of pre-disturbance values
F2 = array of post-disturbance values

The MMMT_once function performs the MMMT for two data sets, F1 and F2, returning the test statistic and p-value for the MMMT. The notation used in MMMT_once is similar to that used in Mowafy et al. (2023). While here this function is called by the function MMMT, MMMT_once could be employed to estimate the test statistic and p-value for the MMMT for any 2 samples.

4) FDR(MMMT_out)

MMMT_out = matrix of output from MMMT, including the first column of p-values for MMMT test at each percentile.

The FDR function performs test of field/global significance using results of the MMMT test across all percentiles using a multiple False Discovery Rate (FDR) Test (see Mowafy et al., 2023). The result is a flag which equals 3 if the p-value of the FDR test < 0.01, 2 if 0.01 < p-value of FDR test < 0.05, 1 if 0.05 < p-value of FDR test < 0.10, and 0 if p-value of FDR test > 0.10. Including this result is an option in the output figure.

5) Plot_func(pre, post, MMMT_out, p, label, logy, ymin, ymax, ex, FDR, FDR_flag, pdfname)

pre = output matrix of Percentiles from running prepro on pre-disturbance record.
post = output matrix of Percentiles from running prepro on post-disturbance record.
MMMT_out = output matrix from running MMMT using pre and post.
p = array of cumulative distribution function values at which percentiles are calculated using the Harrell and Davis (1982) Distribution-Free Quantile Estimator for data in each year between the Start_Month and End_Month.
label = array of 5 characters which contains:
	label[1] = Title of plot (placed on top of figure)
	label[2] = y-axis label
	label[3] = x-axis label
	label[4] = Name/description of pre-disturbance data set (for legend)
	label[5] = Name/description of post-disturbance data set (for legend)
	Note that the user may have to adjust label lengths to properly fit on figure.
ex = logical character that is either TRUE or FALSE. When TRUE, the exceedance probabilities for each percentile will be plotted on the x-axis; when FALSE, the cdfs for each percentile (1 – the exceedance probabilities) will be plotted on the x-axis. 
logy = logical character that is either TRUE or FALSE. When TRUE, the y-axis is plot on a log base-10 scale; when FALSE, the y-axis is plot on a linear scale.
ymin = minimum value on y-axis
ymax = maximum value on y-axis
FDR = logical character that is either TRUE or FALSE. When TRUE, the result of the FDR test are included on the output figure; when FALSE, the result of FDR test is not included.
FDR_out = output from the FDR function
pdfname = character string with name of the pdf file to be created by Plot_func.

The function Plot_func makes a pdf of a plot of the output from the MMMT function. The pre and post matrices are used to plot the percentiles for each year (light lines on the plot) and calculate the medians at each percentile (darker symbols), while the p-values from the MMMT are used to color code the medians for each percentile (red = p-value < 0.01; orange = 0.01 < p-value < 0.05; yellow = 0.05 < p-value < 0.10, and blue = p-value > 0.10). This same color coding is used to represent the result of the FDR test using the color of the text “FDR p-value” when FDR is TRUE.

References:

Fligner, M.A., Rust, S.W. 1982. A modification of Mood’s median test for the generalized Behrens-Fisher problem. Biometrika, 69, 221–226, doi.org/10.1093/biomet/69.1.221.

Harrell, F.E., and Davis, C.E. 1982. A new distribution-free quantile estimator. Biometrika, 69, 635-640, doi.org/10.1093/biomet/69.3.635.

Mowafy, M.H.M.A., Kroll, C.N., and Vogel, R.M. 2023. A graphical non-parametric hydrologic alteration test using flow duration curves. Hydrological Science Journal, doi.org/10.1080/02626667.2023.2248109.
