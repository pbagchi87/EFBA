
Empirical Frequency Band Analysis for Functional Time Series based on Bagchi and Bruce (2020).


## Description:

Instructions for implementing empirical band analysis for nonstationary functional time series

## Dependencies:

Code was developed using R version 4.0.2 ("Taking Off Again"), so code may not function properly on older versions of R. The packages listed in the demo file must also be installed prior to use: fda, fields, viridis, fEBAcpp.

Code has been tested and developed for data of the following dimensions: <= 50000 time points (T), <= 49 tapers (K), <= 1000 observations per block (N). Larger datasets may require more RAM for processing.

## Quick start guide:

Follow the steps below to simulate data, run the search algorithm, and create summary tables and plots of the results. Simulated data is according to the three settings in the paper: white noise (1 band), linear time-varying dynamics (3 bands), and sinusoidal time-varying dynamics (3 bands).

Download the two R files (fEBA_demo.R and fEBA_Rfns.R) and one tar.gz file (fEBAcpp_1.0.tar.gz) into a folder of your choosing. In what follows, we will assume the folder is saved with the following path 'C:\EBAdemo'.

Open R (version 4.0.2 or newer is recommended) and open the file 'C:\EBAdemo\fEBA_demo.R'. This script contains all necessary code to simulate data pertaining to the processes described in the paper and run the search algorithm.

Follow the instructions in the comments of the demo file to simulate data and apply the search algorithm.

## Using fEBA procedure on other data:

The following inputs are necessary to run the fEBA:

**'X'** is a matrix (T x R) containing a realization of the functional time series process you wish to analyze. T is the length of the time series and R is the number of points in the functional domain where the function is observed. Each row corresponds to a realization of the function for a specific point in time.

**'Rsel'** is the number of points in the functional domain to use in computing test statistics.

**'K'** is a number representing how many tapers to use in estimating the approximately stationary local power spectrum using multitaper spectral estimation. Note here that K must be significantly smaller than N in order to achieve reasonable frequency resolution. For example, if you wish to distinguish behavior for frequencies separated by 0.01 or larger (this is the bandwidth) and you had N=1000 observations in each approximately stationary block, then you could have at most K=9 tapers. More generally, using the sine tapers, bw=(K+1)/(N+1).

**'N'** is a number representing how many observations should be contained in each approximately stationary block for the estimation procedure of the local approximately stationary power spectrum. For example, N=1000 means that the time series is broken up into approximately stationary segments each containing 1000 observations. Note that N must be significantly smaller than the total length of the time series.

**'ndraw'** is the number of random draws from the Gaussian process needed to approximate p-values for the test statistic.

**'alpha'** is the significance level for testing each frequency partition. For example, alpha=0.05 corresponds to the 5% significance level or 95% confidence level.

**'std'** is a binary indicator to determine if the variance in each stationary block should be standardized to unit variance (std=TRUE) or not (std=FALSE).

**'blockdiag'** is a binary indicator to determine if the covariance matrix for the Gaussian process should be approximated with a block diagonal structure (blockdiag=TRUE) or not (blockdiag=FALSE).

**'dcap'** is the number of frequencies to test in a single pass.

Once you have created these inputs, you can pass them into the 'fEBA.wrapper' function just as in the demo file for estimation.
