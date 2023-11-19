# Zeta Potential
The file zeta_v17.m takes the time series streaming current data as input and extracts the zeta potential signal from it.
By Siddharth S. Sahu

Detailed steps:

1) Import the data
2) Data selection: possibility to select either the entire data, or only chunks of it, manually excluding noisy portions of the data. The code then further excludes the outliers in the data, and selects the rest for further processing.
3) As the streaming current data are obtained at two different alternating pressure steps (a high and a low pressure), it is necessary to segregate the streaming current data corresponding to these two pressure values. This is done in the code by manually setting a threshold.
4) The set of data points corresponding to each pressure pulse is then averaged.
5) This is followd by a pieace-wise non-linear fit and extrapolation to fill up the missing low and high values for the alternative pressure pulses. This is done so that the the current difference can be calculated for the low and high values of pressure even though it is not possible to measure the currents corresponding to the low and high pressure simultaneously. It was instead done alternatingly with a time difference of 15 seconds.
6) This was followed by the calculation of the zeta potential using the relevant expression known from theory.
7) The code has the option of further filtering the data using Savitzky-Golay filter.
