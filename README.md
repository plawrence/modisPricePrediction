modisPricePrediction
====================

Note: This analysis was performed as a class project.  The code is rough and commented sparsely, but it does the trick.

Data from NRCS Soil Climate Analysis Network (SCAN) sites, MODIS satellite Enhanced Vegetation Index (EVI) pixels, and the price of corn were all tested to determine Granger Causality (GC) of wheat prices, and were simultaneously assessed for forecasting ability. All years from 2007 to 2012 were independently evaluated, and all variables were first-differenced to achieve stationarity. The following variables were found to GC wheat prices, with the associated predictive lags: in 2007, precipitation and MODIS EVI data [lag = 75 days]; in 2008, the mean temperature, the maximum temperature, and potential evapotranspiration [lag = 60 days]; in 2010, the price of corn [lag = 75 days]; in 2011, the mean temperature and maximum temperature [lag = 75 days]; in 2012, soil moisture at 8 inches [lag = 60 days]. 
