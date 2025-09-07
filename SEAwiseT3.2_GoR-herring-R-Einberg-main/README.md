SEAwiseT3.2_GoR-herring-R
SEAwise Analysis and data of Gulf on Riga herring.  Individual growth of the Gulf of Riga (Baltic Sea) spring spawning herring

data.R:

Data sources. Environmental data, i.e. air temperature is obtained from Kihnu weather station (58.09861°N, 23.97028°E), provided by the Estonian Weather Service. 
Zooplankton data from BIOR. Gulf of Riga herring weight and SSB values are obtained from ICES 2023 report: https://doi.org/10.17895/ices.pub.23123768.v1 and Ojaveer et al. 2004.
Time-period: 1961-2020.
Data files (does not contain raw data, except SSB, age-1 abundance and 1-year old fish WAA): updated data ifle GoR_WAA1_2025_16_07

model.R:

Model_01_WAA1_GAM. GAM models of WAA 1 (1961-2020); and GAM models of WAA 1 of separate time-series: TS1 (1961-1988) and TS2 (1989-2020)

output.R:

Output 1. Time-series figures (1961-2020)
Output 2. Model output figures. GAM figures for GoR herring WAA 1 (1961-2020); GAM figures for GoR herring WAA 1 TS1 (1961-1988) and TS2 (1989-2020).
Output 3. Sliding window analysis (1961-2020) figures
