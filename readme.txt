Reference:

Saha, S., & Bradley, J. R. (2024). Incorporating Subsampling into Bayesian Models for High-Dimensional Spatial Data. Bayesian Analysis, 1(1), 1-40.

----------------------------------------------------------

Steps to run the spatial data subset model (SDSM).

Step 1: Run the "DataSimulator.R" code to generate 20M data. The data is saved as "SimulatedData.RData".

Step 2: Run the "SDSM.R" code to fit the SDSM on the simulated dataset "SimulatedData.RData" and to make predictions at the missing locations.

----------------------------------------------------------

Steps to generate dataset with different data size.

Step 1: Go to the "DataSimulator.R" code.

Step 2: Change the value of Nn. Note that, this is the total number of observations (observed + missing).

Step 3: Change the value of m (number of missing locations).

Step 4: Change xnum and ynum depending on the value of Nn for number of grid points in the direction of x-axis and y-axis respectively.

Step 5: Run the script.

----------------------------------------------------------
