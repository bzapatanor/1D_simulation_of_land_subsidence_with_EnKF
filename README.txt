This folder contain:

1) SGS: R Script to generate an ensemble Y(lnK) realizations. Due to the RandomFields package
   being obsolete, provide an alternative method that uses the gstat package.

2) 1D algorithm: the 1D nonlinear subsidence algorithm developed by Neuman et al. (1982),
   modified by Rudolph and Frind (1991), and implemented by Zapata-Norberto, 2019 
   (10.5281/zenodo.13864158).

3) EnKF: The coupled 1D nonlinear subsidence algorithm and ensemble Kalman filter for   
   heterogeneous and highly compressible aquitards, initially developed by Zapata-Norberto, 2019.


Details of workflow can be consulted in Zapata-Norberto, 2019.

Workflow:

1) Generate an ensamble of Y(lnK) realizations with "SGS".

2) Generate one reference realization of Y with "SGS".

3) Generate the error of reference realization (of h and Y) with "SGS".

4) Replace the file "K.txt" by reference realization, generate in (2), and run the "1D algorithm".

5) On desired positions and times, extract the h and/or K observations, joined with their error 
   generated in (2)

6) To simulate the coupled 1D nonlinear subsidence algorithm with ensemble Kalman filter, in EnKF:
  
	a) Replace the file "realizaciones_K.txt" with the ensemble of Y realizations. 

	b) Replace the file "mesures.txt" with the h and/or Y measurement set with their error generated in (5).
