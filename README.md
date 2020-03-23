# SmoothingSplines
Smoothing Splines and lagrangian parameter optimization
This software generates smoothing splines using different algorithms :
\n -> By giving an arbitrary smoothing paramter p
\n -> By caping the global error S tolerated which permits to deduce the optimal smoothing paramter p
\n -> By giving the spline's degree of freedom which permits to deduce the optimal smoothing paramter p

These differents methods are perfectly explained in the following paper (written in french) : 
https://pro.univ-lille.fr/fileadmin/user_upload/pages_pros/francois_boulier/GIS4-CNUM/support.pdf

The entire software is structurized with a Python skeleton code calling Fortran subroutines using f2py3.

1. Compilation with "make"
================ Code editing =====================
2. Set variables pays with any country you want.
   /!\ Data associated to the chosen country must be into the data file fish.csv /!\ 
3. Set any degree of freedom you want, or "max" if the smoothing spline must go through each data point.
=============== Libraries =========================
4. Check that numpy and matplotlib are well installed in your machine.
=============== Launching =========================
5. $ python3 spline.py

Fortran Subroutines : 
-> cholesky.f : Cholesky Decomposition Algorithm, generates a Lower Triangular Matrix
-> df.f : Degrees of freedom computing for all p in [a,b]
-> init_HG.f : H and G matrices initialization
-> init_Q.f : Q matrix initialization
-> init_T.f : T matrix initialization
-> newton_f.f : Computing p with a global error S chosen using Newton-Raphson Algorithm
-> outils.f : Various tools used for a,b,c, and calculation or an algorithm inversing a SPD matrix
-> splines.f : Splines generation 

Python Code : 
-> spline.py : Software's skeleton calling fortran subroutines and ploting generated splines

Data : 
-> example*.txt : Toy examples
-> fish.csv : Fisheries production per country between 1960 and Now

Other : 
-> Images Folder : Contains various images produced by the software on different data sets
-> Makefile : Compiles Fortran subroutines referenced at line 10
-> Rapport.pdf : The report explaining the different algorithms used in the sofwatre (written in french)
