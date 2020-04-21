# covid-logistic


FORTRAN Code to fit a logistic function to a pandemic curve and compute the time lag between infection and death

Programs used to produce results in the paper ‘COVID-19 infection-to-death lag and lockdown-effect lag in Italy and Spain’ by Jorge Buescu, Henrique M. Oliveira and Carlos A. Pires. 

Fortran Program logistic1.f 
Computes the fitting logistic function to a time series correspondent to a pandemic curve, total infected, actively infected or deceased. The parameters are read from the file logistic1.p and their meaning is explained in logistic1.f as comment. The program includes a n iterative subroutine based on the multivariate Newton Raphson algorithm. The ‘linrg’ subroutine from IMSL package is used. If inexistent in the running platform, it must be replaced by another subroutine, be compiled and linked with the appropriate library. 
The computation is performed by minimizing a sum of squared errors (cost function) aling a given period starting at day t1 and ending a day t2. The program makes the fitting for a sequence of t2 values from a smallest to a largest value. 



Fortran Program logistic2.f 
Computes the distance between the curve of daily new actives and the lagged curve of daily new deceased. Values are written in a file and the lag for which the distance is minimal is considered as a measure of the mean lag between infection and death. The parameters are read form file logistic2.p.

The files ‘Italia.acum-actives’ and ‘Spain.acum-actives’ contain data extracted from https://www.worldometers.info/coronavirus/#countries and include:

Number of the day (0=15 February 2020), total number of infected, total number of actively infected , total number of deceased.

