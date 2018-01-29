7.09.2015
author: Dominik Kasperski (matrix elements, intergration, steering file)
e-mail: kasperski.dominik@gmail.com
Author: Joanna Wańczyk (reading data from file, cross section, X^2)
e-mail: joannawanczyk@gmail.com

F2 values are calculated with formulae 3.2N from "Unintegrated Parton Distributions" by M. A. Kimber (2001) 
check if you have following files:
main.cpp, F2.h, F2.cpp, Blad.h, Blad.cpp, Data.h, Data.cpp, Steer.h, Steer.cpp and makefile
you need to have GSL (v. 1.16) library installed (if not you can download it from http://ftp.task.gda.pl/pub/gnu/gsl/gsl-1.16.tar.gz)
installation:
1. unzip downloaded file
2. type 
$./configure
in directory with the unzipped library
3. type
$make
4. type
$sudo make install
if you are in trouble with gsl, please check "README" and "INSTALL" files in gsl directory for help.

after that in directory with a source code of F2 type 
$make
F2.exe is the executable file, so to run it type
$./F2.exe

using "steering_file" you can choose number of iterations and mode: reading data from file (x, y, Q2, F2, cross section, total error) or specifing what points will be calculated.
All F2 values calculated would be in file "F2_values_M.dat", all times of execution in "reportM.txt".
Data about chi^2 etc in "cs_liczone_i_blad_chikw2.dat"

