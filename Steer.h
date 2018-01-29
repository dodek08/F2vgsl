// 7.09.2015
// author: Dominik Kasperski
// kasperski.dominik@gmail.com

#ifndef STER
#define STER

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//klasa obslugi steering file
class Steer {
public:
double xi;
double xf;
double Q2i;
double Q2f;
double nx;
double nQ2;
bool logx;
bool logQ2;
bool vectors;
bool limits;
bool FILE;
int iterations;
vector<double> xvals;
vector<double> Q2vals;
Steer(string); //konstruktor zczytuje z pliku ktorego nazwa jest przekazana jako srting to co trzeba
void printS(); //wypisuje co zczytalo z pliku

};


#endif // STER

///steering_file text in case the real file is lost:
/*
///answer the questions as the example states, otherwise it might not work properly
///logscale options works only if F2 is not for specified points
do you want to read values of x, Q2, F2, y, s_r etc. from file "nce-p.txt"? then other options would be unavalible (y/n):
y
number of iterations:
10000
Do You want to evaluate F2 function for some special points? (y/n)
y
Otherwise, define the domain
x: 1e-7 1e-2
Q2: 1 100
logharitmic x axis? (y/n)
y
logharitmic Q2 axis? (y/n)
y
if x axis is not logharitmic axis, define number of x values between limits xinitial and xfinal:
10
if Q2 axis is not logharitmic axis, define number of Q2 values between limits Q2initial and Q2final:
10
Write here points (x,Q2) in that form:
1e-7 1
1e-7 10
1e-7 100
*/

