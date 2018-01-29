/// gsl!
/// elementy macierzowe z Kimbera
/// Dominik Kasperski
/// kasperski.dominik@gmail.com
/// 25.08.2015
///
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <tuple>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <random>
#include <chrono>
#include <gsl/gsl_rng.h>
#include <utility>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>

using namespace std;

//stale uzywane w obliczeniach
const double als= 0.2;
constexpr double sig0=(2912./100.)/(389./1000.);
constexpr double lambda= 277./1000.;
constexpr double x0=0.41*0.0001;
const double Q02=1.;
const double charge2[]= {4./9., 1./9., 4./9., 1./9.}; //tablica z ladunkami kwarkow, ale niepotrzebna
constexpr double chargesq=10./9.; //uzyte w integrandFL/FT
constexpr double coeff=als/(4*M_PI); //stare (3.*sig0)/(4.*M_PI*M_PI*als);
constexpr double coeff2= 2.0/(M_PI*M_PI)*als*2.*M_PI/2.;
constexpr double coeff3= 1./(4.0*M_PI*M_PI)*als*2.*(M_PI/2.0);
constexpr double phiav=(M_PI/4.);
const double mq2[]={0.,0.,0.,9./4.};

const char name[] = "unintegrated_gluons_proton_BK_lin.dat"; //nazwa pliku z danymi


//obsluga siatki
void wczytaj(); //wczytuje siatke z pliku
double pff(const double&,const double&);//points from file, returns points from file
int find_index_gt(const vector<double>&,const double& val);//finds index of a nelement in a vector which is greater than argument
double fg(const double&,const double&); //interpolated lattice

//funcions eliminating compilation errors due to the some conflict of libraries
vector<double>::iterator max_eleme ( vector<double>::iterator first, vector<double>::iterator last );
vector<double>::iterator min_eleme ( vector<double>::iterator first, vector<double>::iterator last );
bool my_isnan(double f);
void setup_gsl_rng(); //nic chyba nie robi

//F2
//functions from mathematica notebook
double Qs2(const double & x);//not used since we use lattice
double fgg(const double& x, const double& kt); //gluon density only for test prupouse
double D1(const int & i, const double & K, const double & B, const double & Q2);
double J1(const int & i, const double & K, const double & B, const double & Q2);
double J2(const int & i, const double & K, const double & B, const double & Q2, const double & kt);
double J3(const int & i, const double & K, const double & B, const double & Q2, const double & kt);
double J4(const int & i, const double & K, const double & B, const double & Q2);
double J5(const int & i, const double & K, const double & B, const double & Q2, const double & kt);
double J6(const int & i, const double & K, const double & B, const double & Q2, const double & kt);
double invzav(const int & i, const double & K, const double & B, const double & Q2, const double & kt);
double impactFT(const double& x, const double& Q2, const double& kt, const double& K, const double& B, const int& i);
double integrandFT(const double& x, const double& Q2, const double& kt, const double& K, const double& B);
pair<double,double> FL(const double& x, const double& Q2); // (result,error) blad wyznaczenia to srednia wazona bledow
double impactFL(const double& x, const double& Q2, const double& kt, const double& K, const double& B, const int& i);
double integrandFL(const double& x, const double& Q2, const double& kt, const double& K, const double& B);
pair<double,double> FT(const double& x, const double& Q2); // (result,error) blad wyznaczenia to srednia wazona bledow

double f2(double *args, size_t dim, void *params); //F2
double ft(double *args, size_t dim, void *params);
double fl(double *args, size_t dim, void *params);
void set_iterations(double);
void warm_up(); //ustawia generator i strukture do calkowania z gsl
void cool_down(); //uwalnia pamiec
pair<double,double> F2(const double& x, const double& Q2);// (result,error) blad wyznaczenia to srednia wazona bledow
void set_iterations(int); //ustawia liczbe iteracji
size_t liczba_iteracji(void); //informative function about number of iterations
struct pars {double x; double Q2; }; //struktura przechowuje parametry
