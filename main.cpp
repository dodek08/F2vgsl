///gsl!
#include "F2.h"
#include "Data.h"
#include "Steer.h"
#include "Blad.h"


int main()
{

try{


Steer s("steering_file");
s.printS();
set_iterations(s.iterations);
ofstream save,report;
save.open("F2_values_M.dat", ios::out |ios::app );
report.open("reportM.txt", ios::out | ios::app);
const char b1[]={"blad otwarcia pliku zapisu wartosci"};
const char b2[]={"blad otwarcia pliku raportu"};
if(not save.is_open())
    throw Blad(b1);
if(not report.is_open())
    throw Blad(b2);
wczytaj();

cout<<setprecision(6)<<endl;
save<<setprecision(16)<<endl;

vector<double> xx;
vector<double> Q22;
pair<double,double> F2val;

if(s.FILE)
{
    clock_t start = clock();
    Data dane("nce+p.txt");
    dane.read_File();

    vector<double> x = dane.getX();
    vector<double> q2 = dane.getQ2();
    vector<double> y = dane.getY();
    vector<double> cs = dane.getCS();
    vector<double> error_cs = dane.getErrorCS();
    vector<double> f2 = dane.getF2();
    double cs_calc, chi_kw=0, real_error_cs;
    warm_up();
    cout<<"x \t Q2 \t F2 \t F2_error \t cs_calc \t real_error_cs \t chisq_i"<<endl;
    for(unsigned int i=0; i<x.size() ;++i){
        if(f2[i]==9.990) continue;
        else {
            F2val=F2(x[i],q2[i]);
            cs_calc = ( F2val.first - (y[i] * y[i] / (1 + (1 - y[i]) * (1 - y[i]))) * FL(x[i], q2[i]).first);
            real_error_cs = error_cs[i] * cs[i] / 100;
            chi_kw += (cs[i] - cs_calc) * (cs[i] - cs_calc) / (real_error_cs);
            dane.write_File("cs_liczone_i_blad_chikw2.dat", x[i], q2[i], cs_calc, real_error_cs);
            cout<<x[i]<<"\t"<<q2[i]<<"\t"<<F2val.first<<"\t"<<F2val.second<<"\t"<<cs_calc<<"\t"<<real_error_cs<<"\t"<<chi_kw<<endl;
            save<<x[i]<<"\t"<<q2[i]<<"\t"<<F2val.first<<"\t"<<F2val.second<<endl;
        }
    }
    cool_down();

    dane.write_File("cs_liczone_i_blad_chikw.dat",0,0,0,chi_kw);
    cout << "Chi kwadrat = " << chi_kw << endl;
    cout << "Czas wykonywania: " << (double)(clock()-start)/(CLOCKS_PER_SEC) << "sek" << endl;
    report<<"czas realizacji (dane z pliku): "<<(double)(clock()-start)/(CLOCKS_PER_SEC)<<", liczba iteracji: "<<liczba_iteracji()<<endl;
}
else
{
if(s.limits)
{
double xv, xi=s.xi, xf=s.xf;  //granice x
double Q2v, Q2i=s.Q2i, Q2f=s.Q2f; //granice Q2


int imax;
if(s.logx)
    imax=log10(xf)-log10(xi);
else
{
    imax=(xf-xi)/s.nx;
}

int jmax;
if(s.logQ2)
    jmax=log10(Q2f)-log10(Q2i);
else
{
    jmax=(Q2f-Q2i)/s.nQ2;
}

for(int i=0;i<=imax;i++)
    {
    xv=xi*pow(10,i);
    cout<<xv<<endl;
    xx.push_back(xv);
	for (int j=0; j<=jmax;j++)
	{
        Q2v=Q2i*pow(10,j);
        cout<<Q2v<<endl;
        Q22.push_back(Q2v);
        }
        }
}
if(s.vectors)
{
    xx=s.xvals;
    Q22=s.Q2vals;
}


cout<<"(x,Q2)"<<endl;
for(double x:xx )
    {
	for (double Q2 : Q22)
        {
            cout<<"("<<x<<","<<Q2<<")"<<endl;
        }
	}

cout<<"x \t"<<"Q2 \t"<<"F2(x,Q2)"<<endl;
clock_t t;
t=clock();
warm_up();
for(double x:xx )
    {
	for (double Q2 : Q22)
	{
        F2val=F2(x,Q2);
		save<<x<<"\t"<<Q2<<"\t"<<F2val.first<<"\t"<<F2val.second<<"\n";
		cout<<x<<"\t"<<Q2<<F2val.first<<"\t"<<F2val.second<<endl;
    }
    }
cool_down();
t=clock()-t;
report<<"czas realizacji: "<<((double)t)/CLOCKS_PER_SEC<<", liczba iteracji: "<<liczba_iteracji()<<endl;
cout<<"czas realizacji: "<<((double)t)/CLOCKS_PER_SEC<<", liczba iteracji: "<<liczba_iteracji()<<endl;

}//if not file...
save.close();
report.close();
}

catch(Blad B)
{
	B.what();
	cout<<B.a1<<" "<<B.a2<<" "<<B.a3<<" "<<B.a4<<" "<<endl;
}
return 0;
}
