///gsl!

#include "Blad.h"
#include "F2.h"

size_t calls=10000; //jakies przykladowe, zmieniaj przez plik

const gsl_rng_type * T = gsl_rng_taus113; //generator do wyboru z gsl
gsl_rng * r=gsl_rng_alloc(T);

vector<double> X;
vector<double> Kt2; //tu jest kt^2, ale w rzeczywistosci to jest teraz kt
vector<double> Fg;
multimap<tuple<double,double>,double> FG; //przyspiesza wyszukiwanie wartosci FG z pliku
vector<double> Xred; //Xreduced czyli X z unikalnymi wartosciami X
vector<double> Kt2red; //j.w. ale z Kt2
double Xa=0.0001; // X dla ktorego wyliczamy rozklad kt2
double XMAX, XMIN, KT2MAX, KT2MIN, KTRANGE;
double Xindexgt;

gsl_monte_function F; //F2
gsl_monte_function G; //FT
gsl_monte_function H; //FL
gsl_monte_vegas_state *s;

void set_iterations(int d)
{
    calls=(size_t)d;
}

vector<double>::iterator max_eleme ( vector<double>::iterator first, vector<double>::iterator last )
{
  if (first==last) return last;
  vector<double>::iterator largest = first;

  while (++first!=last)
    if (*largest<*first)    // or: if (comp(*largest,*first)) for version (2)
      largest=first;
  return largest;
}


vector<double>::iterator min_eleme ( vector<double>::iterator first, vector<double>::iterator last )
{
  if (first==last) return last;
  vector<double>::iterator smallest = first;

  while (++first!=last)
    if (*first<*smallest)    // or: if (comp(*first,*smallest)) for version (2)
      smallest=first;
  return smallest;
}

inline bool my_isnan(double f)
{
    return (f!=f);
}

//obsluga siatki

void wczytaj()
{
	double tmp,tmp1,ret;
	fstream f;
	int i=-1;
	f.open(name, ios::in | ios::out);
	if (!f.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
	while(!f.eof())
	{
        i++;
		f>>tmp;
		f>>tmp1;
		f>>ret;

		X.push_back(tmp);
		Kt2.push_back(sqrt(tmp1));
		Fg.push_back(ret);

		FG.insert(pair<tuple<double,double>,double>(make_tuple(X.back(),Kt2.back()),Fg.back()));
		if(i==1)
            Xred.push_back(X.front());
		if(i>1)
		{
            if(X.at(i)!=X.at(i-1))
                Xred.push_back(X.back());
        }
	}
	f.close();

	XMAX=(*max_element(X.begin(),X.end()));
	XMIN=(*min_element(X.begin(),X.end()));
	KT2MAX=(*max_element(Kt2.begin(),Kt2.end()));
	KT2MIN=(*min_element(Kt2.begin(),Kt2.end()));
    sort(Kt2.begin(),Kt2.end());
    vector<double>::iterator it=Kt2.begin();
    Kt2red.push_back(*it++);
    while(it!=Kt2.end())
    {
        if(*it!=*(it-1))
                Kt2red.push_back(*it);
        it++;
    }
    KTRANGE=KT2MAX-KT2MIN;
}

//returning only points in file
double pff(const double & x, const double & kt2)
{
		multimap<tuple<double,double>,double,double>::iterator it=FG.find(make_tuple(x,kt2));

	if(it!=FG.end() or (X.back()==x and Kt2.back()==kt2))
	{
		if (it->second<0)
            return 0.0;
        else
            return it->second;
    }
	else
		throw Blad("point out of the lattice");
}

int find_index_gt(const vector<double>& vec, const double& val)
{
    int size = vec.size();
    int right=size+1;
    int left=0;
    int index=(right+left)/2;
    //cout<<val<<" "<<vec.size()<<endl;
    bool not_match=true;
    while(not_match)
    {
        index=(right+left)/2;
        if (vec[index]>=val and vec[index-1]<val)
            {
            not_match=false;
            break;
            }
        else
            if(val<vec[index])
                right=index;
            else
                left=index;
        if(index>size or index<1)
            throw Blad("nie dziala find_index_gt",(double)index, val, vec[index], vec[index-1]);
    }
    return index;
}

double fg(const double& x, const double& kt2)
{
    if(x<XMIN or x>XMAX)
        throw Blad("out of range w fg",x,kt2,0.,0.);
	int tab[3]={0}; //table of indexes
	//unsigned int i=1;

	multimap<tuple<double,double>,double,double>::iterator it=FG.find(make_tuple(x,kt2));

	if(it!=FG.end() or (X.back()==x and Kt2.back()==kt2))
	{
		if (it->second<0)
            return 0.0;
        else
            return it->second;
    }
	else
	{
		tab[0]=find_index_gt(Xred,x);
		tab[1]=tab[0]-1;
		tab[2]=find_index_gt(Kt2red,kt2);
		tab[3]=tab[2]-1;
		//tables of point around the choosen by the user
		double lf[] = {Xred[tab[1]], Kt2red[tab[3]], FG.find(make_tuple(Xred[tab[1]], Kt2red[tab[3]]))->second};
		double rf[] = {Xred[tab[0]], Kt2red[tab[3]], FG.find(make_tuple(Xred[tab[0]], Kt2red[tab[3]]))->second};
		double lr[] = {Xred[tab[1]], Kt2red[tab[2]], FG.find(make_tuple(Xred[tab[1]], Kt2red[tab[2]]))->second};
		double rr[] = {Xred[tab[0]], Kt2red[tab[2]], FG.find(make_tuple(Xred[tab[0]], Kt2red[tab[2]]))->second};
        //interpolation
		double a,b,c,d,f,g; //y=ax+b y1=cx+d y2=fx+g
        a = (lf[2]-rf[2])/(lf[0]-rf[0]);
        c = (lr[2]-rr[2])/(lr[0]-rr[0]);
        if(a==0 and c==0)
            return 0.0;
        b = lf[2]-a*lf[0];
        d = lr[2]-a*lr[0];
		double y = a*x+b;
		//double y1 = c*x+d;
        f = (y-(c*x+d))/(lf[1]-lr[1]);
        g = y-f*lf[1];
		double ret= f*kt2+g;
		if (ret<0)
            return 0.0;
        else
          //  cout<<difftime(stop,start)<<endl;
        return ret;
        }
}

//funkcje to funkcje z mathamatica notebooka
//nie jest wywolywana skoro operujemy siatka
inline double Qs2(const double & x)
{
	double ret= pow((x0/x),lambda);
	//if(my_isnan(ret))
	//	throw Blad("nan w Qs2",x,x0,lambda,0.);
	return ret;
}

inline double fgg(const double& x, const double& kt)
{
    double ret= 3.*sig0/(4*M_PI*M_PI*als)*kt*kt/Qs2(x)*exp(-kt*kt/Qs2(x));
	if(my_isnan(ret))
		throw Blad("nan w fgg",x,kt,Qs2(x),0.);
	return ret;
}

inline double D1(const int & i, const double & K, const double & B, const double & Q2)
{
    return K*K+B*(1.-B)*Q2+mq2[i];
}

inline double J1(const int & i, const double & K, const double & B, const double & Q2)
{
    return K*K/(D1(i,K,B,Q2)*D1(i,K,B,Q2));
}

double J2(const int & i, const double & K, const double & B, const double & Q2, const double & kt)
{
    return  1./sqrt(((D1(i,K,B,Q2) + kt*kt)*(D1(i,K,B,Q2) + kt*kt)) - 4*K*K*kt*kt) - (D1(i,K,B,Q2) +kt*kt)* ((B* (1. - B)* Q2) +mq2[i])/pow(((D1(i,K,B,Q2) + kt*kt)*(D1(i,K,B,Q2) + kt*kt) - 4*K*K*kt*kt),(3./2.));
}

double J3(const int & i, const double & K, const double & B, const double & Q2, const double & kt)
{
    return 1./D1(i,K,B,Q2) *(1. + (K*K - B* (1. - B)* Q2 - mq2[i] - kt*kt)/sqrt((D1(i,K,B,Q2) + kt*kt)*(D1(i,K,B,Q2) + kt*kt) - 4*K*K*kt*kt));
}

inline double J4(const int & i, const double & K, const double & B, const double & Q2)
{
    return 1./(D1(i,K,B,Q2)*D1(i,K,B,Q2));
}

double J5(const int & i, const double & K, const double & B, const double & Q2, const double & kt)
{
    return (D1(i,K,B,Q2) +kt*kt)/pow(((D1(i,K,B,Q2) + kt*kt)*(D1(i,K,B,Q2) + kt*kt) - 4*K*K*kt*kt),(3./2.));
}

double J6(const int & i, const double & K, const double & B, const double & Q2, const double & kt)
{
    return 2./(D1(i,K,B,Q2)*sqrt((D1(i,K,B,Q2) + kt*kt)*(D1(i,K,B,Q2) + kt*kt) - 4*K*K*kt*kt));
}

double invzav(const int & i, const double & K, const double & B, const double & Q2, const double & kt)///UWAGA
{
    return 1./(1. + (K*K *mq2[i])/((1. - B)* Q2) + (K*K + kt*kt -2.*K*kt*cos(phiav) +mq2[i])/(B*Q2));
}

double impactFT(const double & x, const double& Q2, const double& kt, const double& K, const double& B, const int& i)
{
    double xnorm = x/invzav(i,K,B,Q2,kt);
    if(my_isnan(xnorm))
		throw Blad("NAN w impactFT");
    if(((1. - xnorm)<0.) or (xnorm<XMIN) or (xnorm>XMAX) )
        return 0.0;
    else
        return als*Q2/(4.*M_PI)*charge2[i]*4.*kt*K*((B*B + (1. - B)*(1. - B))*(J1(i,K,B,Q2) + J2(i,K,B,Q2,kt) -J3(i,K,B,Q2,kt)) + mq2[i]*(J4(i,K,B,Q2) + J5(i,K,B,Q2,kt) - J6(i,K,B,Q2,kt)))* fg(xnorm, kt)/(kt*kt);
}

double integrandFT(const double& x, const double& Q2, const double& kt, const double& K, const double& B)
{
	double tmp2=0.0;
	for(int i=0;i<4;i++)
        {
            tmp2+=impactFT(x,Q2,kt,K,B,i);
        }
	if(my_isnan(tmp2))
        cout<<"NAN"<<endl;
	//	throw Blad("NAN w integrandFT",kt,K,B,0.);
	return tmp2;
}

double ft(double *args, size_t dim, void *params)
{
    struct pars * fp = (struct pars *)params;
    //prams := {x,Q2}
    //args := {kt,K,B}
    return integrandFT(fp->x,fp->Q2,args[0],args[1],args[2]);
}

pair<double,double> FT(const double& x, const double& Q2)
{
	  if(Q2<=0)
    throw Blad("Q2 out of range");
  if(x>=1)
    throw Blad("x out of range");
    gsl_rng_set(r, chrono::system_clock::now().time_since_epoch().count());
  clock_t t;
  t=clock();
  double result=0., error=0.;		// result and error
  double tmp=0,tmp2=0;
  struct pars pms={x,Q2};
  G.params=&pms;
  //Xindexgt=find_index_gt(Xred,x);
  double Kl[]={.01, 1., 5., 10., 50.};
  double Bl[]={ 0., .1, .5, .7,  1.};
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<4; j++)
    {
        double xl[3] = { KT2MIN, Kl[i], Bl[j] }; //kt,K,B
        double xu[3] = { KT2MAX, Kl[i+1], Bl[j+1] };
        gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/10, r, s, &result, &error);
        //s->stage=1;
        do
        {
            result=0.;
            error=0.;
            gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s, &result, &error);
            //cout<<"chisq<< "<<s->chisq<<endl;
            //s->stage=1;
            if(s->chisq==0)
                //gsl_monte_vegas_init(s);
                break;
        }
        while ((fabs (s->chisq - 1.0) > 0.35) ); //w celu uzyskania najwiekszej dokladnosci
        tmp+=result;
        tmp2+=error*result;
        //cout<<"K["<<Kl[i]<<","<<Kl[i+1]<<"], B["<<Bl[j]<<","<<Bl[j+1]<<"], wynik= "<<tmp<<", blad="<<error<<endl;
    }
  }
  t=clock()-t;
  cout<<"time FT: "<<((double)t)/CLOCKS_PER_SEC<<endl;
  return make_pair(tmp,tmp2/tmp);
}

double impactFL(const double& x, const double& Q2, const double& kt, const double& K, const double& B, const int & i)
{
    double xnorm = x/invzav(i,K,B,Q2,kt);
    if(my_isnan(xnorm))
		throw Blad("NAN w impactFL",K,B,Q2,kt);
    if(((1. - xnorm)<0.) or (xnorm<XMIN) or (xnorm>XMAX) )
        return 0.0;
    else
        return als*Q2/(4.*M_PI)*charge2[i]* 4.*kt*K*(4.*Q2*B*B* (1.-B)*(1.-B)* (J4(i,K,B,Q2) + J5(i,K,B,Q2,kt) -J6(i,K,B,Q2,kt))) *fg(xnorm, kt)/(kt*kt);
}

double integrandFL(const double& x, const double& Q2, const double& kt, const double& K, const double& B)
{
	double tmp1=0.0;
	for(int i=0;i<4;i++)
        {
            tmp1+=impactFL(x,Q2,kt,K,B,i);
        }
	if(my_isnan(tmp1))
        cout<<"NAN"<<endl;
	//	throw Blad("NAN w integrandFL",kt,K,B,0.);
	return tmp1;
}

double fl(double *args, size_t dim, void *params)
{
    struct pars * fp = (struct pars *)params;
    //prams := {x,Q2}
    //args := {kt,K,B}
    return integrandFL(fp->x,fp->Q2,args[0],args[1],args[2]);
}

pair<double,double> FL(const double& x, const double& Q2)
{
	  if(Q2<=0)
    throw Blad("Q2 out of range");
  if(x>=1)
    throw Blad("x out of range");
    gsl_rng_set(r, chrono::system_clock::now().time_since_epoch().count());
  clock_t t;
  t=clock();
  double result=0., error=0.;		// result and error
  double tmp=0,tmp2=0;
  struct pars pms={x,Q2};
  H.params=&pms;
  //Xindexgt=find_index_gt(Xred,x);
  double Kl[]={.01, 1., 5., 10., 50.};
  double Bl[]={ 0., .1, .5, .7,  1.};
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<4; j++)
    {
        double xl[3] = { KT2MIN, Kl[i], Bl[j] }; //kt,K,B
        double xu[3] = { KT2MAX, Kl[i+1], Bl[j+1] };
        gsl_monte_vegas_integrate (&H, xl, xu, 3, calls/10, r, s, &result, &error);
        //s->stage=1;
        do
        {
            result=0.;
            error=0.;
            gsl_monte_vegas_integrate (&H, xl, xu, 3, calls/5, r, s, &result, &error);
            //cout<<"chisq<< "<<s->chisq<<endl;
            //s->stage=1;
            if(s->chisq==0)
                //gsl_monte_vegas_init(s);
                break;
        }
        while ((fabs (s->chisq - 1.0) > 0.35) ); //w celu uzyskania najwiekszej dokladnosci
        tmp+=result;
        tmp2+=error*result;
        //cout<<"K["<<Kl[i]<<","<<Kl[i+1]<<"], B["<<Bl[j]<<","<<Bl[j+1]<<"], wynik= "<<tmp<<", blad="<<error<<endl;
    }
  }
  t=clock()-t;
  cout<<"time FL: "<<((double)t)/CLOCKS_PER_SEC<<endl;
  return make_pair(tmp,tmp2/tmp);
}

double f2(double *args, size_t dim, void *params)
{
    struct pars * fp = (struct pars *)params;
    //prams := {x,Q2}
    //args := {kt,K,B}
    return integrandFL(fp->x,fp->Q2,args[0],args[1],args[2])+integrandFT(fp->x,fp->Q2,args[0],args[1],args[2]);
}


void warm_up()
{
  gsl_rng_env_setup ();
  F.f=&f2;
  G.f=&fl;
  H.f=&ft;
  F.dim=3;
  G.dim=3;
  H.dim=3;
  s = gsl_monte_vegas_alloc(3);
}

void cool_down()
{
    gsl_rng_free(r);
    gsl_monte_vegas_free(s);
}


pair<double,double> F2(const double & x, const double & Q2)
{
  if(Q2<=0)
    throw Blad("Q2 out of range");
  if(x>=1)
    throw Blad("x out of range");
    gsl_rng_set(r, chrono::system_clock::now().time_since_epoch().count());
  clock_t t;
  t=clock();
  double result=0., error=0.;		// result and error
  double tmp=0,tmp2=0;
  struct pars pms={x,Q2};
  F.params=&pms;
  //Xindexgt=find_index_gt(Xred,x);
  double Kl[]={.01, 1., 5., 10., 50.};
  double Bl[]={ 0., .1, .5, .7,  1.};
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<4; j++)
    {
        double xl[3] = { KT2MIN, Kl[i], Bl[j] }; //kt,K,B
        double xu[3] = { KT2MAX, Kl[i+1], Bl[j+1] };
        gsl_monte_vegas_integrate (&F, xl, xu, 3, calls/10, r, s, &result, &error);
        //s->stage=1;
        do
        {
            result=0.;
            error=0.;
            gsl_monte_vegas_integrate (&F, xl, xu, 3, calls/5, r, s, &result, &error);
            //cout<<"chisq<< "<<s->chisq<<endl;
            //s->stage=1;
            if(s->chisq==0)
                //gsl_monte_vegas_init(s);
                break;
        }
        while ((fabs (s->chisq - 1.0) > 0.35) ); //w celu uzyskania najwiekszej dokladnosci
        tmp+=result;
        tmp2+=error*result;
        //cout<<"K["<<Kl[i]<<","<<Kl[i+1]<<"], B["<<Bl[j]<<","<<Bl[j+1]<<"], wynik= "<<tmp<<", blad="<<error<<endl;
    }
  }
  t=clock()-t;
  cout<<"time F2: "<<((double)t)/CLOCKS_PER_SEC<<endl;
  return make_pair(tmp,tmp2/tmp);
}

size_t liczba_iteracji(void)
{
    return calls;
}
