// 7.09.2015
// author: Dominik Kasperski
// kasperski.dominik@gmail.com

#include "Blad.h"

Blad::Blad(const char * init, double q, double w, double e, double r):zaw(init), a1(q), a2(w), a3(e), a4(r) {}

Blad::Blad(const char * init):zaw(init){}

void Blad::what()
{
	cout<<zaw<<endl;
}

