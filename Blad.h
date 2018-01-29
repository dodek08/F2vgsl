// 7.09.2015
// author: Dominik Kasperski
// kasperski.dominik@gmail.com

#ifndef BLAD
#define BLAD

#include <iostream>

using namespace std;

//klasa wyjatku
class Blad {
public:
const char * zaw;
double a1;
double a2;
double a3;
double a4;
Blad(const char*); //konstruktor
Blad(const char* ,double,double,double,double); // konstruktor
void what(); //wypisuje string
};

#endif // BLAD
