//
// Created by Joanna on 02.09.2015.
//
//KLASA OPAKOWUJACA DANE Z PLIKU

#ifndef INTERPOLACJA_FILE_H
#define INTERPOLACJA_FILE_H
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class Data {
private:
    fstream file;
    vector<double> x;                                                       //3 WEKTORY PRZECHOWUJACE 3 KOLEJNE KOLUMNY DANYCH Z PLIKU
    vector<double> q2;
    vector<double> y;
    vector<double> cs;
    vector<double> error_cs;
    vector<double> f2;

public:
    Data(string name);
    void read_File();                                                       //METODA SCZYTUJACA DANE Z PLIKU DO WEKTOROW
    void write_File(string name, double tmp1, double tmp2, double tmp3, double tmp4);
    vector<double> getX();                                                  //METODY ZWRACAJACE POSZCZEGOLNE WEKTORY
    vector<double> getQ2();
    vector<double> getY();
    vector<double> getCS();
    vector<double> getErrorCS();
    vector<double> getF2();
};

#endif //INTERPOLACJA_FILE_H
