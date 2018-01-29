//
// Created by Joanna on 02.09.2015.
//
#include "Blad.h"
#include "Data.h"

Data::Data(string name) {
    file.open(name, ios::in);
}
void Data::read_File() {

    if(file.good() == true){

        cout << "Uzyskano dostep do pliku \n";
        double tmp, tmp_x, tmp_q2, tmp_y, tmp_cs, tmp_error_cs, tmp_f2;
        while(!file.eof()) {
            file >> tmp;
            file >> tmp_q2;
            file >> tmp_x;
            file >> tmp_y;
            file >> tmp_cs;
            file >> tmp_f2;
            file >> tmp;
            file >> tmp;
            file >> tmp;
            file >> tmp_error_cs;
            for(int i=0; i<114; ++i) file >> tmp;
            q2.push_back(tmp_q2);
            x.push_back(tmp_x);
            y.push_back(tmp_y);
            cs.push_back(tmp_cs);
            error_cs.push_back(tmp_error_cs);
            f2.push_back(tmp_f2);
        }
        file.close();

    } else cout << "Dostep do pliku zostal zabroniony! \n";
}
void Data::write_File(string name, double tmp1, double tmp2, double tmp3, double tmp4){

    file.open(name, ios::out | fstream::app);
    if(file.good() == true){

        file << tmp1;
        file << " ";
        file << tmp2;
        file << " ";
        file << tmp3;
        file << " ";
        file << tmp4;
        file << "\n";

        file.close();

    } else cout << "Dostep do pliku zostal zabroniony! \n";
}

vector<double> Data::getX() {
    return x;
}
vector<double> Data::getQ2(){
    return q2;
}
vector<double> Data::getY(){
    return y;
}
vector<double > Data::getCS(){
    return cs;
}
vector<double> Data::getErrorCS() {
    return error_cs;
}
vector<double> Data::getF2(){
    return f2;
}
