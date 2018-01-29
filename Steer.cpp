// 7.09.2015
// author: Dominik Kasperski
// kasperski.dominik@gmail.com

#include "Blad.h"
#include "Steer.h"

Steer::Steer(string s)
{
    fstream file(s);
    char c;
    double tmp;
    string vs;

    while(file.get(c))
    {
        if(c==':')
        {
            file.get(c);
            file.get(c);
            if(c=='y')
            {
                FILE=true;
                vectors=false;
                limits=false;
            }
            else
                FILE=false;
            break;
            }
    }

    while(file.get(c))
    {
        if(c==':')
        {
            file>>iterations;
            break;
        }
    }
if(not FILE)
{
    while(file.get(c))
    {
        if(c==')')
        {
            file.get(c);
            file.get(c);
            if(c=='y')
                vectors=true;
            else
                vectors=false;
            break;
            }
    }

    if(vectors)
    {
        limits=false;
        while(vs!="form:")
        {
            file>>vs;
        }
        while(!file.eof())
        {
            file>>tmp;
            cout<<tmp<<endl;
            xvals.push_back(tmp);
            file>>tmp;
            cout<<tmp<<endl;
            Q2vals.push_back(tmp);
        }
    }
    else
    {
        limits=true;
        while(file.get(c))
        {
        if(c==':')
        {
            file>>xi;
            file>>xf;
            break;
        }
        }

        while(file.get(c))
        {
        if(c==':')
        {
            file>>Q2i;
            file>>Q2f;
            break;
        }
        }

        while(file.get(c))
        {
        if(c==')')
        {
            file.get(c);
            file.get(c);
            if(c=='y')
                logx=true;
            else
                logx=false;
            break;
            }
        }

        while(file.get(c))
        {
        if(c==')')
        {
            file.get(c);
            file.get(c);
            if(c=='y')
                logQ2=true;
            else
                logQ2=false;
            break;
            }
        }

    while(file.get(c))
        {
        if(c==':')
        {
            file>>nx;
            break;
            }
        }
    while(file.get(c))
        {
        if(c==':')
        {
            file>>nQ2;
            break;
            }
        }

    }//else
} //if FILE
    file.close();
    cout<<"steering file red "<<s<<endl;
}

void Steer::printS()
{
    if(FILE)
        cout<<"data from file mode"<<endl;
    else
    {
    if(limits)
    {
    cout<<"limits"<<endl;
    cout<< xi<<" "<< xf<<"\n"<< Q2i<<" "<< Q2f<<"\n";
    if(logx)
        cout<<"x logscale"<<endl;
    if(logQ2)
        cout<<"Q2 logscale"<<endl;
    if(!logx)
    {
        cout<<"number of x points:"<<endl;
        cout<< nx<<endl;
    }
    if(!logQ2)
    {
        cout<<"number of Q2 points:"<<endl;
        cout<< nQ2<<endl;
    }
    }
    else
    {
    cout<<"point values"<<endl;
    vector<double>::iterator it=xvals.begin();
    vector<double>::iterator it2=Q2vals.begin();
    cout<<"x \t Q2"<<endl;
    while(it!=xvals.end())
    {
        cout<<*it++<<"\t"<<*it2++<<endl;
    }
    cout<<"number of iterations: "<<iterations<<endl;
    }
    }
}

