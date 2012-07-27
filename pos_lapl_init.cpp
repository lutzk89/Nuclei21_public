#include<iostream>
#include<time.h>
#include<string>
#include<fstream>
#include<cmath>
#include<stdlib.h>
using namespace std;
int main()
{
	srand(time(NULL));
	ifstream conf("Settings.txt");
	ofstream out("POS.txt"),ext("EXT.txt");
	double x=10,y=10;
	int n=0,D=2;
	string word;
	double rho=1;
	while(!conf.eof())
	{
		conf >> word;
		if(word == "rho")
			conf >> rho;
    if(word == "BoxSize")
			conf >> x >> y;
	}
	n=(int)(x*y*rho);
    double X[2*n],a,b;
    for(int kappa=0;kappa<2*n;kappa++)
        X[kappa]=0;
    double minabs=0.5/sqrt(rho);
    bool allow=false,goon=false;
    while (goon==false)
    {
    int set=0,pool=0;
    goon=true;
	while(set<n)
	{
		a=((double)(rand()%10001)/10000)*x;
        b=((double)(rand()%10001)/10000)*y;
        allow=true;
        for(int kappa=0;kappa<set;kappa++)
        {
            if(pow(X[kappa*D]-a,2)+pow(X[kappa*D+1]-b,2) < pow(minabs,2))
                allow=false;
        }
        pool++;
        if(allow)
        {
            out << a << " " << b << endl;
            set++;
            X[set*D]=a;
            X[set*D+1]=b;
		    ext << "0 1\n";
        }
        else
        {
            pool++;
            srand(pool*time(NULL));
        //    goon=false;
        //    break;
        }
	}
            cout << "allowed nr " << set << endl;
    }
}
