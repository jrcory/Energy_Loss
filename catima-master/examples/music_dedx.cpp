#include "catima/catima.h"
#include <iostream>

using std::cout;
using std::endl;


int main(){
    catima::Material he({ // he material with one atom
        {4,2,0}, // 1H - two atoms
    });
    he.density(1.3e-4).thickness(1.25);//assuming this is in g/cm^3 and cm, might need to check

    catima::Projectile p(30,16); // define projectile, ie 30S

    cout<<"30S->He\n";
    for(double T=0; T<90;T+=1){
        auto result = catima::calculate(p,he,T/12); 
	    cout<<"T = "<<T<<" MeV, dEdx = "<<result.dEdxi<<" MeV/g/cm2"<<endl;
	}

    return 0;
}
