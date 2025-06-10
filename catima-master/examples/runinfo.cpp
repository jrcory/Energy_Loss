
#include "runinfo.h" 
using namespace info;
/*
class runinfo{
    public:
    int runnb,A,Z,gasel;
    int gasA[4]={0};
    int gasZ[4]={0};
    double Ein,density;
    void setinfo(int run);
}; 
*/
void runinfo::setinfo(int nb){
	if (nb>2446 && nb < 2452){
	    A=20;
	    Z=8;
	    gasel=2;
	    gasA[0]=12;
	    gasA[1]=1;
	    gasZ[0]=6; gasZ[1]=1;
	    gasS[0]=1; gasS[1]=4;
	    Ein=54;
	    density = 8.7759e-5;
	    lab = "20O";
	  
	}
	if (nb ==2184 || nb == 2185){
	    A=18;
	    Z=8;
	    gasel=2;
	    gasA[0]=12;
	    gasA[1]=1;
	    gasZ[0]=6; gasZ[1]=1;
	    gasS[0]=1; gasS[1]=4;
   	    Ein=60;
	    density = 1.3e-4;
	    lab = "18O"; 
	}
	if (nb==-1){
	    A=30;
	    Z=16;
	    gasel=1;
	    gasA[0]=4;
	    gasZ[0]=2;
	    gasS[0]=1;
	    Ein=60;
	    density = 1.3e-4;
	    lab = "30S";	
	}
	runnb=nb; 
	return;
}
 
