// This is just a really hacky way to load in all the
// classes to demonstrate how to use the energy loss model.
// Really these should be compiled into a shared library
#include "AtELossTable.h"
#include "AtELossTable.cxx"
#include "AtELossModel.h"
#include "AtELossModel.cxx"
#include "AtSpline.h"
#include "AtSpline.cxx"
#include "AtStringManip.h"
#include "AtStringManip.cxx"

using namespace AtTools;  // All these classes live in the AtTools namespace

void example_music()
{
    //Load in energy loss model for 30S in He gass.
    AtELossTable eModelSrim;
    eModelSrim.LoadSrimTable("./30S_4He_SRIM.txt");

    // Load in energy loss model for 30S in He gas from LISE++
    AtELossTable eModelLise;
    eModelLise.LoadLiseTable("./30S_4He_Lise_um.txt", 30, -1.3e-4, 2); 
    //62 is mass, next is density (only matters if you are chaning 
    // the density of the material while making calculations. Units are g/cm^3.
    // Sign of density sets the units used when exporting the table in LISE.
    // 2 is which column in the table to use (each is a different model). ATIMA-LS is 2.

    double inEn = 74; //Initial energy in MeV
    double distance = 12.5; //Distance in mm

    double E_out_S[21]={0};
    double Eloss_S[20]={0};
    double E_out_L[21]={0};
    double Eloss_L[20]={0};

    E_out_S[0]=E_out_L[0]=inEn;
    for (int i=0; i<20; i++){
	Eloss_S[i]=eModelSrim.GetEnergyLoss(E_out_S[i],distance);
	Eloss_L[i]=eModelLise.GetEnergyLoss(E_out_L[i],distance);
	E_out_S[i+1]=eModelSrim.GetEnergy(E_out_S[i], distance);
	E_out_L[i+1]=eModelLise.GetEnergy(E_out_L[i], distance);

    }

//add in root functionality to plot the results: 
   double x[20]={};
   for (int i=0; i<20; i++){
//	std::cout << "lise: " << Eloss_L[i] << std::endl;
//	std::cout << "srim: " << Eloss_S[i] << std::endl;
        x[i]=i;
   }
   
   TGraph *lise = new TGraph(20,x,Eloss_L);
   TGraph *srim = new TGraph(20,x,Eloss_S);
   TCanvas *c1 = new TCanvas();
   lise->GetHistogram()->SetMaximum(6);
   lise->Draw();
   srim->SetLineColor(2);
   srim->Draw("same");
   auto legend = new TLegend(0.7,0.8,0.9,0.9);
   legend->AddEntry(lise,"Lise");
   legend->AddEntry(srim,"SRIM");
   legend->Draw();

}
