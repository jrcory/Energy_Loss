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

void example_macro()
{
    //Load in energy loss model for Fe26 in He gass.
    AtELossTable eModelSrim;
    eModelSrim.LoadSrimTable("./30S_4He_SRIM.txt");

    // Load in energy loss model for Fe26 in He gas from LISE++
    AtELossTable eModelLise;
    eModelLise.LoadLiseTable("./30S_4He_Lise.txt", 30, 1.3e-4, 2); 
    //62 is mass, next is density (only matters if you are chaning 
    // the density of the material while making calculations. Units are g/cm^3.
    // Sign of density sets the units used when exporting the table in LISE.
    // 2 is which column in the table to use (each is a different model). ATIMA-LS is 2.

    double inEn = 74; //Initial energy in MeV
    double distance = 25; //Distance in mm

    std::cout << "Energy loss from SRIM: " << eModelSrim.GetEnergyLoss(inEn, distance) << std::endl;
    std::cout << "Energy loss from LISE: " << eModelLise.GetEnergyLoss(inEn, distance) << std::endl;

    std::cout <<"Range from SRIM: " << eModelSrim.GetRange(inEn) << std::endl;
    std::cout <<"Range from LISE: " << eModelLise.GetRange(inEn) << std::endl;

    std::cout << "Energy after distance from SRIM: " << eModelSrim.GetEnergy(inEn, distance) << std::endl;
    std::cout << "Energy after distance from LISE: " << eModelLise.GetEnergy(inEn, distance) << std::endl;

}
