#include "AtELossTable.h"
// IWYU pragma: no_include <ext/alloc_traits.h>
#include "AtStringManip.h"//

//#include <FairLogger.h>

#include <algorithm> // for max
#include <cmath>     // for fabs
#include <iostream>

namespace AtTools {

void AtELossTable::LoadTable(const std::vector<double> &energy, const std::vector<double> &dEdX)
{

   if (energy.size() > 0)
      std::cout << "Loading dE/dX from " << energy.front() << " MeV to " << energy.back() << " MeV " << std::endl;
   else
      std::cerr << "Passed table was empty!"<< std::endl;

   std::vector<double> dXdE;
   for (auto &elem : dEdX)
      dXdE.push_back(1 / elem);
   fdXdE = tk::spline(energy, dXdE);
}

AtELossTable::AtELossTable(const std::vector<double> &energy, const std::vector<double> &dEdX, double density)
   : AtELossModel(density)
{
   LoadTable(energy, dEdX);
}

double AtELossTable::GetdEdx(double energy) const
{
   double dedx = 1 / fdXdE(energy);
   return dedx * fdEdxScale;
}

double AtELossTable::GetRange(double energyIni, double energyFin) const
{
   if (energyIni == energyFin)
      return 0;
   if (energyFin < fdXdE.get_x_min()) {
      // LOG(debug) << "Attemting to integrate energy to " << energyFin << " when min energy in table is "
      //            << fdXdE.get_x_min()<< std::endl;
      energyFin = fdXdE.get_x_min();
   }
   return fdXdE.integrate(energyFin, energyIni);
}

double AtELossTable::GetEnergy(double energyIni, double distance) const
{
   if (distance == 0)
      return energyIni;
   if (energyIni < 1e-6 || GetRange(energyIni) < distance)
      return 0.;

   int maxIt = 100;

   double guessEnergy = energyIni - GetdEdx(energyIni) * distance;
   for (int i = 0; i < maxIt; ++i) {

      double range = GetRange(energyIni, guessEnergy);
      if (fabs(range - distance) < fDistErr) {
         // LOG(debug) << "Energy converged in " << i + 1 << " iterations.";
         return guessEnergy;
      }
      // LOG(debug) << "guessE: " << guessEnergy << " d: " << distance << " R: " << range
      //            << " dEdX: " << GetdEdx(guessEnergy) << " dEnergy: " << GetdEdx(guessEnergy) * (distance - range)<< std::endl;
      // Update guess energy using what is essentially Newton's method
      guessEnergy += GetdEdx(guessEnergy) * (range - distance);
   }

   std::cerr << "Energy calculation (" << energyIni << " MeV through " << distance << " mm) failed to converge in "
              << maxIt << " iterations!"<< std::endl;
   return -1;
}

double AtELossTable::GetEnergyOld(double energyIni, double distance) const
{
   if (distance == 0)
      return energyIni;
   if (energyIni < 1e-6)
      return 0.;

   double maxRange = GetRange(energyIni);

   if (distance > maxRange)
      return 0.;

   // Searcg for energyFinal where GetRange(energyIni, energyFinal) == distance

   double lowGuess = 0;
   double highGuess = energyIni;
   double guessEnergy = energyIni / 2.0;
   double range = GetRange(energyIni, guessEnergy);
   double distErr = 1.0e-4; // in mm
   int numIt = 1;

   while (fabs(range - distance) > distErr) {
      numIt++;
      // We went too far
      if (range > distance) {
         lowGuess = guessEnergy;
         guessEnergy = (lowGuess + highGuess) / 2.0;
      } else {
         // We didn't go far enough
         highGuess = guessEnergy;
         guessEnergy = (lowGuess + highGuess) / 2.0;
      }
      range = GetRange(energyIni, guessEnergy);
   }

   std::cout << "Energy converged in " << numIt << " iterations."<< std::endl;
   return guessEnergy;
}

void AtELossTable::LoadSrimTable(std::string fileName)
{
   std::ifstream file(fileName);
   if (!file.is_open()) {
      std::cerr << "Failed to open SRIM file " << fileName<< std::endl;
      exit(1);
   }

   bool atConversion = false;
   double conversion = 0;
   std::vector<double> energy;
   std::vector<double> dEdX;

   while (!file.eof()) {
      // Get the current line
      std::string line;
      std::getline(file, line);
      auto tokens = SplitString(line, ' ');
      // LOG(debug) << "Processing line " << line<< std::endl;

      // If this is the densityt line, grab it and continue
      if (tokens[0] == "Target" && tokens[1] == "Density") {
         std::cout << "Setting target density to: " << tokens[3] << " g/cm^3"<< std::endl;
         SetIniDensity(std::stod(tokens[3]));
         continue;
      }

      // Check if we've reached the conversion part of the table
      atConversion |= (tokens[0] == "Multiply");
      // LOG(debug) << "At converion: " << atConversion;

      // If we're at the conversion stage, look for the correct conversion factor
      try {
         if (atConversion && tokens.at(1) == "MeV" && tokens.at(3) == "mm") {
            conversion = std::stod(tokens.at(0));
            std::cout << "Using conversion factor of " << conversion <<std::endl;
            break;
         }

         //  If this is part of the table, grab it and continue
         double en = std::stod(tokens.at(0)) * GetUnitConversion(tokens.at(1));
         double dedx = std::stod(tokens.at(2)) + std::stod(tokens.at(3));

         // If we made it this far then we have a valid table entry
         energy.push_back(en);
         dEdX.push_back(dedx);

         // LOG(debug) << en << " " << dedx<< std::endl;

      } catch (...) {
      }
   } // end loop over file

   for (auto &dedx : dEdX)
      dedx *= conversion;

   LoadTable(energy, dEdX);
}

void AtELossTable::LoadLiseTable(std::string fileName, double mass, double density, int column)
{
   std::ifstream file(fileName);
   if (!file.is_open()) {
      std::cerr << "Failed to open SRIM file " << fileName<< std::endl;
      exit(1);
   }

   std::vector<double> energy;
   std::vector<double> dEdX;
   SetIniDensity(fabs(density));

   while (!file.eof()) {
      // Get the current line
      std::string line;
      std::getline(file, line);
      auto tokens = SplitString(line, '\t');
      // LOG(debug) << "Processing line " << line;

      try {
         //  If this is part of the table, grab it and continue
         double en = std::stod(tokens.at(0 + column * 2)) * mass;
         double dedx = std::stod(tokens.at(1 + column * 2));

         if (density > 0) // dEdX is in units MeV/(mg/cm^2)
            dedx *= (100 * density);//JC change from 100 to 1000, convert ddensity from g/cm^3 to mg/cm^3 to match units of dedx
         // dEdX is in units MeV/um *JC change - MeV/cm
            dedx *= 1000; //changed from 1000 to .1 to convert to MeV/mm

         // LOG(debug) << "En: " << en << " dEdX " << dedx;
         // If we made it this far then we have a valid table entry
         energy.push_back(en);
         dEdX.push_back(dedx);
      } catch (...) {
      }
   }

   LoadTable(energy, dEdX);
}

double AtELossTable::GetUnitConversion(const std::string &unit)
{
   if (unit == "eV")
      return 1e-6;
   if (unit == "keV")
      return 1e-3;
   if (unit == "MeV")
      return 1e-0;
   if (unit == "GeV")
      return 1e3;
   return 0;
}
} // namespace AtTools
