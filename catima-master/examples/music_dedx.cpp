#include "catima/catima.h"
#include <iostream>
#include <cmath>
#include <TROOT.h>
#include<TCanvas.h>
#include <TGraph.h>
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TFile.h"
#include "THStack.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include "runinfo.h" 

using namespace info;
using std::cout;
using std::endl;
//simulation file that creates traces for beam or fusion events, all events convoluted with gaussian 100kev FWHM for music anode resolution 

//to run: ./music_dedx (anode #)
int main(int argc,char **argv){

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Initillize variables and open files 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (argc<2){cout << "not enough inputs"; return 1;}
//create application to view pop ups in exe
    TApplication *app = new TApplication("app",&argc,argv);
//read in anode input, modeling fusion in this anode
    int anode = atoi(argv[1]);
//set "run number" which will be used to define beam and gas settings (see runinfo.cpp), default to 30S
//quick guide: 2447-2451: 20O, 2184-2185: 18O, -1: 30S prediction 
    int runnb = -1;
    if (argc >2){ runnb = atoi(argv[2]);}
    runinfo run; 
    run.setinfo(runnb);
//create output file 
    std::string fileout = Form("Data/%s/%sSim_Data%i.csv",run.lab,run.lab,anode);  
    std::ofstream cfile;
    cfile.open(fileout);
//open gemini output file to use as input
    TFile *fusion_input = TFile::Open(Form("../../../gemini/Data/%s/%s_fusion%i.root",run.lab,run.lab,anode),"READ");
//intilize tree  
    TTree *t1 = new TTree();
    t1 = (TTree*)fusion_input->Get("GeminiEvent");     
    Int_t mult, A[5], Z[5]; 
    float px[5], py[5], pz[5], E[5];
//set branches  
    t1->SetBranchAddress("mult",&mult);
    t1->SetBranchAddress("px",&px);
    t1->SetBranchAddress("py",&py);
    t1->SetBranchAddress("pz",&pz);
    t1->SetBranchAddress("A",&A);
    t1->SetBranchAddress("Z",&Z);
    t1->SetBranchAddress("E",&E);
    //might need to add theta cm 
    int entries = t1->GetEntries();
    cout << entries << " Gemini events found " << endl;   
//set basic gas and projectile properties
    catima::Material gas;
    for (int i=0; i<run.gasel;i++){gas.add_element(run.gasA[i],run.gasZ[i],run.gasS[i]);}
    gas.density(run.density);
    gas.thickness(run.density*1.25);//might need to change this to 0.25? 
   
    catima::Projectile proj(run.A,run.Z); 
    cout <<Form("A=%i, Z=%i, on A=%i Z=%i \n", run.A, run.Z, run.gasA[0], run.gasZ[0]);
//set scaling factor and fwhm for each anode 
    std::ifstream scales("scales.csv"); 
   // scales.open("scales.csv");
    float scale_fac[20];
    float fwhm[20];
    std::vector<std::string> row;
    std::string line, word, temp;
    while(std::getline(scales,line)){

	row.clear();
	std::stringstream in(line);
	while(std::getline(in,word,',')){
	    row.push_back(word);
	}
	int A_text = std::stoi(row[0]);
	int Z_text = std::stoi(row[1]);
 
	if (A_text == run.A && Z_text == run.Z){
	    for (int i=0; i< 20; i++){
		fwhm[i]=std::stof(row[i+22]);
		scale_fac[i]=std::stof(row[i+2]);	
	    }	
	}
	
    }
scales.close();
   
//random variable for gaussian later
   TRandom3 ran;
//connecting app to the canvas. When window is closed exe will terminate 
    TCanvas *c1 = new TCanvas();
    TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
    rc->Connect("CloseWindow()","TApplication",app,"Terminate()");
//array on anode numbers used for plotting  
    float x[20]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Begin energy loss calculations
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    int first=0; //used for plotting purposes
    for(int j=0; j<entries; j++){	
//check for empty particles in gemini output 
	  
	t1->GetEntry(j);
	float E_out[21]={0};
	float Estart = run.Ein/run.A;
//calculate energy of beam after it travels 18 mm (length of snout before anode 0) 
	auto init = catima::calculate(proj, gas, Estart);
   	E_out[0]= Estart - (init.Eloss/run.A);
        float Eloss[20]={0};
//%%%%%%%Begin eloss prior to fusion %%%%%%%%%%%%%%%%%%%%%%% 
   	 for(int i=0; i<anode; i++){
	     auto result = catima::calculate(proj,gas,E_out[i]);
	       
       	     Eloss[i]=ran.Gaus(result.Eloss,fwhm[i]/2.355); 
	     E_out[i+1]=E_out[i]-(Eloss[i]/run.A);	
	} //end loop over anodes before fusion
//%%%%%%%%%%%Begin eloss for ERs after fusion %%%%%%%%%%%%%%%%%%%%%%%%
        cout << "mult: " << mult << endl;
	for (int k=0;k<mult;k++){
	    if (A[k]<5 || Z[k]==0)continue;//skipping empty particles and neutrons bc theres no eloss 
	   catima::Projectile ER(A[k],Z[k]);
	   cout << "projectile made z=" << Z[k] << " A=" << A[k] << endl;
	   E_out[anode]=E[k]/A[k]; 
	   float theta = atan2(py[k],pz[k]);
	   float h = 1.25/cos(theta);
	   float phi = atan2(px[k],pz[k]);
	   float hf = h/cos(phi);
	   cout << "final length " << hf << " cm" << endl;
	   gas.thickness(8.7759e-5*hf);
	   float Etemp; 
           for (int i=anode;i<20;i++){
	     auto result = catima::calculate(ER,gas,E_out[i]);  
	     Etemp = ran.Gaus(result.Eloss,fwhm[i]/2.355)*scale_fac[i];
       	     Eloss[i]+= Etemp;
	     E_out[i+1]=(E_out[i]-(Etemp/A[k]));	   
 
    	   }//end loop over anodes 

	}//end of loop over mult  
	 if (Eloss[anode+1]==0){continue;}
	 if (Eloss[0]<0.85){continue;}
 	 for (int i=0; i<20; i++){
	     Eloss[i]= Eloss[i];//*scale_fac[i]; 
	    
	 }	
	 
         TGraph *trace = new TGraph(20,x,Eloss);
         trace->GetYaxis()->SetRangeUser(0,5);
         if (first==0){trace->Draw();first=1;} 
         else{ trace->Draw("same");}
	 

        for (int i=0; i<20; i++){
	    if (i==19){cfile << Eloss[i] << endl;}
	    else{cfile << Eloss[i] << ", ";}
	}
    }//end of loop over gemini events


    fusion_input->Close();
    cout << "done " << endl;
    cfile.close();
   app->Run(kTRUE); //canvas main menu ->File ->QUIT ROOT (I wish there was a better way of doing this but ig its fine
    delete app;     
    //app->Terminate();    
    return 0;
}
