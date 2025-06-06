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
#include "TF1.h"
#include "TString.h"
#include "runinfo.h"
using std::cout;
using std::endl;
using namespace info;
//simulation file that creates traces for beam or fusion events, all events convoluted with gaussian 100kev FWHM for music anode resolution 

//to run: ./anodes runnb
int main(int argc,char **argv){
//NOTE::: NEED TO ADD IN DIRECTION OF PARTICLE   
 
    //setting run number and associated variables - see runinfo.cpp for beam, energy, gas, ect.  
    int runnb = 2448;
    if (argc==2){runnb=atoi(argv[1]);} 
    runinfo run; 
    run.setinfo(runnb); 
  
    TApplication *app = new TApplication("app",&argc,argv);//application call to allow figures to pop up in excecutable 
   
    //defniing gas to lose energy in  
    catima::Material gas;
    for (int i=0; i<run.gasel; i++){
	gas.add_element(run.gasA[i],run.gasZ[i],run.gasS[i]);
    }
    gas.density(run.density);
    //thickness defined later bc it changes 


//**********************************
//*******NOTE: CHANGED THICKNES TO 0.25CM 
//*********************************
   
    catima::Projectile proj(run.A,run.Z);
 
    cout <<Form("A=%i, Z=%i on A=%i, Z=%i\n", run.A,run.Z,run.gasA[0],run.gasZ[0]);
    //random variable for gaussian later
    TRandom3 ran;
    //connecting app to the canvas. When window is closed exe will terminate 
   // TCanvas *c1 = new TCanvas();
    //TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
    //rc->Connect("CloseWindow()","TApplication",app,"Terminate()");
    
    float x[20]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};

    float E_out[21]={0};
    //defining a 1.5 anode thickness to represent snout of gas beam travels through
    //energy after going through this is the energy found in anode 0 
    gas.thickness(run.density*1.8);
    float Estart = run.Ein/run.A;
    auto init = catima::calculate(proj,gas,Estart); 
    E_out[0]=Estart - init.Eloss/run.A;
    gas.thickness(run.density*0.25);//put thickness back to 1/5 of an anode 
        
    float Eloss[20]={0};
//******************************************** 
//*******begin loop over anodes *****************
//******************************************* 
    for(int i=0; i<20; i++){
	float Etemp[6];
	Etemp[6]=0;
 	float los_temp[6]={0};
	los_temp[0]=E_out[i];
	//start loop splitting anode into five parts      
	for (int j=0; j<5; j++){
             auto result = catima::calculate(proj,gas,los_temp[j]);  
	     Etemp[j]=result.Eloss;//not conveluted w gaussian yet bc we are comparing to exp data 
	     Etemp[6]+=Etemp[j];
             los_temp[j+1]=los_temp[j]-(Etemp[j]/20);
	}
	Eloss[i]=Etemp[6];
        E_out[i+1]=E_out[i]-(Eloss[i]/20); 
     } //end loop over anodes

       
 
    cout << "done " << endl;
    TString filename = Form("../../../Ana_MSU/MSU_BeamRun%i.root",runnb);
    TFile *input = new TFile(filename.Data());
    cout << "got Tfile " << endl;
    TList *expdata = (TList*)gROOT->FindObject("SummedAnodesCalib");
    
    float scales[20];
    float fwhm[20];
    std::vector< TCanvas*> canvases(20);
    std::vector<TH1F*> graphs(20);
//********************************
//******START LOOP OVER ANODES TO CANVASES
//*******************************
    for (int j=0;j<20;j++){
	TH1F *an0 = (TH1F*)expdata->At(j);
	cout << "got th1F" << endl;
    
	TCanvas *c0 = new TCanvas();
	c0->cd();
	canvases[j]=c0;
	Int_t binmax = an0->GetMaximumBin();
	Double_t x_max = an0->GetBinCenter(binmax);
   
	auto fit = new TF1("fit","gaus",x_max-0.15,x_max+0.15);
	if (j==8 || j==13 || j==17 || j==19){ 
    	    an0->Fit(fit,"R");} 
        else{an0->Fit(fit,"","",x_max-.1,x_max+.1);}
        Double_t mean = fit->GetParameter(1);
        cout << mean << " mean " << endl;
        an0->Draw("same");


        TH1F *h = new TH1F("h","h",550,-0.5,5);
        graphs[j]=h;
        int fill = fit->GetParameter(0)*11;
        cout << fill << std::endl;
        scales[j] = mean/Eloss[j];
        fwhm[j]=fit->GetParameter(2)*2.355; 
        for (int i=1; i<fill;i++){

           h->Fill(ran.Gaus(Eloss[j]*scales[j],(0.1/2.355)));
         }

        h->SetLineColor(kRed);
  
        h->Draw("same");
    }//end loop over anodes
    TCanvas *c2 = new TCanvas();
    c2->cd();
    TGraph *g1 = new TGraph(20,x,scales);
    g1->SetTitle("Scaling factor per anoded");
    g1->Draw();
    TCanvas *c3 = new TCanvas();
    c3->cd();
    TGraph *g2 = new TGraph(20,x,fwhm);
    g2->SetTitle("FWHM per anode");
    g2->Draw();
    cout << "scaling factors per anode : "; 
    for (int i=0;i<19;i++){cout << scales[i] << ",";}
    cout << scales[19] << endl;
    app->Run(kTRUE); //canvas main menu ->File ->QUIT ROOT (I wish there was a better way of doing this but ig its fine
    delete app;     
    input->Close();
    //app->Terminate();    
    return 0;
}
