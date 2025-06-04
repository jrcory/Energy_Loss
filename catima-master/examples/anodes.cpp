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
using std::cout;
using std::endl;
//simulation file that creates traces for beam or fusion events, all events convoluted with gaussian 100kev FWHM for music anode resolution 

//to run: ./music_dedx (anode #)
int main(int argc,char **argv){
//NOTE::: NEED TO ADD IN DIRECTION OF PARTICLE    
  
    int runnb = 2448;
    if (argc==2){runnb=atoi(argv[1]);} 

    double den = 1.3e-4;
    if (argc==3){den==0; }//characterize conversion between pressure in torr to densit  
    TApplication *app = new TApplication("app",&argc,argv);//application call to allow figures to pop up in excecutable 
 
//inistilize basic gas properties for helium run  
    catima::Material he;
    he.add_element(4,2,1);
    he.density(1.3e-4); //density in g/cm^3
    he.thickness(1.625e-4); //thickness in g/cm^2
    
    float p150 = 1.3e-4;
    float p100 = 8.7759e-5;
//initlize basic gas properties for methane run 
    catima::Material methane;
    methane.add_element(1,1,4);
    methane.add_element(12,6,1);
    methane.density(p150);
    methane.thickness(p150*0.25);
//**********************************
//*******NOTE: CHANGED THICKNES TO 0.25CM 
//*********************************
   
//defining commonly used beams 
    catima::Projectile S30(30,16); // define projectile, ie 30S
    catima::Projectile O20(20,8);

  //  cout<<"30S->He\n";
    cout <<"20O->Methane\n";
//random variable for gaussian later
   TRandom3 ran;
//connecting app to the canvas. When window is closed exe will terminate 
   // TCanvas *c1 = new TCanvas();
    //TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
    //rc->Connect("CloseWindow()","TApplication",app,"Terminate()");
    
    float x[20]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    float Eloss_test;
    int first=0; 


	float E_out[21]={0};
   	E_out[0]=2.95;//should be 2.7 for oxygen but we will scale to 2.055
        float Eloss[20]={0};
      //  float E_sig[20];
       // float ang_sig[20];     
//******************************************** 
//*******begin loop over anodes *****************
//******************************************* 
   	 for(int i=0; i<20; i++){
	     float Etemp[6];
	     Etemp[6]=0;
 	     float los_temp[6]={0};
	     los_temp[0]=E_out[i];
	     cout << Etemp[6] << endl;
	     for (int j=0; j<5; j++){
             auto result = catima::calculate(O20,methane,los_temp[j]); 
 	     //Etemp[j]=ran.Gaus(result.Eloss,0.1/2.355);
	     Etemp[j]=result.Eloss;//changed for now bc we are MCing
	     Etemp[6]+=Etemp[j];
             los_temp[j+1]=los_temp[j]-(Etemp[j]/20);
	     }
       	 //    Eloss[i]=ran.Gaus(result.Eloss,0.1/2.355);
      	    //Eloss[i]=result.Eloss;
	    Eloss[i]=Etemp[6];
            cout << "Eloss " << Eloss[i] << endl; 
	    E_out[i+1]=E_out[i]-(Eloss[i]/20); 
           // E_sig[i]=result.sigma_E;
            //ang_sig[i]=result.sigma_a;
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
