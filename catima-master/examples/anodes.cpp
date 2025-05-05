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
using std::cout;
using std::endl;
//simulation file that creates traces for beam or fusion events, all events convoluted with gaussian 100kev FWHM for music anode resolution 

//to run: ./music_dedx (anode #)
int main(int argc,char **argv){
//NOTE::: NEED TO ADD IN DIRECTION OF PARTICLE    
  
    int an_hist = 0;
    if (argc==2){an_hist=atoi(argv[1]);} 
    TApplication *app = new TApplication("app",&argc,argv);//application call to allow figures to pop up in excecutable 
 
//inistilize basic gas properties for helium run  
    catima::Material he;
    he.add_element(4,2,1);
    he.density(1.3e-4); //density in g/cm^3
    he.thickness(1.625e-4); //thickness in g/cm^2

//initlize basic gas properties for methane run 
    catima::Material methane;
    methane.add_element(1,1,4);
    methane.add_element(12,6,1);
    methane.density(9.6534e-5);
    methane.thickness(9.6534e-5*1.25);
   
//defining commonly used beams 
    catima::Projectile S30(30,16); // define projectile, ie 30S
    catima::Projectile O20(20,8);

  //  cout<<"30S->He\n";
    cout <<"20O->Methane\n";
//random variable for gaussian later
   TRandom3 ran;
//connecting app to the canvas. When window is closed exe will terminate 
    TCanvas *c1 = new TCanvas();
    TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
    rc->Connect("CloseWindow()","TApplication",app,"Terminate()");
    
    float x[20]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    float Eloss_test;
    int first=0; 


	float E_out[21]={0};
   	E_out[0]=2.4;//should be 2.7 for oxygen but we will scale to 2.055
        float Eloss[20]={0};
        float E_sig[20];
        float ang_sig[20];      
 
   	 for(int i=0; i<20; i++){
             auto result = catima::calculate(O20,methane,E_out[i]); 
       	 //    Eloss[i]=ran.Gaus(result.Eloss,0.1/2.355);
      	    Eloss[i]=result.Eloss;
         
	    E_out[i+1]=E_out[i]-(Eloss[i]/20); 
            E_sig[i]=result.sigma_E;
            ang_sig[i]=result.sigma_a;
	} //end loop over anodes

       
 
    cout << "done " << endl;
    TFile *input = TFile::Open("../../../Ana_MSU/MSU_Run2451.root");
    cout << "got Tfile " << endl;
    TList *expdata = (TList*)gROOT->FindObject("SummedAnodesCalib");
    float scales[20];
    float fwhm[20];
    std::vector< TCanvas*> canvases(20);
    std::vector<TH1F*> graphs(20);
    for (int j=0;j<3;j++){
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
   for (int i=0;i<19;i++){cout << scales[i] << ",";}
   cout << scales[19] << endl;
    app->Run(kTRUE); //canvas main menu ->File ->QUIT ROOT (I wish there was a better way of doing this but ig its fine
    delete app;     
    input->Close();
    //app->Terminate();    
    return 0;
}
