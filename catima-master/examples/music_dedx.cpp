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
using std::cout;
using std::endl;
//simulation file that creates traces for beam or fusion events, all events convoluted with gaussian 100kev FWHM for music anode resolution 

//to run: ./music_dedx (anode #)
int main(int argc,char **argv){
//NOTE::: NEED TO ADD IN DIRECTION OF PARTICLE    
    if (argc<2){cout << "not enough inputs"; return 1;}
 
    TApplication *app = new TApplication("app",&argc,argv);//application call to allow figures to pop up in excecutable 

    int anode = atoi(argv[1]);//input for what anode fusion starts in, set to 20 if you want a beam event 
    double scales[20]={1.07202,1.03719,1.02841,1.01483,1.02579,1.01567,1.03889,1.0012,1.00787,1.0362,0.985296,1.0186,1.00916,1.008,1.0313,1.00856,0.993695,1.0023,1.01465,0.902206};
//update this to change from 30S to 20O
    TFile *fusion_input = TFile::Open(Form("../../../gemini/20O_fusion%i.root",anode),"READ");//opens gemini output file

    //intilize tree  
    TTree *t1 = new TTree();
    t1 = (TTree*)fusion_input->Get("GeminiEvent");
    
    Int_t mult, A[5], Z[5]; 
    float px[5], py[5], pz[5], E[5];
    
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
//loop over gemini entries, then make beam like trace up to anode, then loop over number of particles post reaction for rest of the anodes
    for(int j=0; j<1000; j++){	

	t1->GetEntry(j);
	float E_out[21]={0};
   	E_out[0]=2.055;//should be 2.7 for oxygen but we will scale
        float Eloss[20]={0};
        float E_sig[20];
        float ang_sig[20];      
 
   	 for(int i=0; i<anode; i++){
             auto result = catima::calculate(O20,methane,E_out[i]); 
       	     Eloss[i]=ran.Gaus(result.Eloss,0.1/2.355)*scales[i];
      	    // Eloss[i]=result.Eloss;
	    E_out[i+1]=E_out[i]-(Eloss[i]/20);
            E_sig[i]=result.sigma_E;
            ang_sig[i]=result.sigma_a;
	} //end loop over anodes before fusion

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
	   methane.thickness(9.6534e-5*hf);
	   float Etemp; 
           for (int i=anode;i<20;i++){
	     auto result = catima::calculate(ER,methane,E_out[i]);  
	     Etemp = ran.Gaus(result.Eloss,0.1/2.355)*1; 
       	     Eloss[i]+= Etemp;
	     E_out[i+1]=(E_out[i]-(Etemp/A[k]));	   
            E_sig[i]=result.sigma_E;
            ang_sig[i]=result.sigma_a;
	   // cout << "finished loop " << i <<  " with e out " << E_out[i+1] << endl;
 
    	   }//end loop over anodes 

	}//end of loop over mult  

  
	 if (Eloss[anode+1]==0){continue;}
         TGraph *trace = new TGraph(20,x,Eloss);
         trace->GetYaxis()->SetRangeUser(0,5);
         if (first==0){trace->Draw();first=1;} 
         else{ trace->Draw("same");}
	 Eloss_test=Eloss[0];
    }//end of loop over gemini events


    fusion_input->Close();
    cout << "done " << endl;
   app->Run(kTRUE); //canvas main menu ->File ->QUIT ROOT (I wish there was a better way of doing this but ig its fine
    delete app;     
    //app->Terminate();    
    return 0;
}
