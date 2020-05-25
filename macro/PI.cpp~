#ifndef __CINT__
#include <iostream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH2.h>
#include <TText.h>
#include <TStyle.h>
#include <TArrow.h>
#include "TArtStoreManager.hh"
#include "TClonesArray.h"
#include "TArtEventStore.hh"
#include "TArtRawEventObject.hh"
#include "TArtRawFADCDataObject.hh"
#include "TArtESPRIParameters.hh"
#include "TArtUserParameters.hh"
#include "TArtCalibSRPPAC.hh"
#include "TArtCalibTDCHit.hh"
#include "TArtCalibPlas.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtSRPPAC.hh"
#include "TArtPlas.hh"
#include "TArtTDCHit.hh"
#include "signal.h"

#endif

bool stoploop = false;
void stop_interrupt(int sig){
  printf("keybprd interrupt\n");
  stoploop = true;
}
double funcZd(double x);
double funcZu(double x);
double funcAZd(double x);
double funcAZu(double x);

double func2Zd(double x);
double func2Zu(double x);
double func2AZd(double x);
double func2AZu(double x);
/*****const value*********/

//Trigger condition of f12N2 
const int tmin_f12n2=-8000; const int tmax_f12n2=-6000;//Sn&Ca
//const int qmin_f12n2=1000;  const int qmax_f12n2=2000;//Sn
const int qmin_f12n2=700;  const int qmax_f12n2=2500;//Ca

const int tmin_f7dia=8300; const int tmax_f7dia=10700;//Ca
const int tmin_f3dia=-1300; const int tmax_f3dia=1000;//Ca

int main(int argc, char * argv[]){
  std::string fname; Long64_t emax; Int_t flugpi=0;//if 1->make pi root file cp
  if(argc == 1){
    fname = "~/quser/work/enyo/rootfile/5056.root";
    //    emax = 5.E+4;
  }else if(argc == 2){
    fname = argv[1];
    // emax = -1;
  }else if(argc == 3){
    fname = argv[1];
    flugpi = atoi(argv[2]);
    // emax = atoi(argv[2]);
    std::cout << fname << std::endl;
  }else{
    std::cout << "Usage: ./reconst 0000.root"<<std::endl;
    return -1;
  }
  std::string outdir = fname.substr(9,2);
  std::string outfile = fname.substr(12);
  std::string runnumber = fname.substr(12,4);

  Int_t Nrun = std::atoi(runnumber.c_str());

  TChain *tespri = new TChain("rtree");
  TString strroot = Form("%s",fname.c_str());
  tespri->Add(strroot);

  //*******Define global variables in variable.h***********
  Double_t vx_f5[2][2];//BLD52(SRPPAC)[51 or 52][X or Y]//Ca run has only 52X  
  Double_t vxe12[2][2]={{-1111,-1111},{-1111,-1111}};//F12Xe raw [left or right] [tdc or adc]
  Double_t vxe12d[2][2]={{-1111,-1111},{-1111,-1111}};//F12XeDelay raw [left or right] [tdc or adc]
  Double_t vxe8[2][2]={{-1111,-1111},{-1111,-1111}};//F8Xe raw [left or right] [tdc or adc]
  Double_t vxe8d[2][2]={{-1111,-1111},{-1111,-1111}};//F8XeDelay raw [left or right] [tdc or adc]
  Double_t vn2[2][2]={{-1111,-999999},{-1111,-999999}};//F12N2 raw [left or right] [tdc or adc]
  Double_t vdia[2][4]={-999999,-999999};//Dia timing raw [f3 or f7][pad 1-4] 
  Double_t vdia_mean[2]={-999999,-999999};//Dia timing raw [f3 or f7] 

  //for analysis
  Double_t vtof_f37=-9999;
  Double_t vxe12_q=-1111;
  Double_t vxe8_q=-1111;
  Double_t vn2_q=-1111;

/******************************************************/

  //*********new TTree for analysis*************/
  TTree *tana = new TTree("tana","tana");

  tana->Branch("bldf5",&vx_f5,"bldf5[2][2]/D");//[BLD51 or 52][X or Y]
  tana->Branch("f12xe",&vxe12,"f12xe[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f12xed",&vxe12d,"f12xed[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f8xe",&vxe8,"f8xe[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f8xed",&vxe8d,"f8xed[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f12n2",&vn2,"f12n2[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("Dia",&vdia,"Dia[2][4]/D");//[f3Dia or f7Dia][pad 1-4]
  tana->Branch("Dia_t",&vdia_mean,"Dia_t[2]/D");//[f3Dia or f7Dia]
  /**************************************************************/

  Double_t cx1=9945.;Double_t cy1=1015.;
  Double_t rx1=56.;Double_t ry1=52.; 
  Double_t cx2=10000.;Double_t cy2=1015.;
  Double_t rx2=380.;Double_t ry2=300.; 
  if(Nrun>=3135){
    cx1=9945.; cy1=975.; rx1=56.; ry1=52.; 
    cx2=10000.; cy2=975.; rx2=380.; ry2=300.; 
  }
  TEllipse *cca48_3133 = new TEllipse(cx1,cy1,rx1,ry1,0.,360.,0);
  TEllipse *call_3133 = new TEllipse(cx2,cy2,rx2,ry2,0.,360.,0);

  /*******histgram********************************/
  Int_t npi=11;//number of hist for PI
  TH2F *hpi[npi];  
  hpi[0]=new TH2F("hpi_00","BLD52X vs f3f7 tof",100,0,100,200,9600,10400);
  hpi[1]=new TH2F("hpi_01","f3f7 tof vs f8XeQ",200,9600,10400,100,200,1500);
  hpi[2]=new TH2F("hpi_02","f3f7tofhsi vs f8XeQ",200,9600,10400,100,200,1500);
  hpi[3]=new TH2F("hpi_03","f3f7 tof vs f12XeQ",200,9600,10400,100,400,1500);
  hpi[4]=new TH2F("hpi_04","f3f7tofhsi vs f12XeQ",200,9600,10400,100,400,1500);
  hpi[5]=new TH2F("hpi_05","f3f7 tof vs f12N2Q",200,9600,10400,100,300,2000);
  hpi[6]=new TH2F("hpi_06","f3f7tofhsi vs f12N2Q",200,9600,10400,100,300,2000);
  hpi[7]=new TH2F("hpi_07","f3f7tofhsi vs f12XeQ GateIN",200,9600,10400,100,400,1500);
  hpi[8]=new TH2F("hpi_08","f3f7tofhsi vs f12XeQ",100,9600,10400,100,400,1500);
  hpi[9]=new TH2F("hpi_09","f3f7tofhsi vs f12XeQ Z=20",200,9600,10400,100,400,1500);
  hpi[10]=new TH2F("hpi_10","f3f7tofhsi vs f12XeQ A/Z=2.4",200,9600,10400,100,400,1500);

  Int_t nsep=2;//number of hist Z, A/Z
  TH1F *hsep[npi];  
  hsep[0]=new TH1F("hsep_00","f3f7tofhosei[A/Z] (Z=20)",200,9600,10400);
  hsep[1]=new TH1F("hsep_01","f12XeQ[Z] (A/Z=2.4)",250,500,1500);
  /**************************************************/

  TClonesArray *srppac_array=0;
  tespri->SetBranchAddress("ESPRISRPPAC",&srppac_array);
  TClonesArray *plas_array=0;
  tespri->SetBranchAddress("ESPRIPLAS",&plas_array);
  TClonesArray *tdc_array=0; 
  tespri->SetBranchAddress("ESPRITDC",&tdc_array);


  Int_t neve = tespri->GetEntries();
  int nbdcok=0; //number of the event of "traking ok" 
  TArtSRPPAC *hit_srppac = 0;
  TArtPlas *hit_pla = 0;
  TArtTDCHit *hit_tdc = 0;

  // std::string outputall = outdir+"/"+"allpi_"+outfile;
  // TFile *anafileall = new TFile(Form("anaresult/%s",outputall.c_str()),"RECREATE");

  /************* Fill data to new tree*******************/
  for(Int_t ne=0;ne<neve;ne++){
    Int_t trigN2[2]={0};//N2 scinti tdc & qdc ok [left&right]
    Int_t trigveto=1;//veto
    tespri->GetEntry(ne);
    Int_t nsrppac = srppac_array->GetEntries();
    Int_t npla = plas_array->GetEntries(); 
    Int_t ntdc = tdc_array->GetEntries();
    Int_t nmd[4]={0};//Dia multiplicity
    Int_t nmn2[2]={0};//N2 multiplicity
   
    for(Int_t j=0;j<2;j++){
      vn2[j][0]=-999999;
    }
    vdia_mean[0]=-999999;vdia_mean[1]=-999999;
    for(Int_t j=0;j<4;j++){
      vdia[0][j]=-999999;vdia[1][j]=-999999;
    }
  
   //    std::cout<<ne<<" "<<nnai<<std::endl;

    /*********  ESPRISRPPAC    *******************************/
    vx_f5[0][0]=-999;vx_f5[0][1]=-999;vx_f5[1][0]=-999;vx_f5[1][1]=-999;
    for(Int_t i=0;i<nsrppac;i++){
      TArtSRPPAC *hit_srppac = (TArtSRPPAC *)srppac_array->At(i);
      Int_t layer = hit_srppac->GetLayer();
      //      std::cout<<ne<<":"<<i<<":"<<layer<<":"<<hit_srppac->GetSRPPACX()<<std::endl;
     
      vx_f5[layer-1][0] = hit_srppac->GetSRPPACX();
      vx_f5[layer-1][1] = hit_srppac->GetSRPPACY();
    }

    /**********qdc from ESPRIPLAS          ****************/  
    for(Int_t i=0;i<2;i++){
      for(Int_t ii=0;ii<2;ii++){
	vxe12[i][ii]=-1111;vxe8[i][ii]=-1111;vn2[i][ii]=-1111;
      }
    }
    for(Int_t i=0;i<npla;i++){
      TArtPlas *hit_pla = (TArtPlas *)plas_array->At(i);
      Int_t id_plane = hit_pla->GetPlaneID();
      Int_t layer = hit_pla->GetLayer();
      Int_t channel = hit_pla->GetCh();
      
      if(layer==3 && id_plane==33){
	if(channel==1){//F12Xe left no Delay
	  vxe12[0][0] = (Double_t)hit_pla->GetTime();
	  vxe12[0][1] = (Double_t)hit_pla->GetEnergy();
	}else if(channel==2){//F12Xe right no Delay
	  vxe12d[0][0] = (Double_t)hit_pla->GetTime();
	  vxe12d[0][1] = (Double_t)hit_pla->GetEnergy();
	}else if(channel==3){//F12Xe right no Delay
	  vxe12[1][0] = (Double_t)hit_pla->GetTime();
	  vxe12[1][1] = (Double_t)hit_pla->GetEnergy();
	}else if(channel==4){//F12Xe right no Delay
	  vxe12d[1][0] = (Double_t)hit_pla->GetTime();
	  vxe12d[1][1] = (Double_t)hit_pla->GetEnergy();
	}
      }
      if(layer==5 && id_plane==35){
      	if(channel==1){//F8Xe left no Delay
	  vxe8[0][0] = (Double_t)hit_pla->GetTime();//no data->dummy
	  vxe8[0][1] = (Double_t)hit_pla->GetEnergy();
	}else if(channel==2){//F8Xe left no Delay
	  vxe8d[0][0] = (Double_t)hit_pla->GetTime();//no data->dummy
	  vxe8d[0][1] = (Double_t)hit_pla->GetEnergy();
      	}else if(channel==3){//F8Xe left no Delay
      	  vxe8[1][0] = (Double_t)hit_pla->GetTime();//no data->dummy
      	  vxe8[1][1] = (Double_t)hit_pla->GetEnergy();
      	}else if(channel==4){//F8Xe left no Delay
      	  vxe8d[1][0] = (Double_t)hit_pla->GetTime();//no data->dummy
      	  vxe8d[1][1] = (Double_t)hit_pla->GetEnergy();
      	}
      }
      if(layer==1 && id_plane==31){
      	if(channel==1){//N2 left no Delay
	  //      	  vn2[0][0] = (Double_t)hit_pla->GetTime();
      	  vn2[0][1] = (Double_t)hit_pla->GetEnergy();
      	}else if(channel==2){//N2 left no Delay
	  //      	  vn2[1][0] = (Double_t)hit_pla->GetTime();
      	  vn2[1][1]= (Double_t)hit_pla->GetEnergy();
      	}
      }

    }//end of the roop npla

    /**********fLTDC from ESPRITDC          ****************/  
    
    for(int i=0;i<ntdc;i++){
      TArtTDCHit *hit_tdc = (TArtTDCHit *)tdc_array->At(i);
      Int_t id_plane = hit_tdc->GetPlaneID();
      Int_t layer = hit_tdc->GetLayer();
      Int_t wireid = hit_tdc->GetWireID();
      Int_t toffset = hit_tdc->GetTzero();
      Int_t fltdc = hit_tdc->GetTDC();

      if(layer==1&&id_plane==31){
	for(Int_t j=0;j<2;j++){
	  if(wireid=j+1&&fltdc>tmin_f12n2&&fltdc<tmax_f12n2){
	    vn2[j][0]=fltdc;
	    //std::cout<<nmd[j]<<":"<<fltdc<<std::endl;
	  }
	}
      }
       

      if(layer==13&&id_plane==97){
      	for(Int_t j=0;j<4;j++){
	  if(wireid==j+1&&fltdc>tmin_f3dia&&fltdc<tmax_f3dia){
    	    vdia[0][j]=fltdc;
    	  }
      	} 
      }
      if(layer==17&&id_plane==98){
      	for(Int_t j=0;j<4;j++){
      	  if(wireid==j+1&&fltdc>tmin_f7dia&&fltdc<tmax_f7dia){
    	    vdia[1][j]=fltdc;
    	  }   
    	}  
      }
    }//end of the roop ntdc
    for(Int_t i=0;i<2;i++){
      if(vdia[i][0]>-999990&&vdia[i][1]>-999990
    	 &&vdia[i][2]>-999990&&vdia[i][3]>-999990){
    	vdia_mean[i]=(vdia[i][0]+vdia[i][1]+vdia[i][2]+vdia[i][3])/4.;
      }
    }

    /**************  analysis here   ******************************/
    vtof_f37=-9999;vxe12_q=-1111;vxe8_q=-1111;vn2_q=-1111;

    /** Trigger & cut setting **/
    for(Int_t i=0;i<2;i++){//N2 trig 
      if(vn2[i][0]>tmin_f12n2 && vn2[i][0]<tmax_f12n2){
	if(vn2[i][1]>qmin_f12n2 && vn2[i][1]<qmax_f12n2){
	  trigN2[i]=1;
	}
      }
    }
    /*-----------------------*/

    vtof_f37 = (vdia_mean[1]-vdia_mean[0]); //ToF between f3 & f7 Dia
    vxe12_q = (vxe12[0][1]-vxe12d[0][1]+vxe12[1][1]-vxe12d[1][1])/2.; //F12 Xe q
    vxe8_q = (vxe8[0][1]-vxe8d[0][1]+vxe8[1][1]-vxe8d[1][1])/2.; //F8 Xe q
    vn2_q = (vn2[0][1]+vn2[1][1])/2.; //F8 Xe q
       
    /*****************************************/
    //FILL to tree!!!
    tana->Fill();//event Fill
    
    Double_t tmpx; Double_t tmpy;
    //Fill to hist when trigger is OK
    if(trigN2[0]*trigN2[1]==1){
      
      hpi[0]->Fill(vx_f5[1][0],vtof_f37);
      hpi[1]->Fill(vtof_f37,vxe8_q);
      hpi[2]->Fill(vtof_f37-(vx_f5[1][0]-60)*4.5,vxe8_q);
      hpi[3]->Fill(vtof_f37,vxe12_q);
      hpi[4]->Fill(vtof_f37-(vx_f5[1][0]-60)*4.5,vxe12_q);
      hpi[5]->Fill(vtof_f37,vn2_q);
      hpi[6]->Fill(vtof_f37-(vx_f5[1][0]-60)*4.5,vn2_q);
      tmpx=vtof_f37-(vx_f5[1][0]-60)*4.5;
      tmpy=vxe12_q;
      
      if((tmpx-cx1)*(tmpx-cx1)/rx1/rx1+(tmpy-cy1)*(tmpy-cy1)/ry1/ry1-1.<=0){
	hpi[7]->Fill(tmpx,tmpy);
      }     
      if((tmpx-cx2)*(tmpx-cx2)/rx2/rx2+(tmpy-cy2)*(tmpy-cy2)/ry2/ry2-1.<=0){
	hpi[8]->Fill(tmpx,tmpy);
      }
      if(Nrun<3135){
	if(tmpy>funcZd(tmpx)&&tmpy<funcZu(tmpx)){
	  hpi[9]->Fill(tmpx,tmpy);
	  hsep[0]->Fill(tmpx);
	}
	if(tmpy>funcAZd(tmpx)&&tmpy<funcAZu(tmpx)){
	  hpi[10]->Fill(tmpx,tmpy);
	  hsep[1]->Fill(tmpy);
	}
      }else{	
	if(tmpy>func2Zd(tmpx)&&tmpy<func2Zu(tmpx)){
	  hpi[9]->Fill(tmpx,tmpy);
	  hsep[0]->Fill(tmpx);
	}
	if(tmpy>func2AZd(tmpx)&&tmpy<func2AZu(tmpx)){
	  hpi[10]->Fill(tmpx,tmpy);
	  hsep[1]->Fill(tmpy);
	}
      }
    } 
    if(ne%100==0){
      printf("Event:%10d/%10d\r"
	     ,ne,neve);
      fflush(stdout);    
    }
    // if(ne%100==50000){
    //   tana->Write();   
    // }
    
  }//end of ne(tespri) roop to neve & fill()
  // anafileall->Close();  
  
  /*************************************************************************/
  /*********test graph**************/
  
  Int_t ng = 20;
  Double_t g1x[ng], g1y[ng], g2x[ng], g2y[ng];
  Double_t g3x[ng], g3y[ng], g4x[ng], g4y[ng];
  for(Int_t i=0;i<ng;i++){
    g1x[i]=9600+40*i; g2x[i]=9600+40*i; g3x[i]=9600+40*i; g4x[i]=9600+40*i;
    if(Nrun<3135){  
      g1y[i]=funcZd(g1x[i]); g2y[i]=funcZu(g2x[i]);
      g3y[i]=funcAZd(g3x[i]); g4y[i]=funcAZu(g3x[i]);
    }else{
      g1y[i]=func2Zd(g1x[i]); g2y[i]=func2Zu(g2x[i]);
      g3y[i]=func2AZd(g3x[i]); g4y[i]=func2AZu(g3x[i]);
    } 
  }

  TGraph *g1 = new TGraph(ng,g1x,g1y); g1->SetLineColor(2);g1->SetLineWidth(2);
  TGraph *g2 = new TGraph(ng,g2x,g2y); g2->SetLineColor(2);g2->SetLineWidth(2);
  TGraph *g3 = new TGraph(ng,g3x,g3y); g3->SetLineColor(2);g3->SetLineWidth(2);
  TGraph *g4 = new TGraph(ng,g4x,g4y); g4->SetLineColor(2);g4->SetLineWidth(2);

  /************* make canvas ****************************/
  cca48_3133->SetFillStyle(0);
  call_3133->SetFillStyle(0);
  TCanvas *ccut = new TCanvas("cut48","cut48",1500,800);
  ccut->Divide(3,2);
  ccut->cd(1);  gPad->SetLogz();
  hpi[4]->Draw("colz");
  cca48_3133->Draw("same");
  ccut->cd(4);  gPad->SetLogz();
  hpi[4]->Draw("colz");
  call_3133->Draw("same");
  g1->Draw("LPLC"); g2->Draw("LPLC");g3->Draw("LPLC"); g4->Draw("LPLC"); 
  ccut->cd(2);  gPad->SetLogz(); hpi[9]->Draw("colz");
  ccut->cd(3);  gPad->SetLogz(); hpi[10]->Draw("colz");
  ccut->cd(5); hsep[0]->Draw();
  ccut->cd(6); hsep[1]->Draw();
  
  TCanvas *csep = new TCanvas("csep","separation",1500,800);
  csep->Divide(2,1);
  TText *t;
  char fpara[100];
  Double_t seprng_AZ[6]={9760.,9870.,9920.,9980.,10020.,10150.};
  Double_t xh_AZ[3]={9650,10000,10100}; 
  Double_t yh_AZ=hsep[0]->GetMaximum();  
  TF1 *ffit_AZ[3];
  for(Int_t i=0;i<3;i++){
    ffit_AZ[i]=new TF1(Form("fitAZ%i",i),"gaus",9600,10200);
  }
  TF1 *f3gaus_AZ=new TF1("f3gausAZ","fitAZ0+fitAZ1+fitAZ2",9600,10200);
  f3gaus_AZ->SetParameters(yh_AZ*0.8,(seprng_AZ[0]+seprng_AZ[1])/2.,20
			   ,yh_AZ*0.9,(seprng_AZ[2]+seprng_AZ[3])/2.,20
 			   ,yh_AZ*0.2,(seprng_AZ[4]+seprng_AZ[5])/2.,20);

  Double_t height_AZ[3];
  Double_t mean_AZ[3];
  Double_t sigma_AZ[3];
  hsep[0]->Fit("f3gausAZ","0","",seprng_AZ[0],seprng_AZ[5]);
  for(Int_t i=0;i<3;i++){
    height_AZ[i]=f3gaus_AZ->GetParameter(3*i); mean_AZ[i]=f3gaus_AZ->GetParameter(3*i+1);
    sigma_AZ[i]=f3gaus_AZ->GetParameter(3*i+2);
  }
  f3gaus_AZ->SetLineColor(46); 

  csep->cd(1);
  hsep[0]->SetStats(false);
  hsep[0]->Draw();
  f3gaus_AZ->Draw("SAME");

  for(Int_t i=0;i<3;i++){
    sprintf(fpara,"h:%6.2f+-%5.2f",height_AZ[i],f3gaus_AZ->GetParError(3*i));
    t = new TText(xh_AZ[i],height_AZ[i],fpara);  t->SetTextSize(0.03); t->Draw();
    sprintf(fpara,"m:%6.3f+-%5.3f",mean_AZ[i],f3gaus_AZ->GetParError(3*i+1));
    t = new TText(xh_AZ[i],height_AZ[i]-yh_AZ/20,fpara);  t->SetTextSize(0.03); t->Draw();
    sprintf(fpara,"sig:%6.3f+-%5.3f",sigma_AZ[i],f3gaus_AZ->GetParError(3*i+2));
    t = new TText(xh_AZ[i],height_AZ[i]-yh_AZ/10,fpara);  t->SetTextSize(0.03); t->Draw();
  }
  TArrow *ar_AZ[2];
  ar_AZ[0] = new TArrow(mean_AZ[0],height_AZ[1]/2.,mean_AZ[1],height_AZ[1]/2.,0.01,"<|>");
  ar_AZ[1] = new TArrow(mean_AZ[1],height_AZ[1]/2.,mean_AZ[2],height_AZ[1]/2.,0.01,"<|>");
  for(Int_t i=0;i<2;i++){
    ar_AZ[i]->SetLineColor(96);ar_AZ[i]->SetFillColor(96); ar_AZ[i]->Draw();
    sprintf(fpara,"%2.2fsig",fabsf(mean_AZ[2*i]-mean_AZ[1])/sigma_AZ[1]);
    t = new TText((mean_AZ[2*i]+mean_AZ[1])/2.+50*(i-1),height_AZ[1]/2.+yh_AZ/20,fpara);
    t->SetTextSize(0.035); t->Draw();
  }
  /***********************************************/
  Double_t seprng_Z[6]={880.,955.,965.,1055.,1070.,1160.};
  //  std::vector<double> seprng_Z{880.,955.,965.,1055.,1070.,1160.};
  if(Nrun>=3135){
    seprng_Z[0]=855.; seprng_Z[1]=915.; seprng_Z[2]=925.;
    seprng_Z[3]=1015.; seprng_Z[4]=1025.; seprng_Z[5]=1110.;
  }
  Double_t xh_Z[3]={600,1050,1160}; 
  Double_t yh_Z=hsep[0]->GetMaximum();  
  TF1 *ffit_Z[3];
  for(Int_t i=0;i<3;i++){
    ffit_Z[i]=new TF1(Form("fitZ%i",i),"gaus",500,1400);
  }
  TF1 *f3gaus_Z=new TF1("f3gausZ","fitZ0+fitZ1+fitZ2",500,1400);
  f3gaus_Z->SetParameters(yh_Z*0.5,(seprng_Z[0]+seprng_Z[1])/2.,15
			  ,yh_Z*0.9,(seprng_Z[2]+seprng_Z[3])/2.,15
			  ,yh_Z*0.5,(seprng_Z[4]+seprng_Z[5])/2.,15);
  Double_t height_Z[3];
  Double_t mean_Z[3];
  Double_t sigma_Z[3];
  hsep[1]->Fit("f3gausZ","0","",seprng_Z[0],seprng_Z[5]);
  for(Int_t i=0;i<3;i++){
    height_Z[i]=f3gaus_Z->GetParameter(3*i); mean_Z[i]=f3gaus_Z->GetParameter(3*i+1);
    sigma_Z[i]=f3gaus_Z->GetParameter(3*i+2);
  }
  f3gaus_Z->SetLineColor(56); 
  csep->cd(2);
  hsep[1]->SetStats(false);
  hsep[1]->Draw();
  f3gaus_Z->Draw("SAME");

  for(Int_t i=0;i<3;i++){
    sprintf(fpara,"h:%6.2f+-%5.2f",height_Z[i],f3gaus_Z->GetParError(3*i));
    t = new TText(xh_Z[i],height_Z[i],fpara);  t->SetTextSize(0.03); t->Draw();
    sprintf(fpara,"m:%6.3f+-%5.3f",mean_Z[i],f3gaus_Z->GetParError(3*i+1));
    t = new TText(xh_Z[i],height_Z[i]-yh_Z/20,fpara);  t->SetTextSize(0.03); t->Draw();
    sprintf(fpara,"sig:%6.3f+-%5.3f",sigma_Z[i],f3gaus_Z->GetParError(3*i+2));
    t = new TText(xh_Z[i],height_Z[i]-yh_Z/10,fpara);  t->SetTextSize(0.03); t->Draw();
  }
  TArrow *ar_Z[2];
  ar_Z[0] = new TArrow(mean_Z[0],height_Z[1]/2.,mean_Z[1],height_Z[1]/2.,0.01,"<|>");
  ar_Z[1] = new TArrow(mean_Z[1],height_Z[1]/2.,mean_Z[2],height_Z[1]/2.,0.01,"<|>");
  for(Int_t i=0;i<2;i++){
    ar_Z[i]->SetLineColor(96);ar_Z[i]->SetFillColor(96); ar_Z[i]->Draw();
    sprintf(fpara,"%2.2fsig",fabsf(mean_Z[2*i]-mean_Z[1])/sigma_Z[1]);
    t = new TText((mean_Z[2*i]+mean_Z[1])/2.+75*(i-1),height_Z[1]/2.+yh_Z/20,fpara);
    t->SetTextSize(0.035); t->Draw();
  }


 /**********************************************/
  std::string output = outdir+"/"+"pi2_"+outfile;
  TFile *anafile = new TFile(Form("anaresult/%s",output.c_str()),"RECREATE");

  tana->Write();
  for(Int_t i=0;i<npi;i++){
    hpi[i]->Write();
  } 
  for(Int_t i=0;i<nsep;i++){
    hsep[i]->Write();
  } 
  ccut->Write();
  csep->Write();
  anafile->Close(); 

  std::cout<<"input = "<<fname<<", evenumber = "<<neve<<std::endl;

  std::cout<<"analysis rootfile is created as anaresult/"<<output<<std::endl;
  return 0;
}

double funcZd(double x){
  double y = 2./11.*x-850;
  return y;
}
double funcZu(double x){
  double y = 2./11.*x-738;
  return y;
}
double funcAZd(double x){
  double y = -11.5/5.*x+23760;
  return y;
}
double funcAZu(double x){
  double y = -11.5/5.*x+24040;
  return y;
}
/**********************/
double func2Zd(double x){
  double y = 2./11.*x-850-30;
  return y;
}
double func2Zu(double x){
  double y = 2./11.*x-738-47;
  return y;
}
double func2AZd(double x){
  double y = -11.5/5.*x+23760-50;
  return y;
}
double func2AZu(double x){
  double y = -11.5/5.*x+24040-50;
  return y;
}
