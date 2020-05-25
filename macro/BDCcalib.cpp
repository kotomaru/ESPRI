#ifndef __CINT__
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TString.h>
#include "TArtStoreManager.hh"
#include "TClonesArray.h"
#include "TArtCalibTDCHit.hh"
#include "TArtCalibRDC.hh"
#include "TArtCalibBDC.hh"
#include "TArtPlas.hh"
#include "TArtBDC.hh"
#include "TArtTDCHit.hh"
#include "TArtDataObject.hh"
#include "signal.h"

#endif


bool stoploop = false;
void stop_interrupt(int sig){
  printf("keybprd interrupt\n");
  stoploop = true;
}

int main(int argc, char * argv[]){
  std::string fname; Long64_t emax;
  if(argc == 1){
    fname = "~/quser/work/enyo/rootfile/5056.root";
    //    emax = 5.E+4;
  }else if(argc == 2){
    fname = argv[1];
    // emax = -1;
  }else if(argc == 3){
    fname = argv[1];
    // emax = atoi(argv[2]);
    std::cout << fname << std::endl;
  }else{
    std::cout << "Usage: blPI 0000.root"<<std::endl;
    return -1;
  }
  std::string outdir = fname.substr(9,2);
  std::string outfile = fname.substr(12);

  TChain *tespri = new TChain("rtree");
  TString strroot = Form("%s",fname.c_str());
  tespri->Add(strroot);

  Int_t tdcval[16]={-1111};
  Int_t nh[16]={0};//reference 1~16 BDC1X->Y->BDC2X->Y
  Double_t parcst[16]={0};//fit parameter : const,mean,sigma
  Double_t parmean[16]={0};
  Double_t parsigma[16]={0};
  
  //*********new TTree for TDC of BDC*************/
  TTree *tree = new TTree("tree","tree");
  tree->Branch("historder",&nh,"historder[16]/I");
  tree->Branch("ParaConst",&parcst,"ParaConst[16]/D");
  tree->Branch("ParaMean",&parmean,"ParaMean[16]/D");
  tree->Branch("ParaSigma",&parsigma,"ParaSigma[16]/D");
  /*******histgram********************************/
  TH1F *h[16];
  Int_t tmin=-2E+3;Int_t tmax=2E+3;Int_t tbin=100;
  h[0] = new TH1F("h01","BDC1 X01-16",tbin,tmin,tmax);
  h[1] = new TH1F("h02","BDC1 X17-32",tbin,tmin,tmax);
  h[2] = new TH1F("h03","BDC1 X19-48",tbin,tmin,tmax);
  h[3] = new TH1F("h04","BDC1 X49-64",tbin,tmin,tmax);
  h[4] = new TH1F("h05","BDC1 Y01-16",tbin,tmin,tmax);
  h[5] = new TH1F("h06","BDC1 Y17-32",tbin,tmin,tmax);
  h[6] = new TH1F("h07","BDC1 Y19-48",tbin,tmin,tmax);
  h[7] = new TH1F("h08","BDC1 Y49-64",tbin,tmin,tmax);
  h[8] = new TH1F("h09","BDC2 X01-16",tbin,tmin,tmax);
  h[9] = new TH1F("h10","BDC2 X17-32",tbin,tmin,tmax);
  h[10] = new TH1F("h11","BDC2 X19-48",tbin,tmin,tmax);
  h[11] = new TH1F("h12","BDC2 X49-64",tbin,tmin,tmax);
  h[12] = new TH1F("h13","BDC2 Y01-16",tbin,tmin,tmax);
  h[13] = new TH1F("h14","BDC2 Y17-32",tbin,tmin,tmax);
  h[14] = new TH1F("h15","BDC2 Y19-48",tbin,tmin,tmax);
  h[15] = new TH1F("h16","BDC2 Y49-64",tbin,tmin,tmax);
  /**************************************************/
  TClonesArray *tdc_array=0;				\
  tespri->SetBranchAddress("ESPRITDC",&tdc_array);
  
  Int_t neve = tespri->GetEntries();
  TArtTDCHit *hit_tdc = 0;
  
  /************* Fill data to new tree*******************/
  for(Int_t ne=0;ne<neve;ne++){
    tespri->GetEntry(ne);
    Int_t ntdc = tdc_array->GetEntries();
    for(Int_t i=0;i<16;i++){tdcval[i]=-9999;} 
    // Int_t nbdc = bdc_array->GetEntries(); 
    /**********tdc from ESPRITDC          ****************/  
    for(Int_t i=0;i<ntdc;i++){
      TArtTDCHit *hit_tdc = (TArtTDCHit *)tdc_array->At(i);
      Int_t id_plane = hit_tdc->GetPlaneID();
      Int_t layer = hit_tdc->GetLayer();
      Int_t toffset = hit_tdc->GetTzero();
      Int_t fltdc = hit_tdc->GetTDC();
      
      for(Int_t j=0;j<8;j++){
	if(layer==j+1){
	  if(id_plane==j+1){
	    tdcval[j]=fltdc;
	    h[j]->Fill(tdcval[j]);
	  //  std::cout<<id_plane<<":"<<layer<<":"<<toffset<<":"<<tdcval[j]<<std::endl;
	  }
	  else if(id_plane==j+9){
	    tdcval[j+8]=fltdc;
	    h[j+8]->Fill(tdcval[j+8]);
	    //	    std::cout<<id_plane<<":"<<layer<<":"<<fltdc<<":"<<tdcval[j+8]<<std::endl;
	  }
	}
      }
    }//end of the roop ntdc
    
    // for(Int_t i=0;i<16;i++){
    //   h[i]->Fill(tdcval[i]);
    //   std::cout<<tdcval[i]<<std::endl;    
    // }
    signal(SIGINT,stop_interrupt); //CTRL + C, interrupt
    
    if(ne%100==0){
      printf("Event:%10d/%10d\r",ne,neve);
      fflush(stdout);    
    }
    
  }//end of neve(tespri) roop & fill()
  /*************************************************************************/
  std::string output = outdir+"/"+"BDCcalib_"+outfile;
  TFile *newfile = new TFile(Form("anaresult/%s",output.c_str()),"RECREATE");
  for(int i=0;i<16;i++){
    h[i]->Write();
  }
  TCanvas *ct = new TCanvas("cBDCfLTDC","cBDCfLTDC",1500,1200);
  ct->Divide(4,4);
  ct->cd(0);
  Int_t fmin = 1.2E+3; Int_t fmax = 2.5E+3;
  for(int i=0;i<16;i++){
    ct->cd(i+1);
    // h[i]->Fit("gaus","","",fmin,fmax);
    h[i]->Draw();
    // TF1 *g = (TF1*)h[i]->GetListOfFunctions()->FindObject("gaus");
    // parcst[i] = g->GetParameter(0);
    // parmean[i] = g->GetParameter(1);
    // parsigma[i] = g->GetParameter(2);
    nh[i] = i+1;
  }
  tree->Fill();
  tree->Write();
  ct->Write();
  ct->Update();
  newfile->Close();  

  std::cout<<"input = "<<fname<<", evenumber = "<<neve<<std::endl;
  std::cout<<"analysis rootfile is created as anaresult/"<<output<<std::endl;
  return 0;
}
