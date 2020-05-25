#ifndef __CINT__
#include <iostream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TText.h>
#include <TString.h>
#include <TGraph.h>
#include <TGraph.h>
#include <TMarker.h>
#include "TArtStoreManager.hh"
#include "TClonesArray.h"
#include "TArtEventStore.hh"
#include "TArtESPRIParameters.hh"
#include "TArtUserParameters.hh"
#include "TArtCalibPlas.hh"
#include "TArtPlas.hh"
#include "TArtCalibTDCHit.hh"
#include "TArtTDCHit.hh"

#endif
#include "signal.h"
#include "include/Divide2.h"

bool stoploop = false;
void stop_interrupt(int sig){
  printf("keybprd interrupt\n");
  stoploop = true;
}
/*****const value*********/
//const int tmin_f12n2=-8000; const int tmax_f12n2=-6000;//Sn&Ca
const Int_t tmin_f12n2[2]={-8000,-7580}; const Int_t tmax_f12n2[2]={-6000,-5580};//Ca(<3135,3135<=)
const int tmin_f7dia=8300; const int tmax_f7dia=10700;//Ca
const int tmin_f3dia=-1300; const int tmax_f3dia=1000;//Ca

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
    std::cout << "Usage: ./Xe 0000.root"<<std::endl;
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
  Int_t tmin_f12n2_r; Int_t tmax_f12n2_r;
  Double_t vxe12[2][2]={{-1111,-1111},{-1111,-1111}};//F12Xe raw [left or right] [tdc or adc]
  Double_t vxe12d[2][2]={{-1111,-1111},{-1111,-1111}};//F12XeDelay raw [left or right] [tdc or adc]
  Double_t vxe8[2][2]={{-1111,-1111},{-1111,-1111}};//F8Xe raw [left or right] [tdc or adc]
  Double_t vxe8d[2][2]={{-1111,-1111},{-1111,-1111}};//F8XeDelay raw [left or right] [tdc or adc]
  Double_t vn2[2][2]={{-1111,-1111},{-1111,-1111}};//F12N2 raw [left or right] [tdc or adc]
  Double_t vxe12_q=-1111;
  Double_t vxe8_q=-1111;
  Int_t trigN2[2]={-1111,-1111};

  if(Nrun>=3135&&Nrun<=3140){
    tmin_f12n2_r=tmin_f12n2[1]; tmax_f12n2_r=tmax_f12n2[1];
  }else{tmin_f12n2_r=tmin_f12n2[0]; tmax_f12n2_r=tmax_f12n2[0];}
/******************************************************/
  Int_t nhist=2;
  TH1F *h[nhist];
  h[0]=new TH1F("h_00","F12XeQ LR mean ",400,500,1500);
  h[1]=new TH1F("h_01","F12XeQ LR mean Calib ",150,17,23);

  //*********new TTree for analysis*************/
  TTree *tana = new TTree("tana","tana");

  tana->Branch("f12xe",&vxe12,"f12xe[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f12xed",&vxe12d,"f12xed[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f8xe",&vxe8,"f8xe[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f8xed",&vxe8d,"f8xed[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f12n2",&vn2,"f12n2[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("trigN2t",&trigN2,"trigN2t[2]/I");//[Left or Right] N2 Timing

  /**************************************************************/

  /**************************************************/

  TClonesArray *plas_array=0;
  tespri->SetBranchAddress("ESPRIPLAS",&plas_array);
  TClonesArray *tdc_array=0;
  tespri->SetBranchAddress("ESPRITDC",&tdc_array);

  
  Int_t neve = tespri->GetEntries();
  TArtPlas *hit_pla = 0;
  TArtTDCHit *hit_tdc = 0;

  /************* Fill data to new tree*******************/
  for(Int_t ne=0;ne<neve;ne++){
    Int_t trigveto=1;//veto
    tespri->GetEntry(ne);
    Int_t npla = plas_array->GetEntries(); 
    Int_t ntdc = tdc_array->GetEntries();
   
    /**********qdc from ESPRIPLAS          ****************/  
    for(Int_t i=0;i<2;i++){
      for(Int_t ii=0;ii<2;ii++){
	vxe12[i][ii]=-999999;vxe8[i][ii]=-999999;
	vxe12d[i][ii]=-999999;vxe8d[i][ii]=-999999;
      }
    }
    for(Int_t i=0;i<2;i++){
      vn2[i][0]=-999999;
    }

    for(Int_t i=0;i<npla;i++){
      TArtPlas *hit_pla = (TArtPlas *)plas_array->At(i);
      Int_t id_plane = hit_pla->GetPlaneID();
      Int_t layer = hit_pla->GetLayer();
      Int_t channel = hit_pla->GetCh();
      
      if(layer==3 && id_plane==33){
	if(channel==1){//F12Xe left no Delay
	  vxe12[0][1] = (Double_t)hit_pla->GetEnergy();
	}else if(channel==2){//F12Xe right with Delay
	  vxe12d[0][1] = (Double_t)hit_pla->GetEnergy();
	}else if(channel==3){//F12Xe right no Delay
       	  vxe12[1][1] = (Double_t)hit_pla->GetEnergy();
	}else if(channel==4){//F12Xe right with Delay
	  vxe12d[1][1] = (Double_t)hit_pla->GetEnergy();
	}
      }
      if(layer==5 && id_plane==35){
      	if(channel==1){//F8Xe left no Delay
	  vxe8[0][1] = (Double_t)hit_pla->GetEnergy();
	}else if(channel==2){//F8Xe left with Delay
	  vxe8d[0][1] = (Double_t)hit_pla->GetEnergy();
      	}else if(channel==3){//F8Xe left no Delay
      	  vxe8[1][1] = (Double_t)hit_pla->GetEnergy();
      	}else if(channel==4){//F8Xe left with Delay
      	  vxe8d[1][1] = (Double_t)hit_pla->GetEnergy();
      	}
      }
    }//end of the roop npla

    /**********fLTDC from ESPRITDC          ****************/  
    

    for(Int_t i=0;i<ntdc;i++){
      TArtTDCHit *hit_tdc = (TArtTDCHit *)tdc_array->At(i);
      Int_t id_plane = hit_tdc->GetPlaneID();
      Int_t layer = hit_tdc->GetLayer();
      Int_t wireid = hit_tdc->GetWireID();
      Int_t toffset = hit_tdc->GetTzero();
      Int_t fltdc = hit_tdc->GetTDC();
      if(layer==3&&id_plane==33){
	if(wireid==1) vxe12[0][0] = fltdc;
	if(wireid==3) vxe12[1][0] = fltdc;
      }
      if(layer==5&&id_plane==35){
	if(wireid==1) vxe8[0][0] = fltdc;
	if(wireid==3) vxe8[1][0] = fltdc;
      }
      if(layer==1&&id_plane==31){
	for(Int_t j=0;j<2;j++){
	  if(wireid=j+1&&fltdc>tmin_f12n2_r&&fltdc<tmax_f12n2_r){
	    vn2[j][0]=fltdc;
	  }
	}
      }

    }//end of the roop ntdc

    /**************  analysis here   ******************************/
    /*-----------------------*/
    for(Int_t i=0;i<2;i++){//N2 trig 
      if(vn2[i][0]>tmin_f12n2_r && vn2[i][0]<tmax_f12n2_r){
	trigN2[i]=1;
      }else{trigN2[i]=0;}
    }
    
    vxe12_q = (vxe12[0][1]-vxe12d[0][1]+vxe12[1][1]-vxe12d[1][1])/2.; //F12 Xe q
    vxe8_q = (vxe8[0][1]-vxe8d[0][1]+vxe8[1][1]-vxe8d[1][1])/2.; //F8 Xe q
       
    /*****************************************/
    //FILL to tree!!!
    tana->Fill();//event Fill

   //Fill to hist when trigger is OK
    h[0]->Fill(vxe12_q);

    if(ne%100==0){
      printf("Event:%10d/%10d\r"
	     ,ne,neve);
      fflush(stdout);    
    }
    if(ne%100==50000){
      tana->Write();   
    }
    
  }//end of ne(tespri) roop to neve & fill()
  
  TCanvas *c12 = new TCanvas("cf12xe","F12Xe",1400,1000);
  TText *t;
  char fpara[200];
  //1:fit by 3or6 separately, 2:fit by 3 gaussian
  Int_t nffnc=3;
  Double_t range[2*nffnc];
  Int_t flugfit=2;
  Double_t xh[nffnc];
  Double_t yh=h[0]->GetMaximum();  
  TF1 *ffit[nffnc+1];
  for(Int_t i=0;i<nffnc;i++){
    ffit[i]=new TF1(Form("fit%i",i),"gaus",0,2000);
  }
  Double_t height[nffnc];
  Double_t mean[nffnc];
  Double_t sigma[nffnc];

  switch(flugfit){
  case 1:
    {
      range[0]=870.;range[1]=960.;range[2]=970.;range[3]=1050.;
      range[4]=1060.;range[5]=1150.;
      //={870.,960.,970.,1050.,1060.,1150};
      //  Double_t range[12]={700.,770.,780.,845.,850.,920.,
      //		      930.,1010.,1020.,1100.,1120.,1200.};
      for(Int_t i=0;i<nffnc;i++) xh[i]=600.;
      //={600,600,600};
      //  Double_t xh[6]={600,600,600,600,600,600}; 

      TF1 *ffit[nffnc+1];
      for(Int_t i=0;i<nffnc;i++){
	ffit[i]->SetParameters(yh*0.9,(range[2*i]+range[2*i+1])/2.,30.);
      }
      
      
      //  c12->Divide(3,2);
      c12->Divide(4,2);
      for(Int_t i=0;i<nffnc;i++){
	c12->cd(i+1);
	h[0]->Draw();
	if(i==0){
	  h[0]->Fit("fit0","0","",range[2*i],range[2*i+1]);
	}else if(i==1){
	  h[0]->Fit("fit1","0","",range[2*i],range[2*i+1]);
	}else if(i==2){
	  h[0]->Fit("fit2","0","",range[2*i],range[2*i+1]);
	}else if(i==3){
	  h[0]->Fit("fit3","0","",range[2*i],range[2*i+1]);
	}else if(i==4){
	  h[0]->Fit("fit4","0","",range[2*i],range[2*i+1]);
	}else if(i==5){
	  h[0]->Fit("fit5","0","",range[2*i],range[2*i+1]);
	}
	ffit[i]->SetLineColor(33+i*10);
	ffit[i]->Draw("SAME");
	height[i]=ffit[i]->GetParameter(0); mean[i]=ffit[i]->GetParameter(1);
	sigma[i]=ffit[i]->GetParameter(2);
	sprintf(fpara,"h:%6.2f+-%5.2f",height[i],ffit[i]->GetParError(0));
	t = new TText(xh[i],yh*0.9,fpara);  t->SetTextSize(0.04); t->Draw();
	sprintf(fpara,"m:%6.3f+-%5.3f",mean[i],ffit[i]->GetParError(1));
	t = new TText(xh[i],yh*0.83,fpara);  t->SetTextSize(0.04); t->Draw();
	sprintf(fpara,"sig:%6.3f+-%5.3f",sigma[i],ffit[i]->GetParError(2));
	t = new TText(xh[i],yh*0.76,fpara);  t->SetTextSize(0.04); t->Draw();
      }
      
      TMarker *mcalib = new TMarker();
      mcalib->SetMarkerStyle(20);
      mcalib->SetMarkerColor(kOrange);
      mcalib->SetMarkerSize(2);
      c12->cd(nffnc+1);
      gPad->SetGrid();
      Double_t Z[nffnc];
      Z[0]=19;Z[1]=20;Z[2]=21;//={19,20,21};
      //  Double_t Z[nffnc]={17,18,19,20,21,22};
      TGraph *gcalib = new TGraph((sizeof(mean)/sizeof(Double_t)),mean,Z);
      gcalib->SetTitle("QDC to Z");
      gcalib->SetFillColor(kYellow); 
      gcalib->Fit("pol2","","",500,1400);
      TF1 *fcalib = gcalib->GetFunction("pol2");
      gcalib->GetXaxis()->SetLimits(500,1400);
      gcalib->GetHistogram()->SetMaximum(14.);
      gcalib->GetHistogram()->SetMaximum(28.);
      gcalib->Draw("AP");
      for(Int_t i=0;i<nffnc;i++){
	mcalib->DrawMarker(mean[i],Z[i]);
      } 
      Double_t pcalib[3];
      pcalib[0] = fcalib->GetParameter(0);
      pcalib[1] = fcalib->GetParameter(1);
      pcalib[2] = fcalib->GetParameter(2);
      
      for(Int_t ne=0;ne<neve;ne++){
	tespri->GetEntry(ne);
	Int_t npla = plas_array->GetEntries(); 
	
	/**********qdc from ESPRIPLAS          ****************/  
	for(Int_t i=0;i<2;i++){
	  for(Int_t ii=0;ii<2;ii++){
	    vxe12[i][ii]=-999999;vxe8[i][ii]=-999999;
	    vxe12d[i][ii]=-999999;vxe8d[i][ii]=-999999;
	  }
	}
	for(Int_t i=0;i<npla;i++){
	  TArtPlas *hit_pla = (TArtPlas *)plas_array->At(i);
	  Int_t channel = hit_pla->GetCh();
	  if(hit_pla->GetLayer()==3 && hit_pla->GetPlaneID()==33){
	    if(channel==1){//F12Xe left no Delay
	      vxe12[0][1] = (Double_t)hit_pla->GetEnergy();
	    }else if(channel==2){//F12Xe right with Delay
	    vxe12d[0][1] = (Double_t)hit_pla->GetEnergy();
	    }else if(channel==3){//F12Xe right no Delay
	      vxe12[1][1] = (Double_t)hit_pla->GetEnergy();
	    }else if(channel==4){//F12Xe right with Delay
	      vxe12d[1][1] = (Double_t)hit_pla->GetEnergy();
	    }
	  }
	}//end of the roop npla
      
	vxe12_q = (vxe12[0][1]-vxe12d[0][1]+vxe12[1][1]-vxe12d[1][1])/2.; //F12 Xe q
	
	//Fill to hist when trigger is OK
	h[1]->Fill(pcalib[0]+pcalib[1]*vxe12_q+pcalib[2]*vxe12_q*vxe12_q);
	
	if(ne%100==0){
	  printf("Event:%10d/%10d\r"
		 ,ne,neve);
	  fflush(stdout);    
	}
      }//end of ne(tespri) roop to neve & fill()
      
      ffit[nffnc]=new TF1("fit6","gaus",10,30);
      ffit[nffnc]->SetParameters(yh*0.98,20.,0.3);
      c12->cd(8);
      h[1]->Draw();
      Double_t yh2=h[1]->GetMaximum();  
      h[1]->Fit("fit6","0","",19.52,20.48);
      ffit[nffnc]->SetLineColor(93);
      ffit[nffnc]->Draw("SAME");
      height[nffnc]=ffit[nffnc]->GetParameter(0); 
      mean[nffnc]=ffit[nffnc]->GetParameter(1);
      sigma[nffnc]=ffit[nffnc]->GetParameter(2);
      sprintf(fpara,"h:%6.2f",height[nffnc]);
      t = new TText(18,yh2*0.9,fpara);  t->SetTextSize(0.04); t->Draw();
      sprintf(fpara,"m:%6.3f",mean[nffnc]);
      t = new TText(18,yh2*0.83,fpara);  t->SetTextSize(0.04); t->Draw();
      sprintf(fpara,"sig:%6.3f",sigma[nffnc]);
      t = new TText(18,yh2*0.76,fpara);  t->SetTextSize(0.04); 
      t->SetTextColor(97); t->Draw();
      
      c12->cd(nffnc+1);
      t = new TText(900,10,runnumber.c_str()); t->SetTextSize(0.1); t->Draw();
    } break;

  /**********************************************/
  case 2://3-gaussian fit
    {
      range[0]=870.;range[1]=960.;range[2]=970.;range[3]=1050.;
      range[4]=1060.;range[5]=1150.;
      //   range[0]=850.;range[1]=920.;range[2]=930.;range[3]=1010.;
      //   range[4]=1020.;range[5]=1100.;
      xh[0]=600.; xh[1]=850.; xh[2]=1100.;
      TF1 *f3gaus;
      f3gaus=new TF1("f3gaus","fit0+fit1+fit2",0,2000);
      f3gaus->SetParameters(yh*0.5,(range[0]+range[1])/2.,30.
			    ,yh*0.9,(range[2]+range[3])/2.,30.
			    ,yh*0.5,(range[4]+range[5])/2.,30.);
      //  Divide2(c12,0.6);
      c12->Divide(1,2);
      c12->cd(1);
      h[0]->SetStats(false); h[0]->Draw();
      h[0]->Fit("f3gaus","0","",range[0],range[5]);
      f3gaus->SetLineColor(53); f3gaus->Draw("SAME");
      for(Int_t i=0;i<nffnc;i++){    
	height[i]=f3gaus->GetParameter(3*i); 
	mean[i]=f3gaus->GetParameter(3*i+1);
	sigma[i]=f3gaus->GetParameter(3*i+2);
	sprintf(fpara,"h:%6.2f+-%5.2f",height[i],f3gaus->GetParError(3*i));
	t = new TText(xh[i],yh*0.9,fpara);  t->SetTextSize(0.04); t->Draw();
	sprintf(fpara,"m:%6.3f+-%5.3f",mean[i],f3gaus->GetParError(3*i+1));
	t = new TText(xh[i],yh*0.83,fpara);  t->SetTextSize(0.04); t->Draw();
	sprintf(fpara,"sig:%6.3f+-%5.3f",sigma[i],f3gaus->GetParError(3*i+2));
	t = new TText(xh[i],yh*0.76,fpara);  t->SetTextSize(0.04); t->Draw();
      }
      std::cout<<mean[0]<<":"<<mean[1]<<":"<<mean[2]<<std::endl;
      c12->cd(2);

      sprintf(fpara,"19-20: %6.2fsig",(mean[1]-mean[0])/sigma[1]);
      t = new TText(0.15,0.8,fpara); t->SetTextSize(0.08); t->Draw();
      sprintf(fpara,"20-21: %6.2fsig",(mean[2]-mean[1])/sigma[1]);
      t = new TText(0.15,0.7,fpara); t->SetTextSize(0.08); t->Draw();
      c12->cd(2);
      t = new TText(0.1,0.9,runnumber.c_str()); t->SetTextSize(0.1); t->Draw();
    } break;
  default:
    std::cout<<"there is no option for fitting."<<std::endl;
    //    break;
  }
  std::stringstream ss;
  ss << flugfit;
  std::string output = outdir+"/"+"Xe"+ss.str()+"_"+outfile;
  TFile *anafile = new TFile(Form("anaresult/%s",output.c_str()),"RECREATE");

  tana->Write();
  for(Int_t i=0;i<nhist;i++){
    h[i]->Write();
  } 
  c12->Write();
  anafile->Close(); 

  std::cout<<"input = "<<fname<<", evenumber = "<<neve<<std::endl;

  std::cout<<"analysis rootfile is created as anaresult/"<<output<<std::endl;
  return 0;
}


