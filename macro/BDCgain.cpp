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
#include <TEllipse.h>
#include <TString.h>
#include <TLatex.h>
#include "TArtStoreManager.hh"
#include "TClonesArray.h"
#include "TArtEventStore.hh"
#include "TArtRawEventObject.hh"
#include "TArtRawFADCDataObject.hh"
#include "TArtESPRIParameters.hh"
#include "TArtUserParameters.hh"
#include "TArtCalibTDCHit.hh"
#include "TArtCalibPlas.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtPlas.hh"
#include "TArtTDCHit.hh"
#include "TArtCalibSRPPAC.hh"
#include "TArtSRPPAC.hh"
#include "signal.h"

#endif

bool stoploop = false;
void stop_interrupt(int sig){
  printf("keybprd interrupt\n");
  stoploop = true;
}
/*****const value*********/

//Trigger condition of f12N2 
const int tmin_f12n2=-8000; const int tmax_f12n2=-6000;//Sn&Ca
//const int qmin_f12n2=1000;  const int qmax_f12n2=2000;//Sn
const int qmin_f12n2=700;  const int qmax_f12n2=2500;//Ca
const int tmin_f7dia=8300; const int tmax_f7dia=10700;//Ca
const int tmin_f3dia=-1300; const int tmax_f3dia=1000;//Ca


int main(int argc, char * argv[]){
  std::string fname; Long64_t emax; Int_t piornot;
  if(argc == 1){
    fname = "~/quser/work/enyo/rootfile/5056.root";
    //    emax = 5.E+4;
    piornot=0;
  }else if(argc == 2){
    fname = argv[1];
    // emax = -1;
    piornot=0;
  }else if(argc == 3){
    fname = argv[1];
    // emax = atoi(argv[2]);
    piornot= atoi(argv[2]);
    std::cout << fname << std::endl;
  }else{
    std::cout << "Usage: ./BDCgain 0000.root"<<std::endl;
    return -1;
  }
  std::string outdir = fname.substr(9,2);
  std::string outfile = fname.substr(12);
  std::string runnumber = fname.substr(12,4);
 
  std::string output = outdir+"/"+"BDCgain_"+outfile;
  //TFile *anafile = new TFile(Form("anaresult/%s",output.c_str()),"RECREATE");


  ///////check!! if your file is before PI just comment out 
  TFile *anafile = new TFile("anaresult/Ca/BDCgain_3133mergedpi.root","RECREATE");

  TChain *tespri = new TChain("rtree");
  TString strroot = Form("%s",fname.c_str());
  tespri->Add(strroot);

  //***for PI
  Double_t vx_f5[2][2];//BLD52(SRPPAC)[51 or 52][X or Y]//Ca run has only 52X  
  Double_t vxe12[2][2]={{-1111,-1111},{-1111,-1111}};//F12Xe raw [left or right] [tdc or adc]
  Double_t vxe12d[2][2]={{-1111,-1111},{-1111,-1111}};//F12XeDelay raw [left or right] [tdc or adc]
  Double_t vn2[2][2]={{-1111,-999999},{-1111,-999999}};//F12N2 raw [left or right] [tdc or adc]
  Double_t vdia[2][4]={-999999,-999999};//Dia timing raw [f3 or f7][pad 1-4] 
  Double_t vdia_mean[2]={-999999,-999999};//Dia timing raw [f3 or f7] 
  //for analysis
  Double_t vtof_f37=-9999;
  Double_t vxe12_q=-1111;

  Double_t cx1=9945.;Double_t cy1=1015.;
  Double_t rx1=56.;Double_t ry1=52.; 
  TEllipse *cca48_3133 = new TEllipse(cx1,cy1,rx1,ry1,0.,360.,0);  
  TH2F *hpi=new TH2F("hpi","f3f7tofhsi vs f12XeQ",200,9600,10400,100,400,1500);
  TCanvas *ccut = new TCanvas("cut48","cut48");

  //*******Define global variables in variable.h***********
  Double_t vn2_rt[2][50];//F12N2 tdc raw [left&right][multiplicity]
  Double_t vbdc_L[16][16][10];//Leading [BDC1X*4,Y*4,BDC2X,Y][wireid][TDCnumber]
  Double_t vbdc_T[16][16][10];//Trailing  [BDC1X*4,Y*4,BDC2X,Y][wireid][TDCnumber]
  Double_t vbdc_tot[16][16][10];//TOT[BDC1X*4,Y*4,BDC2X,Y][wire][TDCnum]
  Double_t vn2_t[2];//In beam bunch
  Double_t vbdc_oL[16][16];//Leading(orderd by TOT)[BDC1X*4,Y*4,BDC2X,Y][number]
  Double_t vbdc_oT[16][16];//Trailing(orderd by TOT)[BDC1X*4,Y*4,BDC2X,Y][number]
  Double_t vbdc_otot[16][16];//TOT(ordered)[BDC1X*4,Y*4,BDC2X,Y][TDCnum]
  
  Int_t nmw[16][16]={0};//TDC number of each wire [BDC1X*4,Y*4,BDC2X,Y][wireid]
  /******************************************************/
  //result
  Int_t nm[16]={-1};//multiplicity[BDC1X*4,Y*4,BDC2X,Y]
  Int_t neff[16][5];//Integral of each multiplicity[BDC1X*4,Y*4,BDC2X,Y][all,1,2,3,4]
  Int_t nespri=0;//Integral of f12N2(LR) is in gate
  Int_t nespripi=0;//Integral of f12N2(LR) is in gate with PI
  char name[100];

  //*********new TTree for analysis*************/
  TTree *tana = new TTree("tana","tana");
  
  tana->Branch("f12n2_rt",&vn2_rt,"f12n2_rt[2][50]/D");//raw TDC[Left or Right]
  tana->Branch("f12n2_t",&vn2_t,"f12n2_t[2]/D");//TDC[Left or Right]

  tana->Branch("bdc_L",&vbdc_L,"bdc_L[16][16][10]/D");
  tana->Branch("bdc_T",&vbdc_T,"bdc_T[16][16][10]/D");
  tana->Branch("bdc_tot",&vbdc_tot,"bdc_tot[16][16][10]/D");
  tana->Branch("bdc_oL",&vbdc_oL,"bdc_oL[16][16]/D");
  tana->Branch("bdc_oT",&vbdc_oT,"bdc_oT[16][16]/D");
  tana->Branch("bdc_otot",&vbdc_otot,"bdc_otot[16][16]/D");
  tana->Branch("bdc_m",&nm,"nm[16]/I");
  /**************************************************************/
  TTree *ttrans = new TTree("ttrans","ttrans");

  ttrans->Branch("nespri",&nespri,"nespri/I");//hit number to f12N2

  /***********************************/  

  /*******histgram********************************/
  TH1F *hm[16]; Int_t nbin_m=12; Int_t nmin_m=-1; Int_t nmax_m=11;
  hm[0]=new TH1F("hm_0","multi BDC1X1",nbin_m,nmin_m,nmax_m);
  hm[1]=new TH1F("hm_1","multi BDC1X2",nbin_m,nmin_m,nmax_m);
  hm[2]=new TH1F("hm_2","multi BDC1X3",nbin_m,nmin_m,nmax_m);
  hm[3]=new TH1F("hm_3","multi BDC1X4",nbin_m,nmin_m,nmax_m);
  hm[4]=new TH1F("hm_4","multi BDC1Y1",nbin_m,nmin_m,nmax_m);
  hm[5]=new TH1F("hm_5","multi BDC1Y2",nbin_m,nmin_m,nmax_m);
  hm[6]=new TH1F("hm_6","multi BDC1Y3",nbin_m,nmin_m,nmax_m);
  hm[7]=new TH1F("hm_7","multi BDC1Y4",nbin_m,nmin_m,nmax_m);
  hm[8]=new TH1F("hm_8","multi BDC2X1",nbin_m,nmin_m,nmax_m);
  hm[9]=new TH1F("hm_9","multi BDC2X2",nbin_m,nmin_m,nmax_m);
  hm[10]=new TH1F("hm_10","multi BDC2X3",nbin_m,nmin_m,nmax_m);
  hm[11]=new TH1F("hm_11","multi BDC2X4",nbin_m,nmin_m,nmax_m);
  hm[12]=new TH1F("hm_12","multi BDC2Y1",nbin_m,nmin_m,nmax_m);
  hm[13]=new TH1F("hm_13","multi BDC2Y2",nbin_m,nmin_m,nmax_m);
  hm[14]=new TH1F("hm_14","multi BDC2Y3",nbin_m,nmin_m,nmax_m);
  hm[15]=new TH1F("hm_15","multi BDC2Y4",nbin_m,nmin_m,nmax_m);
 
  TCanvas *cm1 = new TCanvas("cm1","BDC1 multiplicity",1000,500);
  TCanvas *cm2 = new TCanvas("cm2","BDC2 multiplicity",1000,500);


  /**************************************************/

  TClonesArray *srppac_array=0;
  tespri->SetBranchAddress("ESPRISRPPAC",&srppac_array);
  TClonesArray *tdc_array=0; 
  tespri->SetBranchAddress("ESPRITDC",&tdc_array);
  TClonesArray *plas_array=0;
  tespri->SetBranchAddress("ESPRIPLAS",&plas_array);
  

  Int_t neve = tespri->GetEntries();
  TArtPlas *hit_pla = 0;
  TArtTDCHit *tdc = 0;
  TArtTDCHit *hit_tdc = 0;//for PI & Trigger
  TArtSRPPAC *hit_srppac = 0;

  Int_t tmp1,tmp2,tmp3;  

  for(Int_t i=0;i<16;i++){
    for(Int_t j=0;j<5;j++){
      neff[i][j]=0;
    }
  }
  /************* Fill data to new tree*******************/
 Int_t nentry = tespri->GetEntries();


 for(Int_t ne=0;ne<neve;ne++){
   Int_t trigN2[2]={0};//N2 scinti tdc & qdc ok [left&right]
   Int_t flugpi=0;
   Int_t FLUG_PI=0;//final flug of PI
   tespri->GetEntry(ne);
   Int_t npla = plas_array->GetEntries(); 
   Int_t ntdc = tdc_array->GetEntries();
   Int_t nsrppac = srppac_array->GetEntries();   
   Int_t nmn2[2]={0};//N2 multiplicity
   
   /********ESPRITDC***************************/
   for(Int_t j=0;j<2;j++){
     vn2_t[j]=-999999;
     for(Int_t jj=0;jj<16;jj++){
       vn2_rt[j][jj]=-999999;
     }    
   }
   vdia_mean[0]=-999999;vdia_mean[1]=-999999;
   for(Int_t j=0;j<4;j++){
     vdia[0][j]=-999999;vdia[1][j]=-999999;
   }
 
   for(Int_t i=0;i<ntdc;i++){
     TArtTDCHit *tdc = (TArtTDCHit *)tdc_array->At(i);
     Int_t layer = tdc->GetLayer();
     Int_t ch = tdc->GetWireID();
     Int_t id_plane = tdc->GetPlaneID();
     Int_t wireid = tdc->GetWireID();
     Int_t fltdc = tdc->GetTDC();
     
     if(layer==1&&id_plane==31){
       for(Int_t j=0;j<2;j++){
	 if(wireid=j+1){
	   vn2_rt[j][nmn2[j]]=fltdc;
	   //std::cout<<nmd[j]<<":"<<fltdc<<std::endl;
	   nmn2[j]++; 
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
     
     for(Int_t i=0;i<2;i++){
       if(vdia[i][0]>-999990&&vdia[i][1]>-999990
	  &&vdia[i][2]>-999990&&vdia[i][3]>-999990){
	 vdia_mean[i]=(vdia[i][0]+vdia[i][1]+vdia[i][2]+vdia[i][3])/4.;
       }
     }
   }
   //*******TRIGGER by TDC **************************************/
   for(Int_t i=0;i<2;i++){
     for(Int_t j=0;j<nmn2[i];j++){
       if(vn2_rt[i][j]>tmin_f12n2&&vn2_rt[i][j]<tmax_f12n2){//in beam bunch
	 vn2_t[i]=vn2_rt[i][j];
	 break;
       }
     } 
   }
   /** Trigger & setting **/
   trigN2[0]=0;trigN2[1]=0;
   for(Int_t i=0;i<2;i++){//N2 trig 
     if(vn2_t[i]>tmin_f12n2&&vn2_t[i]<tmax_f12n2){
       // if(vn2_rq[i]>qmin_f12n2 && vn2_rq[i]<qmax_f12n2){
       trigN2[i]=1;
       // 	//	}
     }
   }
   if(trigN2[0]*trigN2[1]==1) nespri++;
   
   /*********  ESPRISRPPAC    *******************************/
   vx_f5[0][0]=-999;vx_f5[0][1]=-999;vx_f5[1][0]=-999;vx_f5[1][1]=-999;
   for(Int_t i=0;i<nsrppac;i++){
     TArtSRPPAC *hit_srppac = (TArtSRPPAC *)srppac_array->At(i);
     Int_t layer = hit_srppac->GetLayer();
     vx_f5[layer-1][0] = hit_srppac->GetSRPPACX();
     vx_f5[layer-1][1] = hit_srppac->GetSRPPACY();
   }
   /********ESPRIPLAS*********/
   for(Int_t i=0;i<npla;i++){
     TArtPlas *hit_pla = (TArtPlas *)plas_array->At(i);
     Int_t id_plane = hit_pla->GetPlaneID();
     Int_t layer = hit_pla->GetLayer();
     Int_t channel = hit_pla->GetCh();
     
     if(layer==3 && id_plane==33){
       if(channel==1){//F12Xe left no Delay
	 vxe12[0][1] = (Double_t)hit_pla->GetEnergy();
       }else if(channel==2){//F12Xe right no Delay
	 vxe12d[0][1] = (Double_t)hit_pla->GetEnergy();
       }else if(channel==3){//F12Xe right no Delay
	 vxe12[1][1] = (Double_t)hit_pla->GetEnergy();
       }else if(channel==4){//F12Xe right no Delay
	 vxe12d[1][1] = (Double_t)hit_pla->GetEnergy();
       }
     }
   } 
   
   Double_t tmpx; Double_t tmpy;
   vtof_f37 = (vdia_mean[1]-vdia_mean[0]); //ToF between f3 & f7 Dia
   vxe12_q = (vxe12[0][1]-vxe12d[0][1]+vxe12[1][1]-vxe12d[1][1])/2.; //F12 Xe q
   tmpx=vtof_f37-(vx_f5[1][0]-60)*4.5;
   tmpy=vxe12_q;
   if(trigN2[0]*trigN2[1]==1){
     if((tmpx-cx1)*(tmpx-cx1)/rx1/rx1+(tmpy-cy1)*(tmpy-cy1)/ry1/ry1-1.<=0){
       hpi->Fill(vtof_f37-(vx_f5[1][0]-60)*4.5,vxe12_q);     
       flugpi=1;
       nespripi++;
     }
   }
   
   if(piornot==1){
     nespri=nespripi;
     if(flugpi==1){FLUG_PI=1;}
     else{FLUG_PI=0;}
   }else{//when dont want to choose Ca
     nespri=nespri; 
     FLUG_PI=1;
   }

   if(FLUG_PI==1){
     //**************BDC gain**********************/
     for(Int_t j=0;j<16;j++){
       nm[j]=-1;
       for(Int_t jj=0;jj<10;jj++){
	 for(Int_t jjj=0;jjj<16;jjj++){
	   vbdc_L[j][jjj][jj]=-999999;
	   vbdc_T[j][jjj][jj]=-999999;
	   vbdc_tot[j][jjj][jj]=-999999;
	 }
       }   
       for(Int_t jj=0;jj<16;jj++){
       nmw[j][jj]=0;
       }
       for(Int_t jj=0;jj<16;jj++){
	 vbdc_oL[j][jj]=-999999;
	 vbdc_oT[j][jj]=-999999;
	 vbdc_otot[j][jj]=-999999;
       }
     }
     
     /*************TDC from ESPRITDC***********/
     
     for(Int_t i=0;i<ntdc;i++){
       TArtTDCHit *hit_tdc = (TArtTDCHit *)tdc_array->At(i);
       Int_t id_plane = hit_tdc->GetPlaneID();
       Int_t layer = hit_tdc->GetLayer();
       Int_t wireid = hit_tdc->GetWireID();
       Int_t toffset = hit_tdc->GetTzero();
       Int_t fltdc = hit_tdc->GetTDC();
       Int_t fttdc = hit_tdc->GetTrailTDC();
       
       for(Int_t ud=0;ud<2;ud++){//BDC1or2
	 for(Int_t xy=0;xy<2;xy++){//XorY
	   for(Int_t k=0;k<4;k++){//layer
	     if(layer==(xy*4+k)+1&&id_plane==(ud*8+xy*4+k)+1){
	       if(wireid>=1&&wireid<=16){
		 if(fltdc>-999990){
		   vbdc_L[ud*8+xy*4+k][wireid-1][nmw[ud*8+xy*4+k][wireid-1]]=fltdc;
		   vbdc_T[ud*8+xy*4+k][wireid-1][nmw[ud*8+xy*4+k][wireid-1]]=fttdc;
		   vbdc_tot[ud*8+xy*4+k][wireid-1][nmw[ud*8+xy*4+k][wireid-1]]=fttdc-fltdc;
		   nmw[ud*8+xy*4+k][wireid-1]++;
		 }
	       }
	     }
	   }
	 }
       }
     }//end of the roop ntdc
     //std::cout<<nmw[2][50]<<std::endl;
     
     /**************  analysis here   ******************************/
     
     /** initialization ***/
     for(Int_t i=0;i<16;i++){//BDC1X*4Y*4,BDC2XY
       for(Int_t ii=0;ii<16;ii++){//wire
	 for(Int_t jm=0;jm<nmw[i][ii];jm++){
	   Double_t tofn2bdc=vbdc_L[i][ii][jm]/10.-(vn2_t[0]+vn2_t[1])/80.;
	   if(tofn2bdc>160.&&tofn2bdc<300.&&trigN2[0]*trigN2[1]==1){//bunch + F12 Trig
	     //if(tofn2bdc>160.&&tofn2bdc<300.){//bunch + F12 Trig	  
	     //first hit at the bunch 
	     vbdc_oL[i][ii]=vbdc_L[i][ii][jm];
	     vbdc_oT[i][ii]=vbdc_T[i][ii][jm];
	     vbdc_otot[i][ii]=vbdc_tot[i][ii][jm];
	     //	    std::cout<<ne<<":"<<i<<":"<<ii<<":"<<jm<<std::endl;	    
	     break;	  
	   }
	 }
       }
     }
     /*** Orderd by TOT **/
     for(Int_t i=0;i<16;i++){//BDC1XY,BDC2XY
       for(Int_t ii=0;ii<16;ii++){//wire
	 for(Int_t j=ii+1;j<16;j++){
	   if(vbdc_otot[i][j]>vbdc_otot[i][ii]&&vbdc_otot[i][j]>0&&vbdc_otot[i][j]<20000){
	     tmp1=vbdc_oL[i][ii];  tmp2=vbdc_oT[i][ii];  tmp3=vbdc_otot[i][ii]; 
	     vbdc_oL[i][ii]=vbdc_oL[i][j]; vbdc_oL[i][j]=tmp1;
	     vbdc_oT[i][ii]=vbdc_oT[i][j]; vbdc_oT[i][j]=tmp2;
	     vbdc_otot[i][ii]=vbdc_otot[i][j]; vbdc_otot[i][j]=tmp3;
	   }
	 }
       }
     }
     //multiplivcity of each event
     for(Int_t i=0;i<16;i++){//BDC1XY,BDC2XY
       for(Int_t ii=0;ii<16;ii++){//number orderd by TOT 
	 if(vbdc_otot[i][ii]>0&&vbdc_otot[i][ii]<2000&&trigN2[0]*trigN2[1]==1){
	   if(nm[i]==-1){
	     nm[i]=nm[i]+2;
	   }
	   else{
	     nm[i]++;
	   }
	 }else{
	   nm[i]=nm[i];
	 }
       }
       if(trigN2[0]*trigN2[1]==1){
	 hm[i]->Fill(nm[i]);
       }
     }
     
     for(Int_t i=0;i<16;i++){
       for(Int_t im=1;im<5;im++){  
	 if(nm[i]==im){
	   neff[i][im]++;
	 }
       }
       for(Int_t im=1;im<=16;im++){  
	 if(nm[i]==im){
	   neff[i][0]++;
	 }
       }
     }
     //std::cout<<neff[0][0]<<std::endl;
     /*****************************************/
     //FILL to tree!!!
     tana->Fill();//event Fill
     
     
     if(ne%100==0){
       printf("Event:%10d/%10d\r",ne,neve);
       fflush(stdout);    
     }
     if(ne%10000==0){
       tana->Write();
     }
   }//end of FLUG_PI
   
   
 }//end of ne(tespri) roop to neve & fill()
 std::cout<<nespri<<":"<<nespri<<":"<<neff[8][0]<<":"<<neff[8][1]<<":"<<neff[8][2]<<std::endl;
 
 std::cout<<neff[0][0]<<":"<<neff[1][0]<<":"<<neff[2][0]<<":"<<neff[3][0]<<":"<<
   neff[4][0]<<":"<<neff[5][0]<<":"<<neff[6][0]<<":"<<neff[7][0]<<":"<<
   neff[8][0]<<":"<<neff[9][0]<<":"<<neff[10][0]<<":"<<neff[11][0]<<":"<<
   neff[12][0]<<":"<<neff[13][0]<<":"<<neff[14][0]<<":"<<neff[15][0]<<std::endl;  
 /************* make canvas ****************************/
 cm1->Divide(4,2);
 TLatex *text[50];
 for(Int_t i=0;i<8;i++){
   sprintf(name, "eff all = %0.4f\n",(Double_t)neff[i][0]/nespri);
   text[i*8]=new TLatex(4,hm[i]->GetMaximum()*0.5,name);
   sprintf(name, "eff m1 = %0.4f\n",(Double_t)neff[i][1]/nespri);
   text[i*8+1]=new TLatex(4,hm[i]->GetMaximum()*0.4,name);
   sprintf(name, "eff m2 = %0.4f\n",(Double_t)neff[i][2]/nespri);
   text[i*8+2]=new TLatex(4,hm[i]->GetMaximum()*0.3,name);
   sprintf(name, "eff m3 = %0.4f\n",(Double_t)neff[i][3]/nespri);
   text[i*8+3]=new TLatex(4,hm[i]->GetMaximum()*0.2,name);
   sprintf(name, "eff m4 = %0.4f\n",(Double_t)neff[i][4]/nespri);
   text[i*8+4]=new TLatex(4,hm[i]->GetMaximum()*0.1,name);
   for(Int_t j=0;j<5;j++){text[i*8+j]->SetTextSize(0.05);  text[i*8+j]->SetTextColor(2);}
 }
 TText *t = new TText(4,hm[0]->GetMaximum()*0.9,runnumber.c_str());t->SetTextSize(0.1);
 for(int i=0;i<8;i++){
   cm1->cd(i+1); hm[i]->Draw();
   if(i==0)t->Draw();
   for(int j=0;j<5;j++){   
     text[i*8+j]->Draw("same");
   } 
 }
 cm2->Divide(4,2);
 TLatex *text2[50];
 for(Int_t i=0;i<8;i++){
   sprintf(name, "eff all = %0.4f\n",(Double_t)neff[i+8][0]/nespri);
   text2[i*8]=new TLatex(4,hm[i+8]->GetMaximum()*0.5,name);
   sprintf(name, "eff m1 = %0.4f\n",(Double_t)neff[i+8][1]/nespri);
   text2[i*8+1]=new TLatex(4,hm[i+8]->GetMaximum()*0.4,name);
   sprintf(name, "eff m2 = %0.4f\n",(Double_t)neff[i+8][2]/nespri);
   text2[i*8+2]=new TLatex(4,hm[i+8]->GetMaximum()*0.3,name);
   sprintf(name, "eff m3 = %0.4f\n",(Double_t)neff[i+8][3]/nespri);
   text2[i*8+3]=new TLatex(4,hm[i+8]->GetMaximum()*0.2,name);
   sprintf(name, "eff m4 = %0.4f\n",(Double_t)neff[i+8][4]/nespri);
   text2[i*8+4]=new TLatex(4,hm[i+8]->GetMaximum()*0.1,name);
   for(Int_t j=0;j<5;j++){text2[i*8+j]->SetTextSize(0.05);  text2[i*8+j]->SetTextColor(2);}
 }
 t = new TText(4,hm[8]->GetMaximum()*0.9,runnumber.c_str());
 t->SetTextSize(0.1);
 for(int i=0;i<8;i++){
   cm2->cd(i+1); hm[i+8]->Draw();
   if(i==0)t->Draw();
   for(int j=0;j<5;j++){   
     text2[i*8+j]->Draw("same");
   } 
 }
 
  
 /**********************************************/

 tana->Write();
 ttrans->Fill();
 ttrans->Write();
 for(int i=0;i<16;i++){
   hm[i]->Write();
 }
 cm1->Write();
 cm2->Write();
 hpi->Write();
 anafile->Close(); 
 
 std::cout<<"input = "<<fname<<", evenumber = "<<neve<<std::endl;
 std::cout<<"analysis rootfile is created as anaresult/"<<output<<std::endl;
 return 0;
}


