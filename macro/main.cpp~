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
#include <TString.h>
#include "TArtStoreManager.hh"
#include "TClonesArray.h"
#include "TArtEventStore.hh"
#include "TArtRawEventObject.hh"
#include "TArtRawFADCDataObject.hh"
#include "TArtESPRIParameters.hh"
#include "TArtUserParameters.hh"
#include "TArtCalibSRPPAC.hh"
#include "TArtCalibRNaI.hh"
#include "TArtCalibTDCHit.hh"
#include "TArtCalibRDC.hh"
#include "TArtCalibBDC.hh"
#include "TArtCalibPlas.hh"
#include "TArtCalibRNaI.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtSRPPAC.hh"
#include "TArtPlas.hh"
#include "TArtBDC.hh"
#include "TArtTDCHit.hh"
#include "TArtRNaI.hh"
#include "TArtRDC.hh"
#include "signal.h"

#include "./include/variable.h"
#include "./include/func.h"

#endif

bool stoploop = false;
void stop_interrupt(int sig){
  printf("keybprd interrupt\n");
  stoploop = true;
}
/*****const value*********/
const double dis_f712[2] = {50.,50.}; //50m for now
const double c = 2.9979E+8; //light speed [m/s]


//NaI pedestal ?
const double vnai_0[28]={98.7,104.3,76.9,83.8,90.2,101.0
			,99.7,90.1,96.2,86.3,75.6,79.3
			,113.1,109.6,110.7,84.9,103.8
			,115.9,102.5,80.2,94.3,105.9,89.5
			,96.9,113.1,79.6,86.0,104.4};
const double gnai_min = 100;

//PdE pedestal
const double vpde_0[4]={132.6,156.6,152.8,114.3};
const double gpde_min = 100;

//E*dE gate
const double gEtdE[2][2]={{5.E+5,15.E+5},{7.E+5,19.E+5}};//proton gate[up,down][min,max]

//Trigger condition of f12N2 
const int tmin_f12n2=-8000; const int tmax_f12n2=-6000;//Sn&Ca
//const int qmin_f12n2=1000;  const int qmax_f12n2=2000;//Sn
const int qmin_f12n2=700;  const int qmax_f12n2=2500;//Ca

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
    std::cout << "Usage: ./reconst 0000.root"<<std::endl;
    return -1;
  }
  std::string outdir = fname.substr(9,2);
  std::string outfile = fname.substr(12);

  TChain *tespri = new TChain("rtree");
  TString strroot = Form("%s",fname.c_str());
  tespri->Add(strroot);

  //*******Define global variables in variable.h***********
  Double_t vx_f5[2][2];//BLD52(SRPPAC)[51 or 52][X or Y]//Ca run has only 52X  
  Double_t vxe12[2][2]={{-1111,-1111},{-1111,-1111}};//F12Xe raw [left or right] [tdc or adc]
  Double_t vxe8[2][2]={{-1111,-1111},{-1111,-1111}};//F8Xe raw [left or right] [tdc or adc]
  Double_t vn2[2][2]={{-1111,-1111},{-1111,-1111}};//F12N2 raw [left or right] [tdc or adc]
  Double_t vveto[2]={-1111,-1111};//VETO raw [tdc or qdc]
  Double_t vcsi[4][2]={{-1111},{-1111},{-1111},{-1111}};//CsI raw [1,2,3,4] [tdc or adc]
  Double_t vpde[4][2]={{-1111},{-1111},{-1111},{-1111}};//PdE [up L&R, down L&R][time & qdc]
  Double_t vdia[2][4]={-999999,-999999};//Dia timing raw [f3 or f7][pad 1-4] 
  Double_t vdia_mean[2]={-999999,-999999};//Dia timing raw [f3 or f7] 
  Double_t vnai[28][2];//NaI[LR * 14][time & energy]

  double vx[2][4];//before correction raw value 
  double vthx=-111,vthy=-111,vxsht=-111,vysht=-111;//result of bdc tracking digree(x,y), diff(cm) at SHT
  double vbdctrk[4]={-111,-111,-111,-111};//result of bdc tracking Theta digree(x,y), diff(cm) at SHT
  Int_t npp_fin[3]={-111,-111,-111};//final number of N2ok, N2&BDCtrkOK,N2&BDC&protonOK

  //for analysis
  Double_t vtof_f37=-9999;
  Double_t vxe12_t=-1111;
  Double_t vxe12_q=-1111;
  Double_t vn2_t=-1111;
  Double_t vtof_f712[2]={-9999,-9999};//[f7Dia->f12Xe,f7Dia->f12N2]
  Double_t vnai_ped[28];//Nai-ped??[LR * 14]
  Double_t vpde_ped[4];//PdE-ped??[LR * 2]
  Double_t vEtdE[14]={};//E times dE [NaI 14hon]
  Double_t vbeta[2]={-1111,-1111}; //beta of the nucleus[f12Xe,f12N2]
  Double_t vgamma=-1111;//gamma of the nucleus

  Int_t npp[3]={0};//number of proton [n2trigOK, n2trig&BDCtrkOK, n2trig&BDCtrk&proton] 
/******************************************************/

  //*********new TTree for analysis*************/
  TTree *tana = new TTree("tana","tana");

  tana->Branch("bldf5",&vx_f5,"bldf5[2][2]/D");//[BLD51 or 52][X or Y]
  tana->Branch("f12xe",&vxe12,"f12xe[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f8xe",&vxe8,"f8xe[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("f12n2",&vn2,"f12n2[2][2]/D");//[Left or Right][Time or ADC]
  tana->Branch("veto",&vveto,"veto[2]/D");//[Time or ADC]
  tana->Branch("csi",&vcsi,"csi[4][2]/D");//[1,2,3,4][Time or ADC]
  tana->Branch("pde",&vpde,"pde[4][2]/D");//[1,2,3,4][Time or ADC]
  tana->Branch("Dia",&vdia,"Dia[2][4]/D");//[f3Dia or f7Dia][pad 1-4]
  tana->Branch("Dia_t",&vdia_mean,"Dia_t[2]/D");//[f3Dia or f7Dia]
  tana->Branch("nai",&vnai,"nai[28][2]/D");//[LR * 14][time or ADC]
  tana->Branch("bdc_raw",&vx,"bdc_raw[2][4]/D");//[BDC1,BDC2][X,A,Y,B]
  tana->Branch("bdctrk",&vbdctrk,"bdctrk[4]/D");//[thx,thy,xsht,ysht]
  /**************************************************************/
  /********* new TTree for result******************************/
  TTree *tres = new TTree("tres","tres");
  tres->Branch("npp",&npp_fin,"npp[3]/I");//[N2ok, N2&BDCtrkOK,N2&BDCtrk&protonOK]
  /*******************************************************/

  /*******histgram********************************/
  Int_t nhraw_t=17;   Int_t nhraw_q=15; 
  TH1F *hraw_t[nhraw_t];
  hraw_t[0] = new TH1F("hraw_t00","f12xel_t",100,-4E+4,13E+4);
  hraw_t[1] = new TH1F("hraw_t01","f12xer_t",100,-4E+4,13E+4);
  hraw_t[2] = new TH1F("hraw_t02","f8xel_t", 100,-4E+4,13E+4);
  hraw_t[3] = new TH1F("hraw_t03","f8xer_t", 100,-4E+4,13E+4);
  hraw_t[4] = new TH1F("hraw_t04","f12n2l_t",1000,-3E+4,10E+4);
  hraw_t[5] = new TH1F("hraw_t05","f12n2r_t",1000,-3E+4,10E+4);
  hraw_t[6] = new TH1F("hraw_t06","veto_t",100,-4E+4,13E+4);
  hraw_t[7] = new TH1F("hraw_t07","csi01_t",100,-4E+4,13E+4);
  hraw_t[8] = new TH1F("hraw_t08","csi02_t", 100,-4E+4,13E+4);
  hraw_t[9] = new TH1F("hraw_t09","csi03_t", 100,-4E+4,13E+4);
  hraw_t[10] = new TH1F("hraw_t10","csi04_t",100,-4E+4,13E+4);
  hraw_t[11] = new TH1F("hraw_t11","pdeUl_t",100,-4E+4,13E+4);
  hraw_t[12] = new TH1F("hraw_t12","pdeUr_t", 100,-4E+4,13E+4);
  hraw_t[13] = new TH1F("hraw_t13","pdeDl_t", 100,-4E+4,13E+4);
  hraw_t[14] = new TH1F("hraw_t14","pdeDr_t",100,-4E+4,13E+4);
  hraw_t[15] = new TH1F("hraw_t15","f3dia_t", 100,-4E+4,13E+4);
  //  hraw_t[16] = new TH1F("hraw_t16","f7dia_t", 100,-4E+4,13E+4);

  TH1F *hraw_q[nhraw_q];
  hraw_q[0] = new TH1F("hraw_q00","f12xel_q",100,0,3000);
  hraw_q[1] = new TH1F("hraw_q01","f12xer_q",100,0,3000);
  hraw_q[2] = new TH1F("hraw_q02","f8xel_q", 100,0,3000);
  hraw_q[3] = new TH1F("hraw_q03","f8xer_q", 100,0,3000);
  hraw_q[4] = new TH1F("hraw_q04","f12n2l_q",100,0,3000);
  hraw_q[5] = new TH1F("hraw_q05","f12n2r_q",100,0,3000);
  hraw_q[6] = new TH1F("hraw_q06","veto_q",100,0,3000);
  hraw_q[7] = new TH1F("hraw_q07","csi01_q",100,0,3000);
  hraw_q[8] = new TH1F("hraw_q08","csi02_q", 100,0,3000);
  hraw_q[9] = new TH1F("hraw_q09","csi03_q", 100,0,3000);
  hraw_q[10] = new TH1F("hraw_q10","csi04l_q",100,0,3000);
  hraw_q[11] = new TH1F("hraw_q11","pdeUl_q",100,0,3000);
  hraw_q[12] = new TH1F("hraw_q12","pdeUr_q",100,0,3000);
  hraw_q[13] = new TH1F("hraw_q13","pdeDl_q",100,0,3000);
  hraw_q[14] = new TH1F("hraw_q14","pdeDr_q", 100,0,3000);

  TH1F *hnai_t[28]; 
  hnai_t[0] = new TH1F("hnai_t00","nai01l_t", 100,-1000,6000);
  hnai_t[1] = new TH1F("hnai_t01","nai01r_t", 100,-1000,6000);
  hnai_t[2] = new TH1F("hnai_t02","nai02l_t", 100,-1000,6000);
  hnai_t[3] = new TH1F("hnai_t03","nai02r_t", 100,-1000,6000);
  hnai_t[4] = new TH1F("hnai_t04","nai03l_t", 100,-1000,6000);
  hnai_t[5] = new TH1F("hnai_t05","nai03r_t", 100,-1000,6000);
  hnai_t[6] = new TH1F("hnai_t06","nai04l_t", 100,-1000,6000);
  hnai_t[7] = new TH1F("hnai_t07","nai04r_t", 100,-1000,6000);
  hnai_t[8] = new TH1F("hnai_t08","nai05l_t", 100,-1000,6000);
  hnai_t[9] = new TH1F("hnai_t09","nai05r_t", 100,-1000,6000);
  hnai_t[10] = new TH1F("hnai_t10","nai06l_t", 100,-1000,6000);
  hnai_t[11] = new TH1F("hnai_t11","nai06r_t", 100,-1000,6000);
  hnai_t[12] = new TH1F("hnai_t12","nai07l_t", 100,-1000,6000);
  hnai_t[13] = new TH1F("hnai_t13","nai07r_t", 100,-1000,6000);
  hnai_t[14] = new TH1F("hnai_t14","nai08l_t", 100,-1000,6000);
  hnai_t[15] = new TH1F("hnai_t15","nai08r_t", 100,-1000,6000);
  hnai_t[16] = new TH1F("hnai_t16","nai09l_t", 100,-1000,6000);
  hnai_t[17] = new TH1F("hnai_t17","nai09r_t", 100,-1000,6000);
  hnai_t[18] = new TH1F("hnai_t18","nai10l_t", 100,-1000,6000);
  hnai_t[19] = new TH1F("hnai_t19","nai10r_t", 100,-1000,6000);
  hnai_t[20] = new TH1F("hnai_t20","nai11l_t", 100,-1000,6000);
  hnai_t[21] = new TH1F("hnai_t21","nai11r_t", 100,-1000,6000);
  hnai_t[22] = new TH1F("hnai_t22","nai12l_t", 100,-1000,6000);
  hnai_t[23] = new TH1F("hnai_t23","nai12r_t", 100,-1000,6000);
  hnai_t[24] = new TH1F("hnai_t24","nai13l_t", 100,-1000,6000);
  hnai_t[25] = new TH1F("hnai_t25","nai13r_t", 100,-1000,6000);
  hnai_t[26] = new TH1F("hnai_t26","nai14l_t", 100,-1000,6000);
  hnai_t[27] = new TH1F("hnai_t27","nai14r_t", 100,-1000,6000);

  TH1F *hnai_q[28];Int_t nbin_nai=500;Int_t nmin_nai=0;Int_t nmax_nai=3000; 
  hnai_q[0] = new TH1F("hnai_q00","nai01l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[1] = new TH1F("hnai_q01","nai01r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[2] = new TH1F("hnai_q02","nai02l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[3] = new TH1F("hnai_q03","nai02r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[4] = new TH1F("hnai_q04","nai03l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[5] = new TH1F("hnai_q05","nai03r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[6] = new TH1F("hnai_q06","nai04l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[7] = new TH1F("hnai_q07","nai04r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[8] = new TH1F("hnai_q08","nai05l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[9] = new TH1F("hnai_q09","nai05r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[10] = new TH1F("hnai_q10","nai06l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[11] = new TH1F("hnai_q11","nai06r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[12] = new TH1F("hnai_q12","nai07l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[13] = new TH1F("hnai_q13","nai07r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[14] = new TH1F("hnai_q14","nai08l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[15] = new TH1F("hnai_q15","nai08r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[16] = new TH1F("hnai_q16","nai09l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[17] = new TH1F("hnai_q17","nai09r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[18] = new TH1F("hnai_q18","nai10l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[19] = new TH1F("hnai_q19","nai10r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[20] = new TH1F("hnai_q20","nai11l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[21] = new TH1F("hnai_q21","nai11r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[22] = new TH1F("hnai_q22","nai12l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[23] = new TH1F("hnai_q23","nai12r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[24] = new TH1F("hnai_q24","nai13l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[25] = new TH1F("hnai_q25","nai13r_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[26] = new TH1F("hnai_q26","nai14l_q",nbin_nai,nmin_nai,nmax_nai);
  hnai_q[27] = new TH1F("hnai_q27","nai14r_q",nbin_nai,nmin_nai,nmax_nai);

  TH1F *hbdcraw[8];Int_t nbin_bdc=100;Int_t nmin_bdc=-10;Int_t nmax_bdc=20; 
  hbdcraw[0] = new TH1F("hbdcraw_00","hraw BDC1 X",nbin_bdc,nmin_bdc,nmax_bdc);
  hbdcraw[1] = new TH1F("hbdcraw_01","hraw BDC1 A",nbin_bdc,nmin_bdc,nmax_bdc);
  hbdcraw[2] = new TH1F("hbdcraw_02","hraw BDC1 Y",nbin_bdc,nmin_bdc,nmax_bdc);
  hbdcraw[3] = new TH1F("hbdcraw_03","hraw BDC1 B",nbin_bdc,nmin_bdc,nmax_bdc);
  hbdcraw[4] = new TH1F("hbdcraw_04","hraw BDC2 X",nbin_bdc,nmin_bdc,nmax_bdc);
  hbdcraw[5] = new TH1F("hbdcraw_05","hraw BDC2 A",nbin_bdc,nmin_bdc,nmax_bdc);
  hbdcraw[6] = new TH1F("hbdcraw_06","hraw BDC2 Y",nbin_bdc,nmin_bdc,nmax_bdc);
  hbdcraw[7] = new TH1F("hbdcraw_07","hraw BDC2 B",nbin_bdc,nmin_bdc,nmax_bdc);

  TH1F *htrk[4];//hist for tracking (all events when trig is ok)
  htrk[0] = new TH1F("htrk_thx","thetaX",100,-10,10);
  htrk[1] = new TH1F("htrk_thy","thetaY",100,-10,10);
  htrk[2] = new TH1F("htrk_x","X@SHT",100,-10,10);
  htrk[3] = new TH1F("htrk_y","Y@SHT",100,-10,10);
  TH2F *hsht[1];
  hsht[0] = new TH2F("hSHT_ion","ion on SHT",200,-10,10,200,-10,10);
  TH1F *htrk_gate[4];//hist for tracking (events when no veto qdc)
  htrk_gate[0] = new TH1F("htrkg_thx","thetaX with Gate",100,-10,10);
  htrk_gate[1] = new TH1F("htrkg_thy","thetaY with Gate",100,-10,10);
  htrk_gate[2] = new TH1F("htrkg_x","X@SHT with Gate",100,-10,10);
  htrk_gate[3] = new TH1F("htrkg_y","Y@SHT with Gate",100,-10,10);
  TH2F *hsht_gate[1];
  hsht_gate[0] = new TH2F("hSHTg_ion","ion on SHT with Gate",200,-10,10,200,-10,10);

  int npi=4;//number of hist for PI
  TH2F *hpi[npi];  
  hpi[0]=new TH2F("hpi_00","f3f7 tof vs f12XeQ",200,40,60,100,0,2000);
  hpi[1]=new TH2F("hpi_01","f7f12 tof vs f12XeQ",100,-530,-490,100,0,2000);
  hpi[2]=new TH2F("hpi_02","f7f12 tof vs f3f7 tof",100,-530,-490,200,40,60);
  hpi[3]=new TH2F("hpi_03","f3f7ToF*beta vs f12XeQ "
		  ,100,10.,18.,200,700,1500);

  TH2F *hEdE[14];
  Int_t nbin_E=100;Int_t nmin_E=-200;Int_t nmax_E=4000; 
  hEdE[0] = new TH2F("hEdE_00","EdE 00",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[1] = new TH2F("hEdE_01","EdE 01",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[2] = new TH2F("hEdE_02","EdE 02",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[3] = new TH2F("hEdE_03","EdE 03",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[4] = new TH2F("hEdE_04","EdE 04",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[5] = new TH2F("hEdE_05","EdE 05",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[6] = new TH2F("hEdE_06","EdE 06",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[7] = new TH2F("hEdE_07","EdE 07",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[8] = new TH2F("hEdE_08","EdE 08",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[9] = new TH2F("hEdE_09","EdE 09",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[10] = new TH2F("hEdE_10","EdE 10",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[11] = new TH2F("hEdE_11","EdE 11",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[12] = new TH2F("hEdE_12","EdE 12",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdE[13] = new TH2F("hEdE_13","EdE 13",nbin_E,nmin_E,nmax_E,100,0,4000);
  TH2F *hEdE_up=new TH2F("hEdE_up","EdE up",nbin_E,nmin_E,nmax_E,100,-200,4000);
  TH2F *hEdE_down=new TH2F("hEdE_down","EdE down",nbin_E,nmin_E,nmax_E,100,-200,4000);
  TH1F *hEtdE_up=new TH1F("hEtdE_up","E*dE up",2000,0,1E+7);
  TH1F *hEtdE_down=new TH1F("hEtdE_down","E*dE down",2000,0,1E+7);

  //hist for EdE prty
  TH2F *hEdEprt1=new TH2F("hEdEprt1","hEdEprt1",nbin_E,nmin_E,nmax_E,100,0,4000);
  TH2F *hEdEprt2=new TH2F("hEdEprt2","hEdEprt2",nbin_E,nmin_E,nmax_E,100,0,4000);
  TH2F *hEdEprt3=new TH2F("hEdEprt3","hEdEprt3",nbin_E,nmin_E,nmax_E,100,0,4000);
  TH2F *hEdEprt4=new TH2F("hEdEprt4","hEdEprt4",nbin_E,nmin_E,nmax_E,100,0,4000);

  TH2F *hEdEfun[14];
  hEdEfun[0]=new TH2F("hEdEfun_00","hEdEfun00",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[1]=new TH2F("hEdEfun_01","hEdEfun01",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[2]=new TH2F("hEdEfun_02","hEdEfun02",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[3]=new TH2F("hEdEfun_03","hEdEfun03",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[4]=new TH2F("hEdEfun_04","hEdEfun04",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[5]=new TH2F("hEdEfun_05","hEdEfun05",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[6]=new TH2F("hEdEfun_06","hEdEfun06",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[7]=new TH2F("hEdEfun_07","hEdEfun07",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[8]=new TH2F("hEdEfun_08","hEdEfun08",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[9]=new TH2F("hEdEfun_09","hEdEfun09",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[10]=new TH2F("hEdEfun_10","hEdEfun10",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[11]=new TH2F("hEdEfun_11","hEdEfun11",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[12]=new TH2F("hEdEfun_12","hEdEfun12",nbin_E,nmin_E,nmax_E,100,0,4000);
  hEdEfun[13]=new TH2F("hEdEfun_13","hEdEfun13",nbin_E,nmin_E,nmax_E,100,0,4000);


  TCanvas *craw_t = new TCanvas("craw_t","raw TDC");
  TCanvas *craw_q = new TCanvas("craw_q","raw QDC");
  TCanvas *cnai_t = new TCanvas("cnai_t","NaI TDC");
  TCanvas *cnai_q = new TCanvas("cnai_q","NaI QDC");
  TCanvas *cbdcraw = new TCanvas("cbdcraw","BDC raw value",900,500);
  TCanvas *ctrk = new TCanvas("ctrk","focus at SHT");
  TCanvas *cEdE = new TCanvas("cEdE","NaI vs pde");
  

  /**************************************************/

  TClonesArray *srppac_array=0;
  tespri->SetBranchAddress("ESPRISRPPAC",&srppac_array);
  TClonesArray *plas_array=0;
  tespri->SetBranchAddress("ESPRIPLAS",&plas_array);
  TClonesArray *bdc_array=0;
  tespri->SetBranchAddress("ESPRIBDC",&bdc_array);
  TClonesArray *tdc_array=0; 
  tespri->SetBranchAddress("ESPRITDC",&tdc_array);
  TClonesArray *nai_array=0; 
  tespri->SetBranchAddress("ESPRINaI",&nai_array);
  TClonesArray *rdc_array=0; 
  tespri->SetBranchAddress("ESPRIRDC",&rdc_array);

  Int_t neve = tespri->GetEntries();
  int nbdcok=0; //number of the event of "traking ok" 
  TArtSRPPAC *hit_srppac = 0;
  TArtPlas *hit_pla = 0;
  TArtBDC *hit_bdc = 0;
  TArtTDCHit *hit_tdc = 0;
  TArtRNaI *hit_nai = 0;
  TArtRDC *hit_rdc = 0;

  /************* Fill data to new tree*******************/
  for(Int_t ne=0;ne<neve;ne++){
    int flugtrk=0;//BDC both ok
    int trigN2[2]={0};//N2 scinti tdc & qdc ok [left&right]
    int trigveto=1;//veto
    tespri->GetEntry(ne);
    Int_t nsrppac = srppac_array->GetEntries();
    Int_t npla = plas_array->GetEntries(); 
    Int_t nbdc = bdc_array->GetEntries(); 
    Int_t ntdc = tdc_array->GetEntries();
    Int_t nnai = nai_array->GetEntries();
    //    std::cout<<ne<<" "<<nnai<<std::endl;
    /*********  ESPRISRPPAC    *******************************/
    vx_f5[0][0]=-999;vx_f5[0][1]=-999;vx_f5[1][0]=-999;vx_f5[1][1]=-999;
    for(Int_t i=0;i<nsrppac;i++){
      TArtSRPPAC *hit_srppac = (TArtSRPPAC *)srppac_array->At(i);
      Int_t layer = hit_srppac->GetLayer();
     
      vx_f5[layer-1][0] = hit_srppac->GetSRPPACX();
      vx_f5[layer-1][1] = hit_srppac->GetSRPPACY();
    }
    /**********qdc from ESPRIPLAS          ****************/  
    for(Int_t i=0;i<npla;i++){
      TArtPlas *hit_pla = (TArtPlas *)plas_array->At(i);
      Int_t id_plane = hit_pla->GetPlaneID();
      Int_t layer = hit_pla->GetLayer();
      Int_t channel = hit_pla->GetCh();
      
      if(layer==3 && id_plane==33){
	if(channel==1){//F12Xe left no Delay
	  vxe12[0][0] = (Double_t)hit_pla->GetTime();
	  vxe12[0][1] = (Double_t)hit_pla->GetEnergy();
	}else if(channel==3){//F12Xe right no Delay
	  vxe12[1][0] = (Double_t)hit_pla->GetTime();
	  vxe12[1][1] = (Double_t)hit_pla->GetEnergy();
	}
	//else if(channel==2){vdia[0] = (Double_t)hit_pla->GetTime();}//f3Dia timing
	//else if(channel==4){vdia[1] = (Double_t)hit_pla->GetTime();}//f7Dia timing
      }
      if(layer==5 && id_plane==35){
      	if(channel==1){//F8Xe left no Delay
	  vxe8[0][0] = (Double_t)hit_pla->GetTime();//no data->dummy
	  vxe8[0][1] = (Double_t)hit_pla->GetEnergy();
      	}else if(channel==3){//F8Xe left no Delay
      	  vxe8[1][0] = (Double_t)hit_pla->GetTime();//no data->dummy
      	  vxe8[1][1] = (Double_t)hit_pla->GetEnergy();
      	}
      }
      if(layer==1 && id_plane==31){
      	if(channel==1){//N2 left no Delay
      	  vn2[0][0] = (Double_t)hit_pla->GetTime();
      	  vn2[0][1] = (Double_t)hit_pla->GetEnergy();
      	}else if(channel==2){//N2 left no Delay
      	  vn2[1][0] = (Double_t)hit_pla->GetTime();
      	  vn2[1][1]= (Double_t)hit_pla->GetEnergy();
      	}
      }
      if(layer==2 && id_plane==32 && channel==1){
	vveto[0] = (Double_t)hit_pla->GetTime();
	vveto[1] = (Double_t)hit_pla->GetEnergy();
      }
      if(layer==4 && id_plane==34){
	for(int i=0;i<4;i++){
	  if(channel==i+1){
	    vcsi[i][0] = (Double_t)hit_pla->GetTime();
	    vcsi[i][1] = (Double_t)hit_pla->GetEnergy();
	  }
	}
      }
      if(layer==6 && id_plane==36){
	//	std::cout<<channel<<std::endl;
	for(int i=0;i<4;i++){
	  if(channel==i+1){
	    vpde[i][0] = (Double_t)hit_pla->GetTime();
	    vpde[i][1] = (Double_t)hit_pla->GetEnergy();
	  }
	}
      }
    }//end of the roop npla
    
    /*************NaI from ESPRIRNaI***********/
    for(Int_t i=0;i<nnai;i++){
      TArtRNaI *hit_nai = (TArtRNaI *)nai_array->At(i);
      Int_t id_plane = hit_nai->GetPlaneID();
      Int_t layer = hit_nai->GetLayer();
      Int_t ch = hit_nai->GetCh();
      Int_t id = hit_nai->GetID();

      // for(Int_t j=0;j<28;j++){//NaI UP
      // 	if(j+1==id){
      // 	  if(layer==1){ 
      // 	    vnai[j][1]=hit_nai->GetEnergy();
      // 	    //	    std::cout<<j<<":"<<vnai[j][1]<<std::endl;
      // 	  }//NaI UP
      // 	  if(layer==2){ vnai[j][1]=hit_nai->GetEnergy();}//NaI Down
      // 	}
      // }
      if(layer==1&&id_plane==38){
      	for(Int_t j=0;j<14;j++){//NaI UP
      	  if(j+1==ch){
      	    vnai[j][1]=hit_nai->GetEnergy();
      	  }
      	}
      }
      if(layer==2&&id_plane==39){
      	for(Int_t j=0;j<14;j++){//NaI Down
      	  if(j+1==ch){
      	    vnai[j+14][1]=hit_nai->GetEnergy();
      	  }
      	}
      }

    }//end of the roop nnai
       
    for(int i=0;i<ntdc;i++){
      // vdia_mean[0]=-999999;vdia_mean[1]=-999999;
      // for(Int_t j=0;j<4;j++){
      // 	vdia[0][j]=-999999;vdia[1][j]=-999999;
      // }
      TArtTDCHit *hit_tdc = (TArtTDCHit *)tdc_array->At(i);
      Int_t id_plane = hit_tdc->GetPlaneID();
      Int_t layer = hit_tdc->GetLayer();
      Int_t wireid = hit_tdc->GetWireID();
      Int_t toffset = hit_tdc->GetTzero();
      Int_t fltdc = hit_tdc->GetTDC();
     
      // if(layer==13&&id_plane==97){
      // 	for(Int_t j=0;j<4;j++){
      // 	  if(wireid==j+1) vdia[0][j]=fltdc;
      // 	  vdia_mean[0]+=vdia[0][j];
      // 	} 
      // 	vdia_mean[0]=vdia_mean[0]/4.;
      // }
      // if(layer==17&&id_plane==98){
      // 	for(Int_t j=0;j<4;j++){
      // 	  if(wireid==j+1) vdia[1][j]=fltdc;
      // 	  vdia_mean[1]+=vdia[1][j];
      // 	}     
      // 	vdia_mean[1]=vdia_mean[1]/4.;
      // }

      if(layer==1&&id_plane==38){
	for(int j=0;j<14;j++){
	  if(wireid==j+1) vnai[j][0]=fltdc;
	}
      }
      if(layer==2&&id_plane==39){
	for(int j=0;j<14;j++){
	  if(wireid==j+1) vnai[j+14][0]=fltdc;
	}
      }
      

    }//end of the roop ntdc
       
    /**********bdc from ESPRIBDC   ****************/  
    for(Int_t i=0;i<nbdc;i++){
      TArtBDC *hit_bdc = (TArtBDC *)bdc_array->At(i);
      Int_t layer = hit_bdc->GetLayer();

      vx[layer-1][0] = hit_bdc->GetBDCX();
      vx[layer-1][1] = hit_bdc->GetBDCA();
      vx[layer-1][2] = hit_bdc->GetBDCY();
      vx[layer-1][3] = hit_bdc->GetBDCB();
      // if(vx[layer-1][0]>=0. && vx[layer-1][0]<=10. && vx[layer-1][2]>=0. && vx[layer-1][2]<=10.){
      // 	std::cout<<ne<<":"<<"laye="<<layer<<":"<<vx[layer-1][0]<<":"<<vx[layer-1][1]
      // 		 <<":"<<vx[layer-1][2]<<":"<<vx[layer-1][3]<<std::endl;
      // }
    }
    /*********************************************/

    /**************  analysis here   ******************************/

    /** Trigger & cut setting **/
    for(int i=0;i<2;i++){//N2 trig 
      if(vn2[i][0]>tmin_f12n2 && vn2[i][0]<tmax_f12n2){
	if(vn2[i][1]>qmin_f12n2 && vn2[i][1]<qmax_f12n2){
	  trigN2[i]=1;
	}
      }
    }
    if(vveto[1]>100.){trigveto=0;}//of veto has qdc, no trig 
    /*-----------------------*/

    vtof_f37 = (vdia_mean[1]-vdia_mean[0])/40.; //ToF between f3 & f7 Dia
    vxe12_t = (vxe12[0][0]+vxe12[1][0])/80.; //F12 Xe t
    vn2_t = (vn2[0][0]+vn2[1][0])/80.; //F12 N2 t
    vtof_f712[0] = vxe12_t - vdia_mean[1]/40.;//ToF f7Dia->f12Xe
    vtof_f712[1] = vn2_t - vdia_mean[1]/40.;//Tof f7Dia->f12N2
    vxe12_q = (vxe12[0][1]+vxe12[1][1])/2.; //F12 Xe q
    vbeta[0] = -dis_f712[0]/(vtof_f712[0]*1.E-9*c); //beta of the nuclei
    vbeta[1] = -dis_f712[1]/(vtof_f712[1]*1.E-9*c); //beta of the nuclei

    // std::cout<<vx[0][0]<<":"<<vx[0][2]<<std::endl;
    if(vx[0][0]>=0. && vx[0][0]<=10. && vx[0][2]>=0. && vx[0][2]<=10.){
      if(vx[1][0]>=0. && vx[1][0]<=10. && vx[1][2]>=0. && vx[1][2]<=10.){
	flugtrk = 1;
	nbdcok = nbdcok+1;
	std::cout<<nbdcok<<std::endl;
	//	std::cout<<vx[0][0]<<":"<<vx[0][2]<<":"<<vx[1][0]<<":"<<vx[1][2]<<std::endl;
      }
    }
    if(flugtrk == 1){
      track(vx,vthx,vthy,vxsht,vysht);
      vbdctrk[0] = vthx;
      vbdctrk[1] = vthy;
      vbdctrk[2] = vxsht;
      vbdctrk[3] = vysht;
      //      std::cout<<vx[0][0]<<":"<<vx[0][2]<<":"<<vx[1][0]<<":"<<vx[1][2]<<std::endl;
      //      std::cout<<vthx<<":"<<vthy<<":"<<vxsht<<":"<<vysht<<std::endl;
    }

    //EdE proton PI
    for(int i=0;i<28;i++){
      vnai_ped[i]=vnai[i][1]-vnai_0[i];
    }
    for(int i=0;i<4;i++){
      vpde_ped[i]=vpde[i][1]-vpde_0[i];
    }
    for(int i=0;i<7;i++){//E times dE = constant UP
      if(vnai_ped[2*i]>gnai_min&&vnai_ped[2*i+1]>gnai_min
	 &&vpde_ped[0]>gpde_min&&vpde_ped[1]>gpde_min){//E times dE
	vEtdE[i]=((vnai_ped[2*i]+vnai_ped[2*i+1])/2.
		  +(vpde_ped[0]+vpde_ped[1])/2.)*(vpde_ped[0]+vpde_ped[1])/2.;
      }
    }
    for(int i=7;i<14;i++){//E times dE = constant DOWN
      if(vnai_ped[2*i]>gnai_min&&vnai_ped[2*i+1]>gnai_min
	 &&vpde_ped[2]>gpde_min&&vpde_ped[3]>gpde_min){//E * dE
	vEtdE[i]=((vnai_ped[2*i]+vnai_ped[2*i+1])/2.
		  +(vpde_ped[2]+vpde_ped[3])/2.)*(vpde_ped[2]+vpde_ped[3])/2.;
      }
    }
    /*********************************************/
    
    /*****************************************/
    //FILL to tree!!!
    tana->Fill();//event Fill
    //raw Hist  
    for(int i=0;i<2;i++){
      hraw_t[i]->Fill(vxe12[i][0]); hraw_q[i]->Fill(vxe12[i][1]); 
      hraw_t[i+2]->Fill(vxe8[i][0]);hraw_q[i+2]->Fill(vxe8[i][1]);
      hraw_t[i+4]->Fill(vn2[i][0]); hraw_q[i+4]->Fill(vn2[i][1]); 
      //  hraw_t[i+15]->Fill(vdia[i]);      
    }
    hraw_t[6]->Fill(vveto[0]);  hraw_q[6]->Fill(vveto[1]);
    for(int i=0;i<4;i++){
	hraw_t[i+7]->Fill(vcsi[i][0]); hraw_q[i+7]->Fill(vcsi[i][1]); 
	hraw_t[i+11]->Fill(vpde[i][0]); hraw_q[i+11]->Fill(vpde[i][1]); 
    }
   for(int i=0;i<28;i++){
	hnai_t[i]->Fill(vnai[i][0]); 
	hnai_q[i]->Fill(vnai[i][1]);
   }
   for(int i=0;i<2;i++){//BDC1,BDC2
     for(int j=0;j<4;j++){//X,A,Y,B
       hbdcraw[i*4+j]->Fill(vx[i][j]);   
     }
   }
   
   //Fill to hist when trigger is OK
   if(trigN2[0]*trigN2[1]==1){
     npp[0]=npp[0]+1;//number of trigN2 OK
     for(int i=0;i<4;i++){
       htrk[i]->Fill(vbdctrk[i]);
     }
     hpi[0]->Fill(vtof_f37,vxe12_q);
     hpi[1]->Fill(vtof_f712[0],vxe12_q);
     hpi[2]->Fill(vtof_f712[0],vtof_f37);
     hpi[3]->Fill(vtof_f37*vbeta[1],vxe12_q);
     
     hsht[0]->Fill(vbdctrk[2],vbdctrk[3]);
     
     // if(trigveto==1){//when veto hit, trigveto==0;    
     //   for(int i=0;i<4;i++){
     // 	 htrk_gate[i]->Fill(vbdctrk[i]);
     //   }
     //   hsht_gate[0]->Fill(vbdctrk[2],vbdctrk[3]);
     // }
     
     if(flugtrk==1){
       npp[1]=npp[1]+1;//number of trigN2&&BDCtrk OK
       for(int i=0;i<7;i++){
	 if(vnai[2*i][1]>0&&vnai[2*i+1][1]>0){
	   hEdE[i]->Fill((vnai[2*i][1]+vnai[2*i+1][1])/2.,(vpde[0][1]+vpde[1][1])/2.);
	   hEdE_up->Fill((vnai_ped[2*i]+vnai_ped[2*i+1])/2.,(vpde_ped[0]+vpde_ped[1])/2.);
	   hEtdE_up->Fill(vEtdE[i]);//pedestal gates have already set above
	   if(vEtdE[i]>gEtdE[0][0]&&vEtdE[i]<gEtdE[0][1])
	     npp[2]=npp[2]+1;//proton +1
	   //	  std::cout<<npp[0]<<":"<<npp[1]<<std::endl;
	 }
       }
       for(int i=7;i<14;i++){
	 if(vnai[2*i][1]>0&&vnai[2*i+1][1]>0){
	   hEdE[i]->Fill((vnai[2*i][1]+vnai[2*i+1][1])/2.,(vpde[2][1]+vpde[3][1])/2.);
	   hEdE_down->Fill((vnai_ped[2*i]+vnai_ped[2*i+1])/2.,(vpde_ped[2]+vpde_ped[3])/2.);
	   hEtdE_down->Fill(vEtdE[i]);
	   if(vEtdE[i]>gEtdE[1][0]&&vEtdE[i]<gEtdE[1][1])
	     npp[2]=npp[2]+1;//proton +1
	   //	  std::cout<<npp[0]<<":"<<npp[1]<<std::endl;
	 }
       }
     }
     signal(SIGINT,stop_interrupt); //CTRL + C, interrupt
   }      
   if(ne%100==0){
     printf("Event:%10d/%10d, N2Ok:%8d,trackEvent:%8d,proton:%8d\r"
	    ,ne,neve,npp[0],npp[1],npp[2]);
     fflush(stdout);    
   }
   
  }//end of ne(tespri) roop to neve & fill()


  /*************************************************************************/
  npp_fin[0]=npp[0];npp_fin[1]=npp[1];npp_fin[2]=npp[2];
  tres->Fill();
  /************* make canvas ****************************/
  craw_t->Divide(5,4);
  craw_q->Divide(5,4);
  for(Int_t i=0;i<nhraw_t;i++){
    craw_t->cd(i+1);hraw_t[i]->Draw();
  }
  for(Int_t i=0;i<nhraw_q;i++){
    craw_q->cd(i+1);hraw_q[i]->Draw();
  }
  cnai_t->Divide(7,4);//TDC canvas 7*4
  cnai_q->Divide(7,4);//QDC canvas 7*4
  for(Int_t i=0;i<2;i++){//UP&DOWN
    for(Int_t j=0;j<7;j++){//
      cnai_t->cd(i*14+j+1); hnai_t[i*14+j*2]->Draw();cnai_t->Modified();
      cnai_t->cd(i*14+j+8); hnai_t[i*14+j*2+1]->Draw();cnai_t->Modified(); 
      cnai_q->cd(i*14+j+1); hnai_q[i*14+j*2]->Draw();cnai_q->Modified();
      cnai_q->cd(i*14+j+8); hnai_q[i*14+j*2+1]->Draw();cnai_q->Modified(); 
    }
  }
  cbdcraw->Divide(4,2);
  for(int i=0;i<8;i++){
    cbdcraw->cd(i+1); hbdcraw[i]->Draw();
  }
  ctrk->Divide(2,2);
  for(int i=0;i<4;i++){
    ctrk->cd(i+1);htrk[i]->Draw();
  }

  int npoint=1000;
  for(int k=0;k<npoint;k++){
    hEdEprt1->Fill(nmax_E/npoint*k,EdEfunc1(nmax_E/npoint*k));
    hEdEprt2->Fill(nmax_E/npoint*k,EdEfunc1(nmax_E/npoint*k)-100);
    hEdEprt3->Fill(nmax_E/npoint*k,EdEfunc2(nmax_E/npoint*k)+200);
    hEdEprt4->Fill(nmax_E/npoint*k,EdEfunc3(nmax_E/npoint*k));
  }
  hEdEfun[0]=hEdEprt4;hEdEfun[1]=hEdEprt2;hEdEfun[2]=hEdEprt1;
  hEdEfun[3]=hEdEprt2;hEdEfun[4]=hEdEprt2;hEdEfun[5]=hEdEprt1;
  hEdEfun[6]=hEdEprt2;hEdEfun[7]=hEdEprt2;hEdEfun[8]=hEdEprt2;
  hEdEfun[9]=hEdEprt1;hEdEfun[10]=hEdEprt3;hEdEfun[11]=hEdEprt1;
  hEdEfun[12]=hEdEprt1;hEdEfun[13]=hEdEprt1;
  cEdE->Divide(7,2);
  for(Int_t i=0;i<14;i++){
    cEdE->cd(i+1);gPad->SetLogz(); 
    hEdE[i]->Draw("colz");
    //    hEdEfun[i]->Draw("colz,same");
  }
  /**********************************************/
  std::string output = outdir+"/"+"ana_"+outfile;
  TFile *anafile = new TFile(Form("anaresult/%s",output.c_str()),"RECREATE");

  tana->Write();
  tres->Write();
  for(int i=0;i<nhraw_t;i++){
    hraw_t[i]->Write();
  }
  for(int i=0;i<nhraw_q;i++){
    hraw_q[i]->Write();
  }
  for(int i=0;i<28;i++){
    hnai_t[i]->Write();hnai_q[i]->Write();
  }
  for(int i=0;i<8;i++){
    hbdcraw[i]->Write();
  }
  for(int i=0;i<4;i++){
    htrk[i]->Write();
    // htrk_gate[i]->Write();
  }
  hsht[0]->Write();
  // hsht_gate[0]->Write();
  for(int i=0;i<npi;i++){
    hpi[i]->Write();
  } 
  for(int i=0;i<14;i++){
    hEdE[i]->Write();
  } 
  hEdE_up->Write();
  hEdE_down->Write();
  hEtdE_up->Write();
  hEtdE_down->Write();

  hEdEprt1->Write();
  hEdEprt2->Write();
  hEdEprt3->Write();
  hEdEprt4->Write();
  craw_t->Write();
  craw_q->Write();  
  cnai_t->Write();
  cnai_q->Write();
  cbdcraw->Write();
  ctrk->Write();
  cEdE->Write();
  anafile->Close(); 

  std::cout<<"input = "<<fname<<", evenumber = "<<neve<<std::endl;
  std::cout<<"Number of Tracking = "<<nbdcok<<std::endl;
  std::cout<<"analysis rootfile is created as anaresult/"<<output<<std::endl;
  return 0;
}


