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
//const int qmin_f12n2=700;  const int qmax_f12n2=2500;//Ca

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
    std::cout << "Usage: ./srppacgain 0000.root"<<std::endl;
    return -1;
  }
  std::string outdir = fname.substr(9,2);
  std::string outfile = fname.substr(12);

  std::string output = outdir+"/"+"srppac_"+outfile;
  TFile *anafile = new TFile(Form("anaresult/%s",output.c_str()),"RECREATE");

  TChain *tespri = new TChain("rtree");
  TString strroot = Form("%s",fname.c_str());
  tespri->Add(strroot);

  //*******Define global variables in variable.h***********
  Double_t vn2_rq[2]={-1111};//F12N2 qdc raw [left&right]
  Double_t vn2_rt[2][50];//F12N2 tdc raw [left&right][multiplicity]
  Double_t vf7dia[4][150]={};//F7Dia LTDC raw [pad1-4] [multi]
  Double_t vbld_L[4][100][10];//Leading [BLD51X,Y,BLD52X,Y][wireid][TDCnumber]
  Double_t vbld_T[4][100][10];//Trailing [BLD51X,Y,BLD52X,Y][wireid][TDCnumber]
  Double_t vbld_tot[4][100][10];//TOT[BLD51X,Y,BLD52X,Y][wire][TDCnum]
  Double_t vn2_t[2];//In beam bunch
  Double_t vf7dia_t[4];//In beam bunch
  Double_t vf7dia_tm;//mean In beam bunch
  Double_t vbld_oL[4][100];//Leading(orderd by TOT) [BLD51X,Y,BLD52X,Y][wireid][TDCnumber]
  Double_t vbld_oT[4][100];//Trailing(orderd by TOT) [BLD51X,Y,BLD52X,Y][wireid][TDCnumber]
  Double_t vbld_otot[4][100];//TOT(ordered)[BLD51X,Y,BLD52X,Y][wire][TDCnum]
  
  Double_t vbld_t0[4][10];//Leading Hit orderd by TOT[BLD51X,Y,BLD52X,Y][TDCnum]
  
  Int_t nmw[4][100]={0};//TDC number of each wire [BLD51X,Y,BLD52X,Y][wireid]
  /******************************************************/
  //result
  Int_t nm[4]={-1};//multiplicity[BLD51X,Y,BLD52X,Y]
  Int_t neff[4][5]={{0},{0},{0},{0}};//Integral of each multiplicity[BDC51XY,BDC52XY][all,1,2,3,4]
  Int_t nespri=0;//Integral of f12N2(LR) is in gate
  Int_t nf7=0;//Integral of f12N2(LR) is in gate
  char name[10];
  
  //*********new TTree for analysis*************/
  TTree *tana = new TTree("tana","tana");
  
  tana->Branch("f12n2_q",&vn2_rq,"f12n2_q[2]/D");//QDC[Left or Right]
  tana->Branch("f12n2_rt",&vn2_rt,"f12n2_rt[2][50]/D");//raw TDC[Left or Right]
  tana->Branch("f12n2_t",&vn2_t,"f12n2_t[2]/D");//TDC[Left or Right]
  tana->Branch("f7dia",&vf7dia,"f7Dia[4][150]/D");//[pad1-4][multi]
  //  tana->Branch("bld5_t",&vbld_t,"bld5_t[4][100][10]/D");//TDC[BLD51X,Y,BLD52X,Y][multi]

  tana->Branch("bld5_L",&vbld_L,"bld5_L[4][100][10]/D");//TDC[BLD51X,Y,BLD52X,Y][multi
  tana->Branch("bld5_T",&vbld_T,"bld5_T[4][100][10]/D");//TDC[BLD51X,Y,BLD52X,Y][multi
  tana->Branch("bld5_tot",&vbld_tot,"bld5_tot[4][100][10]/D");//TDC[BLD51X,Y,BLD52X,Y][multi
  tana->Branch("bld5_oL",&vbld_oL,"bld5_oL[4][100]/D");//TDC[BLD51X,Y,BLD52X,Y][0~93]
  tana->Branch("bld5_oT",&vbld_oT,"bld5_oT[4][100]/D");//TDC[BLD51X,Y,BLD52X,Y][0~93]
  tana->Branch("bld5_otot",&vbld_otot,"bld5_otot[4][100]/D");//TDC[BLD51X,Y,BLD52X,Y][0~93]
  tana->Branch("f7dia_t",&vf7dia_t,"f7Dia_t[4]/D");//first signal in beam bunch[pad]
  tana->Branch("f7dia_tm",&vf7dia_tm,"f7Dia_tm/D");//first signal in beam bunch mean
  tana->Branch("bld5_m",&nm,"nm[4]/I");//multiplicity[BLD51X,Y,BLD52X,Y]
  /**************************************************************/
  TTree *ttrans = new TTree("ttrans","ttrans");

  ttrans->Branch("nespri",&nespri,"nespri/I");//hit number to f12N2
  ttrans->Branch("nf7",&nf7,"nf7/I");//hit number to f12N2
  
  /*******histgram********************************/
  
  
  TH1F *hm[4]; Int_t nbin_m=12; Int_t nmin_m=-1; Int_t nmax_m=11;
  hm[0]=new TH1F("hm_0","multi BLD51X",nbin_m,nmin_m,nmax_m);
  hm[1]=new TH1F("hm_1","multi BLD51Y",nbin_m,nmin_m,nmax_m);
  hm[2]=new TH1F("hm_2","multi BLD52X",nbin_m,nmin_m,nmax_m);
  hm[3]=new TH1F("hm_3","multi BLD52Y",nbin_m,nmin_m,nmax_m);
 
  TCanvas *cm = new TCanvas("cm","multiplicity",1000,500);
  

  /**************************************************/

  TClonesArray *plas_array=0;
  tespri->SetBranchAddress("ESPRIPLAS",&plas_array);
  TClonesArray *tdc_array=0; 
  tespri->SetBranchAddress("ESPRITDC",&tdc_array);
  
  Int_t neve = tespri->GetEntries();
  TArtPlas *hit_pla = 0;
  TArtTDCHit *hit_tdc = 0;

  Int_t tmp1,tmp2,tmp3;  
  /************* Fill data to new tree*******************/
  for(Int_t ne=0;ne<neve;ne++){
    int trigN2[2]={0};//N2 scinti tdc & qdc ok [left&right]
    int trigF7=0;//Dia tdc
    tespri->GetEntry(ne);
    Int_t npla = plas_array->GetEntries(); 
    Int_t ntdc = tdc_array->GetEntries();
    Int_t nmd[4]={0};//Dia multiplicity
    Int_t nmn2[2]={0};//N2 multiplicity
    
    for(Int_t j=0;j<2;j++){
      vn2_t[j]=-999999;
      for(Int_t jj=0;jj<50;jj++){
	vn2_rt[j][jj]=-999999;
      }    
    }
    vf7dia_tm=-999999;
    for(Int_t j=0;j<4;j++){
      for(Int_t jj=0;jj<150;jj++){
	vf7dia[j][jj]=-999999;
      }
      vf7dia_t[j]=-999999;
    }
    for(Int_t j=0;j<4;j++){
      nm[j]=-1;
      for(Int_t jj=0;jj<10;jj++){
	vbld_t0[j][jj]=-999999;
	for(Int_t jjj=0;jjj<100;jjj++){
	  vbld_L[j][jjj][jj]=-999999;
	  vbld_T[j][jjj][jj]=-999999;
	  vbld_tot[j][jjj][jj]=0;
	}
      }   
      for(Int_t jj=0;jj<100;jj++){
	nmw[j][jj]=0;
	vbld_oL[j][jj]=-999999;
	vbld_oT[j][jj]=-999999;
	vbld_otot[j][jj]=-999999;
      }
    }
    
    /**********qdc from ESPRIPLAS          ****************/  
    for(Int_t i=0;i<npla;i++){
      TArtPlas *hit_pla = (TArtPlas *)plas_array->At(i);
      Int_t id_plane = hit_pla->GetPlaneID();
      Int_t layer = hit_pla->GetLayer();
      Int_t channel = hit_pla->GetCh();
      
      if(layer==1 && id_plane==31){
      	if(channel==1){//N2 left no Delay
	  //      	  vn2[0][0] = (Double_t)hit_pla->GetTime();
	  vn2_rq[0] = (Double_t)hit_pla->GetEnergy();
      	}else if(channel==2){//N2 left no Delay
	  //      	  vn2[1][0] = (Double_t)hit_pla->GetTime();
      	  vn2_rq[1]= (Double_t)hit_pla->GetEnergy();
      	}
      }
    }//end of the roop npla
    
    /*************TDC from ESPRITDC***********/
    //    Int_t iok[4] = {0};      
	
    for(Int_t i=0;i<ntdc;i++){
      TArtTDCHit *hit_tdc = (TArtTDCHit *)tdc_array->At(i);
      Int_t id_plane = hit_tdc->GetPlaneID();
      Int_t layer = hit_tdc->GetLayer();
      Int_t wireid = hit_tdc->GetWireID();
      Int_t toffset = hit_tdc->GetTzero();
      Int_t fltdc = hit_tdc->GetTDC();
      Int_t fttdc = hit_tdc->GetTrailTDC();

      if(layer==1&&id_plane==31){
	for(Int_t j=0;j<2;j++){
	  if(wireid=j+1){
	    vn2_rt[j][nmn2[j]]=fltdc;
	    //std::cout<<nmd[j]<<":"<<fltdc<<std::endl;
	    nmn2[j]++; 
	  }
	}
      }
      if(layer==17&&id_plane==98){
	for(Int_t j=0;j<4;j++){
	  if(wireid=j+1){
	    //if(wireid=j+1&&fltdc>8300&&fltdc<10700){
	    vf7dia[j][nmd[j]]=fltdc;
	    nmd[j]++; 	   
	    // vf7dia_t[j]=fltdc;
	    // iok[j]=1;
	  }
	}
      }
      
      for(Int_t ud=0;ud<2;ud++){
    	for(Int_t xy=0;xy<2;xy++){
    	  if(layer==(13+ud*2+xy)&&id_plane==(73+ud*2+xy)){
    	    if(wireid>=1&&wireid<=94){
    	      if(fltdc>-999990){
    		vbld_L[ud*2+xy][wireid-1][nmw[ud*2+xy][wireid-1]]=fltdc;
    		vbld_T[ud*2+xy][wireid-1][nmw[ud*2+xy][wireid-1]]=fttdc;
    		vbld_tot[ud*2+xy][wireid-1][nmw[ud*2+xy][wireid-1]]=fltdc-fttdc;
    		nmw[ud*2+xy][wireid-1]++;
    	      }
    	    }
    	  }
    	}
      }
    }//end of the roop ntdc
    //std::cout<<nmw[2][50]<<std::endl;
    
    for(Int_t i=0;i<4;i++){    
      if(nmd[i]>90){std::cout<<ne<<":"<<i<<":"<<nmd[i]<<std::endl; 
      }
    }
    
    /**************  analysis here   ******************************/
    
    for(Int_t i=0;i<2;i++){
      for(Int_t j=0;j<nmn2[i];j++){
	if(vn2_rt[i][j]>tmin_f12n2&&vn2_rt[i][j]<tmax_f12n2){//in beam bunch
	  vn2_t[i]=vn2_rt[i][j];
	  break;
	}
      } 
    }

    /** Trigger & setting **/
    for(int i=0;i<2;i++){//N2 trig 
      if(vn2_t[i]>tmin_f12n2 && vn2_t[i]<tmax_f12n2){
	//if(vn2[i][1]>qmin_f12n2 && vn2[i][1]<qmax_f12n2){
	trigN2[i]=1;
	//	}
      }
    }
    if(trigN2[0]*trigN2[1]){nespri++;}
    
 
    Int_t iok[4]={0};
    for(Int_t i=0;i<4;i++){
      for(Int_t j=0;j<nmd[i];j++){
    	if(vf7dia[i][j]>8300.&&vf7dia[i][j]<10700.){//in beam bunch
    	  vf7dia_t[i]=vf7dia[i][j];
    	  iok[i]=1;
    	  break;
    	}else{vf7dia_t[i]=-999999;iok[i]=0;}
      } 
    }
    if(iok[0]*iok[1]*iok[2]*iok[3]==1){
      vf7dia_tm=(vf7dia_t[0]+vf7dia_t[1]+vf7dia_t[2]+vf7dia_t[3])/4.;
    }else if(iok[0]+iok[1]+iok[2]+iok[3]==3){
      for(Int_t i=0;i<4;i++){
    	if(iok[0]+iok[1]+iok[2]+iok[3]-iok[i]==3){
    	  vf7dia_tm=(vf7dia_t[0]+vf7dia_t[1]+vf7dia_t[2]+vf7dia_t[3]-vf7dia_t[i])/3.;
    	  break;
    	}
      }
    }else if(iok[0]+iok[1]+iok[2]+iok[3]==2){
      for(Int_t i=0;i<3;i++){
    	for(Int_t j=i+1;j<4;j++){
    	  if(iok[0]+iok[1]+iok[2]+iok[3]-iok[i]-iok[j]==2){
    	    vf7dia_tm=(vf7dia_t[0]+vf7dia_t[1]+vf7dia_t[2]+vf7dia_t[3]
    		       -vf7dia_t[i]-vf7dia_t[j])/2.;
    	    break;
    	  }
    	}
      }
    }else if(iok[0]+iok[1]+iok[2]+iok[3]==1){
      for(Int_t i=0;i<4;i++){
    	if(iok[i]==1){
    	  vf7dia_tm = vf7dia_t[i];
    	  break;
    	}
      }
    }else{vf7dia_tm=-999999;}
    // // if(iok[0]+iok[1]+iok[2]+iok[3]==1){
    // //   std::cout<<ne<<":"<<vf7dia_t[0]<<":"<<vf7dia_t[1]<<":"<<vf7dia_t[2]<<":"<<vf7dia_t[3]<<s
    //    td::endl;
    // }

    if(vf7dia_tm>8300.&&vf7dia_tm<10700.){
      trigF7=1;
      nf7++;
    }

    printf("test\n");
    /** initialization ***/
    for(Int_t i=0;i<4;i++){//BDC51XY,BDC52XY
      for(Int_t ii=0;ii<94;ii++){//wire
	for(Int_t jm=0;jm<nmw[i][ii];jm++){
	  //  Double_t tof_f12f5 = (vn2_t[0]+vn2_t[1])/80. - vbld_L[i][ii][jm]/10.;
	  Double_t tof_f7f5 = vf7dia_tm/40.- vbld_L[i][ii][jm]/10.;
	  // if(tof_f12f5>-224.&&tof_f12f5<-100.&&trigN2[0]*trigN2[1]==1){//bunch+ESPRI Trig
	  //	  if(tof_f7f5>215.&&tof_f7f5<275&&trigF7==1){//bunch + F7 Trig
	  if(tof_f7f5>180.&&tof_f7f5<260&&trigF7==1){//bunch + F7 Trig
	    vbld_oL[i][ii]=vbld_L[i][ii][jm];
	    vbld_oT[i][ii]=vbld_T[i][ii][jm];
	    vbld_otot[i][ii]=vbld_tot[i][ii][jm];
	    //	    std::cout<<ne<<":"<<i<<":"<<ii<<":"<<jm<<std::endl;	   
	    break;
	  }
	}
      }
    }
    /*** Orderd by TOT **/
    for(Int_t i=0;i<4;i++){//BDC51XY,BDC52XY
      for(Int_t ii=0;ii<94;ii++){//wire
	//for(Int_t jm=0;jm<nmw[i][ii];jm++){
	  // if(vbld_t[i][ii-1][jm]>0&&vbld_t[i][ii-1][jm]<700&&
	  //    vbld_t[i][ii-1][jm]<vbld_t0[i][ii-1]){
	for(Int_t j=ii+1;j<94;j++){
	  if(vbld_otot[i][j]>vbld_otot[i][ii]&&vbld_otot[i][j]>0){
	    tmp1=vbld_oL[i][ii];  tmp2=vbld_oT[i][ii];  tmp3=vbld_otot[i][ii]; 
	    vbld_oL[i][ii]=vbld_oL[i][j]; vbld_oL[i][j]=tmp1;
	    vbld_oT[i][ii]=vbld_oT[i][j]; vbld_oT[i][j]=tmp2;
	    vbld_otot[i][ii]=vbld_otot[i][j]; vbld_otot[i][j]=tmp3;
	  }
	}
      }
    }
     
    //    std::cout<<ne<<"    "<<vbld_otot[2][0]<<":"<<vbld_otot[2][1]<<":"<<vbld_otot[2][2]
    //	     <<":"<<vbld_otot[2][3]<<std::endl;

    //multiplivcity of each event
    for(Int_t i=0;i<4;i++){//BDC51XY,BDC52XY
      for(Int_t ii=0;ii<94;ii++){//number orderd by TOT 
	if(vbld_otot[i][ii]>0&&vbld_otot[i][ii]<1000){
	  if(nm[i]==-1){
	    nm[i]=nm[i]+2;
	  }
	  else{
	    nm[i]++;
	  }
	}
      }
    }
    //  std::cout<<nm[0]<<":"<<nm[1]<<":"<<nm[2]<<":"<<nm[3]<<std::endl;
    
    for(Int_t i=0;i<4;i++){
      hm[i]->Fill(nm[i]);
    }
    // }
    for(Int_t i=0;i<4;i++){
      for(Int_t im=1;im<=5;im++){  
	if(nm[i]==im){
	  neff[i][im]++;
	}
      }
      for(Int_t im=1;im<=50;im++){  
	if(nm[i]==im){
	  neff[i][0]++;
	}
      }

    }

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
    
    
  }//end of ne(tespri) roop to neve & fill()
  std::cout<<nf7<<":"<<neff[2][0]<<":"<<neff[2][1]<<":"<<neff[2][2]<<":"<<neff[2][3]<<std::endl;

  /************* make canvas ****************************/
  cm->Divide(2,1);
  TLatex *text[10];
  for(Int_t i=2;i<4;i++){
    sprintf(name, "eff all = %0.4f\n",(Double_t)neff[i][0]/nf7);
    text[(i-2)*5+0]=new TLatex(4,hm[i]->GetMaximum()*0.5,name);
    sprintf(name, "eff m1 = %0.4f\n",(Double_t)neff[i][1]/nf7);
    text[(i-2)*5+1]=new TLatex(4,hm[i]->GetMaximum()*0.4,name);
    sprintf(name, "eff m2 = %0.4f\n",(Double_t)neff[i][2]/nf7);
    text[(i-2)*5+2]=new TLatex(4,hm[i]->GetMaximum()*0.3,name);
    sprintf(name, "eff m3 = %0.4f\n",(Double_t)neff[i][3]/nf7);
    text[(i-2)*5+3]=new TLatex(4,hm[i]->GetMaximum()*0.2,name);
    sprintf(name, "eff m4 = %0.4f\n",(Double_t)neff[i][4]/nf7);
    text[(i-2)*5+4]=new TLatex(4,hm[i]->GetMaximum()*0.1,name);
    for(Int_t j=0;j<5;j++){text[(i-2)*5+j]->SetTextSize(0.05);  text[(i-2)*5+j]->SetTextColor(2);}
  }
  for(int i=0;i<2;i++){
    cm->cd(i+1); hm[i+2]->Draw();
    for(int j=0;j<5;j++){   
      text[i*5+j]->Draw("same");
    } 
  }
  /**********************************************/

  tana->Write();
  ttrans->Fill();
  ttrans->Write();
  for(int i=0;i<4;i++){
    hm[i]->Write();
  }
  cm->Write();
  anafile->Close(); 

  std::cout<<"input = "<<fname<<", evenumber = "<<neve<<std::endl;
  std::cout<<"analysis rootfile is created as anaresult/"<<output<<std::endl;
  return 0;
}


