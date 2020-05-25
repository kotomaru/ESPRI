/* macro to from terashima-san plottrack.C*/
#ifndef __CINT__

#include "Riostream.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEllipse.h"
#include "TBox.h"
#include "TApplication.h"

#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtCalibBDC.hh"
#include "TArtBDC.hh"
#include "TArtCalibTDCHit.hh"
#include "TArtTDCHit.hh"
#include "TArtCalibSRPPAC.hh"
#include "TArtSRPPAC.hh"
#include "TArtCalibPlas.hh"
#include "TArtPlas.hh"

#endif
bool stoploop = false;
void stop_interrupt(int sig){
  printf("keybprd interrupt\n");
  stoploop = true;
}
const int tmin_f12n2=-8000; const int tmax_f12n2=-6000;//Sn&Ca
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
    std::cout << "Usage: ./plottrack_pi 0000.root"<<std::endl;
    return -1;
  }
  std::string outdir = fname.substr(9,2);
  std::string outfile = fname.substr(12);



  Double_t vx_f5[2][2];//BLD52(SRPPAC)[51 or 52][X or Y]//Ca run has only 52X  
  Double_t vxe12[2][2]={{-1111,-1111},{-1111,-1111}};//F12Xe raw [left or right] [tdc or adc]
  Double_t vxe12d[2][2]={{-1111,-1111},{-1111,-1111}};//F12XeDelay raw [left or right] [tdc or adc]
  Double_t vn2[2][2]={{-1111,-999999},{-1111,-999999}};//F12N2 raw [left or right] [tdc or adc]
  Double_t vdia[2][4]={-999999,-999999};//Dia timing raw [f3 or f7][pad 1-4] 
  Double_t vdia_mean[2]={-999999,-999999};//Dia timing raw [f3 or f7] 
  //for analysis
  Double_t vtof_f37=-9999;
  Double_t vxe12_q=-1111;

  //const Double_t geom[2][3]={{4.125,4.125,-141.0},{4.125,4.125,-41.0}};//cm
  //Double_t geom[2][3]={{4.125,4.125,-141.0},{4.125,4.125,-41.0}};//cm
  //  Double_t geom[2][3]={{4.125,4.125,-146.0},{4.125,4.125,-46.0}};//cm //RIBF2016
  Double_t geom[2][3]={{4.125,4.125,-142.72},{4.125,4.125,-42.39}};//cm //RIBF2019
  Double_t xx[2][4];//BDC1, BDC2; x,th,y,ph
  Double_t rx[2][4];//RDCu, RDCd; x,th,y,ph
  Double_t scatt[2];//RDCu, RDCd; x,th,y,ph
  const Double_t z_veto=-21.2;
  const Double_t r2d = 57.2958;
  const Double_t focus = 0.;//focus point default=0.0

  Double_t cx1=9945.;Double_t cy1=1015.;
  Double_t rx1=56.;Double_t ry1=52.; 
  TEllipse *cca48_3133 = new TEllipse(cx1,cy1,rx1,ry1,0.,360.,0);  

  TH2F *hpi=new TH2F("hpi","f3f7tofhsi vs f12XeQ",200,9600,10400,100,400,1500);
  TCanvas *ccut = new TCanvas("cut48","cut48");

  for(Int_t i=0;i<2;i++)
    geom[i][2] += focus;
  Int_t period=1.e4;
  Int_t minx=-170,maxx=+30;
  Int_t miny=-5.,maxy=5.;
  TCanvas *cp = new TCanvas("cp","BDC Profile",0,0,500,750);
  TCanvas *cf = new TCanvas("cf","BDC Focus",500,0,500,500);
  TCanvas *cp_pi = new TCanvas("cp_pi","BDC Profile PI",0,0,500,750);
  TCanvas *cf_pi = new TCanvas("cf_pi","BDC Focus PI",500,0,500,500);

  const Double_t by1 = -4., by2 = +4.;
  const Double_t bx1 = geom[0][2]-10.,bx2 = geom[0][2]+10.;
  const Double_t bx3 = geom[1][2]-10.,bx4 = geom[1][2]+10.;
  TBox *bdc1 = new TBox(bx1,by1,bx2,by2);
  TBox *bdc2 = new TBox(bx3,by1,bx4,by2);
  bdc1->SetFillStyle(3013);bdc1->SetFillColor(2);
  bdc2->SetFillStyle(3013);bdc2->SetFillColor(2);
  TEllipse *tgtxy = new TEllipse(0.,0.,1.25,1.25,0,360,0.);//Size 30 mm
  //    TEllipse *tgtxy = new TEllipse(0.,0.,1.5,1.5,0,360,0.);//Size 30 mm
  TEllipse *vetoxy = new TEllipse(0.1,0.,1.3,1.3,0,360,0.);//Size 26 mm
  TEllipse *tgtx = new TEllipse(0.,0.,0.3,1.25,0,360,-45.);//phi30,45deg
  TEllipse *tgty = new TEllipse(0.,0.,0.3,1.25);//phi30,0deg
  tgtxy->SetFillColor(1);
  tgtxy->SetFillStyle(0);
  vetoxy->SetFillColor(1);
  vetoxy->SetFillStyle(0);

  //TH1F *frame = gPad->DrawFrame(minx, miny, maxx, maxy);
  TH1F *hist1 = new TH1F("hist1","X pos (cm)", 100, -5., 5.);
  TH1F *hist2 = new TH1F("hist2","Theta (deg)", 100, -1.4, 1.4);
  TH1F *hist3 = new TH1F("hist3","Y pos (cm)", 100, -5., 5.);
  TH1F *hist4 = new TH1F("hist4","Phi (deg)", 100, -1.4, 1.4);
  TH2F *hist10 = new TH2F("hist10","X-Y at Tgt", 100, -5., 5., 100, -5.0, 5.0);
  //TH2F *hist20 = new TH2F("hist20","A-B at Tgt", 100, -5., 5., 100, -1.4, 1.4);
  TH2F *hist20 = new TH2F("hist20","X-Y at Veto", 100, -5., 5., 100, -5., 5.);
  //TH2F *hist10 = new TH2F("hist10","X-Y by DS(ID1)", 100, -5., 5., 100, -5.0, 5.0);
  //TH2F *hist20 = new TH2F("hist20","X-Y by Recoil(ID2)", 100, -5., 5., 100, -5., 5.);
  //TH2F *hist11 = new TH2F("hist11","X-A at Tgt", 100, -5., 5., 100, -1.4, 1.4);
  //TH2F *hist21 = new TH2F("hist21","Y-B at Tgt", 100, -5., 5., 100, -1.4, 1.4);
  TH2F *hist11 = new TH2F("hist11","X-A at Tgt(mm mrad)", 100, -40., 40., 100, -40.,40.);
  TH2F *hist21 = new TH2F("hist21","Y-B at Tgt(mm mrad)", 100, -40., 40., 100, -40.,40.);
  TH2F *hist12 = new TH2F("hist12","X Track(cm)", 130, -170.,90.,100,-5., 5.);
  TH2F *hist22 = new TH2F("hist22","Y Track(cm)", 130, -170.,90.,100,-5., 5.);

  TH1F *hist_pi1 = new TH1F("hist_pi1","X pos (cm)", 100, -5., 5.);
  TH1F *hist_pi2 = new TH1F("hist_pi2","Theta (deg)", 100, -1.4, 1.4);
  TH1F *hist_pi3 = new TH1F("hist_pi3","Y pos (cm)", 100, -5., 5.);
  TH1F *hist_pi4 = new TH1F("hist_pi4","Phi (deg)", 100, -1.4, 1.4);
  TH2F *hist_pi10 = new TH2F("hist_pi10","X-Y at Tgt", 100, -5., 5., 100, -5.0, 5.0);
  TH2F *hist_pi20 = new TH2F("hist_pi20","X-Y at Veto", 100, -5., 5., 100, -5., 5.);
  TH2F *hist_pi11 = new TH2F("hist_pi11","X-A at Tgt(mm mrad)", 100, -40., 40., 100, -40.,40.);
  TH2F *hist_pi21 = new TH2F("hist_pi21","Y-B at Tgt(mm mrad)", 100, -40., 40., 100, -40.,40.);
  TH2F *hist_pi12 = new TH2F("hist_pi12","X Track(cm)", 130, -170.,90.,100,-5., 5.);
  TH2F *hist_pi22 = new TH2F("hist_pi22","Y Track(cm)", 130, -170.,90.,100,-5., 5.);


  TChain *tespri = new TChain("rtree");
  TString strroot = Form("%s",fname.c_str());
  tespri->Add(strroot);

  TClonesArray * bdc_array=0;
  tespri->SetBranchAddress("ESPRIBDC",&bdc_array);
  TClonesArray * tdc_array=0;
  tespri->SetBranchAddress("ESPRITDC",&tdc_array);
  TClonesArray *srppac_array=0;
  tespri->SetBranchAddress("ESPRISRPPAC",&srppac_array);
  TClonesArray *plas_array=0;
  tespri->SetBranchAddress("ESPRIPLAS",&plas_array);
  
  Int_t nentry = tespri->GetEntries();

  Int_t total=0;
  for(Int_t ev=0; ev< nentry; ev++){
    Int_t trigN2[2]={0};//N2 scinti tdc & qdc ok [left&right]

    if(ev%1000==0){
      cout<<"Event: "<<ev<<"\t"<<100.*ev/nentry<<" %"<<"\r";
      fflush(stdout); 
    }

    tespri->GetEvent(ev);
    //Coincidence Information 
    Int_t ntdc = tdc_array->GetEntries();
    Int_t coin=0;
   
    for(Int_t j=0;j<2;j++){
      vn2[j][0]=-999999;
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
      if(id_plane==40){//CoinBit
	if(ch==1){//changed at 2016/May/17
	  coin=1;//Reaction
	}else if(ch==2){//changed at 2016/May/17
	  coin=2;//DS
	}else if(ch>2){
	  coin=3;//Both
	}else{
	  coin=0;//0 or not defined
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
      if(layer==1&&id_plane==31){
	for(Int_t j=0;j<2;j++){
	  if(wireid=j+1&&fltdc>tmin_f12n2&&fltdc<tmax_f12n2){
	    vn2[j][0]=fltdc;
	  }
	}
      }
    }
    for(Int_t i=0;i<2;i++){
      if(vdia[i][0]>-999990&&vdia[i][1]>-999990
    	 &&vdia[i][2]>-999990&&vdia[i][3]>-999990){
    	vdia_mean[i]=(vdia[i][0]+vdia[i][1]+vdia[i][2]+vdia[i][3])/4.;
      }
    }
    for(Int_t i=0;i<2;i++){//N2 trig 
      if(vn2[i][0]>tmin_f12n2 && vn2[i][0]<tmax_f12n2){
	  trigN2[i]=1;
      }
    }

    /*********  ESPRISRPPAC    *******************************/
    Int_t nsrppac = srppac_array->GetEntries();   
    vx_f5[0][0]=-999;vx_f5[0][1]=-999;vx_f5[1][0]=-999;vx_f5[1][1]=-999;
    for(Int_t i=0;i<nsrppac;i++){
      TArtSRPPAC *hit_srppac = (TArtSRPPAC *)srppac_array->At(i);
      Int_t layer = hit_srppac->GetLayer();
      //      std::cout<<ne<<":"<<i<<":"<<layer<<":"<<hit_srppac->GetSRPPACX()<<std::endl;
     
      vx_f5[layer-1][0] = hit_srppac->GetSRPPACX();
      vx_f5[layer-1][1] = hit_srppac->GetSRPPACY();
    }

    Int_t npla = plas_array->GetEntries(); 
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
 
    Int_t nbdc = bdc_array->GetEntries();
    //cout<<nbdc<<endl;
    for(Int_t i=0;i<nbdc;i++){
      TArtBDC *bdc = (TArtBDC *)bdc_array->At(i);
      Int_t layer = bdc->GetLayer();
      if((bdc->GetBDCX())<99990&&(bdc->GetBDCY())<99990){
	xx[layer-1][0] = bdc->GetBDCX()-geom[layer-1][0];
	xx[layer-1][1] = bdc->GetBDCA();
	xx[layer-1][2] = bdc->GetBDCY()-geom[layer-1][1];
	xx[layer-1][3] = bdc->GetBDCB();
      }
    }//end of array-chain
    Double_t th = (xx[1][0]-xx[0][0])/(geom[1][2]-geom[0][2]);
    Double_t ph = (xx[1][2]-xx[0][2])/(geom[1][2]-geom[0][2]);
    Double_t x0 = (xx[0][0]*geom[1][2]-xx[1][0]*geom[0][2])/(geom[1][2]-geom[0][2]);
    Double_t y0 = (xx[0][2]*geom[1][2]-xx[1][2]*geom[0][2])/(geom[1][2]-geom[0][2]);
    //cout<<x0<<" "<<y0<<" "<<th*180/3.14<<" "<<ph*180/3.14<<endl;
    hist1->Fill(x0);hist2->Fill(th*r2d);hist3->Fill(y0);hist4->Fill(ph*r2d);
    hist10->Fill(x0,y0);//hist20->Fill(th*r2d,ph*r2d);
    hist20->Fill(x0+th*z_veto,y0+ph*z_veto);
    //if(coin==1)
    //hist10->Fill(x0,y0);
    //if(coin==2)
    //hist20->Fill(x0,y0);
    //hist11->Fill(x0,th*r2d);hist21->Fill(y0,ph*r2d);
    hist11->Fill(x0*10,th*1000);hist21->Fill(y0*10,ph*1000);
    for(Int_t i=0;i<100;i++){
      Double_t zz=minx+(maxx-minx)/100*i;
      hist12->Fill(zz,x0+th*zz);hist22->Fill(zz,y0+ph*zz);
    }

    Double_t tmpx; Double_t tmpy;
    vtof_f37 = (vdia_mean[1]-vdia_mean[0]); //ToF between f3 & f7 Dia
    vxe12_q = (vxe12[0][1]-vxe12d[0][1]+vxe12[1][1]-vxe12d[1][1])/2.; //F12 Xe q
    tmpx=vtof_f37-(vx_f5[1][0]-60)*4.5;
    tmpy=vxe12_q;

   if(trigN2[0]*trigN2[1]==1){
     hpi->Fill(vtof_f37-(vx_f5[1][0]-60)*4.5,vxe12_q);     
     if((tmpx-cx1)*(tmpx-cx1)/rx1/rx1+(tmpy-cy1)*(tmpy-cy1)/ry1/ry1-1.<=0){
       hist_pi1->Fill(x0);hist_pi2->Fill(th*r2d);hist_pi3->Fill(y0);hist_pi4->Fill(ph*r2d);
       hist_pi10->Fill(x0,y0);//hist_pi20->Fill(th*r2d,ph*r2d);
       hist_pi20->Fill(x0+th*z_veto,y0+ph*z_veto);
       hist_pi11->Fill(x0*10,th*1000);hist_pi21->Fill(y0*10,ph*1000);
       for(Int_t i=0;i<100;i++){
	 Double_t zz=minx+(maxx-minx)/100*i;
	 hist_pi12->Fill(zz,x0+th*zz);hist_pi22->Fill(zz,y0+ph*zz);
       }
     }
   }
      // if(ev%period==0){
    //   cp->Clear();
    //   //cp->Divide(2,2);
    //   cp->Divide(2,3);
    //   cp->cd(1);hist1->Draw();cp->cd(2);hist2->Draw();
    //   cp->cd(3);hist3->Draw();cp->cd(4);hist4->Draw();
    //   cp->cd(5);hist10->Draw("col");tgtxy->Draw("same");
    //   cp->cd(6);hist20->Draw("col");vetoxy->Draw("same");
    //   cp->Update();

    //   cf->Clear();
    //   cf->Divide(2,2);
    //   cf->cd(1);hist11->Draw("col");cf->cd(3);hist21->Draw("col");
    //   cf->cd(2);hist12->Draw("col"); bdc1->Draw("same"); bdc2->Draw("same"); tgtx->Draw("same");
    //   cf->cd(4);hist22->Draw("col"); bdc1->Draw("same"); bdc2->Draw("same"); tgty->Draw("same");
    //   cf->Update();
    // }
    total++;
  }//end of event loop

  cout << "Total/Ana: "<<tespri->GetEntries()<<"/"<<total<<" Events"<<endl;
  //fout.close();
  //cp->Update();
  //cf->Update();
  cp->Clear();
  cp->Divide(2,3);
  cp->cd(1);hist1->Draw();cp->cd(2);hist2->Draw();
  cp->cd(3);hist3->Draw();cp->cd(4);hist4->Draw();
  cp->cd(5);hist10->Draw("col");tgtxy->Draw("same");
  cp->cd(6);hist20->Draw("col");vetoxy->Draw("same");
  cp->Update();
  
  cf->Clear();
  cf->Divide(2,2);
  cf->cd(1);hist11->Draw("col");cf->cd(3);hist21->Draw("col");
  cf->cd(2);hist12->Draw("col"); bdc1->Draw("same"); bdc2->Draw("same"); tgtx->Draw("same");
  cf->cd(4);hist22->Draw("col"); bdc1->Draw("same"); bdc2->Draw("same"); tgty->Draw("same");
  cf->Update();

  ccut->Clear();
  cca48_3133->SetFillStyle(0);
  hpi->Draw("colz");
  cca48_3133->Draw("same"); 
  ccut->Update();

  cp_pi->Clear();
  cp_pi->Divide(2,3);
  cp_pi->cd(1);hist_pi1->Draw();cp_pi->cd(2);hist_pi2->Draw();
  cp_pi->cd(3);hist_pi3->Draw();cp_pi->cd(4);hist_pi4->Draw();
  cp_pi->cd(5);hist_pi10->Draw("col");tgtxy->Draw("same");
  cp_pi->cd(6);hist_pi20->Draw("col");vetoxy->Draw("same");
  cp_pi->Update();
  
  cf_pi->Clear();
  cf_pi->Divide(2,2);
  cf_pi->cd(1);hist_pi11->Draw("col");cf_pi->cd(3);hist_pi21->Draw("col");
  cf_pi->cd(2);hist_pi12->Draw("col"); bdc1->Draw("same"); bdc2->Draw("same"); tgtx->Draw("same");
  cf_pi->cd(4);hist_pi22->Draw("col"); bdc1->Draw("same"); bdc2->Draw("same"); tgty->Draw("same");
  cf_pi->Update();

  std::string output = outdir+"/"+"trackpi_"+outfile;
  TFile *anafile = new TFile(Form("anaresult/%s",output.c_str()),"RECREATE");
  cp->Write();
  cf->Write();
  cp_pi->Write();
  cf_pi->Write();
  ccut->Write();
  anafile->Close();
  //  cf->WaitPrimitive();
  std::cout<<"analysis rootfile is created as anaresult/"<<output<<std::endl;
  return 0;
}
