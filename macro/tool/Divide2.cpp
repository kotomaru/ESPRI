#ifndef __CINT__
#include <TString.h>
#include <TPad.h>
#include <TStyle.h>
#include <TList.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include "../include/Divide2.h"
#endif 
 
void Divide2(TPad *c,Double_t ratio){
  //ratio : to set the ratio of upper pad xsize to canvas
 
  //divide TPad horizontal
 
  if (!c->IsEditable()) return;
  TVirtualPad* padsav = (TVirtualPad*)gPad;
  c->cd();
 
  const Char_t* name = c->GetName();
 
  Double_t defLeftMargin   = gStyle->GetPadLeftMargin();
  Double_t defBottomMargin = gStyle->GetPadBottomMargin();
  Double_t dx = 1.0;
 
  Double_t x1low = 0;
  Double_t x1up  = x1low + dx;
  Double_t y1up  = 1.0;
  Double_t y1low = y1up - ratio;
 
  TPad* p1 = new TPad(Form("%s_%d",name,1),Form("%s_%d",name,1),x1low,y1low,x1up,y1up);
  p1->SetLeftMargin(defLeftMargin*ratio);
  p1->SetBottomMargin(0.0);
  p1->SetNumber(1);
  p1->Draw();
 
 
  Double_t x2low = 0;
  Double_t x2up  = x1low + dx;
  Double_t y2up  = 1.0 - ratio;
  Double_t y2low = y2up + ratio - 1.0;
 
  TPad *p2 = new TPad(Form("%s_%d",name,2),Form("%s_%d",name,2),x2low,y2low,x2up,y2up);
  p2->SetTopMargin(0.0);
  p2->SetLeftMargin(defLeftMargin*ratio);
  p2->SetBottomMargin(defBottomMargin/(1.0-ratio));
  p2->SetNumber(2);
  p2->Draw();
 
  c->Modified();
 
#ifndef __CINT__
  if (padsav) padsav->cd();
#endif
 
}
 
void SetDivideHistStyle(TPad *c,Double_t ratio = 0.65,Double_t lfactor = 1.00,Double_t tfactor = 1.00,Double_t tfactor2 = 1.00,Double_t tfactor3 = 1.00){
  //lfactor : to resize label font on all histograms by lfactor
  //tfactor : to resize title font on all histograms by tfactor
  //tfactor2 : to resize Ytitle font on bottom histograms by tfactor2
  //tfactor3 : to resize Xtitle font by tfactor3
 
  Double_t defTitleSizex = gStyle->GetTitleSize("x");
  Double_t defLabelSizex = gStyle->GetLabelSize("x");
  Double_t defTitleSizey = gStyle->GetTitleSize("y");
  Double_t defLabelSizey = gStyle->GetLabelSize("y");
 
  TIter next(c->GetListOfPrimitives());
  TObject* obj;
 
  int count = 0;
  TPad* p[2];
  while ( (obj = next()) && count < 2) {
    if(obj->InheritsFrom("TPad")){
      p[count] = (TPad*)obj;
      count++;
    }
  }
 
  if(!(p[0] && p[1])){
    std::cout << "Error: This pad doesn't have two pads" << std::endl;
    return;
  }
 
  TIter next1(p[0]->GetListOfPrimitives());
  TObject* obj1;
 
  while ((obj1 = next1())) {
    TH1* h = NULL;
    if(obj1->InheritsFrom("TH1")){
      h = (TH1*)obj1;
    }else if(obj1->InheritsFrom("TGraph")){
      h = ((TGraph*)obj1)->GetHistogram();
    }else if(obj1->InheritsFrom("TF1")){
      h = ((TF1*)obj1)->GetHistogram();
    }else if(obj1->InheritsFrom("TMultiGraph")){
      h = ((TMultiGraph*)obj1)->GetHistogram();
    }
    if(h){
      Double_t lsizex  = defLabelSizex/ratio*lfactor;
      Double_t tsizex  = defTitleSizex/ratio*tfactor;
      Double_t lsizey  = defLabelSizey/ratio*lfactor;
      Double_t tsizey  = defTitleSizey/ratio*tfactor;
      Double_t offsety = (0.9*lsizey + 0.3*tsizey) / tsizey;
      Double_t lmargin = 1.5*lsizey + 0.9*tsizey;
 
      h->SetLabelSize(lsizex,"x");
      h->SetTitleSize(tsizex,"x");
      h->SetLabelSize(lsizey,"y");
      h->SetTitleSize(tsizey,"y");
      h->SetTitleOffset(offsety*ratio,"y");
 
      if(ratio <= 0.3){
	h->GetYaxis()->SetNdivisions(504);
      }else if(ratio <= 0.4){
	h->GetYaxis()->SetNdivisions(505);
      }else if(0.65 >= ratio){
	h->GetYaxis()->SetNdivisions(509);
      }
 
 
      p[0]->SetLeftMargin(lmargin*ratio);
      p[0]->Modified();
    }
  }
 
  TIter next2(p[1]->GetListOfPrimitives());
  TObject* obj2;
 
  while ((obj2 = next2())) {
    TH1* h = NULL;
    if(obj2->InheritsFrom("TH1")){
      h = (TH1*)obj2;
    }else if(obj2->InheritsFrom("TGraph")){
      h = ((TGraph*)obj2)->GetHistogram();
    }else if(obj2->InheritsFrom("TF1")){
      h = ((TF1*)obj2)->GetHistogram();
    }else if(obj2->InheritsFrom("TMultiGraph")){
      h = ((TMultiGraph*)obj2)->GetHistogram();
    }
    if(h){
      
      Double_t lsizex  = defLabelSizex/(1.0-ratio)*lfactor;
      Double_t tsizex  = defTitleSizex/(1.0-ratio)*tfactor*tfactor3;
      Double_t lsizey  = defLabelSizey/(1.0-ratio)*lfactor;
      Double_t tsizey  = defTitleSizey/(1.0-ratio)*tfactor*tfactor2;
      Double_t offsetx = (0.8*lsizex + 0.4*tsizex) / tsizex;
      Double_t offsety = (0.9*lsizey + 0.3*tsizey/tfactor2 ) / tsizey;
      Double_t bmargin = 1.3*lsizex + 1.3*tsizex;
      Double_t lmargin = 1.5*lsizey + 0.9*tsizey/tfactor2;
 
      h->SetLabelSize(lsizex,"x");
      h->SetTitleSize(tsizex,"x");
      h->SetLabelSize(lsizey,"y");
      h->SetTitleSize(tsizey,"y");
      h->SetTitleOffset(offsetx,"x");
      h->SetTitleOffset(offsety*(1-ratio),"y");
 
      if(ratio >= 0.65){
	h->GetYaxis()->SetNdivisions(504);
      }else if(ratio >= 0.55){
	h->GetYaxis()->SetNdivisions(505);
      }else if(0.40 <= ratio){
	h->GetYaxis()->SetNdivisions(509);
      }
 
      p[1]->SetLeftMargin(lmargin*(1.0-ratio));
      p[1]->SetBottomMargin(bmargin);
      p[1]->Modified();
 
    }
 
  }
 
  c->Modified();
 
}
