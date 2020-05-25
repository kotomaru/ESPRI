void plotstc(char *name = "BDC", Int_t iflag=0){
  //iflag   1: with fit, 0: wi/o fit
  TCanvas *cc = new TCanvas("cc","",800,600);
  cout<<name<<endl;
  Double_t range=0;//fit max
  if(name=="BDC"){
    range=700.;
    cc->Divide(4,4);  
    for(Int_t i=0;i<8;i++){
      cc->cd(i+1);
      rtree->Draw((TString::Format("(ESPRI%s.fRes[0][%i]+ESPRI%s.fDrf[0][%i])/0.25:ESPRI%s.fTDC[0][%i]>>h%i(200,-100,1000,200,0.0,1.1)","BDC",i,"BDC",i,"BDC",i,i).Data()),"","col");
      //      rtree->Draw((TString::Format("(ESPRI%s.fDrf[0][%i])/0.25:ESPRI%s.fTDC[0][%i]>>h%i(200,-100,1000,200,0.0,1.1)","BDC",i,"BDC",i,i).Data()),"","col");
      cc->cd(i+9);
      rtree->Draw((TString::Format("(ESPRI%s.fRes[1][%i]+ESPRI%s.fDrf[1][%i])/0.25:ESPRI%s.fTDC[1][%i]>>h%i(200,-100,1000,200,0.0,1.1)","BDC",i,"BDC",i,"BDC",i,i+8).Data()),"","col");
      //rtree->Draw((TString::Format("(ESPRI%s.fDrf[1][%i])/0.25:ESPRI%s.fTDC[1][%i]>>h%i(200,-100,1000,200,0.0,1.1)","BDC",i,"BDC",i,i+8).Data()),"","col");
    }
  }else if(name=="RDC"){
    cc->Divide(4,4);  
    for(Int_t i=0;i<7;i++){
      cc->cd(i+1);
      rtree->Draw((TString::Format("(ESPRI%s.fRes[0][%i]+ESPRI%s.fDrf[0][%i])/0.7:ESPRI%s.fTDC[0][%i]>>h%i(200,-100,1000,200,0.0,1.1)","RDC",i,"RDC",i,"RDC",i,i).Data()),"","col");
      cc->cd(i+9);
      rtree->Draw((TString::Format("(ESPRI%s.fRes[1][%i]+ESPRI%s.fDrf[1][%i])/0.7:ESPRI%s.fTDC[1][%i]>>h%i(200,-100,1000,200,0.0,1.1)","RDC",i,"RDC",i,"RDC",i,i+8).Data()),"","col");
    }
  }else{
    cout<<"No Detector"<<endl;
  }
  if(iflag==1){
    Double_t t0=0.,t1=0.,t2=0.,t2=0.,t3=0.,t4=0.;
    TCanvas *ccf = new TCanvas("ccf","",800,600);
    for(Int_t i=0;i<16;i++){
      TH2F *hh = (TH2F*)gROOT->FindObject(TString::Format("h%i",i).Data());
      hh->ProfileX("hh_pfx");
      hh->Draw("col");
      hh_pfx->Fit("pol4","Q","");
      TF1 *fit0 = (TF1 *)hh_pfx->GetFunction("pol4");
      Double_t par[5];
      fit0->GetParameters(par);
      TF1 *fit = new TF1("fit","[0]*sqrt(x)+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",-100,range);
      fit->SetParameter(0,fit0->GetParameter(1));
      hh_pfx->Fit("fit","Q","",0,range);
      Double_t par0,par1,par2,par3,par4;
      par0=fit->GetParameter(0);
      par1=fit->GetParameter(1);
      par2=fit->GetParameter(2);
      par3=fit->GetParameter(3);
      par4=fit->GetParameter(4);
      t0 += par0; t1 += par1; t2 += par2; t3 += par3; t4 += par4;
      cout<<i<<"\t"<<par0<<"\t"<<par1<<"\t"<<par2<<"\t"<<par3<<"\t"<<par4<<endl;
      cc->cd(i+1);
      fit->SetLineColor(4);
      fit->Draw("same");
      ccf->cd();
    }
    cout<<"av:\t"<<t0/16.<<"\t"<<t1/16.<<"\t"<<t2/16.<<"\t"<<t3/16.<<"\t"<<t4/16.<<endl;
    ccf->Close();
  }
}

