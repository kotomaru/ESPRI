void plotres(char *name = "BDC"){
  TCanvas *cc1 = new TCanvas("cc1","",800,600);
  cc1->Divide(4,4);  
  if(name=="BDC"){
    for(Int_t i=0;i<8;i++){
      cc1->cd(i+1);
      rtree->Draw(((TString::Format("ESPRI%s.fDrf[0][%i]>>hh%i(200,-0.3,0.3)","BDC",i,i)).Data()),"","colz");
      cc1->cd(i+9);
      rtree->Draw(((TString::Format("ESPRI%s.fDrf[1][%i]>>hh%i(200,-0.3,0.3)","BDC",i,i+8)).Data()),"","colz");
    }
    
    TCanvas *cc = new TCanvas("cc","",800,600);
    cc->Divide(4,4);
    for(Int_t i=0;i<8;i++){
      cc->cd(i+1);
      rtree->Draw(((TString::Format("ESPRI%s.fRes[0][%i]:ESPRIBDC.fDrf[0][%i]>>h%i(200,-0.3,0.3,200,-0.1,0.1)","BDC",i,i,i)).Data()),"","colz");
      cc->cd(i+9);
      rtree->Draw(((TString::Format("ESPRI%s.fRes[1][%i]:ESPRIBDC.fDrf[1][%i]>>h%i(200,-0.3,0.3,200,-0.1,0.1)","BDC",i,i,i+8)).Data()),"","colz");
    }
  }else{
    cout<<"No detector"<<endl;
  }
}

