//someting wrong -> segmentation violation

void BDCcheck(){
  TCanvas *ct = new TCanvas("ct","fLTDC",1500,900);
  ct->Clear();  
  ct->Divide(4,4); 

  for(Int_t i=0;i<16;i++){ 
    ct->cd(i+1);
    rtree->Draw(((TString::Format("ESPRITDC.fLTDC>>h%i(100,-1000,2000)",i)).Data()),((TString::Format("ESPRITDC.id_plane==%i",i+1)).Data()),"colz",100000);
  }

}
