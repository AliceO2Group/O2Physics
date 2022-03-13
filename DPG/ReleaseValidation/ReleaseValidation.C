#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

TFile* fileOut=0x0;
TFile* fileSummaryOutput=0x0;


TString prefix="";
Int_t correlationCase=0; // at the moment...Later I need to check that 

void ProcessFile(TString fname, TString dirToAnalyse); //needs to be declared
void ProcessMonitorObjectCollection(o2::quality_control::core::MonitorObjectCollection* o2MonObjColl);
void ProcessMonitorObject(o2::quality_control::core::MonitorObject *MonObj);
void WriteHisto(TObject* obj);
void WriteHisto2D(TObject* obj);
void WriteProfile(TObject* obj);
void WriteTEfficiency(TObject* obj);
Bool_t AreIdenticalHistos(TH1* hA, TH1* hB);
Bool_t CompareChiSquareAndBinContent(TH1* hA, TH1* hB, Int_t whichTest);
void CompareHistos(TH1* hA, TH1* hB, TString monobj, Int_t whichTest, Bool_t FinalComparison);
void DrawRatio(TH1* hR);
void DrawRelativeDifference(TH1* hR);
void SelectCriticalHistos(TString whichdir);

//void DrawSummaryCanvas(TFile *file);

void ReleaseValidation(TString filename1="QCpass3", TString filename2="QCpass2", TString ObjectToAnalyse = "", Int_t whichTest = 1, Bool_t isOnGrid = kFALSE, Bool_t SelectCritical = kFALSE){

  if(whichTest != 1 && whichTest != 2 && whichTest != 3) {
    printf("Error, please select which test you want to perform: \n");
    printf("1->Chi-square; 2--> ContBinDiff; 3 --> Chi-square and ContBinDiff");
    return;
  }
  TCanvas* cpdf = new TCanvas("cpdf","cpdf");
  cpdf->Print("plots.pdf[");
  fileSummaryOutput = new TFile("Summary_" + ObjectToAnalyse + ".root","recreate");
  fileSummaryOutput->Close();
 
  TFile* inFile1 = 0x0;
  TFile* inFile2 = 0x0;
  //open files to be read in
  if(isOnGrid) {
    TGrid::Connect("alien://");
     inFile1 = TFile::Open("alien:" + filename1 + ".root", "READ");
     inFile2 = TFile::Open("alien:" + filename2 + ".root", "READ");
     inFile1->ls();
     inFile2->ls();
   } else {
     inFile1 = TFile::Open(filename1 + ".root", "READ");
     inFile2 = TFile::Open(filename2 + ".root", "READ");
   }
   
   inFile1->ls();
   inFile2->ls();

   fileOut=new TFile("newfile1.root","recreate");
   ProcessFile(filename1+ ".root",ObjectToAnalyse);
   fileOut->Close();
   delete fileOut;   // do I need to remove it?? //
   fileOut=new TFile("newfile2.root","recreate");
   ProcessFile(filename2+ ".root",ObjectToAnalyse);
   fileOut->Close();
   delete fileOut;
   
   TFile* fileA=new TFile("newfile1.root");
   TFile* fileB=new TFile("newfile2.root");

   Bool_t IsLastComparison = kFALSE;  // It is kTRUE when the last histogram of the file is considered, on order to properly close the pdf
   
   Int_t nkeys=fileA->GetNkeys();
   TList* lkeys=fileA->GetListOfKeys();
   for(Int_t j=0; j<nkeys; j++){
     if(j ==nkeys-1 ) IsLastComparison = kTRUE;
     TKey* k=(TKey*)lkeys->At(j);
     TString cname=k->GetClassName();
     TString oname=k->GetName();
     if(cname.BeginsWith("TH")){
       TH1* hA=static_cast<TH1*>(fileA->Get(oname.Data()));
       TH1* hB=static_cast<TH1*>(fileB->Get(oname.Data()));
       printf("%s and %s compared \n", hA->GetName(),  hB->GetName());
       if(hA && hB){
	 Bool_t ok=AreIdenticalHistos(hA,hB);
	 if(ok) printf("%s       ---> IDENTICAL\n",oname.Data());
	 else{
	   CompareHistos(hA,hB,ObjectToAnalyse,whichTest, IsLastComparison);
	   //Bool_t ok=CompareHistos(hA,hB, ObjectToAnalyse, IsLastComparison);
	   //if(ok) printf("%s       ---> COMPATIBLE\n",oname.Data());
	   //  else printf("********  %s     ---> BAD *****\n",oname.Data());
	 }
       }else{
	 if(!hA) printf("%s    ---> MISSING in first file\n",hA->GetName());
	 if(!hB) printf("%s    ---> MISSING  in second file\n",hA->GetName());
       }
     }
   }

   //fileSummaryOutput = new TFile("Summary" + ObjectToAnalyse  + ".root","READ");
   if(SelectCritical)  SelectCriticalHistos(ObjectToAnalyse);
   //DrawSummaryCanvas(fileSummaryOutput); 
   return;
   
}

void ProcessFile(TString fname, TString dirToAnalyse){

  TFile *fileBase=TFile::Open(fname.Data());
  Int_t nkeys=fileBase->GetNkeys();
  TList* lkeys=fileBase->GetListOfKeys();
  for(Int_t j=0; j<nkeys; j++){
    prefix="";
    TKey* k=(TKey*)lkeys->At(j);
    TString cname=k->GetClassName();
    TString oname=k->GetName();
    printf("****** KEY %d: %s (class %s)   ******\n",j,oname.Data(),cname.Data());
    if(cname=="o2::quality_control::core::MonitorObjectCollection"){
      o2::quality_control::core::MonitorObjectCollection* o2MonObjColl =(o2::quality_control::core::MonitorObjectCollection*)fileBase->Get(oname.Data());
      TString objName=o2MonObjColl->GetName();
      if(dirToAnalyse.Length()>0 && !objName.Contains(dirToAnalyse.Data())){
	printf("Skip MonitorObjectCollection %s\n",objName.Data());
      }else{
	prefix.Append(Form("%s_",o2MonObjColl->GetName()));
	//	fileSummaryOutput = new TFile(Form("Summary_%s.root",o2MonObjColl->GetName()),"recreate");
	ProcessMonitorObjectCollection(o2MonObjColl);
	//DrawSummaryCanvas(fileSummaryOutput);
      }
    }/*else if(cname=="TList"){
      TList* ls=(TList*)fileBase->Get(oname.Data());
      TString lsName=ls->GetName();
      if(lsName=="cstatsout") ProcessPhysSel(ls);
      if(dirToAnalyse.Length()>0 && !lsName.Contains(dirToAnalyse.Data())){
	printf("Skip TList %s\n",lsName.Data());
      }else{
	prefix.Append(Form("%s_",ls->GetName()));
	ProcessList(ls);
      }
    }else if(cname.BeginsWith("TH")){
      TObject* o=(TObject*)fileBase->Get(oname.Data());
      WriteHisto(o);
    }*/
  }
  return;
}

void ProcessMonitorObjectCollection(o2::quality_control::core::MonitorObjectCollection* o2MonObjColl){
  printf("--- Process o2 Monitor Object Collection %s ---\n",o2MonObjColl->GetName());
  Int_t nent=o2MonObjColl->GetEntries();
  Int_t Counts = 0;
  for(Int_t j=0; j<nent; j++){
    TObject* o=(TObject*)o2MonObjColl->At(j);
    TString clname=o->ClassName();
    TString olname=o->GetName();
    printf("****** %s (class %s)   ******\n",olname.Data(),clname.Data());
    if(clname=="o2::quality_control::core::MonitorObject"){
     o2::quality_control::core::MonitorObject * MonObj=(o2::quality_control::core::MonitorObject*)o2MonObjColl->FindObject(olname.Data());
     //prefix.Append(Form("%s_",MonObj->GetName()));	
     ProcessMonitorObject(MonObj);
     Counts ++;
     if(Counts == 40) break;   // to avoid crushes 
      /*
	}else if(cname=="TObjArray"){
	TObjArray* arr=(TObjArray*)df->Get(oname.Data());
	prefix.Append(Form("%s_",arr->GetName()));	
	//ProcessArray(arr);
	
	}else if(cname.Contains("THnSpars")){
	printf("   --> DIRECTORY %s CONTAINS A THNSPARSE: %s\n",df->GetName(),oname.Data());
	}else if(cname.Contains("TNtuple")){
	printf("   --> DIRECTORY %s CONTAINS A TNtuple: %s\n",df->GetName(),oname.Data());      
	}else if(cname.Contains("TTree")){
	printf("   --> DIRECTORY %s CONTAINS A TTree: %s\n",df->GetName(),oname.Data());
      */      
    }else if(clname.BeginsWith("TH")){
      TObject* o=(TObject*)o2MonObjColl->FindObject(olname.Data());
      WriteHisto(o);
      // }else if(cname.BeginsWith("TProfile")){
      // TObject* o=(TObject*)df->Get(oname.Data());
      //WriteProfile(o);
    }
  }
  printf("%d objects processed \n",Counts);
  return; 
}

void ProcessMonitorObject(o2::quality_control::core::MonitorObject* MonObj){
  printf("------ Process o2 Monitor Object %s ------\n",MonObj->GetName());
  TObject* Obj = (TObject*)MonObj->getObject();
  TString Clname=Obj->ClassName();
  TString Objname=Obj->GetName();
  printf("****** %s (class %s)   ******\n",Objname.Data(),Clname.Data());
  if(Clname.BeginsWith("TH")) { 
    WriteHisto(Obj);
  } //else if(Clname.BeginsWith("TH2")) { 
    //WriteHisto2D(Obj);
  //}
  /*
  else if(Clname=="o2::quality_control::core::MonitorObject"){
    o2::quality_control::core::MonitorObject * SecondMonObj=(o2::quality_control::core::MonitorObject*)MonObj->FindObject(Clname.Data());
    printf("------- more Monitor Objects found---\n");
    prefix.Append(Form("%s_",SecondMonObj->GetName()));
    ProcessMonitorObject(SecondMonObj);
    }*/ else if (Clname.BeginsWith("TProfile")) {
    WriteProfile(Obj);
    //prefix.Append(Form("%s_",Obj->GetName()));
  } else if (Clname == "TEfficiency") {
    WriteTEfficiency(Obj);
    //prefix.Append(Form("%s_",Obj->GetName()));
  } else printf("class %s needs to be analysed \n", Clname.Data());
  return;
}

void WriteHisto(TObject* obj){
  TH1* hA=static_cast<TH1*>(obj);
  TString hAcln=hA->ClassName();
  TDirectory *current=gDirectory;
  TCanvas *cc = new TCanvas();
  if(hAcln.Contains("TH2")) hA->Draw("colz");
  else hA->DrawNormalized();
  cc->SaveAs(Form("%s.png",hA->GetName()));
  fileOut->cd();
  hA->Write(Form("%s%s",prefix.Data(),hA->GetName()));
  current->cd();
  return;
}

void WriteHisto2D(TObject* obj){
  TH2* hA2D=static_cast<TH2*>(obj);
  TDirectory *current=gDirectory;
  new TCanvas();
  hA2D->Draw("colz");
  fileOut->cd();
  hA2D->Write(Form("%s%s",prefix.Data(),hA2D->GetName()));
  current->cd();
  return;
}

void WriteProfile(TObject* obj) { //should I further develop that?
  TProfile* hProf = static_cast<TProfile*>(obj); 
  new TCanvas();
  // hProf->Write(Form("%s%s",prefix.Data(),hProf->GetName()));
  hProf->Draw();
  return;
}

void WriteTEfficiency(TObject* obj) { //should I further develop that?
  TEfficiency* hEff = static_cast<TEfficiency*>(obj); 
  new TCanvas();
  hEff->Write(Form("%s%s",prefix.Data(),hEff->GetName()));
  hEff->Draw("A4");
  return;
}
  
Bool_t AreIdenticalHistos(TH1* hA, TH1* hB){

  Long_t nEventsA=hA->GetEntries(); //temporary: we should use the number of analyzed events
  Long_t nEventsB=hB->GetEntries(); //temporary: we should use the number of analyzed events
  if(nEventsA!=nEventsB) {
    printf(" %s -> Different number of entries: A --> %ld, B --> %ld\n", hA->GetName(), nEventsA, nEventsB);
    return kFALSE;
  }
  for(Int_t ix=1; ix<=hA->GetNbinsX(); ix++){
    for(Int_t iy=1; iy<=hA->GetNbinsY(); iy++){
      for(Int_t iz=1; iz<=hA->GetNbinsZ(); iz++){
	Double_t cA=hA->GetBinContent(ix,iy,iz);
	Double_t cB=hB->GetBinContent(ix,iy,iz);
	if(TMath::Abs(cA-cB)>0.001*TMath::Abs(cA)) return kFALSE;
      }
    }
  }
  return kTRUE;
}
 

void CompareHistos(TH1* hA, TH1* hB, TString monobj, Int_t whichTest, Bool_t FinalComparison){
  
  Int_t nEventsA=hA->GetEntries();
  Int_t nEventsB=hB->GetEntries();
  //if(nEventsA==0 && nEventsB==0) return kTRUE;
  //if(nEventsA==0 && nEventsB!=0) return kFALSE;
  //if(nEventsA!=0 && nEventsB==0) return kFALSE;

  TH1 *hACl = (TH1*)hA->Clone("hACl");   // I will use these two clones of hA and hB to perform other checks later..
  TH1 *hBCl = (TH1*)hB->Clone("hBCl");

  TString outc="";
  Int_t colt=1;
  Bool_t ok = kTRUE;
  if(ok == CompareChiSquareAndBinContent(hA, hB, whichTest)) {
      outc="Compatible";
      colt=kGreen+1;
    } else{
      outc="BAD";
      colt=2;
    }
  
  /*
  Double_t chi2=0;
  Double_t meandiff=0;
  Int_t nBins=0;
  for(Int_t ix=1; ix<=hA->GetNbinsX(); ix++){
    for(Int_t iy=1; iy<=hA->GetNbinsY(); iy++){
      for(Int_t iz=1; iz<=hA->GetNbinsZ(); iz++){
	Double_t cA=hA->GetBinContent(ix,iy,iz);
	Double_t eA=0;
	if(cA<0) printf("Negative counts!!! cA=%f in bin %d %d %d\n",cA,ix,iy,iz);
	else eA=TMath::Sqrt(cA);
	Double_t cB=hB->GetBinContent(ix,iy,iz);
	Double_t eB=0;
	if(cB<0) printf("Negative counts!!! cB=%f in bin %d %d %d\n",cB,ix,iy,iz);
	else eB=TMath::Sqrt(cB);
	Double_t diff=cA-cB;
	if(cA>0 && cB>0){
	  Double_t correl=0.;
	  if(correlationCase==1){
	    // estimate degree of correlation from number of events in histogram
	    // assume that the histogram with less events is a subsample of that with more events
	    correl=TMath::Sqrt(cA/cB);
	    if(cA>cB) correl=1./correl;
	  }
	  Double_t sigma2=eA*eA+eB*eB-2*correl*eA*eB; // to be improved
	  meandiff+=(cA/nEventsA-cB/nEventsB);
	  if(sigma2>0) chi2+=diff*diff/sigma2;
	  nBins++;
	}
      }
    }
  }
  */
  // if(nBins>1){
  //  printf(" -> Different contents: %s  chi2/nBins=%f   meanreldiff=%f \n",hA->GetName(),chi2/nBins,meandiff);
    TCanvas* c = new TCanvas(hA->GetName(),hA->GetName(),1200,600);
    c->Divide(2,1);
    c->cd(1);
    TString hcln=hA->ClassName();
    TString optD="";
    if(hcln.Contains("TH2")) optD="box";
    hA->SetLineColor(1);
    hA->SetMarkerColor(1);
    hA->Scale(1./hA->GetEntries());  // normalize to the number of entries
    TH1* hAc=(TH1*)hA->DrawClone(optD.Data());
    hB->SetLineColor(2);
    hB->SetMarkerColor(2);
    hB->Scale(1./hB->GetEntries()); // normalize to the number of entries
    TH1* hBc=(TH1*)hB->DrawClone(Form("%ssames",optD.Data()));
    //c->Update();
    TPaveStats* stA=(TPaveStats*)hAc->GetListOfFunctions()->FindObject("stats");
    if(stA){
      stA->SetLineColor(1);
      stA->SetTextColor(1);
      stA->SetY1NDC(0.68);
      stA->SetY2NDC(0.88);
    }
    TPaveStats* stB=(TPaveStats*)hBc->GetListOfFunctions()->FindObject("stats");
    if(stB){
      stB->SetLineColor(2);
      stB->SetTextColor(2);
      stB->SetY1NDC(0.45);
      stB->SetY2NDC(0.65);
    }
    
    c->cd(2);
    if(hcln.Contains("TH3")){
      TH1D* hXa=((TH3*)hA)->ProjectionX(Form("%s_xA",hA->GetName()));
      TH1D* hXb=((TH3*)hB)->ProjectionX(Form("%s_xB",hB->GetName()));
      TH1D* hYa=((TH3*)hA)->ProjectionY(Form("%s_yA",hA->GetName()));
      TH1D* hYb=((TH3*)hB)->ProjectionY(Form("%s_yB",hB->GetName()));
      TH1D* hZa=((TH3*)hA)->ProjectionZ(Form("%s_zA",hA->GetName()));
      TH1D* hZb=((TH3*)hB)->ProjectionZ(Form("%s_zB",hB->GetName()));
      hXa->Divide(hXb);
      hYa->Divide(hYb);
      hZa->Divide(hZb);
      TPad* rpad=(TPad*)gPad;
      rpad->Divide(1,3);
      rpad->cd(1);
      DrawRatio(hXa);
      rpad->cd(2);
      DrawRatio(hYa);
      rpad->cd(3);
      DrawRatio(hZa);


    }else{
      TH1 *hArat = (TH1*)hA->Clone("hArat");
      hArat->Divide(hB);
      for(Int_t k=1;k<=hArat->GetNbinsX();k++) hArat->SetBinError(k,0.000000001);
      hArat->SetMinimum(TMath::Max(0.98,0.95*hArat->GetBinContent(hArat->GetMinimumBin())-hArat->GetBinError(hArat->GetMinimumBin())));
      hArat->SetMaximum(TMath::Min(1.02,1.05*hArat->GetBinContent(hArat->GetMaximumBin())+hArat->GetBinError(hArat->GetMaximumBin())));
      hArat->SetStats(0);
      if(hcln.Contains("TH2")) hArat->Draw("colz");
      else if(hcln.Contains("TH1")) DrawRatio(hArat);
      else hArat->Draw();

    }
    c->cd(1);
    /*
    TString outc="";
    Int_t colt=1;
    Bool_t retVal=kTRUE;
    if(chi2/nBins<1.5 && TMath::Abs(meandiff)<0.01){
      outc="Compatible";
      colt=kGreen+1;
    }else{
      outc="BAD";
      colt=2;
      retVal=kFALSE;
    }
    */
    TLatex* toutc=new TLatex(0.2,0.85,outc.Data());
    toutc->SetNDC();
    toutc->SetTextColor(colt);
    toutc->SetTextFont(62);
    toutc->Draw();
    
    c->SaveAs(Form("%s.png",hA->GetName()));
    fileSummaryOutput = new TFile("Summary_" + monobj + ".root","update");
    c->Write(Form("%s%s_Ratio",prefix.Data(),hA->GetName()));
    fileSummaryOutput->ls();
    fileSummaryOutput->Close();
    c->Print("plots.pdf[");
    /*
    if(FinalComparison)
      {
	c->Print("plots.pdf[");
	//c->Print("plots.pdf]");
      }
    else c->Print("plots.pdf[");
    */
    //return retVal;


    // Implement the plotting of the difference between the two histograms, and the relative difference 
    TCanvas* c1 = new TCanvas(hACl->GetName(),hACl->GetName(),1200,600);
    c1->Divide(2,1);
    c1->cd(1);
    
    TString hAClcln=hACl->ClassName();
    TString NoptD="";
    if(hAClcln.Contains("TH2")) NoptD="colz";   //box
    hACl->SetLineColor(1);
    hACl->SetMarkerColor(1);
    hACl->Scale(1./hACl->GetEntries());
    hBCl->Scale(1./hBCl->GetEntries());

    // Subtraction
    TH1 * hDiff = (TH1*)hACl->Clone("hDiff");
    hDiff->Add(hBCl, -1);
    hDiff->DrawClone(NoptD.Data());
    
    TPaveStats* stACl=(TPaveStats*)hACl->GetListOfFunctions()->FindObject("stats");
    if(stACl){
      stACl->SetLineColor(1);
      stACl->SetTextColor(1);
      stACl->SetY1NDC(0.68);
      stACl->SetY2NDC(0.88);
    }
    
    c1->cd(2);
    if(hcln.Contains("TH3")){
      TH1D* hXaCl=((TH3*)hDiff)->ProjectionX(Form("%s_xA",hACl->GetName()));
      TH1D* hXbCl=((TH3*)hBCl)->ProjectionX(Form("%s_xB",hBCl->GetName()));
      TH1D* hYaCl=((TH3*)hDiff)->ProjectionY(Form("%s_yA",hACl->GetName()));
      TH1D* hYbCl=((TH3*)hBCl)->ProjectionY(Form("%s_yB",hBCl->GetName()));
      TH1D* hZaCl=((TH3*)hDiff)->ProjectionZ(Form("%s_zA",hACl->GetName()));
      TH1D* hZbCl=((TH3*)hBCl)->ProjectionZ(Form("%s_zB",hBCl->GetName()));
      hXaCl->Divide(hXbCl);
      hYaCl->Divide(hYbCl);
      hZaCl->Divide(hZbCl);
      TPad* rrpad=(TPad*)gPad;
      rrpad->Divide(1,3);
      rrpad->cd(1);
      DrawRelativeDifference(hXaCl);
      rrpad->cd(2);
      DrawRelativeDifference(hYaCl);
      rrpad->cd(3);
      DrawRelativeDifference(hZaCl);


    }else{
      TH1 * hDiffRel = (TH1*)hDiff->Clone("hDiffRel");
      hDiffRel->Divide(hBCl);
      for(Int_t k=1;k<=hDiffRel->GetNbinsX();k++) hDiffRel->SetBinError(k,0.000000001);
      hDiffRel->SetMinimum(TMath::Max(0.98,0.95*hDiffRel->GetBinContent(hDiffRel->GetMinimumBin())-hDiffRel->GetBinError(hDiffRel->GetMinimumBin())));
      hDiffRel->SetMaximum(TMath::Min(1.02,1.05*hDiffRel->GetBinContent(hDiffRel->GetMaximumBin())+hDiffRel->GetBinError(hDiffRel->GetMaximumBin())));
      hDiffRel->SetStats(0);
      TString hDiffRelcln=hDiffRel->ClassName();
      if(hDiffRelcln.Contains("TH2")) hDiffRel->Draw("colz");
      else if(hDiffRelcln.Contains("TH1")) DrawRelativeDifference(hDiffRel);
      else hDiffRel->Draw();

    }
    /*
    c1->cd(1);
    TString outc="";
    Int_t colt=1;
    Bool_t retVal=kTRUE;
    if(chi2/nBins<1.5 && TMath::Abs(meandiff)<0.01){
      outc="Compatible";
      colt=kGreen+1;
    }else{
      outc="BAD";
      colt=2;
      retVal=kFALSE;
    }
    TLatex* toutc=new TLatex(0.2,0.85,outc.Data());
    toutc->SetNDC();
    toutc->SetTextColor(colt);
    toutc->SetTextFont(62);
    */
    c1->cd(1);
    toutc->Draw();
    c1->SaveAs(Form("%s_Difference.png",hDiff->GetName()));
    fileSummaryOutput = new TFile("Summary_" + monobj + ".root","update");
    c1->Write(Form("%s%s_Difference",prefix.Data(),hA->GetName()));
    fileSummaryOutput->ls();
    fileSummaryOutput->Close();
    if(FinalComparison)
      {
	c1->Print("plots.pdf[");
	c1->Print("plots.pdf]");
      }
    else c1->Print("plots.pdf[");
   
    //return retVal;
    //}
    return;
    //return kTRUE;
}

Bool_t CompareChiSquareAndBinContent(TH1* hA, TH1* hB, Int_t whichTest){
  // implement here some simple checks that the two histograms are statistically compatible
  Int_t nEventsA=hA->GetEntries();
  Int_t nEventsB=hB->GetEntries();
  TString oname=hA->GetEntries();

  if(nEventsA==0 && nEventsB==0) return kTRUE;
  if(nEventsA==0 && nEventsB!=0) return kFALSE;
  if(nEventsA!=0 && nEventsB==0) return kFALSE;
  Double_t chi2=0;
  Double_t meandiff=0;
  Int_t nBins=0;
  for(Int_t ix=1; ix<=hA->GetNbinsX(); ix++){
    for(Int_t iy=1; iy<=hA->GetNbinsY(); iy++){
      for(Int_t iz=1; iz<=hA->GetNbinsZ(); iz++){
	Double_t cA=hA->GetBinContent(ix,iy,iz);
	Double_t eA=0;
	if(cA<0) printf("Negative counts!!! cA=%f in bin %d %d %d\n",cA,ix,iy,iz);
	else eA=TMath::Sqrt(cA);
	Double_t cB=hB->GetBinContent(ix,iy,iz);
	Double_t eB=0;
	if(cB<0) printf("Negative counts!!! cB=%f in bin %d %d %d\n",cB,ix,iy,iz);
	else eB=TMath::Sqrt(cB);
	Double_t diff=cA-cB;
	if(cA>0 && cB>0){
	  Double_t correl=0.;
	  if(correlationCase==1){
	    // estimate degree of correlation from number of events in histogram
	    // assume that the histogram with less events is a subsample of that with more events
	    correl=TMath::Sqrt(cA/cB);
	    if(cA>cB) correl=1./correl;
	  }
	  Double_t sigma2=eA*eA+eB*eB-2*correl*eA*eB; // to be improved
	  meandiff+=(cA/nEventsA-cB/nEventsB);
	  if(sigma2>0) chi2+=diff*diff/sigma2;
	  nBins++;
	}
      }
    }
  }
  if(nBins>1){
    printf(" -> Different contents: %s  chi2/nBins=%f   meanreldiff=%f \n",hA->GetName(),chi2/nBins,meandiff);
    Bool_t retVal = kTRUE;
    switch (whichTest) {
    case 1:
      printf("chi-square test performed. \n");
      if(chi2/nBins<1.5) {
	printf("%s       ---> COMPATIBLE\n",oname.Data());
	retVal = kTRUE;
      }else{
	printf("%s       ---> BAD\n",oname.Data());
	retVal=kFALSE;
      }
      break;
      
    case 2:
      printf("bin-content test performed. \n");
      if(TMath::Abs(meandiff)<0.01) {
	printf("%s       ---> COMPATIBLE\n",oname.Data());
	retVal = kTRUE;
      }else{
	printf("%s       ---> BAD\n",oname.Data());
	retVal=kFALSE;
      }
      break;
      
    case 3:
      printf("chi-square and bin-content test performed. \n");
      if(chi2/nBins<1.5 && TMath::Abs(meandiff)<0.01) {
	printf("%s       ---> COMPATIBLE\n",oname.Data());
	retVal = kTRUE;
      }else{
	printf("%s       ---> BAD\n",oname.Data());
	retVal=kFALSE;
      }
      break;
    }
    
    return retVal; 
  }
  return kTRUE;
}

void DrawRatio(TH1* hR){
  hR->SetMarkerStyle(20);
  hR->SetMarkerSize(0.5);
  hR->SetMinimum(TMath::Max(0.98,0.95*hR->GetBinContent(hR->GetMinimumBin())-hR->GetBinError(hR->GetMinimumBin())));
  hR->SetMaximum(TMath::Min(1.02,1.05*hR->GetBinContent(hR->GetMaximumBin())+hR->GetBinError(hR->GetMaximumBin())));
  hR->SetStats(0);
  hR->GetYaxis()->SetTitle("Ratio");
  hR->Draw("P");
  return;
}

void DrawRelativeDifference(TH1* hR){
  hR->SetMarkerStyle(20);
  hR->SetMarkerSize(0.5);
  hR->SetMinimum(TMath::Max(0.98,0.95*hR->GetBinContent(hR->GetMinimumBin())-hR->GetBinError(hR->GetMinimumBin())));
  hR->SetMaximum(TMath::Min(1.02,1.05*hR->GetBinContent(hR->GetMaximumBin())+hR->GetBinError(hR->GetMaximumBin())));
  hR->SetStats(0);
  hR->GetYaxis()->SetTitle("RelativeDifference");
  hR->Draw("P");
  return;
}
  
void  SelectCriticalHistos(TString whichdir) {
  printf("Select all critical plots..... \n");
  
  vector<string> NamesFromTheList;
  fileSummaryOutput = new TFile("Summary_" + whichdir  + ".root","READ");
  fileSummaryOutput->ls();

  ifstream InputFile;
  InputFile.open("CriticalPlots.txt");
  string string;
  while(!InputFile.eof()) // To get all the lines
    {
      std::getline(InputFile,string);  // Save the names in a string
      NamesFromTheList.push_back(string); // Save the histo names in the string vector
      cout << string << endl;
    }
  InputFile.close();

  // access the string vector elements
  std::cout << "Access the elements of the list of critical..." << std::endl;
  for(Int_t i = 0; i < NamesFromTheList.size(); i++) {
    cout << NamesFromTheList[i] << endl; 
  }
  TCanvas* critic_pdf = new TCanvas("critic_pdf","critic_pdf");
  critic_pdf->Print("critical.pdf[");
  
  Int_t Nkeys=fileSummaryOutput->GetNkeys();
  std::cout << "In the summary file there are " << Nkeys << " plots. \n " << std::endl;
  TList*Lkeys=fileSummaryOutput->GetListOfKeys();
  for(Int_t j=0; j<Nkeys; j++){
    std::cout << "case " << j << std::endl; 
    TKey* k=(TKey*)Lkeys->At(j);
    TString Cname=k->GetClassName();
    TString Oname=k->GetName();
    std::cout << Oname << " " << Cname << std::endl;
    for(Int_t i = 0; i < NamesFromTheList.size(); i++) {
      std::cout << NamesFromTheList[i] << std::endl;
      if(Oname.Contains(NamesFromTheList[i]) && NamesFromTheList[i] != "" ) {
	std::cout << " name file and name from the list: "<< Oname << " e " << NamesFromTheList[i] << std::endl;
	TCanvas* ccc = static_cast<TCanvas*>(fileSummaryOutput->Get(Oname.Data()));
	//ccc->Draw();
	ccc->Print("critical.pdf[");
      }
    }
  }
  critic_pdf->Print("critical.pdf]");
  // access the string vector elements
  //std::cout << "Access the elements of the list of critical..." << std::endl;
  //for(Int_t i = 0; i < NamesFromTheList.size(); i++) {
  //  cout << NamesFromTheList[i] << endl; 
  // }
  
  return;
}


/*void DrawSummaryCanvas(TFile* file){
  file->ls();
  Int_t Nkeys=file->GetNkeys();
  std::cout << Nkeys << std::endl;
  TList*Lkeys=file->GetListOfKeys();
  TCanvas *Sum = new TCanvas("Sum","Sum");
  Sum->Divide(2,Nkeys/2);
  for(Int_t j=0; j<Nkeys; j++){
    TKey* k=(TKey*)Lkeys->At(j);
    TString Cname=k->GetClassName();
    TString Oname=k->GetName();
    Sum->cd(j);
    if(Cname.BeginsWith("TCanvas")){
      TCanvas* ccc = static_cast<TCanvas*>(file->Get(Oname.Data()));
      ccc->Draw();
    }
  }
  Sum->SaveAs("Sum.pdf");
  printf("***************************** \n");
  return;
}
*/
