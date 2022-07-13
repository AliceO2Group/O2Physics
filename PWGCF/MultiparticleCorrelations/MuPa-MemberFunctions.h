// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// void BookBaseList()
// void BookAndNestAllLists()
// void BookControlEventHistograms()
// void BookWeightsHistograms()
// void BookResultsHistograms()
// Bool_t EventCuts(aod::Collision const& collision)
// Bool_t ParticleCuts(aod::Track const& track)

// *) Particle weights:
// void SetWeightsHist(TH1D* const hist, const char *variable)
// TH1D* GetWeightsHist(const char *variable)
// TH1D* GetHistogramWithWeights(const char *filePath, const char *variable)

//============================================================

void BookBaseList()
{
  // ...
  TList* temp = new TList();
  temp->SetOwner(kTRUE);
  fBaseList.setObject(temp);

  // ...
  fBasePro = new TProfile("fBasePro", "flags for the whole analysis", 15, 0., 15.);
  fBasePro->SetStats(kFALSE);
  fBasePro->SetLineColor(COLOR);
  fBasePro->SetFillColor(FILLCOLOR);
  // ...
  fBaseList->Add(fBasePro);

  return;

} // void BookBaseList()

//============================================================

void BookAndNestAllLists()
{
  // *) QA;
  // *) Control event histograms;
  // *) Particle weights;
  // *) Results.

  //if(fVerbose){Green(__PRETTY_FUNCTION__);}

  // *) QA:
  fQAList = new TList();
  fQAList->SetName("QA");
  fQAList->SetOwner(kTRUE);
  fBaseList->Add(fQAList);

  // *) Control event histograms:
  fControlEventHistogramsList = new TList();
  fControlEventHistogramsList->SetName("ControlEventHistograms");
  fControlEventHistogramsList->SetOwner(kTRUE);
  fBaseList->Add(fControlEventHistogramsList);

  // *) Particle weights:
  fWeightsList = new TList();
  fWeightsList->SetName("Weights");
  fWeightsList->SetOwner(kTRUE);
  fBaseList->Add(fWeightsList);

  // *) Results:
  fResultsList = new TList();
  fResultsList->SetName("Results");
  fResultsList->SetOwner(kTRUE);
  fBaseList->Add(fResultsList);

  return;

} // void BookAndNestAllLists()

//============================================================

void BookControlEventHistograms()
{
  // Book all control event histograms.

  // a) Book the profile holding flags;
  // *) ...

  //if(fVerbose){Green(__PRETTY_FUNCTION__);}

  // a) Book the profile holding flags:
  fControlEventHistogramsPro = new TProfile("fControlEventHistogramsPro", "flags for control event histograms", 25, 0., 25.);
  fControlEventHistogramsPro->SetStats(kFALSE);
  fControlEventHistogramsPro->SetLineColor(COLOR);
  fControlEventHistogramsPro->SetFillColor(FILLCOLOR);
  // ...
  fControlEventHistogramsList->Add(fControlEventHistogramsPro);

  Int_t fBeforeAfterColor[2] = {kRed, kGreen}; //! [0 = kRed,1 = kGreen] TBI 20220713 only temporarily here

  // *)
  for (Int_t ba = 0; ba < 2; ba++) // before/after cuts
  {
    //ceh_a.fMultiplicityHist[ba] = new TH1D(Form("fMultiplicityHist[%d]",ba),Form("%s, %s",fRunNumber.Data(),sba[ba].Data()),(Int_t)fMultiplicityBins[0],fMultiplicityBins[1],fMultiplicityBins[2]);
    ceh_a.fMultiplicityHist[ba] = new TH1D(Form("fMultiplicityHist[%d]", ba), "...", 5000, 0., 5000.);
    ceh_a.fMultiplicityHist[ba]->SetStats(kFALSE);
    ceh_a.fMultiplicityHist[ba]->SetLineColor(fBeforeAfterColor[ba]);
    ceh_a.fMultiplicityHist[ba]->SetFillColor(fBeforeAfterColor[ba] - 10);
    ceh_a.fMultiplicityHist[ba]->GetXaxis()->SetTitle("...");
    fControlEventHistogramsList->Add(ceh_a.fMultiplicityHist[ba]);
  }

  return;

} // void BookControlEventHistograms()

//============================================================

void BookWeightsHistograms()
{
  // Book all objects for particle weights.

  // a) Book the profile holding flags;
  // b) Common local labels;
  // c) Histograms.

  //if(fVerbose){Green(__PRETTY_FUNCTION__);}

  // a) Book the profile holding flags:
  fWeightsFlagsPro = new TProfile("fWeightsFlagsPro", "flags for particle weights", 3, 0., 3.);
  fWeightsFlagsPro->SetStats(kFALSE);
  fWeightsFlagsPro->SetLineColor(COLOR);
  fWeightsFlagsPro->SetFillColor(FILLCOLOR);
  fWeightsFlagsPro->GetXaxis()->SetLabelSize(0.05);
  fWeightsFlagsPro->GetXaxis()->SetBinLabel(1, "w_{#varphi}");
  fWeightsFlagsPro->GetXaxis()->SetBinLabel(2, "w_{p_{t}}");
  fWeightsFlagsPro->GetXaxis()->SetBinLabel(3, "w_{#eta}");
  /*
 for(Int_t w=0;w<gWeights;w++) // use weights [phi,pt,eta]
 {
  if(fUseWeights[w])fWeightsFlagsPro->Fill(w+0.5,1.);
 }
 */
  fWeightsList->Add(fWeightsFlagsPro);

  // b) Common local labels: TBI 20220713 book before
  TString sVariable[gWeights] = {"#varphi", "p_{t}", "#eta"}; // [phi,pt,eta,rapidity]
  TString sWeights[gWeights] = {"w_{#varphi}", "w_{p_{t}}", "w_{#eta}"};

  // c) Histograms:
  for (Int_t w = 0; w < gWeights; w++) // use weights [phi,pt,eta]
  {
    //if(!fUseWeights[w]){continue;}
    if (!pw_a.fWeightsHist[w]) // yes, because these histos are cloned from the external ones, see SetWeightsHist(TH1D* const hist, const char *variable)
    {
      //pw_a.fWeightsHist[w] = new TH1D(Form("fWeightsHist[%d]",w),"",(Int_t)fKinematicsBins[w][0],fKinematicsBins[w][1],fKinematicsBins[w][2]);
      pw_a.fWeightsHist[w] = new TH1D(Form("fWeightsHist[%d]", w), "", 200, -100., 100.);
      pw_a.fWeightsHist[w]->SetTitle(Form("Particle weights for %s", sWeights[w].Data()));
      pw_a.fWeightsHist[w]->SetStats(kFALSE);
      pw_a.fWeightsHist[w]->GetXaxis()->SetTitle(sVariable[w].Data());
      pw_a.fWeightsHist[w]->SetFillColor(FILLCOLOR);
      pw_a.fWeightsHist[w]->SetLineColor(COLOR);
    }
    fWeightsList->Add(pw_a.fWeightsHist[w]);
  } // for(Int_t w=0;w<gWeights;w++) // use weights [phi,pt,eta]

} // void BookWeightsHistograms()

//============================================================

void BookResultsHistograms()
{
  // Book all results histograms.

  // a) Book the profile holding flags;
  // *) ...

  //if(fVerbose){Green(__PRETTY_FUNCTION__);}

  // a) Book the profile holding flags:
  fResultsFlagsPro = new TProfile("fResultsFlagsPro", "flags for results histograms", 1, 0., 1.);
  fResultsFlagsPro->SetStats(kFALSE);
  fResultsFlagsPro->SetLineColor(COLOR);
  fResultsFlagsPro->SetFillColor(FILLCOLOR);
  // ...
  fResultsList->Add(fResultsFlagsPro);

  // *)
  fResultsHist = new TH1D("fResultsHist", "...", 10000, -500, 500.);
  fResultsList->Add(fResultsHist);

  return;

} // void BookResultsHistograms()

//============================================================

Bool_t EventCuts(aod::Collision const& collision)
{
  // ...

  if ((collision.posZ() < Vz_min) || (collision.posZ() > Vz_max)) {
    return kFALSE;
  }
  // ...

  return kTRUE;

} // void EventCuts(aod::Collision const& collision)

//============================================================

Bool_t ParticleCuts(aod::Track const& track)
{
  // ...
  if ((track.pt() < pt_min) || (track.pt() > pt_max)) {
    return kFALSE;
  }

  return kTRUE;

} // void ParticleCuts(aod::Track const& tracks)

//============================================================

void SetWeightsHist(TH1D* const hist, const char* variable)
{

  // Copy histogram holding weights from an external file to the corresponding data member.

  // Basic protection:
  if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))) {
    cout << __LINE__ << endl;
    exit(1);
  }

  Int_t ppe = -1;
  if (TString(variable).EqualTo("phi")) {
    ppe = 0;
  }
  if (TString(variable).EqualTo("pt")) {
    ppe = 1;
  }
  if (TString(variable).EqualTo("eta")) {
    ppe = 2;
  }

  // Finally:
  hist->SetDirectory(0);
  pw_a.fWeightsHist[ppe] = (TH1D*)hist->Clone(); // use eventually this line
  //fWeightsHist = (TH1D*)hist->Clone();
  if (!pw_a.fWeightsHist[ppe]) {
    cout << __LINE__ << endl;
    exit(1);
  } // use eventually this line

  // Flag:
  // fUseWeights[ppe] = kTRUE; // use eventually this line

} // void SetWeightsHist(TH1D* const hist, const char *variable)

//============================================================

TH1D* GetWeightsHist(const char* variable)
{

  // The standard getter.

  // Basic protection:
  if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))) {
    cout << __LINE__ << endl;
    exit(1);
  }

  Int_t ppe = -1;
  if (TString(variable).EqualTo("phi")) {
    ppe = 0;
  }
  if (TString(variable).EqualTo("pt")) {
    ppe = 1;
  }
  if (TString(variable).EqualTo("eta")) {
    ppe = 2;
  }

  // Finally:
  return pw_a.fWeightsHist[ppe]; // use eventually this line
                                 //return fWeightsHist;

} // TH1D* GetWeightsHist(const char *variable)

//============================================================

TH1D* GetHistogramWithWeights(const char* filePath, const char* variable)
{
  // ...

  // a) Return value;
  // b) Basic protection for arguments;
  // c) Check if the external ROOT file exists at specified path;
  // d) Access the external ROOT file and fetch the desired histogram with weights;
  // e) Close the external ROOT file.

  // a) Return value:
  TH1D* hist = NULL;

  // b) Basic protection for arguments:
  if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))) {
    cout << __LINE__ << endl;
    exit(1);
  }

  // c) Check if the external ROOT file exists at specified path:
  /*
 if(gSystem->AccessPathName(filePath,kFileExists))
 {
  //Red(Form("if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s",filePath)); // use eventually
  cout<<__LINE__<<endl;
  exit(1);
 }
 */

  // d) Access the external ROOT file and fetch the desired histogram with weights:
  //TFile *weightsFile = TFile::Open(filePath,"READ");
  TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", "");
  if (!alien) {
    cout << __LINE__ << endl;
    exit(1);
  }
  TFile* weightsFile = TFile::Open(Form("alien://%s", filePath), "READ");

  if (!weightsFile) {
    cout << __LINE__ << endl;
    exit(1);
  }

  hist = (TH1D*)(weightsFile->Get("phi_Task=>0.0-5.0_clone_96"));

  //hist = (TH1D*)(weightsFile->Get(Form("%s_%s",variable,fTaskName.Data()))); // 20220712 this was the original line, instead of thew one above, which is there temporarily

  if (!hist) {
    hist = (TH1D*)(weightsFile->Get(Form("%s", variable)));
  } // yes, for some simple tests I can have only histogram named e.g. 'phi'
  //if(!hist){Red(Form("%s_%s",variable,fTaskName.Data())); cout<<__LINE__<<endl;exit(1);}
  hist->SetDirectory(0);
  hist->SetTitle(filePath);

  // e) Close the external ROOT file:
  weightsFile->Close();
  delete weightsFile;
  weightsFile = NULL;

  return hist;

} // TH1D* GetHistogramWithWeights(const char *filePath, const char *variable)

//============================================================
