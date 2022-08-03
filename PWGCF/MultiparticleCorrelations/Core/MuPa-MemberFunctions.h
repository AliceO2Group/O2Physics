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
// void DefaultConfiguration();
// void DefaultBooking();
// void DefaultBinning();
// void DefaultCuts(); // Remark: has to be called after DefaultBinning(), since some default cuts are defined through default binning, to ease bookeeping
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

// *) Utility:
// void Red(const char* text);
// void Green(const char* text);
// void Yellow(const char* text);
// void Blue(const char* text);
// TObject* GetObjectFromList(TList *list, Char_t *objectName); // 20220803 NOT_PORTED_YET
// Int_t NumberOfNonEmptyLines(const char *externalFile); // 20220803 NOT_PORTED_YET
// void Exit(const char* functionName, const Int_t lineNumber, const char* message); // 20220803 NOT_PORTED_YET
// void Warning(const char* functionName, const Int_t lineNumber, const char* message); // 20220803 NOT_PORTED_YET

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
  fBasePro->SetLineColor(eColor);
  fBasePro->SetFillColor(eFillColor);
  // ...
  fBaseList->Add(fBasePro);

} // void BookBaseList()

//============================================================

void DefaultConfiguration()
{
  // Default task configuration.

  // task->SetVerbose(kFALSE);
  fVerbose = kFALSE;

} // void DefaultConfiguration()

//============================================================

void DefaultBooking()
{
  // Set here which histograms are booked by default.

  // a) Event histograms;
  // b) Particle histograms;
  // c) QA;

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Event histograms:
  // Each default setting can be overuled e.g. with: task->SetBookEventHistograms("NumberOfEvents",kFALSE);
  ceh_a.fBookEventHistograms[eNumberOfEvents] = kTRUE;
  ceh_a.fBookEventHistograms[eTotalMultiplicity] = kTRUE;
  ceh_a.fBookEventHistograms[eSelectedParticles] = kTRUE;
  ceh_a.fBookEventHistograms[eCentrality] = kTRUE;
  ceh_a.fBookEventHistograms[eVertex_x] = kTRUE;
  ceh_a.fBookEventHistograms[eVertex_y] = kTRUE;
  ceh_a.fBookEventHistograms[eVertex_z] = kTRUE;

  // b) Particle histograms:
  // ...

  // c) QA:
  // ...

} // void DefaultBooking()

//============================================================

void DefaultBinning()
{
  // Default binning for all histograms.

  // a) Default binning for event histograms;
  // b) Default binning for particle histograms;

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Default binning for event histograms:
  // task->SetEventHistogramsBins("NumberOfEvents",1,0.,1.);
  ceh_a.fEventHistogramsBins[eNumberOfEvents][0] = 1;
  ceh_a.fEventHistogramsBins[eNumberOfEvents][1] = 0.;
  ceh_a.fEventHistogramsBins[eNumberOfEvents][2] = 1.;
  // task->SetEventHistogramsBins("TotalMultiplicity",1000,0.,1000.);
  ceh_a.fEventHistogramsBins[eTotalMultiplicity][0] = 1000;
  ceh_a.fEventHistogramsBins[eTotalMultiplicity][1] = 0.;
  ceh_a.fEventHistogramsBins[eTotalMultiplicity][2] = 10000.;
  // task->SetEventHistogramsBins("SelectedParticles",1000,0.,1000.);
  ceh_a.fEventHistogramsBins[eSelectedParticles][0] = 1000;
  ceh_a.fEventHistogramsBins[eSelectedParticles][1] = 0.;
  ceh_a.fEventHistogramsBins[eSelectedParticles][2] = 10000.;
  // task->SetEventHistogramsBins("Centrality",100,0.,100.);
  ceh_a.fEventHistogramsBins[eCentrality][0] = 100;
  ceh_a.fEventHistogramsBins[eCentrality][1] = 0.;
  ceh_a.fEventHistogramsBins[eCentrality][2] = 100.;
  // task->SetEventHistogramsBins("Vertex_x",1000,-20.,20.);
  ceh_a.fEventHistogramsBins[eVertex_x][0] = 1000;
  ceh_a.fEventHistogramsBins[eVertex_x][1] = -20.;
  ceh_a.fEventHistogramsBins[eVertex_x][2] = 20.;
  // task->SetEventHistogramsBins("Vertex_y",1000,-20.,20.);
  ceh_a.fEventHistogramsBins[eVertex_y][0] = 1000;
  ceh_a.fEventHistogramsBins[eVertex_y][1] = -20.;
  ceh_a.fEventHistogramsBins[eVertex_y][2] = 20.;
  // task->SetEventHistogramsBins("Vertex_z",1000,-20.,20.);
  ceh_a.fEventHistogramsBins[eVertex_z][0] = 1000;
  ceh_a.fEventHistogramsBins[eVertex_z][1] = -20.;
  ceh_a.fEventHistogramsBins[eVertex_z][2] = 20.;

  // b) Default binning for particle histograms:
  // ...

} // void DefaultBinning()

//============================================================

void DefaultCuts()
{
  // ...

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

} // void DefaultCuts()

//============================================================

void BookAndNestAllLists()
{
  // *) QA;
  // *) Control event histograms;
  // *) Particle weights;
  // *) Results.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

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

} // void BookAndNestAllLists()

//============================================================

void BookControlEventHistograms()
{
  // Book all control event histograms.

  // a) Book the profile holding flags;
  // b) Book specific control event histograms.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fControlEventHistogramsPro = new TProfile("fControlEventHistogramsPro", "flags for control event histograms", 25, 0., 25.);
  fControlEventHistogramsPro->SetStats(kFALSE);
  fControlEventHistogramsPro->SetLineColor(eColor);
  fControlEventHistogramsPro->SetFillColor(eFillColor);
  // ...
  fControlEventHistogramsList->Add(fControlEventHistogramsPro);

  Int_t fBeforeAfterColor[2] = {kRed, kGreen}; //! [0 = kRed,1 = kGreen] TBI 20220713 only temporarily here

  // b) Book specific control event histograms:
  TString stype[eEventHistograms_N] = {"NumberOfEvents", "TotalMultiplicity", "SelectedParticles", "Centrality", "Vertex_x", "Vertex_y", "Vertex_z"}; // keep in sync. with enum eEventHistograms
  TString srs[2] = {"rec", "sim"};
  TString sba[2] = {"before", "after"};

  for (Int_t t = 0; t < eEventHistograms_N; t++) // type, see enum eEventHistograms
  {
    if (!ceh_a.fBookEventHistograms[t]) {
      continue;
    }
    for (Int_t rs = 0; rs < 2; rs++) // reco/sim
    {
      for (Int_t ba = 0; ba < 2; ba++) // before/after cuts
      {
        // Skip exceptional cases:
        if (eSelectedParticles == t && eBefore == ba) {
          continue;
        } // Number of selected particles makes sense only after cuts
        // ...
        // Book the rest:
        ceh_a.fEventHistograms[t][rs][ba] = new TH1D(Form("fEventHistograms[%s][%s][%s]", stype[t].Data(), srs[rs].Data(), sba[ba].Data()), Form("%s, %s, %s", stype[t].Data(), srs[rs].Data(), sba[ba].Data()), (Int_t)ceh_a.fEventHistogramsBins[t][0], ceh_a.fEventHistogramsBins[t][1], ceh_a.fEventHistogramsBins[t][2]);
        ceh_a.fEventHistograms[t][rs][ba]->SetLineColor(fBeforeAfterColor[ba]);
        ceh_a.fEventHistograms[t][rs][ba]->SetFillColor(fBeforeAfterColor[ba] - 10);
        fControlEventHistogramsList->Add(ceh_a.fEventHistograms[t][rs][ba]);
      } // for(Int_t ba=0;ba<2;ba++)
    }   // for(Int_t rs=0;rs<2;rs++) // reco/sim
  }     // for(Int_t t=0;t<eEventHistograms_N;t++) // type, see enum eEventHistograms

  // *)
  /*
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
  */

} // void BookControlEventHistograms()

//============================================================

void BookWeightsHistograms()
{
  // Book all objects for particle weights.

  // a) Book the profile holding flags;
  // b) Common local labels;
  // c) Histograms.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fWeightsFlagsPro = new TProfile("fWeightsFlagsPro", "flags for particle weights", 3, 0., 3.);
  fWeightsFlagsPro->SetStats(kFALSE);
  fWeightsFlagsPro->SetLineColor(eColor);
  fWeightsFlagsPro->SetFillColor(eFillColor);
  fWeightsFlagsPro->GetXaxis()->SetLabelSize(0.05);
  fWeightsFlagsPro->GetXaxis()->SetBinLabel(1, "w_{#varphi}");
  fWeightsFlagsPro->GetXaxis()->SetBinLabel(2, "w_{p_{t}}");
  fWeightsFlagsPro->GetXaxis()->SetBinLabel(3, "w_{#eta}");
  /*
 for(Int_t w=0;w<eWeights_N;w++) // use weights [phi,pt,eta]
 {
  if(fUseWeights[w])fWeightsFlagsPro->Fill(w+0.5,1.);
 }
 */
  fWeightsList->Add(fWeightsFlagsPro);

  // b) Common local labels: TBI 20220713 book before
  TString sVariable[eWeights_N] = {"#varphi", "p_{t}", "#eta"}; // [phi,pt,eta,rapidity]
  TString sWeights[eWeights_N] = {"w_{#varphi}", "w_{p_{t}}", "w_{#eta}"};

  // c) Histograms:
  for (Int_t w = 0; w < eWeights_N; w++) // use weights [phi,pt,eta]
  {
    //if(!fUseWeights[w]){continue;}
    if (!pw_a.fWeightsHist[w]) // yes, because these histos are cloned from the external ones, see SetWeightsHist(TH1D* const hist, const char *variable)
    {
      //pw_a.fWeightsHist[w] = new TH1D(Form("fWeightsHist[%d]",w),"",(Int_t)fKinematicsBins[w][0],fKinematicsBins[w][1],fKinematicsBins[w][2]);
      pw_a.fWeightsHist[w] = new TH1D(Form("fWeightsHist[%d]", w), "", 200, -100., 100.);
      pw_a.fWeightsHist[w]->SetTitle(Form("Particle weights for %s", sWeights[w].Data()));
      pw_a.fWeightsHist[w]->SetStats(kFALSE);
      pw_a.fWeightsHist[w]->GetXaxis()->SetTitle(sVariable[w].Data());
      pw_a.fWeightsHist[w]->SetFillColor(eFillColor);
      pw_a.fWeightsHist[w]->SetLineColor(eColor);
    }
    fWeightsList->Add(pw_a.fWeightsHist[w]);
  } // for(Int_t w=0;w<eWeights_N;w++) // use weights [phi,pt,eta]

} // void BookWeightsHistograms()

//============================================================

void BookResultsHistograms()
{
  // Book all results histograms.

  // a) Book the profile holding flags;
  // *) ...

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fResultsFlagsPro = new TProfile("fResultsFlagsPro", "flags for results histograms", 1, 0., 1.);
  fResultsFlagsPro->SetStats(kFALSE);
  fResultsFlagsPro->SetLineColor(eColor);
  fResultsFlagsPro->SetFillColor(eFillColor);
  // ...
  fResultsList->Add(fResultsFlagsPro);

  // *)
  fResultsHist = new TH1D("fResultsHist", "...", 10000, -500, 500.);
  fResultsList->Add(fResultsHist);

} // void BookResultsHistograms()

//============================================================

Bool_t EventCuts(aod::Collision const& collision)
{
  // ...

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

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

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  if ((track.pt() < pt_min) || (track.pt() > pt_max)) {
    return kFALSE;
  }

  return kTRUE;

} // void ParticleCuts(aod::Track const& tracks)

//============================================================

void SetWeightsHist(TH1D* const hist, const char* variable)
{

  // Copy histogram holding weights from an external file to the corresponding data member.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

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

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

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

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

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

void Red(const char* text)
{
  cout << "\n\033[1;31m" << text << "\033[0m\n"
       << endl;
}

//============================================================

void Green(const char* text)
{
  cout << "\n\033[1;32m" << text << "\033[0m\n"
       << endl;
}

//============================================================

void Yellow(const char* text)
{
  cout << "\n\033[1;33m" << text << "\033[0m\n"
       << endl;
}

//============================================================

void Blue(const char* text)
{
  cout << "\n\033[1;34m" << text << "\033[0m\n"
       << endl;
}
