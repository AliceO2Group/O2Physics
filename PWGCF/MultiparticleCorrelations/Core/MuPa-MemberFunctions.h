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

// a) Called directly in init(...);
// b) Called directly in process(...);
// *) Particle weights;
// *) Utility;

// a) Called directly in init(...):
// void BookBaseList()
// void DefaultConfiguration();
// void DefaultBooking();
// void DefaultBinning();
// void DefaultCuts(); // Remark: has to be called after DefaultBinning(), since some default cuts are defined through default binning, to ease bookeeping
// void BookAndNestAllLists()
// void BookEventHistograms()
// void BookParticleHistograms()
// void BookWeightsHistograms()
// void BookResultsHistograms()

// b) Called directly in process(...):
// void FillEventHistograms(aod::Collision const& collision, aod::Tracks const& tracks, const Int_t rs, const Int_t ba); // reco or sim, before or after event cuts
// Bool_t EventCuts(aod::Collision const& collision)
// void FillParticleHistograms(aod::Track const& track, const Int_t rs, const Int_t ba); // reco or sim, before or after particle cuts
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
  // Each default setting can be overuled e.g. with: task->SetBookParticleHistograms("Phi",kFALSE);
  cph_a.fBookParticleHistograms[ePhi] = kTRUE;
  cph_a.fBookParticleHistograms[ePt] = kTRUE;
  cph_a.fBookParticleHistograms[eEta] = kTRUE;

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
  // task->SetParticleHistogramsBins("Phi",360,0.,TMath::TwoPi());
  cph_a.fParticleHistogramsBins[ePhi][0] = 360;
  cph_a.fParticleHistogramsBins[ePhi][1] = 0.;
  cph_a.fParticleHistogramsBins[ePhi][2] = TMath::TwoPi();
  // task->SetParticleHistogramsBins("Pt",1000,0.,20.);
  cph_a.fParticleHistogramsBins[ePt][0] = 1000;
  cph_a.fParticleHistogramsBins[ePt][1] = 0.;
  cph_a.fParticleHistogramsBins[ePt][2] = 20.;
  // task->SetParticleHistogramsBins("Eta",200,-1.,1.);
  cph_a.fParticleHistogramsBins[eEta][0] = 200;
  cph_a.fParticleHistogramsBins[eEta][1] = -1.;
  cph_a.fParticleHistogramsBins[eEta][2] = 1.;

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
  // *) Control particle histograms;
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
  fEventHistogramsList = new TList();
  fEventHistogramsList->SetName("EventHistograms");
  fEventHistogramsList->SetOwner(kTRUE);
  fBaseList->Add(fEventHistogramsList);

  // *) Control particle histograms:
  fParticleHistogramsList = new TList();
  fParticleHistogramsList->SetName("ParticleHistograms");
  fParticleHistogramsList->SetOwner(kTRUE);
  fBaseList->Add(fParticleHistogramsList);

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

void BookEventHistograms()
{
  // Book all event histograms.

  // a) Book the profile holding flags;
  // b) Book specific event histograms.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fEventHistogramsPro = new TProfile("fEventHistogramsPro", "flags for event histograms", 25, 0., 25.);
  fEventHistogramsPro->SetStats(kFALSE);
  fEventHistogramsPro->SetLineColor(eColor);
  fEventHistogramsPro->SetFillColor(eFillColor);
  // ...
  fEventHistogramsList->Add(fEventHistogramsPro);

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
        ceh_a.fEventHistograms[t][rs][ba] = new TH1D(Form("fEventHistograms[%s][%s][%s]", stype[t].Data(), srs[rs].Data(), sba[ba].Data()), Form("%s, %s, %s", stype[t].Data(), srs[rs].Data(), sba[ba].Data()), (Int_t)ceh_a.fEventHistogramsBins[t][0], ceh_a.fEventHistogramsBins[t][1], ceh_a.fEventHistogramsBins[t][2]);
        ceh_a.fEventHistograms[t][rs][ba]->SetLineColor(fBeforeAfterColor[ba]);
        ceh_a.fEventHistograms[t][rs][ba]->SetFillColor(fBeforeAfterColor[ba] - 10);
        fEventHistogramsList->Add(ceh_a.fEventHistograms[t][rs][ba]);
      } // for(Int_t ba=0;ba<2;ba++)
    }   // for(Int_t rs=0;rs<2;rs++) // reco/sim
  }     // for(Int_t t=0;t<eEventHistograms_N;t++) // type, see enum eEventHistograms

} // void BookEventHistograms()

//============================================================

void BookParticleHistograms()
{
  // Book all particle histograms.

  // a) Book the profile holding flags;
  // b) Book specific particle histograms.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fParticleHistogramsPro = new TProfile("fParticleHistogramsPro", "flags for particle histograms", 25, 0., 25.);
  fParticleHistogramsPro->SetStats(kFALSE);
  fParticleHistogramsPro->SetLineColor(eColor);
  fParticleHistogramsPro->SetFillColor(eFillColor);
  // ...
  fParticleHistogramsList->Add(fParticleHistogramsPro);

  Int_t fBeforeAfterColor[2] = {kRed, kGreen}; //! [0 = kRed,1 = kGreen] TBI 20220713 only temporarily here

  // b) Book specific control particle histograms:
  TString stype[eParticleHistograms_N] = {"Phi", "Pt", "Eta"}; // keep in sync. with enum eParticleHistograms
  TString srs[2] = {"rec", "sim"};
  TString sba[2] = {"before", "after"};

  for (Int_t t = 0; t < eParticleHistograms_N; t++) // type, see enum eParticleHistograms
  {
    if (!cph_a.fBookParticleHistograms[t]) {
      continue;
    }
    for (Int_t rs = 0; rs < 2; rs++) // reco/sim
    {
      for (Int_t ba = 0; ba < 2; ba++) // before/after cuts
      {
        cph_a.fParticleHistograms[t][rs][ba] = new TH1D(Form("fParticleHistograms[%s][%s][%s]", stype[t].Data(), srs[rs].Data(), sba[ba].Data()), Form("%s, %s, %s", stype[t].Data(), srs[rs].Data(), sba[ba].Data()), (Int_t)cph_a.fParticleHistogramsBins[t][0], cph_a.fParticleHistogramsBins[t][1], cph_a.fParticleHistogramsBins[t][2]);
        cph_a.fParticleHistograms[t][rs][ba]->SetLineColor(fBeforeAfterColor[ba]);
        cph_a.fParticleHistograms[t][rs][ba]->SetFillColor(fBeforeAfterColor[ba] - 10);
        fParticleHistogramsList->Add(cph_a.fParticleHistograms[t][rs][ba]);
      } // for(Int_t ba=0;ba<2;ba++)
    }   // for(Int_t rs=0;rs<2;rs++) // reco/sim
  }     // for(Int_t t=0;t<eParticleHistograms_N;t++) // type, see enum eParticleHistograms

} // void BookParticleHistograms()

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

void FillEventHistograms(aod::Collision const& collision, aod::Tracks const& tracks, const Int_t rs, const Int_t ba) // reco or sim, before or after event cuts
{
  // Fill all event histograms.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  ceh_a.fEventHistograms[eNumberOfEvents][rs][ba]->Fill(0.5);
  ceh_a.fEventHistograms[eTotalMultiplicity][rs][ba]->Fill(tracks.size());
  ceh_a.fEventHistograms[eVertex_x][rs][ba]->Fill(collision.posX());
  ceh_a.fEventHistograms[eVertex_y][rs][ba]->Fill(collision.posY());
  ceh_a.fEventHistograms[eVertex_z][rs][ba]->Fill(collision.posZ());

  // TBI 20220808 ctd. with centrality, selected tracks

} // void FillEventHistograms(aod::Collision const& collision, aod::Tracks const& tracks, const Int_t rs, const Int_t ba); // reco or sim, before or after event cuts

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

void FillParticleHistograms(aod::Track const& track, const Int_t rs, const Int_t ba) // reco or sim, before or after particle cuts
{
  // Fill all particle histograms.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  cph_a.fParticleHistograms[ePhi][rs][ba]->Fill(track.phi());
  cph_a.fParticleHistograms[ePt][rs][ba]->Fill(track.pt());
  cph_a.fParticleHistograms[eEta][rs][ba]->Fill(track.eta());

  // TBI 20220808 ctd. with other particle histograms, DCA, etc.

} // void FillParticleHistograms(aod::Track const& track, const Int_t rs, const Int_t ba); // reco or sim, before or after particle cuts

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
