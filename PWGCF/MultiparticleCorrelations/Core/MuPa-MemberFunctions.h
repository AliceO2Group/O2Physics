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
// void BookQvectorHistograms()
// void BookCorrelationsHistograms()
// void BookWeightsHistograms()
// void BookResultsHistograms()

// b) Called directly in process(...):
// void ResetEventByEventQuantities();
// void FillEventHistograms(aod::Collision const& collision, aod::Tracks const& tracks, const Int_t rs, const Int_t ba); // reco or sim, before or after event cuts
// Bool_t EventCuts(aod::Collision const& collision)
// void FillParticleHistograms(aod::Track const& track, const Int_t rs, const Int_t ba); // reco or sim, before or after particle cuts
// Bool_t ParticleCuts(aod::Track const& track)
// void CalculateCorrelations();

// *) Q-vectors:
// TComplex Q(Int_t n, Int_t p);
// TComplex One(Int_t n1);
// TComplex Two(Int_t n1, Int_t n2);
// ... TBI 20220809 port the rest ...
// void ResetQ(); // reset the components of generic Q-vectors

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

  // task->SetCalculateQvector(kTRUE);
  fCalculateQvector = kTRUE;

  // task->SetCalculateCorrelations(kTRUE);
  fCalculateCorrelations = kTRUE;

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
  // *) Correlations;
  // *) Q-vectors;
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

  // *) Q-vectors:
  fQvectorList = new TList();
  fQvectorList->SetName("Q-vectors");
  fQvectorList->SetOwner(kTRUE);
  fBaseList->Add(fQvectorList);

  // *) Correlations:
  fCorrelationsList = new TList();
  fCorrelationsList->SetName("Correlations");
  fCorrelationsList->SetOwner(kTRUE);
  fBaseList->Add(fCorrelationsList);

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

void BookQvectorHistograms()
{
  // Book all Q-vector histograms.

  // a) Book the profile holding flags;
  // b) ...

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fQvectorFlagsPro = new TProfile("fQvectorFlagsPro", "flags for Q-vector objects", 3, 0., 3.);
  fQvectorFlagsPro->SetStats(kFALSE);
  fQvectorFlagsPro->SetLineColor(eColor);
  fQvectorFlagsPro->SetFillColor(eFillColor);
  fQvectorFlagsPro->GetXaxis()->SetLabelSize(0.05);
  fQvectorFlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateQvector");
  fQvectorFlagsPro->Fill(0.5, fCalculateQvector);
  fQvectorFlagsPro->GetXaxis()->SetBinLabel(2, "gMaxHarmonic");
  fQvectorFlagsPro->Fill(1.5, gMaxHarmonic);
  fQvectorFlagsPro->GetXaxis()->SetBinLabel(3, "gMaxCorrelator");
  fQvectorFlagsPro->Fill(2.5, gMaxCorrelator);
  fQvectorList->Add(fQvectorFlagsPro);

  // b) ...

} // void BookQvectorHistograms()

//============================================================

void BookCorrelationsHistograms()
{
  // Book all correlations histograms.

  // a) Book the profile holding flags;
  // b) Common local labels;
  // c) Histograms.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fCorrelationsFlagsPro = new TProfile("fCorrelationsFlagsPro", "flags for correlations", 3, 0., 3.);
  fCorrelationsFlagsPro->SetStats(kFALSE);
  fCorrelationsFlagsPro->SetLineColor(eColor);
  fCorrelationsFlagsPro->SetFillColor(eFillColor);
  fCorrelationsFlagsPro->GetXaxis()->SetLabelSize(0.05);
  fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateCorrelations");
  fCorrelationsFlagsPro->Fill(0.5, fCalculateCorrelations);
  // ...
  fCorrelationsList->Add(fCorrelationsFlagsPro);

  if (!fCalculateCorrelations) {
    return;
  }

  // b) Common local labels:
  TString oVariable[4] = {"#varphi_{1}-#varphi_{2}", "#varphi_{1}+#varphi_{2}-#varphi_{3}-#varphi_{4}",
                          "#varphi_{1}+#varphi_{2}+#varphi_{3}-#varphi_{4}-#varphi_{5}-#varphi_{6}",
                          "#varphi_{1}+#varphi_{2}+#varphi_{3}+#varphi_{4}-#varphi_{5}-#varphi_{6}-#varphi_{7}-#varphi_{8}"};

  TString vvVariable[3] = {"int", "mult", "cent"};

  // c) Histograms:
  for (Int_t k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]
  {
    for (Int_t n = 0; n < 6; n++) // harmonic [n=1,n=2,...,n=6]
    {
      for (Int_t v = 0; v < 3; v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
      {
        // ... TBI 20220809 ... port the rest

        // c_a.fCorrelationsPro[k][n][v] = new TProfile(Form("c_a.fCorrelationsPro[%d][%d][%s]",k,n,vvVariable[v].Data()),harmonicArray.Data(),vvvariableNBins[v],vvvariableMinMax[v][0],vvvariableMinMax[v][1]);
        c_a.fCorrelationsPro[k][n][v] = new TProfile(Form("fCorrelationsPro[%d][%d][%s]", k, n, vvVariable[v].Data()), "some title", 2000, 0., 2000.);
        c_a.fCorrelationsPro[k][n][v]->SetStats(kFALSE);
        c_a.fCorrelationsPro[k][n][v]->Sumw2();
        c_a.fCorrelationsPro[k][n][v]->GetXaxis()->SetTitle(vvVariable[v].Data());
        c_a.fCorrelationsPro[k][n][v]->GetYaxis()->SetTitle(Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1), oVariable[k].Data()));
        fCorrelationsList->Add(c_a.fCorrelationsPro[k][n][v]);
      }
    }
  }

} // BookCorrelationsHistograms()

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

void ResetEventByEventQuantities()
{
  // Reset all global event-by-event quantities here:

  // a) Multiplicities;
  // b) Centrality;
  // c) Q-vectors;
  // d) Reset ebe containers for nested loops;
  // e) Fisher-Yates algorithm.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Multiplicities:
  //fMultiplicity = 0.;
  //fSelectedParticles = 0;

  // b) Centrality:
  //fCentrality = 0.;

  // c) Q-vectors:
  if (fCalculateQvector) {
    ResetQ(); // generic Q-vector
    for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
      for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) // weight power
      {
        qv_a.fQvector[h][wp] = TComplex(0., 0.);
      }
    }
  } // if(fCalculateQvector)

  // ... TBI 20220809 port the rest ...

} // void ResetEventByEventQuantities()

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

void CalculateCorrelations()
{
  // Calculate analytically multiparticle correlations from Q-vectors
  // In this method, only isotropic correlations for which all harmonics are the same are evaluated.
  // For the calculus of generic multiparticle correlations, see method CalculateGenericCorrelations()

  // a) Flush 'n' fill the generic Q-vectors;
  // b) Calculate correlations;
  // c) Flush the generic Q-vectors.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  // a) Flush 'n' fill the generic Q-vectors:
  ResetQ();
  for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
    for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) // weight power
    {
      qv_a.fQ[h][wp] = qv_a.fQvector[h][wp];
    }
  }

  // b) Calculate correlations:
  for (Int_t h = 1; h <= gMaxHarmonic; h++) // harmonic
  {
    // 2p:
    //if(fSelectedParticles<2){return;}
    //if(fVerbose){cout<<Form("   => CalculateCorrelations(void), 2p, h = %d .... ",h)<<endl;}
    TComplex two = Two(h, -h);
    Double_t twoC = two.Re(); // cos
    //Double_t twoS = two.Im(); // sin
    Double_t wTwo = Two(0, 0).Re(); // Weight is 'number of combinations' by default TBI 20220809 add support for other weights
    if (wTwo > 0.0) {
      twoC /= wTwo;
    } else {
      return;
    } // ... TBI 20220809 ... use the line below eventually
    // else { Exit(__PRETTY_FUNCTION__,__LINE__,Form("wTwo = %f is not positive. fSelectedParticles = %d",wTwo,fSelectedParticles)); return; }

    /* ... TBI 20220809 ... enable eventually
  if(fCalculateCustomNestedLoop)
  {
   // e-b-e sanity check:
   TArrayI *harmonics = new TArrayI(2);
   harmonics->SetAt(h,0);
   harmonics->SetAt(-h,1);
   Double_t nestedLoopValue = this->CalculateCustomNestedLoop(harmonics);
   if(TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(twoC - nestedLoopValue)>1.e-5)
   {
    Exit(__PRETTY_FUNCTION__,__LINE__,Form("nestedLoopValue = %f is not the same as twoC = %f",nestedLoopValue,twoC)); exit(1);
   }
   else
   {
    cout<<Form("=> e-b-e check with CustomNestedLoop is OK for isotropic 2-p, harmonic %d",h)<<endl;
    //cout<<Form("   value = %f",twoC)<<endl;
   }
   delete harmonics; harmonics = NULL;
  } // if(fCalculateCustomNestedLoop)

  // for on-the-fly and internal validation, rescale results with theoretical value:
  if(fCalculateOnTheFly && fOnTheFlyFlowAmplitudes && fRescaleWithTheoreticalInput &&
     TMath::Abs(fOnTheFlyFlowAmplitudes->GetAt(h-1))>0.){twoC/=pow(fOnTheFlyFlowAmplitudes->GetAt(h-1),2.);}
  else if(fUseInternalValidation && fInternalValidationAmplitudes && fRescaleWithTheoreticalInput &&
          TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){twoC/=pow(fInternalValidationAmplitudes->GetAt(h-1),2.);}

*/

    // integrated:
    if (c_a.fCorrelationsPro[0][h - 1][0]) {
      c_a.fCorrelationsPro[0][h - 1][0]->Fill(0.5, twoC, wTwo);
    }
    // vs. multiplicity:
    //if(c_a.fCorrelationsPro[0][h-1][1]){c_a.fCorrelationsPro[0][h-1][1]->Fill(fSelectedParticles+0.5,twoC,wTwo);}
    // vs. centrality:
    //if(c_a.fCorrelationsPro[0][h-1][2]){c_a.fCorrelationsPro[0][h-1][2]->Fill(fCentrality,twoC,wTwo);}

    // ... TBI 20220809 port the rest ...

  } // for(Int_t h=1;h<=gMaxHarmonic;h++) // harmonic

  // c) Flush the generic Q-vectors:
  ResetQ();

} // void CalculateCorrelations()

//============================================================

TComplex Q(Int_t n, Int_t wp)
{
  // Using the fact that Q{-n,p} = Q{n,p}^*.

  if (n >= 0) {
    return qv_a.fQ[n][wp];
  }
  return TComplex::Conjugate(qv_a.fQ[-n][wp]);

} // TComplex FlowWithMultiparticleCorrelationsTask::Q(Int_t n, Int_t wp)

//============================================================

TComplex One(Int_t n1)
{
  // Generic expression <exp[i(n1*phi1)]>.

  TComplex one = Q(n1, 1);

  return one;

} // TComplex FlowWithMultiparticleCorrelationsTask::One(Int_t n1)

//============================================================

TComplex Two(Int_t n1, Int_t n2)
{
  // Generic two-particle correlation <exp[i(n1*phi1+n2*phi2)]>.

  TComplex two = Q(n1, 1) * Q(n2, 1) - Q(n1 + n2, 2);

  return two;

} // TComplex FlowWithMultiparticleCorrelationsTask::Two(Int_t n1, Int_t n2)

//============================================================

// ... TBI 20220809 ... port the rest

//============================================================

void ResetQ()
{
  // Reset the components of generic Q-vectors. Use it whenever you call the standard functions for correlations, for some custom Q-vectors.

  if (fVerbose) {
    Green(__PRETTY_FUNCTION__);
  }

  for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
    for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) // weight power
    {
      qv_a.fQ[h][wp] = TComplex(0., 0.);
    }
  }

} // void ResetQ()

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
