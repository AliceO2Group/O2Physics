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

#ifndef PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_MEMBERFUNCTIONS_H_
#define PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_MEMBERFUNCTIONS_H_

// a) Called directly in init(...);
// b) Called directly in process(...);
// *) Particle weights;
// *) Utility;

// a) Called directly in init(...):
// void BookBaseList()
// void DefaultConfiguration();
// void DefaultBooking();
// void DefaultBinning();
// void DefaultCuts(); // Remark: has to be called after DefaultBinning(), since
// some default cuts are defined through default binning, to ease bookeeping
// void BookAndNestAllLists()
// void BookEventHistograms()
// void BookParticleHistograms()
// void BookQvectorHistograms()
// void BookCorrelationsHistograms()
// void BookWeightsHistograms()
// void BookNestedLoopsHistograms()
// void BookTest0Histograms()
//  void StoreLabelsInPlaceholder();
//  Bool_t RetrieveCorrelationsLabels();
// void BookResultsHistograms()

// b) Called directly in process(...):
// void ResetEventByEventQuantities();
// void FillEventHistograms(aod::Collision const& collision, aod::Tracks const&
// tracks, const Int_t rs, const Int_t ba); // reco or sim, before or after
// event cuts Bool_t EventCuts(aod::Collision const& collision, aod::Tracks
// const& tracks) void FillParticleHistograms(aod::Track const& track, const
// Int_t rs, const Int_t ba); // reco or sim, before or after particle cuts
// Bool_t ParticleCuts(aod::Track const& track)
// Double_t Weight(const Double_t &value, const char *variable)
// void CalculateCorrelations();
// void CalculateNestedLoops(); // calculate all standard isotropic correlations
// with nested loops Double_t CalculateCustomNestedLoop(TArrayI *harmonics); //
// calculate nested loop for the specified harmonics

// *) Called after all events are processed (former "Terminate()"):
// void ComparisonNestedLoopsVsCorrelations();

// *) Q-vectors:
// TComplex Q(Int_t n, Int_t p)
// TComplex One(Int_t n1)
// TComplex Two(Int_t n1, Int_t n2)
// TComplex Three(Int_t n1, Int_t n2, Int_t n3)
// TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)
// ... TBI 20220809 port the rest ...
// void ResetQ(); // reset the components of generic Q-vectors

// *) Particle weights:
// void SetWeightsHist(TH1D* const hist, const char *variable)
// TH1D* GetWeightsHist(const char *variable)
// TH1D* GetHistogramWithWeights(const char* filePath, const char* runNumber,
// const char* variable)

// *) Test0:
// TObjArray* GetObjArrayWithLabels(const char *filePath)

// *) Utility:
// TObject* GetObjectFromList(TList *list, Char_t *objectName);

//============================================================

void BookBaseList()
{
  // ...

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  TList* temp = new TList();
  temp->SetOwner(kTRUE);
  fBaseList.setObject(temp);
  // fBaseList.object->SetName("4444");

  fBasePro = new TProfile("fBasePro", "flags for the whole analysis",
                          eConfiguration_N, 0., eConfiguration_N);
  fBasePro->SetStats(kFALSE);
  fBasePro->SetLineColor(eColor);
  fBasePro->SetFillColor(eFillColor);

  // Remark: If I want to change the ordering of bin lables, simply change the
  // ordering in enum eConfiguration { ... }, nothing needs to be changed here.
  fBasePro->GetXaxis()->SetBinLabel(eTaskName,
                                    Form("fTaskName = %s", fTaskName.Data()));
  fBasePro->GetXaxis()->SetBinLabel(eRunNumber,
                                    Form("fRunNumber = %s", fRunNumber.Data()));
  fBasePro->GetXaxis()->SetBinLabel(eVerbose, "fVerbose");
  fBasePro->Fill(eVerbose - 0.5, (Int_t)fVerbose);
  fBasePro->GetXaxis()->SetBinLabel(eVerboseForEachParticle,
                                    "fVerboseForEachParticle");
  fBasePro->Fill(eVerboseForEachParticle - 0.5, (Int_t)fVerboseForEachParticle);
  fBasePro->GetXaxis()->SetBinLabel(eUseCCDB, "fUseCCDB");
  fBasePro->Fill(eUseCCDB - 0.5, (Int_t)fUseCCDB);

  fBaseList->Add(fBasePro);

} // void BookBaseList()

//============================================================

void DefaultConfiguration()
{
  // Default task configuration.
  // a) Default values are hardcoded as Configurables in the file
  // MuPa-Configurables.h b) If corresponding fields are available in an
  // external json file at run time, the default values hardcoded here are
  // overwritten with values set in json file.
  //    Remember #1: To take into account configuration from external json file,
  //    use additional flag for executable, e.g.: --configuration
  //    json://my-config.json
  //                 And yes, --configuration json://my-config.json has to be
  //                 supplied to all executables in the pipe chain (workflow),
  //                 not only to mine!
  //    Remember #2: If names of Configurables in the json file are not
  //    identical to the internal definitions in MuPa-Configurables.h, the
  //    settings in json file are silently ignored.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // Configurable<string> cfTaskName{ ... }
  fTaskName = TString(cfTaskName);

  // Configurable<string> cfRunNumber{ ... } // TBI temporarily here, this shall
  // be eventually automatically obtained from e.g. collision.bc().runNumber()
  fRunNumber = TString(cfRunNumber);

  // Configurable<bool> cfVerbose{ ... }
  fVerbose = cfVerbose;

  // Configurable<bool> cfVerboseForEachParticle{ ... }
  fVerboseForEachParticle = cfVerboseForEachParticle;

  // Configurable<bool> cfUseCCDB{ ... }
  fUseCCDB = cfUseCCDB;

  // ...

  // Configurable<bool> cfCalculateTest0{ ... };
  fCalculateTest0 = cfCalculateTest0;

  // Configurable<string> cfFileWithLabels{ ... }
  fFileWithLabels = TString(cfFileWithLabels);

  // Configurable<bool> cfUsePhiWeights{"cfUsePhiWeights", false, "use or not
  // phi weights"};
  pw_a.fUseWeights[wPHI] = cfUsePhiWeights;

  // Configurable<bool> cfUsePtWeights{"cfUsePtWeights", false, "use or not pt
  // weights"};
  pw_a.fUseWeights[wPT] = cfUsePtWeights;

  // Configurable<bool> cfUseEtaWeights{"cfUseEtaWeights", false, "use or not
  // eta weights"};
  pw_a.fUseWeights[wETA] = cfUseEtaWeights;

  // Configurable<string> cfFileWithWeights{ ... }
  fFileWithWeights = TString(cfFileWithWeights);

  // ...

  // task->SetCalculateQvector(kTRUE);
  fCalculateQvector = kTRUE;

  // task->SetCalculateCorrelations(kTRUE);
  fCalculateCorrelations = kTRUE;

  // task->SetCalculateNestedLoops(kFALSE);
  fCalculateNestedLoops = kFALSE;

  // task->SetCalculateCustomNestedLoop(kFALSE); // independent e-b-e
  // cross-check with custom nested loop
  fCalculateCustomNestedLoop = kFALSE;

} // void DefaultConfiguration()

//============================================================

void DefaultBooking()
{
  // Set here which histograms are booked by default.

  // a) Event histograms;
  // b) Particle histograms;
  // c) QA;

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Event histograms:
  // Each default setting can be overuled e.g. with:
  // task->SetBookEventHistograms("NumberOfEvents",kFALSE);
  ceh_a.fBookEventHistograms[eNumberOfEvents] = kTRUE;
  ceh_a.fBookEventHistograms[eTotalMultiplicity] = kTRUE;
  ceh_a.fBookEventHistograms[eSelectedTracks] = kTRUE;
  ceh_a.fBookEventHistograms[eCentrality] = kTRUE;
  ceh_a.fBookEventHistograms[eVertex_x] = kTRUE;
  ceh_a.fBookEventHistograms[eVertex_y] = kTRUE;
  ceh_a.fBookEventHistograms[eVertex_z] = kTRUE;
  ceh_a.fBookEventHistograms[eNContributors] = kTRUE;

  // b) Particle histograms:
  // Each default setting can be overuled e.g. with:
  // task->SetBookParticleHistograms("Phi",kFALSE);
  cph_a.fBookParticleHistograms[ePhi] = kTRUE;
  cph_a.fBookParticleHistograms[ePt] = kTRUE;
  cph_a.fBookParticleHistograms[eEta] = kTRUE;
  cph_a.fBookParticleHistograms[etpcNClsCrossedRows] = kTRUE;
  cph_a.fBookParticleHistograms[eDCA_xy] = kTRUE;
  cph_a.fBookParticleHistograms[eDCA_z] = kTRUE;

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
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Default binning for event histograms:
  // task->SetEventHistogramsBins("NumberOfEvents",1,0.,1.);
  ceh_a.fEventHistogramsBins[eNumberOfEvents][0] = 1;
  ceh_a.fEventHistogramsBins[eNumberOfEvents][1] = 0.;
  ceh_a.fEventHistogramsBins[eNumberOfEvents][2] = 1.;
  // task->SetEventHistogramsBins("TotalMultiplicity",1000,0.,1000.);
  ceh_a.fEventHistogramsBins[eTotalMultiplicity][0] = 10000;
  ceh_a.fEventHistogramsBins[eTotalMultiplicity][1] = 0.;
  ceh_a.fEventHistogramsBins[eTotalMultiplicity][2] = 10000.;
  // task->SetEventHistogramsBins("SelectedTracks",1000,0.,1000.);
  ceh_a.fEventHistogramsBins[eSelectedTracks][0] = 10000;
  ceh_a.fEventHistogramsBins[eSelectedTracks][1] = 0.;
  ceh_a.fEventHistogramsBins[eSelectedTracks][2] = 10000.;
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
  // task->SetEventHistogramsBins("NContributors",1000,0.,1000.);
  ceh_a.fEventHistogramsBins[eNContributors][0] = 1000;
  ceh_a.fEventHistogramsBins[eNContributors][1] = 0.;
  ceh_a.fEventHistogramsBins[eNContributors][2] = 1000.;

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

  // ctd with  etpcNClsCrossedRows, eDCA_xy, eDCA_z

} // void DefaultBinning()

//============================================================

void DefaultCuts()
{
  // Define default cuts. Default cuts are hardwired in MuPa-Configurables.h.

  // a) Default event cuts;
  // b) Default particle cuts.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Default event cuts:
  ceh_a.fEventCuts[eNumberOfEvents][eMin] =
    cNumberOfEvents_min; // Configurable<int>
                         // cNumberOfEvents_min{"cNumberOfEvents_min", ...
  ceh_a.fEventCuts[eNumberOfEvents][eMax] =
    cNumberOfEvents_max; // Configurable<int>
                         // cNumberOfEvents_max{"cNumberOfEvents_max", ...

  ceh_a.fEventCuts[eTotalMultiplicity][eMin] = cTotalMultiplicity_min; // Configurable<int>
                                                                       // cTotalMultiplicity_min{"cTotalMultiplicity_min",
                                                                       // ...
  ceh_a.fEventCuts[eTotalMultiplicity][eMax] = cTotalMultiplicity_max; // Configurable<int>
                                                                       // cTotalMultiplicity_max{"cTotalMultiplicity_max",
                                                                       // ...

  ceh_a.fEventCuts[eSelectedTracks][eMin] =
    cSelectedTracks_min; // Configurable<int>
                         // cSelectedTracks_min{"cSelectedTracks_min", ...
  ceh_a.fEventCuts[eSelectedTracks][eMax] =
    cSelectedTracks_max; // Configurable<int>
                         // cSelectedTracks_max{"cSelectedTracks_max", ...

  ceh_a.fEventCuts[eCentrality][eMin] =
    cCentrality_min; // Configurable<int> cCentrality_min{"cCentrality_min",
                     // ...
  ceh_a.fEventCuts[eCentrality][eMax] =
    cCentrality_max; // Configurable<int> cCentrality_max{"cCentrality_max",
                     // ...

  ceh_a.fEventCuts[eVertex_x][eMin] =
    cVertex_x_min; // Configurable<int> cVertex_x_min{"cVertex_x_min", ...
  ceh_a.fEventCuts[eVertex_x][eMax] =
    cVertex_x_max; // Configurable<int> cVertex_x_max{"cVertex_x_max", ...

  ceh_a.fEventCuts[eVertex_y][eMin] =
    cVertex_y_min; // Configurable<int> cVertex_y_min{"cVertex_y_min", ...
  ceh_a.fEventCuts[eVertex_y][eMax] =
    cVertex_y_max; // Configurable<int> cVertex_y_max{"cVertex_y_max", ...

  ceh_a.fEventCuts[eVertex_z][eMin] =
    cVertex_z_min; // Configurable<int> cVertex_z_min{"cVertex_z_min", ...
  ceh_a.fEventCuts[eVertex_z][eMax] =
    cVertex_z_max; // Configurable<int> cVertex_z_max{"cVertex_z_max", ...

  ceh_a.fEventCuts[eNContributors][eMin] =
    cNContributors_min; // Configurable<int>
                        // cNContributors_min{"cNContributors_min", ...
  ceh_a.fEventCuts[eNContributors][eMax] =
    cNContributors_max; // Configurable<int>
                        // cNContributors_max{"cNContributors_max", ...

  // ...

  // b) Default particle cuts:

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
  // *) Nested loops;
  // *) Test0;
  // *) Results.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
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

  // *) Nested loops:
  fNestedLoopsList = new TList();
  fNestedLoopsList->SetName("NestedLoops");
  fNestedLoopsList->SetOwner(kTRUE);
  fBaseList->Add(fNestedLoopsList);

  // *) Test0:
  fTest0List = new TList();
  fTest0List->SetName("Test0");
  fTest0List->SetOwner(kTRUE);
  fBaseList->Add(fTest0List);

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
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fEventHistogramsPro = new TProfile("fEventHistogramsPro",
                                     "flags for event histograms", 25, 0., 25.);
  fEventHistogramsPro->SetStats(kFALSE);
  fEventHistogramsPro->SetLineColor(eColor);
  fEventHistogramsPro->SetFillColor(eFillColor);
  // ...
  fEventHistogramsList->Add(fEventHistogramsPro);

  Int_t fBeforeAfterColor[2] = {
    kRed,
    kGreen}; //! [0 = kRed,1 = kGreen] TBI 20220713 only temporarily here

  // b) Book specific control event histograms:
  TString stype[eEventHistograms_N] = {
    "NumberOfEvents", "TotalMultiplicity", "SelectedTracks",
    "Centrality", "Vertex_x", "Vertex_y",
    "Vertex_z", "NContributors"}; // keep in sync. with enum
                                  // eEventHistograms
  TString srs[2] = {"rec", "sim"};
  TString sba[2] = {"before", "after"};

  for (Int_t t = 0; t < eEventHistograms_N;
       t++) // type, see enum eEventHistograms
  {
    if (!ceh_a.fBookEventHistograms[t]) {
      continue;
    }
    for (Int_t rs = 0; rs < 2; rs++) // reco/sim
    {
      for (Int_t ba = 0; ba < 2; ba++) // before/after cuts
      {
        ceh_a.fEventHistograms[t][rs][ba] = new TH1D(
          Form("fEventHistograms[%s][%s][%s]", stype[t].Data(),
               srs[rs].Data(), sba[ba].Data()),
          Form("%s, %s, %s", stype[t].Data(), srs[rs].Data(), sba[ba].Data()),
          (Int_t)ceh_a.fEventHistogramsBins[t][0],
          ceh_a.fEventHistogramsBins[t][1], ceh_a.fEventHistogramsBins[t][2]);
        ceh_a.fEventHistograms[t][rs][ba]->SetLineColor(fBeforeAfterColor[ba]);
        ceh_a.fEventHistograms[t][rs][ba]->SetFillColor(fBeforeAfterColor[ba] -
                                                        10);
        fEventHistogramsList->Add(ceh_a.fEventHistograms[t][rs][ba]);
      } // for(Int_t ba=0;ba<2;ba++)
    }   // for(Int_t rs=0;rs<2;rs++) // reco/sim
  }     // for(Int_t t=0;t<eEventHistograms_N;t++) // type, see enum
        // eEventHistograms

} // void BookEventHistograms()

//============================================================

void BookParticleHistograms()
{
  // Book all particle histograms.

  // a) Book the profile holding flags;
  // b) Book specific particle histograms.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fParticleHistogramsPro = new TProfile(
    "fParticleHistogramsPro", "flags for particle histograms", 25, 0., 25.);
  fParticleHistogramsPro->SetStats(kFALSE);
  fParticleHistogramsPro->SetLineColor(eColor);
  fParticleHistogramsPro->SetFillColor(eFillColor);
  // ...
  fParticleHistogramsList->Add(fParticleHistogramsPro);

  Int_t fBeforeAfterColor[2] = {
    kRed,
    kGreen}; //! [0 = kRed,1 = kGreen] TBI 20220713 only temporarily here

  // b) Book specific control particle histograms:
  TString stype[eParticleHistograms_N] = {
    "Phi", "Pt", "Eta", "tpcNClsCrossedRows",
    "DCA_xy", "DCA_z"}; // keep in sync. with enum eParticleHistograms
  TString srs[2] = {"rec", "sim"};
  TString sba[2] = {"before", "after"};

  for (Int_t t = 0; t < eParticleHistograms_N;
       t++) // type, see enum eParticleHistograms
  {
    if (!cph_a.fBookParticleHistograms[t]) {
      continue;
    }
    for (Int_t rs = 0; rs < 2; rs++) // reco/sim
    {
      for (Int_t ba = 0; ba < 2; ba++) // before/after cuts
      {
        cph_a.fParticleHistograms[t][rs][ba] = new TH1D(
          Form("fParticleHistograms[%s][%s][%s]", stype[t].Data(),
               srs[rs].Data(), sba[ba].Data()),
          Form("%s, %s, %s", stype[t].Data(), srs[rs].Data(), sba[ba].Data()),
          (Int_t)cph_a.fParticleHistogramsBins[t][0],
          cph_a.fParticleHistogramsBins[t][1],
          cph_a.fParticleHistogramsBins[t][2]);
        cph_a.fParticleHistograms[t][rs][ba]->SetLineColor(
          fBeforeAfterColor[ba]);
        cph_a.fParticleHistograms[t][rs][ba]->SetFillColor(
          fBeforeAfterColor[ba] - 10);
        fParticleHistogramsList->Add(cph_a.fParticleHistograms[t][rs][ba]);
      } // for(Int_t ba=0;ba<2;ba++)
    }   // for(Int_t rs=0;rs<2;rs++) // reco/sim
  }     // for(Int_t t=0;t<eParticleHistograms_N;t++) // type, see enum
        // eParticleHistograms

} // void BookParticleHistograms()

//============================================================

void BookQvectorHistograms()
{
  // Book all Q-vector histograms.

  // a) Book the profile holding flags;
  // b) ...

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fQvectorFlagsPro =
    new TProfile("fQvectorFlagsPro", "flags for Q-vector objects", 3, 0., 3.);
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
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fCorrelationsFlagsPro = new TProfile("fCorrelationsFlagsPro",
                                       "flags for correlations", 3, 0., 3.);
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
  TString oVariable[4] = {
    "#varphi_{1}-#varphi_{2}",
    "#varphi_{1}+#varphi_{2}-#varphi_{3}-#varphi_{4}",
    "#varphi_{1}+#varphi_{2}+#varphi_{3}-#varphi_{4}-#varphi_{5}-#varphi_{6}",
    "#varphi_{1}+#varphi_{2}+#varphi_{3}+#varphi_{4}-#varphi_{5}-#varphi_{6}-"
    "#varphi_{7}-#varphi_{8}"};

  TString vvVariable[3] = {"int", "mult", "cent"};

  // c) Histograms:
  for (Int_t k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]
  {
    for (Int_t n = 0; n < 6; n++) // harmonic [n=1,n=2,...,n=6]
    {
      for (Int_t v = 0; v < 3;
           v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
      {
        // ... TBI 20220809 ... port the rest

        // c_a.fCorrelationsPro[k][n][v] = new
        // TProfile(Form("c_a.fCorrelationsPro[%d][%d][%s]",k,n,vvVariable[v].Data()),harmonicArray.Data(),vvvariableNBins[v],vvvariableMinMax[v][0],vvvariableMinMax[v][1]);
        c_a.fCorrelationsPro[k][n][v] = new TProfile(
          Form("fCorrelationsPro[%d][%d][%s]", k, n, vvVariable[v].Data()),
          "some title", 2000, 0., 2000.);
        c_a.fCorrelationsPro[k][n][v]->SetStats(kFALSE);
        c_a.fCorrelationsPro[k][n][v]->Sumw2();
        c_a.fCorrelationsPro[k][n][v]->GetXaxis()->SetTitle(
          vvVariable[v].Data());
        c_a.fCorrelationsPro[k][n][v]->GetYaxis()->SetTitle(
          Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1),
               oVariable[k].Data()));
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
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fWeightsFlagsPro =
    new TProfile("fWeightsFlagsPro", "flags for particle weights", 3, 0., 3.);
  fWeightsFlagsPro->SetStats(kFALSE);
  fWeightsFlagsPro->SetLineColor(eColor);
  fWeightsFlagsPro->SetFillColor(eFillColor);
  fWeightsFlagsPro->GetXaxis()->SetLabelSize(0.05);
  fWeightsFlagsPro->GetXaxis()->SetBinLabel(1, "w_{#varphi}");
  fWeightsFlagsPro->GetXaxis()->SetBinLabel(2, "w_{p_{t}}");
  fWeightsFlagsPro->GetXaxis()->SetBinLabel(3, "w_{#eta}");
  for (Int_t w = 0; w < eWeights_N; w++) // use weights [phi,pt,eta]
  {
    if (pw_a.fUseWeights[w])
      fWeightsFlagsPro->Fill(w + 0.5, 1.);
  }
  fWeightsList->Add(fWeightsFlagsPro);

  // b) Common local labels: TBI 20220713 book before
  TString sVariable[eWeights_N] = {"#varphi", "p_{t}", "#eta"}; // [phi,pt,eta]
  TString sWeights[eWeights_N] = {"w_{#varphi}", "w_{p_{t}}", "w_{#eta}"};

  // c) Histograms:
  for (Int_t w = 0; w < eWeights_N; w++) // use weights [phi,pt,eta]
  {
    if (!pw_a.fUseWeights[w]) {
      continue;
    }
    if (!pw_a.fWeightsHist[w]) {
      // yes, because these histos are cloned from the
      // external ones, see SetWeightsHist(TH1D* const
      // hist, const char *variable)

      // pw_a.fWeightsHist[w] = new
      // TH1D(Form("fWeightsHist[%d]",w),"",(Int_t)fKinematicsBins[w][0],fKinematicsBins[w][1],fKinematicsBins[w][2]);
      pw_a.fWeightsHist[w] =
        new TH1D(Form("fWeightsHist[%d]", w), "", 200, -100., 100.);
      pw_a.fWeightsHist[w]->SetTitle(
        Form("Particle weights for %s", sWeights[w].Data()));
      pw_a.fWeightsHist[w]->SetStats(kFALSE);
      pw_a.fWeightsHist[w]->GetXaxis()->SetTitle(sVariable[w].Data());
      pw_a.fWeightsHist[w]->SetFillColor(eFillColor);
      pw_a.fWeightsHist[w]->SetLineColor(eColor);
    }
    fWeightsList->Add(pw_a.fWeightsHist[w]);
  } // for(Int_t w=0;w<eWeights_N;w++) // use weights [phi,pt,eta]

} // void BookWeightsHistograms()

//============================================================

void BookNestedLoopsHistograms()
{
  // Book all nested loops histograms.

  // a) Book the profile holding flags;
  // *) ...

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fNestedLoopsFlagsPro =
    new TProfile("fNestedLoopsFlagsPro", "flags for nested loops", 2, 0., 2.);
  fNestedLoopsFlagsPro->SetStats(kFALSE);
  fNestedLoopsFlagsPro->GetXaxis()->SetLabelSize(0.05);
  fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateNestedLoops");
  fNestedLoopsFlagsPro->Fill(0.5, fCalculateNestedLoops);
  fNestedLoopsFlagsPro->Fill(1.5, fCalculateCustomNestedLoop);
  fNestedLoopsList->Add(fNestedLoopsFlagsPro);

  if (!(fCalculateNestedLoops || fCalculateCustomNestedLoop)) {
    return;
  }

  const Int_t iMaxSize = 2e4;
  nl_a.ftaNestedLoops[0] =
    new TArrayD(iMaxSize); // ebe container for azimuthal angles
  nl_a.ftaNestedLoops[1] = new TArrayD(
    iMaxSize); // ebe container for particle weights (product of all)

  // TBI 20220823 port here if(fCalculatePtCorrelations) { ... } and
  // if(fCalculateEtaCorrelations) { ... }

  if (!fCalculateNestedLoops) {
    return;
  }

  // b) Common local labels (keep 'em in sync with BookCorrelationsHistograms())
  TString oVariable[4] = {
    "#varphi_{1}-#varphi_{2}",
    "#varphi_{1}+#varphi_{2}-#varphi_{3}-#varphi_{4}",
    "#varphi_{1}+#varphi_{2}+#varphi_{3}-#varphi_{4}-#varphi_{5}-#varphi_{6}",
    "#varphi_{1}+#varphi_{2}+#varphi_{3}+#varphi_{4}-#varphi_{5}-#varphi_{6}-"
    "#varphi_{7}-#varphi_{8}"};

  /*
  Int_t vvvariableNBins[5] =
  {1,(Int_t)fMultiplicityBins[0],(Int_t)fCentralityBins[0],
                              fUseCustomKineDependenceBins[PTq] ?
  fKineDependenceBins[PTq]->GetSize()-1 : (Int_t)fKinematicsBins[PT][0],
                              fUseCustomKineDependenceBins[ETAq] ?
  fKineDependenceBins[ETAq]->GetSize()-1 : (Int_t)fKinematicsBins[ETA][0]};
  Double_t vvvariableMinMax[5][2] = { {0.,1.}, // integrated
                                      {fMultiplicityBins[1],fMultiplicityBins[2]},
  // multiplicity {fCentralityBins[1],fCentralityBins[2]}, // centrality
                                      {fKinematicsBins[PT][1],fKinematicsBins[PT][2]},
                                      {fKinematicsBins[ETA][1],fKinematicsBins[ETA][2]}
                                    };
  */

  TString vvVariable[3] = {"integrated", "multiplicity", "centrality"};

  for (Int_t k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]
  {
    for (Int_t n = 0; n < 6; n++) // harmonic [n=1,n=2,...,n=6]
    {
      for (Int_t v = 0; v < 3;
           v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality]
      {

        // if(PTKINE == v  && !fCalculatePtCorrelations){continue;}
        // if(ETAKINE == v  && !fCalculateEtaCorrelations){continue;}

        /*
        // per demand, custom meeting for kine dependence:
        if(PTKINE == v  && fUseCustomKineDependenceBins[PTq])
        {
         fNestedLoopsPro[k][n][v] = new
        TProfile(Form("fNestedLoopsPro[%d][%d][%d]",k,n,v),Form("#LT#LTcos[%s(%s)]#GT#GT",1==n+1?"":Form("%d",n+1),oVariable[k].Data()),fKineDependenceBins[PTq]->GetSize()-1,fKineDependenceBins[PTq]->GetArray());
        }
        else if(ETAKINE == v  && fUseCustomKineDependenceBins[ETAq])
        {
         fNestedLoopsPro[k][n][v] = new
        TProfile(Form("fNestedLoopsPro[%d][%d][%d]",k,n,v),Form("#LT#LTcos[%s(%s)]#GT#GT",1==n+1?"":Form("%d",n+1),oVariable[k].Data()),fKineDependenceBins[ETAq]->GetSize()-1,fKineDependenceBins[ETAq]->GetArray());
        }
        else
        {
        */
        // the default binning:
        // fNestedLoopsPro[k][n][v] = new
        // TProfile(Form("fNestedLoopsPro[%d][%d][%d]",k,n,v),Form("#LT#LTcos[%s(%s)]#GT#GT",1==n+1?"":Form("%d",n+1),oVariable[k].Data()),vvvariableNBins[v],vvvariableMinMax[v][0],vvvariableMinMax[v][1]);
        nl_a.fNestedLoopsPro[k][n][v] = new TProfile(
          Form("fNestedLoopsPro[%d][%d][%d]", k, n, v),
          Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1),
               oVariable[k].Data()),
          2000, 0., 2000.);
        // } // else
        nl_a.fNestedLoopsPro[k][n][v]->SetStats(kFALSE);
        nl_a.fNestedLoopsPro[k][n][v]->Sumw2();
        nl_a.fNestedLoopsPro[k][n][v]->GetXaxis()->SetTitle(
          vvVariable[v].Data());
        // fNestedLoopsPro[k][n][v]->SetFillColor(colorsW[v]-10);
        // fNestedLoopsPro[k][n][v]->SetLineColor(colorsW[v]);
        /*
        if(fUseFixedNumberOfRandomlySelectedTracks && 1==v) // just a warning
        for the meaning of multiplicity in this special case
        {
         nl_a.fNestedLoopsPro[k][n][1]->GetXaxis()->SetTitle("WARNING: for each
        multiplicity, fFixedNumberOfRandomlySelectedTracks is selected randomly
        in Q-vector");
        }
        */

        fNestedLoopsList->Add(nl_a.fNestedLoopsPro[k][n][v]);
      } // for(Int_t v=0;v<5;v++) // variable [0=integrated,1=vs.
        // multiplicity,2=vs. centrality]
    }   // for(Int_t n=0;n<6;n++) // harmonic [n=1,n=2,...,n=6]
  }     // for(Int_t n=0;n<6;n++) // harmonics [n=1,n=2,...,n=6]

} // void BookNestedLoopsHistograms()

//============================================================

void BookTest0Histograms()
{
  // Book all Test0 histograms.

  // a) Book the profile holding flags;
  // b) Book placeholder and make sure all labels are stored in the placeholder;
  // c) Retreive labels from placeholder;
  // d) Book what needs to be booked.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fTest0FlagsPro = new TProfile("fTest0FlagsPro", "flags for Test0", 1, 0., 1.);
  fTest0FlagsPro->SetStats(kFALSE);
  fTest0FlagsPro->GetXaxis()->SetLabelSize(0.04);
  fTest0FlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateTest0");
  fTest0FlagsPro->Fill(0.5, fCalculateTest0);
  fTest0List->Add(fTest0FlagsPro);

  if (!fCalculateTest0) {
    return;
  }

  // b) Book placeholder and make sure all labels are stored in the placeholder:
  StoreLabelsInPlaceholder();
  if (fTest0LabelsPlaceholder) {
    fTest0List->Add(fTest0LabelsPlaceholder);
  }

  // c) Retreive labels from placeholder:
  if (!(this->RetrieveCorrelationsLabels())) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

} // void BookTest0Histograms()

//============================================================

void BookResultsHistograms()
{
  // Book all results histograms.

  // a) Book the profile holding flags;
  // *) ...

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  fResultsFlagsPro = new TProfile("fResultsFlagsPro",
                                  "flags for results histograms", 1, 0., 1.);
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

  // a) Event-by-event quantities;
  // b) Q-vectors;
  // c) Reset ebe containers for nested loops;
  // d) Fisher-Yates algorithm.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Event-by-event quantities:
  fSelectedTracks = 0;
  fCentrality = 0;

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

  // d) Reset ebe containers for nested loops:
  if (fCalculateNestedLoops || fCalculateCustomNestedLoop) {
    if (nl_a.ftaNestedLoops[0]) {
      nl_a.ftaNestedLoops[0]->Reset();
    }
    if (nl_a.ftaNestedLoops[1]) {
      nl_a.ftaNestedLoops[1]->Reset();
    }

    // TBI 20220803 port still if(fCalculatePtCorrelations){...} and
    // if(fCalculateEtaCorrelations){...}

  } // if(fCalculateNestedLoops||fCalculateCustomNestedLoop)

  // ... TBI 20220809 port the rest ...

} // void ResetEventByEventQuantities()

//============================================================

// Bool_t EventCuts(aod::Collision const& collision, soa::Join<aod::Tracks,
// aod::TracksExtra, aod::TracksDCA> const& tracks)
Bool_t EventCuts(aod::Collision const& collision, aod::Tracks const& tracks)
{
  // ...

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // NumberOfEvents: => cut directly in void process( ... )

  // TotalMultiplicity:
  if ((tracks.size() < ceh_a.fEventCuts[eTotalMultiplicity][eMin]) ||
      (tracks.size() > ceh_a.fEventCuts[eTotalMultiplicity][eMax])) {
    return kFALSE;
  }

  // SelectedTracks: => cut directly in void process( ... )

  // Centrality: TBI
  //  if ((  TBI   < ceh_a.fEventCuts[eCentrality][eMin]) || (  TBI   >
  //  ceh_a.fEventCuts[eCentrality][eMax])) {
  //    return kFALSE;
  //  }

  // Vertex_x:
  if ((collision.posX() < ceh_a.fEventCuts[eVertex_x][eMin]) ||
      (collision.posX() > ceh_a.fEventCuts[eVertex_x][eMax])) {
    return kFALSE;
  }

  // Vertex_y:
  if ((collision.posY() < ceh_a.fEventCuts[eVertex_y][eMin]) ||
      (collision.posY() > ceh_a.fEventCuts[eVertex_y][eMax])) {
    return kFALSE;
  }

  // Vertex_z:
  if ((collision.posZ() < ceh_a.fEventCuts[eVertex_z][eMin]) ||
      (collision.posZ() > ceh_a.fEventCuts[eVertex_z][eMax])) {
    return kFALSE;
  }

  // NContributors:
  if ((collision.numContrib() < ceh_a.fEventCuts[eNContributors][eMin]) ||
      (collision.numContrib() > ceh_a.fEventCuts[eNContributors][eMax])) {
    return kFALSE;
  }

  return kTRUE;

} // void EventCuts(aod::Collision const& collision, aod::Tracks const& tracks)

//============================================================

// void FillEventHistograms(aod::Collision const& collision,
// soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& tracks, const
// Int_t rs, const Int_t ba) // reco or sim, before or after event cuts
void FillEventHistograms(
  aod::Collision const& collision, aod::Tracks const& tracks, const Int_t rs,
  const Int_t ba) // reco or sim, before or after event cuts
{
  // Fill all event histograms.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  ceh_a.fEventHistograms[eNumberOfEvents][rs][ba]->Fill(0.5);
  ceh_a.fEventHistograms[eTotalMultiplicity][rs][ba]->Fill(tracks.size());
  // SelectedTracks => filled directly in void process( ... )
  // ceh_a.fEventHistograms[eCentrality][rs][ba]->Fill(  TBI  );  => TBI
  // 20231007 not ready yet
  ceh_a.fEventHistograms[eVertex_x][rs][ba]->Fill(collision.posX());
  ceh_a.fEventHistograms[eVertex_y][rs][ba]->Fill(collision.posY());
  ceh_a.fEventHistograms[eVertex_z][rs][ba]->Fill(collision.posZ());
  ceh_a.fEventHistograms[eNContributors][rs][ba]->Fill(collision.numContrib());

} // void FillEventHistograms(aod::Collision const& collision, aod::Tracks
// const& tracks, const Int_t rs, const Int_t ba); // reco or sim, before or
// after event cuts

//============================================================

Bool_t ParticleCuts(aod::Track const& track)
// Bool_t ParticleCuts(auto const& track)
{
  // ...

  if (fVerboseForEachParticle) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  if ((track.pt() < pt_min) || (track.pt() > pt_max)) {
    return kFALSE;
  }

  return kTRUE;

} // void ParticleCuts(aod::Track const& tracks)

//============================================================

// void FillParticleHistograms(auto const& track, const Int_t rs, const Int_t
// ba) // reco or sim, before or after particle cuts
void FillParticleHistograms(
  aod::Track const& track, const Int_t rs,
  const Int_t ba) // reco or sim, before or after particle cuts
{
  // Fill all particle histograms.

  if (fVerboseForEachParticle) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // From aod::Tracks
  cph_a.fParticleHistograms[ePhi][rs][ba]->Fill(track.phi());
  cph_a.fParticleHistograms[ePt][rs][ba]->Fill(track.pt());
  cph_a.fParticleHistograms[eEta][rs][ba]->Fill(track.eta());

  /*
  // From aod::TracksExtra
  cph_a.fParticleHistograms[etpcNClsCrossedRows][rs][ba]->Fill(track.tpcNClsCrossedRows());

  // From aod::TracksDCA
  cph_a.fParticleHistograms[eDCA_xy][rs][ba]->Fill(track.dcaXY());
  cph_a.fParticleHistograms[eDCA_z][rs][ba]->Fill(track.dcaZ());
  */

} // void FillParticleHistograms(aod::Track const& track, const Int_t rs, const
// Int_t ba); // reco or sim, before or after particle cuts

//============================================================

void CalculateCorrelations()
{
  // Calculate analytically multiparticle correlations from Q-vectors
  // In this method, only isotropic correlations for which all harmonics are the
  // same are evaluated. For the calculus of generic multiparticle correlations,
  // see method CalculateGenericCorrelations()

  // a) Flush 'n' fill the generic Q-vectors;
  // b) Calculate correlations;
  // c) Flush the generic Q-vectors.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
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
    // if(fSelectedTracks<2){return;}
    // if(fVerbose){cout<<Form("   => CalculateCorrelations(void), 2p, h = %d
    // .... ",h)<<endl;}
    TComplex two = Two(h, -h);
    Double_t twoC = two.Re(); // cos
    // Double_t twoS = two.Im(); // sin
    Double_t wTwo =
      Two(0, 0).Re(); // Weight is 'number of combinations' by default TBI
                      // 20220809 add support for other weights
    if (wTwo > 0.0) {
      twoC /= wTwo;
    } else {
      return;
    } // ... TBI 20220809 ... use the line below eventually
    // else { Exit(__PRETTY_FUNCTION__,__LINE__,Form("wTwo = %f is not positive.
    // fSelectedTracks = %d",wTwo,fSelectedTracks)); return; }

    /* ... TBI 20220809 ... enable eventually
  if(fCalculateCustomNestedLoop)
  {
   // e-b-e sanity check:
   TArrayI *harmonics = new TArrayI(2);
   harmonics->SetAt(h,0);
   harmonics->SetAt(-h,1);
   Double_t nestedLoopValue = this->CalculateCustomNestedLoop(harmonics);
   if(TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(twoC -
  nestedLoopValue)>1.e-5)
   {
    Exit(__PRETTY_FUNCTION__,__LINE__,Form("nestedLoopValue = %f is not the same
  as twoC = %f",nestedLoopValue,twoC)); exit(1);
   }
   else
   {
    cout<<Form("=> e-b-e check with CustomNestedLoop is OK for isotropic 2-p,
  harmonic %d",h)<<endl;
    //cout<<Form("   value = %f",twoC)<<endl;
   }
   delete harmonics; harmonics = NULL;
  } // if(fCalculateCustomNestedLoop)

  // for on-the-fly and internal validation, rescale results with theoretical
  value: if(fCalculateOnTheFly && fOnTheFlyFlowAmplitudes &&
  fRescaleWithTheoreticalInput &&
     TMath::Abs(fOnTheFlyFlowAmplitudes->GetAt(h-1))>0.){twoC/=pow(fOnTheFlyFlowAmplitudes->GetAt(h-1),2.);}
  else if(fUseInternalValidation && fInternalValidationAmplitudes &&
  fRescaleWithTheoreticalInput &&
          TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){twoC/=pow(fInternalValidationAmplitudes->GetAt(h-1),2.);}

*/

    // integrated:
    if (c_a.fCorrelationsPro[0][h - 1][0]) {
      c_a.fCorrelationsPro[0][h - 1][0]->Fill(0.5, twoC, wTwo);
    }
    // vs. multiplicity:
    // if(c_a.fCorrelationsPro[0][h-1][1]){c_a.fCorrelationsPro[0][h-1][1]->Fill(fSelectedTracks+0.5,twoC,wTwo);}
    // vs. centrality:
    // if(c_a.fCorrelationsPro[0][h-1][2]){c_a.fCorrelationsPro[0][h-1][2]->Fill(fCentrality,twoC,wTwo);}

    // ... TBI 20220809 port the rest ...

  } // for(Int_t h=1;h<=gMaxHarmonic;h++) // harmonic

  // c) Flush the generic Q-vectors:
  ResetQ();

} // void CalculateCorrelations()

//============================================================

void CalculateNestedLoops()
{
  // Calculate correlations with nested loops.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  LOGF(info, "\033[1;32m fSelectedTracks = %d\033[0m", fSelectedTracks);
  Int_t nParticles = fSelectedTracks;

  /* TBI 20220823 enable the lines below eventually
  if(fUseFixedNumberOfRandomlySelectedTracks)
  {
   nParticles = 0;
   for(Int_t i=0;i<ftaNestedLoops[0]->GetSize();i++)
   {
    if(TMath::Abs(ftaNestedLoops[0]->GetAt(i)) > 0. &&
  TMath::Abs(ftaNestedLoops[1]->GetAt(i)) > 0.){nParticles++;}
   }
  }
   cout<<"nParticles = "<<nParticles<<endl;
  */

  // 2p:
  if (nParticles < 2) {
    return;
  }
  LOGF(info, "\033[1;32m       CalculateNestedLoops(void), 2-p correlations .... \033[0m");
  for (int i1 = 0; i1 < nParticles; i1++) {
    Double_t dPhi1 = nl_a.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl_a.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      Double_t dPhi2 = nl_a.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl_a.ftaNestedLoops[1]->GetAt(i2);
      for (int h = 0; h < 6; h++) // TBI 20220823 hardcoded "h<6"
      {
        // fill cos, 2p, integreated:
        if (nl_a.fNestedLoopsPro[0][h][0]) {
          nl_a.fNestedLoopsPro[0][h][0]->Fill(
            0.5, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
        // fill cos, 2p, vs. M:
        if (nl_a.fNestedLoopsPro[0][h][1]) {
          nl_a.fNestedLoopsPro[0][h][1]->Fill(
            fSelectedTracks + 0.5, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)),
            dW1 * dW2);
        }
        // fill cos, 2p, vs. centrality:
        if (nl_a.fNestedLoopsPro[0][h][2]) {
          nl_a.fNestedLoopsPro[0][h][2]->Fill(
            fCentrality, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
      } // for(int h=1; h<=6; h++)
    }   // for(int i2=0; i2<nParticles; i2++)
  }     // for(int i1=0; i1<nParticles; i1++)

  // TBI port the rest, i.e. 4p, 6p, 8p, etc.

} // void CalculateNestedLoops()

//============================================================

void ComparisonNestedLoopsVsCorrelations()
{
  // Make a ratio fNestedLoopsPro[....]/fCorrelationsPro[....]. If results are
  // the same, these ratios must be 1.

  // a) Integrated comparison;
  // b) Comparison vs. multiplicity;
  // c) Comparison vs. centrality;

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  Int_t nBinsQV = -44;
  Int_t nBinsNL = -44;
  Double_t valueQV = 0.;
  Double_t valueNL = 0.;

  // a) Integrated comparison:
  nBinsQV = c_a.fCorrelationsPro[0][0][0]->GetNbinsX();
  nBinsNL = nl_a.fNestedLoopsPro[0][0][0]->GetNbinsX();
  if (nBinsQV != nBinsNL) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }
  cout << endl;
  cout << "   [0] : integrated" << endl;
  for (Int_t o = 0; o < 4; o++) {
    cout << Form("   ==== <<%d>>-particle correlations ====", 2 * (o + 1))
         << endl;
    for (Int_t h = 0; h < 6; h++) {
      for (Int_t b = 1; b <= nBinsQV; b++) {
        if (c_a.fCorrelationsPro[o][h][0]) {
          valueQV = c_a.fCorrelationsPro[o][h][0]->GetBinContent(b);
        }
        if (nl_a.fNestedLoopsPro[o][h][0]) {
          valueNL = nl_a.fNestedLoopsPro[o][h][0]->GetBinContent(b);
        }
        if (TMath::Abs(valueQV) > 0. && TMath::Abs(valueNL) > 0.) {
          cout << Form("   h=%d, Q-vectors:    ", h + 1) << valueQV << endl;
          cout << Form("   h=%d, Nested loops: ", h + 1) << valueNL << endl;
          if (TMath::Abs(valueQV - valueNL) > 1.e-5) {
            LOGF(info, "\n\033[1;33m[%d][%d][%d] \033[0m\n", o, h, 0);
            LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
                 __PRETTY_FUNCTION__, __LINE__);
          }
        } // if(TMath::Abs(valueQV)>0. && TMath::Abs(valueNL)>0.)
      }   // for(Int_t b=1;b<=nBinsQV;b++)
    }     // for(Int_t h=0;h<6;h++)
    cout << endl;
  } // for(Int_t o=0;o<4;o++)

  // TBI 20230530 port the rest, i.e. vs. multiplicity, centrality, etc.

} // void ComparisonNestedLoopsVsCorrelations()

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

TComplex Three(Int_t n1, Int_t n2, Int_t n3)
{
  // Generic three-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3)]>.

  TComplex three = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) - Q(n1 + n2, 2) * Q(n3, 1) -
                   Q(n2, 1) * Q(n1 + n3, 2) - Q(n1, 1) * Q(n2 + n3, 2) +
                   2. * Q(n1 + n2 + n3, 3);

  return three;

} // TComplex Three(Int_t n1, Int_t n2, Int_t n3)

//============================================================

TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)
{
  // Generic four-particle correlation
  // <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4)]>.

  TComplex four =
    Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) -
    Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) -
    Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) -
    Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) -
    Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) + Q(n2 + n3, 2) * Q(n1 + n4, 2) -
    Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) + Q(n1 + n3, 2) * Q(n2 + n4, 2) +
    2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) +
    Q(n1 + n2, 2) * Q(n3 + n4, 2) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) +
    2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) - 6. * Q(n1 + n2 + n3 + n4, 4);

  return four;

} // TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)

//============================================================

void ResetQ()
{
  // Reset the components of generic Q-vectors. Use it whenever you call the
  // standard functions for correlations, for some custom Q-vectors.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
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
  // Copy histogram holding weights from an external file to the corresponding
  // data member.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // Basic protection:
  if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") ||
        TString(variable).EqualTo("eta"))) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
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
  pw_a.fWeightsHist[ppe] = reinterpret_cast<TH1D*>(hist->Clone());
  if (!pw_a.fWeightsHist[ppe]) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  // Flag:
  pw_a.fUseWeights[ppe] = kTRUE;

} // void SetWeightsHist(TH1D* const hist, const char *variable)

//============================================================

TH1D* GetWeightsHist(const char* variable)
{
  // The standard getter.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // Basic protection:
  if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") ||
        TString(variable).EqualTo("eta"))) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
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
  return pw_a.fWeightsHist[ppe];

} // TH1D* GetWeightsHist(const char *variable)

//============================================================

TH1D* GetHistogramWithWeights(const char* filePath, const char* runNumber,
                              const char* variable)
{
  // ...

  // a) Return value;
  // b) Basic protection for arguments;
  // c) Determine from filePath if the file in on a local machine, or in AliEn,
  // or in CCDB; d) Handle the AliEn case; e) Handle the CCDB case; f) Handle
  // the local case; g) The final touch on histogram with weights.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
    LOGF(info, "\033[1;33m filePath = %s\033[0m", filePath);
    LOGF(info, "\033[1;33m runNumber = %s\033[0m", runNumber);
    LOGF(info, "\033[1;33m variable = %s\033[0m", variable);
    LOGF(info, "\033[1;33m fTaskName = %s\033[0m", fTaskName.Data());
  }

  // a) Return value:
  TH1D* hist = NULL;

  // b) Basic protection for arguments:
  if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") ||
        TString(variable).EqualTo("eta"))) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  // c) Determine from filePath if the file in on a local machine, or in home
  // dir AliEn, or in CCDB:
  //    Algorithm: If filePath begins with "/alice/cern.ch/" then it's in home
  //    dir AliEn.
  //               If filePath begins with "/alice-ccdb.cern.ch/" then it's in
  //               CCDB. Therefore, files in AliEn and CCDB must be specified
  //               with abs path, for local files both abs and relative paths
  //               are just fine.
  Bool_t bFileIsInAliEn = kFALSE;
  Bool_t bFileIsInCCDB = kFALSE;
  if (TString(filePath).BeginsWith("/alice/cern.ch/")) {
    bFileIsInAliEn = kTRUE;
  } else {
    if (TString(filePath).BeginsWith("/alice-ccdb.cern.ch/")) {
      bFileIsInCCDB = kTRUE;
    } // else {
  }   // if (TString(filePath).BeginsWith("/alice/cern.ch/")) {

  if (bFileIsInAliEn) {
    // d) Handle the AliEn case:
    TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", "");
    if (!alien) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    TFile* weightsFile = TFile::Open(Form("alien://%s", filePath), "READ");
    if (!weightsFile) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    TList* baseList = NULL;
    weightsFile->GetObject(
      "ccdb_object", baseList); // TBI 20231008 for simplicity, harwired name
                                // of base TList is "ccdb_object" also for
                                // AliEn case, see if I need to change this
    if (!baseList) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    TList* listWithRuns =
      reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    hist = reinterpret_cast<TH1D*>(GetObjectFromList(
      listWithRuns, Form("%s_%s", variable, fTaskName.Data())));
    if (!hist) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

  } else if (bFileIsInCCDB) {

    // e) Handle the CCDB case: Remember that here I do not access the file,
    //    instead directly object in that file.
    //    My home dir in CCDB: https://alice-ccdb.cern.ch/browse/Users/a/abilandz/
    //    Inspired by:
    //    1. Discussion at:
    //    https://alice-talk.web.cern.ch/t/access-to-lhc-filling-scheme/1073/17
    //    2. See also:
    //    https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/efficiencyGlobal.cxx
    //    https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/efficiencyPerRun.cxx
    //    3. O2 Analysis Tutorial 2.0:
    //    https://indico.cern.ch/event/1267433/timetable/#20230417.detailed

    ccdb->setURL("http://alice-ccdb.cern.ch");
    if (fVerbose) {
      LOGF(info, "\033[1;32mAccessing in CCDB %s\033[0m",
           TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data());
    }

    TList* baseList =
      reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath)
                                                  .ReplaceAll("/alice-ccdb.cern.ch/", "")
                                                  .Data()));

    if (!baseList) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    TList* listWithRuns =
      reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    hist = reinterpret_cast<TH1D*>(GetObjectFromList(
      listWithRuns, Form("%s_%s", variable, fTaskName.Data())));

    if (!hist) {
      LOGF(info, "\033[1;33m%s_%s \033[0m", variable, fTaskName.Data());
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

  } else {

    // f) Handle the local case:
    //    TBI 20231008 In principle, also for the local case in O2, I could
    //    maintain the same local structure of weights as it was in AliPhysics.
    //                 But for simplicity, in O2 I organize local weights in the
    //                 same way as in AliEn or CCDB.

    // Check if the external ROOT file exists at specified path:
    if (gSystem->AccessPathName(filePath, kFileExists)) {
      LOGF(info,
           "\033[1;33m if(gSystem->AccessPathName(filePath,kFileExists)), "
           "filePath = %s \033[0m",
           filePath);
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    TFile* weightsFile = TFile::Open(filePath, "READ");
    if (!weightsFile) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    TList* baseList = NULL;
    weightsFile->GetObject(
      "ccdb_object", baseList); // TBI 20231008 for simplicity, harwired name
                                // of base TList is "ccdb_object" also for
                                // local case, see if I need to change this
    if (!baseList) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    TList* listWithRuns =
      reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    hist = reinterpret_cast<TH1D*>(GetObjectFromList(
      listWithRuns, Form("%s_%s", variable, fTaskName.Data())));
    if (!hist) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

  } // else {

  // g) The final touch on histogram with weights:
  hist->SetDirectory(0);
  hist->SetTitle(Form("%s, %s", filePath, runNumber));

  return hist;

} // TH1D* GetHistogramWithWeights(const char* filePath, const char* runNumber,
// const char* variable)

//============================================================

TObjArray* GetObjArrayWithLabels(const char* filePath)
{
  // This function extracts from an external file TObjArray named "labels", and
  // returns it. External file can be: 1) on a local computer; 2) in home
  // directory AliEn => configurable "cfFileWithLabels" must begin with
  // "/alice/cern.ch/" 3) in CCDB => configurable "cfFileWithLabels" must begin
  // with "/alice-ccdb.cern.ch/" For all CCDB wisdom, see toggle "CCDB" in page
  // "O2"

  // a) Return value;
  // b) Determine from filePath if the file in on a local machine, or in AliEn;
  // c) Handle the AliEn case;
  // d) Handle the CCDB case;
  // e) Handle the local case.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Return value:
  TObjArray* oa;

  // b) Determine from filePath if the file in on a local machine, or in home
  // dir AliEn, or in CCDB:
  //    Algorithm: If filePath begins with "/alice/cern.ch/" then it's in home
  //    dir AliEn.
  //               If filePath begins with "/alice-ccdb.cern.ch/" then it's in
  //               CCDB. Therefore, files in AliEn and CCDB must be specified
  //               with abs path, for local files both abs and relative paths
  //               are just fine.
  Bool_t bFileIsInAliEn = kFALSE;
  Bool_t bFileIsInCCDB = kFALSE;
  if (TString(filePath).BeginsWith("/alice/cern.ch/")) {
    bFileIsInAliEn = kTRUE;
  } else {
    if (TString(filePath).BeginsWith("/alice-ccdb.cern.ch/")) {
      bFileIsInCCDB = kTRUE;
    } // else {
  }   // if (TString(filePath).BeginsWith("/alice/cern.ch/")) {

  TFile* oaFile = NULL; // file holding TObjArray with all labels
  if (bFileIsInAliEn) {
    // c) Handle the AliEn case:
    TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", "");
    if (!alien) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    oaFile = TFile::Open(Form("alien://%s", filePath), "READ");
    if (!oaFile) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    oaFile->GetObject(
      "labels", oa); // TBI 20230530 hardcoded name of TObjArray is "labels",
                     // perhaps I can do this also via Configurables? Keep in
                     // sync with O2::TranslateASCIIintoObjArray
    if (!oa) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
  } else if (bFileIsInCCDB) {

    // d) Handle the CCDB case: Remember that here I do not access the file,
    // instead directly object in that file.
    //    My home dir in CCDB: https://alice-ccdb.cern.ch/browse/Users/a/abilandz/
    //    Inspired by:
    //    1. Discussion at:
    //    https://alice-talk.web.cern.ch/t/access-to-lhc-filling-scheme/1073/17
    //    2. See also:
    //    https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/efficiencyGlobal.cxx
    //    https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/efficiencyPerRun.cxx
    //    3. O2 Analysis Tutorial 2.0:
    //    https://indico.cern.ch/event/1267433/timetable/#20230417.detailed

    ccdb->setURL("http://alice-ccdb.cern.ch");
    oa = reinterpret_cast<TObjArray*>(ccdb->get<TObjArray>(
      TString(filePath)
        .ReplaceAll("/alice-ccdb.cern.ch/", "")
        .Data()));

    if (!oa) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
  } else {

    // e) Handle the local case:
    // Check if the external ROOT file exists at specified path:
    if (gSystem->AccessPathName(filePath, kFileExists)) {
      LOGF(info,
           "\033[1;33m if(gSystem->AccessPathName(filePath,kFileExists)), "
           "filePath = %s \033[0m",
           filePath);
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    oaFile = TFile::Open(filePath, "READ");
    if (!oaFile) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    oaFile->GetObject(
      "labels", oa); // TBI 20230530 hardcoded name of TObjArray is "labels",
                     // perhaps I can do this also via Configurables? Keep in
                     // sync with O2::TranslateASCIIintoObjArray
    if (!oa) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

  } // else {

  return oa;

} // TObjArray* GetObjArrayWithLabels(const char *filePath)

//============================================================

void StoreLabelsInPlaceholder()
{
  // Storal all Test0 labels in the temporary placeholder.

  // a) Initialize all counters;
  // b) Fetch TObjArray with labels from an external file;
  // c) Book the placeholder fTest0LabelsPlaceholder for all labels;
  // d) Finally, store the labels from external source into placeholder.

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Initialize all counters;
  Int_t counter[gMaxCorrelator] = {0}; // is this safe?
  for (Int_t o = 0; o < gMaxCorrelator; o++) {
    counter[o] = 0;
  } // now it's safe :-)

  // b) Fetch TObjArray with labels from an external file:
  TObjArray* oa = GetObjArrayWithLabels(fFileWithLabels.Data());
  if (!oa) {
    LOGF(info, "\033[1;33m fFileWithLabels = %s \033[0m",
         fFileWithLabels.Data());
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  // c) Book the placeholder fTest0LabelsPlaceholder for all labels:
  Int_t nLabels = oa->GetEntries();
  fTest0LabelsPlaceholder =
    new TH1I("fTest0LabelsPlaceholder",
             Form("placeholder for all labels, %d in total", nLabels),
             nLabels, 0, nLabels);
  fTest0LabelsPlaceholder->SetStats(kFALSE);

  // d) Finally, store the labels from external source into placeholder:
  Int_t bin = 1; // used only for fTest0LabelsPlaceholder
  Int_t order = -44;
  for (Int_t e = 0; e < nLabels; e++) {
    TObjArray* temp = TString(oa->At(e)->GetName()).Tokenize(" ");
    if (!temp) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    order = temp->GetEntries();
    delete temp; // yes, otherwise it's a memory leak
    if (0 == order) {
      continue;
    } // empty lines, or the label format which is not supported
    // 1-p => 0, 2-p => 1, etc.:
    t0_a.fTest0Labels[order - 1][counter[order - 1]] =
      new TString(oa->At(e)->GetName()); // okay...
    fTest0LabelsPlaceholder->GetXaxis()->SetBinLabel(
      bin++, t0_a.fTest0Labels[order - 1][counter[order - 1]]->Data());
    // cout<<__LINE__<<":
    // "<<t0_a.fTest0Labels[order-1][counter[order-1]]->Data()<<endl;
    counter[order - 1]++;
    // cout<<TString(line).Data()<<endl;
    // cout<<oa->GetEntries()<<endl;
  } // for(Int_t e=0; e<nLabels; e++)

} // void StoreLabelsInPlaceholder()

//============================================================

Bool_t RetrieveCorrelationsLabels()
{
  // Generate the labels of all correlations of interest, i.e. retrieve them
  // from TH1I *fTest0LabelsPlaceholder

  if (fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  Int_t counter[gMaxCorrelator] = {0}; // is this safe?
  for (Int_t o = 0; o < gMaxCorrelator; o++) {
    counter[o] = 0;
  } // now it's safe :-)

  Int_t nBins = fTest0LabelsPlaceholder->GetXaxis()->GetNbins();

  Int_t order = -44;
  for (Int_t b = 1; b <= nBins; b++) {
    TObjArray* oa = TString(fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b))
                      .Tokenize(" ");
    if (!oa) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    order = oa->GetEntries();
    delete oa; // yes, otherwise it's a memory leak
    if (0 == order) {
      continue;
    } // empty lines, or the label format which is not supported
    // 1-p => 0, 2-p => 1, etc.:
    t0_a.fTest0Labels[order - 1][counter[order - 1]] = new TString(
      fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b)); // okay...
    // cout<<__LINE__<<":
    // "<<fTest0Labels[order-1][counter[order-1]]->Data()<<endl; sleep(1);
    counter[order - 1]++;
  } // for(Int_t b=1;b<=nBins;b++)

  return kTRUE;

} // Bool_t RetrieveCorrelationsLabels()

//============================================================

TObject* GetObjectFromList(TList* list,
                           const Char_t* objectName) // Last update: 20210918
{
  // Get TObject pointer from TList, even if it's in some nested TList. Foreseen
  // to be used to fetch histograms or profiles from files directly. Some ideas
  // taken from TCollection::ls() If you have added histograms directly to files
  // (without TList's), then you can fetch them directly with
  // file->Get("hist-name").

  // Usage: TH1D *hist = (TH1D*)
  // GetObjectFromList("some-valid-TList-pointer","some-object-name");

  // Example 1:
  // GetObjectFromList("some-valid-TList-pointer","some-object-name")->Draw();
  // // yes, for histograms and profiles this is just fine, at least in
  // interpreted code

  // To do:
  // a) Check if I can make it working in compiled mode.
  // b) If I have objects with same name, nested in different TLists, what then?

  // Insanity checks:
  if (!list) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }
  if (!objectName) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }
  if (0 == list->GetEntries()) {
    return NULL;
  }

  // The object is in the current base list:
  TObject* objectFinal =
    list->FindObject(objectName); // final object I am after
  if (objectFinal)
    return objectFinal;

  // Search for object recursively in the nested lists:
  TObject* objectIter; // iterator object in the loop below
  TIter next(list);
  while (
    (objectIter = next())) // double round braces are to silent the warnings
  {
    if (TString(objectIter->ClassName()).EqualTo("TList")) {
      objectFinal = GetObjectFromList(reinterpret_cast<TList*>(objectIter), objectName);
      if (objectFinal)
        return objectFinal;
    }
  } // while(objectIter = next())

  return NULL;

} // TObject* GetObjectFromList(TList *list, Char_t *objectName)

//============================================================

Double_t Weight(const Double_t& value,
                const char* variable) // value, [phi,pt,eta]
{
  // Determine particle weight.

  // Basic protection:
  if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") ||
        TString(variable).EqualTo("eta"))) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  Int_t ppe = 0; // [phi,pt,eta]
  if (TString(variable).EqualTo("pt")) {
    ppe = 1;
  }
  if (TString(variable).EqualTo("eta")) {
    ppe = 2;
  }

  if (!pw_a.fWeightsHist[ppe]) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  Int_t bin = pw_a.fWeightsHist[ppe]->FindBin(value);
  Double_t weight = 0.;
  if (bin > pw_a.fWeightsHist[ppe]->GetNbinsX()) {
    weight = 0.; // we are in the overflow, ignore this particle TBI_20210524 is
                 // this really the correct procedure?
  } else {
    weight = pw_a.fWeightsHist[ppe]->GetBinContent(bin);
  }

  return weight;

} // Weight(const Double_t &value, const char *variable) // value, [phi,pt,eta]

//============================================================

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_MEMBERFUNCTIONS_H_
