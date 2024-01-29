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

//============================================================

void BookBaseList()
{
  // ...

  if (tc.fVerbose) {
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
                                    Form("fTaskName = %s", tc.fTaskName.Data()));

  fBasePro->GetXaxis()->SetBinLabel(eRunNumber,
                                    Form("fRunNumber = %s", tc.fRunNumber.Data()));

  fBasePro->GetXaxis()->SetBinLabel(eVerbose, "fVerbose");
  fBasePro->Fill(eVerbose - 0.5, (Int_t)tc.fVerbose);

  fBasePro->GetXaxis()->SetBinLabel(eVerboseForEachParticle, "fVerboseForEachParticle");
  fBasePro->Fill(eVerboseForEachParticle - 0.5, (Int_t)tc.fVerboseForEachParticle);

  fBasePro->GetXaxis()->SetBinLabel(eDoAdditionalInsanityChecks, "fDoAdditionalInsanityChecks");
  fBasePro->Fill(eDoAdditionalInsanityChecks - 0.5, (Int_t)tc.fDoAdditionalInsanityChecks);

  fBasePro->GetXaxis()->SetBinLabel(eUseCCDB, "fUseCCDB");
  fBasePro->Fill(eUseCCDB - 0.5, (Int_t)tc.fUseCCDB);

  fBasePro->GetXaxis()->SetBinLabel(eProcessRemainingEvents, "fProcessRemainingEvents");
  fBasePro->Fill(eProcessRemainingEvents - 0.5, (Int_t)tc.fProcessRemainingEvents);

  fBasePro->GetXaxis()->SetBinLabel(eWhatToProcess,
                                    Form("WhatToProcess = %s", tc.fWhatToProcess.Data()));

  fBasePro->GetXaxis()->SetBinLabel(eRandomSeed, "fRandomSeed");
  fBasePro->Fill(eRandomSeed - 0.5, (Int_t)tc.fRandomSeed);

  fBasePro->GetXaxis()->SetBinLabel(eUseFisherYates, "fUseFisherYates");
  fBasePro->Fill(eUseFisherYates - 0.5, (Int_t)tc.fUseFisherYates);

  fBasePro->GetXaxis()->SetBinLabel(eFixedNumberOfRandomlySelectedTracks, "fFixedNumberOfRandomlySelectedTracks");
  fBasePro->Fill(eFixedNumberOfRandomlySelectedTracks - 0.5, (Int_t)tc.fFixedNumberOfRandomlySelectedTracks);

  fBaseList->Add(fBasePro);

} // void BookBaseList()

//============================================================

void WhatToProcess()
{
  // Set here what to process: only rec, both rec and sim, only sim.
  // Use in combination with configurable cfWhatToProcess.
  // TBI 20231017 I call this function, but it still has no desired effect, until I can call PROCESS_SWITCH(  ) by passing a variable, instead only literals 'true' or 'false', as it is now

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  if (tc.fWhatToProcess.EqualTo("Rec")) {
    gProcessRec = true;
  } else if (tc.fWhatToProcess.EqualTo("RecSim")) {
    gProcessRecSim = true;
  } else if (tc.fWhatToProcess.EqualTo("Sim")) {
    gProcessSim = true;
  } else {
    LOGF(info, "\033[1;32m This option is not supported! tc.fWhatToProcess = %s \033[0m", tc.fWhatToProcess.Data());
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  // Make sure that only one of these flags is set to kTRUE.
  // This is needed, becase these flags are used in PROCESS_SWITCH,
  // and if 2 or more are kTRUE, then corresponding process function
  // is executed over ALL data, then another process(...) function, etc.
  if ((Int_t)gProcessRec + (Int_t)gProcessRecSim + (Int_t)gProcessSim > 1) {
    LOGF(info, "\033[1;32m Only one flag can be kTRUE: gProcessRec = %d, gProcessRecSim = %d, gProcessSim = %d \033[0m", (Int_t)gProcessRec, (Int_t)gProcessRecSim, (Int_t)gProcessSim);
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

} // WhatToProcess()

//============================================================

void DefaultConfiguration()
{
  // Default task configuration.
  // a) Default values are hardcoded as Configurables in the file
  // MuPa-Configurables.h

  // b) If corresponding fields are available in an external json file at run time, the default values hardcoded here are
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

  // c) Scientific notation is NOT supported in json file. E.g. if you have
  //            "cSelectedTracks_max": "1e3",
  //    that setting and ALL other ones in json are silently ignored.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // Configurable<string> cfTaskName{ ... }
  tc.fTaskName = TString(cfTaskName);

  // Configurable<bool> cfVerbose{ ... }
  tc.fVerbose = cfVerbose;

  // Configurable<bool> cfVerboseForEachParticle{ ... }
  tc.fVerboseForEachParticle = cfVerboseForEachParticle;

  // Configurable<bool> cfDoAdditionalInsanityChecks{ ... }
  tc.fDoAdditionalInsanityChecks = cfDoAdditionalInsanityChecks;

  // Configurable<bool> cfUseCCDB{ ... }
  tc.fUseCCDB = cfUseCCDB;

  // Configurable<string> cfWhatToProcess{ ... )
  tc.fWhatToProcess = TString(cfWhatToProcess);

  // Configurable<unsigned int> cfRandomSeed{ ... )
  tc.fRandomSeed = cfRandomSeed;

  // Configurable<bool> cfUseFisherYates{ ... }
  tc.fUseFisherYates = cfUseFisherYates;

  // Configurable<int> cfFixedNumberOfRandomlySelectedTracks{ ... }
  tc.fFixedNumberOfRandomlySelectedTracks = cfFixedNumberOfRandomlySelectedTracks;

  // ...

  // Configurable<bool> cfCalculateCorrelations{ ... };
  fCalculateCorrelations = cfCalculateCorrelations;

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

  // *) Nested loops:
  // Configurable<string> cfCalculateNestedLoops{ ... }
  fCalculateNestedLoops = cfCalculateNestedLoops;

  // Configurable<string> cfCalculateCustomNestedLoops{ ... }
  fCalculateCustomNestedLoop = cfCalculateCustomNestedLoops;

  // *) TBI 20231108 not ported yet:
  // task->SetCalculateQvector(kTRUE);
  fCalculateQvector = kTRUE;

} // void DefaultConfiguration()

//============================================================

void DefaultBooking()
{
  // Set here which histograms are booked by default.

  // a) Event histograms;
  // b) Particle histograms;
  // c) QA;

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Event histograms:
  // Each default setting can be overuled e.g. with:
  // task->SetBookEventHistograms("NumberOfEvents",kFALSE);
  ceh_a.fBookEventHistograms[eNumberOfEvents] = kTRUE;
  ceh_a.fBookEventHistograms[eTotalMultiplicity] = kTRUE;
  ceh_a.fBookEventHistograms[eSelectedTracks] = kTRUE;
  ceh_a.fBookEventHistograms[eMultFV0M] = kTRUE;
  ceh_a.fBookEventHistograms[eMultFT0M] = kTRUE;
  ceh_a.fBookEventHistograms[eMultTPC] = kTRUE;
  ceh_a.fBookEventHistograms[eMultNTracksPV] = kTRUE;
  ceh_a.fBookEventHistograms[eCentrality] = kTRUE;
  ceh_a.fBookEventHistograms[eVertex_x] = kTRUE;
  ceh_a.fBookEventHistograms[eVertex_y] = kTRUE;
  ceh_a.fBookEventHistograms[eVertex_z] = kTRUE;
  ceh_a.fBookEventHistograms[eNContributors] = kTRUE;
  ceh_a.fBookEventHistograms[eImpactParameter] = kTRUE;

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

  // TBI 20240114 If some of these values are going to change frequently, add support for them in MuPa-Configurables.h,
  // in the same way I did it for DefaultCuts().

  // a) Default binning for event histograms;
  // b) Default binning for particle histograms;
  // c) Default binning for results histograms;
  // d) Variable-length binning set via MuPa-Configurables.h.

  if (tc.fVerbose) {
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
  // ... TBI 20240120
  ceh_a.fEventHistogramsBins[eMultFV0M][0] = 10000;
  ceh_a.fEventHistogramsBins[eMultFV0M][1] = 0.;
  ceh_a.fEventHistogramsBins[eMultFV0M][2] = 10000.;
  // ... TBI 20240120
  ceh_a.fEventHistogramsBins[eMultFT0M][0] = 10000;
  ceh_a.fEventHistogramsBins[eMultFT0M][1] = 0.;
  ceh_a.fEventHistogramsBins[eMultFT0M][2] = 10000.;
  // ... TBI 20240120
  ceh_a.fEventHistogramsBins[eMultTPC][0] = 10000;
  ceh_a.fEventHistogramsBins[eMultTPC][1] = 0.;
  ceh_a.fEventHistogramsBins[eMultTPC][2] = 10000.;
  // ... TBI 20240120
  ceh_a.fEventHistogramsBins[eMultNTracksPV][0] = 10000;
  ceh_a.fEventHistogramsBins[eMultNTracksPV][1] = 0.;
  ceh_a.fEventHistogramsBins[eMultNTracksPV][2] = 10000.;
  // task->SetEventHistogramsBins("Centrality",100,0.,100.);
  ceh_a.fEventHistogramsBins[eCentrality][0] = 110; // intentionally, because if centrality is not determined, it's set to 105.0 at the moment
  ceh_a.fEventHistogramsBins[eCentrality][1] = 0.;
  ceh_a.fEventHistogramsBins[eCentrality][2] = 110.;
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
  // task->SetEventHistogramsBins("ImpactParameter",1000,0.,1000.);
  ceh_a.fEventHistogramsBins[eImpactParameter][0] = 1000000;
  ceh_a.fEventHistogramsBins[eImpactParameter][1] = 0.;
  ceh_a.fEventHistogramsBins[eImpactParameter][2] = 1.; // TBI 20231031 check this, i do not know in which units IP is stored

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

  cph_a.fParticleHistogramsBins[etpcNClsCrossedRows][0] = 200;
  cph_a.fParticleHistogramsBins[etpcNClsCrossedRows][1] = 0.;
  cph_a.fParticleHistogramsBins[etpcNClsCrossedRows][2] = 200.;

  cph_a.fParticleHistogramsBins[eDCA_xy][0] = 2000;
  cph_a.fParticleHistogramsBins[eDCA_xy][1] = -10.;
  cph_a.fParticleHistogramsBins[eDCA_xy][2] = 10.;

  cph_a.fParticleHistogramsBins[eDCA_z][0] = 2000;
  cph_a.fParticleHistogramsBins[eDCA_z][1] = -10.;
  cph_a.fParticleHistogramsBins[eDCA_z][2] = 10.;

  // c) Default binning for results histograms:
  //    Remark: These bins apply to following categories fCorrelationsPro, fNestedLoopsPro, fTest0Pro, and fResultsHist.
  rh_a.fResultsHistogramsFixedLengthBins[AFO_INTEGRATED][0] = 1;
  rh_a.fResultsHistogramsFixedLengthBins[AFO_INTEGRATED][1] = 0.;
  rh_a.fResultsHistogramsFixedLengthBins[AFO_INTEGRATED][2] = 1.;

  rh_a.fResultsHistogramsFixedLengthBins[AFO_MULTIPLICITY][0] = 1000;
  rh_a.fResultsHistogramsFixedLengthBins[AFO_MULTIPLICITY][1] = 0.;
  rh_a.fResultsHistogramsFixedLengthBins[AFO_MULTIPLICITY][2] = 5000.;

  rh_a.fResultsHistogramsFixedLengthBins[AFO_CENTRALITY][0] = 100;
  rh_a.fResultsHistogramsFixedLengthBins[AFO_CENTRALITY][1] = 0.;
  rh_a.fResultsHistogramsFixedLengthBins[AFO_CENTRALITY][2] = 100.;

  rh_a.fResultsHistogramsFixedLengthBins[AFO_PT][0] = cph_a.fParticleHistogramsBins[ePt][0];
  rh_a.fResultsHistogramsFixedLengthBins[AFO_PT][1] = cph_a.fParticleHistogramsBins[ePt][1];
  rh_a.fResultsHistogramsFixedLengthBins[AFO_PT][2] = cph_a.fParticleHistogramsBins[ePt][2];

  rh_a.fResultsHistogramsFixedLengthBins[AFO_ETA][0] = cph_a.fParticleHistogramsBins[eEta][0];
  rh_a.fResultsHistogramsFixedLengthBins[AFO_ETA][1] = cph_a.fParticleHistogramsBins[eEta][1];
  rh_a.fResultsHistogramsFixedLengthBins[AFO_ETA][2] = cph_a.fParticleHistogramsBins[eEta][2];

  // d) Variable-length binning set via MuPa-Configurables.h:
  if (cUseVariableLength_mult_bins) {
    rh_a.fUseResultsHistogramsVariableLengthBins[AFO_MULTIPLICITY] = kTRUE;
    rh_a.fResultsHistogramsVariableLengthBinsString[AFO_MULTIPLICITY] = cVariableLength_mult_bins;
    this->CastStringIntoArray(AFO_MULTIPLICITY);
  }
  if (cUseVariableLength_cent_bins) {
    rh_a.fUseResultsHistogramsVariableLengthBins[AFO_CENTRALITY] = kTRUE;
    rh_a.fResultsHistogramsVariableLengthBinsString[AFO_CENTRALITY] = cVariableLength_cent_bins;
    this->CastStringIntoArray(AFO_CENTRALITY);
  }
  if (cUseVariableLength_pt_bins) {
    rh_a.fUseResultsHistogramsVariableLengthBins[AFO_PT] = kTRUE;
    rh_a.fResultsHistogramsVariableLengthBinsString[AFO_PT] = cVariableLength_pt_bins;
    this->CastStringIntoArray(AFO_PT);
  }
  if (cUseVariableLength_eta_bins) {
    rh_a.fUseResultsHistogramsVariableLengthBins[AFO_ETA] = kTRUE;
    rh_a.fResultsHistogramsVariableLengthBinsString[AFO_ETA] = cVariableLength_eta_bins;
    this->CastStringIntoArray(AFO_ETA);
  }

} // void DefaultBinning()

//============================================================

void CastStringIntoArray(Int_t AFO)
{
  // Temporary function, to be removed eventually. Here temporarily I am casting e.g. a string "1.0,2.0,5.0" into corresponding TArrayD.
  // TBI 20240114 This function is used until I figure out how to pass array directly via configurable.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m Casting a string %s into TArrayD .... \033[0m", rh_a.fResultsHistogramsVariableLengthBinsString[AFO].Data());
  }

  TObjArray* oa = rh_a.fResultsHistogramsVariableLengthBinsString[AFO].Tokenize(",");
  if (!oa) {
    LOGF(fatal, "in function \033[1;31m%s at line %d \n fResultsHistogramsVariableLengthBinsString[AFO] = %s\033[0m", __PRETTY_FUNCTION__, __LINE__, rh_a.fResultsHistogramsVariableLengthBinsString[AFO].Data());
  }
  Int_t nEntries = oa->GetEntries();
  rh_a.fResultsHistogramsVariableLengthBins[AFO] = new TArrayD(nEntries);
  for (Int_t i = 0; i < nEntries; i++) {
    // cout<< TString(oa->At(i)->GetName()).Atof() <<endl;
    rh_a.fResultsHistogramsVariableLengthBins[AFO]->AddAt(TString(oa->At(i)->GetName()).Atof(), i);
  }
  delete oa; // yes, otherwise it's a memory leak

  if (tc.fVerbose) {
    for (Int_t i = 0; i < rh_a.fResultsHistogramsVariableLengthBins[AFO]->GetSize(); i++) {
      LOGF(info, "\033[1;32m [%d] : %f \033[0m", i, rh_a.fResultsHistogramsVariableLengthBins[AFO]->At(i));
    }
  }

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m Done! \033[0m");
  }

} // void CastStringIntoArray(Int_t AFO)

//============================================================

void DefaultCuts()
{
  // Define default cuts. Default cuts are hardwired in MuPa-Configurables.h.

  // a) Default event cuts;
  // b) Default particle cuts.

  if (tc.fVerbose) {
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

  ceh_a.fEventCuts[eImpactParameter][eMin] =
    cImpactParameter_min; // Configurable<int>
                          // cImpactParameter_min{"cImpactParameter_min", ...
  ceh_a.fEventCuts[eImpactParameter][eMax] =
    cImpactParameter_max; // Configurable<int>
                          // cImpactParameter_max{"cImpactParameter_max", ...

  // ...

  // b) Default particle cuts:

} // void DefaultCuts()

//============================================================

void InsanityChecks()
{
  // Do insanity checks on configuration, booking, binning and cuts.

  // *) Insanity checks on configuration;
  // *) Insanity checks on booking;
  // *) Insanity checks on binning;
  // *) Insanity checks on cuts.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // *) Insanity checks on configuration:
  if (tc.fRandomSeed < 0) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && !tc.fUseFisherYates) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  // *) Insanity checks on booking:
  // ...

  // *) Insanity checks on binning:
  // ...

  // *) Insanity checks on cuts:
  // ...

} // void InsanityChecks()

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

  if (tc.fVerbose) {
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

  if (tc.fVerbose) {
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
    "NumberOfEvents", "TotalMultiplicity", "SelectedTracks", "MultFV0M", "MultFT0M", "MultTPC", "MultNTracksPV",
    "Centrality", "Vertex_x", "Vertex_y",
    "Vertex_z", "NContributors", "ImpactParameter"}; // keep in sync. with enum eEventHistograms
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
      if ((gProcessRec && rs == eSim) || (gProcessSim && rs == eRec)) {
        continue; // if I am analyzing only reconstructed data, do not book histos for simulated, and vice versa.
      }
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

  if (tc.fVerbose) {
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
      if ((gProcessRec && rs == eSim) || (gProcessSim && rs == eRec)) {
        continue; // if I am analyzing only reconstructed data, do not book histos for simulated, and vice versa.
      }
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

  if (tc.fVerbose) {
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
  // c) Histograms;
  // d) Few quick insanity checks on booking.

  if (tc.fVerbose) {
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

  // c) Histograms:
  for (Int_t k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]
  {
    for (Int_t n = 0; n < gMaxHarmonic; n++) // harmonic
    {
      for (Int_t v = 0; v < eAsFunctionOf_N;
           v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pt,4=eta]
      {
        // per demand, variable-length binning:
        if (rh_a.fUseResultsHistogramsVariableLengthBins[v]) {
          c_a.fCorrelationsPro[k][n][v] = new TProfile(Form("fCorrelationsPro[%d][%d][%s]", k, n, rh_a.fResultsHistogramsRawName[v].Data()), "some title", rh_a.fResultsHistogramsVariableLengthBins[v]->GetSize() - 1, rh_a.fResultsHistogramsVariableLengthBins[v]->GetArray());
        } else {
          // the default fixed-length binning:
          c_a.fCorrelationsPro[k][n][v] = new TProfile(Form("fCorrelationsPro[%d][%d][%s]", k, n, rh_a.fResultsHistogramsRawName[v].Data()),
                                                       "some title", (Int_t)rh_a.fResultsHistogramsFixedLengthBins[v][0], rh_a.fResultsHistogramsFixedLengthBins[v][1], rh_a.fResultsHistogramsFixedLengthBins[v][2]);
        }

        c_a.fCorrelationsPro[k][n][v]->SetStats(kFALSE);
        c_a.fCorrelationsPro[k][n][v]->Sumw2();
        c_a.fCorrelationsPro[k][n][v]->GetXaxis()->SetTitle(rh_a.fResultsHistogramsXaxisTitle[v].Data());
        c_a.fCorrelationsPro[k][n][v]->GetYaxis()->SetTitle(Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1), oVariable[k].Data()));
        fCorrelationsList->Add(c_a.fCorrelationsPro[k][n][v]);
      }
    } // for (Int_t n = 0; n < gMaxHarmonic; n++) // harmonic
  }   // for (Int_t k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]

  // d) Few quick insanity checks on booking:
  if (c_a.fCorrelationsPro[0][0][AFO_INTEGRATED] && !TString(c_a.fCorrelationsPro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsHistogramsXaxisTitle[eAsFunctionOf_N]
  }
  if (c_a.fCorrelationsPro[0][0][AFO_PT] && !TString(c_a.fCorrelationsPro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("p_{T}")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsHistogramsXaxisTitle[eAsFunctionOf_N]
  }

} // BookCorrelationsHistograms()

//============================================================

void BookWeightsHistograms()
{
  // Book all objects for particle weights.

  // a) Book the profile holding flags;
  // b) Common local labels;
  // c) Histograms.

  if (tc.fVerbose) {
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
  // b) Common local labels (keep 'em in sync with BookCorrelationsHistograms());
  // c) Book what needs to be booked;
  // d) Few quick insanity checks on booking.

  if (tc.fVerbose) {
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

  // c) Book what needs to be booked:
  for (Int_t k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]
  {
    for (Int_t n = 0; n < gMaxHarmonic; n++) // harmonic
    {
      for (Int_t v = 0; v < eAsFunctionOf_N;
           v++) // variable [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pt,4=eta]
      {

        // if(PTKINE == v  && !fCalculatePtCorrelations){continue;}
        // if(ETAKINE == v  && !fCalculateEtaCorrelations){continue;}

        // per demand, variable-length binning:
        if (rh_a.fUseResultsHistogramsVariableLengthBins[v]) {
          nl_a.fNestedLoopsPro[k][n][v] = new TProfile(
            Form("fNestedLoopsPro[%d][%d][%d]", k, n, v),
            Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1),
                 oVariable[k].Data()),
            rh_a.fResultsHistogramsVariableLengthBins[v]->GetSize() - 1, rh_a.fResultsHistogramsVariableLengthBins[v]->GetArray());
        } else {
          // the default fixed-length binning:
          nl_a.fNestedLoopsPro[k][n][v] = new TProfile(
            Form("fNestedLoopsPro[%d][%d][%d]", k, n, v),
            Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1),
                 oVariable[k].Data()),
            (Int_t)rh_a.fResultsHistogramsFixedLengthBins[v][0], rh_a.fResultsHistogramsFixedLengthBins[v][1], rh_a.fResultsHistogramsFixedLengthBins[v][2]);
        } // else

        nl_a.fNestedLoopsPro[k][n][v]->SetStats(kFALSE);
        nl_a.fNestedLoopsPro[k][n][v]->Sumw2();
        nl_a.fNestedLoopsPro[k][n][v]->GetXaxis()->SetTitle(
          rh_a.fResultsHistogramsXaxisTitle[v].Data());

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
    }   // for (Int_t n = 0; n < gMaxHarmonic; n++) // harmonic
  }     // for (Int_t k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]

  // d) Few quick insanity checks on booking:
  if (nl_a.fNestedLoopsPro[0][0][AFO_INTEGRATED] && !TString(nl_a.fNestedLoopsPro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsHistogramsXaxisTitle[eAsFunctionOf_N]
  }
  if (nl_a.fNestedLoopsPro[0][0][AFO_PT] && !TString(nl_a.fNestedLoopsPro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("p_{T}")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsHistogramsXaxisTitle[eAsFunctionOf_N]
  }

} // void BookNestedLoopsHistograms()

//============================================================

void BookTest0Histograms()
{
  // Book all Test0 histograms.

  // a) Book the profile holding flags;
  // b) Book placeholder and make sure all labels are stored in the placeholder;
  // c) Retrieve labels from placeholder;
  // d) Book what needs to be booked;
  // e) Few quick insanity checks on booking.

  if (tc.fVerbose) {
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
  this->StoreLabelsInPlaceholder();
  if (fTest0LabelsPlaceholder) {
    fTest0List->Add(fTest0LabelsPlaceholder);
  }

  // c) Retrieve labels from placeholder:
  if (!(this->RetrieveCorrelationsLabels())) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  // d) Book what needs to be booked:
  for (Int_t mo = 0; mo < gMaxCorrelator; mo++) {
    for (Int_t mi = 0; mi < gMaxIndex; mi++) {
      if (!t0_a.fTest0Labels[mo][mi]) {
        continue;
      }
      {
        for (Int_t v = 0; v < eAsFunctionOf_N; v++) {
          /*
                // decide what is booked, then later valid pointer to fCorrelationsPro[k][n][v] is used as a boolean, in the standard way:
                if(AFO_INTEGRATED == v && fDoNotCalculateCorrelationsAsFunctionOf[AFO_INTEGRATED]){continue;}
                if(AFO_MULTIPLICITY == v && fDoNotCalculateCorrelationsAsFunctionOf[AFO_MULTIPLICITY]){continue;}
                if(AFO_CENTRALITY == v && fDoNotCalculateCorrelationsAsFunctionOf[AFO_CENTRALITY]){continue;}
                if(AFO_PT == v && fDoNotCalculateCorrelationsAsFunctionOf[AFO_PT]){continue;}
                if(AFO_ETA == v && fDoNotCalculateCorrelationsAsFunctionOf[AFO_ETA]){continue;}
          */

          // per demand, variable-length binning:
          if (rh_a.fUseResultsHistogramsVariableLengthBins[v]) {
            t0_a.fTest0Pro[mo][mi][v] = new TProfile(Form("fTest0Pro[%d][%d][%s]", mo, mi, rh_a.fResultsHistogramsRawName[v].Data()), t0_a.fTest0Labels[mo][mi]->Data(), rh_a.fResultsHistogramsVariableLengthBins[v]->GetSize() - 1, rh_a.fResultsHistogramsVariableLengthBins[v]->GetArray());
          } else {
            // the default fixed-length binning:
            t0_a.fTest0Pro[mo][mi][v] = new TProfile(Form("fTest0Pro[%d][%d][%s]", mo, mi, rh_a.fResultsHistogramsRawName[v].Data()), t0_a.fTest0Labels[mo][mi]->Data(), (Int_t)rh_a.fResultsHistogramsFixedLengthBins[v][0], rh_a.fResultsHistogramsFixedLengthBins[v][1], rh_a.fResultsHistogramsFixedLengthBins[v][2]);
          }
          t0_a.fTest0Pro[mo][mi][v]->SetStats(kFALSE);
          t0_a.fTest0Pro[mo][mi][v]->Sumw2();
          t0_a.fTest0Pro[mo][mi][v]->GetXaxis()->SetTitle(rh_a.fResultsHistogramsXaxisTitle[v].Data());
          /*
                if(fUseFixedNumberOfRandomlySelectedParticles && 1==v) // just a warning for the meaning of multiplicity in this special case
                {
                 fTest0Pro[mo][mi][1]->GetXaxis()->SetTitle("WARNING: for each multiplicity, fFixedNumberOfRandomlySelectedParticles is selected randomly in Q-vector");
                }
          */
          fTest0List->Add(t0_a.fTest0Pro[mo][mi][v]); // yes, this has to be here
        }                                             // for(Int_t v=0;v<eAsFunctionOf_N;v++) // variable, see content of enum eAsFunctionOf
      }                                               // if(fTest0Labels[mo][mi])
    }                                                 // for(Int_t mi=0;mi<gMaxIndex;mi++)
  }                                                   // for(Int_t mo=0;mo<gMaxCorrelator;mo++)

  // e) Few quick insanity checks on booking:
  if (t0_a.fTest0Pro[0][0][AFO_INTEGRATED] && !TString(t0_a.fTest0Pro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsHistogramsXaxisTitle[eAsFunctionOf_N]
  }
  if (t0_a.fTest0Pro[0][0][AFO_PT] && !TString(t0_a.fTest0Pro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("p_{T}")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsHistogramsXaxisTitle[eAsFunctionOf_N]
  }

} // void BookTest0Histograms()

//============================================================

void BookResultsHistograms()
{
  // Book all results histograms.

  // a) Book the profile holding flags;
  // *) ...

  if (tc.fVerbose) {
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

template <typename T>
void Preprocess(T const& collision)
{
  // Do all thingies before starting to process data (e.g. count number of events, fetch the run number, etc.).

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // *) If I reached max number of events, ignore the remaining collisions:
  if (MaxNumberOfEvents()) {
    tc.fProcessRemainingEvents = kFALSE;
  }

  // *) Determine and propagate run number info to already booked objects:
  if (!tc.fRunNumberIsDetermined) {
    DetermineAndPropagateRunNumber(collision);
  }
  if (tc.fDoAdditionalInsanityChecks && tc.fRunNumberIsDetermined) {
    CheckCurrentRunNumber(collision);
  }

  // *) Fetch the weights for this particular run number. Do it only once.
  //    TBI 20231012 If eventualy I can access programatically run number in init(...) at run time, this shall go there.
  if (!fParticleWeightsAreFetched) {
    if (pw_a.fUseWeights[wPHI] || pw_a.fUseWeights[wPT] || pw_a.fUseWeights[wETA]) {
      GetParticleWeights();
      fParticleWeightsAreFetched = kTRUE;
    }
  }

} // template <typename T> void Preprocess(T const& collision)

//============================================================

template <typename T>
void DetermineAndPropagateRunNumber(T const& collision)
{
  // Determine and propagate run number info to already booked objects, wherever it's relevant.
  // Make sure in process(...) that this function is called only once.

  // TBI 20231018 At the moment I can access run number info only in process(...) via collision->bc().runNumber(), but not in init(...)
  // Once I can access run number info in init(...), this function shall be called in init(...), not in process(...)

  // a) Determine run number;
  // b) Propagate run number to all booked objects, wherever that info is relevant.

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // a) Determine run number for reconstructed data:
  tc.fRunNumber = Form("%d", collision.bc().runNumber()); // implemented for both aod::Collision and aod::McCollision, so I can use it straight, as long as I have subscribed to aod::BCs
  if (tc.fRunNumber.EqualTo("")) {
    LOGF(error, "\033[1;33m%s fRunNumber is empty, collision->bc().runNumber() failed...\033[0m", __PRETTY_FUNCTION__);
    LOGF(fatal, "collision->bc().runNumber() = %d", collision.bc().runNumber());
  }
  tc.fRunNumberIsDetermined = kTRUE;

  // b) Propagate run number to all booked objects, wherever that info is relevant:
  fBasePro->GetXaxis()->SetBinLabel(eRunNumber, Form("tc.fRunNumber = %s", tc.fRunNumber.Data()));
  // ...

} // template <typename T> void DetermineAndPropagateRunNumber(T const& collision)

//============================================================

template <typename T>
void CheckCurrentRunNumber(T const& collision)
{
  // Insanity check for the current run number.

  if (!tc.fRunNumber.EqualTo(Form("%d", collision.bc().runNumber()))) {
    LOGF(error, "\033[1;33m%s Run number changed within process(). This most likely indicates that a given masterjob is processing 2 or more different runs in one go.\033[0m", __PRETTY_FUNCTION__);
    LOGF(fatal, "tc.fRunNumber = %s, collision.bc().runNumber() = %d", tc.fRunNumber.Data(), collision.bc().runNumber());
  }

} // template <typename T> void CheckCurrentRunNumber(T const& collision)

//============================================================

void ResetEventByEventQuantities()
{
  // Reset all global event-by-event quantities here:

  // a) Event-by-event quantities;
  // b) Q-vectors;
  // c) Reset ebe containers for nested loops;
  // d) Fisher-Yates algorithm.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Event-by-event quantities:
  fSelectedTracks = 0;
  fCentrality = 0;

  // b) Q-vectors:
  if (fCalculateQvector) {
    ResetQ(); // generic Q-vector
    for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
      for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) // weight power
      {
        qv_a.fQvector[h][wp] = TComplex(0., 0.);
      }
    }
  } // if(fCalculateQvector)

  // c) Reset ebe containers for nested loops:
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

  // d) Fisher-Yates algorithm:
  if (tc.fUseFisherYates) {
    delete tc.fRandomIndices;
    tc.fRandomIndices = NULL;
  }

  // ... TBI 20240117 port the rest ...

} // void ResetEventByEventQuantities()

//============================================================

template <eRecSim rs, typename T1, typename T2>
Bool_t EventCuts(T1 const& collision, T2 const& tracks)
{
  // Event cuts on reconstructed and simulated data.

  // a) Event cuts on info available in reconstructed (and corresponding MC truth simulated);
  // b) Event cuts on info available only in simulated data.

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // a) Event cuts on info available in reconstructed ...:
  if constexpr (rs == eRec || rs == eRecAndSim) {
    //   *) NumberOfEvents: => cut directly in void process( ... )
    //   *) TotalMultiplicity:
    if ((tracks.size() < ceh_a.fEventCuts[eTotalMultiplicity][eMin]) ||
        (tracks.size() > ceh_a.fEventCuts[eTotalMultiplicity][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eTotalMultiplicity\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    //   *) SelectedTracks: => cut directly in void process( ... )
    //   *) Centrality: TBI
    //  if ((  TBI   < ceh_a.fEventCuts[eCentrality][eMin]) || (  TBI   >
    //  ceh_a.fEventCuts[eCentrality][eMax])) {
    //    return kFALSE;
    //  }
    //   *) Vertex_x:
    if ((collision.posX() < ceh_a.fEventCuts[eVertex_x][eMin]) ||
        (collision.posX() > ceh_a.fEventCuts[eVertex_x][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eVertex_x\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    //   *) Vertex_y:
    if ((collision.posY() < ceh_a.fEventCuts[eVertex_y][eMin]) ||
        (collision.posY() > ceh_a.fEventCuts[eVertex_y][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eVertex_y\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    //   *) Vertex_z:
    if ((collision.posZ() < ceh_a.fEventCuts[eVertex_z][eMin]) ||
        (collision.posZ() > ceh_a.fEventCuts[eVertex_z][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eVertex_z\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    //   *) NContributors:
    if ((collision.numContrib() < ceh_a.fEventCuts[eNContributors][eMin]) ||
        (collision.numContrib() > ceh_a.fEventCuts[eNContributors][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eNContributors\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    // TBI 20231106 continue here with other event cuts on reconstructed info

    // ... and corresponding MC truth simulated ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip..."); // TBI 20231106 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this collision
        return kFALSE;
      }

      // TBI 20231106 here I cat cut directly on corresponding MC truth simulated, e.g. on collision.mcCollision().posZ(), if necessary

    } // if constexpr (rs == eRecAndSim) {

  } // if constexpr (rs == eRec || rs == eRecAndSim) {

  // b) Event cuts on info available only in simulated data:
  if constexpr (rs == eSim) {
    //   *) Impact parameter:
    if ((collision.impactParameter() < ceh_a.fEventCuts[eImpactParameter][eMin]) ||
        (collision.impactParameter() > ceh_a.fEventCuts[eImpactParameter][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eImpactParameter\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    // ...
  } // if (gProcessSim) {

  return kTRUE;

} // template <eRecSim rs, typename T1, typename T2> Bool_t EventCuts(T1 const& collision, T2 const& tracks)

//============================================================

template <eRecSim rs, typename T1, typename T2>
void FillEventHistograms(T1 const& collision, T2 const& tracks, eBeforeAfter ba)
{
  // Fill all event histograms for reconstructed or simulated data.

  // a) Fill reconstructed (and corresponding MC truth simulated);
  // b) Fill only simulated.

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s eBeforeAfter = %d \033[0m", __FUNCTION__, (Int_t)ba); // just a bare function name
  }

  // a) Fill reconstructed ...:
  if constexpr (rs == eRec || rs == eRecAndSim) {
    ceh_a.fEventHistograms[eNumberOfEvents][eRec][ba]->Fill(0.5);
    ceh_a.fEventHistograms[eVertex_x][eRec][ba]->Fill(collision.posX());
    ceh_a.fEventHistograms[eVertex_y][eRec][ba]->Fill(collision.posY());
    ceh_a.fEventHistograms[eVertex_z][eRec][ba]->Fill(collision.posZ());
    ceh_a.fEventHistograms[eNContributors][eRec][ba]->Fill(collision.numContrib());
    ceh_a.fEventHistograms[eTotalMultiplicity][eRec][ba]->Fill(tracks.size()); // TBI 20231106 check and validate further
    ceh_a.fEventHistograms[eSelectedTracks][eRec][ba]->Fill(fSelectedTracks);  // TBI 20240108 this one makes sense only for eAfter
    ceh_a.fEventHistograms[eMultFT0M][eRec][ba]->Fill(collision.multFT0M());
    ceh_a.fEventHistograms[eMultFV0M][eRec][ba]->Fill(collision.multFV0M());
    ceh_a.fEventHistograms[eMultTPC][eRec][ba]->Fill(collision.multTPC());
    ceh_a.fEventHistograms[eMultNTracksPV][eRec][ba]->Fill(collision.multNTracksPV());
    ceh_a.fEventHistograms[eCentrality][eRec][ba]->Fill(fCentrality); // TBI 20240120 for the time being, I fill only default centrality

    // ... and corresponding MC truth simulated ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }
      ceh_a.fEventHistograms[eNumberOfEvents][eSim][ba]->Fill(0.5);
      ceh_a.fEventHistograms[eVertex_x][eSim][ba]->Fill(collision.mcCollision().posX());
      ceh_a.fEventHistograms[eVertex_y][eSim][ba]->Fill(collision.mcCollision().posY());
      ceh_a.fEventHistograms[eVertex_z][eSim][ba]->Fill(collision.mcCollision().posZ());
      // ceh_a.fEventHistograms[eTotalMultiplicity][eSim][ba]->Fill(tracks.size()); // TBI 20231106 check how to get corresponding MC truth info, and validate further
      // ceh_a.fEventHistograms[eSelectedTracks][eSim][ba]->Fill(fSelectedTracks); // TBI 20240108 this one makes sense only for eAfter + re-think if I really need it here
      // TBI 20240120 eMultFT0M, ..., eMultNTracksPV are not needed here
      // ceh_a.fEventHistograms[eCentrality][eSim][ba]->Fill(fCentrality); // TBI 20240120 this case is still not supported in DetermineCentrality()
    } // if constexpr (rs == eRecAndSim) {
  }   // if constexpr (rs == eRec || rs == eRecAndSim) {

  // b) Fill only simulated:
  if constexpr (rs == eSim) {
    ceh_a.fEventHistograms[eImpactParameter][eSim][ba]->Fill(collision.impactParameter()); // yes, because in this branch 'collision' is always aod::McCollision
    ceh_a.fEventHistograms[eSelectedTracks][eSim][ba]->Fill(fSelectedTracks);              // TBI 20240108 this one makes sense only for eAfter
    // ceh_a.fEventHistograms[eCentrality][eSim][ba]->Fill(fCentrality); // TBI 20240120 this case is still not supported in DetermineCentrality()
    // ceh_a.fEventHistograms[eTotalMultiplicity][eSim][ba]->Fill(tracks.size()); // TBI 20231030 check further how to use the same thing for 'sim'
  } // if constexpr (rs == eSim) {

} // template <eRecSim rs, typename T1, typename T2> void FillEventHistograms(...)

//============================================================

template <eRecSim rs, typename T>
Bool_t ParticleCuts(T const& track)
{
  // Particles cuts.

  // a) Particle cuts on info available in reconstructed (and corresponding MC truth simulated);
  // b) Particle cuts on info available only in simulated data.

  if (tc.fVerboseForEachParticle) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Particle cuts on info available in reconstructed ...:
  if constexpr (rs == eRec || rs == eRecAndSim) {
    if ((track.pt() < pt_min) || (track.pt() > pt_max)) {
      return kFALSE;
    }

    // TBI 20231107 other cuts on reconstructed track ...

    // ... and corresponding MC truth simulated ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for this track, skip...");
        return kFALSE; // TBI 20231107 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this particle
      }
      auto mcparticle = track.mcParticle();                           // corresponding MC truth simulated particle
      if ((mcparticle.pt() < pt_min) || (mcparticle.pt() > pt_max)) { // TBI 20231107 re-thing if i really cut directly on MC truth, keep it in sync with what I did in AliPhysics
        return kFALSE;
      }

      // TBI 20231107 other cuts on corresponding MC truth particle ...

    } // if constexpr (rs == eRecAndSim) {
  }   // if constexpr (rs == eRec || rs == eRecAndSim) {

  // b) Particle cuts on info available only in simulated data:
  if constexpr (rs == eSim) {
    // Remark: in this branch, 'track' is always TracksSim = aod::McParticles
    if ((track.pt() < pt_min) || (track.pt() > pt_max)) {
      return kFALSE;
    }

    // TBI 20231107 other cuts on simulated ...

  } // if constexpr (rs == eSim) {

  return kTRUE;

} // template <typename T> Bool_t ParticleCuts(T const& track)

//============================================================

template <eRecSim rs, typename T>
void FillParticleHistograms(T const& track, eBeforeAfter ba)
{
  // Fill all particle histograms for reconstructed and simulated data.

  // a) Fill reconstructed (and corresponding MC truth simulated);
  // b) Fill only simulated.

  if (tc.fVerboseForEachParticle) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Fill reconstructed ...:
  if constexpr (rs == eRec || rs == eRecAndSim) {
    cph_a.fParticleHistograms[ePhi][eRec][ba]->Fill(track.phi());
    cph_a.fParticleHistograms[ePt][eRec][ba]->Fill(track.pt());
    cph_a.fParticleHistograms[eEta][eRec][ba]->Fill(track.eta());

    // ... and corresponding MC truth simulated ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for this track, skip...");
        return;
      }
      auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle
      cph_a.fParticleHistograms[ePhi][eSim][ba]->Fill(mcparticle.phi());
      cph_a.fParticleHistograms[ePt][eSim][ba]->Fill(mcparticle.pt());
      cph_a.fParticleHistograms[eEta][eSim][ba]->Fill(mcparticle.eta());
    } // if constexpr (rs == eRecAndSim) {
  }   // if constexpr (rs == eRec || rs == eRecAndSim) {

  // b) Fill only simulated:
  if constexpr (rs == eSim) {
    // Remark: in this branch, 'track' is always TracksSim = aod::McParticles
    cph_a.fParticleHistograms[ePhi][eSim][ba]->Fill(track.phi());
    cph_a.fParticleHistograms[ePt][eSim][ba]->Fill(track.pt());
    cph_a.fParticleHistograms[eEta][eSim][ba]->Fill(track.eta());
  } // if constexpr (rs == eSim) {

  /* TBI 20231019 use also these + check further
  // From aod::TracksExtra
  cph_a.fParticleHistograms[etpcNClsCrossedRows][rs][ba]->Fill(track.tpcNClsCrossedRows());

  // From aod::TracksDCA
  cph_a.fParticleHistograms[eDCA_xy][rs][ba]->Fill(track.dcaXY());
  cph_a.fParticleHistograms[eDCA_z][rs][ba]->Fill(track.dcaZ());
  */

} // template <eRecSim rs, typename T> void FillParticleHistograms(...)

//============================================================

void CalculateCorrelations()
{
  // Calculate analytically multiparticle correlations from Q-vectors.
  // In this method, only isotropic correlations for which all harmonics are the
  // same are evaluated.

  // a) Flush 'n' fill the generic Q-vectors;
  // b) Calculate correlations;
  // c) Flush the generic Q-vectors.

  if (tc.fVerbose) {
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
    if (fSelectedTracks < 2) {
      return;
    }
    if (tc.fVerbose) {
      LOGF(info, "\033[1;32m%s => calculating 2-particle correlations....\033[0m", __PRETTY_FUNCTION__);
    }
    TComplex two = Two(h, -h);
    Double_t twoC = two.Re(); // cos
    // Double_t twoS = two.Im(); // sin
    Double_t wTwo = Two(0, 0).Re(); // Weight is 'number of combinations' by default TBI
                                    // 20220809 add support for other weights
    if (wTwo > 0.0) {
      twoC /= wTwo;
    } else {
      LOGF(fatal, "In function \033[1;31m%s at line %d, wTwo = %f <=0. fSelectedTracks = %d\033[0m", __PRETTY_FUNCTION__, __LINE__, wTwo, fSelectedTracks);
    }

    if (fCalculateCustomNestedLoop) {
      // e-b-e sanity check:
      TArrayI* harmonics = new TArrayI(2);
      harmonics->SetAt(h, 0);
      harmonics->SetAt(-h, 1);
      Double_t nestedLoopValue = this->CalculateCustomNestedLoop(harmonics);
      if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(twoC - nestedLoopValue) > 1.e-5) {
        LOGF(fatal, "in function \033[1;31m%s at line %d, nestedLoopValue = %f is not the same as twoC = %f\033[0m", __PRETTY_FUNCTION__, __LINE__, nestedLoopValue, twoC);
      } else {
        LOGF(info, "=> e-b-e check with CustomNestedLoop is OK for isotropic 2-p, harmonic %d", h);
      }
      delete harmonics;
      harmonics = NULL;
    } // if(fCalculateCustomNestedLoop)

    /*
      // for on-the-fly and internal validation, rescale results with theoretical value: if(fCalculateOnTheFly && fOnTheFlyFlowAmplitudes &&
      fRescaleWithTheoreticalInput &&
         TMath::Abs(fOnTheFlyFlowAmplitudes->GetAt(h-1))>0.){twoC/=pow(fOnTheFlyFlowAmplitudes->GetAt(h-1),2.);}
      else if(fUseInternalValidation && fInternalValidationAmplitudes &&
      fRescaleWithTheoreticalInput &&
              TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){twoC/=pow(fInternalValidationAmplitudes->GetAt(h-1),2.);}
    */

    // integrated:
    if (c_a.fCorrelationsPro[0][h - 1][AFO_INTEGRATED]) {
      c_a.fCorrelationsPro[0][h - 1][AFO_INTEGRATED]->Fill(0.5, twoC, wTwo);
    }
    // vs. multiplicity:
    if (c_a.fCorrelationsPro[0][h - 1][AFO_MULTIPLICITY]) {
      c_a.fCorrelationsPro[0][h - 1][AFO_MULTIPLICITY]->Fill(fSelectedTracks + 0.5, twoC, wTwo);
    }
    // vs. centrality:
    if (c_a.fCorrelationsPro[0][h - 1][AFO_CENTRALITY]) {
      c_a.fCorrelationsPro[0][h - 1][AFO_CENTRALITY]->Fill(fCentrality, twoC, wTwo);
    }

    // 4p:
    if (fSelectedTracks < 4) {
      continue;
    } // yes, continue, because I can still calculate 2-p in other harmonics!
    if (tc.fVerbose) {
      LOGF(info, "\033[1;32m%s => calculating 4-particle correlations....\033[0m", __PRETTY_FUNCTION__);
    }
    TComplex four = Four(h, h, -h, -h);
    Double_t fourC = four.Re(); // cos
    // Double_t fourS = four.Im(); // sin
    Double_t wFour = Four(0, 0, 0, 0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
    if (wFour > 0.0) {
      fourC /= wFour;
    } else {
      LOGF(fatal, "In function \033[1;31m%s at line %d, wFour = %f <=0. fSelectedTracks = %d\033[0m", __PRETTY_FUNCTION__, __LINE__, wFour, fSelectedTracks);
      // TBI 20240110 shall I 'continue' here, instead of bailing out?
    }

    if (fCalculateCustomNestedLoop) {
      // e-b-e sanity check:
      TArrayI* harmonics = new TArrayI(4);
      harmonics->SetAt(h, 0);
      harmonics->SetAt(h, 1);
      harmonics->SetAt(-h, 2);
      harmonics->SetAt(-h, 3);
      Double_t nestedLoopValue = this->CalculateCustomNestedLoop(harmonics);
      if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(fourC - nestedLoopValue) > 1.e-5) {
        LOGF(fatal, "in function \033[1;31m%s at line %d, nestedLoopValue = %f is not the same as fourC = %f\033[0m", __PRETTY_FUNCTION__, __LINE__, nestedLoopValue, fourC);
      } else {
        LOGF(info, "=> e-b-e check with CustomNestedLoop is OK for isotropic 4-p, harmonic %d", h);
      }
      delete harmonics;
      harmonics = NULL;
    } // if(fCalculateCustomNestedLoop)

    //    if(fUseInternalValidation && fInternalValidationAmplitudes && fRescaleWithTheoreticalInput &&
    //       TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){fourC/=pow(fInternalValidationAmplitudes->GetAt(h-1),4.);}

    // integrated:
    if (c_a.fCorrelationsPro[1][h - 1][AFO_INTEGRATED]) {
      c_a.fCorrelationsPro[1][h - 1][AFO_INTEGRATED]->Fill(0.5, fourC, wFour);
    }
    // vs. multiplicity:
    if (c_a.fCorrelationsPro[1][h - 1][AFO_MULTIPLICITY]) {
      c_a.fCorrelationsPro[1][h - 1][AFO_MULTIPLICITY]->Fill(fSelectedTracks + 0.5, fourC, wFour);
    }
    // vs. centrality:
    if (c_a.fCorrelationsPro[1][h - 1][AFO_CENTRALITY]) {
      c_a.fCorrelationsPro[1][h - 1][AFO_CENTRALITY]->Fill(fCentrality, fourC, wFour);
    }

    // 6p:
    if (fSelectedTracks < 6) {
      continue;
    } // yes, continue, because I can still calculate 2-p and 4-p in other harmonics!
    if (tc.fVerbose) {
      LOGF(info, "\033[1;32m%s => calculating 6-particle correlations....\033[0m", __PRETTY_FUNCTION__);
    }
    TComplex six = Six(h, h, h, -h, -h, -h);
    Double_t sixC = six.Re(); // cos
    // Double_t sixS = six.Im(); // sin
    Double_t wSix = Six(0, 0, 0, 0, 0, 0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
    if (wSix > 0.0) {
      sixC /= wSix;
    } else {
      LOGF(fatal, "In function \033[1;31m%s at line %d, wSix = %f <=0. fSelectedTracks = %d\033[0m", __PRETTY_FUNCTION__, __LINE__, wSix, fSelectedTracks);
      // TBI 20240110 shall I 'continue' here, instead of bailing out?
    }

    if (fCalculateCustomNestedLoop) {
      // e-b-e sanity check:
      TArrayI* harmonics = new TArrayI(6);
      harmonics->SetAt(h, 0);
      harmonics->SetAt(h, 1);
      harmonics->SetAt(h, 2);
      harmonics->SetAt(-h, 3);
      harmonics->SetAt(-h, 4);
      harmonics->SetAt(-h, 5);
      Double_t nestedLoopValue = this->CalculateCustomNestedLoop(harmonics);
      if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(sixC - nestedLoopValue) > 1.e-5) {
        LOGF(fatal, "in function \033[1;31m%s at line %d, nestedLoopValue = %f is not the same as sixC = %f\033[0m", __PRETTY_FUNCTION__, __LINE__, nestedLoopValue, sixC);
      } else {
        LOGF(info, "=> e-b-e check with CustomNestedLoop is OK for isotropic 6-p, harmonic %d", h);
      }
      delete harmonics;
      harmonics = NULL;
    } // if(fCalculateCustomNestedLoop)

    //    if(fUseInternalValidation && fInternalValidationAmplitudes && fRescaleWithTheoreticalInput &&
    //       TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){sixC/=pow(fInternalValidationAmplitudes->GetAt(h-1),4.);}

    // integrated:
    if (c_a.fCorrelationsPro[2][h - 1][AFO_INTEGRATED]) {
      c_a.fCorrelationsPro[2][h - 1][AFO_INTEGRATED]->Fill(0.5, sixC, wSix);
    }
    // vs. multiplicity:
    if (c_a.fCorrelationsPro[2][h - 1][AFO_MULTIPLICITY]) {
      c_a.fCorrelationsPro[2][h - 1][AFO_MULTIPLICITY]->Fill(fSelectedTracks + 0.5, sixC, wSix);
    }
    // vs. centrality:
    if (c_a.fCorrelationsPro[2][h - 1][AFO_CENTRALITY]) {
      c_a.fCorrelationsPro[2][h - 1][AFO_CENTRALITY]->Fill(fCentrality, sixC, wSix);
    }

    // 8p:
    if (fSelectedTracks < 8) {
      continue;
    } // yes, continue, because I can still calculate 2-p, 4-p and 6-p in other harmonics!
    if (tc.fVerbose) {
      LOGF(info, "\033[1;32m%s => calculating 8-particle correlations....\033[0m", __PRETTY_FUNCTION__);
    }
    TComplex eight = Eight(h, h, h, h, -h, -h, -h, -h);
    Double_t eightC = eight.Re(); // cos
    // Double_t eightS = eight.Im(); // sin
    Double_t wEight = Eight(0, 0, 0, 0, 0, 0, 0, 0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
    if (wEight > 0.0) {
      eightC /= wEight;
    } else {
      LOGF(fatal, "In function \033[1;31m%s at line %d, wEight = %f <=0. fSelectedTracks = %d\033[0m", __PRETTY_FUNCTION__, __LINE__, wEight, fSelectedTracks);
      // TBI 20240110 shall I 'continue' here, instead of bailing out?
    }

    if (fCalculateCustomNestedLoop) {
      // e-b-e sanity check:
      TArrayI* harmonics = new TArrayI(8);
      harmonics->SetAt(h, 0);
      harmonics->SetAt(h, 1);
      harmonics->SetAt(h, 2);
      harmonics->SetAt(h, 3);
      harmonics->SetAt(-h, 4);
      harmonics->SetAt(-h, 5);
      harmonics->SetAt(-h, 6);
      harmonics->SetAt(-h, 7);
      Double_t nestedLoopValue = this->CalculateCustomNestedLoop(harmonics);
      if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(eightC - nestedLoopValue) > 1.e-5) {
        LOGF(fatal, "in function \033[1;31m%s at line %d, nestedLoopValue = %f is not the same as eightC = %f\033[0m", __PRETTY_FUNCTION__, __LINE__, nestedLoopValue, eightC);
      } else {
        LOGF(info, "=> e-b-e check with CustomNestedLoop is OK for isotropic 8-p, harmonic %d", h);
      }
      delete harmonics;
      harmonics = NULL;
    } // if(fCalculateCustomNestedLoop)

    //    if(fUseInternalValidation && fInternalValidationAmplitudes && fRescaleWithTheoreticalInput &&
    //       TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){eightC/=pow(fInternalValidationAmplitudes->GetAt(h-1),4.);}

    // integrated:
    if (c_a.fCorrelationsPro[3][h - 1][AFO_INTEGRATED]) {
      c_a.fCorrelationsPro[3][h - 1][AFO_INTEGRATED]->Fill(0.5, eightC, wEight);
    }
    // vs. multiplicity:
    if (c_a.fCorrelationsPro[3][h - 1][AFO_MULTIPLICITY]) {
      c_a.fCorrelationsPro[3][h - 1][AFO_MULTIPLICITY]->Fill(fSelectedTracks + 0.5, eightC, wEight);
    }
    // vs. centrality:
    if (c_a.fCorrelationsPro[3][h - 1][AFO_CENTRALITY]) {
      c_a.fCorrelationsPro[3][h - 1][AFO_CENTRALITY]->Fill(fCentrality, eightC, wEight);
    }
  } // for(Int_t h=1;h<=gMaxHarmonic;h++) // harmonic

  // c) Flush the generic Q-vectors:
  ResetQ();

} // void CalculateCorrelations()

//============================================================

void CalculateTest0()
{
  // Calculate Test0.

  // a) Flush 'n' fill the generic Q-vectors;
  // b) Calculate correlations;
  // c) Flush the generic Q-vectors.

  if (tc.fVerbose) {
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
  Double_t correlation = 0.; // still has to be divided with 'weight' later, to get average correlation
  Double_t weight = 0.;
  Int_t n[gMaxCorrelator] = {0}; // array holding harmonics

  for (Int_t mo = 0; mo < gMaxCorrelator; mo++) {
    for (Int_t mi = 0; mi < gMaxIndex; mi++) {
      // TBI 20210913 I do not have to loop each time all the way up to gMaxCorrelator and gMaxIndex, but nevermind now, it's not a big efficiency loss.

      // Sanitize the labels (If necessary. Locally this is irrelevant):
      if (!t0_a.fTest0Labels[mo][mi]) // I do not stream them.
      {
        for (Int_t v = 0; v < eAsFunctionOf_N; v++) {
          if (t0_a.fTest0Pro[mo][mi][v]) {
            t0_a.fTest0Labels[mo][mi] = new TString(t0_a.fTest0Pro[mo][mi][v]->GetTitle()); // there is no memory leak here, since this is executed only once due to if(!fTest0Labels[mo][mi])
            break;                                                                          // yes, since for all v they are the same, so I just need to fetch it from one
          }
        }
      } // if(!t0_afTest0Labels[mo][mi])

      if (t0_a.fTest0Labels[mo][mi]) {
        // Extract harmonics from TString, FS is " ":
        for (Int_t h = 0; h <= mo; h++) {
          TObjArray* oa = t0_a.fTest0Labels[mo][mi]->Tokenize(" ");
          if (!oa) {
            LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
                 __PRETTY_FUNCTION__, __LINE__);
          }
          n[h] = TString(oa->At(h)->GetName()).Atoi();
          delete oa; // yes, otherwise it's a memory leak
        }

        switch (mo + 1) // which order? yes, mo+1
        {
          case 1:
            if (fSelectedTracks < 1) {
              return;
            }
            correlation = One(n[0]).Re();
            weight = One(0).Re();
            break;

          case 2:
            if (fSelectedTracks < 2) {
              return;
            }
            correlation = Two(n[0], n[1]).Re();
            weight = Two(0, 0).Re();
            break;

          case 3:
            if (fSelectedTracks < 3) {
              return;
            }
            correlation = Three(n[0], n[1], n[2]).Re();
            weight = Three(0, 0, 0).Re();
            break;

          case 4:
            if (fSelectedTracks < 4) {
              return;
            }
            correlation = Four(n[0], n[1], n[2], n[3]).Re();
            weight = Four(0, 0, 0, 0).Re();
            break;

          case 5:
            if (fSelectedTracks < 5) {
              return;
            }
            correlation = Five(n[0], n[1], n[2], n[3], n[4]).Re();
            weight = Five(0, 0, 0, 0, 0).Re();
            break;

          case 6:
            if (fSelectedTracks < 6) {
              return;
            }
            correlation = Six(n[0], n[1], n[2], n[3], n[4], n[5]).Re();
            weight = Six(0, 0, 0, 0, 0, 0).Re();
            break;

          case 7:
            if (fSelectedTracks < 7) {
              return;
            }
            correlation = Seven(n[0], n[1], n[2], n[3], n[4], n[5], n[6]).Re();
            weight = Seven(0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 8:
            if (fSelectedTracks < 8) {
              return;
            }
            correlation = Eight(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7]).Re();
            weight = Eight(0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 9:
            if (fSelectedTracks < 9) {
              return;
            }
            correlation = Nine(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8]).Re();
            weight = Nine(0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 10:
            if (fSelectedTracks < 10) {
              return;
            }
            correlation = Ten(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9]).Re();
            weight = Ten(0, 0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 11:
            if (fSelectedTracks < 11) {
              return;
            }
            correlation = Eleven(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10]).Re();
            weight = Eleven(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 12:
            if (fSelectedTracks < 12) {
              return;
            }
            correlation = Twelve(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11]).Re();
            weight = Twelve(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          default:
            LOGF(fatal, "in function \033[1;31m%s at line %d\n Not supported yet: %s \n\n\033[0m", __PRETTY_FUNCTION__, __LINE__, t0_a.fTest0Labels[mo][mi]->Data());
        } // switch(mo+1)

        // Insanity check on weight:
        if (!(weight > 0.)) {
          LOGF(fatal, "in function \033[1;31m%s at line %d\n Is perhaps order of correlator bigger than the number of particles?: %s \n\n\033[0m", __PRETTY_FUNCTION__, __LINE__, t0_a.fTest0Labels[mo][mi]->Data());
        }

        // e-b-e sanity check:
        if (fCalculateCustomNestedLoop) {
          TArrayI* harmonics = new TArrayI(mo + 1);
          for (Int_t i = 0; i < mo + 1; i++) {
            harmonics->SetAt(n[i], i);
          }
          Double_t nestedLoopValue = this->CalculateCustomNestedLoop(harmonics);
          if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(correlation / weight - nestedLoopValue) > 1.e-5) {
            LOGF(fatal, "in function \033[1;31m%s at line %d, nestedLoopValue = %f is not the same as correlation/weight = %f, for correlator %s\033[0m", __PRETTY_FUNCTION__, __LINE__, nestedLoopValue, correlation / weight, t0_a.fTest0Labels[mo][mi]->Data());
          } else {
            LOGF(info, "=> e-b-e check with CustomNestedLoop is OK for %d-p Test0 corr. %s", mo + 1, t0_a.fTest0Labels[mo][mi]->Data());
          }
          delete harmonics;
          harmonics = NULL;
        } // if(fCalculateCustomNestedLoop)

        /*
            // To ease comparison, rescale with theoretical value. Now all Test0 results shall be at 1:
            if(fUseInternalValidation && fInternalValidationAmplitudes && fInternalValidationPlanes && fRescaleWithTheoreticalInput)
            {
             TArrayI *harmonics = new TArrayI(mo+1);
             for(Int_t i=0;i<mo+1;i++)
             {
              harmonics->SetAt(n[i],i);
             }
             TComplex theoreticalValue = TheoreticalValue(harmonics,fInternalValidationAmplitudes,fInternalValidationPlanes);
             if(TMath::Abs(theoreticalValue.Re()) > 0.)
             {
              correlation /= theoreticalValue.Re();
             }
             delete harmonics; harmonics = NULL;
            } // if(fUseInternalValidation && fRescaleWithTheoreticalInput)
        */

        // Finally, fill:
        // integrated:
        if (t0_a.fTest0Pro[mo][mi][AFO_INTEGRATED]) {
          t0_a.fTest0Pro[mo][mi][AFO_INTEGRATED]->Fill(0.5, correlation / weight, weight);
        }
        // vs. multiplicity:
        if (t0_a.fTest0Pro[mo][mi][AFO_MULTIPLICITY]) {
          t0_a.fTest0Pro[mo][mi][AFO_MULTIPLICITY]->Fill(fSelectedTracks + 0.5, correlation / weight, weight);
        }
        // vs. centrality:
        if (t0_a.fTest0Pro[mo][mi][AFO_CENTRALITY]) {
          t0_a.fTest0Pro[mo][mi][AFO_CENTRALITY]->Fill(fCentrality, correlation / weight, weight);
        }
      } // if(fTest0Labels[mo][mi])
    }   // for(Int_t mi=0;mi<gMaxIndex;mi++)
  }     // for(Int_t mo=0;mo<gMaxCorrelator;mo++)

  // c) Flush the generic Q-vectors:
  ResetQ();

} // void CalculateTest0()

//============================================================

void CalculateNestedLoops()
{
  // Calculate correlations with nested loops.

  // a) 2-particle nested loops;
  // b) 4-particle nested loops;
  // c) 6-particle nested loops;
  // d) 8-particle nested loops.

  if (tc.fVerbose) {
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

  // a) 2-particle nested loops:
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
      for (int h = 0; h < gMaxHarmonic; h++) {
        // fill cos, 2p, integreated:
        if (nl_a.fNestedLoopsPro[0][h][AFO_INTEGRATED]) {
          nl_a.fNestedLoopsPro[0][h][AFO_INTEGRATED]->Fill(
            0.5, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
        // fill cos, 2p, vs. multiplicity:
        if (nl_a.fNestedLoopsPro[0][h][AFO_MULTIPLICITY]) {
          nl_a.fNestedLoopsPro[0][h][AFO_MULTIPLICITY]->Fill(
            fSelectedTracks + 0.5, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)),
            dW1 * dW2);
        }
        // fill cos, 2p, vs. centrality:
        if (nl_a.fNestedLoopsPro[0][h][AFO_CENTRALITY]) {
          nl_a.fNestedLoopsPro[0][h][AFO_CENTRALITY]->Fill(
            fCentrality, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
      } // for(int h=1; h<=6; h++)
    }   // for(int i2=0; i2<nParticles; i2++)
  }     // for(int i1=0; i1<nParticles; i1++)

  // b) 4-particle nested loops:
  if (nParticles < 4) {
    return;
  }
  LOGF(info, "\033[1;32m       CalculateNestedLoops(void), 4-p correlations .... \033[0m");
  for (int i1 = 0; i1 < nParticles; i1++) {
    Double_t dPhi1 = nl_a.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl_a.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      Double_t dPhi2 = nl_a.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl_a.ftaNestedLoops[1]->GetAt(i2);
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        Double_t dPhi3 = nl_a.ftaNestedLoops[0]->GetAt(i3);
        Double_t dW3 = nl_a.ftaNestedLoops[1]->GetAt(i3);
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          Double_t dPhi4 = nl_a.ftaNestedLoops[0]->GetAt(i4);
          Double_t dW4 = nl_a.ftaNestedLoops[1]->GetAt(i4);
          for (int h = 0; h < gMaxHarmonic; h++) {
            // fill cos, 4p, integreated:
            if (nl_a.fNestedLoopsPro[1][h][AFO_INTEGRATED]) {
              nl_a.fNestedLoopsPro[1][h][AFO_INTEGRATED]->Fill(0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. M:
            if (nl_a.fNestedLoopsPro[1][h][AFO_MULTIPLICITY]) {
              nl_a.fNestedLoopsPro[1][h][AFO_MULTIPLICITY]->Fill(fSelectedTracks + 0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. centrality:
            if (nl_a.fNestedLoopsPro[1][h][AFO_CENTRALITY]) {
              nl_a.fNestedLoopsPro[1][h][AFO_CENTRALITY]->Fill(fCentrality, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
          } // for(int h=0; h<gMaxHarmonic; h++)
        }   // for(int i4=0; i4<nParticles; i4++)
      }     // for(int i3=0; i3<nParticles; i3++)
    }       // for(int i2=0; i2<nTracks; i2++)
  }         // for(int i1=0; i1<nTracks; i1++)

  // c) 6-particle nested loops:
  if (nParticles < 6) {
    return;
  }
  LOGF(info, "\033[1;32m       CalculateNestedLoops(void), 6-p correlations .... \033[0m");
  for (int i1 = 0; i1 < nParticles; i1++) {
    Double_t dPhi1 = nl_a.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl_a.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      Double_t dPhi2 = nl_a.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl_a.ftaNestedLoops[1]->GetAt(i2);
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        Double_t dPhi3 = nl_a.ftaNestedLoops[0]->GetAt(i3);
        Double_t dW3 = nl_a.ftaNestedLoops[1]->GetAt(i3);
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          Double_t dPhi4 = nl_a.ftaNestedLoops[0]->GetAt(i4);
          Double_t dW4 = nl_a.ftaNestedLoops[1]->GetAt(i4);
          for (int i5 = 0; i5 < nParticles; i5++) {
            if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4) {
              continue;
            }
            Double_t dPhi5 = nl_a.ftaNestedLoops[0]->GetAt(i5);
            Double_t dW5 = nl_a.ftaNestedLoops[1]->GetAt(i5);
            for (int i6 = 0; i6 < nParticles; i6++) {
              if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5) {
                continue;
              }
              Double_t dPhi6 = nl_a.ftaNestedLoops[0]->GetAt(i6);
              Double_t dW6 = nl_a.ftaNestedLoops[1]->GetAt(i6);
              for (int h = 0; h < gMaxHarmonic; h++) {
                // fill cos, 6p, integreated:
                if (nl_a.fNestedLoopsPro[2][h][AFO_INTEGRATED]) {
                  nl_a.fNestedLoopsPro[2][h][AFO_INTEGRATED]->Fill(0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. M:
                if (nl_a.fNestedLoopsPro[2][h][AFO_MULTIPLICITY]) {
                  nl_a.fNestedLoopsPro[2][h][AFO_MULTIPLICITY]->Fill(fSelectedTracks + 0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. M:
                if (nl_a.fNestedLoopsPro[2][h][AFO_CENTRALITY]) {
                  nl_a.fNestedLoopsPro[2][h][AFO_CENTRALITY]->Fill(fCentrality, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
              } // for(int h=0; h<gMaxHarmonic; h++)
            }   // if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
          }     // if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
        }       // for(int i4=0; i4<nParticles; i4++)
      }         // for(int i3=0; i3<nParticles; i3++)
    }           // for(int i2=0; i2<nTracks; i2++)
  }             // for(int i1=0; i1<nTracks; i1++)

  // d) 8-particle nested loops:
  if (nParticles < 8) {
    return;
  }
  LOGF(info, "\033[1;32m       CalculateNestedLoops(void), 8-p correlations .... \033[0m");
  for (int i1 = 0; i1 < nParticles; i1++) {
    Double_t dPhi1 = nl_a.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl_a.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      Double_t dPhi2 = nl_a.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl_a.ftaNestedLoops[1]->GetAt(i2);
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        Double_t dPhi3 = nl_a.ftaNestedLoops[0]->GetAt(i3);
        Double_t dW3 = nl_a.ftaNestedLoops[1]->GetAt(i3);
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          Double_t dPhi4 = nl_a.ftaNestedLoops[0]->GetAt(i4);
          Double_t dW4 = nl_a.ftaNestedLoops[1]->GetAt(i4);
          for (int i5 = 0; i5 < nParticles; i5++) {
            if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4) {
              continue;
            }
            Double_t dPhi5 = nl_a.ftaNestedLoops[0]->GetAt(i5);
            Double_t dW5 = nl_a.ftaNestedLoops[1]->GetAt(i5);
            for (int i6 = 0; i6 < nParticles; i6++) {
              if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5) {
                continue;
              }
              Double_t dPhi6 = nl_a.ftaNestedLoops[0]->GetAt(i6);
              Double_t dW6 = nl_a.ftaNestedLoops[1]->GetAt(i6);
              for (int i7 = 0; i7 < nParticles; i7++) {
                if (i7 == i1 || i7 == i2 || i7 == i3 || i7 == i4 || i7 == i5 || i7 == i6) {
                  continue;
                }
                Double_t dPhi7 = nl_a.ftaNestedLoops[0]->GetAt(i7);
                Double_t dW7 = nl_a.ftaNestedLoops[1]->GetAt(i7);
                for (int i8 = 0; i8 < nParticles; i8++) {
                  if (i8 == i1 || i8 == i2 || i8 == i3 || i8 == i4 || i8 == i5 || i8 == i6 || i8 == i7) {
                    continue;
                  }
                  Double_t dPhi8 = nl_a.ftaNestedLoops[0]->GetAt(i8);
                  Double_t dW8 = nl_a.ftaNestedLoops[1]->GetAt(i8);
                  for (int h = 0; h < gMaxHarmonic; h++) {
                    // fill cos, 8p, integreated:
                    if (nl_a.fNestedLoopsPro[3][h][AFO_INTEGRATED]) {
                      nl_a.fNestedLoopsPro[3][h][AFO_INTEGRATED]->Fill(0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. M:
                    if (nl_a.fNestedLoopsPro[3][h][AFO_MULTIPLICITY]) {
                      nl_a.fNestedLoopsPro[3][h][AFO_MULTIPLICITY]->Fill(fSelectedTracks + 0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. M:
                    if (nl_a.fNestedLoopsPro[3][h][AFO_CENTRALITY]) {
                      nl_a.fNestedLoopsPro[3][h][AFO_CENTRALITY]->Fill(fCentrality, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                  } // for(int h=0; h<gMaxHarmonic; h++)
                }   // for(int i8=0; i8<nParticles; i8++)
              }     // for(int i7=0; i7<nParticles; i7++)
            }       // for(int i6=0; i6<nParticles; i6++)
          }         // for(int i5=0; i5<nParticles; i6++)
        }           // for(int i4=0; i4<nParticles; i4++)
      }             // for(int i3=0; i3<nParticles; i3++)
    }               // for(int i2=0; i2<nParticles; i2++)
  }                 // for(int i1=0; i1<nParticles; i1++)

} // void CalculateNestedLoops()

//============================================================

void ComparisonNestedLoopsVsCorrelations()
{
  // Compare analytic results from Q-vectors and brute force results from nested loops.
  // Use only for small multiplicities, when nested loops are still feasible.
  // Results have to be exactly the same in each case.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  Int_t nBinsQV = -44;
  Int_t nBinsNL = -44;
  Double_t valueQV = 0.;
  Double_t valueNL = 0.;

  for (Int_t v = 0; v < 3; v++) { // TBI 20240116 this corresponds to the ordering of variables in enum eAsFunctionOf . Here (for the time being) I compare only int, mult. and cent.
    // a) Integrated comparison:
    nBinsQV = c_a.fCorrelationsPro[0][0][v]->GetNbinsX();
    nBinsNL = nl_a.fNestedLoopsPro[0][0][v]->GetNbinsX();
    if (nBinsQV != nBinsNL) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    LOGF(info, "\033[1;32m   [%d] : %s\033[0m", v, rh_a.fResultsHistogramsXaxisTitle[v].Data());
    for (Int_t o = 0; o < 4; o++) {
      LOGF(info, "\033[1;32m   ==== <<%d>>-particle correlations ====\033[0m", 2 * (o + 1));
      for (Int_t h = 0; h < gMaxHarmonic; h++) {
        for (Int_t b = 1; b <= nBinsQV; b++) {
          if (c_a.fCorrelationsPro[o][h][v]) {
            valueQV = c_a.fCorrelationsPro[o][h][v]->GetBinContent(b);
          }
          if (nl_a.fNestedLoopsPro[o][h][v]) {
            valueNL = nl_a.fNestedLoopsPro[o][h][v]->GetBinContent(b);
          }
          if (TMath::Abs(valueQV) > 0. && TMath::Abs(valueNL) > 0.) {
            LOGF(info, "   bin=%d, h=%d, Q-vectors:    %f", b, h + 1, valueQV);
            LOGF(info, "   bin=%d, h=%d, Nested loops: %f", b, h + 1, valueNL);
            if (TMath::Abs(valueQV - valueNL) > 1.e-5) {
              LOGF(info, "\n\033[1;33m[%d][%d][%d] \033[0m\n", o, h, v);
              LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
                   __PRETTY_FUNCTION__, __LINE__);
            }
          }           // if(TMath::Abs(valueQV)>0. && TMath::Abs(valueNL)>0.)
        }             // for(Int_t b=1;b<=nBinsQV;b++)
      }               // for(Int_t h=0;h<6;h++)
      LOGF(info, ""); // new line
    }                 // for(Int_t o=0;o<4;o++)
  }                   // for (Int_t v = 0; v < 3; v++)

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

TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5)
{
  // Generic five-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5)]>.

  TComplex five = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) - Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) - Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) - Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5, 1) - Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) + Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5, 1) - Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) + Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5, 1) + 2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5, 1) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) + Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5, 1) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5, 1) + 2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5, 1) - 6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5, 1) - Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) + Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5, 2) + Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5, 2) + Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5, 2) - 2. * Q(n2 + n3 + n4, 3) * Q(n1 + n5, 2) - Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) + Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5, 2) + Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5, 2) + Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5, 2) - 2. * Q(n1 + n3 + n4, 3) * Q(n2 + n5, 2) + 2. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5, 3) - 2. * Q(n3 + n4, 2) * Q(n1 + n2 + n5, 3) - Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) + Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5, 2) + Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5, 2) + Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5, 2) - 2. * Q(n1 + n2 + n4, 3) * Q(n3 + n5, 2) + 2. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5, 3) - 2. * Q(n2 + n4, 2) * Q(n1 + n3 + n5, 3) + 2. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5, 3) - 2. * Q(n1 + n4, 2) * Q(n2 + n3 + n5, 3) - 6. * Q(n4, 1) * Q(n1 + n2 + n3 + n5, 4) - Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) + Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5, 2) + Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5, 2) + Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5, 2) - 2. * Q(n1 + n2 + n3, 3) * Q(n4 + n5, 2) + 2. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5, 3) - 2. * Q(n2 + n3, 2) * Q(n1 + n4 + n5, 3) + 2. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5, 3) - 2. * Q(n1 + n3, 2) * Q(n2 + n4 + n5, 3) - 6. * Q(n3, 1) * Q(n1 + n2 + n4 + n5, 4) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5, 3) - 2. * Q(n1 + n2, 2) * Q(n3 + n4 + n5, 3) - 6. * Q(n2, 1) * Q(n1 + n3 + n4 + n5, 4) - 6. * Q(n1, 1) * Q(n2 + n3 + n4 + n5, 4) + 24. * Q(n1 + n2 + n3 + n4 + n5, 5);

  return five;

} // TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5)

//============================================================

TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6)
{
  // Generic six-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6)]>.

  TComplex six = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) - Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) - Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) - Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) - Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n6, 1) + Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n6, 1) - Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n6, 1) + Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n6, 1) + 2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5, 1) * Q(n6, 1) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n6, 1) + Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n6, 1) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5, 1) * Q(n6, 1) + 2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5, 1) * Q(n6, 1) - 6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5, 1) * Q(n6, 1) - Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n6, 1) + Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n6, 1) + Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5, 2) * Q(n6, 1) + Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5, 2) * Q(n6, 1) - 2. * Q(n2 + n3 + n4, 3) * Q(n1 + n5, 2) * Q(n6, 1) - Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n6, 1) + Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n6, 1) + Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5, 2) * Q(n6, 1) + Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5, 2) * Q(n6, 1) - 2. * Q(n1 + n3 + n4, 3) * Q(n2 + n5, 2) * Q(n6, 1) + 2. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5, 3) * Q(n6, 1) - 2. * Q(n3 + n4, 2) * Q(n1 + n2 + n5, 3) * Q(n6, 1) - Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n6, 1) + Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n6, 1) + Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5, 2) * Q(n6, 1) + Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5, 2) * Q(n6, 1) - 2. * Q(n1 + n2 + n4, 3) * Q(n3 + n5, 2) * Q(n6, 1) + 2. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5, 3) * Q(n6, 1) - 2. * Q(n2 + n4, 2) * Q(n1 + n3 + n5, 3) * Q(n6, 1) + 2. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5, 3) * Q(n6, 1) - 2. * Q(n1 + n4, 2) * Q(n2 + n3 + n5, 3) * Q(n6, 1) - 6. * Q(n4, 1) * Q(n1 + n2 + n3 + n5, 4) * Q(n6, 1) - Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n6, 1) + Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n6, 1) + Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5, 2) * Q(n6, 1) + Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5, 2) * Q(n6, 1) - 2. * Q(n1 + n2 + n3, 3) * Q(n4 + n5, 2) * Q(n6, 1) + 2. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5, 3) * Q(n6, 1) - 2. * Q(n2 + n3, 2) * Q(n1 + n4 + n5, 3) * Q(n6, 1) + 2. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5, 3) * Q(n6, 1) - 2. * Q(n1 + n3, 2) * Q(n2 + n4 + n5, 3) * Q(n6, 1) - 6. * Q(n3, 1) * Q(n1 + n2 + n4 + n5, 4) * Q(n6, 1) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5, 3) * Q(n6, 1) - 2. * Q(n1 + n2, 2) * Q(n3 + n4 + n5, 3) * Q(n6, 1) - 6. * Q(n2, 1) * Q(n1 + n3 + n4 + n5, 4) * Q(n6, 1) - 6. * Q(n1, 1) * Q(n2 + n3 + n4 + n5, 4) * Q(n6, 1) + 24. * Q(n1 + n2 + n3 + n4 + n5, 5) * Q(n6, 1) - Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n6, 2) + Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n6, 2) + Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n1 + n6, 2) + Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n1 + n6, 2) - 2. * Q(n2 + n3 + n4, 3) * Q(n5, 1) * Q(n1 + n6, 2) + Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n1 + n6, 2) - Q(n3 + n4, 2) * Q(n2 + n5, 2) * Q(n1 + n6, 2) + Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n1 + n6, 2) - Q(n2 + n4, 2) * Q(n3 + n5, 2) * Q(n1 + n6, 2) - 2. * Q(n4, 1) * Q(n2 + n3 + n5, 3) * Q(n1 + n6, 2) + Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n1 + n6, 2) - Q(n2 + n3, 2) * Q(n4 + n5, 2) * Q(n1 + n6, 2) - 2. * Q(n3, 1) * Q(n2 + n4 + n5, 3) * Q(n1 + n6, 2) - 2. * Q(n2, 1) * Q(n3 + n4 + n5, 3) * Q(n1 + n6, 2) + 6. * Q(n2 + n3 + n4 + n5, 4) * Q(n1 + n6, 2) - Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n6, 2) + Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n6, 2) + Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n2 + n6, 2) + Q(n1, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n2 + n6, 2) - 2. * Q(n1 + n3 + n4, 3) * Q(n5, 1) * Q(n2 + n6, 2) + Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n2 + n6, 2) - Q(n3 + n4, 2) * Q(n1 + n5, 2) * Q(n2 + n6, 2) + Q(n1, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n2 + n6, 2) - Q(n1 + n4, 2) * Q(n3 + n5, 2) * Q(n2 + n6, 2) - 2. * Q(n4, 1) * Q(n1 + n3 + n5, 3) * Q(n2 + n6, 2) + Q(n1, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n2 + n6, 2) - Q(n1 + n3, 2) * Q(n4 + n5, 2) * Q(n2 + n6, 2) - 2. * Q(n3, 1) * Q(n1 + n4 + n5, 3) * Q(n2 + n6, 2) - 2. * Q(n1, 1) * Q(n3 + n4 + n5, 3) * Q(n2 + n6, 2) + 6. * Q(n1 + n3 + n4 + n5, 4) * Q(n2 + n6, 2) + 2. * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n2 + n6, 3) - 2. * Q(n3 + n4, 2) * Q(n5, 1) * Q(n1 + n2 + n6, 3) - 2. * Q(n4, 1) * Q(n3 + n5, 2) * Q(n1 + n2 + n6, 3) - 2. * Q(n3, 1) * Q(n4 + n5, 2) * Q(n1 + n2 + n6, 3) + 4. * Q(n3 + n4 + n5, 3) * Q(n1 + n2 + n6, 3) - Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n5, 1) * Q(n3 + n6, 2) + Q(n1 + n2, 2) * Q(n4, 1) * Q(n5, 1) * Q(n3 + n6, 2) + Q(n2, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n3 + n6, 2) + Q(n1, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n3 + n6, 2) - 2. * Q(n1 + n2 + n4, 3) * Q(n5, 1) * Q(n3 + n6, 2) + Q(n2, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n3 + n6, 2) - Q(n2 + n4, 2) * Q(n1 + n5, 2) * Q(n3 + n6, 2) + Q(n1, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n3 + n6, 2) - Q(n1 + n4, 2) * Q(n2 + n5, 2) * Q(n3 + n6, 2) - 2. * Q(n4, 1) * Q(n1 + n2 + n5, 3) * Q(n3 + n6, 2) + Q(n1, 1) * Q(n2, 1) * Q(n4 + n5, 2) * Q(n3 + n6, 2) - Q(n1 + n2, 2) * Q(n4 + n5, 2) * Q(n3 + n6, 2) - 2. * Q(n2, 1) * Q(n1 + n4 + n5, 3) * Q(n3 + n6, 2) - 2. * Q(n1, 1) * Q(n2 + n4 + n5, 3) * Q(n3 + n6, 2) + 6. * Q(n1 + n2 + n4 + n5, 4) * Q(n3 + n6, 2) + 2. * Q(n2, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n3 + n6, 3) - 2. * Q(n2 + n4, 2) * Q(n5, 1) * Q(n1 + n3 + n6, 3) - 2. * Q(n4, 1) * Q(n2 + n5, 2) * Q(n1 + n3 + n6, 3) - 2. * Q(n2, 1) * Q(n4 + n5, 2) * Q(n1 + n3 + n6, 3) + 4. * Q(n2 + n4 + n5, 3) * Q(n1 + n3 + n6, 3) + 2. * Q(n1, 1) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n3 + n6, 3) - 2. * Q(n1 + n4, 2) * Q(n5, 1) * Q(n2 + n3 + n6, 3) - 2. * Q(n4, 1) * Q(n1 + n5, 2) * Q(n2 + n3 + n6, 3) - 2. * Q(n1, 1) * Q(n4 + n5, 2) * Q(n2 + n3 + n6, 3) + 4. * Q(n1 + n4 + n5, 3) * Q(n2 + n3 + n6, 3) - 6. * Q(n4, 1) * Q(n5, 1) * Q(n1 + n2 + n3 + n6, 4) + 6. * Q(n4 + n5, 2) * Q(n1 + n2 + n3 + n6, 4) - Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n5, 1) * Q(n4 + n6, 2) + Q(n1 + n2, 2) * Q(n3, 1) * Q(n5, 1) * Q(n4 + n6, 2) + Q(n2, 1) * Q(n1 + n3, 2) * Q(n5, 1) * Q(n4 + n6, 2) + Q(n1, 1) * Q(n2 + n3, 2) * Q(n5, 1) * Q(n4 + n6, 2) - 2. * Q(n1 + n2 + n3, 3) * Q(n5, 1) * Q(n4 + n6, 2) + Q(n2, 1) * Q(n3, 1) * Q(n1 + n5, 2) * Q(n4 + n6, 2) - Q(n2 + n3, 2) * Q(n1 + n5, 2) * Q(n4 + n6, 2) + Q(n1, 1) * Q(n3, 1) * Q(n2 + n5, 2) * Q(n4 + n6, 2) - Q(n1 + n3, 2) * Q(n2 + n5, 2) * Q(n4 + n6, 2) - 2. * Q(n3, 1) * Q(n1 + n2 + n5, 3) * Q(n4 + n6, 2) + Q(n1, 1) * Q(n2, 1) * Q(n3 + n5, 2) * Q(n4 + n6, 2) - Q(n1 + n2, 2) * Q(n3 + n5, 2) * Q(n4 + n6, 2) - 2. * Q(n2, 1) * Q(n1 + n3 + n5, 3) * Q(n4 + n6, 2) - 2. * Q(n1, 1) * Q(n2 + n3 + n5, 3) * Q(n4 + n6, 2) + 6. * Q(n1 + n2 + n3 + n5, 4) * Q(n4 + n6, 2) + 2. * Q(n2, 1) * Q(n3, 1) * Q(n5, 1) * Q(n1 + n4 + n6, 3) - 2. * Q(n2 + n3, 2) * Q(n5, 1) * Q(n1 + n4 + n6, 3) - 2. * Q(n3, 1) * Q(n2 + n5, 2) * Q(n1 + n4 + n6, 3) - 2. * Q(n2, 1) * Q(n3 + n5, 2) * Q(n1 + n4 + n6, 3) + 4. * Q(n2 + n3 + n5, 3) * Q(n1 + n4 + n6, 3) + 2. * Q(n1, 1) * Q(n3, 1) * Q(n5, 1) * Q(n2 + n4 + n6, 3) - 2. * Q(n1 + n3, 2) * Q(n5, 1) * Q(n2 + n4 + n6, 3) - 2. * Q(n3, 1) * Q(n1 + n5, 2) * Q(n2 + n4 + n6, 3) - 2. * Q(n1, 1) * Q(n3 + n5, 2) * Q(n2 + n4 + n6, 3) + 4. * Q(n1 + n3 + n5, 3) * Q(n2 + n4 + n6, 3) - 6. * Q(n3, 1) * Q(n5, 1) * Q(n1 + n2 + n4 + n6, 4) + 6. * Q(n3 + n5, 2) * Q(n1 + n2 + n4 + n6, 4) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n5, 1) * Q(n3 + n4 + n6, 3) - 2. * Q(n1 + n2, 2) * Q(n5, 1) * Q(n3 + n4 + n6, 3) - 2. * Q(n2, 1) * Q(n1 + n5, 2) * Q(n3 + n4 + n6, 3) - 2. * Q(n1, 1) * Q(n2 + n5, 2) * Q(n3 + n4 + n6, 3) + 4. * Q(n1 + n2 + n5, 3) * Q(n3 + n4 + n6, 3) - 6. * Q(n2, 1) * Q(n5, 1) * Q(n1 + n3 + n4 + n6, 4) + 6. * Q(n2 + n5, 2) * Q(n1 + n3 + n4 + n6, 4) - 6. * Q(n1, 1) * Q(n5, 1) * Q(n2 + n3 + n4 + n6, 4) + 6. * Q(n1 + n5, 2) * Q(n2 + n3 + n4 + n6, 4) + 24. * Q(n5, 1) * Q(n1 + n2 + n3 + n4 + n6, 5) - Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5 + n6, 2) + Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5 + n6, 2) + Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5 + n6, 2) + Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5 + n6, 2) - 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5 + n6, 2) + Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5 + n6, 2) - Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5 + n6, 2) + Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5 + n6, 2) - Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5 + n6, 2) - 2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5 + n6, 2) + Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5 + n6, 2) - Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5 + n6, 2) - 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5 + n6, 2) - 2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5 + n6, 2) + 6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5 + n6, 2) + 2. * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5 + n6, 3) - 2. * Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5 + n6, 3) - 2. * Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5 + n6, 3) - 2. * Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5 + n6, 3) + 4. * Q(n2 + n3 + n4, 3) * Q(n1 + n5 + n6, 3) + 2. * Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5 + n6, 3) - 2. * Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5 + n6, 3) - 2. * Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5 + n6, 3) - 2. * Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5 + n6, 3) + 4. * Q(n1 + n3 + n4, 3) * Q(n2 + n5 + n6, 3) - 6. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5 + n6, 4) + 6. * Q(n3 + n4, 2) * Q(n1 + n2 + n5 + n6, 4) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5 + n6, 3) - 2. * Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5 + n6, 3) - 2. * Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5 + n6, 3) - 2. * Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5 + n6, 3) + 4. * Q(n1 + n2 + n4, 3) * Q(n3 + n5 + n6, 3) - 6. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5 + n6, 4) + 6. * Q(n2 + n4, 2) * Q(n1 + n3 + n5 + n6, 4) - 6. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5 + n6, 4) + 6. * Q(n1 + n4, 2) * Q(n2 + n3 + n5 + n6, 4) + 24. * Q(n4, 1) * Q(n1 + n2 + n3 + n5 + n6, 5) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5 + n6, 3) - 2. * Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5 + n6, 3) - 2. * Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5 + n6, 3) - 2. * Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5 + n6, 3) + 4. * Q(n1 + n2 + n3, 3) * Q(n4 + n5 + n6, 3) - 6. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5 + n6, 4) + 6. * Q(n2 + n3, 2) * Q(n1 + n4 + n5 + n6, 4) - 6. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5 + n6, 4) + 6. * Q(n1 + n3, 2) * Q(n2 + n4 + n5 + n6, 4) + 24. * Q(n3, 1) * Q(n1 + n2 + n4 + n5 + n6, 5) - 6. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5 + n6, 4) + 6. * Q(n1 + n2, 2) * Q(n3 + n4 + n5 + n6, 4) + 24. * Q(n2, 1) * Q(n1 + n3 + n4 + n5 + n6, 5) + 24. * Q(n1, 1) * Q(n2 + n3 + n4 + n5 + n6, 5) - 120. * Q(n1 + n2 + n3 + n4 + n5 + n6, 6);

  return six;

} // TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6)

//============================================================

TComplex Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)
{
  // Generic seven-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7)]>.

  Int_t harmonic[7] = {n1, n2, n3, n4, n5, n6, n7};

  TComplex seven = Recursion(7, harmonic);

  return seven;

} // end of TComplex Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)

//============================================================

TComplex Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)
{
  // Generic eight-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8)]>.

  Int_t harmonic[8] = {n1, n2, n3, n4, n5, n6, n7, n8};

  TComplex eight = Recursion(8, harmonic);

  return eight;

} // end of Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)

//============================================================

TComplex Nine(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9)
{
  // Generic nine-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9)]>.

  Int_t harmonic[9] = {n1, n2, n3, n4, n5, n6, n7, n8, n9};

  TComplex nine = Recursion(9, harmonic);

  return nine;

} // end of TComplex Nine(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9)

//============================================================

TComplex Ten(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10)
{
  // Generic ten-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9+n10*phi10)]>.

  Int_t harmonic[10] = {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10};

  TComplex ten = Recursion(10, harmonic);

  return ten;

} // end of TComplex Ten(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10)

//============================================================

TComplex Eleven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11)
{
  // Generic eleven-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9+n10*phi10+n11*phi11)]>.

  Int_t harmonic[11] = {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11};

  TComplex eleven = Recursion(11, harmonic);

  return eleven;

} // end of TComplex Eleven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11)

//============================================================

TComplex Twelve(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11, Int_t n12)
{
  // Generic twelve-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9+n10*phi10+n11*phi11+n12*phi12)]>.

  Int_t harmonic[12] = {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12};

  TComplex twelve = Recursion(12, harmonic);

  return twelve;

} // end of TComplex Twelve(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11, Int_t n12)

//============================================================

TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0)
{
  // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by
  // Kristjan Gulbrandsen (gulbrand@nbi.dk).

  Int_t nm1 = n - 1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0)
    return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip)
    return c;

  Int_t multp1 = mult + 1;
  Int_t nm2 = n - 2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n - 3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1)
    return c - c2;
  return c - Double_t(mult) * c2;

} // TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0)

//============================================================

void ResetQ()
{
  // Reset the components of generic Q-vectors. Use it whenever you call the
  // standard functions for correlations, for some custom Q-vectors.

  if (tc.fVerbose) {
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

  if (tc.fVerbose) {
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

  if (tc.fVerbose) {
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

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
    LOGF(info, "\033[1;33m filePath = %s\033[0m", filePath);
    LOGF(info, "\033[1;33m runNumber = %s\033[0m", runNumber);
    LOGF(info, "\033[1;33m variable = %s\033[0m", variable);
    LOGF(info, "\033[1;33m fTaskName = %s\033[0m", tc.fTaskName.Data());
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
  //    dir AliEn. If filePath begins with "/alice-ccdb.cern.ch/" then it's in
  //    CCDB. Therefore, files in AliEn and CCDB must be specified with abs path,
  //    for local files both abs and relative paths are just fine.
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
      listWithRuns, Form("%s_%s", variable, tc.fTaskName.Data())));
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
    if (tc.fVerbose) {
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
      listWithRuns, Form("%s_%s", variable, tc.fTaskName.Data())));

    if (!hist) {
      LOGF(info, "\033[1;33m%s_%s \033[0m", variable, tc.fTaskName.Data());
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
      listWithRuns, Form("%s_%s", variable, tc.fTaskName.Data())));
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

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Return value:
  TObjArray* oa = NULL;

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

    // Fetch TObjArray from external file (keep in sync with local file case below):
    TList* lok = oaFile->GetListOfKeys();
    if (!lok) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    for (Int_t l = 0; l < lok->GetEntries(); l++) {
      oaFile->GetObject(lok->At(l)->GetName(), oa);
      if (oa && TString(oa->ClassName()).EqualTo("TObjArray")) {
        break; // TBI 20231107 the working assumption is that in an external file there is only one TObjArray object,
               // and here I fetch it, whatever its name is. The advantage is that I do not have to do
               // any additional work for TObjArray's name. Since I do not anticipate ever having more than 1
               // TObjArray in an external file, this shall be alright. With the current implementation,
               // if there are multiple TObjArray objects in the same ROOT file, the first one will be fetched.
      }
    } // for(Int_t l=0;l<lok->GetEntries();l++)

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

    // Fetch TObjArray from external file (keep in sync with AliEn file case above):
    TList* lok = oaFile->GetListOfKeys();
    if (!lok) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    for (Int_t l = 0; l < lok->GetEntries(); l++) {
      oaFile->GetObject(lok->At(l)->GetName(), oa);
      if (oa && TString(oa->ClassName()).EqualTo("TObjArray")) {
        break; // TBI 20231107 the working assumption is that in an external file there is only one TObjArray object,
               // and here I fetch it, whatever its name is. The advantage is that I do not have to do
               // any additional work for TObjArray's name. Since I do not anticipate ever having more than 1
               // TObjArray in an external file, this shall be alright. With the current implementation,
               // if there are multiple TObjArray objects in the same ROOT file, the first one will be fetched.
      }
    } // for(Int_t l=0;l<lok->GetEntries();l++)
    if (!oa) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }

  } // else {

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s => Fetched TObjArray named \"%s\" from file %s\033[0m", __PRETTY_FUNCTION__, oa->GetName(), filePath);
  }

  return oa;

} // TObjArray* GetObjArrayWithLabels(const char *filePath)

//============================================================

void StoreLabelsInPlaceholder()
{
  // Storal all Test0 labels in the temporary placeholder.

  // a) Initialize all counters;
  // b) Fetch TObjArray with labels from an external file;
  // c) Book the placeholder fTest0LabelsPlaceholder for all labels;
  // d) Finally, store the labels from external source into placeholder;
  // e) Insantity check on labels.

  if (tc.fVerbose) {
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

  // e) Insantity check on labels:
  //    Here I am merely checking that harmonic larget than gMaxHarmonic was not requested.
  for (Int_t b = 1; b <= fTest0LabelsPlaceholder->GetXaxis()->GetNbins(); b++) {
    TObjArray* temp = TString(fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b)).Tokenize(" ");
    for (Int_t h = 0; h < temp->GetEntries(); h++) {
      if (TMath::Abs(TString(temp->At(h)->GetName()).Atoi()) > gMaxHarmonic) {
        LOGF(info, "\033[1;31m bin = %d, label = %s, gMaxHarmonic = %d\033[0m", b, fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b), (Int_t)gMaxHarmonic);
        LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
      }          // if(TString(temp->At(h)->GetName()).Atoi() > gMaxHarmonic) {
    }            // for(Int_t h = 0; h < temp->GetEntries(); h++) {
    delete temp; // yes, otherwise it's a memory leak
  }              // for(Int_t b = 1; b <= fTest0LabelsPlaceholder->GetXaxis()->GetNbins(); b++) {

} // void StoreLabelsInPlaceholder()

//============================================================

Bool_t RetrieveCorrelationsLabels()
{
  // Generate the labels of all correlations of interest, i.e. retrieve them
  // from TH1I *fTest0LabelsPlaceholder

  if (tc.fVerbose) {
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

void GetParticleWeights()
{
  // Get the particle weights. Call this function only once.

  //    TBI 20231012 Here the current working assumption is that:
  //    a) Corrections do not change within a given run;
  //    b) Hyperloop proceeses the dataset one masterjob per run number.
  //    If any of these 2 assumptions are violated, this code will have to be modified.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  if (pw_a.fUseWeights[wPHI]) {
    TH1D* phiWeights = GetHistogramWithWeights(fFileWithWeights.Data(), tc.fRunNumber.Data(), "phi");
    if (!phiWeights) {
      LOGF(fatal, "in function \033[1;31m%s at line %d, phiWeights is NULL. Check the external file %s with particle weights\033[0m", __PRETTY_FUNCTION__, __LINE__, fFileWithWeights.Data());
    }
    SetWeightsHist(phiWeights, "phi");
  }

  if (pw_a.fUseWeights[wPT]) {
    TH1D* ptWeights = GetHistogramWithWeights(fFileWithWeights.Data(), tc.fRunNumber.Data(), "pt");
    if (!ptWeights) {
      LOGF(fatal, "in function \033[1;31m%s at line %d, ptWeights is NULL. Check the external file %s with particle weights\033[0m", __PRETTY_FUNCTION__, __LINE__, fFileWithWeights.Data());
    }
    SetWeightsHist(ptWeights, "pt");
  }

  if (pw_a.fUseWeights[wETA]) {
    TH1D* etaWeights = GetHistogramWithWeights(fFileWithWeights.Data(), tc.fRunNumber.Data(), "eta");
    if (!etaWeights) {
      LOGF(fatal, "in function \033[1;31m%s at line %d, etaWeights is NULL. Check the external file %s with particle weights\033[0m", __PRETTY_FUNCTION__, __LINE__, fFileWithWeights.Data());
    }
    SetWeightsHist(etaWeights, "eta");
  }

} // void GetParticleWeights()

//============================================================

Bool_t MaxNumberOfEvents()
{
  // Check if max number of events was reached. See also configurable cNumberOfEvents_max.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // *) Return value:
  Bool_t reachedMaxNumberOfEvents = kFALSE;

  // *) Determine from which histogram the relevant info will be taken:
  Int_t rs = -44; // reconstructed or simulated
  if (gProcessRec || gProcessRecSim) {
    rs = eRec;
  } else if (gProcessSim) {
    rs = eSim;
  } else {
    LOGF(fatal, "in function \033[1;31m%s at line %d, not a single flag gProcess* is true \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  // *) Okay, do the thing:
  if (ceh_a.fEventHistograms[eNumberOfEvents][rs][eAfter]->GetBinContent(1) >= ceh_a.fEventCuts[eNumberOfEvents][eMax]) {
    reachedMaxNumberOfEvents = kTRUE;
  }

  // *) Hasta la vista:
  return reachedMaxNumberOfEvents;

} // void MaxNumberOfEvents()

//============================================================

Double_t CalculateCustomNestedLoop(TArrayI* harmonics)
{
  // For the specified harmonics, get the correlation from nested loops.
  // Order of correlator is the number of harmonics, i.e. the number of elements in an array.

  // a) Determine the order of correlator;
  // b) Custom nested loop;
  // c) Return value.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  if (!harmonics) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  Int_t nParticles = fSelectedTracks;
  /* TBI 20231108 enable eventually
  if(fUseFixedNumberOfRandomlySelectedParticles)
  {
   nParticles = 0;
   for(Int_t i=0;i<nl_a.ftaNestedLoops[0]->GetSize();i++)
   {
    if(TMath::Abs(nl_a.ftaNestedLoops[0]->GetAt(i)) > 0. && TMath::Abs(nl_a.ftaNestedLoops[1]->GetAt(i)) > 0.){nParticles++;}
   }
  }
  */

  // a) Determine the order of correlator;
  Int_t order = harmonics->GetSize();
  if (0 == order || order > gMaxCorrelator) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  // b) Custom nested loop:
  TProfile* profile = new TProfile("profile", "", 1, 0., 1.); // helper profile to get all averages automatically
  // profile->Sumw2();
  Double_t value = 0.;  // cos of current multiplet
  Double_t weight = 1.; // weight of current multiplet
  for (int i1 = 0; i1 < nParticles; i1++) {
    Double_t dPhi1 = nl_a.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl_a.ftaNestedLoops[1]->GetAt(i1);
    if (1 == order) {
      value = TMath::Cos(harmonics->GetAt(0) * dPhi1);
      weight = dW1;
      profile->Fill(0.5, value, weight);
      continue;
    }
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      Double_t dPhi2 = nl_a.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl_a.ftaNestedLoops[1]->GetAt(i2);
      if (2 == order) {
        value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2);
        weight = dW1 * dW2;
        profile->Fill(0.5, value, weight);
        continue;
      }
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        Double_t dPhi3 = nl_a.ftaNestedLoops[0]->GetAt(i3);
        Double_t dW3 = nl_a.ftaNestedLoops[1]->GetAt(i3);
        if (3 == order) {
          value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3);
          weight = dW1 * dW2 * dW3;
          profile->Fill(0.5, value, weight);
          continue;
        }
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          Double_t dPhi4 = nl_a.ftaNestedLoops[0]->GetAt(i4);
          Double_t dW4 = nl_a.ftaNestedLoops[1]->GetAt(i4);
          if (4 == order) {
            value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4);
            weight = dW1 * dW2 * dW3 * dW4;
            profile->Fill(0.5, value, weight);
            continue;
          }
          for (int i5 = 0; i5 < nParticles; i5++) {
            if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4) {
              continue;
            }
            Double_t dPhi5 = nl_a.ftaNestedLoops[0]->GetAt(i5);
            Double_t dW5 = nl_a.ftaNestedLoops[1]->GetAt(i5);
            if (5 == order) {
              value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5);
              weight = dW1 * dW2 * dW3 * dW4 * dW5;
              profile->Fill(0.5, value, weight);
              continue;
            }
            for (int i6 = 0; i6 < nParticles; i6++) {
              if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5) {
                continue;
              }
              Double_t dPhi6 = nl_a.ftaNestedLoops[0]->GetAt(i6);
              Double_t dW6 = nl_a.ftaNestedLoops[1]->GetAt(i6);
              if (6 == order) {
                value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5 + harmonics->GetAt(5) * dPhi6);
                weight = dW1 * dW2 * dW3 * dW4 * dW5 * dW6;
                profile->Fill(0.5, value, weight);
                continue;
              }
              for (int i7 = 0; i7 < nParticles; i7++) {
                if (i7 == i1 || i7 == i2 || i7 == i3 || i7 == i4 || i7 == i5 || i7 == i6) {
                  continue;
                }
                Double_t dPhi7 = nl_a.ftaNestedLoops[0]->GetAt(i7);
                Double_t dW7 = nl_a.ftaNestedLoops[1]->GetAt(i7);
                if (7 == order) {
                  value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5 + harmonics->GetAt(5) * dPhi6 + harmonics->GetAt(6) * dPhi7);
                  weight = dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7;
                  profile->Fill(0.5, value, weight);
                  continue;
                }
                for (int i8 = 0; i8 < nParticles; i8++) {
                  if (i8 == i1 || i8 == i2 || i8 == i3 || i8 == i4 || i8 == i5 || i8 == i6 || i8 == i7) {
                    continue;
                  }
                  Double_t dPhi8 = nl_a.ftaNestedLoops[0]->GetAt(i8);
                  Double_t dW8 = nl_a.ftaNestedLoops[1]->GetAt(i8);
                  if (8 == order) {
                    value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5 + harmonics->GetAt(5) * dPhi6 + harmonics->GetAt(6) * dPhi7 + harmonics->GetAt(7) * dPhi8);
                    weight = dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8;
                    profile->Fill(0.5, value, weight);
                    continue;
                  }
                  for (int i9 = 0; i9 < nParticles; i9++) {
                    if (i9 == i1 || i9 == i2 || i9 == i3 || i9 == i4 || i9 == i5 || i9 == i6 || i9 == i7 || i9 == i8) {
                      continue;
                    }
                    Double_t dPhi9 = nl_a.ftaNestedLoops[0]->GetAt(i9);
                    Double_t dW9 = nl_a.ftaNestedLoops[1]->GetAt(i9);
                    if (9 == order) {
                      value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5 + harmonics->GetAt(5) * dPhi6 + harmonics->GetAt(6) * dPhi7 + harmonics->GetAt(7) * dPhi8 + harmonics->GetAt(8) * dPhi9);
                      weight = dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8 * dW9;
                      profile->Fill(0.5, value, weight);
                      continue;
                    }
                    for (int i10 = 0; i10 < nParticles; i10++) {
                      if (i10 == i1 || i10 == i2 || i10 == i3 || i10 == i4 || i10 == i5 || i10 == i6 || i10 == i7 || i10 == i8 || i10 == i9) {
                        continue;
                      }
                      Double_t dPhi10 = nl_a.ftaNestedLoops[0]->GetAt(i10);
                      Double_t dW10 = nl_a.ftaNestedLoops[1]->GetAt(i10);
                      if (10 == order) {
                        value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5 + harmonics->GetAt(5) * dPhi6 + harmonics->GetAt(6) * dPhi7 + harmonics->GetAt(7) * dPhi8 + harmonics->GetAt(8) * dPhi9 + harmonics->GetAt(9) * dPhi10);
                        weight = dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8 * dW9 * dW10;
                        profile->Fill(0.5, value, weight);
                        continue;
                      }
                      for (int i11 = 0; i11 < nParticles; i11++) {
                        if (i11 == i1 || i11 == i2 || i11 == i3 || i11 == i4 || i11 == i5 || i11 == i6 || i11 == i7 || i11 == i8 || i11 == i9 || i11 == i10) {
                          continue;
                        }
                        Double_t dPhi11 = nl_a.ftaNestedLoops[0]->GetAt(i11);
                        Double_t dW11 = nl_a.ftaNestedLoops[1]->GetAt(i11);
                        if (11 == order) {
                          value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5 + harmonics->GetAt(5) * dPhi6 + harmonics->GetAt(6) * dPhi7 + harmonics->GetAt(7) * dPhi8 + harmonics->GetAt(8) * dPhi9 + harmonics->GetAt(9) * dPhi10 + harmonics->GetAt(10) * dPhi11);
                          weight = dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8 * dW9 * dW10 * dW11;
                          profile->Fill(0.5, value, weight);
                          continue;
                        }
                        for (int i12 = 0; i12 < nParticles; i12++) {
                          if (i12 == i1 || i12 == i2 || i12 == i3 || i12 == i4 || i12 == i5 || i12 == i6 || i12 == i7 || i12 == i8 || i12 == i9 || i12 == i10 || i12 == i11) {
                            continue;
                          }
                          Double_t dPhi12 = nl_a.ftaNestedLoops[0]->GetAt(i12);
                          Double_t dW12 = nl_a.ftaNestedLoops[1]->GetAt(i12);
                          if (12 == order) {
                            value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5 + harmonics->GetAt(5) * dPhi6 + harmonics->GetAt(6) * dPhi7 + harmonics->GetAt(7) * dPhi8 + harmonics->GetAt(8) * dPhi9 + harmonics->GetAt(9) * dPhi10 + harmonics->GetAt(10) * dPhi11 + harmonics->GetAt(11) * dPhi12);
                            weight = dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8 * dW9 * dW10 * dW11 * dW12;
                            profile->Fill(0.5, value, weight);
                            continue;
                          }

                          // ... it's easy to continue the above pattern here

                        } // for(int i12=0; i12<nParticles; i12++)
                      }   // for(int i11=0; i11<nParticles; i11++)
                    }     // for(int i10=0; i10<nParticles; i10++)
                  }       // for(int i9=0; i9<nParticles; i9++)
                }         // for(int i8=0; i8<nParticles; i8++)
              }           // for(int i7=0; i7<nParticles; i7++)
            }             // for(int i6=0; i6<nParticles; i6++)
          }               // for(int i5=0; i5<nParticles; i5++)
        }                 // for(int i4=0; i4<nParticles; i4++)
      }                   // for(int i3=0; i3<nParticles; i3++)
    }                     // for(int i2=0; i2<nParticles; i2++)
  }                       // for(int i1=0; i1<nParticles; i1++)

  // c) Return value:
  Double_t finalValue = profile->GetBinContent(1);
  delete profile;
  profile = NULL;
  return finalValue;

} // Double_t CalculateCustomNestedLoop(TArrayI *harmonics)

//============================================================

template <eRecSim rs, typename T>
void DetermineCentrality(T const& collision)
{
  // Determine collision centrality.

  // a) For real data, determine centrality from default centrality estimator;
  // b) For simulated data, determine centrality directly from impact parameter.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // a) For real data, determine centrality from default centrality estimator:
  if constexpr (rs == eRec || rs == eRecAndSim) {
    // fCentrality = gRandom->Uniform(0.,100.);  // collision.centFT0M(); // TBI 20240120 not ready yet, estimators are specific for Run 1,2,3 data processing ...
    fCentrality = collision.centFT0M(); // TBI 20240120 not ready yet, estimators are specific for Run 1,2,3 data processing ...
    // TBI 20240120 I could also here access also corresponding simulated centrality from impact parameter, if available through collision.has_mcCollision()
  }

  // b) For simulated data, determine centrality directly from impact parameter:
  if constexpr (rs == eSim) {
    fCentrality = -44.; // TBI 20240120 add support eventualy
  }                     // if constexpr (rs == eSim) {

  // TBI 20240120 remove this printout eventually:
  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m fCentrality = %f\033[0m", fCentrality);
  }

} // template <eRecSim rs, typename T> void DetermineCentrality(T const& collision)

//============================================================

void RandomIndices(Int_t nTracks)
{
  // Randomize indices using Fisher-Yates algorithm.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  if (nTracks < 1) {
    return;
  }

  // Fisher-Yates algorithm:
  tc.fRandomIndices = new TArrayI(nTracks);
  tc.fRandomIndices->Reset(); // just in case there is some random garbage in memory at init
  for (Int_t i = 0; i < nTracks; i++) {
    tc.fRandomIndices->AddAt(i, i);
  }
  for (Int_t i = nTracks - 1; i >= 1; i--) {
    Int_t j = gRandom->Integer(i + 1);
    Int_t temp = tc.fRandomIndices->GetAt(j);
    tc.fRandomIndices->AddAt(tc.fRandomIndices->GetAt(i), j);
    tc.fRandomIndices->AddAt(temp, i);
  } // end of for(Int_t i=nTracks-1;i>=1;i--)

} // void RandomIndices(Int_t nTracks)

//============================================================

void CalculateEverything()
{
  // Calculate everything for selected events and particles.
  // Remark: Data members for Q-vectors, containers for nested loops, etc., must all be filled when this function is called.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // *) Progress info:
  LOGF(info, "\033[1;32m=> Processing event %d/%d (selected/total) .... \033[0m", (Int_t)ceh_a.fEventHistograms[eNumberOfEvents][eRec][eAfter]->GetBinContent(1), (Int_t)ceh_a.fEventHistograms[eNumberOfEvents][eRec][eBefore]->GetBinContent(1));

  // *) Calculate multiparticle correlations (standard, isotropic, same harmonic):
  if (fCalculateCorrelations) {
    this->CalculateCorrelations();
  }

  // *) Calculate Test0: TBI 20240110 name convention
  if (fCalculateTest0) {
    this->CalculateTest0();
  }

  // *) Calculate nested loops:
  if (fCalculateNestedLoops) {
    this->CalculateNestedLoops();
    this->ComparisonNestedLoopsVsCorrelations(); // I call it here, so comparison is performed cumulatively after each event. The final printout corresponds to all events.
  }

} // void CalculateEverything()

//============================================================

template <eRecSim rs, typename T1, typename T2>
void Steer(T1 const& collision, T2 const& tracks)
{
  // This is the only function to be called in processRec(...), processRecSim(...), and processSim(...).
  // All analysis workflow is defined step-by-step here, via dedicated function calls.
  // The order of function calls obviously matters.

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // *) Do all thingies before starting to process data from this collision (e.g. count number of events, fetch the run number, etc.):
  Preprocess(collision);

  // *) Determine collision centrality:
  DetermineCentrality<rs>(collision);

  // *) Fill event histograms before event cuts:
  FillEventHistograms<rs>(collision, tracks, eBefore);

  // *) Event cuts:
  if (!EventCuts<rs>(collision, tracks)) {
    return;
  }

  // *) Main loop over particles:
  MainLoopOverParticles<rs>(tracks);

  // *) Remaining event cuts which can be applied only after the loop over particles is performed:
  if ((fSelectedTracks < ceh_a.fEventCuts[eSelectedTracks][eMin]) || (fSelectedTracks > ceh_a.fEventCuts[eSelectedTracks][eMax])) {
    if (tc.fVerbose) {
      LOGF(info, "\033[1;31m%s eSelectedTracks\033[0m", __FUNCTION__); // just a bare function name
    }
    ResetEventByEventQuantities();
    return;
  }

  // *) Fill event histograms after event AND particle cuts: // TBI 20240110 not sure still if this one is called here, or it has to be moved above
  FillEventHistograms<rs>(collision, tracks, eAfter);

  // *) Calculate everything for selected events and particles:
  CalculateEverything();

  // *) Reset event-by-event quantities:
  ResetEventByEventQuantities();

} // template <eRecSim rs, typename T1, typename T2> void Steer(T1 const* collision, T2 const* tracks)

//============================================================

template <eRecSim rs, typename T>
void MainLoopOverParticles(T const& tracks)
{
  // This is the main loop over particles, in which Q-vectors and particle histograms are filled, particle cuts applied, etc.

  // Remark:
  // To process only 'rec', set gProcessRec = true via configurable "cfWhatToProcess".
  // To process both 'rec' and 'sim', set gProcessRecSim = true via configurable "cfWhatToProcess".
  // To process only 'sim', set gProcessSim = true via configurable "cfWhatToProcess".

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // *) Local kinematic variables and corresponding particle weights:
  Double_t dPhi = 0., wPhi = 1.; // azimuthal angle and corresponding phi weight
  Double_t dPt = 0., wPt = 1.;   // transverse momentum and corresponding pt weight
  Double_t dEta = 0., wEta = 1.; // pseudorapidity and corresponding eta weight
  Double_t wToPowerP = 1.;       // weight raised to power p
  fSelectedTracks = 0;           // reset number of selected tracks

  // *) If random access of tracks from collection is requested, use Fisher-Yates algorithm to generate random indices:
  if (tc.fUseFisherYates) {
    if (tc.fRandomIndices) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }
    this->RandomIndices(tracks.size());
    if (!tc.fRandomIndices) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }
  }

  // *) Main loop over particles:
  // for (auto& track : tracks) { // default standard way of looping of tracks
  auto track = tracks.iteratorAt(0); // set the type and scope from one instance
  for (int64_t i = 0; i < tracks.size(); i++) {

    // *) Access track sequentially from collection of tracks (default), or randomly using Fisher-Yates algorithm:
    if (!tc.fUseFisherYates) {
      track = tracks.iteratorAt(i);
    } else {
      track = tracks.iteratorAt((int64_t)tc.fRandomIndices->GetAt(i));
    }

    // *) Fill particle histograms before particle cuts:
    FillParticleHistograms<rs>(track, eBefore);

    // *) Particle cuts:
    if (!ParticleCuts<rs>(track)) {
      continue;
    }

    // *) Fill particle histograms after particle cuts:
    FillParticleHistograms<rs>(track, eAfter);

    // *) Fill Q-vectors:
    //  Kinematics (Remark: for "eRecSim" processing, kinematics is taken from reconstructed):
    dPhi = track.phi();
    dPt = track.pt();
    dEta = track.eta();

    // Particle weights:
    if (pw_a.fUseWeights[wPHI]) {
      wPhi = Weight(dPhi, "phi"); // corresponding phi weight
      if (!(wPhi > 0.)) {
        LOGF(error, "\033[1;33m%s wPhi is not positive, skipping this particle for the time being...\033[0m", __PRETTY_FUNCTION__);
        LOGF(error, "dPhi = %f\nwPhi = %f", dPhi, wPhi);
        continue;
      }
    } // if(pw_a.fUseWeights[wPHI])
    if (pw_a.fUseWeights[wPT]) {
      wPt = Weight(dPt, "pt"); // corresponding pt weight
      if (!(wPt > 0.)) {
        LOGF(error, "\033[1;33m%s wPt is not positive, skipping this particle for the time being...\033[0m", __PRETTY_FUNCTION__);
        LOGF(error, "dPt = %f\nwPt = %f", dPt, wPt);
        continue;
      }
    } // if(pw_a.fUseWeights[wPT])
    if (pw_a.fUseWeights[wETA]) {
      wEta = Weight(dEta, "eta"); // corresponding eta weight
      if (!(wEta > 0.)) {
        LOGF(error, "\033[1;33m%s wEta is not positive, skipping this particle for the time being...\033[0m", __PRETTY_FUNCTION__);
        LOGF(error, "dEta = %f\nwEta = %f", dEta, wEta);
        continue;
      }
    } // if(pw_a.fUseWeights[wETA])

    for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
      for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
        if (pw_a.fUseWeights[wPHI] || pw_a.fUseWeights[wPT] || pw_a.fUseWeights[wETA]) {
          wToPowerP = pow(wPhi * wPt * wEta, wp);
        }
        qv_a.fQvector[h][wp] += TComplex(wToPowerP * TMath::Cos(h * dPhi), wToPowerP * TMath::Sin(h * dPhi));
      } // for(Int_t wp=0;wp<gMaxCorrelator+1;wp++)
    }   // for(Int_t h=0;h<gMaxHarmonic*gMaxCorrelator+1;h++)

    // *) Nested loops containers:
    if (fCalculateNestedLoops || fCalculateCustomNestedLoop) {
      if (nl_a.ftaNestedLoops[0]) {
        nl_a.ftaNestedLoops[0]->AddAt(dPhi, fSelectedTracks);
      } // remember that the 2nd argument here must start from 0
      if (nl_a.ftaNestedLoops[1]) {
        nl_a.ftaNestedLoops[1]->AddAt(wPhi * wPt * wEta, fSelectedTracks);
      } // remember that the 2nd argument here must start from 0
    }   // if(fCalculateNestedLoops || fCalculateCustomNestedLoop)

    // *) Counter of selected tracks in the current event:
    fSelectedTracks++;
    if (fSelectedTracks >= cSelectedTracks_max) {
      break;
    }

    // *) Break the loop if fixed number of particles is taken randomly from each event (use always in combination with tc.fUseFisherYates = kTRUE):
    if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && tc.fFixedNumberOfRandomlySelectedTracks == fSelectedTracks) {
      LOGF(info, "\033[1;32mBreaking the loop over particles, since requested fixed number of %d particles was reached\033[0m", tc.fFixedNumberOfRandomlySelectedTracks);
      break;
    }

  } // for (auto& track : tracks)

  // *) Insanity check on fixed number of randomly selected tracks:
  if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && tc.fFixedNumberOfRandomlySelectedTracks < fSelectedTracks) {
    LOGF(fatal, "\033[1;31mIn this event there are too few particles (fSelectedTracks = %d), and requested number of fixed number randomly selected tracks %d couldn't be reached\033[0m", fSelectedTracks, tc.fFixedNumberOfRandomlySelectedTracks);
  }

} // template <eRecSim rs, typename T> void MainLoopOverParticles(T const& tracks) {

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_MEMBERFUNCTIONS_H_
