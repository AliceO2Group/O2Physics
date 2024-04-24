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

// ...
#include <vector>

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

  fBasePro->GetXaxis()->SetBinLabel(eDryRun, "fDryRun");
  fBasePro->Fill(+eDryRun - 0.5, (Int_t)tc.fDryRun);

  fBasePro->GetXaxis()->SetBinLabel(eVerbose, "fVerbose");
  fBasePro->Fill(+eVerbose - 0.5, (Int_t)tc.fVerbose);

  fBasePro->GetXaxis()->SetBinLabel(eVerboseForEachParticle, "fVerboseForEachParticle");
  fBasePro->Fill(+eVerboseForEachParticle - 0.5, (Int_t)tc.fVerboseForEachParticle);

  fBasePro->GetXaxis()->SetBinLabel(eDoAdditionalInsanityChecks, "fDoAdditionalInsanityChecks");
  fBasePro->Fill(+eDoAdditionalInsanityChecks - 0.5, (Int_t)tc.fDoAdditionalInsanityChecks);

  fBasePro->GetXaxis()->SetBinLabel(eUseCCDB, "fUseCCDB");
  fBasePro->Fill(+eUseCCDB - 0.5, (Int_t)tc.fUseCCDB);

  fBasePro->GetXaxis()->SetBinLabel(eWhichProcess,
                                    Form("WhichProcess = %s", tc.fWhichProcess.Data()));

  fBasePro->GetXaxis()->SetBinLabel(eRandomSeed, "fRandomSeed");
  fBasePro->Fill(+eRandomSeed - 0.5, (Int_t)tc.fRandomSeed);

  fBasePro->GetXaxis()->SetBinLabel(eUseFisherYates, "fUseFisherYates");
  fBasePro->Fill(+eUseFisherYates - 0.5, (Int_t)tc.fUseFisherYates);

  fBasePro->GetXaxis()->SetBinLabel(eFixedNumberOfRandomlySelectedTracks, "fFixedNumberOfRandomlySelectedTracks");
  fBasePro->Fill(+eFixedNumberOfRandomlySelectedTracks - 0.5, (Int_t)tc.fFixedNumberOfRandomlySelectedTracks);

  fBasePro->GetXaxis()->SetBinLabel(eUseStopwatch, "fUseStopwatch");
  fBasePro->Fill(+eUseStopwatch - 0.5, (Int_t)tc.fUseStopwatch);

  fBaseList->Add(fBasePro);

} // void BookBaseList()

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

  // d) There are also implicit variables like "doprocessSomeProcessName" within a PROCESS_SWITCH clause. For them, I do not need an entry in Configurables

  // Configurable<string> cfTaskName{ ... }
  tc.fTaskName = TString(cfTaskName);

  // Configurable<bool> cfDryRun{ ... }
  tc.fDryRun = cfDryRun;

  // Configurable<bool> cfVerbose{ ... }
  tc.fVerbose = cfVerbose;

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // yes, here
  }

  // Configurable<bool> cfVerboseForEachParticle{ ... }
  tc.fVerboseForEachParticle = cfVerboseForEachParticle;

  // Configurable<bool> cfDoAdditionalInsanityChecks{ ... }
  tc.fDoAdditionalInsanityChecks = cfDoAdditionalInsanityChecks;

  // Configurable<bool> cfUseCCDB{ ... }
  tc.fUseCCDB = cfUseCCDB;

  // Set automatically what to process, from an implicit variable "doprocessSomeProcessName" within a PROCESS_SWITCH clause:
  // Remark: As of 20240224, I have abandoned Configurable<string> cfWhatToProcess{ ... ), which is now obsolete
  tc.fProcessRec = doprocessRec;
  tc.fProcessRecSim = doprocessRecSim;
  tc.fProcessSim = doprocessSim;
  tc.fProcessRec_Run2 = doprocessRec_Run2;
  tc.fProcessRecSim_Run2 = doprocessRecSim_Run2;
  tc.fProcessSim_Run2 = doprocessSim_Run2;
  tc.fProcessRec_Run1 = doprocessRec_Run1;
  tc.fProcessRecSim_Run1 = doprocessRecSim_Run1;
  tc.fProcessSim_Run1 = doprocessSim_Run1;
  tc.fProcessTest = doprocessTest;

  // Temporarary bailout protection against cases which are not implemented/validated yet:
  if (tc.fProcessSim) {
    LOGF(fatal, "in function \033[1;31m%s at line %d - processSim(...) is not implemented/validated yet \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }
  if (tc.fProcessSim_Run2) {
    LOGF(fatal, "in function \033[1;31m%s at line %d - processSim_Run2(...) is not implemented/validated yet \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }
  if (tc.fProcessRecSim_Run1) {
    LOGF(fatal, "in function \033[1;31m%s at line %d - processRecSim_Run1(...) is not implemented/validated yet \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }
  if (tc.fProcessSim_Run1) {
    LOGF(fatal, "in function \033[1;31m%s at line %d - processSim_Run1(...) is not implemented/validated yet \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  // Set automatically generic flags, from above individual flags:
  tc.fGenericRec = tc.fProcessRec || tc.fProcessRec_Run2 || tc.fProcessRec_Run1 || tc.fProcessTest;
  tc.fGenericRecSim = tc.fProcessRecSim || tc.fProcessRecSim_Run2 || tc.fProcessRecSim_Run1;
  tc.fGenericSim = tc.fProcessSim || tc.fProcessSim_Run2 || tc.fProcessSim_Run1;

  // Set automatically tc.fWhichProcess from above individual flags:
  if (tc.fProcessRec) {
    tc.fWhichProcess = "ProcessRec";
  } else if (tc.fProcessRecSim) {
    tc.fWhichProcess = "ProcessRecSim";
  } else if (tc.fProcessSim) {
    tc.fWhichProcess = "ProcessSim";
  } else if (tc.fProcessRec_Run2) {
    tc.fWhichProcess = "ProcessRec_Run2";
  } else if (tc.fProcessRecSim_Run2) {
    tc.fWhichProcess = "ProcessRecSim_Run2";
  } else if (tc.fProcessSim_Run2) {
    tc.fWhichProcess = "ProcessSim_Run2";
  } else if (tc.fProcessRec_Run1) {
    tc.fWhichProcess = "ProcessRec_Run1";
  } else if (tc.fProcessRecSim_Run1) {
    tc.fWhichProcess = "ProcessRecSim_Run1";
  } else if (tc.fProcessSim_Run1) {
    tc.fWhichProcess = "ProcessSim_Run1";
  } else if (tc.fProcessTest) {
    tc.fWhichProcess = "ProcessTest";
  }

  // Configurable<unsigned int> cfRandomSeed{ ... )
  tc.fRandomSeed = cfRandomSeed;

  // Configurable<bool> cfUseFisherYates{ ... }
  tc.fUseFisherYates = cfUseFisherYates;

  // Configurable<int> cfFixedNumberOfRandomlySelectedTracks{ ... }
  tc.fFixedNumberOfRandomlySelectedTracks = cfFixedNumberOfRandomlySelectedTracks;

  // Configurable<bool> cfUseStopwatch{ ... }
  tc.fUseStopwatch = cfUseStopwatch;

  // *) Event cuts:
  // Configurable<string> cfTrigger{ ... )
  ec.fTrigger = TString(cfTrigger);
  // Configurable<bool> cfUseSel7{ ... }
  ec.fUseSel7 = cfUseSel7;
  // Configurable<bool> cfUseSel8{ ... }
  ec.fUseSel8 = cfUseSel8;
  // Configurable<string> cfCentralityEstimator{ ... )
  ec.fCentralityEstimator = cfCentralityEstimator;
  // ...

  // Configurable<bool> cfCalculateQvectors{ ... }
  qv.fCalculateQvectors = cfCalculateQvectors;

  // Configurable<bool> cfCalculateCorrelations{ ... };
  mupa.fCalculateCorrelations = cfCalculateCorrelations;

  // ...

  // Configurable<bool> cfCalculateTest0{ ... };
  t0.fCalculateTest0 = cfCalculateTest0; // + see below, how it's automatically set via other Test0 flags

  // Configurable<bool> cfCalculateTest0AsFunctionOfIntegrated{ ... } + analogous configurables for other ones:
  t0.fCalculateTest0AsFunctionOf[AFO_INTEGRATED] = cfCalculateTest0AsFunctionOfIntegrated;
  t0.fCalculateTest0AsFunctionOf[AFO_MULTIPLICITY] = cfCalculateTest0AsFunctionOfMultiplicity;
  t0.fCalculateTest0AsFunctionOf[AFO_CENTRALITY] = cfCalculateTest0AsFunctionOfCentrality;
  t0.fCalculateTest0AsFunctionOf[AFO_PT] = cfCalculateTest0AsFunctionOfPt;
  t0.fCalculateTest0AsFunctionOf[AFO_ETA] = cfCalculateTest0AsFunctionOfEta;
  // Use above Test0 flags to automatically set the main Test0 flag:
  for (Int_t v = 0; v < eAsFunctionOf_N; v++) {
    if (t0.fCalculateTest0AsFunctionOf[v]) {
      t0.fCalculateTest0 = true;
      break; // yes, it suffices one diff. flag to be true, for the main Test0 flag to be true
    }
  }

  // Configurable<string> cfFileWithLabels{ ... }
  t0.fFileWithLabels = TString(cfFileWithLabels);

  // Configurable<bool> cfUsePhiWeights{ ... };
  pw.fUseWeights[wPHI] = cfUsePhiWeights;

  // Configurable<bool> cfUsePtWeights{ ... };
  pw.fUseWeights[wPT] = cfUsePtWeights;

  // Configurable<bool> cfUseEtaWeights{ ... };
  pw.fUseWeights[wETA] = cfUseEtaWeights;

  // Configurable<bool> cfUseDiffPhiPtWeights{ ... };
  pw.fUseDiffWeights[wPHIPT] = cfUseDiffPhiPtWeights;

  // Configurable<bool> cfUseDiffPhiEtaWeights{ ... };
  pw.fUseDiffWeights[wPHIETA] = cfUseDiffPhiEtaWeights;

  // Configurable<string> cfFileWithWeights{ ... }
  pw.fFileWithWeights = TString(cfFileWithWeights);

  // ...

  // *) Nested loops:
  // Configurable<string> cfCalculateNestedLoops{ ... }
  nl.fCalculateNestedLoops = cfCalculateNestedLoops;

  // Configurable<string> cfCalculateCustomNestedLoops{ ... }
  nl.fCalculateCustomNestedLoop = cfCalculateCustomNestedLoops;

  // *) Results histograms:
  res.fSaveResultsHistograms = cfSaveResultsHistograms;

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
  // By default all event histograms are booked. If you do not want particular event histogram to be booked,
  // use configurable array cfBookEventHistograms, where you can specify flags 1 (book) or 0 (do not book).
  // Ordering of the flags in that array is interpreted through ordering of enums in enum eEventHistograms. // TBI 20240124 is this safe enough?
  auto lBookEventHistograms = (vector<int>)cfBookEventHistograms; // this is now the local version of that int array from configurable.
  if (lBookEventHistograms.size() != eEventHistograms_N) {
    LOGF(fatal, "in function \033[1;31m%s at line %d Mismatch in the number of flags in configurable cfBookEventHistograms, and number of entries in enum eEventHistograms \n \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  eh.fBookEventHistograms[eNumberOfEvents] = static_cast<bool>(lBookEventHistograms[eNumberOfEvents]);
  eh.fBookEventHistograms[eTotalMultiplicity] = static_cast<bool>(lBookEventHistograms[eTotalMultiplicity]);
  eh.fBookEventHistograms[eSelectedTracks] = static_cast<bool>(lBookEventHistograms[eSelectedTracks]);
  eh.fBookEventHistograms[eMultFV0M] = static_cast<bool>(lBookEventHistograms[eMultFV0M]);
  eh.fBookEventHistograms[eMultFT0M] = static_cast<bool>(lBookEventHistograms[eMultFT0M]);
  eh.fBookEventHistograms[eMultTPC] = static_cast<bool>(lBookEventHistograms[eMultTPC]);
  eh.fBookEventHistograms[eMultNTracksPV] = static_cast<bool>(lBookEventHistograms[eMultNTracksPV]);
  eh.fBookEventHistograms[eCentrality] = static_cast<bool>(lBookEventHistograms[eCentrality]);
  eh.fBookEventHistograms[eVertex_x] = static_cast<bool>(lBookEventHistograms[eVertex_x]);
  eh.fBookEventHistograms[eVertex_y] = static_cast<bool>(lBookEventHistograms[eVertex_y]);
  eh.fBookEventHistograms[eVertex_z] = static_cast<bool>(lBookEventHistograms[eVertex_z]);
  eh.fBookEventHistograms[eNContributors] = static_cast<bool>(lBookEventHistograms[eNContributors]);
  eh.fBookEventHistograms[eImpactParameter] = static_cast<bool>(lBookEventHistograms[eImpactParameter]);

  // b) Particle histograms:
  // By default all particle histograms are booked. If you do not want particular particle histogram to be booked,
  // use configurable array cfBookParticleHistograms, where you can specify flags 1 (book) or 0 (do not book).
  // Ordering of the flags in that array is interpreted through ordering of enums in enum eParticleHistograms. // TBI 20240124 is this safe enough?
  auto lBookParticleHistograms = (vector<int>)cfBookParticleHistograms; // this is now the local version of that int array from configurable. TBI 20240124 why is this casting mandatory?
  if (lBookParticleHistograms.size() != eParticleHistograms_N) {
    LOGF(fatal, "in function \033[1;31m%s at line %d Mismatch in the number of flags in configurable cfBookParticleHistograms, and number of entries in enum eParticleHistograms \n \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  ph.fBookParticleHistograms[ePhi] = static_cast<bool>(lBookParticleHistograms[ePhi]);
  ph.fBookParticleHistograms[ePt] = static_cast<bool>(lBookParticleHistograms[ePt]);
  ph.fBookParticleHistograms[eEta] = static_cast<bool>(lBookParticleHistograms[eEta]);
  ph.fBookParticleHistograms[etpcNClsCrossedRows] = static_cast<bool>(lBookParticleHistograms[etpcNClsCrossedRows]);
  ph.fBookParticleHistograms[eDCA_xy] = static_cast<bool>(lBookParticleHistograms[eDCA_xy]);
  ph.fBookParticleHistograms[eDCA_z] = static_cast<bool>(lBookParticleHistograms[eDCA_z]);

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
  eh.fEventHistogramsBins[eNumberOfEvents][0] = 1;
  eh.fEventHistogramsBins[eNumberOfEvents][1] = 0.;
  eh.fEventHistogramsBins[eNumberOfEvents][2] = 1.;
  // task->SetEventHistogramsBins("TotalMultiplicity",1000,0.,1000.);
  eh.fEventHistogramsBins[eTotalMultiplicity][0] = 10000;
  eh.fEventHistogramsBins[eTotalMultiplicity][1] = 0.;
  eh.fEventHistogramsBins[eTotalMultiplicity][2] = 10000.;
  // task->SetEventHistogramsBins("SelectedTracks",1000,0.,1000.);
  eh.fEventHistogramsBins[eSelectedTracks][0] = 10000;
  eh.fEventHistogramsBins[eSelectedTracks][1] = 0.;
  eh.fEventHistogramsBins[eSelectedTracks][2] = 10000.;
  // ... TBI 20240120
  eh.fEventHistogramsBins[eMultFV0M][0] = 10000;
  eh.fEventHistogramsBins[eMultFV0M][1] = 0.;
  eh.fEventHistogramsBins[eMultFV0M][2] = 10000.;
  // ... TBI 20240120
  eh.fEventHistogramsBins[eMultFT0M][0] = 10000;
  eh.fEventHistogramsBins[eMultFT0M][1] = 0.;
  eh.fEventHistogramsBins[eMultFT0M][2] = 10000.;
  // ... TBI 20240120
  eh.fEventHistogramsBins[eMultTPC][0] = 10000;
  eh.fEventHistogramsBins[eMultTPC][1] = 0.;
  eh.fEventHistogramsBins[eMultTPC][2] = 10000.;
  // ... TBI 20240120
  eh.fEventHistogramsBins[eMultNTracksPV][0] = 10000;
  eh.fEventHistogramsBins[eMultNTracksPV][1] = 0.;
  eh.fEventHistogramsBins[eMultNTracksPV][2] = 10000.;
  // task->SetEventHistogramsBins("Centrality",100,0.,100.);
  eh.fEventHistogramsBins[eCentrality][0] = 110; // intentionally, because if centrality is not determined, it's set to 105.0 at the moment
  eh.fEventHistogramsBins[eCentrality][1] = 0.;
  eh.fEventHistogramsBins[eCentrality][2] = 110.;
  // task->SetEventHistogramsBins("Vertex_x",1000,-20.,20.);
  eh.fEventHistogramsBins[eVertex_x][0] = 1000;
  eh.fEventHistogramsBins[eVertex_x][1] = -20.;
  eh.fEventHistogramsBins[eVertex_x][2] = 20.;
  // task->SetEventHistogramsBins("Vertex_y",1000,-20.,20.);
  eh.fEventHistogramsBins[eVertex_y][0] = 1000;
  eh.fEventHistogramsBins[eVertex_y][1] = -20.;
  eh.fEventHistogramsBins[eVertex_y][2] = 20.;
  // task->SetEventHistogramsBins("Vertex_z",1000,-20.,20.);
  eh.fEventHistogramsBins[eVertex_z][0] = 1000;
  eh.fEventHistogramsBins[eVertex_z][1] = -20.;
  eh.fEventHistogramsBins[eVertex_z][2] = 20.;
  // task->SetEventHistogramsBins("NContributors",1000,0.,1000.);
  eh.fEventHistogramsBins[eNContributors][0] = 1000;
  eh.fEventHistogramsBins[eNContributors][1] = 0.;
  eh.fEventHistogramsBins[eNContributors][2] = 1000.;
  // task->SetEventHistogramsBins("ImpactParameter",1000,0.,1000.);
  eh.fEventHistogramsBins[eImpactParameter][0] = 1000000;
  eh.fEventHistogramsBins[eImpactParameter][1] = 0.;
  eh.fEventHistogramsBins[eImpactParameter][2] = 1.; // TBI 20231031 check this, i do not know in which units IP is stored

  // b) Default binning for particle histograms:
  // task->SetParticleHistogramsBins("Phi",360,0.,TMath::TwoPi());
  ph.fParticleHistogramsBins[ePhi][0] = 360;
  ph.fParticleHistogramsBins[ePhi][1] = 0.;
  ph.fParticleHistogramsBins[ePhi][2] = TMath::TwoPi();
  // task->SetParticleHistogramsBins("Pt",1000,0.,20.);
  ph.fParticleHistogramsBins[ePt][0] = 1000;
  ph.fParticleHistogramsBins[ePt][1] = 0.;
  ph.fParticleHistogramsBins[ePt][2] = 20.;
  // task->SetParticleHistogramsBins("Eta",200,-1.,1.);
  ph.fParticleHistogramsBins[eEta][0] = 200;
  ph.fParticleHistogramsBins[eEta][1] = -2.;
  ph.fParticleHistogramsBins[eEta][2] = 2.;

  ph.fParticleHistogramsBins[etpcNClsCrossedRows][0] = 200;
  ph.fParticleHistogramsBins[etpcNClsCrossedRows][1] = 0.;
  ph.fParticleHistogramsBins[etpcNClsCrossedRows][2] = 200.;

  ph.fParticleHistogramsBins[eDCA_xy][0] = 2000;
  ph.fParticleHistogramsBins[eDCA_xy][1] = -10.;
  ph.fParticleHistogramsBins[eDCA_xy][2] = 10.;

  ph.fParticleHistogramsBins[eDCA_z][0] = 2000;
  ph.fParticleHistogramsBins[eDCA_z][1] = -10.;
  ph.fParticleHistogramsBins[eDCA_z][2] = 10.;

  ph.fParticleHistogramsBins[ePDG][0] = 10000;
  ph.fParticleHistogramsBins[ePDG][1] = -5000.;
  ph.fParticleHistogramsBins[ePDG][2] = 5000.;

  // c) Default binning for results histograms:
  //    Remark: These bins apply to following categories fCorrelationsPro, fNestedLoopsPro, fTest0Pro, and fResultsPro.
  // *) For integrated resullts, binning is always the same:
  res.fResultsProFixedLengthBins[AFO_INTEGRATED][0] = 1;
  res.fResultsProFixedLengthBins[AFO_INTEGRATED][1] = 0.;
  res.fResultsProFixedLengthBins[AFO_INTEGRATED][2] = 1.;

  // *) Fixed-length binning vs. multiplicity:
  auto lFixedLength_mult_bins = (vector<float>)cFixedLength_mult_bins; // this is now the local version of that float array from configurable.
  if (lFixedLength_mult_bins.size() != 3) {
    LOGF(fatal, "in function \033[1;31m%s at line %d => The array cFixedLength_mult_bins must have 3 entries: {nBins, min, max} \n \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }
  res.fResultsProFixedLengthBins[AFO_MULTIPLICITY][0] = lFixedLength_mult_bins[0];
  res.fResultsProFixedLengthBins[AFO_MULTIPLICITY][1] = lFixedLength_mult_bins[1];
  res.fResultsProFixedLengthBins[AFO_MULTIPLICITY][2] = lFixedLength_mult_bins[2];

  // *) Fixed-length binning vs. centrality:
  auto lFixedLength_cent_bins = (vector<float>)cFixedLength_cent_bins; // this is now the local version of that float array from configurable.
  if (lFixedLength_cent_bins.size() != 3) {
    LOGF(fatal, "in function \033[1;31m%s at line %d => The array cFixedLength_cent_bins must have 3 entries: {nBins, min, max} \n \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }
  res.fResultsProFixedLengthBins[AFO_CENTRALITY][0] = lFixedLength_cent_bins[0];
  res.fResultsProFixedLengthBins[AFO_CENTRALITY][1] = lFixedLength_cent_bins[1];
  res.fResultsProFixedLengthBins[AFO_CENTRALITY][2] = lFixedLength_cent_bins[2];

  // *) Fixed-length binning vs. pt:
  auto lFixedLength_pt_bins = (vector<float>)cFixedLength_pt_bins; // this is now the local version of that float array from configurable.
  if (lFixedLength_pt_bins.size() != 3) {
    LOGF(fatal, "in function \033[1;31m%s at line %d => The array cFixedLength_pt_bins must have 3 entries: {nBins, min, max} \n \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }
  res.fResultsProFixedLengthBins[AFO_PT][0] = lFixedLength_pt_bins[0];
  res.fResultsProFixedLengthBins[AFO_PT][1] = lFixedLength_pt_bins[1];
  res.fResultsProFixedLengthBins[AFO_PT][2] = lFixedLength_pt_bins[2];

  // *) Fixed-length binning vs. eta:
  auto lFixedLength_eta_bins = (vector<float>)cFixedLength_eta_bins; // this is now the local version of that float array from configurable.
  if (lFixedLength_eta_bins.size() != 3) {
    LOGF(fatal, "in function \033[1;31m%s at line %d => The array cFixedLength_eta_bins must have 3 entries: {nBins, min, max} \n \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }
  res.fResultsProFixedLengthBins[AFO_ETA][0] = lFixedLength_eta_bins[0];
  res.fResultsProFixedLengthBins[AFO_ETA][1] = lFixedLength_eta_bins[1];
  res.fResultsProFixedLengthBins[AFO_ETA][2] = lFixedLength_eta_bins[2];

  // d) Variable-length binning set via MuPa-Configurables.h:
  // *) Variable-length binning vs. multiplicity:
  if (cUseVariableLength_mult_bins) {
    res.fUseResultsProVariableLengthBins[AFO_MULTIPLICITY] = kTRUE;
    res.fResultsProVariableLengthBinsString[AFO_MULTIPLICITY] = cVariableLength_mult_bins;
    this->CastStringIntoArray(AFO_MULTIPLICITY);
  }
  // *) Variable-length binning vs. centrality:
  if (cUseVariableLength_cent_bins) {
    res.fUseResultsProVariableLengthBins[AFO_CENTRALITY] = kTRUE;
    res.fResultsProVariableLengthBinsString[AFO_CENTRALITY] = cVariableLength_cent_bins;
    this->CastStringIntoArray(AFO_CENTRALITY);
  }
  // *) Variable-length binning vs. pt:
  if (cUseVariableLength_pt_bins) {
    res.fUseResultsProVariableLengthBins[AFO_PT] = kTRUE;
    res.fResultsProVariableLengthBinsString[AFO_PT] = cVariableLength_pt_bins;
    this->CastStringIntoArray(AFO_PT);
  }
  // *) Variable-length binning vs. eta:
  if (cUseVariableLength_eta_bins) {
    res.fUseResultsProVariableLengthBins[AFO_ETA] = kTRUE;
    res.fResultsProVariableLengthBinsString[AFO_ETA] = cVariableLength_eta_bins;
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
    LOGF(info, "\033[1;32m Casting a string %s into TArrayD .... \033[0m", res.fResultsProVariableLengthBinsString[AFO].Data());
  }

  TObjArray* oa = res.fResultsProVariableLengthBinsString[AFO].Tokenize(",");
  if (!oa) {
    LOGF(fatal, "in function \033[1;31m%s at line %d \n fResultsProVariableLengthBinsString[AFO] = %s\033[0m", __PRETTY_FUNCTION__, __LINE__, res.fResultsProVariableLengthBinsString[AFO].Data());
  }
  Int_t nEntries = oa->GetEntries();
  res.fResultsProVariableLengthBins[AFO] = new TArrayD(nEntries);
  for (Int_t i = 0; i < nEntries; i++) {
    // cout<< TString(oa->At(i)->GetName()).Atof() <<endl;
    res.fResultsProVariableLengthBins[AFO]->AddAt(TString(oa->At(i)->GetName()).Atof(), i);
  }
  delete oa; // yes, otherwise it's a memory leak

  if (tc.fVerbose) {
    for (Int_t i = 0; i < res.fResultsProVariableLengthBins[AFO]->GetSize(); i++) {
      LOGF(info, "\033[1;32m [%d] : %f \033[0m", i, res.fResultsProVariableLengthBins[AFO]->At(i));
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
  eh.fEventCuts[eNumberOfEvents][eMin] =
    cNumberOfEvents_min; // Configurable<int>
                         // cNumberOfEvents_min{"cNumberOfEvents_min", ...
  eh.fEventCuts[eNumberOfEvents][eMax] =
    cNumberOfEvents_max; // Configurable<int>
                         // cNumberOfEvents_max{"cNumberOfEvents_max", ...

  eh.fEventCuts[eTotalMultiplicity][eMin] = cTotalMultiplicity_min; // Configurable<int>
                                                                    // cTotalMultiplicity_min{"cTotalMultiplicity_min",
                                                                    // ...
  eh.fEventCuts[eTotalMultiplicity][eMax] = cTotalMultiplicity_max; // Configurable<int>
                                                                    // cTotalMultiplicity_max{"cTotalMultiplicity_max",
                                                                    // ...

  eh.fEventCuts[eSelectedTracks][eMin] =
    cSelectedTracks_min; // Configurable<int>
                         // cSelectedTracks_min{"cSelectedTracks_min", ...
  eh.fEventCuts[eSelectedTracks][eMax] =
    cSelectedTracks_max; // Configurable<int>
                         // cSelectedTracks_max{"cSelectedTracks_max", ...

  eh.fEventCuts[eCentrality][eMin] =
    cCentrality_min; // Configurable<int> cCentrality_min{"cCentrality_min",
                     // ...
  eh.fEventCuts[eCentrality][eMax] =
    cCentrality_max; // Configurable<int> cCentrality_max{"cCentrality_max",
                     // ...

  eh.fEventCuts[eVertex_x][eMin] =
    cVertex_x_min; // Configurable<int> cVertex_x_min{"cVertex_x_min", ...
  eh.fEventCuts[eVertex_x][eMax] =
    cVertex_x_max; // Configurable<int> cVertex_x_max{"cVertex_x_max", ...

  eh.fEventCuts[eVertex_y][eMin] =
    cVertex_y_min; // Configurable<int> cVertex_y_min{"cVertex_y_min", ...
  eh.fEventCuts[eVertex_y][eMax] =
    cVertex_y_max; // Configurable<int> cVertex_y_max{"cVertex_y_max", ...

  eh.fEventCuts[eVertex_z][eMin] =
    cVertex_z_min; // Configurable<int> cVertex_z_min{"cVertex_z_min", ...
  eh.fEventCuts[eVertex_z][eMax] =
    cVertex_z_max; // Configurable<int> cVertex_z_max{"cVertex_z_max", ...

  eh.fEventCuts[eNContributors][eMin] =
    cNContributors_min; // Configurable<int>
                        // cNContributors_min{"cNContributors_min", ...
  eh.fEventCuts[eNContributors][eMax] =
    cNContributors_max; // Configurable<int>
                        // cNContributors_max{"cNContributors_max", ...

  eh.fEventCuts[eImpactParameter][eMin] =
    cImpactParameter_min; // Configurable<int>
                          // cImpactParameter_min{"cImpactParameter_min", ...
  eh.fEventCuts[eImpactParameter][eMax] =
    cImpactParameter_max; // Configurable<int>
                          // cImpactParameter_max{"cImpactParameter_max", ...

  // ...

  // b) Default particle cuts:

} // void DefaultCuts()

//============================================================

void InsanityChecks()
{
  // Do insanity checks on configuration, booking, binning and cuts.

  // a) Insanity checks on configuration;
  // b) Insanity checks on event cuts;
  // c) Insanity checks on booking;
  // d) Insanity checks on binning;
  // e) Insanity checks on cuts.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Insanity checks on configuration:

  // *) Insanity check on individual flags: Make sure that only one process is set to kTRUE.
  //    If 2 or more are kTRUE, then corresponding process function is executed over ALL data, then another process(...) function, etc.
  //    Re-think this if it's possible to run different process(...)'s concurently over the same data.
  if ((Int_t)tc.fProcessRec + (Int_t)tc.fProcessRecSim + (Int_t)tc.fProcessSim + (Int_t)tc.fProcessRec_Run2 + (Int_t)tc.fProcessRecSim_Run2 + (Int_t)tc.fProcessSim_Run2 + (Int_t)tc.fProcessRec_Run1 + (Int_t)tc.fProcessRecSim_Run1 + (Int_t)tc.fProcessSim_Run1 > 1) {
    LOGF(info, "\033[1;31m Only one flag can be kTRUE: tc.fProcessRec = %d, tc.fProcessRecSim = %d, tc.fProcessSim = %d, tc.fProcessRec_Run2 = %d, tc.fProcessRecSim_Run2 = %d, tc.fProcessSim_Run2 = %d, tc.fProcessRec_Run1 = %d, tc.fProcessRecSim_Run1 = %d, tc.fProcessSim_Run1 = %d \033[0m", (Int_t)tc.fProcessRec, (Int_t)tc.fProcessRecSim, (Int_t)tc.fProcessSim, (Int_t)tc.fProcessRec_Run2, (Int_t)tc.fProcessRecSim_Run2, (Int_t)tc.fProcessSim_Run2, (Int_t)tc.fProcessRec_Run1, (Int_t)tc.fProcessRecSim_Run1, (Int_t)tc.fProcessSim_Run1);
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  // *) Seed for random number generator must be non-negative integer:
  // if (tc.fRandomSeed < 0) {
  //   LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  // }

  // *) Insanity checks on event cuts:
  if (tc.fProcessRec_Run2 || tc.fProcessRec_Run1) { // From documentation: Bypass this check if you analyse MC or continuous Run3 data.
    if (!(ec.fTrigger.EqualTo("kINT7"))) {          // TBI 20240223 expand this list with other supported triggers eventually in this category (see if(...) above)
      LOGF(info, "in function \033[1;32m%s at line %d : trigger \"%s\" is not internally supported yet. Add it to the list of supported triggers, if you really want to use that one.\033[0m", __PRETTY_FUNCTION__, __LINE__, ec.fTrigger.Data());
      ec.fUseTrigger = kFALSE;
    } else {
      ec.fUseTrigger = kTRUE; // I am analyzing converted Run 1 or Run 2 real data (not MC!), and the trigger is supported, so let's use it
    }
  }

  if (ec.fUseSel7) { // from doc: for Run 2 data and MC
    if (!(tc.fProcessRec_Run2 || tc.fProcessRecSim_Run2 || tc.fProcessSim_Run2 || tc.fProcessRec_Run1 || tc.fProcessRecSim_Run1 || tc.fProcessSim_Run1)) {
      LOGF(fatal, "in function \033[1;31m%s at line %d use fUseSel7 for Run 2 data and MC\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseSel8) { // from doc: for Run 3 data and MC
    if (!(tc.fProcessRec || tc.fProcessRecSim || tc.fProcessSim)) {
      LOGF(fatal, "in function \033[1;31m%s at line %d use fUseSel8 for Run 3 data and MC\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }
  }

  if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && !tc.fUseFisherYates) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  if (tc.fProcessRec || tc.fProcessRecSim || tc.fProcessSim) {
    // Supported centrality estimators for Run 3 are enlisted here:
    if (!(ec.fCentralityEstimator.EqualTo("centFT0M", TString::kIgnoreCase) ||
          ec.fCentralityEstimator.EqualTo("centFV0A", TString::kIgnoreCase) ||
          ec.fCentralityEstimator.EqualTo("centNTPV", TString::kIgnoreCase))) {
      LOGF(fatal, "in function \033[1;31m%s at line %d. centrality estimator = %d is not supported yet for Run 3 analysis. \033[0m", __PRETTY_FUNCTION__, __LINE__, ec.fCentralityEstimator.Data());
    }
  }

  if (tc.fProcessRec_Run2 || tc.fProcessRecSim_Run2 || tc.fProcessSim_Run2 || tc.fProcessRec_Run1 || tc.fProcessRecSim_Run1 || tc.fProcessSim_Run1) {
    // Supported centrality estimators for Run 3 are enlisted here:
    if (!(ec.fCentralityEstimator.EqualTo("centRun2V0M", TString::kIgnoreCase) ||
          ec.fCentralityEstimator.EqualTo("centRun2SPDTracklets", TString::kIgnoreCase))) {
      LOGF(fatal, "in function \033[1;31m%s at line %d. centrality estimator = %d is not supported yet for converted Run 2 and Run 1 analysis. \033[0m", __PRETTY_FUNCTION__, __LINE__, ec.fCentralityEstimator.Data());
    }
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
  qa.fQAList = new TList();
  qa.fQAList->SetName("QA");
  qa.fQAList->SetOwner(kTRUE);
  fBaseList->Add(qa.fQAList);

  // *) Control event histograms:
  eh.fEventHistogramsList = new TList();
  eh.fEventHistogramsList->SetName("EventHistograms");
  eh.fEventHistogramsList->SetOwner(kTRUE);
  fBaseList->Add(eh.fEventHistogramsList);

  // *) Event cuts:
  ec.fEventCutsList = new TList();
  ec.fEventCutsList->SetName("EventCuts");
  ec.fEventCutsList->SetOwner(kTRUE);
  fBaseList->Add(ec.fEventCutsList);

  // *) Control particle histograms:
  ph.fParticleHistogramsList = new TList();
  ph.fParticleHistogramsList->SetName("ParticleHistograms");
  ph.fParticleHistogramsList->SetOwner(kTRUE);
  fBaseList->Add(ph.fParticleHistogramsList);

  // *) Particle cuts:
  pc.fParticleCutsList = new TList();
  pc.fParticleCutsList->SetName("ParticleCuts");
  pc.fParticleCutsList->SetOwner(kTRUE);
  fBaseList->Add(pc.fParticleCutsList);

  // *) Q-vectors:
  qv.fQvectorList = new TList();
  qv.fQvectorList->SetName("Q-vectors");
  qv.fQvectorList->SetOwner(kTRUE);
  fBaseList->Add(qv.fQvectorList);

  // *) Correlations:
  mupa.fCorrelationsList = new TList();
  mupa.fCorrelationsList->SetName("Correlations");
  mupa.fCorrelationsList->SetOwner(kTRUE);
  fBaseList->Add(mupa.fCorrelationsList);

  // *) Particle weights:
  pw.fWeightsList = new TList();
  pw.fWeightsList->SetName("Weights");
  pw.fWeightsList->SetOwner(kTRUE);
  fBaseList->Add(pw.fWeightsList);

  // *) Nested loops:
  nl.fNestedLoopsList = new TList();
  nl.fNestedLoopsList->SetName("NestedLoops");
  nl.fNestedLoopsList->SetOwner(kTRUE);
  fBaseList->Add(nl.fNestedLoopsList);

  // *) Test0:
  t0.fTest0List = new TList();
  t0.fTest0List->SetName("Test0");
  t0.fTest0List->SetOwner(kTRUE);
  fBaseList->Add(t0.fTest0List);

  // *) Results:
  res.fResultsList = new TList();
  res.fResultsList->SetName("Results");
  res.fResultsList->SetOwner(kTRUE);
  fBaseList->Add(res.fResultsList);

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
  eh.fEventHistogramsPro = new TProfile("fEventHistogramsPro",
                                        "flags for event histograms", 25, 0., 25.);
  eh.fEventHistogramsPro->SetStats(kFALSE);
  eh.fEventHistogramsPro->SetLineColor(eColor);
  eh.fEventHistogramsPro->SetFillColor(eFillColor);
  // ...
  eh.fEventHistogramsList->Add(eh.fEventHistogramsPro);

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
    if (!eh.fBookEventHistograms[t]) {
      continue;
    }
    for (Int_t rs = 0; rs < 2; rs++) // reco/sim
    {
      // If I am analyzing only reconstructed data, do not book histos for simulated, and vice versa.
      // TBI 20240223 tc.fProcessTest is treated as tc.fProcessRec, for the time being
      if ((tc.fGenericRec && rs == eSim) || (tc.fGenericSim && rs == eRec)) {
        continue;
      }

      for (Int_t ba = 0; ba < 2; ba++) // before/after cuts
      {
        eh.fEventHistograms[t][rs][ba] = new TH1D(
          Form("fEventHistograms[%s][%s][%s]", stype[t].Data(),
               srs[rs].Data(), sba[ba].Data()),
          Form("%s, %s, %s", stype[t].Data(), srs[rs].Data(), sba[ba].Data()),
          (Int_t)eh.fEventHistogramsBins[t][0],
          eh.fEventHistogramsBins[t][1], eh.fEventHistogramsBins[t][2]);
        eh.fEventHistograms[t][rs][ba]->SetLineColor(fBeforeAfterColor[ba]);
        eh.fEventHistograms[t][rs][ba]->SetFillColor(fBeforeAfterColor[ba] -
                                                     10);
        eh.fEventHistogramsList->Add(eh.fEventHistograms[t][rs][ba]);
      } // for(Int_t ba=0;ba<2;ba++)
    }   // for(Int_t rs=0;rs<2;rs++) // reco/sim
  }     // for(Int_t t=0;t<eEventHistograms_N;t++) // type, see enum
        // eEventHistograms

} // void BookEventHistograms()

//============================================================

void BookEventCutsHistograms()
{
  // Book all event cuts objects.

  // a) Book the profile holding event cuts flags;

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  ec.fEventCutsPro = new TProfile("fEventCutsPro",
                                  "flags for event cuts", eEventCuts_N - 1, 0., eEventCuts_N - 1);
  ec.fEventCutsPro->SetStats(kFALSE);
  ec.fEventCutsPro->SetLineColor(eColor);
  ec.fEventCutsPro->SetFillColor(eFillColor);

  ec.fEventCutsPro->GetXaxis()->SetBinLabel(eTrigger, Form("fTrigger = %s", ec.fTrigger.Data()));
  ec.fEventCutsPro->GetXaxis()->SetBinLabel(eUseTrigger, "fUseTrigger");
  ec.fEventCutsPro->Fill(+eUseTrigger - 0.5, (Int_t)ec.fUseTrigger);
  ec.fEventCutsPro->GetXaxis()->SetBinLabel(eUseSel7, "fUseSel7");
  ec.fEventCutsPro->Fill(+eUseSel7 - 0.5, (Int_t)ec.fUseSel7);
  ec.fEventCutsPro->GetXaxis()->SetBinLabel(eUseSel8, "fUseSel8");
  ec.fEventCutsPro->Fill(+eUseSel8 - 0.5, (Int_t)ec.fUseSel8);
  ec.fEventCutsPro->GetXaxis()->SetBinLabel(eCentralityEstimator, Form("fCentralityEstimator = %s", ec.fCentralityEstimator.Data()));

  ec.fEventCutsList->Add(ec.fEventCutsPro);

} // void BookEventCutsHistograms()

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
  ph.fParticleHistogramsPro = new TProfile(
    "fParticleHistogramsPro", "flags for particle histograms", 25, 0., 25.);
  ph.fParticleHistogramsPro->SetStats(kFALSE);
  ph.fParticleHistogramsPro->SetLineColor(eColor);
  ph.fParticleHistogramsPro->SetFillColor(eFillColor);
  // ...
  ph.fParticleHistogramsList->Add(ph.fParticleHistogramsPro);

  Int_t fBeforeAfterColor[2] = {
    kRed,
    kGreen}; //! [0 = kRed,1 = kGreen] TBI 20220713 only temporarily here

  // b) Book specific control particle histograms:
  TString stype[eParticleHistograms_N] = {
    "Phi", "Pt", "Eta", "tpcNClsCrossedRows",
    "DCA_xy", "DCA_z", "PDG"}; // keep in sync. with enum eParticleHistograms
  TString srs[2] = {"rec", "sim"};
  TString sba[2] = {"before", "after"};

  for (Int_t t = 0; t < eParticleHistograms_N;
       t++) // type, see enum eParticleHistograms
  {
    if (!ph.fBookParticleHistograms[t]) {
      continue;
    }
    for (Int_t rs = 0; rs < 2; rs++) // reco/sim
    {
      if ((tc.fGenericRec && rs == eSim) || (tc.fGenericSim && rs == eRec)) {
        continue; // if I am analyzing only reconstructed data, do not book histos for simulated, and vice versa.
      }
      for (Int_t ba = 0; ba < 2; ba++) // before/after cuts
      {
        ph.fParticleHistograms[t][rs][ba] = new TH1D(
          Form("fParticleHistograms[%s][%s][%s]", stype[t].Data(),
               srs[rs].Data(), sba[ba].Data()),
          Form("%s, %s, %s", stype[t].Data(), srs[rs].Data(), sba[ba].Data()),
          (Int_t)ph.fParticleHistogramsBins[t][0],
          ph.fParticleHistogramsBins[t][1],
          ph.fParticleHistogramsBins[t][2]);
        ph.fParticleHistograms[t][rs][ba]->SetLineColor(
          fBeforeAfterColor[ba]);
        ph.fParticleHistograms[t][rs][ba]->SetFillColor(
          fBeforeAfterColor[ba] - 10);
        ph.fParticleHistogramsList->Add(ph.fParticleHistograms[t][rs][ba]);
      } // for(Int_t ba=0;ba<2;ba++)
    }   // for(Int_t rs=0;rs<2;rs++) // reco/sim
  }     // for(Int_t t=0;t<eParticleHistograms_N;t++) // type, see enum
        // eParticleHistograms

} // void BookParticleHistograms()

//============================================================

void BookParticleCutsHistograms()
{
  // Book all particle cuts objects.

  // a) Book the profile holding particle cuts flags;

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  pc.fParticleCutsPro = new TProfile("fParticleCutsPro",
                                     "flags for particle cuts", eParticleCuts_N - 1, 0., eParticleCuts_N - 1);
  pc.fParticleCutsPro->SetStats(kFALSE);
  pc.fParticleCutsPro->SetLineColor(eColor);
  pc.fParticleCutsPro->SetFillColor(eFillColor);

  pc.fParticleCutsPro->GetXaxis()->SetBinLabel(eTBI, "TBI");
  // ...

  pc.fParticleCutsList->Add(pc.fParticleCutsPro);

} // void BookParticleCutsHistograms()

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
  qv.fQvectorFlagsPro =
    new TProfile("fQvectorFlagsPro", "flags for Q-vector objects", 3, 0., 3.);
  qv.fQvectorFlagsPro->SetStats(kFALSE);
  qv.fQvectorFlagsPro->SetLineColor(eColor);
  qv.fQvectorFlagsPro->SetFillColor(eFillColor);
  qv.fQvectorFlagsPro->GetXaxis()->SetLabelSize(0.05);
  qv.fQvectorFlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateQvectors");
  qv.fQvectorFlagsPro->Fill(0.5, qv.fCalculateQvectors);
  qv.fQvectorFlagsPro->GetXaxis()->SetBinLabel(2, "gMaxHarmonic");
  qv.fQvectorFlagsPro->Fill(1.5, gMaxHarmonic);
  qv.fQvectorFlagsPro->GetXaxis()->SetBinLabel(3, "gMaxCorrelator");
  qv.fQvectorFlagsPro->Fill(2.5, gMaxCorrelator);
  qv.fQvectorList->Add(qv.fQvectorFlagsPro);

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
  mupa.fCorrelationsFlagsPro = new TProfile("fCorrelationsFlagsPro",
                                            "flags for correlations", 1, 0., 31);
  mupa.fCorrelationsFlagsPro->SetStats(kFALSE);
  mupa.fCorrelationsFlagsPro->SetLineColor(eColor);
  mupa.fCorrelationsFlagsPro->SetFillColor(eFillColor);
  mupa.fCorrelationsFlagsPro->GetXaxis()->SetLabelSize(0.05);
  mupa.fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateCorrelations");
  mupa.fCorrelationsFlagsPro->Fill(0.5, mupa.fCalculateCorrelations);
  // ...
  mupa.fCorrelationsList->Add(mupa.fCorrelationsFlagsPro);

  if (!mupa.fCalculateCorrelations) {
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
        if (!res.fResultsPro[v]) {
          LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
        }
        mupa.fCorrelationsPro[k][n][v] = reinterpret_cast<TProfile*>(res.fResultsPro[v]->Clone(Form("fCorrelationsPro[%d][%d][%s]", k, n, res.fResultsProRawName[v].Data()))); // yes
        mupa.fCorrelationsPro[k][n][v]->SetStats(kFALSE);
        mupa.fCorrelationsPro[k][n][v]->Sumw2();
        mupa.fCorrelationsPro[k][n][v]->GetXaxis()->SetTitle(res.fResultsProXaxisTitle[v].Data());
        mupa.fCorrelationsPro[k][n][v]->GetYaxis()->SetTitle(Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1), oVariable[k].Data()));
        mupa.fCorrelationsList->Add(mupa.fCorrelationsPro[k][n][v]);
      }
    } // for (Int_t n = 0; n < gMaxHarmonic; n++) // harmonic
  }   // for (Int_t k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]

  // d) Few quick insanity checks on booking:
  if (mupa.fCorrelationsPro[0][0][AFO_INTEGRATED] && !TString(mupa.fCorrelationsPro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }
  if (mupa.fCorrelationsPro[0][0][AFO_PT] && !TString(mupa.fCorrelationsPro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("p_{T}")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }

} // BookCorrelationsHistograms()

//============================================================

void BookWeightsHistograms()
{
  // Book all objects for particle weights.

  // a) Book the profile holding flags;
  // b) Histograms;
  // c) Histograms for differential weights.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  pw.fWeightsFlagsPro =
    new TProfile("fWeightsFlagsPro", "flags for particle weights", 5, 0., 5.);
  pw.fWeightsFlagsPro->SetStats(kFALSE);
  pw.fWeightsFlagsPro->SetLineColor(eColor);
  pw.fWeightsFlagsPro->SetFillColor(eFillColor);
  pw.fWeightsFlagsPro->GetXaxis()->SetLabelSize(0.05);
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(1, "w_{#varphi}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(2, "w_{p_{t}}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(3, "w_{#eta}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(4, "w_{#varphi}(p_{t})");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(5, "w_{#varphi}(#eta)");

  for (Int_t w = 0; w < eWeights_N; w++) // use weights [phi,pt,eta]
  {
    if (pw.fUseWeights[w])
      pw.fWeightsFlagsPro->Fill(w + 0.5, 1.);
  }
  for (Int_t w = 0; w < eDiffWeights_N; w++) // use differential weights [phipt,phieta,...]
  {
    if (pw.fUseDiffWeights[w])
      pw.fWeightsFlagsPro->Fill(w + 3.5, 1.); // TBI 20231026 This hadrwired offset of +3.5 will bite me sooner or later, but nevermind now...
  }
  pw.fWeightsList->Add(pw.fWeightsFlagsPro);

  // b) Histograms:
  //    As of 20240216, I have abandoned the idea to generate integrated weights internally, weights
  //    are always fetched and cloned from external files, in any case (local, AliEn, CCDB).
  //    Therefore, add histos with weights to this list only after they are cloned from external files.

  // c) Histograms for differential weights:
  //    Same comment applies as for c) => add histograms to the list, only after they are cloned from external files.

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
  nl.fNestedLoopsFlagsPro =
    new TProfile("fNestedLoopsFlagsPro", "flags for nested loops", 2, 0., 2.);
  nl.fNestedLoopsFlagsPro->SetStats(kFALSE);
  nl.fNestedLoopsFlagsPro->GetXaxis()->SetLabelSize(0.05);
  nl.fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateNestedLoops");
  nl.fNestedLoopsFlagsPro->Fill(0.5, nl.fCalculateNestedLoops);
  nl.fNestedLoopsFlagsPro->Fill(1.5, nl.fCalculateCustomNestedLoop);
  nl.fNestedLoopsList->Add(nl.fNestedLoopsFlagsPro);

  if (!(nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoop)) {
    return;
  }

  const Int_t iMaxSize = 2e4;
  nl.ftaNestedLoops[0] =
    new TArrayD(iMaxSize); // ebe container for azimuthal angles
  nl.ftaNestedLoops[1] = new TArrayD(
    iMaxSize); // ebe container for particle weights (product of all)

  // TBI 20220823 port here if(fCalculatePtCorrelations) { ... } and
  // if(fCalculateEtaCorrelations) { ... }

  if (!nl.fCalculateNestedLoops) {
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

        if (!res.fResultsPro[v]) {
          LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
        }
        nl.fNestedLoopsPro[k][n][v] = reinterpret_cast<TProfile*>(res.fResultsPro[v]->Clone(Form("fNestedLoopsPro[%d][%d][%d]", k, n, v))); // yes
        nl.fNestedLoopsPro[k][n][v]->SetTitle(Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1), oVariable[k].Data()));
        nl.fNestedLoopsPro[k][n][v]->SetStats(kFALSE);
        nl.fNestedLoopsPro[k][n][v]->Sumw2();
        nl.fNestedLoopsPro[k][n][v]->GetXaxis()->SetTitle(
          res.fResultsProXaxisTitle[v].Data());

        /*
        if(fUseFixedNumberOfRandomlySelectedTracks && 1==v) // just a warning
        for the meaning of multiplicity in this special case
        {
         nl.fNestedLoopsPro[k][n][1]->GetXaxis()->SetTitle("WARNING: for each
        multiplicity, fFixedNumberOfRandomlySelectedTracks is selected randomly
        in Q-vector");
        }
        */

        nl.fNestedLoopsList->Add(nl.fNestedLoopsPro[k][n][v]);
      } // for(Int_t v=0;v<5;v++) // variable [0=integrated,1=vs.
        // multiplicity,2=vs. centrality]
    }   // for (Int_t n = 0; n < gMaxHarmonic; n++) // harmonic
  }     // for (Int_t k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]

  // d) Few quick insanity checks on booking:
  if (nl.fNestedLoopsPro[0][0][AFO_INTEGRATED] && !TString(nl.fNestedLoopsPro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }
  if (nl.fNestedLoopsPro[0][0][AFO_PT] && !TString(nl.fNestedLoopsPro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("p_{T}")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
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
  t0.fTest0FlagsPro = new TProfile("fTest0FlagsPro", "flags for Test0", 1, 0., 1.);
  t0.fTest0FlagsPro->SetStats(kFALSE);
  t0.fTest0FlagsPro->GetXaxis()->SetLabelSize(0.04);
  t0.fTest0FlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateTest0");
  t0.fTest0FlagsPro->Fill(0.5, t0.fCalculateTest0);
  t0.fTest0List->Add(t0.fTest0FlagsPro);

  if (!t0.fCalculateTest0) {
    return;
  }

  // b) Book placeholder and make sure all labels are stored in the placeholder:
  this->StoreLabelsInPlaceholder();
  if (t0.fTest0LabelsPlaceholder) {
    t0.fTest0List->Add(t0.fTest0LabelsPlaceholder);
  }

  // c) Retrieve labels from placeholder:
  if (!(this->RetrieveCorrelationsLabels())) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  // d) Book what needs to be booked:
  for (Int_t mo = 0; mo < gMaxCorrelator; mo++) {
    for (Int_t mi = 0; mi < gMaxIndex; mi++) {
      if (!t0.fTest0Labels[mo][mi]) {
        continue;
      }
      {
        for (Int_t v = 0; v < eAsFunctionOf_N; v++) {
          // decide what is booked, then later valid pointer to fCorrelationsPro[k][n][v] is used as a boolean, in the standard way:
          if (AFO_INTEGRATED == v && !t0.fCalculateTest0AsFunctionOf[AFO_INTEGRATED]) {
            continue;
          }
          if (AFO_MULTIPLICITY == v && !t0.fCalculateTest0AsFunctionOf[AFO_MULTIPLICITY]) {
            continue;
          }
          if (AFO_CENTRALITY == v && !t0.fCalculateTest0AsFunctionOf[AFO_CENTRALITY]) {
            continue;
          }
          if (AFO_PT == v && !t0.fCalculateTest0AsFunctionOf[AFO_PT]) {
            continue;
          }
          if (AFO_ETA == v && !t0.fCalculateTest0AsFunctionOf[AFO_ETA]) {
            continue;
          }

          if (!res.fResultsPro[v]) {
            LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
          }

          t0.fTest0Pro[mo][mi][v] = reinterpret_cast<TProfile*>(res.fResultsPro[v]->Clone(Form("fTest0Pro[%d][%d][%s]", mo, mi, res.fResultsProRawName[v].Data()))); // yes
          t0.fTest0Pro[mo][mi][v]->SetStats(kFALSE);
          t0.fTest0Pro[mo][mi][v]->Sumw2();
          t0.fTest0Pro[mo][mi][v]->SetTitle(t0.fTest0Labels[mo][mi]->Data());
          t0.fTest0Pro[mo][mi][v]->GetXaxis()->SetTitle(res.fResultsProXaxisTitle[v].Data());
          /*
                if(fUseFixedNumberOfRandomlySelectedParticles && 1==v) // just a warning for the meaning of multiplicity in this special case
                {
                 fTest0Pro[mo][mi][1]->GetXaxis()->SetTitle("WARNING: for each multiplicity, fFixedNumberOfRandomlySelectedParticles is selected randomly in Q-vector");
                }
          */
          t0.fTest0List->Add(t0.fTest0Pro[mo][mi][v]); // yes, this has to be here
        }                                              // for(Int_t v=0;v<eAsFunctionOf_N;v++) // variable, see content of enum eAsFunctionOf
      }                                                // if(fTest0Labels[mo][mi])
    }                                                  // for(Int_t mi=0;mi<gMaxIndex;mi++)
  }                                                    // for(Int_t mo=0;mo<gMaxCorrelator;mo++)

  // e) Few quick insanity checks on booking:
  if (t0.fTest0Pro[0][0][AFO_INTEGRATED] && !TString(t0.fTest0Pro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }
  if (t0.fTest0Pro[0][0][AFO_PT] && !TString(t0.fTest0Pro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("p_{T}")) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }

} // void BookTest0Histograms()

//============================================================

void BookResultsHistograms()
{
  // Book all results histograms.

  // a) Book the profile holding flags;
  // b) Book results histograms, which in addition act as a sort of "abstract" interface, which defines common binning, etc., for other groups of histograms.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the profile holding flags:
  res.fResultsFlagsPro = new TProfile("fResultsFlagsPro",
                                      "flags for results histograms", 1, 0., 1.);
  res.fResultsFlagsPro->SetStats(kFALSE);
  res.fResultsFlagsPro->SetLineColor(eColor);
  res.fResultsFlagsPro->SetFillColor(eFillColor);
  res.fResultsFlagsPro->GetXaxis()->SetBinLabel(1, "fSaveResultsHistograms");
  res.fResultsFlagsPro->Fill(0.5, res.fSaveResultsHistograms);
  // ...
  res.fResultsList->Add(res.fResultsFlagsPro);

  // b) Book results histograms, which in addition act as a sort of "abstract" interface, which defines common binning, etc., for other groups of histograms:
  for (Int_t v = 0; v < eAsFunctionOf_N; v++) {
    if (res.fUseResultsProVariableLengthBins[v]) {
      // per demand, variable-length binning:
      res.fResultsPro[v] = new TProfile(Form("fResultsPro[%s]", res.fResultsProRawName[v].Data()), "...", res.fResultsProVariableLengthBins[v]->GetSize() - 1, res.fResultsProVariableLengthBins[v]->GetArray());
    } else {
      // the default fixed-length binning:
      res.fResultsPro[v] = new TProfile(Form("fResultsPro[%s]", res.fResultsProRawName[v].Data()), "...", (Int_t)res.fResultsProFixedLengthBins[v][0], res.fResultsProFixedLengthBins[v][1], res.fResultsProFixedLengthBins[v][2]);
    }

    // Optionally, save these histograms. Or just use them as an "abstract" interface for the booking of other group of histograms:
    if (res.fSaveResultsHistograms) {
      res.fResultsList->Add(res.fResultsPro[v]);
    }
  } // for (Int_t v = 0; v < eAsFunctionOf_N; v++) {

} // void BookResultsHistograms()

//============================================================

void BookTheRest()
{
  // Here I book everything not sorted (yes) in specific functions above.

  // a) Book the timer;
  // *) ...

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Book the timer:
  if (tc.fUseStopwatch) {
    tc.fTimer[eGlobal] = new TStopwatch();
    tc.fTimer[eGlobal]->Start();
    tc.fTimer[eLocal] = new TStopwatch();
  }

} // void BookTheRest()

//============================================================

template <typename T>
void Preprocess(T const& collision)
{
  // Do all thingies before starting to process data (e.g. count number of events, fetch the run number, get the weights for this run number, etc.).

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // *) If I reached max number of events, ignore the remaining collisions:
  if (MaxNumberOfEvents()) {
    BailOut();
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
  if (!pw.fParticleWeightsAreFetched) {
    if (pw.fUseWeights[wPHI] || pw.fUseWeights[wPT] || pw.fUseWeights[wETA] || pw.fUseDiffWeights[wPHIPT] || pw.fUseDiffWeights[wPHIETA]) {
      GetParticleWeights();
      pw.fParticleWeightsAreFetched = kTRUE;
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
  ebye.fSelectedTracks = 0;
  ebye.fCentrality = 0;

  // b) Q-vectors:
  if (qv.fCalculateQvectors) {
    // b0) generic Q-vector:
    ResetQ();
    // b1) integrated Q-vector:
    for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
      for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) // weight power
      {
        qv.fQvector[h][wp] = TComplex(0., 0.);
      }
    }
    // b2) diff. Q-vector:
    for (Int_t bin = 1; bin <= gMaxNoBinsKine; bin++) {
      qv.fqVectorEntries[PTq][bin - 1] = 0; // TBI 20240214 shall I loop also over enum's PTq and ETAq? If yes, fix it also below for qv.fqvector[PTq][bin - 1][...
      qv.fqVectorEntries[ETAq][bin - 1] = 0;
      for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
        for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
          qv.fqvector[PTq][bin - 1][h][wp] = TComplex(0., 0.);
          qv.fqvector[ETAq][bin - 1][h][wp] = TComplex(0., 0.);
        } // for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
      }   // for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
    }     // for (Int_t b = 0; b < gMaxNoBinsKine; b++ ) {
  }       // if(qv.fCalculateQvectors)

  // c) Reset ebe containers for nested loops:
  if (nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoop) {
    if (nl.ftaNestedLoops[0]) {
      nl.ftaNestedLoops[0]->Reset();
    }
    if (nl.ftaNestedLoops[1]) {
      nl.ftaNestedLoops[1]->Reset();
    }

    // TBI 20220803 port still if(fCalculatePtCorrelations){...} and
    // if(fCalculateEtaCorrelations){...}

  } // if(nl.fCalculateNestedLoops||nl.fCalculateCustomNestedLoop)

  // *) TBI 20240208 add here support for differential nested loops
  /*
  if(nl.fCalculateCustomNestedLoop)
  {
    ftaNestedLoopsKine[PTq][bin-1][0]->AddAt(dPhi,fqVectorEntries[PTq][bin-1]);
    ftaNestedLoopsKine[PTq][bin-1][1]->AddAt(wPhi*wPt*wEta,fqVectorEntries[PTq][bin-1]);
  }
  */

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

  // *) Offline trigger;
  // *) sel7() and sel8(); TBI 20240223 sort out eventualy;
  // *) Specific direct event cuts on info available in reconstructed (and corresponding MC truth simulated);
  // *) Specific direct event cuts on info available only in simulated data;
  // *) Test case.

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // *) Offline trigger:
  //    From documentation: Bypass this check if you analyse MC or continuous Run3 data.
  //    In addition: remember that I can use it only for process cases where I have joined aod::Collisions with aod::EvSels
  if constexpr (rs == eRec_Run2 || rs == eRec_Run1) {
    if (ec.fUseTrigger) {
      if (ec.fTrigger.EqualTo("kINT7")) {
        if (!collision.alias_bit(kINT7)) {
          if (tc.fVerbose) {
            LOGF(info, "\033[1;31m%s collision.alias_bit(kINT7)\033[0m", __FUNCTION__);                                   // just a bare function name
            LOGF(info, "\033[1;31m%s Bypass this check if you analyse MC or continuous Run3 data.\033[0m", __FUNCTION__); // just a bare function name
          }
          return kFALSE;
        }
      }
    }
  }

  // *) sel7() and sel8(); TBI 20240223 sort out eventualy:
  // if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1 || rs == eSim_Run1) {
  if constexpr (rs == eRec_Run2 || rs == eRec_Run1) { // TBI 20240223 use the line above, after I join aod::Collisions with aod::EvSels also for RecSim and Sim cases for Run 2 and Run 1
    if (ec.fUseSel7) {                                // from doc: for Run 2 data and MC
      if (!collision.sel7()) {
        if (tc.fVerbose) {
          LOGF(info, "\033[1;31m%s collision.sel7()\033[0m", __FUNCTION__); // just a bare function name
        }
        return kFALSE;
      }
    }
  }
  // if constexpr (rs == eRec || rs == eRecAndSim || rs == eSim) {
  if constexpr (rs == eRec || rs == eRecAndSim) { // TBI 20240223 use the line above, after I join aod::Collisions with aod::EvSels also for Sim case for Run 3
    if (ec.fUseSel8) {                            // from doc: for Run 3 data and MC
      if (!collision.sel8()) {
        if (tc.fVerbose) {
          LOGF(info, "\033[1;31m%s collision.sel8()\033[0m", __FUNCTION__); // just a bare function name
        }
        return kFALSE;
      }
    }
  }

  // *) Specific direct event cuts on info available in reconstructed ...:
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    //   *) NumberOfEvents: => cut directly in void process( ... )
    //   *) TotalMultiplicity:
    if ((tracks.size() < eh.fEventCuts[eTotalMultiplicity][eMin]) ||
        (tracks.size() > eh.fEventCuts[eTotalMultiplicity][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eTotalMultiplicity\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    //   *) SelectedTracks: => cut directly in void process( ... )

    //   *) Centrtality:
    if ((ebye.fCentrality < eh.fEventCuts[eCentrality][eMin]) || (ebye.fCentrality > eh.fEventCuts[eCentrality][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eCentrality\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    //   *) Vertex_x:
    if ((collision.posX() < eh.fEventCuts[eVertex_x][eMin]) ||
        (collision.posX() > eh.fEventCuts[eVertex_x][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eVertex_x\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    //   *) Vertex_y:
    if ((collision.posY() < eh.fEventCuts[eVertex_y][eMin]) ||
        (collision.posY() > eh.fEventCuts[eVertex_y][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eVertex_y\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    //   *) Vertex_z:
    if ((collision.posZ() < eh.fEventCuts[eVertex_z][eMin]) ||
        (collision.posZ() > eh.fEventCuts[eVertex_z][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eVertex_z\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    //   *) NContributors:
    if ((collision.numContrib() < eh.fEventCuts[eNContributors][eMin]) ||
        (collision.numContrib() > eh.fEventCuts[eNContributors][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eNContributors\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    // TBI 20231106 continue here with other event cuts on reconstructed info

    // ... and corresponding MC truth simulated ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip..."); // TBI 20231106 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this collision
        return kFALSE;
      }

      // TBI 20231106 here I cat cut directly on corresponding MC truth simulated, e.g. on collision.mcCollision().posZ(), if necessary

    } // if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {

  } // if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // *) Specific direct event cuts on info available only in simulated data:
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {
    //   *) Impact parameter:
    if ((collision.impactParameter() < eh.fEventCuts[eImpactParameter][eMin]) ||
        (collision.impactParameter() > eh.fEventCuts[eImpactParameter][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eImpactParameter\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
    // ...
  } // if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

  // *) Test case:
  if constexpr (rs == eTest) {
    // TBI 20240223 for the time being, eTest cuts only on eRec info:
    // A few example cuts.

    //   *) ...
    if ((tracks.size() < eh.fEventCuts[eTotalMultiplicity][eMin]) ||
        (tracks.size() > eh.fEventCuts[eTotalMultiplicity][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eTotalMultiplicity\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }

    //   *) Vertex_z:
    if ((collision.posZ() < eh.fEventCuts[eVertex_z][eMin]) ||
        (collision.posZ() > eh.fEventCuts[eVertex_z][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eVertex_z\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }

    //   *) Centrality:
    if ((ebye.fCentrality < eh.fEventCuts[eCentrality][eMin]) || (ebye.fCentrality > eh.fEventCuts[eCentrality][eMax])) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eCentrality\033[0m", __FUNCTION__); // just a bare function name
      }
      return kFALSE;
    }
  } // if constexpr (rs == eTest) {

  return kTRUE;

} // template <eRecSim rs, typename T1, typename T2> Bool_t EventCuts(T1 const& collision, T2 const& tracks)

//============================================================

template <eRecSim rs, typename T1, typename T2>
void FillEventHistograms(T1 const& collision, T2 const& tracks, eBeforeAfter ba)
{
  // Fill all event histograms for reconstructed or simulated data.

  // a) Fill reconstructed, and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1);
  // b) Fill only simulated (common to Run 3, Run 2 and Run 1);
  // c) Fill reconstructed (Run 3 specific);
  // d) Fill only simulated (Run 3 specific);
  // e) Fill reconstructed (Run 2 specific);
  // f) Fill only simulated (Run 2 specific);
  // g) Fill reconstructed (Run 1 specific);
  // h) Fill only simulated (Run 1 specific);
  // i) Test case.

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s eBeforeAfter = %d \033[0m", __FUNCTION__, (Int_t)ba); // just a bare function name
  }

  // a) Fill reconstructed ... (common to Run 3, Run 2 and Run 1):
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    !eh.fEventHistograms[eNumberOfEvents][eRec][ba] ? true : eh.fEventHistograms[eNumberOfEvents][eRec][ba]->Fill(0.5); // basically, if histogram is not booked, do nothing. 'true' is a placeholder, for the time being
    !eh.fEventHistograms[eVertex_x][eRec][ba] ? true : eh.fEventHistograms[eVertex_x][eRec][ba]->Fill(collision.posX());
    !eh.fEventHistograms[eVertex_y][eRec][ba] ? true : eh.fEventHistograms[eVertex_y][eRec][ba]->Fill(collision.posY());
    !eh.fEventHistograms[eVertex_z][eRec][ba] ? true : eh.fEventHistograms[eVertex_z][eRec][ba]->Fill(collision.posZ());
    !eh.fEventHistograms[eNContributors][eRec][ba] ? true : eh.fEventHistograms[eNContributors][eRec][ba]->Fill(collision.numContrib());
    !eh.fEventHistograms[eTotalMultiplicity][eRec][ba] ? true : eh.fEventHistograms[eTotalMultiplicity][eRec][ba]->Fill(tracks.size());  // TBI 20231106 check and validate further
    !eh.fEventHistograms[eSelectedTracks][eRec][ba] ? true : eh.fEventHistograms[eSelectedTracks][eRec][ba]->Fill(ebye.fSelectedTracks); // TBI 20240108 this one makes sense only for eAfter
    !eh.fEventHistograms[eMultTPC][eRec][ba] ? true : eh.fEventHistograms[eMultTPC][eRec][ba]->Fill(collision.multTPC());
    !eh.fEventHistograms[eMultNTracksPV][eRec][ba] ? true : eh.fEventHistograms[eMultNTracksPV][eRec][ba]->Fill(collision.multNTracksPV());
    !eh.fEventHistograms[eCentrality][eRec][ba] ? true : eh.fEventHistograms[eCentrality][eRec][ba]->Fill(ebye.fCentrality);

    // ... and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1) ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }
      !eh.fEventHistograms[eNumberOfEvents][eSim][ba] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][ba]->Fill(0.5);
      !eh.fEventHistograms[eVertex_x][eSim][ba] ? true : eh.fEventHistograms[eVertex_x][eSim][ba]->Fill(collision.mcCollision().posX());
      !eh.fEventHistograms[eVertex_y][eSim][ba] ? true : eh.fEventHistograms[eVertex_y][eSim][ba]->Fill(collision.mcCollision().posY());
      !eh.fEventHistograms[eVertex_z][eSim][ba] ? true : eh.fEventHistograms[eVertex_z][eSim][ba]->Fill(collision.mcCollision().posZ());
      // eh.fEventHistograms[eTotalMultiplicity][eSim][ba]->Fill(tracks.size()); // TBI 20231106 check how to get corresponding MC truth info, and validate further
      // eh.fEventHistograms[eSelectedTracks][eSim][ba]->Fill(ebye.fSelectedTracks); // TBI 20240108 this one makes sense only for eAfter + re-think if I really need it here
      // TBI 20240120 eMultFT0M, ..., eMultNTracksPV are not needed here
      // eh.fEventHistograms[eCentrality][eSim][ba]->Fill(ebye.fCentrality); // TBI 20240120 this case is still not supported in DetermineCentrality()
    } // if constexpr (rs == eRecAndSim) {
  }   // if constexpr (rs == eRec || rs == eRecAndSim) {

  // b) Fill only simulated (common to Run 3, Run 2 and Run 1):
  if constexpr (rs == eSim) {
    !eh.fEventHistograms[eImpactParameter][eSim][ba] ? true : eh.fEventHistograms[eImpactParameter][eSim][ba]->Fill(collision.impactParameter()); // yes, because in this branch 'collision' is always aod::McCollision
    !eh.fEventHistograms[eSelectedTracks][eSim][ba] ? true : eh.fEventHistograms[eSelectedTracks][eSim][ba]->Fill(ebye.fSelectedTracks);          // TBI 20240108 this one makes sense only for eAfter
    // eh.fEventHistograms[eCentrality][eSim][ba]->Fill(ebye.fCentrality); // TBI 20240120 this case is still not supported in DetermineCentrality()
    // eh.fEventHistograms[eTotalMultiplicity][eSim][ba]->Fill(tracks.size()); // TBI 20231030 check further how to use the same thing for 'sim'
  } // if constexpr (rs == eSim) {

  // -----------------------------------------------------------------------------

  // c) Fill reconstructed (Run 3 specific):
  if constexpr (rs == eRec || rs == eRecAndSim) {
    !eh.fEventHistograms[eMultFT0M][eRec][ba] ? true : eh.fEventHistograms[eMultFT0M][eRec][ba]->Fill(collision.multFT0M());
    !eh.fEventHistograms[eMultFV0M][eRec][ba] ? true : eh.fEventHistograms[eMultFV0M][eRec][ba]->Fill(collision.multFV0M());

    // ... and corresponding MC truth simulated (Run 3 specific1) ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }
      // !eh.fEventHistograms[eNumberOfEvents][eSim][ba] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][ba]->Fill(0.5);
    } // if constexpr (rs == eRecAndSim) {
  }   // if constexpr (rs == eRec || rs == eRecAndSim) {

  // d) Fill only simulated(Run 3 specific):
  if constexpr (rs == eSim) {
    // !eh.fEventHistograms[eImpactParameter][eSim][ba] ? true : eh.fEventHistograms[eImpactParameter][eSim][ba]->Fill(collision.impactParameter()); // yes, because in this branch 'collision' is always aod::McCollision
  } // if constexpr (rs == eSim) {

  // -----------------------------------------------------------------------------

  // e) Fill reconstructed (Run 2 specific):
  if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2) {
    //! eh.fEventHistograms[eMultFT0M][eRec][ba] ? true : eh.fEventHistograms[eMultFT0M][eRec][ba]->Fill(collision.multFT0M());

    // ... and corresponding MC truth simulated (Run 3 specific1) ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim_Run2) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }
      // !eh.fEventHistograms[eNumberOfEvents][eSim][ba] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][ba]->Fill(0.5);
    } // if constexpr (rs == eRecAndSim_Run2) {
  }   // if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2) {

  // f) Fill only simulated(Run 2 specific):
  if constexpr (rs == eSim_Run2) {
    // !eh.fEventHistograms[eImpactParameter][eSim][ba] ? true : eh.fEventHistograms[eImpactParameter][eSim][ba]->Fill(collision.impactParameter()); // yes, because in this branch 'collision' is always aod::McCollision
  } // if constexpr (rs == eSim_Run2) {

  // -----------------------------------------------------------------------------

  // g) Fill reconstructed (Run 1 specific):
  if constexpr (rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    //! eh.fEventHistograms[eMultFT0M][eRec][ba] ? true : eh.fEventHistograms[eMultFT0M][eRec][ba]->Fill(collision.multFT0M());

    // ... and corresponding MC truth simulated (Run 3 specific1) ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim_Run1) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }
      // !eh.fEventHistograms[eNumberOfEvents][eSim][ba] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][ba]->Fill(0.5);
    } // if constexpr (rs == eRecAndSim_Run1) {
  }   // if constexpr (rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // h) Fill only simulated(Run 1 specific):
  if constexpr (rs == eSim_Run1) {
    // !eh.fEventHistograms[eImpactParameter][eSim][ba] ? true : eh.fEventHistograms[eImpactParameter][eSim][ba]->Fill(collision.impactParameter()); // yes, because in this branch 'collision' is always aod::McCollision
  } // if constexpr (rs == eSim_Run1) {

  // -----------------------------------------------------------------------------

  // i) Test case:
  if constexpr (rs == eTest) {
    // TBI 20240223 for the time being, eTest fills only eRec histos:
    // A few example histograms, just to check if I access corresponding tables:
    !eh.fEventHistograms[eVertex_z][eRec][ba] ? true : eh.fEventHistograms[eVertex_z][eRec][ba]->Fill(collision.posZ());
    !eh.fEventHistograms[eTotalMultiplicity][eRec][ba] ? true : eh.fEventHistograms[eTotalMultiplicity][eRec][ba]->Fill(tracks.size());
    !eh.fEventHistograms[eCentrality][eRec][ba] ? true : eh.fEventHistograms[eCentrality][eRec][ba]->Fill(ebye.fCentrality);
  } // if constexpr (rs == eTest) {

} // template <eRecSim rs, typename T1, typename T2> void FillEventHistograms(...)

//============================================================

template <eRecSim rs, typename T>
Bool_t ValidTrack(T const& track)
{
  // Before I start applying any particle tracks, check if this is a valid track.
  // For instance, Run 2 or Run 1 tracklets are NOT valid tracks, as they carry no pt information, and in this function they are filtered out.

  // See enum TrackTypeEnum in O2/Framework/Core/include/Framework/DataTypes.h for further info.

  // a) Validity checks for tracks in Run 3;
  // b) Validity checks for tracks in Run 2 and 1.

  if (tc.fVerboseForEachParticle) {
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__);
  }

  // a) Validity checks for tracks in Run 3;
  // *) Ensure that I am taking into account propagated tracks (and not e.g. track evaluated at innermost update):
  if constexpr (rs == eRec || rs == eRecAndSim) {
    if (!(track.trackType() == o2::aod::track::TrackTypeEnum::Track)) {
      return kFALSE;
    }
  }

  // b) Validity checks for tracks in Run 2 and 1:
  // *) Ensure that tracklets (no pt information) are skipped:
  // if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
  if constexpr (rs == eRec_Run2 || rs == eRec_Run1) { // TBI 20240226 the check below is crashing for eRecAndSim_Run2 , disabling that case here temporarily. Use eventually the above line
    if (!(track.trackType() == o2::aod::track::TrackTypeEnum::Run2Track)) {
      return kFALSE;
    }
  }

  // *) Temporary here, until I cover also these cases:
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {
    LOGF(fatal, "in function \033[1;31m%s at line %d : add support for TrackTypeEnum here also for cases eSim, eSim_Run2 and eSim_Run1\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  // *) In case there are still some pathological cases with one of kine variable being "nan", treat 'em here:
  // TBI 20240228 disable eventually this check, or add a flag to switch it on/off
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    if (isnan(track.phi()) || isnan(track.pt()) || isnan(track.eta())) {
      LOGF(error, "\033[1;33m%s nan hit, skipping this track \033[0m", __PRETTY_FUNCTION__);
      LOGF(error, "dPhi = %f\ndPt = %f\ndEta = %f\ntrackType() = %d", track.phi(), track.pt(), track.eta(), static_cast<int>(track.trackType())); // TBI 20240226 check if I need this cast here
      return kFALSE;
    }
  }

  // *) All checks above survived, then it's a valid track:
  return kTRUE;

} // template <typename T> Bool_t ValidTrack(T const& track)

//============================================================

template <eRecSim rs, typename T>
Bool_t ParticleCuts(T const& track)
{
  // Particles cuts.

  // a) Particle cuts on info available in reconstructed (and corresponding MC truth simulated);
  // b) Particle cuts on info available only in simulated data;
  // c) Test case.

  // TBI 20240213 at the moment, I take that there is nothing specific for Run 3, 2, 1 here. Otherwise, see what I did in EventCuts

  if (tc.fVerboseForEachParticle) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Particle cuts on info available in reconstructed ...:
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    if ((track.pt() < pt_min) || (track.pt() > pt_max)) {
      return kFALSE;
    }

    // TBI 20231107 other cuts on reconstructed track ...

    // ... and corresponding MC truth simulated ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for this track, skip...");
        return kFALSE; // TBI 20231107 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this particle
      }
      auto mcparticle = track.mcParticle();                           // corresponding MC truth simulated particle
      if ((mcparticle.pt() < pt_min) || (mcparticle.pt() > pt_max)) { // TBI 20231107 re-thing if i really cut directly on MC truth, keep it in sync with what I did in AliPhysics
        return kFALSE;
      }

      // TBI 20231107 other cuts on corresponding MC truth particle ...

    } // if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
  }   // if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // b) Particle cuts on info available only in simulated data:
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {
    // Remark: in this branch, 'track' is always TracksSim = aod::McParticles
    if ((track.pt() < pt_min) || (track.pt() > pt_max)) {
      return kFALSE;
    }

    // TBI 20231107 other cuts on simulated ...

  } // if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

  // c) Test case:
  // TBI 2024023 for the time being, eTest cuts only on eRec info.
  if constexpr (rs == eTest) {
    if ((track.pt() < pt_min) || (track.pt() > pt_max)) {
      return kFALSE;
    }
  } // if constexpr (rs == eTest) {

  return kTRUE;

} // template <typename T> Bool_t ParticleCuts(T const& track)

//============================================================

template <eRecSim rs, typename T>
void FillParticleHistograms(T const& track, eBeforeAfter ba)
{
  // Fill all particle histograms for reconstructed and simulated data.

  // a) Fill reconstructed (and corresponding MC truth simulated);
  // b) Fill only simulated;
  // c) Test case.

  if (tc.fVerboseForEachParticle) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Fill reconstructed ...:
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    !ph.fParticleHistograms[ePhi][eRec][ba] ? true : ph.fParticleHistograms[ePhi][eRec][ba]->Fill(track.phi()); // basically, if hist is not booked, do nothing. 'true' is a placeholder, for the time being
    !ph.fParticleHistograms[ePt][eRec][ba] ? true : ph.fParticleHistograms[ePt][eRec][ba]->Fill(track.pt());
    !ph.fParticleHistograms[eEta][eRec][ba] ? true : ph.fParticleHistograms[eEta][eRec][ba]->Fill(track.eta());

    // ... and corresponding MC truth simulated ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for this track, skip...");
        return;
      }
      auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle
      !ph.fParticleHistograms[ePhi][eSim][ba] ? true : ph.fParticleHistograms[ePhi][eSim][ba]->Fill(mcparticle.phi());
      !ph.fParticleHistograms[ePt][eSim][ba] ? true : ph.fParticleHistograms[ePt][eSim][ba]->Fill(mcparticle.pt());
      !ph.fParticleHistograms[eEta][eSim][ba] ? true : ph.fParticleHistograms[eEta][eSim][ba]->Fill(mcparticle.eta());
      !ph.fParticleHistograms[ePDG][eSim][ba] ? true : ph.fParticleHistograms[ePDG][eSim][ba]->Fill(mcparticle.pdgCode());
    } // if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
  }   // if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // b) Fill only simulated:
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {
    // Remark: in this branch, 'track' is always TracksSim = aod::McParticles
    !ph.fParticleHistograms[ePhi][eSim][ba] ? true : ph.fParticleHistograms[ePhi][eSim][ba]->Fill(track.phi());
    !ph.fParticleHistograms[ePt][eSim][ba] ? true : ph.fParticleHistograms[ePt][eSim][ba]->Fill(track.pt());
    !ph.fParticleHistograms[eEta][eSim][ba] ? true : ph.fParticleHistograms[eEta][eSim][ba]->Fill(track.eta());
    !ph.fParticleHistograms[ePDG][eSim][ba] ? true : ph.fParticleHistograms[ePDG][eSim][ba]->Fill(track.pdgCode());
  } // if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

  /* TBI 20231019 use also these + check further
  // From aod::TracksExtra
  ph.fParticleHistograms[etpcNClsCrossedRows][rs][ba]->Fill(track.tpcNClsCrossedRows());

  // From aod::TracksDCA
  ph.fParticleHistograms[eDCA_xy][rs][ba]->Fill(track.dcaXY());
  ph.fParticleHistograms[eDCA_z][rs][ba]->Fill(track.dcaZ());
  */

  // c) Test case:
  if constexpr (rs == eTest) {
    // TBI 20240223 for the time being, eTest fills eRec histos:
    !ph.fParticleHistograms[ePhi][eRec][ba] ? true : ph.fParticleHistograms[ePhi][eRec][ba]->Fill(track.phi());
    !ph.fParticleHistograms[ePt][eRec][ba] ? true : ph.fParticleHistograms[ePt][eRec][ba]->Fill(track.pt());
    !ph.fParticleHistograms[eEta][eRec][ba] ? true : ph.fParticleHistograms[eEta][eRec][ba]->Fill(track.eta());
  } // if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

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
      qv.fQ[h][wp] = qv.fQvector[h][wp];
    }
  }

  // b) Calculate correlations:
  for (Int_t h = 1; h <= gMaxHarmonic; h++) // harmonic
  {
    // 2p:
    if (ebye.fSelectedTracks < 2) {
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
      LOGF(fatal, "In function \033[1;31m%s at line %d, wTwo = %f <=0. ebye.fSelectedTracks = %d\033[0m", __PRETTY_FUNCTION__, __LINE__, wTwo, ebye.fSelectedTracks);
    }

    if (nl.fCalculateCustomNestedLoop) {
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
    } // if(nl.fCalculateCustomNestedLoop)

    /*
      // for on-the-fly and internal validation, rescale results with theoretical value: if(fCalculateOnTheFly && fOnTheFlyFlowAmplitudes &&
      fRescaleWithTheoreticalInput &&
         TMath::Abs(fOnTheFlyFlowAmplitudes->GetAt(h-1))>0.){twoC/=pow(fOnTheFlyFlowAmplitudes->GetAt(h-1),2.);}
      else if(fUseInternalValidation && fInternalValidationAmplitudes &&
      fRescaleWithTheoreticalInput &&
              TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){twoC/=pow(fInternalValidationAmplitudes->GetAt(h-1),2.);}
    */

    // integrated:
    if (mupa.fCorrelationsPro[0][h - 1][AFO_INTEGRATED]) {
      mupa.fCorrelationsPro[0][h - 1][AFO_INTEGRATED]->Fill(0.5, twoC, wTwo);
    }
    // vs. multiplicity:
    if (mupa.fCorrelationsPro[0][h - 1][AFO_MULTIPLICITY]) {
      mupa.fCorrelationsPro[0][h - 1][AFO_MULTIPLICITY]->Fill(ebye.fSelectedTracks + 0.5, twoC, wTwo);
    }
    // vs. centrality:
    if (mupa.fCorrelationsPro[0][h - 1][AFO_CENTRALITY]) {
      mupa.fCorrelationsPro[0][h - 1][AFO_CENTRALITY]->Fill(ebye.fCentrality, twoC, wTwo);
    }

    // 4p:
    if (ebye.fSelectedTracks < 4) {
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
      LOGF(fatal, "In function \033[1;31m%s at line %d, wFour = %f <=0. ebye.fSelectedTracks = %d\033[0m", __PRETTY_FUNCTION__, __LINE__, wFour, ebye.fSelectedTracks);
      // TBI 20240110 shall I 'continue' here, instead of bailing out?
    }

    if (nl.fCalculateCustomNestedLoop) {
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
    } // if(nl.fCalculateCustomNestedLoop)

    //    if(fUseInternalValidation && fInternalValidationAmplitudes && fRescaleWithTheoreticalInput &&
    //       TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){fourC/=pow(fInternalValidationAmplitudes->GetAt(h-1),4.);}

    // integrated:
    if (mupa.fCorrelationsPro[1][h - 1][AFO_INTEGRATED]) {
      mupa.fCorrelationsPro[1][h - 1][AFO_INTEGRATED]->Fill(0.5, fourC, wFour);
    }
    // vs. multiplicity:
    if (mupa.fCorrelationsPro[1][h - 1][AFO_MULTIPLICITY]) {
      mupa.fCorrelationsPro[1][h - 1][AFO_MULTIPLICITY]->Fill(ebye.fSelectedTracks + 0.5, fourC, wFour);
    }
    // vs. centrality:
    if (mupa.fCorrelationsPro[1][h - 1][AFO_CENTRALITY]) {
      mupa.fCorrelationsPro[1][h - 1][AFO_CENTRALITY]->Fill(ebye.fCentrality, fourC, wFour);
    }

    // 6p:
    if (ebye.fSelectedTracks < 6) {
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
      LOGF(fatal, "In function \033[1;31m%s at line %d, wSix = %f <=0. ebye.fSelectedTracks = %d\033[0m", __PRETTY_FUNCTION__, __LINE__, wSix, ebye.fSelectedTracks);
      // TBI 20240110 shall I 'continue' here, instead of bailing out?
    }

    if (nl.fCalculateCustomNestedLoop) {
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
    } // if(nl.fCalculateCustomNestedLoop)

    //    if(fUseInternalValidation && fInternalValidationAmplitudes && fRescaleWithTheoreticalInput &&
    //       TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){sixC/=pow(fInternalValidationAmplitudes->GetAt(h-1),4.);}

    // integrated:
    if (mupa.fCorrelationsPro[2][h - 1][AFO_INTEGRATED]) {
      mupa.fCorrelationsPro[2][h - 1][AFO_INTEGRATED]->Fill(0.5, sixC, wSix);
    }
    // vs. multiplicity:
    if (mupa.fCorrelationsPro[2][h - 1][AFO_MULTIPLICITY]) {
      mupa.fCorrelationsPro[2][h - 1][AFO_MULTIPLICITY]->Fill(ebye.fSelectedTracks + 0.5, sixC, wSix);
    }
    // vs. centrality:
    if (mupa.fCorrelationsPro[2][h - 1][AFO_CENTRALITY]) {
      mupa.fCorrelationsPro[2][h - 1][AFO_CENTRALITY]->Fill(ebye.fCentrality, sixC, wSix);
    }

    // 8p:
    if (ebye.fSelectedTracks < 8) {
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
      LOGF(fatal, "In function \033[1;31m%s at line %d, wEight = %f <=0. ebye.fSelectedTracks = %d\033[0m", __PRETTY_FUNCTION__, __LINE__, wEight, ebye.fSelectedTracks);
      // TBI 20240110 shall I 'continue' here, instead of bailing out?
    }

    if (nl.fCalculateCustomNestedLoop) {
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
    } // if(nl.fCalculateCustomNestedLoop)

    //    if(fUseInternalValidation && fInternalValidationAmplitudes && fRescaleWithTheoreticalInput &&
    //       TMath::Abs(fInternalValidationAmplitudes->GetAt(h-1))>0.){eightC/=pow(fInternalValidationAmplitudes->GetAt(h-1),4.);}

    // integrated:
    if (mupa.fCorrelationsPro[3][h - 1][AFO_INTEGRATED]) {
      mupa.fCorrelationsPro[3][h - 1][AFO_INTEGRATED]->Fill(0.5, eightC, wEight);
    }
    // vs. multiplicity:
    if (mupa.fCorrelationsPro[3][h - 1][AFO_MULTIPLICITY]) {
      mupa.fCorrelationsPro[3][h - 1][AFO_MULTIPLICITY]->Fill(ebye.fSelectedTracks + 0.5, eightC, wEight);
    }
    // vs. centrality:
    if (mupa.fCorrelationsPro[3][h - 1][AFO_CENTRALITY]) {
      mupa.fCorrelationsPro[3][h - 1][AFO_CENTRALITY]->Fill(ebye.fCentrality, eightC, wEight);
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
      qv.fQ[h][wp] = qv.fQvector[h][wp];
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
      if (!t0.fTest0Labels[mo][mi]) // I do not stream them.
      {
        for (Int_t v = 0; v < eAsFunctionOf_N; v++) {
          if (t0.fTest0Pro[mo][mi][v]) {
            t0.fTest0Labels[mo][mi] = new TString(t0.fTest0Pro[mo][mi][v]->GetTitle()); // there is no memory leak here, since this is executed only once due to if(!fTest0Labels[mo][mi])
            break;                                                                      // yes, since for all v they are the same, so I just need to fetch it from one
          }
        }
      } // if(!t0_afTest0Labels[mo][mi])

      if (t0.fTest0Labels[mo][mi]) {
        // Extract harmonics from TString, FS is " ":
        for (Int_t h = 0; h <= mo; h++) {
          TObjArray* oa = t0.fTest0Labels[mo][mi]->Tokenize(" ");
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
            if (ebye.fSelectedTracks < 1) {
              return;
            }
            correlation = One(n[0]).Re();
            weight = One(0).Re();
            break;

          case 2:
            if (ebye.fSelectedTracks < 2) {
              return;
            }
            correlation = Two(n[0], n[1]).Re();
            weight = Two(0, 0).Re();
            break;

          case 3:
            if (ebye.fSelectedTracks < 3) {
              return;
            }
            correlation = Three(n[0], n[1], n[2]).Re();
            weight = Three(0, 0, 0).Re();
            break;

          case 4:
            if (ebye.fSelectedTracks < 4) {
              return;
            }
            correlation = Four(n[0], n[1], n[2], n[3]).Re();
            weight = Four(0, 0, 0, 0).Re();
            break;

          case 5:
            if (ebye.fSelectedTracks < 5) {
              return;
            }
            correlation = Five(n[0], n[1], n[2], n[3], n[4]).Re();
            weight = Five(0, 0, 0, 0, 0).Re();
            break;

          case 6:
            if (ebye.fSelectedTracks < 6) {
              return;
            }
            correlation = Six(n[0], n[1], n[2], n[3], n[4], n[5]).Re();
            weight = Six(0, 0, 0, 0, 0, 0).Re();
            break;

          case 7:
            if (ebye.fSelectedTracks < 7) {
              return;
            }
            correlation = Seven(n[0], n[1], n[2], n[3], n[4], n[5], n[6]).Re();
            weight = Seven(0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 8:
            if (ebye.fSelectedTracks < 8) {
              return;
            }
            correlation = Eight(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7]).Re();
            weight = Eight(0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 9:
            if (ebye.fSelectedTracks < 9) {
              return;
            }
            correlation = Nine(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8]).Re();
            weight = Nine(0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 10:
            if (ebye.fSelectedTracks < 10) {
              return;
            }
            correlation = Ten(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9]).Re();
            weight = Ten(0, 0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 11:
            if (ebye.fSelectedTracks < 11) {
              return;
            }
            correlation = Eleven(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10]).Re();
            weight = Eleven(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          case 12:
            if (ebye.fSelectedTracks < 12) {
              return;
            }
            correlation = Twelve(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11]).Re();
            weight = Twelve(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
            break;

          default:
            LOGF(fatal, "in function \033[1;31m%s at line %d\n Not supported yet: %s \n\n\033[0m", __PRETTY_FUNCTION__, __LINE__, t0.fTest0Labels[mo][mi]->Data());
        } // switch(mo+1)

        // Insanity check on weight:
        if (!(weight > 0.)) {
          LOGF(fatal, "in function \033[1;31m%s at line %d\n Is perhaps order of correlator bigger than the number of particles? %s \n\n\033[0m", __PRETTY_FUNCTION__, __LINE__, t0.fTest0Labels[mo][mi]->Data());
        }

        // e-b-e sanity check:
        if (nl.fCalculateCustomNestedLoop) {
          TArrayI* harmonics = new TArrayI(mo + 1);
          for (Int_t i = 0; i < mo + 1; i++) {
            harmonics->SetAt(n[i], i);
          }
          Double_t nestedLoopValue = this->CalculateCustomNestedLoop(harmonics);
          if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(correlation / weight - nestedLoopValue) > 1.e-5) {
            LOGF(fatal, "in function \033[1;31m%s at line %d, nestedLoopValue = %f is not the same as correlation/weight = %f, for correlator %s\033[0m", __PRETTY_FUNCTION__, __LINE__, nestedLoopValue, correlation / weight, t0.fTest0Labels[mo][mi]->Data());
          } else {
            LOGF(info, "=> e-b-e check with CustomNestedLoop is OK for %d-p Test0 corr. %s", mo + 1, t0.fTest0Labels[mo][mi]->Data());
          }
          delete harmonics;
          harmonics = NULL;
        } // if(nl.fCalculateCustomNestedLoop)

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
        if (t0.fTest0Pro[mo][mi][AFO_INTEGRATED]) {
          t0.fTest0Pro[mo][mi][AFO_INTEGRATED]->Fill(0.5, correlation / weight, weight);
        }
        // vs. multiplicity:
        if (t0.fTest0Pro[mo][mi][AFO_MULTIPLICITY]) {
          t0.fTest0Pro[mo][mi][AFO_MULTIPLICITY]->Fill(ebye.fSelectedTracks + 0.5, correlation / weight, weight);
        }
        // vs. centrality:
        if (t0.fTest0Pro[mo][mi][AFO_CENTRALITY]) {
          t0.fTest0Pro[mo][mi][AFO_CENTRALITY]->Fill(ebye.fCentrality, correlation / weight, weight);
        }
      } // if(t0.fTest0Labels[mo][mi])
    }   // for(Int_t mi=0;mi<gMaxIndex;mi++)
  }     // for(Int_t mo=0;mo<gMaxCorrelator;mo++)

  // c) Flush the generic Q-vectors:
  ResetQ();

} // void CalculateTest0()

//============================================================

void CalculateKineTest0(const char* kc)
{
  // Calculate analytically kine Test0 from Q-vectors.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // *) ...
  Int_t kb = -1;         // which kine bin
  Int_t qvKine_var = -1; // which eqvectorKine enum
  Int_t nBins = -1;
  if (TString(kc).EqualTo("pt")) {
    kb = AFO_PT;
    qvKine_var = PTq;
    // nBins = fKinematicsBins[PT][0];
    nBins = res.fResultsPro[AFO_PT]->GetNbinsX();
  } else if (TString(kc).EqualTo("eta")) {
    kb = AFO_ETA;
    qvKine_var = ETAq;
    // nBins = fKinematicsBins[ETA][0];
    nBins = res.fResultsPro[AFO_ETA]->GetNbinsX();
  }

  // *) Uniform loop over bin for all kine variables:
  for (Int_t b = 0; b < nBins; b++) {
    // Ensures that in each bin of interest, I have the same cut on number of particles, like in integrated analysis:
    /* TBI 20240221 enable eventually
        if(fUseSelectedTracksCuts) {
          if(fqVectorEntries[qvKine_var][b] < fSelectedTracksCuts[0]){continue;}
          if(fqVectorEntries[qvKine_var][b] > fSelectedTracksCuts[1]){continue;}
        }
    */

    // *) Re-initialize Q-vector to be q-vector in this bin:
    // After that, I can call all standard Q-vector functions again:
    for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
      for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
        qv.fQ[h][wp] = qv.fqvector[qvKine_var][b][h][wp];
      }
    }

    // *) Okay, let's do the differential calculus:
    Double_t correlation = 0.;
    Double_t weight = 0.;
    Int_t n[gMaxCorrelator] = {0}; // array holding harmonics

    for (Int_t mo = 0; mo < gMaxCorrelator; mo++) {
      for (Int_t mi = 0; mi < gMaxIndex; mi++) {
        // TBI 20240221 I do not have to loop each time all the way up to gMaxCorrelator and gMaxIndex, but nevermind now, it's not a big efficiency loss.
        if (t0.fTest0Labels[mo][mi]) {
          // Extract harmonics from TString, FS is " ":
          for (Int_t h = 0; h <= mo; h++) {
            // cout<<Form("h = %d, t0.fTest0Labels[%d][%d] = ",h,mo,mi)<<t0.fTest0Labels[mo][mi]->Data()<<endl;
            TObjArray* oa = t0.fTest0Labels[mo][mi]->Tokenize(" ");
            if (!oa) {
              cout << __LINE__ << endl;
              exit(1);
            }
            n[h] = TString(oa->At(h)->GetName()).Atoi();
            delete oa; // yes, otherwise it's a memory leak
          }

          if (qv.fqVectorEntries[qvKine_var][b] < mo + 1) {
            continue;
          }

          switch (mo + 1) // which order? yes, mo+1
          {
            case 1:
              correlation = One(n[0]).Re();
              weight = One(0).Re();
              break;

            case 2:
              correlation = Two(n[0], n[1]).Re();
              weight = Two(0, 0).Re();
              break;

            case 3:
              correlation = Three(n[0], n[1], n[2]).Re();
              weight = Three(0, 0, 0).Re();
              break;

            case 4:
              correlation = Four(n[0], n[1], n[2], n[3]).Re();
              weight = Four(0, 0, 0, 0).Re();
              break;

            case 5:
              correlation = Five(n[0], n[1], n[2], n[3], n[4]).Re();
              weight = Five(0, 0, 0, 0, 0).Re();
              break;

            case 6:
              correlation = Six(n[0], n[1], n[2], n[3], n[4], n[5]).Re();
              weight = Six(0, 0, 0, 0, 0, 0).Re();
              break;

            case 7:
              correlation = Seven(n[0], n[1], n[2], n[3], n[4], n[5], n[6]).Re();
              weight = Seven(0, 0, 0, 0, 0, 0, 0).Re();
              break;

            case 8:
              correlation = Eight(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7]).Re();
              weight = Eight(0, 0, 0, 0, 0, 0, 0, 0).Re();
              break;

            case 9:
              correlation = Nine(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8]).Re();
              weight = Nine(0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
              break;

            case 10:
              correlation = Ten(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9]).Re();
              weight = Ten(0, 0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
              break;

            case 11:
              correlation = Eleven(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10]).Re();
              weight = Eleven(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
              break;

            case 12:
              correlation = Twelve(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11]).Re();
              weight = Twelve(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).Re();
              break;

            default:
              LOGF(fatal, "in function \033[1;31m%s at line %d\n Not supported yet: %s \n\n\033[0m", __PRETTY_FUNCTION__, __LINE__, t0.fTest0Labels[mo][mi]->Data());
          } // switch(mo+1)

          /*
               // e-b-e sanity check:
               if(fCalculateCustomNestedLoop)
               {
                TArrayI *harmonics = new TArrayI(mo+1);
                for(Int_t i=0;i<mo+1;i++)
                {
                 harmonics->SetAt(n[i],i);
                }
                if(!(weight>0.))
                {
                 cout<<fTest0Labels[mo][mi]->Data()<<endl;
                 Red("Is perhaps order of some requested correlator bigger than the number of particles?");
                 cout<<__LINE__<<endl; exit(1);
                }
                Double_t nestedLoopValue = this->CalculateKineCustomNestedLoop(harmonics,kc,b);
                if(TMath::Abs(nestedLoopValue)>0. && TMath::Abs(correlation/weight - nestedLoopValue)>1.e-5)
                {
                 cout<<fTest0Labels[mo][mi]->Data()<<endl;
                 cout<<"correlation: "<<correlation/weight<<endl;
                 cout<<"custom loop: "<<this->CalculateCustomNestedLoop(harmonics)<<endl;
                 cout<<__LINE__<<endl; exit(1);
                }
                else
                {
                 cout<<Form("=> e-b-e check with CustomNestedLoop is OK for %d-p Test0 corr. %s, bin = %d",mo+1,fTest0Labels[mo][mi]->Data(),b+1)<<endl;
                }
                delete harmonics; harmonics = NULL;
               } // if(fCalculateCustomNestedLoop)
          */

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

          /*
               // Insanity check for the event weight:
               if(!(weight > 0.))
               {
                // If it's negative, that means that sum of particle weights is smaller than "number of particles - 1"
                // In that case, you can simply rescale all particle weights, so that each of them is > 1, basically recalculate weights.root files with such a rescaling.
                cout<<Form("b = %d",b)<<endl;
                cout<<Form("qvKine_var = %d",qvKine_var)<<endl;
                cout<<Form("event weight = %e",weight)<<endl;
                cout<<Form("sum of particle weights = %e",One(0).Re())<<endl;
                cout<<Form("correlation = %f",correlation)<<endl;
                cout<<Form("fTest0Pro[mo][mi][kb]->GetTitle() = %s",fTest0Pro[mo][mi][kb]->GetTitle())<<endl;
                cout<<Form("[mo][mi][kb] = [%d][%d][%d]",mo,mi,kb)<<endl;
                cout<<Form("fSelectedTracks = %d",fSelectedTracks)<<endl;
                cout<<Form("fqVectorEntries[qvKine_var][b] = %d",fqVectorEntries[qvKine_var][b])<<endl;
                cout<<__LINE__<<endl;exit(1);
               }
          */

          // Finally, fill:
          if (t0.fTest0Pro[mo][mi][kb]) {
            t0.fTest0Pro[mo][mi][kb]->Fill(t0.fTest0Pro[mo][mi][kb]->GetXaxis()->GetBinCenter(b + 1), correlation / weight, weight);
          } // fill in the bin center

        } // if(fTest0Labels[mo][mi])
      }   // for(Int_t mi=0;mi<gMaxIndex;mi++)
    }     // for(Int_t mo=0;mo<gMaxCorrelator;mo++)

  } // for(Int_t b=0;b<nBins;b++)

} // CalculateKineTest0(const char* kc)

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

  LOGF(info, "\033[1;32m ebye.fSelectedTracks = %d\033[0m", ebye.fSelectedTracks);
  Int_t nParticles = ebye.fSelectedTracks;

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
    Double_t dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      Double_t dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
      for (int h = 0; h < gMaxHarmonic; h++) {
        // fill cos, 2p, integreated:
        if (nl.fNestedLoopsPro[0][h][AFO_INTEGRATED]) {
          nl.fNestedLoopsPro[0][h][AFO_INTEGRATED]->Fill(
            0.5, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
        // fill cos, 2p, vs. multiplicity:
        if (nl.fNestedLoopsPro[0][h][AFO_MULTIPLICITY]) {
          nl.fNestedLoopsPro[0][h][AFO_MULTIPLICITY]->Fill(
            ebye.fSelectedTracks + 0.5, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)),
            dW1 * dW2);
        }
        // fill cos, 2p, vs. centrality:
        if (nl.fNestedLoopsPro[0][h][AFO_CENTRALITY]) {
          nl.fNestedLoopsPro[0][h][AFO_CENTRALITY]->Fill(
            ebye.fCentrality, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
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
    Double_t dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      Double_t dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        Double_t dPhi3 = nl.ftaNestedLoops[0]->GetAt(i3);
        Double_t dW3 = nl.ftaNestedLoops[1]->GetAt(i3);
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          Double_t dPhi4 = nl.ftaNestedLoops[0]->GetAt(i4);
          Double_t dW4 = nl.ftaNestedLoops[1]->GetAt(i4);
          for (int h = 0; h < gMaxHarmonic; h++) {
            // fill cos, 4p, integreated:
            if (nl.fNestedLoopsPro[1][h][AFO_INTEGRATED]) {
              nl.fNestedLoopsPro[1][h][AFO_INTEGRATED]->Fill(0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. M:
            if (nl.fNestedLoopsPro[1][h][AFO_MULTIPLICITY]) {
              nl.fNestedLoopsPro[1][h][AFO_MULTIPLICITY]->Fill(ebye.fSelectedTracks + 0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. centrality:
            if (nl.fNestedLoopsPro[1][h][AFO_CENTRALITY]) {
              nl.fNestedLoopsPro[1][h][AFO_CENTRALITY]->Fill(ebye.fCentrality, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
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
    Double_t dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      Double_t dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        Double_t dPhi3 = nl.ftaNestedLoops[0]->GetAt(i3);
        Double_t dW3 = nl.ftaNestedLoops[1]->GetAt(i3);
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          Double_t dPhi4 = nl.ftaNestedLoops[0]->GetAt(i4);
          Double_t dW4 = nl.ftaNestedLoops[1]->GetAt(i4);
          for (int i5 = 0; i5 < nParticles; i5++) {
            if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4) {
              continue;
            }
            Double_t dPhi5 = nl.ftaNestedLoops[0]->GetAt(i5);
            Double_t dW5 = nl.ftaNestedLoops[1]->GetAt(i5);
            for (int i6 = 0; i6 < nParticles; i6++) {
              if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5) {
                continue;
              }
              Double_t dPhi6 = nl.ftaNestedLoops[0]->GetAt(i6);
              Double_t dW6 = nl.ftaNestedLoops[1]->GetAt(i6);
              for (int h = 0; h < gMaxHarmonic; h++) {
                // fill cos, 6p, integreated:
                if (nl.fNestedLoopsPro[2][h][AFO_INTEGRATED]) {
                  nl.fNestedLoopsPro[2][h][AFO_INTEGRATED]->Fill(0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. M:
                if (nl.fNestedLoopsPro[2][h][AFO_MULTIPLICITY]) {
                  nl.fNestedLoopsPro[2][h][AFO_MULTIPLICITY]->Fill(ebye.fSelectedTracks + 0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. M:
                if (nl.fNestedLoopsPro[2][h][AFO_CENTRALITY]) {
                  nl.fNestedLoopsPro[2][h][AFO_CENTRALITY]->Fill(ebye.fCentrality, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
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
    Double_t dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      Double_t dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        Double_t dPhi3 = nl.ftaNestedLoops[0]->GetAt(i3);
        Double_t dW3 = nl.ftaNestedLoops[1]->GetAt(i3);
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          Double_t dPhi4 = nl.ftaNestedLoops[0]->GetAt(i4);
          Double_t dW4 = nl.ftaNestedLoops[1]->GetAt(i4);
          for (int i5 = 0; i5 < nParticles; i5++) {
            if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4) {
              continue;
            }
            Double_t dPhi5 = nl.ftaNestedLoops[0]->GetAt(i5);
            Double_t dW5 = nl.ftaNestedLoops[1]->GetAt(i5);
            for (int i6 = 0; i6 < nParticles; i6++) {
              if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5) {
                continue;
              }
              Double_t dPhi6 = nl.ftaNestedLoops[0]->GetAt(i6);
              Double_t dW6 = nl.ftaNestedLoops[1]->GetAt(i6);
              for (int i7 = 0; i7 < nParticles; i7++) {
                if (i7 == i1 || i7 == i2 || i7 == i3 || i7 == i4 || i7 == i5 || i7 == i6) {
                  continue;
                }
                Double_t dPhi7 = nl.ftaNestedLoops[0]->GetAt(i7);
                Double_t dW7 = nl.ftaNestedLoops[1]->GetAt(i7);
                for (int i8 = 0; i8 < nParticles; i8++) {
                  if (i8 == i1 || i8 == i2 || i8 == i3 || i8 == i4 || i8 == i5 || i8 == i6 || i8 == i7) {
                    continue;
                  }
                  Double_t dPhi8 = nl.ftaNestedLoops[0]->GetAt(i8);
                  Double_t dW8 = nl.ftaNestedLoops[1]->GetAt(i8);
                  for (int h = 0; h < gMaxHarmonic; h++) {
                    // fill cos, 8p, integreated:
                    if (nl.fNestedLoopsPro[3][h][AFO_INTEGRATED]) {
                      nl.fNestedLoopsPro[3][h][AFO_INTEGRATED]->Fill(0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. M:
                    if (nl.fNestedLoopsPro[3][h][AFO_MULTIPLICITY]) {
                      nl.fNestedLoopsPro[3][h][AFO_MULTIPLICITY]->Fill(ebye.fSelectedTracks + 0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. M:
                    if (nl.fNestedLoopsPro[3][h][AFO_CENTRALITY]) {
                      nl.fNestedLoopsPro[3][h][AFO_CENTRALITY]->Fill(ebye.fCentrality, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
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
    nBinsQV = mupa.fCorrelationsPro[0][0][v]->GetNbinsX();
    nBinsNL = nl.fNestedLoopsPro[0][0][v]->GetNbinsX();
    if (nBinsQV != nBinsNL) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }
    LOGF(info, "\033[1;32m   [%d] : %s\033[0m", v, res.fResultsProXaxisTitle[v].Data());
    for (Int_t o = 0; o < 4; o++) {
      LOGF(info, "\033[1;32m   ==== <<%d>>-particle correlations ====\033[0m", 2 * (o + 1));
      for (Int_t h = 0; h < gMaxHarmonic; h++) {
        for (Int_t b = 1; b <= nBinsQV; b++) {
          if (mupa.fCorrelationsPro[o][h][v]) {
            valueQV = mupa.fCorrelationsPro[o][h][v]->GetBinContent(b);
          }
          if (nl.fNestedLoopsPro[o][h][v]) {
            valueNL = nl.fNestedLoopsPro[o][h][v]->GetBinContent(b);
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
    return qv.fQ[n][wp];
  }
  return TComplex::Conjugate(qv.fQ[-n][wp]);

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
      qv.fQ[h][wp] = TComplex(0., 0.);
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
  pw.fWeightsHist[ppe] = reinterpret_cast<TH1D*>(hist->Clone());

  if (!pw.fWeightsHist[ppe]) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  // Cosmetics: TBI 20240216 do I really want to overwrite initial cosmetics, perhaps this shall go better into MakeWeights.C ?
  //                         Or I could move all this to GetHistogramWithWeights, where in any case I am setting e.g. histogram title, etc.
  TString sVariable[eWeights_N] = {"#varphi", "p_{t}", "#eta"}; // [phi,pt,eta]
  TString sWeights[eWeights_N] = {"w_{#varphi}", "w_{p_{t}}", "w_{#eta}"};
  pw.fWeightsHist[ppe]->SetStats(kFALSE);
  pw.fWeightsHist[ppe]->GetXaxis()->SetTitle(sVariable[ppe].Data());
  pw.fWeightsHist[ppe]->GetYaxis()->SetTitle(sWeights[ppe].Data());
  pw.fWeightsHist[ppe]->SetFillColor(eFillColor);
  pw.fWeightsHist[ppe]->SetLineColor(eColor);
  pw.fWeightsList->Add(pw.fWeightsHist[ppe]);

  // Flag:
  pw.fUseWeights[ppe] = kTRUE;

} // void SetWeightsHist(TH1D* const hist, const char *variable)

//============================================================

void SetDiffWeightsHist(TH1D* const hist, const char* variable, Int_t bin)
{
  // Copy histogram holding differential weights from an external file to the corresponding data member.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // Basic protection:
  if (!(TString(variable).EqualTo("phipt") || TString(variable).EqualTo("phieta"))) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  Int_t ppe = -1; // TBI 20240215 use enum's instead
  if (TString(variable).EqualTo("phipt")) {
    ppe = 0;
  }
  if (TString(variable).EqualTo("phieta")) {
    ppe = 1;
  }

  // Finally:
  hist->SetDirectory(0);
  pw.fDiffWeightsHist[ppe][bin] = reinterpret_cast<TH1D*>(hist->Clone());

  if (!pw.fDiffWeightsHist[ppe][bin]) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  // Cosmetics: TBI 20240216 do I really want to overwrite initial cosmetics, perhaps this shall go better into MakeWeights.C ?
  //                         Or I could move all this to GetHistogramWithWeights, where in any case I am setting e.g. histogram title, etc.
  TString sVariable[eWeights_N] = {"#varphi", "p_{t}", "#eta"}; // [phi,pt,eta]
  TString sWeights[eWeights_N] = {"w_{#varphi}", "w_{p_{t}}", "w_{#eta}"};
  pw.fDiffWeightsHist[ppe][bin]->SetStats(kFALSE);
  pw.fDiffWeightsHist[ppe][bin]->GetXaxis()->SetTitle(sVariable[ppe].Data());
  pw.fDiffWeightsHist[ppe][bin]->GetYaxis()->SetTitle(sWeights[ppe].Data());
  pw.fDiffWeightsHist[ppe][bin]->SetFillColor(eFillColor);
  pw.fDiffWeightsHist[ppe][bin]->SetLineColor(eColor);
  pw.fWeightsList->Add(pw.fDiffWeightsHist[ppe][bin]);

  // Flag:
  if (!pw.fUseDiffWeights[ppe]) // yes, set it only once to kTRUE, for all bins
  {
    pw.fUseDiffWeights[ppe] = kTRUE;
  }

} // SetDiffWeightsHist(TH1D* const hist, const char *variable, Int_t bin)

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
  return pw.fWeightsHist[ppe];

} // TH1D* GetWeightsHist(const char *variable)

//============================================================

TH1D* GetHistogramWithWeights(const char* filePath, const char* runNumber, const char* variable, Int_t bin = -1)
{
  // Get and return histogram with weights from an external file.
  // If bin > 0, differential weights for that bin are searched for.
  // If bin = -1, integrated weights are searched for, i.e. in this case "bin" variable has no effect.
  // I do it this way, so as to condense GetHistogramWithWeights(...) and GetHistogramWithDiffWeights(...) from MuPa class in
  // one routine here, so that I do not duplicate code related to CCDB access, etc.

  // a) Return value;
  // b) Basic protection for arguments;
  // c) Determine from filePath if the file in on a local machine, or in AliEn, or in CCDB;
  // d) Handle the AliEn case;
  // e) Handle the CCDB case;
  // f) Handle the local case;
  // g) The final touch on histogram with weights.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
    LOGF(info, "\033[1;33m filePath = %s\033[0m", filePath);
    LOGF(info, "\033[1;33m runNumber = %s\033[0m", runNumber);
    LOGF(info, "\033[1;33m variable = %s\033[0m", variable);
    LOGF(info, "\033[1;33m bin = %d\033[0m", bin);
    LOGF(info, "\033[1;33m fTaskName = %s\033[0m", tc.fTaskName.Data());
  }

  // a) Return value:
  TH1D* hist = NULL;
  TList* baseList = NULL;     // base top-level list in the TFile, e.g. named "ccdb_object"
  TList* listWithRuns = NULL; // nested list with run-wise TList's holding run-specific weights

  // b) Basic protection for arguments:
  //    Remark: below I do one more specific check.
  if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta") ||
        TString(variable).EqualTo("phipt") || TString(variable).EqualTo("phieta"))) {
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

    weightsFile->GetObject(
      "ccdb_object", baseList); // TBI 20231008 for simplicity, harwired name
                                // of base TList is "ccdb_object" also for
                                // AliEn case, see if I need to change this
    if (!baseList) {
      // weightsFile->ls();
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
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

    baseList =
      reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath)
                                                  .ReplaceAll("/alice-ccdb.cern.ch/", "")
                                                  .Data()));

    if (!baseList) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    listWithRuns =
      reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
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

    weightsFile->GetObject(
      "ccdb_object", baseList); // TBI 20231008 for simplicity, harwired name
                                // of base TList is "ccdb_object" also for
                                // local case, see if I need to change this
    if (!baseList) {
      // weightsFile->ls();
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

    listWithRuns =
      reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
           __PRETTY_FUNCTION__, __LINE__);
    }

  } // else {

  // g) The final touch on histogram with weights:
  if (-1 == bin) {
    // Integrated weights:
    if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }
    hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, Form("%s_%s", variable, tc.fTaskName.Data())));
    if (!hist) {
      hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, Form("%s", variable))); // yes, for some simple tests I can have only histogram named e.g. 'phi'
    }
    if (!hist) {
      LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
      listWithRuns->ls();
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }
    hist->SetDirectory(0);
    hist->SetTitle(Form("%s, %s", filePath, runNumber));

  } else {
    // Differential weights:
    if (!(TString(variable).EqualTo("phipt") || TString(variable).EqualTo("phieta"))) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }
    hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, Form("%s[%d]_%s", variable, bin, tc.fTaskName.Data())));
    if (!hist) {
      hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, Form("%s[%d]", variable, bin))); // yes, for some simple tests I can have only histogram named e.g. 'phipt[0]'
    }
    if (!hist) {
      LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
      listWithRuns->ls();
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
    }
    hist->SetDirectory(0);

    if (TString(variable).EqualTo("phipt")) {
      hist->SetTitle(Form("%s, %.2f < p_{T} < %.2f", filePath, res.fResultsProVariableLengthBins[AFO_PT]->At(bin), res.fResultsProVariableLengthBins[AFO_PT]->At(bin + 1)));
    }
    if (TString(variable).EqualTo("phieta")) {
      hist->SetTitle(Form("%s, %.2f < #eta < %.2f", filePath, res.fResultsProVariableLengthBins[AFO_ETA]->At(bin), res.fResultsProVariableLengthBins[AFO_ETA]->At(bin + 1)));
    }

  } // else

  return hist;

} // TH1D* GetHistogramWithWeights(const char* filePath, const char* runNumber, const char* variable, Int_t bin = -1)

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
  TObjArray* oa = GetObjArrayWithLabels(t0.fFileWithLabels.Data());
  if (!oa) {
    LOGF(info, "\033[1;33m fFileWithLabels = %s \033[0m",
         t0.fFileWithLabels.Data());
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  // c) Book the placeholder fTest0LabelsPlaceholder for all labels:
  Int_t nLabels = oa->GetEntries();
  t0.fTest0LabelsPlaceholder =
    new TH1I("fTest0LabelsPlaceholder",
             Form("placeholder for all labels, %d in total", nLabels),
             nLabels, 0, nLabels);
  t0.fTest0LabelsPlaceholder->SetStats(kFALSE);

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
    t0.fTest0Labels[order - 1][counter[order - 1]] =
      new TString(oa->At(e)->GetName()); // okay...
    t0.fTest0LabelsPlaceholder->GetXaxis()->SetBinLabel(
      bin++, t0.fTest0Labels[order - 1][counter[order - 1]]->Data());
    // cout<<__LINE__<<":
    // "<<t0.fTest0Labels[order-1][counter[order-1]]->Data()<<endl;
    counter[order - 1]++;
    // cout<<TString(line).Data()<<endl;
    // cout<<oa->GetEntries()<<endl;
  } // for(Int_t e=0; e<nLabels; e++)

  // e) Insantity check on labels:
  //    Here I am merely checking that harmonic larget than gMaxHarmonic was not requested.
  for (Int_t b = 1; b <= t0.fTest0LabelsPlaceholder->GetXaxis()->GetNbins(); b++) {
    TObjArray* temp = TString(t0.fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b)).Tokenize(" ");
    for (Int_t h = 0; h < temp->GetEntries(); h++) {
      if (TMath::Abs(TString(temp->At(h)->GetName()).Atoi()) > gMaxHarmonic) {
        LOGF(info, "\033[1;31m bin = %d, label = %s, gMaxHarmonic = %d\033[0m", b, t0.fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b), (Int_t)gMaxHarmonic);
        LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
      }          // if(TString(temp->At(h)->GetName()).Atoi() > gMaxHarmonic) {
    }            // for(Int_t h = 0; h < temp->GetEntries(); h++) {
    delete temp; // yes, otherwise it's a memory leak
  }              // for(Int_t b = 1; b <= t0.fTest0LabelsPlaceholder->GetXaxis()->GetNbins(); b++) {

} // void StoreLabelsInPlaceholder()

//============================================================

Bool_t RetrieveCorrelationsLabels()
{
  // Generate the labels of all correlations of interest, i.e. retrieve them
  // from TH1I *t0.fTest0LabelsPlaceholder

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  Int_t counter[gMaxCorrelator] = {0}; // is this safe?
  for (Int_t o = 0; o < gMaxCorrelator; o++) {
    counter[o] = 0;
  } // now it's safe :-)

  Int_t nBins = t0.fTest0LabelsPlaceholder->GetXaxis()->GetNbins();

  Int_t order = -44;
  for (Int_t b = 1; b <= nBins; b++) {
    TObjArray* oa = TString(t0.fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b))
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
    t0.fTest0Labels[order - 1][counter[order - 1]] = new TString(
      t0.fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b)); // okay...
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

  if (!pw.fWeightsHist[ppe]) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m",
         __PRETTY_FUNCTION__, __LINE__);
  }

  Int_t bin = pw.fWeightsHist[ppe]->FindBin(value);
  Double_t weight = 0.;
  if (bin > pw.fWeightsHist[ppe]->GetNbinsX()) {
    weight = 0.; // we are in the overflow, ignore this particle TBI_20210524 is
                 // this really the correct procedure?
  } else {
    weight = pw.fWeightsHist[ppe]->GetBinContent(bin);
  }

  return weight;

} // Weight(const Double_t &value, const char *variable) // value, [phi,pt,eta]

//============================================================

Double_t DiffWeight(const Double_t& valueY, const Double_t& valueX, const char* variableX)
{
  // Determine differential particle weight y(x). For the time being, "y = phi" always, but this can be generalized.

  // *) Determine first to which bin the 'valueX' corresponds to.
  //    Based on that, I decide from which histogram I fetch weight for y. See MakeWeights.C

  // TBI 20240208 I need to add support below also for fixed binning case, using fResultsProFixedLengthBins, not only for fResultsProVariableLengthBins
  // TBI 20231026 I do it at the moment this way just to move on, but this can be optimized clearly.

  // *) Mapping between enum's "variableX" on one side, and enum "eAsFunctionOf" on the other:
  Int_t AFO_var = -1;                     // this local variable determines the enum "eAsFunctionOf" which corresponds to enum "eqvectorKine"
  Int_t AFO_weight = -1;                  // this local variable determines the enum "eDiffWeights" which corresponds to enum "eqvectorKine"
  if (TString(variableX).EqualTo("pt")) { // TBI 20240215 check if I can optimize here, i.e. by using strcmp
    AFO_var = AFO_PT;
    AFO_weight = wPHIPT;
  } else if (TString(variableX).EqualTo("eta")) {
    AFO_var = AFO_ETA;
    AFO_weight = wPHIETA;
  }

  // *) Okay, let's do it:
  Int_t binX = 1;
  // if(fInsanityChecksForEachParticle) {
  if (true) { // TBI 20240208 add support to switch on and off this check, which is computationally heavy

    if (valueX < res.fResultsProVariableLengthBins[AFO_var]->At(0)) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
      // underflow. this means that I didn't use the same cuts now, and I was using when making the particle weights. Adjust the cuts.
    }
    if (valueX >= res.fResultsProVariableLengthBins[AFO_var]->At(res.fResultsProVariableLengthBins[AFO_var]->GetSize() - 1)) {
      LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
      // overflow. this means that I didn't use the same cuts now, and I was using when making the particle weights. Adjust the cuts.
    }
  } // fInsanityChecksForEachParticle

  // *) TBI add some comment...
  for (Int_t e = 1; e < res.fResultsProVariableLengthBins[AFO_var]->GetSize(); e++) {
    // Since I set binX = 1, intentionally I skip the first element in the loop, and start from e = 1, instead of e = 0.
    if (valueX < res.fResultsProVariableLengthBins[AFO_var]->At(e)) {
      binX = e;
      break;
    } // gotcha
  }

  // *) Finally, determine weight for y(x):
  if (!pw.fDiffWeightsHist[AFO_weight][binX - 1]) {
    LOGF(info, "\033[1;32mvalueY = %f\033[0m", valueY);
    LOGF(info, "\033[1;32mvalueX = %f\033[0m", valueX);
    LOGF(info, "\033[1;32mvariableX = %s\033[0m", variableX);
    LOGF(info, "\033[1;32mAFO_weight = %d\033[0m", AFO_weight);
    LOGF(info, "\033[1;32mbinX = %d\033[0m", binX);
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  Int_t bin = pw.fDiffWeightsHist[AFO_weight][binX - 1]->FindBin(valueY); // binX - 1, because I histogram for first bin in X is labeled with "[0]", etc.
  Double_t weight = 0.;
  if (bin > pw.fDiffWeightsHist[AFO_weight][binX - 1]->GetNbinsX()) {
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __PRETTY_FUNCTION__, __LINE__); // TBI 20240208 re-think what to do here
  } else {
    weight = pw.fDiffWeightsHist[AFO_weight][binX - 1]->GetBinContent(bin);
  }

  return weight;

} // DiffWeight(const Double_t &valueY, const Double_t &valueX, const char* variableX)

//============================================================

void GetParticleWeights()
{
  // Get the particle weights. Call this function only once.

  //    TBI 20231012 Here the current working assumption is that:
  //    1) Corrections do not change within a given run;
  //    2) Hyperloop proceeses the dataset one masterjob per run number.
  //    If any of these 2 assumptions are violated, this code will have to be modified.

  // a) Integrated weights;
  // b) Differential weights.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // a) Integrated weights:
  // integrated phi weights:
  if (pw.fUseWeights[wPHI]) {
    TH1D* phiWeights = GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "phi");
    if (!phiWeights) {
      LOGF(fatal, "in function \033[1;31m%s at line %d, phiWeights is NULL. Check the external file %s with particle weights\033[0m", __PRETTY_FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
    }
    SetWeightsHist(phiWeights, "phi");
  }

  // integrated pt weights:
  if (pw.fUseWeights[wPT]) {
    TH1D* ptWeights = GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "pt");
    if (!ptWeights) {
      LOGF(fatal, "in function \033[1;31m%s at line %d, ptWeights is NULL. Check the external file %s with particle weights\033[0m", __PRETTY_FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
    }
    SetWeightsHist(ptWeights, "pt");
  }

  // integrated eta weights:
  if (pw.fUseWeights[wETA]) {
    TH1D* etaWeights = GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "eta");
    if (!etaWeights) {
      LOGF(fatal, "in function \033[1;31m%s at line %d, etaWeights is NULL. Check the external file %s with particle weights\033[0m", __PRETTY_FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
    }
    SetWeightsHist(etaWeights, "eta");
  }

  // b) Differential weights:
  // differential phi(pt) weights:
  if (pw.fUseDiffWeights[wPHIPT]) {
    TH1D* phiptWeights = NULL;
    Int_t nPtBins = res.fResultsPro[AFO_PT]->GetXaxis()->GetNbins();
    for (Int_t b = 0; b < nPtBins; b++) {
      phiptWeights = GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "phipt", b);
      if (!phiptWeights) {
        LOGF(fatal, "in function \033[1;31m%s at line %d, phiptWeights is NULL. Check the external file %s with particle weights\033[0m", __PRETTY_FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
      }
      SetDiffWeightsHist(phiptWeights, "phipt", b);
    }
  } // if (pw.fUseDiffWeights[wPHIPT]) {

  // differential phi(eta) weights:
  if (pw.fUseDiffWeights[wPHIETA]) {
    TH1D* phietaWeights = NULL;
    Int_t nEtaBins = res.fResultsPro[AFO_ETA]->GetXaxis()->GetNbins();
    for (Int_t b = 0; b < nEtaBins; b++) {
      GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "phieta", b);
      if (!phietaWeights) {
        LOGF(fatal, "in function \033[1;31m%s at line %d, phietaWeights is NULL. Check the external file %s with particle weights\033[0m", __PRETTY_FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
      }
      SetDiffWeightsHist(phietaWeights, "phieta", b);
    } // for(Int_t b=0; b<nEtaBins; b++) {
  }   // if (pw.fUseDiffWeights[wPHIETA]) {

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
  Int_t rs = -44;                            // reconstructed or simulated
  if (tc.fGenericRec || tc.fGenericRecSim) { // yes, for tc.fGenericRecSim I take info from Rec part
    rs = eRec;
  } else if (tc.fGenericSim) {
    rs = eSim;
  } else {
    LOGF(fatal, "in function \033[1;31m%s at line %d, not a single flag gProcess* is true \033[0m", __PRETTY_FUNCTION__, __LINE__);
  }

  // *) Okay, do the thing:
  if (eh.fEventHistograms[eNumberOfEvents][rs][eAfter] && eh.fEventHistograms[eNumberOfEvents][rs][eAfter]->GetBinContent(1) == eh.fEventCuts[eNumberOfEvents][eMax]) {
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

  Int_t nParticles = ebye.fSelectedTracks;
  /* TBI 20231108 enable eventually
  if(fUseFixedNumberOfRandomlySelectedParticles)
  {
   nParticles = 0;
   for(Int_t i=0;i<nl.ftaNestedLoops[0]->GetSize();i++)
   {
    if(TMath::Abs(nl.ftaNestedLoops[0]->GetAt(i)) > 0. && TMath::Abs(nl.ftaNestedLoops[1]->GetAt(i)) > 0.){nParticles++;}
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
    Double_t dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    Double_t dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
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
      Double_t dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      Double_t dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
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
        Double_t dPhi3 = nl.ftaNestedLoops[0]->GetAt(i3);
        Double_t dW3 = nl.ftaNestedLoops[1]->GetAt(i3);
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
          Double_t dPhi4 = nl.ftaNestedLoops[0]->GetAt(i4);
          Double_t dW4 = nl.ftaNestedLoops[1]->GetAt(i4);
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
            Double_t dPhi5 = nl.ftaNestedLoops[0]->GetAt(i5);
            Double_t dW5 = nl.ftaNestedLoops[1]->GetAt(i5);
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
              Double_t dPhi6 = nl.ftaNestedLoops[0]->GetAt(i6);
              Double_t dW6 = nl.ftaNestedLoops[1]->GetAt(i6);
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
                Double_t dPhi7 = nl.ftaNestedLoops[0]->GetAt(i7);
                Double_t dW7 = nl.ftaNestedLoops[1]->GetAt(i7);
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
                  Double_t dPhi8 = nl.ftaNestedLoops[0]->GetAt(i8);
                  Double_t dW8 = nl.ftaNestedLoops[1]->GetAt(i8);
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
                    Double_t dPhi9 = nl.ftaNestedLoops[0]->GetAt(i9);
                    Double_t dW9 = nl.ftaNestedLoops[1]->GetAt(i9);
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
                      Double_t dPhi10 = nl.ftaNestedLoops[0]->GetAt(i10);
                      Double_t dW10 = nl.ftaNestedLoops[1]->GetAt(i10);
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
                        Double_t dPhi11 = nl.ftaNestedLoops[0]->GetAt(i11);
                        Double_t dW11 = nl.ftaNestedLoops[1]->GetAt(i11);
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
                          Double_t dPhi12 = nl.ftaNestedLoops[0]->GetAt(i12);
                          Double_t dW12 = nl.ftaNestedLoops[1]->GetAt(i12);
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
  // b) For simulated data, determine centrality directly from impact parameter;
  // c) Same as a), just for converted Run 2 data;
  // d) Same as b), just for converted Run 2 data;
  // e) Same as a), just for converted Run 1 data;
  // f) Same as b), just for converted Run 1 data;
  // g) Test case.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // a) For real data, determine centrality from default centrality estimator:
  if constexpr (rs == eRec || rs == eRecAndSim) {
    if (ec.fCentralityEstimator.EqualTo("centFT0M", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centFT0M();
    } else if (ec.fCentralityEstimator.EqualTo("CentFV0A", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centFV0A();
    } else if (ec.fCentralityEstimator.EqualTo("CentNTPV", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centNTPV();
    } else {
      LOGF(fatal, "in function \033[1;31m%s at line %d. centrality estimator = %d is not supported yet. \033[0m", __PRETTY_FUNCTION__, __LINE__, ec.fCentralityEstimator.Data());
    }
    // TBI 20240120 I could also here access also corresponding simulated centrality from impact parameter, if available through collision.has_mcCollision()
  }

  // b) For simulated data, determine centrality directly from impact parameter:
  if constexpr (rs == eSim) {
    ebye.fCentrality = -44.; // TBI 20240120 add support eventualy
  }

  // c) Same as a), just for converted Run 2 data:
  if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2) {
    if (ec.fCentralityEstimator.EqualTo("centRun2V0M", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centRun2V0M();
    } else if (ec.fCentralityEstimator.EqualTo("CentRun2SPDTracklets", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centRun2SPDTracklets();
    } else {
      LOGF(fatal, "in function \033[1;31m%s at line %d. centrality estimator = %d is not supported yet. \033[0m", __PRETTY_FUNCTION__, __LINE__, ec.fCentralityEstimator.Data());
    }
    // TBI 20240120 I could also here access also corresponding simulated centrality from impact parameter, if available through collision.has_mcCollision()
  }

  // d) Same as b), just for converted Run 2 data:
  if constexpr (rs == eSim_Run2) {
    ebye.fCentrality = -44.; // TBI 20240120 add support eventualy
  }

  // e) Same as a), just for converted Run 1 data:
  if constexpr (rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    if (ec.fCentralityEstimator.EqualTo("centRun2V0M", TString::kIgnoreCase)) {
      // ebye.fCentrality = collision.centRun2V0M(); // TBI 20240224 enable when I add support for RecAndSim_Run1
    } else if (ec.fCentralityEstimator.EqualTo("CentRun2SPDTracklets", TString::kIgnoreCase)) {
      // ebye.fCentrality = collision.centRun2SPDTracklets(); // TBI 20240224 enable when I add support for RecAndSim_Run1
    } else {
      LOGF(fatal, "in function \033[1;31m%s at line %d. centrality estimator = %d is not supported yet. \033[0m", __PRETTY_FUNCTION__, __LINE__, ec.fCentralityEstimator.Data());
    }
    // TBI 20240120 I could also here access also corresponding simulated centrality from impact parameter, if available through collision.has_mcCollision()
  }

  // f) Same as b), just for converted Run 1 data:
  if constexpr (rs == eSim_Run1) {
    ebye.fCentrality = -44.; // TBI 20240120 add support eventualy
  }

  // g) Test case:
  if constexpr (rs == eTest) {
    ebye.fCentrality = gRandom->Uniform(0., 100.);
  }

  // *) Print centrality for the audience...:
  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m ebye.fCentrality = %f\033[0m", ebye.fCentrality);
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

void BailOut()
{
  // Use only locally - bail out if maximum number of events was reached, and dump all results by that point in a local ROOT file.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // *) Local variables: TBI 20240130 shall I promote 'em to data members + add support for configurables?
  TString sBailOutFile = "AnalysisResultsBailOut.root";
  TString sDirectoryFile = "multiparticle-correlations-a-b";

  // *) Info message:
  if (eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]) {
    LOGF(info, "\033[1;32m=> Per request, bailing out after %d selected events in the local file %s .\n\033[0m", (Int_t)eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]->GetBinContent(1), sBailOutFile.Data());
  }

  // *) Okay, let's bail out intentionally:
  TFile* f = new TFile(sBailOutFile.Data(), "recreate");
  TDirectoryFile* dirFile = new TDirectoryFile(sDirectoryFile.Data(), sDirectoryFile.Data());
  // TBI 20240130 I cannot add here fBaseList directtly, since that one is declared as OutputObj<TList>
  // Therefore, adding one-by-one nested TList's I want to bail out.
  // Keep in sync with BookAndNestAllLists().
  TList* bailOutList = new TList(); // this is sort of 'fake' fBaseList
  bailOutList->SetOwner(kTRUE);
  bailOutList->SetName(sBaseListName.Data());
  bailOutList->Add(fBasePro); // yes, this one needs a special treatment
  bailOutList->Add(qa.fQAList);
  bailOutList->Add(eh.fEventHistogramsList);
  bailOutList->Add(ec.fEventCutsList);
  bailOutList->Add(ph.fParticleHistogramsList);
  bailOutList->Add(pc.fParticleCutsList);
  bailOutList->Add(qv.fQvectorList);
  bailOutList->Add(mupa.fCorrelationsList);
  bailOutList->Add(pw.fWeightsList);
  bailOutList->Add(nl.fNestedLoopsList);
  bailOutList->Add(t0.fTest0List);
  bailOutList->Add(res.fResultsList);

  // *) Add list with nested list to TDirectoryFile:
  dirFile->Add(bailOutList, kTRUE);
  dirFile->Write(dirFile->GetName(), TObject::kSingleKey + TObject::kOverwrite);
  delete dirFile;
  dirFile = NULL;
  f->Close();

  // *) Hasta la vista:
  LOGF(fatal, "\n\nHasta la vista - bailed out intentionally in function \033[1;31m%s at line %d\n The output file is: %s\n\n\033[0m", __PRETTY_FUNCTION__, __LINE__, sBailOutFile.Data());

} // void BailOut()

//============================================================

void Fillqvector(const Double_t& dPhi, const Double_t& kineVarValue, eqvectorKine kineVarChoice)
{
  // Fill differential q-vector, in generic kinematic variable. Here "kine" originally meant vs. pt or vs. eta, now it's general.
  // Example usage: this->Fillqvector(dPhi, dPt, PTq);

  if (tc.fVerboseForEachParticle) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // *) Mapping between enum's "eqvectorKine" on one side, and "eAsFunctionOf", "eWeights" and "eDiffWeights" on the other:
  //    TBI 20240212 I could promote this also to a member function, if I need it elsewhere. Or I could use TExMap?
  Int_t AFO_var = -1;        // this local variable determines the enum "eAsFunctionOf" which corresponds to enum "eqvectorKine"
  Int_t AFO_weight = -1;     // this local variable determined the enum "eWeights" which corresponds to enum "eqvectorKine"
  Int_t AFO_diffWeight = -1; // this local variable determines the enum "eDiffWeights" which corresponds to enum "eqvectorKine"
  TString AFO_name = "";     // TBI 20240212 most likely, I won't need this one in the final version
  switch (kineVarChoice) {
    case PTq:
      AFO_var = AFO_PT;
      AFO_weight = wPT;
      AFO_diffWeight = wPHIPT;
      AFO_name = "pt";
      break;
    case ETAq:
      AFO_var = AFO_ETA;
      AFO_weight = wETA;
      AFO_diffWeight = wPHIETA;
      AFO_name = "eta";
      break;
    default:
      LOGF(fatal, "in function \033[1;31m%s at line %d. This kineVarChoice = %d is not supported yet. \033[0m", __PRETTY_FUNCTION__, __LINE__, (Int_t)kineVarChoice);
      break;
  } // switch(kineVarChoice)

  // *) Get the desired bin number:
  Int_t bin = -1;
  if (res.fResultsPro[AFO_var]) {
    bin = res.fResultsPro[AFO_var]->FindBin(kineVarValue);
    if (0 >= bin || res.fResultsPro[AFO_var]->GetNbinsX() < bin) { // either underflow or overflow is hit, meaning that histogram is booked in narrower range than cuts
      LOGF(fatal, "in function \033[1;31m%s at line %d => kineVarChoice = %d, bin = %d, kineVarValue = %f \033[0m", __PRETTY_FUNCTION__, __LINE__, (Int_t)kineVarChoice, bin, kineVarValue);
    }
  }

  // *) Get all integrated kinematic weights:
  Double_t wToPowerP = 1.;     // weight raised to power p
  Double_t kineVarWeight = 1.; // e.g. this can be integrated pT or eta weight
  if (pw.fUseWeights[AFO_weight]) {
    kineVarWeight = Weight(kineVarValue, AFO_name.Data()); // corresponding e.g. pt or eta weight
    if (!(kineVarWeight > 0.)) {
      LOGF(fatal, "in function \033[1;31m%s at line %d. kineVarWeight is not positive \033[0m", __PRETTY_FUNCTION__, __LINE__);
      // TBI 20240212 or could I just skip this particle?
    }
  } // if(fUseWeights[AFO_weight]) {

  // *) Get all differential phi-weights for this kinematic variable:
  //    Remark: special treatment is justified for phi-weights, because q-vector is defined in terms of phi-weights.
  Double_t diffPhiWeightsForThisKineVar = 1.;
  if (pw.fUseDiffWeights[AFO_diffWeight]) {
    diffPhiWeightsForThisKineVar = DiffWeight(dPhi, kineVarValue, AFO_name.Data()); // corresponding differential phi weight as a function of e.g. pt or eta
    if (!(diffPhiWeightsForThisKineVar > 0.)) {
      LOGF(fatal, "in function \033[1;31m%s at line %d. diffPhiWeightsForThisKineVar is not positive \033[0m", __PRETTY_FUNCTION__, __LINE__);
      // TBI 20240212 or could I just skip this particle?
    }
  } // if(pw.fUseDiffWeights[AFO_diffWeight]) {

  // *) Finally, fill differential q-vector in that bin:
  for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
    for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
      if (pw.fUseWeights[AFO_weight] || pw.fUseDiffWeights[AFO_diffWeight]) {
        wToPowerP = pow(diffPhiWeightsForThisKineVar * kineVarWeight, wp); // TBI 20240212 supported at the moment: e.g. q-vector vs pt can be weighted only with diff. phi(pt) and integrated pt weights. It cannot be weighted in addition with eta weights, since in any case I anticipate I will do always 1-D analysis, by integrating out all other dependencies
      }
      qv.fqvector[kineVarChoice][bin - 1][h][wp] += TComplex(wToPowerP * TMath::Cos(h * dPhi), wToPowerP * TMath::Sin(h * dPhi));
    } // for(Int_t wp=0;wp<gMaxCorrelator+1;wp++)
  }   // for(Int_t h=0;h<gMaxHarmonic*gMaxCorrelator+1;h++)

  // *) TBI 20240208 add here support for differential nested loops
  /*
  if(nl.fCalculateCustomNestedLoop)
  {
    ftaNestedLoopsKine[PTq][bin-1][0]->AddAt(dPhi,fqVectorEntries[PTq][bin-1]);
    ftaNestedLoopsKine[PTq][bin-1][1]->AddAt(wPhi*wPt*wEta,fqVectorEntries[PTq][bin-1]);
  }
  */

  // *) Multiplicity counter in this bin:
  qv.fqVectorEntries[kineVarChoice][bin - 1]++; // count number of particles in this pt bin in this event

} // void Fillqvector(const Double_t &dPhi, const Double_t &kineVarValue, eqvectorKine kineVarChoice)

//============================================================

void CalculateEverything()
{
  // Calculate everything for selected events and particles.
  // Remark: Data members for Q-vectors, containers for nested loops, etc., must all be filled when this function is called.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__);
  }

  // *) Progress info:
  if (eh.fEventHistograms[eNumberOfEvents][eRec][eBefore] && eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]) {
    LOGF(info, "\033[1;32m=> Processing event %d/%d (selected/total) .... \033[0m", (Int_t)eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]->GetBinContent(1), (Int_t)eh.fEventHistograms[eNumberOfEvents][eRec][eBefore]->GetBinContent(1));
  }

  // *) Calculate multiparticle correlations (standard, isotropic, same harmonic):
  if (mupa.fCalculateCorrelations) {
    this->CalculateCorrelations();
  }

  // *) Calculate Test0: TBI 20240110 name convention
  //    Remark: integrated, vs. M and vs. centrality are all calculated here
  if (t0.fCalculateTest0) {
    this->CalculateTest0();
  }

  // *) Calculate kine Test0: TBI 20240110 name convention
  //    Remark: vs. pt, vs. eta, etc., are all calculated here
  if (t0.fCalculateTest0AsFunctionOf[AFO_PT]) {
    this->CalculateKineTest0("pt");
  }
  if (t0.fCalculateTest0AsFunctionOf[AFO_ETA]) {
    this->CalculateKineTest0("eta");
  }

  // *) Calculate nested loops:
  if (nl.fCalculateNestedLoops) {
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

  if (tc.fDryRun) {
    LOGF(info, "\033[1;32m%s => This is a dry run, bailing out immediately from Steer(...)\033[0m", __FUNCTION__);
    return;
  }

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s starting ...\033[0m", __FUNCTION__); // just a bare function name
  }

  // *) Global timestamp:
  if (tc.fUseStopwatch) {
    LOGF(info, "\033[1;32m\n\n=> Global timer: Steer begins ... %.6f\n\n\033[0m", tc.fTimer[eGlobal]->RealTime());
    tc.fTimer[eGlobal]->Continue(); // yes
  }

  // *) Do all thingies before starting to process data from this collision (e.g. count number of events, fetch the run number, etc.):
  Preprocess(collision);

  // *) Determine collision centrality:
  DetermineCentrality<rs>(collision);

  // *) Fill event histograms before event cuts:
  FillEventHistograms<rs>(collision, tracks, eBefore);

  // *) Print info on the current event number (total, before cuts):
  if (tc.fVerbose) {
    if (eh.fEventHistograms[eNumberOfEvents][eRec][eBefore]) {
      LOGF(info, "\033[1;32m%s : processing event %d\033[0m", __FUNCTION__, (Int_t)eh.fEventHistograms[eNumberOfEvents][eRec][eBefore]->GetBinContent(1));
    } else if (eh.fEventHistograms[eNumberOfEvents][eSim][eBefore]) {
      LOGF(info, "\033[1;32m%s : processing event %d\033[0m", __FUNCTION__, (Int_t)eh.fEventHistograms[eNumberOfEvents][eSim][eBefore]->GetBinContent(1));
    }
  }

  // *) Event cuts:
  if (!EventCuts<rs>(collision, tracks)) {
    return;
  }

  // *) Main loop over particles:
  MainLoopOverParticles<rs>(tracks);

  // *) Remaining event cuts which can be applied only after the loop over particles is performed:
  if ((ebye.fSelectedTracks < eh.fEventCuts[eSelectedTracks][eMin]) || (ebye.fSelectedTracks > eh.fEventCuts[eSelectedTracks][eMax])) {
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

  // *) Print info on the current event number after cuts:
  if (tc.fVerbose) {
    if (eh.fEventHistograms[eNumberOfEvents][eRec][eBefore]) {
      LOGF(info, "\033[1;32m%s : this event passed all cuts %d/%d\033[0m", __FUNCTION__, (Int_t)eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]->GetBinContent(1), (Int_t)eh.fEventHistograms[eNumberOfEvents][eRec][eBefore]->GetBinContent(1));
    } else if (eh.fEventHistograms[eNumberOfEvents][eSim][eBefore]) {
      LOGF(info, "\033[1;32m%s : this event passed all cuts %d/%d\033[0m", __FUNCTION__, (Int_t)eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]->GetBinContent(1), (Int_t)eh.fEventHistograms[eNumberOfEvents][eSim][eBefore]->GetBinContent(1));
    }
  }

  // *) Global timestamp:
  if (tc.fUseStopwatch) {
    LOGF(info, "\033[1;32m\n\n=> Global timer: Steer ends ... %.6f\n\n\033[0m", tc.fTimer[eGlobal]->RealTime());
    tc.fTimer[eGlobal]->Continue(); // yes
  }

} // template <eRecSim rs, typename T1, typename T2> void Steer(T1 const* collision, T2 const* tracks)

//============================================================

template <eRecSim rs, typename T>
void MainLoopOverParticles(T const& tracks)
{
  // This is the main loop over particles, in which Q-vectors and particle histograms are filled, particle cuts applied, etc.

  // Remark #1:
  // *) To process only reconstructed Run 3, use processRec(...), i.e. set field "processRec": "true" in json config
  // *) To process both reconstructed and simulated Run 3, use processRecSim(...), i.e. set field "processRecSim": "true" in json config
  // *) To process only simulated Run 3, use processSim(...), i.e. set field "processSim": "true" in json config

  // Remark #2:
  // *) To process Run 1 and Run 2 converted data, use in the same spirit e.g. processRec_Run2(...), i.e. set field "processRec_Run2": "true" in json config

  // Remark #3:
  // *) There is also processTest(...), to process data with minimum subscription to the tables. To use it, set field "processTest": "true" in json config

  if (tc.fVerbose) {
    // LOGF(info, "\033[1;32m%s\033[0m", __PRETTY_FUNCTION__); // full function signature (including arguments, etc.), too verbose here...
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__); // just a bare function name
  }

  // *) Local kinematic variables and corresponding particle weights:
  Double_t dPhi = 0., wPhi = 1.; // azimuthal angle and corresponding phi weight
  Double_t dPt = 0., wPt = 1.;   // transverse momentum and corresponding pt weight
  Double_t dEta = 0., wEta = 1.; // pseudorapidity and corresponding eta weight
  Double_t wToPowerP = 1.;       // weight raised to power p
  ebye.fSelectedTracks = 0;      // reset number of selected tracks

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

  // *) Local timestamp:
  if (tc.fUseStopwatch) {
    LOGF(info, "\033[1;32m\n\n=> Local timer starts at line %d\n\n\033[0m", __LINE__);
    tc.fTimer[eLocal]->Reset();
    tc.fTimer[eLocal]->Start();
  }

  // *) Main loop over particles:
  // for (auto& track : tracks) { // default standard way of looping of tracks
  auto track = tracks.iteratorAt(0); // set the type and scope from one instance
  for (int64_t i = 0; i < tracks.size(); i++) {

    // *) Skip track objects which are not valid tracks (e.g. Run 2 and 1 tracklets, etc.):
    if (!ValidTrack<rs>(track)) {
      continue;
    }

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
    //  Kinematics (Remark: for "eRecSim" processing, kinematics is taken from "reconstructed"):
    dPhi = track.phi();
    dPt = track.pt();
    dEta = track.eta();

    // Particle weights:
    if (pw.fUseWeights[wPHI]) {
      wPhi = Weight(dPhi, "phi"); // corresponding phi weight
      if (!(wPhi > 0.)) {
        LOGF(error, "\033[1;33m%s wPhi is not positive, skipping this particle for the time being...\033[0m", __PRETTY_FUNCTION__);
        LOGF(error, "dPhi = %f\nwPhi = %f", dPhi, wPhi);
        continue;
      }
    } // if(pw.fUseWeights[wPHI])
    if (pw.fUseWeights[wPT]) {
      wPt = Weight(dPt, "pt"); // corresponding pt weight
      if (!(wPt > 0.)) {
        LOGF(error, "\033[1;33m%s wPt is not positive, skipping this particle for the time being...\033[0m", __PRETTY_FUNCTION__);
        LOGF(error, "dPt = %f\nwPt = %f", dPt, wPt);
        continue;
      }
    } // if(pw.fUseWeights[wPT])
    if (pw.fUseWeights[wETA]) {
      wEta = Weight(dEta, "eta"); // corresponding eta weight
      if (!(wEta > 0.)) {
        LOGF(error, "\033[1;33m%s wEta is not positive, skipping this particle for the time being...\033[0m", __PRETTY_FUNCTION__);
        LOGF(error, "dEta = %f\nwEta = %f", dEta, wEta);
        continue;
      }
    } // if(pw.fUseWeights[wETA])

    if (qv.fCalculateQvectors) {
      for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
        for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
          if (pw.fUseWeights[wPHI] || pw.fUseWeights[wPT] || pw.fUseWeights[wETA]) {
            wToPowerP = pow(wPhi * wPt * wEta, wp);
          }
          qv.fQvector[h][wp] += TComplex(wToPowerP * TMath::Cos(h * dPhi), wToPowerP * TMath::Sin(h * dPhi));
        } // for(Int_t wp=0;wp<gMaxCorrelator+1;wp++)
      }   // for(Int_t h=0;h<gMaxHarmonic*gMaxCorrelator+1;h++)
    }     // if (qv.fCalculateQvectors) {

    // *) Nested loops containers:
    if (nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoop) {
      if (nl.ftaNestedLoops[0]) {
        nl.ftaNestedLoops[0]->AddAt(dPhi, ebye.fSelectedTracks);
      } // remember that the 2nd argument here must start from 0
      if (nl.ftaNestedLoops[1]) {
        nl.ftaNestedLoops[1]->AddAt(wPhi * wPt * wEta, ebye.fSelectedTracks);
      } // remember that the 2nd argument here must start from 0
    }   // if(nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoop)

    // *) Differential q-vectors:
    if (qv.fCalculateQvectors && t0.fCalculateTest0AsFunctionOf[AFO_PT]) {
      this->Fillqvector(dPhi, dPt, PTq); // first 2 arguments are passed by reference, 3rd argument is enum
    }
    if (qv.fCalculateQvectors && t0.fCalculateTest0AsFunctionOf[AFO_ETA]) {
      this->Fillqvector(dPhi, dEta, ETAq); // first 2 arguments are passed by reference, 3rd argument is enum
    }

    // *) Counter of selected tracks in the current event:
    ebye.fSelectedTracks++;
    if (ebye.fSelectedTracks >= cSelectedTracks_max) {
      break;
    }

    // *) Break the loop if fixed number of particles is taken randomly from each event (use always in combination with tc.fUseFisherYates = kTRUE):
    if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && tc.fFixedNumberOfRandomlySelectedTracks == ebye.fSelectedTracks) {
      LOGF(info, "\033[1;32mBreaking the loop over particles, since requested fixed number of %d particles was reached\033[0m", tc.fFixedNumberOfRandomlySelectedTracks);
      break;
    }

  } // for (auto& track : tracks)

  // *) Local timestamp:
  if (tc.fUseStopwatch) {
    LOGF(info, "\033[1;32m\n\n=> Local timer ends at line %d, time elapsed ... %.6f\n\n\033[0m", __LINE__, tc.fTimer[eLocal]->RealTime());
    tc.fTimer[eLocal]->Continue();
  }

  // *) Insanity check on fixed number of randomly selected tracks:
  if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && tc.fFixedNumberOfRandomlySelectedTracks < ebye.fSelectedTracks) {
    LOGF(fatal, "\033[1;31mIn this event there are too few particles (ebye.fSelectedTracks = %d), and requested number of fixed number randomly selected tracks %d couldn't be reached\033[0m", ebye.fSelectedTracks, tc.fFixedNumberOfRandomlySelectedTracks);
  }

} // template <eRecSim rs, typename T> void MainLoopOverParticles(T const& tracks) {

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_MEMBERFUNCTIONS_H_
