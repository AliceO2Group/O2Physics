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
#include <string>

//============================================================

void BookBaseList()
{
  // Book base TList and store task configuration.

  // a) Book base TList;
  // b) Store task configuration.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book base TList:
  TList* temp = new TList();
  temp->SetOwner(true);
  fBaseList.setObject(temp);

  // b) Store task configuration:
  fBasePro = new TProfile("fBasePro", "flags for the whole analysis", eConfiguration_N - 1, 0.5, static_cast<float>(eConfiguration_N) - 0.5);
  // yes, eConfiguration_N - 1 and -0.5, because eConfiguration kicks off from 1
  fBasePro->SetStats(false);
  fBasePro->SetLineColor(eColor);
  fBasePro->SetFillColor(eFillColor);

  // Remark: If I want to change the ordering of bin labels, simply change the
  // ordering in enum eConfiguration { ... }, nothing needs to be changed here.
  fBasePro->GetXaxis()->SetBinLabel(eTaskIsConfiguredFromJson, Form("fTaskIsConfiguredFromJson = %s", tc.fTaskIsConfiguredFromJson.Data()));

  fBasePro->GetXaxis()->SetBinLabel(eTaskName, Form("fTaskName = %s", tc.fTaskName.Data()));

  fBasePro->GetXaxis()->SetBinLabel(eRunNumber, Form("fRunNumber = %s", "__RUN_NUMBER__"));
  // I have to do it this way via placeholder, because run number is available only when i start to process data.
  // Then, I replace placeholder with run number in PropagateRunNumber(...)

  fBasePro->GetXaxis()->SetBinLabel(eDryRun, "fDryRun");
  fBasePro->Fill(eDryRun, static_cast<double>(tc.fDryRun));

  fBasePro->GetXaxis()->SetBinLabel(eVerbose, "fVerbose");
  fBasePro->Fill(eVerbose, static_cast<double>(tc.fVerbose));

  fBasePro->GetXaxis()->SetBinLabel(eVerboseUtility, "fVerboseUtility");
  fBasePro->Fill(eVerboseUtility, static_cast<double>(tc.fVerboseUtility));

  fBasePro->GetXaxis()->SetBinLabel(eVerboseForEachParticle, "fVerboseForEachParticle");
  fBasePro->Fill(eVerboseForEachParticle, static_cast<double>(tc.fVerboseForEachParticle));

  fBasePro->GetXaxis()->SetBinLabel(eVerboseEventCounter, "fVerboseEventCounter");
  fBasePro->Fill(eVerboseEventCounter, static_cast<double>(tc.fVerboseEventCounter));

  fBasePro->GetXaxis()->SetBinLabel(ePlainPrintout, "fPlainPrintout");
  fBasePro->Fill(ePlainPrintout, static_cast<double>(tc.fPlainPrintout));

  fBasePro->GetXaxis()->SetBinLabel(eDoAdditionalInsanityChecks, "fDoAdditionalInsanityChecks");
  fBasePro->Fill(eDoAdditionalInsanityChecks, static_cast<double>(tc.fDoAdditionalInsanityChecks));

  fBasePro->GetXaxis()->SetBinLabel(eInsanityCheckForEachParticle, "fInsanityCheckForEachParticle");
  fBasePro->Fill(eInsanityCheckForEachParticle, static_cast<double>(tc.fInsanityCheckForEachParticle));

  fBasePro->GetXaxis()->SetBinLabel(eWhichProcess, Form("WhichProcess = %s", tc.fWhichProcess.Data()));

  fBasePro->GetXaxis()->SetBinLabel(eRandomSeed, "fRandomSeed");
  fBasePro->Fill(eRandomSeed, static_cast<double>(tc.fRandomSeed));

  fBasePro->GetXaxis()->SetBinLabel(eUseFisherYates, "fUseFisherYates");
  fBasePro->Fill(eUseFisherYates, static_cast<double>(tc.fUseFisherYates));

  fBasePro->GetXaxis()->SetBinLabel(eFixedNumberOfRandomlySelectedTracks, "fFixedNumberOfRandomlySelectedTracks");
  fBasePro->Fill(eFixedNumberOfRandomlySelectedTracks, static_cast<int>(tc.fFixedNumberOfRandomlySelectedTracks));

  fBasePro->GetXaxis()->SetBinLabel(eUseStopwatch, "fUseStopwatch");
  fBasePro->Fill(eUseStopwatch, static_cast<double>(tc.fUseStopwatch));

  fBasePro->GetXaxis()->SetBinLabel(eFloatingPointPrecision, "fFloatingPointPrecision");
  fBasePro->Fill(eFloatingPointPrecision, tc.fFloatingPointPrecision);

  fBasePro->GetXaxis()->SetBinLabel(eSequentialBailout, "fSequentialBailout");
  fBasePro->Fill(eSequentialBailout, static_cast<double>(tc.fSequentialBailout));

  fBasePro->GetXaxis()->SetBinLabel(eUseSpecificCuts, "fUseSpecificCuts");
  fBasePro->Fill(eUseSpecificCuts, static_cast<double>(tc.fUseSpecificCuts));

  fBasePro->GetXaxis()->SetBinLabel(eWhichSpecificCuts, Form("WhichSpecificCuts = %s", tc.fWhichSpecificCuts.Data()));

  fBaseList->Add(fBasePro);

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookBaseList()

//============================================================

void DefaultConfiguration()
{
  // Default task configuration.
  // a) Default values are hardcoded as Configurables in the file MuPa-Configurables.h

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

  // e) Use whenever you can ConfigurableGroup, keep grouping in sync with all struct's in MuPa-DataMembers.h.
  //    This is needed, to avoid this compilation error:
  //        Framework/StructToTuple.h:286:6: error: only 99 names provided for structured binding
  //    *) As far as I can tell, that means that sum of individual data members + struct fields + individual configurables > 99
  //    *) Therefore, wrap up all data members in some struct fields + use in instead of individual configurables ConfigurableGroup whenever possible.
  //    *) Within a given struct field, number of data members do not add to that number. Also, number of enum fields do not add.

  tc.fTaskIsConfiguredFromJson = TString(cf_tc.cfTaskIsConfiguredFromJson);
  tc.fTaskName = TString(cf_tc.cfTaskName);
  tc.fDryRun = cf_tc.cfDryRun;
  tc.fVerbose = cf_tc.cfVerbose;
  if (tc.fVerbose) {
    StartFunction(__FUNCTION__); // yes, here
  }
  tc.fVerboseUtility = cf_tc.cfVerboseUtility;
  tc.fVerboseForEachParticle = cf_tc.cfVerboseForEachParticle;
  tc.fVerboseEventCounter = cf_tc.cfVerboseEventCounter;
  tc.fVerboseEventCut = cf_tc.cfVerboseEventCut;
  tc.fPlainPrintout = cf_tc.cfPlainPrintout;
  tc.fDoAdditionalInsanityChecks = cf_tc.cfDoAdditionalInsanityChecks;
  // Set automatically what to process, from an implicit variable "doprocessSomeProcessName" within a PROCESS_SWITCH clause:
  tc.fProcess[eProcessRec] = doprocessRec;
  tc.fProcess[eProcessRecSim] = doprocessRecSim;
  tc.fProcess[eProcessSim] = doprocessSim;
  tc.fProcess[eProcessRec_Run2] = doprocessRec_Run2;
  tc.fProcess[eProcessRecSim_Run2] = doprocessRecSim_Run2;
  tc.fProcess[eProcessSim_Run2] = doprocessSim_Run2;
  tc.fProcess[eProcessRec_Run1] = doprocessRec_Run1;
  tc.fProcess[eProcessRecSim_Run1] = doprocessRecSim_Run1;
  tc.fProcess[eProcessSim_Run1] = doprocessSim_Run1;
  tc.fProcess[eProcessTest] = doprocessTest;

  // Temporarary bailout protection against cases which are not implemented/validated yet:
  if (tc.fProcess[eProcessSim]) {
    LOGF(fatal, "\033[1;31m%s at line %d - processSim(...) is not implemented/validated yet \033[0m", __FUNCTION__, __LINE__);
    // TBI 20240512 Most likely, for this case i will have to establish a separate workflow. But since I can with the current
    //              workflow run both over Rec and RecSim, this case is not of a high priority
    // TBI 20240512 See also if I need to extand subscription, both in the definition of CollisionSim and TrackSim
  }

  if (tc.fProcess[eProcessSim_Run2]) {
    LOGF(fatal, "\033[1;31m%s at line %d - processSim_Run2(...) is not implemented/validated yet \033[0m", __FUNCTION__, __LINE__);
    // TBI 20240517 see above comments for eProcessSim , most likely they also apply for this case
  }

  if (tc.fProcess[eProcessRecSim_Run1]) {
    LOGF(fatal, "\033[1;31m%s at line %d - processRecSim_Run1(...) is not implemented/validated yet \033[0m", __FUNCTION__, __LINE__);
  }

  if (tc.fProcess[eProcessSim_Run1]) {
    LOGF(fatal, "\033[1;31m%s at line %d - processSim_Run1(...) is not implemented/validated yet \033[0m", __FUNCTION__, __LINE__);
  }

  // Set automatically generic flags, from above individual flags:
  tc.fProcess[eGenericRec] = tc.fProcess[eProcessRec] || tc.fProcess[eProcessRec_Run2] || tc.fProcess[eProcessRec_Run1] || tc.fProcess[eProcessTest];
  tc.fProcess[eGenericRecSim] = tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessRecSim_Run2] || tc.fProcess[eProcessRecSim_Run1];
  tc.fProcess[eGenericSim] = tc.fProcess[eProcessSim] || tc.fProcess[eProcessSim_Run2] || tc.fProcess[eProcessSim_Run1];

  // Set automatically tc.fWhichProcess from above individual flags:
  if (tc.fProcess[eProcessRec]) {
    tc.fWhichProcess = "ProcessRec";
  } else if (tc.fProcess[eProcessRecSim]) {
    tc.fWhichProcess = "ProcessRecSim";
  } else if (tc.fProcess[eProcessSim]) {
    tc.fWhichProcess = "ProcessSim";
  } else if (tc.fProcess[eProcessRec_Run2]) {
    tc.fWhichProcess = "ProcessRec_Run2";
  } else if (tc.fProcess[eProcessRecSim_Run2]) {
    tc.fWhichProcess = "ProcessRecSim_Run2";
  } else if (tc.fProcess[eProcessSim_Run2]) {
    tc.fWhichProcess = "ProcessSim_Run2";
  } else if (tc.fProcess[eProcessRec_Run1]) {
    tc.fWhichProcess = "ProcessRec_Run1";
  } else if (tc.fProcess[eProcessRecSim_Run1]) {
    tc.fWhichProcess = "ProcessRecSim_Run1";
  } else if (tc.fProcess[eProcessSim_Run1]) {
    tc.fWhichProcess = "ProcessSim_Run1";
  } else if (tc.fProcess[eProcessTest]) {
    tc.fWhichProcess = "ProcessTest";
  }

  tc.fRandomSeed = cf_tc.cfRandomSeed;
  tc.fUseFisherYates = cf_tc.cfUseFisherYates;
  tc.fFixedNumberOfRandomlySelectedTracks = cf_tc.cfFixedNumberOfRandomlySelectedTracks;
  tc.fUseStopwatch = cf_tc.cfUseStopwatch;
  tc.fFloatingPointPrecision = cf_tc.cfFloatingPointPrecision;
  tc.fSequentialBailout = cf_tc.cfSequentialBailout;
  tc.fUseSpecificCuts = cf_tc.cfUseSpecificCuts;
  tc.fWhichSpecificCuts = cf_tc.cfWhichSpecificCuts;

  // *) Event histograms (for QA see below):
  eh.fEventHistogramsName[eNumberOfEvents] = "NumberOfEvents";
  eh.fEventHistogramsName[eTotalMultiplicity] = "TotalMultiplicity";
  eh.fEventHistogramsName[eMultiplicity] = "Multiplicity";
  eh.fEventHistogramsName[eReferenceMultiplicity] = "ReferenceMultiplicity";
  eh.fEventHistogramsName[eCentrality] = "Centrality";
  eh.fEventHistogramsName[eVertex_x] = "Vertex_x";
  eh.fEventHistogramsName[eVertex_y] = "Vertex_y";
  eh.fEventHistogramsName[eVertex_z] = "Vertex_z";
  eh.fEventHistogramsName[eNContributors] = "NContributors";
  eh.fEventHistogramsName[eImpactParameter] = "ImpactParameter";
  eh.fEventHistogramsName[eEventPlaneAngle] = "EventPlaneAngle";
  eh.fEventHistogramsName[eOccupancy] = "Occupancy";
  eh.fEventHistogramsName[eInteractionRate] = "InteractionRate";
  eh.fEventHistogramsName[eCurrentRunDuration] = "CurrentRunDuration";
  eh.fEventHistogramsName[eMultMCNParticlesEta08] = "MultMCNParticlesEta08";

  for (int t = 0; t < eEventHistograms_N; t++) {
    if (eh.fEventHistogramsName[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : name of fEventHistogramsName[%d] is not set \033[0m", __FUNCTION__, __LINE__, static_cast<int>(t));
    }
  }

  // *) Event cuts:
  ec.fUseEventCutCounterAbsolute = cf_ec.cfUseEventCutCounterAbsolute;
  ec.fUseEventCutCounterSequential = cf_ec.cfUseEventCutCounterSequential;
  ec.fPrintCutCounterContent = cf_ec.cfPrintCutCounterContent;

  // Set names of all event cuts:
  ec.fEventCutName[eNumberOfEvents] = "NumberOfEvents";
  ec.fEventCutName[eTotalMultiplicity] = "TotalMultiplicity";
  ec.fEventCutName[eMultiplicity] = "Multiplicity";
  ec.fEventCutName[eReferenceMultiplicity] = "ReferenceMultiplicity";
  ec.fEventCutName[eCentrality] = "Centrality";
  ec.fEventCutName[eVertex_x] = "Vertex_x";
  ec.fEventCutName[eVertex_y] = "Vertex_y";
  ec.fEventCutName[eVertex_z] = "Vertex_z";
  ec.fEventCutName[eNContributors] = "NContributors";
  ec.fEventCutName[eImpactParameter] = "ImpactParameter";
  ec.fEventCutName[eEventPlaneAngle] = "EventPlaneAngle";
  ec.fEventCutName[eOccupancy] = "Occupancy";
  ec.fEventCutName[eInteractionRate] = "InteractionRate";
  ec.fEventCutName[eCurrentRunDuration] = "CurrentRunDuration";
  ec.fEventCutName[eMultMCNParticlesEta08] = "MultMCNParticlesEta08";
  ec.fEventCutName[eTrigger] = "Trigger";
  ec.fEventCutName[eSel7] = "Sel7";
  ec.fEventCutName[eSel8] = "Sel8";
  ec.fEventCutName[eCentralityEstimator] = "CentralityEstimator";
  ec.fEventCutName[eMultiplicityEstimator] = "MultiplicityEstimator";
  ec.fEventCutName[eReferenceMultiplicityEstimator] = "ReferenceMultiplicityEstimator";
  ec.fEventCutName[eSelectedEvents] = "SelectedEvents";
  ec.fEventCutName[eNoSameBunchPileup] = "NoSameBunchPileup";
  ec.fEventCutName[eIsGoodZvtxFT0vsPV] = "IsGoodZvtxFT0vsPV";
  ec.fEventCutName[eIsVertexITSTPC] = "IsVertexITSTPC";
  ec.fEventCutName[eIsVertexTOFmatched] = "IsVertexTOFmatched";
  ec.fEventCutName[eIsVertexTRDmatched] = "IsVertexTRDmatched";
  ec.fEventCutName[eNoCollInTimeRangeStrict] = "NoCollInTimeRangeStrict";
  ec.fEventCutName[eNoCollInTimeRangeStandard] = "NoCollInTimeRangeStandard";
  ec.fEventCutName[eNoCollInRofStrict] = "NoCollInRofStrict";
  ec.fEventCutName[eNoCollInRofStandard] = "NoCollInRofStandard";
  ec.fEventCutName[eNoHighMultCollInPrevRof] = "NoHighMultCollInPrevRof";
  ec.fEventCutName[eIsGoodITSLayer3] = "IsGoodITSLayer3";
  ec.fEventCutName[eIsGoodITSLayer0123] = "IsGoodITSLayer0123";
  ec.fEventCutName[eIsGoodITSLayersAll] = "IsGoodITSLayersAll";
  ec.fEventCutName[eOccupancyEstimator] = "OccupancyEstimator";
  ec.fEventCutName[eMinVertexDistanceFromIP] = "MinVertexDistanceFromIP";
  ec.fEventCutName[eCentralityWeights] = "CentralityWeights";
  for (int t = 0; t < eEventCuts_N; t++) {
    if (ec.fEventCutName[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : event cut name is not set for ec.fEventCutName[%d]. The last cut name which was set is \"%s\" \033[0m", __FUNCTION__, __LINE__, t, ec.fEventCutName[t - 1].Data());
    }
  }

  // *) Particle histograms 1D (for QA see below):
  ph.fParticleHistogramsName[ePhi] = "Phi";
  ph.fParticleHistogramsName[ePt] = "Pt";
  ph.fParticleHistogramsName[eEta] = "Eta";
  ph.fParticleHistogramsName[eCharge] = "Charge";
  ph.fParticleHistogramsName[etpcNClsFindable] = "tpcNClsFindable";
  ph.fParticleHistogramsName[etpcNClsShared] = "tpcNClsShared";
  ph.fParticleHistogramsName[eitsChi2NCl] = "itsChi2NCl";
  ph.fParticleHistogramsName[etpcNClsFound] = "tpcNClsFound";
  ph.fParticleHistogramsName[etpcNClsCrossedRows] = "tpcNClsCrossedRows";
  ph.fParticleHistogramsName[eitsNCls] = "itsNCls";
  ph.fParticleHistogramsName[eitsNClsInnerBarrel] = "itsNClsInnerBarrel";
  ph.fParticleHistogramsName[etpcCrossedRowsOverFindableCls] = "tpcCrossedRowsOverFindableCls";
  ph.fParticleHistogramsName[etpcFoundOverFindableCls] = "tpcFoundOverFindableCls";
  ph.fParticleHistogramsName[etpcFractionSharedCls] = "tpcFractionSharedCls";
  ph.fParticleHistogramsName[etpcChi2NCl] = "tpcChi2NCl";
  ph.fParticleHistogramsName[edcaXY] = "dcaXY";
  ph.fParticleHistogramsName[edcaZ] = "dcaZ";
  ph.fParticleHistogramsName[ePDG] = "PDG";
  for (int t = 0; t < eParticleHistograms_N; t++) {
    if (ph.fParticleHistogramsName[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : name of fParticleHistogramsName[%d] is not set \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // *) Particle histograms 2D (for QA see below):
  ph.fParticleHistogramsName2D[ePhiPt] = Form("%s_vs_%s", ph.fParticleHistogramsName[ePhi].Data(), ph.fParticleHistogramsName[ePt].Data()),
  ph.fParticleHistogramsName2D[ePhiEta] = Form("%s_vs_%s", ph.fParticleHistogramsName[ePhi].Data(), ph.fParticleHistogramsName[eEta].Data());
  for (int t = 0; t < eParticleHistograms2D_N; t++) {
    if (ph.fParticleHistogramsName2D[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : name of fParticleHistogramsName2D[%d] is not set \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // *) Particle cuts:
  pc.fUseParticleCutCounterAbsolute = cf_pc.cfUseParticleCutCounterAbsolute;
  pc.fUseParticleCutCounterSequential = cf_pc.cfUseParticleCutCounterSequential;

  // Set names of all particle cuts:
  pc.fParticleCutName[ePhi] = "Phi";
  pc.fParticleCutName[ePt] = "Pt";
  pc.fParticleCutName[eEta] = "Eta";
  pc.fParticleCutName[eCharge] = "Charge";
  pc.fParticleCutName[etpcNClsFindable] = "tpcNClsFindable";
  pc.fParticleCutName[etpcNClsShared] = "tpcNClsShared";
  pc.fParticleCutName[eitsChi2NCl] = "itsChi2NCl";
  pc.fParticleCutName[etpcNClsFound] = "tpcNClsFound";
  pc.fParticleCutName[etpcNClsCrossedRows] = "tpcNClsCrossedRows";
  pc.fParticleCutName[eitsNCls] = "itsNCls";
  pc.fParticleCutName[eitsNClsInnerBarrel] = "itsNClsInnerBarrel";
  pc.fParticleCutName[etpcCrossedRowsOverFindableCls] = "tpcCrossedRowsOverFindableCls";
  pc.fParticleCutName[etpcFoundOverFindableCls] = "tpcFoundOverFindableCls";
  pc.fParticleCutName[etpcFractionSharedCls] = "tpcFractionSharedCls";
  pc.fParticleCutName[etpcChi2NCl] = "tpcChi2NCl";
  pc.fParticleCutName[edcaXY] = "dcaXY";
  pc.fParticleCutName[edcaZ] = "dcaZ";
  pc.fParticleCutName[ePDG] = "PDG";
  pc.fParticleCutName[etrackCutFlag] = "trackCutFlag";
  pc.fParticleCutName[etrackCutFlagFb1] = "trackCutFlagFb1";
  pc.fParticleCutName[etrackCutFlagFb2] = "trackCutFlagFb2";
  pc.fParticleCutName[eisQualityTrack] = "isQualityTrack";
  pc.fParticleCutName[eisPrimaryTrack] = "isPrimaryTrack";
  pc.fParticleCutName[eisInAcceptanceTrack] = "isInAcceptanceTrack";
  pc.fParticleCutName[eisGlobalTrack] = "isGlobalTrack";
  pc.fParticleCutName[eisPVContributor] = "isPVContributor";
  pc.fParticleCutName[ePtDependentDCAxyParameterization] = "PtDependentDCAxyParameterization";
  for (int t = 0; t < eParticleCuts_N; t++) {
    if (pc.fParticleCutName[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : particle cut name is not set for pc.fParticleCutName[%d] \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // *) Q-vectors:
  qv.fCalculateQvectors = cf_qv.cfCalculateQvectors;

  // *) Multiparticle correlations:
  mupa.fCalculateCorrelations = cf_mupa.cfCalculateCorrelations;
  mupa.fCalculateCorrelationsAsFunctionOf[AFO_INTEGRATED] = cf_mupa.cfCalculateCorrelationsAsFunctionOfIntegrated && mupa.fCalculateCorrelations;
  mupa.fCalculateCorrelationsAsFunctionOf[AFO_MULTIPLICITY] = cf_mupa.cfCalculateCorrelationsAsFunctionOfMultiplicity && mupa.fCalculateCorrelations;
  mupa.fCalculateCorrelationsAsFunctionOf[AFO_CENTRALITY] = cf_mupa.cfCalculateCorrelationsAsFunctionOfCentrality && mupa.fCalculateCorrelations;
  mupa.fCalculateCorrelationsAsFunctionOf[AFO_PT] = cf_mupa.cfCalculateCorrelationsAsFunctionOfPt && mupa.fCalculateCorrelations;
  mupa.fCalculateCorrelationsAsFunctionOf[AFO_ETA] = cf_mupa.cfCalculateCorrelationsAsFunctionOfEta && mupa.fCalculateCorrelations;
  mupa.fCalculateCorrelationsAsFunctionOf[AFO_OCCUPANCY] = cf_mupa.cfCalculateCorrelationsAsFunctionOfOccupancy && mupa.fCalculateCorrelations;
  mupa.fCalculateCorrelationsAsFunctionOf[AFO_INTERACTIONRATE] = cf_mupa.cfCalculateCorrelationsAsFunctionOfInteractionRate && mupa.fCalculateCorrelations;
  mupa.fCalculateCorrelationsAsFunctionOf[AFO_CURRENTRUNDURATION] = cf_mupa.cfCalculateCorrelationsAsFunctionOfCurrentRunDuration && mupa.fCalculateCorrelations;
  mupa.fCalculateCorrelationsAsFunctionOf[AFO_VZ] = cf_mupa.cfCalculateCorrelationsAsFunctionOfVz && mupa.fCalculateCorrelations;

  // *) Test0:
  t0.fCalculateTest0 = cf_t0.cfCalculateTest0;
  t0.fCalculateTest0AsFunctionOf[AFO_INTEGRATED] = cf_t0.cfCalculateTest0AsFunctionOfIntegrated && t0.fCalculateTest0;
  t0.fCalculateTest0AsFunctionOf[AFO_MULTIPLICITY] = cf_t0.cfCalculateTest0AsFunctionOfMultiplicity && t0.fCalculateTest0;
  t0.fCalculateTest0AsFunctionOf[AFO_CENTRALITY] = cf_t0.cfCalculateTest0AsFunctionOfCentrality && t0.fCalculateTest0;
  t0.fCalculateTest0AsFunctionOf[AFO_PT] = cf_t0.cfCalculateTest0AsFunctionOfPt && t0.fCalculateTest0;
  t0.fCalculateTest0AsFunctionOf[AFO_ETA] = cf_t0.cfCalculateTest0AsFunctionOfEta && t0.fCalculateTest0;
  t0.fCalculateTest0AsFunctionOf[AFO_OCCUPANCY] = cf_t0.cfCalculateTest0AsFunctionOfOccupancy && t0.fCalculateTest0;
  t0.fCalculateTest0AsFunctionOf[AFO_INTERACTIONRATE] = cf_t0.cfCalculateTest0AsFunctionOfInteractionRate && t0.fCalculateTest0;
  t0.fCalculateTest0AsFunctionOf[AFO_CURRENTRUNDURATION] = cf_t0.cfCalculateTest0AsFunctionOfCurrentRunDuration && t0.fCalculateTest0;
  t0.fCalculateTest0AsFunctionOf[AFO_VZ] = cf_t0.cfCalculateTest0AsFunctionOfVz && t0.fCalculateTest0;

  if (t0.fCalculateTest0) {
    t0.fFileWithLabels = TString(cf_t0.cfFileWithLabels);
    t0.fUseDefaultLabels = cf_t0.cfUseDefaultLabels;
    t0.fWhichDefaultLabels = TString(cf_t0.cfWhichDefaultLabels);
  }

  // *) Particle weights:
  pw.fUseWeights[wPHI] = cf_pw.cfUsePhiWeights;
  pw.fUseWeights[wPT] = cf_pw.cfUsePtWeights;
  pw.fUseWeights[wETA] = cf_pw.cfUseEtaWeights;
  pw.fUseDiffWeights[wPHIPT] = cf_pw.cfUseDiffPhiPtWeights;   // TBI 20250222 obsolete
  pw.fUseDiffWeights[wPHIETA] = cf_pw.cfUseDiffPhiEtaWeights; // TBI 20250222 obsolete

  // **) Differential phi weights:
  auto lWhichDiffPhiWeights = cf_pw.cfWhichDiffPhiWeights.value;
  if (lWhichDiffPhiWeights.size() != eDiffPhiWeights_N) {
    LOGF(info, "\033[1;31m lWhichDiffPhiWeights.size() = %d\033[0m", lWhichDiffPhiWeights.size());
    LOGF(info, "\033[1;31m eDiffPhiWeights_N = %d\033[0m", static_cast<int>(eDiffPhiWeights_N));
    LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfWhichDiffPhiWeights, and number of entries in enum eDiffPhiWeights_N \n \033[0m", __FUNCTION__, __LINE__);
  }
  for (int dpw = 0; dpw < eDiffPhiWeights_N; dpw++) { // "differential phi weight"
    if (TString(lWhichDiffPhiWeights[dpw]).Contains("wPhi")) {
      pw.fUseDiffPhiWeights[wPhiPhiAxis] = Alright(lWhichDiffPhiWeights[dpw]); // if I pass "1-Phi" => true, "0-Phi" => false
    } else if (TString(lWhichDiffPhiWeights[dpw]).Contains("wPt")) {
      pw.fUseDiffPhiWeights[wPhiPtAxis] = Alright(lWhichDiffPhiWeights[dpw]) && pw.fUseDiffPhiWeights[wPhiPhiAxis]; // I chain here with wPhiPhiAxis , so that I can switch off all differential phi weights, if phi itself is not set to true
    } else if (TString(lWhichDiffPhiWeights[dpw]).Contains("wEta")) {
      pw.fUseDiffPhiWeights[wPhiEtaAxis] = Alright(lWhichDiffPhiWeights[dpw]) && pw.fUseDiffPhiWeights[wPhiPhiAxis];
    } else if (TString(lWhichDiffPhiWeights[dpw]).Contains("wCharge")) {
      pw.fUseDiffPhiWeights[wPhiChargeAxis] = Alright(lWhichDiffPhiWeights[dpw]) && pw.fUseDiffPhiWeights[wPhiPhiAxis];
    } else if (TString(lWhichDiffPhiWeights[dpw]).Contains("wCentrality")) {
      pw.fUseDiffPhiWeights[wPhiCentralityAxis] = Alright(lWhichDiffPhiWeights[dpw]) && pw.fUseDiffPhiWeights[wPhiPhiAxis];
    } else if (TString(lWhichDiffPhiWeights[dpw]).Contains("wVertex_z")) {
      pw.fUseDiffPhiWeights[wPhiVertex_zAxis] = Alright(lWhichDiffPhiWeights[dpw]) && pw.fUseDiffPhiWeights[wPhiPhiAxis];
    } else {
      LOGF(fatal, "\033[1;31m%s at line %d : The setting %s in configurable cfWhichDiffPhiWeights is not supported yet. See enum eDiffPhiWeights . \n \033[0m", __FUNCTION__, __LINE__, TString(lWhichDiffPhiWeights[dpw]).Data());
    }
  }

  // **) Differential pt weights:
  auto lWhichDiffPtWeights = cf_pw.cfWhichDiffPtWeights.value;
  if (lWhichDiffPtWeights.size() != eDiffPtWeights_N) {
    LOGF(info, "\033[1;31m lWhichDiffPtWeights.size() = %d\033[0m", lWhichDiffPtWeights.size());
    LOGF(info, "\033[1;31m eDiffPtWeights_N = %d\033[0m", static_cast<int>(eDiffPtWeights_N));
    LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfWhichDiffPtWeights, and number of entries in enum eDiffPtWeights_N \n \033[0m", __FUNCTION__, __LINE__);
  }
  for (int dpw = 0; dpw < eDiffPtWeights_N; dpw++) { // "differential pt weight"
    if (TString(lWhichDiffPtWeights[dpw]).Contains("wPt")) {
      pw.fUseDiffPtWeights[wPtPtAxis] = Alright(lWhichDiffPtWeights[dpw]); // if I pass "1-Pt" => true, "0-Pt" => false
    } else {                                                               // ... TBI 20250222 add support for other dimensions of differential pt weights, in the same spirit i did it for differential phi weights
      LOGF(fatal, "\033[1;31m%s at line %d : The setting %s in configurable cfWhichDiffPtWeights is not supported yet. See enum eDiffPtWeights . \n \033[0m", __FUNCTION__, __LINE__, TString(lWhichDiffPtWeights[dpw]).Data());
    }
  }

  // **) Differential eta weights:
  auto lWhichDiffEtaWeights = cf_pw.cfWhichDiffEtaWeights.value;
  if (lWhichDiffEtaWeights.size() != eDiffEtaWeights_N) {
    LOGF(info, "\033[1;31m lWhichDiffEtaWeights.size() = %d\033[0m", lWhichDiffEtaWeights.size());
    LOGF(info, "\033[1;31m eDiffEtaWeights_N = %d\033[0m", static_cast<int>(eDiffEtaWeights_N));
    LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfWhichDiffEtaWeights, and number of entries in enum eDiffEtaWeights_N \n \033[0m", __FUNCTION__, __LINE__);
  }
  for (int dpw = 0; dpw < eDiffEtaWeights_N; dpw++) { // "differential eta weight"
    if (TString(lWhichDiffEtaWeights[dpw]).Contains("wEta")) {
      pw.fUseDiffEtaWeights[wEtaEtaAxis] = Alright(lWhichDiffEtaWeights[dpw]); // if I pass "1-Eta" => true, "0-Eta" => false
    } else {                                                                   // ... TBI 20250222 add support for other dimensions of differential eta weights, in the same spirit i did it for differential phi weights
      LOGF(fatal, "\033[1;31m%s at line %d : The setting %s in configurable cfWhichDiffEtaWeights is not supported yet. See enum eDiffEtaWeights . \n \033[0m", __FUNCTION__, __LINE__, TString(lWhichDiffEtaWeights[dpw]).Data());
    }
  }

  // **) File holding all particle weights:
  pw.fFileWithWeights = cf_pw.cfFileWithWeights;

  // *) Centrality weights:
  cw.fUseCentralityWeights = cf_cw.cfUseCentralityWeights;
  cw.fFileWithCentralityWeights = cf_cw.cfFileWithCentralityWeights;

  // ...

  // *) Nested loops:
  nl.fCalculateNestedLoops = cf_nl.cfCalculateNestedLoops;
  nl.fCalculateCustomNestedLoops = cf_nl.cfCalculateCustomNestedLoops;
  nl.fCalculateKineCustomNestedLoops = cf_nl.cfCalculateKineCustomNestedLoops;
  nl.fMaxNestedLoop = cf_nl.cfMaxNestedLoop;

  // ...

  // *) Toy NUA:
  auto lApplyNUAPDF = (vector<int>)cf_nua.cfApplyNUAPDF;
  if (lApplyNUAPDF.size() != eNUAPDF_N) {
    LOGF(info, "\033[1;31m lApplyNUAPDF.size() = %d\033[0m", lApplyNUAPDF.size());
    LOGF(info, "\033[1;31m eNUAPDF_N = %d\033[0m", static_cast<int>(eNUAPDF_N));
    LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfApplyNUAPDF, and number of entries in enum eNUAPDF \n \033[0m", __FUNCTION__, __LINE__);
  }
  nua.fApplyNUAPDF[ePhiNUAPDF] = static_cast<bool>(lApplyNUAPDF[ePhiNUAPDF]);
  nua.fApplyNUAPDF[ePtNUAPDF] = static_cast<bool>(lApplyNUAPDF[ePtNUAPDF]);
  nua.fApplyNUAPDF[eEtaNUAPDF] = static_cast<bool>(lApplyNUAPDF[eEtaNUAPDF]);

  // **) Execute the lines below, only if toy NUA (either default or custom) is requested for at least one kine variable:
  if (nua.fApplyNUAPDF[ePhiNUAPDF] || nua.fApplyNUAPDF[ePtNUAPDF] || nua.fApplyNUAPDF[eEtaNUAPDF]) {

    auto lUseDefaultNUAPDF = (vector<int>)cf_nua.cfUseDefaultNUAPDF;
    if (lUseDefaultNUAPDF.size() != eNUAPDF_N) {
      LOGF(info, "\033[1;31m lUseDefaultNUAPDF.size() = %d\033[0m", lUseDefaultNUAPDF.size());
      LOGF(info, "\033[1;31m eNUAPDF_N = %d\033[0m", static_cast<int>(eNUAPDF_N));
      LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfUseDefaultNUAPDF, and number of entries in enum eNUAPDF \n \033[0m", __FUNCTION__, __LINE__);
    }
    nua.fUseDefaultNUAPDF[ePhiNUAPDF] = static_cast<bool>(lUseDefaultNUAPDF[ePhiNUAPDF]);
    nua.fUseDefaultNUAPDF[ePtNUAPDF] = static_cast<bool>(lUseDefaultNUAPDF[ePtNUAPDF]);
    nua.fUseDefaultNUAPDF[eEtaNUAPDF] = static_cast<bool>(lUseDefaultNUAPDF[eEtaNUAPDF]);

    // **) Execute the lines below, only if custom toy NUA is requested in at least one kine variable:
    if (!((nua.fApplyNUAPDF[ePhiNUAPDF] && nua.fUseDefaultNUAPDF[ePhiNUAPDF]) ||
          (nua.fApplyNUAPDF[ePtNUAPDF] && nua.fUseDefaultNUAPDF[ePtNUAPDF]) ||
          (nua.fApplyNUAPDF[eEtaNUAPDF] && nua.fUseDefaultNUAPDF[eEtaNUAPDF]))) {
      // If the above conditon is true, as least one NUA is requested and is not default, i.e. it's custom NUA obtained from external file, which was requested to be used.
      // TBI 20240501 Can I simplify the logic above, it's a bit cryptic...

      // *) external file path with custom NUA histos:
      nua.fFileWithCustomNUA = TString(cf_nua.cfFileWithCustomNUA);

      // *) histogram names with custom NUA distributions in that file + get those histograms immediately here:
      auto lCustomNUAPDFHistNames = (vector<string>)cf_nua.cfCustomNUAPDFHistNames;
      // TBI 20241115 For some reason, the default values of configurable "cfCustomNUAPDFHistNames" are not correctly propagated in the local variables, but I can circumvent that with JSON settings for the time being
      if (lCustomNUAPDFHistNames.size() != eNUAPDF_N) {
        LOGF(info, "\033[1;31m lCustomNUAPDFHistNames.size() = %d\033[0m", lCustomNUAPDFHistNames.size());
        LOGF(info, "\033[1;31m eNUAPDF_N = %d\033[0m", static_cast<int>(eNUAPDF_N));
        LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfCustomNUAPDFHistNames, and number of entries in enum eNUAPDF \n \033[0m", __FUNCTION__, __LINE__);
      }

      if (!nua.fUseDefaultNUAPDF[ePhiNUAPDF]) {
        nua.fCustomNUAPDFHistNames[ePhiNUAPDF] = new TString(lCustomNUAPDFHistNames[ePhiNUAPDF]);
        this->GetHistogramWithCustomNUA(nua.fFileWithCustomNUA.Data(), ePhiNUAPDF);
        if (!nua.fCustomNUAPDF[ePhiNUAPDF]) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }

      if (!nua.fUseDefaultNUAPDF[ePtNUAPDF]) {
        nua.fCustomNUAPDFHistNames[ePtNUAPDF] = new TString(lCustomNUAPDFHistNames[ePtNUAPDF]);
        this->GetHistogramWithCustomNUA(nua.fFileWithCustomNUA.Data(), ePtNUAPDF);
        if (!nua.fCustomNUAPDF[ePtNUAPDF]) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }
      if (!nua.fUseDefaultNUAPDF[eEtaNUAPDF]) {
        nua.fCustomNUAPDFHistNames[eEtaNUAPDF] = new TString(lCustomNUAPDFHistNames[eEtaNUAPDF]);
        this->GetHistogramWithCustomNUA(nua.fFileWithCustomNUA.Data(), eEtaNUAPDF);
        if (!nua.fCustomNUAPDF[eEtaNUAPDF]) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }
    } // if (!(nua.fUseDefaultNUAPDF[ePhiNUAPDF] || nua.fUseDefaultNUAPDF[ePtNUAPDF] || nua.fUseDefaultNUAPDF[eEtaNUAPDF])) {

  } // if ( nua.fApplyNUAPDF[ePhiNUAPDF] || nua.fApplyNUAPDF[ePtNUAPDF] || nua.fApplyNUAPDF[eEtaNUAPDF] ) {

  // *) Internal validation:
  iv.fUseInternalValidation = cf_iv.cfUseInternalValidation;
  iv.fInternalValidationForceBailout = cf_iv.cfInternalValidationForceBailout;
  iv.fnEventsInternalValidation = cf_iv.cfnEventsInternalValidation;
  iv.fRescaleWithTheoreticalInput = cf_iv.cfRescaleWithTheoreticalInput;
  iv.fHarmonicsOptionInternalValidation = new TString(cf_iv.cfHarmonicsOptionInternalValidation);

  // *) Results histograms:
  //    Define axis titles:
  //    Remark: keep ordering in sync with enum eAsFunctionOf
  res.fResultsProXaxisTitle[AFO_INTEGRATED] = "integrated";
  res.fResultsProRawName[AFO_INTEGRATED] = "int"; // this is how it appears simplified in the hist name when saved to the file
  res.fResultsProXaxisTitle[AFO_MULTIPLICITY] = "multiplicity";
  res.fResultsProRawName[AFO_MULTIPLICITY] = "mult";
  res.fResultsProXaxisTitle[AFO_CENTRALITY] = "centrality";
  res.fResultsProRawName[AFO_CENTRALITY] = "cent";
  res.fResultsProXaxisTitle[AFO_PT] = "pt";
  res.fResultsProRawName[AFO_PT] = "pt";
  res.fResultsProXaxisTitle[AFO_ETA] = "eta";
  res.fResultsProRawName[AFO_ETA] = "eta";
  res.fResultsProXaxisTitle[AFO_OCCUPANCY] = "occupancy";
  res.fResultsProRawName[AFO_OCCUPANCY] = "occu";
  res.fResultsProXaxisTitle[AFO_INTERACTIONRATE] = "interaction rate";
  res.fResultsProRawName[AFO_INTERACTIONRATE] = "ir";
  res.fResultsProXaxisTitle[AFO_CURRENTRUNDURATION] = "current run duration";
  res.fResultsProRawName[AFO_CURRENTRUNDURATION] = "crd";
  res.fResultsProXaxisTitle[AFO_VZ] = "vertex z position";
  res.fResultsProRawName[AFO_VZ] = "vz";
  res.fSaveResultsHistograms = cf_res.cfSaveResultsHistograms;

  // *) QA:
  //    Remark: I keep it on the bottom, because here I define some names in temrs of names defined above.
  qa.fCheckUnderflowAndOverflow = cf_qa.cfCheckUnderflowAndOverflow;
  qa.fRebin = cf_qa.cfRebin;

  // **) Reference multiplicity estimators:
  qa.fReferenceMultiplicityEstimatorName[eMultTPC] = "MultTPC";
  qa.fReferenceMultiplicityEstimatorName[eMultFV0M] = "MultFV0M";
  qa.fReferenceMultiplicityEstimatorName[eMultFT0C] = "MultFT0C";
  qa.fReferenceMultiplicityEstimatorName[eMultFT0M] = "MultFT0M";
  qa.fReferenceMultiplicityEstimatorName[eMultNTracksPV] = "MultNTracksPV";
  qa.fReferenceMultiplicityEstimatorName[eMultNTracksGlobal] = "MultNTracksGlobal";
  qa.fReferenceMultiplicityEstimatorName[eMultTracklets] = "MultTracklets";

  // **) Centrality estimators:
  qa.fCentralityEstimatorName[eCentFT0C] = "CentFT0C";
  qa.fCentralityEstimatorName[eCentFT0CVariant1] = "CentFT0CVariant1";
  qa.fCentralityEstimatorName[eCentFT0M] = "CentFT0M";
  qa.fCentralityEstimatorName[eCentFV0A] = "CentFV0A";
  qa.fCentralityEstimatorName[eCentNTPV] = "CentNTPV";
  qa.fCentralityEstimatorName[eCentNGlobal] = "CentNGlobal";
  qa.fCentralityEstimatorName[eCentRun2V0M] = "CentRun2V0M";
  qa.fCentralityEstimatorName[eCentRun2SPDTracklets] = "CentRun2SPDTracklets";

  // **) Occupancy estimators:
  qa.fOccupancyEstimatorName[eTrackOccupancyInTimeRange] = "TrackOccupancyInTimeRange";
  qa.fOccupancyEstimatorName[eFT0COccupancyInTimeRange] = "FT0COccupancyInTimeRange";

  // **) Names of QA 2D event histograms:
  //     Remark: Do NOT use FancyFormatting here, only later in BookQAHistograms() for axis titles!
  qa.fEventHistogramsName2D[eMultiplicity_vs_ReferenceMultiplicity] = Form("%s_vs_%s", eh.fEventHistogramsName[eMultiplicity].Data(), eh.fEventHistogramsName[eReferenceMultiplicity].Data());
  qa.fEventHistogramsName2D[eMultiplicity_vs_NContributors] = Form("%s_vs_%s", eh.fEventHistogramsName[eMultiplicity].Data(), eh.fEventHistogramsName[eNContributors].Data());
  qa.fEventHistogramsName2D[eMultiplicity_vs_Centrality] = Form("%s_vs_%s", eh.fEventHistogramsName[eMultiplicity].Data(), eh.fEventHistogramsName[eCentrality].Data());
  qa.fEventHistogramsName2D[eMultiplicity_vs_Vertex_z] = Form("%s_vs_%s", eh.fEventHistogramsName[eMultiplicity].Data(), eh.fEventHistogramsName[eVertex_z].Data());
  qa.fEventHistogramsName2D[eMultiplicity_vs_Occupancy] = Form("%s_vs_%s", eh.fEventHistogramsName[eMultiplicity].Data(), eh.fEventHistogramsName[eOccupancy].Data());
  qa.fEventHistogramsName2D[eReferenceMultiplicity_vs_NContributors] = Form("%s_vs_%s", eh.fEventHistogramsName[eReferenceMultiplicity].Data(), eh.fEventHistogramsName[eNContributors].Data());
  qa.fEventHistogramsName2D[eReferenceMultiplicity_vs_Centrality] = Form("%s_vs_%s", eh.fEventHistogramsName[eReferenceMultiplicity].Data(), eh.fEventHistogramsName[eCentrality].Data());
  qa.fEventHistogramsName2D[eReferenceMultiplicity_vs_Vertex_z] = Form("%s_vs_%s", eh.fEventHistogramsName[eReferenceMultiplicity].Data(), eh.fEventHistogramsName[eVertex_z].Data());
  qa.fEventHistogramsName2D[eReferenceMultiplicity_vs_Occupancy] = Form("%s_vs_%s", eh.fEventHistogramsName[eReferenceMultiplicity].Data(), eh.fEventHistogramsName[eOccupancy].Data());
  qa.fEventHistogramsName2D[eNContributors_vs_Centrality] = Form("%s_vs_%s", eh.fEventHistogramsName[eNContributors].Data(), eh.fEventHistogramsName[eCentrality].Data());
  qa.fEventHistogramsName2D[eNContributors_vs_Vertex_z] = Form("%s_vs_%s", eh.fEventHistogramsName[eNContributors].Data(), eh.fEventHistogramsName[eVertex_z].Data());
  qa.fEventHistogramsName2D[eNContributors_vs_Occupancy] = Form("%s_vs_%s", eh.fEventHistogramsName[eNContributors].Data(), eh.fEventHistogramsName[eOccupancy].Data());
  qa.fEventHistogramsName2D[eCentrality_vs_Vertex_z] = Form("%s_vs_%s", eh.fEventHistogramsName[eCentrality].Data(), eh.fEventHistogramsName[eVertex_z].Data());
  qa.fEventHistogramsName2D[eCentrality_vs_Occupancy] = Form("%s_vs_%s", eh.fEventHistogramsName[eCentrality].Data(), eh.fEventHistogramsName[eOccupancy].Data());
  qa.fEventHistogramsName2D[eCentrality_vs_ImpactParameter] = Form("%s_vs_%s", eh.fEventHistogramsName[eCentrality].Data(), eh.fEventHistogramsName[eImpactParameter].Data());
  qa.fEventHistogramsName2D[eVertex_z_vs_Occupancy] = Form("%s_vs_%s", eh.fEventHistogramsName[eVertex_z].Data(), eh.fEventHistogramsName[eOccupancy].Data());
  qa.fEventHistogramsName2D[eMultNTracksPV_vs_MultNTracksGlobal] = Form("%s_vs_%s", qa.fReferenceMultiplicityEstimatorName[eMultNTracksPV].Data(), qa.fReferenceMultiplicityEstimatorName[eMultNTracksGlobal].Data());
  qa.fEventHistogramsName2D[eCentFT0C_vs_CentFT0CVariant1] = Form("%s_vs_%s", qa.fCentralityEstimatorName[eCentFT0C].Data(), qa.fCentralityEstimatorName[eCentFT0CVariant1].Data());
  qa.fEventHistogramsName2D[eCentFT0C_vs_CentFT0M] = Form("%s_vs_%s", qa.fCentralityEstimatorName[eCentFT0C].Data(), qa.fCentralityEstimatorName[eCentFT0M].Data());
  qa.fEventHistogramsName2D[eCentFT0C_vs_CentFV0A] = Form("%s_vs_%s", qa.fCentralityEstimatorName[eCentFT0C].Data(), qa.fCentralityEstimatorName[eCentFV0A].Data());
  qa.fEventHistogramsName2D[eCentFT0C_vs_CentNTPV] = Form("%s_vs_%s", qa.fCentralityEstimatorName[eCentFT0C].Data(), qa.fCentralityEstimatorName[eCentNTPV].Data());
  qa.fEventHistogramsName2D[eCentFT0C_vs_CentNGlobal] = Form("%s_vs_%s", qa.fCentralityEstimatorName[eCentFT0C].Data(), qa.fCentralityEstimatorName[eCentNGlobal].Data());
  qa.fEventHistogramsName2D[eCentFT0M_vs_CentNTPV] = Form("%s_vs_%s", qa.fCentralityEstimatorName[eCentFT0M].Data(), qa.fCentralityEstimatorName[eCentNTPV].Data());
  qa.fEventHistogramsName2D[eCentRun2V0M_vs_CentRun2SPDTracklets] = Form("%s_vs_%s", qa.fCentralityEstimatorName[eCentRun2V0M].Data(), qa.fCentralityEstimatorName[eCentRun2SPDTracklets].Data());
  qa.fEventHistogramsName2D[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = Form("%s_vs_%s", qa.fOccupancyEstimatorName[eTrackOccupancyInTimeRange].Data(), qa.fOccupancyEstimatorName[eFT0COccupancyInTimeRange].Data());
  qa.fEventHistogramsName2D[eCurrentRunDuration_vs_InteractionRate] = Form("%s_vs_%s", ec.fEventCutName[eCurrentRunDuration].Data(), ec.fEventCutName[eInteractionRate].Data());

  // ***) Quick insanity check that all names are set:
  for (int t = 0; t < eQAEventHistograms2D_N; t++) {
    if (qa.fEventHistogramsName2D[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : qa.fEventHistogramsName2D[%d] is not set, check corresponding enum eQAEventHistograms2D \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // **) Names of QA 2D particle histograms:
  qa.fParticleHistogramsName2D[ePt_vs_dcaXY] = Form("%s_vs_%s", ph.fParticleHistogramsName[ePt].Data(), ph.fParticleHistogramsName[edcaXY].Data());

  // ***) Quick insanity check that all names are set:
  for (int t = 0; t < eQAParticleHistograms2D_N; t++) {
    if (qa.fParticleHistogramsName2D[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : qa.fParticleHistogramsName2D[%d] is not set, check corresponding enum eQAParticleHistograms2D \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // **) Names of QA 2D particle event histograms:
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_itsNClsEbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), ph.fParticleHistogramsName[eitsNCls].Data()).Data();
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), TString(ph.fParticleHistogramsName[eitsNCls].Data()).Append("NegEtaEbyE").Data()).Data(); // TBI 20241214 time will tell if this Append() is safe enough... Remember that Append works in-place
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), TString(ph.fParticleHistogramsName[eitsNCls].Data()).Append("PosEtaEbyE").Data()).Data(); // TBI 20241214 time will tell if this Append() is safe enough... Remember that Append works in-place
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_Eta0804EbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), TString(ph.fParticleHistogramsName[eEta].Data()).Append("0804EbyE").Data()).Data();             // TBI 20241214 time will tell if this Append() is safe enough... Remember that Append works in-place
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_Eta0400EbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), TString(ph.fParticleHistogramsName[eEta].Data()).Append("0400EbyE").Data()).Data();             // TBI 20241214 time will tell if this Append() is safe enough... Remember that Append works in-place
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_Eta0004EbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), TString(ph.fParticleHistogramsName[eEta].Data()).Append("0004EbyE").Data()).Data();             // TBI 20241214 time will tell if this Append() is safe enough... Remember that Append works in-place
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_Eta0408EbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), TString(ph.fParticleHistogramsName[eEta].Data()).Append("0408EbyE").Data()).Data();             // TBI 20241214 time will tell if this Append() is safe enough... Remember that Append works in-place
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_Pt0005EbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), TString(ph.fParticleHistogramsName[ePt].Data()).Append("0005EbyE").Data()).Data();               // TBI 20241214 time will tell if this Append() is safe enough... Remember that Append works in-place
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_Pt0510EbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), TString(ph.fParticleHistogramsName[ePt].Data()).Append("0510EbyE").Data()).Data();               // TBI 20241214 time will tell if this Append() is safe enough... Remember that Append works in-place
  qa.fQAParticleEventHistogramsName2D[eCurrentRunDuration_vs_Pt1050EbyE] = TString::Format("%s_vs_%s", eh.fEventHistogramsName[eCurrentRunDuration].Data(), TString(ph.fParticleHistogramsName[ePt].Data()).Append("1050EbyE").Data()).Data();               // TBI 20241214 time will tell if this Append() is safe enough... Remember that Append works in-place

  // ***) Quick insanity check that all names are set:
  for (int t = 0; t < eQAParticleEventHistograms2D_N; t++) {
    if (qa.fQAParticleEventHistogramsName2D[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : qa.fQAParticleEventHistogramsName2D[%d] is not set, check corresponding enum eQAParticleEventHistograms2D \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // **) Names of QA 2D "correlations vs." histograms:
  qa.fQACorrelationsVsHistogramsName2D[eCorrelations_vs_Multiplicity] = TString::Format("%s_vs_%s", "Correlations", eh.fEventHistogramsName[eMultiplicity].Data()).Data();
  qa.fQACorrelationsVsHistogramsName2D[eCorrelations_vs_ReferenceMultiplicity] = TString::Format("%s_vs_%s", "Correlations", eh.fEventHistogramsName[eReferenceMultiplicity].Data()).Data();
  qa.fQACorrelationsVsHistogramsName2D[eCorrelations_vs_Centrality] = TString::Format("%s_vs_%s", "Correlations", eh.fEventHistogramsName[eCentrality].Data()).Data();
  // ...
  qa.fQACorrelationsVsHistogramsName2D[eCorrelations_vs_MeanPhi] = TString::Format("%s_vs_%s", "Correlations", ph.fParticleHistogramsName[ePhi].Data()).Data();
  qa.fQACorrelationsVsHistogramsName2D[eCorrelations_vs_MeanPt] = TString::Format("%s_vs_%s", "Correlations", ph.fParticleHistogramsName[ePt].Data()).Data();
  qa.fQACorrelationsVsHistogramsName2D[eCorrelations_vs_MeanEta] = TString::Format("%s_vs_%s", "Correlations", ph.fParticleHistogramsName[eEta].Data()).Data();
  // ...

  // ***) Quick insanity check that all names are set:
  for (int t = 0; t < eQACorrelationsVsHistograms2D_N; t++) {
    if (qa.fQACorrelationsVsHistogramsName2D[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : qa.fQACorrelationsVsHistogramsName2D[%d] is not set, check corresponding enum eQACorrelationsVsHistograms2D \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // **) Names and titles of all categories of sparse histograms:
  ph.fParticleSparseHistogramsName[eDWPhi] = "fParticleSparseHistograms_DWPhi";
  ph.fParticleSparseHistogramsTitle[eDWPhi] = "sparse histogram for differential #phi weights,";

  ph.fParticleSparseHistogramsName[eDWPt] = "fParticleSparseHistograms_DWPt";
  ph.fParticleSparseHistogramsTitle[eDWPt] = "sparse histogram for differential p_{T} weights,";

  ph.fParticleSparseHistogramsName[eDWEta] = "fParticleSparseHistograms_DWEta";
  ph.fParticleSparseHistogramsTitle[eDWEta] = "sparse histogram for differential #eta weights,";

  // ...

  // ** Eta separations:
  es.fCalculateEtaSeparations = cf_es.cfCalculateEtaSeparations;
  es.fCalculateEtaSeparationsAsFunctionOf[AFO_INTEGRATED] = cf_es.cfCalculateEtaSeparationsAsFunctionOfIntegrated && es.fCalculateEtaSeparations;
  es.fCalculateEtaSeparationsAsFunctionOf[AFO_MULTIPLICITY] = cf_es.cfCalculateEtaSeparationsAsFunctionOfMultiplicity && es.fCalculateEtaSeparations;
  es.fCalculateEtaSeparationsAsFunctionOf[AFO_CENTRALITY] = cf_es.cfCalculateEtaSeparationsAsFunctionOfCentrality && es.fCalculateEtaSeparations;
  es.fCalculateEtaSeparationsAsFunctionOf[AFO_PT] = cf_es.cfCalculateEtaSeparationsAsFunctionOfPt && es.fCalculateEtaSeparations;
  es.fCalculateEtaSeparationsAsFunctionOf[AFO_ETA] = false; // this one doesn't make sense in this context, obviously
  es.fCalculateEtaSeparationsAsFunctionOf[AFO_OCCUPANCY] = cf_es.cfCalculateEtaSeparationsAsFunctionOfOccupancy && es.fCalculateEtaSeparations;
  es.fCalculateEtaSeparationsAsFunctionOf[AFO_INTERACTIONRATE] = cf_es.cfCalculateEtaSeparationsAsFunctionOfInteractionRate && es.fCalculateEtaSeparations;
  es.fCalculateEtaSeparationsAsFunctionOf[AFO_CURRENTRUNDURATION] = cf_es.cfCalculateEtaSeparationsAsFunctionOfCurrentRunDuration && es.fCalculateEtaSeparations;
  es.fCalculateEtaSeparationsAsFunctionOf[AFO_VZ] = cf_es.cfCalculateEtaSeparationsAsFunctionOfVz && es.fCalculateEtaSeparations;

  if (es.fCalculateEtaSeparations) {
    auto lEtaSeparationsValues = cf_es.cfEtaSeparationsValues.value;
    if (lEtaSeparationsValues.size() != gMaxNumberEtaSeparations) {
      LOGF(info, "\033[1;31m%s at line %d : lEtaSeparationsValues.size() = %d\n \033[0m", __FUNCTION__, __LINE__, lEtaSeparationsValues.size());
      LOGF(fatal, "\033[1;31m%s at line %d : Provide in configurable cfEtaSeparationsValues precisely %d entries\n \033[0m", __FUNCTION__, __LINE__, static_cast<int>(gMaxNumberEtaSeparations));
    }
    for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
      if (lEtaSeparationsValues[e] < 0.) {
        LOGF(fatal, "\033[1;31m%s at line %d : lEtaSeparationsValues[%d] = %f is not >= 0. \n \033[0m", __FUNCTION__, __LINE__, e, static_cast<float>(lEtaSeparationsValues[e]));
      }
      es.fEtaSeparationsValues[e] = lEtaSeparationsValues[e];
    }

    auto lEtaSeparationsSkipHarmonics = cf_es.cfEtaSeparationsSkipHarmonics.value;
    if (lEtaSeparationsSkipHarmonics.size() != gMaxHarmonic) {
      LOGF(info, "\033[1;31m lEtaSeparationsSkipHarmonics.size() = %d\033[0m", lEtaSeparationsSkipHarmonics.size());
      LOGF(info, "\033[1;31m gMaxHarmonic) = %d\033[0m", static_cast<int>(gMaxHarmonic));
      LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfEtaSeparationsSkipHarmonics, and max number of supported harmonics \n \033[0m", __FUNCTION__, __LINE__);
    }

    for (int h = 0; h < static_cast<int>(lEtaSeparationsSkipHarmonics.size()); h++) {
      es.fEtaSeparationsSkipHarmonics[h] = Alright(lEtaSeparationsSkipHarmonics[h]);
    }

  } // if(es.fCalculateEtaSeparations) {

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void DefaultConfiguration()

//============================================================

bool Alright(TString s)
{
  // Simple utility function, which for a string formatted "0-someName" returns false, and for "1-someName" returns true.

  // a) Insanity check on the format;
  // b) Do the thing.

  if (tc.fVerboseUtility) {
    StartFunction(__FUNCTION__);
    LOGF(info, "\033[1;32m  TString s = %s\033[0m", s.Data());
  }

  bool returnValue = false;

  // a) Insanity check on the format:
  TObjArray* oa = s.Tokenize("-");
  if (!oa) {
    LOGF(fatal, "\033[1;31m%s at line %d : oa is NULL , s = %s\033[0m", __FUNCTION__, __LINE__, s.Data());
  }
  int nEntries = oa->GetEntries();
  if (2 != nEntries) {
    LOGF(fatal, "\033[1;31m%s at line %d : string expected in this function must be formatted as \"someName-0\" or \"someName-1\" => s = %s\033[0m", __FUNCTION__, __LINE__, s.Data());
  }

  // b) Do the thing:
  //    Algorithm: I split "0-someName" or "1-someName" with respect to "-" as a field separator, and check what is in the 1st field.
  if (TString(oa->At(0)->GetName()).EqualTo("0")) {
    delete oa;
    returnValue = false;
  } else if (TString(oa->At(0)->GetName()).EqualTo("1")) {
    delete oa;
    returnValue = true;
  } else {
    LOGF(fatal, "\033[1;31m%s at line %d : string expected in this function must be formatted as \"0-someName\" or \"1-someName\" => s = %s\033[0m", __FUNCTION__, __LINE__, s.Data());
  }

  if (tc.fVerboseUtility) {
    ExitFunction(__FUNCTION__);
  }

  return returnValue;

} // bool Alright(const char* name)

//============================================================

void DefaultBooking()
{
  // Set here which histograms are booked by default.

  // a) Event histograms 1D;
  // b) Event histograms 2D;
  // c) Particle histograms 1D;
  // d) Particle histograms 2D;
  // e) Particle sparse histograms;
  // f) QA.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Event histograms 1D:
  // By default all event histograms are booked. Set this flag to false to switch off booking of all event histograms:
  eh.fFillEventHistograms = cf_eh.cfFillEventHistograms;

  // *) By default all event histograms are booked. If you do not want particular event histogram to be booked,
  // use configurable array cfBookEventHistograms, where you can specify name of the histogram accompanied with flags 1 (book) or 0 (do not book).
  // Supported format: "someName-0" and "someName-1", where "-" is a field separator.
  // Ordering of the flags in that array is interpreted through ordering of enums in enum eEventHistograms.
  auto lBookEventHistograms = cf_eh.cfBookEventHistograms.value; // this is now the local version of that string array from configurable.
  if (lBookEventHistograms.size() != eEventHistograms_N) {
    LOGF(info, "\033[1;31m lBookEventHistograms.size() = %d\033[0m", lBookEventHistograms.size());
    LOGF(info, "\033[1;31m eEventHistograms_N) = %d\033[0m", static_cast<int>(eEventHistograms_N));
    LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfBookEventHistograms, and number of entries in enum eEventHistograms \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering of histogram names in the initialization in configurable cfBookEventHistograms:
  // TBI 20240518 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  for (int name = 0; name < eEventHistograms_N; name++) {
    // TBI 20240518 I could implement even a strickter EqualTo instead of EndsWith, but then I need to tokenize, etc., etc. This shall be safe enough.
    if (!TString(lBookEventHistograms[name]).EndsWith(eh.fEventHistogramsName[name].Data())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfBookEventHistograms => name = %d, lBookEventHistograms[%d] = \"%s\", eh.fEventHistogramsName[%d] = \"%s\" \n Check if you are using an up to date tag. \033[0m", __FUNCTION__, __LINE__, name, name, TString(lBookEventHistograms[name]).Data(), name, eh.fEventHistogramsName[name].Data());
    }
  }

  // I append "&& eh.fFillEventHistograms" below, to switch off booking of all event histograms with one common flag:
  eh.fBookEventHistograms[eNumberOfEvents] = Alright(lBookEventHistograms[eNumberOfEvents]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eTotalMultiplicity] = Alright(lBookEventHistograms[eTotalMultiplicity]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eMultiplicity] = Alright(lBookEventHistograms[eMultiplicity]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eReferenceMultiplicity] = Alright(lBookEventHistograms[eReferenceMultiplicity]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eCentrality] = Alright(lBookEventHistograms[eCentrality]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eVertex_x] = Alright(lBookEventHistograms[eVertex_x]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eVertex_y] = Alright(lBookEventHistograms[eVertex_y]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eVertex_z] = Alright(lBookEventHistograms[eVertex_z]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eNContributors] = Alright(lBookEventHistograms[eNContributors]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eImpactParameter] = Alright(lBookEventHistograms[eImpactParameter]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eEventPlaneAngle] = Alright(lBookEventHistograms[eEventPlaneAngle]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eOccupancy] = Alright(lBookEventHistograms[eOccupancy]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eInteractionRate] = Alright(lBookEventHistograms[eInteractionRate]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eCurrentRunDuration] = Alright(lBookEventHistograms[eCurrentRunDuration]) && eh.fFillEventHistograms;
  eh.fBookEventHistograms[eMultMCNParticlesEta08] = Alright(lBookEventHistograms[eMultMCNParticlesEta08]) && eh.fFillEventHistograms;

  // b) Event histograms 2D:
  // TBI 20240515 Ideally, all 2D shall go to QA group, see below
  // ...

  // c) Particle histograms 1D:
  // By default all 1D particle histograms are booked. Set this flag to false to switch off booking of all 1D particle histograms:
  ph.fFillParticleHistograms = cf_ph.cfFillParticleHistograms;

  // *) If you do not want particular particle histogram to be booked, use configurable array cfBookParticleHistograms, where you can specify flags 1 (book) or 0 (do not book).
  // Ordering of the flags in that array is interpreted through ordering of enums in enum eParticleHistograms. // TBI 20240124 is this safe enough?
  auto lBookParticleHistograms = cf_ph.cfBookParticleHistograms.value; // this is now the local version of that string array from configurable.
  if (lBookParticleHistograms.size() != eParticleHistograms_N) {
    LOGF(info, "\033[1;31m lBookParticleHistograms.size() = %d\033[0m", lBookParticleHistograms.size());
    LOGF(info, "\033[1;31m eParticleHistograms_N) = %d\033[0m", static_cast<int>(eParticleHistograms_N));
    LOGF(fatal, "in function \033[1;31m%s at line %d Mismatch in the number of flags in configurable cfBookParticleHistograms, and number of entries in enum eParticleHistograms \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering of particle histograms in the initialization in configurable cfBookParticleHistograms:
  // TBI 20240518 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  for (int name = 0; name < eParticleHistograms_N; name++) {
    // TBI 20240518 I could implement even a strickter EqualTo instead of EndsWith, but then I need to tokenize, etc., etc. This shall be safe enough.
    if (!TString(lBookParticleHistograms[name]).EndsWith(ph.fParticleHistogramsName[name].Data())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfBookParticleHistograms => name = %d, lBookParticleHistograms[name] = \"%s\", ph.fParticleHistogramsName[name] = \"%s\" \n Check if you are using an up to date tag. \033[0m", __FUNCTION__, __LINE__, name, TString(lBookParticleHistograms[name]).Data(), ph.fParticleHistogramsName[name].Data());
    }
  }

  // I append "&& ph.fFillParticleHistograms" below, to switch off booking of all 1D particle histograms with one common flag:
  ph.fBookParticleHistograms[ePhi] = Alright(lBookParticleHistograms[ePhi]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[ePt] = Alright(lBookParticleHistograms[ePt]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[eEta] = Alright(lBookParticleHistograms[eEta]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[eCharge] = Alright(lBookParticleHistograms[eCharge]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[etpcNClsFindable] = Alright(lBookParticleHistograms[etpcNClsFindable]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[etpcNClsShared] = Alright(lBookParticleHistograms[etpcNClsShared]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[eitsChi2NCl] = Alright(lBookParticleHistograms[eitsChi2NCl]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[etpcNClsFound] = Alright(lBookParticleHistograms[etpcNClsFound]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[etpcNClsCrossedRows] = Alright(lBookParticleHistograms[etpcNClsCrossedRows]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[eitsNCls] = Alright(lBookParticleHistograms[eitsNCls]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[eitsNClsInnerBarrel] = Alright(lBookParticleHistograms[eitsNClsInnerBarrel]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[etpcCrossedRowsOverFindableCls] = Alright(lBookParticleHistograms[etpcCrossedRowsOverFindableCls]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[etpcFoundOverFindableCls] = Alright(lBookParticleHistograms[etpcFoundOverFindableCls]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[etpcFractionSharedCls] = Alright(lBookParticleHistograms[etpcFractionSharedCls]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[etpcChi2NCl] = Alright(lBookParticleHistograms[etpcChi2NCl]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[edcaXY] = Alright(lBookParticleHistograms[edcaXY]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[edcaZ] = Alright(lBookParticleHistograms[edcaZ]) && ph.fFillParticleHistograms;
  ph.fBookParticleHistograms[ePDG] = Alright(lBookParticleHistograms[ePDG]) && ph.fFillParticleHistograms;
  // Remark #1: I do not need here anythig for etrackCutFlagFb1, etrackCutFlagFb2, ... eisGlobalTrack, because they are booleans
  // Remark #2: Nothing special here for ePtDependentDCAxyParameterization, because that is a string.

  // d) Particle histograms 2D:
  // By default all 2D particle histograms are booked. Set this flag to false to switch off booking of all 2D particle histograms:
  ph.fFillParticleHistograms2D = cf_ph.cfFillParticleHistograms2D;

  // If you do not want particular 2D particle histogram to be booked, use configurable array cfBookParticleHistograms2D, where you can specify flags 1 (book) or 0 (do not book).
  // *) Ordering of the flags in that array is interpreted through ordering of enums in enum eParticleHistograms2D.
  auto lBookParticleHistograms2D = cf_ph.cfBookParticleHistograms2D.value; // this is now the local version of that string array from configurable
  // TBI 20241113 For some reason, the default values of configurable "cfBookParticleHistograms2D" are not correctly propagated in the local variables, but I can circumvent that with JSON settings for the time being
  if (lBookParticleHistograms2D.size() != eParticleHistograms2D_N) {
    LOGF(info, "\033[1;31m lBookParticleHistograms2D.size() = %d\033[0m", lBookParticleHistograms2D.size());
    LOGF(info, "\033[1;31m eParticleHistograms2D_N) = %d\033[0m", static_cast<int>(eParticleHistograms2D_N));
    LOGF(fatal, "in function \033[1;31m%s at line %d Mismatch in the number of flags in configurable cfBookParticleHistograms2D, and number of entries in enum eParticleHistograms2D \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering of 2D particle histograms in the initialization in configurable cfBookParticleHistograms2D:
  // TBI 20241109 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  for (int name = 0; name < eParticleHistograms2D_N; name++) {
    // TBI 20241109 I could implement even a strickter EqualTo instead of EndsWith, but then I need to tokenize, etc., etc. This shall be safe enough.
    if (!TString(lBookParticleHistograms2D[name]).EndsWith(ph.fParticleHistogramsName2D[name].Data())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfBookParticleHistograms2D => name = %d, lBookParticleHistograms2D[name] = \"%s\", ph.fParticleHistogramsName2D[name] = \"%s\" \n Check if you are using an up to date tag. \033[0m", __FUNCTION__, __LINE__, name, TString(lBookParticleHistograms2D[name]).Data(), ph.fParticleHistogramsName2D[name].Data());
    }
  }

  // I append "&& ph.fFillParticleHistograms2D" below, to switch off booking of all 2D particle histograms with one common flag:
  ph.fBookParticleHistograms2D[ePhiPt] = Alright(lBookParticleHistograms2D[ePhiPt]) && ph.fFillParticleHistograms2D;
  ph.fBookParticleHistograms2D[ePhiEta] = Alright(lBookParticleHistograms2D[ePhiEta]) && ph.fFillParticleHistograms2D;

  // e) Particle sparse histograms:
  ph.fRebinSparse = cf_ph.cfRebinSparse;

  // *) Categories of sparse histograms:
  auto lBookParticleSparseHistograms = cf_ph.cfBookParticleSparseHistograms.value; // fill or not particulat category of sparse histograms
  if (lBookParticleSparseHistograms.size() != eDiffWeightCategory_N) {
    LOGF(info, "\033[1;31m lBookParticleSparseHistograms.size() = %d\033[0m", lBookParticleSparseHistograms.size());
    LOGF(info, "\033[1;31m eDiffWeightCategory_N) = %d\033[0m", static_cast<int>(eDiffWeightCategory_N));
    LOGF(fatal, "in function \033[1;31m%s at line %d Mismatch in the number of flags in configurable cfBookParticleSparseHistograms, and number of entries in enum eDiffWeightCategory_N \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering in the initialization in configurable cfBookParticleSparseHistograms:
  // TBI 20241109 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  // Algorithm: From [01]-DWPhi I tokenize with respect to "-" the 2nd field, and check if e.g. fParticleSparseHistogramsName_DWPhi ends with it.
  for (int name = 0; name < eDiffWeightCategory_N; name++) {
    TObjArray* oa = TString(lBookParticleSparseHistograms[name]).Tokenize("-");
    if (!oa) {
      LOGF(fatal, "\033[1;31m%s at line %d : name = %s\033[0m", __FUNCTION__, __LINE__, TString(lBookParticleSparseHistograms[name]).Data());
    }
    if (!ph.fParticleSparseHistogramsName[name].EndsWith(oa->At(1)->GetName())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfBookParticleSparseHistograms => name = %d, lBookParticleSparseHistograms[name] = \"%s\", ph.fParticleSparseHistogramsName[name] = \"%s\" \n Check if you are using an up to date tag. \033[0m", __FUNCTION__, __LINE__, name, TString(lBookParticleSparseHistograms[name]).Data(), ph.fParticleSparseHistogramsName[name].Data());
    }
    delete oa;
  }
  // Remark: below exceptionally I do not append the common flag with &&, since each of these flags already stands for one category
  ph.fBookParticleSparseHistograms[eDWPhi] = Alright(lBookParticleSparseHistograms[eDWPhi]);
  ph.fBookParticleSparseHistograms[eDWPt] = Alright(lBookParticleSparseHistograms[eDWPt]);
  ph.fBookParticleSparseHistograms[eDWEta] = Alright(lBookParticleSparseHistograms[eDWEta]);

  // f) QA:

  // **) QA 2D event histograms:
  qa.fFillQAEventHistograms2D = cf_qa.cfFillQAEventHistograms2D;

  // *) If you do not want particular 2D event histogram to be booked, use configurable array cfBookQAEventHistograms2D, where you can specify flags 1 (book) or 0 (do not book).
  // Ordering of the flags in that array is interpreted through ordering of enums in enum eQAEventHistograms2D
  auto lBookQAEventHistograms2D = cf_qa.cfBookQAEventHistograms2D.value; // this is now the local version of that string array from configurable
  // TBI 20241115 For some reason, the default values of configurable "cfBookQAEventHistograms2D" are not correctly propagated in the local variables, but I can circumvent that with JSON settings for the time being
  if (lBookQAEventHistograms2D.size() != eQAEventHistograms2D_N) {
    LOGF(info, "\033[1;31m lBookQAEventHistograms2D.size() = %d\033[0m", lBookQAEventHistograms2D.size());
    LOGF(info, "\033[1;31m eQAEventHistograms2D_N = %d\033[0m", static_cast<int>(eQAEventHistograms2D_N));
    LOGF(fatal, "in function \033[1;31m%s at line %d Mismatch in the number of flags in configurable cfBookQAEventHistograms2D, and number of entries in enum eQAEventHistograms2D \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering of QA 2D event histograms in the initialization in configurable cfBookQAEventHistograms2D:
  // TBI 20240518 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  for (int name = 0; name < eQAEventHistograms2D_N; name++) {
    // TBI 20240518 I could implement even a strickter EqualTo instead of EndsWith, but then I need to tokenize, etc., etc. This shall be safe enough.
    if (!TString(lBookQAEventHistograms2D[name]).EndsWith(qa.fEventHistogramsName2D[name].Data())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfBookQAEventHistograms2D => name = %d, lBookQAEventHistograms2D[name] = \"%s\", qa.fEventHistogramsName2D[name] = \"%s\" \n Check if you are using an up to date tag. \033[0m", __FUNCTION__, __LINE__, name, TString(lBookQAEventHistograms2D[name]).Data(), qa.fEventHistogramsName2D[name].Data());
    }
  }

  // I append "&& qa.fFillQAEventHistograms2D" below, to switch off booking of all 2D event histograms with one common flag:
  qa.fBookQAEventHistograms2D[eMultiplicity_vs_ReferenceMultiplicity] = Alright(lBookQAEventHistograms2D[eMultiplicity_vs_ReferenceMultiplicity]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eMultiplicity_vs_NContributors] = Alright(lBookQAEventHistograms2D[eMultiplicity_vs_NContributors]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eMultiplicity_vs_Centrality] = Alright(lBookQAEventHistograms2D[eMultiplicity_vs_Centrality]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eMultiplicity_vs_Vertex_z] = Alright(lBookQAEventHistograms2D[eMultiplicity_vs_Vertex_z]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eMultiplicity_vs_Occupancy] = Alright(lBookQAEventHistograms2D[eMultiplicity_vs_Occupancy]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eReferenceMultiplicity_vs_NContributors] = Alright(lBookQAEventHistograms2D[eReferenceMultiplicity_vs_NContributors]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eReferenceMultiplicity_vs_Centrality] = Alright(lBookQAEventHistograms2D[eReferenceMultiplicity_vs_Centrality]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eReferenceMultiplicity_vs_Vertex_z] = Alright(lBookQAEventHistograms2D[eReferenceMultiplicity_vs_Vertex_z]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eReferenceMultiplicity_vs_Occupancy] = Alright(lBookQAEventHistograms2D[eReferenceMultiplicity_vs_Occupancy]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eNContributors_vs_Centrality] = Alright(lBookQAEventHistograms2D[eNContributors_vs_Centrality]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eNContributors_vs_Vertex_z] = Alright(lBookQAEventHistograms2D[eNContributors_vs_Vertex_z]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eNContributors_vs_Occupancy] = Alright(lBookQAEventHistograms2D[eNContributors_vs_Occupancy]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentrality_vs_Vertex_z] = Alright(lBookQAEventHistograms2D[eCentrality_vs_Vertex_z]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentrality_vs_Occupancy] = Alright(lBookQAEventHistograms2D[eCentrality_vs_Occupancy]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentrality_vs_ImpactParameter] = Alright(lBookQAEventHistograms2D[eCentrality_vs_ImpactParameter]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eVertex_z_vs_Occupancy] = Alright(lBookQAEventHistograms2D[eVertex_z_vs_Occupancy]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eMultNTracksPV_vs_MultNTracksGlobal] = Alright(lBookQAEventHistograms2D[eMultNTracksPV_vs_MultNTracksGlobal]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentFT0C_vs_CentFT0CVariant1] = Alright(lBookQAEventHistograms2D[eCentFT0C_vs_CentFT0CVariant1]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentFT0C_vs_CentFT0M] = Alright(lBookQAEventHistograms2D[eCentFT0C_vs_CentFT0M]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentFT0C_vs_CentFV0A] = Alright(lBookQAEventHistograms2D[eCentFT0C_vs_CentFV0A]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentFT0C_vs_CentNTPV] = Alright(lBookQAEventHistograms2D[eCentFT0C_vs_CentNTPV]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentFT0C_vs_CentNGlobal] = Alright(lBookQAEventHistograms2D[eCentFT0C_vs_CentNGlobal]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentFT0M_vs_CentNTPV] = Alright(lBookQAEventHistograms2D[eCentFT0M_vs_CentNTPV]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCentRun2V0M_vs_CentRun2SPDTracklets] = Alright(lBookQAEventHistograms2D[eCentRun2V0M_vs_CentRun2SPDTracklets]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = Alright(lBookQAEventHistograms2D[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange]) && qa.fFillQAEventHistograms2D;
  qa.fBookQAEventHistograms2D[eCurrentRunDuration_vs_InteractionRate] = Alright(lBookQAEventHistograms2D[eCurrentRunDuration_vs_InteractionRate]) && qa.fFillQAEventHistograms2D;

  // **) QA 2D particle histograms:
  qa.fFillQAParticleHistograms2D = cf_qa.cfFillQAParticleHistograms2D;

  // *) If you do not want particular 2D particle histogram to be booked, use configurable array cfBookQAParticleHistograms2D, where you can specify flags 1 (book) or 0 (do not book).
  // Ordering of the flags in that array is interpreted through ordering of enums in enum eQAParticleHistograms2D.
  auto lBookQAParticleHistograms2D = cf_qa.cfBookQAParticleHistograms2D.value; // this is now the local version of that string array from configurable
  if (lBookQAParticleHistograms2D.size() != eQAParticleHistograms2D_N) {
    LOGF(info, "\033[1;31m lBookQAParticleHistograms2D.size() = %d\033[0m", lBookQAParticleHistograms2D.size());
    LOGF(info, "\033[1;31m eQAParticleHistograms2D_N = %d\033[0m", static_cast<int>(eQAParticleHistograms2D_N));
    LOGF(fatal, "in function \033[1;31m%s at line %d Mismatch in the number of flags in configurable cfBookQAParticleHistograms2D, and number of entries in enum eParticleHistograms2D \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering of QA 2D particle histograms in the initialization in configurable cfBookQAParticleHistograms2D:
  // TBI 20240518 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  for (int name = 0; name < eQAParticleHistograms2D_N; name++) {
    // TBI 20240518 I could implement even a strickter EqualTo instead of EndsWith, but then I need to tokenize, etc., etc. This shall be safe enough.
    if (!TString(lBookQAParticleHistograms2D[name]).EndsWith(qa.fParticleHistogramsName2D[name].Data())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfBookQAParticleHistograms2D => name = %d, lBookQAParticleHistograms2D[name] = \"%s\", qa.fParticleHistogramsName2D[name] = \"%s\" \n Check if you are using an up to date tag. \033[0m", __FUNCTION__, __LINE__, name, TString(lBookQAParticleHistograms2D[name]).Data(), qa.fParticleHistogramsName2D[name].Data());
    }
  }

  // I append "&& qa.fFillQAParticleHistograms2D" below, to switch off booking of all 2D particle histograms with one common flag:
  qa.fBookQAParticleHistograms2D[ePt_vs_dcaXY] = Alright(lBookQAParticleHistograms2D[ePt_vs_dcaXY]) && qa.fFillQAParticleHistograms2D;

  // **) QA 2D particle event histograms:
  qa.fFillQAParticleEventHistograms2D = cf_qa.cfFillQAParticleEventHistograms2D;

  // *) If you do not want particular 2D particle event histogram to be booked, use configurable array cfBookQAParticleEventHistograms2D, where you can specify flags 1 (book) or 0 (do not book).
  // Ordering of the flags in that array is interpreted through ordering of enums in enum eQAParticleEventHistograms2D.
  auto lBookQAParticleEventHistograms2D = cf_qa.cfBookQAParticleEventHistograms2D.value; // this is now the local version of that string array from configurable
  if (lBookQAParticleEventHistograms2D.size() != eQAParticleEventHistograms2D_N) {
    LOGF(info, "\033[1;31m lBookQAParticleEventHistograms2D.size() = %d\033[0m", lBookQAParticleEventHistograms2D.size());
    LOGF(info, "\033[1;31m eQAParticleEventHistograms2D_N = %d\033[0m", static_cast<int>(eQAParticleEventHistograms2D_N));
    LOGF(fatal, "in function \033[1;31m%s at line %d Mismatch in the number of flags in configurable cfBookQAParticleEventHistograms2D, and number of entries in enum eParticleEventHistograms2D \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering of QA 2D particle event histograms in the initialization in configurable cfBookQAParticleEventHistograms2D:
  // TBI 20240518 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  for (int name = 0; name < eQAParticleEventHistograms2D_N; name++) {
    // TBI 20240518 I could implement even a strickter EqualTo instead of EndsWith, but then I need to tokenize, etc., etc. This shall be safe enough.
    if (!TString(lBookQAParticleEventHistograms2D[name]).EndsWith(qa.fQAParticleEventHistogramsName2D[name].Data())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfBookQAParticleEventHistograms2D => name = %d, lBookQAParticleEventHistograms2D[name] = \"%s\", qa.fParticleEventHistogramsName2D[name] = \"%s\" \n Check if you are using an up to date tag. \033[0m", __FUNCTION__, __LINE__, name, TString(lBookQAParticleEventHistograms2D[name]).Data(), qa.fQAParticleEventHistogramsName2D[name].Data());
    }
  }

  // I append "&& qa.fFillQAParticleEventHistograms2D" below, to switch off booking of all 2D particle event histograms with one common flag:
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsEbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsEbyE]) && qa.fFillQAParticleEventHistograms2D;
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsNegEtaEbyE]) && qa.fFillQAParticleEventHistograms2D;
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsPosEtaEbyE]) && qa.fFillQAParticleEventHistograms2D;
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0804EbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0804EbyE]) && qa.fFillQAParticleEventHistograms2D;
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0400EbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0400EbyE]) && qa.fFillQAParticleEventHistograms2D;
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0004EbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0004EbyE]) && qa.fFillQAParticleEventHistograms2D;
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0408EbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0408EbyE]) && qa.fFillQAParticleEventHistograms2D;
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt0005EbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt0005EbyE]) && qa.fFillQAParticleEventHistograms2D;
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt0510EbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt0510EbyE]) && qa.fFillQAParticleEventHistograms2D;
  qa.fBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt1050EbyE] = Alright(lBookQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt1050EbyE]) && qa.fFillQAParticleEventHistograms2D;

  // **) QA 2D "correlations vs." histograms:
  qa.fFillQACorrelationsVsHistograms2D = cf_qa.cfFillQACorrelationsVsHistograms2D;

  // *) If you do not want particular 2D "correlations vs." histogram to be booked, use configurable array cfBookQACorrelationsVsHistograms2D, where you can specify flags 1 (book) or 0 (do not book).
  // Ordering of the flags in that array is interpreted through ordering of enums in enum eQACorrelationsVsHistograms2D.
  auto lBookQACorrelationsVsHistograms2D = cf_qa.cfBookQACorrelationsVsHistograms2D.value; // this is now the local version of that string array from configurable
  if (lBookQACorrelationsVsHistograms2D.size() != eQACorrelationsVsHistograms2D_N) {
    LOGF(info, "\033[1;31m lBookQACorrelationsVsHistograms2D.size() = %d\033[0m", lBookQACorrelationsVsHistograms2D.size());
    LOGF(info, "\033[1;31m eQACorrelationsVsHistograms2D_N = %d\033[0m", static_cast<int>(eQACorrelationsVsHistograms2D_N));
    LOGF(fatal, "in function \033[1;31m%s at line %d Mismatch in the number of flags in configurable cfBookQACorrelationsVsHistograms2D, and number of entries in enum eCorrelationsVsHistograms2D \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering of QA 2D "correlations vs." histograms in the initialization in configurable cfBookQACorrelationsVsHistograms2D:
  // TBI 20240518 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  for (int name = 0; name < eQACorrelationsVsHistograms2D_N; name++) {
    // TBI 20240518 I could implement even a strickter EqualTo instead of EndsWith, but then I need to tokenize, etc., etc. This shall be safe enough.
    if (!TString(lBookQACorrelationsVsHistograms2D[name]).EndsWith(qa.fQACorrelationsVsHistogramsName2D[name].Data())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfBookQACorrelationsVsHistograms2D => name = %d, lBookQACorrelationsVsHistograms2D[name] = \"%s\", qa.fCorrelationsVsHistogramsName2D[name] = \"%s\" \n Check if you are using an up to date tag. \033[0m", __FUNCTION__, __LINE__, name, TString(lBookQACorrelationsVsHistograms2D[name]).Data(), qa.fQACorrelationsVsHistogramsName2D[name].Data());
    }
  }

  // I append "&& qa.fFillQACorrelationsVsHistograms2D" below, to switch off booking of all 2D "correlations vs." histograms with one common flag:
  qa.fBookQACorrelationsVsHistograms2D[eCorrelations_vs_Multiplicity] = Alright(lBookQACorrelationsVsHistograms2D[eCorrelations_vs_Multiplicity]) && qa.fFillQACorrelationsVsHistograms2D;
  qa.fBookQACorrelationsVsHistograms2D[eCorrelations_vs_ReferenceMultiplicity] = Alright(lBookQACorrelationsVsHistograms2D[eCorrelations_vs_ReferenceMultiplicity]) && qa.fFillQACorrelationsVsHistograms2D;
  qa.fBookQACorrelationsVsHistograms2D[eCorrelations_vs_Centrality] = Alright(lBookQACorrelationsVsHistograms2D[eCorrelations_vs_Centrality]) && qa.fFillQACorrelationsVsHistograms2D;
  // .....
  qa.fBookQACorrelationsVsHistograms2D[eCorrelations_vs_MeanPhi] = Alright(lBookQACorrelationsVsHistograms2D[eCorrelations_vs_MeanPhi]) && qa.fFillQACorrelationsVsHistograms2D;
  qa.fBookQACorrelationsVsHistograms2D[eCorrelations_vs_MeanPt] = Alright(lBookQACorrelationsVsHistograms2D[eCorrelations_vs_MeanPt]) && qa.fFillQACorrelationsVsHistograms2D;
  qa.fBookQACorrelationsVsHistograms2D[eCorrelations_vs_MeanEta] = Alright(lBookQACorrelationsVsHistograms2D[eCorrelations_vs_MeanEta]) && qa.fFillQACorrelationsVsHistograms2D;
  // .....

  // *) min and max harmonics for which this series of histograms will be booked:
  auto lQACorrelationsVsHistogramsMinMaxHarmonic = cf_qa.cfQACorrelationsVsHistogramsMinMaxHarmonic.value;
  qa.fQACorrelationsVsHistogramsMinMaxHarmonic[eMin] = lQACorrelationsVsHistogramsMinMaxHarmonic[eMin];
  qa.fQACorrelationsVsHistogramsMinMaxHarmonic[eMax] = lQACorrelationsVsHistogramsMinMaxHarmonic[eMax];
  // **) insanity check:
  if (!(qa.fQACorrelationsVsHistogramsMinMaxHarmonic[eMin] < qa.fQACorrelationsVsHistogramsMinMaxHarmonic[eMax])) {
    LOGF(fatal, "\033[1;31m%s at line %d : wrong setting for min and max harmonics: min = %d, max = %d \033[0m", __FUNCTION__, __LINE__, qa.fQACorrelationsVsHistogramsMinMaxHarmonic[eMin], qa.fQACorrelationsVsHistogramsMinMaxHarmonic[eMax]);
  }

  // ...

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void DefaultBooking()

//============================================================

void DefaultBinning()
{
  // Default binning for all histograms.

  // TBI 20240114 If some of these values are going to change frequently, add support for them in MuPa-Configurables.h,
  // in the same way I did it for DefaultCuts().
  // At the moment, I added to configurables support only for binning of sparse histograms, because there memory managment is critical.

  // a) Default binning for event histograms;
  // b) Default binning for particle histograms 1D;
  // c) Default binning for particle histograms 2D;
  // d) Default binning for results histograms;
  // e) Variable-length binning for results histograms set via MuPa-Configurables.h.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Default binning for event histograms:
  eh.fEventHistogramsBins[eNumberOfEvents][0] = 1;
  eh.fEventHistogramsBins[eNumberOfEvents][1] = 0.;
  eh.fEventHistogramsBins[eNumberOfEvents][2] = 1.;

  eh.fEventHistogramsBins[eTotalMultiplicity][0] = 10000.;
  eh.fEventHistogramsBins[eTotalMultiplicity][1] = 0.;
  eh.fEventHistogramsBins[eTotalMultiplicity][2] = 100000.;

  eh.fEventHistogramsBins[eMultiplicity][0] = 2000.;
  eh.fEventHistogramsBins[eMultiplicity][1] = 0.;
  eh.fEventHistogramsBins[eMultiplicity][2] = 20000.;

  eh.fEventHistogramsBins[eReferenceMultiplicity][0] = 6000.;
  eh.fEventHistogramsBins[eReferenceMultiplicity][1] = 0.;
  eh.fEventHistogramsBins[eReferenceMultiplicity][2] = 60000.;

  eh.fEventHistogramsBins[eCentrality][0] = 110; // intentionally, because if centrality is not determined, it's set to 105.0 at the moment
  eh.fEventHistogramsBins[eCentrality][1] = 0.;
  eh.fEventHistogramsBins[eCentrality][2] = 110.;

  eh.fEventHistogramsBins[eVertex_x][0] = 1600;
  eh.fEventHistogramsBins[eVertex_x][1] = -0.8;
  eh.fEventHistogramsBins[eVertex_x][2] = 0.8;

  eh.fEventHistogramsBins[eVertex_y][0] = 1600;
  eh.fEventHistogramsBins[eVertex_y][1] = -0.8;
  eh.fEventHistogramsBins[eVertex_y][2] = 0.8;

  eh.fEventHistogramsBins[eVertex_z][0] = 800;
  eh.fEventHistogramsBins[eVertex_z][1] = -40.;
  eh.fEventHistogramsBins[eVertex_z][2] = 40.;

  eh.fEventHistogramsBins[eNContributors][0] = 1000.;
  eh.fEventHistogramsBins[eNContributors][1] = 0.;
  eh.fEventHistogramsBins[eNContributors][2] = 10000.;

  eh.fEventHistogramsBins[eImpactParameter][0] = 1000;
  eh.fEventHistogramsBins[eImpactParameter][1] = 0.;
  eh.fEventHistogramsBins[eImpactParameter][2] = 100.;

  eh.fEventHistogramsBins[eEventPlaneAngle][0] = 720;
  eh.fEventHistogramsBins[eEventPlaneAngle][1] = -o2::constants::math::PI; // just in case somebody uses the convention -Pi < EP < Pi, instead of 0 < EP < 2Pi
  eh.fEventHistogramsBins[eEventPlaneAngle][2] = o2::constants::math::TwoPI;

  eh.fEventHistogramsBins[eOccupancy][0] = 1000;
  eh.fEventHistogramsBins[eOccupancy][1] = 0.;
  eh.fEventHistogramsBins[eOccupancy][2] = 100000.;

  eh.fEventHistogramsBins[eInteractionRate][0] = 1000;
  eh.fEventHistogramsBins[eInteractionRate][1] = 0.;
  eh.fEventHistogramsBins[eInteractionRate][2] = 100.;

  eh.fEventHistogramsBins[eCurrentRunDuration][0] = 10000;
  eh.fEventHistogramsBins[eCurrentRunDuration][1] = 0.;
  eh.fEventHistogramsBins[eCurrentRunDuration][2] = 10000.;

  // b) Default binning for particle histograms 1D:
  ph.fParticleHistogramsBins[ePhi][0] = 360;
  ph.fParticleHistogramsBins[ePhi][1] = 0.;
  ph.fParticleHistogramsBins[ePhi][2] = o2::constants::math::TwoPI;

  ph.fParticleHistogramsBins[ePt][0] = 2000;
  ph.fParticleHistogramsBins[ePt][1] = 0.;
  ph.fParticleHistogramsBins[ePt][2] = 200.;

  ph.fParticleHistogramsBins[eEta][0] = 500;
  ph.fParticleHistogramsBins[eEta][1] = -5.;
  ph.fParticleHistogramsBins[eEta][2] = 5.;

  ph.fParticleHistogramsBins[eCharge][0] = 7;
  ph.fParticleHistogramsBins[eCharge][1] = -3.5; // anticipating I might be storing charge of Delta++, etc.
  ph.fParticleHistogramsBins[eCharge][2] = 3.5;

  ph.fParticleHistogramsBins[etpcNClsFindable][0] = 300;
  ph.fParticleHistogramsBins[etpcNClsFindable][1] = 0.;
  ph.fParticleHistogramsBins[etpcNClsFindable][2] = 300.;

  ph.fParticleHistogramsBins[etpcNClsShared][0] = 200;
  ph.fParticleHistogramsBins[etpcNClsShared][1] = 0.;
  ph.fParticleHistogramsBins[etpcNClsShared][2] = 200.;

  ph.fParticleHistogramsBins[eitsChi2NCl][0] = 200;
  ph.fParticleHistogramsBins[eitsChi2NCl][1] = 0.;
  ph.fParticleHistogramsBins[eitsChi2NCl][2] = 200.;

  ph.fParticleHistogramsBins[etpcNClsFound][0] = 200;
  ph.fParticleHistogramsBins[etpcNClsFound][1] = 0.;
  ph.fParticleHistogramsBins[etpcNClsFound][2] = 200.;

  ph.fParticleHistogramsBins[etpcNClsCrossedRows][0] = 200;
  ph.fParticleHistogramsBins[etpcNClsCrossedRows][1] = 0.;
  ph.fParticleHistogramsBins[etpcNClsCrossedRows][2] = 200.;

  ph.fParticleHistogramsBins[eitsNCls][0] = 10;
  ph.fParticleHistogramsBins[eitsNCls][1] = 0.;
  ph.fParticleHistogramsBins[eitsNCls][2] = 10.;

  ph.fParticleHistogramsBins[eitsNClsInnerBarrel][0] = 10;
  ph.fParticleHistogramsBins[eitsNClsInnerBarrel][1] = 0.;
  ph.fParticleHistogramsBins[eitsNClsInnerBarrel][2] = 10.;

  ph.fParticleHistogramsBins[etpcCrossedRowsOverFindableCls][0] = 1000;
  ph.fParticleHistogramsBins[etpcCrossedRowsOverFindableCls][1] = 0.;
  ph.fParticleHistogramsBins[etpcCrossedRowsOverFindableCls][2] = 10;

  ph.fParticleHistogramsBins[etpcFoundOverFindableCls][0] = 1000;
  ph.fParticleHistogramsBins[etpcFoundOverFindableCls][1] = 0.;
  ph.fParticleHistogramsBins[etpcFoundOverFindableCls][2] = 10.;

  ph.fParticleHistogramsBins[etpcFractionSharedCls][0] = 110;
  ph.fParticleHistogramsBins[etpcFractionSharedCls][1] = -1.; // yes, I saw here entries with negative values TBI 20240507 check what are these values
  ph.fParticleHistogramsBins[etpcFractionSharedCls][2] = 10.;

  ph.fParticleHistogramsBins[etpcChi2NCl][0] = 2500;
  ph.fParticleHistogramsBins[etpcChi2NCl][1] = 0.;
  ph.fParticleHistogramsBins[etpcChi2NCl][2] = 250.;

  ph.fParticleHistogramsBins[edcaXY][0] = 2000;
  ph.fParticleHistogramsBins[edcaXY][1] = -10.;
  ph.fParticleHistogramsBins[edcaXY][2] = 10.;

  ph.fParticleHistogramsBins[edcaZ][0] = 2000;
  ph.fParticleHistogramsBins[edcaZ][1] = -10.;
  ph.fParticleHistogramsBins[edcaZ][2] = 10.;

  ph.fParticleHistogramsBins[ePDG][0] = 2000;
  ph.fParticleHistogramsBins[ePDG][1] = -1000.;
  ph.fParticleHistogramsBins[ePDG][2] = 1000.;

  // c) Default binning for particle histograms 2D:
  //    At the moment, for fixed binning, I just re-use the binning of corresponding 1D histograms.
  //    For variable-length binning, I use binning from fResultsPro[], as for other histograms.
  ph.fParticleHistogramsBins2D[ePhiPt][eX][0] = ph.fParticleHistogramsBins[ePhi][0];
  ph.fParticleHistogramsBins2D[ePhiPt][eX][1] = ph.fParticleHistogramsBins[ePhi][1];
  ph.fParticleHistogramsBins2D[ePhiPt][eX][2] = ph.fParticleHistogramsBins[ePhi][2];
  ph.fParticleHistogramsBins2D[ePhiPt][eY][0] = ph.fParticleHistogramsBins[ePt][0];
  ph.fParticleHistogramsBins2D[ePhiPt][eY][1] = ph.fParticleHistogramsBins[ePt][1];
  ph.fParticleHistogramsBins2D[ePhiPt][eY][2] = ph.fParticleHistogramsBins[ePt][2];

  ph.fParticleHistogramsBins2D[ePhiEta][eX][0] = ph.fParticleHistogramsBins[ePhi][0];
  ph.fParticleHistogramsBins2D[ePhiEta][eX][1] = ph.fParticleHistogramsBins[ePhi][1];
  ph.fParticleHistogramsBins2D[ePhiEta][eX][2] = ph.fParticleHistogramsBins[ePhi][2];
  ph.fParticleHistogramsBins2D[ePhiEta][eY][0] = ph.fParticleHistogramsBins[eEta][0];
  ph.fParticleHistogramsBins2D[ePhiEta][eY][1] = ph.fParticleHistogramsBins[eEta][1];
  ph.fParticleHistogramsBins2D[ePhiEta][eY][2] = ph.fParticleHistogramsBins[eEta][2];

  // d) Default binning for results histograms:
  //    Remark: These bins apply to following categories fCorrelationsPro, fNestedLoopsPro, fTest0Pro, and fResultsPro.
  // *) For integrated resullts, binning is always the same:
  res.fResultsProFixedLengthBins[AFO_INTEGRATED][0] = 1;
  res.fResultsProFixedLengthBins[AFO_INTEGRATED][1] = 0.;
  res.fResultsProFixedLengthBins[AFO_INTEGRATED][2] = 1.;
  // *) Fixed-length binning vs. multiplicity:
  this->InitializeFixedLengthBins(AFO_MULTIPLICITY);
  // *) Fixed-length binning vs. centrality:
  this->InitializeFixedLengthBins(AFO_CENTRALITY);
  // *) Fixed-length binning vs. pt:
  this->InitializeFixedLengthBins(AFO_PT);
  // *) Fixed-length binning vs. eta:
  this->InitializeFixedLengthBins(AFO_ETA);
  // *) Fixed-length binning vs. occupancy:
  this->InitializeFixedLengthBins(AFO_OCCUPANCY);
  // *) Fixed-length binning vs. interaction rate:
  this->InitializeFixedLengthBins(AFO_INTERACTIONRATE);
  // *) Fixed-length binning vs. run duration:
  this->InitializeFixedLengthBins(AFO_CURRENTRUNDURATION);
  // *) Vertex z position:
  this->InitializeFixedLengthBins(AFO_VZ);

  // e) Variable-length binning set via MuPa-Configurables.h:
  // *) Variable-length binning vs. multiplicity:
  if (cf_res.cfUseVariableLength_mult_bins) {
    this->InitializeVariableLengthBins(AFO_MULTIPLICITY);
  }
  // *) Variable-length binning vs. centrality:
  if (cf_res.cfUseVariableLength_cent_bins) {
    this->InitializeVariableLengthBins(AFO_CENTRALITY);
  }
  // *) Variable-length binning vs. pt:
  if (cf_res.cfUseVariableLength_pt_bins) {
    this->InitializeVariableLengthBins(AFO_PT);
  }
  // *) Variable-length binning vs. eta:
  if (cf_res.cfUseVariableLength_eta_bins) {
    this->InitializeVariableLengthBins(AFO_ETA);
  }
  // *) Variable-length binning vs. occupancy:
  if (cf_res.cfUseVariableLength_occu_bins) {
    this->InitializeVariableLengthBins(AFO_OCCUPANCY);
  }
  // *) Variable-length binning vs. interaction rate:
  if (cf_res.cfUseVariableLength_ir_bins) {
    this->InitializeVariableLengthBins(AFO_INTERACTIONRATE);
  }
  // *) Variable-length binning vs. run duration:
  if (cf_res.cfUseVariableLength_crd_bins) {
    this->InitializeVariableLengthBins(AFO_CURRENTRUNDURATION);
  }
  // *) Variable-length binning vs. vertex z position:
  if (cf_res.cfUseVariableLength_vz_bins) {
    this->InitializeVariableLengthBins(AFO_VZ);
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void DefaultBinning()

//============================================================

void InitializeFixedLengthBins(eAsFunctionOf AFO)
{
  // This is a helper function to suppress code bloat in DefaultBinning().
  // It merely initalizes res.fResultsProFixedLengthBins[...] from corresponding configurables + a few other minor thingies.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Common local vector for all fixed-length bins:
  vector<float> lFixedLength_bins;

  switch (AFO) {
    case AFO_MULTIPLICITY:
      lFixedLength_bins = cf_res.cfFixedLength_mult_bins.value;
      break;
    case AFO_CENTRALITY:
      lFixedLength_bins = cf_res.cfFixedLength_cent_bins.value;
      break;
    case AFO_PT:
      lFixedLength_bins = cf_res.cfFixedLength_pt_bins.value;
      break;
    case AFO_ETA:
      lFixedLength_bins = cf_res.cfFixedLength_eta_bins.value;
      break;
    case AFO_OCCUPANCY:
      lFixedLength_bins = cf_res.cfFixedLength_occu_bins.value;
      break;
    case AFO_INTERACTIONRATE:
      lFixedLength_bins = cf_res.cfFixedLength_ir_bins.value;
      break;
    case AFO_CURRENTRUNDURATION:
      lFixedLength_bins = cf_res.cfFixedLength_crd_bins.value;
      break;
    case AFO_VZ:
      lFixedLength_bins = cf_res.cfFixedLength_vz_bins.value;
      break;
    // ...
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : This enum AFO = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(AFO));
      break;
  } // switch(AFO)

  // From this point onward, the code is the same for any AFO variable:
  if (lFixedLength_bins.size() != 3) {
    LOGF(fatal, "in function \033[1;31m%s at line %d => The array cfFixedLength_bins must have have 3 entries: {nBins, min, max} \n \033[0m", __FUNCTION__, __LINE__);
  }
  res.fResultsProFixedLengthBins[AFO][0] = lFixedLength_bins[0];
  res.fResultsProFixedLengthBins[AFO][1] = lFixedLength_bins[1];
  res.fResultsProFixedLengthBins[AFO][2] = lFixedLength_bins[2];

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m %s : fixed-length %s bins \033[0m", __FUNCTION__, res.fResultsProXaxisTitle[AFO].Data());
    LOGF(info, "\033[1;32m [0] : %f \033[0m", res.fResultsProFixedLengthBins[AFO][0]);
    LOGF(info, "\033[1;32m [1] : %f \033[0m", res.fResultsProFixedLengthBins[AFO][1]);
    LOGF(info, "\033[1;32m [2] : %f \033[0m", res.fResultsProFixedLengthBins[AFO][2]);
    ExitFunction(__FUNCTION__);
  }

} // void InitializeFixedLengthBins(eAsFunctionOf AFO)

//============================================================

void InitializeVariableLengthBins(eAsFunctionOf AFO)
{
  // This is a helper function to suppress code bloat in DefaultBinning().
  // It merely initalizes res.fResultsProVariableLengthBins[...] from corresponding configurables + a few other minor thingies.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Common local vector for all variable-length bins:
  vector<float> lVariableLength_bins;

  switch (AFO) {
    case AFO_MULTIPLICITY:
      lVariableLength_bins = cf_res.cfVariableLength_mult_bins.value;
      break;
    case AFO_CENTRALITY:
      lVariableLength_bins = cf_res.cfVariableLength_cent_bins.value;
      break;
    case AFO_PT:
      lVariableLength_bins = cf_res.cfVariableLength_pt_bins.value;
      break;
    case AFO_ETA:
      lVariableLength_bins = cf_res.cfVariableLength_eta_bins.value;
      break;
    case AFO_OCCUPANCY:
      lVariableLength_bins = cf_res.cfVariableLength_occu_bins.value;
      break;
    case AFO_INTERACTIONRATE:
      lVariableLength_bins = cf_res.cfVariableLength_ir_bins.value;
      break;
    case AFO_CURRENTRUNDURATION:
      lVariableLength_bins = cf_res.cfVariableLength_crd_bins.value;
      break;
    case AFO_VZ:
      lVariableLength_bins = cf_res.cfVariableLength_vz_bins.value;
      break;
    // ...
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : This enum AFO = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(AFO));
      break;
  } // switch(AFO)

  // From this point onward, the code is the same for any AFO variable:
  res.fUseResultsProVariableLengthBins[AFO] = true;
  if (lVariableLength_bins.size() < 2) {
    LOGF(fatal, "in function \033[1;31m%s at line %d => The array cfVariableLength_bins must have at least 2 entries \n \033[0m", __FUNCTION__, __LINE__);
  }
  res.fResultsProVariableLengthBins[AFO] = new TArrayF(lVariableLength_bins.size(), lVariableLength_bins.data());
  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m %s : variable-length %s bins \033[0m", __FUNCTION__, res.fResultsProXaxisTitle[AFO].Data());
    for (int i = 0; i < res.fResultsProVariableLengthBins[AFO]->GetSize(); i++) {
      LOGF(info, "\033[1;32m [%d] : %f \033[0m", i, res.fResultsProVariableLengthBins[AFO]->GetAt(i));
    }
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void InitializeVariableLengthBins(eAsFunctionOf AFO)

//============================================================

void CastStringIntoArray(int AFO)
{
  // Temporary function, to be removed eventually. Here temporarily I am casting e.g. a string "1.0,2.0,5.0" into corresponding TArrayD.

  // TBI 20241019 This function is not needed any longer, remove eventually.

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__);
  }

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m Casting a string %s into TArrayD .... \033[0m", res.fResultsProVariableLengthBinsString[AFO].Data());
  }

  TObjArray* oa = res.fResultsProVariableLengthBinsString[AFO].Tokenize(",");
  if (!oa) {
    LOGF(fatal, "in function \033[1;31m%s at line %d \n fResultsProVariableLengthBinsString[AFO] = %s\033[0m", __FUNCTION__, __LINE__, res.fResultsProVariableLengthBinsString[AFO].Data());
  }
  int nEntries = oa->GetEntries();
  res.fResultsProVariableLengthBins[AFO] = new TArrayF(nEntries);
  for (int i = 0; i < nEntries; i++) {
    res.fResultsProVariableLengthBins[AFO]->AddAt(TString(oa->At(i)->GetName()).Atof(), i);
  }
  delete oa; // yes, otherwise it's a memory leak

  if (tc.fVerbose) {
    for (int i = 0; i < res.fResultsProVariableLengthBins[AFO]->GetSize(); i++) {
      LOGF(info, "\033[1;32m [%d] : %f \033[0m", i, res.fResultsProVariableLengthBins[AFO]->At(i));
    }
  }

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m Done! \033[0m");
  }

} // void CastStringIntoArray(int AFO)

//============================================================

void DefaultCuts()
{
  // Define default cuts. Default cuts are hardwired in MuPa-Configurables.h.

  // a) Default event cuts;
  // b) Default particle cuts.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Default event cuts:

  // *) Use or do not use a cut enumerated in eEventHistograms + eEventCuts.
  //    Default cuts are set in configurable cfUseEventCuts
  auto lUseEventCuts = (vector<string>)cf_ec.cfUseEventCuts;
  if (lUseEventCuts.size() != eEventCuts_N) {
    LOGF(info, "\033[1;31m lUseEventCuts.size() = %d\033[0m", lUseEventCuts.size());
    LOGF(info, "\033[1;31m eEventCuts_N = %d\033[0m", static_cast<int>(eEventCuts_N));
    LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfUseEventCuts, and number of entries in enum eEventHistograms + eEventCuts \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering of event cuts in the initialization in configurable cfUseEventCuts:
  // TBI 20240518 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  for (int name = 0; name < eEventCuts_N; name++) {
    // TBI 20240518 I could implement even a strickter EqualTo instead of EndsWith, but then I need to tokenize, etc., etc. This shall be safe enough.
    if (!TString(lUseEventCuts[name]).EndsWith(ec.fEventCutName[name].Data())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfUseEventCuts => name = %d, lUseEventCuts[name] = \"%s\", ec.fEventCutName[name] = \"%s\" \033[0m", __FUNCTION__, __LINE__, name, TString(lUseEventCuts[name]).Data(), ec.fEventCutName[name].Data());
    }
  }

  // *) from enum eEventHistograms:
  ec.fUseEventCuts[eNumberOfEvents] = Alright(lUseEventCuts[eNumberOfEvents]); // total number of events (before event cuts)
  ec.fUseEventCuts[eTotalMultiplicity] = Alright(lUseEventCuts[eTotalMultiplicity]);
  ec.fUseEventCuts[eMultiplicity] = Alright(lUseEventCuts[eMultiplicity]);
  ec.fUseEventCuts[eReferenceMultiplicity] = Alright(lUseEventCuts[eReferenceMultiplicity]);
  ec.fUseEventCuts[eCentrality] = Alright(lUseEventCuts[eCentrality]);
  ec.fUseEventCuts[eVertex_x] = Alright(lUseEventCuts[eVertex_x]);
  ec.fUseEventCuts[eVertex_y] = Alright(lUseEventCuts[eVertex_y]);
  ec.fUseEventCuts[eVertex_z] = Alright(lUseEventCuts[eVertex_z]);
  ec.fUseEventCuts[eNContributors] = Alright(lUseEventCuts[eNContributors]);
  ec.fUseEventCuts[eImpactParameter] = Alright(lUseEventCuts[eImpactParameter]);
  ec.fUseEventCuts[eEventPlaneAngle] = Alright(lUseEventCuts[eEventPlaneAngle]);
  ec.fUseEventCuts[eOccupancy] = Alright(lUseEventCuts[eOccupancy]);
  ec.fUseEventCuts[eInteractionRate] = Alright(lUseEventCuts[eInteractionRate]);
  ec.fUseEventCuts[eCurrentRunDuration] = Alright(lUseEventCuts[eCurrentRunDuration]);
  ec.fUseEventCuts[eMultMCNParticlesEta08] = Alright(lUseEventCuts[eMultMCNParticlesEta08]);

  // *) from enum eEventCuts:
  ec.fUseEventCuts[eTrigger] = Alright(lUseEventCuts[eTrigger]);
  ec.fUseEventCuts[eSel7] = Alright(lUseEventCuts[eSel7]);
  ec.fUseEventCuts[eSel8] = Alright(lUseEventCuts[eSel8]);
  ec.fUseEventCuts[eMultiplicityEstimator] = Alright(lUseEventCuts[eMultiplicityEstimator]);
  ec.fUseEventCuts[eCentralityEstimator] = Alright(lUseEventCuts[eCentralityEstimator]);
  ec.fUseEventCuts[eReferenceMultiplicityEstimator] = Alright(lUseEventCuts[eReferenceMultiplicityEstimator]);
  ec.fUseEventCuts[eSelectedEvents] = Alright(lUseEventCuts[eSelectedEvents]); // selected number of events (after all event cuts)
  ec.fUseEventCuts[eNoSameBunchPileup] = Alright(lUseEventCuts[eNoSameBunchPileup]);
  ec.fUseEventCuts[eIsGoodZvtxFT0vsPV] = Alright(lUseEventCuts[eIsGoodZvtxFT0vsPV]);
  ec.fUseEventCuts[eIsVertexITSTPC] = Alright(lUseEventCuts[eIsVertexITSTPC]);
  ec.fUseEventCuts[eIsVertexTOFmatched] = Alright(lUseEventCuts[eIsVertexTOFmatched]);
  ec.fUseEventCuts[eIsVertexTRDmatched] = Alright(lUseEventCuts[eIsVertexTRDmatched]);
  ec.fUseEventCuts[eNoCollInTimeRangeStrict] = Alright(lUseEventCuts[eNoCollInTimeRangeStrict]);
  ec.fUseEventCuts[eNoCollInTimeRangeStandard] = Alright(lUseEventCuts[eNoCollInTimeRangeStandard]);
  ec.fUseEventCuts[eNoCollInRofStrict] = Alright(lUseEventCuts[eNoCollInRofStrict]);
  ec.fUseEventCuts[eNoCollInRofStandard] = Alright(lUseEventCuts[eNoCollInRofStandard]);
  ec.fUseEventCuts[eNoHighMultCollInPrevRof] = Alright(lUseEventCuts[eNoHighMultCollInPrevRof]);
  ec.fUseEventCuts[eIsGoodITSLayer3] = Alright(lUseEventCuts[eIsGoodITSLayer3]);
  ec.fUseEventCuts[eIsGoodITSLayer0123] = Alright(lUseEventCuts[eIsGoodITSLayer0123]);
  ec.fUseEventCuts[eIsGoodITSLayersAll] = Alright(lUseEventCuts[eIsGoodITSLayersAll]);
  ec.fUseEventCuts[eOccupancyEstimator] = Alright(lUseEventCuts[eOccupancyEstimator]);
  ec.fUseEventCuts[eMinVertexDistanceFromIP] = Alright(lUseEventCuts[eMinVertexDistanceFromIP]);
  ec.fUseEventCuts[eCentralityWeights] = Alright(lUseEventCuts[eCentralityWeights]);

  // **) event cuts defined via booleans:
  ec.fUseEventCuts[eSel7] = ec.fUseEventCuts[eSel7] && cf_ec.cfUseSel7;
  ec.fUseEventCuts[eSel8] = ec.fUseEventCuts[eSel8] && cf_ec.cfUseSel8;
  ec.fUseEventCuts[eNoSameBunchPileup] = ec.fUseEventCuts[eNoSameBunchPileup] && cf_ec.cfUseNoSameBunchPileup;
  ec.fUseEventCuts[eIsGoodZvtxFT0vsPV] = ec.fUseEventCuts[eIsGoodZvtxFT0vsPV] && cf_ec.cfUseIsGoodZvtxFT0vsPV;
  ec.fUseEventCuts[eIsVertexITSTPC] = ec.fUseEventCuts[eIsVertexITSTPC] && cf_ec.cfUseIsVertexITSTPC;
  ec.fUseEventCuts[eIsVertexTOFmatched] = ec.fUseEventCuts[eIsVertexTOFmatched] && cf_ec.cfUseIsVertexTOFmatched;
  ec.fUseEventCuts[eIsVertexTRDmatched] = ec.fUseEventCuts[eIsVertexTRDmatched] && cf_ec.cfUseIsVertexTRDmatched;
  ec.fUseEventCuts[eNoCollInTimeRangeStrict] = ec.fUseEventCuts[eNoCollInTimeRangeStrict] && cf_ec.cfUseNoCollInTimeRangeStrict;
  ec.fUseEventCuts[eNoCollInTimeRangeStandard] = ec.fUseEventCuts[eNoCollInTimeRangeStandard] && cf_ec.cfUseNoCollInTimeRangeStandard;
  ec.fUseEventCuts[eNoCollInRofStrict] = ec.fUseEventCuts[eNoCollInRofStrict] && cf_ec.cfUseNoCollInRofStrict;
  ec.fUseEventCuts[eNoCollInRofStandard] = ec.fUseEventCuts[eNoCollInRofStandard] && cf_ec.cfUseNoCollInRofStandard;
  ec.fUseEventCuts[eNoHighMultCollInPrevRof] = ec.fUseEventCuts[eNoHighMultCollInPrevRof] && cf_ec.cfUseNoHighMultCollInPrevRof;
  ec.fUseEventCuts[eIsGoodITSLayer3] = ec.fUseEventCuts[eIsGoodITSLayer3] && cf_ec.cfUseIsGoodITSLayer3;
  ec.fUseEventCuts[eIsGoodITSLayer0123] = ec.fUseEventCuts[eIsGoodITSLayer0123] && cf_ec.cfUseIsGoodITSLayer0123;
  ec.fUseEventCuts[eIsGoodITSLayersAll] = ec.fUseEventCuts[eIsGoodITSLayersAll] && cf_ec.cfUseIsGoodITSLayersAll;
  ec.fUseEventCuts[eCentralityWeights] = ec.fUseEventCuts[eCentralityWeights] && cf_cw.cfUseCentralityWeights;

  // **) event cuts defined via [min, max):
  //     Remark: I use this one also for events cuts set only via min or via max.
  //             In this case, I set either min or max intentionally to some value which never can be met (see below example for "MinVertexDistanceFromIP")
  auto lNumberOfEvents = (std::vector<int>)cf_ec.cfNumberOfEvents;
  ec.fdEventCuts[eNumberOfEvents][eMin] = lNumberOfEvents[eMin];
  ec.fdEventCuts[eNumberOfEvents][eMax] = lNumberOfEvents[eMax];

  auto lTotalMultiplicity = (std::vector<int>)cf_ec.cfTotalMultiplicity;
  ec.fdEventCuts[eTotalMultiplicity][eMin] = lTotalMultiplicity[eMin];
  ec.fdEventCuts[eTotalMultiplicity][eMax] = lTotalMultiplicity[eMax];

  auto lMultiplicity = (std::vector<float>)cf_ec.cfMultiplicity;
  ec.fdEventCuts[eMultiplicity][eMin] = lMultiplicity[eMin];
  ec.fdEventCuts[eMultiplicity][eMax] = lMultiplicity[eMax];

  auto lReferenceMultiplicity = (std::vector<float>)cf_ec.cfReferenceMultiplicity;
  ec.fdEventCuts[eReferenceMultiplicity][eMin] = lReferenceMultiplicity[eMin];
  ec.fdEventCuts[eReferenceMultiplicity][eMax] = lReferenceMultiplicity[eMax];

  auto lCentrality = (std::vector<float>)cf_ec.cfCentrality;
  ec.fdEventCuts[eCentrality][eMin] = lCentrality[eMin];
  ec.fdEventCuts[eCentrality][eMax] = lCentrality[eMax];

  auto lVertex_x = (std::vector<float>)cf_ec.cfVertex_x;
  ec.fdEventCuts[eVertex_x][eMin] = lVertex_x[eMin];
  ec.fdEventCuts[eVertex_x][eMax] = lVertex_x[eMax];

  auto lVertex_y = (std::vector<float>)cf_ec.cfVertex_y;
  ec.fdEventCuts[eVertex_y][eMin] = lVertex_y[eMin];
  ec.fdEventCuts[eVertex_y][eMax] = lVertex_y[eMax];

  auto lVertex_z = (std::vector<float>)cf_ec.cfVertex_z;
  ec.fdEventCuts[eVertex_z][eMin] = lVertex_z[eMin];
  ec.fdEventCuts[eVertex_z][eMax] = lVertex_z[eMax];

  auto lNContributors = (std::vector<int>)cf_ec.cfNContributors;
  ec.fdEventCuts[eNContributors][eMin] = lNContributors[eMin];
  ec.fdEventCuts[eNContributors][eMax] = lNContributors[eMax];

  auto lImpactParameter = (std::vector<float>)cf_ec.cfImpactParameter;
  ec.fdEventCuts[eImpactParameter][eMin] = lImpactParameter[eMin];
  ec.fdEventCuts[eImpactParameter][eMax] = lImpactParameter[eMax];

  auto lEventPlaneAngle = (std::vector<float>)cf_ec.cfEventPlaneAngle;
  ec.fdEventCuts[eEventPlaneAngle][eMin] = lEventPlaneAngle[eMin];
  ec.fdEventCuts[eEventPlaneAngle][eMax] = lEventPlaneAngle[eMax];

  auto lOccupancy = (std::vector<float>)cf_ec.cfOccupancy;
  ec.fdEventCuts[eOccupancy][eMin] = lOccupancy[eMin];
  ec.fdEventCuts[eOccupancy][eMax] = lOccupancy[eMax];

  auto lInteractionRate = (std::vector<float>)cf_ec.cfInteractionRate;
  ec.fdEventCuts[eInteractionRate][eMin] = lInteractionRate[eMin];
  ec.fdEventCuts[eInteractionRate][eMax] = lInteractionRate[eMax];

  auto lCurrentRunDuration = (std::vector<float>)cf_ec.cfCurrentRunDuration;
  ec.fdEventCuts[eCurrentRunDuration][eMin] = lCurrentRunDuration[eMin];
  ec.fdEventCuts[eCurrentRunDuration][eMax] = lCurrentRunDuration[eMax];

  auto lMultMCNParticlesEta08 = (std::vector<float>)cf_ec.cfMultMCNParticlesEta08;
  ec.fdEventCuts[eMultMCNParticlesEta08][eMin] = lMultMCNParticlesEta08[eMin];
  ec.fdEventCuts[eMultMCNParticlesEta08][eMax] = lMultMCNParticlesEta08[eMax];

  auto lSelectedEvents = (std::vector<int>)cf_ec.cfSelectedEvents;
  ec.fdEventCuts[eSelectedEvents][eMin] = lSelectedEvents[eMin];
  ec.fdEventCuts[eSelectedEvents][eMax] = lSelectedEvents[eMax];

  ec.fdEventCuts[eMinVertexDistanceFromIP][eMin] = cf_ec.cfMinVertexDistanceFromIP; // if vertex is closer to IP than this value, the event is rejected
  ec.fdEventCuts[eMinVertexDistanceFromIP][eMax] = -1;                              // // this value is never checked in any case

  // **) event cuts defined via string:
  ec.fsEventCuts[eMultiplicityEstimator] = cf_ec.cfMultiplicityEstimator;
  ec.fsEventCuts[eReferenceMultiplicityEstimator] = cf_ec.cfReferenceMultiplicityEstimator;
  ec.fsEventCuts[eCentralityEstimator] = cf_ec.cfCentralityEstimator;
  ec.fsEventCuts[eTrigger] = cf_ec.cfTrigger;
  ec.fsEventCuts[eOccupancyEstimator] = cf_ec.cfOccupancyEstimator;

  // ----------------------------------------------------------------------

  // b) Default particle cuts:

  // *) Use or do not use a cut enumerated in eParticleHistograms + eParticleCuts.
  //    Default cuts are set in configurable cfUseParticleCuts
  auto lUseParticleCuts = (std::vector<string>)cf_pc.cfUseParticleCuts;
  if (lUseParticleCuts.size() != eParticleCuts_N) {
    LOGF(info, "\033[1;31m lUseParticleCuts.size() = %d\033[0m", lUseParticleCuts.size());
    LOGF(info, "\033[1;31m eParticleCuts_N = %d\033[0m", static_cast<int>(eParticleCuts_N));
    LOGF(fatal, "\033[1;31m%s at line %d : Mismatch in the number of flags in configurable cfUseParticleCuts, and number of entries in enum eParticleHistograms + eParticleCuts \n \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Insanity check on the content and ordering of particle cuts in the initialization in configurable cfUseParticleCuts:
  // TBI 20240518 I do not need this in fact, I can automate initialization even without ordering in configurable, but it feels with the ordering enforced, it's much safer.
  for (int name = 0; name < eParticleCuts_N; name++) {
    // TBI 20240518 I could implement even a strickter EqualTo instead of EndsWith, but then I need to tokenize, etc., etc. This shall be safe enough.
    if (!TString(lUseParticleCuts[name]).EndsWith(pc.fParticleCutName[name].Data())) {
      LOGF(fatal, "\033[1;31m%s at line %d : Wrong content or ordering of contents in configurable cfUseParticleCuts => name = %d, lUseParticleCuts[name] = \"%s\", pc.fParticleCutName[name] = \"%s\" \033[0m", __FUNCTION__, __LINE__, name, TString(lUseParticleCuts[name]).Data(), pc.fParticleCutName[name].Data());
    }
  }

  // *) from enum eParticleHistograms:
  pc.fUseParticleCuts[ePhi] = Alright(lUseParticleCuts[ePhi]);
  pc.fUseParticleCuts[ePt] = Alright(lUseParticleCuts[ePt]);
  pc.fUseParticleCuts[eEta] = Alright(lUseParticleCuts[eEta]);
  pc.fUseParticleCuts[eCharge] = Alright(lUseParticleCuts[eCharge]);
  pc.fUseParticleCuts[etpcNClsFindable] = Alright(lUseParticleCuts[etpcNClsFindable]);
  pc.fUseParticleCuts[etpcNClsShared] = Alright(lUseParticleCuts[etpcNClsShared]);
  pc.fUseParticleCuts[eitsChi2NCl] = Alright(lUseParticleCuts[eitsChi2NCl]);
  pc.fUseParticleCuts[etpcNClsFound] = Alright(lUseParticleCuts[etpcNClsFound]);
  pc.fUseParticleCuts[etpcNClsCrossedRows] = Alright(lUseParticleCuts[etpcNClsCrossedRows]);
  pc.fUseParticleCuts[eitsNCls] = Alright(lUseParticleCuts[eitsNCls]);
  pc.fUseParticleCuts[eitsNClsInnerBarrel] = Alright(lUseParticleCuts[eitsNClsInnerBarrel]);
  pc.fUseParticleCuts[etpcCrossedRowsOverFindableCls] = Alright(lUseParticleCuts[etpcCrossedRowsOverFindableCls]);
  pc.fUseParticleCuts[etpcFoundOverFindableCls] = Alright(lUseParticleCuts[etpcFoundOverFindableCls]);
  pc.fUseParticleCuts[etpcFractionSharedCls] = Alright(lUseParticleCuts[etpcFractionSharedCls]);
  pc.fUseParticleCuts[etpcChi2NCl] = Alright(lUseParticleCuts[etpcChi2NCl]);
  pc.fUseParticleCuts[edcaXY] = Alright(lUseParticleCuts[edcaXY]);
  pc.fUseParticleCuts[edcaZ] = Alright(lUseParticleCuts[edcaZ]);
  pc.fUseParticleCuts[ePDG] = Alright(lUseParticleCuts[ePDG]);
  pc.fUseParticleCuts[etrackCutFlag] = Alright(lUseParticleCuts[etrackCutFlag]);
  pc.fUseParticleCuts[etrackCutFlagFb1] = Alright(lUseParticleCuts[etrackCutFlagFb1]);
  pc.fUseParticleCuts[etrackCutFlagFb2] = Alright(lUseParticleCuts[etrackCutFlagFb2]);
  pc.fUseParticleCuts[eisQualityTrack] = Alright(lUseParticleCuts[eisQualityTrack]);
  pc.fUseParticleCuts[eisPrimaryTrack] = Alright(lUseParticleCuts[eisPrimaryTrack]);
  pc.fUseParticleCuts[eisInAcceptanceTrack] = Alright(lUseParticleCuts[eisInAcceptanceTrack]);
  pc.fUseParticleCuts[eisGlobalTrack] = Alright(lUseParticleCuts[eisGlobalTrack]);
  pc.fUseParticleCuts[eisPVContributor] = Alright(lUseParticleCuts[eisPVContributor]);
  pc.fUseParticleCuts[ePtDependentDCAxyParameterization] = Alright(lUseParticleCuts[ePtDependentDCAxyParameterization]);

  // **) particles cuts defined via booleans:
  pc.fUseParticleCuts[etrackCutFlag] = pc.fUseParticleCuts[etrackCutFlag] && cf_pc.cftrackCutFlag;
  pc.fUseParticleCuts[etrackCutFlagFb1] = pc.fUseParticleCuts[etrackCutFlagFb1] && cf_pc.cftrackCutFlagFb1;
  pc.fUseParticleCuts[etrackCutFlagFb2] = pc.fUseParticleCuts[etrackCutFlagFb2] && cf_pc.cftrackCutFlagFb2;
  pc.fUseParticleCuts[eisQualityTrack] = pc.fUseParticleCuts[eisQualityTrack] && cf_pc.cfisQualityTrack;
  pc.fUseParticleCuts[eisPrimaryTrack] = pc.fUseParticleCuts[eisPrimaryTrack] && cf_pc.cfisPrimaryTrack;
  pc.fUseParticleCuts[eisInAcceptanceTrack] = pc.fUseParticleCuts[eisInAcceptanceTrack] && cf_pc.cfisInAcceptanceTrack;
  pc.fUseParticleCuts[eisGlobalTrack] = pc.fUseParticleCuts[eisGlobalTrack] && cf_pc.cfisGlobalTrack;
  pc.fUseParticleCuts[eisPVContributor] = pc.fUseParticleCuts[eisPVContributor] && cf_pc.cfisPVContributor;

  // **) particles cuts defined via [min, max):
  auto lPhi = (std::vector<float>)cf_pc.cfPhi;
  pc.fdParticleCuts[ePhi][eMin] = lPhi[eMin];
  pc.fdParticleCuts[ePhi][eMax] = lPhi[eMax];

  auto lPt = (std::vector<float>)cf_pc.cfPt;
  pc.fdParticleCuts[ePt][eMin] = lPt[eMin];
  pc.fdParticleCuts[ePt][eMax] = lPt[eMax];

  auto lEta = (std::vector<float>)cf_pc.cfEta;
  pc.fdParticleCuts[eEta][eMin] = lEta[eMin];
  pc.fdParticleCuts[eEta][eMax] = lEta[eMax];

  auto lCharge = (std::vector<float>)cf_pc.cfCharge;
  pc.fdParticleCuts[eCharge][eMin] = lCharge[eMin];
  pc.fdParticleCuts[eCharge][eMax] = lCharge[eMax];

  auto ltpcNClsFindable = (std::vector<float>)cf_pc.cftpcNClsFindable;
  pc.fdParticleCuts[etpcNClsFindable][eMin] = ltpcNClsFindable[eMin];
  pc.fdParticleCuts[etpcNClsFindable][eMax] = ltpcNClsFindable[eMax];

  auto ltpcNClsShared = (std::vector<float>)cf_pc.cftpcNClsShared;
  pc.fdParticleCuts[etpcNClsShared][eMin] = ltpcNClsShared[eMin];
  pc.fdParticleCuts[etpcNClsShared][eMax] = ltpcNClsShared[eMax];

  auto litsChi2NCl = (std::vector<float>)cf_pc.cfitsChi2NCl;
  pc.fdParticleCuts[eitsChi2NCl][eMin] = litsChi2NCl[eMin];
  pc.fdParticleCuts[eitsChi2NCl][eMax] = litsChi2NCl[eMax];

  auto ltpcNClsFound = (std::vector<float>)cf_pc.cftpcNClsFound;
  pc.fdParticleCuts[etpcNClsFound][eMin] = ltpcNClsFound[eMin];
  pc.fdParticleCuts[etpcNClsFound][eMax] = ltpcNClsFound[eMax];

  auto ltpcNClsCrossedRows = (std::vector<float>)cf_pc.cftpcNClsCrossedRows;
  pc.fdParticleCuts[etpcNClsCrossedRows][eMin] = ltpcNClsCrossedRows[eMin];
  pc.fdParticleCuts[etpcNClsCrossedRows][eMax] = ltpcNClsCrossedRows[eMax];

  auto litsNCls = (std::vector<float>)cf_pc.cfitsNCls;
  pc.fdParticleCuts[eitsNCls][eMin] = litsNCls[eMin];
  pc.fdParticleCuts[eitsNCls][eMax] = litsNCls[eMax];

  auto litsNClsInnerBarrel = (std::vector<float>)cf_pc.cfitsNClsInnerBarrel;
  pc.fdParticleCuts[eitsNClsInnerBarrel][eMin] = litsNClsInnerBarrel[eMin];
  pc.fdParticleCuts[eitsNClsInnerBarrel][eMax] = litsNClsInnerBarrel[eMax];

  auto ltpcCrossedRowsOverFindableCls = (std::vector<float>)cf_pc.cftpcCrossedRowsOverFindableCls;
  pc.fdParticleCuts[etpcCrossedRowsOverFindableCls][eMin] = ltpcCrossedRowsOverFindableCls[eMin];
  pc.fdParticleCuts[etpcCrossedRowsOverFindableCls][eMax] = ltpcCrossedRowsOverFindableCls[eMax];

  auto ltpcFoundOverFindableCls = (std::vector<float>)cf_pc.cftpcFoundOverFindableCls;
  pc.fdParticleCuts[etpcFoundOverFindableCls][eMin] = ltpcFoundOverFindableCls[eMin];
  pc.fdParticleCuts[etpcFoundOverFindableCls][eMax] = ltpcFoundOverFindableCls[eMax];

  auto ltpcFractionSharedCls = (std::vector<float>)cf_pc.cftpcFractionSharedCls;
  pc.fdParticleCuts[etpcFractionSharedCls][eMin] = ltpcFractionSharedCls[eMin];
  pc.fdParticleCuts[etpcFractionSharedCls][eMax] = ltpcFractionSharedCls[eMax];

  auto ltpcChi2NCl = (std::vector<float>)cf_pc.cftpcChi2NCl;
  pc.fdParticleCuts[etpcChi2NCl][eMin] = ltpcChi2NCl[eMin];
  pc.fdParticleCuts[etpcChi2NCl][eMax] = ltpcChi2NCl[eMax];

  auto ldcaXY = (std::vector<float>)cf_pc.cfdcaXY;
  pc.fdParticleCuts[edcaXY][eMin] = ldcaXY[eMin];
  pc.fdParticleCuts[edcaXY][eMax] = ldcaXY[eMax];

  auto ldcaZ = (std::vector<float>)cf_pc.cfdcaZ;
  pc.fdParticleCuts[edcaZ][eMin] = ldcaZ[eMin];
  pc.fdParticleCuts[edcaZ][eMax] = ldcaZ[eMax];

  auto lPDG = (std::vector<float>)cf_pc.cfPDG;
  pc.fdParticleCuts[ePDG][eMin] = lPDG[eMin];
  pc.fdParticleCuts[ePDG][eMax] = lPDG[eMax];

  // **) particles cuts defined via string:
  pc.fsParticleCuts[ePtDependentDCAxyParameterization] = cf_pc.cfPtDependentDCAxyParameterization;

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void DefaultCuts()

//============================================================

void SpecificCuts(TString whichSpecificCuts)
{
  // After default cuts are applied, on top of them apply analysis-specific cuts. Has to be called after DefaultBinning() and DefaultCuts().
  // Here I hardwire defalt cuts and settings for a given period which will overwrite whatever is set in configurables.
  // When I do systematic checks, this option shall NOT be used, because values for some cuts which I plan to vary, are also hardwired here.
  // Both event and particle cuts are hardwired here. As well as some other settings.
  // For the time being, all specific cuts are defaulted and tuned for the latest reconstruction pass.

  // a) Mapping;
  // b) Implementation of analysis-specific cuts.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Mapping:
  //    I need this mapping, for the switch statement below. TBI 20241120 I could introduce a utility function for this as well...
  eSpecificCuts specificCuts = eSpecificCuts_N;
  if (whichSpecificCuts.EqualTo("LHC23zzh")) {
    specificCuts = eLHC23zzh;
  } else if (whichSpecificCuts.EqualTo("LHC24ar")) {
    specificCuts = eLHC24ar;
  } else if (whichSpecificCuts.EqualTo("LHC24as")) {
    specificCuts = eLHC24as;
  } else if (whichSpecificCuts.EqualTo("LHC15o")) {
    specificCuts = eLHC15o;
  } else {
    LOGF(fatal, "\033[1;31m%s at line %d : whichSpecificCuts = %s is not supported \033[0m", __FUNCTION__, __LINE__, whichSpecificCuts.Data());
  }

  // b) Implementation of analysis-specific cuts:
  //    Remark #1: Whichever cuts start to repeat below across different case statements, promote them into DefaultCuts(), i.e. hardwire those values in configurables.
  //               The idea is to keep here cuts only which are specific for particular analysis, and which are unlikely ever to change as a default cut for that particular analysis.
  //    Remark #2: Remember that the values for the cuts hardwired here overwrite the ones set as default values in configurables.
  //               If you want to reconfigure all cuts below manually via configurables, simply do not call SpecificCuts, i.e. set in JSON "cfUseSpecificCuts": "false"
  //               Therefore, if I want to vary some of these cuts via configurables as a part of systematics, I must set in JSON "cfUseSpecificCuts": "false"
  //    Remark #3: Most up-to-date documentation of each cut is in enum file.
  switch (specificCuts) {

    case eLHC23zzh:

      // In this branch I implement default cuts and settings for PbPb Run 3 datasets collected in 2023.
      // If I change some cut here, keep in sync. with other branches (e.g. for 2024 data).

      // Event cuts:
      ec.fUseEventCuts[eSel7] = false;
      ec.fUseEventCuts[eSel8] = true;
      ec.fUseEventCuts[eNoSameBunchPileup] = true;
      ec.fUseEventCuts[eIsVertexITSTPC] = true;
      ec.fUseEventCuts[eNoCollInTimeRangeStandard] = true;
      ec.fUseEventCuts[eNoCollInTimeRangeStrict] = false;
      ec.fUseEventCuts[eNoCollInRofStandard] = true;
      ec.fUseEventCuts[eNoCollInRofStrict] = false;

      // Particle cuts:
      pc.fUseParticleCuts[eitsNCls] = true;
      pc.fdParticleCuts[eitsNCls][eMin] = 5.;
      pc.fdParticleCuts[eitsNCls][eMax] = 1000.;

      pc.fUseParticleCuts[etpcNClsFound] = true;
      pc.fdParticleCuts[etpcNClsFound][eMin] = 70.;
      pc.fdParticleCuts[etpcNClsFound][eMax] = 1000.;

      pc.fUseParticleCuts[etpcNClsCrossedRows] = true;
      pc.fdParticleCuts[etpcNClsCrossedRows][eMin] = 70.;
      pc.fdParticleCuts[etpcNClsCrossedRows][eMax] = 1000.;

      pc.fUseParticleCuts[etpcCrossedRowsOverFindableCls] = true;
      pc.fdParticleCuts[etpcCrossedRowsOverFindableCls][eMin] = 0.8;
      pc.fdParticleCuts[etpcCrossedRowsOverFindableCls][eMax] = 1000.;

      pc.fUseParticleCuts[etpcFoundOverFindableCls] = true;
      pc.fdParticleCuts[etpcFoundOverFindableCls][eMin] = 0.8;
      pc.fdParticleCuts[etpcFoundOverFindableCls][eMax] = 1000.;

      pc.fUseParticleCuts[etpcFractionSharedCls] = true;
      pc.fdParticleCuts[etpcFractionSharedCls][eMin] = -1000.;
      pc.fdParticleCuts[etpcFractionSharedCls][eMax] = 0.4;

      pc.fUseParticleCuts[etpcChi2NCl] = true;
      pc.fdParticleCuts[etpcChi2NCl][eMin] = -1000.;
      pc.fdParticleCuts[etpcChi2NCl][eMax] = 4.0;

      pc.fUseParticleCuts[edcaXY] = true;
      pc.fdParticleCuts[edcaXY][eMin] = -2.4;
      pc.fdParticleCuts[edcaXY][eMax] = 2.4;

      pc.fUseParticleCuts[edcaZ] = true;
      pc.fdParticleCuts[edcaZ][eMin] = -3.2;
      pc.fdParticleCuts[edcaZ][eMax] = 3.2;

      pc.fUseParticleCuts[eisInAcceptanceTrack] = false; // see enum
      pc.fUseParticleCuts[eisGlobalTrack] = false;       // only for Run 2
      pc.fUseParticleCuts[eisPVContributor] = true;

      break;

    case eLHC24ar:
    case eLHC24as:

      // In this branch I implement default cuts and settings for PbPb Run 3 datasets collected in 2024:
      // If I change some cut here, keep in sync. with other branches (e.g. for 2023 data).
      // As of 20250207, all cuts are the same as for 2023, expect that here I do NOT use eIsGoodZvtxFT0vsPV and eNoHighMultCollInPrevRof

      // Event cuts:
      ec.fUseEventCuts[eSel7] = false;
      ec.fUseEventCuts[eSel8] = true;
      ec.fUseEventCuts[eNoSameBunchPileup] = true;
      ec.fUseEventCuts[eIsVertexITSTPC] = true;
      ec.fUseEventCuts[eNoCollInTimeRangeStandard] = true;
      ec.fUseEventCuts[eNoCollInTimeRangeStrict] = false;
      ec.fUseEventCuts[eNoCollInRofStandard] = true;
      ec.fUseEventCuts[eNoCollInRofStrict] = false;
      ec.fUseEventCuts[eIsGoodZvtxFT0vsPV] = false;       // diff commpared to 2023
      ec.fUseEventCuts[eNoHighMultCollInPrevRof] = false; // diff commpared to 2023

      // Particle cuts:
      pc.fUseParticleCuts[eitsNCls] = true;
      pc.fdParticleCuts[eitsNCls][eMin] = 5.;
      pc.fdParticleCuts[eitsNCls][eMax] = 1000.;

      pc.fUseParticleCuts[etpcNClsFound] = true;
      pc.fdParticleCuts[etpcNClsFound][eMin] = 70.;
      pc.fdParticleCuts[etpcNClsFound][eMax] = 1000.;

      pc.fUseParticleCuts[etpcNClsCrossedRows] = true;
      pc.fdParticleCuts[etpcNClsCrossedRows][eMin] = 70.;
      pc.fdParticleCuts[etpcNClsCrossedRows][eMax] = 1000.;

      pc.fUseParticleCuts[etpcCrossedRowsOverFindableCls] = true;
      pc.fdParticleCuts[etpcCrossedRowsOverFindableCls][eMin] = 0.8;
      pc.fdParticleCuts[etpcCrossedRowsOverFindableCls][eMax] = 1000.;

      pc.fUseParticleCuts[etpcFoundOverFindableCls] = true;
      pc.fdParticleCuts[etpcFoundOverFindableCls][eMin] = 0.8;
      pc.fdParticleCuts[etpcFoundOverFindableCls][eMax] = 1000.;

      pc.fUseParticleCuts[etpcFractionSharedCls] = true;
      pc.fdParticleCuts[etpcFractionSharedCls][eMin] = -1000.;
      pc.fdParticleCuts[etpcFractionSharedCls][eMax] = 0.4;

      pc.fUseParticleCuts[etpcChi2NCl] = true;
      pc.fdParticleCuts[etpcChi2NCl][eMin] = -1000.;
      pc.fdParticleCuts[etpcChi2NCl][eMax] = 4.0;

      pc.fUseParticleCuts[edcaXY] = true;
      pc.fdParticleCuts[edcaXY][eMin] = -2.4;
      pc.fdParticleCuts[edcaXY][eMax] = 2.4;

      pc.fUseParticleCuts[edcaZ] = true;
      pc.fdParticleCuts[edcaZ][eMin] = -3.2;
      pc.fdParticleCuts[edcaZ][eMax] = 3.2;

      pc.fUseParticleCuts[eisInAcceptanceTrack] = false; // see enum
      pc.fUseParticleCuts[eisGlobalTrack] = false;       // only for Run 2
      pc.fUseParticleCuts[eisPVContributor] = true;

      break;

    case eLHC15o:

      // In this branch I implement default cuts and settings for Run 2 datasets:

      // Event cuts:
      // ec.fUseEventCuts[eSel7] = true; // TBI 20250115 ehen i procees in "Rec" some converted Run 2 MC, it removes 99% of events, see enum
      ec.fUseEventCuts[eSel8] = false;
      ec.fUseEventCuts[eNoSameBunchPileup] = false;
      ec.fUseEventCuts[eIsGoodZvtxFT0vsPV] = false;
      ec.fUseEventCuts[eIsVertexITSTPC] = false;
      ec.fUseEventCuts[eNoCollInTimeRangeStrict] = false;
      ec.fUseEventCuts[eNoCollInRofStrict] = false;
      ec.fUseEventCuts[eNoHighMultCollInPrevRof] = false;
      ec.fUseEventCuts[eNoCollInTimeRangeStandard] = false;
      ec.fUseEventCuts[eNoCollInRofStrict] = false;
      ec.fUseEventCuts[eNoCollInRofStandard] = false;
      ec.fUseEventCuts[eNoCollInRofStandard] = false;
      ec.fUseEventCuts[eIsGoodITSLayer3] = false;
      ec.fUseEventCuts[eIsGoodITSLayer0123] = false;
      ec.fUseEventCuts[eIsGoodITSLayersAll] = false;

      // ec.fUseEventCuts[eTrigger] = true;
      // ec.fsEventCuts[eTrigger] = "kINT7"; // TBI 20250115 cannot be used when i procees in "Rec" some converted Run 2 MC, see enum

      // ...

      // Particle cuts:
      pc.fUseParticleCuts[eisInAcceptanceTrack] = false; // see enum
      pc.fUseParticleCuts[etrackCutFlagFb1] = false;     // only for Run 3
      pc.fUseParticleCuts[etrackCutFlagFb2] = false;     // only for Run 3
      pc.fUseParticleCuts[eisPVContributor] = false;     // only for Run 3

      // ...

      // The rest:
      mupa.fCalculateCorrelationsAsFunctionOf[AFO_OCCUPANCY] = false;
      mupa.fCalculateCorrelationsAsFunctionOf[AFO_INTERACTIONRATE] = false;
      mupa.fCalculateCorrelationsAsFunctionOf[AFO_CURRENTRUNDURATION] = false;

      t0.fCalculateTest0AsFunctionOf[AFO_OCCUPANCY] = false;
      t0.fCalculateTest0AsFunctionOf[AFO_INTERACTIONRATE] = false;
      t0.fCalculateTest0AsFunctionOf[AFO_CURRENTRUNDURATION] = false;

      es.fCalculateEtaSeparationsAsFunctionOf[AFO_OCCUPANCY] = false;
      es.fCalculateEtaSeparationsAsFunctionOf[AFO_INTERACTIONRATE] = false;
      es.fCalculateEtaSeparationsAsFunctionOf[AFO_CURRENTRUNDURATION] = false;

      // ...

      break;

      // ...

    default:
      LOGF(fatal, "\033[1;31m%s at line %d : specificCuts = %d is not supported yet \033[0m", __FUNCTION__, __LINE__, static_cast<int>(specificCuts));
      break;
  } // switch((specificCuts))

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void SpecificCuts(const char* specificCutsName)

//============================================================

void InsanityChecksOnDefinitionsOfConfigurables()
{
  // Do insanity checks on values obtained from configurables before using them in the remaining function.
  // This is really important, because one misconfigured configurable (e.g. boolean set to string), causes the whole json config to die silently, and
  // only default values from MuPa-Configurables.h are used.
  // Here I only check if configurables are correctly defined, I do NOT here initialize local variables with configurables, that is done later.
  // Example misconfiguration in JSON:
  //       "var": "true",   =>   var = 1 + other configurables are processed correctly
  //       "var": "truee",  =>   var = 0 + all settings in JSON for configurables are ingored silently

  // TBI 20241127 finalize this function eventually. This is not urgent, though, as only a check below on cfTaskIsConfiguredFromJson covers most cases already.

  // Remark: Ordering below reflects the ordering in Configurables.h, not in DataMembers.h
  // a) Task configuration;
  // b) QA;
  // c) Event histograms;
  // d) Event cuts;
  // e) Particle histograms;
  // f) Particle cuts;
  // g) Q-vectors;
  // h) Multiparticle correlations;
  // i) Test0;
  // j) Particle weights;
  // k) Centrality weights;
  // l) Nested loops;
  // m) Toy NUA;
  // n) Internal validation;
  // o) Results histograms.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Task configuration:
  if (!TString(cf_tc.cfTaskIsConfiguredFromJson).EqualTo("yes")) {
    LOGF(fatal, "\033[1;31m%s at line %d : configurable cfTaskIsConfiguredFromJson = \"%s\", but it has to be set to \"yes\" in JSON => most likely some other configurable is misconfigured and all remaining settings in JSON are ignored silently\033[0m", __FUNCTION__, __LINE__, TString(cf_tc.cfTaskIsConfiguredFromJson).Data());
  }

  // b) QA:
  // ...

  // c) Event histograms:
  // ...

  // d) Event cuts:
  // ...

  // e) Particle histograms:
  // ...

  // f) Particle cuts:
  // ...

  // g) Q-vectors:
  // ...

  // h) Multiparticle correlations:
  // ...

  // i) Test0:
  // ...

  // j) Particle weights:
  // ...

  // k) Centrality weights:
  // ...

  // l) Nested loops:
  // ...

  // m) Toy NUA:
  // ...

  // n) Internal validation:
  // ...

  // o) Results histograms:
  // ...

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // InsanityChecksOnDefinitionsOfConfigurables()

//============================================================

void InsanityChecksBeforeBooking()
{
  // Do insanity checks on configuration, binning and cuts. Values obtained from configurables are checked before being used in InsanityChecksOnDefinitionsOfConfigurables().
  // Remember that here I cannot do insanity checks on local histograms, etc., because they are not booked yet.
  // For those additional checks, use InsanityChecksAfterBooking().

  // a) Insanity checks on configuration;
  // b) Ensure that Run 1/2 specific cuts and flags are used only in Run 1/2 (both data and sim);
  // c) Ensure that Run 3 specific cuts and flags are used only in Run 3 (both data and sim);
  // d) Insanity checks on binning;
  // e) Insanity checks on events cuts;
  // f) Insanity checks on Toy NUA;
  // g) Insanity checks on internal validation;
  // h) Insanity checks on results histograms.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Insanity checks on configuration:

  // **) Dry run and internal validation are not meant to be run together:
  if (tc.fDryRun && iv.fUseInternalValidation) {
    LOGF(fatal, "\033[1;31m%s at line %d : Dry run and internal validation are not meant to be run together\033[0m", __FUNCTION__, __LINE__);
  }

  // **) Cannot calculate multiparticle correlations, in case Q-vectors are not filled:
  if (mupa.fCalculateCorrelations && !qv.fCalculateQvectors) {
    LOGF(fatal, "\033[1;31m%s at line %d : Cannot calculate multiparticle correlations, in case Q-vectors are not calculated \033[0m", __FUNCTION__, __LINE__);
  }

  // **) If some differential "correlations" flag is set to true, but the main fCalculateCorrelations is false, only print the warning that that differential correlations won't be calculated.
  //     This is not fatal, because this way I can turn off all differential "correlations" flags, just by setting fCalculateCorrelations to false, e.g. when I want to fill only control histograms.
  for (int v = 0; v < eAsFunctionOf_N; v++) {
    if (mupa.fCalculateCorrelationsAsFunctionOf[v] && !mupa.fCalculateCorrelations) {
      LOGF(warning, "\033[1;33m%s at line %d : mupa.fCalculateCorrelationsAsFunctionOf[%d] is true, but mupa.fCalculateCorrelations is false. This differential correlations won't be calculated.\033[0m", __FUNCTION__, __LINE__, v);
    }
  }

  // **) Cannot calculate Test0, in case Q-vectors are not filled:
  if (t0.fCalculateTest0 && !qv.fCalculateQvectors) {
    LOGF(fatal, "\033[1;31m%s at line %d : Cannot calculate Test0, in case Q-vectors are not filled \033[0m", __FUNCTION__, __LINE__);
  }

  // **) If some differential Test0 flag is set to true, but the main fCalculateTest0 is false, only print the warning that that differential Test0 won't be calculated.
  //     This is not fatal, because this way I can turn off all differential Test0 flags, just by setting fCalculateTest0 to false, e.g. when I want to fill only control histograms.
  for (int v = 0; v < eAsFunctionOf_N; v++) {
    if (t0.fCalculateTest0AsFunctionOf[v] && !t0.fCalculateTest0) {
      LOGF(warning, "\033[1;33m%s at line %d : t0.fCalculateTest0AsFunctionOf[%d] is true, but t0.fCalculateTest0 is false. This differential Test0 won't be calculated.\033[0m", __FUNCTION__, __LINE__, v);
    }
  }

  // **) Enforce that if fixed number of randomly selected tracks is used in the analysis, that Fisher-Yates algorithm is enabled:
  if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && !tc.fUseFisherYates) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // **) When it comes to DCAxy cut, ensure that either flat or pt-dependent cut is used, but not both:
  if (pc.fUseParticleCuts[edcaXY] && pc.fUseParticleCuts[ePtDependentDCAxyParameterization]) {
    LOGF(fatal, "\033[1;31m%s at line %d : use either flat or pt-dependent DCAxy cut, but not both \033[0m", __FUNCTION__, __LINE__);
  }

  // **) Insanity check on individual flags: Make sure that only one process is set to true.
  //     If 2 or more are true, then corresponding process function is executed over ALL data, then another process(...) function, etc.
  //     Re-think this if it's possible to run different process(...)'s concurently over the same data.
  if (static_cast<int>(tc.fProcess[eProcessRec]) + static_cast<int>(tc.fProcess[eProcessRecSim]) + static_cast<int>(tc.fProcess[eProcessSim]) + static_cast<int>(tc.fProcess[eProcessRec_Run2]) + static_cast<int>(tc.fProcess[eProcessRecSim_Run2]) + static_cast<int>(tc.fProcess[eProcessSim_Run2]) + static_cast<int>(tc.fProcess[eProcessRec_Run1]) + static_cast<int>(tc.fProcess[eProcessRecSim_Run1]) + static_cast<int>(tc.fProcess[eProcessSim_Run1]) > 1) {
    LOGF(info, "\033[1;31m Only one flag can be true: tc.fProcess[eProcessRec] = %d, tc.fProcess[eProcessRecSim] = %d, tc.fProcess[eProcessSim] = %d, tc.fProcess[eProcessRec_Run2] = %d, tc.fProcess[eProcessRecSim_Run2] = %d, tc.fProcess[eProcessSim_Run2] = %d, tc.fProcess[eProcessRec_Run1] = %d, tc.fProcess[eProcessRecSim_Run1] = %d, tc.fProcess[eProcessSim_Run1] = %d \033[0m", static_cast<int>(tc.fProcess[eProcessRec]), static_cast<int>(tc.fProcess[eProcessRecSim]), static_cast<int>(tc.fProcess[eProcessSim]), static_cast<int>(tc.fProcess[eProcessRec_Run2]), static_cast<int>(tc.fProcess[eProcessRecSim_Run2]), static_cast<int>(tc.fProcess[eProcessSim_Run2]), static_cast<int>(tc.fProcess[eProcessRec_Run1]), static_cast<int>(tc.fProcess[eProcessRecSim_Run1]), static_cast<int>(tc.fProcess[eProcessSim_Run1]));
    LOGF(fatal, "in function \033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // **) Insanity checks on event cuts:

  // **) This check is meant to prevent the case when I want to bailout for max number of events, but I do not fill event histograms:
  if (ec.fdEventCuts[eNumberOfEvents][eMax] < 1e6) { // TBI 20241011 Do I need to tune 1000000000
    // If I do not want to bail out when max number of events is reached, then in the configurable I have e.g. cfNumberOfEvents{"cfNumberOfEvents", {-1, 1000000000}
    // So if the upper limit is set to some number < 1e6, I want to bail out for that number of events.
    // TBI 20241011 this is a bit shaky, but nevermind now...
    if (!eh.fBookEventHistograms[eNumberOfEvents]) {
      LOGF(fatal, "\033[1;31m%s at line %d : Bailout for max number of events cannot be done, unless eh.fBookEventHistograms[eNumberOfEvents] is true.\033[0m", __FUNCTION__, __LINE__);
    }
  }

  // **) Check if the trigger makes sense or was validated for this dataset:
  if (ec.fUseEventCuts[eTrigger]) {

    // Validated and supported Run 3 triggers:
    if (tc.fProcess[eProcessRec]) {
      if (!ec.fsEventCuts[eTrigger].EqualTo("kTVXinTRD")) {
        LOGF(fatal, "\033[1;31m%s at line %d : trigger \"%s\" is not internally validated or supported for Run 3. Add it to the list of supported triggers, if you really want to use that one.\033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eTrigger].Data());
      }
    }

    // Validated and supported Run 2 triggers:
    if (tc.fProcess[eProcessRec_Run2]) {
      if (!ec.fsEventCuts[eTrigger].EqualTo("kINT7")) {
        //        LOGF(fatal, "\033[1;31m%s at line %d : trigger \"%s\" is not internally validated/supported yet for Run 2. Add it to the list of supported triggers, if you really want to use that one.\033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eTrigger].Data());
      }
    }

    // Validated and supported Run 1 triggers:
    //  ...

  } // if (ec.fUseEventCuts[eTrigger]) {

  // **) Check if the cut on MinVertexDistanceFromIP makes sense:
  if (ec.fUseEventCuts[eMinVertexDistanceFromIP]) {
    if (!(ec.fdEventCuts[eMinVertexDistanceFromIP][eMin] > 0.)) {
      LOGF(fatal, "\033[1;31m%s at line %d : trigger ec.fdEventCuts[eMinVertexDistanceFromIP][eMin] = %f must be positive. Check the setting of configurable cfMinVertexDistanceFromIP\033[0m", __FUNCTION__, __LINE__, ec.fdEventCuts[eMinVertexDistanceFromIP][eMin]);
    }
  }

  // **) Enforce the usage of particular trigger for this dataset:
  if (tc.fProcess[eProcessRec_Run2]) {
    // TBI 20250115 Not really sure I need this - if I want to run only "Rec" over Monte Carlo, then obviously the condition below is pointless.
    //              Also here I need to be able automaticaly to determine whether I am processing real data or Monte Carlo, from the dataset itself.
    // TBI 20240517 for the time being, here I am enforcing that "kINT7" is mandatory for Run 2
    // TBI 20241209 I still have to validate it for Run 1 converted real data => then expand if(...) statement above

    /* commented out temporariy, see TBI 20250115 above
    if (!(ec.fUseEventCuts[eTrigger] && ec.fsEventCuts[eTrigger].EqualTo("kINT7"))) {
      LOGF(fatal, "\033[1;31m%s at line %d : trigger \"%s\" is not internally validated/supported yet. Add it to the list of supported triggers, if you really want to use that one.\033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eTrigger].Data());
    } else {
      LOGF(info, "\033[1;32m%s at line %d : WARNING => trigger \"%s\" can be used only on real converted Run 2 and Run 1 data. For MC converted Run 2 and Run 1 data, this trigger shouldn't be used.\033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eTrigger].Data());
      // TBI 20240517 I need here programmatic access to "event-selection-task" flags "isMC and "isRunMC" . Then I can directly bail out.
    }
    */
  }

  // **) Ensure that fFloatingPointPrecision makes sense:
  if (!(tc.fFloatingPointPrecision > 0.)) {
    LOGF(fatal, "\033[1;31m%s at line %d : set fFloatingPointPrecision = %f to some small positive value, which will determine if two floats are the same \033[0m", __FUNCTION__, __LINE__, tc.fFloatingPointPrecision);
  }

  // **) Ensure that fSequentialBailout makes sense:
  if (!(tc.fSequentialBailout >= 0)) {
    LOGF(fatal, "\033[1;31m%s at line %d : set fSequentialBailout = %d either to 0 (not used), or to positive integer.\033[0m", __FUNCTION__, __LINE__, tc.fSequentialBailout);
  }

  // **) Ensure that I do not spill over with number of dimensions in sparse histograms:
  if (eDiffPhiWeights_N > gMaxNumberSparseDimensions) {
    LOGF(fatal, "\033[1;31m%s at line %d : set eDiffPhiWeights_N = %d is bigger than gMaxNumberSparseDimensions = %d\033[0m", __FUNCTION__, __LINE__, static_cast<int>(eDiffPhiWeights_N), gMaxNumberSparseDimensions);
  }
  if (eDiffPtWeights_N > gMaxNumberSparseDimensions) {
    LOGF(fatal, "\033[1;31m%s at line %d : set eDiffPtWeights_N = %d is bigger than gMaxNumberSparseDimensions = %d\033[0m", __FUNCTION__, __LINE__, static_cast<int>(eDiffPtWeights_N), gMaxNumberSparseDimensions);
  }
  if (eDiffEtaWeights_N > gMaxNumberSparseDimensions) {
    LOGF(fatal, "\033[1;31m%s at line %d : set eDiffEtaWeights_N = %d is bigger than gMaxNumberSparseDimensions = %d\033[0m", __FUNCTION__, __LINE__, static_cast<int>(eDiffEtaWeights_N), gMaxNumberSparseDimensions);
  }

  // b) Ensure that Run 1/2 specific cuts and flags are used only in Run 1/2 (both data and sim):
  // **) Ensure that eSel7 is used only for converted Run 2 and Run 1 (both data and sim):
  if (ec.fUseEventCuts[eSel7]) {
    if (!(tc.fProcess[eProcessRec_Run2] || tc.fProcess[eProcessRecSim_Run2] || tc.fProcess[eProcessSim_Run2] || tc.fProcess[eProcessRec_Run1] || tc.fProcess[eProcessRecSim_Run1] || tc.fProcess[eProcessSim_Run1])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use fSel7 for Run 2 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  // **) Supported reference multiplicity estimators for Run 1 and 2 are enlisted here:
  if (tc.fProcess[eProcessRec_Run2] || tc.fProcess[eProcessRecSim_Run2] || tc.fProcess[eProcessRec_Run1] || tc.fProcess[eProcessRecSim_Run1]) {
    if (!(ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("MultTracklets", TString::kIgnoreCase))) {
      LOGF(fatal, "\033[1;31m%s at line %d : reference multiplicity  estimator = %s is not supported yet for Run 1 and 2 analysis.\nUse \"MultTracklets\"\033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eReferenceMultiplicityEstimator].Data());
    }
  } else if (tc.fProcess[eProcessSim_Run2] || tc.fProcess[eProcessSim_Run1]) {
    LOGF(fatal, "\033[1;31m%s at line %d : eProcessSim is not validated yet \033[0m", __FUNCTION__, __LINE__);
  }

  // **) Supported centrality estimators for Run 1 and 2 are enlisted here:
  if (tc.fProcess[eProcessRec_Run2] || tc.fProcess[eProcessRecSim_Run2] || tc.fProcess[eProcessSim_Run2] || tc.fProcess[eProcessRec_Run1] || tc.fProcess[eProcessRecSim_Run1] || tc.fProcess[eProcessSim_Run1]) {
    if (!(ec.fsEventCuts[eCentralityEstimator].EqualTo("centRun2V0M", TString::kIgnoreCase) ||
          ec.fsEventCuts[eCentralityEstimator].EqualTo("centRun2SPDTracklets", TString::kIgnoreCase))) {
      LOGF(fatal, "\033[1;31m%s at line %d : centrality estimator = %s is not supported yet for converted Run 2 and Run 1 analysis.\nUse either \"centRun2V0M\" or \"centRun2SPDTracklets\" (case sensitive!) \033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eCentralityEstimator].Data());
    }
  }

  // **) Protection against particle cuts which are available, but not yet validated, or are meaningless, in Run 2 and 1:
  if (tc.fProcess[eProcessRec_Run2] || tc.fProcess[eProcessRecSim_Run2] || tc.fProcess[eProcessSim_Run2] || tc.fProcess[eProcessRec_Run1] || tc.fProcess[eProcessRecSim_Run1] || tc.fProcess[eProcessSim_Run1]) {
    if (pc.fUseParticleCuts[etrackCutFlag]) {
      LOGF(fatal, "\033[1;31m%s at line %d : particle cut etrackCutFlag is not validated, as of 20250113 it has no effect in Run 2 and Run 1 \033[0m", __FUNCTION__, __LINE__);
    }
    if (pc.fUseParticleCuts[etrackCutFlagFb1]) {
      LOGF(fatal, "\033[1;31m%s at line %d : particle cut etrackCutFlagFb1 is not validated, as of 20250113 it kills all reconstructed tracks in Run 2 and Run 1 \033[0m", __FUNCTION__, __LINE__);
    }
    if (pc.fUseParticleCuts[etrackCutFlagFb2]) {
      LOGF(fatal, "\033[1;31m%s at line %d : particle cut etrackCutFlagFb2 is not validated, as of 20250113 it kills all reconstructed tracks in Run 2 and Run 1 \033[0m", __FUNCTION__, __LINE__);
    }
  }

  // ...

  // c) Ensure that Run 3 specific cuts and flags are used only in Run 3 (both data and sim):
  // **) Ensure that eSel8 is used only in Run 3 (both data and sim):
  if (ec.fUseEventCuts[eSel8]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eSel8 only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eNoSameBunchPileup]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eNoSameBunchPileup only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eIsGoodZvtxFT0vsPV]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eIsGoodZvtxFT0vsPV only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eIsVertexITSTPC]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eIsVertexITSTPC only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eIsVertexTOFmatched]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eIsVertexTOFmatched only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eIsVertexTRDmatched]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eIsVertexTRDmatched only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eNoCollInTimeRangeStrict]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eNoCollInTimeRangeStrict only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eNoCollInTimeRangeStandard]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eNoCollInTimeRangeStandard only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eNoCollInRofStrict]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eNoCollInRofStrict only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eNoCollInRofStandard]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eNoCollInRofStandard only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eNoHighMultCollInPrevRof]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eNoHighMultCollInPrevRof only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eIsGoodITSLayer3]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eIsGoodITSLayer3 only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eIsGoodITSLayer0123]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eIsGoodITSLayer0123 only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (ec.fUseEventCuts[eIsGoodITSLayersAll]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : use eIsGoodITSLayersAll only for Run 3 data and MC\033[0m", __FUNCTION__, __LINE__);
    }
  }

  // **) Supported reference multiplicity estimators for Run 3 are enlisted here:
  if (tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim]) {
    if (!(ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("MultTPC", TString::kIgnoreCase) ||
          ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("MultFV0M", TString::kIgnoreCase) ||
          ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("MultFT0C", TString::kIgnoreCase) ||
          ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("MultFT0M", TString::kIgnoreCase) ||
          ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("MultNTracksPV", TString::kIgnoreCase))) {
      LOGF(fatal, "\033[1;31m%s at line %d : reference multiplicity  estimator = %s is not supported yet for Run 3 analysis.\nUse \"MultTPC\", \"MultFV0M\", \"MultFT0C\", \"MultFT0M\" or \"MultNTracksPV\"\033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eReferenceMultiplicityEstimator].Data());
    }
  } else if (tc.fProcess[eProcessSim]) {
    LOGF(fatal, "\033[1;31m%s at line %d : eProcessSim is not validated yet \033[0m", __FUNCTION__, __LINE__);
  }

  // **) Supported centrality estimators for Run 3 are enlisted here:
  if (tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim]) {
    if (!(ec.fsEventCuts[eCentralityEstimator].EqualTo("centFT0C", TString::kIgnoreCase) ||
          ec.fsEventCuts[eCentralityEstimator].EqualTo("centFT0CVariant1", TString::kIgnoreCase) ||
          ec.fsEventCuts[eCentralityEstimator].EqualTo("centFT0M", TString::kIgnoreCase) ||
          ec.fsEventCuts[eCentralityEstimator].EqualTo("centFV0A", TString::kIgnoreCase) ||
          ec.fsEventCuts[eCentralityEstimator].EqualTo("centNTPV", TString::kIgnoreCase) ||
          ec.fsEventCuts[eCentralityEstimator].EqualTo("centNGlobal", TString::kIgnoreCase))) {
      LOGF(fatal, "\033[1;31m%s at line %d : centrality estimator = %s is not supported yet for Run 3 analysis.\nUse \"centFT0C\", \"centFT0CVariant1\", \"centFT0M\", \"centFV0A\", \"centNTPV\", pr , \"centNGlobal\"\033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eCentralityEstimator].Data());
    }
  } else if (tc.fProcess[eProcessSim]) {
    LOGF(fatal, "\033[1;31m%s at line %d : eProcessSim is not validated yet \033[0m", __FUNCTION__, __LINE__);
  }

  // **) Supported occupancy estimators for Run 3 are enlisted here:
  if (tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim]) {
    if (!(ec.fsEventCuts[eOccupancyEstimator].EqualTo("TrackOccupancyInTimeRange", TString::kIgnoreCase) ||
          ec.fsEventCuts[eOccupancyEstimator].EqualTo("FT0COccupancyInTimeRange", TString::kIgnoreCase))) {
      LOGF(fatal, "\033[1;31m%s at line %d : occupancy estimator = %s is not supported yet for Run 3 analysis. \033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eOccupancyEstimator].Data());
    }
  }

  // **) Protection against particle cuts which are available, but not yet validated, or are meaningless, in Run 3:
  if (tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim]) {
    if (pc.fUseParticleCuts[etrackCutFlag]) {
      LOGF(fatal, "\033[1;31m%s at line %d : particle cut trackCutFlag is not validated in Run 3 as of 20250113 => it has no effect\033[0m", __FUNCTION__, __LINE__);
    }
    if (pc.fUseParticleCuts[eisQualityTrack]) {
      LOGF(fatal, "\033[1;31m%s at line %d : particle cut isQualityTrack is not validated in Run 3 as of 20250113 => it kills all reconstructed tracks \033[0m", __FUNCTION__, __LINE__);
    }
    if (pc.fUseParticleCuts[eisGlobalTrack]) {
      LOGF(fatal, "\033[1;31m%s at line %d : particle cut isGlobalTrack cannot be used in Run 3 => it kills all reconstructed tracks.\n To select global track in Run 3, use etrackCutFlagFb1 or etrackCutFlagFb2, see documentation in enum\033[0m", __FUNCTION__, __LINE__);
    }
  }

  // **) Protection on particle cuts which can be used only in Run 3:
  // trackCutFlag, trackCutFlagFb1, trackCutFlagFb2 => use only one at the time
  if (static_cast<int>(pc.fUseParticleCuts[etrackCutFlag]) + static_cast<int>(pc.fUseParticleCuts[etrackCutFlagFb1]) + static_cast<int>(pc.fUseParticleCuts[etrackCutFlagFb2]) >= 2) {
    LOGF(fatal, "\033[1;31m%s at line %d : use only one of trackCutFlag, trackCutFlagFb1, trackCutFlagFb2 at time. \033[0m", __FUNCTION__, __LINE__);
  }

  // isPVContributor:
  if (pc.fUseParticleCuts[eisPVContributor]) {
    if (!(tc.fProcess[eProcessRec] || tc.fProcess[eProcessRecSim] || tc.fProcess[eProcessSim])) {
      LOGF(fatal, "\033[1;31m%s at line %d : particle cut isPVContributor can be used only in Run 3\033[0m", __FUNCTION__, __LINE__);
    }
  }

  // ...

  // d) Insanity checks on binning:
  // ...

  // e) Insanity checks on events cuts:
  if (ec.fsEventCuts[eMultiplicityEstimator].EqualTo("ReferenceMultiplicity", TString::kIgnoreCase) && ec.fUseEventCuts[eMultiplicity]) {
    LOGF(fatal, "\033[1;31m%s at line %d : use ec.fUseEventCuts[eMultiplicity] only when eMultiplicityEstimator = SelectedTracks. Otherwise, things can happen... \033[0m", __FUNCTION__, __LINE__);
  }
  // ...

  // f) Insanity checks on Toy NUA:
  // ...

  // g) Insanity checks on internal validation:
  //    Remark: I check here only in the settings I could define in DefaultConfiguration().
  //            The other insanity checks are in BookInternalValidationHistograms() or in InsanityChecksAfterBooking()
  if (iv.fUseInternalValidation) {
    if (iv.fnEventsInternalValidation <= 0) {
      LOGF(fatal, "\033[1;31m%s at line %d : iv.fnEventsInternalValidation <= 0 => Set number of events to positive integer\033[0m", __FUNCTION__, __LINE__);
    }

    if (!(iv.fHarmonicsOptionInternalValidation->EqualTo("constant", TString::kIgnoreCase) ||
          iv.fHarmonicsOptionInternalValidation->EqualTo("correlated", TString::kIgnoreCase) ||
          iv.fHarmonicsOptionInternalValidation->EqualTo("persistent", TString::kIgnoreCase))) {
      LOGF(fatal, "\033[1;31m%s at line %d : fHarmonicsOptionInternalValidation = %s is not supported. \033[0m", __FUNCTION__, __LINE__, iv.fHarmonicsOptionInternalValidation->Data());
    }

    if (iv.fRescaleWithTheoreticalInput && (nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoops || nl.fCalculateKineCustomNestedLoops)) {
      LOGF(fatal, "\033[1;31m%s at line %d : rescaling with theoretical input is not supported when cross-check is done with nested loops. \033[0m", __FUNCTION__, __LINE__);
    }

    if (ec.fsEventCuts[eMultiplicityEstimator].EqualTo("ReferenceMultiplicity", TString::kIgnoreCase)) {
      LOGF(fatal, "\033[1;31m%s at line %d : in IV eMultiplicityEstimator cannot be set to \"ReferenceMultiplicity\" (yet) \033[0m", __FUNCTION__, __LINE__);
    }

  } // if (iv.fUseInternalValidation) {

  // h) Insanity checks on results histograms:
  //  **) Check if all arrays are initialized until the end:
  for (int afo = 0; afo < eAsFunctionOf_N; afo++) {
    if (res.fResultsProXaxisTitle[afo].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : res.fResultsProXaxisTitle[%d] is empty.\033[0m", __FUNCTION__, __LINE__, afo);
    }

    if (res.fResultsProRawName[afo].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : res.fResultsProRawName[%d] is empty.\033[0m", __FUNCTION__, __LINE__, afo);
    }
  } // for(int afo = 0; afo < eAsFunctionOf_N; afo++) {

  // ...

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void InsanityChecksBeforeBooking()

//============================================================

void InsanityChecksAfterBooking()
{
  // Do insanity checks on all booked histograms, etc.,
  // Configuration, binning and cuts are checked already before booking in InsanityChecksBeforeBooking().

  // a) Insanity checks on booking;
  // b) Insanity checks on internal validation;
  // ...

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Insanity checks on booking:

  // **) Check that the last bin is not empty in fBasePro, and that there is no underflow or overflow bins:
  if (TString(fBasePro->GetXaxis()->GetBinLabel(eConfiguration_N - 1)).EqualTo("")) {
    LOGF(fatal, "\033[1;31m%s at line %d : The very last bin of \"fBasePro\" doesn't have the title, check the booking of this hostogram. \033[0m", __FUNCTION__, __LINE__);
  }
  if (TMath::Abs(fBasePro->GetBinContent(0)) > 0.) {
    LOGF(fatal, "\033[1;31m%s at line %d : In \"fBasePro\" something was filled in the underflow, check the booking of this hostogram. \033[0m", __FUNCTION__, __LINE__);
  }
  if (TMath::Abs(fBasePro->GetBinContent(eConfiguration_N)) > 0.) {
    LOGF(fatal, "\033[1;31m%s at line %d : In \"fBasePro\" something was filled in the overflow, check the booking of this hostogram. \033[0m", __FUNCTION__, __LINE__);
  }

  // ...

  // b) Insanity checks on internal validation:
  if (iv.fUseInternalValidation) {

    // **) Check that rescaling is used only when it makes sense:
    if (iv.fRescaleWithTheoreticalInput && iv.fHarmonicsOptionInternalValidation->EqualTo("correlated")) {
      LOGF(fatal, "\033[1;31m%s at line %d : rescaling with theoretical input doesn't make sanse for fHarmonicsOptionInternalValidation = \"correlated\". \033[0m", __FUNCTION__, __LINE__);
    }
    if (iv.fRescaleWithTheoreticalInput && iv.fHarmonicsOptionInternalValidation->EqualTo("persistent")) {
      LOGF(fatal, "\033[1;31m%s at line %d : rescaling with theoretical input doesn't make sanse for fHarmonicsOptionInternalValidation = \"persistent\". \033[0m", __FUNCTION__, __LINE__);
    }

    // **) Print a warning if this histogram is not booked:
    if (!eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]) {
      LOGF(warning, "\033[1;31m%s at line %d : eh.fEventHistograms[eNumberOfEvents][eSim][eAfter] is not booked => no info on the total number of events in internal validation can be provided \033[0m", __FUNCTION__, __LINE__);
    }

  } // end of if (iv.fUseInternalValidation) {

  // ...

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void InsanityChecksAfterBooking()

//============================================================

bool Skip(int recOrSim)
{
  // Decide here whether a certain histogram, etc., will be booked and used both for eRec and eSim.
  // Same for cuts.

  if (tc.fVerboseUtility) {
    StartFunction(__FUNCTION__);
  }

  // *) Insanity check:
  if (!(recOrSim == eRec || recOrSim == eSim)) {
    LOGF(fatal, "\033[1;31m%s at line %d : recOrSim = %d \033[0m", __FUNCTION__, __LINE__, recOrSim);
  }

  // *) If I am doing internal validation, I book and fill only eSim:
  if (iv.fUseInternalValidation) {
    if (recOrSim == eRec) {
      return true; // yes, skip
    } else {
      return false; // this is eSim, do not skip
    }
  }

  // *) If I am analyzing only reconstructed data, do not book histos for simulated, and vice versa.
  //    TBI 20240223 tc.fProcess[eProcessTest] is treated as tc.fProcess[eProcessRec], for the time being
  if ((tc.fProcess[eGenericRec] && recOrSim == eSim) || (tc.fProcess[eGenericSim] && recOrSim == eRec)) {
    return true; // yes, skip
  }

  return false; // by default, I do not skip anything

} // bool Skip(int recOrSim)

//============================================================

void BookAndNestAllLists()
{
  // *) QA;
  //  **) QA event histograms;
  //  **) QA particle histograms;
  //  **) QA particle event histograms;
  // *) Control event histograms;
  // *) Control particle histograms;
  // *) Correlations;
  // *) Q-vectors;
  // *) Particle weights;
  // *) Centrality weights;
  // *) Nested loops;
  // *) Toy NUA;
  // *) Internal validation;
  // *) Test0;
  // *) Eta separations;
  // *) Results.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) QA:
  qa.fQAList = new TList();
  qa.fQAList->SetName("QA");
  qa.fQAList->SetOwner(true);
  fBaseList->Add(qa.fQAList);

  //  **) QA event histograms;
  if (qa.fFillQAEventHistograms2D) {
    qa.fQAEventList = new TList();
    qa.fQAEventList->SetName("QAEvent");
    qa.fQAEventList->SetOwner(true);
    qa.fQAList->Add(qa.fQAEventList); // yes, this one is nested within base QA TList
  }

  //  **) QA particle histograms;
  if (qa.fFillQAParticleHistograms2D) {
    qa.fQAParticleList = new TList();
    qa.fQAParticleList->SetName("QAParticle");
    qa.fQAParticleList->SetOwner(true);
    qa.fQAList->Add(qa.fQAParticleList); // yes, this one is nested within base QA TList
  }

  //  **) QA particle event histograms;
  if (qa.fFillQAParticleEventHistograms2D) {
    qa.fQAParticleEventList = new TList();
    qa.fQAParticleEventList->SetName("QAParticleEvent");
    qa.fQAParticleEventList->SetOwner(true);
    qa.fQAList->Add(qa.fQAParticleEventList); // yes, this one is nested within base QA TList
  }

  //  **) QA "correlations vs." histograms;
  if (qa.fFillQACorrelationsVsHistograms2D) {
    qa.fQACorrelationsVsList = new TList();
    qa.fQACorrelationsVsList->SetName("QACorrelationsVs");
    qa.fQACorrelationsVsList->SetOwner(true);
    qa.fQAList->Add(qa.fQACorrelationsVsList); // yes, this one is nested within base QA TList
  }

  // *) Event cuts:
  ec.fEventCutsList = new TList();
  ec.fEventCutsList->SetName("EventCuts");
  ec.fEventCutsList->SetOwner(true);
  fBaseList->Add(ec.fEventCutsList);

  // *) Control event histograms:
  eh.fEventHistogramsList = new TList();
  eh.fEventHistogramsList->SetName("EventHistograms");
  eh.fEventHistogramsList->SetOwner(true);
  fBaseList->Add(eh.fEventHistogramsList);

  // *) Particle cuts:
  pc.fParticleCutsList = new TList();
  pc.fParticleCutsList->SetName("ParticleCuts");
  pc.fParticleCutsList->SetOwner(true);
  fBaseList->Add(pc.fParticleCutsList);

  // *) Control particle histograms:
  ph.fParticleHistogramsList = new TList();
  ph.fParticleHistogramsList->SetName("ParticleHistograms");
  ph.fParticleHistogramsList->SetOwner(true);
  fBaseList->Add(ph.fParticleHistogramsList);

  // *) Q-vectors:
  qv.fQvectorList = new TList();
  qv.fQvectorList->SetName("Q-vectors");
  qv.fQvectorList->SetOwner(true);
  fBaseList->Add(qv.fQvectorList);

  // *) Correlations:
  mupa.fCorrelationsList = new TList();
  mupa.fCorrelationsList->SetName("Correlations");
  mupa.fCorrelationsList->SetOwner(true);
  fBaseList->Add(mupa.fCorrelationsList);

  // *) Particle weights:
  pw.fWeightsList = new TList();
  pw.fWeightsList->SetName("Weights");
  pw.fWeightsList->SetOwner(true);
  fBaseList->Add(pw.fWeightsList);

  // *) Centrality weights:
  cw.fCentralityWeightsList = new TList();
  cw.fCentralityWeightsList->SetName("CentralityWeights");
  cw.fCentralityWeightsList->SetOwner(true);
  fBaseList->Add(cw.fCentralityWeightsList);

  // *) Nested loops:
  nl.fNestedLoopsList = new TList();
  nl.fNestedLoopsList->SetName("NestedLoops");
  nl.fNestedLoopsList->SetOwner(true);
  fBaseList->Add(nl.fNestedLoopsList);

  // *) Toy NUA:
  nua.fNUAList = new TList();
  nua.fNUAList->SetName("ToyNUA");
  nua.fNUAList->SetOwner(true);
  fBaseList->Add(nua.fNUAList);

  // *) Internal validation:
  iv.fInternalValidationList = new TList();
  iv.fInternalValidationList->SetName("InternalValidation");
  iv.fInternalValidationList->SetOwner(true);
  fBaseList->Add(iv.fInternalValidationList);

  // *) Test0:
  t0.fTest0List = new TList();
  t0.fTest0List->SetName("Test0");
  t0.fTest0List->SetOwner(true);
  fBaseList->Add(t0.fTest0List);

  // *) Eta separations:
  es.fEtaSeparationsList = new TList();
  es.fEtaSeparationsList->SetName("EtaSeparations");
  es.fEtaSeparationsList->SetOwner(true);
  fBaseList->Add(es.fEtaSeparationsList);

  // *) Results:
  res.fResultsList = new TList();
  res.fResultsList->SetName("Results");
  res.fResultsList->SetOwner(true);
  fBaseList->Add(res.fResultsList);

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookAndNestAllLists()

//============================================================

void BookQAHistograms()
{
  // Book all QA histograms and other related objects.

  // TBI 20240520 There is a bit of code bloat here - I could introduce a new enum eEventParticle, and then use eEvent = 0 and eParticle = 1

  // a) Book the profile holding flags;
  // b) Common local variables;
  // c) Book specific QA 2D event histograms;
  // d) Book specific QA 2D particle histograms;
  // e) Book specific QA 2D particle event histograms;
  // f) Book specific QA 2D "correlations vs." histograms.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Print the warning message, because with too many 2D histograms with double precision, the code crashes in terminate, due to:
  /*
  [1450742:multiparticle-correlations-a-b]: [13:30:27][STATE] Exiting FairMQ state machine
  [1450742:multiparticle-correlations-a-b]: [13:30:27][FATAL] error while setting up workflow in o2-analysis-cf-multiparticle-correlations-ab: shmem: could not create a message of size 1282720912, alignment: 64, free memory: 1358639296
  [1450742:multiparticle-correlations-a-b]: terminate called after throwing an instance of 'o2::framework::RuntimeErrorRef'
  [1450742:multiparticle-correlations-a-b]: *** Program crashed (Aborted)
  [1450742:multiparticle-correlations-a-b]: Backtrace by DPL:
  */
  if (tc.fVerbose) {
    LOGF(info, "\033[1;33m%s: !!!! WARNING !!!! With too many 2D histograms with double precision, the code will crash in terminate (\"... shmem: could not create a message of size ...\") . Locally, you can circumvent this while testing by calling Bailout() explicitly. !!!! WARNING !!!! \033[0m", __FUNCTION__);
  }

  // a) Book the profile holding flags:
  qa.fQAHistogramsPro = new TProfile("fQAHistogramsPro", "flags for QA histograms", 6, 0., 6.);
  qa.fQAHistogramsPro->SetStats(false);
  qa.fQAHistogramsPro->SetLineColor(eColor);
  qa.fQAHistogramsPro->SetFillColor(eFillColor);
  qa.fQAHistogramsPro->GetXaxis()->SetBinLabel(1, "fCheckUnderflowAndOverflow");
  qa.fQAHistogramsPro->Fill(0.5, static_cast<int>(qa.fCheckUnderflowAndOverflow));
  qa.fQAHistogramsPro->GetXaxis()->SetBinLabel(2, "fFillQAEventHistograms2D");
  qa.fQAHistogramsPro->Fill(1.5, static_cast<int>(qa.fFillQAEventHistograms2D));
  qa.fQAHistogramsPro->GetXaxis()->SetBinLabel(3, "fFillQAParticleHistograms2D");
  qa.fQAHistogramsPro->Fill(2.5, static_cast<int>(qa.fFillQAParticleHistograms2D));
  qa.fQAHistogramsPro->GetXaxis()->SetBinLabel(4, "fFillQAParticleEventHistograms2D");
  qa.fQAHistogramsPro->Fill(3.5, static_cast<int>(qa.fFillQAParticleEventHistograms2D));
  qa.fQAHistogramsPro->GetXaxis()->SetBinLabel(5, "fFillQACorrelationsVsHistograms2D");
  qa.fQAHistogramsPro->Fill(4.5, static_cast<int>(qa.fFillQACorrelationsVsHistograms2D));
  qa.fQAHistogramsPro->GetXaxis()->SetBinLabel(6, "fRebin");
  qa.fQAHistogramsPro->Fill(5.5, static_cast<int>(qa.fRebin));

  // ...

  qa.fQAList->Add(qa.fQAHistogramsPro);

  // b) Common local variables:
  // ...

  // c) Book specific QA 2D event histograms:
  // Binning of 2D event histos: TBI 20240503 see if you can automate all this
  int nBins_x_Event[eQAEventHistograms2D_N] = {0};
  double min_x_Event[eQAEventHistograms2D_N] = {0.};
  double max_x_Event[eQAEventHistograms2D_N] = {0.};
  TString title_x_Event[eQAEventHistograms2D_N] = {""};
  int nBins_y_Event[eQAEventHistograms2D_N] = {0};
  double min_y_Event[eQAEventHistograms2D_N] = {0.};
  double max_y_Event[eQAEventHistograms2D_N] = {0.};
  TString title_y_Event[eQAEventHistograms2D_N] = {""};

  // *) "Multiplicity_vs_ReferenceMultiplicity":
  nBins_x_Event[eMultiplicity_vs_ReferenceMultiplicity] = static_cast<int>(eh.fEventHistogramsBins[eMultiplicity][0] / qa.fRebin);
  min_x_Event[eMultiplicity_vs_ReferenceMultiplicity] = eh.fEventHistogramsBins[eMultiplicity][1];
  max_x_Event[eMultiplicity_vs_ReferenceMultiplicity] = eh.fEventHistogramsBins[eMultiplicity][2];
  title_x_Event[eMultiplicity_vs_ReferenceMultiplicity] = FancyFormatting(eh.fEventHistogramsName[eMultiplicity].Data());
  nBins_y_Event[eMultiplicity_vs_ReferenceMultiplicity] = static_cast<int>(eh.fEventHistogramsBins[eReferenceMultiplicity][0] / qa.fRebin);
  min_y_Event[eMultiplicity_vs_ReferenceMultiplicity] = eh.fEventHistogramsBins[eReferenceMultiplicity][1];
  max_y_Event[eMultiplicity_vs_ReferenceMultiplicity] = eh.fEventHistogramsBins[eReferenceMultiplicity][2];
  title_y_Event[eMultiplicity_vs_ReferenceMultiplicity] = FancyFormatting(eh.fEventHistogramsName[eReferenceMultiplicity].Data());

  // *) "Multiplicity_vs_NContributors":
  nBins_x_Event[eMultiplicity_vs_NContributors] = static_cast<int>(eh.fEventHistogramsBins[eMultiplicity][0] / qa.fRebin);
  min_x_Event[eMultiplicity_vs_NContributors] = eh.fEventHistogramsBins[eMultiplicity][1];
  max_x_Event[eMultiplicity_vs_NContributors] = eh.fEventHistogramsBins[eMultiplicity][2];
  title_x_Event[eMultiplicity_vs_NContributors] = FancyFormatting(eh.fEventHistogramsName[eMultiplicity].Data());
  nBins_y_Event[eMultiplicity_vs_NContributors] = static_cast<int>(eh.fEventHistogramsBins[eNContributors][0] / qa.fRebin);
  min_y_Event[eMultiplicity_vs_NContributors] = eh.fEventHistogramsBins[eNContributors][1];
  max_y_Event[eMultiplicity_vs_NContributors] = eh.fEventHistogramsBins[eNContributors][2];
  title_y_Event[eMultiplicity_vs_NContributors] = FancyFormatting(eh.fEventHistogramsName[eNContributors].Data());

  // *) "Multiplicity_vs_Centrality":
  nBins_x_Event[eMultiplicity_vs_Centrality] = static_cast<int>(eh.fEventHistogramsBins[eMultiplicity][0] / qa.fRebin);
  min_x_Event[eMultiplicity_vs_Centrality] = eh.fEventHistogramsBins[eMultiplicity][1];
  max_x_Event[eMultiplicity_vs_Centrality] = eh.fEventHistogramsBins[eMultiplicity][2];
  title_x_Event[eMultiplicity_vs_Centrality] = FancyFormatting(eh.fEventHistogramsName[eMultiplicity].Data());
  nBins_y_Event[eMultiplicity_vs_Centrality] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eMultiplicity_vs_Centrality] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eMultiplicity_vs_Centrality] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eMultiplicity_vs_Centrality] = FancyFormatting(eh.fEventHistogramsName[eCentrality].Data());

  // *) "Multiplicity_vs_Vertex_z":
  nBins_x_Event[eMultiplicity_vs_Vertex_z] = static_cast<int>(eh.fEventHistogramsBins[eMultiplicity][0] / qa.fRebin);
  min_x_Event[eMultiplicity_vs_Vertex_z] = eh.fEventHistogramsBins[eMultiplicity][1];
  max_x_Event[eMultiplicity_vs_Vertex_z] = eh.fEventHistogramsBins[eMultiplicity][2];
  title_x_Event[eMultiplicity_vs_Vertex_z] = FancyFormatting(eh.fEventHistogramsName[eMultiplicity].Data());
  nBins_y_Event[eMultiplicity_vs_Vertex_z] = static_cast<int>(eh.fEventHistogramsBins[eVertex_z][0]);
  min_y_Event[eMultiplicity_vs_Vertex_z] = eh.fEventHistogramsBins[eVertex_z][1];
  max_y_Event[eMultiplicity_vs_Vertex_z] = eh.fEventHistogramsBins[eVertex_z][2];
  title_y_Event[eMultiplicity_vs_Vertex_z] = FancyFormatting(eh.fEventHistogramsName[eVertex_z].Data());

  // *) "Multiplicity_vs_Occupancy":
  nBins_x_Event[eMultiplicity_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eMultiplicity][0] / qa.fRebin);
  min_x_Event[eMultiplicity_vs_Occupancy] = eh.fEventHistogramsBins[eMultiplicity][1];
  max_x_Event[eMultiplicity_vs_Occupancy] = eh.fEventHistogramsBins[eMultiplicity][2];
  title_x_Event[eMultiplicity_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eMultiplicity].Data());
  nBins_y_Event[eMultiplicity_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eOccupancy][0] / qa.fRebin);
  min_y_Event[eMultiplicity_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][1];
  max_y_Event[eMultiplicity_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][2];
  title_y_Event[eMultiplicity_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eOccupancy].Data());

  // *) "ReferenceMultiplicity_vs_NContributors":
  nBins_x_Event[eReferenceMultiplicity_vs_NContributors] = static_cast<int>(eh.fEventHistogramsBins[eReferenceMultiplicity][0] / qa.fRebin);
  min_x_Event[eReferenceMultiplicity_vs_NContributors] = eh.fEventHistogramsBins[eReferenceMultiplicity][1];
  max_x_Event[eReferenceMultiplicity_vs_NContributors] = eh.fEventHistogramsBins[eReferenceMultiplicity][2];
  title_x_Event[eReferenceMultiplicity_vs_NContributors] = FancyFormatting(eh.fEventHistogramsName[eReferenceMultiplicity].Data());
  nBins_y_Event[eReferenceMultiplicity_vs_NContributors] = static_cast<int>(eh.fEventHistogramsBins[eNContributors][0] / qa.fRebin);
  min_y_Event[eReferenceMultiplicity_vs_NContributors] = eh.fEventHistogramsBins[eNContributors][1];
  max_y_Event[eReferenceMultiplicity_vs_NContributors] = eh.fEventHistogramsBins[eNContributors][2];
  title_y_Event[eReferenceMultiplicity_vs_NContributors] = FancyFormatting(eh.fEventHistogramsName[eNContributors].Data());

  // *) "ReferenceMultiplicity_vs_Centrality":
  nBins_x_Event[eReferenceMultiplicity_vs_Centrality] = static_cast<int>(eh.fEventHistogramsBins[eReferenceMultiplicity][0] / qa.fRebin);
  min_x_Event[eReferenceMultiplicity_vs_Centrality] = eh.fEventHistogramsBins[eReferenceMultiplicity][1];
  max_x_Event[eReferenceMultiplicity_vs_Centrality] = eh.fEventHistogramsBins[eReferenceMultiplicity][2];
  title_x_Event[eReferenceMultiplicity_vs_Centrality] = FancyFormatting(eh.fEventHistogramsName[eReferenceMultiplicity].Data());
  nBins_y_Event[eReferenceMultiplicity_vs_Centrality] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eReferenceMultiplicity_vs_Centrality] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eReferenceMultiplicity_vs_Centrality] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eReferenceMultiplicity_vs_Centrality] = FancyFormatting(eh.fEventHistogramsName[eCentrality].Data());

  // *) "ReferenceMultiplicity_vs_Vertex_z":
  nBins_x_Event[eReferenceMultiplicity_vs_Vertex_z] = static_cast<int>(eh.fEventHistogramsBins[eReferenceMultiplicity][0] / qa.fRebin);
  min_x_Event[eReferenceMultiplicity_vs_Vertex_z] = eh.fEventHistogramsBins[eReferenceMultiplicity][1];
  max_x_Event[eReferenceMultiplicity_vs_Vertex_z] = eh.fEventHistogramsBins[eReferenceMultiplicity][2];
  title_x_Event[eReferenceMultiplicity_vs_Vertex_z] = FancyFormatting(eh.fEventHistogramsName[eReferenceMultiplicity].Data());
  nBins_y_Event[eReferenceMultiplicity_vs_Vertex_z] = static_cast<int>(eh.fEventHistogramsBins[eVertex_z][0]);
  min_y_Event[eReferenceMultiplicity_vs_Vertex_z] = eh.fEventHistogramsBins[eVertex_z][1];
  max_y_Event[eReferenceMultiplicity_vs_Vertex_z] = eh.fEventHistogramsBins[eVertex_z][2];
  title_y_Event[eReferenceMultiplicity_vs_Vertex_z] = FancyFormatting(eh.fEventHistogramsName[eVertex_z].Data());

  // *) "ReferenceMultiplicity_vs_Occupancy":
  nBins_x_Event[eReferenceMultiplicity_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eReferenceMultiplicity][0] / qa.fRebin);
  min_x_Event[eReferenceMultiplicity_vs_Occupancy] = eh.fEventHistogramsBins[eReferenceMultiplicity][1];
  max_x_Event[eReferenceMultiplicity_vs_Occupancy] = eh.fEventHistogramsBins[eReferenceMultiplicity][2];
  title_x_Event[eReferenceMultiplicity_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eReferenceMultiplicity].Data());
  nBins_y_Event[eReferenceMultiplicity_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eOccupancy][0] / qa.fRebin);
  min_y_Event[eReferenceMultiplicity_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][1];
  max_y_Event[eReferenceMultiplicity_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][2];
  title_y_Event[eReferenceMultiplicity_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eOccupancy].Data());

  // *) "NContributors_vs_Centrality":
  nBins_x_Event[eNContributors_vs_Centrality] = static_cast<int>(eh.fEventHistogramsBins[eNContributors][0] / qa.fRebin);
  min_x_Event[eNContributors_vs_Centrality] = eh.fEventHistogramsBins[eNContributors][1];
  max_x_Event[eNContributors_vs_Centrality] = eh.fEventHistogramsBins[eNContributors][2];
  title_x_Event[eNContributors_vs_Centrality] = FancyFormatting(eh.fEventHistogramsName[eNContributors].Data());
  nBins_y_Event[eNContributors_vs_Centrality] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eNContributors_vs_Centrality] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eNContributors_vs_Centrality] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eNContributors_vs_Centrality] = FancyFormatting(eh.fEventHistogramsName[eCentrality].Data());

  // *) "NContributors_vs_Vertex_z":
  nBins_x_Event[eNContributors_vs_Vertex_z] = static_cast<int>(eh.fEventHistogramsBins[eNContributors][0] / qa.fRebin);
  min_x_Event[eNContributors_vs_Vertex_z] = eh.fEventHistogramsBins[eNContributors][1];
  max_x_Event[eNContributors_vs_Vertex_z] = eh.fEventHistogramsBins[eNContributors][2];
  title_x_Event[eNContributors_vs_Vertex_z] = FancyFormatting(eh.fEventHistogramsName[eNContributors].Data());
  nBins_y_Event[eNContributors_vs_Vertex_z] = static_cast<int>(eh.fEventHistogramsBins[eVertex_z][0]);
  min_y_Event[eNContributors_vs_Vertex_z] = eh.fEventHistogramsBins[eVertex_z][1];
  max_y_Event[eNContributors_vs_Vertex_z] = eh.fEventHistogramsBins[eVertex_z][2];
  title_y_Event[eNContributors_vs_Vertex_z] = FancyFormatting(eh.fEventHistogramsName[eVertex_z].Data());

  // *) "NContributors_vs_Occupancy":
  nBins_x_Event[eNContributors_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eNContributors][0] / qa.fRebin);
  min_x_Event[eNContributors_vs_Occupancy] = eh.fEventHistogramsBins[eNContributors][1];
  max_x_Event[eNContributors_vs_Occupancy] = eh.fEventHistogramsBins[eNContributors][2];
  title_x_Event[eNContributors_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eNContributors].Data());
  nBins_y_Event[eNContributors_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eOccupancy][0] / qa.fRebin);
  min_y_Event[eNContributors_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][1];
  max_y_Event[eNContributors_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][2];
  title_y_Event[eNContributors_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eOccupancy].Data());

  // *) "Centrality_vs_Vertex_z":
  nBins_x_Event[eCentrality_vs_Vertex_z] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentrality_vs_Vertex_z] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentrality_vs_Vertex_z] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentrality_vs_Vertex_z] = FancyFormatting(eh.fEventHistogramsName[eCentrality].Data());
  nBins_y_Event[eCentrality_vs_Vertex_z] = static_cast<int>(eh.fEventHistogramsBins[eVertex_z][0]);
  min_y_Event[eCentrality_vs_Vertex_z] = eh.fEventHistogramsBins[eVertex_z][1];
  max_y_Event[eCentrality_vs_Vertex_z] = eh.fEventHistogramsBins[eVertex_z][2];
  title_y_Event[eCentrality_vs_Vertex_z] = FancyFormatting(eh.fEventHistogramsName[eVertex_z].Data());

  // *) "Centrality_vs_Occupancy":
  nBins_x_Event[eCentrality_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentrality_vs_Occupancy] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentrality_vs_Occupancy] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentrality_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eCentrality].Data());
  nBins_y_Event[eCentrality_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eOccupancy][0] / qa.fRebin);
  min_y_Event[eCentrality_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][1];
  max_y_Event[eCentrality_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][2];
  title_y_Event[eCentrality_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eOccupancy].Data());

  // *) "Centrality_vs_ImpactParameter":
  nBins_x_Event[eCentrality_vs_ImpactParameter] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentrality_vs_ImpactParameter] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentrality_vs_ImpactParameter] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentrality_vs_ImpactParameter] = FancyFormatting(eh.fEventHistogramsName[eCentrality].Data());
  nBins_y_Event[eCentrality_vs_ImpactParameter] = static_cast<int>(eh.fEventHistogramsBins[eImpactParameter][0] / qa.fRebin);
  min_y_Event[eCentrality_vs_ImpactParameter] = eh.fEventHistogramsBins[eImpactParameter][1];
  max_y_Event[eCentrality_vs_ImpactParameter] = eh.fEventHistogramsBins[eImpactParameter][2];
  title_y_Event[eCentrality_vs_ImpactParameter] = FancyFormatting(eh.fEventHistogramsName[eImpactParameter].Data());

  // *) "Vertex_z_vs_Occupancy":
  nBins_x_Event[eVertex_z_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eVertex_z][0]);
  min_x_Event[eVertex_z_vs_Occupancy] = eh.fEventHistogramsBins[eVertex_z][1];
  max_x_Event[eVertex_z_vs_Occupancy] = eh.fEventHistogramsBins[eVertex_z][2];
  title_x_Event[eVertex_z_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eVertex_z].Data());
  nBins_y_Event[eVertex_z_vs_Occupancy] = static_cast<int>(eh.fEventHistogramsBins[eOccupancy][0] / qa.fRebin);
  min_y_Event[eVertex_z_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][1];
  max_y_Event[eVertex_z_vs_Occupancy] = eh.fEventHistogramsBins[eOccupancy][2];
  title_y_Event[eVertex_z_vs_Occupancy] = FancyFormatting(eh.fEventHistogramsName[eOccupancy].Data());

  // *) "eMultNTracksPV_vs_MultNTracksGlobal":
  nBins_x_Event[eMultNTracksPV_vs_MultNTracksGlobal] = static_cast<int>(eh.fEventHistogramsBins[eMultiplicity][0] / qa.fRebin);
  min_x_Event[eMultNTracksPV_vs_MultNTracksGlobal] = eh.fEventHistogramsBins[eMultiplicity][1];
  max_x_Event[eMultNTracksPV_vs_MultNTracksGlobal] = eh.fEventHistogramsBins[eMultiplicity][2];
  title_x_Event[eMultNTracksPV_vs_MultNTracksGlobal] = FancyFormatting(qa.fReferenceMultiplicityEstimatorName[eMultNTracksPV].Data());
  nBins_y_Event[eMultNTracksPV_vs_MultNTracksGlobal] = static_cast<int>(eh.fEventHistogramsBins[eMultiplicity][0] / qa.fRebin);
  min_y_Event[eMultNTracksPV_vs_MultNTracksGlobal] = eh.fEventHistogramsBins[eMultiplicity][1];
  max_y_Event[eMultNTracksPV_vs_MultNTracksGlobal] = eh.fEventHistogramsBins[eMultiplicity][2];
  title_y_Event[eMultNTracksPV_vs_MultNTracksGlobal] = FancyFormatting(qa.fReferenceMultiplicityEstimatorName[eMultNTracksGlobal].Data());

  // *) "eCentFT0C_vs_CentFT0CVariant1":
  nBins_x_Event[eCentFT0C_vs_CentFT0CVariant1] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentFT0C_vs_CentFT0CVariant1] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentFT0C_vs_CentFT0CVariant1] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentFT0C_vs_CentFT0CVariant1] = FancyFormatting(qa.fCentralityEstimatorName[eCentFT0C].Data());
  nBins_y_Event[eCentFT0C_vs_CentFT0CVariant1] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eCentFT0C_vs_CentFT0CVariant1] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eCentFT0C_vs_CentFT0CVariant1] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eCentFT0C_vs_CentFT0CVariant1] = FancyFormatting(qa.fCentralityEstimatorName[eCentFT0CVariant1].Data());

  // *) "eCentFT0C_vs_CentFT0M":
  nBins_x_Event[eCentFT0C_vs_CentFT0M] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentFT0C_vs_CentFT0M] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentFT0C_vs_CentFT0M] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentFT0C_vs_CentFT0M] = FancyFormatting(qa.fCentralityEstimatorName[eCentFT0C].Data());
  nBins_y_Event[eCentFT0C_vs_CentFT0M] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eCentFT0C_vs_CentFT0M] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eCentFT0C_vs_CentFT0M] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eCentFT0C_vs_CentFT0M] = FancyFormatting(qa.fCentralityEstimatorName[eCentFT0M].Data());

  // *) "eCentFT0C_vs_CentFV0A":
  nBins_x_Event[eCentFT0C_vs_CentFV0A] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentFT0C_vs_CentFV0A] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentFT0C_vs_CentFV0A] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentFT0C_vs_CentFV0A] = FancyFormatting(qa.fCentralityEstimatorName[eCentFT0C].Data());
  nBins_y_Event[eCentFT0C_vs_CentFV0A] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eCentFT0C_vs_CentFV0A] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eCentFT0C_vs_CentFV0A] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eCentFT0C_vs_CentFV0A] = FancyFormatting(qa.fCentralityEstimatorName[eCentFV0A].Data());

  // *) "eCentFT0C_vs_CentNTPV":
  nBins_x_Event[eCentFT0C_vs_CentNTPV] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentFT0C_vs_CentNTPV] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentFT0C_vs_CentNTPV] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentFT0C_vs_CentNTPV] = FancyFormatting(qa.fCentralityEstimatorName[eCentFT0C].Data());
  nBins_y_Event[eCentFT0C_vs_CentNTPV] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eCentFT0C_vs_CentNTPV] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eCentFT0C_vs_CentNTPV] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eCentFT0C_vs_CentNTPV] = FancyFormatting(qa.fCentralityEstimatorName[eCentNTPV].Data());

  // *) "eCentFT0C_vs_CentNGlobal":
  nBins_x_Event[eCentFT0C_vs_CentNGlobal] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentFT0C_vs_CentNGlobal] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentFT0C_vs_CentNGlobal] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentFT0C_vs_CentNGlobal] = FancyFormatting(qa.fCentralityEstimatorName[eCentFT0C].Data());
  nBins_y_Event[eCentFT0C_vs_CentNGlobal] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eCentFT0C_vs_CentNGlobal] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eCentFT0C_vs_CentNGlobal] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eCentFT0C_vs_CentNGlobal] = FancyFormatting(qa.fCentralityEstimatorName[eCentNGlobal].Data());

  // *) "eCentFT0M_vs_CentNTPV":
  nBins_x_Event[eCentFT0M_vs_CentNTPV] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentFT0M_vs_CentNTPV] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentFT0M_vs_CentNTPV] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentFT0M_vs_CentNTPV] = FancyFormatting(qa.fCentralityEstimatorName[eCentFT0M].Data());
  nBins_y_Event[eCentFT0M_vs_CentNTPV] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eCentFT0M_vs_CentNTPV] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eCentFT0M_vs_CentNTPV] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eCentFT0M_vs_CentNTPV] = FancyFormatting(qa.fCentralityEstimatorName[eCentNTPV].Data());

  // *) "eCentRun2V0M_vs_CentRun2SPDTracklets":
  nBins_x_Event[eCentRun2V0M_vs_CentRun2SPDTracklets] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_x_Event[eCentRun2V0M_vs_CentRun2SPDTracklets] = eh.fEventHistogramsBins[eCentrality][1];
  max_x_Event[eCentRun2V0M_vs_CentRun2SPDTracklets] = eh.fEventHistogramsBins[eCentrality][2];
  title_x_Event[eCentRun2V0M_vs_CentRun2SPDTracklets] = FancyFormatting(qa.fCentralityEstimatorName[eCentRun2V0M].Data());
  nBins_y_Event[eCentRun2V0M_vs_CentRun2SPDTracklets] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_Event[eCentRun2V0M_vs_CentRun2SPDTracklets] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_Event[eCentRun2V0M_vs_CentRun2SPDTracklets] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_Event[eCentRun2V0M_vs_CentRun2SPDTracklets] = FancyFormatting(qa.fCentralityEstimatorName[eCentRun2SPDTracklets].Data());

  // *) "eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange":
  nBins_x_Event[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = static_cast<int>(eh.fEventHistogramsBins[eOccupancy][0] / qa.fRebin);
  min_x_Event[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = eh.fEventHistogramsBins[eOccupancy][1];
  max_x_Event[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = eh.fEventHistogramsBins[eOccupancy][2];
  title_x_Event[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = FancyFormatting(qa.fOccupancyEstimatorName[eTrackOccupancyInTimeRange].Data());
  nBins_y_Event[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = static_cast<int>(eh.fEventHistogramsBins[eOccupancy][0] / qa.fRebin);
  min_y_Event[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = eh.fEventHistogramsBins[eOccupancy][1];
  max_y_Event[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = eh.fEventHistogramsBins[eOccupancy][2];
  title_y_Event[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange] = FancyFormatting(qa.fOccupancyEstimatorName[eFT0COccupancyInTimeRange].Data());

  // *) "eCurrentRunDuration_vs_InteractionRate":
  nBins_x_Event[eCurrentRunDuration_vs_InteractionRate] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0] / qa.fRebin);
  min_x_Event[eCurrentRunDuration_vs_InteractionRate] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_Event[eCurrentRunDuration_vs_InteractionRate] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_Event[eCurrentRunDuration_vs_InteractionRate] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_Event[eCurrentRunDuration_vs_InteractionRate] = static_cast<int>(eh.fEventHistogramsBins[eInteractionRate][0] / qa.fRebin);
  min_y_Event[eCurrentRunDuration_vs_InteractionRate] = eh.fEventHistogramsBins[eInteractionRate][1];
  max_y_Event[eCurrentRunDuration_vs_InteractionRate] = eh.fEventHistogramsBins[eInteractionRate][2];
  title_y_Event[eCurrentRunDuration_vs_InteractionRate] = FancyFormatting(eh.fEventHistogramsName[eInteractionRate].Data());

  // ...

  // *) Quick insanity check on title_x_Event and title_y_Event:
  for (int t = 0; t < eQAEventHistograms2D_N; t++) {

    // **) title_x_Event:
    if (tc.fVerbose) {
      LOGF(info, "\033[1;32m title_x_Event[%d] = %s \033[0m", t, title_x_Event[t].Data());
    }
    if (title_x_Event[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : title_x_Event[%d] is not set, check corresponding enum \033[0m", __FUNCTION__, __LINE__, t);
    }

    // **) title_y_Event:
    if (tc.fVerbose) {
      LOGF(info, "\033[1;32m title_y_Event[%d] = %s \033[0m", t, title_y_Event[t].Data());
    }
    if (title_y_Event[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : title_y_Event[%d] is not set, check corresponding enum  \033[0m", __FUNCTION__, __LINE__, t);
    }

  } // for (int t = 0; t < eQAEventHistograms2D_N; t++) {

  // Okay, let's book 'em all:
  for (int t = 0; t < eQAEventHistograms2D_N; t++) // type, see enum eQAEventHistograms2D
  {
    if (!qa.fBookQAEventHistograms2D[t]) {
      continue;
    }
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {

      if (Skip(rs)) {
        continue;
      }

      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {

        // Special treatment for eMultiplicity => I will never fill this one before the cuts, if Multiplicity = SelectedTracks, obviously:
        if (ba == eBefore && (title_x_Event[t].BeginsWith("Multiplicity") || title_y_Event[t].BeginsWith("Multiplicity")) && ec.fsEventCuts[eMultiplicityEstimator].EqualTo("SelectedTracks", TString::kIgnoreCase)) {
          // TBI 20241123 what remains ill-defined is the case when Multiplicity != SelectedTracks , check that further
          // TBI 20241123 not sure if checking with BeginsWith(...) x2 is robust enough
          // TBI 20241123 just like I have Skip(rs), introduce the same thingie for "ba" counter + propagate to other member functions
          // TBI 20241124 there is a corner case when eMultiplicityEstimator itself is "ReferenceMultiplicity" => all 2D QA booked both before and after cuts,
          //              but it's filled trivially before the cuts, because Multiplicity is always 0. Re-think this at some point.
          continue;
        }

        qa.fQAEventHistograms2D[t][rs][ba] = new TH2F(
          TString::Format("fQAEventHistograms2D[%s][%s][%s]", qa.fEventHistogramsName2D[t].Data(), gc.srs[rs].Data(), gc.sba[ba].Data()),
          TString::Format("%s, %s, %s", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.sba_long[ba].Data()), // __RUN_NUMBER__ is handled in PropagateRunNumber(...)
          nBins_x_Event[t], min_x_Event[t], max_x_Event[t], nBins_y_Event[t], min_y_Event[t], max_y_Event[t]);
        qa.fQAEventHistograms2D[t][rs][ba]->GetXaxis()->SetTitle(title_x_Event[t].Data());
        qa.fQAEventHistograms2D[t][rs][ba]->GetYaxis()->SetTitle(title_y_Event[t].Data());
        qa.fQAEventHistograms2D[t][rs][ba]->SetLineColor(ec.fBeforeAfterColor[ba]);
        qa.fQAEventHistograms2D[t][rs][ba]->SetFillColor(ec.fBeforeAfterColor[ba] - 10);
        qa.fQAEventHistograms2D[t][rs][ba]->SetOption("col");
        qa.fQAEventList->Add(qa.fQAEventHistograms2D[t][rs][ba]);
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for(int t=0;t<eQAEventHistograms2D_N;t++) // type, see enum eEventHistograms2D

  // ...

  // c) Book specific QA 2D particle histograms:
  // Binning of 2D particle histos: TBI 20240503 see if you can automate all this
  int nBins_x_Particle[eQAParticleHistograms2D_N] = {0};
  double min_x_Particle[eQAParticleHistograms2D_N] = {0.};
  double max_x_Particle[eQAParticleHistograms2D_N] = {0.};
  TString title_x_Particle[eQAParticleHistograms2D_N] = {""};
  int nBins_y_Particle[eQAParticleHistograms2D_N] = {0};
  double min_y_Particle[eQAParticleHistograms2D_N] = {0.};
  double max_y_Particle[eQAParticleHistograms2D_N] = {0.};
  TString title_y_Particle[eQAParticleHistograms2D_N] = {""};

  // *) "pt_vs_dcaXY":
  nBins_x_Particle[ePt_vs_dcaXY] = static_cast<int>(ph.fParticleHistogramsBins[ePt][0]); // TBI 20240702 add support for rebinning
  min_x_Particle[ePt_vs_dcaXY] = ph.fParticleHistogramsBins[ePt][1];
  max_x_Particle[ePt_vs_dcaXY] = ph.fParticleHistogramsBins[ePt][2];
  title_x_Particle[ePt_vs_dcaXY] = FancyFormatting(ph.fParticleHistogramsName[ePt].Data());
  nBins_y_Particle[ePt_vs_dcaXY] = static_cast<int>(ph.fParticleHistogramsBins[edcaXY][0]); // TBI 20240702 add support for rebinning
  min_y_Particle[ePt_vs_dcaXY] = ph.fParticleHistogramsBins[edcaXY][1];
  max_y_Particle[ePt_vs_dcaXY] = ph.fParticleHistogramsBins[edcaXY][2];
  title_y_Particle[ePt_vs_dcaXY] = FancyFormatting(ph.fParticleHistogramsName[edcaXY].Data());

  // ...

  // *) Quick insanity check on title_x_Particle and title_y_Particle:
  for (int t = 0; t < eQAParticleHistograms2D_N; t++) {
    if (title_x_Particle[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : title_x_Particle[%d] is not set, check corresponding enum \033[0m", __FUNCTION__, __LINE__, t);
    }
    if (title_y_Particle[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : title_y_Particle[%d] is not set, check corresponding enum  \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // Okay, let's book 'em all:
  for (int t = 0; t < eQAParticleHistograms2D_N; t++) // type, see enum eQAParticleHistograms2D
  {
    if (!qa.fBookQAParticleHistograms2D[t]) {
      continue;
    }
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {

      if (Skip(rs)) {
        continue;
      }

      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        qa.fQAParticleHistograms2D[t][rs][ba] = new TH2F(
          TString::Format("fQAParticleHistograms2D[%s][%s][%s]", qa.fParticleHistogramsName2D[t].Data(), gc.srs[rs].Data(), gc.sba[ba].Data()),
          TString::Format("%s, %s, %s", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.sba_long[ba].Data()), // __RUN_NUMBER__ is handled in PropagateRunNumber(...)
          nBins_x_Particle[t], min_x_Particle[t], max_x_Particle[t], nBins_y_Particle[t], min_y_Particle[t], max_y_Particle[t]);

        qa.fQAParticleHistograms2D[t][rs][ba]->GetXaxis()->SetTitle(title_x_Particle[t].Data());
        qa.fQAParticleHistograms2D[t][rs][ba]->GetYaxis()->SetTitle(title_y_Particle[t].Data());
        qa.fQAParticleHistograms2D[t][rs][ba]->SetLineColor(ec.fBeforeAfterColor[ba]);
        qa.fQAParticleHistograms2D[t][rs][ba]->SetFillColor(ec.fBeforeAfterColor[ba] - 10);
        qa.fQAParticleHistograms2D[t][rs][ba]->SetOption("col");
        qa.fQAParticleList->Add(qa.fQAParticleHistograms2D[t][rs][ba]);
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for(int t=0;t<eQAParticleHistograms2D_N;t++) // type, see enum eParticleHistograms2D

  // e) Book specific QA 2D particle event histograms:

  // Binning of 2D particle event histos:
  int nBins_x_ParticleEvent[eQAParticleEventHistograms2D_N] = {0};
  double min_x_ParticleEvent[eQAParticleEventHistograms2D_N] = {0.};
  double max_x_ParticleEvent[eQAParticleEventHistograms2D_N] = {0.};
  TString title_x_ParticleEvent[eQAParticleEventHistograms2D_N] = {""};
  int nBins_y_ParticleEvent[eQAParticleEventHistograms2D_N] = {0};
  double min_y_ParticleEvent[eQAParticleEventHistograms2D_N] = {0.};
  double max_y_ParticleEvent[eQAParticleEventHistograms2D_N] = {0.};
  TString title_y_ParticleEvent[eQAParticleEventHistograms2D_N] = {""};

  // *) "eCurrentRunDuration_vs_itsNClsEbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsEbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsEbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsEbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsEbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  // nBins_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsEbyE] = static_cast<int>(ph.fParticleHistogramsBins[eitsNCls][0]);
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsEbyE] = 100;
  min_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsEbyE] = ph.fParticleHistogramsBins[eitsNCls][1];
  max_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsEbyE] = ph.fParticleHistogramsBins[eitsNCls][2];
  title_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsEbyE] = TString::Format("#LT%s#GT", FancyFormatting(ph.fParticleHistogramsName[eitsNCls].Data()));

  // *) "eCurrentRunDuration_vs_itsNClsNegEtaEbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = 100;
  min_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = ph.fParticleHistogramsBins[eitsNCls][1];
  max_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = ph.fParticleHistogramsBins[eitsNCls][2];
  title_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsNegEtaEbyE] = TString::Format("#LT%s#GT, #eta < 0", FancyFormatting(ph.fParticleHistogramsName[eitsNCls].Data()));

  // *) "eCurrentRunDuration_vs_itsNClsPosEtaEbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = 100;
  min_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = ph.fParticleHistogramsBins[eitsNCls][1];
  max_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = ph.fParticleHistogramsBins[eitsNCls][2];
  title_y_ParticleEvent[eCurrentRunDuration_vs_itsNClsPosEtaEbyE] = TString::Format("#LT%s#GT, #eta > 0", FancyFormatting(ph.fParticleHistogramsName[eitsNCls].Data()));

  // *) "eCurrentRunDuration_vs_Eta0804EbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_Eta0804EbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_Eta0804EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_Eta0804EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_Eta0804EbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_Eta0804EbyE] = 80;
  min_y_ParticleEvent[eCurrentRunDuration_vs_Eta0804EbyE] = -1.0; // TBI 20241214 intentionally temporarily overshooting, to trace down overflow and underflow, if any
  max_y_ParticleEvent[eCurrentRunDuration_vs_Eta0804EbyE] = 1.0;  // TBI 20241214 intentionally temporarily overshooting, to trace down overflow and underflow, if any
  title_y_ParticleEvent[eCurrentRunDuration_vs_Eta0804EbyE] = TString::Format("#LT%s#GT, -0.8 < #eta < -0.4", FancyFormatting(ph.fParticleHistogramsName[eEta].Data()));

  // *) "eCurrentRunDuration_vs_Eta0400EbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_Eta0400EbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_Eta0400EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_Eta0400EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_Eta0400EbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_Eta0400EbyE] = 80;
  min_y_ParticleEvent[eCurrentRunDuration_vs_Eta0400EbyE] = -1.0;
  max_y_ParticleEvent[eCurrentRunDuration_vs_Eta0400EbyE] = 1.0;
  title_y_ParticleEvent[eCurrentRunDuration_vs_Eta0400EbyE] = TString::Format("#LT%s#GT, -0.4 < #eta < 0.0", FancyFormatting(ph.fParticleHistogramsName[eEta].Data()));

  // *) "eCurrentRunDuration_vs_Eta0004EbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_Eta0004EbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_Eta0004EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_Eta0004EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_Eta0004EbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_Eta0004EbyE] = 80;
  min_y_ParticleEvent[eCurrentRunDuration_vs_Eta0004EbyE] = -1.0;
  max_y_ParticleEvent[eCurrentRunDuration_vs_Eta0004EbyE] = 1.0;
  title_y_ParticleEvent[eCurrentRunDuration_vs_Eta0004EbyE] = TString::Format("#LT%s#GT, 0.0 < #eta < 0.4", FancyFormatting(ph.fParticleHistogramsName[eEta].Data()));

  // *) "eCurrentRunDuration_vs_Eta0408EbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_Eta0408EbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_Eta0408EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_Eta0408EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_Eta0408EbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_Eta0408EbyE] = 80;
  min_y_ParticleEvent[eCurrentRunDuration_vs_Eta0408EbyE] = -1.0;
  max_y_ParticleEvent[eCurrentRunDuration_vs_Eta0408EbyE] = 1.0;
  title_y_ParticleEvent[eCurrentRunDuration_vs_Eta0408EbyE] = TString::Format("#LT%s#GT, 0.4 < #eta < 0.8", FancyFormatting(ph.fParticleHistogramsName[eEta].Data()));

  // *) "eCurrentRunDuration_vs_Pt0005EbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_Pt0005EbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_Pt0005EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_Pt0005EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_Pt0005EbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_Pt0005EbyE] = 400;
  min_y_ParticleEvent[eCurrentRunDuration_vs_Pt0005EbyE] = 0.0;  // TBI 20241214 intentionally temporarilyovershooting, to trace down overflow and underflow, if any
  max_y_ParticleEvent[eCurrentRunDuration_vs_Pt0005EbyE] = 10.0; // TBI 20241214 intentionally temporarily overshooting, to trace down overflow and underflow, if any
  title_y_ParticleEvent[eCurrentRunDuration_vs_Pt0005EbyE] = TString::Format("#LT%s#GT, 0.0 < p_{T} < 0.5 GeV/c", FancyFormatting(ph.fParticleHistogramsName[ePt].Data()));

  // *) "eCurrentRunDuration_vs_Pt0510EbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_Pt0510EbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_Pt0510EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_Pt0510EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_Pt0510EbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_Pt0510EbyE] = 400;
  min_y_ParticleEvent[eCurrentRunDuration_vs_Pt0510EbyE] = 0.0;
  max_y_ParticleEvent[eCurrentRunDuration_vs_Pt0510EbyE] = 10.0;
  title_y_ParticleEvent[eCurrentRunDuration_vs_Pt0510EbyE] = TString::Format("#LT%s#GT, 0.5 < p_{T} < 1.0 GeV/c", FancyFormatting(ph.fParticleHistogramsName[ePt].Data()));

  // *) "eCurrentRunDuration_vs_Pt1050EbyE":
  nBins_x_ParticleEvent[eCurrentRunDuration_vs_Pt1050EbyE] = static_cast<int>(eh.fEventHistogramsBins[eCurrentRunDuration][0]);
  min_x_ParticleEvent[eCurrentRunDuration_vs_Pt1050EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][1];
  max_x_ParticleEvent[eCurrentRunDuration_vs_Pt1050EbyE] = eh.fEventHistogramsBins[eCurrentRunDuration][2];
  title_x_ParticleEvent[eCurrentRunDuration_vs_Pt1050EbyE] = FancyFormatting(eh.fEventHistogramsName[eCurrentRunDuration].Data());
  nBins_y_ParticleEvent[eCurrentRunDuration_vs_Pt1050EbyE] = 400;
  min_y_ParticleEvent[eCurrentRunDuration_vs_Pt1050EbyE] = 0.0;
  max_y_ParticleEvent[eCurrentRunDuration_vs_Pt1050EbyE] = 10.0;
  title_y_ParticleEvent[eCurrentRunDuration_vs_Pt1050EbyE] = TString::Format("#LT%s#GT, 1.0 < p_{T} < 5.0 GeV/c", FancyFormatting(ph.fParticleHistogramsName[ePt].Data()));

  // ...

  // *) Quick insanity check on title_x_ParticleEvent and title_y_ParticleEvent:
  for (int t = 0; t < eQAParticleEventHistograms2D_N; t++) {
    if (title_x_ParticleEvent[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : title_x_ParticleEvent[%d] is not set, check corresponding enum \033[0m", __FUNCTION__, __LINE__, t);
    }
    if (title_y_ParticleEvent[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : title_y_ParticleEvent[%d] is not set, check corresponding enum  \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // Okay, let's book 'em all:
  for (int t = 0; t < eQAParticleEventHistograms2D_N; t++) // type, see enum eQAParticleEventHistograms2D
  {
    if (!qa.fBookQAParticleEventHistograms2D[t]) {
      continue;
    }
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {

      if (Skip(rs)) {
        continue;
      }

      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {

        if (ba == eBefore) { // TBI 20241214 re-think if I need these additional QA particle event histos before cuts
          continue;
        }

        qa.fQAParticleEventHistograms2D[t][rs][ba] = new TH2F(
          TString::Format("fQAParticleEventHistograms2D[%s][%s][%s]", qa.fQAParticleEventHistogramsName2D[t].Data(), gc.srs[rs].Data(), gc.sba[ba].Data()),
          TString::Format("%s, %s, %s", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.sba_long[ba].Data()), // __RUN_NUMBER__ is handled in PropagateRunNumber(...)
          nBins_x_ParticleEvent[t], min_x_ParticleEvent[t], max_x_ParticleEvent[t], nBins_y_ParticleEvent[t], min_y_ParticleEvent[t], max_y_ParticleEvent[t]);

        qa.fQAParticleEventHistograms2D[t][rs][ba]->GetXaxis()->SetTitle(title_x_ParticleEvent[t].Data());
        qa.fQAParticleEventHistograms2D[t][rs][ba]->GetYaxis()->SetTitle(title_y_ParticleEvent[t].Data());
        qa.fQAParticleEventHistograms2D[t][rs][ba]->SetLineColor(ec.fBeforeAfterColor[ba]);
        qa.fQAParticleEventHistograms2D[t][rs][ba]->SetFillColor(ec.fBeforeAfterColor[ba] - 10);
        qa.fQAParticleEventHistograms2D[t][rs][ba]->SetOption("col");
        qa.fQAParticleEventList->Add(qa.fQAParticleEventHistograms2D[t][rs][ba]);
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for(int t=0;t<eQAParticleEventHistograms2D_N;t++) // type, see enum eParticleEventHistograms2D

  // Helper profile ...
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {

    if (Skip(rs)) {
      continue;
    }

    for (int ba = 0; ba < 2; ba++) // before/after cuts
    {

      if (ba == eBefore) { // TBI 20241214 re-think if I need these additional QA particle event histos before cuts
        continue;
      }

      qa.fQAParticleEventProEbyE[rs][ba] = new TProfile(
        TString::Format("fParticleEventProEbyE[%s][%s]", gc.srs[rs].Data(), gc.sba[ba].Data()),
        TString::Format("%s, %s", gc.srs_long[rs].Data(), gc.sba_long[ba].Data()),
        eQAParticleEventProEbyE_N, 0., eQAParticleEventProEbyE_N);
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(eitsNClsEbyE, "#LTitsNCls#GT"); // TBI 20241214 this bin labeling is not really needed, as I never save this TProfile persistently
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(eitsNClsNegEtaEbyE, "#LTitsNClsNegEta#GT");
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(eitsNClsPosEtaEbyE, "#LTitsNClsPosEta#GT");
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(eEta0804EbyE, "#LTEta0804EbyE#GT");
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(eEta0400EbyE, "#LTEta0400EbyE#GT");
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(eEta0004EbyE, "#LTEta0004EbyE#GT");
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(eEta0408EbyE, "#LTEta0408EbyE#GT");
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(ePt0005EbyE, "#LTPt0005EbyE#GT");
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(ePt0510EbyE, "#LTPt0510EbyE#GT");
      qa.fQAParticleEventProEbyE[rs][ba]->GetXaxis()->SetBinLabel(ePt1050EbyE, "#LTPt1050EbyE#GT");
    }
  }

  // f) Book specific QA 2D "correlations vs." histograms:

  // Binning of 2D "correlations vs." histos:
  // Remark: I use the same binning for x axis for all 2D QA histos in this category, therefore here implementation in shorter.
  int nBins_x_CorrelationsVs = 2000;
  double min_x_CorrelationsVs = -1.;
  double max_x_CorrelationsVs = 1.;
  TString title_x_CorrelationsVs = "#LT2#GT"; // harmonic I store elsewhere TBI-today document here where
  int nBins_y_CorrelationsVs[eQACorrelationsVsHistograms2D_N] = {0};
  double min_y_CorrelationsVs[eQACorrelationsVsHistograms2D_N] = {0.};
  double max_y_CorrelationsVs[eQACorrelationsVsHistograms2D_N] = {0.};
  TString title_y_CorrelationsVs[eQACorrelationsVsHistograms2D_N] = {""};

  // *) "eCorrelations_vs_Multiplicity":
  nBins_y_CorrelationsVs[eCorrelations_vs_Multiplicity] = static_cast<int>(eh.fEventHistogramsBins[eMultiplicity][0] / qa.fRebin);
  min_y_CorrelationsVs[eCorrelations_vs_Multiplicity] = eh.fEventHistogramsBins[eMultiplicity][1];
  max_y_CorrelationsVs[eCorrelations_vs_Multiplicity] = eh.fEventHistogramsBins[eMultiplicity][2];
  title_y_CorrelationsVs[eCorrelations_vs_Multiplicity] = FancyFormatting(eh.fEventHistogramsName[eMultiplicity].Data());

  // *) "eCorrelations_vs_ReferenceMultiplicity":
  nBins_y_CorrelationsVs[eCorrelations_vs_ReferenceMultiplicity] = static_cast<int>(eh.fEventHistogramsBins[eReferenceMultiplicity][0] / qa.fRebin);
  min_y_CorrelationsVs[eCorrelations_vs_ReferenceMultiplicity] = eh.fEventHistogramsBins[eReferenceMultiplicity][1];
  max_y_CorrelationsVs[eCorrelations_vs_ReferenceMultiplicity] = eh.fEventHistogramsBins[eReferenceMultiplicity][2];
  title_y_CorrelationsVs[eCorrelations_vs_ReferenceMultiplicity] = FancyFormatting(eh.fEventHistogramsName[eReferenceMultiplicity].Data());

  // *) "eCorrelations_vs_Centrality":
  nBins_y_CorrelationsVs[eCorrelations_vs_Centrality] = static_cast<int>(eh.fEventHistogramsBins[eCentrality][0]);
  min_y_CorrelationsVs[eCorrelations_vs_Centrality] = eh.fEventHistogramsBins[eCentrality][1];
  max_y_CorrelationsVs[eCorrelations_vs_Centrality] = eh.fEventHistogramsBins[eCentrality][2];
  title_y_CorrelationsVs[eCorrelations_vs_Centrality] = FancyFormatting(eh.fEventHistogramsName[eCentrality].Data());

  // .....

  // *) "eCorrelations_vs_MeanPhi":
  nBins_y_CorrelationsVs[eCorrelations_vs_MeanPhi] = 200;
  min_y_CorrelationsVs[eCorrelations_vs_MeanPhi] = 2.;
  max_y_CorrelationsVs[eCorrelations_vs_MeanPhi] = 4.;
  title_y_CorrelationsVs[eCorrelations_vs_MeanPhi] = FancyFormatting(ph.fParticleHistogramsName[ePhi].Data());

  // *) "eCorrelations_vs_MeanPt":
  nBins_y_CorrelationsVs[eCorrelations_vs_MeanPt] = 200;
  min_y_CorrelationsVs[eCorrelations_vs_MeanPt] = 0.0;
  max_y_CorrelationsVs[eCorrelations_vs_MeanPt] = 2.0;
  title_y_CorrelationsVs[eCorrelations_vs_MeanPt] = FancyFormatting(ph.fParticleHistogramsName[ePt].Data());

  // *) "eCorrelations_vs_MeanEta":
  nBins_y_CorrelationsVs[eCorrelations_vs_MeanEta] = 600;
  min_y_CorrelationsVs[eCorrelations_vs_MeanEta] = -0.3;
  max_y_CorrelationsVs[eCorrelations_vs_MeanEta] = 0.3;
  title_y_CorrelationsVs[eCorrelations_vs_MeanEta] = FancyFormatting(ph.fParticleHistogramsName[eEta].Data());

  // .....

  // *) Quick insanity check on title_x_CorrelationsVs and title_y_CorrelationsVs:
  if (title_x_CorrelationsVs.EqualTo("")) {
    LOGF(fatal, "\033[1;31m%s at line %d : title_x_CorrelationsVs is not set, check corresponding enum \033[0m", __FUNCTION__, __LINE__);
  }
  for (int t = 0; t < eQACorrelationsVsHistograms2D_N; t++) {
    if (title_y_CorrelationsVs[t].EqualTo("")) {
      LOGF(fatal, "\033[1;31m%s at line %d : title_y_CorrelationsVs[%d] is not set, check corresponding enum  \033[0m", __FUNCTION__, __LINE__, t);
    }
  }

  // Okay, let's book 'em all:
  for (int t = 0; t < eQACorrelationsVsHistograms2D_N; t++) // type, see enum eQACorrelationsVsHistograms2D
  {
    if (!qa.fBookQACorrelationsVsHistograms2D[t]) {
      continue;
    }

    for (int h = 0; h < gMaxHarmonic; h++) {

      if (h + 1 < qa.fQACorrelationsVsHistogramsMinMaxHarmonic[eMin] || h + 1 >= qa.fQACorrelationsVsHistogramsMinMaxHarmonic[eMax]) {
        continue;
      }

      for (int rs = 0; rs < 2; rs++) // reco/sim
      {

        if (Skip(rs)) {
          continue;
        }

        qa.fQACorrelationsVsHistograms2D[t][h][rs] = new TH2F(
          TString::Format("fQACorrelationsVsHistograms2D[%s][%d][%s]", qa.fQACorrelationsVsHistogramsName2D[t].Data(), h, gc.srs[rs].Data()),
          TString::Format("%s, %s", "__RUN_NUMBER__", gc.srs_long[rs].Data()), // __RUN_NUMBER__ is handled in PropagateRunNumber(...)
          nBins_x_CorrelationsVs, min_x_CorrelationsVs, max_x_CorrelationsVs, nBins_y_CorrelationsVs[t], min_y_CorrelationsVs[t], max_y_CorrelationsVs[t]);

        qa.fQACorrelationsVsHistograms2D[t][h][rs]->GetXaxis()->SetTitle(TString::Format("%s (harmonic = %d)", title_x_CorrelationsVs.Data(), h + 1));
        qa.fQACorrelationsVsHistograms2D[t][h][rs]->GetYaxis()->SetTitle(title_y_CorrelationsVs[t].Data());
        qa.fQACorrelationsVsHistograms2D[t][h][rs]->SetLineColor(ec.fBeforeAfterColor[eAfter]);
        qa.fQACorrelationsVsHistograms2D[t][h][rs]->SetFillColor(ec.fBeforeAfterColor[eAfter] - 10);
        qa.fQACorrelationsVsHistograms2D[t][h][rs]->SetOption("col");
        qa.fQACorrelationsVsList->Add(qa.fQACorrelationsVsHistograms2D[t][h][rs]);
      } // for(int rs=0;rs<2;rs++) // reco/sim

    } // for (int h = 0; h < gMaxHarmonic; h++) {

  } // for(int t=0;t<eQACorrelationsVsHistograms2D_N;t++) // type, see enum eCorrelationsVsHistograms2D

  // ...

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookQAHistograms()

//============================================================

void BookEventHistograms()
{
  // Book all event histograms.

  // a) Book the profile holding flags;
  // b) Book specific event histograms 1D;
  // c) Book specific event histograms 2D.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  eh.fEventHistogramsPro = new TProfile("fEventHistogramsPro", "flags for event histograms", 1, 0., 1.);
  eh.fEventHistogramsPro->SetStats(false);
  eh.fEventHistogramsPro->SetLineColor(eColor);
  eh.fEventHistogramsPro->SetFillColor(eFillColor);
  eh.fEventHistogramsPro->GetXaxis()->SetBinLabel(1, "fFillEventHistograms");
  eh.fEventHistogramsPro->Fill(0.5, static_cast<int>(eh.fFillEventHistograms));
  // ...
  eh.fEventHistogramsList->Add(eh.fEventHistogramsPro);

  // b) Book specific control event histograms 1D:
  // ...

  for (int t = 0; t < eEventHistograms_N; t++) // type, see enum eEventHistograms
  {
    if (!eh.fBookEventHistograms[t]) {
      continue;
    }

    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      if (Skip(rs)) {
        continue;
      }

      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {

        // Special treatment for eMultiplicity => I will never fill this one before the cuts, if Multiplicity = SelectedTracks, obviously:
        if (ba == eBefore && eh.fEventHistogramsName[t].EqualTo("Multiplicity") && ec.fsEventCuts[eMultiplicityEstimator].EqualTo("SelectedTracks", TString::kIgnoreCase)) {
          // TBI 20241123 what remains ill-defined is the case when Multiplicity != SelectedTracks , check that further
          continue;
        }
        eh.fEventHistograms[t][rs][ba] = new TH1F(
          TString::Format("fEventHistograms[%s][%s][%s]", eh.fEventHistogramsName[t].Data(), gc.srs[rs].Data(), gc.sba[ba].Data()),
          TString::Format("%s, %s, %s", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.sba_long[ba].Data()), // __RUN_NUMBER__ is handled in PropagateRunNumber(...)
          static_cast<int>(eh.fEventHistogramsBins[t][0]),
          eh.fEventHistogramsBins[t][1], eh.fEventHistogramsBins[t][2]);
        eh.fEventHistograms[t][rs][ba]->GetXaxis()->SetTitle(FancyFormatting(eh.fEventHistogramsName[t].Data()));
        eh.fEventHistograms[t][rs][ba]->SetLineColor(ec.fBeforeAfterColor[ba]);
        eh.fEventHistograms[t][rs][ba]->SetFillColor(ec.fBeforeAfterColor[ba] - 10);
        eh.fEventHistogramsList->Add(eh.fEventHistograms[t][rs][ba]);
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for(int t=0;t<eEventHistograms_N;t++) // type, see enum eEventHistograms

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookEventHistograms()

//============================================================

void BookEventCutsHistograms()
{
  // Book all event cuts objects.

  // a) Book the profile holding event cuts flags;
  // b) Book event cut counter maps;
  // c) Book event cut counter histograms.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  ec.fEventCutsPro = new TProfile("fEventCutsPro", "flags for event cuts", eEventCuts_N, -0.5, static_cast<float>(eEventCuts_N) - 0.5);
  if (tc.fUseSpecificCuts) {
    ec.fEventCutsPro->SetTitle(TString::Format("%s (hardwired analysis-specific cuts = %s)", ec.fEventCutsPro->GetTitle(), tc.fWhichSpecificCuts.Data()).Data());
  } else {
    ec.fEventCutsPro->SetTitle(TString::Format("%s (hardwired analysis-specific cuts not used)", ec.fEventCutsPro->GetTitle()).Data());
  }
  ec.fEventCutsPro->SetStats(false);
  ec.fEventCutsPro->SetLineColor(eColor);
  ec.fEventCutsPro->SetFillColor(eFillColor);
  ec.fEventCutsPro->GetXaxis()->SetLabelSize(0.025);
  for (int cut = 0; cut < eEventCuts_N; cut++) {
    ec.fEventCutsPro->GetXaxis()->SetBinLabel(1 + cut, ec.fEventCutName[cut].Data()); // Remark: check always if bin labels here correspond to ordering in enum eEventCuts
    ec.fEventCutsPro->Fill(cut, static_cast<int>(ec.fUseEventCuts[cut]));
  }
  ec.fEventCutsList->Add(ec.fEventCutsPro);

  // b) Book event cut counter maps:
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {
    // If I am analyzing only reconstructed data, do not book maps for simulated, and vice versa.
    if ((tc.fProcess[eGenericRec] && rs == eSim) || (tc.fProcess[eGenericSim] && rs == eRec)) {
      continue;
    }
    ec.fEventCutCounterMap[rs] = new TExMap();
    ec.fEventCutCounterMapInverse[rs] = new TExMap();
  }

  // c) Book event cut counter histograms:
  // ...
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {

    if (Skip(rs)) {
      continue;
    }

    for (int cc = 0; cc < eCutCounter_N; cc++) // cut counter
    {

      if ((!ec.fUseEventCutCounterAbsolute && cc == eAbsolute) || (!ec.fUseEventCutCounterSequential && cc == eSequential)) {
        continue;
      }

      ec.fEventCutCounterHist[rs][cc] = new TH1I(TString::Format("fEventCutCounterHist[%s][%s]", gc.srs[rs].Data(), gc.scc[cc].Data()), TString::Format("%s, %s, event cut counter (%s)", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.scc_long[cc].Data()), eEventCuts_N, 0.5, static_cast<double>(eEventCuts_N) + 0.5); // I cast in double the last argument, because that's what this particular TH1I constructor expects. And yes, +0.5, because eEventCuts kicks off from 0
      ec.fEventCutCounterHist[rs][cc]->SetStats(false);
      ec.fEventCutCounterHist[rs][cc]->SetLineColor(eColor);
      ec.fEventCutCounterHist[rs][cc]->SetFillColor(eFillColor);
      ec.fEventCutCounterHist[rs][cc]->GetXaxis()->SetLabelSize(0.025);

      // Remark: Bin labels are set later in a dry call to EventCuts, to accomodate sequential event cut counting
      ec.fEventCutsList->Add(ec.fEventCutCounterHist[rs][cc]);

    } // for (int cc = 0; cc < eCutCounter_N; cc++) // enum eCutCounter

  } // for (int rs = 0; rs < 2; rs++) // reco/sim

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookEventCutsHistograms()

//============================================================

void BookParticleHistograms()
{
  // Book all particle histograms.

  // a) Book the profile holding flags;
  // b) Book specific particle histograms 1D;
  // c) Book specific particle histograms 2D;
  // e) Default binning for particle sparse histograms (yes, here, see comments below);
  // d) Book specific particle sparse histograms (n-dimensions).

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  ph.fParticleHistogramsPro = new TProfile(
    "fParticleHistogramsPro", "flags for particle histograms", 1, 0., 1.);
  ph.fParticleHistogramsPro->SetStats(false);
  ph.fParticleHistogramsPro->SetLineColor(eColor);
  ph.fParticleHistogramsPro->SetFillColor(eFillColor);
  ph.fParticleHistogramsPro->GetXaxis()->SetLabelSize(0.025);
  ph.fParticleHistogramsPro->GetXaxis()->SetBinLabel(1, "fFillParticleHistograms");
  ph.fParticleHistogramsPro->Fill(0.5, static_cast<int>(ph.fFillParticleHistograms));
  // ...
  ph.fParticleHistogramsList->Add(ph.fParticleHistogramsPro);

  // b) Book specific particle histograms 1D:
  // ...
  for (int t = 0; t < eParticleHistograms_N; t++) // type, see enum eParticleHistograms
  {
    if (!ph.fBookParticleHistograms[t]) {
      continue;
    }
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {

      if (Skip(rs)) {
        continue;
      }

      // **) PDG makes sense only for Sim:
      if ((tc.fProcess[eGenericRec] || tc.fProcess[eGenericRecSim]) && rs == eRec) {
        if (t == ePDG) {
          continue;
        }
      }

      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        ph.fParticleHistograms[t][rs][ba] = new TH1F(TString::Format("fParticleHistograms[%s][%s][%s]", ph.fParticleHistogramsName[t].Data(), gc.srs[rs].Data(), gc.sba[ba].Data()),
                                                     TString::Format("%s, %s, %s", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.sba_long[ba].Data()),
                                                     static_cast<int>(ph.fParticleHistogramsBins[t][0]), ph.fParticleHistogramsBins[t][1], ph.fParticleHistogramsBins[t][2]);
        ph.fParticleHistograms[t][rs][ba]->SetLineColor(ec.fBeforeAfterColor[ba]);
        ph.fParticleHistograms[t][rs][ba]->SetFillColor(ec.fBeforeAfterColor[ba] - 10);
        ph.fParticleHistograms[t][rs][ba]->GetXaxis()->SetTitle(FancyFormatting(ph.fParticleHistogramsName[t].Data()));
        ph.fParticleHistograms[t][rs][ba]->SetMinimum(1.e-4); // so that I can switch to log scale, even if some bins are empty
        // Remark: For empty histograms, when plotting interactively, because of this line, I will get
        //   E-TCanvas::Range: illegal world coordinates range ....
        // But it's harmless, because in any case I do not care about the content of empty histogram...
        ph.fParticleHistograms[t][rs][ba]->SetOption("hist"); // do not plot marker and error (see BanishmentLoopOverParticles why errors are not reliable) for each bin, only content + filled area.
        ph.fParticleHistogramsList->Add(ph.fParticleHistograms[t][rs][ba]);
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for(int t=0;t<eParticleHistograms_N;t++) // type, see enum eParticleHistograms

  // c) Book specific particle histograms 2D:
  // keep ordering in sync. with enum eParticleHistograms2D
  TString stitleX2D[] = {FancyFormatting(ph.fParticleHistogramsName[ePhi].Data()), FancyFormatting(ph.fParticleHistogramsName[ePhi].Data())};
  TString stitleY2D[] = {FancyFormatting(ph.fParticleHistogramsName[ePt].Data()), FancyFormatting(ph.fParticleHistogramsName[eEta].Data())};

  if (sizeof(stitleX2D) / sizeof(stitleX2D[0]) != eParticleHistograms2D_N) {
    LOGF(info, "\033[1;31m mismatch - add same number of names for 2D particle histograms as you have data members \033[0m");
    LOGF(info, "\033[1;31m sizeof(stitleX2D)/sizeof(stitleX2D[0]) = %d \033[0m", sizeof(stitleX2D) / sizeof(stitleX2D[0]));
    LOGF(info, "\033[1;31m eParticleHistograms2D_N = %d \033[0m", static_cast<int>(eParticleHistograms2D_N));
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }
  if (sizeof(stitleY2D) / sizeof(stitleY2D[0]) != eParticleHistograms2D_N) {
    LOGF(info, "\033[1;31m mismatch - add same number of names for 2D particle histograms as you have data members \033[0m");
    LOGF(info, "\033[1;31m sizeof(stitleY2D)/sizeof(stitleY2D[0]) = %d \033[0m", sizeof(stitleY2D) / sizeof(stitleY2D[0]));
    LOGF(info, "\033[1;31m eParticleHistograms2D_N = %d \033[0m", static_cast<int>(eParticleHistograms2D_N));
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  for (int t = 0; t < eParticleHistograms2D_N; t++) // type, see enum eParticleHistograms2D
  {
    if (!ph.fBookParticleHistograms2D[t]) {
      continue;
    }
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {

      if (Skip(rs)) {
        continue;
      }

      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {

        // optional variable-length binning for y-axis (for supported observables):
        if (ph.fParticleHistogramsName2D[t].EqualTo("Phi_vs_Pt") && res.fUseResultsProVariableLengthBins[AFO_PT]) {

          // Remark: placeholder __RUN_NUMBER__ is handled in PropagateRunNumber(...)

          // *) variable-length binning for phi vs pt, but only in pt axis:
          ph.fParticleHistograms2D[t][rs][ba] = new TH2D(TString::Format("fParticleHistograms2D[%s][%s][%s]", ph.fParticleHistogramsName2D[t].Data(), gc.srs[rs].Data(), gc.sba[ba].Data()),
                                                         TString::Format("%s, %s, %s", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.sba_long[ba].Data()),
                                                         static_cast<int>(ph.fParticleHistogramsBins2D[t][eX][0]), ph.fParticleHistogramsBins2D[t][eX][1], ph.fParticleHistogramsBins2D[t][eX][2],
                                                         res.fResultsPro[AFO_PT]->GetXaxis()->GetXbins()->GetSize() - 1, res.fResultsPro[AFO_PT]->GetXaxis()->GetXbins()->GetArray()); // yes, x-axis of "results vs pt" hist is y-axis here for 2D.
        } else if (ph.fParticleHistogramsName2D[t].EqualTo("Phi_vs_Eta") && res.fUseResultsProVariableLengthBins[AFO_ETA]) {

          // *) variable-length binning for phi vs eta, but only in eta axis:
          ph.fParticleHistograms2D[t][rs][ba] = new TH2D(TString::Format("fParticleHistograms2D[%s][%s][%s]", ph.fParticleHistogramsName2D[t].Data(), gc.srs[rs].Data(), gc.sba[ba].Data()),
                                                         TString::Format("%s, %s, %s", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.sba_long[ba].Data()),
                                                         static_cast<int>(ph.fParticleHistogramsBins2D[t][eX][0]), ph.fParticleHistogramsBins2D[t][eX][1], ph.fParticleHistogramsBins2D[t][eX][2],
                                                         res.fResultsPro[AFO_ETA]->GetXaxis()->GetXbins()->GetSize() - 1, res.fResultsPro[AFO_ETA]->GetXaxis()->GetXbins()->GetArray()); // yes, x-axis of "results vs pt" hist is y-axis here for 2D
        } else {
          // default fixed-length binning:
          // Remark: Remember that I cannot use here  GetXaxis()->GetXbins()->GetArray() as for variable-width case, because for fixed-width case, this is always 0
          //         See https://root-forum.cern.ch/t/get-bin-array/7276/9
          ph.fParticleHistograms2D[t][rs][ba] = new TH2D(TString::Format("fParticleHistograms2D[%s][%s][%s]", ph.fParticleHistogramsName2D[t].Data(), gc.srs[rs].Data(), gc.sba[ba].Data()),
                                                         TString::Format("%s, %s, %s", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.sba_long[ba].Data()),
                                                         static_cast<int>(ph.fParticleHistogramsBins2D[t][eX][0]), ph.fParticleHistogramsBins2D[t][eX][1], ph.fParticleHistogramsBins2D[t][eX][2],
                                                         static_cast<int>(ph.fParticleHistogramsBins2D[t][eY][0]), ph.fParticleHistogramsBins2D[t][eY][1], ph.fParticleHistogramsBins2D[t][eY][2]);
        }
        ph.fParticleHistograms2D[t][rs][ba]->SetLineColor(ec.fBeforeAfterColor[ba]);
        ph.fParticleHistograms2D[t][rs][ba]->SetFillColor(ec.fBeforeAfterColor[ba] - 10);
        ph.fParticleHistograms2D[t][rs][ba]->GetXaxis()->SetTitle(stitleX2D[t].Data());
        ph.fParticleHistograms2D[t][rs][ba]->GetYaxis()->SetTitle(stitleY2D[t].Data());
        ph.fParticleHistogramsList->Add(ph.fParticleHistograms2D[t][rs][ba]);
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for(int t=0;t<eParticleHistograms_N;t++) // type, see enum
    // eParticleHistograms

  // d) Default binning for particle sparse histograms:
  //    Remark 0: This requires the special treatment, because I re-use in some cases bins from results histograns.
  //              Therefore, I can do all this only after BookResultsHistograms() was already called.
  //    Remark 1: I anticipate I will need them only when I need to calculate differential weights, therefore I couple them intentionally
  //              with enum's for differential weights from very beginning.
  //    Remark 2: Whenever possible, I re-use binning from results histograms.
  //    Remark 3: For variable-length binning, for each dimension of THnSparse, I have to call SetBinEdges (see below).
  //              Therefore, to facilitate the whole procedure, fixed-length bins which I implemented directly (e.g. for phi dimension, which I do not have in results histograms),
  //              I convert also in arrays. For fixed-length bins in results histograms I do NOT have to do that, because for that case I call GetArray() in any case, which is
  //              doing such conversion automatically.

  // **) eDiffWeightCategory = eDWPhi:

  TAxis* lAxis = NULL; // local helper TAxis, to convert in one line the booking of fixed-length array into array of corresponding bin edges

  // ***) phi-axis for diff phi weights: at the moment I support only fixed-length binning, which optionally can be made finer or coarser with ph.fRebinSparse configurable:
  ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiPhiAxis] = static_cast<int>(180. / ph.fRebinSparse);
  lAxis = new TAxis(ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiPhiAxis], 0., o2::constants::math::TwoPI);
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiPhiAxis] = new TArrayD(1 + ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiPhiAxis]);
  for (int bin = 1; bin <= lAxis->GetNbins(); bin++) {
    ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiPhiAxis]->AddAt(lAxis->GetBinLowEdge(bin), bin - 1);
  }
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiPhiAxis]->AddAt(lAxis->GetBinLowEdge(1 + lAxis->GetNbins()), lAxis->GetNbins()); // special treatment for last bin
  delete lAxis;
  ph.fParticleSparseHistogramsAxisTitle[eDWPhi][wPhiPhiAxis] = FancyFormatting("Phi");

  // ***) pt-axis for diff phi weights: I re-use binning from results histograms
  ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiPtAxis] = res.fResultsPro[AFO_PT]->GetNbinsX();
  lAxis = res.fResultsPro[AFO_PT]->GetXaxis();
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiPtAxis] = new TArrayD(1 + ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiPtAxis]);
  for (int bin = 1; bin <= lAxis->GetNbins(); bin++) {
    ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiPtAxis]->AddAt(lAxis->GetBinLowEdge(bin), bin - 1);
  }
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiPtAxis]->AddAt(lAxis->GetBinLowEdge(1 + lAxis->GetNbins()), lAxis->GetNbins()); // special treatment for last bin
  // delete lAxis; // I do not need to delete here, only when new TAxis(...)
  ph.fParticleSparseHistogramsAxisTitle[eDWPhi][wPhiPtAxis] = FancyFormatting("Pt");

  // ***) eta-axis for diff phi weights: I re-use binning from results histograms
  ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiEtaAxis] = res.fResultsPro[AFO_ETA]->GetNbinsX();
  lAxis = res.fResultsPro[AFO_ETA]->GetXaxis();
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiEtaAxis] = new TArrayD(1 + ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiEtaAxis]);
  for (int bin = 1; bin <= lAxis->GetNbins(); bin++) {
    ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiEtaAxis]->AddAt(lAxis->GetBinLowEdge(bin), bin - 1);
  }
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiEtaAxis]->AddAt(lAxis->GetBinLowEdge(1 + lAxis->GetNbins()), lAxis->GetNbins()); // special treatment for last bin
  // delete lAxis; // I do not need to delete here, only when new TAxis(...)
  ph.fParticleSparseHistogramsAxisTitle[eDWPhi][wPhiEtaAxis] = FancyFormatting("Eta");

  // ***) charge-axis for diff phi weights: I support only fixed-length binning, nothing really to ever change here:
  ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiChargeAxis] = 2;
  lAxis = new TAxis(ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiChargeAxis], -1.5, 1.5);
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiChargeAxis] = new TArrayD(1 + ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiChargeAxis]);
  for (int bin = 1; bin <= lAxis->GetNbins(); bin++) {
    ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiChargeAxis]->AddAt(lAxis->GetBinLowEdge(bin), bin - 1);
  }
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiChargeAxis]->AddAt(lAxis->GetBinLowEdge(1 + lAxis->GetNbins()), lAxis->GetNbins()); // special treatment for last bin
  // delete lAxis; // I do not need to delete here, only when new TAxis(...)
  ph.fParticleSparseHistogramsAxisTitle[eDWPhi][wPhiChargeAxis] = FancyFormatting("Charge");

  // ***) centrality-axis for diff phi weights: I re-use binning from results histograms
  ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiCentralityAxis] = res.fResultsPro[AFO_CENTRALITY]->GetNbinsX();
  lAxis = res.fResultsPro[AFO_CENTRALITY]->GetXaxis();
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiCentralityAxis] = new TArrayD(1 + ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiCentralityAxis]);
  for (int bin = 1; bin <= lAxis->GetNbins(); bin++) {
    ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiCentralityAxis]->AddAt(lAxis->GetBinLowEdge(bin), bin - 1);
  }
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiCentralityAxis]->AddAt(lAxis->GetBinLowEdge(1 + lAxis->GetNbins()), lAxis->GetNbins()); // special treatment for last bin
  // delete lAxis; // I do not need to delete here, only when new TAxis(...)
  ph.fParticleSparseHistogramsAxisTitle[eDWPhi][wPhiCentralityAxis] = "Centrality"; // TBI 20250222 I cannot call here FancyFormatting for "Centrality", because ec.fsEventCuts[eCentralityEstimator] is still not fetched and set from configurable. Re-think how to proceed for this specific case.

  // ***) vertex_z-axis for diff phi weights: I re-use binning from results histograms
  ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiVertex_zAxis] = res.fResultsPro[AFO_VZ]->GetNbinsX();
  lAxis = res.fResultsPro[AFO_VZ]->GetXaxis();
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiVertex_zAxis] = new TArrayD(1 + ph.fParticleSparseHistogramsNBins[eDWPhi][wPhiVertex_zAxis]);
  for (int bin = 1; bin <= lAxis->GetNbins(); bin++) {
    ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiVertex_zAxis]->AddAt(lAxis->GetBinLowEdge(bin), bin - 1);
  }
  ph.fParticleSparseHistogramsBinEdges[eDWPhi][wPhiVertex_zAxis]->AddAt(lAxis->GetBinLowEdge(1 + lAxis->GetNbins()), lAxis->GetNbins()); // special treatment for last bin
  // delete lAxis; // I do not need to delete here, only when new TAxis(...)
  ph.fParticleSparseHistogramsAxisTitle[eDWPhi][wPhiVertex_zAxis] = "Vertex_z"; // TBI 20250222 I cannot call here FancyFormatting for "Centrality", because ec.fsEventCuts[eCentralityEstimator]

  // ...

  // **) eDiffWeightCategory = eDWPt:

  // ***) pt-axis for diff pt weights: I re-use binning from results histograms
  ph.fParticleSparseHistogramsNBins[eDWPt][wPtPtAxis] = res.fResultsPro[AFO_PT]->GetNbinsX();
  lAxis = res.fResultsPro[AFO_PT]->GetXaxis();
  ph.fParticleSparseHistogramsBinEdges[eDWPt][wPtPtAxis] = new TArrayD(1 + ph.fParticleSparseHistogramsNBins[eDWPt][wPtPtAxis]);
  for (int bin = 1; bin <= lAxis->GetNbins(); bin++) {
    ph.fParticleSparseHistogramsBinEdges[eDWPt][wPtPtAxis]->AddAt(lAxis->GetBinLowEdge(bin), bin - 1);
  }
  ph.fParticleSparseHistogramsBinEdges[eDWPt][wPtPtAxis]->AddAt(lAxis->GetBinLowEdge(1 + lAxis->GetNbins()), lAxis->GetNbins()); // special treatment for last bin
  // delete lAxis; // I do not need to delete here, only when new TAxis(...)
  ph.fParticleSparseHistogramsAxisTitle[eDWPt][wPtPtAxis] = FancyFormatting("Pt");

  // ...

  // **) eDiffWeightCategory = eDWEta:

  // ***) eta-axis for diff eta weights: I re-use binning from results histograms
  ph.fParticleSparseHistogramsNBins[eDWEta][wEtaEtaAxis] = res.fResultsPro[AFO_ETA]->GetNbinsX();
  lAxis = res.fResultsPro[AFO_ETA]->GetXaxis();
  ph.fParticleSparseHistogramsBinEdges[eDWEta][wEtaEtaAxis] = new TArrayD(1 + ph.fParticleSparseHistogramsNBins[eDWEta][wEtaEtaAxis]);
  for (int bin = 1; bin <= lAxis->GetNbins(); bin++) {
    ph.fParticleSparseHistogramsBinEdges[eDWEta][wEtaEtaAxis]->AddAt(lAxis->GetBinLowEdge(bin), bin - 1);
  }
  ph.fParticleSparseHistogramsBinEdges[eDWEta][wEtaEtaAxis]->AddAt(lAxis->GetBinLowEdge(1 + lAxis->GetNbins()), lAxis->GetNbins()); // special treatment for last bin
  // delete lAxis; // I do not need to delete here, only when new TAxis(...)
  ph.fParticleSparseHistogramsAxisTitle[eDWEta][wEtaEtaAxis] = FancyFormatting("Eta");

  // ...

  // e) Book specific particle sparse histograms (n-dimensions):
  if (ph.fBookParticleSparseHistograms[eDWPhi]) {
    BookParticleSparseHistograms(eDWPhi);
  }

  if (ph.fBookParticleSparseHistograms[eDWPt]) {
    BookParticleSparseHistograms(eDWPt);
  }

  if (ph.fBookParticleSparseHistograms[eDWEta]) {
    BookParticleSparseHistograms(eDWEta);
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookParticleHistograms()

//============================================================

void BookParticleSparseHistograms(eDiffWeightCategory dwc)
{
  // This is a helper function for BookParticleHistograms(), merely to reduce code bloat.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Determine number of dimensions for sparse histogram for this differential weight category:
  int nDimensions = -1;
  switch (dwc) {
    case eDWPhi: {
      nDimensions = static_cast<int>(eDiffPhiWeights_N);
      break;
    }
    case eDWPt: {
      nDimensions = static_cast<int>(eDiffPtWeights_N);
      break;
    }
    case eDWEta: {
      nDimensions = static_cast<int>(eDiffEtaWeights_N);
      break;
    }
    default: {
      LOGF(fatal, "\033[1;31m%s at line %d : This differential weight category, dwc = %d, is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(dwc));
      break;
    }
  } // switch(dwc)

  // *) Determine binning for all dimensions:
  TArrayI* nBins = new TArrayI(nDimensions);
  for (int d = 0; d < nDimensions; d++) {
    nBins->AddAt(static_cast<int>(ph.fParticleSparseHistogramsNBins[dwc][d]), d);
  }

  // *) Book THnSparse with correct number of bins for each dimension, but void bin edges:
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {
    if (Skip(rs)) {
      continue;
    }
    // Remark: Here I have a bit unusual convention for the name and title, but okay...
    ph.fParticleSparseHistograms[dwc][rs] = new THnSparseF(TString::Format("%s[%s]", ph.fParticleSparseHistogramsName[dwc].Data(), gc.srs[rs].Data()), TString::Format("__RUN_NUMBER__, %s, %s", gc.srs_long[rs].Data(), ph.fParticleSparseHistogramsTitle[dwc].Data()), nDimensions, nBins->GetArray(), NULL, NULL);

    // *) For each dimension set bin edges, axis title, etc.:
    for (int d = 0; d < nDimensions; d++) {
      ph.fParticleSparseHistograms[dwc][rs]->SetBinEdges(d, ph.fParticleSparseHistogramsBinEdges[dwc][d]->GetArray());
      ph.fParticleSparseHistograms[dwc][rs]->GetAxis(d)->SetTitle(ph.fParticleSparseHistogramsAxisTitle[dwc][d].Data());
    }

    // *) Finally, add the fully configured THnSparse to its TList:
    ph.fParticleHistogramsList->Add(ph.fParticleSparseHistograms[dwc][rs]);
  } // for (int rs = 0; rs < 2; rs++) // reco/sim

  // *) Clean up:
  delete nBins;

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookParticleSparseHistograms()

//============================================================

void BookParticleCutsHistograms()
{
  // Book all particle cuts objects.

  // a) Book the profile holding flags;
  // b) Book particle cut counter maps;
  // c) Book the particle cut counter (absolute);
  // d) Book the formula for pt-dependent DCAxy cut.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  pc.fParticleCutsPro = new TProfile("fParticleCutsPro", "flags for particle cuts", eParticleCuts_N, -0.5, static_cast<float>(eParticleCuts_N) - 0.5);
  if (tc.fUseSpecificCuts) {
    pc.fParticleCutsPro->SetTitle(TString::Format("%s (hardwired analysis-specific cuts = %s)", pc.fParticleCutsPro->GetTitle(), tc.fWhichSpecificCuts.Data()).Data());
  } else {
    pc.fParticleCutsPro->SetTitle(TString::Format("%s (hardwired analysis-specific cuts not used)", pc.fParticleCutsPro->GetTitle()).Data());
  }
  pc.fParticleCutsPro->SetStats(false);
  pc.fParticleCutsPro->SetLineColor(eColor);
  pc.fParticleCutsPro->SetFillColor(eFillColor);
  for (int cut = 0; cut < eParticleCuts_N; cut++) {
    pc.fParticleCutsPro->GetXaxis()->SetBinLabel(1 + cut, pc.fParticleCutName[cut].Data()); // Remark: check always if bin labels here correspond to ordering in enum eParticleCuts
    pc.fParticleCutsPro->Fill(cut, static_cast<int>(pc.fUseParticleCuts[cut]));
  }
  pc.fParticleCutsList->Add(pc.fParticleCutsPro);

  // b) Book particle cut counter maps:
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {
    // If I am analyzing only reconstructed data, do not book maps for simulated, and vice versa.
    if ((tc.fProcess[eGenericRec] && rs == eSim) || (tc.fProcess[eGenericSim] && rs == eRec)) {
      continue;
    }
    pc.fParticleCutCounterMap[rs] = new TExMap();
    pc.fParticleCutCounterMapInverse[rs] = new TExMap();
  }

  // c) Book the particle cut counter (absolute):
  // ...
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {

    if (Skip(rs)) {
      continue;
    }

    for (int cc = 0; cc < eCutCounter_N; cc++) // cut counter
    {

      if ((!pc.fUseParticleCutCounterAbsolute && cc == eAbsolute) || (!pc.fUseParticleCutCounterSequential && cc == eSequential)) {
        continue;
      }

      pc.fParticleCutCounterHist[rs][cc] = new TH1I(TString::Format("fParticleCutCounterHist[%s][%s]", gc.srs[rs].Data(), gc.scc[cc].Data()), TString::Format("%s, %s, particle cut counter (%s)", "__RUN_NUMBER__", gc.srs_long[rs].Data(), gc.scc_long[cc].Data()), eParticleCuts_N, 0.5, static_cast<double>(eParticleCuts_N) + 0.5);
      // I cast in double the last argument, because that's what this particular TH1I constructor expects
      // Yes, +0.5, because eParticleCuts kicks off from 0
      pc.fParticleCutCounterHist[rs][cc]->SetStats(false);
      pc.fParticleCutCounterHist[rs][cc]->SetLineColor(eColor);
      pc.fParticleCutCounterHist[rs][cc]->SetFillColor(eFillColor);
      pc.fParticleCutCounterHist[rs][cc]->GetXaxis()->SetLabelSize(0.025);
      // Remark: Bin labels are set later in a dry call to ParticleCuts, to accomodate sequential particle cut counting
      pc.fParticleCutsList->Add(pc.fParticleCutCounterHist[rs][cc]);

    } // for (int cc = 0; cc < eCutCounter_N; cc++) // enum eCutCounter

  } // for (int rs = 0; rs < 2; rs++) // reco/sim

  // d) Book the formula for pt-dependent DCAxy cut:
  if (pc.fUseParticleCuts[ePtDependentDCAxyParameterization]) {
    pc.fPtDependentDCAxyFormula = new TFormula("fPtDependentDCAxyFormula", pc.fsParticleCuts[ePtDependentDCAxyParameterization].Data());
    // As a quick insanity check, try immediately to evaluate something from this formula:
    if (std::isnan(pc.fPtDependentDCAxyFormula->Eval(1.44))) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
  } // if(pc.fUseParticleCuts[ePtDependentDCAxyParameterization]) {

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookParticleCutsHistograms()

//============================================================

void BookQvectorHistograms()
{
  // Book all Q-vector histograms.

  // a) Book the profile holding flags;
  // b) Book multiplicity distributions in A and B, for each eta separation;
  // c) ...

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  qv.fQvectorFlagsPro =
    new TProfile("fQvectorFlagsPro", "flags for Q-vector objects", 3, 0., 3.);
  qv.fQvectorFlagsPro->SetStats(false);
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

  // b) Book multiplicity distributions in A and B, for each eta separation:
  if (es.fCalculateEtaSeparations) {
    TString sEtaSep[2] = {"A", "B"}; // A <=> -eta , B <=> + eta
    TString sEtaSep_long[2] = {TString::Format("%.2f < #eta <", pc.fdParticleCuts[eEta][eMin]), TString::Format("< #eta < %.2f", pc.fdParticleCuts[eEta][eMax])};
    // yes, here I define first the part of intervals as etaCutMin < eta < "subevent boundary", and "subevent" boundary < eta < etaCutMax
    // Then below in the loop, I inject for "subevent boundary" the corresponding fEtaSeparationsValues (devided by 2, becaus it's symmetric round 0)
    for (int ab = 0; ab < 2; ab++) {   // ab = 0 <=> -eta , ab = 1 <=> + eta
      for (int rs = 0; rs < 2; rs++) { // reco/sim
        if (Skip(rs)) {
          continue;
        }
        for (int ba = 0; ba < 2; ba++) { // before/after cuts
          if (eBefore == ba) {
            continue; // it make sense to fill these histos only for "eAfter", because Q-vectors are not filled for "eBefore"
          }
          for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
            qv.fMabDist[ab][rs][ba][e] = new TH1F(Form("fMabDist[%s][%s][%s][%d]", sEtaSep[ab].Data(), gc.srs[rs].Data(), gc.sba[ba].Data(), e),
                                                  Form("%s, %s, %s, %s", "__RUN_NUMBER__",
                                                       0 == ab ? Form("%s -%.2f", sEtaSep_long[ab].Data(), es.fEtaSeparationsValues[e] / 2.) : Form("%.2f %s", es.fEtaSeparationsValues[e] / 2., sEtaSep_long[ab].Data()), gc.srs_long[rs].Data(), gc.sba_long[ba].Data()),
                                                  static_cast<int>(eh.fEventHistogramsBins[eMultiplicity][0]), eh.fEventHistogramsBins[eMultiplicity][1], eh.fEventHistogramsBins[eMultiplicity][2]); // TBI 20241207 I have hardwired in this constructor "0 == ab", this can backfire...
            qv.fMabDist[ab][rs][ba][e]->SetLineColor(ec.fBeforeAfterColor[ba]);
            qv.fMabDist[ab][rs][ba][e]->SetFillColor(ec.fBeforeAfterColor[ba] - 10);
            qv.fMabDist[ab][rs][ba][e]->GetXaxis()->SetTitle("subevent multiplicity (sum of particle weights)");
            qv.fMabDist[ab][rs][ba][e]->SetMinimum(1.e-4); // so that I can switch to log scale, even if some bins are empty
            // Remark: For empty histograms, when plotting interactively, because of this line, I will get
            //   E-TCanvas::Range: illegal world coordinates range ....
            // But it's harmless, because in any case I do not care about the content of empty histogram...
            qv.fMabDist[ab][rs][ba][e]->SetOption("hist"); // do not plot marker and error (see BanishmentLoopOverParticles why errors are not reliable) for each bin, only content + filled area.
            qv.fQvectorList->Add(qv.fMabDist[ab][rs][ba][e]);
          }
        }
      }
    }
  } // if (es.fCalculateEtaSeparations) {

  // c) ...

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

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
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  mupa.fCorrelationsFlagsPro = new TProfile("fCorrelationsFlagsPro",
                                            "flags for correlations", 1, 0., 1.);
  mupa.fCorrelationsFlagsPro->SetStats(false);
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
  for (int k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]
  {
    for (int n = 0; n < gMaxHarmonic; n++) // harmonic
    {
      for (int v = 0; v < eAsFunctionOf_N; v++) {

        // decide what is booked, then later valid pointer to fCorrelationsPro[k][n][v] is used as a boolean, in the standard way:
        if (AFO_INTEGRATED == v && !mupa.fCalculateCorrelationsAsFunctionOf[AFO_INTEGRATED]) {
          continue;
        }
        if (AFO_MULTIPLICITY == v && !mupa.fCalculateCorrelationsAsFunctionOf[AFO_MULTIPLICITY]) {
          continue;
        }
        if (AFO_CENTRALITY == v && !mupa.fCalculateCorrelationsAsFunctionOf[AFO_CENTRALITY]) {
          continue;
        }
        if (AFO_PT == v && !mupa.fCalculateCorrelationsAsFunctionOf[AFO_PT]) {
          continue;
        }
        if (AFO_ETA == v && !mupa.fCalculateCorrelationsAsFunctionOf[AFO_ETA]) {
          continue;
        }
        if (AFO_OCCUPANCY == v && !mupa.fCalculateCorrelationsAsFunctionOf[AFO_OCCUPANCY]) {
          continue;
        }
        if (AFO_INTERACTIONRATE == v && !mupa.fCalculateCorrelationsAsFunctionOf[AFO_INTERACTIONRATE]) {
          continue;
        }
        if (AFO_CURRENTRUNDURATION == v && !mupa.fCalculateCorrelationsAsFunctionOf[AFO_CURRENTRUNDURATION]) {
          continue;
        }
        if (AFO_VZ == v && !mupa.fCalculateCorrelationsAsFunctionOf[AFO_VZ]) {
          continue;
        }

        if (!res.fResultsPro[v]) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
        mupa.fCorrelationsPro[k][n][v] = reinterpret_cast<TProfile*>(res.fResultsPro[v]->Clone(Form("fCorrelationsPro[%d][%d][%s]", k, n, res.fResultsProRawName[v].Data()))); // yes
        mupa.fCorrelationsPro[k][n][v]->SetStats(false);
        mupa.fCorrelationsPro[k][n][v]->Sumw2();
        mupa.fCorrelationsPro[k][n][v]->GetXaxis()->SetTitle(FancyFormatting(res.fResultsProXaxisTitle[v].Data()));
        mupa.fCorrelationsPro[k][n][v]->GetYaxis()->SetTitle(Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1), oVariable[k].Data()));
        mupa.fCorrelationsList->Add(mupa.fCorrelationsPro[k][n][v]);
      }
    } // for (int n = 0; n < gMaxHarmonic; n++) // harmonic
  } // for (int k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]

  // d) Few quick insanity checks on booking:
  if (mupa.fCorrelationsPro[0][0][AFO_INTEGRATED] && !TString(mupa.fCorrelationsPro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }
  if (mupa.fCorrelationsPro[0][0][AFO_PT] && !TString(mupa.fCorrelationsPro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("p_{T}")) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // BookCorrelationsHistograms()

//============================================================

void BookWeightsHistograms()
{
  // Book all objects for particle weights.

  // a) Book the profile holding flags;
  // b) Histograms for integrated weights;
  // c) Histograms for differential weights;
  // d) Sparse histograms for differential phi weights.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  pw.fWeightsFlagsPro = new TProfile("fWeightsFlagsPro", "flags for particle weights", 13, 0., 13.);
  pw.fWeightsFlagsPro->SetStats(false);
  pw.fWeightsFlagsPro->SetLineColor(eColor);
  pw.fWeightsFlagsPro->SetFillColor(eFillColor);
  pw.fWeightsFlagsPro->GetXaxis()->SetLabelSize(0.035);
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(1, "w_{#varphi}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(2, "w_{p_{t}}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(3, "w_{#eta}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(4, "(w_{#varphi})_{| p_{T}}"); // TBI 20241019 not sure if this is the final notation, keep in sync with void SetDiffWeightsHist(...)
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(5, "(w_{#varphi})_{| #eta}");  // TBI 20241019 not sure if this is the final notation, keep in sync with void SetDiffWeightsHist(...)

  // **) differential phi weights using sparse:
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(6, "(w_{#varphi})_{phi axis (sparse)}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(7, "(w_{#varphi})_{p_{T} axis (sparse)}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(8, "(w_{#varphi})_{#eta axis (sparse)}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(9, "(w_{#varphi})_{charge axis (sparse)}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(10, "(w_{#varphi})_{centrality axis (sparse)}");
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(11, "(w_{#varphi})_{vertex_z axis (sparse)}");

  // **) differential pt weights using sparse:
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(12, "(w_{p_{T}})_{pt axis (sparse)}");

  // **) differential eta weights using sparse:
  pw.fWeightsFlagsPro->GetXaxis()->SetBinLabel(13, "(w_{#eta})_{eta axis (sparse)}");

  for (int w = 0; w < eWeights_N; w++) // use weights [phi,pt,eta]
  {
    if (pw.fUseWeights[w]) {
      pw.fWeightsFlagsPro->Fill(w + 0.5, 1.);
    }
  }

  // **) use differential weights [phipt,phieta,...] TBI 20250514 obsolete
  if (pw.fUseDiffWeights[wPHIPT]) {
    pw.fWeightsFlagsPro->Fill(3.5, 1.);
  }
  if (pw.fUseDiffWeights[wPHIETA]) {
    pw.fWeightsFlagsPro->Fill(4.5, 1.);
  }

  // **) differential phi weights using sparse:
  if (pw.fUseDiffPhiWeights[wPhiPhiAxis]) {
    pw.fWeightsFlagsPro->Fill(5.5, 1.);
  }
  if (pw.fUseDiffPhiWeights[wPhiPtAxis]) {
    pw.fWeightsFlagsPro->Fill(6.5, 1.);
  }
  if (pw.fUseDiffPhiWeights[wPhiEtaAxis]) {
    pw.fWeightsFlagsPro->Fill(7.5, 1.);
  }
  if (pw.fUseDiffPhiWeights[wPhiChargeAxis]) {
    pw.fWeightsFlagsPro->Fill(8.5, 1.);
  }
  if (pw.fUseDiffPhiWeights[wPhiCentralityAxis]) {
    pw.fWeightsFlagsPro->Fill(9.5, 1.);
  }
  if (pw.fUseDiffPhiWeights[wPhiVertex_zAxis]) {
    pw.fWeightsFlagsPro->Fill(10.5, 1.);
  }

  // **) differential pt weights using sparse:
  if (pw.fUseDiffPtWeights[wPtPtAxis]) {
    pw.fWeightsFlagsPro->Fill(11.5, 1.);
  }

  // **) differential eta weights using sparse:
  if (pw.fUseDiffEtaWeights[wEtaEtaAxis]) {
    pw.fWeightsFlagsPro->Fill(12.5, 1.);
  }

  pw.fWeightsList->Add(pw.fWeightsFlagsPro);

  // b) Histograms for integrated weights:
  //    As of 20240216, I have abandoned the idea to generate integrated weights internally, weights
  //    are always fetched and cloned from external files, in any case (local, AliEn, CCDB).
  //    Therefore, add histos with weights to this list only after they are cloned from external files.

  // c) Histograms for differential weights:
  //    Same comment applies as for b) => add histograms to the list, only after they are cloned from external files.

  // d) Sparse histograms for differential phi weights:
  //    Same comment applies as for b) => add sparse histograms to the list, only after they are cloned from external files.

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookWeightsHistograms()

//============================================================

void BookCentralityWeightsHistograms()
{
  // Book all objects for centrality weights.

  // a) Book the profile holding flags;
  // b) Histograms for centrality weights.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  cw.fCentralityWeightsFlagsPro =
    new TProfile("fCentralityWeightsFlagsPro", "flags for centrality weights", 1, 0., 1.);
  cw.fCentralityWeightsFlagsPro->SetStats(false);
  cw.fCentralityWeightsFlagsPro->SetLineColor(eColor);
  cw.fCentralityWeightsFlagsPro->SetFillColor(eFillColor);
  cw.fCentralityWeightsFlagsPro->GetXaxis()->SetLabelSize(0.05);
  cw.fCentralityWeightsFlagsPro->GetXaxis()->SetBinLabel(1, TString::Format("Use centrality weights for estimator %s", ec.fsEventCuts[eCentralityEstimator].Data()));
  if (cw.fUseCentralityWeights) {
    cw.fCentralityWeightsFlagsPro->Fill(0.5, 1.); // TBI 20241118 shall I automate this?
  }
  cw.fCentralityWeightsList->Add(cw.fCentralityWeightsFlagsPro);

  // b) Histograms for centrality weights:
  //    As of 20240216, I have abandoned the idea to generate centrality weights internally, centrality weights
  //    are always fetched and cloned from external files, in any case (local, AliEn, CCDB).
  //    Therefore, add histos with centrality weights to this list only after they are cloned from external files.

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookCentralityWeightsHistograms()

//============================================================

void BookNestedLoopsHistograms()
{
  // Book all nested loops histograms.

  // a) Book the profile holding flags;
  // b) Common local labels (keep 'em in sync with BookCorrelationsHistograms());
  // c) Book what needs to be booked;
  // d) Few quick insanity checks on booking.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  nl.fNestedLoopsFlagsPro =
    new TProfile("fNestedLoopsFlagsPro", "flags for nested loops", 4, 0., 4.);
  nl.fNestedLoopsFlagsPro->SetStats(false);
  nl.fNestedLoopsFlagsPro->SetLineColor(eColor);
  nl.fNestedLoopsFlagsPro->SetFillColor(eFillColor);
  nl.fNestedLoopsFlagsPro->GetXaxis()->SetLabelSize(0.03);
  nl.fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateNestedLoops");
  nl.fNestedLoopsFlagsPro->Fill(0.5, nl.fCalculateNestedLoops);
  nl.fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(2, "fCalculateCustomNestedLoops");
  nl.fNestedLoopsFlagsPro->Fill(1.5, nl.fCalculateCustomNestedLoops);
  nl.fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(3, "fCalculateKineCustomNestedLoops");
  nl.fNestedLoopsFlagsPro->Fill(2.5, nl.fCalculateKineCustomNestedLoops);
  nl.fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(4, "fMaxNestedLoop");
  nl.fNestedLoopsFlagsPro->Fill(3.5, nl.fMaxNestedLoop);
  nl.fNestedLoopsList->Add(nl.fNestedLoopsFlagsPro);

  if (!(nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoops || nl.fCalculateKineCustomNestedLoops)) {
    // TBI 20240326 I have to keep all 3 flags above, because for instance TArrayD* ftaNestedLoops[2] is used as a storage both for nl.fCalculateNestedLoops and nl.fCalculateCustomNestedLoops
    return;
  }

  // *) Book containers for integrated nested loops:
  if (nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoops) {
    const int iMaxSize = 2e4;
    nl.ftaNestedLoops[0] = new TArrayD(iMaxSize); // ebe container for azimuthal angles
    nl.ftaNestedLoops[1] = new TArrayD(iMaxSize); // ebe container for particle weights (product of all)
  }

  // *) Book containers for differential nested loops:
  if (nl.fCalculateKineCustomNestedLoops) {
    const int iMaxSize = 2e4;
    for (int b = 0; b < res.fResultsPro[AFO_PT]->GetNbinsX(); b++) {
      nl.ftaNestedLoopsKine[PTq][b][0] = new TArrayD(iMaxSize);
      nl.ftaNestedLoopsKine[PTq][b][1] = new TArrayD(iMaxSize);
    }
    for (int b = 0; b < res.fResultsPro[AFO_ETA]->GetNbinsX(); b++) {
      nl.ftaNestedLoopsKine[ETAq][b][0] = new TArrayD(iMaxSize);
      nl.ftaNestedLoopsKine[ETAq][b][1] = new TArrayD(iMaxSize);
    }
  }

  // b) Common local labels (keep 'em in sync with BookCorrelationsHistograms())
  TString oVariable[4] = {
    "#varphi_{1}-#varphi_{2}",
    "#varphi_{1}+#varphi_{2}-#varphi_{3}-#varphi_{4}",
    "#varphi_{1}+#varphi_{2}+#varphi_{3}-#varphi_{4}-#varphi_{5}-#varphi_{6}",
    "#varphi_{1}+#varphi_{2}+#varphi_{3}+#varphi_{4}-#varphi_{5}-#varphi_{6}-"
    "#varphi_{7}-#varphi_{8}"};

  // c) Book what needs to be booked:
  if (!(nl.fCalculateNestedLoops)) { // TBI 20240404 for the time being, I can keep it here, but eventualy it will have to go elsewhere
    return;
  }
  for (int k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]
  {
    // TBI 20240405 I could break here, with respect to what nl.fMaxNestedLoop was set to

    for (int n = 0; n < gMaxHarmonic; n++) // harmonic
    {
      for (int v = 0; v < eAsFunctionOf_N; v++) {

        if (!res.fResultsPro[v]) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
        nl.fNestedLoopsPro[k][n][v] = reinterpret_cast<TProfile*>(res.fResultsPro[v]->Clone(Form("fNestedLoopsPro[%d][%d][%d]", k, n, v))); // yes
        nl.fNestedLoopsPro[k][n][v]->SetTitle(Form("#LT#LTcos[%s(%s)]#GT#GT", 1 == n + 1 ? "" : Form("%d", n + 1), oVariable[k].Data()));
        nl.fNestedLoopsPro[k][n][v]->SetStats(false);
        nl.fNestedLoopsPro[k][n][v]->Sumw2();
        nl.fNestedLoopsPro[k][n][v]->GetXaxis()->SetTitle(res.fResultsProXaxisTitle[v].Data());

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
      } // for(int v=0;v<5;v++) // variable [0=integrated,1=vs.
        // multiplicity,2=vs. centrality]
    } // for (int n = 0; n < gMaxHarmonic; n++) // harmonic
  } // for (int k = 0; k < 4; k++) // order [2p=0,4p=1,6p=2,8p=3]

  // d) Few quick insanity checks on booking:
  if (nl.fNestedLoopsPro[0][0][AFO_INTEGRATED] && !TString(nl.fNestedLoopsPro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(info, "\033[1;33mnl.fNestedLoopsPro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle() = %s \033[0m", nl.fNestedLoopsPro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle());
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }
  if (nl.fNestedLoopsPro[0][0][AFO_PT] && !TString(nl.fNestedLoopsPro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("pt")) { // I do not need here fancy formatting
    LOGF(info, "\033[1;33mnl.fNestedLoopsPro[0][0][AFO_PT]->GetXaxis()->GetTitle() = %s \033[0m", nl.fNestedLoopsPro[0][0][AFO_PT]->GetXaxis()->GetTitle());
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookNestedLoopsHistograms()

//============================================================

void BookNUAHistograms()
{
  // Book all objects for Toy NUA.

  // a) Book the profile holding flags;
  // b) Common local labels;
  // c) Histograms.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  nua.fNUAFlagsPro = new TProfile("fNUAFlagsPro", "flags for Toy NUA", 6, 0.5, 6.5);
  nua.fNUAFlagsPro->SetStats(false);
  nua.fNUAFlagsPro->SetLineColor(eColor);
  nua.fNUAFlagsPro->SetFillColor(eFillColor);
  nua.fNUAFlagsPro->GetXaxis()->SetLabelSize(0.03);
  // TBI 20240429 the binning below is a bit fragile, but ok...
  nua.fNUAFlagsPro->GetXaxis()->SetBinLabel(static_cast<int>(1 + ePhiNUAPDF), "fApplyNUAPDF[phi]");
  nua.fNUAFlagsPro->GetXaxis()->SetBinLabel(static_cast<int>(1 + ePtNUAPDF), "fApplyNUAPDF[pt]");
  nua.fNUAFlagsPro->GetXaxis()->SetBinLabel(static_cast<int>(1 + eEtaNUAPDF), "fApplyNUAPDF[eta]");
  if (nua.fApplyNUAPDF[ePhiNUAPDF]) {
    nua.fNUAFlagsPro->Fill(static_cast<int>(1 + ePhiNUAPDF), 1.);
  }
  if (nua.fApplyNUAPDF[ePtNUAPDF]) {
    nua.fNUAFlagsPro->Fill(static_cast<int>(1 + ePtNUAPDF), 1.);
  }
  if (nua.fApplyNUAPDF[eEtaNUAPDF]) {
    nua.fNUAFlagsPro->Fill(static_cast<int>(1 + eEtaNUAPDF), 1.);
  }
  nua.fNUAFlagsPro->GetXaxis()->SetBinLabel(static_cast<int>(4 + ePhiNUAPDF), "fUseDefaultNUAPDF[phi]");
  nua.fNUAFlagsPro->GetXaxis()->SetBinLabel(static_cast<int>(4 + ePtNUAPDF), "fUseDefaultNUAPDF[pt]");
  nua.fNUAFlagsPro->GetXaxis()->SetBinLabel(static_cast<int>(4 + eEtaNUAPDF), "fUseDefaultNUAPDF[eta]");
  if (nua.fUseDefaultNUAPDF[ePhiNUAPDF]) {
    nua.fNUAFlagsPro->Fill(static_cast<int>(4 + ePhiNUAPDF), 1.);
  }
  if (nua.fUseDefaultNUAPDF[ePtNUAPDF]) {
    nua.fNUAFlagsPro->Fill(static_cast<int>(4 + ePtNUAPDF), 1.);
  }
  if (nua.fUseDefaultNUAPDF[eEtaNUAPDF]) {
    nua.fNUAFlagsPro->Fill(static_cast<int>(4 + eEtaNUAPDF), 1.);
  }
  nua.fNUAList->Add(nua.fNUAFlagsPro);

  if (!(nua.fApplyNUAPDF[ePhiNUAPDF] || nua.fApplyNUAPDF[ePtNUAPDF] || nua.fApplyNUAPDF[eEtaNUAPDF])) {
    return;
  }

  // b) Common local labels:
  TString sVariable[eNUAPDF_N] = {"#varphi", "p_{t}", "#eta"}; // has to be in sync with the ordering of enum eNUAPDF

  // c) Histograms:
  for (int pdf = 0; pdf < eNUAPDF_N; pdf++) // use pdfs for NUA in (phi, pt, eta, ...)
  {
    if (!nua.fCustomNUAPDF[pdf]) // yes, because these histos are cloned from the external ones, see void SetNUAPDF(TH1D* const hist, const char* variable);
    {
      // otherwise, book here TF1 objects with default pdfs for NUA:

      // *) default NUA for azimuthal angle pdf:
      if (sVariable[pdf].EqualTo("#varphi")) {
        if (!nua.fApplyNUAPDF[ePhiNUAPDF]) {
          continue;
        }
        // Define default detector acceptance in azimuthal angle: Two sectors, with different probabilities.
        double dFirstSector[2] = {-(3. / 4.) * o2::constants::math::PI, -(1. / 4.) * o2::constants::math::PI}; // first sector is defined as [-3Pi/4,Pi/4]
        double dSecondSector[2] = {(1. / 3.) * o2::constants::math::PI, (2. / 3.) * o2::constants::math::PI};  // second sector is defined as [Pi/3,2Pi/3]
        double dProbability[2] = {0.3, 0.5};                                                                   // probabilities
        nua.fDefaultNUAPDF[ePhiNUAPDF] = new TF1(TString::Format("fDefaultNUAPDF[%d]", ePhiNUAPDF), "1.-(x>=[0])*(1.-[4]) + (x>=[1])*(1.-[4]) - (x>=[2])*(1.-[5]) + (x>=[3])*(1.-[5]) ",
                                                 ph.fParticleHistogramsBins[ePhi][1], ph.fParticleHistogramsBins[ePhi][2]);
        nua.fDefaultNUAPDF[ePhiNUAPDF]->SetParameter(0, dFirstSector[0]);
        nua.fDefaultNUAPDF[ePhiNUAPDF]->SetParameter(1, dFirstSector[1]);
        nua.fDefaultNUAPDF[ePhiNUAPDF]->SetParameter(2, dSecondSector[0]);
        nua.fDefaultNUAPDF[ePhiNUAPDF]->SetParameter(3, dSecondSector[1]);
        nua.fDefaultNUAPDF[ePhiNUAPDF]->SetParameter(4, dProbability[0]);
        nua.fDefaultNUAPDF[ePhiNUAPDF]->SetParameter(5, dProbability[1]);
        nua.fNUAList->Add(nua.fDefaultNUAPDF[ePhiNUAPDF]);

      } else if (sVariable[pdf].EqualTo("p_{t}")) {

        // *) default NUA for transverse momentum pdf:
        if (!nua.fApplyNUAPDF[ePtNUAPDF]) {
          continue;
        }
        // Define default detector acceptance in transverse momentum: One sectors, with probability < 1.
        double dSector[2] = {0.4, 0.8}; // sector is defined as 0.8 < pT < 1.2
        double dProbability = 0.3;      // probability, so after being set this way, only 30% of particles in that sector are reconstructed
        nua.fDefaultNUAPDF[ePtNUAPDF] = new TF1(TString::Format("fDefaultNUAPDF[%d]", ePtNUAPDF), "1.-(x>=[0])*(1.-[2]) + (x>=[1])*(1.-[2])",
                                                ph.fParticleHistogramsBins[ePt][1], ph.fParticleHistogramsBins[ePt][2]);
        nua.fDefaultNUAPDF[ePtNUAPDF]->SetParameter(0, dSector[0]);
        nua.fDefaultNUAPDF[ePtNUAPDF]->SetParameter(1, dSector[1]);
        nua.fDefaultNUAPDF[ePtNUAPDF]->SetParameter(2, dProbability);
        nua.fNUAList->Add(nua.fDefaultNUAPDF[ePtNUAPDF]);

      } else if (sVariable[pdf].EqualTo("#eta")) {

        // *) default NUA for pseudorapidity pdf:
        if (!nua.fApplyNUAPDF[eEtaNUAPDF]) {
          continue;
        }
        // Define default detector acceptance in pseudorapidity: One sectors, with probability < 1.
        double dSector[2] = {2.0, 2.5}; // sector is defined as 0.5 < eta < 1.0
        double dProbability = 0.5;      // probability, so after being set this way, only 50% of particles in that sector are reconstructed
        nua.fDefaultNUAPDF[eEtaNUAPDF] = new TF1(TString::Format("fDefaultNUAPDF[%d]", eEtaNUAPDF), "1.-(x>=[0])*(1.-[2]) + (x>=[1])*(1.-[2])",
                                                 ph.fParticleHistogramsBins[eEta][1], ph.fParticleHistogramsBins[eEta][2]);
        nua.fDefaultNUAPDF[eEtaNUAPDF]->SetParameter(0, dSector[0]);
        nua.fDefaultNUAPDF[eEtaNUAPDF]->SetParameter(1, dSector[1]);
        nua.fDefaultNUAPDF[eEtaNUAPDF]->SetParameter(2, dProbability);
        nua.fNUAList->Add(nua.fDefaultNUAPDF[eEtaNUAPDF]);
      } else {
        LOGF(fatal, "\033[1;31m%s at line %d : pdf = %s is not supported (yet)\n \033[0m", __FUNCTION__, __LINE__, sVariable[pdf].Data());
      }

    } else { // if(!nua.fCustomNUAPDF[pdf])
      // generic cosmetics for custom user-supplied pdfs via histograms:
      nua.fCustomNUAPDF[pdf]->SetTitle(TString::Format("Custom user-provided NUA for %s", sVariable[pdf].Data()));
      nua.fCustomNUAPDF[pdf]->SetStats(false);
      nua.fCustomNUAPDF[pdf]->GetXaxis()->SetTitle(sVariable[pdf].Data());
      nua.fCustomNUAPDF[pdf]->SetFillColor(eFillColor);
      nua.fCustomNUAPDF[pdf]->SetLineColor(eColor);
      nua.fNUAList->Add(nua.fCustomNUAPDF[pdf]);
    } // if(!nua.fCustomNUAPDF[pdf])

    // Get the max values of pdfs, so that later in Accept(...) there is no loss of efficiency, when would need to calculate the same thing for each particle:
    if (!nua.fUseDefaultNUAPDF[pdf] && nua.fCustomNUAPDF[pdf]) { // pdf is a loop variable
      // custom, user-provided pdf via TH1D object:
      nua.fMaxValuePDF[pdf] = nua.fCustomNUAPDF[pdf]->GetMaximum();
    } else if (nua.fUseDefaultNUAPDF[pdf] && nua.fDefaultNUAPDF[pdf]) {
      // default pdf implemented as TF1 object:
      nua.fMaxValuePDF[pdf] = nua.fDefaultNUAPDF[pdf]->GetMaximum(ph.fParticleHistogramsBins[pdf][1], ph.fParticleHistogramsBins[pdf][2]);
    }

  } // for(int pdf=0;pdf<eNUAPDF_N;pdf++) // use pdfs for NUA in (phi, pt, eta, ...).

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookNUAHistograms()

//============================================================

void BookInternalValidationHistograms()
{
  // Book all internal validation histograms.

  // a) Book the profile holding flags;
  // b) Book and fill container for vn amplitudes;
  // c) Book and fill container for Psin planes;
  // d) Handle multiplicity for internal validation.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  iv.fInternalValidationFlagsPro = new TProfile("fInternalValidationFlagsPro", "flags for internal validation", 4, 0., 4.);
  iv.fInternalValidationFlagsPro->SetStats(false);
  iv.fInternalValidationFlagsPro->SetLineColor(eColor);
  iv.fInternalValidationFlagsPro->SetFillColor(eFillColor);
  iv.fInternalValidationFlagsPro->GetXaxis()->SetLabelSize(0.04);
  iv.fInternalValidationList->Add(iv.fInternalValidationFlagsPro);

  iv.fInternalValidationFlagsPro->GetXaxis()->SetBinLabel(1, "fUseInternalValidation");
  iv.fInternalValidationFlagsPro->Fill(0.5, iv.fUseInternalValidation);
  iv.fInternalValidationFlagsPro->GetXaxis()->SetBinLabel(2, "fnEventsInternalValidation");
  iv.fInternalValidationFlagsPro->Fill(1.5, iv.fnEventsInternalValidation);
  iv.fInternalValidationFlagsPro->GetXaxis()->SetBinLabel(3, "fRescaleWithTheoreticalInput");
  iv.fInternalValidationFlagsPro->Fill(2.5, iv.fRescaleWithTheoreticalInput);
  iv.fInternalValidationFlagsPro->GetXaxis()->SetBinLabel(4, TString::Format("option = %s", iv.fHarmonicsOptionInternalValidation->Data()));
  iv.fInternalValidationFlagsPro->Fill(3.5, 1); // redundant, because here I only care about bin label, but preserves symmetry in this code snippet...

  // *) Book object beyond this line only if internal validation was requested:
  if (!iv.fUseInternalValidation) {
    return;
  }

  // b) Book and fill container vn amplitudes:
  iv.fInternalValidationVnPsin[eVn] = new TArrayD(gMaxHarmonic);
  auto lInternalValidationAmplitudes = (vector<float>)cf_iv.cfInternalValidationAmplitudes; // this is now the local version of that array from configurable
  if (lInternalValidationAmplitudes.size() < 1) {
    LOGF(fatal, "\033[1;31m%s at line %d : set at least one vn amplitude in array cfInternalValidationAmplitudes\n \033[0m", __FUNCTION__, __LINE__);
  }
  if (lInternalValidationAmplitudes.size() > gMaxHarmonic) {
    LOGF(fatal, "\033[1;31m%s at line %d : lInternalValidationAmplitudes.size() > gMaxHarmonic \n \033[0m", __FUNCTION__, __LINE__);
  }
  for (int i = 0; i < static_cast<int>(lInternalValidationAmplitudes.size()); i++) {
    iv.fInternalValidationVnPsin[eVn]->SetAt(lInternalValidationAmplitudes[i], i);
  }

  // c) Book and fill container for Psin planes:
  iv.fInternalValidationVnPsin[ePsin] = new TArrayD(gMaxHarmonic);
  auto lInternalValidationPlanes = (vector<float>)cf_iv.cfInternalValidationPlanes;
  if (lInternalValidationPlanes.size() < 1) {
    LOGF(fatal, "\033[1;31m%s at line : %d set at least one Psi plane in array cfInternalValidationPlanes\n \033[0m", __FUNCTION__, __LINE__);
  }
  if (lInternalValidationPlanes.size() > gMaxHarmonic) {
    LOGF(fatal, "\033[1;31m%s at line %d : lInternalValidationPlanes.size() > gMaxHarmonic \n \033[0m", __FUNCTION__, __LINE__);
  }
  if (lInternalValidationAmplitudes.size() != lInternalValidationPlanes.size()) {
    LOGF(fatal, "\033[1;31m%s at line %d : lInternalValidationAmplitudes.size() != lInternalValidationPlanes.size() \n \033[0m", __FUNCTION__, __LINE__);
  }
  for (int i = 0; i < static_cast<int>(lInternalValidationPlanes.size()); i++) {
    iv.fInternalValidationVnPsin[ePsin]->SetAt(lInternalValidationPlanes[i], i);
  }

  // d) Handle multiplicity for internal validation:
  auto lMultRangeInternalValidation = (vector<int>)cf_iv.cfMultRangeInternalValidation;
  iv.fMultRangeInternalValidation[eMin] = lMultRangeInternalValidation[eMin];
  iv.fMultRangeInternalValidation[eMax] = lMultRangeInternalValidation[eMax];
  if (iv.fMultRangeInternalValidation[eMin] >= iv.fMultRangeInternalValidation[eMax]) {
    LOGF(fatal, "\033[1;31m%s at line %d : iv.fMultRangeInternalValidation[eMin] >= iv.fMultRangeInternalValidation[eMax] \n \033[0m", __FUNCTION__, __LINE__);
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // BookInternalValidationHistograms()

//============================================================

TComplex TheoreticalValue(TArrayI* harmonics, TArrayD* amplitudes, TArrayD* planes)
{
  // For the specified harmonics, from available amplitudes and symmetry planes, return the theoretical value of correlator.
  // See Eq. (2) in MVC, originally derived in R. S. Bhalerao, M. Luzum, and J.-Y. Ollitrault, Phys. Rev. C 84, 034910 (2011), arXiv:1104.4740 [nucl-th].

  // a) Insanity checks;
  // b) Main calculus;
  // c) Return value.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Insanity checks:
  if (!harmonics) {
    LOGF(fatal, "\033[1;31m%s at line %d : !harmonics \n \033[0m", __FUNCTION__, __LINE__);
  }
  if (!amplitudes) {
    LOGF(fatal, "\033[1;31m%s at line %d : !amplitudes \n \033[0m", __FUNCTION__, __LINE__);
  }
  if (!planes) {
    LOGF(fatal, "\033[1;31m%s at line %d : !planes \n \033[0m", __FUNCTION__, __LINE__);
  }
  if (amplitudes->GetSize() != planes->GetSize()) {
    LOGF(fatal, "\033[1;31m%s at line %d : amplitudes->GetSize() != planes->GetSize() \n \033[0m", __FUNCTION__, __LINE__);
  }

  // b) Main calculus:
  TComplex value = TComplex(1., 0., true); // yes, polar representation
  for (int h = 0; h < harmonics->GetSize(); h++) {
    //  Using polar form of TComplex (double re, double im=0, bool polar=false):
    value *= TComplex(amplitudes->GetAt(TMath::Abs(harmonics->GetAt(h)) - 1), 1. * harmonics->GetAt(h) * planes->GetAt(TMath::Abs(harmonics->GetAt(h)) - 1), true);
  } // for(int h=0;h<harmonics->GetSize();h++)

  // c) Return value:
  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }
  return value;

} // TComplex TheoreticalValue(TArrayI *harmonics, TArrayD *amplitudes, TArrayD *planes)

//============================================================

void InternalValidation()
{
  // Internal validation against theoretical values in on-the-fly study for all implemented correlators.

  // Last update: 20250121

  // To do:
  // 20250121 At the moment, I do not support here differential phi weights. If I decide to add that feature, basically I need to generalize Accept() for 2D case,
  //          where e.g. phi(pt) weights will be given with some toy 2D pdf.

  // *) Set and propagate some fake run number;
  // *) Fetch the weights for this particular run number. Do it only once;
  // a) Fourier like p.d.f. for azimuthal angles and flow amplitudes;
  // b) Loop over on-the-fly events.
  //    b0) Reset ebye quantities;
  //    b1) Determine multiplicity, centrality, reaction plane and configure p.d.f. for azimuthal angles if harmonics are not constant e-by-e;
  //    b2) Fill event histograms before cuts;
  //    b3) Loop over particles;
  //    b4) Fill event histograms after cuts;
  //    b5) Calculate everything for selected events and particles;
  // c) Delete persistent objects.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Set and propagate some fake run number:
  tc.fRunNumber = "123456";
  PropagateRunNumber();

  // *) Fetch the weights for this particular run number. Do it only once.
  //    TBI 20231012 If eventualy I can access programatically run number in init(...) at run time, this shall go there.
  if (!pw.fParticleWeightsAreFetched) {
    if (pw.fUseWeights[wPHI] || pw.fUseWeights[wPT] || pw.fUseWeights[wETA] || pw.fUseDiffWeights[wPHIPT] || pw.fUseDiffWeights[wPHIETA]) {
      GetParticleWeights();
      pw.fParticleWeightsAreFetched = true;
    }
    // differential phi weights:
    if (pw.fUseDiffPhiWeights[wPhiPhiAxis]) { // Yes, I check only the first flag. This way, I can switch off all differential phi weights by setting 0-wPhi in config.
                                              // On the other hand, it doesn't make sense to calculate differential phi weights without having phi axis.
                                              // At any point I shall be able to fall back to integrated phi weights, that corresponds to the case wheh "1-wPhi" and all others are "0-w..."
      GetParticleWeights();
      pw.fParticleWeightsAreFetched = true;
    }
  } // if (!pw.fParticleWeightsAreFetched) {

  // a) Fourier like p.d.f. for azimuthal angles and flow amplitudes:
  TF1* fPhiPDF = NULL;
  TF3* fvnPDF = NULL;

  if (iv.fHarmonicsOptionInternalValidation->EqualTo("constant")) {
    // For this option, vn's and psin's are constant for all simulated events, therefore I can configure fPhiPDF outside of loop over events.
    // Remark: The last parameter [18] is a random reaction plane, keep in sync with fPhiPDF->SetParameter(18,fReactionPlane); below
    //         Keep also in sync with const int gMaxHarmonic = 9; in *GlobalConstants.h
    fPhiPDF = new TF1("fPhiPDF", "1 + 2.*[0]*TMath::Cos(x-[1]-[18]) + 2.*[2]*TMath::Cos(2.*(x-[3]-[18])) + 2.*[4]*TMath::Cos(3.*(x-[5]-[18])) + 2.*[6]*TMath::Cos(4.*(x-[7]-[18])) + 2.*[8]*TMath::Cos(5.*(x-[9]-[18])) + 2.*[10]*TMath::Cos(6.*(x-[11]-[18])) + 2.*[12]*TMath::Cos(7.*(x-[13]-[18])) + 2.*[14]*TMath::Cos(8.*(x-[15]-[18])) + 2.*[16]*TMath::Cos(9.*(x-[17]-[18]))", 0., o2::constants::math::TwoPI);
    for (int h = 0; h < gMaxHarmonic; h++) {
      fPhiPDF->SetParName(2 * h, TString::Format("v_{%d}", h + 1));       // set name v_n
      fPhiPDF->SetParName(2 * h + 1, TString::Format("Psi_{%d}", h + 1)); // set name psi_n
      // initialize v_n:
      if (iv.fInternalValidationVnPsin[eVn] && h + 1 <= iv.fInternalValidationVnPsin[eVn]->GetSize()) {
        fPhiPDF->SetParameter(2 * h, iv.fInternalValidationVnPsin[eVn]->GetAt(h));
      } else {
        fPhiPDF->SetParameter(2 * h, 0.);
      }
      // initialize psi_n:
      if (iv.fInternalValidationVnPsin[ePsin] && h + 1 <= iv.fInternalValidationVnPsin[ePsin]->GetSize()) {
        fPhiPDF->SetParameter(2 * h + 1, iv.fInternalValidationVnPsin[ePsin]->GetAt(h));
      } else {
        fPhiPDF->SetParameter(2 * h + 1, 0.);
      }
    } // for(int h=0;h<gMaxHarmonic;h++)
    // cross-check set vn's and psin's:

    if (tc.fVerbose) {
      LOGF(info, "=> This is initial configuration for p.d.f. used in internal validation:");
      for (int h = 0; h < 2 * gMaxHarmonic; h++) {
        LOGF(info, Form("%d %s = %f", h, fPhiPDF->GetParName(h), fPhiPDF->GetParameter(h)));
      }
      LOGF(info, "Remark: Parameter [18] at the moment is reaction plane.\n");
    } // if (tc.fVerbose) {

  } else if (iv.fHarmonicsOptionInternalValidation->EqualTo("correlated")) { // if(iv.fHarmonicsOptionInternalValidation->EqualTo("constant"))
    // For this option, three selected vn's (v1,v2,v3) are correlated, and all psin's are set to zero, for simplicity.
    // Remark: The last parameter [3] is a random reaction plane, keep in sync with fPhiPDF->SetParameter(3,fReactionPlane); below
    //         Keep also in sync with const int gMaxHarmonic = 9; in *GlobalConstants.h

    // Azimuthal angles are sampled from this pdf:
    fPhiPDF = new TF1("fPhiPDF", "1 + 2.*[0]*TMath::Cos(x-[3]) + 2.*[1]*TMath::Cos(2.*(x-[3])) + 2.*[2]*TMath::Cos(3.*(x-[3]))", 0., o2::constants::math::TwoPI);
    // With this parameterization, I have:
    //  [0] => v1
    //  [1] => v2
    //  [2] => v3
    //  [3] => RP
    fPhiPDF->SetParName(0, "v_{1}");
    fPhiPDF->SetParName(1, "v_{2}");
    fPhiPDF->SetParName(2, "v_{3}");
    fPhiPDF->SetParName(3, "RP");

    // vn amplitudes are sampled e-b-e from this pdf:
    fvnPDF = new TF3("fvnPDF", "x + 2.*y - 3.*z", 0.07, 0.08, 0.06, 0.07, 0.05, 0.06); // v1 \in [0.07,0.08], v2 \in [0.06,0.07], v3 \in [0.05,0.06]
    // check for example message 'W-TF3::GetRandom3: function:fvnPDF has 27000 negative values: abs assumed' in the log file
    // All the amplitudes v1, v2 and v3, and RP are determined e-b-e, and then set in fPhiPDF below

  } else if (iv.fHarmonicsOptionInternalValidation->EqualTo("persistent")) { // if(iv.fHarmonicsOptionInternalValidation->EqualTo("persistent"))
    // For this option, three selected vn's (v1,v2,v3) are correlated in the same way as in "correlated" case, but in addition, the persistent
    // non-vanishing correlation among SPCs Psi1, Psi2 and Psi3 is introduced, in the same way as in arXiv:1901.06968, Sec. II D.

    // Remark: In this example, there is no Reaction Plane, instead Psi1 and Psi2 are sampled uniformly, and the equation for Psi3 is hardwired,
    //         to introduce strong and persistent SPC correlation, see arXiv:1901.06968, Sec. II D.
    //         Keep also in sync with const int gMaxHarmonic = 9; in *GlobalConstants.h

    // Azimuthal angles are sampled from this pdf:
    fPhiPDF = new TF1("fPhiPDF", "1 + 2.*[0]*TMath::Cos(x-[3]) + 2.*[1]*TMath::Cos(2.*(x-[4])) + 2.*[2]*TMath::Cos(3.*(x-[5]))", 0., o2::constants::math::TwoPI);
    // With this parameterization, I have:
    //  [0] => v1
    //  [1] => v2
    //  [2] => v3
    //  [3] => Psi1
    //  [4] => Psi2
    //  [5] => Psi3
    fPhiPDF->SetParName(0, "v_{1}");
    fPhiPDF->SetParName(1, "v_{2}");
    fPhiPDF->SetParName(2, "v_{3}");
    fPhiPDF->SetParName(3, "Psi_{1}");
    fPhiPDF->SetParName(4, "Psi_{2}");
    fPhiPDF->SetParName(5, "Psi_{3}");

    // vn amplitudes are sampled e-b-e from this pdf (yes, for simplicity, I keep it the same as in "correlated" case):
    fvnPDF = new TF3("fvnPDF", "x + 2.*y - 3.*z", 0.07, 0.08, 0.06, 0.07, 0.05, 0.06); // v1 \in [0.07,0.08], v2 \in [0.06,0.07], v3 \in [0.05,0.06]
    // check for example message 'W-TF3::GetRandom3: function:fvnPDF has 27000 negative values: abs assumed' in the log file
    // All the amplitudes v1, v2 and v3, and symmetry planes Psi_{1}, Psi_{2} and Psi_{3} are determined e-b-e, and then set in fPhiPDF below
  } // else if(fHarmonicsOptionInternalValidation->EqualTo("persistent"))

  // b) Loop over on-the-fly events:
  double v1 = 0., v2 = 0., v3 = 0.;
  for (int e = 0; e < static_cast<int>(iv.fnEventsInternalValidation); e++) {

    // b0) Reset ebye quantities:
    ResetEventByEventQuantities();

    // b1) Determine multiplicity, centrality, reaction plane and configure p.d.f. for azimuthal angles if harmonics are not constant e-by-e:
    int nMult = static_cast<int>(gRandom->Uniform(iv.fMultRangeInternalValidation[eMin], iv.fMultRangeInternalValidation[eMax]));

    double fReactionPlane = gRandom->Uniform(0., o2::constants::math::TwoPI); // no cast is needed, since Uniform(...) returns double
    if (iv.fHarmonicsOptionInternalValidation->EqualTo("constant")) {
      fPhiPDF->SetParameter(18, fReactionPlane);
    } else if (iv.fHarmonicsOptionInternalValidation->EqualTo("correlated")) {
      fPhiPDF->SetParameter(3, fReactionPlane);
    } // Remark: I do not need here anything for option "persistent", because RP is not used for that case. See below how 3 symmetry planes are introduced with persistent correlation

    ebye.fCentrality = static_cast<float>(gRandom->Uniform(0., 100.));           // this is perfectly fine for this exercise
    ebye.fOccupancy = static_cast<float>(gRandom->Uniform(0., 10000.));          // this is perfectly fine for this exercise
    ebye.fInteractionRate = static_cast<float>(gRandom->Uniform(0., 10000.));    // this is perfectly fine for this exercise
    ebye.fCurrentRunDuration = static_cast<float>(gRandom->Uniform(0., 86400.)); // this is perfectly fine for this exercise
    ebye.fVz = static_cast<float>(gRandom->Uniform(-20., 20.));                  // this is perfectly fine for this exercise

    //    b2) Fill event histograms before cuts:
    if (eh.fFillEventHistograms) {
      !eh.fEventHistograms[eNumberOfEvents][eSim][eBefore] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][eBefore]->Fill(0.5);
      !eh.fEventHistograms[eTotalMultiplicity][eSim][eBefore] ? true : eh.fEventHistograms[eTotalMultiplicity][eSim][eBefore]->Fill(nMult);
      !eh.fEventHistograms[eCentrality][eSim][eBefore] ? true : eh.fEventHistograms[eCentrality][eSim][eBefore]->Fill(ebye.fCentrality);
      !eh.fEventHistograms[eOccupancy][eSim][eBefore] ? true : eh.fEventHistograms[eOccupancy][eSim][eBefore]->Fill(ebye.fOccupancy);
      !eh.fEventHistograms[eInteractionRate][eSim][eBefore] ? true : eh.fEventHistograms[eInteractionRate][eSim][eBefore]->Fill(ebye.fInteractionRate);
      !eh.fEventHistograms[eCurrentRunDuration][eSim][eBefore] ? true : eh.fEventHistograms[eCurrentRunDuration][eSim][eBefore]->Fill(ebye.fCurrentRunDuration);
      !eh.fEventHistograms[eVertex_z][eSim][eBefore] ? true : eh.fEventHistograms[eVertex_z][eSim][eBefore]->Fill(ebye.fVz);
      !eh.fEventHistograms[eEventPlaneAngle][eSim][eBefore] ? true : eh.fEventHistograms[eEventPlaneAngle][eSim][eBefore]->Fill(fReactionPlane);
    }

    // ... here I could implement some event cuts, if necessary ...

    // configure p.d.f. for azimuthal angles if harmonics are not constant e-by-e, for option "correlated":
    if (iv.fHarmonicsOptionInternalValidation->EqualTo("correlated")) {
      // Sample 3 correlated vn's from TF3 fvnPDF, and with them initialize fPhiPDF:
      fvnPDF->GetRandom3(v1, v2, v3);
      fPhiPDF->SetParameter(0, v1);
      fPhiPDF->SetParameter(1, v2);
      fPhiPDF->SetParameter(2, v3);
      // reaction plane is set above already
    } // if(fHarmonicsOptionInternalValidation->EqualTo("correlated"))

    // configure p.d.f. for azimuthal angles if harmonics are not constant e-by-e, for option "persistent":
    if (iv.fHarmonicsOptionInternalValidation->EqualTo("persistent")) {

      // Sample 3 correlated vn's from TF3 fvnPDF, and with them initialize fPhiPDF:
      fvnPDF->GetRandom3(v1, v2, v3);
      fPhiPDF->SetParameter(0, v1);
      fPhiPDF->SetParameter(1, v2);
      fPhiPDF->SetParameter(2, v3);

      // Persistent symmetry plane correlation:
      double Psi1 = gRandom->Uniform(0., o2::constants::math::TwoPI);
      double Psi2 = gRandom->Uniform(0., o2::constants::math::TwoPI);
      double Psi3 = (1. / 3.) * ((o2::constants::math::PI / 4.) + 2. * Psi2 + Psi1); // see arXiv:1901.06968, Sec. II D.
      fPhiPDF->SetParameter(3, Psi1);
      fPhiPDF->SetParameter(4, Psi2);
      fPhiPDF->SetParameter(5, Psi3);

      // Remark: reaction plane is not needed for case "persistent"

    } // if(fHarmonicsOptionInternalValidation->EqualTo("persistent"))

    // b3) Loop over particles:
    double dPhi = 0.;
    double dPt = 0.;
    double dEta = 0.;

    // *) Define min and max ranges for sampling:
    double dPt_min = res.fResultsPro[AFO_PT]->GetXaxis()->GetBinLowEdge(1);                                           // yes, low edge of first bin is pt min
    double dPt_max = res.fResultsPro[AFO_PT]->GetXaxis()->GetBinLowEdge(1 + res.fResultsPro[AFO_PT]->GetNbinsX());    // yes, low edge of overflow bin is max pt
    double dEta_min = res.fResultsPro[AFO_ETA]->GetXaxis()->GetBinLowEdge(1);                                         // yes, low edge of first bin is eta min
    double dEta_max = res.fResultsPro[AFO_ETA]->GetXaxis()->GetBinLowEdge(1 + res.fResultsPro[AFO_ETA]->GetNbinsX()); // yes, low edge of overflow bin is max eta

    for (int p = 0; p < nMult; p++) {
      // Particle angle:
      dPhi = fPhiPDF->GetRandom();

      // *) To increase performance, sample pt or eta only if requested:
      if (mupa.fCalculateCorrelationsAsFunctionOf[AFO_PT] || t0.fCalculateTest0AsFunctionOf[AFO_PT] || es.fCalculateEtaSeparationsAsFunctionOf[AFO_PT]) {
        dPt = gRandom->Uniform(dPt_min, dPt_max);
      }

      if (mupa.fCalculateCorrelationsAsFunctionOf[AFO_ETA] || t0.fCalculateTest0AsFunctionOf[AFO_ETA] || es.fCalculateEtaSeparations) {
        // Yes, I have to use here es.fCalculateEtaSeparations , and not some differential flag, like for pt case above
        dEta = gRandom->Uniform(dEta_min, dEta_max);
      }

      // *) Fill few selected particle histograms before cuts here directly:
      // Remark: I do not call FillParticleHistograms<rs>(track, eBefore), as I do not want to bother to make here full 'track' object, etc., just to fill simple kine info:
      if (ph.fFillParticleHistograms || ph.fFillParticleHistograms2D) {
        // 1D:
        !ph.fParticleHistograms[ePhi][eSim][eBefore] ? true : ph.fParticleHistograms[ePhi][eSim][eBefore]->Fill(dPhi);
        !ph.fParticleHistograms[ePt][eSim][eBefore] ? true : ph.fParticleHistograms[ePt][eSim][eBefore]->Fill(dPt);
        !ph.fParticleHistograms[eEta][eSim][eBefore] ? true : ph.fParticleHistograms[eEta][eSim][eBefore]->Fill(dEta);
        // 2D:
        !ph.fParticleHistograms2D[ePhiPt][eSim][eBefore] ? true : ph.fParticleHistograms2D[ePhiPt][eSim][eBefore]->Fill(dPhi, dPt);
        !ph.fParticleHistograms2D[ePhiEta][eSim][eBefore] ? true : ph.fParticleHistograms2D[ePhiEta][eSim][eBefore]->Fill(dPhi, dEta);
      }

      // *) Particle cuts (only support for Toy NUA is provided, for the time being):
      //    NUA:
      if (nua.fApplyNUAPDF[ePhiNUAPDF] && !Accept(dPhi, ePhiNUAPDF)) {
        continue;
      }
      if (nua.fApplyNUAPDF[ePtNUAPDF] && !Accept(dPt, ePtNUAPDF)) {
        continue;
      }
      if (nua.fApplyNUAPDF[eEtaNUAPDF] && !Accept(dEta, eEtaNUAPDF)) {
        continue;
      }

      // *) Fill few selected particle histograms after cuts here directly here:
      // Remark: I do not call FillParticleHistograms<rs>(track, eAfter), as I do not want to bother to make here full 'track' object, etc., just to fill simple kine info:
      if (ph.fFillParticleHistograms || ph.fFillParticleHistograms2D) {
        // 1D:
        !ph.fParticleHistograms[ePhi][eSim][eAfter] ? true : ph.fParticleHistograms[ePhi][eSim][eAfter]->Fill(dPhi);
        !ph.fParticleHistograms[ePt][eSim][eAfter] ? true : ph.fParticleHistograms[ePt][eSim][eAfter]->Fill(dPt);
        !ph.fParticleHistograms[eEta][eSim][eAfter] ? true : ph.fParticleHistograms[eEta][eSim][eAfter]->Fill(dEta);
        // 2D:
        !ph.fParticleHistograms2D[ePhiPt][eSim][eAfter] ? true : ph.fParticleHistograms2D[ePhiPt][eSim][eAfter]->Fill(dPhi, dPt);
        !ph.fParticleHistograms2D[ePhiEta][eSim][eAfter] ? true : ph.fParticleHistograms2D[ePhiEta][eSim][eAfter]->Fill(dPhi, dEta);
      }

      // Remark: Keep in sync all calls and flags below with the ones in MainLoopOverParticles().
      // *) Integrated Q-vectors:
      if (qv.fCalculateQvectors || es.fCalculateEtaSeparations) {
        this->FillQvector(dPhi, dPt, dEta); // all 3 arguments are passed by reference
      }

      // *) Differential q-vectors:
      // **) pt-dependence:
      if (qv.fCalculateQvectors && (mupa.fCalculateCorrelationsAsFunctionOf[AFO_PT] || t0.fCalculateTest0AsFunctionOf[AFO_PT]) && !es.fCalculateEtaSeparations) {
        // In this branch I do not need eta separation, so the lighter call can be executed:
        this->Fillqvector(dPhi, dPt, PTq); // first 2 arguments are passed by reference, 3rd argument is enum
      } else if (es.fCalculateEtaSeparations && es.fCalculateEtaSeparationsAsFunctionOf[AFO_PT]) {
        // In this branch I do need eta separation, so the heavier call must be executed:
        // Remark: Within Fillqvector() I check again all the relevant flags.
        this->Fillqvector(dPhi, dPt, PTq, dEta); // first 2 arguments and the last one are passed by reference, 3rd argument is enum. "kine" variable is the 2nd argument
      }
      // **) eta-dependence:
      if (qv.fCalculateQvectors && (mupa.fCalculateCorrelationsAsFunctionOf[AFO_ETA] || t0.fCalculateTest0AsFunctionOf[AFO_ETA])) {
        // Remark: For eta dependence I do not consider es.fCalculateEtaSeparations, because in this context that calculation is meaningless.
        this->Fillqvector(dPhi, dEta, ETAq); // first 2 arguments are passed by reference, 3rd argument is enum
      }

      // *) Fill nested loops containers:
      if (nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoops) {
        this->FillNestedLoopsContainers(ebye.fSelectedTracks, dPhi, dPt, dEta); // all 4 arguments are passed by reference
      }

      // *) Counter of selected tracks in the current event:
      //    Remark: This has to go after FillNestedLoopsContainers(...), because ebye.fSelectedTracks is used as a particle index there.
      ebye.fSelectedTracks++;
      if (ebye.fSelectedTracks >= ec.fdEventCuts[eMultiplicity][eMax]) {
        break;
      }

    } // for(int p=0;p<nMult;p++)

    // *) Determine multiplicity of this event, for all "vs. mult" results:
    DetermineMultiplicity();

    //    b4) Fill event histograms after all event and particle cuts:
    if (eh.fFillEventHistograms) {
      !eh.fEventHistograms[eNumberOfEvents][eSim][eAfter] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]->Fill(0.5);
      !eh.fEventHistograms[eTotalMultiplicity][eSim][eAfter] ? true : eh.fEventHistograms[eTotalMultiplicity][eSim][eAfter]->Fill(nMult);
      !eh.fEventHistograms[eMultiplicity][eSim][eAfter] ? true : eh.fEventHistograms[eMultiplicity][eSim][eAfter]->Fill(ebye.fMultiplicity);
      !eh.fEventHistograms[eCentrality][eSim][eAfter] ? true : eh.fEventHistograms[eCentrality][eSim][eAfter]->Fill(ebye.fCentrality);
      !eh.fEventHistograms[eOccupancy][eSim][eAfter] ? true : eh.fEventHistograms[eOccupancy][eSim][eAfter]->Fill(ebye.fOccupancy);
      !eh.fEventHistograms[eInteractionRate][eSim][eAfter] ? true : eh.fEventHistograms[eCentrality][eSim][eAfter]->Fill(ebye.fInteractionRate);
      !eh.fEventHistograms[eCurrentRunDuration][eSim][eAfter] ? true : eh.fEventHistograms[eCurrentRunDuration][eSim][eAfter]->Fill(ebye.fCurrentRunDuration);
      !eh.fEventHistograms[eVertex_z][eSim][eAfter] ? true : eh.fEventHistograms[eVertex_z][eSim][eAfter]->Fill(ebye.fVz);
      !eh.fEventHistograms[eEventPlaneAngle][eSim][eAfter] ? true : eh.fEventHistograms[eEventPlaneAngle][eSim][eAfter]->Fill(fReactionPlane);
    }

    // *) Fill subevent multiplicities:
    //    Remark: I can call this one only after Qa and Qb vectors are filled:
    if (es.fCalculateEtaSeparations) {
      FillSubeventMultiplicities<eSim>();
    }

    // b5) Calculate everything for selected events and particles:
    CalculateEverything();

    // *) Reset event-by-event quantities:
    ResetEventByEventQuantities();

    // *) Print info on the current event number (within current real event):
    LOGF(info, "   Event # %d/%d (within current real event) ....", e + 1, static_cast<int>(iv.fnEventsInternalValidation));

    // *) Determine all event counters:
    DetermineEventCounters();

    // *) Sequential bailout: After each tc.fSequentialBailout events, I bail out:
    if (iv.fInternalValidationForceBailout && tc.fSequentialBailout > 0 && eh.fEventCounter[eProcessed] > 0 && 0 == eh.fEventCounter[eProcessed] % tc.fSequentialBailout) {
      BailOut();
    }

    // *) If I reached max number of events, ignore the remaining collisions:
    if (MaxNumberOfEvents(eAfter)) {
      if (iv.fInternalValidationForceBailout) {
        BailOut(true);
      }
    }

  } // for(int e=0;e<e<static_cast<int>(iv.fnEventsInternalValidation);e++)

  // *) Print info on the current event number (total):
  if (tc.fVerboseEventCounter) {
    PrintEventCounter(eAfter);
  }

  // c) Delete persistent objects:
  if (fPhiPDF) {
    delete fPhiPDF;
  }
  if (fvnPDF) {
    delete fvnPDF;
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void InternalValidation()

//============================================================

bool Accept(const double& value, int var)
{
  // Given the acceptance profile for this observable, accept or not that observable for the analysis.
  // Use in Toy NUA studies.

  // Remark: var corresponds to the field in enum eNUAPDF { ePhiNUAPDF, ePtNUAPDF, eEtaNUAPDF };
  //         Therefore, always call this function as e.g. Accept(someAngle, ePhiNUAPDF) or Accept(somePt, ePtNUAPDF)

  if (tc.fVerboseForEachParticle) {
    LOGF(info, "\033[1;32m%s\033[0m", __FUNCTION__);
  }

  // Basic protection:
  if (nua.fUseDefaultNUAPDF[var] && !nua.fDefaultNUAPDF[var]) {
    LOGF(info, "\033[1;33m%s var = %d\033[0m", static_cast<int>(var));
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  } else if (!nua.fUseDefaultNUAPDF[var] && !nua.fCustomNUAPDF[var]) {
    LOGF(info, "\033[1;33m%s var = %d\033[0m", static_cast<int>(var));
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  bool bAccept = true; // return value

  double acceptanceProbability = 1.;
  double correspondingAcceptance = -44.;
  if (!nua.fUseDefaultNUAPDF[var]) {
    correspondingAcceptance = nua.fCustomNUAPDF[var]->GetBinContent(nua.fCustomNUAPDF[var]->FindBin(value));
  } else {
    correspondingAcceptance = nua.fDefaultNUAPDF[var]->Eval(value);
  }

  // Probability to accept:
  acceptanceProbability = 1. - (nua.fMaxValuePDF[var] - correspondingAcceptance) / nua.fMaxValuePDF[var];

  // Accept or not:
  (gRandom->Uniform(0, 1) < acceptanceProbability) ? bAccept = true : bAccept = false;

  return bAccept;

} // bool Accept(const double &value, int var)

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
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  t0.fTest0FlagsPro = new TProfile("fTest0FlagsPro", "flags for Test0", 1, 0., 1.);
  t0.fTest0FlagsPro->SetStats(false);
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
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // d) Book what needs to be booked:
  for (int mo = 0; mo < gMaxCorrelator; mo++) {
    for (int mi = 0; mi < gMaxIndex; mi++) {
      if (!t0.fTest0Labels[mo][mi]) {
        continue;
      }
      {
        for (int v = 0; v < eAsFunctionOf_N; v++) {
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
          if (AFO_OCCUPANCY == v && !t0.fCalculateTest0AsFunctionOf[AFO_OCCUPANCY]) {
            continue;
          }
          if (AFO_INTERACTIONRATE == v && !t0.fCalculateTest0AsFunctionOf[AFO_INTERACTIONRATE]) {
            continue;
          }
          if (AFO_CURRENTRUNDURATION == v && !t0.fCalculateTest0AsFunctionOf[AFO_CURRENTRUNDURATION]) {
            continue;
          }
          if (AFO_VZ == v && !t0.fCalculateTest0AsFunctionOf[AFO_VZ]) {
            continue;
          }

          if (!res.fResultsPro[v]) {
            LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
          }

          t0.fTest0Pro[mo][mi][v] = reinterpret_cast<TProfile*>(res.fResultsPro[v]->Clone(Form("fTest0Pro[%d][%d][%s]", mo, mi, res.fResultsProRawName[v].Data()))); // yes
          t0.fTest0Pro[mo][mi][v]->SetStats(false);
          t0.fTest0Pro[mo][mi][v]->Sumw2();
          t0.fTest0Pro[mo][mi][v]->SetTitle(t0.fTest0Labels[mo][mi]->Data());
          t0.fTest0Pro[mo][mi][v]->GetXaxis()->SetTitle(FancyFormatting(res.fResultsProXaxisTitle[v].Data()));
          /*
                if(fUseFixedNumberOfRandomlySelectedParticles && 1==v) // just a warning for the meaning of multiplicity in this special case
                {
                 fTest0Pro[mo][mi][1]->GetXaxis()->SetTitle("WARNING: for each multiplicity, fFixedNumberOfRandomlySelectedParticles is selected randomly in Q-vector");
                }
          */
          t0.fTest0List->Add(t0.fTest0Pro[mo][mi][v]); // yes, this has to be here
        } // for(int v=0;v<eAsFunctionOf_N;v++) // variable, see content of enum eAsFunctionOf
      } // if(fTest0Labels[mo][mi])
    } // for(int mi=0;mi<gMaxIndex;mi++)
  } // for(int mo=0;mo<gMaxCorrelator;mo++)

  // e) Few quick insanity checks on booking:
  if (t0.fTest0Pro[0][0][AFO_INTEGRATED] && !TString(t0.fTest0Pro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }
  if (t0.fTest0Pro[0][0][AFO_PT] && !TString(t0.fTest0Pro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("p_{T}")) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookTest0Histograms()

//============================================================

void BookEtaSeparationsHistograms()
{
  // Book all eta separations histograms.

  // a) Book the profile holding flags;
  // b) Book what needs to be booked;
  // c) Few quick insanity checks on booking.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  es.fEtaSeparationsFlagsPro = new TProfile("fEtaSeparationsFlagsPro", "flags for eta separations", 1, 0., 1.);
  es.fEtaSeparationsFlagsPro->SetStats(false);
  es.fEtaSeparationsFlagsPro->SetLineColor(eColor);
  es.fEtaSeparationsFlagsPro->SetFillColor(eFillColor);
  es.fEtaSeparationsFlagsPro->GetXaxis()->SetLabelSize(0.04);
  es.fEtaSeparationsFlagsPro->GetXaxis()->SetBinLabel(1, "fCalculateEtaSeparations");
  es.fEtaSeparationsFlagsPro->Fill(0.5, es.fCalculateEtaSeparations);
  es.fEtaSeparationsList->Add(es.fEtaSeparationsFlagsPro);

  if (!es.fCalculateEtaSeparations) {
    return;
  }

  // b) Book what needs to be booked:
  for (int h = 0; h < gMaxHarmonic; h++) {
    if (es.fEtaSeparationsSkipHarmonics[h]) {
      continue;
    }
    for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
      for (int v = 0; v < eAsFunctionOf_N; v++) {
        // decide what is booked, then later valid pointer to es.fEtaSeparationsPro[h][e][v] is used as a boolean, in the standard way:
        if (AFO_INTEGRATED == v && !es.fCalculateEtaSeparationsAsFunctionOf[AFO_INTEGRATED]) {
          continue;
        }
        if (AFO_MULTIPLICITY == v && !es.fCalculateEtaSeparationsAsFunctionOf[AFO_MULTIPLICITY]) {
          continue;
        }
        if (AFO_CENTRALITY == v && !es.fCalculateEtaSeparationsAsFunctionOf[AFO_CENTRALITY]) {
          continue;
        }
        if (AFO_PT == v && !es.fCalculateEtaSeparationsAsFunctionOf[AFO_PT]) {
          continue;
        }
        if (AFO_ETA == v && !es.fCalculateEtaSeparationsAsFunctionOf[AFO_ETA]) {
          continue;
        }
        if (AFO_OCCUPANCY == v && !es.fCalculateEtaSeparationsAsFunctionOf[AFO_OCCUPANCY]) {
          continue;
        }
        if (AFO_INTERACTIONRATE == v && !es.fCalculateEtaSeparationsAsFunctionOf[AFO_INTERACTIONRATE]) {
          continue;
        }
        if (AFO_CURRENTRUNDURATION == v && !es.fCalculateEtaSeparationsAsFunctionOf[AFO_CURRENTRUNDURATION]) {
          continue;
        }
        if (AFO_VZ == v && !es.fCalculateEtaSeparationsAsFunctionOf[AFO_VZ]) {
          continue;
        }

        if (!res.fResultsPro[v]) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }

        es.fEtaSeparationsPro[h][e][v] = reinterpret_cast<TProfile*>(res.fResultsPro[v]->Clone(Form("fEtaSeparationsPro[%d][%d][%s]", h, e, res.fResultsProRawName[v].Data()))); // yes
        es.fEtaSeparationsPro[h][e][v]->SetStats(false);
        es.fEtaSeparationsPro[h][e][v]->Sumw2();
        es.fEtaSeparationsPro[h][e][v]->SetTitle(Form("%d -%d, |#Delta#eta| > %.2f", h + 1, h + 1, es.fEtaSeparationsValues[e]));
        es.fEtaSeparationsPro[h][e][v]->GetXaxis()->SetTitle(FancyFormatting(res.fResultsProXaxisTitle[v].Data()));
        es.fEtaSeparationsList->Add(es.fEtaSeparationsPro[h][e][v]); // yes, this has to be here
      } // for(int v=0;v<eAsFunctionOf_N;v++) // variable, see content of enum eAsFunctionOf
    } // for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
  } // for (int h = 0; h < gMaxHarmonic; h++) {

  // c) Few quick insanity checks on booking:
  if (es.fEtaSeparationsPro[0][0][AFO_INTEGRATED] && !TString(es.fEtaSeparationsPro[0][0][AFO_INTEGRATED]->GetXaxis()->GetTitle()).EqualTo("integrated")) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }
  if (es.fEtaSeparationsPro[0][0][AFO_PT] && !TString(es.fEtaSeparationsPro[0][0][AFO_PT]->GetXaxis()->GetTitle()).EqualTo("p_{T}")) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__); // ordering in enum eAsFunctionOf is not the same as in TString fResultsProXaxisTitle[eAsFunctionOf_N]
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookEtaSeparationsHistograms()

//============================================================

void BookResultsHistograms()
{
  // Book all results histograms.

  // a) Book the profile holding flags;
  // b) Book results histograms, which in addition act as a sort of "abstract" interface, which defines common binning, etc., for other groups of histograms.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the profile holding flags:
  res.fResultsFlagsPro = new TProfile("fResultsFlagsPro", "flags for results histograms", 1, 0., 1.);
  res.fResultsFlagsPro->SetStats(false);
  res.fResultsFlagsPro->SetLineColor(eColor);
  res.fResultsFlagsPro->SetFillColor(eFillColor);
  res.fResultsFlagsPro->GetXaxis()->SetBinLabel(1, "fSaveResultsHistograms");
  res.fResultsFlagsPro->Fill(0.5, res.fSaveResultsHistograms);
  // ...
  res.fResultsList->Add(res.fResultsFlagsPro);

  // b) Book results histograms, which in addition act as a sort of "abstract" interface, which defines common binning, etc., for other groups of histograms:
  for (int v = 0; v < eAsFunctionOf_N; v++) {
    if (res.fUseResultsProVariableLengthBins[v]) {
      // per demand, variable-length binning:
      res.fResultsPro[v] = new TProfile(Form("fResultsPro[%s]", res.fResultsProRawName[v].Data()), "...", res.fResultsProVariableLengthBins[v]->GetSize() - 1, res.fResultsProVariableLengthBins[v]->GetArray());
    } else {
      // the default fixed-length binning:
      res.fResultsPro[v] = new TProfile(Form("fResultsPro[%s]", res.fResultsProRawName[v].Data()), "...", static_cast<int>(res.fResultsProFixedLengthBins[v][0]), res.fResultsProFixedLengthBins[v][1], res.fResultsProFixedLengthBins[v][2]);
    }

    // Optionally, save these histograms. Or just use them as an "abstract" interface for the booking of other group of histograms:
    if (res.fSaveResultsHistograms) {
      res.fResultsList->Add(res.fResultsPro[v]);
    }
  } // for (int v = 0; v < eAsFunctionOf_N; v++) {

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookResultsHistograms()

//============================================================

void BookTheRest()
{
  // Here I book everything not sorted (yes) in specific functions above.

  // a) Book the timer;
  // *) ...

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Book the timer:
  if (tc.fUseStopwatch) {
    tc.fTimer[eGlobal] = new TStopwatch();
    tc.fTimer[eGlobal]->Start();
    tc.fTimer[eLocal] = new TStopwatch();
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void BookTheRest()

//============================================================

template <eRecSim rs, typename T1, typename T2>
void Preprocess(T1 const& collision, T2 const& bcs)
{
  // Do all thingies before starting to process data (e.g. count number of events, fetch the run number, get the weights for this run number, etc.).

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Determine all event counters:
  DetermineEventCounters();

  // *) Sequential bailout: After each tc.fSequentialBailout events, I bail out:
  if (tc.fSequentialBailout > 0 && eh.fEventCounter[eProcessed] > 0 && 0 == eh.fEventCounter[eProcessed] % tc.fSequentialBailout) {
    BailOut();
  }

  // *) If I reached max number of events, ignore the remaining collisions:
  if (MaxNumberOfEvents(eAfter) || MaxNumberOfEvents(eBefore)) { // TBI 20240510 this is a bit confusing, implemented this way. Shall I split off?
    BailOut(true);
  }

  // *) Determine and propagate run number info to already booked objects:
  if (!tc.fRunNumberIsDetermined) {
    DetermineRunNumber<rs>(collision, bcs);
    PropagateRunNumber();
  }
  if (tc.fDoAdditionalInsanityChecks && tc.fRunNumberIsDetermined) {
    CheckCurrentRunNumber<rs>(collision, bcs);
  }

  // *) Fetch the weights for this particular run number. Do it only once.
  //    TBI 20231012 If eventualy I can access programatically run number in init(...) at run time, this shall go there.
  if (!pw.fParticleWeightsAreFetched) {

    // integrated weights and differentials weights without sparse histograms (the latter is becoming obsolete):
    if (pw.fUseWeights[wPHI] || pw.fUseWeights[wPT] || pw.fUseWeights[wETA] || pw.fUseDiffWeights[wPHIPT] || pw.fUseDiffWeights[wPHIETA]) {
      pw.fParticleWeightsAreFetched = true;
    }

    // differential particle weights using sparse histogreams:
    if (pw.fUseDiffPhiWeights[wPhiPhiAxis] || pw.fUseDiffPtWeights[wPtPtAxis] || pw.fUseDiffPtWeights[wEtaEtaAxis]) {
      // Yes, I check only the first flag. This way, I can e.g. switch off all differential phi weights by setting 0-wPhi in config.
      // On the other hand, it doesn't make sense to calculate differential phi weights without having phi axis.
      // At any point I shall be able to fall back e.g. to integrated phi weights, that corresponds to the case wheh "1-wPhi" and all others are "0-w..."
      // Same for differential pt or eta weights.
      GetParticleWeights();
      pw.fParticleWeightsAreFetched = true;
    }

  } // if (!pw.fParticleWeightsAreFetched) {

  // *) Fetch the centrality weights for this particular run number. Do it only once.
  //    TBI 20231012 If eventualy I can access programatically run number in init(...) at run time, this shall go there.
  if (!cw.fCentralityWeightsAreFetched) {
    if (cw.fUseCentralityWeights) {
      GetCentralityWeights();
      cw.fCentralityWeightsAreFetched = true;
    }
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T1, typename T2> void Preprocess(T1 const& collision, T2 const& bcs)

//============================================================

template <eRecSim rs, typename T1, typename T2>
void DetermineRunNumber(T1 const& collision, T2 const&)
{
  // Determine run number and all related thingies.
  // Make sure in process(...) that this function is called only once.

  // TBI 20231018 At the moment I can access run number info only in process(...), but not in init(...)
  // Once I can access run number info in init(...), this function shall be called in init(...), not in process(...)

  // a) Determine run number for Run 3 and Run 2 real data;
  // b) Determine run number for the rest. TBI 20241126 differentiate this support as well, e.g. for eRecSim and eSim.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Determine run number for Run 3 real data:
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2) {

    // **) Determine run number:
    // Get start timestamp and end timemstamp for this run in miliseconds, and convert both of them in seconds:
    // o see O2/CCDB/src/BasicCCDBManager.cxx, O2/CCDB/include/CCDB/BasicCCDBManager.h
    // o example usage in O2Physics/PWGLF/TableProducer/Common/zdcSP.cxx
    auto bc = collision.template foundBC_as<T2>(); // I have the same code snippet at other places, keep in sync.
    tc.fRunNumber = Form("%d", bc.runNumber());
    if (tc.fRunNumber.EqualTo("")) {
      LOGF(error, "\033[1;33m%fRunNumber is empty, bc.runNumber() failed...\033[0m");
      LOGF(fatal, "\033[1;31m%s at line %d : bc.runNumber() = %d \033[0m", __FUNCTION__, __LINE__, bc.runNumber());
    }

    // **) Determine SoR, EoR, and run duration:
    auto runDuration = ccdb->getRunDuration(bc.runNumber());                         // this is total run duration, not the current one (see below)
    tc.fRunTime[eStartOfRun] = std::floor(runDuration.first * 0.001);                // in seconds since Unix epoch
    tc.fRunTime[eEndOfRun] = std::ceil(runDuration.second * 0.001);                  // in seconds since Unix epoch
    tc.fRunTime[eDurationInSec] = tc.fRunTime[eEndOfRun] - tc.fRunTime[eStartOfRun]; // yes, this is now run duration in seconds

    if (!(tc.fRunTime[eStartOfRun] > 0)) {
      LOGF(fatal, "\033[1;31m%s at line %d : tc.fRunTime[eStartOfRun] = %d is not positive\033[0m", __FUNCTION__, __LINE__, tc.fRunTime[eStartOfRun]);
    }
    if (!(tc.fRunTime[eEndOfRun] > 0)) {
      LOGF(fatal, "\033[1;31m%s at line %d : tc.fRunTime[eEndOfRun] = %d is not positive\033[0m", __FUNCTION__, __LINE__, tc.fRunTime[eEndOfRun]);
    }
    if (!(tc.fRunTime[eDurationInSec] > 0)) {
      LOGF(fatal, "\033[1;31m%s at line %d : tc.fRunTime[eDurationInSec] = %d is not positive\033[0m", __FUNCTION__, __LINE__, tc.fRunTime[eDurationInSec]);
    }

  } else {
    // b) Determine run number for the rest. TBI 20241126 differentiate this support as well, e.g. for eRecSim and eSim.
    LOGF(fatal, "\033[1;31m%s at line %d : bc.runNumber() is not validated yet for this case\033[0m", __FUNCTION__, __LINE__);
  }
  tc.fRunNumberIsDetermined = true;

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T1, typename T2> void DetermineRunNumber(T1 const& collision, T2 const&)

//============================================================

void PropagateRunNumber()
{
  // Propagate run number info to already booked objects, wherever it's relevant.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Do some local insanity checks:
  if (tc.fRunNumber.EqualTo("")) {
    LOGF(fatal, "\033[1;31m%s at line %d : tc.fRunNumber is empty \033[0m", __FUNCTION__, __LINE__);
  }

  // *) base:
  fBasePro->GetXaxis()->SetBinLabel(eRunNumber, Form("fRunNumber = %s", tc.fRunNumber.Data()));

  // *) common title var:
  TString histTitle = "";

  // *) event cuts:
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {
    for (int cc = 0; cc < eCutCounter_N; cc++) // enum eCutCounter
    {
      if (!ec.fEventCutCounterHist[rs][cc]) {
        continue;
      }
      histTitle = ec.fEventCutCounterHist[rs][cc]->GetTitle();
      if (histTitle.Contains("__RUN_NUMBER__")) {
        histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
        ec.fEventCutCounterHist[rs][cc]->SetTitle(histTitle.Data());
      }
    }
  }

  // *) event histograms 1D:
  for (int t = 0; t < eEventHistograms_N; t++) // type, see enum eEventHistograms
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!eh.fEventHistograms[t][rs][ba]) {
          continue;
        }
        histTitle = eh.fEventHistograms[t][rs][ba]->GetTitle();
        if (histTitle.Contains("__RUN_NUMBER__")) {
          histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
          eh.fEventHistograms[t][rs][ba]->SetTitle(histTitle.Data());
        }
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for(int t=0;t<eEventHistograms_N;t++) // type, see enum        // eEventHistograms

  // *) event histograms 2D:
  for (int t = 0; t < eQAEventHistograms2D_N; t++) // type, see enum eEventHistograms2D
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!qa.fQAEventHistograms2D[t][rs][ba]) {
          continue;
        }
        histTitle = qa.fQAEventHistograms2D[t][rs][ba]->GetTitle();
        if (histTitle.Contains("__RUN_NUMBER__")) {
          histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
          qa.fQAEventHistograms2D[t][rs][ba]->SetTitle(histTitle.Data());
        }
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for (int t = 0; t < eQAEventHistograms2D_N; t++) // type, see enum eEventHistograms2D

  // *) particle histograms 2D:
  for (int t = 0; t < eQAParticleHistograms2D_N; t++) // type, see enum eParticleHistograms2D
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!qa.fQAParticleHistograms2D[t][rs][ba]) {
          continue;
        }
        histTitle = qa.fQAParticleHistograms2D[t][rs][ba]->GetTitle();
        if (histTitle.Contains("__RUN_NUMBER__")) {
          histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
          qa.fQAParticleHistograms2D[t][rs][ba]->SetTitle(histTitle.Data());
        }
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for (int t = 0; t < eQAParticleHistograms2D_N; t++) // type, see enum eParticleHistograms2D

  // *) particle event histograms 2D:
  for (int t = 0; t < eQAParticleEventHistograms2D_N; t++) // type, see enum eParticleEventHistograms2D
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!qa.fQAParticleEventHistograms2D[t][rs][ba]) {
          continue;
        }
        histTitle = qa.fQAParticleEventHistograms2D[t][rs][ba]->GetTitle();
        if (histTitle.Contains("__RUN_NUMBER__")) {
          histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
          qa.fQAParticleEventHistograms2D[t][rs][ba]->SetTitle(histTitle.Data());
        }
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for (int t = 0; t < eQAParticleEventHistograms2D_N; t++) // type, see enum eParticleEventHistograms2D

  // *) particle sparse histograms:
  for (int t = 0; t < eDiffWeightCategory_N; t++) // category, see enum eDiffWeightCategory
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      if (!ph.fParticleSparseHistograms[t][rs]) {
        continue;
      }
      histTitle = ph.fParticleSparseHistograms[t][rs]->GetTitle();
      if (histTitle.Contains("__RUN_NUMBER__")) {
        histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
        ph.fParticleSparseHistograms[t][rs]->SetTitle(histTitle.Data());
      }
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for (int t = 0; t < eDiffWeightCategory; t++) // category, see enum eDiffWeightCategory

  // *) "correlations vs." histograms 2D:
  for (int t = 0; t < eQACorrelationsVsHistograms2D_N; t++) // type, see enum eCorrelationsVsHistograms2D
  {
    for (int h = 0; h < gMaxHarmonic; h++) {
      for (int rs = 0; rs < 2; rs++) // reco/sim
      {
        if (!qa.fQACorrelationsVsHistograms2D[t][h][rs]) {
          continue;
        }
        histTitle = qa.fQACorrelationsVsHistograms2D[t][h][rs]->GetTitle();
        if (histTitle.Contains("__RUN_NUMBER__")) {
          histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
          qa.fQACorrelationsVsHistograms2D[t][h][rs]->SetTitle(histTitle.Data());
        }
      } // for(int rs=0;rs<2;rs++) // reco/sim
    } // for (int h = 0; h < gMaxHarmonic; h++)
  } // for (int t = 0; t < eQACorrelationsVsHistograms2D_N; t++) // type, see enum eCorrelationsVsHistograms2D

  // *) particle cuts:
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {
    for (int cc = 0; cc < eCutCounter_N; cc++) // enum eCutCounter
    {
      if (!pc.fParticleCutCounterHist[rs][cc]) {
        continue;
      }
      histTitle = pc.fParticleCutCounterHist[rs][cc]->GetTitle();
      if (histTitle.Contains("__RUN_NUMBER__")) {
        histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
        pc.fParticleCutCounterHist[rs][cc]->SetTitle(histTitle.Data());
      }
    }
  }

  // *) particle histograms 1D:
  for (int t = 0; t < eParticleHistograms_N; t++) // type, see enum eParticleHistograms
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!ph.fParticleHistograms[t][rs][ba]) {
          continue;
        }
        histTitle = ph.fParticleHistograms[t][rs][ba]->GetTitle();
        if (histTitle.Contains("__RUN_NUMBER__")) {
          histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
          ph.fParticleHistograms[t][rs][ba]->SetTitle(histTitle.Data());
        }
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for(int t=0;t<eParticleHistograms_N;t++) // type, see enum  eParticleHistograms

  // *) particle histograms 2D:
  for (int t = 0; t < eParticleHistograms2D_N; t++) // type, see enum eParticleHistograms2D
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!ph.fParticleHistograms2D[t][rs][ba]) {
          continue;
        }
        histTitle = ph.fParticleHistograms2D[t][rs][ba]->GetTitle();
        if (histTitle.Contains("__RUN_NUMBER__")) {
          histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
          ph.fParticleHistograms2D[t][rs][ba]->SetTitle(histTitle.Data());
        }
      } // for(int ba=0;ba<2;ba++)
    } // for(int rs=0;rs<2;rs++) // reco/sim
  } // for(int t=0;t<eParticleHistograms_N;t++) // type, see enum eParticleHistograms2D

  // *) eta separations:
  for (int ab = 0; ab < 2; ab++) {                           // ab = 0 <=> -eta , ab = 1 <=> + eta
    for (int rs = 0; rs < 2; rs++) {                         // reco/sim
      for (int ba = 0; ba < 2; ba++) {                       // before/after cuts
        for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
          if (!qv.fMabDist[ab][rs][ba][e]) {
            continue;
          }
          histTitle = qv.fMabDist[ab][rs][ba][e]->GetTitle();
          if (histTitle.Contains("__RUN_NUMBER__")) {
            histTitle.ReplaceAll("__RUN_NUMBER__", tc.fRunNumber.Data()); // it replaces in-place
            qv.fMabDist[ab][rs][ba][e]->SetTitle(histTitle.Data());
          }
        }
      }
    }
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // PropagateRunNumber()

//============================================================

template <eRecSim rs, typename T1, typename T2>
void CheckCurrentRunNumber(T1 const& collision, T2 const&)
{
  // Insanity check for the current run number and related thingies.
  // Used only during validation.

  // a) Support for Run 3 and Run 2 real data;
  // b) The rest. TBI 20241126 differentiate this support as well, e.g. for eRecSim and eSim. But Run 2 and Run 1 most likely will stay as before

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Support for Run 3 and Run 2 real data:
  //    TBI 20250112 enable other cases, after validating them
  //    TBI 20250112 Remember that I can get total run duration in converted data, but not current run duration.
  if constexpr (rs == eRec || rs == eRec_Run2) {

    // **) Check run number:
    auto bc = collision.template foundBC_as<T2>(); // I have the same code snippet at other places, keep in sync.
    if (!tc.fRunNumber.EqualTo(Form("%d", bc.runNumber()))) {
      LOGF(error, "\033[1;33m%s Run number changed within process(). This most likely indicates that a given masterjob is processing 2 or more different runs in one go.\033[0m", __FUNCTION__);
      LOGF(fatal, "tc.fRunNumber = %s, bc.runNumber() = %d", tc.fRunNumber.Data(), bc.runNumber());
    }

    // **) Check SoR, EoR, and run duration:
    auto runDuration = ccdb->getRunDuration(bc.runNumber());    // this is total run duration, not the current one (see below)
    int64_t startOfRun = std::floor(runDuration.first * 0.001); // in seconds since Unix epoch
    int64_t endOfRun = std::ceil(runDuration.second * 0.001);   // in seconds since Unix epoch
    int64_t durationInSec = endOfRun - startOfRun;              // yes, this is now run duration in seconds

    // **) Insanity check on SoR:
    if (!(tc.fRunTime[eStartOfRun] == startOfRun)) {
      LOGF(error, "\033[1;33m%s tc.fRunTime[eStartOfRun] changed within process(). This most likely indicates that a given masterjob is processing 2 or more different runs in one go.\033[0m", __FUNCTION__);
      LOGF(fatal, "tc.fRunTime[eStartOfRun] = %d, startOfRun = %d", tc.fRunTime[eStartOfRun], startOfRun);
    }

    // **) Insanity check on EoR:
    if (!(tc.fRunTime[eEndOfRun] == endOfRun)) {
      LOGF(error, "\033[1;33m%s tc.fRunTime[eEndOfRun] changed within process(). This most likely indicates that a given masterjob is processing 2 or more different runs in one go.\033[0m", __FUNCTION__);
      LOGF(fatal, "tc.fRunTime[eEndOfRun] = %d, endOfRun = %d", tc.fRunTime[eEndOfRun], endOfRun);
    }

    // **) Insanity check on run duration:
    if (!(tc.fRunTime[eDurationInSec] == durationInSec)) {
      LOGF(error, "\033[1;33m%s tc.fRunTime[eDurationInSec] changed within process(). This most likely indicates that a given masterjob is processing 2 or more different runs in one go.\033[0m", __FUNCTION__);
      LOGF(fatal, "tc.fRunTime[eDurationInSec] = %d, durationInSec = %d", tc.fRunTime[eDurationInSec], durationInSec);
    }

  } else {
    // b) The rest:

    if (!tc.fRunNumber.EqualTo(Form("%d", collision.bc().runNumber()))) {
      LOGF(error, "\033[1;33m%s Run number changed within process(). This most likely indicates that a given masterjob is processing 2 or more different runs in one go.\033[0m", __FUNCTION__);
      LOGF(fatal, "tc.fRunNumber = %s, collision.bc().runNumber() = %d", tc.fRunNumber.Data(), collision.bc().runNumber());
    }

  } // to else

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T1, typename T2> void CheckCurrentRunNumber(T1 const& collision, T2 const&)

//============================================================

void ResetEventByEventQuantities()
{
  // Reset all global event-by-event quantities here:

  // a) Event-by-event quantities;
  // b) Q-vectors;
  // c) Reset ebe containers for nested loops;
  // d) Fisher-Yates algorithm;
  // e) QA.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Event-by-event quantities:
  ebye.fSelectedTracks = 0;
  ebye.fMultiplicity = 0.;
  ebye.fReferenceMultiplicity = 0.;
  ebye.fCentrality = 0.;
  ebye.fOccupancy = 0.;
  ebye.fInteractionRate = 0.;
  ebye.fCurrentRunDuration = 0.;
  ebye.fVz = 0.;

  // b) Q-vectors:
  if (qv.fCalculateQvectors) {

    // b0) generic Q-vector:
    ResetQ();
    // b1) integrated Q-vector:
    for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
      for (int wp = 0; wp < gMaxCorrelator + 1; wp++) // weight power
      {
        qv.fQvector[h][wp] = TComplex(0., 0.);
      }
    }

    // b2) diff. Q-vector:
    for (int bin = 1; bin <= gMaxNoBinsKine; bin++) {
      qv.fqVectorEntries[PTq][bin - 1] = 0; // TBI 20240214 shall I loop also over enum's PTq and ETAq? If yes, fix it also below for qv.fqvector[PTq][bin - 1][...
      qv.fqVectorEntries[ETAq][bin - 1] = 0;
      for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
        for (int wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
          qv.fqvector[PTq][bin - 1][h][wp] = TComplex(0., 0.);
          qv.fqvector[ETAq][bin - 1][h][wp] = TComplex(0., 0.);
        } // for (int wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
      } // for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
    } // for (int b = 0; b < gMaxNoBinsKine; b++ ) {
  } // if(qv.fCalculateQvectors)

  // b3) integrated Q-vector needed for calculations with eta separations:
  if (es.fCalculateEtaSeparations) {
    for (int ab = 0; ab < 2; ab++) { // ab = 0 <=> -eta , ab = 1 <=> + eta
      for (int h = 0; h < gMaxHarmonic; h++) {
        if (es.fEtaSeparationsSkipHarmonics[h]) {
          continue;
        }
        for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
          qv.fQabVector[ab][h][e] = TComplex(0., 0.);
        }
      }
    }
    for (int ab = 0; ab < 2; ab++) {                       // ab = 0 <=> -eta , ab = 1 <=> + eta
      for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
        qv.fMab[ab][e] = 0.;
      }
    }
  }

  // b4) diff. q-vector in pt needed for calculations with eta separations:
  if (es.fCalculateEtaSeparationsAsFunctionOf[AFO_PT]) { // yes, for the time being, only as a function of pt makes sense if eta separation is used
    for (int ab = 0; ab < 2; ab++) {                     // ab = 0 <=> -eta , ab = 1 <=> + eta
      for (int bin = 1; bin <= gMaxNoBinsKine; bin++) {
        for (int h = 0; h < gMaxHarmonic; h++) {
          if (es.fEtaSeparationsSkipHarmonics[h]) {
            continue;
          }
          for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
            qv.fqabVector[ab][bin - 1][h][e] = TComplex(0., 0.); // yes, bin - 1 here
          }
        }
      }
    }

    for (int ab = 0; ab < 2; ab++) { // ab = 0 <=> -eta , ab = 1 <=> + eta
      for (int bin = 1; bin <= gMaxNoBinsKine; bin++) {
        for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
          qv.fmab[ab][bin - 1][e] = 0.; // yes, bin - 1 here
        }
      }
    }
  }

  // c) Reset ebe containers for nested loops:
  if (nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoops) {
    if (nl.ftaNestedLoops[0]) {
      nl.ftaNestedLoops[0]->Reset();
    }
    if (nl.ftaNestedLoops[1]) {
      nl.ftaNestedLoops[1]->Reset();
    }

  } // if(nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoops)

  if (nl.fCalculateKineCustomNestedLoops) {
    for (int b = 0; b < res.fResultsPro[AFO_PT]->GetNbinsX(); b++) {
      nl.ftaNestedLoopsKine[PTq][b][0]->Reset();
      nl.ftaNestedLoopsKine[PTq][b][1]->Reset();
    }
    for (int b = 0; b < res.fResultsPro[AFO_ETA]->GetNbinsX(); b++) {
      nl.ftaNestedLoopsKine[ETAq][b][0]->Reset();
      nl.ftaNestedLoopsKine[ETAq][b][1]->Reset();
    }
  } // if(nl.fCalculateKineCustomNestedLoops) {

  // d) Fisher-Yates algorithm:
  if (tc.fUseFisherYates) {
    delete tc.fRandomIndices;
    tc.fRandomIndices = NULL;
  }

  // e) QA:
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {
    for (int ba = 0; ba < 2; ba++) // before/after cuts
    {
      if (qa.fQAParticleEventProEbyE[rs][ba]) {
        qa.fQAParticleEventProEbyE[rs][ba]->Reset();
      }
    }
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void ResetEventByEventQuantities()

//============================================================

template <eRecSim rs, typename T1, typename T2>
void EventCutsCounters(T1 const& collision, T2 const& tracks)
{
  // Use this function to fill absolute and sequential event cut counters. Use only during QA, as this is computationally heavy.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Establish ordering of binning in event cut counters histograms, which resembles ordering of event cuts implementation:
  if (!ec.fEventCutCounterBinLabelingIsDone) {
    ec.fEventCutCounterBinNumber[eRec] = 1; // remember that I cannot use 'rs' here as an index, because the enum eRecSim covers separately Run 1 and Run 2, etc.
    ec.fEventCutCounterBinNumber[eSim] = 1;
    EventCuts<rs>(collision, tracks, eCutCounterBinning); // dry call, to establish the map fEventCutCounterMap and its inverse

    // **) Special treatment for event cuts implemented outside of EventCuts(), like eMultiplicity:
    //     Algorithm: I simply add eMultiplicity at the end of what was esatablished by now in the above call EventCuts<rs>(collision, tracks, eCutCounterBinning)
    //     unless proven it shall be done some other way.
    if (ec.fEventCutCounterMap[eRec]) { // TBI 20240414 also here have to hardcode 'eRec', because 'rs' spans over all enums in eRecSim => I definitely need 'generic Rec' case, perhaps via TExMap ?
                                        //              But I have already tc.fProcess[eGenericRec] and tc.fProcess[eGenericRecSim], available, shall I simply re-use them?
      ec.fEventCutCounterMap[eRec]->Add(ec.fEventCutCounterBinNumber[eRec], eMultiplicity);
      ec.fEventCutCounterMapInverse[eRec]->Add(eMultiplicity, ec.fEventCutCounterBinNumber[eRec]);
      ec.fEventCutCounterBinNumber[eRec]++; // yes
    }
    if (ec.fEventCutCounterMap[eSim]) { // TBI 20240414 also here have to hardcode 'eSim', because 'rs' spans over all enums in eRecSim => I definitely need 'generic Rec' case, perhaps via TExMap ?
                                        //              But I have already tc.fProcess[eGenericRec] and tc.fProcess[eGenericRecSim], available, shall I simply re-use them?
      ec.fEventCutCounterMap[eSim]->Add(ec.fEventCutCounterBinNumber[eSim], eMultiplicity);
      ec.fEventCutCounterMapInverse[eSim]->Add(eMultiplicity, ec.fEventCutCounterBinNumber[eSim]);
      ec.fEventCutCounterBinNumber[eSim]++; // yes
    }

    // **) Map this ordering into bin labels of actual histograms for event cut counters:
    for (int rec_sim = 0; rec_sim < 2; rec_sim++) // reco/sim => I use here exceptionally different var 'rec_sim', not the shadow 'rs' in the template parameter
    {
      for (int cc = 0; cc < eCutCounter_N; cc++) // enum eCutCounter
      {
        if (!ec.fEventCutCounterHist[rec_sim][cc]) {
          continue;
        }
        for (int bin = 1; bin < ec.fEventCutCounterBinNumber[rec_sim]; bin++) // implemented and used cuts in this analysis
        {
          ec.fEventCutCounterHist[rec_sim][cc]->GetXaxis()->SetBinLabel(bin, FancyFormatting(ec.fEventCutName[ec.fEventCutCounterMap[rec_sim]->GetValue(bin)].Data()));
        }
        for (int bin = ec.fEventCutCounterBinNumber[rec_sim]; bin <= eEventCuts_N; bin++) // implemented, but unused cuts in this analysis
        {
          ec.fEventCutCounterHist[rec_sim][cc]->GetXaxis()->SetBinLabel(bin, Form("binNo = %d (unused cut)", bin));
          // Remark: I have to write here something concrete as a bin label, if I leave "TBI" for all bin labels here for cuts which were not used,
          // I get this harmless but annoying warning during merging:
          //   Warning in <TH1Merger::CheckForDuplicateLabels>: Histogram fEventCutCounterHist[rec][seq] has duplicate labels in the x axis. Bin contents will be merged in a single bin
          // TBI 20241130 as a better solution, I shall re-define this histogram with the narower range on x-axis...
        }
        // All cuts which were implemeted, but not used I simply do not show (i can always UnZoom x-axis in TBrowser, if I want to see 'em):
        ec.fEventCutCounterHist[rec_sim][cc]->GetXaxis()->SetRangeUser(ec.fEventCutCounterHist[rec_sim][cc]->GetBinLowEdge(1), ec.fEventCutCounterHist[rec_sim][cc]->GetBinLowEdge(ec.fEventCutCounterBinNumber[rec_sim]));
      }
    }

    ec.fEventCutCounterBinLabelingIsDone = true; // this flag ensures that this specific binning is performed only once, for the first processed event
    // delete ec.fEventCutCounterMap[eRec]; // TBI 20240508 if i do not need them later, I could delete here
    // delete ec.fEventCutCounterMap[eSim];
    // delete ec.fEventCutCounterMapInverse[eRec];
    // delete ec.fEventCutCounterMapInverse[eSim];
  } // if (!ec.fEventCutCounterBinLabelingIsDone) {

  // *) Event cut counter (absolute):
  if (ec.fUseEventCutCounterAbsolute) {
    ec.fEventCutCounterBinNumber[eRec] = 1;
    ec.fEventCutCounterBinNumber[eSim] = 1;
    EventCuts<rs>(collision, tracks, eCutCounterAbsolute);

    // **) Special treatments:
    //     a) eMultiplicity: It doesn't make sense to treat this one in eCutCounterAbsolute
  }

  // *) Event cut counter (sequential):
  if (ec.fUseEventCutCounterSequential) {
    ec.fEventCutCounterBinNumber[eRec] = 1;
    ec.fEventCutCounterBinNumber[eSim] = 1;
    EventCuts<rs>(collision, tracks, eCutCounterSequential);

    // **) Special treatments:
    //     a) eMultiplicity: Since cut on eMultiplicity is implenented outside of EventCuts
    //                         I call EventCut(rs, eMultiplicity, eCutCounterSequential) directly where its implemented.
    //                         Add same treatment for other special cases, but do not forget above to expand **) Special treatment for event cuts ...
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T1, typename T2> void EventCutsCounters(T1 const& collision, T2 const& tracks, eCutModus cutModus)

//============================================================

template <eRecSim rs, typename T1, typename T2>
bool EventCuts(T1 const& collision, T2 const& tracks, eCutModus cutModus)
{
  // Event cuts on reconstructed and simulated data. Supports event cut counters, both absolute and sequential.
  // There is also a related enum eEventCuts.
  // Remark: I have added to all if statemets below which deals with floats, e.g. TMath::Abs(ebye.fCentrality - ec.fdEventCuts[eCentrality][eMax]) < tc.fFloatingPointPrecision ,
  //           to enforce the ROOT convention: "lower boundary included, upper boundary excluded"

  // a) Event cuts on reconstructed, and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1);
  // b) Event cuts only on simulated (common to Run 3, Run 2 and Run 1);
  // c) Event cuts on reconstructed, and corresponding MC truth simulated (Run 3 specific);
  // d) Event cuts on simulated (Run 3 specific);
  // e) Event cuts on reconstructed, and corresponding MC truth simulated (Run 1 and 2 specific); // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  // f) Event cuts on simulated (Run 1 and 2 specific); // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  // *) Event cuts on Test case.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Event cuts on reconstructed, and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1) ...
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

    //   *) NumberOfEvents: => this event cut is implemented directly in Steer(...)

    //   *) SelectedEvents: => this event cut is implemented directly in Steer(...)

    //   *) Offline trigger:
    //      Remark from documentation: Bypass this check if you analyse MC or continuous Run3 data.
    //      Documentation:
    //        a) O2Physics/Common/CCDB/TriggerAliases.h => available trigger aliases
    //        b) O2Physics/Common/CCDB/macros/upload_trigger_aliases.C => definitions of each trigger alias
    //      In addition: remember that I can use it only for process cases where I have joined aod::Collisions with aod::EvSels
    //      TBI 20240517 I didn't validate this trigger on Run 1, in fact, I have added protection against its usage in InsanityChecks.
    if (ec.fUseEventCuts[eTrigger]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eTrigger, eCutCounterBinning);
      } else if (ec.fsEventCuts[eTrigger].EqualTo("kINT7") && !collision.alias_bit(kINT7)) { // Validated only for Run 2
        if (!EventCut(eRec, eTrigger, cutModus)) {
          return false;
        }
      } else if (ec.fsEventCuts[eTrigger].EqualTo("kTVXinTRD") && !collision.alias_bit(kTVXinTRD)) { // Validated only for Run 3
        if (!EventCut(eRec, eTrigger, cutModus)) {
          return false;
        }
      }
      // ...
    }

    //  collision.alias_bit(kTVXinTRD);

    //   *) Sel8: // see definition in Common/TableProducer/eventSelection.cxx
    if (ec.fUseEventCuts[eSel8]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eSel8, eCutCounterBinning);
      } else if (!collision.sel8()) {
        if (!EventCut(eRec, eSel8, cutModus)) {
          return false;
        }
      }
    }

    //   *) TotalMultiplicity:
    if (ec.fUseEventCuts[eTotalMultiplicity]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eTotalMultiplicity, eCutCounterBinning);
      } else if (tracks.size() < ec.fdEventCuts[eTotalMultiplicity][eMin] || tracks.size() > ec.fdEventCuts[eTotalMultiplicity][eMax] || TMath::Abs(tracks.size() - ec.fdEventCuts[eTotalMultiplicity][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eTotalMultiplicity, cutModus)) {
          return false;
        }
      }
    }

    //   *) Multiplicity:
    //      Remark:  This cut is implemented directly in Steer(...), because I allow the possibility that ebye.fMultiplicity = ebye.fSelectedTracks .
    //               In fact, that will be true in most cases of practical interest.

    //   *) Reference multiplicity:
    //      Remark: In this member function, reference multiplicity is just a number, and any specific setting for Run 3, 2, or 1 is already done in DetermineReferenceMultiplicity(...)
    if (ec.fUseEventCuts[eReferenceMultiplicity]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eReferenceMultiplicity, eCutCounterBinning);
      } else if (ebye.fReferenceMultiplicity < ec.fdEventCuts[eReferenceMultiplicity][eMin] || ebye.fReferenceMultiplicity > ec.fdEventCuts[eReferenceMultiplicity][eMax] || TMath::Abs(ebye.fReferenceMultiplicity - ec.fdEventCuts[eReferenceMultiplicity][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eReferenceMultiplicity, cutModus)) {
          return false;
        }
      }
    }

    //   *) Centrality:
    //      Remark: In this member function, centrality is just a number, and any specific setting for Run 3, 2, or 1 is already done in DetermineCentrality(...)
    if (ec.fUseEventCuts[eCentrality]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eCentrality, eCutCounterBinning);
      } else if (ebye.fCentrality < ec.fdEventCuts[eCentrality][eMin] || ebye.fCentrality > ec.fdEventCuts[eCentrality][eMax] || TMath::Abs(ebye.fCentrality - ec.fdEventCuts[eCentrality][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eCentrality, cutModus)) {
          return false;
        }
      }
    }

    //   *) Vertex_x:
    if (ec.fUseEventCuts[eVertex_x]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eVertex_x, eCutCounterBinning);
      } else if (collision.posX() < ec.fdEventCuts[eVertex_x][eMin] || collision.posX() > ec.fdEventCuts[eVertex_x][eMax] || TMath::Abs(collision.posX() - ec.fdEventCuts[eVertex_x][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eVertex_x, cutModus)) {
          return false;
        }
      }
    }

    //   *) Vertex_y:
    if (ec.fUseEventCuts[eVertex_y]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eVertex_y, eCutCounterBinning);
      } else if (collision.posY() < ec.fdEventCuts[eVertex_y][eMin] || collision.posY() > ec.fdEventCuts[eVertex_y][eMax] || TMath::Abs(collision.posY() - ec.fdEventCuts[eVertex_y][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eVertex_y, cutModus)) {
          return false;
        }
      }
    }

    //   *) Vertex_z:
    if (ec.fUseEventCuts[eVertex_z]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eVertex_z, eCutCounterBinning);
      } else if (collision.posZ() < ec.fdEventCuts[eVertex_z][eMin] || collision.posZ() > ec.fdEventCuts[eVertex_z][eMax] || TMath::Abs(collision.posZ() - ec.fdEventCuts[eVertex_z][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eVertex_z, cutModus)) {
          return false;
        }
      }
    }

    //   *) MinVertexDistanceFromIP (minimal vertex distance from nominal Interaction Point). If vertex is closer that this value, this event is rejected:
    if (ec.fUseEventCuts[eMinVertexDistanceFromIP]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eMinVertexDistanceFromIP, eCutCounterBinning);
      } else if (sqrt(pow(collision.posX(), 2.) + pow(collision.posY(), 2.) + pow(collision.posZ(), 2.)) < ec.fdEventCuts[eMinVertexDistanceFromIP][eMin]) {
        if (!EventCut(eRec, eMinVertexDistanceFromIP, cutModus)) {
          return false;
        }
      }
    }

    //   *) NContributors:
    if (ec.fUseEventCuts[eNContributors]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eNContributors, eCutCounterBinning);
      } else if (collision.numContrib() < ec.fdEventCuts[eNContributors][eMin] || collision.numContrib() > ec.fdEventCuts[eNContributors][eMax] || TMath::Abs(collision.numContrib() - ec.fdEventCuts[eNContributors][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eNContributors, cutModus)) {
          return false;
        }
      }
    }

    // ...

    // ... and corresponding MC truth simulated:
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
    // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
    if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip..."); // TBI 20231106 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this collision
        return false;
      }

      // In this branch I can cut additionally and directly on corresponding MC truth simulated, e.g. on collision.mcCollision().posZ().
      // In case I implement something here, remember to switch from eRec to eSim when calling e.g. EventCut(...)

      // ...

    } // if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {

  } // if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // -------------------------------------------------------------------------

  // b) Event cuts only on simulated (common to Run 3, Run 2 and Run 1):
  //    Remark #1: This branch is relevant when processing ONLY simulated data at generator level.
  //    Remark #2: In this branch 'collision' is always o2::aod::McCollision, see https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

    //   *) NumberOfEvents: => this event cut is implemented directly in Steer(...)

    //   *) Impact parameter:
    if (ec.fUseEventCuts[eImpactParameter]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eSim, eImpactParameter, eCutCounterBinning);
      } else if (collision.impactParameter() < ec.fdEventCuts[eImpactParameter][eMin] || collision.impactParameter() > ec.fdEventCuts[eImpactParameter][eMax] || TMath::Abs(collision.impactParameter() - ec.fdEventCuts[eImpactParameter][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eSim, eImpactParameter, cutModus)) {
          return false;
        }
      }
    }

    //   *) Event plane angle:
    if (ec.fUseEventCuts[eEventPlaneAngle]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eSim, eEventPlaneAngle, eCutCounterBinning);
      } else if (collision.eventPlaneAngle() < ec.fdEventCuts[eEventPlaneAngle][eMin] || collision.eventPlaneAngle() > ec.fdEventCuts[eEventPlaneAngle][eMax] || TMath::Abs(collision.eventPlaneAngle() - ec.fdEventCuts[eEventPlaneAngle][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eSim, eEventPlaneAngle, cutModus)) {
          return false;
        }
      }
    }

    //   *) TotalMultiplicity:
    //      TBI 20240509 check what is the Monte Carlo analogy for tracks.size()

    //   *) Multiplicity:
    //      Remark: This cut is implemented directly in Steer(...) TBI 20240508 check how to implement this one with the current re-write

    //   *) Centrality: this is related to eImpactParameter. TBI 20240509 How do I proceed here? Shall i calculate it in void DetermineCentrality( ... ), from IP, and store it in ebye.fCentrality?

    //   *) Vertex_x:
    if (ec.fUseEventCuts[eVertex_x]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eSim, eVertex_x, eCutCounterBinning);
      } else if (collision.posX() < ec.fdEventCuts[eVertex_x][eMin] || collision.posX() > ec.fdEventCuts[eVertex_x][eMax] || TMath::Abs(collision.posX() - ec.fdEventCuts[eVertex_x][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eSim, eVertex_x, cutModus)) {
          return false;
        }
      }
    }

    //   *) Vertex_y:
    if (ec.fUseEventCuts[eVertex_y]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eSim, eVertex_y, eCutCounterBinning);
      } else if (collision.posY() < ec.fdEventCuts[eVertex_y][eMin] || collision.posY() > ec.fdEventCuts[eVertex_y][eMax] || TMath::Abs(collision.posY() - ec.fdEventCuts[eVertex_y][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eSim, eVertex_y, cutModus)) {
          return false;
        }
      }
    }

    //   *) Vertex_z:
    if (ec.fUseEventCuts[eVertex_z]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eSim, eVertex_z, eCutCounterBinning);
      } else if (collision.posZ() < ec.fdEventCuts[eVertex_z][eMin] || collision.posZ() > ec.fdEventCuts[eVertex_z][eMax] || TMath::Abs(collision.posZ() - ec.fdEventCuts[eVertex_z][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eSim, eVertex_z, cutModus)) {
          return false;
        }
      }
    }

    //   *) MinVertexDistanceFromIP (minimal vertex distance from nominal Interaction Point). If vertex is closer that this value, this event is rejected:
    if (ec.fUseEventCuts[eMinVertexDistanceFromIP]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eMinVertexDistanceFromIP, eCutCounterBinning);
      } else if (sqrt(pow(collision.posX(), 2.) + pow(collision.posY(), 2.) + pow(collision.posZ(), 2.)) < ec.fdEventCuts[eMinVertexDistanceFromIP][eMin]) {
        if (!EventCut(eRec, eMinVertexDistanceFromIP, cutModus)) {
          return false;
        }
      }
    }

    //   *) Sel8: TBI 20240509

    //   *) SelectedEvents: => this event cut is implemented directly in Steer(...)

    // ...

  } // if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

  // -------------------------------------------------------------------------

  // c) Event cuts on reconstructed, and corresponding MC truth simulated (Run 3 specific):
  //    Remark: I implement here only the event cuts which are not already in group a) above, and which make sense only for Run 3 data.
  if constexpr (rs == eRec || rs == eRecAndSim) {

    // For Run 3 multiplicities, I subscribe to o2::aod::Mults
    // See how it is defined as Joined table at https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html#o2-analysis-multiplicity-table
    // Therefore, I need always a header Common/DataModel/Multiplicity.h and o2-analysis-multiplicity-table in the workflow
    // TBI 20240509 check also o2::aod::MultExtra

    //   *) Occupancy:
    if (ec.fUseEventCuts[eOccupancy]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eOccupancy, eCutCounterBinning);
      } else if (ebye.fOccupancy < ec.fdEventCuts[eOccupancy][eMin] || ebye.fOccupancy > ec.fdEventCuts[eOccupancy][eMax] || TMath::Abs(ebye.fOccupancy - ec.fdEventCuts[eOccupancy][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eOccupancy, cutModus)) {
          return false;
        }
      }
    }

    //   *) InteractionRate:
    if (ec.fUseEventCuts[eInteractionRate]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eInteractionRate, eCutCounterBinning);
      } else if (ebye.fInteractionRate < ec.fdEventCuts[eInteractionRate][eMin] || ebye.fInteractionRate > ec.fdEventCuts[eInteractionRate][eMax] || TMath::Abs(ebye.fInteractionRate - ec.fdEventCuts[eInteractionRate][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eInteractionRate, cutModus)) {
          return false;
        }
      }
    }

    //   *) CurrentRunDuration: // TBI 20241128 check if I can use this one also on Run 2 and Run 1, most likely not
    if (ec.fUseEventCuts[eCurrentRunDuration]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eCurrentRunDuration, eCutCounterBinning);
      } else if (ebye.fCurrentRunDuration < ec.fdEventCuts[eCurrentRunDuration][eMin] || ebye.fCurrentRunDuration > ec.fdEventCuts[eCurrentRunDuration][eMax] || TMath::Abs(ebye.fCurrentRunDuration - ec.fdEventCuts[eCurrentRunDuration][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eCurrentRunDuration, cutModus)) {
          return false;
        }
      }
    }

    //   *) NoSameBunchPileup: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eNoSameBunchPileup]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eNoSameBunchPileup, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        if (!EventCut(eRec, eNoSameBunchPileup, cutModus)) {
          return false;
        }
      }
    }

    //   *) IsGoodZvtxFT0vsPV: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eIsGoodZvtxFT0vsPV]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eIsGoodZvtxFT0vsPV, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        if (!EventCut(eRec, eIsGoodZvtxFT0vsPV, cutModus)) {
          return false;
        }
      }
    }

    //   *) IsVertexITSTPC: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eIsVertexITSTPC]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eIsVertexITSTPC, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        if (!EventCut(eRec, eIsVertexITSTPC, cutModus)) {
          return false;
        }
      }
    }

    //   *) IsVertexTOFmatched: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eIsVertexTOFmatched]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eIsVertexTOFmatched, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
        if (!EventCut(eRec, eIsVertexTOFmatched, cutModus)) {
          return false;
        }
      }
    }

    //   *) IsVertexTRDmatched: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eIsVertexTRDmatched]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eIsVertexTRDmatched, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
        if (!EventCut(eRec, eIsVertexTRDmatched, cutModus)) {
          return false;
        }
      }
    }

    //   *) NoCollInTimeRangeStrict: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eNoCollInTimeRangeStrict]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eNoCollInTimeRangeStrict, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
        if (!EventCut(eRec, eNoCollInTimeRangeStrict, cutModus)) {
          return false;
        }
      }
    }

    //   *) NoCollInTimeRangeStandard: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eNoCollInTimeRangeStandard]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eNoCollInTimeRangeStandard, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        if (!EventCut(eRec, eNoCollInTimeRangeStandard, cutModus)) {
          return false;
        }
      }
    }

    //   *) NoCollInRofStrict: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eNoCollInRofStrict]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eNoCollInRofStrict, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
        if (!EventCut(eRec, eNoCollInRofStrict, cutModus)) {
          return false;
        }
      }
    }

    //   *) NoCollInRofStandard: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eNoCollInRofStandard]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eNoCollInRofStandard, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        if (!EventCut(eRec, eNoCollInRofStandard, cutModus)) {
          return false;
        }
      }
    }

    //   *) NoHighMultCollInPrevRof: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eNoHighMultCollInPrevRof]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eNoHighMultCollInPrevRof, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
        if (!EventCut(eRec, eNoHighMultCollInPrevRof, cutModus)) {
          return false;
        }
      }
    }

    //   *) IsGoodITSLayer3: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eIsGoodITSLayer3]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eIsGoodITSLayer3, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)) {
        if (!EventCut(eRec, eIsGoodITSLayer3, cutModus)) {
          return false;
        }
      }
    }

    //   *) IsGoodITSLayer0123: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eIsGoodITSLayer0123]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eIsGoodITSLayer0123, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
        if (!EventCut(eRec, eIsGoodITSLayer0123, cutModus)) {
          return false;
        }
      }
    }

    //   *) IsGoodITSLayersAll: // see O2Physics/Common/CCDB/EventSelectionParams.cxx
    if (ec.fUseEventCuts[eIsGoodITSLayersAll]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eIsGoodITSLayersAll, eCutCounterBinning);
      } else if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        if (!EventCut(eRec, eIsGoodITSLayersAll, cutModus)) {
          return false;
        }
      }
    }

    // ...

    //  *) Centrality weights (flattening):
    // Remark 1: Since I am getting centrality weights from centrality distribution AFTER all the events cuts, flattening must be applied here after all other event cuts:
    // Remark 2: Whatever I change here, change also in the corresponding branch for Run 2 and Run 1.
    //           Yes, I have to replicate for this special event cut the same code, since in each case it has to be applied at the very end.
    if (ec.fUseEventCuts[eCentralityWeights]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eCentralityWeights, eCutCounterBinning);
      } else if (gRandom->Uniform(0, 1) > CentralityWeight(ebye.fCentrality)) { // yes, since centralityWeight is normalized probability (see CentralityWeight(...))
        if (!EventCut(eRec, eCentralityWeights, cutModus)) {
          return false;
        }
      }
    }

    // Remark: If I need any further event cut, implement it BEFORE event cut "Centrality weights (flattening)", which must be implemented last.

    // ... and corresponding MC truth simulated (Run 3 specific):
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
    // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
    if constexpr (rs == eRecAndSim) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip..."); // TBI 20231106 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this collision
        return false;
      }

      // In this branch I can cut additionally and directly on corresponding MC truth simulated.
      // Remark: I implement here only the event cuts which are not already in group a) above, and which make sense only for Run 3 data.
      // In case I implement something here, remember to switch from eRec to eSim when calling e.g. EventCut(...)

      // ...

    } // if constexpr (rs == eRecAndSim) {

  } // if constexpr (rs == eRec || rs == eRecAndSim) {

  // -------------------------------------------------------------------------

  // d) Event cuts on simulated (Run 3 specific):
  //    Remark #1: I implement here only the event cuts which are not already in group b) above, and which make sense only for Run 3 data.
  //    Remark #2: In this branch 'collision' is always o2::aod::McCollision, see https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
  //               See how I handled the case b) above.
  if constexpr (rs == eSim) {

    // ...

  } // if constexpr (rs == eSim) {

  // -------------------------------------------------------------------------

  // e) Event cuts on reconstructed, and corresponding MC truth simulated (Run 1 and 2 specific):
  //    Remark: I implement here only the event cuts which are not already in group a) above, and which make sense only for Run 1 and 2 data.
  if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

    //   *) Sel7:
    if (ec.fUseEventCuts[eSel7]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eSel7, eCutCounterBinning);
      } else if (!collision.sel7()) {
        if (!EventCut(eRec, eSel7, cutModus)) {
          return false;
        }
      }
    }

    // ...

    //  *) Centrality weights (flattening):
    // Remark 1: Since I am getting centrality weights from centrality distribution AFTER all the events cuts, flattening must be applied here after all other event cuts:
    // Remark 2: Whatever I change here, change also in the corresponding branch for Run 3.
    //           Yes, I have to replicate for this special event cut the same code, since in each case it has to be applied at the very end.
    if (ec.fUseEventCuts[eCentralityWeights]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eCentralityWeights, eCutCounterBinning);
      } else if (gRandom->Uniform(0, 1) > CentralityWeight(ebye.fCentrality)) { // yes, since centralityWeight is normalized probability (see CentralityWeight(...))
        if (!EventCut(eRec, eCentralityWeights, cutModus)) {
          return false;
        }
      }
    }

    // Remark: If I need any further event cut, implement it BEFORE event cut "Centrality weights (flattening)", which must be implemented last.

    // ... and corresponding MC truth simulated:
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
    // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
    if constexpr (rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip..."); // TBI 20231106 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this collision
        return false;
      }

      // In this branch I can cut additionally and directly on corresponding MC truth simulated.
      // Remark: I implement here only the event cuts which are not already in group a) above, and which make sense only for Run 1 and 2 data.
      // In case I implement something here, remember to switch from eRec to eSim when calling e.g. EventCut(...)

      // ...

    } // if constexpr (rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {

  } // if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // -------------------------------------------------------------------------

  // f) Event cuts on simulated (Run 1 and 2 specific)
  //    Remark #1: I implement here only the event cuts which are not already in group b) above, and which make sense only for Run 1 and 2 data.
  //    Remark #2: In this branch 'collision' is always o2::aod::McCollision, see https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
  //               See how I handled the case b) above.
  if constexpr (rs == eSim_Run2 || rs == eSim_Run1) {

    // ...

  } // if constexpr (rs == eSim_Run2 || rs == eSim_Run1) {

  // -------------------------------------------------------------------------

  // *) Test case:
  if constexpr (rs == eTest) {
    // This branch corresponds to process with minimal subscription - I implement just a few example cuts, just for testing purposes.
    // Only eRec is support in Test for the time being.

    //   *) TotalMultiplicity:
    if (ec.fUseEventCuts[eTotalMultiplicity]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eTotalMultiplicity, eCutCounterBinning);
      } else if (tracks.size() < ec.fdEventCuts[eTotalMultiplicity][eMin] || tracks.size() > ec.fdEventCuts[eTotalMultiplicity][eMax] || TMath::Abs(tracks.size() - ec.fdEventCuts[eTotalMultiplicity][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eTotalMultiplicity, cutModus)) {
          return false;
        }
      }
    }

    //   *) Vertex_z:
    if (ec.fUseEventCuts[eVertex_z]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eSim, eVertex_z, eCutCounterBinning);
      } else if (collision.posZ() < ec.fdEventCuts[eVertex_z][eMin] || collision.posZ() > ec.fdEventCuts[eVertex_z][eMax] || TMath::Abs(collision.posZ() - ec.fdEventCuts[eVertex_z][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eSim, eVertex_z, cutModus)) {
          return false;
        }
      }
    }

    //   *) Centrality:
    if (ec.fUseEventCuts[eCentrality]) {
      if (cutModus == eCutCounterBinning) {
        EventCut(eRec, eCentrality, eCutCounterBinning);
      } else if (ebye.fCentrality < ec.fdEventCuts[eCentrality][eMin] || ebye.fCentrality > ec.fdEventCuts[eCentrality][eMax] || TMath::Abs(ebye.fCentrality - ec.fdEventCuts[eCentrality][eMax]) < tc.fFloatingPointPrecision) {
        if (!EventCut(eRec, eCentrality, cutModus)) {
          return false;
        }
      }
    }

    // ...

  } // if constexpr (rs == eTest) {

  return true;

} // template <eRecSim rs, typename T1, typename T2> bool EventCuts(T1 const& collision, T2 const& tracks)

//============================================================

bool EventCut(int rs, int eventCut, eCutModus cutModus)
{
  // Helper function to reduce code bloat in EventCuts(). It's meant to be used only in EventCuts().
  // It can be used also in exceptional cases outside of EventCuts(), like for eMultiplicity, but use with care.
  // For instance, I can call EventCut(eRec, eCentrality, eCutCounterSequential) directly, only if I have checked that
  // fUseEventCutCounterSequential is true, etc.

  // Remark: Remember that as a second argument I cannot use enum eEventCuts, because here in one go I take both enum eEventCuts and enum eEventHistograms .

  // *) Insanity checks on arguments:
  if (!(0 == rs || 1 == rs)) {
    LOGF(fatal, "\033[1;31m%s at line %d : 'rs' must be generic Rec or Sim index, rs = %d \033[0m", __FUNCTION__, __LINE__, rs);
  }
  if (eventCut >= eEventCuts_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : eventCut >= eEventCuts_N, eventCut = %d , eEventCuts_N = %d \033[0m", __FUNCTION__, __LINE__, eventCut, static_cast<int>(eEventCuts_N));
  }

  // *) Do the thing:
  switch (cutModus) {
    case eCut:
      if (tc.fVerboseEventCut) {
        LOGF(info, "\033[1;31mEvent didn't survive the cut: %s\033[0m", ec.fEventCutName[eventCut].Data());
      }
      return false;
      break;
    case eCutCounterBinning:
      ec.fEventCutCounterMap[rs]->Add(ec.fEventCutCounterBinNumber[rs], eventCut);
      ec.fEventCutCounterMapInverse[rs]->Add(eventCut, ec.fEventCutCounterBinNumber[rs]);
      ec.fEventCutCounterBinNumber[rs]++; // yes
      return true;
      break;
    case eCutCounterAbsolute:
      ec.fEventCutCounterHist[rs][eAbsolute]->Fill(ec.fEventCutCounterMapInverse[rs]->GetValue(eventCut));
      return true; // yes, so that I can proceed with another cut in EventCuts
      break;
    case eCutCounterSequential:
      ec.fEventCutCounterHist[rs][eSequential]->Fill(ec.fEventCutCounterMapInverse[rs]->GetValue(eventCut));
      return false; // yes, so that I bail out from EventCuts
      break;
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : This cutModus = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(cutModus));
      break;
  } // switch(cutModus)

  return false; // obsolete, but it suppresses the warning...

} // bool EventCut(int rs, int eventCut, eCutModus cutModus)

//============================================================

bool RemainingEventCuts()
{
  // Remaining event cuts which can be applied ONLY after the main loop over particles.
  // For instance, cut on total number of selected particles (eMultiplicity).
  // Remark #1: Whichever cut I implement here, update EventCutsCounters(...) for that cut (like I did for eMultiplicity, as a sort of template).
  // Remark #2: I do not have here templated arguments like in EventCuts(), because I do not anticipate using any getter from the framework directly here.
  // Remark #3: With the current implementation, I support here only eCutCounterSequential, i.e. eCutCounterAbsolute is not supported for cuts applied here.

  // a) Determine if this function was called for generic rec or generic sim:
  // *) eMultiplicity;
  // ...

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Determine if this function was called for generic rec or generic sim:
  //    Remark: I can do it in this simplified way, because I do not anticipate I will call here any getters from the framework.
  int rs = -1;
  if (tc.fProcess[eGenericRec] || tc.fProcess[eGenericRecSim]) {
    rs = eRec; // yes, I do not count in RecSim mode separately particles and rec and sim level which survived particle cuts
  } else if (tc.fProcess[eGenericSim]) {
    rs = eSim;
  }

  // *) Multiplicity: (see documentation for ebye.fMultiplicity for its definition)
  if (ec.fUseEventCuts[eMultiplicity]) {
    if (ebye.fMultiplicity < ec.fdEventCuts[eMultiplicity][eMin] || ebye.fMultiplicity > ec.fdEventCuts[eMultiplicity][eMax] || TMath::Abs(ebye.fMultiplicity - ec.fdEventCuts[eMultiplicity][eMax]) < tc.fFloatingPointPrecision) {
      // Remark: I have to implement RemainingEventCuts() in a slightly different way as EventCuts()
      EventCut(rs, eMultiplicity, eCut);      // just a printout that this event didn't survive this cut
      if (ec.fUseEventCutCounterSequential) { // yes, this is important. Otherwise fEventCutCounterHist can be used in EventCut(...), even though it's NULL
        EventCut(rs, eMultiplicity, eCutCounterSequential);
      }
      return false;
    }
  }

  return true;

} // bool RemainingEventCuts()

//============================================================

template <eRecSim rs>
void FillSubeventMultiplicities()
{
  // Fill subevent (defined via eta separation) multiplicities.

  // a) Fill reconstructed (common to Run 3, Run 2 and Run 1 + Test mode);
  // b) Fill only simulated (common to Run 3, Run 2 and Run 1).

  // Remark: This function has to be called after Q-vectors are filled. It makes sense to fill these histograms only for "eAfter",
  //         becase Q-vectors are not filled before the event cuts.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Fill reconstructed (common to Run 3, Run 2 and Run 1 + Test mode):
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1 || rs == eTest) {
    for (int ab = 0; ab < 2; ab++) {                       // ab = 0 <=> -eta , ab = 1 <=> + eta
      for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
        !qv.fMabDist[ab][eRec][eAfter][e] ? true : qv.fMabDist[ab][eRec][eAfter][e]->Fill(qv.fMab[ab][e]);
      }
    }
  }

  // b) Fill only simulated (common to Run 3, Run 2 and Run 1):
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {
    for (int ab = 0; ab < 2; ab++) {                       // ab = 0 <=> -eta , ab = 1 <=> + eta
      for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
        !qv.fMabDist[ab][eSim][eAfter][e] ? true : qv.fMabDist[ab][eSim][eAfter][e]->Fill(qv.fMab[ab][e]);
      }
    }
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void FillSubeventMultiplicities()

//============================================================

template <eRecSim rs, typename T1, typename T2>
void FillEventHistograms(T1 const& collision, T2 const& tracks, eBeforeAfter ba)
{
  // Fill all event histograms for reconstructed or simulated data. QA event histograms are also filled here.

  // a) Fill reconstructed, and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1);
  // b) Fill only simulated (common to Run 3, Run 2 and Run 1);
  // c) Fill reconstructed (Run 3 specific);
  // d) Fill only simulated (Run 3 specific);
  // e) Fill reconstructed (Run 1 and 2 specific); // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  // f) Fill only simulated (Run 1 and 2 specific); // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  // g) Test case.

  // Remark: in most cases, all histogram which depend on eMultiplicity are booked only for "after", because by default, Multiplicity = SelectedTracks.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Fill reconstructed ... (common to Run 3, Run 2 and Run 1):
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    if (eh.fFillEventHistograms) {
      !eh.fEventHistograms[eNumberOfEvents][eRec][ba] ? true : eh.fEventHistograms[eNumberOfEvents][eRec][ba]->Fill(0.5); // basically, if histogram is not booked, do nothing. 'true' is a placeholder, for the time being
      !eh.fEventHistograms[eVertex_x][eRec][ba] ? true : eh.fEventHistograms[eVertex_x][eRec][ba]->Fill(collision.posX());
      !eh.fEventHistograms[eVertex_y][eRec][ba] ? true : eh.fEventHistograms[eVertex_y][eRec][ba]->Fill(collision.posY());
      !eh.fEventHistograms[eVertex_z][eRec][ba] ? true : eh.fEventHistograms[eVertex_z][eRec][ba]->Fill(collision.posZ());
      !eh.fEventHistograms[eNContributors][eRec][ba] ? true : eh.fEventHistograms[eNContributors][eRec][ba]->Fill(collision.numContrib());
      !eh.fEventHistograms[eTotalMultiplicity][eRec][ba] ? true : eh.fEventHistograms[eTotalMultiplicity][eRec][ba]->Fill(tracks.size()); // TBI 20231106 check and validate further
      !eh.fEventHistograms[eMultiplicity][eRec][ba] ? true : eh.fEventHistograms[eMultiplicity][eRec][ba]->Fill(ebye.fMultiplicity);
      !eh.fEventHistograms[eReferenceMultiplicity][eRec][ba] ? true : eh.fEventHistograms[eReferenceMultiplicity][eRec][ba]->Fill(ebye.fReferenceMultiplicity);
      !eh.fEventHistograms[eCentrality][eRec][ba] ? true : eh.fEventHistograms[eCentrality][eRec][ba]->Fill(ebye.fCentrality);
    }

    // QA:
    if (qa.fFillQAEventHistograms2D) {
      !qa.fQAEventHistograms2D[eMultiplicity_vs_ReferenceMultiplicity][eRec][ba] ? true : qa.fQAEventHistograms2D[eMultiplicity_vs_ReferenceMultiplicity][eRec][ba]->Fill(ebye.fMultiplicity, ebye.fReferenceMultiplicity);
      !qa.fQAEventHistograms2D[eMultiplicity_vs_NContributors][eRec][ba] ? true : qa.fQAEventHistograms2D[eMultiplicity_vs_NContributors][eRec][ba]->Fill(ebye.fMultiplicity, collision.numContrib());
      !qa.fQAEventHistograms2D[eMultiplicity_vs_Centrality][eRec][ba] ? true : qa.fQAEventHistograms2D[eMultiplicity_vs_Centrality][eRec][ba]->Fill(ebye.fMultiplicity, ebye.fCentrality);
      !qa.fQAEventHistograms2D[eMultiplicity_vs_Vertex_z][eRec][ba] ? true : qa.fQAEventHistograms2D[eMultiplicity_vs_Vertex_z][eRec][ba]->Fill(ebye.fMultiplicity, collision.posZ());
      !qa.fQAEventHistograms2D[eReferenceMultiplicity_vs_NContributors][eRec][ba] ? true : qa.fQAEventHistograms2D[eReferenceMultiplicity_vs_NContributors][eRec][ba]->Fill(ebye.fReferenceMultiplicity, collision.numContrib());
      !qa.fQAEventHistograms2D[eReferenceMultiplicity_vs_Centrality][eRec][ba] ? true : qa.fQAEventHistograms2D[eReferenceMultiplicity_vs_Centrality][eRec][ba]->Fill(ebye.fReferenceMultiplicity, ebye.fCentrality);
      !qa.fQAEventHistograms2D[eReferenceMultiplicity_vs_Vertex_z][eRec][ba] ? true : qa.fQAEventHistograms2D[eReferenceMultiplicity_vs_Vertex_z][eRec][ba]->Fill(ebye.fReferenceMultiplicity, collision.posZ());
      !qa.fQAEventHistograms2D[eNContributors_vs_Centrality][eRec][ba] ? true : qa.fQAEventHistograms2D[eNContributors_vs_Centrality][eRec][ba]->Fill(collision.numContrib(), ebye.fCentrality);
      !qa.fQAEventHistograms2D[eNContributors_vs_Vertex_z][eRec][ba] ? true : qa.fQAEventHistograms2D[eNContributors_vs_Vertex_z][eRec][ba]->Fill(collision.numContrib(), collision.posZ());
      !qa.fQAEventHistograms2D[eCentrality_vs_Vertex_z][eRec][ba] ? true : qa.fQAEventHistograms2D[eCentrality_vs_Vertex_z][eRec][ba]->Fill(ebye.fCentrality, collision.posZ());
    }

    if (qa.fFillQACorrelationsVsHistograms2D && qa.fQAParticleEventProEbyE[eRec][ba] && ba == eAfter) { // fill only for eAfter, because I do not calculate Q-vectors before cuts

      // Calculate quickly 2-p correlation in harmonic h for this event: TBI 20250114 shall I add this also to some EbyE variable? There is no really much of a code bloat for the time being...

      // Flush 'n' fill the generic Q-vectors:
      ResetQ();
      int lMaxCorrelator = 2; // used only here locally
      for (int h = 0; h < gMaxHarmonic * lMaxCorrelator + 1; h++) {
        for (int wp = 0; wp < lMaxCorrelator + 1; wp++) // weight power
        {
          qv.fQ[h][wp] = qv.fQvector[h][wp];
        }
      }

      for (int h = 1; h <= gMaxHarmonic; h++) {
        TComplex two = Two(h, -h);
        double twoC = two.Re(); // cos
        // double twoS = two.Im(); // sin
        double wTwo = Two(0, 0).Re(); // Weight is 'number of combinations' by default TBI
                                      // 20220809 add support for other weights
        if (wTwo > 0.0) {
          twoC /= wTwo;
        } else {
          LOGF(fatal, "In function \033[1;31m%s at line %d, wTwo = %f <=0. ebye.fSelectedTracks = %d\033[0m", __FUNCTION__, __LINE__, wTwo, ebye.fSelectedTracks);
        }
        !qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_Multiplicity][h - 1][eRec] ? true : qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_Multiplicity][h - 1][eRec]->Fill(twoC, ebye.fMultiplicity);
        !qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_ReferenceMultiplicity][h - 1][eRec] ? true : qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_ReferenceMultiplicity][h - 1][eRec]->Fill(twoC, ebye.fReferenceMultiplicity);
        !qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_Centrality][h - 1][eRec] ? true : qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_Centrality][h - 1][eRec]->Fill(twoC, ebye.fCentrality);
        // .....
        !qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_MeanPhi][h - 1][eRec] ? true : qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_MeanPhi][h - 1][eRec]->Fill(twoC, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eMeanPhi));
        !qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_MeanPt][h - 1][eRec] ? true : qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_MeanPt][h - 1][eRec]->Fill(twoC, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eMeanPt));
        !qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_MeanEta][h - 1][eRec] ? true : qa.fQACorrelationsVsHistograms2D[eCorrelations_vs_MeanEta][h - 1][eRec]->Fill(twoC, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eMeanEta));
        // .....
      }

      // Flush the generic Q-vectors:
      ResetQ();

    } // if (qa.fFillQACorrelationsVsHistograms2D && qa.fQAParticleEventProEbyE[eRec][ba] && ba == eAfter) {

    // ... and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1) ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }
      if (eh.fFillEventHistograms) {
        !eh.fEventHistograms[eNumberOfEvents][eSim][ba] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][ba]->Fill(0.5);
        !eh.fEventHistograms[eVertex_x][eSim][ba] ? true : eh.fEventHistograms[eVertex_x][eSim][ba]->Fill(collision.mcCollision().posX());
        !eh.fEventHistograms[eVertex_y][eSim][ba] ? true : eh.fEventHistograms[eVertex_y][eSim][ba]->Fill(collision.mcCollision().posY());
        !eh.fEventHistograms[eVertex_z][eSim][ba] ? true : eh.fEventHistograms[eVertex_z][eSim][ba]->Fill(collision.mcCollision().posZ());
        !eh.fEventHistograms[eImpactParameter][eSim][ba] ? true : eh.fEventHistograms[eImpactParameter][eSim][ba]->Fill(collision.mcCollision().impactParameter());
        !eh.fEventHistograms[eEventPlaneAngle][eSim][ba] ? true : eh.fEventHistograms[eEventPlaneAngle][eSim][ba]->Fill(collision.mcCollision().eventPlaneAngle());
        // eh.fEventHistograms[eTotalMultiplicity][eSim][ba]->Fill(tracks.size()); // TBI 20231106 check how to get corresponding MC truth info, and validate further
        // eh.fEventHistograms[eMultiplicity][eSim][ba]->Fill(ebye.fMultiplicity); // TBI 20241123 re-think if I really need it here. If yes, most likely I will have to
        //              generalize fSelectedTracks to an array, to counter separately selected sim particles
        // eh.fEventHistograms[eCentrality][eSim][ba]->Fill(ebye.fCentrality); // TBI 20240120 this case is still not supported in DetermineCentrality()
      }

      // QA:
      if (qa.fFillQAEventHistograms2D) {
        !qa.fQAEventHistograms2D[eCentrality_vs_ImpactParameter][eSim][ba] ? true : qa.fQAEventHistograms2D[eCentrality_vs_ImpactParameter][eSim][ba]->Fill(ebye.fCentrality, collision.mcCollision().impactParameter());
        // ...
      }
    } // if constexpr (rs == eRecAndSim) {
  } // if constexpr (rs == eRec || rs == eRecAndSim) {

  // -----------------------------------------------------------------------------

  // b) Fill only simulated (common to Run 3, Run 2 and Run 1):
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {
    if (eh.fFillEventHistograms) {
      !eh.fEventHistograms[eImpactParameter][eSim][ba] ? true : eh.fEventHistograms[eImpactParameter][eSim][ba]->Fill(collision.impactParameter()); // yes, because in this branch 'collision' is always aod::McCollision
      !eh.fEventHistograms[eEventPlaneAngle][eSim][ba] ? true : eh.fEventHistograms[eEventPlaneAngle][eSim][ba]->Fill(collision.eventPlaneAngle()); // yes, because in this branch 'collision' is always aod::McCollision
      !eh.fEventHistograms[eMultiplicity][eSim][ba] ? true : eh.fEventHistograms[eMultiplicity][eSim][ba]->Fill(ebye.fMultiplicity);
      // eh.fEventHistograms[eCentrality][eSim][ba]->Fill(ebye.fCentrality); // TBI 20240120 this case is still not supported in DetermineCentrality()
      // eh.fEventHistograms[eReferenceMultiplicity][eSim][ba]->Fill(ebye.fReferenceMultiplicity); // TBI 20241123 this case is still not supported in DetermineReferenceMultiplicity()
      // eh.fEventHistograms[eTotalMultiplicity][eSim][ba]->Fill(tracks.size()); // TBI 20231030 check further how to use the same thing for 'sim'
    }

    // Eta separations:
    if (es.fCalculateEtaSeparations) {
      for (int ab = 0; ab < 2; ab++) {                       // ab = 0 <=> -eta , ab = 1 <=> + eta
        for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
          !qv.fMabDist[ab][eSim][ba][e] ? true : qv.fMabDist[ab][eSim][ba][e]->Fill(qv.fMab[ab][e]);
        }
      }
    }
  }

  // -----------------------------------------------------------------------------

  // c) Fill reconstructed (Run 3 specific):
  if constexpr (rs == eRec || rs == eRecAndSim) {
    if (eh.fFillEventHistograms) {
      !eh.fEventHistograms[eOccupancy][eRec][ba] ? true : eh.fEventHistograms[eOccupancy][eRec][ba]->Fill(ebye.fOccupancy);
      !eh.fEventHistograms[eInteractionRate][eRec][ba] ? true : eh.fEventHistograms[eInteractionRate][eRec][ba]->Fill(ebye.fInteractionRate);
      !eh.fEventHistograms[eCurrentRunDuration][eRec][ba] ? true : eh.fEventHistograms[eCurrentRunDuration][eRec][ba]->Fill(ebye.fCurrentRunDuration);
    }
    // QA:
    if (qa.fFillQAEventHistograms2D) {
      // General (estimators can be chosen via configurables):
      !qa.fQAEventHistograms2D[eMultiplicity_vs_Occupancy][eRec][ba] ? true : qa.fQAEventHistograms2D[eMultiplicity_vs_Occupancy][eRec][ba]->Fill(ebye.fMultiplicity, ebye.fOccupancy);
      !qa.fQAEventHistograms2D[eReferenceMultiplicity_vs_Occupancy][eRec][ba] ? true : qa.fQAEventHistograms2D[eReferenceMultiplicity_vs_Occupancy][eRec][ba]->Fill(ebye.fReferenceMultiplicity, ebye.fOccupancy);
      !qa.fQAEventHistograms2D[eNContributors_vs_Occupancy][eRec][ba] ? true : qa.fQAEventHistograms2D[eNContributors_vs_Occupancy][eRec][ba]->Fill(collision.numContrib(), ebye.fOccupancy);
      !qa.fQAEventHistograms2D[eCentrality_vs_Occupancy][eRec][ba] ? true : qa.fQAEventHistograms2D[eCentrality_vs_Occupancy][eRec][ba]->Fill(ebye.fCentrality, ebye.fOccupancy);
      !qa.fQAEventHistograms2D[eVertex_z_vs_Occupancy][eRec][ba] ? true : qa.fQAEventHistograms2D[eVertex_z_vs_Occupancy][eRec][ba]->Fill(collision.posZ(), ebye.fOccupancy);
      // ...

      // Specific (estimators are hardwired):
      !qa.fQAEventHistograms2D[eMultNTracksPV_vs_MultNTracksGlobal][eRec][ba] ? true : qa.fQAEventHistograms2D[eMultNTracksPV_vs_MultNTracksGlobal][eRec][ba]->Fill(qa.fReferenceMultiplicity[eMultNTracksPV], qa.fReferenceMultiplicity[eMultNTracksGlobal]); // TBI 20241209 check if I can use this one for Run 2 and 1
      !qa.fQAEventHistograms2D[eCentFT0C_vs_CentFT0CVariant1][eRec][ba] ? true : qa.fQAEventHistograms2D[eCentFT0C_vs_CentFT0CVariant1][eRec][ba]->Fill(qa.fCentrality[eCentFT0C], qa.fCentrality[eCentFT0CVariant1]);
      !qa.fQAEventHistograms2D[eCentFT0C_vs_CentFT0M][eRec][ba] ? true : qa.fQAEventHistograms2D[eCentFT0C_vs_CentFT0M][eRec][ba]->Fill(qa.fCentrality[eCentFT0C], qa.fCentrality[eCentFT0M]);
      !qa.fQAEventHistograms2D[eCentFT0C_vs_CentFV0A][eRec][ba] ? true : qa.fQAEventHistograms2D[eCentFT0C_vs_CentFV0A][eRec][ba]->Fill(qa.fCentrality[eCentFT0C], qa.fCentrality[eCentFV0A]);
      !qa.fQAEventHistograms2D[eCentFT0C_vs_CentNTPV][eRec][ba] ? true : qa.fQAEventHistograms2D[eCentFT0C_vs_CentNTPV][eRec][ba]->Fill(qa.fCentrality[eCentFT0C], qa.fCentrality[eCentNTPV]);
      !qa.fQAEventHistograms2D[eCentFT0C_vs_CentNGlobal][eRec][ba] ? true : qa.fQAEventHistograms2D[eCentFT0C_vs_CentNGlobal][eRec][ba]->Fill(qa.fCentrality[eCentFT0C], qa.fCentrality[eCentNGlobal]);
      !qa.fQAEventHistograms2D[eCentFT0M_vs_CentNTPV][eRec][ba] ? true : qa.fQAEventHistograms2D[eCentFT0M_vs_CentNTPV][eRec][ba]->Fill(qa.fCentrality[eCentFT0M], qa.fCentrality[eCentNTPV]);
      !qa.fQAEventHistograms2D[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange][eRec][ba] ? true : qa.fQAEventHistograms2D[eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange][eRec][ba]->Fill(collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());
      !qa.fQAEventHistograms2D[eCurrentRunDuration_vs_InteractionRate][eRec][ba] ? true : qa.fQAEventHistograms2D[eCurrentRunDuration_vs_InteractionRate][eRec][ba]->Fill(ebye.fCurrentRunDuration, ebye.fInteractionRate);
    }

    if (qa.fFillQAParticleEventHistograms2D && qa.fQAParticleEventProEbyE[eRec][ba]) {
      // This is a special category, where I do correlation <some-particle-property> vs. some-event-property.
      // I use 'number of combinations' as a weight, which here reduces simply to the 'number of entries' weight.
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsEbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsEbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eitsNClsEbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(eitsNClsEbyE));
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsNegEtaEbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsNegEtaEbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eitsNClsNegEtaEbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(eitsNClsNegEtaEbyE));
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsPosEtaEbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_itsNClsPosEtaEbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eitsNClsPosEtaEbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(eitsNClsPosEtaEbyE));
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0804EbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0804EbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eEta0804EbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(eEta0804EbyE));
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0400EbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0400EbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eEta0400EbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(eEta0400EbyE));
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0004EbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0004EbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eEta0004EbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(eEta0004EbyE));
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0408EbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Eta0408EbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(eEta0408EbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(eEta0408EbyE));
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt0005EbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt0005EbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(ePt0005EbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(ePt0005EbyE));
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt0510EbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt0510EbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(ePt0510EbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(ePt0510EbyE));
      !qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt1050EbyE][eRec][ba] ? true : qa.fQAParticleEventHistograms2D[eCurrentRunDuration_vs_Pt1050EbyE][eRec][ba]->Fill(ebye.fCurrentRunDuration, qa.fQAParticleEventProEbyE[eRec][ba]->GetBinContent(ePt1050EbyE), qa.fQAParticleEventProEbyE[eRec][ba]->GetBinEntries(ePt1050EbyE));

      // ...

    } // if (qa.fFillQAParticleEventHistograms2D && qa.fQAParticleEventProEbyE[eRec][ba]) {

    // ... and corresponding MC truth simulated (Run 3 specific)
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
    if constexpr (rs == eRecAndSim) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }

      // !eh.fEventHistograms[eMultMCNParticlesEta08][eSim][ba] ? true : eh.fEventHistograms[eMultMCNParticlesEta08][eSim][ba]->Fill(collision.multMCNParticlesEta08());

      // !eh.fEventHistograms[eNumberOfEvents][eSim][ba] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][ba]->Fill(0.5);
    } // if constexpr (rs == eRecAndSim) {
  } // if constexpr (rs == eRec || rs == eRecAndSim) {

  // -----------------------------------------------------------------------------

  // d) Fill only simulated (Run 3 specific):
  if constexpr (rs == eSim) {
    // !eh.fEventHistograms[eImpactParameter][eSim][ba] ? true : eh.fEventHistograms[eImpactParameter][eSim][ba]->Fill(collision.impactParameter()); // yes, because in this branch 'collision' is always aod::McCollision
    // !eh.fEventHistograms[eEventPlaneAngle][eSim][ba] ? true : eh.fEventHistograms[eEventPlaneAngle][eSim][ba]->Fill(collision.eventPlaneAngle()); // yes, because in this branch 'collision' is always aod::McCollision
  } // if constexpr (rs == eSim) {

  // -----------------------------------------------------------------------------

  // e) Fill reconstructed (Run 1 and 2 specific): // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    if (eh.fFillEventHistograms) {
      // ...
    }
    // QA:
    if (qa.fFillQAEventHistograms2D) {
      !qa.fQAEventHistograms2D[eCentRun2V0M_vs_CentRun2SPDTracklets][eRec][ba] ? true : qa.fQAEventHistograms2D[eCentRun2V0M_vs_CentRun2SPDTracklets][eRec][ba]->Fill(qa.fCentrality[eCentRun2V0M], qa.fCentrality[eCentRun2SPDTracklets]);
    }

    // ... and corresponding MC truth simulated (Run 1 and Run 2 specific):
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
    if constexpr (rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }
      // !eh.fEventHistograms[eNumberOfEvents][eSim][ba] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][ba]->Fill(0.5);

    } // if constexpr (rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
  } // if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // -----------------------------------------------------------------------------

  // f) Fill only simulated (Run 1 and 2 specific): // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  if constexpr (rs == eSim_Run2 || rs == eSim_Run1) {
    // !eh.fEventHistograms[eImpactParameter][eSim][ba] ? true : eh.fEventHistograms[eImpactParameter][eSim][ba]->Fill(collision.impactParameter()); // yes, because in this branch 'collision' is always aod::McCollision
    // !eh.fEventHistograms[eEventPlaneAngle][eSim][ba] ? true : eh.fEventHistograms[eEventPlaneAngle][eSim][ba]->Fill(collision.eventPlaneAngle()); // yes, because in this branch 'collision' is always aod::McCollision
  } // if constexpr (rs == eSim_Run2 || rs == eSim_Run1) {

  // -----------------------------------------------------------------------------

  // g) Test case:
  if constexpr (rs == eTest) {
    // TBI 20240223 for the time being, eTest fills only eRec histos:
    // A few example histograms, just to check if I access corresponding tables:
    if (eh.fFillEventHistograms) {
      !eh.fEventHistograms[eVertex_z][eRec][ba] ? true : eh.fEventHistograms[eVertex_z][eRec][ba]->Fill(collision.posZ());
      !eh.fEventHistograms[eTotalMultiplicity][eRec][ba] ? true : eh.fEventHistograms[eTotalMultiplicity][eRec][ba]->Fill(tracks.size());
      !eh.fEventHistograms[eCentrality][eRec][ba] ? true : eh.fEventHistograms[eCentrality][eRec][ba]->Fill(ebye.fCentrality);
    }
  } // if constexpr (rs == eTest) {

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T1, typename T2> void FillEventHistograms(...)

//============================================================

void CheckUnderflowAndOverflow()
{
  // Check and bail out if in event and particle histograms there are entries which went to underflow or overflow bins.

  // a) Event histograms 1D;
  // b) Event histograms 2D;
  // c) Particle histograms 1D;
  // d) Particle histograms 2D;
  // e) QA Event histograms 2D;
  // f) QA Particle histograms 2D;
  // g) QA Particle event histograms 2D.

  if (tc.fVerboseForEachParticle) {
    StartFunction(__FUNCTION__);
  }

  // a) Event histograms 1D:
  for (int t = 0; t < eEventHistograms_N; t++) // type, see enum eEventHistograms
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!eh.fEventHistograms[t][rs][ba]) {
          continue;
        } else if (eh.fEventHistograms[t][rs][ba]->GetBinContent(0) > 0) {
          LOGF(fatal, "\033[1;31m%s at line %d : underflow in fEventHistograms[%d][%d][%d] => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba);
        } else if (eh.fEventHistograms[t][rs][ba]->GetBinContent(eh.fEventHistograms[t][rs][ba]->GetNbinsX() + 1) > 0) {
          LOGF(fatal, "\033[1;31m%s at line %d : overflow in fEventHistograms[%d][%d][%d] => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba);
        }
      }
    }
  }

  // b) Event histograms 2D:
  // ...

  // c) Particle histograms 1D:
  for (int t = 0; t < eParticleHistograms_N; t++) // type, see enum eParticleHistograms
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!ph.fParticleHistograms[t][rs][ba]) {
          continue;
        } else if (ph.fParticleHistograms[t][rs][ba]->GetBinContent(0) > 0) {
          LOGF(fatal, "\033[1;31m%s at line %d : underflow in fParticleHistograms[%d][%d][%d] => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba);
        } else if (ph.fParticleHistograms[t][rs][ba]->GetBinContent(ph.fParticleHistograms[t][rs][ba]->GetNbinsX() + 1) > 0) {
          LOGF(fatal, "\033[1;31m%s at line %d : overflow in fParticleHistograms[%d][%d][%d] => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba);
        }
      }
    }
  }

  // d) Particle histograms 2D:
  for (int t = 0; t < eParticleHistograms2D_N; t++) // type, see enum eParticleHistograms2D
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!ph.fParticleHistograms2D[t][rs][ba]) {
          continue;
        }

        // Underflow and overflow in x:
        for (int binY = 0; binY <= ph.fParticleHistograms2D[t][rs][ba]->GetNbinsY(); binY++) {
          if (ph.fParticleHistograms2D[t][rs][ba]->GetBinContent(ph.fParticleHistograms2D[t][rs][ba]->GetBin(0, binY)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : underflow in x variable in fParticleHistograms2D[%d][%d][%d], for binY = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binY);
          }
          if (ph.fParticleHistograms2D[t][rs][ba]->GetBinContent(ph.fParticleHistograms2D[t][rs][ba]->GetBin(ph.fParticleHistograms2D[t][rs][ba]->GetNbinsX() + 1, binY)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : overflow in x variable in fParticleHistograms2D[%d][%d][%d], for binY = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binY);
          }
        } // for (int binY = 0; binY <= ph.fParticleHistograms2D[t][rs][ba]->GetNbinsY(); binY++) {

        // Underflow and overflow in y:
        for (int binX = 0; binX <= ph.fParticleHistograms2D[t][rs][ba]->GetNbinsX(); binX++) {
          if (ph.fParticleHistograms2D[t][rs][ba]->GetBinContent(ph.fParticleHistograms2D[t][rs][ba]->GetBin(binX, 0)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : underflow in y variable in fParticleHistograms2D[%d][%d][%d], for binX = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binX);
          }

          if (ph.fParticleHistograms2D[t][rs][ba]->GetBinContent(ph.fParticleHistograms2D[t][rs][ba]->GetBin(binX, ph.fParticleHistograms2D[t][rs][ba]->GetNbinsY() + 1)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : overflow in y variable in fParticleHistograms2D[%d][%d][%d], for binX = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binX);
          }
        } // for (int binX = 0; binX <= ph.fParticleHistograms2D[t][rs][ba]->GetNbinsX(); binX++) {
      } // for (int ba = 0; ba < 2; ba++) // before/after cuts
    } // for (int rs = 0; rs < 2; rs++) // reco/sim
  } // for (int t = 0; t < eParticleHistograms2D_N; t++) // type, see enum eParticleHistograms2D

  // e) QA Event histograms 2D:
  for (int t = 0; t < eQAEventHistograms2D_N; t++) // type, see enum eQAEventHistograms2D
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!qa.fQAEventHistograms2D[t][rs][ba]) {
          continue;
        }

        // Underflow and overflow in x:
        for (int binY = 0; binY <= qa.fQAEventHistograms2D[t][rs][ba]->GetNbinsY(); binY++) {
          if (qa.fQAEventHistograms2D[t][rs][ba]->GetBinContent(qa.fQAEventHistograms2D[t][rs][ba]->GetBin(0, binY)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : underflow in x variable in fEventHistograms2D[%d][%d][%d], for binY = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binY);
          }
          if (qa.fQAEventHistograms2D[t][rs][ba]->GetBinContent(qa.fQAEventHistograms2D[t][rs][ba]->GetBin(qa.fQAEventHistograms2D[t][rs][ba]->GetNbinsX() + 1, binY)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : overflow in x variable in fEventHistograms2D[%d][%d][%d], for binY = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binY);
          }
        } // for (int binY = 0; binY <= qa.fQAEventHistograms2D[t][rs][ba]->GetNbinsY(); binY++) {

        // Underflow and overflow in y:
        for (int binX = 0; binX <= qa.fQAEventHistograms2D[t][rs][ba]->GetNbinsX(); binX++) {
          if (qa.fQAEventHistograms2D[t][rs][ba]->GetBinContent(qa.fQAEventHistograms2D[t][rs][ba]->GetBin(binX, 0)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : underflow in y variable in fEventHistograms2D[%d][%d][%d], for binX = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binX);
          }

          if (qa.fQAEventHistograms2D[t][rs][ba]->GetBinContent(qa.fQAEventHistograms2D[t][rs][ba]->GetBin(binX, qa.fQAEventHistograms2D[t][rs][ba]->GetNbinsY() + 1)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : overflow in y variable in fEventHistograms2D[%d][%d][%d], for binX = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binX);
          }
        } // for (int binX = 0; binX <= qa.fQAEventHistograms2D[t][rs][ba]->GetNbinsX(); binX++) {
      } // for (int ba = 0; ba < 2; ba++) // before/after cuts
    } // for (int rs = 0; rs < 2; rs++) // reco/sim
  } // for (int t = 0; t < eQAEventHistograms2D_N; t++) // type, see enum eQAEventHistograms2D

  // f) QA Particle histograms 2D:
  for (int t = 0; t < eQAParticleHistograms2D_N; t++) // type, see enum eQAParticleHistograms2D
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!qa.fQAParticleHistograms2D[t][rs][ba]) {
          continue;
        }

        // Underflow and overflow in x:
        for (int binY = 0; binY <= qa.fQAParticleHistograms2D[t][rs][ba]->GetNbinsY(); binY++) {
          if (qa.fQAParticleHistograms2D[t][rs][ba]->GetBinContent(qa.fQAParticleHistograms2D[t][rs][ba]->GetBin(0, binY)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : underflow in x variable in fParticleHistograms2D[%d][%d][%d], for binY = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binY);
          }
          if (qa.fQAParticleHistograms2D[t][rs][ba]->GetBinContent(qa.fQAParticleHistograms2D[t][rs][ba]->GetBin(qa.fQAParticleHistograms2D[t][rs][ba]->GetNbinsX() + 1, binY)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : overflow in x variable in fParticleHistograms2D[%d][%d][%d], for binY = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binY);
          }
        } // for (int binY = 0; binY <= qa.fQAParticleHistograms2D[t][rs][ba]->GetNbinsY(); binY++) {

        // Underflow and overflow in y:
        for (int binX = 0; binX <= qa.fQAParticleHistograms2D[t][rs][ba]->GetNbinsX(); binX++) {
          if (qa.fQAParticleHistograms2D[t][rs][ba]->GetBinContent(qa.fQAParticleHistograms2D[t][rs][ba]->GetBin(binX, 0)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : underflow in y variable in fParticleHistograms2D[%d][%d][%d], for binX = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binX);
          }

          if (qa.fQAParticleHistograms2D[t][rs][ba]->GetBinContent(qa.fQAParticleHistograms2D[t][rs][ba]->GetBin(binX, qa.fQAParticleHistograms2D[t][rs][ba]->GetNbinsY() + 1)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : overflow in y variable in fParticleHistograms2D[%d][%d][%d], for binX = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binX);
          }
        } // for (int binX = 0; binX <= qa.fQAParticleHistograms2D[t][rs][ba]->GetNbinsX(); binX++) {
      } // for (int ba = 0; ba < 2; ba++) // before/after cuts
    } // for (int rs = 0; rs < 2; rs++) // reco/sim
  } // for (int t = 0; t < eQAParticleHistograms2D_N; t++) // type, see enum eParticleHistograms2D

  // g) QA Particle event histograms 2D:
  //    TBI 20241212 I never validated this code block
  for (int t = 0; t < eQAParticleEventHistograms2D_N; t++) // type, see enum eQAParticleEventHistograms2D
  {
    for (int rs = 0; rs < 2; rs++) // reco/sim
    {
      for (int ba = 0; ba < 2; ba++) // before/after cuts
      {
        if (!qa.fQAParticleEventHistograms2D[t][rs][ba]) {
          continue;
        }

        // Underflow and overflow in x:
        for (int binY = 0; binY <= qa.fQAParticleEventHistograms2D[t][rs][ba]->GetNbinsY(); binY++) {
          if (qa.fQAParticleEventHistograms2D[t][rs][ba]->GetBinContent(qa.fQAParticleEventHistograms2D[t][rs][ba]->GetBin(0, binY)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : underflow in x variable in fParticleEventHistograms2D[%d][%d][%d], for binY = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binY);
          }
          if (qa.fQAParticleEventHistograms2D[t][rs][ba]->GetBinContent(qa.fQAParticleEventHistograms2D[t][rs][ba]->GetBin(qa.fQAParticleEventHistograms2D[t][rs][ba]->GetNbinsX() + 1, binY)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : overflow in x variable in fParticleEventHistograms2D[%d][%d][%d], for binY = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binY);
          }
        } // for (int binY = 0; binY <= qa.fQAParticleEventHistograms2D[t][rs][ba]->GetNbinsY(); binY++) {

        // Underflow and overflow in y:
        for (int binX = 0; binX <= qa.fQAParticleEventHistograms2D[t][rs][ba]->GetNbinsX(); binX++) {
          if (qa.fQAParticleEventHistograms2D[t][rs][ba]->GetBinContent(qa.fQAParticleEventHistograms2D[t][rs][ba]->GetBin(binX, 0)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : underflow in y variable in fParticleEventHistograms2D[%d][%d][%d], for binX = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binX);
          }

          if (qa.fQAParticleEventHistograms2D[t][rs][ba]->GetBinContent(qa.fQAParticleEventHistograms2D[t][rs][ba]->GetBin(binX, qa.fQAParticleEventHistograms2D[t][rs][ba]->GetNbinsY() + 1)) > 0) {
            LOGF(fatal, "\033[1;31m%s at line %d : overflow in y variable in fParticleEventHistograms2D[%d][%d][%d], for binX = %d  => optimize default binning for this histogram\033[0m", __FUNCTION__, __LINE__, t, rs, ba, binX);
          }
        } // for (int binX = 0; binX <= qa.fQAParticleEventHistograms2D[t][rs][ba]->GetNbinsX(); binX++) {
      } // for (int ba = 0; ba < 2; ba++) // before/after cuts
    } // for (int rs = 0; rs < 2; rs++) // reco/sim
  } // for (int t = 0; t < eQAParticleEventHistograms2D_N; t++) // type, see enum eParticleEventHistograms2D

  if (tc.fVerboseForEachParticle) {
    ExitFunction(__FUNCTION__);
  }

} // void CheckUnderflowAndOverflow()

//============================================================

template <eRecSim rs, typename T>
bool ValidTrack(T const& track)
{
  // Before I start applying any track cuts, check if this is a valid track.
  // For instance, Run 2 or Run 1 tracklets are NOT valid tracks, as they carry no pt information, and in this function they are filtered out.

  // See enum TrackTypeEnum in O2/Framework/Core/include/Framework/DataTypes.h for further info.

  // a) Validity checks for tracks in Run 3;
  // b) Validity checks for tracks in Run 2 and 1;
  // c) Additional validity checks for all tracks (in Run 3, 2 and 1), use only during debugging.

  if (tc.fVerboseForEachParticle) {
    StartFunction(__FUNCTION__);
    LOGF(info, "  track.phi() = %f", track.phi());
    LOGF(info, "  track.pt()  = %f", track.pt());
    LOGF(info, "  track.eta() = %f", track.eta());
    // LOGF(info, "track.trackType() = %d", static_cast<int>(track.trackType())); TBI 20240404 this is not supported for MC particles
  }

  // a) Validity checks for tracks in Run 3:
  // *) Ensure that I am taking into account propagated tracks (and not e.g. track evaluated at innermost update):
  if constexpr (rs == eRec || rs == eRecAndSim) {
    if (!(track.trackType() == o2::aod::track::TrackTypeEnum::Track)) {
      if (tc.fVerboseForEachParticle) {
        LOGF(info, "\033[1;31m%s track.trackType() == o2::aod::track::TrackTypeEnum::Trac\033[0m", __FUNCTION__);
      }
      return false;
    }
  }

  // b) Validity checks for tracks in Run 2 and 1:
  // *) Ensure that tracklets (no pt information) are skipped:
  if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    if (!(track.trackType() == o2::aod::track::TrackTypeEnum::Run2Track)) {
      if (tc.fVerboseForEachParticle) {
        LOGF(info, "\033[1;31m%s track.trackType() == o2::aod::track::TrackTypeEnum::Run2Track\033[0m", __FUNCTION__);
      }
      return false;
    }
  }

  // *) Temporary here, until I cover also these cases:
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {
    LOGF(fatal, "\033[1;31m%s at line %d : add support for TrackTypeEnum here also for cases eSim, eSim_Run2 and eSim_Run1\033[0m", __FUNCTION__, __LINE__);
  }

  // c) Additional validity checks for all tracks (in Run 3, 2 and 1), use only during debugging:
  if (tc.fInsanityCheckForEachParticle) {

    // *) std::isnan() check (remember that 'nan' is 0./0., inf-inf, etc. However 'inf' itself is NOT a 'nan', therefore std::isnan(1./0.) is false, std::isnan(0./0.) is true, etc.):
    if (std::isnan(track.phi()) || std::isnan(track.pt()) || std::isnan(track.eta())) {
      if (tc.fVerboseForEachParticle) {
        LOGF(info, "\033[1;31m%s std::isnan(track.phi()) || std::isnan(track.pt()) || std::isnan(track.eta())\033[0m", __FUNCTION__);
        LOGF(error, "track.phi() = %f\ntrack.pt() = %f\ntrack.eta() = %f", track.phi(), track.pt(), track.eta());
      }
      return false;
    }

    // *) ...
    // ...

  } // if(tc.fInsanityCheckForEachParticle) {

  // *) All checks above survived, then it's a valid track:
  return true;

} // template <eRecSim rs, typename T> bool ValidTrack(T const& track)

//============================================================

template <eRecSim rs, typename T>
void ParticleCutsCounters(T const& track)
{
  // Use this function to fill absolute and sequential particle cut counters. Use only during QA, as this is computationally heavy (I mean really).

  if (tc.fVerboseForEachParticle) {
    StartFunction(__FUNCTION__);
  }

  // *) Establish ordering of binning in particle cut counters histograms, which resembles ordering of particle cuts implementation:
  if (!pc.fParticleCutCounterBinLabelingIsDone) {
    pc.fParticleCutCounterBinNumber[eRec] = 1; // remember that I cannot use 'rs' here as an index, because the enum eRecSim covers separately Run 1 and Run 2, etc.
    pc.fParticleCutCounterBinNumber[eSim] = 1;
    ParticleCuts<rs>(track, eCutCounterBinning); // dry call, to establish the map fParticleCutCounterMap and its inverse

    // **) Map this ordering into bin labels of actual histograms for particle cut counters:
    for (int rec_sim = 0; rec_sim < 2; rec_sim++) // reco/sim => I use here exceptionally different var 'rec_sim', not the shadow 'rs' in the template parameter
    {
      for (int cc = 0; cc < eCutCounter_N; cc++) // enum eCutCounter
      {
        if (!pc.fParticleCutCounterHist[rec_sim][cc]) {
          continue;
        }
        for (int bin = 1; bin < pc.fParticleCutCounterBinNumber[rec_sim]; bin++) // implemented and used particle cuts in this analysis
        {
          pc.fParticleCutCounterHist[rec_sim][cc]->GetXaxis()->SetBinLabel(bin, FancyFormatting(pc.fParticleCutName[pc.fParticleCutCounterMap[rec_sim]->GetValue(bin)].Data()));
        }
        for (int bin = pc.fParticleCutCounterBinNumber[rec_sim]; bin <= eParticleCuts_N; bin++) // implemented, but unused particle cuts in this analysis
        {
          pc.fParticleCutCounterHist[rec_sim][cc]->GetXaxis()->SetBinLabel(bin, Form("binNo = %d (unused cut)", bin));
          // Remark: I have to write here something concrete as a bin label, if I leave "TBI" for all bin labels here for cuts which were not used,
          // I get this harmless but annoying warning during merging:
          //   Warning in <TH1Merger::CheckForDuplicateLabels>: Histogram fParticleCutCounterHist[rec][seq] has duplicate labels in the x axis. Bin contents will be merged in a single bin
          // TBI 20241130 as a better solution, I shall re-define this histogram with the narower range on x-axis...
        }
        // All cuts which were implemeted, but not used I simply do not show (i can always UnZoom x-axis in TBrowser, if I want to see 'em).
        pc.fParticleCutCounterHist[rec_sim][cc]->GetXaxis()->SetRangeUser(pc.fParticleCutCounterHist[rec_sim][cc]->GetBinLowEdge(1), pc.fParticleCutCounterHist[rec_sim][cc]->GetBinLowEdge(pc.fParticleCutCounterBinNumber[rec_sim]));
      }
    }
    pc.fParticleCutCounterBinLabelingIsDone = true; // this flag ensures that this specific binning is performed only once, for the first processed particle
    // delete pc.fParticleCutCounterMap[eRec]; // TBI 20240508 if i do not need them later, I could delete here
    // delete pc.fParticleCutCounterMap[eSim];
    // delete pc.fParticleCutCounterMapInverse[eRec];
    // delete pc.fParticleCutCounterMapInverse[eSim];
  } // if (!pc.fParticleCutCounterBinLabelingIsDone) {

  // *) Particle cut counter (absolute):
  if (pc.fUseParticleCutCounterAbsolute) {
    pc.fParticleCutCounterBinNumber[eRec] = 1;
    pc.fParticleCutCounterBinNumber[eSim] = 1;
    ParticleCuts<rs>(track, eCutCounterAbsolute);
  }

  // *) Particle cut counter (sequential):
  if (pc.fUseParticleCutCounterSequential) {
    pc.fParticleCutCounterBinNumber[eRec] = 1;
    pc.fParticleCutCounterBinNumber[eSim] = 1;
    ParticleCuts<rs>(track, eCutCounterSequential);
  }

} // template <eRecSim rs, typename T> void ParticleCutsCounters(T const& track)

//============================================================

template <eRecSim rs, typename T>
bool ParticleCuts(T const& track, eCutModus cutModus)
{
  // Particle cuts on reconstructed and simulated data. Supports particle cut counters, both absolute and sequential.
  // There is also a related enum eParticleCuts.
  // Remark: I have added to all if statemets below which deals with floats, e.g. TMath::Abs(track.eta() - pc.fdParticleCuts[eEta][eMax]) < tc.fFloatingPointPrecision ,
  //           to enforce the ROOT convention: "lower boundary included, upper boundary excluded"

  // a) Particle cuts on reconstructed, and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1);
  // b) Particle cuts only on simulated (common to Run 3, Run 2 and Run 1);
  // c) Particle cuts on reconstructed, and corresponding MC truth simulated (Run 3 specific);
  // d) Particle cuts on simulated (Run 3 specific);
  // e) Particle cuts on reconstructed, and corresponding MC truth simulated (Run 1 and 2 specific); // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  // f) Particle cuts on simulated (Run 1 and 2 specific); // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  // *) Particle cuts on Test case;
  // *) Toy NUA.

  if (tc.fVerboseForEachParticle) {
    StartFunction(__FUNCTION__);
  }

  // a) Particle cuts on reconstructed, and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1) ...
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

    // *) Phi:
    if (pc.fUseParticleCuts[ePhi]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, ePhi, eCutCounterBinning);
      } else if (track.phi() < pc.fdParticleCuts[ePhi][eMin] || track.phi() > pc.fdParticleCuts[ePhi][eMax] || TMath::Abs(track.phi() - pc.fdParticleCuts[ePhi][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, ePhi, cutModus)) {
          return false;
        }
      }
    }

    // *) Pt:
    if (pc.fUseParticleCuts[ePt]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, ePt, eCutCounterBinning);
      } else if (track.pt() < pc.fdParticleCuts[ePt][eMin] || track.pt() > pc.fdParticleCuts[ePt][eMax] || TMath::Abs(track.pt() - pc.fdParticleCuts[ePt][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, ePt, cutModus)) {
          return false;
        }
      }
    }

    // *) Eta:
    if (pc.fUseParticleCuts[eEta]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eEta, eCutCounterBinning);
      } else if (track.eta() < pc.fdParticleCuts[eEta][eMin] || track.eta() > pc.fdParticleCuts[eEta][eMax] || TMath::Abs(track.eta() - pc.fdParticleCuts[eEta][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, eEta, cutModus)) {
          return false;
        }
      }
    }

    // *) Charge:
    if (pc.fUseParticleCuts[eCharge]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eCharge, eCutCounterBinning);
      } else if (0 == track.sign() || track.sign() < pc.fdParticleCuts[eCharge][eMin] || track.sign() > pc.fdParticleCuts[eCharge][eMax]) {
        // With first condition, I always throw away neutral particles.
        // I can use safely == here, because track.sign() returns short int.
        if (!ParticleCut(eRec, eCharge, cutModus)) {
          return false;
        }
      }
    }

    // *) tpcNClsFindable:
    if (pc.fUseParticleCuts[etpcNClsFindable]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etpcNClsFindable, eCutCounterBinning);
      } else if (track.tpcNClsFindable() < pc.fdParticleCuts[etpcNClsFindable][eMin] || track.tpcNClsFindable() > pc.fdParticleCuts[etpcNClsFindable][eMax]) {
        if (!ParticleCut(eRec, etpcNClsFindable, cutModus)) {
          return false;
        }
      }
    }

    // *) tpcNClsShared:
    if (pc.fUseParticleCuts[etpcNClsShared]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etpcNClsShared, eCutCounterBinning);
      } else if (track.tpcNClsShared() < pc.fdParticleCuts[etpcNClsShared][eMin] || track.tpcNClsShared() > pc.fdParticleCuts[etpcNClsShared][eMax]) {
        if (!ParticleCut(eRec, etpcNClsShared, cutModus)) {
          return false;
        }
      }
    }

    // *) itsChi2NCl
    if (pc.fUseParticleCuts[eitsChi2NCl]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eitsChi2NCl, eCutCounterBinning);
      } else if (track.itsChi2NCl() < pc.fdParticleCuts[eitsChi2NCl][eMin] || track.itsChi2NCl() > pc.fdParticleCuts[eitsChi2NCl][eMax]) {
        if (!ParticleCut(eRec, eitsChi2NCl, cutModus)) {
          return false;
        }
      }
    }

    // *) tpcNClsFound:
    if (pc.fUseParticleCuts[etpcNClsFound]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etpcNClsFound, eCutCounterBinning);
      } else if (track.tpcNClsFound() < pc.fdParticleCuts[etpcNClsFound][eMin] || track.tpcNClsFound() > pc.fdParticleCuts[etpcNClsFound][eMax]) {
        if (!ParticleCut(eRec, etpcNClsFound, cutModus)) {
          return false;
        }
      }
    }

    // *) tpcNClsCrossedRows:
    if (pc.fUseParticleCuts[etpcNClsCrossedRows]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etpcNClsCrossedRows, eCutCounterBinning);
      } else if (track.tpcNClsCrossedRows() < pc.fdParticleCuts[etpcNClsCrossedRows][eMin] || track.tpcNClsCrossedRows() > pc.fdParticleCuts[etpcNClsCrossedRows][eMax]) {
        if (!ParticleCut(eRec, etpcNClsCrossedRows, cutModus)) {
          return false;
        }
      }
    }

    // *) itsNCls:
    if (pc.fUseParticleCuts[eitsNCls]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eitsNCls, eCutCounterBinning);
      } else if (track.itsNCls() < pc.fdParticleCuts[eitsNCls][eMin] || track.itsNCls() > pc.fdParticleCuts[eitsNCls][eMax]) {
        if (!ParticleCut(eRec, eitsNCls, cutModus)) {
          return false;
        }
      }
    }

    // *) itsNClsInnerBarrel:
    if (pc.fUseParticleCuts[eitsNClsInnerBarrel]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eitsNClsInnerBarrel, eCutCounterBinning);
      } else if (track.itsNClsInnerBarrel() < pc.fdParticleCuts[eitsNClsInnerBarrel][eMin] || track.itsNClsInnerBarrel() > pc.fdParticleCuts[eitsNClsInnerBarrel][eMax]) {
        if (!ParticleCut(eRec, eitsNClsInnerBarrel, cutModus)) {
          return false;
        }
      }
    }

    // *) tpcCrossedRowsOverFindableCls:
    if (pc.fUseParticleCuts[etpcCrossedRowsOverFindableCls]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etpcCrossedRowsOverFindableCls, eCutCounterBinning);
      } else if (track.tpcCrossedRowsOverFindableCls() < pc.fdParticleCuts[etpcCrossedRowsOverFindableCls][eMin] || track.tpcCrossedRowsOverFindableCls() > pc.fdParticleCuts[etpcCrossedRowsOverFindableCls][eMax] || TMath::Abs(track.tpcCrossedRowsOverFindableCls() - pc.fdParticleCuts[etpcCrossedRowsOverFindableCls][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, etpcCrossedRowsOverFindableCls, cutModus)) {
          return false;
        }
      }
    }

    // *) tpcFoundOverFindableCls:
    if (pc.fUseParticleCuts[etpcFoundOverFindableCls]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etpcFoundOverFindableCls, eCutCounterBinning);
      } else if (track.tpcFoundOverFindableCls() < pc.fdParticleCuts[etpcFoundOverFindableCls][eMin] || track.tpcFoundOverFindableCls() > pc.fdParticleCuts[etpcFoundOverFindableCls][eMax] || TMath::Abs(track.tpcFoundOverFindableCls() - pc.fdParticleCuts[etpcFoundOverFindableCls][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, etpcFoundOverFindableCls, cutModus)) {
          return false;
        }
      }
    }

    // *) tpcFractionSharedCls:
    if (pc.fUseParticleCuts[etpcFractionSharedCls]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etpcFractionSharedCls, eCutCounterBinning);
      } else if (track.tpcFractionSharedCls() < pc.fdParticleCuts[etpcFractionSharedCls][eMin] || track.tpcFractionSharedCls() > pc.fdParticleCuts[etpcFractionSharedCls][eMax] || TMath::Abs(track.tpcFractionSharedCls() - pc.fdParticleCuts[etpcFractionSharedCls][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, etpcFractionSharedCls, cutModus)) {
          return false;
        }
      }
    }

    // *) tpcChi2NCl:
    if (pc.fUseParticleCuts[etpcChi2NCl]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etpcChi2NCl, eCutCounterBinning);
      } else if (track.tpcChi2NCl() < pc.fdParticleCuts[etpcChi2NCl][eMin] || track.tpcChi2NCl() > pc.fdParticleCuts[etpcChi2NCl][eMax] || TMath::Abs(track.tpcChi2NCl() - pc.fdParticleCuts[etpcChi2NCl][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, etpcChi2NCl, cutModus)) {
          return false;
        }
      }
    }

    // *) dcaXY:
    if (pc.fUseParticleCuts[edcaXY]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, edcaXY, eCutCounterBinning);
      } else if (track.dcaXY() < pc.fdParticleCuts[edcaXY][eMin] || track.dcaXY() > pc.fdParticleCuts[edcaXY][eMax] || TMath::Abs(track.dcaXY() - pc.fdParticleCuts[edcaXY][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, edcaXY, cutModus)) {
          return false;
        }
      }
    }

    // *) dcaZ:
    if (pc.fUseParticleCuts[edcaZ]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, edcaZ, eCutCounterBinning);
      } else if (track.dcaZ() < pc.fdParticleCuts[edcaZ][eMin] || track.dcaZ() > pc.fdParticleCuts[edcaZ][eMax] || TMath::Abs(track.dcaZ() - pc.fdParticleCuts[edcaZ][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, edcaZ, cutModus)) {
          return false;
        }
      }
    }

    // *) trackCutFlag:
    if (pc.fUseParticleCuts[etrackCutFlag]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etrackCutFlag, eCutCounterBinning);
      } else if (!track.trackCutFlag()) {
        if (!ParticleCut(eRec, etrackCutFlag, cutModus)) {
          return false;
        }
      }
    }

    // *) trackCutFlagFb1:
    if (pc.fUseParticleCuts[etrackCutFlagFb1]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etrackCutFlagFb1, eCutCounterBinning);
      } else if (!track.trackCutFlagFb1()) {
        if (!ParticleCut(eRec, etrackCutFlagFb1, cutModus)) {
          return false;
        }
      }
    }

    // *) trackCutFlagFb2:
    if (pc.fUseParticleCuts[etrackCutFlagFb2]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, etrackCutFlagFb2, eCutCounterBinning);
      } else if (!track.trackCutFlagFb2()) {
        if (!ParticleCut(eRec, etrackCutFlagFb2, cutModus)) {
          return false;
        }
      }
    }

    // *) isQualityTrack:
    if (pc.fUseParticleCuts[eisQualityTrack]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eisQualityTrack, eCutCounterBinning);
      } else if (!track.isQualityTrack()) {
        if (!ParticleCut(eRec, eisQualityTrack, cutModus)) {
          return false;
        }
      }
    }

    // *) isPrimaryTrack:
    if (pc.fUseParticleCuts[eisPrimaryTrack]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eisPrimaryTrack, eCutCounterBinning);
      } else if (!track.isPrimaryTrack()) {
        if (!ParticleCut(eRec, eisPrimaryTrack, cutModus)) {
          return false;
        }
      }
    }

    // *) isInAcceptanceTrack:
    if (pc.fUseParticleCuts[eisInAcceptanceTrack]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eisInAcceptanceTrack, eCutCounterBinning);
      } else if (!track.isInAcceptanceTrack()) {
        if (!ParticleCut(eRec, eisInAcceptanceTrack, cutModus)) {
          return false;
        }
      }
    }

    // *) isGlobalTrack:
    if (pc.fUseParticleCuts[eisGlobalTrack]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eisGlobalTrack, eCutCounterBinning);
      } else if (!track.isGlobalTrack()) {
        if (!ParticleCut(eRec, eisGlobalTrack, cutModus)) {
          return false;
        }
      }
    }

    // *) isPVContributor:
    if (pc.fUseParticleCuts[eisPVContributor]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eisPVContributor, eCutCounterBinning);
      } else if (!track.isPVContributor()) {
        if (!ParticleCut(eRec, eisPVContributor, cutModus)) {
          return false;
        }
      }
    }

    // *) PtDependentDCAxyParameterization:
    if (pc.fUseParticleCuts[ePtDependentDCAxyParameterization]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, ePtDependentDCAxyParameterization, eCutCounterBinning);
      } else if (TMath::Abs(track.dcaXY()) > pc.fPtDependentDCAxyFormula->Eval(track.pt())) {
        if (!ParticleCut(eRec, ePtDependentDCAxyParameterization, cutModus)) {
          return false;
        }
      }
    }

    // ...

    // ... and corresponding MC truth simulated:
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
    // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
    if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {

      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for this track, skip...");
        return false; // TBI 20231107 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this track
      }
      // auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle

      // In this branch I can cut additionally and directly on corresponding MC truth simulated, e.g. on mcparticle.pt()
      // In case I implement something here, remember to switch from eRec to eSim when calling e.g. ParticleCut(...)

      /*
          // *) Phi: TBI 2024-511 re-think if i really cut directly on MC truth kine and other info and keep it in sync with what I did in AliPhysics
          if (pc.fUseParticleCuts[ePhi]) {
            if (cutModus == eCutCounterBinning) {
              ParticleCut(eSim, ePhi, eCutCounterBinning);
            } else if (mcparticle.phi() < pc.fdParticleCuts[ePhi][eMin] || mcparticle.phi() > pc.fdParticleCuts[ePhi][eMax]) {
              if (!ParticleCut(eSim, ePhi, cutModus)) {
                return false;
              }
            }
          }
      */
      // *) Charge: TBI 20240511 mcparticle.sign() doesn't exist, here most likely i need to cut on the signature of mcparticle.pdg() but check further, because e is negative charge, but PDG is 11, etc.
      /*
            if (pc.fUseParticleCuts[eCharge]) {
              if (cutModus == eCutCounterBinning) {
                ParticleCut(eSim, eCharge, eCutCounterBinning);
              } else if (0 == mcparticle.sign() || mcparticle.sign() < pc.fdParticleCuts[eCharge][eMin] || mcparticle.sign() > pc.fdParticleCuts[eCharge][eMax]) {
                // TBI 20240511 with first condition, I always throw away neutral particles, so for the time being that is hardcoded
                if (!ParticleCut(eSim, eCharge, cutModus)) {
                  return false;
                }
              }
            }
      */
      // TBI 20240511 add cut on PDG

      // ...

    } // if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {

  } // if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // -------------------------------------------------------------------------

  // b) Particle cuts only on simulated (common to Run 3, Run 2 and Run 1):
  //    Remark #1: This branch is relevant when processing ONLY simulated data at generator level.
  //    Remark #2: In this branch, 'track' is always TracksSim = aod::McParticles, see https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

    // *) Phi:
    if (pc.fUseParticleCuts[ePhi]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eSim, ePhi, eCutCounterBinning);
      } else if (track.phi() < pc.fdParticleCuts[ePhi][eMin] || track.phi() > pc.fdParticleCuts[ePhi][eMax] || TMath::Abs(track.phi() - pc.fdParticleCuts[ePhi][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eSim, ePhi, cutModus)) {
          return false;
        }
      }
    }

    // *) Pt:
    if (pc.fUseParticleCuts[ePt]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eSim, ePt, eCutCounterBinning);
      } else if (track.pt() < pc.fdParticleCuts[ePt][eMin] || track.pt() > pc.fdParticleCuts[ePt][eMax] || TMath::Abs(track.pt() - pc.fdParticleCuts[ePt][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eSim, ePt, cutModus)) {
          return false;
        }
      }
    }

    // *) Eta:
    if (pc.fUseParticleCuts[eEta]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eSim, eEta, eCutCounterBinning);
      } else if (track.eta() < pc.fdParticleCuts[eEta][eMin] || track.eta() > pc.fdParticleCuts[eEta][eMax] || TMath::Abs(track.eta() - pc.fdParticleCuts[eEta][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eSim, eEta, cutModus)) {
          return false;
        }
      }
    }

    /*
        // *) Charge:
        if (pc.fUseParticleCuts[eCharge]) {
          if (cutModus == eCutCounterBinning) {
            ParticleCut(eSim, eCharge, eCutCounterBinning);
          } else if (0 == track.sign() || track.sign() < pc.fdParticleCuts[eCharge][eMin] || track.sign() > pc.fdParticleCuts[eCharge][eMax]) {
            // TBI 20240511 with first condition, I always throw away neutral particles, so for the time being that is hardcoded
            if (!ParticleCut(eSim, eCharge, cutModus)) {
              return false;
            }
          }
        }
    */
    // TBI 20240511 add cut on PDG

    // ...

  } // if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

  // -------------------------------------------------------------------------

  // c) Particle cuts on reconstructed, and corresponding MC truth simulated (Run 3 specific):
  //    Remark: I implement here only the particle cuts which are not already in group a) above, and which make sense only for Run 3 data.
  if constexpr (rs == eRec || rs == eRecAndSim) {

    // ...

    // ... and corresponding MC truth simulated (Run 3 specific):
    if constexpr (rs == eRecAndSim) {

      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for this track, skip...");
        return false; // TBI 20231107 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this track
      }
      // auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle

      // ...

    } // if constexpr (rs == eRecAndSim) {

  } // if constexpr (rs == eRec || rs == eRecAndSim) {

  // -------------------------------------------------------------------------

  // d) Particle cuts on simulated (Run 3 specific):
  //    Remark #1: I implement here only the particle cuts which are not already in group b) above, and which make sense only for Run 3 data.
  //    Remark #2: In this branch, 'track' is always TracksSim = aod::McParticles, see https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
  //               See how I handled the case b) above.
  if constexpr (rs == eSim) {

    // ...

  } // if constexpr (rs == eSim) {

  // -------------------------------------------------------------------------

  // e) Particle cuts on reconstructed, and corresponding MC truth simulated (Run 1 and 2 specific):
  //    Remark: I implement here only the particle cuts which are not already in group a) above, and which make sense only for Run 1 and 2 data.
  if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

    // ...

    // ... and corresponding MC truth simulated:
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
    // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
    if constexpr (rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {

      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for this track, skip...");
        return false; // TBI 20231107 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this track
      }
      // auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle

      // ...

    } // if constexpr (rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {

  } // if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // -------------------------------------------------------------------------

  // f) Event cuts on simulated (Run 1 and 2 specific)
  //    Remark #1: I implement here only the event cuts which are not already in group b) above, and which make sense only for Run 1 and 2 data.
  //    Remark #2: In this branch, 'track' is always TracksSim = aod::McParticles, see https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
  //               See how I handled the case b) above.
  if constexpr (rs == eSim_Run2 || rs == eSim_Run1) {

    // ...

  } // if constexpr (rs == eSim_Run2 || rs == eSim_Run1) {

  // -------------------------------------------------------------------------

  // *) Test case:
  if constexpr (rs == eTest) {
    // This branch corresponds to process with minimal subscription - I implement just a few example cuts, just for testing purposes.
    // Only eRec is support in Test for the time being.

    // *) Phi:
    if (pc.fUseParticleCuts[ePhi]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, ePhi, eCutCounterBinning);
      } else if (track.phi() < pc.fdParticleCuts[ePhi][eMin] || track.phi() > pc.fdParticleCuts[ePhi][eMax] || TMath::Abs(track.phi() - pc.fdParticleCuts[ePhi][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, ePhi, cutModus)) {
          return false;
        }
      }
    }

    // *) Pt:
    if (pc.fUseParticleCuts[ePt]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, ePt, eCutCounterBinning);
      } else if (track.pt() < pc.fdParticleCuts[ePt][eMin] || track.pt() > pc.fdParticleCuts[ePt][eMax] || TMath::Abs(track.pt() - pc.fdParticleCuts[ePt][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, ePt, cutModus)) {
          return false;
        }
      }
    }

    // *) Eta:
    if (pc.fUseParticleCuts[eEta]) {
      if (cutModus == eCutCounterBinning) {
        ParticleCut(eRec, eEta, eCutCounterBinning);
      } else if (track.eta() < pc.fdParticleCuts[eEta][eMin] || track.eta() > pc.fdParticleCuts[eEta][eMax] || TMath::Abs(track.eta() - pc.fdParticleCuts[eEta][eMax]) < tc.fFloatingPointPrecision) {
        if (!ParticleCut(eRec, eEta, cutModus)) {
          return false;
        }
      }
    }

    // ...

  } // if constexpr (rs == eTest) {

  // -------------------------------------------------------------------------

  // *) Toy NUA:
  if (nua.fApplyNUAPDF[ePhiNUAPDF] || nua.fApplyNUAPDF[ePtNUAPDF] || nua.fApplyNUAPDF[eEtaNUAPDF]) {

    // Remark: I do not for the time being add Toy NUA cuts to particle cut counters, since in this case I can inspect direcly from phi, pt and eta distributions.

    // Local kine variables on which support for Toy NUA is implemented and applied:
    double dPhi = 0.;
    double dPt = 0.;
    double dEta = 0.;

    // *) Apply Toy NUA on info available in reconstructed (and the corresponding MC truth simulated track);
    if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
      dPhi = track.phi();
      dPt = track.pt();
      dEta = track.eta();

      // Apply NUA on these kine variables:
      if (nua.fApplyNUAPDF[ePhiNUAPDF] && !Accept(dPhi, ePhiNUAPDF)) {
        return false;
      }
      if (nua.fApplyNUAPDF[ePtNUAPDF] && !Accept(dPt, ePtNUAPDF)) {
        return false;
      }
      if (nua.fApplyNUAPDF[eEtaNUAPDF] && !Accept(dEta, eEtaNUAPDF)) {
        return false;
      }

      // ... and corresponding MC truth simulated ( see https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
      if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC particle for this track, skip...");
          return false; // TBI 20231107 re-think. I shouldn't probably get to this point, if MC truth info doesn't exist for this particle
        }
        auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle
        dPhi = mcparticle.phi();
        dPt = mcparticle.pt();
        dEta = mcparticle.eta();

        // Apply NUA on these kine variables:
        if (nua.fApplyNUAPDF[ePhiNUAPDF] && !Accept(dPhi, ePhiNUAPDF)) {
          return false;
        }
        if (nua.fApplyNUAPDF[ePtNUAPDF] && !Accept(dPt, ePtNUAPDF)) {
          return false;
        }
        if (nua.fApplyNUAPDF[eEtaNUAPDF] && !Accept(dEta, eEtaNUAPDF)) {
          return false;
        }

      } // if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
    } // if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

    // *) Apply Toy NUA on info available only in simulated data:
    if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {
      // Remark: in this branch, 'track' is always TracksSim = aod::McParticles
      dPhi = track.phi();
      dPt = track.pt();
      dEta = track.eta();

      // Apply NUA on these kine variables:
      if (nua.fApplyNUAPDF[ePhiNUAPDF] && !Accept(dPhi, ePhiNUAPDF)) {
        return false;
      }
      if (nua.fApplyNUAPDF[ePtNUAPDF] && !Accept(dPt, ePtNUAPDF)) {
        return false;
      }
      if (nua.fApplyNUAPDF[eEtaNUAPDF] && !Accept(dEta, eEtaNUAPDF)) {
        return false;
      }
    } // if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

  } // if(nua.fApplyNUAPDF[ePhiNUAPDF] || nua.fApplyNUAPDF[ePtNUAPDF] || nua.fApplyNUAPDF[eEtaNUAPDF]) {

  return true;

} // template <eRecSim rs, typename T> bool ParticleCuts(T const& track, eCutModus cutModus)

//============================================================

bool ParticleCut(int rs, int particleCut, eCutModus cutModus)
{
  // Helper function to reduce code bloat in ParticleCuts(). It's meant to be used only in ParticleCuts().

  // Remark: Remember that as a second argument I cannot use enum eParticleCuts, because here in one go I take both enum eParticleCuts and enum eParticleHistograms .

  switch (cutModus) {
    case eCut:
      if (tc.fVerboseForEachParticle) {
        LOGF(info, "\033[1;31mParticle didn't pass the cut: %s\033[0m", pc.fParticleCutName[particleCut].Data());
      }
      return false;
      break;
    case eCutCounterBinning:
      pc.fParticleCutCounterMap[rs]->Add(pc.fParticleCutCounterBinNumber[rs], particleCut);
      pc.fParticleCutCounterMapInverse[rs]->Add(particleCut, pc.fParticleCutCounterBinNumber[rs]);
      pc.fParticleCutCounterBinNumber[rs]++; // yes
      return true;
      break;
    case eCutCounterAbsolute:
      pc.fParticleCutCounterHist[rs][eAbsolute]->Fill(pc.fParticleCutCounterMapInverse[rs]->GetValue(particleCut));
      return true; // yes, so that I can proceed with another cut in ParticleCuts
      break;
    case eCutCounterSequential:
      pc.fParticleCutCounterHist[rs][eSequential]->Fill(pc.fParticleCutCounterMapInverse[rs]->GetValue(particleCut));
      return false; // yes, so that I bail out from ParticleCuts
      break;
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : This cutModus = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(cutModus));
      break;
  } // switch(cutModus)

  return false; // obsolete, but it suppresses the warning...

} // bool ParticleCut(int rs, int particleCut, eCutModus cutModus)

//============================================================

template <eRecSim rs, typename T>
void FillParticleHistograms(T const& track, eBeforeAfter ba, int weight = 1)
{
  // Fill all particle histograms for reconstructed and simulated data.

  // a) Fill reconstructed, and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1);
  // b) Fill only simulated (common to Run 3, Run 2 and Run 1);
  // c) Fill reconstructed, and corresponding MC truth simulated (Run 3 specific);
  // d) Fill only simulated (Run 3 specific);
  // e) Fill reconstructed, and corresponding MC truth simulated (Run 1 and 2 specific); // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  // f) Fill only simulated (Run 1 and 2 specific); // In case there is some corner case between Run 1 and Run 2, simply branch further this one
  // i) Test case.

  // Remark #1: Why weight is introduced as 3rd argument, see explanation in BanishmentLoopOverParticles.
  //            The efficiency loss for default fill between weight = 1 and previous implementation with no weight, is negligible.
  // Remark #2: After calling BanishmentLoopOverParticles, the GetMean(), GetRMS(), GetStdDev(), skewness, kurtosis, etc., of histogram are unaffected.
  //            But GetMeanErorr(), GetRMSError(), etc. are affected.
  //            For instance, GetMean() remains the same, because (x+y)/(1+1) = (x+y+z-z)(1+1+1-1). But whenever weight in the formula is taken directly to some higher power,
  //            like in the calculation of GetMeanError(), this idea with BanishmentLoopOverParticles is not applicable (also when I enable Setw2() in histograms).
  //            Since from particle histograms I only care about the number of entries, I rarely need even GetMean(), and basically never GetMeanError(),
  //            I use BanishmentLoopOverParticles . Alternatively, I would need new set of histograms, fill them separately, etc.

  if (tc.fVerboseForEachParticle) {
    StartFunction(__FUNCTION__);
  }

  if (tc.fInsanityCheckForEachParticle) {
    if (1 != TMath::Abs(weight)) {
      LOGF(fatal, "\033[1;31m%s at line %d : in the current implementation, weight for particle histograms can be only +1 or -1, weight = %d\033[0m", __FUNCTION__, __LINE__, weight);
    }
  }

  // a) Fill reconstructed ... (common to Run 3, Run 2 and Run 1):
  if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    // Remark: Remember to use only eRec and eSim as array indices in histos, also for rs == eRecAndSim, etc. TBI 20240504 shall I introduce generic enum egRec and egSim for this sake?
    // TBI 20240414 also here have to hardcode 'eRec', because 'rs' spans over all enums in eRecSim => I definitely need 'generic Rec' case, perhaps via TExMap ?
    //              But I have already tc.fProcess[eGenericRec] and tc.fProcess[eGenericRecSim], available, shall I simply re-use them?

    // 1D:
    if (ph.fFillParticleHistograms) {
      // From o2::aod::Tracks
      !ph.fParticleHistograms[ePhi][eRec][ba] ? true : ph.fParticleHistograms[ePhi][eRec][ba]->Fill(track.phi(), weight);
      !ph.fParticleHistograms[ePt][eRec][ba] ? true : ph.fParticleHistograms[ePt][eRec][ba]->Fill(track.pt(), weight);
      !ph.fParticleHistograms[eEta][eRec][ba] ? true : ph.fParticleHistograms[eEta][eRec][ba]->Fill(track.eta(), weight);
      !ph.fParticleHistograms[eCharge][eRec][ba] ? true : ph.fParticleHistograms[eCharge][eRec][ba]->Fill(track.sign(), weight);

      // From o2::aod::TracksExtra_001
      !ph.fParticleHistograms[etpcNClsFindable][eRec][ba] ? true : ph.fParticleHistograms[etpcNClsFindable][eRec][ba]->Fill(track.tpcNClsFindable(), weight);
      !ph.fParticleHistograms[etpcNClsShared][eRec][ba] ? true : ph.fParticleHistograms[etpcNClsShared][eRec][ba]->Fill(track.tpcNClsShared(), weight);
      !ph.fParticleHistograms[eitsChi2NCl][eRec][ba] ? true : ph.fParticleHistograms[eitsChi2NCl][eRec][ba]->Fill(track.itsChi2NCl(), weight);
      !ph.fParticleHistograms[etpcNClsFound][eRec][ba] ? true : ph.fParticleHistograms[etpcNClsFound][eRec][ba]->Fill(track.tpcNClsFound(), weight);
      !ph.fParticleHistograms[etpcNClsCrossedRows][eRec][ba] ? true : ph.fParticleHistograms[etpcNClsCrossedRows][eRec][ba]->Fill(track.tpcNClsCrossedRows(), weight);
      !ph.fParticleHistograms[eitsNCls][eRec][ba] ? true : ph.fParticleHistograms[eitsNCls][eRec][ba]->Fill(track.itsNCls(), weight);
      !ph.fParticleHistograms[eitsNClsInnerBarrel][eRec][ba] ? true : ph.fParticleHistograms[eitsNClsInnerBarrel][eRec][ba]->Fill(track.itsNClsInnerBarrel(), weight);
      !ph.fParticleHistograms[etpcCrossedRowsOverFindableCls][eRec][ba] ? true : ph.fParticleHistograms[etpcCrossedRowsOverFindableCls][eRec][ba]->Fill(track.tpcCrossedRowsOverFindableCls(), weight);
      !ph.fParticleHistograms[etpcFoundOverFindableCls][eRec][ba] ? true : ph.fParticleHistograms[etpcFoundOverFindableCls][eRec][ba]->Fill(track.tpcFoundOverFindableCls(), weight);
      !ph.fParticleHistograms[etpcFractionSharedCls][eRec][ba] ? true : ph.fParticleHistograms[etpcFractionSharedCls][eRec][ba]->Fill(track.tpcFractionSharedCls(), weight);
      !ph.fParticleHistograms[etpcChi2NCl][eRec][ba] ? true : ph.fParticleHistograms[etpcChi2NCl][eRec][ba]->Fill(track.tpcChi2NCl(), weight);

      // From o2::aod::TracksDCA
      // Remark: For this one, in Run 3 workflow I need helper task o2-analysis-track-propagation, while in Run 2 and 1 I need o2-analysis-trackextension .
      !ph.fParticleHistograms[edcaXY][eRec][ba] ? true : ph.fParticleHistograms[edcaXY][eRec][ba]->Fill(track.dcaXY(), weight);
      !ph.fParticleHistograms[edcaZ][eRec][ba] ? true : ph.fParticleHistograms[edcaZ][eRec][ba]->Fill(track.dcaZ(), weight);
    }

    // 2D:
    if (ph.fFillParticleHistograms2D) {
      !ph.fParticleHistograms2D[ePhiPt][eRec][ba] ? true : ph.fParticleHistograms2D[ePhiPt][eRec][ba]->Fill(track.phi(), track.pt(), weight);
      !ph.fParticleHistograms2D[ePhiEta][eRec][ba] ? true : ph.fParticleHistograms2D[ePhiEta][eRec][ba]->Fill(track.phi(), track.eta(), weight);
    } // if (ph.fFillParticleHistograms2D) {

    // nD (THnSparse):
    if (ba == eAfter) { // yes, I feel sparse histograms only AFTER cuts for the time being
      // **) eDWPhi : here the fundamental 0-th axis never to be projected out is "phi"
      if (ph.fBookParticleSparseHistograms[eDWPhi]) {
        // Remark: It is mandatory that ordering in initialization here resembles the ordering in enum eDiffPhiWeights
        double vector[eDiffPhiWeights_N] = {track.phi(), track.pt(), track.eta(), static_cast<double>(track.sign()), ebye.fCentrality, ebye.fVz};
        ph.fParticleSparseHistograms[eDWPhi][eRec]->Fill(vector, weight);
      }
      // **) eDWPt : here the fundamental 0-th axis never to be projected out is "pt"
      if (ph.fBookParticleSparseHistograms[eDWPt]) {
        // Remark: It is mandatory that ordering in initialization here resembles the ordering in enum eDiffPtWeights
        double vector[eDiffPtWeights_N] = {track.pt()};
        ph.fParticleSparseHistograms[eDWPt][eRec]->Fill(vector, weight);
      }
      // **) eDWEta : here the fundamental 0-th axis never to be projected out is "eta"
      if (ph.fBookParticleSparseHistograms[eDWEta]) {
        // Remark: It is mandatory that ordering in initialization here resembles the ordering in enum eDiffEtaWeights
        double vector[eDiffEtaWeights_N] = {track.eta()};
        ph.fParticleSparseHistograms[eDWEta][eRec]->Fill(vector, weight);
      }
    } // if (ba == eAfter) {

    // QA:
    if (qa.fFillQAParticleHistograms2D) {
      !qa.fQAParticleHistograms2D[ePt_vs_dcaXY][eRec][ba] ? true : qa.fQAParticleHistograms2D[ePt_vs_dcaXY][eRec][ba]->Fill(track.pt(), track.dcaXY(), weight);
    }
    if ((qa.fFillQAParticleEventHistograms2D || qa.fFillQACorrelationsVsHistograms2D) && qa.fQAParticleEventProEbyE[eRec][ba]) {
      // Here I only fill the helper profile to get average of requested particle variable for current event:
      qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eitsNClsEbyE) - 0.5, track.itsNCls(), weight);

      if (track.eta() < 0.) {
        qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eitsNClsNegEtaEbyE) - 0.5, track.itsNCls(), weight);
      } else if (track.eta() > 0.) { // TBI 20241214 for the time being, I do not care about the corner case eta = 0.
        qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eitsNClsPosEtaEbyE) - 0.5, track.itsNCls(), weight);
      }

      if (-0.8 < track.eta() && track.eta() < -0.4) {
        qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eEta0804EbyE) - 0.5, track.eta(), weight);
      } else if (-0.4 < track.eta() && track.eta() < 0.0) {
        qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eEta0400EbyE) - 0.5, track.eta(), weight);
      } else if (0.0 < track.eta() && track.eta() < 0.4) {
        qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eEta0004EbyE) - 0.5, track.eta(), weight);
      } else if (0.4 < track.eta() && track.eta() < 0.8) {
        qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eEta0408EbyE) - 0.5, track.eta(), weight);
      }

      if (0.0 < track.pt() && track.pt() < 0.5) {
        qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(ePt0005EbyE) - 0.5, track.pt(), weight);
      } else if (0.5 < track.pt() && track.pt() < 1.0) {
        qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(ePt0510EbyE) - 0.5, track.pt(), weight);
      } else if (1.0 < track.pt() && track.pt() < 5.0) {
        qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(ePt1050EbyE) - 0.5, track.pt(), weight);
      }

      // eMeanPhi:
      qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eMeanPhi) - 0.5, track.phi(), weight);

      // eMeanPt:
      qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eMeanPt) - 0.5, track.pt(), weight);

      // eMeanEta:
      qa.fQAParticleEventProEbyE[eRec][ba]->Fill(static_cast<float>(eMeanEta) - 0.5, track.eta(), weight);

      // ...

    } // if ...

    // ... and corresponding MC truth simulated (common to Run 3, Run 2 and Run 1)
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
    // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo

    if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "  No MC particle for this track, skip...");
        return;
      }
      auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle

      // 1D:
      if (ph.fFillParticleHistograms) {
        !ph.fParticleHistograms[ePhi][eSim][ba] ? true : ph.fParticleHistograms[ePhi][eSim][ba]->Fill(mcparticle.phi(), weight);
        !ph.fParticleHistograms[ePt][eSim][ba] ? true : ph.fParticleHistograms[ePt][eSim][ba]->Fill(mcparticle.pt(), weight);
        !ph.fParticleHistograms[eEta][eSim][ba] ? true : ph.fParticleHistograms[eEta][eSim][ba]->Fill(mcparticle.eta(), weight);
        // !ph.fParticleHistograms[eCharge][eSim][ba] ? true : ph.fParticleHistograms[eCharge][eSim][ba]->Fill( ... ); // TBI 20240511 there is no mcparticle.sign())
        !ph.fParticleHistograms[ePDG][eSim][ba] ? true : ph.fParticleHistograms[ePDG][eSim][ba]->Fill(mcparticle.pdgCode(), weight); // TBI 20240512 this one gets filles correctly, deduce from it charge signature
      }

      // 2D:
      if (ph.fFillParticleHistograms2D) {
        !ph.fParticleHistograms2D[ePhiPt][eSim][ba] ? true : ph.fParticleHistograms2D[ePhiPt][eSim][ba]->Fill(mcparticle.phi(), mcparticle.pt(), weight);
        !ph.fParticleHistograms2D[ePhiEta][eSim][ba] ? true : ph.fParticleHistograms2D[ePhiEta][eSim][ba]->Fill(mcparticle.phi(), mcparticle.eta(), weight);
      } // if(ph.fFillParticleHistograms2D) {

      // nD (THnSparse):
      if (ba == eAfter) { // yes, I feel sparse histograms only AFTER cuts for the time being
        // **) eDWPhi : here the fundamental 0-th axis never to be projected out is "phi"
        if (ph.fBookParticleSparseHistograms[eDWPhi]) {
          // Remark: It is mandatory that ordering in initialization here resembles the ordering in enum eDiffPhiWeights
          double vector[eDiffPhiWeights_N] = {mcparticle.phi(), mcparticle.pt(), mcparticle.eta(), 0., 0., 0.};
          // TBI 20250223 I do not have access to particle charge signature here => I set it to 0 temporarily.
          //              Then, I did not calculate and store centrality for "sim" => I set it to 0 temporarily.
          //              Same for vertex z, I could trivially extend ebye.fVz also for "sim" dimension => I set it to 0 temporarily here, until that's done.
          ph.fParticleSparseHistograms[eDWPhi][eSim]->Fill(vector, weight);
        }
        // **) eDWPt : here the fundamental 0-th axis never to be projected out is "pt"
        if (ph.fBookParticleSparseHistograms[eDWPt]) {
          // Remark: It is mandatory that ordering in initialization here resembles the ordering in enum eDiffPtWeights
          double vector[eDiffPtWeights_N] = {mcparticle.pt()};
          ph.fParticleSparseHistograms[eDWPt][eSim]->Fill(vector, weight);
        }
        // **) eDWEta : here the fundamental 0-th axis never to be projected out is "eta"
        if (ph.fBookParticleSparseHistograms[eDWEta]) {
          // Remark: It is mandatory that ordering in initialization here resembles the ordering in enum eDiffEtaWeights
          double vector[eDiffEtaWeights_N] = {mcparticle.eta()};
          ph.fParticleSparseHistograms[eDWEta][eSim]->Fill(vector, weight);
        }
      } // if (ba == eAfter) {

    } // if constexpr (rs == eRecAndSim || rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
  } // if constexpr (rs == eRec || rs == eRecAndSim || rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // -----------------------------------------------------------------------------

  // b) Fill only simulated (common to Run 3, Run 2 and Run 1):
  //    Remark #1: This branch is relevant when processing ONLY simulated data at generator level.
  //    Remark #2: In this branch, 'track' is always TracksSim = aod::McParticles, see https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
  if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {
    // 1D:
    if (ph.fFillParticleHistograms) {
      !ph.fParticleHistograms[ePhi][eSim][ba] ? true : ph.fParticleHistograms[ePhi][eSim][ba]->Fill(track.phi(), weight);
      !ph.fParticleHistograms[ePt][eSim][ba] ? true : ph.fParticleHistograms[ePt][eSim][ba]->Fill(track.pt(), weight);
      !ph.fParticleHistograms[eEta][eSim][ba] ? true : ph.fParticleHistograms[eEta][eSim][ba]->Fill(track.eta(), weight);
      // !ph.fParticleHistograms[eCharge][eSim][ba] ? true : ph.fParticleHistograms[eCharge][eSim][ba]->Fill( ... ); // TBI 20240511 there is no mcparticle.sign())
      !ph.fParticleHistograms[ePDG][eSim][ba] ? true : ph.fParticleHistograms[ePDG][eSim][ba]->Fill(track.pdgCode(), weight);
    } // if(ph.fFillParticleHistograms) {

    // 2D:
    if (ph.fFillParticleHistograms2D) {
      !ph.fParticleHistograms2D[ePhiPt][eSim][ba] ? true : ph.fParticleHistograms2D[ePhiPt][eSim][ba]->Fill(track.phi(), track.pt(), weight);
      !ph.fParticleHistograms2D[ePhiEta][eSim][ba] ? true : ph.fParticleHistograms2D[ePhiEta][eSim][ba]->Fill(track.phi(), track.eta(), weight);
    } // if(ph.fFillParticleHistograms2D) {
  } // if constexpr (rs == eSim || rs == eSim_Run2 || rs == eSim_Run1) {

  // -----------------------------------------------------------------------------

  // c) Fill reconstructed ... (Run 3 specific):
  if constexpr (rs == eRec || rs == eRecAndSim) {
    // TBI 20240511 check If I can use them for Run 2 and Run 1, but extending TracksRecSim_Run2 to Tracks_extra, etc.
    // Remark: Remember to use only eRec and eSim as array indices in histos, also for rs == eRecAndSim, etc. TBI 20240504 shall I introduce generic enum egRec and egSim for this sake?

    // ...

    // ... and corresponding MC truth simulated (Run 3 specific):
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "  No MC particle for this track, skip...");
        return;
      }

      // auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle

      // ...

    } // if constexpr (rs == eRecAndSim) {
  } // if constexpr (rs == eRec || rs == eRecAndSim) {

  // -----------------------------------------------------------------------------

  // d) Fill only simulated (Run 3 specific):
  //    Remark #1: I fill here only the histograms which are not already filled in group b) above, and which make sense only for Run 3 data.
  //    Remark #2: In this branch 'track' is always TracksSim = aod::McParticles, see https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
  //               See how I handled the case b) above.
  if constexpr (rs == eSim) {

    // ...

  } // if constexpr (rs == eSim) {

  // -----------------------------------------------------------------------------

  // e) Fill reconstructed ... (Run 1 and 2 specific):
  //    Remark: I fill here only the histograms which are not already filled in group a) above, and which make sense only for Run 1 and 2 data.
  if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

    // ...

    // ... and corresponding MC truth simulated (Run 1 and 2 specific):
    // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx ):
    if constexpr (rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "  No MC particle for this track, skip...");
        return;
      }

      // auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle

      // ...

    } // if constexpr (rs == eRecAndSim_Run2 || rs == eRecAndSim_Run1) {
  } // if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {

  // -----------------------------------------------------------------------------

  // f) Fill only simulated (Run 1 and 2 specific):
  //    Remark #1: I fill here only histograms which are not already filled in group b) above, and which make sense only for Run 1 and 2 data.
  //    Remark #2: In this branch In this branch 'track' is always TracksSim = aod::McParticles, see https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
  //               See how I handled the case b) above.
  if constexpr (rs == eSim_Run2 || rs == eSim_Run1) {

    // ...

  } // if constexpr (rs == eSim_Run2 || rs == eSim_Run1) {

  // -----------------------------------------------------------------------------

  // *) Test case:
  if constexpr (rs == eTest) {
    // This branch corresponds to process with minimal subscription - I implement just a few example cuts, just for testing purposes.
    // Only eRec is support in Test for the time being.
    // 1D:
    if (ph.fFillParticleHistograms) {
      !ph.fParticleHistograms[ePhi][eRec][ba] ? true : ph.fParticleHistograms[ePhi][eRec][ba]->Fill(track.phi(), weight);
      !ph.fParticleHistograms[ePt][eRec][ba] ? true : ph.fParticleHistograms[ePt][eRec][ba]->Fill(track.pt(), weight);
      !ph.fParticleHistograms[eEta][eRec][ba] ? true : ph.fParticleHistograms[eEta][eRec][ba]->Fill(track.eta(), weight);
    }
    // 2D:
    if (ph.fFillParticleHistograms2D) {
      !ph.fParticleHistograms2D[ePhiPt][eRec][ba] ? true : ph.fParticleHistograms2D[ePhiPt][eRec][ba]->Fill(track.phi(), track.pt(), weight);
      !ph.fParticleHistograms2D[ePhiEta][eRec][ba] ? true : ph.fParticleHistograms2D[ePhiEta][eRec][ba]->Fill(track.phi(), track.eta(), weight);
    }
  } // if constexpr (rs == eTest) {

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
    StartFunction(__FUNCTION__);
  }

  // a) Flush 'n' fill the generic Q-vectors:
  ResetQ();
  for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
    for (int wp = 0; wp < gMaxCorrelator + 1; wp++) // weight power
    {
      qv.fQ[h][wp] = qv.fQvector[h][wp];
    }
  }

  // b) Calculate correlations:
  for (int h = 1; h <= gMaxHarmonic; h++) // harmonic
  {
    // 2p:
    if (ebye.fSelectedTracks < 2) {
      return;
    }
    if (tc.fVerbose) {
      LOGF(info, "  calculating 2-particle correlations ....");
    }
    TComplex two = Two(h, -h);
    double twoC = two.Re(); // cos
    // double twoS = two.Im(); // sin
    double wTwo = Two(0, 0).Re(); // Weight is 'number of combinations' by default TBI
                                  // 20220809 add support for other weights
    if (wTwo > 0.0) {
      twoC /= wTwo;
    } else {
      LOGF(fatal, "In function \033[1;31m%s at line %d, wTwo = %f <=0. ebye.fSelectedTracks = %d\033[0m", __FUNCTION__, __LINE__, wTwo, ebye.fSelectedTracks);
    }

    if (nl.fCalculateCustomNestedLoops) {
      // e-b-e sanity check:
      TArrayI* harmonics = new TArrayI(2);
      harmonics->SetAt(h, 0);
      harmonics->SetAt(-h, 1);
      double nestedLoopValue = this->CalculateCustomNestedLoops(harmonics);
      if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(twoC - nestedLoopValue) > tc.fFloatingPointPrecision) {
        LOGF(fatal, "\033[1;31m%s at line %d : nestedLoopValue = %f is not the same as twoC = %f\033[0m", __FUNCTION__, __LINE__, nestedLoopValue, twoC);
      } else {
        LOGF(info, "\033[1;32m ebye check (integrated) with CustomNestedLoops is OK for isotropic 2-p, harmonic %d\033[0m", h);
      }
      delete harmonics;
      harmonics = NULL;
    } // if(nl.fCalculateCustomNestedLoops)

    // for on-the-fly and internal validation, rescale results with theoretical value:
    if (iv.fUseInternalValidation && iv.fRescaleWithTheoreticalInput && iv.fInternalValidationVnPsin[eVn] && TMath::Abs(iv.fInternalValidationVnPsin[eVn]->GetAt(h - 1)) > 0.) {
      twoC /= pow(iv.fInternalValidationVnPsin[eVn]->GetAt(h - 1), 2.);
    }

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
    // vs. occupancy:
    if (mupa.fCorrelationsPro[0][h - 1][AFO_OCCUPANCY]) {
      mupa.fCorrelationsPro[0][h - 1][AFO_OCCUPANCY]->Fill(ebye.fOccupancy, twoC, wTwo);
    }
    // vs. interaction rate:
    if (mupa.fCorrelationsPro[0][h - 1][AFO_INTERACTIONRATE]) {
      mupa.fCorrelationsPro[0][h - 1][AFO_INTERACTIONRATE]->Fill(ebye.fInteractionRate, twoC, wTwo);
    }
    // vs. current run duration:
    if (mupa.fCorrelationsPro[0][h - 1][AFO_CURRENTRUNDURATION]) {
      mupa.fCorrelationsPro[0][h - 1][AFO_CURRENTRUNDURATION]->Fill(ebye.fCurrentRunDuration, twoC, wTwo);
    }
    // vs. vertex z position:
    if (mupa.fCorrelationsPro[0][h - 1][AFO_VZ]) {
      mupa.fCorrelationsPro[0][h - 1][AFO_VZ]->Fill(ebye.fVz, twoC, wTwo);
    }

    // 4p:
    if (ebye.fSelectedTracks < 4) {
      continue;
    } // yes, continue, because I can still calculate 2-p in other harmonics!
    if (tc.fVerbose) {
      LOGF(info, "  calculating 4-particle correlations ....");
    }
    TComplex four = Four(h, h, -h, -h);
    double fourC = four.Re(); // cos
    // double fourS = four.Im(); // sin
    double wFour = Four(0, 0, 0, 0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
    if (wFour > 0.0) {
      fourC /= wFour;
    } else {
      LOGF(fatal, "In function \033[1;31m%s at line %d, wFour = %f <=0. ebye.fSelectedTracks = %d\033[0m", __FUNCTION__, __LINE__, wFour, ebye.fSelectedTracks);
      // TBI 20240110 shall I 'continue' here, instead of bailing out?
    }

    if (nl.fCalculateCustomNestedLoops) {
      // e-b-e sanity check:
      TArrayI* harmonics = new TArrayI(4);
      harmonics->SetAt(h, 0);
      harmonics->SetAt(h, 1);
      harmonics->SetAt(-h, 2);
      harmonics->SetAt(-h, 3);
      double nestedLoopValue = this->CalculateCustomNestedLoops(harmonics);
      if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(fourC - nestedLoopValue) > tc.fFloatingPointPrecision) {
        LOGF(fatal, "\033[1;31m%s at line %d : nestedLoopValue = %f is not the same as fourC = %f\033[0m", __FUNCTION__, __LINE__, nestedLoopValue, fourC);
      } else {
        LOGF(info, "\033[1;32m ebye check (integrated) with CustomNestedLoops is OK for isotropic 4-p, harmonic %d\033[0m", h);
      }
      delete harmonics;
      harmonics = NULL;
    } // if(nl.fCalculateCustomNestedLoops)

    // for on-the-fly and internal validation, rescale results with theoretical value:
    if (iv.fUseInternalValidation && iv.fRescaleWithTheoreticalInput && iv.fInternalValidationVnPsin[eVn] && TMath::Abs(iv.fInternalValidationVnPsin[eVn]->GetAt(h - 1)) > 0.) {
      fourC /= pow(iv.fInternalValidationVnPsin[eVn]->GetAt(h - 1), 4.);
    }

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
    // vs. occupancy:
    if (mupa.fCorrelationsPro[1][h - 1][AFO_OCCUPANCY]) {
      mupa.fCorrelationsPro[1][h - 1][AFO_OCCUPANCY]->Fill(ebye.fOccupancy, fourC, wFour);
    }
    // vs. interaction rate:
    if (mupa.fCorrelationsPro[1][h - 1][AFO_INTERACTIONRATE]) {
      mupa.fCorrelationsPro[1][h - 1][AFO_INTERACTIONRATE]->Fill(ebye.fInteractionRate, fourC, wFour);
    }
    // vs. current run duration:
    if (mupa.fCorrelationsPro[1][h - 1][AFO_CURRENTRUNDURATION]) {
      mupa.fCorrelationsPro[1][h - 1][AFO_CURRENTRUNDURATION]->Fill(ebye.fCurrentRunDuration, fourC, wFour);
    }
    // vs. vertex z position:
    if (mupa.fCorrelationsPro[1][h - 1][AFO_VZ]) {
      mupa.fCorrelationsPro[1][h - 1][AFO_VZ]->Fill(ebye.fVz, fourC, wFour);
    }

    // 6p:
    if (ebye.fSelectedTracks < 6) {
      continue;
    } // yes, continue, because I can still calculate 2-p and 4-p in other harmonics!
    if (tc.fVerbose) {
      LOGF(info, "  calculating 6-particle correlations ....");
    }
    TComplex six = Six(h, h, h, -h, -h, -h);
    double sixC = six.Re(); // cos
    // double sixS = six.Im(); // sin
    double wSix = Six(0, 0, 0, 0, 0, 0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
    if (wSix > 0.0) {
      sixC /= wSix;
    } else {
      LOGF(fatal, "In function \033[1;31m%s at line %d, wSix = %f <=0. ebye.fSelectedTracks = %d\033[0m", __FUNCTION__, __LINE__, wSix, ebye.fSelectedTracks);
      // TBI 20240110 shall I 'continue' here, instead of bailing out?
    }

    if (nl.fCalculateCustomNestedLoops) {
      // e-b-e sanity check:
      TArrayI* harmonics = new TArrayI(6);
      harmonics->SetAt(h, 0);
      harmonics->SetAt(h, 1);
      harmonics->SetAt(h, 2);
      harmonics->SetAt(-h, 3);
      harmonics->SetAt(-h, 4);
      harmonics->SetAt(-h, 5);
      double nestedLoopValue = this->CalculateCustomNestedLoops(harmonics);
      if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(sixC - nestedLoopValue) > tc.fFloatingPointPrecision) {
        LOGF(fatal, "\033[1;31m%s at line %d : nestedLoopValue = %f is not the same as sixC = %f\033[0m", __FUNCTION__, __LINE__, nestedLoopValue, sixC);
      } else {
        LOGF(info, "\033[1;32m ebye check (integrated) with CustomNestedLoops is OK for isotropic 6-p, harmonic %d\033[0m", h);
      }
      delete harmonics;
      harmonics = NULL;
    } // if(nl.fCalculateCustomNestedLoops)

    // for on-the-fly and internal validation, rescale results with theoretical value:
    if (iv.fUseInternalValidation && iv.fRescaleWithTheoreticalInput && iv.fInternalValidationVnPsin[eVn] && TMath::Abs(iv.fInternalValidationVnPsin[eVn]->GetAt(h - 1)) > 0.) {
      sixC /= pow(iv.fInternalValidationVnPsin[eVn]->GetAt(h - 1), 6.);
    }

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
    // vs. occupancy:
    if (mupa.fCorrelationsPro[2][h - 1][AFO_OCCUPANCY]) {
      mupa.fCorrelationsPro[2][h - 1][AFO_OCCUPANCY]->Fill(ebye.fOccupancy, sixC, wSix);
    }
    // vs. interaction rate:
    if (mupa.fCorrelationsPro[2][h - 1][AFO_INTERACTIONRATE]) {
      mupa.fCorrelationsPro[2][h - 1][AFO_INTERACTIONRATE]->Fill(ebye.fInteractionRate, sixC, wSix);
    }
    // vs. current run duration:
    if (mupa.fCorrelationsPro[2][h - 1][AFO_CURRENTRUNDURATION]) {
      mupa.fCorrelationsPro[2][h - 1][AFO_CURRENTRUNDURATION]->Fill(ebye.fCurrentRunDuration, sixC, wSix);
    }
    // vs. vertex z position:
    if (mupa.fCorrelationsPro[2][h - 1][AFO_VZ]) {
      mupa.fCorrelationsPro[2][h - 1][AFO_VZ]->Fill(ebye.fVz, sixC, wSix);
    }

    // 8p:
    if (ebye.fSelectedTracks < 8) {
      continue;
    } // yes, continue, because I can still calculate 2-p, 4-p and 6-p in other harmonics!
    if (tc.fVerbose) {
      LOGF(info, "  calculating 8-particle correlations ....");
    }
    TComplex eight = Eight(h, h, h, h, -h, -h, -h, -h);
    double eightC = eight.Re(); // cos
    // double eightS = eight.Im(); // sin
    double wEight = Eight(0, 0, 0, 0, 0, 0, 0, 0).Re(); // Weight is 'number of combinations' by default TBI_20210515 add support for other weights
    if (wEight > 0.0) {
      eightC /= wEight;
    } else {
      LOGF(fatal, "In function \033[1;31m%s at line %d, wEight = %f <=0. ebye.fSelectedTracks = %d\033[0m", __FUNCTION__, __LINE__, wEight, ebye.fSelectedTracks);
      // TBI 20240110 shall I 'continue' here, instead of bailing out?
    }

    if (nl.fCalculateCustomNestedLoops) {
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
      double nestedLoopValue = this->CalculateCustomNestedLoops(harmonics);
      if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(eightC - nestedLoopValue) > tc.fFloatingPointPrecision) {
        LOGF(fatal, "\033[1;31m%s at line %d : nestedLoopValue = %f is not the same as eightC = %f\033[0m", __FUNCTION__, __LINE__, nestedLoopValue, eightC);
      } else {
        LOGF(info, "\033[1;32m ebye check (integrated) with CustomNestedLoops is OK for isotropic 8-p, harmonic %d\033[0m", h);
      }
      delete harmonics;
      harmonics = NULL;
    } // if(nl.fCalculateCustomNestedLoops)

    // for on-the-fly and internal validation, rescale results with theoretical value:
    if (iv.fUseInternalValidation && iv.fRescaleWithTheoreticalInput && iv.fInternalValidationVnPsin[eVn] && TMath::Abs(iv.fInternalValidationVnPsin[eVn]->GetAt(h - 1)) > 0.) {
      eightC /= pow(iv.fInternalValidationVnPsin[eVn]->GetAt(h - 1), 8.);
    }

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
    // vs. occupancy:
    if (mupa.fCorrelationsPro[3][h - 1][AFO_OCCUPANCY]) {
      mupa.fCorrelationsPro[3][h - 1][AFO_OCCUPANCY]->Fill(ebye.fOccupancy, eightC, wEight);
    }
    // vs. interaction rate:
    if (mupa.fCorrelationsPro[3][h - 1][AFO_INTERACTIONRATE]) {
      mupa.fCorrelationsPro[3][h - 1][AFO_INTERACTIONRATE]->Fill(ebye.fInteractionRate, eightC, wEight);
    }
    // vs. current run duration:
    if (mupa.fCorrelationsPro[3][h - 1][AFO_CURRENTRUNDURATION]) {
      mupa.fCorrelationsPro[3][h - 1][AFO_CURRENTRUNDURATION]->Fill(ebye.fCurrentRunDuration, eightC, wEight);
    }
    // vs. vertex z position:
    if (mupa.fCorrelationsPro[3][h - 1][AFO_VZ]) {
      mupa.fCorrelationsPro[3][h - 1][AFO_VZ]->Fill(ebye.fVz, eightC, wEight);
    }
  } // for(int h=1;h<=gMaxHarmonic;h++) // harmonic

  // c) Flush the generic Q-vectors:
  ResetQ();

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void CalculateCorrelations()

//============================================================

void CalculateKineCorrelations(eAsFunctionOf AFO_variable)
{
  // Calculate analytically differential multiparticle correlations from Q-vectors.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) ...
  eqvectorKine qvKine = eqvectorKine_N; // which eqvectorKine enum
  // int nBins = -1; // TBI 20241111 temporarily commented out just to suppress warnings

  switch (AFO_variable) {
    case AFO_PT:
      qvKine = PTq;
      // nBins = res.fResultsPro[AFO_PT]->GetNbinsX(); // TBI 20241111 temporarily commented out just to suppress warnings
      break;
    case AFO_ETA:
      qvKine = ETAq;
      // nBins = res.fResultsPro[AFO_ETA]->GetNbinsX(); // TBI 20241111 temporarily commented out just to suppress warnings
      break;
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : This AFO_variable = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(AFO_variable));
      break;
  } // switch(AFO_variable)

  // *) Insanity checks on above settings:
  if (qvKine == eqvectorKine_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : qvKine == eqvectorKine_N => add some more entries to the case statement \033[0m", __FUNCTION__, __LINE__);
  }

  // ...

  LOGF(warning, "\033[1;33m%s at line %d : Not implemented yet, this is just a placeholder for future implementation.\033[0m", __FUNCTION__, __LINE__);

  // ...

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void CalculateKineCorrelations(eAsFunctionOf AFO_variable)

//============================================================

void CalculateTest0()
{
  // Calculate Test0.

  // a) Flush 'n' fill the generic Q-vectors;
  // b) Calculate correlations;
  // c) Flush the generic Q-vectors.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Flush 'n' fill the generic Q-vectors:
  ResetQ();
  for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
    for (int wp = 0; wp < gMaxCorrelator + 1; wp++) // weight power
    {
      qv.fQ[h][wp] = qv.fQvector[h][wp];
    }
  }

  // b) Calculate correlations:
  double correlation = 0.; // still has to be divided with 'weight' later, to get average correlation
  double weight = 0.;
  int n[gMaxCorrelator] = {0}; // array holding harmonics

  for (int mo = 0; mo < gMaxCorrelator; mo++) {
    for (int mi = 0; mi < gMaxIndex; mi++) {
      // TBI 20210913 I do not have to loop each time all the way up to gMaxCorrelator and gMaxIndex, but nevermind now, it's not a big efficiency loss.

      // Sanitize the labels (If necessary. Locally this is irrelevant):
      if (!t0.fTest0Labels[mo][mi]) // I do not stream them.
      {
        for (int v = 0; v < eAsFunctionOf_N; v++) {
          if (t0.fTest0Pro[mo][mi][v]) {
            t0.fTest0Labels[mo][mi] = new TString(t0.fTest0Pro[mo][mi][v]->GetTitle()); // there is no memory leak here, since this is executed only once due to if(!fTest0Labels[mo][mi])
            break;                                                                      // yes, since for all v they are the same, so I just need to fetch it from one
          }
        }
      } // if(!t0_afTest0Labels[mo][mi])

      if (t0.fTest0Labels[mo][mi]) {
        // Extract harmonics from TString, FS is " ":
        for (int h = 0; h <= mo; h++) {
          TObjArray* oa = t0.fTest0Labels[mo][mi]->Tokenize(" ");
          if (!oa) {
            LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
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
            LOGF(fatal, "\033[1;31m%s at line %d : Not supported yet: t0.fTest0Labels[mo][mi]->Data() = %s\033[0m", __FUNCTION__, __LINE__, t0.fTest0Labels[mo][mi]->Data());
        } // switch(mo+1)

        // Insanity check on weight:
        if (!(weight > 0.)) {
          LOGF(fatal, "\033[1;31m%s at line %d : weight = %f => Is perhaps order of correlator bigger than the number of particles? t0.fTest0Labels[mo][mi]->Data() = %s \033[0m", __FUNCTION__, __LINE__, weight, t0.fTest0Labels[mo][mi]->Data());
        }

        // e-b-e sanity check:
        if (nl.fCalculateCustomNestedLoops) {
          TArrayI* harmonics = new TArrayI(mo + 1);
          for (int i = 0; i < mo + 1; i++) {
            harmonics->SetAt(n[i], i);
          }
          double nestedLoopValue = this->CalculateCustomNestedLoops(harmonics);
          if (!(TMath::Abs(nestedLoopValue) > 0.)) {
            LOGF(info, "  ebye check (integrated) with CustomNestedLoops was NOT calculated for %d-p Test0 corr. %s", mo + 1, t0.fTest0Labels[mo][mi]->Data());
          } else if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(correlation / weight - nestedLoopValue) > tc.fFloatingPointPrecision) {
            LOGF(fatal, "\033[1;31m%s at line %d : nestedLoopValue = %f is not the same as correlation/weight = %f, for correlator %s\033[0m", __FUNCTION__, __LINE__, nestedLoopValue, correlation / weight, t0.fTest0Labels[mo][mi]->Data());
          } else {
            LOGF(info, "\033[1;32m ebye check (integrated) with CustomNestedLoops is OK for %d-p Test0 corr. %s\033[0m", mo + 1, t0.fTest0Labels[mo][mi]->Data());
          }
          delete harmonics;
          harmonics = NULL;
        } // if(nl.fCalculateCustomNestedLoops)

        // To ease comparison, rescale with theoretical value. Now all Test0 results shall be at 1. Remember that contribution from symmetry planes is here also relevant (in general):
        if (iv.fUseInternalValidation && iv.fRescaleWithTheoreticalInput && iv.fInternalValidationVnPsin[eVn] && iv.fInternalValidationVnPsin[ePsin]) {
          TArrayI* harmonics = new TArrayI(mo + 1);
          for (int i = 0; i < mo + 1; i++) {
            harmonics->SetAt(n[i], i);
          }
          TComplex theoreticalValue = this->TheoreticalValue(harmonics, iv.fInternalValidationVnPsin[eVn], iv.fInternalValidationVnPsin[ePsin]);
          if (TMath::Abs(theoreticalValue.Re()) > 0.) {
            correlation /= theoreticalValue.Re();
          }
          // TBI 20240424 for the time being, I do not do anything with imaginary part, but I could eventually...
          delete harmonics;
          harmonics = NULL;
        } // if(fUseInternalValidation && fRescaleWithTheoreticalInput)

        // Finally, fill:
        // integrated:
        if (t0.fTest0Pro[mo][mi][AFO_INTEGRATED]) {
          t0.fTest0Pro[mo][mi][AFO_INTEGRATED]->Fill(0.5, correlation / weight, weight);
        }
        // vs. multiplicity:
        if (t0.fTest0Pro[mo][mi][AFO_MULTIPLICITY]) {
          t0.fTest0Pro[mo][mi][AFO_MULTIPLICITY]->Fill(ebye.fMultiplicity + 0.5, correlation / weight, weight);
        }
        // vs. centrality:
        if (t0.fTest0Pro[mo][mi][AFO_CENTRALITY]) {
          t0.fTest0Pro[mo][mi][AFO_CENTRALITY]->Fill(ebye.fCentrality, correlation / weight, weight);
        }
        // vs. occupancy:
        if (t0.fTest0Pro[mo][mi][AFO_OCCUPANCY]) {
          t0.fTest0Pro[mo][mi][AFO_OCCUPANCY]->Fill(ebye.fOccupancy, correlation / weight, weight);
        }
        // vs. interaction rate:
        if (t0.fTest0Pro[mo][mi][AFO_INTERACTIONRATE]) {
          t0.fTest0Pro[mo][mi][AFO_INTERACTIONRATE]->Fill(ebye.fInteractionRate, correlation / weight, weight);
        }
        // vs. current run duration:
        if (t0.fTest0Pro[mo][mi][AFO_CURRENTRUNDURATION]) {
          t0.fTest0Pro[mo][mi][AFO_CURRENTRUNDURATION]->Fill(ebye.fCurrentRunDuration, correlation / weight, weight);
        }
        // vs. vertex z position:
        if (t0.fTest0Pro[mo][mi][AFO_VZ]) {
          t0.fTest0Pro[mo][mi][AFO_VZ]->Fill(ebye.fVz, correlation / weight, weight);
        }
      } // if(t0.fTest0Labels[mo][mi])
    } // for(int mi=0;mi<gMaxIndex;mi++)
  } // for(int mo=0;mo<gMaxCorrelator;mo++)

  // c) Flush the generic Q-vectors:
  ResetQ();

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void CalculateTest0()

//============================================================

void CalculateKineTest0(eAsFunctionOf AFO_variable)
{
  // Calculate analytically kine Test0 from Q-vectors.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) ...
  eqvectorKine qvKine = eqvectorKine_N; // which eqvectorKine enum
  int nBins = -1;

  switch (AFO_variable) {
    case AFO_PT:
      qvKine = PTq;
      nBins = res.fResultsPro[AFO_PT]->GetNbinsX();
      break;
    case AFO_ETA:
      qvKine = ETAq;
      nBins = res.fResultsPro[AFO_ETA]->GetNbinsX();
      break;
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : This AFO_variable = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(AFO_variable));
      break;
  } // switch(AFO_variable)

  // *) Insanity checks on above settings:
  if (qvKine == eqvectorKine_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : qvKine == eqvectorKine_N => add some more entries to the case statement \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Uniform loop over bin for all kine variables:
  for (int b = 0; b < nBins; b++) {

    // *) Ensures that in each bin of interest, I have the same cut on number of particles, like in integrated analysis:
    if ((qv.fqVectorEntries[qvKine][b] < ec.fdEventCuts[eMultiplicity][eMin]) || (qv.fqVectorEntries[qvKine][b] > ec.fdEventCuts[eMultiplicity][eMax] || TMath::Abs(qv.fqVectorEntries[qvKine][b] - ec.fdEventCuts[eMultiplicity][eMax]) < tc.fFloatingPointPrecision)) {
      if (tc.fVerbose) {
        LOGF(info, "\033[1;31m%s eMultiplicity cut in bin = %d, for qvKine = %d\033[0m", __FUNCTION__, b, static_cast<int>(qvKine));
      }
    }

    // *) Re-initialize Q-vector to be q-vector in this bin:
    // After that, I can call all standard Q-vector functions again:
    for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
      for (int wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
        qv.fQ[h][wp] = qv.fqvector[qvKine][b][h][wp];
      }
    }

    // *) Okay, let's do the differential calculus:
    double correlation = 0.;
    double weight = 0.;
    int n[gMaxCorrelator] = {0}; // array holding harmonics

    for (int mo = 0; mo < gMaxCorrelator; mo++) {
      for (int mi = 0; mi < gMaxIndex; mi++) {
        // TBI 20240221 I do not have to loop each time all the way up to gMaxCorrelator and gMaxIndex, but nevermind now, it's not a big efficiency loss.
        if (t0.fTest0Labels[mo][mi]) {
          // Extract harmonics from TString, FS is " ":
          for (int h = 0; h <= mo; h++) {
            // cout<<Form("h = %d, t0.fTest0Labels[%d][%d] = ",h,mo,mi)<<t0.fTest0Labels[mo][mi]->Data()<<endl;
            TObjArray* oa = t0.fTest0Labels[mo][mi]->Tokenize(" ");
            if (!oa) {
              LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
            }
            n[h] = TString(oa->At(h)->GetName()).Atoi();
            delete oa; // yes, otherwise it's a memory leak
          }

          if (qv.fqVectorEntries[qvKine][b] < mo + 1) {
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
              LOGF(fatal, "\033[1;31m%s at line %d : not supported yet: %s \n\n\033[0m", __FUNCTION__, __LINE__, t0.fTest0Labels[mo][mi]->Data());
          } // switch(mo+1)

          // *) e-b-e sanity check:
          if (nl.fCalculateKineCustomNestedLoops) {
            TArrayI* harmonics = new TArrayI(mo + 1);
            for (int i = 0; i < mo + 1; i++) {
              harmonics->SetAt(n[i], i);
            }
            if (!(weight > 0.)) {
              LOGF(fatal, "\033[1;31m%s at line %d : is perhaps order of some requested correlator bigger than the number of particles? Correlator = %s \033[0m", __FUNCTION__, __LINE__, t0.fTest0Labels[mo][mi]->Data());
            }
            double nestedLoopValue = this->CalculateKineCustomNestedLoops(harmonics, AFO_variable, b);
            if (!(TMath::Abs(nestedLoopValue) > 0.)) {
              LOGF(info, "  e-b-e check with CalculateKineCustomNestedLoops was NOT calculated for %d-p Test0 corr. %s, bin = %d", mo + 1, t0.fTest0Labels[mo][mi]->Data(), b + 1);
            } else if (TMath::Abs(nestedLoopValue) > 0. && TMath::Abs(correlation / weight - nestedLoopValue) > tc.fFloatingPointPrecision) {
              LOGF(fatal, "\033[1;31m%s at line %d : correlator: %s \n correlation: %f \n custom loop: %f \033[0m", __FUNCTION__, __LINE__, t0.fTest0Labels[mo][mi]->Data(), correlation / weight, nestedLoopValue);
            } else {
              LOGF(info, "\033[1;32m ebye check (differential) with CalculateKineCustomNestedLoops is OK for %d-p Test0 corr. %s, bin = %d\033[0m", mo + 1, t0.fTest0Labels[mo][mi]->Data(), b + 1);
            }
            delete harmonics;
            harmonics = NULL;
          } // if(nl.fCalculateKineCustomNestedLoops)

          // To ease comparison, rescale with theoretical value. Now all Test0 results shall be at 1:
          if (iv.fUseInternalValidation && iv.fRescaleWithTheoreticalInput && iv.fInternalValidationVnPsin[eVn] && iv.fInternalValidationVnPsin[ePsin]) {
            TArrayI* harmonics = new TArrayI(mo + 1);
            for (int i = 0; i < mo + 1; i++) {
              harmonics->SetAt(n[i], i);
            }
            TComplex theoreticalValue = TheoreticalValue(harmonics, iv.fInternalValidationVnPsin[eVn], iv.fInternalValidationVnPsin[ePsin]);
            if (TMath::Abs(theoreticalValue.Re()) > 0.) {
              correlation /= theoreticalValue.Re();
            }
            // TBI 20240424 for the time being, I do not do anything with imaginary part, but I could eventually...
            delete harmonics;
            harmonics = NULL;
          } // if(fUseInternalValidation && fRescaleWithTheoreticalInput)

          // Insanity check for the event weight:
          if (!(weight > 0.)) {
            // If it's negative, that means that sum of particle weights is smaller than "number of particles - 1"
            // In that case, you can simply rescale all particle weights, so that each of them is > 1, basically recalculate weights.root files with such a rescaling.
            LOGF(info, "\n\033[1;33m b = %d \033[0m\n", b);
            LOGF(info, "\n\033[1;33m qvKine = %d \033[0m\n", static_cast<int>(qvKine));
            LOGF(info, "\n\033[1;33m event weight = %e \033[0m\n", weight);
            LOGF(info, "\n\033[1;33m sum of particle weights = %e \033[0m\n", One(0).Re());
            LOGF(info, "\n\033[1;33m correlation = %f \033[0m\n", correlation);
            LOGF(info, "\n\033[1;33m t0.fTest0Pro[mo][mi][AFO_variable]->GetTitle() = %s \033[0m\n", t0.fTest0Pro[mo][mi][AFO_variable]->GetTitle());
            LOGF(info, "\n\033[1;33m [mo][mi][AFO_variable] = [%d][%d][%d] \033[0m\n", mo, mi, static_cast<int>(AFO_variable));
            LOGF(info, "\n\033[1;33m ebye.fSelectedTracks = %d \033[0m\n", ebye.fSelectedTracks);
            LOGF(info, "\n\033[1;33m qv.fqVectorEntries[qvKine][b] = %d \033[0m\n", qv.fqVectorEntries[qvKine][b]);
            LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
          }

          // Finally, fill:
          if (t0.fTest0Pro[mo][mi][AFO_variable]) {
            t0.fTest0Pro[mo][mi][AFO_variable]->Fill(t0.fTest0Pro[mo][mi][AFO_variable]->GetXaxis()->GetBinCenter(b + 1), correlation / weight, weight);
          } // fill in the bin center

        } // if(fTest0Labels[mo][mi])
      } // for(int mi=0;mi<gMaxIndex;mi++)
    } // for(int mo=0;mo<gMaxCorrelator;mo++)

  } // for(int b=0;b<nBins;b++)

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // CalculateKineTest0(eAsFunctionOf AFO_variable)

//============================================================

void CalculateEtaSeparations()
{
  // Calculate correlations with pseudorapidity separations.

  // Remark: this is a port and generalization of void AliFlowAnalysisWithMultiparticleCorrelations::CalculateEtaGaps(AliFlowEventSimple *anEvent)

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Calculate 2-p correlations with eta separations from Qa (-eta, index [0]) and Qb (+eta, index [1]) vectors:
  double correlation = 0.;
  double weight = 0.;
  for (int h = 0; h < gMaxHarmonic; h++) {
    if (es.fEtaSeparationsSkipHarmonics[h]) {
      continue;
    }
    for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
      if (!(qv.fQabVector[0][h][e].Rho() > 0. && qv.fQabVector[1][h][e].Rho() > 0.)) {
        continue;
      }
      if (!(qv.fMab[0][e] > 0. && qv.fMab[1][e] > 0.)) {
        continue;
      }

      // calculate correlation and weights with particular eta separation:
      correlation = TComplex(qv.fQabVector[0][h][e] * TComplex::Conjugate(qv.fQabVector[1][h][e])).Re();
      weight = qv.fMab[0][e] * qv.fMab[1][e];

      // for on-the-fly and internal validation, rescale results with theoretical value:
      if (iv.fUseInternalValidation && iv.fRescaleWithTheoreticalInput && iv.fInternalValidationVnPsin[eVn] && TMath::Abs(iv.fInternalValidationVnPsin[eVn]->GetAt(h)) > 0.) {
        correlation /= pow(iv.fInternalValidationVnPsin[eVn]->GetAt(h), 2.);
      }

      // integrated:
      if (es.fEtaSeparationsPro[h][e][AFO_INTEGRATED]) {
        es.fEtaSeparationsPro[h][e][AFO_INTEGRATED]->Fill(0.5, correlation / weight, weight);
      }

      // vs. multiplicity:
      if (es.fEtaSeparationsPro[h][e][AFO_MULTIPLICITY]) {
        es.fEtaSeparationsPro[h][e][AFO_MULTIPLICITY]->Fill(ebye.fMultiplicity + 0.5, correlation / weight, weight);
      }

      // vs. centrality:
      if (es.fEtaSeparationsPro[h][e][AFO_CENTRALITY]) {
        es.fEtaSeparationsPro[h][e][AFO_CENTRALITY]->Fill(ebye.fCentrality, correlation / weight, weight);
      }

      // vs. occupancy:
      if (es.fEtaSeparationsPro[h][e][AFO_OCCUPANCY]) {
        es.fEtaSeparationsPro[h][e][AFO_OCCUPANCY]->Fill(ebye.fOccupancy, correlation / weight, weight);
      }

      // vs. interaction rate:
      if (es.fEtaSeparationsPro[h][e][AFO_INTERACTIONRATE]) {
        es.fEtaSeparationsPro[h][e][AFO_INTERACTIONRATE]->Fill(ebye.fInteractionRate, correlation / weight, weight);
      }

      // vs. current run duration:
      if (es.fEtaSeparationsPro[h][e][AFO_CURRENTRUNDURATION]) {
        es.fEtaSeparationsPro[h][e][AFO_CURRENTRUNDURATION]->Fill(ebye.fCurrentRunDuration, correlation / weight, weight);
      }

      // vs. vertex z position:
      if (es.fEtaSeparationsPro[h][e][AFO_VZ]) {
        es.fEtaSeparationsPro[h][e][AFO_VZ]->Fill(ebye.fVz, correlation / weight, weight);
      }

    } //  for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
  } // for (int h = 0; h < gMaxHarmonic; h++) {

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void CalculateEtaSeparations()

//============================================================

void CalculateKineEtaSeparations(eAsFunctionOf AFO_variable)
{
  // Calculate differential correlations with pseudorapidity separations.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) ...
  eqvectorKine qvKine = eqvectorKine_N; // which eqvectorKine enum
  int nBins = -1;

  switch (AFO_variable) {
    case AFO_PT:
      qvKine = PTq;
      nBins = res.fResultsPro[AFO_PT]->GetNbinsX();
      break;
    case AFO_ETA:
      LOGF(fatal, "\033[1;31m%s at line %d : It doesn't make sense (i.e. AFO_ETA cannot be used here). \033[0m", __FUNCTION__, __LINE__, static_cast<int>(AFO_variable));
      break; // obsolete, but it supresses the warning
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : This AFO_variable = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(AFO_variable));
      break;
  } // switch(AFO_variable)

  // *) Insanity checks on above settings:
  if (qvKine == eqvectorKine_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : qvKine == eqvectorKine_N => add some more entries to the case statement \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Uniform loop over bin for all kine variables:
  for (int b = 0; b < nBins; b++) {

    /* TBI 20241206 Do I need to adapt and apply this cut, also for Qa and Qb? If so, most likely I would need to apply it on sum, i.e. on entries in Qa + Qb

        // *) Ensures that in each bin of interest, I have the same cut on number of particles, like in integrated analysis:
        if ((qv.fqVectorEntries[qvKine][b] < ec.fdEventCuts[eMultiplicity][eMin]) || (qv.fqVectorEntries[qvKine][b] > ec.fdEventCuts[eMultiplicity][eMax] || TMath::Abs(qv.fqVectorEntries[qvKine][b] - ec.fdEventCuts[eMultiplicity][eMax]) < tc.fFloatingPointPrecision)) {
          if (tc.fVerbose) {
            LOGF(info, "\033[1;31m%s eMultiplicity cut in bin = %d, for qvKine = %d\033[0m", __FUNCTION__, b, static_cast<int>(qvKine));
          }
        }
    */

    // Calculate differential 2-p correlations with eta separations from Qa (-eta, index [0]) and Qb (+eta, index [1]) vectors:
    double correlation = 0.;
    double weight = 0.;
    for (int h = 0; h < gMaxHarmonic; h++) {
      if (es.fEtaSeparationsSkipHarmonics[h]) {
        continue;
      }

      for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
        if (!(qv.fqabVector[0][b][h][e].Rho() > 0. && qv.fqabVector[1][b][h][e].Rho() > 0.)) {
          continue;
        }
        if (!(qv.fmab[0][b][e] > 0. && qv.fmab[1][b][e] > 0.)) {
          continue;
        }

        // calculate correlation and weights with particular eta separation:
        correlation = TComplex(qv.fqabVector[0][b][h][e] * TComplex::Conjugate(qv.fqabVector[1][b][h][e])).Re();
        weight = qv.fmab[0][b][e] * qv.fmab[1][b][e];

        // for on-the-fly and internal validation, rescale results with theoretical value:
        if (iv.fUseInternalValidation && iv.fRescaleWithTheoreticalInput && iv.fInternalValidationVnPsin[eVn] && TMath::Abs(iv.fInternalValidationVnPsin[eVn]->GetAt(h)) > 0.) {
          correlation /= pow(iv.fInternalValidationVnPsin[eVn]->GetAt(h), 2.);
        }

        // finally, fill:
        if (es.fEtaSeparationsPro[h][e][AFO_variable]) {
          es.fEtaSeparationsPro[h][e][AFO_variable]->Fill(es.fEtaSeparationsPro[h][e][AFO_variable]->GetXaxis()->GetBinCenter(b + 1), correlation / weight, weight);
        }
      }
    }
  } // for (int b = 0; b < nBins; b++) {

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void CalculateKineEtaSeparations()

//============================================================

void FillNestedLoopsContainers(const int& particleIndex, const double& dPhi, const double& dPt, const double& dEta)
{
  // Fill into the nested loop containers the current particle.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  if (tc.fInsanityCheckForEachParticle) {
    if (particleIndex < 0) {
      LOGF(fatal, "\033[1;31m%s at line %d : particleIndex = %d\033[0m", __FUNCTION__, __LINE__, particleIndex);
    }
    if (!(TMath::Abs(nl.ftaNestedLoops[0]->GetAt(particleIndex - 1)) > 0.)) {
      LOGF(fatal, "\033[1;31m%s at line %d : there are empty elements in nl.ftaNestedLoops[0] \033[0m", __FUNCTION__, __LINE__);
      // I need this protection, to ensure that all array entries are filled. If not, most likely a particle passed all
      // selection criteria, and it wasn't added to the nested loops containers
    }
  }

  // *) Fill container for angles:
  if (nl.ftaNestedLoops[0]) {
    nl.ftaNestedLoops[0]->AddAt(dPhi, particleIndex); // remember that the 2nd argument here must start from 0
  }

  // *) Fill container for weights:
  if (nl.ftaNestedLoops[1]) {
    // TBI 20240501 there is a bit of efficiency loss here, because I access Weight() again here.
    // But it doesn't matter really, in any case I evaluate nested loops only for small M during debugging.
    // Otherwise, just promote weights to data members, and initialize them only once for a given particle.
    double wPhi = 1.;
    double wPt = 1.;
    double wEta = 1.;
    if (pw.fUseWeights[wPHI]) {
      wPhi = Weight(dPhi, wPHI);
    }
    if (pw.fUseWeights[wPT]) {
      wPt = Weight(dPt, wPT);
    }
    if (pw.fUseWeights[wETA]) {
      wEta = Weight(dEta, wETA);
    }
    nl.ftaNestedLoops[1]->AddAt(wPhi * wPt * wEta, particleIndex); // remember that the 2nd argument here must start from 0
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void FillNestedLoopsContainers(const int& particleIndex, const double& dPhi, const double& dPt, const double& dEta)

//============================================================

void CalculateNestedLoops()
{
  // Calculate correlations with nested loops.

  // a) 2-particle nested loops;
  // b) 4-particle nested loops;
  // c) 6-particle nested loops;
  // d) 8-particle nested loops.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  LOGF(info, "  ebye.fSelectedTracks = %d", ebye.fSelectedTracks);
  int nParticles = ebye.fSelectedTracks;

  /* TBI 20220823 enable the lines below eventually
  if(fUseFixedNumberOfRandomlySelectedTracks)
  {
   nParticles = 0;
   for(int i=0;i<ftaNestedLoops[0]->GetSize();i++)
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
  if (nl.fMaxNestedLoop > 0 && nl.fMaxNestedLoop < 2) {
    return;
  }
  LOGF(info, "  Calculating 2-p correlations with nested loops .... ");
  for (int i1 = 0; i1 < nParticles; i1++) {
    double dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    double dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      double dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      double dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
      for (int h = 0; h < gMaxHarmonic; h++) {
        // fill cos, 2p, integreated:
        if (nl.fNestedLoopsPro[0][h][AFO_INTEGRATED]) {
          nl.fNestedLoopsPro[0][h][AFO_INTEGRATED]->Fill(
            0.5, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
        // fill cos, 2p, vs. multiplicity:
        if (nl.fNestedLoopsPro[0][h][AFO_MULTIPLICITY]) {
          nl.fNestedLoopsPro[0][h][AFO_MULTIPLICITY]->Fill(
            ebye.fMultiplicity + 0.5, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)),
            dW1 * dW2);
        }
        // fill cos, 2p, vs. centrality:
        if (nl.fNestedLoopsPro[0][h][AFO_CENTRALITY]) {
          nl.fNestedLoopsPro[0][h][AFO_CENTRALITY]->Fill(
            ebye.fCentrality, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
        // fill cos, 2p, vs. occupancy:
        if (nl.fNestedLoopsPro[0][h][AFO_OCCUPANCY]) {
          nl.fNestedLoopsPro[0][h][AFO_OCCUPANCY]->Fill(
            ebye.fOccupancy, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
        // fill cos, 2p, vs. interaction rate:
        if (nl.fNestedLoopsPro[0][h][AFO_INTERACTIONRATE]) {
          nl.fNestedLoopsPro[0][h][AFO_INTERACTIONRATE]->Fill(
            ebye.fInteractionRate, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
        // fill cos, 2p, vs. current run duration:
        if (nl.fNestedLoopsPro[0][h][AFO_CURRENTRUNDURATION]) {
          nl.fNestedLoopsPro[0][h][AFO_CURRENTRUNDURATION]->Fill(
            ebye.fCurrentRunDuration, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }
        // fill cos, 2p, vs. vertex z position:
        if (nl.fNestedLoopsPro[0][h][AFO_VZ]) {
          nl.fNestedLoopsPro[0][h][AFO_VZ]->Fill(
            ebye.fVz, TMath::Cos((h + 1.) * (dPhi1 - dPhi2)), dW1 * dW2);
        }

      } // for(int h=1; h<=6; h++)
    } // for(int i2=0; i2<nParticles; i2++)
  } // for(int i1=0; i1<nParticles; i1++)
  LOGF(info, "  Done! ");

  // b) 4-particle nested loops:
  if (nParticles < 4) {
    return;
  }
  if (nl.fMaxNestedLoop > 0 && nl.fMaxNestedLoop < 4) {
    return;
  }
  LOGF(info, "  Calculating 4-p correlations with nested loops .... ");
  for (int i1 = 0; i1 < nParticles; i1++) {
    double dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    double dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      double dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      double dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        double dPhi3 = nl.ftaNestedLoops[0]->GetAt(i3);
        double dW3 = nl.ftaNestedLoops[1]->GetAt(i3);
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          double dPhi4 = nl.ftaNestedLoops[0]->GetAt(i4);
          double dW4 = nl.ftaNestedLoops[1]->GetAt(i4);
          for (int h = 0; h < gMaxHarmonic; h++) {
            // fill cos, 4p, integreated:
            if (nl.fNestedLoopsPro[1][h][AFO_INTEGRATED]) {
              nl.fNestedLoopsPro[1][h][AFO_INTEGRATED]->Fill(0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. M:
            if (nl.fNestedLoopsPro[1][h][AFO_MULTIPLICITY]) {
              nl.fNestedLoopsPro[1][h][AFO_MULTIPLICITY]->Fill(ebye.fMultiplicity + 0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. centrality:
            if (nl.fNestedLoopsPro[1][h][AFO_CENTRALITY]) {
              nl.fNestedLoopsPro[1][h][AFO_CENTRALITY]->Fill(ebye.fCentrality, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. occupancy:
            if (nl.fNestedLoopsPro[1][h][AFO_OCCUPANCY]) {
              nl.fNestedLoopsPro[1][h][AFO_OCCUPANCY]->Fill(ebye.fOccupancy, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. interaction rate:
            if (nl.fNestedLoopsPro[1][h][AFO_INTERACTIONRATE]) {
              nl.fNestedLoopsPro[1][h][AFO_INTERACTIONRATE]->Fill(ebye.fInteractionRate, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. current run duratione:
            if (nl.fNestedLoopsPro[1][h][AFO_CURRENTRUNDURATION]) {
              nl.fNestedLoopsPro[1][h][AFO_CURRENTRUNDURATION]->Fill(ebye.fCurrentRunDuration, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
            // fill cos, 4p, all harmonics, vs. vertex z position:
            if (nl.fNestedLoopsPro[1][h][AFO_VZ]) {
              nl.fNestedLoopsPro[1][h][AFO_VZ]->Fill(ebye.fVz, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 - dPhi3 - dPhi4)), dW1 * dW2 * dW3 * dW4);
            }
          } // for(int h=0; h<gMaxHarmonic; h++)
        } // for(int i4=0; i4<nParticles; i4++)
      } // for(int i3=0; i3<nParticles; i3++)
    } // for(int i2=0; i2<nTracks; i2++)
  } // for(int i1=0; i1<nTracks; i1++)
  LOGF(info, "  Done! ");

  // c) 6-particle nested loops:
  if (nParticles < 6) {
    return;
  }
  if (nl.fMaxNestedLoop > 0 && nl.fMaxNestedLoop < 6) {
    return;
  }
  LOGF(info, "  Calculating 6-p correlations with nested loops .... ");
  for (int i1 = 0; i1 < nParticles; i1++) {
    double dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    double dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      double dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      double dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        double dPhi3 = nl.ftaNestedLoops[0]->GetAt(i3);
        double dW3 = nl.ftaNestedLoops[1]->GetAt(i3);
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          double dPhi4 = nl.ftaNestedLoops[0]->GetAt(i4);
          double dW4 = nl.ftaNestedLoops[1]->GetAt(i4);
          for (int i5 = 0; i5 < nParticles; i5++) {
            if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4) {
              continue;
            }
            double dPhi5 = nl.ftaNestedLoops[0]->GetAt(i5);
            double dW5 = nl.ftaNestedLoops[1]->GetAt(i5);
            for (int i6 = 0; i6 < nParticles; i6++) {
              if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5) {
                continue;
              }
              double dPhi6 = nl.ftaNestedLoops[0]->GetAt(i6);
              double dW6 = nl.ftaNestedLoops[1]->GetAt(i6);
              for (int h = 0; h < gMaxHarmonic; h++) {
                // fill cos, 6p, integreated:
                if (nl.fNestedLoopsPro[2][h][AFO_INTEGRATED]) {
                  nl.fNestedLoopsPro[2][h][AFO_INTEGRATED]->Fill(0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. M:
                if (nl.fNestedLoopsPro[2][h][AFO_MULTIPLICITY]) {
                  nl.fNestedLoopsPro[2][h][AFO_MULTIPLICITY]->Fill(ebye.fMultiplicity + 0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. centrality:
                if (nl.fNestedLoopsPro[2][h][AFO_CENTRALITY]) {
                  nl.fNestedLoopsPro[2][h][AFO_CENTRALITY]->Fill(ebye.fCentrality, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. occupancy:
                if (nl.fNestedLoopsPro[2][h][AFO_OCCUPANCY]) {
                  nl.fNestedLoopsPro[2][h][AFO_OCCUPANCY]->Fill(ebye.fOccupancy, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. interaction rate:
                if (nl.fNestedLoopsPro[2][h][AFO_INTERACTIONRATE]) {
                  nl.fNestedLoopsPro[2][h][AFO_INTERACTIONRATE]->Fill(ebye.fInteractionRate, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. current run duration:
                if (nl.fNestedLoopsPro[2][h][AFO_CURRENTRUNDURATION]) {
                  nl.fNestedLoopsPro[2][h][AFO_CURRENTRUNDURATION]->Fill(ebye.fCurrentRunDuration, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
                // fill cos, 6p, all harmonics, vs. vertex z position:
                if (nl.fNestedLoopsPro[2][h][AFO_VZ]) {
                  nl.fNestedLoopsPro[2][h][AFO_VZ]->Fill(ebye.fVz, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 - dPhi4 - dPhi5 - dPhi6)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6);
                }
              } // for(int h=0; h<gMaxHarmonic; h++)
            } // if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
          } // if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
        } // for(int i4=0; i4<nParticles; i4++)
      } // for(int i3=0; i3<nParticles; i3++)
    } // for(int i2=0; i2<nTracks; i2++)
  } // for(int i1=0; i1<nTracks; i1++)
  LOGF(info, "  Done! ");

  // d) 8-particle nested loops:
  if (nParticles < 8) {
    return;
  }
  if (nl.fMaxNestedLoop > 0 && nl.fMaxNestedLoop < 8) {
    return;
  }
  LOGF(info, "  Calculating 8-p correlations with nested loops .... ");
  for (int i1 = 0; i1 < nParticles; i1++) {
    double dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    double dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
    for (int i2 = 0; i2 < nParticles; i2++) {
      if (i2 == i1) {
        continue;
      }
      double dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      double dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
      for (int i3 = 0; i3 < nParticles; i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        }
        double dPhi3 = nl.ftaNestedLoops[0]->GetAt(i3);
        double dW3 = nl.ftaNestedLoops[1]->GetAt(i3);
        for (int i4 = 0; i4 < nParticles; i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          }
          double dPhi4 = nl.ftaNestedLoops[0]->GetAt(i4);
          double dW4 = nl.ftaNestedLoops[1]->GetAt(i4);
          for (int i5 = 0; i5 < nParticles; i5++) {
            if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4) {
              continue;
            }
            double dPhi5 = nl.ftaNestedLoops[0]->GetAt(i5);
            double dW5 = nl.ftaNestedLoops[1]->GetAt(i5);
            for (int i6 = 0; i6 < nParticles; i6++) {
              if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5) {
                continue;
              }
              double dPhi6 = nl.ftaNestedLoops[0]->GetAt(i6);
              double dW6 = nl.ftaNestedLoops[1]->GetAt(i6);
              for (int i7 = 0; i7 < nParticles; i7++) {
                if (i7 == i1 || i7 == i2 || i7 == i3 || i7 == i4 || i7 == i5 || i7 == i6) {
                  continue;
                }
                double dPhi7 = nl.ftaNestedLoops[0]->GetAt(i7);
                double dW7 = nl.ftaNestedLoops[1]->GetAt(i7);
                for (int i8 = 0; i8 < nParticles; i8++) {
                  if (i8 == i1 || i8 == i2 || i8 == i3 || i8 == i4 || i8 == i5 || i8 == i6 || i8 == i7) {
                    continue;
                  }
                  double dPhi8 = nl.ftaNestedLoops[0]->GetAt(i8);
                  double dW8 = nl.ftaNestedLoops[1]->GetAt(i8);
                  for (int h = 0; h < gMaxHarmonic; h++) {
                    // fill cos, 8p, integreated:
                    if (nl.fNestedLoopsPro[3][h][AFO_INTEGRATED]) {
                      nl.fNestedLoopsPro[3][h][AFO_INTEGRATED]->Fill(0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. M:
                    if (nl.fNestedLoopsPro[3][h][AFO_MULTIPLICITY]) {
                      nl.fNestedLoopsPro[3][h][AFO_MULTIPLICITY]->Fill(ebye.fMultiplicity + 0.5, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. centrality:
                    if (nl.fNestedLoopsPro[3][h][AFO_CENTRALITY]) {
                      nl.fNestedLoopsPro[3][h][AFO_CENTRALITY]->Fill(ebye.fCentrality, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. occupancy:
                    if (nl.fNestedLoopsPro[3][h][AFO_OCCUPANCY]) {
                      nl.fNestedLoopsPro[3][h][AFO_OCCUPANCY]->Fill(ebye.fOccupancy, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. interaction rate:
                    if (nl.fNestedLoopsPro[3][h][AFO_INTERACTIONRATE]) {
                      nl.fNestedLoopsPro[3][h][AFO_INTERACTIONRATE]->Fill(ebye.fInteractionRate, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. current run duration:
                    if (nl.fNestedLoopsPro[3][h][AFO_CURRENTRUNDURATION]) {
                      nl.fNestedLoopsPro[3][h][AFO_CURRENTRUNDURATION]->Fill(ebye.fCurrentRunDuration, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                    // fill cos, 8p, all harmonics, vs. vertex z position:
                    if (nl.fNestedLoopsPro[3][h][AFO_VZ]) {
                      nl.fNestedLoopsPro[3][h][AFO_VZ]->Fill(ebye.fVz, TMath::Cos((h + 1.) * (dPhi1 + dPhi2 + dPhi3 + dPhi4 - dPhi5 - dPhi6 - dPhi7 - dPhi8)), dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8);
                    }
                  } // for(int h=0; h<gMaxHarmonic; h++)
                } // for(int i8=0; i8<nParticles; i8++)
              } // for(int i7=0; i7<nParticles; i7++)
            } // for(int i6=0; i6<nParticles; i6++)
          } // for(int i5=0; i5<nParticles; i6++)
        } // for(int i4=0; i4<nParticles; i4++)
      } // for(int i3=0; i3<nParticles; i3++)
    } // for(int i2=0; i2<nParticles; i2++)
  } // for(int i1=0; i1<nParticles; i1++)

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void CalculateNestedLoops()

//============================================================

void ComparisonNestedLoopsVsCorrelations()
{
  // Compare analytic results from Q-vectors and brute force results from nested loops.
  // Use only for small multiplicities, when nested loops are still feasible.
  // Results have to be exactly the same in each case.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  int nBinsQV = -44;
  int nBinsNL = -44;
  double valueQV = 0.;
  double valueNL = 0.;

  for (int v = 0; v < eAsFunctionOf_N; v++) { // This corresponds to the ordering of variables in enum eAsFunctionOf . Here (for the time being) I compare only int, mult, cent and occu.
    if (v == AFO_PT || v == AFO_ETA) {
      continue; // TBI 20241112 correlations vs pt and vs eta are not implemented yet
    }
    nBinsQV = mupa.fCorrelationsPro[0][0][v]->GetNbinsX();
    nBinsNL = nl.fNestedLoopsPro[0][0][v]->GetNbinsX();
    if (nBinsQV != nBinsNL) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    LOGF(info, "\033[1;32m   [%d] : %s\033[0m", v, res.fResultsProXaxisTitle[v].Data());
    for (int o = 0; o < 4; o++) {
      LOGF(info, "\033[1;32m   ==== <<%d>>-particle correlations ====\033[0m", 2 * (o + 1));
      for (int h = 0; h < gMaxHarmonic; h++) {
        for (int b = 1; b <= nBinsQV; b++) {
          if (mupa.fCorrelationsPro[o][h][v]) {
            valueQV = mupa.fCorrelationsPro[o][h][v]->GetBinContent(b);
          }
          if (nl.fNestedLoopsPro[o][h][v]) {
            valueNL = nl.fNestedLoopsPro[o][h][v]->GetBinContent(b);
          }
          if (TMath::Abs(valueQV) > 0. && TMath::Abs(valueNL) > 0.) {
            LOGF(info, "   bin=%d, h=%d, Q-vectors:    %f", b, h + 1, valueQV);
            LOGF(info, "   bin=%d, h=%d, Nested loops: %f", b, h + 1, valueNL);
            if (TMath::Abs(valueQV - valueNL) > tc.fFloatingPointPrecision) {
              LOGF(info, "\n\033[1;33m[%d][%d][%d] \033[0m\n", o, h, v);
              LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
            }
          } // if(TMath::Abs(valueQV)>0. && TMath::Abs(valueNL)>0.)
        } // for(int b=1;b<=nBinsQV;b++)
      } // for (int h = 0; h < gMaxHarmonic; h++) {
      LOGF(info, ""); // new line
    } // for(int o=0;o<4;o++)
  } // for (int v = 0; v < 3; v++)

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void ComparisonNestedLoopsVsCorrelations()

//============================================================

TComplex Q(int n, int wp)
{
  // Using the fact that Q{-n,p} = Q{n,p}^*.

  if (n >= 0) {
    return qv.fQ[n][wp];
  }
  return TComplex::Conjugate(qv.fQ[-n][wp]);

} // TComplex FlowWithMultiparticleCorrelationsTask::Q(int n, int wp)

//============================================================

TComplex One(int n1)
{
  // Generic expression <exp[i(n1*phi1)]>.

  TComplex one = Q(n1, 1);

  return one;

} // TComplex FlowWithMultiparticleCorrelationsTask::One(int n1)

//============================================================

TComplex Two(int n1, int n2)
{
  // Generic two-particle correlation <exp[i(n1*phi1+n2*phi2)]>.

  TComplex two = Q(n1, 1) * Q(n2, 1) - Q(n1 + n2, 2);

  return two;

} // TComplex FlowWithMultiparticleCorrelationsTask::Two(int n1, int n2)

//============================================================

TComplex Three(int n1, int n2, int n3)
{
  // Generic three-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3)]>.

  TComplex three = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) - Q(n1 + n2, 2) * Q(n3, 1) -
                   Q(n2, 1) * Q(n1 + n3, 2) - Q(n1, 1) * Q(n2 + n3, 2) +
                   2. * Q(n1 + n2 + n3, 3);

  return three;

} // TComplex Three(int n1, int n2, int n3)

//============================================================

TComplex Four(int n1, int n2, int n3, int n4)
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

} // TComplex Four(int n1, int n2, int n3, int n4)

//============================================================

TComplex Five(int n1, int n2, int n3, int n4, int n5)
{
  // Generic five-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5)]>.

  TComplex five = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) - Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) - Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) - Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5, 1) - Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) + Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5, 1) - Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) + Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5, 1) + 2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5, 1) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) + Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5, 1) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5, 1) + 2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5, 1) - 6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5, 1) - Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) + Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5, 2) + Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5, 2) + Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5, 2) - 2. * Q(n2 + n3 + n4, 3) * Q(n1 + n5, 2) - Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) + Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5, 2) + Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5, 2) + Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5, 2) - 2. * Q(n1 + n3 + n4, 3) * Q(n2 + n5, 2) + 2. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5, 3) - 2. * Q(n3 + n4, 2) * Q(n1 + n2 + n5, 3) - Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) + Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5, 2) + Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5, 2) + Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5, 2) - 2. * Q(n1 + n2 + n4, 3) * Q(n3 + n5, 2) + 2. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5, 3) - 2. * Q(n2 + n4, 2) * Q(n1 + n3 + n5, 3) + 2. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5, 3) - 2. * Q(n1 + n4, 2) * Q(n2 + n3 + n5, 3) - 6. * Q(n4, 1) * Q(n1 + n2 + n3 + n5, 4) - Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) + Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5, 2) + Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5, 2) + Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5, 2) - 2. * Q(n1 + n2 + n3, 3) * Q(n4 + n5, 2) + 2. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5, 3) - 2. * Q(n2 + n3, 2) * Q(n1 + n4 + n5, 3) + 2. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5, 3) - 2. * Q(n1 + n3, 2) * Q(n2 + n4 + n5, 3) - 6. * Q(n3, 1) * Q(n1 + n2 + n4 + n5, 4) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5, 3) - 2. * Q(n1 + n2, 2) * Q(n3 + n4 + n5, 3) - 6. * Q(n2, 1) * Q(n1 + n3 + n4 + n5, 4) - 6. * Q(n1, 1) * Q(n2 + n3 + n4 + n5, 4) + 24. * Q(n1 + n2 + n3 + n4 + n5, 5);

  return five;

} // TComplex Five(int n1, int n2, int n3, int n4, int n5)

//============================================================

TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6)
{
  // Generic six-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6)]>.

  TComplex six = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) - Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) - Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) - Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) - Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n6, 1) + Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n6, 1) - Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n6, 1) + Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n6, 1) + 2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5, 1) * Q(n6, 1) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n6, 1) + Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n6, 1) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5, 1) * Q(n6, 1) + 2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5, 1) * Q(n6, 1) - 6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5, 1) * Q(n6, 1) - Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n6, 1) + Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n6, 1) + Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5, 2) * Q(n6, 1) + Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5, 2) * Q(n6, 1) - 2. * Q(n2 + n3 + n4, 3) * Q(n1 + n5, 2) * Q(n6, 1) - Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n6, 1) + Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n6, 1) + Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5, 2) * Q(n6, 1) + Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5, 2) * Q(n6, 1) - 2. * Q(n1 + n3 + n4, 3) * Q(n2 + n5, 2) * Q(n6, 1) + 2. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5, 3) * Q(n6, 1) - 2. * Q(n3 + n4, 2) * Q(n1 + n2 + n5, 3) * Q(n6, 1) - Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n6, 1) + Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n6, 1) + Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5, 2) * Q(n6, 1) + Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5, 2) * Q(n6, 1) - 2. * Q(n1 + n2 + n4, 3) * Q(n3 + n5, 2) * Q(n6, 1) + 2. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5, 3) * Q(n6, 1) - 2. * Q(n2 + n4, 2) * Q(n1 + n3 + n5, 3) * Q(n6, 1) + 2. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5, 3) * Q(n6, 1) - 2. * Q(n1 + n4, 2) * Q(n2 + n3 + n5, 3) * Q(n6, 1) - 6. * Q(n4, 1) * Q(n1 + n2 + n3 + n5, 4) * Q(n6, 1) - Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n6, 1) + Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n6, 1) + Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5, 2) * Q(n6, 1) + Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5, 2) * Q(n6, 1) - 2. * Q(n1 + n2 + n3, 3) * Q(n4 + n5, 2) * Q(n6, 1) + 2. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5, 3) * Q(n6, 1) - 2. * Q(n2 + n3, 2) * Q(n1 + n4 + n5, 3) * Q(n6, 1) + 2. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5, 3) * Q(n6, 1) - 2. * Q(n1 + n3, 2) * Q(n2 + n4 + n5, 3) * Q(n6, 1) - 6. * Q(n3, 1) * Q(n1 + n2 + n4 + n5, 4) * Q(n6, 1) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5, 3) * Q(n6, 1) - 2. * Q(n1 + n2, 2) * Q(n3 + n4 + n5, 3) * Q(n6, 1) - 6. * Q(n2, 1) * Q(n1 + n3 + n4 + n5, 4) * Q(n6, 1) - 6. * Q(n1, 1) * Q(n2 + n3 + n4 + n5, 4) * Q(n6, 1) + 24. * Q(n1 + n2 + n3 + n4 + n5, 5) * Q(n6, 1) - Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n6, 2) + Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n6, 2) + Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n1 + n6, 2) + Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n1 + n6, 2) - 2. * Q(n2 + n3 + n4, 3) * Q(n5, 1) * Q(n1 + n6, 2) + Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n1 + n6, 2) - Q(n3 + n4, 2) * Q(n2 + n5, 2) * Q(n1 + n6, 2) + Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n1 + n6, 2) - Q(n2 + n4, 2) * Q(n3 + n5, 2) * Q(n1 + n6, 2) - 2. * Q(n4, 1) * Q(n2 + n3 + n5, 3) * Q(n1 + n6, 2) + Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n1 + n6, 2) - Q(n2 + n3, 2) * Q(n4 + n5, 2) * Q(n1 + n6, 2) - 2. * Q(n3, 1) * Q(n2 + n4 + n5, 3) * Q(n1 + n6, 2) - 2. * Q(n2, 1) * Q(n3 + n4 + n5, 3) * Q(n1 + n6, 2) + 6. * Q(n2 + n3 + n4 + n5, 4) * Q(n1 + n6, 2) - Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n6, 2) + Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n6, 2) + Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n2 + n6, 2) + Q(n1, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n2 + n6, 2) - 2. * Q(n1 + n3 + n4, 3) * Q(n5, 1) * Q(n2 + n6, 2) + Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n2 + n6, 2) - Q(n3 + n4, 2) * Q(n1 + n5, 2) * Q(n2 + n6, 2) + Q(n1, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n2 + n6, 2) - Q(n1 + n4, 2) * Q(n3 + n5, 2) * Q(n2 + n6, 2) - 2. * Q(n4, 1) * Q(n1 + n3 + n5, 3) * Q(n2 + n6, 2) + Q(n1, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n2 + n6, 2) - Q(n1 + n3, 2) * Q(n4 + n5, 2) * Q(n2 + n6, 2) - 2. * Q(n3, 1) * Q(n1 + n4 + n5, 3) * Q(n2 + n6, 2) - 2. * Q(n1, 1) * Q(n3 + n4 + n5, 3) * Q(n2 + n6, 2) + 6. * Q(n1 + n3 + n4 + n5, 4) * Q(n2 + n6, 2) + 2. * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n2 + n6, 3) - 2. * Q(n3 + n4, 2) * Q(n5, 1) * Q(n1 + n2 + n6, 3) - 2. * Q(n4, 1) * Q(n3 + n5, 2) * Q(n1 + n2 + n6, 3) - 2. * Q(n3, 1) * Q(n4 + n5, 2) * Q(n1 + n2 + n6, 3) + 4. * Q(n3 + n4 + n5, 3) * Q(n1 + n2 + n6, 3) - Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n5, 1) * Q(n3 + n6, 2) + Q(n1 + n2, 2) * Q(n4, 1) * Q(n5, 1) * Q(n3 + n6, 2) + Q(n2, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n3 + n6, 2) + Q(n1, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n3 + n6, 2) - 2. * Q(n1 + n2 + n4, 3) * Q(n5, 1) * Q(n3 + n6, 2) + Q(n2, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n3 + n6, 2) - Q(n2 + n4, 2) * Q(n1 + n5, 2) * Q(n3 + n6, 2) + Q(n1, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n3 + n6, 2) - Q(n1 + n4, 2) * Q(n2 + n5, 2) * Q(n3 + n6, 2) - 2. * Q(n4, 1) * Q(n1 + n2 + n5, 3) * Q(n3 + n6, 2) + Q(n1, 1) * Q(n2, 1) * Q(n4 + n5, 2) * Q(n3 + n6, 2) - Q(n1 + n2, 2) * Q(n4 + n5, 2) * Q(n3 + n6, 2) - 2. * Q(n2, 1) * Q(n1 + n4 + n5, 3) * Q(n3 + n6, 2) - 2. * Q(n1, 1) * Q(n2 + n4 + n5, 3) * Q(n3 + n6, 2) + 6. * Q(n1 + n2 + n4 + n5, 4) * Q(n3 + n6, 2) + 2. * Q(n2, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n3 + n6, 3) - 2. * Q(n2 + n4, 2) * Q(n5, 1) * Q(n1 + n3 + n6, 3) - 2. * Q(n4, 1) * Q(n2 + n5, 2) * Q(n1 + n3 + n6, 3) - 2. * Q(n2, 1) * Q(n4 + n5, 2) * Q(n1 + n3 + n6, 3) + 4. * Q(n2 + n4 + n5, 3) * Q(n1 + n3 + n6, 3) + 2. * Q(n1, 1) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n3 + n6, 3) - 2. * Q(n1 + n4, 2) * Q(n5, 1) * Q(n2 + n3 + n6, 3) - 2. * Q(n4, 1) * Q(n1 + n5, 2) * Q(n2 + n3 + n6, 3) - 2. * Q(n1, 1) * Q(n4 + n5, 2) * Q(n2 + n3 + n6, 3) + 4. * Q(n1 + n4 + n5, 3) * Q(n2 + n3 + n6, 3) - 6. * Q(n4, 1) * Q(n5, 1) * Q(n1 + n2 + n3 + n6, 4) + 6. * Q(n4 + n5, 2) * Q(n1 + n2 + n3 + n6, 4) - Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n5, 1) * Q(n4 + n6, 2) + Q(n1 + n2, 2) * Q(n3, 1) * Q(n5, 1) * Q(n4 + n6, 2) + Q(n2, 1) * Q(n1 + n3, 2) * Q(n5, 1) * Q(n4 + n6, 2) + Q(n1, 1) * Q(n2 + n3, 2) * Q(n5, 1) * Q(n4 + n6, 2) - 2. * Q(n1 + n2 + n3, 3) * Q(n5, 1) * Q(n4 + n6, 2) + Q(n2, 1) * Q(n3, 1) * Q(n1 + n5, 2) * Q(n4 + n6, 2) - Q(n2 + n3, 2) * Q(n1 + n5, 2) * Q(n4 + n6, 2) + Q(n1, 1) * Q(n3, 1) * Q(n2 + n5, 2) * Q(n4 + n6, 2) - Q(n1 + n3, 2) * Q(n2 + n5, 2) * Q(n4 + n6, 2) - 2. * Q(n3, 1) * Q(n1 + n2 + n5, 3) * Q(n4 + n6, 2) + Q(n1, 1) * Q(n2, 1) * Q(n3 + n5, 2) * Q(n4 + n6, 2) - Q(n1 + n2, 2) * Q(n3 + n5, 2) * Q(n4 + n6, 2) - 2. * Q(n2, 1) * Q(n1 + n3 + n5, 3) * Q(n4 + n6, 2) - 2. * Q(n1, 1) * Q(n2 + n3 + n5, 3) * Q(n4 + n6, 2) + 6. * Q(n1 + n2 + n3 + n5, 4) * Q(n4 + n6, 2) + 2. * Q(n2, 1) * Q(n3, 1) * Q(n5, 1) * Q(n1 + n4 + n6, 3) - 2. * Q(n2 + n3, 2) * Q(n5, 1) * Q(n1 + n4 + n6, 3) - 2. * Q(n3, 1) * Q(n2 + n5, 2) * Q(n1 + n4 + n6, 3) - 2. * Q(n2, 1) * Q(n3 + n5, 2) * Q(n1 + n4 + n6, 3) + 4. * Q(n2 + n3 + n5, 3) * Q(n1 + n4 + n6, 3) + 2. * Q(n1, 1) * Q(n3, 1) * Q(n5, 1) * Q(n2 + n4 + n6, 3) - 2. * Q(n1 + n3, 2) * Q(n5, 1) * Q(n2 + n4 + n6, 3) - 2. * Q(n3, 1) * Q(n1 + n5, 2) * Q(n2 + n4 + n6, 3) - 2. * Q(n1, 1) * Q(n3 + n5, 2) * Q(n2 + n4 + n6, 3) + 4. * Q(n1 + n3 + n5, 3) * Q(n2 + n4 + n6, 3) - 6. * Q(n3, 1) * Q(n5, 1) * Q(n1 + n2 + n4 + n6, 4) + 6. * Q(n3 + n5, 2) * Q(n1 + n2 + n4 + n6, 4) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n5, 1) * Q(n3 + n4 + n6, 3) - 2. * Q(n1 + n2, 2) * Q(n5, 1) * Q(n3 + n4 + n6, 3) - 2. * Q(n2, 1) * Q(n1 + n5, 2) * Q(n3 + n4 + n6, 3) - 2. * Q(n1, 1) * Q(n2 + n5, 2) * Q(n3 + n4 + n6, 3) + 4. * Q(n1 + n2 + n5, 3) * Q(n3 + n4 + n6, 3) - 6. * Q(n2, 1) * Q(n5, 1) * Q(n1 + n3 + n4 + n6, 4) + 6. * Q(n2 + n5, 2) * Q(n1 + n3 + n4 + n6, 4) - 6. * Q(n1, 1) * Q(n5, 1) * Q(n2 + n3 + n4 + n6, 4) + 6. * Q(n1 + n5, 2) * Q(n2 + n3 + n4 + n6, 4) + 24. * Q(n5, 1) * Q(n1 + n2 + n3 + n4 + n6, 5) - Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5 + n6, 2) + Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5 + n6, 2) + Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5 + n6, 2) + Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5 + n6, 2) - 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5 + n6, 2) + Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5 + n6, 2) - Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5 + n6, 2) + Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5 + n6, 2) - Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5 + n6, 2) - 2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5 + n6, 2) + Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5 + n6, 2) - Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5 + n6, 2) - 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5 + n6, 2) - 2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5 + n6, 2) + 6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5 + n6, 2) + 2. * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5 + n6, 3) - 2. * Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5 + n6, 3) - 2. * Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5 + n6, 3) - 2. * Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5 + n6, 3) + 4. * Q(n2 + n3 + n4, 3) * Q(n1 + n5 + n6, 3) + 2. * Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5 + n6, 3) - 2. * Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5 + n6, 3) - 2. * Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5 + n6, 3) - 2. * Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5 + n6, 3) + 4. * Q(n1 + n3 + n4, 3) * Q(n2 + n5 + n6, 3) - 6. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5 + n6, 4) + 6. * Q(n3 + n4, 2) * Q(n1 + n2 + n5 + n6, 4) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5 + n6, 3) - 2. * Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5 + n6, 3) - 2. * Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5 + n6, 3) - 2. * Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5 + n6, 3) + 4. * Q(n1 + n2 + n4, 3) * Q(n3 + n5 + n6, 3) - 6. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5 + n6, 4) + 6. * Q(n2 + n4, 2) * Q(n1 + n3 + n5 + n6, 4) - 6. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5 + n6, 4) + 6. * Q(n1 + n4, 2) * Q(n2 + n3 + n5 + n6, 4) + 24. * Q(n4, 1) * Q(n1 + n2 + n3 + n5 + n6, 5) + 2. * Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5 + n6, 3) - 2. * Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5 + n6, 3) - 2. * Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5 + n6, 3) - 2. * Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5 + n6, 3) + 4. * Q(n1 + n2 + n3, 3) * Q(n4 + n5 + n6, 3) - 6. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5 + n6, 4) + 6. * Q(n2 + n3, 2) * Q(n1 + n4 + n5 + n6, 4) - 6. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5 + n6, 4) + 6. * Q(n1 + n3, 2) * Q(n2 + n4 + n5 + n6, 4) + 24. * Q(n3, 1) * Q(n1 + n2 + n4 + n5 + n6, 5) - 6. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5 + n6, 4) + 6. * Q(n1 + n2, 2) * Q(n3 + n4 + n5 + n6, 4) + 24. * Q(n2, 1) * Q(n1 + n3 + n4 + n5 + n6, 5) + 24. * Q(n1, 1) * Q(n2 + n3 + n4 + n5 + n6, 5) - 120. * Q(n1 + n2 + n3 + n4 + n5 + n6, 6);

  return six;

} // TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6)

//============================================================

TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{
  // Generic seven-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7)]>.

  int harmonic[7] = {n1, n2, n3, n4, n5, n6, n7};

  TComplex seven = Recursion(7, harmonic);

  return seven;

} // end of TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)

//============================================================

TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
  // Generic eight-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8)]>.

  int harmonic[8] = {n1, n2, n3, n4, n5, n6, n7, n8};

  TComplex eight = Recursion(8, harmonic);

  return eight;

} // end of Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)

//============================================================

TComplex Nine(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9)
{
  // Generic nine-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9)]>.

  int harmonic[9] = {n1, n2, n3, n4, n5, n6, n7, n8, n9};

  TComplex nine = Recursion(9, harmonic);

  return nine;

} // end of TComplex Nine(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9)

//============================================================

TComplex Ten(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10)
{
  // Generic ten-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9+n10*phi10)]>.

  int harmonic[10] = {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10};

  TComplex ten = Recursion(10, harmonic);

  return ten;

} // end of TComplex Ten(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10)

//============================================================

TComplex Eleven(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10, int n11)
{
  // Generic eleven-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9+n10*phi10+n11*phi11)]>.

  int harmonic[11] = {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11};

  TComplex eleven = Recursion(11, harmonic);

  return eleven;

} // end of TComplex Eleven(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10, int n11)

//============================================================

TComplex Twelve(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10, int n11, int n12)
{
  // Generic twelve-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8+n9*phi9+n10*phi10+n11*phi11+n12*phi12)]>.

  int harmonic[12] = {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12};

  TComplex twelve = Recursion(12, harmonic);

  return twelve;

} // end of TComplex Twelve(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10, int n11, int n12)

//============================================================

TComplex Recursion(int n, int* harmonic, int mult = 1, int skip = 0)
{
  // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by
  // Kristjan Gulbrandsen (gulbrand@nbi.dk).

  int nm1 = n - 1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0)
    return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip)
    return c;

  int multp1 = mult + 1;
  int nm2 = n - 2;
  int counter1 = 0;
  int hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  int counter2 = n - 3;
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
  return c - static_cast<double>(mult) * c2;

} // TComplex Recursion(int n, int* harmonic, int mult = 1, int skip = 0)

//============================================================

void ResetQ()
{
  // Reset the components of generic Q-vectors. Use it whenever you call the
  // standard functions for correlations, for some custom Q-vectors.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
    for (int wp = 0; wp < gMaxCorrelator + 1; wp++) // weight power
    {
      qv.fQ[h][wp] = TComplex(0., 0.);
    }
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void ResetQ()

//============================================================

void SetWeightsHist(TH1D* const hist, eWeights whichWeight)
{
  // Copy histogram holding weights from an external file to the corresponding data member.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Finally:
  hist->SetDirectory(0);
  pw.fWeightsHist[whichWeight] = reinterpret_cast<TH1D*>(hist->Clone());

  if (!pw.fWeightsHist[whichWeight]) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // Cosmetics: TBI 20240216 do I really want to overwrite initial cosmetics, perhaps this shall go better into MakeWeights.C ?
  //                         Or I could move all this to GetHistogramWithWeights, where in any case I am setting e.g. histogram title, etc.
  TString sVariable[eWeights_N] = {"#varphi", "p_{t}", "#eta"}; // [phi,pt,eta]
  TString sWeights[eWeights_N] = {"w_{#varphi}", "w_{p_{t}}", "w_{#eta}"};
  pw.fWeightsHist[whichWeight]->SetStats(false);
  pw.fWeightsHist[whichWeight]->GetXaxis()->SetTitle(sVariable[whichWeight].Data());
  pw.fWeightsHist[whichWeight]->GetYaxis()->SetTitle(sWeights[whichWeight].Data());
  pw.fWeightsHist[whichWeight]->SetFillColor(eFillColor);
  pw.fWeightsHist[whichWeight]->SetLineColor(eColor);
  if (!pw.fWeightsList) {
    LOGF(fatal, "\033[1;31m%s at line %d: fWeightsList is NULL. That means that you have called SetWeightsHist(...) in init(), before this TList was booked.\033[0m", __FUNCTION__, __LINE__);
  }
  pw.fWeightsList->Add(pw.fWeightsHist[whichWeight]); // This is working at the moment, because I am fetching all weights in Preprocess(), which is called after init()
                                                      // But if eventually it will be possible to fetch run number programatically in init(), I will have to re-think this line.

  // Flag:
  pw.fUseWeights[whichWeight] = true;

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void SetWeightsHist(TH1D* const hist, eWeights whichWeight)

//============================================================

void SetDiffWeightsHist(TH1D* const hist, eDiffWeights whichDiffWeight, int bin)
{
  // Copy histogram holding differential weights from an external file to the corresponding data member.

  // Remark: Do not edit histogram title here, because that's done in GetHistogramWithWeights(), because I have "filePath" info there locally.
  //         Only if I promote "filePath" to data members, re-think the design of this function, and what goes where.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Finally:
  hist->SetDirectory(0);
  pw.fDiffWeightsHist[whichDiffWeight][bin] = reinterpret_cast<TH1D*>(hist->Clone());

  if (!pw.fDiffWeightsHist[whichDiffWeight][bin]) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // Cosmetics: TBI 20240216 do I really want to overwrite initial cosmetics, perhaps this shall go better into MakeWeights.C ?
  //                         Or I could move all this to GetHistogramWithWeights, where in any case I am setting e.g. histogram title, etc.
  TString sVariable[eDiffWeights_N] = {"#varphi", "#varphi"}; // yes, for the time being, x-axis is always phi
  TString sWeights[eDiffWeights_N] = {"(w_{#varphi})_{| p_{T}}", "(w_{#varphi})_{| #eta}"};
  pw.fDiffWeightsHist[whichDiffWeight][bin]->SetStats(false);
  pw.fDiffWeightsHist[whichDiffWeight][bin]->GetXaxis()->SetTitle(sVariable[whichDiffWeight].Data());
  pw.fDiffWeightsHist[whichDiffWeight][bin]->GetYaxis()->SetTitle(sWeights[whichDiffWeight].Data());
  pw.fDiffWeightsHist[whichDiffWeight][bin]->SetFillColor(eFillColor);
  pw.fDiffWeightsHist[whichDiffWeight][bin]->SetLineColor(eColor);
  pw.fWeightsList->Add(pw.fDiffWeightsHist[whichDiffWeight][bin]); // This is working at the moment, because I am fetching all weights in Preprocess(), which is called after init()
                                                                   // But if eventually it will be possible to fetch run number programatically in init(), I will have to re-think this line.

  // Flag:
  if (!pw.fUseDiffWeights[whichDiffWeight]) // yes, set it only once to true, for all bins
  {
    pw.fUseDiffWeights[whichDiffWeight] = true;
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // SetDiffWeightsHist(TH1D* const hist, const char *variable, int bin)

//============================================================

void SetDiffWeightsSparse(THnSparseF* const sparse, eDiffWeightCategory dwc)
{
  // Copy sparse histogram holding differential phi, pt, eta, etc., weights from an external file to the corresponding data member.

  // Remark: Do not edit sparse histogram title here, because that's done in GetHistogramWithWeights(), because I have "filePath" info there locally.
  //         Only if I promote "filePath" to data members, re-think the design of this function, and what goes where.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Finally:
  // sparse->SetDirectory(0); I cannot use this for sparse
  pw.fDiffWeightsSparse[dwc] = reinterpret_cast<THnSparseF*>(sparse->Clone());

  if (!pw.fDiffWeightsSparse[dwc]) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // Within current analysis the dimension of weight for each category won't change, therefore I store it permanently:
  pw.fDWdimension[dwc] = pw.fDiffWeightsSparse[dwc]->GetNdimensions();

  // I book here immediately vectors needed to fetch the weight from the right bin of THnSparse:
  pw.fFindBinVector[dwc] = new TArrayD(pw.fDWdimension[dwc]);

  // Finally, add to corresponding TList:
  pw.fWeightsList->Add(pw.fDiffWeightsSparse[dwc]);

  /*

TBI-today

    // Cosmetics: TBI 20240216 do I really want to overwrite initial cosmetics, perhaps this shall go better into MakeWeights.C ?
    //                         Or I could move all this to GetHistogramWithWeights, where in any case I am setting e.g. histogram title, etc.
    TString sVariable[eDiffWeights_N] = {"#varphi", "#varphi"}; // yes, for the time being, x-axis is always phi
    TString sWeights[eDiffWeights_N] = {"(w_{#varphi})_{| p_{T}}", "(w_{#varphi})_{| #eta}"};
    pw.fDiffWeightsSparse[whichDiffWeight][bin]->SetStats(false);
    pw.fDiffWeightsSparse[whichDiffWeight][bin]->GetXaxis()->SetTitle(sVariable[whichDiffWeight].Data());
    pw.fDiffWeightsSparse[whichDiffWeight][bin]->GetYaxis()->SetTitle(sWeights[whichDiffWeight].Data());
    pw.fDiffWeightsSparse[whichDiffWeight][bin]->SetFillColor(eFillColor);
    pw.fDiffWeightsSparse[whichDiffWeight][bin]->SetLineColor(eColor);
    pw.fWeightsList->Add(pw.fDiffWeightsSparse[whichDiffWeight][bin]); // This is working at the moment, because I am fetching all weights in Preprocess(), which is called after init()
                                                                     // But if eventually it will be possible to fetch run number programatically in init(), I will have to re-think this line.

    // Flag:
    if (!pw.fUseDiffWeights[whichDiffWeight]) // yes, set it only once to true, for all bins
    {
      pw.fUseDiffWeights[whichDiffWeight] = true;
    }

    if (tc.fVerbose) {
      ExitFunction(__FUNCTION__);
    }

  */

} // void SetDiffWeightsSparse(THnSparseF* const sparse)

//============================================================

void SetCentralityWeightsHist(TH1D* const hist)
{
  // Copy histogram holding weights from an external file to the corresponding data member.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Finally:
  hist->SetDirectory(0);
  cw.fCentralityWeightsHist = reinterpret_cast<TH1D*>(hist->Clone());

  if (!cw.fCentralityWeightsHist) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // Cosmetics: TBI 20240216 do I really want to overwrite initial cosmetics, perhaps this shall go better into MakeCentralityWeights.C ?
  //                         Or I could move all this to GetHistogramWithCentralityWeights, where in any case I am setting e.g. histogram title, etc.
  cw.fCentralityWeightsHist->SetStats(false);
  cw.fCentralityWeightsHist->GetXaxis()->SetTitle("Centrality percentile");
  cw.fCentralityWeightsHist->GetYaxis()->SetTitle(Form("Centrality weight (%s)", ec.fsEventCuts[eCentralityEstimator].Data()));
  cw.fCentralityWeightsHist->SetFillColor(eFillColor);
  cw.fCentralityWeightsHist->SetLineColor(eColor);
  if (!cw.fCentralityWeightsList) {
    LOGF(fatal, "\033[1;31m%s at line %d: fCentralityWeightsList is NULL. That means that you have called SetCentralityWeightsHist(...) in init(), before this TList was booked.\033[0m", __FUNCTION__, __LINE__);
  }
  cw.fCentralityWeightsList->Add(cw.fCentralityWeightsHist); // This is working at the moment, because I am fetching all centrality weights in Preprocess(), which is called after init()
                                                             // But if eventually it will be possible to fetch run number programatically in init(), I will have to re-think this line.

  // Flag:
  cw.fUseCentralityWeights = true;

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void SetCentralityWeightsHist(TH1D* const hist)

//============================================================

TH1D* GetWeightsHist(eWeights whichWeight)
{
  // The standard getter.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // ...

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  // Finally:
  return pw.fWeightsHist[whichWeight];

} // TH1D* GetWeightsHist(eWeights whichWeigh)

//============================================================

TH1D* GetHistogramWithWeights(const char* filePath, const char* runNumber, const char* variable, int bin = -1)
{
  // Get and return histogram with weights from an external file.
  // If bin > 0, differential weights for that bin are searched for.
  // If bin = -1, integrated weights are searched for, i.e. in this case "bin" variable has no effect.
  // I do it this way, so as to condense GetHistogramWithWeights(...) and GetHistogramWithDiffWeights(...) from MuPa class in
  // one routine here, so that I do not duplicate code related to CCDB access, etc.

  // TBI 20240504: Here I can keep const char* variable , i.e. no need to switch to enums, because this function is called only once, at init.
  //               Nevertheless, I could switch to enums and make it more general, i.e. I could introduce additional data members and configurables,
  //               for the names of histograms with weights. Like I did it in void GetHistogramWithCustomNUA(const char* filePath, eNUAPDF variable)

  // TBI 20241021 Strictly speaking, I do not need to pass here first 2 arguments, "filePath" and "runNumber", because they are initialized at call from data members.
  //              But since this function is called only once, it's not an important performance loss. But re-think the design here eventually.
  //              If I decide to promote filePath to data member, implement it as an array, to allow possibility that different catagories of weights are fetched from different external files.

  // a) Return value;
  // b) Basic protection for arguments;
  // c) Determine from filePath if the file in on a local machine, or in AliEn, or in CCDB;
  // d) Handle the AliEn case;
  // e) Handle the CCDB case;
  // f) Handle the local case;
  // g) The final touch on histogram with weights.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
    LOGF(info, "\033[1;33m filePath = %s\033[0m", filePath);
    LOGF(info, "\033[1;33m runNumber = %s\033[0m", runNumber);
    LOGF(info, "\033[1;33m variable = %s\033[0m", variable);
    LOGF(info, "\033[1;33m bin = %d (if bin = -1, integrated weights are searched for)\033[0m", bin);
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
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // c) Determine from filePath if the file in on a local machine, or in home
  // dir AliEn, or in CCDB:
  //    Algorithm: If filePath begins with "/alice/cern.ch/" then it's in home
  //    dir AliEn. If filePath begins with "/alice-ccdb.cern.ch/" then it's in
  //    CCDB. Therefore, files in AliEn and CCDB must be specified with abs path,
  //    for local files both abs and relative paths are just fine.
  bool bFileIsInAliEn = false;
  bool bFileIsInCCDB = false;
  if (TString(filePath).BeginsWith("/alice/cern.ch/")) {
    bFileIsInAliEn = true;
  } else {
    if (TString(filePath).BeginsWith("/alice-ccdb.cern.ch/")) {
      bFileIsInCCDB = true;
    } // else {
  } // if (TString(filePath).BeginsWith("/alice/cern.ch/")) {

  if (bFileIsInAliEn) {
    // d) Handle the AliEn case:
    TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", "");
    if (!alien) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    TFile* weightsFile = TFile::Open(Form("alien://%s", filePath), "READ");
    if (!weightsFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    weightsFile->GetObject(
      "ccdb_object", baseList); // TBI 20231008 for simplicity, hardwired name
                                // of base TList is "ccdb_object" also for
                                // AliEn case, see if I need to change this
    if (!baseList) {
      // weightsFile->ls();
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
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
      LOGF(info, "\033[1;32mAccessing in CCDB %s\033[0m", TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data());
    }

    baseList = reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));

    if (!baseList) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
    }

  } else {

    // f) Handle the local case:
    //    TBI 20231008 In principle, also for the local case in O2, I could
    //    maintain the same local structure of weights as it was in AliPhysics.
    //                 But for simplicity, in O2 I organize local weights in the
    //                 same way as in AliEn or CCDB.

    // Check if the external ROOT file exists at specified path:
    if (gSystem->AccessPathName(filePath, kFileExists)) {
      LOGF(info, "\033[1;33m if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s \033[0m", filePath);
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    TFile* weightsFile = TFile::Open(filePath, "READ");
    if (!weightsFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    weightsFile->GetObject("ccdb_object", baseList); // TBI 20231008 for simplicity, hardwired name
                                                     // of base TList is "ccdb_object" also for
                                                     // local case, see if I need to change this
    if (!baseList) {
      // weightsFile->ls();
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        // baseList->ls();
        LOGF(fatal, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current run number = %s\033[0m", __FUNCTION__, __LINE__, tc.fRunNumber.Data());
      }
    }

  } // else {

  // g) The final touch on histogram with weights:
  TString histName = "";
  if (-1 == bin) {
    // Integrated weights:
    if (!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    // fetch histogram directly from this list:
    histName = TString::Format("%s_%s", variable, tc.fTaskName.Data());
    LOGF(info, "\033[1;33m%s at line %d : fetching directly hist with name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
    hist = reinterpret_cast<TH1D*>(listWithRuns->FindObject(histName.Data()));
    // if the previous search failed, descend recursively also into the nested lists:
    if (!hist) {
      LOGF(info, "\033[1;33m%s at line %d : previous attempt failed, fetching instead recursively hist with name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
      hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, histName.Data()));
    }
    if (!hist) {
      histName = TString::Format("%s", variable); // yes, for some simple tests I can have only histogram named e.g. 'phi'
      LOGF(info, "\033[1;33m%s at line %d : last attempt, fetching instead hist with trivial name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
      hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, histName.Data()));
    }
    if (!hist) {
      listWithRuns->ls();
      LOGF(fatal, "\033[1;31m%s at line %d : couldn't fetch hist with name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
    }
    hist->SetDirectory(0);
    hist->SetTitle(Form("%s, %s", filePath, runNumber)); // I have to do it here, because only here I have "filePath" av

  } else {
    // Differential weights:
    if (!(TString(variable).EqualTo("phipt") || TString(variable).EqualTo("phieta"))) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    // fetch histogram directly from this list:
    histName = TString::Format("%s[%d]_%s", variable, bin, tc.fTaskName.Data());
    LOGF(info, "\033[1;33m%s at line %d : fetching directly hist with name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
    hist = reinterpret_cast<TH1D*>(listWithRuns->FindObject(histName.Data()));
    // if the previous search failed, descend recursively also into the nested lists:
    if (!hist) {
      LOGF(info, "\033[1;33m%s at line %d : previous attempt failed, fetching instead recursively hist with name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
      hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, Form("%s[%d]_%s", variable, bin, tc.fTaskName.Data())));
    }
    if (!hist) {
      histName = TString::Format("%s[%d]", variable, bin); // yes, for some simple tests I can have only histogram named e.g. 'phipt[0]'
      LOGF(info, "\033[1;33m%s at line %d : last attempt, fetching instead hist with trivial name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
      hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, histName.Data()));
    }
    if (!hist) {
      listWithRuns->ls();
      LOGF(fatal, "\033[1;31m%s at line %d : couldn't fetch hist with name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
    }

    // *) insanity check for differential weights => check if boundaries of current bin are the same as bin boundaries for which these weights were calculated.
    //    This way I ensure that weights correspond to same kinematic cuts and binning as in current analysis.
    //      Current example format which was set in MakeWeights.C: someString(s), min < kinematic-variable-name < max
    //      Algorithm: IFS is " " and I take (N-1)th and (N-5)th entry:
    TObjArray* oa = TString(hist->GetTitle()).Tokenize(" ");
    if (!oa) {
      LOGF(fatal, "in function \033[1;31m%s at line %d \n hist->GetTitle() = %s\033[0m", __FUNCTION__, __LINE__, hist->GetTitle());
    }
    int nEntries = oa->GetEntries();

    // I need to figure out corresponding variable from results histograms and its formatting:
    eAsFunctionOf AFO = eAsFunctionOf_N;
    const char* lVariableName = "";
    if (TString(variable).EqualTo("phipt")) {
      AFO = AFO_PT;
      lVariableName = FancyFormatting("Pt");
    } else if (TString(variable).EqualTo("phieta")) {
      AFO = AFO_ETA;
      lVariableName = FancyFormatting("Eta");
    } else {
      LOGF(fatal, "\033[1;31m%s at line %d : name = %s is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<string>(variable));
    }

    // Get min and max value for bin, stored locally:
    float min = res.fResultsPro[AFO]->GetBinLowEdge(bin + 1);
    float max = res.fResultsPro[AFO]->GetBinLowEdge(bin + 2);
    if (min > max) {
      LOGF(fatal, "\033[1;33m min = %f, max = %f, res.fResultsPro[AFO]->GetName() = %s\033[0m", min, max, res.fResultsPro[AFO]->GetName());
    }

    // Compare with min and max value stored in external weights.root file using MakeWeights.C:
    if (!(TMath::Abs(TString(oa->At(nEntries - 1)->GetName()).Atof() - max) < tc.fFloatingPointPrecision)) {
      LOGF(info, "\033[1;33m hist->GetTitle() = %s, res.fResultsPro[AFO]->GetName() = %s\033[0m", hist->GetTitle(), res.fResultsPro[AFO]->GetName());
      LOGF(fatal, "in function \033[1;31m%s at line %d : mismatch in upper bin boundaries \n from title = %f , local = %f\033[0m", __FUNCTION__, __LINE__, TString(oa->At(nEntries - 1)->GetName()).Atof(), max);
    }
    if (!(TMath::Abs(TString(oa->At(nEntries - 5)->GetName()).Atof() - min) < tc.fFloatingPointPrecision)) {
      LOGF(info, "\033[1;33m hist->GetTitle() = %s, res.fResultsPro[AFO]->GetName() = %s\033[0m", hist->GetTitle(), res.fResultsPro[AFO]->GetName());
      LOGF(fatal, "in function \033[1;31m%s at line %d : mismatch in lower bin boundaries \n from title = %f , local = %f\033[0m", __FUNCTION__, __LINE__, TString(oa->At(nEntries - 5)->GetName()).Atof(), min);
    }
    delete oa; // yes, otherwise it's a memory leak

    // *) final settings and cosmetics:
    hist->SetDirectory(0);
    hist->SetTitle(Form("%s, %.2f < %s < %.2f", filePath, min, lVariableName, max));

  } // else

  // TBI 20241021 if I need to split hist title across two lines, use this technique:
  // hist->SetTitle(Form("#splitline{#scale[0.6]{%s}}{#scale[0.4]{%s}}",hist->GetTitle(),filePath));

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return hist;

} // TH1D* GetHistogramWithWeights(const char* filePath, const char* runNumber, const char* variable, int bin = -1)

//============================================================

THnSparseF* GetSparseHistogramWithWeights(const char* filePath, const char* runNumber, const char* whichCategory, const char* whichDimensions)
{
  // Get and return sparse histogram with weights from an external file.

  // Remark 1: "whichCategory" always indicates the default x-axis (0th dimension), for instance for "differential phi weights" it's "phi"

  // Remark 2: "whichDimensions" is formatted as follows: <dim1>_<dim_2>_..., for instance "pt_cent", if weights are calculated differentially as a function of pt and centrality
  //           If empty, that is also fine, I am fetching integrated <whichCategory> weights, for instance integrated phi-weights.

  // Remark 3: The nameing convention for sparse histogram in the output file is: <whichCategory>_<whichDimensions>_multiparticle-correlations-a-b_<tc.fTaskName>
  //           a) I allow possibility that "multiparticle-correlations-a-b_" is not present in the name
  //           b) In HL, fTaskName is typically subwagon name. Therefoere, it's mandatory that for a given subwagon in HL, BOTH subwagon name and fTaskName are set to the same name
  //              TBI 20250215 If I can get within my task at run time subwagon name, I can automate this step. Check if that is possible

  // TBI 20240504: Here I can keep const char* variable , i.e. no need to switch to enums, because this function is called only once, at init.
  //               Nevertheless, I could switch to enums and make it more general, i.e. I could introduce additional data members and configurables,
  //               for the names of histograms with weights. Like I did it in void GetHistogramWithCustomNUA(const char* filePath, eNUAPDF variable)

  // TBI 20241021 Strictly speaking, I do not need to pass here first 2 arguments, "filePath" and "runNumber", because they are initialized at call from data members.
  //              But since this function is called only once, it's not an important performance loss. But re-think the design here eventually.
  //              If I decide to promote filePath to data member, implement it as an array, to allow possibility that different catagories of weights are fetched from different external files.

  // a) Return value;
  // b) Basic protection for arguments;
  // c) Determine from filePath if the file in on a local machine, or in AliEn, or in CCDB;
  // d) Handle the AliEn case;
  // e) Handle the CCDB case;
  // f) Handle the local case;
  // g) The final touch on histogram with weights.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
    LOGF(info, "\033[1;33m filePath = %s\033[0m", filePath);
    LOGF(info, "\033[1;33m runNumber = %s\033[0m", runNumber);
    LOGF(info, "\033[1;33m whichDimensions = %s\033[0m", whichDimensions);
  }

  // a) Return value:
  THnSparseF* sparseHist = NULL;
  TList* baseList = NULL;     // base top-level list in the TFile, e.g. named "ccdb_object"
  TList* listWithRuns = NULL; // nested list with run-wise TList's holding run-specific weights

  // b) Basic protection for arguments:
  //    Remark: below I do one more specific check.
  if (!(TString(whichCategory).EqualTo("phi"))) { // TBI 20250215 I could in the future extend support to differential pT weights, etc.
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }
  if (TString(whichDimensions).EqualTo("")) {
    LOGF(warning, "\033[1;33m%s at line %d : whichDimensions is empty, accessing only integrated %s weights\033[0m", __FUNCTION__, __LINE__, whichCategory);
  }

  // c) Determine from filePath if the file in on a local machine, or in home
  // dir AliEn, or in CCDB:
  //    Algorithm: If filePath begins with "/alice/cern.ch/" then it's in home
  //    dir AliEn. If filePath begins with "/alice-ccdb.cern.ch/" then it's in
  //    CCDB. Therefore, files in AliEn and CCDB must be specified with abs path,
  //    for local files both abs and relative paths are just fine.
  bool bFileIsInAliEn = false;
  bool bFileIsInCCDB = false;
  if (TString(filePath).BeginsWith("/alice/cern.ch/")) {
    bFileIsInAliEn = true;
  } else {
    if (TString(filePath).BeginsWith("/alice-ccdb.cern.ch/")) {
      bFileIsInCCDB = true;
    } // else {
  } // if (TString(filePath).BeginsWith("/alice/cern.ch/")) {

  if (bFileIsInAliEn) {
    // d) Handle the AliEn case:
    TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", "");
    if (!alien) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    TFile* weightsFile = TFile::Open(Form("alien://%s", filePath), "READ");
    if (!weightsFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    weightsFile->GetObject("ccdb_object", baseList); // TBI 20231008 for simplicity, hardwired name
                                                     // of base TList is "ccdb_object" also for
                                                     // AliEn case, see if I need to change this
    if (!baseList) {
      // weightsFile->ls();
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
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
      LOGF(info, "\033[1;32mAccessing in CCDB %s\033[0m", TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data());
    }

    baseList = reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));

    if (!baseList) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
    }

  } else {

    // f) Handle the local case:
    //    TBI 20231008 In principle, also for the local case in O2, I could
    //    maintain the same local structure of weights as it was in AliPhysics.
    //                 But for simplicity, in O2 I organize local weights in the
    //                 same way as in AliEn or CCDB.

    // Check if the external ROOT file exists at specified path:
    if (gSystem->AccessPathName(filePath, kFileExists)) {
      LOGF(info, "\033[1;33m if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s \033[0m", filePath);
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    TFile* weightsFile = TFile::Open(filePath, "READ");
    if (!weightsFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    weightsFile->GetObject("ccdb_object", baseList); // TBI 20231008 for simplicity, hardwired name
                                                     // of base TList is "ccdb_object" also for
                                                     // local case, see if I need to change this
    if (!baseList) {
      // weightsFile->ls();
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        // baseList->ls();
        LOGF(fatal, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current run number = %s\033[0m", __FUNCTION__, __LINE__, tc.fRunNumber.Data());
      }
    }

  } // else {

  // g) The final touch on sparse histogram with weights:
  TString sparseHistName = "";
  if (TString(whichDimensions).EqualTo("")) {
    sparseHistName = TString::Format("%s_multiparticle-correlations-a-b", whichCategory);
  } else if (TString(whichDimensions).BeginsWith("_")) { // TBI 20250215 alternativelly, I can remove leading "_" before calling this function
    sparseHistName = TString::Format("%s%s_multiparticle-correlations-a-b", whichCategory, whichDimensions);
  } else {
    sparseHistName = TString::Format("%s_%s_multiparticle-correlations-a-b", whichCategory, whichDimensions);
  }
  // *) If not empty, I still need to appent TaskName (i.e. the cut name):
  if (!TString(tc.fTaskName).EqualTo("")) {
    sparseHistName += tc.fTaskName.Data();
  }

  // 1. fetch histogram directly from this list: const char* whichCategory, const char* whichDimensions
  LOGF(info, "\033[1;33m%s at line %d : fetching directly from the list sparse histogram with name = %s\033[0m", __FUNCTION__, __LINE__, sparseHistName.Data());
  sparseHist = reinterpret_cast<THnSparseF*>(listWithRuns->FindObject(sparseHistName.Data()));
  if (!sparseHist) {
    // try once again by chopping off "multiparticle-correlations-a-b_" from name:
    TString tmp = sparseHistName; // yes, because "ReplaceAll" below replaces in-place, and I will need sparseHistName unmodified still later
    sparseHist = reinterpret_cast<THnSparseF*>(listWithRuns->FindObject(tmp.ReplaceAll("multiparticle-correlations-a-b_", "")));
  }

  // 2. if the previous search failed, descend recursively into the nested lists:
  if (!sparseHist) {
    LOGF(info, "\033[1;33m%s at line %d : previous attempt failed, fetching instead recursively sparse histogram with name = %s\033[0m", __FUNCTION__, __LINE__, sparseHistName.Data());
    sparseHist = reinterpret_cast<THnSparseF*>(GetObjectFromList(listWithRuns, sparseHistName.Data()));
    if (!sparseHist) {
      // try once again by chopping off "multiparticle-correlations-a-b_" from name:
      TString tmp = sparseHistName; // yes, because "ReplaceAll" below replaces in-place, and I will need sparseHistName unmodified still later
      sparseHist = reinterpret_cast<THnSparseF*>(GetObjectFromList(listWithRuns, tmp.ReplaceAll("multiparticle-correlations-a-b_", "")));
    }
  }

  if (!sparseHist) {
    listWithRuns->ls();
    LOGF(fatal, "\033[1;31m%s at line %d : couldn't fetch sparse histogram with name = %s from this list\033[0m", __FUNCTION__, __LINE__, sparseHistName.Data());
  }

  sparseHist->SetTitle(Form("%s, %s", filePath, runNumber)); // I have to do it here, because only here I have "filePath" available

  //    hist->SetTitle(Form("%s, %.2f < %s < %.2f", filePath, min, lVariableName, max));

  /*
      // *) insanity check for differential weights => check if boundaries of current bin are the same as bin boundaries for which these weights were calculated.
      //    This way I ensure that weights correspond to same kinematic cuts and binning as in current analysis.
      //      Current example format which was set in MakeWeights.C: someString(s), min < kinematic-variable-name < max
      //      Algorithm: IFS is " " and I take (N-1)th and (N-5)th entry:
      TObjArray* oa = TString(hist->GetTitle()).Tokenize(" ");
      if (!oa) {
        LOGF(fatal, "in function \033[1;31m%s at line %d \n hist->GetTitle() = %s\033[0m", __FUNCTION__, __LINE__, hist->GetTitle());
      }
      int nEntries = oa->GetEntries();

      // I need to figure out corresponding variable from results histograms and its formatting:
      eAsFunctionOf AFO = eAsFunctionOf_N;
      const char* lVariableName = "";
      if (TString(variable).EqualTo("phipt")) {
        AFO = AFO_PT;
        lVariableName = FancyFormatting("Pt");
      } else if (TString(variable).EqualTo("phieta")) {
        AFO = AFO_ETA;
        lVariableName = FancyFormatting("Eta");
      } else {
        LOGF(fatal, "\033[1;31m%s at line %d : name = %s is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<string>(variable));
      }

      // Get min and max value for bin, stored locally:
      float min = res.fResultsPro[AFO]->GetBinLowEdge(bin + 1);
      float max = res.fResultsPro[AFO]->GetBinLowEdge(bin + 2);
      if (min > max) {
        LOGF(fatal, "\033[1;33m min = %f, max = %f, res.fResultsPro[AFO]->GetName() = %s\033[0m", min, max, res.fResultsPro[AFO]->GetName());
      }

      // Compare with min and max value stored in external weights.root file using MakeWeights.C:
      if (!(TMath::Abs(TString(oa->At(nEntries - 1)->GetName()).Atof() - max) < tc.fFloatingPointPrecision)) {
        LOGF(info, "\033[1;33m hist->GetTitle() = %s, res.fResultsPro[AFO]->GetName() = %s\033[0m", hist->GetTitle(), res.fResultsPro[AFO]->GetName());
        LOGF(fatal, "in function \033[1;31m%s at line %d : mismatch in upper bin boundaries \n from title = %f , local = %f\033[0m", __FUNCTION__, __LINE__, TString(oa->At(nEntries - 1)->GetName()).Atof(), max);
      }
      if (!(TMath::Abs(TString(oa->At(nEntries - 5)->GetName()).Atof() - min) < tc.fFloatingPointPrecision)) {
        LOGF(info, "\033[1;33m hist->GetTitle() = %s, res.fResultsPro[AFO]->GetName() = %s\033[0m", hist->GetTitle(), res.fResultsPro[AFO]->GetName());
        LOGF(fatal, "in function \033[1;31m%s at line %d : mismatch in lower bin boundaries \n from title = %f , local = %f\033[0m", __FUNCTION__, __LINE__, TString(oa->At(nEntries - 5)->GetName()).Atof(), min);
      }
      delete oa; // yes, otherwise it's a memory leak

      // *) final settings and cosmetics:
      hist->SetDirectory(0);

  */

  // TBI 20241021 if I need to split hist title across two lines, use this technique:
  // hist->SetTitle(Form("#splitline{#scale[0.6]{%s}}{#scale[0.4]{%s}}",hist->GetTitle(),filePath));

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return sparseHist;

} // THnSparseF* GetSparseHistogramWithWeights(const char* filePath, const char* runNumber, const char* whichCategory, const char* whichDimensions)

//============================================================

TH1D* GetHistogramWithCentralityWeights(const char* filePath, const char* runNumber)
{
  // Get and return histogram with centrality weights from an external file.

  // TBI 20241118 Shall I merge this function with GetHistogramWithWeights(...) as there is a bit of code bloat?

  // TBI 20241021 Strictly speaking, I do not need to pass here 2 arguments, "filePath" and "runNumber", because they are initialized at call from data members.
  //              But since this function is called only once, it's not an important performance loss. But re-think the design here eventually.

  // a) Return value;
  // b) Basic protection for arguments;
  // c) Determine from filePath if the file in on a local machine, or in AliEn, or in CCDB;
  // d) Handle the AliEn case;
  // e) Handle the CCDB case;
  // f) Handle the local case;
  // g) The final touch on histogram with centrality weights.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
    LOGF(info, "\033[1;33m filePath = %s\033[0m", filePath);
    LOGF(info, "\033[1;33m runNumber = %s\033[0m", runNumber);
    LOGF(info, "\033[1;33m fTaskName = %s\033[0m", tc.fTaskName.Data());
  }

  // a) Return value:
  TH1D* hist = NULL;
  TList* baseList = NULL;     // base top-level list in the TFile, e.g. named "ccdb_object"
  TList* listWithRuns = NULL; // nested list with run-wise TList's holding run-specific weights

  // b) Basic protection for arguments:
  // ...

  // c) Determine from filePath if the file in on a local machine, or in home dir AliEn, or in CCDB:
  //    Algorithm: If filePath begins with "/alice/cern.ch/" then it's in home
  //    dir AliEn. If filePath begins with "/alice-ccdb.cern.ch/" then it's in
  //    CCDB. Therefore, files in AliEn and CCDB must be specified with abs path,
  //    for local files both abs and relative paths are just fine.
  bool bFileIsInAliEn = false;
  bool bFileIsInCCDB = false;
  if (TString(filePath).BeginsWith("/alice/cern.ch/")) {
    bFileIsInAliEn = true;
  } else {
    if (TString(filePath).BeginsWith("/alice-ccdb.cern.ch/")) {
      bFileIsInCCDB = true;
    } // else {
  } // if (TString(filePath).BeginsWith("/alice/cern.ch/")) {

  if (bFileIsInAliEn) {
    // d) Handle the AliEn case:
    TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", "");
    if (!alien) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    TFile* centralityWeightsFile = TFile::Open(Form("alien://%s", filePath), "READ");
    if (!centralityWeightsFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    centralityWeightsFile->GetObject("ccdb_object", baseList); // TBI 20231008 for simplicity, hardwired name
                                                               // of base TList is "ccdb_object" also for
                                                               // AliEn case, see if I need to change this
    if (!baseList) {
      // centralityWeightsFile->ls();
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
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
      LOGF(info, "\033[1;32mAccessing in CCDB %s\033[0m", TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data());
    }

    baseList = reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));

    if (!baseList) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
    }

  } else {

    // f) Handle the local case:
    //    TBI 20231008 In principle, also for the local case in O2, I could
    //    maintain the same local structure of weights as it was in AliPhysics.
    //    But for simplicity, in O2 I organize local weights in the same way as in AliEn or CCDB.

    // Check if the external ROOT file exists at specified path:
    if (gSystem->AccessPathName(filePath, kFileExists)) {
      LOGF(info, "\033[1;33m if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s \033[0m", filePath);
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    TFile* centralityWeightsFile = TFile::Open(filePath, "READ");
    if (!centralityWeightsFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    /*
        // xxxxxxxxxxxx  TBI 20241124 remove this code

        hist = reinterpret_cast<TH1D*>(centralityWeightsFile->Get("FT0C_Default list name")); // TBI 20241122 temporary workaround
        if (!hist) {
          Exit();
        }
        hist->SetDirectory(0);
        hist->SetTitle(Form("%s, %s", filePath, runNumber)); // I have to do it here, because only here I have "filePath" available
        return hist;

        // xxxxxxxxxxxx
    */

    centralityWeightsFile->GetObject("ccdb_object", baseList); // TBI 20231008 for simplicity, hardwired name
                                                               // of base TList is "ccdb_object" also for
                                                               // local case, see if I need to change this

    if (!baseList) {
      // centralityWeightsFile->ls();
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
    }

  } // else {

  // g) The final touch on histogram with centrality weights:
  TString histName = "";

  // fetch histogram directly from this list:
  // Remark: histName must be formated as e.g. "FT0C_multiparticle-correlations-a-b" for default analysis, or "FT0C_multiparticle-correlations-a-b_someCut"

  // Isolate short centrality estimator name, e.g. "FT0C" or "V0M" TBI 20250122 move this to utility function, because I have the same code in FancyFormatting()
  TString tmp = ec.fsEventCuts[eCentralityEstimator]; // I have to introduce local TString tmp, because ReplaceAll replaces in-place
  if (tmp.BeginsWith("CentRun2")) {
    tmp.ReplaceAll("CentRun2", ""); // "CentRun2V0M" => "V0M"
  } else if (tmp.BeginsWith("Cent")) {
    tmp.ReplaceAll("Cent", ""); // "CentFT0C" => FT0C"
  } else {
    LOGF(fatal, "\033[1;31m%s at line %d : the case tmp = %s is not supported yet\033[0m", __FUNCTION__, __LINE__, tmp.Data());
  }

  histName = TString::Format("%s_multiparticle-correlations-a-b", tmp.Data()); // I can hardwire here the name, as long as my main task name is struct MultiparticleCorrelationsAB
  if (!tc.fTaskName.EqualTo("")) {
    // for non-default analysis (e.g. in subwagon), append still "_someName", where "someName" is subwagon = taskname
    histName += "_";
    histName += tc.fTaskName;
  }

  LOGF(info, "\033[1;33m%s at line %d : fetching directly hist with name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
  hist = reinterpret_cast<TH1D*>(listWithRuns->FindObject(histName.Data()));
  // if the previous search failed, descend recursively also into the nested lists:
  if (!hist) {
    LOGF(info, "\033[1;33m%s at line %d : previous attempt failed, fetching instead recursively hist with name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
    hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, histName.Data()));
  }
  if (!hist) {
    histName = tmp; // yes, for some simple tests I can have only histogram named e.g. 'FT0C'
    LOGF(info, "\033[1;33m%s at line %d : last attempt, fetching instead hist with trivial name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
    hist = reinterpret_cast<TH1D*>(GetObjectFromList(listWithRuns, histName.Data()));
  }
  if (!hist) {
    listWithRuns->ls();
    LOGF(fatal, "\033[1;31m%s at line %d : couldn't fetch hist with name = %s\033[0m", __FUNCTION__, __LINE__, histName.Data());
  }
  hist->SetDirectory(0);
  hist->SetTitle(Form("%s, %s", filePath, runNumber)); // I have to do it here, because only here I have "filePath" av

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return hist;

} // TH1D* GetHistogramWithCentralityWeights(const char* filePath, const char* runNumber)

//============================================================

TObjArray* GetDefaultObjArrayWithLabels(const char* whichDefaultLabels)
{
  // To speed up testing, I hardwire here some labels and use them directly as they are.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Define TObjArray:
  TObjArray* arr = new TObjArray();
  arr->SetOwner();

  // Define some labels, depending on the chosen option for whichDefaultLabels:
  if (TString(whichDefaultLabels).EqualTo("trivial")) {
    const int nLabels = 1;
    TString labels[nLabels] = {"2 -2"};
    for (int l = 0; l < nLabels; l++) {
      TObjString* objstr = new TObjString(labels[l].Data());
      arr->Add(objstr);
    }
  } else if (TString(whichDefaultLabels).EqualTo("scan2p")) {
    const int nLabels = 6;
    TString labels[nLabels] = {"1 -1", "2 -2", "3 -3", "4 -4", "5 -5", "6 -6"};
    for (int l = 0; l < nLabels; l++) {
      TObjString* objstr = new TObjString(labels[l].Data());
      arr->Add(objstr);
    }
  } else if (TString(whichDefaultLabels).EqualTo("standard")) {
    const int nLabels = 7;
    TString labels[nLabels] = {"1 -1", "2 -2", "3 -3", "2 1 -1 -2", "3 1 -1 -3", "3 2 -2 -3", "3 2 1 -1 -2 -3"};
    for (int l = 0; l < nLabels; l++) {
      TObjString* objstr = new TObjString(labels[l].Data());
      arr->Add(objstr);
    }
  } else if (TString(whichDefaultLabels).EqualTo("isotropic")) {
    const int nLabels = 8;
    TString labels[nLabels] = {"1 -1", "2 -2", "3 -3", "4 -4", "1 1 -1 -1", "2 2 -2 -2", "3 3 -3 -3", "4 4 -4 -4"};
    for (int l = 0; l < nLabels; l++) {
      TObjString* objstr = new TObjString(labels[l].Data());
      arr->Add(objstr);
    }
  } else if (TString(whichDefaultLabels).EqualTo("upto8th")) {
    const int nLabels = 7; // yes, because I do not care about 1-p
    TString labels[nLabels] = {"1 -1", "1 1 -1", "1 1 -1 -1", "1 1 -1 -1 -1", "1 1 1 -1 -1 -1", "1 1 1 1 -1 -1 -1", "1 1 1 1 -1 -1 -1 -1"};
    for (int l = 0; l < nLabels; l++) {
      TObjString* objstr = new TObjString(labels[l].Data());
      arr->Add(objstr);
    }
  } else if (TString(whichDefaultLabels).EqualTo("upto10th")) {
    const int nLabels = 9; // yes, because I do not care about 1-p
    TString labels[nLabels] = {"1 -1", "1 1 -1", "1 1 -1 -1", "1 1 -1 -1 -1", "1 1 1 -1 -1 -1", "1 1 1 1 -1 -1 -1", "1 1 1 1 -1 -1 -1 -1", "1 1 1 1 -1 -1 -1 -1 -1", "1 1 1 1 1 -1 -1 -1 -1 -1"};
    for (int l = 0; l < nLabels; l++) {
      TObjString* objstr = new TObjString(labels[l].Data());
      arr->Add(objstr);
    }
  } else if (TString(whichDefaultLabels).EqualTo("upto12th")) {
    const int nLabels = 11; // yes, because I do not care about 1-p
    TString labels[nLabels] = {"1 -1", "1 1 -1", "1 1 -1 -1", "1 1 -1 -1 -1", "1 1 1 -1 -1 -1", "1 1 1 1 -1 -1 -1", "1 1 1 1 -1 -1 -1 -1", "1 1 1 1 -1 -1 -1 -1 -1", "1 1 1 1 1 -1 -1 -1 -1 -1", "1 1 1 1 1 1 -1 -1 -1 -1 -1", "1 1 1 1 1 1 -1 -1 -1 -1 -1 -1"};
    for (int l = 0; l < nLabels; l++) {
      TObjString* objstr = new TObjString(labels[l].Data());
      arr->Add(objstr);
    }
  } else {
    LOGF(fatal, "\033[1;31m%s at line %d : whichDefaultLabels = %s is not supported yet \033[0m", __FUNCTION__, __LINE__, whichDefaultLabels);
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return arr;

} // TObjArray* GetDefaultObjArrayWithLabels(const char* whichDefaultLabels)

//============================================================

TObjArray* GetObjArrayWithLabels(const char* filePath)
{
  // This function extracts from an external file TObjArray named "labels", and
  // returns it. External file can be:
  //  1) on a local computer;
  //  2) in home directory AliEn => configurable "cfFileWithLabels" must begin with "/alice/cern.ch/"
  //  3) in CCDB => configurable "cfFileWithLabels" must begin with "/alice-ccdb.cern.ch/"
  // For all CCDB wisdom, see toggle "CCDB" in page "O2"

  // a) Return value;
  // b) Determine from filePath if the file in on a local machine, or in AliEn;
  // c) Handle the AliEn case;
  // d) Handle the CCDB case;
  // e) Handle the local case.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
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
  bool bFileIsInAliEn = false;
  bool bFileIsInCCDB = false;
  if (TString(filePath).BeginsWith("/alice/cern.ch/")) {
    bFileIsInAliEn = true;
  } else {
    if (TString(filePath).BeginsWith("/alice-ccdb.cern.ch/")) {
      bFileIsInCCDB = true;
    } // else {
  } // if (TString(filePath).BeginsWith("/alice/cern.ch/")) {

  TFile* oaFile = NULL; // file holding TObjArray with all labels
  if (bFileIsInAliEn) {
    // c) Handle the AliEn case:
    TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", "");
    if (!alien) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    oaFile = TFile::Open(Form("alien://%s", filePath), "READ");
    if (!oaFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    // Fetch TObjArray from external file (keep in sync with local file case below):
    TList* lok = oaFile->GetListOfKeys();
    if (!lok) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    for (int l = 0; l < lok->GetEntries(); l++) {
      oaFile->GetObject(lok->At(l)->GetName(), oa);
      if (oa && TString(oa->ClassName()).EqualTo("TObjArray")) {
        break; // TBI 20231107 the working assumption is that in an external file there is only one TObjArray object,
               // and here I fetch it, whatever its name is. The advantage is that I do not have to do
               // any additional work for TObjArray's name. Since I do not anticipate ever having more than 1
               // TObjArray in an external file, this shall be alright. With the current implementation,
               // if there are multiple TObjArray objects in the same ROOT file, the first one will be fetched.
      }
    } // for(int l=0;l<lok->GetEntries();l++)

    if (!oa) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
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
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
  } else {

    // e) Handle the local case:
    // Check if the external ROOT file exists at specified path:
    if (gSystem->AccessPathName(filePath, kFileExists)) {
      LOGF(info, "\033[1;33m if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s \033[0m", filePath);
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    oaFile = TFile::Open(filePath, "READ");
    if (!oaFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    // Fetch TObjArray from external file (keep in sync with AliEn file case above):
    TList* lok = oaFile->GetListOfKeys();
    if (!lok) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    for (int l = 0; l < lok->GetEntries(); l++) {
      oaFile->GetObject(lok->At(l)->GetName(), oa);
      if (oa && TString(oa->ClassName()).EqualTo("TObjArray")) {
        break; // TBI 20231107 the working assumption is that in an external file there is only one TObjArray object,
               // and here I fetch it, whatever its name is. The advantage is that I do not have to do
               // any additional work for TObjArray's name. Since I do not anticipate ever having more than 1
               // TObjArray in an external file, this shall be alright. With the current implementation,
               // if there are multiple TObjArray objects in the same ROOT file, the first one will be fetched.
      }
    } // for(int l=0;l<lok->GetEntries();l++)
    if (!oa) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

  } // else {

  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m%s => Fetched TObjArray named \"%s\" from file %s\033[0m", __FUNCTION__, oa->GetName(), filePath);
    ExitFunction(__FUNCTION__);
  }

  return oa;

} // TObjArray* GetObjArrayWithLabels(const char *filePath)

//============================================================

void GetHistogramWithCustomNUA(const char* filePath, eNUAPDF variable)
{
  // Get and set immediately histogram with custom NUA for specified variable, from an external file.
  // This structure of an external file is mandatory at the moment:
  //  *) There is a TList named "ccdb_object", which holds all histograms for custom NUA;
  //  *) These histograms are TH1D objects.
  // This is a port of void FlowWithMultiparticleCorrelationsTask::SetNUAPDF(TH1D* const hist, const char* variable)

  // Remark: Unlike for weights and labels, here I am trying to do everythign in the same function, that's why here return type is void, instead of TH1D*.

  // TBI 20240501 there is a bit of code repetition below, add with analogous functions GetHistogramWithWeights(...) and GetObjArrayWithLabels(...)

  // a) Local objects;
  // b) Basic protection for arguments;
  // c) Determine from filePath if the file in on a local machine, or in AliEn, or in CCDB;
  // d) Handle the AliEn case;
  // e) Handle the CCDB case;
  // f) Handle the local case;
  // g) The final touch.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
    LOGF(info, "\033[1;32m filePath = %s\033[0m", filePath);
    LOGF(info, "\033[1;32m variable = %d\033[0m", static_cast<int>(variable));
    LOGF(info, "\033[1;32m nua.fCustomNUAPDFHistNames[variable]->Data() = %s\033[0m", nua.fCustomNUAPDFHistNames[variable]->Data());
    LOGF(info, "\033[1;32m fTaskName = %s\033[0m", tc.fTaskName.Data());
  }

  // *) Basic protection:
  if (nua.fCustomNUAPDFHistNames[variable] && nua.fCustomNUAPDFHistNames[variable]->EqualTo("")) {
    LOGF(info, "\033[1;32m empty TString, variable = %d\033[0m", static_cast<int>(variable));
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // a) Local objects:
  TH1D* hist = NULL;
  TList* baseList = NULL; // base top-level list in the TFile, always named "ccdb_object"

  // b) Basic protection for arguments:
  // ...

  // c) Determine from filePath if the file in on a local machine, or in home dir AliEn, or in CCDB:
  //    Algorithm:
  //     *) If filePath begins with "/alice/cern.ch/" then it's in home dir AliEn.
  //     *) If filePath begins with "/alice-ccdb.cern.ch/" then it's in  CCDB.
  //     *) It's a local file otherwise.
  //    Therefore, files in AliEn and CCDB must be specified with abs path, for local files both abs and relative paths are just fine.
  bool bFileIsInAliEn = false;
  bool bFileIsInCCDB = false;
  if (TString(filePath).BeginsWith("/alice/cern.ch/")) {
    bFileIsInAliEn = true;
  } else {
    if (TString(filePath).BeginsWith("/alice-ccdb.cern.ch/")) {
      bFileIsInCCDB = true;
    } // else {
  } // if (TString(filePath).BeginsWith("/alice/cern.ch/")) {

  if (bFileIsInAliEn) {
    // d) Handle the AliEn case:
    TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", "");
    if (!alien) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    TFile* customNUAFile = TFile::Open(Form("alien://%s", filePath), "READ");
    if (!customNUAFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    customNUAFile->GetObject("ccdb_object", baseList);
    if (!baseList) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

  } else if (bFileIsInCCDB) {

    // e) Handle the CCDB case: Remember that here I do not access the file, instead directly object in that file.
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
      LOGF(info, "\033[1;32mAccessing in CCDB %s\033[0m", TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data());
    }

    baseList = reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));

    if (!baseList) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

  } else {

    // f) Handle the local case:

    // Check if the external ROOT file exists at specified path:
    if (gSystem->AccessPathName(filePath, kFileExists)) {
      LOGF(info, "\033[1;33m if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s \033[0m", filePath);
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    TFile* customNUAFile = TFile::Open(filePath, "READ");
    if (!customNUAFile) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    customNUAFile->GetObject("ccdb_object", baseList);
    if (!baseList) {
      // customNUAFile->ls();
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    // hist = reinterpret_cast<TH1D*>(GetObjectFromList(baseList, nua.fCustomNUAPDFHistNames[variable]->Data()));

  } // else {

  // g) The final touch:
  hist = reinterpret_cast<TH1D*>(GetObjectFromList(baseList, nua.fCustomNUAPDFHistNames[variable]->Data()));
  if (!hist) {
    LOGF(info, "\033[1;31m hist is NULL \033[0m");
    LOGF(info, "\033[1;31m variable = %d\033[0m", static_cast<int>(variable));
    LOGF(info, "\033[1;31m nua.fCustomNUAPDFHistNames[variable]->Data() = %s\033[0m", nua.fCustomNUAPDFHistNames[variable]->Data());
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }
  hist->SetDirectory(0);
  nua.fCustomNUAPDF[variable] = reinterpret_cast<TH1D*>(hist->Clone());
  nua.fCustomNUAPDF[variable]->SetTitle(Form("%s", filePath));

  // TBI 20240501 if additional cosmetics is needed, it can be implemented here

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void GetHistogramWithCustomNUA(const char* filePath, eNUAPDF variable)

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
    StartFunction(__FUNCTION__);
  }

  // a) Initialize all counters;
  int counter[gMaxCorrelator] = {0}; // is this safe?
  for (int o = 0; o < gMaxCorrelator; o++) {
    counter[o] = 0;
  } // now it's safe :-)

  // b) Fetch TObjArray with labels from an external file:
  TObjArray* oa = NULL;
  if (t0.fUseDefaultLabels) {
    oa = GetDefaultObjArrayWithLabels(t0.fWhichDefaultLabels.Data());
  } else {
    oa = GetObjArrayWithLabels(t0.fFileWithLabels.Data());
  }
  if (!oa) {
    LOGF(info, "\033[1;33m fFileWithLabels = %s \033[0m",
         t0.fFileWithLabels.Data());
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // c) Book the placeholder fTest0LabelsPlaceholder for all labels:
  int nLabels = oa->GetEntries();
  t0.fTest0LabelsPlaceholder =
    new TH1I("fTest0LabelsPlaceholder",
             Form("placeholder for all labels, %d in total", nLabels),
             nLabels, 0, nLabels);
  t0.fTest0LabelsPlaceholder->SetStats(false);

  // d) Finally, store the labels from external source into placeholder:
  int bin = 1; // used only for fTest0LabelsPlaceholder
  int order = -44;
  for (int e = 0; e < nLabels; e++) {
    TObjArray* temp = TString(oa->At(e)->GetName()).Tokenize(" ");
    if (!temp) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
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
  } // for(int e=0; e<nLabels; e++)

  // e) Insantity check on labels:
  //    Here I am merely checking that harmonic larget than gMaxHarmonic was not requested.
  for (int b = 1; b <= t0.fTest0LabelsPlaceholder->GetXaxis()->GetNbins(); b++) {
    TObjArray* temp = TString(t0.fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b)).Tokenize(" ");
    for (int h = 0; h < temp->GetEntries(); h++) {
      if (TMath::Abs(TString(temp->At(h)->GetName()).Atoi()) > gMaxHarmonic) {
        LOGF(info, "\033[1;31m bin = %d, label = %s, gMaxHarmonic = %d\033[0m", b, t0.fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b), static_cast<int>(gMaxHarmonic));
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      } // if(TString(temp->At(h)->GetName()).Atoi() > gMaxHarmonic) {
    } // for(int h = 0; h < temp->GetEntries(); h++) {
    delete temp; // yes, otherwise it's a memory leak
  } // for(int b = 1; b <= t0.fTest0LabelsPlaceholder->GetXaxis()->GetNbins(); b++) {

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void StoreLabelsInPlaceholder()

//============================================================

bool RetrieveCorrelationsLabels()
{
  // Generate the labels of all correlations of interest, i.e. retrieve them
  // from TH1I *t0.fTest0LabelsPlaceholder

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  int counter[gMaxCorrelator] = {0}; // is this safe?
  for (int o = 0; o < gMaxCorrelator; o++) {
    counter[o] = 0;
  } // now it's safe :-)

  int nBins = t0.fTest0LabelsPlaceholder->GetXaxis()->GetNbins();

  int order = -44;
  for (int b = 1; b <= nBins; b++) {
    TObjArray* oa = TString(t0.fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b))
                      .Tokenize(" ");
    if (!oa) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    order = oa->GetEntries();
    delete oa; // yes, otherwise it's a memory leak
    if (0 == order) {
      continue;
    } // empty lines, or the label format which is not supported
    // 1-p => 0, 2-p => 1, etc.:
    t0.fTest0Labels[order - 1][counter[order - 1]] = new TString(t0.fTest0LabelsPlaceholder->GetXaxis()->GetBinLabel(b)); // okay...
    counter[order - 1]++;
  } // for(int b=1;b<=nBins;b++)

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return true;

} // bool RetrieveCorrelationsLabels()

//============================================================

TObject* GetObjectFromList(TList* list, const char* objectName) // Last update: 20210918
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

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Insanity checks:
  if (!list) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }
  if (!objectName) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
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

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return NULL;

} // TObject* GetObjectFromList(TList *list, char *objectName)

//============================================================

double Weight(const double& value, eWeights whichWeight) // value, integrated [phi,pt,eta] weight
{
  // Determine particle weight.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
    LOGF(info, "\033[1;32m value = %f\033[0m", value);
    LOGF(info, "\033[1;32m variable = %d\033[0m", static_cast<int>(whichWeight));
  }

  if (!pw.fWeightsHist[whichWeight]) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  int bin = pw.fWeightsHist[whichWeight]->FindBin(value);
  double weight = 0.;
  if (bin > pw.fWeightsHist[whichWeight]->GetNbinsX()) {
    weight = 0.; // we are in the overflow, ignore this particle TBI_20210524 is
                 // this really the correct procedure?
  } else {
    weight = pw.fWeightsHist[whichWeight]->GetBinContent(bin);
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return weight;

} // Weight(const double &value, eWeights whichWeight) // value, integrated [phi,pt,eta] weight

//============================================================

double WeightFromSparse(const double& dPhi, const double& dPt, const double& dEta, const double& dCharge, eDiffWeightCategory dwc)
{
  // Determine differential multidimensional particle weight using sparse histograms.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Reduce dimensionality is possible, i.e. look up only the dimensions in THnSparse which were requested in this analysis:
  Int_t dim = 1; // yes, because dimension 0 is always reserved for each category
  switch (dwc) {
    case eDWPhi: {
      // Remember that ordering here has to resemble ordering in eDiffPhiWeights
      pw.fFindBinVector[dwc]->AddAt(dPhi, 0); // special treatment for phi in eDWPhi category
      if (pw.fUseDiffPhiWeights[wPhiPtAxis]) {
        pw.fFindBinVector[dwc]->AddAt(dPt, dim++);
      }
      if (pw.fUseDiffPhiWeights[wPhiEtaAxis]) {
        pw.fFindBinVector[dwc]->AddAt(dEta, dim++);
      }
      if (pw.fUseDiffPhiWeights[wPhiChargeAxis]) {
        pw.fFindBinVector[dwc]->AddAt(dCharge, dim++);
      }
      if (pw.fUseDiffPhiWeights[wPhiCentralityAxis]) {
        pw.fFindBinVector[dwc]->AddAt(ebye.fCentrality, dim++);
      }
      if (pw.fUseDiffPhiWeights[wPhiVertex_zAxis]) {
        pw.fFindBinVector[dwc]->AddAt(ebye.fVz, dim++);
      }
      // ...
      break;
    }
    case eDWPt: {
      pw.fFindBinVector[dwc]->AddAt(dPt, 0); // special treatment for pt in eDWPt category
      // Remember that ordering here has to resemble ordering in eDiffPtWeights
      // if(pw.fUseDiffPtWeights[...]) {
      //   pw.fFindBinVector[dwc]->AddAt(..., dim++); // skeleton for next dimension
      // }
      // ...
      break;
    }
    case eDWEta: {
      pw.fFindBinVector[dwc]->AddAt(dEta, 0); // special treatment for eta in eDWEta category
      // Remember that ordering here has to resemble ordering in eDiffEtaWeights
      // if(pw.fUseDiffEtaWeights[...]) {
      //   pw.fFindBinVector[dwc]->AddAt(..., dim++); // skeleton for next dimension
      // }
      // ...
      break;
    }
    default: {
      LOGF(fatal, "\033[1;31m%s at line %d : This differential weight category, dwc = %d, is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(dwc));
      break;
    }
  } // switch(dwc)

  // *) Insanity check:
  // **) ...
  if (!pw.fDiffWeightsSparse[dwc]) {
    LOGF(fatal, "\033[1;31m dwc = %d\033[0m", static_cast<int>(dwc));
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // **) Check that dimensions of vector I will use to fetch the right bin content and sparse histogram with weights do match:
  if (tc.fInsanityCheckForEachParticle) {
    if (pw.fFindBinVector[dwc]->GetSize() != pw.fDiffWeightsSparse[dwc]->GetNdimensions()) {
      LOGF(fatal, "\033[1;31m dwc = %d\033[0m", static_cast<int>(dwc));
      LOGF(fatal, "\033[1;31m pw.fFindBinVector[dwc]->GetSize() = %d\033[0m", pw.fFindBinVector[dwc]->GetSize());
      LOGF(fatal, "\033[1;31m pw.fDiffWeightsSparse[dwc]->GetNdimensions() = %d\033[0m", pw.fDiffWeightsSparse[dwc]->GetNdimensions());
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
  } // if(tc.fInsanityCheckForEachParticle)

  // *) okay, let's fetch the weight:
  int bin = pw.fDiffWeightsSparse[dwc]->GetBin(pw.fFindBinVector[dwc]->GetArray()); // this is the general bin, corresponding to the actual multidimensional bin
  // TBI 20250224 do I need some insanity check here, e.g. that bin is neither in overflow nor in underflow?
  double weight = pw.fDiffWeightsSparse[dwc]->GetBinContent(bin);

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return weight;

} // double WeightFromSparse(...)

//============================================================

double DiffWeight(const double& valueY, const double& valueX, eqvectorKine variableX)
{
  // Determine differential particle weight y(x). For the time being, "y = phi" always, but this can be generalized.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Determine first to which bin the 'valueX' corresponds to.
  //    Based on that, I decide from which histogram I fetch weight for y. See MakeWeights.C

  // *) Mapping between enums "eqvectorKine" on one side, and enums "eAsFunctionOf" and "eDiffWeights" on the other:
  eAsFunctionOf AFO_var = eAsFunctionOf_N;      // this local variable determines the enum "eAsFunctionOf" which corresponds to enum "eqvectorKine"
  eDiffWeights AFO_diffWeight = eDiffWeights_N; // this local variable determines the enum "eDiffWeights" which corresponds to enum "eqvectorKine"
  if (variableX == PTq) {
    AFO_var = AFO_PT;
    AFO_diffWeight = wPHIPT;
  } else if (variableX == ETAq) {
    AFO_var = AFO_ETA;
    AFO_diffWeight = wPHIETA;
  }

  // *) Insanity checks on above settings:
  if (AFO_var == eAsFunctionOf_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : AFO_var == eAsFunctionOf_N => add some more entries to the case statement \033[0m", __FUNCTION__, __LINE__);
  }
  if (AFO_diffWeight == eDiffWeights_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : AFO_diffWeight == eDiffWeights_N => add some more entries to the case statement \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Determine first to which bin the 'valueX' corresponds to.
  //    Based on that, I decide from which histogram I fetch weight for y. See MakeWeights.C
  int binX = res.fResultsPro[AFO_var]->FindBin(valueX);
  if (tc.fInsanityCheckForEachParticle) // enable only during debugging, as this check is computationally heavy.
  {
    if (binX < 1) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      // underflow. this means that I didn't use the same cuts now, and I was using when making the particle weights. Adjust the cuts.
    }
    if (binX > res.fResultsPro[AFO_var]->GetNbinsX()) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      // overflow. this means that I didn't use the same cuts now, and I was using when making the particle weights. Adjust the cuts.
    }
  } // if(tc.fInsanityCheckForEachParticle)

  // *) Finally, determine weight for y(x):
  if (!pw.fDiffWeightsHist[AFO_diffWeight][binX - 1]) {
    LOGF(info, "\033[1;32mvalueY = %f\033[0m", valueY);
    LOGF(info, "\033[1;32mvalueX = %f\033[0m", valueX);
    LOGF(info, "\033[1;32mvariableX = %d\033[0m", static_cast<int>(variableX));
    LOGF(info, "\033[1;32mAFO_diffWeight = %d\033[0m", static_cast<int>(AFO_diffWeight));
    LOGF(info, "\033[1;32mbinX = %d\033[0m", binX);
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  int bin = pw.fDiffWeightsHist[AFO_diffWeight][binX - 1]->FindBin(valueY); // binX - 1, because I histogram for first bin in X is labeled with "[0]", etc.
  if (tc.fInsanityCheckForEachParticle)                                     // enable only during debugging, as this check is computationally heavy.
  {
    if (bin < 1) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      // underflow. this means that I didn't use the same cuts now, and I was using when making the particle weights. Adjust the cuts.
    }
    if (bin > pw.fDiffWeightsHist[AFO_diffWeight][binX - 1]->GetNbinsX()) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      // overflow. this means that I didn't use the same cuts now, and I was using when making the particle weights. Adjust the cuts.
    }
  } // if(tc.fInsanityCheckForEachParticle)

  double diffWeight = pw.fDiffWeightsHist[AFO_diffWeight][binX - 1]->GetBinContent(bin);
  if (tc.fInsanityCheckForEachParticle) // enable only during debugging, as this check is computationally heavy.
  {
    if (diffWeight < 0.) { // or <= 0 ? TBI 20240324 rethink
      LOGF(fatal, "\033[1;31m%s at line %d : diffWeight < 0\033[0m", __FUNCTION__, __LINE__);
    }
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return diffWeight;

} // DiffWeight(const double &valueY, const double &valueX, eqvectorKine variableX)

//============================================================

void GetParticleWeights()
{
  // Get the particle weights. Call this function only once.

  //    TBI 20231012 Here the current working assumption is that:
  //    1) Corrections do not change within a given run;
  //    2) Hyperloop proceeses the dataset one masterjob per run number.
  //    If any of these 2 assumptions are violated, this code will have to be modified.

  // a) Integrated weights;
  // b) Differential weights; => TBI 20250225 this is now obsolete and superseeded with c), where I use more general approach with sparse histograms
  // c) Differential phi weights using sparse histograms;
  // d) Differential pt weights using sparse histograms;
  // e) Differential eta weights using sparse histograms.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Integrated weights:
  // integrated phi weights:
  if (pw.fUseWeights[wPHI]) {
    TH1D* phiWeights = GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "phi");
    if (!phiWeights) {
      LOGF(fatal, "in function \033[1;31m%s at line %d, phiWeights is NULL. Check the external file %s with particle weights\033[0m", __FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
    }
    SetWeightsHist(phiWeights, wPHI);
  }

  // integrated pt weights:
  if (pw.fUseWeights[wPT]) {
    TH1D* ptWeights = GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "pt");
    if (!ptWeights) {
      LOGF(fatal, "\033[1;31m%s at line %d : ptWeights is NULL. Check the external file %s with particle weights\033[0m", __FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
    }
    SetWeightsHist(ptWeights, wPT);
  }

  // integrated eta weights:
  if (pw.fUseWeights[wETA]) {
    TH1D* etaWeights = GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "eta");
    if (!etaWeights) {
      LOGF(fatal, "\033[1;31m%s at line %d : etaWeights is NULL. Check the external file %s with particle weights\033[0m", __FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
    }
    SetWeightsHist(etaWeights, wETA);
  }

  // b) Differential weights:
  // differential phi(pt) weights:
  if (pw.fUseDiffWeights[wPHIPT]) {
    TH1D* phiptWeights = NULL;
    int nPtBins = res.fResultsPro[AFO_PT]->GetXaxis()->GetNbins();
    for (int b = 0; b < nPtBins; b++) {

      // *) check if particles in this pt bin survive particle cuts in pt. If not, skip this bin, because for that pt bin weights are simply not available:
      if (!(res.fResultsPro[AFO_PT]->GetBinLowEdge(b + 2) > pc.fdParticleCuts[ePt][eMin])) {
        // this branch protects against the case when I am e.g. in pt bin [0.0,0.2], and pt cut is 0.2 < pt < 5.0
        LOGF(info, "\033[1;33m%s at line %d : you are requesting phi(pt) weight for pt bin = %d from (%f,%f), which is outside (below) pt phase space = (%f,%f). Skipping this bin. \033[0m", __FUNCTION__, __LINE__, b, res.fResultsPro[AFO_PT]->GetBinLowEdge(b + 1), res.fResultsPro[AFO_PT]->GetBinLowEdge(b + 2), pc.fdParticleCuts[ePt][eMin], pc.fdParticleCuts[ePt][eMax]);
        continue;
      }
      if (!(res.fResultsPro[AFO_PT]->GetBinLowEdge(b + 1) < pc.fdParticleCuts[ePt][eMax])) {
        // this branch protects against the case when I am e.g. in pt bin [5.0,10.0], and pt cut is 0.2 < pt < 5.0
        LOGF(info, "\033[1;33m%s at line %d : you are requesting phi(pt) weight for pt bin = %d from (%f,%f), which is outside (above) pt phase space = (%f,%f). Skipping this bin. \033[0m", __FUNCTION__, __LINE__, b, res.fResultsPro[AFO_PT]->GetBinLowEdge(b + 1), res.fResultsPro[AFO_PT]->GetBinLowEdge(b + 2), pc.fdParticleCuts[ePt][eMin], pc.fdParticleCuts[ePt][eMax]);
        continue;
      }

      // *) okay, this pt bin is within pt phase-space window, defined by pt cut:
      phiptWeights = GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "phipt", b);
      if (!phiptWeights) {
        LOGF(fatal, "\033[1;31m%s at line %d : phiptWeights is NULL. Check the external file %s with particle weights\033[0m", __FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
      }

      // *) okay, just use this histogram with weights:
      SetDiffWeightsHist(phiptWeights, wPHIPT, b);
    }
  } // if (pw.fUseDiffWeights[wPHIPT]) {

  // differential phi(eta) weights:
  if (pw.fUseDiffWeights[wPHIETA]) {
    TH1D* phietaWeights = NULL;
    int nEtaBins = res.fResultsPro[AFO_ETA]->GetXaxis()->GetNbins();
    for (int b = 0; b < nEtaBins; b++) {

      // *) check if particles in this eta bin survive particle cuts in eta. If not, skip this bin, because for that eta bin weights are simply not available:
      if (!(res.fResultsPro[AFO_ETA]->GetBinLowEdge(b + 2) > pc.fdParticleCuts[eEta][eMin])) {
        // this branch protects against the case when I am e.g. in eta bin [-1.0,-0.8], and eta cut is -0.8 < eta < 0.8
        LOGF(info, "\033[1;33m%s at line %d : you are requesting phi(eta) weight for eta bin = %d from (%f,%f), which is outside (below) eta phase space = (%f,%f). Skipping this bin. \033[0m", __FUNCTION__, __LINE__, b, res.fResultsPro[AFO_ETA]->GetBinLowEdge(b + 1), res.fResultsPro[AFO_ETA]->GetBinLowEdge(b + 2), pc.fdParticleCuts[eEta][eMin], pc.fdParticleCuts[eEta][eMax]);
        continue;
      }
      if (!(res.fResultsPro[AFO_ETA]->GetBinLowEdge(b + 1) < pc.fdParticleCuts[eEta][eMax])) {
        // this branch protects against the case when I am e.g. in eta bin [0.8,1.0], and eta cut is 0.8 < eta < 1.0
        LOGF(info, "\033[1;33m%s at line %d : you are requesting phi(eta) weight for eta bin = %d from (%f,%f), which is outside (above) eta phase space = (%f,%f). Skipping this bin. \033[0m", __FUNCTION__, __LINE__, b, res.fResultsPro[AFO_ETA]->GetBinLowEdge(b + 1), res.fResultsPro[AFO_ETA]->GetBinLowEdge(b + 2), pc.fdParticleCuts[eEta][eMin], pc.fdParticleCuts[eEta][eMax]);
        continue;
      }

      // *) okay, this eta bin is within eta phase-space window, defined by eta cut:
      phietaWeights = GetHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), "phieta", b);
      if (!phietaWeights) {
        LOGF(fatal, "\033[1;31m%s at line %d : phietaWeights is NULL. Check the external file %s with particle weights\033[0m", __FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
      }

      // *) okay, just use this histogram with weights:
      SetDiffWeightsHist(phietaWeights, wPHIETA, b);
    } // for(int b=0; b<nEtaBins; b++) {
  } // if (pw.fUseDiffWeights[wPHIETA]) {

  // c) Differential phi weights using sparse histograms:
  if (pw.fUseDiffPhiWeights[wPhiPhiAxis]) { // yes, remember that flag for phi axis serves also as a common boolean to switch off all differential phi weights

    TString whichCategory = "phi"; // differential phi weights

    TString whichDimensions = ""; // differential phi weights as a function of particular dimension
    // Remark: the naming convention hardwired here for axes dimensions have to be in sync with what I have in the macro to make these weights
    if (pw.fUseDiffPhiWeights[wPhiPtAxis]) {
      whichDimensions += "_pt";
    }
    if (pw.fUseDiffPhiWeights[wPhiEtaAxis]) {
      whichDimensions += "_eta";
    }
    if (pw.fUseDiffPhiWeights[wPhiChargeAxis]) {
      whichDimensions += "_charge";
    }
    if (pw.fUseDiffPhiWeights[wPhiCentralityAxis]) {
      whichDimensions += "_centrality";
    }
    if (pw.fUseDiffPhiWeights[wPhiVertex_zAxis]) {
      whichDimensions += "_vertex_z";
    }
    // ...

    // TBI-today ... check if particles weights are avaiable for the phase window I have selected for each dimension with cuts

    THnSparseF* diffWeightsSparse = GetSparseHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), whichCategory.Data(), whichDimensions.Data());
    if (!diffWeightsSparse) {
      LOGF(fatal, "\033[1;31m%s at line %d : diffWeightsSparse  for category \"phi\" is NULL. Check the external file %s with particle weights\033[0m", __FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
    }

    // okay, just use this sparse histogram with weights:
    SetDiffWeightsSparse(diffWeightsSparse, eDWPhi);

  } // if (pw.fUseDiffPhiWeights[wPhiPhiAxis]) {

  // d) Differential pt weights using sparse histograms:
  if (pw.fUseDiffPtWeights[wPtPtAxis]) { // yes, remember that flag for pt axis serves also as a common boolean to switch off all differential pt weights

    TString whichCategory = "pt"; // differential pt weights

    TString whichDimensions = ""; // differential pt weights as a function of particular dimension
    // Remark: the naming convention hardwired here for axes dimensions have to be in sync with what I have in the macro to make these weights
    // ... TBI 20250222 proceed here in the same way as above for phi weights

    // TBI-today ... check if particles weights are avaiable for the phase window I have selected for each dimension with cuts

    THnSparseF* diffWeightsSparse = GetSparseHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), whichCategory.Data(), whichDimensions.Data());
    if (!diffWeightsSparse) {
      LOGF(fatal, "\033[1;31m%s at line %d : diffWeightsSparse for category \"pt\" is NULL. Check the external file %s with particle weights\033[0m", __FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
    }

    // okay, just use this sparse histogram with weights:
    SetDiffWeightsSparse(diffWeightsSparse, eDWPt);

  } // if (pw.fUseDiffPtWeights[wPtPtAxis]) {

  // e) Differential eta weights using sparse histograms:
  if (pw.fUseDiffEtaWeights[wEtaEtaAxis]) { // yes, remember that flag for eta axis serves also as a common boolean to switch off all differential eta weights

    TString whichCategory = "eta"; // differential eta weights

    TString whichDimensions = ""; // differential eta weights as a function of particular dimension
    // Remark: the naming convention hardwired here for axes dimensions have to be in sync with what I have in the macro to make these weights
    // ... TBI 20250222 proceed here in the same way as above for phi weights

    // TBI-today ... check if particles weights are avaiable for the phase window I have selected for each dimension with cuts

    THnSparseF* diffWeightsSparse = GetSparseHistogramWithWeights(pw.fFileWithWeights.Data(), tc.fRunNumber.Data(), whichCategory.Data(), whichDimensions.Data());
    if (!diffWeightsSparse) {
      LOGF(fatal, "\033[1;31m%s at line %d : diffWeightsSparse for category \"pt\" is NULL. Check the external file %s with particle weights\033[0m", __FUNCTION__, __LINE__, pw.fFileWithWeights.Data());
    }

    // okay, just use this sparse histogram with weights:
    SetDiffWeightsSparse(diffWeightsSparse, eDWEta);

  } // if (pw.fUseDiffEtaWeights[wEtaEtaAxis]) {

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void GetParticleWeights()

//============================================================

void GetCentralityWeights()
{
  // Get the centrality weights. Call this function only once.

  //    TBI 20231012 Here the current working assumption is that:
  //    1) Corrections do not change within a given run;
  //    2) Hyperloop proceeses the dataset one masterjob per run number.
  //    If any of these 2 assumptions are violated, this code will have to be modified.

  // a) Centrality weights;
  // b) ...

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Centrality weights:
  if (cw.fUseCentralityWeights) {
    TH1D* centralityWeights = GetHistogramWithCentralityWeights(cw.fFileWithCentralityWeights.Data(), tc.fRunNumber.Data());
    if (!centralityWeights) {
      LOGF(fatal, "in function \033[1;31m%s at line %d : centralityWeights is NULL. Check the external file %s with centrality weights\033[0m", __FUNCTION__, __LINE__, cw.fFileWithCentralityWeights.Data());
    }
    SetCentralityWeightsHist(centralityWeights);
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void GetCentralityWeights()

//============================================================

double CentralityWeight(const double& value) // centrality value
{
  // Determine centrality weight.

  // Ported directly from double AliAnalysisTaskMuPa::CentralityWeight(const double &value)

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  if (!cw.fCentralityWeightsHist) {
    LOGF(fatal, "\033[1;31m%s at line %d : cw.fCentralityWeightsHist is NULL. \033[0m", __FUNCTION__, __LINE__);
  }

  int bin = cw.fCentralityWeightsHist->FindBin(value);
  double weight = 0.;
  if (bin > cw.fCentralityWeightsHist->GetNbinsX()) {
    weight = 0.; // we are in the overflow, ignore this case
  } else {
    weight = cw.fCentralityWeightsHist->GetBinContent(bin) * cw.fCentralityWeightsHist->GetBinWidth(bin); // yes, since fCentralityWeightsHist is normalized p.d.f.
                                                                                                          // (I ensure that with the macro which makes centrality weights)
  }

  // In this context, it is assumed that centrality weight is a normalized probability (I ensure that with the macro which makes centrality weights):
  if (weight < 0. || weight > 1.) {
    LOGF(fatal, "\033[1;31m%s at line %d : weight = %f \033[0m", __FUNCTION__, __LINE__, weight);
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

  return weight;

} // double CentralityWeight(const double& value)

//============================================================

bool MaxNumberOfEvents(eBeforeAfter ba)
{
  // Check if max number of events was reached. Can be used for cut eNumberOfEvents (= total events, with ba = eBefore), and eSelectedEvents (ba = eAfter).

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Return value:
  bool reachedMaxNumberOfEvents = false;

  // *) Internal validation case (special treatment):
  if (iv.fUseInternalValidation) {
    if (eh.fEventHistograms[eNumberOfEvents][eSim][eAfter] && (eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]->GetBinContent(1) == ec.fdEventCuts[eNumberOfEvents][eMax] || eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]->GetBinContent(1) == ec.fdEventCuts[eSelectedEvents][eMax])) {
      return true;
    } else {
      return false;
    }
  }

  // *) Determine from which histogram the relevant info will be taken:
  int rs = -44;                                                  // reconstructed or simulated
  if (tc.fProcess[eGenericRec] || tc.fProcess[eGenericRecSim]) { // yes, for tc.fProcess[eGenericRecSim] I take info from Rec part
    rs = eRec;
  } else if (tc.fProcess[eGenericSim]) {
    rs = eSim;
  } else {
    LOGF(fatal, "\033[1;31m%s at line %d : not a single flag gProcess* is true \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Okay, do the thing:
  switch (ba) {
    case eBefore:
      if (eh.fEventHistograms[eNumberOfEvents][rs][eBefore] && eh.fEventHistograms[eNumberOfEvents][rs][eBefore]->GetBinContent(1) == ec.fdEventCuts[eNumberOfEvents][eMax]) {
        reachedMaxNumberOfEvents = true;
      }
      break;
    case eAfter:
      if (eh.fEventHistograms[eNumberOfEvents][rs][eAfter] && eh.fEventHistograms[eNumberOfEvents][rs][eAfter]->GetBinContent(1) == ec.fdEventCuts[eSelectedEvents][eMax]) {
        reachedMaxNumberOfEvents = true;
      }
      break;
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : enum ba = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(ba));
      break;
  }

  // *) Hasta la vista:
  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }
  return reachedMaxNumberOfEvents;

} // void MaxNumberOfEvents(eBeforeAfter ba)

//============================================================

void PrintEventCounter(eBeforeAfter ba)
{
  // Print how many events were processed by now.
  // Remark: If I am processing RecSim, the counter is corresponding to Rec.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Print or die:
  switch (ba) {
    case eBefore:
      if (!tc.fPlainPrintout) {
        LOGF(info, "\033[1;32m%s : processing event %d ....\033[0m", __FUNCTION__, eh.fEventCounter[eTotal]);
      } else {
        LOGF(info, "%s : processing event %d ....", __FUNCTION__, eh.fEventCounter[eTotal]);
      }
      break;
    case eAfter:
      if (!tc.fPlainPrintout) {
        LOGF(info, "\033[1;32m%s : event passed all cuts %d/%d\033[0m", __FUNCTION__, eh.fEventCounter[eProcessed], eh.fEventCounter[eTotal]);
      } else {
        LOGF(info, "%s : event passed all cuts %d/%d", __FUNCTION__, eh.fEventCounter[eProcessed], eh.fEventCounter[eTotal]);
      }
      break;
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : enum ba = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(ba));
      break;
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void PrintEventCounter(eBeforeAfter ba)

//============================================================

void EventCounterForDryRun(eEventCounterForDryRun eVar)
{
  // Simple utility function which either fills histogram with event count, or prints its current content.
  // Remark: Use only in combination with tc.fDryRun = true, otherwise I might be filling the same histogram in different member functions, there is a protection below.
  // It fills or prints per call, therefore I do not have to pass 'collision' objects, etc.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  if (!tc.fDryRun) {
    LOGF(fatal, "\033[1;31m%s at line %d : for the time being, function EventCounterForDryRun(...) can be safely used only for tc.fDryRun = true \033[0m", __FUNCTION__, __LINE__);
  }

  switch (eVar) {
    case eFill:
      // Fill event counter:
      !eh.fEventHistograms[eNumberOfEvents][eRec][eAfter] ? true : eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]->Fill(0.5);
      !eh.fEventHistograms[eNumberOfEvents][eSim][eAfter] ? true : eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]->Fill(0.5);
      break;
    case ePrint:
      // Print current status of event counter:
      // Remark: if I am processing RecSim, the counter is corresponding to Rec.
      if (eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]) {
        LOGF(info, "Processing event %d (dry run) ....", static_cast<int>(eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]->GetBinContent(1)));
      } else if (eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]) {
        LOGF(info, "Processing event %d (dry run) ....", static_cast<int>(eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]->GetBinContent(1)));
      }
      break;
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : enum eVar = %d is not supported yet in eEventCounter. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(eVar));
      break;
  } // switch(eVar)

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void EventCounterForDryRun(eEventCounterForDryRun eVar)

//============================================================

const char* FancyFormatting(const char* name)
{
  // Simple utility function to convert ordinary name into fancier formatting.

  // Examples:
  //  1. use LaTeX syntax (as supported by ROOT!), for the case when it's possible (e.g. "Phi" => "#{varphi}");
  //  2. add additional information to defalt name (e.g. "Centrality" => "Centrality (V0M)"
  //  3. ...

  if (tc.fVerboseUtility) {
    StartFunction(__FUNCTION__);
    LOGF(info, "\033[1;32m  const char* name = %s\033[0m", name);
  }

  // By default, do nothing and return the same thing:
  const char* fancyFormatting = name;

  // Special cases supported by now:
  if (TString(name).EqualTo("Phi", TString::kIgnoreCase)) {
    fancyFormatting = "#varphi";
  } else if (TString(name).EqualTo("Pt", TString::kIgnoreCase)) {
    fancyFormatting = "p_{T}";
  } else if (TString(name).EqualTo("Eta", TString::kIgnoreCase)) {
    fancyFormatting = "#eta";
  } else if (TString(name).EqualTo("Vertex_x")) {
    fancyFormatting = "V_{x}";
  } else if (TString(name).EqualTo("Vertex_y")) {
    fancyFormatting = "V_{y}";
  } else if (TString(name).EqualTo("Vertex_z")) {
    fancyFormatting = "V_{z}";
  } else if (TString(name).EqualTo("TotalMultiplicity")) {
    fancyFormatting = "TotalMultiplicity (tracks.size())";
  } else if (TString(name).EqualTo("Multiplicity", TString::kIgnoreCase)) {
    fancyFormatting = Form("Multiplicity (%s)", ec.fsEventCuts[eMultiplicityEstimator].Data());
  } else if (TString(name).EqualTo("ReferenceMultiplicity")) {
    fancyFormatting = Form("ReferenceMultiplicity (%s)", ec.fsEventCuts[eReferenceMultiplicityEstimator].Data());
  } else if (TString(name).EqualTo("Centrality", TString::kIgnoreCase)) {
    TString tmp = ec.fsEventCuts[eCentralityEstimator]; // I have to introduce local TString tmp, because ReplaceAll replaces in-place
    if (tmp.BeginsWith("CentRun2")) {
      fancyFormatting = Form("Centrality (%s)", tmp.ReplaceAll("CentRun2", "").Data()); // "CentRun2V0M" => "Centrality (V0M)"
    } else if (tmp.BeginsWith("Cent")) {
      fancyFormatting = Form("Centrality (%s)", tmp.ReplaceAll("Cent", "").Data()); // "CentFT0C" => "Centrality (FT0C)"
    } else {
      LOGF(fatal, "\033[1;31m%s at line %d : the case tmp = \"%s\" is not supported yet\033[0m", __FUNCTION__, __LINE__, tmp.Data());
    }
  } else if (TString(name).EqualTo("Trigger")) {
    fancyFormatting = Form("Trigger (%s)", ec.fsEventCuts[eTrigger].Data());
  } else if (TString(name).EqualTo("TrackOccupancyInTimeRange")) {
    fancyFormatting = "trackOccupancyInTimeRange()";
  } else if (TString(name).EqualTo("FT0COccupancyInTimeRange")) {
    fancyFormatting = "ft0cOccupancyInTimeRange()";
  } else if (TString(name).EqualTo("Occupancy", TString::kIgnoreCase)) {
    fancyFormatting = Form("Occupancy (%s)", ec.fsEventCuts[eOccupancyEstimator].Data());
  } else if (TString(name).EqualTo("InteractionRate", TString::kIgnoreCase) || TString(name).EqualTo("Interaction Rate", TString::kIgnoreCase)) {
    fancyFormatting = "Interaction Rate [kHz]"; // TBI 20241127 do I leave kHz hardwired here?
  } else if (TString(name).EqualTo("CurrentRunDuration", TString::kIgnoreCase) || TString(name).EqualTo("Current Run Duration", TString::kIgnoreCase)) {
    fancyFormatting = "Current run duration [s] (i.e. time in seconds since start of run)";
  }

  if (tc.fVerboseUtility) {
    ExitFunction(__FUNCTION__);
  }

  return fancyFormatting;

} // const char* FancyFormatting(const char *name)

//============================================================

double CalculateCustomNestedLoops(TArrayI* harmonics)
{
  // For the specified harmonics, get the correlation from nested loops.
  // Order of correlator is the number of harmonics, i.e. the number of elements in an array.

  // a) Determine the order of correlator;
  // b) Custom nested loop;
  // c) Return value.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  if (!harmonics) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  int nParticles = ebye.fSelectedTracks;
  /* TBI 20231108 enable eventually
  if(fUseFixedNumberOfRandomlySelectedParticles)
  {
   nParticles = 0;
   for(int i=0;i<nl.ftaNestedLoops[0]->GetSize();i++)
   {
    if(TMath::Abs(nl.ftaNestedLoops[0]->GetAt(i)) > 0. && TMath::Abs(nl.ftaNestedLoops[1]->GetAt(i)) > 0.){nParticles++;}
   }
  }
  */

  // a) Determine the order of correlator;
  int order = harmonics->GetSize();
  if (0 == order || order > gMaxCorrelator) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }
  if (nl.fMaxNestedLoop > 0 && nl.fMaxNestedLoop < order) {
    LOGF(info, "  nl.fMaxNestedLoop > 0 && nl.fMaxNestedLoop < order, where nl.fMaxNestedLoop = %d, order = %d", nl.fMaxNestedLoop, order);
    return 0.; // TBI 20240405 Is this really safe here? Re-think...
  }

  // b) Custom nested loop:
  TProfile* profile = new TProfile("profile", "", 1, 0., 1.); // helper profile to get all averages automatically
  // profile->Sumw2();
  double value = 0.;  // cos of current multiplet
  double weight = 1.; // weight of current multiplet
  for (int i1 = 0; i1 < nParticles; i1++) {
    double dPhi1 = nl.ftaNestedLoops[0]->GetAt(i1);
    double dW1 = nl.ftaNestedLoops[1]->GetAt(i1);
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
      double dPhi2 = nl.ftaNestedLoops[0]->GetAt(i2);
      double dW2 = nl.ftaNestedLoops[1]->GetAt(i2);
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
        double dPhi3 = nl.ftaNestedLoops[0]->GetAt(i3);
        double dW3 = nl.ftaNestedLoops[1]->GetAt(i3);
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
          double dPhi4 = nl.ftaNestedLoops[0]->GetAt(i4);
          double dW4 = nl.ftaNestedLoops[1]->GetAt(i4);
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
            double dPhi5 = nl.ftaNestedLoops[0]->GetAt(i5);
            double dW5 = nl.ftaNestedLoops[1]->GetAt(i5);
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
              double dPhi6 = nl.ftaNestedLoops[0]->GetAt(i6);
              double dW6 = nl.ftaNestedLoops[1]->GetAt(i6);
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
                double dPhi7 = nl.ftaNestedLoops[0]->GetAt(i7);
                double dW7 = nl.ftaNestedLoops[1]->GetAt(i7);
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
                  double dPhi8 = nl.ftaNestedLoops[0]->GetAt(i8);
                  double dW8 = nl.ftaNestedLoops[1]->GetAt(i8);
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
                    double dPhi9 = nl.ftaNestedLoops[0]->GetAt(i9);
                    double dW9 = nl.ftaNestedLoops[1]->GetAt(i9);
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
                      double dPhi10 = nl.ftaNestedLoops[0]->GetAt(i10);
                      double dW10 = nl.ftaNestedLoops[1]->GetAt(i10);
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
                        double dPhi11 = nl.ftaNestedLoops[0]->GetAt(i11);
                        double dW11 = nl.ftaNestedLoops[1]->GetAt(i11);
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
                          double dPhi12 = nl.ftaNestedLoops[0]->GetAt(i12);
                          double dW12 = nl.ftaNestedLoops[1]->GetAt(i12);
                          if (12 == order) {
                            value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5 + harmonics->GetAt(5) * dPhi6 + harmonics->GetAt(6) * dPhi7 + harmonics->GetAt(7) * dPhi8 + harmonics->GetAt(8) * dPhi9 + harmonics->GetAt(9) * dPhi10 + harmonics->GetAt(10) * dPhi11 + harmonics->GetAt(11) * dPhi12);
                            weight = dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8 * dW9 * dW10 * dW11 * dW12;
                            profile->Fill(0.5, value, weight);
                            continue;
                          }

                          // ... it's easy to continue the above pattern here

                        } // for(int i12=0; i12<nParticles; i12++)
                      } // for(int i11=0; i11<nParticles; i11++)
                    } // for(int i10=0; i10<nParticles; i10++)
                  } // for(int i9=0; i9<nParticles; i9++)
                } // for(int i8=0; i8<nParticles; i8++)
              } // for(int i7=0; i7<nParticles; i7++)
            } // for(int i6=0; i6<nParticles; i6++)
          } // for(int i5=0; i5<nParticles; i5++)
        } // for(int i4=0; i4<nParticles; i4++)
      } // for(int i3=0; i3<nParticles; i3++)
    } // for(int i2=0; i2<nParticles; i2++)
  } // for(int i1=0; i1<nParticles; i1++)

  // c) Return value:
  double finalValue = profile->GetBinContent(1);
  delete profile;
  profile = NULL;
  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }
  return finalValue;

} // double CalculateCustomNestedLoops(TArrayI *harmonics)

//============================================================

double CalculateKineCustomNestedLoops(TArrayI* harmonics, eAsFunctionOf AFO_variable, int bin)
{
  // For the specified harmonics, kine variable, and bin, get the correlation from nested loops.
  // Order of correlator is the number of harmonics, i.e. the number of elements in an array.

  // a) Determine the order of correlator;
  // b) Custom nested loop;
  // c) Return value.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  if (!harmonics) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // *) ...
  eqvectorKine qvKine = eqvectorKine_N; // which component of q-vector
  TString kineVarName = "";
  switch (AFO_variable) {
    case AFO_PT:
      qvKine = PTq;
      kineVarName = "pt";
      break;
    case AFO_ETA:
      qvKine = ETAq;
      kineVarName = "eta";
      break;
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : This AFO_variable = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(AFO_variable));
      break;
  } // switch(AFO_variable)

  // *) Insanity checks on above settings:
  if (qvKine == eqvectorKine_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : qvKine == eqvectorKine_N => add some more entries to the case statement \033[0m", __FUNCTION__, __LINE__);
  }

  if (0 > bin || res.fResultsPro[AFO_variable]->GetNbinsX() < bin) { // this 'bin' starts from 0, i.e. this is an array bin
    // either underflow or overflow is hit, meaning that histogram is booked in narrower range than cuts
    LOGF(fatal, "\033[1;31m%s at line %d => AFO_variable = %d, bin = %d\033[0m", __FUNCTION__, __LINE__, static_cast<int>(AFO_variable), bin);
  }

  // Get the number of particles in this kine bin:
  int nParticles = 0;
  for (int i = 0; i < nl.ftaNestedLoopsKine[qvKine][bin][0]->GetSize(); i++) {
    if (TMath::Abs(nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i)) > 0.) {
      nParticles++;
    }
  }

  // 'qvKine' is enum eqvectorKine:
  if (!res.fResultsPro[AFO_variable]) {
    LOGF(fatal, "\033[1;31m%s at line %d : AFO_variable = %d, bin = %d \033[0m", __FUNCTION__, __LINE__, static_cast<int>(AFO_variable), bin);
  }

  LOGF(info, " Processing qvKine = %d (vs. %s), nParticles in this kine bin = %d, bin range = [%f,%f) ....", static_cast<int>(qvKine), kineVarName.Data(), nParticles, res.fResultsPro[AFO_variable]->GetBinLowEdge(bin + 1), res.fResultsPro[AFO_variable]->GetBinLowEdge(bin + 2));

  // a) Determine the order of correlator;
  int order = harmonics->GetSize();
  if (0 == order || order > gMaxCorrelator) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }
  if (order > nParticles) {
    LOGF(info, "  There is no enough particles in this bin to calculate the requested correlator");
    return 0.; // TBI 20240405 Is this really safe here? Re-think...
  }
  if (nl.fMaxNestedLoop > 0 && nl.fMaxNestedLoop < order) {
    LOGF(info, "  nl.fMaxNestedLoop > 0 && nl.fMaxNestedLoop < order, where nl.fMaxNestedLoop = %d, order = %d", nl.fMaxNestedLoop, order);
    return 0.; // TBI 20240405 Is this really safe here? Re-think...
  }

  // b) Custom nested loop:
  TProfile* profile = new TProfile("profile", "", 1, 0., 1.); // helper profile to get all averages automatically
  // profile->Sumw2();
  double value = 0.;  // cos of current multiplet
  double weight = 1.; // weight of current multiplet
  for (int i1 = 0; i1 < nParticles; i1++) {
    double dPhi1 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i1);
    double dW1 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i1);
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
      double dPhi2 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i2);
      double dW2 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i2);
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
        double dPhi3 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i3);
        double dW3 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i3);
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
          double dPhi4 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i4);
          double dW4 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i4);
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
            double dPhi5 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i5);
            double dW5 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i5);
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
              double dPhi6 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i6);
              double dW6 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i6);
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
                double dPhi7 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i7);
                double dW7 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i7);
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
                  double dPhi8 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i8);
                  double dW8 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i8);
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
                    double dPhi9 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i9);
                    double dW9 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i9);
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
                      double dPhi10 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i10);
                      double dW10 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i10);
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
                        double dPhi11 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i11);
                        double dW11 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i11);
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
                          double dPhi12 = nl.ftaNestedLoopsKine[qvKine][bin][0]->GetAt(i12);
                          double dW12 = nl.ftaNestedLoopsKine[qvKine][bin][1]->GetAt(i12);
                          if (12 == order) {
                            value = TMath::Cos(harmonics->GetAt(0) * dPhi1 + harmonics->GetAt(1) * dPhi2 + harmonics->GetAt(2) * dPhi3 + harmonics->GetAt(3) * dPhi4 + harmonics->GetAt(4) * dPhi5 + harmonics->GetAt(5) * dPhi6 + harmonics->GetAt(6) * dPhi7 + harmonics->GetAt(7) * dPhi8 + harmonics->GetAt(8) * dPhi9 + harmonics->GetAt(9) * dPhi10 + harmonics->GetAt(10) * dPhi11 + harmonics->GetAt(11) * dPhi12);
                            weight = dW1 * dW2 * dW3 * dW4 * dW5 * dW6 * dW7 * dW8 * dW9 * dW10 * dW11 * dW12;
                            profile->Fill(0.5, value, weight);
                            continue;
                          }

                          // ... it's easy to continue the above pattern here

                        } // for(int i12=0; i12<nParticles; i12++)
                      } // for(int i11=0; i11<nParticles; i11++)
                    } // for(int i10=0; i10<nParticles; i10++)
                  } // for(int i9=0; i9<nParticles; i9++)
                } // for(int i8=0; i8<nParticles; i8++)
              } // for(int i7=0; i7<nParticles; i7++)
            } // for(int i6=0; i6<nParticles; i6++)
          } // for(int i5=0; i5<nParticles; i5++)
        } // for(int i4=0; i4<nParticles; i4++)
      } // for(int i3=0; i3<nParticles; i3++)
    } // for(int i2=0; i2<nParticles; i2++)
  } // for(int i1=0; i1<nParticles; i1++)

  // c) Return value:
  double finalValue = profile->GetBinContent(1);
  delete profile;
  profile = NULL;
  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }
  return finalValue;

} // double CalculateKineCustomNestedLoops(TArrayI *harmonics, eAsFunctionOf AFO_variable, int bin)

//============================================================

void DetermineMultiplicity()
{
  // Determine multiplicity for "vs. mult" results.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  if (ec.fsEventCuts[eMultiplicityEstimator].EqualTo("SelectedTracks", TString::kIgnoreCase)) {
    ebye.fMultiplicity = static_cast<float>(ebye.fSelectedTracks);
  } else if (ec.fsEventCuts[eMultiplicityEstimator].EqualTo("ReferenceMultiplicity", TString::kIgnoreCase)) {
    ebye.fMultiplicity = ebye.fReferenceMultiplicity;
  } else {
    LOGF(fatal, "\033[1;31m%s at line %d : multiplicity estimator = %s is not supported yet. \033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eMultiplicityEstimator].Data());
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void DetermineMultiplicity()

//============================================================

template <eRecSim rs, typename T>
void DetermineReferenceMultiplicity(T const& collision)
{
  // Determine collision reference multiplicity.

  // a) Determine reference multiplicity for real Run 3 data;
  // b) Determine reference multiplicity for simulated Run 3 data;
  // c) Same as a), just for converted Run 2 and Run 1 data;
  // d) Same as b), just for converted Run 2 and Run 1 data;
  // e) Test case;
  // f) Print reference multiplicity  for the audience...

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Determine reference multiplicity for real Run 3 data:
  if constexpr (rs == eRec || rs == eRecAndSim) {
    // Local convention for name of reference multiplicity estimator: use the same name as the getter, case insensitive.
    if (ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("multTPC", TString::kIgnoreCase)) {
      ebye.fReferenceMultiplicity = collision.multTPC();
    } else if (ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("multFV0M", TString::kIgnoreCase)) {
      ebye.fReferenceMultiplicity = collision.multFV0M();
    } else if (ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("multFT0C", TString::kIgnoreCase)) {
      ebye.fReferenceMultiplicity = collision.multFT0C();
    } else if (ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("multFT0M", TString::kIgnoreCase)) {
      ebye.fReferenceMultiplicity = collision.multFT0M();
    } else if (ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("multNTracksPV", TString::kIgnoreCase)) {
      ebye.fReferenceMultiplicity = collision.multNTracksPV();
    } else if (ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("multNTracksGlobal", TString::kIgnoreCase)) {
      // ebye.fReferenceMultiplicity = collision.multNTracksGlobal(); // TBI 20241209 not validated yet
    } else {
      LOGF(fatal, "\033[1;31m%s at line %d : reference multiplicity estimator = %d is not supported yet for Run 3. \033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eReferenceMultiplicityEstimator].Data());
    }
    // QA:
    if (qa.fFillQAEventHistograms2D) { // TBI 20240515 this flag is too general here, I need to make it more specific
      qa.fReferenceMultiplicity[eMultTPC] = collision.multTPC();
      qa.fReferenceMultiplicity[eMultFV0M] = collision.multFV0M();
      qa.fReferenceMultiplicity[eMultFT0C] = collision.multFT0C();
      qa.fReferenceMultiplicity[eMultFT0M] = collision.multFT0M();
      qa.fReferenceMultiplicity[eMultNTracksPV] = collision.multNTracksPV();
      // qa.fReferenceMultiplicity[eMultNTracksGlobal] = collision.multNTracksGlobal(); // TBI 20241209 not validated yet
    }

    // TBI 20241123 check if corresponding simulated ref. mult. is available through collision.has_mcCollision()
    // ...
  }

  // b) Determine reference multiplicity for simulated Run 3 data:
  if constexpr (rs == eSim) {
    ebye.fReferenceMultiplicity = -44.; // TBI 20241123 check what to use here and add support eventualy
  }

  // c) Same as a), just for converted Run 2 and Run 1 data:
  if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2 || rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    if (ec.fsEventCuts[eReferenceMultiplicityEstimator].EqualTo("multTracklets", TString::kIgnoreCase)) {
      ebye.fReferenceMultiplicity = collision.multTracklets();
    } else {
      LOGF(fatal, "\033[1;31m%s at line %d : reference multiplicity estimator = %d is not supported yet for Run 2. \033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eReferenceMultiplicityEstimator].Data());
    }
    // QA:
    if (qa.fFillQAEventHistograms2D) { // TBI 20240515 this flag is too general here, I need to make it more specific
      // ...
    }

    // TBI 20241123 check if corresponding simulated ref. mult. is available through collision.has_mcCollision()
    // ...
  }

  // d) Same as b), just for converted Run 2 and Run 1 data:
  if constexpr (rs == eSim_Run2 || rs == eSim_Run1) {
    ebye.fReferenceMultiplicity = -44.; // TBI 20241123 check what to use here and add support eventualy
  }

  // e) Test case:
  if constexpr (rs == eTest) {
    ebye.fReferenceMultiplicity = static_cast<float>(gRandom->Uniform(0., 5000.)); // TBI 20241123 I could implement here a getter, if there is one available both for Run 3 and Run 2/1
  }

  // f) Print centrality for the audience...:
  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m ebye.fReferenceMultiplicity = %f\033[0m", ebye.fReferenceMultiplicity);
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T> void DetermineReferenceMultiplicity(T const& collision)

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
  // g) Test case;
  // h) Print centrality for the audience...

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) For real data, determine centrality from default centrality estimator:
  if constexpr (rs == eRec || rs == eRecAndSim) {
    // Local convention for name of centrality estimator: use the same name as the getter, case insensitive.
    if (ec.fsEventCuts[eCentralityEstimator].EqualTo("centFT0C", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centFT0C();
    } else if (ec.fsEventCuts[eCentralityEstimator].EqualTo("centFT0CVariant1", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centFT0CVariant1();
    } else if (ec.fsEventCuts[eCentralityEstimator].EqualTo("centFT0M", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centFT0M();
    } else if (ec.fsEventCuts[eCentralityEstimator].EqualTo("centFV0A", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centFV0A();
    } else if (ec.fsEventCuts[eCentralityEstimator].EqualTo("centNTPV", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centNTPV();
    } else if (ec.fsEventCuts[eCentralityEstimator].EqualTo("centNGlobal", TString::kIgnoreCase)) {
      // ebye.fCentrality = collision.centNGlobal(); // TBI 20250128 enable eventually
    } else {
      LOGF(fatal, "\033[1;31m%s at line %d : centrality estimator = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eCentralityEstimator].Data());
    }
    // QA:
    if (qa.fFillQAEventHistograms2D) {
      qa.fCentrality[eCentFT0C] = collision.centFT0C();
      qa.fCentrality[eCentFT0CVariant1] = collision.centFT0CVariant1();
      qa.fCentrality[eCentFT0M] = collision.centFT0M();
      qa.fCentrality[eCentFV0A] = collision.centFV0A();
      qa.fCentrality[eCentNTPV] = collision.centNTPV();
      // qa.fCentrality[eCentNGlobal] = collision.centNGlobal(); // TBI 20250128 enable eventually
    }

    // TBI 20240120 I could also here access also corresponding simulated centrality from impact parameter, if available through collision.has_mcCollision()
  }

  // b) For simulated data, determine centrality directly from impact parameter:
  if constexpr (rs == eSim) {
    ebye.fCentrality = -44.; // TBI 20240120 add support eventualy
  }

  // c) Same as a), just for converted Run 2 data:
  if constexpr (rs == eRec_Run2 || rs == eRecAndSim_Run2) {
    if (ec.fsEventCuts[eCentralityEstimator].EqualTo("centRun2V0M", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centRun2V0M();
    } else if (ec.fsEventCuts[eCentralityEstimator].EqualTo("centRun2SPDTracklets", TString::kIgnoreCase)) {
      ebye.fCentrality = collision.centRun2SPDTracklets();
    } else {
      LOGF(fatal, "\033[1;31m%s at line %d : centrality estimator = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eCentralityEstimator].Data());
    }
    // QA:
    if (qa.fFillQAEventHistograms2D) { // TBI 20240515 this flag is too general here, I need to make it more specific
      qa.fCentrality[eCentRun2V0M] = collision.centRun2V0M();
      qa.fCentrality[eCentRun2SPDTracklets] = collision.centRun2SPDTracklets();
    }

    // TBI 20240120 I could also here access also corresponding simulated centrality from impact parameter, if available through collision.has_mcCollision()
  }

  // d) Same as b), just for converted Run 2 data:
  if constexpr (rs == eSim_Run2) {
    ebye.fCentrality = -44.; // TBI 20240120 add support eventualy
  }

  // e) Same as a), just for converted Run 1 data:
  if constexpr (rs == eRec_Run1 || rs == eRecAndSim_Run1) {
    if (ec.fsEventCuts[eCentralityEstimator].EqualTo("centRun2V0M", TString::kIgnoreCase)) {
      // ebye.fCentrality = collision.centRun2V0M(); // TBI 20240224 enable when I add support for RecAndSim_Run1
    } else if (ec.fsEventCuts[eCentralityEstimator].EqualTo("CentRun2SPDTracklets", TString::kIgnoreCase)) {
      // ebye.fCentrality = collision.centRun2SPDTracklets(); // TBI 20240224 enable when I add support for RecAndSim_Run1
    } else {
      LOGF(fatal, "\033[1;31m%s at line %d : centrality estimator = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eCentralityEstimator].Data());
    }
    // TBI 20240120 I could also here access also corresponding simulated centrality from impact parameter, if available through collision.has_mcCollision()
  }

  // f) Same as b), just for converted Run 1 data:
  if constexpr (rs == eSim_Run1) {
    ebye.fCentrality = -44.; // TBI 20240515 add support eventualy, or merge with Run 2 branch. It seems that in converted Run 1 there is no centrality.
  }

  // g) Test case:
  if constexpr (rs == eTest) {
    ebye.fCentrality = static_cast<float>(gRandom->Uniform(0., 100.));
  }

  // h) Print centrality for the audience...:
  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m ebye.fCentrality = %f\033[0m", ebye.fCentrality);
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T> void DetermineCentrality(T const& collision)

//============================================================

template <eRecSim rs, typename T>
void DetermineOccupancy(T const& collision)
{
  // Determine collision occupancy.

  // a) Determine occupancy from default occupancy estimator, only for eRec and eRecAndSim;
  // b) For all other cases, set occupancy to -1 (not defined).
  // c) Print occupancy for the audience...

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Determine occupancy from default occupancy estimator, only for eRec and eRecAndSim:
  if constexpr (rs == eRec || rs == eRecAndSim) {
    if (ec.fsEventCuts[eOccupancyEstimator].EqualTo("TrackOccupancyInTimeRange", TString::kIgnoreCase)) {
      ebye.fOccupancy = collision.trackOccupancyInTimeRange();
    } else if (ec.fsEventCuts[eOccupancyEstimator].EqualTo("FT0COccupancyInTimeRange", TString::kIgnoreCase)) {
      ebye.fOccupancy = collision.ft0cOccupancyInTimeRange();
    } else {
      LOGF(fatal, "\033[1;31m%s at line %d : occupancy estimator = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, ec.fsEventCuts[eOccupancyEstimator].Data());
    }
    // QA:
    if (qa.fFillQAEventHistograms2D) { // TBI 20240515 this flag is too general here, I need to make it more specific
      qa.fOccupancy[eTrackOccupancyInTimeRange] = collision.trackOccupancyInTimeRange();
      qa.fOccupancy[eFT0COccupancyInTimeRange] = collision.ft0cOccupancyInTimeRange();
    }
  } else {
    // b) For all other cases, set occupancy to -1 (not defined):
    ebye.fOccupancy = -1.;
    // QA:
    if (qa.fFillQAEventHistograms2D) { // TBI 20240515 this flag is too general here, I need to make it more specific
      for (int oe = 0; oe < eOccupancyEstimators_N; oe++) {
        qa.fOccupancy[oe] = -1.;
      }
    }
  }

  // c) Print occupancy for the audience...:
  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m ebye.fOccupancy = %f\033[0m", ebye.fOccupancy);
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T> void DetermineOccupancy(T const& collision)

//============================================================

template <eRecSim rs, typename T1, typename T2>
void DetermineInteractionRateAndCurrentRunDuration(T1 const& collision, T2 const&)
{
  // Determine interaction rate and current run duration in Run 3.

  // Cannot be used in converted Run 2 and Run 1, because mRateFetcher.fetch... line below crashes with example line:
  //    [228607:multiparticle-correlations-a-b]: [10:02:38][ERROR] Requested resource does not exist: http://alice-ccdb.cern.ch//GLO/Config/GRPLHCIF/1449947476529/
  //    [228607:multiparticle-correlations-a-b]: [10:02:38][FATAL] Got nullptr from CCDB for path GLO/Config/GRPLHCIF and timestamp 1449947476529

  // a) Determine interaction rate and current run duration only for eRec;
  // b) For all other cases, set interaction rate to -1 for the time being;
  // c) Print interaction rate and current run duration for the audience...

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a1) Determine interaction rate only for eRec:
  if constexpr (rs == eRec) {                      // TBI 20250112 check still eRecSim mode here
    auto bc = collision.template foundBC_as<T2>(); // I have the same code snippet at other places, keep in sync.
    double hadronicRate = mRateFetcher.fetch(ccdb.service, static_cast<uint64_t>(bc.timestamp()), static_cast<int>(bc.runNumber()), "ZNC hadronic") * 1.e-3;
    if (hadronicRate > 0.) {
      ebye.fInteractionRate = static_cast<float>(hadronicRate);
    } else {
      LOGF(warning, "\033[1;31m%s at line %d : hadronicRate = %f is meaningless \033[0m", __FUNCTION__, __LINE__, hadronicRate);
      // I hit indeed at negative hadronic rate in LHC24ar/559545/apass1 dataset. But I do not really need to bail out here, because that collision in
      // any case will not pass a cut in configurable cfInteractionRate . Therefore, I print a warning, and then can grep it from the log, if necessary.
    }

    // a2) Determine the current run duration:
    // TBI 20250107 I could move this to a separate function?
    ebye.fCurrentRunDuration = std::floor(bc.timestamp() * 0.001) - tc.fRunTime[eStartOfRun];
    if (ebye.fCurrentRunDuration > tc.fRunTime[eDurationInSec]) {
      LOGF(fatal, "\033[1;31m%s at line %d : ebye.fCurrentRunDuration = %d is bigger than tc.fRunTime[eDurationInSec] = %d, which is meaningless \033[0m", __FUNCTION__, __LINE__, static_cast<int>(ebye.fCurrentRunDuration), static_cast<int>(tc.fRunTime[eDurationInSec]));
    }
  } else {
    // b) For all other cases, set interaction rate to -1:
    ebye.fInteractionRate = -1.;
    ebye.fCurrentRunDuration = -1.;
  }

  // c) Print interaction rate and run duration for the audience...:
  if (tc.fVerbose) {
    LOGF(info, "\033[1;32m ebye.fInteractionRate = %f kHz\033[0m", ebye.fInteractionRate);
    if (qa.fBookQAEventHistograms2D[eCurrentRunDuration_vs_InteractionRate]) { // TBI 20241127 do I check this flag, or pointer, like in FillEventHistograms(...) ?
      LOGF(info, "\033[1;32m ebye.fCurrentRunDuration = %f s (in seconds after SOR)\033[0m", ebye.fCurrentRunDuration);
    }
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T1, typename T2> void DetermineInteractionRateAndCurrentRunDuration(T1 const& collision, T2 const& bcs)

//============================================================

template <eRecSim rs, typename T>
void DetermineVertexZ(T const& collision)
{
  // Determine vetex z position.

  // TBI 20250108 I could use ebye.fVz determined here to fill event histograms, but it's not a big deal to fetch it there also via collision.posZ()

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  ebye.fVz = collision.posZ();

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void DetermineVertexZ(T const& collision)

//============================================================

void DetermineEventCounters()
{
  // Determine all event counters.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // Remark: For "RecSim", the total number of events is taken from eRec.
  if (eh.fEventHistograms[eNumberOfEvents][eRec][eBefore] && eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]) {
    eh.fEventCounter[eTotal] = static_cast<int>(eh.fEventHistograms[eNumberOfEvents][eRec][eBefore]->GetBinContent(1));
    eh.fEventCounter[eProcessed] = static_cast<int>(eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]->GetBinContent(1));
  } else if (eh.fEventHistograms[eNumberOfEvents][eSim][eBefore] && eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]) {
    // Remark: This branch covers automatically also internal validation, because I book and fill there only eSim.
    eh.fEventCounter[eTotal] = static_cast<int>(eh.fEventHistograms[eNumberOfEvents][eSim][eBefore]->GetBinContent(1));
    eh.fEventCounter[eProcessed] = static_cast<int>(eh.fEventHistograms[eNumberOfEvents][eSim][eAfter]->GetBinContent(1));
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void DetermineEventCounters()

//============================================================

void RandomIndices(int nTracks)
{
  // Randomize indices using Fisher-Yates algorithm.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  if (nTracks < 1) {
    return;
  }

  // Fisher-Yates algorithm:
  tc.fRandomIndices = new TArrayI(nTracks);
  tc.fRandomIndices->Reset(); // just in case there is some random garbage in memory at init
  for (int i = 0; i < nTracks; i++) {
    tc.fRandomIndices->AddAt(i, i);
  }
  for (int i = nTracks - 1; i >= 1; i--) {
    int j = gRandom->Integer(i + 1);
    int temp = tc.fRandomIndices->GetAt(j);
    tc.fRandomIndices->AddAt(tc.fRandomIndices->GetAt(i), j);
    tc.fRandomIndices->AddAt(temp, i);
  } // end of for(int i=nTracks-1;i>=1;i--)

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void RandomIndices(int nTracks)

//============================================================

template <eRecSim rs, typename T>
void BanishmentLoopOverParticles(T const& tracks)
{
  // This is the quick banishment loop over particles, as a support for eSelectedTracks cut (used through eMultiplicity, see comments for ebye.fMultiplicity).
  // This is particularly relevant to get all efficiency corrections right.
  // The underlying problem is that particle histograms got filled before eSelectedTracks could be applied in Steer.
  // Therefore, particle histograms got filled even for events which were rejected by eSelectedTracks cut.
  // In this loop, for those few specific events (typically low-multiplicity outliers), particle histograms are re-filled again with weight -1,
  // which in effect cancels the previos fill in the MainLoopOverParticles.

  // Remark: I have to use here all additional checks, like ValidTrack, as in the MainLoopOverParticles.
  //         Therefore, it's important to have local variable lSelectedTracks, so that I can cross-compare
  //         at the end with central data member ebye.fSelectedTracks

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) If random access of tracks from collection is requested, use Fisher-Yates algorithm to generate random indices:
  //    Remark: It is very important that I use exactly the same random sequence from FY already generated in the MainLoop
  if (tc.fUseFisherYates) {
    if (!tc.fRandomIndices) {
      LOGF(fatal, "\033[1;31m%s at line %d : I have to use here exactly the same random sequence from FY already generated in the MainLoopOverParticles, but it's not available \033[0m", __FUNCTION__, __LINE__);
    }
  }

  // *) Counter of selected tracks in the current event:
  int lSelectedTracks = 0; // I could reset and reuse here ebye.fSelectedTracks, but it's safer to use separate local variable, as I can do additional insanity checks here

  // *) Banishment loop over particles:
  // for (auto& track : tracks) { // default standard way of looping of tracks
  auto track = tracks.iteratorAt(0); // set the type and scope from one instance
  for (int64_t i = 0; i < tracks.size(); i++) {

    // *) Access track sequentially from collection of tracks (default), or randomly using Fisher-Yates algorithm:
    if (!tc.fUseFisherYates) {
      track = tracks.iteratorAt(i);
    } else {
      track = tracks.iteratorAt(static_cast<int64_t>(tc.fRandomIndices->GetAt(i)));
    }

    // *) Skip track objects which are not valid tracks (e.g. Run 2 and 1 tracklets, etc.):
    if (!ValidTrack<rs>(track)) {
      continue;
    }

    // *) Banish particle histograms before particle cuts:
    // TBI 20240515 I banish for the time being only particle histograms AFTER particle cuts.
    //              If I start to banish here also particle histograms BEFORE particle cuts, then see if I have to do it also for event histograms BEFORE cuts.
    //              Event histograms AFTER cuts are not affected.
    //    if (ph.fFillParticleHistograms || ph.fFillParticleHistograms2D) {
    //      FillParticleHistograms<rs>(track, eBefore, -1);
    //    }

    // *) Particle cuts:
    if (!ParticleCuts<rs>(track, eCut)) { // Main call for particle cuts.
      continue;                           // not return!!
    }

    // *) Banish particle histograms after particle cuts:
    if (ph.fFillParticleHistograms || ph.fFillParticleHistograms2D || qa.fFillQAParticleHistograms2D) {
      FillParticleHistograms<rs>(track, eAfter, -1); // with negative weight -1, I effectively remove the previous fill for this track
    }

    // *) Increase the local selected particle counter:
    lSelectedTracks++;
    if (lSelectedTracks >= ec.fdEventCuts[eMultiplicity][eMax]) {
      break;
    }

    // *) Break the loop if fixed number of particles is taken randomly from each event (use always in combination with tc.fUseFisherYates = true):
    if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && tc.fFixedNumberOfRandomlySelectedTracks == lSelectedTracks) {
      LOGF(info, "%s : Breaking the loop over particles, since requested fixed number of %d particles was reached", __FUNCTION__, tc.fFixedNumberOfRandomlySelectedTracks);
      break;
    }

  } // for (auto& track : tracks)

  // *) Quick insanity checks (mandatory!):
  if (lSelectedTracks != ebye.fSelectedTracks) {
    LOGF(fatal, "\033[1;31m%s at line %d : lSelectedTracks != ebye.fSelectedTracks , lSelectedTracks = %d, ebye.fSelectedTracks = %d \033[0m", __FUNCTION__, __LINE__, lSelectedTracks, ebye.fSelectedTracks);
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T> void BanishmentLoopOverParticles(T const& tracks) {

//============================================================

void PrintCutCounterContent()
{
  // Prints on the screen content of fEventCutCounterHist[][] (all which were booked).

  // a) Insanity checks;
  // b) Print or die.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // a) Insanity checks:
  if (!(ec.fUseEventCutCounterAbsolute || ec.fUseEventCutCounterSequential)) {
    LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
  }

  // b) Print or die:
  for (int rs = 0; rs < 2; rs++) // reco/sim
  {
    for (int cc = 0; cc < eCutCounter_N; cc++) // enum eCutCounter
    {
      if (!(ec.fEventCutCounterHist[rs][cc])) {
        continue;
      }
      LOGF(info, "\033[1;32m\nPrinting the content of event cut counter histogram %s\033[0m", ec.fEventCutCounterHist[rs][cc]->GetName());
      for (int bin = 1; bin <= ec.fEventCutCounterHist[rs][cc]->GetNbinsX(); bin++) {
        if (TString(ec.fEventCutCounterHist[rs][cc]->GetXaxis()->GetBinLabel(bin)).EqualTo("TBI")) { // TBI 20240514 temporary workaround, "TBI" can't persist here
          continue;
        }
        LOGF(info, "bin = %d => %s : %d", bin, ec.fEventCutCounterHist[rs][cc]->GetXaxis()->GetBinLabel(bin), static_cast<int>(ec.fEventCutCounterHist[rs][cc]->GetBinContent(bin)));
      }
    } // for (int cc = 0; cc < eCutCounter_N; cc++) // enum eCutCounter
  } // for (int rs = 0; rs < 2; rs++) // reco/sim

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void PrintCutCounterContent()

//============================================================

void Trace(const char* functionName, int lineNumber)
{
  // A simple utility wrapper. Use only during debugging, sprinkle calls to this function here and there, as follows
  //    Trace(__FUNCTION__, __LINE__);

  LOGF(info, "\033[1;32m%s .... line %d\033[0m", functionName, lineNumber);

} // void Trace(const char* functionName, int lineNumber)

//============================================================

void Exit()
{
  // A simple utility wrapper. Used only during debugging.
  // Use directly as:  Exit();
  // Line number, function name, formatting, etc, are determinad automatically.

  LOGF(info, "\n\n\n\n\n\n\n\n\n\n");
  exit(1);

} // void Exit()

//============================================================

void StartFunction(const char* functionName)
{
  // A simple utility wrapper, used when tc.fVerbose = true. It merely ensures uniform formatting of notification when the function starts.

  LOGF(info, "\033[1;32mStart %s\033[0m", functionName); // prints in green

} // void StartFunction(const char* functionName)

//============================================================

void ExitFunction(const char* functionName)
{
  // A simple utility wrapper, used when tc.fVerbose = true. It merely ensures uniform formatting of notification when the function exits.

  LOGF(info, "\033[1;32mExit %s\033[0m", functionName); // prints in green

} // void ExitFunction(const char* functionName)

//============================================================

void BailOut(bool finalBailout = false)
{
  // Use only locally - bail out if maximum number of events was reached, and dump all results by that point in a local ROOT file.
  // If fSequentialBailout > 0, bail out is performed each fSequentialBailout events, each time in a new local ROOT file.
  // For sequential bailout, the naming scheme of ROOT files is AnalysisResultsBailOut_eh.fEventCounter[eProcessed].root .
  // If ROOT file with the same name already exists, BailOut is not performed, since the argument is that
  // it's pointless to perform Bailout for same eh.fEventCounter[eProcessed], even if eh.fEventCounter[eTotal] changed.
  // Only if finalBailout = true, I will overwrite the existing file with the same name.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Local variables: TBI 20240130 shall I promote 'em to data members + add support for configurables?
  TString sBailOutFile = "AnalysisResultsBailOut.root";
  TString sDirectoryFile = "multiparticle-correlations-a-b";

  // *) For sequential bailout, I need to adapt the ROOT file name each time this function is called:
  if (tc.fSequentialBailout > 0) {
    sBailOutFile.ReplaceAll(".root", Form("_%d.root", eh.fEventCounter[eProcessed])); // replaces in-place
    // basically, at 1st call "AnalysisResultsBailOut.root" => "AnalysisResultsBailOut_1*eh.fEventCounter[eProcessed].root",
    //            at 2nd call "AnalysisResultsBailOut.root" => "AnalysisResultsBailOut_2*eh.fEventCounter[eProcessed].root", etc.
    if (!finalBailout && !gSystem->AccessPathName(sBailOutFile.Data(), kFileExists)) { // only for finalBailout = true, I will overwrite the existing file with the same name.
      LOGF(info, "\033[1;33m\nsBailOutFile = %s already exits, that means that eh.fEventCounter[eProcessed] is the same as in the previous call of BailOut.\nJust skipping and waiting more events to pass selection criteria... \033[0m", sBailOutFile.Data());
      return;
    }
  }

  // *) Info message:
  if (eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]) {
    LOGF(info, "\033[1;32m=> Per request, bailing out after %d selected events in the local file %s .\n\033[0m", static_cast<int>(eh.fEventHistograms[eNumberOfEvents][eRec][eAfter]->GetBinContent(1)), sBailOutFile.Data());
  }

  // *) Okay, let's bail out intentionally:
  TFile* f = new TFile(sBailOutFile.Data(), "recreate");
  TDirectoryFile* dirFile = new TDirectoryFile(sDirectoryFile.Data(), sDirectoryFile.Data());
  // TBI 20240130 I cannot add here fBaseList directly, since that one is declared as OutputObj<TList>
  // Therefore, adding one-by-one nested TList's I want to bail out.
  // Keep in sync with BookAndNestAllLists().
  TList* bailOutList = new TList(); // this is sort of 'fake' fBaseList
  bailOutList->SetOwner(false);     // yes, beacause for sequential bailout, with SetOwner(true) the code is crashing after 1st sequential bailout is done
  bailOutList->SetName(sBaseListName.Data());
  bailOutList->Add(fBasePro); // yes, this one needs a special treatment
  bailOutList->Add(qa.fQAList);
  bailOutList->Add(ec.fEventCutsList);
  bailOutList->Add(eh.fEventHistogramsList);
  bailOutList->Add(pc.fParticleCutsList);
  bailOutList->Add(ph.fParticleHistogramsList);
  bailOutList->Add(qv.fQvectorList);
  bailOutList->Add(mupa.fCorrelationsList);
  bailOutList->Add(pw.fWeightsList);
  bailOutList->Add(cw.fCentralityWeightsList);
  bailOutList->Add(nl.fNestedLoopsList);
  bailOutList->Add(nua.fNUAList);
  bailOutList->Add(iv.fInternalValidationList);
  bailOutList->Add(t0.fTest0List);
  bailOutList->Add(es.fEtaSeparationsList);
  bailOutList->Add(res.fResultsList);

  // *) Add list with nested list to TDirectoryFile:
  dirFile->Add(bailOutList, true);
  dirFile->Write(dirFile->GetName(), TObject::kSingleKey + TObject::kOverwrite);
  delete dirFile;
  dirFile = NULL;
  f->Close();

  if (tc.fVerbose && !(tc.fSequentialBailout > 0)) { // then it will be called only once, for the only and permanent bailout
    ExitFunction(__FUNCTION__);
  }

  // *) Hasta la vista:
  if (finalBailout) {
    LOGF(fatal, "\033[1;31mHasta la vista - bailed out permanently in function %s at line %d\n The output file is: %s\n\n\033[0m", __FUNCTION__, __LINE__, sBailOutFile.Data());
  } else {
    LOGF(info, "\033[1;32mBailed out sequentially in function %s at line %d\n The output file is: %s\n\n\033[0m", __FUNCTION__, __LINE__, sBailOutFile.Data());
    if (tc.fVerbose) {
      ExitFunction(__FUNCTION__);
    }
  }

} // void BailOut(bool finalBailout = false)

//============================================================

void FillQvector(const double& dPhi, const double& dPt, const double& dEta)
{
  // Fill integrated Q-vector.
  // Example usage: this->FillQvector(dPhi, dPt, dEta);

  // TBI 20240430 I could optimize further, and have a bare version of this function when weights are NOT used.
  //              But since usage of weights amounts to checking a few simple booleans here, I do not anticipate any big gain in efficiency...

  if (tc.fVerboseForEachParticle) {
    StartFunction(__FUNCTION__);
    LOGF(info, "\033[1;32m dPhi = %f\033[0m", dPhi);
    LOGF(info, "\033[1;32m dPt  = %f\033[0m", dPt);
    LOGF(info, "\033[1;32m dEta = %f\033[0m", dEta);
  }

  // Particle weights:
  double wPhi = 1.;      // integrated phi weight
  double wPt = 1.;       // integrated pt weight
  double wEta = 1.;      // integrated eta weight
  double wToPowerP = 1.; // weight raised to power p

  if (pw.fUseWeights[wPHI]) {
    wPhi = Weight(dPhi, wPHI);
    if (!(wPhi > 0.)) {
      LOGF(error, "\033[1;33m%s wPhi is not positive\033[0m", __FUNCTION__);
      LOGF(fatal, "dPhi = %f\nwPhi = %f", dPhi, wPhi);
    }
  } // if(pw.fUseWeights[wPHI])

  if (pw.fUseWeights[wPT]) {
    wPt = Weight(dPt, wPT); // corresponding pt weight
    if (!(wPt > 0.)) {
      LOGF(error, "\033[1;33m%s wPt is not positive\033[0m", __FUNCTION__);
      LOGF(fatal, "dPt = %f\nwPt = %f", dPt, wPt);
    }
  } // if(pw.fUseWeights[wPT])

  if (pw.fUseWeights[wETA]) {
    wEta = Weight(dEta, wETA); // corresponding eta weight
    if (!(wEta > 0.)) {
      LOGF(error, "\033[1;33m%s wEta is not positive\033[0m", __FUNCTION__);
      LOGF(fatal, "dEta = %f\nwEta = %f", dEta, wEta);
    }
  } // if(pw.fUseWeights[wETA])

  if (qv.fCalculateQvectors) {
    for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
      for (int wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
        if (pw.fUseWeights[wPHI] || pw.fUseWeights[wPT] || pw.fUseWeights[wETA]) {
          wToPowerP = pow(wPhi * wPt * wEta, wp);
          qv.fQvector[h][wp] += TComplex(wToPowerP * TMath::Cos(h * dPhi), wToPowerP * TMath::Sin(h * dPhi)); // Q-vector with weights
        } else {
          qv.fQvector[h][wp] += TComplex(TMath::Cos(h * dPhi), TMath::Sin(h * dPhi)); // bare Q-vector without weights
        }
      } // for(int wp=0;wp<gMaxCorrelator+1;wp++)
    } // for(int h=0;h<gMaxHarmonic*gMaxCorrelator+1;h++)
  } // if (qv.fCalculateQvectors) {

  if (es.fCalculateEtaSeparations) { // yes, I can decouple this one from if (qv.fCalculateQvectors)
    if (dEta < 0.) {
      for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
        if (dEta < -1. * es.fEtaSeparationsValues[e] / 2.) { // yes, if eta separation is 0.2, then separation interval runs from -0.1 to 0.1
          qv.fMab[0][e] += wPhi * wPt * wEta;
          for (int h = 0; h < gMaxHarmonic; h++) {
            if (es.fEtaSeparationsSkipHarmonics[h]) {
              continue;
            }
            qv.fQabVector[0][h][e] += TComplex(wPhi * wPt * wEta * TMath::Cos((h + 1) * dPhi), wPhi * wPt * wEta * TMath::Sin((h + 1) * dPhi));
          }
        } // for (int h = 0; h < gMaxHarmonic; h++) {
      } // for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
    } else if (dEta > 0.) {
      for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
        if (dEta > es.fEtaSeparationsValues[e] / 2.) { // yes, if eta separation is 0.2, then separation interval runs from -0.1 to 0.1
          qv.fMab[1][e] += wPhi * wPt * wEta;
          for (int h = 0; h < gMaxHarmonic; h++) {
            {
              if (es.fEtaSeparationsSkipHarmonics[h]) {
                continue;
              }
              qv.fQabVector[1][h][e] += TComplex(wPhi * wPt * wEta * TMath::Cos((h + 1) * dPhi), wPhi * wPt * wEta * TMath::Sin((h + 1) * dPhi));
            }
          } // for (int h = 0; h < gMaxHarmonic; h++) {
        } // for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
      }
    }
  } // if(es.fCalculateEtaSeparations) {

  if (tc.fVerboseForEachParticle) {
    ExitFunction(__FUNCTION__);
  }

} // void FillQvector(const double& dPhi, const double& dPt, const double& dEta)

//============================================================

void FillQvectorFromSparse(const double& dPhi, const double& dPt, const double& dEta, const double& dCharge)
{
  // Fill integrated Q-vector using sparse histograms.

  // Remark: I pass by reference particle quantities, while event quantities (centrality, vertex z, ...) I fetch from data members (or from global variables in a macro).

  // To do:
  // 20250224 do I need to switch to this function also in InternalValidation()? I still use simple FillQvector() there.
  //          That would really make sense only after I add support for usage of particle weights in InternalValidation()

  if (tc.fVerboseForEachParticle) {
    StartFunction(__FUNCTION__);
    LOGF(info, "\033[1;32m dPhi    = %f\033[0m", dPhi);
    LOGF(info, "\033[1;32m dPt     = %f\033[0m", dPt);
    LOGF(info, "\033[1;32m dEta    = %f\033[0m", dEta);
    LOGF(info, "\033[1;32m dCharge = %f\033[0m", dCharge);
  }

  // Particle weights from sparse histograms:
  double wPhi = 1.;      // differential multidimensional phi weight, its dimensions are defined via enum eDiffPhiWeights
  double wPt = 1.;       // differential multidimensional pt weight, its dimensions are defined via enum eDiffPtWeights
  double wEta = 1.;      // differential multidimensional eta weight, its dimensions are defined via enum eDiffEtaWeights
  double wToPowerP = 1.; // weight raised to power p

  // *) Multidimensional phi weights:
  if (pw.fUseDiffPhiWeights[wPhiPhiAxis]) { // yes, 0th axis serves as a comon boolean for this category
    wPhi = WeightFromSparse(dPhi, dPt, dEta, dCharge, eDWPhi);
    // last argument is enum eDiffWeightCategory. Event quantities, e.g. centraliy and vz, I do not need to pass, because
    // for them I have ebye data members
    if (!(wPhi > 0.)) {
      LOGF(error, "\033[1;33m%s wPhi is not positive\033[0m", __FUNCTION__);
      LOGF(error, "dPhi = %f", dPhi);
      if (pw.fUseDiffPhiWeights[wPhiPtAxis]) {
        LOGF(fatal, "dPt = %f", dPt);
      }
      if (pw.fUseDiffPhiWeights[wPhiEtaAxis]) {
        LOGF(fatal, "dEta = %f", dEta);
      }
      if (pw.fUseDiffPhiWeights[wPhiChargeAxis]) {
        LOGF(fatal, "dCharge = %f", dCharge);
      }
      if (pw.fUseDiffPhiWeights[wPhiCentralityAxis]) {
        LOGF(fatal, "ebye.Centrality = %f", ebye.fCentrality);
      }
      if (pw.fUseDiffPhiWeights[wPhiVertex_zAxis]) {
        LOGF(fatal, "ebye.Vz = %f", ebye.fVz);
      }
      LOGF(fatal, "Multidimensional weight for enabled dimensions is wPhi = %f", wPhi);
    }
  } // if(pw.fUseDiffPhiWeights[wPhiPhiAxis])

  // *) Multidimensional pt weights:
  if (pw.fUseDiffPtWeights[wPtPtAxis]) {                     // yes, 0th axis serves as a comon boolean for this category
    wPt = WeightFromSparse(dPhi, dPt, dEta, dCharge, eDWPt); // TBI 20250224 not sure if this is the right/best approach
    // last argument is enum eDiffWeightCategory. Event quantities, e.g. centraliy and vz, I do not need to pass, because
    // for them I have ebye data members
    if (!(wPt > 0.)) {
      LOGF(error, "\033[1;33m%s wPt is not positive\033[0m", __FUNCTION__);
      LOGF(error, "dPt = %f", dPt);
      if (pw.fUseDiffPtWeights[wPtPtAxis]) {
        LOGF(fatal, "dPt = %f", dPt);
      }
      LOGF(fatal, "Multidimensional weight for enabled dimensions is wPt = %f", wPt);
    }
  } // if(pw.fUseDiffPtWeights[wPtPtAxis])

  // *) Multidimensional eta weights:
  if (pw.fUseDiffEtaWeights[wEtaEtaAxis]) {                    // yes, 0th axis serves as a comon boolean for this category
    wEta = WeightFromSparse(dPhi, dPt, dEta, dCharge, eDWEta); // TBI 20250224 not sure if this is the right/best approach
    // last argument is enum eDiffWeightCategory. Event quantities, e.g. centraliy and vz, I do not need to pass, because
    // for them I have ebye data members
    if (!(wEta > 0.)) {
      LOGF(error, "\033[1;33m%s wEta is not positive\033[0m", __FUNCTION__);
      LOGF(error, "dEta = %f", dEta);
      if (pw.fUseDiffEtaWeights[wEtaEtaAxis]) {
        LOGF(fatal, "dEta = %f", dEta);
      }
      LOGF(fatal, "Multidimensional weight for enabled dimensions is wEta = %f", wEta);
    }
  } // if(pw.fUseDiffEtaWeights[wEtaEtaAxis])

  if (qv.fCalculateQvectors) {
    for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
      for (int wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
        if (pw.fUseDiffPhiWeights[wPhiPhiAxis] || pw.fUseDiffPtWeights[wPtPtAxis] || pw.fUseDiffEtaWeights[wEtaEtaAxis]) {
          wToPowerP = pow(wPhi * wPt * wEta, wp);
          qv.fQvector[h][wp] += TComplex(wToPowerP * TMath::Cos(h * dPhi), wToPowerP * TMath::Sin(h * dPhi)); // Q-vector with weights
        } else {
          qv.fQvector[h][wp] += TComplex(TMath::Cos(h * dPhi), TMath::Sin(h * dPhi)); // bare Q-vector without weights
        }
      } // for(int wp=0;wp<gMaxCorrelator+1;wp++)
    } // for(int h=0;h<gMaxHarmonic*gMaxCorrelator+1;h++)
  } // if (qv.fCalculateQvectors) {

  if (es.fCalculateEtaSeparations) { // yes, I can decouple this one from if (qv.fCalculateQvectors)
    if (dEta < 0.) {
      for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
        if (dEta < -1. * es.fEtaSeparationsValues[e] / 2.) { // yes, if eta separation is 0.2, then separation interval runs from -0.1 to 0.1
          qv.fMab[0][e] += wPhi * wPt * wEta;
          for (int h = 0; h < gMaxHarmonic; h++) {
            if (es.fEtaSeparationsSkipHarmonics[h]) {
              continue;
            }
            qv.fQabVector[0][h][e] += TComplex(wPhi * wPt * wEta * TMath::Cos((h + 1) * dPhi), wPhi * wPt * wEta * TMath::Sin((h + 1) * dPhi));
          }
        } // for (int h = 0; h < gMaxHarmonic; h++) {
      } // for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
    } else if (dEta > 0.) {
      for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
        if (dEta > es.fEtaSeparationsValues[e] / 2.) { // yes, if eta separation is 0.2, then separation interval runs from -0.1 to 0.1
          qv.fMab[1][e] += wPhi * wPt * wEta;
          for (int h = 0; h < gMaxHarmonic; h++) {
            {
              if (es.fEtaSeparationsSkipHarmonics[h]) {
                continue;
              }
              qv.fQabVector[1][h][e] += TComplex(wPhi * wPt * wEta * TMath::Cos((h + 1) * dPhi), wPhi * wPt * wEta * TMath::Sin((h + 1) * dPhi));
            }
          } // for (int h = 0; h < gMaxHarmonic; h++) {
        } // for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
      }
    }
  } // if(es.fCalculateEtaSeparations) {

  if (tc.fVerboseForEachParticle) {
    ExitFunction(__FUNCTION__);
  }

} // void FillQvectorFromSparse(const double& dPhi, const double& dPt, const double& dEta, const double& dCharge)

//============================================================

void Fillqvector(const double& dPhi, const double& kineVarValue, eqvectorKine kineVarChoice, const double& dEta = 0.)
{
  // Fill differential q-vector, in generic kinematic variable. Here "kine" originally meant vs. pt or vs. eta, now it's general.
  // Example usage #1: this->Fillqvector(dPhi, dPt, PTq); // differential q-vectors without using eta separations
  // Example usage #2: this->Fillqvector(dPhi, dPt, PTq, dEta); // differential q-vectors with using eta separations (I need dEta of particle to decide whether particle is added to qa or qb)

  if (tc.fVerboseForEachParticle) {
    StartFunction(__FUNCTION__);
  }

  // *) Mapping between enum's "eqvectorKine" on one side, and "eAsFunctionOf", "eWeights" and "eDiffWeights" on the other:
  //    TBI 20240212 I could promote this also to a member function, if I need it elsewhere. Or I could use TExMap?
  eAsFunctionOf AFO_var = eAsFunctionOf_N;      // this local variable determines the enum "eAsFunctionOf" which corresponds to enum "eqvectorKine"
  eWeights AFO_weight = eWeights_N;             // this local variable determines the enum "eWeights" which corresponds to enum "eqvectorKine"
  eDiffWeights AFO_diffWeight = eDiffWeights_N; // this local variable determines the enum "eDiffWeights" which corresponds to enum "eqvectorKine"
  switch (kineVarChoice) {
    case PTq:
      AFO_var = AFO_PT;
      AFO_weight = wPT;
      AFO_diffWeight = wPHIPT;
      break;
    case ETAq:
      AFO_var = AFO_ETA;
      AFO_weight = wETA;
      AFO_diffWeight = wPHIETA;
      break;
    default:
      LOGF(fatal, "\033[1;31m%s at line %d : this kineVarChoice = %d is not supported yet. \033[0m", __FUNCTION__, __LINE__, static_cast<int>(kineVarChoice));
      break;
  } // switch(kineVarChoice)

  // *) Insanity checks on above settings:
  if (AFO_var == eAsFunctionOf_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : AFO_var == eAsFunctionOf_N => add some more entries to the case statement \033[0m", __FUNCTION__, __LINE__);
  }
  if (AFO_weight == eWeights_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : AFO_weight == eWeights_N => add some more entries to the case statement \033[0m", __FUNCTION__, __LINE__);
  }
  if (AFO_diffWeight == eDiffWeights_N) {
    LOGF(fatal, "\033[1;31m%s at line %d : AFO_diffWeight == eDiffWeights_N => add some more entries to the case statement \033[0m", __FUNCTION__, __LINE__);
  }

  // *) Get the desired bin number:
  int bin = -1;
  if (res.fResultsPro[AFO_var]) {
    bin = res.fResultsPro[AFO_var]->FindBin(kineVarValue);         // this 'bin' starts from 1, i.e. this is genuine histogram bin
    if (0 >= bin || res.fResultsPro[AFO_var]->GetNbinsX() < bin) { // either underflow or overflow is hit, meaning that histogram is booked in narrower range than cuts
      LOGF(fatal, "\033[1;31m%s at line %d : kineVarChoice = %d, bin = %d, kineVarValue = %f \033[0m", __FUNCTION__, __LINE__, static_cast<int>(kineVarChoice), bin, kineVarValue);
    }
  }

  // *) Get all integrated kinematic weights:
  double wToPowerP = 1.;     // weight raised to power p
  double kineVarWeight = 1.; // e.g. this can be integrated pT or eta weight
  if (pw.fUseWeights[AFO_weight]) {
    kineVarWeight = Weight(kineVarValue, AFO_weight); // corresponding e.g. pt or eta weight
    if (!(kineVarWeight > 0.)) {
      LOGF(fatal, "\033[1;31m%s at line %d : kineVarWeight is not positive \033[0m", __FUNCTION__, __LINE__);
      // TBI 20240212 or could I just skip this particle?
    }
  } // if(fUseWeights[AFO_weight]) {

  // *) Get all differential phi-weights for this kinematic variable:
  //    Remark: special treatment is justified for phi-weights, because q-vector is defined in terms of phi-weights.
  double diffPhiWeightsForThisKineVar = 1.;
  if (pw.fUseDiffWeights[AFO_diffWeight]) {
    diffPhiWeightsForThisKineVar = DiffWeight(dPhi, kineVarValue, kineVarChoice); // corresponding differential phi weight as a function of e.g. pt or eta
    if (!(diffPhiWeightsForThisKineVar > 0.)) {
      LOGF(fatal, "\033[1;31m%s at line %d : diffPhiWeightsForThisKineVar is not positive \033[0m", __FUNCTION__, __LINE__);
      // TBI 20240212 or could I just skip this particle?
    }
  } // if(pw.fUseDiffWeights[AFO_diffWeight]) {

  // *) Finally, fill differential q-vector in that bin:
  for (int h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
    for (int wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
      if (pw.fUseWeights[AFO_weight] || pw.fUseDiffWeights[AFO_diffWeight]) {
        // TBI 20240212 supported at the moment: e.g. q-vector vs pt can be weighted only with diff. phi(pt) and integrated pt weights.
        // It cannot be weighted in addition with eta weights, since in any case I anticipate I will do always 1-D analysis, by integrating out all other dependencies
        wToPowerP = pow(diffPhiWeightsForThisKineVar * kineVarWeight, wp);
        qv.fqvector[kineVarChoice][bin - 1][h][wp] += TComplex(wToPowerP * TMath::Cos(h * dPhi), wToPowerP * TMath::Sin(h * dPhi)); // q-vector with weights
      } else {
        qv.fqvector[kineVarChoice][bin - 1][h][wp] += TComplex(TMath::Cos(h * dPhi), TMath::Sin(h * dPhi)); // bare q-vector without weights
      }
    } // for(int wp=0;wp<gMaxCorrelator+1;wp++)
  } // for(int h=0;h<gMaxHarmonic*gMaxCorrelator+1;h++)

  // *) Differential nested loops:
  if (nl.fCalculateKineCustomNestedLoops) {
    nl.ftaNestedLoopsKine[kineVarChoice][bin - 1][0]->AddAt(dPhi, qv.fqVectorEntries[kineVarChoice][bin - 1]);
    nl.ftaNestedLoopsKine[kineVarChoice][bin - 1][1]->AddAt(diffPhiWeightsForThisKineVar * kineVarWeight, qv.fqVectorEntries[kineVarChoice][bin - 1]);
  }

  // *) Multiplicity counter in this bin:
  qv.fqVectorEntries[kineVarChoice][bin - 1]++; // count number of particles in this pt bin in this event

  // *) Usage of eta separations in differential correlations:
  if (es.fCalculateEtaSeparations && es.fCalculateEtaSeparationsAsFunctionOf[AFO_var]) { // yes, I can decouple this one from if (qv.fCalculateQvectors)

    if (AFO_var == AFO_ETA) {
      LOGF(fatal, "\033[1;31m%s at line %d : AFO_var == AFO_ETA . This doesn't make any sense in this context. \033[0m", __FUNCTION__, __LINE__);
    }

    if (dEta < 0.) {
      for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
        if (dEta < -1. * es.fEtaSeparationsValues[e] / 2.) {                      // yes, if eta separation is 0.2, then separation interval runs from -0.1 to 0.1
          qv.fmab[0][bin - 1][e] += diffPhiWeightsForThisKineVar * kineVarWeight; // Remark: I can hardwire linear weight like this only for 2-p correlation
          for (int h = 0; h < gMaxHarmonic; h++) {
            if (es.fEtaSeparationsSkipHarmonics[h]) {
              continue;
            }
            qv.fqabVector[0][bin - 1][h][e] += TComplex(diffPhiWeightsForThisKineVar * kineVarWeight * TMath::Cos((h + 1) * dPhi), diffPhiWeightsForThisKineVar * kineVarWeight * TMath::Sin((h + 1) * dPhi)); // Remark: I can hardwire linear weight like this only for 2-p correlation
          }
        } // for (int h = 0; h < gMaxHarmonic; h++) {
      } // for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
    } else if (dEta > 0.) {
      for (int e = 0; e < gMaxNumberEtaSeparations; e++) {
        if (dEta > es.fEtaSeparationsValues[e] / 2.) {                            // yes, if eta separation is 0.2, then separation interval runs from -0.1 to 0.1
          qv.fmab[1][bin - 1][e] += diffPhiWeightsForThisKineVar * kineVarWeight; // Remark: I can hardwire linear weight like this only for 2-p correlation
          for (int h = 0; h < gMaxHarmonic; h++) {
            {
              if (es.fEtaSeparationsSkipHarmonics[h]) {
                continue;
              }
              qv.fqabVector[1][bin - 1][h][e] += TComplex(diffPhiWeightsForThisKineVar * kineVarWeight * TMath::Cos((h + 1) * dPhi), diffPhiWeightsForThisKineVar * kineVarWeight * TMath::Sin((h + 1) * dPhi)); // Remark: I can hardwire linear weight like this only for 2-p correlation
            }
          } // for (int h = 0; h < gMaxHarmonic; h++) {
        } // for (int e = 0; e < gMaxNumberEtaSeparations; e++) { // eta separation
      }
    }
  } // if(es.fCalculateEtaSeparations) {

  if (tc.fVerboseForEachParticle) {
    ExitFunction(__FUNCTION__);
  }

} // void Fillqvector(const double& dPhi, const double& kineVarValue, eqvectorKine kineVarChoice)

//============================================================

void CalculateEverything()
{
  // Calculate everything for selected events and particles.
  // Remark: Data members for Q-vectors, containers for nested loops, etc., must all be filled when this function is called.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Calculate multiparticle correlations (standard, isotropic, same harmonic):
  if (qv.fCalculateQvectors && mupa.fCalculateCorrelations) {
    this->CalculateCorrelations();
  }

  // *) Calculate differential ("kine") multiparticle correlations:
  //    Remark: vs. pt, vs. eta, etc., are all calculated here
  if (qv.fCalculateQvectors && mupa.fCalculateCorrelationsAsFunctionOf[AFO_PT]) {
    this->CalculateKineCorrelations(AFO_PT);
  }
  if (qv.fCalculateQvectors && mupa.fCalculateCorrelationsAsFunctionOf[AFO_ETA]) {
    this->CalculateKineCorrelations(AFO_ETA);
  }

  // *) Calculate Test0: TBI 20240110 name convention
  //    Remark: integrated, vs. M and vs. centrality are all calculated here
  if (qv.fCalculateQvectors && t0.fCalculateTest0) {
    this->CalculateTest0();
  }

  // *) Calculate kine Test0: TBI 20240110 name convention
  //    Remark: vs. pt, vs. eta, etc., are all calculated here
  if (qv.fCalculateQvectors && t0.fCalculateTest0AsFunctionOf[AFO_PT]) {
    this->CalculateKineTest0(AFO_PT);
  }
  if (qv.fCalculateQvectors && t0.fCalculateTest0AsFunctionOf[AFO_ETA]) {
    this->CalculateKineTest0(AFO_ETA);
  }

  // *) Calculate nested loops:
  if (nl.fCalculateNestedLoops) {
    this->CalculateNestedLoops();
    if (mupa.fCalculateCorrelations) {
      // I do not have option here for Test0, because in Test0 I cross-check either e-by-e with CustomNestedLoops or
      // for all events  with IV + fRescaleWithTheoreticalInput = true
      this->ComparisonNestedLoopsVsCorrelations(); // I call it here, so comparison is performed cumulatively after each event. The final printout corresponds to all events.
    }
  }

  // *) Calculate correlations with eta separations:
  if (es.fCalculateEtaSeparations) {
    this->CalculateEtaSeparations();
    if (es.fCalculateEtaSeparationsAsFunctionOf[AFO_PT]) {
      this->CalculateKineEtaSeparations(AFO_PT); // The implementation of CalculateKineEtaSeparations( ... ) is generic and can be used for any other "kine" variable, for which it makes sense
    }
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // void CalculateEverything()

//============================================================

template <eRecSim rs, typename T>
void MainLoopOverParticles(T const& tracks)
{
  // This is the main loop over particles, in which Q-vectors (both integrated and differential) and particle histograms are filled, particle cuts applied, etc.

  // Remark #1:
  // *) To process only reconstructed Run 3, use processRec(...), i.e. set field "processRec": "true" in json config
  // *) To process both reconstructed and simulated Run 3, use processRecSim(...), i.e. set field "processRecSim": "true" in json config
  // *) To process only simulated Run 3, use processSim(...), i.e. set field "processSim": "true" in json config

  // Remark #2:
  // *) To process Run 1 and Run 2 converted data, use in the same spirit e.g. processRec_Run2(...), i.e. set field "processRec_Run2": "true" in json config

  // Remark #3:
  // *) There is also processTest(...), to process data with minimum subscription to the tables. To use it, set field "processTest": "true" in json config

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Declare local kinematic variables:
  double dPhi = 0.; // azimuthal angle
  double dPt = 0.;  // transverse momentum
  double dEta = 0.; // pseudorapidity

  // *) If random access of tracks from collection is requested, use Fisher-Yates algorithm to generate random indices:
  if (tc.fUseFisherYates) {
    if (tc.fRandomIndices) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    this->RandomIndices(tracks.size());
    if (!tc.fRandomIndices) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
  }

  // *) Local timestamp:
  if (tc.fUseStopwatch && tc.fVerboseUtility) {
    LOGF(info, "  Local timer starts at line %d", __LINE__);
    tc.fTimer[eLocal]->Reset();
    tc.fTimer[eLocal]->Start();
  }

  // *) Main loop over particles:
  // for (auto& track : tracks) { // default standard way of looping of tracks
  auto track = tracks.iteratorAt(0); // set the type and scope from one instance
  for (int64_t i = 0; i < tracks.size(); i++) {

    // *) Access track sequentially from collection of tracks (default), or randomly using Fisher-Yates algorithm:
    if (!tc.fUseFisherYates) {
      track = tracks.iteratorAt(i);
    } else {
      track = tracks.iteratorAt(static_cast<int64_t>(tc.fRandomIndices->GetAt(i)));
    }

    // *) Skip track objects which are not valid tracks (e.g. Run 2 and 1 tracklets, etc.):
    if (!ValidTrack<rs>(track)) {
      continue;
    }

    // *) Fill particle histograms before particle cuts:
    if (ph.fFillParticleHistograms || ph.fFillParticleHistograms2D || qa.fFillQAParticleHistograms2D) {
      FillParticleHistograms<rs>(track, eBefore);
    }

    // *) Particle cuts counters (use only during QA, as this is computationally heavy):
    if (pc.fUseParticleCutCounterAbsolute || pc.fUseParticleCutCounterSequential) {
      ParticleCutsCounters<rs>(track);
    }

    // *) Particle cuts:
    if (!ParticleCuts<rs>(track, eCut)) { // Main call for event cuts.
      continue;                           // not return!!
    }

    // *) Fill particle histograms after particle cuts:
    if (ph.fFillParticleHistograms || ph.fFillParticleHistograms2D || qa.fFillQAParticleHistograms2D) {
      FillParticleHistograms<rs>(track, eAfter);
    }

    // *) Intitialize local kinematic variables:
    //    Remark: for "eRecSim" processing, kinematics is taken from "reconstructed".
    dPhi = track.phi();
    dPt = track.pt();
    dEta = track.eta();

    // Remark: Keep in sync all calls and flags below with the ones in InternalValidation().
    // *) Integrated Q-vectors:
    if (qv.fCalculateQvectors || es.fCalculateEtaSeparations) {
      if (!(pw.fUseDiffPhiWeights[wPhiPhiAxis] || pw.fUseDiffPtWeights[wPtPtAxis] || pw.fUseDiffPtWeights[wEtaEtaAxis])) {
        // legacy integrated weights:
        this->FillQvector(dPhi, dPt, dEta); // all 3 arguments are passed by reference
      } else {
        // this is now the new approach, with sparse histograms:
        this->FillQvectorFromSparse(dPhi, dPt, dEta, track.sign()); // particle arguments are passed by reference.
                                                                    // Event observables (centrality, vertex z, ...), I do not need to pass as arguments,
                                                                    // as I have data members for them (ebye.fCentrality, ebye.Vz, ...)
      }
    }

    // *) Differential q-vectors:
    // **) pt-dependence:
    if (qv.fCalculateQvectors && (mupa.fCalculateCorrelationsAsFunctionOf[AFO_PT] || t0.fCalculateTest0AsFunctionOf[AFO_PT]) && !es.fCalculateEtaSeparations) {
      // In this branch I do not need eta separation, so the lighter call can be executed:
      this->Fillqvector(dPhi, dPt, PTq); // first 2 arguments are passed by reference, 3rd argument is enum
    } else if (es.fCalculateEtaSeparations && es.fCalculateEtaSeparationsAsFunctionOf[AFO_PT]) {
      // In this branch I do need eta separation, so the heavier call must be executed:
      // Remark: Within Fillqvector() I check again all the relevant flags.
      this->Fillqvector(dPhi, dPt, PTq, dEta); // first 2 arguments and the last one are passed by reference, 3rd argument is enum. "kine" variable is the 2nd argument
    }
    // **) eta-dependence:
    if (qv.fCalculateQvectors && (mupa.fCalculateCorrelationsAsFunctionOf[AFO_ETA] || t0.fCalculateTest0AsFunctionOf[AFO_ETA])) {
      // Remark: For eta dependence I do not consider es.fCalculateEtaSeparations, because in this context that calculation is meaningless.
      this->Fillqvector(dPhi, dEta, ETAq); // first 2 arguments are passed by reference, 3rd argument is enum
    }

    // *) Fill nested loops containers:
    if (nl.fCalculateNestedLoops || nl.fCalculateCustomNestedLoops) {
      this->FillNestedLoopsContainers(ebye.fSelectedTracks, dPhi, dPt, dEta); // all 4 arguments are passed by reference
    }

    // *) Counter of selected tracks in the current event:
    ebye.fSelectedTracks++;
    if (ebye.fSelectedTracks >= ec.fdEventCuts[eMultiplicity][eMax]) {
      break;
    }

    // *) Break the loop if fixed number of particles is taken randomly from each event (use always in combination with tc.fUseFisherYates = true):
    if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && tc.fFixedNumberOfRandomlySelectedTracks == ebye.fSelectedTracks) {
      LOGF(info, "%s : Breaking the loop over particles, since requested fixed number of %d particles was reached", __FUNCTION__, tc.fFixedNumberOfRandomlySelectedTracks);
      break;
    }

  } // for (auto& track : tracks)

  // *) Local timestamp:
  if (tc.fUseStopwatch && tc.fVerboseUtility) {
    LOGF(info, "  Local timer ends at line %d, time elapsed ... %.6f", __LINE__, tc.fTimer[eLocal]->RealTime());
    tc.fTimer[eLocal]->Continue();
  }

  // *) Insanity check on fixed number of randomly selected tracks:
  if (tc.fFixedNumberOfRandomlySelectedTracks > 0 && tc.fFixedNumberOfRandomlySelectedTracks < ebye.fSelectedTracks) {
    LOGF(fatal, "\033[1;31mIn this event there are too few particles (ebye.fSelectedTracks = %d), and requested number of fixed number randomly selected tracks %d couldn't be reached\033[0m", ebye.fSelectedTracks, tc.fFixedNumberOfRandomlySelectedTracks);
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }
} // template <eRecSim rs, typename T> void MainLoopOverParticles(T const& tracks) {

//============================================================

template <eRecSim rs, typename T1, typename T2, typename T3>
void Steer(T1 const& collision, T2 const& bcs, T3 const& tracks)
{
  // This is the only function to be called in processRec(...), processRecSim(...), and processSim(...).
  // All analysis workflow is defined step-by-step here, via dedicated function calls.
  // The order of function calls obviously matters.

  if (tc.fVerbose) {
    StartFunction(__FUNCTION__);
  }

  // *) Dry run:
  if (tc.fDryRun) {
    EventCounterForDryRun(eFill);
    EventCounterForDryRun(ePrint);
    Preprocess<rs>(collision, bcs); // yes, so that e.g. I can only test if the particle and centrality weights were correctly fetched from external file and initialized locally into data members
    return;
  }

  // *) Reset event-by-event quantities: TBI 20240430 I do not need this call also here really, but it doesn't hurt either...
  ResetEventByEventQuantities();

  // *) Only do internal validation for all implemented correlators against the theoretical values:
  if (iv.fUseInternalValidation) {
    InternalValidation();
    return;
  }

  // *) Global timestamp:
  if (tc.fUseStopwatch) {
    LOGF(info, "\033[1;32m=> Global timer: Steer begins ... %.6f\033[0m", tc.fTimer[eGlobal]->RealTime());
    tc.fTimer[eGlobal]->Continue(); // yes
  }

  // *) Do all thingies before starting to process data from this collision (e.g. cut on number of events (both total and selected), fetch the run number, etc.):
  Preprocess<rs>(collision, bcs);

  // *) Determine collision reference multiplicity:
  DetermineReferenceMultiplicity<rs>(collision);

  // *) Determine collision centrality:
  DetermineCentrality<rs>(collision);

  // *) Determine collision occupancy:
  DetermineOccupancy<rs>(collision);

  // *) Determine collision interaction rate and current run duration:
  DetermineInteractionRateAndCurrentRunDuration<rs>(collision, bcs);

  // *) Determine vertex z position:
  DetermineVertexZ<rs>(collision);

  // *) Fill event histograms before event cuts:
  if (eh.fFillEventHistograms || qa.fFillQAEventHistograms2D || qa.fFillQAParticleEventHistograms2D) {
    // Remark: I do not above the flag fFillQACorrelationsVsHistograms2D, because as a part of QA I calculate <2> only after cuts in any case
    FillEventHistograms<rs>(collision, tracks, eBefore);
  }

  // *) Print info on the current event number (total, before cuts):
  if (tc.fVerboseEventCounter) {
    PrintEventCounter(eBefore);
  }

  // *) Event cuts counters (use only during QA, as this is computationally heavy):
  if (ec.fUseEventCutCounterAbsolute || ec.fUseEventCutCounterSequential) {
    EventCutsCounters<rs>(collision, tracks);
  }

  // *) Event cuts:
  if (!EventCuts<rs>(collision, tracks, eCut)) { // Main call for event cuts
    return;
  }

  // *) Main loop over particles:
  MainLoopOverParticles<rs>(tracks);

  // *) Determine multiplicity of this event, for all "vs. mult" results:
  DetermineMultiplicity();

  // *) Remaining event cuts which can be applied only after the loop over particles is performed:
  if (!RemainingEventCuts()) {
    // yes, I need to remove particles from ParticleHistograms, which were filled in the MainLoopOverParticles also for events which didn't survive RemainingEventCuts
    BanishmentLoopOverParticles<rs>(tracks);
    ResetEventByEventQuantities();
    return;
  }

  // *) Fill event histograms after event AND particle cuts:
  if (eh.fFillEventHistograms || qa.fFillQAEventHistograms2D || qa.fFillQAParticleEventHistograms2D || qa.fFillQACorrelationsVsHistograms2D) {
    FillEventHistograms<rs>(collision, tracks, eAfter);
  }

  // *) Fill subevent multiplicities:
  //    Remark: I can call this one only after Qa and Qb vectors are filled, and after all particle and event cuts:
  if (es.fCalculateEtaSeparations) {
    FillSubeventMultiplicities<rs>();
  }

  // *) Calculate everything for selected events and particles:
  CalculateEverything();

  // *) Reset event-by-event quantities:
  ResetEventByEventQuantities();

  // *) QA:
  if (qa.fCheckUnderflowAndOverflow) { // TBI 20240507 introduce eventualy common function QA(), within which I will call all specific QA functions
    CheckUnderflowAndOverflow();
  }

  // *) Print info on the current event number after cuts:
  if (tc.fVerboseEventCounter) {
    PrintEventCounter(eAfter);
  }

  // *) Per request, print content of event cut counters:
  if (ec.fPrintCutCounterContent) {
    PrintCutCounterContent();
  }

  // *) Global timestamp:
  if (tc.fUseStopwatch) {
    LOGF(info, "\033[1;32m=> Global timer: Steer ends ... %.6f\033[0m\n", tc.fTimer[eGlobal]->RealTime());
    tc.fTimer[eGlobal]->Continue(); // yes
  }

  if (tc.fVerbose) {
    ExitFunction(__FUNCTION__);
  }

} // template <eRecSim rs, typename T1, typename T2> void Steer(T1 const* collision, T2 const* tracks)

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_MEMBERFUNCTIONS_H_
