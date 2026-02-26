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

/// \file MuPa-Enums.h
/// \brief ... TBI 20250425
/// \author Ante.Bilandzic@cern.ch

#ifndef PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_
#define PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_

enum eConfiguration {
  eTaskIsConfiguredFromJson = 1, // here I start from 1 exceptionally, because these enums are used as bin contents, and ROOT starts counting bins from 1
  eTaskName,
  eRunNumber,
  eDryRun,
  eVerbose,
  eVerboseUtility,
  eVerboseForEachParticle,
  eVerboseEventCounter,
  ePlainPrintout,
  eDoAdditionalInsanityChecks,
  eInsanityCheckForEachParticle,
  eWhichProcess,
  eRandomSeed,
  eUseFisherYates,
  eFixedNumberOfRandomlySelectedTracks,
  eUseStopwatch,
  eFloatingPointPrecision,
  eSequentialBailout,
  eUseSpecificCuts,
  eWhichSpecificCuts,
  eSkipTheseRuns,
  eUseSetBinLabel, // Use or not ->SetBinLabel(...)
  eUseClone,       // Use or not ->Clone()
  eUseFormula,     // Use or not class TFormula
  eUseDatabasePDG, // Use or not class TDatabasePDG
  eConfiguration_N
};

enum eProcess {
  eProcessRec = 0,     // Run 3, only reconstructed
  eProcessRecSim,      // Run 3, both reconstructed and simulated
  eProcessSim,         // Run 3, only simulated
  eProcessRec_Run2,    // Run 2, only reconstructed
  eProcessRecSim_Run2, // Run 2, both reconstructed and simulated
  eProcessSim_Run2,    // Run 2, only simulated
  eProcessRec_Run1,    // Run 1, only reconstructed
  eProcessRecSim_Run1, // Run 1, both reconstructed and simulated
  eProcessSim_Run1,    // Run 1, only simulated
  eProcessTest,        // minimum subscription to the tables, for testing purposes
  eProcessQA,          // maximum subscription to the tables, for QA purposes. Basically: eProcessRec + otherwise unnecessary tables
  eProcessHepMChi,     // special subscription when I extract info from the table HepMCHeavyIons TBI 20250429 merge this case eventually with RecSim cases
  // Generic flags, calculated and set from individual flags above in DefaultConfiguration(), AFTER process switch was taken into account:
  eGenericRec,    // generic "Rec" case, eTest is treated for the time being as "Rec". eQA is also in this category
  eGenericRecSim, // generic "RecSim" case
  eGenericSim,    // generic "Sim" case
  eProcess_N
};

enum eRecSim { eRec = 0,
               eSim,
               eRecAndSim,
               eRec_Run2, // converted Run 2 data
               eSim_Run2,
               eRecAndSim_Run2,
               eRec_Run1, // converted Run 1 data
               eSim_Run1,
               eRecAndSim_Run1,
               eTest, // remember that as of 20250315 I can use "cfWhichSpecificCuts": "Test" in JSON, to configure quickly all cuts for this case
               eQA };

enum eBeforeAfter { eBefore = 0, // use this one for cuts
                    eAfter = 1 };

enum eMinMax { eMin = 0,
               eMax = 1 };

enum eXYZ { eX = 0,
            eY = 1,
            eZ = 2 };

enum eDefaultColors { eColor = kBlack,
                      eFillColor = kGray };

enum eWeights { wPHI = 0,
                wPT = 1,
                wETA = 2,
                wCHARGE = 3,
                eWeights_N };

enum eDiffWeights { // TBI 20250215 this is now obsolete, superseeded with more general implementation, see enums eDiffWeightCategory, eDiffPhiWeights, etc.
  wPHIPT = 0,
  wPHIETA,
  wPHICHARGE,
  eDiffWeights_N
};

enum eDiffWeightCategory {
  eDWPhi = 0, // corresponds to eDiffPhiWeights structure, here the fundamental 0-th axis never to be projected out is "phi"
  eDWPt,      // corresponds to eDiffPtWeights structure, here the fundamental 0-th axis never to be projected out is "pt"
  eDWEta,     // corresponds to eDiffEtaWeights structure, here the fundamental 0-th axis never to be projected out is "eta"
  // ...
  eDiffWeightCategory_N
};

enum eDiffPhiWeights {
  wPhiPhiAxis = 0, // this is the main axis in this category, the only axis which shall never be projected out. If I project out all remaining axes, I shall recover the standard integrated phi weights
  wPhiPtAxis,
  wPhiEtaAxis,
  wPhiChargeAxis,
  wPhiCentralityAxis,
  wPhiVertexZAxis,
  eDiffPhiWeights_N
};

enum eDiffPtWeights {
  wPtPtAxis = 0,
  wPtChargeAxis,
  wPtCentralityAxis,
  eDiffPtWeights_N
};

enum eDiffEtaWeights {
  wEtaEtaAxis = 0,
  wEtaChargeAxis,
  wEtaCentralityAxis,
  eDiffEtaWeights_N
};

enum eVnPsin { eVn = 0,
               ePsin = 1 };

enum eEventHistograms {
  eNumberOfEvents = 0,    // Total events = eNumberOfEvents + eBefore, Selected events = eNumberOfEvents + eAfter
  eTotalMultiplicity,     // TBI 20241123 I define it as tracks.size(), but most likely I do not need this
  eMultiplicity,          // see documentation for ebye.fMultiplicity
  eReferenceMultiplicity, // see documentation for ebye.fReferenceMultiplicity
  eCentrality,            // default centrality estimator
  eVertexX,
  eVertexY,
  eVertexZ,
  eNContributors, // number of tracks used for the vertex
  eImpactParameter,
  eEventPlaneAngle,
  eOccupancy,             // from helper task o2-analysis-event-selection, see also IA's presentation in https://indico.cern.ch/event/1464946, slide 38. Use specific occupancy estimator via eOccupancyEstimator
  eInteractionRate,       // from utility ctpRateFetcher .
  eCurrentRunDuration,    // calculated with utility ctpRateFetcher
  eMultMCNParticlesEta08, // from helper task table o2::aod::MultMCExtras
  eEventHistograms_N
};

enum eEventCuts {
  // a) For available event selection bits, check https://github.com/AliceO2Group/O2Physics/blob/master/Common/CCDB/EventSelectionParams.cxx (and .h for better documentation)
  // b) Some settings are configurable, check: https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/eventSelection.cxx
  eTrigger = eEventHistograms_N,   // Implemented and validated so far:
                                   // a) Run 3: "kTVXinTRD" (use optionally for systematics, and only in real data)
                                   // b) Run 2: "kINT7" (at the moment there is a warning if not used in real data.). Not validated in Monte Carlo.
                                   // c) Run 1: TBI 20241209 check if I can use kINT7 also for Run 1
  eSel7,                           // See def. of sel7 in Ref. b) above. Event selection decision based on V0A & V0C => use only in Run 2 and Run 1.
                                   // TBI 20250115 it removes 99% of events in MC LHC21i6a, check this further
  eSel8,                           // See def. of sel8 in Ref. b) above. Event selection decision based on TVX => use only in Run 3, both for data and MC
                                   // *) As of 20240410, kNoITSROFrameBorder (only in MC) and kNoTimeFrameBorder event selection cuts are part of Sel8
                                   //    See also email from EK from 2024041
  eMultiplicityEstimator,          // see documentation for ebye.fMultiplicity
  eReferenceMultiplicityEstimator, // see documentation for ebye.fReferenceMultiplicity
  eCentralityEstimator,            // the default centrality estimator, set via configurable. All supported centrality estimators, for QA, etc, are in enum eCentralityEstimators
  eSelectedEvents,                 // selected events = eNumberOfEvents + eAfter => therefore I do not need a special histogram for it
  eNoSameBunchPileup,              // reject collisions in case of pileup with another collision in the same foundBC (emails from IA on 20240404 and EK on 20240410)
  eIsGoodZvtxFT0vsPV,              // small difference between z-vertex from PV and from FT0 (emails from IA on 20240404 and EK on 20240410)
                                   // Avoid using kIsGoodZvtxFT0vsPV selection bit for Pb-Pb 2024 apass1, see IA email from 20250115.
                                   // Therefore, until further notice, use this one in LHC23zzh, but not in LHC24ar and LHC24as
  eIsVertexITSTPC,                 // at least one ITS-TPC track (reject vertices built from ITS-only tracks) (emails from IA on 20240404 and EK on 20240410
  eIsVertexTOFmatched,             // at least one of vertex contributors is matched to TOF
  eIsVertexTRDmatched,             // at least one of vertex contributors is matched to TRD
  eNoCollInTimeRangeStrict,        // rejects a collision if there are other events in dtime +/- 10 μs, see IA Slide 39 in https://indico.cern.ch/event/1462154/
                                   // 20250122 Per feedback from IA, use this one only as a part of systematic check, and use eNoCollInTimeRangeStandard by default
  eNoCollInTimeRangeStandard,      // rejects a collision if there are other events in dtime +/- 2 μs + additional cuts on multiplicity, see IA Slide 39 in https://indico.cern.ch/event/1462154/
  eNoCollInRofStrict,              // rejects a collision if there are other events within the same ROF (in-ROF pileup), ROF = "ITS Readout Frames",
                                   // see IA Slide 39 in https://indico.cern.ch/event/1462154/
                                   // 20250122 Per feedback from IA, use this one only as a part of systematic check, and use eNoCollInRofStandard by default
  eNoCollInRofStandard,            // same as previous + additional cuts on multiplicity, see IA Slide 39 in https://indico.cern.ch/event/1462154/
  eNoHighMultCollInPrevRof,        // veto an event if FT0C amplitude in previous ITS ROF is above threshold (default is >5000 a.e. by FT0C), see IA Slide 39 in https://indico.cern.ch/event/1462154/
                                   // 20250122 Per feedback from IA, use it only in 2023 PbPb data (e.g. eLHC23zzh), in 2024 PbPb data this one has no effect (do not use in eLHC24ar and eLHC24as)
  eIsGoodITSLayer3,                // number of inactive chips on ITS layer 3 is below maximum allowed value
  eIsGoodITSLayer0123,             // numbers of inactive chips on ITS layers 0-3 are below maximum allowed values
  eIsGoodITSLayersAll,             // numbers of inactive chips on all ITS layers are below maximum allowed values
  eOccupancyEstimator,             // the default Occupancy estimator, set via configurable. All supported centrality estimators, for QA, etc, are in enum eOccupancyEstimators
  eMinVertexDistanceFromIP,        // if sqrt(vx^2+vy^2+vz^2) < MinVertexDistanceFromIP, the event is rejected. This way, I remove suspicious events with |vertex| = 0.
  eNoPileupTPC,                    // no pileup in TPC
  eNoPileupFromSPD,                // no pileup according to SPD vertexer
  eNoSPDOnVsOfPileup,              // no out-of-bunch pileup according to online-vs-offline SPD correlation
  eRefMultVsNContrUp,              // formula for upper boundary cut in eReferenceMultiplicity_vs_NContributors (remember that I use naming convention "x_vs_y")
  eRefMultVsNContrLow,             // formula for lower boundary cut in eReferenceMultiplicity_vs_NContributors (remember that I use naming convention "x_vs_y")
  eCentralityCorrelationsCut,      // port of void SetCentralityCorrelationsCuts(...) from MuPa class. Example format: "CentFT0C_CentFT0M", so IFS is "_", until proven otherwise

  // RCT flags, see https://indico.cern.ch/event/1545907/ + up-to-date code in Common/CCDB/RCTSelectionFlags.h
  // Remark 1: For the time being, I support here differentially 6 flags used to define the combined "CBT" flag, see if (label == "CBT") in Common/CCDB/RCTSelectionFlags.h
  // Remark 2: If I want to use directly the combined "CBT" flag, see how it can be done using RCTFlagsChecker in
  //           https://github.com/AliceO2Group/O2Physics/blob/master/DPG/Tasks/AOTEvent/timeDependentQa.cxx#L115
  //           But check before the memory status after RCTFlagsChecker is used.
  eFT0Bad,
  eITSBad,
  eITSLimAccMCRepr,
  eTPCBadTracking,
  eTPCLimAccMCRepr,
  eTPCBadPID,
  // ...
  eCentralityWeights, // used for centrality flattening. Remember that this event cut must be implemented very last,
                      // therefore I have it separately implemented for Run 3,2,1 in EventCuts() at the very end in each case.
                      // Use only for small non-uniformity in centrality distribution (e.g. of the biggest dip in distribution is up to 20% compared to uniform part of cent. distribution),
                      // otherwise this flattening is too costly in terms of statistics.
  eEventCuts_N
};

enum eEventCutsFormulas { // special treatment for all event cuts defined via mathematical formula, because for them I have to do one additional layer of booking using TFormula
  eRefMultVsNContrUp_Formula = 0,
  eRefMultVsNContrLow_Formula,
  eEventCutsFormulas_N
};

enum eParticleHistograms {

  // from o2::aod::Tracks (Track parameters at their point closest to the collision vertex)
  ePhi = 0,
  ePt,
  eEta,
  eCharge, // Charge: positive: 1, negative: -1

  // from o2::aod::TracksExtra_001 - I keep the ordering here the same as in the TracksExtra_001 table
  etpcNClsFindable,
  etpcNClsShared,
  eitsChi2NCl, // TBI 20250110 I see for this one [478682:track-selection]: [15:35:00][INFO] Track selection, set max chi2 per cluster ITS: 36
               //              But even with open particle cuts, this distribution doesn't cross 30... There is a sudden drop round 22, but when I apply other cuts
               //              that tail is gone already.
  etpcNClsFound,
  etpcNClsCrossedRows,
  eitsNCls,
  eitsNClsInnerBarrel,
  etpcCrossedRowsOverFindableCls,
  etpcFoundOverFindableCls, // TBI 20250110 I keep this one in sync with values for etpcCrossedRowsOverFindableCls
  etpcFractionSharedCls,
  etpcChi2NCl, // TBI 20250110 this one shall resemble aodTrack->GetTPCchi2()/aodTrack->GetTPCNcls(), but cross-check with the experts. Particles with tpcChi2NCl > 4. I reject now by default.
               //              See what I documented in AliPhysics below // task->SetParticleCuts("TPCChi2perNDF",4.,-44); // VAL
               // 20250123 in some Run 2 analysis, 2.5 was used as a default. Check that value as a part of systematics
  // from o2::aod::TracksDCA
  edcaXY,
  edcaZ,

  // the rest:
  ePDG,

  // counter:
  eParticleHistograms_N
};

enum eParticleHistograms2D { // All 2D histograms are first implemented in eQAParticleHistograms2D, the ones which I need regularly, are then promoted to this category.
  ePhiPt = 0,
  ePhiEta,
  eParticleHistograms2D_N
};

enum eParticleCuts {

  // from o2::aod::TrackSelection (https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html#o2-analysis-trackselection)
  // See also O2Physics/Common/DataModel/TrackSelectionTables.h
  etrackCutFlag = eParticleHistograms_N, // General selection, with centrally tuned particle cuts for tpcNClsFound, itsNCls, etc.
                                         // As of 20250113, this cut still has not effect, neither in Run 3 nor in converted Run 2. Use instead trackCutFlagFb1 and/or trackCutFlagFb2 below.
  etrackCutFlagFb1,                      // Global tracks in Run 3. Closest possible match to global track definition in Run 2, which are selected with eisGlobalTrack.
                                         // For the definition, see:
                                         //  a) "filtbit1" in https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/trackselection.cxx#L128
                                         //  b) "getGlobalTrackSelectionRun3ITSMatch" in https://github.com/AliceO2Group/O2Physics/blob/master/Common/Core/TrackSelectionDefaults.cxx#L43
                                         // When I use this flag, make sure I do NOT cut on something on which this cut is already cutting by default (e.g. pt-dependent DCA xy cut)
  etrackCutFlagFb2,                      // Global tracks in Run 3, similar as etrackCutFlagFb1, but more stringent (since 2 points in ITS are required in inner barrel (IB)).
                                         // Unlike etrackCutFlagFb1 (1 ITS point is required), it produces a 20% dip in azimuthal acceptance for 1.2 < phi < 1.6, in LHC24ar/559545
                                         // DCAxy and z are significantly further depleted, when compared to etrackCutFlagFb1
  eisQualityTrack,                       // Do not use in Run 3, but it can be used in Run 2 and Run 1
                                         // In Run 2, it is already requested in isGlobalTrack, through definition kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks
                                         // See O2Physics/Common/DataModel/TrackSelectionTables.h for further details.
                                         // Therefore, vary in Run 2 only is isGlobalTrack is NOT requested by default.
  eisPrimaryTrack,                       // Validated in Run 3. See also isPVContributor
                                         // In Run 2, it is already requested in isGlobalTrack, through definition kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks
                                         // See O2Physics/Common/DataModel/TrackSelectionTables.h for further details.
                                         // Therefore, vary in Run 2 only is isGlobalTrack is NOT requested by default.
  eisInAcceptanceTrack,                  // kInAcceptanceTracks = kPtRange | kEtaRange . Pt is open, and |eta| < 0.8.
                                         // But after I already cut directly on 0.2 < pt < 5.0 and |eta| < 0.8, it has no effect.
                                         // Can be used both in Run 3 and Run 2.
                                         // TBI 20250113 remove this cut eventually from the code, because I cut direcly on  0.2 < pt < 5.0 and |eta| < 0.8 in any case.
  eisGlobalTrack,                        // Do not use in Run 3, it can be used directly only in Run 2 and Run 1, see definition in:
                                         //   https://github.com/AliceO2Group/O2Physics/blob/master/Common/Core/TrackSelectionDefaults.cxx#L23
                                         // For Run 3 global tracks, I need to use TrackSelection getGlobalTrackSelectionRun3ITSMatch(int matching, int passFlag) from
                                         //   https://github.com/AliceO2Group/O2Physics/blob/master/Common/Core/TrackSelectionDefaults.cxx#L43
                                         // That is precisely definition of filtBit1 in https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/trackselection.cxx
                                         // So: etrackCutFlagFb1 in Run 3 is a closest match to eisGlobalTrack in Run 2
  eisPVContributor,                      // Run 3: Has this track contributed to the collision vertex fit
                                         // Tracks used in vertex fit are flagged as "contributors" (-> track.isPVContributor() for AO2D tracks). Such a track is
                                         // allowed to contribute to only one PV. => See further details in RS presentation https://indico.cern.ch/event/1453901/timetable/#12-track-reconstruction
                                         // This cut affects significantly distributions of other tracking parameters. Most notably, after using this cut, DCAz distribution is reduced to ~ 1mm range.
                                         // But for global tracks in any case we request very stringent DCA cut.
                                         // pt and eta distributions are only mildly affected.
                                         // Do not use in Run 2 and Run 1.
                                         // It's not the same as isPrimaryTrack cut, albeit there is an overlap.
  // special treatment:
  ePtDependentDCAxyParameterization,

  eParticleCuts_N
};

enum eAsFunctionOf { // this is a specific enum only for 1D dependence
  // 1D:
  AFO_INTEGRATED = 0,
  AFO_MULTIPLICITY, // vs. default multiplicity, which is (at the moment) fSelectedTracks, i.e. number of tracks in Q-vector
  AFO_CENTRALITY,   // vs. default centrality estimator, see how it's calculated in DetermineCentrality(...)
  AFO_PT,
  AFO_ETA,
  AFO_OCCUPANCY,          // vs. default "occupancy" variable which is (at the moment) "FT0COccupancyInTimeRange" (alternative is "TrackOccupancyInTimeRange")
  AFO_INTERACTIONRATE,    // vs. "interation rate"
  AFO_CURRENTRUNDURATION, // vs. "current run duration", i.e. vs "seconds since start of run"
  AFO_VZ,                 // vs. "vertex z position"
  AFO_CHARGE,             // vs. "particle charge"
  // ...
  eAsFunctionOf_N
}; // prefix is needed, to avoid conflict with enum eKinematics

enum eAsFunctionOf2D { // this is a specific enum only for 2D dependence
  // 2D:
  AFO_CENTRALITY_PT = 0,
  AFO_CENTRALITY_ETA,
  AFO_CENTRALITY_CHARGE,
  AFO_CENTRALITY_VZ,
  AFO_PT_ETA,
  AFO_PT_CHARGE,
  AFO_ETA_CHARGE,
  // ...
  eAsFunctionOf2D_N
};

enum eAsFunctionOf3D { // this is a specific enum only for 3D dependence
  // 3D:
  AFO_CENTRALITY_PT_ETA = 0,
  AFO_CENTRALITY_PT_CHARGE,
  AFO_CENTRALITY_PT_VZ,
  AFO_CENTRALITY_ETA_VZ,
  AFO_CENTRALITY_ETA_CHARGE,
  AFO_CENTRALITY_VZ_CHARGE,
  AFO_PT_ETA_CHARGE,
  // ...
  eAsFunctionOf3D_N
};

enum eNUAPDF {
  ePhiNUAPDF = 0,
  ePtNUAPDF,
  eEtaNUAPDF,
  eNUAPDF_N
};

enum eqvectorKine { // Here "kine" originally meant "kinematic", i.e. vs. pt or vs. eta, now it's general.
  // 1D:
  PTq = 0,
  ETAq,
  CHARGEq,
  // ...

  // 2D: // Yes, I linearize 2D case, in an analogy with "global bin" structure for multidimensional histograms.
  PT_ETAq,
  PT_CHARGEq,
  ETA_CHARGEq,
  // ...

  // 3D: // Yes, I linearize 3D case, in an analogy with "global bin" structure for multidimensional histograms.
  PT_ETA_CHARGEq,
  // ...
  eqvectorKine_N
};

enum eTimer {
  eGlobal = 0,
  eLocal,
  eTimer_N
};

enum eEventCounterForDryRun {
  eFill = 0,
  ePrint
};

enum eCutModus {
  eCut = 0,           // standard, i.e. no cut counters are used
  eCutCounterBinning, // dry call to EventCuts and ParticleCuts, just to establish order of binning in CutCountets, which resembles order of cut implementation
  eCutCounterAbsolute,
  eCutCounterSequential
};

enum eCutCounter {
  eAbsolute = 0,
  eSequential,
  eCutCounter_N
};

enum eQAEventHistograms2D {
  // General (estimators can be chosen via configurables):
  eMultiplicity_vs_ReferenceMultiplicity = 0, // multiplicity is x, reference multiplicity is y. I can swap offline if needed: histOriginal->GetBinContent(x,y); histSwapped->Fill(y,x);
  eMultiplicity_vs_NContributors,
  eMultiplicity_vs_Centrality,
  eMultiplicity_vs_VertexZ,
  eMultiplicity_vs_Occupancy,
  eMultiplicity_vs_InteractionRate, // TBI 20250331 I ctd. below with more histos in category eMultiplicity_vs_... - re-organize at some point bookkeping here
  eReferenceMultiplicity_vs_NContributors,
  eReferenceMultiplicity_vs_Centrality,
  eReferenceMultiplicity_vs_VertexZ,
  eReferenceMultiplicity_vs_Occupancy,
  eReferenceMultiplicity_vs_InteractionRate,
  eNContributors_vs_Centrality,
  eNContributors_vs_VertexZ,
  eNContributors_vs_Occupancy,
  eNContributors_vs_InteractionRate,
  eCentrality_vs_VertexZ,
  eCentrality_vs_Occupancy,
  eCentrality_vs_ImpactParameter, // [sim] = reconstructed centrality vs. simulated impact parameter. [rec] = ... TBI 20241210
  eCentrality_vs_InteractionRate,
  eVertexZ_vs_Occupancy,
  eVertexZ_vs_InteractionRate,
  // ...
  // Specific (everything is hardwired):
  eMultNTracksPV_vs_MultNTracksGlobal,  // Run 3 multiplicity
  eCentFT0C_vs_CentFT0CVariant1,        // Run 3 centrality
  eCentFT0C_vs_CentFT0M,                // Run 3 centrality
  eCentFT0C_vs_CentFV0A,                // Run 3 centrality
  eCentFT0C_vs_CentNTPV,                // Run 3 centrality
  eCentFT0C_vs_CentNGlobal,             // Run 3 centrality
  eCentFT0M_vs_CentNTPV,                // Run 3 centrality
  eCentRun2V0M_vs_CentRun2SPDTracklets, // Run 2 centrality (do not use in Run 1 converted, because there is no centrality information)
  eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange,
  eCurrentRunDuration_vs_InteractionRate, // ...
  // ...
  // Unsorted category: TBI 20250331 not sure if I will keep these ones permanently:
  eMultiplicity_vs_FT0CAmplitudeOnFoundBC,
  eCentFT0C_vs_FT0CAmplitudeOnFoundBC,
  eCentrality_vs_CentralitySim, // correlation between centrality determined from reconstructed data, and centrality determined at generated level (from IP). I save it as [eSim], not as [eRec]
  // ...
  eQAEventHistograms2D_N
};

enum eQAParticleHistograms2D {
  ePt_vs_dcaXY,
  eQAParticleHistograms2D_N
};

enum eQAParticleEventHistograms2D {
  // In this category I do correlation <some-particle-property> vs. some-event-property.
  // The < ... > goes over all particles in that event.
  // All < ... > over particles are calculated with helper TProfile fQAParticleEventProEbyE
  // For instance: <nITScls> vs. current run duration
  eCurrentRunDuration_vs_itsNClsEbyE,
  eCurrentRunDuration_vs_itsNClsNegEtaEbyE,
  eCurrentRunDuration_vs_itsNClsPosEtaEbyE,
  eCurrentRunDuration_vs_Eta0804EbyE,
  eCurrentRunDuration_vs_Eta0400EbyE,
  eCurrentRunDuration_vs_Eta0004EbyE,
  eCurrentRunDuration_vs_Eta0408EbyE,
  eCurrentRunDuration_vs_Pt0005EbyE,
  eCurrentRunDuration_vs_Pt0510EbyE,
  eCurrentRunDuration_vs_Pt1050EbyE,
  eQAParticleEventHistograms2D_N
};

enum eQAParticleEventProEbyE {
  eitsNClsEbyE = 1,                   // Labels average <itsNCls> in a given event (therefore "EbyE" is appended). Yes, from one, because it runs over bin content and entries in TProfile for most of the time.
  eitsNClsNegEtaEbyE,                 // <itsNCls> in a given event for eta < 0
  eitsNClsPosEtaEbyE,                 // <itsNCls> in a given event for eta > 0
  eEta0804EbyE,                       // <eta> in a given event for -0.8 < eta < -0.4
  eEta0400EbyE,                       // <eta> in a given event for -0.4 < eta <  0.0
  eEta0004EbyE,                       // <eta> in a given event for  0.0 < eta <  0.4
  eEta0408EbyE,                       // <eta> in a given event for  0.4 < eta <  0.8
  ePt0005EbyE,                        // <pt> in a given event for  0.0 < pt < 0.5
  ePt0510EbyE,                        // <pt> in a given event for  0.5 < pt < 1.0
  ePt1050EbyE,                        // <pt> in a given event for  1.0 < pt < 5.0
  eMeanPhi,                           // <phi> in an event TBI 20250214 I need to unify naming convention for <> with previous enums in above in the series, but okay...
  eMeanPt,                            // <pt> in an event
  eMeanEta,                           // <eta> in an event
  eMeanCharge,                        // <charge> in an event
  eMeantpcNClsFindable,               // <tpcNClsFindable> in an event
  eMeantpcNClsShared,                 // <tpcNClsShared> in an event
  eMeanitsChi2NCl,                    // <itsChi2NCl> in an event
  eMeantpcNClsFound,                  // <tpcNClsFound> in an event
  eMeantpcNClsCrossedRows,            // <tpcNClsCrossedRows> in an event
  eMeanitsNCls,                       // <itsNCls> in an event
  eMeanitsNClsInnerBarrel,            // <itsNClsInnerBarrel> in an event
  eMeantpcCrossedRowsOverFindableCls, // <tpcCrossedRowsOverFindableCls> in an event
  eMeantpcFoundOverFindableCls,       // <tpcFoundOverFindableCls> in an event
  eMeantpcFractionSharedCls,          // <tpcFractionSharedCls> in an event
  eMeantpcChi2NCl,                    // <tpcChi2NCl> in an event
  eMeandcaXY,                         // <dcaXY> in an event
  eMeandcaZ,                          // <dcaZ> in an event
  eQAParticleEventProEbyE_N
};

enum eQACorrelationsVsHistograms2D {
  // In this category I correlate <2> vs. some-event-property.
  // For instance: <2> vs. ref. mult
  //               <2> vs. <pT>, where <pT> is calculated from all particles in that event (so in this sense, it's an event property as well)
  // Remark 1: If I would ever need the same thingie for <4>, <6>, etc., just introduce new dimension in 2D histogram
  // Remark 2: All < ... > over particles are calculated with helper TProfile fQAParticleEventProEbyE
  eCorrelations_vs_Multiplicity = 0,
  eCorrelations_vs_ReferenceMultiplicity,
  eCorrelations_vs_Centrality,
  // ...
  eCorrelations_vs_MeanPhi,
  eCorrelations_vs_MeanPt,
  eCorrelations_vs_MeanEta,
  eCorrelations_vs_MeanCharge,
  eCorrelations_vs_MeantpcNClsFindable,
  eCorrelations_vs_MeantpcNClsShared,
  eCorrelations_vs_MeanitsChi2NCl,
  eCorrelations_vs_MeantpcNClsFound,
  eCorrelations_vs_MeantpcNClsCrossedRows,
  eCorrelations_vs_MeanitsNCls,
  eCorrelations_vs_MeanitsNClsInnerBarrel,
  eCorrelations_vs_MeantpcCrossedRowsOverFindableCls,
  eCorrelations_vs_MeantpcFoundOverFindableCls,
  eCorrelations_vs_MeantpcFractionSharedCls,
  eCorrelations_vs_MeantpcChi2NCl,
  eCorrelations_vs_MeandcaXY,
  eCorrelations_vs_MeandcaZ,
  // ...
  eQACorrelationsVsHistograms2D_N
};

enum eQACorrelationsVsInteractionRateVsProfiles2D_N {
  // In this category I fill <2> in 2D profile spanned by IR vs. some-other-observable (IR is always x axis)
  // For instance: <2> is filled in TProfile2D spanned by IR vs. CurrentRunDuration (crd)
  //               <2> is filled in TProfile2D spanned by IR vs. <pT>, where <pT> is calculated from all particles in that event
  // Remark 1: If I would ever need the same thingie for <4>, <6>, etc., just introduce new dimension in 2D profile
  // Remark 2: All < ... > over particles are calculated with helper TProfile fQAParticleEventProEbyE
  eCorrelationsVsInteractionRate_vs_CurrentRunDuration = 0,
  eCorrelationsVsInteractionRate_vs_Multiplicity,
  eCorrelationsVsInteractionRate_vs_ReferenceMultiplicity,
  eCorrelationsVsInteractionRate_vs_Centrality,
  // ...
  eCorrelationsVsInteractionRate_vs_MeanPhi,
  eCorrelationsVsInteractionRate_vs_SigmaMeanPhi,
  eCorrelationsVsInteractionRate_vs_MeanPt,
  eCorrelationsVsInteractionRate_vs_SigmaMeanPt,
  eCorrelationsVsInteractionRate_vs_MeanEta,
  eCorrelationsVsInteractionRate_vs_SigmaMeanEta,
  // ...
  eQACorrelationsVsInteractionRateVsProfiles2D_N
};

enum eReferenceMultiplicityEstimators {
  // Run 3:
  eMultTPC = 0,
  eMultFV0M,          // ref. mult from helper task o2-analysis-multiplicity-table
  eMultFT0C,          // ref. mult from helper task o2-analysis-multiplicity-table
  eMultFT0M,          // ref. mult from helper task o2-analysis-multiplicity-table
  eMultNTracksPV,     // ref. mult from helper task o2-analysis-multiplicity-table
  eMultNTracksGlobal, // ref. mult from helper task o2-analysis-multiplicity-table
  // Run 2:
  eMultTracklets, // ref. mult from helper task o2-analysis-multiplicity-table, use only for Run 2
  eReferenceMultiplicityEstimators_N
};

enum eCentralityEstimators {
  // Run 3:
  eCentFT0C = 0,
  eCentFT0CVariant1,
  eCentFT0M,
  eCentFV0A,
  eCentNTPV,
  eCentNGlobal,
  // Run 2:
  eCentRun2V0M,
  eCentRun2SPDTracklets,
  eCentralityEstimators_N
};

enum eOccupancyEstimators {
  eTrackOccupancyInTimeRange, // from helper task o2-analysis-event-selection, see also IA's presentation in https://indico.cern.ch/event/1464946, slide 38
  eFT0COccupancyInTimeRange,  // from helper task o2-analysis-event-selection
  eOccupancyEstimators_N
};

enum eEventCounter {
  eTotal,     // total number of events, before any cuts are applied
  eProcessed, // number of processed events, i.e. number of events which survived cuts and on which analysis have been performed
  eEventCounter_N
};

enum eSpecificCuts {
  // Run 3:
  eLHC23zzh,
  eLHC24ar,
  eLHC24as,
  // Run 2:
  eLHC15o,
  // Run 1:
  // ...
  // Cuts for minimal subscription, "processTest": "true in JSON
  eTestCuts,
  eSpecificCuts_N
};

enum eRunTime {
  eStartOfRun = 0, // in abs. seconds since Unix epoch
  eEndOfRun,       // in abs. seconds since Unix epoch
  eDurationInSec,  // in seconds
  eRunTime_N
};

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_
