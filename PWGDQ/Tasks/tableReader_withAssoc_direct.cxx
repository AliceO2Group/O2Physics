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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//   Configurable workflow for running several DQ or other PWG analyses

#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/PID/PIDTOFParamService.h"
#include "Common/Core/TableHelper.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonUtils/ConfigurableParam.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <DetectorsVertexing/PVertexer.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH2.h>
#include <THashList.h>
#include <TList.h>
#include <TMathBase.h>
#include <TString.h>

#include <RtypesCore.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::common::core;

Zorro zorro;

// Some definitions
namespace o2::aod
{
namespace dqanalysisflags
{
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);                                     //! Hash used in event mixing
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 32);                     //! Event decision
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelected, isBarrelSelected, 32);                   //! Barrel track decisions
DECLARE_SOA_COLUMN(BarrelAmbiguityInBunch, barrelAmbiguityInBunch, int8_t);          //! Barrel track in-bunch ambiguity
DECLARE_SOA_COLUMN(BarrelAmbiguityOutOfBunch, barrelAmbiguityOutOfBunch, int8_t);    //! Barrel track out of bunch ambiguity
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, 32); //! Barrel prefilter decisions

// Bcandidate columns
DECLARE_SOA_COLUMN(RunNumber, runNumber, uint64_t);
DECLARE_SOA_COLUMN(EventIdx, eventIdx, uint64_t);
DECLARE_SOA_COLUMN(EventTimestamp, eventTimestamp, uint64_t);
DECLARE_SOA_COLUMN(massBcandidate, MBcandidate, float);
DECLARE_SOA_COLUMN(MassDileptonCandidate, massDileptonCandidate, float);
DECLARE_SOA_COLUMN(deltaMassBcandidate, deltaMBcandidate, float);
DECLARE_SOA_COLUMN(pTBcandidate, PtBcandidate, float);
DECLARE_SOA_COLUMN(EtaBcandidate, etaBcandidate, float);
DECLARE_SOA_COLUMN(PhiBcandidate, phiBcandidate, float);
DECLARE_SOA_COLUMN(RapBcandidate, rapBcandidate, float);
DECLARE_SOA_COLUMN(LxyBcandidate, lxyBcandidate, float);
DECLARE_SOA_COLUMN(LxyBcandidateErr, lxyBcandidateErr, float);
DECLARE_SOA_COLUMN(LxyzBcandidate, lxyzBcandidate, float);
DECLARE_SOA_COLUMN(LxyzBcandidateErr, lxyzBcandidateErr, float);
DECLARE_SOA_COLUMN(LzBcandidate, lzBcandidate, float);
DECLARE_SOA_COLUMN(LzBcandidateErr, lzBcandidateErr, float);
DECLARE_SOA_COLUMN(TauxyBcandidate, tauxyBcandidate, float);
DECLARE_SOA_COLUMN(TauxyBcandidateErr, tauxyBcandidateErr, float);
DECLARE_SOA_COLUMN(TauzBcandidate, tauzBcandidate, float);
DECLARE_SOA_COLUMN(TauzBcandidateErr, tauzBcandidateErr, float);
DECLARE_SOA_COLUMN(CosPBcandidate, cosPBcandidate, float);
DECLARE_SOA_COLUMN(Chi2Bcandidate, chi2Bcandidate, float);
DECLARE_SOA_COLUMN(GlobalIndexassoc, globalIndexassoc, uint64_t);
DECLARE_SOA_COLUMN(GlobalIndexleg1, globalIndexleg1, uint64_t);
DECLARE_SOA_COLUMN(GlobalIndexleg2, globalIndexleg2, uint64_t);
DECLARE_SOA_COLUMN(Ptassoc, ptassoc, float);
DECLARE_SOA_COLUMN(PINassoc, pINassoc, float);
DECLARE_SOA_COLUMN(Etaassoc, etaassoc, float);
DECLARE_SOA_COLUMN(Phiassoc, phiassoc, float);
DECLARE_SOA_COLUMN(Ptpair, ptpair, float);
DECLARE_SOA_COLUMN(Etapair, etapair, float);
DECLARE_SOA_COLUMN(Ptleg1, ptleg1, float);
DECLARE_SOA_COLUMN(PINleg1, pINleg1, float);
DECLARE_SOA_COLUMN(Etaleg1, etaleg1, float);
DECLARE_SOA_COLUMN(Phileg1, phileg1, float);
DECLARE_SOA_COLUMN(Ptleg2, ptleg2, float);
DECLARE_SOA_COLUMN(PINleg2, pINleg2, float);
DECLARE_SOA_COLUMN(Etaleg2, etaleg2, float);
DECLARE_SOA_COLUMN(Phileg2, phileg2, float);
DECLARE_SOA_COLUMN(TPCnsigmaKaassoc, tpcnsigmaKaassoc, float);
DECLARE_SOA_COLUMN(TPCnsigmaPiassoc, tpcnsigmaPiassoc, float);
DECLARE_SOA_COLUMN(TPCnsigmaPrassoc, tpcnsigmaPrassoc, float);
DECLARE_SOA_COLUMN(TOFnsigmaKaassoc, tofnsigmaKaassoc, float);
DECLARE_SOA_COLUMN(TPCnsigmaElleg1, tpcnsigmaElleg1, float);
DECLARE_SOA_COLUMN(TPCnsigmaPileg1, tpcnsigmaPileg1, float);
DECLARE_SOA_COLUMN(TPCnsigmaPrleg1, tpcnsigmaPrleg1, float);
DECLARE_SOA_COLUMN(TPCnsigmaElleg2, tpcnsigmaElleg2, float);
DECLARE_SOA_COLUMN(TPCnsigmaPileg2, tpcnsigmaPileg2, float);
DECLARE_SOA_COLUMN(TPCnsigmaPrleg2, tpcnsigmaPrleg2, float);
DECLARE_SOA_COLUMN(ITSClusterMapassoc, itsClusterMapassoc, uint8_t);
DECLARE_SOA_COLUMN(ITSClusterMapleg1, itsClusterMapleg1, uint8_t);
DECLARE_SOA_COLUMN(ITSClusterMapleg2, itsClusterMapleg2, uint8_t);
DECLARE_SOA_COLUMN(ITSChi2assoc, itsChi2assoc, float);
DECLARE_SOA_COLUMN(ITSChi2leg1, itsChi2leg1, float);
DECLARE_SOA_COLUMN(ITSChi2leg2, itsChi2leg2, float);
DECLARE_SOA_COLUMN(TPCNclsassoc, tpcNclsassoc, float);
DECLARE_SOA_COLUMN(TPCNclsleg1, tpcNclsleg1, float);
DECLARE_SOA_COLUMN(TPCNclsleg2, tpcNclsleg2, float);
DECLARE_SOA_COLUMN(TPCChi2assoc, tpcChi2assoc, float);
DECLARE_SOA_COLUMN(TPCChi2leg1, tpcChi2leg1, float);
DECLARE_SOA_COLUMN(TPCChi2leg2, tpcChi2leg2, float);
DECLARE_SOA_BITMAP_COLUMN(IsJpsiFromBSelected, isJpsiFromBSelected, 32);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);

DECLARE_SOA_COLUMN(Massee, massee, float);
DECLARE_SOA_COLUMN(Etaee, etaee, float);
DECLARE_SOA_COLUMN(Rapee, rapee, float);
DECLARE_SOA_COLUMN(Phiee, phiee, float);
DECLARE_SOA_COLUMN(Ptee, ptee, float);
DECLARE_SOA_COLUMN(Lxyee, lxyee, float);
DECLARE_SOA_COLUMN(LxyeePoleMass, lxyeepolemass, float);
DECLARE_SOA_COLUMN(Lzee, lzee, float);
DECLARE_SOA_COLUMN(LxyeePoleMassPVrecomputed, lxyeePoleMassPVrecomputed, float);
DECLARE_SOA_COLUMN(MultiplicityFT0A, multiplicityFT0AJPsi2ee, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0CJPsi2ee, float);
DECLARE_SOA_COLUMN(PercentileFT0M, percentileFT0MJPsi2ee, float);
DECLARE_SOA_COLUMN(MultiplicityNContrib, multiplicityNContribJPsi2ee, float);
DECLARE_SOA_COLUMN(AmbiguousInBunchPairs, AmbiguousJpsiPairsInBunch, bool);
DECLARE_SOA_COLUMN(AmbiguousOutOfBunchPairs, AmbiguousJpsiPairsOutOfBunch, bool);
DECLARE_SOA_COLUMN(Corrassoc, corrassoc, bool);
DECLARE_SOA_BITMAP_COLUMN(IsMuonSelected, isMuonSelected, 32);                //! Muon track decisions (joinable to FwdTrackAssoc)
DECLARE_SOA_COLUMN(MuonAmbiguityInBunch, muonAmbiguityInBunch, int8_t);       //! Muon track in-bunch ambiguity
DECLARE_SOA_COLUMN(MuonAmbiguityOutOfBunch, muonAmbiguityOutOfBunch, int8_t); //! Muon track out of bunch ambiguity
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASHA", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected);
DECLARE_SOA_TABLE(BarrelAmbiguities, "AOD", "DQBARRELAMB", dqanalysisflags::BarrelAmbiguityInBunch, dqanalysisflags::BarrelAmbiguityOutOfBunch);
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsBarrelSelectedPrefilter);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);                                               //! joinable to FwdTrackAssoc
DECLARE_SOA_TABLE(MuonAmbiguities, "AOD", "DQMUONAMB", dqanalysisflags::MuonAmbiguityInBunch, dqanalysisflags::MuonAmbiguityOutOfBunch); //! joinable to FwdTracks

DECLARE_SOA_TABLE(JPsieeCandidates, "AOD", "DQPSEUDOPROPER",
                  dqanalysisflags::Massee, dqanalysisflags::Ptee, dqanalysisflags::Etaee, dqanalysisflags::Rapee,
                  dqanalysisflags::Phiee, dqanalysisflags::Lxyee, dqanalysisflags::LxyeePoleMass, dqanalysisflags::Lzee, dqanalysisflags::LxyeePoleMassPVrecomputed,
                  dqanalysisflags::AmbiguousInBunchPairs, dqanalysisflags::AmbiguousOutOfBunchPairs,
                  dqanalysisflags::MultiplicityFT0A, dqanalysisflags::MultiplicityFT0C, dqanalysisflags::PercentileFT0M, dqanalysisflags::MultiplicityNContrib);

DECLARE_SOA_TABLE(BmesonCandidates, "AOD", "DQBMESONS",
                  dqanalysisflags::RunNumber, dqanalysisflags::EventIdx, dqanalysisflags::EventTimestamp,
                  dqanalysisflags::massBcandidate, dqanalysisflags::MassDileptonCandidate, dqanalysisflags::deltaMassBcandidate,
                  dqanalysisflags::pTBcandidate, dqanalysisflags::EtaBcandidate, dqanalysisflags::PhiBcandidate, dqanalysisflags::RapBcandidate,
                  dqanalysisflags::LxyBcandidate, dqanalysisflags::LxyBcandidateErr, dqanalysisflags::LxyzBcandidate, dqanalysisflags::LxyzBcandidateErr,
                  dqanalysisflags::LzBcandidate, dqanalysisflags::LzBcandidateErr, dqanalysisflags::TauxyBcandidate, dqanalysisflags::TauxyBcandidateErr,
                  dqanalysisflags::TauzBcandidate, dqanalysisflags::TauzBcandidateErr, dqanalysisflags::CosPBcandidate, dqanalysisflags::Chi2Bcandidate,
                  dqanalysisflags::GlobalIndexassoc, dqanalysisflags::GlobalIndexleg1, dqanalysisflags::GlobalIndexleg2,
                  dqanalysisflags::PINassoc, dqanalysisflags::Etaassoc, dqanalysisflags::Ptpair, dqanalysisflags::Etapair,
                  dqanalysisflags::PINleg1, dqanalysisflags::Etaleg1, dqanalysisflags::PINleg2, dqanalysisflags::Etaleg2,
                  dqanalysisflags::TPCnsigmaKaassoc, dqanalysisflags::TPCnsigmaPiassoc, dqanalysisflags::TPCnsigmaPrassoc, dqanalysisflags::TOFnsigmaKaassoc,
                  dqanalysisflags::TPCnsigmaElleg1, dqanalysisflags::TPCnsigmaPileg1, dqanalysisflags::TPCnsigmaPrleg1,
                  dqanalysisflags::TPCnsigmaElleg2, dqanalysisflags::TPCnsigmaPileg2, dqanalysisflags::TPCnsigmaPrleg2,
                  dqanalysisflags::ITSClusterMapassoc, dqanalysisflags::ITSClusterMapleg1, dqanalysisflags::ITSClusterMapleg2,
                  dqanalysisflags::ITSChi2assoc, dqanalysisflags::ITSChi2leg1, dqanalysisflags::ITSChi2leg2,
                  dqanalysisflags::TPCNclsassoc, dqanalysisflags::TPCNclsleg1, dqanalysisflags::TPCNclsleg2,
                  dqanalysisflags::TPCChi2assoc, dqanalysisflags::TPCChi2leg1, dqanalysisflags::TPCChi2leg2,
                  dqanalysisflags::IsJpsiFromBSelected, dqanalysisflags::IsBarrelSelected);
} // namespace o2::aod

// Using definitions (data-only)
using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra>;
using MyEventsSelected = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::EventCuts, aod::MixingHashes>;
using MyEventsWithDqFilter = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::DQEventFilter>;

using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA,
                                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksWithCovNoTOF = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA,
                                             aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                             aod::pidTPCFullKa, aod::pidTPCFullPr>;
using MyBarrelTracksWithCovWithAmbiguities = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA,
                                                       aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                                       aod::pidTPCFullKa, aod::pidTPCFullPr, aod::BarrelAmbiguities>;

using MyMuonTracksWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksDCA, aod::FwdTrackCovFwd>;
using MyMuonTracksWithCovWithAmbiguities = soa::Join<aod::FwdTracks, aod::FwdTracksDCA, aod::FwdTrackCovFwd, aod::MuonAmbiguities>;

using MyDielectronCandidates = soa::Join<aod::Dielectrons, aod::DielectronsExtra>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMapWithMults = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionMult | VarManager::ObjTypes::CollisionMultExtra;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkTrackFillMapWithCovNoTOF = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackTPCPID | VarManager::ObjTypes::TrackTOFService;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups); // defines histograms for all tasks

// Enum containing the ordering of statistics histograms to be written in the QA file
enum ZorroStatHist {
  kStatsZorroInfo = 0,
  kStatsZorroSel
};

struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TList> fStatsList{"Statistics"};
  Configurable<std::string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<std::string> fConfigMixingVariablesJson{"cfgMixingVarsJSON", "", "Mixing configs in JSON format"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigEventCutsJSON{"cfgEventCutsJSON", "", "Additional event cuts specified in JSON format"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Add event histograms defined via JSON formatting (see HistogramsLibrary)"};

  Configurable<float> fConfigSplitCollisionsDeltaZ{"cfgSplitCollisionsDeltaZ", 1.0, "maximum delta-z (cm) between two collisions to consider them as split candidates"};
  Configurable<unsigned int> fConfigSplitCollisionsDeltaBC{"cfgSplitCollisionsDeltaBC", 100, "maximum delta-BC between two collisions to consider them as split candidates; do not apply if value is negative"};
  Configurable<bool> fConfigCheckSplitCollisions{"cfgCheckSplitCollisions", false, "If true, run the split collision check and fill histograms"};

  // Zorro selection
  struct : ConfigurableGroup {
    Configurable<bool> fConfigRunZorro{"cfgRunZorro", false, "Enable event selection with zorro"};
    Configurable<std::string> fConfigZorroTrigMask{"cfgZorroTriggerMask", "fDiMuon", "DQ Trigger masks: fSingleE,fLMeeIMR,fLMeeHMR,fDiElectron,fSingleMuLow,fSingleMuHigh,fDiMuon"};
    Configurable<bool> fConfigRunZorroSel{"cfgRunZorroSel", false, "Select events with trigger mask"};
    Configurable<uint64_t> fBcTolerance{"cfgBcTolerance", 100, "Number of BCs of margin for software triggers"};
    Configurable<std::string> fConfigCcdbPathZorro{"ccdb-path-zorro", "/Users/m/mpuccio/EventFiltering/OTS/Chunked/", "base path to the ccdb object for zorro"};
  } fConfigZorro;

  // RCT selection
  struct : ConfigurableGroup {
    Configurable<bool> fConfigUseRCT{"cfgUseRCT", false, "Enable event selection with RCT flags"};
    Configurable<bool> fCheckZDC{"cfgCheckZDC", false, "Check ZDC quality in the RCT flag checker"};
    Configurable<std::string> fConfigRCTLabel{"cfgRCTLabel", "CBT", "RCT flag labels : CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo"};
  } fConfigRCT;

  // CCDB connection configurables
  struct : ConfigurableGroup {
    Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/z/zhxiong/TPCPID/PostCalib", "base path to the ccdb object"};
    Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
    Configurable<std::string> fConfigGeoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> fConfigGrpMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> fZShiftPath{"zShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
    Configurable<std::string> fConfigGrpMagPathRun2{"grpmagPathRun2", "GLO/GRP/GRP", "CCDB path of the GRPObject (Usage for Run 2)"};
  } fConfigCCDB;

  // TPC postcalibration related options
  struct : ConfigurableGroup {
    Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas(electrons, pions, protons)"};
    Configurable<int> fConfigTPCpostCalibType{"cfgTPCpostCalibType", 1, "1: (TPCncls,pIN,eta) calibration typically for pp, 2: (eta,nPV,nLong,tLong) calibration typically for PbPb"};
    Configurable<bool> fConfigTPCuseInterpolatedCalib{"cfgTPCpostCalibUseInterpolation", true, "If true, use interpolated calibration values (default: true)"};
    Configurable<bool> fConfigComputeTPCpostCalibKaon{"cfgTPCpostCalibKaon", false, "If true, compute TPC post-calibrated n-sigmas for kaons"};
    Configurable<bool> fConfigIsOnlyforMaps{"cfgIsforMaps", false, "If true, run for postcalibration maps only"};
    Configurable<bool> fConfigSaveElectronSample{"cfgSaveElectronSample", false, "If true, only save electron sample"};
  } fConfigPostCalibTPC;

  Configurable<bool> fIsRun2{"cfgIsRun2", false, "Whether we analyze Run-2 or Run-3 data"};

  // RCT flag checker
  o2::aod::rctsel::RCTFlagsChecker rctChecker{"CBT"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;

  AnalysisCompositeCut* fEventCut;

  o2::parameters::GRPObject* fGrpMagRun2 = nullptr; // for run 2, we access the GRPObject from GLO/GRP/GRP
  o2::parameters::GRPMagField* fGrpMag = nullptr;   // for run 3, we access GRPMagField from GLO/Config/GRPMagField

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  std::map<int64_t, bool> fSelMap;                     // key: reduced event global index, value: event selection decision
  std::map<uint64_t, std::vector<int64_t>> fBCCollMap; // key: global BC, value: vector of reduced event global indices
  int fCurrentRun;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    if (eventCutStr != "") {
      AnalysisCut* cut = dqcuts::GetAnalysisCut(eventCutStr.Data());
      if (cut != nullptr) {
        fEventCut->AddCut(cut);
      }
    }

    // Additional cuts via JSON
    TString eventCutJSONStr = fConfigEventCutsJSON.value;
    if (eventCutJSONStr != "") {
      std::vector<AnalysisCut*> jsonCuts = dqcuts::GetCutsFromJSON(eventCutJSONStr.Data());
      for (auto& cutIt : jsonCuts) {
        fEventCut->AddCut(cutIt);
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(true);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "TimeFrameStats;Event_BeforeCuts;Event_AfterCuts;", fConfigAddEventHistogram.value.data());
      if (fConfigCheckSplitCollisions) {
        DefineHistograms(fHistMan, "OutOfBunchCorrelations;SameBunchCorrelations;", "");
      }
      dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str());
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    // Zorro information: kStatsZorroInfo
    // Zorro trigger selection: kStatsZorroSel
    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);
    TH2D* histZorroInfo = new TH2D("ZorroInfo", "Zorro information", 1, -0.5, 0.5, 1, -0.5, 0.5);
    fStatsList->AddAt(histZorroInfo, kStatsZorroInfo);

    TH2D* histZorroSel = new TH2D("ZorroSel", "trigger of interested", 1, -0.5, 0.5, 1, -0.5, 0.5);
    fStatsList->AddAt(histZorroSel, kStatsZorroSel);

    TString mixVarsString = fConfigMixingVariables.value;
    TString mixVarsJsonString = fConfigMixingVariablesJson.value;
    std::unique_ptr<TObjArray> objArray(mixVarsString.Tokenize(","));
    if (objArray->GetEntries() > 0 || mixVarsJsonString != "") {
      fMixHandler = new MixingHandler("mixingHandler", "mixing handler");
      fMixHandler->Init();
      if (objArray->GetEntries() > 0) {
        for (int iVar = 0; iVar < objArray->GetEntries(); ++iVar) {
          dqmixing::SetUpMixing(fMixHandler, objArray->At(iVar)->GetName());
        }
      }
      if (mixVarsJsonString != "") {
        dqmixing::SetUpMixingFromJSON(fMixHandler, mixVarsJsonString.Data());
      }
    }

    fCurrentRun = -1;
    fCCDB->setURL(fConfigCCDB.fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigCCDB.fConfigNoLaterThan.value);
    fCCDBApi.init(fConfigCCDB.fConfigCcdbUrl.value);

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      fCCDB->get<TGeoManager>(fConfigCCDB.fConfigGeoPath.value);
    }

    if (fConfigRCT.fConfigUseRCT.value) {
      rctChecker.init(fConfigRCT.fConfigRCTLabel, fConfigRCT.fCheckZDC.value);
    }

  }

  template <uint32_t TEventFillMap, typename TEvents>
  void runEventSelection(TEvents const& events, BCsWithTimestamps const& bcs)
  {

    if (bcs.size() > 0 && fCurrentRun != bcs.begin().runNumber()) {
      if (fConfigPostCalibTPC.fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCCDB.fConfigCcdbPathTPC.value, bcs.begin().timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
        if (fConfigPostCalibTPC.fConfigComputeTPCpostCalibKaon) {
          VarManager::SetCalibrationObject(VarManager::kTPCKaonMean, calibList->FindObject("mean_map_kaon"));
          VarManager::SetCalibrationObject(VarManager::kTPCKaonSigma, calibList->FindObject("sigma_map_kaon"));
        }
        if (fConfigPostCalibTPC.fConfigTPCpostCalibType == 2) {
          VarManager::SetCalibrationObject(VarManager::kTPCElectronStatus, calibList->FindObject("status_map_electron"));
          VarManager::SetCalibrationObject(VarManager::kTPCPionStatus, calibList->FindObject("status_map_pion"));
          VarManager::SetCalibrationObject(VarManager::kTPCProtonStatus, calibList->FindObject("status_map_proton"));
          if (fConfigPostCalibTPC.fConfigComputeTPCpostCalibKaon) {
            VarManager::SetCalibrationObject(VarManager::kTPCKaonStatus, calibList->FindObject("status_map_kaon"));
          }
        }
        VarManager::SetCalibrationType(fConfigPostCalibTPC.fConfigTPCpostCalibType, fConfigPostCalibTPC.fConfigTPCuseInterpolatedCalib);
      }
      if (fIsRun2 == true) {
        fGrpMagRun2 = fCCDB->getForTimeStamp<o2::parameters::GRPObject>(fConfigCCDB.fConfigGrpMagPathRun2, bcs.begin().timestamp());
        if (fGrpMagRun2 != nullptr) {
          o2::base::Propagator::initFieldFromGRP(fGrpMagRun2);
        }
      } else {
        fGrpMag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigCCDB.fConfigGrpMagPath, bcs.begin().timestamp());
        auto* fZShift = fCCDB->getForTimeStamp<std::vector<float>>(fConfigCCDB.fZShiftPath, bcs.begin().timestamp());
        if (fGrpMag != nullptr) {
          o2::base::Propagator::initFieldFromGRP(fGrpMag);
          VarManager::SetMagneticField(fGrpMag->getNominalL3Field());
        }
        if (fZShift != nullptr && !fZShift->empty()) {
          VarManager::SetZShift((*fZShift)[0]);
        }
        /*if (fConfigVariousOptions.fPropMuon) {
          VarManager::SetupMuonMagField();
        }*/
      }
      std::map<std::string, std::string> metadataRCT, header;
      header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", bcs.begin().runNumber()), metadataRCT, -1);
      uint64_t sor = std::atol(header["SOR"].c_str());
      uint64_t eor = std::atol(header["EOR"].c_str());
      VarManager::SetSORandEOR(sor, eor);

      fCurrentRun = bcs.begin().runNumber();
    } // end updating the CCDB quantities at change of run

    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillTimeFrame(bcs);
    VarManager::FillTimeFrame(events);
    if (fConfigQA) {
      fHistMan->FillHistClass("TimeFrameStats", VarManager::fgValues);
    }

    fSelMap.clear();
    fBCCollMap.clear();

    for (auto& event : events) {
      auto bc = event.template bc_as<BCsWithTimestamps>();

      VarManager::ResetValues(VarManager::kNTFWiseVariables, VarManager::kNEventWiseVariables);
      VarManager::FillBC(bc);
      VarManager::FillEvent<TEventFillMap>(event);

      bool decision = false;
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }

      if (fConfigZorro.fConfigRunZorro) {
        zorro.setBaseCCDBPath(fConfigZorro.fConfigCcdbPathZorro.value);
        zorro.setBCtolerance(fConfigZorro.fBcTolerance);
        zorro.initCCDB(fCCDB.service, fCurrentRun, bc.timestamp(), fConfigZorro.fConfigZorroTrigMask.value);
        zorro.populateExternalHists(fCurrentRun, reinterpret_cast<TH2D*>(fStatsList->At(kStatsZorroInfo)), reinterpret_cast<TH2D*>(fStatsList->At(kStatsZorroSel)));

        if (!fEventCut->IsSelected(VarManager::fgValues) || (fConfigRCT.fConfigUseRCT.value && !rctChecker(event))) {
          continue;
        }

        bool zorroSel = zorro.isSelected(bc.globalBC(), fConfigZorro.fBcTolerance, reinterpret_cast<TH2D*>(fStatsList->At(kStatsZorroSel)));
        if (fConfigZorro.fConfigRunZorroSel && (!zorroSel)) {
          continue;
        }
      } else {

        if (!fEventCut->IsSelected(VarManager::fgValues) || (fConfigRCT.fConfigUseRCT.value && !rctChecker(event))) {
          continue;
        }
      }

      decision = true;
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }

      fSelMap[event.globalIndex()] = decision;
      if (fBCCollMap.find(bc.globalBC()) == fBCCollMap.end()) {
        std::vector<int64_t> evIndices = {event.globalIndex()};
        fBCCollMap[bc.globalBC()] = evIndices;
      } else {
        auto& evIndices = fBCCollMap[bc.globalBC()];
        evIndices.push_back(event.globalIndex());
      }
      if (fMixHandler != nullptr) {
        int hh = fMixHandler->FindEventCategory(VarManager::fgValues);
        hash(hh);
      }
    }

  }

  // Variant of runEventSelection that first checks the DqFilters EMu prefilter bit.
  // Events not passing the EMu filter bit are skipped entirely, reducing track/muon
  // propagation and PID computation for the majority of collisions.
  template <uint32_t TEventFillMap, typename TEvents>
  void runEventSelectionWithFilter(TEvents const& events, BCsWithTimestamps const& bcs)
  {

    if (bcs.size() > 0 && fCurrentRun != bcs.begin().runNumber()) {
      if (fConfigPostCalibTPC.fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCCDB.fConfigCcdbPathTPC.value, bcs.begin().timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
        if (fConfigPostCalibTPC.fConfigComputeTPCpostCalibKaon) {
          VarManager::SetCalibrationObject(VarManager::kTPCKaonMean, calibList->FindObject("mean_map_kaon"));
          VarManager::SetCalibrationObject(VarManager::kTPCKaonSigma, calibList->FindObject("sigma_map_kaon"));
        }
        if (fConfigPostCalibTPC.fConfigTPCpostCalibType == 2) {
          VarManager::SetCalibrationObject(VarManager::kTPCElectronStatus, calibList->FindObject("status_map_electron"));
          VarManager::SetCalibrationObject(VarManager::kTPCPionStatus, calibList->FindObject("status_map_pion"));
          VarManager::SetCalibrationObject(VarManager::kTPCProtonStatus, calibList->FindObject("status_map_proton"));
          if (fConfigPostCalibTPC.fConfigComputeTPCpostCalibKaon) {
            VarManager::SetCalibrationObject(VarManager::kTPCKaonStatus, calibList->FindObject("status_map_kaon"));
          }
        }
        VarManager::SetCalibrationType(fConfigPostCalibTPC.fConfigTPCpostCalibType, fConfigPostCalibTPC.fConfigTPCuseInterpolatedCalib);
      }
      if (fIsRun2 == true) {
        fGrpMagRun2 = fCCDB->getForTimeStamp<o2::parameters::GRPObject>(fConfigCCDB.fConfigGrpMagPathRun2, bcs.begin().timestamp());
        if (fGrpMagRun2 != nullptr) {
          o2::base::Propagator::initFieldFromGRP(fGrpMagRun2);
        }
      } else {
        fGrpMag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigCCDB.fConfigGrpMagPath, bcs.begin().timestamp());
        auto* fZShift = fCCDB->getForTimeStamp<std::vector<float>>(fConfigCCDB.fZShiftPath, bcs.begin().timestamp());
        if (fGrpMag != nullptr) {
          o2::base::Propagator::initFieldFromGRP(fGrpMag);
          VarManager::SetMagneticField(fGrpMag->getNominalL3Field());
        }
        if (fZShift != nullptr && !fZShift->empty()) {
          VarManager::SetZShift((*fZShift)[0]);
        }
      }
      std::map<std::string, std::string> metadataRCT, header;
      header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", bcs.begin().runNumber()), metadataRCT, -1);
      uint64_t sor = std::atol(header["SOR"].c_str());
      uint64_t eor = std::atol(header["EOR"].c_str());
      VarManager::SetSORandEOR(sor, eor);

      fCurrentRun = bcs.begin().runNumber();
    } // end updating the CCDB quantities at change of run

    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillTimeFrame(bcs);
    VarManager::FillTimeFrame(events);
    if (fConfigQA) {
      fHistMan->FillHistClass("TimeFrameStats", VarManager::fgValues);
    }

    fSelMap.clear();
    fBCCollMap.clear();

    for (auto& event : events) {
      // Skip events that did not pass any filterPP selection.
      // The bit position depends on filterPP config (fNBarrelCuts + fNMuonCuts + emu_index),
      // so check eventFilter != 0 rather than a hardcoded bit.
      if (event.eventFilter() == 0) {
        continue;
      }

      auto bc = event.template bc_as<BCsWithTimestamps>();

      VarManager::ResetValues(VarManager::kNTFWiseVariables, VarManager::kNEventWiseVariables);
      VarManager::FillBC(bc);
      VarManager::FillEvent<TEventFillMap>(event);

      bool decision = false;
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }

      if (fConfigZorro.fConfigRunZorro) {
        zorro.setBaseCCDBPath(fConfigZorro.fConfigCcdbPathZorro.value);
        zorro.setBCtolerance(fConfigZorro.fBcTolerance);
        zorro.initCCDB(fCCDB.service, fCurrentRun, bc.timestamp(), fConfigZorro.fConfigZorroTrigMask.value);
        zorro.populateExternalHists(fCurrentRun, reinterpret_cast<TH2D*>(fStatsList->At(kStatsZorroInfo)), reinterpret_cast<TH2D*>(fStatsList->At(kStatsZorroSel)));

        if (!fEventCut->IsSelected(VarManager::fgValues) || (fConfigRCT.fConfigUseRCT.value && !rctChecker(event))) {
          continue;
        }

        bool zorroSel = zorro.isSelected(bc.globalBC(), fConfigZorro.fBcTolerance, reinterpret_cast<TH2D*>(fStatsList->At(kStatsZorroSel)));
        if (fConfigZorro.fConfigRunZorroSel && (!zorroSel)) {
          continue;
        }
      } else {

        if (!fEventCut->IsSelected(VarManager::fgValues) || (fConfigRCT.fConfigUseRCT.value && !rctChecker(event))) {
          continue;
        }
      }

      decision = true;
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }

      fSelMap[event.globalIndex()] = decision;
      if (fBCCollMap.find(bc.globalBC()) == fBCCollMap.end()) {
        std::vector<int64_t> evIndices = {event.globalIndex()};
        fBCCollMap[bc.globalBC()] = evIndices;
      } else {
        auto& evIndices = fBCCollMap[bc.globalBC()];
        evIndices.push_back(event.globalIndex());
      }
      if (fMixHandler != nullptr) {
        int hh = fMixHandler->FindEventCategory(VarManager::fgValues);
        hash(hh);
      }
    }

  }

  template <uint32_t TEventFillMap, typename TEvents>
  void publishSelections(TEvents const& events)
  {
    std::map<int64_t, bool> collisionSplittingMap; // key: event global index, value: whether pileup event is a possible splitting

    // Reset the fValues array and fill event observables
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);

    // loop over the BC map, get the collision vectors and make in-bunch and out of bunch 2-event correlations
    for (auto bc1It = fBCCollMap.begin(); bc1It != fBCCollMap.end(); ++bc1It) {
      uint64_t bc1 = bc1It->first;
      auto const& bc1Events = bc1It->second;

      // same bunch event correlations, if more than 1 collisions in this bunch
      if (bc1Events.size() > 1) {
        for (auto ev1It = bc1Events.begin(); ev1It != bc1Events.end(); ++ev1It) {
          auto ev1 = events.rawIteratorAt(*ev1It);
          for (auto ev2It = std::next(ev1It); ev2It != bc1Events.end(); ++ev2It) {
            auto ev2 = events.rawIteratorAt(*ev2It);
            // compute 2-event quantities and mark the candidate split collisions
            VarManager::FillTwoEvents(ev1, ev2);
            if (TMath::Abs(VarManager::fgValues[VarManager::kTwoEvDeltaZ]) < fConfigSplitCollisionsDeltaZ) { // this is a possible collision split
              collisionSplittingMap[*ev1It] = true;
              collisionSplittingMap[*ev2It] = true;
            }
            if (fConfigQA) {
              fHistMan->FillHistClass("SameBunchCorrelations", VarManager::fgValues);
            }
          } // end second event loop
        } // end first event loop
      } // end if BC1 events > 1

      // loop over the following BCs in the TF
      for (auto bc2It = std::next(bc1It); bc2It != fBCCollMap.end(); ++bc2It) {
        uint64_t bc2 = bc2It->first;
        if ((bc2 > bc1 ? bc2 - bc1 : bc1 - bc2) > fConfigSplitCollisionsDeltaBC) {
          break;
        }
        auto const& bc2Events = bc2It->second;

        // loop over events in the first BC
        for (auto ev1It : bc1Events) {
          auto ev1 = events.rawIteratorAt(ev1It);
          // loop over events in the second BC
          for (auto ev2It : bc2Events) {
            auto ev2 = events.rawIteratorAt(ev2It);
            // compute 2-event quantities and mark the candidate split collisions
            VarManager::FillTwoEvents(ev1, ev2);
            if (TMath::Abs(VarManager::fgValues[VarManager::kTwoEvDeltaZ]) < fConfigSplitCollisionsDeltaZ) { // this is a possible collision split
              collisionSplittingMap[ev1It] = true;
              collisionSplittingMap[ev2It] = true;
            }
            if (fConfigQA) {
              fHistMan->FillHistClass("OutOfBunchCorrelations", VarManager::fgValues);
            }
          }
        }
      }
    }

    // publish the table
    uint32_t evSel = static_cast<uint32_t>(0);
    for (auto& event : events) {
      evSel = 0;
      if (fSelMap[event.globalIndex()]) { // event passed the user cuts
        evSel |= (static_cast<uint32_t>(1) << 0);
      }
      auto bc = event.template bc_as<BCsWithTimestamps>();
      std::vector<int64_t> sameBunchEvents = fBCCollMap[bc.globalBC()];
      if (sameBunchEvents.size() > 1) { // event with in-bunch pileup
        evSel |= (static_cast<uint32_t>(1) << 1);
      }
      if (collisionSplittingMap.find(event.globalIndex()) != collisionSplittingMap.end()) { // event with possible fake in-bunch pileup (collision splitting)
        evSel |= (static_cast<uint32_t>(1) << 2);
      }
      eventSel(evSel);
    }
  }

  void processDirect(MyEvents const& events, BCsWithTimestamps const& bcs)
  {
    runEventSelection<gkEventFillMapWithMults>(events, bcs);
    publishSelections<gkEventFillMapWithMults>(events);
  }

  void processDirectWithFilter(MyEventsWithDqFilter const& events, BCsWithTimestamps const& bcs)
  {
    runEventSelectionWithFilter<gkEventFillMapWithMults>(events, bcs);
    publishSelections<gkEventFillMapWithMults>(events);
  }

  void processDummy(aod::Collisions&) {}

  PROCESS_SWITCH(AnalysisEventSelection, processDirect, "Run event selection on framework AO2Ds", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDirectWithFilter, "Run event selection on framework AO2Ds with DqFilters EMu prefilter", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy function", true);
};

struct AnalysisTrackSelection {

  Produces<aod::BarrelTrackCuts> trackSel;
  Produces<aod::BarrelAmbiguities> trackAmbiguities;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<std::string> fConfigCuts{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigCutsJSON{"cfgBarrelTrackCutsJSON", "", "Additional list of barrel track cuts in JSON format"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
  Configurable<bool> fConfigPublishAmbiguity{"cfgPublishAmbiguity", true, "If true, publish ambiguity table and fill QA histograms"};
  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/z/zhxiong/TPCPID/PostCalib", "base path to the ccdb object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  Service<o2::pid::tof::TOFResponse> fTofResponse;

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut*> fTrackCuts;
  std::vector<TString> fHistNamesReco;

  int fCurrentRun; // current run (needed to detect run changes for loading CCDB parameters)

  std::map<int64_t, std::vector<int64_t>> fNAssocsInBunch;    // key: track global index, value: vector of global index for events associated in-bunch (events that have in-bunch pileup or splitting)
  std::map<int64_t, std::vector<int64_t>> fNAssocsOutOfBunch; // key: track global index, value: vector of global index for events associated out-of-bunch (events that have no in-bunch pileup)

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy"))
      return;

    VarManager::SetDefaultVarNames();
    fCurrentRun = 0;

    // Setup track cuts
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // Extra cuts from JSON
    TString addTrackCutsStr = fConfigCutsJSON.value;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (auto& t : addTrackCuts) {
        fTrackCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    // Setup histogram manager
    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // Configure histogram classes for each track cut;
      TString histClasses = "TimeFrameStats;AssocsBarrel_BeforeCuts;";
      for (auto& cut : fTrackCuts) {
        TString nameStr = Form("AssocsBarrel_%s", cut->GetName());
        fHistNamesReco.push_back(nameStr);
        histClasses += Form("%s;", nameStr.Data());
      }

      DefineHistograms(fHistMan, histClasses.Data(), fConfigAddTrackHistogram.value.data());
      if (fConfigPublishAmbiguity) {
        DefineHistograms(fHistMan, "TrackBarrel_AmbiguityInBunch;TrackBarrel_AmbiguityOutOfBunch;", "ambiguity");
      }

      dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                       // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    fTofResponse->initSetup(fCCDB, context);

  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTracks>
  void runTrackSelection(TrackAssoc const& assocs, BCsWithTimestamps const& bcs, TEvents const& events, TTracks const& tracks)
  {
    // determine if TEvents table contains aod::Collisions
    // bool hasCollisions = std::is_same<typename TEvents::BaseType, aod::Collisions>::value;

    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillTimeFrame(events);
    VarManager::FillTimeFrame(tracks);
    if (fConfigQA)
      fHistMan->FillHistClass("TimeFrameStats", VarManager::fgValues);

    if (bcs.size() > 0 && fCurrentRun != bcs.begin().runNumber()) {
      if (fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, bcs.begin().timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
      }

      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bcs.begin().timestamp());
      if (grpmag != nullptr) {
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", bcs.begin().timestamp());
      }
      fCurrentRun = bcs.begin().runNumber();
    }

    trackSel.reserve(assocs.size());
    trackAmbiguities.reserve(tracks.size());


    for (auto& assoc : assocs) {
      auto event = assoc.template collision_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        trackSel(0);
        continue;
      }

      ///
      auto track = tracks.rawIteratorAt(assoc.trackId());
      auto evFromTrack = events.rawIteratorAt(track.collisionId());
      if (!evFromTrack.isEventSelected_bit(0)) {
        trackSel(0);
        continue;
      }

      VarManager::ResetValues(VarManager::kNTFWiseVariables, VarManager::kNBarrelTrackVariables);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);

      VarManager::FillTrack<TTrackFillMap>(track);
      // compute quantities which depend on the associated collision, such as DCA
      if (track.collisionId() != event.globalIndex())
        VarManager::FillTrackCollision<TTrackFillMap>(track, event);

      if (fConfigQA)
        fHistMan->FillHistClass("AssocsBarrel_BeforeCuts", VarManager::fgValues);

      int iCut = 0;
      uint32_t filterMap = static_cast<uint32_t>(0);
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, iCut++) {
        if ((*cut)->IsSelected(VarManager::fgValues)) {
          filterMap |= (static_cast<uint32_t>(1) << iCut);
          if (fConfigQA) {
            fHistMan->FillHistClass(fHistNamesReco[iCut], VarManager::fgValues);
          }
        }
      } // end loop over cuts
      trackSel(filterMap);

      // count the number of associations per track
      if (fConfigPublishAmbiguity && filterMap > 0) {
        if (event.isEventSelected_bit(1)) {
          // for this track, count the number of associated collisions with in-bunch pileup and out of bunch associations
          if (fNAssocsInBunch.find(track.globalIndex()) == fNAssocsInBunch.end()) {
            std::vector<int64_t> evVector = {event.globalIndex()};
            fNAssocsInBunch[track.globalIndex()] = evVector;
          } else {
            auto& evVector = fNAssocsInBunch[track.globalIndex()];
            evVector.push_back(event.globalIndex());
          }
        } else {
          if (fNAssocsOutOfBunch.find(track.globalIndex()) == fNAssocsOutOfBunch.end()) {
            std::vector<int64_t> evVector = {event.globalIndex()};
            fNAssocsOutOfBunch[track.globalIndex()] = evVector;
          } else {
            auto& evVector = fNAssocsOutOfBunch[track.globalIndex()];
            evVector.push_back(event.globalIndex());
          }
        }
      }
    } // end loop over associations

    //  QA the collision-track associations
    //  TODO: some tracks can be associated to both collisions that have in bunch pileup and collisions from different bunches
    //        So one could QA these tracks separately
    if (fConfigPublishAmbiguity) {
      if (fConfigQA) {
        for (auto& [trackIdx, evIndices] : fNAssocsInBunch) {
          if (evIndices.size() == 1) {
            continue;
          }
          auto track = tracks.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
          VarManager::FillTrack<TTrackFillMap>(track);
          VarManager::fgValues[VarManager::kBarrelNAssocsInBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackBarrel_AmbiguityInBunch", VarManager::fgValues);
        } // end loop over in-bunch ambiguous tracks

        for (auto& [trackIdx, evIndices] : fNAssocsOutOfBunch) {
          if (evIndices.size() == 1) {
            continue;
          }
          auto track = tracks.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
          VarManager::FillTrack<TTrackFillMap>(track);
          VarManager::fgValues[VarManager::kBarrelNAssocsOutOfBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackBarrel_AmbiguityOutOfBunch", VarManager::fgValues);
        } // end loop over out-of-bunch ambiguous tracks
      }

      // publish the ambiguity table
      for (auto& track : tracks) {
        int8_t nInBunch = 0;
        if (fNAssocsInBunch.find(track.globalIndex()) != fNAssocsInBunch.end()) {
          nInBunch = fNAssocsInBunch[track.globalIndex()].size();
        }
        int8_t nOutOfBunch = 0;
        if (fNAssocsOutOfBunch.find(track.globalIndex()) != fNAssocsOutOfBunch.end()) {
          nOutOfBunch = fNAssocsOutOfBunch[track.globalIndex()].size();
        }
        trackAmbiguities(nInBunch, nOutOfBunch);
      }
    }

  }

  void processWithCov(TrackAssoc const& assocs, BCsWithTimestamps const& bcs, MyEventsSelected const& events, MyBarrelTracksWithCov const& tracks)
  {
    runTrackSelection<gkEventFillMapWithMults, gkTrackFillMapWithCov>(assocs, bcs, events, tracks);
  }

  void processWithCovTOFService(TrackAssoc const& assocs, BCsWithTimestamps const& bcs, MyEventsSelected const& events, MyBarrelTracksWithCovNoTOF const& tracks)
  {
    fTofResponse->processSetup(bcs.iteratorAt(0));
    auto tracksWithTOFservice = soa::Attach<MyBarrelTracksWithCovNoTOF, o2::aod::TOFNSigmaDynEl, o2::aod::TOFNSigmaDynPi,
                                            o2::aod::TOFNSigmaDynKa, o2::aod::TOFNSigmaDynPr>(tracks);
    runTrackSelection<gkEventFillMapWithMults, gkTrackFillMapWithCovNoTOF>(assocs, bcs, events, tracksWithTOFservice);
  }

  void processDummy(MyEvents&) {}

  PROCESS_SWITCH(AnalysisTrackSelection, processWithCov, "Run barrel track selection on DQ skimmed tracks w/ cov matrix associations", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processWithCovTOFService, "Run barrel track selection on DQ skimmed tracks w/ cov matrix associations, with TOF service", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy function", true);
};

struct AnalysisPrefilterSelection {
  Produces<aod::Prefilter> prefilter; // joinable with TracksAssoc

  // Configurables
  Configurable<std::string> fConfigPrefilterTrackCut{"cfgPrefilterTrackCut", "", "Prefilter track cut"};
  Configurable<std::string> fConfigPrefilterPairCut{"cfgPrefilterPairCut", "", "Prefilter pair cut"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Track cuts for which to run the prefilter"};
  // Track related options
  Configurable<bool> fPropTrack{"cfgPropTrack", false, "Propagate tracks to associated collision to recalculate DCA and momentum vector"};

  std::map<uint32_t, uint32_t> fPrefilterMap;
  AnalysisCompositeCut* fPairCut;
  uint32_t fPrefilterMask;
  int fPrefilterCutBit;

  Preslice<aod::TrackAssoc> trackAssocsPerCollision = aod::track_association::collisionId;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    bool runPrefilter = true;
    // get the list of track cuts to be prefiltered
    TString trackCutsStr = fConfigTrackCuts.value;
    TObjArray* objArrayTrackCuts = nullptr;
    if (!trackCutsStr.IsNull()) {
      objArrayTrackCuts = trackCutsStr.Tokenize(",");
      if (objArrayTrackCuts == nullptr) {
        runPrefilter = false;
      }
    } else {
      LOG(warn) << " No track cuts to prefilter! Prefilter will not be run";
      runPrefilter = false;
    }
    // get the cut to be used as loose selection
    TString prefilterTrackCutStr = fConfigPrefilterTrackCut.value;
    if (prefilterTrackCutStr.IsNull()) {
      LOG(warn) << " No prefilter loose selection specified! Prefilter will not be run";
      runPrefilter = false;
    }

    fPrefilterMask = 0;
    fPrefilterCutBit = -1;
    if (runPrefilter) {
      // get the list of cuts that were computed in the barrel track-selection task and create a bit mask
      //  to mark just the ones we want to apply a prefilter on
      string trackCuts;
      getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", trackCuts, false);
      TString allTrackCutsStr = trackCuts;
      // check also the cuts added via JSON and add them to the string of cuts
      getTaskOptionValue<string>(context, "analysis-track-selection", "cfgBarrelTrackCutsJSON", trackCuts, false);
      TString addTrackCutsStr = trackCuts;
      if (addTrackCutsStr != "") {
        std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
        for (auto& t : addTrackCuts) {
          allTrackCutsStr += Form(",%s", t->GetName());
        }
      }

      std::unique_ptr<TObjArray> objArray(allTrackCutsStr.Tokenize(","));
      if (objArray == nullptr) {
        LOG(fatal) << " Not getting any track cuts from the barrel-track-selection ";
      }
      if (objArray->FindObject(prefilterTrackCutStr.Data()) == nullptr) {
        LOG(fatal) << " Prefilter track cut not among the cuts calculated by the track-selection task! ";
      }
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        if (objArrayTrackCuts->FindObject(tempStr.Data()) != nullptr) {
          fPrefilterMask |= (static_cast<uint32_t>(1) << icut);
        }
        if (tempStr.CompareTo(fConfigPrefilterTrackCut.value) == 0) {
          fPrefilterCutBit = icut;
        }
      }
      // setup the prefilter pair cut
      fPairCut = new AnalysisCompositeCut(true);
      TString pairCutStr = fConfigPrefilterPairCut.value;
      if (!pairCutStr.IsNull()) {
        fPairCut = dqcuts::GetCompositeCut(pairCutStr.Data());
      }
    }
    if (fPrefilterMask == static_cast<uint32_t>(0) || fPrefilterCutBit < 0) {
      LOG(warn) << "No specified loose cut or track cuts for prefiltering. This task will do nothing.";
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();

    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  template <typename T>
  void runPrefilter(MyEvents::iterator const& event, soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts> const& assocs, T const& /*tracks*/)
  {
    // cout << "AnalysisPrefilterSelection::runPrefilter() called for event " << event.globalIndex() << " with " << assocs.size() << " track associations" << endl;
    if (fPrefilterCutBit < 0 || fPrefilterMask == 0) {
      return;
    }

    for (auto& [assoc1, assoc2] : o2::soa::combinations(assocs, assocs)) {
      auto track1 = assoc1.template track_as<T>();
      auto track2 = assoc2.template track_as<T>();

      // NOTE: here we restrict to just pairs of opposite sign (conversions), but in principle this can be made
      // a configurable and check also same-sign pairs (track splitting)
      if (track1.sign() * track2.sign() > 0) {
        continue;
      }

      // here we check the cuts fulfilled by both tracks, for both the tight and loose selections
      uint32_t track1Candidate = (assoc1.isBarrelSelected_raw() & fPrefilterMask);
      uint32_t track2Candidate = (assoc2.isBarrelSelected_raw() & fPrefilterMask);
      bool track1Loose = assoc1.isBarrelSelected_bit(fPrefilterCutBit);
      bool track2Loose = assoc2.isBarrelSelected_bit(fPrefilterCutBit);

      if (!((track1Candidate > 0 && track2Loose) || (track2Candidate > 0 && track1Loose))) {
        continue;
      }

      // compute pair quantities
      VarManager::FillPair<VarManager::kDecayToEE, gkTrackFillMapWithCov>(track1, track2);
      if (fPropTrack) {
        VarManager::FillPairCollision<VarManager::kDecayToEE, gkTrackFillMapWithCov>(event, track1, track2);
      }
      // if the pair fullfils the criteria, add an entry into the prefilter map for the two tracks
      if (fPairCut->IsSelected(VarManager::fgValues)) {
        if (fPrefilterMap.find(track1.globalIndex()) == fPrefilterMap.end() && track1Candidate > 0) {
          fPrefilterMap[track1.globalIndex()] = track1Candidate;
        }
        if (fPrefilterMap.find(track2.globalIndex()) == fPrefilterMap.end() && track2Candidate > 0) {
          fPrefilterMap[track2.globalIndex()] = track2Candidate;
        }
      }
    } // end loop over combinations
    // cout << "AnalysisPrefilterSelection::runPrefilter() completed for event " << event.globalIndex() << endl;
  }

  void processBarrel(MyEvents const& events, soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts> const& assocs, MyBarrelTracksWithCov const& tracks)
  {
    fPrefilterMap.clear();

    for (auto& event : events) {
      auto groupedAssocs = assocs.sliceBy(trackAssocsPerCollision, event.globalIndex());
      groupedAssocs.bindInternalIndicesTo(&assocs);

      if (groupedAssocs.size() > 1) {
        runPrefilter(event, groupedAssocs, tracks);
      }
    }

    uint32_t mymap = -1;
    // If cuts were not configured, then produce a map with all 1's and publish it for all associations
    if (fPrefilterCutBit < 0 || fPrefilterMask == 0) {
      for (int i = 0; i < assocs.size(); ++i) {
        prefilter(mymap);
      }
    } else {
      for (auto& assoc : assocs) {
        // TODO: just use the index from the assoc (no need to cast the whole track)
        // auto track = assoc.template track_as<MyBarrelTracksWithCov>();
        mymap = -1;
        // if (fPrefilterMap.find(track.globalIndex()) != fPrefilterMap.end()) {
        if (fPrefilterMap.find(assoc.trackId()) != fPrefilterMap.end()) {
          // NOTE: publish the bitwise negated bits (~), so there will be zeroes for cuts that failed the prefiltering and 1 everywhere else
          // mymap = ~fPrefilterMap[track.globalIndex()];
          mymap = ~fPrefilterMap[assoc.trackId()];
          prefilter(mymap);
        } else {
          prefilter(mymap); // track did not pass the prefilter selections, so publish just 1's
        }
      }
    }
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisPrefilterSelection, processBarrel, "Run Prefilter selection on barrel tracks", false);
  PROCESS_SWITCH(AnalysisPrefilterSelection, processDummy, "Do nothing", true);
};

// Produces a table with muon decisions (joinable to FwdTrackAssoc)
struct AnalysisMuonSelection {
  Produces<aod::MuonTrackCuts> muonSel;
  Produces<aod::MuonAmbiguities> muonAmbiguities;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<std::string> fConfigCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<std::string> fConfigCutsJSON{"cfgMuonCutsJSON", "", "Additional list of muon cuts in JSON format"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
  Configurable<bool> fConfigPublishAmbiguity{"cfgPublishAmbiguity", true, "If true, publish ambiguity table and fill QA histograms"};

  struct : ConfigurableGroup {
    Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<int64_t> noLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } fConfigCCDB;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan = nullptr;
  std::vector<AnalysisCompositeCut*> fMuonCuts;

  int fCurrentRun = 0;

  // key: FwdTrack global index, value: vector of collision global indices
  std::map<int64_t, std::vector<int64_t>> fNAssocsInBunch;
  std::map<int64_t, std::vector<int64_t>> fNAssocsOutOfBunch;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fMuonCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // extra cuts from JSON
    TString addCutsStr = fConfigCutsJSON.value;
    if (addCutsStr != "") {
      std::vector<AnalysisCut*> addCuts = dqcuts::GetCutsFromJSON(addCutsStr.Data());
      for (auto& t : addCuts) {
        fMuonCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars);

    if (fConfigQA) {
      if (fHistMan == nullptr) {
        fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
        fHistMan->SetUseDefaultVariableNames(kTRUE);
        fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

        TString histDirNames = "TrackMuon_BeforeCuts;";
        for (auto& cut : fMuonCuts) {
          histDirNames += Form("TrackMuon_%s;", cut->GetName());
        }
        if (fConfigPublishAmbiguity) {
          histDirNames += "TrackMuon_AmbiguityInBunch;TrackMuon_AmbiguityOutOfBunch;";
        }
        DefineHistograms(fHistMan, histDirNames.Data(), fConfigAddMuonHistogram.value.data());
        dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str());
        VarManager::SetUseVars(fHistMan->GetUsedVars());
        fOutputList.setObject(fHistMan->GetMainHistogramList());
      }
    }

    fCCDB->setURL(fConfigCCDB.url.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigCCDB.noLaterThan.value);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      fCCDB->get<TGeoManager>(fConfigCCDB.geoPath);
    }
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvents, typename TMuons>
  void runMuonSelection(BCsWithTimestamps const& bcs,
                        aod::FwdTrackAssoc const& assocs,
                        TEvents const& /*events*/, TMuons const& muons)
  {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

    if (bcs.size() > 0 && fCurrentRun != bcs.begin().runNumber()) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigCCDB.grpMagPath, bcs.begin().timestamp());
      if (grpmag != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", bcs.begin().timestamp());
      }
      fCurrentRun = bcs.begin().runNumber();
    }

    muonSel.reserve(assocs.size());
    if (fConfigPublishAmbiguity) {
      muonAmbiguities.reserve(muons.size());
    }
    uint32_t filterMap = static_cast<uint32_t>(0);
    int iCut = 0;

    for (auto& assoc : assocs) {
      auto event = assoc.template collision_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        muonSel(0);
        continue;
      }
      VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
      VarManager::FillEvent<TEventFillMap>(event);

      auto track = assoc.template fwdtrack_as<TMuons>();
      filterMap = static_cast<uint32_t>(0);
      VarManager::FillTrack<TMuonFillMap>(track);
      if (fConfigQA) {
        fHistMan->FillHistClass("TrackMuon_BeforeCuts", VarManager::fgValues);
      }
      iCut = 0;
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, iCut++) {
        if ((*cut)->IsSelected(VarManager::fgValues)) {
          filterMap |= (static_cast<uint32_t>(1) << iCut);
          if (fConfigQA) {
            fHistMan->FillHistClass(Form("TrackMuon_%s", (*cut)->GetName()), VarManager::fgValues);
          }
        }
      }
      muonSel(filterMap);

      if (fConfigPublishAmbiguity && filterMap > 0) {
        if (event.isEventSelected_bit(1)) { // in-bunch pileup flag
          if (fNAssocsInBunch.find(track.globalIndex()) == fNAssocsInBunch.end()) {
            fNAssocsInBunch[track.globalIndex()] = {event.globalIndex()};
          } else {
            fNAssocsInBunch[track.globalIndex()].push_back(event.globalIndex());
          }
        } else {
          if (fNAssocsOutOfBunch.find(track.globalIndex()) == fNAssocsOutOfBunch.end()) {
            fNAssocsOutOfBunch[track.globalIndex()] = {event.globalIndex()};
          } else {
            fNAssocsOutOfBunch[track.globalIndex()].push_back(event.globalIndex());
          }
        }
      }
    } // end loop over assocs

    if (fConfigPublishAmbiguity) {
      if (fConfigQA) {
        for (auto& [trackIdx, evIndices] : fNAssocsInBunch) {
          if (evIndices.size() <= 1)
            continue;
          auto track = muons.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
          VarManager::FillTrack<TMuonFillMap>(track);
          VarManager::fgValues[VarManager::kMuonNAssocsInBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackMuon_AmbiguityInBunch", VarManager::fgValues);
        }
        for (auto& [trackIdx, evIndices] : fNAssocsOutOfBunch) {
          if (evIndices.size() <= 1)
            continue;
          auto track = muons.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
          VarManager::FillTrack<TMuonFillMap>(track);
          VarManager::fgValues[VarManager::kMuonNAssocsOutOfBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackMuon_AmbiguityOutOfBunch", VarManager::fgValues);
        }
      }
      // publish ambiguity table (one row per FwdTrack)
      for (auto& track : muons) {
        int8_t nInBunch = 0;
        if (fNAssocsInBunch.find(track.globalIndex()) != fNAssocsInBunch.end()) {
          nInBunch = static_cast<int8_t>(fNAssocsInBunch[track.globalIndex()].size());
        }
        int8_t nOutOfBunch = 0;
        if (fNAssocsOutOfBunch.find(track.globalIndex()) != fNAssocsOutOfBunch.end()) {
          nOutOfBunch = static_cast<int8_t>(fNAssocsOutOfBunch[track.globalIndex()].size());
        }
        muonAmbiguities(nInBunch, nOutOfBunch);
      }
    }
  }

  void processDirect(MyEventsSelected const& events, BCsWithTimestamps const& bcs,
                     aod::FwdTrackAssoc const& assocs,
                     MyMuonTracksWithCov const& muons)
  {
    runMuonSelection<gkEventFillMapWithMults, gkMuonFillMapWithCov>(bcs, assocs, events, muons);
  }

  void processDummy(MyEvents&) { /* do nothing */ }

  PROCESS_SWITCH(AnalysisMuonSelection, processDirect, "Run muon selection on AO2D FwdTracks", false);
  PROCESS_SWITCH(AnalysisMuonSelection, processDummy, "Dummy function", true);
};

struct AnalysisSameEventPairing {
  Produces<aod::Dielectrons> dielectronList;
  Produces<aod::DielectronsExtra> dielectronsExtraList;
  Produces<aod::DielectronsInfo> dielectronInfoList;
  Produces<aod::DielectronsAll> dielectronAllList;
  Produces<aod::DileptonInfo> dileptonInfoList;
  Produces<aod::JPsieeCandidates> PromptNonPromptSepTable;
  Produces<aod::ElectronMuons> electronmuonList;

  o2::base::MatLayerCylSet* fLUT = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  OutputObj<THashList> fOutputList{"output"};

  // Histogram manager
  HistogramManager* fHistMan = nullptr;

  // Config options
  // ConfigOptions fConfigOptions;

  struct : ConfigurableGroup {
    Configurable<std::string> track{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
    Configurable<std::string> muon{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
    Configurable<std::string> pair{"cfgPairCuts", "", "Comma separated list of pair cuts, !!! Use only if you know what you are doing, otherwise leave empty"};
    Configurable<bool> collSplitting{"cfgRemoveCollSplittingCandidates", false, "If true, remove collision splitting candidates as determined by the event selection task upstream"};
    Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
    Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
    Configurable<bool> useRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
    Configurable<float> magField{"cfgMagField", 5.0f, "Manually set magnetic field"};
    Configurable<bool> flatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
    Configurable<bool> useKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
    Configurable<bool> recomputePV{"cfgRecomputePV", false, "Recompute primary vertex using PVertexer to calculate unbiased pseudoproperDL"};
    Configurable<bool> removeDiamondConstrPV{"cfgRemoveDiamondPV", false, "remove diamond constrain for PV recomputation"};
    Configurable<bool> useAbsDCA{"cfgUseAbsDCA", false, "Use absolute DCA minimization instead of chi^2 minimization in secondary vertexing"};
    Configurable<bool> propToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
    Configurable<bool> corrFullGeo{"cfgCorrFullGeo", false, "Use full geometry to correct for MCS effects in track propagation"};
    Configurable<bool> noCorr{"cfgNoCorrFwdProp", false, "Do not correct for MCS effects in track propagation"};
    Configurable<std::string> collisionSystem{"syst", "pp", "Collision system, pp or PbPb"};
    Configurable<float> centerMassEnergy{"energy", 13600, "Center of mass energy in GeV"};
    Configurable<bool> fPropTrack{"cfgPropTrack", true, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};
    Configurable<bool> fConfigMiniTree{"cfgMiniTree", false, "Produce a single flat table with minimal information for analysis"};
    Configurable<float> fConfigMiniTreeMinMass{"cfgMiniTreeMinMass", 2, "Min. mass cut for minitree"};
    Configurable<float> fConfigMiniTreeMaxMass{"cfgMiniTreeMaxMass", 5, "Max. mass cut for minitree"};
  } fConfigOptions;

  struct : ConfigurableGroup {
    Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } fConfigCCDB;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  // vectors needed for PV recomputation
  std::vector<int64_t> pvContribGlobIDs;
  std::vector<o2::track::TrackParCov> pvContribTrackPars;
  std::vector<bool> vec_useTrk_PVrefit;

  // keep histogram class names in maps, so we don't have to buld their names in the pair loops
  std::map<int, std::vector<TString>> fTrackHistNames;
  std::map<int, std::vector<TString>> fMuonHistNames;
  std::map<int, std::vector<TString>> fTrackMuonHistNames; // for electron-muon pairs: key = iTrack * fNCutsMuon + iMuon

  std::vector<AnalysisCompositeCut> fPairCuts;
  AnalysisCompositeCut fMCGenAccCut;
  // bool fUseMCGenAccCut = false;

  uint32_t fTrackFilterMask; // mask for the track cuts required in this task to be applied on the barrel cuts produced upstream
  uint32_t fMuonFilterMask;  // mask for the muon cuts required in this task to be applied on the muon cuts produced upstream
  int fNCutsBarrel;
  int fNCutsMuon;
  int fNPairCuts;
  bool fHasTwoProngGenMCsignals = false;

  bool fEnableBarrelHistos;
  bool fEnableBarrelMuonHistos;

  std::vector<TString> fTrackCuts; // barrel cut names, used in EMu histogram filling
  std::vector<TString> fMuonCuts;  // muon cut names, used in EMu histogram filling

  Preslice<soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts, aod::Prefilter>> trackAssocsPerCollision = aod::track_association::collisionId;
  Preslice<soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts>> trackEmuAssocsPerCollision = aod::track_association::collisionId;
  Preslice<soa::Join<aod::FwdTrackAssoc, aod::MuonTrackCuts>> muonAssocsPerCollision = aod::track_association::collisionId;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    fEnableBarrelHistos = context.mOptions.get<bool>("processBarrelOnly");
    fEnableBarrelMuonHistos = context.mOptions.get<bool>("processElectronMuonDirect");

    // Keep track of all the histogram class names to avoid composing strings in the pairing loop
    TString histNames = "";
    TString cutNamesStr = fConfigOptions.pair.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // get the list of cuts for tracks/muons, check that they were played by the barrel/muon selection tasks
    //   and make a mask for active cuts (barrel and muon selection tasks may run more cuts, needed for other analyses)
    TString trackCutsStr = fConfigOptions.track.value;
    TObjArray* objArrayTrackCuts = nullptr;
    if (!trackCutsStr.IsNull()) {
      objArrayTrackCuts = trackCutsStr.Tokenize(",");
    }
    TString muonCutsStr = fConfigOptions.muon.value;
    TObjArray* objArrayMuonCuts = nullptr;
    if (!muonCutsStr.IsNull()) {
      objArrayMuonCuts = muonCutsStr.Tokenize(",");
    }

    // get the barrel track selection cuts
    string tempCuts;
    getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", tempCuts, false);
    TString tempCutsStr = tempCuts;
    // check also the cuts added via JSON and add them to the string of cuts
    getTaskOptionValue<string>(context, "analysis-track-selection", "cfgBarrelTrackCutsJSON", tempCuts, false);
    TString addTrackCutsStr = tempCuts;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (auto& t : addTrackCuts) {
        tempCutsStr += Form(",%s", t->GetName());
      }
    }

    // check that the barrel track cuts array required in this task is not empty
    if (!trackCutsStr.IsNull()) {
      // tokenize and loop over the barrel cuts produced by the barrel track selection task
      std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
      fNCutsBarrel = objArray->GetEntries();
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        // if the current barrel selection cut is required in this task, then switch on the corresponding bit in the mask
        // and assign histogram directories
        if (objArrayTrackCuts->FindObject(tempStr.Data()) != nullptr) {
          fTrackFilterMask |= (static_cast<uint32_t>(1) << icut);
          fTrackCuts.push_back(tempStr);

          if (fEnableBarrelHistos) {
            // assign the pair hist directories for the current cut
            std::vector<TString> names = {
              Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
            if (fConfigOptions.fConfigQA) {
              // assign separate hist directories for ambiguous tracks
              names.push_back(Form("PairsBarrelSEPM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEPP_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEMM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEPM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEPP_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelSEMM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
            }
            for (auto& n : names) {
              histNames += Form("%s;", n.Data());
            }
            fTrackHistNames[icut] = names;

            // if there are pair cuts specified, assign hist directories for each barrel cut - pair cut combination
            // NOTE: This could possibly lead to large histogram outputs. It is strongly advised to use pair cuts only
            //   if you know what you are doing.
            TString cutNamesStr = fConfigOptions.pair.value;
            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              fNPairCuts = objArrayPair->GetEntries();
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                names = {
                  Form("PairsBarrelSEPM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsBarrelSEPP_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsBarrelSEMM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
                histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
                // NOTE: In the numbering scheme for the map key, we use the number of barrel cuts in the barrel-track selection task
                fTrackHistNames[fNCutsBarrel + icut * fNPairCuts + iPairCut] = names;
              } // end loop (pair cuts)
            } // end if (pair cuts)
          } // end if enableBarrelHistos
        }
      }
    }

    // get the muon track selection cuts (from analysis-muon-selection task)
    getTaskOptionValue<string>(context, "analysis-muon-selection", "cfgMuonCuts", tempCuts, false);
    tempCutsStr = tempCuts;
    // check also the cuts added via JSON
    getTaskOptionValue<string>(context, "analysis-muon-selection", "cfgMuonCutsJSON", tempCuts, false);
    TString addMuonCutsStr = tempCuts;
    if (addMuonCutsStr != "") {
      std::vector<AnalysisCut*> addMuonCuts = dqcuts::GetCutsFromJSON(addMuonCutsStr.Data());
      for (auto& t : addMuonCuts) {
        tempCutsStr += Form(",%s", t->GetName());
      }
    }

    // build fMuonFilterMask and, if needed, fTrackMuonHistNames for EMu pairing
    if (!muonCutsStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
      fNCutsMuon = objArray->GetEntries();
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        if (objArrayMuonCuts->FindObject(tempStr.Data()) != nullptr) {
          fMuonFilterMask |= (static_cast<uint32_t>(1) << icut);
          fMuonCuts.push_back(tempStr);

          if (fEnableBarrelMuonHistos) {
            // assign PairsEleMu histogram directories for each (barrel cut, muon cut) combination
            int seqTrackIdx = 0; // sequential index into fTrackCuts (which contains only required cuts)
            for (int iTrack = 0; iTrack < fNCutsBarrel; ++iTrack) {
              // skip barrel cuts not required in this task
              if (!(fTrackFilterMask & (static_cast<uint32_t>(1) << iTrack)))
                continue;
              TString trackCutName = fTrackCuts[seqTrackIdx];
              seqTrackIdx++;
              std::vector<TString> names = {
                Form("PairsEleMuSEPM_%s_%s", trackCutName.Data(), tempStr.Data()),
                Form("PairsEleMuSEPP_%s_%s", trackCutName.Data(), tempStr.Data()),
                Form("PairsEleMuSEMM_%s_%s", trackCutName.Data(), tempStr.Data())};
              histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());

              // pair-cut variants
              TString pairCutsStr = fConfigOptions.pair.value;
              if (!pairCutsStr.IsNull()) {
                std::unique_ptr<TObjArray> objArrayPair(pairCutsStr.Tokenize(","));
                int nPairCuts = objArrayPair->GetEntries();
                for (int iPairCut = 0; iPairCut < nPairCuts; ++iPairCut) {
                  names = {
                    Form("PairsEleMuSEPM_%s_%s_%s", trackCutName.Data(), tempStr.Data(), objArrayPair->At(iPairCut)->GetName()),
                    Form("PairsEleMuSEPP_%s_%s_%s", trackCutName.Data(), tempStr.Data(), objArrayPair->At(iPairCut)->GetName()),
                    Form("PairsEleMuSEMM_%s_%s_%s", trackCutName.Data(), tempStr.Data(), objArrayPair->At(iPairCut)->GetName())};
                  histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
                  int index = iTrack * (fNCutsMuon * nPairCuts) + icut * nPairCuts + iPairCut;
                  fTrackMuonHistNames[index] = names;
                }
              } else {
                int index = iTrack * fNCutsMuon + icut;
                fTrackMuonHistNames[index] = names;
              }
            } // end loop barrel cuts
          } // end if fEnableBarrelMuonHistos
        }
      } // end loop muon cuts
    } // end if (muonCutsStr)

    fCurrentRun = 0;

    fCCDB->setURL(fConfigCCDB.url.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();

    if (fConfigOptions.noCorr) {
      VarManager::SetupFwdDCAFitterNoCorr();
    } else if (fConfigOptions.corrFullGeo || (fConfigOptions.useKFVertexing && fConfigOptions.propToPCA)) {
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        fCCDB->get<TGeoManager>(fConfigCCDB.geoPath);
      }
    } else {
      fLUT = o2::base::MatLayerCylSet::rectifyPtrFromFile(fCCDB->get<o2::base::MatLayerCylSet>(fConfigCCDB.lutPath));
      VarManager::SetupMatLUTFwdDCAFitter(fLUT);
    }

    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    VarManager::SetCollisionSystem((TString)fConfigOptions.collisionSystem, fConfigOptions.centerMassEnergy); // set collision system and center of mass energy

    DefineHistograms(fHistMan, histNames.Data(), fConfigOptions.fConfigAddSEPHistogram.value.data());     // define all histograms
    dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigOptions.fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                                      // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

  }

  void initParamsFromCCDB(uint64_t timestamp, bool withTwoProngFitter = true)
  {
    if (fConfigOptions.useRemoteField.value) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigCCDB.grpMagPath, timestamp);
      o2::base::MatLayerCylSet* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(fCCDB->get<o2::base::MatLayerCylSet>(fConfigCCDB.lutPath));
      float magField = 0.0;
      if (grpmag != nullptr) {
        magField = grpmag->getNominalL3Field();
        o2::base::Propagator::initFieldFromGRP(grpmag);
        o2::base::Propagator::Instance()->setMatLUT(lut);
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", timestamp);
      }
      if (withTwoProngFitter) {
        if (fConfigOptions.useKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(magField);
        } else {
          VarManager::SetupTwoProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigOptions.useAbsDCA.value); // TODO: get these parameters from Configurables
          VarManager::SetupTwoProngFwdDCAFitter(magField, true, 200.0f, 1.0e-3f, 0.9f, fConfigOptions.useAbsDCA.value);
        }
      } else {
        VarManager::SetupTwoProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigOptions.useAbsDCA.value); // needed because take in varmanager Bz from fgFitterTwoProngBarrel for PhiV calculations
      }
    } else {
      if (withTwoProngFitter) {
        if (fConfigOptions.useKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(fConfigOptions.magField.value);
        } else {
          VarManager::SetupTwoProngDCAFitter(fConfigOptions.magField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigOptions.useAbsDCA.value); // TODO: get these parameters from Configurables
          VarManager::SetupTwoProngFwdDCAFitter(fConfigOptions.magField.value, true, 200.0f, 1.0e-3f, 0.9f, fConfigOptions.useAbsDCA.value);
        }
      } else {
        VarManager::SetupTwoProngDCAFitter(fConfigOptions.magField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigOptions.useAbsDCA.value); // needed because take in varmanager Bz from fgFitterTwoProngBarrel for PhiV calculations
      }
    }
  }

  template <typename Events, typename TTracks, typename Tracks>
  bool refitPVWithPVertexer(Events const& collision, TTracks const& tracks, Tracks const& t1, Tracks const& t2, o2::dataformats::VertexBase& pvRefitted)
  {
    // --- build PV contributor list ---
    pvContribGlobIDs.clear();
    pvContribTrackPars.clear();
    // int nMyPVContrib = 0; int nMyPVContribOrig = 0;
    for (auto const& trk : tracks) {
      // check if it is PV contributor
      if (!trk.isPVContributor())
        continue;
      // check if it contributes to the vtx of this collision
      if (trk.collisionId() != collision.globalIndex())
        continue;
      // nMyPVContribOrig++;
      // --- remove t1 and t2 if they are PV contributors ---
      if (trk.globalIndex() == t1.globalIndex() || trk.globalIndex() == t2.globalIndex())
        continue;
      // add tracks and parameters to the list
      pvContribGlobIDs.push_back(trk.globalIndex());
      pvContribTrackPars.push_back(getTrackParCov(trk));
      // nMyPVContrib++;
    }

    // cout << "contributors from collision: " << collision.numContrib() << " - from refitting: before -> " <<  nMyPVContribOrig << " after -> " << nMyPVContrib << endl;
    vec_useTrk_PVrefit.assign(pvContribGlobIDs.size(), true);
    // --- build VertexBase from event collision ---
    o2::dataformats::VertexBase Pvtx;
    Pvtx.setX(collision.posX());
    Pvtx.setY(collision.posY());
    Pvtx.setZ(collision.posZ());
    Pvtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(),
                collision.covXZ(), collision.covYZ(), collision.covZZ());

    // --- configure vertexer ---
    o2::vertexing::PVertexer vertexer;
    if (fConfigOptions.removeDiamondConstrPV) {
      o2::conf::ConfigurableParam::updateFromString("pvertexer.useMeanVertexConstraint=false");
    }
    vertexer.init();

    bool PVrefit_doable = vertexer.prepareVertexRefit(pvContribTrackPars, Pvtx);
    if (!PVrefit_doable)
      return false;

    // --- do the refit ---
    pvRefitted = vertexer.refitVertex(vec_useTrk_PVrefit, Pvtx);

    return true;
  }

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <bool TTwoProngFitter, int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTracks>
  void runSameEventPairing(TEvents const& events, BCsWithTimestamps const& bcs, Preslice<soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts, aod::Prefilter>>& preslice, soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& assocs, TTracks const& tracks)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    if (fCurrentRun != bcs.begin().runNumber()) {
      initParamsFromCCDB(bcs.begin().timestamp(), TTwoProngFitter);
      fCurrentRun = bcs.begin().runNumber();
    }

    TString cutNames = fConfigOptions.track.value;
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    int ncuts = fNCutsBarrel;

    uint32_t twoTrackFilter = 0;
    int sign1 = 0;
    int sign2 = 0;

    dielectronList.reserve(1);
    // dimuonList.reserve(1);
    dielectronsExtraList.reserve(1);
    // dimuonsExtraList.reserve(1);
    dielectronInfoList.reserve(1);
    dileptonInfoList.reserve(1);

    if (fConfigOptions.flatTables.value) {
      dielectronAllList.reserve(1);
      // dimuonAllList.reserve(1);
    }

    constexpr bool eventHasQvector = ((TEventFillMap & VarManager::ObjTypes::CollisionQvect) > 0);
    constexpr bool trackHasCov = ((TTrackFillMap & VarManager::ObjTypes::TrackCov) > 0);

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0))
        continue;

      if (fConfigOptions.collSplitting && event.isEventSelected_bit(2)) {
        continue;
      }

      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event, VarManager::fgValues);

      auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());
      if (groupedAssocs.size() == 0)
        continue;

      for (auto& [a1, a2] : o2::soa::combinations(groupedAssocs, groupedAssocs)) {

        if constexpr (TPairType == VarManager::kDecayToEE) {

          twoTrackFilter = a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isBarrelSelectedPrefilter_raw() & fTrackFilterMask;
          if (!twoTrackFilter)
            continue;

          auto t1 = a1.template track_as<TTracks>();
          auto t2 = a2.template track_as<TTracks>();
          sign1 = t1.sign();
          sign2 = t2.sign();

          // store the ambiguity number of the two dilepton legs in the last 4 digits of the two-track filter
          if (t1.barrelAmbiguityInBunch() > 1) {
            twoTrackFilter |= (static_cast<uint32_t>(1) << 28);
          }
          if (t2.barrelAmbiguityInBunch() > 1) {
            twoTrackFilter |= (static_cast<uint32_t>(1) << 29);
          }
          if (t1.barrelAmbiguityOutOfBunch() > 1) {
            twoTrackFilter |= (static_cast<uint32_t>(1) << 30);
          }
          if (t2.barrelAmbiguityOutOfBunch() > 1) {
            twoTrackFilter |= (static_cast<uint32_t>(1) << 31);
          }

          VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
          // fill tables
          dielectronList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                         VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                         t1.sign() + t2.sign(), twoTrackFilter, -1);

          dielectronInfoList(event.globalIndex(), t1.globalIndex(), t2.globalIndex());
          dileptonInfoList(event.globalIndex(), event.posX(), event.posY(), event.posZ());

          if (fConfigOptions.fPropTrack) {
            VarManager::FillPairCollision<TPairType, TTrackFillMap>(event, t1, t2);
          }

          if constexpr (TTwoProngFitter) {
            VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigOptions.propToPCA);
            o2::dataformats::VertexBase pvRefit;
            if (fConfigOptions.recomputePV) {
              VarManager::SetPVrecalculationKF(false);
              VarManager::ResetValues(VarManager::kVertexingLxyProjectedRecalculatePV, VarManager::kVertexingLxyProjectedRecalculatePV + 1);
              VarManager::ResetValues(VarManager::kVertexingTauxyProjectedPoleJPsiMassRecalculatePV, VarManager::kVertexingTauxyProjectedPoleJPsiMassRecalculatePV + 1);
              // cout << "primary vertex (before): x -> " << event.posX() << " y -> " << event.posY() << " z -> " << event.posZ() << endl;
              o2::dataformats::VertexBase pvRefit;
              bool ok = refitPVWithPVertexer(event, tracks, t1, t2, pvRefit);
              if (ok)
                VarManager::FillPairVertexingRecomputePV<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, pvRefit);
              // cout << "primary vertex (after): ok -> " << ok << " x -> " << pvRefit.getX() << " y -> " << pvRefit.getY() << " z -> " << pvRefit.getZ() << endl;
            }
          }

          if constexpr (trackHasCov && TTwoProngFitter) {
            dielectronsExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
            if (fConfigOptions.flatTables.value) {
              dielectronAllList(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), twoTrackFilter, -1,
                                // t1.pt(), t1.eta(), t1.phi(), t1.itsClusterMap(), t1.itsChi2NCl(), t1.tpcNClsCrossedRows(), t1.tpcNClsFound(), t1.tpcChi2NCl(), t1.dcaXY(), t1.dcaZ(), t1.tpcSignal(), t1.tpcNSigmaEl(), t1.tpcNSigmaPi(), t1.tpcNSigmaPr(), t1.beta(), t1.tofNSigmaEl(), t1.tofNSigmaPi(), t1.tofNSigmaPr(),
                                t1.pt(), t1.eta(), t1.phi(), t1.itsClusterMap(), t1.itsChi2NCl(), t1.tpcNClsCrossedRows(), t1.tpcNClsFound(), t1.tpcChi2NCl(), t1.dcaXY(), t1.dcaZ(), t1.tpcSignal(), t1.tpcNSigmaEl(), t1.tpcNSigmaPi(), t1.tpcNSigmaPr(), -999.0, -999.0, -999.0, -999.0,
                                // t2.pt(), t2.eta(), t2.phi(), t2.itsClusterMap(), t2.itsChi2NCl(), t2.tpcNClsCrossedRows(), t2.tpcNClsFound(), t2.tpcChi2NCl(), t2.dcaXY(), t2.dcaZ(), t2.tpcSignal(), t2.tpcNSigmaEl(), t2.tpcNSigmaPi(), t2.tpcNSigmaPr(), t2.beta(), t2.tofNSigmaEl(), t2.tofNSigmaPi(), t2.tofNSigmaPr(),
                                t2.pt(), t2.eta(), t2.phi(), t2.itsClusterMap(), t2.itsChi2NCl(), t2.tpcNClsCrossedRows(), t2.tpcNClsFound(), t2.tpcChi2NCl(), t2.dcaXY(), t2.dcaZ(), t2.tpcSignal(), t2.tpcNSigmaEl(), t2.tpcNSigmaPi(), t2.tpcNSigmaPr(), -999.0, -999.0, -999.0, -999.0,
                                VarManager::fgValues[VarManager::kKFTrack0DCAxyz], VarManager::fgValues[VarManager::kKFTrack1DCAxyz], VarManager::fgValues[VarManager::kKFDCAxyzBetweenProngs], VarManager::fgValues[VarManager::kKFTrack0DCAxy], VarManager::fgValues[VarManager::kKFTrack1DCAxy], VarManager::fgValues[VarManager::kKFDCAxyBetweenProngs],
                                VarManager::fgValues[VarManager::kKFTrack0DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack0DeviationxyFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationxyFromPV],
                                VarManager::fgValues[VarManager::kKFMass], VarManager::fgValues[VarManager::kKFChi2OverNDFGeo], VarManager::fgValues[VarManager::kVertexingLxyz], VarManager::fgValues[VarManager::kVertexingLxyzOverErr], VarManager::fgValues[VarManager::kVertexingLxy], VarManager::fgValues[VarManager::kVertexingLxyOverErr], VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr], VarManager::fgValues[VarManager::kKFCosPA], VarManager::fgValues[VarManager::kKFJpsiDCAxyz], VarManager::fgValues[VarManager::kKFJpsiDCAxy],
                                VarManager::fgValues[VarManager::kKFPairDeviationFromPV], VarManager::fgValues[VarManager::kKFPairDeviationxyFromPV],
                                VarManager::fgValues[VarManager::kKFMassGeoTop], VarManager::fgValues[VarManager::kKFChi2OverNDFGeoTop],
                                VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingTauxyProjected],
                                VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
            }
          }

          if constexpr (eventHasQvector) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }
        }
        // Fill normal histograms
        bool isAmbiInBunch = (twoTrackFilter & (1 << 28)) || (twoTrackFilter & (1 << 29));
        bool isAmbiOutOfBunch = (twoTrackFilter & (1 << 30)) || (twoTrackFilter & (1 << 31));

        for (int icut = 0; icut < ncuts; icut++) { // loop over cut definitions
          if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
            if (sign1 * sign2 < 0) { // opposite sign pairs
              fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
              PromptNonPromptSepTable(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kRap], VarManager::fgValues[VarManager::kPhi], VarManager::fgValues[VarManager::kVertexingTauxyProjected], VarManager::fgValues[VarManager::kVertexingTauxyProjectedPoleJPsiMass], VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingTauxyProjectedPoleJPsiMassRecalculatePV], isAmbiInBunch, isAmbiOutOfBunch, VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);
              if (fConfigOptions.fConfigQA) {
                if (isAmbiInBunch) {
                  fHistMan->FillHistClass(histNames[icut][3].Data(), VarManager::fgValues);
                }
                if (isAmbiOutOfBunch) {
                  fHistMan->FillHistClass(histNames[icut][3 + 3].Data(), VarManager::fgValues);
                }
              }
            } else if (sign1 > 0) { // ++ pairs
              fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
              if (fConfigOptions.fConfigQA) {
                if (isAmbiInBunch) {
                  fHistMan->FillHistClass(histNames[icut][4].Data(), VarManager::fgValues);
                }
                if (isAmbiOutOfBunch) {
                  fHistMan->FillHistClass(histNames[icut][4 + 3].Data(), VarManager::fgValues);
                }
              }
            } else { // -- pairs
              fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
              if (fConfigOptions.fConfigQA) {
                if (isAmbiInBunch) {
                  fHistMan->FillHistClass(histNames[icut][5].Data(), VarManager::fgValues);
                }
                if (isAmbiOutOfBunch) {
                  fHistMan->FillHistClass(histNames[icut][5 + 3].Data(), VarManager::fgValues);
                }
              }
            }

            // Pair cuts
            for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++) {
              AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
              if (!cut.IsSelected(VarManager::fgValues))
                continue;              // apply pair cuts
              if (sign1 * sign2 < 0) { // opposite sign pairs
                fHistMan->FillHistClass(histNames[ncuts + icut * fPairCuts.size() + iPairCut][0].Data(), VarManager::fgValues);
              } else if (sign1 > 0) { // ++ pairs
                fHistMan->FillHistClass(histNames[ncuts + icut * fPairCuts.size() + iPairCut][1].Data(), VarManager::fgValues);
              } else { // -- pairs
                fHistMan->FillHistClass(histNames[ncuts + icut * fPairCuts.size() + iPairCut][2].Data(), VarManager::fgValues);
              }
            } // end loop (pair cuts)
          }
        } // end loop (cuts)
      } // end loop over pairs of track associations
    } // end loop over events

  }

  // Template function for electron-muon same-event pairing (barrel x muon, full index policy)
  template <bool TTwoProngFitter, int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap,
            typename TEvents, typename TTrackAssocs, typename TTracks, typename TMuonAssocs, typename TMuons>
  void runEmuSameEventPairing(TEvents const& events, BCsWithTimestamps const& bcs,
                              Preslice<TTrackAssocs>& preslice1, TTrackAssocs const& assocs1, TTracks const& /*tracks1*/,
                              Preslice<TMuonAssocs>& preslice2, TMuonAssocs const& assocs2, TMuons const& /*tracks2*/)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    if (fCurrentRun != bcs.begin().runNumber()) {
      initParamsFromCCDB(bcs.begin().timestamp(), TTwoProngFitter);
      fCurrentRun = bcs.begin().runNumber();
    }

    const auto& histNames = fTrackMuonHistNames;
    int nPairCuts = (fPairCuts.size() > 0) ? static_cast<int>(fPairCuts.size()) : 1;

    electronmuonList.reserve(1);

    uint32_t twoTrackFilter = 0;
    int sign1 = 0;
    int sign2 = 0;

    constexpr bool eventHasQvector = ((TEventFillMap & VarManager::ObjTypes::CollisionQvect) > 0);

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0))
        continue;
      if (fConfigOptions.collSplitting && event.isEventSelected_bit(2))
        continue;

      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event, VarManager::fgValues);

      auto groupedAssocs1 = assocs1.sliceBy(preslice1, event.globalIndex());
      if (groupedAssocs1.size() == 0)
        continue;
      auto groupedAssocs2 = assocs2.sliceBy(preslice2, event.globalIndex());
      if (groupedAssocs2.size() == 0)
        continue;

      for (auto& [a1, a2] : o2::soa::combinations(soa::CombinationsFullIndexPolicy(groupedAssocs1, groupedAssocs2))) {
        if (!(a1.isBarrelSelected_raw() & fTrackFilterMask))
          continue;
        if (!(a2.isMuonSelected_raw() & fMuonFilterMask))
          continue;

        auto t1 = a1.template track_as<TTracks>();
        auto t2 = a2.template fwdtrack_as<TMuons>();
        sign1 = t1.sign();
        sign2 = t2.sign();

        twoTrackFilter = 0;
        int minCuts = std::min(fNCutsBarrel, fNCutsMuon);
        for (int i = 0; i < minCuts; ++i) {
          if ((a1.isBarrelSelected_raw() & (1u << i)) && (a2.isMuonSelected_raw() & (1u << i))) {
            twoTrackFilter |= (1u << i);
          }
        }
        // store ambiguity flags in bits 28-31
        if (t1.barrelAmbiguityInBunch() > 1)
          twoTrackFilter |= (1u << 28);
        if (t2.muonAmbiguityInBunch() > 1)
          twoTrackFilter |= (1u << 29);
        if (t1.barrelAmbiguityOutOfBunch() > 1)
          twoTrackFilter |= (1u << 30);
        if (t2.muonAmbiguityOutOfBunch() > 1)
          twoTrackFilter |= (1u << 31);

        VarManager::FillPair<VarManager::kElectronMuon, TTrackFillMap>(t1, t2);
        if (fConfigOptions.fPropTrack) {
          VarManager::FillPairCollision<VarManager::kElectronMuon, TTrackFillMap>(event, t1, t2);
        }
        if constexpr (eventHasQvector) {
          VarManager::FillPairVn<VarManager::kElectronMuon>(t1, t2);
        }

        electronmuonList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                         VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta],
                         VarManager::fgValues[VarManager::kPhi],
                         t1.sign() + t2.sign(), twoTrackFilter, 0);

        for (int iTrack = 0; iTrack < fNCutsBarrel; ++iTrack) {
          if (!(a1.isBarrelSelected_raw() & (1u << iTrack)))
            continue;
          for (int iMuon = 0; iMuon < fNCutsMuon; ++iMuon) {
            if (!(a2.isMuonSelected_raw() & (1u << iMuon)))
              continue;
            for (unsigned int iPairCut = 0; iPairCut < (fPairCuts.empty() ? 1u : static_cast<unsigned int>(fPairCuts.size())); iPairCut++) {
              if (!fPairCuts.empty()) {
                AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
                if (!cut.IsSelected(VarManager::fgValues))
                  continue;
              }
              int index = iTrack * (fNCutsMuon * nPairCuts) + iMuon * nPairCuts + static_cast<int>(iPairCut);
              auto itHist = histNames.find(index);
              if (itHist == histNames.end())
                continue;
              if (sign1 * sign2 < 0) {
                fHistMan->FillHistClass(itHist->second[0].Data(), VarManager::fgValues);
              } else if (sign1 > 0) {
                fHistMan->FillHistClass(itHist->second[1].Data(), VarManager::fgValues);
              } else {
                fHistMan->FillHistClass(itHist->second[2].Data(), VarManager::fgValues);
              }
            } // end pair cut loop
          } // end muon cut loop
        } // end barrel cut loop

      } // end combinations loop
    } // end event loop
  }

  void processBarrelOnly(MyEventsSelected const& events, BCsWithTimestamps const& bcs,
                         soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                         MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithMults, gkTrackFillMapWithCov>(events, bcs, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processElectronMuonDirect(
    MyEventsSelected const& events, BCsWithTimestamps const& bcs,
    soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
    MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
    soa::Join<aod::FwdTrackAssoc, aod::MuonTrackCuts> const& muonAssocs,
    MyMuonTracksWithCovWithAmbiguities const& muons)
  {
    runEmuSameEventPairing<true, VarManager::kElectronMuon,
                           gkEventFillMapWithMults, gkTrackFillMapWithCov, gkMuonFillMapWithCov>(
      events, bcs,
      trackEmuAssocsPerCollision, barrelAssocs, barrelTracks,
      muonAssocsPerCollision, muonAssocs, muons);
  }

  void processDummy(MyEvents&) { /* do nothing */ }

  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnly, "Run barrel only pairing", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processElectronMuonDirect, "Run electron-muon pairing on AO2D tracks/fwd-tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Initialize metadata for TOF response
  o2::pid::tof::TOFResponseImpl::metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisMuonSelection>(cfgc),
    adaptAnalysisTask<AnalysisPrefilterSelection>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc)};
  // adaptAnalysisTask<AnalysisDileptonTrack>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //  The histogram classes are provided in the histClasses string, separated by semicolon ";"
  //  The histogram classes and their components histograms are defined below depending on the name of the histogram class
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    TString histName = histGroups;
    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("TimeFrameStats")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "timeframe");
    }
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", histName);
    }

    if (classStr.Contains("SameBunchCorrelations") || classStr.Contains("OutOfBunchCorrelations")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "two-collisions", histName);
    }

    if ((classStr.Contains("Track") || classStr.Contains("Assoc")) && !classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
        if (classStr.Contains("PIDCalibElectron")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_electron");
        }
        if (classStr.Contains("PIDCalibPion")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_pion");
        }
        if (classStr.Contains("PIDCalibProton")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_proton");
        }
        if (classStr.Contains("Ambiguity")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "ambiguity");
        }
      }
    }
    if (classStr.Contains("Muon") && !classStr.Contains("Pairs")) {
      if (!classStr.Contains("Ambiguity")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
      } else {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon-ambiguity");
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("Triplets")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "barrel,vertexing");
    }

    if (classStr.Contains("DileptonTrack") && !classStr.Contains("ME")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-track", histName);
    }

    if (classStr.Contains("DileptonTrackME")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-track", "mixedevent");
    }

    if (classStr.Contains("HadronsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
    }

    if (classStr.Contains("DileptonHadronInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-mass");
    }

    if (classStr.Contains("DileptonHadronCorrelation")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-correlation");
    }

    if (classStr.Contains("MCTruthEenergyCorrelators")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "energy-correlator-gen");
    }
  } // end loop over histogram classes
}

/*
struct AnalysisDileptonTrack {
  int fCurrentRun = -1;

  // Preslice per associazioni e dileptoni
  Preslice<soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts>> trackAssocsPerCollision = aod::track_association::collisionId;
  Preslice<MyDielectronCandidates> dielectronsPerCollision = aod::reducedpair::reducedeventId;

  // Configurazioni e istogrammi
  // ConfigOptions fConfigOptions;
  struct : ConfigurableGroup {
    Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "kaonPID", "Comma separated list of track cuts to be correlated with the dileptons"};
    Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 2.8, "Low mass cut for the dileptons used in analysis"};
    Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 3.2, "High mass cut for the dileptons used in analysis"};
    Configurable<float> fConfigDileptonLxyCut{"cfgDileptonLxyCut", 0.0, "Lxy cut for dileptons used in the triplet vertexing"};
    Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
    Configurable<float> fConfigDileptonLowpTCut{"cfgDileptonLowpTCut", 0.0, "Low pT cut for dileptons used in the triplet vertexing"};
    Configurable<float> fConfigDileptonHighpTCut{"cfgDileptonHighpTCut", 1E5, "High pT cut for dileptons used in the triplet vertexing"};
    Configurable<float> fConfigDileptonRapCutAbs{"cfgDileptonRapCutAbs", 1.0, "Rap cut for dileptons used in the triplet vertexing"};
    Configurable<std::string> fConfigHistogramSubgroups{"cfgDileptonTrackHistogramsSubgroups", "invmass,vertexing", "Comma separated list of dilepton-track histogram subgroups"};
    Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
    Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 5, "Event mixing pool depth"};
    Configurable<bool> fConfigPublishTripletTable{"cfgPublishTripletTable", false, "Publish the triplet tables, BmesonCandidates"};
    Configurable<bool> fConfigApplyMassEC{"cfgApplyMassEC", false, "Apply fit mass for sideband for the energy correlator study"};
    Configurable<std::vector<int>> fConfigSavelessevents{"cfgSavelessevents", std::vector<int>{1, 0}, "Save less events for the energy correlator study"};
  } fConfigOptions;

  struct : ConfigurableGroup {
    Configurable<bool> fConfigUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
    Configurable<std::string> fConfigGRPmagPath{"cfgGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};
    Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
    Configurable<std::string> fConfigGeoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } fConfigCCDBOptions;

  HistogramManager* fHistMan = nullptr;

  void processBarrel(soa::Filtered<MyEventsSelected> const& events, BCsWithTimestamps const& bcs,
                     soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts> const& assocs,
                     MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {

    if (events.size() == 0) return;

    if (fCurrentRun != bcs.begin().runNumber()) {
      initParamsFromCCDB(bcs.begin().timestamp());
      fCurrentRun = bcs.begin().runNumber();
    }

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) continue;

      std::vector<int> fSavelessevents = fConfigOptions.fConfigSavelessevents;
      if (fSavelessevents[0] > 1 && event.globalIndex() % fSavelessevents[0] == fSavelessevents[1]) continue;

      auto groupedBarrelAssocs = assocs.sliceBy(trackAssocsPerCollision, event.globalIndex());
      auto groupedDielectrons = dileptons.sliceBy(dielectronsPerCollision, event.globalIndex());

      runDileptonHadron<VarManager::kBtoJpsiEEK, gkEventFillMapWithMults, gkTrackFillMapWithCov>(
          event, bcs, groupedBarrelAssocs, tracks, groupedDielectrons);
    }

  }

  void processDummy(MyEvents&) {
    // funzione dummy, non fa nulla
  }

  void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups) {
    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      histMan->AddHistClass(classStr.Data());

      if (classStr.Contains("TimeFrameStats")) {
        dqhistograms::DefineHistograms(histMan, classStr, "timeframe");
      }
      if (classStr.Contains("Event")) {
        dqhistograms::DefineHistograms(histMan, classStr, "event", histGroups);
      }
      if (classStr.Contains("Track") && !classStr.Contains("Pairs")) {
        if (classStr.Contains("Barrel")) {
          dqhistograms::DefineHistograms(histMan, classStr, "track", histGroups);
        }
      }
      if (classStr.Contains("Pairs")) {
        dqhistograms::DefineHistograms(histMan, classStr, "pair", histGroups);
      }
      if (classStr.Contains("DileptonTrack") && !classStr.Contains("ME")) {
        dqhistograms::DefineHistograms(histMan, classStr, "dilepton-track", histGroups);
      }
      if (classStr.Contains("DileptonTrackME")) {
        dqhistograms::DefineHistograms(histMan, classStr, "dilepton-track", "mixedevent");
      }
      if (classStr.Contains("DileptonHadronInvMass")) {
        dqhistograms::DefineHistograms(histMan, classStr, "dilepton-hadron-mass");
      }
      if (classStr.Contains("DileptonHadronCorrelation")) {
        dqhistograms::DefineHistograms(histMan, classStr, "dilepton-hadron-correlation");
      }
    }
  }

  PROCESS_SWITCH(AnalysisDileptonTrack, processBarrel, "Run barrel dilepton-track pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processDummy, "Dummy function", true);
}; */
