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

#include "Common/Core/PID/PIDTOFParamService.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <Common/Core/trackUtilities.h>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <DetectorsVertexing/PVertexer.h>
#include <ReconstructionDataFormats/Track.h>

#include "TGeoGlobalMagField.h"
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TObjString.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

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
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASHA", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected);
DECLARE_SOA_TABLE(BarrelAmbiguities, "AOD", "DQBARRELAMB", dqanalysisflags::BarrelAmbiguityInBunch, dqanalysisflags::BarrelAmbiguityOutOfBunch);
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsBarrelSelectedPrefilter);

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

using MyDielectronCandidates = soa::Join<aod::Dielectrons, aod::DielectronsExtra>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMapWithMults = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionMult | VarManager::ObjTypes::CollisionMultExtra;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkTrackFillMapWithCovNoTOF = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackTPCPID | VarManager::ObjTypes::TrackTOFService;
constexpr static uint32_t gkDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups); // defines histograms for all tasks

template <typename TMap>
void PrintBitMap(TMap map, int nbits)
{
  for (int i = 0; i < nbits; i++) {
    cout << ((map & (TMap(1) << i)) > 0 ? "1" : "0");
  }
}

struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigEventCutsJSON{"cfgEventCutsJSON", "", "Additional event cuts specified in JSON format"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Add event histograms defined via JSON formatting (see HistogramsLibrary)"};

  Configurable<float> fConfigSplitCollisionsDeltaZ{"cfgSplitCollisionsDeltaZ", 1.0, "maximum delta-z (cm) between two collisions to consider them as split candidates"};
  Configurable<unsigned int> fConfigSplitCollisionsDeltaBC{"cfgSplitCollisionsDeltaBC", 100, "maximum delta-BC between two collisions to consider them as split candidates; do not apply if value is negative"};
  Configurable<bool> fConfigCheckSplitCollisions{"cfgCheckSplitCollisions", false, "If true, run the split collision check and fill histograms"};

  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;

  AnalysisCompositeCut* fEventCut;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  std::map<int64_t, bool> fSelMap;                     // key: reduced event global index, value: event selection decision
  std::map<uint64_t, std::vector<int64_t>> fBCCollMap; // key: global BC, value: vector of reduced event global indices
  int fCurrentRun;

  void init(o2::framework::InitContext& context)
  {
    cout << "AnalysisEventSelection::init() called" << endl;
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

    TString mixVarsString = fConfigMixingVariables.value;
    std::unique_ptr<TObjArray> objArray(mixVarsString.Tokenize(","));
    if (objArray->GetEntries() > 0) {
      fMixHandler = new MixingHandler("mixingHandler", "mixing handler");
      fMixHandler->Init();
      for (int iVar = 0; iVar < objArray->GetEntries(); ++iVar) {
        dqmixing::SetUpMixing(fMixHandler, objArray->At(iVar)->GetName());
      }
    }

    fCurrentRun = -1;
    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    fCCDBApi.init(fConfigCcdbUrl.value);
    cout << "AnalysisEventSelection::init() completed" << endl;
  }

  template <uint32_t TEventFillMap, typename TEvents>
  void runEventSelection(TEvents const& events, BCsWithTimestamps const& bcs)
  {
    cout << "AnalysisEventSelection::runEventSelection() called with " << events.size() << " events and " << bcs.size() << " BCs" << endl;
    if (bcs.size() > 0 && bcs.begin().runNumber() != fCurrentRun) {
      std::map<std::string, std::string> metadataRCT, header;
      header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", bcs.begin().runNumber()), metadataRCT, -1);
      uint64_t sor = std::atol(header["SOR"].c_str());
      uint64_t eor = std::atol(header["EOR"].c_str());
      VarManager::SetSORandEOR(sor, eor);
    }

    cout << "Filling TimeFrame statistics histograms" << endl;
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillTimeFrame(bcs);
    VarManager::FillTimeFrame(events);
    if (fConfigQA) {
      fHistMan->FillHistClass("TimeFrameStats", VarManager::fgValues);
    }

    fSelMap.clear();
    fBCCollMap.clear();

    cout << "Starting event loop for event selection" << endl;
    for (auto& event : events) {
      auto bc = event.template bc_as<BCsWithTimestamps>();

      VarManager::ResetValues(VarManager::kNTFWiseVariables, VarManager::kNEventWiseVariables);
      VarManager::FillBC(bc);
      VarManager::FillEvent<TEventFillMap>(event);

      bool decision = false;
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }
      if (fEventCut->IsSelected(VarManager::fgValues)) {
        if (fConfigQA) {
          fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
        }
        decision = true;
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

    cout << "AnalysisEventSelection::runEventSelection() completed" << endl;
  }

  template <uint32_t TEventFillMap, typename TEvents>
  void publishSelections(TEvents const& events)
  {
    cout << "AnalysisEventSelection::publishSelections() called" << endl;
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
    cout << "AnalysisEventSelection::publishSelections() completed" << endl;
  }

  void processDirect(MyEvents const& events, BCsWithTimestamps const& bcs)
  {
    cout << "AnalysisEventSelection::processDirect() called" << endl;
    runEventSelection<gkEventFillMapWithMults>(events, bcs);
    publishSelections<gkEventFillMapWithMults>(events);
    cout << "AnalysisEventSelection::processDirect() completed" << endl;
  }

  void processDummy(aod::Collisions&) {}

  PROCESS_SWITCH(AnalysisEventSelection, processDirect, "Run event selection on framework AO2Ds", false);
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
    cout << "AnalysisTrackSelection::init() called" << endl;
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

    cout << "AnalysisTrackSelection::init() completed" << endl;
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTracks>
  void runTrackSelection(TrackAssoc const& assocs, BCsWithTimestamps const& bcs, TEvents const& events, TTracks const& tracks)
  {
    cout << "AnalysisTrackSelection::runTrackSelection() called" << endl;
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

    cout << "Starting loop over track associations" << endl;

    for (auto& assoc : assocs) {
      auto event = assoc.template collision_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        trackSel(0);
        continue;
      }

      VarManager::ResetValues(VarManager::kNTFWiseVariables, VarManager::kNBarrelTrackVariables);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);

      auto track = tracks.rawIteratorAt(assoc.trackId());
      VarManager::FillTrack<TTrackFillMap>(track);
      // compute quantities which depend on the associated collision, such as DCA
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

    cout << "AnalysisTrackSelection::runTrackSelection() completed" << endl;
  }

  void processWithCov(TrackAssoc const& assocs, BCsWithTimestamps const& bcs, MyEventsSelected const& events, MyBarrelTracksWithCov const& tracks)
  {
    cout << "AnalysisTrackSelection::processWithCov() called" << endl;
    runTrackSelection<gkEventFillMapWithMults, gkTrackFillMapWithCov>(assocs, bcs, events, tracks);
    cout << "AnalysisTrackSelection::processWithCov() completed" << endl;
  }

  void processWithCovTOFService(TrackAssoc const& assocs, BCsWithTimestamps const& bcs, MyEventsSelected const& events, MyBarrelTracksWithCovNoTOF const& tracks)
  {
    cout << "AnalysisTrackSelection::processWithCovTOFService() called" << endl;
    fTofResponse->processSetup(bcs.iteratorAt(0));
    auto tracksWithTOFservice = soa::Attach<MyBarrelTracksWithCovNoTOF, o2::aod::TOFNSigmaDynEl, o2::aod::TOFNSigmaDynPi,
                                            o2::aod::TOFNSigmaDynKa, o2::aod::TOFNSigmaDynPr>(tracks);
    runTrackSelection<gkEventFillMapWithMults, gkTrackFillMapWithCovNoTOF>(assocs, bcs, events, tracksWithTOFservice);
    cout << "AnalysisTrackSelection::processWithCovTOFService() completed" << endl;
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
    cout << "AnalysisPrefilterSelection::init() called" << endl;
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
    cout << "AnalysisPrefilterSelection::init() completed" << endl;
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
    cout << "AnalysisPrefilterSelection::processBarrel() called" << endl;
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
    cout << "AnalysisPrefilterSelection::processBarrel() completed" << endl;
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisPrefilterSelection, processBarrel, "Run Prefilter selection on barrel tracks", false);
  PROCESS_SWITCH(AnalysisPrefilterSelection, processDummy, "Do nothing", true);
};

struct AnalysisSameEventPairing {
  Produces<aod::Dielectrons> dielectronList;
  Produces<aod::DielectronsExtra> dielectronsExtraList;
  Produces<aod::DielectronsInfo> dielectronInfoList;
  Produces<aod::DielectronsAll> dielectronAllList;
  Produces<aod::DileptonInfo> dileptonInfoList;
  Produces<aod::JPsieeCandidates> PromptNonPromptSepTable;

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

  Preslice<soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts, aod::Prefilter>> trackAssocsPerCollision = aod::track_association::collisionId;
  // Preslice<soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

  void init(o2::framework::InitContext& context)
  {
    cout << "AnalysisSameEventPairing::init() called" << endl;
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    fEnableBarrelHistos = context.mOptions.get<bool>("processBarrelOnly");
    // fEnableMuonHistos = context.mOptions.get<bool>("processMuonOnlySkimmed");

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
    /*TString muonCutsStr = fConfigOptions.muon.value;
    TObjArray* objArrayMuonCuts = nullptr;
    if (!muonCutsStr.IsNull()) {
      objArrayMuonCuts = muonCutsStr.Tokenize(",");
    }*/

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

    /*
    // get the muon track selection cuts
    getTaskOptionValue<string>(context, "analysis-muon-selection", "cfgMuonCuts", tempCuts, false);
    tempCutsStr = tempCuts;
    // check also the cuts added via JSON and add them to the string of cuts
    getTaskOptionValue<string>(context, "analysis-muon-selection", "cfgMuonCutsJSON", tempCuts, false);
    TString addMuonCutsStr = tempCuts;
    if (addMuonCutsStr != "") {
      std::vector<AnalysisCut*> addMuonCuts = dqcuts::GetCutsFromJSON(addMuonCutsStr.Data());
      for (auto& t : addMuonCuts) {
        tempCutsStr += Form(",%s", t->GetName());
      }
    }

    // check that in this task we have specified muon cuts
    if (!muonCutsStr.IsNull()) {
      // loop over the muon cuts computed by the muon selection task and build a filter mask for those required in this task
      std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
      fNCutsMuon = objArray->GetEntries();
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        if (objArrayMuonCuts->FindObject(tempStr.Data()) != nullptr) {
          // update the filter mask
          fMuonFilterMask |= (static_cast<uint32_t>(1) << icut);

          if (fEnableMuonHistos) {
            // assign pair hist directories for each required muon cut
            std::vector<TString> names = {
              Form("PairsMuonSEPM_%s", objArray->At(icut)->GetName()),
              Form("PairsMuonSEPP_%s", objArray->At(icut)->GetName()),
              Form("PairsMuonSEMM_%s", objArray->At(icut)->GetName())};
            if (fConfigOptions.fConfigQA) {
              // assign separate hist directories for ambiguous tracks
              names.push_back(Form("PairsMuonSEPM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEPP_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEMM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEPM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEPP_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEMM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
            }
            for (auto& n : names) {
              histNames += Form("%s;", n.Data());
            }
            fMuonHistNames[icut] = names;

            // if there are specified pair cuts, assign hist dirs for each muon cut - pair cut combination
            TString cutNamesStr = fConfigOptions.pair.value;
            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              fNPairCuts = objArrayPair->GetEntries();
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                names = {
                  Form("PairsMuonSEPM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsMuonSEPP_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsMuonSEMM_%s_%s", objArray->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
                histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
                fMuonHistNames[fNCutsMuon + icut * fNCutsMuon + iPairCut] = names;
              } // end loop (pair cuts)
            } // end if (pair cuts)

            // assign hist directories for pairs matched to MC signals for each (muon cut, MCrec signal) combination
            if (!sigNamesStr.IsNull()) {
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
                auto sig = fRecMCSignals.at(isig);
                names = {
                  Form("PairsMuonSEPM_%s_%s", objArray->At(icut)->GetName(), sig->GetName()),
                  Form("PairsMuonSEPP_%s_%s", objArray->At(icut)->GetName(), sig->GetName()),
                  Form("PairsMuonSEMM_%s_%s", objArray->At(icut)->GetName(), sig->GetName()),
                };
                if (fConfigOptions.fConfigQA) {
                  names.push_back(Form("PairsMuonSEPMCorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsMuonSEPMIncorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsMuonSEPM_ambiguousInBunch_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsMuonSEPM_ambiguousInBunchCorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsMuonSEPM_ambiguousInBunchIncorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsMuonSEPM_ambiguousOutOfBunch_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsMuonSEPM_ambiguousOutOfBunchCorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsMuonSEPM_ambiguousOutOfBunchIncorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                }
                for (auto& n : names) {
                  histNames += Form("%s;", n.Data());
                }
                fMuonHistNamesMCmatched.try_emplace(icut * fRecMCSignals.size() + isig, names);
              } // end loop over MC signals
            }
          }
        }
      } // end loop over cuts
    } // end if (muonCutsStr)
*/

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

    cout << "AnalysisSameEventPairing::init() completed" << endl;
  }

  void initParamsFromCCDB(uint64_t timestamp, bool withTwoProngFitter = true)
  {
    cout << "AnalysisSameEventPairing::initParamsFromCCDB() called for timestamp " << timestamp << endl;
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
    cout << "AnalysisSameEventPairing::initParamsFromCCDB() completed" << endl;
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
    cout << "AnalysisSameEventPairing::runSameEventPairing() called" << endl;
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

    cout << "AnalysisSameEventPairing::runSameEventPairing() completed" << endl;
  }

  void processBarrelOnly(MyEventsSelected const& events, BCsWithTimestamps const& bcs,
                         soa::Join<aod::TrackAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                         MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    cout << "AnalysisSameEventPairing::processBarrelOnly() called" << endl;
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithMults, gkTrackFillMapWithCov>(events, bcs, trackAssocsPerCollision, barrelAssocs, barrelTracks);
    cout << "AnalysisSameEventPairing::processBarrelOnly() completed" << endl;
  }

  void processDummy(MyEvents&) { /* do nothing */ }

  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnly, "Run barrel only pairing", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Initialize metadata for TOF response
  o2::pid::tof::TOFResponseImpl::metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
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
    std::cout << "AnalysisDileptonTrack::processBarrel() called" << std::endl;

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

    std::cout << "AnalysisDileptonTrack::processBarrel() completed" << std::endl;
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
