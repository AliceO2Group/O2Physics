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
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/Core/TableHelper.h"

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
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);                                     //! Hash used in event mixing //need to understand
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 32);                     //! Event decision
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelected, isBarrelSelected, 32);                   //! Barrel track decisions
DECLARE_SOA_COLUMN(BarrelAmbiguityInBunch, barrelAmbiguityInBunch, int8_t);          //! Barrel track in-bunch ambiguity
DECLARE_SOA_COLUMN(BarrelAmbiguityOutOfBunch, barrelAmbiguityOutOfBunch, int8_t);    //! Barrel track out of bunch ambiguity
DECLARE_SOA_BITMAP_COLUMN(IsMuonSelected, isMuonSelected, 32);                       //! Muon track decisions (joinable to ReducedMuonsAssoc)
DECLARE_SOA_COLUMN(MuonAmbiguityInBunch, muonAmbiguityInBunch, int8_t);              //! Muon track in-bunch ambiguity
DECLARE_SOA_COLUMN(MuonAmbiguityOutOfBunch, muonAmbiguityOutOfBunch, int8_t);        //! Muon track out of bunch ambiguity
DECLARE_SOA_BITMAP_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, 32); //! Barrel prefilter decisions (joinable to ReducedTracksAssoc)
// Bcandidate columns for ML analysis of B->Jpsi+K
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
DECLARE_SOA_COLUMN(MCLxyBcandidate, MClxyBcandidate, float);
DECLARE_SOA_COLUMN(MCLxyzBcandidate, MClxyzBcandidate, float);
DECLARE_SOA_COLUMN(MCLzBcandidate, MClzBcandidate, float);
DECLARE_SOA_COLUMN(MCTauxyBcandidate, MCtauxyBcandidate, float);
DECLARE_SOA_COLUMN(MCTauzBcandidate, MCtauzBcandidate, float);
DECLARE_SOA_COLUMN(CosPBcandidate, cosPBcandidate, float);
DECLARE_SOA_COLUMN(MCCosPBcandidate, MCcosPBcandidate, float);
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
DECLARE_SOA_COLUMN(McFlag, mcFlag, int8_t);
DECLARE_SOA_BITMAP_COLUMN(IsJpsiFromBSelected, isJpsiFromBSelected, 32);
// Candidate columns for prompt-non-prompt JPsi separation
DECLARE_SOA_COLUMN(Massee, massee, float);
DECLARE_SOA_COLUMN(Etaee, etaee, float);
DECLARE_SOA_COLUMN(Rapee, rapee, float);
DECLARE_SOA_COLUMN(Phiee, phiee, float);
DECLARE_SOA_COLUMN(Ptee, ptee, float);
DECLARE_SOA_COLUMN(Lxyee, lxyee, float);
DECLARE_SOA_COLUMN(LxyeePoleMass, lxyeepolemass, float);
DECLARE_SOA_COLUMN(Lzee, lzee, float);
DECLARE_SOA_COLUMN(MultiplicityFT0A, multiplicityFT0AJPsi2ee, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0CJPsi2ee, float);
DECLARE_SOA_COLUMN(PercentileFT0M, percentileFT0MJPsi2ee, float);
DECLARE_SOA_COLUMN(MultiplicityNContrib, multiplicityNContribJPsi2ee, float);
DECLARE_SOA_COLUMN(AmbiguousInBunchPairs, AmbiguousJpsiPairsInBunch, bool);
DECLARE_SOA_COLUMN(AmbiguousOutOfBunchPairs, AmbiguousJpsiPairsOutOfBunch, bool);
DECLARE_SOA_COLUMN(Corrassoc, corrassoc, bool);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
// Candidate columns efficiency calculation for prompt-non-prompt JPsi separation
DECLARE_SOA_COLUMN(OniaPt, oniaPt, float);
DECLARE_SOA_COLUMN(OniaY, oniaY, float);
DECLARE_SOA_COLUMN(OniaEta, oniaEta, float);
DECLARE_SOA_COLUMN(OniaPhi, oniaPhi, float);
DECLARE_SOA_COLUMN(OniaVz, oniaVz, float);
DECLARE_SOA_COLUMN(OniaVtxZ, oniaVtxZ, float);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);                                                            //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASHA", dqanalysisflags::MixingHash);                                                            //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected);                                                    //!  joinable to ReducedTracksAssoc
DECLARE_SOA_TABLE(BarrelAmbiguities, "AOD", "DQBARRELAMB", dqanalysisflags::BarrelAmbiguityInBunch, dqanalysisflags::BarrelAmbiguityOutOfBunch); //!  joinable to ReducedBarrelTracks
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);                                                       //!  joinable to ReducedMuonsAssoc
DECLARE_SOA_TABLE(MuonAmbiguities, "AOD", "DQMUONAMB", dqanalysisflags::MuonAmbiguityInBunch, dqanalysisflags::MuonAmbiguityOutOfBunch);         //!  joinable to ReducedMuonTracks
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsBarrelSelectedPrefilter);                                                  //!  joinable to ReducedTracksAssoc
DECLARE_SOA_TABLE(BmesonCandidates, "AOD", "DQBMESONS",
                  dqanalysisflags::RunNumber, dqanalysisflags::EventIdx, dqanalysisflags::EventTimestamp,
                  dqanalysisflags::massBcandidate, dqanalysisflags::MassDileptonCandidate, dqanalysisflags::deltaMassBcandidate, dqanalysisflags::pTBcandidate, dqanalysisflags::EtaBcandidate, dqanalysisflags::PhiBcandidate, dqanalysisflags::RapBcandidate,
                  dqanalysisflags::LxyBcandidate, dqanalysisflags::LxyBcandidateErr, dqanalysisflags::LxyzBcandidate, dqanalysisflags::LxyzBcandidateErr, dqanalysisflags::LzBcandidate, dqanalysisflags::LzBcandidateErr,
                  dqanalysisflags::TauxyBcandidate, dqanalysisflags::TauxyBcandidateErr, dqanalysisflags::TauzBcandidate, dqanalysisflags::TauzBcandidateErr, dqanalysisflags::CosPBcandidate, dqanalysisflags::Chi2Bcandidate,
                  dqanalysisflags::MCLxyBcandidate, dqanalysisflags::MCLxyzBcandidate, dqanalysisflags::MCLzBcandidate,
                  dqanalysisflags::MCTauxyBcandidate, dqanalysisflags::MCTauzBcandidate, dqanalysisflags::MCCosPBcandidate,
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
                  dqanalysisflags::IsJpsiFromBSelected, dqanalysisflags::IsBarrelSelected, dqanalysisflags::McFlag);
DECLARE_SOA_TABLE(JPsiMuonCandidates, "AOD", "DQJPSIMUONA",
                  dqanalysisflags::DeltaEta, dqanalysisflags::DeltaPhi,
                  dqanalysisflags::MassDileptonCandidate, dqanalysisflags::Ptpair, dqanalysisflags::Etapair, dqanalysisflags::Ptassoc, dqanalysisflags::Etaassoc, dqanalysisflags::Phiassoc,
                  dqanalysisflags::Ptleg1, dqanalysisflags::Etaleg1, dqanalysisflags::Phileg1, dqanalysisflags::Ptleg2, dqanalysisflags::Etaleg2, dqanalysisflags::Phileg2,
                  dqanalysisflags::McFlag);
DECLARE_SOA_TABLE(JPsieeCandidates, "AOD", "DQPSEUDOPROPER", dqanalysisflags::Massee, dqanalysisflags::Ptee, dqanalysisflags::Etaee, dqanalysisflags::Rapee, dqanalysisflags::Phiee, dqanalysisflags::Lxyee, dqanalysisflags::LxyeePoleMass, dqanalysisflags::Lzee, dqanalysisflags::AmbiguousInBunchPairs, dqanalysisflags::AmbiguousOutOfBunchPairs, dqanalysisflags::Corrassoc, dqanalysisflags::MultiplicityFT0A, dqanalysisflags::MultiplicityFT0C, dqanalysisflags::PercentileFT0M, dqanalysisflags::MultiplicityNContrib);
DECLARE_SOA_TABLE(OniaMCTruth, "AOD", "MCTRUTHONIA", dqanalysisflags::OniaPt, dqanalysisflags::OniaEta, dqanalysisflags::OniaY, dqanalysisflags::OniaPhi, dqanalysisflags::OniaVz, dqanalysisflags::OniaVtxZ, dqanalysisflags::MultiplicityFT0A, dqanalysisflags::MultiplicityFT0C, dqanalysisflags::PercentileFT0M, dqanalysisflags::MultiplicityNContrib);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedMCEventLabels>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::ReducedMCEventLabels>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::ReducedMCEventLabels>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedMCEventLabels>;
using MyEventsVtxCovSelectedMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll, aod::ReducedMCEventLabels>;
using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector>;
using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
using MyEventsVtxCovHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedMCEventLabels, aod::MixingHashes>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksWithAmbiguities = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksWithCovWithAmbiguities = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities, aod::ReducedTracksBarrelLabels>;
using MyBarrelTracksWithCovWithAmbiguitiesWithColl = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities, aod::ReducedTracksBarrelLabels, aod::ReducedTracksBarrelInfo>;
using MyDielectronCandidates = soa::Join<aod::Dielectrons, aod::DielectronsExtra>;
using MyDitrackCandidates = soa::Join<aod::Ditracks, aod::DitracksExtra>;
using MyDimuonCandidates = soa::Join<aod::Dimuons, aod::DimuonsExtra>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsLabels>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::ReducedMuonsLabels>;
using MyMuonTracksWithCovWithAmbiguities = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonAmbiguities, aod::ReducedMuonsLabels>;
// using MyMuonTracksSelectedWithColl = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
constexpr static uint32_t gkEventFillMapWithCovMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventMultExtra;
// constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
// constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCovWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;
// constexpr static uint32_t gkTrackFillMapWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;

constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;
// constexpr static uint32_t gkMuonFillMapWithColl = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCollInfo;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups); // defines histograms for all tasks

template <typename TMap>
void PrintBitMap(TMap map, int nbits)
{
  for (int i = 0; i < nbits; i++) {
    cout << ((map & (TMap(1) << i)) > 0 ? "1" : "0");
  }
}

// Analysis task that produces event decisions and the Hash table used in event mixing
struct AnalysisEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  Configurable<std::string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigEventCutsJSON{"cfgEventCutsJSON", "", "Additional event cuts specified in JSON format"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddEventMCHistogram{"cfgAddEventMCHistogram", "generator", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Add event histograms defined via JSON formatting (see HistogramsLibrary)"};

  Configurable<float> fConfigSplitCollisionsDeltaZ{"splitCollisionsDeltaZ", 1.0, "maximum delta-z (cm) between two collisions to consider them as split candidates"};
  Configurable<unsigned int> fConfigSplitCollisionsDeltaBC{"splitCollisionsDeltaBC", 100, "maximum delta-BC between two collisions to consider them as split candidates; do not apply if value is negative"};
  Configurable<bool> fConfigCheckSplitCollisions{"checkSplitCollisions", false, "If true, run the split collision check and fill histograms"};

  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;

  AnalysisCompositeCut* fEventCut;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  std::map<int64_t, bool> fSelMap;                     // key: reduced event global index, value: event selection decision
  std::map<uint64_t, std::vector<int64_t>> fBCCollMap; // key: global BC, value: vector of reduced event global indices
  std::map<string, string> fMetadataRCT, fHeader;
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
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;", fConfigAddEventHistogram.value.data());
      if (fConfigCheckSplitCollisions) {
        DefineHistograms(fHistMan, "OutOfBunchCorrelations;SameBunchCorrelations;", "");
      }
      DefineHistograms(fHistMan, "EventsMC", fConfigAddEventMCHistogram.value.data());
      dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // aditional histograms via JSON
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
  }

  template <uint32_t TEventFillMap, typename TEvents, typename TEventsMC>
  void runEventSelection(TEvents const& events, TEventsMC const& mcEvents)
  {
    if (events.size() > 0 && events.begin().runNumber() != fCurrentRun) {
      fHeader = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", events.begin().runNumber()), fMetadataRCT, -1);
      uint64_t sor = std::atol(fHeader["SOR"].c_str());
      uint64_t eor = std::atol(fHeader["EOR"].c_str());
      VarManager::SetSORandEOR(sor, eor);
    }

    fSelMap.clear();
    fBCCollMap.clear();

    for (auto& event : events) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<TEventFillMap>(event);
      if (event.has_reducedMCevent()) {
        VarManager::FillEvent<VarManager::ObjTypes::ReducedEventMC>(event.reducedMCevent());
      }

      bool decision = false;
      // if QA is requested fill histograms before event selections
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues); // automatically fill all the histograms in the class Event
      }
      if (fEventCut->IsSelected(VarManager::fgValues)) {
        if (fConfigQA) {
          fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
        }
        decision = true;
      }
      fSelMap[event.globalIndex()] = decision;
      if (fBCCollMap.find(event.globalBC()) == fBCCollMap.end()) {
        std::vector<int64_t> evIndices = {event.globalIndex()};
        fBCCollMap[event.globalBC()] = evIndices;
      } else {
        auto& evIndices = fBCCollMap[event.globalBC()];
        evIndices.push_back(event.globalIndex());
      }
      if (fMixHandler != nullptr) {
        int hh = fMixHandler->FindEventCategory(VarManager::fgValues);
        hash(hh);
      }
    }

    for (auto& event : mcEvents) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<VarManager::ObjTypes::ReducedEventMC>(event);
      if (fConfigQA) {
        fHistMan->FillHistClass("EventsMC", VarManager::fgValues);
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
            fHistMan->FillHistClass("SameBunchCorrelations", VarManager::fgValues);
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
            fHistMan->FillHistClass("OutOfBunchCorrelations", VarManager::fgValues);
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
      std::vector<int64_t> sameBunchEvents = fBCCollMap[event.globalBC()];
      if (sameBunchEvents.size() > 1) { // event with in-bunch pileup
        evSel |= (static_cast<uint32_t>(1) << 1);
      }
      if (collisionSplittingMap.find(event.globalIndex()) != collisionSplittingMap.end()) { // event with possible fake in-bunch pileup (collision splitting)
        evSel |= (static_cast<uint32_t>(1) << 2);
      }
      eventSel(evSel);
    }
  }

  void processSkimmed(MyEvents const& events, aod::ReducedMCEvents const& mcEvents)
  {
    runEventSelection<gkEventFillMap>(events, mcEvents);
    publishSelections<gkEventFillMap>(events);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processDummy, "Dummy function", true);
};

// Produces a table with barrel track decisions (joinable to the ReducedTracksAssociations)
// Here one should add all the track cuts needed through the workflow (e.g. cuts for same-even pairing, electron prefiltering, track for dilepton-track correlations)
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

  Configurable<std::string> fConfigMCSignals{"cfgTrackMCSignals", "", "Comma separated list of MC signals"};
  Configurable<std::string> fConfigMCSignalsJSON{"cfgTrackMCsignalsJSON", "", "Additional list of MC signals via JSON"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut*> fTrackCuts;
  std::vector<MCSignal*> fMCSignals; // list of signals to be checked
  std::vector<TString> fHistNamesReco;
  std::vector<TString> fHistNamesMCMatched;

  int fCurrentRun; // current run (needed to detect run changes for loading CCDB parameters)

  std::map<int64_t, std::vector<int64_t>> fNAssocsInBunch;    // key: track global index, value: vector of global index for events associated in-bunch (events that have in-bunch pileup or splitting)
  std::map<int64_t, std::vector<int64_t>> fNAssocsOutOfBunch; // key: track global index, value: vector of global index for events associated out-of-bunch (events that have no in-bunch pileup)

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    fCurrentRun = 0;
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // add extra cuts from JSON
    TString addTrackCutsStr = fConfigCutsJSON.value;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (auto& t : addTrackCuts) {
        fTrackCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    TString configSigNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> sigNamesArray(configSigNamesStr.Tokenize(","));
    // Setting the MC signals
    for (int isig = 0; isig < sigNamesArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(sigNamesArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        fMCSignals.push_back(sig);
      }
    }
    // Add the MCSignals from the JSON config
    TString addMCSignalsStr = fConfigMCSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        fMCSignals.push_back(mcIt);
      }
    }

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // Configure histogram classes for each track cut;
      // Add histogram classes for each track cut and for each requested MC signal (reconstructed tracks with MC truth)
      TString histClasses = "AssocsBarrel_BeforeCuts;";
      for (auto& cut : fTrackCuts) {
        TString nameStr = Form("AssocsBarrel_%s", cut->GetName());
        fHistNamesReco.push_back(nameStr);
        histClasses += Form("%s;", nameStr.Data());
        for (auto& sig : fMCSignals) {
          TString nameStr2 = Form("AssocsCorrectBarrel_%s_%s", cut->GetName(), sig->GetName());
          fHistNamesMCMatched.push_back(nameStr2);
          histClasses += Form("%s;", nameStr2.Data());
          nameStr2 = Form("AssocsIncorrectBarrel_%s_%s", cut->GetName(), sig->GetName());
          fHistNamesMCMatched.push_back(nameStr2);
          histClasses += Form("%s;", nameStr2.Data());
        }
      }

      DefineHistograms(fHistMan, histClasses.Data(), fConfigAddTrackHistogram.value.data());
      if (fConfigPublishAmbiguity) {
        DefineHistograms(fHistMan, "TrackBarrel_AmbiguityInBunch;TrackBarrel_AmbiguityOutOfBunch;", "ambiguity");
      }
      dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                       // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    if (fConfigComputeTPCpostCalib) {
      fCCDB->setURL(fConfigCcdbUrl.value);
      fCCDB->setCaching(true);
      fCCDB->setLocalObjectValidityChecking();
      fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTracks, typename TEventsMC, typename TTracksMC>
  void runTrackSelection(ReducedTracksAssoc const& assocs, TEvents const& events, TTracks const& tracks, TEventsMC const& /*eventsMC*/, TTracksMC const& tracksMC)
  {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

    // TODO: Check if postcalibration needed for MC
    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      if (fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, events.begin().timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
      }

      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
      }

      fCurrentRun = events.begin().runNumber();
    }

    trackSel.reserve(assocs.size());
    trackAmbiguities.reserve(tracks.size());

    // Loop over associations
    for (auto& assoc : assocs) {
      auto event = assoc.template reducedevent_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        trackSel(0);
        continue;
      }

      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);
      if (event.has_reducedMCevent()) {
        VarManager::FillEvent<VarManager::ObjTypes::ReducedEventMC>(event.reducedMCevent());
      }

      auto track = assoc.template reducedtrack_as<TTracks>();
      VarManager::FillTrack<TTrackFillMap>(track);
      // compute quantities which depend on the associated collision, such as DCA
      VarManager::FillTrackCollision<TTrackFillMap>(track, event);

      bool isCorrectAssoc = false;
      if (track.has_reducedMCTrack()) {
        auto trackMC = track.reducedMCTrack();
        auto eventMCfromTrack = trackMC.reducedMCevent();
        if (event.has_reducedMCevent()) {
          isCorrectAssoc = (eventMCfromTrack.globalIndex() == event.reducedMCevent().globalIndex());
        }
        VarManager::FillTrackMC(tracksMC, trackMC);
      }

      if (fConfigQA) {
        fHistMan->FillHistClass("AssocsBarrel_BeforeCuts", VarManager::fgValues);
      }

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

      // compute MC matching decisions and fill histograms for matched associations
      int isig = 0;
      if (filterMap > 0 && track.has_reducedMCTrack() && fConfigQA) {
        // loop over all MC signals
        for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
          // check if this MC signal is matched
          if ((*sig)->CheckSignal(true, track.reducedMCTrack())) {
            // mcDecision |= (static_cast<uint32_t>(1) << isig);
            //  loop over cuts and fill histograms for the cuts that are fulfilled
            for (unsigned int icut = 0; icut < fTrackCuts.size(); icut++) {
              if (filterMap & (static_cast<uint32_t>(1) << icut)) {
                if (isCorrectAssoc) {
                  fHistMan->FillHistClass(fHistNamesMCMatched[icut * 2 * fMCSignals.size() + 2 * isig].Data(), VarManager::fgValues);
                } else {
                  fHistMan->FillHistClass(fHistNamesMCMatched[icut * 2 * fMCSignals.size() + 2 * isig + 1].Data(), VarManager::fgValues);
                }
              }
            } // end loop over cuts
          }
        }

        // fill histograms
        /*for (unsigned int i = 0; i < fMCSignals.size(); i++) {
          if (!(mcDecision & (static_cast<uint32_t>(1) << i))) {
            continue;
          }
          for (unsigned int j = 0; j < fTrackCuts.size(); j++) {
            if (filterMap & (static_cast<uint32_t>(1) << j)) {
              if (isCorrectAssoc) {
                fHistMan->FillHistClass(fHistNamesMCMatched[j * fMCSignals.size() + 2 * i].Data(), VarManager::fgValues);
              } else {
                fHistMan->FillHistClass(fHistNamesMCMatched[j * fMCSignals.size() + 2 * i + 1].Data(), VarManager::fgValues);
              }
            }
          } // end loop over cuts
        } // end loop over MC signals
        */
      } // end if (filterMap > 0)

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

    // QA the collision-track associations
    // TODO: some tracks can be associated to both collisions that have in bunch pileup and collisions from different bunches
    //       So one could QA these tracks separately
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
  } // end runTrackSelection()

  void processSkimmed(ReducedTracksAssoc const& assocs, MyEventsSelected const& events, MyBarrelTracks const& tracks, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runTrackSelection<gkEventFillMap, gkTrackFillMap>(assocs, events, tracks, eventsMC, tracksMC);
  }
  void processSkimmedWithCov(ReducedTracksAssoc const& assocs, MyEventsVtxCovSelected const& events, MyBarrelTracksWithCov const& tracks, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runTrackSelection<gkEventFillMapWithCov, gkTrackFillMapWithCov>(assocs, events, tracks, eventsMC, tracksMC);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run barrel track selection on DQ skimmed track associations", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmedWithCov, "Run barrel track selection on DQ skimmed tracks w/ cov matrix associations", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processDummy, "Dummy function", true);
};

// Produces a table with muon decisions (joinable to the ReducedMuonsAssociations)
// Here one should add all the track cuts needed through the workflow (e.g. cuts for same-event pairing, track for dilepton-track correlations)
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
  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};
  Configurable<std::string> fConfigGeoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<std::string> fConfigMCSignals{"cfgMuonMCSignals", "", "Comma separated list of MC signals"};
  Configurable<std::string> fConfigMCSignalsJSON{"cfgMuonMCsignalsJSON", "", "Additional list of MC signals via JSON"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut*> fMuonCuts;
  std::vector<TString> fHistNamesReco;
  std::vector<TString> fHistNamesMCMatched;
  std::vector<MCSignal*> fMCSignals; // list of signals to be checked

  int fCurrentRun; // current run kept to detect run changes and trigger loading params from CCDB

  std::map<int64_t, std::vector<int64_t>> fNAssocsInBunch;    // key: track global index, value: vector of global index for events associated in-bunch (events that have in-bunch pileup or splitting)
  std::map<int64_t, std::vector<int64_t>> fNAssocsOutOfBunch; // key: track global index, value: vector of global index for events associated out-of-bunch (events that have no in-bunch pileup)

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    fCurrentRun = 0;
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fMuonCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // Add cuts configured via JSON
    TString addCutsStr = fConfigCutsJSON.value;
    if (addCutsStr != "") {
      std::vector<AnalysisCut*> addCuts = dqcuts::GetCutsFromJSON(addCutsStr.Data());
      for (auto& t : addCuts) {
        fMuonCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    TString configSigNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> sigNamesArray(configSigNamesStr.Tokenize(","));
    // Setting the MC signals
    for (int isig = 0; isig < sigNamesArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(sigNamesArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        fMCSignals.push_back(sig);
      }
    }
    // Add the MCSignals from the JSON config
    TString addMCSignalsStr = fConfigMCSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != 1) { // NOTE: only 1 prong signals
          continue;
        }
        fMCSignals.push_back(mcIt);
      }
    }

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // Configure histogram classes for each track cut;
      // Add histogram classes for each track cut and for each requested MC signal (reconstructed tracks with MC truth)
      TString histClasses = "AssocsMuon_BeforeCuts;";
      for (auto& cut : fMuonCuts) {
        TString nameStr = Form("AssocsMuon_%s", cut->GetName());
        fHistNamesReco.push_back(nameStr);
        histClasses += Form("%s;", nameStr.Data());
        for (auto& sig : fMCSignals) {
          TString nameStr2 = Form("AssocsCorrectMuon_%s_%s", cut->GetName(), sig->GetName());
          fHistNamesMCMatched.push_back(nameStr2);
          histClasses += Form("%s;", nameStr2.Data());
          nameStr2 = Form("AssocsIncorrectMuon_%s_%s", cut->GetName(), sig->GetName());
          fHistNamesMCMatched.push_back(nameStr2);
          histClasses += Form("%s;", nameStr2.Data());
        }
      }
      if (fConfigPublishAmbiguity) {
        histClasses += "Muon_AmbiguityInBunch;Muon_AmbiguityOutOfBunch;";
      }

      DefineHistograms(fHistMan, histClasses.Data(), fConfigAddMuonHistogram.value.data());  // define all histograms
      dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                       // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      fCCDB->get<TGeoManager>(fConfigGeoPath);
    }
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvents, typename TMuons>
  void runMuonSelection(ReducedMuonsAssoc const& assocs, TEvents const& events, TMuons const& muons, ReducedMCEvents const& /*eventsMC*/, ReducedMCTracks const& muonsMC)
  {
    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        // LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
        //  If the magnetic field is not found it is configured by had by the user
        VarManager::SetMagneticField(fConfigMagField.value);
      }
      fCurrentRun = events.begin().runNumber();
    }

    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();
    muonSel.reserve(assocs.size());

    for (auto& assoc : assocs) {
      auto event = assoc.template reducedevent_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        muonSel(0);
        continue;
      }
      VarManager::ResetValues(0, VarManager::kNVars);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);
      VarManager::FillEvent<VarManager::ObjTypes::ReducedEventMC>(event.reducedMCevent());

      auto track = assoc.template reducedmuon_as<TMuons>();
      VarManager::FillTrack<TMuonFillMap>(track);

      bool isCorrectAssoc = false;
      if (track.has_reducedMCTrack()) {
        auto trackMC = track.reducedMCTrack();
        auto eventMCfromTrack = trackMC.reducedMCevent();
        isCorrectAssoc = (eventMCfromTrack.globalIndex() == event.reducedMCevent().globalIndex());
        VarManager::FillTrackMC(muonsMC, trackMC);
      }

      if (fConfigQA) {
        fHistMan->FillHistClass("AssocsMuon_BeforeCuts", VarManager::fgValues);
      }

      int iCut = 0;
      uint32_t filterMap = static_cast<uint32_t>(0);
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, iCut++) {
        if ((*cut)->IsSelected(VarManager::fgValues)) {
          filterMap |= (static_cast<uint32_t>(1) << iCut);
          if (fConfigQA) {
            fHistMan->FillHistClass(fHistNamesReco[iCut].Data(), VarManager::fgValues);
          }
        }
      } // end loop over cuts
      muonSel(filterMap);
      muonAmbiguities.reserve(muons.size());

      // if no filter fulfilled, continue
      if (!filterMap) {
        continue;
      }

      // everything below is related to filling QA histograms
      if (!fConfigQA) {
        continue;
      }

      // compute MC matching decisions
      uint32_t mcDecision = static_cast<uint32_t>(0);
      int isig = 0;
      for (auto sig = fMCSignals.begin(); sig != fMCSignals.end(); sig++, isig++) {
        if constexpr ((TMuonFillMap & VarManager::ObjTypes::ReducedMuon) > 0) {
          if (track.has_reducedMCTrack()) {
            if ((*sig)->CheckSignal(true, track.reducedMCTrack())) {
              mcDecision |= (static_cast<uint32_t>(1) << isig);
            }
          }
        }
      }

      // fill histograms
      for (unsigned int i = 0; i < fMCSignals.size(); i++) {
        if (!(mcDecision & (static_cast<uint32_t>(1) << i))) {
          continue;
        }
        for (unsigned int j = 0; j < fMuonCuts.size(); j++) {
          if (filterMap & (static_cast<uint32_t>(1) << j)) {
            if (isCorrectAssoc) {
              fHistMan->FillHistClass(fHistNamesMCMatched[j * 2 * fMCSignals.size() + 2 * i].Data(), VarManager::fgValues);
            } else {
              fHistMan->FillHistClass(fHistNamesMCMatched[j * 2 * fMCSignals.size() + 2 * i + 1].Data(), VarManager::fgValues);
            }
          }
        } // end loop over cuts
      } // end loop over MC signals

      // count the number of associations per track
      if (fConfigPublishAmbiguity && filterMap > 0) {
        if (event.isEventSelected_bit(1)) {
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
    } // end loop over assocs

    // QA the collision-track associations
    // TODO: some tracks can be associated to both collisions that have in bunch pileup and collisions from different bunches
    //       So one could QA these tracks separately
    if (fConfigPublishAmbiguity) {
      for (auto& [trackIdx, evIndices] : fNAssocsInBunch) {
        if (evIndices.size() == 1) {
          continue;
        }
        auto track = muons.rawIteratorAt(trackIdx);
        VarManager::ResetValues(0, VarManager::kNVars);
        VarManager::FillTrack<TMuonFillMap>(track);
        VarManager::fgValues[VarManager::kMuonNAssocsInBunch] = static_cast<float>(evIndices.size());
        fHistMan->FillHistClass("Muon_AmbiguityInBunch", VarManager::fgValues);
      } // end loop over in-bunch ambiguous tracks

      for (auto& [trackIdx, evIndices] : fNAssocsOutOfBunch) {
        if (evIndices.size() == 1) {
          continue;
        }
        auto track = muons.rawIteratorAt(trackIdx);
        VarManager::ResetValues(0, VarManager::kNVars);
        VarManager::FillTrack<TMuonFillMap>(track);
        VarManager::fgValues[VarManager::kMuonNAssocsOutOfBunch] = static_cast<float>(evIndices.size());
        fHistMan->FillHistClass("Muon_AmbiguityOutOfBunch", VarManager::fgValues);
      } // end loop over out-of-bunch ambiguous tracks

      // publish the ambiguity table
      for (auto& track : muons) {
        int8_t nInBunch = 0;
        if (fNAssocsInBunch.find(track.globalIndex()) != fNAssocsInBunch.end()) {
          nInBunch = fNAssocsInBunch[track.globalIndex()].size();
        }
        int8_t nOutOfBunch = 0;
        if (fNAssocsOutOfBunch.find(track.globalIndex()) != fNAssocsOutOfBunch.end()) {
          nOutOfBunch = fNAssocsOutOfBunch[track.globalIndex()].size();
        }
        muonAmbiguities(nInBunch, nOutOfBunch);
      }
    }
  }

  void processSkimmed(ReducedMuonsAssoc const& assocs, MyEventsSelected const& events, MyMuonTracks const& muons, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runMuonSelection<gkEventFillMap, gkMuonFillMap>(assocs, events, muons, eventsMC, tracksMC);
  }
  void processSkimmedWithCov(ReducedMuonsAssoc const& assocs, MyEventsVtxCovSelected const& events, MyMuonTracksWithCov const& muons, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    runMuonSelection<gkEventFillMapWithCov, gkMuonFillMapWithCov>(assocs, events, muons, eventsMC, tracksMC);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisMuonSelection, processSkimmed, "Run muon selection on DQ skimmed muons", false);
  PROCESS_SWITCH(AnalysisMuonSelection, processSkimmedWithCov, "Run muon selection on DQ skimmed muons, with event and track covariances", false);
  PROCESS_SWITCH(AnalysisMuonSelection, processDummy, "Dummy function", true);
};

// Run the prefilter selection (e.g. electron prefiltering for photon conversions)
// This takes uses a sample of tracks selected with loose cuts (fConfigPrefilterTrackCut) and combines them
//  with the sample of tracks to be used in downstream analysis (fConfigTrackCuts). If a pair is found to pass
//  the pair prefilter cut (cfgPrefilterPairCut), the analysis track is tagged to be removed from analysis.
// TODO: Add optional QA histograms
struct AnalysisPrefilterSelection {
  Produces<aod::Prefilter> prefilter; // joinable with ReducedTracksAssoc

  // Configurables
  Configurable<std::string> fConfigPrefilterTrackCut{"cfgPrefilterTrackCut", "", "Prefilter track cut"};
  Configurable<std::string> fConfigPrefilterPairCut{"cfgPrefilterPairCut", "", "Prefilter pair cut"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Track cuts for which to run the prefilter"};
  // Track related options
  Configurable<bool> fPropTrack{"cfgPropTrack", false, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};

  std::map<uint32_t, uint32_t> fPrefilterMap;
  AnalysisCompositeCut* fPairCut;
  uint32_t fPrefilterMask;
  int fPrefilterCutBit;

  Preslice<aod::ReducedTracksAssoc> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

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

  template <uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runPrefilter(TEvent const& event, soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& assocs, TTracks const& /*tracks*/)
  {
    if (fPrefilterCutBit < 0 || fPrefilterMask == 0) {
      return;
    }

    for (auto& [assoc1, assoc2] : o2::soa::combinations(assocs, assocs)) {
      auto track1 = assoc1.template reducedtrack_as<TTracks>();
      auto track2 = assoc2.template reducedtrack_as<TTracks>();

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
      VarManager::FillPair<VarManager::kDecayToEE, TTrackFillMap>(track1, track2);
      if (fPropTrack) {
        VarManager::FillPairCollision<VarManager::kDecayToEE, TTrackFillMap>(event, track1, track2);
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
  }

  void processBarrelSkimmed(MyEvents const& events, soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& assocs, MyBarrelTracks const& tracks)
  {

    fPrefilterMap.clear();

    for (auto& event : events) {
      auto groupedAssocs = assocs.sliceBy(trackAssocsPerCollision, event.globalIndex());
      if (groupedAssocs.size() > 1) {
        runPrefilter<gkTrackFillMap>(event, groupedAssocs, tracks);
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
        auto track = assoc.template reducedtrack_as<MyBarrelTracks>();
        mymap = -1;
        if (fPrefilterMap.find(track.globalIndex()) != fPrefilterMap.end()) {
          // NOTE: publish the bitwise negated bits (~), so there will be zeroes for cuts that failed the prefiltering and 1 everywhere else
          mymap = ~fPrefilterMap[track.globalIndex()];
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

  PROCESS_SWITCH(AnalysisPrefilterSelection, processBarrelSkimmed, "Run Prefilter selection on reduced tracks", false);
  PROCESS_SWITCH(AnalysisPrefilterSelection, processDummy, "Do nothing", true);
};

// Run the same-event pairing
// This task assumes that both legs of the resonance fulfill the same cuts (symmetric decay channel)
//  Runs combinatorics for barrel-barrel, muon-muon and barrel-muon combinations
// TODO: implement properly the barrel-muon combinations
// The task implements also process functions for running event mixing
struct AnalysisSameEventPairing {

  Produces<aod::Dielectrons> dielectronList;
  Produces<aod::Dimuons> dimuonList;
  Produces<aod::DielectronsExtra> dielectronsExtraList;
  Produces<aod::DielectronsInfo> dielectronInfoList;
  Produces<aod::DielectronsAll> dielectronAllList;
  Produces<aod::DimuonsExtra> dimuonsExtraList;
  Produces<aod::DimuonsAll> dimuonAllList;
  Produces<aod::DileptonsMiniTreeGen> dileptonMiniTreeGen;
  Produces<aod::DileptonsMiniTreeRec> dileptonMiniTreeRec;
  Produces<aod::DileptonsInfo> dileptonInfoList;
  Produces<aod::JPsieeCandidates> PromptNonPromptSepTable;
  Produces<aod::OniaMCTruth> MCTruthTableEffi;
  Produces<aod::DileptonPolarization> dileptonPolarList;

  o2::base::MatLayerCylSet* fLUT = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  OutputObj<THashList> fOutputList{"output"};

  struct : ConfigurableGroup {
    Configurable<std::string> track{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
    Configurable<std::string> muon{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
    Configurable<std::string> pair{"cfgPairCuts", "", "Comma separated list of pair cuts, !!! Use only if you know what you are doing, otherwise leave empty"};
    Configurable<std::string> MCgenAcc{"cfgMCGenAccCuts", "", "Comma separated list of MC generated particles acceptance cuts, !!! Use only if you know what you are doing, otherwise leave empty"};
    // TODO: Add pair cuts via JSON
  } fConfigCuts;

  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  struct : ConfigurableGroup {
    Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } fConfigCCDB;

  struct : ConfigurableGroup {
    Configurable<bool> useRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
    Configurable<float> magField{"cfgMagField", 5.0f, "Manually set magnetic field"};
    Configurable<bool> flatTables{"cfgFlatTables", false, "Produce a single flat tables with all relevant information of the pairs and single tracks"};
    Configurable<bool> polarTables{"cfgPolarTables", false, "Produce tables with dilepton polarization information"};
    Configurable<bool> useKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
    Configurable<bool> useAbsDCA{"cfgUseAbsDCA", false, "Use absolute DCA minimization instead of chi^2 minimization in secondary vertexing"};
    Configurable<bool> propToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
    Configurable<bool> corrFullGeo{"cfgCorrFullGeo", false, "Use full geometry to correct for MCS effects in track propagation"};
    Configurable<bool> noCorr{"cfgNoCorrFwdProp", false, "Do not correct for MCS effects in track propagation"};
    Configurable<std::string> collisionSystem{"syst", "pp", "Collision system, pp or PbPb"};
    Configurable<float> centerMassEnergy{"energy", 13600, "Center of mass energy in GeV"};
  } fConfigOptions;

  struct : ConfigurableGroup {
    Configurable<std::string> genSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};
    Configurable<std::string> genSignalsJSON{"cfgMCGenSignalsJSON", "", "Additional list of MC signals (generated) via JSON"};
    Configurable<std::string> recSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
    Configurable<std::string> recSignalsJSON{"cfgMCRecSignalsJSON", "", "Comma separated list of MC signals (reconstructed) via JSON"};
    Configurable<bool> skimSignalOnly{"cfgSkimSignalOnly", false, "Configurable to select only matched candidates"};
  } fConfigMC;

  struct : ConfigurableGroup {
    Configurable<bool> fConfigMiniTree{"useMiniTree.cfgMiniTree", false, "Produce a single flat table with minimal information for analysis"};
    Configurable<float> fConfigMiniTreeMinMass{"useMiniTree.cfgMiniTreeMinMass", 2, "Min. mass cut for minitree"};
    Configurable<float> fConfigMiniTreeMaxMass{"useMiniTree.cfgMiniTreeMaxMass", 5, "Max. mass cut for minitree"};
  } useMiniTree;

  // Track related options
  Configurable<bool> fPropTrack{"cfgPropTrack", true, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  // Filter filterEventSelected = aod::dqanalysisflags::isEventSelected & uint32_t(1);
  Filter eventFilter = aod::dqanalysisflags::isEventSelected > static_cast<uint32_t>(0);

  HistogramManager* fHistMan;

  // keep histogram class names in maps, so we don't have to buld their names in the pair loops
  std::map<int, std::vector<TString>> fTrackHistNames;
  std::map<int, std::vector<TString>> fBarrelHistNamesMCmatched;
  std::map<int, std::vector<TString>> fMuonHistNames;
  std::map<int, std::vector<TString>> fMuonHistNamesMCmatched;
  std::vector<MCSignal*> fRecMCSignals;
  std::vector<MCSignal*> fGenMCSignals;

  std::vector<AnalysisCompositeCut> fPairCuts;
  std::vector<AnalysisCut*> fMCGenAccCuts;
  bool fUseMCGenAccCut = false;

  uint32_t fTrackFilterMask; // mask for the track cuts required in this task to be applied on the barrel cuts produced upstream
  uint32_t fMuonFilterMask;  // mask for the muon cuts required in this task to be applied on the muon cuts produced upstream
  int fNCutsBarrel;
  int fNCutsMuon;
  int fNPairCuts;
  bool fHasTwoProngGenMCsignals = false;

  bool fEnableBarrelHistos;
  bool fEnableMuonHistos;

  Preslice<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter>> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    bool isMCGen = context.mOptions.get<bool>("processMCGen") || context.mOptions.get<bool>("processMCGenWithGrouping") || context.mOptions.get<bool>("processBarrelOnlySkimmed");
    VarManager::SetDefaultVarNames();

    fEnableBarrelHistos = context.mOptions.get<bool>("processAllSkimmed") || context.mOptions.get<bool>("processBarrelOnlySkimmed") || context.mOptions.get<bool>("processBarrelOnlyWithCollSkimmed");
    fEnableMuonHistos = context.mOptions.get<bool>("processAllSkimmed") || context.mOptions.get<bool>("processMuonOnlySkimmed");

    // Keep track of all the histogram class names to avoid composing strings in the pairing loop
    TString histNames = "";
    TString cutNamesStr = fConfigCuts.pair.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // get the list of cuts for tracks/muons, check that they were played by the barrel/muon selection tasks
    //   and make a mask for active cuts (barrel and muon selection tasks may run more cuts, needed for other analyses)
    TString trackCutsStr = fConfigCuts.track.value;
    TObjArray* objArrayTrackCuts = nullptr;
    if (!trackCutsStr.IsNull()) {
      objArrayTrackCuts = trackCutsStr.Tokenize(",");
    }
    TString muonCutsStr = fConfigCuts.muon.value;
    TObjArray* objArrayMuonCuts = nullptr;
    if (!muonCutsStr.IsNull()) {
      objArrayMuonCuts = muonCutsStr.Tokenize(",");
    }

    // Setting the MC rec signal names
    TString sigNamesStr = fConfigMC.recSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required
          continue;
        }
        fRecMCSignals.push_back(sig);
      }
    }

    // Add the MCSignals from the JSON config
    TString addMCSignalsStr = fConfigMC.recSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != 2) { // NOTE: only 2 prong signals
          continue;
        }
        fRecMCSignals.push_back(mcIt);
      }
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

          if (fEnableBarrelHistos) {
            // assign the pair hist directories for the current cut
            std::vector<TString> names = {
              Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
            if (fConfigQA) {
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
            TString cutNamesStr = fConfigCuts.pair.value;
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

            // assign hist directories for the MC matched pairs for each (track cut,MCsignal) combination
            if (!sigNamesStr.IsNull()) {
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
                auto sig = fRecMCSignals.at(isig);
                names = {
                  Form("PairsBarrelSEPM_%s_%s", objArray->At(icut)->GetName(), sig->GetName()),
                  Form("PairsBarrelSEPP_%s_%s", objArray->At(icut)->GetName(), sig->GetName()),
                  Form("PairsBarrelSEMM_%s_%s", objArray->At(icut)->GetName(), sig->GetName())};
                if (fConfigQA) {
                  names.push_back(Form("PairsBarrelSEPMCorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPMIncorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousInBunch_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousInBunchCorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousInBunchIncorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousOutOfBunch_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousOutOfBunchCorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                  names.push_back(Form("PairsBarrelSEPM_ambiguousOutOfBunchIncorrectAssoc_%s_%s", objArray->At(icut)->GetName(), sig->GetName()));
                }
                for (auto& n : names) {
                  histNames += Form("%s;", n.Data());
                }
                fBarrelHistNamesMCmatched.try_emplace(icut * fRecMCSignals.size() + isig, names);
              } // end loop over MC signals
            }
          } // end if enableBarrelHistos
        }
      }
    }

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
            if (fConfigQA) {
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
            TString cutNamesStr = fConfigCuts.pair.value;
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
                if (fConfigQA) {
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

    // get the mc generated acceptance cuts
    TString mcgenCutsStr = fConfigCuts.MCgenAcc.value;
    if (!mcgenCutsStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(mcgenCutsStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        AnalysisCut* cut = dqcuts::GetAnalysisCut(objArray->At(icut)->GetName());
        if (cut != nullptr) {
          fMCGenAccCuts.push_back(cut);
        }
      }
      fUseMCGenAccCut = true;
    }

    // Add histogram classes for each specified MCsignal at the generator level
    // TODO: create a std::vector of hist classes to be used at Fill time, to avoid using Form in the process function
    TString sigGenNamesStr = fConfigMC.genSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        fGenMCSignals.push_back(sig);
      }
    }

    // Add the MCSignals from the JSON config
    TString addMCSignalsGenStr = fConfigMC.genSignalsJSON.value;
    if (addMCSignalsGenStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsGenStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() > 2) { // NOTE: only 2 prong signals
          continue;
        }
        fGenMCSignals.push_back(mcIt);
      }
    }

    if (isMCGen) {
      for (auto& sig : fGenMCSignals) {
        if (sig->GetNProngs() == 1) {
          histNames += Form("MCTruthGen_%s;", sig->GetName()); // TODO: Add these names to a std::vector to avoid using Form in the process function
          histNames += Form("MCTruthGenSel_%s;", sig->GetName());
        } else if (sig->GetNProngs() == 2) {
          histNames += Form("MCTruthGenPair_%s;", sig->GetName());
          histNames += Form("MCTruthGenPairSel_%s;", sig->GetName());
          fHasTwoProngGenMCsignals = true;
          // for these pair level signals, also add histograms for each MCgenAcc cut if specified
          if (fUseMCGenAccCut) {
            for (auto& cut : fMCGenAccCuts) {
              histNames += Form("MCTruthGenPairSel_%s_%s;", cut->GetName(), sig->GetName());
            }
          }
        }
      }
    }

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

    DefineHistograms(fHistMan, histNames.Data(), fConfigAddSEPHistogram.value.data());     // define all histograms
    dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                       // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void initParamsFromCCDB(uint64_t timestamp, bool withTwoProngFitter = true)
  {
    if (fConfigOptions.useRemoteField.value) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigCCDB.grpMagPath, timestamp);
      float magField = 0.0;
      if (grpmag != nullptr) {
        magField = grpmag->getNominalL3Field();
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

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <bool TTwoProngFitter, int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTrackAssocs, typename TTracks>
  void runSameEventPairing(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& assocs, TTracks const& /*tracks*/, ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& /*mcTracks*/)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    if (fCurrentRun != events.begin().runNumber()) {
      initParamsFromCCDB(events.begin().timestamp(), TTwoProngFitter);
      fCurrentRun = events.begin().runNumber();
    }

    TString cutNames = fConfigCuts.track.value;
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    std::map<int, std::vector<TString>> histNamesMC = fBarrelHistNamesMCmatched;
    int ncuts = fNCutsBarrel;
    if constexpr (TPairType == VarManager::kDecayToMuMu) {
      cutNames = fConfigCuts.muon.value;
      histNames = fMuonHistNames;
      histNamesMC = fMuonHistNamesMCmatched;
      ncuts = fNCutsMuon;
    }

    uint32_t twoTrackFilter = static_cast<uint32_t>(0);
    int sign1 = 0;
    int sign2 = 0;
    uint32_t mcDecision = static_cast<uint32_t>(0);
    bool isCorrectAssoc_leg1 = false;
    bool isCorrectAssoc_leg2 = false;
    dielectronList.reserve(1);
    dimuonList.reserve(1);
    dielectronsExtraList.reserve(1);
    dimuonsExtraList.reserve(1);
    dielectronInfoList.reserve(1);
    dileptonInfoList.reserve(1);
    if (fConfigOptions.flatTables.value) {
      dielectronAllList.reserve(1);
      dimuonAllList.reserve(1);
    }
    if (useMiniTree.fConfigMiniTree) {
      dileptonMiniTreeGen.reserve(1);
      dileptonMiniTreeRec.reserve(1);
    }
    if (fConfigOptions.polarTables.value) {
      dileptonPolarList.reserve(1);
    }
    constexpr bool eventHasQvector = ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0);
    constexpr bool trackHasCov = ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelCov) > 0);

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      uint8_t evSel = event.isEventSelected_raw();
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
      VarManager::FillEvent<VarManager::ObjTypes::ReducedEventMC>(event.reducedMCevent(), VarManager::fgValues);

      auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());
      if (groupedAssocs.size() == 0) {
        continue;
      }

      for (auto& [a1, a2] : o2::soa::combinations(groupedAssocs, groupedAssocs)) {

        if constexpr (TPairType == VarManager::kDecayToEE) {
          twoTrackFilter = a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isBarrelSelectedPrefilter_raw() & fTrackFilterMask;

          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }

          auto t1 = a1.template reducedtrack_as<TTracks>();
          auto t2 = a2.template reducedtrack_as<TTracks>();
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

          // run MC matching for this pair
          int isig = 0;
          mcDecision = 0;
          for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
            if (t1.has_reducedMCTrack() && t2.has_reducedMCTrack()) {
              if ((*sig)->CheckSignal(true, t1.reducedMCTrack(), t2.reducedMCTrack())) {
                mcDecision |= (static_cast<uint32_t>(1) << isig);
              }
            }
          } // end loop over MC signals
          if (t1.has_reducedMCTrack() && t2.has_reducedMCTrack()) {
            isCorrectAssoc_leg1 = (t1.reducedMCTrack().reducedMCevent() == event.reducedMCevent());
            isCorrectAssoc_leg2 = (t2.reducedMCTrack().reducedMCevent() == event.reducedMCevent());
          }

          VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
          if (fPropTrack) {
            VarManager::FillPairCollision<TPairType, TTrackFillMap>(event, t1, t2);
          }
          if constexpr (TTwoProngFitter) {
            VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigOptions.propToPCA);
          }
          if constexpr (eventHasQvector) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }
          if (!fConfigMC.skimSignalOnly || (fConfigMC.skimSignalOnly && mcDecision > 0)) {
            dielectronList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                           VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                           t1.sign() + t2.sign(), twoTrackFilter, mcDecision);

            if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackCollInfo) > 0) {
              dielectronInfoList(t1.collisionId(), t1.trackId(), t2.trackId());
              dileptonInfoList(t1.collisionId(), event.posX(), event.posY(), event.posZ());
            }
            if constexpr (trackHasCov && TTwoProngFitter) {
              dielectronsExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
              if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackCollInfo) > 0) {
                if (fConfigOptions.flatTables.value && t1.has_reducedMCTrack() && t2.has_reducedMCTrack()) {
                  dielectronAllList(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), twoTrackFilter, mcDecision,
                                    t1.pt(), t1.eta(), t1.phi(), t1.itsClusterMap(), t1.itsChi2NCl(), t1.tpcNClsCrossedRows(), t1.tpcNClsFound(), t1.tpcChi2NCl(), t1.dcaXY(), t1.dcaZ(), t1.tpcSignal(), t1.tpcNSigmaEl(), t1.tpcNSigmaPi(), t1.tpcNSigmaPr(), t1.beta(), t1.tofNSigmaEl(), t1.tofNSigmaPi(), t1.tofNSigmaPr(),
                                    t2.pt(), t2.eta(), t2.phi(), t2.itsClusterMap(), t2.itsChi2NCl(), t2.tpcNClsCrossedRows(), t2.tpcNClsFound(), t2.tpcChi2NCl(), t2.dcaXY(), t2.dcaZ(), t2.tpcSignal(), t2.tpcNSigmaEl(), t2.tpcNSigmaPi(), t2.tpcNSigmaPr(), t2.beta(), t2.tofNSigmaEl(), t2.tofNSigmaPi(), t2.tofNSigmaPr(),
                                    VarManager::fgValues[VarManager::kKFTrack0DCAxyz], VarManager::fgValues[VarManager::kKFTrack1DCAxyz], VarManager::fgValues[VarManager::kKFDCAxyzBetweenProngs], VarManager::fgValues[VarManager::kKFTrack0DCAxy], VarManager::fgValues[VarManager::kKFTrack1DCAxy], VarManager::fgValues[VarManager::kKFDCAxyBetweenProngs],
                                    VarManager::fgValues[VarManager::kKFTrack0DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack0DeviationxyFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationxyFromPV],
                                    VarManager::fgValues[VarManager::kKFMass], VarManager::fgValues[VarManager::kKFChi2OverNDFGeo], VarManager::fgValues[VarManager::kVertexingLxyz], VarManager::fgValues[VarManager::kVertexingLxyzOverErr], VarManager::fgValues[VarManager::kVertexingLxy], VarManager::fgValues[VarManager::kVertexingLxyOverErr], VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr], VarManager::fgValues[VarManager::kKFCosPA], VarManager::fgValues[VarManager::kKFJpsiDCAxyz], VarManager::fgValues[VarManager::kKFJpsiDCAxy],
                                    VarManager::fgValues[VarManager::kKFPairDeviationFromPV], VarManager::fgValues[VarManager::kKFPairDeviationxyFromPV],
                                    VarManager::fgValues[VarManager::kKFMassGeoTop], VarManager::fgValues[VarManager::kKFChi2OverNDFGeoTop],
                                    VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingTauxyProjected],
                                    VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
                }
                if (fConfigOptions.polarTables.value && t1.has_reducedMCTrack() && t2.has_reducedMCTrack()) {
                  dileptonPolarList(VarManager::fgValues[VarManager::kCosThetaHE], VarManager::fgValues[VarManager::kPhiHE], VarManager::fgValues[VarManager::kPhiTildeHE],
                                    VarManager::fgValues[VarManager::kCosThetaCS], VarManager::fgValues[VarManager::kPhiCS], VarManager::fgValues[VarManager::kPhiTildeCS],
                                    VarManager::fgValues[VarManager::kCosThetaPP], VarManager::fgValues[VarManager::kPhiPP], VarManager::fgValues[VarManager::kPhiTildePP],
                                    VarManager::fgValues[VarManager::kCosThetaRM],
                                    VarManager::fgValues[VarManager::kCosThetaStarTPC], VarManager::fgValues[VarManager::kCosThetaStarFT0A], VarManager::fgValues[VarManager::kCosThetaStarFT0C]);
                }
              }
            }
          }
        }

        if constexpr (TPairType == VarManager::kDecayToMuMu) {
          twoTrackFilter = a1.isMuonSelected_raw() & a2.isMuonSelected_raw() & fMuonFilterMask;
          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }
          auto t1 = a1.template reducedmuon_as<TTracks>();
          auto t2 = a2.template reducedmuon_as<TTracks>();
          if (t1.matchMCHTrackId() == t2.matchMCHTrackId() && t1.matchMCHTrackId() >= 0)
            continue;
          if (t1.matchMFTTrackId() == t2.matchMFTTrackId() && t1.matchMFTTrackId() >= 0)
            continue;
          sign1 = t1.sign();
          sign2 = t2.sign();
          // store the ambiguity number of the two dilepton legs in the last 4 digits of the two-track filter
          if (t1.muonAmbiguityInBunch() > 1) {
            twoTrackFilter |= (static_cast<uint32_t>(1) << 28);
          }
          if (t2.muonAmbiguityInBunch() > 1) {
            twoTrackFilter |= (static_cast<uint32_t>(1) << 29);
          }
          if (t1.muonAmbiguityOutOfBunch() > 1) {
            twoTrackFilter |= (static_cast<uint32_t>(1) << 30);
          }
          if (t2.muonAmbiguityOutOfBunch() > 1) {
            twoTrackFilter |= (static_cast<uint32_t>(1) << 31);
          }

          // run MC matching for this pair
          int isig = 0;
          mcDecision = 0;
          for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
            if (t1.has_reducedMCTrack() && t2.has_reducedMCTrack()) {
              if ((*sig)->CheckSignal(true, t1.reducedMCTrack(), t2.reducedMCTrack())) {
                mcDecision |= (static_cast<uint32_t>(1) << isig);
              }
            }
          } // end loop over MC signals

          if (t1.has_reducedMCTrack() && t2.has_reducedMCTrack()) {
            isCorrectAssoc_leg1 = (t1.reducedMCTrack().reducedMCevent() == event.reducedMCevent());
            isCorrectAssoc_leg2 = (t2.reducedMCTrack().reducedMCevent() == event.reducedMCevent());
          }

          VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
          if (fPropTrack) {
            VarManager::FillPairCollision<TPairType, TTrackFillMap>(event, t1, t2);
          }
          if constexpr (TTwoProngFitter) {
            VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigOptions.propToPCA);
          }
          if constexpr (eventHasQvector) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }

          dimuonList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                     VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                     t1.sign() + t2.sign(), twoTrackFilter, mcDecision);
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedMuonCollInfo) > 0) {
            dileptonInfoList(t1.collisionId(), event.posX(), event.posY(), event.posZ());
          }

          if constexpr (TTwoProngFitter) {
            dimuonsExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingLz], VarManager::fgValues[VarManager::kVertexingLxy]);
            if (fConfigOptions.flatTables.value && t1.has_reducedMCTrack() && t2.has_reducedMCTrack()) {
              dimuonAllList(event.posX(), event.posY(), event.posZ(), event.numContrib(),
                            event.selection_raw(), evSel,
                            event.reducedMCevent().mcPosX(), event.reducedMCevent().mcPosY(), event.reducedMCevent().mcPosZ(),
                            VarManager::fgValues[VarManager::kMass],
                            mcDecision,
                            VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), VarManager::fgValues[VarManager::kVertexingChi2PCA],
                            VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingTauzErr],
                            VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr],
                            VarManager::fgValues[VarManager::kCosPointingAngle],
                            VarManager::fgValues[VarManager::kPt1], VarManager::fgValues[VarManager::kEta1], VarManager::fgValues[VarManager::kPhi1], t1.sign(),
                            VarManager::fgValues[VarManager::kPt2], VarManager::fgValues[VarManager::kEta2], VarManager::fgValues[VarManager::kPhi2], t2.sign(),
                            t1.fwdDcaX(), t1.fwdDcaY(), t2.fwdDcaX(), t2.fwdDcaY(),
                            t1.mcMask(), t2.mcMask(),
                            t1.chi2MatchMCHMID(), t2.chi2MatchMCHMID(),
                            t1.chi2MatchMCHMFT(), t2.chi2MatchMCHMFT(),
                            t1.chi2(), t2.chi2(),
                            t1.reducedMCTrack().pt(), t1.reducedMCTrack().eta(), t1.reducedMCTrack().phi(), t1.reducedMCTrack().e(),
                            t2.reducedMCTrack().pt(), t2.reducedMCTrack().eta(), t2.reducedMCTrack().phi(), t2.reducedMCTrack().e(),
                            t1.reducedMCTrack().vx(), t1.reducedMCTrack().vy(), t1.reducedMCTrack().vz(), t1.reducedMCTrack().vt(),
                            t2.reducedMCTrack().vx(), t2.reducedMCTrack().vy(), t2.reducedMCTrack().vz(), t2.reducedMCTrack().vt(),
                            (twoTrackFilter & (static_cast<uint32_t>(1) << 28)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 29)), (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31)),
                            -999.0, -999.0, -999.0, -999.0, -999.0,
                            -999.0, -999.0, -999.0, -999.0, -999.0,
                            -999.0, VarManager::fgValues[VarManager::kMultDimuons],
                            VarManager::fgValues[VarManager::kVertexingPz], VarManager::fgValues[VarManager::kVertexingSV]);
            }
          }
        }
        // TODO: the model for the electron-muon combination has to be thought through
        /*if constexpr (TPairType == VarManager::kElectronMuon) {
          twoTrackFilter = a1.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isMuonSelected_raw() & fTwoTrackFilterMask;
        }*/

        // Fill histograms
        bool isAmbiInBunch = false;
        bool isAmbiOutOfBunch = false;
        bool isCorrect_pair = false;
        if (isCorrectAssoc_leg1 && isCorrectAssoc_leg2)
          isCorrect_pair = true;

        for (int icut = 0; icut < ncuts; icut++) {
          if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
            isAmbiInBunch = (twoTrackFilter & (static_cast<uint32_t>(1) << 28)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 29));
            isAmbiOutOfBunch = (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31));
            if (sign1 * sign2 < 0) {                                                    // +- pairs
              fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues); // reconstructed, unmatched
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {        // loop over MC signals
                if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                  PromptNonPromptSepTable(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kRap], VarManager::fgValues[VarManager::kPhi], VarManager::fgValues[VarManager::kVertexingTauxyProjected], VarManager::fgValues[VarManager::kVertexingTauxyProjectedPoleJPsiMass], VarManager::fgValues[VarManager::kVertexingTauzProjected], isAmbiInBunch, isAmbiOutOfBunch, isCorrect_pair, VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);
                  fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][0].Data(), VarManager::fgValues); // matched signal
                  if (useMiniTree.fConfigMiniTree) {
                    if constexpr (TPairType == VarManager::kDecayToMuMu) {
                      twoTrackFilter = a1.isMuonSelected_raw() & a2.isMuonSelected_raw() & fMuonFilterMask;
                      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
                        continue;
                      }
                      auto t1 = a1.template reducedmuon_as<TTracks>();
                      auto t2 = a2.template reducedmuon_as<TTracks>();

                      float dileptonMass = VarManager::fgValues[VarManager::kMass];
                      if (dileptonMass > useMiniTree.fConfigMiniTreeMinMass && dileptonMass < useMiniTree.fConfigMiniTreeMaxMass) {
                        // In the miniTree the positive daughter is positioned as first
                        if (t1.sign() > 0) {
                          dileptonMiniTreeRec(mcDecision,
                                              VarManager::fgValues[VarManager::kMass],
                                              VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], VarManager::fgValues[VarManager::kCentFT0C],
                                              t1.reducedMCTrack().pt(), t1.reducedMCTrack().eta(), t1.reducedMCTrack().phi(),
                                              t2.reducedMCTrack().pt(), t2.reducedMCTrack().eta(), t2.reducedMCTrack().phi(),
                                              t1.pt(), t1.eta(), t1.phi(),
                                              t2.pt(), t2.eta(), t2.phi());
                        } else {
                          dileptonMiniTreeRec(mcDecision,
                                              VarManager::fgValues[VarManager::kMass],
                                              VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], VarManager::fgValues[VarManager::kCentFT0C],
                                              t2.reducedMCTrack().pt(), t2.reducedMCTrack().eta(), t2.reducedMCTrack().phi(),
                                              t1.reducedMCTrack().pt(), t1.reducedMCTrack().eta(), t1.reducedMCTrack().phi(),
                                              t2.pt(), t2.eta(), t2.phi(),
                                              t1.pt(), t1.eta(), t1.phi());
                        }
                      }
                    }
                  }
                  if (fConfigQA) {
                    if (isCorrectAssoc_leg1 && isCorrectAssoc_leg2) { // correct track-collision association
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][3].Data(), VarManager::fgValues);
                    } else { // incorrect track-collision association
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][4].Data(), VarManager::fgValues);
                    }
                    if (isAmbiInBunch) { // ambiguous in bunch
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][5].Data(), VarManager::fgValues);
                      if (isCorrectAssoc_leg1 && isCorrectAssoc_leg2) {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][6].Data(), VarManager::fgValues);
                      } else {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][7].Data(), VarManager::fgValues);
                      }
                    }
                    if (isAmbiOutOfBunch) { // ambiguous out of bunch
                      fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][8].Data(), VarManager::fgValues);
                      if (isCorrectAssoc_leg1 && isCorrectAssoc_leg2) {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][9].Data(), VarManager::fgValues);
                      } else {
                        fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][10].Data(), VarManager::fgValues);
                      }
                    }
                  }
                }
                if (fConfigQA) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][3].Data(), VarManager::fgValues);
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][3 + 3].Data(), VarManager::fgValues);
                  }
                }
              }
            } else {
              if (sign1 > 0) { // ++ pairs
                fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
                for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                  if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                    fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][1].Data(), VarManager::fgValues);
                  }
                }
                if (fConfigQA) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][4].Data(), VarManager::fgValues);
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][4 + 3].Data(), VarManager::fgValues);
                  }
                }
              } else { // -- pairs
                fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
                for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                  if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                    fHistMan->FillHistClass(histNamesMC[icut * fRecMCSignals.size() + isig][2].Data(), VarManager::fgValues);
                  }
                }
                if (fConfigQA) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][5].Data(), VarManager::fgValues);
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][5 + 3].Data(), VarManager::fgValues);
                  }
                }
              }
            }
            for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++) {
              AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
              if (!(cut.IsSelected(VarManager::fgValues))) // apply pair cuts
                continue;
              if (sign1 * sign2 < 0) {
                fHistMan->FillHistClass(histNames[ncuts + icut * fPairCuts.size() + iPairCut][0].Data(), VarManager::fgValues);
              } else {
                if (sign1 > 0) {
                  fHistMan->FillHistClass(histNames[ncuts + icut * fPairCuts.size() + iPairCut][1].Data(), VarManager::fgValues);
                } else {
                  fHistMan->FillHistClass(histNames[ncuts + icut * fPairCuts.size() + iPairCut][2].Data(), VarManager::fgValues);
                }
              }
            } // end loop (pair cuts)
          }
        } // end loop (cuts)
      } // end loop over pairs of track associations
    } // end loop over events
  }

  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;

  template <int TPairType>
  void runMCGenWithGrouping(MyEventsVtxCovSelected const& events, ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    uint32_t mcDecision = 0;
    int isig = 0;

    for (auto& mctrack : mcTracks) {
      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }
    // Fill Generated histograms taking into account selected collisions
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }

      for (auto& track : mcTracks) {
        if (track.reducedMCeventId() != event.reducedMCeventId()) {
          continue;
        }
        VarManager::FillTrackMC(mcTracks, track);
        auto track_raw = mcTracks.rawIteratorAt(track.globalIndex());
        mcDecision = 0;
        isig = 0;
        for (auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, track_raw)) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
            MCTruthTableEffi(VarManager::fgValues[VarManager::kMCPt], VarManager::fgValues[VarManager::kMCEta], VarManager::fgValues[VarManager::kMCY], VarManager::fgValues[VarManager::kMCPhi], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCVtxZ], VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);

            if (useMiniTree.fConfigMiniTree) {
              auto mcEvent = mcEvents.rawIteratorAt(track_raw.reducedMCeventId());
              dileptonMiniTreeGen(mcDecision, mcEvent.impactParameter(), track_raw.pt(), track_raw.eta(), track_raw.phi(), -999, -999, -999);
            }
          }
          isig++;
        }
      }
    } // end loop over reconstructed events

    if (fHasTwoProngGenMCsignals) {
      for (auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
        auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
        auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
        if (t1_raw.reducedMCeventId() == t2_raw.reducedMCeventId()) {
          for (auto& sig : fGenMCSignals) {
            if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required here
              continue;
            }
            if (sig->CheckSignal(true, t1_raw, t2_raw)) {
              VarManager::FillPairMC<TPairType>(t1, t2);
              fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig->GetName()), VarManager::fgValues);
            }
          }
        }
      }
    }
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }
      // CURRENTLY ONLY FOR 1-GENERATION 2-PRONG SIGNALS
      if (fHasTwoProngGenMCsignals) {
        auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reducedMCeventId());
        groupedMCTracks.bindInternalIndicesTo(&mcTracks);
        for (auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {
          auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
          auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
          if (t1_raw.reducedMCeventId() == t2_raw.reducedMCeventId()) {
            mcDecision = 0;
            isig = 0;
            for (auto& sig : fGenMCSignals) {
              if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required here
                continue;
              }
              if (sig->CheckSignal(true, t1_raw, t2_raw)) {
                mcDecision |= (static_cast<uint32_t>(1) << isig);
                VarManager::FillPairMC<TPairType>(t1, t2);
                fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s", sig->GetName()), VarManager::fgValues);
                // Fill also acceptance cut histograms if requested
                if (fUseMCGenAccCut) {
                  for (auto& accCut : fMCGenAccCuts) {
                    if (accCut.IsSelected(VarManager::fgValues)) {
                      fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s_%s", sig->GetName(), accCut.GetName()), VarManager::fgValues);
                    }
                  }
                }
                if (useMiniTree.fConfigMiniTree) {
                  // WARNING! To be checked
                  dileptonMiniTreeGen(mcDecision, -999, t1.pt(), t1.eta(), t1.phi(), t2.pt(), t2.eta(), t2.phi());
                }
              }
              isig++;
            }
          }
        }
      } // end loop over reconstructed events
    }
  }

  void processAllSkimmed(MyEventsVtxCovSelected const& events,
                         soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs, MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                         soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons,
                         ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(events, muonAssocsPerCollision, muonAssocs, muons, mcEvents, mcTracks);
    // Feature replaced by processMCGen
    /*if (fConfigMC.runMCGenPair) {
      runMCGen<VarManager::kDecayToEE>(mcEvents, mcTracks);
    }*/
    // runSameEventPairing<true, VarManager::kElectronMuon, gkEventFillMap, gkTrackFillMap>(event, tracks, muons);
  }

  void processBarrelOnlySkimmed(MyEventsVtxCovSelected const& events,
                                soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                MyBarrelTracksWithCovWithAmbiguities const& barrelTracks, ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks);
    runMCGenWithGrouping<VarManager::kDecayToEE>(events, mcEvents, mcTracks);
  }

  void processBarrelOnlyWithCollSkimmed(MyEventsVtxCovSelected const& events,
                                        soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                        MyBarrelTracksWithCovWithAmbiguitiesWithColl const& barrelTracks, ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCovWithColl>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks);
    runMCGenWithGrouping<VarManager::kDecayToEE>(events, mcEvents, mcTracks);
  }

  void processMuonOnlySkimmed(MyEventsVtxCovSelected const& events,
                              soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons, ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(events, muonAssocsPerCollision, muonAssocs, muons, mcEvents, mcTracks);
    // Feature replaced by processMCGen
    /* if (fConfigMC.runMCGenPair) {
       runMCGen<VarManager::kDecayToMuMu>(mcEvents, mcTracks);
     }*/
  }

  PresliceUnsorted<ReducedMCTracks> perReducedMcGenEvent = aod::reducedtrackMC::reducedMCeventId;

  void processMCGen(soa::Filtered<MyEventsVtxCovSelected> const& events, ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    // Fill Generated histograms taking into account all generated tracks
    uint32_t mcDecision = 0;
    int isig = 0;

    for (auto& mctrack : mcTracks) {
      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }

    // Fill Generated histograms taking into account selected collisions
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }
      VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
      VarManager::FillEvent<VarManager::ObjTypes::ReducedEventMC>(event.reducedMCevent(), VarManager::fgValues);
      // auto groupedMCTracks = mcTracks.sliceBy(perReducedMcGenEvent, event.reducedMCeventId());
      // groupedMCTracks.bindInternalIndicesTo(&mcTracks);
      // for (auto& track : groupedMCTracks) {
      for (auto& track : mcTracks) {
        if (track.reducedMCeventId() != event.reducedMCeventId()) {
          continue;
        }
        VarManager::FillTrackMC(mcTracks, track);
        auto track_raw = mcTracks.rawIteratorAt(track.globalIndex());
        // auto track_raw = groupedMCTracks.rawIteratorAt(track.globalIndex());
        mcDecision = 0;
        isig = 0;
        for (auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, track_raw)) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
            MCTruthTableEffi(VarManager::fgValues[VarManager::kMCPt], VarManager::fgValues[VarManager::kMCEta], VarManager::fgValues[VarManager::kMCY], VarManager::fgValues[VarManager::kMCPhi], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCVtxZ], VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);

            if (useMiniTree.fConfigMiniTree) {
              auto mcEvent = mcEvents.rawIteratorAt(track_raw.reducedMCeventId());
              dileptonMiniTreeGen(mcDecision, mcEvent.impactParameter(), track_raw.pt(), track_raw.eta(), track_raw.phi(), -999, -999, -999);
            }
          }
          isig++;
        }
      }
    } // end loop over reconstructed events
    if (fHasTwoProngGenMCsignals) {
      for (auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
        auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
        auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
        if (t1_raw.reducedMCeventId() == t2_raw.reducedMCeventId()) {
          for (auto& sig : fGenMCSignals) {
            if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required here
              continue;
            }
            if (sig->CheckSignal(true, t1_raw, t2_raw)) {
              VarManager::FillPairMC<VarManager::kDecayToMuMu>(t1, t2); // NOTE: This feature will only work for muons
              fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig->GetName()), VarManager::fgValues);
            }
          }
        }
      }
    }
    // Fill Generated PAIR histograms taking into account selected collisions
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }

      if (fHasTwoProngGenMCsignals) {
        for (auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
          if (t1.reducedMCeventId() != event.reducedMCeventId()) {
            continue;
          }
          if (t2.reducedMCeventId() != event.reducedMCeventId()) {
            continue;
          }
          auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
          auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
          if (t1_raw.reducedMCeventId() == t2_raw.reducedMCeventId()) {
            mcDecision = 0;
            isig = 0;
            for (auto& sig : fGenMCSignals) {
              if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required here
                continue;
              }
              if (sig->CheckSignal(true, t1_raw, t2_raw)) {
                mcDecision |= (static_cast<uint32_t>(1) << isig);
                VarManager::FillPairMC<VarManager::kDecayToMuMu>(t1, t2); // NOTE: This feature will only work for muons
                fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s", sig->GetName()), VarManager::fgValues);
                if (useMiniTree.fConfigMiniTree) {
                  // WARNING! To be checked
                  dileptonMiniTreeGen(mcDecision, -999, t1.pt(), t1.eta(), t1.phi(), t2.pt(), t2.eta(), t2.phi());
                }
              }
              isig++;
            }
          }
        }
      }
    } // end loop over reconstructed events
  }

  void processMCGenWithGrouping(soa::Filtered<MyEventsVtxCovSelected> const& events, ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    uint32_t mcDecision = 0;
    int isig = 0;

    for (auto& mctrack : mcTracks) {
      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }
    // Fill Generated histograms taking into account selected collisions
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }

      for (auto& track : mcTracks) {
        if (track.reducedMCeventId() != event.reducedMCeventId()) {
          continue;
        }
        VarManager::FillTrackMC(mcTracks, track);
        auto track_raw = mcTracks.rawIteratorAt(track.globalIndex());
        mcDecision = 0;
        isig = 0;
        for (auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, track_raw)) {
            mcDecision |= (static_cast<uint32_t>(1) << isig);
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
            MCTruthTableEffi(VarManager::fgValues[VarManager::kMCPt], VarManager::fgValues[VarManager::kMCEta], VarManager::fgValues[VarManager::kMCY], VarManager::fgValues[VarManager::kMCPhi], VarManager::fgValues[VarManager::kMCVz], VarManager::fgValues[VarManager::kMCVtxZ], VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);

            if (useMiniTree.fConfigMiniTree) {
              auto mcEvent = mcEvents.rawIteratorAt(track_raw.reducedMCeventId());
              dileptonMiniTreeGen(mcDecision, mcEvent.impactParameter(), track_raw.pt(), track_raw.eta(), track_raw.phi(), -999, -999, -999);
            }
          }
          isig++;
        }
      }
    } // end loop over reconstructed events
    if (fHasTwoProngGenMCsignals) {
      for (auto& [t1, t2] : combinations(mcTracks, mcTracks)) {
        auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
        auto t2_raw = mcTracks.rawIteratorAt(t2.globalIndex());
        if (t1_raw.reducedMCeventId() == t2_raw.reducedMCeventId()) {
          for (auto& sig : fGenMCSignals) {
            if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required here
              continue;
            }
            if (sig->CheckSignal(true, t1_raw, t2_raw)) {
              VarManager::FillPairMC<VarManager::kDecayToMuMu>(t1, t2); // NOTE: This feature will only work for muons
              fHistMan->FillHistClass(Form("MCTruthGenPair_%s", sig->GetName()), VarManager::fgValues);
            }
          }
        }
      }
    }
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }
      // CURRENTLY ONLY FOR 1-GENERATION 2-PRONG SIGNALS
      if (fHasTwoProngGenMCsignals) {
        auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reducedMCeventId());
        groupedMCTracks.bindInternalIndicesTo(&mcTracks);
        for (auto& [t1, t2] : combinations(groupedMCTracks, groupedMCTracks)) {
          auto t1_raw = groupedMCTracks.rawIteratorAt(t1.globalIndex());
          auto t2_raw = groupedMCTracks.rawIteratorAt(t2.globalIndex());
          if (t1_raw.reducedMCeventId() == t2_raw.reducedMCeventId()) {
            mcDecision = 0;
            isig = 0;
            for (auto& sig : fGenMCSignals) {
              if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required here
                continue;
              }
              if (sig->CheckSignal(true, t1_raw, t2_raw)) {
                mcDecision |= (static_cast<uint32_t>(1) << isig);
                VarManager::FillPairMC<VarManager::kDecayToMuMu>(t1, t2); // NOTE: This feature will only work for muons
                fHistMan->FillHistClass(Form("MCTruthGenPairSel_%s", sig->GetName()), VarManager::fgValues);
                if (useMiniTree.fConfigMiniTree) {
                  // WARNING! To be checked
                  dileptonMiniTreeGen(mcDecision, -999, t1.pt(), t1.eta(), t1.phi(), t2.pt(), t2.eta(), t2.phi());
                }
              }
              isig++;
            }
          }
        }
      } // end loop over reconstructed events
    }
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processAllSkimmed, "Run all types of pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlySkimmed, "Run barrel only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlyWithCollSkimmed, "Run barrel only pairing, with skimmed tracks and with collision information", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMuonOnlySkimmed, "Run muon only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMCGen, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMCGenWithGrouping, "Loop over MC particle stack (grouped MCTracks) and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", true);
};

// Run pairing for resonance with legs fulfilling separate cuts (asymmetric decay channel)
struct AnalysisAsymmetricPairing {

  Produces<aod::Ditracks> ditrackList;
  Produces<aod::DitracksExtra> ditrackExtraList;

  o2::base::MatLayerCylSet* fLUT = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  // Output objects
  OutputObj<THashList> fOutputList{"output"};

  // Configurables
  Configurable<std::string> fConfigLegCuts{"cfgLegCuts", "", "<leg-A-1>:<leg-B-1>[:<leg-C-1>],[<leg-A-2>:<leg-B-2>[:<leg-C-1>],...]"};
  Configurable<uint32_t> fConfigLegAFilterMask{"cfgLegAFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<uint32_t> fConfigLegBFilterMask{"cfgLegBFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<uint32_t> fConfigLegCFilterMask{"cfgLegCFilterMask", 0, "Filter mask corresponding to cuts in event-selection"};
  Configurable<std::string> fConfigCommonTrackCuts{"cfgCommonTrackCuts", "", "Comma separated list of cuts to be applied to all legs"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts"};
  Configurable<std::string> fConfigPairCutsJSON{"cfgPairCutsJSON", "", "Additional list of pair cuts in JSON format"};
  Configurable<bool> fConfigSkipAmbiguousIdCombinations{"cfgSkipAmbiguousIdCombinations", true, "Choose whether to skip pairs/triples which pass a stricter combination of cuts, e.g. KKPi triplets for D+ -> KPiPi"};

  Configurable<std::string> fConfigHistogramSubgroups{"cfgAsymmetricPairingHistogramsSubgroups", "barrel,vertexing", "Comma separated list of asymmetric-pairing histogram subgroups"};
  Configurable<bool> fConfigSameSignHistograms{"cfgSameSignHistograms", false, "Include same sign pair histograms for 2-prong decays"};
  // Configurable<bool> fConfigAmbiguousHistograms{"cfgAmbiguousHistograms", false, "Include separate histograms for pairs/triplets with ambiguous tracks"};
  Configurable<bool> fConfigReflectedHistograms{"cfgReflectedHistograms", false, "Include separate histograms for pairs which are reflections of previously counted pairs"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> fConfigGRPMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> fConfigUseRemoteField{"cfgUseRemoteField", false, "Choose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};

  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fConfigUseAbsDCA{"cfgUseAbsDCA", false, "Use absolute DCA minimization instead of chi^2 minimization in secondary vertexing"};
  Configurable<bool> fConfigPropToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
  Configurable<std::string> fConfigLutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};

  Configurable<std::string> fConfigMCGenSignals{"cfgBarrelMCGenSignals", "", "Comma separated list of MC signals (generated)"};
  Configurable<std::string> fConfigMCRecSignals{"cfgBarrelMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
  Configurable<std::string> fConfigMCRecSignalsJSON{"cfgMCRecSignalsJSON", "", "Additional list of MC signals (reconstructed) via JSON"};
  Configurable<std::string> fConfigMCGenSignalsJSON{"cfgMCGenSignalsJSON", "", "Comma separated list of MC signals (generated) via JSON"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;

  std::vector<AnalysisCompositeCut*> fPairCuts;
  int fNPairHistPrefixes;

  std::vector<MCSignal*> fRecMCSignals;
  std::vector<MCSignal*> fGenMCSignals;

  // Filter masks to find legs in BarrelTrackCuts table
  uint32_t fLegAFilterMask;
  uint32_t fLegBFilterMask;
  uint32_t fLegCFilterMask;
  // Maps tracking which combination of leg cuts the track cuts participate in
  std::map<int, uint32_t> fConstructedLegAFilterMasksMap;
  std::map<int, uint32_t> fConstructedLegBFilterMasksMap;
  std::map<int, uint32_t> fConstructedLegCFilterMasksMap;
  // Filter map for common track cuts
  uint32_t fCommonTrackCutMask;
  // Map tracking which common track cut the track cuts correspond to
  std::map<int, uint32_t> fCommonTrackCutFilterMasks;

  int fNLegCuts;
  int fNPairCuts = 0;
  int fNCommonTrackCuts;
  // vectors for cut names and signal names, for easy access when calling FillHistogramList()
  std::vector<TString> fLegCutNames;
  std::vector<TString> fPairCutNames;
  std::vector<TString> fCommonCutNames;
  std::vector<TString> fRecMCSignalNames;

  Filter eventFilter = aod::dqanalysisflags::isEventSelected > static_cast<uint32_t>(0);

  Preslice<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

  // Partitions for triplets and asymmetric pairs
  Partition<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> legACandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegAFilterMask) > static_cast<uint32_t>(0);
  Partition<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> legBCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegBFilterMask) > static_cast<uint32_t>(0);
  Partition<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> legCCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegCFilterMask) > static_cast<uint32_t>(0);

  // Map to track how many times a pair of tracks has been encountered
  std::map<std::pair<int32_t, int32_t>, int8_t> fPairCount;

  void init(o2::framework::InitContext& context)
  {
    bool isMCGen = context.mOptions.get<bool>("processMCGen") || context.mOptions.get<bool>("processMCGenWithEventSelection");
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Get the leg cut filter masks
    fLegAFilterMask = fConfigLegAFilterMask.value;
    fLegBFilterMask = fConfigLegBFilterMask.value;
    fLegCFilterMask = fConfigLegCFilterMask.value;

    // Get the pair cuts
    TString cutNamesStr = fConfigPairCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // Extra pair cuts via JSON
    TString addPairCutsStr = fConfigPairCutsJSON.value;
    if (addPairCutsStr != "") {
      std::vector<AnalysisCut*> addPairCuts = dqcuts::GetCutsFromJSON(addPairCutsStr.Data());
      for (auto& t : addPairCuts) {
        fPairCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
        cutNamesStr += Form(",%s", t->GetName());
      }
    }
    std::unique_ptr<TObjArray> objArrayPairCuts(cutNamesStr.Tokenize(","));
    fNPairCuts = objArrayPairCuts->GetEntries();
    for (int j = 0; j < fNPairCuts; j++) {
      fPairCutNames.push_back(objArrayPairCuts->At(j)->GetName());
    }

    // Setting the MC rec signal names
    TString sigNamesStr = fConfigMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        fRecMCSignals.push_back(sig);
      }
    }
    // Add the reco MCSignals from the JSON config
    TString addMCSignalsStr = fConfigMCRecSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != 2 && mcIt->GetNProngs() != 3) {
          LOG(fatal) << "Signal at reconstructed level requested (" << mcIt->GetName() << ") " << "does not have 2 or 3 prongs! Fix it";
        }
        fRecMCSignals.push_back(mcIt);
        sigNamesStr += Form(",%s", mcIt->GetName());
      }
    }
    // Put all the reco MCSignal names in the vector for histogram naming
    std::unique_ptr<TObjArray> objArrayRecMCSignals(sigNamesStr.Tokenize(","));
    for (int i = 0; i < objArrayRecMCSignals->GetEntries(); i++) {
      fRecMCSignalNames.push_back(objArrayRecMCSignals->At(i)->GetName());
    }

    // Get the barrel track selection cuts
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
    std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
    // Get the common leg cuts
    int commonCutIdx;
    TString commonNamesStr = fConfigCommonTrackCuts.value;
    if (!commonNamesStr.IsNull()) { // if common track cuts
      std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
      fNCommonTrackCuts = objArrayCommon->GetEntries();
      for (int icut = 0; icut < fNCommonTrackCuts; ++icut) {
        commonCutIdx = objArray->IndexOf(objArrayCommon->At(icut));
        if (commonCutIdx >= 0) {
          fCommonTrackCutMask |= static_cast<uint32_t>(1) << objArray->IndexOf(objArrayCommon->At(icut));
          fCommonTrackCutFilterMasks[icut] = static_cast<uint32_t>(1) << objArray->IndexOf(objArrayCommon->At(icut));
          fCommonCutNames.push_back(objArrayCommon->At(icut)->GetName());
        } else {
          LOGF(fatal, "Common track cut %s was not calculated upstream. Check the config!", objArrayCommon->At(icut)->GetName());
        }
      }
    }
    // Check that the leg cut masks make sense
    if (static_cast<int>(std::floor(TMath::Log2(fLegAFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "fConfigLegAFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(TMath::Log2(fLegAFilterMask))) + 1, objArray->GetEntries());
    }
    if (static_cast<int>(std::floor(TMath::Log2(fLegBFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "fConfigLegBFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(TMath::Log2(fLegBFilterMask))) + 1, objArray->GetEntries());
    }
    if (static_cast<int>(std::floor(TMath::Log2(fLegCFilterMask))) + 1 > objArray->GetEntries()) {
      LOGF(fatal, "fConfigLegCFilterMask has highest bit at position %d, but track-selection only has %d cuts!", static_cast<int>(std::floor(TMath::Log2(fLegCFilterMask))) + 1, objArray->GetEntries());
    }

    // Get the cuts defining the legs
    uint32_t fConstructedLegAFilterMask = 0;
    uint32_t fConstructedLegBFilterMask = 0;
    uint32_t fConstructedLegCFilterMask = 0;
    TString legCutsStr = fConfigLegCuts.value;
    std::unique_ptr<TObjArray> objArrayLegs(legCutsStr.Tokenize(","));
    if (objArrayLegs->GetEntries() == 0 && !isMCGen) {
      LOG(fatal) << "No cuts defining legs. Check the config!";
    }
    fNLegCuts = objArrayLegs->GetEntries();
    std::vector<bool> isThreeProng;
    int legAIdx;
    int legBIdx;
    int legCIdx;
    // Loop over leg defining cuts
    for (int icut = 0; icut < fNLegCuts; ++icut) {
      TString legsStr = objArrayLegs->At(icut)->GetName();
      std::unique_ptr<TObjArray> legs(legsStr.Tokenize(":"));
      if (legs->GetEntries() == 3) {
        isThreeProng.push_back(true);
      } else if (legs->GetEntries() == 2) {
        isThreeProng.push_back(false);
      } else {
        LOGF(fatal, "Leg cuts %s has the wrong format and could not be parsed!", legsStr.Data());
        continue;
      }
      // Find leg cuts in the track selection cuts
      legAIdx = objArray->IndexOf(legs->At(0));
      if (legAIdx >= 0) {
        fConstructedLegAFilterMask |= static_cast<uint32_t>(1) << legAIdx;
        fConstructedLegAFilterMasksMap[icut] |= static_cast<uint32_t>(1) << legAIdx;
      } else {
        LOGF(fatal, "Leg A cut %s was not calculated upstream. Check the config!", legs->At(0)->GetName());
        continue;
      }
      legBIdx = objArray->IndexOf(legs->At(1));
      if (legBIdx >= 0) {
        fConstructedLegBFilterMask |= static_cast<uint32_t>(1) << legBIdx;
        fConstructedLegBFilterMasksMap[icut] |= static_cast<uint32_t>(1) << legBIdx;
      } else {
        LOGF(fatal, "Leg B cut %s was not calculated upstream. Check the config!", legs->At(1)->GetName());
        continue;
      }
      if (isThreeProng[icut]) {
        legCIdx = objArray->IndexOf(legs->At(2));
        if (legCIdx >= 0) {
          fConstructedLegCFilterMask |= static_cast<uint32_t>(1) << legCIdx;
          fConstructedLegCFilterMasksMap[icut] |= static_cast<uint32_t>(1) << legCIdx;
        } else {
          LOGF(fatal, "Leg C cut %s was not calculated upstream. Check the config!", legs->At(2)->GetName());
          continue;
        }
      }
      // Leg cut config is fine, store the leg cut name in a vector
      fLegCutNames.push_back(legsStr);

      // Define histogram and histogram directory names
      if (isThreeProng[icut]) {
        DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s", legsStr.Data()), fConfigHistogramSubgroups.value.data());
        if (fConfigQA) {
          DefineHistograms(fHistMan, Form("TripletsBarrelSE_ambiguous_%s", legsStr.Data()), fConfigHistogramSubgroups.value.data());
        }

        std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
          DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName()), fConfigHistogramSubgroups.value.data());
        }

        TString cutNamesStr = fConfigPairCuts.value;
        if (!cutNamesStr.IsNull()) { // if pair cuts
          std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
          fNPairCuts = objArrayPair->GetEntries();
          for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
            DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), objArrayPair->At(iPairCut)->GetName()), fConfigHistogramSubgroups.value.data());
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName()), fConfigHistogramSubgroups.value.data());
            } // end loop (common cuts)
          } // end loop (pair cuts)
        } // end if (pair cuts)

        // TODO: assign hist directories for the MC matched triplets for each (leg cut combo,MCsignal) combination
        if (!sigNamesStr.IsNull()) {
          for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
            auto sig = fRecMCSignals.at(isig);
            DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());

            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
            }

            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s", legsStr.Data(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
                for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                  DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
                } // end loop (common cuts)
              } // end loop (pair cuts)
            } // end if (pair cuts)
          } // end loop over MC signals
        } // end if (MC signals)
      } else {
        std::vector<TString> pairHistPrefixes = {"PairsBarrelSEPM"};
        if (fConfigSameSignHistograms.value) {
          pairHistPrefixes.push_back("PairsBarrelSEPP");
          pairHistPrefixes.push_back("PairsBarrelSEMM");
        }
        fNPairHistPrefixes = pairHistPrefixes.size();

        for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
          DefineHistograms(fHistMan, Form("%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()), fConfigHistogramSubgroups.value.data());
        }
        if (fConfigQA) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            DefineHistograms(fHistMan, Form("%s_ambiguous_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()), fConfigHistogramSubgroups.value.data());
          }
        }
        if (fConfigReflectedHistograms.value) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            DefineHistograms(fHistMan, Form("%s_reflected_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data()), fConfigHistogramSubgroups.value.data());
          }
        }

        std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
          for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
            DefineHistograms(fHistMan, Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName()), fConfigHistogramSubgroups.value.data());
          }
        }

        if (!cutNamesStr.IsNull()) { // if pair cuts
          std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
          fNPairCuts = objArrayPair->GetEntries();
          for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
            for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
              DefineHistograms(fHistMan, Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayPair->At(iPairCut)->GetName()), fConfigHistogramSubgroups.value.data());
            }
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                DefineHistograms(fHistMan, Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName()), fConfigHistogramSubgroups.value.data());
              }
            } // end loop (common cuts)
          } // end loop (pair cuts)
        } // end if (pair cuts)

        // assign hist directories for the MC matched triplets for each (leg cut combo,MCsignal) combination
        if (!sigNamesStr.IsNull()) {
          for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) {
            auto sig = fRecMCSignals.at(isig);
            for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
              DefineHistograms(fHistMan, Form("%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());
            }
            if (fConfigReflectedHistograms.value) {
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                DefineHistograms(fHistMan, Form("%s_reflected_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());
              }
            }

            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                DefineHistograms(fHistMan, Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
              }
            }

            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) { // loop over pair cuts
                for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                  DefineHistograms(fHistMan, Form("%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
                }
                for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                  for (int iPrefix = 0; iPrefix < fNPairHistPrefixes; ++iPrefix) {
                    DefineHistograms(fHistMan, Form("%s_%s_%s_%s_%s", pairHistPrefixes[iPrefix].Data(), legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName(), objArrayPair->At(iPairCut)->GetName(), sig->GetName()), fConfigHistogramSubgroups.value.data());
                  }
                } // end loop (common cuts)
              } // end loop (pair cuts)
            } // end if (pair cuts)
          } // end loop over MC signals
        } // end if (MC signals)
      }
    }

    // Add histogram classes for each specified MCsignal at the generator level
    // TODO: create a std::vector of hist classes to be used at Fill time, to avoid using Form in the process function
    TString sigGenNamesStr = fConfigMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 1) { // NOTE: 1-prong signals required
          fGenMCSignals.push_back(sig);
          DefineHistograms(fHistMan, Form("MCTruthGen_%s;", sig->GetName()), fConfigHistogramSubgroups.value.data());    // TODO: Add these names to a std::vector to avoid using Form in the process function
          DefineHistograms(fHistMan, Form("MCTruthGenSel_%s;", sig->GetName()), fConfigHistogramSubgroups.value.data()); // TODO: Add these names to a std::vector to avoid using Form in the process function
        }
      }
    }

    // Add the gen MCSignals from the JSON config
    addMCSignalsStr = fConfigMCGenSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() == 1) {
          fGenMCSignals.push_back(mcIt);
          DefineHistograms(fHistMan, Form("MCTruthGen_%s;", mcIt->GetName()), fConfigHistogramSubgroups.value.data());    // TODO: Add these names to a std::vector to avoid using Form in the process function
          DefineHistograms(fHistMan, Form("MCTruthGenSel_%s;", mcIt->GetName()), fConfigHistogramSubgroups.value.data()); // TODO: Add these names to a std::vector to avoid using Form in the process function
        }
      }
    }

    // Make sure the leg cuts are covered by the configured filter masks
    if (fLegAFilterMask != fConstructedLegAFilterMask) {
      LOGF(fatal, "cfgLegAFilterMask (%d) is not equal to the mask constructed by the cuts specified in cfgLegCuts (%d)!", fLegAFilterMask, fConstructedLegAFilterMask);
    }
    if (fLegBFilterMask != fConstructedLegBFilterMask) {
      LOGF(fatal, "cfgLegBFilterMask (%d) is not equal to the mask constructed by the cuts specified in cfgLegCuts (%d)!", fLegBFilterMask, fConstructedLegBFilterMask);
    }
    if (fLegCFilterMask != fConstructedLegCFilterMask) {
      LOGF(fatal, "cfgLegCFilterMask (%d) is not equal to the mask constructed by the cuts specified in cfgLegCuts (%d)!", fLegCFilterMask, fConstructedLegCFilterMask);
    }
    // Make sure only pairs or only triplets of leg cuts were given
    int tripletCheckSum = std::count(isThreeProng.begin(), isThreeProng.end(), true);
    if (tripletCheckSum != 0 && tripletCheckSum != fNLegCuts) {
      LOGF(fatal, "A mix of pairs and triplets was given as leg cuts. Check your config!");
    }

    fCurrentRun = 0;

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();

    fLUT = o2::base::MatLayerCylSet::rectifyPtrFromFile(fCCDB->get<o2::base::MatLayerCylSet>(fConfigLutPath));
    VarManager::SetupMatLUTFwdDCAFitter(fLUT);

    dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON
    VarManager::SetUseVars(fHistMan->GetUsedVars());                                       // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void initParamsFromCCDB(uint64_t timestamp, bool isTriplets)
  {
    if (fConfigUseRemoteField.value) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigGRPMagPath, timestamp);
      float magField = 0.0;
      if (grpmag != nullptr) {
        magField = grpmag->getNominalL3Field();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", timestamp);
      }
      if (isTriplets) {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupThreeProngKFParticle(magField);
        } else {
          VarManager::SetupThreeProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value);
        }
      } else {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(magField);
        } else {
          VarManager::SetupTwoProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value); // TODO: get these parameters from Configurables
        }
      }
    } else {
      if (isTriplets) {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupThreeProngKFParticle(fConfigMagField.value);
        } else {
          VarManager::SetupThreeProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value);
        }
      } else {
        if (fConfigUseKFVertexing.value) {
          VarManager::SetupTwoProngKFParticle(fConfigMagField.value);
        } else {
          VarManager::SetupTwoProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, fConfigUseAbsDCA.value); // TODO: get these parameters from Configurables
        }
      }
    }
  }

  // Template function to run same event pairing with asymmetric pairs (e.g. kaon-pion)
  template <bool TTwoProngFitter, int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTrackAssocs, typename TTracks>
  void runAsymmetricPairing(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& /*assocs*/, TTracks const& /*tracks*/, ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& /*mcTracks*/)
  {
    fPairCount.clear();

    if (events.size() > 0) { // Additional protection to avoid crashing of events.begin().runNumber()
      if (fCurrentRun != events.begin().runNumber()) {
        initParamsFromCCDB(events.begin().timestamp(), false);
        fCurrentRun = events.begin().runNumber();
      }
    }

    int sign1 = 0;
    int sign2 = 0;
    uint32_t mcDecision = 0;
    ditrackList.reserve(1);
    ditrackExtraList.reserve(1);

    constexpr bool trackHasCov = ((TTrackFillMap & VarManager::ObjTypes::TrackCov) > 0 || (TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelCov) > 0);

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event, VarManager::fgValues);

      auto groupedLegAAssocs = legACandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegAAssocs.size() == 0) {
        continue;
      }
      auto groupedLegBAssocs = legBCandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegBAssocs.size() == 0) {
        continue;
      }

      for (auto& [a1, a2] : combinations(soa::CombinationsFullIndexPolicy(groupedLegAAssocs, groupedLegBAssocs))) {

        uint32_t twoTrackFilter = 0;
        uint32_t twoTrackCommonFilter = 0;
        uint32_t pairFilter = 0;
        bool isPairIdWrong = false;
        for (int icut = 0; icut < fNLegCuts; ++icut) {
          // Find leg pair definitions both candidates participate in
          if ((a1.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) && (a2.isBarrelSelected_raw() & fConstructedLegBFilterMasksMap[icut])) {
            twoTrackFilter |= static_cast<uint32_t>(1) << icut;
            // If the supposed pion passes a kaon cut, this is a K+K-. Skip it.
            if (TPairType == VarManager::kDecayToKPi && fConfigSkipAmbiguousIdCombinations.value) {
              if (a2.isBarrelSelected_raw() & fLegAFilterMask) {
                isPairIdWrong = true;
              }
            }
          }
        }

        if (!twoTrackFilter || isPairIdWrong) {
          continue;
        }

        // Find common track cuts both candidates pass
        twoTrackCommonFilter |= a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & fCommonTrackCutMask;

        auto t1 = a1.template reducedtrack_as<TTracks>();
        auto t2 = a2.template reducedtrack_as<TTracks>();

        // Avoid self-pairs
        if (t1.globalIndex() == t2.globalIndex()) {
          continue;
        }

        bool isReflected = false;
        std::pair<int32_t, int32_t> trackIds(t1.globalIndex(), t2.globalIndex());
        if (fPairCount.find(trackIds) != fPairCount.end()) {
          // Double counting is possible due to track-collision ambiguity. Skip pairs which were counted before
          fPairCount[trackIds] += 1;
          continue;
        }
        if (fPairCount.find(std::pair(trackIds.second, trackIds.first)) != fPairCount.end()) {
          isReflected = true;
        }
        fPairCount[trackIds] += 1;

        sign1 = t1.sign();
        sign2 = t2.sign();
        // store the ambiguity number of the two dilepton legs in the last 4 digits of the two-track filter
        if (t1.barrelAmbiguityInBunch() > 1 || t1.barrelAmbiguityOutOfBunch() > 1) {
          twoTrackFilter |= static_cast<uint32_t>(1) << 30;
        }
        if (t2.barrelAmbiguityInBunch() > 1 || t2.barrelAmbiguityOutOfBunch() > 1) {
          twoTrackFilter |= static_cast<uint32_t>(1) << 31;
        }

        // run MC matching for this pair
        int isig = 0;
        mcDecision = 0;
        for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
          if (t1.has_reducedMCTrack() && t2.has_reducedMCTrack()) {
            VarManager::FillPairMC<TPairType>(t1.reducedMCTrack(), t2.reducedMCTrack());
            if ((*sig)->CheckSignal(true, t1.reducedMCTrack(), t2.reducedMCTrack())) {
              mcDecision |= static_cast<uint32_t>(1) << isig;
            }
          }
        } // end loop over MC signals

        VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
        if constexpr (TTwoProngFitter) {
          VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigPropToPCA);
        }

        // Fill histograms
        bool isAmbi = false;
        for (int icut = 0; icut < fNLegCuts; icut++) {
          if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
            isAmbi = (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31));
            if (sign1 * sign2 < 0) {                                                                                // +- pairs
              fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s", fLegCutNames[icut].Data()), VarManager::fgValues); // reconstructed, unmatched
              if (isAmbi && fConfigQA) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
              }
              if (isReflected && fConfigReflectedHistograms.value) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_reflected_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
              }
            } else if (fConfigSameSignHistograms.value) {
              if (sign1 > 0) { // ++ pairs
                fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                if (isAmbi && fConfigQA) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
                if (isReflected && fConfigReflectedHistograms.value) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_reflected_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
              } else { // -- pairs
                fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                if (isAmbi && fConfigQA) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
                if (isReflected && fConfigReflectedHistograms) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_reflected_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
              }
            }
            for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
              if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                if (sign1 * sign2 < 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                  if (isReflected && fConfigReflectedHistograms.value) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPM_reflected_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                  }
                } else if (fConfigSameSignHistograms.value) {
                  if (sign1 > 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    if (isReflected && fConfigReflectedHistograms.value) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPP_reflected_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    }
                  } else {
                    fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    if (isReflected && fConfigReflectedHistograms.value) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEMM_reflected_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    }
                  }
                }
              }
            }
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
              if (twoTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
                if (sign1 * sign2 < 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), VarManager::fgValues);
                } else if (fConfigSameSignHistograms.value) {
                  if (sign1 > 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), VarManager::fgValues);
                  } else {
                    fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), VarManager::fgValues);
                  }
                }
                for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                  if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                    if (sign1 * sign2 < 0) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    } else if (fConfigSameSignHistograms.value) {
                      if (sign1 > 0) {
                        fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                      } else {
                        fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                      }
                    }
                  }
                }
              }
            } // end loop (common cuts)
            int iPairCut = 0;
            for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, iPairCut++) {
              if (!((*cut)->IsSelected(VarManager::fgValues))) // apply pair cuts
                continue;
              pairFilter |= (static_cast<uint32_t>(1) << iPairCut);
              // Histograms with pair cuts
              if (sign1 * sign2 < 0) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
              } else if (fConfigSameSignHistograms.value) {
                if (sign1 > 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                } else {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                }
              }
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                  if (sign1 * sign2 < 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                  } else if (fConfigSameSignHistograms.value) {
                    if (sign1 > 0) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    } else {
                      fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                    }
                  }
                }
              }
              // Histograms with pair cuts and common track cuts
              for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                if (twoTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
                  if (sign1 * sign2 < 0) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                  } else if (fConfigSameSignHistograms.value) {
                    if (sign1 > 0) {
                      fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                    } else {
                      fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
                    }
                  }
                  for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                    if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                      if (sign1 * sign2 < 0) {
                        fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                      } else if (fConfigSameSignHistograms.value) {
                        if (sign1 > 0) {
                          fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                        } else {
                          fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues);
                        }
                      }
                    }
                  }
                }
              }
            } // end loop (pair cuts)
          }
        } // end loop (cuts)
        ditrackList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                    VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                    t1.sign() + t2.sign(), twoTrackFilter, pairFilter, twoTrackCommonFilter);
        if constexpr (trackHasCov && TTwoProngFitter) {
          ditrackExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
        }
      } // end inner assoc loop (leg A)
    } // end event loop
  }

  // Template function to run same event triplets (e.g. D+->K-pi+pi+)
  template <bool TThreeProngFitter, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTrackAssocs, typename TTracks>
  void runThreeProng(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& /*assocs*/, TTracks const& tracks, ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& /*mcTracks*/, VarManager::PairCandidateType tripletType)
  {
    if (events.size() > 0) { // Additional protection to avoid crashing of events.begin().runNumber()
      if (fCurrentRun != events.begin().runNumber()) {
        initParamsFromCCDB(events.begin().timestamp(), true);
        fCurrentRun = events.begin().runNumber();
      }
    }

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event, VarManager::fgValues);

      auto groupedLegAAssocs = legACandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegAAssocs.size() == 0) {
        continue;
      }
      auto groupedLegBAssocs = legBCandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegBAssocs.size() == 0) {
        continue;
      }
      auto groupedLegCAssocs = legCCandidateAssocs.sliceBy(preslice, event.globalIndex());
      if (groupedLegCAssocs.size() == 0) {
        continue;
      }

      // Based on triplet type, make suitable combinations of the partitions
      if (tripletType == VarManager::kTripleCandidateToPKPi) {
        for (auto& [a1, a2, a3] : combinations(soa::CombinationsFullIndexPolicy(groupedLegAAssocs, groupedLegBAssocs, groupedLegCAssocs))) {
          readTriplet<TThreeProngFitter, TEventFillMap, TTrackFillMap>(a1, a2, a3, tracks, event, tripletType);
        }
      } else if (tripletType == VarManager::kTripleCandidateToKPiPi) {
        for (auto& a1 : groupedLegAAssocs) {
          for (auto& [a2, a3] : combinations(groupedLegBAssocs, groupedLegCAssocs)) {
            readTriplet<TThreeProngFitter, TEventFillMap, TTrackFillMap>(a1, a2, a3, tracks, event, tripletType);
          }
        }
      } else {
        LOG(fatal) << "Given tripletType not recognized. Don't know how to make combinations!" << endl;
      }
    } // end event loop
  }

  // Helper function to process triplet
  template <bool TThreeProngFitter, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TTrackAssoc, typename TTracks, typename TEvent>
  void readTriplet(TTrackAssoc const& a1, TTrackAssoc const& a2, TTrackAssoc const& a3, TTracks const& /*tracks*/, TEvent const& event, VarManager::PairCandidateType tripletType)
  {
    uint32_t mcDecision = 0;

    uint32_t threeTrackFilter = 0;
    uint32_t threeTrackCommonFilter = 0;
    for (int icut = 0; icut < fNLegCuts; ++icut) {
      // Find out which leg cut combinations the triplet passes
      if ((a1.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) && (a2.isBarrelSelected_raw() & fConstructedLegBFilterMasksMap[icut]) && (a3.isBarrelSelected_raw() & fConstructedLegCFilterMasksMap[icut])) {
        threeTrackFilter |= (static_cast<uint32_t>(1) << icut);
        if (tripletType == VarManager::kTripleCandidateToPKPi && fConfigSkipAmbiguousIdCombinations.value) {
          // Check if the supposed pion passes as a proton or kaon, if so, skip this triplet. It is pKp or pKK.
          if ((a3.isBarrelSelected_raw() & fLegAFilterMask) || (a3.isBarrelSelected_raw() & fLegBFilterMask)) {
            return;
          }
          // Check if the supposed kaon passes as a proton, if so, skip this triplet. It is ppPi.
          if (a2.isBarrelSelected_raw() & fLegAFilterMask) {
            return;
          }
        }
        if (tripletType == VarManager::kTripleCandidateToKPiPi && fConfigSkipAmbiguousIdCombinations.value) {
          // Check if one of the supposed pions pass as a kaon, if so, skip this triplet. It is KKPi.
          if ((a2.isBarrelSelected_raw() & fLegAFilterMask) || (a3.isBarrelSelected_raw() & fLegAFilterMask)) {
            return;
          }
        }
      }
    }
    if (!threeTrackFilter) {
      return;
    }

    // Find common track cuts all candidates pass
    threeTrackCommonFilter |= a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & a3.isBarrelSelected_raw() & fCommonTrackCutMask;

    auto t1 = a1.template reducedtrack_as<TTracks>();
    auto t2 = a2.template reducedtrack_as<TTracks>();
    auto t3 = a3.template reducedtrack_as<TTracks>();

    // Avoid self-pairs
    if (t1 == t2 || t1 == t3 || t2 == t3) {
      return;
    }
    // Check charge
    if (tripletType == VarManager::kTripleCandidateToKPiPi) {
      if (!((t1.sign() == -1 && t2.sign() == 1 && t3.sign() == 1) || (t1.sign() == 1 && t2.sign() == -1 && t3.sign() == -1))) {
        return;
      }
    }
    if (tripletType == VarManager::kTripleCandidateToPKPi) {
      if (!((t1.sign() == 1 && t2.sign() == -1 && t3.sign() == 1) || (t1.sign() == -1 && t2.sign() == 1 && t3.sign() == -1))) {
        return;
      }
    }

    // store the ambiguity of the three legs in the last 3 digits of the two-track filter
    if (t1.barrelAmbiguityInBunch() > 1 || t1.barrelAmbiguityOutOfBunch() > 1) {
      threeTrackFilter |= (static_cast<uint32_t>(1) << 29);
    }
    if (t2.barrelAmbiguityInBunch() > 1 || t2.barrelAmbiguityOutOfBunch() > 1) {
      threeTrackFilter |= (static_cast<uint32_t>(1) << 30);
    }
    if (t3.barrelAmbiguityInBunch() > 1 || t3.barrelAmbiguityOutOfBunch() > 1) {
      threeTrackFilter |= (static_cast<uint32_t>(1) << 31);
    }

    // run MC matching for this triplet
    int isig = 0;
    mcDecision = 0;
    for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
      if (t1.has_reducedMCTrack() && t2.has_reducedMCTrack() && t3.has_reducedMCTrack()) {
        if ((*sig)->CheckSignal(true, t1.reducedMCTrack(), t2.reducedMCTrack(), t3.reducedMCTrack())) {
          mcDecision |= (static_cast<uint32_t>(1) << isig);
        }
      }
    } // end loop over MC signals

    VarManager::FillTriple(t1, t2, t3, VarManager::fgValues, tripletType);
    if constexpr (TThreeProngFitter) {
      VarManager::FillTripletVertexing<TEventFillMap, TTrackFillMap>(event, t1, t2, t3, tripletType);
    }

    // Fill histograms
    bool isAmbi = false;
    for (int icut = 0; icut < fNLegCuts; icut++) {
      isAmbi = (threeTrackFilter & (static_cast<uint32_t>(1) << 29)) || (threeTrackFilter & (static_cast<uint32_t>(1) << 30)) || (threeTrackFilter & (static_cast<uint32_t>(1) << 31));
      if (threeTrackFilter & (static_cast<uint32_t>(1) << icut)) {
        fHistMan->FillHistClass(Form("TripletsBarrelSE_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
        for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
          if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
            fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues); // matched signal
          }
        } // end loop (MC signals)
        if (fConfigQA && isAmbi) {
          fHistMan->FillHistClass(Form("TripletsBarrelSE_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
        }
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
          if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
            fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), VarManager::fgValues);
            for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
              if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues); // matched signal
              }
            } // end loop (MC signals)
          }
        } // end loop (common cuts)
        int iPairCut = 0;
        for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, iPairCut++) {
          if (!((*cut)->IsSelected(VarManager::fgValues))) // apply pair cuts
            continue;
          // Histograms with pair cuts
          fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
          for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
            if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
              fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues); // matched signal
            }
          } // end loop (MC signals)
          // Histograms with pair cuts and common track cuts
          for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
            if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
              fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
              for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
                if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                  fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fRecMCSignalNames[isig].Data()), VarManager::fgValues); // matched signal
                }
              } // end loop (MC signals)
            }
          }
        } // end loop (pair cuts)
      }
    } // end loop (cuts)
  }

  void processKaonPionSkimmed(MyEventsVtxCovSelected const& events,
                              soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                              MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                              ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    runAsymmetricPairing<true, VarManager::kDecayToKPi, gkEventFillMapWithCov, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks);
  }

  void processKaonPionSkimmedMultExtra(MyEventsVtxCovSelectedMultExtra const& events,
                                       soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                                       MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                                       ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    runAsymmetricPairing<true, VarManager::kDecayToKPi, gkEventFillMapWithCovMultExtra, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks);
  }

  void processKaonPionPionSkimmed(MyEventsVtxCovSelected const& events,
                                  soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                                  MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                                  ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    runThreeProng<true, gkEventFillMapWithCov, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, mcEvents, mcTracks, VarManager::kTripleCandidateToKPiPi);
  }

  void processMCGen(ReducedMCTracks const& mcTracks)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // group all the MC tracks which belong to the MC event corresponding to the current reconstructed event
    // auto groupedMCTracks = tracksMC.sliceBy(aod::reducedtrackMC::reducedMCeventId, event.reducedMCevent().globalIndex());
    for (auto& mctrack : mcTracks) {

      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }
  }

  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;

  void processMCGenWithEventSelection(soa::Filtered<MyEventsVtxCovSelected> const& events,
                                      ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& mcTracks)
  {
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }

      auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reducedMCeventId());
      groupedMCTracks.bindInternalIndicesTo(&mcTracks);
      for (auto& track : groupedMCTracks) {

        VarManager::FillTrackMC(mcTracks, track);

        auto track_raw = groupedMCTracks.rawIteratorAt(track.globalIndex());
        for (auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, track_raw)) {
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
          }
        }
      }
    } // end loop over reconstructed events
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionSkimmed, "Run kaon pion pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionSkimmedMultExtra, "Run kaon pion pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionPionSkimmed, "Run kaon pion pion triplets, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processMCGen, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processMCGenWithEventSelection, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", true);
};

// Combines dileptons with barrel or muon tracks for either resonance or correlation analyses
// Dileptons produced with all the selection cuts specified in the same-event pairing task are combined with the
//   tracks passing the fConfigTrackCut cut. The dileptons cuts from the same-event pairing task are auto-detected
struct AnalysisDileptonTrack {
  Produces<aod::BmesonCandidates> BmesonsTable;
  Produces<aod::JPsiMuonCandidates> DileptonTrackTable;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "kaonPID", "Comma separated list of track cuts to be correlated with the dileptons"};
  Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 2.8, "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 3.2, "High mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonpTCut{"cfgDileptonpTCut", 0.0, "pT cut for dileptons used in the triplet vertexing"};
  Configurable<float> fConfigDileptonLxyCut{"cfgDileptonLxyCut", 0.0, "Lxy cut for dileptons used in the triplet vertexing"};
  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<float> fConfigDileptonLowpTCut{"cfgDileptonLowpTCut", 0.0, "Low pT cut for dileptons used in the triplet vertexing"};
  Configurable<float> fConfigDileptonHighpTCut{"cfgDileptonHighpTCut", 1E5, "High pT cut for dileptons used in the triplet vertexing"};
  Configurable<float> fConfigDileptonRapCutAbs{"cfgDileptonRapCutAbs", 1.0, "Rap cut for dileptons used in the triplet vertexing"};
  Configurable<std::string> fConfigHistogramSubgroups{"cfgDileptonTrackHistogramsSubgroups", "invmass,vertexing", "Comma separated list of dilepton-track histogram subgroups"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  Configurable<bool> fConfigUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<std::string> fConfigGRPmagPath{"cfgGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};
  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> fConfigGeoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<std::string> fConfigMCRecSignals{"cfgMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};
  Configurable<std::string> fConfigMCGenSignals{"cfgMCGenSignals", "", "Comma separated list of MC signals (generated)"};
  Configurable<std::string> fConfigMCRecSignalsJSON{"cfgMCRecSignalsJSON", "", "Additional list of MC signals (reconstructed) via JSON"};
  Configurable<std::string> fConfigMCGenSignalsJSON{"cfgMCGenSignalsJSON", "", "Comma separated list of MC signals (generated) via JSON"};
  Configurable<int> fConfigMCGenSignalDileptonLegPos{"cfgMCGenSignalDileptonLegPos", 0, "generator level positive dilepton leg signal (bit number according to table-maker)"};
  Configurable<int> fConfigMCGenSignalDileptonLegNeg{"cfgMCGenSignalDileptonLegNeg", 0, "generator level negative dilepton leg signal (bit number according to table-maker)"};
  Configurable<int> fConfigMCGenSignalHadron{"cfgMCGenSignalHadron", 0, "generator level associated hadron signal (bit number according to table-maker)"};
  Configurable<float> fConfigMCGenDileptonLegPtMin{"cfgMCGenDileptonLegPtMin", 1.0f, "minimum pt for the dilepton leg"};
  Configurable<float> fConfigMCGenHadronPtMin{"cfgMCGenHadronPtMin", 1.0f, "minimum pt for the hadron"};
  Configurable<float> fConfigMCGenDileptonLegEtaAbs{"cfgMCGenDileptonLegEtaAbs", 0.9f, "eta abs range for the dilepton leg"};
  Configurable<float> fConfigMCGenHadronEtaAbs{"cfgMCGenHadronEtaAbs", 0.9f, "eta abs range for the hadron"};
  Configurable<bool> fConfigUseMCRapcut{"cfgUseMCRapcut", false, "Use Rap cut for dileptons used in the triplet vertexing(reconstructed)"};
  Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 5, "Event mixing pool depth"};

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.
  int fNCuts;
  int fNLegCuts;
  int fNPairCuts;
  int fNCommonTrackCuts;
  std::map<int, int> fCommonTrackCutMap;
  uint32_t fTrackCutBitMap; // track cut bit mask to be used in the selection of tracks associated with dileptons
  // vector for single-lepton and track cut names for easy access when calling FillHistogramList()
  std::vector<TString> fTrackCutNames;
  std::vector<TString> fLegCutNames;
  // vector for pair cut names, used mainly for pairs built via the asymmetric pairing task
  std::vector<TString> fPairCutNames;
  std::vector<TString> fCommonPairCutNames;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  // TODO: The filter expressions seem to always use the default value of configurables, not the values from the actual configuration file
  Filter eventFilter = aod::dqanalysisflags::isEventSelected > static_cast<uint32_t>(0);
  Filter dileptonFilter = aod::reducedpair::pt > fConfigDileptonLowpTCut&& aod::reducedpair::pt<fConfigDileptonHighpTCut && aod::reducedpair::mass> fConfigDileptonLowMass&& aod::reducedpair::mass<fConfigDileptonHighMass && aod::reducedpair::sign == 0 && aod::reducedpair::lxy> fConfigDileptonLxyCut;
  Filter filterBarrel = aod::dqanalysisflags::isBarrelSelected > static_cast<uint32_t>(0);
  Filter filterMuon = aod::dqanalysisflags::isMuonSelected > static_cast<uint32_t>(0);

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesHadron;
  HistogramManager* fHistMan;

  std::vector<MCSignal*> fRecMCSignals;
  std::vector<MCSignal*> fGenMCSignals;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> fHashBin;

  void init(o2::framework::InitContext& context)
  {
    bool isBarrel = context.mOptions.get<bool>("processBarrelSkimmed");
    bool isBarrelAsymmetric = context.mOptions.get<bool>("processDstarToD0Pi");
    bool isMuon = context.mOptions.get<bool>("processMuonSkimmed");
    bool isMCGen = context.mOptions.get<bool>("processMCGen") || context.mOptions.get<bool>("processMCGenWithEventSelection");
    bool isDummy = context.mOptions.get<bool>("processDummy");
    bool isMCGen_energycorrelators = context.mOptions.get<bool>("processMCGenEnergyCorrelators") || context.mOptions.get<bool>("processMCGenEnergyCorrelatorsPion");
    bool isMCGen_energycorrelatorsME = context.mOptions.get<bool>("processMCGenEnergyCorrelatorsME") || context.mOptions.get<bool>("processMCGenEnergyCorrelatorsPionME");

    if (isDummy) {
      if (isBarrel || isMuon || isBarrelAsymmetric || isMCGen) {
        LOG(fatal) << "Dummy function is enabled even if there are normal process functions running! Fix your config!" << endl;
      } else {
        LOG(info) << "Dummy function is enabled. Skipping the rest of the init function" << endl;
        return;
      }
    }

    fCurrentRun = 0;

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      fCCDB->get<TGeoManager>(fConfigGeoPath);
    }

    fValuesDilepton = new float[VarManager::kNVars];
    fValuesHadron = new float[VarManager::kNVars];
    fTrackCutBitMap = 0;
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(true);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    TString sigNamesStr = fConfigMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    if (!sigNamesStr.IsNull()) {
      for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
        MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
        if (sig) {
          if (sig->GetNProngs() != 3) {
            LOG(fatal) << "Signal at reconstructed level requested (" << sig->GetName() << ") " << "does not have 3 prongs! Fix it";
          }
          fRecMCSignals.push_back(sig);
        } else {
          LOG(fatal) << "Signal at reconstructed level requested (" << objRecSigArray->At(isig)->GetName() << ") " << "could not be retrieved from the library! -> skipped";
        }
      }
    }

    // Add the reco MCSignals from the JSON config
    TString addMCSignalsStr = fConfigMCRecSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != 3) {
          LOG(fatal) << "Signal at reconstructed level requested (" << mcIt->GetName() << ") " << "does not have 3 prongs! Fix it";
        }
        fRecMCSignals.push_back(mcIt);
      }
    }

    // Add histogram classes for each specified MCsignal at the generator level
    // TODO: create a std::vector of hist classes to be used at Fill time, to avoid using Form in the process function
    TString sigGenNamesStr = fConfigMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 1) { // NOTE: 1-prong signals required
          fGenMCSignals.push_back(sig);
        }
      }
    }

    // Add the gen MCSignals from the JSON config
    addMCSignalsStr = fConfigMCGenSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() == 1) {
          fGenMCSignals.push_back(mcIt);
        }
      }
    }

    // For each track/muon selection used to produce dileptons, create a separate histogram directory using the
    // name of the track/muon cut.

    // Get the list of single track and muon cuts computed in the dedicated tasks upstream
    // We need this to know the order in which they were computed, and also to make sure that in this task we do not ask
    //   for cuts which were not computed (in which case this will trigger a fatal)
    string cfgTrackSelection_TrackCuts;
    if (isBarrel || isBarrelAsymmetric) {
      getTaskOptionValue<string>(context, "analysis-track-selection", "cfgTrackCuts", cfgTrackSelection_TrackCuts, false);
    } else {
      getTaskOptionValue<string>(context, "analysis-muon-selection", "cfgMuonCuts", cfgTrackSelection_TrackCuts, false);
    }
    TObjArray* cfgTrackSelection_objArrayTrackCuts = nullptr;
    if (!cfgTrackSelection_TrackCuts.empty()) {
      cfgTrackSelection_objArrayTrackCuts = TString(cfgTrackSelection_TrackCuts).Tokenize(",");
    }
    // get also the list of cuts specified via the JSON parameters
    if (isBarrel || isBarrelAsymmetric) {
      getTaskOptionValue<string>(context, "analysis-track-selection", "cfgBarrelTrackCutsJSON", cfgTrackSelection_TrackCuts, false);
    } else {
      getTaskOptionValue<string>(context, "analysis-muon-selection", "cfgMuonCutsJSON", cfgTrackSelection_TrackCuts, false);
    }
    if (!cfgTrackSelection_TrackCuts.empty()) {
      if (cfgTrackSelection_objArrayTrackCuts == nullptr) {
        cfgTrackSelection_objArrayTrackCuts = new TObjArray();
      }
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(cfgTrackSelection_TrackCuts.data());
      for (auto& t : addTrackCuts) {
        TObjString* tempObjStr = new TObjString(t->GetName());
        cfgTrackSelection_objArrayTrackCuts->Add(tempObjStr);
      }
    }
    if (cfgTrackSelection_objArrayTrackCuts->GetEntries() == 0) {
      LOG(fatal) << " No track cuts found in the barrel or muon upstream tasks";
    }
    // store all the computed track cut names in a vector
    for (int icut = 0; icut < cfgTrackSelection_objArrayTrackCuts->GetEntries(); icut++) {
      fTrackCutNames.push_back(cfgTrackSelection_objArrayTrackCuts->At(icut)->GetName());
    }
    // get the list of associated track cuts to be combined with the dileptons,
    //   check that these were computed upstream, and create a bit mask
    TObjArray* cfgDileptonTrack_objArrayTrackCuts = nullptr;
    if (!fConfigTrackCuts.value.empty()) {
      cfgDileptonTrack_objArrayTrackCuts = TString(fConfigTrackCuts.value).Tokenize(",");
    } else {
      LOG(fatal) << " No track cuts specified! Check it out!";
    }
    // loop over these cuts and check they were computed upstream (otherwise trigger a fatal)
    for (int icut = 0; icut < cfgDileptonTrack_objArrayTrackCuts->GetEntries(); icut++) {
      if (!cfgTrackSelection_objArrayTrackCuts->FindObject(cfgDileptonTrack_objArrayTrackCuts->At(icut)->GetName())) {
        LOG(fatal) << "Specified track cut (" << cfgDileptonTrack_objArrayTrackCuts->At(icut)->GetName() << ") not found in the list of computed cuts by the single barrel / muon selection tasks";
      }
    }
    // loop over all the upstream cuts and make a bit mask for the track cuts specified in this task
    for (int icut = 0; icut < cfgTrackSelection_objArrayTrackCuts->GetEntries(); icut++) {
      if (cfgDileptonTrack_objArrayTrackCuts->FindObject(cfgTrackSelection_objArrayTrackCuts->At(icut)->GetName())) {
        fTrackCutBitMap |= (static_cast<uint32_t>(1) << icut);
      }
    }
    // finally, store the total number of upstream tasks, for easy access
    fNCuts = fTrackCutNames.size();

    // get the cuts employed for same-event pairing
    // NOTE: The track/muon cuts in analysis-same-event-pairing are used to select electrons/muons to build dielectrons/dimuons
    // NOTE: The cfgPairCuts in analysis-same-event-pairing are used to apply an additional selection on top of the already produced dileptons
    //        but this is only used for histograms, not for the produced dilepton tables
    string cfgPairing_TrackCuts;
    string cfgPairing_PairCuts;
    string cfgPairing_PairCutsJSON;
    string cfgPairing_CommonTrackCuts;
    if (isBarrel) {
      getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgTrackCuts", cfgPairing_TrackCuts, false);
      getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgPairCuts", cfgPairing_PairCuts, false);
    } else if (isMuon) {
      getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgMuonCuts", cfgPairing_TrackCuts, false);
      getTaskOptionValue<string>(context, "analysis-same-event-pairing", "cfgPairCuts", cfgPairing_PairCuts, false);
    } else if (isBarrelAsymmetric) {
      getTaskOptionValue<string>(context, "analysis-asymmetric-pairing", "cfgLegCuts", cfgPairing_TrackCuts, false);
      getTaskOptionValue<string>(context, "analysis-asymmetric-pairing", "cfgPairCuts", cfgPairing_PairCuts, false);
      getTaskOptionValue<string>(context, "analysis-asymmetric-pairing", "cfgCommonTrackCuts", cfgPairing_CommonTrackCuts, false);
    }
    if (cfgPairing_TrackCuts.empty()) {
      LOG(fatal) << "There are no dilepton cuts specified in the upstream in the same-event-pairing or asymmetric-pairing";
    }

    // If asymmetric pair is used, it may have common track cuts
    TString cfgPairing_strCommonTrackCuts = cfgPairing_CommonTrackCuts;
    if (!cfgPairing_strCommonTrackCuts.IsNull()) { // if common track cuts
      std::unique_ptr<TObjArray> objArrayCommon(cfgPairing_strCommonTrackCuts.Tokenize(","));
      fNCommonTrackCuts = objArrayCommon->GetEntries();
      for (int icut = 0; icut < fNCommonTrackCuts; ++icut) {
        for (int iicut = 0; iicut < cfgTrackSelection_objArrayTrackCuts->GetEntries(); iicut++) {
          if (std::strcmp(cfgTrackSelection_objArrayTrackCuts->At(iicut)->GetName(), objArrayCommon->At(icut)->GetName()) == 0) {
            fCommonTrackCutMap[icut] = iicut;
            fCommonPairCutNames.push_back(objArrayCommon->At(icut)->GetName());
          }
        }
      }
    } // end if (common cuts)

    // Get also the pair cuts specified via the JSON parameters
    if (isBarrelAsymmetric) {
      getTaskOptionValue<string>(context, "analysis-asymmetric-pairing", "cfgPairCutsJSON", cfgPairing_PairCutsJSON, false);
      TString addPairCutsStr = cfgPairing_PairCutsJSON;
      if (addPairCutsStr != "") {
        std::vector<AnalysisCut*> addPairCuts = dqcuts::GetCutsFromJSON(addPairCutsStr.Data());
        for (auto& t : addPairCuts) {
          cfgPairing_PairCuts += Form(",%s", t->GetName());
        }
      }
    }

    std::unique_ptr<TObjArray> objArrayPairCuts(TString(cfgPairing_PairCuts).Tokenize(","));
    fNPairCuts = objArrayPairCuts->GetEntries();
    for (int j = 0; j < fNPairCuts; j++) {
      fPairCutNames.push_back(objArrayPairCuts->At(j)->GetName());
    }

    // array of single lepton cuts specified in the same-analysis-pairing task
    std::unique_ptr<TObjArray> cfgPairing_objArrayTrackCuts(TString(cfgPairing_TrackCuts).Tokenize(","));
    // If asymmetric pairs are used, the number of cuts should come from the asymmetric-pairing task
    if (isBarrelAsymmetric) {
      fNLegCuts = cfgPairing_objArrayTrackCuts->GetEntries();
    } else {
      fNLegCuts = fNCuts;
    }

    // loop over single lepton cuts
    if (isBarrel || isBarrelAsymmetric || isMuon) {
      for (int icut = 0; icut < fNLegCuts; ++icut) {

        TString pairLegCutName;

        // here we check that this cut is one of those used for building the dileptons
        if (!isBarrelAsymmetric) {
          if (!cfgPairing_objArrayTrackCuts->FindObject(fTrackCutNames[icut].Data())) {
            continue;
          }
          pairLegCutName = fTrackCutNames[icut].Data();
        } else {
          // For asymmetric pairs we access the leg cuts instead
          pairLegCutName = static_cast<TObjString*>(cfgPairing_objArrayTrackCuts->At(icut))->GetString();
        }
        fLegCutNames.push_back(pairLegCutName);

        // define dilepton histograms
        DefineHistograms(fHistMan, Form("DileptonsSelected_%s", pairLegCutName.Data()), "barrel,vertexing");
        // loop over track cuts and create dilepton - track histogram directories
        for (int iCutTrack = 0; iCutTrack < fNCuts; iCutTrack++) {

          // here we check that this track cut is one of those required to associate with the dileptons
          if (!(fTrackCutBitMap & (static_cast<uint32_t>(1) << iCutTrack))) {
            continue;
          }

          DefineHistograms(fHistMan, Form("DileptonTrack_%s_%s", pairLegCutName.Data(), fTrackCutNames[iCutTrack].Data()), fConfigHistogramSubgroups.value.data());
          for (auto& sig : fRecMCSignals) {
            DefineHistograms(fHistMan, Form("DileptonTrackMCMatched_%s_%s_%s", pairLegCutName.Data(), fTrackCutNames[iCutTrack].Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());
          }

          if (!cfgPairing_strCommonTrackCuts.IsNull()) {
            std::unique_ptr<TObjArray> objArrayCommon(cfgPairing_strCommonTrackCuts.Tokenize(","));
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              DefineHistograms(fHistMan, Form("DileptonsSelected_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data()), "barrel,vertexing");
              DefineHistograms(fHistMan, Form("DileptonTrack_%s_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data(), fTrackCutNames[iCutTrack].Data()), fConfigHistogramSubgroups.value.data());
              for (auto& sig : fRecMCSignals) {
                DefineHistograms(fHistMan, Form("DileptonTrackMCMatched_%s_%s_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data(), fTrackCutNames[iCutTrack].Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());
              }
            }
          }

          if (fNPairCuts != 0) {

            for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) {
              DefineHistograms(fHistMan, Form("DileptonsSelected_%s_%s", pairLegCutName.Data(), fPairCutNames[iPairCut].Data()), "barrel,vertexing");
              DefineHistograms(fHistMan, Form("DileptonTrack_%s_%s_%s", pairLegCutName.Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iCutTrack].Data()), fConfigHistogramSubgroups.value.data());
              for (auto& sig : fRecMCSignals) {
                DefineHistograms(fHistMan, Form("DileptonTrackMCMatched_%s_%s_%s_%s", pairLegCutName.Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iCutTrack].Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());
              }

              if (!cfgPairing_strCommonTrackCuts.IsNull()) {
                std::unique_ptr<TObjArray> objArrayCommon(cfgPairing_strCommonTrackCuts.Tokenize(","));
                for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                  DefineHistograms(fHistMan, Form("DileptonsSelected_%s_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), "barrel,vertexing");
                  DefineHistograms(fHistMan, Form("DileptonTrack_%s_%s_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iCutTrack].Data()), fConfigHistogramSubgroups.value.data());
                  for (auto& sig : fRecMCSignals) {
                    DefineHistograms(fHistMan, Form("DileptonTrack_%s_%s_%s_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iCutTrack].Data(), sig->GetName()), fConfigHistogramSubgroups.value.data());
                  }
                }
              }
            }
          }
        } // end loop over track cuts to be combined with dileptons / di-tracks
      } // end loop over pair leg track cuts
    } // end if (isBarrel || isBarrelAsymmetric || isMuon)

    if (isMCGen) {
      for (auto& sig : fGenMCSignals) {
        DefineHistograms(fHistMan, Form("MCTruthGen_%s", sig->GetName()), "");
        DefineHistograms(fHistMan, Form("MCTruthGenSel_%s", sig->GetName()), "");
      }
      for (auto& sig : fRecMCSignals) {
        DefineHistograms(fHistMan, Form("MCTruthGenSelBR_%s", sig->GetName()), "");
        DefineHistograms(fHistMan, Form("MCTruthGenSelBRAccepted_%s", sig->GetName()), "");
      }
    }

    if (isMCGen_energycorrelators) {
      for (auto& sig : fGenMCSignals) {
        DefineHistograms(fHistMan, Form("MCTruthEenergyCorrelators_%s", sig->GetName()), "");
      }
    }

    if (isMCGen_energycorrelatorsME) {
      for (auto& sig : fGenMCSignals) {
        DefineHistograms(fHistMan, Form("MCTruthEenergyCorrelatorsME_%s", sig->GetName()), "");
      }
    }

    TString addHistsStr = fConfigAddJSONHistograms.value;
    if (addHistsStr != "") {
      dqhistograms::AddHistogramsFromJSON(fHistMan, addHistsStr.Data());
    }
    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // init parameters from CCDB
  void initParamsFromCCDB(uint64_t timestamp)
  {
    if (fConfigUseRemoteField.value) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigGRPmagPath.value, timestamp);
      float magField = 0.0;
      if (grpmag != nullptr) {
        magField = grpmag->getNominalL3Field();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", timestamp);
      }
      if (fConfigUseKFVertexing.value) {
        VarManager::SetupThreeProngKFParticle(magField);
      } else {
        VarManager::SetupThreeProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false); // TODO: get these parameters from Configurables
      }
    } else {
      if (fConfigUseKFVertexing.value) {
        VarManager::SetupThreeProngKFParticle(fConfigMagField.value);
      } else {
        VarManager::SetupThreeProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false); // TODO: get these parameters from Configurables
      }
    }
  }

  // Template function to run pair - hadron combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks, typename TTrackAssocs, typename TDileptons>
  void runDileptonHadron(TEvent const& event, TTrackAssocs const& assocs, TTracks const& tracks, TDileptons const& dileptons, ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& mcTracks)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesHadron);
    VarManager::FillEvent<VarManager::ObjTypes::ReducedEventMC>(event.reducedMCevent(), fValuesHadron);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);
    VarManager::FillEvent<VarManager::ObjTypes::ReducedEventMC>(event.reducedMCevent(), fValuesDilepton);

    uint32_t mcDecision = static_cast<uint32_t>(0);
    size_t isig = 0;

    for (auto dilepton : dileptons) {
      // get full track info of tracks based on the index
      auto lepton1 = tracks.rawIteratorAt(dilepton.index0Id());
      auto lepton2 = tracks.rawIteratorAt(dilepton.index1Id());
      if (!lepton1.has_reducedMCTrack() || !lepton2.has_reducedMCTrack()) {
        continue;
      }
      auto lepton1MC = lepton1.reducedMCTrack();
      auto lepton2MC = lepton2.reducedMCTrack();
      // Check that the dilepton has zero charge
      if (dilepton.sign() != 0) {
        continue;
      }
      // dilepton rap cut
      float rap = dilepton.rap();
      if (fConfigUseMCRapcut && abs(rap) > fConfigDileptonRapCutAbs)
        continue;

      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton);

      // fill selected dilepton histograms for each specified selection
      for (int icut = 0; icut < fNLegCuts; icut++) {

        if (!dilepton.filterMap_bit(icut)) {
          continue;
        }

        // regular dileptons
        fHistMan->FillHistClass(Form("DileptonsSelected_%s", fLegCutNames[icut].Data()), fValuesDilepton);

        // other pairs, e.g.: D0s
        if constexpr (TCandidateType == VarManager::kDstarToD0KPiPi) { // Dielectrons and Dimuons don't have the PairFilterMap column
          for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
            if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
              fHistMan->FillHistClass(Form("DileptonsSelected_%s_%s", fLegCutNames[icut].Data(), fCommonPairCutNames[iCommonCut].Data()), fValuesDilepton);
            }
          }
          for (int iPairCut = 0; iPairCut < fNPairCuts; iPairCut++) {
            if (dilepton.pairFilterMap_bit(iPairCut)) {
              fHistMan->FillHistClass(Form("DileptonsSelected_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[icut].Data()), fValuesDilepton);
              for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
                if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
                  fHistMan->FillHistClass(Form("DileptonsSelected_%s_%s_%s", fLegCutNames[icut].Data(), fCommonPairCutNames[iCommonCut].Data(), fPairCutNames[icut].Data()), fValuesDilepton);
                }
              }
            }
          }
        }
      }

      // loop over track associations
      for (auto& assoc : assocs) {
        VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
        VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);

        uint32_t trackSelection = 0;
        if constexpr (TCandidateType == VarManager::kBtoJpsiEEK) {
          // check the cuts fulfilled by this candidate track; if none just continue
          trackSelection = (assoc.isBarrelSelected_raw() & fTrackCutBitMap);
          if (!trackSelection) {
            continue;
          }
          // get the track from this association
          auto track = assoc.template reducedtrack_as<TTracks>();
          // check that this track is not included in the current dilepton
          if (track.globalIndex() == dilepton.index0Id() || track.globalIndex() == dilepton.index1Id()) {
            continue;
          }
          // compute needed quantities
          VarManager::FillDileptonHadron(dilepton, track, fValuesHadron);
          VarManager::FillDileptonTrackVertexing<TCandidateType, TEventFillMap, TTrackFillMap>(event, lepton1, lepton2, track, fValuesHadron);

          auto trackMC = track.reducedMCTrack();
          mcDecision = 0;
          isig = 0;
          for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
            if ((*sig)->CheckSignal(true, lepton1MC, lepton2MC, trackMC)) {
              mcDecision |= (static_cast<uint32_t>(1) << isig);
            }
          }

          // fill MC truth values for the B hadron
          auto currentMCParticle = trackMC;
          if (mcDecision > 0) {
            while (true) {
              if (currentMCParticle.has_mothers()) {
                currentMCParticle = currentMCParticle.template mothers_first_as<ReducedMCTracks>();
                if (std::abs(currentMCParticle.pdgCode()) > 500 && std::abs(currentMCParticle.pdgCode()) < 549) { // nb! hardcoded pdgcodes
                  VarManager::FillTrackMC(mcTracks, currentMCParticle, fValuesHadron);
                  break;
                }
              } else {
                break;
              }
            }
            // fill mc truth vertexing (for the associated track as this will have a displaced vertex, while the B hadron is produced in the PV)
            VarManager::FillTrackCollisionMC<VarManager::kBtoJpsiEEK, TEventFillMap>(trackMC, currentMCParticle, event.reducedMCevent(), fValuesHadron);
          }

          // table to be written out for ML analysis
          BmesonsTable(event.runNumber(), event.globalIndex(), event.timestamp(),
                       fValuesHadron[VarManager::kPairMass], dilepton.mass(), fValuesHadron[VarManager::kDeltaMass], fValuesHadron[VarManager::kPairPt], fValuesHadron[VarManager::kPairEta], fValuesHadron[VarManager::kPairPhi], fValuesHadron[VarManager::kPairRap],
                       fValuesHadron[VarManager::kVertexingLxy], fValuesHadron[VarManager::kVertexingLxyErr], fValuesHadron[VarManager::kVertexingLxyz], fValuesHadron[VarManager::kVertexingLxyzErr], fValuesHadron[VarManager::kVertexingLz], fValuesHadron[VarManager::kVertexingLzErr],
                       fValuesHadron[VarManager::kVertexingTauxy], fValuesHadron[VarManager::kVertexingTauxyErr], fValuesHadron[VarManager::kVertexingTauz], fValuesHadron[VarManager::kVertexingTauzErr], fValuesHadron[VarManager::kCosPointingAngle],
                       fValuesHadron[VarManager::kVertexingChi2PCA],
                       fValuesHadron[VarManager::kMCVertexingLxy], fValuesHadron[VarManager::kMCVertexingLxyz], fValuesHadron[VarManager::kMCVertexingLz],
                       fValuesHadron[VarManager::kMCVertexingTauxy], fValuesHadron[VarManager::kMCVertexingTauz], fValuesHadron[VarManager::kMCCosPointingAngle],
                       track.globalIndex(), lepton1.globalIndex(), lepton2.globalIndex(),
                       track.tpcInnerParam(), track.eta(), dilepton.pt(), dilepton.eta(), lepton1.tpcInnerParam(), lepton1.eta(), lepton2.tpcInnerParam(), lepton2.eta(),
                       track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tpcNSigmaPr(), track.tofNSigmaKa(),
                       lepton1.tpcNSigmaEl(), lepton1.tpcNSigmaPi(), lepton1.tpcNSigmaPr(),
                       lepton2.tpcNSigmaEl(), lepton2.tpcNSigmaPi(), lepton2.tpcNSigmaPr(),
                       track.itsClusterMap(), lepton1.itsClusterMap(), lepton2.itsClusterMap(),
                       track.itsChi2NCl(), lepton1.itsChi2NCl(), lepton2.itsChi2NCl(),
                       track.tpcNClsFound(), lepton1.tpcNClsFound(), lepton2.tpcNClsFound(),
                       track.tpcChi2NCl(), lepton1.tpcChi2NCl(), lepton2.tpcChi2NCl(),
                       dilepton.filterMap_raw(), trackSelection, mcDecision);
        }

        if constexpr (TCandidateType == VarManager::kDstarToD0KPiPi) {
          trackSelection = (assoc.isBarrelSelected_raw() & fTrackCutBitMap);
          if (!trackSelection) {
            continue;
          }

          auto track = assoc.template reducedtrack_as<TTracks>();
          if (track.globalIndex() == dilepton.index0Id() || track.globalIndex() == dilepton.index1Id()) {
            continue;
          }
          // Check that the charge combination makes sense for D*+ -> D0 pi+ or D*- -> D0bar pi-
          if (!((track.sign() == 1 && lepton1.sign() == -1 && lepton2.sign() == 1) || (track.sign() == -1 && lepton1.sign() == 1 && lepton2.sign() == -1))) {
            continue;
          }
          VarManager::FillDileptonHadron(dilepton, track, fValuesHadron);
          VarManager::FillDileptonTrackVertexing<TCandidateType, TEventFillMap, TTrackFillMap>(event, lepton1, lepton2, track, fValuesHadron);

          auto trackMC = track.reducedMCTrack();
          mcDecision = 0;
          isig = 0;
          for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
            if ((*sig)->CheckSignal(true, lepton1MC, lepton2MC, trackMC)) {
              mcDecision |= (static_cast<uint32_t>(1) << isig);
            }
          }
        }

        if constexpr (TCandidateType == VarManager::kBcToThreeMuons) {
          trackSelection = (assoc.isMuonSelected_raw() & fTrackCutBitMap);
          if (!trackSelection) {
            continue;
          }

          auto track = assoc.template reducedmuon_as<TTracks>();
          if (track.globalIndex() == dilepton.index0Id() || track.globalIndex() == dilepton.index1Id()) {
            continue;
          }

          VarManager::FillDileptonHadron(dilepton, track, fValuesHadron);
          VarManager::FillDileptonTrackVertexing<TCandidateType, TEventFillMap, TTrackFillMap>(event, lepton1, lepton2, track, fValuesHadron);

          if (!track.has_reducedMCTrack()) {
            continue;
          }
          auto trackMC = track.reducedMCTrack();
          mcDecision = 0;
          isig = 0;
          for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
            if ((*sig)->CheckSignal(true, lepton1MC, lepton2MC, trackMC)) {
              mcDecision |= (static_cast<uint32_t>(1) << isig);
            }
          }
          // Fill table for correlation analysis
          DileptonTrackTable(fValuesHadron[VarManager::kDeltaEta], fValuesHadron[VarManager::kDeltaPhi],
                             dilepton.mass(), dilepton.pt(), dilepton.eta(), track.pt(), track.eta(), track.phi(),
                             lepton1.pt(), lepton1.eta(), lepton1.phi(), lepton2.pt(), lepton2.eta(), lepton2.phi(),
                             mcDecision);
        }

        // Fill histograms for the triplets
        // loop over dilepton / ditrack cuts and MC signals
        for (int icut = 0; icut < fNCuts; icut++) {

          if (!dilepton.filterMap_bit(icut)) {
            continue;
          }

          // loop over specified track cuts (the tracks to be combined with the dileptons)
          for (int iTrackCut = 0; iTrackCut < fNCuts; iTrackCut++) {

            if (!(trackSelection & (static_cast<uint32_t>(1) << iTrackCut))) {
              continue;
            }

            fHistMan->FillHistClass(Form("DileptonTrack_%s_%s", fLegCutNames[icut].Data(), fTrackCutNames[iTrackCut].Data()), fValuesHadron);
            for (uint32_t isig = 0; isig < fRecMCSignals.size(); isig++) {
              if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                fHistMan->FillHistClass(Form("DileptonTrackMCMatched_%s_%s_%s", fLegCutNames[icut].Data(), fTrackCutNames[iTrackCut].Data(), fRecMCSignals[isig]->GetName()), fValuesHadron);
              }
            }

            if constexpr (TCandidateType == VarManager::kDstarToD0KPiPi) { // Dielectrons and Dimuons don't have the PairFilterMap column
              for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
                if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
                  fHistMan->FillHistClass(Form("DileptonTrack_%s_%s_%s", fLegCutNames[icut].Data(), fCommonPairCutNames[iCommonCut].Data(), fTrackCutNames[iTrackCut].Data()), fValuesHadron);
                  for (uint32_t isig = 0; isig < fRecMCSignals.size(); isig++) {
                    if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                      fHistMan->FillHistClass(Form("DileptonTrackMCMatched_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonPairCutNames[iCommonCut].Data(), fTrackCutNames[iTrackCut].Data(), fRecMCSignals[isig]->GetName()), fValuesHadron);
                    }
                  }
                }
              }
              for (int iPairCut = 0; iPairCut < fNPairCuts; iPairCut++) {
                if (dilepton.pairFilterMap_bit(iPairCut)) {
                  fHistMan->FillHistClass(Form("DileptonTrack_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iTrackCut].Data()), fValuesHadron);
                  for (uint32_t isig = 0; isig < fRecMCSignals.size(); isig++) {
                    if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                      fHistMan->FillHistClass(Form("DileptonTrackMCMatched_%s_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iTrackCut].Data(), fRecMCSignals[isig]->GetName()), fValuesHadron);
                    }
                  }
                  for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
                    if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
                      fHistMan->FillHistClass(Form("DileptonTrack_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonPairCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iTrackCut].Data()), fValuesHadron);
                      for (uint32_t isig = 0; isig < fRecMCSignals.size(); isig++) {
                        if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
                          fHistMan->FillHistClass(Form("DileptonTrackMCMatched_%s_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonPairCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iTrackCut].Data(), fRecMCSignals[isig]->GetName()), fValuesHadron);
                        }
                      }
                    }
                  }
                }
              }
            }
          } // end loop over track cuts
        } // end loop over dilepton cuts
      } // end loop over associations
    } // end loop over dileptons
  }

  Preslice<aod::ReducedTracksAssoc> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<MyDielectronCandidates> dielectronsPerCollision = aod::reducedpair::reducedeventId;
  Preslice<MyDitrackCandidates> ditracksPerCollision = aod::reducedpair::reducedeventId;

  void processBarrelSkimmed(soa::Filtered<MyEventsVtxCovSelected> const& events,
                            soa::Filtered<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> const& assocs,
                            MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons,
                            ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    // set up KF or DCAfitter
    if (events.size() == 0) {
      return;
    }
    if (fCurrentRun != events.begin().runNumber()) { // start: runNumber
      initParamsFromCCDB(events.begin().timestamp());
      fCurrentRun = events.begin().runNumber();
    } // end: runNumber
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      auto groupedBarrelAssocs = assocs.sliceBy(trackAssocsPerCollision, event.globalIndex());
      auto groupedDielectrons = dileptons.sliceBy(dielectronsPerCollision, event.globalIndex());
      runDileptonHadron<VarManager::kBtoJpsiEEK, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, groupedBarrelAssocs, tracks, groupedDielectrons, mcEvents, mcTracks);
    }
  }

  void processDstarToD0Pi(soa::Filtered<MyEventsVtxCovSelected> const& events,
                          soa::Filtered<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> const& assocs,
                          MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDitrackCandidates> const& ditracks,
                          ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    // set up KF or DCAfitter
    if (events.size() == 0) {
      return;
    }
    if (fCurrentRun != events.begin().runNumber()) { // start: runNumber
      initParamsFromCCDB(events.begin().timestamp());
      fCurrentRun = events.begin().runNumber();
    } // end: runNumber
    for (auto& event : events) {
      auto groupedBarrelAssocs = assocs.sliceBy(trackAssocsPerCollision, event.globalIndex());
      auto groupedDitracks = ditracks.sliceBy(ditracksPerCollision, event.globalIndex());
      runDileptonHadron<VarManager::kDstarToD0KPiPi, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, groupedBarrelAssocs, tracks, groupedDitracks, mcEvents, mcTracks);
    }
  }

  Preslice<aod::ReducedMuonsAssoc> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<MyDimuonCandidates> dimuonsPerCollision = aod::reducedpair::reducedeventId;

  void processMuonSkimmed(soa::Filtered<MyEventsVtxCovSelected> const& events,
                          soa::Filtered<soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>> const& assocs,
                          MyMuonTracksWithCov const& tracks, soa::Filtered<MyDimuonCandidates> const& dileptons,
                          ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    // set up KF or DCAfitter
    if (events.size() == 0) {
      return;
    }
    if (fCurrentRun != events.begin().runNumber()) { // start: runNumber
      initParamsFromCCDB(events.begin().timestamp());
      fCurrentRun = events.begin().runNumber();
    } // end: runNumber
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      auto groupedMuonAssocs = assocs.sliceBy(muonAssocsPerCollision, event.globalIndex());
      auto groupedDimuons = dileptons.sliceBy(dimuonsPerCollision, event.globalIndex());
      runDileptonHadron<VarManager::kBcToThreeMuons, gkEventFillMapWithCov, gkMuonFillMapWithCov>(event, groupedMuonAssocs, tracks, groupedDimuons, mcEvents, mcTracks);
    }
  }

  void processMCGen(ReducedMCTracks const& mcTracks)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // group all the MC tracks which belong to the MC event corresponding to the current reconstructed event
    // auto groupedMCTracks = tracksMC.sliceBy(aod::reducedtrackMC::reducedMCeventId, event.reducedMCevent().globalIndex());
    for (auto& mctrack : mcTracks) {

      if ((std::abs(mctrack.pdgCode()) > 400 && std::abs(mctrack.pdgCode()) < 599) ||
          (std::abs(mctrack.pdgCode()) > 4000 && std::abs(mctrack.pdgCode()) < 5999) ||
          mctrack.mcReducedFlags() > 0) {
        /*cout << ">>>>>>>>>>>>>>>>>>>>>>> track idx / pdg / selections: " << mctrack.globalIndex() << " / " << mctrack.pdgCode() << " / ";
        PrintBitMap(mctrack.mcReducedFlags(), 16);
        cout << endl;
        if (mctrack.has_mothers()) {
          for (auto& m : mctrack.mothersIds()) {
            if (m < mcTracks.size()) { // protect against bad mother indices
              auto aMother = mcTracks.rawIteratorAt(m);
              cout << "<<<<<< mother idx / pdg: " << m << " / " << aMother.pdgCode() << endl;
            }
          }
        }

        if (mctrack.has_daughters()) {
          for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
            if (d < mcTracks.size()) { // protect against bad daughter indices
              auto aDaughter = mcTracks.rawIteratorAt(d);
              cout << "<<<<<< daughter idx / pdg: " << d << " / " << aDaughter.pdgCode() << endl;
            }
          }
        }*/
      }

      VarManager::FillTrackMC(mcTracks, mctrack);
      // NOTE: Signals are checked here mostly based on the skimmed MC stack, so depending on the requested signal, the stack could be incomplete.
      // NOTE: However, the working model is that the decisions on MC signals are precomputed during skimming and are stored in the mcReducedFlags member.
      // TODO:  Use the mcReducedFlags to select signals
      for (auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, mctrack)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }
  }

  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;

  void processMCGenWithEventSelection(soa::Filtered<MyEventsVtxCovSelected> const& events,
                                      ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& mcTracks)
  {
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }

      auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reducedMCeventId());
      groupedMCTracks.bindInternalIndicesTo(&mcTracks);
      for (auto& track : groupedMCTracks) {

        VarManager::FillTrackMC(mcTracks, track);

        auto track_raw = groupedMCTracks.rawIteratorAt(track.globalIndex());
        for (auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, track_raw)) {
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
          }
        }
      }

      // make a list of all MC tracks in the MC collision corresponding to the current reconstructed event
      std::vector<size_t> mcTrackIndices;
      for (auto& t : groupedMCTracks) {
        mcTrackIndices.push_back(t.globalIndex());
      }

      // make a three nested for loop over all MC tracks in the vector
      for (auto t1 : mcTrackIndices) {
        auto track1 = mcTracks.rawIteratorAt(*(&t1));
        for (auto t2 : mcTrackIndices) {
          if (t1 == t2)
            continue;
          // if (t2 < t1) continue;
          auto track2 = mcTracks.rawIteratorAt(*(&t2));
          for (auto t3 : mcTrackIndices) {
            if (t3 == t1)
              continue;
            if (t3 == t2)
              continue;
            auto track3 = mcTracks.rawIteratorAt(*(&t3));

            for (auto& sig : fRecMCSignals) {
              if (sig->CheckSignal(true, track1, track2, track3)) {

                VarManager::FillTripleMC<VarManager::kBtoJpsiEEK>(track1, track2, track3, VarManager::fgValues); // nb! hardcoded for jpsiK

                fHistMan->FillHistClass(Form("MCTruthGenSelBR_%s", sig->GetName()), VarManager::fgValues);

                // apply kinematic cuts
                if (track1.pt() < fConfigMCGenDileptonLegPtMin.value || std::abs(track1.eta()) > fConfigMCGenDileptonLegEtaAbs.value) {
                  continue;
                }
                if (track2.pt() < fConfigMCGenDileptonLegPtMin.value || std::abs(track2.eta()) > fConfigMCGenDileptonLegEtaAbs.value) {
                  continue;
                }
                if (track3.pt() < fConfigMCGenHadronPtMin.value || std::abs(track3.eta()) > fConfigMCGenHadronEtaAbs.value) {
                  continue;
                }
                fHistMan->FillHistClass(Form("MCTruthGenSelBRAccepted_%s", sig->GetName()), VarManager::fgValues);
              }
            }
          }
        }
      }
    } // end loop over reconstructed events
  }

  template <int THadronMassType, typename TEvent, typename TMCTracks>
  void runEnergyCorrelators(TEvent const& event, TMCTracks const& mcTracks)
  {
    auto groupedMCTracks = mcTracks.sliceBy(perReducedMcEvent, event.reducedMCeventId());
    groupedMCTracks.bindInternalIndicesTo(&mcTracks);
    for (auto& t1 : groupedMCTracks) {
      auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
      // apply kinematic cuts for signal
      if ((t1_raw.pt() < fConfigDileptonLowpTCut || t1_raw.pt() > fConfigDileptonHighpTCut))
        continue;
      if (abs(t1_raw.y()) > fConfigDileptonRapCutAbs)
        continue;
      // for the energy correlators
      for (auto& t2 : groupedMCTracks) {
        auto t2_raw = groupedMCTracks.rawIteratorAt(t2.globalIndex());
        if (TMath::Abs(t2_raw.pdgCode()) == 443 || TMath::Abs(t2_raw.pdgCode()) == 11 || TMath::Abs(t2_raw.pdgCode()) == 22)
          continue;
        if (t2_raw.pt() < fConfigMCGenHadronPtMin.value || std::abs(t2_raw.eta()) > fConfigMCGenHadronEtaAbs.value)
          continue;
        if (t2_raw.getGenStatusCode() <= 0)
          continue;
        VarManager::FillEnergyCorrelatorsMC<THadronMassType>(t1_raw, t2_raw, VarManager::fgValues);
        for (auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, t1_raw)) {
            fHistMan->FillHistClass(Form("MCTruthEenergyCorrelators_%s", sig->GetName()), VarManager::fgValues);
          }
        }
      }
    }
  }
  void processMCGenEnergyCorrelators(soa::Filtered<MyEventsVtxCovSelected> const& events,
                                     ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& mcTracks)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }
      runEnergyCorrelators<VarManager::kJpsiHadronMass>(event, mcTracks);
    }
  }

  void processMCGenEnergyCorrelatorsPion(soa::Filtered<MyEventsVtxCovSelected> const& events,
                                         ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& mcTracks)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (!event.has_reducedMCevent()) {
        continue;
      }
      runEnergyCorrelators<VarManager::kJpsiPionMass>(event, mcTracks);
    }
  }

  template <int THadronMassType, typename TEvent, typename TMCTracks>
  void runEnergyCorrelatorsMixedEvent(TEvent const& event1, TEvent const& event2, TMCTracks const& mcTracks)
  {
    auto groupedMCTracks1 = mcTracks.sliceBy(perReducedMcEvent, event1.reducedMCeventId());
    auto groupedMCTracks2 = mcTracks.sliceBy(perReducedMcEvent, event2.reducedMCeventId());
    groupedMCTracks1.bindInternalIndicesTo(&mcTracks);
    groupedMCTracks2.bindInternalIndicesTo(&mcTracks);
    for (auto& t1 : groupedMCTracks1) {
      auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
      // apply kinematic cuts for signal
      if ((t1_raw.pt() < fConfigDileptonLowpTCut || t1_raw.pt() > fConfigDileptonHighpTCut))
        continue;
      if (abs(t1_raw.y()) > fConfigDileptonRapCutAbs)
        continue;
      // for the energy correlators
      for (auto& t2 : groupedMCTracks2) {
        auto t2_raw = groupedMCTracks2.rawIteratorAt(t2.globalIndex());
        if (TMath::Abs(t2_raw.pdgCode()) == 443 || TMath::Abs(t2_raw.pdgCode()) == 11 || TMath::Abs(t2_raw.pdgCode()) == 22)
          continue;
        if (t2_raw.pt() < fConfigMCGenHadronPtMin.value || std::abs(t2_raw.eta()) > fConfigMCGenHadronEtaAbs.value)
          continue;
        if (t2_raw.getGenStatusCode() <= 0)
          continue;
        VarManager::FillEnergyCorrelatorsMC<THadronMassType>(t1_raw, t2_raw, VarManager::fgValues);
        for (auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, t1_raw)) {
            fHistMan->FillHistClass(Form("MCTruthEenergyCorrelatorsME_%s", sig->GetName()), VarManager::fgValues);
          }
        }
      }
    }
  }

  void processMCGenEnergyCorrelatorsME(soa::Filtered<MyEventsVtxCovHashSelected> const& events,
                                       ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& mcTracks)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    // loop over two event comibnations
    for (auto& [event1, event2] : selfCombinations(fHashBin, fConfigMixingDepth.value, -1, events, events)) {
      if (!event1.isEventSelected_bit(0) || !event2.isEventSelected_bit(0)) {
        continue;
      }
      if (!event1.has_reducedMCevent() || !event2.has_reducedMCevent()) {
        continue;
      }
      runEnergyCorrelatorsMixedEvent<VarManager::kJpsiHadronMass>(event1, event2, mcTracks);
    }
  }

  void processMCGenEnergyCorrelatorsPionME(soa::Filtered<MyEventsVtxCovHashSelected> const& events,
                                           ReducedMCEvents const& /*mcEvents*/, ReducedMCTracks const& mcTracks)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    // loop over two event comibnations
    for (auto& [event1, event2] : selfCombinations(fHashBin, fConfigMixingDepth.value, -1, events, events)) {
      if (!event1.isEventSelected_bit(0) || !event2.isEventSelected_bit(0)) {
        continue;
      }
      if (!event1.has_reducedMCevent() || !event2.has_reducedMCevent()) {
        continue;
      }
      runEnergyCorrelatorsMixedEvent<VarManager::kJpsiPionMass>(event1, event2, mcTracks);
    }
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisDileptonTrack, processBarrelSkimmed, "Run barrel dilepton-track pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processDstarToD0Pi, "Run barrel pairing of D0 daughters with pion candidate, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMuonSkimmed, "Run muon dilepton-track pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMCGen, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMCGenWithEventSelection, "Loop over MC particle stack and fill generator level histograms", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMCGenEnergyCorrelators, "Loop over MC particle stack and fill generator level histograms(energy correlators)", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMCGenEnergyCorrelatorsPion, "Loop over MC particle stack and fill generator level histograms(energy correlators)", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMCGenEnergyCorrelatorsME, "Loop over MC particle stack and fill generator level histograms(energy correlators)", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMCGenEnergyCorrelatorsPionME, "Loop over MC particle stack and fill generator level histograms(energy correlators)", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEventSelection>(cfgc),
    adaptAnalysisTask<AnalysisTrackSelection>(cfgc),
    adaptAnalysisTask<AnalysisMuonSelection>(cfgc),
    adaptAnalysisTask<AnalysisPrefilterSelection>(cfgc),
    adaptAnalysisTask<AnalysisSameEventPairing>(cfgc),
    adaptAnalysisTask<AnalysisAsymmetricPairing>(cfgc),
    adaptAnalysisTask<AnalysisDileptonTrack>(cfgc)};
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

    if (classStr.Contains("MCTruthGenPair")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_pair", histName);
    }

    if (classStr.Contains("MCTruthGenSelBR")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_triple");
    } else if (classStr.Contains("MCTruthGen")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_track");
    }

    // if (classStr.Contains("MCTruthGen")) {
    //   dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_track");
    // }

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
