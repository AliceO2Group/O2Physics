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
#include "PWGDQ/Core/DQMlResponse.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/CCDB/EventSelectionParams.h"
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
#include "Framework/Configurable.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"
#include "ITSMFTBase/DPLAlpideParam.h"

#include "TGeoGlobalMagField.h"
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TObjString.h>
#include <TString.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
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
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 8);                      //! Event decision
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
DECLARE_SOA_COLUMN(deltamassBcandidate, deltaMBcandidate, float);
DECLARE_SOA_COLUMN(pTBcandidate, PtBcandidate, float);
DECLARE_SOA_COLUMN(EtaBcandidate, etaBcandidate, float);
DECLARE_SOA_COLUMN(PhiBcandidate, phiBcandidate, float);
DECLARE_SOA_COLUMN(RapBcandidate, rapBcandidate, float);
DECLARE_SOA_COLUMN(LxyBcandidate, lxyBcandidate, float);
DECLARE_SOA_COLUMN(LxyzBcandidate, lxyzBcandidate, float);
DECLARE_SOA_COLUMN(LzBcandidate, lzBcandidate, float);
DECLARE_SOA_COLUMN(TauxyBcandidate, tauxyBcandidate, float);
DECLARE_SOA_COLUMN(TauzBcandidate, tauzBcandidate, float);
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
// Candidate columns for prompt-non-prompt JPsi separation
DECLARE_SOA_COLUMN(Massee, massJPsi2ee, float);
DECLARE_SOA_COLUMN(Ptee, ptJPsi2ee, float);
DECLARE_SOA_COLUMN(Lxyee, lxyJPsi2ee, float);
DECLARE_SOA_COLUMN(LxyeePoleMass, lxyJPsi2eePoleMass, float);
DECLARE_SOA_COLUMN(Lzee, lzJPsi2ee, float);
DECLARE_SOA_COLUMN(AmbiguousInBunchPairs, AmbiguousJpsiPairsInBunch, bool);
DECLARE_SOA_COLUMN(AmbiguousOutOfBunchPairs, AmbiguousJpsiPairsOutOfBunch, bool);
DECLARE_SOA_COLUMN(MultiplicityFT0A, multiplicityFT0AJPsi2ee, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0CJPsi2ee, float);
DECLARE_SOA_COLUMN(PercentileFT0M, percentileFT0MJPsi2ee, float);
DECLARE_SOA_COLUMN(MultiplicityNContrib, multiplicityNContribJPsi2ee, float);
// Candidate columns for JPsi/muon correlations
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTSA", dqanalysisflags::IsEventSelected);                                                            //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASHA", dqanalysisflags::MixingHash);                                                             //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTSA", dqanalysisflags::IsBarrelSelected);                                                    //!  joinable to ReducedTracksAssoc
DECLARE_SOA_TABLE(BarrelAmbiguities, "AOD", "DQBARRELAMBA", dqanalysisflags::BarrelAmbiguityInBunch, dqanalysisflags::BarrelAmbiguityOutOfBunch); //!  joinable to ReducedBarrelTracks
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTSA", dqanalysisflags::IsMuonSelected);                                                       //!  joinable to ReducedMuonsAssoc
DECLARE_SOA_TABLE(MuonAmbiguities, "AOD", "DQMUONAMBA", dqanalysisflags::MuonAmbiguityInBunch, dqanalysisflags::MuonAmbiguityOutOfBunch);         //!  joinable to ReducedMuonTracks
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTERA", dqanalysisflags::IsBarrelSelectedPrefilter);                                                  //!  joinable to ReducedTracksAssoc
DECLARE_SOA_TABLE(BmesonCandidates, "AOD", "DQBMESONSA",
                  dqanalysisflags::RunNumber, dqanalysisflags::EventIdx, dqanalysisflags::EventTimestamp,
                  dqanalysisflags::massBcandidate, dqanalysisflags::MassDileptonCandidate, dqanalysisflags::deltamassBcandidate, dqanalysisflags::pTBcandidate, dqanalysisflags::EtaBcandidate, dqanalysisflags::PhiBcandidate, dqanalysisflags::RapBcandidate,
                  dqanalysisflags::LxyBcandidate, dqanalysisflags::LxyzBcandidate, dqanalysisflags::LzBcandidate,
                  dqanalysisflags::TauxyBcandidate, dqanalysisflags::TauzBcandidate, dqanalysisflags::CosPBcandidate, dqanalysisflags::Chi2Bcandidate,
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
DECLARE_SOA_TABLE(JPsiMuonCandidates, "AOD", "DQJPSIMUONA",
                  dqanalysisflags::DeltaEta, dqanalysisflags::DeltaPhi,
                  dqanalysisflags::MassDileptonCandidate, dqanalysisflags::Ptpair, dqanalysisflags::Etapair, dqanalysisflags::Ptassoc, dqanalysisflags::Etaassoc, dqanalysisflags::Phiassoc,
                  dqanalysisflags::Ptleg1, dqanalysisflags::Etaleg1, dqanalysisflags::Phileg1, dqanalysisflags::Ptleg2, dqanalysisflags::Etaleg2, dqanalysisflags::Phileg2);
DECLARE_SOA_TABLE(JPsieeCandidates, "AOD", "DQPSEUDOPROPER", dqanalysisflags::Massee, dqanalysisflags::Ptee, dqanalysisflags::Lxyee, dqanalysisflags::LxyeePoleMass, dqanalysisflags::Lzee, dqanalysisflags::AmbiguousInBunchPairs, dqanalysisflags::AmbiguousOutOfBunchPairs, dqanalysisflags::MultiplicityFT0A, dqanalysisflags::MultiplicityFT0C, dqanalysisflags::PercentileFT0M, dqanalysisflags::MultiplicityNContrib);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsMultExtraQVector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll, aod::ReducedEventsQvectorCentr, aod::ReducedEventsQvectorCentrExtra>;
using MyEventsZdc = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedZdcs>;
using MyEventsMultExtraZdc = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll, aod::ReducedZdcs>;
using MyEventsMultExtraZdcFit = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll, aod::ReducedZdcs, aod::ReducedFITs>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsMultExtraSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll, aod::EventCuts>;
using MyEventsVtxCovSelectedMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;
using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvectorCentr, aod::ReducedEventsQvectorCentrExtra>;
using MyEventsVtxCovSelectedQvectorWithHash = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvectorCentr, aod::ReducedEventsQvectorCentrExtra, aod::MixingHashes>;
using MyEventsVtxCovZdcFitSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::ReducedZdcs, aod::ReducedFITs, aod::EventCuts>;
using MyEventsVtxCovZdcFitSelectedMultExtra = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::ReducedZdcs, aod::ReducedFITs, aod::EventCuts, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
using MyEventsHashSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector>;
using MyEventsQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvectorCentr, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;
using MyEventsQvectorCentrSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvectorCentr, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll, aod::EventCuts>;
using MyEventsHashSelectedQvectorCentr = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvectorCentr, aod::ReducedEventsMultPV, aod::ReducedEventsMultAll>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithAmbiguities = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithCovWithAmbiguities = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities>;
using MyBarrelTracksWithCovWithAmbiguitiesWithColl = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelAmbiguities, aod::ReducedTracksBarrelInfo>;
using MyDielectronCandidates = soa::Join<aod::Dielectrons, aod::DielectronsExtra>;
using MyDitrackCandidates = soa::Join<aod::Ditracks, aod::DitracksExtra>;
using MyDimuonCandidates = soa::Join<aod::Dimuons, aod::DimuonsExtra>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksWithCovWithAmbiguities = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonAmbiguities>;
using MyMuonTracksSelectedWithColl = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkEventFillMapWithZdc = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ReducedZdc;
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
// constexpr static uint32_t gkEventFillMapWithCovFlow = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkEventFillMapWithCovZdcFit = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ReducedZdc | VarManager::ReducedFit;
constexpr static uint32_t gkEventFillMapWithMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventMultExtra;
// New fillmap
constexpr static uint32_t gkEventFillMapWithMultExtraWithQVector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventMultExtra | VarManager::ObjTypes::CollisionQvect;
constexpr static uint32_t gkEventFillMapWithMultExtraZdc = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventMultExtra | VarManager::ReducedZdc;
constexpr static uint32_t gkEventFillMapWithMultExtraZdcFit = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventMultExtra | VarManager::ReducedZdc | VarManager::ReducedFit;
constexpr static uint32_t gkEventFillMapWithCovZdcFitMultExtra = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ReducedZdc | VarManager::ReducedFit | VarManager::ReducedEventMultExtra;
constexpr static uint32_t gkEventFillMapWithQvectorCentr = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::CollisionQvect | VarManager::ObjTypes::ReducedEventMultExtra;
// constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
// constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkTrackFillMapWithCovWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;

// constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;
// constexpr static uint32_t gkMuonFillMapWithColl = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCollInfo;

constexpr static int pairTypeEE = VarManager::kDecayToEE;
constexpr static int pairTypeMuMu = VarManager::kDecayToMuMu;
// constexpr static int pairTypeEMu = VarManager::kElectronMuon;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups); // defines histograms for all tasks

template <typename TMap>
void PrintBitMap(TMap map, int nbits)
{
  for (int i = 0; i < nbits; i++) {
    cout << ((map & (TMap(1) << i)) > 0 ? "1" : "0");
  }
}

// Analysis task that produces event decisions (analysis cut, in bunch pileup and split collision check) and the Hash table used in event mixing
struct AnalysisEventSelection {

  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};

  // TODO: Provide the mixing variables and binning directly via configurables (e.g. vectors of float)
  Configurable<std::string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigEventCutsJSON{"cfgEventCutsJSON", "", "Additional event cuts specified in JSON format"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Add event histograms defined via JSON formatting (see HistogramsLibrary)"};
  Configurable<bool> fConfigQA{"cfgQA", true, "If true, QA histograms will be created and filled"};

  Configurable<int> fConfigITSROFrameStartBorderMargin{"cfgITSROFrameStartBorderMargin", -1, "Number of bcs at the start of ITS RO Frame border. Take from CCDB if -1"};
  Configurable<int> fConfigITSROFrameEndBorderMargin{"cfgITSROFrameEndBorderMargin", -1, "Number of bcs at the end of ITS RO Frame border. Take from CCDB if -1"};

  Configurable<float> fConfigSplitCollisionsDeltaZ{"cfgSplitCollisionsDeltaZ", 1.0, "maximum delta-z (cm) between two collisions to consider them as split candidates"};
  Configurable<unsigned int> fConfigSplitCollisionsDeltaBC{"cfgSplitCollisionsDeltaBC", 100, "maximum delta-BC between two collisions to consider them as split candidates; do not apply if value is negative"};
  Configurable<bool> fConfigCheckSplitCollisions{"cfgCheckSplitCollisions", false, "If true, run the split collision check and fill histograms"};

  Configurable<bool> fConfigRunZorro{"cfgRunZorro", false, "Enable event selection with zorro [WARNING: under debug, do not enable!]"};
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

    bool isAnyProcessEnabled = context.mOptions.get<bool>("processSkimmed") || context.mOptions.get<bool>("processSkimmedWithZdc") || context.mOptions.get<bool>("processSkimmedWithMultExtra") || context.mOptions.get<bool>("processSkimmedWithMultExtraZdc") || context.mOptions.get<bool>("processSkimmedWithMultExtraZdcFit") || context.mOptions.get<bool>("processSkimmedWithQvectorCentr");
    bool isDummyEnabled = context.mOptions.get<bool>("processDummy");

    if (isDummyEnabled) {
      if (isAnyProcessEnabled) {
        LOG(fatal) << "You have enabled both the dummy process and at least one normal process function! Check your config!";
      }
      return;
    }
    if (!isAnyProcessEnabled) {
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
        DefineHistograms(fHistMan, "SameBunchCorrelations;OutOfBunchCorrelations;", "");
      }
      // Additional histograms via JSON
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
  }

  template <uint32_t TEventFillMap, typename TEvents>
  void runEventSelection(TEvents const& events)
  {
    if (events.size() > 0 && events.begin().runNumber() != fCurrentRun) {
      std::map<std::string, std::string> metadataRCT, header;
      header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", events.begin().runNumber()), metadataRCT, -1);
      uint64_t sor = std::atol(header["SOR"].c_str());
      uint64_t eor = std::atol(header["EOR"].c_str());
      VarManager::SetSORandEOR(sor, eor);

      auto alppar = fCCDB->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", events.begin().timestamp());
      EventSelectionParams* par = fCCDB->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", events.begin().timestamp());
      int itsROFrameStartBorderMargin = fConfigITSROFrameStartBorderMargin < 0 ? par->fITSROFrameStartBorderMargin : fConfigITSROFrameStartBorderMargin;
      int itsROFrameEndBorderMargin = fConfigITSROFrameEndBorderMargin < 0 ? par->fITSROFrameEndBorderMargin : fConfigITSROFrameEndBorderMargin;
      VarManager::SetITSROFBorderselection(alppar->roFrameBiasInBC, alppar->roFrameLengthInBC, itsROFrameStartBorderMargin, itsROFrameEndBorderMargin);
      fCurrentRun = events.begin().runNumber();
    }

    fSelMap.clear();
    fBCCollMap.clear();

    for (auto& event : events) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<TEventFillMap>(event);

      bool decision = false;
      // Fill histograms in the class Event, before cuts
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }

      // Apply event cuts and fill histograms after selection
      if (fEventCut->IsSelected(VarManager::fgValues)) {
        if (fConfigRunZorro) {
          if (event.tag_bit(56)) { // This is the bit used for the software trigger event selections [TO BE DONE: find a more clear way to use it]
            if (fConfigQA) {
              fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
            }
            decision = true;
          }
        } else {
          if (fConfigQA) {
            fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
          }
          decision = true;
        }
      }

      // fill the event decision map
      fSelMap[event.globalIndex()] = decision;

      // Fill the BC map of events
      if (fBCCollMap.find(event.globalBC()) == fBCCollMap.end()) {
        std::vector<int64_t> evIndices = {event.globalIndex()};
        fBCCollMap[event.globalBC()] = evIndices;
      } else {
        auto& evIndices = fBCCollMap[event.globalBC()];
        evIndices.push_back(event.globalIndex());
      }

      // create the mixing hash and publish it into the hash table
      if (fMixHandler != nullptr) {
        int hh = fMixHandler->FindEventCategory(VarManager::fgValues);
        hash(hh);
      }
    }
  }

  template <uint32_t TEventFillMap, typename TEvents>
  void publishSelections(TEvents const& events)
  {
    // Create a map for collisions which are candidate of being split
    // key: event global index, value: whether pileup event is a possible splitting
    std::map<int64_t, bool> collisionSplittingMap;

    if (fConfigCheckSplitCollisions) {
      // Reset the fValues array and fill event observables
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      // loop over the BC map, get the collision vectors and make in-bunch and out of bunch 2-event correlations
      for (auto bc1It = fBCCollMap.begin(); bc1It != fBCCollMap.end(); ++bc1It) {
        uint64_t bc1 = bc1It->first;
        auto& bc1Events = bc1It->second;

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
          auto& bc2Events = bc2It->second;

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
    }

    // publish the table
    uint8_t evSel = static_cast<uint8_t>(0);
    for (auto& event : events) {
      evSel = 0;
      if (fSelMap[event.globalIndex()]) { // event passed the user cuts
        evSel |= (static_cast<uint8_t>(1) << 0);
      }
      std::vector<int64_t> sameBunchEvents = fBCCollMap[event.globalBC()];
      if (sameBunchEvents.size() > 1) { // event with in-bunch pileup
        evSel |= (static_cast<uint8_t>(1) << 1);
      }
      if (collisionSplittingMap.find(event.globalIndex()) != collisionSplittingMap.end()) { // event with possible fake in-bunch pileup (collision splitting)
        evSel |= (static_cast<uint8_t>(1) << 2);
      }
      eventSel(evSel);
    }
  }

  void processSkimmed(MyEvents const& events)
  {
    runEventSelection<gkEventFillMap>(events);
    publishSelections<gkEventFillMap>(events);
  }
  void processSkimmedWithZdc(MyEventsZdc const& events)
  {
    runEventSelection<gkEventFillMapWithZdc>(events);
    publishSelections<gkEventFillMapWithZdc>(events);
  }
  void processSkimmedWithMultExtra(MyEventsMultExtra const& events)
  {
    runEventSelection<gkEventFillMapWithMultExtra>(events);
    publishSelections<gkEventFillMapWithMultExtra>(events);
  }
  void processSkimmedWithMultExtraZdc(MyEventsMultExtraZdc const& events)
  {
    runEventSelection<gkEventFillMapWithMultExtraZdc>(events);
    publishSelections<gkEventFillMapWithMultExtraZdc>(events);
  }
  void processSkimmedWithMultExtraZdcFit(MyEventsMultExtraZdcFit const& events)
  {
    runEventSelection<gkEventFillMapWithMultExtraZdcFit>(events);
    publishSelections<gkEventFillMapWithMultExtraZdcFit>(events);
  }
  void processSkimmedWithQvectorCentr(MyEventsQvectorCentr const& events)
  {
    runEventSelection<gkEventFillMapWithQvectorCentr>(events);
    publishSelections<gkEventFillMapWithQvectorCentr>(events);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisEventSelection, processSkimmed, "Run event selection on DQ skimmed events", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedWithZdc, "Run event selection on DQ skimmed events, with ZDC", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedWithMultExtra, "Run event selection on DQ skimmed events, with mult extra", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedWithMultExtraZdc, "Run event selection on DQ skimmed events, with mult extra and ZDC", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedWithMultExtraZdcFit, "Run event selection on DQ skimmed events, with mult extra, ZDC and FIT", false);
  PROCESS_SWITCH(AnalysisEventSelection, processSkimmedWithQvectorCentr, "Run event selection on DQ skimmed events, with Q-vector", false);
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
  Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<bool> fConfigPublishAmbiguity{"cfgPublishAmbiguity", true, "If true, publish ambiguity table and fill QA histograms"};

  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/z/zhxiong/TPCPID/PostCalib", "base path to the ccdb object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas"};
  Configurable<int> fConfigTPCpostCalibType{"cfgTPCpostCalibType", 1, "1: (TPCncls,pIN,eta) calibration typically for pp, 2: (eta,nPV,nLong,tLong) calibration typically for PbPb"};
  Configurable<bool> fConfigTPCuseInterpolatedCalib{"cfgTPCpostCalibUseInterpolation", true, "If true, use interpolated calibration values (default: true)"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // Track related options
  Configurable<bool> fPropTrack{"cfgPropTrack", true, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut*> fTrackCuts;

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
        fTrackCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // Extra cuts via JSON
    TString addTrackCutsStr = fConfigCutsJSON.value;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (auto& t : addTrackCuts) {
        fTrackCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // set one histogram directory for each defined track cut
      TString histDirNames = "TrackBarrel_BeforeCuts;";
      for (auto& cut : fTrackCuts) {
        histDirNames += Form("TrackBarrel_%s;", cut->GetName());
      }
      if (fConfigPublishAmbiguity) {
        histDirNames += "TrackBarrel_AmbiguityInBunch;TrackBarrel_AmbiguityOutOfBunch;";
      }

      DefineHistograms(fHistMan, histDirNames.Data(), fConfigAddTrackHistogram.value.data()); // define all histograms
      dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str());  // ad-hoc histograms via JSON
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                        // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    fCCDBApi.init(fConfigCcdbUrl.value);
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTracks>
  void runTrackSelection(ReducedTracksAssoc const& assocs, TEvents const& events, TTracks const& tracks)
  {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      if (fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, events.begin().timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
        if (fConfigTPCpostCalibType == 2) {
          VarManager::SetCalibrationObject(VarManager::kTPCElectronStatus, calibList->FindObject("status_map_electron"));
          VarManager::SetCalibrationObject(VarManager::kTPCPionStatus, calibList->FindObject("status_map_pion"));
          VarManager::SetCalibrationObject(VarManager::kTPCProtonStatus, calibList->FindObject("status_map_proton"));
        }
        VarManager::SetCalibrationType(fConfigTPCpostCalibType, fConfigTPCuseInterpolatedCalib);
      }

      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
      }

      std::map<std::string, std::string> metadataRCT, header;
      header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", events.begin().runNumber()), metadataRCT, -1);
      uint64_t sor = std::atol(header["SOR"].c_str());
      uint64_t eor = std::atol(header["EOR"].c_str());
      VarManager::SetSORandEOR(sor, eor);

      fCurrentRun = events.begin().runNumber();
    }

    trackSel.reserve(assocs.size());
    trackAmbiguities.reserve(tracks.size());
    uint32_t filterMap = static_cast<uint32_t>(0);
    int iCut = 0;

    for (auto& assoc : assocs) {

      // if the event from this association is not selected, reject also the association
      auto event = assoc.template reducedevent_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        trackSel(0);
        continue;
      }
      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);

      auto track = assoc.template reducedtrack_as<TTracks>();
      filterMap = static_cast<uint32_t>(0);
      VarManager::FillTrack<TTrackFillMap>(track);
      // compute quantities which depend on the associated collision, such as DCA
      if (fPropTrack) {
        VarManager::FillTrackCollision<TTrackFillMap>(track, event);
      }
      if (fConfigQA) {
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }
      iCut = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, iCut++) {
        if ((*cut)->IsSelected(VarManager::fgValues)) {
          filterMap |= (static_cast<uint32_t>(1) << iCut);
          if (fConfigQA) {
            fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut)->GetName()), VarManager::fgValues);
          }
        }
      } // end loop over cuts

      // publish the decisions
      trackSel(filterMap);

      // count the number of associations per track
      if (fConfigPublishAmbiguity && filterMap > 0) {
        // for this track, count the number of associated collisions with in-bunch pileup and out of bunch associations
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
    } // end loop over associations

    if (fConfigPublishAmbiguity) {
      // QA the collision-track associations
      if (fConfigQA) {
        for (auto& [trackIdx, evIndices] : fNAssocsInBunch) {
          if (evIndices.size() == 1) {
            continue;
          }
          auto track = tracks.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
          VarManager::FillTrack<TTrackFillMap>(track);
          // Exceptionally, set the VarManager ambiguity number here, to be used in histograms
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
          // Exceptionally, set the VarManager ambiguity number here
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
    } // end if (fConfigPublishAmbiguity)

  } // end runTrackSelection()

  void processSkimmed(ReducedTracksAssoc const& assocs, MyEventsSelected const& events, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkEventFillMap, gkTrackFillMap>(assocs, events, tracks);
  }
  void processSkimmedWithMultExtra(ReducedTracksAssoc const& assocs, MyEventsMultExtraSelected const& events, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkEventFillMapWithMultExtra, gkTrackFillMap>(assocs, events, tracks);
  }
  void processSkimmedWithCov(ReducedTracksAssoc const& assocs, MyEventsVtxCovSelected const& events, MyBarrelTracksWithCov const& tracks)
  {
    runTrackSelection<gkEventFillMapWithCov, gkTrackFillMapWithCov>(assocs, events, tracks);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmed, "Run barrel track selection on DQ skimmed track associations", false);
  PROCESS_SWITCH(AnalysisTrackSelection, processSkimmedWithMultExtra, "Run barrel track selection on DQ skimmed track associations, with extra multiplicity tables", false);
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
  Configurable<std::string> fConfigGeoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut*> fMuonCuts;

  int fCurrentRun; // current run kept to detect run changes and trigger loading params from CCDB

  std::map<int64_t, std::vector<int64_t>> fNAssocsInBunch;    // key: muon global index, value: vector of global index for events associated in-bunch (events that have in-bunch pileup or splitting)
  std::map<int64_t, std::vector<int64_t>> fNAssocsOutOfBunch; // key: muon global index, value: vector of global index for events associated out-of-bunch (events that have no in-bunch pileup)

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
    // extra cuts from JSON
    TString addCutsStr = fConfigCutsJSON.value;
    if (addCutsStr != "") {
      std::vector<AnalysisCut*> addCuts = dqcuts::GetCutsFromJSON(addCutsStr.Data());
      for (auto& t : addCuts) {
        fMuonCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // set one histogram directory for each defined track cut
      TString histDirNames = "TrackMuon_BeforeCuts;";
      for (auto& cut : fMuonCuts) {
        histDirNames += Form("TrackMuon_%s;", cut->GetName());
      }
      if (fConfigPublishAmbiguity) {
        histDirNames += "TrackMuon_AmbiguityInBunch;TrackMuon_AmbiguityOutOfBunch;";
      }

      DefineHistograms(fHistMan, histDirNames.Data(), fConfigAddMuonHistogram.value.data()); // define all histograms
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
  void runMuonSelection(ReducedMuonsAssoc const& assocs, TEvents const& events, TMuons const& muons)
  {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();

    if (events.size() > 0 && fCurrentRun != events.begin().runNumber()) {
      o2::parameters::GRPMagField* grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, events.begin().timestamp());
      if (grpmag != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
        VarManager::SetMagneticField(grpmag->getNominalL3Field());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", events.begin().timestamp());
      }
      fCurrentRun = events.begin().runNumber();
    }

    muonSel.reserve(assocs.size());
    muonAmbiguities.reserve(muons.size());
    uint32_t filterMap = static_cast<uint32_t>(0);
    int iCut = 0;

    for (auto& assoc : assocs) {
      auto event = assoc.template reducedevent_as<TEvents>();
      if (!event.isEventSelected_bit(0)) {
        muonSel(0);
        continue;
      }
      VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
      // fill event information which might be needed in histograms/cuts that combine track and event properties
      VarManager::FillEvent<TEventFillMap>(event);

      auto track = assoc.template reducedmuon_as<TMuons>();
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
      } // end loop over cuts
      muonSel(filterMap);

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
      } // end if (fConfigPublishAmbiguity)
    } // end loop over assocs

    if (fConfigPublishAmbiguity) {
      // QA the collision-track associations
      if (fConfigQA) {
        for (auto& [trackIdx, evIndices] : fNAssocsInBunch) {
          if (evIndices.size() == 1) {
            continue;
          }
          auto track = muons.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
          VarManager::FillTrack<TMuonFillMap>(track);
          VarManager::fgValues[VarManager::kMuonNAssocsInBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackMuon_AmbiguityInBunch", VarManager::fgValues);
        } // end loop over in-bunch ambiguous tracks

        for (auto& [trackIdx, evIndices] : fNAssocsOutOfBunch) {
          if (evIndices.size() == 1) {
            continue;
          }
          auto track = muons.rawIteratorAt(trackIdx);
          VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
          VarManager::FillTrack<TMuonFillMap>(track);
          VarManager::fgValues[VarManager::kMuonNAssocsOutOfBunch] = static_cast<float>(evIndices.size());
          fHistMan->FillHistClass("TrackMuon_AmbiguityOutOfBunch", VarManager::fgValues);
        } // end loop over out-of-bunch ambiguous tracks
      }

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

  void processSkimmed(ReducedMuonsAssoc const& assocs, MyEventsVtxCovSelected const& events, MyMuonTracksWithCov const& muons)
  {
    runMuonSelection<gkEventFillMapWithCov, gkMuonFillMapWithCov>(assocs, events, muons);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisMuonSelection, processSkimmed, "Run muon selection on DQ skimmed muons", false);
  PROCESS_SWITCH(AnalysisMuonSelection, processDummy, "Dummy function", true);
};

// Run the prefilter selection (e.g. electron prefiltering for photon conversions)
// This takes uses a sample of tracks selected with loose cuts (fConfigPrefilterTrackCut) and combines them
//  with the sample of tracks to be used in downstream analysis (fConfigTrackCuts). If a pair is found to pass
//  the pair prefilter cut (cfgPrefilterPairCut), the analysis track is tagged to be removed from analysis.
struct AnalysisPrefilterSelection {
  Produces<aod::Prefilter> prefilter; // joinable with ReducedTracksAssoc

  // Configurables
  Configurable<std::string> fConfigPrefilterTrackCut{"cfgPrefilterTrackCut", "", "Prefilter track cut"};
  Configurable<std::string> fConfigPrefilterPairCut{"cfgPrefilterPairCut", "", "Prefilter pair cut"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Track cuts for which to run the prefilter"};

  // TODO: Add prefilter pair cut via JSON

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
    VarManager::SetDefaultVarNames();

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

    fPrefilterMask = static_cast<uint32_t>(0);
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

    VarManager::SetUseVars(AnalysisCut::fgUsedVars);                                   // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    // VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  template <uint32_t TTrackFillMap, typename TTracks>
  void runPrefilter(soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& assocs, TTracks const& /*tracks*/)
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

      // Here we check if we should apply the prefilter for this pair
      if (!((track1Candidate > 0 && track2Loose) || (track2Candidate > 0 && track1Loose))) {
        continue;
      }

      // compute pair quantities
      VarManager::FillPair<VarManager::kDecayToEE, TTrackFillMap>(track1, track2);
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
        runPrefilter<gkTrackFillMap>(groupedAssocs, tracks);
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
  Produces<aod::DimuonsExtra> dimuonsExtraList;
  Produces<aod::DielectronsAll> dielectronAllList;
  Produces<aod::DimuonsAll> dimuonAllList;
  Produces<aod::DileptonFlow> dileptonFlowList;
  Produces<aod::DileptonsInfo> dileptonInfoList;
  Produces<aod::JPsieeCandidates> PromptNonPromptSepTable;
  Produces<aod::DileptonPolarization> dileptonPolarList;

  o2::base::MatLayerCylSet* fLUT = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  OutputObj<THashList> fOutputList{"output"};

  struct : ConfigurableGroup {
    Configurable<std::string> track{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
    Configurable<std::string> muon{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
    Configurable<std::string> pair{"cfgPairCuts", "", "Comma separated list of pair cuts"};
    Configurable<bool> event{"cfgRemoveCollSplittingCandidates", false, "If true, remove collision splitting candidates as determined by the event selection task upstream"};
    // TODO: Add pair cuts via JSON
  } fConfigCuts;

  Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 100, "Number of Events stored for event mixing"};
  // Configurable<std::string> fConfigAddEventMixingHistogram{"cfgAddEventMixingHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
  Configurable<bool> fConfigQA{"cfgQA", true, "If true, fill output histograms"};
  Configurable<bool> fConfigAmbiguousMuonHistograms{"cfgAmbiguousMuonHistograms", true, "If true, fill ambiguous histograms"};

  struct : ConfigurableGroup {
    Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> GrpLhcIfPath{"grplhcif", "GLO/Config/GRPLHCIF", "Path on the CCDB for the GRPLHCIF object"};
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
    Configurable<bool> propTrack{"cfgPropTrack", true, "Propgate tracks to associated collision to recalculate DCA and momentum vector"};
    Configurable<bool> useRemoteCollisionInfo{"cfgUseRemoteCollisionInfo", false, "Use remote collision information from CCDB"};
  } fConfigOptions;
  struct : ConfigurableGroup {
    Configurable<bool> applyBDT{"applyBDT", false, "Flag to apply ML selections"};
    Configurable<std::string> fConfigBdtCutsJSON{"fConfigBdtCutsJSON", "", "Additional list of BDT cuts in JSON format"};

    Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/j/jseo/ML/PbPbPsi/default/"}, "Paths of models on CCDB"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  } fConfigML;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  Filter filterEventSelected = aod::dqanalysisflags::isEventSelected > static_cast<uint8_t>(0);

  HistogramManager* fHistMan;

  o2::analysis::DQMlResponse<float> fDQMlResponse;
  std::vector<float> fOutputMlPsi2ee = {}; // TODO: check this is needed or not

  // keep histogram class names in maps, so we don't have to buld their names in the pair loops
  std::map<int, std::vector<TString>> fTrackHistNames;
  std::map<int, std::vector<TString>> fMuonHistNames;
  std::map<int, std::vector<TString>> fTrackMuonHistNames;
  std::vector<AnalysisCompositeCut> fPairCuts;
  std::vector<TString> fTrackCuts;
  std::map<std::pair<uint32_t, uint32_t>, uint32_t> fAmbiguousPairs;

  uint32_t fTrackFilterMask; // mask for the track cuts required in this task to be applied on the barrel cuts produced upstream
  uint32_t fMuonFilterMask;  // mask for the muon cuts required in this task to be applied on the muon cuts produced upstream
  int fNCutsBarrel;
  int fNCutsMuon;
  int fNPairCuts;

  bool fEnableBarrelMixingHistos;
  bool fEnableBarrelHistos;
  bool fEnableMuonHistos;
  bool fEnableMuonMixingHistos;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> hashBin;

  Preslice<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter>> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

  void init(o2::framework::InitContext& context)
  {
    fEnableBarrelHistos = context.mOptions.get<bool>("processAllSkimmed") || context.mOptions.get<bool>("processBarrelOnlySkimmed") || context.mOptions.get<bool>("processBarrelOnlyWithCollSkimmed") || context.mOptions.get<bool>("processBarrelOnlySkimmedNoCov") || context.mOptions.get<bool>("processBarrelOnlySkimmedNoCovWithMultExtra") || context.mOptions.get<bool>("processBarrelOnlyWithQvectorCentrSkimmedNoCov");
    fEnableBarrelMixingHistos = context.mOptions.get<bool>("processMixingAllSkimmed") || context.mOptions.get<bool>("processMixingBarrelSkimmed") || context.mOptions.get<bool>("processMixingBarrelSkimmedFlow") || context.mOptions.get<bool>("processMixingBarrelWithQvectorCentrSkimmedNoCov");
    fEnableMuonHistos = context.mOptions.get<bool>("processAllSkimmed") || context.mOptions.get<bool>("processMuonOnlySkimmed") || context.mOptions.get<bool>("processMuonOnlySkimmedMultExtra") || context.mOptions.get<bool>("processMuonOnlySkimmedFlow") || context.mOptions.get<bool>("processMixingMuonSkimmed");
    fEnableMuonMixingHistos = context.mOptions.get<bool>("processMixingAllSkimmed") || context.mOptions.get<bool>("processMixingMuonSkimmed");

    if (context.mOptions.get<bool>("processDummy")) {
      if (fEnableBarrelHistos || fEnableBarrelMixingHistos || fEnableMuonHistos || fEnableMuonMixingHistos) {
        LOG(fatal) << "No other processing tasks should be enabled if the processDummy is enabled!!";
      }
      return;
    }
    VarManager::SetDefaultVarNames();

    // Keep track of all the histogram class names to avoid composing strings in the pairing loop
    TString histNames = "";
    std::vector<TString> names;
    fTrackCuts.clear();

    // NOTE: Pair cuts are only applied on the histogram output. The produced pair tables do not have these cuts applied
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

    if (fConfigML.applyBDT) {
      // BDT cuts via JSON
      std::vector<double> binsMl;
      o2::framework::LabeledArray<double> cutsMl;
      std::vector<int> cutDirMl;
      int nClassesMl = 1;
      std::vector<std::string> namesInputFeatures;
      std::vector<std::string> onnxFileNames;

      auto config = o2::aod::dqmlcuts::GetBdtScoreCutsAndConfigFromJSON(fConfigML.fConfigBdtCutsJSON.value.c_str());

      if (std::holds_alternative<dqmlcuts::BinaryBdtScoreConfig>(config)) {
        auto& cfg = std::get<dqmlcuts::BinaryBdtScoreConfig>(config);
        binsMl = cfg.binsMl;
        nClassesMl = 1;
        cutsMl = cfg.cutsMl;
        cutDirMl = cfg.cutDirs;
        namesInputFeatures = cfg.inputFeatures;
        onnxFileNames = cfg.onnxFiles;
        fDQMlResponse.setBinsCent(cfg.binsCent);
        fDQMlResponse.setBinsPt(cfg.binsPt);
        fDQMlResponse.setCentType(cfg.centType);
        LOG(info) << "Using BDT cuts for binary classification";
      } else {
        auto& cfg = std::get<dqmlcuts::MultiClassBdtScoreConfig>(config);
        binsMl = cfg.binsMl;
        nClassesMl = 3;
        cutsMl = cfg.cutsMl;
        cutDirMl = cfg.cutDirs;
        namesInputFeatures = cfg.inputFeatures;
        onnxFileNames = cfg.onnxFiles;
        fDQMlResponse.setBinsCent(cfg.binsCent);
        fDQMlResponse.setBinsPt(cfg.binsPt);
        fDQMlResponse.setCentType(cfg.centType);
        LOG(info) << "Using BDT cuts for multiclass classification";
      }

      fDQMlResponse.configure(binsMl, cutsMl, cutDirMl, nClassesMl);
      if (fConfigML.loadModelsFromCCDB) {
        fCCDBApi.init(fConfigCCDB.url);
        fDQMlResponse.setModelPathsCCDB(onnxFileNames, fCCDBApi, fConfigML.modelPathsCCDB, fConfigML.timestampCCDB);
      } else {
        fDQMlResponse.setModelPathsLocal(onnxFileNames);
      }
      fDQMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      fDQMlResponse.init();
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
    if (!trackCutsStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
      fNCutsBarrel = objArray->GetEntries();
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        fTrackCuts.push_back(tempStr);
        if (objArrayTrackCuts->FindObject(tempStr.Data()) != nullptr) {
          fTrackFilterMask |= (static_cast<uint32_t>(1) << icut);

          if (fEnableBarrelHistos) {
            names = {
              Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
              Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
            names.push_back(Form("PairsBarrelSEPM_ambiguousextra_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsBarrelSEPP_ambiguousextra_%s", objArray->At(icut)->GetName()));
            names.push_back(Form("PairsBarrelSEMM_ambiguousextra_%s", objArray->At(icut)->GetName()));
            histNames += Form("%s;%s;%s;", names[3].Data(), names[4].Data(), names[5].Data());
            if (fEnableBarrelMixingHistos) {
              names.push_back(Form("PairsBarrelMEPM_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelMEPP_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsBarrelMEMM_%s", objArray->At(icut)->GetName()));
              histNames += Form("%s;%s;%s;", names[6].Data(), names[7].Data(), names[8].Data());
            }
            fTrackHistNames[icut] = names;

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
                fTrackHistNames[fNCutsBarrel + icut * fNPairCuts + iPairCut] = names;
              } // end loop (pair cuts)
            } // end if (pair cuts)
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

    if (!muonCutsStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(tempCutsStr.Tokenize(","));
      fNCutsMuon = objArray->GetEntries();
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        if (objArrayMuonCuts->FindObject(tempStr.Data()) != nullptr) {
          fMuonFilterMask |= (static_cast<uint32_t>(1) << icut);

          if (fEnableMuonHistos) {
            // no pair cuts
            names = {
              Form("PairsMuonSEPM_%s", objArray->At(icut)->GetName()),
              Form("PairsMuonSEPP_%s", objArray->At(icut)->GetName()),
              Form("PairsMuonSEMM_%s", objArray->At(icut)->GetName())};
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
            if (fEnableMuonMixingHistos) {
              names.push_back(Form("PairsMuonMEPM_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonMEPP_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonMEMM_%s", objArray->At(icut)->GetName()));
              histNames += Form("%s;%s;%s;", names[3].Data(), names[4].Data(), names[5].Data());
            }
            if (fConfigAmbiguousMuonHistograms) {
              names.push_back(Form("PairsMuonSEPM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEPP_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEMM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEPM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEPP_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEMM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEPM_unambiguous_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEPP_unambiguous_%s", objArray->At(icut)->GetName()));
              names.push_back(Form("PairsMuonSEMM_unambiguous_%s", objArray->At(icut)->GetName()));
              histNames += Form("%s;%s;%s;", names[(fEnableMuonMixingHistos ? 6 : 3)].Data(), names[(fEnableMuonMixingHistos ? 7 : 4)].Data(), names[(fEnableMuonMixingHistos ? 8 : 5)].Data());
              histNames += Form("%s;%s;%s;", names[(fEnableMuonMixingHistos ? 9 : 6)].Data(), names[(fEnableMuonMixingHistos ? 10 : 7)].Data(), names[(fEnableMuonMixingHistos ? 11 : 8)].Data());
              histNames += Form("%s;%s;%s;", names[(fEnableMuonMixingHistos ? 12 : 9)].Data(), names[(fEnableMuonMixingHistos ? 13 : 10)].Data(), names[(fEnableMuonMixingHistos ? 14 : 11)].Data());
              if (fEnableMuonMixingHistos) {
                names.push_back(Form("PairsMuonMEPM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
                names.push_back(Form("PairsMuonMEPP_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
                names.push_back(Form("PairsMuonMEMM_ambiguousInBunch_%s", objArray->At(icut)->GetName()));
                names.push_back(Form("PairsMuonMEPM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
                names.push_back(Form("PairsMuonMEPP_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
                names.push_back(Form("PairsMuonMEMM_ambiguousOutOfBunch_%s", objArray->At(icut)->GetName()));
                names.push_back(Form("PairsMuonMEPM_unambiguous_%s", objArray->At(icut)->GetName()));
                names.push_back(Form("PairsMuonMEPP_unambiguous_%s", objArray->At(icut)->GetName()));
                names.push_back(Form("PairsMuonMEMM_unambiguous_%s", objArray->At(icut)->GetName()));
                histNames += Form("%s;%s;%s;", names[15].Data(), names[16].Data(), names[17].Data());
                histNames += Form("%s;%s;%s;", names[18].Data(), names[19].Data(), names[20].Data());
                histNames += Form("%s;%s;%s;", names[21].Data(), names[22].Data(), names[23].Data());
              }
            }
            fMuonHistNames[icut] = names;

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
          }
        }
      }
    }

    fCurrentRun = 0;

    fCCDB->setURL(fConfigCCDB.url.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDBApi.init(fConfigCCDB.url.value);

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

    /*if (context.mOptions.get<bool>("processElectronMuonSkimmed") || context.mOptions.get<bool>("processAllSkimmed")) {
      TString cutNamesBarrel = fConfigCuts.track.value;
      TString cutNamesMuon = fConfigCuts.muon.value;
      if (!cutNamesBarrel.IsNull() && !cutNamesMuon.IsNull()) {
        std::unique_ptr<TObjArray> objArrayBarrel(cutNamesBarrel.Tokenize(","));
        std::unique_ptr<TObjArray> objArrayMuon(cutNamesMuon.Tokenize(","));
        if (objArrayBarrel->GetEntries() == objArrayMuon->GetEntries()) {   // one must specify equal number of barrel and muon cuts
          for (int icut = 0; icut < objArrayBarrel->GetEntries(); ++icut) { // loop over track cuts
            // no pair cuts
            names = {
              Form("PairsEleMuSEPM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuSEPP_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName()),
              Form("PairsEleMuSEMM_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName())};
            histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
            fTrackMuonHistNames.push_back(names);

            TString cutNamesStr = fConfigCuts.pair.value;
            if (!cutNamesStr.IsNull()) { // if pair cuts
              std::unique_ptr<TObjArray> objArrayPair(cutNamesStr.Tokenize(","));
              for (int iPairCut = 0; iPairCut < objArrayPair->GetEntries(); ++iPairCut) { // loop over pair cuts
                std::vector<TString> names = {
                  Form("PairsEleMuSEPM_%s_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsEleMuSEPP_%s_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName()),
                  Form("PairsEleMuSEMM_%s_%s_%s", objArrayBarrel->At(icut)->GetName(), objArrayMuon->At(icut)->GetName(), objArrayPair->At(iPairCut)->GetName())};
                histNames += Form("%s;%s;%s;", names[0].Data(), names[1].Data(), names[2].Data());
                fTrackMuonHistNames.push_back(names);
              } // end loop (pair cuts)
            }   // end if (pair cuts)
          }     // end loop (track cuts)
        }       // end if (equal number of cuts)
      }         // end if (track cuts)
    }*/

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(true);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      VarManager::SetCollisionSystem((TString)fConfigOptions.collisionSystem, fConfigOptions.centerMassEnergy); // set collision system and center of mass energy
      DefineHistograms(fHistMan, histNames.Data(), fConfigAddSEPHistogram.value.data());                        // define all histograms
      dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str());                    // ad-hoc histograms via JSON
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                                          // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  void initParamsFromCCDB(uint64_t timestamp, int runNumber, bool withTwoProngFitter = true)
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

    std::map<string, string> metadataRCT, header;
    header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", runNumber), metadataRCT, -1);
    uint64_t sor = std::atol(header["SOR"].c_str());
    uint64_t eor = std::atol(header["EOR"].c_str());
    VarManager::SetSORandEOR(sor, eor);

    if (fConfigOptions.useRemoteCollisionInfo) {
      o2::parameters::GRPLHCIFData* grpo = fCCDB->getForTimeStamp<o2::parameters::GRPLHCIFData>(fConfigCCDB.GrpLhcIfPath, timestamp);
      VarManager::SetCollisionSystem(grpo);
    }
  }

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <bool TTwoProngFitter, int TPairType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvents, typename TTrackAssocs, typename TTracks>
  void runSameEventPairing(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& assocs, TTracks const& /*tracks*/)
  {
    if (events.size() > 0) { // Additional protection to avoid crashing of events.begin().runNumber()
      if (fCurrentRun != events.begin().runNumber()) {
        initParamsFromCCDB(events.begin().timestamp(), events.begin().runNumber(), TTwoProngFitter);
        fCurrentRun = events.begin().runNumber();
      }
    }

    TString cutNames = fConfigCuts.track.value;
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    int ncuts = fNCutsBarrel;
    int histIdxOffset = 0;
    if constexpr (TPairType == pairTypeMuMu) {
      cutNames = fConfigCuts.muon.value;
      histNames = fMuonHistNames;
      ncuts = fNCutsMuon;
      if (fEnableMuonMixingHistos) {
        histIdxOffset = 3;
      }
    }
    if constexpr (TPairType == pairTypeEE) {
      if (fEnableBarrelMixingHistos) {
        histIdxOffset = 3;
      }
    }
    /*if constexpr (TPairType == pairTypeEMu) {
      cutNames = fConfigCuts.muon.value;
      histNames = fTrackMuonHistNames;
    }*/

    uint32_t twoTrackFilter = static_cast<uint32_t>(0);
    uint32_t dileptonMcDecision = static_cast<uint32_t>(0); // placeholder, copy of the dqEfficiency.cxx one
    int sign1 = 0;
    int sign2 = 0;
    dielectronList.reserve(1);
    dimuonList.reserve(1);
    dielectronsExtraList.reserve(1);
    dielectronInfoList.reserve(1);
    dimuonsExtraList.reserve(1);
    dileptonInfoList.reserve(1);
    dileptonFlowList.reserve(1);
    if (fConfigOptions.flatTables.value) {
      dielectronAllList.reserve(1);
      dimuonAllList.reserve(1);
    }
    if (fConfigOptions.polarTables.value) {
      dileptonPolarList.reserve(1);
    }
    fAmbiguousPairs.clear();
    constexpr bool eventHasQvector = ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0);
    constexpr bool eventHasQvectorCentr = ((TEventFillMap & VarManager::ObjTypes::CollisionQvect) > 0);
    constexpr bool trackHasCov = ((TTrackFillMap & VarManager::ObjTypes::TrackCov) > 0 || (TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelCov) > 0);
    bool isSelectedBDT = false;

    for (auto& event : events) {
      if (!event.isEventSelected_bit(0)) {
        continue;
      }
      if (fConfigCuts.event && event.isEventSelected_bit(2)) {
        continue;
      }
      uint8_t evSel = event.isEventSelected_raw();
      // Reset the fValues array
      VarManager::ResetValues(0, VarManager::kNVars);
      // VarManager::FillEvent<gkEventFillMap>(event, VarManager::fgValues);
      VarManager::FillEvent<TEventFillMap>(event, VarManager::fgValues);

      auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());
      if (groupedAssocs.size() == 0) {
        continue;
      }

      bool isFirst = true;
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
          // TODO: Make sure that we do not work with more than 28 track bits
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
          // compute quantities which depend on the associated collision, such as DCA
          if (fConfigOptions.propTrack) {
            VarManager::FillPairCollision<TPairType, TTrackFillMap>(event, t1, t2);
          }
          if constexpr (TTwoProngFitter) {
            VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigOptions.propToPCA);
          }
          if constexpr (eventHasQvector) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }
          if constexpr (eventHasQvectorCentr) {
            VarManager::FillPairVn<TEventFillMap, TPairType>(t1, t2);
          }

          dielectronList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                         VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                         t1.sign() + t2.sign(), twoTrackFilter, 0);

          if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackCollInfo) > 0) {
            dielectronInfoList(t1.collisionId(), t1.trackId(), t2.trackId());
            dileptonInfoList(t1.collisionId(), event.posX(), event.posY(), event.posZ());
          }
          if (fConfigOptions.polarTables.value) {
            dileptonPolarList(VarManager::fgValues[VarManager::kCosThetaHE], VarManager::fgValues[VarManager::kPhiHE], VarManager::fgValues[VarManager::kPhiTildeHE],
                              VarManager::fgValues[VarManager::kCosThetaCS], VarManager::fgValues[VarManager::kPhiCS], VarManager::fgValues[VarManager::kPhiTildeCS],
                              VarManager::fgValues[VarManager::kCosThetaPP], VarManager::fgValues[VarManager::kPhiPP], VarManager::fgValues[VarManager::kPhiTildePP],
                              VarManager::fgValues[VarManager::kCosThetaRM],
                              VarManager::fgValues[VarManager::kCosThetaStarTPC], VarManager::fgValues[VarManager::kCosThetaStarFT0A], VarManager::fgValues[VarManager::kCosThetaStarFT0C]);
          }
          if constexpr (trackHasCov && TTwoProngFitter) {
            dielectronsExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
            if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelPID) > 0) {
              if (fConfigML.applyBDT) {
                std::vector<float> dqInputFeatures = fDQMlResponse.getInputFeatures(t1, t2, VarManager::fgValues);

                if (dqInputFeatures.empty()) {
                  LOG(fatal) << "Input features for ML selection are empty! Please check your configuration.";
                  return;
                }

                int modelIndex = -1;
                const auto& binsCent = fDQMlResponse.getBinsCent();
                const auto& binsPt = fDQMlResponse.getBinsPt();
                const std::string& centType = fDQMlResponse.getCentType();

                if ("kCentFT0C" == centType) {
                  modelIndex = o2::aod::dqmlcuts::getMlBinIndex(VarManager::fgValues[VarManager::kCentFT0C], VarManager::fgValues[VarManager::kPt], binsCent, binsPt);
                } else if ("kCentFT0A" == centType) {
                  modelIndex = o2::aod::dqmlcuts::getMlBinIndex(VarManager::fgValues[VarManager::kCentFT0A], VarManager::fgValues[VarManager::kPt], binsCent, binsPt);
                } else if ("kCentFT0M" == centType) {
                  modelIndex = o2::aod::dqmlcuts::getMlBinIndex(VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kPt], binsCent, binsPt);
                } else {
                  LOG(fatal) << "Unknown centrality estimation type: " << centType;
                  return;
                }

                if (modelIndex < 0) {
                  LOG(info) << "Ml index is negative! This means that the centrality/pt is not in the range of the model bins.";
                  continue;
                }

                LOG(debug) << "Model index: " << modelIndex << ", pT: " << VarManager::fgValues[VarManager::kPt] << ", centrality (kCentFT0C): " << VarManager::fgValues[VarManager::kCentFT0C];
                isSelectedBDT = fDQMlResponse.isSelectedMl(dqInputFeatures, modelIndex, fOutputMlPsi2ee);
                VarManager::FillBdtScore(fOutputMlPsi2ee); // TODO: check if this is needed or not
              }

              if (fConfigML.applyBDT && !isSelectedBDT)
                continue;

              if (fConfigOptions.flatTables.value) {
                dielectronAllList(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), twoTrackFilter, dileptonMcDecision,
                                  t1.pt(), t1.eta(), t1.phi(), t1.itsClusterMap(), t1.itsChi2NCl(), t1.tpcNClsCrossedRows(), t1.tpcNClsFound(), t1.tpcChi2NCl(), t1.dcaXY(), t1.dcaZ(), t1.tpcSignal(), t1.tpcNSigmaEl(), t1.tpcNSigmaPi(), t1.tpcNSigmaPr(), t1.beta(), t1.tofNSigmaEl(), t1.tofNSigmaPi(), t1.tofNSigmaPr(),
                                  t2.pt(), t2.eta(), t2.phi(), t2.itsClusterMap(), t2.itsChi2NCl(), t2.tpcNClsCrossedRows(), t2.tpcNClsFound(), t2.tpcChi2NCl(), t2.dcaXY(), t2.dcaZ(), t2.tpcSignal(), t2.tpcNSigmaEl(), t2.tpcNSigmaPi(), t2.tpcNSigmaPr(), t2.beta(), t2.tofNSigmaEl(), t2.tofNSigmaPi(), t2.tofNSigmaPr(),
                                  VarManager::fgValues[VarManager::kKFTrack0DCAxyz], VarManager::fgValues[VarManager::kKFTrack1DCAxyz], VarManager::fgValues[VarManager::kKFDCAxyzBetweenProngs], VarManager::fgValues[VarManager::kKFTrack0DCAxy], VarManager::fgValues[VarManager::kKFTrack1DCAxy], VarManager::fgValues[VarManager::kKFDCAxyBetweenProngs],
                                  VarManager::fgValues[VarManager::kKFTrack0DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationFromPV], VarManager::fgValues[VarManager::kKFTrack0DeviationxyFromPV], VarManager::fgValues[VarManager::kKFTrack1DeviationxyFromPV],
                                  VarManager::fgValues[VarManager::kKFMass], VarManager::fgValues[VarManager::kKFChi2OverNDFGeo], VarManager::fgValues[VarManager::kVertexingLxyz], VarManager::fgValues[VarManager::kVertexingLxyzOverErr], VarManager::fgValues[VarManager::kVertexingLxy], VarManager::fgValues[VarManager::kVertexingLxyOverErr], VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr], VarManager::fgValues[VarManager::kKFCosPA], VarManager::fgValues[VarManager::kKFJpsiDCAxyz], VarManager::fgValues[VarManager::kKFJpsiDCAxy],
                                  VarManager::fgValues[VarManager::kKFPairDeviationFromPV], VarManager::fgValues[VarManager::kKFPairDeviationxyFromPV],
                                  VarManager::fgValues[VarManager::kKFMassGeoTop],
                                  VarManager::fgValues[VarManager::kKFChi2OverNDFGeoTop], VarManager::fgValues[VarManager::kVertexingTauzProjected], VarManager::fgValues[VarManager::kVertexingTauxyProjected], VarManager::fgValues[VarManager::kVertexingLzProjected], VarManager::fgValues[VarManager::kVertexingLxyProjected]);
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

          VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
          // compute quantities which depend on the associated collision, such as DCA
          if (fConfigOptions.propTrack) {
            VarManager::FillPairCollision<TPairType, TTrackFillMap>(event, t1, t2);
          }
          if constexpr (TTwoProngFitter) {
            VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigOptions.propToPCA);
          }
          if constexpr (eventHasQvector) {
            VarManager::FillPairVn<TPairType>(t1, t2);
          }
          if constexpr (eventHasQvectorCentr) {
            VarManager::FillPairVn<TEventFillMap, TPairType>(t1, t2);
          }

          dimuonList(event.globalIndex(), VarManager::fgValues[VarManager::kMass],
                     VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi],
                     t1.sign() + t2.sign(), twoTrackFilter, 0);
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedMuonCollInfo) > 0) {
            dileptonInfoList(t1.collisionId(), event.posX(), event.posY(), event.posZ());
          }

          if constexpr (TTwoProngFitter) {
            dimuonsExtraList(t1.globalIndex(), t2.globalIndex(), VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingLz], VarManager::fgValues[VarManager::kVertexingLxy]);
            if (fConfigOptions.flatTables.value) {
              dimuonAllList(event.posX(), event.posY(), event.posZ(), event.numContrib(),
                            event.selection_raw(), evSel,
                            -999., -999., -999.,
                            VarManager::fgValues[VarManager::kMass],
                            false,
                            VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), VarManager::fgValues[VarManager::kVertexingChi2PCA],
                            VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingTauzErr],
                            VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr],
                            VarManager::fgValues[VarManager::kCosPointingAngle],
                            VarManager::fgValues[VarManager::kPt1], VarManager::fgValues[VarManager::kEta1], VarManager::fgValues[VarManager::kPhi1], t1.sign(),
                            VarManager::fgValues[VarManager::kPt2], VarManager::fgValues[VarManager::kEta2], VarManager::fgValues[VarManager::kPhi2], t2.sign(),
                            t1.fwdDcaX(), t1.fwdDcaY(), t2.fwdDcaX(), t2.fwdDcaY(),
                            0., 0.,
                            t1.chi2MatchMCHMID(), t2.chi2MatchMCHMID(),
                            t1.chi2MatchMCHMFT(), t2.chi2MatchMCHMFT(),
                            t1.chi2(), t2.chi2(),
                            -999., -999., -999., -999.,
                            -999., -999., -999., -999.,
                            -999., -999., -999., -999.,
                            -999., -999., -999., -999.,
                            (twoTrackFilter & (static_cast<uint32_t>(1) << 28)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 29)), (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31)),
                            VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kU3Q3],
                            VarManager::fgValues[VarManager::kR2EP_AB], VarManager::fgValues[VarManager::kR2SP_AB], VarManager::fgValues[VarManager::kCentFT0C],
                            VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kCos3DeltaPhi],
                            VarManager::fgValues[VarManager::kCORR2POI], VarManager::fgValues[VarManager::kCORR4POI], VarManager::fgValues[VarManager::kM01POI], VarManager::fgValues[VarManager::kM0111POI], VarManager::fgValues[VarManager::kMultDimuons],
                            VarManager::fgValues[VarManager::kVertexingPz], VarManager::fgValues[VarManager::kVertexingSV]);
            }
            if constexpr ((TTrackFillMap & VarManager::ObjTypes::ReducedMuonCollInfo) > 0) {
              if constexpr (eventHasQvector == true || eventHasQvectorCentr == true) {
                dileptonFlowList(t1.collisionId(), VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kCentFT0C],
                                 VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), isFirst,
                                 VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kR2SP_AB], VarManager::fgValues[VarManager::kR2SP_AC], VarManager::fgValues[VarManager::kR2SP_BC],
                                 VarManager::fgValues[VarManager::kU3Q3], VarManager::fgValues[VarManager::kR3SP],
                                 VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kR2EP_AB], VarManager::fgValues[VarManager::kR2EP_AC], VarManager::fgValues[VarManager::kR2EP_BC],
                                 VarManager::fgValues[VarManager::kCos3DeltaPhi], VarManager::fgValues[VarManager::kR3EP],
                                 VarManager::fgValues[VarManager::kCORR2POI], VarManager::fgValues[VarManager::kCORR4POI], VarManager::fgValues[VarManager::kM01POI], VarManager::fgValues[VarManager::kM0111POI],
                                 VarManager::fgValues[VarManager::kCORR2REF], VarManager::fgValues[VarManager::kCORR4REF], VarManager::fgValues[VarManager::kM11REF], VarManager::fgValues[VarManager::kM1111REF],
                                 VarManager::fgValues[VarManager::kMultDimuons], VarManager::fgValues[VarManager::kMultA]);
              }
            }
          }
          if (t1.sign() != t2.sign()) {
            isFirst = false;
          }
        }
        // TODO: the model for the electron-muon combination has to be thought through
        /*if constexpr (TPairType == VarManager::kElectronMuon) {
          twoTrackFilter = a1.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isMuonSelected_raw() & fTwoTrackFilterMask;
        }*/

        // Fill histograms
        bool isAmbiInBunch = false;
        bool isAmbiOutOfBunch = false;
        bool isUnambiguous = false;
        bool isLeg1Ambi = false;
        bool isLeg2Ambi = false;
        bool isAmbiExtra = false;

        if (fConfigML.applyBDT && !isSelectedBDT)
          continue;

        for (int icut = 0; icut < ncuts; icut++) {
          if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
            isAmbiInBunch = (twoTrackFilter & (static_cast<uint32_t>(1) << 28)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 29));
            isAmbiOutOfBunch = (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31));
            isUnambiguous = !(isAmbiInBunch || isAmbiOutOfBunch);
            isLeg1Ambi = (twoTrackFilter & (static_cast<uint32_t>(1) << 28) || (twoTrackFilter & (static_cast<uint32_t>(1) << 30)));
            isLeg2Ambi = (twoTrackFilter & (static_cast<uint32_t>(1) << 29) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31)));
            if constexpr (TPairType == VarManager::kDecayToEE) {
              if (isLeg1Ambi && isLeg2Ambi) {
                std::pair<uint32_t, uint32_t> iPair(a1.reducedtrackId(), a2.reducedtrackId());
                if (fAmbiguousPairs.find(iPair) != fAmbiguousPairs.end()) {
                  if (fAmbiguousPairs[iPair] & (static_cast<uint32_t>(1) << icut)) { // if this pair is already stored with this cut
                    isAmbiExtra = true;
                  } else {
                    fAmbiguousPairs[iPair] |= static_cast<uint32_t>(1) << icut;
                  }
                } else {
                  fAmbiguousPairs[iPair] = static_cast<uint32_t>(1) << icut;
                }
              }
            }
            if (sign1 * sign2 < 0) {
              PromptNonPromptSepTable(VarManager::fgValues[VarManager::kMass], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kVertexingTauxyProjected], VarManager::fgValues[VarManager::kVertexingTauxyProjectedPoleJPsiMass], VarManager::fgValues[VarManager::kVertexingTauzProjected], isAmbiInBunch, isAmbiOutOfBunch, VarManager::fgValues[VarManager::kMultFT0A], VarManager::fgValues[VarManager::kMultFT0C], VarManager::fgValues[VarManager::kCentFT0M], VarManager::fgValues[VarManager::kVtxNcontribReal]);
              if constexpr (TPairType == VarManager::kDecayToMuMu) {
                fHistMan->FillHistClass(histNames[icut][0].Data(), VarManager::fgValues);
		if (fConfigAmbiguousMuonHistograms) {
                if (isAmbiInBunch) {
                  fHistMan->FillHistClass(histNames[icut][3 + histIdxOffset].Data(), VarManager::fgValues);
                }
                if (isAmbiOutOfBunch) {
                  fHistMan->FillHistClass(histNames[icut][3 + histIdxOffset + 3].Data(), VarManager::fgValues);
                }
                if (isUnambiguous) {
                  fHistMan->FillHistClass(histNames[icut][3 + histIdxOffset + 6].Data(), VarManager::fgValues);
                }
		}
              }
              if constexpr (TPairType == VarManager::kDecayToEE) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s", fTrackCuts[icut].Data()), VarManager::fgValues);
                if (isAmbiExtra) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPM_ambiguousextra_%s", fTrackCuts[icut].Data()), VarManager::fgValues);
                }
              }
            } else {
              if (sign1 > 0) {
                if constexpr (TPairType == VarManager::kDecayToMuMu) {
                  fHistMan->FillHistClass(histNames[icut][1].Data(), VarManager::fgValues);
		  if (fConfigAmbiguousMuonHistograms) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][4 + histIdxOffset].Data(), VarManager::fgValues);
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][4 + histIdxOffset + 3].Data(), VarManager::fgValues);
                  }
                  if (isUnambiguous) {
                    fHistMan->FillHistClass(histNames[icut][4 + histIdxOffset + 6].Data(), VarManager::fgValues);
                  }
		  }
                }
                if constexpr (TPairType == VarManager::kDecayToEE) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s", fTrackCuts[icut].Data()), VarManager::fgValues);
                  if (isAmbiExtra) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEPP_ambiguousextra_%s", fTrackCuts[icut].Data()), VarManager::fgValues);
                  }
                }
              } else {
                if constexpr (TPairType == VarManager::kDecayToMuMu) {
                  fHistMan->FillHistClass(histNames[icut][2].Data(), VarManager::fgValues);
		  if (fConfigAmbiguousMuonHistograms) {
                  if (isAmbiInBunch) {
                    fHistMan->FillHistClass(histNames[icut][5 + histIdxOffset].Data(), VarManager::fgValues);
                  }
                  if (isAmbiOutOfBunch) {
                    fHistMan->FillHistClass(histNames[icut][5 + histIdxOffset + 3].Data(), VarManager::fgValues);
                  }
                  if (isUnambiguous) {
                    fHistMan->FillHistClass(histNames[icut][5 + histIdxOffset + 6].Data(), VarManager::fgValues);
                  }
		  }
                }
                if constexpr (TPairType == VarManager::kDecayToEE) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s", fTrackCuts[icut].Data()), VarManager::fgValues);
                  if (isAmbiExtra) {
                    fHistMan->FillHistClass(Form("PairsBarrelSEMM_ambiguousextra_%s", fTrackCuts[icut].Data()), VarManager::fgValues);
                  }
                }
              }
            }
            for (unsigned int iPairCut = 0; iPairCut < fPairCuts.size(); iPairCut++) {
              AnalysisCompositeCut cut = fPairCuts.at(iPairCut);
              if (!(cut.IsSelected(VarManager::fgValues))) // apply pair cuts
                continue;
              if (sign1 * sign2 < 0) {
                fHistMan->FillHistClass(histNames[ncuts + icut * ncuts + iPairCut][0].Data(), VarManager::fgValues);
              } else {
                if (sign1 > 0) {
                  fHistMan->FillHistClass(histNames[ncuts + icut * ncuts + iPairCut][1].Data(), VarManager::fgValues);
                } else {
                  fHistMan->FillHistClass(histNames[ncuts + icut * ncuts + iPairCut][2].Data(), VarManager::fgValues);
                }
              }
            } // end loop (pair cuts)
          }
        } // end loop (cuts)
      } // end loop over pairs of track associations
    } // end loop over events
  }

  template <int TPairType, uint32_t TEventFillMap, typename TAssoc1, typename TAssoc2, typename TTracks1, typename TTracks2>
  void runMixedPairing(TAssoc1 const& assocs1, TAssoc2 const& assocs2, TTracks1 const& /*tracks1*/, TTracks2 const& /*tracks2*/)
  {
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    int pairSign = 0;
    int ncuts = 0;
    uint32_t twoTrackFilter = static_cast<uint32_t>(0);
    for (auto& a1 : assocs1) {
      for (auto& a2 : assocs2) {
        if constexpr (TPairType == VarManager::kDecayToEE) {
          twoTrackFilter = a1.isBarrelSelected_raw() & a2.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isBarrelSelectedPrefilter_raw() & fTrackFilterMask;
          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }
          auto t1 = a1.template reducedtrack_as<TTracks1>();
          auto t2 = a2.template reducedtrack_as<TTracks2>();
          VarManager::FillPairME<TEventFillMap, TPairType>(t1, t2);
          if constexpr ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0) {
            VarManager::FillPairVn<TEventFillMap, TPairType>(t1, t2);
          }
          if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionQvect) > 0) {
            VarManager::FillPairVn<TEventFillMap, TPairType>(t1, t2);
          }
          pairSign = t1.sign() + t2.sign();
          ncuts = fNCutsBarrel;
        }
        if constexpr (TPairType == VarManager::kDecayToMuMu) {
          twoTrackFilter = a1.isMuonSelected_raw() & a2.isMuonSelected_raw() & fMuonFilterMask;
          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }
          auto t1 = a1.template reducedmuon_as<TTracks1>();
          auto t2 = a2.template reducedmuon_as<TTracks2>();
          if (t1.matchMCHTrackId() == t2.matchMCHTrackId() && t1.matchMCHTrackId() >= 0)
            continue;
          if (t1.matchMFTTrackId() == t2.matchMFTTrackId() && t1.matchMCHTrackId() >= 0)
            continue;
          VarManager::FillPairME<TEventFillMap, TPairType>(t1, t2);
          if constexpr ((TEventFillMap & VarManager::ObjTypes::ReducedEventQvector) > 0) {
            VarManager::FillPairVn<TEventFillMap, TPairType>(t1, t2);
          }
          pairSign = t1.sign() + t2.sign();
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
          ncuts = fNCutsMuon;
          histNames = fMuonHistNames;

          if (fConfigOptions.flatTables.value) {
            dimuonAllList(-999., -999., -999., -999.,
                          0, 0,
                          -999., -999., -999.,
                          VarManager::fgValues[VarManager::kMass],
                          false,
                          VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], t1.sign() + t2.sign(), VarManager::fgValues[VarManager::kVertexingChi2PCA],
                          VarManager::fgValues[VarManager::kVertexingTauz], VarManager::fgValues[VarManager::kVertexingTauzErr],
                          VarManager::fgValues[VarManager::kVertexingTauxy], VarManager::fgValues[VarManager::kVertexingTauxyErr],
                          VarManager::fgValues[VarManager::kCosPointingAngle],
                          t1.pt(), t1.eta(), t1.phi(), t1.sign(),
                          t2.pt(), t2.eta(), t2.phi(), t2.sign(),
                          t1.fwdDcaX(), t1.fwdDcaY(), t2.fwdDcaX(), t2.fwdDcaY(),
                          0., 0.,
                          t1.chi2MatchMCHMID(), t2.chi2MatchMCHMID(),
                          t1.chi2MatchMCHMFT(), t2.chi2MatchMCHMFT(),
                          t1.chi2(), t2.chi2(),
                          -999., -999., -999., -999.,
                          -999., -999., -999., -999.,
                          -999., -999., -999., -999.,
                          -999., -999., -999., -999.,
                          (twoTrackFilter & (static_cast<uint32_t>(1) << 28)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 29)), (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31)),
                          VarManager::fgValues[VarManager::kU2Q2], VarManager::fgValues[VarManager::kU3Q3],
                          VarManager::fgValues[VarManager::kR2EP_AB], VarManager::fgValues[VarManager::kR2SP_AB], VarManager::fgValues[VarManager::kCentFT0C],
                          VarManager::fgValues[VarManager::kCos2DeltaPhi], VarManager::fgValues[VarManager::kCos3DeltaPhi],
                          VarManager::fgValues[VarManager::kCORR2POI], VarManager::fgValues[VarManager::kCORR4POI], VarManager::fgValues[VarManager::kM01POI], VarManager::fgValues[VarManager::kM0111POI], VarManager::fgValues[VarManager::kMultDimuons],
                          VarManager::fgValues[VarManager::kVertexingPz], VarManager::fgValues[VarManager::kVertexingSV]);
          }
        }
        /*if constexpr (TPairType == VarManager::kElectronMuon) {
          twoTrackFilter = a1.isBarrelSelected_raw() & a1.isBarrelSelectedPrefilter_raw() & a2.isMuonSelected_raw() & fTrackFilterMask;
        }*/

        bool isAmbiInBunch = false;
        bool isAmbiOutOfBunch = false;
        bool isUnambiguous = false;
        for (int icut = 0; icut < ncuts; icut++) {
          if (!(twoTrackFilter & (static_cast<uint32_t>(1) << icut))) {
            continue; // cut not passed
          }
          isAmbiInBunch = (twoTrackFilter & (static_cast<uint32_t>(1) << 28)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 29));
          isAmbiOutOfBunch = (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31));
          isUnambiguous = !((twoTrackFilter & (static_cast<uint32_t>(1) << 28)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 29)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31)));
          if (pairSign == 0) {
            if constexpr (TPairType == VarManager::kDecayToMuMu) {
              fHistMan->FillHistClass(histNames[icut][3].Data(), VarManager::fgValues);
	      if (fConfigAmbiguousMuonHistograms) {
              if (isAmbiInBunch) {
                fHistMan->FillHistClass(histNames[icut][15].Data(), VarManager::fgValues);
              }
              if (isAmbiOutOfBunch) {
                fHistMan->FillHistClass(histNames[icut][18].Data(), VarManager::fgValues);
              }
              if (isUnambiguous) {
                fHistMan->FillHistClass(histNames[icut][21].Data(), VarManager::fgValues);
              }
	      }
            }
            if constexpr (TPairType == VarManager::kDecayToEE) {
              fHistMan->FillHistClass(Form("PairsBarrelMEPM_%s", fTrackCuts[icut].Data()), VarManager::fgValues);
            }
          } else {
            if (pairSign > 0) {
              if constexpr (TPairType == VarManager::kDecayToMuMu) {
                fHistMan->FillHistClass(histNames[icut][4].Data(), VarManager::fgValues);
		if (fConfigAmbiguousMuonHistograms) {
                if (isAmbiInBunch) {
                  fHistMan->FillHistClass(histNames[icut][16].Data(), VarManager::fgValues);
                }
                if (isAmbiOutOfBunch) {
                  fHistMan->FillHistClass(histNames[icut][19].Data(), VarManager::fgValues);
                }
                if (isUnambiguous) {
                  fHistMan->FillHistClass(histNames[icut][22].Data(), VarManager::fgValues);
                }
		}
              }
              if constexpr (TPairType == VarManager::kDecayToEE) {
                fHistMan->FillHistClass(Form("PairsBarrelMEPP_%s", fTrackCuts[icut].Data()), VarManager::fgValues);
              }
            } else {
              if constexpr (TPairType == VarManager::kDecayToMuMu) {
                fHistMan->FillHistClass(histNames[icut][5].Data(), VarManager::fgValues);
		if (fConfigAmbiguousMuonHistograms) {
                if (isAmbiInBunch) {
                  fHistMan->FillHistClass(histNames[icut][17].Data(), VarManager::fgValues);
                }
                if (isAmbiOutOfBunch) {
                  fHistMan->FillHistClass(histNames[icut][20].Data(), VarManager::fgValues);
                }
                if (isUnambiguous) {
                  fHistMan->FillHistClass(histNames[icut][23].Data(), VarManager::fgValues);
                }
		}
              }
              if constexpr (TPairType == VarManager::kDecayToEE) {
                fHistMan->FillHistClass(Form("PairsBarrelMEMM_%s", fTrackCuts[icut].Data()), VarManager::fgValues);
              }
            }
          }
        } // end for (cuts)
      } // end for (track2)
    } // end for (track1)
  }

  // barrel-barrel and muon-muon event mixing
  template <int TPairType, uint32_t TEventFillMap, typename TEvents, typename TAssocs, typename TTracks>
  void runSameSideMixing(TEvents& events, TAssocs const& assocs, TTracks const& tracks, Preslice<TAssocs>& preSlice)
  {
    events.bindExternalIndices(&assocs);
    int mixingDepth = fConfigMixingDepth.value;
    fAmbiguousPairs.clear();
    for (auto& [event1, event2] : selfCombinations(hashBin, mixingDepth, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<TEventFillMap>(event1, VarManager::fgValues);

      auto assocs1 = assocs.sliceBy(preSlice, event1.globalIndex());
      assocs1.bindExternalIndices(&events);

      auto assocs2 = assocs.sliceBy(preSlice, event2.globalIndex());
      assocs2.bindExternalIndices(&events);

      runMixedPairing<TPairType, TEventFillMap>(assocs1, assocs2, tracks, tracks);
    } // end event loop
  }

  void processAllSkimmed(MyEventsVtxCovSelected const& events,
                         soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs, MyBarrelTracksWithCovWithAmbiguities const& barrelTracks,
                         soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(events, muonAssocsPerCollision, muonAssocs, muons);
    // runSameEventPairing<true, VarManager::kElectronMuon, gkEventFillMap, gkTrackFillMap>(event, tracks, muons);
  }

  void processBarrelOnlySkimmed(MyEventsVtxCovSelected const& events,
                                soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processBarrelOnlySkimmedFlow(MyEventsVtxCovSelectedQvector const& events,
                                    soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                    MyBarrelTracksWithAmbiguities const& barrelTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithMultExtraWithQVector, gkTrackFillMap>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processBarrelOnlySkimmedNoCov(MyEventsSelected const& events,
                                     soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                     MyBarrelTracksWithAmbiguities const& barrelTracks)
  {
    runSameEventPairing<false, VarManager::kDecayToEE, gkEventFillMap, gkTrackFillMap>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processBarrelOnlySkimmedNoCovWithMultExtra(MyEventsMultExtraSelected const& events,
                                                  soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                                  MyBarrelTracksWithAmbiguities const& barrelTracks)
  {
    runSameEventPairing<false, VarManager::kDecayToEE, gkEventFillMapWithMultExtra, gkTrackFillMap>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processBarrelOnlyWithCollSkimmed(MyEventsVtxCovSelected const& events,
                                        soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                        MyBarrelTracksWithCovWithAmbiguitiesWithColl const& barrelTracks)
  {
    runSameEventPairing<true, VarManager::kDecayToEE, gkEventFillMapWithCov, gkTrackFillMapWithCovWithColl>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processBarrelOnlyWithQvectorCentrSkimmedNoCov(MyEventsQvectorCentrSelected const& events,
                                                     soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& barrelAssocs,
                                                     MyBarrelTracksWithAmbiguities const& barrelTracks)
  {
    runSameEventPairing<false, VarManager::kDecayToEE, gkEventFillMapWithQvectorCentr, gkTrackFillMap>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processMuonOnlySkimmed(MyEventsVtxCovSelected const& events,
                              soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons)
  {
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(events, muonAssocsPerCollision, muonAssocs, muons);
  }

  void processMuonOnlySkimmedMultExtra(MyEventsVtxCovSelectedMultExtra const& events,
                                       soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons)
  {
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithMultExtra, gkMuonFillMapWithCov>(events, muonAssocsPerCollision, muonAssocs, muons);
  }

  void processMuonOnlySkimmedFlow(MyEventsQvectorCentrSelected const& events,
                                  soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons)
  {
    runSameEventPairing<true, VarManager::kDecayToMuMu, gkEventFillMapWithMultExtraWithQVector, gkMuonFillMapWithCov>(events, muonAssocsPerCollision, muonAssocs, muons);
  }

  void processMixingAllSkimmed(soa::Filtered<MyEventsHashSelected>& events,
                               soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& trackAssocs, MyBarrelTracksWithCov const& tracks,
                               soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons)
  {
    runSameSideMixing<pairTypeEE, gkEventFillMap>(events, trackAssocs, tracks, trackAssocsPerCollision);
    runSameSideMixing<pairTypeMuMu, gkEventFillMap>(events, muonAssocs, muons, muonAssocsPerCollision);
  }

  void processMixingBarrelSkimmed(soa::Filtered<MyEventsHashSelected>& events,
                                  soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& trackAssocs, aod::ReducedTracks const& tracks)
  {
    runSameSideMixing<pairTypeEE, gkEventFillMap>(events, trackAssocs, tracks, trackAssocsPerCollision);
  }

  void processMixingBarrelSkimmedFlow(soa::Filtered<MyEventsVtxCovSelectedQvectorWithHash>& events,
                                      soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& trackAssocs, aod::ReducedTracks const& tracks)
  {
    runSameSideMixing<pairTypeEE, gkEventFillMapWithMultExtraWithQVector>(events, trackAssocs, tracks, trackAssocsPerCollision);
  }

  void processMixingBarrelWithQvectorCentrSkimmedNoCov(soa::Filtered<MyEventsHashSelectedQvectorCentr>& events,
                                                       soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts, aod::Prefilter> const& trackAssocs, MyBarrelTracksWithAmbiguities const& tracks)
  {
    runSameSideMixing<pairTypeEE, gkEventFillMapWithQvectorCentr>(events, trackAssocs, tracks, trackAssocsPerCollision);
  }

  void processMixingMuonSkimmed(soa::Filtered<MyEventsHashSelected>& events,
                                soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts> const& muonAssocs, MyMuonTracksWithCovWithAmbiguities const& muons)
  {
    runSameSideMixing<pairTypeMuMu, gkEventFillMap>(events, muonAssocs, muons, muonAssocsPerCollision);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisSameEventPairing, processAllSkimmed, "Run all types of pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlySkimmed, "Run barrel only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlyWithCollSkimmed, "Run barrel only pairing, with skimmed tracks and with collision information", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlySkimmedNoCov, "Run barrel only pairing (no covariances), with skimmed tracks and with collision information", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlySkimmedNoCovWithMultExtra, "Run barrel only pairing (no covariances), with skimmed tracks, with collision information, with MultsExtra", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlyWithQvectorCentrSkimmedNoCov, "Run barrel only pairing (no covariances), with skimmed tracks, with Qvector from central framework", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processBarrelOnlySkimmedFlow, "Run barrel only pairing, with skimmed tracks and with flow", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMuonOnlySkimmed, "Run muon only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMuonOnlySkimmedMultExtra, "Run muon only pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMuonOnlySkimmedFlow, "Run muon only pairing, with skimmed tracks and flow", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMixingAllSkimmed, "Run all types of mixed pairing, with skimmed tracks/muons", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMixingBarrelSkimmed, "Run barrel type mixing pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMixingBarrelSkimmedFlow, "Run barrel type mixing pairing, with flow, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMixingBarrelWithQvectorCentrSkimmedNoCov, "Run barrel type mixing pairing, with skimmed tracks and with Qvector from central framework", false);
  PROCESS_SWITCH(AnalysisSameEventPairing, processMixingMuonSkimmed, "Run muon type mixing pairing, with skimmed muons", false);
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
  Configurable<uint32_t> fConfigLegAFilterMask{"cfgLegAFilterMask", 0, "Filter mask corresponding to cuts in track-selection"};
  Configurable<uint32_t> fConfigLegBFilterMask{"cfgLegBFilterMask", 0, "Filter mask corresponding to cuts in track-selection"};
  Configurable<uint32_t> fConfigLegCFilterMask{"cfgLegCFilterMask", 0, "Filter mask corresponding to cuts in track-selection"};
  Configurable<std::string> fConfigCommonTrackCuts{"cfgCommonTrackCuts", "", "Comma separated list of cuts to be applied to all legs"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts"};
  Configurable<std::string> fConfigPairCutsJSON{"cfgPairCutsJSON", "", "Additional list of pair cuts in JSON format"};
  Configurable<bool> fConfigSkipAmbiguousIdCombinations{"cfgSkipAmbiguousIdCombinations", true, "Choose whether to skip pairs/triples which pass a stricter combination of cuts, e.g. KKPi triplets for D+ -> KPiPi"};

  Configurable<std::string> fConfigHistogramSubgroups{"cfgAsymmetricPairingHistogramsSubgroups", "barrel,vertexing", "Comma separated list of asymmetric-pairing histogram subgroups"};
  Configurable<bool> fConfigSameSignHistograms{"cfgSameSignHistograms", false, "Include same sign pair histograms for 2-prong decays"};
  Configurable<bool> fConfigAmbiguousHistograms{"cfgAmbiguousHistograms", false, "Include separate histograms for pairs/triplets with ambiguous tracks"};
  Configurable<bool> fConfigReflectedHistograms{"cfgReflectedHistograms", false, "Include separate histograms for pairs which are reflections of previously counted pairs"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};

  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> fConfigGRPMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> fConfigUseRemoteField{"cfgUseRemoteField", false, "Choose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};

  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fConfigUseAbsDCA{"cfgUseAbsDCA", false, "Use absolute DCA minimization instead of chi^2 minimization in secondary vertexing"};
  Configurable<bool> fConfigPropToPCA{"cfgPropToPCA", false, "Propagate tracks to secondary vertex"};
  Configurable<std::string> fConfigLutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  HistogramManager* fHistMan;

  std::vector<AnalysisCompositeCut*> fPairCuts;
  int fNPairHistPrefixes;

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
  int fNPairCuts;
  int fNCommonTrackCuts;
  // vectors for cut names and signal names, for easy access when calling FillHistogramList()
  std::vector<TString> fLegCutNames;
  std::vector<TString> fPairCutNames;
  std::vector<TString> fCommonCutNames;

  Preslice<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

  // Partitions for triplets and asymmetric pairs
  Partition<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> legACandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegAFilterMask) > static_cast<uint32_t>(0);
  Partition<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> legBCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegBFilterMask) > static_cast<uint32_t>(0);
  Partition<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> legCCandidateAssocs = (o2::aod::dqanalysisflags::isBarrelSelected & fConfigLegCFilterMask) > static_cast<uint32_t>(0);

  // Map to track how many times a pair of tracks has been encountered
  std::map<std::pair<int32_t, int32_t>, int8_t> fPairCount;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Get the leg cut filter maps
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
    if (objArrayLegs->GetEntries() == 0) {
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
        fConstructedLegAFilterMask |= (static_cast<uint32_t>(1) << legAIdx);
        fConstructedLegAFilterMasksMap[icut] |= static_cast<uint32_t>(1) << legAIdx;
      } else {
        LOGF(fatal, "Leg A cut %s was not calculated upstream. Check the config!", legs->At(0)->GetName());
        continue;
      }
      legBIdx = objArray->IndexOf(legs->At(1));
      if (legBIdx >= 0) {
        fConstructedLegBFilterMask |= (static_cast<uint32_t>(1) << legBIdx);
        fConstructedLegBFilterMasksMap[icut] |= static_cast<uint32_t>(1) << legBIdx;
      } else {
        LOGF(fatal, "Leg B cut %s was not calculated upstream. Check the config!", legs->At(1)->GetName());
        continue;
      }
      if (isThreeProng[icut]) {
        legCIdx = objArray->IndexOf(legs->At(2));
        if (legCIdx >= 0) {
          fConstructedLegCFilterMask |= (static_cast<uint32_t>(1) << legCIdx);
          fConstructedLegCFilterMasksMap[icut] |= static_cast<uint32_t>(1) << legCIdx;
        } else {
          LOGF(fatal, "Leg C cut %s was not calculated upstream. Check the config!", legs->At(2)->GetName());
          continue;
        }
      }
      // Leg cut config is fine, store the leg cut name in a vector
      fLegCutNames.push_back(legsStr);

      if (isThreeProng[icut]) {
        DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s", legsStr.Data()), fConfigHistogramSubgroups.value.data());
        if (fConfigAmbiguousHistograms.value) {
          DefineHistograms(fHistMan, Form("TripletsBarrelSE_ambiguous_%s", legsStr.Data()), fConfigHistogramSubgroups.value.data());
        }

        std::unique_ptr<TObjArray> objArrayCommon(commonNamesStr.Tokenize(","));
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
          DefineHistograms(fHistMan, Form("TripletsBarrelSE_%s_%s", legsStr.Data(), objArrayCommon->At(iCommonCut)->GetName()), fConfigHistogramSubgroups.value.data());
        }

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
        if (fConfigAmbiguousHistograms.value) {
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
  void runAsymmetricPairing(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& /*assocs*/, TTracks const& /*tracks*/)
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

        uint32_t twoTrackFilter = static_cast<uint32_t>(0);
        uint32_t twoTrackCommonFilter = static_cast<uint32_t>(0);
        uint32_t pairFilter = static_cast<uint32_t>(0);
        for (int icut = 0; icut < fNLegCuts; ++icut) {
          // Find leg pair definitions both candidates participate in
          if ((a1.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) && (a2.isBarrelSelected_raw() & fConstructedLegBFilterMasksMap[icut])) {
            twoTrackFilter |= (static_cast<uint32_t>(1) << icut);
            // If the supposed pion passes a kaon cut, this is a K+K-. Skip it.
            if (TPairType == VarManager::kDecayToKPi && fConfigSkipAmbiguousIdCombinations.value) {
              if (a2.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) {
                twoTrackFilter &= ~(static_cast<uint32_t>(1) << icut);
              }
            }
          }
        }
        if (!twoTrackFilter) {
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
          twoTrackFilter |= (static_cast<uint32_t>(1) << 30);
        }
        if (t2.barrelAmbiguityInBunch() > 1 || t2.barrelAmbiguityOutOfBunch() > 1) {
          twoTrackFilter |= (static_cast<uint32_t>(1) << 31);
        }

        VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);
        if constexpr (TTwoProngFitter) {
          VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, t1, t2, fConfigPropToPCA);
        }

        // Fill histograms
        bool isAmbi = false;
        for (int icut = 0; icut < fNLegCuts; icut++) {
          if (twoTrackFilter & (static_cast<uint32_t>(1) << icut)) {
            isAmbi = (twoTrackFilter & (static_cast<uint32_t>(1) << 30)) || (twoTrackFilter & (static_cast<uint32_t>(1) << 31));
            if (sign1 * sign2 < 0) {
              fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s", fLegCutNames[icut].Data()), VarManager::fgValues); // reconstructed, unmatched
              if (isAmbi && fConfigAmbiguousHistograms.value) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
              }
              if (isReflected && fConfigReflectedHistograms.value) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPM_reflected_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
              }
            } else if (fConfigSameSignHistograms.value) {
              if (sign1 > 0) {
                fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                if (isAmbi && fConfigAmbiguousHistograms.value) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
                if (isReflected && fConfigReflectedHistograms.value) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEPP_reflected_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
              } else {
                fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                if (isAmbi && fConfigAmbiguousHistograms.value) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
                }
                if (isReflected && fConfigReflectedHistograms) {
                  fHistMan->FillHistClass(Form("PairsBarrelSEMM_reflected_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
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
  void runThreeProng(TEvents const& events, Preslice<TTrackAssocs>& preslice, TTrackAssocs const& /*assocs*/, TTracks const& tracks, VarManager::PairCandidateType tripletType)
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
    uint32_t threeTrackFilter = static_cast<uint32_t>(0);
    uint32_t threeTrackCommonFilter = static_cast<uint32_t>(0);
    for (int icut = 0; icut < fNLegCuts; ++icut) {
      // Find out which leg cut combinations the triplet passes
      if ((a1.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) && (a2.isBarrelSelected_raw() & fConstructedLegBFilterMasksMap[icut]) && (a3.isBarrelSelected_raw() & fConstructedLegCFilterMasksMap[icut])) {
        threeTrackFilter |= (static_cast<uint32_t>(1) << icut);
        if (tripletType == VarManager::kTripleCandidateToPKPi && fConfigSkipAmbiguousIdCombinations.value) {
          // Check if the supposed pion passes as a proton or kaon, if so, skip this triplet. It is pKp or pKK.
          if ((a3.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) || (a3.isBarrelSelected_raw() & fConstructedLegBFilterMasksMap[icut])) {
            return;
          }
          // Check if the supposed kaon passes as a proton, if so, skip this triplet. It is ppPi.
          if (a2.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) {
            return;
          }
        }
        if (tripletType == VarManager::kTripleCandidateToKPiPi && fConfigSkipAmbiguousIdCombinations.value) {
          // Check if one of the supposed pions pass as a kaon, if so, skip this triplet. It is KKPi.
          if ((a2.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut]) || (a3.isBarrelSelected_raw() & fConstructedLegAFilterMasksMap[icut])) {
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

    VarManager::FillTriple(t1, t2, t3, VarManager::fgValues, tripletType);
    if constexpr (TThreeProngFitter) {
      VarManager::FillTripletVertexing<TEventFillMap, TTrackFillMap>(event, t1, t2, t3, tripletType);
    }

    // Fill histograms
    for (int icut = 0; icut < fNLegCuts; icut++) {
      if (threeTrackFilter & (static_cast<uint32_t>(1) << icut)) {
        fHistMan->FillHistClass(Form("TripletsBarrelSE_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
        if (fConfigAmbiguousHistograms.value && ((threeTrackFilter & (static_cast<uint32_t>(1) << 29)) || (threeTrackFilter & (static_cast<uint32_t>(1) << 30)) || (threeTrackFilter & (static_cast<uint32_t>(1) << 31)))) {
          fHistMan->FillHistClass(Form("TripletsBarrelSE_ambiguous_%s", fLegCutNames[icut].Data()), VarManager::fgValues);
        }
        for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
          if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
            fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data()), VarManager::fgValues);
          }
        } // end loop (common cuts)
        int iPairCut = 0;
        for (auto cut = fPairCuts.begin(); cut != fPairCuts.end(); cut++, iPairCut++) {
          if (!((*cut)->IsSelected(VarManager::fgValues))) // apply pair cuts
            continue;
          // Histograms with pair cuts
          fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
          // Histograms with pair cuts and common track cuts
          for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
            if (threeTrackCommonFilter & fCommonTrackCutFilterMasks[iCommonCut]) {
              fHistMan->FillHistClass(Form("TripletsBarrelSE_%s_%s_%s", fLegCutNames[icut].Data(), fCommonCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), VarManager::fgValues);
            }
          }
        } // end loop (pair cuts)
      }
    } // end loop (cuts)
  }

  void processKaonPionSkimmed(MyEventsVtxCovZdcFitSelected const& events,
                              soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                              MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runAsymmetricPairing<true, VarManager::kDecayToKPi, gkEventFillMapWithCovZdcFit, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processKaonPionSkimmedMultExtra(MyEventsVtxCovZdcFitSelectedMultExtra const& events,
                                       soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                                       MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runAsymmetricPairing<true, VarManager::kDecayToKPi, gkEventFillMapWithCovZdcFitMultExtra, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks);
  }

  void processKaonPionPionSkimmed(MyEventsVtxCovZdcFitSelected const& events,
                                  soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                                  MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runThreeProng<true, gkEventFillMapWithCovZdcFit, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, VarManager::kTripleCandidateToKPiPi);
  }

  void processKaonPionPionSkimmedMultExtra(MyEventsVtxCovZdcFitSelectedMultExtra const& events,
                                           soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                                           MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runThreeProng<true, gkEventFillMapWithCovZdcFitMultExtra, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, VarManager::kTripleCandidateToKPiPi);
  }

  void processProtonKaonPionSkimmed(MyEventsVtxCovZdcFitSelected const& events,
                                    soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts> const& barrelAssocs,
                                    MyBarrelTracksWithCovWithAmbiguities const& barrelTracks)
  {
    runThreeProng<true, gkEventFillMapWithCovZdcFit, gkTrackFillMapWithCov>(events, trackAssocsPerCollision, barrelAssocs, barrelTracks, VarManager::kTripleCandidateToPKPi);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionSkimmed, "Run kaon pion pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionPionSkimmed, "Run kaon pion pion triplets, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionSkimmedMultExtra, "Run kaon pion pairing, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processKaonPionPionSkimmedMultExtra, "Run kaon pion pion triplets, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processProtonKaonPionSkimmed, "Run proton kaon pion triplets, with skimmed tracks", false);
  PROCESS_SWITCH(AnalysisAsymmetricPairing, processDummy, "Dummy function, enabled only if none of the others are enabled", true);
};

// Combines dileptons with barrel or muon tracks for either resonance or correlation analyses
// Dileptons produced with all the selection cuts specified in the same-event pairing task are combined with the
//   tracks passing the fConfigTrackCut cut. The dileptons cuts from the same-event pairing task are auto-detected
struct AnalysisDileptonTrack {
  Produces<aod::BmesonCandidates> BmesonsTable;
  Produces<aod::JPsiMuonCandidates> DileptonTrackTable;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "kaonPID", "Comma separated list of cuts for the track to be correlated with the dileptons"};
  Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 2.8, "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 3.2, "High mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonLowpTCut{"cfgDileptonLowpTCut", 0.0, "Low pT cut for dileptons used in the triplet vertexing"};
  Configurable<float> fConfigDileptonHighpTCut{"cfgDileptonHighpTCut", 1E5, "High pT cut for dileptons used in the triplet vertexing"};
  Configurable<float> fConfigDileptonRapCutAbs{"cfgDileptonRapCutAbs", 1.0, "Rap cut for dileptons used in the triplet vertexing"};
  Configurable<float> fConfigDileptonLxyCut{"cfgDileptonLxyCut", 0.0, "Lxy cut for dileptons used in the triplet vertexing"};
  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};

  Configurable<std::string> fConfigHistogramSubgroups{"cfgDileptonTrackHistogramsSubgroups", "invmass,vertexing", "Comma separated list of dilepton-track histogram subgroups"};
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
  Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 5, "Event mixing pool depth"};

  Configurable<bool> fConfigUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<std::string> fConfigGRPmagPath{"cfgGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};

  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> fConfigGeoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<bool> fConfigUseRapcut{"cfgUseMCRapcut", false, "Use Rap cut for dileptons used in the triplet vertexing"};
  Configurable<bool> fConfigEnergycorrelator{"cfgEnergycorrelator", false, "Add some hist for energy correlator study"};

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.
  int fNCuts;      // number of dilepton leg cuts
  int fNLegCuts;
  int fNPairCuts; // number of pair cuts
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
  Filter eventFilter = aod::dqanalysisflags::isEventSelected > static_cast<uint8_t>(0);
  Filter dileptonFilter = aod::reducedpair::pt > fConfigDileptonLowpTCut&& aod::reducedpair::pt<fConfigDileptonHighpTCut && aod::reducedpair::mass> fConfigDileptonLowMass&& aod::reducedpair::mass<fConfigDileptonHighMass && aod::reducedpair::sign == 0 && aod::reducedpair::lxy> fConfigDileptonLxyCut;
  Filter filterBarrel = aod::dqanalysisflags::isBarrelSelected > static_cast<uint32_t>(0);
  Filter filterMuon = aod::dqanalysisflags::isMuonSelected > static_cast<uint32_t>(0);

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesHadron;
  HistogramManager* fHistMan;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> fHashBin;

  void init(o2::framework::InitContext& context)
  {
    bool isBarrel = context.mOptions.get<bool>("processBarrelSkimmed");
    bool isBarrelME = context.mOptions.get<bool>("processBarrelMixedEvent");
    bool isBarrelAsymmetric = context.mOptions.get<bool>("processDstarToD0Pi");
    bool isMuon = context.mOptions.get<bool>("processMuonSkimmed");
    bool isMuonME = context.mOptions.get<bool>("processMuonMixedEvent");

    // If the dummy process is enabled, skip the entire init
    if (context.mOptions.get<bool>("processDummy")) {
      if (isBarrel || isBarrelME || isBarrelAsymmetric || isMuon || isMuonME) {
        LOG(fatal) << "If processDummy is enabled, no other process functions should be enabled! Or switch off the processDummy!";
      }
      return;
    }

    fCurrentRun = 0;
    fTrackCutBitMap = 0;
    fValuesDilepton = new float[VarManager::kNVars];
    fValuesHadron = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(true);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // For each track/muon selection used to produce dileptons, create a separate histogram directory using the
    // name of the track/muon cut.
    if (isBarrel || isMuon || isBarrelAsymmetric) {
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
      // get the list of associated track cuts to be combined with the dileptons and
      //   check that these were computed upstream and create a bit mask
      TObjArray* cfgDileptonTrack_objArrayTrackCuts = nullptr;
      if (!fConfigTrackCuts.value.empty()) {
        cfgDileptonTrack_objArrayTrackCuts = TString(fConfigTrackCuts.value).Tokenize(",");
      } else {
        LOG(fatal) << " No track cuts specified! Check it out!";
      }
      for (int icut = 0; icut < cfgDileptonTrack_objArrayTrackCuts->GetEntries(); icut++) {
        if (!cfgTrackSelection_objArrayTrackCuts->FindObject(cfgDileptonTrack_objArrayTrackCuts->At(icut)->GetName())) {
          LOG(fatal) << "Specified track cut (" << cfgDileptonTrack_objArrayTrackCuts->At(icut)->GetName() << ") not found in the list of computed cuts by the single barrel / muon selection tasks";
        }
      }
      for (int icut = 0; icut < cfgTrackSelection_objArrayTrackCuts->GetEntries(); icut++) {
        if (cfgDileptonTrack_objArrayTrackCuts->FindObject(cfgTrackSelection_objArrayTrackCuts->At(icut)->GetName())) {
          fTrackCutBitMap |= (static_cast<uint32_t>(1) << icut);
        }
      }
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

          if (!cfgPairing_strCommonTrackCuts.IsNull()) {
            std::unique_ptr<TObjArray> objArrayCommon(cfgPairing_strCommonTrackCuts.Tokenize(","));
            for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
              DefineHistograms(fHistMan, Form("DileptonsSelected_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data()), "barrel,vertexing");
              DefineHistograms(fHistMan, Form("DileptonTrack_%s_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data(), fTrackCutNames[iCutTrack].Data()), fConfigHistogramSubgroups.value.data());
            }
          }

          if (fNPairCuts != 0) {

            for (int iPairCut = 0; iPairCut < fNPairCuts; ++iPairCut) {
              DefineHistograms(fHistMan, Form("DileptonsSelected_%s_%s", pairLegCutName.Data(), fPairCutNames[iPairCut].Data()), "barrel,vertexing");
              DefineHistograms(fHistMan, Form("DileptonTrack_%s_%s_%s", pairLegCutName.Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iCutTrack].Data()), fConfigHistogramSubgroups.value.data());

              if (!cfgPairing_strCommonTrackCuts.IsNull()) {
                std::unique_ptr<TObjArray> objArrayCommon(cfgPairing_strCommonTrackCuts.Tokenize(","));
                for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; ++iCommonCut) {
                  DefineHistograms(fHistMan, Form("DileptonsSelected_%s_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data()), "barrel,vertexing");
                  DefineHistograms(fHistMan, Form("DileptonTrack_%s_%s_%s_%s", pairLegCutName.Data(), fCommonPairCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iCutTrack].Data()), fConfigHistogramSubgroups.value.data());
                }
              }
            }
          }

          if (isBarrelME || isMuonME) {
            DefineHistograms(fHistMan, Form("DileptonTrackME_%s_%s", pairLegCutName.Data(), fTrackCutNames[iCutTrack].Data()), "dilepton-hadron-array-correlation"); // define ME histograms
            if (fConfigEnergycorrelator) {
              DefineHistograms(fHistMan, Form("DileptonTrackECME_%s_%s", pairLegCutName.Data(), fTrackCutNames[iCutTrack].Data()), "energy-correlator"); // define ME histograms
            }
          }
        } // end loop over track cuts to be combined with dileptons / di-tracks
      } // end loop over pair leg track cuts
    }

    dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // ad-hoc histograms via JSON

    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOG(info) << "Loading geometry from CCDB in dilepton-track task";
      fCCDB->get<TGeoManager>(fConfigGeoPath);
    }
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
  void runDileptonHadron(TEvent const& event, TTrackAssocs const& assocs, TTracks const& tracks, TDileptons const& dileptons)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesHadron);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);

    for (auto dilepton : dileptons) {
      // get full track info of tracks based on the index
      auto lepton1 = tracks.rawIteratorAt(dilepton.index0Id());
      auto lepton2 = tracks.rawIteratorAt(dilepton.index1Id());
      // Check that the dilepton has zero charge
      if (dilepton.sign() != 0) {
        continue;
      }
      // dilepton rap cut
      float rap = dilepton.rap();
      if (fConfigUseRapcut && abs(rap) > fConfigDileptonRapCutAbs)
        continue;

      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton);

      // loop over existing dilepton leg cuts (e.g. electron1, electron2, etc)
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
      } // end loop over single lepton selections

      // loop over hadrons
      for (auto& assoc : assocs) {

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
          // table to be written out for ML analysis
          BmesonsTable(event.runNumber(), event.globalIndex(), event.timestamp(), fValuesHadron[VarManager::kPairMass], dilepton.mass(), fValuesHadron[VarManager::kDeltaMass], fValuesHadron[VarManager::kPairPt], fValuesHadron[VarManager::kPairEta], fValuesHadron[VarManager::kPairPhi], fValuesHadron[VarManager::kPairRap],
                       fValuesHadron[VarManager::kVertexingLxy], fValuesHadron[VarManager::kVertexingLxyz], fValuesHadron[VarManager::kVertexingLz],
                       fValuesHadron[VarManager::kVertexingTauxy], fValuesHadron[VarManager::kVertexingTauz], fValuesHadron[VarManager::kCosPointingAngle],
                       fValuesHadron[VarManager::kVertexingChi2PCA],
                       track.globalIndex(), lepton1.globalIndex(), lepton2.globalIndex(),
                       track.tpcInnerParam(), track.eta(), dilepton.pt(), dilepton.eta(), lepton1.tpcInnerParam(), lepton1.eta(), lepton2.tpcInnerParam(), lepton2.eta(),
                       track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tpcNSigmaPr(), track.tofNSigmaKa(),
                       lepton1.tpcNSigmaEl(), lepton1.tpcNSigmaPi(), lepton1.tpcNSigmaPr(),
                       lepton2.tpcNSigmaEl(), lepton2.tpcNSigmaPi(), lepton2.tpcNSigmaPr(),
                       track.itsClusterMap(), lepton1.itsClusterMap(), lepton2.itsClusterMap(),
                       track.itsChi2NCl(), lepton1.itsChi2NCl(), lepton2.itsChi2NCl(),
                       track.tpcNClsFound(), lepton1.tpcNClsFound(), lepton2.tpcNClsFound(),
                       track.tpcChi2NCl(), lepton1.tpcChi2NCl(), lepton2.tpcChi2NCl(),
                       dilepton.filterMap_raw(), trackSelection);
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
          // Fill table for correlation analysis
          DileptonTrackTable(fValuesHadron[VarManager::kDeltaEta], fValuesHadron[VarManager::kDeltaPhi],
                             dilepton.mass(), dilepton.pt(), dilepton.eta(), track.pt(), track.eta(), track.phi(),
                             lepton1.pt(), lepton1.eta(), lepton1.phi(), lepton2.pt(), lepton2.eta(), lepton2.phi());
        }

        // Fill histograms for the triplets
        // loop over dilepton / ditrack cuts
        for (int icut = 0; icut < fNLegCuts; icut++) {

          if (!dilepton.filterMap_bit(icut)) {
            continue;
          }

          // loop over specified track cuts (the tracks to be combined with the dileptons)
          for (int iTrackCut = 0; iTrackCut < fNCuts; iTrackCut++) {

            if (!(trackSelection & (static_cast<uint32_t>(1) << iTrackCut))) {
              continue;
            }
            fHistMan->FillHistClass(Form("DileptonTrack_%s_%s", fLegCutNames[icut].Data(), fTrackCutNames[iTrackCut].Data()), fValuesHadron);

            if constexpr (TCandidateType == VarManager::kDstarToD0KPiPi) { // Dielectrons and Dimuons don't have the PairFilterMap column
              for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
                if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
                  fHistMan->FillHistClass(Form("DileptonTrack_%s_%s_%s", fLegCutNames[icut].Data(), fCommonPairCutNames[iCommonCut].Data(), fTrackCutNames[iTrackCut].Data()), fValuesHadron);
                }
              }
              for (int iPairCut = 0; iPairCut < fNPairCuts; iPairCut++) {
                if (dilepton.pairFilterMap_bit(iPairCut)) {
                  fHistMan->FillHistClass(Form("DileptonTrack_%s_%s_%s", fLegCutNames[icut].Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iTrackCut].Data()), fValuesHadron);
                  for (int iCommonCut = 0; iCommonCut < fNCommonTrackCuts; iCommonCut++) {
                    if (dilepton.commonFilterMap_bit(fCommonTrackCutMap[iCommonCut])) {
                      fHistMan->FillHistClass(Form("DileptonTrack_%s_%s_%s_%s", fLegCutNames[icut].Data(), fCommonPairCutNames[iCommonCut].Data(), fPairCutNames[iPairCut].Data(), fTrackCutNames[iTrackCut].Data()), fValuesHadron);
                    }
                  }
                }
              }
            }
          } // end loop over track cuts
        } // end loop over dilepton cuts
      }
    }
  }

  Preslice<aod::ReducedTracksAssoc> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<MyDielectronCandidates> dielectronsPerCollision = aod::reducedpair::reducedeventId;
  Preslice<MyDitrackCandidates> ditracksPerCollision = aod::reducedpair::reducedeventId;

  void processBarrelSkimmed(soa::Filtered<MyEventsVtxCovSelected> const& events,
                            soa::Filtered<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> const& assocs,
                            MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
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
      auto groupedDielectrons = dileptons.sliceBy(dielectronsPerCollision, event.globalIndex());
      runDileptonHadron<VarManager::kBtoJpsiEEK, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, groupedBarrelAssocs, tracks, groupedDielectrons);
    }
  }

  void processDstarToD0Pi(soa::Filtered<MyEventsVtxCovZdcFitSelected> const& events,
                          soa::Filtered<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> const& assocs,
                          MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDitrackCandidates> const& ditracks)
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
      runDileptonHadron<VarManager::kDstarToD0KPiPi, gkEventFillMapWithCovZdcFit, gkTrackFillMapWithCov>(event, groupedBarrelAssocs, tracks, groupedDitracks);
    }
  }

  Preslice<aod::ReducedMuonsAssoc> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<MyDimuonCandidates> dimuonsPerCollision = aod::reducedpair::reducedeventId;

  void processMuonSkimmed(soa::Filtered<MyEventsVtxCovSelected> const& events,
                          soa::Filtered<soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>> const& assocs,
                          MyMuonTracksWithCov const& tracks, soa::Filtered<MyDimuonCandidates> const& dileptons)
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
      auto groupedMuonAssocs = assocs.sliceBy(muonAssocsPerCollision, event.globalIndex());
      auto groupedDimuons = dileptons.sliceBy(dimuonsPerCollision, event.globalIndex());
      runDileptonHadron<VarManager::kBcToThreeMuons, gkEventFillMapWithCov, gkMuonFillMapWithCov>(event, groupedMuonAssocs, tracks, groupedDimuons);
    }
  }

  void processBarrelMixedEvent(soa::Filtered<MyEventsHashSelected>& events,
                               soa::Filtered<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> const& assocs,
                               MyBarrelTracksWithCov const&, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    if (events.size() == 0) {
      return;
    }
    events.bindExternalIndices(&dileptons);
    events.bindExternalIndices(&assocs);

    // loop over two event comibnations
    for (auto& [event1, event2] : selfCombinations(fHashBin, fConfigMixingDepth.value, -1, events, events)) {
      // fill event quantities
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMap>(event1, VarManager::fgValues);

      // get the dilepton slice for event1
      auto evDileptons = dileptons.sliceBy(dielectronsPerCollision, event1.globalIndex());
      evDileptons.bindExternalIndices(&events);

      // get the track associations slice for event2
      auto evAssocs = assocs.sliceBy(trackAssocsPerCollision, event2.globalIndex());
      evAssocs.bindExternalIndices(&events);

      // loop over associations
      for (auto& assoc : evAssocs) {

        // check that this track fulfills at least one of the specified cuts
        uint32_t trackSelection = (assoc.isBarrelSelected_raw() & fTrackCutBitMap);
        if (!trackSelection) {
          continue;
        }

        // get the track from this association
        auto track = assoc.template reducedtrack_as<MyBarrelTracksWithCov>();

        // loop over dileptons
        for (auto dilepton : evDileptons) {

          // compute dilepton - track quantities
          VarManager::FillDileptonHadron(dilepton, track, VarManager::fgValues);

          // loop over dilepton leg cuts and track cuts and fill histograms separately for each combination
          for (int icut = 0; icut < fNCuts; icut++) {
            if (!dilepton.filterMap_bit(icut)) {
              continue;
            }
            for (uint32_t iTrackCut = 0; iTrackCut < fTrackCutNames.size(); iTrackCut++) {
              if (trackSelection & (static_cast<uint32_t>(1) << iTrackCut)) {
                fHistMan->FillHistClass(Form("DileptonTrackME_%s_%s", fTrackCutNames[icut].Data(), fTrackCutNames[iTrackCut].Data()), VarManager::fgValues);
                if (fConfigEnergycorrelator) {
                  fHistMan->FillHistClass(Form("DileptonTrackECME_%s_%s", fTrackCutNames[icut].Data(), fTrackCutNames[iTrackCut].Data()), VarManager::fgValues);
                }
              }
            }
          }
        } // end for (dileptons)
      } // end for (assocs)
    } // end event loop
  }

  void processMuonMixedEvent(soa::Filtered<MyEventsHashSelected>& events,
                             soa::Filtered<soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>> const& assocs,
                             MyMuonTracksWithCov const&, soa::Filtered<MyDimuonCandidates> const& dileptons)
  {
    if (events.size() == 0) {
      return;
    }
    events.bindExternalIndices(&dileptons);
    events.bindExternalIndices(&assocs);

    for (auto& [event1, event2] : selfCombinations(fHashBin, fConfigMixingDepth.value, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMap>(event1, VarManager::fgValues);

      auto evDileptons = dileptons.sliceBy(dimuonsPerCollision, event1.globalIndex());
      evDileptons.bindExternalIndices(&events);

      auto evAssocs = assocs.sliceBy(muonAssocsPerCollision, event2.globalIndex());
      evAssocs.bindExternalIndices(&events);

      for (auto& assoc : evAssocs) {

        uint32_t muonSelection = assoc.isMuonSelected_raw() & fTrackCutBitMap;
        if (!muonSelection) {
          continue;
        }
        auto track = assoc.template reducedmuon_as<MyMuonTracksWithCov>();

        for (auto dilepton : evDileptons) {
          VarManager::FillDileptonHadron(dilepton, track, VarManager::fgValues);
          for (int icut = 0; icut < fNCuts; icut++) {
            if (!dilepton.filterMap_bit(icut)) {
              continue;
            }
            for (uint32_t iTrackCut = 0; iTrackCut < fTrackCutNames.size(); iTrackCut++) {
              if (muonSelection & (static_cast<uint32_t>(1) << iTrackCut)) {
                fHistMan->FillHistClass(Form("DileptonTrackME_%s_%s", fTrackCutNames[icut].Data(), fTrackCutNames[iTrackCut].Data()), VarManager::fgValues);
              }
            }
          }
        } // end for (dileptons)
      } // end for (assocs)
    } // end event loop
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisDileptonTrack, processBarrelSkimmed, "Run barrel dilepton-track pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processDstarToD0Pi, "Run barrel pairing of D0 daughters with pion candidate, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMuonSkimmed, "Run muon dilepton-track pairing, using skimmed data", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processBarrelMixedEvent, "Run barrel dilepton-hadron mixed event pairing", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processMuonMixedEvent, "Run muon dilepton-hadron mixed event pairing", false);
  PROCESS_SWITCH(AnalysisDileptonTrack, processDummy, "Dummy function", true);
};

struct AnalysisDileptonTrackTrack {
  OutputObj<THashList> fOutputList{"output"};

  Configurable<std::string> fConfigTrackCut1{"cfgTrackCut1", "pionPIDCut1", "track1 cut"}; // used for select the tracks from SelectedTracks
  Configurable<std::string> fConfigTrackCut2{"cfgTrackCut2", "pionPIDCut2", "track2 cut"}; // used for select the tracks from SelectedTracks
  Configurable<std::string> fConfigDileptonCut{"cfgDiLeptonCut", "pairJpsi2", "Dilepton cut"};
  Configurable<std::string> fConfigQuadrupletCuts{"cfgQuadrupletCuts", "pairX3872Cut1", "Comma separated list of Dilepton-Track-Track cut"};
  Configurable<std::string> fConfigAddDileptonHistogram{"cfgAddDileptonHistogram", "barrel", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddQuadrupletHistogram{"cfgAddQuadrupletHistogram", "xtojpsipipi", "Comma separated list of histograms"};

  Configurable<bool> fConfigSetupFourProngFitter{"cfgSetupFourProngFitter", false, "Use DCA for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fConfigUseKFVertexing{"cfgUseKFVertexing", false, "Use KF Particle for secondary vertex reconstruction (DCAFitter is used by default)"};
  Configurable<bool> fConfigUseRemoteField{"cfgUseRemoteField", false, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<std::string> fConfigGRPmagPath{"cfgGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};

  Produces<aod::DileptonTrackTrackCandidates> DileptonTrackTrackTable;

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.
  // uint32_t fTrackCutBitMap; // track cut bit mask to be used in the selection of tracks associated with dileptons
  // cut name setting
  TString fTrackCutName1;
  TString fTrackCutName2;
  bool fIsSameTrackCut = false;
  AnalysisCompositeCut fDileptonCut;
  std::vector<TString> fQuadrupletCutNames;
  std::vector<AnalysisCompositeCut> fQuadrupletCuts;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  Filter eventFilter = aod::dqanalysisflags::isEventSelected > static_cast<uint8_t>(0);
  Filter dileptonFilter = aod::reducedpair::sign == 0;
  Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > static_cast<uint32_t>(0);

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map

  // use some values array to avoid mixing up the quantities
  float* fValuesQuadruplet;
  HistogramManager* fHistMan;

  void init(o2::framework::InitContext& context)
  {
    // bool isBarrel = context.mOptions.get<bool>("processJpsiPiPi");

    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    fCurrentRun = 0;

    fValuesQuadruplet = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // define cuts
    fTrackCutName1 = fConfigTrackCut1.value;
    fTrackCutName2 = fConfigTrackCut2.value;
    if (fTrackCutName1 == fTrackCutName2) {
      fIsSameTrackCut = true;
    }
    TString configDileptonCutNamesStr = fConfigDileptonCut.value;
    fDileptonCut = *dqcuts::GetCompositeCut(configDileptonCutNamesStr.Data());
    TString configQuadruletCutNamesStr = fConfigQuadrupletCuts.value;
    std::unique_ptr<TObjArray> objArray(configQuadruletCutNamesStr.Tokenize(","));
    for (Int_t icut = 0; icut < objArray->GetEntries(); ++icut) {
      TString cutName = objArray->At(icut)->GetName();
      fQuadrupletCutNames.push_back(cutName);
      fQuadrupletCuts.push_back(*dqcuts::GetCompositeCut(cutName.Data()));
    }

    if (!context.mOptions.get<bool>("processDummy")) {
      DefineHistograms(fHistMan, Form("Pairs_%s", configDileptonCutNamesStr.Data()), fConfigAddDileptonHistogram.value.data());
      if (!configQuadruletCutNamesStr.IsNull()) {
        for (std::size_t icut = 0; icut < fQuadrupletCutNames.size(); ++icut) {
          if (fIsSameTrackCut) {
            DefineHistograms(fHistMan, Form("QuadrupletSEPM_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram.value.data());
          } else {
            DefineHistograms(fHistMan, Form("QuadrupletSEPM_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram.value.data());
            DefineHistograms(fHistMan, Form("QuadrupletSEMP_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram.value.data());
          }
          DefineHistograms(fHistMan, Form("QuadrupletSEPP_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram.value.data());
          DefineHistograms(fHistMan, Form("QuadrupletSEMM_%s", fQuadrupletCutNames[icut].Data()), fConfigAddQuadrupletHistogram.value.data());
        }
      }
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
        VarManager::SetupTwoProngKFParticle(magField);
        VarManager::SetupFourProngKFParticle(magField);
      } else if (fConfigSetupFourProngFitter.value) {
        VarManager::SetupTwoProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false);  // TODO: get these parameters from Configurables
        VarManager::SetupFourProngDCAFitter(magField, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false); // TODO: get these parameters from Configurables
      }
    } else {
      if (fConfigUseKFVertexing.value) {
        VarManager::SetupTwoProngKFParticle(fConfigMagField.value);
        VarManager::SetupFourProngKFParticle(fConfigMagField.value);
      } else if (fConfigSetupFourProngFitter.value) {
        LOGP(info, "Setting up DCA fitter for two and four prong candidates");
        VarManager::SetupTwoProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false);  // TODO: get these parameters from Configurables
        VarManager::SetupFourProngDCAFitter(fConfigMagField.value, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, false); // TODO: get these parameters from Configurables
      }
    }
  }

  // Template function to run pair - hadron combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks, typename TTrackAssocs, typename TDileptons>
  void runDileptonTrackTrack(TEvent const& event, TTrackAssocs const& assocs, TTracks const& tracks, TDileptons const& dileptons)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesQuadruplet);
    VarManager::FillEvent<TEventFillMap>(event, fValuesQuadruplet);
    // int indexOffset = -999;
    // std::vector<int> trackGlobalIndexes;

    for (auto dilepton : dileptons) {
      // get full track info of tracks based on the index

      int indexLepton1 = dilepton.index0Id();
      int indexLepton2 = dilepton.index1Id();
      auto lepton1 = tracks.rawIteratorAt(dilepton.index0Id());
      auto lepton2 = tracks.rawIteratorAt(dilepton.index1Id());
      // Check that the dilepton has zero charge
      if (dilepton.sign() != 0) {
        continue;
      }
      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesQuadruplet);
      // LOGP(info, "is dilepton selected: {}", fDileptonCut.IsSelected(fValuesQuadruplet));

      // apply the dilepton cut
      if (!fDileptonCut.IsSelected(fValuesQuadruplet))
        continue;

      fHistMan->FillHistClass(Form("Pairs_%s", fDileptonCut.GetName()), fValuesQuadruplet);

      // loop over hadrons pairs
      for (auto& [a1, a2] : o2::soa::combinations(assocs, assocs)) {
        uint32_t trackSelection = 0;
        if (fIsSameTrackCut) {
          trackSelection = ((a1.isBarrelSelected_raw() & (static_cast<uint32_t>(1) << 1)) || (a2.isBarrelSelected_raw() & (static_cast<uint32_t>(1) << 1)));
        } else {
          trackSelection = ((a1.isBarrelSelected_raw() & (static_cast<uint32_t>(1) << 1)) && (a2.isBarrelSelected_raw() & (static_cast<uint32_t>(1) << 2)));
        }
        // LOGP(info, "trackSelection: {}, a1: {}, a2: {}", trackSelection, a1.isBarrelSelected_raw(), a2.isBarrelSelected_raw());
        if (!trackSelection) {
          continue;
        }

        // get the track from this association
        auto track1 = a1.template reducedtrack_as<TTracks>();
        auto track2 = a2.template reducedtrack_as<TTracks>();
        // avoid self combinations
        if (track1.globalIndex() == indexLepton1 || track1.globalIndex() == indexLepton2 || track2.globalIndex() == indexLepton1 || track2.globalIndex() == indexLepton2) {
          continue;
        }

        // fill variables
        VarManager::FillDileptonTrackTrack<TCandidateType>(dilepton, track1, track2, fValuesQuadruplet);
        if (fConfigSetupFourProngFitter || fConfigUseKFVertexing) {
          // LOGP(info, "Using KF or DCA fitter for secondary vertexing");
          VarManager::FillDileptonTrackTrackVertexing<TCandidateType, TEventFillMap, TTrackFillMap>(event, lepton1, lepton2, track1, track2, fValuesQuadruplet);
        }

        int iCut = 0;
        uint32_t CutDecision = 0;
        for (auto cutname = fQuadrupletCutNames.begin(); cutname != fQuadrupletCutNames.end(); cutname++, iCut++) {
          // apply dilepton-track-track cut
          if (fQuadrupletCuts[iCut].IsSelected(fValuesQuadruplet)) {
            CutDecision |= (1 << iCut);
            if (fIsSameTrackCut) {
              if (track1.sign() * track2.sign() < 0) {
                fHistMan->FillHistClass(Form("QuadrupletSEPM_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
              }
            } else {
              if ((track1.sign() < 0) && (track2.sign() > 0)) {
                fHistMan->FillHistClass(Form("QuadrupletSEMP_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
              } else if ((track1.sign() > 0) && (track2.sign() < 0)) {
                fHistMan->FillHistClass(Form("QuadrupletSEPM_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
              }
            }
            if ((track1.sign() > 0) && (track2.sign() > 0)) {
              fHistMan->FillHistClass(Form("QuadrupletSEPP_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
            } else if ((track1.sign() < 0) && (track2.sign() < 0)) {
              fHistMan->FillHistClass(Form("QuadrupletSEMM_%s", fQuadrupletCutNames[iCut].Data()), fValuesQuadruplet);
            }
          }
        } // loop over dilepton-track-track cuts

        // fill table
        if (!CutDecision)
          continue;
        DileptonTrackTrackTable(fValuesQuadruplet[VarManager::kQuadMass], fValuesQuadruplet[VarManager::kQuadPt], fValuesQuadruplet[VarManager::kQuadEta], fValuesQuadruplet[VarManager::kQuadPhi], fValuesQuadruplet[VarManager::kRap],
                                fValuesQuadruplet[VarManager::kQ], fValuesQuadruplet[VarManager::kDeltaR1], fValuesQuadruplet[VarManager::kDeltaR2], fValuesQuadruplet[VarManager::kDeltaR],
                                dilepton.mass(), dilepton.pt(), dilepton.eta(), dilepton.phi(), dilepton.sign(),
                                lepton1.tpcNSigmaEl(), lepton1.tpcNSigmaPi(), lepton1.tpcNSigmaPr(), lepton1.tpcNClsFound(),
                                lepton2.tpcNSigmaEl(), lepton2.tpcNSigmaPi(), lepton2.tpcNSigmaPr(), lepton2.tpcNClsFound(),
                                fValuesQuadruplet[VarManager::kDitrackMass], fValuesQuadruplet[VarManager::kDitrackPt], track1.pt(), track2.pt(), track1.eta(), track2.eta(), track1.phi(), track2.phi(), track1.sign(), track2.sign(), track1.tpcNSigmaPi(), track2.tpcNSigmaPi(), track1.tpcNSigmaKa(), track2.tpcNSigmaKa(), track1.tpcNSigmaPr(), track1.tpcNSigmaPr(), track1.tpcNClsFound(), track2.tpcNClsFound(),
                                fValuesQuadruplet[VarManager::kKFMass], fValuesQuadruplet[VarManager::kVertexingProcCode], fValuesQuadruplet[VarManager::kVertexingChi2PCA], fValuesQuadruplet[VarManager::kCosPointingAngle], fValuesQuadruplet[VarManager::kKFDCAxyzBetweenProngs], fValuesQuadruplet[VarManager::kKFChi2OverNDFGeo],
                                fValuesQuadruplet[VarManager::kVertexingLz], fValuesQuadruplet[VarManager::kVertexingLxy], fValuesQuadruplet[VarManager::kVertexingLxyz], fValuesQuadruplet[VarManager::kVertexingTauz], fValuesQuadruplet[VarManager::kVertexingTauxy], fValuesQuadruplet[VarManager::kVertexingLzErr], fValuesQuadruplet[VarManager::kVertexingLxyzErr],
                                fValuesQuadruplet[VarManager::kVertexingTauzErr], fValuesQuadruplet[VarManager::kVertexingLzProjected], fValuesQuadruplet[VarManager::kVertexingLxyProjected], fValuesQuadruplet[VarManager::kVertexingLxyzProjected], fValuesQuadruplet[VarManager::kVertexingTauzProjected], fValuesQuadruplet[VarManager::kVertexingTauxyProjected]);
      }
    }
  }

  Preslice<aod::ReducedTracksAssoc> trackAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
  Preslice<MyDielectronCandidates> dielectronsPerCollision = aod::reducedpair::reducedeventId;
  Preslice<MyDitrackCandidates> ditracksPerCollision = aod::reducedpair::reducedeventId;

  void processJpsiPiPi(soa::Filtered<MyEventsVtxCovSelected> const& events,
                       soa::Filtered<soa::Join<aod::ReducedTracksAssoc, aod::BarrelTrackCuts>> const& assocs,
                       MyBarrelTracksWithCov const& tracks, soa::Filtered<MyDielectronCandidates> const& dileptons)
  {
    if (events.size() == 0) {
      return;
    }
    if (fCurrentRun != events.begin().runNumber()) { // start: runNumber
      initParamsFromCCDB(events.begin().timestamp());
      fCurrentRun = events.begin().runNumber();
    } // end: runNumber
    for (auto& event : events) {
      auto groupedBarrelAssocs = assocs.sliceBy(trackAssocsPerCollision, event.globalIndex());
      auto groupedDielectrons = dileptons.sliceBy(dielectronsPerCollision, event.globalIndex());
      runDileptonTrackTrack<VarManager::kXtoJpsiPiPi, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, groupedBarrelAssocs, tracks, groupedDielectrons);
    }
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(AnalysisDileptonTrackTrack, processJpsiPiPi, "Run barrel pairing of J/psi with pion candidate", false);
  PROCESS_SWITCH(AnalysisDileptonTrackTrack, processDummy, "Dummy function", true);
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
    adaptAnalysisTask<AnalysisDileptonTrack>(cfgc),
    adaptAnalysisTask<AnalysisDileptonTrackTrack>(cfgc)};
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

    if (classStr.Contains("Track") && !classStr.Contains("Pairs")) {
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
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
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
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-track", "dilepton-hadron-array-correlation");
    }

    if (classStr.Contains("DileptonTrackECME")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-track", "energy-correlator");
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

    if (classStr.Contains("Quadruplet")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-dihadron", histName);
    }
  } // end loop over histogram classes
}
