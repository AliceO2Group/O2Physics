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
/// @author  Kaare Endrup Iversen
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/StringHelpers.h>
#include <Framework/runDataProcessing.h>

#include <sys/types.h>

#include <RtypesCore.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

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
// DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 8);                      //! Event decision
// DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, int);
DECLARE_SOA_BITMAP_COLUMN(IsMuonSelected, isMuonSelected, 32);                       //! Muon track decisions (joinable to ReducedMuonsAssoc)
}

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTSA", dqanalysisflags::IsEventSelected);                                                            //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTSA", dqanalysisflags::IsMuonSelected);                                                       //!  joinable to ReducedMuonsAssoc
} 

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;

using MyPairCandidatesSelected = soa::Join<aod::Dimuons, aod::DimuonsExtra>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::MuonTrackCuts>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksSelectedWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;

struct DqJPsiMuonCorrelations {

  // Configurables for the dilepton and dilepton cuts
  Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 2., "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 4., "High mass cut for the dileptons used in analysis"};

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define the filter for events
  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;

  // Define the filter for the dileptons
  Filter dileptonFilter = aod::reducedpair::sign == 0;

  constexpr static uint32_t fgDimuonsFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::Pair;   // fill map

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesMuon;

  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    fValuesDilepton = new float[VarManager::kNVars];
    fValuesMuon = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();

    LOG(info) << "Initialized DqJPsiMuonCorrelations task" << std::endl;
  }

  // Template function to run pair - muon combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuons, typename TDileptons>
  void runDileptonMuon(TEvent const& event, TMuons const& muons, TDileptons const& dileptons)
  {
    // VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesMuon);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesMuon);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);

    if (dileptons.size() > 0) {
      LOG(info) << "Dileptons in event: " << dileptons.size() << std::endl;
    } else {
      LOG(info) << "No dileptons in event" << std::endl;
    }
  }

  void processSkimmedDimuon(MyEventsSelected::iterator const& event, MyMuonTracksSelected const& muons, soa::Filtered<MyPairCandidatesSelected> const& dileptons)
  {
    LOG(info) << "Processing event " << event.globalIndex() << " with " << muons.size() << " muons and " << dileptons.size() << " dileptons" << std::endl;
    runDileptonMuon<VarManager::kDecayToMuMu, gkEventFillMap, gkMuonFillMap>(event, muons, dileptons);
  }
  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(DqJPsiMuonCorrelations, processSkimmedDimuon, "Run dilepton-muon pairing, using skimmed data", false);
  PROCESS_SWITCH(DqJPsiMuonCorrelations, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DqJPsiMuonCorrelations>(cfgc)};
}
