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
// Contact: luca.micheletti@to.infn.it, INFN
//
#include <iostream>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;
using namespace o2::aod;

// Some definitions
namespace o2::aod
{

namespace dqanalysisflags
{
// TODO: the barrel amd muon selection columns are bit maps so unsigned types should be used, however, for now this is not supported in Filter expressions
// TODO: For now in the tasks we just statically convert from unsigned int to int, which should be fine as long as we do
//      not use a large number of bits (>=30)
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, int);
DECLARE_SOA_COLUMN(IsPrefilterVetoed, isPrefilterVetoed, int);
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASH", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected, dqanalysisflags::IsBarrelSelectedPrefilter);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsPrefilterVetoed);
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;
using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector>;
using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
using MyEventsHashSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithPrefilter = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::Prefilter>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyPairCandidatesSelected = soa::Join<aod::Dileptons, aod::DileptonsExtra>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::MuonTrackCuts>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksSelectedWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonTrackCuts>;
using MyD0CandidatesSelected = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;

// Global function used to define needed histogram classes
void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar); // defines histograms for all tasks

struct taskJpsiHf {
  //
  // This task combines dilepton candidates with a heavy flavors and could be used for example
  //

  OutputObj<THashList> fOutputList{"output"};
  // TODO: For now this is only used to determine the position in the filter bit map for the hadron cut
  Configurable<string> fConfigTrackCuts{"cfgLeptonCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  // comment: add list of subgroups (must define subgroups under )
  Configurable<std::string> fConfigAddDileptonHadHistogram{"cfgAddDileptonHadHistogram", "", "Comma separated list of histograms"};

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for HF flagged candidates"};
  Configurable<int> selectionTopol{"selectionTopol", 1, "Selection Flag for topologically selected candidates"};
  Configurable<int> selectionCand{"selectionCand", 1, "Selection Flag for conj. topol. selected candidates"};
  Configurable<int> selectionPid{"selectionPid", 1, "Selection Flag for reco PID candidates"};

  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;
  Filter dileptonFilter = aod::reducedpair::mass > 2.92f && aod::reducedpair::mass < 3.16f && aod::reducedpair::sign == 0;

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map
  constexpr static uint32_t fgDmesonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair;   // fill map

  Preslice<MyD0CandidatesSelected> perCollision = aod::track::collisionId;

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesHadron;
  float* fValuesDmeson;
  HistogramManager* fHistMan;

  // NOTE: the barrel track filter is shared between the filters for dilepton electron candidates (first n-bits)
  //       and the associated hadrons (n+1 bit) --> see the barrel track selection task
  //      The current condition should be replaced when bitwise operators will become available in Filter expressions

  int fNHadronCutBit;

  void init(o2::framework::InitContext& context)
  {
    fValuesDilepton = new float[VarManager::kNVars];
    fValuesHadron = new float[VarManager::kNVars];
    fValuesDmeson = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // TODO: Create separate histogram directories for each selection used in the creation of the dileptons
    // TODO: Implement possibly multiple selections for the associated track ?
    if (context.mOptions.get<bool>("processSkimmedJpsiD0")) {
      DefineHistograms(fHistMan, "DimuonsSelected;DmesonsSelected;DileptonHadronInvMass", fConfigAddDileptonHadHistogram); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    TString configCutNamesStr = fConfigTrackCuts.value;
    if (!configCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(configCutNamesStr.Tokenize(","));
      fNHadronCutBit = objArray->GetEntries() - 1;
    } else {
      fNHadronCutBit = 0;
    }
  }

  // Template function to run pair - hadron combinations
  template <uint32_t TEventFillMap, typename TEvent, typename TDqTrack, typename THfTrack>
  void runDileptonDmeson(TEvent const& event, TDqTrack const& dileptons, THfTrack const& dmesons)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDmeson);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDmeson);

    // std::cout << "---> Event global index: " << event.globalIndex() << std::endl;

    // loop over D mesons
    for (auto& dmeson : dmesons) {
      // std::cout << "------> D meson collision index: " << dmeson.collisionId() << std::endl;
      if (!(dmeson.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(yD0(dmeson)) > yCandMax) {
        continue;
      }

      if (dmeson.isSelD0() >= selectionFlagD0) {
        VarManager::FillHadron(dmeson, fValuesDmeson, invMassD0ToPiK(dmeson));
        fHistMan->FillHistClass("DmesonsSelected", fValuesDmeson);
      }

      if (dmeson.isSelD0bar() >= selectionFlagD0bar) {
        VarManager::FillHadron(dmeson, fValuesDmeson, invMassD0barToKPi(dmeson));
        fHistMan->FillHistClass("DmesonsSelected", fValuesDmeson);
      }
    }

    // loop over dileptons
    for (auto dilepton : dileptons) {
      // std::cout << "------> Dilepton collision index: " << dilepton.reducedeventId() << std::endl;
      // Check that the dilepton has zero charge
      if (dilepton.sign() != 0) {
        continue;
      }

      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton);
      fHistMan->FillHistClass("DimuonsSelected", fValuesDilepton);

      // loop over D mesons
      for (auto& dmeson : dmesons) {
        if (dmeson.isSelD0() >= selectionFlagD0) {
          VarManager::FillDileptonHadron(dilepton, dmeson, fValuesDmeson, invMassD0ToPiK(dmeson));
          fHistMan->FillHistClass("DileptonHadronInvMass", fValuesDmeson);
        }
        if (dmeson.isSelD0bar() >= selectionFlagD0bar) {
          VarManager::FillDileptonHadron(dilepton, dmeson, fValuesDmeson, invMassD0barToKPi(dmeson));
          fHistMan->FillHistClass("DileptonHadronInvMass", fValuesDmeson);
        }
      }
    }
    // std::cout << "------------------------" << std::endl;
  }

  void processSkimmedJpsiD0(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyPairCandidatesSelected const& dileptons, MyD0CandidatesSelected const& dmesons)
  {
    auto groupedDmesonCandidates = dmesons.sliceBy(perCollision, event.globalIndex());
    runDileptonDmeson<gkEventFillMapWithCov>(event, dileptons, groupedDmesonCandidates);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(taskJpsiHf, processSkimmedJpsiD0, "Run dilepton-D meson pairing, using skimmed data", false);
  PROCESS_SWITCH(taskJpsiHf, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taskJpsiHf>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar)
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

    TString histName = configVar.value;
    // NOTE: The level of detail for histogramming can be controlled via configurables

    if (classStr.Contains("DimuonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "dimuon");
    }

    if (classStr.Contains("DmesonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "dmeson");
    }

    if (classStr.Contains("DileptonHadronInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-mass");
    }
  } // end loop over histogram classes
}
