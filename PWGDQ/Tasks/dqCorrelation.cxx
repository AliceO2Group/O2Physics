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
//
#include <iostream>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>
#include <TRandom3.h>
#include <cmath>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
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

#include "PWGCF/GenericFramework/Core/GFWPowerArray.h"
#include "TProfile.h"

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
// using MyEventsVtxCovSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts, aod::ReducedEventsQvector>;
// using MyEventsQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsQvector>;
// using MyEventsHashSelectedQvector = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes, aod::ReducedEventsQvector>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithPrefilter = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::Prefilter>;
using MyBarrelTracksSelectedWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;
using MyBarrelTracksSelectedWithColl = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts, aod::ReducedTracksBarrelInfo>;
// using MyPairCandidatesSelected = soa::Join<aod::Dileptons, aod::DileptonsExtra>;
using MyPairCandidatesSelected = aod::Dileptons;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo, aod::MuonTrackCuts>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo, aod::ReducedMuonsCov>;
using MyMuonTracksSelectedWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo, aod::ReducedMuonsCov, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
// constexpr static uint32_t gkEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
//  constexpr static uint32_t gkEventFillMapWithQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventQvector;
//  constexpr static uint32_t gkEventFillMapWithCovQvector = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov | VarManager::ObjTypes::ReducedEventQvector;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
// constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
//  constexpr static uint32_t gkTrackFillMapWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;

constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
//  constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;

// constexpr static int pairTypeEE = VarManager::kDecayToEE;
// constexpr static int pairTypeMuMu = VarManager::kDecayToMuMu;
// constexpr static int pairTypeEMu = VarManager::kElectronMuon;

void DefineHistograms(HistogramManager* histMan, TString histClasses);
void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar);

// Fill the FlowContainer
void FillFC(GFW* fGFW, OutputObj<FlowContainer> fFC, const GFW::CorrConfig& corrconf, const double& cent, const double& rndm, bool fillflag);

struct DqCumulantFlow {

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -6.0, 1.5}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {100, 0, 20}, "pt axis for histograms"};
  ConfigurableAxis axisMass{"axisMass", {40, 2, 4}, "mass axis for histograms"};
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Configurable<uint8_t> fConfigNPow{"cfgNPow", 0, "Power of weights for Q vector"};
  // Configurables for the reference flow
  Configurable<string> fConfigTrackCuts{"cfgLeptonCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  Configurable<float> fConfigCutPtMin{"cfgCutPtMin", 1.0f, "Minimal pT for tracks"};
  Configurable<float> fConfigCutPtMax{"cfgCutPtMax", 12.0f, "Maximal pT for tracks"};
  Configurable<float> fConfigCutEtaMin{"cfgCutEtaMin", -0.8f, "Eta min range for tracks"};
  Configurable<float> fConfigCutEtaMax{"cfgCutEtaMax", 0.8f, "Eta max range for tracks"};
  Configurable<float> fConfigEtaLimitMin{"cfgEtaLimitMin", -0.4f, "Eta gap min separation, only if using subEvents"};
  Configurable<float> fConfigEtaLimitMax{"cfgEtaLimitMax", 0.4f, "Eta gap max separation, only if using subEvents"};

  /////////////////////

  Configurable<float> fConfigMuonCutEtaMin{"cfgCutEtaMinMuon", -4.f, "Eta min range for muons"};
  Configurable<float> fConfigMuonCutEtaMax{"cfgCutEtaMaxMuon", -2.5f, "Eta max range for muons"};

  Configurable<float> fConfigMultiplicityMin{"cfgMultiplicityMin", 3500.0f, "Multiplicity for 30%"};
  Configurable<float> fConfigMultiplicityMax{"cfgMultiplicityMax", 2000.0f, "Multiplicity for 60%"};

  /////////////////////

  // Configurables for the dilepton and dilepton cuts
  Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 2., "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 4., "High mass cut for the dileptons used in analysis"};

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // // // CCDB access
  // Service<ccdb::BasicCCDBManager> ccdb;
  // // CCDB access for the efficiencies and acceptances
  // Configurable<std::string> fConfigEfficiency{"ccdb-path-efficiency", "Users/r/rcaron/efficiency", "CCDB path to efficiency object"};
  // Configurable<std::string> fConfigAcceptance{"ccdb-path-acceptance", "", "CCDB path to acceptance or GFWWeights object"};
  // Configurable<std::string> fConfigURL{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // Configurables for histograms
  Configurable<std::string> fConfigAddDileptonHadHistogram{"cfgAddDileptonHadHistogram", "", "Comma separated list of histograms"};

  // Configurables for FlowContainer (e.g charged particles pt-differential v2{2}, v2{3}, ...)
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1}, "multiplicity / centrality axis for histograms"};

  ////////////////////////////////////////////////////////////////////////////////
  // Define the filter for events
  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;

  // Define the filter for barrel tracks and forward tracks
  Filter trackFilter = (nabs(aod::track::eta) <= fConfigCutEtaMax) && (aod::track::pt > fConfigCutPtMin) && (aod::track::pt < fConfigCutPtMax);
  Filter fwdFilter = (nabs(aod::fwdtrack::eta) <= fConfigMuonCutEtaMax) && (nabs(aod::fwdtrack::eta) >= fConfigMuonCutEtaMin);
  // Filter fwdFilter = (aod::fwdtrack::eta < -2.f) && (aod::fwdtrack::eta > -4.5f);

  // Define the filter for the dileptons
  // Filter dileptonFilter = aod::reducedpair::mass > fConfigDileptonLowMass.value&& aod::reducedpair::mass < fConfigDileptonHighMass.value&& aod::reducedpair::sign == 0;
  Filter dileptonFilter = aod::reducedpair::sign == 0;
  // Filter filterBarrelTrackSelected = aod::dqanalysisflags::isBarrelSelected > 0;
  // Filter filterMuonSelected = aod::dqanalysisflags::isMuonlSelected > 0;

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair; // fill map
  constexpr static uint32_t fgDimuonsFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::Pair;   // fill map

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesHadron;

  // NOTE: the barrel track filter is shared between the filters for dilepton electron candidates (first n-bits)
  //       and the associated hadrons (n+1 bit) --> see the barrel track selection task
  //      The current condition should be replaced when bitwise operators will become available in Filter expressions
  int fNHadronCutBit;

  NoBinningPolicy<aod::dqanalysisflags::MixingHash> hashBin;

  // Histograms used for optionnal efficiency and non-uniform acceptance corrections
  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance = nullptr;
  } cfg;

  // Define output
  HistogramManager* fHistMan = nullptr;
  OutputObj<FlowContainer> fFC{FlowContainer("flowContainer")};
  OutputObj<THashList> fOutputList{"outputQA"};
  HistogramRegistry registry{"registry"};

  // Define global variables for generic framework
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;

  void init(o2::framework::InitContext& context)
  {

    /////////////////////////////////////////////////////////////////////////////////////

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Add some output objects to the histogram registry
    registry.add("hPhi_barrel", "", {HistType::kTH1D, {axisPhi}});
    registry.add("hEta_barrel_muon", "", {HistType::kTH1D, {axisEta}});
    registry.add("hPt_muon", "", {HistType::kTH1D, {axisPt}});
    registry.add("hmass_muon", "", {HistType::kTH1D, {{300, 0.0, 5.0}}});
    registry.add("hmass_muon_cent", "", {HistType::kTH1D, {{300, 0.0, 5.0}}});
    registry.add("hmass_muon_pt", "", {HistType::kTH1D, {{300, 2.0, 5.0}}});

    // Add histograms for other correlator configurations
    registry.add("c22_mass_muon", "c22_mass_muon", {HistType::kTProfile, {axisMass}});
    registry.add("c22_mass_muon_cent", "c22_mass_muon_cent", {HistType::kTProfile, {axisMass}});
    registry.add("c24_mass_muon", "c24_mass_muon", {HistType::kTProfile, {axisMass}});
    registry.add("c24_mass_muon_cent", "c24_mass_muon_cent", {HistType::kTProfile, {axisMass}});
    registry.add("c22_cent_muon", "c22_cent_muon", {HistType::kTProfile, {{100, 0, 100}}});
    registry.add("hCent_muon", "", {HistType::kTH1D, {{10, 0.0, 100.0}}});
    registry.add("hCent_muon2", "", {HistType::kTH1D, {{10, 0.0, 100.0}}});

    registry.add("dimuon_mass", "", {HistType::kTH1D, {{300, 0.0, 5.0}}});
    /////////////////////////////////////////////////////////////////////////////////////

    fValuesDilepton = new float[VarManager::kNVars];
    fValuesHadron = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // TODO: Create separate histogram directories for each selection used in the creation of the dileptons
    // TODO: Implement possibly multiple selections for the associated track ?
    if (context.mOptions.get<bool>("processSkimmed") || context.mOptions.get<bool>("processSkimmedDimuon")) {
      DefineHistograms(fHistMan, "DileptonsSelected", fConfigAddDileptonHadHistogram); // define all histograms
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
    TString configCutNamesStr = fConfigTrackCuts.value;
    if (!configCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(configCutNamesStr.Tokenize(","));
      fNHadronCutBit = objArray->GetEntries() - 1;
    } else {
      fNHadronCutBit = 0;
    }

    // ccdb->setURL(fConfigURL.value);
    // ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // ccdb->setCreatedNotAfter(fConfigNoLaterThan.value);

    VarManager::SetDefaultVarNames();

    // // Global effiencies
    // if (fConfigEfficiency.value.empty() == false) {
    //   cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(fConfigEfficiency.value, fConfigNoLaterThan.value);
    // }

    // // Global acceptance or GFWeights to correct for NUA in the track loop
    // if (fConfigAcceptance.value.empty() == false) {
    //   cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(fConfigAcceptance.value, fConfigNoLaterThan.value);
    // }

    // Reference flow
    TObjArray* oba = new TObjArray();
    // REF name declaration
    oba->Add(new TNamed("ChGap22", "ChGap22"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap24", "ChGap24"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChFull22", "ChFull22")); // no-gap case
    oba->Add(new TNamed("ChFull24", "ChFull24")); // no-gap case
    oba->Add(new TNamed("ChGap32", "ChGap32"));   // gap-case
    // POIs name declaration
    oba->Add(new TNamed("CdhFull22", "CdhFull22"));               // no-gap case
    oba->Add(new TNamed("CdhFull24Dimuons", "CdhFull24Dimuons")); // no-gap case

    fFC->SetName("flowContainer");
    fFC->Initialize(oba, axisMultiplicity, 10);
    delete oba;

    // Define regions of positive and negative eta in order to create gaps
    fGFW->AddRegion("refN", fConfigCutEtaMin, fConfigEtaLimitMin, 1, 1);
    fGFW->AddRegion("refP", fConfigEtaLimitMax, fConfigCutEtaMax, 1, 1);
    fGFW->AddRegion("full", fConfigCutEtaMin, fConfigCutEtaMax, 1, 2);
    fGFW->AddRegion("dilepton", fConfigCutEtaMin, fConfigCutEtaMax, 1, 4);
    fGFW->AddRegion("dimuon", fConfigMuonCutEtaMin, fConfigMuonCutEtaMax, 1, 5);

    // Defined the different charged particle correlations
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2} refN {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2 2} refN {-2 -2}", "ChGap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {3} refN {-3}", "ChGap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("dilepton {2} full {-2}", "CdhFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("dimuon {2} full {-2}", "CdhFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("dimuon {2} full {2 -2 -2}", "CdhFull24Dimuons", kFALSE));

    fGFW->CreateRegions();
  }

  template <char... chars>
  void FillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& axis)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        registry.fill(tarName, axis, val, dnx);
      return;
    }
    return;
  }

  // Template function to run pair - hadron combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runDileptonHadron(TEvent const& event, TTracks const& tracks, soa::Filtered<MyPairCandidatesSelected> const& dileptons)
  {
    fGFW->Clear();
    VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesHadron);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);

    // Set the global index offset to find the proper lepton
    // TO DO: remove it once the issue with lepton index is solved
    // int indexOffset = -999;
    std::vector<int> trackGlobalIndexes;

    float weff = 1.0, wacc = 1.0;

    if (dileptons.size() > 0) {
      for (auto track : tracks) {
        trackGlobalIndexes.push_back(track.globalIndex());
        // std::cout << track.index() << " " << track.globalIndex() << std::endl;
      }

      for (auto& track : tracks) {
        if (cfg.mEfficiency) {
          weff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(track.pt()));
        } else {
          weff = 1.0;
        }
        if (weff == 0) {
          continue;
        }
        weff = 1. / weff;
        if (cfg.mAcceptance) {
          wacc = cfg.mAcceptance->GetNUA(track.phi(), track.eta(), event.posZ());
        } else {
          wacc = 1.0;
        }

        if (track.eta() >= -0.9 && track.eta() <= 0.9) {

          // Fill the GFW for each track to compute Q vector and correction using weights
          fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, 3); // using default values for ptin = 0 and mask = 3

          registry.fill(HIST("hPhi_barrel"), track.phi());
          registry.fill(HIST("hEta_barrel_muon"), track.eta());
        }
      }

      VarManager::fgValues[VarManager::kCentFT0C];

      for (auto dilepton : dileptons) {
        registry.fill(HIST("dimuon_mass"), dilepton.mass());

        // VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton);
        VarManager::FillTrack<fgDimuonsFillMap>(dilepton, fValuesDilepton);
        fHistMan->FillHistClass("DileptonsSelected", fValuesDilepton);

        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::ReducedTrack)) {
          fGFW->Fill(dilepton.eta(), 0, dilepton.phi(), wacc * weff, 4);
        } else if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::ReducedMuon)) {

          if (dilepton.eta() >= -4.0 && dilepton.eta() <= -2.5) {
            fGFW->Fill(dilepton.eta(), 0, dilepton.phi(), wacc * weff, 5);
            registry.fill(HIST("hEta_barrel_muon"), dilepton.eta());
            registry.fill(HIST("hPt_muon"), dilepton.pt());
            registry.fill(HIST("hmass_muon"), dilepton.mass());
            registry.fill(HIST("hCent_muon"), fValuesDilepton[VarManager::kCentFT0C]);

            FillProfile(corrconfigs.at(6), HIST("c22_mass_muon"), dilepton.mass());
            FillProfile(corrconfigs.at(7), HIST("c24_mass_muon"), dilepton.mass());

            if ((fValuesDilepton[VarManager::kCentFT0C] < 50) && ((fValuesDilepton[VarManager::kCentFT0C] > 30))) {

              FillProfile(corrconfigs.at(6), HIST("c22_cent_muon"), fValuesDilepton[VarManager::kCentFT0C]);

              FillProfile(corrconfigs.at(6), HIST("c22_mass_muon_cent"), dilepton.mass());

              FillProfile(corrconfigs.at(7), HIST("c24_mass_muon_cent"), dilepton.mass());
              registry.fill(HIST("hmass_muon_cent"), dilepton.mass());

              registry.fill(HIST("hCent_muon2"), fValuesDilepton[VarManager::kCentFT0C]);
            }

            if ((dilepton.pt() > 5) && (dilepton.pt() < 1000)) {

              registry.fill(HIST("hmass_muon_pt"), dilepton.mass());
            }
          }
        }
      }
    }

    double l_Random = fRndm->Rndm(); // used only to compute correlators
    bool fillFlag = kFALSE;          // could be used later
    for (uint64_t l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      // FillFC(corrconfigs.at(l_ind), VarManager::fgValues[VarManager::kCentFT0C], l_Random, fillFlag);
      // Because of the absence of the multipilcity, the multiplicity is set to 5 temporarily.
      // FillFC(fGFW, fFC, corrconfigs.at(l_ind), 5, l_Random, fillFlag);
      FillFC(fGFW, fFC, corrconfigs.at(l_ind), VarManager::fgValues[VarManager::kCentFT0C], l_Random, fillFlag);
    }
  }

  void processSkimmed(MyEventsSelected::iterator const& event, MyBarrelTracksSelected const& tracks, soa::Filtered<MyPairCandidatesSelected> const& dileptons)
  {
    runDileptonHadron<VarManager::kBtoJpsiEEK, gkEventFillMap, gkTrackFillMap>(event, tracks, dileptons);
  }
  void processSkimmedDimuon(MyEventsSelected::iterator const& event, MyBarrelTracksSelected const& tracks, soa::Filtered<MyPairCandidatesSelected> const& dileptons)
  {
    runDileptonHadron<VarManager::kDecayToMuMu, gkEventFillMap, gkMuonFillMap>(event, tracks, dileptons);
  }
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  // same for Covariance...
  //    void processSkimmed(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelectedWithCov const& tracks, soa::Filtered<MyPairCandidatesSelected> const& dileptons)
  //  {
  //    runDileptonHadron<VarManager::kBtoJpsiEEK, gkEventFillMapWithCov, gkTrackFillMapWithCov>(event, tracks, dileptons);
  //  }
  //  void processSkimmedDimuon(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelectedWithCov  const& tracks, soa::Filtered<MyPairCandidatesSelected> const& dileptons)
  //  {
  //    runDileptonHadron<VarManager::kDecayToMuMu, gkEventFillMapWithCov, gkMuonFillMapWithCov>(event, tracks, dileptons);
  //  }
  //  void processDummy(MyEvents&)
  //  {
  //    // do nothing
  //  }

  PROCESS_SWITCH(DqCumulantFlow, processSkimmed, "Run dilepton-hadron pairing, using skimmed data", false);
  PROCESS_SWITCH(DqCumulantFlow, processSkimmedDimuon, "Run dilepton-hadron pairing, using skimmed data", false);
  PROCESS_SWITCH(DqCumulantFlow, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DqCumulantFlow>(cfgc)};
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
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", histName);
    }

    if (classStr.Contains("Track") && !classStr.Contains("Pairs")) {

      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "barrel");
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
  } // end loop over histogram classes
}

void FillFC(GFW* fGFW, OutputObj<FlowContainer> fFC, const GFW::CorrConfig& corrconf, const double& cent, const double& rndm, bool fillflag)
{
  // Calculate the correlations from the GFW
  double dnx, dny, valx;
  dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
  dny = fGFW->Calculate(corrconf, 0, kTRUE).imag();
  if (dnx == 0) {
    return;
  }

  if (!corrconf.pTDif) {
    valx = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
    if (TMath::Abs(valx) < 1) {
      fFC->FillProfile(corrconf.Head.c_str(), cent, valx, 1, rndm);
      if (dny == 0) {
        return;
      }
    }
    return;
  }
  uint8_t nAxisPtBins = 31;
  for (int i = 1; i <= nAxisPtBins; i++) {
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0) {
      return;
    }
    valx = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
    if (TMath::Abs(valx) < 1) {
      // Fill the charged particle correlation vs pT profiles
      fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, valx, 1., rndm);
    }
    return;
  }
}
