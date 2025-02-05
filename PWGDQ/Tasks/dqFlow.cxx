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
///       Task to compute Q vectors and other quantities related from the generic framework
///       Generic framework O2 version is a port of the AliPhysics version
///       To be used in the DQ analyses aiming for flow measurements
///       Run the standalone task with:
///       o2-analysis-timestamp --aod-file AO2D.root -b | o2-analysis-event-selection -b | o2-analysis-multiplicity-table -b | o2-analysis-centrality-table -b | o2-analysis-fdd-converter -b | o2-analysis-trackselection -b | o2-analysis-trackextension -b | o2-analysis-pid-tpc-full -b | o2-analysis-pid-tof-full -b | o2-analysis-pid-tof-base -b | o2-analysis-pid-tof-beta -b | o2-analysis-dq-flow -b
///       tested (June 2, 2022) on AO2D.root files from train production 242

#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Qvectors.h"

using std::complex;
using std::cout;
using std::endl;
using std::pow;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::analysis;

// Declarations of various short names
using MyBcs = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
using MyEventsWithCentRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
// using MyEventsWithCentQvectRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
// using MyEventsWithCentQvectRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorFV0AVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs, aod::QvectorTPCallVecs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
using MyEventsWithCentQvectRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorTPCposs, aod::QvectorTPCnegs, aod::QvectorTPCalls, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;
using MyMuons = aod::FwdTracks;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCentRun2;
constexpr static uint32_t gkEventFillMapRun3 = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent;
constexpr static uint32_t gkEventFillMapRun3Qvect = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent | VarManager::ObjTypes::CollisionQvectCentr;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct DQEventQvector {

  std::vector<double> ptbinning = {
    0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
    0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
    1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
    2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};

  Produces<ReducedEventsQvector> eventQvector;
  Produces<ReducedEventsQvectorExtra> eventQvectorExtra;
  Produces<ReducedEventsQvectorCentr> eventQvectorCentr;
  Produces<ReducedEventsQvectorCentrExtra> eventQvectorCentrExtra;
  Produces<ReducedEventsRefFlow> eventRefFlow;
  Produces<ReducedEventsQvectorZN> eventQvectorZN;
  Produces<ReducedZdc> eventReducedZdc;
  Produces<ReducedZdcExtra> eventReducedZdcExtra;

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", true, "If true, fill QA histograms"};
  Configurable<bool> fConfigFlow{"cfgFlow", true, "If true, fill Flow histograms"};
  Configurable<bool> fConfigFillWeights{"cfgFillWeights", false, "If true, fill histogram for weights"};
  Configurable<float> fConfigCutPtMin{"cfgCutPtMin", 0.2f, "Minimal pT for tracks"};
  Configurable<float> fConfigCutPtMax{"cfgCutPtMax", 12.0f, "Maximal pT for tracks"};
  Configurable<float> fConfigCutEtaMin{"cfgCutEtaMin", -0.8f, "Eta min range for tracks"};
  Configurable<float> fConfigCutEtaMax{"cfgCutEtaMax", 0.8f, "Eta max range for tracks"};
  Configurable<int> fConfigCutTPCNClMin{"cfgCutTPCNclMin", 0, "Min requirement for number of TPC clusters"};
  Configurable<float> fConfigEtaLimitMin{"cfgEtaLimitMin", -0.4f, "Eta gap min separation, only if using subEvents"};
  Configurable<float> fConfigEtaLimitMax{"cfgEtaLimitMax", 0.4f, "Eta gap max separation, only if using subEvents"};
  // Configurable<uint> fConfigNPow{"cfgNPow", 0, "Power of weights for Q vector"};
  // Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 3.0, {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}}, "Configuration for binning"};

  // Access to the efficiencies and acceptances from CCDB
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> fConfigEfficiency{"ccdb-path-efficiency", "Users/r/rcaron/efficiency", "CCDB path to efficiency object"};
  Configurable<std::string> fConfigAcceptance{"ccdb-path-acceptance", "", "CCDB path to acceptance or GFWWeights object"};
  Configurable<std::string> fConfigURL{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // Configurables for FlowContainer (e.g charged particles pt-differential v2{2}, v2{3}, ...)
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1}, "multiplicity / centrality axis for histograms"};

  // Define the filter for barrel tracks and forward tracks
  Filter trackFilter = (requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true);
  Filter fwdFilter = (aod::fwdtrack::eta < -2.45f) && (aod::fwdtrack::eta > -3.6f);

  // Histograms used for optionnal efficiency and non-uniform acceptance corrections
  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance = nullptr;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  HistogramManager* fHistMan = nullptr;
  AnalysisCompositeCut* fEventCut;
  OutputObj<THashList> fOutputList{"outputQA"};
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};

  // Define global variables for generic framework
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;

  // Initialize CCDB, efficiencies and acceptances from CCDB, histograms, GFW, FlowContainer
  void init(o2::framework::InitContext&)
  {
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut = new AnalysisCompositeCut(true);
    if (!eventCutStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(eventCutStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fEventCut->AddCut(dqcuts::GetAnalysisCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    ccdb->setURL(fConfigURL.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(fConfigNoLaterThan.value);

    VarManager::SetDefaultVarNames();

    if (fConfigQA) {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      DefineHistograms(fHistMan, "Event_BeforeCuts_GFW;Event_AfterCuts_GFW;Event_BeforeCuts_centralFW;Event_AfterCuts_centralFW"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    // Weights for acceptance
    int ptbins = ptbinning.size() - 1;
    fPtAxis = new TAxis(ptbins, &ptbinning[0]);
    if (fConfigFillWeights) {
      // fWeights->SetPtBins(ptbins, &ptbinning[0]); // in the default case, it will accept everything
      fWeights->init(true, false); // true for data, false for MC
    }

    // Reference flow
    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("ChGap22", "ChGap22"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap24", "ChGap24"));   // for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChFull22", "ChFull22")); // no-gap case
    oba->Add(new TNamed("ChFull24", "ChFull24")); // no-gap case
    oba->Add(new TNamed("ChGap32", "ChGap32"));   // gap-case
    fFC->SetName("FlowContainer");
    // fFC->SetXAxis(fPtAxis); // For Robin's first version, there was no configuration of axis for flow container
    fFC->Initialize(oba, axisMultiplicity, 10);
    delete oba;

    // Define regions of positive and negative eta in order to create gaps
    fGFW->AddRegion("refN", fConfigCutEtaMin, fConfigEtaLimitMin, 1, 1);
    fGFW->AddRegion("refP", fConfigEtaLimitMax, fConfigCutEtaMax, 1, 1);
    fGFW->AddRegion("full", fConfigCutEtaMin, fConfigCutEtaMax, 1, 2);
    // Defined the different charged particle correlations
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2} refN {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2 2} refN {-2 -2}", "ChGap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {3 -3}", "ChFull32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {3} refN {-3}", "ChGap32", kFALSE));

    fGFW->CreateRegions();
  }

  template <uint32_t TTrackFillMap, typename TTrack>
  bool isTrackSelected(TTrack const& track)
  {
    constexpr bool isBarrelTrack = ((TTrackFillMap & VarManager::ObjTypes::Track) > 0);
    if constexpr (isBarrelTrack) {
      if (!(track.pt() > fConfigCutPtMin && track.pt() < fConfigCutPtMax)) {
        return false;
      }
      if (!(track.eta() > fConfigCutEtaMin && track.eta() < fConfigCutEtaMax)) {
        return false;
      }
      if (track.tpcNClsFound() <= fConfigCutTPCNClMin) {
        return false;
      }
      if (!track.passedITSNCls()) {
        return false;
      }
      if (!track.passedITSChi2NDF()) {
        return false;
      }
      if (!track.passedITSHits()) {
        return false;
      }
      if (!track.passedTPCCrossedRowsOverNCls()) {
        return false;
      }
      if (!track.passedTPCChi2NDF()) {
        return false;
      }
      if (!track.passedDCAxy()) {
        return false;
      }
      if (!track.passedDCAz()) {
        return false;
      }
    }
    return true;
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (cfg.correctionsLoaded) {
      return;
    }
    if (fConfigAcceptance.value.empty() == false) {
      cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(fConfigAcceptance, timestamp);
      if (cfg.mAcceptance) {
        LOGF(info, "Loaded acceptance weights from %s (%p)", fConfigAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      } else {
        LOGF(warning, "Could not load acceptance weights from %s (%p)", fConfigAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      }
    }
    if (fConfigEfficiency.value.empty() == false) {
      cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(fConfigEfficiency, timestamp);
      if (cfg.mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", fConfigAcceptance.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", fConfigAcceptance.value.c_str(), (void*)cfg.mEfficiency);
    }
    cfg.correctionsLoaded = true;
  }

  // Fill the FlowContainer
  void FillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm, bool /*fillflag*/)
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

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, typename TEvent>
  void runFillQvectorFromCentralFW(TEvent const& collision)
  {
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<TEventFillMap>(collision);

    VarManager::FillQVectorFromCentralFW(collision);

    if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts_centralFW", VarManager::fgValues);
      if (fEventCut->IsSelected(VarManager::fgValues)) {
        fHistMan->FillHistClass("Event_AfterCuts_centralFW", VarManager::fgValues);
      }
    }

    // Fill the tree for the reduced event table with Q vector quantities
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      eventQvectorCentr(collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.qvecFT0CRe(), collision.qvecFT0CIm(), collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.qvecFV0ARe(), collision.qvecFV0AIm(), collision.qvecTPCposRe(), collision.qvecTPCposIm(), collision.qvecTPCnegRe(), collision.qvecTPCnegIm(),
                        collision.sumAmplFT0A(), collision.sumAmplFT0C(), collision.sumAmplFT0M(), collision.sumAmplFV0A(), collision.nTrkTPCpos(), collision.nTrkTPCneg());
      eventQvectorCentrExtra(collision.qvecTPCallRe(), collision.qvecTPCallIm(), collision.nTrkTPCall());
    }
  }

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runFillQvector(TEvent const& collision, MyBcs const&, TTracks const& tracks1, aod::Zdcs const&)
  {
    // Fill the event properties within the VarManager
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<TEventFillMap>(collision);

    // TODO: bc that could be used later to get timestamp for acceptance/GFWWeights
    // auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    // if (fConfigAcceptance.value.empty() == false) { cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(fConfigAcceptance.value, bc.timestamp());}

    fGFW->Clear();
    float centrality;
    // auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    auto bc = collision.template bc_as<MyBcs>();
    // Load weights from CCDB
    loadCorrections(bc.timestamp());
    constexpr bool eventHasCentRun2 = ((TEventFillMap & VarManager::ObjTypes::CollisionCentRun2) > 0);
    constexpr bool eventHasCentRun3 = ((TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0);

    // Get centrality for current collision
    if constexpr (eventHasCentRun2) {
      centrality = VarManager::fgValues[VarManager::kCentVZERO];
    }
    if constexpr (eventHasCentRun3) {
      centrality = VarManager::fgValues[VarManager::kCentFT0C];
    }

    // Acceptance and efficiency weights
    float weff = 1.0, wacc = 1.0;

    // Fill the GFW object in the track loop
    for (auto& track : tracks1) {

      // Selections for barrel tracks
      if (!isTrackSelected<TTrackFillMap>(track)) {
        continue;
      }

      // Fill weights for Q-vector correction: this should be enabled for a first run to get weights
      if (fConfigFillWeights) {
        fWeights->fill(track.phi(), track.eta(), collision.posZ(), track.pt(), centrality, 0);
      }

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
        wacc = cfg.mAcceptance->getNUA(track.phi(), track.eta(), collision.posZ());
      } else {
        wacc = 1.0;
      }
      // Fill the GFW for each track to compute Q vector and correction using weights
      fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, 3); // using default values for ptin = 0 and mask = 3
    }

    float l_Random = fRndm->Rndm(); // used only to compute correlators
    bool fillFlag = kFALSE;         // could be used later
    for (uint64_t l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      if constexpr (eventHasCentRun2) {
        FillFC(corrconfigs.at(l_ind), VarManager::fgValues[VarManager::kCentVZERO], l_Random, fillFlag);
      }
      if constexpr (eventHasCentRun3) {
        FillFC(corrconfigs.at(l_ind), VarManager::fgValues[VarManager::kCentFT0C], l_Random, fillFlag);
      }
    }

    // Define quantities needed for the different eta regions using weighed Q-vectors
    double S10N = 0.;
    double S10P = 0.;
    double S10Full = 0.;
    double S11N = 0.;
    double S11P = 0.;
    double S11Full = 0.;
    double S12Full = 0.;
    double S13Full = 0.;
    double S21Full = 0.;
    double S22Full = 0.;
    double S14Full = 0.;
    double S31Full = 0.;
    double S41Full = 0.;
    complex<double> Q11N(0., 0.);
    complex<double> Q11P(0., 0.);
    complex<double> Q11Full(0., 0.);
    complex<double> Q21N(0., 0.);
    complex<double> Q21P(0., 0.);
    complex<double> Q21Full(0., 0.);
    complex<double> Q23Full(0., 0.);
    complex<double> Q31N(0., 0.);
    complex<double> Q31P(0., 0.);
    complex<double> Q31Full(0., 0.);
    complex<double> Q41N(0., 0.);
    complex<double> Q41P(0., 0.);
    complex<double> Q41Full(0., 0.);
    complex<double> Q42Full(0., 0.);

    if (fGFW && (tracks1.size() > 0)) {
      // Obtain the GFWCumulant where Q is calculated (index=region, with different eta gaps)
      GFWCumulant gfwCumN = fGFW->GetCumulant(0);
      GFWCumulant gfwCumP = fGFW->GetCumulant(1);
      GFWCumulant gfwCumFull = fGFW->GetCumulant(2);

      // S(1,0) for event multiplicity
      S10N = gfwCumN.Vec(0, 0).real();
      S10P = gfwCumP.Vec(0, 0).real();
      S10Full = gfwCumFull.Vec(0, 0).real();

      // S(1,1) for weighted sum
      S11N = gfwCumN.Vec(0, 1).real();
      S11P = gfwCumP.Vec(0, 1).real();
      S11Full = gfwCumFull.Vec(0, 1).real();

      // Extra S(n,p) for cumulants
      S12Full = gfwCumFull.Vec(0, 2).real();
      S13Full = gfwCumFull.Vec(0, 3).real();
      S14Full = gfwCumFull.Vec(0, 4).real();
      S21Full = pow(gfwCumFull.Vec(0, 1).real(), 2);
      S22Full = pow(gfwCumFull.Vec(0, 2).real(), 2);
      S31Full = pow(gfwCumFull.Vec(0, 1).real(), 3);
      S41Full = pow(gfwCumFull.Vec(0, 1).real(), 4);

      // Get the Q vector for selected harmonic, power (for minPt=0)
      Q11N = gfwCumN.Vec(1, 1);
      Q11P = gfwCumP.Vec(1, 1);
      Q11Full = gfwCumFull.Vec(1, 1);

      Q21N = gfwCumN.Vec(2, 1);
      Q21P = gfwCumP.Vec(2, 1);
      Q21Full = gfwCumFull.Vec(2, 1);

      Q31N = gfwCumN.Vec(3, 1);
      Q31P = gfwCumP.Vec(3, 1);
      Q31Full = gfwCumFull.Vec(3, 1);

      Q41N = gfwCumN.Vec(4, 1);
      Q41P = gfwCumP.Vec(4, 1);
      Q41Full = gfwCumFull.Vec(4, 1);

      // Extra Q-vectors for cumulants
      Q23Full = gfwCumFull.Vec(2, 3);
      Q42Full = gfwCumFull.Vec(4, 2);
    }

    // Fill the VarManager::fgValues with the Q vector quantities
    VarManager::FillQVectorFromGFW(collision, Q11Full, Q11N, Q11P, Q21Full, Q21N, Q21P, Q31Full, Q31N, Q31P, Q41Full, Q41N, Q41P, Q23Full, Q42Full, S10Full, S10N, S10P, S11Full, S11N, S11P, S12Full, S13Full, S14Full, S21Full, S22Full, S31Full, S41Full);

    if (fConfigQA) {
      if ((tracks1.size() > 0) && (S10N * S10P * S10Full != 0.0)) {
        fHistMan->FillHistClass("Event_BeforeCuts_GFW", VarManager::fgValues);
        if (fEventCut->IsSelected(VarManager::fgValues)) {
          fHistMan->FillHistClass("Event_AfterCuts_GFW", VarManager::fgValues);
        }
      }
    }

    // Fill the tree for the reduced event table with Q vector quantities
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      eventQvector(VarManager::fgValues[VarManager::kQ1X0A], VarManager::fgValues[VarManager::kQ1Y0A], VarManager::fgValues[VarManager::kQ1X0B], VarManager::fgValues[VarManager::kQ1Y0B], VarManager::fgValues[VarManager::kQ1X0C], VarManager::fgValues[VarManager::kQ1Y0C], VarManager::fgValues[VarManager::kQ2X0A], VarManager::fgValues[VarManager::kQ2Y0A], VarManager::fgValues[VarManager::kQ2X0B], VarManager::fgValues[VarManager::kQ2Y0B], VarManager::fgValues[VarManager::kQ2X0C], VarManager::fgValues[VarManager::kQ2Y0C], VarManager::fgValues[VarManager::kMultA], VarManager::fgValues[VarManager::kMultB], VarManager::fgValues[VarManager::kMultC], VarManager::fgValues[VarManager::kQ3X0A], VarManager::fgValues[VarManager::kQ3Y0A], VarManager::fgValues[VarManager::kQ3X0B], VarManager::fgValues[VarManager::kQ3Y0B], VarManager::fgValues[VarManager::kQ3X0C], VarManager::fgValues[VarManager::kQ3Y0C], VarManager::fgValues[VarManager::kQ4X0A], VarManager::fgValues[VarManager::kQ4Y0A], VarManager::fgValues[VarManager::kQ4X0B], VarManager::fgValues[VarManager::kQ4Y0B], VarManager::fgValues[VarManager::kQ4X0C], VarManager::fgValues[VarManager::kQ4Y0C]);
      eventQvectorExtra(VarManager::fgValues[VarManager::kQ42XA], VarManager::fgValues[VarManager::kQ42YA], VarManager::fgValues[VarManager::kQ23XA], VarManager::fgValues[VarManager::kQ23YA], VarManager::fgValues[VarManager::kS11A], VarManager::fgValues[VarManager::kS12A], VarManager::fgValues[VarManager::kS13A], VarManager::fgValues[VarManager::kS31A]);
      eventRefFlow(VarManager::fgValues[VarManager::kM11REF], VarManager::fgValues[VarManager::kM11REFetagap], VarManager::fgValues[VarManager::kM1111REF], VarManager::fgValues[VarManager::kCORR2REF], VarManager::fgValues[VarManager::kCORR2REFetagap], VarManager::fgValues[VarManager::kCORR4REF], centrality);
    }

    if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionQvectCentr) > 0) {
      VarManager::FillQVectorFromCentralFW(collision);

      if (bc.has_zdc()) {
        auto zdc = bc.zdc();
        VarManager::FillSpectatorPlane(zdc);
      }

      if ((tracks1.size() > 0) && (VarManager::fgValues[VarManager::kMultA] * VarManager::fgValues[VarManager::kMultB] * VarManager::fgValues[VarManager::kMultC] != 0.0)) {
        if (fConfigQA) {
          fHistMan->FillHistClass("Event_BeforeCuts_centralFW", VarManager::fgValues);
          if (fEventCut->IsSelected(VarManager::fgValues)) {
            fHistMan->FillHistClass("Event_AfterCuts_centralFW", VarManager::fgValues);
          }
        }
      }
      if (fEventCut->IsSelected(VarManager::fgValues)) {
        eventQvectorCentr(collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.qvecFT0CRe(), collision.qvecFT0CIm(), collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.qvecFV0ARe(), collision.qvecFV0AIm(), collision.qvecTPCposRe(), collision.qvecTPCposIm(), collision.qvecTPCnegRe(), collision.qvecTPCnegIm(),
                          collision.sumAmplFT0A(), collision.sumAmplFT0C(), collision.sumAmplFT0M(), collision.sumAmplFV0A(), collision.nTrkTPCpos(), collision.nTrkTPCneg());
        eventQvectorCentrExtra(collision.qvecTPCallRe(), collision.qvecTPCallIm(), collision.nTrkTPCall());
        if (bc.has_zdc()) {
          eventQvectorZN(VarManager::fgValues[VarManager::kQ1ZNAX], VarManager::fgValues[VarManager::kQ1ZNAY], VarManager::fgValues[VarManager::kQ1ZNCX], VarManager::fgValues[VarManager::kQ1ZNCY]);
          eventReducedZdc(VarManager::fgValues[VarManager::kEnergyCommonZNA], VarManager::fgValues[VarManager::kEnergyCommonZNC], VarManager::fgValues[VarManager::kEnergyCommonZPA], VarManager::fgValues[VarManager::kEnergyCommonZPC],
                          VarManager::fgValues[VarManager::kTimeZNA], VarManager::fgValues[VarManager::kTimeZNC], VarManager::fgValues[VarManager::kTimeZPA], VarManager::fgValues[VarManager::kTimeZPC]);
          eventReducedZdcExtra(VarManager::fgValues[VarManager::kEnergyZNA1], VarManager::fgValues[VarManager::kEnergyZNA2], VarManager::fgValues[VarManager::kEnergyZNA3], VarManager::fgValues[VarManager::kEnergyZNA4],
                               VarManager::fgValues[VarManager::kEnergyZNC1], VarManager::fgValues[VarManager::kEnergyZNC2], VarManager::fgValues[VarManager::kEnergyZNC3], VarManager::fgValues[VarManager::kEnergyZNC4]);
        } else {
          eventQvectorZN(-999, -999, -999, -999);
          eventReducedZdc(-999, -999, -999, -999, -999, -999, -999, -999);
          eventReducedZdcExtra(-999, -999, -999, -999, -999, -999, -999, -999);
        }
      }
    }
  }

  // Process to fill Q vector using barrel tracks in a reduced event table for barrel/muon tracks flow related analyses Run 2
  void processBarrelQvectorRun2(MyEventsWithCent::iterator const& collisions, MyBcs const& bcs, soa::Filtered<MyBarrelTracks> const& tracks, aod::Zdcs const& zdcs)
  {
    runFillQvector<gkEventFillMap, gkTrackFillMap>(collisions, bcs, tracks, zdcs);
  }

  // Process to fill Q vector using barrel tracks in a reduced event table for barrel/muon tracks flow related analyses Run 3
  void processBarrelQvector(MyEventsWithCentRun3::iterator const& collisions, MyBcs const& bcs, soa::Filtered<MyBarrelTracks> const& tracks, aod::Zdcs const& zdcs)
  {
    runFillQvector<gkEventFillMapRun3, gkTrackFillMap>(collisions, bcs, tracks, zdcs);
  }

  // Process to fill Q vector using barrel tracks in a reduced event table for barrel/muon tracks flow related analyses Run 3
  void processAllQvector(MyEventsWithCentQvectRun3::iterator const& collisions, MyBcs const& bcs, soa::Filtered<MyBarrelTracks> const& tracks, aod::Zdcs const& zdcs)
  {
    runFillQvector<gkEventFillMapRun3Qvect, gkTrackFillMap>(collisions, bcs, tracks, zdcs);
  }

  // Process to fill Q vector using barrel tracks in a reduced event table for barrel/muon tracks flow related analyses Run 3
  void processCentralQvector(MyEventsWithCentQvectRun3::iterator const& collisions)
  {
    runFillQvectorFromCentralFW<gkEventFillMapRun3>(collisions);
  }

  // Process to fill Q vector using forward tracks in a reduced event table for barrel/muon tracks flow related analyses Run 3
  void processForwardQvector(MyEventsWithCentRun3::iterator const& collisions, MyBcs const& bcs, soa::Filtered<aod::MFTTracks> const& tracks, aod::Zdcs const& zdcs)
  {
    runFillQvector<gkEventFillMapRun3, 0u>(collisions, bcs, tracks, zdcs);
  }

  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQEventQvector, processBarrelQvectorRun2, "Run q-vector task on barrel tracks for Run2", false);
  PROCESS_SWITCH(DQEventQvector, processBarrelQvector, "Run q-vector task on barrel tracks for Run3", false);
  PROCESS_SWITCH(DQEventQvector, processAllQvector, "Run q-vector task on barrel tracks for Run3 and using central q-vector", false);
  PROCESS_SWITCH(DQEventQvector, processCentralQvector, "Run q-vector task using central q-vector", false);
  PROCESS_SWITCH(DQEventQvector, processForwardQvector, "Run q-vector task on forward tracks for Run3", false);
  PROCESS_SWITCH(DQEventQvector, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQEventQvector>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "qvector,cross,trigger,cent,res");
    }

    if (classStr.Contains("Track")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "cent,kine,time,map,its,clssize,tpc");
    }
  }
}
