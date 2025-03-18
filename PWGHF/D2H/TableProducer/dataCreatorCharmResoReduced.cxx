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

/// \file dataCreatorCharmResoReduced.cxx
/// \brief Creation of D-V0 pairs
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, UniTO Turin
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "DetectorsBase/Propagator.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

// event types
enum Event : uint8_t {
  Processed = 0,
  NoDV0Selected,
  DV0Selected,
  kNEvent
};

enum DecayChannel : uint8_t {
  DstarV0 = 0,
  DplusV0,
  DstarTrack
};

enum BachelorType : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda,
  Track
};

enum DType : uint8_t {
  Dplus = 1,
  Dstar
};

enum DecayTypeMc : uint8_t {
  Ds1ToDStarK0ToD0PiK0s = 1,
  Ds2StarToDplusK0sToPiKaPiPiPi,
  Ds1ToDStarK0ToDPlusPi0K0s,
  Ds1ToDStarK0ToD0PiK0sPart,
  Ds1ToDStarK0ToD0NoPiK0sPart,
  Ds1ToDStarK0ToD0PiK0sOneMu,
  Ds2StarToDplusK0sOneMu
};

enum PartialMatchMc : uint8_t {
  K0Matched = 0,
  D0Matched,
  DStarMatched,
  DPlusMatched,
  K0MuMatched,
  DStarMuMatched
};

/// Creation of D-V0 pairs
struct HfDataCreatorCharmResoReduced {

  // Produces AOD tables to store collision information
  Produces<aod::HfRedCollisions> hfReducedCollision; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfOrigColCounts> hfCollisionCounter; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  // tracks, V0 and D candidates reduced tables
  Produces<aod::HfRedVzeros> hfCandV0;            // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfRedTrkNoParams> hfTrackNoParam; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfRed3PrNoTrks> hfCandD;          // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  // ML optional Tables
  Produces<aod::HfRed3ProngsMl> hfCandDMl; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  // MC Tables
  Produces<aod::HfMcRecRedDV0s> rowHfDV0McRecReduced;
  Produces<aod::HfMcGenRedResos> rowHfResoMcGenReduced;

  // CCDB configuration
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Configurable<bool> propagateV0toPV{"propagateV0toPV", false, "Enable or disable V0 propagation to V0"};
  Configurable<bool> doMcRecQa{"doMcRecQa", true, "Fill QA histograms for Mc matching"};

  int runNumber{0}; // needed to detect if the run changed and trigger update of calibrations etc.
  // selection D
  struct : ConfigurableGroup {
    std::string prefix = "dmesons";
    Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for D"};
    Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 & Pi"};
  } cfgDmesCuts;

  // selection V0
  struct : ConfigurableGroup {
    std::string prefix = "v0s";
    Configurable<float> deltaMassK0s{"deltaMassK0s", 0.02, "delta mass cut for K0S"};
    Configurable<float> deltaMassLambda{"deltaMassLambda", 0.01, "delta mass cut for Lambda"};
    Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};
    Configurable<float> etaMaxDau{"etaMaxDau", 5.f, "maximum eta V0 daughters"};
    Configurable<float> trackNclusItsCut{"trackNclusItsCut", 0, "Minimum number of ITS clusters for V0 daughter"};
    Configurable<int> trackNCrossedRowsTpc{"trackNCrossedRowsTpc", 50, "Minimum TPC crossed rows"};
    Configurable<float> trackNsharedClusTpc{"trackNsharedClusTpc", 1000, "Maximum number of shared TPC clusters for V0 daughter"};
    Configurable<float> dcaDau{"dcaDau", 1.f, "DCA V0 daughters"};
    Configurable<float> dcaMaxDauToPv{"dcaMaxDauToPv", 0.1f, "Maximum daughter's DCA to PV"};
    Configurable<float> dcaPv{"dcaPv", 1.f, "DCA V0 to PV"};
    Configurable<double> cosPa{"cosPa", 0.99f, "V0 CosPA"};
    Configurable<float> radiusMin{"radiusMin", 0.9f, "Minimum v0 radius accepted"};
    Configurable<float> nSigmaTpc{"nSigmaTpc", 4.f, "Nsigmatpc"};
    Configurable<float> nSigmaTofPr{"nSigmaTofPr", 4.f, "N sigma TOF for protons only"};
  } cfgV0Cuts;

  // selection single tracks
  struct : ConfigurableGroup {
    std::string prefix = "single_tracks";
    Configurable<int> setTrackSelections{"setTrackSelections", 2, "flag to apply track selections: 0=none; 1=global track w/o DCA selection; 2=global track; 3=only ITS quality"};
    Configurable<float> maxEta{"maxEta", 0.8, "maximum pseudorapidity for single tracks to be paired with D mesons"};
    Configurable<float> minPt{"minPt", 0.1, "minimum pT for single tracks to be paired with D mesons"};
    Configurable<float> maxNsigmaTpcPi{"maxNsigmaTpcPi", -1., "maximum pion NSigma in TPC for single tracks to be paired with D mesons; set negative to reject"};
    Configurable<float> maxNsigmaTpcKa{"maxNsigmaTpcKa", -1., "maximum kaon NSigma in TPC for single tracks to be paired with D mesons; set negative to reject"};
    Configurable<float> maxNsigmaTpcPr{"maxNsigmaTpcPr", 3., "maximum proton NSigma in TPC for single tracks to be paired with D mesons; set negative to reject"};
  } cfgSingleTrackCuts;

  // other configurables
  Configurable<bool> rejectPairsWithCommonDaughter{"rejectPairsWithCommonDaughter", true, "flag to reject already at this stage the pairs that share a daughter track"};

  // material correction for track propagation
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  HfHelper hfHelper;
  o2::hf_evsel::HfEventSelection hfEvSel;
  o2::vertexing::DCAFitterN<2> fitter;
  double bz{0.};
  Service<o2::framework::O2DatabasePDG> pdg;

  // bool isHfCandResoConfigFilled = false;

  // Helper struct to pass V0 informations
  struct {
    std::array<float, 3> pos;
    std::array<float, 3> mom;
    std::array<float, 3> momPos;
    std::array<float, 3> momNeg;
    float pT;
    float cosPA;
    float dcaV0ToPv;
    float dcaDau;
    float alpha;
    float eta;
    float radius;
    float mK0Short;
    float mLambda;
    uint8_t v0Type;
  } candidateV0;

  struct {
    float invMassD;
    float ptD;
    float invMassDdau;
    float invMassKPiPiV0;
    float ptReso;
  } varUtils;
  using CandsDplusFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandsDplusFilteredWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandDstarFiltered = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi>>;
  using CandDstarFilteredWithMl = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>>;
  using TracksWithPID = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using TracksIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTOFFullPi, aod::pidTPCPr, aod::pidTOFFullPr>;
  using TracksIUWithPIDAndMC = soa::Join<TracksIUWithPID, aod::McTrackLabels>;

  Filter filterSelectDplus = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= cfgDmesCuts.selectionFlagDplus);
  Filter filterSelectedCandDstar = (aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == cfgDmesCuts.selectionFlagDstarToD0Pi);

  Preslice<CandsDplusFiltered> candsDplusPerCollision = aod::hf_cand::collisionId;
  Preslice<CandsDplusFilteredWithMl> candsDplusPerCollisionWithMl = aod::hf_cand::collisionId;
  Preslice<CandDstarFiltered> candsDstarPerCollision = aod::hf_cand::collisionId;
  Preslice<CandDstarFilteredWithMl> candsDstarPerCollisionWithMl = aod::hf_cand::collisionId;
  Preslice<aod::V0s> candsV0PerCollision = aod::v0::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // histograms
    constexpr int kNBinsEvents = kNEvent;
    std::string labels[kNBinsEvents];
    labels[Event::Processed] = "processed";
    labels[Event::NoDV0Selected] = "without DV0 pairs";
    labels[Event::DV0Selected] = "with DV0 pairs";
    static const AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    const AxisSpec axisPt{50, 0.f, 50.f, ""};
    const AxisSpec axisP{100, 0.f, 10.f, ""};
    const AxisSpec axisDeDx{500, 0.f, 1000.f, ""};
    const AxisSpec axisMassDplus{200, 1.7f, 2.1f, ""};
    const AxisSpec axisMassDstar{200, 0.139f, 0.179f, ""};
    const AxisSpec axisMassLambda{100, 1.05f, 1.35f, ""};
    const AxisSpec axisMassKzero{100, 0.35f, 0.65f, ""};
    const AxisSpec axisMassDsj{400, 0.49f, 0.89f, ""};

    registry.add("hMassVsPtDplusAll", "Dplus candidates (all, regardless the pairing with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDplus}});
    registry.add("hMassVsPtDstarAll", "Dstar candidates (all, regardless the pairing with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDstar}});
    registry.add("hMassVsPtDplusPaired", "Dplus candidates (paired with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDplus}});
    registry.add("hMassVsPtDstarPaired", "Dstar candidates (paired with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDstar}});

    registry.add("hMassVsPtK0s", "K0^{s} candidates;#it{p}_{T} (GeV/#it{c});inv. mass (#pi^{#plus}#pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassKzero}});
    registry.add("hMassVsPtLambda", "Lambda candidates;#it{p}_{T} (GeV/#it{c});inv. mass (p #pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassLambda}});
    registry.add("hdEdxVsP", "Tracks;#it{p} (GeV/#it{c});d#it{E}/d#it{x};entries", {HistType::kTH2F, {axisP, axisDeDx}});

    registry.add("hMassDs1", "Ds1 candidates;m_{Ds1} - m_{D^{*}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassDsj}});
    registry.add("hMassDsStar2", "Ds^{*}2 candidates; Ds^{*}2 - m_{D^{#plus}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassDsj}});
    registry.add("hMassXcRes", "XcRes candidates; XcRes - m_{D^{#plus}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{300, 1.1, 1.4}}});
    registry.add("hMassDstarProton", "D^{*}-proton candidates;m_{D^{*}p} - m_{D^{*}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0.9, 1.4}}});
    registry.add("hDType", "D selection flag", {HistType::kTH1F, {{5, -2.5, 2.5}}});

    registry.add("hMCRecCounter", "Number of Reconstructed MC Matched candidates per channel", {HistType::kTH1F, {{17, -8.5, 8.5}}});
    registry.add("hMCRecDebug", "Debug of MC Reco", {HistType::kTH1F, {{16, -0.5, 15.5}}});
    registry.add("hMCRecOrigin", "Origin of Matched particles", {HistType::kTH1F, {{3, -0.5, 2.5}}});

    registry.add("hMCGenCounter", "Number of Generated particles; Decay Channel Flag; pT [GeV/c]", {HistType::kTH2F, {{17, -8.5, 8.5}, {100, 0, 50}}});
    registry.add("hMCSignCounter", "Sign of Generated particles", {HistType::kTH1F, {{3, -1.5, 1.5}}});
    registry.add("hMCGenOrigin", "Origin of Generated particles", {HistType::kTH1F, {{3, -0.5, 2.5}}});
    registry.add("hMCOriginCounterWrongDecay", "Origin of Generated particles in Wrong decay", {HistType::kTH1F, {{3, -0.5, 2.5}}});

    if (doMcRecQa) {
      registry.add("hMassVsPtK0Matched", "K0s candidates Matched ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassKzero}});
      registry.add("hMassVsPtD0Matched", "D0 candidates Matched ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDplus}});
      registry.add("hMassVsPtDstarMatched", "Dstar candidates Matched ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDstar}});
      registry.add("hMassVsPtDplusMatched", "Dplus candidates Matched ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDplus}});
      registry.add("hMassVsPtDs1Matched", "Ds1 candidates Matched ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDsj}});
      registry.add("hMassVsPtDs2StarMatched", "Ds2Star candidates Matched ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDsj}});
      registry.add("hMassVsPtK0MatchedPiToMu", "K0s candidates Matched with PiToMu decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassKzero}});
      registry.add("hMassVsPtD0MatchedPiToMu", "D0 candidates Matched with PiToMu decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDplus}});
      registry.add("hMassVsPtDstarMatchedPiToMu", "Dstar candidates Matched with PiToMu decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDstar}});
      registry.add("hMassVsPtDplusMatchedPiToMu", "Dplus candidates Matched with PiToMu decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDplus}});
      registry.add("hMassVsPtDs1MatchedPiToMu", "Ds1 candidates Matched with PiToMu decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDsj}});
      registry.add("hMassVsPtDs2StarMatchedPiToMu", "Ds2Star candidates Matched with PiToMu decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDsj}});
      registry.add("hMassVsPtD0MatchedKaToPi", "D0 candidates Matched with KaToPi decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDplus}});
      registry.add("hMassVsPtDstarMatchedKaToPi", "Dstar candidates Matched with KaToPi decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDstar}});
      registry.add("hMassVsPtDplusMatchedKaToPi", "Dplus candidates Matched with KaToPi decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDplus}});
      registry.add("hMassVsPtDs1MatchedKaToPi", "Ds1 candidates Matched with KaToPi decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDsj}});
      registry.add("hMassVsPtDs2StarMatchedKaToPi", "Ds2Star candidates Matched with KaToPi decay ;#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPt, axisMassDsj}});
    }

    // Configure CCDB access
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(url);
    runNumber = 0;
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    // Configure DCA fitter
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxDXYIni(4);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);
  }

  /// Basic track quality selections for V0 daughters
  /// \param Tr is a track
  /// \param dDaughtersIds are the IDs of the D meson daughter tracks
  template <typename Tr>
  bool selectV0Daughter(Tr const& track, const std::array<int, 3>& dDaughtersIds)
  {
    // acceptance selection
    if (std::abs(track.eta()) > cfgV0Cuts.etaMaxDau) {
      return false;
    }
    // Tpc Refit
    if (!(track.hasTPC())) {
      return false;
    }
    // track quality selection
    if (track.itsNCls() < cfgV0Cuts.trackNclusItsCut ||
        track.tpcNClsFound() < cfgV0Cuts.trackNCrossedRowsTpc ||
        track.tpcNClsCrossedRows() < cfgV0Cuts.trackNCrossedRowsTpc ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcNClsShared() > cfgV0Cuts.trackNsharedClusTpc) {
      return false;
    }
    // rejection of tracks that share a daughter with the D meson
    if (rejectPairsWithCommonDaughter && std::find(dDaughtersIds.begin(), dDaughtersIds.end(), track.globalIndex()) != dDaughtersIds.end()) {
      return false;
    }
    return true;
  }

  // Utility to find which v0 daughter carries the largest fraction of the mother longitudinal momentum
  float alphaAP(std::array<float, 3> const& momA, std::array<float, 3> const& momB, std::array<float, 3> const& momC)
  {
    float momTot = std::sqrt(std::pow(momA[0], 2.) + std::pow(momA[1], 2.) + std::pow(momA[2], 2.));
    float lQlPos = (momB[0] * momA[0] + momB[1] * momA[1] + momB[2] * momA[2]) / momTot;
    float lQlNeg = (momC[0] * momA[0] + momC[1] * momA[1] + momC[2] * momA[2]) / momTot;
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
  }
  // Utility to find DCA of V0 to Primary vertex
  float calculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
  {
    return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
  }
  /// Basic selection of V0 candidates
  /// \param collision is the current collision
  /// \param dauTracks are the v0 daughter tracks
  /// \param dDaughtersIds are the IDs of the D meson daughter tracks
  /// \return a bitmap with mass hypotesis if passes all cuts
  template <typename Coll, typename Tr>
  bool buildAndSelectV0(const Coll& collision, const std::array<int, 3>& dDaughtersIds, const std::array<Tr, 2>& dauTracks)
  {
    auto trackPos = dauTracks[0];
    auto trackNeg = dauTracks[1];

    // single-tracks selection
    if (!selectV0Daughter(trackPos, dDaughtersIds) || !selectV0Daughter(trackNeg, dDaughtersIds))
      return false;
    // daughters DCA to V0's collision primary vertex
    gpu::gpustd::array<float, 2> dcaInfo;
    auto trackPosPar = getTrackPar(trackPos);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPosPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
    auto trackPosDcaXY = dcaInfo[0];
    auto trackNegPar = getTrackPar(trackNeg);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackNegPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
    auto trackNegDcaXY = dcaInfo[0];
    if (std::fabs(trackPosDcaXY) < cfgV0Cuts.dcaMaxDauToPv || std::fabs(trackNegDcaXY) < cfgV0Cuts.dcaMaxDauToPv) {
      return false;
    }
    // vertex reconstruction
    auto trackPosCov = getTrackParCov(trackPos);
    auto trackNegCov = getTrackParCov(trackNeg);
    int nCand = 0;
    try {
      nCand = fitter.process(trackPosCov, trackNegCov);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call!";
      return false;
    }
    if (nCand == 0) {
      return false;
    }
    // compute candidate momentum from tracks propagated to decay vertex
    auto& trackPosProp = fitter.getTrack(0);
    auto& trackNegProp = fitter.getTrack(1);
    trackPosProp.getPxPyPzGlo(candidateV0.momPos);
    trackNegProp.getPxPyPzGlo(candidateV0.momNeg);
    for (int i = 0; i < 3; ++i) {
      candidateV0.mom[i] = candidateV0.momPos[i] + candidateV0.momNeg[i];
    }
    candidateV0.pT = std::hypot(candidateV0.mom[0], candidateV0.mom[1]);
    // topological selections:
    // v0 eta
    candidateV0.eta = RecoDecay::eta(candidateV0.mom);
    if (std::abs(candidateV0.eta) > cfgV0Cuts.etaMax) {
      return false;
    }
    // daughters DCA
    candidateV0.dcaDau = std::sqrt(fitter.getChi2AtPCACandidate());
    if (candidateV0.dcaDau > cfgV0Cuts.dcaDau) {
      return false;
    }
    // v0 radius
    const auto& vtx = fitter.getPCACandidate();
    candidateV0.radius = std::hypot(vtx[0], vtx[1]);
    if (candidateV0.radius < cfgV0Cuts.radiusMin) {
      return false;
    }
    for (int i = 0; i < 3; i++) {
      candidateV0.pos[i] = vtx[i];
    }
    // v0 DCA to primary vertex
    candidateV0.dcaV0ToPv = calculateDCAStraightToPV(
      vtx[0], vtx[1], vtx[2],
      candidateV0.momPos[0] + candidateV0.momNeg[0],
      candidateV0.momPos[1] + candidateV0.momNeg[1],
      candidateV0.momPos[2] + candidateV0.momNeg[2],
      collision.posX(), collision.posY(), collision.posZ());
    if (std::abs(candidateV0.dcaV0ToPv) > cfgV0Cuts.dcaPv) {
      return false;
    }
    // v0 cosine of pointing angle
    std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
    candidateV0.cosPA = RecoDecay::cpa(primVtx, vtx, candidateV0.mom);
    if (candidateV0.cosPA < cfgV0Cuts.cosPa) {
      return false;
    }

    // distinguish between K0s, and Lambda hypotesys
    candidateV0.v0Type = {BIT(K0s) | BIT(Lambda) | BIT(AntiLambda)};
    // for lambda hypotesys define if its lambda or anti-lambda
    candidateV0.alpha = alphaAP(candidateV0.mom, candidateV0.momPos, candidateV0.momNeg);
    bool matter = candidateV0.alpha > 0;
    CLRBIT(candidateV0.v0Type, matter ? AntiLambda : Lambda);
    auto massPos = matter ? o2::constants::physics::MassProton : o2::constants::physics::MassPionCharged;
    auto massNeg = matter ? o2::constants::physics::MassPionCharged : o2::constants::physics::MassProton;
    // mass hypotesis
    candidateV0.mLambda = RecoDecay::m(std::array{candidateV0.momPos, candidateV0.momNeg}, std::array{massPos, massNeg});
    candidateV0.mK0Short = RecoDecay::m(std::array{candidateV0.momPos, candidateV0.momNeg}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
    if (std::fabs(candidateV0.mK0Short - MassK0) > cfgV0Cuts.deltaMassK0s) {
      CLRBIT(candidateV0.v0Type, K0s);
    }
    if (std::fabs(candidateV0.mLambda - MassLambda0) > cfgV0Cuts.deltaMassLambda) {
      CLRBIT(candidateV0.v0Type, Lambda);
      CLRBIT(candidateV0.v0Type, AntiLambda);
    }
    // PID
    if (TESTBIT(candidateV0.v0Type, K0s)) {
      if ((trackPos.hasTPC() && std::fabs(trackPos.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc) ||
          (trackNeg.hasTPC() && std::fabs(trackNeg.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc))
        CLRBIT(candidateV0.v0Type, K0s);
    }
    if (TESTBIT(candidateV0.v0Type, Lambda)) {
      if ((trackPos.hasTPC() && std::fabs(trackPos.tpcNSigmaPr()) > cfgV0Cuts.nSigmaTpc) ||
          (trackPos.hasTOF() && std::fabs(trackPos.tofNSigmaPr()) > cfgV0Cuts.nSigmaTofPr) ||
          (trackNeg.hasTPC() && std::fabs(trackNeg.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc))
        CLRBIT(candidateV0.v0Type, Lambda);
    }
    if (TESTBIT(candidateV0.v0Type, AntiLambda)) {
      if ((trackPos.hasTPC() && std::fabs(trackPos.tpcNSigmaPi()) > cfgV0Cuts.nSigmaTpc) ||
          (trackNeg.hasTPC() && std::fabs(trackNeg.tpcNSigmaPr()) > cfgV0Cuts.nSigmaTpc) ||
          (trackNeg.hasTOF() && std::fabs(trackNeg.tofNSigmaPr()) > cfgV0Cuts.nSigmaTofPr))
        CLRBIT(candidateV0.v0Type, AntiLambda);
    }
    if (candidateV0.v0Type == 0) {
      return false;
    }
    return true;
  }

  /// Basic selection of tracks
  /// \param track is the track
  /// \param dDaughtersIds are the IDs of the D meson daughter tracks
  /// \return true if passes all cuts
  template <typename Tr>
  bool isTrackSelected(const Tr& track, const std::array<int, 3>& dDaughtersIds)
  {

    if (rejectPairsWithCommonDaughter && std::find(dDaughtersIds.begin(), dDaughtersIds.end(), track.globalIndex()) != dDaughtersIds.end()) {
      return false;
    }

    switch (cfgSingleTrackCuts.setTrackSelections) {
      case 1:
        if (!track.isGlobalTrackWoDCA()) {
          return false;
        }
        break;
      case 2:
        if (!track.isGlobalTrack()) {
          return false;
        }
        break;
      case 3:
        if (!track.isQualityTrackITS()) {
          return false;
        }
        break;
    }

    if (track.pt() < cfgSingleTrackCuts.minPt) {
      return false;
    }

    if (std::abs(track.eta()) > cfgSingleTrackCuts.maxEta) {
      return false;
    }

    if (!track.hasTPC()) {
      return false;
    }

    bool isPion = std::abs(track.tpcNSigmaPi()) < cfgSingleTrackCuts.maxNsigmaTpcPi;
    bool isKaon = std::abs(track.tpcNSigmaKa()) < cfgSingleTrackCuts.maxNsigmaTpcKa;
    bool isProton = std::abs(track.tpcNSigmaPr()) < cfgSingleTrackCuts.maxNsigmaTpcPr;

    if (!isPion && !isKaon && !isProton) { // we keep the track if is it compatible with at least one of the PID hypotheses selected
      return false;
    }

    return true;
  }

  /// Function for filling MC reco information in the tables
  /// \param particlesMc is the table with MC particles
  /// \param vecDaughtersReso is the vector with all daughter tracks (bachelor pion in last position)
  /// \param indexHfCandCharm is the index of the charm-hadron bachelor in the reduced table
  /// \param indexCandV0 is the index of the v0 bachelor in the reduced table
  template <uint8_t decChannel, typename PParticles, typename TTrack>
  void fillMcRecoInfo(const PParticles& particlesMc,
                      const std::vector<TTrack>& vecDaughtersReso,
                      int& indexHfCandCharm,
                      int& indexCandV0)
  {

    // we check the MC matching to be stored
    int8_t sign{0};
    int8_t signDStar{0};
    int8_t signDPlus{0};
    int8_t signD0{0};
    int8_t signV0{0};
    int8_t flag{0};
    int8_t debug{0};
    int8_t origin{0};
    int8_t nPiToMuReso{0}, nPiToMuV0, nPiToMuD0{0}, nPiToMuDstar{0}, nPiToMuDplus{0};
    int8_t nKaToPiReso{0}, nKaToPiV0, nKaToPiD0{0}, nKaToPiDstar{0}, nKaToPiDplus{0};
    std::vector<int> idxBhadMothers{};
    float motherPt{-1.f};
    int indexRecReso{-1}, indexRecDstar{-1}, indexRecDplus{-1}, indexRecD0{-1}, indexRecK0{-1}, indexRecResoPartReco{-1};

    if constexpr (decChannel == DecayChannel::DstarV0) {
      // Ds1 → D* K0 → (D0 π+) K0s → ((K-π+) π+)(π+π-)
      indexRecD0 = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1]}, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &signD0, 1, &nPiToMuD0, &nKaToPiD0);
      indexRecK0 = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, std::array{vecDaughtersReso[3], vecDaughtersReso[4]}, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signV0, 2, &nPiToMuV0, &nKaToPiV0);
      if (indexRecD0 > -1) {
        indexRecDstar = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2]}, Pdg::kDStar, std::array{-kKPlus, +kPiPlus, +kPiPlus}, true, &signDStar, 2, &nPiToMuDstar, &nKaToPiDstar);
      }
      if (indexRecD0 > -1 && indexRecDstar > -1 && indexRecK0 > -1) {
        indexRecReso = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]}, Pdg::kDS1, std::array{+kPiPlus, -kKPlus, +kPiPlus, +kPiPlus, -kPiPlus}, true, &sign, 3, &nPiToMuReso, &nKaToPiReso);
        if (indexRecReso > -1 && nPiToMuReso == 0 && nKaToPiReso == 0) {
          flag = sign * DecayTypeMc::Ds1ToDStarK0ToD0PiK0s;
        } else if (indexRecReso > -1 && nPiToMuReso >= 1 && nKaToPiReso == 0) {
          flag = sign * DecayTypeMc::Ds1ToDStarK0ToD0PiK0sOneMu;
        }
      }

      // Ds1+ not matched: we check if it is partially reco
      if (indexRecReso < 0) {
        indexRecResoPartReco = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]}, Pdg::kDS1, std::array{+kPiPlus, -kKPlus, +kPiPlus, +kPiPlus, -kPiPlus}, true, &sign, 3);
        indexRecDplus = RecoDecay::getMatchedMCRec(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2]}, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &signDPlus, 2);
        if (indexRecResoPartReco > -1) { // we look for decays of D* or D0 with more daughters
          if (indexRecDstar < 0 && indexRecK0 > -1) {
            auto indexRecDstarPartReco = RecoDecay::getMatchedMCRec<false, false, true, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2]}, Pdg::kDStar, std::array{-kKPlus, +kPiPlus, +kPiPlus}, true, &signDStar, 3);
            if (indexRecDstarPartReco > -1) {
              if (indexRecDplus > -1) { // Ds1 -> D* K0s -> D+  π0 K0s
                flag = sign * DecayTypeMc::Ds1ToDStarK0ToDPlusPi0K0s;
              } else {
                auto indexRecDzeroPartReco = RecoDecay::getMatchedMCRec<false, false, true, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1]}, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &signD0, 2);
                if (indexRecDzeroPartReco > -1) { // Ds1 -> D* K0s -> D0 π+ K0s -> K- π+ π0 π+ K0s
                  flag = sign * DecayTypeMc::Ds1ToDStarK0ToD0PiK0sPart;
                }
              }
            }
          }
        } else { // we look for D* not matched, but all the other ones yes, we check if we only lost the soft pion
          if (indexRecD0 > -1 && indexRecK0 > -1 && indexRecDstar < 0) {
            indexRecResoPartReco = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[3], vecDaughtersReso[4]}, Pdg::kDS1, std::array{+kPiPlus, -kKPlus, +kPiPlus, -kPiPlus}, true, &sign, 3);
            if (indexRecResoPartReco > -1) {
              flag = sign * DecayTypeMc::Ds1ToDStarK0ToD0NoPiK0sPart;
            }
          }
        }
      }
      if (flag != 0) {
        int indexParticle{-1};
        if (indexRecReso > -1) {
          indexParticle = indexRecReso;
        } else if (indexRecResoPartReco > -1) {
          indexParticle = indexRecResoPartReco;
        }
        auto particleReso = particlesMc.iteratorAt(indexParticle);
        origin = RecoDecay::getCharmHadronOrigin(particlesMc, particleReso, false, &idxBhadMothers);
        motherPt = particleReso.pt();
      }
      if (doMcRecQa) {
        if (indexRecReso > -1) {
          if (nPiToMuReso == 0 && nKaToPiReso == 0) {
            registry.fill(HIST("hMassVsPtDs1Matched"), varUtils.ptD, varUtils.invMassKPiPiV0 - varUtils.invMassD);
          }
          if (nPiToMuReso >= 1) {
            registry.fill(HIST("hMassVsPtDs1MatchedPiToMu"), varUtils.ptD, varUtils.invMassKPiPiV0 - varUtils.invMassD);
          }
          if (nKaToPiReso >= 1) {
            registry.fill(HIST("hMassVsPtDs1MatchedKaToPi"), varUtils.ptD, varUtils.invMassKPiPiV0 - varUtils.invMassD);
          }
        }
        if (indexRecD0 > -1) {
          if (nPiToMuD0 == 0 && nKaToPiD0 == 0) {
            registry.fill(HIST("hMassVsPtD0Matched"), varUtils.ptD, varUtils.invMassDdau);
          }
          if (nPiToMuD0 >= 1) {
            registry.fill(HIST("hMassVsPtD0MatchedPiToMu"), varUtils.ptD, varUtils.invMassDdau);
          }
          if (nKaToPiD0 >= 1) {
            registry.fill(HIST("hMassVsPtD0MatchedKaToPi"), varUtils.ptD, varUtils.invMassDdau);
          }
        }
        if (indexRecDstar > -1) {
          if (nPiToMuDstar == 0 && nKaToPiDstar == 0) {
            registry.fill(HIST("hMassVsPtDstarMatched"), varUtils.ptD, varUtils.invMassD - varUtils.invMassDdau);
          }
          if (nPiToMuDstar >= 1) {
            registry.fill(HIST("hMassVsPtDstarMatchedPiToMu"), varUtils.ptD, varUtils.invMassD - varUtils.invMassDdau);
          }
          if (nKaToPiDstar >= 1) {
            registry.fill(HIST("hMassVsPtDstarMatchedKaToPi"), varUtils.ptD, varUtils.invMassD - varUtils.invMassDdau);
          }
        }
        if (indexRecK0 > -1) {
          if (nPiToMuV0 == 0 && nKaToPiV0 == 0) {
            registry.fill(HIST("hMassVsPtK0Matched"), candidateV0.pT, candidateV0.mK0Short);
          }
          if (nPiToMuV0 >= 1) {
            registry.fill(HIST("hMassVsPtK0MatchedPiToMu"), candidateV0.pT, candidateV0.mK0Short);
          }
          if (nKaToPiV0 >= 1) {
            registry.fill(HIST("hMassVsPtK0MatchedKaToPi"), candidateV0.pT, candidateV0.mK0Short);
          }
        }
      }
    } else if constexpr (decChannel == DecayChannel::DplusV0) {
      // Ds2Star → D+ K0 → (π+K-π+) K0s → (π+K-π+)(π+π-)
      indexRecK0 = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particlesMc, std::array{vecDaughtersReso[3], vecDaughtersReso[4]}, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signV0, 2, &nPiToMuV0, &nKaToPiV0);
      indexRecDplus = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2]}, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &signDPlus, 2, &nPiToMuDplus, &nKaToPiDplus);
      if (indexRecK0 > -1 && indexRecDplus > -1) {
        indexRecReso = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]}, Pdg::kDS2Star, std::array{+kPiPlus, -kKPlus, +kPiPlus, +kPiPlus, -kPiPlus}, true, &sign, 3, &nPiToMuReso, &nKaToPiReso);
        if (indexRecReso > -1 && nPiToMuReso == 0 && nKaToPiReso == 0) {
          flag = sign * DecayTypeMc::Ds2StarToDplusK0sToPiKaPiPiPi;
        } else if (indexRecReso > -1 && nPiToMuReso >= 1 && nKaToPiReso == 0) {
          flag = sign * DecayTypeMc::Ds2StarToDplusK0sOneMu;
        } else if (indexRecReso < 0) {
          // Verify partly reconstructed decay Ds1 -> D* K0s -> D+  π0 K0s
          indexRecDstar = RecoDecay::getMatchedMCRec<false, false, true, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2]}, Pdg::kDStar, std::array{-kKPlus, +kPiPlus, +kPiPlus}, true, &signDStar, 2);
          if (indexRecDstar > -1) {
            indexRecReso = RecoDecay::getMatchedMCRec<false, true, true, true, true>(particlesMc, std::array{vecDaughtersReso[0], vecDaughtersReso[1], vecDaughtersReso[2], vecDaughtersReso[3], vecDaughtersReso[4]}, Pdg::kDS1, std::array{+kPiPlus, -kKPlus, +kPiPlus, +kPiPlus, -kPiPlus}, true, &sign, 3);
            if (indexRecReso > -1) {
              flag = sign * DecayTypeMc::Ds1ToDStarK0ToDPlusPi0K0s;
            }
          }
        }
      }
      if (flag != 0) {
        auto particleReso = particlesMc.iteratorAt(indexRecReso);
        origin = RecoDecay::getCharmHadronOrigin(particlesMc, particleReso, false, &idxBhadMothers);
        motherPt = particleReso.pt();
      }
      if (doMcRecQa) {
        if (indexRecReso > -1) {
          if (nPiToMuReso == 0 && nKaToPiReso == 0) {
            registry.fill(HIST("hMassVsPtDs2StarMatched"), varUtils.ptD, varUtils.invMassKPiPiV0 - varUtils.invMassD);
          }
          if (nPiToMuReso >= 1) {
            registry.fill(HIST("hMassVsPtDs2StarMatchedPiToMu"), varUtils.ptD, varUtils.invMassKPiPiV0 - varUtils.invMassD);
          }
          if (nKaToPiReso >= 1) {
            registry.fill(HIST("hMassVsPtDs2StarMatchedKaToPi"), varUtils.ptD, varUtils.invMassKPiPiV0 - varUtils.invMassD);
          }
        }
        if (indexRecDplus > -1) {
          if (nPiToMuDplus == 0 && nKaToPiDplus == 0) {
            registry.fill(HIST("hMassVsPtDplusMatched"), varUtils.ptD, varUtils.invMassD);
          }
          if (nPiToMuDplus >= 1) {
            registry.fill(HIST("hMassVsPtDplusMatchedPiToMu"), varUtils.ptD, varUtils.invMassD);
          }
          if (nKaToPiDplus >= 1) {
            registry.fill(HIST("hMassVsPtDplusMatchedKaToPi"), varUtils.ptD, varUtils.invMassD);
          }
        }
        if (indexRecK0 > -1) {
          if (nPiToMuV0 == 0 && nKaToPiV0 == 0) {
            registry.fill(HIST("hMassVsPtK0Matched"), candidateV0.pT, candidateV0.mK0Short);
          }
          if (nPiToMuV0 >= 1) {
            registry.fill(HIST("hMassVsPtK0MatchedPiToMu"), candidateV0.pT, candidateV0.mK0Short);
          }
          if (nKaToPiV0 >= 1) {
            registry.fill(HIST("hMassVsPtK0MatchedKaToPi"), candidateV0.pT, candidateV0.mK0Short);
          }
        }
      }
    } // DecayChannel::DplusV0
    if (flag != 0) {
      registry.fill(HIST("hMCRecCounter"), flag);
      registry.fill(HIST("hMCRecOrigin"), origin);
    } else {
      if (indexRecK0 > -1) {
        SETBIT(debug, PartialMatchMc::K0Matched);
      }
      if (indexRecD0 > -1) {
        SETBIT(debug, PartialMatchMc::D0Matched);
      }
      if (indexRecDstar > -1) {
        SETBIT(debug, PartialMatchMc::DStarMatched);
      }
      if (indexRecDplus > -1) {
        SETBIT(debug, PartialMatchMc::DPlusMatched);
      }
      registry.fill(HIST("hMCRecDebug"), debug);
    }
    rowHfDV0McRecReduced(indexHfCandCharm, indexCandV0, flag, debug, origin, signD0, motherPt);
  }

  template <bool withMl, bool doMc, uint8_t DecayChannel, typename Coll, typename CCands, typename Tr, typename PParticles, typename BBach>
  void runDataCreation(Coll const& collision,
                       CCands const& candsD,
                       BBach const& bachelors,
                       Tr const&,
                       PParticles const& particlesMc,
                       aod::BCsWithTimestamps const&)
  {
    // helpers for ReducedTables filling
    int indexHfReducedCollision = hfReducedCollision.lastIndex() + 1;
    // std::map where the key is the V0.globalIndex() and
    // the value is the V0 index in the table of the selected v0s
    std::map<int64_t, int64_t> selectedV0s;
    std::map<int64_t, int64_t> selectedTracks;
    bool fillHfReducedCollision = false;
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
    }
    fitter.setBz(bz);
    // loop on D candidates
    for (const auto& candD : candsD) {
      // initialize variables depending on decay channel
      bool fillHfCandD = false;
      std::array<float, 3> pVecD;
      std::array<float, 3> pVecProng2;
      std::array<float, 3> secondaryVertexD;
      std::array<int, 3> prongIdsD;
      int8_t dtype;
      std::array<float, 3> bdtScores;
      std::vector<typename Tr::iterator> charmHadDauTracks{};
      varUtils.ptD = candD.pt();
      if constexpr (DecayChannel == DecayChannel::DstarV0 || DecayChannel == DecayChannel::DstarTrack) {
        if (candD.signSoftPi() > 0) {
          varUtils.invMassD = candD.invMassDstar();
          varUtils.invMassDdau = candD.invMassD0();
        } else {
          varUtils.invMassD = candD.invMassAntiDstar();
          varUtils.invMassDdau = candD.invMassD0Bar();
        }
        pVecD = candD.pVector();
        secondaryVertexD[0] = candD.xSecondaryVertexD0();
        secondaryVertexD[1] = candD.ySecondaryVertexD0();
        secondaryVertexD[2] = candD.zSecondaryVertexD0();
        prongIdsD[0] = candD.prong0Id();
        prongIdsD[1] = candD.prong1Id();
        prongIdsD[2] = candD.prongPiId();
        pVecProng2 = candD.pVecSoftPi();
        charmHadDauTracks.push_back(candD.template prong0_as<Tr>());
        charmHadDauTracks.push_back(candD.template prong1_as<Tr>());
        charmHadDauTracks.push_back(candD.template prongPi_as<Tr>());
        dtype = candD.signSoftPi() * DType::Dstar;
        if constexpr (withMl) {
          std::copy(candD.mlProbDstarToD0Pi().begin(), candD.mlProbDstarToD0Pi().end(), bdtScores.begin());
        }
        registry.fill(HIST("hMassVsPtDstarAll"), candD.pt(), varUtils.invMassD - varUtils.invMassDdau);
      } else if constexpr (DecayChannel == DecayChannel::DplusV0) {
        auto prong0 = candD.template prong0_as<Tr>();
        varUtils.invMassD = hfHelper.invMassDplusToPiKPi(candD);
        pVecD = candD.pVector();
        secondaryVertexD[0] = candD.xSecondaryVertex();
        secondaryVertexD[1] = candD.ySecondaryVertex();
        secondaryVertexD[2] = candD.zSecondaryVertex();
        prongIdsD[0] = candD.prong0Id();
        prongIdsD[1] = candD.prong1Id();
        prongIdsD[2] = candD.prong2Id();
        pVecProng2 = candD.pVectorProng2();
        dtype = static_cast<int8_t>(prong0.sign() * DType::Dplus);
        charmHadDauTracks.push_back(candD.template prong0_as<Tr>());
        charmHadDauTracks.push_back(candD.template prong1_as<Tr>());
        charmHadDauTracks.push_back(candD.template prong2_as<Tr>());
        if constexpr (withMl) {
          std::copy(candD.mlProbDplusToPiKPi().begin(), candD.mlProbDplusToPiKPi().end(), bdtScores.begin());
        }
        registry.fill(HIST("hMassVsPtDplusAll"), candD.pt(), varUtils.invMassD);
      } // else if

      // Get single track variables
      float chi2TpcDauMax = -1.f;
      int nItsClsDauMin = 8, nTpcCrossRowsDauMin = 200;
      for (const auto& charmHadTrack : charmHadDauTracks) {
        if (charmHadTrack.itsNCls() < nItsClsDauMin) {
          nItsClsDauMin = charmHadTrack.itsNCls();
        }
        if (charmHadTrack.tpcNClsCrossedRows() < nTpcCrossRowsDauMin) {
          nTpcCrossRowsDauMin = charmHadTrack.tpcNClsCrossedRows();
        }
        if (charmHadTrack.tpcChi2NCl() > chi2TpcDauMax) {
          chi2TpcDauMax = charmHadTrack.tpcChi2NCl();
        }
      }

      if constexpr (DecayChannel == DecayChannel::DplusV0 || DecayChannel == DecayChannel::DstarV0) {
        // Loop on V0 candidates
        for (const auto& v0 : bachelors) {
          auto trackPos = v0.template posTrack_as<Tr>();
          auto trackNeg = v0.template negTrack_as<Tr>();
          // Apply selsection
          auto v0DauTracks = std::array{trackPos, trackNeg};
          if (!buildAndSelectV0(collision, prongIdsD, v0DauTracks)) {
            continue;
          }
          // Get single track variables
          float chi2TpcDauV0Max = -1.f;
          int nItsClsDauV0Min = 8, nTpcCrossRowsDauV0Min = 200;
          for (const auto& v0Track : v0DauTracks) {
            if (v0Track.itsNCls() < nItsClsDauV0Min) {
              nItsClsDauV0Min = v0Track.itsNCls();
            }
            if (v0Track.tpcNClsCrossedRows() < nTpcCrossRowsDauV0Min) {
              nTpcCrossRowsDauV0Min = v0Track.tpcNClsCrossedRows();
            }
            if (v0Track.tpcChi2NCl() > chi2TpcDauV0Max) {
              chi2TpcDauV0Max = v0Track.tpcChi2NCl();
            }
          }
          // propagate V0 to primary vertex (if enabled)
          if (propagateV0toPV) {
            std::array<float, 3> pVecV0Orig = {candidateV0.mom[0], candidateV0.mom[1], candidateV0.mom[2]};
            gpu::gpustd::array<float, 2> dcaInfo;
            auto trackParK0 = o2::track::TrackPar(candidateV0.pos, pVecV0Orig, 0, true);
            trackParK0.setPID(o2::track::PID::K0);
            trackParK0.setAbsCharge(0);
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParK0, 2.f, matCorr, &dcaInfo);
            getPxPyPz(trackParK0, candidateV0.mom);
          }
          if (TESTBIT(candidateV0.v0Type, K0s)) {
            if constexpr (DecayChannel == DecayChannel::DplusV0) {
              varUtils.invMassKPiPiV0 = RecoDecay::m(std::array{candD.pVectorProng0(), candD.pVectorProng1(), candD.pVectorProng2(), candidateV0.mom}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0});
              // varUtils.ptReso = RecoDecay::pt(std::array{candD.pVectorProng0(), candD.pVectorProng1(), candD.pVectorProng2(), candidateV0.mom});
            } else if (DecayChannel == DecayChannel::DstarV0) {
              // varUtils.ptReso = RecoDecay::pt(std::array{candD.pVectorProng0(), candD.pVectorProng1(), candD.pVecSoftPi(), candidateV0.mom});
              if (candD.signSoftPi() > 0) {
                varUtils.invMassKPiPiV0 = RecoDecay::m(std::array{candD.pVectorProng0(), candD.pVectorProng1(), candD.pVecSoftPi(), candidateV0.mom}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0});
              } else {
                varUtils.invMassKPiPiV0 = RecoDecay::m(std::array{candD.pVectorProng1(), candD.pVectorProng0(), candD.pVecSoftPi(), candidateV0.mom}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0});
              }
            }
            registry.fill(HIST("hMassVsPtK0s"), candidateV0.pT, candidateV0.mK0Short);
            if constexpr (DecayChannel == DecayChannel::DstarV0) {
              registry.fill(HIST("hMassDs1"), varUtils.invMassKPiPiV0 - varUtils.invMassD);
            } else if constexpr (DecayChannel == DecayChannel::DplusV0) {
              registry.fill(HIST("hMassDsStar2"), varUtils.invMassKPiPiV0 - varUtils.invMassD);
            }
          }
          bool isLambda = TESTBIT(candidateV0.v0Type, Lambda);
          bool isAntiLambda = TESTBIT(candidateV0.v0Type, AntiLambda);
          if (isLambda || isAntiLambda) {
            if constexpr (DecayChannel == DecayChannel::DplusV0) {
              varUtils.invMassKPiPiV0 = RecoDecay::m(std::array{candD.pVectorProng0(), candD.pVectorProng1(), candD.pVectorProng2(), candidateV0.mom}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassLambda0});
            } else if (DecayChannel == DecayChannel::DstarV0) {
              if (candD.signSoftPi() > 0) {
                varUtils.invMassKPiPiV0 = RecoDecay::m(std::array{candD.pVectorProng0(), candD.pVectorProng1(), candD.pVecSoftPi(), candidateV0.mom}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassLambda0});
              } else {
                varUtils.invMassKPiPiV0 = RecoDecay::m(std::array{candD.pVectorProng1(), candD.pVectorProng0(), candD.pVecSoftPi(), candidateV0.mom}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassLambda0});
              }
            }
            if (isLambda || isAntiLambda) {
              registry.fill(HIST("hMassVsPtLambda"), candidateV0.pT, candidateV0.mLambda);
            }
            if constexpr (DecayChannel == DecayChannel::DplusV0) {
              registry.fill(HIST("hMassXcRes"), varUtils.invMassKPiPiV0 - varUtils.invMassD);
            }
          }
          // fill V0 table
          // if information on V0 already stored, go to next V0
          if (!selectedV0s.count(v0.globalIndex())) {
            hfCandV0(trackPos.globalIndex(), trackNeg.globalIndex(),
                     indexHfReducedCollision,
                     candidateV0.pos[0], candidateV0.pos[1], candidateV0.pos[2],
                     candidateV0.momPos[0], candidateV0.momPos[1], candidateV0.momPos[2],
                     candidateV0.momNeg[0], candidateV0.momNeg[1], candidateV0.momNeg[2],
                     candidateV0.cosPA,
                     candidateV0.dcaV0ToPv,
                     nItsClsDauV0Min, nTpcCrossRowsDauV0Min, chi2TpcDauV0Max,
                     candidateV0.v0Type);
            selectedV0s[v0.globalIndex()] = hfCandV0.lastIndex();
          }
          fillHfCandD = true;
          // Optional filling of MC Rec table
          if constexpr (doMc) {
            std::vector<typename Tr::iterator> charmResoDauTracks{};
            for (const auto& track : charmHadDauTracks) {
              charmResoDauTracks.push_back(track);
            }
            charmResoDauTracks.push_back(trackPos);
            charmResoDauTracks.push_back(trackNeg);
            int indexHfCandCharm = hfCandD.lastIndex() + 1;
            int indexHfCandV0 = hfCandV0.lastIndex();
            fillMcRecoInfo<DecayChannel>(particlesMc, charmResoDauTracks, indexHfCandCharm, indexHfCandV0);
          }
        } // V0 loop
      } else if constexpr (DecayChannel == DecayChannel::DstarTrack) {
        for (const auto& trackIndex : bachelors) {
          auto track = trackIndex.template track_as<Tr>();
          if (!isTrackSelected(track, prongIdsD)) {
            continue;
          }

          // if the track has been reassociated, re-propagate it to PV (minor difference)
          auto trackParCovTrack = getTrackParCov(track);
          o2::gpu::gpustd::array<float, 2> dcaTrack{track.dcaXY(), track.dcaZ()};
          std::array<float, 3> pVecTrack = track.pVector();
          if (track.collisionId() != collision.globalIndex()) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovTrack, 2.f, matCorr, &dcaTrack);
            getPxPyPz(trackParCovTrack, pVecTrack);
          }

          registry.fill(HIST("hdEdxVsP"), track.p(), track.tpcSignal());
          float invMassKPiPiP{0.f};
          if (candD.signSoftPi() > 0) {
            invMassKPiPiP = RecoDecay::m(std::array{candD.pVectorProng0(), candD.pVectorProng1(), candD.pVecSoftPi(), pVecTrack}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassProton});
          } else {
            invMassKPiPiP = RecoDecay::m(std::array{candD.pVectorProng1(), candD.pVectorProng0(), candD.pVecSoftPi(), pVecTrack}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassProton});
          }
          registry.fill(HIST("hMassDstarProton"), invMassKPiPiP - varUtils.invMassD);
          if (!selectedTracks.count(track.globalIndex())) {
            hfTrackNoParam(indexHfReducedCollision,
                           track.px(), track.py(), track.pz(), track.sign(),
                           track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                           track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                           track.hasTOF(), track.hasTPC(), track.itsNCls(), track.tpcNClsCrossedRows(), track.tpcChi2NCl());
            selectedTracks[track.globalIndex()] = hfTrackNoParam.lastIndex();
          }
          fillHfCandD = true;
        } // track loop
      }

      if (fillHfCandD) { // fill candDplus table only once per D candidate, only if at least one V0 is found
        hfCandD(prongIdsD[0], prongIdsD[1], prongIdsD[2],
                indexHfReducedCollision,
                secondaryVertexD[0], secondaryVertexD[1], secondaryVertexD[2],
                candD.pxProng0(), candD.pyProng0(), candD.pzProng0(),
                candD.pxProng1(), candD.pyProng1(), candD.pzProng1(),
                pVecProng2[0], pVecProng2[1], pVecProng2[2],
                nItsClsDauMin, nTpcCrossRowsDauMin, chi2TpcDauMax, dtype);
        if constexpr (withMl) {
          hfCandDMl(bdtScores[0], bdtScores[1], bdtScores[2], -1., -1., -1.);
        }
        fillHfReducedCollision = true;
        if constexpr (DecayChannel == DecayChannel::DstarV0 || DecayChannel == DecayChannel::DstarTrack) {
          registry.fill(HIST("hMassVsPtDstarPaired"), candD.pt(), varUtils.invMassD - varUtils.invMassDdau);
        } else if constexpr (DecayChannel == DecayChannel::DplusV0) {
          registry.fill(HIST("hMassVsPtDplusPaired"), candD.pt(), varUtils.invMassD);
        }
        registry.fill(HIST("hDType"), dtype);
      }
    } // candsD loop
    registry.fill(HIST("hEvents"), 1 + Event::Processed);
    if (!fillHfReducedCollision) {
      registry.fill(HIST("hEvents"), 1 + Event::NoDV0Selected);
      return;
    }
    registry.fill(HIST("hEvents"), 1 + Event::DV0Selected);
    float centrality = -1.f;
    uint16_t hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
    // fill collision table if it contains a DPi pair a minima
    hfReducedCollision(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), hfRejMap, bz);
  } // run data creation

  template <uint8_t decayChannel>
  void runMcGen(aod::McParticles const& particlesMc)
  {
    // Match generated particles.
    for (const auto& particle : particlesMc) {
      int8_t sign{0};
      int8_t flag{0};
      int8_t signDStar{0};
      int8_t signDPlus{0};
      int8_t signV0{0};
      int8_t origin = 0;
      std::vector<int> idxBhadMothers{};

      if constexpr (decayChannel == DecayChannel::DstarV0) {
        // Ds1 → D* K0
        if (RecoDecay::isMatchedMCGen<false, true>(particlesMc, particle, Pdg::kDS1, std::array{static_cast<int>(Pdg::kDStar), +kK0}, true, &sign, 1)) {
          registry.fill(HIST("hMCSignCounter"), sign);
          origin = RecoDecay::getCharmHadronOrigin(particlesMc, particle, false, &idxBhadMothers);
          registry.fill(HIST("hMCGenOrigin"), origin);
          auto candV0MC = particlesMc.rawIteratorAt(particle.daughtersIds().back());
          auto candDStarMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          // K0 -> K0s -> π+π-
          if (RecoDecay::isMatchedMCGen<false, true>(particlesMc, candV0MC, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signV0, 2)) {
            // D* -> D0 π+ -> K-π+π+
            if (RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{static_cast<int>(Pdg::kD0), +static_cast<int>(kPiPlus)}, true, &signDStar, 1)) {
              auto candD0MC = particlesMc.rawIteratorAt(candDStarMC.daughtersIds().front());
              if (RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{-kKPlus, +kPiPlus, +kPiPlus}, true, &signDStar, 2)) {
                flag = signDStar * DecayTypeMc::Ds1ToDStarK0ToD0PiK0s;
              } else if (RecoDecay::isMatchedMCGen(particlesMc, candD0MC, Pdg::kD0, std::array{-kKPlus, +kPiPlus, +kPiPlus, +kPi0}, true, &signDStar, 2) ||
                         RecoDecay::isMatchedMCGen(particlesMc, candD0MC, Pdg::kD0, std::array{-kKPlus, +kPiPlus, +kPiPlus, -kPi0}, true, &signDStar, 2)) {
                flag = signDStar * DecayTypeMc::Ds1ToDStarK0ToD0PiK0sPart;
              }
            } else if (RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{static_cast<int>(Pdg::kDPlus), static_cast<int>(kGamma)}, true, &signDStar, 1) ||
                       RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{static_cast<int>(Pdg::kDPlus), -static_cast<int>(kGamma)}, true, &signDStar, 1) ||
                       RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{static_cast<int>(Pdg::kDPlus), static_cast<int>(kPi0)}, true, &signDStar, 1) ||
                       RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{static_cast<int>(Pdg::kDPlus), -static_cast<int>(kPi0)}, true, &signDStar, 1)) {
              auto candDPlusMC = particlesMc.rawIteratorAt(candDStarMC.daughtersIds().front());
              if (RecoDecay::isMatchedMCGen(particlesMc, candDPlusMC, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &signDPlus, 2))
                flag = sign * DecayTypeMc::Ds1ToDStarK0ToDPlusPi0K0s;
            }
          }
        } else {
          if (std::abs(particle.pdgCode()) == Pdg::kDS1) {
            origin = RecoDecay::getCharmHadronOrigin(particlesMc, particle, false, &idxBhadMothers);
            registry.fill(HIST("hMCOriginCounterWrongDecay"), origin);
          }
        }
        // save information for task
        if (flag == 0) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), MassDS1);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs;
        std::array<float, 2> yProngs;
        std::array<float, 2> etaProngs;
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        registry.fill(HIST("hMCGenCounter"), flag, ptParticle);
        rowHfResoMcGenReduced(flag, origin, ptParticle, yParticle, etaParticle,
                              ptProngs[0], yProngs[0], etaProngs[0],
                              ptProngs[1], yProngs[1], etaProngs[1]);
      } else if constexpr (decayChannel == DecayChannel::DplusV0) { // Ds2Star → D+ K0
        if (RecoDecay::isMatchedMCGen<false, true>(particlesMc, particle, Pdg::kDS2Star, std::array{static_cast<int>(Pdg::kDPlus), +kK0}, true, &sign, 1)) {
          registry.fill(HIST("hMCSignCounter"), sign);
          origin = RecoDecay::getCharmHadronOrigin(particlesMc, particle, false, &idxBhadMothers);
          registry.fill(HIST("hMCGenOrigin"), origin);
          auto candV0MC = particlesMc.rawIteratorAt(particle.daughtersIds().back());
          auto candDPlusMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          // K0 -> K0s -> π+π-
          if (RecoDecay::isMatchedMCGen<false, true>(particlesMc, candV0MC, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signV0, 2)) {
            // D* -> D0 π+ -> K-π+π+
            if (RecoDecay::isMatchedMCGen(particlesMc, candDPlusMC, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &signDPlus, 2)) {
              flag = sign * DecayTypeMc::Ds2StarToDplusK0sToPiKaPiPiPi;
            }
          }
        } else if (RecoDecay::isMatchedMCGen<false, true>(particlesMc, particle, Pdg::kDS1, std::array{static_cast<int>(Pdg::kDStar), +kK0}, true, &sign, 1)) {
          auto candV0MC = particlesMc.rawIteratorAt(particle.daughtersIds().back());
          // K0 -> K0s -> π+π-
          if (RecoDecay::isMatchedMCGen<false, true>(particlesMc, candV0MC, kK0, std::array{+kPiPlus, -kPiPlus}, true, &signV0, 2)) {
            auto candDStarMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
            // D* -> D+ π0/γ ->π+K-π+ π0/γ
            if (RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{static_cast<int>(Pdg::kDPlus), static_cast<int>(kGamma)}, true, &signDStar, 1) ||
                RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{static_cast<int>(Pdg::kDPlus), -static_cast<int>(kGamma)}, true, &signDStar, 1) ||
                RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{static_cast<int>(Pdg::kDPlus), static_cast<int>(kPi0)}, true, &signDStar, 1) ||
                RecoDecay::isMatchedMCGen(particlesMc, candDStarMC, Pdg::kDStar, std::array{static_cast<int>(Pdg::kDPlus), -static_cast<int>(kPi0)}, true, &signDStar, 1)) {
              auto candDPlusMC = particlesMc.rawIteratorAt(candDStarMC.daughtersIds().front());
              if (RecoDecay::isMatchedMCGen(particlesMc, candDPlusMC, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &signDPlus, 2))
                flag = sign * DecayTypeMc::Ds1ToDStarK0ToDPlusPi0K0s;
            }
          }
        } else {
          if (std::abs(particle.pdgCode()) == Pdg::kDS2Star) {
            origin = RecoDecay::getCharmHadronOrigin(particlesMc, particle, false, &idxBhadMothers);
            // LOGF(info, "Found DS2Star that decays into %d, %d", particlesMc.rawIteratorAt(particle.daughtersIds().front()).pdgCode(),particlesMc.rawIteratorAt(particle.daughtersIds().back()).pdgCode());
            registry.fill(HIST("hMCOriginCounterWrongDecay"), origin);
          }
        }
        // save information for task
        if (flag == 0) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), MassDS2Star);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs;
        std::array<float, 2> yProngs;
        std::array<float, 2> etaProngs;
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        registry.fill(HIST("hMCGenCounter"), flag, ptParticle);
        rowHfResoMcGenReduced(flag, origin, ptParticle, yParticle, etaParticle,
                              ptProngs[0], yProngs[0], etaProngs[0],
                              ptProngs[1], yProngs[1], etaProngs[1]);
      } // Dplus V0
    } // for loop
  } // gen

  void processDplusV0(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                      CandsDplusFiltered const& candsDplus,
                      aod::V0s const& V0s,
                      TracksIUWithPID const& tracks,
                      aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::DplusV0>(collision, candsDThisColl, V0sThisColl, tracks, tracks, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDplusV0, "Process Dplus candidates paired with V0s without MC info and without ML info", true);

  void processDplusV0MC(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                        CandsDplusFiltered const& candsDplus,
                        aod::V0s const& V0s,
                        TracksIUWithPIDAndMC const& tracks,
                        aod::McParticles const& particlesMc,
                        aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::DplusV0>(collision, candsDThisColl, V0sThisColl, tracks, particlesMc, bcs);
    }
    runMcGen<DecayChannel::DplusV0>(particlesMc);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDplusV0MC, "Process DPlus candidates paired with V0s with MC matching and without ML info", false);

  void processDplusV0WithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                            CandsDplusFilteredWithMl const& candsDplus,
                            aod::V0s const& V0s,
                            TracksIUWithPID const& tracks,
                            aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<true, false, DecayChannel::DplusV0>(collision, candsDThisColl, V0sThisColl, tracks, tracks, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDplusV0WithMl, "Process Dplus candidates paired with V0s with ML info", false);

  void processDplusV0MCWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                              CandsDplusFilteredWithMl const& candsDplus,
                              aod::V0s const& V0s,
                              TracksIUWithPIDAndMC const& tracks,
                              aod::McParticles const& particlesMc,
                              aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<true, true, DecayChannel::DplusV0>(collision, candsDThisColl, V0sThisColl, tracks, particlesMc, bcs);
    }
    runMcGen<DecayChannel::DplusV0>(particlesMc);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDplusV0MCWithMl, "Process DPlus candidates paired with V0s with MC matching and with ML info", false);

  void processDstarV0(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                      CandDstarFiltered const& candsDstar,
                      aod::V0s const& V0s,
                      TracksIUWithPID const& tracks,
                      aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::DstarV0>(collision, candsDThisColl, V0sThisColl, tracks, tracks, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDstarV0, "Process DStar candidates paired with V0s without MC info and without ML info", false);

  void processDstarV0MC(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                        CandDstarFiltered const& candsDstar,
                        aod::V0s const& V0s,
                        TracksIUWithPIDAndMC const& tracks,
                        aod::McParticles const& particlesMc,
                        aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<false, true, DecayChannel::DstarV0>(collision, candsDThisColl, V0sThisColl, tracks, particlesMc, bcs);
    }
    runMcGen<DecayChannel::DstarV0>(particlesMc);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDstarV0MC, "Process DStar candidates paired with V0s with MC matching and without ML info", false);

  void processDstarV0WithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                            CandDstarFilteredWithMl const& candsDstar,
                            aod::V0s const& V0s,
                            TracksIUWithPID const& tracks,
                            aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollisionWithMl, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<true, false, DecayChannel::DstarV0>(collision, candsDThisColl, V0sThisColl, tracks, tracks, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDstarV0WithMl, "Process DStar candidates paired with V0s with ML info", false);

  void processDstarV0MCWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                              CandDstarFilteredWithMl const& candsDstar,
                              aod::V0s const& V0s,
                              TracksIUWithPIDAndMC const& tracks,
                              aod::McParticles const& particlesMc,
                              aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<true, true, DecayChannel::DstarV0>(collision, candsDThisColl, V0sThisColl, tracks, particlesMc, bcs);
    }
    runMcGen<DecayChannel::DstarV0>(particlesMc);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDstarV0MCWithMl, "Process MC DStar candidates paired with V0s with ML info", false);

  void processDstarTrack(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                         CandDstarFiltered const& candsDstar,
                         aod::TrackAssoc const& trackIndices,
                         TracksWithPID const& tracks,
                         aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DecayChannel::DstarTrack>(collision, candsDThisColl, trackIdsThisColl, tracks, tracks, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDstarTrack, "Process DStar candidates paired with tracks without MC info and without ML info", false);

  void processDstarTrackWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                               CandDstarFilteredWithMl const& candsDstar,
                               aod::TrackAssoc const& trackIndices,
                               TracksWithPID const& tracks,
                               aod::BCsWithTimestamps const& bcs)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollisionWithMl, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<true, false, DecayChannel::DstarTrack>(collision, candsDThisColl, trackIdsThisColl, tracks, tracks, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoReduced, processDstarTrackWithMl, "Process DStar candidates paired with tracks with ML info", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorCharmResoReduced>(cfgc)};
}
