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

/// \file piKpRAA.cxx

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/ZDCConstants.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h" // required for Filter op.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include "TPDGCode.h"
#include "TVector3.h"
#include <TString.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <string_view>
#include <typeinfo>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;
using namespace o2::constants::math;
using namespace o2::framework::expressions;

using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::FT0MultZeqs, o2::aod::CentFT0Cs, aod::TPCMults, o2::aod::BarrelMults>;
using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection, aod::TracksCovIU, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCEl, aod::pidTOFEl, aod::pidTPCMu, aod::pidTOFMu>;

struct piKpRAA {

  static constexpr float kZero{0.};
  static constexpr float kOne{1.};
  static constexpr float kMinCharge{3.f};

  static constexpr float DefaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {DefaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  struct : ConfigurableGroup {
    Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};

    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
    Configurable<float> minEtaDaughter{"minEtaDaughter", -0.8, "Daughter minimum-eta selection"};
    Configurable<float> maxEtaDaughter{"maxEtaDaughter", +0.8, "Daughter maximum-eta selection"};
    Configurable<float> minPt{"minPt", 0.15, "minimum pt of the tracks"};
    Configurable<float> maxPt{"maxPt", 20.0, "maximum pt of the tracks"};

    // Standard 5 topological criteria
    Configurable<float> v0cospa{"v0cospa", 0.995, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcanegtopv{"dcanegtopv", .1, "min DCA Neg To PV (cm)"};
    Configurable<float> dcapostopv{"dcapostopv", .1, "min DCA Pos To PV (cm)"};
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};

    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

    // Selection
    Configurable<bool> applyInvMassSel{"applyInvMassSel", true, "Select V0s close to the Inv. mass value"};
    Configurable<float> minMassK0s{"minMassK0", 0.4f, "Min. Inv. Mass K0"};
    Configurable<float> maxMassK0s{"maxMassK0", 0.6f, "Max. Inv. Mass K0"};
    Configurable<float> minMassLambda{"minMassLambda", 1.1f, "Min. Inv. Mass Lambda"};
    Configurable<float> maxMassLambda{"maxMassLambda", 1.2f, "Max. Inv. Mass Lambda"};
    Configurable<float> minMassGamma{"minMassGamma", 0.000922f, "Min. Inv. Mass Gamma"};
    Configurable<float> maxMassGamma{"maxMassGamma", 0.002022f, "Max. Inv. Mass Gamma"};

    // PID (TPC/TOF)
    Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};

  } v0Selections;

  // Configurables Event Selection
  Configurable<bool> isNoCollInTimeRangeStrict{"isNoCollInTimeRangeStrict", true, "use isNoCollInTimeRangeStrict?"};
  Configurable<bool> isNoCollInTimeRangeStandard{"isNoCollInTimeRangeStandard", false, "use isNoCollInTimeRangeStandard?"};
  Configurable<bool> isNoCollInRofStrict{"isNoCollInRofStrict", true, "use isNoCollInRofStrict?"};
  Configurable<bool> isNoCollInRofStandard{"isNoCollInRofStandard", false, "use isNoCollInRofStandard?"};
  Configurable<bool> isNoHighMultCollInPrevRof{"isNoHighMultCollInPrevRof", true, "use isNoHighMultCollInPrevRof?"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", false, "use isNoCollInTimeRangeNarrow?"};
  Configurable<bool> isOccupancyCut{"isOccupancyCut", true, "Occupancy cut?"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", false, "T0C Occu cut"};
  Configurable<bool> useMidRapNchSel{"useMidRapNchSel", true, "Use mid-rapidit Nch selection"};
  Configurable<bool> skipRecoColGTOne{"skipRecoColGTOne", true, "Remove collisions if reconstructed more than once"};
  Configurable<std::string> detector4Calibration{"detector4Calibration", "T0M", "Detector for nSigma-Nch rejection"};

  // Event selection
  Configurable<float> posZcut{"posZcut", +10.0, "z-vertex position cut"};
  Configurable<float> minT0CcentCut{"minT0CcentCut", 0.0, "Min T0C Cent. cut"};
  Configurable<float> maxT0CcentCut{"maxT0CcentCut", 100.0, "Max T0C Cent. cut"};
  Configurable<float> minOccCut{"minOccCut", 0., "min Occu cut"};
  Configurable<float> maxOccCut{"maxOccCut", 500., "max Occu cut"};

  ConfigurableAxis binsPtPhiCut{"binsPtPhiCut", {VARIABLE_WIDTH, 0.0, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 45.0, 50.0}, "pT"};
  ConfigurableAxis binsPtV0s{"binsPtV0s", {VARIABLE_WIDTH, 0.0, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 9.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0}, "pT"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12}, "pT binning"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "T0C binning"};
  ConfigurableAxis axisArmAlpha{"axisArmAlpha", {200, -1.0, 1.0}, "Armenteros alpha"};
  ConfigurableAxis axisArmqT{"axisArmqT", {200, 0.0f, 0.3f}, "Armenteros qT"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, "Mass K0Short"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, "Mass Lambda"};
  ConfigurableAxis axisGammaMass{"axisGammaMass", {200, 0.0f, 0.5f}, "Mass Gamma"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10.0f, 10.0f}, "N sigma TPC"};
  ConfigurableAxis axisdEdx{"axisdEdx", {140, 20.0, 160.0}, "dEdx binning"};

  // CCDB paths
  Configurable<std::string> paTHEff{"paTHEff", "Users/o/omvazque/MCcorrection/perTimeStamp/TrackingEff", "base path to the ccdb object"};
  Configurable<std::string> paTHFD{"paTHFD", "Users/o/omvazque/MCcorrection/perTimeStamp/FeedDown", "base path to the ccdb object"};
  Configurable<std::string> paTHmeanNch{"paTHmeanNch", "Users/o/omvazque/FitMeanNch_9May2025", "base path to the ccdb object"};
  Configurable<std::string> paTHsigmaNch{"paTHsigmaNch", "Users/o/omvazque/FitSigmaNch_9May2025", "base path to the ccdb object"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  enum EvCutLabel {
    All = 1,
    SelEigth,
    NoSameBunchPileup,
    IsGoodZvtxFT0vsPV,
    NoCollInTimeRangeStrict,
    NoCollInTimeRangeStandard,
    NoCollInRofStrict,
    NoCollInRofStandard,
    NoHighMultCollInPrevRof,
    NoCollInTimeRangeNarrow,
    OccuCut,
    Centrality,
    VtxZ,
    Zdc,
    TZero,
    Tdc,
    Zem
  };

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH2F* hEfficiency = nullptr;
    TH2F* hFeedDown = nullptr;
    bool correctionsLoaded = false;
  } cfg;

  struct NchConfig {
    TH1F* hMeanNch = nullptr;
    TH1F* hSigmaNch = nullptr;
    bool calibrationsLoaded = false;
  } cfgNch;

  TrackSelection trkSelDaugthers;
  TrackSelection trkSelGlobal;
  TrackSelection trkSelDaugthersV0s()
  {
    TrackSelection selectedTracks;
    selectedTracks.SetEtaRange(-0.8f, 0.8f);
    selectedTracks.SetMinNCrossedRowsTPC(70);
    return selectedTracks;
  }

  void
    init(InitContext const&)
  {

    trkSelDaugthers = trkSelDaugthersV0s();
    trkSelGlobal = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default);

    // define axes you want to use
    const AxisSpec axisZpos{48, -12., 12., "Vtx_{z} (cm)"};
    const AxisSpec axisEvent{14, 0.5, 14.5, ""};
    const AxisSpec axisEta{40, -1., +1., "#eta"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtPhiCut{binsPtPhiCut, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtV0s{binsPtV0s, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisCent{binsCent, "T0C centrality"};

    registry.add("EventCounter", ";;Events", kTH1F, {axisEvent});

    auto hstat = registry.get<TH1>(HIST("EventCounter"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "SelEigth");
    x->SetBinLabel(3, "NoSameBunchPileup");
    x->SetBinLabel(4, "GoodZvtxFT0vsPV");
    x->SetBinLabel(5, "NoCollInTimeRangeStrict");
    x->SetBinLabel(6, "NoCollInTimeRangeStandard");
    x->SetBinLabel(7, "NoCollInRofStrict");
    x->SetBinLabel(8, "NoCollInRofStandard");
    x->SetBinLabel(9, "NoHighMultCollInPrevRof");
    x->SetBinLabel(10, "NoCollInTimeRangeNarrow");
    x->SetBinLabel(11, "Occupancy Cut");
    x->SetBinLabel(12, "Cent. Sel.");
    x->SetBinLabel(13, "VtxZ cut");

    if (doprocessCalibrationAndV0s) {
      registry.add("zPos", ";;Entries;", kTH1F, {axisZpos});
      registry.add("T0Ccent", ";;Entries", kTH1F, {axisCent});

      registry.add("dcaDauVsPt", ";V0 #it{p}_{T} (GeV/#it{c});DCA_{xy} (cm) daughters;", kTH2F, {{{axisPtV0s}, {200, -10., 10.}}});
      registry.add("nSigmaPiFromK0", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigmaPiFromLambda", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigmaPrFromLambda", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigmaPiFromAntiLambda", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigmaPrFromAntiLambda", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigmaElFromGammaConv", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("ArmAll", "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmAfterTopoSel", "Armenteros-Podolanski anfter topological selection;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmK0NOSel", "Armenteros-Podolanski WITH OUT 5 #times q_{T} > #alpha selection;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmK0", "Armenteros-Podolanski WITH 5 #times q_{T} > #alpha selection;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmLambda", "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmAntiLambda", "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmGamma", "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("MassK0ShortVsPt", ";;Inv. Mass (GeV/#it{c}^{2});", kTH2F, {axisPtV0s, axisK0Mass});
      registry.add("MassLambdaVsPt", ";;Inv. Mass (GeV/#it{c}^{2});", kTH2F, {axisPtV0s, axisLambdaMass});
      registry.add("MassAntiLambdaVsPt", ";;Inv. Mass (GeV/#it{c}^{2});", kTH2F, {axisPtV0s, axisLambdaMass});
      registry.add("MassGammaVsPt", ";;Inv. Mass (GeV/#it{c}^{2});", kTH2F, {axisPtV0s, axisGammaMass});
      registry.add("dEdxPiPos", ";Momentum (GeV/#it{c}); dE/dx", kTH2F, {axisPtV0s, axisdEdx});
      registry.add("dEdxPiNeg", ";Momentum (GeV/#it{c}); dE/dx", kTH2F, {axisPtV0s, axisdEdx});
      registry.add("dEdxPrPos", ";Momentum (GeV/#it{c}); dE/dx", kTH2F, {axisPtV0s, axisdEdx});
      registry.add("dEdxPrNeg", ";Momentum (GeV/#it{c}); dE/dx", kTH2F, {axisPtV0s, axisdEdx});
      registry.add("dEdxEl", ";Momentum (GeV/#it{c}); dE/dx", kTH2F, {axisPtV0s, axisdEdx});

      registry.add("NclFindableVsPt", ";#it{p}_{T} (GeV/#it{c}); Ncl TPC", kTH2F, {{{axisPtV0s}, {160, 0, 160}}});
      registry.add("NclFindableVsPtp", ";#it{p}_{T} (GeV/#it{c}); #LTNcl#GT TPC", kTProfile, {axisPtV0s});
      registry.add("NclFoundVsPt", ";#it{p}_{T} (GeV/#it{c}); Ncl TPC", kTH2F, {{{axisPtV0s}, {160, 0, 160}}});
      registry.add("NclFoundVsPtp", ";#it{p}_{T} (GeV/#it{c}); #LTNcl#GT TPC", kTProfile, {axisPtV0s});
      registry.add("NclFindableMinusFoundVsPt", ";#it{p}_{T} (GeV/#it{c}); Ncl TPC", kTH2F, {{{axisPtV0s}, {21, -0.5, 20.5}}});
      registry.add("NclFindableMinusFoundVsPtp", ";#it{p}_{T} (GeV/#it{c}); #LTNcl#GT TPC", kTProfile, {axisPtV0s});
      registry.add("NclFoundVsPtVsPhip", "", kTProfile2D, {{{axisPtPhiCut}, {350, 0.0, 0.35}}});
      registry.add("NclVsEta", ";#eta; Ncl TPC", kTH2F, {{{axisEta}, {160, 0, 160}}});
      registry.add("dEdxVsMomentum", ";Momentum (GeV/#it{c}); dE/dx", kTH2F, {axisPtV0s, axisdEdx});
      registry.add("dEdxVsEtaPiMIP", ";#eta; dE/dx MIP Pions (0.4 < #it{p} < 0.6 GeV/#it{c})", kTH2F, {{{axisEta}, {50, 20, 70}}});
      registry.add("dEdxVsEtaPiMIPp", ";#eta; #LTdE/dx#GT MIP Pions", kTProfile, {axisEta});
    }

    LOG(info) << "\tccdbNoLaterThan=" << ccdbNoLaterThan.value;
    LOG(info) << "\tpaTHEff=" << paTHEff.value;
    LOG(info) << "\tpaTHFD=" << paTHFD.value;
    LOG(info) << "\tuseMidRapNchSel=" << useMidRapNchSel.value;
    LOG(info) << "\tdetector4Calibration=" << detector4Calibration.value;
    LOG(info) << "\tpaTHmeanNch=" << paTHmeanNch.value;
    LOG(info) << "\tpaTHsigmaNch=" << paTHsigmaNch.value;
    LOG(info) << "\tminPt=" << v0Selections.minPt;
    LOG(info) << "\tmaxPt=" << v0Selections.maxPt;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
  }

  void processCalibrationAndV0s(ColEvSels::iterator const& collision, BCsRun3 const& /**/, aod::V0Datas const& v0s, aod::FV0As const& /**/, aod::FT0s const& /**/, TracksFull const& tracks)
  {
    // LOG(info) << " Collisions size: " << collisions.size() << "
    // Table's size: " << collisions.tableSize() << "\n";
    // LOG(info) << "Run number: " << foundBC.runNumber() << "\n";
    if (!isEventSelected(collision)) {
      return;
    }

    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    uint64_t timeStamp{foundBC.timestamp()};
    const int magField{getMagneticField(timeStamp)};

    registry.fill(HIST("zPos"), collision.posZ());
    registry.fill(HIST("T0Ccent"), collision.centFT0C());

    for (const auto& track : tracks) {

      if (!trkSelGlobal.IsSelected(track))
        continue;

      if (track.pt() < v0Selections.minPt || track.pt() > v0Selections.maxPt)
        continue;

      const float momentum{track.p()};
      const float pt{track.pt()};
      const float phi{track.phi()};
      const float eta{track.eta()};
      const float dedx{track.tpcSignal()};

      if (momentum > 0.4 && momentum < 0.6 && dedx < 65.0) {
        registry.fill(HIST("dEdxVsEtaPiMIP"), eta, dedx);
        registry.fill(HIST("dEdxVsEtaPiMIPp"), eta, dedx);
      }

      registry.fill(HIST("NclFindableVsPt"), pt, track.tpcNClsFindable());
      registry.fill(HIST("NclFindableVsPtp"), pt, track.tpcNClsFindable());
      registry.fill(HIST("NclFoundVsPt"), pt, track.tpcNClsFound());
      registry.fill(HIST("NclFoundVsPtp"), pt, track.tpcNClsFound());
      registry.fill(HIST("NclFindableMinusFoundVsPt"), pt, track.tpcNClsFindableMinusFound());
      registry.fill(HIST("NclFindableMinusFoundVsPtp"), pt, track.tpcNClsFindableMinusFound());
      registry.fill(HIST("NclVsEta"), eta, track.tpcNClsFound());

      float phiPrime{phi};
      const int charge{track.sign()};
      PhiPrime(phiPrime, magField, charge);
      registry.fill(HIST("NclFoundVsPtVsPhip"), pt, phiPrime, track.tpcNClsFound());

      registry.fill(HIST("dEdxVsMomentum"), momentum, dedx);
    }

    for (const auto& v0 : v0s) {

      // Select V0 type
      if (v0.v0Type() != v0Selections.v0TypeSelection)
        continue;

      // Positive-(negative-)charged tracks (daughters)
      const auto& posTrack = v0.posTrack_as<TracksFull>();
      const auto& negTrack = v0.negTrack_as<TracksFull>();

      if (posTrack.sign() == negTrack.sign())
        continue;

      if (!(passesTrackSelectionDaughters(posTrack) && passesTrackSelectionDaughters(negTrack)))
        continue;

      const TVector3 ppos(posTrack.px(), posTrack.py(), posTrack.pz());
      const TVector3 pneg(negTrack.px(), negTrack.py(), negTrack.pz());
      double alpha, qT;

      GetArmeterosVariables(ppos, pneg, alpha, qT);
      registry.fill(HIST("ArmAll"), alpha, qT);

      if (!passesV0TopologicalSelection(v0))
        continue;

      registry.fill(HIST("ArmAfterTopoSel"), alpha, qT);
      registry.fill(HIST("dcaDauVsPt"), v0.pt(), v0.dcapostopv());
      registry.fill(HIST("dcaDauVsPt"), v0.pt(), v0.dcanegtopv());

      if (passesK0Selection(collision, v0)) {
        registry.fill(HIST("nSigmaPiFromK0"), posTrack.pt(), posTrack.tpcNSigmaPi());
        registry.fill(HIST("nSigmaPiFromK0"), negTrack.pt(), negTrack.tpcNSigmaPi());
        registry.fill(HIST("ArmK0NOSel"), alpha, qT);
        if (v0Selections.armPodCut * qT > std::abs(alpha)) {
          if (v0Selections.applyInvMassSel) {
            if (!(v0.mK0Short() > v0Selections.minMassK0s && v0.mK0Short() < v0Selections.maxMassK0s))
              continue;
          }
          registry.fill(HIST("ArmK0"), alpha, qT);
          registry.fill(HIST("MassK0ShortVsPt"), v0.pt(), v0.mK0Short());

          registry.fill(HIST("dEdxPiPos"), posTrack.p(), posTrack.tpcSignal());
          registry.fill(HIST("dEdxPiNeg"), negTrack.p(), negTrack.tpcSignal());
        }
      }

      if (passesLambdaSelection(collision, v0)) {
        if (v0Selections.applyInvMassSel) {
          if (!(v0.mLambda() > v0Selections.minMassLambda && v0.mLambda() < v0Selections.maxMassLambda))
            continue;
        }
        registry.fill(HIST("nSigmaPrFromLambda"), posTrack.pt(), posTrack.tpcNSigmaPr());
        registry.fill(HIST("nSigmaPiFromLambda"), negTrack.pt(), negTrack.tpcNSigmaPi());
        registry.fill(HIST("ArmLambda"), alpha, qT);
        registry.fill(HIST("MassLambdaVsPt"), v0.pt(), v0.mLambda());
        registry.fill(HIST("dEdxPiNeg"), negTrack.p(), negTrack.tpcSignal());
        registry.fill(HIST("dEdxPrPos"), posTrack.p(), posTrack.tpcSignal());
        // std::cout << "pos charge = " << posTrack.sign() << " | neg charge = " << negTrack.sign() << '\n';
      }

      if (passesAntiLambdaSelection(collision, v0)) {
        if (v0Selections.applyInvMassSel) {
          if (!(v0.mAntiLambda() > v0Selections.minMassLambda && v0.mAntiLambda() < v0Selections.maxMassLambda))
            continue;
        }
        registry.fill(HIST("nSigmaPrFromAntiLambda"), posTrack.pt(), posTrack.tpcNSigmaPi());
        registry.fill(HIST("nSigmaPiFromAntiLambda"), negTrack.pt(), negTrack.tpcNSigmaPr());
        registry.fill(HIST("ArmAntiLambda"), alpha, qT);
        registry.fill(HIST("MassAntiLambdaVsPt"), v0.pt(), v0.mAntiLambda());

        registry.fill(HIST("dEdxPiPos"), posTrack.p(), posTrack.tpcSignal());
        registry.fill(HIST("dEdxPrNeg"), negTrack.p(), negTrack.tpcSignal());
      }
      if (passesGammaSelection(collision, v0)) {

        // const double ePos{static_cast<double>(posTrack.energy(o2::constants::physics::MassElectron))};
        // const double eEle{static_cast<double>(negTrack.energy(o2::constants::physics::MassElectron))};
        // const float px[2] = {posTrack.px(), negTrack.px()};
        // const float py[2] = {posTrack.py(), negTrack.py()};
        // const float pz[2] = {posTrack.pz(), negTrack.pz()};
        // const double invMass = std::sqrt(std::pow(ePos + eEle, 2.0) - (std::pow(px[0] + px[1], 2.0) + std::pow(py[0] + py[1], 2.0) + std::pow(pz[0] + pz[1], 2.0)));

        if (v0Selections.applyInvMassSel) {
          if (!(v0.mGamma() > v0Selections.minMassGamma && v0.mGamma() < v0Selections.maxMassGamma))
            continue;
        }

        registry.fill(HIST("MassGammaVsPt"), v0.pt(), v0.mGamma());
        registry.fill(HIST("dEdxEl"), posTrack.p(), posTrack.tpcSignal());
        registry.fill(HIST("dEdxEl"), negTrack.p(), negTrack.tpcSignal());

        registry.fill(HIST("nSigmaElFromGammaConv"), negTrack.pt(), negTrack.tpcNSigmaEl());
        registry.fill(HIST("nSigmaElFromGammaConv"), posTrack.pt(), posTrack.tpcNSigmaEl());
        registry.fill(HIST("ArmGamma"), alpha, qT);
      }
    }
  }
  PROCESS_SWITCH(piKpRAA, processCalibrationAndV0s, "Process QA", true);

  template <typename T, typename U>
  void GetArmeterosVariables(const T& ppos, const T& pneg, U& alpha, U& qT)
  {

    alpha = 0., qT = 0.;
    TVector3 pV0 = ppos + pneg;
    double pV0mag = pV0.Mag();
    if (pV0mag < 1e-9)
      return; // protect against zero momentum

    const TVector3 u = pV0 * (1.0 / pV0mag);

    double pLpos = ppos.Dot(u);
    double pLneg = pneg.Dot(u);

    // qT: transverse momentum of the + track w.r.t. V0 direction
    TVector3 pTpos = ppos - pLpos * u;
    qT = pTpos.Mag();

    // α: longitudinal asymmetry (uses + and − labels by charge)
    double denom = pLpos + pLneg;
    if (std::abs(denom) < 1e-9)
      return; // avoid 0 division (unphysical for V0s)

    alpha = (pLpos - pLneg) / denom; // equivalently / pV0mag
  }

  // Daughters DCA selection
  template <typename T>
  bool passesDCASelectionDaughters(const T& v0)
  {

    bool isSelected{false};
    // const double ptPos{std::sqrt(std::pow(v0.pxpos(), 2.0) + std::pow(v0.pypos(), 2.0))};
    // const double ptNeg{std::sqrt(std::pow(v0.pxneg(), 2.0) + std::pow(v0.pyneg(), 2.0))};
    // const double dcaPtDepPos{0.0105 + 0.035 * std::pow(ptPos, -1.1)};
    // const double dcaPtDepNeg{0.0105 + 0.035 * std::pow(ptNeg, -1.1)};
    // const double dcaSelPos = std::max(0.1, dcaPtDepPos);
    // const double dcaSelNeg = std::max(0.1, dcaPtDepNeg);

    const double dcaPos{std::fabs(v0.dcapostopv())};
    const double dcaNeg{std::fabs(v0.dcanegtopv())};

    isSelected = dcaPos > v0Selections.dcapostopv && dcaNeg > v0Selections.dcanegtopv ? true : false;
    return isSelected;
  }

  template <typename T>
  bool passesTrackSelectionDaughters(const T& track)
  {

    bool isSelected = trkSelDaugthers.IsSelected(track) ? true : false;
    return isSelected;
  }

  // V0 topological selection
  template <typename T>
  bool passesV0TopologicalSelection(const T& v0)
  {

    bool isSelected = v0.v0radius() > v0Selections.v0radius && v0.v0radius() < v0Selections.v0radiusMax && passesDCASelectionDaughters(v0) && v0.v0cosPA() > v0Selections.v0cospa && v0.dcaV0daughters() < v0Selections.dcav0dau ? true : false;

    return isSelected;
  }

  template <typename C, typename T>
  bool passesK0Selection(const C& collision, const T& v0)
  {
    // Selection on rapiditty, proper lifetime, and Nsigma Pion

    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPi())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPi())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaPi())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaPi())};
    const double posRadiusNsigma{std::sqrt(std::pow(posTPCNsigma, 2.) + std::pow(posTOFNsigma, 2.))};
    const double negRadiusNsigma{std::sqrt(std::pow(negTPCNsigma, 2.) + std::pow(negTOFNsigma, 2.))};

    // bool isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S") && std::abs(v0.yK0Short()) < v0Selections.rapidityCut && posTPCNsigma < v0Selections.tpcPidNsigmaCut && negTPCNsigma < v0Selections.tpcPidNsigmaCut ? true : false;
    bool isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S") && std::abs(v0.yK0Short()) < v0Selections.rapidityCut && posTrack.hasTOF() && negTrack.hasTOF() && posRadiusNsigma < v0Selections.tpcPidNsigmaCut && negRadiusNsigma < v0Selections.tpcPidNsigmaCut ? true : false;

    return isSelected;
  }

  template <typename C, typename T>
  bool passesLambdaSelection(const C& collision, const T& v0)
  {
    // Selection on rapiditty, proper lifetime, and Nsigma Pion

    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPr())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPi())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaPr())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaPi())};
    const double posRadiusNsigma{std::sqrt(std::pow(posTPCNsigma, 2.) + std::pow(posTOFNsigma, 2.))};
    const double negRadiusNsigma{std::sqrt(std::pow(negTPCNsigma, 2.) + std::pow(negTOFNsigma, 2.))};

    bool isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda") && std::abs(v0.yLambda()) < v0Selections.rapidityCut && posTrack.hasTOF() && negTrack.hasTOF() && posRadiusNsigma < v0Selections.tpcPidNsigmaCut && negRadiusNsigma < v0Selections.tpcPidNsigmaCut ? true : false;

    return isSelected;
  }

  template <typename C, typename T>
  bool passesAntiLambdaSelection(const C& collision, const T& v0)
  {
    // Selection on rapiditty, proper lifetime, and Nsigma Pion

    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPi())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPr())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaPi())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaPr())};
    const double posRadiusNsigma{std::sqrt(std::pow(posTPCNsigma, 2.) + std::pow(posTOFNsigma, 2.))};
    const double negRadiusNsigma{std::sqrt(std::pow(negTPCNsigma, 2.) + std::pow(negTOFNsigma, 2.))};

    bool isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda") && std::abs(v0.yLambda()) < v0Selections.rapidityCut && posTrack.hasTOF() && negTrack.hasTOF() && posRadiusNsigma < v0Selections.tpcPidNsigmaCut && negRadiusNsigma < v0Selections.tpcPidNsigmaCut ? true : false;

    return isSelected;
  }

  template <typename C, typename T>
  bool passesGammaSelection(const C& /*collision*/, const T& v0)
  {
    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaEl())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaEl())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaEl())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaEl())};
    const double posRadiusNsigma{std::sqrt(std::pow(posTPCNsigma, 2.) + std::pow(posTOFNsigma, 2.))};
    const double negRadiusNsigma{std::sqrt(std::pow(negTPCNsigma, 2.) + std::pow(negTOFNsigma, 2.))};
    const float yGamma = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma);

    if (!(std::abs(yGamma) < v0Selections.rapidityCut))
      return false;

    if (!(posTrack.hasTOF() && negTrack.hasTOF()))
      return false;

    bool isSelected = posRadiusNsigma < v0Selections.tpcPidNsigmaCut && negRadiusNsigma < v0Selections.tpcPidNsigmaCut ? true : false;
    // bool isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda") && std::abs(v0.yLambda()) < v0Selections.rapidityCut && posTrack.hasTOF() && negTrack.hasTOF() && posRadiusNsigma < v0Selections.tpcPidNsigmaCut && negRadiusNsigma < v0Selections.tpcPidNsigmaCut ? true : false;

    return isSelected;
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  bool passesGeometricalCut()
  {
    return true;
  }

  void PhiPrime(float& phi, const int& magField, const int& charge)
  {

    if (magField < 0) // for negatve polarity field
      phi = o2::constants::math::TwoPI - phi;
    if (charge < 0) // for negatve charge
      phi = o2::constants::math::TwoPI - phi;

    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    registry.fill(HIST("EventCounter"), EvCutLabel::All);
    if (!col.sel8()) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::SelEigth);

    if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::NoSameBunchPileup);

    if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::IsGoodZvtxFT0vsPV);

    if (isNoCollInTimeRangeStrict) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInTimeRangeStrict);
    }

    if (isNoCollInTimeRangeStandard) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInTimeRangeStandard);
    }

    if (isNoCollInRofStrict) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInRofStrict);
    }

    if (isNoCollInRofStandard) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInRofStandard);
    }

    if (isNoHighMultCollInPrevRof) {
      if (!col.selection_bit(
            o2::aod::evsel::kNoHighMultCollInPrevRof)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoHighMultCollInPrevRof);
    }

    // To be used in combination with FT0C-based occupancy
    if (isNoCollInTimeRangeNarrow) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInTimeRangeNarrow);
    }

    if (isOccupancyCut) {
      auto occuValue{isApplyFT0CbasedOccupancy ? col.ft0cOccupancyInTimeRange() : col.trackOccupancyInTimeRange()};
      if (occuValue < minOccCut || occuValue > maxOccCut) {
        return false;
      }
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::OccuCut);

    if (col.centFT0C() < minT0CcentCut || col.centFT0C() > maxT0CcentCut) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::Centrality);

    // Z-vertex position cut
    if (std::fabs(col.posZ()) > posZcut) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::VtxZ);

    return true;
  }

  template <typename T, typename U>
  void getPTpowers(const T& pTs, const T& vecEff, const T& vecFD, U& pOne,
                   U& wOne, U& pTwo, U& wTwo, U& pThree, U& wThree,
                   U& pFour, U& wFour)
  {
    pOne = wOne = pTwo = wTwo = pThree = wThree = pFour = wFour = 0.;
    for (std::size_t i = 0; i < pTs.size(); ++i) {
      const double pTi{pTs.at(i)};
      const double eFFi{vecEff.at(i)};
      const double fDi{vecFD.at(i)};
      const double wEighti{std::pow(eFFi, -1.) * fDi};
      pOne += wEighti * pTi;
      wOne += wEighti;
      pTwo += std::pow(wEighti * pTi, 2.);
      wTwo += std::pow(wEighti, 2.);
      pThree += std::pow(wEighti * pTi, 3.);
      wThree += std::pow(wEighti, 3.);
      pFour += std::pow(wEighti * pTi, 4.);
      wFour += std::pow(wEighti, 4.);
    }
  }

  void loadCorrections(uint64_t timeStamp)
  {
    //        if (cfg.correctionsLoaded) return;

    if (paTHEff.value.empty() == false) {
      cfg.hEfficiency =
        ccdb->getForTimeStamp<TH2F>(paTHEff, timeStamp);
      if (cfg.hEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s",
             paTHEff.value.c_str());
      }
    }

    if (paTHFD.value.empty() == false) {
      cfg.hFeedDown = ccdb->getForTimeStamp<TH2F>(paTHFD, timeStamp);
      if (cfg.hFeedDown == nullptr) {
        LOGF(fatal, "Could not load feed down histogram from %s",
             paTHFD.value.c_str());
      }
    }
    cfg.correctionsLoaded = true;
  }

  void loadNchCalibrations(uint64_t timeStamp)
  {
    if (paTHmeanNch.value.empty() == false) {
      cfgNch.hMeanNch =
        ccdb->getForTimeStamp<TH1F>(paTHmeanNch, timeStamp);
      if (cfgNch.hMeanNch == nullptr) {
        LOGF(fatal, "Could not load hMeanNch histogram from %s",
             paTHmeanNch.value.c_str());
      }
    }

    if (paTHsigmaNch.value.empty() == false) {
      cfgNch.hSigmaNch =
        ccdb->getForTimeStamp<TH1F>(paTHsigmaNch, timeStamp);
      if (cfgNch.hSigmaNch == nullptr) {
        LOGF(fatal, "Could not load hSigmaNch histogram from %s",
             paTHsigmaNch.value.c_str());
      }
    }
    cfgNch.calibrationsLoaded = true;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<piKpRAA>(cfgc)};
}
