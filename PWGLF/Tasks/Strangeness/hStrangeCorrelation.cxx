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
///
/// \brief This task serves to do hadron-(strange hadron) correlation studies.
///  The yield will be calculated using the two-particle correlation method.
///  Trigger particle : Hadrons
///  Associated Particles : V0s or Cascades
///  this task requires the hStrangeCorrelationFilter to have been run before.
///
/// \author Kai Cui (kaicui@mails.ccnu.edu.cn)
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)
/// \author David Dobrigkeit Chinellato (david.dobrigkeit.chinellato@cern.ch)
/// \author Zhongbao Yin (Zhong-Bao.Yin@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFHStrangeCorrelationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Framework/StaticFor.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

// simple checkers
#define bitset(var, nbit) ((var) |= (1 << (nbit)))
#define bitcheck(var, nbit) ((var) & (1 << (nbit)))

using TracksComplete = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
using V0DatasWithoutTrackX = soa::Join<aod::V0Indices, aod::V0Cores>;

struct correlateStrangeness {
  // for efficiency corrections if requested
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Service<o2::framework::O2DatabasePDG> pdgDB;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> doCorrelationK0Short{"doCorrelationK0Short", true, "do K0Short correlation"};
  Configurable<bool> doCorrelationLambda{"doCorrelationLambda", false, "do Lambda correlation"};
  Configurable<bool> doCorrelationAntiLambda{"doCorrelationAntiLambda", false, "do AntiLambda correlation"};
  Configurable<bool> doCorrelationXiMinus{"doCorrelationXiMinus", false, "do XiMinus correlation"};
  Configurable<bool> doCorrelationXiPlus{"doCorrelationXiPlus", false, "do XiPlus correlation"};
  Configurable<bool> doCorrelationOmegaMinus{"doCorrelationOmegaMinus", false, "do OmegaMinus correlation"};
  Configurable<bool> doCorrelationOmegaPlus{"doCorrelationOmegaPlus", false, "do OmegaPlus correlation"};
  Configurable<bool> doCorrelationPion{"doCorrelationPion", false, "do Pion correlation"};
  Configurable<bool> doGenEventSelection{"doGenEventSelection", true, "use event selections when performing closure test for the gen events"};
  Configurable<bool> selectINELgtZERO{"selectINELgtZERO", true, "select INEL>0 events"};
  Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
  Configurable<bool> skipUnderOverflowInTHn{"skipUnderOverflowInTHn", false, "skip under/overflow in THns"};
  Configurable<int> mixingParameter{"mixingParameter", 10, "how many events are mixed"};
  Configurable<bool> doMCassociation{"doMCassociation", false, "fill everything only for MC associated"};
  Configurable<bool> doTriggPhysicalPrimary{"doTriggPhysicalPrimary", false, "require physical primary for trigger particles"};
  Configurable<bool> doAssocPhysicalPrimary{"doAssocPhysicalPrimary", false, "require physical primary for associated particles"};
  Configurable<bool> doLambdaPrimary{"doLambdaPrimary", false, "do primary selection for lambda"};
  Configurable<bool> doAutocorrelationRejection{"doAutocorrelationRejection", true, "reject pairs where trigger Id is the same as daughter particle Id"};

  // Axes - configurable for smaller sizes
  ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0, 2 * M_PI}, "#phi"};
  ConfigurableAxis axisEta{"axisEta", {80, -0.8, +0.8}, "#eta"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta #varphi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {50, -1.6, 1.6}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.0, 1.0, 2.0, 3.0, 100}, "pt associated axis for histograms"};
  ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisMultCount{"axisMultCount", {VARIABLE_WIDTH, 0, 200, 400, 600, 800, 1000, 1400, 1800, 2300, 2800, 3300, 4000, 5000, 6000}, "Mixing bins - multiplicity"};
  ConfigurableAxis axisMassNSigma{"axisMassNSigma", {40, -2, 2}, "Axis for mass Nsigma"};

  // for topo var QA
  struct : ConfigurableGroup {
    Configurable<float> maxPeakNSigma{"maxPeakNSigma", 5, "Peak region edge definition (in sigma)"};
    Configurable<float> minBgNSigma{"minBgNSigma", 5, "Bg region edge closest to peak (in sigma)"};
    Configurable<float> maxBgNSigma{"maxBgNSigma", 10, "Bg region edge furthest to peak (in sigma)"};
  } massWindowConfigurations; // allows for gap between peak and bg in case someone wants to

  // Implementation of on-the-spot efficiency correction
  Configurable<bool> applyEfficiencyCorrection{"applyEfficiencyCorrection", false, "apply efficiency correction"};
  Configurable<bool> applyEfficiencyForTrigger{"applyEfficiencyForTrigger", false, "apply efficiency correction for the trigger particle"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> efficiencyCCDBPath{"efficiencyCCDBPath", "GLO/Config/GeometryAligned", "Path of the efficiency corrections"};

  // Configurables for doing subwagon systematics
  // Group all settings necessary for systematics in a specific ConfigurableGroup
  struct : ConfigurableGroup {
    std::string prefix = "systematics";
    // --- Track quality variations (single track, both trigger and assoc daughters)
    Configurable<int> minTPCNCrossedRowsTrigger{"minTPCNCrossedRowsTrigger", 70, "Minimum TPC crossed rows (trigger)"};
    Configurable<int> minTPCNCrossedRowsAssociated{"minTPCNCrossedRowsAssociated", 70, "Minimum TPC crossed rows (associated)"};
    Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};
    Configurable<int> triggerMaxTPCSharedClusters{"triggerMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive)"};
    Configurable<bool> triggerRequireL0{"triggerRequireL0", false, "require ITS L0 cluster for trigger"};

    // --- Trigger: DCA variation from basic formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstant{"dcaXYconstant", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};

    // --- Associated: topological variable variation (OK to vary all-at-once, at least for first study)
    Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> dcaV0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
    Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
    Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
    Configurable<float> v0RadiusMin{"v0radiusmin", 0.5, "v0radius"};
    Configurable<float> v0RadiusMax{"v0radiusmax", 200, "v0radius"};

    // cascade selections
    Configurable<double> casc_cospa{"casc_cospa", 0.95, "casc_cospa"};
    Configurable<float> casc_dcacascdau{"casc_dcacascdau", 1.0, "casc_dcacascdau"};
    Configurable<float> casc_dcabachtopv{"casc_dcabachtopv", 0.1, "casc_dcabachtopv"};
    Configurable<float> casc_cascradius{"casc_cascradius", 0.5, "casc_cascradius"};
    Configurable<float> casc_v0masswindow{"casc_v0masswindow", 0.01, "casc_v0masswindow"};
    Configurable<float> casc_mindcav0topv{"casc_mindcav0topv", 0.01, "casc_mindcav0topv"};

    // dE/dx for associated daughters
    Configurable<int> dEdxCompatibility{"dEdxCompatibility", 1, "0: loose, 1: normal, 2: tight. Defined in hStrangeCorrelationFilter"};

    // (N.B.: sources that can be investigated in post are not listed!)
  } systCuts;

  // objects to use for efficiency corrections
  TH2F* hEfficiencyTrigger;
  TH2F* hEfficiencyPion;
  TH2F* hEfficiencyK0Short;
  TH2F* hEfficiencyLambda;
  TH2F* hEfficiencyAntiLambda;
  TH2F* hEfficiencyXiMinus;
  TH2F* hEfficiencyXiPlus;
  TH2F* hEfficiencyOmegaMinus;
  TH2F* hEfficiencyOmegaPlus;

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningType colBinning{{axisVtxZ, axisMult}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.

  // collision slicing for mixed events
  Preslice<aod::TriggerTracks> collisionSliceTracks = aod::triggerTracks::collisionId;
  Preslice<aod::AssocV0s> collisionSliceV0s = aod::assocV0s::collisionId;
  Preslice<aod::AssocCascades> collisionSliceCascades = aod::assocCascades::collisionId;
  Preslice<aod::AssocPions> collisionSlicePions = aod::assocPions::collisionId;
  Preslice<aod::McParticles> perCollision = aod::mcparticle::mcCollisionId;

  static constexpr std::string_view v0names[] = {"K0Short", "Lambda", "AntiLambda"};
  static constexpr std::string_view cascadenames[] = {"XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
  static constexpr std::string_view particlenames[] = {"Pion", "K0Short", "Lambda", "AntiLambda", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};

  uint8_t doCorrelation;
  int mRunNumber;

  std::vector<std::vector<float>> axisRanges;

  /// Function to aid in calculating delta-phi
  /// \param phi1 first phi value
  /// \param phi2 second phi value
  Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2)
  {
    Double_t deltaPhi = phi1 - phi2;
    if (deltaPhi < -TMath::Pi() / 2.) {
      deltaPhi += 2. * TMath::Pi();
    }
    if (deltaPhi > 3 * TMath::Pi() / 2.) {
      deltaPhi -= 2. * TMath::Pi();
    }
    return deltaPhi;
  }

  /// Function to load efficiencies to memory from CCDB
  /// \param bc provided such that the run number can be used
  void initEfficiencyFromCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOG(info) << "Loading efficiencies from CCDB for run " << mRunNumber << " now...";
    auto timeStamp = bc.timestamp();

    TList* listEfficiencies = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, timeStamp);

    if (!listEfficiencies) {
      LOG(fatal) << "Problem getting TList object with efficiencies!";
    }

    hEfficiencyTrigger = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyTrigger"));
    hEfficiencyK0Short = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyK0Short"));
    hEfficiencyLambda = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyLambda"));
    hEfficiencyAntiLambda = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyAntiLambda"));
    hEfficiencyXiMinus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyXiMinus"));
    hEfficiencyXiPlus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyXiPlus"));
    hEfficiencyOmegaMinus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyOmegaMinus"));
    hEfficiencyOmegaPlus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyOmegaPlus"));
    LOG(info) << "Efficiencies now loaded for " << mRunNumber;
  }
  template <class TTrack>
  bool isValidTrigger(TTrack track)
  {
    if (track.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsTrigger) {
      return false; // crossed rows
    }
    if (!track.hasITS() && systCuts.triggerRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > systCuts.triggerMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(bitcheck(track.itsClusterMap(), 0)) && systCuts.triggerRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    // systematic variations: trigger DCAxy
    if (std::abs(track.dcaXY()) > systCuts.dcaXYconstant + systCuts.dcaXYpTdep * std::abs(track.signed1Pt())) {
      return false;
    }
    if (track.pt() > axisRanges[3][1] || track.pt() < axisRanges[3][0]) {
      return false;
    }
    return true;
  }
  void fillCorrelationsV0(aod::TriggerTracks const& triggers, aod::AssocV0s const& assocs, bool mixing, float pvz, float mult)
  {
    for (auto& triggerTrack : triggers) {
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg))
        continue;

      if (!mixing) {
        float efficiency = 1.0f;
        if (applyEfficiencyForTrigger) {
          efficiency = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        }
        float weight = (applyEfficiencyForTrigger) ? 1. / efficiency : 1.0f;
        histos.fill(HIST("sameEvent/TriggerParticlesV0"), trigg.pt(), mult, weight);
      }

      for (auto& assocCandidate : assocs) {
        auto assoc = assocCandidate.v0Core_as<V0DatasWithoutTrackX>();

        //---] syst cuts [---
        if (assoc.v0radius() < systCuts.v0RadiusMin || assoc.v0radius() > systCuts.v0RadiusMax ||
            std::abs(assoc.dcapostopv()) < systCuts.dcapostopv || std::abs(assoc.dcanegtopv()) < systCuts.dcanegtopv ||
            assoc.v0cosPA() < systCuts.v0cospa || assoc.dcaV0daughters() > systCuts.dcaV0dau)
          continue;

        //---] removing autocorrelations [---
        auto postrack = assoc.posTrack_as<TracksComplete>();
        auto negtrack = assoc.negTrack_as<TracksComplete>();
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == postrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsV0"), 0.5);
            continue;
          }
          if (trigg.globalIndex() == negtrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsV0"), 0.5);
            continue;
          }
        }

        //---] track quality check [---
        if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
          continue;

        float deltaphi = ComputeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
          continue;
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
          continue;
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
          continue;

        TH2F* hEfficiencyV0[3];
        hEfficiencyV0[0] = hEfficiencyK0Short;
        hEfficiencyV0[1] = hEfficiencyLambda;
        hEfficiencyV0[2] = hEfficiencyAntiLambda;
        static_for<0, 2>([&](auto i) {
          constexpr int index = i.value;
          float efficiency = 1.0f;
          if (applyEfficiencyCorrection) {
            efficiency = hEfficiencyV0[index]->Interpolate(ptassoc, assoc.eta());
          }
          if (applyEfficiencyForTrigger) {
            efficiency = efficiency * hEfficiencyTrigger->Interpolate(pttrigger, trigg.eta());
          }

          float weight = (applyEfficiencyCorrection || applyEfficiencyForTrigger) ? 1. / efficiency : 1.0f;
          if (bitcheck(doCorrelation, index) && (!applyEfficiencyCorrection || efficiency != 0)) {
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < -massWindowConfigurations.minBgNSigma)
              histos.fill(HIST("sameEvent/LeftBg/") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < +massWindowConfigurations.maxPeakNSigma)
              histos.fill(HIST("sameEvent/Signal/") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < +massWindowConfigurations.maxBgNSigma)
              histos.fill(HIST("sameEvent/RightBg/") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < -massWindowConfigurations.minBgNSigma)
              histos.fill(HIST("mixedEvent/LeftBg/") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < +massWindowConfigurations.maxPeakNSigma)
              histos.fill(HIST("mixedEvent/Signal/") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < +massWindowConfigurations.maxBgNSigma)
              histos.fill(HIST("mixedEvent/RightBg/") + HIST(v0names[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
          }
        });
      }
    }
  }

  void fillCorrelationsCascade(aod::TriggerTracks const& triggers, aod::AssocCascades const& assocs, bool mixing, float pvx, float pvy, float pvz, float mult)
  {
    for (auto& triggerTrack : triggers) {
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg))
        continue;

      if (!mixing) {
        float efficiency = 1.0f;
        if (applyEfficiencyForTrigger) {
          efficiency = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        }
        float weight = (applyEfficiencyForTrigger) ? 1. / efficiency : 1.0f;
        histos.fill(HIST("sameEvent/TriggerParticlesCascade"), trigg.pt(), mult, weight);
      }
      for (auto& assocCandidate : assocs) {
        auto assoc = assocCandidate.cascData();

        //---] syst cuts [---
        if (std::abs(assoc.dcapostopv()) < systCuts.dcapostopv ||
            std::abs(assoc.dcanegtopv()) < systCuts.dcanegtopv ||
            assoc.dcabachtopv() < systCuts.casc_dcabachtopv ||
            assoc.dcaV0daughters() > systCuts.dcaV0dau ||
            assoc.dcacascdaughters() > systCuts.casc_dcacascdau ||
            assoc.v0cosPA(pvx, pvy, pvz) < systCuts.v0cospa ||
            assoc.casccosPA(pvx, pvy, pvz) < systCuts.casc_cospa ||
            assoc.cascradius() < systCuts.casc_cascradius ||
            std::abs(assoc.dcav0topv(pvx, pvy, pvz)) < systCuts.casc_mindcav0topv ||
            std::abs(assoc.mLambda() - pdgDB->Mass(3122)) > systCuts.casc_v0masswindow)
          continue;

        //---] removing autocorrelations [---
        auto postrack = assoc.posTrack_as<TracksComplete>();
        auto negtrack = assoc.negTrack_as<TracksComplete>();
        auto bachtrack = assoc.bachelor_as<TracksComplete>();
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == postrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsCascades"), 0.5);
            continue;
          }
          if (trigg.globalIndex() == negtrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsCascades"), 0.5);
            continue;
          }
          if (trigg.globalIndex() == bachtrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsCascades"), 0.5);
            continue;
          }
        }

        //---] track quality check [---
        if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || bachtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
          continue;

        float deltaphi = ComputeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
          continue;
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
          continue;
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
          continue;

        TH2F* hEfficiencyCascade[4];
        hEfficiencyCascade[0] = hEfficiencyXiMinus;
        hEfficiencyCascade[1] = hEfficiencyXiPlus;
        hEfficiencyCascade[2] = hEfficiencyOmegaMinus;
        hEfficiencyCascade[3] = hEfficiencyOmegaPlus;

        static_for<0, 3>([&](auto i) {
          constexpr int index = i.value;
          float efficiency = 1.0f;
          if (applyEfficiencyCorrection) {
            efficiency = hEfficiencyCascade[index]->GetBinContent(hEfficiencyCascade[index]->GetXaxis()->FindBin(ptassoc), hEfficiencyCascade[index]->GetYaxis()->FindBin(assoc.eta()));
          }
          if (applyEfficiencyForTrigger) {
            efficiency = efficiency * hEfficiencyTrigger->Interpolate(pttrigger, trigg.eta());
          }
          float weight = (applyEfficiencyCorrection || applyEfficiencyForTrigger) ? 1. / efficiency : 1.0f;
          if (bitcheck(doCorrelation, index + 3) && (!applyEfficiencyCorrection || efficiency != 0)) {
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < -massWindowConfigurations.minBgNSigma)
              histos.fill(HIST("sameEvent/LeftBg/") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < +massWindowConfigurations.maxPeakNSigma)
              histos.fill(HIST("sameEvent/Signal/") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < +massWindowConfigurations.maxBgNSigma)
              histos.fill(HIST("sameEvent/RightBg/") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < -massWindowConfigurations.minBgNSigma)
              histos.fill(HIST("mixedEvent/LeftBg/") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < +massWindowConfigurations.maxPeakNSigma)
              histos.fill(HIST("mixedEvent/Signal/") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
            if (assocCandidate.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(index) && assocCandidate.invMassNSigma(index) < +massWindowConfigurations.maxBgNSigma)
              histos.fill(HIST("mixedEvent/RightBg/") + HIST(cascadenames[index]), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult, weight);
          }
        });
      }
    }
  }

  void fillCorrelationsPion(aod::TriggerTracks const& triggers, aod::AssocPions const& assocs, bool mixing, float pvz, float mult)
  {

    for (auto& triggerTrack : triggers) {
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg))
        continue;

      if (!mixing)
        histos.fill(HIST("sameEvent/TriggerParticlesPion"), trigg.pt(), mult);
      for (auto& assocTrack : assocs) {
        auto assoc = assocTrack.track_as<TracksComplete>();

        //---] removing autocorrelations [---
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == assoc.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsPions"), 0.5);
            continue;
          }
        }

        //---] track quality check [---
        if (assoc.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
          continue;

        float deltaphi = ComputeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
          continue;
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
          continue;
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
          continue;

        if (!mixing)
          histos.fill(HIST("sameEvent/Pion"), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);
        else
          histos.fill(HIST("mixedEvent/Pion"), deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);
      }
    }
  }

  void init(InitContext const&)
  {
    mRunNumber = 0;
    hEfficiencyPion = 0x0;
    hEfficiencyK0Short = 0x0;
    hEfficiencyLambda = 0x0;
    hEfficiencyAntiLambda = 0x0;
    hEfficiencyXiMinus = 0x0;
    hEfficiencyXiPlus = 0x0;
    hEfficiencyOmegaMinus = 0x0;
    hEfficiencyOmegaPlus = 0x0;

    // set bitmap for convenience
    doCorrelation = 0;
    if (doCorrelationK0Short)
      bitset(doCorrelation, 0);
    if (doCorrelationLambda)
      bitset(doCorrelation, 1);
    if (doCorrelationAntiLambda)
      bitset(doCorrelation, 2);
    if (doCorrelationXiMinus)
      bitset(doCorrelation, 3);
    if (doCorrelationXiPlus)
      bitset(doCorrelation, 4);
    if (doCorrelationOmegaMinus)
      bitset(doCorrelation, 5);
    if (doCorrelationOmegaPlus)
      bitset(doCorrelation, 6);
    if (doCorrelationPion)
      bitset(doCorrelation, 7);

    // Store axis ranges to prevent spurious filling
    // axis status:
    // --- Delta-phi is safe -> math forbids insanity
    // --- Delta-eta depends on pre-filter -> check
    // --- pT assoc depends on binning -> check
    // --- vertex Z is safe -> skipped at evsel level
    // --- multiplicity -> check

    // grab axis edge from ConfigurableAxes
    const AxisSpec preAxisDeltaPhi{axisDeltaPhi, "#Delta#varphi"};
    const AxisSpec preAxisDeltaEta{axisDeltaEta, "#Delta#eta"};
    const AxisSpec preAxisPtAssoc{axisPtAssoc, "#it{p}_{T}^{assoc} (GeV/c)"};
    const AxisSpec preAxisPtTrigger{axisPtTrigger, "#it{p}_{T}^{trigger} (GeV/c)"};
    const AxisSpec preAxisVtxZ{axisVtxZ, "vertex Z (cm)"};
    const AxisSpec preAxisMult{axisMult, "mult percentile"};

    // store the original axes in specific TH1Cs for completeness
    histos.add("axes/hDeltaPhiAxis", "", kTH1C, {preAxisDeltaPhi});
    histos.add("axes/hDeltaEtaAxis", "", kTH1C, {preAxisDeltaEta});
    histos.add("axes/hPtAssocAxis", "", kTH1C, {preAxisPtAssoc});
    histos.add("axes/hPtTriggerAxis", "", kTH1C, {preAxisPtTrigger});
    histos.add("axes/hVertexZAxis", "", kTH1C, {preAxisVtxZ});
    histos.add("axes/hMultAxis", "", kTH1C, {preAxisMult});

    std::vector<double> edgesDeltaPhiOrig = preAxisDeltaPhi.binEdges;
    std::vector<double> edgesDeltaEtaOrig = preAxisDeltaEta.binEdges;
    std::vector<double> edgesPtAssocOrig = preAxisPtAssoc.binEdges;
    std::vector<double> edgesPtTriggerOrig = preAxisPtTrigger.binEdges;
    std::vector<double> edgesVtxZOrig = preAxisVtxZ.binEdges;
    std::vector<double> edgesMultOrig = preAxisMult.binEdges;

    std::vector<float> rangesDeltaPhi = {static_cast<float>(edgesDeltaPhiOrig[0]), static_cast<float>(edgesDeltaPhiOrig[edgesDeltaPhiOrig.size() - 1])};
    std::vector<float> rangesDeltaEta = {static_cast<float>(edgesDeltaEtaOrig[0]), static_cast<float>(edgesDeltaEtaOrig[edgesDeltaEtaOrig.size() - 1])};
    std::vector<float> rangesPtAssoc = {static_cast<float>(edgesPtAssocOrig[0]), static_cast<float>(edgesPtAssocOrig[edgesPtAssocOrig.size() - 1])};
    std::vector<float> rangesPtTrigger = {static_cast<float>(edgesPtTriggerOrig[0]), static_cast<float>(edgesPtTriggerOrig[edgesPtTriggerOrig.size() - 1])};
    std::vector<float> rangesVtxZ = {static_cast<float>(edgesVtxZOrig[0]), static_cast<float>(edgesVtxZOrig[edgesVtxZOrig.size() - 1])};
    std::vector<float> rangesMult = {static_cast<float>(edgesMultOrig[0]), static_cast<float>(edgesMultOrig[edgesMultOrig.size() - 1])};

    axisRanges.emplace_back(rangesDeltaPhi);
    axisRanges.emplace_back(rangesDeltaEta);
    axisRanges.emplace_back(rangesPtAssoc);
    axisRanges.emplace_back(rangesPtTrigger);
    axisRanges.emplace_back(rangesVtxZ);
    axisRanges.emplace_back(rangesMult);

    std::vector<double> edgesDeltaPhi;
    std::vector<double> edgesDeltaEta;
    std::vector<double> edgesPtAssoc;
    std::vector<double> edgesPtTrigger;
    std::vector<double> edgesVtxZ;
    std::vector<double> edgesMult;

    // v--- skipUnderOverflowInTHn ---v
    //
    // if enabled, this will change the axes such that they will solely cover the interval from
    // edge[1] to edge[n-1]; this will mean that the bin 1 and bin N will be stored in
    // under / overflow bins and will have to be manually unpacked. Do not forget to do the manual
    // unpacking a posteriori!
    //
    // this feature is meant to save memory conveniently.
    // it should actually be implemented centrally in ROOT but ok, this will do it for now.

    int offset = skipUnderOverflowInTHn ? 1 : 0;
    // ===] delta-phi [===
    if (!preAxisDeltaPhi.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesDeltaPhiOrig.size()) - offset; i++)
        edgesDeltaPhi.emplace_back(edgesDeltaPhiOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesDeltaPhiOrig[0];
      double delta = (edgesDeltaPhiOrig[1] - edgesDeltaPhiOrig[0]) / preAxisDeltaPhi.nBins.value();
      for (int i = offset; i < preAxisDeltaPhi.nBins.value() + 1 - offset; i++)
        edgesDeltaPhi.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] delta-eta [===
    if (!preAxisDeltaEta.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesDeltaEtaOrig.size()) - offset; i++)
        edgesDeltaEta.emplace_back(edgesDeltaEtaOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesDeltaEtaOrig[0];
      double delta = (edgesDeltaEtaOrig[1] - edgesDeltaEtaOrig[0]) / preAxisDeltaEta.nBins.value();
      for (int i = offset; i < preAxisDeltaEta.nBins.value() + 1 - offset; i++)
        edgesDeltaEta.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] pt assoc [===
    if (!preAxisPtAssoc.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesPtAssocOrig.size()) - offset; i++)
        edgesPtAssoc.emplace_back(edgesPtAssocOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesPtAssocOrig[0];
      double delta = (edgesPtAssocOrig[1] - edgesPtAssocOrig[0]) / preAxisPtAssoc.nBins.value();
      for (int i = offset; i < preAxisPtAssoc.nBins.value() + 1 - offset; i++)
        edgesPtAssoc.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] pt trigger [===
    if (!preAxisPtTrigger.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesPtTriggerOrig.size()) - offset; i++)
        edgesPtTrigger.emplace_back(edgesPtTriggerOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesPtTriggerOrig[0];
      double delta = (edgesPtTriggerOrig[1] - edgesPtTriggerOrig[0]) / preAxisPtTrigger.nBins.value();
      for (int i = offset; i < preAxisPtTrigger.nBins.value() + 1 - offset; i++)
        edgesPtTrigger.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] vtx Z [===
    if (!preAxisVtxZ.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesVtxZOrig.size()) - offset; i++)
        edgesVtxZ.emplace_back(edgesVtxZOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesVtxZOrig[0];
      double delta = (edgesVtxZOrig[1] - edgesVtxZOrig[0]) / preAxisVtxZ.nBins.value();
      for (int i = offset; i < preAxisVtxZ.nBins.value() + 1 - offset; i++)
        edgesVtxZ.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] mult percentile [===
    if (!preAxisMult.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesMultOrig.size()) - offset; i++)
        edgesMult.emplace_back(edgesMultOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesMultOrig[0];
      double delta = (edgesMultOrig[1] - edgesMultOrig[0]) / preAxisMult.nBins.value();
      for (int i = offset; i < preAxisMult.nBins.value() + 1 - offset; i++)
        edgesMult.emplace_back(min + static_cast<double>(i) * delta);
    }

    LOGF(info, "Initialized THnF axis delta-phi with %i bins.", edgesDeltaPhi.size() - 1);
    LOGF(info, "Initialized THnF axis delta-eta with %i bins.", edgesDeltaEta.size() - 1);
    LOGF(info, "Initialized THnF axis pTassoc with %i bins.", edgesPtAssoc.size() - 1);
    LOGF(info, "Initialized THnF axis pTtrigger with %i bins.", edgesPtTrigger.size() - 1);
    LOGF(info, "Initialized THnF axis vertex-Z with %i bins.", edgesVtxZ.size() - 1);
    LOGF(info, "Initialized THnF axis multiplicity with %i bins.", edgesMult.size() - 1);

    const AxisSpec axisDeltaPhiNDim{edgesDeltaPhi, "#Delta#varphi"};
    const AxisSpec axisDeltaEtaNDim{edgesDeltaEta, "#Delta#eta"};
    const AxisSpec axisPtAssocNDim{edgesPtAssoc, "#it{p}_{T}^{assoc} (GeV/c)"};
    const AxisSpec axisPtTriggerNDim{edgesPtTrigger, "#it{p}_{T}^{trigger} (GeV/c)"};
    const AxisSpec axisVtxZNDim{edgesVtxZ, "vertex Z (cm)"};
    const AxisSpec axisMultNDim{edgesMult, "mult percentile"};

    if (bitcheck(doCorrelation, 0)) {
      histos.add("h3dK0ShortSpectrum", "h3dK0ShortSpectrum", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("h3dK0ShortSpectrumY", "h3dK0ShortSpectrumY", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("sameEvent/Signal/K0Short", "K0Short", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
    }
    if (bitcheck(doCorrelation, 1)) {
      histos.add("h3dLambdaSpectrum", "h3dLambdaSpectrum", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("h3dLambdaSpectrumY", "h3dLambdaSpectrumY", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("sameEvent/Signal/Lambda", "Lambda", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
    }
    if (bitcheck(doCorrelation, 2)) {
      histos.add("h3dAntiLambdaSpectrum", "h3dAntiLambdaSpectrum", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("h3dAntiLambdaSpectrumY", "h3dAntiLambdaSpectrumY", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("sameEvent/Signal/AntiLambda", "AntiLambda", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
    }
    if (bitcheck(doCorrelation, 3)) {
      histos.add("h3dXiMinusSpectrum", "h3dXiMinusSpectrum", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("h3dXiMinusSpectrumY", "h3dXiMinusSpectrumY", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("sameEvent/Signal/XiMinus", "XiMinus", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
    }
    if (bitcheck(doCorrelation, 4)) {
      histos.add("h3dXiPlusSpectrum", "h3dXiPlusSpectrum", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("h3dXiPlusSpectrumY", "h3dXiPlusSpectrumY", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("sameEvent/Signal/XiPlus", "XiPlus", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
    }
    if (bitcheck(doCorrelation, 5)) {
      histos.add("h3dOmegaMinusSpectrum", "h3dOmegaMinusSpectrum", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("h3dOmegaMinusSpectrumY", "h3dOmegaMinusSpectrumY", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("sameEvent/Signal/OmegaMinus", "OmegaMinus", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
    }
    if (bitcheck(doCorrelation, 6)) {
      histos.add("h3dOmegaPlusSpectrum", "h3dOmegaPlusSpectrum", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("h3dOmegaPlusSpectrumY", "h3dOmegaPlusSpectrumY", kTH3F, {axisPtQA, axisMult, axisMassNSigma});
      histos.add("sameEvent/Signal/OmegaPlus", "OmegaPlus", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
    }
    if (bitcheck(doCorrelation, 7)) {
      histos.add("sameEvent/Pion", "Pion", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
    }
    LOGF(info, "Init THnFs done");
    if (doCorrelationK0Short || doCorrelationLambda || doCorrelationAntiLambda || doCorrelationXiMinus || doCorrelationXiPlus || doCorrelationOmegaMinus || doCorrelationOmegaPlus) {
      histos.addClone("sameEvent/Signal/", "sameEvent/LeftBg/");
      histos.addClone("sameEvent/Signal/", "sameEvent/RightBg/");
    }

    // mixed-event correlation functions
    if (doprocessMixedEventHV0s || doprocessMixedEventHCascades || doprocessMixedEventHPions) {
      histos.addClone("sameEvent/", "mixedEvent/");
    }

    // Some QA plots
    histos.add("hGeneratedQAPtTrigger", "hGeneratedQAPtTrigger", kTH2F, {axisPtQA, {5, -0.5f, 4.5f}});
    histos.add("hGeneratedQAPtAssociatedK0", "hGeneratedQAPtAssociatedK0", kTH2F, {axisPtQA, {5, -0.5f, 4.5f}});
    histos.add("hClosureQAPtTrigger", "hClosureQAPtTrigger", kTH2F, {axisPtQA, {5, -0.5f, 4.5f}});
    histos.add("hClosureQAPtAssociatedK0", "hClosureQAPtAssociatedK0", kTH2F, {axisPtQA, {5, -0.5f, 4.5f}});

    histos.add("hTrackEtaVsPtVsPhi", "hTrackEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hK0ShortEtaVsPtVsPhi", "hK0ShortEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hK0ShortEtaVsPtVsPhiBg", "hK0ShortEtaVsPtVsPhiBg", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hLambdaEtaVsPtVsPhi", "hLambdaEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hLambdaEtaVsPtVsPhiBg", "hLambdaEtaVsPtVsPhiBg", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hAntiLambdaEtaVsPtVsPhi", "hAntiLambdaEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hAntiLambdaEtaVsPtVsPhiBg", "hAntiLambdaEtaVsPtVsPhiBg", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hXiMinusEtaVsPtVsPhi", "hXiMinusEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hXiMinusEtaVsPtVsPhiBg", "hXiMinusEtaVsPtVsPhiBg", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hXiPlusEtaVsPtVsPhi", "hXiPlusEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hXiPlusEtaVsPtVsPhiBg", "hXiPlusEtaVsPtVsPhiBg", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hOmegaMinusEtaVsPtVsPhi", "hOmegaMinusEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hOmegaMinusEtaVsPtVsPhiBg", "hOmegaMinusEtaVsPtVsPhiBg", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hOmegaPlusEtaVsPtVsPhi", "hOmegaPlusEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hOmegaPlusEtaVsPtVsPhiBg", "hOmegaPlusEtaVsPtVsPhiBg", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hPionEtaVsPtVsPhi", "hPionEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hTriggerPrimaryEtaVsPt", "hTriggerPrimaryEtaVsPt", kTH3F, {axisPtQA, axisEta, axisMult});
    histos.add("hTriggerAllSelectedEtaVsPt", "hTriggerAllSelectedEtaVsPt", kTH3F, {axisPtQA, axisEta, axisMult});

    histos.add("hClosureTestEventCounter", "hClosureTestEventCounter", kTH1F, {{10, 0, 10}});

    histos.add("hNumberOfRejectedPairsV0", "hNumberOfRejectedPairsV0", kTH1F, {{1, 0, 1}});
    histos.add("hNumberOfRejectedPairsCascades", "hNumberOfRejectedPairsCascades", kTH1F, {{1, 0, 1}});
    histos.add("hNumberOfRejectedPairsPions", "hNumberOfRejectedPairsPions", kTH1F, {{1, 0, 1}});

    histos.add("sameEvent/TriggerParticlesV0", "TriggersV0", kTH2F, {axisPtQA, axisMult});
    histos.add("sameEvent/TriggerParticlesCascade", "TriggersCascade", kTH2F, {axisPtQA, axisMult});
    histos.add("sameEvent/TriggerParticlesPion", "TriggersPion", kTH2F, {axisPtQA, axisMult});

    // mixing QA
    histos.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{140, -0.5, 139.5}});
    histos.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{140, -0.5, 139.5}});
    histos.add("MixingQA/hMEpvz1", ";pvz;Entries", kTH1F, {{30, -15, 15}});
    histos.add("MixingQA/hMEpvz2", ";pvz;Entries", kTH1F, {{30, -15, 15}});

    // Event QA
    histos.add("EventQA/hMixingQA", "mixing QA", kTH1F, {{2, -0.5, 1.5}});
    histos.add("EventQA/hMult", "Multiplicity", kTH1F, {axisMult});
    histos.add("EventQA/hPvz", ";pvz;Entries", kTH1F, {{30, -15, 15}});
    histos.add("EventQA/hMultFT0vsTPC", ";centFT0M;multNTracksPVeta1", kTH2F, {{100, 0, 100}, {300, 0, 300}});

    // MC generated plots
    if (doprocessMCGenerated) {
      histos.add("Generated/hTrigger", "", kTH2F, {axisPtQA, axisEta});
      histos.add("Generated/hPion", "", kTH2F, {axisPtQA, axisEta});
      histos.add("Generated/hK0Short", "", kTH2F, {axisPtQA, axisEta});
      histos.add("Generated/hLambda", "", kTH2F, {axisPtQA, axisEta});
      histos.add("Generated/hAntiLambda", "", kTH2F, {axisPtQA, axisEta});
      histos.add("Generated/hXiMinus", "", kTH2F, {axisPtQA, axisEta});
      histos.add("Generated/hXiPlus", "", kTH2F, {axisPtQA, axisEta});
      histos.add("Generated/hOmegaMinus", "", kTH2F, {axisPtQA, axisEta});
      histos.add("Generated/hOmegaPlus", "", kTH2F, {axisPtQA, axisEta});

      histos.addClone("Generated/", "GeneratedWithPV/");

      // histograms within |y|<0.5, vs multiplicity
      histos.add("GeneratedWithPV/hPion_MidYVsMult", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hK0Short_MidYVsMult", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hLambda_MidYVsMult", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hAntiLambda_MidYVsMult", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hXiMinus_MidYVsMult", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hXiPlus_MidYVsMult", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hOmegaMinus_MidYVsMult", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hOmegaPlus_MidYVsMult", "", kTH2F, {axisPtQA, axisMult});

      histos.add("GeneratedWithPV/hPion_MidYVsMult_TwoPVsOrMore", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hK0Short_MidYVsMult_TwoPVsOrMore", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hLambda_MidYVsMult_TwoPVsOrMore", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hAntiLambda_MidYVsMult_TwoPVsOrMore", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hXiMinus_MidYVsMult_TwoPVsOrMore", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hXiPlus_MidYVsMult_TwoPVsOrMore", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hOmegaMinus_MidYVsMult_TwoPVsOrMore", "", kTH2F, {axisPtQA, axisMult});
      histos.add("GeneratedWithPV/hOmegaPlus_MidYVsMult_TwoPVsOrMore", "", kTH2F, {axisPtQA, axisMult});
    }
    if (doprocessClosureTest) {
      if (doCorrelationPion) {
        histos.add("ClosureTest/sameEvent/Pion", "Pion", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/hPion", "", kTH3F, {axisPtQA, axisEta, axisPhi});
      }
      if (doCorrelationK0Short) {
        histos.add("ClosureTest/sameEvent/K0Short", "K0Short", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/hK0Short", "", kTH3F, {axisPtQA, axisEta, axisPhi});
      }
      if (doCorrelationLambda) {
        histos.add("ClosureTest/sameEvent/Lambda", "Lambda", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/hLambda", "", kTH3F, {axisPtQA, axisEta, axisPhi});
      }
      if (doCorrelationAntiLambda) {
        histos.add("ClosureTest/sameEvent/AntiLambda", "AntiLambda", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/hAntiLambda", "", kTH3F, {axisPtQA, axisEta, axisPhi});
      }
      if (doCorrelationXiMinus) {
        histos.add("ClosureTest/sameEvent/XiMinus", "XiMinus", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/hXiMinus", "", kTH3F, {axisPtQA, axisEta, axisPhi});
      }
      if (doCorrelationXiPlus) {
        histos.add("ClosureTest/sameEvent/XiPlus", "XiPlus", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/XiPlus", "", kTH3F, {axisPtQA, axisEta, axisPhi});
      }
      if (doCorrelationOmegaMinus) {
        histos.add("ClosureTest/sameEvent/OmegaMinus", "OmegaMinus", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/hOmegaMinus", "", kTH3F, {axisPtQA, axisEta, axisPhi});
      }
      if (doCorrelationOmegaPlus) {
        histos.add("ClosureTest/sameEvent/OmegaPlus", "OmegaPlus", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/hOmegaPlus", "", kTH3F, {axisPtQA, axisEta, axisPhi});
      }
      histos.add("ClosureTest/hTrigger", "Trigger Tracks", kTH3F, {axisPtQA, axisEta, axisMult});
    }

    // visual inspection of sizes
    histos.print();

    // initialize CCDB *only* if efficiency correction requested
    // skip if not requested, saves a bit of time
    if (applyEfficiencyCorrection) {
      ccdb->setURL(ccdburl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
    }
  }

  void processSameEventHV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision,
                            aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks,
                            V0DatasWithoutTrackX const&, aod::V0sLinked const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    // ________________________________________________
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    if (TMath::Abs(collision.posZ()) > zVertexCut) {
      return;
    }
    if (collision.centFT0M() > axisRanges[5][1] || collision.centFT0M() < axisRanges[5][0]) {
      return;
    }
    if (!collision.isInelGt0() && selectINELgtZERO) {
      return;
    }
    // ________________________________________________
    if (!doprocessSameEventHCascades) {
      histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), collision.centFT0M()}));
      histos.fill(HIST("EventQA/hMult"), collision.centFT0M());
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
      histos.fill(HIST("EventQA/hMultFT0vsTPC"), collision.centFT0M(), collision.multNTracksPVeta1());
    }
    // Do basic QA

    if (applyEfficiencyCorrection) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initEfficiencyFromCCDB(bc);
    }
    TH2F* hEfficiencyV0[3];
    hEfficiencyV0[0] = hEfficiencyK0Short;
    hEfficiencyV0[1] = hEfficiencyLambda;
    hEfficiencyV0[2] = hEfficiencyAntiLambda;

    for (auto const& v0 : associatedV0s) {
      auto v0Data = v0.v0Core_as<V0DatasWithoutTrackX>();

      //---] track quality check [---
      auto postrack = v0Data.posTrack_as<TracksComplete>();
      auto negtrack = v0Data.negTrack_as<TracksComplete>();
      if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
        continue;

      //---] syst cuts [---
      if (v0Data.v0radius() < systCuts.v0RadiusMin || v0Data.v0radius() > systCuts.v0RadiusMax ||
          std::abs(v0Data.dcapostopv()) < systCuts.dcapostopv || std::abs(v0Data.dcanegtopv()) < systCuts.dcanegtopv ||
          v0Data.v0cosPA() < systCuts.v0cospa || v0Data.dcaV0daughters() > systCuts.dcaV0dau)
        continue;

      static_for<0, 2>([&](auto i) {
        constexpr int index = i.value;
        float efficiency = 1.0f;
        if (applyEfficiencyCorrection) {
          efficiency = hEfficiencyV0[index]->Interpolate(v0Data.pt(), v0Data.eta());
        }
        float weight = applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        if (v0.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || v0.mcTrue(index)) && (!doAssocPhysicalPrimary || v0.mcPhysicalPrimary()) && (!applyEfficiencyCorrection || efficiency != 0)) {
          if (-massWindowConfigurations.maxPeakNSigma < v0.invMassNSigma(index) && v0.invMassNSigma(index) < +massWindowConfigurations.maxPeakNSigma)
            histos.fill(HIST("h") + HIST(v0names[index]) + HIST("EtaVsPtVsPhi"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
          if ((-massWindowConfigurations.maxBgNSigma < v0.invMassNSigma(index) && v0.invMassNSigma(index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < v0.invMassNSigma(index) && v0.invMassNSigma(index) < +massWindowConfigurations.maxBgNSigma))
            histos.fill(HIST("h") + HIST(v0names[index]) + HIST("EtaVsPtVsPhiBg"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
          if (bitcheck(doCorrelation, index)) {
            histos.fill(HIST("h3d") + HIST(v0names[index]) + HIST("Spectrum"), v0Data.pt(), collision.centFT0M(), v0.invMassNSigma(index), weight);
            if (std::abs(v0Data.rapidity(index)) < 0.5) {
              histos.fill(HIST("h3d") + HIST(v0names[index]) + HIST("SpectrumY"), v0Data.pt(), collision.centFT0M(), v0.invMassNSigma(index), weight);
            }
          }
        }
      });
    }
    if (!doprocessSameEventHCascades) {
      for (auto const& triggerTrack : triggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track))
          continue;
        histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
          continue;
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
      }
    }

    // ________________________________________________
    // Do hadron - V0 correlations
    fillCorrelationsV0(triggerTracks, associatedV0s, false, collision.posZ(), collision.centFT0M());
  }

  void processSameEventHCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision,
                                 aod::AssocV0s const&, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                                 V0DatasWithoutTrackX const&, aod::V0sLinked const&, aod::CascDatas const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    // ________________________________________________
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    if (TMath::Abs(collision.posZ()) > zVertexCut) {
      return;
    }
    if (collision.centFT0M() > axisRanges[5][1] || collision.centFT0M() < axisRanges[5][0]) {
      return;
    }
    if (!collision.isInelGt0() && selectINELgtZERO) {
      return;
    }
    // ________________________________________________
    histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), collision.centFT0M()}));
    histos.fill(HIST("EventQA/hMult"), collision.centFT0M());
    histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    // Do basic QA
    if (applyEfficiencyCorrection) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initEfficiencyFromCCDB(bc);
    }
    TH2F* hEfficiencyCascade[4];
    hEfficiencyCascade[0] = hEfficiencyXiMinus;
    hEfficiencyCascade[1] = hEfficiencyXiPlus;
    hEfficiencyCascade[2] = hEfficiencyOmegaMinus;
    hEfficiencyCascade[3] = hEfficiencyOmegaPlus;

    for (auto const& casc : associatedCascades) {
      auto cascData = casc.cascData();

      //---] syst cuts [---
      if (std::abs(cascData.dcapostopv()) < systCuts.dcapostopv ||
          std::abs(cascData.dcanegtopv()) < systCuts.dcanegtopv ||
          cascData.dcabachtopv() < systCuts.casc_dcabachtopv ||
          cascData.dcaV0daughters() > systCuts.dcaV0dau ||
          cascData.dcacascdaughters() > systCuts.casc_dcacascdau ||
          cascData.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < systCuts.v0cospa ||
          cascData.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < systCuts.casc_cospa ||
          cascData.cascradius() < systCuts.casc_cascradius ||
          std::abs(cascData.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < systCuts.casc_mindcav0topv ||
          std::abs(cascData.mLambda() - pdgDB->Mass(3122)) > systCuts.casc_v0masswindow)
        continue;

      //---] track quality check [---
      auto postrack = cascData.posTrack_as<TracksComplete>();
      auto negtrack = cascData.negTrack_as<TracksComplete>();
      auto bachtrack = cascData.bachelor_as<TracksComplete>();
      if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || bachtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
        continue;

      static_for<0, 3>([&](auto i) {
        constexpr int index = i.value;
        float efficiency = 1.0f;
        if (applyEfficiencyCorrection) {
          efficiency = hEfficiencyCascade[index]->Interpolate(cascData.pt(), cascData.eta());
        }
        float weight = applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        if (casc.compatible(index, systCuts.dEdxCompatibility) && (!doMCassociation || casc.mcTrue(index)) && (!doAssocPhysicalPrimary || casc.mcPhysicalPrimary()) && (!applyEfficiencyCorrection || efficiency != 0)) {
          if (bitcheck(doCorrelation, index + 3)) {
            histos.fill(HIST("h3d") + HIST(cascadenames[index]) + HIST("Spectrum"), cascData.pt(), collision.centFT0M(), casc.invMassNSigma(index), weight);
            if (std::abs(cascData.rapidity(index)) < 0.5) {
              histos.fill(HIST("h3d") + HIST(cascadenames[index]) + HIST("SpectrumY"), cascData.pt(), collision.centFT0M(), casc.invMassNSigma(index), weight);
            }
          }
          if (-massWindowConfigurations.maxPeakNSigma < casc.invMassNSigma(index) && casc.invMassNSigma(index) < +massWindowConfigurations.maxPeakNSigma)
            histos.fill(HIST("h") + HIST(cascadenames[index]) + HIST("EtaVsPtVsPhi"), cascData.pt(), cascData.eta(), cascData.phi(), weight);
          if ((-massWindowConfigurations.maxBgNSigma < casc.invMassNSigma(index) && casc.invMassNSigma(index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < casc.invMassNSigma(index) && casc.invMassNSigma(index) < +massWindowConfigurations.maxBgNSigma))
            histos.fill(HIST("h") + HIST(cascadenames[index]) + HIST("EtaVsPtVsPhiBg"), cascData.pt(), cascData.eta(), cascData.phi(), weight);
        }
      });
    }
    for (auto const& triggerTrack : triggerTracks) {
      auto track = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(track))
        continue;
      histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
      histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
    }

    // ________________________________________________
    // Do hadron - cascade correlations
    fillCorrelationsCascade(triggerTracks, associatedCascades, false, collision.posX(), collision.posY(), collision.posZ(), collision.centFT0M());
  }
  void processSameEventHPions(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision,
                              aod::AssocPions const& associatedPions, aod::TriggerTracks const& triggerTracks,
                              TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    // ________________________________________________
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    if (TMath::Abs(collision.posZ()) > zVertexCut) {
      return;
    }
    if (collision.centFT0M() > axisRanges[5][1] || collision.centFT0M() < axisRanges[5][0]) {
      return;
    }
    if (!collision.isInelGt0() && selectINELgtZERO) {
      return;
    }
    // ________________________________________________
    if (!doprocessSameEventHCascades && !doprocessSameEventHV0s) {
      histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), collision.centFT0M()}));
      histos.fill(HIST("EventQA/hMult"), collision.centFT0M());
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    }
    // Do basic QA
    for (auto const& pion : associatedPions) {
      auto pionTrack = pion.track_as<TracksComplete>();
      histos.fill(HIST("hPionEtaVsPtVsPhi"), pionTrack.pt(), pionTrack.eta(), pionTrack.phi());
    }
    if (!doprocessSameEventHCascades && !doprocessSameEventHV0s) {
      for (auto const& triggerTrack : triggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track))
          continue;
        histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
          continue;
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
      }
    }

    // ________________________________________________
    // Do hadron - Pion correlations
    fillCorrelationsPion(triggerTracks, associatedPions, false, collision.posZ(), collision.centFT0M());
  }
  void processMixedEventHV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions,
                             aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks,
                             V0DatasWithoutTrackX const&, aod::V0sLinked const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, mixingParameter, -1, collisions, collisions)) {
      // ________________________________________________
      if (applyEfficiencyCorrection) {
        auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
        initEfficiencyFromCCDB(bc);
      }
      // Perform basic event selection on both collisions
      if (!collision1.sel8() || !collision2.sel8())
        continue;
      if (TMath::Abs(collision1.posZ()) > zVertexCut || TMath::Abs(collision2.posZ()) > zVertexCut)
        continue;
      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0])
        continue;
      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0])
        continue;
      if ((!collision1.isInelGt0() || !collision2.isInelGt0()) && selectINELgtZERO)
        continue;

      if (!doprocessMixedEventHCascades) {
        if (collision1.globalIndex() == collision2.globalIndex()) {
          histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
        }
        histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
        histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
        histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0M()}));
      }
      // ________________________________________________
      // Do slicing
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
      auto slicedAssocV0s = associatedV0s.sliceBy(collisionSliceV0s, collision2.globalIndex());
      // ________________________________________________
      // Do hadron - V0 correlations
      fillCorrelationsV0(slicedTriggerTracks, slicedAssocV0s, true, collision1.posZ(), collision1.centFT0M());
    }
  }
  void processMixedEventHCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions,
                                  aod::AssocV0s const&, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                                  V0DatasWithoutTrackX const&, aod::V0sLinked const&, aod::CascDatas const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, mixingParameter, -1, collisions, collisions)) {
      // ________________________________________________
      if (applyEfficiencyCorrection) {
        auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
        initEfficiencyFromCCDB(bc);
      }
      // Perform basic event selection on both collisions
      if (!collision1.sel8() || !collision2.sel8())
        continue;
      if (TMath::Abs(collision1.posZ()) > zVertexCut || TMath::Abs(collision2.posZ()) > zVertexCut)
        continue;
      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0])
        continue;
      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0])
        continue;
      if ((!collision1.isInelGt0() || !collision2.isInelGt0()) && selectINELgtZERO)
        continue;

      if (collision1.globalIndex() == collision2.globalIndex()) {
        histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
      }

      histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
      histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
      histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0M()}));
      // ________________________________________________
      // Do slicing
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
      auto slicedAssocCascades = associatedCascades.sliceBy(collisionSliceCascades, collision2.globalIndex());
      // ________________________________________________
      // Do hadron - cascade correlations
      fillCorrelationsCascade(slicedTriggerTracks, slicedAssocCascades, true, collision1.posX(), collision1.posY(), collision1.posZ(), collision1.centFT0M());
    }
  }
  void processMixedEventHPions(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions,
                               aod::AssocPions const& assocPions, aod::TriggerTracks const& triggerTracks,
                               TracksComplete const&)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, mixingParameter, -1, collisions, collisions)) {
      // ________________________________________________
      // Perform basic event selection on both collisions
      if (!collision1.sel8() || !collision2.sel8())
        continue;
      if (TMath::Abs(collision1.posZ()) > zVertexCut || TMath::Abs(collision2.posZ()) > zVertexCut)
        continue;
      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0])
        continue;
      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0])
        continue;
      if ((!collision1.isInelGt0() || !collision2.isInelGt0()) && selectINELgtZERO)
        continue;

      if (collision1.globalIndex() == collision2.globalIndex()) {
        histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
      }

      histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
      histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
      histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0M()}));
      // ________________________________________________
      // Do slicing
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
      auto slicedAssocPions = assocPions.sliceBy(collisionSlicePions, collision2.globalIndex());
      // ________________________________________________
      // Do hadron - cascade correlations
      fillCorrelationsPion(slicedTriggerTracks, slicedAssocPions, true, collision1.posZ(), collision1.centFT0M());
    }
  }

  void processMCGenerated(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>> const& collisions, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("hClosureTestEventCounter"), 2.5f);

    for (auto const& mcParticle : mcParticles) {
      Double_t geta = mcParticle.eta();
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      Double_t gpt = mcParticle.pt();
      if (abs(mcParticle.pdgCode()) == 211 || abs(mcParticle.pdgCode()) == 321 || abs(mcParticle.pdgCode()) == 2212 || abs(mcParticle.pdgCode()) == 11 || abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), gpt, 0.0f); // step 1: before all selections
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (abs(mcParticle.pdgCode()) == 310 && doCorrelationK0Short) {
          histos.fill(HIST("hGeneratedQAPtAssociatedK0"), gpt, 0.0f); // step 1: before all selections
        }
      }
    }

    for (auto const& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary())
        continue;
      if (abs(mcParticle.pdgCode()) == 211)
        histos.fill(HIST("Generated/hPion"), mcParticle.pt(), mcParticle.eta());
      if (abs(mcParticle.pdgCode()) == 310)
        histos.fill(HIST("Generated/hK0Short"), mcParticle.pt(), mcParticle.eta());
      if (mcParticle.pdgCode() == 3122)
        histos.fill(HIST("Generated/hLambda"), mcParticle.pt(), mcParticle.eta());
      if (mcParticle.pdgCode() == -3122)
        histos.fill(HIST("Generated/hAntiLambda"), mcParticle.pt(), mcParticle.eta());
      if (mcParticle.pdgCode() == 3312)
        histos.fill(HIST("Generated/hXiMinus"), mcParticle.pt(), mcParticle.eta());
      if (mcParticle.pdgCode() == -3312)
        histos.fill(HIST("Generated/hXiPlus"), mcParticle.pt(), mcParticle.eta());
      if (mcParticle.pdgCode() == 3334)
        histos.fill(HIST("Generated/hOmegaMinus"), mcParticle.pt(), mcParticle.eta());
      if (mcParticle.pdgCode() == -3334)
        histos.fill(HIST("Generated/hOmegaPlus"), mcParticle.pt(), mcParticle.eta());
    }
    if (collisions.size() < 1)
      return;

    // determine best collision properties
    int biggestNContribs = -1;
    int bestCollisionFT0Mpercentile = -1;
    float bestCollisionVtxZ = 0.0f;
    bool bestCollisionSel8 = false;
    bool bestCollisionINELgtZERO = false;
    for (auto& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCollisionFT0Mpercentile = collision.centFT0M();
        bestCollisionSel8 = collision.sel8();
        bestCollisionVtxZ = collision.posZ();
        bestCollisionINELgtZERO = collision.isInelGt0();
      }
    }

    if (collisions.size() > 1) {
      for (auto const& mcParticle : mcParticles) {
        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (abs(mcParticle.y()) < 0.5) {
          if (abs(mcParticle.pdgCode()) == 211)
            histos.fill(HIST("GeneratedWithPV/hPion_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
          if (abs(mcParticle.pdgCode()) == 310)
            histos.fill(HIST("GeneratedWithPV/hK0Short_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
          if (mcParticle.pdgCode() == 3122)
            histos.fill(HIST("GeneratedWithPV/hLambda_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
          if (mcParticle.pdgCode() == -3122)
            histos.fill(HIST("GeneratedWithPV/hAntiLambda_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
          if (mcParticle.pdgCode() == 3312)
            histos.fill(HIST("GeneratedWithPV/hXiMinus_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
          if (mcParticle.pdgCode() == -3312)
            histos.fill(HIST("GeneratedWithPV/hXiPlus_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
          if (mcParticle.pdgCode() == 3334)
            histos.fill(HIST("GeneratedWithPV/hOmegaMinus_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
          if (mcParticle.pdgCode() == -3334)
            histos.fill(HIST("GeneratedWithPV/hOmegaPlus_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
        }
      }
    }

    // do selections on best collision
    // WARNING: if 2 PV case large, this will not necessarily be fine!
    //          caution advised!
    if (!bestCollisionSel8)
      return;
    if (std::abs(bestCollisionVtxZ) > 10.0f)
      return;
    if (!bestCollisionINELgtZERO)
      return;

    histos.fill(HIST("hClosureTestEventCounter"), 3.5f);

    for (auto const& mcParticle : mcParticles) {
      Double_t geta = mcParticle.eta();
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      Double_t gpt = mcParticle.pt();
      if (abs(mcParticle.pdgCode()) == 211 || abs(mcParticle.pdgCode()) == 321 || abs(mcParticle.pdgCode()) == 2212 || abs(mcParticle.pdgCode()) == 11 || abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), gpt, 1.0f); // step 2: after event selection
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (abs(mcParticle.pdgCode()) == 310 && doCorrelationK0Short) {
          histos.fill(HIST("hGeneratedQAPtAssociatedK0"), gpt, 1.0f); // step 2: before all selections
        }
      }
    }

    for (auto const& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      Double_t geta = mcParticle.eta();
      Double_t gpt = mcParticle.pt();
      if (abs(mcParticle.pdgCode()) == 211 || abs(mcParticle.pdgCode()) == 321 || abs(mcParticle.pdgCode()) == 2212 || abs(mcParticle.pdgCode()) == 11 || abs(mcParticle.pdgCode()) == 13)
        histos.fill(HIST("GeneratedWithPV/hTrigger"), gpt, geta);
      if (abs(mcParticle.pdgCode()) == 211)
        histos.fill(HIST("GeneratedWithPV/hPion"), gpt, geta);
      if (abs(mcParticle.pdgCode()) == 310)
        histos.fill(HIST("GeneratedWithPV/hK0Short"), gpt, geta);
      if (mcParticle.pdgCode() == 3122)
        histos.fill(HIST("GeneratedWithPV/hLambda"), gpt, geta);
      if (mcParticle.pdgCode() == -3122)
        histos.fill(HIST("GeneratedWithPV/hAntiLambda"), gpt, geta);
      if (mcParticle.pdgCode() == 3312)
        histos.fill(HIST("GeneratedWithPV/hXiMinus"), gpt, geta);
      if (mcParticle.pdgCode() == -3312)
        histos.fill(HIST("GeneratedWithPV/hXiPlus"), gpt, geta);
      if (mcParticle.pdgCode() == 3334)
        histos.fill(HIST("GeneratedWithPV/hOmegaMinus"), gpt, geta);
      if (mcParticle.pdgCode() == -3334)
        histos.fill(HIST("GeneratedWithPV/hOmegaPlus"), gpt, geta);

      if (abs(mcParticle.y()) < 0.5) {
        if (abs(mcParticle.pdgCode()) == 211)
          histos.fill(HIST("GeneratedWithPV/hPion_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
        if (abs(mcParticle.pdgCode()) == 310)
          histos.fill(HIST("GeneratedWithPV/hK0Short_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
        if (mcParticle.pdgCode() == 3122)
          histos.fill(HIST("GeneratedWithPV/hLambda_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
        if (mcParticle.pdgCode() == -3122)
          histos.fill(HIST("GeneratedWithPV/hAntiLambda_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
        if (mcParticle.pdgCode() == 3312)
          histos.fill(HIST("GeneratedWithPV/hXiMinus_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
        if (mcParticle.pdgCode() == -3312)
          histos.fill(HIST("GeneratedWithPV/hXiPlus_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
        if (mcParticle.pdgCode() == 3334)
          histos.fill(HIST("GeneratedWithPV/hOmegaMinus_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
        if (mcParticle.pdgCode() == -3334)
          histos.fill(HIST("GeneratedWithPV/hOmegaPlus_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
      }
    }
  }
  void processClosureTest(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>> const& recCollisions, aod::McParticles const& mcParticles)
  {

    std::vector<uint32_t> triggerIndices;
    std::vector<std::vector<uint32_t>> associatedIndices;
    std::vector<uint32_t> piIndices;
    std::vector<uint32_t> k0ShortIndices;
    std::vector<uint32_t> lambdaIndices;
    std::vector<uint32_t> antiLambdaIndices;
    std::vector<uint32_t> xiMinusIndices;
    std::vector<uint32_t> xiPlusIndices;
    std::vector<uint32_t> omegaMinusIndices;
    std::vector<uint32_t> omegaPlusIndices;

    for (auto const& mcParticle : mcParticles) {
      Double_t geta = mcParticle.eta();
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      Double_t gpt = mcParticle.pt();
      if (abs(mcParticle.pdgCode()) == 211 || abs(mcParticle.pdgCode()) == 321 || abs(mcParticle.pdgCode()) == 2212 || abs(mcParticle.pdgCode()) == 11 || abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 0.0f); // step 1: no event selection whatsoever
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (abs(mcParticle.pdgCode()) == 310 && doCorrelationK0Short) {
          histos.fill(HIST("hClosureQAPtAssociatedK0"), gpt, 0.0f); // step 1: no event selection whatsoever
        }
      }
    }

    histos.fill(HIST("hClosureTestEventCounter"), 0.5f);

    int bestCollisionFT0Mpercentile = -1;
    float bestCollisionVtxZ = 0.0f;
    bool bestCollisionSel8 = false;
    bool bestCollisionINELgtZERO = false;
    int biggestNContribs = -1;
    for (auto& recCollision : recCollisions) {
      if (biggestNContribs < recCollision.numContrib()) {
        biggestNContribs = recCollision.numContrib();
        bestCollisionFT0Mpercentile = recCollision.centFT0M();
        bestCollisionSel8 = recCollision.sel8();
        bestCollisionVtxZ = recCollision.posZ();
        bestCollisionINELgtZERO = recCollision.isInelGt0();
      }
    }
    if (doGenEventSelection) {
      if (!bestCollisionSel8)
        return;
      if (std::abs(bestCollisionVtxZ) > zVertexCut)
        return;
      if (!bestCollisionINELgtZERO)
        return;
      if (bestCollisionFT0Mpercentile > axisRanges[5][1] || bestCollisionFT0Mpercentile < axisRanges[5][0]) {
        return;
      }
    }

    histos.fill(HIST("hClosureTestEventCounter"), 1.5f);

    for (auto const& mcParticle : mcParticles) {
      Double_t geta = mcParticle.eta();
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      Double_t gpt = mcParticle.pt();
      if (abs(mcParticle.pdgCode()) == 211 || abs(mcParticle.pdgCode()) == 321 || abs(mcParticle.pdgCode()) == 2212 || abs(mcParticle.pdgCode()) == 11 || abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 1.0f); // step 2: after event selection
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (abs(mcParticle.pdgCode()) == 310 && doCorrelationK0Short) {
          histos.fill(HIST("hClosureQAPtAssociatedK0"), gpt, 1.0f); // step 2: after event selection
        }
      }
    }

    int iteratorNum = -1;
    for (auto const& mcParticle : mcParticles) {
      iteratorNum = iteratorNum + 1;
      Double_t geta = mcParticle.eta();
      Double_t gpt = mcParticle.pt();
      Double_t gphi = mcParticle.phi();
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      if (abs(mcParticle.pdgCode()) == 211 || abs(mcParticle.pdgCode()) == 321 || abs(mcParticle.pdgCode()) == 2212 || abs(mcParticle.pdgCode()) == 11 || abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          triggerIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hTrigger"), gpt, geta, bestCollisionFT0Mpercentile);
        }
      }
      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (abs(mcParticle.pdgCode()) == 211 && doCorrelationPion) {
          piIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hPion"), gpt, geta, gphi);
        }
        if (abs(mcParticle.pdgCode()) == 310 && doCorrelationK0Short) {
          k0ShortIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hK0Short"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == 3122 && doCorrelationLambda) {
          lambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == -3122 && doCorrelationAntiLambda) {
          antiLambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hAntiLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == 3312 && doCorrelationXiMinus) {
          xiMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hXiMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == -3312 && doCorrelationXiPlus) {
          xiPlusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hXiPlus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == 3334 && doCorrelationOmegaMinus) {
          omegaMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hOmegaMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == -3334 && doCorrelationOmegaPlus) {
          omegaPlusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hOmegaPlus"), gpt, geta, gphi);
        }
      }
    }
    associatedIndices.emplace_back(piIndices);
    associatedIndices.emplace_back(k0ShortIndices);
    associatedIndices.emplace_back(lambdaIndices);
    associatedIndices.emplace_back(antiLambdaIndices);
    associatedIndices.emplace_back(xiMinusIndices);
    associatedIndices.emplace_back(xiPlusIndices);
    associatedIndices.emplace_back(omegaMinusIndices);
    associatedIndices.emplace_back(omegaPlusIndices);

    for (Int_t iTrigger = 0; iTrigger < triggerIndices.size(); iTrigger++) {
      auto triggerParticle = mcParticles.iteratorAt(triggerIndices[iTrigger]);
      // check range of trigger particle
      if (triggerParticle.pt() > axisRanges[3][1] || triggerParticle.pt() < axisRanges[3][0]) {
        continue;
      }
      Double_t getatrigger = triggerParticle.eta();
      Double_t gphitrigger = triggerParticle.phi();
      Double_t pttrigger = triggerParticle.pt();
      auto const& mother = triggerParticle.mothers_first_as<aod::McParticles>();
      Int_t globalIndex = mother.globalIndex();
      static_for<0, 7>([&](auto i) { // associated loop
        constexpr int index = i.value;
        for (Int_t iassoc = 0; iassoc < associatedIndices[index].size(); iassoc++) {
          auto assocParticle = mcParticles.iteratorAt(associatedIndices[index][iassoc]);
          if (triggerIndices[iTrigger] != associatedIndices[index][iassoc] && globalIndex != assocParticle.globalIndex()) { // avoid self
            Double_t getaassoc = assocParticle.eta();
            Double_t gphiassoc = assocParticle.phi();
            Double_t ptassoc = assocParticle.pt();
            Double_t deltaphi = ComputeDeltaPhi(gphitrigger, gphiassoc);
            Double_t deltaeta = getatrigger - getaassoc;

            // skip if basic ranges not met
            if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
              continue;
            if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
              continue;
            if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
              continue;
            histos.fill(HIST("ClosureTest/sameEvent/") + HIST(particlenames[index]), ComputeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, bestCollisionVtxZ, bestCollisionFT0Mpercentile);
          }
        }
      });
    }
  }

  PROCESS_SWITCH(correlateStrangeness, processSameEventHV0s, "Process same events, h-V0s", true);
  PROCESS_SWITCH(correlateStrangeness, processSameEventHCascades, "Process same events, h-Cascades", true);
  PROCESS_SWITCH(correlateStrangeness, processSameEventHPions, "Process same events, h-Pion", true);
  PROCESS_SWITCH(correlateStrangeness, processMixedEventHV0s, "Process mixed events, h-V0s", true);
  PROCESS_SWITCH(correlateStrangeness, processMixedEventHCascades, "Process mixed events, h-Cascades", true);
  PROCESS_SWITCH(correlateStrangeness, processMixedEventHPions, "Process mixed events, h-Pion", true);
  PROCESS_SWITCH(correlateStrangeness, processMCGenerated, "Process MC generated", true);
  PROCESS_SWITCH(correlateStrangeness, processClosureTest, "Process Closure Test", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<correlateStrangeness>(cfgc)};
}
