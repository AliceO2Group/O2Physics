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
/// \brief this is a starting point for the Resonances tutorial
/// \author sourav kundu
/// \since 02/11/2023

#include "PWGLF/DataModel/ReducedF1ProtonTables.h"

#include "Common/Core/trackUtilities.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjString.h>

#include <fairlogger/Logger.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct f1protoncorrelation {

  double bz = 0.;
  double bz2 = 0.;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // PID selection
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> typeofCombined{"typeofCombined", 1, "type of combined"};
  // PID selection
  Configurable<bool> fillSparse{"fillSparse", 1, "Fill Sparse"};
  Configurable<bool> doold{"doold", 0, "Do old calculation"};
  Configurable<bool> fillRotation{"fillRotation", 1, "Fill rotation"};
  Configurable<bool> pdepPID{"pdepPID", 1, "Momentum dependent pi, k PID"};
  Configurable<int> strategyPIDPion{"strategyPIDPion", 0, "PID strategy Pion"};
  Configurable<int> strategyPIDKaon{"strategyPIDKaon", 0, "PID strategy Kaon"};
  Configurable<float> maxKKS0Mass{"maxKKS0Mass", 1.025, "Maximum kaon kshort mass"};
  Configurable<float> maxMomentumPion{"maxMomentumPion", 4.0, "Maximum momentum Pion"};
  Configurable<float> maxMomentumKaon{"maxMomentumKaon", 4.0, "Maximum momentum Kaon"};
  Configurable<float> momentumTOFPionMin{"momentumTOFPionMin", 0.8, "Pion momentum TOF Min"};
  Configurable<float> momentumTOFKaonMin{"momentumTOFKaonMin", 0.5, "Kaon momentum TOF Min"};
  Configurable<float> momentumTOFPionMax{"momentumTOFPionMax", 1.2, "Pion momentum TOF Max"};
  Configurable<float> momentumTOFKaonMax{"momentumTOFKaonMax", 0.9, "Kaon momentum TOF Max"};
  Configurable<float> momentumTOFProton{"momentumTOFProton", 0.7, "Proton momentum TOF"};
  Configurable<float> momentumProtonMax{"momentumProtonMax", 3.0, "Maximum proton momentum"};
  Configurable<float> momentumProtonMin{"momentumProtonMin", 0.1, "Minimum proton momentum"};
  Configurable<float> lowPtF1{"lowPtF1", 1.0, "PT cut F1"};
  Configurable<int> nRot{"nRot", 4, "Number of rotational bkg"};
  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  Configurable<int> nEvtMixingBkg{"nEvtMixingBkg", 5, "Number of events to mix for background reconstruction"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 10.0}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0, 40.0, 80.0, 500.0}, "Mixing bins - number of contributor"};

  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {100, 1.0, 1.4}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisMultiplicity{"configThnAxisMultiplicity", {100, 0.0, 500}, "multiplicity"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisMt{"configThnAxisMt", {100, 0.0, 10.}, "#it{M}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisKstar{"configThnAxisKstar", {100, 0.0, 1.0}, "#it{k}^{*} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisPtProton{"configThnAxisPtProton", {20, 0.0, 4.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisNsigma{"configThnAxisNsigma", {90, -9.0, 9.0}, "NsigmaCombined"};
  ConfigurableAxis configThnAxisCharge{"configThnAxisCharge", {5, -2.5, 2.5}, "Charge"};

  // -------------------------
  // Systematics control
  // -------------------------
  Configurable<int> nSysRand{"nSysRand", 500, "Number of random systematic combinations (in addition to default sysId=0)"};
  Configurable<uint32_t> sysSeed{"sysSeed", 12345u, "Seed for systematics (reproducible)"};

  struct SysCuts {
    // Primary tracks (π/K/p)
    float maxDcaxy;
    float maxDcaz;
    int minTPCCrossedRows;
    int minTPCClusters;

    // V0 (K0s)
    float minCPA;
    float minRadius;
    float maxDcaDaughters;
    float maxDcaV0;
    float maxLifetime;
    float minDcaD1;
    float minDcaD2;

    // PID modes
    int pidPi; // 0/1/2 depending on what you use
    int pidK;  // 0/1/2
    int pidP;  // 0/1/2
  };

  std::vector<SysCuts> sysCuts; // sysCuts[0] = default; sysCuts[1..] random unique
  int nSysTotal = 1;

  // -------------------------
  // Build sysCuts (default + N random unique)
  // -------------------------
  void buildSystematicCuts()
  {
    // choices from your table/picture
    const std::array<float, 2> optDcaxy{0.01f, 0.03f};
    const std::array<float, 2> optDcaz{0.01f, 0.03f};
    const std::array<int, 2> optNcrs{90, 100};
    const std::array<int, 2> optNcls{90, 100};

    const std::array<float, 2> optCPA{0.99f, 0.995f};
    const std::array<float, 2> optRad{0.8f, 1.0f};
    const std::array<float, 2> optDcaDD{0.9f, 0.8f};
    const std::array<float, 2> optDcaV0{0.2f, 0.15f};
    const std::array<float, 2> optLife{16.f, 18.f};
    const std::array<float, 2> optDcaD1{0.06f, 0.08f};
    const std::array<float, 2> optDcaD2{0.06f, 0.08f};

    // bit layout:
    // primary: 4 bits  (0..3)
    // v0:      7 bits  (4..10)
    // pid:     3 bits  (11..13)
    // total: 14 bits -> 16384 combos
    auto buildFromMask = [&](uint32_t m) -> SysCuts {
      SysCuts c{};
      c.maxDcaxy = optDcaxy[(m >> 0) & 1u];
      c.maxDcaz = optDcaz[(m >> 1) & 1u];
      c.minTPCCrossedRows = optNcrs[(m >> 2) & 1u];
      c.minTPCClusters = optNcls[(m >> 3) & 1u];

      c.minCPA = optCPA[(m >> 4) & 1u];
      c.minRadius = optRad[(m >> 5) & 1u];
      c.maxDcaDaughters = optDcaDD[(m >> 6) & 1u];
      c.maxDcaV0 = optDcaV0[(m >> 7) & 1u];
      c.maxLifetime = optLife[(m >> 8) & 1u];
      c.minDcaD1 = optDcaD1[(m >> 9) & 1u];
      c.minDcaD2 = optDcaD2[(m >> 10) & 1u];

      c.pidPi = (m >> 11) & 1u;
      c.pidK = (m >> 12) & 1u;
      c.pidP = (m >> 13) & 1u;
      return c;
    };

    // DEFAULT (sysId=0): set to your baseline
    SysCuts def{};
    def.maxDcaxy = 0.05f;
    def.maxDcaz = 0.05f;
    def.minTPCCrossedRows = 80;
    def.minTPCClusters = 80;

    def.minCPA = 0.985f;
    def.minRadius = 0.5f;
    def.maxDcaDaughters = 1.0f;
    def.maxDcaV0 = 0.3f;
    def.maxLifetime = 20.f;
    def.minDcaD1 = 0.05f;
    def.minDcaD2 = 0.05f;

    def.pidPi = 0;
    def.pidK = 0;
    def.pidP = 0;

    std::vector<uint32_t> masks;
    masks.reserve(16384);
    for (uint32_t m = 0; m < 16384; ++m) {
      masks.push_back(m);
    }

    std::mt19937 rng(sysSeed);
    std::shuffle(masks.begin(), masks.end(), rng);

    sysCuts.clear();
    sysCuts.reserve(1 + nSysRand);
    sysCuts.push_back(def);

    const int nPick = std::min<int>(nSysRand, (int)masks.size());
    for (int i = 0, picked = 0; picked < nPick && i < (int)masks.size(); ++i) {
      if (masks[i] == 0u) { // avoid trivial mask that can reproduce default
        continue;
      }
      sysCuts.push_back(buildFromMask(masks[i]));
      ++picked;
    }
    nSysTotal = (int)sysCuts.size();
  }

  // -------------------------
  // Helpers for systematics selections
  // -------------------------
  inline bool passPrimary(float dcaxy, float dcaz, float ncrs, float ncls, const SysCuts& c) const
  {
    if (std::abs(dcaxy) > c.maxDcaxy)
      return false;
    if (std::abs(dcaz) > c.maxDcaz)
      return false;
    if ((int)ncrs < c.minTPCCrossedRows)
      return false;
    if ((int)ncls < c.minTPCClusters)
      return false;
    return true;
  }

  inline bool passV0(const aod::F1Tracks::iterator& t, const SysCuts& c) const
  {
    if (t.k0Cpa() < c.minCPA)
      return false;
    if (t.k0Radius() < c.minRadius)
      return false;
    if (t.k0DcaDaughters() > c.maxDcaDaughters)
      return false;
    if (t.k0Dca() > c.maxDcaV0)
      return false;
    if (t.k0LifeTime() > c.maxLifetime)
      return false;
    if (std::abs(t.k0D1Dcaxy()) < c.minDcaD1)
      return false;
    if (std::abs(t.k0D2Dcaxy()) < c.minDcaD2)
      return false;

    return true;
  }

  // PID dispatchers
  inline bool passPionPID(int pidMode,
                          const aod::F1Tracks::iterator& f1track,
                          const TLorentzVector& pion,
                          float maxMomPi) const
  {
    // pidMode:
    // 0 = default:   no TOF -> |TPC|<2.5 ; with TOF -> sqrt(TPC^2 + TOF^2)<2.5
    // 1 = syst-1:    same but use 2.0 instead of 2.5 everywhere
    // 2 = syst-2:    no TOF -> |TPC|<2.0 ; with TOF -> sqrt(TPC^2 + TOF^2)<2.5

    if (pion.Pt() > maxMomPi) {
      return false;
    }

    const float nsTPC = f1track.f1d1TPC();
    const float nsTOF = f1track.pionTOF(); // TOF nσ stored in reduced table

    const float cutNoTOF = (pidMode == 0) ? 2.5f : 2.0f;   // mode0:2.5, mode1:2.0, mode2:2.0
    const float cutWithTOF = (pidMode == 1) ? 2.0f : 2.5f; // mode0:2.5, mode1:2.0, mode2:2.5

    if (f1track.f1d1TOFHit() != 1) {
      return (std::abs(nsTPC) < cutNoTOF);
    }

    const float comb = std::sqrt(nsTPC * nsTPC + nsTOF * nsTOF);
    return (comb < cutWithTOF);
  }
  inline bool passKaonPID(int pidMode,
                          const aod::F1Tracks::iterator& f1track,
                          const TLorentzVector& kaon,
                          float maxMomK) const
  {
    // pidMode:
    // 0 = default:   no TOF -> momentum-dependent TPC cut (original ladder) with TOF -> circular cut sqrt(TPC^2 + TOF^2) < 2.5
    // 1 = syst-1:    same, but use 2.0 instead of 2.5 everywhere (including ladder upper edges)
    // 2 = syst-2:    no TOF -> ladder with "tight" = 2.0 (upper edges also 2.0)
    //                with TOF -> circular cut < 2.5

    if (kaon.Pt() > maxMomK) {
      return false;
    }

    const float nsTPC = f1track.f1d2TPC();
    const float nsTOF = f1track.kaonTOF(); // TOF nσ stored in reduced table

    const float cutNoTOFBase = (pidMode == 0) ? 2.5f : 2.0f; // mode0:2.5, mode1/2:2.0
    const float cutWithTOF = (pidMode == 1) ? 2.0f : 2.5f;   // mode0:2.5, mode1:2.0, mode2:2.5

    // --- no TOF: momentum-dependent TPC ladder
    if (f1track.f1d2TOFHit() != 1) {
      const float pt = kaon.Pt();

      if (pt <= 0.5f) {
        if (nsTPC < -cutNoTOFBase || nsTPC > cutNoTOFBase)
          return false;
      } else if (pt <= 0.7f) {
        const float low = -1.5f;
        const float high = cutNoTOFBase; // original high=+2.5
        if (nsTPC < low || nsTPC > high)
          return false;
      } else if (pt <= 1.0f) {
        const float low = -1.0f;
        const float high = cutNoTOFBase; // original high=+2.5
        if (nsTPC < low || nsTPC > high)
          return false;
      } else {
        if (nsTPC < -cutNoTOFBase || nsTPC > cutNoTOFBase)
          return false;
      }
      return true;
    }

    // --- TOF available: circular cut in (TPC,TOF) nσ plane
    const float comb = std::sqrt(nsTPC * nsTPC + nsTOF * nsTOF);
    return (comb < cutWithTOF);
  }
  inline bool passProtonPID(int pidMode,
                            const aod::ProtonTracks::iterator& ptrack,
                            const TLorentzVector& proton,
                            float pMin,
                            float pMax,
                            float pTofThr) const
  {
    // pidMode:
    // 0 = default:  p < thr -> |TPC| < 2.5   ;  p >= thr -> TOF mandatory AND circular(TPC,TOF) < 2.0
    // 1 = syst-1:   p < thr -> |TPC| < 2.0   ;  p >= thr -> TOF mandatory AND circular(TPC,TOF) < 2.0
    // 2 = syst-2:   p < thr -> |TPC| < 2.5   ;  p >= thr -> TOF mandatory AND circular(TPC,TOF) < 2.5

    if (proton.Pt() > pMax || proton.Pt() < pMin) {
      return false;
    }

    const float cutTPC =
      (pidMode == 1) ? 2.0f : 2.5f; // mode0:2.5, mode1:2.0, mode2:2.5
    const float cutCircle =
      (pidMode == 2) ? 2.0f : 2.5f; // mode0:2.5, mode1:2.5, mode2:2.0

    // below threshold: TPC-only
    if (proton.P() < pTofThr) {
      return (std::abs(ptrack.protonNsigmaTPC()) < cutTPC);
    }

    // above threshold: TOF must exist
    if (ptrack.protonTOFHit() != 1) {
      return false;
    }

    // circular cut in (TPC,TOF)
    const float nsTPC = ptrack.protonNsigmaTPC();
    const float nsTOF = ptrack.protonNsigmaTOF();
    const float comb = std::sqrt(nsTPC * nsTPC + nsTOF * nsTOF);
    return (comb < cutCircle);
  }

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    buildSystematicCuts();

    nSysTotal = (int)sysCuts.size(); // or whatever vector you fill
    LOGF(info, "sysCuts.size()=%zu  nSysTotal=%d", sysCuts.size(), nSysTotal);
    const AxisSpec thnAxisSys{nSysTotal, -0.5f, float(nSysTotal) - 0.5f, "sysId"};
    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisMt{configThnAxisMt, "#it{M}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtProton{configThnAxisPtProton, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisKstar{configThnAxisKstar, "#it{k}^{*} (GeV/#it{c})"};
    const AxisSpec thnAxisNsigma{configThnAxisNsigma, "NsigmaCombined"};
    const AxisSpec thnAxisCharge{configThnAxisCharge, "Charge"};
    const AxisSpec thnAxisMultiplicity{configThnAxisMultiplicity, "Multiplicity"};

    // register histograms
    histos.add("hPhaseSpaceProtonKaonSame", "hPhaseSpaceProtonKaonSame", kTH3F, {{100, -0.3f, 0.3f}, {100, -0.3, 0.3}, {10, 0.0, 1.0}});
    histos.add("hPhaseSpaceProtonPionSame", "hPhaseSpaceProtonPionSame", kTH3F, {{100, -0.3f, 0.3f}, {100, -0.3, 0.3}, {10, 0.0, 1.0}});
    histos.add("hPhaseSpaceProtonKaonMix", "hPhaseSpaceProtonKaonMix", kTH3F, {{100, -0.3f, 0.3f}, {100, -0.3, 0.3}, {10, 0.0, 1.0}});
    histos.add("hPhaseSpaceProtonPionMix", "hPhaseSpaceProtonPionMix", kTH3F, {{100, -0.3f, 0.3f}, {100, -0.3, 0.3}, {10, 0.0, 1.0}});

    histos.add("hNsigmaProtonTPC", "Nsigma Proton TPC distribution", kTH3F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}, {100, 0.0f, 10.0f}});
    histos.add("hNsigmaKaonTPC", "Nsigma Kaon TPC distribution", kTH2F, {{100, -5.0f, 5.0f}, {100, 0.0f, 10.0f}});
    histos.add("hNsigmaPionTPC", "Nsigma Pion TPC distribution", kTH2F, {{100, -5.0f, 5.0f}, {100, 0.0f, 10.0f}});
    histos.add("hNsigmaPionKaonTPC", "Nsigma Pion Kaon TPC correlation", kTH2F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}});
    histos.add("h2SameEventPtCorrelation", "Pt correlation of F1 and proton", kTH3F, {{100, 0.0f, 1.0f}, {100, 0.0, 10.0}, {100, 0.0, 10.0}});
    histos.add("h2SameEventf1pptCorrelation", "h2SameEventf1pptCorrelation", kTH3F, {{100, 1.0, 1.4}, {300, 0.0f, 3.0f}, {100, 0.0, 10.0}});
    histos.add("h2SameEventInvariantMassUnlike_mass_SYS", "SE unlike (SYS)", kTHnSparseF, {thnAxisSys, thnAxisKstar, thnAxisMt, thnAxisInvMass, thnAxisMultiplicity});
    histos.add("h2SameEventInvariantMassLike_mass_SYS", "SE like (SYS)", kTHnSparseF, {thnAxisSys, thnAxisKstar, thnAxisMt, thnAxisInvMass, thnAxisMultiplicity});
    histos.add("h2MixEventInvariantMassUnlike_mass_SYS", "ME unlike (SYS)", kTHnSparseF, {thnAxisSys, thnAxisKstar, thnAxisMt, thnAxisInvMass, thnAxisMultiplicity});
    histos.add("h2MixEventInvariantMassLike_mass_SYS", "ME like (SYS)", kTHnSparseF, {thnAxisSys, thnAxisKstar, thnAxisMt, thnAxisInvMass, thnAxisMultiplicity});

    if (doold) {
      histos.add("h2SameEventInvariantMassUnlike_mass", "Unlike Sign Invariant mass of f1 same event", kTHnSparseF, {thnAxisKstar, thnAxisPt, thnAxisInvMass, thnAxisCharge, thnAxisMultiplicity});
      histos.add("h2SameEventInvariantMassLike_mass", "Like Sign Invariant mass of f1 same event", kTHnSparseF, {thnAxisKstar, thnAxisPt, thnAxisInvMass, thnAxisCharge, thnAxisMultiplicity});
      histos.add("h2SameEventInvariantMassRot_mass", "Rotational Invariant mass of f1 same event", kTHnSparseF, {thnAxisKstar, thnAxisPt, thnAxisInvMass, thnAxisCharge});

      histos.add("h2MixEventInvariantMassUnlike_mass", "Unlike Sign Invariant mass of f1 mix event", kTHnSparseF, {thnAxisKstar, thnAxisPt, thnAxisInvMass, thnAxisCharge, thnAxisMultiplicity});
      histos.add("h2MixEventInvariantMassLike_mass", "Like Sign Invariant mass of f1 mix event", kTHnSparseF, {thnAxisKstar, thnAxisPt, thnAxisInvMass, thnAxisCharge, thnAxisMultiplicity});
      histos.add("h2MixEventInvariantMassRot_mass", "Rotational Sign Invariant mass of f1 mix event", kTHnSparseF, {thnAxisKstar, thnAxisPt, thnAxisInvMass, thnAxisCharge});

      histos.add("h2MixEventInvariantMassUnlike_mass_SEFP", "Unlike-sign invariant mass of f1 mix event (SE-F1P: π mixed, p same event)", kTHnSparseF, {thnAxisKstar, thnAxisPt, thnAxisInvMass, thnAxisCharge});
      histos.add("h2MixEventInvariantMassUnlike_mass_DEFP", "Unlike-sign invariant mass of f1 mix event (DE-F1P: π + p mixed)", kTHnSparseF, {thnAxisKstar, thnAxisPt, thnAxisInvMass, thnAxisCharge});
    }
    if (fillSparse) {
      histos.add("SEMassUnlike", "SEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPtProton, thnAxisKstar, thnAxisNsigma, thnAxisCharge});
      histos.add("SEMassLike", "SEMassLike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPtProton, thnAxisKstar, thnAxisNsigma, thnAxisCharge});
      histos.add("SEMassRot", "SEMassRot", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPtProton, thnAxisKstar, thnAxisNsigma, thnAxisCharge});

      histos.add("MEMassUnlike", "MEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPtProton, thnAxisKstar, thnAxisNsigma, thnAxisCharge});
      histos.add("MEMassLike", "MEMassLike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPtProton, thnAxisKstar, thnAxisNsigma, thnAxisCharge});
      histos.add("MEMassRot", "MEMassRot", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPtProton, thnAxisKstar, thnAxisNsigma, thnAxisCharge});
    }
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  int getMagneticField(uint64_t timestamp)
  {
    // Get the magnetic field
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("/GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %2.2f kG", timestamp, 0.1 * grpo->getNominalL3Field());
    }
    return 0.1 * grpo->getNominalL3Field();
  }

  /// Magnetic field to be provided in Tesla
  static constexpr float tmpRadiiTPC[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};
  float PhiAtSpecificRadiiTPC(const TLorentzVector part1, const TLorentzVector part2, float charge1 = 0, int charge2 = 0, float magfield1 = 0.0, float magfield2 = 0.0)
  {
    float pt1 = part1.Pt();
    float phi1 = part1.Phi();
    float value1 = 0.0;
    float count1 = 0.0;
    for (size_t i = 0; i < 9; i++) {
      auto arg1 = 0.3 * charge1 * magfield1 * tmpRadiiTPC[i] * 0.01 / (2. * pt1);
      if (std::fabs(arg1) < 1) {
        value1 = value1 + (phi1 - std::asin(arg1));
        count1 = count1 + 1.0;
      }
    }
    value1 = value1 / count1;

    float pt2 = part2.Pt();
    float phi2 = part2.Phi();
    float value2 = 0.0;
    float count2 = 0.0;
    for (size_t i = 0; i < 9; i++) {
      auto arg2 = 0.3 * charge2 * magfield2 * tmpRadiiTPC[i] * 0.01 / (2. * pt2);
      if (std::fabs(arg2) < 1) {
        value2 = value2 + (phi2 - std::asin(arg2));
        count2 = count2 + 1.0;
      }
    }
    value2 = value2 / count2;
    return value1 - value2;
  }

  // get kstar
  TLorentzVector trackSum, PartOneCMS, PartTwoCMS, trackRelK;
  float getkstar(const TLorentzVector part1,
                 const TLorentzVector part2)
  {
    // const TLorentzVector trackSum = part1 + part2;
    trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    // TLorentzVector PartOneCMS(part1);
    // TLorentzVector PartTwoCMS(part2);
    PartOneCMS.SetXYZM(part1.Px(), part1.Py(), part1.Pz(), part1.M());
    PartTwoCMS.SetXYZM(part2.Px(), part2.Py(), part2.Pz(), part2.M());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    // const TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;
    trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  float getmT(const TLorentzVector part1, const TLorentzVector part2)
  {
    trackSum = part1 + part2;
    float kT = 0.5 * trackSum.Pt();
    return std::sqrt((kT * kT) + 0.25 * (part1.M() + part2.M()) * (part1.M() + part2.M()));
  }

  float combinedTPC;
  TLorentzVector F1, Proton, F1ProtonPair, Pion, Kaon, Kshort;
  TLorentzVector F1Rot, PionRot, KaonKshortPair;
  // Process the data in same event

  int currentRunNumber = -999;
  int lastRunNumber = -999;

  void processSE(aod::RedF1PEvents::iterator const& collision, aod::F1Tracks const& f1tracks, aod::ProtonTracks const& protontracks)
  {

    // auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    // currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    currentRunNumber = collision.runNumber();
    if (currentRunNumber != lastRunNumber) {
      bz = getMagneticField(collision.timestamp());
    }
    lastRunNumber = currentRunNumber;

    for (auto f1track : f1tracks) {
      if (f1track.f1MassKaonKshort() > maxKKS0Mass) {
        continue;
      }
      F1.SetXYZM(f1track.f1Px(), f1track.f1Py(), f1track.f1Pz(), f1track.f1Mass());
      Pion.SetXYZM(f1track.f1d1Px(), f1track.f1d1Py(), f1track.f1d1Pz(), 0.139);
      Kaon.SetXYZM(f1track.f1d2Px(), f1track.f1d2Py(), f1track.f1d2Pz(), 0.493);
      Kshort.SetXYZM(f1track.f1d3Px(), f1track.f1d3Py(), f1track.f1d3Pz(), 0.497);
      KaonKshortPair = Kaon + Kshort;
      if (Pion.Pt() > maxMomentumPion || Kaon.Pt() > maxMomentumKaon) {
        continue;
      }
      if (pdepPID) {
        if (Kaon.Pt() <= 0.5 && (f1track.f1d2TPC() < -2.5 || f1track.f1d2TPC() > 2.5)) {
          continue;
        }
        if (Kaon.Pt() > 0.5 && Kaon.Pt() <= 0.7 && (f1track.f1d2TPC() < -1.5 || f1track.f1d2TPC() > 2.5)) {
          continue;
        }
        if (Kaon.Pt() > 0.7 && Kaon.Pt() <= 1.0 && (f1track.f1d2TPC() < -1.0 || f1track.f1d2TPC() > 2.5)) {
          continue;
        }
        if (Kaon.Pt() > 1.0 && (f1track.f1d2TPC() < -2.5 || f1track.f1d2TPC() > 2.5)) {
          continue;
        }
        if (Pion.Pt() < 2.0 && (f1track.f1d1TPC() < -2.5 || f1track.f1d1TPC() > 2.5)) {
          continue;
        }
        if (Pion.Pt() > 2.0 && (f1track.f1d1TPC() < -2.5 || f1track.f1d1TPC() > 2.5)) {
          continue;
        }
      }
      if (strategyPIDPion == 1 && Pion.Pt() > momentumTOFPionMin && Pion.Pt() <= momentumTOFPionMax && f1track.f1d1TOFHit() != 1) {
        continue;
      }
      if (strategyPIDKaon == 1 && Kaon.Pt() > momentumTOFKaonMin && Kaon.Pt() <= momentumTOFKaonMax && f1track.f1d2TOFHit() != 1) {
        continue;
      }
      if (strategyPIDKaon == 2 && Kaon.Pt() > momentumTOFKaonMin && Kaon.Pt() <= momentumTOFKaonMax && f1track.f1d2TPC() < -1.0 && f1track.f1d2TOFHit() != 1) {
        continue;
      }
      histos.fill(HIST("hNsigmaKaonTPC"), f1track.f1d2TPC(), Kaon.Pt());
      histos.fill(HIST("hNsigmaPionTPC"), f1track.f1d1TPC(), Pion.Pt());
      histos.fill(HIST("hNsigmaPionKaonTPC"), f1track.f1d1TPC(), f1track.f1d2TPC());
      if (typeofCombined == 0) {
        combinedTPC = TMath::Sqrt(f1track.f1d1TPC() * f1track.f1d1TPC() + f1track.f1d2TPC() * f1track.f1d2TPC());
      }
      if (typeofCombined == 1) {
        combinedTPC = (f1track.f1d1TPC() - f1track.f1d2TPC()) / (f1track.f1d1TPC() + f1track.f1d2TPC());
      }
      for (auto protontrack : protontracks) {
        Proton.SetXYZM(protontrack.protonPx(), protontrack.protonPy(), protontrack.protonPz(), 0.938);
        if (Proton.Pt() > momentumProtonMax || Proton.Pt() < momentumProtonMin) {
          continue;
        }
        if (Proton.P() < momentumTOFProton && TMath::Abs(protontrack.protonNsigmaTPC()) > 2.5) {
          continue;
        }
        if (Proton.P() >= momentumTOFProton && (protontrack.protonTOFHit() != 1 || TMath::Abs(protontrack.protonNsigmaTOF()) > 2.5)) {
          continue;
        }
        if ((f1track.f1PionIndex() == protontrack.f1ProtonIndex()) || (f1track.f1KaonIndex() == protontrack.f1ProtonIndex()) || (f1track.f1KshortPositiveIndex() == protontrack.f1ProtonIndex()) || (f1track.f1KshortNegativeIndex() == protontrack.f1ProtonIndex())) {
          continue;
        }
        auto relative_momentum = getkstar(F1, Proton);
        if (relative_momentum <= 0.5) {
          histos.fill(HIST("hNsigmaProtonTPC"), protontrack.protonNsigmaTPC(), Proton.Pt());
        }
        histos.fill(HIST("h2SameEventPtCorrelation"), relative_momentum, F1.Pt(), Proton.Pt());

        if (f1track.f1SignalStat() > 0) {
          // check charge
          float pairCharge = f1track.f1SignalStat() * protontrack.protonCharge();
          int f1Charge = f1track.f1SignalStat();
          int pionCharge = -1;
          int kaonCharge = 1;
          if (f1Charge == 2) {
            pionCharge = 1;
            kaonCharge = -1;
          }
          if (kaonCharge == protontrack.protonCharge())
            histos.fill(HIST("hPhaseSpaceProtonKaonSame"), Proton.Eta() - Kaon.Eta(), PhiAtSpecificRadiiTPC(Proton, Kaon, protontrack.protonCharge(), kaonCharge, bz, bz), relative_momentum); // Phase Space Proton kaon
          if (pionCharge == protontrack.protonCharge())
            histos.fill(HIST("hPhaseSpaceProtonPionSame"), Proton.Eta() - Pion.Eta(), PhiAtSpecificRadiiTPC(Proton, Pion, protontrack.protonCharge(), pionCharge, bz, bz), relative_momentum); // Phase Space Proton Pionsyst
          histos.fill(HIST("h2SameEventInvariantMassUnlike_mass"), relative_momentum, F1.Pt(), F1.M(), pairCharge, collision.numContrib());                                                  // F1 sign = 1 unlike, F1 sign = -1 like
          if (fillSparse) {
            histos.fill(HIST("SEMassUnlike"), F1.M(), F1.Pt(), Proton.Pt(), relative_momentum, combinedTPC, pairCharge);
          }
          if (fillRotation) {
            for (int nrotbkg = 0; nrotbkg < nRot; nrotbkg++) {
              auto anglestart = 5.0 * TMath::Pi() / 6.0;
              auto angleend = 7.0 * TMath::Pi() / 6.0;
              auto anglestep = (angleend - anglestart) / (1.0 * (9.0 - 1.0));
              auto rotangle = anglestart + nrotbkg * anglestep;
              auto rotPionPx = Pion.Px() * std::cos(rotangle) - Pion.Py() * std::sin(rotangle);
              auto rotPionPy = Pion.Px() * std::sin(rotangle) + Pion.Py() * std::cos(rotangle);
              PionRot.SetXYZM(rotPionPx, rotPionPy, Pion.Pz(), Pion.M());
              F1Rot = PionRot + KaonKshortPair;
              if (F1Rot.Pt() < 1.0) {
                continue;
              }
              auto relative_momentum_rot = getkstar(F1Rot, Proton);
              histos.fill(HIST("h2SameEventInvariantMassRot_mass"), relative_momentum_rot, F1Rot.Pt(), F1Rot.M(), pairCharge);
              if (fillSparse) {
                histos.fill(HIST("SEMassRot"), F1Rot.M(), F1Rot.Pt(), Proton.Pt(), relative_momentum_rot, combinedTPC, pairCharge);
              }
            }
          }
        }
        if (f1track.f1SignalStat() == -1) {
          histos.fill(HIST("h2SameEventInvariantMassLike_mass"), relative_momentum, F1.Pt(), F1.M(), protontrack.protonCharge(), collision.numContrib());
          if (fillSparse) {
            histos.fill(HIST("SEMassLike"), F1.M(), F1.Pt(), Proton.Pt(), relative_momentum, combinedTPC, protontrack.protonCharge());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(f1protoncorrelation, processSE, "Process same event", false);
  // Processing Event Mixing
  SliceCache cache;
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  BinningType colBinning{{CfgVtxBins, CfgMultBins}, true};
  Preslice<aod::F1Tracks> tracksPerCollisionPresliceF1 = aod::f1protondaughter::redF1PEventId;
  Preslice<aod::ProtonTracks> tracksPerCollisionPresliceP = aod::f1protondaughter::redF1PEventId;
  void processME(aod::RedF1PEvents& collisions,
                 aod::F1Tracks& f1tracks,
                 aod::ProtonTracks& protontracks)
  {
    for (auto const& [collision1, collision2] :
         selfCombinations(colBinning, nEvtMixingBkg, -1, collisions, collisions)) {
      if (collision1.index() == collision2.index())
        continue;

      // Preslices
      auto f1_c1 = f1tracks.sliceBy(tracksPerCollisionPresliceF1, collision1.globalIndex());
      auto f1_c2 = f1tracks.sliceBy(tracksPerCollisionPresliceF1, collision2.globalIndex());
      auto p_c1 = protontracks.sliceBy(tracksPerCollisionPresliceP, collision1.globalIndex());
      auto p_c2 = protontracks.sliceBy(tracksPerCollisionPresliceP, collision2.globalIndex());

      // -------------------------------
      // CASE 1: SE-F1P  (π mixed from c2, K+K0s from c1, proton from c1)
      // -------------------------------
      for (auto const& t1 : f1_c1) {
        if (t1.f1MassKaonKshort() > maxKKS0Mass)
          continue;

        Kaon.SetXYZM(t1.f1d2Px(), t1.f1d2Py(), t1.f1d2Pz(), 0.493);
        Kshort.SetXYZM(t1.f1d3Px(), t1.f1d3Py(), t1.f1d3Pz(), 0.497);
        KaonKshortPair = Kaon + Kshort;

        if (Kaon.Pt() > maxMomentumKaon)
          continue;
        if (pdepPID) {
          if (Kaon.Pt() <= 0.5 && (t1.f1d2TPC() < -2.5 || t1.f1d2TPC() > 2.5))
            continue;
          if (Kaon.Pt() > 0.5 && Kaon.Pt() <= 0.7 && (t1.f1d2TPC() < -1.5 || t1.f1d2TPC() > 2.5))
            continue;
          if (Kaon.Pt() > 0.7 && Kaon.Pt() <= 1.0 && (t1.f1d2TPC() < -1.0 || t1.f1d2TPC() > 2.5))
            continue;
          if (Kaon.Pt() > 1.0 && (t1.f1d2TPC() < -2.5 || t1.f1d2TPC() > 2.5))
            continue;
        }
        if (strategyPIDKaon == 1 &&
            Kaon.Pt() > momentumTOFKaonMin && Kaon.Pt() <= momentumTOFKaonMax &&
            t1.f1d2TOFHit() != 1)
          continue;

        for (auto const& t2 : p_c1) { // proton from c1
          Proton.SetXYZM(t2.protonPx(), t2.protonPy(), t2.protonPz(), 0.938);
          if (Proton.Pt() > momentumProtonMax || Proton.Pt() < momentumProtonMin)
            continue;
          if (Proton.P() < momentumTOFProton && TMath::Abs(t2.protonNsigmaTPC()) > 2.5)
            continue;
          if (Proton.P() >= momentumTOFProton && (t2.protonTOFHit() != 1 || TMath::Abs(t2.protonNsigmaTOF()) > 2.5))
            continue;

          for (auto const& t3 : f1_c2) { // pion source from c2
            Pion.SetXYZM(t3.f1d1Px(), t3.f1d1Py(), t3.f1d1Pz(), 0.139);
            if (Pion.Pt() > maxMomentumPion)
              continue;
            if (pdepPID) {
              if (Pion.Pt() < 2.0 && (t3.f1d1TPC() < -2.5 || t3.f1d1TPC() > 2.5))
                continue;
              if (Pion.Pt() >= 2.0 && (t3.f1d1TPC() < -2.5 || t3.f1d1TPC() > 2.5))
                continue;
            }
            if (strategyPIDPion == 1 &&
                Pion.Pt() > momentumTOFPionMin && Pion.Pt() <= momentumTOFPionMax &&
                t3.f1d1TOFHit() != 1)
              continue;

            // Fake f1: π(mixed) + (K+K0s from c1)
            F1 = Pion + KaonKshortPair;

            // keep only unlike-sign branch
            if (t1.f1SignalStat() <= 0)
              continue;

            int f1Charge = t1.f1SignalStat();
            float pairQ = f1Charge * t2.protonCharge();

            auto kstar = getkstar(F1, Proton);

            histos.fill(HIST("h2MixEventInvariantMassUnlike_mass_SEFP"),
                        kstar, F1.Pt(), F1.M(), pairQ);
          }
        }
      }

      // -------------------------------
      // CASE 2: DE-F1P  (π mixed from c2, K+K0s from c1, proton from c2)
      // -------------------------------
      for (auto const& t1 : f1_c1) {
        if (t1.f1MassKaonKshort() > maxKKS0Mass)
          continue;

        Kaon.SetXYZM(t1.f1d2Px(), t1.f1d2Py(), t1.f1d2Pz(), 0.493);
        Kshort.SetXYZM(t1.f1d3Px(), t1.f1d3Py(), t1.f1d3Pz(), 0.497);
        KaonKshortPair = Kaon + Kshort;

        if (Kaon.Pt() > maxMomentumKaon)
          continue;
        if (pdepPID) {
          if (Kaon.Pt() <= 0.5 && (t1.f1d2TPC() < -2.5 || t1.f1d2TPC() > 2.5))
            continue;
          if (Kaon.Pt() > 0.5 && Kaon.Pt() <= 0.7 && (t1.f1d2TPC() < -1.5 || t1.f1d2TPC() > 2.5))
            continue;
          if (Kaon.Pt() > 0.7 && Kaon.Pt() <= 1.0 && (t1.f1d2TPC() < -1.0 || t1.f1d2TPC() > 2.5))
            continue;
          if (Kaon.Pt() > 1.0 && (t1.f1d2TPC() < -2.5 || t1.f1d2TPC() > 2.5))
            continue;
        }
        if (strategyPIDKaon == 1 &&
            Kaon.Pt() > momentumTOFKaonMin && Kaon.Pt() <= momentumTOFKaonMax &&
            t1.f1d2TOFHit() != 1)
          continue;

        for (auto const& t2 : p_c2) { // proton from c2
          Proton.SetXYZM(t2.protonPx(), t2.protonPy(), t2.protonPz(), 0.938);
          if (Proton.Pt() > momentumProtonMax || Proton.Pt() < momentumProtonMin)
            continue;
          if (Proton.P() < momentumTOFProton && TMath::Abs(t2.protonNsigmaTPC()) > 2.5)
            continue;
          if (Proton.P() >= momentumTOFProton && (t2.protonTOFHit() != 1 || TMath::Abs(t2.protonNsigmaTOF()) > 2.5))
            continue;

          for (auto const& t3 : f1_c2) { // pion from c2
            Pion.SetXYZM(t3.f1d1Px(), t3.f1d1Py(), t3.f1d1Pz(), 0.139);
            if (Pion.Pt() > maxMomentumPion)
              continue;
            if (pdepPID) {
              if (Pion.Pt() < 2.0 && (t3.f1d1TPC() < -2.5 || t3.f1d1TPC() > 2.5))
                continue;
              if (Pion.Pt() >= 2.0 && (t3.f1d1TPC() < -2.5 || t3.f1d1TPC() > 2.5))
                continue;
            }
            if (strategyPIDPion == 1 &&
                Pion.Pt() > momentumTOFPionMin && Pion.Pt() <= momentumTOFPionMax &&
                t3.f1d1TOFHit() != 1)
              continue;

            F1 = Pion + KaonKshortPair;

            if (t1.f1SignalStat() <= 0)
              continue;

            int f1Charge = t1.f1SignalStat();
            float pairQ = f1Charge * t2.protonCharge();

            auto kstar = getkstar(F1, Proton);

            histos.fill(HIST("h2MixEventInvariantMassUnlike_mass_DEFP"),
                        kstar, F1.Pt(), F1.M(), pairQ);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(f1protoncorrelation, processME, "Process EventMixing for combinatorial background (SE-F1P & DE-F1P, minimal)", false);

  void processMEOpti(aod::RedF1PEvents& collisions, aod::F1Tracks& f1tracks, aod::ProtonTracks& protontracks)
  {
    // for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockFullIndexPolicy(colBinningFemto, nEvtMixing, -1, collisions, collisions))){
    for (auto const& [collision1, collision2] : selfCombinations(colBinning, nEvtMixing, -1, collisions, collisions)) {
      // LOGF(info, "Mixed event collisions: (%d, %d)", collision1.index(), collision2.index());
      if (collision1.index() == collision2.index()) {
        continue;
      }
      currentRunNumber = collision1.runNumber();
      if (currentRunNumber != lastRunNumber) {
        bz = getMagneticField(collision1.timestamp());
        bz2 = getMagneticField(collision2.timestamp());
      }
      lastRunNumber = currentRunNumber;
      auto groupF1 = f1tracks.sliceBy(tracksPerCollisionPresliceF1, collision1.globalIndex());
      auto groupProton = protontracks.sliceBy(tracksPerCollisionPresliceP, collision2.globalIndex());
      // auto groupF1 = f1tracks.sliceByCached(aod::f1protondaughter::redF1PEventId, collision1.globalIndex(), cache);
      // auto groupProton = protontracks.sliceByCached(aod::f1protondaughter::redF1PEventId, collision2.globalIndex(), cache);
      for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(groupF1, groupProton))) {
        if (t1.f1MassKaonKshort() > maxKKS0Mass) {
          continue;
        }
        F1.SetXYZM(t1.f1Px(), t1.f1Py(), t1.f1Pz(), t1.f1Mass());
        Pion.SetXYZM(t1.f1d1Px(), t1.f1d1Py(), t1.f1d1Pz(), 0.139);
        Kaon.SetXYZM(t1.f1d2Px(), t1.f1d2Py(), t1.f1d2Pz(), 0.493);
        Kshort.SetXYZM(t1.f1d3Px(), t1.f1d3Py(), t1.f1d3Pz(), 0.497);
        KaonKshortPair = Kaon + Kshort;
        if (Pion.Pt() > maxMomentumPion || Kaon.Pt() > maxMomentumKaon) {
          continue;
        }
        if (pdepPID) {
          if (Kaon.Pt() <= 0.5 && (t1.f1d2TPC() < -2.5 || t1.f1d2TPC() > 2.5)) {
            continue;
          }
          if (Kaon.Pt() > 0.5 && Kaon.Pt() <= 0.7 && (t1.f1d2TPC() < -1.5 || t1.f1d2TPC() > 2.5)) {
            continue;
          }
          if (Kaon.Pt() > 0.7 && Kaon.Pt() <= 1.0 && (t1.f1d2TPC() < -1.0 || t1.f1d2TPC() > 2.5)) {
            continue;
          }
          if (Kaon.Pt() > 1.0 && (t1.f1d2TPC() < -2.5 || t1.f1d2TPC() > 2.5)) {
            continue;
          }
          if (Pion.Pt() < 2.0 && (t1.f1d1TPC() < -2.5 || t1.f1d1TPC() > 2.5)) {
            continue;
          }
          if (Pion.Pt() > 2.0 && (t1.f1d1TPC() < -2.5 || t1.f1d1TPC() > 2.5)) {
            continue;
          }
        }
        if (strategyPIDPion == 1 && Pion.Pt() > momentumTOFPionMin && Pion.Pt() <= momentumTOFPionMax && t1.f1d1TOFHit() != 1) {
          continue;
        }
        if (strategyPIDKaon == 1 && Kaon.Pt() > momentumTOFKaonMin && Kaon.Pt() <= momentumTOFKaonMax && t1.f1d2TOFHit() != 1) {
          continue;
        }
        if (typeofCombined == 0) {
          combinedTPC = TMath::Sqrt(t1.f1d1TPC() * t1.f1d1TPC() + t1.f1d2TPC() * t1.f1d2TPC());
        }
        if (typeofCombined == 1) {
          combinedTPC = (t1.f1d1TPC() - t1.f1d2TPC()) / (t1.f1d1TPC() + t1.f1d2TPC());
        }
        Proton.SetXYZM(t2.protonPx(), t2.protonPy(), t2.protonPz(), 0.938);
        if (Proton.Pt() > momentumProtonMax || Proton.Pt() < momentumProtonMin) {
          continue;
        }
        if (Proton.P() < momentumTOFProton && TMath::Abs(t2.protonNsigmaTPC()) > 2.5) {
          continue;
        }
        if (Proton.P() >= momentumTOFProton && (t2.protonTOFHit() != 1 || TMath::Abs(t2.protonNsigmaTOF()) > 2.5)) {
          continue;
        }
        auto relative_momentum = getkstar(F1, Proton);
        if (t1.f1SignalStat() > 0) {
          float pairCharge = t1.f1SignalStat() * t2.protonCharge();
          int f1Charge = t1.f1SignalStat();
          int pionCharge = -1;
          int kaonCharge = 1;
          if (f1Charge == 2) {
            pionCharge = 1;
            kaonCharge = -1;
          }
          histos.fill(HIST("h2MixEventInvariantMassUnlike_mass"), relative_momentum, F1.Pt(), F1.M(), pairCharge, collision1.numContrib());                                         // F1 sign = 1 unlike, F1 sign = -1 like
          histos.fill(HIST("hPhaseSpaceProtonKaonMix"), Proton.Eta() - Kaon.Eta(), PhiAtSpecificRadiiTPC(Proton, Kaon, t2.protonCharge(), kaonCharge, bz, bz2), relative_momentum); // Phase Space Proton kaon
          histos.fill(HIST("hPhaseSpaceProtonPionMix"), Proton.Eta() - Pion.Eta(), PhiAtSpecificRadiiTPC(Proton, Pion, t2.protonCharge(), pionCharge, bz, bz2), relative_momentum); // Phase Space Proton Pion
          if (fillSparse) {
            histos.fill(HIST("MEMassUnlike"), F1.M(), F1.Pt(), Proton.Pt(), relative_momentum, combinedTPC, pairCharge);
          }

          if (fillRotation) {
            for (int nrotbkg = 0; nrotbkg < nRot; nrotbkg++) {
              auto anglestart = 5.0 * TMath::Pi() / 6.0;
              auto angleend = 7.0 * TMath::Pi() / 6.0;
              auto anglestep = (angleend - anglestart) / (1.0 * (9.0 - 1.0));
              auto rotangle = anglestart + nrotbkg * anglestep;
              auto rotPionPx = Pion.Px() * std::cos(rotangle) - Pion.Py() * std::sin(rotangle);
              auto rotPionPy = Pion.Px() * std::sin(rotangle) + Pion.Py() * std::cos(rotangle);
              PionRot.SetXYZM(rotPionPx, rotPionPy, Pion.Pz(), Pion.M());
              F1Rot = PionRot + KaonKshortPair;
              if (F1Rot.Pt() < 1.0) {
                continue;
              }
              auto relative_momentum_rot = getkstar(F1Rot, Proton);
              if (t1.f1SignalStat() > 0) {
                histos.fill(HIST("h2MixEventInvariantMassRot_mass"), relative_momentum_rot, F1Rot.Pt(), F1Rot.M(), pairCharge);
                if (fillSparse) {
                  histos.fill(HIST("MEMassRot"), F1Rot.M(), F1Rot.Pt(), Proton.Pt(), relative_momentum_rot, combinedTPC, pairCharge);
                }
              }
            }
          }
        }
        if (t1.f1SignalStat() == -1) {
          histos.fill(HIST("h2MixEventInvariantMassLike_mass"), relative_momentum, F1.Pt(), F1.M(), t2.protonCharge(), collision1.numContrib());
          if (fillSparse) {
            histos.fill(HIST("MEMassLike"), F1.M(), F1.Pt(), Proton.Pt(), relative_momentum, combinedTPC, t2.protonCharge());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(f1protoncorrelation, processMEOpti, "Process EventMixing for combinatorial background Optimal", false);

  void processSys(aod::RedF1PEvents::iterator const& collision,
                  aod::F1Tracks const& f1tracks,
                  aod::ProtonTracks const& protontracks)
  {
    const float maxMomPi = maxMomentumPion;
    const float maxMomK = maxMomentumKaon;

    const float pTofPiMin = momentumTOFPionMin;
    const float pTofPiMax = momentumTOFPionMax;
    const float pTofKMin = momentumTOFKaonMin;
    const float pTofKMax = momentumTOFKaonMax;

    const float pTofP = momentumTOFProton;
    const float pMaxP = momentumProtonMax;
    const float pMinP = momentumProtonMin;

    currentRunNumber = collision.runNumber();
    if (currentRunNumber != lastRunNumber) {
      bz = getMagneticField(collision.timestamp());
    }
    lastRunNumber = currentRunNumber;
    auto countf1 = 0;

    for (auto f1track : f1tracks) {

      if (f1track.f1MassKaonKshort() > maxKKS0Mass)
        continue;

      F1.SetXYZM(f1track.f1Px(), f1track.f1Py(), f1track.f1Pz(), f1track.f1Mass());
      Pion.SetXYZM(f1track.f1d1Px(), f1track.f1d1Py(), f1track.f1d1Pz(), 0.139);
      Kaon.SetXYZM(f1track.f1d2Px(), f1track.f1d2Py(), f1track.f1d2Pz(), 0.493);
      Kshort.SetXYZM(f1track.f1d3Px(), f1track.f1d3Py(), f1track.f1d3Pz(), 0.497);
      KaonKshortPair = Kaon + Kshort;

      std::vector<int> activeSys;
      activeSys.reserve((size_t)nSysTotal);

      for (int sysId = 0; sysId < nSysTotal; ++sysId) {
        const auto& sc = sysCuts[sysId];
        // Primary π
        if (!passPrimary(f1track.pionDcaxy(), f1track.pionDcaz(),
                         f1track.pionTPCNcrs(), f1track.pionTPCNcls(), sc))
          continue;
        // Primary K
        if (!passPrimary(f1track.kaonDcaxy(), f1track.kaonDcaz(),
                         f1track.kaonTPCNcrs(), f1track.kaonTPCNcls(), sc))
          continue;
        // V0 (K0s)
        if (!passV0(f1track, sc))
          continue;

        if (!passPionPID(sc.pidPi, f1track, Pion, maxMomPi))
          continue;
        if (!passKaonPID(sc.pidK, f1track, Kaon, maxMomK))
          continue;
        if (sysId == 0) {
          histos.fill(HIST("hNsigmaKaonTPC"), f1track.f1d2TPC(), Kaon.Pt());
          histos.fill(HIST("hNsigmaPionTPC"), f1track.f1d1TPC(), Pion.Pt());
          histos.fill(HIST("hNsigmaPionKaonTPC"), f1track.f1d1TPC(), f1track.f1d2TPC());
          countf1 = countf1 + 1;
        }
        activeSys.push_back(sysId);
      }

      if (activeSys.empty())
        continue;

      // Proton loop
      for (auto protontrack : protontracks) {
        Proton.SetXYZM(protontrack.protonPx(), protontrack.protonPy(), protontrack.protonPz(), 0.938);

        if ((f1track.f1PionIndex() == protontrack.f1ProtonIndex()) ||
            (f1track.f1KaonIndex() == protontrack.f1ProtonIndex()) ||
            (f1track.f1KshortPositiveIndex() == protontrack.f1ProtonIndex()) ||
            (f1track.f1KshortNegativeIndex() == protontrack.f1ProtonIndex()))
          continue;

        const auto& sc0 = sysCuts[0];

        if (countf1 && passPrimary(protontrack.protonDcaxy(), protontrack.protonDcaz(), protontrack.protonTPCNcrs(), protontrack.protonTPCNcls(), sc0)) {
          histos.fill(HIST("hNsigmaProtonTPC"), protontrack.protonNsigmaTPC(), protontrack.protonNsigmaTOF(), Proton.Pt());
        }

        // physics variables
        auto relative_momentum = getkstar(F1, Proton);
        auto mT = getmT(F1, Proton);

        std::vector<int> activePair;
        activePair.reserve(activeSys.size());

        for (int sysId : activeSys) {
          const auto& sc = sysCuts[sysId];

          if (!passPrimary(protontrack.protonDcaxy(), protontrack.protonDcaz(),
                           protontrack.protonTPCNcrs(), protontrack.protonTPCNcls(), sc))
            continue;

          if (!passProtonPID(sc.pidP, protontrack, Proton, pMinP, pMaxP, pTofP))
            continue;
          if (sysId == 0) {
            int f1Charge = f1track.f1SignalStat();
            int pionCharge = -1;
            int kaonCharge = 1;
            if (f1Charge == 2) {
              pionCharge = 1;
              kaonCharge = -1;
            }
            if (kaonCharge == protontrack.protonCharge())
              histos.fill(HIST("hPhaseSpaceProtonKaonSame"), Proton.Eta() - Kaon.Eta(), PhiAtSpecificRadiiTPC(Proton, Kaon, protontrack.protonCharge(), kaonCharge, bz, bz), relative_momentum); // Phase Space Proton kaon
            if (pionCharge == protontrack.protonCharge())
              histos.fill(HIST("hPhaseSpaceProtonPionSame"), Proton.Eta() - Pion.Eta(), PhiAtSpecificRadiiTPC(Proton, Pion, protontrack.protonCharge(), pionCharge, bz, bz), relative_momentum); // Phase Space Proton Pion
            histos.fill(HIST("h2SameEventf1pptCorrelation"), F1.M(), relative_momentum, Proton.Pt());
          }
          activePair.push_back(sysId);
        }
        if (activePair.empty())
          continue;
        if (f1track.f1SignalStat() > 0) {
          for (int sysId : activePair) {
            histos.fill(HIST("h2SameEventInvariantMassUnlike_mass_SYS"), sysId, relative_momentum, mT, F1.M(), collision.numContrib());
          }
        }
        if (f1track.f1SignalStat() == -1) {
          for (int sysId : activePair) {
            histos.fill(HIST("h2SameEventInvariantMassLike_mass_SYS"), sysId, relative_momentum, mT, F1.M(), collision.numContrib());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(f1protoncorrelation, processSys, "Process SAME-EVENT with systematics (sysId axis, default + random)", false);

  void processMESysOpti(aod::RedF1PEvents& collisions,
                        aod::F1Tracks& f1tracks,
                        aod::ProtonTracks& protontracks)
  {
    const float maxMomPi = maxMomentumPion;
    const float maxMomK = maxMomentumKaon;
    const float pTofPiMin = momentumTOFPionMin;
    const float pTofPiMax = momentumTOFPionMax;
    const float pTofKMin = momentumTOFKaonMin;
    const float pTofKMax = momentumTOFKaonMax;
    const float pTofP = momentumTOFProton;
    const float pMaxP = momentumProtonMax;
    const float pMinP = momentumProtonMin;
    // const bool  doPdep    = pdepPID;

    for (auto const& [collision1, collision2] : selfCombinations(colBinning, nEvtMixing, -1, collisions, collisions)) {
      if (collision1.index() == collision2.index())
        continue;

      currentRunNumber = collision1.runNumber();
      if (currentRunNumber != lastRunNumber) {
        bz = getMagneticField(collision1.timestamp());
        bz2 = getMagneticField(collision2.timestamp());
      }
      lastRunNumber = currentRunNumber;

      auto groupF1 = f1tracks.sliceBy(tracksPerCollisionPresliceF1, collision1.globalIndex());
      auto groupProton = protontracks.sliceBy(tracksPerCollisionPresliceP, collision2.globalIndex());

      for (auto& [t1, t2] : soa::combinations(o2::soa::CombinationsFullIndexPolicy(groupF1, groupProton))) {

        if (t1.f1MassKaonKshort() > maxKKS0Mass)
          continue;

        // 4-vectors
        F1.SetXYZM(t1.f1Px(), t1.f1Py(), t1.f1Pz(), t1.f1Mass());
        Pion.SetXYZM(t1.f1d1Px(), t1.f1d1Py(), t1.f1d1Pz(), 0.139);
        Kaon.SetXYZM(t1.f1d2Px(), t1.f1d2Py(), t1.f1d2Pz(), 0.493);
        Kshort.SetXYZM(t1.f1d3Px(), t1.f1d3Py(), t1.f1d3Pz(), 0.497);
        KaonKshortPair = Kaon + Kshort;
        Proton.SetXYZM(t2.protonPx(), t2.protonPy(), t2.protonPz(), 0.938);
        auto relative_momentum = getkstar(F1, Proton);
        auto mT = getmT(F1, Proton);
        // sys list for this (F1, p) pair
        std::vector<int> activePair;
        activePair.reserve((size_t)nSysTotal);

        for (int sysId = 0; sysId < nSysTotal; ++sysId) {
          const auto& sc = sysCuts[sysId];

          // π/K primary
          if (!passPrimary(t1.pionDcaxy(), t1.pionDcaz(), t1.pionTPCNcrs(), t1.pionTPCNcls(), sc))
            continue;
          if (!passPrimary(t1.kaonDcaxy(), t1.kaonDcaz(), t1.kaonTPCNcrs(), t1.kaonTPCNcls(), sc))
            continue;

          // V0 cuts
          if (!passV0(t1, sc))
            continue;

          // PID (F1 side)
          if (!passPionPID(sc.pidPi, t1, Pion, maxMomPi))
            continue;
          if (!passKaonPID(sc.pidK, t1, Kaon, maxMomK))
            continue;

          // proton primary + PID
          if (!passPrimary(t2.protonDcaxy(), t2.protonDcaz(), t2.protonTPCNcrs(), t2.protonTPCNcls(), sc))
            continue;
          if (!passProtonPID(sc.pidP, t2, Proton, pMinP, pMaxP, pTofP))
            continue;
          activePair.push_back(sysId);
          if (sysId == 0) {
            int f1Charge = t1.f1SignalStat();
            int pionCharge = -1;
            int kaonCharge = 1;
            if (f1Charge == 2) {
              pionCharge = 1;
              kaonCharge = -1;
            }
            if (kaonCharge == t2.protonCharge())
              histos.fill(HIST("hPhaseSpaceProtonKaonMix"), Proton.Eta() - Kaon.Eta(), PhiAtSpecificRadiiTPC(Proton, Kaon, t2.protonCharge(), kaonCharge, bz, bz), relative_momentum); // Phase Space Proton kaon
            if (pionCharge == t2.protonCharge())
              histos.fill(HIST("hPhaseSpaceProtonPionMix"), Proton.Eta() - Pion.Eta(), PhiAtSpecificRadiiTPC(Proton, Pion, t2.protonCharge(), pionCharge, bz, bz), relative_momentum); // Phase Space Proton Pion
          }
        }
        if (activePair.empty())
          continue;
        if (t1.f1SignalStat() > 0) {
          for (int sysId : activePair) {
            histos.fill(HIST("h2MixEventInvariantMassUnlike_mass_SYS"), sysId, relative_momentum, mT, F1.M(), collision1.numContrib());
          }
        }
        if (t1.f1SignalStat() == -1) {
          for (int sysId : activePair) {
            histos.fill(HIST("h2MixEventInvariantMassLike_mass_SYS"), sysId, relative_momentum, mT, F1.M(), collision1.numContrib());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(f1protoncorrelation, processMESysOpti, "Process MIX-EVENT OPTI with systematics (sysId axis, default + random)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<f1protoncorrelation>(cfgc)}; }
