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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGUD/DataModel/UDTables.h"

#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandAnalyzer {
  bool fIsMC = true;

  Preslice<o2::aod::UDMcParticles> perMcCollision = o2::aod::udmcparticle::udMcCollisionId;

  float ft0DummyTime = 32.767f;
  float fT0CBBlower = -1.0; // ns
  float fT0CBBupper = 1.0;  // ns

  std::unordered_map<int32_t, float> pdgsMass;

  enum pdgs {
    kPdgElectron = 11,
    kPdgMuon = 13,
    kPdgTau = 15,
    kPdgPion = 211
  };

  // selection flags for processCandidate()
  enum selections {
    kSelIdealPID = 0, // MC is used: PID using PDG codes
    kSelIsNotFake,    // MC is used: check MC particle IDs (background particles are not stored => mcID < 0)
    kSelUnlikeSign,
    kSelNoFT0,
    kSelNoFV0A,
    kSelNoFDD,
    kSelPID,
    kSelDDCA,
    kSelPt,
    kNSelectors
  };

  std::vector<std::vector<std::vector<float>>> fMeansSigmas;

  Configurable<int32_t> fPrimaryPdg{"primaryPdg", 13, "Set 'primary' PDG code: e.g. 15 for ditau production"};
  Configurable<int32_t> fTargetPdg{"targetPdg", 11, "Target particle PDG for 'histSwitch' p_T distributions: e.g. electrons (11) from tau decays"};
  Configurable<int32_t> fTPCPIDSwitch{"tpcPIDSwitch", 0, "TPC PID switch: 0 -- two muons, 1 -- two pions, 2 -- two electrons, 3 -- electron + muon/pion"};
  Configurable<int32_t> fHistSwitch{"histSwitch", 0, "What information to collect: 0 -- pair mass, 1 -- p_T of target particle, 2 -- both"};

  float fMinPt = 0.;
  float fMaxPt = 1.;

  static constexpr int32_t nBinsMass = 500;
  static constexpr float minMass = 0;
  static constexpr float maxMass = 5;

  static constexpr int32_t nBinsPt = 500;
  static constexpr float minPt = 0;
  static constexpr float maxPt = 5;

  HistogramRegistry registry{
    "registry",
    {{"MC/PairMass", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"MC/Eta", ";#eta;", {HistType::kTH1D, {{100, -6., 6.}}}},
     // separate selectors stored in "SelCounter" (see init())
     {"Selection/PairMass/All", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/UnlikeSign", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/NoFT0", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/NoFV0A", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/NoFDD", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/DDCA", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/DDCA_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/PID", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/PIDSig", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/PID_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/PIDSig_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/IdealPID", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/IdealPIDSig", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/IdealPID_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/IdealPIDSig_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     //
     {"Selection/TargetPt/All", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/UnlikeSign", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/NoFT0", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/NoFV0A", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/NoFDD", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/DDCA", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/DDCA_MC", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/PID", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/PIDSig", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/PID_MC", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/PIDSig_MC", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/IdealPID", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/IdealPIDSig", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/IdealPID_MC", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/IdealPIDSig_MC", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     //
     {"Selection/TPCSignals/All", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"Selection/TPCSignals/IdealPID", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"Selection/TPCSignals/PID", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"Selection/TPCSignals/IdealPIDSig", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"Selection/TPCSignals/PIDSig", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"Selection/TPCSignals/DDCA", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"Selection/TPCSignals/DDCASig", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     //
     {"Selection/DDCA/All", ";DDCA_z; DDCA_xy;", {HistType::kTH2D, {{500, -5., 5.}, {500, -5., 5.}}}},
     {"Selection/DDCA/Sig", ";DDCA_z; DDCA_xy;", {HistType::kTH2D, {{500, -5., 5.}, {500, -5., 5.}}}},
     {"Selection/DDCA/MvsDCAZAll", ";#it{m}, GeV; DDCA_z;", {HistType::kTH2D, {{500, 0., 5.}, {500, -5., 5.}}}},
     {"Selection/DDCA/MvsDCAXYAll", ";#it{m}, GeV; DDCA_xy;", {HistType::kTH2D, {{500, 0., 5.}, {500, -5., 5.}}}},
     //
     {"Selection/TPC_nSigmaEl", ";#sigma;", {HistType::kTH1D, {{200, -10., 10.}}}},
     {"Selection/TPC_nSigmaPi", ";#sigma;", {HistType::kTH1D, {{200, -10., 10.}}}}}};

  using Candidates = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSels>;
  using BarrelTracks = soa::Join<o2::aod::UDTracks, o2::aod::UDTracksExtra, o2::aod::UDTracksDCA, o2::aod::UDTracksPID>;
  using FwdTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;

  void init(InitContext&)
  {
    const AxisSpec axisSel{kNSelectors, 0., double(kNSelectors), ""};
    registry.add("Selection/SelCounter", "", kTH1F, {axisSel});
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelIdealPID + 1, "kSelIdealPID");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelIsNotFake + 1, "kSelIsNotFake");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelUnlikeSign + 1, "kSelUnlikeSign");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelNoFT0 + 1, "kSelNoFT0");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelNoFV0A + 1, "kSelNoFV0A");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelNoFDD + 1, "kSelNoFDD");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelPID + 1, "kSelPID");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelDDCA + 1, "kSelDDCA");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelPt + 1, "kSelPt");

    // populate "pdg->particle mass" map
    pdgsMass[kPdgElectron] = 0.000511;
    pdgsMass[kPdgMuon] = 0.10566;
    pdgsMass[kPdgTau] = 1.777;
    pdgsMass[kPdgPion] = 0.13957;
  }

  void processMCParts(o2::aod::UDMcCollisions const& mcCollisions, o2::aod::UDMcParticles const& mcParticles)
  {
    int32_t nMCEvents = mcCollisions.size();
    // collect MC distributions
    for (int32_t i = 0; i < nMCEvents; i++) {
      auto mcPartsGroup = mcParticles.sliceBy(perMcCollision, i);
      std::vector<aod::UDMcParticle> mcPartsFiltered;
      for (const auto& mcPart : mcPartsGroup) {
        if (std::abs(mcPart.pdgCode()) != fPrimaryPdg) {
          continue;
        }
        mcPartsFiltered.emplace_back(mcPart);
      }
      // sanity check
      if (mcPartsFiltered.size() != 2) {
        continue;
      }
      const auto& part1 = mcPartsFiltered[0];
      const auto& part2 = mcPartsFiltered[1];
      TLorentzVector p1, p2, p;
      p1.SetXYZM(part1.px(), part1.py(), part1.pz(), pdgsMass[fPrimaryPdg]);
      p2.SetXYZM(part2.px(), part2.py(), part2.pz(), pdgsMass[fPrimaryPdg]);
      p = p1 + p2;
      registry.fill(HIST("MC/PairMass"), p.M());
      registry.fill(HIST("MC/Eta"), p1.Eta());
      registry.fill(HIST("MC/Eta"), p2.Eta());
      mcPartsFiltered.clear();
    }
  }

  template <typename TTrack>
  bool checkTPCPID(TTrack& tr1, TTrack& tr2, float& m1, float& m2, bool* pidFlags, int32_t /*pdg1*/ = -1, int32_t /*pdg2*/ = -1)
  {
    bool pass = false;

    float nTPCSigmaEl1 = tr1.tpcNSigmaEl();
    float nTPCSigmaPi1 = tr1.tpcNSigmaPi();
    float nTPCSigmaEl2 = tr2.tpcNSigmaEl();
    float nTPCSigmaPi2 = tr2.tpcNSigmaPi();

    // "soft" selection
    bool isEl1 = std::abs(nTPCSigmaEl1) < 3.0f && std::abs(nTPCSigmaEl1) < std::abs(nTPCSigmaPi1);
    bool isEl2 = std::abs(nTPCSigmaEl2) < 3.0f && std::abs(nTPCSigmaEl2) < std::abs(nTPCSigmaPi2);
    bool isPi1 = tr1.tpcSignal() < 77; // std::abs(nTPCSigmaPi1) < 3.0f && std::abs(nTPCSigmaEl1) > std::abs(nTPCSigmaPi1);
    bool isPi2 = tr2.tpcSignal() < 77; // std::abs(nTPCSigmaPi2) < 3.0f && std::abs(nTPCSigmaEl2) > std::abs(nTPCSigmaPi2);

    // PID flags to check sigmas later
    pidFlags[0] = isEl1;
    pidFlags[1] = isEl2;
    pidFlags[2] = isPi1;
    pidFlags[3] = isPi2;

    // QA histograms
    if (isEl1)
      registry.fill(HIST("Selection/TPC_nSigmaEl"), nTPCSigmaEl1);
    if (isEl2)
      registry.fill(HIST("Selection/TPC_nSigmaEl"), nTPCSigmaEl2);
    if (isPi1)
      registry.fill(HIST("Selection/TPC_nSigmaPi"), nTPCSigmaPi1);
    if (isPi2)
      registry.fill(HIST("Selection/TPC_nSigmaPi"), nTPCSigmaPi2);

    // two muons/pions
    if ((fTPCPIDSwitch == 0 || fTPCPIDSwitch == 1) && isPi1 && isPi2) {
      pass = true;
      m1 = pdgsMass[kPdgPion];
      m2 = pdgsMass[kPdgPion];
    }

    // two electrons
    if (fTPCPIDSwitch == 2 && isEl1 && isEl2) {
      pass = true;
      m1 = pdgsMass[kPdgElectron];
      m2 = pdgsMass[kPdgElectron];
    }

    // electron + muon/pion
    if (fTPCPIDSwitch == 3) {
      if (isEl1 && isPi2 && !isEl2 && !isPi1) {
        pass = true;
        m1 = pdgsMass[kPdgElectron];
        m2 = pdgsMass[kPdgPion];
      }
      if (!isEl1 && !isPi2 && isEl2 && isPi1) {
        pass = true;
        m1 = pdgsMass[kPdgPion];
        m2 = pdgsMass[kPdgElectron];
      }
    }

    return pass;
  }

  bool checkIdealPID(int32_t pdg1, int32_t pdg2)
  {
    bool pass = false;
    // two muons
    if (fTPCPIDSwitch == 0) {
      bool trueSig1 = std::abs(pdg1) == kPdgMuon;
      bool trueSig2 = std::abs(pdg2) == kPdgMuon;
      pass = trueSig1 && trueSig2;
    }
    // two pions
    if (fTPCPIDSwitch == 1) {
      bool trueSig1 = std::abs(pdg1) == kPdgPion;
      bool trueSig2 = std::abs(pdg2) == kPdgPion;
      pass = trueSig1 && trueSig2;
    }
    // two electrons
    if (fTPCPIDSwitch == 2) {
      bool trueSig1 = std::abs(pdg1) == kPdgElectron;
      bool trueSig2 = std::abs(pdg2) == kPdgElectron;
      pass = trueSig1 && trueSig2;
    }
    // electron + muon/pion
    if (fTPCPIDSwitch == 3) {
      bool trueSig1 = std::abs(pdg1) == kPdgElectron && (std::abs(pdg2) == kPdgMuon || std::abs(pdg2) == kPdgPion);
      bool trueSig2 = std::abs(pdg2) == kPdgElectron && (std::abs(pdg1) == kPdgMuon || std::abs(pdg1) == kPdgPion);
      pass = trueSig1 || trueSig2;
    }
    return pass;
  }

  void fillMassDistr(float m, float mmc, std::unordered_map<int32_t, bool>& selFlags)
  {
    if (mmc < 0) { // just fill reco mass, if not using MC
      mmc = m;
    }
    registry.fill(HIST("Selection/PairMass/All"), m);
    // unlike-sign
    bool selector = selFlags[kSelPt] && selFlags[kSelUnlikeSign];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/UnlikeSign"), m);
    }
    // unlike sign + no FT0
    selector = selector && selFlags[kSelNoFT0];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/NoFT0"), m);
    }
    selector = selector && selFlags[kSelNoFV0A];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/NoFV0A"), m);
    }
    selector = selector && selFlags[kSelNoFDD];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/NoFDD"), m);
    }
    // unlike sign + no FT0 + delta(DCA)
    bool selectorDCA = selector && selFlags[kSelDDCA];
    if (selectorDCA) {
      registry.fill(HIST("Selection/PairMass/DDCA"), m);
      registry.fill(HIST("Selection/PairMass/DDCA_MC"), mmc);
    }
    // unlike sign + no FT0 + [ideal]PID
    bool selectorIdeal = selector && selFlags[kSelIdealPID];
    selector = selector && selFlags[kSelPID];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/PID"), m);
      registry.fill(HIST("Selection/PairMass/PID_MC"), mmc);
    }
    if (selectorIdeal) {
      registry.fill(HIST("Selection/PairMass/IdealPID"), m);
      registry.fill(HIST("Selection/PairMass/IdealPID_MC"), mmc);
    }
    if (selector && selFlags[kSelIsNotFake]) {
      registry.fill(HIST("Selection/PairMass/PIDSig"), m);
      registry.fill(HIST("Selection/PairMass/PIDSig_MC"), mmc);
    }
    if (selectorIdeal && selFlags[kSelIsNotFake]) {
      registry.fill(HIST("Selection/PairMass/IdealPIDSig"), m);
      registry.fill(HIST("Selection/PairMass/IdealPIDSig_MC"), mmc);
    }
  }

  void fillPtDistr(float pt, float ptmc, std::unordered_map<int32_t, bool>& selFlags)
  {
    if (ptmc < 0) { // just fill reco pt, if not using MC
      ptmc = pt;
    }
    registry.fill(HIST("Selection/TargetPt/All"), pt);
    // unlike-sign
    bool selector = selFlags[kSelPt] && selFlags[kSelUnlikeSign];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/UnlikeSign"), pt);
    }
    // unlike sign + no FT0
    selector = selector && selFlags[kSelNoFT0];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/NoFT0"), pt);
    }
    selector = selector && selFlags[kSelNoFV0A];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/NoFV0A"), pt);
    }
    selector = selector && selFlags[kSelNoFDD];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/NoFDD"), pt);
    }
    // unlike sign + no FT0 + delta(DCA)
    bool selectorDCA = selector && selFlags[kSelDDCA];
    if (selectorDCA) {
      registry.fill(HIST("Selection/TargetPt/DDCA"), pt);
      registry.fill(HIST("Selection/TargetPt/DDCA_MC"), ptmc);
    }
    // unlike sign + no FT0 + [ideal]PID
    bool selectorIdeal = selector && selFlags[kSelIdealPID];
    selector = selector && selFlags[kSelPID];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/PID"), pt);
      registry.fill(HIST("Selection/TargetPt/PID_MC"), ptmc);
    }
    if (selectorIdeal) {
      registry.fill(HIST("Selection/TargetPt/IdealPID"), pt);
      registry.fill(HIST("Selection/TargetPt/IdealPID_MC"), ptmc);
    }
    if (selector && selFlags[kSelIsNotFake]) {
      registry.fill(HIST("Selection/TargetPt/PIDSig"), pt);
      registry.fill(HIST("Selection/TargetPt/PIDSig_MC"), ptmc);
    }
    if (selectorIdeal && selFlags[kSelIsNotFake]) {
      registry.fill(HIST("Selection/TargetPt/IdealPIDSig"), pt);
      registry.fill(HIST("Selection/TargetPt/IdealPIDSig_MC"), ptmc);
    }
  }

  // naive DCA check
  template <typename TTrack>
  bool checkDDCA(TTrack& tr1, TTrack& tr2, TLorentzVector& p, bool isNotFake = false)
  {
    float dDcaZ = tr1.dcaZ() - tr2.dcaZ();
    float dDcaXY = tr1.dcaXY() - tr2.dcaXY();
    registry.fill(HIST("Selection/DDCA/All"), dDcaZ, dDcaXY);
    registry.fill(HIST("Selection/DDCA/MvsDCAZAll"), p.M(), dDcaZ);
    registry.fill(HIST("Selection/DDCA/MvsDCAXYAll"), p.M(), dDcaXY);
    if (isNotFake) {
      registry.fill(HIST("Selection/DDCA/Sig"), dDcaZ, dDcaXY);
    }
    bool pass = std::abs(dDcaZ) < 4.f;
    return pass;
  }

  template <int32_t processSwitch, typename TTrack1, typename TTrack2>
  void processCandidate(Candidates::iterator const& cand, TTrack1& tr1, TTrack2& tr2,
                        o2::aod::UDMcParticles* mcParticles, o2::aod::UDMcTrackLabels* mcTrackLabels, o2::aod::UDMcFwdTrackLabels* mcFwdTrackLabels)
  {
    std::unordered_map<int32_t, bool> selFlags; // holder of selection flags
    float mmc = -1;
    float pt1mc = -1;
    float pt2mc = -1;
    int32_t pdg1 = -1;
    int32_t pdg2 = -1;
    //
    if (fIsMC) {
      TLorentzVector mcP1, mcP2;
      int32_t mcPartId1;
      int32_t mcPartId2;
      // forward
      if constexpr (processSwitch == 0) {
        mcPartId1 = mcFwdTrackLabels->iteratorAt(tr1.globalIndex()).udMcParticleId();
        mcPartId2 = mcFwdTrackLabels->iteratorAt(tr2.globalIndex()).udMcParticleId();
      }
      // semi-forward
      if constexpr (processSwitch == 1) {
        mcPartId1 = mcFwdTrackLabels->iteratorAt(tr1.globalIndex()).udMcParticleId();
        mcPartId2 = mcTrackLabels->iteratorAt(tr2.globalIndex()).udMcParticleId();
      }
      // central
      if constexpr (processSwitch == 2) {
        mcPartId1 = mcTrackLabels->iteratorAt(tr1.globalIndex()).udMcParticleId();
        mcPartId2 = mcTrackLabels->iteratorAt(tr2.globalIndex()).udMcParticleId();
      }
      const auto& mcPart1 = mcParticles->iteratorAt(mcPartId1);
      const auto& mcPart2 = mcParticles->iteratorAt(mcPartId2);
      pdg1 = mcPart1.pdgCode();
      pdg2 = mcPart2.pdgCode();
      float m1mc = pdgsMass[pdg1];
      float m2mc = pdgsMass[pdg2];
      mcP1.SetXYZM(mcPart1.px(), mcPart1.py(), mcPart1.pz(), m1mc);
      mcP2.SetXYZM(mcPart2.px(), mcPart2.py(), mcPart2.pz(), m2mc);
      pt1mc = mcP1.Pt();
      pt2mc = mcP2.Pt();
      mmc = (mcP1 + mcP2).M();
      selFlags[kSelIdealPID] = checkIdealPID(pdg1, pdg2);
      selFlags[kSelIsNotFake] = mcPartId1 != -1 && mcPartId2 != -1;
    }
    // checking "realistic" PID and DCA for central barrel
    float m1 = pdgsMass[fTargetPdg];
    float m2 = pdgsMass[fTargetPdg];
    bool pidFlags[4] = {false};
    if constexpr (processSwitch == 2) {
      selFlags[kSelPID] = checkTPCPID(tr1, tr2, m1, m2, pidFlags, pdg1, pdg2);
    }
    // unlike-sign tracks requirement
    selFlags[kSelUnlikeSign] = (tr1.sign() * tr2.sign()) < 0;
    // check FT0 signal
    bool hasNoFT0 = true;
    if constexpr (processSwitch == 2) {
      bool isBB = cand.bbFT0A() || cand.bbFT0C();
      bool isBG = cand.bgFT0A() || cand.bgFT0C();
      hasNoFT0 = !isBB && !isBG;
    } else {
      // if there is a signal, candidate passes if timeA is dummy
      // and timeC is between +/- 1 ns
      bool checkA = std::abs(cand.timeFT0A() - ft0DummyTime) < 1e-3;
      bool checkC = cand.timeFT0C() > fT0CBBlower && cand.timeFT0C() < fT0CBBupper;
      hasNoFT0 = checkA && checkC;
    }
    selFlags[kSelNoFT0] = hasNoFT0;
    // check FV0 signal
    bool hasNoFV0A = true;
    if constexpr (processSwitch == 2) {
      bool isBB = cand.bbFV0A();
      bool isBG = cand.bgFV0A();
      hasNoFV0A = !isBB && !isBG;
    }
    selFlags[kSelNoFV0A] = hasNoFV0A;
    // check FDD signal
    bool hasNoFDD = true;
    if constexpr (processSwitch == 2) {
      bool isBB = cand.bbFDDA() || cand.bbFDDC();
      bool isBG = cand.bgFDDA() || cand.bgFDDC();
      hasNoFDD = !isBB && !isBG;
    }
    selFlags[kSelNoFDD] = hasNoFDD;
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m1);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m2);
    TLorentzVector p = p1 + p2;
    selFlags[kSelPt] = p.Pt() > fMinPt && p.Pt() < fMaxPt;
    if constexpr (processSwitch == 2) {
      selFlags[kSelDDCA] = checkDDCA(tr1, tr2, p, selFlags[kSelIsNotFake]);
    }
    // selection counters
    if (selFlags[kSelIdealPID]) { // has real meaning for MC only
      registry.fill(HIST("Selection/SelCounter"), kSelIdealPID, 1);
    }
    if (selFlags[kSelIsNotFake]) { // has real meaning for MC only
      registry.fill(HIST("Selection/SelCounter"), kSelIsNotFake, 1);
    }
    if (selFlags[kSelUnlikeSign]) {
      registry.fill(HIST("Selection/SelCounter"), kSelUnlikeSign, 1);
    }
    if (selFlags[kSelNoFT0]) {
      registry.fill(HIST("Selection/SelCounter"), kSelNoFT0, 1);
    }
    if (selFlags[kSelNoFV0A]) {
      registry.fill(HIST("Selection/SelCounter"), kSelNoFV0A, 1);
    }
    if (selFlags[kSelNoFDD]) {
      registry.fill(HIST("Selection/SelCounter"), kSelNoFDD, 1);
    }
    if (selFlags[kSelPID]) {
      registry.fill(HIST("Selection/SelCounter"), kSelPID, 1);
    }
    if (selFlags[kSelDDCA]) {
      registry.fill(HIST("Selection/SelCounter"), kSelDDCA, 1);
    }
    if (selFlags[kSelPt]) {
      registry.fill(HIST("Selection/SelCounter"), kSelPt, 1);
    }
    // collect mass distributions if needed
    if (fHistSwitch == 0 || fHistSwitch == 2) {
      float m = p.M();
      fillMassDistr(m, mmc, selFlags);
    }
    // collect pt distributions if needed
    if (fHistSwitch == 1 || fHistSwitch == 2) {
      float pt1 = p1.Pt();
      float pt2 = p2.Pt();
      bool fill1 = false;
      bool fill2 = false;
      if (fIsMC) { // if this is MC, we can just use real PDG
        fill1 = pdg1 == fTargetPdg;
        fill2 = pdg2 == fTargetPdg;
      } else { // otherwise, using PID hypotheses based on TPC signals
        if (fTargetPdg == kPdgElectron) {
          fill1 = pidFlags[0];
          fill2 = pidFlags[1];
        }
        if (fTargetPdg == kPdgMuon || fTargetPdg == kPdgPion) {
          fill1 = pidFlags[2];
          fill2 = pidFlags[3];
        }
      }
      if (fill1) {
        fillPtDistr(pt1, pt1mc, selFlags);
      }
      if (fill2) {
        fillPtDistr(pt2, pt2mc, selFlags);
      }
    }
    // collect TPC signals if Central Barrel
    if constexpr (processSwitch == 2) {
      registry.fill(HIST("Selection/TPCSignals/All"), tr1.tpcSignal(), tr2.tpcSignal());
      if (selFlags[kSelIdealPID]) {
        registry.fill(HIST("Selection/TPCSignals/IdealPID"), tr1.tpcSignal(), tr2.tpcSignal());
      }
      if (selFlags[kSelIdealPID] && selFlags[kSelIsNotFake]) { // for MC
        registry.fill(HIST("Selection/TPCSignals/IdealPIDSig"), tr1.tpcSignal(), tr2.tpcSignal());
      }
      if (selFlags[kSelPID]) {
        registry.fill(HIST("Selection/TPCSignals/PID"), tr1.tpcSignal(), tr2.tpcSignal());
      }
      if (selFlags[kSelPID] && selFlags[kSelIsNotFake]) { // for MC
        registry.fill(HIST("Selection/TPCSignals/PIDSig"), tr1.tpcSignal(), tr2.tpcSignal());
      }
      if (selFlags[kSelDDCA]) {
        registry.fill(HIST("Selection/TPCSignals/DDCA"), tr1.tpcSignal(), tr2.tpcSignal());
      }
      if (selFlags[kSelDDCA] && selFlags[kSelIsNotFake]) { // for MC
        registry.fill(HIST("Selection/TPCSignals/DDCASig"), tr1.tpcSignal(), tr2.tpcSignal());
      }
    }
  }

  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  // process candidates with 2 muon tracks
  void processFwdMC(Candidates const& eventCandidates,
                    FwdTracks const& fwdTracks,
                    o2::aod::UDMcCollisions const& mcCollisions,
                    o2::aod::UDMcParticles& mcParticles,
                    o2::aod::UDMcFwdTrackLabels& mcFwdTrackLabels)
  {
    fIsMC = true;

    processMCParts(mcCollisions, mcParticles);

    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    // assuming that candidates have exatly 2 muon tracks and 0 barrel tracks
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      const auto& cand = eventCandidates.iteratorAt(candID);
      const auto& tr1 = fwdTracks.iteratorAt(trId1);
      const auto& tr2 = fwdTracks.iteratorAt(trId2);
      processCandidate<0>(cand, tr1, tr2, &mcParticles, (o2::aod::UDMcTrackLabels*)nullptr, &mcFwdTrackLabels);
    }
  }

  // process candidates with 2 muon tracks
  void processFwd(Candidates const& eventCandidates,
                  FwdTracks const& fwdTracks)
  {
    fIsMC = false;

    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    // assuming that candidates have exatly 2 muon tracks and 0 barrel tracks
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      const auto& cand = eventCandidates.iteratorAt(candID);
      const auto& tr1 = fwdTracks.iteratorAt(trId1);
      const auto& tr2 = fwdTracks.iteratorAt(trId2);
      processCandidate<0>(cand, tr1, tr2, (o2::aod::UDMcParticles*)nullptr, (o2::aod::UDMcTrackLabels*)nullptr, (o2::aod::UDMcFwdTrackLabels*)nullptr);
    }
  }

  // process candidates with 1 muon and 1 barrel tracks
  void processSemiFwdMC(Candidates const& eventCandidates,
                        FwdTracks const& fwdTracks,
                        BarrelTracks const& barTracks,
                        o2::aod::UDMcCollisions const& mcCollisions,
                        o2::aod::UDMcParticles& mcParticles,
                        o2::aod::UDMcFwdTrackLabels& mcFwdTrackLabels,
                        o2::aod::UDMcTrackLabels& mcTrackLabels)
  {
    fIsMC = true;

    processMCParts(mcCollisions, mcParticles);

    // "value" = vectors of track IDs
    //  first track -> forward
    //  second track -> central barrel
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);
    collectCandIDs(tracksPerCand, barTracks);

    // assuming that candidates have exatly 1 muon track and 1 barrel track
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      const auto& cand = eventCandidates.iteratorAt(candID);
      const auto& tr1 = fwdTracks.iteratorAt(trId1);
      const auto& tr2 = barTracks.iteratorAt(trId2);
      processCandidate<1>(cand, tr1, tr2, &mcParticles, &mcTrackLabels, &mcFwdTrackLabels);
    }
  }

  // process candidates with 1 muon and 1 barrel tracks
  void processSemiFwd(Candidates const& eventCandidates,
                      FwdTracks const& fwdTracks,
                      BarrelTracks const& barTracks)
  {
    fIsMC = false;

    // "value" = vectors of track IDs
    //  first track -> forward
    //  second track -> central barrel
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);
    collectCandIDs(tracksPerCand, barTracks);

    // assuming that candidates have exatly 1 muon track and 1 barrel track
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      const auto& cand = eventCandidates.iteratorAt(candID);
      const auto& tr1 = fwdTracks.iteratorAt(trId1);
      const auto& tr2 = barTracks.iteratorAt(trId2);
      processCandidate<1>(cand, tr1, tr2, (o2::aod::UDMcParticles*)nullptr, (o2::aod::UDMcTrackLabels*)nullptr, (o2::aod::UDMcFwdTrackLabels*)nullptr);
    }
  }

  // process candidates with 2 central barrel tracks
  void processCentralMC(Candidates const& eventCandidates,
                        BarrelTracks const& barTracks,
                        o2::aod::UDMcCollisions const& mcCollisions,
                        o2::aod::UDMcParticles& mcParticles,
                        o2::aod::UDMcTrackLabels& mcTrackLabels)
  {
    fIsMC = true;

    processMCParts(mcCollisions, mcParticles);

    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, barTracks);

    // assuming that candidates have exatly 2 central barrel tracks
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      const auto& cand = eventCandidates.iteratorAt(candID);
      const auto& tr1 = barTracks.iteratorAt(trId1);
      const auto& tr2 = barTracks.iteratorAt(trId2);
      processCandidate<2>(cand, tr1, tr2, &mcParticles, &mcTrackLabels, (o2::aod::UDMcFwdTrackLabels*)nullptr);
    }
  }

  // process candidates with 2 central barrel tracks
  void processCentral(Candidates const& eventCandidates,
                      BarrelTracks const& barTracks)
  {
    fIsMC = false;

    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, barTracks);

    // assuming that candidates have exatly 2 central barrel tracks
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      const auto& cand = eventCandidates.iteratorAt(candID);
      const auto& tr1 = barTracks.iteratorAt(trId1);
      const auto& tr2 = barTracks.iteratorAt(trId2);
      processCandidate<2>(cand, tr1, tr2, (o2::aod::UDMcParticles*)nullptr, (o2::aod::UDMcTrackLabels*)nullptr, (o2::aod::UDMcFwdTrackLabels*)nullptr);
    }
  }

  PROCESS_SWITCH(UpcCandAnalyzer, processFwdMC, "Analyse forward candidates with MC information", false);
  PROCESS_SWITCH(UpcCandAnalyzer, processSemiFwdMC, "Analyse semiforward candidates with MC information", false);
  PROCESS_SWITCH(UpcCandAnalyzer, processCentralMC, "Analyse central candidates with MC information", false);
  PROCESS_SWITCH(UpcCandAnalyzer, processFwd, "Analyse forward candidates", false);
  PROCESS_SWITCH(UpcCandAnalyzer, processSemiFwd, "Analyse semiforward candidates", false);
  PROCESS_SWITCH(UpcCandAnalyzer, processCentral, "Analyse central candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcCandAnalyzer>(cfgc)};
}
