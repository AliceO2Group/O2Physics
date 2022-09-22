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
  std::map<int32_t, float> pdgsMass;

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
    kSelPID,
    kSelDDCA,
    kNSelectors
  };

  std::vector<std::vector<std::vector<float>>> fMeansSigmas;

  Configurable<int32_t> fPrimaryPdg{"primaryPdg", 13, "Set 'primary' PDG code: e.g. 15 for ditau production"};
  Configurable<int32_t> fTargetPdg{"targetPdg", 11, "Target particle PDG for 'histSwitch' p_T distributions: e.g. electrons (11) from tau decays"};
  Configurable<int32_t> fTPCPIDSwitch{"tpcPIDSwitch", 0, "TPC PID switch: 0 -- two muons, 1 -- two pions, 2 -- two electrons, 3 -- electron + muon/pion"};
  Configurable<int32_t> fHistSwitch{"histSwitch", 0, "What information to collect: 0 -- pair mass, 1 -- p_T of target particle, 2 -- both"};

  HistogramRegistry registry{
    "registry",
    {{"CollectMC/PairMass", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"CollectMC/Eta", ";#eta;", {HistType::kTH1D, {{100, -6., 6.}}}},
     // separate selectors stored in "SelCounter" (see init())
     {"ProcessCandidate/PairMass/All", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/UnlikeSign", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/NoFT0", ";#eta;", {HistType::kTH1D, {{100, -6., 6.}}}},
     {"ProcessCandidate/PairMass/PID", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/DDCA", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/DDCA_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/PIDSig", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/PID_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/PIDSig_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/IdealPID", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/IdealPIDSig", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/IdealPID_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/PairMass/IdealPIDSig_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"ProcessCandidate/TargetPt/All", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/TargetPt/UnlikeSign", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/TargetPt/NoFT0", ";#eta;", {HistType::kTH1D, {{100, -6., 6.}}}},
     {"ProcessCandidate/TargetPt/PID", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/TargetPt/IdealPID", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/TargetPt/PID_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/TargetPt/IdealPID_MC", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"ProcessCandidate/TargetPt/IsNotFake", ";#it{m}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"ProcessCandidate/TPCSignals/All", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"ProcessCandidate/TPCSignals/IdealPID", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"ProcessCandidate/TPCSignals/PID", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"ProcessCandidate/TPCSignals/IdealPIDSig", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     {"ProcessCandidate/TPCSignals/PIDSig", ";TPC signal 1; TPC signal 2;", {HistType::kTH2D, {{200, 0., 200.}, {200, 0., 200.}}}},
     //
     {"PID_TPC/nSigmaEl", ";#sigma;", {HistType::kTH1D, {{200, -10., 10.}}}},
     {"PID_TPC/nSigmaPi", ";#sigma;", {HistType::kTH1D, {{200, -10., 10.}}}}}};

  using BarrelTracks = soa::Join<o2::aod::UDTracks, o2::aod::UDTracksCov, o2::aod::UDTracksExtra, o2::aod::UDTracksDCA, o2::aod::UDTracksPID, o2::aod::UDMcTrackLabels, o2::aod::UDTrackCollisionIDs>;
  using FwdTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTrackCollisionIDs, o2::aod::UDFwdTracksExtra, o2::aod::UDMcFwdTrackLabels>;

  void init(InitContext&)
  {
    const AxisSpec axisSel{kNSelectors, 0., double(kNSelectors), ""};
    registry.add("ProcessCandidate/SelCounter", "", kTH1F, {axisSel});
    registry.get<TH1>(HIST("ProcessCandidate/SelCounter"))->GetXaxis()->SetBinLabel(kSelIdealPID + 1, "kSelIdealPID");
    registry.get<TH1>(HIST("ProcessCandidate/SelCounter"))->GetXaxis()->SetBinLabel(kSelIsNotFake + 1, "kSelIsNotFake");
    registry.get<TH1>(HIST("ProcessCandidate/SelCounter"))->GetXaxis()->SetBinLabel(kSelUnlikeSign + 1, "kSelUnlikeSign");
    registry.get<TH1>(HIST("ProcessCandidate/SelCounter"))->GetXaxis()->SetBinLabel(kSelNoFT0 + 1, "kSelNoFT0");
    registry.get<TH1>(HIST("ProcessCandidate/SelCounter"))->GetXaxis()->SetBinLabel(kSelPID + 1, "kSelPID");
    registry.get<TH1>(HIST("ProcessCandidate/SelCounter"))->GetXaxis()->SetBinLabel(kSelDDCA + 1, "kSelDDCA");

    // populate "pdg->particle mass" map
    pdgsMass[11] = 0.000511;
    pdgsMass[13] = 0.10566;
    pdgsMass[15] = 1.777;
    pdgsMass[211] = 0.13957;

    // parameters from manual fitting based on full simulations
    // todo: make automatic ?
    fMeansSigmas.resize(5);
    // electrons
    fMeansSigmas[0] = {{70.523, 10.193},
                       {71.996, 10.088},
                       {72.648, 10.213},
                       {72.912, 10.210},
                       {72.894, 10.389},
                       {73.046, 10.516},
                       {73.234, 10.491},
                       {73.609, 10.624},
                       {73.705, 10.819},
                       {74.243, 10.517},
                       {74.209, 10.837}};
    // pions
    fMeansSigmas[1] = {{40.108, 6.045},
                       {39.629, 5.966},
                       {39.694, 6.052},
                       {40.024, 6.133},
                       {40.462, 6.242},
                       {40.967, 6.291},
                       {41.474, 6.460},
                       {43.257, 6.839},
                       {47.087, 7.478},
                       {49.761, 7.884},
                       {52.781, 8.626}};
    // muons
    fMeansSigmas[2] = {{42.628, 6.774},
                       {43.057, 6.881},
                       {43.608, 7.055},
                       {44.398, 7.109},
                       {45.040, 7.290},
                       {45.637, 7.410},
                       {46.207, 7.483},
                       {48.678, 8.132},
                       {52.879, 8.540},
                       {55.512, 9.224},
                       {58.778, 9.619}};
    // kaons
    fMeansSigmas[3] = {{94.405, 18.034},
                       {75.569, 14.678},
                       {61.561, 10.732},
                       {53.871, 8.798},
                       {49.074, 7.879},
                       {45.916, 7.310},
                       {43.852, 6.934},
                       {40.836, 6.508},
                       {40.044, 6.371},
                       {41.214, 6.563},
                       {43.406, 7.178}};
    // protons
    fMeansSigmas[4] = {{232.635, 25.993},
                       {169.174, 21.339},
                       {134.202, 19.744},
                       {112.447, 18.291},
                       {94.572, 15.720},
                       {80.135, 12.959},
                       {70.425, 11.239},
                       {52.880, 10.203},
                       {41.196, 6.562},
                       {39.712, 6.282},
                       {40.418, 6.629}};
  }

  void collectMC(o2::aod::UDMcCollisions const& mcCollisions, o2::aod::UDMcParticles const& mcParticles)
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
      registry.fill(HIST("PairsMC/Mass"), p.M());
      registry.fill(HIST("PairsMC/Eta"), p1.Eta());
      registry.fill(HIST("PairsMC/Eta"), p2.Eta());
      mcPartsFiltered.clear();
    }
  }

  template <typename TTrack>
  bool checkTPCPID(TTrack& tr1, TTrack& tr2, float& m1, float& m2, bool* pidFlags, int32_t pdg1 = -1, int32_t pdg2 = -1)
  {
    bool pass = false;

    float tpcSignal1 = tr1.tpcSignal();
    float tpcSignal2 = tr2.tpcSignal();
    float p1 = std::sqrt(tr1.px() * tr1.px() + tr1.py() * tr1.py() + tr1.pz() * tr1.pz());
    float p2 = std::sqrt(tr2.px() * tr2.px() + tr2.py() * tr2.py() + tr2.pz() * tr2.pz());

    // using manually fitted parameters
    const int32_t nbins = 11;
    double p_l[nbins] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2., 3., 4.};
    double p_h[nbins] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2., 3., 4., 10.};

    int32_t bin1;
    int32_t bin2;

    for (int32_t ib = 0; ib < nbins; ib++) {
      bin1 = ib; // if p > 10. we will get the last bin
      if (p1 > p_l[ib] && p1 < p_h[ib])
        break;
    }

    for (int32_t ib = 0; ib < nbins; ib++) {
      bin2 = ib; // if p > 10. we will get the last bin
      if (p2 > p_l[ib] && p2 < p_h[ib])
        break;
    }

    float nTPCSigmaEl1 = (fMeansSigmas[0][bin1][0] - tpcSignal1) / fMeansSigmas[0][bin1][1];
    float nTPCSigmaPi1 = (fMeansSigmas[1][bin1][0] - tpcSignal1) / fMeansSigmas[1][bin1][1];

    float nTPCSigmaEl2 = (fMeansSigmas[0][bin2][0] - tpcSignal2) / fMeansSigmas[0][bin2][1];
    float nTPCSigmaPi2 = (fMeansSigmas[1][bin2][0] - tpcSignal2) / fMeansSigmas[1][bin2][1];

    // QA histograms for MC
    if (std::abs(pdg1) == 11)
      registry.fill(HIST("PID_TPC/nSigmaEl"), nTPCSigmaEl1);
    if (std::abs(pdg2) == 11)
      registry.fill(HIST("PID_TPC/nSigmaEl"), nTPCSigmaEl2);
    if (std::abs(pdg1) == 211)
      registry.fill(HIST("PID_TPC/nSigmaPi"), nTPCSigmaPi1);
    if (std::abs(pdg2) == 211)
      registry.fill(HIST("PID_TPC/nSigmaPi"), nTPCSigmaPi2);

    // "soft" selection
    bool isEl1 = std::abs(nTPCSigmaEl1) < 3.0f && std::abs(nTPCSigmaEl1) < std::abs(nTPCSigmaPi1);
    bool isEl2 = std::abs(nTPCSigmaEl2) < 3.0f && std::abs(nTPCSigmaEl2) < std::abs(nTPCSigmaPi2);
    bool isPi1 = std::abs(nTPCSigmaPi1) < 3.0f && std::abs(nTPCSigmaEl1) > std::abs(nTPCSigmaPi1);
    bool isPi2 = std::abs(nTPCSigmaPi2) < 3.0f && std::abs(nTPCSigmaEl2) > std::abs(nTPCSigmaPi2);

    // PID flags to check sigmas later
    pidFlags[0] = isEl1;
    pidFlags[1] = isEl2;
    pidFlags[2] = isPi1;
    pidFlags[3] = isPi2;

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

  void fillMassDistr(float m, float mmc, std::map<int32_t, bool>& selFlags)
  {
    if (mmc < 0) { // just fill reco mass, if not using MC
      mmc = m;
    }
    registry.fill(HIST("ProcessCandidate/PairMass/All"), m);
    // unlike-sign
    bool selector = selFlags[kSelUnlikeSign];
    if (selector) {
      registry.fill(HIST("ProcessCandidate/PairMass/UnlikeSign"), m);
    }
    // unlike sign + no FT0
    selector = selector && selFlags[kSelNoFT0];
    if (selector) {
      registry.fill(HIST("ProcessCandidate/PairMass/NoFT0"), m);
    }
    // unlike sign + no FT0 + delta(DCA)
    bool selectorDCA = selector && selFlags[kSelDDCA];
    if (selectorDCA) {
      registry.fill(HIST("ProcessCandidate/PairMass/DDCA"), m);
      registry.fill(HIST("ProcessCandidate/PairMass/DDCA_MC"), mmc);
    }
    // unlike sign + no FT0 + [ideal]PID
    bool selectorIdeal = selector && selFlags[kSelIdealPID];
    selector = selector && selFlags[kSelPID];
    if (selector) {
      registry.fill(HIST("ProcessCandidate/PairMass/PID"), m);
      registry.fill(HIST("ProcessCandidate/PairMass/PID_MC"), mmc);
    }
    if (selectorIdeal) {
      registry.fill(HIST("ProcessCandidate/PairMass/IdealPID"), m);
      registry.fill(HIST("ProcessCandidate/PairMass/IdealPID_MC"), mmc);
    }
    if (selector && selFlags[kSelIsNotFake]) {
      registry.fill(HIST("ProcessCandidate/PairMass/PIDSig"), m);
      registry.fill(HIST("ProcessCandidate/PairMass/PIDSig_MC"), mmc);
    }
    if (selectorIdeal && selFlags[kSelIsNotFake]) {
      registry.fill(HIST("ProcessCandidate/PairMass/IdealPIDSig"), m);
      registry.fill(HIST("ProcessCandidate/PairMass/IdealPIDSig_MC"), mmc);
    }
  }

  void fillPtDistr(float pt, float ptmc, std::map<int32_t, bool>& selFlags)
  {
    if (ptmc < 0) { // just fill reco pt, if not using MC
      ptmc = pt;
    }
    registry.fill(HIST("ProcessCandidate/TargetPt/All"), pt);
    // unlike-sign
    bool selector = selFlags[kSelUnlikeSign];
    if (selector) {
      registry.fill(HIST("ProcessCandidate/TargetPt/UnlikeSign"), pt);
    }
    // unlike sign + no FT0
    selector = selector && selFlags[kSelNoFT0];
    if (selector) {
      registry.fill(HIST("ProcessCandidate/TargetPt/NoFT0"), pt);
    }
    // unlike sign + no FT0 + delta(DCA)
    bool selectorDCA = selector && selFlags[kSelDDCA];
    if (selectorDCA) {
      registry.fill(HIST("ProcessCandidate/TargetPt/DDCA"), pt);
      registry.fill(HIST("ProcessCandidate/TargetPt/DDCA_MC"), ptmc);
    }
    // unlike sign + no FT0 + [ideal]PID
    bool selectorIdeal = selector && selFlags[kSelIdealPID];
    selector = selector && selFlags[kSelPID];
    if (selector) {
      registry.fill(HIST("ProcessCandidate/TargetPt/PID"), pt);
      registry.fill(HIST("ProcessCandidate/TargetPt/PID_MC"), ptmc);
    }
    if (selectorIdeal) {
      registry.fill(HIST("ProcessCandidate/TargetPt/IdealPID"), pt);
      registry.fill(HIST("ProcessCandidate/TargetPt/IdealPID_MC"), ptmc);
    }
    if (selector && selFlags[kSelIsNotFake]) {
      registry.fill(HIST("ProcessCandidate/TargetPt/PIDSig"), pt);
      registry.fill(HIST("ProcessCandidate/TargetPt/PIDSig_MC"), ptmc);
    }
    if (selectorIdeal && selFlags[kSelIsNotFake]) {
      registry.fill(HIST("ProcessCandidate/TargetPt/IdealPIDSig"), pt);
      registry.fill(HIST("ProcessCandidate/TargetPt/IdealPIDSig_MC"), ptmc);
    }
  }

  // naive DCA check
  template <typename TTrack>
  bool checkDDCA(TTrack& tr1, TTrack& tr2)
  {
    float dDcaZ = tr1.dcaZ() - tr2.dcaZ();
    float dDcaXY = tr1.dcaXY() - tr2.dcaXY();
    float r = std::sqrt(dDcaZ * dDcaZ + dDcaXY * dDcaXY);
    bool pass = r < 1.5f;
    return pass;
  }

  template <bool isCentral, typename TTrack1, typename TTrack2>
  void processCandidate(o2::aod::UDCollision const& cand, TTrack1& tr1, TTrack2& tr2, o2::aod::UDMcParticles* mcParticles)
  {
    std::map<int32_t, bool> selFlags; // holder of selection flags
    float mmc = -1;
    float pt1mc = -1;
    float pt2mc = -1;
    int32_t pdg1 = -1;
    int32_t pdg2 = -1;
    //
    if (fIsMC) {
      TLorentzVector mcP1, mcP2;
      int32_t mcPartId1 = tr1.udMcParticleId();
      int32_t mcPartId2 = tr2.udMcParticleId();
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
    if constexpr (isCentral) {
      selFlags[kSelPID] = checkTPCPID(tr1, tr2, m1, m2, pidFlags, pdg1, pdg2);
      selFlags[kSelDDCA] = checkDDCA(tr1, tr2);
    }
    // unlike-sign tracks requirement
    selFlags[kSelUnlikeSign] = tr1.sign() * tr2.sign() < 0;
    // check FT0 signal
    if (cand.hasFT0()) {
      bool hasNoFT0 = true;
      if constexpr (isCentral) {
        hasNoFT0 = false;
      } else {
        // if there is a signal, candidate passes if timeA is dummy
        // and timeC is between +/- 1 ns
        bool checkA = std::abs(cand.timeFT0A() - ft0DummyTime) < 1e-3;
        bool checkC = cand.timeFT0C() > -1. && cand.timeFT0C() < 1.;
        if (!(checkA && checkC)) {
          hasNoFT0 = false;
        }
      }
      selFlags[kSelNoFT0] = hasNoFT0;
    }
    // selection counters
    if (selFlags[kSelIdealPID]) { // has real meaning for MC only
      registry.fill(HIST("ProcessCandidate/SelCounter"), kSelIdealPID, 1);
    }
    if (selFlags[kSelIsNotFake]) { // has real meaning for MC only
      registry.fill(HIST("ProcessCandidate/SelCounter"), kSelIsNotFake, 1);
    }
    if (selFlags[kSelUnlikeSign]) {
      registry.fill(HIST("ProcessCandidate/SelCounter"), kSelUnlikeSign, 1);
    }
    if (selFlags[kSelNoFT0]) {
      registry.fill(HIST("ProcessCandidate/SelCounter"), kSelNoFT0, 1);
    }
    if (selFlags[kSelPID]) {
      registry.fill(HIST("ProcessCandidate/SelCounter"), kSelPID, 1);
    }
    if (selFlags[kSelDDCA]) {
      registry.fill(HIST("ProcessCandidate/SelCounter"), kSelDDCA, 1);
    }
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m1);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m2);
    // collect mass distributions if needed
    if (fHistSwitch == 0 || fHistSwitch == 2) {
      float m = (p1 + p2).M();
      fillMassDistr(m, mmc, selFlags);
    }
    // collect pt distributions if needed
    if (fHistSwitch == 1 || fHistSwitch == 2) {
      float pt1 = p1.Pt();
      float pt2 = p2.Pt();
      bool fill1, fill2;
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
    if constexpr (isCentral) {
      registry.fill(HIST("ProcessCandidate/TPCSignals/All"), tr1.tpcSignal(), tr2.tpcSignal());
      if (selFlags[kSelIdealPID]) {
        registry.fill(HIST("ProcessCandidate/TPCSignals/IdealPID"), tr1.tpcSignal(), tr2.tpcSignal());
      }
      if (selFlags[kSelIdealPID] && selFlags[kSelIsNotFake]) { // for MC
        registry.fill(HIST("ProcessCandidate/TPCSignals/IdealPIDSig"), tr1.tpcSignal(), tr2.tpcSignal());
      }
      if (selFlags[kSelPID]) {
        registry.fill(HIST("ProcessCandidate/TPCSignals/PID"), tr1.tpcSignal(), tr2.tpcSignal());
      }
      if (selFlags[kSelPID] && selFlags[kSelIsNotFake]) { // for MC
        registry.fill(HIST("ProcessCandidate/TPCSignals/PIDSig"), tr1.tpcSignal(), tr2.tpcSignal());
      }
    }
  }

  // process candidates with 2 muon tracks
  void processFwdMC(o2::aod::UDCollisions const& eventCandidates,
                    FwdTracks const& fwdTracks,
                    o2::aod::UDMcCollisions const& mcCollisions,
                    o2::aod::UDMcParticles& mcParticles)
  {
    fIsMC = true;

    collectMC(mcCollisions, mcParticles);

    // assuming that candidates have exatly 2 muon tracks and 0 barrel tracks
    for (const auto& cand : eventCandidates) {
      auto candID = cand.globalIndex();
      auto fwdTracksPerCand = fwdTracks.select(o2::aod::udfwdtrack::udCollisionId == candID);
      const auto& tr1 = fwdTracksPerCand.iteratorAt(0);
      const auto& tr2 = fwdTracksPerCand.iteratorAt(1);
      processCandidate<false>(cand, tr1, tr2, &mcParticles);
    }
  }

  // process candidates with 2 muon tracks
  void processFwd(o2::aod::UDCollisions const& eventCandidates,
                  FwdTracks const& fwdTracks,
                  o2::aod::UDMcCollisions const& mcCollisions)
  {
    fIsMC = false;

    // assuming that candidates have exatly 2 muon tracks and 0 barrel tracks
    for (const auto& cand : eventCandidates) {
      auto candID = cand.globalIndex();
      auto fwdTracksPerCand = fwdTracks.select(o2::aod::udfwdtrack::udCollisionId == candID);
      const auto& tr1 = fwdTracksPerCand.iteratorAt(0);
      const auto& tr2 = fwdTracksPerCand.iteratorAt(1);
      processCandidate<false>(cand, tr1, tr2, (o2::aod::UDMcParticles*)nullptr);
    }
  }

  // process candidates with 1 muon and 1 barrel tracks
  void processSemiFwdMC(o2::aod::UDCollisions const& eventCandidates,
                        FwdTracks const& fwdTracks,
                        BarrelTracks const& barTracks,
                        o2::aod::UDMcCollisions const& mcCollisions,
                        o2::aod::UDMcParticles& mcParticles)
  {
    fIsMC = true;

    collectMC(mcCollisions, mcParticles);

    // assuming that candidates have exatly 1 muon track and 1 barrel track
    for (const auto& cand : eventCandidates) {
      auto candID = cand.globalIndex();
      auto fwdTracksPerCand = fwdTracks.select(o2::aod::udfwdtrack::udCollisionId == candID);
      auto barTracksPerCand = barTracks.select(o2::aod::udtrack::udCollisionId == candID);
      const auto& tr1 = fwdTracksPerCand.iteratorAt(0);
      const auto& tr2 = barTracksPerCand.iteratorAt(0);
      processCandidate<false>(cand, tr1, tr2, &mcParticles);
    }
  }

  // process candidates with 1 muon and 1 barrel tracks
  void processSemiFwd(o2::aod::UDCollisions const& eventCandidates,
                      FwdTracks const& fwdTracks,
                      BarrelTracks const& barTracks,
                      o2::aod::UDMcCollisions const& mcCollisions)
  {
    fIsMC = false;

    // assuming that candidates have exatly 1 muon track and 1 barrel track
    for (const auto& cand : eventCandidates) {
      auto candID = cand.globalIndex();
      auto fwdTracksPerCand = fwdTracks.select(o2::aod::udfwdtrack::udCollisionId == candID);
      auto barTracksPerCand = barTracks.select(o2::aod::udtrack::udCollisionId == candID);
      const auto& tr1 = fwdTracksPerCand.iteratorAt(0);
      const auto& tr2 = barTracksPerCand.iteratorAt(0);
      processCandidate<false>(cand, tr1, tr2, (o2::aod::UDMcParticles*)nullptr);
    }
  }

  // process candidates with 2 central barrel tracks
  void processCentralMC(o2::aod::UDCollisions const& eventCandidates,
                        BarrelTracks const& barTracks,
                        o2::aod::UDMcCollisions const& mcCollisions,
                        o2::aod::UDMcParticles& mcParticles)
  {
    fIsMC = true;

    collectMC(mcCollisions, mcParticles);

    // assuming that candidates have exatly 2 central barrel tracks
    for (const auto& cand : eventCandidates) {
      auto candID = cand.globalIndex();
      auto barTracksPerCand = barTracks.select(o2::aod::udtrack::udCollisionId == candID);
      const auto& tr1 = barTracksPerCand.iteratorAt(0);
      const auto& tr2 = barTracksPerCand.iteratorAt(1);
      processCandidate<true>(cand, tr1, tr2, &mcParticles);
    }
  }

  // process candidates with 2 central barrel tracks
  void processCentral(o2::aod::UDCollisions const& eventCandidates,
                      BarrelTracks const& barTracks,
                      o2::aod::UDMcCollisions const& mcCollisions)
  {
    fIsMC = false;

    // assuming that candidates have exatly 2 central barrel tracks
    for (const auto& cand : eventCandidates) {
      auto candID = cand.globalIndex();
      auto barTracksPerCand = barTracks.select(o2::aod::udtrack::udCollisionId == candID);
      const auto& tr1 = barTracksPerCand.iteratorAt(0);
      const auto& tr2 = barTracksPerCand.iteratorAt(1);
      processCandidate<true>(cand, tr1, tr2, (o2::aod::UDMcParticles*)nullptr);
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
