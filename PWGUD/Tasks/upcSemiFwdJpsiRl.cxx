// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file upcSemiFwdJpsiRl.cxx
/// \brief UPC semi-forward J/psi -> mu+mu- analysis pairing one forward MCH-MID muon
///        with one central-barrel track. Consumes UD tables produced by
///        upcCandProducerSemiFwd.
/// \author Roman Lavicka, roman.lavicka@cern.ch
/// \since  10.05.2026
///
/// executable name o2-analysis-ud-upc-semi-fwd-jpsi-rl

#include "PWGUD/DataModel/UDTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <TRandom3.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <vector>

// table for saving tree with semi-forward dimuon info
namespace jpsisfw
{
// dimuon
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(PhiAv, phiAv, float);
DECLARE_SOA_COLUMN(PhiCh, phiCh, float);
// positive (p) and negative (n) tracks
DECLARE_SOA_COLUMN(EnergyP, energyP, float);
DECLARE_SOA_COLUMN(Pxp, pxp, float);
DECLARE_SOA_COLUMN(Pyp, pyp, float);
DECLARE_SOA_COLUMN(Pzp, pzp, float);
DECLARE_SOA_COLUMN(Ptp, ptp, float);
DECLARE_SOA_COLUMN(Etap, etap, float);
DECLARE_SOA_COLUMN(Phip, phip, float);
DECLARE_SOA_COLUMN(IsFwdP, isFwdP, int); // 1 if the positive track is forward, 0 if barrel
DECLARE_SOA_COLUMN(EnergyN, energyN, float);
DECLARE_SOA_COLUMN(Pxn, pxn, float);
DECLARE_SOA_COLUMN(Pyn, pyn, float);
DECLARE_SOA_COLUMN(Pzn, pzn, float);
DECLARE_SOA_COLUMN(Ptn, ptn, float);
DECLARE_SOA_COLUMN(Etan, etan, float);
DECLARE_SOA_COLUMN(Phin, phin, float);
DECLARE_SOA_COLUMN(IsFwdN, isFwdN, int);
// zn
DECLARE_SOA_COLUMN(Tzna, tzna, float);
DECLARE_SOA_COLUMN(Ezna, ezna, float);
DECLARE_SOA_COLUMN(Tznc, tznc, float);
DECLARE_SOA_COLUMN(Eznc, eznc, float);
DECLARE_SOA_COLUMN(Nclass, nclass, int);
} // namespace jpsisfw

namespace o2::aod
{
DECLARE_SOA_TABLE(JpsiSemiFwdRL, "AOD", "JPSISFW",
                  jpsisfw::RunNumber,
                  jpsisfw::M, jpsisfw::Energy, jpsisfw::Px, jpsisfw::Py, jpsisfw::Pz, jpsisfw::Pt, jpsisfw::Rap, jpsisfw::Phi,
                  jpsisfw::PhiAv, jpsisfw::PhiCh,
                  jpsisfw::EnergyP, jpsisfw::Pxp, jpsisfw::Pyp, jpsisfw::Pzp, jpsisfw::Ptp, jpsisfw::Etap, jpsisfw::Phip, jpsisfw::IsFwdP,
                  jpsisfw::EnergyN, jpsisfw::Pxn, jpsisfw::Pyn, jpsisfw::Pzn, jpsisfw::Ptn, jpsisfw::Etan, jpsisfw::Phin, jpsisfw::IsFwdN,
                  jpsisfw::Tzna, jpsisfw::Ezna, jpsisfw::Tznc, jpsisfw::Eznc, jpsisfw::Nclass);
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// constants used in the forward-muon selection (mirror upcFwdJpsiRl)
const float kRAbsMin = 17.6;
const float kRAbsMid = 26.5;
const float kRAbsMax = 89.5;
const float kPDca1 = 200.;
const float kPDca2 = 200.;
const float kFwdEtaMin = -4.0;
const float kFwdEtaMax = -2.5;
const float kFwdPtMin = 0.;

const float kMaxAmpV0A = 100.;
const float kMaxZDCTime = 2.;
const float kMaxZDCTimeHisto = 10.;
const float kInvalidFloat = -999.;
const int kMaxRelBCsV0A = 1;

struct UpcSemiFwdJpsiRl {

  using CandidatesSemiFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using BarrelTracks = soa::Join<o2::aod::UDTracks, o2::aod::UDTracksExtra, o2::aod::UDTracksDCA,
                                 o2::aod::UDTracksFlags, o2::aod::UDTracksPID>;

  Produces<o2::aod::JpsiSemiFwdRL> dimuSel;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry reg0n0n{"reg0n0n", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry regXn0n{"regXn0n", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry regXnXn{"regXnXn", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // CONFIGURABLES
  static constexpr double Pi = o2::constants::math::PI;
  // pair kinematics
  Configurable<int> nBinsPt{"nBinsPt", 250, "N bins in pT histo"};
  Configurable<float> lowPt{"lowPt", 0., "lower limit in pT histo"};
  Configurable<float> highPt{"highPt", 2., "upper limit in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 500, "N bins in mass histo"};
  Configurable<float> lowMass{"lowMass", 0., "lower limit in mass histo"};
  Configurable<float> highMass{"highMass", 10., "upper limit in mass histo"};
  Configurable<int> nBinsEta{"nBinsEta", 600, "N bins in eta histo"};
  Configurable<float> lowEta{"lowEta", -10., "lower limit in eta histo"};
  Configurable<float> highEta{"highEta", 5., "upper limit in eta histo"};
  Configurable<int> nBinsRapidity{"nBinsRapidity", 500, "N bins in rapidity histo"};
  Configurable<float> lowRapidity{"lowRapidity", -4.5, "lower limit in rapidity histo"};
  Configurable<float> highRapidity{"highRapidity", 0., "upper limit in rapidity histo"};
  Configurable<int> nBinsPhi{"nBinsPhi", 600, "N bins in phi histo"};
  Configurable<float> lowPhi{"lowPhi", -Pi, "lower limit in phi histo"};
  Configurable<float> highPhi{"highPhi", Pi, "upper limit in phi histo"};
  // single-track histos
  Configurable<int> nBinsPtSingle{"nBinsPtSingle", 500, "N bins in pT histo single track"};
  Configurable<float> lowPtSingle{"lowPtSingle", 0., "lower limit in pT histo single track"};
  Configurable<float> highPtSingle{"highPtSingle", 5., "upper limit in pT histo single track"};
  Configurable<int> nBinsEtaSingle{"nBinsEtaSingle", 250, "N bins in eta histo single track"};
  Configurable<float> lowEtaSingle{"lowEtaSingle", -4.5, "lower limit in eta histo single track"};
  Configurable<float> highEtaSingle{"highEtaSingle", 1.5, "upper limit in eta histo single track"};
  Configurable<int> nBinsPhiSingle{"nBinsPhiSingle", 600, "N bins in phi histo single track"};
  Configurable<float> lowPhiSingle{"lowPhiSingle", -Pi, "lower limit in phi histo single track"};
  Configurable<float> highPhiSingle{"highPhiSingle", Pi, "upper limit in phi histo single track"};
  // ZDC
  Configurable<int> nBinsZDCen{"nBinsZDCen", 200, "N bins in ZN energy"};
  Configurable<float> lowEnZN{"lowEnZN", -50., "lower limit in ZN energy histo"};
  Configurable<float> highEnZN{"highEnZN", 250., "upper limit in ZN energy histo"};

  // central-barrel-muon selection
  Configurable<float> barEtaMin{"barEtaMin", -0.9, "barrel min eta"};
  Configurable<float> barEtaMax{"barEtaMax", 0.9, "barrel max eta"};
  Configurable<float> barPtMin{"barPtMin", 0.7, "barrel min pT (GeV/c)"};
  Configurable<int> barMinITSNCls{"barMinITSNCls", 4, "barrel min number of ITS clusters"};
  Configurable<int> barMinTPCNClsCR{"barMinTPCNClsCR", 70, "barrel min number of TPC crossed rows"};
  Configurable<float> barMaxTpcNSigmaMu{"barMaxTpcNSigmaMu", 4., "barrel max |TPC nSigma muon|"};
  Configurable<bool> barRequireTOF{"barRequireTOF", false, "require TOF on the barrel track"};
  Configurable<float> barMaxTofNSigmaMu{"barMaxTofNSigmaMu", 3., "barrel max |TOF nSigma muon| when TOF is present"};

  void init(InitContext&)
  {
    std::vector<double> ptFitBinning = {
      0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
      0.11, 0.12, 0.13, 0.14, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.50,
      0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.50,
      3.00, 3.50};

    const AxisSpec axisPt{nBinsPt, lowPt, highPt, "#it{p}_{T} GeV/#it{c}"};
    const AxisSpec axisPtFit = {ptFitBinning, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axisMass{nBinsMass, lowMass, highMass, "m_{#mu#mu} GeV/#it{c}^{2}"};
    const AxisSpec axisEta{nBinsEta, lowEta, highEta, "#eta"};
    const AxisSpec axisRapidity{nBinsRapidity, lowRapidity, highRapidity, "Rapidity"};
    const AxisSpec axisPhi{nBinsPhi, lowPhi, highPhi, "#varphi"};
    const AxisSpec axisPtSingle{nBinsPtSingle, lowPtSingle, highPtSingle, "#it{p}_{T}_{ trk} GeV/#it{c}"};
    const AxisSpec axisTimeZN{200, -10, 10, "ZDC time (ns)"};
    const AxisSpec axisEnergyZNA{nBinsZDCen, lowEnZN, highEnZN, "ZNA energy (TeV)"};
    const AxisSpec axisEnergyZNC{nBinsZDCen, lowEnZN, highEnZN, "ZNC energy (TeV)"};
    const AxisSpec axisEtaSingle{nBinsEtaSingle, lowEtaSingle, highEtaSingle, "#eta_{trk}"};
    const AxisSpec axisPhiSingle{nBinsPhiSingle, lowPhiSingle, highPhiSingle, "#varphi_{trk}"};

    // pair histos
    registry.add("hMass", "Invariant mass of mu pairs;;#counts", kTH1D, {axisMass});
    registry.add("hPt", "Transverse momentum of mu pairs;;#counts", kTH1D, {axisPt});
    registry.add("hPtFit", "Transverse momentum of mu pairs;;#counts", kTH1D, {axisPtFit});
    registry.add("hEta", "Pseudorapidity of mu pairs;;#counts", kTH1D, {axisEta});
    registry.add("hRapidity", "Rapidity of mu pairs;;#counts", kTH1D, {axisRapidity});
    registry.add("hPhi", "#varphi of mu pairs;;#counts", kTH1D, {axisPhi});
    registry.add("hCharge", "Charge;;;#counts", kTH1D, {{5, -2.5, 2.5}});
    registry.add("hContrib", "hContrib;;#counts", kTH1D, {{20, -0.5, 19.5}});
    registry.add("hEvSign", "Sum of the charges of all the tracks in each event;;#counts", kTH1D, {{5, -2.5, 2.5}});
    registry.add("hSameSign", "hSameSign;;#counts", kTH1D, {{20, -0.5, 19.5}});
    registry.add("hPhiCharge", "#phi #it{charge}", kTH1D, {axisPhi});
    registry.add("hPhiAverage", "#phi #it{average}", kTH1D, {axisPhi});

    // forward / barrel single-track histos
    registry.add("hPtTrkFwd", "Pt of forward muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hPtTrkBar", "Pt of barrel muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hEtaTrkFwd", "#eta of forward muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hEtaTrkBar", "#eta of barrel muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hPhiTrkFwd", "#varphi of forward muons;;#counts", kTH1D, {axisPhiSingle});
    registry.add("hPhiTrkBar", "#varphi of barrel muons;;#counts", kTH1D, {axisPhiSingle});

    // ZDC
    registry.add("hTimeZNA", "ZNA Times;;#counts", kTH1D, {axisTimeZN});
    registry.add("hTimeZNC", "ZNC Times;;#counts", kTH1D, {axisTimeZN});
    registry.add("hEnergyZN", "ZNA vs ZNC energy", kTH2D, {axisEnergyZNA, axisEnergyZNC});

    // neutron classes
    reg0n0n.add("hMass", "Invariant mass of mu pairs - 0n0n;;#counts", kTH1D, {axisMass});
    reg0n0n.add("hPt", "Transverse momentum of mu pairs - 0n0n;;#counts", kTH1D, {axisPt});
    reg0n0n.add("hEta", "Pseudorapidity of mu pairs - 0n0n;;#counts", kTH1D, {axisEta});
    reg0n0n.add("hRapidity", "Rapidity of mu pairs - 0n0n;;#counts", kTH1D, {axisRapidity});
    reg0n0n.add("hPtFit", "Transverse momentum of mu pairs - 0n0n;;#counts", kTH1D, {axisPtFit});

    regXn0n.add("hMass", "Invariant mass of mu pairs - Xn0n;;#counts", kTH1D, {axisMass});
    regXn0n.add("hPt", "Transverse momentum of mu pairs - Xn0n;;#counts", kTH1D, {axisPt});
    regXn0n.add("hEta", "Pseudorapidity of mu pairs - Xn0n;;#counts", kTH1D, {axisEta});
    regXn0n.add("hRapidity", "Rapidity of mu pairs - Xn0n;;#counts", kTH1D, {axisRapidity});
    regXn0n.add("hPtFit", "Transverse momentum of mu pairs - Xn0n;;#counts", kTH1D, {axisPtFit});

    regXnXn.add("hMass", "Invariant mass of mu pairs - XnXn;;#counts", kTH1D, {axisMass});
    regXnXn.add("hPt", "Transverse momentum of mu pairs - XnXn;;#counts", kTH1D, {axisPt});
    regXnXn.add("hEta", "Pseudorapidity of mu pairs - XnXn;;#counts", kTH1D, {axisEta});
    regXnXn.add("hRapidity", "Rapidity of mu pairs - XnXn;;#counts", kTH1D, {axisRapidity});
    regXnXn.add("hPtFit", "Transverse momentum of mu pairs - XnXn;;#counts", kTH1D, {axisPtFit});
  }

  // FUNCTIONS
  using LorentzVec = ROOT::Math::PxPyPzMVector;

  // collect tracks per candidate into a map
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

  struct ZDCinfo {
    float timeA;
    float timeC;
    float enA;
    float enC;
  };

  void collectCandZDCInfo(std::unordered_map<int32_t, ZDCinfo>& zdcPerCand, o2::aod::UDZdcsReduced const& ZDCs)
  {
    for (const auto& zdc : ZDCs) {
      int32_t candId = zdc.udCollisionId();
      if (candId < 0) {
        continue;
      }
      zdcPerCand[candId].timeA = zdc.timeZNA();
      zdcPerCand[candId].timeC = zdc.timeZNC();
      zdcPerCand[candId].enA = zdc.energyCommonZNA();
      zdcPerCand[candId].enC = zdc.energyCommonZNC();

      if (std::isinf(zdcPerCand[candId].timeA))
        zdcPerCand[candId].timeA = kInvalidFloat;
      if (std::isinf(zdcPerCand[candId].timeC))
        zdcPerCand[candId].timeC = kInvalidFloat;
      if (std::isinf(zdcPerCand[candId].enA))
        zdcPerCand[candId].enA = kInvalidFloat;
      if (std::isinf(zdcPerCand[candId].enC))
        zdcPerCand[candId].enC = kInvalidFloat;
    }
  }

  // forward-muon selection (MCH-MID muon, like upcFwdJpsiRl)
  template <typename TTrack>
  bool isFwdMuonSelected(const TTrack& fwdTrack)
  {
    float rAbs = fwdTrack.rAtAbsorberEnd();
    float pDca = fwdTrack.pDca();
    auto mMu = o2::constants::physics::MassMuon;
    LorentzVec p(fwdTrack.px(), fwdTrack.py(), fwdTrack.pz(), mMu);
    float eta = p.Eta();
    float pt = p.Pt();
    float pDcaMax = rAbs < kRAbsMid ? kPDca1 : kPDca2;

    if (eta < kFwdEtaMin || eta > kFwdEtaMax)
      return false;
    if (pt < kFwdPtMin)
      return false;
    if (rAbs < kRAbsMin || rAbs > kRAbsMax)
      return false;
    if (pDca > pDcaMax)
      return false;
    if (fwdTrack.chi2MatchMCHMID() <= 0) // require MID match
      return false;
    return true;
  }

  // barrel-muon selection
  template <typename TTrack>
  bool isBarrelMuonSelected(const TTrack& barTrack)
  {
    auto mMu = o2::constants::physics::MassMuon;
    LorentzVec p(barTrack.px(), barTrack.py(), barTrack.pz(), mMu);
    float eta = p.Eta();
    float pt = p.Pt();

    if (pt < barPtMin)
      return false;
    if (eta < barEtaMin || eta > barEtaMax)
      return false;
    if (barTrack.itsNCls() < static_cast<uint8_t>(barMinITSNCls))
      return false;
    if (barTrack.tpcNClsCrossedRows() < static_cast<int16_t>(barMinTPCNClsCR))
      return false;
    if (std::abs(barTrack.tpcNSigmaMu()) > barMaxTpcNSigmaMu)
      return false;
    if (barRequireTOF && !barTrack.hasTOF())
      return false;
    if (barTrack.hasTOF() && std::abs(barTrack.tofNSigmaMu()) > barMaxTofNSigmaMu)
      return false;
    return true;
  }

  // azimuth anisotropy phi
  void computePhiAnis(LorentzVec p1, LorentzVec p2, int sign1, float& phiAverage, float& phiCharge)
  {
    auto tSum = p1 + p2;
    float halfUnity = 0.5;
    decltype(tSum) tDiffCh, tDiffAv;
    if (sign1 > 0) {
      tDiffCh = p1 - p2;
      if (gRandom->Rndm() > halfUnity)
        tDiffAv = p1 - p2;
      else
        tDiffAv = p2 - p1;
    } else {
      tDiffCh = p2 - p1;
      if (gRandom->Rndm() > halfUnity)
        tDiffAv = p2 - p1;
      else
        tDiffAv = p1 - p2;
    }
    phiAverage = ROOT::Math::VectorUtil::DeltaPhi(tSum, tDiffAv);
    phiCharge = ROOT::Math::VectorUtil::DeltaPhi(tSum, tDiffCh);
  }

  // process a single (forward, barrel) pair
  void processPair(CandidatesSemiFwd::iterator const& cand,
                   ForwardTracks::iterator const& fwd, BarrelTracks::iterator const& bar,
                   ZDCinfo const& zdc)
  {
    // V0 selection
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= kMaxRelBCsV0A) {
        if (ampsV0A[i] > kMaxAmpV0A)
          return;
      }
    }

    // opposite charge only
    int pairCharge = fwd.sign() + bar.sign();
    if (pairCharge != 0) {
      registry.fill(HIST("hSameSign"), cand.numContrib());
      return;
    }

    // track selection
    if (!isFwdMuonSelected(fwd))
      return;
    if (!isBarrelMuonSelected(bar))
      return;

    // form Lorentz vectors
    auto mMu = o2::constants::physics::MassMuon;
    LorentzVec pFwd(fwd.px(), fwd.py(), fwd.pz(), mMu);
    LorentzVec pBar(bar.px(), bar.py(), bar.pz(), mMu);
    auto p = pFwd + pBar;

    // pair-level cuts
    if (p.M() < lowMass || p.M() > highMass)
      return;
    if (p.Pt() < lowPt || p.Pt() > highPt)
      return;
    if (p.Rapidity() < lowRapidity || p.Rapidity() > highRapidity)
      return;

    // azimuth anisotropy
    float phiAverage = 0;
    float phiCharge = 0;
    computePhiAnis(pFwd, pBar, fwd.sign(), phiAverage, phiCharge);

    // ZDC info histograms
    if (std::abs(zdc.timeA) < kMaxZDCTimeHisto)
      registry.fill(HIST("hTimeZNA"), zdc.timeA);
    if (std::abs(zdc.timeC) < kMaxZDCTimeHisto)
      registry.fill(HIST("hTimeZNC"), zdc.timeC);
    registry.fill(HIST("hEnergyZN"), zdc.enA, zdc.enC);

    // neutron classes
    bool neutronA = std::abs(zdc.timeA) < kMaxZDCTime && !std::isinf(zdc.timeA);
    bool neutronC = std::abs(zdc.timeC) < kMaxZDCTime && !std::isinf(zdc.timeC);
    int znClass = -1;

    if (!neutronC && !neutronA) {
      znClass = 1;
      reg0n0n.fill(HIST("hMass"), p.M());
      reg0n0n.fill(HIST("hPt"), p.Pt());
      reg0n0n.fill(HIST("hPtFit"), p.Pt());
      reg0n0n.fill(HIST("hEta"), p.Eta());
      reg0n0n.fill(HIST("hRapidity"), p.Rapidity());
    } else if (neutronA ^ neutronC) {
      znClass = neutronA ? 2 : 3;
      regXn0n.fill(HIST("hMass"), p.M());
      regXn0n.fill(HIST("hPt"), p.Pt());
      regXn0n.fill(HIST("hPtFit"), p.Pt());
      regXn0n.fill(HIST("hEta"), p.Eta());
      regXn0n.fill(HIST("hRapidity"), p.Rapidity());
    } else if (neutronA && neutronC) {
      znClass = 4;
      regXnXn.fill(HIST("hMass"), p.M());
      regXnXn.fill(HIST("hPt"), p.Pt());
      regXnXn.fill(HIST("hPtFit"), p.Pt());
      regXnXn.fill(HIST("hEta"), p.Eta());
      regXnXn.fill(HIST("hRapidity"), p.Rapidity());
    }

    // single-track histos by side
    registry.fill(HIST("hPtTrkFwd"), pFwd.Pt());
    registry.fill(HIST("hPtTrkBar"), pBar.Pt());
    registry.fill(HIST("hEtaTrkFwd"), pFwd.Eta());
    registry.fill(HIST("hEtaTrkBar"), pBar.Eta());
    registry.fill(HIST("hPhiTrkFwd"), pFwd.Phi());
    registry.fill(HIST("hPhiTrkBar"), pBar.Phi());

    // pair histos
    registry.fill(HIST("hContrib"), cand.numContrib());
    registry.fill(HIST("hEvSign"), cand.netCharge());
    registry.fill(HIST("hMass"), p.M());
    registry.fill(HIST("hPt"), p.Pt());
    registry.fill(HIST("hPtFit"), p.Pt());
    registry.fill(HIST("hEta"), p.Eta());
    registry.fill(HIST("hRapidity"), p.Rapidity());
    registry.fill(HIST("hPhi"), p.Phi());
    registry.fill(HIST("hCharge"), fwd.sign());
    registry.fill(HIST("hCharge"), bar.sign());
    registry.fill(HIST("hPhiAverage"), phiAverage);
    registry.fill(HIST("hPhiCharge"), phiCharge);

    // assign positive/negative leg
    bool fwdIsPos = fwd.sign() > 0;
    const LorentzVec& pPos = fwdIsPos ? pFwd : pBar;
    const LorentzVec& pNeg = fwdIsPos ? pBar : pFwd;
    int isFwdPos = fwdIsPos ? 1 : 0;
    int isFwdNeg = fwdIsPos ? 0 : 1;

    dimuSel(cand.runNumber(),
            p.M(), p.E(), p.Px(), p.Py(), p.Pz(), p.Pt(), p.Rapidity(), p.Phi(),
            phiAverage, phiCharge,
            pPos.E(), pPos.Px(), pPos.Py(), pPos.Pz(), pPos.Pt(), pPos.Eta(), pPos.Phi(), isFwdPos,
            pNeg.E(), pNeg.Px(), pNeg.Py(), pNeg.Pz(), pNeg.Pt(), pNeg.Eta(), pNeg.Phi(), isFwdNeg,
            zdc.timeA, zdc.enA, zdc.timeC, zdc.enC, znClass);
  }

  // PROCESS FUNCTION
  void processData(CandidatesSemiFwd const& eventCandidates,
                   o2::aod::UDZdcsReduced const& ZDCs,
                   ForwardTracks const& fwdTracks,
                   BarrelTracks const& barrelTracks)
  {
    std::unordered_map<int32_t, std::vector<int32_t>> fwdTracksPerCand;
    collectCandIDs(fwdTracksPerCand, fwdTracks);

    std::unordered_map<int32_t, std::vector<int32_t>> barTracksPerCand;
    collectCandIDs(barTracksPerCand, barrelTracks);

    std::unordered_map<int32_t, ZDCinfo> zdcPerCand;
    collectCandZDCInfo(zdcPerCand, ZDCs);

    // loop over candidates that have at least one forward track
    for (const auto& [candID, fwdIds] : fwdTracksPerCand) {
      auto itBar = barTracksPerCand.find(candID);
      if (itBar == barTracksPerCand.end())
        continue; // need at least one barrel companion
      const auto& barIds = itBar->second;
      auto cand = eventCandidates.iteratorAt(candID);

      ZDCinfo zdc;
      auto itZdc = zdcPerCand.find(candID);
      if (itZdc != zdcPerCand.end()) {
        zdc = itZdc->second;
      } else {
        zdc.timeA = kInvalidFloat;
        zdc.timeC = kInvalidFloat;
        zdc.enA = kInvalidFloat;
        zdc.enC = kInvalidFloat;
      }

      for (const auto& fwdIdx : fwdIds) {
        for (const auto& barIdx : barIds) {
          auto fwd = fwdTracks.iteratorAt(fwdIdx);
          auto bar = barrelTracks.iteratorAt(barIdx);
          processPair(cand, fwd, bar, zdc);
        }
      }
    }
  }

  PROCESS_SWITCH(UpcSemiFwdJpsiRl, processData, "Process semi-forward J/psi candidates", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcSemiFwdJpsiRl>(cfgc),
  };
}
