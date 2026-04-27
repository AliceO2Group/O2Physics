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

/// \file upcFwdJpsiRL.cxx
/// \brief UPC forward J/psi analysis with configurable track-type candidate selection
///
/// Candidate types (configured via candidateType):
///   0 = GlobalMuon-GlobalMuon: both tracks are GlobalMuonTrack (enum 0) or GlobalMuonTrackOtherMatch (enum 1)
///   1 = Mixed: at least one GlobalMuonTrack (enum 0,1) and the other MuonStandaloneTrack (enum 3)
///   2 = MuonStandalone-MuonStandalone: both tracks are MuonStandaloneTrack (enum 3)

/// executable name o2-analysis-ud-upc-fwd-jpsi-rl

/// \author Roman Lavicka <roman.lavicka@cern.ch>

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
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::fwdtrack;

// constants used in the track selection
const float kRAbsMin = 17.6;
const float kRAbsMid = 26.5;
const float kRAbsMax = 89.5;
const float kPDca1 = 200.;
const float kPDca2 = 200.;
const float kEtaMin = -4.0;
const float kEtaMax = -2.5;
const float kPtMin = 0.;

const float kMaxAmpV0A = 100.;
const float kMaxZDCTime = 2.;
const float kMaxZDCTimeHisto = 10.;

struct UpcFwdJpsiRL {

  Service<o2::framework::O2DatabasePDG> pdg;

  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;

  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry reg0n0n{"reg0n0n", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry regXn0n{"regXn0n", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry regXnXn{"regXnXn", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // CONFIGURABLES
  static constexpr double Pi = o2::constants::math::PI;
  // candidate type selection
  // 0 = GlobalMuon-GlobalMuon (both tracks enum 0 or 1)
  // 1 = Mixed (one GlobalMuon enum 0,1 + one MuonStandalone enum 3)
  // 2 = MuonStandalone-MuonStandalone (both tracks enum 3)
  Configurable<int> candidateType{"candidateType", 0, "Candidate type: 0=GlobalMuon-GlobalMuon, 1=Mixed, 2=MuonStandalone-MuonStandalone"};
  // pT of muon pairs
  Configurable<int> nBinsPt{"nBinsPt", 250, "N bins in pT histo"};
  Configurable<float> lowPt{"lowPt", 0., "lower limit in pT histo"};
  Configurable<float> highPt{"highPt", 2, "upper limit in pT histo"};
  // mass of muon pairs
  Configurable<int> nBinsMass{"nBinsMass", 500, "N bins in mass histo"};
  Configurable<float> lowMass{"lowMass", 0., "lower limit in mass histo"};
  Configurable<float> highMass{"highMass", 10., "upper limit in mass histo"};
  // eta of muon pairs
  Configurable<int> nBinsEta{"nBinsEta", 600, "N bins in eta histo"};
  Configurable<float> lowEta{"lowEta", -10., "lower limit in eta histo"};
  Configurable<float> highEta{"highEta", -2., "upper limit in eta histo"};
  // rapidity of muon pairs
  Configurable<int> nBinsRapidity{"nBinsRapidity", 250, "N bins in rapidity histo"};
  Configurable<float> lowRapidity{"lowRapidity", -4.5, "lower limit in rapidity histo"};
  Configurable<float> highRapidity{"highRapidity", -2., "upper limit in rapidity histo"};
  // phi of muon pairs
  Configurable<int> nBinsPhi{"nBinsPhi", 600, "N bins in phi histo"};
  Configurable<float> lowPhi{"lowPhi", -Pi, "lower limit in phi histo"};
  Configurable<float> highPhi{"highPhi", Pi, "upper limit in phi histo"};
  // pT of single muons
  Configurable<int> nBinsPtSingle{"nBinsPtSingle", 500, "N bins in pT histo single muon"};
  Configurable<float> lowPtSingle{"lowPtSingle", 0., "lower limit in pT histo single muon"};
  Configurable<float> highPtSingle{"highPtSingle", 2., "upper limit in pT histo single muon"};
  // eta of single muons
  Configurable<int> nBinsEtaSingle{"nBinsEtaSingle", 250, "N bins in eta histo single muon"};
  Configurable<float> lowEtaSingle{"lowEtaSingle", -4.5, "lower limit in eta histo single muon"};
  Configurable<float> highEtaSingle{"highEtaSingle", -2., "upper limit in eta histo single muon"};
  // phi of single muons
  Configurable<int> nBinsPhiSingle{"nBinsPhiSingle", 600, "N bins in phi histo single muon"};
  Configurable<float> lowPhiSingle{"lowPhiSingle", -Pi, "lower limit in phi histo single muon"};
  Configurable<float> highPhiSingle{"highPhiSingle", Pi, "upper limit in phi histo single muon"};
  // ZDC
  Configurable<int> nBinsZDCen{"nBinsZDCen", 200, "N bins in ZN energy"};
  Configurable<float> lowEnZN{"lowEnZN", -50., "lower limit in ZN energy histo"};
  Configurable<float> highEnZN{"highEnZN", 250., "upper limit in ZN energy histo"};

  void init(InitContext&)
  {
    // binning of pT axis for fit
    std::vector<double> ptFitBinning = {
      0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
      0.11, 0.12, 0.13, 0.14, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.50,
      0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.50,
      3.00, 3.50};

    std::vector<double> ptFitBinningHalfWidth = {
      0.00, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05,
      0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10,
      0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15,
      0.1625, 0.175, 0.1875, 0.20, 0.225, 0.25, 0.275, 0.30, 0.35, 0.40,
      0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
      1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.25,
      2.50, 2.75, 3.00, 3.25, 3.50};

    // axes
    const AxisSpec axisPt{nBinsPt, lowPt, highPt, "#it{p}_{T} GeV/#it{c}"};
    const AxisSpec axisPtFit = {ptFitBinning, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axisPtFit2 = {ptFitBinningHalfWidth, "#it{p}_{T} (GeV/c)"};
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

    // histos
    registry.add("hMass", "Invariant mass of muon pairs;;#counts", kTH1D, {axisMass});
    registry.add("hPt", "Transverse momentum of muon pairs;;#counts", kTH1D, {axisPt});
    registry.add("hPtFit", "Transverse momentum of muon pairs;;#counts", kTH1D, {axisPtFit});
    registry.add("hPtFit2", "Transverse momentum of muon pairs;;#counts", kTH1D, {axisPtFit2});
    registry.add("hEta", "Pseudorapidty of muon pairs;;#counts", kTH1D, {axisEta});
    registry.add("hRapidity", "Rapidty of muon pairs;;#counts", kTH1D, {axisRapidity});
    registry.add("hPhi", "#varphi of muon pairs;;#counts", kTH1D, {axisPhi});
    registry.add("hCharge", "Charge;;;#counts", kTH1D, {{5, -2.5, 2.5}});
    registry.add("hContrib", "hContrib;;#counts", kTH1D, {{6, -0.5, 5.5}});
    registry.add("hEvSign", "Sum of the charges of all the tracks in each event;;#counts", kTH1D, {{5, -2.5, 2.5}});
    registry.add("hPtTrkPos", "Pt of positive muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hPtTrkNeg", "Pt of negative muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hEtaTrkPos", "#eta of positive muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hEtaTrkNeg", "#eta of negative muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hPhiTrkPos", "#varphi of positive muons;;#counts", kTH1D, {axisPhiSingle});
    registry.add("hPhiTrkNeg", "#varphi of negative muons;;#counts", kTH1D, {axisPhiSingle});
    registry.add("hSameSign", "hSameSign;;#counts", kTH1D, {{6, -0.5, 5.5}});
    registry.add("hPhiCharge", "#phi #it{charge}", kTH1D, {axisPhi});
    registry.add("hPhiAverage", "#phi #it{average}", kTH1D, {axisPhi});

    // ZDC
    registry.add("hTimeZNA", "ZNA Times;;#counts", kTH1D, {axisTimeZN});
    registry.add("hTimeZNC", "ZNC Times;;#counts", kTH1D, {axisTimeZN});
    registry.add("hEnergyZN", "ZNA vs ZNC energy", kTH2D, {axisEnergyZNA, axisEnergyZNC});

    // neutron classes
    reg0n0n.add("hMass", "Invariant mass of muon pairs - 0n0n;;#counts", kTH1D, {axisMass});
    reg0n0n.add("hPt", "Transverse momentum of muon pairs - 0n0n;;#counts", kTH1D, {axisPt});
    reg0n0n.add("hEta", "Pseudorapidty of muon pairs - 0n0n;;#counts", kTH1D, {axisEta});
    reg0n0n.add("hRapidity", "Rapidty of muon pairs - 0n0n;;#counts", kTH1D, {axisRapidity});
    reg0n0n.add("hPtFit", "Transverse momentum of muon pairs - 0n0n;;#counts", kTH1D, {axisPtFit});

    regXn0n.add("hMass", "Invariant mass of muon pairs - Xn0n;;#counts", kTH1D, {axisMass});
    regXn0n.add("hPt", "Transverse momentum of muon pairs - Xn0n;;#counts", kTH1D, {axisPt});
    regXn0n.add("hEta", "Pseudorapidty of muon pairs - Xn0n;;#counts", kTH1D, {axisEta});
    regXn0n.add("hRapidity", "Rapidty of muon pairs - Xn0n;;#counts", kTH1D, {axisRapidity});
    regXn0n.add("hPtFit", "Transverse momentum of muon pairs - Xn0n;;#counts", kTH1D, {axisPtFit});

    regXnXn.add("hMass", "Invariant mass of muon pairs - XnXn;;#counts", kTH1D, {axisMass});
    regXnXn.add("hPt", "Transverse momentum of muon pairs - XnXn;;#counts", kTH1D, {axisPt});
    regXnXn.add("hEta", "Pseudorapidty of muon pairs - XnXn;;#counts", kTH1D, {axisEta});
    regXnXn.add("hRapidity", "Rapidty of muon pairs - XnXn;;#counts", kTH1D, {axisRapidity});
    regXnXn.add("hPtFit", "Transverse momentum of muon pairs - XnXn;;#counts", kTH1D, {axisPtFit});
  }

  // FUNCTIONS

  // retrieve particle mass from TDatabasePDG
  float particleMass(int pid)
  {
    auto mass = pdg->Mass(pid);
    return mass;
  }

  // check if track is a GlobalMuonTrack (enum 0 or 1)
  template <typename TTrack>
  bool isGlobalMuon(const TTrack& track)
  {
    return track.trackType() == ForwardTrackTypeEnum::GlobalMuonTrack ||
           track.trackType() == ForwardTrackTypeEnum::GlobalMuonTrackOtherMatch;
  }

  // check if track is a MuonStandaloneTrack (enum 3)
  template <typename TTrack>
  bool isMuonStandalone(const TTrack& track)
  {
    return track.trackType() == ForwardTrackTypeEnum::MuonStandaloneTrack;
  }

  // check if the track pair matches the selected candidate type
  template <typename TTrack>
  bool passCandidateType(const TTrack& tr1, const TTrack& tr2)
  {
    bool isGlobal1 = isGlobalMuon(tr1);
    bool isGlobal2 = isGlobalMuon(tr2);
    bool isStandalone1 = isMuonStandalone(tr1);
    bool isStandalone2 = isMuonStandalone(tr2);

    switch (static_cast<int>(candidateType)) {
      case 0: // both GlobalMuon
        return isGlobal1 && isGlobal2;
      case 1: // mixed: one GlobalMuon + one MuonStandalone
        return (isGlobal1 && isStandalone2) || (isStandalone1 && isGlobal2);
      case 2: // both MuonStandalone
        return isStandalone1 && isStandalone2;
      default:
        return false;
    }
  }

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

  // struct used to store the ZDC info
  struct ZDCinfo {
    float timeA;
    float timeC;
    float enA;
    float enC;
  };

  // collect ZDC info per candidate
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
        zdcPerCand[candId].timeA = -999;
      if (std::isinf(zdcPerCand[candId].timeC))
        zdcPerCand[candId].timeC = -999;
      if (std::isinf(zdcPerCand[candId].enA))
        zdcPerCand[candId].enA = -999;
      if (std::isinf(zdcPerCand[candId].enC))
        zdcPerCand[candId].enC = -999;
    }
  }

  // muon track selection
  template <typename TTrack>
  bool isMuonSelected(const TTrack& fwdTrack)
  {
    float rAbs = fwdTrack.rAtAbsorberEnd();
    float pDca = fwdTrack.pDca();
    TLorentzVector p;
    auto mMu = o2::constants::physics::MassMuon;
    p.SetXYZM(fwdTrack.px(), fwdTrack.py(), fwdTrack.pz(), mMu);
    float eta = p.Eta();
    float pt = p.Pt();
    float pDcaMax = rAbs < kRAbsMid ? kPDca1 : kPDca2;

    if (eta < kEtaMin || eta > kEtaMax)
      return false;
    if (pt < kPtMin)
      return false;
    if (rAbs < kRAbsMin || rAbs > kRAbsMax)
      return false;
    if (pDca > pDcaMax)
      return false;
    return true;
  }

  // compute phi for azimuth anisotropy
  void computePhiAnis(TLorentzVector p1, TLorentzVector p2, int sign1, float& phiAverage, float& phiCharge)
  {
    TLorentzVector tSum, tDiffAv, tDiffCh;
    tSum = p1 + p2;
    float halfUnity = 0.5;
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

    phiAverage = tSum.DeltaPhi(tDiffAv);
    phiCharge = tSum.DeltaPhi(tDiffCh);
  }

  // process a single candidate
  void processCand(CandidatesFwd::iterator const& cand,
                   ForwardTracks::iterator const& tr1, ForwardTracks::iterator const& tr2,
                   ZDCinfo const& zdc)
  {
    // candidate type selection
    if (!passCandidateType(tr1, tr2))
      return;

    // V0 selection
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 1) {
        if (ampsV0A[i] > kMaxAmpV0A)
          return;
      }
    }

    // select opposite charge events only
    if (cand.netCharge() != 0) {
      registry.fill(HIST("hSameSign"), cand.numContrib());
      return;
    }

    // MCH-MID match selection: both tracks must have MCH-MID match
    int nMIDs = 0;
    if (tr1.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (tr2.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (nMIDs != 2)
      return;

    // track selection
    if (!isMuonSelected(tr1))
      return;
    if (!isMuonSelected(tr2))
      return;

    // form Lorentz vectors
    TLorentzVector p1, p2;
    auto mMu = o2::constants::physics::MassMuon;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
    TLorentzVector p = p1 + p2;

    // cut on pair kinematics
    if (p.M() < lowMass || p.M() > highMass)
      return;
    if (p.Pt() < lowPt || p.Pt() > highPt)
      return;
    if (p.Rapidity() < lowRapidity || p.Rapidity() > highRapidity)
      return;

    // compute phi for azimuth anisotropy
    float phiAverage = 0;
    float phiCharge = 0;
    computePhiAnis(p1, p2, tr1.sign(), phiAverage, phiCharge);

    // zdc info
    if (std::abs(zdc.timeA) < kMaxZDCTimeHisto)
      registry.fill(HIST("hTimeZNA"), zdc.timeA);
    if (std::abs(zdc.timeC) < kMaxZDCTimeHisto)
      registry.fill(HIST("hTimeZNC"), zdc.timeC);
    registry.fill(HIST("hEnergyZN"), zdc.enA, zdc.enC);

    // neutron classes
    bool neutronA = std::abs(zdc.timeA) < kMaxZDCTime && !std::isinf(zdc.timeA);
    bool neutronC = std::abs(zdc.timeC) < kMaxZDCTime && !std::isinf(zdc.timeC);

    if (!neutronC && !neutronA) {
      reg0n0n.fill(HIST("hMass"), p.M());
      reg0n0n.fill(HIST("hPt"), p.Pt());
      reg0n0n.fill(HIST("hPtFit"), p.Pt());
      reg0n0n.fill(HIST("hEta"), p.Eta());
      reg0n0n.fill(HIST("hRapidity"), p.Rapidity());
    } else if (neutronA ^ neutronC) {
      regXn0n.fill(HIST("hMass"), p.M());
      regXn0n.fill(HIST("hPt"), p.Pt());
      regXn0n.fill(HIST("hPtFit"), p.Pt());
      regXn0n.fill(HIST("hEta"), p.Eta());
      regXn0n.fill(HIST("hRapidity"), p.Rapidity());
    } else if (neutronA && neutronC) {
      regXnXn.fill(HIST("hMass"), p.M());
      regXnXn.fill(HIST("hPt"), p.Pt());
      regXnXn.fill(HIST("hPtFit"), p.Pt());
      regXnXn.fill(HIST("hEta"), p.Eta());
      regXnXn.fill(HIST("hRapidity"), p.Rapidity());
    }

    // fill inclusive histos
    registry.fill(HIST("hContrib"), cand.numContrib());
    registry.fill(HIST("hPtTrkPos"), p1.Pt());
    registry.fill(HIST("hPtTrkNeg"), p2.Pt());
    registry.fill(HIST("hEtaTrkPos"), p1.Eta());
    registry.fill(HIST("hEtaTrkNeg"), p2.Eta());
    registry.fill(HIST("hPhiTrkPos"), p1.Phi());
    registry.fill(HIST("hPhiTrkNeg"), p2.Phi());
    registry.fill(HIST("hEvSign"), cand.netCharge());
    registry.fill(HIST("hMass"), p.M());
    registry.fill(HIST("hPt"), p.Pt());
    registry.fill(HIST("hPtFit"), p.Pt());
    registry.fill(HIST("hPtFit2"), p.Pt());
    registry.fill(HIST("hEta"), p.Eta());
    registry.fill(HIST("hRapidity"), p.Rapidity());
    registry.fill(HIST("hPhi"), p.Phi());
    registry.fill(HIST("hCharge"), tr1.sign());
    registry.fill(HIST("hCharge"), tr2.sign());
    registry.fill(HIST("hPhiAverage"), phiAverage);
    registry.fill(HIST("hPhiCharge"), phiCharge);
  }

  // PROCESS FUNCTION
  void processData(CandidatesFwd const& eventCandidates,
                   o2::aod::UDZdcsReduced const& ZDCs,
                   ForwardTracks const& fwdTracks)
  {
    // map with the tracks
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    // map with the ZDC info
    std::unordered_map<int32_t, ZDCinfo> zdcPerCand;
    collectCandZDCInfo(zdcPerCand, ZDCs);

    // loop over the candidates
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      auto cand = eventCandidates.iteratorAt(candID);
      auto tr1 = fwdTracks.iteratorAt(trId1);
      auto tr2 = fwdTracks.iteratorAt(trId2);

      ZDCinfo zdc;
      if (zdcPerCand.count(candID) != 0) {
        zdc = zdcPerCand.at(candID);
      } else {
        zdc.timeA = -999;
        zdc.timeC = -999;
        zdc.enA = -999;
        zdc.enC = -999;
      }

      processCand(cand, tr1, tr2, zdc);
    }
  }

  PROCESS_SWITCH(UpcFwdJpsiRL, processData, "", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcFwdJpsiRL>(cfgc),
  };
}
