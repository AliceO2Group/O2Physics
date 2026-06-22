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

/// \file fwdMuonsUpc.cxx
/// \brief perform some selections on fwd events and saves the results

/// executable name o2-analysis-ud-fwd-muon-upc

/// \author Andrea Giovanni Riffero <andrea.giovanni.riffero@cern.ch>

#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
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

#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <vector>

// table for saving tree with info on data
namespace dimu
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
// tracks positive (p) and negative (n)
DECLARE_SOA_COLUMN(EnergyP, energyP, float);
DECLARE_SOA_COLUMN(Pxp, pxp, float);
DECLARE_SOA_COLUMN(Pyp, pyp, float);
DECLARE_SOA_COLUMN(Pzp, pzp, float);
DECLARE_SOA_COLUMN(Ptp, ptp, float);
DECLARE_SOA_COLUMN(Etap, etap, float);
DECLARE_SOA_COLUMN(Phip, phip, float);
DECLARE_SOA_COLUMN(TrackTypep, trackTypep, int);
DECLARE_SOA_COLUMN(EnergyN, energyN, float);
DECLARE_SOA_COLUMN(Pxn, pxn, float);
DECLARE_SOA_COLUMN(Pyn, pyn, float);
DECLARE_SOA_COLUMN(Pzn, pzn, float);
DECLARE_SOA_COLUMN(Ptn, ptn, float);
DECLARE_SOA_COLUMN(Etan, etan, float);
DECLARE_SOA_COLUMN(Phin, phin, float);
DECLARE_SOA_COLUMN(TrackTypen, trackTypen, int);
// zn
DECLARE_SOA_COLUMN(Tzna, tzna, float);
DECLARE_SOA_COLUMN(Ezna, ezna, float);
DECLARE_SOA_COLUMN(Tznc, tznc, float);
DECLARE_SOA_COLUMN(Eznc, eznc, float);
DECLARE_SOA_COLUMN(Nclass, nclass, int);
// other info
DECLARE_SOA_COLUMN(PDCAp, pDCAp, float);
DECLARE_SOA_COLUMN(PDCAn, pDCAn, float);
DECLARE_SOA_COLUMN(AmpV0A, ampV0A, std::vector<float>);
} // namespace dimu

namespace o2::aod
{
DECLARE_SOA_TABLE(DiMu, "AOD", "DIMU",
                  dimu::RunNumber,
                  dimu::M, dimu::Energy, dimu::Px, dimu::Py, dimu::Pz, dimu::Pt, dimu::Rap, dimu::Phi,
                  dimu::PhiAv, dimu::PhiCh,
                  dimu::EnergyP, dimu::Pxp, dimu::Pyp, dimu::Pzp, dimu::Ptp, dimu::Etap, dimu::Phip, dimu::TrackTypep,
                  dimu::EnergyN, dimu::Pxn, dimu::Pyn, dimu::Pzn, dimu::Ptn, dimu::Etan, dimu::Phin, dimu::TrackTypen,
                  dimu::Tzna, dimu::Ezna, dimu::Tznc, dimu::Eznc, dimu::Nclass,
                  dimu::PDCAp, dimu::PDCAn, dimu::AmpV0A);
} // namespace o2::aod

// for saving tree with info on gen MC
namespace gendimu
{
// dimuon
DECLARE_SOA_COLUMN(GenM, genM, float);
DECLARE_SOA_COLUMN(GenPt, genPt, float);
DECLARE_SOA_COLUMN(GenRap, genRap, float);
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);
DECLARE_SOA_COLUMN(GenPhiAv, genPhiAv, float);
DECLARE_SOA_COLUMN(GenPhiCh, genPhiCh, float);
// tracks positive (p) and negative (n)
DECLARE_SOA_COLUMN(GenPtp, genPtp, float);
DECLARE_SOA_COLUMN(GenEtap, genEtap, float);
DECLARE_SOA_COLUMN(GenPhip, genPhip, float);
DECLARE_SOA_COLUMN(GenPtn, genPtn, float);
DECLARE_SOA_COLUMN(GenEtan, genEtan, float);
DECLARE_SOA_COLUMN(GenPhin, genPhin, float);
} // namespace gendimu

namespace o2::aod
{
DECLARE_SOA_TABLE(GenDimu, "AOD", "GENDIMU",
                  gendimu::GenM, gendimu::GenPt, gendimu::GenRap, gendimu::GenPhi,
                  gendimu::GenPhiAv, gendimu::GenPhiCh,
                  gendimu::GenPtp, gendimu::GenEtap, gendimu::GenPhip,
                  gendimu::GenPtn, gendimu::GenEtan, gendimu::GenPhin);
} // namespace o2::aod

// for saving tree with info on reco MC
namespace recodimu
{
// dimuon
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(PhiAv, phiAv, float);
DECLARE_SOA_COLUMN(PhiCh, phiCh, float);
// tracks positive (p) and negative (n)
DECLARE_SOA_COLUMN(Ptp, ptp, float);
DECLARE_SOA_COLUMN(Etap, etap, float);
DECLARE_SOA_COLUMN(Phip, phip, float);
DECLARE_SOA_COLUMN(TrackTypep, trackTypep, int);
DECLARE_SOA_COLUMN(Ptn, ptn, float);
DECLARE_SOA_COLUMN(Etan, etan, float);
DECLARE_SOA_COLUMN(Phin, phin, float);
DECLARE_SOA_COLUMN(TrackTypen, trackTypen, int);
// gen info dimuon
DECLARE_SOA_COLUMN(GenPt, genPt, float);
DECLARE_SOA_COLUMN(GenRap, genRap, float);
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);
// gen info trks
DECLARE_SOA_COLUMN(GenPtp, genPtp, float);
DECLARE_SOA_COLUMN(GenEtap, genEtap, float);
DECLARE_SOA_COLUMN(GenPhip, genPhip, float);
DECLARE_SOA_COLUMN(GenPtn, genPtn, float);
DECLARE_SOA_COLUMN(GenEtan, genEtan, float);
DECLARE_SOA_COLUMN(GenPhin, genPhin, float);
// other info
DECLARE_SOA_COLUMN(PDCAp, pDCAp, float);
DECLARE_SOA_COLUMN(PDCAn, pDCAn, float);
DECLARE_SOA_COLUMN(AmpV0A, ampV0A, std::vector<float>);
} // namespace recodimu

namespace o2::aod
{
DECLARE_SOA_TABLE(RecoDimu, "AOD", "RECODIMU",
                  recodimu::RunNumber,
                  recodimu::M, recodimu::Pt, recodimu::Rap, recodimu::Phi,
                  recodimu::PhiAv, recodimu::PhiCh,
                  recodimu::Ptp, recodimu::Etap, recodimu::Phip, recodimu::TrackTypep,
                  recodimu::Ptn, recodimu::Etan, recodimu::Phin, recodimu::TrackTypen,
                  recodimu::GenPt, recodimu::GenRap, recodimu::GenPhi,
                  recodimu::GenPtp, recodimu::GenEtap, recodimu::GenPhip,
                  recodimu::GenPtn, recodimu::GenEtan, recodimu::GenPhin,
                  recodimu::PDCAp, recodimu::PDCAn, recodimu::AmpV0A);
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// constants used in the track selection
const float kRAbsMin = 17.6;
const float kRAbsMid = 26.5;
const float kRAbsMax = 89.5;
float kEtaMin = -4.0;
float kEtaMax = -2.5;
const float kPtMin = 0.;

const int kReqMatchMIDTracks = 2;
const int kReqMatchMFTTracks = 2;
const int kMaxChi2MFTMatch = 30;
const float kMaxZDCTime = 2.;
const int k2Tracks = 2;
const int k4Tracks = 4;

struct FwdMuonsUpc {

  // a PDG object
  Service<o2::framework::O2DatabasePDG> pdg;

  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using CompleteFwdTracks = soa::Join<ForwardTracks, o2::aod::UDMcFwdTrackLabels>;

  Produces<o2::aod::DiMu> dimuSel;
  Produces<o2::aod::GenDimu> dimuGen;
  Produces<o2::aod::RecoDimu> dimuReco;

  // defining histograms using histogram registry: different histos for the different process functions
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // CONFIGURABLES
  // pT of muon pairs
  ConfigurableAxis axisPt{"axisPt", {250, 0.0f, 2.0f}, "#it{p}_{T}^{#mu#mu} (GeV/#it{c})"};
  // pT fit
  ConfigurableAxis axisPtFit{"axisPtFit", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.50, 3.00, 3.50}, "#it{p}_{T} (GeV/c)"};
  // mass of muon pairs
  ConfigurableAxis axisMass{"axisMass", {500, 0.0f, 10.0f}, "m_{#mu#mu} (GeV/#it{c}^{2})"};
  // rapidity of muon pairs
  ConfigurableAxis axisRapidity{"axisRapidity", {250, -4.5f, -2.0f}, "y_{#mu#mu}"};

  // cuts on pair kinematics
  Configurable<float> lowPt{"lowPt", 0., "Low pT cut"};
  Configurable<float> highPt{"highPt", 2, "High pT cut"};
  Configurable<float> lowMass{"lowMass", 0., "Low mass cut"};
  Configurable<float> highMass{"highMass", 10., "High mass cut"};
  Configurable<float> lowRapidity{"lowRapidity", -4.5, "Low rapidity cut"};
  Configurable<float> highRapidity{"highRapidity", -2., "High rapidity cut"};

  // cuts on pDCA
  Configurable<float> pDCA1{"pDCA1", 200., "pDCA cut for tracks with rAbs < 26.5 cm"};
  Configurable<float> pDCA2{"pDCA2", 200., "pDCA cut for tracks with rAbs > 26.5 cm"};

  // cut on V0A amplitude
  Configurable<float> maxAmpV0A{"maxAmpV0A", 100., "Max amplitude in V0A"};

  // my track type
  // 0 = MCH-MID-MFT
  // 1 = MCH-MID
  Configurable<int> myTrackType{"myTrackType", 3, "My track type"};

  void init(InitContext&)
  {
    // axis
    const AxisSpec ptAxis{axisPt, "#it{p}_{T}^{#mu#mu} (GeV/#it{c})", "ptAxis"};
    const AxisSpec ptFitAxis{axisPtFit, "#it{p}_{T} (GeV/c)", "ptFitAxis"};
    const AxisSpec massAxis{axisMass, "m_{#mu#mu} (GeV/#it{c}^{2})", "massAxis"};
    const AxisSpec rapidityAxis{axisRapidity, "y_{#mu#mu}", "rapidityAxis"};

    // histos
    registry.add("hMass", "Invariant mass of muon pairs;;#counts", kTH1D, {massAxis});
    registry.add("hPt", "Transverse momentum of muon pairs;;#counts", kTH1D, {ptAxis});
    registry.add("hPtFit", "Transverse momentum of muon pairs;;#counts", kTH1D, {ptFitAxis});
    registry.add("hRapidity", "Rapidty of muon pairs;;#counts", kTH1D, {rapidityAxis});
  }

  // FUNCTIONS

  // template function that fills a map with the collision id of each udcollision as key
  // and a vector with the tracks
  // map == (key, element) == (udCollisionId, vector of trks)
  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks const& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  // template function that fills a map with the collision id of each udmccollision as key
  // and a vector with the tracks
  // map == (key, element) == (udMcCollisionId, vector of mc particles)
  template <typename TTracks>
  void collectMcCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks const& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udMcCollisionId();
      if (candId < 0) {
        continue;
      }
      if (std::abs(tr.pdgCode()) != PDG_t::kMuonMinus) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  // template function that fills a map with the collision id of each udcollision as key
  // and a vector with the tracks and corresponding geneated particles
  // map == (key, element) == (udCollisionId, vector(track1, mcPart1, track2, mcPart2))
  template <typename TTracks>
  void collectRecoCandID(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId < 0)
        continue;

      if (!tr.has_udMcParticle()) {
        // LOGF(debug,"tr does not have mc part");
        continue;
      }
      // retrieve mc particle from the reco track
      auto mcPart = tr.udMcParticle();

      tracksPerCand[candId].push_back(tr.globalIndex());
      tracksPerCand[candId].push_back(mcPart.globalIndex());
    }
  }

  // struct used to store the ZDC info in a map
  struct ZDCinfo {
    float timeA = 0.f;
    float timeC = 0.f;
    float enA = 0.f;
    float enC = 0.f;
    int32_t id = -1;
  };

  // function that fills a map with the collision id of each udcollision as key
  // and a ZDCinfo struct with the ZDC information
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

      // take care of the infinity
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

  // function to select muon tracks
  template <typename TTracks>
  bool isMuonSelected(const TTracks& fwdTrack)
  {
    float rAbs = fwdTrack.rAtAbsorberEnd();
    float pDca = fwdTrack.pDca();

    std::array<float, 3> trackMomentum{fwdTrack.px(), fwdTrack.py(), fwdTrack.pz()};
    float eta = RecoDecay::eta(trackMomentum);
    float pt = RecoDecay::pt(trackMomentum);
    float pDcaMax = rAbs < kRAbsMid ? pDCA1 : pDCA2;

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

  // function to compute phi for azimuth anisotropy
  void computePhiAnis(ROOT::Math::PxPyPzMVector p1, ROOT::Math::PxPyPzMVector p2, int sign1, float& phiAverage, float& phiCharge)
  {
    ROOT::Math::PxPyPzMVector tSum, tDiffAv, tDiffCh;
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

    // average
    phiAverage = ROOT::Math::VectorUtil::DeltaPhi(tSum, tDiffAv);
    // charge
    phiCharge = ROOT::Math::VectorUtil::DeltaPhi(tSum, tDiffCh);
  }

  // function that processes the candidates:
  // it applies V0 selection, trk selection, kine selection, and fills the histograms
  // it also divides the data in neutron classes
  // used for real data
  void processCand(CandidatesFwd::iterator const& cand,
                   ForwardTracks::iterator const& tr1, ForwardTracks::iterator const& tr2,
                   ZDCinfo const& zdc)
  {
    // V0 selection
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 0) {
        if (ampsV0A[i] > maxAmpV0A)
          return;
      }
    }

    // select events with exactly 2 forward tracks
    if (cand.numContrib() != k2Tracks) {
      return;
    }

    // select opposite charge events only
    if (cand.netCharge() != 0) {
      return;
    }

    // MCH-MID match selection
    int nMIDs = 0;
    if (tr1.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (tr2.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (nMIDs != kReqMatchMIDTracks)
      return;

    // MFT-MID match selection (if MFT is requested by the trackType)
    if (myTrackType == 0) {
      // if MFT is requested check that the tracks is inside the MFT acceptance
      kEtaMin = -3.6;
      kEtaMax = -2.5;

      int nMFT = 0;
      if (tr1.chi2MatchMCHMFT() > 0 && tr1.chi2MatchMCHMFT() < kMaxChi2MFTMatch)
        nMFT++;
      if (tr2.chi2MatchMCHMFT() > 0 && tr2.chi2MatchMCHMFT() < kMaxChi2MFTMatch)
        nMFT++;
      if (nMFT != kReqMatchMFTTracks)
        return;
    }

    // track selection
    if (!isMuonSelected(*tr1))
      return;
    if (!isMuonSelected(*tr2))
      return;

    // form Lorentz vectors
    ROOT::Math::PxPyPzMVector p1{tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassMuon};
    ROOT::Math::PxPyPzMVector p2{tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassMuon};
    ROOT::Math::PxPyPzMVector p = p1 + p2;

    // cut on pair kinematics
    // select mass
    if (p.M() < lowMass)
      return;
    if (p.M() > highMass)
      return;
    // select pt
    if (p.Pt() < lowPt)
      return;
    if (p.Pt() > highPt)
      return;
    // select rapidity
    if (p.Rapidity() < lowRapidity)
      return;
    if (p.Rapidity() > highRapidity)
      return;

    // compute phi for azimuth anisotropy
    float phiAverage = 0;
    float phiCharge = 0;
    computePhiAnis(p1, p2, tr1.sign(), phiAverage, phiCharge);

    // divide the events in neutron classes
    bool neutronA = false;
    bool neutronC = false;
    int znClass = -1;

    if (std::abs(zdc.timeA) < kMaxZDCTime)
      neutronA = true;
    if (std::abs(zdc.timeC) < kMaxZDCTime)
      neutronC = true;

    // assign neutron class label
    // 0n0n
    if (neutronC == false && neutronA == false) {
      znClass = 1;
    } else if (neutronA ^ neutronC) { // Xn0n + 0nXn
      if (neutronA)
        znClass = 2;
      else if (neutronC)
        znClass = 3;
    } else if (neutronA && neutronC) { // XnXn
      znClass = 4;
    }

    // fill the histos
    registry.fill(HIST("hMass"), p.M());
    registry.fill(HIST("hPt"), p.Pt());
    registry.fill(HIST("hPtFit"), p.Pt());
    registry.fill(HIST("hRapidity"), p.Rapidity());

    // store the event to save it into a tree
    if (tr1.sign() > 0) {
      dimuSel(cand.runNumber(),
              p.M(), p.E(), p.Px(), p.Py(), p.Pz(), p.Pt(), p.Rapidity(), p.Phi(),
              phiAverage, phiCharge,
              p1.E(), p1.Px(), p1.Py(), p1.Pz(), p1.Pt(), p1.Eta(), p1.Phi(), static_cast<int>(myTrackType),
              p2.E(), p2.Px(), p2.Py(), p2.Pz(), p2.Pt(), p2.Eta(), p2.Phi(), static_cast<int>(myTrackType),
              zdc.timeA, zdc.enA, zdc.timeC, zdc.enC, znClass,
              tr1.pDca(), tr2.pDca(), ampsV0A);
    } else {
      dimuSel(cand.runNumber(),
              p.M(), p.E(), p.Px(), p.Py(), p.Pz(), p.Pt(), p.Rapidity(), p.Phi(),
              phiAverage, phiCharge,
              p2.E(), p2.Px(), p2.Py(), p2.Pz(), p2.Pt(), p2.Eta(), p2.Phi(), static_cast<int>(myTrackType),
              p1.E(), p1.Px(), p1.Py(), p1.Pz(), p1.Pt(), p1.Eta(), p1.Phi(), static_cast<int>(myTrackType),
              zdc.timeA, zdc.enA, zdc.timeC, zdc.enC, znClass,
              tr2.pDca(), tr1.pDca(), ampsV0A);
    }
  }

  // function that processes the MC gen candidates:
  // it applies some kinematics cut and fills the histograms
  void processMcGenCand(aod::UDMcCollisions::iterator const& /*mcCand*/,
                        aod::UDMcParticles::iterator const& McPart1, aod::UDMcParticles::iterator const& McPart2)
  {

    // check that all pairs are mu+mu-
    if (std::abs(McPart1.pdgCode()) != PDG_t::kMuonMinus || std::abs(McPart2.pdgCode()) != PDG_t::kMuonMinus) {
      LOGF(debug, "PDG codes: %d | %d", McPart1.pdgCode(), McPart2.pdgCode());
      return;
    }
    if (McPart1.pdgCode() + McPart2.pdgCode() != 0) {
      return;
    }

    // create Lorentz vectors
    ROOT::Math::PxPyPzMVector p1{McPart1.px(), McPart1.py(), McPart1.pz(), o2::constants::physics::MassMuon};
    ROOT::Math::PxPyPzMVector p2{McPart2.px(), McPart2.py(), McPart2.pz(), o2::constants::physics::MassMuon};
    ROOT::Math::PxPyPzMVector p = p1 + p2;

    // cut on pair kinematics
    // select mass
    if (p.M() < lowMass)
      return;
    if (p.M() > highMass)
      return;
    // select pt
    if (p.Pt() < lowPt)
      return;
    if (p.Pt() > highPt)
      return;
    // select rapidity
    if (p.Rapidity() < lowRapidity)
      return;
    if (p.Rapidity() > highRapidity)
      return;

    // compute phi for azimuth anisotropy
    float phiAverage = 0;
    float phiCharge = 0;
    computePhiAnis(p1, p2, -McPart1.pdgCode(), phiAverage, phiCharge);

    // fill the histos
    registry.fill(HIST("hMass"), p.M());
    registry.fill(HIST("hPt"), p.Pt());
    registry.fill(HIST("hPtFit"), p.Pt());
    registry.fill(HIST("hRapidity"), p.Rapidity());

    // store the event to save it into a tree
    if (McPart1.pdgCode() < 0) {
      dimuGen(p.M(), p.Pt(), p.Rapidity(), p.Phi(),
              phiAverage, phiCharge,
              p1.Pt(), p1.Eta(), p1.Phi(),
              p2.Pt(), p2.Eta(), p2.Phi());
    } else {
      dimuGen(p.M(), p.Pt(), p.Rapidity(), p.Phi(),
              phiAverage, phiCharge,
              p2.Pt(), p2.Eta(), p2.Phi(),
              p1.Pt(), p1.Eta(), p1.Phi());
    }
  }

  // function that processes MC reco candidates
  // it applies V0 selection, trk selection, kine selection, and fills the histograms
  void processMcRecoCand(CandidatesFwd::iterator const& cand,
                         CompleteFwdTracks::iterator const& tr1, aod::UDMcParticles::iterator const& McPart1,
                         CompleteFwdTracks::iterator const& tr2, aod::UDMcParticles::iterator const& McPart2)
  {

    // check that all pairs are mu+mu-
    if (std::abs(McPart1.pdgCode()) != PDG_t::kMuonMinus || std::abs(McPart2.pdgCode()) != PDG_t::kMuonMinus) {
      LOGF(debug, "PDG codes: %d | %d", McPart1.pdgCode(), McPart2.pdgCode());
      return;
    }

    // V0 selection
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 0) {
        if (ampsV0A[i] > maxAmpV0A)
          return;
      }
    }

    // select events with exactly 2 forward tracks
    if (cand.numContrib() != k2Tracks) {
      return;
    }

    // select opposite charge events only
    if (cand.netCharge() != 0) {
      return;
    }

    // MCH-MID match selection
    int nMIDs = 0;
    if (tr1.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (tr2.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (nMIDs != kReqMatchMIDTracks)
      return;

    // MFT-MID match selection (if MFT is requested by the trackType)
    if (myTrackType == 0) {
      // if MFT is requested check that the tracks is inside the MFT acceptance
      kEtaMin = -3.6;
      kEtaMax = -2.5;

      int nMFT = 0;
      if (tr1.chi2MatchMCHMFT() > 0 && tr1.chi2MatchMCHMFT() < kMaxChi2MFTMatch)
        nMFT++;
      if (tr2.chi2MatchMCHMFT() > 0 && tr2.chi2MatchMCHMFT() < kMaxChi2MFTMatch)
        nMFT++;
      if (nMFT != kReqMatchMFTTracks)
        return;
    }

    // track selection
    if (!isMuonSelected(*tr1))
      return;
    if (!isMuonSelected(*tr2))
      return;

    // form Lorentz vectors
    ROOT::Math::PxPyPzMVector p1{tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassMuon};
    ROOT::Math::PxPyPzMVector p2{tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassMuon};
    ROOT::Math::PxPyPzMVector p = p1 + p2;

    // cut on pair kinematics (reco candidates)
    // select mass
    if (p.M() < lowMass)
      return;
    if (p.M() > highMass)
      return;
    // select pt
    if (p.Pt() < lowPt)
      return;
    if (p.Pt() > highPt)
      return;
    // select rapidity
    if (p.Rapidity() < lowRapidity)
      return;
    if (p.Rapidity() > highRapidity)
      return;

    // compute phi for azimuth anisotropy
    float phiAverage = 0;
    float phiCharge = 0;
    computePhiAnis(p1, p2, tr1.sign(), phiAverage, phiCharge);

    ROOT::Math::PxPyPzMVector p1Mc{McPart1.px(), McPart1.py(), McPart1.pz(), o2::constants::physics::MassMuon};
    ROOT::Math::PxPyPzMVector p2Mc{McPart2.px(), McPart2.py(), McPart2.pz(), o2::constants::physics::MassMuon};
    ROOT::Math::PxPyPzMVector pMc = p1Mc + p2Mc;

    // compute gen phi for azimuth anisotropy
    float phiGenAverage = 0;
    float phiGenCharge = 0;
    computePhiAnis(p1Mc, p2Mc, -McPart1.pdgCode(), phiGenAverage, phiGenCharge);

    // print info in case of problems
    if (tr1.sign() * McPart1.pdgCode() > 0 || tr2.sign() * McPart2.pdgCode() > 0) {
      LOGF(debug, "Problem: ");
      LOGF(debug, "real: %d | %d", (int)tr1.sign(), (int)tr2.sign());
      LOGF(debug, "mc  : %i | %i", (int)McPart1.pdgCode(), (int)McPart2.pdgCode());
      LOGF(debug, "contrib: %d", (int)cand.numContrib());
    }

    // fill the histos
    registry.fill(HIST("hMass"), p.M());
    registry.fill(HIST("hPt"), p.Pt());
    registry.fill(HIST("hPtFit"), p.Pt());
    registry.fill(HIST("hRapidity"), p.Rapidity());

    // store the event to save it into a tree
    if (tr1.sign() > 0) {
      dimuReco(cand.runNumber(),
               p.M(), p.Pt(), p.Rapidity(), p.Phi(),
               phiAverage, phiCharge,
               p1.Pt(), p1.Eta(), p1.Phi(), static_cast<int>(myTrackType),
               p2.Pt(), p2.Eta(), p2.Phi(), static_cast<int>(myTrackType),
               // gen info
               pMc.Pt(), pMc.Rapidity(), pMc.Phi(),
               p1Mc.Pt(), p1Mc.Eta(), p1Mc.Phi(),
               p2Mc.Pt(), p2Mc.Eta(), p2Mc.Phi(),
               tr1.pDca(), tr2.pDca(), ampsV0A);
    } else {
      dimuReco(cand.runNumber(),
               p.M(), p.Pt(), p.Rapidity(), p.Phi(),
               phiAverage, phiCharge,
               p2.Pt(), p2.Eta(), p2.Phi(), static_cast<int>(myTrackType),
               p1.Pt(), p1.Eta(), p1.Phi(), static_cast<int>(myTrackType),
               // gen info
               pMc.Pt(), pMc.Rapidity(), pMc.Phi(),
               p2Mc.Pt(), p2Mc.Eta(), p2Mc.Phi(),
               p1Mc.Pt(), p1Mc.Eta(), p1Mc.Phi(),
               tr2.pDca(), tr1.pDca(), ampsV0A);
    }
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
      if (item.second.size() != k2Tracks) {
        LOGF(debug, "number track = %d", item.second.size());
        continue;
      }
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

  PROCESS_SWITCH(FwdMuonsUpc, processData, "", true);

  // process MC Truth
  void processMcGen(aod::UDMcCollisions const& mccollisions, aod::UDMcParticles const& McParts)
  {
    // map with the tracks
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectMcCandIDs(tracksPerCand, McParts);

    // loop over the candidates
    for (const auto& item : tracksPerCand) {
      if (item.second.size() != k2Tracks) {
        LOGF(debug, "mc parts = %d", item.second.size());
        for (const auto& id : item.second) {
          auto p = McParts.iteratorAt(id);
          LOGF(debug,
               "  part %d: pdg=%d status=%d has_mothers=%d has_daughters=%d",
               id, p.pdgCode(), p.statusCode(),
               p.has_mothers(), p.has_daughters());
        }
        continue;
      }
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      auto cand = mccollisions.iteratorAt(candID);
      auto tr1 = McParts.iteratorAt(trId1);
      auto tr2 = McParts.iteratorAt(trId2);

      processMcGenCand(cand, tr1, tr2);
    }
  }
  PROCESS_SWITCH(FwdMuonsUpc, processMcGen, "", false);

  // process reco MC (gen info included)
  void processMcReco(CandidatesFwd const& eventCandidates,
                     CompleteFwdTracks const& fwdTracks,
                     aod::UDMcCollisions const& /*mccollisions*/,
                     aod::UDMcParticles const& McParts)
  {
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCandAll;
    collectRecoCandID(tracksPerCandAll, fwdTracks);

    // loop over the candidates
    for (const auto& item : tracksPerCandAll) {
      if (item.second.size() != k4Tracks) {
        LOGF(debug, "number track (reco + gen) = %d", item.second.size());
        continue;
      }

      // get the reco candidate
      auto cand = eventCandidates.iteratorAt(item.first);

      // get reco tracks and corresponding gen particles
      auto tr1 = fwdTracks.iteratorAt(item.second[0]);
      auto trMc1 = McParts.iteratorAt(item.second[1]);
      auto tr2 = fwdTracks.iteratorAt(item.second[2]);
      auto trMc2 = McParts.iteratorAt(item.second[3]);

      // check that the method used here gets the the MC particles
      // as the one used by Nazar
      auto nzTrMc1 = McParts.iteratorAt(tr1.udMcParticleId());
      auto nzTrMc2 = McParts.iteratorAt(tr2.udMcParticleId());

      if (nzTrMc1 != trMc1)
        LOGF(debug, "diff wrt Nazar!");
      if (nzTrMc2 != trMc2)
        LOGF(debug, "diff wrt Nazar!");

      processMcRecoCand(cand, tr1, trMc1, tr2, trMc2);
    }
  }
  PROCESS_SWITCH(FwdMuonsUpc, processMcReco, "", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FwdMuonsUpc>(cfgc),
  };
}
