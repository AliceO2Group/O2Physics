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

/// \file   lambda1405analysis.cxx
/// \brief Analysis task for lambda1405 via sigma kink decay
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>

#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFFullPi, aod::pidTOFFullPr>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

struct lambda1405candidate {
  // Columns for Lambda(1405) candidate
  bool isSigmaPlus = false;  // True if compatible with Sigma+
  bool isSigmaMinus = false; // True if compatible with Sigma-
  float mass = -1;           // Invariant mass of the Lambda(1405) candidate
  float sigmaMinusMass = -1; // Invariant mass of the Sigma- candidate
  float sigmaPlusMass = -1;  // Invariant mass of the Sigma+ candidate
  float pt = -1;             // pT of the Lambda(1405) candidate
  int sigmaSign = 0;         // Sign of the Sigma candidate: 1 for matter, -1 for antimatter
  float sigmaPt = -1;        // pT of the Sigma daughter
  float piPt = -1;           // pT of the pion daughter
  float nSigmaTPCPi = -1;    // Number of sigmas for the pion candidate
  float nSigmaTOFPi = -1;    // Number of sigmas for the pion candidate using TOF
  int kinkDauID = 0;         // ID of the pion from Sigma decay in MC
  int sigmaID = 0;           // ID of the Sigma candidate in MC
  int piID = 0;              // ID of the pion candidate in MC
};

struct lambda1405analysis {
  int lambda1405PdgCode = 102132;     // PDG code for Lambda(1405)
  lambda1405candidate lambda1405Cand; // Lambda(1405) candidate structure
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda1405{"lambda1405", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // Configurable for event selection
  Configurable<float> cutzvertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutEtaDaught{"cutEtaDaughter", 0.8f, "Eta cut for daughter tracks"};
  Configurable<float> cutDCAtoPVSigma{"cutDCAtoPVSigma", 0.1f, "Max DCA to primary vertex for Sigma candidates (cm)"};
  Configurable<float> cutDCAtoPVPiFromSigma{"cutDCAtoPVPiFromSigma", 2., "Min DCA to primary vertex for pion from Sigma candidates (cm)"};

  Configurable<float> cutSigmaRadius{"cutSigmaRadius", 20.f, "Minimum radius for Sigma candidates (cm)"};
  Configurable<float> cutSigmaMass{"cutSigmaMass", 0.1, "Sigma mass window (MeV/c^2)"};
  Configurable<float> cutNITSClusPi{"cutNITSClusPi", 5, "Minimum number of ITS clusters for pion candidate"};
  Configurable<float> cutNTPCClusPi{"cutNTPCClusPi", 90, "Minimum number of TPC clusters for pion candidate"};
  Configurable<float> cutNSigmaTPC{"cutNSigmaTPC", 3, "NSigmaTPCPion"};
  Configurable<float> cutNSigmaTOF{"cutNSigmaTOF", 3, "NSigmaTOFPion"};

  Configurable<bool> doLSBkg{"doLikeSignBkg", false, "Use like-sign background"};
  Configurable<bool> useTOF{"useTOF", false, "Use TOF for PID for pion candidates"};

  Preslice<aod::KinkCands> mKinkPerCol = aod::track::collisionId;
  Preslice<aod::TracksIU> mPerColTracks = aod::track::collisionId;

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec ptAxis{100, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptPiAxis{100, -5, 5, "#it{p}_{T}^{#pi} (GeV/#it{c})"};
    const AxisSpec ptResolutionAxis{100, -0.5, 0.5, "#it{p}_{T}^{rec} - #it{p}_{T}^{gen} (GeV/#it{c})"};
    const AxisSpec massAxis{100, 1.3, 1.6, "m (GeV/#it{c}^{2})"};
    const AxisSpec massResolutionAxis{100, -0.1, 0.1, "m_{rec} - m_{gen} (GeV/#it{c}^{2})"};
    const AxisSpec nSigmaPiAxis{100, -5, 5, "n#sigma_{#pi}"};
    const AxisSpec sigmaMassAxis{100, 1.1, 1.4, "m (GeV/#it{c}^{2})"};
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // lambda1405 to sigmaminus
    rLambda1405.add("h2PtMass_0", "h2PtMass_0", {HistType::kTH2F, {ptAxis, massAxis}});
    rLambda1405.add("h2PtMassSigmaBeforeCuts_0", "h2PtMassSigmaBeforeCuts_0", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rLambda1405.add("h2PtMassSigma_0", "h2PtMassSigma_0", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rLambda1405.add("h2SigmaMassVsMass_0", "h2SigmaMassVsMass_0", {HistType::kTH2F, {massAxis, sigmaMassAxis}});
    rLambda1405.add("h2PtPiNSigma_0", "h2PtPiNSigma_0", {HistType::kTH2F, {ptPiAxis, nSigmaPiAxis}});
    rLambda1405.add("h2PtPiNSigmaTOF_0", "h2PtPiNSigmaTOF_0", {HistType::kTH2F, {ptPiAxis, nSigmaPiAxis}});
    // lambda1405 to sigmaplus
    rLambda1405.add("h2PtMass_1", "h2PtMass_1", {HistType::kTH2F, {ptAxis, massAxis}});
    rLambda1405.add("h2PtMassSigmaBeforeCuts_1", "h2PtMassSigmaBeforeCuts_1", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rLambda1405.add("h2PtMassSigma_1", "h2PtMassSigma_1", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rLambda1405.add("h2SigmaMassVsMass_1", "h2SigmaMassVsMass_1", {HistType::kTH2F, {massAxis, sigmaMassAxis}});
    rLambda1405.add("h2PtPiNSigma_1", "h2PtPiNSigma_1", {HistType::kTH2F, {ptPiAxis, nSigmaPiAxis}});
    rLambda1405.add("h2PtPiNSigmaTOF_1", "h2PtPiNSigmaTOF_1", {HistType::kTH2F, {ptPiAxis, nSigmaPiAxis}});

    if (doprocessMC) {
      // Add MC histograms if needed, to sigmaminus
      rLambda1405.add("h2MassResolution_0", "h2MassResolution_0", {HistType::kTH2F, {massAxis, massResolutionAxis}});
      rLambda1405.add("h2PtResolution_0", "h2PtResolution_0", {HistType::kTH2F, {ptAxis, ptResolutionAxis}});
      rLambda1405.add("h2PtMassMC_0", "h2PtMassMC_0", {HistType::kTH2F, {ptAxis, massAxis}});
      // Add MC histograms if needed, to sigmaplus
      rLambda1405.add("h2MassResolution_1", "h2MassResolution_1", {HistType::kTH2F, {massAxis, massResolutionAxis}});
      rLambda1405.add("h2PtResolution_1", "h2PtResolution_1", {HistType::kTH2F, {ptAxis, ptResolutionAxis}});
      rLambda1405.add("h2PtMassMC_1", "h2PtMassMC_1", {HistType::kTH2F, {ptAxis, massAxis}});
    }
  }

  template <typename Ttrack>
  bool selectPiTrack(const Ttrack& candidate, bool piFromSigma)
  {
    if (std::abs(candidate.tpcNSigmaPi()) > cutNSigmaTPC || candidate.tpcNClsFound() < cutNTPCClusPi || std::abs(candidate.eta()) > cutEtaDaught) {
      return false;
    }
    if (piFromSigma) {
      return true;
    }

    if (candidate.itsNCls() < cutNITSClusPi) {
      return false;
    }

    if (useTOF && !candidate.hasTOF()) {
      return false;
    }

    if (useTOF && std::abs(candidate.tofNSigmaPi()) > cutNSigmaTOF) {
      return false;
    }

    return true; // Track is selected
  }

  template <typename Ttrack>
  bool selectProTrack(const Ttrack& candidate, bool prFromSigma)
  {
    if (std::abs(candidate.tpcNSigmaPr()) > cutNSigmaTPC || candidate.tpcNClsFound() < cutNTPCClusPi || std::abs(candidate.eta()) > cutEtaDaught) {
      return false;
    }
    if (prFromSigma) {
      return true;
    }
    if (candidate.itsNCls() < cutNITSClusPi) {
      return false;
    }
    if (useTOF && !candidate.hasTOF()) {
      return false;
    }
    if (useTOF && std::abs(candidate.tofNSigmaPr()) > cutNSigmaTOF) {
      return false;
    }
    return true; // Track is selected
  }

  bool selectCandidate(aod::KinkCands::iterator const& sigmaCand, TracksFull const& tracks)
  {
    auto kinkDauTrack = sigmaCand.trackDaug_as<TracksFull>();
    bool isPiKink = selectPiTrack(kinkDauTrack, true);
    bool isProKink = selectProTrack(kinkDauTrack, true);
    if (!isPiKink && !isProKink) {
      return false;
    }

    if (isPiKink) {
      rLambda1405.fill(HIST("h2PtMassSigmaBeforeCuts_0"), sigmaCand.mothSign() * sigmaCand.ptMoth(), sigmaCand.mSigmaMinus());
      rLambda1405.fill(HIST("h2PtPiNSigma_0"), sigmaCand.mothSign() * kinkDauTrack.pt(), kinkDauTrack.tpcNSigmaPi());
    }
    if (isProKink) {
      rLambda1405.fill(HIST("h2PtMassSigmaBeforeCuts_1"), sigmaCand.mothSign() * sigmaCand.ptMoth(), sigmaCand.mSigmaPlus());
      rLambda1405.fill(HIST("h2PtPiNSigma_1"), sigmaCand.mothSign() * kinkDauTrack.pt(), kinkDauTrack.tpcNSigmaPr());
    }

    lambda1405Cand.isSigmaPlus = isProKink && (sigmaCand.mSigmaPlus() > o2::constants::physics::MassSigmaPlus - cutSigmaMass && sigmaCand.mSigmaPlus() < o2::constants::physics::MassSigmaPlus + cutSigmaMass);
    lambda1405Cand.isSigmaMinus = isPiKink && (sigmaCand.mSigmaMinus() > o2::constants::physics::MassSigmaMinus - cutSigmaMass && sigmaCand.mSigmaMinus() < o2::constants::physics::MassSigmaMinus + cutSigmaMass);
    if (!lambda1405Cand.isSigmaPlus && !lambda1405Cand.isSigmaMinus) {
      return false;
    }
    float sigmaRad = std::hypot(sigmaCand.xDecVtx(), sigmaCand.yDecVtx());
    if (std::abs(sigmaCand.dcaMothPv()) > cutDCAtoPVSigma || std::abs(sigmaCand.dcaDaugPv()) < cutDCAtoPVPiFromSigma || sigmaRad < cutSigmaRadius) {
      return false;
    }

    for (const auto& piTrack : tracks) {
      if (!doLSBkg) {
        if (piTrack.sign() == sigmaCand.mothSign()) {
          continue;
        }
      } else {
        if (piTrack.sign() != sigmaCand.mothSign()) {
          continue;
        }
      }
      if (!selectPiTrack(piTrack, false)) {
        continue;
      }
      auto sigmaMom = std::array{sigmaCand.pxMoth(), sigmaCand.pyMoth(), sigmaCand.pzMoth()};
      auto piMom = std::array{piTrack.px(), piTrack.py(), piTrack.pz()};
      float pt = std::hypot(sigmaMom[0] + piMom[0], sigmaMom[1] + piMom[1]);
      double massSigma = lambda1405Cand.isSigmaMinus ? sigmaCand.mSigmaMinus() : sigmaCand.mSigmaPlus();
      float invMass = RecoDecay::m(std::array{sigmaMom, piMom}, std::array{massSigma, o2::constants::physics::MassPiPlus});
      if (invMass < 1.3 || invMass > 1.5) {
        continue;
      }
      lambda1405Cand.kinkDauID = kinkDauTrack.globalIndex();
      lambda1405Cand.sigmaID = sigmaCand.globalIndex();
      lambda1405Cand.piID = piTrack.globalIndex();
      lambda1405Cand.mass = invMass;
      lambda1405Cand.sigmaMinusMass = sigmaCand.mSigmaMinus();
      lambda1405Cand.sigmaPlusMass = sigmaCand.mSigmaPlus();
      lambda1405Cand.sigmaSign = sigmaCand.mothSign();
      lambda1405Cand.pt = pt;
      lambda1405Cand.sigmaPt = sigmaCand.ptMoth();
      lambda1405Cand.piPt = piTrack.pt();
      lambda1405Cand.nSigmaTPCPi = piTrack.tpcNSigmaPi();
      if (useTOF) {
        lambda1405Cand.nSigmaTOFPi = piTrack.tofNSigmaPi();
      } else {
        lambda1405Cand.nSigmaTOFPi = -999; // Not used if TOF is not enabled
      }
      return true; // Candidate is selected
    }
    return false; // No valid pion track found
  }

  void processData(CollisionsFull::iterator const& collision, aod::KinkCands const& kinkCands, TracksFull const& tracks)
  {
    if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
      return;
    }
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    for (const auto& sigmaCand : kinkCands) {
      if (selectCandidate(sigmaCand, tracks)) {
        if (lambda1405Cand.isSigmaMinus) {
          rLambda1405.fill(HIST("h2PtMass_0"), lambda1405Cand.sigmaSign * lambda1405Cand.pt, lambda1405Cand.mass);
          rLambda1405.fill(HIST("h2PtMassSigma_0"), lambda1405Cand.sigmaSign * lambda1405Cand.sigmaPt, lambda1405Cand.sigmaMinusMass);
          rLambda1405.fill(HIST("h2SigmaMassVsMass_0"), lambda1405Cand.mass, lambda1405Cand.sigmaMinusMass);
          rLambda1405.fill(HIST("h2PtPiNSigmaTOF_0"), lambda1405Cand.sigmaSign * lambda1405Cand.piPt, lambda1405Cand.nSigmaTOFPi);
        }
        if (lambda1405Cand.isSigmaPlus) {
          rLambda1405.fill(HIST("h2PtMass_1"), lambda1405Cand.sigmaSign * lambda1405Cand.pt, lambda1405Cand.mass);
          rLambda1405.fill(HIST("h2PtMassSigma_1"), lambda1405Cand.sigmaSign * lambda1405Cand.sigmaPt, lambda1405Cand.sigmaPlusMass);
          rLambda1405.fill(HIST("h2SigmaMassVsMass_1"), lambda1405Cand.mass, lambda1405Cand.sigmaPlusMass);
          rLambda1405.fill(HIST("h2PtPiNSigmaTOF_1"), lambda1405Cand.sigmaSign * lambda1405Cand.piPt, lambda1405Cand.nSigmaTOFPi);
        }
      }
    }
  }
  PROCESS_SWITCH(lambda1405analysis, processData, "Data processing", true);

  void processMC(CollisionsFullMC const& collisions, aod::KinkCands const& kinkCands, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC, TracksFull const& tracks)
  {
    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
        continue;
      }
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      auto sigmaCandsPerCol = kinkCands.sliceBy(mKinkPerCol, collision.globalIndex());
      auto tracksPerCol = tracks.sliceBy(mPerColTracks, collision.globalIndex());
      for (const auto& sigmaCand : sigmaCandsPerCol) {
        if (selectCandidate(sigmaCand, tracksPerCol)) {
          // Do MC association
          auto mcLabPiKink = trackLabelsMC.rawIteratorAt(lambda1405Cand.kinkDauID);
          auto mcLabSigma = trackLabelsMC.rawIteratorAt(lambda1405Cand.sigmaID);
          auto mcLabPi = trackLabelsMC.rawIteratorAt(lambda1405Cand.piID);
          if (!mcLabSigma.has_mcParticle() || mcLabPiKink.has_mcParticle() || mcLabPi.has_mcParticle()) {
            continue; // Skip if no valid MC association
          }
          auto mcTrackPiKink = mcLabPiKink.mcParticle_as<aod::McParticles>();
          auto mcTrackSigma = mcLabSigma.mcParticle_as<aod::McParticles>();
          auto mcTrackPi = mcLabPi.mcParticle_as<aod::McParticles>();
          if (std::abs(mcTrackPiKink.pdgCode()) != 211 || std::abs(mcTrackSigma.pdgCode()) != 3122 || std::abs(mcTrackPi.pdgCode()) != 211) {
            continue; // Skip if not a valid pion or Sigma candidate
          }
          if (!mcTrackPiKink.has_mothers() || !mcTrackSigma.has_mothers() || !mcTrackPi.has_mothers()) {
            continue; // Skip if no mothers found
          }
          // check if kink pi comes from the sigma
          bool isPiFromSigma = false;
          for (const auto& piMother : mcTrackPiKink.mothers_as<aod::McParticles>()) {
            if (piMother.globalIndex() == mcTrackSigma.globalIndex()) {
              isPiFromSigma = true;
              break; // Found the mother, exit loop
            }
          }
          if (!isPiFromSigma) {
            continue; // Skip if the pion does not come from the Sigma
          }
          // check that labpi and labsigma have the same mother (a lambda1405 candidate)
          int lambda1405Id = -1;
          for (const auto& piMother : mcTrackPi.mothers_as<aod::McParticles>()) {
            for (const auto& sigmaMother : mcTrackSigma.mothers_as<aod::McParticles>()) {
              if (piMother.globalIndex() == sigmaMother.globalIndex() && std::abs(piMother.pdgCode()) == lambda1405PdgCode) {
                lambda1405Id = piMother.globalIndex();
                break; // Found the mother, exit loop
              }
            }
          }
          if (lambda1405Id == -1) {
            continue; // Skip if the Sigma and pion do not share the same lambda1405 candidate
          }
          auto lambda1405Mother = particlesMC.rawIteratorAt(lambda1405Id);
          LOG(info) << "Particle selected!";
          float lambda1405Mass = std::sqrt(lambda1405Mother.e() * lambda1405Mother.e() - lambda1405Mother.p() * lambda1405Mother.p());
          if (lambda1405Cand.isSigmaMinus) {
            rLambda1405.fill(HIST("h2PtMass_0"), lambda1405Cand.sigmaSign * lambda1405Cand.pt, lambda1405Cand.mass);
            rLambda1405.fill(HIST("h2PtMassSigma_0"), lambda1405Cand.sigmaSign * lambda1405Cand.sigmaPt, lambda1405Cand.sigmaMinusMass);
            rLambda1405.fill(HIST("h2SigmaMassVsMass_0"), lambda1405Cand.mass, lambda1405Cand.sigmaMinusMass);
            rLambda1405.fill(HIST("h2PtPiNSigma_0"), lambda1405Cand.piPt, lambda1405Cand.nSigmaTPCPi);
            rLambda1405.fill(HIST("h2MassResolution_0"), lambda1405Mass, lambda1405Mass - lambda1405Cand.mass);
            rLambda1405.fill(HIST("h2PtResolution_0"), lambda1405Cand.pt, lambda1405Cand.pt - lambda1405Mother.pt());
          }
          if (lambda1405Cand.isSigmaPlus) {
            rLambda1405.fill(HIST("h2PtMass_1"), lambda1405Cand.sigmaSign * lambda1405Cand.pt, lambda1405Cand.mass);
            rLambda1405.fill(HIST("h2PtMassSigma_1"), lambda1405Cand.sigmaSign * lambda1405Cand.sigmaPt, lambda1405Cand.sigmaPlusMass);
            rLambda1405.fill(HIST("h2SigmaMassVsMass_1"), lambda1405Cand.mass, lambda1405Cand.sigmaPlusMass);
            rLambda1405.fill(HIST("h2PtPiNSigma_1"), lambda1405Cand.piPt, lambda1405Cand.nSigmaTPCPi);
            rLambda1405.fill(HIST("h2MassResolution_1"), lambda1405Mass, lambda1405Mass - lambda1405Cand.mass);
            rLambda1405.fill(HIST("h2PtResolution_1"), lambda1405Cand.pt, lambda1405Cand.pt - lambda1405Mother.pt());
          }
        }
      }
    }
    // Loop over generated particles to fill MC histograms
    for (const auto& mcPart : particlesMC) {
      if (std::abs(mcPart.pdgCode()) != lambda1405PdgCode) {
        continue; // Only consider Lambda(1405) candidates
      }

      if (!mcPart.has_daughters()) {
        continue; // Skip if no daughters
      }
      // Check if the Lambda(1405) has a Sigma daughter
      bool hasSigmaDaughter = false;
      int dauPdgCode = 0;
      for (const auto& daughter : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(daughter.pdgCode()) == 3122 || std::abs(daughter.pdgCode()) == 3222) { // Sigma PDG code
          hasSigmaDaughter = true;
          dauPdgCode = daughter.pdgCode();
          break; // Found a Sigma daughter, exit loop
        }
      }
      if (!hasSigmaDaughter) {
        continue; // Skip if no Sigma daughter found
      }

      float mcMass = std::sqrt(mcPart.e() * mcPart.e() - mcPart.p() * mcPart.p());
      dauPdgCode ? rLambda1405.fill(HIST("h2PtMassMC_0"), mcPart.pt(), mcMass) : rLambda1405.fill(HIST("h2PtMassMC_1"), mcPart.pt(), mcMass);
    }
  }
  PROCESS_SWITCH(lambda1405analysis, processMC, "MC processing", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambda1405analysis>(cfgc)};
}
