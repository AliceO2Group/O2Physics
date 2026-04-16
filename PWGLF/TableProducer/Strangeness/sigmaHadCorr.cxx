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

/// \file   sigmaHadCorr.cxx
/// \brief Analysis task for sigma-hadron correlations
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>

#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFSigmaHadTables.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullPr>;
using TracksFullMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullPr, aod::McTrackLabels>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

struct sigmaHadCand {

  float ptHad() const
  {
    return std::hypot(pxHad, pyHad);
  }
  float sigmaPt() const
  {
    return std::hypot(sigmaPx, sigmaPy);
  }

  int sigmaCharge = 0;       // Charge of the sigma candidate
  float sigmaMass = -1.f;    // Mass of the Sigma candidate
  float sigmaPx = -1;        // Px of the Sigma candidate
  float sigmaPy = -1;        // Py of the Sigma candidate
  float sigmaPz = -1;        // Pz of the Sigma candidate
  float sigmaDauPx = -1;     // Px of the daughter track from Sigma decay
  float sigmaDauPy = -1;     // Py of the daughter track from Sigma decay
  float sigmaDauPz = -1;     // Pz of the daughter track from Sigma decay
  float sigmaDecRadius = -1; // Decay radius of the Sigma candidate
  float sigmaCosPA = -1;     // Cosine of pointing angle of the Sigma candidate

  int chargeHad = 0; // Charge of the hadron candidate
  float pxHad = -1;  // Px of the hadron candidate
  float pyHad = -1;  // Py of the hadron candidate
  float pzHad = -1;  // Pz of the hadron candidate

  float nSigmaTPCHad = -1; // Number of sigmas for the hadron candidate
  float nSigmaTOFHad = -1; // Number of sigmas for the hadron candidate using TOF

  int kinkDauID = -1; // ID of the pion from Sigma decay in MC
  int sigmaID = -1;   // ID of the Sigma candidate in MC
  int hadID = -1;     // ID of the hadron candidate in MC
};

struct sigmaHadCorrTask {

  std::vector<sigmaHadCand> sigmaHadCandidates;        // Vector to store Sigma-hadron candidates
  Produces<aod::SigmaProtonCands> outputDataTable;     // Output table for Sigma-hadron candidates
  Produces<aod::SigmaProtonMCCands> outputDataTableMC; // Output table for Sigma-hadron candidates in MC
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rSigmaHad{"sigmaHad", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // Configurable for event selection
  Configurable<float> cutzvertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};

  Configurable<bool> doSigmaPion{"doSigmaPion", false, "If true, pair Sigma with pions instead of protons"};
  Configurable<bool> doSigmaMinus{"doSigmaMinus", true, "If true, pair Sigma- candidates, else Sigma+"};
  Configurable<float> cutMaxKStar{"cutMaxKStar", 1.5, "Maximum k* for Sigma-hadron pairs (GeV/c)"};

  Configurable<float> minPtSigma{"minPtSigma", 1.f, "Minimum pT for Sigma candidates (GeV/c)"};
  Configurable<bool> useRecalculatedSigmaMomentum{"useRecalculatedSigmaMomentum", true, "If true, compute k* using Sigma momentum recalculated from daughter kinematics"};
  Configurable<float> cutDCAtoPVSigma{"cutDCAtoPVSigma", 0.1f, "Max DCA to primary vertex for Sigma candidates (cm)"};
  Configurable<float> cutSigmaRadius{"cutSigmaRadius", 20.f, "Minimum radius for Sigma candidates (cm)"};
  Configurable<float> cutSigmaMass{"cutSigmaMass", 0.3, "Sigma mass window (GeV/c^2)"};
  Configurable<float> alphaAPCut{"alphaAPCut", 0., "Alpha AP cut for Sigma candidates"};
  Configurable<float> qtAPCutLow{"qtAPCutLow", 0.15, "Lower qT AP cut for Sigma candidates (GeV/c)"};
  Configurable<float> qtAPCutHigh{"qtAPCutHigh", 0.2, "Upper qT AP cut for Sigma candidates (GeV/c)"};
  Configurable<float> cutEtaDaught{"cutEtaDaughter", 0.8f, "Eta cut for daughter tracks"};
  Configurable<float> ptMinTOFKinkDau{"ptMinTOFKinkDau", 0.75f, "Minimum pT to require TOF for kink daughter PID (GeV/c)"};
  Configurable<bool> applyTOFPIDKinkDaughter{"applyTOFPIDKinkDaughter", false, "If true, apply TOF PID cut to the kink daughter track"};

  Configurable<float> cutNITSClusHad{"cutNITSClusHad", 5, "Minimum number of ITS clusters for hadron track"};
  Configurable<float> cutNTPCClusHad{"cutNTPCClusHad", 90, "Minimum number of TPC clusters for hadron track"};
  Configurable<float> ptMinHad{"ptMinHad", 0.2f, "Minimum pT for hadron track (GeV/c)"};
  Configurable<float> ptMinTOFHad{"ptMinTOFHad", 0.75f, "Minimum pT to require TOF for hadron PID (GeV/c)"};
  Configurable<float> cutNSigmaTPC{"cutNSigmaTPC", 3, "TPC nSigma cut for hadron track"};
  Configurable<float> cutNSigmaTOF{"cutNSigmaTOF", 3, "TOF nSigma cut for hadron track"};

  Configurable<bool> fillOutputTree{"fillOutputTree", true, "If true, fill the output tree with Sigma-hadron candidates"};
  Configurable<bool> fillSparseInvMassKstar{"fillSparseInvMassKstar", false, "If true, fill THn with invmass, k*, sigma charge, proton charge, sigma decay radius, cosPA, sigma pt"};

  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0, 40.0, 80.0, 500.0}, "Mixing bins - number of contributor"};
  Configurable<int> nEvtMixingBkg{"nEvtMixingBkg", 5, "Number of events to mix for background reconstruction"};

  Preslice<aod::KinkCands> kinkCandsPerCollisionPreslice = aod::kinkcand::collisionId;
  Preslice<TracksFull> tracksPerCollisionPreslice = aod::track::collisionId;
  Preslice<TracksFullMC> tracksMCPerCollisionPreslice = aod::track::collisionId;

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec ptAxis{100, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec massResolutionAxis{100, -0.1, 0.1, "m_{rec} - m_{gen} (GeV/#it{c}^{2})"};
    const AxisSpec nSigmaHadAxis{100, -5, 5, "n#sigma_{had}"};
    const AxisSpec sigmaMassAxis{50, 1.1, 1.3, "m (GeV/#it{c}^{2})"};
    const AxisSpec kStarAxis{200, 0.0, 2., "k* (GeV/#it{c})"};
    const AxisSpec ptHadAxis{100, 0.0, 6.0, "#it{p}_{T,had} (GeV/#it{c})"};
    const AxisSpec sigmaPtAxis{100, 0.0, 6.0, "#it{p}_{T,#Sigma} (GeV/#it{c})"};
    const AxisSpec sigmaPtAxisCoarse{30, 0.0, 6.0, "#it{p}_{T,#Sigma} (GeV/#it{c})"};
    const AxisSpec sigmaChargeAxis{2, -1.5, 1.5, "#Sigma charge"};
    const AxisSpec hadronChargeAxis{2, -1.5, 1.5, "Hadron charge"};
    const AxisSpec sigmaDecRadiusAxis{25, 14.5, 40.5, "#Sigma decay radius (cm)"};
    const AxisSpec sigmaDecRadiusAxisCoarse{5, 14.5, 40.5, "#Sigma decay radius (cm)"};
    const AxisSpec cosPAAxis{50, 0.9, 1.0, "cos(PA)"};
    const AxisSpec cosPAAxisCoarse{5, 0.9, 1.0, "cos(PA)"};
    const AxisSpec alphaAPAxis{100, -1.0, 1.0, "#alpha_{AP}"};
    const AxisSpec qtAPAxis{100, 0.0, 0.5, "q_{T,AP} (GeV/#it{c})"};
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};

    // qa histograms
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    // Dedicated QA folder
    rSigmaHad.add("QA/hHadronPt", "Hadron #it{p}_{T}", {HistType::kTH1F, {ptHadAxis}});
    rSigmaHad.add("QA/h2TPCNSigmaHadVsPtHad", "TPC n#sigma_{had} vs #it{p}_{T,had}", {HistType::kTH2F, {ptHadAxis, nSigmaHadAxis}});
    rSigmaHad.add("QA/h2TOFNSigmaHadVsPtHad", "TOF n#sigma_{had} vs #it{p}_{T,had}", {HistType::kTH2F, {ptHadAxis, nSigmaHadAxis}});
    rSigmaHad.add("QA/hSigmaPt", "#Sigma #it{p}_{T}", {HistType::kTH1F, {sigmaPtAxis}});
    rSigmaHad.add("QA/hSigmaPtRecal", "#Sigma #it{p}_{T} recalculated", {HistType::kTH1F, {sigmaPtAxis}});
    rSigmaHad.add("QA/h2InvMassVsPtSigma", "m_{#Sigma} vs #it{p}_{T,#Sigma}", {HistType::kTH2F, {sigmaPtAxis, sigmaMassAxis}});
    rSigmaHad.add("QA/h2QtAPvsAlphaAP", "q_{T,AP} vs #alpha_{AP}", {HistType::kTH2F, {alphaAPAxis, qtAPAxis}});

    if (fillSparseInvMassKstar) {
      rSigmaHad.add("hSparseSigmaHad",
                    "7D THnSparse: invmass, k*, sigma charge, hadron charge, sigma decay radius, cosPA, sigma pt",
                    {HistType::kTHnSparseF, {sigmaMassAxis, kStarAxis, sigmaChargeAxis, hadronChargeAxis, sigmaDecRadiusAxisCoarse, cosPAAxisCoarse, sigmaPtAxisCoarse}});
      rSigmaHad.add("hSparseSigmaHadMC",
                    "8D THnSparse (MC): invmass, k*, sigma charge, hadron charge, sigma decay radius, cosPA, sigma pt, k* gen",
                    {HistType::kTHnSparseF, {sigmaMassAxis, kStarAxis, sigmaChargeAxis, hadronChargeAxis, sigmaDecRadiusAxisCoarse, cosPAAxisCoarse, sigmaPtAxisCoarse, kStarAxis}});
    }

    LOG(info) << "Sigma-hadron correlation task initialized";
    LOG(info) << "Process SE enabled: " << doprocessSameEvent;
    LOG(info) << "Process ME enabled: " << doprocessMixedEvent;
    LOG(info) << "Process SE MC enabled: " << doprocessSameEventMC;
    LOG(info) << "Process ME MC enabled: " << doprocessMixedEventMC;
    LOG(info) << "Pairing mode: " << (doSigmaPion ? "Sigma-pion" : "Sigma-proton");
  }

  float getAlphaAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    std::array<float, 3> momMissing = {momMother[0] - momKink[0], momMother[1] - momKink[1], momMother[2] - momKink[2]};
    float lQlP = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float lQlN = std::inner_product(momMother.begin(), momMother.end(), momMissing.begin(), 0.f);
    return (lQlP - lQlN) / (lQlP + lQlN);
  }

  float getQtAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    float dp = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float p2V0 = std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f);
    float p2A = std::inner_product(momKink.begin(), momKink.end(), momKink.begin(), 0.f);
    return std::sqrt(p2A - dp * dp / p2V0);
  }

  float getCosPA(const std::array<float, 3>& momMother, const std::array<float, 3>& decayVertex, const std::array<float, 3>& primaryVertex)
  {
    std::array<float, 3> decayVec = {decayVertex[0] - primaryVertex[0], decayVertex[1] - primaryVertex[1], decayVertex[2] - primaryVertex[2]};
    float dotProduct = std::inner_product(momMother.begin(), momMother.end(), decayVec.begin(), 0.f);
    float momMotherMag = std::sqrt(std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f));
    float decayVecMag = std::sqrt(std::inner_product(decayVec.begin(), decayVec.end(), decayVec.begin(), 0.f));
    return dotProduct / (momMotherMag * decayVecMag);
  }

  float getRecalculatedSigmaMomentum(float sigmaPx, float sigmaPy, float sigmaPz, float sigmaDauPx, float sigmaDauPy, float sigmaDauPz)
  {
    // Sigma- -> n + pi-  (charged daughter = pion, neutral daughter = neutron)
    // Sigma+ -> p + pi0  (charged daughter = proton, neutral daughter = pi0)
    float massChargedDau = doSigmaMinus ? o2::constants::physics::MassPionCharged : o2::constants::physics::MassProton;
    float massNeutralDau = doSigmaMinus ? o2::constants::physics::MassNeutron : o2::constants::physics::MassPionNeutral;
    float massSigma = doSigmaMinus ? o2::constants::physics::MassSigmaMinus : o2::constants::physics::MassSigmaPlus;

    float pMother = std::sqrt(sigmaPx * sigmaPx + sigmaPy * sigmaPy + sigmaPz * sigmaPz);
    if (pMother < 1e-6f) {
      return -999.f;
    }
    float versorX = sigmaPx / pMother;
    float versorY = sigmaPy / pMother;
    float versorZ = sigmaPz / pMother;
    float eChDau = std::sqrt(massChargedDau * massChargedDau + sigmaDauPx * sigmaDauPx + sigmaDauPy * sigmaDauPy + sigmaDauPz * sigmaDauPz);
    float a = versorX * sigmaDauPx + versorY * sigmaDauPy + versorZ * sigmaDauPz;
    float K = massSigma * massSigma + massChargedDau * massChargedDau - massNeutralDau * massNeutralDau;
    float A = 4.f * (eChDau * eChDau - a * a);
    float B = -4.f * a * K;
    float C = 4.f * eChDau * eChDau * massSigma * massSigma - K * K;
    if (std::abs(A) < 1e-6f) {
      return -999.f;
    }
    float D = B * B - 4.f * A * C;
    if (D < 0.f) {
      return -999.f;
    }
    float sqrtD = std::sqrt(D);
    float P1 = (-B + sqrtD) / (2.f * A);
    float P2 = (-B - sqrtD) / (2.f * A);
    if (P2 < 0.f && P1 < 0.f) {
      return -999.f;
    }
    if (P2 < 0.f) {
      return P1;
    }
    float p1Diff = std::abs(P1 - pMother);
    float p2Diff = std::abs(P2 - pMother);
    return (p1Diff < p2Diff) ? P1 : P2;
  }

  std::array<float, 3> getSigmaMomentumForKstar(float sigmaPx, float sigmaPy, float sigmaPz, float sigmaDauPx, float sigmaDauPy, float sigmaDauPz)
  {
    std::array<float, 3> sigmaMomentum = {sigmaPx, sigmaPy, sigmaPz};
    if (!useRecalculatedSigmaMomentum) {
      return sigmaMomentum;
    }

    float pNew = getRecalculatedSigmaMomentum(sigmaPx, sigmaPy, sigmaPz, sigmaDauPx, sigmaDauPy, sigmaDauPz);
    if (pNew <= 0.f) {
      return sigmaMomentum;
    }

    float pOld = std::sqrt(sigmaPx * sigmaPx + sigmaPy * sigmaPy + sigmaPz * sigmaPz);
    if (pOld <= 0.f) {
      return sigmaMomentum;
    }

    float scale = pNew / pOld;
    sigmaMomentum[0] *= scale;
    sigmaMomentum[1] *= scale;
    sigmaMomentum[2] *= scale;
    return sigmaMomentum;
  }

  float getHadTrackMass()
  {
    return doSigmaPion ? o2::constants::physics::MassPionCharged : o2::constants::physics::MassProton;
  }

  float getSigmaMassForKstar()
  {
    return doSigmaMinus ? o2::constants::physics::MassSigmaMinus : o2::constants::physics::MassSigmaPlus;
  }

  template <typename Ttrack>
  float getTPCNSigmaHad(const Ttrack& track)
  {
    return doSigmaPion ? track.tpcNSigmaPi() : track.tpcNSigmaPr();
  }

  template <typename Ttrack>
  float getTOFNSigmaHad(const Ttrack& track)
  {
    return doSigmaPion ? track.tofNSigmaPi() : track.tofNSigmaPr();
  }

  TLorentzVector trackSum, PartOneCMS, PartTwoCMS, trackRelK;
  float getKStar(float sigmaPx, float sigmaPy, float sigmaPz, float pxHad, float pyHad, float pzHad)
  {
    TLorentzVector part1; // Sigma
    TLorentzVector part2; // Hadron track (proton/pion)
    part1.SetXYZM(sigmaPx, sigmaPy, sigmaPz, getSigmaMassForKstar());
    part2.SetXYZM(pxHad, pyHad, pzHad, getHadTrackMass());
    trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    PartOneCMS.SetXYZM(part1.Px(), part1.Py(), part1.Pz(), part1.M());
    PartTwoCMS.SetXYZM(part2.Px(), part2.Py(), part2.Pz(), part2.M());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  template <typename Ttrack>
  bool selectHadTrack(const Ttrack& candidate)
  {
    if (candidate.pt() < ptMinHad) {
      return false;
    }
    if (std::abs(getTPCNSigmaHad(candidate)) > cutNSigmaTPC || candidate.tpcNClsFound() < cutNTPCClusHad || std::abs(candidate.eta()) > cutEtaDaught) {
      return false;
    }

    if (candidate.itsNCls() < cutNITSClusHad) {
      return false;
    }

    if (candidate.pt() >= ptMinTOFHad) {
      if (!candidate.hasTOF()) {
        return false;
      }
      if (std::abs(getTOFNSigmaHad(candidate)) > cutNSigmaTOF) {
        return false;
      }
    }
    return true; // Track is selected
  }

  template <typename Ttrack>
  bool selectSigma(aod::KinkCands::iterator const& sigmaCand, Ttrack const& kinkDauTrack)
  {
    float mass = doSigmaMinus ? sigmaCand.mSigmaMinus() : sigmaCand.mSigmaPlus();
    std::array<float, 3> momMoth = {sigmaCand.pxMoth(), sigmaCand.pyMoth(), sigmaCand.pzMoth()};
    std::array<float, 3> momDaug = {sigmaCand.pxDaug(), sigmaCand.pyDaug(), sigmaCand.pzDaug()};
    float alphaAP = getAlphaAP(momMoth, momDaug);
    float qtAP = getQtAP(momMoth, momDaug);

    if (alphaAP > alphaAPCut || (qtAP < qtAPCutLow || qtAP > qtAPCutHigh)) {
      return false;
    }
    float decRad = std::hypot(sigmaCand.xDecVtx(), sigmaCand.yDecVtx());
    if (decRad < cutSigmaRadius) {
      return false;
    }

    if (doSigmaMinus) {
      if (mass < o2::constants::physics::MassSigmaMinus - cutSigmaMass || mass > o2::constants::physics::MassSigmaMinus + cutSigmaMass) {
        return false;
      }
      if (std::abs(kinkDauTrack.tpcNSigmaPi()) > cutNSigmaTPC) {
        return false;
      }
    } else {
      if (mass < o2::constants::physics::MassSigmaPlus - cutSigmaMass || mass > o2::constants::physics::MassSigmaPlus + cutSigmaMass) {
        return false;
      }
      if (std::abs(kinkDauTrack.tpcNSigmaPr()) > cutNSigmaTPC) {
        return false;
      }
    }

    if (applyTOFPIDKinkDaughter) {
      if (kinkDauTrack.pt() >= ptMinTOFKinkDau) {
        if (!kinkDauTrack.hasTOF()) {
          return false;
        }
        float kinkDauTOFNSigma = doSigmaMinus ? kinkDauTrack.tofNSigmaPi() : kinkDauTrack.tofNSigmaPr();
        if (std::abs(kinkDauTOFNSigma) > cutNSigmaTOF) {
          return false;
        }
      }
    }

    if (std::abs(sigmaCand.dcaMothPv()) > cutDCAtoPVSigma) {
      return false;
    }
    return true;
  }

  template <typename Ttrack, typename Tcollision>
  void fillTreeAndHistograms(aod::KinkCands const& kinkCands, Ttrack const& tracksDauSigma, Ttrack const& tracks, Tcollision const& collision, bool isMC)
  {
    for (const auto& sigmaCand : kinkCands) {
      auto kinkDauTrack = tracksDauSigma.rawIteratorAt(sigmaCand.trackDaugId());
      if (!selectSigma(sigmaCand, kinkDauTrack)) {
        continue;
      }

      std::array<float, 3> momMothAll = {sigmaCand.pxMoth(), sigmaCand.pyMoth(), sigmaCand.pzMoth()};
      std::array<float, 3> momDaugAll = {sigmaCand.pxDaug(), sigmaCand.pyDaug(), sigmaCand.pzDaug()};
      float alphaAP = getAlphaAP(momMothAll, momDaugAll);
      float qtAP = getQtAP(momMothAll, momDaugAll);
      rSigmaHad.fill(HIST("QA/h2QtAPvsAlphaAP"), alphaAP, qtAP);

      auto sigmaPRecal = getSigmaMomentumForKstar(sigmaCand.pxMoth(), sigmaCand.pyMoth(), sigmaCand.pzMoth(),
                                                  sigmaCand.pxDaug(), sigmaCand.pyDaug(), sigmaCand.pzDaug());
      float sigmaPtRecal = std::hypot(sigmaPRecal[0], sigmaPRecal[1]);
      float sigmaMassForQa = doSigmaMinus ? sigmaCand.mSigmaMinus() : sigmaCand.mSigmaPlus();

      if (sigmaPtRecal < minPtSigma) {
        continue;
      }

      rSigmaHad.fill(HIST("QA/hSigmaPt"), sigmaCand.ptMoth());
      rSigmaHad.fill(HIST("QA/hSigmaPtRecal"), sigmaPtRecal);
      rSigmaHad.fill(HIST("QA/h2InvMassVsPtSigma"), sigmaPtRecal, sigmaMassForQa);

      for (const auto& hadTrack : tracks) {
        if (hadTrack.globalIndex() == sigmaCand.trackDaugId()) {
          continue;
        }

        if (!selectHadTrack(hadTrack)) {
          continue;
        }

        sigmaHadCand candidate;
        candidate.sigmaCharge = sigmaCand.mothSign();
        candidate.sigmaPx = sigmaCand.pxMoth();
        candidate.sigmaPy = sigmaCand.pyMoth();
        candidate.sigmaPz = sigmaCand.pzMoth();
        candidate.sigmaDauPx = sigmaCand.pxDaug();
        candidate.sigmaDauPy = sigmaCand.pyDaug();
        candidate.sigmaDauPz = sigmaCand.pzDaug();
        candidate.sigmaDecRadius = std::hypot(sigmaCand.xDecVtx(), sigmaCand.yDecVtx());

        std::array<float, 3> momMoth = {sigmaCand.pxMoth(), sigmaCand.pyMoth(), sigmaCand.pzMoth()};
        std::array<float, 3> decayVtx = {sigmaCand.xDecVtx(), sigmaCand.yDecVtx(), sigmaCand.zDecVtx()};
        std::array<float, 3> primaryVtx = {collision.posX(), collision.posY(), collision.posZ()};
        candidate.sigmaCosPA = getCosPA(momMoth, decayVtx, primaryVtx);

        candidate.chargeHad = hadTrack.sign();
        candidate.pxHad = hadTrack.px();
        candidate.pyHad = hadTrack.py();
        candidate.pzHad = hadTrack.pz();
        candidate.nSigmaTPCHad = getTPCNSigmaHad(hadTrack);
        candidate.nSigmaTOFHad = getTOFNSigmaHad(hadTrack);
        candidate.sigmaMass = doSigmaMinus ? sigmaCand.mSigmaMinus() : sigmaCand.mSigmaPlus();

        candidate.sigmaID = sigmaCand.trackMothId();
        candidate.kinkDauID = sigmaCand.trackDaugId();
        candidate.hadID = hadTrack.globalIndex();

        float kStar = getKStar(sigmaPRecal[0], sigmaPRecal[1], sigmaPRecal[2], candidate.pxHad, candidate.pyHad, candidate.pzHad);
        if (kStar > cutMaxKStar) {
          continue;
        }

        rSigmaHad.fill(HIST("QA/hHadronPt"), candidate.ptHad());
        rSigmaHad.fill(HIST("QA/h2TPCNSigmaHadVsPtHad"), candidate.ptHad(), candidate.nSigmaTPCHad);
        if (hadTrack.hasTOF()) {
          rSigmaHad.fill(HIST("QA/h2TOFNSigmaHadVsPtHad"), candidate.ptHad(), candidate.nSigmaTOFHad);
        }
        if (fillSparseInvMassKstar && !isMC) {
          rSigmaHad.fill(HIST("hSparseSigmaHad"),
                         candidate.sigmaMass,
                         kStar,
                         candidate.sigmaCharge,
                         candidate.chargeHad,
                         candidate.sigmaDecRadius,
                         candidate.sigmaCosPA,
                         sigmaPtRecal);
        }
        sigmaHadCandidates.push_back(candidate);
      }
    }
  }

  void processSameEvent(CollisionsFull const& collisions, aod::KinkCands const& kinkCands, TracksFull const& tracks)
  {
    for (auto const& collision : collisions) {

      sigmaHadCandidates.clear();
      auto kinkCands_c = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision.globalIndex());
      auto tracks_c = tracks.sliceBy(tracksPerCollisionPreslice, collision.globalIndex());
      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
        continue;
      }
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      fillTreeAndHistograms(kinkCands_c, tracks, tracks_c, collision, false);
      if (fillOutputTree) {
        // Fill output table
        for (const auto& candidate : sigmaHadCandidates) {
          outputDataTable(candidate.sigmaCharge,
                          candidate.sigmaPx,
                          candidate.sigmaPy,
                          candidate.sigmaPz,
                          candidate.sigmaDauPx,
                          candidate.sigmaDauPy,
                          candidate.sigmaDauPz,
                          candidate.sigmaDecRadius,
                          candidate.sigmaCosPA,
                          candidate.chargeHad,
                          candidate.pxHad,
                          candidate.pyHad,
                          candidate.pzHad,
                          candidate.nSigmaTPCHad,
                          candidate.nSigmaTOFHad);
        }
      }
    }
  }
  PROCESS_SWITCH(sigmaHadCorrTask, processSameEvent, "Process Same event", true);

  // Processing Event Mixing
  SliceCache cache;
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  BinningType colBinning{{CfgVtxBins, CfgMultBins}, true};

  void processMixedEvent(const CollisionsFull& collisions, const aod::KinkCands& kinkCands, const TracksFull& tracks)
  {
    for (auto const& [collision1, collision2] :
         selfCombinations(colBinning, nEvtMixingBkg, -1, collisions, collisions)) {
      if (collision1.index() == collision2.index())
        continue;

      sigmaHadCandidates.clear();
      if (std::abs(collision1.posZ()) > cutzvertex || !collision1.sel8()) {
        continue;
      }
      if (std::abs(collision2.posZ()) > cutzvertex || !collision2.sel8()) {
        continue;
      }
      auto kinkCands_c1 = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision1.globalIndex());
      auto tracks_c1 = tracks.sliceBy(tracksPerCollisionPreslice, collision1.globalIndex());
      auto tracks_c2 = tracks.sliceBy(tracksPerCollisionPreslice, collision2.globalIndex());
      fillTreeAndHistograms(kinkCands_c1, tracks, tracks_c2, collision1, false);

      if (fillOutputTree) {
        // Fill output table
        for (const auto& candidate : sigmaHadCandidates) {
          outputDataTable(candidate.sigmaCharge,
                          candidate.sigmaPx,
                          candidate.sigmaPy,
                          candidate.sigmaPz,
                          candidate.sigmaDauPx,
                          candidate.sigmaDauPy,
                          candidate.sigmaDauPz,
                          candidate.sigmaDecRadius,
                          candidate.sigmaCosPA,
                          candidate.chargeHad,
                          candidate.pxHad,
                          candidate.pyHad,
                          candidate.pzHad,
                          candidate.nSigmaTPCHad,
                          candidate.nSigmaTOFHad);
        }
      }
    }
    LOG(debug) << "Processing mixed event";
  }
  PROCESS_SWITCH(sigmaHadCorrTask, processMixedEvent, "Process Mixed event", false);

  void processSameEventMC(CollisionsFullMC const& collisions, aod::KinkCands const& kinkCands, TracksFullMC const& tracks, aod::McParticles const&)
  {
    for (auto const& collision : collisions) {

      sigmaHadCandidates.clear();
      auto kinkCands_c = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision.globalIndex());
      auto tracks_c = tracks.sliceBy(tracksMCPerCollisionPreslice, collision.globalIndex());

      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
        continue;
      }
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      fillTreeAndHistograms(kinkCands_c, tracks, tracks_c, collision, true);
      for (const auto& candidate : sigmaHadCandidates) {
        auto mcLabelSigma = tracks.rawIteratorAt(candidate.sigmaID);
        auto mcLabelSigmaDau = tracks.rawIteratorAt(candidate.kinkDauID);
        auto mcLabelHad = tracks.rawIteratorAt(candidate.hadID);

        if (!mcLabelSigma.has_mcParticle() || !mcLabelSigmaDau.has_mcParticle() || !mcLabelHad.has_mcParticle()) {
          continue; // Skip candidates where MC truth is not available
        }

        auto mcPartSigma = mcLabelSigma.mcParticle_as<aod::McParticles>();
        auto mcPartSigmaDau = mcLabelSigmaDau.mcParticle_as<aod::McParticles>();
        auto mcPartHad = mcLabelHad.mcParticle_as<aod::McParticles>();
        auto pdgSigma = mcPartSigma.pdgCode();
        auto pdgSigmaDau = mcLabelSigmaDau.has_mcParticle() ? mcPartSigmaDau.pdgCode() : -999;
        auto pdgHad = mcLabelHad.has_mcParticle() ? mcPartHad.pdgCode() : -999;

        float sigmaPtGen = std::hypot(mcPartSigma.px(), mcPartSigma.py());
        float hadPtGen = std::hypot(mcPartHad.px(), mcPartHad.py());
        float kStarGen = getKStar(mcPartSigma.px(), mcPartSigma.py(), mcPartSigma.pz(), mcPartHad.px(), mcPartHad.py(), mcPartHad.pz());

        if (fillSparseInvMassKstar) {
          auto sigmaMomForKstar = getSigmaMomentumForKstar(candidate.sigmaPx, candidate.sigmaPy, candidate.sigmaPz,
                                                           candidate.sigmaDauPx, candidate.sigmaDauPy, candidate.sigmaDauPz);
          float kStarRec = getKStar(sigmaMomForKstar[0], sigmaMomForKstar[1], sigmaMomForKstar[2], candidate.pxHad, candidate.pyHad, candidate.pzHad);
          float sigmaPtUsed = std::hypot(sigmaMomForKstar[0], sigmaMomForKstar[1]);
          rSigmaHad.fill(HIST("hSparseSigmaHadMC"),
                         candidate.sigmaMass,
                         kStarRec,
                         candidate.sigmaCharge,
                         candidate.chargeHad,
                         candidate.sigmaDecRadius,
                         candidate.sigmaCosPA,
                         sigmaPtUsed,
                         kStarGen);
        }

        if (fillOutputTree) {
          outputDataTableMC(candidate.sigmaCharge,
                            candidate.sigmaPx,
                            candidate.sigmaPy,
                            candidate.sigmaPz,
                            candidate.sigmaDauPx,
                            candidate.sigmaDauPy,
                            candidate.sigmaDauPz,
                            candidate.sigmaDecRadius,
                            candidate.sigmaCosPA,
                            candidate.chargeHad,
                            candidate.pxHad,
                            candidate.pyHad,
                            candidate.pzHad,
                            candidate.nSigmaTPCHad,
                            candidate.nSigmaTOFHad,
                            pdgSigma,
                            pdgSigmaDau,
                            pdgHad,
                            sigmaPtGen,
                            hadPtGen,
                            kStarGen);
        }
      }
    }
  }
  PROCESS_SWITCH(sigmaHadCorrTask, processSameEventMC, "Process Same event MC", false);

  void processMixedEventMC(const CollisionsFullMC& collisions, const aod::KinkCands& kinkCands, const TracksFullMC& tracks, const aod::McParticles&)
  {
    for (auto const& [collision1, collision2] :
         selfCombinations(colBinning, nEvtMixingBkg, -1, collisions, collisions)) {
      if (collision1.index() == collision2.index())
        continue;

      sigmaHadCandidates.clear();
      if (std::abs(collision1.posZ()) > cutzvertex || !collision1.sel8()) {
        continue;
      }
      if (std::abs(collision2.posZ()) > cutzvertex || !collision2.sel8()) {
        continue;
      }
      auto kinkCands_c1 = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision1.globalIndex());
      auto tracks_c1 = tracks.sliceBy(tracksPerCollisionPreslice, collision1.globalIndex());
      auto tracks_c2 = tracks.sliceBy(tracksPerCollisionPreslice, collision2.globalIndex());
      fillTreeAndHistograms(kinkCands_c1, tracks, tracks_c2, collision1, true);

      for (const auto& candidate : sigmaHadCandidates) {
        auto mcLabelSigma = tracks.rawIteratorAt(candidate.sigmaID);
        auto mcLabelSigmaDau = tracks.rawIteratorAt(candidate.kinkDauID);
        auto mcLabelHad = tracks.rawIteratorAt(candidate.hadID);

        if (!mcLabelSigma.has_mcParticle() || !mcLabelSigmaDau.has_mcParticle() || !mcLabelHad.has_mcParticle()) {
          continue; // Skip candidates where MC truth is not available
        }

        auto mcPartSigma = mcLabelSigma.mcParticle_as<aod::McParticles>();
        auto mcPartSigmaDau = mcLabelSigmaDau.mcParticle_as<aod::McParticles>();
        auto mcPartHad = mcLabelHad.mcParticle_as<aod::McParticles>();
        auto pdgSigma = mcPartSigma.pdgCode();
        auto pdgSigmaDau = mcLabelSigmaDau.has_mcParticle() ? mcPartSigmaDau.pdgCode() : -999;
        auto pdgHad = mcLabelHad.has_mcParticle() ? mcPartHad.pdgCode() : -999;
        float sigmaPtGen = std::hypot(mcPartSigma.px(), mcPartSigma.py());
        float hadPtGen = std::hypot(mcPartHad.px(), mcPartHad.py());
        float kStarGen = getKStar(mcPartSigma.px(), mcPartSigma.py(), mcPartSigma.pz(), mcPartHad.px(), mcPartHad.py(), mcPartHad.pz());

        if (fillSparseInvMassKstar) {
          auto sigmaMomForKstar = getSigmaMomentumForKstar(candidate.sigmaPx, candidate.sigmaPy, candidate.sigmaPz,
                                                           candidate.sigmaDauPx, candidate.sigmaDauPy, candidate.sigmaDauPz);
          float kStarRec = getKStar(sigmaMomForKstar[0], sigmaMomForKstar[1], sigmaMomForKstar[2], candidate.pxHad, candidate.pyHad, candidate.pzHad);
          float sigmaPtUsed = std::hypot(sigmaMomForKstar[0], sigmaMomForKstar[1]);
          rSigmaHad.fill(HIST("hSparseSigmaHadMC"),
                         candidate.sigmaMass,
                         kStarRec,
                         candidate.sigmaCharge,
                         candidate.chargeHad,
                         candidate.sigmaDecRadius,
                         candidate.sigmaCosPA,
                         sigmaPtUsed,
                         kStarGen);
        }

        if (fillOutputTree) {
          outputDataTableMC(candidate.sigmaCharge,
                            candidate.sigmaPx,
                            candidate.sigmaPy,
                            candidate.sigmaPz,
                            candidate.sigmaDauPx,
                            candidate.sigmaDauPy,
                            candidate.sigmaDauPz,
                            candidate.sigmaDecRadius,
                            candidate.sigmaCosPA,
                            candidate.chargeHad,
                            candidate.pxHad,
                            candidate.pyHad,
                            candidate.pzHad,
                            candidate.nSigmaTPCHad,
                            candidate.nSigmaTOFHad,
                            pdgSigma,
                            pdgSigmaDau,
                            pdgHad,
                            sigmaPtGen,
                            hadPtGen,
                            kStarGen);
        }
      }
    }
    LOG(debug) << "Processing mixed event MC";
  }
  PROCESS_SWITCH(sigmaHadCorrTask, processMixedEventMC, "Process Mixed event MC", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<sigmaHadCorrTask>(cfgc)};
}
