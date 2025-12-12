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

/// \file   sigmaProtonCorr.cxx
/// \brief Analysis task for sigma-proton correlations
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>

#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFSigmaProtonTables.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullPr>;
using TracksFullMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullPr, aod::McTrackLabels>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

struct sigmaProtonCand {

  float ptPr() const
  {
    return std::hypot(pxPr, pyPr);
  }
  float sigmaPt() const
  {
    return std::hypot(sigmaPx, sigmaPy);
  }

  int sigmaCharge = 0;       // Charge of the sigma candidate
  int sigmaMass = -1;        // Mass of the Sigma candidate
  float sigmaPx = -1;        // Px of the Sigma candidate
  float sigmaPy = -1;        // Py of the Sigma candidate
  float sigmaPz = -1;        // Pz of the Sigma candidate
  float sigmaDauPx = -1;     // Px of the daughter track from Sigma decay
  float sigmaDauPy = -1;     // Py of the daughter track from Sigma decay
  float sigmaDauPz = -1;     // Pz of the daughter track from Sigma decay
  float sigmaDecRadius = -1; // Decay radius of the Sigma candidate

  int chargePr = 0; // Charge of the proton candidate
  float pxPr = -1;  // Px of the proton candidate
  float pyPr = -1;  // Py of the proton candidate
  float pzPr = -1;  // Pz of the proton candidate

  float nSigmaTPCPr = -1; // Number of sigmas for the proton candidate
  float nSigmaTOFPr = -1; // Number of sigmas for the proton candidate using TOF

  int kinkDauID = -1; // ID of the pion from Sigma decay in MC
  int sigmaID = -1;   // ID of the Sigma candidate in MC
  int prID = -1;      // ID of the proton candidate in MC
};

struct sigmaProtonCorrTask {

  std::vector<sigmaProtonCand> sigmaProtonCandidates;  // Vector to store Sigma-Proton candidates
  Produces<aod::SigmaProtonCands> outputDataTable;     // Output table for Sigma-Proton candidates
  Produces<aod::SigmaProtonMCCands> outputDataTableMC; // Output table for Sigma-Proton candidates in MC
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rSigmaProton{"sigmaProton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // Configurable for event selection
  Configurable<float> cutzvertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutEtaDaught{"cutEtaDaughter", 0.8f, "Eta cut for daughter tracks"};
  Configurable<float> cutDCAtoPVSigma{"cutDCAtoPVSigma", 0.1f, "Max DCA to primary vertex for Sigma candidates (cm)"};
  Configurable<bool> doSigmaMinus{"doSigmaMinus", true, "If true, pair Sigma- candidates, else Sigma+"};
  Configurable<float> cutSigmaRadius{"cutSigmaRadius", 20.f, "Minimum radius for Sigma candidates (cm)"};
  Configurable<float> cutSigmaMass{"cutSigmaMass", 0.3, "Sigma mass window (MeV/c^2)"};
  Configurable<float> alphaAPCut{"alphaAPCut", 0., "Alpha AP cut for Sigma candidates"};
  Configurable<float> qtAPCutLow{"qtAPCutLow", 0.15, "Lower qT AP cut for Sigma candidates (GeV/c)"};
  Configurable<float> qtAPCutHigh{"qtAPCutHigh", 0.2, "Upper qT AP cut for Sigma candidates (GeV/c)"};

  Configurable<float> cutNITSClusPr{"cutNITSClusPr", 5, "Minimum number of ITS clusters for proton candidate"};
  Configurable<float> cutNTPCClusPr{"cutNTPCClusPr", 90, "Minimum number of TPC clusters for proton candidate"};
  Configurable<float> cutNSigmaTPC{"cutNSigmaTPC", 3, "NSigmaTPCPr"};
  Configurable<float> cutNSigmaTOF{"cutNSigmaTOF", 3, "NSigmaTOFPr"};

  Configurable<bool> fillOutputTree{"fillOutputTree", true, "If true, fill the output tree with Sigma-Proton candidates"};

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
    const AxisSpec nSigmaPrAxis{100, -5, 5, "n#sigma_{#pr}"};
    const AxisSpec sigmaMassAxis{100, 1.1, 1.4, "m (GeV/#it{c}^{2})"};
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};

    // qa histograms
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rSigmaProton.add("h2PtMassSigma", "h2PtMassSigma", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rSigmaProton.add("h2PtPrNSigma", "h2PtPrNSigma", {HistType::kTH2F, {ptAxis, nSigmaPrAxis}});
    rSigmaProton.add("h2PtPrNSigmaTOF", "h2PtPrNSigmaTOF", {HistType::kTH2F, {ptAxis, nSigmaPrAxis}});

    LOG(info) << "Sigma-Proton correlation task initialized";
    LOG(info) << "Process SE enabled: " << doprocessSameEvent;
    LOG(info) << "Process ME enabled: " << doprocessMixedEvent;
    LOG(info) << "Process SE MC enabled: " << doprocessSameEventMC;
    LOG(info) << "Process ME MC enabled: " << doprocessMixedEventMC;
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

  template <typename Ttrack>
  bool selectPrTrack(const Ttrack& candidate)
  {
    if (std::abs(candidate.tpcNSigmaPr()) > cutNSigmaTPC || candidate.tpcNClsFound() < cutNTPCClusPr || std::abs(candidate.eta()) > cutEtaDaught) {
      return false;
    }

    if (candidate.itsNCls() < cutNITSClusPr) {
      return false;
    }

    float ptMinTOF = 0.75f; // Minimum pT to use TOF for proton PID
    if (candidate.pt() < ptMinTOF) {
      return true; // No TOF cut for low pT
    }

    if (!candidate.hasTOF()) {
      return false;
    }
    if (std::abs(candidate.tofNSigmaPr()) > cutNSigmaTOF) {
      return false;
    }
    return true; // Track is selected
  }

  template <typename Ttrack>
  bool selectSigma(aod::KinkCands::iterator const& sigmaCand, Ttrack const&)
  {

    auto kinkDauTrack = sigmaCand.trackDaug_as<Ttrack>();
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
    if (std::abs(sigmaCand.dcaMothPv()) > cutDCAtoPVSigma) {
      return false;
    }
    return true;
  }

  template <typename Ttrack>
  void fillTreeAndHistograms(aod::KinkCands const& kinkCands, Ttrack const& tracks, Ttrack const& tracksDauSigma)
  {
    for (const auto& sigmaCand : kinkCands) {
      if (selectSigma(sigmaCand, tracksDauSigma)) {
        if (doSigmaMinus) {
          rSigmaProton.fill(HIST("h2PtMassSigma"), sigmaCand.mothSign() * sigmaCand.ptMoth(), sigmaCand.mSigmaMinus());
        } else {
          rSigmaProton.fill(HIST("h2PtMassSigma"), sigmaCand.mothSign() * sigmaCand.ptMoth(), sigmaCand.mSigmaPlus());
        }

        for (const auto& prTrack : tracks) {
          if (!selectPrTrack(prTrack)) {
            continue;
          }

          sigmaProtonCand candidate;
          candidate.sigmaCharge = sigmaCand.mothSign();
          candidate.sigmaPx = sigmaCand.pxMoth();
          candidate.sigmaPy = sigmaCand.pyMoth();
          candidate.sigmaPz = sigmaCand.pzMoth();
          candidate.sigmaDauPx = sigmaCand.pxDaug();
          candidate.sigmaDauPy = sigmaCand.pyDaug();
          candidate.sigmaDauPz = sigmaCand.pzDaug();
          candidate.sigmaDecRadius = std::hypot(sigmaCand.xDecVtx(), sigmaCand.yDecVtx());

          candidate.chargePr = prTrack.sign();
          candidate.pxPr = prTrack.px();
          candidate.pyPr = prTrack.py();
          candidate.pzPr = prTrack.pz();
          candidate.nSigmaTPCPr = prTrack.tpcNSigmaPr();
          candidate.nSigmaTOFPr = prTrack.tofNSigmaPr();
          candidate.sigmaMass = doSigmaMinus ? sigmaCand.mSigmaMinus() : sigmaCand.mSigmaPlus();

          candidate.sigmaID = sigmaCand.trackMothId();
          candidate.kinkDauID = sigmaCand.trackDaugId();
          candidate.prID = prTrack.globalIndex();

          rSigmaProton.fill(HIST("h2PtPrNSigma"), candidate.ptPr(), candidate.nSigmaTPCPr);
          if (prTrack.hasTOF()) {
            rSigmaProton.fill(HIST("h2PtPrNSigmaTOF"), candidate.ptPr(), candidate.nSigmaTOFPr);
          }
          sigmaProtonCandidates.push_back(candidate);
        }
      }
    }
  }

  void processSameEvent(CollisionsFull const& collisions, aod::KinkCands const& kinkCands, TracksFull const& tracks)
  {
    for (auto const& collision : collisions) {

      sigmaProtonCandidates.clear();
      auto kinkCands_c = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision.globalIndex());
      auto tracks_c = tracks.sliceBy(tracksPerCollisionPreslice, collision.globalIndex());
      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
        continue;
      }
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      fillTreeAndHistograms(kinkCands_c, tracks_c, tracks_c);
      if (fillOutputTree) {
        // Fill output table
        for (const auto& candidate : sigmaProtonCandidates) {
          outputDataTable(candidate.sigmaCharge,
                          candidate.sigmaPx,
                          candidate.sigmaPy,
                          candidate.sigmaPz,
                          candidate.sigmaDauPx,
                          candidate.sigmaDauPy,
                          candidate.sigmaDauPz,
                          candidate.sigmaDecRadius,
                          candidate.chargePr,
                          candidate.pxPr,
                          candidate.pyPr,
                          candidate.pzPr,
                          candidate.nSigmaTPCPr,
                          candidate.nSigmaTOFPr);
        }
      }
    }
  }
  PROCESS_SWITCH(sigmaProtonCorrTask, processSameEvent, "Process Same event", true);

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

      sigmaProtonCandidates.clear();
      if (std::abs(collision1.posZ()) > cutzvertex || !collision1.sel8()) {
        continue;
      }
      if (std::abs(collision2.posZ()) > cutzvertex || !collision2.sel8()) {
        continue;
      }
      auto kinkCands_c1 = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision1.globalIndex());
      auto tracks_c1 = tracks.sliceBy(tracksPerCollisionPreslice, collision1.globalIndex());
      auto tracks_c2 = tracks.sliceBy(tracksPerCollisionPreslice, collision2.globalIndex());
      fillTreeAndHistograms(kinkCands_c1, tracks_c1, tracks_c2);

      auto kinkCands_c2 = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision2.globalIndex());
      fillTreeAndHistograms(kinkCands_c2, tracks_c2, tracks_c1);

      if (fillOutputTree) {
        // Fill output table
        for (const auto& candidate : sigmaProtonCandidates) {
          outputDataTable(candidate.sigmaCharge,
                          candidate.sigmaPx,
                          candidate.sigmaPy,
                          candidate.sigmaPz,
                          candidate.sigmaDauPx,
                          candidate.sigmaDauPy,
                          candidate.sigmaDauPz,
                          candidate.sigmaDecRadius,
                          candidate.chargePr,
                          candidate.pxPr,
                          candidate.pyPr,
                          candidate.pzPr,
                          candidate.nSigmaTPCPr,
                          candidate.nSigmaTOFPr);
        }
      }
    }
    LOG(debug) << "Processing mixed event";
  }
  PROCESS_SWITCH(sigmaProtonCorrTask, processMixedEvent, "Process Mixed event", false);

  void processSameEventMC(CollisionsFullMC const& collisions, aod::KinkCands const& kinkCands, TracksFullMC const& tracks, aod::McParticles const&)
  {
    for (auto const& collision : collisions) {

      sigmaProtonCandidates.clear();
      auto kinkCands_c = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision.globalIndex());
      auto tracks_c = tracks.sliceBy(tracksMCPerCollisionPreslice, collision.globalIndex());

      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
        continue;
      }
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      fillTreeAndHistograms(kinkCands_c, tracks_c, tracks_c);
      if (fillOutputTree) {
        // Fill output table
        for (const auto& candidate : sigmaProtonCandidates) {
          auto mcLabelSigma = tracks.rawIteratorAt(candidate.sigmaID);
          auto mcLabelSigmaDau = tracks.rawIteratorAt(candidate.kinkDauID);
          auto mcLabelPr = tracks.rawIteratorAt(candidate.prID);
          auto pdgSigma = mcLabelSigma.has_mcParticle() ? mcLabelSigma.mcParticle_as<aod::McParticles>().pdgCode() : -999;
          auto pdgSigmaDau = mcLabelSigmaDau.has_mcParticle() ? mcLabelSigmaDau.mcParticle_as<aod::McParticles>().pdgCode() : -999;
          auto pdgPr = mcLabelPr.has_mcParticle() ? mcLabelPr.mcParticle_as<aod::McParticles>().pdgCode() : -999;
          outputDataTableMC(candidate.sigmaCharge,
                            candidate.sigmaPx,
                            candidate.sigmaPy,
                            candidate.sigmaPz,
                            candidate.sigmaDauPx,
                            candidate.sigmaDauPy,
                            candidate.sigmaDauPz,
                            candidate.sigmaDecRadius,
                            candidate.chargePr,
                            candidate.pxPr,
                            candidate.pyPr,
                            candidate.pzPr,
                            candidate.nSigmaTPCPr,
                            candidate.nSigmaTOFPr,
                            pdgSigma,
                            pdgSigmaDau,
                            pdgPr);
        }
      }
    }
  }
  PROCESS_SWITCH(sigmaProtonCorrTask, processSameEventMC, "Process Same event MC", false);

  void processMixedEventMC(const CollisionsFullMC& collisions, const aod::KinkCands& kinkCands, const TracksFullMC& tracks, const aod::McParticles&)
  {
    for (auto const& [collision1, collision2] :
         selfCombinations(colBinning, nEvtMixingBkg, -1, collisions, collisions)) {
      if (collision1.index() == collision2.index())
        continue;

      sigmaProtonCandidates.clear();
      if (std::abs(collision1.posZ()) > cutzvertex || !collision1.sel8()) {
        continue;
      }
      if (std::abs(collision2.posZ()) > cutzvertex || !collision2.sel8()) {
        continue;
      }
      auto kinkCands_c1 = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision1.globalIndex());
      auto tracks_c1 = tracks.sliceBy(tracksPerCollisionPreslice, collision1.globalIndex());
      auto tracks_c2 = tracks.sliceBy(tracksPerCollisionPreslice, collision2.globalIndex());
      fillTreeAndHistograms(kinkCands_c1, tracks_c1, tracks_c2);

      auto kinkCands_c2 = kinkCands.sliceBy(kinkCandsPerCollisionPreslice, collision2.globalIndex());
      fillTreeAndHistograms(kinkCands_c2, tracks_c2, tracks_c1);

      if (fillOutputTree) {
        // Fill output table
        for (const auto& candidate : sigmaProtonCandidates) {
          auto mcLabelSigma = tracks.rawIteratorAt(candidate.sigmaID);
          auto mcLabelSigmaDau = tracks.rawIteratorAt(candidate.kinkDauID);
          auto mcLabelPr = tracks.rawIteratorAt(candidate.prID);
          auto pdgSigma = mcLabelSigma.has_mcParticle() ? mcLabelSigma.mcParticle_as<aod::McParticles>().pdgCode() : -999;
          auto pdgSigmaDau = mcLabelSigmaDau.has_mcParticle() ? mcLabelSigmaDau.mcParticle_as<aod::McParticles>().pdgCode() : -999;
          auto pdgPr = mcLabelPr.has_mcParticle() ? mcLabelPr.mcParticle_as<aod::McParticles>().pdgCode() : -999;
          outputDataTableMC(candidate.sigmaCharge,
                            candidate.sigmaPx,
                            candidate.sigmaPy,
                            candidate.sigmaPz,
                            candidate.sigmaDauPx,
                            candidate.sigmaDauPy,
                            candidate.sigmaDauPz,
                            candidate.sigmaDecRadius,
                            candidate.chargePr,
                            candidate.pxPr,
                            candidate.pyPr,
                            candidate.pzPr,
                            candidate.nSigmaTPCPr,
                            candidate.nSigmaTOFPr,
                            pdgSigma,
                            pdgSigmaDau,
                            pdgPr);
        }
      }
    }
    LOG(debug) << "Processing mixed event MC";
  }
  PROCESS_SWITCH(sigmaProtonCorrTask, processMixedEventMC, "Process Mixed event MC", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<sigmaProtonCorrTask>(cfgc)};
}
