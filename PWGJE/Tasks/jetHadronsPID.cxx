// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include "Framework/runDataProcessing.h"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/PseudoJet.hh>

#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::constants::math;

// Define convenient aliases for commonly used table joins
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;


using PionTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCPi, aod::pidTOFPi>;

struct PIDHadronsInJets {

  // Histogram registry for data
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Parameters for ppRef analysis
  Configurable<bool> isppRefAnalysis{"isppRefAnalysis", false, "Is ppRef analysis"};
  Configurable<double> cfgAreaFrac{"cfgAreaFrac", 0.6, "fraction of jet area"};
  Configurable<double> cfgEtaJetMax{"cfgEtaJetMax", 0.5, "max jet eta"};
  Configurable<double> cfgMinPtTrack{"cfgMinPtTrack", 0.1, "minimum pt of tracks for jet reconstruction"};

  // Event selection criteria
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "Reject events near the ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "Reject events near the TF border"};
  Configurable<bool> requireVtxITSTPC{"requireVtxITSTPC", true, "Require at least one ITS-TPC matched track"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "Reject events with same-bunch pileup collisions"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "Require consistent FT0 vs PV z-vertex"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "Require at least one vertex track matched to TOF"};

  // Jet selection parameters
  Configurable<double> minJetPt{"minJetPt", 10.0, "Minimum pt of the jet after bkg subtraction"};
  Configurable<double> maxJetPt{"maxJetPt", 1e+06, "Maximum pt of the jet after bkg subtraction"};
  Configurable<double> rJet{"rJet", 0.4, "Jet resolution parameter R"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<bool> applyAreaCut{"applyAreaCut", true, "apply area cut"};
  Configurable<double> maxNormalizedJetArea{"maxNormalizedJetArea", 1.0, "area cut"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from the edge"};

  // Track quality parameters
  Configurable<bool> requirePvContributor{"requirePvContributor", false, "require that the track is a PV contributor"};
  Configurable<int> minItsNclusters{"minItsNclusters", 5, "minimum number of ITS clusters"};
  Configurable<int> minTpcNcrossedRows{"minTpcNcrossedRows", 100, "minimum number of TPC crossed pad rows"};
  Configurable<double> minChiSquareTpc{"minChiSquareTpc", 0.0, "minimum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0, "maximum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0, "maximum ITS chi^2/Ncls"};
  Configurable<double> minPt{"minPt", 0.3, "minimum pt of the tracks"};
  Configurable<double> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<double> maxEta{"maxEta", +0.8, "maximum eta"};
  Configurable<double> maxDcaxy{"maxDcaxy", 0.05, "Maximum DCAxy"};
  Configurable<double> maxDcaz{"maxDcaz", 0.05, "Maximum DCAz"};
  
  Configurable<bool> setMCDefaultItsParams{"setMCDefaultItsParams", true, "set MC default parameters"};

  // Utility object for jet background subtraction methods
  JetBkgSubUtils backgroundSub;

  // Initialize ITS PID Response object
  o2::aod::ITSResponse itsResponse;

  void init(InitContext const&)
  {
    if (setMCDefaultItsParams) {
      itsResponse.setMCDefaultParameters();
    }

    // pid of pions
    registryData.add("pion_jet_tpc", "TPC Pion PID in Jets", HistType::kTH2F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}, {100, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("pion_jet_tof", "TOF Pion PID in Jets", HistType::kTH2F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}, {100, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("pion_jet_its", "ITS Pion PID in Jets", HistType::kTH2F, {{120, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"}, {100, -3.0, 3.0, "n#sigma_{ITS}"}});
  }

  // ITS hit helper
  template <typename TrackIts>
  bool hasITSHit(const TrackIts& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  // Single-track selection for jet reconstruction
  template <typename JetTrack>
  bool passedTrackSelectionForJetReconstruction(const JetTrack& track)
  {
    static constexpr int MinTpcCr = 70;
    static constexpr double MaxChi2Tpc = 4.0;
    static constexpr double MaxChi2Its = 36.0;
    static constexpr double DcaxyMaxTrackPar0 = 0.0105;
    static constexpr double DcaxyMaxTrackPar1 = 0.035;
    static constexpr double DcaxyMaxTrackPar2 = 1.1;
    static constexpr double DcazMaxTrack = 2.0;

    if (!track.hasITS() || !track.hasTPC()) return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3))) return false;
    if (track.tpcNClsCrossedRows() < MinTpcCr) return false;
    if (track.tpcChi2NCl() > MaxChi2Tpc) return false;
    if (track.itsChi2NCl() > MaxChi2Its) return false;
    if (std::fabs(track.eta()) > maxEta) return false;
    if (track.pt() < cfgMinPtTrack) return false;
    if (std::fabs(track.dcaXY()) > (DcaxyMaxTrackPar0 + DcaxyMaxTrackPar1 / std::pow(track.pt(), DcaxyMaxTrackPar2))) return false;
    if (std::fabs(track.dcaZ()) > DcazMaxTrack) return false;
    return true;
  }

  // Single-track selection for constituents
  template <typename PionTrack>
  bool passedTrackSelection(const PionTrack& track)
  {
    if (requirePvContributor && !(track.isPVContributor())) return false;
    if (!track.hasITS() || !track.hasTPC()) return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3))) return false;
    if (track.itsNCls() < minItsNclusters) return false;
    if (track.tpcNClsCrossedRows() < minTpcNcrossedRows) return false;
    if (track.tpcChi2NCl() < minChiSquareTpc || track.tpcChi2NCl() > maxChiSquareTpc) return false;
    if (track.itsChi2NCl() > maxChiSquareIts) return false;
    if (track.eta() < minEta || track.eta() > maxEta) return false;
    if (track.pt() < minPt) return false;
    return true;
  }

  // Process Data
  void process(SelectedCollisions::iterator const& collision, PionTracks const& tracks)
  {
    // Apply standard event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx) return;
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) return;
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) return;
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) return;
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) return;
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) return;
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) return;

    // Build FastJet particles
    int id(-1);
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {
      id++;
      if (!passedTrackSelectionForJetReconstruction(track)) continue;

      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
      fourMomentum.set_user_index(id);
      fjParticles.emplace_back(fourMomentum);
    }

    if (fjParticles.empty()) return;

    // Cluster particles
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

    // Loop over reconstructed jets
    for (const auto& jet : jets) {

      if (!isppRefAnalysis && ((std::fabs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))) continue;
      if (isppRefAnalysis && std::fabs(jet.eta()) > cfgEtaJetMax) continue;

      auto jetForSub = jet;
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
      
      if (isppRefAnalysis && (jet.pt() < minJetPt || jet.pt() > maxJetPt)) continue;
      if (!isppRefAnalysis && (jetMinusBkg.pt() < minJetPt || jetMinusBkg.pt() > maxJetPt)) continue;

      double normalizedJetArea = jet.area() / (PI * rJet * rJet);
      if (applyAreaCut && (!isppRefAnalysis) && normalizedJetArea > maxNormalizedJetArea) continue;
      if (isppRefAnalysis && (jet.area() < cfgAreaFrac * PI * rJet * rJet)) continue;

      std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();

      // loop to fill historgrams
      for (const auto& particle : jetConstituents) {

        auto const& track = tracks.iteratorAt(particle.user_index());
        
        // Constituent Track Selection (includes DCA checks)
        if (!passedTrackSelection(track)) continue;
        if (std::fabs(track.dcaXY()) > maxDcaxy || std::fabs(track.dcaZ()) > maxDcaz) continue;

        double pt = track.pt();
        // double nSigmaITSPi = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Pion>(track));


        // Check Pion PID (+/- 3 Sigma)
        double nsigmaTPCPi = track.tpcNSigmaPi();
        
        // Fill TPC
        if (std::abs(nsigmaTPCPi) <= 3.0) {
            registryData.fill(HIST("pion_jet_tpc"), pt, nsigmaTPCPi);
        }
        
        // Fill TOF
        if (track.hasTOF()) {
            double nsigmaTOFPi = track.tofNSigmaPi();
            if (std::abs(nsigmaTOFPi) <= 3.0) {
                registryData.fill(HIST("pion_jet_tof"), pt, nsigmaTOFPi);
            }
        }

        // // Fill ITS
        // if (std::abs(nSigmaITSPi) <= 3.0) {
        //     registryData.fill(HIST("pion_jet_its"), pt, nSigmaITSPi);
        // }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PIDHadronsInJets>(cfgc)};
}