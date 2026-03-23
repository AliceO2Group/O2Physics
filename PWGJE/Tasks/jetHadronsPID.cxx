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

#include <TVector2.h>
#include <TVector3.h>

#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::constants::math;

// Define convenient aliases for commonly used table joins
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;


using HadronTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, 
                               aod::pidTPCFullPi, aod::pidTOFFullPi, 
                               aod::pidTPCFullKa, aod::pidTOFFullKa, 
                               aod::pidTPCFullPr, aod::pidTOFFullPr>;

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
  Configurable<double> maxPt{"maxPt", 4.0, "maximum pt of the tracks for PID analysis"}; // pt cut at 4.0
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

    // Pions, Kaons, Protons in Jet
    registryData.add("pion_jet_tpc", "TPC Pion PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("pion_jet_tof", "TOF Pion PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("kaon_jet_tpc", "TPC Kaon PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("kaon_jet_tof", "TOF Kaon PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("proton_jet_tpc", "TPC Proton PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("proton_jet_tof", "TOF Proton PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});

    // Pions, Kaons, Protons Underlying Event
    registryData.add("pion_ue_tpc", "TPC Pion PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("pion_ue_tof", "TOF Pion PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("kaon_ue_tpc", "TPC Kaon PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("kaon_ue_tof", "TOF Kaon PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("proton_ue_tpc", "TPC Proton PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("proton_ue_tof", "TOF Proton PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
  }

  // Compute two transverse directions orthogonal to vector p
  void getPerpendicularDirections(const TVector3& p, TVector3& u1, TVector3& u2)
  {
    double px = p.X(), py = p.Y(), pz = p.Z();
    double px2 = px * px, py2 = py * py, pz2 = pz * pz;
    double pz4 = pz2 * pz2;

    if (px == 0 && py == 0) { u1.SetXYZ(0, 0, 0); u2.SetXYZ(0, 0, 0); return; }
    if (px == 0 && py != 0) {
      double ux = std::sqrt(py2 - pz4 / py2);
      double uy = -pz2 / py;
      u1.SetXYZ(ux, uy, pz); u2.SetXYZ(-ux, uy, pz); return;
    }
    if (py == 0 && px != 0) {
      double ux = -pz2 / px;
      double uy = std::sqrt(px2 - pz4 / px2);
      u1.SetXYZ(ux, uy, pz); u2.SetXYZ(ux, -uy, pz); return;
    }

    double a = px2 + py2;
    double b = 2.0 * px * pz2;
    double c = pz4 - py2 * py2 - px2 * py2;
    double delta = b * b - 4.0 * a * c;

    if (delta < 0 || a == 0) { u1.SetXYZ(0, 0, 0); u2.SetXYZ(0, 0, 0); return; }
    double u1x = (-b + std::sqrt(delta)) / (2.0 * a);
    u1.SetXYZ(u1x, (-pz2 - px * u1x) / py, pz);
    double u2x = (-b - std::sqrt(delta)) / (2.0 * a);
    u2.SetXYZ(u2x, (-pz2 - px * u2x) / py, pz);
  }

  double getDeltaPhi(double a1, double a2)
  {
    double deltaPhi(0);
    double phi1 = TVector2::Phi_0_2pi(a1);
    double phi2 = TVector2::Phi_0_2pi(a2);
    double diff = std::fabs(phi1 - phi2);

    if (diff <= PI) deltaPhi = diff;
    if (diff > PI) deltaPhi = TwoPI - diff;
    return deltaPhi;
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
    if (track.pt() < minPt || track.pt() > maxPt) return false;
    return true;
  }

  // Process Data
  void process(SelectedCollisions::iterator const& collision, HadronTracks const& tracks)
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


    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    
    if (jets.empty()) return;

    auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

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


      double coneRadius = std::sqrt(jet.area() / PI);
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
      getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
      
      if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) continue;

      std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();
      for (const auto& particle : jetConstituents) {

        auto const& track = tracks.iteratorAt(particle.user_index());
        
        if (!passedTrackSelection(track)) continue;
        if (std::fabs(track.dcaXY()) > maxDcaxy || std::fabs(track.dcaZ()) > maxDcaz) continue;

        double pt = track.pt();

        if (std::abs(track.tpcNSigmaPi()) <= 3.0) {
            registryData.fill(HIST("pion_jet_tpc"), pt, track.tpcNSigmaPi());
        }
        if (track.hasTOF() && std::abs(track.tofNSigmaPi()) <= 3.0) {
            registryData.fill(HIST("pion_jet_tof"), pt, track.tofNSigmaPi());
        }

        if (std::abs(track.tpcNSigmaKa()) <= 3.0) {
            registryData.fill(HIST("kaon_jet_tpc"), pt, track.tpcNSigmaKa());
        }
        if (track.hasTOF() && std::abs(track.tofNSigmaKa()) <= 3.0) {
            registryData.fill(HIST("kaon_jet_tof"), pt, track.tofNSigmaKa());
        }

        if (std::abs(track.tpcNSigmaPr()) <= 3.0) {
            registryData.fill(HIST("proton_jet_tpc"), pt, track.tpcNSigmaPr());
        }
        if (track.hasTOF() && std::abs(track.tofNSigmaPr()) <= 3.0) {
            registryData.fill(HIST("proton_jet_tof"), pt, track.tofNSigmaPr());
        }
      }

      for (auto const& track : tracks) {

        if (!passedTrackSelection(track)) continue;
        if (std::fabs(track.dcaXY()) > maxDcaxy || std::fabs(track.dcaZ()) > maxDcaz) continue;

        double deltaEtaUe1 = track.eta() - ueAxis1.Eta();
        double deltaPhiUe1 = getDeltaPhi(track.phi(), ueAxis1.Phi());
        double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
        
        double deltaEtaUe2 = track.eta() - ueAxis2.Eta();
        double deltaPhiUe2 = getDeltaPhi(track.phi(), ueAxis2.Phi());
        double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

        double maxConeRadius = coneRadius;
        if (applyAreaCut) {
          maxConeRadius = std::sqrt(maxNormalizedJetArea) * rJet;
        }

        if (deltaRUe1 > maxConeRadius && deltaRUe2 > maxConeRadius) continue;

        double pt = track.pt();

        if (std::abs(track.tpcNSigmaPi()) <= 3.0) {
            registryData.fill(HIST("pion_ue_tpc"), pt, track.tpcNSigmaPi());
        }
        if (track.hasTOF() && std::abs(track.tofNSigmaPi()) <= 3.0) {
            registryData.fill(HIST("pion_ue_tof"), pt, track.tofNSigmaPi());
        }

        if (std::abs(track.tpcNSigmaKa()) <= 3.0) {
            registryData.fill(HIST("kaon_ue_tpc"), pt, track.tpcNSigmaKa());
        }
        if (track.hasTOF() && std::abs(track.tofNSigmaKa()) <= 3.0) {
            registryData.fill(HIST("kaon_ue_tof"), pt, track.tofNSigmaKa());
        }

        if (std::abs(track.tpcNSigmaPr()) <= 3.0) {
            registryData.fill(HIST("proton_ue_tpc"), pt, track.tpcNSigmaPr());
        }
        if (track.hasTOF() && std::abs(track.tofNSigmaPr()) <= 3.0) {
            registryData.fill(HIST("proton_ue_tof"), pt, track.tofNSigmaPr());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PIDHadronsInJets>(cfgc)};
}