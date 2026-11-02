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
///
/// \file antinucleiInJetsCDF.cxx
///
/// \brief Analysis task for antinuclei in jets using the CDF technique
/// \author Alberto Caliva (alberto.caliva@cern.ch)
/// \since August 22, 2024

#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>
#include <TVector3.h>

#include <cmath>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

// Define convenient aliases for commonly used table joins
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using RecCollisionsMc = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels>;
using GenCollisionsMc = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
using AntiNucleiTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPr, aod::pidTOFFullDe>;
using AntiNucleiTracksMc = soa::Join<AntiNucleiTracks, aod::McTrackLabels>;

// Lightweight particle container
struct ReducedParticle {
  double px;
  double py;
  double pz;
  double y;
  double nsigmaTPC;
  double nsigmaTOF;
  bool hasTOF;
  bool isPhysPrim;

  // Transverse Momentum
  double pt() const
  {
    return std::sqrt(px * px + py * py);
  }
};

struct AntinucleiInJetsCDF {

  // Histogram registries for real-data and MC analyses
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Event selection criteria
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "Reject events near the ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "Reject events near the TF border"};
  Configurable<bool> requireVtxITSTPC{"requireVtxITSTPC", true, "Require at least one ITS-TPC matched track"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "Reject events with same-bunch pileup collisions"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "Require consistent FT0 vs PV z-vertex"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "Require vtx track matched to TOF"};

  // Global Parameters
  Configurable<double> ptLeadingMin{"ptLeadingMin", 5.0, "pt Leading Min"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};

  // Track quality, kinematic, and PID selection parameters
  Configurable<bool> requirePvContributor{"requirePvContributor", false, "require that the track is a PV contributor"};
  Configurable<int> minItsNclusters{"minItsNclusters", 5, "minimum number of ITS clusters"};
  Configurable<int> minTpcNcrossedRows{"minTpcNcrossedRows", 100, "minimum number of TPC crossed pad rows"};
  Configurable<double> minChiSquareTpc{"minChiSquareTpc", 0.0, "minimum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0, "maximum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0, "maximum ITS chi^2/Ncls"};
  Configurable<double> minPt{"minPt", 0.3, "minimum pt of the tracks"};
  Configurable<double> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<double> maxEta{"maxEta", +0.8, "maximum eta"};
  Configurable<double> minY{"minY", -0.5, "minimum rapidity"};
  Configurable<double> maxY{"maxY", +0.5, "maximum rapidity"};
  Configurable<double> maxDcaxy{"maxDcaxy", 0.05, "Maximum DCAxy"};
  Configurable<double> maxDcaz{"maxDcaz", 0.05, "Maximum DCAz"};
  Configurable<double> minNsigmaTpc{"minNsigmaTpc", -3.0, "Minimum nsigma TPC"};
  Configurable<double> maxNsigmaTpc{"maxNsigmaTpc", +3.0, "Maximum nsigma TPC"};
  Configurable<double> minNsigmaTof{"minNsigmaTof", -3.0, "Minimum nsigma TOF"};
  Configurable<double> maxNsigmaTof{"maxNsigmaTof", +3.5, "Maximum nsigma TOF"};

  void init(InitContext const&)
  {
    // Binning
    double min = 0.0;
    double max = 6.0;
    int nbins = 120;

    // Process real data
    if (doprocessData) {

      // Event counter
      registryData.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{20, 0, 20, "counter"}});

      // Antiproton histograms
      registryData.add("antip_toward_tpc", "antip_toward_tpc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antip_toward_tof", "antip_toward_tof", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
      registryData.add("antip_transv_tpc", "antip_transv_tpc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antip_transv_tof", "antip_transv_tof", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});

      // Antideuteron histograms
      registryData.add("antid_toward_tpc", "antid_toward_tpc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antid_toward_tof", "antid_toward_tof", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
      registryData.add("antid_transv_tpc", "antid_transv_tpc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antid_transv_tof", "antid_transv_tof", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    }

    // Process generated MC
    if (doprocessGenMC) {

      // Event counter
      registryMC.add("genEvents", "number of generated events in mc", HistType::kTH1F, {{10, 0, 10, "counter"}});

      // Histograms for generated particles in the full event
      registryMC.add("antip_gen", "antip_gen", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antid_gen", "antid_gen", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      // Histograms for generated particles in different azimuthal regions
      registryMC.add("antip_toward_gen", "antip_toward_gen", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_transv_gen", "antip_transv_gen", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antid_toward_gen", "antid_toward_gen", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antid_transv_gen", "antid_transv_gen", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    }

    // Process reconstructed MC
    if (doprocessRecMC) {

      // Event counter
      registryMC.add("recEvents", "number of reconstructed events in mc", HistType::kTH1F, {{10, 0, 10, "counter"}});

      // Histograms for reconstructed particles in the full event
      registryMC.add("antip_rec_tpc", "antip_rec_tpc", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_rec_tof", "antip_rec_tof", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antid_rec_tpc", "antid_rec_tpc", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antid_rec_tof", "antid_rec_tof", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      // Histograms for reconstructed particles in different azimuthal regions
      registryMC.add("antip_toward_rec_tpc", "antip_toward_rec_tpc", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_toward_rec_tof", "antip_toward_rec_tof", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_transv_rec_tpc", "antip_transv_rec_tpc", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_transv_rec_tof", "antip_transv_rec_tof", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antid_toward_rec_tpc", "antid_toward_rec_tpc", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antid_toward_rec_tof", "antid_toward_rec_tof", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antid_transv_rec_tpc", "antid_transv_rec_tpc", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antid_transv_rec_tof", "antid_transv_rec_tof", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      // Histograms for secondary antiprotons
      registryMC.add("antip_toward_prim", "antip_toward_prim", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_toward_all", "antip_toward_all", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_transv_prim", "antip_transv_prim", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_transv_all", "antip_transv_all", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    }
  }

  // ITS hit
  template <typename TrackIts>
  bool hasITSHit(const TrackIts& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  // Single-Track Selection for leading track
  template <typename JetTrack>
  bool passedTrackSelectionLeadingTrack(const JetTrack& track)
  {
    static constexpr int MinTpcCr = 70;
    static constexpr double MaxChi2Tpc = 4.0;
    static constexpr double MaxChi2Its = 36.0;
    static constexpr double MaxPseudorap = 0.8;
    static constexpr double MinPtTrack = 0.3;
    static constexpr double DcaxyMaxTrackPar0 = 0.0105;
    static constexpr double DcaxyMaxTrackPar1 = 0.035;
    static constexpr double DcaxyMaxTrackPar2 = 1.1;
    static constexpr double DcazMaxTrack = 2.0;

    if (!track.hasITS())
      return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < MinTpcCr)
      return false;
    if (track.tpcChi2NCl() > MaxChi2Tpc)
      return false;
    if (track.itsChi2NCl() > MaxChi2Its)
      return false;
    if (std::abs(track.eta()) > MaxPseudorap)
      return false;
    if (track.pt() < MinPtTrack)
      return false;
    if (std::abs(track.dcaXY()) > (DcaxyMaxTrackPar0 + DcaxyMaxTrackPar1 / std::pow(track.pt(), DcaxyMaxTrackPar2)))
      return false;
    if (std::abs(track.dcaZ()) > DcazMaxTrack)
      return false;
    return true;
  }

  // Single-track selection for antinuclei
  template <typename AntinucleusTrack>
  bool passedTrackSelection(const AntinucleusTrack& track)
  {
    if (track.sign() >= 0)
      return false;
    if (requirePvContributor && !(track.isPVContributor()))
      return false;
    if (!track.hasITS())
      return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    if (track.itsNCls() < minItsNclusters)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minTpcNcrossedRows)
      return false;
    if (track.tpcChi2NCl() < minChiSquareTpc)
      return false;
    if (track.tpcChi2NCl() > maxChiSquareTpc)
      return false;
    if (track.itsChi2NCl() > maxChiSquareIts)
      return false;
    if (track.eta() < minEta || track.eta() > maxEta)
      return false;
    if (track.pt() < minPt)
      return false;

    return true;
  }

  // Rapidity
  double getRapidity(double px, double py, double pz, double mass)
  {
    const double energy = std::sqrt(px * px + py * py + pz * pz + mass * mass);
    ROOT::Math::PxPyPzEVector lorentzVect(px, py, pz, energy);
    return lorentzVect.Rapidity();
  }

  // Check if a particle is in the toward region with respect to the leading particle
  bool isParticleInTowardRegion(const TVector3& particle, const TVector3& leadingParticle)
  {
    const double deltaPhi = std::remainder(particle.Phi() - leadingParticle.Phi(), o2::constants::math::TwoPI);
    return std::abs(deltaPhi) < o2::constants::math::PIThird;
  }

  // Check if a particle is in the transverse region with respect to the leading particle
  bool isParticleInTransverseRegion(const TVector3& particle, const TVector3& leadingParticle)
  {
    static constexpr double Two = 2.0;
    const double deltaPhi = std::remainder(particle.Phi() - leadingParticle.Phi(), o2::constants::math::TwoPI);
    return (std::abs(deltaPhi) >= o2::constants::math::PIThird && std::abs(deltaPhi) < Two * o2::constants::math::PIThird);
  }

  // Check if particle is a physical primary or a decay product of a heavy-flavor hadron
  bool isPhysicalPrimaryOrFromHF(aod::McParticle const& particle, aod::McParticles const& mcParticles)
  {
    // Keep only pi, K, p, d, e, mu
    int pdg = std::abs(particle.pdgCode());
    if (!(pdg == PDG_t::kPiPlus || pdg == PDG_t::kKPlus || pdg == PDG_t::kProton || pdg == o2::constants::physics::Pdg::kDeuteron || pdg == PDG_t::kElectron || pdg == PDG_t::kMuonMinus))
      return false;

    // Constants for identifying heavy-flavor (charm and bottom) content from PDG codes
    static constexpr int CharmQuark = 4;
    static constexpr int BottomQuark = 5;
    static constexpr int Hundreds = 100;
    static constexpr int Thousands = 1000;

    // Check if particle is from heavy-flavor decay
    bool fromHF = false;
    if (particle.has_mothers()) {
      auto mother = mcParticles.iteratorAt(particle.mothersIds()[0]);
      int motherPdg = std::abs(mother.pdgCode());
      fromHF = (motherPdg / Hundreds == CharmQuark || motherPdg / Hundreds == BottomQuark || motherPdg / Thousands == CharmQuark || motherPdg / Thousands == BottomQuark);
    }

    // Select only physical primary particles or from heavy-flavor
    return (particle.isPhysicalPrimary() || fromHF);
  }

  // Process data
  void processData(SelectedCollisions::iterator const& collision, AntiNucleiTracks const& tracks)
  {
    // Event counter: before event selection
    registryData.fill(HIST("number_of_events_data"), 0.5);

    // Standard event selection
    if (!collision.sel8() || std::abs(collision.posZ()) > zVtx)
      return;
    registryData.fill(HIST("number_of_events_data"), 1.5);

    // Reject events near the ITS Read-Out Frame border
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return;
    registryData.fill(HIST("number_of_events_data"), 2.5);

    // Reject events at the time frame border
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);

    // Require at least one ITS-TPC matched track
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return;
    registryData.fill(HIST("number_of_events_data"), 4.5);

    // Reject events with same-bunch pileup
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    registryData.fill(HIST("number_of_events_data"), 5.5);

    // Require consistent FT0 vs PV z-vertex
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    registryData.fill(HIST("number_of_events_data"), 6.5);

    // Require TOF match for at least one vertex track
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      return;
    registryData.fill(HIST("number_of_events_data"), 7.5);

    // Leading track pt and momentum
    double ptMax(0);
    TVector3 leadingTrackVec;

    // Loop over reconstructed tracks
    for (auto const& track : tracks) {
      if (passedTrackSelectionLeadingTrack(track) && track.pt() > ptMax) {
        leadingTrackVec.SetXYZ(track.px(), track.py(), track.pz());
        ptMax = track.pt();
      }
    }

    // Event counter: skip events with pt < pt_leading_min
    if (ptMax < ptLeadingMin)
      return;
    registryData.fill(HIST("number_of_events_data"), 8.5);

    // Loop over reconstructed tracks
    for (auto const& track : tracks) {

      // Apply track and DCA selections
      if (!passedTrackSelection(track))
        continue;
      if (std::abs(track.dcaXY()) > maxDcaxy)
        continue;
      if (std::abs(track.dcaZ()) > maxDcaz)
        continue;

      // Variables
      const double nsigmaTPCPr = track.tpcNSigmaPr();
      const double nsigmaTOFPr = track.tofNSigmaPr();
      const double nsigmaTPCDe = track.tpcNSigmaDe();
      const double nsigmaTOFDe = track.tofNSigmaDe();
      const double yp = getRapidity(track.px(), track.py(), track.pz(), o2::constants::physics::MassProton);
      const double yd = getRapidity(track.px(), track.py(), track.pz(), o2::constants::physics::MassDeuteron);
      const TVector3 trackVec(track.px(), track.py(), track.pz());

      // Find azimuthal region
      const bool isToward = isParticleInTowardRegion(trackVec, leadingTrackVec);
      const bool isTransv = isParticleInTransverseRegion(trackVec, leadingTrackVec);

      // Antiprotons in the toward region
      if (isToward && yp >= minY && yp <= maxY) {
        registryData.fill(HIST("antip_toward_tpc"), track.pt(), nsigmaTPCPr);
        if (track.hasTOF() && nsigmaTPCPr >= minNsigmaTpc && nsigmaTPCPr <= maxNsigmaTpc)
          registryData.fill(HIST("antip_toward_tof"), track.pt(), nsigmaTOFPr);
      }

      // Antiprotons in the transverse region
      if (isTransv && yp >= minY && yp <= maxY) {
        registryData.fill(HIST("antip_transv_tpc"), track.pt(), nsigmaTPCPr);
        if (track.hasTOF() && nsigmaTPCPr >= minNsigmaTpc && nsigmaTPCPr <= maxNsigmaTpc)
          registryData.fill(HIST("antip_transv_tof"), track.pt(), nsigmaTOFPr);
      }

      // Antideuterons in the toward region
      if (isToward && yd >= minY && yd <= maxY) {
        registryData.fill(HIST("antid_toward_tpc"), track.pt(), nsigmaTPCDe);
        if (track.hasTOF() && nsigmaTPCDe >= minNsigmaTpc && nsigmaTPCDe <= maxNsigmaTpc)
          registryData.fill(HIST("antid_toward_tof"), track.pt(), nsigmaTOFDe);
      }

      // Antideuterons in the transverse region
      if (isTransv && yd >= minY && yd <= maxY) {
        registryData.fill(HIST("antid_transv_tpc"), track.pt(), nsigmaTPCDe);
        if (track.hasTOF() && nsigmaTPCDe >= minNsigmaTpc && nsigmaTPCDe <= maxNsigmaTpc)
          registryData.fill(HIST("antid_transv_tof"), track.pt(), nsigmaTOFDe);
      }
    }
  }
  PROCESS_SWITCH(AntinucleiInJetsCDF, processData, "Process real-data analysis", true);

  // Define preslice to group MC particles by their associated MC collision
  Preslice<aod::McParticles> mcParticlesPerMcCollision = o2::aod::mcparticle::mcCollisionId;

  // Process generated MC
  void processGenMC(GenCollisionsMc const& collisions, aod::McParticles const& mcParticles)
  {
    // Define per-event particle containers
    std::vector<TVector3> antip;
    std::vector<TVector3> antid;

    // Loop over generated collisions
    for (const auto& collision : collisions) {

      // Clear containers at the start of the event loop
      antip.clear();
      antid.clear();

      // Event counter: before event selection
      registryMC.fill(HIST("genEvents"), 0.5);

      // Apply event selection: require vertex position to be within the allowed z range
      if (std::abs(collision.posZ()) > zVtx)
        continue;
      registryMC.fill(HIST("genEvents"), 1.5);

      // Get particles in this MC collision
      const auto mcParticlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, collision.globalIndex());

      // Leading particle momentum
      double ptMax(0);
      TVector3 leadingParticleVec;

      // Loop over MC particles
      for (const auto& particle : mcParticlesThisMcColl) {

        // Select physical primary particles or HF decay products
        if (!isPhysicalPrimaryOrFromHF(particle, mcParticles))
          continue;

        // Select particles within acceptance
        if (particle.eta() < minEta || particle.eta() > maxEta)
          continue;

        // Store leading particle momentum
        if (particle.pt() > ptMax) {
          ptMax = particle.pt();
          leadingParticleVec.SetXYZ(particle.px(), particle.py(), particle.pz());
        }

        // Physical primary seletion
        if (!particle.isPhysicalPrimary())
          continue;

        // Rapidity selection
        if (particle.y() < minY || particle.y() > maxY)
          continue;

        // Momentum vector
        TVector3 pVec(particle.px(), particle.py(), particle.pz());

        // Fill histogram for generated antiprotons and store momentum
        if (particle.pdgCode() == PDG_t::kProtonBar) {
          registryMC.fill(HIST("antip_gen"), particle.pt());
          antip.emplace_back(pVec);
        }

        // Fill histogram for generated antideuterons and store momentum
        if (particle.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
          registryMC.fill(HIST("antid_gen"), particle.pt());
          antid.emplace_back(pVec);
        }
      }

      // Event Counter: skip events with pt < pt_leading_min
      if (ptMax < ptLeadingMin)
        continue;
      registryMC.fill(HIST("genEvents"), 2.5);

      // Loop over antiprotons
      for (const auto& antipVec : antip) {
        if (isParticleInTowardRegion(antipVec, leadingParticleVec))
          registryMC.fill(HIST("antip_toward_gen"), antipVec.Pt());
        if (isParticleInTransverseRegion(antipVec, leadingParticleVec))
          registryMC.fill(HIST("antip_transv_gen"), antipVec.Pt());
      }

      // Loop over antideuterons
      for (const auto& antidVec : antid) {
        if (isParticleInTowardRegion(antidVec, leadingParticleVec))
          registryMC.fill(HIST("antid_toward_gen"), antidVec.Pt());
        if (isParticleInTransverseRegion(antidVec, leadingParticleVec))
          registryMC.fill(HIST("antid_transv_gen"), antidVec.Pt());
      }
    }
  }
  PROCESS_SWITCH(AntinucleiInJetsCDF, processGenMC, "Process gen MC", false);

  // Define preslice to group reconstructed MC tracks by collision
  Preslice<AntiNucleiTracksMc> mcTracksPerMcCollision = o2::aod::track::collisionId;

  // Process reconstructed MC
  void processRecMC(RecCollisionsMc const& collisions, AntiNucleiTracksMc const& mcTracks, McParticles const&)
  {
    // Define per-event particle containers
    std::vector<ReducedParticle> antip;
    std::vector<ReducedParticle> antid;

    // Loop over all reconstructed collisions
    for (const auto& collision : collisions) {

      // Clear containers at the start of the event loop
      antip.clear();
      antid.clear();

      // Event counter: before event selection
      registryMC.fill(HIST("recEvents"), 0.5);

      // Standard event selection
      if (!collision.sel8() || std::abs(collision.posZ()) > zVtx)
        continue;
      registryMC.fill(HIST("recEvents"), 1.5);

      // Reject events near the ITS read-out frame border
      if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
        continue;
      registryMC.fill(HIST("recEvents"), 2.5);

      // Reject events at the time frame border
      if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
        continue;
      registryMC.fill(HIST("recEvents"), 3.5);

      // Require at least one ITS-TPC matched track
      if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
        continue;
      registryMC.fill(HIST("recEvents"), 4.5);

      // Reject events with same-bunch pileup
      if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
        continue;
      registryMC.fill(HIST("recEvents"), 5.5);

      // Require consistent FT0 vs PV z-vertex
      if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      registryMC.fill(HIST("recEvents"), 6.5);

      // Require TOF match for at least one vertex track
      if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
        continue;
      registryMC.fill(HIST("recEvents"), 7.5);

      // Get tracks in this MC collision
      const auto mcTracksThisMcColl = mcTracks.sliceBy(mcTracksPerMcCollision, collision.globalIndex());

      // Leading-track momentum
      TVector3 leadingTrackVec;
      double ptMax(0);

      // Loop over reconstructed tracks
      for (auto const& track : mcTracksThisMcColl) {

        // Select leading particle
        if (passedTrackSelectionLeadingTrack(track) && track.pt() > ptMax) {
          leadingTrackVec.SetXYZ(track.px(), track.py(), track.pz());
          ptMax = track.pt();
        }

        // Get corresponding MC particle
        if (!track.has_mcParticle())
          continue;
        const auto mcparticle = track.mcParticle();

        // Track and DCA selections
        if (!passedTrackSelection(track))
          continue;
        if (std::abs(track.dcaXY()) > maxDcaxy)
          continue;
        if (std::abs(track.dcaZ()) > maxDcaz)
          continue;

        // Variables
        const double nsigmaTPCPr = track.tpcNSigmaPr();
        const double nsigmaTOFPr = track.tofNSigmaPr();
        const double nsigmaTPCDe = track.tpcNSigmaDe();
        const double nsigmaTOFDe = track.tofNSigmaDe();
        const bool hasTOF = track.hasTOF();
        const bool isAntip = mcparticle.pdgCode() == PDG_t::kProtonBar;
        const bool isAntid = mcparticle.pdgCode() == -o2::constants::physics::Pdg::kDeuteron;
        const bool isPhysPrim = mcparticle.isPhysicalPrimary();
        const double yp = getRapidity(track.px(), track.py(), track.pz(), o2::constants::physics::MassProton);
        const double yd = getRapidity(track.px(), track.py(), track.pz(), o2::constants::physics::MassDeuteron);

        // Fill antiproton and antideuteron vectors
        if (isAntip)
          antip.push_back({track.px(), track.py(), track.pz(), yp, nsigmaTPCPr, nsigmaTOFPr, hasTOF, isPhysPrim});
        if (isAntid)
          antid.push_back({track.px(), track.py(), track.pz(), yd, nsigmaTPCDe, nsigmaTOFDe, hasTOF, isPhysPrim});

        // Selection of physical primary
        if (!isPhysPrim)
          continue;

        // Fill antiprotons for full-event efficiency
        if (isAntip && yp >= minY && yp <= maxY && nsigmaTPCPr >= minNsigmaTpc && nsigmaTPCPr <= maxNsigmaTpc) {
          registryMC.fill(HIST("antip_rec_tpc"), track.pt());
          if (hasTOF && nsigmaTOFPr >= minNsigmaTof && nsigmaTOFPr <= maxNsigmaTof) {
            registryMC.fill(HIST("antip_rec_tof"), track.pt());
          }
        }

        // Fill antideuterons for full-event efficiency
        if (isAntid && yd >= minY && yd <= maxY && nsigmaTPCDe >= minNsigmaTpc && nsigmaTPCDe <= maxNsigmaTpc) {
          registryMC.fill(HIST("antid_rec_tpc"), track.pt());
          if (hasTOF && nsigmaTOFDe >= minNsigmaTof && nsigmaTOFDe <= maxNsigmaTof) {
            registryMC.fill(HIST("antid_rec_tof"), track.pt());
          }
        }
      }

      // Event Counter: skip events with pt < pt_leading_min
      if (ptMax < ptLeadingMin)
        continue;
      registryMC.fill(HIST("recEvents"), 8.5);

      // Loop over antiprotons
      for (const auto& antipVec : antip) {

        TVector3 trackVec(antipVec.px, antipVec.py, antipVec.pz);
        const bool isToward = isParticleInTowardRegion(trackVec, leadingTrackVec);
        const bool isTransv = isParticleInTransverseRegion(trackVec, leadingTrackVec);

        // Antiprotons in Toward Region
        if (isToward && antipVec.y >= minY && antipVec.y <= maxY) {
          registryMC.fill(HIST("antip_toward_all"), antipVec.pt());
          if (antipVec.isPhysPrim)
            registryMC.fill(HIST("antip_toward_prim"), antipVec.pt());
        }

        // Antiprotons in transverse Region
        if (isTransv && antipVec.y >= minY && antipVec.y <= maxY) {
          registryMC.fill(HIST("antip_transv_all"), antipVec.pt());
          if (antipVec.isPhysPrim)
            registryMC.fill(HIST("antip_transv_prim"), antipVec.pt());
        }

        // Select physical primary
        if (!antipVec.isPhysPrim)
          continue;

        // Select antiprotons inside rapidity window
        if (antipVec.y < minY || antipVec.y > maxY)
          continue;

        // Antiprotons in Toward Region
        if (isToward && antipVec.nsigmaTPC >= minNsigmaTpc && antipVec.nsigmaTPC <= maxNsigmaTpc) {
          registryMC.fill(HIST("antip_toward_rec_tpc"), antipVec.pt());
          if (antipVec.hasTOF && antipVec.nsigmaTOF >= minNsigmaTof && antipVec.nsigmaTOF <= maxNsigmaTof)
            registryMC.fill(HIST("antip_toward_rec_tof"), antipVec.pt());
        }

        // Antiprotons in transverse Region
        if (isTransv && antipVec.nsigmaTPC >= minNsigmaTpc && antipVec.nsigmaTPC <= maxNsigmaTpc) {
          registryMC.fill(HIST("antip_transv_rec_tpc"), antipVec.pt());
          if (antipVec.hasTOF && antipVec.nsigmaTOF >= minNsigmaTof && antipVec.nsigmaTOF <= maxNsigmaTof)
            registryMC.fill(HIST("antip_transv_rec_tof"), antipVec.pt());
        }
      }

      // Loop over antideuterons
      for (const auto& antidVec : antid) {

        TVector3 trackVec(antidVec.px, antidVec.py, antidVec.pz);
        const bool isToward = isParticleInTowardRegion(trackVec, leadingTrackVec);
        const bool isTransv = isParticleInTransverseRegion(trackVec, leadingTrackVec);

        // Select physical primary
        if (!antidVec.isPhysPrim)
          continue;

        // Select antideuterons inside rapiditry window
        if (antidVec.y < minY || antidVec.y > maxY)
          continue;

        // Antideuterons in Toward Region
        if (isToward && antidVec.nsigmaTPC >= minNsigmaTpc && antidVec.nsigmaTPC <= maxNsigmaTpc) {
          registryMC.fill(HIST("antid_toward_rec_tpc"), antidVec.pt());
          if (antidVec.hasTOF && antidVec.nsigmaTOF >= minNsigmaTof && antidVec.nsigmaTOF <= maxNsigmaTof)
            registryMC.fill(HIST("antid_toward_rec_tof"), antidVec.pt());
        }

        // Antideuterons in transverse Region
        if (isTransv && antidVec.nsigmaTPC >= minNsigmaTpc && antidVec.nsigmaTPC <= maxNsigmaTpc) {
          registryMC.fill(HIST("antid_transv_rec_tpc"), antidVec.pt());
          if (antidVec.hasTOF && antidVec.nsigmaTOF >= minNsigmaTof && antidVec.nsigmaTOF <= maxNsigmaTof)
            registryMC.fill(HIST("antid_transv_rec_tof"), antidVec.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(AntinucleiInJetsCDF, processRecMC, "Process reconstructed MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AntinucleiInJetsCDF>(cfgc)};
}
