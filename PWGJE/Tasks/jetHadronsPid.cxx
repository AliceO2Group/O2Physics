// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file jetHadronsPid.cxx
/// \brief Analysis of hadrons in jets
/// \author Aleksandra Mulewicz, WUT Warsaw, aleksandra.mulewicz@cern.ch
/// \author Leonard Lorenc, WUT Warsaw, leonard.lorenc@cern.ch
/// \author Hubert Zalewski, WUT Warsaw, hubert.zalewski@cern.ch
/// \author Małgorzata Janik malgorzata.janik@cern.ch
/// \author Daniela Ruggiano daniela.ruggiano@cern.ch

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
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TVector2.h>
#include <TVector3.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/PseudoJet.hh>

#include <cmath>
#include <set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::constants::math;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;

using HadronTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA,
                               aod::pidTPCFullPi, aod::pidTOFFullPi,
                               aod::pidTPCFullKa, aod::pidTOFFullKa,
                               aod::pidTPCFullPr, aod::pidTOFFullPr>;

using HadronTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA,
                                 aod::pidTPCFullPi, aod::pidTOFFullPi,
                                 aod::pidTPCFullKa, aod::pidTOFFullKa,
                                 aod::pidTPCFullPr, aod::pidTOFFullPr,
                                 aod::McTrackLabels>;

struct jetHadronsPid {

  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<bool> isppRefAnalysis{"isppRefAnalysis", false, "Is ppRef analysis"};
  Configurable<double> cfgAreaFrac{"cfgAreaFrac", 0.6, "fraction of jet area"};
  Configurable<double> cfgEtaJetMax{"cfgEtaJetMax", 0.5, "max jet eta"};
  Configurable<double> cfgMinPtTrack{"cfgMinPtTrack", 0.1, "minimum pt of tracks for jet reconstruction"};

  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "Reject events near the ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "Reject events near the TF border"};
  Configurable<bool> requireVtxITSTPC{"requireVtxITSTPC", true, "Require at least one ITS-TPC matched track"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "Reject events with same-bunch pileup collisions"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "Require consistent FT0 vs PV z-vertex"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "Require at least one vertex track matched to TOF"};

  Configurable<double> minJetPt{"minJetPt", 10.0, "Minimum pt of the jet after bkg subtraction"};
  Configurable<double> maxJetPt{"maxJetPt", 1e+06, "Maximum pt of the jet after bkg subtraction"};
  Configurable<double> rJet{"rJet", 0.4, "Jet resolution parameter R"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<bool> applyAreaCut{"applyAreaCut", true, "apply area cut"};
  Configurable<double> maxNormalizedJetArea{"maxNormalizedJetArea", 1.0, "area cut"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from the edge"};

  Configurable<bool> requirePvContributor{"requirePvContributor", false, "require that the track is a PV contributor"};
  Configurable<int> minItsNclusters{"minItsNclusters", 5, "minimum number of ITS clusters"};
  Configurable<int> minTpcNcrossedRows{"minTpcNcrossedRows", 70, "minimum number of TPC crossed pad rows"};
  Configurable<double> minChiSquareTpc{"minChiSquareTpc", 0.0, "minimum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0, "maximum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0, "maximum ITS chi^2/Ncls"};
  Configurable<double> minPt{"minPt", 0.3, "minimum pt of the tracks"};
  Configurable<double> maxPt{"maxPt", 4.0, "maximum pt of the tracks for PID analysis"};
  Configurable<double> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<double> maxEta{"maxEta", +0.8, "maximum eta"};
  Configurable<double> maxDcaxy{"maxDcaxy", 0.2, "Maximum DCAxy"};
  Configurable<double> maxDcaz{"maxDcaz", 0.1, "Maximum DCAz"};

  Configurable<bool> setMCDefaultItsParams{"setMCDefaultItsParams", true, "set MC default parameters"};

  JetBkgSubUtils backgroundSub;
  o2::aod::ITSResponse itsResponse;

  void init(InitContext const&)
  {
    if (setMCDefaultItsParams) {
      itsResponse.setMCDefaultParameters();
    }

    registryData.add("n_events", "Event counter", HistType::kTH1F, {{1, 0.5, 1.5, "N_{events}"}});
    registryData.add("n_events_raw", "All events", HistType::kTH1F, {{1, 0.5, 1.5, ""}});

    registryData.add("jet_pt_raw", "Jet pT raw", HistType::kTH1F, {{200, 0.0, 200.0, "#it{p}_{T}^{raw} (GeV/#it{c})"}});
    registryData.add("jet_pt_subtracted", "Jet pT subtracted", HistType::kTH1F, {{200, 0.0, 200.0, "#it{p}_{T}^{sub} (GeV/#it{c})"}});
    registryData.add("jet_pt_raw_vs_sub", "Raw vs sub jet pT", HistType::kTH2F, {{200, 0, 200}, {200, 0, 200}});
    registryData.add("jet_eta", "Jet eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta_{jet}"}});
    registryData.add("jet_phi", "Jet phi", HistType::kTH1F, {{100, 0.0, TwoPI, "#phi_{jet}"}});
    registryData.add("jet_area", "Jet area", HistType::kTH1F, {{100, 0.0, 1.5, "Area"}});
    registryData.add("jet_n_constituents", "Jet multiplicity", HistType::kTH1I, {{100, 0, 100, "N_{constituents}"}});

    registryData.add("pion_pure_tpc", "TPC Pion PID", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("pion_pure_tof", "TOF Pion PID", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("pion_pure_pt", "Pion pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("pion_pure_eta", "Pion Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("pion_pure_dcaz", "Pion DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("pion_pure_dcaxy", "Pion DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});

    registryData.add("kaon_pure_tpc", "TPC Kaon PID", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("kaon_pure_tof", "TOF Kaon PID", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("kaon_pure_pt", "Kaon pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("kaon_pure_eta", "Kaon Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("kaon_pure_dcaz", "Kaon DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("kaon_pure_dcaxy", "Kaon DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});

    registryData.add("proton_pure_tpc", "TPC Proton PID", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("proton_pure_tof", "TOF Proton PID", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("proton_pure_pt", "Proton pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("proton_pure_eta", "Proton Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("proton_pure_dcaz", "Proton DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("proton_pure_dcaxy", "Proton DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});

    registryData.add("z_vtx", "Z-Vertex Distribution", HistType::kTH1F, {{200, -20.0, 20.0, "Z-Vertex (cm)"}});
    registryData.add("tracks_in_jets", "Number of Tracks Inside Jets", HistType::kTH1F, {{100, 0, 100, "N_{tracks}"}});
    registryData.add("tracks_outside_jets", "Number of Tracks Outside Jets", HistType::kTH1F, {{500, 0, 500, "N_{tracks}"}});

    registryData.add("pion_jet_tpc", "TPC Pion PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("pion_jet_tof", "TOF Pion PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("kaon_jet_tpc", "TPC Kaon PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("kaon_jet_tof", "TOF Kaon PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("proton_jet_tpc", "TPC Proton PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("proton_jet_tof", "TOF Proton PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});

    registryData.add("pion_jet_pt", "Pion pT in Jets", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("pion_jet_eta", "Pion Eta in Jets", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("pion_jet_dcaxy", "Pion DCAxy in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("pion_jet_dcaz", "Pion DCAz in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("kaon_jet_pt", "Kaon pT in Jets", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("kaon_jet_eta", "Kaon Eta in Jets", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("kaon_jet_dcaxy", "Kaon DCAxy in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("kaon_jet_dcaz", "Kaon DCAz in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("proton_jet_pt", "Proton pT in Jets", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("proton_jet_eta", "Proton Eta in Jets", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("proton_jet_dcaxy", "Proton DCAxy in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("proton_jet_dcaz", "Proton DCAz in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("pion_ue_tpc", "TPC Pion PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("pion_ue_tof", "TOF Pion PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("pion_ue_pt", "Pion pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("pion_ue_eta", "Pion Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("pion_ue_dcaxy", "Pion DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("pion_ue_dcaz", "Pion DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("kaon_ue_tpc", "TPC Kaon PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("kaon_ue_tof", "TOF Kaon PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("kaon_ue_pt", "Kaon pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("kaon_ue_eta", "Kaon Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("kaon_ue_dcaxy", "Kaon DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("kaon_ue_dcaz", "Kaon DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("proton_ue_tpc", "TPC Proton PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("proton_ue_tof", "TOF Proton PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("proton_ue_pt", "Proton pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("proton_ue_eta", "Proton Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("proton_ue_dcaxy", "Proton DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("proton_ue_dcaz", "Proton DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("mc_gen_pion_pt", "Generated Pion pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc_rec_pion_pt", "Reconstructed Pion pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc_gen_kaon_pt", "Generated Kaon pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc_rec_kaon_pt", "Reconstructed Kaon pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc_gen_proton_pt", "Generated Proton pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc_rec_proton_pt", "Reconstructed Proton pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});

    registryData.add("rec_pion_all", "All identified pions", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T}"}});
    registryData.add("rec_kaon_all", "All identified kaons", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T}"}});
    registryData.add("rec_proton_all", "All identified protons", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T}"}});
    registryData.add("contamination_matrix_pion", "Pion Contamination Matrix;True absolute PDG;Rec #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{3000, -0.5, 2999.5}, {120, 0.0, 4.0}});
    registryData.add("contamination_matrix_kaon", "Kaon Contamination Matrix;True absolute PDG;Rec #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{3000, -0.5, 2999.5}, {120, 0.0, 4.0}});
    registryData.add("contamination_matrix_proton", "Proton Contamination Matrix;True absolute PDG;Rec #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{3000, -0.5, 2999.5}, {120, 0.0, 4.0}});

    registryData.add("mc_sec_pion_pt", "Reconstructed Pion Secondaries", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc_sec_kaon_pt", "Reconstructed Kaon Secondaries", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc_sec_proton_pt", "Reconstructed Proton Secondaries", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
  }

  void getPerpendicularDirections(const TVector3& p, TVector3& u1, TVector3& u2)
  {
    double px = p.X(), py = p.Y(), pz = p.Z();
    double px2 = px * px, py2 = py * py, pz2 = pz * pz;
    double pz4 = pz2 * pz2;

    if (px == 0 && py == 0) {
      u1.SetXYZ(0, 0, 0);
      u2.SetXYZ(0, 0, 0);
      return;
    }
    if (px == 0 && py != 0) {
      double ux = std::sqrt(py2 - pz4 / py2);
      double uy = -pz2 / py;
      u1.SetXYZ(ux, uy, pz);
      u2.SetXYZ(-ux, uy, pz);
      return;
    }
    if (py == 0 && px != 0) {
      double ux = -pz2 / px;
      double uy = std::sqrt(px2 - pz4 / px2);
      u1.SetXYZ(ux, uy, pz);
      u2.SetXYZ(ux, -uy, pz);
      return;
    }

    double a = px2 + py2;
    double b = 2.0 * px * pz2;
    double c = pz4 - py2 * py2 - px2 * py2;
    double delta = b * b - 4.0 * a * c;

    if (delta < 0 || a == 0) {
      u1.SetXYZ(0, 0, 0);
      u2.SetXYZ(0, 0, 0);
      return;
    }
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
    double diff = std::abs(phi1 - phi2);

    if (diff <= PI)
      deltaPhi = diff;
    if (diff > PI)
      deltaPhi = TwoPI - diff;
    return deltaPhi;
  }

  template <typename TrackIts>
  bool hasITSHit(const TrackIts& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

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

    if (!track.hasITS() || !track.hasTPC())
      return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    if (track.tpcNClsCrossedRows() < MinTpcCr)
      return false;
    if (track.tpcChi2NCl() > MaxChi2Tpc)
      return false;
    if (track.itsChi2NCl() > MaxChi2Its)
      return false;
    if (std::abs(track.eta()) > maxEta)
      return false;
    if (track.pt() < cfgMinPtTrack)
      return false;
    if (std::abs(track.dcaXY()) > (DcaxyMaxTrackPar0 + DcaxyMaxTrackPar1 / std::pow(track.pt(), DcaxyMaxTrackPar2)))
      return false;
    if (std::abs(track.dcaZ()) > DcazMaxTrack)
      return false;
    return true;
  }

  template <typename PionTrack>
  bool passedTrackSelection(const PionTrack& track)
  {
    if (requirePvContributor && !(track.isPVContributor()))
      return false;
    if (!track.hasITS() || !track.hasTPC())
      return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    if (track.itsNCls() < minItsNclusters)
      return false;
    if (track.tpcNClsCrossedRows() < minTpcNcrossedRows)
      return false;
    if (track.tpcChi2NCl() < minChiSquareTpc || track.tpcChi2NCl() > maxChiSquareTpc)
      return false;
    if (track.itsChi2NCl() > maxChiSquareIts)
      return false;
    if (track.eta() < minEta || track.eta() > maxEta)
      return false;
    if (track.pt() < minPt || track.pt() > maxPt)
      return false;
    return true;
  }

  void processForJets(SelectedCollisions::iterator const& collision, HadronTracks const& tracks)
  {
    registryData.fill(HIST("n_events_raw"), 1);

    if (!collision.sel8() || std::abs(collision.posZ()) > zVtx)
      return;
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return;
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      return;
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return;
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      return;

    registryData.fill(HIST("n_events"), 1);
    registryData.fill(HIST("z_vtx"), collision.posZ());

    for (auto const& track : tracks) {
      if (!passedTrackSelection(track))
        continue;
      if (std::abs(track.dcaXY()) > maxDcaxy || std::abs(track.dcaZ()) > maxDcaz)
        continue;

      double pt = track.pt();
      double eta = track.eta();
      double dcaxy = track.dcaXY();
      double dcaz = track.dcaZ();

      bool passTpcPi = (std::abs(track.tpcNSigmaPi()) <= 3.0);
      bool passTpcKa = (std::abs(track.tpcNSigmaKa()) <= 3.0);
      bool passTpcPr = (std::abs(track.tpcNSigmaPr()) <= 3.0);

      bool passTofPi = track.hasTOF() && (std::abs(track.tofNSigmaPi()) <= 3.0);
      bool passTofKa = track.hasTOF() && (std::abs(track.tofNSigmaKa()) <= 3.0);
      bool passTofPr = track.hasTOF() && (std::abs(track.tofNSigmaPr()) <= 3.0);

      const double ptThreshold = 0.8;

      if (passTpcPi) {
        registryData.fill(HIST("pion_pure_tpc"), pt, track.tpcNSigmaPi());

        if (passTofPi) {
          registryData.fill(HIST("pion_pure_tof"), pt, track.tofNSigmaPi());
        }

        if (pt < ptThreshold || passTofPi) {
          registryData.fill(HIST("pion_pure_pt"), pt);
          registryData.fill(HIST("pion_pure_eta"), eta);
          registryData.fill(HIST("pion_pure_dcaxy"), dcaxy);
          registryData.fill(HIST("pion_pure_dcaz"), dcaz);
        }
      }

      if (passTpcKa) {
        registryData.fill(HIST("kaon_pure_tpc"), pt, track.tpcNSigmaKa());

        if (passTofKa) {
          registryData.fill(HIST("kaon_pure_tof"), pt, track.tofNSigmaKa());
        }

        if (pt < ptThreshold || passTofKa) {
          registryData.fill(HIST("kaon_pure_pt"), pt);
          registryData.fill(HIST("kaon_pure_eta"), eta);
          registryData.fill(HIST("kaon_pure_dcaxy"), dcaxy);
          registryData.fill(HIST("kaon_pure_dcaz"), dcaz);
        }
      }

      if (passTpcPr) {
        registryData.fill(HIST("proton_pure_tpc"), pt, track.tpcNSigmaPr());

        if (passTofPr) {
          registryData.fill(HIST("proton_pure_tof"), pt, track.tofNSigmaPr());
        }

        if (pt < ptThreshold || passTofPr) {
          registryData.fill(HIST("proton_pure_pt"), pt);
          registryData.fill(HIST("proton_pure_eta"), eta);
          registryData.fill(HIST("proton_pure_dcaxy"), dcaxy);
          registryData.fill(HIST("proton_pure_dcaz"), dcaz);
        }
      }
    }

    int id(-1);
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {
      id++;
      if (!passedTrackSelectionForJetReconstruction(track))
        continue;

      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
      fourMomentum.set_user_index(id);
      fjParticles.emplace_back(fourMomentum);
    }

    if (fjParticles.empty())
      return;

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

    if (jets.empty())
      return;

    auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

    std::set<int> tracksInJetsSet;

    for (const auto& jet : jets) {
      registryData.fill(HIST("jet_pt_raw"), jet.pt());
      registryData.fill(HIST("jet_eta"), jet.eta());
      registryData.fill(HIST("jet_phi"), jet.phi());
      registryData.fill(HIST("jet_area"), jet.area());
      registryData.fill(HIST("jet_n_constituents"), jet.constituents().size());

      if (!isppRefAnalysis && ((std::abs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge)))
        continue;
      if (isppRefAnalysis && std::abs(jet.eta()) > cfgEtaJetMax)
        continue;

      auto jetForSub = jet;
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
      registryData.fill(HIST("jet_pt_subtracted"), jetMinusBkg.pt());
      registryData.fill(HIST("jet_pt_raw_vs_sub"), jet.pt(), jetMinusBkg.pt());

      if (isppRefAnalysis && (jet.pt() < minJetPt || jet.pt() > maxJetPt))
        continue;
      if (!isppRefAnalysis && (jetMinusBkg.pt() < minJetPt || jetMinusBkg.pt() > maxJetPt))
        continue;

      double normalizedJetArea = jet.area() / (PI * rJet * rJet);
      if (applyAreaCut && (!isppRefAnalysis) && normalizedJetArea > maxNormalizedJetArea)
        continue;
      if (isppRefAnalysis && (jet.area() < cfgAreaFrac * PI * rJet * rJet))
        continue;

      double coneRadius = std::sqrt(jet.area() / PI);
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
      getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);

      if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0)
        continue;

      std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();

      for (const auto& particle : jetConstituents) {

        int trackIdx = particle.user_index();
        auto const& track = tracks.iteratorAt(trackIdx);

        tracksInJetsSet.insert(trackIdx);

        if (!passedTrackSelection(track))
          continue;
        if (std::abs(track.dcaXY()) > maxDcaxy || std::abs(track.dcaZ()) > maxDcaz)
          continue;

        double pt = track.pt();
        double eta = track.eta();
        double dcaxy = track.dcaXY();
        double dcaz = track.dcaZ();

        bool passTpcPi = (std::abs(track.tpcNSigmaPi()) <= 3.0);
        bool passTpcKa = (std::abs(track.tpcNSigmaKa()) <= 3.0);
        bool passTpcPr = (std::abs(track.tpcNSigmaPr()) <= 3.0);

        bool passTofPi = track.hasTOF() && (std::abs(track.tofNSigmaPi()) <= 3.0);
        bool passTofKa = track.hasTOF() && (std::abs(track.tofNSigmaKa()) <= 3.0);
        bool passTofPr = track.hasTOF() && (std::abs(track.tofNSigmaPr()) <= 3.0);

        const double ptThreshold = 0.8;

        if (passTpcPi) {
          registryData.fill(HIST("pion_jet_tpc"), pt, track.tpcNSigmaPi());

          if (passTofPi) {
            registryData.fill(HIST("pion_jet_tof"), pt, track.tofNSigmaPi());
          }

          if (pt < ptThreshold || passTofPi) {
            registryData.fill(HIST("pion_jet_pt"), pt);
            registryData.fill(HIST("pion_jet_eta"), eta);
            registryData.fill(HIST("pion_jet_dcaxy"), dcaxy);
            registryData.fill(HIST("pion_jet_dcaz"), dcaz);
          }
        }

        if (passTpcKa) {
          registryData.fill(HIST("kaon_jet_tpc"), pt, track.tpcNSigmaKa());

          if (passTofKa) {
            registryData.fill(HIST("kaon_jet_tof"), pt, track.tofNSigmaKa());
          }

          if (pt < ptThreshold || passTofKa) {
            registryData.fill(HIST("kaon_jet_pt"), pt);
            registryData.fill(HIST("kaon_jet_eta"), eta);
            registryData.fill(HIST("kaon_jet_dcaxy"), dcaxy);
            registryData.fill(HIST("kaon_jet_dcaz"), dcaz);
          }
        }

        if (passTpcPr) {
          registryData.fill(HIST("proton_jet_tpc"), pt, track.tpcNSigmaPr());

          if (passTofPr) {
            registryData.fill(HIST("proton_jet_tof"), pt, track.tofNSigmaPr());
          }

          if (pt < ptThreshold || passTofPr) {
            registryData.fill(HIST("proton_jet_pt"), pt);
            registryData.fill(HIST("proton_jet_eta"), eta);
            registryData.fill(HIST("proton_jet_dcaxy"), dcaxy);
            registryData.fill(HIST("proton_jet_dcaz"), dcaz);
          }
        }
      }

      for (auto const& track : tracks) {

        if (!passedTrackSelection(track))
          continue;
        if (std::abs(track.dcaXY()) > maxDcaxy || std::abs(track.dcaZ()) > maxDcaz)
          continue;

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

        if (deltaRUe1 > maxConeRadius && deltaRUe2 > maxConeRadius)
          continue;

        double pt = track.pt();
        double eta = track.eta();
        double dcaxy = track.dcaXY();
        double dcaz = track.dcaZ();

        bool passTpcPi = (std::abs(track.tpcNSigmaPi()) <= 3.0);
        bool passTpcKa = (std::abs(track.tpcNSigmaKa()) <= 3.0);
        bool passTpcPr = (std::abs(track.tpcNSigmaPr()) <= 3.0);

        bool passTofPi = track.hasTOF() && (std::abs(track.tofNSigmaPi()) <= 3.0);
        bool passTofKa = track.hasTOF() && (std::abs(track.tofNSigmaKa()) <= 3.0);
        bool passTofPr = track.hasTOF() && (std::abs(track.tofNSigmaPr()) <= 3.0);

        const double ptThreshold = 0.8;

        if (passTpcPi) {
          registryData.fill(HIST("pion_ue_tpc"), pt, track.tpcNSigmaPi());
          if (passTofPi) {
            registryData.fill(HIST("pion_ue_tof"), pt, track.tofNSigmaPi());
          }
          if (pt < ptThreshold || passTofPi) {
            registryData.fill(HIST("pion_ue_pt"), pt);
            registryData.fill(HIST("pion_ue_eta"), eta);
            registryData.fill(HIST("pion_ue_dcaxy"), dcaxy);
            registryData.fill(HIST("pion_ue_dcaz"), dcaz);
          }
        }

        if (passTpcKa) {
          registryData.fill(HIST("kaon_ue_tpc"), pt, track.tpcNSigmaKa());
          if (passTofKa) {
            registryData.fill(HIST("kaon_ue_tof"), pt, track.tofNSigmaKa());
          }
          if (pt < ptThreshold || passTofKa) {
            registryData.fill(HIST("kaon_ue_pt"), pt);
            registryData.fill(HIST("kaon_ue_eta"), eta);
            registryData.fill(HIST("kaon_ue_dcaxy"), dcaxy);
            registryData.fill(HIST("kaon_ue_dcaz"), dcaz);
          }
        }

        if (passTpcPr) {
          registryData.fill(HIST("proton_ue_tpc"), pt, track.tpcNSigmaPr());
          if (passTofPr) {
            registryData.fill(HIST("proton_ue_tof"), pt, track.tofNSigmaPr());
          }
          if (pt < ptThreshold || passTofPr) {
            registryData.fill(HIST("proton_ue_pt"), pt);
            registryData.fill(HIST("proton_ue_eta"), eta);
            registryData.fill(HIST("proton_ue_dcaxy"), dcaxy);
            registryData.fill(HIST("proton_ue_dcaz"), dcaz);
          }
        }
      }
    }

    registryData.fill(HIST("tracks_in_jets"), tracksInJetsSet.size());

    int nTracksOut = 0;
    int trackIteratorIdx = 0;
    for (auto const& track : tracks) {
      if (passedTrackSelection(track) && tracksInJetsSet.find(trackIteratorIdx) == tracksInJetsSet.end()) {
        nTracksOut++;
      }
      trackIteratorIdx++;
    }
    registryData.fill(HIST("tracks_outside_jets"), nTracksOut);
  }

  PROCESS_SWITCH(jetHadronsPid, processForJets, "Pid Analysis in and outside jets", true);

  void processMC(SelectedCollisions::iterator const& collision, aod::McParticles const& mcParticles, HadronTracksMC const& tracks)
  {
    if (!collision.sel8() || std::abs(collision.posZ()) > zVtx)
      return;
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return;
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      return;
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return;
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      return;

    for (auto const& mcpart : mcParticles) {
      if (!mcpart.isPhysicalPrimary())
        continue;
      if (std::abs(mcpart.eta()) > 0.8)
        continue;

      int pdg = std::abs(mcpart.pdgCode());
      double pt = mcpart.pt();

      if (pdg == 211) {
        registryData.fill(HIST("mc_gen_pion_pt"), pt);
      } else if (pdg == 321) {
        registryData.fill(HIST("mc_gen_kaon_pt"), pt);
      } else if (pdg == 2212) {
        registryData.fill(HIST("mc_gen_proton_pt"), pt);
      }
    }

    const double ptThreshold = 0.8;

    for (auto const& track : tracks) {
      if (!passedTrackSelection(track))
        continue;
      if (std::abs(track.dcaXY()) > maxDcaxy || std::abs(track.dcaZ()) > maxDcaz)
        continue;

      double pt = track.pt();

      bool passTpcPi = (std::abs(track.tpcNSigmaPi()) <= 3.0);
      bool passTofPi = track.hasTOF() && (std::abs(track.tofNSigmaPi()) <= 3.0);

      bool passTpcKa = (std::abs(track.tpcNSigmaKa()) <= 3.0);
      bool passTofKa = track.hasTOF() && (std::abs(track.tofNSigmaKa()) <= 3.0);

      bool passTpcPr = (std::abs(track.tpcNSigmaPr()) <= 3.0);
      bool passTofPr = track.hasTOF() && (std::abs(track.tofNSigmaPr()) <= 3.0);

      if (!track.has_mcParticle())
        continue;
      auto const& trueParticle = track.mcParticle();

      int pdg = std::abs(trueParticle.pdgCode());
      bool isPrimary = trueParticle.isPhysicalPrimary();

      if (passTpcPi && (pt < ptThreshold || passTofPi)) {
        registryData.fill(HIST("rec_pion_all"), pt);
        registryData.fill(HIST("contamination_matrix_pion"), pdg, pt);
        if (isPrimary) {
          if (pdg == 211)
            registryData.fill(HIST("mc_rec_pion_pt"), pt);
        } else {
          registryData.fill(HIST("mc_sec_pion_pt"), pt);
        }
      }

      if (passTpcKa && (pt < ptThreshold || passTofKa)) {
        registryData.fill(HIST("rec_kaon_all"), pt);
        registryData.fill(HIST("contamination_matrix_kaon"), pdg, pt);
        if (isPrimary) {
          if (pdg == 321)
            registryData.fill(HIST("mc_rec_kaon_pt"), pt);
        } else {
          registryData.fill(HIST("mc_sec_kaon_pt"), pt);
        }
      }

      if (passTpcPr && (pt < ptThreshold || passTofPr)) {
        registryData.fill(HIST("rec_proton_all"), pt);
        registryData.fill(HIST("contamination_matrix_proton"), pdg, pt);
        if (isPrimary) {
          if (pdg == 2212)
            registryData.fill(HIST("mc_rec_proton_pt"), pt);
        } else {
          registryData.fill(HIST("mc_sec_proton_pt"), pt);
        }
      }
    }
  }
  PROCESS_SWITCH(jetHadronsPid, processMC, "Run on Monte Carlo", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<jetHadronsPid>(cfgc)};
}
