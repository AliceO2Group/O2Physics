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
/// \author Alberto Caliva (alberto.caliva@cern.ch), Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since May 31, 2024

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>
#include <TVector3.h>
#include <cmath>
#include <vector>

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/DCA.h"
#define mXi 1.32171
#define mOmega 1.67245

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

using FullTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

struct tracked_cascade_properties {

  // QC Histograms
  HistogramRegistry registryQC{
    "registryQC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Analysis Histograms: Data
  HistogramRegistry registryData{
    "registryData",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Global Parameters
  Configurable<float> zVtx{"zVtx", 10.0f, "z vertex cut"};

  // Cascade Parameters
  Configurable<float> minimumCascRadius{"minimumCascRadius", 5.0f, "Minimum Cascade Radius"};
  Configurable<float> maximumCascRadius{"maximumCascRadius", 18.0f, "Maximum Cascade Radius"};

  // Mass Cuts
  Configurable<float> mMin_xi{"mMin_xi", 1.315f, "mMin Xi"};
  Configurable<float> mMax_xi{"mMax_xi", 1.328f, "mMax Xi"};
  Configurable<float> mMin_omega{"mMin_omega", 1.665f, "mMin Omega"};
  Configurable<float> mMax_omega{"mMax_omega", 1.680f, "mMax Omega"};

  void init(InitContext const&)
  {
    registryQC.add("matchingChi2", "matching Chi2", HistType::kTH1F, {{200, 0, 1000, "#chi^{2}_{matching}"}});
    registryQC.add("topologyChi2", "topology Chi2", HistType::kTH1F, {{500, 0, 0.5, "#chi^{2}_{topology}"}});
    registryQC.add("nITScls_vs_p_xi", "nITS Xi", HistType::kTH2F, {{100, 0, 10, "#it{p} (GeV/#it{c})"}, {8, 0, 8, "n_{ITS}^{cls}"}});
    registryQC.add("nITScls_vs_p_omega", "nITS Omega", HistType::kTH2F, {{100, 0, 10, "#it{p} (GeV/#it{c})"}, {8, 0, 8, "n_{ITS}^{cls}"}});
    registryQC.add("decayXY", "decayXY", HistType::kTH2F, {{500, -50, 50, "x"}, {500, -50, 50, "y"}});
    registryQC.add("deltaClsSize", "deltaClsSize", HistType::kTH1F, {{40, -20, 20, "#DeltaClsSize"}});
    registryQC.add("deltaP", "deltaP", HistType::kTH1F, {{1000, -1, 1, "#Deltap"}});
    registryQC.add("deltaEta", "deltaEta", HistType::kTH1F, {{200, -0.5, 0.5, "#Delta#eta"}});
    registryQC.add("deltaNclsITS", "deltaNclsITS", HistType::kTH1F, {{20, -10, 10, "#DeltaN"}});
    registryQC.add("deltaNclsITS_track", "deltaNclsITS_track", HistType::kTH1F, {{20, -10, 10, "#DeltaN"}});
    registryQC.add("deltaNclsITS_itstrack", "deltaNclsITS_itstrack", HistType::kTH1F, {{20, -10, 10, "#DeltaN"}});

    registryData.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{5, 0, 5, "Event Cuts"}});
    registryData.add("xi_pos_avgclustersize", "xi_pos_avgclustersize", HistType::kTH3F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("xi_neg_avgclustersize", "xi_neg_avgclustersize", HistType::kTH3F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("omega_pos_avgclustersize", "omega_pos_avgclustersize", HistType::kTH3F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("omega_neg_avgclustersize", "omega_neg_avgclustersize", HistType::kTH3F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {16, -0.8, 0.8, "#eta"}});

    registryData.add("xi_pos_avgclustersize_cosL", "xi_pos_avgclustersize_cosL", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("xi_neg_avgclustersize_cosL", "xi_neg_avgclustersize_cosL", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_pos_avgclustersize_cosL", "omega_pos_avgclustersize_cosL", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_neg_avgclustersize_cosL", "omega_neg_avgclustersize_cosL", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});

    registryData.add("xi_pos_avgclustersize_cosL_vs_betagamma", "xi_pos_avgclustersize_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("xi_neg_avgclustersize_cosL_vs_betagamma", "xi_neg_avgclustersize_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_pos_avgclustersize_cosL_vs_betagamma", "omega_pos_avgclustersize_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_neg_avgclustersize_cosL_vs_betagamma", "omega_neg_avgclustersize_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});

    registryData.add("xi_mass_pos", "xi_mass_pos", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("xi_mass_neg", "xi_mass_neg", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("omega_mass_pos", "omega_mass_pos", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
    registryData.add("omega_mass_neg", "omega_mass_neg", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
  }

  double track_inclination(double eta)
  {
    double lambda(0);
    double theta = 2.0 * atan(exp(-eta));
    if (theta <= TMath::Pi() / 2.0)
      lambda = 0.5 * TMath::Pi() - theta;
    if (theta > TMath::Pi() / 2.0)
      lambda = theta - 0.5 * TMath::Pi();
    return lambda;
  }

  void processData(SelectedCollisions::iterator const& collision, aod::AssignedTrackedCascades const& trackedCascades,
                   aod::Cascades const&, FullTracks const&)
  {
    registryData.fill(HIST("number_of_events_data"), 0.5);
    if (!collision.sel8())
      return;

    registryData.fill(HIST("number_of_events_data"), 1.5);
    if (abs(collision.posZ()) > zVtx)
      return;

    registryData.fill(HIST("number_of_events_data"), 2.5);

    for (const auto& trackedCascade : trackedCascades) {

      const auto track = trackedCascade.track_as<FullTracks>();
      const auto trackITS = trackedCascade.itsTrack_as<FullTracks>();

      // Comparison between track and ITStrack
      registryQC.fill(HIST("deltaP"), track.p() - trackITS.p());
      registryQC.fill(HIST("deltaEta"), track.eta() - trackITS.eta());
      registryQC.fill(HIST("deltaNclsITS"), track.itsNCls() - trackITS.itsNCls());
      for (int i = 0; i < 7; i++) {
        registryQC.fill(HIST("deltaClsSize"), track.itsClsSizeInLayer(i) - trackITS.itsClsSizeInLayer(i));
      }

      const auto& casc = trackedCascade.cascade();
      const auto& btrack = casc.bachelor_as<FullTracks>();
      double dx = trackedCascade.decayX();
      double dy = trackedCascade.decayY();
      double r = sqrt(dx * dx + dy * dy);
      if (r < minimumCascRadius || r > maximumCascRadius)
        continue;

      registryQC.fill(HIST("matchingChi2"), trackedCascade.matchingChi2());
      registryQC.fill(HIST("topologyChi2"), trackedCascade.topologyChi2());
      registryQC.fill(HIST("decayXY"), dx, dy);

      // Calculate (Average) Cluster Size
      double averageClusterSize(0);
      int nCls(0);
      for (int i = 0; i < 7; i++) {
        int clusterSize = trackITS.itsClsSizeInLayer(i);
        averageClusterSize += static_cast<double>(clusterSize);
        if (clusterSize > 0)
          nCls++;
      }
      averageClusterSize = averageClusterSize / static_cast<double>(nCls);

      registryQC.fill(HIST("deltaNclsITS_track"), nCls - track.itsNCls());
      registryQC.fill(HIST("deltaNclsITS_itstrack"), nCls - trackITS.itsNCls());

      // Xi Mass
      if (btrack.sign() > 0) {
        registryData.fill(HIST("xi_mass_pos"), track.p(), trackedCascade.xiMass());
      }
      if (btrack.sign() < 0) {
        registryData.fill(HIST("xi_mass_neg"), track.p(), trackedCascade.xiMass());
      }

      // Track Inclination
      double lambda = track_inclination(track.eta());

      // Xi
      if (trackedCascade.xiMass() > mMin_xi && trackedCascade.xiMass() < mMax_xi) {
        registryQC.fill(HIST("nITScls_vs_p_xi"), track.p(), trackITS.itsNCls());
        if (btrack.sign() > 0) {
          registryData.fill(HIST("xi_pos_avgclustersize"), track.p(), averageClusterSize, track.eta());
          registryData.fill(HIST("xi_pos_avgclustersize_cosL"), track.p(), averageClusterSize * cos(lambda));
          registryData.fill(HIST("xi_pos_avgclustersize_cosL_vs_betagamma"), track.p() / mXi, averageClusterSize * cos(lambda));
        }
        if (btrack.sign() < 0) {
          registryData.fill(HIST("xi_neg_avgclustersize"), track.p(), averageClusterSize, track.eta());
          registryData.fill(HIST("xi_neg_avgclustersize_cosL"), track.p(), averageClusterSize * cos(lambda));
          registryData.fill(HIST("xi_neg_avgclustersize_cosL_vs_betagamma"), track.p() / mXi, averageClusterSize * cos(lambda));
        }
        continue;
      }

      // Omega Mass
      if (btrack.sign() > 0) {
        registryData.fill(HIST("omega_mass_pos"), track.p(), trackedCascade.omegaMass());
      }
      if (btrack.sign() < 0) {
        registryData.fill(HIST("omega_mass_neg"), track.p(), trackedCascade.omegaMass());
      }

      // Omega
      if (trackedCascade.omegaMass() > mMin_omega && trackedCascade.omegaMass() < mMax_omega) {
        registryQC.fill(HIST("nITScls_vs_p_omega"), track.p(), trackITS.itsNCls());
        if (btrack.sign() > 0) {
          registryData.fill(HIST("omega_pos_avgclustersize"), track.p(), averageClusterSize, track.eta());
          registryData.fill(HIST("omega_pos_avgclustersize_cosL"), track.p(), averageClusterSize * cos(lambda));
          registryData.fill(HIST("omega_pos_avgclustersize_cosL_vs_betagamma"), track.p() / mOmega, averageClusterSize * cos(lambda));
        }
        if (btrack.sign() < 0) {
          registryData.fill(HIST("omega_neg_avgclustersize"), track.p(), averageClusterSize, track.eta());
          registryData.fill(HIST("omega_neg_avgclustersize_cosL"), track.p(), averageClusterSize * cos(lambda));
          registryData.fill(HIST("omega_neg_avgclustersize_cosL_vs_betagamma"), track.p() / mOmega, averageClusterSize * cos(lambda));
        }
      }
    }
  }
  PROCESS_SWITCH(tracked_cascade_properties, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tracked_cascade_properties>(cfgc)};
}
