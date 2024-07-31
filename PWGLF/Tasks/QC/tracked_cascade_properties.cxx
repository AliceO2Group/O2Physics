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
  Configurable<float> mMin_xi{"mMin_xi", 1.31f, "mMin Xi"};
  Configurable<float> mMax_xi{"mMax_xi", 1.33f, "mMax Xi"};
  Configurable<float> mMin_omega{"mMin_omega", 1.66f, "mMin Omega"};
  Configurable<float> mMax_omega{"mMax_omega", 1.68f, "mMax Omega"};

  void init(InitContext const&)
  {
    registryQC.add("matchingChi2", "matching Chi2", HistType::kTH1F, {{200, 0, 400, "#chi^{2}_{matching}"}});
    registryQC.add("topologyChi2", "topology Chi2", HistType::kTH1F, {{200, 0, 20, "#chi^{2}_{topology}"}});
    registryQC.add("nITS_Xi", "nITS Xi", HistType::kTH2F, {{100, 0, 10, "#it{p} (GeV/#it{c})"}, {8, 0, 8, "n_{ITS}^{cls}"}});
    registryQC.add("nITS_Omega", "nITS Omega", HistType::kTH2F, {{100, 0, 10, "#it{p} (GeV/#it{c})"}, {8, 0, 8, "n_{ITS}^{cls}"}});
    registryQC.add("tgl_Distr", "tgl_Distr", HistType::kTH1F, {{200, 0, 200, "tan (#lambda)"}});
    registryQC.add("nITSclusters", "nITSclusters", HistType::kTH1F, {{7, 0, 7, "n_{cls}^{ITS}"}});
    registryQC.add("clusterSize", "clusterSize", HistType::kTH1F, {{200, 0, 200, "cluster size ITS"}});
    registryQC.add("decayXY", "decayXY", HistType::kTH2F, {{200, -50, 50, "x"}, {200, -50, 50, "y"}});
    registryQC.add("signBachelor", "signBachelor", HistType::kTH1F, {{4, -2, 2, "sign"}});
    registryQC.add("ITSclusterSizeTrkCasc", "ITSclusterSizeTrkCasc", HistType::kTH1F, {{200, 0, 20, "ITS cluster Size"}});
    registryQC.add("DeltaLambda", "DeltaLambda", HistType::kTH1F, {{200, -0.1, 0.1, "#lambda - #lambda_{1}"}});

    registryData.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{5, 0, 5, "Event Cuts"}});
    registryData.add("xi_pos_clustersize", "xi_pos_clustersize", HistType::kTH3F, {{100, 0.0, 100, "ITS cluster size"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("xi_neg_clustersize", "xi_neg_clustersize", HistType::kTH3F, {{100, 0.0, 100, "ITS cluster size"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("omega_pos_clustersize", "omega_pos_clustersize", HistType::kTH3F, {{100, 0.0, 100, "ITS cluster size"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("omega_neg_clustersize", "omega_neg_clustersize", HistType::kTH3F, {{100, 0.0, 100, "ITS cluster size"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});

    registryData.add("xi_pos_avgclustersize", "xi_pos_avgclustersize", HistType::kTH3F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("xi_neg_avgclustersize", "xi_neg_avgclustersize", HistType::kTH3F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("omega_pos_avgclustersize", "omega_pos_avgclustersize", HistType::kTH3F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("omega_neg_avgclustersize", "omega_neg_avgclustersize", HistType::kTH3F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});

    registryData.add("xi_pos_avgclustersize_cosL", "xi_pos_avgclustersize_cosL", HistType::kTH2F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}});
    registryData.add("xi_neg_avgclustersize_cosL", "xi_neg_avgclustersize_cosL", HistType::kTH2F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}});
    registryData.add("omega_pos_avgclustersize_cosL", "omega_pos_avgclustersize_cosL", HistType::kTH2F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}});
    registryData.add("omega_neg_avgclustersize_cosL", "omega_neg_avgclustersize_cosL", HistType::kTH2F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}});

    registryData.add("xi_pos_avgclustersize_sinL", "xi_pos_avgclustersize_sinL", HistType::kTH2F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT sin(#lambda)"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}});
    registryData.add("xi_neg_avgclustersize_sinL", "xi_neg_avgclustersize_sinL", HistType::kTH2F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT sin(#lambda)"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}});
    registryData.add("omega_pos_avgclustersize_sinL", "omega_pos_avgclustersize_sinL", HistType::kTH2F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT sin(#lambda)"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}});
    registryData.add("omega_neg_avgclustersize_sinL", "omega_neg_avgclustersize_sinL", HistType::kTH2F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT sin(#lambda)"}, {100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}});

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

      const auto track = trackedCascade.itsTrack_as<FullTracks>();
      const auto& casc = trackedCascade.cascade();
      const auto& btrack = casc.bachelor_as<FullTracks>();
      double dx = trackedCascade.decayX();
      double dy = trackedCascade.decayY();
      double r = sqrt(dx * dx + dy * dy);
      if (r < minimumCascRadius || r > maximumCascRadius)
        continue;

      registryQC.fill(HIST("matchingChi2"), trackedCascade.matchingChi2());
      registryQC.fill(HIST("topologyChi2"), trackedCascade.topologyChi2());
      registryQC.fill(HIST("nITSclusters"), track.itsNCls());
      registryQC.fill(HIST("ITSclusterSizeTrkCasc"), trackedCascade.itsClsSize());
      registryQC.fill(HIST("decayXY"), dx, dy);
      registryQC.fill(HIST("signBachelor"), btrack.sign());

      // Calculate (Average) Cluster Size
      int clusterSize[7];
      double averageClusterSize(0);
      for (int i = 0; i < 7; i++) {
        clusterSize[i] = track.itsClsSizeInLayer(i);
        registryQC.fill(HIST("clusterSize"), clusterSize[i]);
        averageClusterSize += static_cast<double>(clusterSize[i]);
      }
      averageClusterSize = averageClusterSize / static_cast<double>(track.itsNCls());

      // Track Inclination
      registryQC.fill(HIST("tgl_Distr"), track.tgl());
      double lambda = track_inclination(track.eta());
      double lambda1 = atan(track.tgl());
      double cosL = cos(lambda);
      double sinL = sin(lambda);
      registryQC.fill(HIST("DeltaLambda"), lambda - lambda1);

      // Xi
      if (trackedCascade.xiMass() > mMin_xi && trackedCascade.xiMass() < mMax_xi) {
        registryQC.fill(HIST("nITS_Xi"), track.p(), track.itsNCls());
        if (btrack.sign() > 0) {
          for (int i = 0; i < 7; i++) {
            registryData.fill(HIST("xi_pos_clustersize"), clusterSize[i], track.p(), track.eta());
          }
          registryData.fill(HIST("xi_pos_avgclustersize"), averageClusterSize, track.p(), track.eta());
          registryData.fill(HIST("xi_pos_avgclustersize_cosL"), averageClusterSize * cosL, track.p());
          registryData.fill(HIST("xi_pos_avgclustersize_sinL"), averageClusterSize * sinL, track.p());
          registryData.fill(HIST("xi_mass_pos"), track.p(), trackedCascade.xiMass());
        }
        if (btrack.sign() < 0) {
          for (int i = 0; i < 7; i++) {
            registryData.fill(HIST("xi_neg_clustersize"), clusterSize[i], track.p(), track.eta());
          }
          registryData.fill(HIST("xi_neg_avgclustersize"), averageClusterSize, track.p(), track.eta());
          registryData.fill(HIST("xi_neg_avgclustersize_cosL"), averageClusterSize * cosL, track.p());
          registryData.fill(HIST("xi_neg_avgclustersize_sinL"), averageClusterSize * sinL, track.p());
          registryData.fill(HIST("xi_mass_neg"), track.p(), trackedCascade.xiMass());
        }
      }

      // Omega
      if (trackedCascade.omegaMass() > mMin_omega && trackedCascade.omegaMass() < mMax_omega) {
        registryQC.fill(HIST("nITS_Omega"), track.p(), track.itsNCls());
        if (btrack.sign() > 0) {
          for (int i = 0; i < 7; i++) {
            registryData.fill(HIST("omega_pos_clustersize"), clusterSize[i], track.p(), track.eta());
          }
          registryData.fill(HIST("omega_pos_avgclustersize"), averageClusterSize, track.p(), track.eta());
          registryData.fill(HIST("omega_pos_avgclustersize_cosL"), averageClusterSize * cosL, track.p());
          registryData.fill(HIST("omega_pos_avgclustersize_sinL"), averageClusterSize * sinL, track.p());
          registryData.fill(HIST("omega_mass_pos"), track.p(), trackedCascade.omegaMass());
        }
        if (btrack.sign() < 0) {
          for (int i = 0; i < 7; i++) {
            registryData.fill(HIST("omega_neg_clustersize"), clusterSize[i], track.p(), track.eta());
          }
          registryData.fill(HIST("omega_neg_avgclustersize"), averageClusterSize, track.p(), track.eta());
          registryData.fill(HIST("omega_neg_avgclustersize_cosL"), averageClusterSize * cosL, track.p());
          registryData.fill(HIST("omega_neg_avgclustersize_sinL"), averageClusterSize * sinL, track.p());
          registryData.fill(HIST("omega_mass_neg"), track.p(), trackedCascade.omegaMass());
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
