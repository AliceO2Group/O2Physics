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
//
/// \file trackedCascadeProperties.cxx
///
/// \brief task to study the average cluster size of tracked cascades
///
/// \author Alberto Caliva (alberto.caliva@cern.ch), Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since May 31, 2024

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"

#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using std::array;

// Define type aliases for joined tables
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using FullTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

struct TrackedCascadeProperties {

  // Instantiate the CCDB manager service and API interface
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  // Instantiate the main Zorro processing object and define an output to store summary information
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Histogram registry for quality control
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Histogram registry for data
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Global Parameters
  Configurable<float> zVtx{"zVtx", 10.0f, "z vertex cut"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
  Configurable<std::string> triggerList{"triggerList", "fTrackedOmega, fTrackedXi, fOmegaLargeRadius, fDoubleOmega, fOmegaHighMult, fSingleXiYN, fQuadrupleXi, fDoubleXi, fhadronOmega, fOmegaXi, fTripleXi, fOmega", "Trigger list"};

  // Analysis Selections
  Configurable<int> minItsClustersCasc{"minItsClustersCasc", 4, "min ITS Clusters"};
  Configurable<float> massMinXi{"massMinXi", 1.315f, "mMin Xi"};
  Configurable<float> massMaxXi{"massMaxXi", 1.328f, "mMax Xi"};
  Configurable<float> massMinOmega{"massMinOmega", 1.665f, "mMin Omega"};
  Configurable<float> massMaxOmega{"massMaxOmega", 1.680f, "mMax Omega"};

  // Initialize CCDB access and histogram registry for Zorro processing
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
      zorro.populateHistRegistry(registryData, bc.runNumber());
    }
  }

  void init(InitContext const&)
  {
    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    // Quality Control Histograms
    registryQC.add("matchingChi2", "matching Chi2", HistType::kTH1F, {{200, 0, 1000, "#chi^{2}_{matching}"}});
    registryQC.add("topologyChi2", "topology Chi2", HistType::kTH1F, {{500, 0, 0.5, "#chi^{2}_{topology}"}});
    registryQC.add("nITScls_vs_p_xi", "nITS Xi", HistType::kTH2F, {{100, 0, 10, "#it{p} (GeV/#it{c})"}, {8, 0, 8, "n_{ITS}^{cls}"}});
    registryQC.add("nITScls_vs_p_omega", "nITS Omega", HistType::kTH2F, {{100, 0, 10, "#it{p} (GeV/#it{c})"}, {8, 0, 8, "n_{ITS}^{cls}"}});
    registryQC.add("decayXY", "decayXY", HistType::kTH2F, {{500, -50, 50, "x"}, {500, -50, 50, "y"}});

    // Event Counter
    registryData.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{5, 0, 5, "Event Cuts"}});

    // Average cluster size vs momentum
    registryData.add("xi_pos_avgclustersize_cosL_vs_p", "xi_pos_avgclustersize_cosL_vs_p", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("xi_neg_avgclustersize_cosL_vs_p", "xi_neg_avgclustersize_cosL_vs_p", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_pos_avgclustersize_cosL_vs_p", "omega_pos_avgclustersize_cosL_vs_p", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_neg_avgclustersize_cosL_vs_p", "omega_neg_avgclustersize_cosL_vs_p", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});

    // Average cluster size vs betagamma
    registryData.add("xi_pos_avgclustersize_cosL_vs_betagamma", "xi_pos_avgclustersize_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("xi_neg_avgclustersize_cosL_vs_betagamma", "xi_neg_avgclustersize_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_pos_avgclustersize_cosL_vs_betagamma", "omega_pos_avgclustersize_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_neg_avgclustersize_cosL_vs_betagamma", "omega_neg_avgclustersize_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});

    // Cluster size using truncated mean vs momentum
    registryData.add("xi_pos_avgclustersize_trunc_cosL_vs_p", "xi_pos_avgclustersize_trunc_cosL_vs_p", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("xi_neg_avgclustersize_trunc_cosL_vs_p", "xi_neg_avgclustersize_trunc_cosL_vs_p", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_pos_avgclustersize_trunc_cosL_vs_p", "omega_pos_avgclustersize_trunc_cosL_vs_p", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_neg_avgclustersize_trunc_cosL_vs_p", "omega_neg_avgclustersize_trunc_cosL_vs_p", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});

    // Cluster size using truncated mean vs betagamma
    registryData.add("xi_pos_avgclustersize_trunc_cosL_vs_betagamma", "xi_pos_avgclustersize_trunc_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("xi_neg_avgclustersize_trunc_cosL_vs_betagamma", "xi_neg_avgclustersize_trunc_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_pos_avgclustersize_trunc_cosL_vs_betagamma", "omega_pos_avgclustersize_trunc_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});
    registryData.add("omega_neg_avgclustersize_trunc_cosL_vs_betagamma", "omega_neg_avgclustersize_trunc_cosL_vs_betagamma", HistType::kTH2F, {{200, 0.0, 10.0, "#beta#gamma"}, {100, 0.0, 20.0, "#LT ITS cluster size #GT cos(#lambda)"}});

    // mass histograms
    registryData.add("xi_mass_pos", "xi_mass_pos", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("xi_mass_neg", "xi_mass_neg", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("omega_mass_pos", "omega_mass_pos", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
    registryData.add("omega_mass_neg", "omega_mass_neg", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
  }

  double trackInclination(double eta)
  {
    double lambda(0);
    double theta = 2.0 * std::atan(std::exp(-eta));
    if (theta <= PIHalf)
      lambda = PIHalf - theta;
    if (theta > PIHalf)
      lambda = theta - PIHalf;
    return lambda;
  }

  int findBin(const std::vector<double>& edges, double value)
  {
    auto it = std::upper_bound(edges.begin(), edges.end(), value);
    int index = static_cast<int>(it - edges.begin()) - 1;
    if (index < 0 || index >= static_cast<int>(edges.size()) - 1) {
      return -1; // value is out of bounds
    }
    return index;
  }

  void processData(SelectedCollisions::iterator const& collision, aod::AssignedTrackedCascades const& trackedCascades,
                   aod::Cascades const&, FullTracks const&, aod::BCsWithTimestamps const&)
  {
    // Number of events before any selection
    registryData.fill(HIST("number_of_events_data"), 0.5);

    // Retrieve the bunch crossing information with timestamps from the collision
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    // If skimmed processing is enabled, apply Zorro trigger selection
    if (cfgSkimmedProcessing && !zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC())) {
      return;
    }

    // Number of events after skimming selection
    registryData.fill(HIST("number_of_events_data"), 1.5);
    if (!collision.sel8())
      return;

    // Number of events after sel8
    registryData.fill(HIST("number_of_events_data"), 2.5);
    if (std::abs(collision.posZ()) > zVtx)
      return;

    // Number of events after zVtx cut
    registryData.fill(HIST("number_of_events_data"), 3.5);

    // radii of ITS layers
    std::vector<double> edgesItsLayers = {0.0, 2.2, 2.8, 3.6, 20.0, 22.0, 37.0, 39.0, 100.0};

    // Loop over tracked cascades
    for (const auto& trackedCascade : trackedCascades) {

      // Get tracked cascade
      const auto track = trackedCascade.track_as<FullTracks>();
      const auto trackITS = trackedCascade.itsTrack_as<FullTracks>();
      const auto& casc = trackedCascade.cascade();
      const auto& btrack = casc.bachelor_as<FullTracks>();
      double dx = trackedCascade.decayX();
      double dy = trackedCascade.decayY();
      double r = std::sqrt(dx * dx + dy * dy);
      int nItsLayersCrossed = findBin(edgesItsLayers, r);

      // Fill QC histograms
      registryQC.fill(HIST("matchingChi2"), trackedCascade.matchingChi2());
      registryQC.fill(HIST("topologyChi2"), trackedCascade.topologyChi2());
      registryQC.fill(HIST("decayXY"), dx, dy);

      // Compute average cluster size and truncated mean
      double sumClusterSize = 0.0;
      double sumClusterSizeTrunc = 0.0;
      double maxClusterSize = 0.0;
      double averageClusterSize = 0.0;
      double averageClusterSizeTrunc = 0.0;
      int nCls = 0;

      for (int i = 0; i < nItsLayersCrossed; i++) {
        double clusterSize = static_cast<double>(trackITS.itsClsSizeInLayer(i));
        if (clusterSize > 0) {
          sumClusterSize += clusterSize;
          sumClusterSizeTrunc += clusterSize;
          nCls++;
          if (clusterSize > maxClusterSize) {
            maxClusterSize = clusterSize;
          }
        }
      }
      if (nCls > 0) {
        averageClusterSize = sumClusterSize / static_cast<double>(nCls);
      }
      if (nCls > 1) {
        averageClusterSizeTrunc = (sumClusterSizeTrunc - maxClusterSize) / static_cast<double>(nCls - 1);
      }

      // Apply selection on number of ITS clusters
      if (nCls < minItsClustersCasc)
        continue;

      // Xi Mass
      if (btrack.sign() > 0) {
        registryData.fill(HIST("xi_mass_pos"), track.p(), trackedCascade.xiMass());
      }
      if (btrack.sign() < 0) {
        registryData.fill(HIST("xi_mass_neg"), track.p(), trackedCascade.xiMass());
      }

      // Variables
      double lambda = trackInclination(track.eta());
      double clsSizeCosL = averageClusterSize * std::cos(lambda);
      double clsSizeCosLtrunc = averageClusterSizeTrunc * std::cos(lambda);
      double bgXi = track.p() / MassXiPlusBar;
      double bgOmega = track.p() / MassOmegaPlusBar;

      // Xi
      if (trackedCascade.xiMass() > massMinXi && trackedCascade.xiMass() < massMaxXi) {
        registryQC.fill(HIST("nITScls_vs_p_xi"), track.p(), trackITS.itsNCls());
        if (btrack.sign() > 0) {
          registryData.fill(HIST("xi_pos_avgclustersize_cosL_vs_p"), track.p(), clsSizeCosL);
          registryData.fill(HIST("xi_pos_avgclustersize_cosL_vs_betagamma"), bgXi, clsSizeCosL);
          registryData.fill(HIST("xi_pos_avgclustersize_trunc_cosL_vs_p"), track.p(), clsSizeCosLtrunc);
          registryData.fill(HIST("xi_pos_avgclustersize_trunc_cosL_vs_betagamma"), bgXi, clsSizeCosLtrunc);
        }
        if (btrack.sign() < 0) {
          registryData.fill(HIST("xi_neg_avgclustersize_cosL_vs_p"), track.p(), clsSizeCosL);
          registryData.fill(HIST("xi_neg_avgclustersize_cosL_vs_betagamma"), bgXi, clsSizeCosL);
          registryData.fill(HIST("xi_neg_avgclustersize_trunc_cosL_vs_p"), track.p(), clsSizeCosLtrunc);
          registryData.fill(HIST("xi_neg_avgclustersize_trunc_cosL_vs_betagamma"), bgXi, clsSizeCosLtrunc);
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
      if (trackedCascade.omegaMass() > massMinOmega && trackedCascade.omegaMass() < massMaxOmega) {
        registryQC.fill(HIST("nITScls_vs_p_omega"), track.p(), trackITS.itsNCls());
        if (btrack.sign() > 0) {
          registryData.fill(HIST("omega_pos_avgclustersize_cosL_vs_p"), track.p(), clsSizeCosL);
          registryData.fill(HIST("omega_pos_avgclustersize_cosL_vs_betagamma"), bgOmega, clsSizeCosL);
          registryData.fill(HIST("omega_pos_avgclustersize_trunc_cosL_vs_p"), track.p(), clsSizeCosLtrunc);
          registryData.fill(HIST("omega_pos_avgclustersize_trunc_cosL_vs_betagamma"), bgOmega, clsSizeCosLtrunc);
        }
        if (btrack.sign() < 0) {
          registryData.fill(HIST("omega_neg_avgclustersize_cosL_vs_p"), track.p(), clsSizeCosL);
          registryData.fill(HIST("omega_neg_avgclustersize_cosL_vs_betagamma"), bgOmega, clsSizeCosL);
          registryData.fill(HIST("omega_neg_avgclustersize_trunc_cosL_vs_p"), track.p(), clsSizeCosLtrunc);
          registryData.fill(HIST("omega_neg_avgclustersize_trunc_cosL_vs_betagamma"), bgOmega, clsSizeCosLtrunc);
        }
      }
    }
  }
  PROCESS_SWITCH(TrackedCascadeProperties, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrackedCascadeProperties>(cfgc)};
}
