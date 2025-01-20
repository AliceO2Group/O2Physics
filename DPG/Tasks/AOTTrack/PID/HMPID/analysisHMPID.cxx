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

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <TTree.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace variables_table // declaration of columns to create
{
DECLARE_SOA_COLUMN(ChAngle, chAngle, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(MomentumHMPID, momentumHMPID, float);
DECLARE_SOA_COLUMN(MomentumTrack, momentumTrack, float);
DECLARE_SOA_COLUMN(Xtrack, xtrack, float);
DECLARE_SOA_COLUMN(Ytrack, ytrack, float);
DECLARE_SOA_COLUMN(Xmip, xmip, float);
DECLARE_SOA_COLUMN(Ymip, ymip, float);
DECLARE_SOA_COLUMN(Nphotons, nphotons, float);
DECLARE_SOA_COLUMN(ChargeMIP, chargeMIP, float);
DECLARE_SOA_COLUMN(ClusterSize, clustersize, float);
DECLARE_SOA_COLUMN(Chamber, chamber, float);
DECLARE_SOA_COLUMN(Photons_charge, photons_charge, float[10]);

DECLARE_SOA_COLUMN(EtaTrack, etatrack, float);
DECLARE_SOA_COLUMN(PhiTrack, phitrack, float);

DECLARE_SOA_COLUMN(ITSNcluster, itsNcluster, float);
DECLARE_SOA_COLUMN(TPCNcluster, tpcNcluster, float);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, float);
DECLARE_SOA_COLUMN(TPCchi2, tpcChi2, float);
DECLARE_SOA_COLUMN(ITSchi2, itsChi2, float);

DECLARE_SOA_COLUMN(DCAxy, dcaxy, float);
DECLARE_SOA_COLUMN(DCAz, dcaz, float);

DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNsigmaPi, float);
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNsigmaPi, float);
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNsigmaKa, float);
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNsigmaKa, float);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNsigmaPr, float);
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNsigmaPr, float);
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNsigmaDe, float);
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNsigmaDe, float);

} // namespace variables_table

DECLARE_SOA_TABLE(HMPID_analysis, "AOD", "HMPIDANALYSIS",
                  variables_table::ChAngle, variables_table::Phi, variables_table::Eta, variables_table::MomentumHMPID,
                  variables_table::MomentumTrack, variables_table::Xtrack, variables_table::Ytrack, variables_table::Xmip,
                  variables_table::Ymip, variables_table::Nphotons, variables_table::ChargeMIP, variables_table::ClusterSize,
                  variables_table::Chamber, variables_table::Photons_charge, variables_table::EtaTrack, variables_table::PhiTrack,
                  variables_table::ITSNcluster, variables_table::TPCNcluster, variables_table::TPCNClsCrossedRows,
                  variables_table::TPCchi2, variables_table::ITSchi2, variables_table::DCAxy, variables_table::DCAz,
                  variables_table::TPCNSigmaPi, variables_table::TOFNSigmaPi, variables_table::TPCNSigmaKa, variables_table::TOFNSigmaKa,
                  variables_table::TPCNSigmaPr, variables_table::TOFNSigmaPr, variables_table::TPCNSigmaDe, variables_table::TOFNSigmaDe);
} // namespace o2::aod

struct pidHmpidAnalysis {

  Produces<aod::HMPID_analysis> HMPID_analysis;

  // using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                    aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;

  void process(const aod::HMPIDs& hmpids,
               TrackCandidates const&,
               CollisionCandidates const&)
  {

    for (const auto& t : hmpids) {
      if (t.track_as<TrackCandidates>().isGlobalTrack() != (uint8_t) true) {
        continue;
      }

      const auto& track = t.track_as<TrackCandidates>();

      if (!track.hasITS() || !track.hasTPC() || !track.hasTOF()) {
        continue;
      }

      float hmpidPhotsCharge2[10];

      for (int i = 0; i < 10; i++) {
        hmpidPhotsCharge2[i] = t.hmpidPhotsCharge()[i];
      }

      /////FILL TABLE
      HMPID_analysis(t.hmpidSignal(), t.track_as<TrackCandidates>().phi(), t.track_as<TrackCandidates>().eta(), t.hmpidMom(),
                     track.p(), t.hmpidXTrack(), t.hmpidYTrack(), t.hmpidXMip(),
                     t.hmpidYMip(), t.hmpidNPhotons(), t.hmpidQMip(), (t.hmpidClusSize() % 1000000) / 1000, t.hmpidClusSize() / 1000000,
                     hmpidPhotsCharge2, track.eta(), track.phi(), track.itsNCls(), track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                     track.tpcChi2NCl(), track.itsChi2NCl(), track.dcaXY(), track.dcaZ(),
                     track.tpcNSigmaPi(), track.tofNSigmaPi(), track.tpcNSigmaKa(), track.tofNSigmaKa(),
                     track.tpcNSigmaPr(), track.tofNSigmaPr(), track.tpcNSigmaDe(), track.tofNSigmaDe());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) { return WorkflowSpec{adaptAnalysisTask<pidHmpidAnalysis>(cfg)}; }
