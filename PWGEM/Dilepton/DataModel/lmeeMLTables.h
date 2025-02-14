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

#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#ifndef PWGEM_DILEPTON_DATAMODEL_LMEEMLTABLES_H_
#define PWGEM_DILEPTON_DATAMODEL_LMEEMLTABLES_H_

namespace o2::aod
{

namespace pwgem::dilepton
{
enum class PID_Label : int {
  kUnDef = -1,
  kElectron = 0,
  kMuon = 1,
  kPion = 2,
  kKaon = 3,
  kProton = 4,
}; // this can be used for eID.

enum class Track_Type : int {
  kPrimary = 0,
  kSecondary = 1,
}; // this can be used for selecting electron from primary or photon conversion.

} // namespace pwgem::dilepton

namespace emprimarytrack
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);               //!
DECLARE_SOA_COLUMN(PIDLabel, pidlabel, int);                     //!
DECLARE_SOA_COLUMN(TrackType, tracktype, int);                   //!
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int);             //!
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int); //!
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITS, meanClusterSizeITS, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
} // namespace emprimarytrack

// reconstructed track information
DECLARE_SOA_TABLE(EMPrimaryTracks, "AOD", "EMPTRACK", //!
                  o2::soa::Index<>, emprimarytrack::CollisionId, collision::PosZ, collision::NumContrib,
                  track::Pt, track::Eta, track::Phi, track::Tgl, track::Signed1Pt,
                  track::DcaXY, track::DcaZ, track::CYY, track::CZZ, track::CZY,
                  track::TPCNClsFindable, emprimarytrack::TPCNClsFound, emprimarytrack::TPCNClsCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::TOFChi2, track::DetectorMap, emprimarytrack::PIDLabel, emprimarytrack::TrackType,

                  // dynamic column
                  emprimarytrack::P<track::Pt, track::Eta>,
                  emprimarytrack::MeanClusterSizeITS<track::ITSClusterSizes>);

// iterators
using EMPrimaryTrack = EMPrimaryTracks::iterator;

} // namespace o2::aod

#endif // PWGEM_DILEPTON_DATAMODEL_LMEEMLTABLES_H_
