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
#include <TMath.h>

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
};
} // namespace pwgem::dilepton

namespace emprimarytrack
{
DECLARE_SOA_COLUMN(PIDLabel, pidlabel, int);                     //!
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int);             //!
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int); //!
} // namespace emprimarytrack

// reconstructed track information
DECLARE_SOA_TABLE(EMPrimaryTracks, "AOD", "EMPTRACK", //!
                  o2::soa::Index<>, collision::PosZ, collision::NumContrib,
                  track::Pt, track::Eta, track::Phi, track::Tgl, track::Signed1Pt,
                  track::DcaXY, track::DcaZ, track::CYY, track::CZZ, track::CZY,
                  track::TPCNClsFindable, emprimarytrack::TPCNClsFound, emprimarytrack::TPCNClsCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap, emprimarytrack::PIDLabel);

// iterators
using EMPrimaryTrack = EMPrimaryTracks::iterator;

} // namespace o2::aod

#endif // PWGEM_DILEPTON_DATAMODEL_LMEEMLTABLES_H_
