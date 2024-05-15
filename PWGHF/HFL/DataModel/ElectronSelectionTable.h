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

/// \file ElectronSelectionTable.h
/// \author Rashi Gupta <rashi.gupta@cern.ch>, IIT Indore
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore

#ifndef PWGHF_DATAMODEL_ELECTRONSELECTIONTABLE_H_
#define PWGHF_DATAMODEL_ELECTRONSELECTIONTABLE_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
// definition of columns and tables forElectron Selection
namespace hf_sel_electron
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);          //! CollisionID of Electrons
DECLARE_SOA_INDEX_COLUMN(Track, track);                  //! Track ID of Electrons
DECLARE_SOA_COLUMN(TrackEta, tracketa, float);           //!  Track Eta of Electrons
DECLARE_SOA_COLUMN(TrackPhi, trackphi, float);           //! Track Phi of Electrons
DECLARE_SOA_COLUMN(TrackPt, trackpt, float);             //! Track Transverse momentum of Electron
DECLARE_SOA_COLUMN(TrackP, trackp, float);               //!  Track momentum of Electron
DECLARE_SOA_COLUMN(TrackRapidity, trackrapidity, float); //!  Rapadity of Electron
DECLARE_SOA_COLUMN(TrackDcaXY, trackdcaXY, float);       //!  dca of Electron in XY direction
DECLARE_SOA_COLUMN(TrackDcaZ, trackdcaZ, float);         //!  dca of Electron in z direction

DECLARE_SOA_COLUMN(TrackTPCNSigmaE, trackTPCNSigmaE, float); //! NSigma electron (TPC PID)
DECLARE_SOA_COLUMN(TrackTOFNSigmaE, trackTOFNSigmaE, float); //! NSigma electron (TOF PID)

// cluster values
DECLARE_SOA_COLUMN(ClusterEnergy, clusterE, float);        //! cluster energy (GeV)
DECLARE_SOA_COLUMN(ClusterEta, clusterEta, float);         //! cluster pseudorapidity (calculated using vertex)
DECLARE_SOA_COLUMN(ClusterPhi, clusterPhi, float);         //! cluster azimuthal angle (calculated using vertex)
DECLARE_SOA_COLUMN(M02, m02, float);                       //! shower shape long axis
DECLARE_SOA_COLUMN(M20, m20, float);                       //! shower shape short axis
DECLARE_SOA_COLUMN(ClusterNCells, clusterNCells, uint8_t); //! number of cells in cluster
DECLARE_SOA_COLUMN(ClusterTime, clusterTime, float);       //! cluster time (ns)

DECLARE_SOA_COLUMN(DeltaEtaMatch, deltaEtaMatch, float); //! dEta  matched track  to calorimeter
DECLARE_SOA_COLUMN(DeltaPhiMatch, deltaPhiMatch, float); //! dPhi  matched track to calorimeter
DECLARE_SOA_COLUMN(ISEMcal, isEMcal, bool);

} // namespace hf_sel_electron
DECLARE_SOA_TABLE(HfSelEl, "AOD", "HFSELEL", //! Electron Informations
                  o2::soa::Index<>,
                  hf_sel_electron::CollisionId,
                  hf_sel_electron::TrackId,
                  hf_sel_electron::TrackEta,
                  hf_sel_electron::TrackPhi,
                  hf_sel_electron::TrackPt,
                  hf_sel_electron::TrackP,
                  hf_sel_electron::TrackRapidity,
                  hf_sel_electron::TrackDcaXY,
                  hf_sel_electron::TrackDcaZ,
                  hf_sel_electron::TrackTPCNSigmaE,
                  hf_sel_electron::TrackTOFNSigmaE,
                  hf_sel_electron::ClusterEnergy,
                  hf_sel_electron::ClusterEta,
                  hf_sel_electron::ClusterPhi,
                  hf_sel_electron::M02,
                  hf_sel_electron::M20,
                  hf_sel_electron::ClusterNCells,
                  hf_sel_electron::ClusterTime,
                  hf_sel_electron::DeltaEtaMatch,
                  hf_sel_electron::DeltaPhiMatch,
                  hf_sel_electron::ISEMcal);

} // namespace o2::aod

#endif // PWGHF_DATAMODEL_ELECTRONSELECTIONTABLE_H_

