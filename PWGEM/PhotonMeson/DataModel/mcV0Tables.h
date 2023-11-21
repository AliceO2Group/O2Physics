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

/// \brief Tables for investigating the MC V0 phase-space
/// \author felix.schlepper@cern.ch

#ifndef PWGEM_PHOTONMESON_DATAMODEL_MCV0TABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_MCV0TABLES_H_

#include "Framework/ASoA.h"

namespace o2::aod
{

namespace mcV0
{
enum TrackTypes : uint8_t {
  ITSTPC_ITSTPC = 0,
  TPConly_TPConly,
  ITSonly_ITSonly,
  ITSTPC_TPConly,
  ITSTPC_ITSonly,
  TPConly_ITSonly,
};

// V0
DECLARE_SOA_COLUMN(V0X, v0X, float);   //! Reconstructed x position of V0
DECLARE_SOA_COLUMN(V0Y, v0Y, float);   //! Reconstructed y position of V0
DECLARE_SOA_COLUMN(V0Z, v0z, float);   //! Reconstructed z position of V0
DECLARE_SOA_COLUMN(V0Px, v0Px, float); //! Reconstructed Px of V0
DECLARE_SOA_COLUMN(V0Py, v0Py, float); //! Reconstructed Py of V0
DECLARE_SOA_COLUMN(V0Pz, v0Pz, float); //! Reconstructed Pz of V0
// Track Pos
DECLARE_SOA_COLUMN(RecoPosPt, recoPosPt, float);        //! Reconstructed pt of pos track
DECLARE_SOA_COLUMN(RecoPosEta, recoPosEta, float);      //! Reconstructed eta of pos track
DECLARE_SOA_COLUMN(RecoPosTgl, recoPosTgl, float);      //! Reconstructed tgl of pos track
DECLARE_SOA_COLUMN(RecoPosX, recoPosX, float);          //! Reconstructed X of pos track
DECLARE_SOA_COLUMN(RecoPosY, recoPosY, float);          //! Reconstructed X of pos track
DECLARE_SOA_COLUMN(RecoPosZ, recoPosZ, float);          //! Reconstructed X of pos track
DECLARE_SOA_COLUMN(RecoPosT, recoPosT, float);          //! Estimated time of the pos track in ns wrt collision()
DECLARE_SOA_COLUMN(RecoPosSign, recoPosSign, int16_t);  //! Sign of pos track
DECLARE_SOA_COLUMN(RecoPosHasITS, recoPosHasITS, bool); //! if pos track has ITS
DECLARE_SOA_COLUMN(RecoPosHasTPC, recoPosHasTPC, bool); //! if pos track has TPC
DECLARE_SOA_COLUMN(RecoPosHasTRD, recoPosHasTRD, bool); //! if pos track has TRD
DECLARE_SOA_COLUMN(RecoPosHasTOF, recoPosHasTOF, bool); //! if pos track has TOF
// Track Ele
DECLARE_SOA_COLUMN(RecoElePt, recoElePt, float);        //! Reconstructed pt of ele track
DECLARE_SOA_COLUMN(RecoEleTgl, recoEleTgl, float);      //! Reconstructed tgl of ele track
DECLARE_SOA_COLUMN(RecoEleEta, recoEleEta, float);      //! Reconstructed eta of ele track
DECLARE_SOA_COLUMN(RecoEleX, recoEleX, float);          //! Reconstructed X of ele track
DECLARE_SOA_COLUMN(RecoEleY, recoEleY, float);          //! Reconstructed X of ele track
DECLARE_SOA_COLUMN(RecoEleZ, recoEleZ, float);          //! Reconstructed X of ele track
DECLARE_SOA_COLUMN(RecoEleT, recoEleT, float);          //! Estimated time of the ele track in ns wrt collision()
DECLARE_SOA_COLUMN(RecoEleSign, recoEleSign, int16_t);  //! Sign of ele track
DECLARE_SOA_COLUMN(RecoEleHasITS, recoEleHasITS, bool); //! If ele track has ITS
DECLARE_SOA_COLUMN(RecoEleHasTPC, recoEleHasTPC, bool); //! If ele track has TPC
DECLARE_SOA_COLUMN(RecoEleHasTRD, recoEleHasTRD, bool); //! If ele track has TRD
DECLARE_SOA_COLUMN(RecoEleHasTOF, recoEleHasTOF, bool); //! If ele track has TOF
// MC particle pos
DECLARE_SOA_COLUMN(MCPosPt, mcPosPt, float);   //! Generated pt of pos particle
DECLARE_SOA_COLUMN(MCPosEta, mcPosEta, float); //! Generated eta of pos particle
DECLARE_SOA_COLUMN(MCPosX, mcPosX, float);     //! Generated x of pos particle
DECLARE_SOA_COLUMN(MCPosY, mcPosY, float);     //! Generated y of pos particle
DECLARE_SOA_COLUMN(MCPosZ, mcPosZ, float);     //! Generated z of pos particle
// MC particle ele
DECLARE_SOA_COLUMN(MCElePt, mcElePt, float);   //! Generated pt of ele particle
DECLARE_SOA_COLUMN(MCEleEta, mcEleEta, float); //! Generated eta of ele particle
DECLARE_SOA_COLUMN(MCEleX, mcEleX, float);     //! Generated x of ele particle
DECLARE_SOA_COLUMN(MCEleY, mcEleY, float);     //! Generated y of ele particle
DECLARE_SOA_COLUMN(MCEleZ, mcEleZ, float);     //! Generated z of ele particle
// Propagated MC pos
DECLARE_SOA_COLUMN(MCPosPropX, mcPosPropX, float);     //! Propagated generated x of pos particle
DECLARE_SOA_COLUMN(MCPosPropY, mcPosPropY, float);     //! Propagated generated y of pos particle
DECLARE_SOA_COLUMN(MCPosPropZ, mcPosPropZ, float);     //! Propagated generated z of pos particle
DECLARE_SOA_COLUMN(MCPosPropEta, mcPosPropEta, float); //! Propagated generated eta of pos particle
DECLARE_SOA_COLUMN(MCPosPropTgl, mcPosPropTgl, float); //! Propagated generated tgl of pos particle
DECLARE_SOA_COLUMN(MCPosPropPt, mcPosPropPt, float);   //! Propagated generated Pt of pos particle
// Propagated MC ele
DECLARE_SOA_COLUMN(MCElePropX, mcElePropX, float);     //! Propagated generated x of ele particle
DECLARE_SOA_COLUMN(MCElePropY, mcElePropY, float);     //! Propagated generated y of ele particle
DECLARE_SOA_COLUMN(MCElePropZ, mcElePropZ, float);     //! Propagated generated z of ele particle
DECLARE_SOA_COLUMN(MCElePropEta, mcElePropEta, float); //! Propagated generated eta of ele particle
DECLARE_SOA_COLUMN(MCElePropTgl, mcElePropTgl, float); //! Propagated generated tgl of ele particle
DECLARE_SOA_COLUMN(MCElePropPt, mcElePropPt, float);   //! Propagated generated Pt of ele particle
// MC mother particle
DECLARE_SOA_COLUMN(MCMotherIsPrimary, mcMotherIsPrimary, bool);     //! If MC mother is primary
DECLARE_SOA_COLUMN(MCMotherIsGenerator, mcMotherIsGenerator, bool); //! If MC mother is generated
DECLARE_SOA_COLUMN(MCMotherPt, mcMotherPt, float);                  //! Mother particle Pt
DECLARE_SOA_COLUMN(MCMotherEta, mcMotherEta, float);                //! Mother particle eta
DECLARE_SOA_COLUMN(MCMotherX, mcMotherX, float);                    //! Mother particle vertex x
DECLARE_SOA_COLUMN(MCMotherY, mcMotherY, float);                    //! Mother particle vertex y
DECLARE_SOA_COLUMN(MCMotherZ, mcMotherZ, float);                    //! Mother particle vertex z
// MC Primary Vertex
DECLARE_SOA_COLUMN(MCPVtxX, mcPVtxX, float); //! MC Primary vertex X
DECLARE_SOA_COLUMN(MCPVtxY, mcPVtxY, float); //! MC Primary vertex Y
DECLARE_SOA_COLUMN(MCPVtxZ, mcPVtxZ, float); //! MC Primary vertex Z
// TrackType
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t); //! Track type, see enum TrackTypes
} // namespace mcV0

DECLARE_SOA_TABLE(MCV0, "AOD", "MCV0", o2::soa::Index<>,
                  // V0
                  mcV0::V0X, mcV0::V0Y, mcV0::V0Z,
                  mcV0::V0Px, mcV0::V0Py, mcV0::V0Pz,
                  /* mcV0::V0Pt<mcV0::V0Px, mcV0::V0Py>, */
                  // Track Pos
                  mcV0::RecoPosPt, mcV0::RecoPosEta, mcV0::RecoPosTgl,
                  mcV0::RecoPosX, mcV0::RecoPosY, mcV0::RecoPosZ,
                  mcV0::RecoPosT, mcV0::RecoPosSign,
                  mcV0::RecoPosHasITS, mcV0::RecoPosHasTPC, mcV0::RecoPosHasTRD, mcV0::RecoPosHasTOF,
                  // Track Ele
                  mcV0::RecoElePt, mcV0::RecoEleEta, mcV0::RecoEleTgl,
                  mcV0::RecoEleX, mcV0::RecoEleY, mcV0::RecoEleZ,
                  mcV0::RecoEleT, mcV0::RecoEleSign,
                  mcV0::RecoEleHasITS, mcV0::RecoEleHasTPC, mcV0::RecoEleHasTRD, mcV0::RecoEleHasTOF,
                  // MC particle pos
                  mcV0::MCPosPt, mcV0::MCPosEta,
                  mcV0::MCPosX, mcV0::MCPosY, mcV0::MCPosZ,
                  // MC particle ele
                  mcV0::MCElePt, mcV0::MCEleEta,
                  mcV0::MCEleX, mcV0::MCEleY, mcV0::MCEleZ,
                  // Propagated MC pos
                  mcV0::MCPosPropX, mcV0::MCPosPropY, mcV0::MCPosPropZ,
                  mcV0::MCPosPropEta, mcV0::MCPosPropTgl, mcV0::MCPosPropPt,
                  // Propagated MC pos
                  mcV0::MCElePropX, mcV0::MCElePropY, mcV0::MCElePropZ,
                  mcV0::MCElePropEta, mcV0::MCElePropTgl, mcV0::MCElePropPt,
                  // MC mother particle
                  mcV0::MCMotherIsPrimary, mcV0::MCMotherIsGenerator, mcV0::MCMotherPt, mcV0::MCMotherEta,
                  mcV0::MCMotherX, mcV0::MCMotherY, mcV0::MCMotherZ,
                  // Primary Vertex
                  mcV0::MCPVtxX, mcV0::MCPVtxY, mcV0::MCPVtxZ,
                  // Track type
                  mcV0::TrackType);
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_MCV0TABLES_H_
