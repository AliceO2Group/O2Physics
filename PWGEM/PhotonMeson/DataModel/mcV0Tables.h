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

#include "Framework/AnalysisDataModel.h"

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

// Track Pos
DECLARE_SOA_COLUMN(RecoPosPt, recoPosPt, float);   //! Reconstructed pt of positive sign track
DECLARE_SOA_COLUMN(RecoPosEta, recoPosEta, float); //! Reconstructed eta of positive sign track
DECLARE_SOA_COLUMN(RecoPosTgl, recoPosTgl, float); //! Reconstructed tgl of positive sign track
DECLARE_SOA_COLUMN(RecoPosX, recoPosX, float);     //! Reconstructed X of positive sign track
DECLARE_SOA_COLUMN(RecoPosY, recoPosY, float);     //! Reconstructed X of positive sign track
DECLARE_SOA_COLUMN(RecoPosZ, recoPosZ, float);     //! Reconstructed X of positive sign track
// Track Ele
DECLARE_SOA_COLUMN(RecoElePt, recoElePt, float);   //! Reconstructed pt of negative sign track
DECLARE_SOA_COLUMN(RecoEleTgl, recoEleTgl, float); //! Reconstructed tgl of negative sign track
DECLARE_SOA_COLUMN(RecoEleEta, recoEleEta, float); //! Reconstructed eta of negative sign track
DECLARE_SOA_COLUMN(RecoEleX, recoEleX, float);     //! Reconstructed X of negative sign track
DECLARE_SOA_COLUMN(RecoEleY, recoEleY, float);     //! Reconstructed X of negative sign track
DECLARE_SOA_COLUMN(RecoEleZ, recoEleZ, float);     //! Reconstructed X of negative sign track
// MC particle pos
DECLARE_SOA_COLUMN(MCPosPt, mcPosPt, float);   //! Generated pt of positive sign particle
DECLARE_SOA_COLUMN(MCPosEta, mcPosEta, float); //! Generated eta of positive sign particle
DECLARE_SOA_COLUMN(MCPosX, mcPosX, float);     //! Generated x of positive sign particle
DECLARE_SOA_COLUMN(MCPosY, mcPosY, float);     //! Generated y of positive sign particle
DECLARE_SOA_COLUMN(MCPosZ, mcPosZ, float);     //! Generated z of positive sign particle
// MC particle ele
DECLARE_SOA_COLUMN(MCElePt, mcElePt, float);   //! Generated pt of negative sign particle
DECLARE_SOA_COLUMN(MCEleEta, mcEleEta, float); //! Generated eta of negative sign particle
DECLARE_SOA_COLUMN(MCEleX, mcEleX, float);     //! Generated x of negative sign particle
DECLARE_SOA_COLUMN(MCEleY, mcEleY, float);     //! Generated y of negative sign particle
DECLARE_SOA_COLUMN(MCEleZ, mcEleZ, float);     //! Generated z of negative sign particle
// Propagated MC position
DECLARE_SOA_COLUMN(MCPosPropX, mcPosPropX, float); //! Propagated generated x of positive sign particle
DECLARE_SOA_COLUMN(MCPosPropY, mcPosPropY, float); //! Propagated generated y of positive sign particle
DECLARE_SOA_COLUMN(MCPosPropZ, mcPosPropZ, float); //! Propagated generated z of positive sign particle
DECLARE_SOA_COLUMN(MCElePropX, mcElePropX, float); //! Propagated generated x of negative sign particle
DECLARE_SOA_COLUMN(MCElePropY, mcElePropY, float); //! Propagated generated y of negative sign particle
DECLARE_SOA_COLUMN(MCElePropZ, mcElePropZ, float); //! Propagated generated z of negative sign particle
// TrackType
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t); //! Track type, see enum TrackTypes
} // namespace mcV0

DECLARE_SOA_TABLE(MCV0, "AOD", "MCV0", o2::soa::Index<>,
                  // Track Pos
                  mcV0::RecoPosPt, mcV0::RecoPosEta, mcV0::RecoPosTgl,
                  mcV0::RecoPosX, mcV0::RecoPosY, mcV0::RecoPosZ,
                  // Track Ele
                  mcV0::RecoElePt, mcV0::RecoEleEta, mcV0::RecoEleTgl,
                  mcV0::RecoEleX, mcV0::RecoEleY, mcV0::RecoEleZ,
                  // MC particle pos
                  mcV0::MCPosPt, mcV0::MCPosEta,
                  mcV0::MCPosX, mcV0::MCPosY, mcV0::MCPosZ,
                  // MC particle ele
                  mcV0::MCElePt, mcV0::MCEleEta,
                  mcV0::MCEleX, mcV0::MCEleY, mcV0::MCEleZ,
                  // Propagated MC position
                  mcV0::MCPosPropX, mcV0::MCPosPropY, mcV0::MCPosPropZ,
                  mcV0::MCElePropX, mcV0::MCElePropY, mcV0::MCElePropZ,
                  // Track type
                  mcV0::TrackType);
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_MCV0TABLES_H_
