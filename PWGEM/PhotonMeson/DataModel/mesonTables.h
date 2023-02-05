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

#include "Framework/AnalysisDataModel.h"

#ifndef PWGEM_PHOTONMESON_DATAMODEL_MESONTABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_MESONTABLES_H_

namespace o2::aod
{
namespace calomeson
{
DECLARE_SOA_INDEX_COLUMN(Collision,
                         collision); //! Collision to which this meson belongs
DECLARE_SOA_INDEX_COLUMN_FULL(DaugtherPhotonOne, daugtherPhotonOne, int,
                              SkimEMCClusters,
                              "_One"); //! 1st daughter if this meson
DECLARE_SOA_INDEX_COLUMN_FULL(DaugtherPhotonTwo, daugtherPhotonTwo, int,
                              SkimEMCClusters,
                              "_Two"); //! 2nd daughter of this meson
DECLARE_SOA_COLUMN(Oa, oa, float);     //! opening angle between the two daugthers
DECLARE_SOA_COLUMN(Px, px, float);     //! px
DECLARE_SOA_COLUMN(Py, py, float);     //! py
DECLARE_SOA_COLUMN(Pz, pz, float);     //! pz
DECLARE_SOA_COLUMN(E, e, float);       //! E
DECLARE_SOA_COLUMN(Alpha,
                   alpha, float);      //! energy asymmertry of the two daughter particles
DECLARE_SOA_COLUMN(Minv, minv, float); //! invariant mass of the meson
DECLARE_SOA_COLUMN(Eta, eta, float);   //! pseudorapidity of the meson
DECLARE_SOA_COLUMN(Phi, phi, float);   //! phi angle of the meson
DECLARE_SOA_COLUMN(Pt, pt, float);     //! pT of the meson
} // namespace calomeson
DECLARE_SOA_TABLE(CaloMeson, "AOD", "CALOMESON", //!
                  o2::soa::Index<>, calomeson::CollisionId, calomeson::DaugtherPhotonOneId, calomeson::DaugtherPhotonTwoId,
                  calomeson::Oa, calomeson::Px, calomeson::Py, calomeson::Pz, calomeson::E, calomeson::Alpha,
                  calomeson::Minv, calomeson::Eta, calomeson::Phi, calomeson::Pt);
} // namespace o2::aod
#endif // PWGEM_PHOTONMESON_DATAMODEL_MESONTABLES_H_
