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
#ifndef TUTORIALS_SRC_TRAININGTREE_H_
#define TUTORIALS_SRC_TRAININGTREE_H_
#include "Framework/AnalysisDataModel.h"

/// Definition of output table to store the flat tree for re-weighting model
/// training
namespace o2::aod
{
namespace vars
{
DECLARE_SOA_COLUMN(Zvtx, vertexZ, float);           //! Vertex Z
DECLARE_SOA_COLUMN(Xvtx, vertexX, float);           //! Vertex X
DECLARE_SOA_COLUMN(Yvtx, vertexY, float);           //! Vertex Y
DECLARE_SOA_COLUMN(MeanPt, meanPt, float);          //! per-collision mean track Pt
DECLARE_SOA_COLUMN(BarrelMult, centralMult, float); //! Central multiplicity of a collison
DECLARE_SOA_COLUMN(FMT0, forwardMultT0, float);     //! Forward multiplicity of a collision (T0)
DECLARE_SOA_COLUMN(FMV0, forwardMultV0, float);     //! Forward multiplicity of a collision (V0)
} // namespace vars

DECLARE_SOA_TABLE(TrainingTree, "AOD", "TRTR", o2::soa::Index<>,
                  vars::Zvtx,
                  vars::Xvtx,
                  vars::Yvtx,
                  vars::MeanPt,
                  vars::BarrelMult,
                  vars::FMT0,
                  vars::FMV0);
} // namespace o2::aod

namespace o2::analysis
{
template <typename Tracks>
auto meanPt(Tracks const& tracks)
{
  auto apt = 0.f;
  auto npt = 0;
  for (auto& track : tracks) {
    if (isfinite(track.pt()) && (std::abs(track.pt()) > 1e-3)) {
      ++npt;
      apt += track.pt();
    }
  }
  return apt / static_cast<float>(npt);
}
} // namespace o2::analysis

#endif // TUTORIALS_SRC_TRAININGTREE_H_
