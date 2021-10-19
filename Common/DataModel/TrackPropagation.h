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

#ifndef O2_ANALYSIS_TRACKPROPAGATION_H_
#define O2_ANALYSIS_TRACKPROPAGATION_H_

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "MathUtils/Utils.h"
#include <cmath>
#include "Framework/DataTypes.h"

namespace o2::aod
{
namespace trackpropagated //
{

// Columns to store user needed track params
DECLARE_SOA_COLUMN(Sign, sign, int8_t); //! sign of the track's charge
DECLARE_SOA_COLUMN(Pt, pt, float);      //! track transverse momentum GeV/c
DECLARE_SOA_COLUMN(Phi, phi, float);    //! track phi
DECLARE_SOA_COLUMN(Eta, eta, float);    //! track eta

DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! x-component of the track momentum GeV/c
                           [](float pt, float phi) -> float {
                             return pt * std::cos(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! y-component of the track momentum GeV/c
                           [](float pt, float phi) -> float {
                             return pt * std::sin(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! z-component of the track momentum GeV/c
                           [](float pt, float eta) -> float {
                             return pt * std::sinh(eta);
                           });

} // namespace trackpropagated

DECLARE_SOA_TABLE(tracksPropagated, "AOD", "TRACKPROPAG", //! commonly used track parameters, propagated to the primary vertex
                  o2::soa::Index<>,
                  track::CollisionId,
                  track::TrackType,
                  trackpropagated::Sign,
                  trackpropagated::Pt,
                  trackpropagated::Phi,
                  trackpropagated::Eta,
                  trackpropagated::Px<trackpropagated::Pt, trackpropagated::Phi>,
                  trackpropagated::Py<trackpropagated::Pt, trackpropagated::Phi>,
                  trackpropagated::Pz<trackpropagated::Pt, trackpropagated::Eta>);

using trackPropagated = tracksPropagated::iterator;

DECLARE_SOA_TABLE(tracksParPropagated, "AOD", "TRACKPARPROPAG", //! additional track parameters, propagated to the primary vertex
                  track::X, track::Alpha,
                  track::Y, track::Z, track::Snp, track::Tgl,
                  track::Signed1Pt,
                  track::Px<track::Signed1Pt, track::Snp, track::Alpha>,
                  track::Py<track::Signed1Pt, track::Snp, track::Alpha>,
                  track::Pz<track::Signed1Pt, track::Tgl>,
                  track::Sign<track::Signed1Pt>);

using trackParPropagated = tracksParPropagated::iterator;
// TODO: replace this table with dynamical columns in the tracksPropagated table

} // namespace o2::aod

#endif // O2_ANALYSIS_TRACKPROPAGATION_H_
