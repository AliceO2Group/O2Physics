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

///
/// \brief Table definitions for background subtraction
///
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_DATAMODEL_JETSUBTRACTION_H_
#define PWGJE_DATAMODEL_JETSUBTRACTION_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

namespace o2::aod
{

namespace bkgrho
{
DECLARE_SOA_COLUMN(Rho, rho, float);   //!
DECLARE_SOA_COLUMN(RhoM, rhoM, float); //!
} // namespace bkgrho

namespace bkgcharged
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
} // namespace bkgcharged

namespace bkgd0
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfD0Bases, "_0");
} // namespace bkgd0

namespace bkglc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfCand3Prong, "_0");
} // namespace bkglc

namespace bkgbplus
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfCandBplus, "_0");
} // namespace bkgbplus

DECLARE_SOA_TABLE(BkgChargedRhos, "AOD", "BkgCRho",
                  o2::soa::Index<>,
                  bkgcharged::JCollisionId,
                  bkgrho::Rho,
                  bkgrho::RhoM);

DECLARE_SOA_TABLE(BkgD0Rhos, "AOD", "BkgD0Rho",
                  o2::soa::Index<>,
                  bkgd0::CandidateId,
                  bkgrho::Rho,
                  bkgrho::RhoM);

DECLARE_SOA_TABLE(BkgLcRhos, "AOD", "BkgLcRho",
                  o2::soa::Index<>,
                  bkglc::CandidateId,
                  bkgrho::Rho,
                  bkgrho::RhoM);

DECLARE_SOA_TABLE(BkgBplusRhos, "AOD", "BkgBPlRho",
                  o2::soa::Index<>,
                  bkgbplus::CandidateId,
                  bkgrho::Rho,
                  bkgrho::RhoM);

namespace jtracksub
{

DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
DECLARE_SOA_COLUMN(TrackSel, trackSel, uint8_t);
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py,
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz,
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p,
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
} // namespace jtracksub

DECLARE_SOA_TABLE(JTrackSubs, "AOD", "JTrackSubs",
                  o2::soa::Index<>,
                  bkgcharged::JCollisionId,
                  jtracksub::Pt,
                  jtracksub::Eta,
                  jtracksub::Phi,
                  jtracksub::Energy,
                  jtracksub::TrackSel,
                  jtracksub::Px<jtracksub::Pt, jtracksub::Phi>,
                  jtracksub::Py<jtracksub::Pt, jtracksub::Phi>,
                  jtracksub::Pz<jtracksub::Pt, jtracksub::Eta>,
                  jtracksub::P<jtracksub::Pt, jtracksub::Eta>);

using JTrackSub = JTrackSubs::iterator;

DECLARE_SOA_TABLE(JTrackD0Subs, "AOD", "JTrackD0Subs",
                  o2::soa::Index<>,
                  bkgd0::CandidateId,
                  jtracksub::Pt,
                  jtracksub::Eta,
                  jtracksub::Phi,
                  jtracksub::Energy,
                  jtracksub::TrackSel,
                  jtracksub::Px<jtracksub::Pt, jtracksub::Phi>,
                  jtracksub::Py<jtracksub::Pt, jtracksub::Phi>,
                  jtracksub::Pz<jtracksub::Pt, jtracksub::Eta>,
                  jtracksub::P<jtracksub::Pt, jtracksub::Eta>);

using JTrackD0Sub = JTrackD0Subs::iterator;

DECLARE_SOA_TABLE(JTrackLcSubs, "AOD", "JTrackLcSubs",
                  o2::soa::Index<>,
                  bkglc::CandidateId,
                  jtracksub::Pt,
                  jtracksub::Eta,
                  jtracksub::Phi,
                  jtracksub::Energy,
                  jtracksub::TrackSel,
                  jtracksub::Px<jtracksub::Pt, jtracksub::Phi>,
                  jtracksub::Py<jtracksub::Pt, jtracksub::Phi>,
                  jtracksub::Pz<jtracksub::Pt, jtracksub::Eta>,
                  jtracksub::P<jtracksub::Pt, jtracksub::Eta>);

using JTrackLcSub = JTrackLcSubs::iterator;

DECLARE_SOA_TABLE(JTrackBplusSubs, "AOD", "JTrackBPlSubs",
                  o2::soa::Index<>,
                  bkgbplus::CandidateId,
                  jtracksub::Pt,
                  jtracksub::Eta,
                  jtracksub::Phi,
                  jtracksub::Energy,
                  jtracksub::TrackSel,
                  jtracksub::Px<jtracksub::Pt, jtracksub::Phi>,
                  jtracksub::Py<jtracksub::Pt, jtracksub::Phi>,
                  jtracksub::Pz<jtracksub::Pt, jtracksub::Eta>,
                  jtracksub::P<jtracksub::Pt, jtracksub::Eta>);

using JTrackBplusSub = JTrackBplusSubs::iterator;

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETSUBTRACTION_H_
