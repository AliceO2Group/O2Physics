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

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataDQ.h"

#include <Framework/ASoA.h>

#include <cmath>

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
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
} // namespace bkgcharged

namespace bkgd0
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfD0Bases, "_0");
} // namespace bkgd0

namespace bkgd0mc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfD0PBases, "_0");
} // namespace bkgd0mc

namespace bkgdplus
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfDplusBases, "_0");
} // namespace bkgdplus

namespace bkgdplusmc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfDplusPBases, "_0");
} // namespace bkgdplusmc

namespace bkgds
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfDsBases, "_0");
} // namespace bkgds

namespace bkgdsmc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfDsPBases, "_0");
} // namespace bkgdsmc

namespace bkgdstar
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfDstarBases, "_0");
} // namespace bkgdstar

namespace bkgdstarmc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfDstarPBases, "_0");
} // namespace bkgdstarmc

namespace bkglc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfLcBases, "_0");
} // namespace bkglc

namespace bkglcmc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfLcPBases, "_0");
} // namespace bkglcmc

namespace bkgb0
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfB0Bases, "_0");
} // namespace bkgb0

namespace bkgb0mc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfB0PBases, "_0");
} // namespace bkgb0mc

namespace bkgbplus
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfBplusBases, "_0");
} // namespace bkgbplus

namespace bkgbplusmc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfBplusPBases, "_0");
} // namespace bkgbplusmc

namespace bkgxictoxipipi
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfXicToXiPiPiBases, "_0");
} // namespace bkgxictoxipipi

namespace bkgxictoxipipimc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfXicToXiPiPiPBases, "_0");
} // namespace bkgxictoxipipimc

namespace bkgdielectron
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, Dielectrons, "_0");
} // namespace bkgdielectron

namespace bkgdielectronmc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, JDielectronMcs, "_0");
} // namespace bkgdielectronmc

DECLARE_SOA_TABLE(BkgChargedRhos, "AOD", "BkgCRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM);

DECLARE_SOA_TABLE(BkgChargedMcRhos, "AOD", "BkgCMcRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(BkgD0Rhos, "AOD", "BkgD0Rho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<2>);

DECLARE_SOA_TABLE(BkgD0McRhos, "AOD", "BkgD0McRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<3>);

DECLARE_SOA_TABLE(BkgDplusRhos, "AOD", "BkgDPRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<4>);

DECLARE_SOA_TABLE(BkgDplusMcRhos, "AOD", "BkgDPMcRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<5>);

DECLARE_SOA_TABLE(BkgDsRhos, "AOD", "BkgDSRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<6>);

DECLARE_SOA_TABLE(BkgDsMcRhos, "AOD", "BkgDSMcRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<7>);

DECLARE_SOA_TABLE(BkgDstarRhos, "AOD", "BkgDSTRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<8>);

DECLARE_SOA_TABLE(BkgDstarMcRhos, "AOD", "BkgDSTMcRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<9>);

DECLARE_SOA_TABLE(BkgLcRhos, "AOD", "BkgLCRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<10>);

DECLARE_SOA_TABLE(BkgLcMcRhos, "AOD", "BkgLCMcRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<11>);

DECLARE_SOA_TABLE(BkgB0Rhos, "AOD", "BkgB0Rho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<12>);

DECLARE_SOA_TABLE(BkgB0McRhos, "AOD", "BkgB0McRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<13>);

DECLARE_SOA_TABLE(BkgBplusRhos, "AOD", "BkgBPRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<14>);

DECLARE_SOA_TABLE(BkgBplusMcRhos, "AOD", "BkgBPMcRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<15>);

DECLARE_SOA_TABLE(BkgXicToXiPiPiRhos, "AOD", "BkgXICXPPRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<16>);

DECLARE_SOA_TABLE(BkgXicToXiPiPiMcRhos, "AOD", "BkgXICXPPMcRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<17>);

DECLARE_SOA_TABLE(BkgDielectronRhos, "AOD", "BkgDIELRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<18>);

DECLARE_SOA_TABLE(BkgDielectronMcRhos, "AOD", "BkgDIELMcRho",
                  o2::soa::Index<>,
                  bkgrho::Rho,
                  bkgrho::RhoM,
                  o2::soa::Marker<19>);

DECLARE_SOA_TABLE(JTrackSubs, "AOD", "JTrackSubs",
                  o2::soa::Index<>,
                  bkgcharged::JCollisionId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackSub = JTrackSubs::iterator;

DECLARE_SOA_TABLE(JMcParticleSubs, "AOD", "JMcPartSubs",
                  o2::soa::Index<>,
                  bkgcharged::JMcCollisionId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleSub = JMcParticleSubs::iterator;

DECLARE_SOA_TABLE(JTrackD0Subs, "AOD", "JTrackD0Subs",
                  o2::soa::Index<>,
                  bkgd0::CandidateId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackD0Sub = JTrackD0Subs::iterator;

DECLARE_SOA_TABLE(JMcParticleD0Subs, "AOD", "JMcPartD0Subs",
                  o2::soa::Index<>,
                  bkgd0mc::CandidateId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleD0Sub = JMcParticleD0Subs::iterator;

DECLARE_SOA_TABLE(JTrackDplusSubs, "AOD", "JTrackDPSubs",
                  o2::soa::Index<>,
                  bkgdplus::CandidateId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackDplusSub = JTrackDplusSubs::iterator;

DECLARE_SOA_TABLE(JMcParticleDplusSubs, "AOD", "JMcPartDPSubs",
                  o2::soa::Index<>,
                  bkgdplusmc::CandidateId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleDplusSub = JMcParticleDplusSubs::iterator;

DECLARE_SOA_TABLE(JTrackDsSubs, "AOD", "JTrackDSSubs",
                  o2::soa::Index<>,
                  bkgds::CandidateId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackDsSub = JTrackDsSubs::iterator;

DECLARE_SOA_TABLE(JMcParticleDsSubs, "AOD", "JMcPartDSSubs",
                  o2::soa::Index<>,
                  bkgdsmc::CandidateId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleDsSub = JMcParticleDsSubs::iterator;

DECLARE_SOA_TABLE(JTrackDstarSubs, "AOD", "JTrackDSTSubs",
                  o2::soa::Index<>,
                  bkgdstar::CandidateId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackDstarSub = JTrackDstarSubs::iterator;

DECLARE_SOA_TABLE(JMcParticleDstarSubs, "AOD", "JMcPartDSTSubs",
                  o2::soa::Index<>,
                  bkgdstarmc::CandidateId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleDstarSub = JMcParticleDstarSubs::iterator;

DECLARE_SOA_TABLE(JTrackLcSubs, "AOD", "JTrackLCSubs",
                  o2::soa::Index<>,
                  bkglc::CandidateId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackLcSub = JTrackLcSubs::iterator;

DECLARE_SOA_TABLE(JMcParticleLcSubs, "AOD", "JMcPartLCSubs",
                  o2::soa::Index<>,
                  bkglcmc::CandidateId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleLcSub = JMcParticleLcSubs::iterator;

DECLARE_SOA_TABLE(JTrackB0Subs, "AOD", "JTrackB0Subs",
                  o2::soa::Index<>,
                  bkgb0::CandidateId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackB0Sub = JTrackB0Subs::iterator;

DECLARE_SOA_TABLE(JMcParticleB0Subs, "AOD", "JMcPartB0Subs",
                  o2::soa::Index<>,
                  bkgb0mc::CandidateId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleB0Sub = JMcParticleB0Subs::iterator;

DECLARE_SOA_TABLE(JTrackBplusSubs, "AOD", "JTrackBPSubs",
                  o2::soa::Index<>,
                  bkgbplus::CandidateId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackBplusSub = JTrackBplusSubs::iterator;

DECLARE_SOA_TABLE(JMcParticleBplusSubs, "AOD", "JMcPartBPSubs",
                  o2::soa::Index<>,
                  bkgbplusmc::CandidateId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleBplusSub = JMcParticleBplusSubs::iterator;

DECLARE_SOA_TABLE(JTrackXicToXiPiPiSubs, "AOD", "JTrackXICXPPCSubs",
                  o2::soa::Index<>,
                  bkgxictoxipipi::CandidateId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackXicToXiPiPiSub = JTrackXicToXiPiPiSubs::iterator;

DECLARE_SOA_TABLE(JMcParticleXicToXiPiPiSubs, "AOD", "JMcPartXICXPPCSubs",
                  o2::soa::Index<>,
                  bkgxictoxipipimc::CandidateId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleXicToXiPiPiSub = JMcParticleXicToXiPiPiSubs::iterator;

DECLARE_SOA_TABLE(JTrackDielectronSubs, "AOD", "JTrackDIELSubs",
                  o2::soa::Index<>,
                  bkgdielectron::CandidateId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>);

using JTrackDielectronSub = JTrackDielectronSubs::iterator;

DECLARE_SOA_TABLE(JMcParticleDielectronSubs, "AOD", "JMcPartDIELSubs",
                  o2::soa::Index<>,
                  bkgdielectronmc::CandidateId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
                  jmcparticle::PdgCode,
                  jmcparticle::StatusCode,
                  jmcparticle::Flags,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::Energy<jmcparticle::E>,
                  jmcparticle::ProducedByGenerator<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::FromBackgroundEvent<jmcparticle::Flags>,                         // this will give a nonsensical value and should not be used
                  jmcparticle::GetProcess<jmcparticle::Flags, jmcparticle::StatusCode>,         // this will give a nonsensical value and should not be used
                  jmcparticle::GetGenStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>,   // this will give a nonsensical value and should not be used
                  jmcparticle::GetHepMCStatusCode<jmcparticle::Flags, jmcparticle::StatusCode>, // this will give a nonsensical value and should not be used
                  jmcparticle::IsPhysicalPrimary<jmcparticle::Flags>);

using JMcParticleDielectronSub = JMcParticleDielectronSubs::iterator;

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETSUBTRACTION_H_
