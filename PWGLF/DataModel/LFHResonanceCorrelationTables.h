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

/// \brief This task serves to do hadron-resonance correlation studies.
///  The yield will be calculated using the two-particle correlation method.
///  Trigger particle : Hadrons
///  Associated Particles : Phi, K*0
///
/// \author Hirak Kumar Koley (hirak.koley@cern.ch)

#ifndef PWGLF_DATAMODEL_LFHRESONANCECORRELATIONTABLES_H_
#define PWGLF_DATAMODEL_LFHRESONANCECORRELATIONTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

// Simple checker
#define bitcheck(var, nbit) ((var) & (1 << (nbit)))

namespace o2::aod
{
/// _________________________________________
/// Table for storing trigger track indices
namespace triggerTracks
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                       //!
DECLARE_SOA_COLUMN(MCPhysicalPrimary, mcPhysicalPrimary, bool);       // true physical primary flag
DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_Trigger"); //!
DECLARE_SOA_COLUMN(MCOriginalPt, mcOriginalPt, float);                // true generated pt
DECLARE_SOA_COLUMN(IsLeading, isLeading, bool);                       // is leading track in the event
} // namespace triggerTracks
DECLARE_SOA_TABLE(TriggerTracks, "AOD", "TRIGGERTRACKS", o2::soa::Index<>, triggerTracks::CollisionId, triggerTracks::MCPhysicalPrimary, triggerTracks::TrackId, triggerTracks::MCOriginalPt, triggerTracks::IsLeading);
namespace triggerTrackExtras
{
DECLARE_SOA_COLUMN(Extra, extra, int); // true physical primary flag
} // namespace triggerTrackExtras
DECLARE_SOA_TABLE(TriggerTrackExtras, "AOD", "TRIGGERTRACKEXTRAs", triggerTrackExtras::Extra);

/// _________________________________________
/// Table for storing assoc track indices
namespace assocHadrons
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                     //!
DECLARE_SOA_COLUMN(MCPhysicalPrimary, mcPhysicalPrimary, bool);     // true physical primary flag
DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_Assoc"); //!
DECLARE_SOA_COLUMN(MCOriginalPt, mcOriginalPt, float);              // true generated pt
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);                          // pdg code of the MC particle
} // namespace assocHadrons
DECLARE_SOA_TABLE(AssocHadrons, "AOD", "ASSOCHADRONS", o2::soa::Index<>, assocHadrons::CollisionId, assocHadrons::MCPhysicalPrimary, assocHadrons::TrackId, assocHadrons::MCOriginalPt, assocHadrons::PDGCode);
/// _________________________________________
/// Table for storing assoc track PID
namespace assocPID
{
DECLARE_SOA_COLUMN(NSigmaTPCPi, nSigmaTPCPi, float);
DECLARE_SOA_COLUMN(NSigmaTPCKa, nSigmaTPCKa, float);
DECLARE_SOA_COLUMN(NSigmaTPCPr, nSigmaTPCPr, float);
DECLARE_SOA_COLUMN(NSigmaTPCEl, nSigmaTPCEl, float);
DECLARE_SOA_COLUMN(NSigmaTOFPi, nSigmaTOFPi, float);
DECLARE_SOA_COLUMN(NSigmaTOFKa, nSigmaTOFKa, float);
DECLARE_SOA_COLUMN(NSigmaTOFPr, nSigmaTOFPr, float);
DECLARE_SOA_COLUMN(NSigmaTOFEl, nSigmaTOFEl, float);
} // namespace assocPID
DECLARE_SOA_TABLE(AssocPID, "AOD", "ASSOCPID", assocPID::NSigmaTPCPi, assocPID::NSigmaTPCKa, assocPID::NSigmaTPCPr, assocPID::NSigmaTPCEl, assocPID::NSigmaTOFPi, assocPID::NSigmaTOFKa, assocPID::NSigmaTOFPr, assocPID::NSigmaTOFEl);

/// _________________________________________
/// Table for storing associated phi candidate indices
namespace assocPhis
{

DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!

// MC information
DECLARE_SOA_COLUMN(MCTruePhi, mcTruePhi, bool);
DECLARE_SOA_COLUMN(MCPhysicalPrimary, mcPhysicalPrimary, bool);

// Quality variable
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Mass, mass, float);

DECLARE_SOA_INDEX_COLUMN_FULL(PosTrackDaughter, posTrackDaughter, int, Tracks, "_PhiPosDaughter");
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrackDaughter, negTrackDaughter, int, Tracks, "_PhiNegDaughter");

} // namespace assocPhis

DECLARE_SOA_TABLE(
  AssocPhis,
  "AOD",
  "ASSOCPHIS",
  o2::soa::Index<>,
  assocPhis::CollisionId,
  assocPhis::MCTruePhi,
  assocPhis::MCPhysicalPrimary,
  assocPhis::Pt,
  assocPhis::Eta,
  assocPhis::Phi,
  assocPhis::Mass,
  assocPhis::PosTrackDaughterId,
  assocPhis::NegTrackDaughterId);

/// _________________________________________
/// Table for storing associated K*0 candidate indices
/// K*0 -> K+ pi- and anti-K*0 -> K- pi+ are stored together (sign-agnostic,
/// same convention as AssocPhis): PosTrackDaughter/NegTrackDaughter identify
/// the positively/negatively charged daughter, irrespective of which one is
/// the kaon and which one is the pion for that candidate.
namespace assocKstars
{

DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!

// MC information
DECLARE_SOA_COLUMN(MCTrueKstar, mcTrueKstar, bool);
DECLARE_SOA_COLUMN(MCPhysicalPrimary, mcPhysicalPrimary, bool);

// Quality variable
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Mass, mass, float);

DECLARE_SOA_INDEX_COLUMN_FULL(PosTrackDaughter, posTrackDaughter, int, Tracks, "_KstarPosDaughter");
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrackDaughter, negTrackDaughter, int, Tracks, "_KstarNegDaughter");

} // namespace assocKstars

DECLARE_SOA_TABLE(
  AssocKstars,
  "AOD",
  "ASSOCKSTARS",
  o2::soa::Index<>,
  assocKstars::CollisionId,
  assocKstars::MCTrueKstar,
  assocKstars::MCPhysicalPrimary,
  assocKstars::Pt,
  assocKstars::Eta,
  assocKstars::Phi,
  assocKstars::Mass,
  assocKstars::PosTrackDaughterId,
  assocKstars::NegTrackDaughterId);
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFHRESONANCECORRELATIONTABLES_H_
