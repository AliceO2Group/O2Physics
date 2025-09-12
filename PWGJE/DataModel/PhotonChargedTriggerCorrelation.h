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
/// \brief table definitions for photon-hadron correlation analyses
///
/// \author Julius Kinner
/// \file PhotonChargedTriggerCorrelation.h

#ifndef PWGJE_DATAMODEL_PHOTONCHARGEDTRIGGERCORRELATION_H_
#define PWGJE_DATAMODEL_PHOTONCHARGEDTRIGGERCORRELATION_H_

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{

// basic correlation particle columns
namespace corr_particle
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetCollision, jetCollision, int, JCollisions, "");
DECLARE_SOA_INDEX_COLUMN_FULL(JetMcCollision, jetMcCollision, int, JMcCollisions, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
} // namespace corr_particle

// reco

// collision extension
namespace collision_extra_corr
{
DECLARE_SOA_COLUMN(SelEv, selEv, bool);
DECLARE_SOA_COLUMN(TrigEv, trigEv, bool);
DECLARE_SOA_COLUMN(NGlobalTracks, nGlobalTracks, int);
} // namespace collision_extra_corr
DECLARE_SOA_TABLE(CollisionsExtraCorr, "AOD", "COLLISIONSEXTRACORR",
                  collision_extra_corr::SelEv, collision_extra_corr::TrigEv, collision_extra_corr::NGlobalTracks);

// trigger
namespace trigger
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetTrack, jetTrack, int, JetTracks, "");
} // namespace trigger
DECLARE_SOA_TABLE(Triggers, "AOD", "TRIGGERS",
                  o2::soa::Index<>, corr_particle::JetCollisionId, trigger::JetTrackId,
                  corr_particle::Pt, corr_particle::Phi, corr_particle::Eta);
using Trigger = Triggers::iterator;

// hadrons (global tracks)
namespace hadron
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetTrack, jetTrack, int, JetTracks, "");
} // namespace hadron
DECLARE_SOA_TABLE(Hadrons, "AOD", "HADRONS",
                  o2::soa::Index<>, corr_particle::JetCollisionId, hadron::JetTrackId,
                  corr_particle::Pt, corr_particle::Phi, corr_particle::Eta);
using Hadron = Hadrons::iterator;

// pipm
namespace pipm
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetTrack, jetTrack, int, JetTracks, "");
} // namespace pipm
DECLARE_SOA_TABLE(Pipms, "AOD", "PIPMS",
                  o2::soa::Index<>, corr_particle::JetCollisionId, pipm::JetTrackId,
                  corr_particle::Pt, corr_particle::Phi, corr_particle::Eta);
using Pipm = Pipms::iterator;

// photonPCM
namespace photon_pcm
{
DECLARE_SOA_INDEX_COLUMN_FULL(V0PhotonKF, v0PhotonKF, int, V0PhotonsKF, "");
DECLARE_SOA_COLUMN(PosTrackId, posTrackId, int);
DECLARE_SOA_COLUMN(NegTrackId, negTrackId, int);
} // namespace photon_pcm
DECLARE_SOA_TABLE(PhotonPCMs, "AOD", "PHOTONPCMS",
                  o2::soa::Index<>, corr_particle::JetCollisionId, photon_pcm::V0PhotonKFId,
                  photon_pcm::PosTrackId, photon_pcm::NegTrackId,
                  corr_particle::Pt, corr_particle::Phi, corr_particle::Eta);
using PhotonPCM = PhotonPCMs::iterator;

// photonPCM pairs (pi0)
namespace photon_pcm_pair
{
DECLARE_SOA_INDEX_COLUMN_FULL(V0PhotonKF1, v0PhotonKF1, int, V0PhotonsKF, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(V0PhotonKF2, v0PhotonKF2, int, V0PhotonsKF, "_2");
DECLARE_SOA_COLUMN(PosTrack1Id, posTrack1Id, int);
DECLARE_SOA_COLUMN(NegTrack1Id, negTrack1Id, int);
DECLARE_SOA_COLUMN(PosTrack2Id, posTrack2Id, int);
DECLARE_SOA_COLUMN(NegTrack2Id, negTrack2Id, int);
DECLARE_SOA_COLUMN(Mgg, mgg, float);
} // namespace photon_pcm_pair
DECLARE_SOA_TABLE(PhotonPCMPairs, "AOD", "PHOTONPCMPAIRS",
                  o2::soa::Index<>, corr_particle::JetCollisionId, photon_pcm_pair::V0PhotonKF1Id, photon_pcm_pair::V0PhotonKF2Id,
                  photon_pcm_pair::PosTrack1Id, photon_pcm_pair::NegTrack1Id, photon_pcm_pair::PosTrack2Id, photon_pcm_pair::NegTrack2Id,
                  corr_particle::Pt, corr_particle::Phi, corr_particle::Eta, photon_pcm_pair::Mgg);
using PhotonPCMPair = PhotonPCMPairs::iterator;

// mc

// mcCollision extension
namespace mc_collision_extra_corr
{
DECLARE_SOA_COLUMN(TrigEv, trigEv, bool);
DECLARE_SOA_COLUMN(NChargedInEtaRange, nChargedInEtaRange, int);
} // namespace mc_collision_extra_corr
DECLARE_SOA_TABLE(McCollisionsExtraCorr, "AOD", "MCCOLLISIONSEXTRACORR",
                  mc_collision_extra_corr::TrigEv, mc_collision_extra_corr::NChargedInEtaRange);

// trigger
namespace trigger_particle
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetMcParticle, jetMcParticle, int, JetParticles, "");
} // namespace trigger_particle
DECLARE_SOA_TABLE(TriggerParticles, "AOD", "TRIGGERPARTICLES",
                  o2::soa::Index<>, corr_particle::JetMcCollisionId, trigger_particle::JetMcParticleId,
                  corr_particle::Pt, corr_particle::Phi, corr_particle::Eta);
using TriggerParticle = TriggerParticles::iterator;
} // namespace o2::aod

#endif // PWGJE_DATAMODEL_PHOTONCHARGEDTRIGGERCORRELATION_H_
