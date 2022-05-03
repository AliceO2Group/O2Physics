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
#ifndef O2_ANALYSIS_DPTDPTFILTER_H
#define O2_ANALYSIS_DPTDPTFILTER_H

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include <TDatabasePDG.h>

namespace o2
{
namespace aod
{
using CollisionsEvSelCent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentRun2V0Ms>;
using CollisionEvSelCent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentRun2V0Ms>::iterator;
using CollisionsEvSel = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
using CollisionEvSel = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>::iterator;
using TrackData = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection>::iterator;
} // namespace aod
namespace analysis
{
namespace dptdptfilter
{
/// \enum SystemType
/// \brief The type of the system under analysis
enum SystemType {
  kNoSystem = 0, ///< no system defined
  kpp,           ///< **p-p** system
  kpPb,          ///< **p-Pb** system
  kPbp,          ///< **Pb-p** system
  kPbPb,         ///< **Pb-Pb** system
  kXeXe,         ///< **Xe-Xe** system
  kppRun3,
  knSystems ///< number of handled systems
};

/// \enum DataType
/// \brief Which kind of data is the task addressing
enum DataType {
  kData = 0,     ///< actual data, not generated
  kMC,           ///< Generator level and detector level
  kFastMC,       ///< Gererator level but stored dataset
  kOnTheFly,     ///< On the fly generator level data
  kDataNoEvtSel, ///< actual data but not event selection available yet
  knGenData      ///< number of different generator data types
};

/// \enum CentMultEstimatorType
/// \brief The detector used to estimate centrality/multiplicity
enum CentMultEstimatorType {
  kNOCM = 0,           ///< do not use centrality/multiplicity estimator
  kV0M,                ///< V0M centrality/multiplicity estimator
  kV0A,                ///< V0A centrality/multiplicity estimator
  kV0C,                ///< V0C centrality/multiplicity estimator
  kCL0,                ///< CL0 centrality/multiplicity estimator
  kCL1,                ///< CL1 centrality/multiplicity estimator
  knCentMultEstimators ///< number of centrality/mutiplicity estimator
};

/// \enum TriggerSelectionType
/// \brief The type of trigger to apply for event selection
enum TriggerSelectionType {
  kNONE = 0,       ///< do not use trigger selection
  kMB,             ///< Minimum bias trigger
  knEventSelection ///< number of triggers for event selection
};

//============================================================================================
// The DptDptFilter configuration objects
//============================================================================================
int ptbins = 18;
float ptlow = 0.2, ptup = 2.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int zvtxbins = 40;
float zvtxlow = -10.0, zvtxup = 10.0;
int phibins = 72;
float philow = 0.0;
float phiup = constants::math::TwoPI;

int tracktype = 1;
int trackonecharge = 1;
int tracktwocharge = -1;

TrackSelection* globalRun3 = nullptr;
TrackSelection* globalSDDRun3 = nullptr;

inline void initializeTrackSelection()
{
  switch (tracktype) {
    case 3: /* Run3 track */
      globalRun3 = new TrackSelection(getGlobalTrackSelection());
      globalRun3->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      globalRun3->ResetITSRequirements();
      globalRun3->SetRequireHitsInITSLayers(1, {0, 1, 2});
      globalSDDRun3 = new TrackSelection(getGlobalTrackSelection());
      globalSDDRun3->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      globalSDDRun3->ResetITSRequirements();
      globalSDDRun3->SetRequireNoHitsInITSLayers({0, 1, 2});
      globalSDDRun3->SetRequireHitsInITSLayers(1, {3});
      break;
    default:
      break;
  }
}

SystemType fSystem = kNoSystem;
DataType fDataType = kData;
CentMultEstimatorType fCentMultEstimator = kV0M;
TriggerSelectionType fTriggerSelection = kMB;

/* adaptations for the pp nightly checks */
analysis::CheckRangeCfg traceDCAOutliers;
bool traceOutOfSpeciesParticles = false;
int recoIdMethod = 0;
bool useOwnTrackSelection = false;
TrackSelection ownTrackSelection = getGlobalTrackSelection();
bool useOwnParticleSelection = false;
float particleMaxDCAxy = 999.9f;
float particleMaxDCAZ = 999.9f;
bool traceCollId0 = false;

TDatabasePDG* fPDG = nullptr;

inline TriggerSelectionType getTriggerSelection(std::string const& triggstr)
{
  if (triggstr.empty() or triggstr == "MB") {
    return kMB;
  } else if (triggstr == "None") {
    return kNONE;
  } else {
    LOGF(fatal, "Wrong trigger selection: %s", triggstr.c_str());
    return kNONE;
  }
}

inline SystemType getSystemType(std::string const& sysstr)
{
  /* we have to figure out how extract the system type */
  if (sysstr.empty() or (sysstr == "PbPb")) {
    return kPbPb;
  } else if (sysstr == "pp") {
    return kpp;
  } else if (sysstr == "pPb") {
    return kpPb;
  } else if (sysstr == "Pbp") {
    return kPbp;
  } else if (sysstr == "pPb") {
    return kpPb;
  } else if (sysstr == "XeXe") {
    return kXeXe;
  } else if (sysstr == "ppRun3") {
    return kppRun3;
  } else {
    LOGF(fatal, "DptDptCorrelations::getSystemType(). Wrong system type: %s", sysstr.c_str());
  }
  return kPbPb;
}

/// \brief Type of data according to the configuration string
/// \param datastr The data type configuration string
/// \return Internal code for the passed kind of data string
inline DataType getDataType(std::string const& datastr)
{
  /* we have to figure out how extract the type of data*/
  if (datastr.empty() or (datastr == "data")) {
    return kData;
  } else if (datastr == "datanoevsel") {
    return kDataNoEvtSel;
  } else if (datastr == "MC") {
    return kMC;
  } else if (datastr == "FastMC") {
    return kFastMC;
  } else if (datastr == "OnTheFlyMC") {
    return kOnTheFly;
  } else {
    LOGF(fatal, "DptDptCorrelations::getDataType(). Wrong type of dat: %d", datastr.c_str());
  }
  return kData;
}

inline CentMultEstimatorType getCentMultEstimator(std::string const& datastr)
{
  if (datastr == "V0M") {
    return kV0M;
  } else if (datastr == "NOCM") {
    return kNOCM;
  } else {
    LOGF(fatal, "Centrality/Multiplicity estimator %s not supported yet", datastr.c_str());
  }
  return kNOCM;
}

//////////////////////////////////////////////////////////////////////////////////
/// Trigger selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief Trigger selection for reconstructed and detector level collision tables
template <typename CollisionObject>
inline bool triggerSelectionReco(CollisionObject const& collision)
{
  bool trigsel = false;
  switch (fSystem) {
    case kpp:
    case kpPb:
    case kPbp:
    case kPbPb:
    case kXeXe:
      switch (fTriggerSelection) {
        case kMB:
          switch (fDataType) {
            case kData:
              if (collision.alias()[kINT7]) {
                if (collision.sel7()) {
                  trigsel = true;
                }
              }
              break;
            case kMC:
              if (collision.sel7()) {
                trigsel = true;
              }
              break;
            default:
              trigsel = true;
              break;
          }
          break;
        case kNONE:
          trigsel = true;
          break;
        default:
          break;
      }
      break;
    case kppRun3:
      switch (fTriggerSelection) {
        case kMB:
          if (collision.sel8()) {
            trigsel = true;
          }
          break;
        case kNONE:
          trigsel = true;
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
  return trigsel;
}

/// \brief Trigger selection by default: unknow subscribed collision table
template <typename CollisionObject>
inline bool triggerSelection(CollisionObject const&)
{
  LOGF(fatal, "Trigger selection not implemented for this kind of collisions");
  return false;
}

/// \brief Trigger selection for reconstructed collision tables without centrality/multiplicity
template <>
inline bool triggerSelection<aod::CollisionEvSel>(aod::CollisionEvSel const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for reconstructed collision tables with centrality/multiplicity
template <>
inline bool triggerSelection<aod::CollisionEvSelCent>(aod::CollisionEvSelCent const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for detector level collision tables without centrality/multiplicity
template <>
inline bool triggerSelection<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for detector level collision tables with centrality/multiplicity
template <>
inline bool triggerSelection<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for generator level collison table
template <>
inline bool triggerSelection<aod::McCollision>(aod::McCollision const&)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////////
/// Multiplicity extraction
//////////////////////////////////////////////////////////////////////////////////

/// \brief Extract the collision multiplicity from the event selection information
template <typename CollisionObject>
inline float extractMultiplicity(CollisionObject const& collision)
{
  float mult = 0.0;
  switch (fSystem) {
    case kpp:
    case kpPb:
    case kPbp:
    case kPbPb:
    case kXeXe:
      /* for the time being let's extract V0M */
      mult = collision.multFV0M();
      break;
    case kppRun3:
      /* for the time being let's extract T0M */
      mult = collision.multFT0M();
      break;
    default:
      break;
  }
  return mult;
}

//////////////////////////////////////////////////////////////////////////////////
/// Centrality selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief Centrality selection when there is centrality/multiplicity information
template <typename CollisionObject>
inline bool centralitySelectionMult(CollisionObject collision, float& centmult)
{
  bool centmultsel = false;
  switch (fCentMultEstimator) {
    case kV0M:
      if (collision.centRun2V0M() < 100 and 0 < collision.centRun2V0M()) {
        centmult = collision.centRun2V0M();
        centmultsel = true;
      }
      break;
    default:
      break;
  }
  return centmultsel;
}

/// \brief Centrality selection when there is not centrality/multiplicity information
template <typename CollisionObject>
inline bool centralitySelectionNoMult(CollisionObject const&, float& centmult)
{
  bool centmultsel = false;
  switch (fCentMultEstimator) {
    case kNOCM:
      centmult = 50.0;
      centmultsel = true;
      break;
    default:
      break;
  }
  return centmultsel;
}

/// \brief Centrality selection by default: unknown subscribed collision table
template <typename CollisionObject>
inline bool centralitySelection(CollisionObject const&, float&)
{
  LOGF(fatal, "Centrality selection not implemented for this kind of collisions");
  return false;
}

/// \brief Centrality selection for reconstructed and detector level collision tables with centrality/multiplicity information
template <>
inline bool centralitySelection<aod::CollisionEvSelCent>(aod::CollisionEvSelCent const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for reconstructed and detector level collision tables without centrality/multiplicity information
template <>
inline bool centralitySelection<aod::CollisionEvSel>(aod::CollisionEvSel const& collision, float& centmult)
{
  return centralitySelectionNoMult(collision, centmult);
}

/// \brief Centrality selection for detector level collision tables without centrality/multiplicity
template <>
inline bool centralitySelection<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator const& collision, float& centmult)
{
  return centralitySelectionNoMult(collision, centmult);
}

/// \brief Centrality selection for detector level collision tables with centrality/multiplicity
template <>
inline bool centralitySelection<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for generator level collision table
template <>
inline bool centralitySelection<aod::McCollision>(aod::McCollision const&, float& centmult)
{
  if (centmult < 100 and 0 < centmult) {
    return true;
  } else {
    return false;
  }
}

//////////////////////////////////////////////////////////////////////////////////
/// Event selection
//////////////////////////////////////////////////////////////////////////////////

template <typename CollisionObject>
inline bool IsEvtSelected(CollisionObject const& collision, float& centormult)
{
  bool trigsel = triggerSelection(collision);

  bool zvtxsel = false;
  /* TODO: vertex quality checks */
  if (zvtxlow < collision.posZ() and collision.posZ() < zvtxup) {
    zvtxsel = true;
  }

  bool centmultsel = centralitySelection(collision, centormult);

  return trigsel and zvtxsel and centmultsel;
}

//////////////////////////////////////////////////////////////////////////////////
/// Track selection
//////////////////////////////////////////////////////////////////////////////////

template <typename TrackObject>
inline bool matchTrackType(TrackObject const& track)
{
  if (useOwnTrackSelection) {
    return ownTrackSelection.IsSelected(track);
  } else {
    switch (tracktype) {
      case 1:
        if (track.isGlobalTrack() != 0 || track.isGlobalTrackSDD() != 0) {
          return true;
        } else {
          return false;
        }
        break;
      case 3: /* Run3 track */
        if (globalRun3->IsSelected(track) || globalSDDRun3->IsSelected(track)) {
          return true;
        } else {
          return false;
        }
        break;
      default:
        return false;
    }
  }
}

template <typename TrackObject>
inline void AcceptTrack(TrackObject const& track, uint8_t& asone, uint8_t& astwo)
{
  asone = uint8_t(false);
  astwo = uint8_t(false);

  /* TODO: incorporate a mask in the scanned tracks table for the rejecting track reason */
  if (matchTrackType(track)) {
    if (ptlow < track.pt() and track.pt() < ptup and etalow < track.eta() and track.eta() < etaup) {
      if (((track.sign() > 0) and (trackonecharge > 0)) or ((track.sign() < 0) and (trackonecharge < 0))) {
        asone = uint8_t(true);
      }
      if (((track.sign() > 0) and (tracktwocharge > 0)) or ((track.sign() < 0) and (tracktwocharge < 0))) {
        astwo = uint8_t(true);
      }
    }
  }
}

template <typename ParticleObject, typename MCCollisionObject>
void exploreMothers(ParticleObject& particle, MCCollisionObject& collision)
{
  for (auto& m : particle.template mothers_as<aod::McParticles>()) {
    LOGF(info, "   mother index: %d", m.globalIndex());
    LOGF(info, "   Tracking back mother");
    LOGF(info, "   assigned collision Id: %d, looping on collision Id: %d", m.mcCollisionId(), collision.globalIndex());
    LOGF(info, "   index: %d, pdg code: %d", m.globalIndex(), m.pdgCode());
    LOGF(info, "   Passed  isPhysicalPrimary(): %s", m.isPhysicalPrimary() ? "YES" : "NO");

    exploreMothers(m, collision);
  }
}

template <typename ParticleObject, typename MCCollisionObject>
inline void AcceptParticle(ParticleObject& particle, MCCollisionObject const& collision, uint8_t& asone, uint8_t& astwo)
{
  asone = uint8_t(false);
  astwo = uint8_t(false);

  float charge = (fPDG->GetParticle(particle.pdgCode())->Charge() / 3 >= 1) ? 1.0 : ((fPDG->GetParticle(particle.pdgCode())->Charge() / 3 <= -1) ? -1.0 : 0.0);

  if (particle.isPhysicalPrimary()) {
    if ((particle.mcCollisionId() == 0) and traceCollId0) {
      LOGF(info, "Particle %d passed isPhysicalPrimary", particle.globalIndex());
    }
    if (useOwnParticleSelection) {
      float dcaxy = TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
      float dcaz = TMath::Abs(particle.vz() - collision.posZ());
      if (not((dcaxy < particleMaxDCAxy) and (dcaz < particleMaxDCAZ))) {
        if ((particle.mcCollisionId() == 0) and traceCollId0) {
          LOGF(info, "Rejecting particle with dcaxy: %.2f and dcaz: %.2f", dcaxy, dcaz);
          LOGF(info, "   assigned collision Id: %d, looping on collision Id: %d", particle.mcCollisionId(), collision.globalIndex());
          LOGF(info, "   Collision x: %.5f, y: %.5f, z: %.5f", collision.posX(), collision.posY(), collision.posZ());
          LOGF(info, "   Particle x: %.5f, y: %.5f, z: %.5f", particle.vx(), particle.vy(), particle.vz());
          LOGF(info, "   index: %d, pdg code: %d", particle.globalIndex(), particle.pdgCode());

          exploreMothers(particle, collision);
        }
        return;
      }
    }
    if (ptlow < particle.pt() and particle.pt() < ptup and etalow < particle.eta() and particle.eta() < etaup) {
      if (((charge > 0) and (trackonecharge > 0)) or ((charge < 0) and (trackonecharge < 0))) {
        asone = uint8_t(true);
      }
      if (((charge > 0) and (tracktwocharge > 0)) or ((charge < 0) and (tracktwocharge < 0))) {
        astwo = uint8_t(true);
      }
    }
  } else {
    if ((particle.mcCollisionId() == 0) and traceCollId0) {
      LOGF(info, "Particle %d NOT passed isPhysicalPrimary", particle.globalIndex());
    }
  }
}

} // namespace dptdptfilter
} // namespace analysis
} // namespace o2

#endif // O2_ANALYSIS_DPTDPTFILTER_H
