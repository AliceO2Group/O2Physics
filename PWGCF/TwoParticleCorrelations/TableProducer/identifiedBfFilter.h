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

#ifndef PWGCF_TWOPARTICLECORRELATIONS_TABLEPRODUCER_IDENTIFIEDBFFILTER_H_
#define PWGCF_TWOPARTICLECORRELATIONS_TABLEPRODUCER_IDENTIFIEDBFFILTER_H_

#include <vector>
#include <string>

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
using CollisionsEvSelCent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>;
using CollisionEvSelCent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>::iterator;
using CollisionsEvSelRun2Cent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s>;
using CollisionEvSelRun2Cent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s>::iterator;
using CollisionsEvSel = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
using CollisionEvSel = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>::iterator;
using TrackData = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>::iterator;
} // namespace aod
namespace analysis
{
namespace identifiedbffilter
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
  kppRun3,       ///< **p-p Run 3** system
  kPbPbRun3,     ///< **Pb-Pb Run 3** system
  knSystems      ///< number of handled systems
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
  kV0M,                ///< V0M centrality/multiplicity estimator Run 1/2
  kCL0,                ///< CL0 centrality/multiplicity estimator Run 1/2
  kCL1,                ///< CL1 centrality/multiplicity estimator Run 1/2
  kFV0A,               ///< FV0A centrality/multiplicity estimator Run 3
  kFT0M,               ///< FT0M centrality/multiplicity estimator Run 3
  kFT0A,               ///< FT0A centrality/multiplicity estimator Run 3
  kFT0C,               ///< FT0C centrality/multiplicity estimator Run 3
  kNTPV,               ///< NTPV centrality/multiplicity estimator Run 3
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
// The IdentifiedBfFilter configuration objects
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

std::vector<TrackSelection*> trackFilters = {};
bool dca2Dcut = false;
float maxDCAz = 1e6f;
float maxDCAxy = 1e6f;

inline void initializeTrackSelection()
{
  switch (tracktype) {
    case 1: { /* Run2 global track */
      TrackSelection* globalRun2 = new TrackSelection(getGlobalTrackSelection());
      globalRun2->SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
      globalRun2->SetMaxChi2PerClusterTPC(2.5f);
      TrackSelection* globalSDDRun2 = new TrackSelection(getGlobalTrackSelectionSDD());
      globalSDDRun2->SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
      globalSDDRun2->SetMaxChi2PerClusterTPC(2.5f);
      trackFilters.push_back(globalRun2);
      trackFilters.push_back(globalSDDRun2);
    } break;
    case 3: { /* Run3 track */
      TrackSelection* globalRun3 = new TrackSelection(getGlobalTrackSelection());
      globalRun3->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      globalRun3->ResetITSRequirements();
      globalRun3->SetRequireHitsInITSLayers(1, {0, 1, 2});
      TrackSelection* globalSDDRun3 = new TrackSelection(getGlobalTrackSelection());
      globalSDDRun3->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      globalSDDRun3->ResetITSRequirements();
      globalSDDRun3->SetRequireNoHitsInITSLayers({0, 1, 2});
      globalSDDRun3->SetRequireHitsInITSLayers(1, {3});
      trackFilters.push_back(globalRun3);
      trackFilters.push_back(globalSDDRun3);
    } break;
    case 5: { /* Run2 TPC only track */
      TrackSelection* tpcOnly = new TrackSelection;
      tpcOnly->SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
      tpcOnly->SetMinNClustersTPC(50);
      tpcOnly->SetMaxChi2PerClusterTPC(4);
      tpcOnly->SetMaxDcaZ(3.2f);
      maxDCAz = 3.2;
      tpcOnly->SetMaxDcaXY(2.4f);
      maxDCAxy = 2.4;
      dca2Dcut = true;
      trackFilters.push_back(tpcOnly);
    } break;
    case 7: { /* Run3 TPC only track */
      TrackSelection* tpcOnly = new TrackSelection;
      tpcOnly->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      tpcOnly->SetMinNClustersTPC(50);
      tpcOnly->SetMaxChi2PerClusterTPC(4);
      tpcOnly->SetMaxDcaZ(3.2f);
      maxDCAz = 3.2;
      tpcOnly->SetMaxDcaXY(2.4f);
      maxDCAxy = 2.4;
      dca2Dcut = true;
      trackFilters.push_back(tpcOnly);
    } break;
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
  if (triggstr.empty() || triggstr == "MB") {
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
  if (sysstr.empty() || (sysstr == "PbPb")) {
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
  } else if (sysstr == "PbPbRun3") {
    return kPbPbRun3;
  } else {
    LOGF(fatal, "IdentifiedBfCorrelations::getSystemType(). Wrong system type: %s", sysstr.c_str());
  }
  return kPbPb;
}

/// \brief Type of data according to the configuration string
/// \param datastr The data type configuration string
/// \return Internal code for the passed kind of data string
inline DataType getDataType(std::string const& datastr)
{
  /* we have to figure out how extract the type of data*/
  if (datastr.empty() || (datastr == "data")) {
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
    LOGF(fatal, "IdentifiedBfCorrelations::getDataType(). Wrong type of dat: %d", datastr.c_str());
  }
  return kData;
}

inline CentMultEstimatorType getCentMultEstimator(std::string const& datastr)
{
  if (datastr == "V0M") {
    return kV0M;
  } else if (datastr == "CL0") {
    return kCL0;
  } else if (datastr == "CL1") {
    return kCL1;
  } else if (datastr == "FV0A") {
    return kFV0A;
  } else if (datastr == "FT0M") {
    return kFT0M;
  } else if (datastr == "FT0A") {
    return kFT0A;
  } else if (datastr == "FT0C") {
    return kFT0C;
  } else if (datastr == "NTPV") {
    return kNTPV;
  } else if (datastr == "NOCM") {
    return kNOCM;
  } else {
    LOGF(fatal, "Centrality/Multiplicity estimator %s not supported yet", datastr.c_str());
  }
  return kNOCM;
}

inline std::string getCentMultEstimatorName(CentMultEstimatorType est)
{
  switch (est) {
    case kV0M:
      return "V0M";
      break;
    case kCL0:
      return "CL0";
      break;
    case kCL1:
      return "CL1";
      break;
    case kFV0A:
      return "FV0A";
      break;
    case kFT0M:
      return "FT0M";
      break;
    case kFT0A:
      return "FT0A";
      break;
    case kFT0C:
      return "FT0C";
      break;
    case kNTPV:
      return "NTPV";
      break;
    case kNOCM:
      return "NOCM";
      break;
    default:
      LOGF(fatal, "Centrality/Multiplicity estimator %d not supported yet", (int)est);
      return "WRONG";
      break;
  }
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
              if (collision.alias_bit(kINT7)) {
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
    case kPbPbRun3:
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

/// \brief Trigger selection for reconstructed collision tables with Run 2 centrality/multiplicity
template <>
inline bool triggerSelection<aod::CollisionEvSelRun2Cent>(aod::CollisionEvSelRun2Cent const& collision)
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

/// \brief Trigger selection for detector level collision tables with Run 2 centrality/multiplicity
template <>
inline bool triggerSelection<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator const& collision)
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
inline float extractMultiplicity(CollisionObject const& collision, CentMultEstimatorType est)
{
  switch (est) {
    case kV0M:
      return collision.multFV0M();
      break;
    case kCL0:
      return collision.multTracklets();
      break;
    case kCL1:
      return collision.multTracklets();
      break;
    case kFV0A:
      return collision.multFV0A();
      break;
    case kFT0M:
      return collision.multFT0M();
      break;
    case kFT0A:
      return collision.multFT0A();
      break;
    case kFT0C:
      return collision.multFT0M();
      break;
    case kNTPV:
      return collision.multNTracksPV();
      break;
    case kNOCM:
      return collision.multFT0M();
      break;
    default:
      LOGF(fatal, "Centrality/Multiplicity estimator %d not supported yet", (int)est);
      return collision.multFT0M();
      break;
  }
}

//////////////////////////////////////////////////////////////////////////////////
/// Centrality selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief Centrality/multiplicity percentile
template <typename CollisionObject>
float getCentMultPercentile(CollisionObject collision)
{
  if constexpr (framework::has_type_v<aod::cent::CentRun2V0M, typename CollisionObject::all_columns> ||
                framework::has_type_v<aod::cent::CentRun2CL0, typename CollisionObject::all_columns> ||
                framework::has_type_v<aod::cent::CentRun2CL1, typename CollisionObject::all_columns>) {
    switch (fCentMultEstimator) {
      case kV0M:
        return collision.centRun2V0M();
        break;
      case kCL0:
        return collision.centRun2CL0();
        break;
      case kCL1:
        return collision.centRun2CL1();
        break;
      default:
        return 105.0;
        break;
    }
  }
  if constexpr (framework::has_type_v<aod::cent::CentFV0A, typename CollisionObject::all_columns> ||
                framework::has_type_v<aod::cent::CentFT0M, typename CollisionObject::all_columns> ||
                framework::has_type_v<aod::cent::CentFT0A, typename CollisionObject::all_columns> ||
                framework::has_type_v<aod::cent::CentFT0C, typename CollisionObject::all_columns> ||
                framework::has_type_v<aod::cent::CentNTPV, typename CollisionObject::all_columns>) {
    switch (fCentMultEstimator) {
      case kFV0A:
        return collision.centFV0A();
        break;
      case kFT0M:
        return collision.centFT0M();
        break;
      case kFT0A:
        return collision.centFT0A();
        break;
      case kFT0C:
        return collision.centFT0C();
        break;
      case kNTPV:
        return collision.centNTPV();
        break;
      default:
        return 105.0;
        break;
    }
  }
}

/// \brief Centrality selection when there is centrality/multiplicity information
template <typename CollisionObject>
inline bool centralitySelectionMult(CollisionObject collision, float& centmult)
{
  float mult = getCentMultPercentile(collision);
  if (mult < 100 && 0 < mult) {
    centmult = mult;
    return true;
  }
  return false;
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

/// \brief Centrality selection for reconstructed and detector level collision tables with Run 2 centrality/multiplicity information
template <>
inline bool centralitySelection<aod::CollisionEvSelRun2Cent>(aod::CollisionEvSelRun2Cent const& collision, float& centmult)
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

/// \brief Centrality selection for detector level collision tables with Run 2 centrality/multiplicity
template <>
inline bool centralitySelection<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for generator level collision table
template <>
inline bool centralitySelection<aod::McCollision>(aod::McCollision const&, float& centmult)
{
  if (centmult < 100 && 0 < centmult) {
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
  if (zvtxlow < collision.posZ() && collision.posZ() < zvtxup) {
    zvtxsel = true;
  }

  bool centmultsel = centralitySelection(collision, centormult);

  return trigsel && zvtxsel && centmultsel;
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
    for (auto filter : trackFilters) {
      if (filter->IsSelected(track)) {
        if (dca2Dcut) {
          if (track.dcaXY() * track.dcaXY() / maxDCAxy / maxDCAxy + track.dcaZ() * track.dcaZ() / maxDCAz / maxDCAz > 1) {
            return false;
          } else {
            return true;
          }
        } else {
          return true;
        }
      }
    }
    return false;
  }
}

/// \brief Accepts or not the passed track
/// \param track the track of interest
/// \return the internal track id, -1 if not accepted
/// TODO: the PID implementation
/// For the time being we keep the convention
/// - positive track pid even
/// - negative track pid odd
/// - charged hadron 0/1
template <typename TrackObject>
inline int8_t AcceptTrack(TrackObject const& track)
{
  /* TODO: incorporate a mask in the scanned tracks table for the rejecting track reason */
  if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
    if (track.mcParticleId() < 0) {
      return -1;
    }
  }

  if (matchTrackType(track)) {
    if (ptlow < track.pt() && track.pt() < ptup && etalow < track.eta() && track.eta() < etaup) {
      if (track.sign() > 0) {
        return 0;
      }
      if (track.sign() < 0) {
        return 1;
      }
    }
  }
  return -1;
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

/// \brief Accepts or not the passed generated particle
/// \param track the particle of interest
/// \return the internal particle id, -1 if not accepted
/// TODO: the PID implementation
/// For the time being we keep the convention
/// - positive particle pid even
/// - negative particle pid odd
/// - charged hadron 0/1
template <typename ParticleObject, typename MCCollisionObject>
inline int8_t AcceptParticle(ParticleObject& particle, MCCollisionObject const& collision)
{
  float charge = (fPDG->GetParticle(particle.pdgCode())->Charge() / 3 >= 1) ? 1.0 : ((fPDG->GetParticle(particle.pdgCode())->Charge() / 3 <= -1) ? -1.0 : 0.0);

  if (particle.isPhysicalPrimary()) {
    if ((particle.mcCollisionId() == 0) && traceCollId0) {
      LOGF(info, "Particle %d passed isPhysicalPrimary", particle.globalIndex());
    }
    if (useOwnParticleSelection) {
      float dcaxy = TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
      float dcaz = TMath::Abs(particle.vz() - collision.posZ());
      if (!((dcaxy < particleMaxDCAxy) && (dcaz < particleMaxDCAZ))) {
        if ((particle.mcCollisionId() == 0) && traceCollId0) {
          LOGF(info, "Rejecting particle with dcaxy: %.2f and dcaz: %.2f", dcaxy, dcaz);
          LOGF(info, "   assigned collision Id: %d, looping on collision Id: %d", particle.mcCollisionId(), collision.globalIndex());
          LOGF(info, "   Collision x: %.5f, y: %.5f, z: %.5f", collision.posX(), collision.posY(), collision.posZ());
          LOGF(info, "   Particle x: %.5f, y: %.5f, z: %.5f", particle.vx(), particle.vy(), particle.vz());
          LOGF(info, "   index: %d, pdg code: %d", particle.globalIndex(), particle.pdgCode());

          exploreMothers(particle, collision);
        }
        return -1;
      }
    }
    if (ptlow < particle.pt() && particle.pt() < ptup && etalow < particle.eta() && particle.eta() < etaup) {
      if (charge > 0) {
        return 0;
      }
      if (charge < 0) {
        return 1;
      }
    }
  } else {
    if ((particle.mcCollisionId() == 0) && traceCollId0) {
      LOGF(info, "Particle %d NOT passed isPhysicalPrimary", particle.globalIndex());
    }
  }
  return -1;
}

} // namespace identifiedbffilter
} // namespace analysis
} // namespace o2

#endif // PWGCF_TWOPARTICLECORRELATIONS_TABLEPRODUCER_IDENTIFIEDBFFILTER_H_
