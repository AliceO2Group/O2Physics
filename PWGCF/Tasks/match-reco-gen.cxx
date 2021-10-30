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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/MC.h"
#include "Common/Core/PID/PIDResponse.h"
#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TList.h>
#include <TDirectory.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile3D.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;

#define MATCHRECGENLOGCOLLISIONS debug
#define MATCHRECGENLOGTRACKS debug

#include "Framework/runDataProcessing.h"

namespace o2
{
namespace aod
{
/* we have to change from int to bool when bool columns work properly */
namespace dptdptcorrelations
{
DECLARE_SOA_COLUMN(EventAccepted, eventaccepted, uint8_t); //! If the collision/event has been accepted or not
DECLARE_SOA_COLUMN(EventCentMult, centmult, float);        //! The centrality/multiplicity pecentile
} // namespace dptdptcorrelations
DECLARE_SOA_TABLE(AcceptedEvents, "AOD", "ACCEPTEDEVENTS", //! Accepted reconstructed collisions/events filtered table
                  o2::soa::Index<>,
                  collision::BCId,
                  collision::PosZ,
                  dptdptcorrelations::EventAccepted,
                  dptdptcorrelations::EventCentMult);
using AcceptedEvent = AcceptedEvents::iterator;
DECLARE_SOA_TABLE(AcceptedTrueEvents, "AOD", "ACCTRUEEVENTS", //! Accepted generated collisions/events filtered table
                  o2::soa::Index<>,
                  collision::BCId,
                  mccollision::PosZ,
                  dptdptcorrelations::EventAccepted,
                  dptdptcorrelations::EventCentMult);
using AcceptedTrueEvent = AcceptedTrueEvents::iterator;
namespace dptdptcorrelations
{
DECLARE_SOA_INDEX_COLUMN(AcceptedEvent, event);                      //! Reconstructed collision/event
DECLARE_SOA_INDEX_COLUMN(AcceptedTrueEvent, mcevent);                //! Generated collision/event
DECLARE_SOA_COLUMN(TrackacceptedAsOne, trackacceptedasone, uint8_t); //! Track accepted as type one
DECLARE_SOA_COLUMN(TrackacceptedAsTwo, trackacceptedastwo, uint8_t); //! Track accepted as type two
DECLARE_SOA_COLUMN(Pt, pt, float);                                   //! The track transverse momentum
DECLARE_SOA_COLUMN(Eta, eta, float);                                 //! The track pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                                 //! The track azimuthal angle
} // namespace dptdptcorrelations
DECLARE_SOA_TABLE(ScannedTracks, "AOD", "SCANNEDTRACKS", //! The reconstructed tracks filtered table
                  dptdptcorrelations::AcceptedEventId,
                  dptdptcorrelations::TrackacceptedAsOne,
                  dptdptcorrelations::TrackacceptedAsTwo,
                  dptdptcorrelations::Pt,
                  dptdptcorrelations::Eta,
                  dptdptcorrelations::Phi);
DECLARE_SOA_TABLE(ScannedTrueTracks, "AOD", "SCANTRUETRACKS", //! The generated particles filtered table
                  dptdptcorrelations::AcceptedTrueEventId,
                  dptdptcorrelations::TrackacceptedAsOne,
                  dptdptcorrelations::TrackacceptedAsTwo,
                  dptdptcorrelations::Pt,
                  dptdptcorrelations::Eta,
                  dptdptcorrelations::Phi);

using CollisionsEvSelCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>;
using CollisionEvSelCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator;
using CollisionsEvSel = soa::Join<aod::Collisions, aod::EvSels>;
using CollisionEvSel = soa::Join<aod::Collisions, aod::EvSels>::iterator;
using TrackData = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection>::iterator;
using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::ScannedTracks>>;
using FilteredTrackData = Partition<aod::FilteredTracks>::filtered_iterator;
} // namespace aod
} // namespace o2

namespace dptdptcorrelations
{
/* all this is made configurable */
int ptbins = 18;
float ptlow = 0.2, ptup = 2.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int zvtxbins = 40;
float zvtxlow = -10.0, zvtxup = 10.0;
int phibins = 72;
float philow = 0.0;
float phiup = M_PI * 2;
float phibinshift = 0.5;
float etabinwidth = (etaup - etalow) / float(etabins);
float phibinwidth = (phiup - philow) / float(phibins);
int deltaetabins = etabins * 2 - 1;
float deltaetalow = etalow - etaup, deltaetaup = etaup - etalow;
float deltaetabinwidth = (deltaetaup - deltaetalow) / float(deltaetabins);
int deltaphibins = phibins;
float deltaphibinwidth = M_PI * 2 / deltaphibins;
float deltaphilow = 0.0 - deltaphibinwidth / 2.0;
float deltaphiup = M_PI * 2 - deltaphibinwidth / 2.0;

int tracktype = 1;
int trackonecharge = 1;
int tracktwocharge = -1;
bool processpairs = false;
std::string fTaskConfigurationString = "PendingToConfigure";

/// \enum SystemType
/// \brief The type of the system under analysis
enum SystemType {
  kNoSystem = 0, ///< no system defined
  kpp,           ///< **p-p** system
  kpPb,          ///< **p-Pb** system
  kPbp,          ///< **Pb-p** system
  kPbPb,         ///< **Pb-Pb** system
  kXeXe,         ///< **Xe-Xe** system
  knSystems      ///< number of handled systems
};

/// \enum GeneratorType
/// \brief Which kid of generator data is the task addressing
enum GenType {
  kData = 0, ///< actual data, not generated
  kMC,       ///< Generator level and detector level
  kFastMC,   ///< Gererator level but stored dataset
  kOnTheFly, ///< On the fly generator level data
  knGenData  ///< number of different generator data types
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

/// \enum MatchRecoGenSpecies
/// \brief The species considered by the matching tast
enum MatchRecoGenSpecies {
  kCharged = 0, ///< charged particle/track
  kElectron,    ///< electron
  kMuon,        ///< muon
  kPion,        ///< pion
  kKaon,        ///< kaon
  kProton,      ///< proton
  kNoOfSpecies, ///< the number of considered species
  kWrongSpecies = -1
};

const char* speciesName[kNoOfSpecies] = {"h", "e", "mu", "pi", "ka", "p"};

const char* speciesTitle[kNoOfSpecies] = {"", "e", "#mu", "#pi", "K", "p"};

namespace filteranalysistask
{
//============================================================================================
// The DptDptCorrelationsFilterAnalysisTask output objects
//============================================================================================
SystemType fSystem = kNoSystem;
GenType fDataType = kData;
CentMultEstimatorType fCentMultEstimator = kV0M;
analysis::CheckRangeCfg traceDCAOutliers;
bool traceOutOfSpeciesParticles = false;
int recoIdMethod = 0;
bool useOwnTrackSelection = false;
TrackSelection ownTrackSelection = getGlobalTrackSelection();
bool useOwnParticleSelection = false;
bool particleMaxDCAxy = 999.9;
bool particleMaxDCAZ = 999.9;
bool traceCollId0 = false;

TDatabasePDG* fPDG = nullptr;
TH1F* fhCentMultB = nullptr;
TH1F* fhCentMultA = nullptr;
TH1F* fhVertexZB = nullptr;
TH1F* fhVertexZA = nullptr;
TH1F* fhPB = nullptr;
TH1F* fhPA[kNoOfSpecies] = {nullptr};
TH1F* fhPtB = nullptr;
TH1F* fhPtA[kNoOfSpecies] = {nullptr};
TH1F* fhPtPosB = nullptr;
TH1F* fhPtPosA[kNoOfSpecies] = {nullptr};
TH1F* fhPtNegB = nullptr;
TH1F* fhPtNegA[kNoOfSpecies] = {nullptr};

TH1F* fhEtaB = nullptr;
TH1F* fhEtaA = nullptr;

TH1F* fhPhiB = nullptr;
TH1F* fhPhiA = nullptr;

TH1F* fhDCAxyB = nullptr;
TH1F* fhDCAxyA = nullptr;
TH1F* fhFineDCAxyA = nullptr;
TH1F* fhDCAzB = nullptr;
TH1F* fhDCAzA = nullptr;
TH1F* fhFineDCAzA = nullptr;

TH1F* fhTrueCentMultB = nullptr;
TH1F* fhTrueCentMultA = nullptr;
TH1F* fhTrueVertexZB = nullptr;
TH1F* fhTrueVertexZA = nullptr;
TH1F* fhTruePB = nullptr;
TH1F* fhTruePA[kNoOfSpecies] = {nullptr};
TH1F* fhTruePtB = nullptr;
TH1F* fhTruePtA[kNoOfSpecies] = {nullptr};
TH1F* fhTruePtPosB = nullptr;
TH1F* fhTruePtPosA[kNoOfSpecies] = {nullptr};
TH1F* fhTruePtNegB = nullptr;
TH1F* fhTruePtNegA[kNoOfSpecies] = {nullptr};

TH1F* fhTrueEtaB = nullptr;
TH1F* fhTrueEtaA = nullptr;

TH1F* fhTruePhiB = nullptr;
TH1F* fhTruePhiA = nullptr;

TH1F* fhTrueDCAxyB = nullptr;
TH1F* fhTrueDCAxyA = nullptr;
TH1F* fhTrueDCAzB = nullptr;
TH1F* fhTrueDCAxyBid = nullptr;
TH1F* fhTrueDCAzA = nullptr;
} // namespace filteranalysistask

namespace filteranalysistaskqa
{
TH1F* fhTracksOne = nullptr;
TH1F* fhTracksTwo = nullptr;
TH1F* fhTracksOneAndTwo = nullptr;
TH1F* fhTracksNone = nullptr;
TH1F* fhTracksOneUnsel = nullptr;
TH1F* fhTracksTwoUnsel = nullptr;
TH1F* fhTracksOneAndTwoUnsel = nullptr;
TH1F* fhTracksNoneUnsel = nullptr;
TH1F* fhSelectedEvents = nullptr;
} // namespace filteranalysistaskqa

namespace correlationstask
{

/// \enum TrackPairs
/// \brief The track combinations hadled by the class
enum TrackPairs {
  kOO = 0,    ///< one-one pairs
  kOT,        ///< one-two pairs
  kTO,        ///< two-one pairs
  kTT,        ///< two-two pairs
  nTrackPairs ///< the number of track pairs
};
} // namespace correlationstask

/// \brief System type according to configuration string
/// \param sysstr The system configuration string
/// \return The internal code for the passed system string
SystemType getSystemType(std::string const& sysstr)
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
  } else {
    LOGF(fatal, "DptDptCorrelations::getSystemType(). Wrong system type: %d", sysstr.c_str());
  }
  return kPbPb;
}

/// \brief Type of data according to the configuration string
/// \param datastr The data type configuration string
/// \return Internal code for the passed kind of data string
GenType getGenType(std::string const& datastr)
{
  /* we have to figure out how extract the type of data*/
  if (datastr.empty() or (datastr == "data")) {
    return kData;
  } else if (datastr == "MC") {
    return kMC;
  } else if (datastr == "FastMC") {
    return kFastMC;
  } else if (datastr == "OnTheFlyMC") {
    return kOnTheFly;
  } else {
    LOGF(fatal, "DptDptCorrelations::getGenType(). Wrong type of dat: %d", datastr.c_str());
  }
  return kData;
}

template <typename CollisionObject>
bool IsEvtSelected(CollisionObject const& collision, float& centormult)
{
  using namespace filteranalysistask;
  using namespace dptdptcorrelations;

  bool trigsel = false;
  if (fDataType != kData) {
    trigsel = true;
  } else if (collision.alias()[kINT7]) {
    if (collision.sel7()) {
      trigsel = true;
    }
  }

  bool zvtxsel = false;
  /* TODO: vertex quality checks */
  if (zvtxlow < collision.posZ() and collision.posZ() < zvtxup) {
    zvtxsel = true;
  }

  bool centmultsel = false;
  switch (fCentMultEstimator) {
    case kV0M:
      if (collision.centV0M() < 100 and 0 < collision.centV0M()) {
        centormult = collision.centV0M();
        centmultsel = true;
      }
      break;
    default:
      break;
  }
  return trigsel and zvtxsel and centmultsel;
}

template <typename CollisionObject>
bool IsEvtSelectedNoCentMult(CollisionObject const& collision, float& centormult)
{
  using namespace filteranalysistask;
  using namespace dptdptcorrelations;

  bool trigsel = false;
  if (fDataType != kData) {
    trigsel = true;
  } else if (collision.alias()[kINT7]) {
    if (collision.sel7() or collision.sel8()) {
      trigsel = true;
    }
  }

  bool zvtxsel = false;
  /* TODO: vertex quality checks */
  if (zvtxlow < collision.posZ() and collision.posZ() < zvtxup) {
    zvtxsel = true;
  }

  bool centmultsel = false;
  switch (fCentMultEstimator) {
    case kNOCM:
      centormult = 50.0;
      centmultsel = true;
      break;
    default:
      break;
  }
  return trigsel and zvtxsel and centmultsel;
}

template <typename CollisionObject>
bool IsTrueEvtSelected(CollisionObject const& collision, float centormult)
{
  using namespace filteranalysistask;
  using namespace dptdptcorrelations;

  bool zvtxsel = false;
  /* TODO: vertex quality checks */
  if (zvtxlow < collision.posZ() and collision.posZ() < zvtxup) {
    zvtxsel = true;
  }

  bool centmultsel = false;
  if (centormult < 100 and 0 < centormult) {
    centmultsel = true;
  }

  return zvtxsel and centmultsel;
}

template <typename TrackObject>
bool matchTrackType(TrackObject const& track)
{
  using namespace filteranalysistask;

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
        if (track.isGlobalTrack() != 0 || track.isGlobalTrackSDD() != 0) {
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
inline void AcceptTrack(TrackObject const& track, bool& asone, bool& astwo)
{
  using namespace filteranalysistask;

  asone = false;
  astwo = false;

  /* TODO: incorporate a mask in the scanned tracks table for the rejecting track reason */
  if (matchTrackType(track)) {
    if (ptlow < track.pt() and track.pt() < ptup and etalow < track.eta() and track.eta() < etaup) {
      if (((track.sign() > 0) and (trackonecharge > 0)) or ((track.sign() < 0) and (trackonecharge < 0))) {
        asone = true;
      }
      if (((track.sign() > 0) and (tracktwocharge > 0)) or ((track.sign() < 0) and (tracktwocharge < 0))) {
        astwo = true;
      }
    }
  }
}

template <typename ParticleObject, typename MCCollisionObject>
inline void AcceptParticle(ParticleObject& particle, MCCollisionObject const& collision, bool& asone, bool& astwo)
{
  using namespace filteranalysistask;

  asone = false;
  astwo = false;

  float charge = (fPDG->GetParticle(particle.pdgCode())->Charge() / 3 >= 1) ? 1.0 : ((fPDG->GetParticle(particle.pdgCode())->Charge() / 3 <= -1) ? -1.0 : 0.0);

  if (MC::isPhysicalPrimary(particle)) {
    if ((particle.mcCollisionId() == 0) and traceCollId0) {
      LOGF(info, "Particle %d passed isPhysicalPrimary", particle.globalIndex());
    }
    if (useOwnParticleSelection) {
      float dcaxy = TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
      float dcaz = TMath::Abs(particle.vz() - collision.posZ());
      if (not((dcaxy < particleMaxDCAxy) and (dcaz < particleMaxDCAZ))) {
        if ((particle.mcCollisionId() == 0) and traceCollId0) {
          auto currparticle = particle;
          LOGF(info, "Rejecting particle with dcaxy: %.2f and dcaz: %.2f", dcaxy, dcaz);
          LOGF(info, "   assigned collision Id: %d, looping on collision Id: %d", currparticle.mcCollisionId(), collision.globalIndex());
          LOGF(info, "   Collision x: %.5f, y: %.5f, z: %.5f", collision.posX(), collision.posY(), collision.posZ());
          LOGF(info, "   Particle x: %.5f, y: %.5f, z: %.5f", particle.vx(), particle.vy(), particle.vz());
          LOGF(info, "   index: %d, pdg code: %d", currparticle.globalIndex(), currparticle.pdgCode());
          while (currparticle.has_mother0()) {
            LOGF(info, "   mother0 index: %d, mother1 index: %d", currparticle.mother0Id(), currparticle.mother1Id());
            LOGF(info, "  Tracking back mother0 index");
            auto newcurrparticle = currparticle.template mother0_as<aod::McParticles>();
            LOGF(info, "   assigned collision Id: %d, looping on collision Id: %d", newcurrparticle.mcCollisionId(), collision.globalIndex());
            LOGF(info, "   index: %d, pdg code: %d", newcurrparticle.globalIndex(), newcurrparticle.pdgCode());
            LOGF(info, "   Passed  isPhysicalPrimary(): %s", MC::isPhysicalPrimary(newcurrparticle) ? "YES" : "NO");
            currparticle = newcurrparticle;
          }
        }
        return;
      }
    }
    if (ptlow < particle.pt() and particle.pt() < ptup and etalow < particle.eta() and particle.eta() < etaup) {
      if (((charge > 0) and (trackonecharge > 0)) or ((charge < 0) and (trackonecharge < 0))) {
        asone = true;
      }
      if (((charge > 0) and (tracktwocharge > 0)) or ((charge < 0) and (tracktwocharge < 0))) {
        astwo = true;
      }
    }
  } else {
    if ((particle.mcCollisionId() == 0) and traceCollId0) {
      LOGF(info, "Particle %d NOT passed isPhysicalPrimary", particle.globalIndex());
    }
  }
}

template <typename TrackObject>
void fillTrackHistosBeforeSelection(TrackObject const& track)
{
  using namespace filteranalysistask;

  fhPB->Fill(track.p());
  fhPtB->Fill(track.pt());
  fhEtaB->Fill(track.eta());
  fhPhiB->Fill(track.phi());
  if (track.sign() > 0) {
    fhPtPosB->Fill(track.pt());
  } else {
    fhPtNegB->Fill(track.pt());
  }
  fhDCAxyB->Fill(track.dcaXY());
  fhDCAzB->Fill(track.dcaZ());
}

template <typename TrackObject>
void fillTrackHistosAfterSelection(TrackObject const& track, MatchRecoGenSpecies sp)
{
  using namespace filteranalysistask;

  /* the charged species should have been called first so avoid double counting */
  if (sp == kCharged) {
    fhEtaA->Fill(track.eta());
    fhPhiA->Fill(track.phi());
    fhDCAxyA->Fill(track.dcaXY());
    fhDCAzA->Fill(track.dcaZ());
    if (track.dcaXY() < 1.0) {
      fhFineDCAxyA->Fill(track.dcaXY());
    }
    if (track.dcaZ() < 1.0) {
      fhFineDCAzA->Fill(track.dcaZ());
    }
  }
  fhPA[sp]->Fill(track.p());
  fhPtA[sp]->Fill(track.pt());
  if (track.sign() > 0) {
    fhPtPosA[sp]->Fill(track.pt());
  } else {
    fhPtNegA[sp]->Fill(track.pt());
  }
}

template <typename ParticleObject, typename MCCollisionObject>
void fillParticleHistosBeforeSelection(ParticleObject const& particle, MCCollisionObject const& collision, float charge)
{
  using namespace filteranalysistask;

  fhTruePB->Fill(particle.p());
  fhTruePtB->Fill(particle.pt());
  fhTrueEtaB->Fill(particle.eta());
  fhTruePhiB->Fill(particle.phi());
  if (charge > 0) {
    fhTruePtPosB->Fill(particle.pt());
  } else if (charge < 0) {
    fhTruePtNegB->Fill(particle.pt());
  }

  float dcaxy = TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                            (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
  if (traceDCAOutliers.mDoIt and (traceDCAOutliers.mLowValue < dcaxy) and (dcaxy < traceDCAOutliers.mUpValue)) {
    fhTrueDCAxyBid->Fill(TString::Format("%d", particle.pdgCode()).Data(), 1.0);
  }

  fhTrueDCAxyB->Fill(TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                 (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())));
  fhTrueDCAzB->Fill((particle.vz() - collision.posZ()));
}

template <typename ParticleObject, typename MCCollisionObject>
void fillParticleHistosAfterSelection(ParticleObject const& particle, MCCollisionObject const& collision, float charge, MatchRecoGenSpecies sp)
{
  using namespace filteranalysistask;

  /* the charged species should have been called first so avoid double counting */
  if (sp == kCharged) {
    fhTrueEtaA->Fill(particle.eta());
    fhTruePhiA->Fill(particle.phi());
    float dcaxy = TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                              (particle.vy() - collision.posY()) * (particle.vy() - collision.posY()));
    if (traceDCAOutliers.mDoIt and (traceDCAOutliers.mLowValue < dcaxy) and (dcaxy < traceDCAOutliers.mUpValue)) {
      LOGF(info, "DCAxy outlier: Particle with index %d and pdg code %d assigned to MC collision %d, pT: %f, phi: %f, eta: %f",
           particle.globalIndex(), particle.pdgCode(), particle.mcCollisionId(), particle.pt(), particle.phi(), particle.eta());
      LOGF(info, "               With status %d and flags %0X", particle.statusCode(), particle.flags());
    }

    fhTrueDCAxyA->Fill(TMath::Sqrt((particle.vx() - collision.posX()) * (particle.vx() - collision.posX()) +
                                   (particle.vy() - collision.posY()) * (particle.vy() - collision.posY())));
    fhTrueDCAzA->Fill((particle.vz() - collision.posZ()));
  }
  fhTruePA[sp]->Fill(particle.p());
  fhTruePtA[sp]->Fill(particle.pt());
  if (charge > 0) {
    fhTruePtPosA[sp]->Fill(particle.pt());
  } else {
    fhTruePtNegA[sp]->Fill(particle.pt());
  }
}

template <typename TrackObject>
inline MatchRecoGenSpecies IdentifyTrack(TrackObject const& track)
{
  float nsigmas[kNoOfSpecies];
  if (track.p() < 0.8) {
    nsigmas[kCharged] = 999.0f;
    nsigmas[kElectron] = track.tpcNSigmaEl();
    nsigmas[kMuon] = track.tpcNSigmaMu();
    nsigmas[kPion] = track.tpcNSigmaPi();
    nsigmas[kKaon] = track.tpcNSigmaKa();
    nsigmas[kProton] = track.tpcNSigmaPr();
  } else {
    /* introduce require TOF flag */
    if (track.hasTOF()) {
      nsigmas[kCharged] = 999.0f;
      nsigmas[kElectron] = sqrtf(track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl());
      nsigmas[kMuon] = sqrtf(track.tpcNSigmaMu() * track.tpcNSigmaMu() + track.tofNSigmaMu() * track.tofNSigmaMu());
      nsigmas[kPion] = sqrtf(track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi());
      nsigmas[kKaon] = sqrtf(track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa());
      nsigmas[kProton] = sqrtf(track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr());
    } else {
      nsigmas[kCharged] = 999.0f;
      nsigmas[kElectron] = track.tpcNSigmaEl();
      nsigmas[kMuon] = track.tpcNSigmaMu();
      nsigmas[kPion] = track.tpcNSigmaPi();
      nsigmas[kKaon] = track.tpcNSigmaKa();
      nsigmas[kProton] = track.tpcNSigmaPr();
    }
  }
  float min_nsigma = 999.0f;
  MatchRecoGenSpecies sp_min_nsigma;
  for (int sp = 0; sp < kNoOfSpecies; ++sp) {
    if (nsigmas[sp] < min_nsigma) {
      min_nsigma = nsigmas[sp];
      sp_min_nsigma = MatchRecoGenSpecies(sp);
    }
  }
  bool doublematch = false;
  if (min_nsigma < 3.0) {
    for (int sp = 0; (sp < kNoOfSpecies) and not doublematch; ++sp) {
      if (sp != sp_min_nsigma) {
        if (nsigmas[sp] < 3.0) {
          doublematch = true;
        }
      }
    }
    if (doublematch) {
      return kWrongSpecies;
    } else {
      return sp_min_nsigma;
    }
  } else {
    return kWrongSpecies;
  }
}

template <typename ParticleObject>
inline MatchRecoGenSpecies IdentifyParticle(ParticleObject const& particle)
{
  using namespace filteranalysistask;
  constexpr int pdgcodeEl = 11L;
  constexpr int pdgcodeMu = 13L;
  constexpr int pdgcodePi = 211L;
  constexpr int pdgcodeKa = 321L;
  constexpr int pdgcodePr = 2212L;

  int pdgcode = abs(particle.pdgCode());

  switch (pdgcode) {
    case pdgcodeEl:
      return kElectron;
      break;
    case pdgcodeMu:
      return kMuon;
      break;
    case pdgcodePi:
      return kPion;
      break;
    case pdgcodeKa:
      return kKaon;
      break;
    case pdgcodePr:
      return kProton;
      break;

    default:
      if (traceOutOfSpeciesParticles) {
        LOGF(info, "Wrong particle passed selection cuts. PDG code: %d", pdgcode);
      }
      return kWrongSpecies;
      break;
  }
}

} /* end namespace dptdptcorrelations */

// Task for <dpt,dpt> correlations analysis
// FIXME: this should really inherit from AnalysisTask but
//        we need GCC 7.4+ for that

using namespace dptdptcorrelations;

struct DptDptCorrelationsFilterAnalysisTask {
  Configurable<int> cfgTrackType{"trktype", 1, "Type of selected tracks: 0 = no selection, 1 = global tracks FB96, 3 = Run3 tracks. Default 1"};
  Configurable<std::string> cfgCentMultEstimator{"centmultestimator", "V0M", "Centrality/multiplicity estimator detector: V0M, NOCM: none. Default V0M"};
  Configurable<std::string> cfgSystem{"syst", "PbPb", "System: pp, PbPb, Pbp, pPb, XeXe. Default PbPb"};
  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"binning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<o2::analysis::CheckRangeCfg> cfgTraceDCAOutliers{"trackdcaoutliers", {false, 0.0, 0.0}, "Track the generator level DCAxy outliers: false/true, low dcaxy, up dcaxy. Default {false,0.0,0.0}"};
  Configurable<float> cfgTraceOutOfSpeciesParticles{"trackoutparticles", false, "Track the particles which are not e,mu,pi,K,p: false/true. Default false"};
  Configurable<int> cfgRecoIdMethod{"recoidmethod", 0, "Method for identifying reconstructed tracks: 0 PID, 1 mcparticle. Default 0"};
  Configurable<o2::analysis::TrackSelectionCfg> cfgTrackSelection{"tracksel", {false, false, 0, 70, 0.8, 2.4, 3.2}, "Track selection: {useit: true/false, ongen: true/false, tpccls, tpcxrws, tpcxrfc, dcaxy, dcaz}. Default {false,0.70.0.8,2.4,3.2}"};
  Configurable<bool> cfgTraceCollId0{"tracecollid0", false, "Trace particles in collisions id 0. Default false"};

  OutputObj<TList> fOutput{"MatchingRecoGenGlobalInfo", OutputObjHandlingPolicy::AnalysisObject};

  Produces<aod::AcceptedEvents> acceptedevents;
  Produces<aod::ScannedTracks> scannedtracks;
  Produces<aod::AcceptedTrueEvents> acceptedtrueevents;
  Produces<aod::ScannedTrueTracks> scannedtruetracks;

  template <typename TrackListObject, typename CollisionIndex>
  void filterDetectorLevelTracks(TrackListObject const& ftracks, CollisionIndex colix, aod::McParticles const&)
  {
    using namespace filteranalysistask;

    int acceptedtracks = 0;

    for (auto& track : ftracks) {
      bool asone = false;
      bool astwo = false;
      if (not(track.mcParticleId() < 0)) {
        /* correctly reconstructed track */
        /* before track selection */
        fillTrackHistosBeforeSelection(track);

        /* track selection */
        /* tricky because the boolean columns issue */
        AcceptTrack(track, asone, astwo);
        if (asone or astwo) {
          /* the track has been accepted */
          /* fill the charged tracks histograms */
          fillTrackHistosAfterSelection(track, kCharged);
          /* let's identify it */
          MatchRecoGenSpecies sp = kWrongSpecies;
          if (cfgRecoIdMethod == 0) {
            sp = IdentifyTrack(track);
          } else if (cfgRecoIdMethod == 1) {
            sp = IdentifyParticle(track.mcParticle());
          }
          if (sp != kWrongSpecies) {
            /* fill the species histograms */
            fillTrackHistosAfterSelection(track, sp);
          }
          acceptedtracks++;
          scannedtracks(colix, (uint8_t)asone, (uint8_t)astwo, track.pt(), track.eta(), track.phi());
          LOGF(MATCHRECGENLOGTRACKS, "Accepted track with global Id %d and with assigned collision Id %d", track.globalIndex(), track.collisionId());
        }
      }
    }
    LOGF(MATCHRECGENLOGCOLLISIONS, "Accepted %d reconstructed tracks", acceptedtracks);
  }

  template <typename ParticleListObject, typename MCCollisionObject, typename CollisionIndex>
  void filterParticles(ParticleListObject const& particles, MCCollisionObject const& mccollision, CollisionIndex colix)
  {
    using namespace filteranalysistask;

    int acceptedparticles = 0;

    for (auto& particle : particles) {
      float charge = 0.0;
      TParticlePDG* pdgparticle = fPDG->GetParticle(particle.pdgCode());
      if (pdgparticle != nullptr) {
        charge = (pdgparticle->Charge() / 3 >= 1) ? 1.0 : ((pdgparticle->Charge() / 3 <= -1) ? -1.0 : 0.0);
      }

      bool asone = false;
      bool astwo = false;
      if (charge != 0) {
        /* before particle selection */
        fillParticleHistosBeforeSelection(particle, mccollision, charge);

        /* track selection */
        /* tricky because the boolean columns issue */
        AcceptParticle(particle, mccollision, asone, astwo);
        if (asone or astwo) {
          /* the track has been accepted */
          /* fill the charged particle histograms */
          fillParticleHistosAfterSelection(particle, mccollision, charge, kCharged);
          /* let's identify the particle */
          MatchRecoGenSpecies sp = IdentifyParticle(particle);
          if (sp != kWrongSpecies) {
            fillParticleHistosAfterSelection(particle, mccollision, charge, sp);
          }
          acceptedparticles++;
          scannedtruetracks(colix, (uint8_t)asone, (uint8_t)astwo, particle.pt(), particle.eta(), particle.phi());
        }
      } else {
        if ((particle.mcCollisionId() == 0) and traceCollId0) {
          LOGF(info, "Particle %d with fractional charge or equal to zero", particle.globalIndex());
        }
      }
    }
    LOGF(MATCHRECGENLOGCOLLISIONS, "Accepted %d generated particles", acceptedparticles);
  }

  void init(InitContext const&)
  {
    using namespace filteranalysistask;

    LOGF(info, "FilterAnalysisTask::init()");

    /* update with the configurable values */
    /* the binning */
    ptbins = cfgBinning->mPTbins;
    ptlow = cfgBinning->mPTmin;
    ptup = cfgBinning->mPTmax;
    etabins = cfgBinning->mEtabins;
    etalow = cfgBinning->mEtamin;
    etaup = cfgBinning->mEtamax;
    zvtxbins = cfgBinning->mZVtxbins;
    zvtxlow = cfgBinning->mZVtxmin;
    zvtxup = cfgBinning->mZVtxmax;
    /* the track types and combinations */
    tracktype = cfgTrackType.value;
    /* the centrality/multiplicity estimation */
    if (cfgCentMultEstimator->compare("V0M") == 0) {
      fCentMultEstimator = kV0M;
    } else if (cfgCentMultEstimator->compare("NOCM") == 0) {
      fCentMultEstimator = kNOCM;
    } else {
      LOGF(fatal, "Centrality/Multiplicity estimator %s not supported yet", cfgCentMultEstimator->c_str());
    }
    traceDCAOutliers = cfgTraceDCAOutliers;
    traceOutOfSpeciesParticles = cfgTraceOutOfSpeciesParticles;
    recoIdMethod = cfgRecoIdMethod;
    if (cfgTrackSelection->mUseIt) {
      useOwnTrackSelection = true;
      if (cfgTrackSelection->mOnGen) {
        useOwnParticleSelection = true;
        particleMaxDCAxy = cfgTrackSelection->mDCAxy;
        particleMaxDCAZ = cfgTrackSelection->mDCAz;
      }
      ownTrackSelection.SetMinNClustersTPC(cfgTrackSelection->mTPCclusters);
      ownTrackSelection.SetMinNCrossedRowsTPC(cfgTrackSelection->mTPCxRows);
      ownTrackSelection.SetMinNCrossedRowsOverFindableClustersTPC(cfgTrackSelection->mTPCXRoFClusters);
      ownTrackSelection.SetMaxDcaXYPtDep(std::function<float(float)>{});
      ownTrackSelection.SetMaxDcaXY(cfgTrackSelection->mDCAxy);
      ownTrackSelection.SetMaxDcaZ(cfgTrackSelection->mDCAz);
      o2::aod::track::TrackTypeEnum ttype;
      switch (tracktype) {
        case 1:
          ttype = o2::aod::track::Run2Track;
          break;
        case 3:
          ttype = o2::aod::track::Track;
          break;
        default:
          ttype = o2::aod::track::Track;
          break;
      }
      ownTrackSelection.SetTrackType(ttype);
    } else {
      useOwnTrackSelection = false;
    }
    traceCollId0 = cfgTraceCollId0;

    /* if the system type is not known at this time, we have to put the initalization somewhere else */
    fSystem = getSystemType(cfgSystem);
    fDataType = kMC;
    fPDG = TDatabasePDG::Instance();

    /* create the output list which will own the task histograms */
    TList* fOutputList = new TList();
    fOutputList->SetOwner(true);
    fOutput.setObject(fOutputList);

    /* incorporate configuration parameters to the output */
    fOutputList->Add(new TParameter<Int_t>("TrackType", cfgTrackType, 'f'));
    fOutputList->Add(new TParameter<Int_t>("TrackOneCharge", trackonecharge, 'f'));
    fOutputList->Add(new TParameter<Int_t>("TrackTwoCharge", tracktwocharge, 'f'));

    if ((fDataType == kData) or (fDataType == kMC)) {
      /* create the reconstructed data histograms */
      if (fSystem > kPbp) {
        fhCentMultB = new TH1F("CentralityB", "Centrality before cut; centrality (%)", 100, 0, 100);
        fhCentMultA = new TH1F("CentralityA", "Centrality; centrality (%)", 100, 0, 100);
      } else {
        /* for pp, pPb and Pbp systems use multiplicity instead */
        fhCentMultB = new TH1F("MultiplicityB", "Multiplicity before cut; multiplicity (%)", 100, 0, 100);
        fhCentMultA = new TH1F("MultiplicityA", "Multiplicity; multiplicity (%)", 100, 0, 100);
      }

      fhVertexZB = new TH1F("VertexZB", "Vertex Z; z_{vtx}", 60, -15, 15);
      fhVertexZA = new TH1F("VertexZA", "Vertex Z; z_{vtx}", zvtxbins, zvtxlow, zvtxup);

      fhPB = new TH1F("fHistPB", "p distribution for reconstructed before;p (GeV/c);dN/dp (c/GeV)", 100, 0.0, 15.0);
      fhPtB = new TH1F("fHistPtB", "p_{T} distribution for reconstructed before;p_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhPtPosB = new TH1F("fHistPtPosB", "P_{T} distribution for reconstructed (#plus) before;P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhPtNegB = new TH1F("fHistPtNegB", "P_{T} distribution for reconstructed (#minus) before;P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhEtaB = new TH1F("fHistEtaB", "#eta distribution for reconstructed before;#eta;counts", 40, -2.0, 2.0);
      fhEtaA = new TH1F("fHistEtaA", "#eta distribution for reconstructed;#eta;counts", etabins, etalow, etaup);
      fhPhiB = new TH1F("fHistPhiB", "#phi distribution for reconstructed before;#phi;counts", 360, 0.0, 2 * M_PI);
      fhPhiA = new TH1F("fHistPhiA", "#phi distribution for reconstructed;#phi;counts", 360, 0.0, 2 * M_PI);
      fhDCAxyB = new TH1F("DCAxyB", "DCA_{xy} distribution for reconstructed before;DCA_{xy} (cm);counts", 1000, -4.0, 4.0);
      fhDCAxyA = new TH1F("DCAxyA", "DCA_{xy} distribution for reconstructed;DCA_{xy} (cm);counts", 1000, -4., 4.0);
      fhFineDCAxyA = new TH1F("FineDCAxyA", "DCA_{xy} distribution for reconstructed;DCA_{xy} (cm);counts", 4000, -1.0, 1.0);
      fhDCAzB = new TH1F("DCAzB", "DCA_{z} distribution for reconstructed before;DCA_{z} (cm);counts", 1000, -4.0, 4.0);
      fhDCAzA = new TH1F("DCAzA", "DCA_{z} distribution for reconstructed;DCA_{z} (cm);counts", 1000, -4.0, 4.0);
      fhFineDCAzA = new TH1F("FineDCAzA", "DCA_{z} distribution for reconstructed;DCA_{z} (cm);counts", 4000, -1.0, 1.0);

      for (int sp = 0; sp < kNoOfSpecies; ++sp) {
        fhPA[sp] = new TH1F(TString::Format("fHistPA_%s", speciesName[sp]).Data(),
                            TString::Format("p distribution for reconstructed %s;p (GeV/c);dN/dp (c/GeV)", speciesTitle[sp]).Data(),
                            ptbins, ptlow, ptup);
        fhPtA[sp] = new TH1F(TString::Format("fHistPtA_%s", speciesName[sp]),
                             TString::Format("p_{T} distribution for reconstructed %s;p_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                             ptbins, ptlow, ptup);
        fhPtPosA[sp] = new TH1F(TString::Format("fHistPtPosA_%s", speciesName[sp]),
                                TString::Format("P_{T} distribution for reconstructed  %s^{#plus};P_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                ptbins, ptlow, ptup);
        fhPtNegA[sp] = new TH1F(TString::Format("fHistPtNegA_%s", speciesName[sp]),
                                TString::Format("P_{T} distribution for reconstructed  %s^{#minus};P_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                ptbins, ptlow, ptup);
      }

      /* add the hstograms to the output list */
      fOutputList->Add(fhCentMultB);
      fOutputList->Add(fhCentMultA);
      fOutputList->Add(fhVertexZB);
      fOutputList->Add(fhVertexZA);
      fOutputList->Add(fhPB);
      fOutputList->Add(fhPtB);
      fOutputList->Add(fhPtPosB);
      fOutputList->Add(fhPtNegB);
      fOutputList->Add(fhEtaB);
      fOutputList->Add(fhEtaA);
      fOutputList->Add(fhPhiB);
      fOutputList->Add(fhPhiA);
      fOutputList->Add(fhDCAxyB);
      fOutputList->Add(fhDCAxyA);
      fOutputList->Add(fhFineDCAxyA);
      fOutputList->Add(fhDCAzB);
      fOutputList->Add(fhDCAzA);
      fOutputList->Add(fhFineDCAzA);

      for (int sp = 0; sp < kNoOfSpecies; ++sp) {
        fOutputList->Add(fhPA[sp]);
        fOutputList->Add(fhPtA[sp]);
        fOutputList->Add(fhPtPosA[sp]);
        fOutputList->Add(fhPtNegA[sp]);
      }
    }

    if (fDataType != kData) {
      /* create the true data histograms */
      if (fSystem > kPbp) {
        fhTrueCentMultB = new TH1F("TrueCentralityB", "Centrality before (truth); centrality (%)", 100, 0, 100);
        fhTrueCentMultA = new TH1F("TrueCentralityA", "Centrality (truth); centrality (%)", 100, 0, 100);
      } else {
        /* for pp, pPb and Pbp systems use multiplicity instead */
        fhTrueCentMultB = new TH1F("TrueMultiplicityB", "Multiplicity before (truth); multiplicity (%)", 100, 0, 100);
        fhTrueCentMultA = new TH1F("TrueMultiplicityA", "Multiplicity (truth); multiplicity (%)", 100, 0, 100);
      }

      fhTrueVertexZB = new TH1F("TrueVertexZB", "Vertex Z before (truth); z_{vtx}", 60, -15, 15);
      fhTrueVertexZA = new TH1F("TrueVertexZA", "Vertex Z (truth); z_{vtx}", zvtxbins, zvtxlow, zvtxup);

      fhTruePB = new TH1F("fTrueHistPB", "p distribution before (truth);p (GeV/c);dN/dp (c/GeV)", 100, 0.0, 15.0);
      fhTruePtB = new TH1F("fTrueHistPtB", "p_{T} distribution before (truth);p_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhTruePtPosB = new TH1F("fTrueHistPtPosB", "P_{T} distribution (#plus) before (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhTruePtNegB = new TH1F("fTrueHistPtNegB", "P_{T} distribution (#minus) before (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", 100, 0.0, 15.0);
      fhTrueEtaB = new TH1F("fTrueHistEtaB", "#eta distribution before (truth);#eta;counts", 40, -2.0, 2.0);
      fhTrueEtaA = new TH1F("fTrueHistEtaA", "#eta distribution (truth);#eta;counts", etabins, etalow, etaup);
      fhTruePhiB = new TH1F("fTrueHistPhiB", "#phi distribution before (truth);#phi;counts", 360, 0.0, 2 * M_PI);
      fhTruePhiA = new TH1F("fTrueHistPhiA", "#phi distribution (truth);#phi;counts", 360, 0.0, 2 * M_PI);
      fhTrueDCAxyB = new TH1F("TrueDCAxyB", "DCA_{xy} distribution for generated before;DCA_{xy} (cm);counts", 1000, -4.0, 4.0);
      if (traceDCAOutliers.mDoIt) {
        fhTrueDCAxyBid = new TH1F("PDGCodeDCAxyB",
                                  TString::Format("PDG code within %.2f<|DCA_{#it{xy}}|<%.2f; PDG code", traceDCAOutliers.mLowValue, traceDCAOutliers.mUpValue).Data(),
                                  100, 0.5, 100.5);
      }
      fhTrueDCAxyA = new TH1F("TrueDCAxyA", "DCA_{xy} distribution for generated;DCA_{xy};counts (cm)", 1000, -4., 4.0);
      fhTrueDCAzB = new TH1F("TrueDCAzB", "DCA_{z} distribution for generated before;DCA_{z} (cm);counts", 1000, -4.0, 4.0);
      fhTrueDCAzA = new TH1F("TrueDCAzA", "DCA_{z} distribution for generated;DCA_{z} (cm);counts", 1000, -4.0, 4.0);

      for (int sp = 0; sp < kNoOfSpecies; ++sp) {
        fhTruePA[sp] = new TH1F(TString::Format("fTrueHistPA_%s", speciesName[sp]).Data(),
                                TString::Format("p distribution %s (truth);p (GeV/c);dN/dp (c/GeV)", speciesTitle[sp]).Data(),
                                ptbins, ptlow, ptup);
        fhTruePtA[sp] = new TH1F(TString::Format("fTrueHistPtA_%s", speciesName[sp]),
                                 TString::Format("p_{T} distribution %s (truth);p_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                 ptbins, ptlow, ptup);
        fhTruePtPosA[sp] = new TH1F(TString::Format("fTrueHistPtPosA_%s", speciesName[sp]),
                                    TString::Format("P_{T} distribution %s^{#plus} (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                    ptbins, ptlow, ptup);
        fhTruePtNegA[sp] = new TH1F(TString::Format("fTrueHistPtNegA_%s", speciesName[sp]),
                                    TString::Format("P_{T} distribution %s^{#minus} (truth);P_{T} (GeV/c);dN/dP_{T} (c/GeV)", speciesTitle[sp]).Data(),
                                    ptbins, ptlow, ptup);
      }

      /* add the hstograms to the output list */
      fOutputList->Add(fhTrueCentMultB);
      fOutputList->Add(fhTrueCentMultA);
      fOutputList->Add(fhTrueVertexZB);
      fOutputList->Add(fhTrueVertexZA);
      fOutputList->Add(fhTruePB);
      fOutputList->Add(fhTruePtB);
      fOutputList->Add(fhTruePtPosB);
      fOutputList->Add(fhTruePtNegB);
      fOutputList->Add(fhTrueEtaB);
      fOutputList->Add(fhTrueEtaA);
      fOutputList->Add(fhTruePhiB);
      fOutputList->Add(fhTruePhiA);
      fOutputList->Add(fhTrueDCAxyB);
      if (traceDCAOutliers.mDoIt) {
        fOutputList->Add(fhTrueDCAxyBid);
      }
      fOutputList->Add(fhTrueDCAxyA);
      fOutputList->Add(fhTrueDCAzB);
      fOutputList->Add(fhTrueDCAzA);

      for (int sp = 0; sp < kNoOfSpecies; ++sp) {
        fOutputList->Add(fhTruePA[sp]);
        fOutputList->Add(fhTruePtA[sp]);
        fOutputList->Add(fhTruePtPosA[sp]);
        fOutputList->Add(fhTruePtNegA[sp]);
      }
    }
  }

  using FullTracksPID = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksCov, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;

  void processWithCentDetectorLevel(aod::CollisionEvSelCent const& collision, FullTracksPID const& ftracks, aod::McParticles const& particles)
  {
    using namespace filteranalysistask;

    LOGF(MATCHRECGENLOGCOLLISIONS, "FilterAnalysisTask::processWithCentDetectorLevel(). New collision with %d tracks", ftracks.size());

    fhCentMultB->Fill(collision.centV0M());
    fhVertexZB->Fill(collision.posZ());
    bool acceptedevent = false;
    float centormult = -100.0;
    if (IsEvtSelected(collision, centormult)) {
      acceptedevent = true;
      fhCentMultA->Fill(collision.centV0M());
      fhVertexZA->Fill(collision.posZ());
      acceptedevents(collision.bcId(), collision.posZ(), (uint8_t)acceptedevent, centormult);

      filterDetectorLevelTracks(ftracks, acceptedevents.lastIndex(), particles);
    } else {
      acceptedevents(collision.bcId(), collision.posZ(), (uint8_t)acceptedevent, centormult);
      for (auto& track : ftracks) {
        scannedtracks(acceptedevents.lastIndex(), (uint8_t) false, (uint8_t) false, track.pt(), track.eta(), track.phi());
      }
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsFilterAnalysisTask, processWithCentDetectorLevel, "Process MC detector level with centrality", false);

  void processWithoutCentDetectorLevel(aod::CollisionEvSel const& collision, FullTracksPID const& ftracks, aod::McParticles const& particles)
  {
    using namespace filteranalysistask;

    LOGF(MATCHRECGENLOGCOLLISIONS, "FilterAnalysisTask::processWithoutCentDetectorLevel(). New collision with collision id %d and with %d tracks", collision.bcId(), ftracks.size());

    /* the task does not have access to either centrality nor multiplicity
       classes information, so it has to live without it.
       For the time being we assign a value of 50% */
    fhCentMultB->Fill(50.0);
    fhVertexZB->Fill(collision.posZ());
    bool acceptedevent = false;
    float centormult = -100.0;
    if (IsEvtSelectedNoCentMult(collision, centormult)) {
      acceptedevent = true;
      fhCentMultA->Fill(50.0);
      fhVertexZA->Fill(collision.posZ());
      acceptedevents(collision.bcId(), collision.posZ(), (uint8_t)acceptedevent, centormult);

      LOGF(MATCHRECGENLOGCOLLISIONS, "Accepted collision with BC id %d and collision Id %d", collision.bcId(), collision.globalIndex());
      filterDetectorLevelTracks(ftracks, acceptedevents.lastIndex(), particles);
    } else {
      acceptedevents(collision.bcId(), collision.posZ(), (uint8_t)acceptedevent, centormult);
      for (auto& track : ftracks) {
        scannedtracks(acceptedevents.lastIndex(), (uint8_t) false, (uint8_t) false, track.pt(), track.eta(), track.phi());
      }
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsFilterAnalysisTask, processWithoutCentDetectorLevel, "Process MC detector level without centrality", false);

  void processWithCentGeneratorLevel(aod::McCollision const& mccollision,
                                     soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentV0Ms> const& collisions,
                                     aod::McParticles const& mcparticles)
  {
    using namespace filteranalysistask;

    LOGF(MATCHRECGENLOGCOLLISIONS, "FilterAnalysisTask::processWithCentGeneratorLevel(). New generated collision %d reconstructed collisions and %d particles", collisions.size(), mcparticles.size());

    /* TODO: in here we have to decide what to do in the following cases
       - On the fly production -> clearly we will need a different process
       - reconstructed collisions without generated associated -> we need a different process or a different signature
       - multiplicity/centrality classes extracted from the reconstructed collision but then
       - generated collision without associated reconstructed collision: how to extract mutliplicity/centrality classes?
       - generated collision with several associated reconstructed collisions: from which to extract multiplicity/centrality classes?
    */
    if (collisions.size() > 1) {
      LOGF(error, "FilterAnalysisTask::processWithCentGeneratorLevel(). Generated collision with more than one reconstructed collisions. Processing only the first for centrality/multiplicity classes extraction");
    }

    for (auto& collision : collisions) {
      float cent = collision.centV0M();
      fhTrueCentMultB->Fill(cent);
      fhTrueVertexZB->Fill(mccollision.posZ());

      bool acceptedevent = false;
      if (IsTrueEvtSelected(mccollision, cent)) {
        acceptedevent = true;
        fhTrueCentMultA->Fill(cent);
        fhTrueVertexZA->Fill(mccollision.posZ());
        acceptedtrueevents(mccollision.bcId(), mccollision.posZ(), (uint8_t)acceptedevent, cent);

        filterParticles(mcparticles, collision, acceptedtrueevents.lastIndex());
      } else {
        acceptedtrueevents(mccollision.bcId(), mccollision.posZ(), (uint8_t)acceptedevent, cent);
        for (auto& particle : mcparticles) {
          scannedtruetracks(acceptedtrueevents.lastIndex(), (uint8_t) false, (uint8_t) false, particle.pt(), particle.eta(), particle.phi());
        }
      }
      break; /* TODO: only processing the first reconstructed collision for centrality/multiplicity class estimation */
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsFilterAnalysisTask, processWithCentGeneratorLevel, "Process generated with centrality", false);

  void processWithoutCentGeneratorLevel(aod::McCollision const& mccollision,
                                        aod::McParticles const& mcparticles)
  {
    using namespace filteranalysistask;

    LOGF(MATCHRECGENLOGCOLLISIONS, "FilterAnalysisTask::processWithoutCentGeneratorLevel(). New generated collision with %d particles", mcparticles.size());

    /* the task does not have access to either centrality nor multiplicity
       classes information, so it has to live without it.
       For the time being we assign a value of 50% */
    fhTrueCentMultB->Fill(50.0);
    fhTrueVertexZB->Fill(mccollision.posZ());

    bool acceptedevent = false;
    float centormult = 50.0;
    if (IsTrueEvtSelected(mccollision, centormult)) {
      acceptedevent = true;
      fhTrueCentMultA->Fill(centormult);
      fhTrueVertexZA->Fill(mccollision.posZ());
      acceptedtrueevents(mccollision.bcId(), mccollision.posZ(), (uint8_t)acceptedevent, centormult);

      filterParticles(mcparticles, mccollision, acceptedtrueevents.lastIndex());
    } else {
      acceptedtrueevents(mccollision.bcId(), mccollision.posZ(), (uint8_t)acceptedevent, centormult);
      for (auto& particle : mcparticles) {
        scannedtruetracks(acceptedtrueevents.lastIndex(), (uint8_t) false, (uint8_t) false, particle.pt(), particle.eta(), particle.phi());
      }
    }
  }
  PROCESS_SWITCH(DptDptCorrelationsFilterAnalysisTask, processWithoutCentGeneratorLevel, "Process generated without centrality", false);
};

// Checking the filtered tables
/* it seems we cannot use a base class task */
// struct TracksAndEventClassificationQABase {

void initQATask(InitContext const&, TList* outlst)
{
  using namespace filteranalysistaskqa;

  fhTracksOne = new TH1F("TracksOne", "Tracks as track one;number of tracks;events", 1500, 0.0, 1500.0);
  fhTracksTwo = new TH1F("TracksTwo", "Tracks as track two;number of tracks;events", 1500, 0.0, 1500.0);
  fhTracksOneAndTwo = new TH1F("TracksOneAndTwo", "Tracks as track one and as track two;number of tracks;events", 1500, 0.0, 1500.0);
  fhTracksNone = new TH1F("TracksNone", "Not selected tracks;number of tracks;events", 1500, 0.0, 1500.0);
  fhTracksOneUnsel = new TH1F("TracksOneUnsel", "Tracks as track one;number of tracks;events", 1500, 0.0, 1500.0);
  fhTracksTwoUnsel = new TH1F("TracksTwoUnsel", "Tracks as track two;number of tracks;events", 1500, 0.0, 1500.0);
  fhTracksOneAndTwoUnsel = new TH1F("TracksOneAndTwoUnsel", "Tracks as track one and as track two;number of tracks;events", 1500, 0.0, 1500.0);
  fhTracksNoneUnsel = new TH1F("TracksNoneUnsel", "Not selected tracks;number of tracks;events", 1500, 0.0, 1500.0);
  fhSelectedEvents = new TH1F("SelectedEvents", "Selected events;;events", 2, 0.0, 2.0);
  fhSelectedEvents->GetXaxis()->SetBinLabel(1, "Not selected events");
  fhSelectedEvents->GetXaxis()->SetBinLabel(2, "Selected events");

  outlst->Add(fhTracksOne);
  outlst->Add(fhTracksTwo);
  outlst->Add(fhTracksOneAndTwo);
  outlst->Add(fhTracksNone);
  outlst->Add(fhTracksOneUnsel);
  outlst->Add(fhTracksTwoUnsel);
  outlst->Add(fhTracksOneAndTwoUnsel);
  outlst->Add(fhTracksNoneUnsel);
  outlst->Add(fhSelectedEvents);
}

template <typename FilteredCollision, typename FilteredTracks>
void processQATask(FilteredCollision const& collision,
                   FilteredTracks const& tracks)
{
  using namespace filteranalysistaskqa;

  if (collision.eventaccepted() != (uint8_t) true) {
    fhSelectedEvents->Fill(0.5);
  } else {
    fhSelectedEvents->Fill(1.5);
  }

  int ntracks_one = 0;
  int ntracks_two = 0;
  int ntracks_one_and_two = 0;
  int ntracks_none = 0;
  for (auto& track : tracks) {
    if ((track.trackacceptedasone() != (uint8_t) true) and (track.trackacceptedastwo() != (uint8_t) true)) {
      ntracks_none++;
    }
    if ((track.trackacceptedasone() == (uint8_t) true) and (track.trackacceptedastwo() == (uint8_t) true)) {
      ntracks_one_and_two++;
    }
    if (track.trackacceptedasone() == (uint8_t) true) {
      ntracks_one++;
    }
    if (track.trackacceptedastwo() == (uint8_t) true) {
      ntracks_two++;
    }
  }
  if (collision.eventaccepted() != (uint8_t) true) {
    /* control for non selected events */
    fhTracksOneUnsel->Fill(ntracks_one);
    fhTracksTwoUnsel->Fill(ntracks_two);
    fhTracksNoneUnsel->Fill(ntracks_none);
    fhTracksOneAndTwoUnsel->Fill(ntracks_one_and_two);
  } else {
    fhTracksOne->Fill(ntracks_one);
    fhTracksTwo->Fill(ntracks_two);
    fhTracksNone->Fill(ntracks_none);
    fhTracksOneAndTwo->Fill(ntracks_one_and_two);
  }
}
// };

/* it seems we cannot use a base class task */
// struct TracksAndEventClassificationQARec : TracksAndEventClassificationQABase {
struct TracksAndEventClassificationQARec {
  OutputObj<TList> fOutput{"FliterTaskRecoQA", OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const& context)
  {
    TList* fOutputList = new TList();
    fOutputList->SetName("FilterTaskRecoQA");
    fOutputList->SetOwner(true);
    fOutput.setObject(fOutputList);

    initQATask(context, fOutputList);
  }

  Filter onlyacceptedevents = (aod::dptdptcorrelations::eventaccepted == (uint8_t) true);
  Filter onlyacceptedtracks = ((aod::dptdptcorrelations::trackacceptedasone == (uint8_t) true) or (aod::dptdptcorrelations::trackacceptedastwo == (uint8_t) true));

  void process(soa::Filtered<aod::AcceptedEvents>::iterator const& collision, soa::Filtered<aod::ScannedTracks> const& tracks)
  {
    LOGF(MATCHRECGENLOGCOLLISIONS, "New filtered collision with BC id %d and with %d accepted tracks", collision.bcId(), tracks.size());
    processQATask(collision, tracks);
  }
};

/* it seems we cannot use a base class task */
// struct TracksAndEventClassificationQAGen : TracksAndEventClassificationQABase {
struct TracksAndEventClassificationQAGen {
  OutputObj<TList> fOutput{"FliterTaskGenQA", OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const& context)
  {
    TList* fOutputList = new TList();
    fOutputList->SetName("FilterTaskGenQA");
    fOutputList->SetOwner(true);
    fOutput.setObject(fOutputList);

    initQATask(context, fOutputList);
  }

  Filter onlyacceptedevents = (aod::dptdptcorrelations::eventaccepted == (uint8_t) true);
  Filter onlyacceptedtracks = ((aod::dptdptcorrelations::trackacceptedasone == (uint8_t) true) or (aod::dptdptcorrelations::trackacceptedastwo == (uint8_t) true));

  void process(soa::Filtered<aod::AcceptedTrueEvents>::iterator const& collision, soa::Filtered<aod::ScannedTrueTracks> const& tracks)
  {
    LOGF(MATCHRECGENLOGCOLLISIONS, "New filtered generated collision with BC id %d and with %d accepted tracks", collision.bcId(), tracks.size());
    processQATask(collision, tracks);
  }
};

namespace recogenmap
{
std::vector<std::vector<int64_t>> mclabelpos[2];
std::vector<std::vector<int64_t>> mclabelneg[2];
} // namespace recogenmap

/// \brief Checks the correspondence generator level <=> detector level
struct CheckGeneratorLevelVsDetectorLevel {
  Configurable<int> cfgTrackType{"trktype", 1, "Type of selected tracks: 0 = no selection, 1 = global tracks FB96"};
  Configurable<std::string> cfgCentMultEstimator{"centmultestimator", "V0M", "Centrality/multiplicity estimator detector:  V0M, NOCM: none. Default V0M"};
  Configurable<std::string> cfgDataType{"datatype", "data", "Data type: data, MC, FastMC, OnTheFlyMC. Default data"};
  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"binning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<bool> cfgTrackMultiRec{"trackmultirec", false, "Track muli-reconstructed particles: true, false. Default false"};
  Configurable<bool> cfgTrackCollAssoc{"trackcollassoc", false, "Track collision id association, track-mcparticle-mccollision vs. track-collision-mccollision: true, false. Default false"};

  HistogramRegistry histos{"RecoGenHistograms", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  TDatabasePDG* fPDG = nullptr;
  typedef enum { kBEFORE = 0,
                 kAFTER } beforeafterselection;
  typedef enum { kPOSITIVE = 0,
                 kNEGATIVE } colllabelsign;
  enum { kMATCH = 0,
         kDONTMATCH };

  void init(InitContext const& context)
  {
    using namespace dptdptcorrelations;
    using namespace dptdptcorrelations::filteranalysistask;

    /* update with the configurable values */
    /* the binning */
    ptbins = cfgBinning->mPTbins;
    ptlow = cfgBinning->mPTmin;
    ptup = cfgBinning->mPTmax;
    etabins = cfgBinning->mEtabins;
    etalow = cfgBinning->mEtamin;
    etaup = cfgBinning->mEtamax;
    zvtxbins = cfgBinning->mZVtxbins;
    zvtxlow = cfgBinning->mZVtxmin;
    zvtxup = cfgBinning->mZVtxmax;
    /* the track types and combinations */
    tracktype = cfgTrackType.value;
    /* the centrality/multiplicity estimation */
    if (cfgCentMultEstimator->compare("V0M") == 0) {
      fCentMultEstimator = kV0M;
    } else if (cfgCentMultEstimator->compare("NOCM") == 0) {
      fCentMultEstimator = kNOCM;
    } else {
      LOGF(fatal, "Centrality/Multiplicity estimator %s not supported yet", cfgCentMultEstimator->c_str());
    }
    fDataType = getGenType(cfgDataType);

    constexpr float TWOPI = 2.0f * static_cast<float>(M_PI);
    fPDG = TDatabasePDG::Instance();

    AxisSpec deltaEta = {100, -2, 2, "#Delta#eta"};
    AxisSpec deltaPhi = {100, 0, TWOPI, "#Delta#varphi (rad)"};
    AxisSpec deltaPt = {1000, 0, 4, "#Delta#it{p}_{T} (GeV/#it{c})"};
    AxisSpec mrectimes = {11, -0.5f, 10.5f, "##/particle"};
    AxisSpec detectors = {32, -0.5, 31.5, "Detectors"};
    std::vector<std::string> detectorlbls = {"", "ITS", "TPC", "ITS+TPC", "TRD", "ITS+TRD", "TPC+TRD", "ITS+TPC+TRD",
                                             "TOF", "ITS+TOF", "TPC+TOF", "ITS+TPC+TOF", "TRD+TOF", "ITS+TRD+TOF", "TPC+TRD+TOF", "ITS+TPC+TRD+TOF",
                                             "UNKN", "ITS+UNKN", "TPC+UNKN", "ITS+TPC+UNKN", "TRD+UNKN", "ITS+TRD+UNKN", "TPC+TRD+UNKN", "ITS+TPC+TRD+UNKN",
                                             "TOF+UNKN", "ITS+TOF+UNKN", "TPC+TOF+UNKN", "ITS+TPC+TOF+UNKN", "TRD+TOF+UNKN", "ITS+TRD+TOF+UNKN", "TPC+TRD+TOF+UNKN", "ITS+TPC+TRD+TOF+UNKN"};
    std::vector<std::string> matchlbs = {"match", "don't match"};

    histos.add("before/positivecolid/mrDeltaEta", "#Delta#eta multirec tracks", kTH1F, {deltaEta});
    histos.add("before/positivecolid/mrDeltaPhi", "#Delta#varphi multirec tracks", kTH1F, {deltaPhi});
    histos.add("before/positivecolid/mrDeltaPt", "#Delta#it{p}_{T} multirec tracks", kTH1F, {deltaPt});
    histos.add("before/positivecolid/multirec", "Multiple reconstruction", kTH1F, {mrectimes});
    histos.add("before/positivecolid/genrecoeta", "#eta Generated vs reconstructed", kTH2F, {{100, -1.0, 1.0, "#eta reco"}, {100, -1.0, 1.0, "#eta gen"}});
    histos.add("before/positivecolid/genrecophi", "#varphi Generated vs reconstructed", kTH2F, {{100, 0, TWOPI, "#varphi (rad) reco"}, {100, 0, TWOPI, "#varphi (rad) gen"}});
    histos.add("before/positivecolid/genrecopt", "#it{p}_{T} Generated vs reconstructed", kTH2F, {{1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c}) reco"}, {1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c}) gen"}});
    histos.add("before/positivecolid/detectormap", "Active detectors", kTH1F, {detectors});
    histos.add("before/positivecolid/matchcollid", "particle MC coll Id <=> track coll MC coll Id", kTH1F, {{2, 0.0, 2.0}});
    histos.add("before/positivecolid/genrecomreta", "#eta Generated vs reconstructed (mr)", kTH2F, {{100, -1.0, 1.0, "#eta reco"}, {100, -1.0, 1.0, "#eta gen"}});
    histos.add("before/positivecolid/genrecomrphi", "#varphi Generated vs reconstructed (mr)", kTH2F, {{100, 0, TWOPI, "#varphi (rad) reco"}, {100, 0, TWOPI, "#varphi (rad) gen"}});
    histos.add("before/positivecolid/genrecomrpt", "#it{p}_{T} Generated vs reconstructed (mr)", kTH2F, {{1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c}) reco"}, {1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c}) gen"}});
    histos.add("before/positivecolid/recomreta", "#eta Reconstructed (mr)", kTH1F, {{100, -1.0, 1.0, "#eta"}});
    histos.add("before/positivecolid/recomrphi", "#varphi Reconstructed (mr)", kTH1F, {{100, 0, TWOPI, "#varphi (rad)"}});
    histos.add("before/positivecolid/recomrpt", "#it{p}_{T} Reconstructed (mr)", kTH1F, {{1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("before/positivecolid/detectormapmr", "Active detectors (mr)", kTH1F, {detectors});
    histos.add("before/positivecolid/matchcollidmr", "particle MC coll Id <=> track coll MC coll Id (mr)", kTH1F, {{2, 0.0, 2.0}});
    histos.add("before/positivecolid/dcaxy", "DCA_{xy} Reconstructed", kTH1F, {{1000, -4.0, 4.0, "DCA_{xy} (cm)"}});
    histos.add("before/positivecolid/dcaz", "DCA_{z} Reconstructed", kTH1F, {{1000, -4.0, 4.0, "DCA_{z} (cm)"}});
    histos.add("before/positivecolid/finedcaxy", "DCA_{xy} Reconstructed", kTH1F, {{2000, -1.0, 1.0, "DCA_{xy} (cm)"}});
    histos.add("before/positivecolid/finedcaz", "DCA_{z} Reconstructed", kTH1F, {{2000, -1.0, 1.0, "DCA_{z} (cm)"}});
    histos.add("before/positivecolid/dcaxymr", "DCA_{xy} Reconstructed (mr)", kTH1F, {{1000, -4.0, 4.0, "DCA_{xy} (cm)"}});
    histos.add("before/positivecolid/dcazmr", "DCA_{z} Reconstructed (mr)", kTH1F, {{1000, -4.0, 4.0, "DCA_{z} (cm)"}});
    histos.add("before/positivecolid/finedcaxymr", "DCA_{xy} Reconstructed (mr)", kTH1F, {{2000, -1.0, 1.0, "DCA_{xy} (cm)"}});
    histos.add("before/positivecolid/finedcazmr", "DCA_{z} Reconstructed (mr)", kTH1F, {{2000, -1.0, 1.0, "DCA_{z} (cm)"}});
    for (int i = 0; i < detectorlbls.size(); ++i) {
      histos.get<TH1>(HIST("before/positivecolid/detectormap"))->GetXaxis()->SetBinLabel(i + 1, detectorlbls[i].c_str());
      histos.get<TH1>(HIST("before/positivecolid/detectormapmr"))->GetXaxis()->SetBinLabel(i + 1, detectorlbls[i].c_str());
    }
    for (int i = 0; i < matchlbs.size(); ++i) {
      histos.get<TH1>(HIST("before/positivecolid/matchcollid"))->GetXaxis()->SetBinLabel(i + 1, matchlbs[i].c_str());
      histos.get<TH1>(HIST("before/positivecolid/matchcollidmr"))->GetXaxis()->SetBinLabel(i + 1, matchlbs[i].c_str());
    }

    /* clone the set for the other cases */
    histos.addClone("before/positivecolid/", "after/positivecolid/");
    histos.addClone("before/positivecolid/", "before/negativecolid/");
    histos.addClone("before/positivecolid/", "after/negativecolid/");
    histos.add("after/positivecolid/pdgcodemr", "PDG code x-collision multi-reconstructed", kTH1F, {{100, 0.5, 100.5, "PDG code"}});
  }

  template <beforeafterselection ba, colllabelsign collsign, typename TracskListObject, typename ParticlesListObject, typename CollisionsListObject>
  void collectData(TracskListObject const& tracks, ParticlesListObject const& mcParticles, CollisionsListObject const& colls)
  {
    using namespace recogenmap;

    static constexpr std::string_view dir[] = {"before/", "after/"};
    static constexpr std::string_view colldir[] = {"positivecolid/", "negativecolid/"};
    constexpr float TWOPI = 2.0F * static_cast<float>(M_PI);

    int nrec_poslabel = 0;
    int nrec_neglabel = 0;
    int nrec_poslabel_crosscoll = 0;

    for (int ixpart = 0; ixpart < mcParticles.size(); ++ixpart) {
      auto particle = mcParticles.iteratorAt(ixpart);
      /* multireconstructed tracks only for positive labels */
      int nrec = mclabelpos[collsign][ixpart].size();
      nrec_poslabel += mclabelpos[collsign][ixpart].size();
      nrec_neglabel += mclabelneg[collsign][ixpart].size();

      if (nrec > 1) {
        /* multireconstruction only from positive labels */
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("multirec"), nrec);

        if (collsign == kPOSITIVE) {
          /* check the cross collision reconstruction */
          bool crosscollfound = false;
          for (int i = 0; (i < mclabelpos[collsign][ixpart].size()) and not crosscollfound; ++i) {
            for (int j = i + 1; (j < mclabelpos[collsign][ixpart].size()) and not crosscollfound; ++j) {
              auto track1 = tracks.iteratorAt(mclabelpos[collsign][ixpart][i]);
              auto track2 = tracks.iteratorAt(mclabelpos[collsign][ixpart][j]);

              if (track1.collisionId() != track2.collisionId()) {
                nrec_poslabel_crosscoll++;
                crosscollfound = true;
              }
            }
          }
          if (crosscollfound and (ba == kAFTER)) {
            if (cfgTrackMultiRec) {
              LOGF(info, "BEGIN multi-reconstructed: ==================================================================");
              LOGF(info, "Particle with index %d and pdg code %d assigned to MC collision %d, pT: %f, phi: %f, eta: %f",
                   particle.globalIndex(), particle.pdgCode(), particle.mcCollisionId(), particle.pt(), particle.phi(), particle.eta());
              LOGF(info, "With status %d and flags %0X and multi-reconstructed as: ==================================", particle.statusCode(), particle.flags());
              for (int i = 0; i < mclabelpos[collsign][ixpart].size(); ++i) {
                auto track = tracks.iteratorAt(mclabelpos[collsign][ixpart][i]);
                auto coll = colls.iteratorAt(track.collisionId());
                LOGF(info, "Track with index %d and label %d assigned to collision %d, with associated MC collision %d",
                     track.globalIndex(), ixpart, track.collisionId(), coll.mcCollisionId());
              }
              LOGF(info, "END multi-reconstructed:   ==================================================================");
            }
            histos.get<TH1>(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("pdgcodemr"))->Fill(TString::Format("%d", particle.pdgCode()).Data(), 1.0);
          }
        }

        for (int i = 0; i < mclabelpos[collsign][ixpart].size(); ++i) {
          auto track1 = tracks.iteratorAt(mclabelpos[collsign][ixpart][i]);
          for (int j = i + 1; j < mclabelpos[collsign][ixpart].size(); ++j) {
            auto track2 = tracks.iteratorAt(mclabelpos[collsign][ixpart][j]);

            float deltaeta = track1.eta() - track2.eta();
            float deltaphi = track1.phi() - track2.phi();
            if (deltaphi < 0) {
              deltaphi += TWOPI;
            }
            if (deltaphi > TWOPI) {
              deltaphi -= TWOPI;
            }
            float deltapt = (track1.pt() > track2.pt()) ? track1.pt() - track2.pt() : track2.pt() - track1.pt();

            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("mrDeltaEta"), deltaeta);
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("mrDeltaPhi"), deltaphi);
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("mrDeltaPt"), deltapt);
          }
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("recomreta"), track1.eta());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("recomrphi"), track1.phi());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("recomrpt"), track1.pt());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("detectormapmr"), track1.detectorMap());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("dcaxymr"), track1.dcaXY());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("dcazmr"), track1.dcaZ());
          if (track1.dcaXY() < 1.0) {
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("finedcaxymr"), track1.dcaXY());
          }
          if (track1.dcaZ() < 1.0) {
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("finedcazmr"), track1.dcaZ());
          }
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecomreta"), track1.eta(), particle.eta());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecomrphi"), track1.phi(), particle.phi());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecomrpt"), track1.pt(), particle.pt());
          if (particle.mcCollisionId() != colls.iteratorAt(track1.collisionId()).mcCollisionId()) {
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("matchcollidmr"), kDONTMATCH + 0.5f);
          } else {
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("matchcollidmr"), kMATCH + 0.5f);
          }
        }
      } else if (nrec > 0) {
        auto track = tracks.iteratorAt(mclabelpos[collsign][ixpart][0]);
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecoeta"), track.eta(), particle.eta());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecophi"), track.phi(), particle.phi());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecopt"), track.pt(), particle.pt());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("detectormap"), track.detectorMap());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("dcaxy"), track.dcaXY());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("dcaz"), track.dcaZ());
        if (track.dcaXY() < 1.0) {
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("finedcaxy"), track.dcaXY());
        }
        if (track.dcaZ() < 1.0) {
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("finedcaz"), track.dcaZ());
        }
        if (particle.mcCollisionId() != colls.iteratorAt(track.collisionId()).mcCollisionId()) {
          if ((ba == kAFTER) and (collsign == kPOSITIVE) and cfgTrackCollAssoc) {
            LOGF(info, "Particle with index %d and pdg code %d assigned to MC collision %d, pT: %f, phi: %f, eta: %f",
                 particle.globalIndex(), particle.pdgCode(), particle.mcCollisionId(), particle.pt(), particle.phi(), particle.eta());
            LOGF(info, "        with status %d and flags %0X and", particle.statusCode(), particle.flags());
            LOGF(info, "        associated to track with index %d and label %d assigned to collision %d, with associated MC collision %d",
                 track.globalIndex(), ixpart, track.collisionId(), colls.iteratorAt(track.collisionId()).mcCollisionId());
          }
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("matchcollid"), kDONTMATCH + 0.5f);
        } else {
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("matchcollid"), kMATCH + 0.5f);
        }
      }
    }

    if (collsign == kPOSITIVE) {
      LOGF(info, "Reconstructed tracks (%s) with positive collision ID: %d with positive label, %d with negative label, %d with cross collision",
           ba == kAFTER ? "after" : "before", nrec_poslabel, nrec_neglabel, nrec_poslabel_crosscoll);
    } else {
      LOGF(info, "Reconstructed tracks (%s) with negative collision ID: %d with positive label, %d with negative label",
           ba == kAFTER ? "after" : "before", nrec_poslabel, nrec_neglabel);
    }
  }

  void
    processMapChecksBeforeCuts(soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::McTrackLabels> const& tracks,
                               soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
                               aod::McParticles const& mcParticles)
  {
    using namespace recogenmap;

    for (int i = 0; i < 2; ++i) {
      mclabelpos[i].clear();
      mclabelneg[i].clear();
      mclabelpos[i].resize(mcParticles.size());
      mclabelneg[i].resize(mcParticles.size());
    }

    size_t nreco = tracks.size();
    size_t ngen = 0;

    for (auto& part : mcParticles) {
      auto pdgpart = fPDG->GetParticle(part.pdgCode());
      if (pdgpart != nullptr) {
        float charge = (pdgpart->Charge() >= 3) ? 1.0 : ((pdgpart->Charge() <= -3) ? -1.0 : 0.0);
        if (charge != 0.0) {
          ngen++;
        }
      }
    }

    // Let's go through the reco-gen mapping to detect multi-reconstructed particles
    // For the time being we are only interested in the information based on the reconstructed tracks
    LOGF(info, "New dataframe (DF) with %d generated charged particles and %d reconstructed tracks", ngen, nreco);

    for (auto& track : tracks) {
      int64_t recix = track.globalIndex();
      int32_t label = track.mcParticleId();

      LOGF(MATCHRECGENLOGTRACKS, "Track with global Id %d and collision Id %d has label %d associated to MC collision %d", recix, track.collisionId(), label, track.mcParticle().mcCollisionId());
      if (track.collisionId() < 0) {
        if (label >= 0) {
          mclabelpos[kNEGATIVE][label].push_back(recix);
        } else {
          mclabelneg[kNEGATIVE][-label].push_back(recix);
        }
      } else {
        if (label >= 0) {
          mclabelpos[kPOSITIVE][label].push_back(recix);
        } else {
          mclabelneg[kPOSITIVE][-label].push_back(recix);
        }
      }
    }

    collectData<kBEFORE, kPOSITIVE>(tracks, mcParticles, collisions);
    collectData<kBEFORE, kNEGATIVE>(tracks, mcParticles, collisions);
  }
  PROCESS_SWITCH(CheckGeneratorLevelVsDetectorLevel, processMapChecksBeforeCuts, "Process detector <=> generator levels mapping checks before selection cuts", false);

  void processMapChecksCutsWithCent(soa::Join<aod::FullTracks, aod::TracksExtended, aod::TrackSelection, aod::McTrackLabels> const& tracks,
                                    soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels> const& collisions,
                                    aod::McParticles const& mcParticles)
  {
    using namespace recogenmap;

    for (int i = 0; i < 2; ++i) {
      mclabelpos[i].clear();
      mclabelneg[i].clear();
      mclabelpos[i].resize(mcParticles.size());
      mclabelneg[i].resize(mcParticles.size());
    }

    size_t nreco = 0;
    size_t ngen = 0;

    for (auto& part : mcParticles) {
      auto pdgpart = fPDG->GetParticle(part.pdgCode());
      if (pdgpart != nullptr) {
        float charge = (pdgpart->Charge() >= 3) ? 1.0 : ((pdgpart->Charge() <= -3) ? -1.0 : 0.0);
        if (charge != 0.0) {
          ngen++;
        }
      }
    }

    // Let's go through the reco-gen mapping to detect multi-reconstructed particles
    for (auto& track : tracks) {
      int64_t recix = track.globalIndex();
      int32_t label = track.mcParticleId();
      if (not(track.collisionId() < 0)) {
        auto coll = collisions.iteratorAt(track.collisionId());
        float centormult = -100.0f;
        if (IsEvtSelected(coll, centormult)) {
          if (not(label < 0)) {
            bool asone = false;
            bool astwo = false;

            AcceptTrack(track, asone, astwo);
            if (asone or astwo) {
              /* the track has been accepted */
              nreco++;
              LOGF(MATCHRECGENLOGTRACKS, "Accepted track with global Id %d and collision Id %d has label %d associated to MC collision %d", recix, track.collisionId(), label, track.mcParticle().mcCollisionId());
              mclabelpos[kPOSITIVE][label].push_back(recix);
            }
          }
        }
      }
    }
    LOGF(info, "New dataframe (DF) with %d generated charged particles and %d reconstructed accepted tracks", ngen, nreco);

    collectData<kAFTER, kPOSITIVE>(tracks, mcParticles, collisions);
  }
  PROCESS_SWITCH(CheckGeneratorLevelVsDetectorLevel, processMapChecksCutsWithCent, "Process detector <=> generator levels mapping checks after selection cuts", false);

  void processMapChecksCutsWithoutCent(soa::Join<aod::FullTracks, aod::TracksExtended, aod::TrackSelection, aod::McTrackLabels> const& tracks,
                                       soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels> const& collisions,
                                       aod::McParticles const& mcParticles)
  {
    using namespace recogenmap;

    for (int i = 0; i < 2; ++i) {
      mclabelpos[i].clear();
      mclabelneg[i].clear();
      mclabelpos[i].resize(mcParticles.size());
      mclabelneg[i].resize(mcParticles.size());
    }

    size_t nreco = 0;
    size_t ngen = 0;

    for (auto& part : mcParticles) {
      auto pdgpart = fPDG->GetParticle(part.pdgCode());
      if (pdgpart != nullptr) {
        float charge = (pdgpart->Charge() >= 3) ? 1.0 : ((pdgpart->Charge() <= -3) ? -1.0 : 0.0);
        if (charge != 0.0) {
          ngen++;
        }
      }
    }

    // Let's go through the reco-gen mapping to detect multi-reconstructed particles
    for (auto& track : tracks) {
      int64_t recix = track.globalIndex();
      int32_t label = track.mcParticleId();
      if (not(label < 0)) {
        if (not(track.collisionId() < 0)) {
          auto coll = collisions.iteratorAt(track.collisionId());
          float centormult = -100.0f;
          if (IsEvtSelectedNoCentMult(coll, centormult)) {
            bool asone = false;
            bool astwo = false;

            AcceptTrack(track, asone, astwo);
            if (asone or astwo) {
              /* the track has been accepted */
              nreco++;
              LOGF(MATCHRECGENLOGTRACKS, "Accepted track with global Id %d and collision Id %d has label %d associated to MC collision %d", recix, track.collisionId(), label, track.mcParticle().mcCollisionId());
              mclabelpos[kPOSITIVE][label].push_back(recix);
            }
          }
        }
      }
    }
    LOGF(info, "New dataframe (DF) with %d generated charged particles and %d reconstructed accepted tracks", ngen, nreco);

    collectData<kAFTER, kPOSITIVE>(tracks, mcParticles, collisions);
  }
  PROCESS_SWITCH(CheckGeneratorLevelVsDetectorLevel, processMapChecksCutsWithoutCent, "Process detector <=> generator levels mapping checks after selection cuts", false);

  void processDummy(aod::Collisions const& colls)
  {
    LOGF(MATCHRECGENLOGCOLLISIONS, "Got a new set of %d collisions", colls.size());
  }
  PROCESS_SWITCH(CheckGeneratorLevelVsDetectorLevel, processDummy, "Dummy process of detector <=> generator levels mapping checks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<DptDptCorrelationsFilterAnalysisTask>(cfgc, SetDefaultProcesses{{{"processWithoutCent", true}, {"processWithoutCentMC", true}}}),
    adaptAnalysisTask<TracksAndEventClassificationQARec>(cfgc),
    adaptAnalysisTask<TracksAndEventClassificationQAGen>(cfgc),
    adaptAnalysisTask<CheckGeneratorLevelVsDetectorLevel>(cfgc, SetDefaultProcesses{{{"processMapChecksCutsWithCent", false}}})};
  return workflow;
}
