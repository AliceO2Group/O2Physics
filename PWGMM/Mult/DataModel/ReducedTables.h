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

#ifndef PWGMM_MULT_DATAMODEL_REDUCEDTABLES_H_
#define PWGMM_MULT_DATAMODEL_REDUCEDTABLES_H_
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"

namespace o2::aod
{

#define BCcols o2::soa::Index<>, \
               bc::RunNumber

// Reduced BCs as a root index
DECLARE_SOA_TABLE(RBCs, "AOD", "RBC",
                  BCcols,
                  soa::Marker<1>);
DECLARE_SOA_TABLE(StoredRBCs, "AOD1", "RBC",
                  BCcols,
                  soa::Marker<2>);

namespace rcol
{
DECLARE_SOA_INDEX_COLUMN(RBC, rbc);
}

#define Ccols o2::soa::Index<>,            \
              rcol::RBCId,                 \
              collision::PosX,             \
              collision::PosY,             \
              collision::PosZ,             \
              collision::CollisionTimeRes, \
              mult::MultFT0A,              \
              mult::MultFT0C,              \
              mult::MultFDDA,              \
              mult::MultFDDC,              \
              mult::MultZNA,               \
              mult::MultZNC,               \
              mult::MultNTracksPV,         \
              mult::MultNTracksPVeta1,     \
              mult::MultNTracksPVetaHalf

#define CCcols cent::CentFV0A, \
               cent::CentFT0M, \
               cent::CentFT0A, \
               cent::CentFT0C, \
               cent::CentFDDM, \
               cent::CentNTPV

// Reduced Collisions
DECLARE_SOA_TABLE(RCollisions, "AOD", "RCOLLISION",
                  Ccols,
                  soa::Marker<1>)
DECLARE_SOA_TABLE(StoredRCollisions, "AOD1", "RCOLLISION",
                  Ccols,
                  soa::Marker<2>)

DECLARE_SOA_TABLE(RCents, "AOD", "RCENTS",
                  CCcols,
                  soa::Marker<1>)
DECLARE_SOA_TABLE(StoredRCents, "AOD1", "RCENTS",
                  CCcols,
                  soa::Marker<2>)

// Reduced tracks (is this needed?)
namespace rtrack
{
DECLARE_SOA_INDEX_COLUMN(RCollision, rcollision);
DECLARE_SOA_COLUMN(Weight, weight, float);
} // namespace rtrack
#define Tcols o2::soa::Index<>,     \
              rtrack::RCollisionId, \
              rtrack::Weight,       \
              track::Pt,            \
              track::P,             \
              track::Eta,           \
              track::Phi,           \
              track::DcaXY,         \
              track::DcaZ

DECLARE_SOA_TABLE(RTracks, "AOD", "RTRACK",
                  Tcols,
                  soa::Marker<1>)
DECLARE_SOA_TABLE(StoredRTracks, "AOD1", "RTRACK",
                  Tcols,
                  soa::Marker<2>)

#define TFcols o2::soa::Index<>,     \
               rtrack::RCollisionId, \
               rtrack::Weight,       \
               fwdtrack::Pt,         \
               fwdtrack::P,          \
               fwdtrack::Eta,        \
               fwdtrack::Phi,        \
               fwdtrack::FwdDcaX,    \
               fwdtrack::FwdDcaY

DECLARE_SOA_TABLE(RFTracks, "AOD", "RFTRACK",
                  TFcols,
                  soa::Marker<1>)
DECLARE_SOA_TABLE(StoredRFTracks, "AOD1", "RFTRACK",
                  TFcols,
                  soa::Marker<2>)

// Reduced MC collisions
namespace rmccol
{
DECLARE_SOA_COLUMN(Weight, weight, float);
}
#define MCCcols o2::soa::Index<>,             \
                rcol::RBCId,                  \
                rmccol::Weight,               \
                mccollision::PosX,            \
                mccollision::PosY,            \
                mccollision::PosZ,            \
                mccollision::ImpactParameter, \
                mult::MultMCFT0A,             \
                mult::MultMCFT0C,             \
                mult::MultMCNParticlesEta05,  \
                mult::MultMCNParticlesEta10

DECLARE_SOA_TABLE(RMCCollisions, "AOD", "RMCCOLLISION",
                  MCCcols,
                  soa::Marker<1>)
DECLARE_SOA_TABLE(StoredRMCCollisions, "AOD1", "RMCCOLLISION",
                  MCCcols,
                  soa::Marker<2>)

// Extra MC tables
namespace rhepmc
{
DECLARE_SOA_INDEX_COLUMN(RMCCollision, rmccollison);
}
#define HMCcols rhepmc::RMCCollisionId,   \
                hepmcxsection::XsectGen,  \
                hepmcxsection::PtHard,    \
                hepmcxsection::NMPI,      \
                hepmcxsection::ProcessId, \
                hepmcpdfinfo::Id1,        \
                hepmcpdfinfo::Id2,        \
                hepmcpdfinfo::PdfId1,     \
                hepmcpdfinfo::PdfId2,     \
                hepmcpdfinfo::X1,         \
                hepmcpdfinfo::X2,         \
                hepmcpdfinfo::ScalePdf,   \
                hepmcpdfinfo::Pdf1,       \
                hepmcpdfinfo::Pdf2

DECLARE_SOA_TABLE(RHepMCinfos, "AOD", "RHEPMCINFO",
                  HMCcols,
                  soa::Marker<1>);
DECLARE_SOA_TABLE(StoredRHepMCinfos, "AOD1", "RHEPMCINFO",
                  HMCcols,
                  soa::Marker<2>);

#define HMCHIcols rhepmc::RMCCollisionId,         \
                  hepmcheavyion::NcollHard,       \
                  hepmcheavyion::NpartProj,       \
                  hepmcheavyion::NpartTarg,       \
                  hepmcheavyion::Ncoll,           \
                  hepmcheavyion::ImpactParameter, \
                  hepmcheavyion::EventPlaneAngle, \
                  hepmcheavyion::SigmaInelNN,     \
                  hepmcheavyion::Centrality

DECLARE_SOA_TABLE(RHepMCHIs, "AOD", "RHEPMCHI",
                  HMCHIcols,
                  soa::Marker<1>);
DECLARE_SOA_TABLE(StoredRHepMCHIs, "AOD1", "RHEPMCHI",
                  HMCHIcols,
                  soa::Marker<2>);

namespace rparticle
{
DECLARE_SOA_INDEX_COLUMN(RMCCollision, rmccollision);
}

// Reduced MC particles (is this needed?)
#define RMCPcols o2::soa::Index<>,             \
                 rparticle::RMCCollisionId,    \
                 mcparticle::PdgCode,          \
                 mcparticle::MothersIds,       \
                 mcparticle::DaughtersIdSlice, \
                 mcparticle::Pt,               \
                 mcparticle::P,                \
                 mcparticle::Eta,              \
                 mcparticle::Y,                \
                 mcparticle::Phi,              \
                 mcparticle::E,                \
                 mcparticle::Weight

DECLARE_SOA_TABLE(RMCParticles, "AOD", "RMCPARTICLE",
                  RMCPcols,
                  soa::Marker<1>)
DECLARE_SOA_TABLE(StoredRMCParticles, "AOD1", "RMCPARTICLE",
                  RMCPcols,
                  soa::Marker<2>)

// label tables
namespace rlabels
{
DECLARE_SOA_INDEX_COLUMN(RMCCollision, rmccollision);
DECLARE_SOA_INDEX_COLUMN(RMCParticle, rmcparticle);
} // namespace rlabels
DECLARE_SOA_TABLE(RMCTrackLabels, "AOD", "RMCTRKLABEL",
                  rlabels::RMCParticleId, soa::Marker<1>)
DECLARE_SOA_TABLE(StoredRMCTrackLabels, "AOD1", "RMCTRKLABEL",
                  rlabels::RMCParticleId, soa::Marker<2>)

DECLARE_SOA_TABLE(RMCColLabels, "AOD", "RMCCOLLABEL",
                  rlabels::RMCCollisionId, soa::Marker<1>)
DECLARE_SOA_TABLE(StoredRMCColLabels, "AOD1", "RMCCOLLABEL",
                  rlabels::RMCCollisionId, soa::Marker<2>)

namespace features
{
DECLARE_SOA_COLUMN(GeneratedCentralMultiplicity, t, int);
DECLARE_SOA_COLUMN(ReconstructedCentralMultiplicity, m, int);
DECLARE_SOA_COLUMN(GeneratedVertexX, vtX, float);
DECLARE_SOA_COLUMN(GeneratedVertexY, vtY, float);
DECLARE_SOA_COLUMN(GeneratedVertexZ, vtZ, float);
DECLARE_SOA_COLUMN(ReconstructedVertexX, vmX, float);
DECLARE_SOA_COLUMN(ReconstructedVertexY, vmY, float);
DECLARE_SOA_COLUMN(ReconstructedVertexZ, vmZ, float);
DECLARE_SOA_COLUMN(TimeRes, tres, float);
DECLARE_SOA_COLUMN(GeneratedForwardMultiplicityA, tfA, int);
DECLARE_SOA_COLUMN(GeneratedForwardMultiplicityC, tfC, int);
DECLARE_SOA_COLUMN(ReconstructedForwardMultiplicityA, mfA, float);
DECLARE_SOA_COLUMN(ReconstructedForwardMultiplicityC, mfC, float);
DECLARE_SOA_COLUMN(ImpactParameter, b, float);
DECLARE_SOA_COLUMN(ProcessId, id, int);
} // namespace features

DECLARE_SOA_TABLE(RFeatMins, "AOD", "RFEATMIN",
                  soa::Index<>,
                  features::GeneratedCentralMultiplicity,
                  features::ReconstructedCentralMultiplicity,
                  features::ProcessId,
                  features::ImpactParameter,
                  features::ReconstructedVertexX,
                  features::ReconstructedVertexY,
                  features::ReconstructedVertexZ,
                  features::TimeRes,
                  features::ReconstructedForwardMultiplicityA,
                  features::ReconstructedForwardMultiplicityC);

} // namespace o2::aod
namespace o2::soa
{
DECLARE_EQUIVALENT_FOR_INDEX(aod::RBCs, aod::StoredRBCs);
DECLARE_EQUIVALENT_FOR_INDEX(aod::RCollisions, aod::StoredRCollisions);
DECLARE_EQUIVALENT_FOR_INDEX(aod::RMCCollisions, aod::StoredRMCCollisions);
DECLARE_EQUIVALENT_FOR_INDEX(aod::RMCParticles, aod::StoredRMCParticles);
} // namespace o2::soa
#endif // PWGMM_MULT_DATAMODEL_REDUCEDTABLES_H_
