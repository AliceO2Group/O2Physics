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
/// \brief QA task for Cascade analysis using derived data
///
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)

#ifndef PWGLF_DATAMODEL_CASCQAANALYSIS_H_
#define PWGLF_DATAMODEL_CASCQAANALYSIS_H_

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "TRandom.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace mycascades
{

enum EvFlags : uint8_t {
  EvINEL = 0x1,    // INEL Event
  EvINELgt0 = 0x2, // Event with at least 1 PV contributors from the |eta| < 1
  EvINELgt1 = 0x4  // Event with at least 2 PV contributors from the |eta| < 1
};

DECLARE_SOA_COLUMN(CollisionZ, zcoll, float);
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);
DECLARE_SOA_COLUMN(MultFV0A, multFV0A, float);
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(RapXi, rapxi, float);
DECLARE_SOA_COLUMN(RapOmega, rapomega, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(MassXi, massxi, float);
DECLARE_SOA_COLUMN(MassOmega, massomega, float);
DECLARE_SOA_COLUMN(MassLambdaDau, masslambdadau, float);
DECLARE_SOA_COLUMN(CascRadius, cascradius, float);
DECLARE_SOA_COLUMN(V0Radius, v0radius, float);
DECLARE_SOA_COLUMN(CascCosPA, casccospa, float);
DECLARE_SOA_COLUMN(V0CosPA, v0cospa, float);
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);
DECLARE_SOA_COLUMN(DCABachToPV, dcabachtopv, float);
DECLARE_SOA_COLUMN(DCACascDaughters, dcacascdaughters, float);
DECLARE_SOA_COLUMN(DCAV0Daughters, dcav0daughters, float);
DECLARE_SOA_COLUMN(DCAV0ToPV, dcav0topv, float);
DECLARE_SOA_COLUMN(PosEta, poseta, float);
DECLARE_SOA_COLUMN(NegEta, negeta, float);
DECLARE_SOA_COLUMN(BachEta, bacheta, float);
DECLARE_SOA_COLUMN(PosITSHits, positshits, float);
DECLARE_SOA_COLUMN(NegITSHits, negitshits, float);
DECLARE_SOA_COLUMN(BachITSHits, bachitshits, float);
DECLARE_SOA_COLUMN(CtauXi, ctauxi, float);
DECLARE_SOA_COLUMN(CtauOmega, ctauomega, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPr, ntpcsigmanegpr, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPr, ntpcsigmapospr, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPi, ntpcsigmanegpi, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPi, ntpcsigmapospi, float);
DECLARE_SOA_COLUMN(NTPCSigmaBachPi, ntpcsigmabachpi, float);
DECLARE_SOA_COLUMN(NTPCSigmaBachKa, ntpcsigmabachka, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPr, ntofsigmanegpr, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPr, ntofsigmapospr, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPi, ntofsigmanegpi, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPi, ntofsigmapospi, float);
DECLARE_SOA_COLUMN(NTOFSigmaBachPi, ntofsigmabachpi, float);
DECLARE_SOA_COLUMN(NTOFSigmaBachKa, ntofsigmabachka, float);
DECLARE_SOA_COLUMN(PosNTPCClusters, posntpcscls, float);
DECLARE_SOA_COLUMN(NegNTPCClusters, negntpcscls, float);
DECLARE_SOA_COLUMN(BachNTPCClusters, bachntpcscls, float);
DECLARE_SOA_COLUMN(PosNTPCCrossedRows, posntpccrrows, float);
DECLARE_SOA_COLUMN(NegNTPCCrossedRows, negntpccrrows, float);
DECLARE_SOA_COLUMN(BachNTPCCrossedRows, bachntpccrrows, float);
DECLARE_SOA_COLUMN(PosHasTOF, poshastof, float);
DECLARE_SOA_COLUMN(NegHasTOF, neghastof, float);
DECLARE_SOA_COLUMN(BachHasTOF, bachhastof, float);
DECLARE_SOA_COLUMN(PosPt, pospt, float);
DECLARE_SOA_COLUMN(NegPt, negpt, float);
DECLARE_SOA_COLUMN(BachPt, bachpt, float);
DECLARE_SOA_COLUMN(McPdgCode, mcPdgCode, float);                     //! -1 unknown
DECLARE_SOA_COLUMN(IsPrimary, isPrimary, float);                     //! -1 unknown, 0 not primary, 1 primary
DECLARE_SOA_COLUMN(BachBaryonCosPA, bachBaryonCosPA, float);         //! avoid bach-baryon correlated inv mass structure in analysis
DECLARE_SOA_COLUMN(BachBaryonDCAxyToPV, bachBaryonDCAxyToPV, float); //! avoid bach-baryon correlated inv mass structure in analysis
DECLARE_SOA_COLUMN(EventSelFilterBitMask, eventSelFilterBitMask, uint8_t);

DECLARE_SOA_DYNAMIC_COLUMN(IsINEL, isINEL, //! True if the Event belongs to the INEL event class
                           [](uint8_t flags) -> bool { return (flags & EvFlags::EvINEL) == EvFlags::EvINEL; });
DECLARE_SOA_DYNAMIC_COLUMN(IsINELgt0, isINELgt0, //! True if the Event belongs to the INELgt0 event class
                           [](uint8_t flags) -> bool { return (flags & EvFlags::EvINELgt0) == EvFlags::EvINELgt0; });
DECLARE_SOA_DYNAMIC_COLUMN(IsINELgt1, isINELgt1, //! True if the Event belongs to the INELgt1 event class
                           [](uint8_t flags) -> bool { return (flags & EvFlags::EvINELgt1) == EvFlags::EvINELgt1; });
} // namespace mycascades

DECLARE_SOA_TABLE(MyCascades, "AOD", "MYCASCADES", o2::soa::Index<>,
                  mycascades::CollisionZ, mycascades::MultFT0M, mycascades::MultFV0A, mycascades::Sign, mycascades::Pt, mycascades::RapXi, mycascades::RapOmega, mycascades::Eta, mycascades::MassXi, mycascades::MassOmega, mycascades::MassLambdaDau, mycascades::CascRadius, mycascades::V0Radius,
                  mycascades::CascCosPA, mycascades::V0CosPA, mycascades::DCAPosToPV, mycascades::DCANegToPV,
                  mycascades::DCABachToPV, mycascades::DCACascDaughters, mycascades::DCAV0Daughters, mycascades::DCAV0ToPV, mycascades::PosEta, mycascades::NegEta,
                  mycascades::BachEta, mycascades::PosITSHits, mycascades::NegITSHits, mycascades::BachITSHits,
                  mycascades::CtauXi, mycascades::CtauOmega,
                  mycascades::NTPCSigmaNegPr, mycascades::NTPCSigmaPosPr, mycascades::NTPCSigmaNegPi, mycascades::NTPCSigmaPosPi, mycascades::NTPCSigmaBachPi, mycascades::NTPCSigmaBachKa,
                  mycascades::NTOFSigmaNegPr, mycascades::NTOFSigmaPosPr, mycascades::NTOFSigmaNegPi,
                  mycascades::NTOFSigmaPosPi, mycascades::NTOFSigmaBachPi, mycascades::NTOFSigmaBachKa,
                  mycascades::PosNTPCClusters, mycascades::NegNTPCClusters, mycascades::BachNTPCClusters,
                  mycascades::PosNTPCCrossedRows, mycascades::NegNTPCCrossedRows, mycascades::BachNTPCCrossedRows,
                  mycascades::PosHasTOF, mycascades::NegHasTOF, mycascades::BachHasTOF,
                  mycascades::PosPt, mycascades::NegPt, mycascades::BachPt,
                  mycascades::McPdgCode, mycascades::IsPrimary,
                  mycascades::BachBaryonCosPA, mycascades::BachBaryonDCAxyToPV,
                  mycascades::EventSelFilterBitMask,
                  mycascades::IsINEL<mycascades::EventSelFilterBitMask>,
                  mycascades::IsINELgt0<mycascades::EventSelFilterBitMask>,
                  mycascades::IsINELgt1<mycascades::EventSelFilterBitMask>);

DECLARE_SOA_TABLE(CascTraining, "AOD", "CascTraining", o2::soa::Index<>,
                  mycascades::MultFT0M, mycascades::Sign, mycascades::Pt, mycascades::Eta, mycascades::MassXi, mycascades::MassOmega, mycascades::MassLambdaDau, mycascades::CascRadius, mycascades::V0Radius, mycascades::CascCosPA, mycascades::V0CosPA, mycascades::DCAPosToPV, mycascades::DCANegToPV,
                  mycascades::DCABachToPV, mycascades::DCACascDaughters, mycascades::DCAV0Daughters, mycascades::DCAV0ToPV, mycascades::BachBaryonCosPA, mycascades::BachBaryonDCAxyToPV, mycascades::McPdgCode);

namespace myMCcascades
{

enum EvFlags : uint8_t {
  EvINEL = 0x1,    // INEL Event
  EvINELgt0 = 0x2, // Event with at least 1 PV contributors from the |eta| < 1
  EvINELgt1 = 0x4  // Event with at least 2 PV contributors from the |eta| < 1
};

DECLARE_SOA_COLUMN(CollisionZ, zcoll, float);
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(IsPrimary, isPrimary, bool);
DECLARE_SOA_COLUMN(NAssocColl, nAssocColl, int); // Number of reconstructed collisions assoceated to the generated one of this cascade
DECLARE_SOA_COLUMN(NChInFT0M, nChInFT0M, float); // Number of charged particles in FT0M acceptance
DECLARE_SOA_COLUMN(AssCollisionTypeFilterBitMask, assCollisionTypeFilterBitMask, uint8_t);
DECLARE_SOA_COLUMN(McCollisionTypeFilterBitMask, mcCollisionTypeFilterBitMask, uint8_t);

DECLARE_SOA_DYNAMIC_COLUMN(IsINELassoc, isINELassoc, //! True if there's at least 1 reconstructed INEL event for the generated one of this cascade
                           [](uint8_t flags) -> bool { return (flags & EvFlags::EvINEL) == EvFlags::EvINEL; });
DECLARE_SOA_DYNAMIC_COLUMN(IsINELgt0assoc, isINELgt0assoc, //! True if there's at least 1 reconstructed INEL>0 event for the generated one of this cascade
                           [](uint8_t flags) -> bool { return (flags & EvFlags::EvINELgt0) == EvFlags::EvINELgt0; });
DECLARE_SOA_DYNAMIC_COLUMN(IsINELgt1assoc, isINELgt1assoc, //! True if there's at least 1 reconstructed INEL>1 event for the generated one of this cascade
                           [](uint8_t flags) -> bool { return (flags & EvFlags::EvINELgt1) == EvFlags::EvINELgt1; });

DECLARE_SOA_DYNAMIC_COLUMN(IsINEL, isINEL, //! True if the Event belongs to the INEL event class
                           [](uint8_t flags) -> bool { return (flags & EvFlags::EvINEL) == EvFlags::EvINEL; });
DECLARE_SOA_DYNAMIC_COLUMN(IsINELgt0, isINELgt0, //! True if the Event belongs to the INELgt0 event class
                           [](uint8_t flags) -> bool { return (flags & EvFlags::EvINELgt0) == EvFlags::EvINELgt0; });
DECLARE_SOA_DYNAMIC_COLUMN(IsINELgt1, isINELgt1, //! True if the Event belongs to the INELgt1 event class
                           [](uint8_t flags) -> bool { return (flags & EvFlags::EvINELgt1) == EvFlags::EvINELgt1; });
} // namespace myMCcascades

DECLARE_SOA_TABLE(MyMCCascades, "AOD", "MYMCCASCADES", o2::soa::Index<>,
                  myMCcascades::CollisionZ, myMCcascades::Sign, myMCcascades::PdgCode,
                  myMCcascades::Y, myMCcascades::Eta, myMCcascades::Phi, myMCcascades::Pt,
                  myMCcascades::IsPrimary, myMCcascades::NAssocColl, myMCcascades::NChInFT0M,
                  myMCcascades::AssCollisionTypeFilterBitMask,
                  myMCcascades::McCollisionTypeFilterBitMask,
                  myMCcascades::IsINELassoc<myMCcascades::AssCollisionTypeFilterBitMask>,
                  myMCcascades::IsINELgt0assoc<myMCcascades::AssCollisionTypeFilterBitMask>,
                  myMCcascades::IsINELgt1assoc<myMCcascades::AssCollisionTypeFilterBitMask>,
                  myMCcascades::IsINEL<myMCcascades::McCollisionTypeFilterBitMask>,
                  myMCcascades::IsINELgt0<myMCcascades::McCollisionTypeFilterBitMask>,
                  myMCcascades::IsINELgt1<myMCcascades::McCollisionTypeFilterBitMask>);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_CASCQAANALYSIS_H_
