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
/// \brief QA task for lambda polarization induced by jet analysis using derived data
///
/// \author Youpeng Su (yousu@cern.ch)

#ifndef PWGLF_DATAMODEL_LAMBDAJETPOLARIZATION_H_
#define PWGLF_DATAMODEL_LAMBDAJETPOLARIZATION_H_

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Boost.h"
#include "Math/Vector4D.h"
#include "TRandom.h"

namespace o2::aod
{
DECLARE_SOA_TABLE(MyCollisions, "AOD", "MYCOLLISION", //! vertex information of collision
                  o2::soa::Index<>, collision::PosZ);
using MyCollision = MyCollisions::iterator;

DECLARE_SOA_TABLE(MyCollisionsV0, "AOD", "MYCOLLISIONV0", //! vertex information of collision
                  o2::soa::Index<>, collision::PosX);
using MyCollisionV0s = MyCollisionsV0::iterator;

namespace myTable
{
DECLARE_SOA_INDEX_COLUMN(MyCollision, mycollision);
DECLARE_SOA_COLUMN(MyCollisionV0, mycollisionv0, Int_t);
DECLARE_SOA_COLUMN(V0px, v0px, Float_t);
DECLARE_SOA_COLUMN(V0py, v0py, Float_t);
DECLARE_SOA_COLUMN(V0pz, v0pz, Float_t);
DECLARE_SOA_COLUMN(V0pT, v0pt, Float_t);
DECLARE_SOA_COLUMN(V0Lambdamass, v0Lambdamass, Float_t);
DECLARE_SOA_COLUMN(V0protonpx, v0protonpx, Float_t);
DECLARE_SOA_COLUMN(V0protonpy, v0protonpy, Float_t);
DECLARE_SOA_COLUMN(V0protonpz, v0protonpz, Float_t);
DECLARE_SOA_COLUMN(MyCollisionJet, mycollisionjet, Int_t);
DECLARE_SOA_COLUMN(Jetpx, jetpx, Float_t);
DECLARE_SOA_COLUMN(Jetpy, jetpy, Float_t);
DECLARE_SOA_COLUMN(Jetpz, jetpz, Float_t);
DECLARE_SOA_COLUMN(JetpT, jetpt, Float_t);
DECLARE_SOA_COLUMN(MyCollisionLeadingJet, mycollisionleadingjet, Int_t);
DECLARE_SOA_COLUMN(LeadingJetpx, leadingjetpx, Float_t);
DECLARE_SOA_COLUMN(LeadingJetpy, leadingjetpy, Float_t);
DECLARE_SOA_COLUMN(LeadingJetpz, leadingjetpz, Float_t);
DECLARE_SOA_COLUMN(LeadingJetpT, leadingjetpt, Float_t);

} // namespace myTable

DECLARE_SOA_TABLE(MyTable, "AOD", "MYTABLE", o2::soa::Index<>,
                  myTable::MyCollisionId, myTable::MyCollisionV0, myTable::V0px, myTable::V0py, myTable::V0pz, myTable::V0pT, myTable::V0Lambdamass,
                  myTable::V0protonpx, myTable::V0protonpy, myTable::V0protonpz);

DECLARE_SOA_TABLE(MyTableAnti, "AOD", "MYTABLEAnti", o2::soa::Index<>,
                  myTable::MyCollisionId, myTable::MyCollisionV0, myTable::V0px, myTable::V0py, myTable::V0pz, myTable::V0pT, myTable::V0Lambdamass,
                  myTable::V0protonpx, myTable::V0protonpy, myTable::V0protonpz);

DECLARE_SOA_TABLE(MyTableJet, "AOD", "MYTABLEJet", o2::soa::Index<>,
                  myTable::MyCollisionId, myTable::MyCollisionJet, myTable::Jetpx, myTable::Jetpy, myTable::Jetpz, myTable::JetpT);

DECLARE_SOA_TABLE(MyTableLeadingJet, "AOD", "LeadingJet", o2::soa::Index<>, myTable::MyCollisionId, myTable::MyCollisionLeadingJet, myTable::LeadingJetpx, myTable::LeadingJetpy, myTable::LeadingJetpz, myTable::LeadingJetpT);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LAMBDAJETPOLARIZATION_H_
