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
/// \file LFSlimNucleiTables.h
/// \brief Slim nuclei tables
///

#ifndef PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_
#define PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#include <Math/Vector4D.h>

namespace o2::aod
{
namespace lfv0he3
{
DECLARE_SOA_COLUMN(Z, z, float);
DECLARE_SOA_COLUMN(CentT0C, centT0C, float);
} // namespace lfv0he3
DECLARE_SOA_TABLE(LFEvents, "AOD", "LFEVENT", o2::soa::Index<>, lfv0he3::Z, lfv0he3::CentT0C);

namespace lfv0he3
{
DECLARE_SOA_INDEX_COLUMN(LFEvent, lfEvent); // Collision ID for the event
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(DCAxy, dcaXY, float);
DECLARE_SOA_COLUMN(DCAz, dcaZ, float);
DECLARE_SOA_COLUMN(TPCnCls, tpcNCls, int);
DECLARE_SOA_COLUMN(TPCnClsPID, tpcNClsPID, int);
DECLARE_SOA_COLUMN(ITSClusterSizes, itsClusterSizes, uint32_t);
DECLARE_SOA_COLUMN(NsigmaTPCPion, nSigmaTPCPion, float);
DECLARE_SOA_COLUMN(NsigmaTPCProton, nSigmaTPCProton, float);
DECLARE_SOA_COLUMN(NsigmaTPC, nSigmaTPC, float);
DECLARE_SOA_COLUMN(DCAdaughters, dcaDaughters, float);
DECLARE_SOA_COLUMN(DCAPVProton, dcaPVProton, float);
DECLARE_SOA_COLUMN(DCAPVPion, dcaPVPion, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
} // namespace lfv0he3
DECLARE_SOA_TABLE_VERSIONED(LFHe3_000, "AOD", "LFHE3V0", 0, lfv0he3::LFEventId, lfv0he3::Pt, lfv0he3::Eta, lfv0he3::Phi, lfv0he3::DCAxy, lfv0he3::DCAz, lfv0he3::TPCnCls, lfv0he3::ITSClusterSizes, lfv0he3::NsigmaTPC, lfv0he3::Sign);
DECLARE_SOA_TABLE_VERSIONED(LFLambda_000, "AOD", "LFLAMBDA", 0, lfv0he3::LFEventId, lfv0he3::Pt, lfv0he3::Eta, lfv0he3::Phi, lfv0he3::Mass, lfv0he3::CosPA, lfv0he3::DCAdaughters, lfv0he3::DCAPVProton, lfv0he3::DCAPVPion, lfv0he3::V0Radius, lfv0he3::Sign);

DECLARE_SOA_TABLE_VERSIONED(LFHe3_001, "AOD", "LFHE3V0", 1, lfv0he3::LFEventId, lfv0he3::Pt, lfv0he3::Eta, lfv0he3::Phi, lfv0he3::DCAxy, lfv0he3::DCAz, lfv0he3::TPCnCls, lfv0he3::TPCnClsPID, lfv0he3::ITSClusterSizes, lfv0he3::NsigmaTPC, lfv0he3::Sign);
DECLARE_SOA_TABLE_VERSIONED(LFLambda_001, "AOD", "LFLAMBDA", 1, lfv0he3::LFEventId, lfv0he3::Pt, lfv0he3::Eta, lfv0he3::Phi, lfv0he3::Mass, lfv0he3::CosPA, lfv0he3::DCAdaughters, lfv0he3::DCAPVProton, lfv0he3::DCAPVPion, lfv0he3::V0Radius, lfv0he3::NsigmaTPCProton, lfv0he3::NsigmaTPCPion, lfv0he3::Sign);
} // namespace o2::aod

struct he3Candidate {
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> momentum; // 4-momentum of the He3 candidate
  float nSigmaTPC = -999.f;                                            // TPC nSigma for He3
  float dcaXY = -999.f;
  float dcaZ = -999.f;
  int tpcNClsFound = 0;         // Number of TPC clusters found
  int tpcNClsPID = 0;           // Number of TPC clusters used for PID
  int itsNCls = 0;              // Number of ITS clusters
  uint32_t itsClusterSizes = 0; // ITS cluster sizes
  int8_t sign = 0;              // Charge sign of the He3 candidate
};

struct lambdaCandidate {
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> momentum;
  float mass = -1.f;              // Lambda mass
  float cosPA = -2.f;             // Cosine of pointing angle
  float dcaV0Daughters = -999.f;  // DCA between V0 daughters
  float dcaProtonToPV = -999.f;   // DCA of the proton to primary vertex
  float dcaPionToPV = -999.f;     // DCA of the pion to primary vertex
  float v0Radius = -1.f;          // V0 radius
  float protonNSigmaTPC = -999.f; // Proton TPC nSigma
  float pionNSigmaTPC = -999.f;   // Pion TPC nSigma
  int8_t sign = 0; // Charge sign of the Lambda candidate
};

#endif // PWGLF_DATAMODEL_LFSLIMNUCLEITABLES_H_
