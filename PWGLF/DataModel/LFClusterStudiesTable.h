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
//
// Author: Giorgio Alberto Lucia

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGLF_DATAMODEL_LFCLUSTERSTUDIESTABLE_H_
#define PWGLF_DATAMODEL_LFCLUSTERSTUDIESTABLE_H_

namespace o2::aod
{

namespace LFClusterStudiesTables
{
DECLARE_SOA_COLUMN(PMother, pMother, float);
DECLARE_SOA_COLUMN(PtMother, ptMother, float);
DECLARE_SOA_COLUMN(EtaMother, etaMother, float);
DECLARE_SOA_COLUMN(PhiMother, phiMother, float);
DECLARE_SOA_COLUMN(MassMother, massMother, float);
DECLARE_SOA_COLUMN(PdgCodeMother, pdgCodeMother, int);
DECLARE_SOA_COLUMN(RadiusMother, radiusMother, float);
DECLARE_SOA_COLUMN(DcaMotherPV, dcaMotherPV, float);
DECLARE_SOA_COLUMN(CosPAMother, cosPAMother, float);
DECLARE_SOA_COLUMN(AlphaAPMother, alphaAPMother, float);
DECLARE_SOA_COLUMN(QtAPMother, qtAPMother, float);
DECLARE_SOA_COLUMN(McPdgCodeMother, mcPdgCodeMother, int);

/**
 * PartID:
 * 0: e
 * 1: #pi
 * 2: K
 * 3: p
 * 4: d
 * 5: ^{3}He
 */
DECLARE_SOA_COLUMN(PartID, partID, uint8_t);
DECLARE_SOA_COLUMN(PartIDMc, partIDMc, int);
DECLARE_SOA_COLUMN(IsPositive, isPositive, bool);

DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(PTPC, pTPC, float);
DECLARE_SOA_COLUMN(PIDinTrk, pidInTrk, uint32_t);
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);
DECLARE_SOA_COLUMN(DcaToPV, dcaToPV, float);
DECLARE_SOA_COLUMN(ItsClusterSize, itsClusterSize, uint32_t);
DECLARE_SOA_COLUMN(TpcSignal, tpcSignal, float);
DECLARE_SOA_COLUMN(TpcNcls, tpcNcls, uint8_t);
DECLARE_SOA_COLUMN(TpcNSigma, tpcNSigma, float);
DECLARE_SOA_COLUMN(TofNSigma, tofNSigma, float);
DECLARE_SOA_COLUMN(TofMass, tofMass, float);
DECLARE_SOA_COLUMN(Chi2its, chi2its, float);
DECLARE_SOA_COLUMN(Chi2tpc, chi2tpc, float);
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);
DECLARE_SOA_COLUMN(McPdgCode, mcPdgCode, int);

DECLARE_SOA_COLUMN(RunNumber, runNumber, int);

} // namespace LFClusterStudiesTables

DECLARE_SOA_TABLE(
  ClStTable, "AOD", "CLSTTABLE",
  LFClusterStudiesTables::P,
  LFClusterStudiesTables::Eta,
  LFClusterStudiesTables::Phi,
  LFClusterStudiesTables::ItsClusterSize,
  LFClusterStudiesTables::PartID);

DECLARE_SOA_TABLE(
  ClStTableMc, "AOD", "CLSTTABLEMC",
  LFClusterStudiesTables::PartIDMc);

DECLARE_SOA_TABLE(
  ClStTableExtra, "AOD", "CLSTTABLEEXTRA",
  LFClusterStudiesTables::PTPC,
  LFClusterStudiesTables::PIDinTrk,
  LFClusterStudiesTables::TpcNSigma,
  LFClusterStudiesTables::TofNSigma,
  LFClusterStudiesTables::TofMass,
  LFClusterStudiesTables::CosPAMother,
  LFClusterStudiesTables::MassMother);

DECLARE_SOA_TABLE(
  ClStTableColl, "AOD", "CLSTTABLECOLL",
  LFClusterStudiesTables::RunNumber);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFCLUSTERSTUDIESTABLE_H_
