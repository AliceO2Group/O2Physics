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

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

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

DECLARE_SOA_COLUMN(PDaughterA, pDaughterA, float);
DECLARE_SOA_COLUMN(PtDaughterA, ptDaughterA, float);
DECLARE_SOA_COLUMN(EtaDaughterA, etaDaughterA, float);
DECLARE_SOA_COLUMN(PhiDaughterA, phiDaughterA, float);
DECLARE_SOA_COLUMN(PTPCDaughterA, pTPCDaughterA, float);
DECLARE_SOA_COLUMN(PDGCodeDaughterA, pdgCodeDaughterA, int);
DECLARE_SOA_COLUMN(DcaToPVDaughterA, dcaToPVDaughterA, float);
DECLARE_SOA_COLUMN(ItsClusterSizeDaughterA, itsClusterSizeDaughterA, uint32_t);
DECLARE_SOA_COLUMN(TpcSignalDaughterA, tpcSignalDaughterA, float);
DECLARE_SOA_COLUMN(TpcNclsDaughterA, tpcNclsDaughterA, uint8_t);
DECLARE_SOA_COLUMN(TpcNsigmaElDaughterA, tpcNsigmaElDaughterA, float);
DECLARE_SOA_COLUMN(TpcNsigmaPiDaughterA, tpcNsigmaPiDaughterA, float);
DECLARE_SOA_COLUMN(TpcNsigmaKaDaughterA, tpcNsigmaKaDaughterA, float);
DECLARE_SOA_COLUMN(TpcNsigmaPrDaughterA, tpcNsigmaPrDaughterA, float);
DECLARE_SOA_COLUMN(TpcNsigmaDeDaughterA, tpcNsigmaDeDaughterA, float);
DECLARE_SOA_COLUMN(TpcNsigmaHeDaughterA, tpcNsigmaHeDaughterA, float);
DECLARE_SOA_COLUMN(TofNsigmaElDaughterA, tofNsigmaElDaughterA, float);
DECLARE_SOA_COLUMN(TofNsigmaPiDaughterA, tofNsigmaPiDaughterA, float);
DECLARE_SOA_COLUMN(TofNsigmaKaDaughterA, tofNsigmaKaDaughterA, float);
DECLARE_SOA_COLUMN(TofNsigmaPrDaughterA, tofNsigmaPrDaughterA, float);
DECLARE_SOA_COLUMN(TofNsigmaDeDaughterA, tofNsigmaDeDaughterA, float);
DECLARE_SOA_COLUMN(TofNsigmaHeDaughterA, tofNsigmaHeDaughterA, float);
DECLARE_SOA_COLUMN(Chi2itsDaughterA, chi2itsDaughterA, float);
DECLARE_SOA_COLUMN(Chi2tpcDaughterA, chi2tpcDaughterA, float);
DECLARE_SOA_COLUMN(HasTPCDaughterA, hasTPCDaughterA, bool);

DECLARE_SOA_COLUMN(PDaughterB, pDaughterB, float);
DECLARE_SOA_COLUMN(PtDaughterB, ptDaughterB, float);
DECLARE_SOA_COLUMN(EtaDaughterB, etaDaughterB, float);
DECLARE_SOA_COLUMN(PhiDaughterB, phiDaughterB, float);
DECLARE_SOA_COLUMN(PTPCDaughterB, pTPCDaughterB, float);
DECLARE_SOA_COLUMN(PDGCodeDaughterB, pdgCodeDaughterB, int);
DECLARE_SOA_COLUMN(DcaToPVDaughterB, dcaToPVDaughterB, float);
DECLARE_SOA_COLUMN(ItsClusterSizeDaughterB, itsClusterSizeDaughterB, uint32_t);
DECLARE_SOA_COLUMN(TpcSignalDaughterB, tpcSignalDaughterB, float);
DECLARE_SOA_COLUMN(TpcNclsDaughterB, tpcNclsDaughterB, uint8_t);
DECLARE_SOA_COLUMN(TpcNsigmaElDaughterB, tpcNsigmaElDaughterB, float);
DECLARE_SOA_COLUMN(TpcNsigmaPiDaughterB, tpcNsigmaPiDaughterB, float);
DECLARE_SOA_COLUMN(TpcNsigmaKaDaughterB, tpcNsigmaKaDaughterB, float);
DECLARE_SOA_COLUMN(TpcNsigmaPrDaughterB, tpcNsigmaPrDaughterB, float);
DECLARE_SOA_COLUMN(TpcNsigmaDeDaughterB, tpcNsigmaDeDaughterB, float);
DECLARE_SOA_COLUMN(TpcNsigmaHeDaughterB, tpcNsigmaHeDaughterB, float);
DECLARE_SOA_COLUMN(TofNsigmaElDaughterB, tofNsigmaElDaughterB, float);
DECLARE_SOA_COLUMN(TofNsigmaPiDaughterB, tofNsigmaPiDaughterB, float);
DECLARE_SOA_COLUMN(TofNsigmaKaDaughterB, tofNsigmaKaDaughterB, float);
DECLARE_SOA_COLUMN(TofNsigmaPrDaughterB, tofNsigmaPrDaughterB, float);
DECLARE_SOA_COLUMN(TofNsigmaDeDaughterB, tofNsigmaDeDaughterB, float);
DECLARE_SOA_COLUMN(TofNsigmaHeDaughterB, tofNsigmaHeDaughterB, float);
DECLARE_SOA_COLUMN(Chi2itsDaughterB, chi2itsDaughterB, float);
DECLARE_SOA_COLUMN(Chi2tpcDaughterB, chi2tpcDaughterB, float);
DECLARE_SOA_COLUMN(HasTPCDaughterB, hasTPCDaughterB, bool);
} // namespace LFClusterStudiesTables

DECLARE_SOA_TABLE(
  ClStV0Table, "AOD", "CLSTV0TABLE",
  LFClusterStudiesTables::PMother,
  // LFClusterStudiesTables::PtMother,
  LFClusterStudiesTables::EtaMother,
  LFClusterStudiesTables::PhiMother,
  // LFClusterStudiesTables::MassMother,
  LFClusterStudiesTables::PdgCodeMother,
  LFClusterStudiesTables::RadiusMother,
  LFClusterStudiesTables::DcaMotherPV,
  LFClusterStudiesTables::CosPAMother,
  LFClusterStudiesTables::AlphaAPMother,
  LFClusterStudiesTables::QtAPMother,
  LFClusterStudiesTables::PDaughterA,
  // LFClusterStudiesTables::PtDaughterA,
  LFClusterStudiesTables::EtaDaughterA,
  LFClusterStudiesTables::PhiDaughterA,
  LFClusterStudiesTables::PTPCDaughterA,
  // LFClusterStudiesTables::PDGCodeDaughterA,
  LFClusterStudiesTables::DcaToPVDaughterA,
  LFClusterStudiesTables::ItsClusterSizeDaughterA,
  LFClusterStudiesTables::TpcSignalDaughterA,
  LFClusterStudiesTables::TpcNclsDaughterA,
  LFClusterStudiesTables::TpcNsigmaElDaughterA,
  LFClusterStudiesTables::TpcNsigmaPiDaughterA,
  LFClusterStudiesTables::TpcNsigmaKaDaughterA,
  LFClusterStudiesTables::TpcNsigmaPrDaughterA,
  LFClusterStudiesTables::TpcNsigmaDeDaughterA,
  LFClusterStudiesTables::TpcNsigmaHeDaughterA,
  LFClusterStudiesTables::TofNsigmaElDaughterA,
  LFClusterStudiesTables::TofNsigmaPiDaughterA,
  LFClusterStudiesTables::TofNsigmaKaDaughterA,
  LFClusterStudiesTables::TofNsigmaPrDaughterA,
  LFClusterStudiesTables::TofNsigmaDeDaughterA,
  LFClusterStudiesTables::TofNsigmaHeDaughterA,
  LFClusterStudiesTables::Chi2itsDaughterA,
  LFClusterStudiesTables::Chi2tpcDaughterA,
  LFClusterStudiesTables::HasTPCDaughterA,
  LFClusterStudiesTables::PDaughterB,
  // LFClusterStudiesTables::PtDaughterB,
  LFClusterStudiesTables::EtaDaughterB,
  LFClusterStudiesTables::PhiDaughterB,
  LFClusterStudiesTables::PTPCDaughterB,
  // LFClusterStudiesTables::PDGCodeDaughterB,
  LFClusterStudiesTables::DcaToPVDaughterB,
  LFClusterStudiesTables::ItsClusterSizeDaughterB,
  LFClusterStudiesTables::TpcSignalDaughterB,
  LFClusterStudiesTables::TpcNclsDaughterB,
  LFClusterStudiesTables::TpcNsigmaElDaughterB,
  LFClusterStudiesTables::TpcNsigmaPiDaughterB,
  LFClusterStudiesTables::TpcNsigmaKaDaughterB,
  LFClusterStudiesTables::TpcNsigmaPrDaughterB,
  LFClusterStudiesTables::TpcNsigmaDeDaughterB,
  LFClusterStudiesTables::TpcNsigmaHeDaughterB,
  LFClusterStudiesTables::TofNsigmaElDaughterB,
  LFClusterStudiesTables::TofNsigmaPiDaughterB,
  LFClusterStudiesTables::TofNsigmaKaDaughterB,
  LFClusterStudiesTables::TofNsigmaPrDaughterB,
  LFClusterStudiesTables::TofNsigmaDeDaughterB,
  LFClusterStudiesTables::TofNsigmaHeDaughterB,
  LFClusterStudiesTables::Chi2itsDaughterB,
  LFClusterStudiesTables::Chi2tpcDaughterB,
  LFClusterStudiesTables::HasTPCDaughterB);

DECLARE_SOA_TABLE(
  ClStCascTable, "AOD", "CLSTCASCTABLE",
  LFClusterStudiesTables::PMother,
  // LFClusterStudiesTables::PtMother,
  LFClusterStudiesTables::EtaMother,
  LFClusterStudiesTables::PhiMother,
  LFClusterStudiesTables::MassMother,
  // LFClusterStudiesTables::PdgCodeMother,
  LFClusterStudiesTables::RadiusMother,
  LFClusterStudiesTables::DcaMotherPV,
  LFClusterStudiesTables::CosPAMother,
  LFClusterStudiesTables::PDaughterA,
  // LFClusterStudiesTables::PtDaughterA,
  LFClusterStudiesTables::EtaDaughterA,
  LFClusterStudiesTables::PhiDaughterA,
  LFClusterStudiesTables::PTPCDaughterA,
  LFClusterStudiesTables::PDGCodeDaughterA,
  LFClusterStudiesTables::DcaToPVDaughterA,
  LFClusterStudiesTables::ItsClusterSizeDaughterA,
  LFClusterStudiesTables::TpcSignalDaughterA,
  LFClusterStudiesTables::TpcNclsDaughterA,
  LFClusterStudiesTables::TpcNsigmaElDaughterA,
  LFClusterStudiesTables::TpcNsigmaPiDaughterA,
  LFClusterStudiesTables::TpcNsigmaKaDaughterA,
  LFClusterStudiesTables::TpcNsigmaPrDaughterA,
  LFClusterStudiesTables::TpcNsigmaDeDaughterA,
  LFClusterStudiesTables::TpcNsigmaHeDaughterA,
  LFClusterStudiesTables::TofNsigmaElDaughterA,
  LFClusterStudiesTables::TofNsigmaPiDaughterA,
  LFClusterStudiesTables::TofNsigmaKaDaughterA,
  LFClusterStudiesTables::TofNsigmaPrDaughterA,
  LFClusterStudiesTables::TofNsigmaDeDaughterA,
  LFClusterStudiesTables::TofNsigmaHeDaughterA,
  LFClusterStudiesTables::Chi2itsDaughterA,
  LFClusterStudiesTables::Chi2tpcDaughterA,
  LFClusterStudiesTables::HasTPCDaughterA);

DECLARE_SOA_TABLE(
  ClStNucTable, "AOD", "CLSTNUCTABLE",
  LFClusterStudiesTables::PDaughterA,
  // LFClusterStudiesTables::PtDaughterA,
  LFClusterStudiesTables::EtaDaughterA,
  LFClusterStudiesTables::PhiDaughterA,
  LFClusterStudiesTables::PTPCDaughterA,
  LFClusterStudiesTables::PDGCodeDaughterA,
  LFClusterStudiesTables::DcaToPVDaughterA,
  LFClusterStudiesTables::ItsClusterSizeDaughterA,
  LFClusterStudiesTables::TpcSignalDaughterA,
  LFClusterStudiesTables::TpcNclsDaughterA,
  LFClusterStudiesTables::TpcNsigmaElDaughterA,
  LFClusterStudiesTables::TpcNsigmaPiDaughterA,
  LFClusterStudiesTables::TpcNsigmaKaDaughterA,
  LFClusterStudiesTables::TpcNsigmaPrDaughterA,
  LFClusterStudiesTables::TpcNsigmaDeDaughterA,
  LFClusterStudiesTables::TpcNsigmaHeDaughterA,
  LFClusterStudiesTables::TofNsigmaElDaughterA,
  LFClusterStudiesTables::TofNsigmaPiDaughterA,
  LFClusterStudiesTables::TofNsigmaKaDaughterA,
  LFClusterStudiesTables::TofNsigmaPrDaughterA,
  LFClusterStudiesTables::TofNsigmaDeDaughterA,
  LFClusterStudiesTables::TofNsigmaHeDaughterA,
  LFClusterStudiesTables::Chi2itsDaughterA,
  LFClusterStudiesTables::Chi2tpcDaughterA,
  LFClusterStudiesTables::HasTPCDaughterA);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFCLUSTERSTUDIESTABLE_H_
