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

/// \file LFSpinAlignVecMesonTables.h
/// \brief Derived-data tables for K*0 and phi spin-alignment studies

#ifndef PWGLF_DATAMODEL_LFSPINALIGNVECMESONTABLES_H_
#define PWGLF_DATAMODEL_LFSPINALIGNVECMESONTABLES_H_

#include "PWGLF/DataModel/EPCalibrationTables.h"

#include "Common/DataModel/Centrality.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
DECLARE_SOA_TABLE(SpinAlignEvents, "AOD", "SAVMEVENT",
                  o2::soa::Index<>,
                  bc::RunNumber,
                  collision::PosZ,
                  cent::CentFT0M,
                  cent::CentFT0C,
                  epcalibrationtable::PsiFT0C,
                  epcalibrationtable::PsiFT0A,
                  epcalibrationtable::PsiTPC,
                  epcalibrationtable::QFT0C,
                  epcalibrationtable::QFT0A,
                  epcalibrationtable::QTPC);
using SpinAlignEvent = SpinAlignEvents::iterator;

namespace savmcand
{
DECLARE_SOA_INDEX_COLUMN(SpinAlignEvent, spinAlignEvent);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(PosPx, posPx, float);
DECLARE_SOA_COLUMN(PosPy, posPy, float);
DECLARE_SOA_COLUMN(PosPz, posPz, float);
DECLARE_SOA_COLUMN(PosTPCNSigma, posTPCNSigma, float);
DECLARE_SOA_COLUMN(PosTOFNSigma, posTOFNSigma, float);
DECLARE_SOA_COLUMN(PosHasTOF, posHasTOF, bool);
DECLARE_SOA_COLUMN(NegPx, negPx, float);
DECLARE_SOA_COLUMN(NegPy, negPy, float);
DECLARE_SOA_COLUMN(NegPz, negPz, float);
DECLARE_SOA_COLUMN(NegTPCNSigma, negTPCNSigma, float);
DECLARE_SOA_COLUMN(NegTOFNSigma, negTOFNSigma, float);
DECLARE_SOA_COLUMN(NegHasTOF, negHasTOF, bool);
DECLARE_SOA_COLUMN(Flavour, flavour, int8_t);
} // namespace savmcand

DECLARE_SOA_TABLE(SAKstarCands, "AOD", "SAKSTARCAND",
                  o2::soa::Index<>,
                  savmcand::SpinAlignEventId,
                  savmcand::Flavour,
                  savmcand::Mass,
                  savmcand::Px,
                  savmcand::Py,
                  savmcand::Pz,
                  savmcand::PosPx,
                  savmcand::PosPy,
                  savmcand::PosPz,
                  savmcand::NegPx,
                  savmcand::NegPy,
                  savmcand::NegPz,
                  savmcand::PosTPCNSigma,
                  savmcand::PosTOFNSigma,
                  savmcand::PosHasTOF,
                  savmcand::NegTPCNSigma,
                  savmcand::NegTOFNSigma,
                  savmcand::NegHasTOF);
using SAKstarCand = SAKstarCands::iterator;

DECLARE_SOA_TABLE(SAPhiCands, "AOD", "SAPHICAND",
                  o2::soa::Index<>,
                  savmcand::SpinAlignEventId,
                  savmcand::Mass,
                  savmcand::Px,
                  savmcand::Py,
                  savmcand::Pz,
                  savmcand::PosPx,
                  savmcand::PosPy,
                  savmcand::PosPz,
                  savmcand::NegPx,
                  savmcand::NegPy,
                  savmcand::NegPz,
                  savmcand::PosTPCNSigma,
                  savmcand::PosTOFNSigma,
                  savmcand::PosHasTOF,
                  savmcand::NegTPCNSigma,
                  savmcand::NegTOFNSigma,
                  savmcand::NegHasTOF);
using SAPhiCand = SAPhiCands::iterator;

namespace savmtrack
{
DECLARE_SOA_INDEX_COLUMN(SpinAlignEvent, spinAlignEvent);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float);
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNSigmaKa, float);
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNSigmaPi, float);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);
DECLARE_SOA_COLUMN(IsKaon, isKaon, bool);
DECLARE_SOA_COLUMN(IsPion, isPion, bool);
} // namespace savmtrack

DECLARE_SOA_TABLE(SAVMTracks, "AOD", "SAVMTRACK",
                  o2::soa::Index<>,
                  savmtrack::SpinAlignEventId,
                  savmtrack::Px,
                  savmtrack::Py,
                  savmtrack::Pz,
                  savmtrack::Sign,
                  savmtrack::TPCNSigmaKa,
                  savmtrack::TPCNSigmaPi,
                  savmtrack::TOFNSigmaKa,
                  savmtrack::TOFNSigmaPi,
                  savmtrack::HasTOF,
                  savmtrack::IsKaon,
                  savmtrack::IsPion);
using SAVMTrack = SAVMTracks::iterator;
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSPINALIGNVECMESONTABLES_H_
