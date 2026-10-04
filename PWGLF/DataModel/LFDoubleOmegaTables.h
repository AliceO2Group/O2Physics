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

#ifndef PWGLF_DATAMODEL_LFDOUBLEOMEGATABLES_H_
#define PWGLF_DATAMODEL_LFDOUBLEOMEGATABLES_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{

namespace DoubleOmegaTables
{
DECLARE_SOA_COLUMN(PtDoubleOmega, ptDoubleOmega, float);
DECLARE_SOA_COLUMN(EtaDoubleOmega, etaDoubleOmega, float);
DECLARE_SOA_COLUMN(PhiDoubleOmega, phiDoubleOmega, float);
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);
DECLARE_SOA_COLUMN(PrimaryVtxX, primaryVtxX, float);
DECLARE_SOA_COLUMN(PrimaryVtxY, primaryVtxY, float);
DECLARE_SOA_COLUMN(PrimaryVtxZ, primaryVtxZ, float);
DECLARE_SOA_COLUMN(CosPAOmega, cosPAOmega, float);
DECLARE_SOA_COLUMN(CosPADirectLambda, cosPADirectLambda, float);
DECLARE_SOA_COLUMN(CosPADoubleOmega, cosPADoubleOmega, float);
DECLARE_SOA_COLUMN(CosPAOmegaToPV, cosPAOmegaToPV, float);
DECLARE_SOA_COLUMN(CosPADirectLambdaToPV, cosPADirectLambdaToPV, float);
DECLARE_SOA_COLUMN(DCAxyOmegaToPV, dcaXYOmegaToPV, float);
DECLARE_SOA_COLUMN(DCAxyDirectLambdaToPV, dcaXYDirectLambdaToPV, float);
DECLARE_SOA_COLUMN(DCAxyDirectKaonToPV, dcaXYDirectKaonToPV, float);
// Transverse DCA to the primary vertex of the final-state tracks, ordered as
// the detector map: Omega kaon, Omega proton, Omega pion, direct proton, direct pion.
DECLARE_SOA_COLUMN(DCAxyOmegaKaonToPV, dcaXYOmegaKaonToPV, float);
DECLARE_SOA_COLUMN(DCAxyOmegaProtonToPV, dcaXYOmegaProtonToPV, float);
DECLARE_SOA_COLUMN(DCAxyOmegaPionToPV, dcaXYOmegaPionToPV, float);
DECLARE_SOA_COLUMN(DCAxyDirectProtonToPV, dcaXYDirectProtonToPV, float);
DECLARE_SOA_COLUMN(DCAxyDirectPionToPV, dcaXYDirectPionToPV, float);
DECLARE_SOA_COLUMN(TPCNSigmaDirectKaon, tpcNSigmaDirectKaon, float); // -999 when TPC is unavailable or the candidate is not reconstructed.
DECLARE_SOA_COLUMN(TOFNSigmaDirectKaon, tofNSigmaDirectKaon, float); // -999 when TOF is unavailable or the candidate is not reconstructed.
// Three bits per daughter (ITS, TPC, TOF), ordered as: Omega kaon, Omega
// proton, Omega pion, direct proton, direct pion, direct kaon (also for c.c.).
// Bits 18-31 are zero; generated-only rows have a zero bitmap.
DECLARE_SOA_COLUMN(DaughterDetectorMap, daughterDetectorMap, uint32_t);
DECLARE_SOA_COLUMN(MassDoubleOmega, massDoubleOmega, float);
DECLARE_SOA_COLUMN(MassOmega, massOmega, float);
DECLARE_SOA_COLUMN(MassXi, massXi, float);

DECLARE_SOA_COLUMN(GenPt, genPt, float);
DECLARE_SOA_COLUMN(GenEta, genEta, float);
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);
DECLARE_SOA_COLUMN(GenDecayLength, genDecayLength, float);
DECLARE_SOA_COLUMN(PdgDoubleOmega, pdgDoubleOmega, int);
DECLARE_SOA_COLUMN(GenMotherId, genMotherId, int64_t); // MC particle index within the input dataframe.
DECLARE_SOA_COLUMN(IsReco, isReco, bool);
} // namespace DoubleOmegaTables

#define DOUBLE_OMEGA_RECO_COLUMNS             \
  DoubleOmegaTables::PtDoubleOmega,           \
    DoubleOmegaTables::EtaDoubleOmega,        \
    DoubleOmegaTables::PhiDoubleOmega,        \
    DoubleOmegaTables::DecayVtxX,             \
    DoubleOmegaTables::DecayVtxY,             \
    DoubleOmegaTables::DecayVtxZ,             \
    DoubleOmegaTables::PrimaryVtxX,           \
    DoubleOmegaTables::PrimaryVtxY,           \
    DoubleOmegaTables::PrimaryVtxZ,           \
    DoubleOmegaTables::CosPAOmega,            \
    DoubleOmegaTables::CosPADirectLambda,     \
    DoubleOmegaTables::CosPADoubleOmega,      \
    DoubleOmegaTables::CosPAOmegaToPV,        \
    DoubleOmegaTables::CosPADirectLambdaToPV, \
    DoubleOmegaTables::DCAxyOmegaToPV,        \
    DoubleOmegaTables::DCAxyDirectLambdaToPV, \
    DoubleOmegaTables::DCAxyDirectKaonToPV,   \
    DoubleOmegaTables::DCAxyOmegaKaonToPV,    \
    DoubleOmegaTables::DCAxyOmegaProtonToPV,  \
    DoubleOmegaTables::DCAxyOmegaPionToPV,    \
    DoubleOmegaTables::DCAxyDirectProtonToPV, \
    DoubleOmegaTables::DCAxyDirectPionToPV,   \
    DoubleOmegaTables::TPCNSigmaDirectKaon,   \
    DoubleOmegaTables::TOFNSigmaDirectKaon,   \
    DoubleOmegaTables::DaughterDetectorMap,   \
    DoubleOmegaTables::MassDoubleOmega,       \
    DoubleOmegaTables::MassOmega,             \
    DoubleOmegaTables::MassXi

DECLARE_SOA_TABLE(DoubleOmegaTable, "AOD", "DOUBLEOMEGATABLE",
                  DOUBLE_OMEGA_RECO_COLUMNS);

DECLARE_SOA_TABLE(DoubleOmegaTableMC, "AOD", "DBLOMEGAMCTABLE",
                  DOUBLE_OMEGA_RECO_COLUMNS,
                  DoubleOmegaTables::GenPt,
                  DoubleOmegaTables::GenEta,
                  DoubleOmegaTables::GenPhi,
                  DoubleOmegaTables::GenDecayLength,
                  DoubleOmegaTables::PdgDoubleOmega,
                  DoubleOmegaTables::GenMotherId,
                  DoubleOmegaTables::IsReco);

#undef DOUBLE_OMEGA_RECO_COLUMNS

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFDOUBLEOMEGATABLES_H_
