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
DECLARE_SOA_COLUMN(PtDoubleOmega, ptDoubleOmega, float); // Signed transverse momentum of the double-Omega candidate, signed by the charge of the direct kaon.
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

DECLARE_SOA_COLUMN(GenPt, genPt, float); // Signed transverse momentum of the generated double-Omega candidate, signed by the charge of the direct kaon.
DECLARE_SOA_COLUMN(GenEta, genEta, float);
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);
DECLARE_SOA_COLUMN(GenDecayLength, genDecayLength, float);
DECLARE_SOA_COLUMN(PdgDoubleOmega, pdgDoubleOmega, int);
DECLARE_SOA_COLUMN(GenMotherId, genMotherId, int64_t); // MC particle index within the input dataframe.
DECLARE_SOA_COLUMN(IsReco, isReco, bool);

DECLARE_SOA_COLUMN(MixingType, mixingType, uint8_t); // 1: Omega, 2: Lambda, 3: kaon from another collision; 4: Omega, Lambda and kaon from three different collisions.
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

DECLARE_SOA_TABLE(DoubleOmegaMixedTable, "AOD", "DBLOMEGAMIXTABLE",
                  DOUBLE_OMEGA_RECO_COLUMNS,
                  DoubleOmegaTables::MixingType);

#undef DOUBLE_OMEGA_RECO_COLUMNS

// Derived data: Omega, Lambda and kaon candidates of the collisions with at least one Omega, one
// Lambda and one kaon candidate not sharing tracks with it. The tables only index each other, so
// they are self-contained and dataframes can be merged.
namespace domega
{
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float); //! FT0M multiplicity of the collision
DECLARE_SOA_COLUMN(Bz, bz, float);             //! Nominal magnetic field (kG)
} // namespace domega

DECLARE_SOA_TABLE(DOmegaColls, "AOD", "DOMEGACOLL", //! Collisions with double-Omega building blocks
                  o2::soa::Index<>,
                  collision::PosX, collision::PosY, collision::PosZ,
                  domega::MultFT0M,
                  bc::RunNumber,
                  domega::Bz);

namespace domega
{
DECLARE_SOA_INDEX_COLUMN_FULL(DOmegaColl, dOmegaColl, int, DOmegaColls, ""); //! Collision index

// Track parametrisation covariance, in the order expected by getTrackParCov.
DECLARE_SOA_COLUMN(CYY, cYY, float);
DECLARE_SOA_COLUMN(CZY, cZY, float);
DECLARE_SOA_COLUMN(CZZ, cZZ, float);
DECLARE_SOA_COLUMN(CSnpY, cSnpY, float);
DECLARE_SOA_COLUMN(CSnpZ, cSnpZ, float);
DECLARE_SOA_COLUMN(CSnpSnp, cSnpSnp, float);
DECLARE_SOA_COLUMN(CTglY, cTglY, float);
DECLARE_SOA_COLUMN(CTglZ, cTglZ, float);
DECLARE_SOA_COLUMN(CTglSnp, cTglSnp, float);
DECLARE_SOA_COLUMN(CTglTgl, cTglTgl, float);
DECLARE_SOA_COLUMN(C1PtY, c1PtY, float);
DECLARE_SOA_COLUMN(C1PtZ, c1PtZ, float);
DECLARE_SOA_COLUMN(C1PtSnp, c1PtSnp, float);
DECLARE_SOA_COLUMN(C1PtTgl, c1PtTgl, float);
DECLARE_SOA_COLUMN(C1Pt21Pt2, c1Pt21Pt2, float);

DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float);   //! -999 when TPC is unavailable
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNSigmaKa, float);   //! -999 when TOF is unavailable
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t); //! ITS, TPC, TOF bits
DECLARE_SOA_COLUMN(DcaXYToPV, dcaXYToPV, float);       //! Transverse DCA to the primary vertex of its collision
DECLARE_SOA_COLUMN(IsKaon, isKaon, bool);              //! Selected as a direct-kaon candidate
} // namespace domega

#define DOMEGA_TRACK_PAR_COV_COLUMNS                                                      \
  track::X, track::Alpha, track::Y, track::Z, track::Snp, track::Tgl, track::Signed1Pt,   \
    domega::CYY, domega::CZY, domega::CZZ, domega::CSnpY, domega::CSnpZ, domega::CSnpSnp, \
    domega::CTglY, domega::CTglZ, domega::CTglSnp, domega::CTglTgl,                       \
    domega::C1PtY, domega::C1PtZ, domega::C1PtSnp, domega::C1PtTgl, domega::C1Pt21Pt2

DECLARE_SOA_TABLE(DOmegaTracks, "AOD", "DOMEGATRACK", //! Direct-kaon candidates and Omega/Lambda daughters, at the innermost update
                  o2::soa::Index<>,
                  domega::DOmegaCollId,
                  DOMEGA_TRACK_PAR_COV_COLUMNS,
                  domega::TPCNSigmaKa,
                  domega::TOFNSigmaKa,
                  domega::DetectorMap,
                  domega::DcaXYToPV,
                  domega::IsKaon,
                  track::Sign<track::Signed1Pt>);

namespace domega
{
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, DOmegaTracks, "_Pos");  //! Positive daughter
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, DOmegaTracks, "_Neg");  //! Negative daughter
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int, DOmegaTracks, "_Bach"); //! Bachelor kaon of the Omega

DECLARE_SOA_COLUMN(DecayX, decayX, float);
DECLARE_SOA_COLUMN(DecayY, decayY, float);
DECLARE_SOA_COLUMN(DecayZ, decayZ, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(MassOmega, massOmega, float);
DECLARE_SOA_COLUMN(MassXi, massXi, float);
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);
DECLARE_SOA_COLUMN(DcaDaughters, dcaDaughters, float); //! DCA between the daughters
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                //! Charge of the Omega (Omega table) or of the Omega the Lambda can be paired with (Lambda table)
} // namespace domega

DECLARE_SOA_TABLE(DOmegaOmegas, "AOD", "DOMEGAOMEGA", //! Omega candidates; track parameters at the Omega decay vertex
                  o2::soa::Index<>,
                  domega::DOmegaCollId,
                  domega::PosTrackId, domega::NegTrackId, domega::BachelorId,
                  DOMEGA_TRACK_PAR_COV_COLUMNS,
                  domega::DecayX, domega::DecayY, domega::DecayZ,
                  domega::Px, domega::Py, domega::Pz,
                  domega::MassOmega, domega::MassXi,
                  domega::DcaDaughters,
                  domega::Sign);

DECLARE_SOA_TABLE(DOmegaLambdas, "AOD", "DOMEGALAMBDA", //! Lambda candidates; neutral track parameters at the Lambda decay vertex
                  o2::soa::Index<>,
                  domega::DOmegaCollId,
                  domega::PosTrackId, domega::NegTrackId,
                  DOMEGA_TRACK_PAR_COV_COLUMNS,
                  domega::DecayX, domega::DecayY, domega::DecayZ,
                  domega::Px, domega::Py, domega::Pz,
                  domega::MassLambda,
                  domega::DcaDaughters,
                  domega::Sign);

#undef DOMEGA_TRACK_PAR_COV_COLUMNS

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFDOUBLEOMEGATABLES_H_
