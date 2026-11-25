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

/// \file LFDoubleOmegaTrackTables.h
/// \brief Tables of doubleOmegaTrackTreeCreator: the training sample of the double-Omega track selection
///        (DOmegaTrackTrain) and the V0 and cascade candidates built from the tracks selected by the track
///        BDT (DOmegaV0s, DOmegaCascades), analogous to aod::V0s and aod::Cascades (same index getters:
///        collision, posTrack, negTrack, v0, bachelor).

#ifndef PWGLF_DATAMODEL_LFDOUBLEOMEGATRACKTABLES_H_
#define PWGLF_DATAMODEL_LFDOUBLEOMEGATRACKTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace domegatrack
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(TpcInnerParam, tpcInnerParam, float); //! Momentum at the TPC inner wall (PID)
DECLARE_SOA_COLUMN(TpcSignal, tpcSignal, float);
DECLARE_SOA_COLUMN(TpcNSigmaPr, tpcNSigmaPr, float);
DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TpcNSigmaKa, tpcNSigmaKa, float);
DECLARE_SOA_COLUMN(TofNSigmaPr, tofNSigmaPr, float); //! -999 without TOF
DECLARE_SOA_COLUMN(TofNSigmaPi, tofNSigmaPi, float); //! -999 without TOF
DECLARE_SOA_COLUMN(TofNSigmaKa, tofNSigmaKa, float); //! -999 without TOF
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);             //! Signed transverse DCA to the primary vertex (cm), -999 if the propagation failed
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);               //! Longitudinal DCA to the primary vertex (cm)
DECLARE_SOA_COLUMN(SigmaDcaXY, sigmaDcaXY, float);   //! Uncertainty of DcaXY, track and vertex (cm)
DECLARE_SOA_COLUMN(SigmaDcaZ, sigmaDcaZ, float);     //! Uncertainty of DcaZ, track and vertex (cm)
DECLARE_SOA_COLUMN(ItsClusterMap, itsClusterMap, uint8_t);
DECLARE_SOA_COLUMN(ItsClusterSizes, itsClusterSizes, uint32_t);
DECLARE_SOA_COLUMN(ItsChi2NCl, itsChi2NCl, float);
DECLARE_SOA_COLUMN(TpcNClsFound, tpcNClsFound, int16_t);
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, int16_t);
DECLARE_SOA_COLUMN(TpcNClsShared, tpcNClsShared, uint8_t);
DECLARE_SOA_COLUMN(TpcChi2NCl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t); //! ITS, TPC, TRD, TOF bits
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);
DECLARE_SOA_COLUMN(NumContrib, numContrib, uint16_t); //! Primary-vertex contributors of the collision (multiplicity proxy)
DECLARE_SOA_COLUMN(PosZ, posZ, float);                //! Primary-vertex z of the collision (cm)
DECLARE_SOA_COLUMN(Role, role, int8_t);               //! See doubleOmegaTrackTreeCreator.cxx
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);            //! Of the MC particle, 0 without MC particle
DECLARE_SOA_COLUMN(MotherPdgCode, motherPdgCode, int);
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(PtMC, ptMC, float);
DECLARE_SOA_COLUMN(ProductionRadius, productionRadius, float); //! Transverse radius of the MC production vertex (cm)
DECLARE_SOA_COLUMN(McKey, mcKey, int64_t);                     //! MC particle, unique across dataframes (dataframe << 32 | index), -1 without MC particle
DECLARE_SOA_COLUMN(AncestorPdgCode, ancestorPdgCode, int);     //! First particle of the decay chain (see doubleOmegaTrackTreeCreator.cxx), 0 without MC particle
DECLARE_SOA_COLUMN(AncestorPtMC, ancestorPtMC, float);         //! Its generated pT (GeV/c), -999 without MC particle
} // namespace domegatrack

DECLARE_SOA_TABLE(DOmegaTrackTrain, "AOD", "DOMEGATRKTRAIN", //! Tracks for the training of the double-Omega track selection
                  domegatrack::Pt, domegatrack::Eta, domegatrack::Phi, domegatrack::Sign,
                  domegatrack::TpcInnerParam, domegatrack::TpcSignal,
                  domegatrack::TpcNSigmaPr, domegatrack::TpcNSigmaPi, domegatrack::TpcNSigmaKa,
                  domegatrack::TofNSigmaPr, domegatrack::TofNSigmaPi, domegatrack::TofNSigmaKa,
                  domegatrack::DcaXY, domegatrack::DcaZ, domegatrack::SigmaDcaXY, domegatrack::SigmaDcaZ,
                  domegatrack::ItsClusterMap, domegatrack::ItsClusterSizes, domegatrack::ItsChi2NCl,
                  domegatrack::TpcNClsFound, domegatrack::TpcNClsCrossedRows,
                  domegatrack::TpcNClsShared, domegatrack::TpcChi2NCl,
                  domegatrack::DetectorMap, domegatrack::IsPVContributor,
                  domegatrack::NumContrib, domegatrack::PosZ,
                  domegatrack::Role, domegatrack::PdgCode, domegatrack::MotherPdgCode,
                  domegatrack::IsPhysicalPrimary, domegatrack::PtMC, domegatrack::ProductionRadius, domegatrack::McKey,
                  domegatrack::AncestorPdgCode, domegatrack::AncestorPtMC);

// V0 and cascade candidates from the tracks selected by the track BDT. The track indices point to the Tracks
// of the AO2D the candidates were built from: like aod::V0s and aod::Cascades, the tables are used together
// with it.
namespace domegav0
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //! Collision index
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //! Positive track (proton of a Lambda, pion of an anti-Lambda)
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //! Negative track (pion of a Lambda, antiproton of an anti-Lambda)
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                 //! +1 Lambda, -1 anti-Lambda
} // namespace domegav0

DECLARE_SOA_TABLE(DOmegaV0s, "AOD", "DOMEGAV0", //! Lambda candidates, analogous to aod::V0s
                  o2::soa::Index<>, domegav0::CollisionId, domegav0::PosTrackId, domegav0::NegTrackId, domegav0::Sign);

namespace domegacascade
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                     //! Collision index
DECLARE_SOA_INDEX_COLUMN_FULL(V0, v0, int, DOmegaV0s, "");          //! Lambda candidate
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int, Tracks, ""); //! Bachelor kaon
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                             //! Charge of the Omega: -1 Omega-, +1 Omega+
} // namespace domegacascade

DECLARE_SOA_TABLE(DOmegaCascades, "AOD", "DOMEGACASCADE", //! Omega candidates, analogous to aod::Cascades
                  o2::soa::Index<>, domegacascade::CollisionId, domegacascade::V0Id, domegacascade::BachelorId, domegacascade::Sign);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFDOUBLEOMEGATRACKTABLES_H_
