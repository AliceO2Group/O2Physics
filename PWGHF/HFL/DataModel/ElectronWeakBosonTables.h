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

/// \file ElectronWeakBosonTables.h
/// \brief Derived-data tables for the W/Z -> e (+ hadron) analysis (taskElectronWeakBoson)
///
/// Three tables:
///  - HfWbColls     : one row per stored (triggered, skimmed) collision
///  - HfWbTracks    : charged tracks of that collision (associated hadrons + every electron-candidate track)
///  - HfWbElectrons : high-pT tracks only, with what the electron / Z reconstruction needs on top of HfWbTracks:
///                    track parameters + covariance (KFParticle), track isolation sums,
///                    closest EMCal cluster and the cluster-isolation sum.
///
/// The getter names deliberately follow the AO2D ones (centFT0M(), multNTracksPV(), isGlobalTrackWoPtEta(),
/// tpcNClsCrossedRows(), x(), alpha(), cYY(), ...) so that the helpers already used on AO2D
/// (o2::hf_centrality::getCentralityColl, getTrackParCov, isSelectedTrackKine of the task) work unchanged.
///
/// \author S. Sakai & S. Ito (Univ. of Tsukuba)

#ifndef PWGHF_HFL_DATAMODEL_ELECTRONWEAKBOSONTABLES_H_
#define PWGHF_HFL_DATAMODEL_ELECTRONWEAKBOSONTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cmath>
#include <cstdint>

namespace o2::aod
{
// ---------------------------------------------------------------------------
// Collisions
// ---------------------------------------------------------------------------
namespace hf_wb_coll
{
DECLARE_SOA_COLUMN(PosZ, posZ, float);                 //! primary-vertex z (cm)
DECLARE_SOA_COLUMN(Sel8, sel8, bool);                  //! sel8 decision (stored, not applied at production)
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);         //! FT0A centrality percentile
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);         //! FT0C centrality percentile
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);         //! FT0M centrality percentile
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);         //! FV0A centrality percentile
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);         //! FT0M raw multiplicity
DECLARE_SOA_COLUMN(MultNTracksPV, multNTracksPV, int); //! number of PV contributors
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);         //! run number
DECLARE_SOA_COLUMN(MagField, magField, float);         //! nominal Bz (kG), for KFParticle::SetField
DECLARE_SOA_COLUMN(RhoUE, rhoUE, float);               //! UE density from random cones (GeV / unit area)
} // namespace hf_wb_coll

DECLARE_SOA_TABLE(HfWbColls, "AOD", "HFWBCOLL", //! collisions for the W/Z -> e analysis
                  o2::soa::Index<>,
                  hf_wb_coll::PosZ,
                  hf_wb_coll::Sel8,
                  hf_wb_coll::CentFT0A,
                  hf_wb_coll::CentFT0C,
                  hf_wb_coll::CentFT0M,
                  hf_wb_coll::CentFV0A,
                  hf_wb_coll::MultFT0M,
                  hf_wb_coll::MultNTracksPV,
                  hf_wb_coll::RunNumber,
                  hf_wb_coll::MagField,
                  hf_wb_coll::RhoUE);
using HfWbColl = HfWbColls::iterator;

// ---------------------------------------------------------------------------
// Tracks
// ---------------------------------------------------------------------------
namespace hf_wb_track
{
DECLARE_SOA_INDEX_COLUMN(HfWbColl, hfWbColl);                                                       //! derived collision
DECLARE_SOA_COLUMN(Pt, pt, float);                                                                  //! pT (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                                                                //! pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                                                                //! azimuth (0..2pi)
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                                             //! charge sign
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);                                                            //! DCA xy (cm)
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                                                              //! DCA z (cm)
DECLARE_SOA_COLUMN(TpcNClsFound, tpcNClsFound, int16_t);                                            //! TPC found clusters
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, int16_t);                                //! TPC crossed rows
DECLARE_SOA_COLUMN(ItsNCls, itsNCls, uint8_t);                                                      //! ITS clusters
DECLARE_SOA_COLUMN(TpcChi2NCl, tpcChi2NCl, float);                                                  //! TPC chi2 / cluster
DECLARE_SOA_COLUMN(ItsChi2NCl, itsChi2NCl, float);                                                  //! ITS chi2 / cluster
DECLARE_SOA_COLUMN(TpcNSigmaEl, tpcNSigmaEl, float);                                                //! TPC n sigma electron
DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, float);                                                //! TPC n sigma pion
DECLARE_SOA_COLUMN(IsGlobalTrackWoPtEta, isGlobalTrackWoPtEta, bool);                               //! global track without pT/eta cut
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * std::cosh(eta); }); //! momentum
} // namespace hf_wb_track

DECLARE_SOA_TABLE(HfWbTracks, "AOD", "HFWBTRACK", //! tracks for the W/Z -> e analysis
                  o2::soa::Index<>,
                  hf_wb_track::HfWbCollId,
                  hf_wb_track::Pt,
                  hf_wb_track::Eta,
                  hf_wb_track::Phi,
                  hf_wb_track::Sign,
                  hf_wb_track::DcaXY,
                  hf_wb_track::DcaZ,
                  hf_wb_track::TpcNClsFound,
                  hf_wb_track::TpcNClsCrossedRows,
                  hf_wb_track::ItsNCls,
                  hf_wb_track::TpcChi2NCl,
                  hf_wb_track::ItsChi2NCl,
                  hf_wb_track::TpcNSigmaEl,
                  hf_wb_track::TpcNSigmaPi,
                  hf_wb_track::IsGlobalTrackWoPtEta,
                  hf_wb_track::P<hf_wb_track::Pt, hf_wb_track::Eta>);
using HfWbTrack = HfWbTracks::iterator;

// ---------------------------------------------------------------------------
// Electron candidates (high-pT tracks)
// ---------------------------------------------------------------------------
namespace hf_wb_ele
{
DECLARE_SOA_INDEX_COLUMN(HfWbTrack, hfWbTrack); //! row in HfWbTracks (kinematics, PID, quality)

// track parametrisation at the reference point (same getters as aod::TracksIU/TracksCov -> getTrackParCov works)
DECLARE_SOA_COLUMN(X, x, float);
DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Z, z, float);
DECLARE_SOA_COLUMN(Snp, snp, float);
DECLARE_SOA_COLUMN(Tgl, tgl, float);
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);
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

// track isolation, summed over ALL tracks of the collision within rIsolation (production value), self included
DECLARE_SOA_COLUMN(PSumCone, pSumCone, float);   //! sum of p in the cone (GeV/c)
DECLARE_SOA_COLUMN(NTrkCone, nTrkCone, int16_t); //! number of tracks in the cone

// closest EMCal cluster among the matched clusters passing the production time window
// (no rMatchMax cut applied: the analysis cuts on dEta/dPhi)
DECLARE_SOA_COLUMN(HasEmcMatch, hasEmcMatch, bool);      //! a cluster candidate exists
DECLARE_SOA_COLUMN(EmcDEta, emcDEta, float);             //! eta(track at EMCal) - eta(cluster)
DECLARE_SOA_COLUMN(EmcDPhi, emcDPhi, float);             //! phi(track at EMCal) - phi(cluster), in [-pi, pi)
DECLARE_SOA_COLUMN(EmcEnergy, emcEnergy, float);         //! cluster energy (GeV)
DECLARE_SOA_COLUMN(EmcEta, emcEta, float);               //! cluster eta
DECLARE_SOA_COLUMN(EmcPhi, emcPhi, float);               //! cluster phi
DECLARE_SOA_COLUMN(EmcM02, emcM02, float);               //! cluster M02
DECLARE_SOA_COLUMN(EmcTime, emcTime, float);             //! cluster time (ns)
DECLARE_SOA_COLUMN(EmcNCells, emcNCells, int16_t);       //! cluster number of cells
DECLARE_SOA_COLUMN(EmcNMatchInR, emcNMatchInR, uint8_t); //! clusters within the production rMatchMax (QA)
DECLARE_SOA_COLUMN(EmcConeEnergy, emcConeEnergy, float); //! cluster energy in the cone, the cluster itself EXCLUDED, no UE subtraction
} // namespace hf_wb_ele

DECLARE_SOA_TABLE(HfWbElectrons, "AOD", "HFWBELE", //! electron candidates for the W/Z -> e analysis
                  o2::soa::Index<>,
                  hf_wb_track::HfWbCollId,
                  hf_wb_ele::HfWbTrackId,
                  hf_wb_ele::X, hf_wb_ele::Alpha,
                  hf_wb_ele::Y, hf_wb_ele::Z, hf_wb_ele::Snp, hf_wb_ele::Tgl, hf_wb_ele::Signed1Pt,
                  hf_wb_ele::CYY, hf_wb_ele::CZY, hf_wb_ele::CZZ,
                  hf_wb_ele::CSnpY, hf_wb_ele::CSnpZ, hf_wb_ele::CSnpSnp,
                  hf_wb_ele::CTglY, hf_wb_ele::CTglZ, hf_wb_ele::CTglSnp, hf_wb_ele::CTglTgl,
                  hf_wb_ele::C1PtY, hf_wb_ele::C1PtZ, hf_wb_ele::C1PtSnp, hf_wb_ele::C1PtTgl, hf_wb_ele::C1Pt21Pt2,
                  hf_wb_ele::PSumCone,
                  hf_wb_ele::NTrkCone,
                  hf_wb_ele::HasEmcMatch,
                  hf_wb_ele::EmcDEta,
                  hf_wb_ele::EmcDPhi,
                  hf_wb_ele::EmcEnergy,
                  hf_wb_ele::EmcEta,
                  hf_wb_ele::EmcPhi,
                  hf_wb_ele::EmcM02,
                  hf_wb_ele::EmcTime,
                  hf_wb_ele::EmcNCells,
                  hf_wb_ele::EmcNMatchInR,
                  hf_wb_ele::EmcConeEnergy);
using HfWbElectron = HfWbElectrons::iterator;
} // namespace o2::aod

#endif // PWGHF_HFL_DATAMODEL_ELECTRONWEAKBOSONTABLES_H_
