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

/// \file pidFeatureExtractor.h
/// \brief Data model for the PID feature extractor: PidFeaturesData /
///        PidFeaturesMc table definitions. Included with a plain quoted
///        filename ("pidFeatureExtractor.h"), which the compiler resolves
///        relative to the including .cxx's own directory first - so this
///        header and every .cxx that uses it must stay in the same folder;
///        no repo-root-relative path to get wrong.
///
/// \author Robert Forynski

#ifndef PID_FEATURE_EXTRACTOR_H_
#define PID_FEATURE_EXTRACTOR_H_

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace pidfeat
{
// Kinematics not already stored as static columns upstream (Pt, TrackType,
// DcaXY, DcaZ, TPCSignal, TRDSignal, TRDPattern, TrackEtaEMCAL,
// TrackPhiEMCAL, ITSClusterSizes, ITSChi2NCl, TPCNClsFound are reused
// directly from aod::track:: in the table definitions below).
DECLARE_SOA_COLUMN(P, p, float);       //! Track momentum magnitude (GeV/c)
DECLARE_SOA_COLUMN(Px, px, float);     //! Track x-momentum (GeV/c)
DECLARE_SOA_COLUMN(Py, py, float);     //! Track y-momentum (GeV/c)
DECLARE_SOA_COLUMN(Pz, pz, float);     //! Track z-momentum (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);   //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);   //! Azimuthal angle
DECLARE_SOA_COLUMN(Sign, sign, float); //! Track charge sign

// Event level, duplicated per track for a self-contained flat table.
DECLARE_SOA_COLUMN(Vz, vz, float); //! Collision vertex z (cm)

// TOF mass: not exposed as a static column upstream.
DECLARE_SOA_COLUMN(TofMass, tofMass, float); //! TOF-reconstructed mass, NaN if !hasTOF()

// Detector-presence flags. Kept explicit (rather than inferred solely from
// NaN sentinels) because the downstream FSE/attention model conditions on a
// per-track detector mask.
DECLARE_SOA_COLUMN(HasTPC, hasTPC, uint8_t);     //!
DECLARE_SOA_COLUMN(HasTOF, hasTOF, uint8_t);     //!
DECLARE_SOA_COLUMN(HasTRD, hasTRD, uint8_t);     //! trdPattern() > 0
DECLARE_SOA_COLUMN(HasEMCal, hasEMCal, uint8_t); //! trackEtaEmcal() in acceptance
DECLARE_SOA_COLUMN(HasHMPID, hasHMPID, uint8_t); //! matched in the sparse HMPID table

// HMPID: sparse table, matched by hand per collision (see buildHmpidMap() in
// pidFeatureExtractor.cxx).
DECLARE_SOA_COLUMN(HmpidSignal, hmpidSignal, float);   //! Cherenkov angle (rad), NaN if !hasHMPID
DECLARE_SOA_COLUMN(HmpidQMip, hmpidQMip, float);       //!
DECLARE_SOA_COLUMN(HmpidNPhotons, hmpidNPhotons, int); //!
DECLARE_SOA_COLUMN(HmpidClusSize, hmpidClusSize, int); //!
DECLARE_SOA_COLUMN(HmpidMom, hmpidMom, float);         //!

// Bayesian PID posteriors. NaN when computeBayesianPid is false, or when the
// track doesn't have TPC (see computeBayesianProbs() in
// pidFeatureExtractor.cxx) - never a value that could be mistaken for a real
// posterior.
DECLARE_SOA_COLUMN(BayesProbPi, bayesProbPi, float); //!
DECLARE_SOA_COLUMN(BayesProbKa, bayesProbKa, float); //!
DECLARE_SOA_COLUMN(BayesProbPr, bayesProbPr, float); //!
DECLARE_SOA_COLUMN(BayesProbEl, bayesProbEl, float); //!

// MC truth.
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, uint8_t); //!
} // namespace pidfeat

// Real/raw data: reconstructed features only.
DECLARE_SOA_TABLE(PidFeaturesData, "AOD", "PIDFEATDATA", //!
                  o2::soa::Index<>,
                  pidfeat::P, aod::track::Pt, pidfeat::Px, pidfeat::Py, pidfeat::Pz,
                  pidfeat::Eta, pidfeat::Phi, pidfeat::Sign, aod::track::TrackType,
                  pidfeat::Vz, aod::cent::CentFT0C,
                  aod::track::DcaXY, aod::track::DcaZ,
                  pidfeat::HasTPC, aod::track::TPCSignal,
                  pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr, pidtpc::TPCNSigmaEl,
                  aod::track::TPCNClsFound, aod::track::TPCChi2NCl,
                  pidfeat::HasTOF, pidfeat::TofMass, aod::pidtofbeta::Beta,
                  pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr, pidtof::TOFNSigmaEl,
                  pidfeat::HasTRD, aod::track::TRDSignal, aod::track::TRDChi2, aod::track::TRDPattern,
                  aod::track::ITSClusterSizes, aod::track::ITSChi2NCl,
                  pidfeat::HasEMCal, aod::track::TrackEtaEMCAL, aod::track::TrackPhiEMCAL,
                  pidfeat::HasHMPID, pidfeat::HmpidSignal, pidfeat::HmpidQMip,
                  pidfeat::HmpidNPhotons, pidfeat::HmpidClusSize, pidfeat::HmpidMom,
                  pidfeat::BayesProbPi, pidfeat::BayesProbKa, pidfeat::BayesProbPr, pidfeat::BayesProbEl);

// MC: reconstructed features + truth. Same reconstructed columns as
// PidFeaturesData plus PdgCode/IsPhysicalPrimary - kept as a distinct table
// (rather than optional columns on one table) because O2 tables have a
// fixed schema.
DECLARE_SOA_TABLE(PidFeaturesMc, "AOD", "PIDFEATMC", //!
                  o2::soa::Index<>,
                  pidfeat::P, aod::track::Pt, pidfeat::Px, pidfeat::Py, pidfeat::Pz,
                  pidfeat::Eta, pidfeat::Phi, pidfeat::Sign, aod::track::TrackType,
                  pidfeat::Vz, aod::cent::CentFT0C,
                  aod::track::DcaXY, aod::track::DcaZ,
                  pidfeat::HasTPC, aod::track::TPCSignal,
                  pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr, pidtpc::TPCNSigmaEl,
                  aod::track::TPCNClsFound, aod::track::TPCChi2NCl,
                  pidfeat::HasTOF, pidfeat::TofMass, aod::pidtofbeta::Beta,
                  pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr, pidtof::TOFNSigmaEl,
                  pidfeat::HasTRD, aod::track::TRDSignal, aod::track::TRDChi2, aod::track::TRDPattern,
                  aod::track::ITSClusterSizes, aod::track::ITSChi2NCl,
                  pidfeat::HasEMCal, aod::track::TrackEtaEMCAL, aod::track::TrackPhiEMCAL,
                  pidfeat::HasHMPID, pidfeat::HmpidSignal, pidfeat::HmpidQMip,
                  pidfeat::HmpidNPhotons, pidfeat::HmpidClusSize, pidfeat::HmpidMom,
                  pidfeat::BayesProbPi, pidfeat::BayesProbKa, pidfeat::BayesProbPr, pidfeat::BayesProbEl,
                  aod::mcparticle::PdgCode, pidfeat::IsPhysicalPrimary);
} // namespace o2::aod

#endif // PID_FEATURE_EXTRACTOR_H_
