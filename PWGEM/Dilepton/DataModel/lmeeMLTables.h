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

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"

#include <vector>

#ifndef PWGEM_DILEPTON_DATAMODEL_LMEEMLTABLES_H_
#define PWGEM_DILEPTON_DATAMODEL_LMEEMLTABLES_H_

namespace o2::aod
{

namespace pwgem::dilepton::ml
{
enum class PID_Label : uint8_t {
  kElectron = 0,
  kMuon = 1,
  kPion = 2,
  kKaon = 3,
  kProton = 4,
}; // this can be used for eID.

enum class Track_Type : uint8_t {
  kPrimary = 0,
  kSecondary = 1,
}; // this can be used for selecting electron from primary or photon conversion.

} // namespace pwgem::dilepton::ml

namespace emmltrack
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                   //!
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, float);               //!
DECLARE_SOA_COLUMN(PIDLabel, pidlabel, uint8_t);                     //!
DECLARE_SOA_COLUMN(TrackType, tracktype, uint8_t);                   //!
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, uint8_t);             //!
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t); //!
DECLARE_SOA_COLUMN(TPCNClsPID, tpcNClsPID, uint8_t);                 //!
DECLARE_SOA_COLUMN(IsForValidation, isForValidation, bool);          //!
DECLARE_SOA_COLUMN(Sign, sign, short);                               //!
DECLARE_SOA_COLUMN(P, p, float);                                     //!
// DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITS, meanClusterSizeITS, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSob, meanClusterSizeITSob, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 3; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
} // namespace emmltrack

// reconstructed track information
DECLARE_SOA_TABLE(EMTracksForMLPID, "AOD", "EMTRACKMLPID", //!
                  o2::soa::Index<>, collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, emmltrack::HadronicRate,
                  emmltrack::P, track::Tgl, emmltrack::Sign,
                  track::TPCNClsFindable, emmltrack::TPCNClsFound, emmltrack::TPCNClsCrossedRows, emmltrack::TPCNClsPID,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal,
                  pidtofbeta::Beta,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::TOFChi2, track::DetectorMap, emmltrack::PIDLabel,

                  // dynamic column
                  emmltrack::MeanClusterSizeITS<track::ITSClusterSizes>,
                  emmltrack::MeanClusterSizeITSob<track::ITSClusterSizes>);

DECLARE_SOA_TABLE(EMPIDsEl, "AOD", "EMPIDEL", pidtpc::TPCNSigmaEl, pidtof::TOFNSigmaEl); // Joinable with EMTracksForMLPID
DECLARE_SOA_TABLE(EMPIDsPi, "AOD", "EMPIDPI", pidtpc::TPCNSigmaPi, pidtof::TOFNSigmaPi); // Joinable with EMTracksForMLPID
DECLARE_SOA_TABLE(EMPIDsKa, "AOD", "EMPIDKA", pidtpc::TPCNSigmaKa, pidtof::TOFNSigmaKa); // Joinable with EMTracksForMLPID
DECLARE_SOA_TABLE(EMPIDsPr, "AOD", "EMPIDPR", pidtpc::TPCNSigmaPr, pidtof::TOFNSigmaPr); // Joinable with EMTracksForMLPID

// iterators
using EMTrackForMLPID = EMTracksForMLPID::iterator;
using EMPIDEl = EMPIDsEl::iterator;
using EMPIDPi = EMPIDsPi::iterator;
using EMPIDKa = EMPIDsKa::iterator;
using EMPIDPr = EMPIDsPr::iterator;

namespace emmlfwdtrack
{
DECLARE_SOA_COLUMN(PtMatchedMCHMID, ptMatchedMCHMID, float);          //! pt of MCH-MID track in MFT-MCH-MID track at PV
DECLARE_SOA_COLUMN(EtaMatchedMCHMID, etaMatchedMCHMID, float);        //! eta of MCH-MID track in MFT-MCH-MID track at PV
DECLARE_SOA_COLUMN(PhiMatchedMCHMID, phiMatchedMCHMID, float);        //! phi of MCH-MID track in MFT-MCH-MID track at PV
DECLARE_SOA_COLUMN(XMatchedMCHMIDatMP, xMatchedMCHMIDatMP, float);    //! x of MCH-MID track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(YMatchedMCHMIDatMP, yMatchedMCHMIDatMP, float);    //! y of MCH-MID track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(XMatchedMFTatMP, xMatchedMFTatMP, float);          //! x of MFT track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(YMatchedMFTatMP, yMatchedMFTatMP, float);          //! y of MFT track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                               //!
DECLARE_SOA_COLUMN(Chi2MFT, chi2MFT, float);                          //! chi2 of MFT standalone track
DECLARE_SOA_COLUMN(NClustersMFT, nClustersMFT, uint8_t);              //!
DECLARE_SOA_COLUMN(IsPrimary, isPrimary, bool);                       //!
DECLARE_SOA_COLUMN(IsCorrectMatchMFTMCH, isCorrectMatchMFTMCH, bool); //!
} // namespace emmlfwdtrack

DECLARE_SOA_TABLE_VERSIONED(EMFwdTracksForML_000, "AOD", "EMFWDTRKML", 0, //!
                            o2::soa::Index<>, collision::PosZ, collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, emmltrack::HadronicRate,
                            fwdtrack::TrackType, fwdtrack::Pt, fwdtrack::Eta, fwdtrack::Phi, emmlfwdtrack::Sign,
                            fwdtrack::FwdDcaX, fwdtrack::FwdDcaY,
                            emmlfwdtrack::PtMatchedMCHMID, emmlfwdtrack::EtaMatchedMCHMID, emmlfwdtrack::PhiMatchedMCHMID,
                            emmlfwdtrack::XMatchedMCHMIDatMP, emmlfwdtrack::YMatchedMCHMIDatMP,
                            emmlfwdtrack::XMatchedMFTatMP, emmlfwdtrack::YMatchedMFTatMP,
                            fwdtrack::NClusters, fwdtrack::PDca, fwdtrack::RAtAbsorberEnd,
                            fwdtrack::Chi2, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT,
                            // fwdtrack::MCHBitMap, fwdtrack::MIDBitMap, fwdtrack::MIDBoards,
                            fwdtrack::MFTClusterSizesAndTrackFlags, emmlfwdtrack::Chi2MFT, emmlfwdtrack::NClustersMFT, mcparticle::PdgCode, emmlfwdtrack::IsPrimary, emmlfwdtrack::IsCorrectMatchMFTMCH);

using EMFwdTracksForML = EMFwdTracksForML_000;
// iterators
using EMFwdTrackForML = EMFwdTracksForML::iterator;

} // namespace o2::aod

#endif // PWGEM_DILEPTON_DATAMODEL_LMEEMLTABLES_H_
