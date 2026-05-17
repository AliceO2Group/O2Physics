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

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

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
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                             //!
DECLARE_SOA_COLUMN(EtaGen, etaGen, float);                           //!
DECLARE_SOA_COLUMN(PhiGen, phiGen, float);                           //!

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
// DECLARE_SOA_COLUMN(Signed1PtMFTMCHMIDatPV, signed1PtMFTMCHMIDatPV, float); //! pt of MCH-MID track in MFT-MCH-MID track at PV
// DECLARE_SOA_COLUMN(EtaMFTMCHMIDatPV, etaMFTMCHMIDatPV, float);             //! eta of MCH-MID track in MFT-MCH-MID track at PV
// DECLARE_SOA_COLUMN(PhiMFTMCHMIDatPV, phiMFTMCHMIDatPV, float);             //! phi of MCH-MID track in MFT-MCH-MID track at PV
//
// DECLARE_SOA_COLUMN(Signed1PtMCHMIDatPV, signed1PtMCHMIDatPV, float); //! pt of MCH-MID track in MFT-MCH-MID track at PV
// DECLARE_SOA_COLUMN(EtaMCHMIDatPV, etaMCHMIDatPV, float);             //! eta of MCH-MID track in MFT-MCH-MID track at PV
// DECLARE_SOA_COLUMN(PhiMCHMIDatPV, phiMCHMIDatPV, float);             //! phi of MCH-MID track in MFT-MCH-MID track at PV

DECLARE_SOA_COLUMN(Signed1PtMCHMIDatMP, signed1PtMCHMIDatMP, float); //! pt of MCH-MID track in MFT-MCH-MID track at MP
DECLARE_SOA_COLUMN(TglMCHMIDatMP, tglMCHMIDatMP, float);             //! tgl of MCH-MID track in MFT-MCH-MID track at MP
DECLARE_SOA_COLUMN(PhiMCHMIDatMP, phiMCHMIDatMP, float);             //! phi of MCH-MID track in MFT-MCH-MID track at MP

DECLARE_SOA_COLUMN(Signed1PtErrMCHMIDatMP, signed1PtErrMCHMIDatMP, float); //! pt of MCH-MID track in MFT-MCH-MID track at MP
DECLARE_SOA_COLUMN(TglErrMCHMIDatMP, tglErrMCHMIDatMP, float);             //! tgl of MCH-MID track in MFT-MCH-MID track at MP
DECLARE_SOA_COLUMN(PhiErrMCHMIDatMP, phiErrMCHMIDatMP, float);             //! phi of MCH-MID track in MFT-MCH-MID track at MP

DECLARE_SOA_COLUMN(Signed1PtMFTatMP, signed1PtMFTatMP, float); //! pt of MCH-MID track in MFT-MCH-MID track at MP
DECLARE_SOA_COLUMN(TglMFTatMP, tglMFTatMP, float);             //! tgl of MCH-MID track in MFT-MCH-MID track at MP
DECLARE_SOA_COLUMN(PhiMFTatMP, phiMFTatMP, float);             //! phi of MCH-MID track in MFT-MCH-MID track at MP

DECLARE_SOA_COLUMN(Signed1PtErrMFTatMP, signed1PtErrMFTatMP, float); //! pt of MCH-MID track in MFT-MCH-MID track at MP
DECLARE_SOA_COLUMN(TglErrMFTatMP, tglErrMFTatMP, float);             //! tgl of MCH-MID track in MFT-MCH-MID track at MP
DECLARE_SOA_COLUMN(PhiErrMFTatMP, phiErrMFTatMP, float);             //! phi of MCH-MID track in MFT-MCH-MID track at MP

DECLARE_SOA_COLUMN(XMCHMIDatMP, xMCHMIDatMP, float);       //! x of MCH-MID track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(YMCHMIDatMP, yMCHMIDatMP, float);       //! y of MCH-MID track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(XErrMCHMIDatMP, xErrMCHMIDatMP, float); //! x error of MCH-MID track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(YErrMCHMIDatMP, yErrMCHMIDatMP, float); //! y error of MCH-MID track in MFT-MCH-MID track at matching plane

DECLARE_SOA_COLUMN(XMFTatMP, xMFTatMP, float);       //! x of MFT track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(YMFTatMP, yMFTatMP, float);       //! y of MFT track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(XErrMFTatMP, xErrMFTatMP, float); //! x error of MFT track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(YErrMFTatMP, yErrMFTatMP, float); //! y error of MFT track in MFT-MCH-MID track at matching plane

DECLARE_SOA_COLUMN(Chi2MFT, chi2MFT, float);                //! chi2/ndf of MFT track
DECLARE_SOA_COLUMN(Chi2MCHMID, chi2MCHMID, float);          //! chi2/ndf of MCH-MID track
DECLARE_SOA_COLUMN(Chi2MFTMCHMID, chi2MFTMCHMID, float);    //! chi2/ndf of MFT-MCH-MID track
DECLARE_SOA_COLUMN(NClustersMFT, nClustersMFT, uint8_t);    //!
DECLARE_SOA_COLUMN(IsPrimaryMFT, isPrimaryMFT, bool);       //!
DECLARE_SOA_COLUMN(IsPrimaryMCHMID, isPrimaryMCHMID, bool); //!
DECLARE_SOA_COLUMN(IsCorrectMatch, isCorrectMatch, bool);   //!
DECLARE_SOA_COLUMN(PdgCodeMFT, pdgCodeMFT, int);            //!
DECLARE_SOA_COLUMN(PdgCodeMCHMID, pdgCodeMCHMID, int);      //!
DECLARE_SOA_COLUMN(MatchMCHTrackId, mchTrackId, int);       //!
DECLARE_SOA_COLUMN(DFId, dfId, uint64_t);                   //!

DECLARE_SOA_COLUMN(MultMFT, multMFT, uint16_t); //! number of MFTsa tracks per collision
} // namespace emmlfwdtrack

DECLARE_SOA_TABLE(EMFwdTracksForML, "AOD", "EMFWDTRKML", //!
                  o2::soa::Index<>, collision::PosZ, /*collision::NumContrib,*/ mult::MultFT0C, /*evsel::NumTracksInTimeRange,*/ evsel::SumAmpFT0CInTimeRange, emmltrack::HadronicRate, emmlfwdtrack::MultMFT,

                  emmlfwdtrack::Signed1PtMFTatMP, emmlfwdtrack::TglMFTatMP, emmlfwdtrack::PhiMFTatMP,
                  emmlfwdtrack::XMFTatMP, emmlfwdtrack::YMFTatMP,

                  emmlfwdtrack::Signed1PtMCHMIDatMP, emmlfwdtrack::TglMCHMIDatMP, emmlfwdtrack::PhiMCHMIDatMP,
                  emmlfwdtrack::XMCHMIDatMP, emmlfwdtrack::YMCHMIDatMP,

                  fwdtrack::Chi2MatchMCHMFT,
                  emmlfwdtrack::PdgCodeMFT, emmlfwdtrack::IsPrimaryMFT,
                  emmlfwdtrack::PdgCodeMCHMID, emmlfwdtrack::IsPrimaryMCHMID,
                  emmlfwdtrack::IsCorrectMatch, emmlfwdtrack::MatchMCHTrackId, emmlfwdtrack::DFId);

// iterators
using EMFwdTrackForML = EMFwdTracksForML::iterator;

DECLARE_SOA_TABLE(EMFwdTrackErrsForML, "AOD", "EMFWDTRKERRML", //! Joinable with EMFwdTracksForML
                  /*emmlfwdtrack::Signed1PtErrMFTatMP,*/ emmlfwdtrack::TglErrMFTatMP, emmlfwdtrack::PhiErrMFTatMP,
                  emmlfwdtrack::XErrMFTatMP, emmlfwdtrack::YErrMFTatMP,
                  /*emmlfwdtrack::Signed1PtErrMCHMIDatMP,*/ emmlfwdtrack::TglErrMCHMIDatMP, emmlfwdtrack::PhiErrMCHMIDatMP,
                  emmlfwdtrack::XErrMCHMIDatMP, emmlfwdtrack::YErrMCHMIDatMP);

// iterators
using EMFwdTrackErrForML = EMFwdTrackErrsForML::iterator;

// for SemiCharmTag at midrapidity, only electrons
namespace emmlevent
{
DECLARE_SOA_COLUMN(SubGeneratorId, subGeneratorId, int); //! sub generator Id of mc collision
} // namespace emmlevent
namespace emmltrack
{
DECLARE_SOA_COLUMN(IsMotherFromBeauty, isMotherFromBeauty, bool); //! is b quark included in decay history
DECLARE_SOA_COLUMN(Signed1PtL, signedPtL, float);                 //! sign/pT of lepton
DECLARE_SOA_COLUMN(EtaL, etaL, float);                            //! eta of lepton
DECLARE_SOA_COLUMN(PhiL, phiL, float);                            //! phi of lepton
DECLARE_SOA_COLUMN(ImpParXYL, impParXYL, float);                  //! impact parameter for lepton in XY plane
DECLARE_SOA_COLUMN(ImpParZL, impParZL, float);                    //! impact parameter for lepton in Z plane
DECLARE_SOA_COLUMN(ImpParCYYL, impParCYYL, float);                //! sigma of impact parameter for lepton in XY
DECLARE_SOA_COLUMN(ImpParCZYL, impParCZYL, float);                //! sigma of impact parameter for lepton, correlaion term
DECLARE_SOA_COLUMN(ImpParCZZL, impParCZZL, float);                //! sigma of impact parameter for lepton in Z
DECLARE_SOA_COLUMN(PdgCodeMother, pdgCodeMother, int);            //! pdg code of mother of lepton
DECLARE_SOA_COLUMN(IsCorrectCollision, isCorrectCollision, bool); //! LH pair is associated to correct collision.
} // namespace emmltrack

DECLARE_SOA_TABLE(EMMLLeptons, "AOD", "EMMLLEPTON", //!
                  o2::soa::Index<>, collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, emmlevent::SubGeneratorId,
                  emmltrack::Signed1PtL, emmltrack::EtaL,
                  emmltrack::ImpParXYL, emmltrack::ImpParZL, emmltrack::ImpParCYYL, emmltrack::ImpParCZYL, emmltrack::ImpParCZZL,
                  emmltrack::IsMotherFromBeauty, emmltrack::PdgCodeMother, emmltrack::IsCorrectCollision);
// iterators
using EMMLLepton = EMMLLeptons::iterator;

namespace emmllhpair
{
DECLARE_SOA_INDEX_COLUMN(EMMLLepton, emmllepton); //! most propable emeventId

DECLARE_SOA_COLUMN(Signed1PtH, signedPtH, float);  //! sign/pT of associated hadron
DECLARE_SOA_COLUMN(PtH, ptH, float);               //! pT of associated hadron
DECLARE_SOA_COLUMN(EtaH, etaH, float);             //! eta of associated hadron
DECLARE_SOA_COLUMN(PhiH, phiH, float);             //! phi of associated hadron
DECLARE_SOA_COLUMN(RapidityV0, rapidityV0, float); //! rapidity of associated V0
DECLARE_SOA_COLUMN(RapidityC, rapidityC, float);   //! rapidity of associated Cascade

DECLARE_SOA_COLUMN(ImpParXYH, impParXYH, float);   //! impact parameter for V0/Cascade in XY plane
DECLARE_SOA_COLUMN(ImpParZH, impParZH, float);     //! impact parameter for V0/Cascade in Z plane
DECLARE_SOA_COLUMN(ImpParCYYH, impParCYYH, float); //! sigma of impact parameter for V0/Cascade in XY
DECLARE_SOA_COLUMN(ImpParCZYH, impParCZYH, float); //! sigma of impact parameter for V0/Cascade, correlaion term
DECLARE_SOA_COLUMN(ImpParCZZH, impParCZZH, float); //! sigma of impact parameter for V0/Cascade in Z

DECLARE_SOA_COLUMN(V0CPA, v0cpa, float);     //! cosPA of V0
DECLARE_SOA_COLUMN(V0CPAXY, v0cpaXY, float); //! cosPA of V0 in XY plane
DECLARE_SOA_COLUMN(V0CPARZ, v0cpaRZ, float); //! cosPA of V0 in XY plane

DECLARE_SOA_COLUMN(CascCPA, casccpa, float);     //! cosPA of Cascade
DECLARE_SOA_COLUMN(CascCPAXY, casccpaXY, float); //! cosPA of Cascade in XY plane
DECLARE_SOA_COLUMN(CascCPARZ, casccpaRZ, float); //! cosPA of Cascade in RZ plane

DECLARE_SOA_COLUMN(V0Type, v0Type, uint8_t);           //! v0 type, 0 = K0S, 1 = Lambda, 2 = AntiLambda
DECLARE_SOA_COLUMN(CascadeType, cascadeType, uint8_t); //! cascade type, 0 = XiMunus, 1 = OmegaMinus

// LH pair variables
DECLARE_SOA_COLUMN(MassLH, massLH, float); //! invariant mass of LH assuming pion
DECLARE_SOA_COLUMN(PtLH, ptLH, float);     //! pT of LH pair
DECLARE_SOA_COLUMN(DcaLH, dcalh, float);   //! DCA between lepton and hadron
DECLARE_SOA_COLUMN(CPA, cpa, float);       //! cosine of pointing angle of LH pair
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);   //! cosine of pointing angle of LH pair in XY
DECLARE_SOA_COLUMN(CPARZ, cpaRZ, float);   //! cosine of pointing angle of LH pair in RZ

DECLARE_SOA_COLUMN(Lxy, lxy, float);         //! decay length of LH pair
DECLARE_SOA_COLUMN(Lz, lz, float);           //! decay length of LH pair
DECLARE_SOA_COLUMN(Lxyz, lxyz, float);       //! decay length of LH pair
DECLARE_SOA_COLUMN(LxyErr, lxyErr, float);   //! decay length resolution of LH pair
DECLARE_SOA_COLUMN(LzErr, lzErr, float);     //! decay length resolution of LH pair
DECLARE_SOA_COLUMN(LxyzErr, lxyzErr, float); //! decay length resolution of LH pair

DECLARE_SOA_COLUMN(ImpParXY, impParXY, float);   //! impact parameter for LH in XY plane
DECLARE_SOA_COLUMN(ImpParZ, impParZ, float);     //! impact parameter for LH in Z plane
DECLARE_SOA_COLUMN(ImpParCYY, impParCYY, float); //! sigma of impact parameter for LH in XY
DECLARE_SOA_COLUMN(ImpParCZY, impParCZY, float); //! sigma of impact parameter for LH, correlation term
DECLARE_SOA_COLUMN(ImpParCZZ, impParCZZ, float); //! sigma of impact parameter for LH in Z

DECLARE_SOA_COLUMN(PdgCodeH, pdgCodeH, int);                    //! pdg code of associated hadron
DECLARE_SOA_COLUMN(PdgCodeIM, pdgCodeIM, int);                  //! pdg code of intermediate hadron from HF hadrons. e.g K*, D*
DECLARE_SOA_COLUMN(FoundCommonMother, foundCommonMother, bool); //! decay length resolution of LH pair
} // namespace emmllhpair

DECLARE_SOA_TABLE(EMMLLHPairs, "AOD", "EMMLLHPAIR", //!
                  emmllhpair::EMMLLeptonId,
                  emmllhpair::Signed1PtH, emmllhpair::EtaH,
                  emmllhpair::ImpParXYH, emmllhpair::ImpParZH, emmllhpair::ImpParCYYH, emmllhpair::ImpParCZYH, emmllhpair::ImpParCZZH,
                  pidtpc::TPCNSigmaPi, pidtof::TOFNSigmaPi,
                  pidtpc::TPCNSigmaKa, pidtof::TOFNSigmaKa,
                  pidtpc::TPCNSigmaPr, pidtof::TOFNSigmaPr,
                  emmllhpair::MassLH, emmllhpair::PtLH, emmllhpair::DcaLH, emmllhpair::CPA, emmllhpair::CPAXY, emmllhpair::CPARZ,
                  emmllhpair::Lxy, emmllhpair::Lz, emmllhpair::Lxyz, emmllhpair::LxyErr, emmllhpair::LzErr, emmllhpair::LxyzErr,
                  emmllhpair::ImpParXY, emmllhpair::ImpParZ, emmllhpair::ImpParCYY, emmllhpair::ImpParCZY, emmllhpair::ImpParCZZ,
                  emmllhpair::PdgCodeH, emmllhpair::PdgCodeIM, emmllhpair::FoundCommonMother);
// iterators
using EMMLLHPair = EMMLLHPairs::iterator;

DECLARE_SOA_TABLE(EMMLLV0Pairs, "AOD", "EMMLLV0PAIR", //!
                  emmllhpair::EMMLLeptonId, emmllhpair::V0Type,
                  emmllhpair::PtH, emmllhpair::RapidityV0,
                  emmllhpair::V0CPA, emmllhpair::V0CPAXY, emmllhpair::V0CPARZ,
                  emmllhpair::ImpParXYH, emmllhpair::ImpParZH, emmllhpair::ImpParCYYH, emmllhpair::ImpParCZYH, emmllhpair::ImpParCZZH,
                  emmllhpair::MassLH, emmllhpair::PtLH, emmllhpair::DcaLH, emmllhpair::CPA, emmllhpair::CPAXY, emmllhpair::CPARZ,
                  emmllhpair::Lxy, emmllhpair::Lz, emmllhpair::Lxyz, emmllhpair::LxyErr, emmllhpair::LzErr, emmllhpair::LxyzErr,
                  emmllhpair::ImpParXY, emmllhpair::ImpParZ, emmllhpair::ImpParCYY, emmllhpair::ImpParCZY, emmllhpair::ImpParCZZ,
                  emmllhpair::PdgCodeH, emmllhpair::PdgCodeIM, emmllhpair::FoundCommonMother);
// iterators
using EMMLLV0Pair = EMMLLV0Pairs::iterator;

DECLARE_SOA_TABLE(EMMLLCascPairs, "AOD", "EMMLLCAPAIR", //!
                  emmllhpair::EMMLLeptonId, emmllhpair::CascadeType,
                  emmllhpair::Signed1PtH, emmllhpair::RapidityC,
                  emmllhpair::CascCPA, emmllhpair::CascCPAXY, emmllhpair::CascCPARZ,
                  emmllhpair::ImpParXYH, emmllhpair::ImpParZH, emmllhpair::ImpParCYYH, emmllhpair::ImpParCZYH, emmllhpair::ImpParCZZH,
                  emmllhpair::MassLH, emmllhpair::PtLH, emmllhpair::DcaLH, emmllhpair::CPA, emmllhpair::CPAXY, emmllhpair::CPARZ,
                  emmllhpair::Lxy, emmllhpair::Lz, emmllhpair::Lxyz, emmllhpair::LxyErr, emmllhpair::LzErr, emmllhpair::LxyzErr,
                  emmllhpair::ImpParXY, emmllhpair::ImpParZ, emmllhpair::ImpParCYY, emmllhpair::ImpParCZY, emmllhpair::ImpParCZZ,
                  emmllhpair::PdgCodeH, emmllhpair::PdgCodeIM, emmllhpair::FoundCommonMother);
// iterators
using EMMLLCascPair = EMMLLCascPairs::iterator;

} // namespace o2::aod

#endif // PWGEM_DILEPTON_DATAMODEL_LMEEMLTABLES_H_
