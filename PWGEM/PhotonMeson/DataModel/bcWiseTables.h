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
///
/// \file bcWiseTables.h
///
/// \brief This header provides the table definitions to store very lightweight EMCal clusters per BC
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
///

#ifndef PWGEM_PHOTONMESON_DATAMODEL_BCWISETABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_BCWISETABLES_H_

#include "Framework/AnalysisDataModel.h"

#include <limits>

namespace o2::aod
{

namespace emdownscaling
{
enum Observable {
  kDefinition,
  kEnergy,
  kEta,
  kPhi,
  kNCells,
  kM02,
  kTime,
  kCent,
  kZVtx,
  kFT0Amp,
  kpT,
  kMu,
  kTimeSinceSOF,
  nObservables
};

// Values in tables are stored in downscaled format to save disk space
const float downscalingFactors[nObservables]{
  1E0,   // Cluster definition
  1E3,   // Cluster energy
  1E4,   // Cluster eta
  1E4,   // Cluster phi
  1E0,   // Number of cells
  1E4,   // M02
  1E2,   // Cluster time
  2E0,   // FT0 centrality
  1E3,   // Z-vertex position
  1E-1,  // FT0M amplitude
  1E3,   // MC pi0 pt
  1E5,   // Mu
  5E-1}; // Time since start of fill (since ADJUST decleration)
} // namespace emdownscaling

namespace bcwisebc
{
DECLARE_SOA_COLUMN(HasFT0, hasFT0, bool);                                //! has_foundFT0()
DECLARE_SOA_COLUMN(HasTVX, hasTVX, bool);                                //! has the TVX trigger flag
DECLARE_SOA_COLUMN(HaskTVXinEMC, haskTVXinEMC, bool);                    //! kTVXinEMC
DECLARE_SOA_COLUMN(HasEMCCell, hasEMCCell, bool);                        //! at least one EMCal cell in the BC
DECLARE_SOA_COLUMN(StoredFT0CCentrality, storedFt0CCentrality, uint8_t); //! FT0C centrality (0-100) (x2)
DECLARE_SOA_COLUMN(StoredFT0MCentrality, storedFT0MCentrality, uint8_t); //! FT0M centrality (0-100) (x2)
DECLARE_SOA_COLUMN(StoredFT0MAmplitude, storedFT0MAmplitude, uint16_t);  //! ft0a+c amplitude
DECLARE_SOA_COLUMN(StoredMu, storedMu, uint16_t);                        //! probability of TVX collision per BC (x1000)
DECLARE_SOA_COLUMN(StoredTimeSinceSOF, storedTimeSinceSOF, uint16_t);    //! time since decreation of ADJUST in seconds (x2)

DECLARE_SOA_DYNAMIC_COLUMN(FT0CCentrality, ft0cCentrality, [](uint8_t storedft0ccentrality) -> float { return std::nextafter(storedft0ccentrality / emdownscaling::downscalingFactors[emdownscaling::kCent], std::numeric_limits<float>::infinity()); });  //! Centrality (0-100)
DECLARE_SOA_DYNAMIC_COLUMN(FT0MCentrality, ft0mCentrality, [](uint8_t storedft0mcentrality) -> float { return std::nextafter(storedft0mcentrality / emdownscaling::downscalingFactors[emdownscaling::kCent], std::numeric_limits<float>::infinity()); });  //! Centrality (0-100)
DECLARE_SOA_DYNAMIC_COLUMN(FT0MAmplitude, ft0Amplitude, [](uint16_t storedFT0MAmplitude) -> float { return std::nextafter(storedFT0MAmplitude / emdownscaling::downscalingFactors[emdownscaling::kFT0Amp], std::numeric_limits<float>::infinity()); });    //! FT0M amplitude
DECLARE_SOA_DYNAMIC_COLUMN(Mu, mu, [](uint16_t storedMu) -> float { return std::nextafter(storedMu / emdownscaling::downscalingFactors[emdownscaling::kMu], std::numeric_limits<float>::infinity()); });                                                   //! probability of TVX collision per BC
DECLARE_SOA_DYNAMIC_COLUMN(TimeSinceSOF, timeSinceSOF, [](uint16_t storedTimeSinceSOF) -> float { return std::nextafter(storedTimeSinceSOF / emdownscaling::downscalingFactors[emdownscaling::kTimeSinceSOF], std::numeric_limits<float>::infinity()); }); //! probability of TVX collision per BC
} // namespace bcwisebc
DECLARE_SOA_TABLE(BCWiseBCs, "AOD", "BCWISEBC", //! table of bc wise centrality estimation and event selection input
                  o2::soa::Index<>, bcwisebc::HasFT0, bcwisebc::HasTVX, bcwisebc::HaskTVXinEMC, bcwisebc::HasEMCCell, bcwisebc::StoredFT0CCentrality, bcwisebc::StoredFT0MCentrality,
                  bcwisebc::StoredFT0MAmplitude, bcwisebc::StoredMu, bcwisebc::StoredTimeSinceSOF, bcwisebc::FT0CCentrality<bcwisebc::StoredFT0CCentrality>, bcwisebc::FT0MCentrality<bcwisebc::StoredFT0MCentrality>, bcwisebc::FT0MAmplitude<bcwisebc::StoredFT0MAmplitude>, bcwisebc::Mu<bcwisebc::StoredMu>, bcwisebc::TimeSinceSOF<bcwisebc::StoredTimeSinceSOF>);

DECLARE_SOA_INDEX_COLUMN(BCWiseBC, bcWiseBC); //! bunch crossing ID used as index

namespace bcwisecollision
{
DECLARE_SOA_COLUMN(StoredZVtx, storedZVtx, int16_t); //! Z-vertex position (x1000)

DECLARE_SOA_DYNAMIC_COLUMN(ZVtx, zVtx, [](int16_t storedzvtx) -> float { return storedzvtx / emdownscaling::downscalingFactors[emdownscaling::kZVtx]; }); //! Z-Vertex
} // namespace bcwisecollision
DECLARE_SOA_TABLE(BCWiseCollisions, "AOD", "BCWISECOLL", //! table of skimmed EMCal clusters
                  o2::soa::Index<>, BCWiseBCId, bcwisecollision::StoredZVtx, bcwisecollision::ZVtx<bcwisecollision::StoredZVtx>);

namespace bcwisecluster
{
DECLARE_SOA_COLUMN(StoredDefinition, storedDefinition, int8_t); //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_COLUMN(StoredE, storedE, int16_t);                  //! cluster energy (1 MeV -> Maximum cluster energy of ~32 GeV)
DECLARE_SOA_COLUMN(StoredEta, storedEta, int16_t);              //! cluster pseudorapidity (x10,000)
DECLARE_SOA_COLUMN(StoredPhi, storedPhi, uint16_t);             //! cluster azimuthal angle (x10 000) from 0 to 2pi
DECLARE_SOA_COLUMN(StoredNCells, storedNCells, int8_t);         //! number of cells in cluster
DECLARE_SOA_COLUMN(StoredM02, storedM02, int16_t);              //! shower shape long axis (x10 000)
DECLARE_SOA_COLUMN(StoredTime, storedTime, int16_t);            //! cluster time (10 ps resolution)
DECLARE_SOA_COLUMN(StoredIsExotic, storedIsExotic, bool);       //! flag to mark cluster as exotic

DECLARE_SOA_DYNAMIC_COLUMN(Definition, definition, [](int8_t storedDefinition) -> int8_t { return storedDefinition; });                                                                                           //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_DYNAMIC_COLUMN(E, e, [](int16_t storedE) -> float { return std::nextafter(storedE / emdownscaling::downscalingFactors[emdownscaling::kEnergy], std::numeric_limits<float>::infinity()); });           //! cluster energy (GeV)
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](int16_t storedEta) -> float { return std::nextafter(storedEta / emdownscaling::downscalingFactors[emdownscaling::kEta], std::numeric_limits<float>::infinity()); });      //! cluster pseudorapidity
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](uint16_t storedPhi) -> float { return std::nextafter(storedPhi / emdownscaling::downscalingFactors[emdownscaling::kPhi], std::numeric_limits<float>::infinity()); });     //! cluster azimuthal angle (0 to 2pi)
DECLARE_SOA_DYNAMIC_COLUMN(NCells, nCells, [](int16_t storedNCells) -> int16_t { return storedNCells; });                                                                                                         //! number of cells in cluster
DECLARE_SOA_DYNAMIC_COLUMN(M02, m02, [](int16_t storedM02) -> float { return std::nextafter(storedM02 / emdownscaling::downscalingFactors[emdownscaling::kM02], std::numeric_limits<float>::infinity()); });      //! shower shape long axis
DECLARE_SOA_DYNAMIC_COLUMN(Time, time, [](int16_t storedTime) -> float { return std::nextafter(storedTime / emdownscaling::downscalingFactors[emdownscaling::kTime], std::numeric_limits<float>::infinity()); }); //! cluster time (ns)
DECLARE_SOA_DYNAMIC_COLUMN(IsExotic, isExotic, [](bool storedIsExotic) -> bool { return storedIsExotic; });                                                                                                       //! flag to mark cluster as exotic

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float storedE, float storedEta) -> float { return storedE / emdownscaling::downscalingFactors[emdownscaling::kEnergy] / std::cosh(storedEta / emdownscaling::downscalingFactors[emdownscaling::kEta]); }); //! cluster pt, assuming m=0 (photons)
} // namespace bcwisecluster

DECLARE_SOA_TABLE(BCWiseClusters, "AOD", "BCWISECLUSTER", //! table of skimmed EMCal clusters
                  o2::soa::Index<>, BCWiseBCId, bcwisecluster::StoredDefinition, bcwisecluster::StoredE, bcwisecluster::StoredEta, bcwisecluster::StoredPhi, bcwisecluster::StoredNCells, bcwisecluster::StoredM02, bcwisecluster::StoredTime, bcwisecluster::StoredIsExotic,
                  bcwisecluster::Definition<bcwisecluster::StoredDefinition>, bcwisecluster::E<bcwisecluster::StoredE>, bcwisecluster::Eta<bcwisecluster::StoredEta>, bcwisecluster::Phi<bcwisecluster::StoredPhi>, bcwisecluster::NCells<bcwisecluster::StoredNCells>, bcwisecluster::M02<bcwisecluster::StoredM02>, bcwisecluster::Time<bcwisecluster::StoredTime>, bcwisecluster::IsExotic<bcwisecluster::StoredIsExotic>,
                  bcwisecluster::Pt<bcwisecluster::StoredE, bcwisecluster::StoredEta>);

namespace bcwisemcmesons
{
DECLARE_SOA_COLUMN(StoredPt, storedPt, uint16_t); //! Transverse momentum of generated pi0 (1 MeV -> Maximum pi0 pT of ~65 GeV)
DECLARE_SOA_COLUMN(IsAccepted, isAccepted, bool); //! Both decay photons are within the EMCal acceptance
DECLARE_SOA_COLUMN(IsPrimary, isPrimary, bool);   //! mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator()
DECLARE_SOA_COLUMN(IsFromWD, isFromWD, bool);     //! Pi0 from a weak decay according to pwgem::photonmeson::utils::mcutil::IsFromWD

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](uint16_t storedpt) -> float { return std::nextafter(storedpt / emdownscaling::downscalingFactors[emdownscaling::kpT], std::numeric_limits<float>::infinity()); }); //! pT of pi0 (GeV)
} // namespace bcwisemcmesons

DECLARE_SOA_TABLE(BCWiseMCPi0s, "AOD", "BCWISEMCPI0", //! table of pi0s on MC level
                  o2::soa::Index<>, BCWiseBCId, bcwisemcmesons::StoredPt, bcwisemcmesons::IsAccepted, bcwisemcmesons::IsPrimary, bcwisemcmesons::IsFromWD,
                  bcwisemcmesons::Pt<bcwisemcmesons::StoredPt>);
DECLARE_SOA_TABLE(BCWiseMCEtas, "AOD", "BCWISEMCETA", //! table of eta mesons on MC level
                  o2::soa::Index<>, BCWiseBCId, bcwisemcmesons::StoredPt, bcwisemcmesons::IsAccepted, bcwisemcmesons::IsPrimary, bcwisemcmesons::IsFromWD,
                  bcwisemcmesons::Pt<bcwisemcmesons::StoredPt>);

namespace bcwisemccluster
{
DECLARE_SOA_COLUMN(MesonID, mesonID, int32_t);          //! Index of the mother mesom (-1 if not from a pi0 or eta)
DECLARE_SOA_COLUMN(IsEta, isEta, bool);                 //! Boolean flag to indicate if the cluster is from an eta meson, otherwise it is from a pi0
DECLARE_SOA_COLUMN(StoredTrueE, storedTrueE, uint16_t); //! energy of cluster inducing particle (1 MeV -> Maximum cluster energy of ~65 GeV)

DECLARE_SOA_DYNAMIC_COLUMN(TrueE, trueE, [](uint16_t storedTrueE) -> float { return std::nextafter(storedTrueE / emdownscaling::downscalingFactors[emdownscaling::kEnergy], std::numeric_limits<float>::infinity()); }); //! energy of cluster inducing particle (GeV)
} // namespace bcwisemccluster

DECLARE_SOA_TABLE(BCWiseMCClusters, "AOD", "BCWISEMCCLS", //! table of MC information for clusters -> To be joined with the cluster table
                  o2::soa::Index<>, BCWiseBCId, bcwisemccluster::MesonID, bcwisemccluster::IsEta, bcwisemccluster::StoredTrueE,
                  bcwisemccluster::TrueE<bcwisemccluster::StoredTrueE>);

} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_BCWISETABLES_H_
