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

namespace o2::aod
{

namespace emdownscaling
{
enum Observables {
  kDefinition,
  kEnergy,
  kEta,
  kPhi,
  kNCells,
  kM02,
  kTime,
  kFT0MCent,
  kZVtx,
  kFT0Amp,
  kpT,
  nObservables
};

// Values in tables are stored in downscaled format to save disk space
const float downscalingFactors[nObservables]{
  1E0,  // Cluster definition
  1E3,  // Cluster energy
  1E4,  // Cluster eta
  1E4,  // Cluster phi
  1E0,  // Number of cells
  1E4,  // M02
  1E2,  // Cluster time
  2E0,  // FT0M centrality
  1E3,  // Z-vertex position
  1E-1, // FT0M amplitude
  1E3}; // MC pi0 pt
} // namespace emdownscaling

namespace bcwisebc
{
DECLARE_SOA_COLUMN(HasFT0, hasFT0, bool);                               //! has_foundFT0()
DECLARE_SOA_COLUMN(HasTVX, hasTVX, bool);                               //! has the TVX trigger flag
DECLARE_SOA_COLUMN(HaskTVXinEMC, haskTVXinEMC, bool);                   //! kTVXinEMC
DECLARE_SOA_COLUMN(HasEMCCell, hasEMCCell, bool);                       //! at least one EMCal cell in the BC
DECLARE_SOA_COLUMN(HasNoTFROFBorder, hasNoTFROFBorder, bool);           //! not in the TF border or ITS ROF border region
DECLARE_SOA_COLUMN(StoredFT0MAmplitude, storedFT0MAmplitude, uint16_t); //! ft0a+c amplitude

DECLARE_SOA_DYNAMIC_COLUMN(FT0MAmplitude, ft0Amplitude, [](uint16_t storedFT0MAmplitude) -> float { return storedFT0MAmplitude / emdownscaling::downscalingFactors[emdownscaling::kFT0Amp]; }); //! FT0M amplitude
} // namespace bcwisebc
DECLARE_SOA_TABLE(BCWiseBCs, "AOD", "BCWISEBC", //! table of bc wise centrality estimation and event selection input
                  o2::soa::Index<>, bcwisebc::HasFT0, bcwisebc::HasTVX, bcwisebc::HaskTVXinEMC, bcwisebc::HasEMCCell, bcwisebc::HasNoTFROFBorder,
                  bcwisebc::StoredFT0MAmplitude, bcwisebc::FT0MAmplitude<bcwisebc::StoredFT0MAmplitude>);

DECLARE_SOA_INDEX_COLUMN(BCWiseBC, bcWiseBC); //! bunch crossing ID used as index

namespace bcwisecollision
{
DECLARE_SOA_COLUMN(StoredCentrality, storedCentrality, uint8_t); //! FT0M centrality (0-100) (x2)
DECLARE_SOA_COLUMN(StoredZVtx, storedZVtx, int16_t);             //! Z-vertex position (x1000)

DECLARE_SOA_DYNAMIC_COLUMN(Centrality, centrality, [](uint8_t storedcentrality) -> float { return storedcentrality / emdownscaling::downscalingFactors[emdownscaling::kFT0MCent]; }); //! Centrality (0-100)
DECLARE_SOA_DYNAMIC_COLUMN(ZVtx, zVtx, [](uint8_t storedzvtx) -> float { return storedzvtx / emdownscaling::downscalingFactors[emdownscaling::kZVtx]; });                             //! Centrality (0-100)
} // namespace bcwisecollision
DECLARE_SOA_TABLE(BCWiseCollisions, "AOD", "BCWISECOLL", //! table of skimmed EMCal clusters
                  o2::soa::Index<>, BCWiseBCId, bcwisecollision::StoredCentrality, bcwisecollision::StoredZVtx,
                  bcwisecollision::Centrality<bcwisecollision::StoredCentrality>, bcwisecollision::ZVtx<bcwisecollision::StoredZVtx>);

namespace bcwisecluster
{
DECLARE_SOA_COLUMN(StoredDefinition, storedDefinition, uint8_t); //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_COLUMN(StoredE, storedE, uint16_t);                  //! cluster energy (1 MeV -> Maximum cluster energy of ~65 GeV)
DECLARE_SOA_COLUMN(StoredEta, storedEta, int16_t);               //! cluster pseudorapidity (x10,000)
DECLARE_SOA_COLUMN(StoredPhi, storedPhi, uint16_t);              //! cluster azimuthal angle (x10 000) from 0 to 2pi
DECLARE_SOA_COLUMN(StoredNCells, storedNCells, uint8_t);         //! number of cells in cluster
DECLARE_SOA_COLUMN(StoredM02, storedM02, uint16_t);              //! shower shape long axis (x10 000)
DECLARE_SOA_COLUMN(StoredTime, storedTime, int16_t);             //! cluster time (10 ps resolution)
DECLARE_SOA_COLUMN(StoredIsExotic, storedIsExotic, bool);        //! flag to mark cluster as exotic

DECLARE_SOA_DYNAMIC_COLUMN(Definition, definition, [](uint8_t storedDefinition) -> uint8_t { return storedDefinition; });                                 //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_DYNAMIC_COLUMN(E, e, [](uint16_t storedE) -> float { return storedE / emdownscaling::downscalingFactors[emdownscaling::kEnergy]; });          //! cluster energy (GeV)
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](int16_t storedEta) -> float { return storedEta / emdownscaling::downscalingFactors[emdownscaling::kEta]; });      //! cluster pseudorapidity
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](uint16_t storedPhi) -> float { return storedPhi / emdownscaling::downscalingFactors[emdownscaling::kPhi]; });     //! cluster azimuthal angle (0 to 2pi)
DECLARE_SOA_DYNAMIC_COLUMN(NCells, nCells, [](uint16_t storedNCells) -> uint16_t { return storedNCells; });                                               //! number of cells in cluster
DECLARE_SOA_DYNAMIC_COLUMN(M02, m02, [](uint16_t storedM02) -> float { return storedM02 / emdownscaling::downscalingFactors[emdownscaling::kM02]; });     //! shower shape long axis
DECLARE_SOA_DYNAMIC_COLUMN(Time, time, [](int16_t storedTime) -> float { return storedTime / emdownscaling::downscalingFactors[emdownscaling::kTime]; }); //! cluster time (ns)
DECLARE_SOA_DYNAMIC_COLUMN(IsExotic, isExotic, [](bool storedIsExotic) -> bool { return storedIsExotic; });                                               //! flag to mark cluster as exotic

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float storedE, float storedEta) -> float { return storedE / emdownscaling::downscalingFactors[emdownscaling::kEnergy] / std::cosh(storedEta / emdownscaling::downscalingFactors[emdownscaling::kEta]); }); //! cluster pt, assuming m=0 (photons)
} // namespace bcwisecluster

DECLARE_SOA_TABLE(BCWiseClusters, "AOD", "BCWISECLUSTER", //! table of skimmed EMCal clusters
                  o2::soa::Index<>, BCWiseBCId, bcwisecluster::StoredDefinition, bcwisecluster::StoredE, bcwisecluster::StoredEta, bcwisecluster::StoredPhi, bcwisecluster::StoredNCells, bcwisecluster::StoredM02, bcwisecluster::StoredTime, bcwisecluster::StoredIsExotic,
                  bcwisecluster::Definition<bcwisecluster::StoredDefinition>, bcwisecluster::E<bcwisecluster::StoredE>, bcwisecluster::Eta<bcwisecluster::StoredEta>, bcwisecluster::Phi<bcwisecluster::StoredPhi>, bcwisecluster::NCells<bcwisecluster::StoredNCells>, bcwisecluster::M02<bcwisecluster::StoredM02>, bcwisecluster::Time<bcwisecluster::StoredTime>, bcwisecluster::IsExotic<bcwisecluster::StoredIsExotic>,
                  bcwisecluster::Pt<bcwisecluster::StoredE, bcwisecluster::StoredEta>);

namespace bcwisemcpi0s
{
DECLARE_SOA_COLUMN(StoredPt, storedPt, uint16_t); //! Transverse momentum of generated pi0 (1 MeV -> Maximum pi0 pT of ~65 GeV)
DECLARE_SOA_COLUMN(IsAccepted, isAccepted, bool); //! Both decay photons are within the EMCal acceptance
DECLARE_SOA_COLUMN(IsPrimary, isPrimary, bool);   //! mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator()
DECLARE_SOA_COLUMN(IsFromWD, isFromWD, bool);     //! Pi0 from a weak decay according to pwgem::photonmeson::utils::mcutil::IsFromWD

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](uint16_t storedpt) -> float { return storedpt / emdownscaling::downscalingFactors[emdownscaling::kpT]; }); //! pT of pi0 (GeV)
} // namespace bcwisemcpi0s

DECLARE_SOA_TABLE(BCWiseMCPi0s, "AOD", "BCWISEMCPI0", //! table of pi0s on MC level
                  o2::soa::Index<>, BCWiseBCId, bcwisemcpi0s::StoredPt, bcwisemcpi0s::IsAccepted, bcwisemcpi0s::IsPrimary, bcwisemcpi0s::IsFromWD,
                  bcwisemcpi0s::Pt<bcwisemcpi0s::StoredPt>);

namespace bcwisemccluster
{
DECLARE_SOA_COLUMN(Pi0ID, pi0ID, int32_t);              //! Index of the mother pi0 (-1 if not from pi0)
DECLARE_SOA_COLUMN(StoredTrueE, storedTrueE, uint16_t); //! energy of cluster inducing particle (1 MeV -> Maximum cluster energy of ~65 GeV)

DECLARE_SOA_DYNAMIC_COLUMN(TrueE, trueE, [](uint16_t storedTrueE) -> float { return storedTrueE / emdownscaling::downscalingFactors[emdownscaling::kEnergy]; }); //! energy of cluster inducing particle (GeV)
} // namespace bcwisemccluster

DECLARE_SOA_TABLE(BCWiseMCClusters, "AOD", "BCWISEMCCLS", //! table of MC information for clusters -> To be joined with the cluster table
                  o2::soa::Index<>, BCWiseBCId, bcwisemccluster::Pi0ID, bcwisemccluster::StoredTrueE,
                  bcwisemccluster::TrueE<bcwisemccluster::StoredTrueE>);

} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_BCWISETABLES_H_
