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
/// \brief Table definitions for gamma-jet analyses
///
/// \author Florian Jonas <florian.jonas@cern.ch>

#ifndef PWGJE_DATAMODEL_GAMMAJETANALYSISTREE_H_
#define PWGJE_DATAMODEL_GAMMAJETANALYSISTREE_H_

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/AnalysisDataModel.h"

namespace o2::aod::gjanalysis
{
enum class ClusterOrigin {
  kUnknown = 0,
  kPhoton, // dominant amount of energy from the cluster is from a photon
  kPromptPhoton,
  kDirectPromptPhoton,
  kFragmentationPhoton,
  kDecayPhoton,     // the particle that produced the cluster is a decay product
  kDecayPhotonPi0,  // the cluster was produced by a pi0 decay
  kDecayPhotonEta,  // the cluster was produced by a eta decay
  kMergedPi0,       // the cluster was produced by a merged pi0, i.e. two photons contribute to the cluster that both come from pi0 decay
  kMergedEta,       // the cluster was produced by a merged eta, i.e. two photons contribute to the cluster that both come from eta decay
  kConvertedPhoton, // the cluster was produced by a converted photon, i.e. a photon that converted to an electron-positron pair and one of the electrons was detected in the cluster
};
enum class ParticleOrigin {
  kUnknown = 0,
  kPromptPhoton,
  kDirectPromptPhoton,
  kFragmentationPhoton,
  kDecayPhoton,
  kDecayPhotonPi0,
  kDecayPhotonEta,
  kDecayPhotonOther,
  kPi0
};
} // namespace o2::aod::gjanalysis
namespace o2::aod
{

// Collision level information
namespace gjevent
{ //! event index
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(Rho, rho, float);
DECLARE_SOA_COLUMN(EventSel, eventSel, uint8_t);
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);
DECLARE_SOA_BITMAP_COLUMN(Alias, alias, 32);
} // namespace gjevent
DECLARE_SOA_TABLE(GjEvents, "AOD", "GJEVENT", o2::soa::Index<>, gjevent::Multiplicity, gjevent::Centrality, gjevent::Rho, gjevent::EventSel, gjevent::Occupancy, gjevent::Alias)

using GjEvent = GjEvents::iterator;

// Information about the MC collision that was matched to the reconstructed collision
namespace gjmcevent
{
DECLARE_SOA_INDEX_COLUMN(GjEvent, gjevent);
DECLARE_SOA_COLUMN(Weight, weight, double);
DECLARE_SOA_COLUMN(Rho, rho, float);                              // gen level rho
DECLARE_SOA_COLUMN(IsMultipleAssigned, isMultipleAssigned, bool); // if the corresponding MC collision matched to this rec collision was also matched to other rec collisions (allows to skip those on analysis level   )
} // namespace gjmcevent
DECLARE_SOA_TABLE(GjMCEvents, "AOD", "GJMCEVENT", gjmcevent::GjEventId, gjmcevent::Weight, gjmcevent::Rho, gjmcevent::IsMultipleAssigned)
// Information about EMCal clusters
namespace gjgamma
{
DECLARE_SOA_INDEX_COLUMN(GjEvent, gjevent);                            //! event index
DECLARE_SOA_COLUMN(Energy, energy, float);                             //! cluster energy (GeV)
DECLARE_SOA_COLUMN(Definition, definition, int);                       //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_COLUMN(Eta, eta, float);                                   //! cluster pseudorapidity (calculated using vertex)
DECLARE_SOA_COLUMN(Phi, phi, float);                                   //! cluster azimuthal angle (calculated using vertex)
DECLARE_SOA_COLUMN(M02, m02, float);                                   //! shower shape long axis
DECLARE_SOA_COLUMN(M20, m20, float);                                   //! shower shape short axis
DECLARE_SOA_COLUMN(NCells, nCells, ushort);                            //! number of cells in cluster
DECLARE_SOA_COLUMN(Time, time, float);                                 //! cluster time (ns)
DECLARE_SOA_COLUMN(IsExotic, isExotic, bool);                          //! flag to mark cluster as exotic
DECLARE_SOA_COLUMN(DistanceToBadChannel, distanceToBadChannel, float); //! distance to bad channel
DECLARE_SOA_COLUMN(NLM, nlm, ushort);                                  //! number of local maxima
DECLARE_SOA_COLUMN(IsoRaw, isoraw, float);                             //! isolation in cone not corrected for Rho
DECLARE_SOA_COLUMN(PerpConeRho, perpconerho, float);                   //! rho in perpendicular cone
DECLARE_SOA_COLUMN(TMdeltaPhi, tmdeltaphi, float);                     //! delta phi between cluster and closest match
DECLARE_SOA_COLUMN(TMdeltaEta, tmdeltaeta, float);                     //! delta eta between cluster and closest match
DECLARE_SOA_COLUMN(TMtrackP, tmtrackp, float);                         //! track momentum of closest match, -1 if no match found
} // namespace gjgamma
DECLARE_SOA_TABLE(GjGammas, "AOD", "GJGAMMA",
                  gjgamma::GjEventId, gjgamma::Energy, gjgamma::Definition, gjgamma::Eta, gjgamma::Phi, gjgamma::M02, gjgamma::M20, gjgamma::NCells, gjgamma::Time, gjgamma::IsExotic, gjgamma::DistanceToBadChannel, gjgamma::NLM, gjgamma::IsoRaw, gjgamma::PerpConeRho, gjgamma::TMdeltaPhi, gjgamma::TMdeltaEta, gjgamma::TMtrackP)

// MC information for reconstructed EMCal clusters
namespace gjgammamcinfo
{
DECLARE_SOA_COLUMN(Origin, origin, uint16_t);
DECLARE_SOA_COLUMN(LeadingEnergyFraction, leadingEnergyFraction, float); // fraction of energy from the leading MC particle
} // namespace gjgammamcinfo
DECLARE_SOA_TABLE(GjGammaMCInfos, "AOD", "GJGAMMAMCINFO", gjgamma::GjEventId, gjgammamcinfo::Origin, gjgammamcinfo::LeadingEnergyFraction)

// Generator level particle information from the MC collision that was matched to the reconstructed collision
namespace gjmcparticle
{
DECLARE_SOA_INDEX_COLUMN(GjEvent, gjevent);
DECLARE_SOA_COLUMN(Energy, energy, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, ushort);        // TODO also add smoe origin of particle? maybe only save original pi0 and eta and photon (not decay photons)
DECLARE_SOA_COLUMN(MCIsolation, mcIsolation, float); // isolation in cone on mc gen level
DECLARE_SOA_COLUMN(Origin, origin, uint16_t);        // origin of particle
} // namespace gjmcparticle
DECLARE_SOA_TABLE(GjMCParticles, "AOD", "GJMCPARTICLE", gjmcparticle::GjEventId, gjmcparticle::Energy, gjmcparticle::Eta, gjmcparticle::Phi, gjmcparticle::Pt, gjmcparticle::PdgCode, gjmcparticle::MCIsolation, gjmcparticle::Origin)

// Reconstructed charged jet information
namespace gjchjet
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Radius, radius, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Area, area, float);
DECLARE_SOA_COLUMN(LeadingTrackPt, leadingtrackpt, float);
DECLARE_SOA_COLUMN(PerpConeRho, perpconerho, float);
DECLARE_SOA_COLUMN(NConstituents, nConstituents, ushort);
} // namespace gjchjet
DECLARE_SOA_TABLE(GjChargedJets, "AOD", "GJCHJET", gjgamma::GjEventId, gjchjet::Pt, gjchjet::Eta, gjchjet::Phi, gjchjet::Radius, gjchjet::Energy, gjchjet::Mass, gjchjet::Area, gjchjet::LeadingTrackPt, gjchjet::PerpConeRho, gjchjet::NConstituents)

using GjChargedJet = GjChargedJets::iterator;

// Jet substructure information (vectors stored per jet)
namespace gjjetsubstructure
{
DECLARE_SOA_COLUMN(EnergyMother, energyMother, std::vector<float>); //! energy of mother subjet at each splitting
DECLARE_SOA_COLUMN(PtLeading, ptLeading, std::vector<float>);       //! pt of leading subjet at each splitting
DECLARE_SOA_COLUMN(PtSubLeading, ptSubLeading, std::vector<float>); //! pt of subleading subjet at each splitting
DECLARE_SOA_COLUMN(Theta, theta, std::vector<float>);               //! opening angle theta at each splitting
} // namespace gjjetsubstructure
DECLARE_SOA_TABLE(GjJetSubstructures, "AOD", "GJJETSUBSTR",
                  gjgamma::GjEventId,
                  gjjetsubstructure::EnergyMother,
                  gjjetsubstructure::PtLeading,
                  gjjetsubstructure::PtSubLeading,
                  gjjetsubstructure::Theta)

// MC information for reconstructed charged jet
namespace gjchjetmcinfo
{
DECLARE_SOA_COLUMN(MatchedJetIndexGeo, matchedJetIndexGeo, int);
DECLARE_SOA_COLUMN(MatchedJetIndexPt, matchedJetIndexPt, int);
} // namespace gjchjetmcinfo
DECLARE_SOA_TABLE(GjChJetMCInfos, "AOD", "GJCHJETMCINFO", gjgamma::GjEventId, gjchjetmcinfo::MatchedJetIndexGeo, gjchjetmcinfo::MatchedJetIndexPt)

// MC information for generator level jets of associated MC collision
namespace gjmcjet
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Radius, radius, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Area, area, float);
DECLARE_SOA_COLUMN(PerpConeRho, perpconerho, float);
} // namespace gjmcjet
DECLARE_SOA_TABLE(GjMCJets, "AOD", "GJMCJET", gjgamma::GjEventId, gjmcjet::Pt, gjmcjet::Eta, gjmcjet::Phi, gjmcjet::Radius, gjmcjet::Energy, gjmcjet::Mass, gjmcjet::Area, gjmcjet::PerpConeRho)

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_GAMMAJETANALYSISTREE_H_
