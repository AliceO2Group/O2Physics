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

#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

namespace o2::aod
{

namespace gjevent
{ // TODO add rho                          //! event index
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(Rho, rho, float);
DECLARE_SOA_COLUMN(EventSel, eventSel, uint8_t);
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);
DECLARE_SOA_BITMAP_COLUMN(Alias, alias, 32);
} // namespace gjevent
DECLARE_SOA_TABLE(GjEvents, "AOD", "GJEVENT", o2::soa::Index<>, gjevent::Multiplicity, gjevent::Centrality, gjevent::Rho, gjevent::EventSel, gjevent::Occupancy, gjevent::Alias)

using GjEvent = GjEvents::iterator;

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
DECLARE_SOA_COLUMN(IsoRaw, isoraw, ushort);                            //! isolation in cone not corrected for Rho
DECLARE_SOA_COLUMN(PerpConeRho, perpconerho, float);                   //! rho in perpendicular cone
DECLARE_SOA_COLUMN(TMdeltaPhi, tmdeltaphi, float);                     //! delta phi between cluster and closest match
DECLARE_SOA_COLUMN(TMdeltaEta, tmdeltaeta, float);                     //! delta eta between cluster and closest match
DECLARE_SOA_COLUMN(TMtrackP, tmtrackp, float);                         //! track momentum of closest match, -1 if no match found
} // namespace gjgamma
DECLARE_SOA_TABLE(GjGammas, "AOD", "GJGAMMA",
                  gjgamma::GjEventId, gjgamma::Energy, gjgamma::Definition, gjgamma::Eta, gjgamma::Phi, gjgamma::M02, gjgamma::M20, gjgamma::NCells, gjgamma::Time, gjgamma::IsExotic, gjgamma::DistanceToBadChannel, gjgamma::NLM, gjgamma::IsoRaw, gjgamma::PerpConeRho, gjgamma::TMdeltaPhi, gjgamma::TMdeltaEta, gjgamma::TMtrackP)
namespace gjchjet
{
DECLARE_SOA_INDEX_COLUMN(GjEvent, gjevent);
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
DECLARE_SOA_TABLE(GjChargedJets, "AOD", "GJCHJET", gjchjet::GjEventId, gjchjet::Pt, gjchjet::Eta, gjchjet::Phi, gjchjet::Radius, gjchjet::Energy, gjchjet::Mass, gjchjet::Area, gjchjet::LeadingTrackPt, gjchjet::PerpConeRho, gjchjet::NConstituents)
} // namespace o2::aod

#endif // PWGJE_DATAMODEL_GAMMAJETANALYSISTREE_H_
