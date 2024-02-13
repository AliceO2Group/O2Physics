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

#ifndef PWGMM_MULT_CORE_INCLUDE_HISTOGRAMS_H_
#define PWGMM_MULT_CORE_INCLUDE_HISTOGRAMS_H_
#include "TPDGCode.h"
#include <string_view>

namespace pwgmm::mult
{
// particle species to consider for tracking efficiency
static constexpr std::string_view species[] = {"pi", "p", "e", "K"};
static constexpr std::array<int, 4> speciesIds{kPiPlus, kProton, kElectron, kKPlus};
static constexpr std::string_view prefix = "Tracks/Control/";
static constexpr std::string_view PtGenSuff = "/PtGen";
static constexpr std::string_view PtGenIdxSuff = "/PtGenI";
static constexpr std::string_view PtEffSuff = "/PtEfficiency";
static constexpr std::string_view PtEffIdxSuff = "/PtEfficiencyI";

// histogram registry labels
static constexpr std::string_view BinnedPrefix = "Centrality";
static constexpr std::string_view InclusivePrefix = "Inclusive";

namespace histograms
{
// events and collisions
static constexpr std::string_view BCSelection = "Events/BCSelection";             // BC selection categories
static constexpr std::string_view EventSelection = "Events/Selection";            // Collision selection categories
static constexpr std::string_view NtrkZvtx = "Events/NtrkZvtx";                   // N tracks vs vtx Z for selected collisions
static constexpr std::string_view NtrkZvtxGen = "Events/NtrkZvtxGen";             // -- for selected simulated collisions
static constexpr std::string_view NtrkZvtxGen_t = "Events/NtrkZvtxGen_t";         // N particles vs vtx Z for generated events
static constexpr std::string_view Efficiency = "Events/Efficiency";               // simulated event selection efficiency
static constexpr std::string_view EfficiencyMult = "Events/EfficiencyMult";       // simulated event selection efficiency vs generated multiplicity
static constexpr std::string_view NotFoundZvtx = "Events/NotFoundEventZvtx";      // vtx Z distribution of events without reconstructed collisions
static constexpr std::string_view Response = "Events/Response";                   // simulated multiplicity response (N tracks vs N particles)
static constexpr std::string_view MultiResponse = "Events/MultiResponse";         // -- multi-estimator
static constexpr std::string_view SplitMult = "Events/SplitMult";                 // split reconstructed events vs generated multiplicity
static constexpr std::string_view EventChi2 = "Events/Control/Chi2";              // collisions chi2 distribution
static constexpr std::string_view EventTimeRes = "Events/Control/TimeResolution"; // collisions time resolution distribution

// particles and tracks
static constexpr std::string_view EtaZvtx = "Tracks/EtaZvtx";                                 // eta vs vtx Z distribution of tracks
static constexpr std::string_view EtaZvtx_gt0 = "Tracks/EtaZvtx_gt0";                         // -- for INEL>0 collisions
static constexpr std::string_view EtaZvtx_PVgt0 = "Tracks/EtaZvtx_PVgt0";                     // -- for INEL>0 (PV)
static constexpr std::string_view EtaZvtxGen = "Tracks/EtaZvtxGen";                           // eta vs vtx Z distribution of simulated tracks
static constexpr std::string_view EtaZvtxGen_gt0 = "Tracks/EtaZvtxGen_gt0";                   // -- for INEL>0 collisions
static constexpr std::string_view EtaZvtxGen_PVgt0 = "Tracks/EtaZvtxGen_PVgt0";               // -- for INEL>0 (PV)
static constexpr std::string_view EtaZvtxGen_t = "Tracks/EtaZvtxGen_t";                       // -- of particles
static constexpr std::string_view EtaZvtxGen_gt0t = "Tracks/EtaZvtxGen_gt0t";                 // -- of particles for INEL>0 events
static constexpr std::string_view ReassignedEtaZvtx = "Tracks/Control/ReassignedEtaZvtx";     // -- of reassigned ambiguous tracks
static constexpr std::string_view ExtraEtaZvtx = "Tracks/Control/ExtraEtaZvtx";               // -- of adopted orphan tracks
static constexpr std::string_view PhiEta = "Tracks/PhiEta";                                   // eta vs phi distribution of tracks
static constexpr std::string_view PhiEtaDuplicates = "Tracks/Control/PhiEtaDuplicates";       // -- of tracks belonging to the same particle
static constexpr std::string_view PhiEtaGen = "Tracks/PhiEtaGen";                             // -- of simulated tracks
static constexpr std::string_view PhiEtaGenDuplicates = "Tracks/Control/PhiEtaGenDuplicates"; // -- of multi-reconstructed particles
static constexpr std::string_view ReassignedPhiEta = "Tracks/Control/ReassignedPhiEta";       // -- of reassigned ambiguous tracks
static constexpr std::string_view ExtraPhiEta = "Tracks/Control/ExtraPhiEta";                 // -- of adopted orphaned tracks
static constexpr std::string_view PtEta = "Tracks/Control/PtEta";                             // Pt vs eta distribution of tracks
static constexpr std::string_view PtEtaGen = "Tracks/Control/PtEtaGen";                       // Pt vs eta distribution of simulated tracks
static constexpr std::string_view DCAXYPt = "Tracks/Control/DCAXYPt";                         // transversal DCA vs Pt distribution of tracks
static constexpr std::string_view ReassignedDCAXYPt = "Tracks/Control/ReassignedDCAXYPt";     // -- of reassigned ambiguous tracks
static constexpr std::string_view ExtraDCAXYPt = "Tracks/Control/ExtraDCAXYPt";               // -- of adopted orphan tracks
static constexpr std::string_view DCAZPt = "Tracks/Control/DCAZPt";                           // longitudal DCA vs Pt distribution of tracks
static constexpr std::string_view ReassignedDCAZPt = "Tracks/Control/ReassignedDCAZPt";       // -- of reassigned ambiguous tracks
static constexpr std::string_view ExtraDCAZPt = "Tracks/Control/ExtraDCAZPt";                 // -- of adopted orphan tracks
static constexpr std::string_view ReassignedZvtxCorr = "Tracks/Control/ReassignedZvtxCorr";   // original vs reassigned vtx Z correlation for reassigned ambiguous tracks

// efficiencies
static constexpr std::string_view PtGen = "Tracks/Control/PtGen";                                               // pt distribution of particles
static constexpr std::string_view PtGenF = "Tracks/Control/{}/PtGen";                                           // -- format placeholder
static constexpr std::string_view PtGenIdx = "Tracks/Control/PtGenI";                                           // -- for the indexed efficiency
static constexpr std::string_view PtGenIdxF = "Tracks/Control/{}/PtGenI";                                       // -- format placeholder
static constexpr std::string_view PtGenNoEtaCut = "Tracks/Control/PtGenNoEtaCut";                               // -- with no eta restriction
static constexpr std::string_view PtGenIdxNoEtaCut = "Tracks/Control/PtGenINoEtaCut";                           // -- for the indexed eff. with no eta restriction
static constexpr std::string_view PtEfficiency = "Tracks/Control/PtEfficiency";                                 // generator-level pt distribution of selected tracks
static constexpr std::string_view PtEfficiencyF = "Tracks/Control/{}/PtEfficiency";                             // -- format placeholder
static constexpr std::string_view PtEfficiencyIdx = "Tracks/Control/PtEfficiencyI";                             // -- for the indexed efficiency
static constexpr std::string_view PtEfficiencyIdxF = "Tracks/Control/{}/PtEfficiencyI";                         // -- format placeholder
static constexpr std::string_view PtEfficiencyNoEtaCut = "Tracks/Control/PtEfficiencyNoEtaCut";                 // -- with no eta restriction
static constexpr std::string_view PtEfficiencyIdxNoEtaCut = "Tracks/Control/PtEfficiencyINoEtaCut";             // -- for the indexed eff. with no eta restriction
static constexpr std::string_view PtEfficiencyFakes = "Tracks/Control/PtFakes";                                 // pt distribution of fake tracks
static constexpr std::string_view PtEfficiencySecondariesIdx = "Tracks/Control/PtSecondariesI";                 // generator-level pt distribution of secondary particles
static constexpr std::string_view PtEfficiencySecondariesIdxNoEtaCut = "Tracks/Control/PtSecondariesINoEtaCut"; // -- for the indexed efficiency

// misc.
static constexpr std::string_view Mask = "Tracks/Control/Mask";           // reco status bitmask
static constexpr std::string_view ITSlayers = "Tracks/Control/ITSLayers"; // ITS layers hit distribution
} // namespace histograms
} // namespace pwgmm::mult

#endif // PWGMM_MULT_CORE_INCLUDE_HISTOGRAMS_H_
