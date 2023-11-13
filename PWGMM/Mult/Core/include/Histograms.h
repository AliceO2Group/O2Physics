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

#ifndef PWGMM_MULT_HISTOGRAMS_H
#define PWGMM_MULT_HISTOGRAMS_H
#include <string_view>

namespace pwgmm::mult {
static constexpr std::array<std::array<std::string_view, 2>, 2> categories{
  {
   {"Tracks", "Events"},                       //
   {"Tracks/Centrality", "Events/Centrality"}  //
  }                                            //
};

namespace histograms {
// events and collisions
static constexpr std::string_view NtrkZvtx = "NtrkZvtx";                              // N tracks vs vtx Z for selected collisions
static constexpr std::string_view NtrkZvtxGen = "NtrkZvtxGen";                        // -- for selected simulated collisions
static constexpr std::string_view NtrkZvtxGen_t = "NtrkZvtxGen_t";                    // N particles vs vtx Z for generated events
static constexpr std::string_view Efficiency = "Efficiency";                          // simulated event selection efficiency
static constexpr std::string_view EfficiencyMult = "EfficiencyMult";                  // simulated event selection efficiency vs generated multiplicity
static constexpr std::string_view NotFoundZvtx = "NotFoundEventZvtx";                 // vtx Z distribution of events without reconstructed collisions
static constexpr std::string_view Response = "Response";                              // simulated multiplicity response (N tracks vs N particles)
static constexpr std::string_view MultiResponse = "MultiResponse";                    // -- multi-estimator
static constexpr std::string_view SplitMult = "SplitMult";                            // split reconstructed events vs generated multiplicity

// particles and tracks
static constexpr std::string_view EtaZvtx = "EtaZvtx";                                 // eta vs vtx Z distribution of tracks
static constexpr std::string_view EtaZvtx_gt0 = "EtaZvtx_gt0";                         // -- for INEL>0 collisions
static constexpr std::string_view EtaZvtx_PVgt0 = "EtaZvtx_PVgt0";                     // -- for INEL>0 (PV)
static constexpr std::string_view EtaZvtxGen = "EtaZvtxGen";                           // eta vs vtx Z distribution of simulated tracks
static constexpr std::string_view EtaZvtxGen_gt0 = "EtaZvtxGen_gt0";                   // -- for INEL>0 collisions
static constexpr std::string_view EtaZvtxGen_PVgt0 = "EtaZvtxGen_PVgt0";               // -- for INEL>0 (PV)
static constexpr std::string_view EtaZvtxGen_gt0t = "EtaZvtxGen_gt0t";                 // -- of particles for INEL>0 events
static constexpr std::string_view ReassignedEtaZvtx = "Control/ReassignedEtaZvtx";     // -- of reassigned ambiguous tracks
static constexpr std::string_view PhiEta = "PhiEta";                                   // eta vs phi distribution of tracks
static constexpr std::string_view PhiEtaGen = "PhiEtaGen";                             // eta vs phi distribution of simulated tracks
static constexpr std::string_view ReassignedPhiEta = "Control/ReassignedPhiEta";       // -- of reassigned ambiguous tracks
static constexpr std::string_view PtEta = "Control/PtEta";                             // Pt vs eta distribution of tracks
static constexpr std::string_view PtEtaGen = "Control/PtEtaGen";                       // Pt vs eta distribution of simulated tracks
static constexpr std::string_view DCAXYPt = "Control/DCAXYPt";                         // transversal DCA vs Pt distribution of tracks
static constexpr std::string_view ReassignedDCAXYPt = "Control/ReassignedDCAXYPt";     // -- of reassigned ambiguous tracks
static constexpr std::string_view DCAZPt = "Control/DCAZPt";                           // longitudal DCA vs Pt distribution of tracks
static constexpr std::string_view ReassignedDCAZPt = "Control/ReassignedDCAZPt";       // -- of reassigned ambiguous tracks
static constexpr std::string_view ReassignedZvtxCorr = "Control/ReassignedZvtxCorr";   // original vs reassigned vtx Z correlation for reassigned ambiguous tracks

}
}

#endif // PWGMM_MULT_HISTOGRAMS_H
