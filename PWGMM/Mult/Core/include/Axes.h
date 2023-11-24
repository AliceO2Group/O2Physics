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

#ifndef PWGMM_MULT_CORE_INCLUDE_AXES_H_
#define PWGMM_MULT_CORE_INCLUDE_AXES_H_
#include "Framework/HistogramSpec.h"
#include "CommonConstants/MathConstants.h"

namespace pwgmm::mult
{
using namespace o2::framework;
static constexpr std::string_view ptAxisName = "p_{T} (GeV/c)";
// common axis definitions
AxisSpec ZAxis = {301, -30.1, 30.1, "Z_{vtx} (cm)"};    // Z vertex in cm
AxisSpec DeltaZAxis = {61, -6.1, 6.1, "#Delta Z (cm)"}; // Z vertex difference in cm
AxisSpec DCAAxis = {601, -3.01, 3.01};                  // DCA in cm
AxisSpec EtaAxis = {22, -2.2, 2.2, "#eta"};             // Eta

AxisSpec PhiAxis = {629, 0, o2::constants::math::TwoPI, "#phi"}; // Phi (azimuthal angle)
AxisSpec PtAxis = {2401, -0.005, 24.005, ptAxisName.data()};     // Large fine-binned Pt
// Large wide-binned Pt (for efficiency)
AxisSpec PtAxisEff = {{0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                       1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0},
                      ptAxisName.data()};
AxisSpec PtAxis_wide = {1041, -0.05, 104.05, ptAxisName.data()};         // Smaller wider-binned Pt
AxisSpec FT0CAxis = {1001, -0.5, 1000.5, "FT0C amplitude (arb. units)"}; // FT0C amplitudes
AxisSpec FT0AAxis = {3001, -0.5, 3000.5, "FT0A amplitude (arb. units)"}; // FT0A amplitudes
AxisSpec FDDAxis = {3001, -0.5, 3000.5, "FDD amplitude (arb. units)"};   // FDD amplitudes
AxisSpec RapidityAxis = {102, -10.2, 10.2, "Y"};                         // Rapidity
AxisSpec ScaleAxis = {121, -0.5, 120.5, "Event scale (GeV)"};            // Event scale
AxisSpec MPIAxis = {51, -0.5, 50.5, "N_{MPI}"};                          // N_{MPI}
AxisSpec ProcAxis = {21, 89.5, 110.5};                                   // Process flag

// event selection/efficiency binning
enum struct EvSelBins : int {
  kAll = 1,
  kSelected = 2,
  kSelectedgt0 = 3,
  kSelectedPVgt0 = 4,
  kRejected = 5
};

// labels for event selection axis
std::array<std::string_view, static_cast<size_t>(EvSelBins::kRejected) + 1> EvSelBinLabels{
  "dummy",
  "All",
  "Selected",
  "Selected INEL>0",
  "Selected INEL>0 (PV)",
  "Rejected"};

enum struct EvEffBins : int {
  kGen = 1,
  kGengt0 = 2,
  kRec = 3,
  kSelected = 4,
  kSelectedgt0 = 5,
  kSelectedPVgt0 = 6
};

// labels for event efficiency axis
std::array<std::string_view, static_cast<size_t>(EvEffBins::kSelectedPVgt0) + 1> EvEffBinLabels{
  "dummy",
  "Generated",
  "Generated INEL>0",
  "Reconstructed",
  "Selected",
  "Selected INEL>0",
  "Selected INEL>0 (PV)"};
} // namespace pwgmm::mult
#endif // PWGMM_MULT_CORE_INCLUDE_AXES_H_
