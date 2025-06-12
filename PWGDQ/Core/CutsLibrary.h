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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//

#ifndef PWGDQ_CORE_CUTSLIBRARY_H_
#define PWGDQ_CORE_CUTSLIBRARY_H_

#include <string>
#include <vector>
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/VarManager.h"

// ///////////////////////////////////////////////
//   These are the Cuts used in the CEFP Task   //
//   to select tracks in the event selection    //
// ///////////////////////////////////////////////
//
// Electron 2022 cuts :
//    - Single e pT > 1.0 GeV/c
//    - Single e eta = [-0.9 ; 0.9]
//    - Is SPD Any : yes
//    - TPC chi2 < 4.0
//    - TPC N Clusters = [70 ; 160]
//    - n-sigma_e = [-4.0 ; 4.0]
//    - n-sigma_pi > 2.5
//    - n-sigma_pr > 2.5
//    - PID post-calibration : Yes
//    - Track-collision association : No
//  For the dielectron Cut : pairNoCut
//    mee > 0 GeV/c2
//
// Electron 2023 & 2024 cuts :
//    - Single e pT > 1.0 GeV/c
//    - Single e eta = [-0.9 ; 0.9]
//    - Is SPD Any : yes
//    - TPC chi2 < 4.0
//    - TPC N Clusters = [70 ; 160]
//    - n-sigma_e = [-4.0 ; 4.0]
//    - n-sigma_pi > 2.5
//    - n-sigma_pr > 2.5
//    - PID post-calibration : No
//    - Track-collision association : Yes
//  For the dielectron Cut : pairMassLow5
//    mee > 1.8 GeV/c2
//
// Low Mass electrons 2023 & 2024 cuts :
//    - Single e pT > 0.4 GeV/c
//    - Single e eta = [-0.8 ; 0.8]
//    - Is SPD Any : yes
//    - TPC chi2 < 4.0
//    - ITS chi2 < 6.0
//    - TPC N Clusters = [70 ; 170]
//    - ITS N Clusters = [3.5 ; 7.5]
//    - n-sigma_e = [-4.0 ; 4.0]
//    - n-sigma_pi > 3.5 for 0.0 < pIN < 2.0
//    - n-sigma_pr > 2.5 for 2.0 < pIN < 1e+10
//    - n-sigma_e TOF = [-4.0 ; 4.0] for 0.3 < pIN < 1e+10
//    - PID post-calibration : No
//    - Track-collision association : Yes
//  For the dielectron Cut :
//    - Intermediate Mass Range ee trigger : mee > 1.3 GeV/c2 (pairMass1_3)
//    - High Mass Range ee trigger : mee > 3.5 GeV/c2 (pairMassLow12)
//
//    Muons Cuts 2022 :
//     - Single mu Low (High) pT > 0.7 (4.0) GeV/c
//     - Single mu eta = [-4.0 ; -2.5]
//     - Rabs = [17.6 ; 89.5]
//     - p x DCA = ~6 sigma_[p x DCA]
//     - Matching MCH-MID : No
//     - Track-collision association : No
//    For the dimuon cut
//     m_mumu > 1.8 GeV/c2
//
//    Muons Cuts 2023 & 2024 :
//     - Single mu Low (High) pT > 0.7 (20.0) GeV/c
//     - Single mu eta = [-4.0 ; -2.5]
//     - Rabs = [17.6 ; 89.5]
//     - p x DCA = ~10 sigma_[p x DCA]
//     - Matching MCH-MID : Yes
//     - Track-collision association : Yes
//    For the dimuon cut
//     m_mumu > 1.8 GeV/c2
//
//
// ///////////////////////////////////////////////
//           End of Cuts for CEFP               //
// ///////////////////////////////////////////////

#include "rapidjson/document.h"

namespace o2::aod
{
namespace dqcuts
{
AnalysisCompositeCut* GetCompositeCut(const char* cutName);
AnalysisCut* GetAnalysisCut(const char* cutName);

std::vector<AnalysisCut*> GetCutsFromJSON(const char* json);
// AnalysisCut** GetCutsFromJSON(const char* json);
template <typename T>
bool ValidateJSONAnalysisCut(T cut);
template <typename T>
bool ValidateJSONAddCut(T cut, bool isSimple);
template <typename T>
AnalysisCut* ParseJSONAnalysisCut(T cut, const char* cutName);
template <typename T>
bool ValidateJSONAnalysisCompositeCut(T cut);
template <typename T>
AnalysisCompositeCut* ParseJSONAnalysisCompositeCut(T key, const char* cutName);
} // namespace dqcuts
} // namespace o2::aod

AnalysisCompositeCut* o2::aod::dqcuts::GetCompositeCut(const char* cutName);
#endif // PWGDQ_CORE_CUTSLIBRARY_H_
