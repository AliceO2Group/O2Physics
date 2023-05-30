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
#ifndef PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
#define PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_

#include <string>

// Task configuration:
Configurable<string> cfTaskName{"cfTaskName", "Default task name", "set task name - use eventually to determine weights for this task"};
Configurable<bool> cfVerbose{"cfVerbose", true, "run or not in verbose mode"};

// Test0:
Configurable<bool> cfCalculateTest0{"cfCalculateTest0", false, "calculate or not Test0"};
Configurable<string> cfLabels{"cfLabels", "/home/abilandz/DatasetsO2/labels.root", "TBI description"};

// Event cuts:
Configurable<float> Vz_min{"Vz_min", -10.0, "minimum vertex z range [cm]"};
Configurable<float> Vz_max{"Vz_max", 10.0, "maximum vertex z range [cm]"};

// Particle cuts:
Configurable<float> pt_min{"pt_min", 0.2, "minimum track pt value [GeV/c]"};
Configurable<float> pt_max{"pt_max", 5.0, "maximum track pt value [GeV/c]"};

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
