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
Configurable<string> cfTaskName{
  "cfTaskName", "Default task name",
  "set task name - use eventually to determine weights for this task"};
Configurable<string> cfRunNumber{
  "cfRunNumber", "Some run number",
  "set run number - TBI temporarily here, this shall be eventually "
  "automatically obtained from e.g. collision.bc().runNumber()"};
Configurable<bool> cfVerbose{
  "cfVerbose", false,
  "run or not in verbose mode (but not for function calls per particle)"};
Configurable<bool> cfVerboseForEachParticle{
  "cfVerboseForEachParticle", false,
  "run or not in verbose mode (also for function calls per particle)"};
Configurable<bool> cfUseCCDB{
  "cfUseCCDB", true,
  "access personal files from CCDB or from home dir in AliEn"};

// Test0:
Configurable<bool> cfCalculateTest0{"cfCalculateTest0", false,
                                    "calculate or not Test0"};
Configurable<string> cfFileWithLabels{"cfFileWithLabels",
                                      "/home/abilandz/DatasetsO2/labels.root", "path to external ROOT file which specifies all labels"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"

// Particle weights:
Configurable<bool> cfUsePhiWeights{"cfUsePhiWeights", false,
                                   "use or not phi weights"};
Configurable<bool> cfUsePtWeights{"cfUsePtWeights", false,
                                  "use or not pt weights"};
Configurable<bool> cfUseEtaWeights{"cfUseEtaWeights", false,
                                   "use or not eta weights"};
Configurable<string> cfFileWithWeights{"cfFileWithWeights",
                                       "/home/abilandz/DatasetsO2/weights.root", "path to external ROOT file which holds all particle weights in O2 format"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"

// Event cuts:
Configurable<int> cNumberOfEvents_min{
  "cNumberOfEvents_min", -1,
  "minimum number of events to process (set to -1 to ignore)"};
Configurable<int> cNumberOfEvents_max{"cNumberOfEvents_max", 1000000000,
                                      "maximum number of events to process"};
Configurable<int> cTotalMultiplicity_min{
  "cTotalMultiplicity_min", -1,
  "minimum total multiplicity of an event to be processed (set to -1 to "
  "ignore)"};
Configurable<int> cTotalMultiplicity_max{
  "cTotalMultiplicity_max", 1000000000,
  "maximum total multiplicity of an event to be processed"};
Configurable<int> cSelectedTracks_min{
  "cSelectedTracks_min", -1,
  "minimum number of selected tracks (set to -1 to ignore)"};
Configurable<int> cSelectedTracks_max{"cSelectedTracks_max", 1000000000,
                                      "maximum number of selected tracks"};
Configurable<float> cCentrality_min{"cCentrality_min", 0.,
                                    "minimum centrality"};
Configurable<float> cCentrality_max{"cCentrality_max", 100.,
                                    "maximum centrality"};
Configurable<float> cVertex_x_min{"cVertex_x_min", -10.0,
                                  "minimum vertex x range [cm]"};
Configurable<float> cVertex_x_max{"cVertex_x_max", 10.0,
                                  "maximum vertex x range [cm]"};
Configurable<float> cVertex_y_min{"cVertex_y_min", -10.0,
                                  "minimum vertex y range [cm]"};
Configurable<float> cVertex_y_max{"cVertex_y_max", 10.0,
                                  "maximum vertex y range [cm]"};
Configurable<float> cVertex_z_min{"cVertex_z_min", -10.0,
                                  "minimum vertex z range [cm]"};
Configurable<float> cVertex_z_max{"cVertex_z_max", 10.0,
                                  "maximum vertex z range [cm]"};
Configurable<int> cNContributors_min{
  "cNContributors_min", -1,
  "minimum number of vertex contributors (set to -1 to ignore)"};
Configurable<int> cNContributors_max{"cNContributors_max", 1000000000,
                                     "maximum number of vertex contributors"};

// Particle cuts:
Configurable<float> pt_min{"pt_min", 0.2, "minimum track pt value [GeV/c]"};
Configurable<float> pt_max{"pt_max", 5.0, "maximum track pt value [GeV/c]"};

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
