# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

add_subdirectory(Skimming)
add_subdirectory(OpenData)
add_subdirectory(ML)
add_subdirectory(PWGEM)
add_subdirectory(PWGUD)
add_subdirectory(PWGLF)
add_subdirectory(PWGCF)
add_subdirectory(PWGMM)
add_subdirectory(PWGHF)

o2physics_add_dpl_workflow(histogram-track-selection
                    SOURCES src/histogramTrackSelection.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(track-iteration
                  SOURCES src/trackIteration.cxx
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(pid
                  SOURCES src/pidTpcTof.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(full-track-iteration
                  SOURCES src/fullTrackIteration.cxx
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(collision-tracks-iteration
                  SOURCES src/collisionTracksIteration.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(new-collections
                  SOURCES src/newCollections.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(associated-derived
                  SOURCES src/associatedExample.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(dynamic-columns
                  SOURCES src/dynamicColumns.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(histograms
                  SOURCES src/histograms.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(histograms-full-tracks
                  SOURCES src/histogramsFullTracks.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(filters
                  SOURCES src/filters.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(partitions
                  SOURCES src/partitions.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(tableio-write
                  SOURCES src/tableIOout.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(tableio-read
                  SOURCES src/tableIOin.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(jet-analysis
                  SOURCES src/jetAnalysis.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(tracks-combinations
                  SOURCES src/tracksCombinations.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(event-mixing
                  SOURCES src/eventMixing.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(event-mixing-validation
                  SOURCES src/eventMixingValidation.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(mc-histograms
                  SOURCES src/mcHistograms.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(mc-only
                  SOURCES src/mcOnly.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(ccdbaccess
                  SOURCES src/ccdbaccess.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::CCDB O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(ccdbtableaccess
                  SOURCES src/ccdbtableaccess.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::CCDB O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(weak-decay-iteration
                  SOURCES src/weakDecayIteration.cxx
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(extended-columns
                  SOURCES src/extendedColumns.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(extended-tables
                  SOURCES src/extendedTables.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(zdc-vzero-iteration
                  SOURCES src/ZDCVZeroIteration.cxx
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(histogram-registry
                  SOURCES src/histogramRegistry.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_library(ConfigurableCut
               SOURCES  src/configurableCut.cxx
               PUBLIC_INCLUDE_DIRECTORIES include
               PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore)

o2physics_target_root_dictionary(ConfigurableCut
                          HEADERS include/configurableCut.h
                          LINKDEF src/ConfigurableCutLinkDef.h)

o2physics_add_dpl_workflow(configurable-objects
                    SOURCES src/configurableObjects.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::ConfigurableCut
                    COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(muon-iteration
                  SOURCES src/muonIteration.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(schema-evolution
                  SOURCES src/schemaEvolution.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(efficiency-per-run
                  SOURCES src/efficiencyPerRun.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(efficiency-global
                  SOURCES src/efficiencyGlobal.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(compatible-bcs
                  SOURCES src/compatibleBCs.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(multiplicity-event-track-selection
                  SOURCES src/multiplicityEventTrackSelection.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(multiprocess-example
                  SOURCES src/multiProcess.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(conditional-expressions
                  SOURCES src/conditionalExpressions.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(using-pdg
                  SOURCES src/usingPDGService.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(propagated-tracks
                  SOURCES src/propagatedTracks.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(tree-creator
                  SOURCES src/treeCreator.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(reweighting
                  SOURCES src/reweighting.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ONNXRuntime::ONNXRuntime
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(hep-mc
                  SOURCES src/hepMC.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(preprocess
                  SOURCES src/preprocess.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(consume
                  SOURCES src/consume.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)
