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

add_subdirectory(Matching)

if(FastJet_FOUND)

o2physics_add_dpl_workflow(jet-deriveddata-producer
                    SOURCES derivedDataProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::EventFilteringUtils O2Physics::SGCutParHolder
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-deriveddata-trigger-producer
                    SOURCES derivedDataTriggerProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-deriveddata-selector
                    SOURCES derivedDataSelector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-deriveddata-writer
                    SOURCES derivedDataWriter.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-luminosity-producer
                    SOURCES luminosityProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-luminosity-calculator
                    SOURCES luminosityCalculator.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-eventweight-mcd
                    SOURCES jetEventWeightMCD.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-eventweight-mcp
                    SOURCES jetEventWeightMCP.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-outlier-rejector
                    SOURCES mcOutlierRejector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-track-derived
                    SOURCES jetTrackDerived.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-hf-definition
                    SOURCES heavyFlavourDefinition.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-taggerhf
                    SOURCES jetTaggerHF.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(estimator-rho
                    SOURCES rhoEstimator.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(subtractor-eventwiseconstituent
                    SOURCES eventwiseConstituentSubtractor.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-sv-reconstruction
                    SOURCES secondaryVertexReconstruction.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

endif()


o2physics_add_dpl_workflow(emcal-correction-task
                    SOURCES emcalCorrectionTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2::EMCALBase O2::EMCALReconstruction O2::EMCALCalibration O2Physics::PWGJECore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(emcal-matchedtracks-writer
                    SOURCES emcalMatchedTracksTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2::EMCALBase O2::EMCALReconstruction
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(emcal-cluster-hadronic-correction-task
                    SOURCES emcalClusterHadronicCorrectionTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2::EMCALBase O2::EMCALReconstruction
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(slim-tables-producer
                    SOURCES slimTablesProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
