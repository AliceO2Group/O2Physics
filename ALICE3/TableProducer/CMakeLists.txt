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

add_subdirectory(OTF)

o2physics_add_dpl_workflow(alice3-trackselection
                    SOURCES alice3-trackselection.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-trackextension
                    SOURCES alice3-trackextension.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-centrality
                    SOURCES alice3-centrality.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-decaypreselector
                    SOURCES alice3-decaypreselector.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-decayfinder
                    SOURCES alice3-decayfinder.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-multicharm-table
                    SOURCES alice3-multicharmTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-correlatorddbar
                    SOURCES alice3-correlatorDDbar.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-hf-selector-3prong
                    SOURCES alice3HfSelector3Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-hf-tree-creator-3prong
                    SOURCES alice3HfTreeCreator3Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-tracking-translator
                    SOURCES alice3TrackingTranslator.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-dq-table-maker
                    SOURCES alice3-dq-table-maker.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2Physics::AnalysisCCDB O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)
