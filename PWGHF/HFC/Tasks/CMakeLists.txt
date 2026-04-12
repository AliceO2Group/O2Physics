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

o2physics_add_dpl_workflow(task-charm-hadrons-track-femto-dream
                    SOURCES taskCharmHadronsTrackFemtoDream.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-charm-hadrons-v0-femto-dream
                    SOURCES taskCharmHadronsV0FemtoDream.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-correlation-d0-hadrons
                    SOURCES taskCorrelationD0Hadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-correlation-d-dbar
                    SOURCES taskCorrelationDDbar.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-correlation-d-meson-pairs
                    SOURCES taskCorrelationDMesonPairs.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-correlation-dplus-dplus-reduced
                    SOURCES taskCorrelationDplusDplusReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-correlation-dplus-hadrons
                    SOURCES taskCorrelationDplusHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-correlation-ds-hadrons
                    SOURCES taskCorrelationDsHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-correlation-dstar-hadrons
                    SOURCES taskCorrelationDstarHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-correlation-hfe-hadrons
                    SOURCES taskCorrelationHfeHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-correlation-lc-hadrons
                    SOURCES taskCorrelationLcHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-flow
                    SOURCES taskFlow.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/PWGCF/Core)
