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

o2physics_add_dpl_workflow(correlator-d0-d0bar
                    SOURCES correlatorD0D0bar.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-d0-d0bar-barrel-full-pid
                    SOURCES correlatorD0D0barBarrelFullPid.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-d0-hadrons
                    SOURCES correlatorD0Hadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-d-meson-pairs
                    SOURCES correlatorDMesonPairs.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-dplus-dminus
                    SOURCES correlatorDplusDminus.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-dplus-dplus-reduced
                    SOURCES correlatorDplusDplusReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-dplus-hadrons
                    SOURCES correlatorDplusHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-ds-hadrons
                    SOURCES correlatorDsHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derived-data-creator-correlations-reduced
                    SOURCES derivedDataCreatorCorrelationsReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-flow-charm-hadrons-reduced
                    SOURCES correlatorFlowCharmHadronsReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-ds-hadrons-reduced
                    SOURCES correlatorDsHadronsReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-dstar-hadrons
                    SOURCES correlatorDstarHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-hfe-hadrons
                    SOURCES correlatorHfeHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-lc-hadrons
                    SOURCES correlatorLcHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(correlator-lc-sc-hadrons
                    SOURCES correlatorLcScHadrons.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(producer-charm-hadrons-track-femto-dream
                    SOURCES producerCharmHadronsTrackFemtoDream.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(producer-charm-hadrons-v0-femto-dream
                    SOURCES producerCharmHadronsV0FemtoDream.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)
