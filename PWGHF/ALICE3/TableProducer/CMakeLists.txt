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

# Candidate creators

o2physics_add_dpl_workflow(candidate-creator-chic
                    SOURCES candidateCreatorChic.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-x
                    SOURCES candidateCreatorX.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-xicc
                    SOURCES candidateCreatorXicc.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

# Candidate selectors

o2physics_add_dpl_workflow(candidate-selector-chic-to-jpsi-gamma
                    SOURCES candidateSelectorChicToJpsiGamma.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-d0-alice3-barrel
                    SOURCES candidateSelectorD0Alice3Barrel.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-d0-alice3-forward
                    SOURCES candidateSelectorD0Alice3Forward.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-d0-parametrized-pid
                    SOURCES candidateSelectorD0ParametrizedPid.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-jpsi
                    SOURCES candidateSelectorJpsi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-alice3
                    SOURCES candidateSelectorLcAlice3.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lc-parametrized-pid
                    SOURCES candidateSelectorLcParametrizedPid.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-x-to-jpsi-pi-pi
                    SOURCES candidateSelectorXToJpsiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-xicc-to-p-k-pi-pi
                    SOURCES candidateSelectorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Tree creators

o2physics_add_dpl_workflow(tree-creator-chic-to-jpsi-gamma
                    SOURCES treeCreatorChicToJpsiGamma.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-x-to-jpsi-pi-pi
                    SOURCES treeCreatorXToJpsiPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tree-creator-xicc-to-p-k-pi-pi
                    SOURCES treeCreatorXiccToPKPiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
