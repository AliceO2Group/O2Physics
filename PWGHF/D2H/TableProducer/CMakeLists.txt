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

o2physics_add_dpl_workflow(candidate-creator-b-to-jpsi-reduced
                    SOURCES candidateCreatorBToJpsiReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-b0-reduced
                    SOURCES candidateCreatorB0Reduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-bplus-reduced
                    SOURCES candidateCreatorBplusReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-bs-reduced
                    SOURCES candidateCreatorBsReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-charm-reso-reduced
                    SOURCES candidateCreatorCharmResoReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-creator-lb-reduced
                    SOURCES candidateCreatorLbReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

# Candidate selectors

o2physics_add_dpl_workflow(candidate-selector-b0-to-d-pi-reduced
                    SOURCES candidateSelectorB0ToDPiReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-bplus-to-d0-pi-reduced
                    SOURCES candidateSelectorBplusToD0PiReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-bs-to-ds-pi-reduced
                    SOURCES candidateSelectorBsToDsPiReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(candidate-selector-lb-to-lc-pi-reduced
                    SOURCES candidateSelectorLbToLcPiReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

# Data creators

o2physics_add_dpl_workflow(data-creator-charm-had-pi-reduced
                    SOURCES dataCreatorCharmHadPiReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(data-creator-charm-reso-to-dstar-reduced
                    SOURCES dataCreatorCharmResoToDstarReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils KFParticle::KFParticle
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(data-creator-charm-reso-to-dplus-reduced
                    SOURCES dataCreatorCharmResoToDplusReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils KFParticle::KFParticle
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(data-creator-charm-reso-to-d0-reduced
                    SOURCES dataCreatorCharmResoToD0Reduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils KFParticle::KFParticle
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(data-creator-jpsi-had-reduced
                    SOURCES dataCreatorJpsiHadReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

# Converters

o2physics_add_dpl_workflow(converter-reduced-3-prongs-ml
                    SOURCES converterReduced3ProngsMl.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

 o2physics_add_dpl_workflow(converter-reduced-hadron-daus-pid
                     SOURCES converterReducedHadronDausPid.cxx
                     PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                     COMPONENT_NAME Analysis)
