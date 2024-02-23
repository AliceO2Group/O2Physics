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

o2physics_add_dpl_workflow(twopartcorr-register-skim
                           SOURCES twoParticleCorrelationsRegisterSkimming.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(twopartcorr-skim
                           SOURCES twoParticleCorrelationsFullSkimming.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(twopartcorr-notstored-skim
                           SOURCES twoParticleCorrelationsNotStoredSkimming.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(twopartcorr-filter
                           SOURCES twoParticleCorrelationsFiltering.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(identifiedbf-filter
                           SOURCES identifiedBfFilter.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)
