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

o2physics_add_dpl_workflow(two-particle-correlations-register-skimming
                           SOURCES twoParticleCorrelationsRegisterSkimming.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(two-particle-correlations-full-skimming
                           SOURCES twoParticleCorrelationsFullSkimming.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(two-particle-correlations-not-stored-skimming
                           SOURCES twoParticleCorrelationsNotStoredSkimming.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(two-particle-correlations-filtering
                           SOURCES twoParticleCorrelationsFiltering.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(longrange-maker
                           SOURCES longrangeMaker.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore O2Physics::SGCutParHolder
                           COMPONENT_NAME Analysis)
