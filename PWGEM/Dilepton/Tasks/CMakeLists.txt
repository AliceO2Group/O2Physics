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


o2physics_add_dpl_workflow(efficiency-ee
                    SOURCES emEfficiencyEE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_executable(lmee-lf-cocktail
                  SOURCES lmeeLFCocktail.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2::SimulationDataFormat O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lmee-hf-cocktail
                  SOURCES lmeeHFCocktail.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-templates
                    SOURCES MCtemplates.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(smearing
                    SOURCES smearing.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(table-reader-barrel
                    SOURCES tableReaderBarrel.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dielectron-qc
                    SOURCES dielectronQC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2::DCAFitter O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dielectron-qc-mc
                    SOURCES dielectronQCMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2::DCAFitter O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dimuon-qc
                    SOURCES dimuonQC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2::DCAFitter O2::GlobalTracking O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dimuon-qc-mc
                    SOURCES dimuonQCMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2::DCAFitter O2::GlobalTracking O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(emu-correlation
                    SOURCES emuCorrelation.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2::GlobalTracking O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(emu-correlation-mc
                    SOURCES emuCorrelationMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::MLCore O2::GlobalTracking O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(photon-hbt-pcmpcm
                    SOURCES PhotonHBTPCMPCM.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(photon-hbt-eeee
                    SOURCES PhotonHBTEEEE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(photon-hbt-pcmee
                    SOURCES PhotonHBTPCMEE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

