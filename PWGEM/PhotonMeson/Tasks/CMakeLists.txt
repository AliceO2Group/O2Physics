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

o2physics_add_dpl_workflow(gammaconversions
                    SOURCES gammaConversions.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(gammaconversionstruthonlymc
                    SOURCES gammaConversionsTruthOnlyMc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(emc-pi0-qc
                    SOURCES emcalPi0QC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pcm-qc
                    SOURCES pcmQC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pcm-qc-mc
                    SOURCES pcmQCMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dalitz-ee-qc
                    SOURCES dalitzEEQC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dalitz-ee-qc-mc
                    SOURCES dalitzEEQCMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore O2Physics::PWGEMDileptonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(emcal-qc
                    SOURCES emcalQC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phos-qc
                    SOURCES phosQC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(single-photon
                    SOURCES SinglePhoton.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(single-photon-mc
                    SOURCES SinglePhotonMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pi0eta-to-gammagamma-pcmpcm
                    SOURCES Pi0EtaToGammaGammaPCMPCM.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pi0eta-to-gammagamma-pcmdalitzee
                    SOURCES Pi0EtaToGammaGammaPCMDalitzEE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pi0eta-to-gammagamma-emcemc
                    SOURCES Pi0EtaToGammaGammaEMCEMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pi0eta-to-gammagamma-mc-pcmpcm
                    SOURCES Pi0EtaToGammaGammaMCPCMPCM.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pi0eta-to-gammagamma-mc-pcmdalitzee
                    SOURCES Pi0EtaToGammaGammaMCPCMDalitzEE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pi0eta-to-gammagamma-mc-emcemc
                    SOURCES Pi0EtaToGammaGammaMCEMCEMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tagging-pi0
                    SOURCES TaggingPi0.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tagging-pi0-mc
                    SOURCES TaggingPi0MC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tag-and-probe
                    SOURCES TagAndProbe.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(material-budget
                    SOURCES MaterialBudget.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(material-budget-mc
                    SOURCES MaterialBudgetMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(check-mc-v0
                    SOURCES CheckMCV0.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pi0-flow-emc
                    SOURCES taskPi0FlowEMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(prefilter-photon
                    SOURCES prefilterPhoton.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)

