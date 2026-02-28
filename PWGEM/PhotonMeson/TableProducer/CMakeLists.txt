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

add_subdirectory(Converters)

o2physics_add_dpl_workflow(skimmer-gamma-conversion
                    SOURCES skimmerGammaConversion.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DCAFitter O2Physics::AnalysisCore KFParticle::KFParticle
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(skimmer-gamma-conversion-truthonlymc
                    SOURCES skimmerGammaConversionTruthOnlyMc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(photon-conversion-builder
                    SOURCES photonconversionbuilder.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DCAFitter O2Physics::AnalysisCore O2Physics::MLCore KFParticle::KFParticle O2Physics::TPCDriftManager
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(create-pcm
                    SOURCES createPCM.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(create-emevent-photon
                    SOURCES createEMEventPhoton.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGJECore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(associate-mc-info-photon
                    SOURCES associateMCinfoPhoton.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(skimmer-gamma-calo
                    SOURCES skimmerGammaCalo.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(bc-wise-cluster-skimmer
                    SOURCES bcWiseClusterSkimmer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2Physics::AnalysisCCDB O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(skimmer-phos
                    SOURCES skimmerPHOS.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::PHOSBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(skimmer-primary-electron-from-dalitzee
                    SOURCES skimmerPrimaryElectronFromDalitzEE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(skimmer-dalitz-ee
                    SOURCES skimmerDalitzEE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(non-lin-producer
                    SOURCES nonLinProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(material-budget-weights
                    SOURCES materialBudgetWeights.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::CCDB ROOT::Hist ROOT::Core
                    COMPONENT_NAME Analysis)
