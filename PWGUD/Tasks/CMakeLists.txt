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

o2physics_add_dpl_workflow(upc
                           SOURCES upcAnalysis.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-spectra
                           SOURCES SGSpectraAnalyzer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-fourpi
                           SOURCES SG_FourPi_Analyzer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-d0
                           SOURCES SG_D0_Analyzer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-fit
                           SOURCES SG_FIT_Analyzer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-forward
                           SOURCES upcForward.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(diff-mcdatascanner
                           SOURCES diffMCDataScanner.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(diff-mcqa
                           SOURCES diffMCQA.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGCutparHolder
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(diff-qa
                           SOURCES diffQA.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGCutparHolder
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dgcand-analyzer
                           SOURCES DGCandAnalyzer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGCutparHolder O2Physics::UDGoodRunSelector O2Physics::DGPIDSelector O2Physics::UDFSParser
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upccand-analyzer
                           SOURCES UPCCandidateAnalyzer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upccand-producer-qa
                           SOURCES UPCCandidateProducerQa.cpp
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-mft
                           SOURCES upcMft.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsBase O2::DetectorsCommonDataFormats
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-veto
                           SOURCES upcVetoAnalysis.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tautau13topo
                           SOURCES upcTauTau13topo.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(upc-tau-rl
                           SOURCES upcTauCentralBarrelRL.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsBase O2::DetectorsCommonDataFormats
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(polarisation-rho
                           SOURCES PolarisationRho.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)
