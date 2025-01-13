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

o2physics_add_dpl_workflow(twopartcorr
                           SOURCES twoParticleCorrelations.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(r2p24id
                           SOURCES r2p2-4-id.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdacorr
                           SOURCES lambdaR2Correlation.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(identifiedbf
                           SOURCES identifiedbf.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(identifiedbf-filter-qa
                           SOURCES identifiedbfqa.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dpt-dpt-efficiency-and-qc
                           SOURCES dptDptEfficiencyAndQc.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(twopartcorr-per-run-qc
                           SOURCES perRunQc.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore O2Physics::AnalysisCCDB
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(twopartcorr-per-run-extraqc
                           SOURCES perRunExtraQc.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore O2Physics::AnalysisCCDB
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(corr-sparse
                           SOURCES corrSparse.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(neutron-proton-corr-zdc
                           SOURCES neutronProtonCorrZdc.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)




