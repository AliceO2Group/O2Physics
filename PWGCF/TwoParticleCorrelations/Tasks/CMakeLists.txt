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

o2physics_add_dpl_workflow(two-particle-correlations
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

o2physics_add_dpl_workflow(lambda-spin-polarization
                           SOURCES lambdaSpinPolarization.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dpt-dpt-efficiency-and-qc
                           SOURCES dptDptEfficiencyAndQc.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dpt-dpt-per-run-qc
                           SOURCES dptDptPerRunQc.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore O2Physics::AnalysisCCDB O2::DataFormatsITSMFT
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dpt-dpt-per-run-extra-qc
                           SOURCES dptDptPerRunExtraQc.cxx
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

o2physics_add_dpl_workflow(di-hadron-cor
                           SOURCES diHadronCor.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-di-hadron
                           SOURCES pidDiHadron.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(longrange-correlation
                           SOURCES longrangeCorrelation.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(longrangecorr-derived
                           SOURCES longrangecorrDerived.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(long-range-dihadron-cor
                           SOURCES longRangeDihadronCor.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(two-particle-correlation-pp
                           SOURCES twoParticleCorrelationPp.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-decorrelation
                           SOURCES flowDecorrelation.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdacascadecorrelation
    SOURCES Lambdacascadecorrelation.cxx
    PUBLIC_LINK_LIBRARIES 
        O2Physics::AnalysisCore 
        O2::Framework
    COMPONENT_NAME Analysis)