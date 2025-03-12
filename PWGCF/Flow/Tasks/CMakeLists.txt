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

o2physics_add_dpl_workflow(flow-pt-efficiency
                    SOURCES flowPtEfficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-task
                    SOURCES flowTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-runby-run
                    SOURCES flowRunbyRun.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-mc
                    SOURCES flowMc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-qa
                    SOURCES flowQa.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-gfw-task
                    SOURCES flowGfwTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-zdc-task
                    SOURCES flowZdcTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-gf
                    SOURCES flowAnalysisGF.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-pbpb-pikp
                    SOURCES flowPbpbPikp.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-gfw-omega-xi
                    SOURCES flowGfwOmegaXi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-pid-cme
                    SOURCES pidcme.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-sp
                    SOURCES flowSP.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(resonances-gfw-flow
                    SOURCES resonancesGfwFlow.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-efficiency-casc
                    SOURCES flowEfficiencyCasc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)
