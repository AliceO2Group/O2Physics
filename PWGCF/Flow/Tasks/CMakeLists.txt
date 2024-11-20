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
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-task
                    SOURCES FlowTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-runbyrun
                    SOURCES FlowRunbyRun.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-gfw-pbpb
                    SOURCES FlowGFWPbPb.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-zdc-task
                    SOURCES FlowZDCtask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-gf
                    SOURCES flowAnalysisGF.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-pbpb-pikp-task
                    SOURCES FlowPbPbpikp.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-gfw-omegaxi
                    SOURCES flowGFWOmegaXi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)
