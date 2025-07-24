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

# General purpose
o2physics_add_dpl_workflow(epvector
    SOURCES epvector.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(spvector
    SOURCES spvector.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tpcpid
    SOURCES lfTPCPID.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(zdcsp
    SOURCES zdcSP.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-centrality
    SOURCES mcCentrality.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(kink-builder
    SOURCES kinkBuilder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
    COMPONENT_NAME Analysis)
