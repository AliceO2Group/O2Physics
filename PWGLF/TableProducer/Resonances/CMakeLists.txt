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

# Resonances
o2physics_add_dpl_workflow(f1protoninitializer
    SOURCES f1protonInitializer.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(f1protonreducedtable
    SOURCES f1protonreducedtable.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(doublephitable
    SOURCES doublephitable.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(filterf1proton
    SOURCES filterf1proton.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(resonance-initializer
    SOURCES resonanceInitializer.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(resonance-module-initializer
    SOURCES resonanceModuleInitializer.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(resonance-merge-df
    SOURCES resonanceMergeDF.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(heptaquarktable
    SOURCES HeptaQuarktable.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cksspinalignment
    SOURCES cksspinalignment.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(resonance-tree-creator
    SOURCES resonanceTreeCreator.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)
