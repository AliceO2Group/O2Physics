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
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(gammaconversionstruthonlymc
                    SOURCES gammaConversionsTruthOnlyMc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phoscellqa
                    SOURCES phosCellQA.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::DetectorsVertexing O2::PHOSBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phoscluqa
                    SOURCES phosCluQA.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::DetectorsVertexing O2::PHOSBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(skimmergammaconversions
                    SOURCES skimmerGammaConversions.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(skimmergammaconversionstruthonlymc
                    SOURCES skimmerGammaConversionsTruthOnlyMc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)
