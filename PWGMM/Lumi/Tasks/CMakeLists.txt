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

o2physics_add_dpl_workflow(lumi
                        SOURCES lumi.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                        O2::ReconstructionDataFormats
                        O2::DetectorsBase
                        O2::DetectorsCommonDataFormats
                        O2::DetectorsVertexing
                        COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(lumifddft0
                        SOURCES LumiFDDFT0.cxx
                        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                        O2::ReconstructionDataFormats
                        O2::DetectorsBase
                        O2::DetectorsCommonDataFormats
                        O2::DetectorsVertexing
                        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fitvdm
        SOURCES fitLumi.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
        O2::ReconstructionDataFormats
        O2::DetectorsBase
        O2::DetectorsCommonDataFormats
        O2::DetectorsVertexing
        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lumistab
        SOURCES lumiStability.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
        O2::ReconstructionDataFormats
        COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lumi-stability-light-ions
        SOURCES lumiStabilityLightIons.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCCDB O2Physics::AnalysisCore
        COMPONENT_NAME Analysis)