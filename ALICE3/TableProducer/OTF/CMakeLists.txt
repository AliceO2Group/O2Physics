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

o2physics_add_dpl_workflow(onthefly-tracker
                    SOURCES onTheFlyTracker.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsCommonDataFormats O2::DetectorsVertexing O2::DCAFitter O2Physics::ALICE3Core O2Physics::FastTracker
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(onthefly-tofpid
                    SOURCES onTheFlyTofPid.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsCommonDataFormats O2Physics::ALICE3Core
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(onthefly-richpid
                    SOURCES onTheFlyRichPid.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsCommonDataFormats O2Physics::ALICE3Core
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(on-the-fly-tracker-pid
                    SOURCES onTheFlyTrackerPid.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsCommonDataFormats O2Physics::ALICE3Core
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(on-the-fly-detector-geometry-provider
                    SOURCES onTheFlyDetectorGeometryProvider.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::ALICE3Core O2Physics::FastTracker
                    COMPONENT_NAME Analysis)
