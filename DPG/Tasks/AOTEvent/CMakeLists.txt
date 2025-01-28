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

o2physics_add_dpl_workflow(event-selection-qa
                           SOURCES eventSelectionQa.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lumi-qa
                           SOURCES lumiQa.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(time-dependent-qa
                           SOURCES timeDependentQa.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsBase O2::TPCCalibration
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(detector-occupancy-qa
                           SOURCES detectorOccupancyQa.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(rof-border-qa
                           SOURCES rofBorderQa.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(matching-qa
                           SOURCES matchingQa.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(rof-occupancy-qa
                           SOURCES rofOccupancyQa.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::DetectorsBase
                           COMPONENT_NAME Analysis)
