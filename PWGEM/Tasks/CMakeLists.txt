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

o2physics_add_dpl_workflow(phos-cell-q-a
                    SOURCES phosCellQA.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::PHOSBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phos-clu-q-a
                    SOURCES phosCluQA.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::PHOSBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phos-trig-q-a
                    SOURCES phosTrigQA.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::PHOSBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phos-pi0
                    SOURCES phosPi0.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::PHOSBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phos-calibration
                    SOURCES phosCalibration.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::PHOSBase O2::PHOSReconstruction
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phos-align
                    SOURCES phosAlign.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::PHOSBase O2::PHOSReconstruction
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phos-nbar
                    SOURCES phosNbar.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::PHOSBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phos-nonlin
                    SOURCES phosNonlin.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::PHOSBase O2::PHOSReconstruction
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phos-el-id
                    SOURCES phosElId.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2::PHOSBase
                    COMPONENT_NAME Analysis)
