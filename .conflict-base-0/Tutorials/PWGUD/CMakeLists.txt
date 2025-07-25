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

o2physics_add_dpl_workflow(udtutorial-01
                           SOURCES UDTutorial_01.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udtutorial-02a
                           SOURCES UDTutorial_02a.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udtutorial-02b
                           SOURCES UDTutorial_02b.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udtutorial-03a
                           SOURCES UDTutorial_03a.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udtutorial-03b
                           SOURCES UDTutorial_03b.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udtutorial-04
                           SOURCES UDTutorial_04.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udtutorial-05
                           SOURCES UDTutorial_05.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udtutorial-06
                           SOURCES UDTutorial_06.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(udtutorial-07
                           SOURCES UDTutorial_07.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)
