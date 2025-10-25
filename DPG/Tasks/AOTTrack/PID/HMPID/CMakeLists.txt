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

# HMPID
o2physics_add_dpl_workflow(hmpid-qa
    SOURCES hmpidQa.cxx
    PUBLIC_LINK_LIBRARIES O2::HMPIDBase O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hmpid-table-producer
    SOURCES hmpidTableProducer.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hmpid-deuteron
    SOURCES hmpidDeuteron.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::ReconstructionDataFormats O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)
