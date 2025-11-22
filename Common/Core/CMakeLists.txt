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

o2physics_add_library(AnalysisCore
        SOURCES TrackSelection.cxx
        OrbitRange.cxx
        PID/ParamBase.cxx
        PID/PIDTOF.cxx
        PID/PIDTOFParamService.cxx
        CollisionAssociation.cxx
        TrackSelectionDefaults.cxx
        EventPlaneHelper.cxx
        TableHelper.cxx
        MetadataHelper.cxx
        CollisionTypeHelper.cxx
        FFitWeights.cxx
        PUBLIC_LINK_LIBRARIES O2::Framework O2::DataFormatsParameters ROOT::EG O2::CCDB ROOT::Physics O2::FT0Base O2::FV0Base O2::DataFormatsParamTOF)

o2physics_target_root_dictionary(AnalysisCore
        HEADERS TrackSelection.h
        TrackSelectionDefaults.h
        OrbitRange.h
        TableHelper.h
        MetadataHelper.h
        EventPlaneHelper.h
        PID/ParamBase.h
        PID/DetectorResponse.h
        PID/PIDTOF.h
        PID/TPCPIDResponse.h
        PID/PIDTOFParamService.h
        CollisionTypeHelper.h
        FFitWeights.h
        LINKDEF AnalysisCoreLinkDef.h)

o2physics_add_library(EventFilteringUtils
    SOURCES Zorro.cxx ZorroSummary.cxx
    INSTALL_HEADERS ZorroHelper.h ZorroSummary.h
    PUBLIC_LINK_LIBRARIES O2::Framework O2::CCDB ROOT::EG O2::CCDB ROOT::Physics Arrow::arrow_shared)

o2physics_target_root_dictionary(EventFilteringUtils
   HEADERS ZorroHelper.h ZorroSummary.h
   LINKDEF EventFilteringUtilsLinkDef.h)

o2physics_add_header_only_library(TPCDriftManager
        HEADERS TPCVDriftManager.h
        INTERFACE_LINK_LIBRARIES O2::DataFormatsTPC)
