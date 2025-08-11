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

o2physics_add_library(ALICE3Core
                      SOURCES TOFResoALICE3.cxx
                              TrackUtilities.cxx
                              DelphesO2TrackSmearer.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework
                                            O2Physics::AnalysisCore)

o2physics_target_root_dictionary(ALICE3Core
                      HEADERS TOFResoALICE3.h
                              TrackUtilities.h
                              DelphesO2TrackSmearer.h
                      LINKDEF ALICE3CoreLinkDef.h)

o2physics_add_library(FastTracker
                      SOURCES FastTracker.cxx
                              DetLayer.cxx
                              DelphesO2LutWriter.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework
                                            O2Physics::AnalysisCore
                                            O2Physics::ALICE3Core)

o2physics_target_root_dictionary(FastTracker
                      HEADERS FastTracker.h
                              DetLayer.h
                              DelphesO2LutWriter.h
                      LINKDEF FastTrackerLinkDef.h)
