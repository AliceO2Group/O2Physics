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
                       CollisionAssociation.cxx
                       TrackSelectionDefaults.cxx
               PUBLIC_LINK_LIBRARIES O2::Framework ROOT::EG O2::CCDB ROOT::Physics)

o2physics_target_root_dictionary(AnalysisCore
              HEADERS   TrackSelection.h
                        TrackSelectionDefaults.h
                        OrbitRange.h
                        PID/ParamBase.h
                        PID/DetectorResponse.h
                        PID/PIDTOF.h
                        PID/TPCPIDResponse.h
              LINKDEF AnalysisCoreLinkDef.h)

o2physics_add_executable(merger
              COMPONENT_NAME aod
              SOURCES aodMerger.cxx
              PUBLIC_LINK_LIBRARIES ROOT::Core ROOT::Net)

o2physics_add_executable(thinner
              COMPONENT_NAME aod
              SOURCES aodThinner.cxx
              PUBLIC_LINK_LIBRARIES ROOT::Core ROOT::Net)
