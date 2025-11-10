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

# Core event filtering utilities library
o2physics_add_library(EventFilteringUtils
    SOURCES Zorro.cxx ZorroSummary.cxx
    INSTALL_HEADERS ../ZorroHelper.h ../ZorroSummary.h
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore Arrow::arrow_shared)

o2physics_target_root_dictionary(EventFilteringUtils
   HEADERS ../ZorroHelper.h ../ZorroSummary.h
   LINKDEF EventFilteringUtilsLinkDef.h)
