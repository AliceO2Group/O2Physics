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

include_guard()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})

include(FeatureSummary)

find_package(O2 CONFIG)
set_package_properties(O2 PROPERTIES TYPE REQUIRED)

find_package(KFParticle)
set_package_properties(KFParticle PROPERTIES TYPE REQUIRED)

find_package(fjcontrib)
set_package_properties(fjcontrib PROPERTIES TYPE REQUIRED)

find_package(ONNXRuntime)

feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)
