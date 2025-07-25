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

function(o2physics_define_options)

  option(BUILD_SHARED_LIBS "Build shared libs" ON)

  option(ENABLE_CASSERT "Enable asserts" OFF)

  option(ENABLE_UPGRADES "Enable detectors for upgrades" OFF)

  option(O2PHYSICS_WARNINGS_AS_ERRORS "Treat selected warnings as errors" OFF)
endfunction()
