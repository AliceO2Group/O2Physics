# Copyright 2019-2099 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

if (TARGET Rivet)
  o2physics_add_dpl_workflow(rivet
    SOURCES rivet.cxx
    PUBLIC_LINK_LIBRARIES
    O2::Framework
    O2::FrameworkFoundation
    O2::SimulationDataFormat
    O2::Headers
    O2::Steer
    # O2::Mergers
    O2::Generators
    Rivet::Rivet
    Yoda::Yoda
    O2Physics::RivetOutput
    COMPONENT_NAME Analysis)
endif(TARGET Rivet)
#
# EOF
#
