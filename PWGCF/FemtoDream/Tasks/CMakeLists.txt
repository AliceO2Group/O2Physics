# Copyright 2019-2024 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

o2physics_add_dpl_workflow(femtodream-pair-track-track
          SOURCES femtoDreamPairTaskTrackTrack.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtodream-triplet-track-track-track
          SOURCES femtoDreamTripletTaskTrackTrackTrack.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtodream-pair-track-v0
          SOURCES femtoDreamPairTaskTrackV0.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtodream-pair-track-cascade
          SOURCES femtoDreamPairTaskTrackCascade.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtodream-triplet-track-track-v0
          SOURCES femtoDreamTripletTaskTrackTrackV0.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtodream-debug-track
          SOURCES femtoDreamDebugTrack.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtodream-debug-v0
          SOURCES femtoDreamDebugV0.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtodream-debug-cascade
          SOURCES femtoDreamDebugCascade.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtodream-collision-masker
          SOURCES femtoDreamCollisionMasker.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femtodream-hash
          SOURCES femtoDreamHashTask.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)
