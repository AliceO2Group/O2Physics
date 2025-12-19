# Copyright 2019-2025 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

o2physics_add_dpl_workflow(femto-track-qa
	  SOURCES femtoTrackQa.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-twotrackresonance-qa
	  SOURCES femtoTwotrackresonanceQa.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-v0-qa
	  SOURCES femtoV0Qa.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-kink-qa
	  SOURCES femtoKinkQa.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-cascade-qa
	  SOURCES femtoCascadeQa.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-pair-track-track
	  SOURCES femtoPairTrackTrack.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-pair-track-v0
	  SOURCES femtoPairTrackV0.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-pair-track-two-track-resonance
	  SOURCES femtoPairTrackTwoTrackResonance.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-pair-track-cascade
	  SOURCES femtoPairTrackCascade.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-pair-track-kink
	  SOURCES femtoPairTrackKink.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-pair-v0-v0
	  SOURCES femtoPairV0V0.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)
