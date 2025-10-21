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

add_subdirectory(Duplicates)
add_subdirectory(Substructure)

if(FastJet_FOUND)

o2physics_add_dpl_workflow(jet-matching-mc-ch
                    SOURCES jetMatchingMCCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-full
                    SOURCES jetMatchingMCFull.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-neutral
                    SOURCES jetMatchingMCNeutral.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-d0-ch
                    SOURCES jetMatchingMCD0Charged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-dplus-ch
                    SOURCES jetMatchingMCDplusCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-ds-ch
                    SOURCES jetMatchingMCDsCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-dstar-ch
                    SOURCES jetMatchingMCDstarCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-lc-ch
                    SOURCES jetMatchingMCLcCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-b0-ch
                    SOURCES jetMatchingMCB0Charged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-bplus-ch
                    SOURCES jetMatchingMCBplusCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-xictoxipipi-ch
                    SOURCES jetMatchingMCXicToXiPiPiCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-dielectron-ch
                    SOURCES jetMatchingMCDielectronCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-v0-ch
                    SOURCES jetMatchingMCV0Charged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-ch
                    SOURCES jetMatchingMCSubCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-d0-ch
                    SOURCES jetMatchingMCSubD0Charged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-dplus-ch
                    SOURCES jetMatchingMCSubDplusCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-ds-ch
                    SOURCES jetMatchingMCSubDsCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-dstar-ch
                    SOURCES jetMatchingMCSubDstarCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-lc-ch
                    SOURCES jetMatchingMCSubLcCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-b0-ch
                    SOURCES jetMatchingMCSubB0Charged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-bplus-ch
                    SOURCES jetMatchingMCSubBplusCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-xictoxipipi-ch
                    SOURCES jetMatchingMCSubXicToXiPiPiCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-mc-sub-dielectron-ch
                    SOURCES jetMatchingMCSubDielectronCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-ch
                    SOURCES jetMatchingSubCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-d0-ch
                    SOURCES jetMatchingSubD0Charged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-dplus-ch
                    SOURCES jetMatchingSubDplusCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-ds-ch
                    SOURCES jetMatchingSubDsCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-dstar-ch
                    SOURCES jetMatchingSubDstarCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-lc-ch
                    SOURCES jetMatchingSubLcCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-bplus-ch
                    SOURCES jetMatchingSubBplusCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-b0-ch
                    SOURCES jetMatchingSubB0Charged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-xictoxipipi-ch
                    SOURCES jetMatchingSubXicToXiPiPiCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-matching-sub-dielectron-ch
                    SOURCES jetMatchingSubDielectronCharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

endif()
