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

if(FastJet_FOUND)

o2physics_add_dpl_workflow(jet-finder-data-charged
                    SOURCES jetFinderDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-mcd-charged
                    SOURCES jetFinderMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-mcp-charged
                    SOURCES jetFinderMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-data-full
                    SOURCES jetFinderDataFull.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-mcd-full
                    SOURCES jetFinderMCDFull.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-mcp-full
                    SOURCES jetFinderMCPFull.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-data-neutral
                    SOURCES jetFinderDataNeutral.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-mcd-neutral
                    SOURCES jetFinderMCDNeutral.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-mcp-neutral
                    SOURCES jetFinderMCPNeutral.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)


o2physics_add_dpl_workflow(jet-finder-d0-data-charged
                    SOURCES jetFinderD0DataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-d0-mcd-charged
                    SOURCES jetFinderD0MCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-d0-mcp-charged
                    SOURCES jetFinderD0MCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dplus-data-charged
                    SOURCES jetFinderDplusDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dplus-mcd-charged
                    SOURCES jetFinderDplusMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dplus-mcp-charged
                    SOURCES jetFinderDplusMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-ds-data-charged
                    SOURCES jetFinderDsDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-ds-mcd-charged
                    SOURCES jetFinderDsMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-ds-mcp-charged
                    SOURCES jetFinderDsMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dstar-data-charged
                    SOURCES jetFinderDstarDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dstar-mcd-charged
                    SOURCES jetFinderDstarMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dstar-mcp-charged
                    SOURCES jetFinderDstarMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-lc-data-charged
                    SOURCES jetFinderLcDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-lc-mcd-charged
                    SOURCES jetFinderLcMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-lc-mcp-charged
                    SOURCES jetFinderLcMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-b0-data-charged
                    SOURCES jetFinderB0DataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-b0-mcd-charged
                    SOURCES jetFinderB0MCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-b0-mcp-charged
                    SOURCES jetFinderB0MCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-bplus-data-charged
                    SOURCES jetFinderBplusDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-bplus-mcd-charged
                    SOURCES jetFinderBplusMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-bplus-mcp-charged
                    SOURCES jetFinderBplusMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-xictoxipipi-data-charged
                    SOURCES jetFinderXicToXiPiPiDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-xictoxipipi-mcd-charged
                    SOURCES jetFinderXicToXiPiPiMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-xictoxipipi-mcp-charged
                    SOURCES jetFinderXicToXiPiPiMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-v0-data-charged
                    SOURCES jetFinderV0DataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-v0-mcd-charged
                    SOURCES jetFinderV0MCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-v0-mcp-charged
                    SOURCES jetFinderV0MCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dielectron-data-charged
                    SOURCES jetFinderDielectronDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dielectron-mcd-charged
                    SOURCES jetFinderDielectronMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dielectron-mcp-charged
                    SOURCES jetFinderDielectronMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-d0d0bar-data-charged
                    SOURCES jetFinderD0D0BarDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-d0d0bar-mcd-charged
                    SOURCES jetFinderD0D0BarMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-d0d0bar-mcp-charged
                    SOURCES jetFinderD0D0BarMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dplusdminus-data-charged
                    SOURCES jetFinderDplusDminusDataCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dplusdminus-mcd-charged
                    SOURCES jetFinderDplusDminusMCDCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jet-finder-dplusdminus-mcp-charged
                    SOURCES jetFinderDplusDminusMCPCharged.cxx
                    REUSE_FROM JetDataModelPCH
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::PWGJECore O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

endif()
