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

o2physics_add_library(UDGoodRunSelector
                      SOURCES UDGoodRunSelector.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ROOT::EG RapidJSON::RapidJSON)

o2physics_target_root_dictionary(UDGoodRunSelector
                                 HEADERS UDGoodRunSelector.h
                                 LINKDEF UDGoodRunSelectorLinkDef.h)

o2physics_add_library(UDFSParser
                      SOURCES UDFSParser.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(UDFSParser
                                 HEADERS UDFSParser.h
                                 LINKDEF UDFSParserLinkDef.h)

o2physics_add_library(DGCutparHolder
                      SOURCES DGCutparHolder.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(DGCutparHolder
                                 HEADERS DGCutparHolder.h
                                 LINKDEF DGCutparHolderLinkDef.h)

o2physics_add_library(SGCutParHolder
                      SOURCES SGCutParHolder.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(SGCutParHolder
                      HEADERS SGCutParHolder.h
                      LINKDEF SGCutParHolderLinkDef.h)

o2physics_add_library(DGPIDSelector
                      SOURCES DGPIDSelector.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ROOT::EG)

o2physics_target_root_dictionary(DGPIDSelector
                                 HEADERS DGPIDSelector.h
                                 LINKDEF DGPIDSelectorLinkDef.h)

o2physics_add_library(UPCCutparHolder
                      SOURCES UPCCutparHolder.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(UPCCutparHolder
                                 HEADERS UPCCutparHolder.h
                                 LINKDEF UPCCutparHolderLinkDef.h)

o2physics_add_library(decayTree
                      SOURCES decayTree.cxx
                      PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ROOT::EG RapidJSON::RapidJSON)

o2physics_target_root_dictionary(decayTree
                                 HEADERS decayTree.h
                                 LINKDEF decayTreeLinkDef.h)
