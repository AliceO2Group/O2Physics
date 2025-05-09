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

o2physics_add_library(PWGDQCore
               SOURCES  VarManager.cxx
                        HistogramManager.cxx
                        HistogramsLibrary.cxx
                        CutsLibrary.cxx
                        MixingLibrary.cxx
                        MCSignalLibrary.cxx
                        MixingHandler.cxx
                        AnalysisCut.cxx
                        AnalysisCompositeCut.cxx
                        MCProng.cxx
                        MCSignal.cxx
               PUBLIC_LINK_LIBRARIES O2::Framework O2::DCAFitter O2::GlobalTracking O2Physics::AnalysisCore KFParticle::KFParticle)

o2physics_target_root_dictionary(PWGDQCore
                          HEADERS   AnalysisCut.h
                                    AnalysisCompositeCut.h
                                    VarManager.h
                                    HistogramManager.h
                                    CutsLibrary.h
                                    MixingHandler.h
                                    MixingLibrary.h
                                    HistogramsLibrary.h
                                    MCProng.h
                                    MCSignal.h
                                    MCSignalLibrary.h
                          LINKDEF PWGDQCoreLinkDef.h)
