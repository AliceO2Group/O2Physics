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

o2physics_add_library(JCorran
                    SOURCES JFFlucAnalysis.cxx
                    JFFlucAnalysisO2Hist.cxx
                    FlowJSPCAnalysis.cxx
                    FlowJHistManager.cxx
                    JEPFlowAnalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore)

o2physics_target_root_dictionary(JCorran
                    HEADERS JFFlucAnalysis.h
                    JFFlucAnalysisO2Hist.h
                    FlowJHistManager.h
                    FlowJSPCAnalysis.h
                    FlowJSPCObservables.h
                    JEPFlowAnalysis.h
                    LINKDEF JCORRANLinkDef.h)
