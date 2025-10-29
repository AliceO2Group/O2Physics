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

# ITS
o2physics_add_dpl_workflow(pid-its
    SOURCES pidITS.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

# TOF
o2physics_add_dpl_workflow(pid-tof-base
                    SOURCES pidTOFBase.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::TOFBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tof
                    SOURCES pidTOF.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::TOFWorkflowUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tof-merge
                    SOURCES pidTOFMerge.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::TOFWorkflowUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tof-beta
                    SOURCES pidTOFbeta.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::TOFWorkflowUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tof-full
                    SOURCES pidTOFFull.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::TOFWorkflowUtils
                    COMPONENT_NAME Analysis)

# TPC

o2physics_add_dpl_workflow(pid-tpc-service
                    SOURCES pidTPCService.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tpc-service-run2
                    SOURCES pidTPCServiceRun2.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tpc-base
                    SOURCES pidTPCBase.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tpc
                    SOURCES pidTPC.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase O2Physics::AnalysisCCDB O2Physics::MLCore
                    COMPONENT_NAME Analysis)

# HMPID

# BAYES

o2physics_add_dpl_workflow(pid-bayes
                    SOURCES pidBayes.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
