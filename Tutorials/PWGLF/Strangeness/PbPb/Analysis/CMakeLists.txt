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

o2physics_add_dpl_workflow(strangeness-pbpb-skeleton
                  SOURCES strangeness_pbpb_skeleton.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(strangeness-pbpb-step0
                  SOURCES strangeness_pbpb_step0.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(strangeness-pbpb-step1
                  SOURCES strangeness_pbpb_step1.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(strangeness-pbpb-step2
                  SOURCES strangeness_pbpb_step2.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(strangeness-pbpb-step3
                  SOURCES strangeness_pbpb_step3.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(strangeness-pbpb-step4
                  SOURCES strangeness_pbpb_step4.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

