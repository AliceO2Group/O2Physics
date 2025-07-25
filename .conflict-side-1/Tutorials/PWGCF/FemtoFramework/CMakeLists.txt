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

o2physics_add_dpl_workflow(femtodream-tutorial-1
                    SOURCES src/CFTutorialTask1.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(femtodream-tutorial-2
                    SOURCES src/CFTutorialTask2.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(femtodream-tutorial-3
                    SOURCES src/CFTutorialTask3.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(femtodream-tutorial-4
                    SOURCES src/CFTutorialTask4.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME AnalysisTutorial)

o2physics_add_dpl_workflow(femtodream-tutorial-5
                    SOURCES src/CFTutorialTask5.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME AnalysisTutorial)
