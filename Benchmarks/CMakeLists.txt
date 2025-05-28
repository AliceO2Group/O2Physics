# Copyright 2019-2020 CERN and copyright holders of ALICE O2. See
# https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public License
# v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

o2physics_add_dpl_workflow(produce-dynamic-extension
                           SOURCES src/produce-dynamic-extension.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME AnalysisBenchmark)

o2physics_add_dpl_workflow(produce-dynamic-extension-func
                           SOURCES src/produce-dynamic-extension-func.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME AnalysisBenchmark)

o2physics_add_dpl_workflow(produce-expression
                           SOURCES src/produce-expression.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME AnalysisBenchmark)

o2physics_add_dpl_workflow(produce-dynamic-extension-calib
                           SOURCES src/produce-dynamic-extension-calib.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                           COMPONENT_NAME AnalysisBenchmark)

o2physics_add_dpl_workflow(produce-expression-calib
                           SOURCES src/produce-expression-calib.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                           COMPONENT_NAME AnalysisBenchmark)
