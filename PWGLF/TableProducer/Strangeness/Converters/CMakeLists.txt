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

o2physics_add_dpl_workflow(stradautrackstofpidconverter
                    SOURCES stradautrackstofpidconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stradautracksextraconverter2
                    SOURCES stradautracksextraconverter2.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strarawcentsconverter
                    SOURCES strarawcentsconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strarawcentsconverter2v4
                    SOURCES strarawcentsconverter2v4.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter
                    SOURCES straevselsconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter2
                    SOURCES straevselsconverter2.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::ITStracking
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter3
                    SOURCES straevselsconverter3.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::ITStracking
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter4
                    SOURCES straevselsconverter4.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::ITStracking
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter2rawcents
                    SOURCES straevselsconverter2rawcents.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter2rawcents2
                    SOURCES straevselsconverter2rawcents2.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter2rawcents3
                    SOURCES straevselsconverter2rawcents3.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0coresconverter
                    SOURCES v0coresconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0coresconverter2
                    SOURCES v0coresconverter2.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0mlscoresconverter
                    SOURCES v0mlscoresconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stramccollisionconverter
                    SOURCES stramccollisionconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strastampsconverter
                    SOURCES strastampsconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stracentconverter
                    SOURCES stracentconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stramccollmultconverter
                    SOURCES stramccollmultconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
