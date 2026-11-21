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

o2physics_add_dpl_workflow(stradautracksconverter
                    SOURCES stradautracksconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stradautrackstpcpidconverter
                    SOURCES stradautrackstpcpidconverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stradautrackstofpidconverter
                    SOURCES stradautrackstofpidconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stradautrackstofpidconverter2
                    SOURCES stradautrackstofpidconverter2.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stradautrackstofpidconverter3
                    SOURCES stradautrackstofpidconverter3.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stradautrackstofpidconverter4
                    SOURCES stradautrackstofpidconverter4.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stradautracksextraconverter
                    SOURCES stradautracksextraconverter.cxx
                    REUSE_FROM LfStrangenessPCH
		                PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stradautracksextraconverter2
                    SOURCES stradautracksextraconverter2.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stradautracksextraconverter3
                    SOURCES stradautracksextraconverter3.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strarawcentsconverter
                    SOURCES strarawcentsconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strarawcentsconverter2v4
                    SOURCES strarawcentsconverter2v4.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter
                    SOURCES straevselsconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter2
                    SOURCES straevselsconverter2.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter3
                    SOURCES straevselsconverter3.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter4
                    SOURCES straevselsconverter4.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter5
                    SOURCES straevselsconverter5.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter6
                    SOURCES straevselsconverter6.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselextrasconverter
                    SOURCES straevselextrasconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselextrasconverter2
                    SOURCES straevselextrasconverter2.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter2rawcents
                    SOURCES straevselsconverter2rawcents.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter2rawcents2
                    SOURCES straevselsconverter2rawcents2.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(straevselsconverter2rawcents3
                    SOURCES straevselsconverter2rawcents3.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0coresconverter
                    SOURCES v0coresconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0coresconverter2
                    SOURCES v0coresconverter2.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0mlscoresconverter
                    SOURCES v0mlscoresconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stramccollisionconverter
                    SOURCES stramccollisionconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stramccollisionconverter2
                    SOURCES stramccollisionconverter2.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strastampsconverter
                    SOURCES strastampsconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stracentconverter
                    SOURCES stracentconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stracentconverter2
                    SOURCES stracentconverter2.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stramccollmultconverter
                    SOURCES stramccollmultconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stramccollmultconverter2
                    SOURCES stramccollmultconverter2.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(zdcneutronsconverter
                    SOURCES zdcneutronsconverter.cxx
                    REUSE_FROM LfStrangenessPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
