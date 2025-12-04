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

o2physics_add_dpl_workflow(task-electron-weak-boson
                    SOURCES taskElectronWeakBoson.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore KFParticle::KFParticle O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-muon-charm-beauty-separation
                    SOURCES taskMuonCharmBeautySeparation.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-single-electron
                    SOURCES taskSingleElectron.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-single-muon
                    SOURCES taskSingleMuon.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-single-muon-mult
                    SOURCES taskSingleMuonMult.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-single-muon-reader
                    SOURCES taskSingleMuonReader.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-single-muon-reader-assoc
                    SOURCES taskSingleMuonReaderAssoc.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-single-muon-source
                    SOURCES taskSingleMuonSource.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
