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


o2physics_add_dpl_workflow(event-converter2
                  SOURCES eventConverter2.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-converter3
                  SOURCES eventConverter3.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-converter4
                  SOURCES eventConverter4.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-converter5
                  SOURCES eventConverter5.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(electron-converter2
                  SOURCES electronConverter2.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(electron-converter3
                  SOURCES electronConverter3.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(electron-converter4
                  SOURCES electronConverter4.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(electron-converter5
                  SOURCES electronConverter5.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(track-converter1
                  SOURCES trackConverter1.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(muon-converter1
                  SOURCES muonConverter1.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(muon-converter2
                  SOURCES muonConverter2.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(muon-selfid-converter1
                  SOURCES muonSelfIdConverter1.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mcparticle-converter1
                  SOURCES mcParticleConverter1.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qvec-converter
                    SOURCES qvecConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-norm-converter1
                  SOURCES eventNormConverter1.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-norm-converter2
                  SOURCES eventNormConverter2.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-mult-converter1
                  SOURCES eventMultConverter1.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-cent-converter1
                  SOURCES eventCentConverter1.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qvec2-converter0
                  SOURCES qvec2Converter0.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qvec3-converter0
                  SOURCES qvec3Converter0.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(bc-converter1
                  SOURCES bcConverter1.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

