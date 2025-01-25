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

o2physics_add_dpl_workflow(v0converter
                    SOURCES v0converter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mc-converter
                    SOURCES mcConverter.cxx
                    PUBLIC_LINK_LIBRARIES
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tracks-extra-converter
                    SOURCES tracksExtraConverter.cxx
                    PUBLIC_LINK_LIBRARIES
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tracks-extra-v002-converter
                    SOURCES tracksExtraV002Converter.cxx
                    PUBLIC_LINK_LIBRARIES
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mft-tracks-converter
                    SOURCES mftTracksConverter.cxx
                    PUBLIC_LINK_LIBRARIES
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fdd-converter
                    SOURCES fddConverter.cxx
                    PUBLIC_LINK_LIBRARIES
                COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(collision-converter
                    SOURCES collisionConverter.cxx
                    PUBLIC_LINK_LIBRARIES
                COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mccollision-converter
                SOURCES mcCollisionConverter.cxx
                PUBLIC_LINK_LIBRARIES
            COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(bc-converter
                    SOURCES bcConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(bc-flags-creator
                    SOURCES bcFlagsCreator.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(calo-label-converter
                    SOURCES caloLabelConverter.cxx
                    PUBLIC_LINK_LIBRARIES
                COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(zdc-converter
                    SOURCES zdcConverter.cxx
                    PUBLIC_LINK_LIBRARIES
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hmpid-converter
                    SOURCES hmpConverter.cxx
                    PUBLIC_LINK_LIBRARIES
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multsextra-converter
                    SOURCES multsExtraConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multmcextras-converter
                    SOURCES multMCExtrasConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(trackqa-converter
                    SOURCES trackQAConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(trackqa-converter-002
                    SOURCES trackQA002Converter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(run2bcinfos-converter
                    SOURCES run2bcinfosConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(run2-tracks-extra-converter
                    SOURCES run2tracksExtraConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(run2-tiny-to-full-pid
                    SOURCES run2TinyToFullPID.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
