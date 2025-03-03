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

o2physics_add_dpl_workflow(meanpt-fluctuations
                    SOURCES MeanptFluctuations.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mean-pt-fluc-id
                    SOURCES meanPtFlucId.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(netproton-cumulants
                    SOURCES NetProtonCumulants.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(netproton-cumulants-mc
                    SOURCES netprotonCumulantsMc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(identified-meanpt-fluctuations
                    SOURCES IdentifiedMeanPtFluctuations.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(robust-fluctuation-observables
                    SOURCES RobustFluctuationObservables.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0pt-pp-task
                    SOURCES Diff_pT_fluct_PID.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(factorial-moments
                    SOURCES FactorialMomentsTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(kaon-isospin-fluctuations
                    SOURCES kaonIsospinFluctuations.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(netcharge-fluctuations
                    SOURCES netchargeFluctuations.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nch-cumulants-id
                    SOURCES nchCumulantsId.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)
