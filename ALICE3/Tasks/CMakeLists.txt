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

o2physics_add_dpl_workflow(alice3-qa-multiplicity
                    SOURCES alice3-qa-multiplicity.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-qa-singleparticle
                    SOURCES alice3-qa-singleparticle.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-lutmaker
                    SOURCES alice3-lutmaker.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework FairRoot::Base O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-pid-rich-qa
                    SOURCES pidRICHqa.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-ecal-qa
                    SOURCES ECALqa.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-pid-ftof-qa
                    SOURCES pidFTOFqa.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-pid-separation-power
                    SOURCES alice3SeparationPower.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-cdeuteron
                    SOURCES alice3-cdeuteron.cxx
                    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-dilepton
                    SOURCES alice3-dilepton.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::FrameworkPhysicsSupport
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-taskcorrelationddbar
                    SOURCES alice3-taskcorrelationDDbar.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-multicharm
                    SOURCES alice3-multicharm.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-hf-task-3prong
                    SOURCES alice3HfTask3Prong.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-efficiency
                    SOURCES alice3Efficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-tracking-performance
                    SOURCES alice3TrackingPerformance.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-pid-evaluation
                    SOURCES alice3PidEvaluation.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-strangeness
                    SOURCES alice3-strangeness.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(alice3-dq-efficiency
                    SOURCES alice3-dq-efficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)