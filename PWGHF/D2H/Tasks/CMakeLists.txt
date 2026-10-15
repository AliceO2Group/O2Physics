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

o2physics_add_dpl_workflow(task-b0
                    SOURCES taskB0.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-b0-reduced
                    SOURCES taskB0Reduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-b0-to-jpsi-k0-star-reduced
                    SOURCES taskB0ToJpsiK0StarReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-bplus
                    SOURCES taskBplus.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-bplus-reduced
                    SOURCES taskBplusReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-bplus-to-jpsi-k-reduced
                    SOURCES taskBplusToJpsiKReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-bs-reduced
                    SOURCES taskBsReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-bs-to-jpsi-phi-reduced
                    SOURCES taskBsToJpsiPhiReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-bs
                    SOURCES taskBs.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-cd
                    SOURCES taskCd.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-charm-polarisation
                    SOURCES taskCharmPolarisation.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-charm-reso-to-d-v0-reduced
                    SOURCES taskCharmResoToDV0Reduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-charm-reso-to-d-trk-reduced
                    SOURCES taskCharmResoToDTrkReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-d0
                    SOURCES taskD0.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::SGCutParHolder O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-directed-flow-charm-hadrons
                    SOURCES taskDirectedFlowCharmHadrons.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-dplus
                    SOURCES taskDplus.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::SGCutParHolder O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-ds
                    SOURCES taskDs.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-dstar-to-d0-pi
                    SOURCES taskDstarToD0Pi.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-flattenicity-d0-lc
                    SOURCES taskFlattenicityD0Lc.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::SGCutParHolder O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-flow-charm-hadrons
                    SOURCES taskFlowCharmHadrons.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-lb
                    SOURCES taskLb.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-lb-reduced
                    SOURCES taskLbReduced.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-lc
                    SOURCES taskLc.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::SGCutParHolder O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-upc-lc
                    SOURCES taskUpcLc.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::SGCutParHolder O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-lc-to-k0s-p
                    SOURCES taskLcToK0sP.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-net-charm-fluctuations
                    SOURCES taskNetCharmFluctuations.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-omegac0-to-omega-pi
                    SOURCES taskOmegac0ToOmegaPi.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-omegac0-to-omega-pi-qa
                    SOURCES taskOmegac0ToOmegaPiQa.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-pt-fluc-charm-hadrons
                    SOURCES taskPtFlucCharmHadrons.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-sigmac
                    SOURCES taskSigmac.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-sigmac-to-cascade
                    SOURCES taskSigmacToCascade.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-xic
                    SOURCES taskXic.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-xic-to-xi-pi-pi
                    SOURCES taskXicToXiPiPi.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-xic0-to-xi-pi
                    SOURCES taskXic0ToXiPi.cxx
                    REUSE_FROM HfCandidatesPCH
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-hidden-charm
                    SOURCES taskHiddenCharm.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-deuteron-from-lb
                    SOURCES taskDeuteronFromLb.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)
