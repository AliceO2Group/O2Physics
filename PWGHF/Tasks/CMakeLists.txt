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

o2physics_add_dpl_workflow(task-charm-had-impact-par
                    SOURCES taskCharmHadImpactPar.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# o2physics_add_dpl_workflow(task-lc-centrality
#                     SOURCES taskLcCentrality.cxx
#                     PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
#                     COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-mc-efficiency
                    SOURCES taskMcEfficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-mc-efficiency-to-xi-pi
                    SOURCES taskMcEfficiencyToXiPi.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-mc-gen-pt-rap-shapes
                    SOURCES taskMcGenPtRapShapes.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-mc-validation
                    SOURCES taskMcValidation.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-multiplicity-estimator-correlation
                    SOURCES taskMultiplicityEstimatorCorrelation.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-pid-studies
                    SOURCES taskPidStudies.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# o2physics_add_dpl_workflow(task-sel-optimisation
#                     SOURCES taskSelOptimisation.cxx
#                     PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
#                     COMPONENT_NAME Analysis)
