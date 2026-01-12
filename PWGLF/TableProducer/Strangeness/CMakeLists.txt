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

add_subdirectory(Converters)

o2physics_add_dpl_workflow(cascadefinder
    SOURCES cascadefinder.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeflow
    SOURCES cascadeflow.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascademcbuilder
    SOURCES cascademcbuilder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascademcfinder
    SOURCES cascademcfinder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangenesstofpid
    SOURCES strangenesstofpid.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadespawner
    SOURCES cascadespawner.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascqaanalysis
    SOURCES cascqaanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hstrangecorrelationfilter
    SOURCES hStrangeCorrelationFilter.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(double-casc-tree-creator
    SOURCES doubleCascTreeCreator.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzerofinder
    SOURCES lambdakzerofinder.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeromcbuilder
    SOURCES lambdakzeromcbuilder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeromcfinder
    SOURCES lambdakzeromcfinder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzerospawner
    SOURCES lambdakzerospawner.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sigmaminus-task
    SOURCES sigmaminustask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strange-tree-creator
    SOURCES strangeTreeCreator.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(st-coll-ids
    SOURCES stCollIds.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangederivedbuilder
    SOURCES strangederivedbuilder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangenessbuilder
    SOURCES strangenessbuilder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle O2Physics::TPCDriftManager O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0-selector
    SOURCES v0selector.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0qaanalysis
    SOURCES v0qaanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdalambdatable
    SOURCES LambdaLambdatable.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

# ML selection
o2physics_add_dpl_workflow(lambdakzeromlselectiontreecreator
    SOURCES lambdakzeroMLSelectionTreeCreator.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascademlselectiontreecreator
    SOURCES cascadeMLSelectionTreeCreator.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeromlselection
    SOURCES lambdakzeromlselection.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascademlselection
    SOURCES cascademlselection.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sigma0builder
    SOURCES sigma0builder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdajetpolarizationbuilder
    SOURCES lambdaJetpolarizationbuilder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stracents
    SOURCES stracents.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdaspincorrelation
    SOURCES lambdaspincorrelation.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

if(FastJet_FOUND)
    o2physics_add_dpl_workflow(lfinjets
        SOURCES lfinjets.cxx
        PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib O2Physics::EventFilteringUtils
        COMPONENT_NAME Analysis)
endif()

o2physics_add_dpl_workflow(phi-strange-correlator
    SOURCES phiStrangeCorrelator.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)
