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

o2physics_add_dpl_workflow(hyperon-reco-test
    SOURCES hyperon-reco-test.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derivedlambdakzeroanalysis
    SOURCES derivedlambdakzeroanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeroanalysis-mc
    SOURCES lambdakzeroanalysisMC.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeanalysis
    SOURCES cascadeanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeanalysismc
    SOURCES cascadeanalysisMC.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0postprocessing
    SOURCES v0postprocessing.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadecorrelations
    SOURCES cascadecorrelations.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(non-prompt-cascade
    SOURCES nonPromptCascade.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2::ReconstructionDataFormats O2Physics::AnalysisCore O2::DetectorsBase O2::DetectorsVertexing O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(k0mixedevents
    SOURCES k0_mixed_events.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(xi-lambda-corr-task
    SOURCES xiLambdaCorr.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(vzero-cascade-absorption
    SOURCES vzero_cascade_absorption.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derivedcascadeanalysis
    SOURCES derivedcascadeanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascpostprocessing
    SOURCES cascpostprocessing.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hstrangecorrelation
    SOURCES hStrangeCorrelation.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sigmaanalysis
    SOURCES sigmaanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phik0shortanalysis
    SOURCES phik0shortanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phi-strange-correlation
    SOURCES phiStrangeCorrelation.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdapolarization
    SOURCES lambdapolarization.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdapolsp
    SOURCES lambdapolsp.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdak0seff
    SOURCES lambdak0seff.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-lambda-spin-corr
    SOURCES taskLambdaSpinCorr.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascpolsp
    SOURCES cascpolsp.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

if(FastJet_FOUND)
o2physics_add_dpl_workflow(strangeness-in-jets
    SOURCES strangenessInJets.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangeness-in-jets-ions
    SOURCES strangenessInJetsIons.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sjet-tree-creator
    SOURCES sjetTreeCreator.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib
    COMPONENT_NAME Analysis)
endif()

o2physics_add_dpl_workflow(v0topologicalcuts
    SOURCES v0topologicalcuts.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0ptinvmassplots
    SOURCES v0ptinvmassplots.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derivedupcanalysis
    SOURCES derivedupcanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdak0sflattenicity
    SOURCES lambdak0sflattenicity.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdalambda
    SOURCES lambdalambda.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdajetpolarization
    SOURCES lambdaJetpolarization.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdaspincorrderived
    SOURCES lambdaspincorrderived.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangenessderivedbinnedinfo
    SOURCES strangenessderivedbinnedinfo.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdatwopartpolarization
    SOURCES lambdaTwoPartPolarization.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeanalysislightions
    SOURCES cascadeAnalysisLightIonsDerivedData.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strange-casc-track
    SOURCES strangeCascTrack.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)
