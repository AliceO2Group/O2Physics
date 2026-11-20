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
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derivedlambdakzeroanalysis
    SOURCES derivedlambdakzeroanalysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeroanalysis-mc
    SOURCES lambdakzeroanalysisMC.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeanalysis
    SOURCES cascadeanalysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeanalysismc
    SOURCES cascadeanalysisMC.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0postprocessing
    SOURCES v0postprocessing.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadecorrelations
    SOURCES cascadecorrelations.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(non-prompt-cascade
    SOURCES nonPromptCascade.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::Framework O2::ReconstructionDataFormats O2Physics::AnalysisCore O2::DetectorsBase O2::DetectorsVertexing O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(k0mixedevents
    SOURCES k0_mixed_events.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(xi-lambda-corr-task
    SOURCES xiLambdaCorr.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(vzero-cascade-absorption
    SOURCES vzero_cascade_absorption.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derivedcascadeanalysis
    SOURCES derivedcascadeanalysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derivedxi0analysis
    SOURCES derivedXi0analysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascpostprocessing
    SOURCES cascpostprocessing.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hstrangecorrelation
    SOURCES hStrangeCorrelation.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sigmaanalysis
    SOURCES sigmaanalysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phik0shortanalysis
    SOURCES phik0shortanalysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phi-strange-correlation
    SOURCES phiStrangeCorrelation.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdapolarization
    SOURCES lambdapolarization.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdapolsp
    SOURCES lambdapolsp.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdak0seff
    SOURCES lambdak0seff.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(task-lambda-spin-corr
    SOURCES taskLambdaSpinCorr.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascpolsp
    SOURCES cascpolsp.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

if(FastJet_FOUND)
o2physics_add_dpl_workflow(strangeness-in-jets
    SOURCES strangenessInJets.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangeness-in-jets-ions
    SOURCES strangenessInJetsIons.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sjet-tree-creator
    SOURCES sjetTreeCreator.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib
    COMPONENT_NAME Analysis)
endif()

o2physics_add_dpl_workflow(v0topologicalcuts
    SOURCES v0topologicalcuts.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0ptinvmassplots
    SOURCES v0ptinvmassplots.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(derivedupcanalysis
    SOURCES derivedupcanalysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdak0sflattenicity
    SOURCES lambdak0sflattenicity.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdalambda
    SOURCES lambdalambda.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdajetpolarization
    SOURCES lambdaJetpolarization.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdajetpolarizationionsderived
    SOURCES lambdaJetPolarizationIonsDerived.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdaspincorrderived
    SOURCES lambdaspincorrderived.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangenessderivedbinnedinfo
    SOURCES strangenessderivedbinnedinfo.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambda-two-part-polarization
    SOURCES lambdaTwoPartPolarization.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeanalysislightions
    SOURCES cascadeAnalysisLightIonsDerivedData.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strange-casc-track
    SOURCES strangeCascTrack.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(zdccalderived
    SOURCES zdccalderived.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(forwardlambdakzeroanalysis
    SOURCES forwardlambdakzeroanalysis.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2::GlobalTracking O2::DCAFitter
    COMPONENT_NAME Analysis)
