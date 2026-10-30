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
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadeflow
    SOURCES cascadeflow.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascademcbuilder
    SOURCES cascademcbuilder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascademcfinder
    SOURCES cascademcfinder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangenesstofpid
    SOURCES strangenesstofpid.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadespawner
    SOURCES cascadespawner.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascqaanalysis
    SOURCES cascqaanalysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascderivedqaanalysis
    SOURCES cascderivedqaanalysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hstrangecorrelationfilter # o2-linter: disable=name/o2-workflow (Keep historical workflow name for backward compatibility.)
    SOURCES hStrangeCorrelationFilter.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(double-omega-tree-creator
    SOURCES doubleOmegaTreeCreator.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(double-omega-builder
    SOURCES doubleOmegaBuilder.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzerofinder
    SOURCES lambdakzerofinder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeromcbuilder
    SOURCES lambdakzeromcbuilder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeromcfinder
    SOURCES lambdakzeromcfinder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzerospawner
    SOURCES lambdakzerospawner.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sigmaminus-task # o2-linter: disable=name/o2-workflow (Keep historical workflow name for backward compatibility.)
    SOURCES sigmaminustask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sigmaplusbuilder
    SOURCES sigmaplusbuilder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strange-tree-creator
    SOURCES strangeTreeCreator.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(st-coll-ids
    SOURCES stCollIds.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangederivedbuilder
    SOURCES strangederivedbuilder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangenessbuilder
    SOURCES strangenessbuilder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter KFParticle::KFParticle O2Physics::TPCDriftManager O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0-selector # o2-linter: disable=name/o2-workflow (Keep historical workflow name for backward compatibility.)
    SOURCES v0selector.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0qaanalysis
    SOURCES v0qaanalysis.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdalambdatable # o2-linter: disable=name/o2-workflow (Keep historical workflow name for backward compatibility.)
    SOURCES LambdaLambdatable.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
    COMPONENT_NAME Analysis)

# ML selection
o2physics_add_dpl_workflow(lambdakzeromlselectiontreecreator # o2-linter: disable=name/o2-workflow (Keep historical workflow name for backward compatibility.)
    SOURCES lambdakzeroMLSelectionTreeCreator.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascademlselectiontreecreator # o2-linter: disable=name/o2-workflow (Keep historical workflow name for backward compatibility.)
    SOURCES cascadeMLSelectionTreeCreator.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeromlselection
    SOURCES lambdakzeromlselection.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascademlselection
    SOURCES cascademlselection.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sigma0builder
    SOURCES sigma0builder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2Physics::AnalysisCore O2Physics::MLCore O2Physics::AnalysisCCDB O2::DCAFitter
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdajetpolarizationbuilder # o2-linter: disable=name/o2-workflow (Keep historical workflow name for backward compatibility.)
    SOURCES lambdaJetpolarizationbuilder.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdajetpolarizationions # o2-linter: disable=name/o2-workflow (Keep historical workflow name for backward compatibility.)
    SOURCES lambdaJetPolarizationIons.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet O2Physics::AnalysisCCDB
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(stracents
    SOURCES stracents.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdaspincorrelation
    SOURCES lambdaspincorrelation.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sigmahadcorr # o2-linter: disable=name/o2-workflow (Keep historical workflow name for backward compatibility.)
    SOURCES sigmaHadCorr.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

if(FastJet_FOUND)
    o2physics_add_dpl_workflow(lfinjets
        SOURCES lfinjets.cxx
        REUSE_FROM LfStrangenessPCH
        PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib O2Physics::EventFilteringUtils
        COMPONENT_NAME Analysis)
endif()

o2physics_add_dpl_workflow(phi-strange-correlator
    SOURCES phiStrangeCorrelator.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(strangeness-mc-centrality
    SOURCES strangenessMcCentrality.cxx
    REUSE_FROM LfStrangenessPCH
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)
