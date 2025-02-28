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

o2physics_add_dpl_workflow(decay3bodybuilder
    SOURCES decay3bodybuilder.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter KFParticle::KFParticle O2Physics::AnalysisCore O2::TOFBase O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hyhefour-builder
    SOURCES hyhe4builder.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertriton-reco-task
    SOURCES hyperRecoTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lnn-reco-task
    SOURCES lnnRecoTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertriton3bodyfinder
    SOURCES hypertriton3bodyfinder.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nucleustreecreator
    SOURCES LFTreeCreatorNuclei.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(he3hadronfemto
    SOURCES he3HadronFemto.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-spectra
    SOURCES nucleiSpectra.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(spectra-derived
    SOURCES spectraDerivedMaker.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(threebodymcfinder
    SOURCES threebodymcfinder.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(threebody-reco-task
    SOURCES threebodyRecoTask.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore O2Physics::EventFilteringUtils
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(ebye-maker
    SOURCES ebyeMaker.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cluster-studies-tree-creator
    SOURCES LFTreeCreatorClusterStudies.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pidtof-generic
    SOURCES pidTOFGeneric.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::TOFBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(threebody-kf-task
    SOURCES threebodyKFTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypernuclei-kf-reco-task
    SOURCES hypKfRecoTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore KFParticle::KFParticle
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypernuclei-kf-tree-creator
    SOURCES hypKfTreeCreator.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tr-he-analysis
    SOURCES trHeAnalysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(reduced3body-creator
    SOURCES reduced3bodyCreator.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore KFParticle::KFParticle O2Physics::EventFilteringUtils O2::TOFBase
    COMPONENT_NAME Analysis)
