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
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hyhefour-builder
    SOURCES hyhe4builder.cxx
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hyper-kink-reco-task
    SOURCES hyperKinkRecoTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertriton-reco-task
    SOURCES hyperRecoTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
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

o2physics_add_dpl_workflow(lithium4analysis
    SOURCES lithium4analysis.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-spectra
    SOURCES nucleiSpectra.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
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
    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)
