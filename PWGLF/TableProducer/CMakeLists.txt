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

# General purpose
o2physics_add_dpl_workflow(tpcpid
                    SOURCES lfTPCPID.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Strangeness
o2physics_add_dpl_workflow(lambdakzerobuilder
                    SOURCES lambdakzerobuilder.cxx
                    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzerolabelbuilder
                    SOURCES lambdakzerolabelbuilder.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzerofinder
                    SOURCES lambdakzerofinder.cxx
                    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzeromcfinder
                    SOURCES lambdakzeromcfinder.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lambdakzerotofpid
                    SOURCES lambdakzerotofpid.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadebuilder
                    SOURCES cascadebuilder.cxx
                    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore KFParticle::KFParticle
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadelabelbuilder
                    SOURCES cascadelabelbuilder.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascadefinder
                    SOURCES cascadefinder.cxx
                    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascademcfinder
                    SOURCES cascademcfinder.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hstrangecorrelationfilter
                    SOURCES hStrangeCorrelationFilter.cxx
                    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(st-coll-ids
                    SOURCES stCollIds.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Nuclei
o2physics_add_dpl_workflow(nucleustreecreator
                    SOURCES LFTreeCreatorNuclei.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-spectra
                    SOURCES nucleiSpectra.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lithium4analysis
                    SOURCES lithium4analysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Hypernuclei
o2physics_add_dpl_workflow(hypertriton-reco-task
                    SOURCES hyperRecoTask.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DCAFitter
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hyhefour-builder
                    SOURCES hyhe4builder.cxx
                    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(threebodymcfinder
                    SOURCES threebodymcfinder.cxx
                    PUBLIC_LINK_LIBRARIES O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertriton3bodybuilder
                    SOURCES hypertriton3bodybuilder.cxx
                    PUBLIC_LINK_LIBRARIES  O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hypertriton3bodyfinder
                    SOURCES hypertriton3bodyfinder.cxx
                    PUBLIC_LINK_LIBRARIES  O2::DCAFitter O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Resonances
o2physics_add_dpl_workflow(reso2initializer
                    SOURCES LFResonanceInitializer.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(f1protoninitializer
                    SOURCES f1protonInitializer.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(filterf1proton
                    SOURCES filterf1proton.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(v0qaanalysis
                    SOURCES v0qaanalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cascqaanalysis
                    SOURCES cascqaanalysis.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

# Spectra
o2physics_add_dpl_workflow(spectra-derived
                    SOURCES spectraDerivedMaker.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
