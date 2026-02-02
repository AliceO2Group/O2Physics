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

o2physics_add_dpl_workflow(central-event-filter-task
    SOURCES cefpTask.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DataFormatsCTP
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(selected-bc-range-task
    SOURCES selectBCRange.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(nuclei-filter
    SOURCES PWGLF/nucleiFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::TOFBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(diffraction-filter
    SOURCES PWGUD/diffractionFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::DGCutparHolder
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(diffraction-bcfilter
    SOURCES PWGUD/diffractionBCFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::DGCutparHolder
    COMPONENT_NAME Analysis)

o2physics_add_header_only_library(HFFilterHelpers
    HEADERS PWGHF/HFFilterHelpers.h)

o2physics_add_dpl_workflow(hf-filter
    SOURCES PWGHF/HFFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::HFFilterHelpers O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-filter-charm-hadron-signals
    SOURCES PWGHF/HFFilterCharmHadronSignals.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::HFFilterHelpers
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hf-filter-prepare-ml-samples
    SOURCES PWGHF/HFFilterPrepareMLSamples.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::HFFilterHelpers O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(cf-filter
    SOURCES PWGCF/CFFilterAll.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore KFParticle::KFParticle O2::ReconstructionDataFormats O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(je-filter
    SOURCES PWGJE/jetFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(je-hf-filter
    SOURCES PWGJE/jetHFFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fje-filter
    SOURCES PWGJE/fullJetFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGJECore FastJet::FastJet FastJet::Contrib Boost::system
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lf-strangeness-filter
    SOURCES PWGLF/strangenessFilter.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore KFParticle::KFParticle O2::ReconstructionDataFormats O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mult-filter
    SOURCES PWGMM/multFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(em-global-dimuon-filter
    SOURCES PWGEM/globalDimuonFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::GlobalTracking
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(em-photon-filter
    SOURCES PWGEM/EMPhotonFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase O2Physics::PWGEMPhotonMesonCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(em-photon-filter-qc
    SOURCES PWGEM/EMPhotonFilterQC.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase O2Physics::PWGEMPhotonMesonCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(heavy-neutral-meson-filter
    SOURCES PWGEM/HeavyNeutralMesonFilter.cxx
    PUBLIC_LINK_LIBRARIES O2::Framework O2::EMCALBase O2::EMCALCalib O2Physics::AnalysisCore O2Physics::PWGEMPhotonMesonCore
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(lf-f1proton-filter
    SOURCES PWGLF/filterf1proton.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore KFParticle::KFParticle O2::DetectorsBase
    COMPONENT_NAME Analysis)

 o2physics_add_dpl_workflow(lf-doublephi-filter
    SOURCES PWGLF/filterdoublephi.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(h2fromlbfilter
    SOURCES PWGHF/H2fromLbFilter.cxx
    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
    COMPONENT_NAME Analysis)
