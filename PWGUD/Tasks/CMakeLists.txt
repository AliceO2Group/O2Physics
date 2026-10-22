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

o2physics_add_dpl_workflow(upc-polarisation-jpsi-incoh
                           SOURCES upcPolarisationJpsiIncoh.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-track-vertexing-qa
                           SOURCES upcTrackVertexingQA.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-analysis
                           SOURCES upcAnalysis.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-spectra
                           SOURCES sgSpectraAnalyzer.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-pid-spectra-table
                           SOURCES sgPIDSpectraTable.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-pid-analyzer
                           SOURCES sgPIDAnalyzer.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-pid-spectra
                           SOURCES sgPIDSpectra.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-incl-jpsi
                           SOURCES sgInclJpsi.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-excl-omega
                           SOURCES sgExclOmega.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-excl-universe
                           SOURCES sgExcUniverse.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-fourpi
                           SOURCES sgFourPiAnalyzer.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-sixpi
                           SOURCES sgSixPiAnalyzer.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-twopi
                           SOURCES sgTwoPiAnalyzer.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-d0
                           SOURCES sgD0Analyzer.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-fit
                           SOURCES sgFITAnalyzer.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-forward
                           SOURCES upcForward.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(diff-mcdatascanner
                           SOURCES diffMCDataScanner.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(diff-mcqa
                           SOURCES diffMCQA.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGCutparHolder
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(diff-qa
                           SOURCES diffQA.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGCutparHolder
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(dg-cand-analyzer
                           SOURCES dgCandAnalyzer.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGCutparHolder O2Physics::UDGoodRunSelector O2Physics::DGPIDSelector O2Physics::UDFSParser
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-cand-analyzer
                           SOURCES upcCandidateAnalyzer.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-cand-producer-qa
                           SOURCES upcCandidateProducerQa.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-mft
                           SOURCES upcMft.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsBase O2::DetectorsCommonDataFormats
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-veto-analysis
                           SOURCES upcVetoAnalysis.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(tautau13topo
                           SOURCES upcTauTau13topo.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-tau-rl
                           SOURCES upcTauRl.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsBase O2::DetectorsCommonDataFormats
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(polarisation-rho
                           SOURCES polarisationRho.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-jpsi-corr
                           SOURCES upcJpsiCorr.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(exclusive-phi
                           SOURCES exclusivePhi.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-photonuclear-jmg
                           SOURCES upcPhotonuclearAnalysisJMG.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore O2::ReconstructionDataFormats O2::DetectorsBase O2::DetectorsCommonDataFormats
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(exclusive-two-protons
                           SOURCES exclusiveTwoProtons.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(exclusive-two-protons-sg
                           SOURCES exclusiveTwoProtonsSG.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(exclusive-phi-leptons
                           SOURCES exclusivePhiLeptons.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(exclusive-phi-leptons-trees
                           SOURCES exclusivePhiLeptonsTrees.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(exclusive-pentaquark
                           SOURCES exclusivePentaquark.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(inclusive-phikstar-sd
                           SOURCES sginclusivePhiKstarSD.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-exclusive-phi
                           SOURCES sgExclusivePhi.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-exclusive-phi-its
                           SOURCES sgExclusivePhiITSselections.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-pion-analysis
                           SOURCES upcPionAnalysis.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fwd-muons-upc
                           SOURCES fwdMuonsUpc.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-fwd-jpsi-rl
                           SOURCES upcFwdJpsiRl.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-semi-fwd-jpsi-rl
                           SOURCES upcSemiFwdJpsiRl.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-event-its-rof-counter
                           SOURCES upcEventItsRofCounter.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-rho-analysis
                           SOURCES upcRhoAnalysis.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-rho-prime-analysis
                           SOURCES upcRhoPrimeAnalysis.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-by-event
                           SOURCES eventByevent.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(exclusive-rho-to-four-pi
                           SOURCES exclusiveRhoTo4Pi.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::DGPIDSelector
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-quarkonia-central-barrel
                           SOURCES upcQuarkoniaCentralBarrel.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::DetectorsBase
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-test-mc-std-tabs-rl
                           SOURCES upcTestMcStdTabsRl.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::ReconstructionDataFormats O2::DetectorsBase O2::DetectorsCommonDataFormats
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(decaytree-analyzer
                           SOURCES decayTreeAnalyzer.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::UDGoodRunSelector O2Physics::decayTree
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-cumulants-upc
                           SOURCES flowCumulantsUpc.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-correlations-upc
                           SOURCES flowCorrelationsUpc.cxx
                           REUSE_FROM UDTablesPCH
			               PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::GFWCore O2Physics::PWGCFCore O2Physics::TwoPartCorrCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-mc-upc
                           SOURCES flowMcUpc.cxx
			               PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::GFWCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(flow-long-range-corr-upc
                           SOURCES flowLongRangeCorrUpc.cxx
                           REUSE_FROM UDTablesPCH
			               PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::AnalysisCCDB O2Physics::PWGCFCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(analysis-mc-dpm-jet-sg-v3
                           SOURCES analysisMCDPMJetSGv3.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-test-rct-tables
                           SOURCES upcTestRctTables.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-exclusive-jpsi-midrapidity
                           SOURCES sgExclusiveJpsiMidrapidity.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-vm-rof
                           SOURCES upcVmRof.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(upc-test-fit-bit-mapping
                           SOURCES upcTestFitBitMapping.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pt-spectra-inclusive-upc
                           SOURCES ptSpectraInclusiveUpc.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(sg-deuteron-spectra
                           SOURCES sgDeuteronSpectra.cxx
                           REUSE_FROM UDTablesPCH
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)
