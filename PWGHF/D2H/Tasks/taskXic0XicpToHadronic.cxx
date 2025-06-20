// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskXic0XicpToHadronic.cxx
/// \brief analysis task for Xic0 -> Xi Pi & Xicp -> Xi Pi Pi  
///
/// \author Jinhyun Park <jinhyun.park@cern.ch>, Pusan National University
/// \author Krista Smith <krista.lizbeth.smith@cern.ch>, Pusan National University

// Mandatory headers
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"

// Datamodels
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

// AOB
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Xic0Xicp analysis task
struct HfTaskXic0XicpToHadronic {

	struct : ConfigurableGroup {
		Configurable<int>selectionFlagXic0Xicp{"selectionFlagXic0Xicp", 1, "Selection flag for Xic0"};
		Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle's rapidity"};
		Configurable<float> yCandRecoMax{"yCandRecoMax", 0.8, "max. reco particle's rapidity"};
		Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track's pseudo-rapidity"};
		Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track's pT"};
		Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic0_xicp_to_hadronic::vecBinsPt}, "pT bin limits"};
		// MC checks
		Configurable<bool> checkDecayTypeMc{"checkDecayTypeMc", false, "Flag to enable Decaytype histogram"};
		// THnSparse for ML selection check
		Configurable<bool> enableTHn{"enableTHn", false, "Fill THnSparse for Xic"};
		
		// Axis
		ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {400, 0., 40.}, ""};
		ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {300, 1.8, 3.0}, ""};
		ConfigurableAxis thnConfigAxisPtProng{"thnConfigAxisPtProng", {300, 0., 30.}, ""};
		ConfigurableAxis thnConfigAxisChi2PCA{"thnConfigAxisChi2PCA", {200, 0., 20}, ""};
		ConfigurableAxis thnConfigAxisDecLength{"thnConfigAxisDecLength", {200, 0., 0.5}, ""};
		ConfigurableAxis thnConfigAxisDecLengthXY{"thnConfigAxisDecLengthXY", {200, 0., 0.5}, ""};
		ConfigurableAxis thnConfigAxisCPA{"thnConfigAxisCPA", {110, -1.1, 1.1}, ""};
		ConfigurableAxis thnConfigAxisBdtScoreBkg{"thnConfigAxisBdtScoreBkg", {100, 0., 1.}, ""};
		ConfigurableAxis thnConfigAxisBdtScorePrompt{"thnConfigAxisBdtScorePrompt", {100, 0., 1.}, ""};
		ConfigurableAxis thnConfigAxisBdtScoreNonPrompt{"thnConfigAxisBdtScoreNonPrompt", {100, 0., 1.}, ""};
		ConfigurableAxis binsDecayLength{"binsDecayLength", {200, 0., 0.5}, ""};
		ConfigurableAxis binsErrDecayLength{"binsErrDecayLength", {100, 0., 1.}, ""};
		ConfigurableAxis binsDCA{"binsDCA", {100, -0.05, 0.05}, ""};
		ConfigurableAxis binsImpParErr{"binsImpParErr", {200, -0.1, 0.1}, ""};
		ConfigurableAxis binsSV{"binsSV", {200, -5., 5.}, ""};
		ConfigurableAxis binsChi2{"binsChi2", {200, 0., 0.1}, ""};
	} configs;

	Service<o2::framework::O2DatabasePDG> pdg;
	HistogramRegistry registry{"registry"};

	void init(InitContext const&)
	{
		// array of enabled process functions
		std::array<bool, 8> doprocess{doprocessXic0WithDCAFitter,
									  doprocessXic0WithKFParticle,
									  doprocessXicpWithDCAFitter,
									  doprocessXicpWithKFParticle,
									  doprocessMcXic0WithDCAFitter,
									  doprocessMcXic0WithKFParticle,
									  doprocessMcXicpWithDCAFitter,
									  doprocessMcXicpWithKFParticle};

		// Check process function is enalbled.
		if ((std::accumulate(doprocess.begin(), doprocess.end(), 0))==0) {
			LOGP(fatal, "!!!No process function enabled!!!");
		}

		// Check if more than one sv reconstruction method is enabled
		if ((doprocess[0] || doprocess[2] || doprocess[4] || doprocess[6]) && (doprocess[1] || doprocess[3] || doprocess[5] || doprocess[7])) {
			LOGP(fatal, "!!!Cannot enable DCAFitter and KFParticle at the same time");	
		}

		static const AxisSpec axisMassXic0Xicp		= {300, 1.8, 3.0, "inv. mass (GeV/#it{c}^{2}"};
		static const AxisSpec axisMassXicpRes	= {300, 1.8, 3.0, "inv. mass (GeV/#it{c}^{2}"}; // -> This is for Xicp's resonance decay...
		static const AxisSpec axisPt				= {(std::vector<double>)configs.binsPt, "#it{p}_{T} (GeV/#it{c})"};
		static const AxisSpec axisDecayLength		= {configs.binsDecayLength};
		static const AxisSpec axisErrDecayLength	= {configs.binsErrDecayLength};
		static const AxisSpec axisDCA				= {configs.binsDCA};
		static const AxisSpec axisImpParErr			= {configs.binsImpParErr};
		static const AxisSpec axisSV				= {configs.binsSV};
		static const AxisSpec axisChi2				= {configs.binsChi2};

		// Defnie histograms for data analysis
		if (doprocess[0] || doprocess[1] || doprocess[2] || doprocess[3]) {
			// Common properties of Xic0 and Xicp
			registry.add("hMass", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; mass #Xi^{#mp} #pi^{#pm}; entries", {HistType::kTH2F, {axisMassXic0Xicp, axisPt}});
			registry.add("hPt", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{400, 0., 40.}}});
			registry.add("hEta", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; #it{#eta}; entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
			registry.add("hRapidity", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; #it{y}; entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
			registry.add("hCPA", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; Cosine of pointing angle; entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
			registry.add("hCPAXY", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; Cosine of pointing angle xy; entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
			registry.add("hDecayLength", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; decay length (cm); entries", {HistType::kTH2F, {axisDecayLength, axisPt}});
			registry.add("hErrDecayLength", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; decay length error (cm); entries", {HistType::kTH2F, {axisErrDecayLength, axisPt}});
			registry.add("hDecayLengthXY", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; decay length xy (cm); entries", {HistType::kTH2F, {axisDecayLength, axisPt}});
			registry.add("hErrDecayLengthXY", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; decay length xy error (cm); entries", {HistType::kTH2F, {axisErrDecayLength, axisPt}});
			registry.add("hSVx", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; secondary vertex position x (cm); entries", {HistType::kTH2F, {axisSV, axisPt}});
			registry.add("hSVy", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; secondary vertey position y (cm); entries", {HistType::kTH2F, {axisSV, axisPt}});
			registry.add("hSVz", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; secondary vertez position z (cm); entries", {HistType::kTH2F, {axisSV, axisPt}});
			registry.add("hCPAXi", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates;#Xi^{#minus} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
			registry.add("hCPAxyXi", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates;#Xi^{#minus} candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
			registry.add("hCPALambda", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates;#Lambda candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
			registry.add("hCPAxyLambda", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates;#Lambda candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
			registry.add("hPtProng0", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; #Xi^{#mp} #it{p}_{T}; entries", {HistType::kTH1F, {{200, 0., 20.}}});
			registry.add("hPtProng1", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; #pi^{#pm} #it{p}_{T}; entries", {HistType::kTH1F, {{200, 0., 20.}}});
			registry.add("hPtProng0VsPt", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; #Xi^{#mp} #it{p}_{T}; entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
			registry.add("hPtProng1VsPt", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates; #pi^{#pm} #it{p}_{T}; entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});      
			registry.add("hd0Prong0", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates;prong 0 (#Xi^{#mp}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
			registry.add("hd0Prong1", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates;prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
			registry.add("hImpParErr", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates;prongs impact parameter error (cm);entries", {HistType::kTH2F, {axisImpParErr, axisPt}});
			registry.add("hChi2PCA", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.5}, axisPt}});
			
			if (doprocess[2] || doprocess[3]) { // Xicp specific histograms -> Will be suported in the future!
				registry.add("hMass", "#Xi^{0}_{#plus} candidates; mass #Xi^{#mp} #pi^{#pm} #pi^{#pm}; entries", {HistType::kTH2F, {axisMassXic0Xicp, axisPt}});
				registry.add("hPtProng0", "#Xi^{0}_{c} candidates; #Xi^{#mp} #it{p}_{T}; entries", {HistType::kTH1F, {{200, 0., 20.}}});
				registry.add("hPtProng1", "#Xi^{0}_{c} candidates; #pi^{#pm} #it{p}_{T}; entries", {HistType::kTH1F, {{200, 0., 20.}}});
				registry.add("hPtProng2", "#Xi^{0}_{c} candidates; #pi^{#pm} #it{p}_{T}; entries", {HistType::kTH1F, {{200, 0., 20.}}});
				registry.add("hPtProng0VsPt", "#Xi^{0}_{c} candidates; #Xi^{#mp} #it{p}_{T}; entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
				registry.add("hPtProng1VsPt", "#Xi^{0}_{c} candidates; #pi^{#pm} #it{p}_{T}; entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
				registry.add("hPtProng2VsPt", "#Xi^{0}_{c} candidates; #pi^{#pm} #it{p}_{T}; entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});     
				registry.add("hd0Prong0", "#Xi^{#plus}_{c} candidates;prong 0 (#Xi^{#mp}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
				registry.add("hd0Prong1", "#Xi^{#plus}_{c} candidates;prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
				registry.add("hd0Prong2", "#Xi^{#plus}_{c} candidates;prong 2 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
				registry.add("hImpParErr", "#Xi^{#plus}_{c} candidates;prongs impact parameter error (cm);entries", {HistType::kTH2F, {axisImpParErr, axisPt}});
				registry.add("hChi2PCA", "#Xi^{#plus}_{c} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.5}, axisPt}});
			}
			
			if (doprocess[1] || doprocess[3]) { // If KF specific histograms
				registry.add("hChi2geoXi", "#Xi^{#plus}_{c} candidates;#Xi^{#mp} #chi^{2}_{geo};entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2geoLam", "#Xi^{#plus}_{c} candidates;#Lambda #chi^{2}_{geo};entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2topoToPV", "#Xi^{#plus}_{c} candidates;#Xi^{#plus}_{c} candidate #chi^{2}_{topo} to PV;entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2topoXiToXic0Xicp", "#Xi^{#plus}_{c} candidates;#Xi^{#mp} candidate #chi^{2}_{topo} to #Xi^{#plus}_{c};entries", {HistType::kTH2F, {axisChi2, axisPt}});
			}
		} // end Data

		// Define histograms for MC analysis
		if (doprocess[4] || doprocess[5] || doprocess[6] || doprocess[7]) {
			// MC reconstructed - Common features
			registry.add("hPtGenSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (gen+rec);candidate #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      		registry.add("hPtRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      		registry.add("hPtRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      		registry.add("hPtProng0RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 0 (#Xi^{#mp}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      		registry.add("hPtProng0RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 0 (#Xi^{#mp}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      		registry.add("hPtProng1RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 1 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
      		registry.add("hPtProng1RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 1 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
      		registry.add("hPtProng0VsPtRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{#mp} #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
      		registry.add("hPtProng0VsPtRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{#mp} #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
      		registry.add("hPtProng1VsPtRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 1 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
      		registry.add("hPtProng1VsPtRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 1 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
      		registry.add("hEtaRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      		registry.add("hEtaRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      		registry.add("hRapidityRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      		registry.add("hRapidityRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      		registry.add("hSVxRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate secondary vertex position x (cm);entries", {HistType::kTH2F, {axisSV, axisPt}});
      		registry.add("hSVxRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate secondary vertex position x (cm);entries", {HistType::kTH2F, {axisSV, axisPt}});
      		registry.add("hSVyRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate secondary vertex position y (cm);entries", {HistType::kTH2F, {axisSV, axisPt}});
      		registry.add("hSVyRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate secondary vertex position y (cm);entries", {HistType::kTH2F, {axisSV, axisPt}});
      		registry.add("hSVzRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate secondary vertex position z (cm);entries", {HistType::kTH2F, {axisSV, axisPt}});
      		registry.add("hSVzRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate secondary vertex position z (cm);entries", {HistType::kTH2F, {axisSV, axisPt}});
      		registry.add("hCPARecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      		registry.add("hCPARecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      		registry.add("hCPAxyRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate CPAxy;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      		registry.add("hCPAxyRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate CPAxy;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      		registry.add("hMassRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);inv. mass  #Xi^{#mp} #pi^{#pm} #pi^{#pm} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.8, 3.0}, axisPt}});
      		registry.add("hMassRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);inv. mass  #Xi^{#mp} #pi^{#pm} #pi^{#pm} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.8, 3.0}, axisPt}});
      		registry.add("hDecLengthRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate decay length (cm);entries", {HistType::kTH2F, {axisDecayLength, axisPt}});
      		registry.add("hDecLengthRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate decay length (cm);entries", {HistType::kTH2F, {axisDecayLength, axisPt}});
      		registry.add("hErrDecLengthRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate decay length (cm);entries", {HistType::kTH2F, {axisErrDecayLength, axisPt}});
      		registry.add("hErrDecLengthRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate decay length (cm);entries", {HistType::kTH2F, {axisErrDecayLength, axisPt}});
      		registry.add("hDecLengthXYRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate decay length xy (cm);entries", {HistType::kTH2F, {axisDecayLength, axisPt}});
      		registry.add("hDecLengthXYRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate decay length xy(cm);entries", {HistType::kTH2F, {axisDecayLength, axisPt}});
      		registry.add("hErrDecLengthXYRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate decay length xy (cm);entries", {HistType::kTH2F, {axisErrDecayLength, axisPt}});
      		registry.add("hErrDecLengthXYRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{0}_{c}, #Xi^{#plus}_{c} candidate decay length xy(cm);entries", {HistType::kTH2F, {axisErrDecayLength, axisPt}});
      		registry.add("hd0Prong0RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 0 (#Xi^{#mp}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
      		registry.add("hd0Prong0RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 0 (#Xi^{#mp}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
      		registry.add("hd0Prong1RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
		  	registry.add("hd0Prong1RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
		  	registry.add("hImpParErrRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prongs impact parameter error (cm);entries", {HistType::kTH2F, {axisImpParErr, axisPt}});
		  	registry.add("hImpParErrRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prongs impact parameter error (cm);entries", {HistType::kTH2F, {axisImpParErr, axisPt}});
		  	registry.add("hChi2PCARecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, axisPt}});
		  	registry.add("hChi2PCARecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, axisPt}});
		  	registry.add("hCPAXiRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{#minus} cosine of pointing angle;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
		  	registry.add("hCPAXiRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{#minus} cosine of pointing angle;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
		  	registry.add("hCPAxyXiRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Xi^{#minus} candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
		  	registry.add("hCPAxyXiRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Xi^{#minus} candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
		  	registry.add("hCPALambdaRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Lambda candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
		  	registry.add("hCPALambdaRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Lambda candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
		  	registry.add("hCPAxyLambdaRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);#Lambda candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
		  	registry.add("hCPAxyLambdaRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);#Lambda candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});

			// Xicp spedific histograms
			if(doprocess[6] || doprocess[7]) {
				// MC Reco
				registry.add("hPtProng2RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
      			registry.add("hPtProng2RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
				registry.add("hPtProng2vsPtRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
      			registry.add("hPtProng2vsPtRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
				registry.add("hd0Prong2RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 2 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
		  		registry.add("hd0Prong2RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 2 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
				registry.add("hMassXiPi1RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 1) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
				registry.add("hMassXiPi1RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 1) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
		  		registry.add("hMassXiPi2RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 2) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
		  		registry.add("hMassXiPi2RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 2) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});

				// MC Generated
				registry.add("hPtProng2Gen", "MC particles (generated);prong 2 (#pi^{#pm}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
				registry.add("hYProng2Gen", "MC particles (generated);prong 2 (#pi^{#pm}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
				registry.add("hYProng2Gen", "MC particles (generated);prong 2 (#pi^{#pm}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
			}
      
			// MC reconstructed with KFParticle
			if (doprocess[5] || doprocess[7]) {
				registry.add("hChi2topoToPVRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#plus}_{c} candidate #chi^{2}_{topo} to PV;entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2topoToPVRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#plus}_{c} candidate #chi^{2}_{topo} to PV;entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2geoXiRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#mp} #chi^{2}_{geo};entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2geoXiRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#mp} #chi^{2}_{geo};entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2geoLamRecSig", "#Xi^{#plus}_{c} candidates (matched);#Lambda #chi^{2}_{geo};entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2geoLamRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Lambda #chi^{2}_{geo};entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2topoXiToXic0XicpRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#mp} candidate #chi^{2}_{topo} to #Xi^{#plus}_{c};entries", {HistType::kTH2F, {axisChi2, axisPt}});
				registry.add("hChi2topoXiToXic0XicpRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#mp} candidate #chi^{2}_{topo} to #Xi^{#plus}_{c};entries", {HistType::kTH2F, {axisChi2, axisPt}});
			}
			// MC generated
			registry.add("hPtProng0Gen", "MC particles (generated);prong 0 (#Xi^{#mp}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{300, 0., 30.}, axisPt}});
			registry.add("hPtProng1Gen", "MC particles (generated);prong 1 (#pi^{#pm}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
			registry.add("hEtaProng0Gen", "MC particles (generated);prong 0 (#Xi^{#mp}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
			registry.add("hEtaProng1Gen", "MC particles (generated);prong 1 (#pi^{#pm}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
			registry.add("hYProng0Gen", "MC particles (generated);prong 0 (#Xi^{#mp}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
			registry.add("hYProng1Gen", "MC particles (generated);prong 1 (#pi^{#pm}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
			registry.add("hPtGen", "MC particles (generated);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
			registry.add("hEtaGen", "MC particles (generated);#Xi^{#plus}_{c} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
			registry.add("hYGen", "MC particles (generated);#Xi^{#plus}_{c} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
			registry.add("hSVxGen", "#Xi^{#plus}_{c} candidates (generated);#Xi^{#plus}_{c} candidate secondary vertex position x (cm);entries", {HistType::kTH2F, {axisSV, axisPt}});
			registry.add("hSVyGen", "#Xi^{#plus}_{c} candidates (generated);#Xi^{#plus}_{c} candidate secondary vertex position y (cm);entries", {HistType::kTH2F, {axisSV, axisPt}});
			registry.add("hSVzGen", "#Xi^{#plus}_{c} candidates (generated);#Xi^{#plus}_{c} candidate secondary vertex position z (cm);entries", {HistType::kTH2F, {axisSV, axisPt}});
			registry.add("hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
			registry.add("hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);#Xi^{#plus}_{c} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
			registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);#Xi^{#plus}_{c} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});

			if(doprocess[6] || doprocess[7]) { // Xicp specific columns
				// MC Reco
				registry.add("hPtProng2RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
      			registry.add("hPtProng2RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
				registry.add("hPtProng2vsPtRecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
      			registry.add("hPtProng2vsPtRecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
				registry.add("hd0Prong2RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);prong 2 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
		  		registry.add("hd0Prong2RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);prong 2 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisDCA, axisPt}});
				registry.add("hMassXiPi1RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 1) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
				registry.add("hMassXiPi1RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 1) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
		  		registry.add("hMassXiPi2RecSig", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (matched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 2) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
		  		registry.add("hMassXiPi2RecBg", "#Xi^{0}_{c}, #Xi^{#plus}_{c} candidates (unmatched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 2) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});

				// MC Generated
				registry.add("hPtProng2Gen", "MC particles (generated);prong 2 (#pi^{#pm}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
				registry.add("hYProng2Gen", "MC particles (generated);prong 2 (#pi^{#pm}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
				registry.add("hYProng2Gen", "MC particles (generated);prong 2 (#pi^{#pm}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
			}

		} // end MC
		  
		if (configs.enableTHn) { // -> Will be implemented in the future
		} // end THn

	}// end init

	/// Selection of Xic0(Xicp) daughter in geometrical acceptance
	/// \param etaProng is the pseudo rapidity of prong
	/// \param pTProng is the pT of prong
	/// \return true if prong is in geometrical acceptance
	template <typename T = float>
	bool isProngInAcceptance(const T& etaProng, const T& ptProng)
	{
		return std::abs(etaProng) <= configs.etaTrackMax && ptProng >= configs.ptTrackMin;
	}
	
	// Function to fill histogram
	template <bool doXicp, bool useKfParticle, bool useMl, typename TCandTable>
	void fillHistograms(TCandTable const& candidates)
	{
		for (const auto& candidate : candidates) {

			auto yCandXic0Xicp = (doXicp) ? candidate.y(o2::constants::physics::MassXiCPlus) : candidate.y(o2::constants::physics::MassXiC0);
			if (configs.yCandRecoMax >= 0. && std::abs(yCandXic0Xicp) > configs.yCandRecoMax) {
				continue;
			}
			
			auto ptCandXic0Xicp = candidate.pt();

			registry.fill(HIST("hPt"), ptCandXic0Xicp);
			registry.fill(HIST("hPtProng0"), candidate.ptProng0());
			registry.fill(HIST("hPtProng1"), candidate.ptProng1());
			registry.fill(HIST("hEta"), candidate.eta(), ptCandXic0Xicp);
			registry.fill(HIST("hRapidity"), yCandXic0Xicp, ptCandXic0Xicp);
			registry.fill(HIST("hCPA"), candidate.cpa(), ptCandXic0Xicp);
			registry.fill(HIST("hCPAXY"), candidate.cpaXY(), ptCandXic0Xicp);
			registry.fill(HIST("hMass"), candidate.invMassXic0(), ptCandXic0Xicp);
			registry.fill(HIST("hDecayLength"), candidate.decayLength(), ptCandXic0Xicp); 
			registry.fill(HIST("hErrDecayLength"), candidate.errorDecayLength(), ptCandXic0Xicp); 
			registry.fill(HIST("hErrDecayLengthXY"), candidate.errorDecayLengthXY(), ptCandXic0Xicp); 
			registry.fill(HIST("hSVx"), candidate.xSecondaryVertex(), ptCandXic0Xicp);
			registry.fill(HIST("hSVy"), candidate.ySecondaryVertex(), ptCandXic0Xicp);
			registry.fill(HIST("hSVz"), candidate.zSecondaryVertex(), ptCandXic0Xicp);
			registry.fill(HIST("hPtProng0VsPt"), candidate.ptProng0(), ptCandXic0Xicp);
			registry.fill(HIST("hPtProng1VsPt"), candidate.ptProng1(), ptCandXic0Xicp);
			registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), ptCandXic0Xicp);
			registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), ptCandXic0Xicp);
			registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), ptCandXic0Xicp);
			registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), ptCandXic0Xicp);
			registry.fill(HIST("hChi2PCA"), candidate.chi2PCA(), ptCandXic0Xicp);
			registry.fill(HIST("hCPAXi"), candidate.cosPaXi(), ptCandXic0Xicp);
			registry.fill(HIST("hCPAxyXi"), candidate.cosPaXYXi(), ptCandXic0Xicp);
			registry.fill(HIST("hCPALambda"), candidate.cosPaLambda(), ptCandXic0Xicp);
			registry.fill(HIST("hCPAxyLambda"), candidate.cosPaXYLambda(), ptCandXic0Xicp);

			if constexpr(!doXicp && useKfParticle) {
				registry.fill(HIST("hChi2topoToPV"), candidate.chi2TopoXic0ToPV(), ptCandXic0Xicp);
				registry.fill(HIST("hChi2topoXiToXic0Xicp"), candidate.chi2TopoXiToXic0(), ptCandXic0Xicp);
				registry.fill(HIST("hChi2geoXi"), candidate.kfCascadeChi2(), ptCandXic0Xicp);
				registry.fill(HIST("hChi2geoLam"), candidate.kfV0Chi2(), ptCandXic0Xicp);
			}

			if constexpr (doXicp) {
				// -> Will be implemented in the future. Currently, the derived table structure for Xicp not decided
			}
		}

	}

	// Function for MC analysis and histogram filling
	// Now, only analysis for xic0 is implemented 
	template <bool doXicp, bool useKfParticle, bool useMl, typename TCandTable>
	void fillHistogramsMc(TCandTable const& candidates,
						  soa::Join<aod::McParticles, aod::HfCandXic0McGen> const& mcParticles,
						  aod::TracksWMc const&)
	{
		std::vector<int> arrDaughIdx;

		// MC rec
		for (const auto& candidate : candidates) {
			auto yCandXic0 = candidate.y(o2::constants::physics::MassXiC0);
			if (configs.yCandRecoMax >= 0. && std::abs(yCandXic0) > configs.yCandRecoMax) {
				continue;
			}

			auto ptCandXic0 = candidate.pt();
			int flagMatchRecXic0 = std::abs(candidate.flagMcMatchRec());

			if (TESTBIT(flagMatchRecXic0, hf_cand_xic0_xicp_to_hadronic::DecayType::Xic0ToXiPi)) { // process signal
				auto idxMother = RecoDecay::getMother(mcParticles, candidate.template pi_as<aod::TracksWMc>().template mcParticle_as<soa::Join<aod::McParticles, aod::HfCandXic0McGen>>(), o2::constants::physics::Pdg::kXiC0, true);
				auto particleMother = mcParticles.rawIteratorAt(idxMother);
				        
				registry.fill(HIST("hPtGenSig"), particleMother.pt());
				registry.fill(HIST("hPtRecSig"), ptCandXic0);
				registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0());
				registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1());
				registry.fill(HIST("hPtProng0VsPtRecSig"), candidate.ptProng0(), ptCandXic0);
				registry.fill(HIST("hPtProng1VsPtRecSig"), candidate.ptProng1(), ptCandXic0);
				registry.fill(HIST("hEtaRecSig"), candidate.eta(), ptCandXic0);
				registry.fill(HIST("hRapidityRecSig"), yCandXic0, ptCandXic0);
				registry.fill(HIST("hSVxRecSig"), candidate.xSecondaryVertex(), ptCandXic0);
				registry.fill(HIST("hSVyRecSig"), candidate.ySecondaryVertex(), ptCandXic0);
				registry.fill(HIST("hSVzRecSig"), candidate.zSecondaryVertex(), ptCandXic0);
				registry.fill(HIST("hCPARecSig"), candidate.cpa(), ptCandXic0);
				registry.fill(HIST("hCPAxyRecSig"), candidate.cpaXY(), ptCandXic0);
				registry.fill(HIST("hMassRecSig"), candidate.invMassXic0(), ptCandXic0);
				registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), ptCandXic0);
				registry.fill(HIST("hErrDecLengthRecSig"), candidate.errorDecayLength(), ptCandXic0);
				registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), ptCandXic0);
				registry.fill(HIST("hErrDecLengthXYRecSig"), candidate.errorDecayLengthXY(), ptCandXic0);
				registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), ptCandXic0);
				registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), ptCandXic0);
				registry.fill(HIST("hImpParErrRecSig"), candidate.errorImpactParameter0(), ptCandXic0);
				registry.fill(HIST("hImpParErrRecSig"), candidate.errorImpactParameter1(), ptCandXic0);
				registry.fill(HIST("hChi2PCARecSig"), candidate.chi2PCA(), ptCandXic0);
				registry.fill(HIST("hCPAXiRecSig"), candidate.cosPaXi(), ptCandXic0);
				registry.fill(HIST("hCPAxyXiRecSig"), candidate.cosPaXYXi(), ptCandXic0);
				registry.fill(HIST("hCPALambdaRecSig"), candidate.cosPaLambda(), ptCandXic0);
				registry.fill(HIST("hCPAxyLambdaRecSig"), candidate.cosPaXYLambda(), ptCandXic0);

				if constexpr (useKfParticle) {
					registry.fill(HIST("hChi2topoToPVRecSig"), candidate.chi2TopoXic0ToPV(), ptCandXic0);
					registry.fill(HIST("hChi2topoXiToXic0XicpRecSig"), candidate.chi2TopoXiToXic0(), ptCandXic0);
					registry.fill(HIST("hChi2geoXiRecSig"), candidate.kfCascadeChi2(), ptCandXic0);
					registry.fill(HIST("hChi2geoLamRecSig"), candidate.kfV0Chi2(), ptCandXic0);
				}

			} else { // process background
				
				registry.fill(HIST("hPtRecBg"), ptCandXic0);
				registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0());
				registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1());
				registry.fill(HIST("hPtProng0VsPtRecBg"), candidate.ptProng0(), ptCandXic0);
				registry.fill(HIST("hPtProng1VsPtRecBg"), candidate.ptProng1(), ptCandXic0);
				registry.fill(HIST("hEtaRecBg"), candidate.eta(), ptCandXic0);
				registry.fill(HIST("hRapidityRecBg"), yCandXic0, ptCandXic0);
				registry.fill(HIST("hSVxRecBg"), candidate.xSecondaryVertex(), ptCandXic0);
				registry.fill(HIST("hSVyRecBg"), candidate.ySecondaryVertex(), ptCandXic0);
				registry.fill(HIST("hSVzRecBg"), candidate.zSecondaryVertex(), ptCandXic0);
				registry.fill(HIST("hCPARecBg"), candidate.cpa(), ptCandXic0);
				registry.fill(HIST("hCPAxyRecBg"), candidate.cpaXY(), ptCandXic0);
				registry.fill(HIST("hMassRecBg"), candidate.invMassXic0(), ptCandXic0);
				registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), ptCandXic0);
				registry.fill(HIST("hErrDecLengthRecBg"), candidate.errorDecayLength(), ptCandXic0);
				registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), ptCandXic0);
				registry.fill(HIST("hErrDecLengthXYRecBg"), candidate.errorDecayLengthXY(), ptCandXic0);
				registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), ptCandXic0);
				registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), ptCandXic0);
				registry.fill(HIST("hImpParErrRecBg"), candidate.errorImpactParameter0(), ptCandXic0);
				registry.fill(HIST("hImpParErrRecBg"), candidate.errorImpactParameter1(), ptCandXic0);
				registry.fill(HIST("hChi2PCARecBg"), candidate.chi2PCA(), ptCandXic0);
				registry.fill(HIST("hCPAXiRecBg"), candidate.cosPaXi(), ptCandXic0);
				registry.fill(HIST("hCPAxyXiRecBg"), candidate.cosPaXYXi(), ptCandXic0);
				registry.fill(HIST("hCPALambdaRecBg"), candidate.cosPaLambda(), ptCandXic0);
				registry.fill(HIST("hCPAxyLambdaRecBg"), candidate.cosPaXYLambda(), ptCandXic0);
				
				if constexpr (useKfParticle) {
					registry.fill(HIST("hChi2topoToPVRecBg"), candidate.chi2TopoXic0ToPV(), ptCandXic0);
					registry.fill(HIST("hChi2topoXiToXic0XicpRecBg"), candidate.chi2TopoXiToXic0(), ptCandXic0);
					registry.fill(HIST("hChi2geoXiRecBg"), candidate.kfCascadeChi2(), ptCandXic0);
					registry.fill(HIST("hChi2geoLamRecBg"), candidate.kfV0Chi2(), ptCandXic0);
				}
			}
		} // end MC rec

		// MC gen
		for (const auto& particle : mcParticles) {

		if (TESTBIT(std::abs(particle.flagMcMatchGen()), hf_cand_xic0_xicp_to_hadronic::DecayType::Xic0ToXiPi)) {
			arrDaughIdx.clear();

			auto ptParticle = particle.pt();
			auto yParticle = RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiC0);
			if (configs.yCandGenMax >= 0. && std::abs(yParticle) > configs.yCandGenMax) {
				continue;
			}

			// Get Kinematic variables of Xi- Pi+
			std::array<float, 3> ptProngs;
			std::array<float, 3> yProngs;
			std::array<float, 3> etaProngs;
			std::array<float, 3> prodVtxXProngs;
			std::array<float, 3> prodVtxYProngs;
			std::array<float, 3> prodVtxZProngs;
			int counter = 0;
			RecoDecay::getDaughters(particle, &arrDaughIdx, std::array{+kXiMinus, +kPiPlus}, 1); // -> Max depth of search...1 or 2?

			for (auto iProng = 0u; iProng < arrDaughIdx.size(); ++iProng) {
				auto daughI = mcParticles.rawIteratorAt(arrDaughIdx[iProng]);
				ptProngs[counter] = daughI.pt();
				etaProngs[counter] = daughI.pt();
				yProngs[counter] = daughI.pt();
				prodVtxXProngs[counter] = daughI.vx();
				prodVtxYProngs[counter] = daughI.vy();
				prodVtxZProngs[counter] = daughI.vz();
				counter++;
			}

			registry.fill(HIST("hPtProng0Gen"), ptProngs[0], ptParticle);
			registry.fill(HIST("hPtProng1Gen"), ptProngs[1], ptParticle);
			registry.fill(HIST("hEtaProng0Gen"), etaProngs[0], ptParticle);
			registry.fill(HIST("hEtaProng1Gen"), etaProngs[1], ptParticle);
			registry.fill(HIST("hYProng0Gen"), yProngs[0], ptParticle);
			registry.fill(HIST("hYProng1Gen"), yProngs[1], ptParticle);
			registry.fill(HIST("hPtGen"), ptParticle);
			registry.fill(HIST("hYGen"), yParticle, ptParticle);
			registry.fill(HIST("hEtaGen"), particle.eta(), ptParticle);
			registry.fill(HIST("hSVxGen"), prodVtxXProngs[0], ptParticle);
			registry.fill(HIST("hSVyGen"), prodVtxYProngs[0], ptParticle);
			registry.fill(HIST("hSVzGen"), prodVtxZProngs[0], ptParticle);

			if (!isProngInAcceptance(etaProngs[0], ptProngs[0]) || isProngInAcceptance(etaProngs[1], ptProngs[1]) || isProngInAcceptance(etaProngs[2], ptProngs[2])) {
				continue;
			}
			registry.fill(HIST("hPtGenWithProngsInAcceptance"), ptParticle);
			registry.fill(HIST("hEtaGenWithProngsInAcceptance"), particle.eta(), ptParticle);
			registry.fill(HIST("hYGenWithProngsInAcceptance"), yParticle, ptParticle);
		}

		}// end MC gen
	}

	////////////////////////////////////////////
	//										  //
	//    Data analysis and fill histograms   // 
	//										  //
	////////////////////////////////////////////
	
	Filter filterSelectCandidates = (aod::hf_sel_xic0_xicp_to_hadronic::isSelXic0ToXiPi >= configs.selectionFlagXic0Xicp);

	void processXic0WithDCAFitter(soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfSelXic0ToXiPi>> const& candidates)
	{
		fillHistograms<false, false, false>(candidates);
	}
	PROCESS_SWITCH(HfTaskXic0XicpToHadronic, processXic0WithDCAFitter, "Process data(Xic0) with DCAFitter", true);

	void processXic0WithKFParticle(soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfCandXic0KF, aod::HfSelXic0ToXiPi>> const& candidates)
	{
		fillHistograms<false, true, false>(candidates);
	}
	PROCESS_SWITCH(HfTaskXic0XicpToHadronic, processXic0WithKFParticle, "Process data(Xic0) with DCAFitter", false);

	void processXicpWithDCAFitter(aod::Collisions const&)
	{

	}
	PROCESS_SWITCH(HfTaskXic0XicpToHadronic, processXicpWithDCAFitter, "Process data with(Xicp) DCAFitter", false);

	void processXicpWithKFParticle(aod::Collisions const&)
	{

	}
	PROCESS_SWITCH(HfTaskXic0XicpToHadronic, processXicpWithKFParticle, "Process data with(Xicp) DCAFitter", false);

	//////////////////////////////////////////
	//										//
	//    MC analysis and fill histograms   // 
	//										//
	//////////////////////////////////////////
	
	void processMcXic0WithDCAFitter(soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfSelXic0ToXiPi, aod::HfCandXic0McRec>> const& candidates,
									soa::Join<aod::McParticles, aod::HfCandXic0McGen> const& mcParticles,
									aod::TracksWMc const& tracksWMc)
	{
		fillHistogramsMc<false, false, false>(candidates, mcParticles, tracksWMc);
	}
	PROCESS_SWITCH(HfTaskXic0XicpToHadronic, processMcXic0WithDCAFitter, "Process mc(Xic0) with DCAFitter", false);

	void processMcXic0WithKFParticle(soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfCandXic0KF, aod::HfSelXic0ToXiPi, aod::HfCandXic0McRec>> const& candidates,
									 soa::Join<aod::McParticles, aod::HfCandXic0McGen> const& mcParticles,
									 aod::TracksWMc const& tracksWMc)
	{
		fillHistogramsMc<false, true, false>(candidates, mcParticles, tracksWMc);

	}
	PROCESS_SWITCH(HfTaskXic0XicpToHadronic, processMcXic0WithKFParticle, "Process mc(Xic0) with DCAFitter", false);

	void processMcXicpWithDCAFitter(aod::Collisions const&)
	{

	}
	PROCESS_SWITCH(HfTaskXic0XicpToHadronic, processMcXicpWithDCAFitter, "Process mc(Xicp) with DCAFitter", false);

	void processMcXicpWithKFParticle(aod::Collisions const&)
	{

	}
	PROCESS_SWITCH(HfTaskXic0XicpToHadronic, processMcXicpWithKFParticle, "Process mc(Xicp) with DCAFitter", false);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return WorkflowSpec
	{
		adaptAnalysisTask<HfTaskXic0XicpToHadronic>(cfgc)
	};
}
