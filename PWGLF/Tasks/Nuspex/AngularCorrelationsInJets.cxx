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
//
// Authors: Lars Jörgensen
// Date: 29.05.2024

#include <cmath>
#include <vector>
#include <TMath.h>
#include <TObjArray.h>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "TVector2.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA,aod::pidTPCLfFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>;

struct AngularCorrelationsInJets {

	// HistogramRegistry registryName{"folderTitle", {}, OutputObjHandlingPolicy::AnalysisObject, <sortHistograms:bool>, <createDir:bool>};
	HistogramRegistry registryData{"jetOutput", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
	// OutputObj<Type> histogramName{Type("histogramName", "histogramTitle;Axis", nbins,minbin,maxbin)};

	void init(InitContext&) {

		// std::vector<double> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};

		// AxisSpec specName = {binningInfo, "axisLabel"};
		AxisSpec ptAxis 		  = {1000,0,100, "#it{p}_{T} [GeV/#it{c}]"};
		AxisSpec particleTypeAxis = {4,1,5, "[p, ap, d, ad]"};
		AxisSpec nsigmapTAxis 	  = {1000,-50,50, "#it{p}_{T} [GeV/#it{c}]"};
		AxisSpec nsigmaAxis 	  = {1000,-5,5, "n#sigma"};
		AxisSpec dcazAxis		  = {200,-3,3, "DCA_{z} [cm]"};
		AxisSpec dcaxyAxis		  = {200,-2,2, "DCA_{xy} [cm]"};
		AxisSpec angDistPhiAxis	  = {1000,-2,5, "#Delta#varphi"};
		AxisSpec angDistEtaAxis   = {1000,-2,2, "#Delta#eta"};

		// registryName.add("histogramName", "histogramTitle", HistType::Type, {{binningInfo}});

		// Counters
		registryData.add("hNumberOfEvents", "Number of events", HistType::kTH1I, {{1,0,1}});
		registryData.add("hNumberOfJets", "Total number of jets", HistType::kTH1I, {{1,0,1}});
		registryData.add("hEventProtocol", "Event protocol", HistType::kTH1I, {{20,0,20}});
		registryData.add("hTrackProtocol", "Track protocol", HistType::kTH1I, {{20,0,20}});
		registryData.add("hNumPartInJet", "Number of particles in a jet", HistType::kTH1I, {{200,0,200}});

		// (Pseudo)Rapidity
		registryData.add("hEtaFullEvent", "Particle pseudorapidity;#eta", HistType::kTH1F, {{200,-1,1}});
		registryData.add("hJetRapidity", "Jet rapidity;#it{y}", HistType::kTH1F, {{200,-1,1}});

		// pT
		registryData.add("hPtFullEvent", "p_{T} after basic cuts", HistType::kTH1F, {ptAxis});
		registryData.add("hPtJetParticle", "p_{T} of particles in jets", HistType::kTH1F, {ptAxis});
		registryData.add("hPtSubtractedJet", "Subtracted jet p_{T}", HistType::kTH1F, {ptAxis});
		registryData.add("hPtJetProtonDeuteron", "p_{T} of (anti)p, (anti)d", HistType::kTH2F, {particleTypeAxis, ptAxis});
		registryData.add("hPtTotalJet", "p_{T} of entire jet;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1F, {{2000,0,500}});
		registryData.add("hPtDiff", "pT difference PseudoJet/original track;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1D, {{100,-0.0000005,0.0000005}});

		// nSigma
		registryData.add("hTPCnsigmaProton", "TPC n#sigma for (anti)proton", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
		registryData.add("hTOFnsigmaProton", "TOF n#sigma for (anti)proton", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
		// registryData.add("hITSnsigmaProton", "ITS n#sigma for (anti)proton", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
		registryData.add("hTPCnsigmaDeuteron", "TPC n#sigma for (anti)deuteron", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
		registryData.add("hTOFnsigmaDeuteron", "TOF n#sigma for (anti)deuteron", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
		// registryData.add("hITSnsigmaDeuteron", "ITS n#sigma for (anti)deuteron", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});

		// DCA
		registryData.add("hDCAxyFullEvent", "DCA_{xy} of full event", HistType::kTH1F, {dcaxyAxis});
		registryData.add("hDCAzFullEvent", "DCA_{z} of full event", HistType::kTH1F, {dcazAxis});
		registryData.add("hDCAzJetParticle", "DCA_{z} of jet particles", HistType::kTH1F, {dcazAxis});

		// Angular Distributions
		registryData.add("hDeltaPhiSE", "#Delta#varphi of (anti)p, (anti)d in single event", HistType::kTH2D, {particleTypeAxis, angDistPhiAxis});
		registryData.add("hDeltaPhiME", "#Delta#varphi of (anti)p, (anti)d in mixed events", HistType::kTH2D, {particleTypeAxis, angDistPhiAxis});
		registryData.add("hDeltaPhiEtaSE", "#Delta#varphi vs #Delta#eta of (anti)p, (anti)d in single event", HistType::kTH3D, {particleTypeAxis, angDistPhiAxis, angDistEtaAxis});
		registryData.add("hDeltaPhiEtaME", "#Delta#varphi vs #Delta#eta of (anti)p, (anti)d in mixed events", HistType::kTH3D, {particleTypeAxis, angDistPhiAxis, angDistEtaAxis});
	}

	// Configurable<Type> cfgName{"nameOnHyperloopAndJson", value, "Flag shown on hyperloop"};
	// Preliminary Cuts
	// Configurable<int> fTPCRefit{"TPCRefit", 0, "Require TPC refit"};
	// Configurable<int> fITSRefit{"ITSRefit", 0, "Require ITS refit"};
	Configurable<float> fMinNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
	Configurable<float> fMinReqClusterITS{"minReqClusterITS", 2.0, "min number of clusters required in ITS"};
	Configurable<float> fMinRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.7f, "min ratio of crossed rows over findable clusters TPC"};
	Configurable<float> fMaxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
	Configurable<float> fMaxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
	Configurable<float> fMaxDCAxy{"maxDCA_xy", 0.5f, "max DCA to vertex xy"};
	Configurable<float> fMaxDCAz{"maxDCA_z", 2.4f, "max DCA to vertex z"};
	Configurable<float> fMaxEta{"maxEta", 0.8, "max pseudorapidity"};

	// Jet Cuts
	Configurable<float> fJetRadius{"jetRadius", 0.4, "jet radius R"};
	Configurable<float> fMinJetPt{"minJetPt", 10.0, "minimum total pT to accept jet"};
	Configurable<float> fMinJetParticlePt{"minJetParticlePt", 0.0, "minimum pT to accept jet particle"};
	Configurable<float> fMinLeadingPt{"minLeadingPt", 5.0, "minimum pT for leading track"};

	// Proton Cuts
	Configurable<float> fProtonDCAxy{"protonDCAxy", 0.5, "[proton] DCAxy cut"};
	Configurable<float> fProtonDCAz{"protonDCAz", 1.0, "[proton] DCAz cut"};
	Configurable<float> fProtonITSTOFpT{"protonITSTOFswitchpT", 0.7, "[proton] pT for switch between ITS+TPC/TPC+TOF"};
	Configurable<float> fProtonITSnsig{"protonITSnsigma", 3.0, "[proton] max ITS nsigma"};
	Configurable<float> fProtonTPCnsigITS{"protonTPCnsigmaIfITS", 5.0, "[proton] max TPC nsigma with ITS"};
	Configurable<float> fProtonTPCnsigTOF{"protonTPCnsigmaIfTOF", 4.0, "[proton] max TPC nsigma with TOF"};
	Configurable<float> fProtonTOFnsig{"protonTOFnsigma", 10.0, "[proton] max TOF nsigma"};

	// Antiproton Cuts
	Configurable<float> fAntiprotonDCAxy{"antiprotonDCAxy", 0.5, "[antiproton] DCAxy cut"};
	Configurable<float> fAntiprotonDCAz{"antiprotonDCAz", 1.0, "[antiproton] DCAz cut"};
	Configurable<float> fAntiprotonITSTOFpT{"antiprotonITSTOFswitchpT", 0.7, "[antiproton] pT for switch between ITS+TPC/TPC+TOF"};
	Configurable<float> fAntiprotonITSnsig{"antiprotonITSnsigma", 3.0, "[antiproton] max ITS nsigma"};
	Configurable<float> fAntiprotonTPCnsigITS{"antiprotonTPCnsigmaIfITS", 5.0, "[antiproton] max TPC nsigma with ITS"};
	Configurable<float> fAntiprotonTPCnsigTOF{"antiprotonTPCnsigmaIfTOF", 4.0, "[antiproton] max TPC nsigma with TOF"};
	Configurable<float> fAntiprotonTOFnsig{"antiprotonTOFnsigma", 10.0, "[antiproton] max TOF nsigma"};

	// Deuteron Cuts
	Configurable<float> fDeuteronDCAxy{"deuteronDCAxy", 0.5, "[deuteron] DCAxy cut"};
	Configurable<float> fDeuteronDCAz{"deuteronDCAz", 1.0, "[deuteron] DCAz cut"};
	Configurable<float> fDeuteronITSTOFpT{"deuteronITSTOFswitchpT", 0.7, "[deuteron] pT for switch between ITS+TPC/TPC+TOF"};
	Configurable<float> fDeuteronITSnsig{"deuteronITSnsigma", 3.0, "[deuteron] max ITS nsigma"};
	Configurable<float> fDeuteronTPCnsigITS{"deuteronTPCnsigmaIfITS", 5.0, "[deuteron] max TPC nsigma with ITS"};
	Configurable<float> fDeuteronTPCnsigTOF{"deuteronTPCnsigmaIfTOF", 4.0, "[deuteron] max TPC nsigma with TOF"};
	Configurable<float> fDeuteronTOFnsig{"deuteronTOFnsigma", 10.0, "[deuteron] max TOF nsigma"};

	// Antideuteron Cuts
	Configurable<float> fAntideuteronDCAxy{"antideuteronDCAxy", 0.5, "[antideuteron] DCAxy cut"};
	Configurable<float> fAntideuteronDCAz{"antideuteronDCAz", 1.0, "[antideuteron] DCAz cut"};
	Configurable<float> fAntideuteronITSTOFpT{"antideuteronITSTOFswitchpT", 0.7, "[antideuteron] pT for switch between ITS+TPC/TPC+TOF"};
	Configurable<float> fAntideuteronITSnsig{"antideuteronITSnsigma", 3.0, "[antideuteron] max ITS nsigma"};
	Configurable<float> fAntideuteronTPCnsigITS{"antideuteronTPCnsigmaIfITS", 5.0, "[antideuteron] max TPC nsigma with ITS"};
	Configurable<float> fAntideuteronTPCnsigTOF{"antideuteronTPCnsigmaIfTOF", 4.0, "[antideuteron] max TPC nsigma with TOF"};
	Configurable<float> fAntideuteronTOFnsig{"antideuteronTOFnsigma", 10.0, "[antideuteron] max TOF nsigma"};

	Configurable<int> fTrackBufferSize{"trackBufferSize", 2000, "Number of mixed-event tracks being stored"};

	std::vector<typename FullTracks::iterator> fTrackBufferProton;
	std::vector<typename FullTracks::iterator> fTrackBufferAntiproton;
	std::vector<typename FullTracks::iterator> fTrackBufferDeuteron;
	std::vector<typename FullTracks::iterator> fTrackBufferAntideuteron;

	//****************************************************************************************************

	template <typename CollisionType, typename TracksType>
	void jetReco(const CollisionType& collision, const TracksType& tracks) { // this turns TracksType into a pointer, no? How then can we access attributes with . instead of ->?
		registryData.fill(HIST("hNumberOfEvents"), 0);

		std::vector<fastjet::PseudoJet> jetInput;
		std::vector<fastjet::PseudoJet> jets;
		std::vector<fastjet::PseudoJet> constituents;
		fastjet::PseudoJet hardestJet(0.,0.,0.,0.);
		fastjet::PseudoJet subtractedJet(0.,0.,0.,0.);
		jetInput.clear();
		jets.clear();
		constituents.clear();

		for (const auto& track : tracks) {
			registryData.fill(HIST("hTrackProtocol"), 0);
			registryData.fill(HIST("hPtFullEvent"), track.pt());

			fastjet::PseudoJet inputPseudoJet(track.px(), track.py(), track.pz(), track.energy(track.mass())); // is this fine? just TOF for invariant mass?
			inputPseudoJet.set_user_index(track.globalIndex());
			jetInput.emplace_back(inputPseudoJet);
		}

		double ghost_maxrap = 1.0;
		double ghost_area = 0.005;
		int ghost_repeat = 1;
		fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, fJetRadius);
		fastjet::JetDefinition jetDefBkg(fastjet::kt_algorithm, 0.5);
		fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(ghost_maxrap, ghost_repeat, ghost_area));
		fastjet::AreaDefinition areaDefBkg(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghost_maxrap));
		fastjet::ClusterSequenceArea clusterSeq(jetInput, jetDef, areaDef); // or CSActiveArea?
		jets = sorted_by_pt(clusterSeq.inclusive_jets());
		if (jets.size() == 0) return;
		
		hardestJet = jets[0];

		if(hardestJet.pt() < fMinJetPt) return;
		if(hardestJet.constituents().size() < 2) return;
		registryData.fill(HIST("hNumberOfJets"), 0);

		constituents = hardestJet.constituents();

		for (int i=0; i<(int)constituents.size(); i++) {
			registryData.fill(HIST("hPtJetParticle"), constituents[i].pt());
		}

		fastjet::Selector selector = fastjet::SelectorAbsEtaMax(1.0) * (!fastjet::SelectorNHardest(2));
		fastjet::JetMedianBackgroundEstimator bkgEst(selector, jetDefBkg, areaDefBkg);
		fastjet::Subtractor subtractor(&bkgEst);
		subtractor.set_use_rho_m(true);
		bkgEst.set_particles(jetInput);

		subtractedJet = subtractor(hardestJet);
		if (subtractedJet.has_constituents()) {
			for (int i=0; i<(int)subtractedJet.constituents().size(); i++) {
				registryData.fill(HIST("hPtSubtractedJet"), subtractedJet.constituents()[i].pt());
			}
		}

		std::vector<typename TracksType::iterator> jetProtons;
		std::vector<typename TracksType::iterator> jetAntiprotons;
		std::vector<typename TracksType::iterator> jetDeuterons;
		std::vector<typename TracksType::iterator> jetAntideuterons;

		for (int i=0; i<(int)constituents.size(); i++) {
			fastjet::PseudoJet pseudoParticle = constituents[i];
        	int id = pseudoParticle.user_index();
			typename TracksType::iterator jetParticle = tracks.iteratorAt(id);
			double ptDiff = pseudoParticle.pt() - jetParticle.pt();
			registryData.fill(HIST("hPtDiff"), ptDiff);
			int particleType = 0;

			registryData.fill(HIST("hDCAzJetParticle"), jetParticle.dcaZ());
			if (jetParticle.pt()<fMinJetParticlePt) continue;
			if (IsProton(jetParticle) || IsAntiproton(jetParticle)) {// collect (anti)protons in jet
				registryData.fill(HIST("hTPCnsigmaProton"), jetParticle.pt()*jetParticle.sign(), jetParticle.tpcNSigmaPr());
				registryData.fill(HIST("hTOFnsigmaProton"), jetParticle.pt()*jetParticle.sign(), jetParticle.tofNSigmaPr());
				// registryData.fill(HIST("hITSnsigmaProton"),track.pt()*track.sign(), track.itsNSigmaPr());
				if (IsProton(jetParticle)) {
					particleType = 1;
					registryData.fill(HIST("hTrackProtocol"), 6.5); // # protons
					registryData.fill(HIST("hPtJetProtonDeuteron"), particleType, jetParticle.pt());
					jetProtons.emplace_back(jetParticle);
				} else {
					particleType = 2;
					registryData.fill(HIST("hTrackProtocol"), 7.5); // # antiprotons
					registryData.fill(HIST("hPtJetProtonDeuteron"), particleType, jetParticle.pt());
					jetAntiprotons.emplace_back(jetParticle);
				}
			} else if (IsDeuteron(jetParticle) || IsAntideuteron(jetParticle)) {// collect (anti)deuterons in jet
				registryData.fill(HIST("hTPCnsigmaDeuteron"), jetParticle.pt()*jetParticle.sign(), jetParticle.tpcNSigmaDe());
				registryData.fill(HIST("hTOFnsigmaDeuteron"), jetParticle.pt()*jetParticle.sign(), jetParticle.tofNSigmaDe());
				// registryData.fill(HIST("hITSnsigmaDeuteron"),track.pt()*track.sign(), track.itsNSigmaDe());
				if (IsDeuteron(jetParticle)) {
					particleType = 3;
					registryData.fill(HIST("hTrackProtocol"), 8.5); // # deuterons
					registryData.fill(HIST("hPtJetProtonDeuteron"), particleType, jetParticle.pt());
					jetDeuterons.emplace_back(jetParticle);
				} else {
					particleType = 4;
					registryData.fill(HIST("hTrackProtocol"), 9.5); // # antideuterons
					registryData.fill(HIST("hPtJetProtonDeuteron"), particleType, jetParticle.pt());
					jetAntideuterons.emplace_back(jetParticle);
				}
			}
		} // for (int i=0; i<(int)constituents.size(); i++)

		if ((jetProtons.size()<2) && (jetAntiprotons.size()<2) && (jetDeuterons.size()<2) && (jetAntideuterons.size()<2)) return;
		registryData.fill(HIST("hEventProtocol"), 5.5);

		if (jetProtons.size()>1) {
			for (int i=0; i<(int)jetProtons.size(); i++) {
				for (int j=i+1; j<(int)jetProtons.size(); j++) {
					double DeltaPhi = TVector2::Phi_0_2pi(jetProtons[i].phi() - jetProtons[j].phi());
					double DeltaEta = TMath::Abs(jetProtons[i].eta() - jetProtons[j].eta());
					if (DeltaPhi>(1.5*TMath::Pi())) {
						DeltaPhi = DeltaPhi - 2*TMath::Pi();
					}
					registryData.fill(HIST("hDeltaPhiSE"), 1, DeltaPhi);
					registryData.fill(HIST("hDeltaPhiEtaSE"), 1, DeltaPhi, DeltaEta);
				}
				FillMixedEventDeltas(jetProtons[i], 1);
				SetTrackBuffer(jetProtons[i], 1);
			}
		}
		if (jetAntiprotons.size()>1) {
			for (int i=0; i<(int)jetAntiprotons.size(); i++) {
				for (int j=i+1; j<(int)jetAntiprotons.size(); j++) {
					double DeltaPhi = TVector2::Phi_0_2pi(jetAntiprotons[i].phi() - jetAntiprotons[j].phi());
					double DeltaEta = TMath::Abs(jetAntiprotons[i].eta() - jetAntiprotons[j].eta());
					if (DeltaPhi>(1.5*TMath::Pi())) {
						DeltaPhi = DeltaPhi - 2*TMath::Pi();
					}
					registryData.fill(HIST("hDeltaPhiSE"), 2, DeltaPhi);
					registryData.fill(HIST("hDeltaPhiEtaSE"), 2, DeltaPhi, DeltaEta);
				}
				FillMixedEventDeltas(jetAntiprotons[i], 2);
				SetTrackBuffer(jetAntiprotons[i], 2);
			}
		}
		if (jetDeuterons.size()>1) {
			for (int i=0; i<(int)jetDeuterons.size(); i++) {
				for (int j=i+1; j<(int)jetDeuterons.size(); j++) {
					double DeltaPhi = TVector2::Phi_0_2pi(jetDeuterons[i].phi() - jetDeuterons[j].phi());
					double DeltaEta = TMath::Abs(jetDeuterons[i].eta() - jetDeuterons[j].eta());
					if (DeltaPhi>(1.5*TMath::Pi())) {
						DeltaPhi = DeltaPhi - 2*TMath::Pi();
					}
					registryData.fill(HIST("hDeltaPhiSE"), 3, DeltaPhi);
					registryData.fill(HIST("hDeltaPhiEtaSE"), 3, DeltaPhi, DeltaEta);
				}
				FillMixedEventDeltas(jetDeuterons[i], 3);
				SetTrackBuffer(jetDeuterons[i], 3);
			}
		}
		if (jetAntideuterons.size()>1) {
			for (int i=0; i<(int)jetAntideuterons.size(); i++) {
				for (int j=i+1; j<(int)jetAntideuterons.size(); j++) {
					double DeltaPhi = TVector2::Phi_0_2pi(jetAntideuterons[i].phi() - jetAntideuterons[j].phi());
					double DeltaEta = TMath::Abs(jetAntideuterons[i].eta() - jetAntideuterons[j].eta());
					if (DeltaPhi>(1.5*TMath::Pi())) {
						DeltaPhi = DeltaPhi - 2*TMath::Pi();
					}
					registryData.fill(HIST("hDeltaPhiSE"), 4, DeltaPhi);
					registryData.fill(HIST("hDeltaPhiEtaSE"), 4, DeltaPhi, DeltaEta);
				}
				FillMixedEventDeltas(jetAntideuterons[i], 4);
				SetTrackBuffer(jetAntideuterons[i], 4);
			}
		}
	}

	template <typename T1>
	bool IsProton(const T1& track) {
		bool isProton = false;
		if (track.sign()<0) return isProton;

		if (TMath::Abs(track.dcaXY()) > fProtonDCAxy) return isProton;
		if (TMath::Abs(track.dcaZ()) > fProtonDCAz)   return isProton;

		if (track.pt()<fProtonITSTOFpT/*  && TMath::Abs(track.itsNSigmaPr())<fProtonITSnsig */ && TMath::Abs(track.tpcNSigmaPr())<fProtonTPCnsigITS) isProton = true;
		if (track.pt()>fProtonITSTOFpT && TMath::Abs(track.tpcNSigmaPr())<fProtonTPCnsigTOF && TMath::Abs(track.tofNSigmaPr())<fProtonTOFnsig) isProton = true;

		return isProton;
	}

	template <typename T2>
	bool IsAntiproton(const T2& track) {
		bool isAntiproton = false;
		if (track.sign()<0) return isAntiproton;

		if (TMath::Abs(track.dcaXY()) > fAntiprotonDCAxy) return isAntiproton;
		if (TMath::Abs(track.dcaZ()) > fAntiprotonDCAz)   return isAntiproton;

		if (track.pt()<fAntiprotonITSTOFpT/*  && TMath::Abs(track.itsNSigmaPr())<fAntiprotonITSnsig */ && TMath::Abs(track.tpcNSigmaPr())<fAntiprotonTPCnsigITS) isAntiproton = true;
		if (track.pt()>fAntiprotonITSTOFpT && TMath::Abs(track.tpcNSigmaPr())<fAntiprotonTPCnsigTOF && TMath::Abs(track.tofNSigmaPr())<fAntiprotonTOFnsig) isAntiproton = true;

		return isAntiproton;
	}

	template <typename T3>
	bool IsDeuteron(const T3& track) {
		bool isDeuteron = false;
		if (track.sign()<0) return isDeuteron;

		if (TMath::Abs(track.dcaXY()) > fDeuteronDCAxy) return isDeuteron;
		if (TMath::Abs(track.dcaZ()) > fDeuteronDCAz)   return isDeuteron;

		if (track.pt()<fDeuteronITSTOFpT/*  && TMath::Abs(track.itsNSigmaPr())<fDeuteronITSnsig */ && TMath::Abs(track.tpcNSigmaDe())<fDeuteronTPCnsigITS) isDeuteron = true;
		if (track.pt()>fDeuteronITSTOFpT && TMath::Abs(track.tpcNSigmaPr())<fDeuteronTPCnsigTOF && TMath::Abs(track.tofNSigmaDe())<fDeuteronTOFnsig) isDeuteron = true;

		return isDeuteron;
	}

	template <typename T4>
	bool IsAntideuteron(const T4& track) {
		bool isAntideuteron = false;
		if (track.sign()<0) return isAntideuteron;

		if (TMath::Abs(track.dcaXY()) > fAntideuteronDCAxy) return isAntideuteron;
		if (TMath::Abs(track.dcaZ()) > fAntideuteronDCAz)   return isAntideuteron;

		if (track.pt()<fAntideuteronITSTOFpT/*  && TMath::Abs(track.itsNSigmaDe())<fAntideuteronITSnsig */ && TMath::Abs(track.tpcNSigmaDe())<fAntideuteronTPCnsigITS) isAntideuteron = true;
		if (track.pt()>fAntideuteronITSTOFpT && TMath::Abs(track.tpcNSigmaDe())<fAntideuteronTPCnsigTOF && TMath::Abs(track.tofNSigmaDe())<fAntideuteronTOFnsig) isAntideuteron = true;

		return isAntideuteron;
	}

	template <typename T5>
	void SetTrackBuffer(const T5& track, int particleType) {
		switch (particleType) {
			case 1:
				if ((int)fTrackBufferProton.size()==fTrackBufferSize) {
					fTrackBufferProton.insert(fTrackBufferProton.begin(), track);
					fTrackBufferProton.resize(fTrackBufferSize);
				} else if ((int)fTrackBufferProton.size()<fTrackBufferSize) {
					fTrackBufferProton.emplace_back(track);
				}
			break;
			case 2:
				if ((int)fTrackBufferAntiproton.size()==fTrackBufferSize) {
					fTrackBufferAntiproton.insert(fTrackBufferAntiproton.begin(), track);
					fTrackBufferAntiproton.resize(fTrackBufferSize);
				} else if ((int)fTrackBufferAntiproton.size()<fTrackBufferSize) {
					fTrackBufferAntiproton.emplace_back(track);
				}
			break;
			case 3:
				if ((int)fTrackBufferDeuteron.size()==fTrackBufferSize) {
					fTrackBufferDeuteron.insert(fTrackBufferDeuteron.begin(), track);
					fTrackBufferDeuteron.resize(fTrackBufferSize);
				} else if ((int)fTrackBufferDeuteron.size()<fTrackBufferSize) {
					fTrackBufferDeuteron.emplace_back(track);
				}
			break;
			case 4:
				if ((int)fTrackBufferAntideuteron.size()==fTrackBufferSize) {
					fTrackBufferAntideuteron.insert(fTrackBufferAntideuteron.begin(), track);
					fTrackBufferAntideuteron.resize(fTrackBufferSize);
				} else if ((int)fTrackBufferAntideuteron.size()<fTrackBufferSize) {
					fTrackBufferAntideuteron.emplace_back(track);
				}
			break;
			default:
			LOG(warn) << "SetTrackBuffer: invalid particle ID!";
		}
	}

	template <typename T6>
	void FillMixedEventDeltas(const T6& track, int particleType) {
		switch (particleType) {
			case 1:
				if (fTrackBufferProton.size() == 0) return;
				for (int i=0; i<(int)fTrackBufferProton.size(); i++) { // can I do this even if the trackbuffer isn't even full yet?
					double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - fTrackBufferProton[i].phi());
					double DeltaEta = TMath::Abs(track.eta() - fTrackBufferProton[i].eta());
					registryData.fill(HIST("hDeltaPhiME"), particleType, DeltaPhi);
					registryData.fill(HIST("hDeltaPhiEtaME"), particleType, DeltaPhi, DeltaEta);
				}
			break;
			case 2:
				if (fTrackBufferAntiproton.size() == 0) return;
				for (int i=0; i<(int)fTrackBufferAntiproton.size(); i++) {
					double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - fTrackBufferAntiproton[i].phi());
					double DeltaEta = TMath::Abs(track.eta() - fTrackBufferAntiproton[i].eta());
					registryData.fill(HIST("hDeltaPhiME"), particleType, DeltaPhi);
					registryData.fill(HIST("hDeltaPhiEtaME"), particleType, DeltaPhi, DeltaEta);
				}
			break;
			case 3:
				if (fTrackBufferDeuteron.size() == 0) return;
				for (int i=0; i<(int)fTrackBufferDeuteron.size(); i++) {
					double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - fTrackBufferDeuteron[i].phi());
					double DeltaEta = TMath::Abs(track.eta() - fTrackBufferDeuteron[i].eta());
					registryData.fill(HIST("hDeltaPhiME"), particleType, DeltaPhi);
					registryData.fill(HIST("hDeltaPhiEtaME"), particleType, DeltaPhi, DeltaEta);
				}
			break;
			case 4:
				if (fTrackBufferAntideuteron.size() == 0) return;
				for (int i=0; i<(int)fTrackBufferAntideuteron.size(); i++) {
					double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - fTrackBufferAntideuteron[i].phi());
					double DeltaEta = TMath::Abs(track.eta() - fTrackBufferAntideuteron[i].eta());
					registryData.fill(HIST("hDeltaPhiME"), particleType, DeltaPhi);
					registryData.fill(HIST("hDeltaPhiEtaME"), particleType, DeltaPhi, DeltaEta);
				}
			break;
			default:
			LOG(warn) << "FillMixedEventDeltas: invalid particle ID!";
		}
	}

  //****************************************************************************************************

	Filter prelimTrackCuts = (/* aod::track::TPCrefit == fTPCRefit && aod::track::ITSrefit == fITSRefit && */aod::track::itsChi2NCl<fMaxChi2ITS && aod::track::tpcChi2NCl<fMaxChi2TPC && nabs(aod::track::dcaXY) < fMaxDCAxy && nabs(aod::track::dcaZ) < fMaxDCAz && nabs(aod::track::eta) < fMaxEta);

  void process_jet_reco(aod::Collision const& collision, soa::Filtered<FullTracks> const& tracks)
  {
    jetReco(collision, tracks);
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, process_jet_reco, "reconstruct jets", true);
};

//****************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AngularCorrelationsInJets>(cfgc, TaskName{"angular-correlations-in-jets"})};
}
