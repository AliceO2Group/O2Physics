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

/// \file taskCorrelationHfeHadrons.cxx
/// \brief task for azimuthal correlation of heavy flavour hadron decay electron with charge particles
/// \author Rashi gupta <rashi.gupta@cern.ch>, IIT Indore 
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore 

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "Framework/AnalysisDataModel.h"
#include "THnSparse.h"
#include "Common/Core/RecoDecay.h"


using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;
using namespace o2::soa;

using GTrks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra,o2::aod::pidTPCFullEl, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension>;
using GTrk = GTrks::iterator;

struct HfCorrelatorElHadrons{

	Configurable<bool> Isrun3{"Isrun3",true, "Data is from Run3 or Run2"};
	Configurable<bool> fclusterinfo{"fclusterinfo",true, "EMCal cluster info before and after track match"};
	Configurable<int> fNDelPhiBins{"fNDelPhiBins",32, "Bins for #Delta#varphi bins"};
	Configurable<float> etaTrackMin{"etaTrackMin", -0.6f, "Eta range for electron tracks"};
	Configurable<float> etaTrackMinDcright{"etaTrackMinDcright", 0.22f, "Eta range for electron tracks"};
	Configurable<float> etaTrackMinDcleft{"etaTrackMinDcleft", -0.6f, "Eta range for electron tracks"};
	Configurable<float> ptTrackMin{"ptTrackMin", 2.0f , "Transverse MOmentum  range for  electron tracks"};
	Configurable<float> phiTrackMaxEM{"phiTrackMaxEM", 5.708f, "phi range for electron tracks associated Emcal"};
	Configurable<float> phiTrackMinEM{"phiTrackMinEM", 4.5355f, "phi range for electron tracks associated Emcal"};
	Configurable<float> phiTrackMaxDC{"phiTrackMaxDC",3.3621f, "phi range for electron tracks associated Dcal"};
	Configurable<float> phiTrackMinDc{"phiTrackMinDc", 1.3955f, "phi range for electron tracks associated Dcal"};
	Configurable<float> etaTrackMax{"etaTrackMax", 0.6f, "Eta range for electron tracks"};
	Configurable<float> etaTrackMaxDcright{"etaTrackMaxDcright", 0.6f, "Eta range for electron Dcal tracks"};
	Configurable<float> etaTrackMaxDcleft{"etaTrackMaxDcleft", -0.22f, "Eta range for electron Dcal tracks"};
	Configurable<float> eopElMin{"eopElMin", 0.8f, "Minimum E/p for electron tracks"};
	Configurable<float> eopElMax{"eopElMax", 1.2f, "Maximum E/p for electron tracks"};
	Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 0.5f, "DCA XY cut"};
	Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1.0f, "DCA Z cut"};
	Configurable<float> tpcnSigElMin{"tpcnSigElMin", -1.0f, "min Electron TPCnsigma"};
	Configurable<float> tpcnSigElMax{"tpcnSigElMax", 3.0f, "max Electron TPCnsigma"};
	Configurable<float> m02ElMin{"m02ElMin", 0.02f, "min Electron M02"};
	Configurable<float> m02ElMax{"m02ElMax", 0.9f, "max Electron M02"};
	Configurable<float> m20ElMin{"m20ElMin", 0.0f, "min Electron M20"};
	Configurable<float> m20ElMax{"m20ElMax", 1000.f, "max Electron M20"};
	Configurable<float> ptAssocMin{"ptAssocMin",0.2f, "Minimum pT for associated track"};
	Configurable<float> etaAssocTrack{"etaAssocTrack", 0.8f, "Eta range for Associated or partner electron tracks"};
	Configurable<float> deltaphiMinEMC{"deltaphiMinEMC", 0.01f, "Phi distance of EMCAL cluster to its closest track"};
	Configurable<float> deltaetaMinEMC{"deltaetaMinEMC", 0.01f, "Eta distance of EMCAL cluster to its closest track"};
	Configurable<float> timeEMCClsCut{"timeEMCClsCut", 100.f, "Cluster time"};
	Configurable<bool> pTcondition{"pTcondition",true, "Electron pT should be greater than Associate particle pT"};

	PresliceUnsorted<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalmatchedtrack::trackId; 

	HistogramConfigSpec hCorrelSpec{HistType::kTHnSparseD,{{30, 0., 30.}, {20,0.,20.}, {fNDelPhiBins, -TMath::Pi()/2, 3*TMath::Pi()/2}, {50, -1.8, 1.8}}};
	HistogramConfigSpec hTrackInfoSpec{HistType::kTHnSparseD,{{500, 0, 160},{300,-15,15},{500,0.,50.},{500,0.,50.},{100,-1.5,1.5},{100, 0, 7}}};
	HistogramConfigSpec hTrackallInfoSpec{HistType::kTHnSparseD,{{500, 0, 160},{300,-15,15},{500,0.,50.},{500,0.,50.},{100,-1.5,1.5},{100, 0, 7},{2000, -1, 1},{2000 ,-1,1}, {3,0,3}}};
	HistogramConfigSpec hClusterInfoSpec{HistType::kTHnSparseD,{{300, 0.0, 30.0},{100,-0.9,0.9},{200,0,6.3},{50,0,50},{1800,-900,900}}};
	HistogramConfigSpec hPIDSpec{HistType::kTHnSparseD,{{500, 0.0, 50.0},{500,0.,50.},{300,-15,15},{300, 0.0, 30.0},{400,0,2},{400,0,2}}};
	HistogramConfigSpec hDphiDetaclustertrackSpec{HistType::kTH3F, {{400,-0.2,0.2},{400,-0.2,0.2},{280,0,70}}};


	HistogramRegistry registry{
		"registry",
			{{"hNevents", "hNevents", {HistType::kTH1F, {{3,1,4}}}},
				{"zvertex", "z vertex", {HistType::kTH1F, {{100, -20, 20}}}},

				{"fSprsHHCorrl","Sparse for Dphi and Deta with Inclusive electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",hCorrelSpec},
				{"fSprsInclusiveEHCorrl","Sparse for Dphi and Deta with Inclusive electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;", hCorrelSpec},
				{"fSprsdEdxnSigmaPt","Sparse TPC info; dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi; DcaXY;Dcaz;passEMcal; ", hTrackallInfoSpec},
				{"fSprshadroninformation","Sparse hadron info; dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi;",hTrackInfoSpec},
				{"fSprsClusterInfoBe","Cluster Info before match; Energy();#eta;#varphi", hClusterInfoSpec},
				{"fSprsClusterInfoAf","Cluster Info after match; Energy(GeV);#eta;#varphi", hClusterInfoSpec},
				{"fSprsPIDafterMatch","PID Info after match;dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi;",hTrackInfoSpec},
				{"fSprsPIDafterPIDcuts","PID Info after PID cuts; #it{p}(GeV#it{/c});#it{p}_{T}(GeV#it{/c});n_{#sigma}^{e};GeV;M02;M20", hPIDSpec},
				{"ClsTrkEtaPhiDiffTime", "ClsTrkEtaPhiDiffTime;#Delta#eta;#Delta#varphi;Sec;",hDphiDetaclustertrackSpec}}};


	void init(o2::framework::InitContext&)
	{
		registry.get<THnSparse>(HIST("fSprsInclusiveEHCorrl"))->Sumw2();
		registry.get<THnSparse>(HIST("fSprsHHCorrl"))->Sumw2();
		registry.get<THnSparse>(HIST("fSprsdEdxnSigmaPt"))->Sumw2();
		registry.get<THnSparse>(HIST("fSprsClusterInfoBe"))->Sumw2();
		registry.get<THnSparse>(HIST("fSprsClusterInfoAf"))->Sumw2();
		registry.get<THnSparse>(HIST("fSprsPIDafterMatch"))->Sumw2();
		registry.get<THnSparse>(HIST("fSprsPIDafterPIDcuts"))->Sumw2();
		registry.get<THnSparse>(HIST("fSprshadroninformation"))->Sumw2();

	}

	// correlation function for electron hadron
	void Correlation(GTrk const& eTrack,GTrks const& assotrks, int SparseNo = -1)
	{
		if(SparseNo == -1) {cout<<"Error: pass sparse value from '0' to 'N' "<<endl; return;}
		std::shared_ptr<THnSparse> hEHcorrArray[4] = {
			registry.get<THnSparse>(HIST("fSprsHHCorrl")), 
			registry.get<THnSparse>(HIST("fSprsInclusiveEHCorrl"))};

		//Construct Deta Phi between electrons and hadrons

		Double_t  Dphi = -999;
		Double_t  Deta = -999;
		Double_t ptHad =  -999;
		Double_t ptEle =  -999;
		Double_t phiEle =  -999;
		Double_t phiHad =  -999;
		Double_t etaEle =  -999;
		Double_t etaHad =  -999;

		for(const auto& aTrack:assotrks){
			if(aTrack.globalIndex() == eTrack.globalIndex()) continue;

			ptHad = aTrack.pt();
			ptEle = eTrack.pt();
			phiEle = eTrack.phi();
			phiHad = aTrack.phi();
			etaEle = eTrack.eta();
			etaHad = aTrack.eta();

			//Apply Hadron cuts
			if(abs(etaHad) > etaAssocTrack)  continue; 
			if(ptHad < ptAssocMin)  continue;  
			if(abs(aTrack.dcaXY()) > dcaXYTrackMax || abs(eTrack.dcaZ())> dcaZTrackMax) continue;
			if(!aTrack.isGlobalTrackWoDCA())  continue;


			registry.fill(HIST("fSprshadroninformation"), aTrack.tpcSignal(),aTrack.tpcNSigmaEl(),aTrack.p(), ptHad , etaHad, phiHad);

			if(pTcondition && (ptEle > ptHad)) continue;

			Dphi  = RecoDecay::constrainAngle( phiEle - phiHad, -o2::constants::math::PIHalf);
			Deta = etaEle - etaHad;
			if(SparseNo >= 0) hEHcorrArray[SparseNo]->Fill(ptEle, ptHad, Dphi, Deta);
		}
	}

	Filter CollisionFilter = nabs(aod::collision::posZ) < 10.f && aod::collision::numContrib > (uint16_t) 1;
	using aodCollisions = soa::Filtered<soa::Join<aod::Collisions,aod::Mults, aod::EvSels>>::iterator;
	void process(aodCollisions const& collision , aod::EMCALClusters const& mAnalysisClusters, o2::aod::EMCALMatchedTracks const& MatchedTracks, GTrks  const& tracks)                      
	{ 
		if (!(Isrun3 ? collision.sel8() : (collision.sel7() && collision.alias_bit(kINT7)))) return;

		registry.fill(HIST("hNevents"), 1);
		registry.fill(HIST("zvertex"), collision.posZ());

		/////////////////////////////////
		// cluster info before match  ///
		///////////////////////////////
		if(fclusterinfo){
			for (const auto& clusterbf:mAnalysisClusters)
			{
				registry.fill(HIST("fSprsClusterInfoBe"),clusterbf.energy(),clusterbf.eta(),clusterbf.phi(),clusterbf.nCells(),clusterbf.time());
			}
		}
		int PassEMCal;
		Double_t phiTrack = -999;
		Double_t etaTrack = -999;
		Double_t pTrack =  -999;
		Double_t ptTrack =  -999;
		Double_t dcaxyTrack=  -999;
		Double_t dcazTrack =  -999;
		Double_t tpcnsigmaTrack =  -999;

		for (auto& Track : tracks) {

			phiTrack = Track.phi();
			etaTrack = Track.eta();
			pTrack =  Track.p();
			ptTrack =  Track.pt();
			dcaxyTrack=   Track.dcaXY();
			dcazTrack =  Track.dcaZ();
			tpcnsigmaTrack =  Track.tpcNSigmaEl();

			if(!Track.isGlobalTrackWoDCA())  continue;
			PassEMCal = 0;
			if((phiTrack > phiTrackMinEM && phiTrack  < phiTrackMaxEM) &&  (etaTrack > etaTrackMin && etaTrack < etaTrackMax))   PassEMCal = 1; //EMcal acceptance passed
			if ((phiTrack  > phiTrackMinDc && phiTrack  < phiTrackMaxDC) && ((etaTrack > etaTrackMinDcright && etaTrack < etaTrackMaxDcright) || (etaTrack > etaTrackMinDcleft && etaTrack < etaTrackMaxDcleft))) PassEMCal = 2; //Dcal acceptance passed

			registry.fill(HIST("fSprsdEdxnSigmaPt"), Track.tpcSignal(),tpcnsigmaTrack ,pTrack , ptTrack ,etaTrack ,phiTrack , dcaxyTrack , dcazTrack , PassEMCal); //track infor after filter bit

			// Apply Track cut
			if(ptTrack < ptTrackMin)  continue;  
			if(abs(dcaxyTrack) > dcaXYTrackMax || abs(dcazTrack)> dcaZTrackMax) continue;

			auto tracksofcluster = MatchedTracks.sliceBy(perClusterMatchedTracks, Track.globalIndex());
			Double_t phiMatchTrack = -999;
			Double_t etaMatchTrack = -999;
			Double_t pMatchTrack =  -999;
			Double_t ptMatchTrack =  -999;
			Double_t tpcnsigmaMatchTrack =  -999;
			Double_t phiMatchCluster = -999;
			Double_t etaMatchCluster = -999;
			Double_t EMatchCluster =  -999;
			Double_t m02MatchCluster = -999;
			Double_t m20MatchCluster =  -999;
			Double_t timeMatchCluster   =  -999;
			for(const auto& emtrack : tracksofcluster){

				if(Track.globalIndex()!=emtrack.trackId()) continue;
				auto mtrack = emtrack.track_as<GTrks>(); 
				auto cluster = emtrack.emcalcluster_as<aod::EMCALClusters>();

				phiMatchTrack = mtrack.phi();
				etaMatchTrack = mtrack.eta();
				pMatchTrack =  mtrack.p();
				ptMatchTrack =  mtrack.pt();
				tpcnsigmaMatchTrack = mtrack.tpcNSigmaEl();
				phiMatchCluster = cluster.phi();
				etaMatchCluster = cluster.eta();
				EMatchCluster =  cluster.energy();
				m02MatchCluster = cluster.m02();
				m20MatchCluster =  cluster.m20();
				timeMatchCluster  =  cluster.time();
				Correlation(mtrack,tracks,0); //"0" stands for filling Di-hadron

				if(etaMatchTrack < etaTrackMin || etaMatchTrack > etaTrackMax) continue;
				double  deltaPhiEMT = -999.;
				double  deltaEtaEMT = -999.;

				deltaPhiEMT = mtrack.trackPhiEmcal() - phiMatchCluster;
				deltaEtaEMT = mtrack.trackEtaEmcal() - etaMatchCluster;

				registry.fill(HIST("ClsTrkEtaPhiDiffTime"), deltaEtaEMT , deltaPhiEMT, timeMatchCluster  );

				// Track and cluster Matching 
				if(TMath::Abs(deltaPhiEMT) > deltaphiMinEMC || TMath::Abs(deltaEtaEMT )> deltaetaMinEMC) continue;
				if(timeMatchCluster   > timeEMCClsCut) continue;

				if(fclusterinfo) registry.fill(HIST("fSprsClusterInfoAf"), EMatchCluster ,etaMatchCluster, phiMatchCluster ,cluster.nCells(),timeMatchCluster);

				registry.fill(HIST("fSprsPIDafterMatch") , mtrack.tpcSignal(),  tpcnsigmaMatchTrack, pMatchTrack ,ptMatchTrack, etaMatchTrack , phiMatchTrack);

				double  eop = EMatchCluster/ pMatchTrack ;

				//Apply Electron Identification cuts

				if(( tpcnsigmaMatchTrack < tpcnSigElMin || tpcnsigmaMatchTrack  > tpcnSigElMax) || (eop < eopElMin || eop > eopElMax) || (m02MatchCluster < m02ElMin || m02MatchCluster > m02ElMax) || (m20MatchCluster < m20ElMin ||    m20MatchCluster > m20ElMax)) continue;

				registry.fill(HIST("fSprsPIDafterPIDcuts"), pMatchTrack , ptMatchTrack,  tpcnsigmaMatchTrack , EMatchCluster , m02MatchCluster , m20MatchCluster);

				Correlation(mtrack,tracks,1); //"1" stands for filling Electron-hadron      


			}

		}
	}

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return WorkflowSpec{

		adaptAnalysisTask<HfCorrelatorElHadrons>(cfgc),

	};
}
