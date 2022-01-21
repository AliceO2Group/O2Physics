
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>

//centrality
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

////TODO: remove useless:
//#include "ReconstructionDataFormats/Track.h"
//#include "Framework/AnalysisDataModel.h"

#include "Framework/HistogramRegistry.h"

#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
//#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
////

//#include <Math/GenVector/PxPyPzE4D.h>
#include <Math/Vector4D.h>
#include <Math/LorentzVector.h>
#include <TRandom.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using namespace ROOT;
using namespace ROOT::Math;

/*namespace o2::aod
{
}*/

//typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> LorentzVectorD;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

class JFlucCatalyst{
public:
	O2_DEFINE_CONFIGURABLE(zvertex,double,8.0,"Accepted z-vertex range");
	O2_DEFINE_CONFIGURABLE(ptmin,double,0.2,"Minimal pT for tracks");
	O2_DEFINE_CONFIGURABLE(ptmax,double,5.0,"Maximal pT for tracks");
	O2_DEFINE_CONFIGURABLE(charge,int,0,"Particle charge: 0 = all; 1 = positive; -1 = negative");
	O2_DEFINE_CONFIGURABLE(trackingMode,int,0,"Tracking mode: 0 = global; 1 = hybrid");
	Configurable<bool> cutOutliers{"cutOutliers",false,"Cut outlier events"};

	Service<ccdb::BasicCCDBManager> ccdb;
	Configurable<std::string> url{"ccdb-url","http://ccdb-test.cern.ch:8080","CCDB repository URL"};
	O2_DEFINE_CONFIGURABLE(mapCentFlattening,std::string,"","CCDB path to centrality flattening map");
	O2_DEFINE_CONFIGURABLE(mapNUACorrection,std::string,"","CCDB path to NUA correction map");
	Configurable<long> nolaterthan{"ccdb-no-later-than",std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),"Latest acceptable timestamp of creation for the object"};

	TH1D *pcentFlatteningMap = 0;
	TList *pnuaMapList = 0;

	//XXX this will be virtual
	Int_t GetCentBin(Double_t cent){
		return 0;
	}

	void init(InitContext const &ic){
		LOGF(info,"JFlucCatalyst init()");

		ccdb->setURL(url.value);
		ccdb->setCaching(true);
		ccdb->setCreatedNotAfter(nolaterthan.value);

		if(!mapCentFlattening.value.empty()){
			pcentFlatteningMap = ccdb->getForTimeStamp<TH1D>(mapCentFlattening.value,nolaterthan.value);
			if(pcentFlatteningMap)
				LOGF(info,"Centrality flattening enabled. Loaded %s.",mapCentFlattening.value.c_str());
			else LOGF(info,"Failed to load centrality flattening histogram %s. Flattening will be disabled.",mapCentFlattening.value.c_str());
		}
		if(!mapNUACorrection.value.empty()){
			pnuaMapList = ccdb->getForTimeStamp<TList>(mapNUACorrection.value,nolaterthan.value);
			if(pnuaMapList)
				LOGF(info,"NUA correction enabled. Loaded %s.",mapNUACorrection.value.c_str());
			else LOGF(info,"Failed to load NUA correction catalog %s. Correction will be disabled.",mapNUACorrection.value.c_str());
		}

		//fInputList = new TClonesArray("AliJBaseTrack",2500);
		//fInputList = new TClonesArray("PtEtaPhiEVector",2500);
		//fInputList->SetOwner(kTRUE);
		fInputList = new std::vector<PtEtaPhiEVector>();
		fInputList->reserve(2500);

		collisionId = 0;
	}

	//void process(aod::Collision const& collision, aod::Tracks const& tracks){
	void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms, aod::CentRun2CL0s>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TOFSignal> const& tracks){
		if(!collision.alias()[kINT7] || !collision.sel7())
			return;
		if(std::abs(collision.posZ()) > zvertex)
			return;

		Double_t cent = collision.centV0M();
		if(pcentFlatteningMap){
			Int_t bin = pcentFlatteningMap->GetXaxis()->FindBin(cent);
			if(gRandom->Uniform(0,1) > pcentFlatteningMap->GetBinContent(bin))
				return;
		}

		//TODO: outlier cutting
		UInt_t FB32Tracks = 0;
		UInt_t FB32TOFTracks = 0;
		for(auto &track : tracks){
			//fb
			Double_t tofTime = track.tofSignal();
			//track.tofDz();
			//if(track.hasTOF() && std::abs
			//track.isGlobalTrackSDD(); //hybrid
		}

		if(cutOutliers.value){
			//--
			/*
			const AliVVertex* vtTrc = event->GetPrimaryVertex();
			const AliVVertex* vtSPD = event->GetPrimaryVertexSPD();
			double covTrc[6],covSPD[6];
			vtTrc->GetCovarianceMatrix(covTrc);
			vtSPD->GetCovarianceMatrix(covSPD);
			double dz = vtTrc->GetZ()-vtSPD->GetZ();
			double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
			double errTrc = TMath::Sqrt(covTrc[5]);
			double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
			if(TMath::Abs(dz) > 0.2 || nsigTot > 10 || nsigTrc > 20)
				return kFALSE;

			AliMultSelection *pms = (AliMultSelection*)event->FindListObject("MultSelection");
			if(!pms){
				AliError("MultSelection unavailable.");
				return kFALSE;
			}*/

			//TODO: how to get SPD primary vertex?
			collision.posX();
			collision.covXX();//[5] = covYZ
			auto primaryVertex = getPrimaryVertex(collision);
			primaryVertex.getZ();
			auto covMatrix = primaryVertex.getCov();
			//auto primaryVertexSPD = getPrimaryVertexSPD(collision);

			double centCL0 = collision.centRun2CL0();
			double center = 0.973488*centCL0+0.0157497;
			double sigma = 0.673612+centCL0*(0.0290718+centCL0*(-0.000546728+centCL0*5.82749e-06));
			if(cent < center-5.0*sigma || cent > center+5.5*sigma || cent < 0.0 || cent > 60.0)
				return;
		}

		Int_t cbin = GetCentBin(cent);
		if(cbin < 0)
			return;

		auto bc = collision.bc_as<aod::BCsWithTimestamps>();
		TH1 *pweightMap = pnuaMapList?
			(TH1*)pnuaMapList->FindObject(Form("PhiWeights_%u_%02u",bc.runNumber(),cbin)):0;

		
		fInputList->clear();
		for(auto &track : tracks){
			//if(!track.hasTOF())
			//	continue;

			if(trackingMode == 0 && !track.isGlobalTrack())
				continue;
			else
			if(trackingMode == 1 && !track.isGlobalTrackSDD())
				continue;

			Double_t pt = track.pt();
			if(pt < ptmin || pt > ptmax)
				continue;

			Int_t ch = track.sign();
			if(charge != 0 && charge*ch < 0)
				continue;

			Double_t eta = track.eta();
			Double_t phi = track.phi();

			Double_t phiWeight = 1.0;
			if(pweightMap){
				Int_t bin = pweightMap->FindBin(phi,eta,zvertex);
				phiWeight = pweightMap->GetBinContent(bin);
			}

			//AliJBaseTrack *ptrack = new ((*fInputList)[ntrack++])AliJBaseTrack;
			//PtEtaPhiEVector *ptrack = ((*fInputList)[ntrack++]) PtEtaPhiEVector;
			fInputList->emplace_back();
			PtEtaPhiEVector t = fInputList->back();
			
			//ptrack->SetLabel(ntrack);
			//ptrack->SetParticleType(0);
			t.SetPt(pt);
			t.SetEta(eta);
			t.SetPhi(phi);
			//ptrack->SetE(track.E());
			//ptrack->SetCharge(ch);
			//ptrack->SetTrackEff(1.0);
			//ptrack->SetWeight(1.0); // phi weight
		}
		++collisionId;
		LOGF(info,"event %u processed with %u tracks.",collisionId,fInputList->size());
	}

	std::vector<PtEtaPhiEVector> *fInputList;
	int collisionId;

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
	return WorkflowSpec{
		adaptAnalysisTask<JFlucCatalyst>(cfgc)
	};
}

