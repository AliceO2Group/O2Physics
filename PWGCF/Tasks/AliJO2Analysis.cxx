
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"


#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "ReconstructionDataFormats/V0.h"

#include "AliJO2Catalyst.h"
#include "AliJFFlucAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

class JFlucAnalysis{
public:
	~JFlucAnalysis(){
		delete pcf;
	}

	O2_DEFINE_CONFIGURABLE(etamin,double,0.4,"Minimal eta for tracks");
	O2_DEFINE_CONFIGURABLE(etamax,double,0.8,"Maximal eta for tracks");

	//OutputObj<TDirectory> output{TDirectory("jflucO2","jflucO2","",0)};
	OutputObj<TDirectory> output{"jflucO2"};

	void init(InitContext const &ic){
		//
		fInputList.reserve(2500);

		pcf = new AliJFFlucAnalysis("jflucAnalysis");
		pcf->SetNumBins(sizeof(jflucCentBins)/sizeof(jflucCentBins[0]));
		pcf->AddFlags(AliJFFlucAnalysis::FLUC_EBE_WEIGHTING);

		output->cd();
		pcf->UserCreateOutputObjects();
	}

	//void process(aod::Collision const& collision, aod::ParticleTrack const& tracks){
	//void process(soa::Join<aod::Collisions, aod::CentV0Ms>::iterator const& collision, aod::ParticleTrack const& tracks){
	void process(soa::Join<aod::Collisions, aod::CollisionData>::iterator const& collision, aod::ParticleTrack const& tracks){
		if(tracks.size() == 0)
			return; //rejected event
		fInputList.clear();

		for(auto &track : tracks){
			fInputList.emplace_back();
			PtEtaPhiEVector &t = fInputList.back();
			
			//ptrack->SetLabel(ntrack);
			//ptrack->SetParticleType(0);
			t.SetPt(track.pt());
			t.SetEta(track.eta());
			t.SetPhi(track.phi());
		}

		const double fVertex[3] = {collision.posX(),collision.posY(),collision.posZ()};

		pcf->Init();
		pcf->SetInputList(&fInputList);
		pcf->SetEventCentralityAndBin(collision.cent(),collision.cbin());
		pcf->SetEventVertex(fVertex);
		pcf->SetEtaRange(etamin,etamax);
		pcf->UserExec("");

	}
	std::vector<PtEtaPhiEVector> fInputList;
	AliJFFlucAnalysis *pcf;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
	return WorkflowSpec{
		adaptAnalysisTask<JFlucAnalysis>(cfgc)
	};
}

