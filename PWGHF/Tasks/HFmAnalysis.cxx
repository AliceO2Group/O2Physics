#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Common/DataModel/EventSelection.h"
#include <TH1F.h>
#include <TH2I.h>
#include <THashList.h>
#include <TString.h>
#include <vector>
#include <memory>
#include <cstring>
#include <THnSparse.h>
#include "TDatabasePDG.h"
//=============================================================================
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;

//TODO
// Explicit name space

//namespace o2::framework::expressions;

namespace o2::aod {
  namespace dqanalysisflags {
    DECLARE_SOA_COLUMN(IsEventSelected,isEventSelected,int);
  }

  DECLARE_SOA_TABLE(EventCuts,"AOD","EVENTCUTS",dqanalysisflags::IsEventSelected);
}
//=============================================================================

using theCollisions    = o2::soa::Join<aod::Collisions,aod::EvSels>;
using theMcCollisions  = o2::soa::Join<aod::Collisions,aod::McCollisionLabels,aod::EvSels>;
using MyEventsSelected = o2::soa::Join<aod::Collisions,aod::EvSels,o2::aod::EventCuts>;
using McEventsSelected = o2::soa::Join<aod::Collisions,aod::McCollisionLabels,aod::EvSels,o2::aod::EventCuts>;
using MyMuons          = o2::soa::Join<aod::FullFwdTracks,aod::McFwdTrackLabels>; 
//=============================================================================

constexpr static uint32_t gkMuonFillMap  (VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov);
constexpr static uint32_t gkEventFillMap (VarManager::ObjTypes::BC   | VarManager::ObjTypes::Collision);
constexpr static uint32_t gkTrackMCFillMap (VarManager::ObjTypes::ParticleMC);

void DefineHistograms(HistogramManager *, TString);
//=============================================================================

struct EventSelection {
  Configurable<bool>  fConfigTriggerUsed     { "cfgtriggerUsed", false, "whether apply the software trigger" };
  Configurable<bool>  fConfigVtxZUsed        { "cfgVtxZUsed", false, "whether apply the VtxZ trigger" };
  Configurable<float> fConfigSoftwareTrigger { "cfgsoftwaretrigger", VarManager::kIsMuonSingleLowPt7, "softerware tigger flag" };

  Configurable<float> fConfigVtxZMIN { "cfgVtxZMIN", -10., "min. z of primary vertex [cm]" };
  Configurable<float> fConfigVtxZMAX { "cfgVtxZMAX",  10., "max. z of primary vertex [cm]" };
//=============================================================================

  HistogramManager *fHistMan;
  OutputObj<THashList> fOutputList{ "output" };
//=============================================================================

  float *fValues;
  AnalysisCompositeCut *fEventCut;

  void init(o2::framework::InitContext &) {
    VarManager::SetDefaultVarNames();
    fValues = new float[VarManager::kNVars];
    fEventCut = new AnalysisCompositeCut("event selection","Event Selection",true);

    AnalysisCut fCut;
    if (fConfigVtxZUsed) fCut.AddCut(VarManager::kVtxZ,fConfigVtxZMIN,fConfigVtxZMAX,false);
    if (fConfigTriggerUsed) fCut.AddCut(fConfigSoftwareTrigger,0.5,1.5,false);
    fEventCut->AddCut(&fCut);

    fHistMan = new HistogramManager("analysisHistos","aa",VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames,VarManager::fgVariableUnits);

    DefineHistograms(fHistMan,"EventBeforeCuts;EventAfterCuts;");
    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }
//=============================================================================

  Produces<aod::EventCuts> eventSel;
  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSel(TEvent const &event, aod::BCs const &bcs) {
    VarManager::ResetValues(0,VarManager::kNEventWiseVariables);
    VarManager::FillEvent<TEventFillMap>(event,fValues);
    fHistMan->FillHistClass("EventBeforeCuts",fValues);

    if (fEventCut->IsSelected(fValues)) {
      fHistMan->FillHistClass("EventAfterCuts",fValues);
      eventSel(1);
    } else {
      eventSel(0);
    }
  }
//=============================================================================

  void processEventDataAO2D(theCollisions::iterator const &event, aod::BCs const &bcs) {
    runEventSel<gkEventFillMap>(event, bcs);
  }
//=============================================================================

  void processEventMCAO2D(theMcCollisions::iterator const &event, aod::BCs const &bcs) {
    runEventSel<gkEventFillMap>(event, bcs);
  }
//=============================================================================

  PROCESS_SWITCH(EventSelection, processEventDataAO2D, "run event selection with Data AO2D file", true);
  PROCESS_SWITCH(EventSelection, processEventMCAO2D, "run event selection with MC AO2D file", false);
};
//=============================================================================

struct MuonSelection {
  Configurable<std::string> fConfigCuts { "cfgMuonCuts", "muonQualityCuts", "muon selection" };
//=============================================================================

  float *fValues;
  AnalysisCompositeCut* fTrackCut;

  o2::framework::HistogramRegistry registry {
    "registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true };
//=============================================================================

  void init(o2::framework::InitContext &) {
    AxisSpec apT  { 200,  0., 200., "#it{p}_{T} (GeV/#it{c})" };
    AxisSpec aEta { 500, -5.,   5., "#eta" };
    AxisSpec aDCA { 500,  0.,   5., "DCA (cm)" };
    AxisSpec aSign {  5, -2.5,   2.5, "Charge"};
    AxisSpec aP	   {  500,  0., 500., "p (GeV/#it{c})"};
    AxisSpec aVtxZ {  400,-20.,  20., "Vertex Z (cm)"};

    AxisSpec aTrkType {   5,  -0.5,   4.5, "TrackType" };
    AxisSpec aMCmask  { 200,  -0.5, 199.5, "mcMask"    };
   
    AxisSpec apTGen  { 200,  0., 200., "#it{p}_{T} Truth (GeV/#it{c})" };
    AxisSpec aEtaGen { 500, -5.,   5., "#eta_{Truth}" };
    AxisSpec aPGen   { 500,  0., 500., "p_{Truth} (GeV/#it{c})"};
    AxisSpec apTDif  { 200,  -2., 2., "#it{p}_{T} diff (GeV/#it{c})" };
    AxisSpec aEtaDif { 200, -2.,   2., "#eta diff" };
    AxisSpec aPDif   { 200,  -2., 2., "p diff (GeV/#it{c})"};

    HistogramConfigSpec cfgTHnMu { HistType::kTHnSparseD, { apT, aEta, aDCA, aSign, aP, aVtxZ, aTrkType, aMCmask}, 8 };
    HistogramConfigSpec cfgTHnPt { HistType::kTHnSparseD, { apT, apTGen, apTDif, aTrkType}, 4 };
    HistogramConfigSpec cfgTHnEta { HistType::kTHnSparseD, { aEta, aEtaGen, aEtaDif, aTrkType}, 4 };
    HistogramConfigSpec cfgTHnP { HistType::kTHnSparseD, { aP, aPGen, aPDif, aTrkType}, 4 };

    registry.add("muBcuts", "", cfgTHnMu);
    registry.add("muAcuts", "", cfgTHnMu);
    registry.add("pTBcuts", "", cfgTHnPt);
    registry.add("pTAcuts", "", cfgTHnPt);
    registry.add("etaBcuts", "", cfgTHnEta);
    registry.add("etaAcuts", "", cfgTHnEta);
    registry.add("pBcuts", "", cfgTHnP);
    registry.add("pAcuts", "", cfgTHnP);

    VarManager::SetDefaultVarNames();
    fValues = new float[VarManager::kNVars];

    fTrackCut = new AnalysisCompositeCut(true);
    TString selectStr = fConfigCuts.value;
    fTrackCut->AddCut(dqcuts::GetAnalysisCut(selectStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);
  }
//=============================================================================
  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuons>
  void runDataMuonSel(TEvent const &event, aod::BCs const &bcs, TMuons const &tracks) {
    if(event.isEventSelected()==0) return;

    VarManager::ResetValues(0,VarManager::kNMuonTrackVariables,fValues);
    VarManager::FillEvent<gkEventFillMap>(event,fValues);

    for(auto &track: tracks) {
      VarManager::FillTrack<TMuonFillMap>(track,fValues);

//compute DCAXY
      std::vector<double> vcovs {
        track.cXX(),   track.cXY(),   track.cYY(),
        track.cPhiX(), track.cPhiY(), track.cPhiPhi(),
        track.cTglX(), track.cTglY(), track.cTglPhi(), track.cTglTgl(),
        track.c1PtX(), track.c1PtY(), track.c1PtPhi(), track.c1PtTgl(), track.c1Pt21Pt2() };

      SMatrix55 tcovs(vcovs.begin(),vcovs.end());
      SMatrix5  tpars(track.x(),track.y(),track.phi(),track.tgl(),track.signed1Pt());

      o2::track::TrackParCovFwd tparcov{track.z(),tpars,tcovs,track.chi2()};
      tparcov.propagateToZlinear(event.posZ());

      const auto dcaX(tparcov.getX()-event.posX());
      const auto dcaY(tparcov.getY()-event.posY());
      const auto dcaXY(std::sqrt(dcaX*dcaX+dcaY*dcaY));

      registry.fill(HIST("muBcuts"),
        fValues[VarManager::kPt],
        fValues[VarManager::kEta], dcaXY,
	fValues[VarManager::kCharge], track.p(),
	fValues[VarManager::kVtxZ],
        fValues[VarManager::kMuonTrackType],0);

      if (fTrackCut->IsSelected(fValues)) registry.fill(HIST("muAcuts"),
        fValues[VarManager::kPt],
        fValues[VarManager::kEta], dcaXY,
        fValues[VarManager::kCharge], track.p(),
	fValues[VarManager::kVtxZ],
        fValues[VarManager::kMuonTrackType],0);	
	}
  }
//=============================================================================


  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, uint32_t TTrackMCFillMap, typename TEvent, typename TMuons, typename TMC>
  void runMCMuonSel(TEvent const &event, aod::BCs const &bcs, TMuons const &tracks, TMC const &mc) {
    if(event.isEventSelected()==0) return;

    VarManager::ResetValues(0,VarManager::kNMuonTrackVariables,fValues);
    VarManager::FillEvent<gkEventFillMap>(event,fValues);

    for(auto &track: tracks) {
      VarManager::FillTrack<TMuonFillMap>(track,fValues);

	if(!track.has_mcParticle()){continue;}
	auto mcp=track.mcParticle();

	VarManager::FillTrack<TTrackMCFillMap>(mcp,fValues);
	
//compute DCAXY
      std::vector<double> vcovs {
        track.cXX(),   track.cXY(),   track.cYY(),
        track.cPhiX(), track.cPhiY(), track.cPhiPhi(),
        track.cTglX(), track.cTglY(), track.cTglPhi(), track.cTglTgl(),
        track.c1PtX(), track.c1PtY(), track.c1PtPhi(), track.c1PtTgl(), track.c1Pt21Pt2() };

      SMatrix55 tcovs(vcovs.begin(),vcovs.end());
      SMatrix5  tpars(track.x(),track.y(),track.phi(),track.tgl(),track.signed1Pt());

      o2::track::TrackParCovFwd tparcov{track.z(),tpars,tcovs,track.chi2()};
      tparcov.propagateToZlinear(event.posZ());

      const auto dcaX(tparcov.getX()-event.posX());
      const auto dcaY(tparcov.getY()-event.posY());
      const auto dcaXY(std::sqrt(dcaX*dcaX+dcaY*dcaY));

      registry.fill(HIST("muBcuts"),
        fValues[VarManager::kPt],
        fValues[VarManager::kEta], dcaXY,
	fValues[VarManager::kCharge], track.p(),
	fValues[VarManager::kVtxZ],
        fValues[VarManager::kMuonTrackType],
        track.mcMask());

      registry.fill(HIST("pTBcuts"),
	fValues[VarManager::kPt], fValues[VarManager::kMCPt], fValues[VarManager::kMCPt]-fValues[VarManager::kPt], fValues[VarManager::kMuonTrackType]);
      registry.fill(HIST("etaBcuts"),
	fValues[VarManager::kEta], fValues[VarManager::kMCEta], fValues[VarManager::kMCEta]-fValues[VarManager::kEta], fValues[VarManager::kMuonTrackType]);
      registry.fill(HIST("pBcuts"),
	fValues[VarManager::kP], fValues[VarManager::kMCPt]*std::cosh(fValues[VarManager::kMCEta]), fValues[VarManager::kMCPt]*std::cosh(fValues[VarManager::kMCEta])-fValues[VarManager::kP], fValues[VarManager::kMuonTrackType]);


      if (fTrackCut->IsSelected(fValues)){ registry.fill(HIST("muAcuts"),
        fValues[VarManager::kPt],
        fValues[VarManager::kEta], dcaXY,
        fValues[VarManager::kCharge], track.p(),
	fValues[VarManager::kVtxZ],
        fValues[VarManager::kMuonTrackType],
        track.mcMask());
	
      registry.fill(HIST("pTAcuts"),
	fValues[VarManager::kPt], fValues[VarManager::kMCPt], fValues[VarManager::kMCPt]-fValues[VarManager::kPt], fValues[VarManager::kMuonTrackType]);
      registry.fill(HIST("etaAcuts"),
	fValues[VarManager::kEta], fValues[VarManager::kMCEta], fValues[VarManager::kMCEta]-fValues[VarManager::kEta], fValues[VarManager::kMuonTrackType]);
      registry.fill(HIST("pAcuts"),
	fValues[VarManager::kP], fValues[VarManager::kMCPt]*std::cosh(fValues[VarManager::kMCEta]), fValues[VarManager::kMCPt]*std::cosh(fValues[VarManager::kMCEta])-fValues[VarManager::kP], fValues[VarManager::kMuonTrackType]);
		}
	}
  }
//=============================================================================

  void processMuonDataAO2D(MyEventsSelected::iterator const &event, aod::BCs const &bcs,
  aod::FullFwdTracks const &tracks) {
	runDataMuonSel<gkEventFillMap,gkMuonFillMap>(event, bcs, tracks);
  }
//=============================================================================
  void processMuonMCAO2D(McEventsSelected::iterator const &event, aod::BCs const &bcs,
  MyMuons const &tracks, aod::McParticles_001 const &mc) {
	runMCMuonSel<gkEventFillMap,gkMuonFillMap,gkTrackMCFillMap>(event, bcs, tracks, mc);
  }
//=============================================================================

  PROCESS_SWITCH(MuonSelection,processMuonDataAO2D,"run muonselection with Data AO2D file",true);
  PROCESS_SWITCH(MuonSelection,processMuonMCAO2D,"run muonselection with MC AO2D file",false);
};
//=============================================================================

WorkflowSpec defineDataProcessing(ConfigContext const &cfgc) {
  return WorkflowSpec {
    adaptAnalysisTask<EventSelection>(cfgc),
    adaptAnalysisTask<MuonSelection>(cfgc),
  };
}
//=============================================================================

void DefineHistograms(HistogramManager *histMan, TString histClasses) {
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));

  for (Int_t iclass=0; iclass<objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    if (classStr.Contains("Event")) dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,all");
  }
}
//=============================================================================
