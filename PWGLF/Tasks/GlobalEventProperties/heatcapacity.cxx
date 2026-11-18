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
///
/// \brief 
/// \author everyone

#include <Framework/AnalysisDataModel.h> // tpcNClsCrossedRows
#include <Framework/AnalysisTask.h> // 
#include <Framework/runDataProcessing.h> // 
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/Configurable.h>
#include <Framework/AnalysisHelpers.h>

#include <Common/DataModel/TrackSelectionTables.h> // TracksDCA
#include <Common/DataModel/EventSelection.h> // ft0cOccupancyInTimeRange,trackOccupancyInTimeRange
#include <Common/DataModel/Centrality.h> // 
#include <Common/DataModel/Multiplicity.h> // multNTracksPVeta1, multNTracksGlobal

#include "Common/CCDB/EventSelectionParams.h" //kIsTriggerTVX, kNoSameBunchPileup, kIsGoodZvtxFT0vsPV

#include <string>
#include <vector>
#include <cmath>


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions; 

using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::PVMults, aod::MultsGlobal, aod::Mults>; 
using TrackCandidates = soa::Join<aod::Tracks, aod::TracksDCA, aod::TracksExtra>; 

struct CentralityClassificationPbPb
{

  struct EventSelections : ConfigurableGroup {
    std::string prefix = "EventSelections-Configurations"; // JSON group name
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC (Run 3 only)"};
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference (Run 3 only)"};
    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", true, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};

    Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
    Configurable<float> maxCent{"maxCent", 100.0f, "max centrality percentil"};
    Configurable<float> minCent{"minCent", 0.0f, "min centrality percentil"};
    
  } eventSelections;

  
  struct TrackSelections : ConfigurableGroup {
    std::string prefix = "TrackSelections-Configurations"; // JSON group name
    Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
  
    Configurable<bool> requireTPCRefit{"RequireTPCRefit", true, "Require TPC Refit"};
    Configurable<bool> requireITSRefit{"RequireITSRefit", true, "Require ITS Refit"};
    Configurable<bool> requireITSibAny{"requireITSibAny", true, "Require ITS inner barrel (ib) any hits"};

    //Configurable<float> dcaZpos{"dcaZpos", 0.2, "DCA Pos To PV (cm)"};
    Configurable<float> dcaZconst{"dcaZconst", 0.0028f, "Constant term for DCAz cut (cm)"};
    Configurable<float> dcaZpt{"dcaZpt", 0.0175f, "pT-dependent coefficient for DCAz cut (cm)"};
    Configurable<float> dcaZexp{"dcaZexp", -0.88f, "Exponent for DCAz cut"};
    Configurable<float> dcaXYconst{"dcaXYconst", 0.0035f, "Constant term for DCAxy cut (cm)"};
    Configurable<float> dcaXYpt{"dcaXYpt", 0.0145f, "pT-dependent coefficient for DCAxy cut (cm)"};
    Configurable<float> dcaXYexp{"dcaXYexp", 0.91f, "Exponent for DCAxy cut"};
    Configurable<float> maxeta{"maxeta", 2.0f, "pseudorapidity cut"};
    Configurable<float> minpT{"minpT", 0.15f, "Min pT cut"};
    Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "Max Chi2 per cluster TPC "};
    Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "Max Chi2 per cluster ITS "};
    Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  } trackSelections;
  

  struct AxisConfigurations : ConfigurableGroup {
    std::string prefix = "axis-Configurations"; // JSON group name
    //ConfigurableAxis axisPt{"axispT", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisEta{"axisEta", {100, -2.0f, 2.0f}, "Eta"};
    ConfigurableAxis axisPhi{"axisPhi", {100, 0, 2 * M_PI}, "  Phi"};
    ConfigurableAxis axispT{"axispT", {100, 0.15f, 10}, "  pT"};
    ConfigurableAxis axisdcaZ{"axisdcaZ",{2000, -3.0f, 3.0f}, "dcaZ"};
    ConfigurableAxis axisdcaXY{"axisdcaXY",{2000, -1.0f, 1.0f}, "dcaXY"};
    ConfigurableAxis axisEvSelection{"axisEvSelection",{10, 0.5f, 10.5f}, ""};
    ConfigurableAxis axisTrackSelection{"axisTrackSelection",{8, 0.5f, 8.5f}, ""};
    ConfigurableAxis axisCent{"axisCent",{100, 0, 100}, ""};
    //ConfigurableAxis axisUltCent{"axisUltCent",{VARIABLE_WIDTH, 0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0}, ""};
    ConfigurableAxis axisMult{"axisMult",{500, 0.0f, 5000.0f}, "Multiplicity"};
    //ConfigurableAxis axisAmplitude{"axisAmplitude", {VARIABLE_WIDTH, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140,160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000,2500, 3000, 3500,4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000, 12000, 14000, 16000, 18000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,60000, 70000}, "FT0C Amplitude"};
    ConfigurableAxis axisAmplitude{"axisAmplitude",{1400, 0.0f, 70000.0f},"FT0C Amplitude"};
    ConfigurableAxis axisVz{"axisVz", {200, -10.0f, 10.0f}, "V_{z} (cm)"}; //nuevo
  } axisConfigurations;



  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};


  void init(InitContext const&)
  {
    histos.add("hEventSelection", "hEventSelection", kTH1D, {axisConfigurations.axisEvSelection});
    histos.add("hTrackSelection", "hTrackSelection", kTH1D, {axisConfigurations.axisTrackSelection});

    histos.add("hEventNchCorrelationBefCuts", "hEventNchCorrelationBefCuts", kTH2F, {{5000, 0, 5000}, {2500, 0, 2500}});
    histos.add("hEventNchCorrelation", "hEventNchCorrelation", kTH2F, {{5000, 0, 5000}, {2500, 0, 2500}});

    histos.add("hEventPVcontributorsVsCentralityBefCuts", "hEventPVcontributorsVsCentralityBefCuts", kTH2F, {{100, 0, 100}, {5000, 0, 5000}});
    histos.add("hEventPVcontributorsVsCentrality", "hEventPVcontributorsVsCentrality", kTH2F, {{100, 0, 100}, {5000, 0, 5000}});
    
    histos.add("hEventGlobalTracksVsCentralityBefCuts", "hEventGlobalTracksVsCentralityBefCuts", kTH2F, {{100, 0, 100}, {2500, 0, 2500}});
    histos.add("hEventGlobalTracksVsCentrality", "hEventGlobalTracksVsCentrality", kTH2F, {{100, 0, 100}, {2500, 0, 2500}});
  
    histos.add("hNTracksVSetaBefCuts", "hNTracksVSetaBefCuts", kTH1D, {axisConfigurations.axisEta});
    histos.add("hNTracksVSeta", "hNTracksVSeta", kTH1D, {axisConfigurations.axisEta});

    histos.add("hNTracksVSphiBefCuts", "hNTracksVSphiBefCuts", kTH1D, {axisConfigurations.axisPhi}); 
    histos.add("hNTracksVSphi", "hNTracksVSphi", kTH1D, {axisConfigurations.axisPhi}); 

    histos.add("hdcaZVSptBefCuts", "hdcaZVSptBefCuts", kTH2F, {axisConfigurations.axisdcaZ, axisConfigurations.axispT});
    histos.add("hdcaZVSpt", "hdcaZVSpt", kTH2F, {axisConfigurations.axisdcaZ, axisConfigurations.axispT});

    histos.add("hdcaXYVSptBefCuts", "hdcaXYVSptBefCut", kTH2F, {axisConfigurations.axisdcaXY, axisConfigurations.axispT});
    histos.add("hdcaXYVSpt", "hdcaXYVSpt", kTH2F, {axisConfigurations.axisdcaXY, axisConfigurations.axispT});

    histos.add("hFT0CAmplitude", "hFT0CAmplitude", kTH1D, {axisConfigurations.axisAmplitude}); 
    //histos.add("hUltCent", "hUltCent", kTH1F, {axisConfigurations.axisUltCent});
    histos.add("hFT0CAmplitudeVSMult", "hFT0CAmplitudeVSMult", kTH2F, {axisConfigurations.axisAmplitude, axisConfigurations.axisMult});
    histos.add("hFT0CAmplitudeVSpT", "hFT0CAmplitudeVSpT", kTH2F, {axisConfigurations.axisAmplitude, axisConfigurations.axispT});
    histos.add("hFT0CAmplitudeVsVz", "hFT0CAmplitudeVsVz", kTH2F, {axisConfigurations.axisVz, axisConfigurations.axisAmplitude});
    
    histos.add("hMultSelectedTracks", "Multiplicity of selected tracks;N_{tracks};Events", kTH1D, {axisConfigurations.axisMult});
    histos.add("hMultGlobalTracks", "Global tracks multiplicity;multNTracksGlobal;Events", kTH1D, {axisConfigurations.axisMult});
    histos.add("hMultPVeta1", "PV contributors (|#eta| < 1);multNTracksPVeta1;Events", kTH1D, {axisConfigurations.axisMult});


    auto hEvSelCount = histos.get<TH1>(HIST("hEventSelection"));
    hEvSelCount->GetXaxis()->SetBinLabel(1, "All events");
    hEvSelCount->GetXaxis()->SetBinLabel(2, "sel8");
    hEvSelCount->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
    hEvSelCount->GetXaxis()->SetBinLabel(4, "Vertex Z position");
    hEvSelCount->GetXaxis()->SetBinLabel(5, "kIsGoodZvtxFT0vsPV");
    hEvSelCount->GetXaxis()->SetBinLabel(6, "kNoSameBunchPileup");
    hEvSelCount->GetXaxis()->SetBinLabel(7, "Below min occup."); 
    hEvSelCount->GetXaxis()->SetBinLabel(8, "Above max occup.");
    hEvSelCount->GetXaxis()->SetBinLabel(9, "Centrality");
    hEvSelCount->GetXaxis()->SetBinLabel(10, "Selected events");

    auto hTrackSelCount = histos.get<TH1>(HIST("hTrackSelection"));
    hTrackSelCount->GetXaxis()->SetBinLabel(1, "All tracks");
    hTrackSelCount->GetXaxis()->SetBinLabel(2, "TPC tracks");
    hTrackSelCount->GetXaxis()->SetBinLabel(3, "ITS tracks");
    hTrackSelCount->GetXaxis()->SetBinLabel(4, "eta selection");
    hTrackSelCount->GetXaxis()->SetBinLabel(5, "min pT selection");
    hTrackSelCount->GetXaxis()->SetBinLabel(6, "dcaXY");
    hTrackSelCount->GetXaxis()->SetBinLabel(7, "dcaZ");
    hTrackSelCount->GetXaxis()->SetBinLabel(8, "Selected tracks");

  }

  template <typename TCollision>
  bool isEventAccepted(const TCollision& collision, bool fillHists)
  {
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 1);
      histos.fill(HIST("hEventNchCorrelationBefCuts"), collision.multNTracksPVeta1(), collision.multNTracksGlobal());
      histos.fill(HIST("hEventPVcontributorsVsCentralityBefCuts"), collision.centFT0C(), collision.multNTracksPVeta1());
      histos.fill(HIST("hEventGlobalTracksVsCentralityBefCuts"), collision.centFT0C(), collision.multNTracksGlobal());
    }

    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 2);
    }

    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 3 );
    }

    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 4);
    }

    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 5 );
    }

    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 6);
    }

    float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 7);
    }

    if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 8);
    }
    /*  PONER CONFIGURABLES */
    if (collision.centFT0C() < eventSelections.minCent || collision.centFT0C() > eventSelections.maxCent) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 9);
    }
      

    if (fillHists) {
      //Arfer Event Selection
      histos.fill(HIST("hEventSelection"), 10);
      histos.fill(HIST("hEventNchCorrelation"), collision.multNTracksPVeta1(), collision.multNTracksGlobal());
      histos.fill(HIST("hEventPVcontributorsVsCentrality"), collision.centFT0C(), collision.multNTracksPVeta1());
      histos.fill(HIST("hEventGlobalTracksVsCentrality"), collision.centFT0C(), collision.multNTracksGlobal());
    }
  
    return true;
  }


  template <typename CheckTrack>
  bool IsTrackSelected(CheckTrack const& track, bool fillHists)
  {
    if (fillHists) {
      histos.fill(HIST("hTrackSelection"), 1);
    }

    //TPC cuts
    if (trackSelections.requireTPCRefit && !track.hasTPC()) {
      return false;
    } 
    if (track.hasTPC() && (track.tpcNClsCrossedRows() < trackSelections.minTPCrows || track.tpcChi2NCl() > trackSelections.maxChi2TPC || track.tpcCrossedRowsOverFindableCls() < trackSelections.minNCrossedRowsOverFindableClustersTPC)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hTrackSelection"), 2);
    }

    //ITS cuts
    if (trackSelections.requireITSRefit && !track.hasITS()) {
      return false;
    }
    if (trackSelections.requireITSibAny && !(track.itsClusterMap() & 0x7)) {
      return false;
    }
    if (track.hasITS() && track.itsChi2NCl() > trackSelections.maxChi2ITS) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hTrackSelection"), 3);
    }

    // eta cut
    if(fabs(track.eta()) > trackSelections.maxeta){
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hTrackSelection"), 4);
    }

    // pT cut
    if(track.pt() < trackSelections.minpT){
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hTrackSelection"), 5);
    }

    // dcaXY cut
    float maxDCAxy = trackSelections.dcaXYconst + trackSelections.dcaXYpt * std::pow(track.pt(), trackSelections.dcaXYexp);
    if(fabs(track.dcaXY()) > maxDCAxy){
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hTrackSelection"), 6);
    }

    //dcaZ VS pt before max dcaZ=0.2 selection, but after dcaXY
    if (fillHists) {
      histos.fill(HIST("hdcaZVSptBefCuts"), track.dcaZ(), track.pt());
    }

     // dcaZ cut
    float maxDCAz = trackSelections.dcaZconst + trackSelections.dcaZpt * std::pow(track.pt(), trackSelections.dcaZexp);
    if(fabs(track.dcaZ()) > maxDCAz){
      return false;
    } 
    if (fillHists) {
      histos.fill(HIST("hTrackSelection"), 7);
    }

    //dcaZ VS pt after max dcaZ=0.2 selection 
    if (fillHists) {
      histos.fill(HIST("hdcaZVSpt"), track.dcaZ(), track.pt());
      histos.fill(HIST("hTrackSelection"), 8);
    }

    return true;
  }


  void process(CollisionCandidates::iterator const& collision, TrackCandidates const& tracks, aod::FT0s const& )
  {  
   
    if (!isEventAccepted(collision, true)) {
      return;
    }

    //Looping over tracks before selection
    for (int i = 0; i < tracks.size(); i++)  {
      auto track = tracks.iteratorAt(i);
      histos.fill(HIST("hNTracksVSetaBefCuts"), track.eta());
      histos.fill(HIST("hNTracksVSphiBefCuts"), track.phi()); 
      histos.fill(HIST("hdcaXYVSptBefCuts"), track.dcaXY(), track.pt()); 
    }
     //Collision selection level 
    if (!collision.has_foundFT0()) {
      return;
    }
    const auto& ft0 = collision.foundFT0();
    float ft0cAmp = ft0.sumAmpC();
    int multiplicity = collision.multNTracksPVeta1();
    histos.fill(HIST("hFT0CAmplitude"),ft0cAmp);
    histos.fill(HIST("hFT0CAmplitudeVsVz"), collision.posZ(), ft0cAmp);
    histos.fill(HIST("hFT0CAmplitudeVSMult"), ft0cAmp, multiplicity);


    //Looping over selected tracks
    int nSelTracks = 0;
    for (int i = 0; i < tracks.size(); i++)  {
      auto track = tracks.iteratorAt(i);
      if (!IsTrackSelected(track, true)) {
        continue;
      }
      nSelTracks++;
      histos.fill(HIST("hNTracksVSeta"), track.eta());
      histos.fill(HIST("hNTracksVSphi"), track.phi());
      histos.fill(HIST("hdcaXYVSpt"), track.dcaXY(), track.pt());
      histos.fill(HIST("hFT0CAmplitudeVSpT"), ft0cAmp, track.pt());
    }    

    //histograms for mult
    histos.fill(HIST("hMultSelectedTracks"), nSelTracks);
    histos.fill(HIST("hMultGlobalTracks"), collision.multNTracksGlobal());
    histos.fill(HIST("hMultPVeta1"), collision.multNTracksPVeta1());

  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityClassificationPbPb>(cfgc)
  };
}
