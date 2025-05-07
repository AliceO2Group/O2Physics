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
/// \file pmdqa.cxx
///
/// \brief QA task to check PMD info on Run 2 converted data
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since February 19, 2025

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "CCDB/BasicCCDBManager.h"
#include "TH1F.h"
#include "TH2F.h"
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::aod::run2;
using namespace o2::framework;
using namespace o2::aod::evsel;
using namespace o2::framework::expressions;

namespace o2::aod
{
  namespace pmdtrack
  {
    DECLARE_SOA_INDEX_COLUMN(Collision, collision);
    DECLARE_SOA_ARRAY_INDEX_COLUMN(Collision, collisions);
    DECLARE_SOA_INDEX_COLUMN(BC, bc);
    DECLARE_SOA_SLICE_INDEX_COLUMN(Pmd,pmd);
  } // namespace pmdtrack
  
  DECLARE_SOA_INDEX_TABLE_USER(PMDTracksIndex, BCs, "PMDTRKIDX", pmdtrack::CollisionId, pmdtrack::BCId, pmdtrack::PmdIdSlice);
}

struct BuiltPmdIndex {
  // build the index table PMDTracksIndex
  Builds<aod::PMDTracksIndex> idx;
  void init(InitContext const&){};
};

struct PmdQa {

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis axisEventBin{"axisEventBin", {4, 0.5, 4.5}, ""};
  ConfigurableAxis axisVtxZBin{"axisVtxZBin", {40, -20, 20}, ""};
  ConfigurableAxis axisNPMDtracksBin{"axisNPMDtracksBin", {500, 0, 500}, "Number of pmdtracks"};
  ConfigurableAxis axisClsxyBin{"axisClsxyBin", {200,-100,100}, ""};
  ConfigurableAxis axisAdcBin{"axisAdcBin", {200,0,2000}, ""};
  ConfigurableAxis axisEtaBin{"axisEtaBin", {10,2.1,4.1}, ""};
  ConfigurableAxis axisNcellBin{"axisNcellBin", {50,-0.5,49.5}, ""};
  Configurable<int> fMipCut{"fMipCut", 432, "fMipCut"};
    
  void init(InitContext&)
  {

    AxisSpec axisEvent = {axisEventBin, "Event", "EventAxis"};
    AxisSpec axisVtxZ = {axisVtxZBin, "VtxZ", "VtxZAxis"};
    AxisSpec axisNPMDtracks = {axisNPMDtracksBin, "NPMDtracks", "NPMDtracksAxis"};
    AxisSpec axisClsxy = {axisClsxyBin, "Clsxy", "ClsxyAxis"};
    AxisSpec axisAdc = {axisAdcBin, "Adc", "AdcAxis"};
    AxisSpec axisEta = {axisEtaBin, "Eta", "EtaAxis"};
    AxisSpec axisNcell = {axisNcellBin, "Ncell", "NcellAxis"};
    
    histos.add("hEventHist", "hEventHist", kTH1F, {axisEvent});
    histos.add("hVtxZHist", "hVtxZHist", kTH1F, {axisVtxZ});
    histos.add("hNPMDtracks", "Number of pmdtracks", kTH1F, {axisNPMDtracks});
    histos.add("hClusXY", "hClusXY", kTH2F, {axisClsxy, axisClsxy});
    histos.add("hClusAdc", "hClusAdc", kTH1F, {axisAdc});
    histos.add("hetacls", "hetacls", kTH1F, {axisEta});
    histos.add("hclsncell", "hclsncell", kTH1F, {axisNcell});
  }

  using coltable = soa::Join<aod::Collisions, aod::PMDTracksIndex>;
  using colevsel = soa::Join<coltable, aod::EvSels>;
    
  void process(colevsel::iterator const& collision, aod::Pmds const&)
  {
    histos.fill(HIST("hEventHist"), 1);
    if (collision.sel7()) {
      return;
    }
    histos.fill(HIST("hEventHist"), 2);
    if (std::abs(collision.posZ()) >= 10.) {
      return;
    }
    histos.fill(HIST("hEventHist"), 3);
    histos.fill(HIST("hVtxZHist"), collision.posZ());
    
    if (collision.has_pmd()) {
      histos.fill(HIST("hEventHist"), 4);
      auto tracks = collision.pmd();
      histos.fill(HIST("hNPMDtracks"), tracks.size());
      for (const auto& track : tracks) {
	if (track.pmddet() == 1) {
	  return;
	}
	if(track.pmdclsz() == 0) {
	  return;
	}
	if (!track.pmdmodule()) {
	  return;
	}
	histos.fill(HIST("hClusXY"), track.pmdclsx(), track.pmdclsy());
	histos.fill(HIST("hClusAdc"), track.pmdclsadc());
	float rdist = TMath::Sqrt(track.pmdclsx()*track.pmdclsx() + track.pmdclsy()*track.pmdclsy());
	float theta = TMath::ATan2(rdist,track.pmdclsz());
	float etacls  = -TMath::Log(TMath::Tan(0.5*theta));
	if (track.pmdclsadc() > fMipCut && track.pmdncell() > 2) {
	  if(etacls>2.3 && etacls<3.9){
	    histos.fill(HIST("hetacls"), etacls);
	    histos.fill(HIST("hclsncell"), track.pmdncell());
	  }
	}
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PmdQa>(cfgc),
    adaptAnalysisTask<BuiltPmdIndex>(cfgc),
  };
}
