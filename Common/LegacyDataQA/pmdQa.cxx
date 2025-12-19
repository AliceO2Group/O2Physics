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
/// \file pmdQa.cxx
///
/// \brief QA task to check PMD info on Run 2 converted data
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since May 17, 2025

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PmdTable.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::aod::run2;
using namespace o2::framework;
using namespace o2::aod::evsel;
using namespace o2::framework::expressions;

struct BuiltPmdIndex {
  // build the index table PMDTracksIndex
  Builds<aod::PMDTracksIndex> idx;
  void init(InitContext const&) {}
};

struct PmdQa {

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis axisEventBin{"axisEventBin", {4, 0.5, 4.5}, ""};
  ConfigurableAxis axisVtxZBin{"axisVtxZBin", {40, -20, 20}, ""};
  ConfigurableAxis axisNPMDtracksBin{"axisNPMDtracksBin", {500, 0, 500}, "Number of pmdtracks"};
  ConfigurableAxis axisClsxyBin{"axisClsxyBin", {200, -100, 100}, ""};
  ConfigurableAxis axisAdcBin{"axisAdcBin", {200, 0, 2000}, ""};
  ConfigurableAxis axisEtaBin{"axisEtaBin", {10, 2.1, 4.1}, ""};
  ConfigurableAxis axisNcellBin{"axisNcellBin", {50, -0.5, 49.5}, ""};
  Configurable<int> fMipCut{"fMipCut", 432, "fMipCut"};
  Configurable<int> fNcellCut{"fNcellCut", 2, "fNcellCut"};
  Configurable<float> fEtalow{"fEtalow", 2.3, "fEtalow"};
  Configurable<float> fEtahigh{"fEtahigh", 3.9, "fEtahigh"};
  Configurable<float> fVtxCut{"fVtxCut", 10.0, "fVtxCut"};

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

  using ColTable = soa::Join<aod::Collisions, aod::PMDTracksIndex>;
  using ColevSel = soa::Join<ColTable, aod::EvSels>;

  void process(ColevSel::iterator const& collision, aod::Pmds const&)
  {
    histos.fill(HIST("hEventHist"), 1);
    if (collision.sel7()) {
      return;
    }
    histos.fill(HIST("hEventHist"), 2);
    if (std::abs(collision.posZ()) >= fVtxCut) {
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
        if (track.pmdclsz() == 0) {
          return;
        }
        if (!track.pmdmodule()) {
          return;
        }
        histos.fill(HIST("hClusXY"), track.pmdclsx(), track.pmdclsy());
        histos.fill(HIST("hClusAdc"), track.pmdclsadc());
        float rdist = std::sqrt(track.pmdclsx() * track.pmdclsx() + track.pmdclsy() * track.pmdclsy());
        float theta = std::atan2(rdist, track.pmdclsz());
        float etacls = -std::log(std::tan(0.5 * theta));
        if (track.pmdclsadc() > fMipCut && track.pmdncell() > fNcellCut) {
          if (etacls > fEtalow && etacls < fEtahigh) {
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
