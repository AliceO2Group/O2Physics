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
/// \file taskSingleElectron.cxx
/// \brief task for electrons from heavy-flavour hadron decays
/// \author Jonghan Park (Jeonbuk National University)

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskSingleElectron {

  // Produces

  // Configurable
  Configurable<float> ptTrackMax{"ptTrackMax", 10., "max pt cut"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.5, "min pt cut"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "eta cut"};
  Configurable<int> ptcNCrossedRowMax{"ptcNCrossedRowMax", 70, "max of TPC n cluster crossed rows"};
  Configurable<float> tpcNClsFoundOverFindableMin{"tpcNClsFoundOverFindableMin", 0.8, "min # of TPC found/findable clusters"};
  Configurable<float> tpcChi2perNClMax{"tpcChi2perNClMax", 4., "min # of tpc chi2 per clusters"};
  Configurable<int> itsIBClsMin{"itsIBClsMin", 3, "min # of its clusters in IB"};
  Configurable<float> dcaxyMax{"dcaxyMax", 1., "max of track dca in xy"};
  Configurable<float> dcazMax{"dcazMax", 2., "max of track dca in z"};
  Configurable<float> tofNSigmaMax{"tofNSigmaMax", 3., "max of tof nsigma"};
  Configurable<float> tpcNSigmaMin{"tpcNSigmaMin", -1., "min of tpc nsigma"};
  Configurable<float> tpcNSigmaMax{"tpcNSigmaMax", 3., "max of tpc nsigma"};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsDeltaPhi{"nBinsDeltaPhi", 100, "N bins in DeltaPhi histo"};
  Configurable<int> nBinsDeltaEta{"nBinsDeltaEta", 100, "N bins in DeltaEta histo"};

  // SliceCache
  SliceCache cache;

  // using declarations
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using TracksEl = soa::Join<aod::Tracks, aod::TrackExtra, aod::TracksDCA, aod::pidTOFFullEl, aod::pidTPCFullEl>;

  // Filter
  Filter collZFilter = nabs(aod::collision::posZ) < 10.0f;

  // Partition

  // ConfigurableAxis
  ConfigurableAxis axisPtEl{"axisPtEl", {VARIABLE_WIDTH, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.75f, 2.0f, 2.25f, 2.5f, 2.75f, 3.f, 3.5f, 4.0f, 5.0f, 6.0f, 8.0f, 10.0f}, "electron pt bins"};

  // AxisSpec
  const AxisSpec axisEvt{4, 0., 4., "nEvents"};
  const AxisSpec axisNCont{100, 0., 100., "nCont"};
  const AxisSpec axisPosZ{600, -30., 30., "Z_{pos}"};
  const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
  const AxisSpec axisPt{nBinsPt, 0., 15., "p_{T}"};
  const AxisSpec axisNsig{800, -20., 20.};
  const AxisSpec axisTrackIp{4000, -0.2, 0.2, "dca"};
  const AxisSpec axisDeltaPhi{nBinsDeltaPhi, -PIHalf, PI + PIHalf, "#Delta#varphi"};
  const AxisSpec axisDeltaEta{nBinsDeltaEta, -1.0, +1.0, "#Delta#eta"};

  // Histogram registry
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // create histograms
    histos.add("hEventCounter", "hEventCounter", kTH1F, {axisEvt});
    histos.add("nEvents", "Number of events", kTH1F, {{1, 0., 1.}});
    histos.add("VtxZ", "VtxZ; cm; entries", kTH1F, {axisPosZ});
    histos.add("etaTrack", "etaTrack; #eta; entries", kTH1F, {axisEta});
    histos.add("ptTrack", "#it{p}_{T} distribution of selected tracks; #it{p}_{T} (GeV/#it{c}); entries", kTH1F, {axisPt});

    // QA plots for trigger track selection
    histos.add("tpcNClsTrack", "tpcNClsTrack", kTH1F, {{200, 0, 200}});
    histos.add("tpcFoundFindableTrack", "", kTH1F, {{10, 0, 1}});
    histos.add("tpcChi2Track", "", kTH1F, {{100, 0, 10}});
    histos.add("itsIBClsTrack", "", kTH1F, {{10, 0, 10}});
    histos.add("dcaXYTrack", "", kTH1F, {{600, -3, 3}});
    histos.add("dcaZTrack", "", kTH1F, {{600, -3, 3}});

    // QA plots for associated track selection
    histos.add("tpcNClsAsso", "tpcNClsAsso", kTH1F, {{200, 0, 200}});
    histos.add("tpcFoundFindableAsso", "", kTH1F, {{10, 0, 1}});
    histos.add("tpcChi2Asso", "", kTH1F, {{100, 0, 10}});
    histos.add("itsIBClsAsso", "", kTH1F, {{10, 0, 10}});
    histos.add("dcaXYAsso", "", kTH1F, {{600, -3, 3}});
    histos.add("dcaZAsso", "", kTH1F, {{600, -3, 3}});

    histos.add("correlationFunction", "correlationFunction", kTH1F, {axisDeltaPhi});
    histos.add("correlationFunction2d", "correlationFunction2d", kTH2F, {axisDeltaPhi, axisDeltaEta});

    // pid
    histos.add("tofNSigPt", "", kTH2F, {{axisPtEl}, {axisNsig}});
    histos.add("tofNSigPtQA", "", kTH2F, {{axisPtEl}, {axisNsig}});
    histos.add("tpcNSigPt", "", kTH2F, {{axisPtEl}, {axisNsig}});
    histos.add("tpcNSigPtAfterTofCut", "", kTH2F, {{axisPtEl}, {axisNsig}});
    histos.add("tpcNSigPtQA", "", kTH2F, {{axisPtEl}, {axisNsig}});

    // track impact parameter
    histos.add("dcaTrack", "", kTH2F, {{axisPtEl}, {axisTrackIp}});
  }

  template <typename TrackType>
  bool TrackSel(TrackType track)
  {
    if (track.pt() > ptTrackMax || track.pt() < ptTrackMin)
      return false;
    if (std::abs(track.eta()) > etaTrackMax)
      return false;

    int tpcNClsFound = track.tpcNClsCrossedRows();
    int tpcNClsFindable = track.tpcNClsFindable();
    float tpcFoundOverFindable = (tpcNClsFindable ? static_cast<float>(tpcNClsFound) / static_cast<float>(tpcNClsFindable) : 0);

    if (tpcNClsFound < ptcNCrossedRowMax)
      return false;
    if (tpcFoundOverFindable < tpcNClsFoundOverFindableMin)
      return false;
    if (track.tpcChi2NCl() > tpcChi2perNClMax)
      return false;

    if (!(track.itsNClsInnerBarrel() == itsIBClsMin))
      return false;

    if (std::abs(track.dcaXY()) > dcaxyMax)
      return false;
    if (std::abs(track.dcaZ()) > dcazMax)
      return false;

    histos.fill(HIST("etaTrack"), track.eta());
    histos.fill(HIST("ptTrack"), track.pt());

    histos.fill(HIST("tpcNClsTrack"), tpcNClsFound);
    histos.fill(HIST("tpcFoundFindableTrack"), tpcFoundOverFindable);
    histos.fill(HIST("tpcChi2Track"), track.tpcChi2NCl());
    histos.fill(HIST("itsIBClsTrack"), track.itsNClsInnerBarrel());
    histos.fill(HIST("dcaXYTrack"), track.dcaXY());
    histos.fill(HIST("dcaZTrack"), track.dcaZ());

    return true;
  }

  template <typename TrackType>
  bool AssoTrackSel(TrackType track)
  {

    if (std::abs(track.eta()) > etaTrackMax)
      return false;
    if (track.pt() < 4.0f || track.pt() > 6.0f)
      return false;

    int tpcNClsFound = track.tpcNClsCrossedRows();
    int tpcNClsFindable = track.tpcNClsFindable();
    float tpcFoundOverFindable = (tpcNClsFindable ? static_cast<float>(tpcNClsFound) / static_cast<float>(tpcNClsFindable) : 0);

    if (tpcNClsFound < 60)
      return false;
    if (track.tpcChi2NCl() > tpcChi2perNClMax)
      return false;

    if (std::abs(track.dcaXY()) > dcaxyMax)
      return false;
    if (std::abs(track.dcaZ()) > dcazMax)
      return false;

    histos.fill(HIST("tpcNClsAsso"), tpcNClsFound);
    histos.fill(HIST("tpcFoundFindableAsso"), tpcFoundOverFindable);
    histos.fill(HIST("tpcChi2Asso"), track.tpcChi2NCl());
    histos.fill(HIST("itsIBClsAsso"), track.itsNClsInnerBarrel());
    histos.fill(HIST("dcaXYAsso"), track.dcaXY());
    histos.fill(HIST("dcaZAsso"), track.dcaZ());

    return true;
  }

  double ComputeDeltaPhi(double phi1, double phi2)
  {
    double deltaPhi = phi1 - phi2;
    if (deltaPhi < -PIHalf)
      deltaPhi += 2. * PI;
    if (deltaPhi > 3. * PIHalf)
      deltaPhi -= 2. * PI;
    return deltaPhi;
  }

  void process(soa::Filtered<MyCollisions>::iterator const& collision, TracksEl const& tracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);

    if (!collision.sel8())
      return;
    histos.fill(HIST("hEventCounter"), 1.5);

    if (collision.numContrib() < 2)
      return;
    histos.fill(HIST("hEventCounter"), 2.5);

    histos.fill(HIST("VtxZ"), collision.posZ());
    histos.fill(HIST("hEventCounter"), 3.5);
    histos.fill(HIST("nEvents"), 0.5);

    for (auto& track : tracks) {

      if (!TrackSel(track))
        continue;

      histos.fill(HIST("tofNSigPt"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("tpcNSigPt"), track.pt(), track.tpcNSigmaEl());

      if (std::abs(track.tofNSigmaEl()) > tofNSigmaMax)
        continue;
      histos.fill(HIST("tofNSigPtQA"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("tpcNSigPtAfterTofCut"), track.pt(), track.tpcNSigmaEl());

      if (track.tpcNSigmaEl() < tpcNSigmaMin || track.tpcNSigmaEl() > tpcNSigmaMax)
        continue;
      histos.fill(HIST("tpcNSigPtQA"), track.pt(), track.tpcNSigmaEl());

      histos.fill(HIST("dcaTrack"), track.pt(), track.dcaXY());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskSingleElectron>(cfgc)};
}
