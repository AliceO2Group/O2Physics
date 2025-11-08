// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamPairEfficiency.cxx
/// \brief Task to produce dNdeta for pair triggered events
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// heaviyly inspiered by phik0shortanalysis.cxx

#include "PWGCF/FemtoDream/Core/femtoDreamMath.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "TPDGCode.h"

#include "fairlogger/Logger.h"

#include <algorithm>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::aod::track;

struct FemtoDreamPairEfficiency {

  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;
  using RecoCollisions = soa::Join<Collisions, aod::McCollisionLabels>;
  using GenCollisions = soa::Join<aod::McCollisions, aod::MultMCExtras, aod::McCentFT0Ms>;

  // Defining the type of the tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                               pidTPCFullEl, pidTPCFullPi, pidTPCFullKa,
                               pidTPCFullPr, pidTPCFullDe, pidTPCFullTr, pidTPCFullHe,
                               pidTOFFullEl, pidTOFFullPi, pidTOFFullKa,
                               pidTOFFullPr, pidTOFFullDe, pidTOFFullTr, pidTOFFullHe>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;
  // filter for tracks
  static constexpr TrackSelectionFlags::flagtype TrackSelectionITS = TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF | TrackSelectionFlags::kITSHits;
  static constexpr TrackSelectionFlags::flagtype TrackSelectionTPC = TrackSelectionFlags::kTPCNCls | TrackSelectionFlags::kTPCCrossedRowsOverNCls | TrackSelectionFlags::kTPCChi2NDF;
  static constexpr TrackSelectionFlags::flagtype TrackSelectionDCA = TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

  Filter trackFilter = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                       ncheckbit(aod::track::trackCutFlag, TrackSelectionITS) &&
                       ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC), ncheckbit(aod::track::trackCutFlag, TrackSelectionTPC), true) &&
                       ncheckbit(aod::track::trackCutFlag, TrackSelectionDCA) &&
                       ncheckbit(aod::track::trackCutFlag, TrackSelectionFlags::kInAcceptanceTracks);

  using FilteredFullTracks = soa::Filtered<FullTracks>;
  using FilteredFullMCTracks = soa::Filtered<FullMCTracks>;

  Service<o2::framework::O2DatabasePDG> pdgDB;

  SliceCache cache;
  Preslice<McParticles> perColMc = mcparticle::mcCollisionId;

  struct : ConfigurableGroup {
    std::string prefix = std::string("Mode");
    Configurable<bool> countPairs{"countPairs", true, "Count Pairs"};
    Configurable<bool> countTriplets{"countTriplets", false, "Count Triplets"};
  } mode;

  // Event Selections
  struct : ConfigurableGroup {
    std::string prefix = std::string("EventSelection");
    Configurable<float> zvtxAbsMax{"zvtxAbsMax", 10.f, "|z-Vertex| max"};
    Configurable<bool> offlineCheck{"offlineCheck", true, "Check for Sel8"};
    Configurable<float> etaAbsMax{"etaAbsMax", 0.8f, "Common eta cut for particles in dNdEta distribution"};
    Configurable<float> kstarMax{"kstarMax", 999.f, "Cut on kstar"};
    Configurable<float> q3Max{"q3Max", 999.f, "Cut on Q3"};
  } eventSelection;

  struct : ConfigurableGroup {
    std::string prefix = std::string("EventBinning");
    ConfigurableAxis zvtxBinning{"zvtxBinning", {240, -12, 12}, "z-vertex Binning"};
    ConfigurableAxis centBinning{"centBinning", {120, 0, 120}, "Centrality Binning"};
    ConfigurableAxis multBinning{"multBinning", {300, 0, 300}, "Multiplicity Binning"};
    ConfigurableAxis dNdetaBinning{"dNdetaBinning", {200, -1, 1}, "dNdeta Binning"};
    ConfigurableAxis kstarBinning{"kstarBinning", {200, 0, 2}, "dNdeta Binning"};
    ConfigurableAxis pairBinning{"pairBinning", {50, 0, 50}, "pair Binning"};
  } eventBinning;

  // Track 1 selections
  struct : ConfigurableGroup {
    std::string prefix = std::string("SelectionTrack1");
    Configurable<int> sign{"sign", 1, "Sign of charge"};
    Configurable<float> ptMin{"ptMin", 0.0f, "pt min"};
    Configurable<float> ptMax{"ptMax", 3.0f, "pt max"};
    Configurable<float> etaAbsMax{"etaAbsMax", 0.8f, "|eta| max"};
    Configurable<float> dcazAbsMax{"dcazAbsMax", 0.1f, "|dca_z| max"};
    Configurable<bool> useDcaxyPtDepCut{"useDcaxyPtDepCut", true, "|dca_z| max"};
    Configurable<float> tpcClusterMin{"tpcClusterMin", 80.f, "TPC clusters min"};
    Configurable<float> tpcCrossedOverClusterMin{"tpcCrossedOverClusterMin", 0.83f, "TPC clusters/TPC crossed rows min"};
    Configurable<float> tpcCrossedMin{"tpcCrossedMin", 70.f, "TPC crossed rows min"};
    Configurable<float> tpcSharedMax{"tpcSharedMax", 160.f, "TPC shared clusters max"};
    Configurable<float> itsClusterMin{"itsClusterMin", 0.f, "ITS clusters min"};
    Configurable<float> itsIbClusterMin{"itsIbClusterMin", 0.f, "ITS inner barrle min"};
    Configurable<int> pdgCode{"pdgCode", 2212, "PDG code"};
    Configurable<float> pidThreshold{"pidThreshold", 0.75f, "Momentum threshold for PID"};
    Configurable<float> itsNsigmaMax{"itsNsigmaMax", 99.f, "its nsigma max"};
    Configurable<float> tpcNsigmaMax{"tpcNsigmaMax", 3.f, "TPC nsigma max"};
    Configurable<float> tpctofNsigmaMax{"tpctofNsigmaMax", 3.f, "TPCTOC nsigma max"};
  } trackCuts1;

  // Track 2 selections
  struct : ConfigurableGroup {
    std::string prefix = std::string("SelectionTrack2");
    Configurable<int> sign{"sign", 1, "Sign of charge"};
    Configurable<float> ptMin{"ptMin", 0.0f, "pt min"};
    Configurable<float> ptMax{"ptMax", 3.0f, "pt max"};
    Configurable<float> etaAbsMax{"etaAbsMax", 0.8f, "|eta| max"};
    Configurable<float> dcazAbsMax{"dcazAbsMax", 0.1f, "|dca_z| max"};
    Configurable<bool> useDcaxyPtDepCut{"useDcaxyPtDepCut", true, "|dca_z| max"};
    Configurable<float> tpcClusterMin{"tpcClusterMin", 80.f, "TPC clusters min"};
    Configurable<float> tpcCrossedOverClusterMin{"tpcCrossedOverClusterMin", 0.83f, "TPC clusters/TPC crossed rows min"};
    Configurable<float> tpcCrossedMin{"tpcCrossedMin", 70.f, "TPC crossed rows min"};
    Configurable<float> tpcSharedMax{"tpcSharedMax", 160.f, "TPC shared clusters max"};
    Configurable<float> itsClusterMin{"itsClusterMin", 0.f, "ITS clusters min"};
    Configurable<float> itsIbClusterMin{"itsIbClusterMin", 0.f, "ITS inner barrle min"};
    Configurable<int> pdgCode{"pdgCode", 2212, "PDG code"};
    Configurable<float> pidThreshold{"pidThreshold", 0.75f, "Momentum threshold for PID"};
    Configurable<float> itsNsigmaMax{"itsNsigmaMax", 99.f, "its nsigma max"};
    Configurable<float> tpcNsigmaMax{"tpcNsigmaMax", 3.f, "TPC nsigma max"};
    Configurable<float> tpctofNsigmaMax{"tpctofNsigmaMax", 3.f, "TPCTOC nsigma max"};
  } trackCuts2;

  struct : ConfigurableGroup {
    std::string prefix = std::string("SelectionTrack3");
    Configurable<int> sign{"sign", 1, "Sign of charge"};
    Configurable<float> ptMin{"ptMin", 0.0f, "pt min"};
    Configurable<float> ptMax{"ptMax", 3.0f, "pt max"};
    Configurable<float> etaAbsMax{"etaAbsMax", 0.8f, "|eta| max"};
    Configurable<float> dcazAbsMax{"dcazAbsMax", 0.1f, "|dca_z| max"};
    Configurable<bool> useDcaxyPtDepCut{"useDcaxyPtDepCut", true, "|dca_z| max"};
    Configurable<float> tpcClusterMin{"tpcClusterMin", 80.f, "TPC clusters min"};
    Configurable<float> tpcCrossedOverClusterMin{"tpcCrossedOverClusterMin", 0.83f, "TPC clusters/TPC crossed rows min"};
    Configurable<float> tpcCrossedMin{"tpcCrossedMin", 70.f, "TPC crossed rows min"};
    Configurable<float> tpcSharedMax{"tpcSharedMax", 160.f, "TPC shared clusters max"};
    Configurable<float> itsClusterMin{"itsClusterMin", 0.f, "ITS clusters min"};
    Configurable<float> itsIbClusterMin{"itsIbClusterMin", 0.f, "ITS inner barrle min"};
    Configurable<int> pdgCode{"pdgCode", 2212, "PDG code"};
    Configurable<float> pidThreshold{"pidThreshold", 0.75f, "Momentum threshold for PID"};
    Configurable<float> itsNsigmaMax{"itsNsigmaMax", 99.f, "its nsigma max"};
    Configurable<float> tpcNsigmaMax{"tpcNsigmaMax", 3.f, "TPC nsigma max"};
    Configurable<float> tpctofNsigmaMax{"tpctofNsigmaMax", 3.f, "TPCTOC nsigma max"};
  } trackCuts3;

  HistogramRegistry dataEventHist{"dataEventHist", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry mcEventHist{"mcEventHist", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  HistogramRegistry dataTrackHist{"dataTrackHist", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  void init(InitContext&)
  {

    if (mode.countPairs.value && mode.countTriplets.value) {
      LOG(fatal) << "We cannot count pairs and triplets at the same time. Turn one of the off.";
    }

    dataEventHist.add("hDataEventSelection", "hDataEventSelection", kTH1F, {{6, -0.5f, 5.5f}});
    dataEventHist.get<TH1>(HIST("hDataEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    dataEventHist.get<TH1>(HIST("hDataEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    dataEventHist.get<TH1>(HIST("hDataEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    dataEventHist.get<TH1>(HIST("hDataEventSelection"))->GetXaxis()->SetBinLabel(4, "INEL>0 cut");
    dataEventHist.get<TH1>(HIST("hDataEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a pair");
    dataEventHist.get<TH1>(HIST("hDataEventSelection"))->GetXaxis()->SetBinLabel(6, "With at least a lowkstar pair");

    // Event information
    dataEventHist.add("hDataVertexZ", "hVertexZ", kTH1F, {eventBinning.zvtxBinning});
    dataEventHist.add("hDataMultiplicity", "hMultiplicity", kTH1F, {eventBinning.multBinning});
    dataEventHist.add("hDataCentrality", "Centrality", kTH1F, {eventBinning.centBinning});

    // Eta distribution for dN/deta values estimation in Data
    dataEventHist.add("hDataCentralityVsEtaDistribtion", "Eta vs multiplicity in Data", kTH2F, {eventBinning.centBinning, eventBinning.dNdetaBinning});

    dataTrackHist.add("Track1/pt", "Track 1 pt", kTH1F, {{600, 0, 6}});
    dataTrackHist.add("Track1/eta", "Track 1 eta", kTH1F, {{200, -1, 1}});
    dataTrackHist.add("Track1/phi", "Track 1 phi", kTH1F, {{720, 0, o2::constants::math::TwoPI}});
    dataTrackHist.add("Track1/tpcCluster", "Track 1 cluster", kTH1F, {{160, 0, 160}});
    dataTrackHist.add("Track1/tpcCrossed", "Track 1 crossed rows", kTH1F, {{160, 0, 160}});
    dataTrackHist.add("Track1/tpcShared", "Track 1 cluster shared", kTH1F, {{160, 0, 160}});
    dataTrackHist.add("Track1/itsCluster", "Track 1 its cluster", kTH1F, {{7, 0, 7}});
    dataTrackHist.add("Track1/itsIbCluster", "Track 1 its cluster inner barrel", kTH1F, {{3, 0, 3}});
    dataTrackHist.add("Track1/itsNsigma", "Track 1 nsigma its", kTH2F, {{600, 0, 6}, {1000, -5, 5}});
    dataTrackHist.add("Track1/tpcNsigma", "Track 1 nsigma tpc", kTH2F, {{600, 0, 6}, {1000, -5, 5}});
    dataTrackHist.add("Track1/tpctofNsigma", "Track 1 nsigma tpctof", kTH2F, {{600, 0, 6}, {500, 0, 5}});

    dataTrackHist.add("Track2/pt", "Track 2 pt", kTH1F, {{600, 0, 6}});
    dataTrackHist.add("Track2/eta", "Track 2 eta", kTH1F, {{200, -1, 1}});
    dataTrackHist.add("Track2/phi", "Track 2 phi", kTH1F, {{720, 0, o2::constants::math::TwoPI}});
    dataTrackHist.add("Track2/tpcCluster", "Track 2 cluster", kTH1F, {{160, 0, 160}});
    dataTrackHist.add("Track2/tpcCrossed", "Track 2 crossed rows", kTH1F, {{160, 0, 160}});
    dataTrackHist.add("Track2/tpcShared", "Track 2 cluster shared", kTH1F, {{160, 0, 160}});
    dataTrackHist.add("Track2/itsCluster", "Track 2 its cluster", kTH1F, {{0, 0, 7}});
    dataTrackHist.add("Track2/itsIbCluster", "Track 2 its cluster inner barrel", kTH1F, {{3, 0, 3}});
    dataTrackHist.add("Track2/itsNsigma", "Track 2 nsigma its", kTH2F, {{600, 0, 6}, {1000, -5, 5}});
    dataTrackHist.add("Track2/tpcNsigma", "Track 2 nsigma tpc", kTH2F, {{600, 0, 6}, {1000, -5, 5}});
    dataTrackHist.add("Track2/tpctofNsigma", "Track 2 nsigma tpctof", kTH2F, {{600, 0, 6}, {500, 0, 5}});

    if (mode.countTriplets) {
      dataTrackHist.add("Track3/pt", "Track 3 pt", kTH1F, {{600, 0, 6}});
      dataTrackHist.add("Track3/eta", "Track 3 eta", kTH1F, {{200, -1, 1}});
      dataTrackHist.add("Track3/phi", "Track 3 phi", kTH1F, {{720, 0, o2::constants::math::TwoPI}});
      dataTrackHist.add("Track3/tpcCluster", "Track 3 cluster", kTH1F, {{160, 0, 160}});
      dataTrackHist.add("Track3/tpcCrossed", "Track 3 crossed rows", kTH1F, {{160, 0, 160}});
      dataTrackHist.add("Track3/tpcShared", "Track 3 cluster shared", kTH1F, {{160, 0, 160}});
      dataTrackHist.add("Track3/itsCluster", "Track 3 its cluster", kTH1F, {{0, 0, 7}});
      dataTrackHist.add("Track3/itsIbCluster", "Track 3 its cluster inner barrel", kTH1F, {{3, 0, 3}});
      dataTrackHist.add("Track3/itsNsigma", "Track 3 nsigma its", kTH2F, {{600, 0, 6}, {1000, -5, 5}});
      dataTrackHist.add("Track3/tpcNsigma", "Track 3 nsigma tpc", kTH2F, {{600, 0, 6}, {1000, -5, 5}});
      dataTrackHist.add("Track3/tpctofNsigma", "Track 3 nsigma tpctof", kTH2F, {{600, 0, 6}, {500, 0, 5}});
    }

    mcEventHist.add("hRecoEventSelection", "hRecoEventSelection", kTH1F, {{7, -0.5f, 6.5f}});
    mcEventHist.get<TH1>(HIST("hRecoEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    mcEventHist.get<TH1>(HIST("hRecoEventSelection"))->GetXaxis()->SetBinLabel(2, "Sel8 cut");
    mcEventHist.get<TH1>(HIST("hRecoEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    mcEventHist.get<TH1>(HIST("hRecoEventSelection"))->GetXaxis()->SetBinLabel(4, "INEL>0 cut");
    mcEventHist.get<TH1>(HIST("hRecoEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a gen coll");
    mcEventHist.get<TH1>(HIST("hRecoEventSelection"))->GetXaxis()->SetBinLabel(6, "With at least a pair/triplet");
    mcEventHist.get<TH1>(HIST("hRecoEventSelection"))->GetXaxis()->SetBinLabel(7, "With at least a lowkstar pair/lowq3 triplet");

    mcEventHist.add("hGenEventSelection", "hGenEventSelection", kTH1F, {{6, -0.5f, 5.5f}});
    mcEventHist.get<TH1>(HIST("hGenEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    mcEventHist.get<TH1>(HIST("hGenEventSelection"))->GetXaxis()->SetBinLabel(2, "posZ cut");
    mcEventHist.get<TH1>(HIST("hGenEventSelection"))->GetXaxis()->SetBinLabel(3, "INEL>0 cut");
    mcEventHist.get<TH1>(HIST("hGenEventSelection"))->GetXaxis()->SetBinLabel(4, "With at least a pair/triplet");
    mcEventHist.get<TH1>(HIST("hGenEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a lowkstar pair/lowq3 triplet");
    mcEventHist.get<TH1>(HIST("hGenEventSelection"))->GetXaxis()->SetBinLabel(6, "With at least a reco coll");

    // MC Event information for Rec and Gen
    mcEventHist.add("hRecoMcVertexZ", "McVertexZ of reconstructed collision", kTH1F, {eventBinning.zvtxBinning});
    mcEventHist.add("hRecoMcMultiplicity", "McMultiplicity of reconstructed collision", kTH1F, {eventBinning.multBinning});
    mcEventHist.add("hRecoMcCentrality", "McCentrality of reconstructed collision", kTH1F, {eventBinning.centBinning});
    mcEventHist.add("hRecoMcKstar", "Mck* in reconstructed collision", kTH1F, {eventBinning.kstarBinning});
    mcEventHist.add("hRecoMcNumberOfPair", "McNumber of pairs in reconstructed collision", kTH1F, {eventBinning.pairBinning});

    mcEventHist.add("hRecoMcCentralityVsEtaDistribution", "McCentrality vs Eta in Reco", kTH2F, {eventBinning.centBinning, eventBinning.dNdetaBinning});
    mcEventHist.add("hRecoMcCentralityVsMcEtaDistribution", "McCentrality vs McEta in Reco (Check)", kTH2F, {eventBinning.centBinning, eventBinning.dNdetaBinning});

    mcEventHist.add("hGenMcVertexZ", "McVertex of generated collision", kTH1F, {eventBinning.zvtxBinning});
    mcEventHist.add("hGenMcMultiplicity", "McMultiplicity fo generated collision", kTH1F, {eventBinning.multBinning});
    mcEventHist.add("hGenMcCentrality", "McCentrality of generated collision", kTH1F, {eventBinning.centBinning});
    mcEventHist.add("hGenMcKstar", "Mckk* of generated collision", kTH1F, {eventBinning.kstarBinning});
    mcEventHist.add("hGenMcNumberOfPair", "McNumber of pairs in generated collision", kTH1F, {eventBinning.pairBinning});

    mcEventHist.add("hGenMcCentralityWithAllAssoc", "McCentrality in all associated collision", kTH1F, {eventBinning.centBinning});
    mcEventHist.add("hGenMcCentralityWithAssoc", "McCentrality in associated collision", kTH1F, {eventBinning.centBinning});
    mcEventHist.add("hGenMcCentralityVsMcEtaDistribution", "McCentrality vs McdNdEta in generated collision", kTH2F, {eventBinning.centBinning, eventBinning.dNdetaBinning});
    mcEventHist.add("hGenMcCentralityVsMcEtaDistributionWithAllAssoc", "McCentrality vs McdNdEta in all associated collisions", kTH2F, {eventBinning.centBinning, eventBinning.dNdetaBinning});
    mcEventHist.add("hGenMcCentralityVsMcEtaDistributionWithAssoc", "McCentrality vs McdNdEta in associated collisions", kTH2F, {eventBinning.centBinning, eventBinning.dNdetaBinning});
  }

  int pairsFound = 0;
  std::vector<float> vecKstar{};

  int tripletsFound = 0;
  std::vector<float> vecq3{};

  template <typename T1, typename T2>
  bool checkRecoTrackSelections(const T1& track, const T2& sel)
  {
    if (track.sign() != sel.sign.value) {
      return false;
    }
    if (track.pt() < sel.ptMin.value || track.pt() > sel.ptMax.value) {
      return false;
    }
    if (std::fabs(track.eta()) > sel.etaAbsMax.value) {
      return false;
    }
    if (sel.useDcaxyPtDepCut.value && std::fabs(track.dcaXY()) > (0.0105 + (0.035 / std::powf(track.pt(), 1.1f)))) {
      return false;
    }
    if (std::fabs(track.dcaZ()) > sel.dcazAbsMax.value) {
      return false;
    }
    if (track.tpcNClsFound() < sel.tpcClusterMin.value) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < sel.tpcCrossedMin.value) {
      return false;
    }
    if (track.tpcNClsShared() > sel.tpcSharedMax.value) {
      return false;
    }
    if (track.itsNCls() < sel.itsClusterMin.value) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < sel.itsIbClusterMin.value) {
      return false;
    }
    return true;
  }

  template <typename T>
  std::array<float, 3> getNSigmaValues(const T& track, int pdgCode)
  {
    std::array<float, 3> nsigma;
    switch (std::abs(pdgCode)) {
      case kPiPlus:
        nsigma[0] = track.itsNSigmaPi();
        nsigma[1] = track.tpcNSigmaPi();
        nsigma[2] = track.tofNSigmaPi();
        break;
      case kKPlus:
        nsigma[0] = track.itsNSigmaKa();
        nsigma[1] = track.tpcNSigmaKa();
        nsigma[2] = track.tofNSigmaKa();
        break;
      case kProton:
        nsigma[0] = track.itsNSigmaPr();
        nsigma[1] = track.tpcNSigmaPr();
        nsigma[2] = track.tofNSigmaPr();
        break;
      case constants::physics::kDeuteron:
        nsigma[0] = track.itsNSigmaDe();
        nsigma[1] = track.tpcNSigmaDe();
        nsigma[2] = track.tofNSigmaDe();
        break;
      default:
        LOG(fatal) << "PDG code " << pdgCode << " is not supported";
    }
    return nsigma;
  }

  template <typename T1, typename T2>
  bool checkRecoTrackPidSelections(const T1& track, const T2& sel)
  {
    auto nsigma = getNSigmaValues(track, sel.pdgCode.value);
    if (track.p() < sel.pidThreshold.value) {
      if (std::fabs(nsigma[1]) > sel.tpcNsigmaMax.value || std::fabs(nsigma[0]) > sel.itsNsigmaMax.value) {
        return false;
      }
    } else {
      if (std::hypot(nsigma[1], nsigma[2]) > sel.tpctofNsigmaMax.value) {
        return false;
      }
    }
    return true;
  }

  template <bool isMc, typename T>
  bool checkEventSelections(const T& col, bool qa)
  {
    if constexpr (!isMc) {
      if (qa) {
        dataEventHist.fill(HIST("hDataEventSelection"), 0); // all collisins
      }
      if (!col.sel8()) {
        return false;
      }
      if (qa) {
        dataEventHist.fill(HIST("hDataEventSelection"), 1); // sel8 collisions
      }
      if (std::fabs(col.posZ()) > eventSelection.zvtxAbsMax.value) {
        return false;
      }
      if (qa) {
        dataEventHist.fill(HIST("hDataEventSelection"), 2); // vertex-Z selected
      }
      if (!col.isInelGt0()) {
        return false;
      }
      if (qa) {
        dataEventHist.fill(HIST("hDataEventSelection"), 3); // INEL>0 collisions
      }
    } else {
      if (qa) {
        mcEventHist.fill(HIST("hRecoEventSelection"), 0); // all collisions
      }
      if (!col.sel8()) {
        return false;
      }
      if (qa) {
        mcEventHist.fill(HIST("hRecoEventSelection"), 1); // sel8 collisions
      }
      if (std::fabs(col.posZ()) > eventSelection.zvtxAbsMax.value) {
        return false;
      }
      if (qa) {
        mcEventHist.fill(HIST("hRecoEventSelection"), 2); // vertex-Z selected
      }
      if (!col.isInelGt0()) {
        return false;
      }
      if (qa) {
        mcEventHist.fill(HIST("hRecoEventSelection"), 3); // INEL>0 collisions
      }
    }
    return true;
  }

  template <typename T>
  int countRecoPairs(const T& tracks)
  {
    int pairs = 0;
    for (auto track1 = tracks.begin(); track1 != tracks.end(); track1++) {
      if (!checkRecoTrackSelections(track1, trackCuts1) || !checkRecoTrackPidSelections(track1, trackCuts1)) {
        continue;
      }
      for (auto track2 = track1 + 1; track2 != tracks.end(); track2++) {
        if (!checkRecoTrackSelections(track2, trackCuts2) || !checkRecoTrackPidSelections(track2, trackCuts2)) {
          continue;
        }

        auto nsigma1 = getNSigmaValues(track1, trackCuts1.pdgCode.value);
        dataTrackHist.fill(HIST("Track1/pt"), track1.pt());
        dataTrackHist.fill(HIST("Track1/eta"), track1.eta());
        dataTrackHist.fill(HIST("Track1/phi"), track1.phi());
        dataTrackHist.fill(HIST("Track1/tpcCluster"), track1.tpcNClsFound());
        dataTrackHist.fill(HIST("Track1/tpcCrossed"), track1.tpcNClsCrossedRows());
        dataTrackHist.fill(HIST("Track1/tpcShared"), track1.tpcNClsShared());
        dataTrackHist.fill(HIST("Track1/itsCluster"), track1.itsNCls());
        dataTrackHist.fill(HIST("Track1/itsIbCluster"), track1.itsNClsInnerBarrel());
        dataTrackHist.fill(HIST("Track1/itsNsigma"), track1.p(), nsigma1[0]);
        dataTrackHist.fill(HIST("Track1/tpcNsigma"), track1.p(), nsigma1[1]);
        dataTrackHist.fill(HIST("Track1/tpctofNsigma"), track1.p(), std::hypot(nsigma1[1], nsigma1[2]));

        auto nsigma2 = getNSigmaValues(track2, trackCuts2.pdgCode.value);
        dataTrackHist.fill(HIST("Track2/pt"), track2.pt());
        dataTrackHist.fill(HIST("Track2/eta"), track2.eta());
        dataTrackHist.fill(HIST("Track2/phi"), track2.phi());
        dataTrackHist.fill(HIST("Track2/tpcCluster"), track2.tpcNClsFound());
        dataTrackHist.fill(HIST("Track2/tpcCrossed"), track2.tpcNClsCrossedRows());
        dataTrackHist.fill(HIST("Track2/tpcShared"), track2.tpcNClsShared());
        dataTrackHist.fill(HIST("Track2/itsCluster"), track2.itsNCls());
        dataTrackHist.fill(HIST("Track2/itsIbCluster"), track2.itsNClsInnerBarrel());
        dataTrackHist.fill(HIST("Track2/itsNsigma"), track2.p(), nsigma2[0]);
        dataTrackHist.fill(HIST("Track2/tpcNsigma"), track2.p(), nsigma2[1]);
        dataTrackHist.fill(HIST("Track2/tpctofNsigma"), track1.p(), std::hypot(nsigma2[1], nsigma2[2]));

        pairs++;

        vecKstar.push_back(femtoDream::FemtoDreamMath::getkstar(track1, femtoDream::getMass(trackCuts1.pdgCode.value), track2, femtoDream::getMass(trackCuts2.pdgCode.value)));
      }
    }
    return pairs;
  }

  template <typename T>
  int countRecoTriplets(const T& tracks)
  {
    int triplets = 0;
    for (auto track1 = tracks.begin(); track1 != tracks.end(); track1++) {
      if (!checkRecoTrackSelections(track1, trackCuts1) || !checkRecoTrackPidSelections(track1, trackCuts1)) {
        continue;
      }
      for (auto track2 = track1 + 1; track2 != tracks.end(); track2++) {
        if (!checkRecoTrackSelections(track2, trackCuts2) || !checkRecoTrackPidSelections(track2, trackCuts2)) {
          continue;
        }

        for (auto track3 = track2 + 1; track3 != tracks.end(); track3++) {
          if (!checkRecoTrackSelections(track3, trackCuts3) || !checkRecoTrackPidSelections(track3, trackCuts3)) {
            continue;
          }

          auto nsigma1 = getNSigmaValues(track1, trackCuts1.pdgCode.value);
          dataTrackHist.fill(HIST("Track1/pt"), track1.pt());
          dataTrackHist.fill(HIST("Track1/eta"), track1.eta());
          dataTrackHist.fill(HIST("Track1/phi"), track1.phi());
          dataTrackHist.fill(HIST("Track1/tpcCluster"), track1.tpcNClsFound());
          dataTrackHist.fill(HIST("Track1/tpcCrossed"), track1.tpcNClsCrossedRows());
          dataTrackHist.fill(HIST("Track1/tpcShared"), track1.tpcNClsShared());
          dataTrackHist.fill(HIST("Track1/itsCluster"), track1.itsNCls());
          dataTrackHist.fill(HIST("Track1/itsIbCluster"), track1.itsNClsInnerBarrel());
          dataTrackHist.fill(HIST("Track1/itsNsigma"), track1.p(), nsigma1[0]);
          dataTrackHist.fill(HIST("Track1/tpcNsigma"), track1.p(), nsigma1[1]);
          dataTrackHist.fill(HIST("Track1/tpctofNsigma"), track1.p(), std::hypot(nsigma1[1], nsigma1[2]));

          auto nsigma2 = getNSigmaValues(track2, trackCuts2.pdgCode.value);
          dataTrackHist.fill(HIST("Track2/pt"), track2.pt());
          dataTrackHist.fill(HIST("Track2/eta"), track2.eta());
          dataTrackHist.fill(HIST("Track2/phi"), track2.phi());
          dataTrackHist.fill(HIST("Track2/tpcCluster"), track2.tpcNClsFound());
          dataTrackHist.fill(HIST("Track2/tpcCrossed"), track2.tpcNClsCrossedRows());
          dataTrackHist.fill(HIST("Track2/tpcShared"), track2.tpcNClsShared());
          dataTrackHist.fill(HIST("Track2/itsCluster"), track2.itsNCls());
          dataTrackHist.fill(HIST("Track2/itsIbCluster"), track2.itsNClsInnerBarrel());
          dataTrackHist.fill(HIST("Track2/itsNsigma"), track2.p(), nsigma2[0]);
          dataTrackHist.fill(HIST("Track2/tpcNsigma"), track2.p(), nsigma2[1]);
          dataTrackHist.fill(HIST("Track2/tpctofNsigma"), track2.p(), std::hypot(nsigma2[1], nsigma2[2]));

          auto nsigma3 = getNSigmaValues(track3, trackCuts3.pdgCode.value);
          dataTrackHist.fill(HIST("Track3/pt"), track3.pt());
          dataTrackHist.fill(HIST("Track3/eta"), track3.eta());
          dataTrackHist.fill(HIST("Track3/phi"), track3.phi());
          dataTrackHist.fill(HIST("Track3/tpcCluster"), track3.tpcNClsFound());
          dataTrackHist.fill(HIST("Track3/tpcCrossed"), track3.tpcNClsCrossedRows());
          dataTrackHist.fill(HIST("Track3/tpcShared"), track3.tpcNClsShared());
          dataTrackHist.fill(HIST("Track3/itsCluster"), track3.itsNCls());
          dataTrackHist.fill(HIST("Track3/itsIbCluster"), track3.itsNClsInnerBarrel());
          dataTrackHist.fill(HIST("Track3/itsNsigma"), track3.p(), nsigma2[0]);
          dataTrackHist.fill(HIST("Track3/tpcNsigma"), track3.p(), nsigma2[1]);
          dataTrackHist.fill(HIST("Track3/tpctofNsigma"), track3.p(), std::hypot(nsigma3[1], nsigma3[2]));

          triplets++;

          vecq3.push_back(femtoDream::FemtoDreamMath::getQ3(track1, femtoDream::getMass(trackCuts1.pdgCode.value), track2, femtoDream::getMass(trackCuts2.pdgCode.value), track3, femtoDream::getMass(trackCuts3.pdgCode.value)));
        }
      }
    }
    return triplets;
  }

  template <typename T>
  int countGenPairs(const T& tracks)
  {
    int pairs = 0;
    for (auto track1 = tracks.begin(); track1 != tracks.end(); track1++) {
      if (!track1.isPhysicalPrimary() ||
          track1.pdgCode() != trackCuts1.pdgCode.value ||
          track1.pt() < trackCuts1.ptMin.value ||
          track1.pt() > trackCuts1.ptMax.value ||
          std::fabs(track1.eta()) > trackCuts1.etaAbsMax.value) {
        continue;
      }
      for (auto track2 = track1 + 1; track2 != tracks.end(); track2++) {
        if (!track2.isPhysicalPrimary() ||
            track2.pdgCode() != trackCuts2.pdgCode.value ||
            track2.pt() < trackCuts2.ptMin.value ||
            track2.pt() > trackCuts2.ptMax.value ||
            std::fabs(track2.eta()) > trackCuts2.etaAbsMax.value) {
          continue;
        }
        vecKstar.push_back(femtoDream::FemtoDreamMath::getkstar(track1, femtoDream::getMass(trackCuts1.pdgCode.value), track2, femtoDream::getMass(trackCuts2.pdgCode.value)));
        pairs++;
      }
    }
    return pairs;
  }

  template <typename T>
  int countGenTriplets(const T& tracks)
  {
    int triplets = 0;
    for (auto track1 = tracks.begin(); track1 != tracks.end(); track1++) {
      if (!track1.isPhysicalPrimary() ||
          track1.pdgCode() != trackCuts1.pdgCode.value ||
          track1.pt() < trackCuts1.ptMin.value ||
          track1.pt() > trackCuts1.ptMax.value ||
          std::fabs(track1.eta()) > trackCuts1.etaAbsMax.value) {
        continue;
      }
      for (auto track2 = track1 + 1; track2 != tracks.end(); track2++) {
        if (!track2.isPhysicalPrimary() ||
            track2.pdgCode() != trackCuts2.pdgCode.value ||
            track2.pt() < trackCuts2.ptMin.value ||
            track2.pt() > trackCuts2.ptMax.value ||
            std::fabs(track2.eta()) > trackCuts2.etaAbsMax.value) {
          continue;
        }
        for (auto track3 = track2 + 1; track3 != tracks.end(); track3++) {
          if (!track3.isPhysicalPrimary() ||
              track3.pdgCode() != trackCuts3.pdgCode.value ||
              track3.pt() < trackCuts3.ptMin.value ||
              track3.pt() > trackCuts3.ptMax.value ||
              std::fabs(track3.eta()) > trackCuts3.etaAbsMax.value) {
            continue;
          }
          vecq3.push_back(femtoDream::FemtoDreamMath::getQ3(track1, femtoDream::getMass(trackCuts1.pdgCode.value), track2, femtoDream::getMass(trackCuts2.pdgCode.value), track3, femtoDream::getMass(trackCuts3.pdgCode.value)));
          triplets++;
        }
      }
    }
    return triplets;
  }

  void processdNdetaData(Collisions::iterator const& collision, FilteredFullTracks const& tracks)
  {
    if (!checkEventSelections<false>(collision, true))
      return;

    auto tracksWithItsPid = soa::Attach<FilteredFullTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);

    if (mode.countPairs.value) {
      vecKstar.clear();
      pairsFound = countRecoPairs(tracksWithItsPid);
      if (pairsFound == 0) {
        return;
      }
    }
    if (mode.countTriplets.value) {
      vecq3.clear();
      tripletsFound = countRecoTriplets(tracksWithItsPid);
      if (tripletsFound == 0) {
        return;
      }
    }

    dataEventHist.fill(HIST("hDataEventSelection"), 4); // found a pair

    if (mode.countPairs.value) {
      if (*std::min_element(vecKstar.begin(), vecKstar.end()) > eventSelection.kstarMax) {
        return;
      }
    }

    if (mode.countTriplets.value) {
      if (*std::min_element(vecq3.begin(), vecq3.end()) > eventSelection.q3Max) {
        return;
      }
    }

    dataEventHist.fill(HIST("hDataEventSelection"), 5); // found a lowkstar pair

    float centrality = collision.centFT0M();

    dataEventHist.fill(HIST("hDataVertexZ"), collision.posZ());
    dataEventHist.fill(HIST("hDataMultiplicity"), collision.multNTracksPV());
    dataEventHist.fill(HIST("hDataCentrality"), centrality);

    for (const auto& track : tracks)
      dataEventHist.fill(HIST("hDataCentralityVsEtaDistribtion"), centrality, track.eta());
  }

  PROCESS_SWITCH(FemtoDreamPairEfficiency, processdNdetaData, "Process function for dN/deta values in MCReco", true);

  void processdNdetaMCReco(RecoCollisions::iterator const& collision, FilteredFullMCTracks const& McTracks, GenCollisions const&, McParticles const& mcParticles)
  {
    if (!checkEventSelections<true>(collision, true)) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    mcEventHist.fill(HIST("hRecoEventSelection"), 4);

    const auto& mcCollision = collision.mcCollision_as<GenCollisions>();
    auto mcParticlesThisColl = mcParticles.sliceByCached(mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    if (mode.countPairs.value) {
      vecKstar.clear();
      pairsFound = countGenPairs(mcParticlesThisColl);
      if (pairsFound == 0) {
        return;
      }
    }
    if (mode.countTriplets.value) {
      vecq3.clear();
      tripletsFound = countGenTriplets(mcParticlesThisColl);
      if (tripletsFound == 0) {
        return;
      }
    }

    mcEventHist.fill(HIST("hRecoEventSelection"), 5);

    if (mode.countPairs.value) {
      if (*std::min_element(vecKstar.begin(), vecKstar.end()) > eventSelection.kstarMax) {
        return;
      }
    }

    if (mode.countTriplets.value) {
      if (*std::min_element(vecq3.begin(), vecq3.end()) > eventSelection.q3Max) {
        return;
      }
    }

    mcEventHist.fill(HIST("hRecoEventSelection"), 6);

    float genCentrality = mcCollision.centFT0M();

    mcEventHist.fill(HIST("hRecoMcVertexZ"), mcCollision.posZ());
    mcEventHist.fill(HIST("hRecoMcMultiplicity"), mcCollision.multMCNParticlesEta08());
    mcEventHist.fill(HIST("hRecoMcCentrality"), mcCollision.centFT0M());
    mcEventHist.fill(HIST("hRecoMcNumberOfPair"), vecKstar.size());
    for (const auto& kstar : vecKstar) {
      mcEventHist.fill(HIST("hRecoMcKstar"), kstar);
    }

    for (const auto& track : McTracks) {
      if (!track.has_mcParticle())
        continue;
      auto mcTrack = track.mcParticle_as<McParticles>();
      if (!mcTrack.isPhysicalPrimary() || std::fabs(mcTrack.eta()) > eventSelection.etaAbsMax.value)
        continue;

      mcEventHist.fill(HIST("hRecoMcCentralityVsEtaDistribution"), genCentrality, mcTrack.eta());
    }

    for (const auto& mcParticle : mcParticlesThisColl) {
      if (!mcParticle.isPhysicalPrimary() || std::fabs(mcParticle.eta()) > eventSelection.etaAbsMax.value) {
        continue;
      }
      auto pdgTrack = pdgDB->GetParticle(mcParticle.pdgCode());
      if (pdgTrack == nullptr) {
        continue;
      }
      if (pdgTrack->Charge() == 0) {
        continue;
      }
      mcEventHist.fill(HIST("hRecoMcCentralityVsMcEtaDistribution"), genCentrality, mcParticle.eta());
    }
  }
  PROCESS_SWITCH(FemtoDreamPairEfficiency, processdNdetaMCReco, "Process function for dN/deta values in MCReco", true);

  void processdNdetaMCGen(GenCollisions::iterator const& mcCollision, soa::SmallGroups<RecoCollisions> const& collisions, McParticles const& mcParticles)
  {
    mcEventHist.fill(HIST("hGenEventSelection"), 0);
    if (std::fabs(mcCollision.posZ()) > eventSelection.zvtxAbsMax.value) {
      return;
    }
    mcEventHist.fill(HIST("hGenEventSelection"), 1);
    if (!pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) {
      return;
    }
    mcEventHist.fill(HIST("hGenEventSelection"), 2);

    if (mode.countPairs.value) {
      vecKstar.clear();
      pairsFound = countGenPairs(mcParticles);
      if (pairsFound == 0) {
        return;
      }
    }
    if (mode.countTriplets.value) {
      vecq3.clear();
      tripletsFound = countGenTriplets(mcParticles);
      if (tripletsFound == 0) {
        return;
      }
    }

    mcEventHist.fill(HIST("hGenEventSelection"), 3);

    if (mode.countPairs.value) {
      if (*std::min_element(vecKstar.begin(), vecKstar.end()) > eventSelection.kstarMax) {
        return;
      }
    }

    if (mode.countTriplets.value) {
      if (*std::min_element(vecq3.begin(), vecq3.end()) > eventSelection.q3Max) {
        return;
      }
    }

    mcEventHist.fill(HIST("hGenEventSelection"), 4);
    mcEventHist.fill(HIST("hGenMcVertexZ"), mcCollision.posZ());
    mcEventHist.fill(HIST("hGenMcMultiplicity"), mcCollision.multMCNParticlesEta08());
    mcEventHist.fill(HIST("hGenMcCentrality"), mcCollision.centFT0M());
    mcEventHist.fill(HIST("hGenMcNumberOfPair"), vecKstar.size());

    for (const auto& kstar : vecKstar) {
      mcEventHist.fill(HIST("hGenMcKstar"), kstar);
    }

    float genCentrality = mcCollision.centFT0M();
    uint64_t numberAssocCollisions = 0;

    for (const auto& collision : collisions) {
      if (checkEventSelections<true>(collision, false)) {
        mcEventHist.fill(HIST("hGenMcCentralityWithAllAssoc"), genCentrality);
        for (const auto& mcParticle : mcParticles) {
          if (!mcParticle.isPhysicalPrimary() || std::abs(mcParticle.eta()) > eventSelection.etaAbsMax.value)
            continue;
          auto pdgTrack = pdgDB->GetParticle(mcParticle.pdgCode());
          if (pdgTrack == nullptr)
            continue;
          if (pdgTrack->Charge() == 0)
            continue;
          mcEventHist.fill(HIST("hGenMcCentralityVsMcEtaDistributionWithAllAssoc"), genCentrality, mcParticle.eta());
        }
        numberAssocCollisions++;
      }
    }

    if (numberAssocCollisions > 0) {
      mcEventHist.fill(HIST("hGenMcCentralityWithAssoc"), genCentrality);
      mcEventHist.fill(HIST("hGenEventSelection"), 5);
    }

    for (const auto& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary() || std::abs(mcParticle.eta()) > eventSelection.etaAbsMax.value)
        continue;

      auto pdgTrack = pdgDB->GetParticle(mcParticle.pdgCode());
      if (pdgTrack == nullptr) {
        continue;
      }
      if (pdgTrack->Charge() == 0) {
        continue;
      }

      mcEventHist.fill(HIST("hGenMcCentralityVsMcEtaDistribution"), genCentrality, mcParticle.eta());
      if (numberAssocCollisions > 0) {
        mcEventHist.fill(HIST("hGenMcCentralityVsMcEtaDistributionWithAssoc"), genCentrality, mcParticle.eta());
      }
    }
  }
  PROCESS_SWITCH(FemtoDreamPairEfficiency, processdNdetaMCGen, "Process function for dN/deta values in MCReco", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoDreamPairEfficiency>(cfgc)};
  return workflow;
}
