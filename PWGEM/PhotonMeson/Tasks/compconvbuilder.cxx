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

/// \file compconvbuilder.cxx
/// \brief QA task for photons in the EM and LF builder
/// \author S. Mrozinski, smrozins@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TPDGCode.h>

#include <cstdlib>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::constants;

using namespace o2::aod::pwgem::dilepton::utils::mcutil;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

using MyCollisions = soa::Join<aod::EMEvents_004, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::BinnedGenPts>;
using MyMCCollision = MyMCCollisions::iterator;

using MyStraCollisions = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels>;
using MyStraCollision = MyStraCollisions::iterator;

using MyTracksIUMC = soa::Join<aod::TracksIU, aod::McTrackLabels>;

using V0DerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0Indices, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0GammaMLScores>;

using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

struct Compconvbuilder {
  HistogramRegistry registry{"Compconvbuilder"};

  enum ConversionBuilderID {
    EMBuilder = 0,
    LFBuilder = 1,
    EMOnly = 2,
    LFOnly = 3,
    Common = 4,
    NConversionBuilder
  };

  static constexpr std::string_view kConversionBuilder[NConversionBuilder] = {"EMBuilder/", "LFBuilder/", "EMOnly/", "LFOnly/", "Common/"};
  static constexpr std::string_view kEventTypes[2] = {"before/", "after/"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
  } eventcuts;

  void defineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
  }

  // Link V0-photons to their collision
  Preslice<MyV0Photons> perV0PhotonCollision = aod::v0photonkf::emphotoneventId;

  void init(InitContext const& /*ctx*/)
  {

    defineEMEventCut();

    for (int i = 0; i < NConversionBuilder; ++i) {

      registry.add<TH1>(string(kConversionBuilder[i]) + "hPt", ";p_{T} (GeV/c); Counts", kTH1F, {{1000, 0., 10.}});
      registry.add<TH1>(string(kConversionBuilder[i]) + "hR", ";R_{conv} (cm); Counts", kTH1F, {{100, 0., 100.}});

      registry.add<TH1>(string(kConversionBuilder[i]) + "hEta", ";#eta; Counts", kTH1F, {{200, -1.0f, 1.0f}});

      registry.add<TH1>(string(kConversionBuilder[i]) + "hcosPA", ";R_{conv} (cm); Counts", kTH1F, {{100, 0.99f, 1.0f}});

      registry.add<TH1>(string(kConversionBuilder[i]) + "MatchedDeltaRec", ";#Delta colID_{rec};Counts", kTH1F, {{21, -10.5, 10.5}});

      registry.add<TH1>(string(kConversionBuilder[i]) + "hZ", ";z (cm);Counts", kTH1F, {{200, -100, 100}});

      registry.add<TH2>(string(kConversionBuilder[i]) + "hZR", "conversion point in RZ;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100, 100}, {200, 0.0f, 100.0f}});

      registry.add<TH2>(string(kConversionBuilder[i]) + "hAP", "AP plot;#alpha;q_{T} (GeV/c)", kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}});

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionGen/Z_res", "Conversion radius resolution;z_{conv, gen} (cm);R_{conv, gen} (cm);#varphi_{gen} (rad.);#eta_{gen};p_{T, gen};z_{conv, rec} - z_{conv, gen} (cm);",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {120, -30, 30}

                              },
                              false);

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionGen/R_res", "Conversion radius resolution;z_{conv, gen} (cm);R_{conv, gen} (cm);#varphi_{gen} (rad.);#eta_{gen};p_{T, gen};R_{conv, rec} - R_{conv, gen} (cm);",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {120, -30, 30}

                              },
                              false);

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionGen/Phi_res", "#varphi resolution;z_{conv, gen} (cm);R_{conv, gen} (cm);#varphi_{gen} (rad.);#eta_{gen};p_{T, gen};#varphi_{conv, rec} - #varphi_{conv, gen} (cm);",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {100, -0.2f, 0.2f}

                              },
                              false);

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionGen/Pt_res", "Conversion radius resolution;z_{conv, gen} (cm);R_{conv, gen} (cm);#varphi_{gen} (rad.);#eta_{gen};p_{T, gen};p_{T, rec} - p_{T, gen}/p_{T, gen};",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {200, -1.0f, 1.0f}

                              },
                              false);

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionGen/Eta_res", "Conversion radius resolution;z_{conv, gen} (cm);R_{conv, gen} (cm);#varphi_{gen} (rad.);#eta_{gen};p_{T, gen};#eta_{conv, rec} - #eta_{conv, gen} (cm);",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {100, -0.5f, 0.5f}

                              },
                              false);

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionRec/Z_res", "Conversion radius resolution;z_{conv, rec} (cm);R_{conv, rec} (cm);#varphi_{rec} (rad.);#eta_{rec};p_{T, rec};z_{conv, rec} - z_{conv, gen} (cm);",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {120, -30, 30}

                              },
                              false);

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionRec/R_res", "Conversion radius resolution;z_{conv, rec} (cm);R_{conv, rec} (cm);#varphi_{rec} (rad.);#eta_{rec};p_{T, rec};R_{conv, rec} - R_{conv, gen} (cm);",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {120, -30, 30}

                              },
                              false);

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionRec/Phi_res", "#varphi resolution;z_{conv, rec} (cm);R_{conv, rec} (cm);#varphi_{rec} (rad.);#eta_{rec};p_{T, rec};#varphi_{conv, rec} - #varphi_{conv, gen} (cm);",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {100, -0.2f, 0.2f}

                              },
                              false);

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionRec/Pt_res", "Conversion radius resolution;z_{conv, rec} (cm);R_{conv, rec} (cm);#varphi_{rec} (rad.);#eta_{rec};p_{T, rec};p_{T, rec} - p_{T, gen}/p_{T, gen};",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {200, -1.0f, 1.0f}

                              },
                              false);

      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ResolutionRec/Eta_res", "Conversion radius resolution;z_{conv, rec} (cm);R_{conv, rec} (cm);#varphi_{rec} (rad.);#eta_{rec};p_{T, rec};#eta_{conv, rec} - #eta_{conv, gen} (cm);",
                              kTHnSparseF,
                              {{200, -100, 100},
                               {200, 0, 100},
                               {90, 0, math::TwoPI},
                               {200, -1.0f, 1.0f},
                               {500, 0, 10},
                               {100, -0.5f, 0.5f}

                              },
                              false);
      registry.add<THnSparse>(string(kConversionBuilder[i]) + "ConvInfo", "Conversion radius resolution;x_{conv} (cm);y_{conv} (cm);z_{conv} (cm);R_{conv} (cm);#varphi (rad.);#eta;p_{T, gen};",
                              kTHnSparseF,
                              {
                                {200, -100, 100},
                                {200, -100, 100},
                                {200, -100, 100},
                                {200, 0, 100},
                                {90, 0, math::TwoPI},
                                {200, -1.0f, 1.0f},
                                {500, 0, 10},

                              },
                              false);

      registry.add<TH1>(string(kConversionBuilder[i]) + "V0Leg/Asymmetry", "", kTH1F, {{100, 0, 1}});

      auto hCollisionCounter = registry.add<TH1>(string(kConversionBuilder[i]) + "Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1F, {{10, 0.5, 10.5}}, false);
      hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
      hCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
      hCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
      hCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
      hCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
      hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
      hCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
      hCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
      hCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
      hCollisionCounter->GetXaxis()->SetBinLabel(10, "accepted");

      registry.add(string(kConversionBuilder[i]) + "Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
      registry.add(string(kConversionBuilder[i]) + "Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
      registry.add(string(kConversionBuilder[i]) + "Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
      registry.add(string(kConversionBuilder[i]) + "Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{300, 0, 6000}, {300, 0, 6000}}, false);
      registry.add(string(kConversionBuilder[i]) + "Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
      registry.add(string(kConversionBuilder[i]) + "Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
      registry.add(string(kConversionBuilder[i]) + "Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
      registry.add(string(kConversionBuilder[i]) + "Event/before/hCentFT0MvsMultNTracksPV", "hCentFT0MvsMultNTracksPV;centrality FT0M (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {600, 0, 6000}}, false);
      registry.add(string(kConversionBuilder[i]) + "Event/before/hMultFT0MvsMultNTracksPV", "hMultFT0MvsMultNTracksPV;mult. FT0M;N_{track} to PV", kTH2F, {{600, 0, 6000}, {600, 0, 6000}}, false);
      registry.addClone(string(kConversionBuilder[i]) + "Event/before/", string(kConversionBuilder[i]) + "Event/after/");
    }

    registry.add<TH1>("truePhotons/hPt_Converted", "Converted Photons; p_{T} (GeV/c); Counts", kTH1F, {{100, 0., 10.}});
    registry.add<TH1>("truePhotons/hR_Converted", "Converted Photons; R (cm); Counts", kTH1F, {{100, 0., 100.}});

    registry.add<THnSparse>("truePhotons/Sparse_Converted", "Conversion radius resolution;x_{conv} (cm);z_{conv} (cm);y_{conv} (cm);R_{conv} (cm);#varphi (rad.);#eta;p_{T, gen};",
                            kTHnSparseF,
                            {{200, -100, 100},
                             {200, -100, 100},
                             {200, -100, 100},
                             {200, 0, 100},
                             {90, 0, math::TwoPI},
                             {200, -1.0f, 1.0f},
                             {500, 0, 10}},
                            false);

    auto h = registry.add<TH1>("EMBuilder/hV0SignType", "Crosscheck", kTH1F, {{3, 0.5, 4.5}}, false);
    h->GetXaxis()->SetBinLabel(1, "Same-sign");
    h->GetXaxis()->SetBinLabel(2, "Opposite-sign");
    h->GetXaxis()->SetBinLabel(3, "Zero-sign");

    auto h2 = registry.add<TH1>("EMBuilder/hV0ElectronPositronTrue", "pair in MC truth;;counts", kTH1F, {{2, -0.5, 1.5}});
    h2->GetXaxis()->SetBinLabel(1, "Mismatch");
    h2->GetXaxis()->SetBinLabel(2, "Good");
  }

  template <const int ev_id, int type, typename TCollision>
  void fillEventInfo(TCollision const& collision, const float /*weight*/ = 1.f)
  {
    registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 1.0);

    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 2.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 3.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 4.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 5.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 6.0);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 7.0);
    }

    if (collision.sel8()) {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 8.0);
    }
    if (std::fabs(collision.posZ()) < eventcuts.cfgZvtxMax) {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hCollisionCounter"), 9.0);
    }

    registry.fill(HIST(kConversionBuilder[type]) + HIST("Event/") + HIST(kEventTypes[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
  }

  template <int type>
  void fillLegInfo(auto& v0, auto& posleg, auto& negleg)
  {
    if constexpr (type == 0) {

      float ptPos = posleg.pt();
      float ptNeg = negleg.pt();
      float asym = (ptPos - ptNeg) / (ptPos + ptNeg);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("V0Leg/Asymmetry"), asym);

    } else {

      float ptPos = v0.postrackpt();
      float ptNeg = negleg.negtrackpt();
      float asym = (ptPos - ptNeg) / (ptPos + ptNeg);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("V0Leg/Asymmetry"), asym);
    }
  }

  template <int type>
  void fillV0Info(auto& v0, auto& v0MC, auto& mcleg)
  {
    registry.fill(HIST(kConversionBuilder[type]) + HIST("hPt"), v0.pt());
    registry.fill(HIST(kConversionBuilder[type]) + HIST("hEta"), v0.eta());
    registry.fill(HIST(kConversionBuilder[type]) + HIST("hAP"), v0.alpha(), v0.qtarm());

    if constexpr (type == EMBuilder || type == EMOnly) {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("hZ"), v0.vz());
      registry.fill(HIST(kConversionBuilder[type]) + HIST("hcosPA"), v0.cospa());
      registry.fill(HIST(kConversionBuilder[type]) + HIST("hZR"), v0.vz(), v0.v0radius());

      float deltapT = v0.pt() - v0MC.pt();
      float deltaZ = v0.vz() - mcleg.vz();
      float deltaPhi = v0.phi() - v0MC.phi();
      float deltaEta = v0.eta() - v0MC.eta();
      float deltaR = v0.v0radius() - std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2));

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/Z_res"),
                    mcleg.vz(),
                    std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)),
                    v0MC.phi(),
                    v0MC.eta(),
                    v0MC.pt(),
                    deltaZ);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/R_res"),
                    mcleg.vz(),
                    std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)),
                    v0MC.phi(),
                    v0MC.eta(),
                    v0MC.pt(),
                    deltaR);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/Phi_res"),
                    mcleg.vz(),
                    std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)),
                    v0MC.phi(),
                    v0MC.eta(),
                    v0MC.pt(),
                    deltaPhi);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/Pt_res"),
                    mcleg.vz(),
                    std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)),
                    v0MC.phi(),
                    v0MC.eta(),
                    v0MC.pt(),
                    deltapT);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/Eta_res"),
                    mcleg.vz(),
                    std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)),
                    v0MC.phi(),
                    v0MC.eta(),
                    v0MC.pt(),
                    deltaEta);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/Z_res"),
                    v0.vz(),
                    v0.v0radius(),
                    v0.phi(),
                    v0.eta(),
                    v0.pt(),
                    deltaZ);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/R_res"),
                    v0.vz(),
                    v0.v0radius(),
                    v0.phi(),
                    v0.eta(),
                    v0.pt(),
                    deltaR);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/Phi_res"),
                    v0.vz(),
                    v0.v0radius(),
                    v0.phi(),
                    v0.eta(),
                    v0.pt(),
                    deltaPhi);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/Pt_res"),
                    v0.vz(),
                    v0.v0radius(),
                    v0.phi(),
                    v0.eta(),
                    v0.pt(),
                    deltapT);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/Eta_res"),
                    v0.vz(),
                    v0.v0radius(),
                    v0.phi(),
                    v0.eta(),
                    v0.pt(),
                    deltaEta);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ConvInfo"),
                    v0.vx(),       // 0
                    v0.vy(),       // 1
                    v0.vz(),       // 2
                    v0.v0radius(), // 3
                    v0.phi(),      // 4
                    v0.eta(),      // 5
                    v0.pt());      // 6

    } else {
      registry.fill(HIST(kConversionBuilder[type]) + HIST("hZ"), v0.z());
      registry.fill(HIST(kConversionBuilder[type]) + HIST("hcosPA"), v0.v0cosPA());
      registry.fill(HIST(kConversionBuilder[type]) + HIST("hZR"), v0.z(), v0.v0radius());

      float deltaR = v0.v0radius() - std::hypot(v0MC.xMC(), v0MC.yMC());

      float phiRec = v0.phi();
      RecoDecay::constrainAngle(phiRec);

      float phiMC = std::atan2(v0MC.pyMC(), v0MC.pxMC());
      RecoDecay::constrainAngle(phiMC);

      float deltaPhi = phiRec - phiMC;
      RecoDecay::constrainAngle(deltaPhi);

      float etaGen = 0.5f * std::log((std::hypot(v0MC.pxMC(), v0MC.pyMC(), v0MC.pzMC()) + v0MC.pzMC()) / (std::hypot(v0MC.pxMC(), v0MC.pyMC(), v0MC.pzMC()) - v0MC.pzMC()));

      float etaRec = 0.5f * std::log((std::hypot(v0.px(), v0.py(), v0.pz()) + v0.pz()) / (std::hypot(v0.px(), v0.py(), v0.pz()) - v0.pz()));

      float deltapT = v0.pt() - v0MC.ptMC();
      float deltaEta = etaRec - etaGen;
      float deltaZ = v0.z() - v0MC.zMC();

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/Z_res"),
                    v0MC.zMC(),
                    std::hypot(v0MC.xMC(), v0MC.yMC()),
                    phiMC,
                    etaGen,
                    v0MC.ptMC(),
                    deltaZ);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/R_res"),
                    v0MC.zMC(),
                    std::hypot(v0MC.xMC(), v0MC.yMC()),
                    phiMC,
                    etaGen,
                    v0MC.ptMC(),
                    deltaR);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/Phi_res"),
                    v0MC.zMC(),
                    std::hypot(v0MC.xMC(), v0MC.yMC()),
                    phiMC,
                    etaGen,
                    v0MC.ptMC(),
                    deltaPhi);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/Pt_res"),
                    v0MC.zMC(),
                    std::hypot(v0MC.xMC(), v0MC.yMC()),
                    phiMC,
                    etaGen,
                    v0MC.ptMC(),
                    deltapT);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionGen/Eta_res"),
                    v0MC.zMC(),
                    std::hypot(v0MC.xMC(), v0MC.yMC()),
                    phiMC,
                    etaGen,
                    v0MC.ptMC(),
                    deltaEta);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/Z_res"),
                    v0.z(),
                    v0.v0radius(),
                    phiRec,
                    v0.eta(),
                    v0.pt(),
                    deltaZ);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/R_res"),
                    v0.z(),
                    v0.v0radius(),
                    phiRec,
                    v0.eta(),
                    v0.pt(),
                    deltaR);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/Phi_res"),
                    v0.z(),
                    v0.v0radius(),
                    phiRec,
                    v0.eta(),
                    v0.pt(),
                    deltaPhi);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/Eta_res"),
                    v0.z(),
                    v0.v0radius(),
                    phiRec,
                    v0.eta(),
                    v0.pt(),
                    deltaEta);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ResolutionRec/Pt_res"),
                    v0.z(),
                    v0.v0radius(),
                    phiRec,
                    v0.eta(),
                    v0.pt(),
                    deltapT);

      registry.fill(HIST(kConversionBuilder[type]) + HIST("ConvInfo"),
                    v0.x(),
                    v0.y(),
                    v0.z(),
                    v0.v0radius(),
                    phiRec,
                    v0.eta(),
                    v0.pt());
    }

    registry.fill(HIST(kConversionBuilder[type]) + HIST("hR"), v0.v0radius());
  }

  Preslice<V0DerivedMCDatas> perCollisionMCDerived = o2::aod::v0data::straCollisionId;

  void processLFV0sMC(MyStraCollisions const& stracollisions,
                      soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&,
                      V0DerivedMCDatas const& strangeV0s,
                      DauTracks const&)
  {

    for (const auto& collision : stracollisions) {

      fillEventInfo<0, LFBuilder>(collision);

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      fillEventInfo<1, LFBuilder>(collision);

      registry.fill(HIST((kConversionBuilder[1])) + HIST("Event/before/hCollisionCounter"), 10.0);
      registry.fill(HIST((kConversionBuilder[1])) + HIST("Event/after/hCollisionCounter"), 10.0); // accepted

      auto myV0s = strangeV0s.sliceBy(perCollisionMCDerived, collision.globalIndex());

      for (auto const& v0 : myV0s) {
        if (!v0.has_v0MCCore()) {
          continue;
        }

        auto v0MC = v0.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

        if (v0MC.pdgCode() != kGamma || !v0MC.isPhysicalPrimary()) {
          continue;
        }

        auto posTrack = v0.template posTrackExtra_as<DauTracks>();

        fillV0Info<LFBuilder>(v0, v0MC, posTrack);
      }
    }
  }

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emphotoneventId;

  void processEMV0sMC(MyV0Photons const& v0s, aod::EMMCParticles const& mcparticles, MyMCV0Legs const&, MyCollisions const& collisions)
  {

    for (const auto& collision : collisions) {

      fillEventInfo<0, EMBuilder>(collision);

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      fillEventInfo<1, EMBuilder>(collision);

      registry.fill(HIST((kConversionBuilder[0])) + HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      registry.fill(HIST((kConversionBuilder[0])) + HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      auto v0PhotonsColl = v0s.sliceBy(perCollision, collision.globalIndex());

      for (auto const& v0 : v0PhotonsColl) {

        auto pos = v0.posTrack_as<MyMCV0Legs>();
        auto ele = v0.negTrack_as<MyMCV0Legs>();
        auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
        auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

        int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, kPositron, kElectron, kGamma, mcparticles);

        auto mcphoton = mcparticles.iteratorAt(photonid);

        if (mcphoton.isPhysicalPrimary()) {

          fillV0Info<EMBuilder>(v0, mcphoton, elemc);
        }

        if (pos.sign() * ele.sign() > 0) {
          registry.fill(HIST("EMBuilder/hV0SignType"), 1); // same-sign
        } else if (pos.sign() * ele.sign() < 0) {
          registry.fill(HIST("EMBuilder/hV0SignType"), 2); // opposite-sign
        } else {
          registry.fill(HIST("EMBuilder/hV0SignType"), 3); // zero or undefined
        }

        if ((posmc.pdgCode() == kElectron && elemc.pdgCode() == kPositron) || (posmc.pdgCode() == kPositron && elemc.pdgCode() == kElectron)) {
          registry.fill(HIST("EMBuilder/hV0ElectronPositronTrue"), 1); // good
        } else {
          registry.fill(HIST("EMBuilder/hV0ElectronPositronTrue"), 0); // mismatch
        }
      }
    }
  }

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;

  void processConvV0s(MyCollisions const& collisions,
                      MyMCCollisions const&,
                      aod::EMMCParticles const& mcparticles,
                      MyTracksIUMC const& tracks)
  {
    for (const auto& collision : collisions) {

      const float minR = 5.0f;
      const float maxR = 90.f;

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto mccollision = collision.template emmcevent_as<MyMCCollisions>();

      auto mcstack = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      auto mcTracksColl = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());

      std::unordered_map<int, int> mc2trk;
      for (const auto& trk : tracks) {
        if (trk.mcParticleId() >= 0) {
          mc2trk[trk.mcParticleId()] = trk.globalIndex();
        }
      }

      for (const auto& mc : mcTracksColl) {
        if (mc.pdgCode() != kGamma || !mc.isPhysicalPrimary()) {
          continue;
        }

        auto daughters = mc.daughtersIds();
        if (daughters.size() != 2) { // o2-linter: disable=magic-number (this is pretty clear and not magic)
          continue;
        }

        auto d1 = mcparticles.iteratorAt(daughters[0]);
        auto d2 = mcparticles.iteratorAt(daughters[1]);

        if (std::abs(d1.pdgCode()) != kElectron || std::abs(d2.pdgCode()) != kElectron) {
          continue;
        }

        float r = std::hypot(d1.vx(), d1.vy());
        if (r < minR || r > maxR) {
          continue;
        }

        registry.fill(HIST("truePhotons/hPt_Converted"), mc.pt());
        registry.fill(HIST("truePhotons/hR_Converted"), r);

        registry.fill(HIST("truePhotons/Sparse_Converted"), d1.vx(), mc.y(), d1.vz(), r, mc.phi(), mc.eta(), mc.pt());

        int id1 = mc2trk.count(d1.globalIndex()) ? mc2trk[d1.globalIndex()] : -1;
        int id2 = mc2trk.count(d2.globalIndex()) ? mc2trk[d2.globalIndex()] : -1;
        if (id1 < 0 || id2 < 0) {
          continue;
        }
      }
    }
  }

  void processMatchCategories(
    MyCollisions const& collisions,
    aod::EMMCEvents const&,
    MyTracksIUMC const& tracksgen,
    MyV0Photons const& emV0s,
    soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&,
    V0DerivedMCDatas const& lfV0s,
    MyMCV0Legs const&,
    aod::EMMCParticles const&,
    aod::McParticles const& mcparticles,
    DauTracks const&)
  {
    std::unordered_map<int, int> trackToMcLabel;
    for (auto const& t : tracksgen) {
      int label = t.mcParticleId();
      if (label >= 0) {
        trackToMcLabel[t.globalIndex()] = label;
      }
    }

    for (const auto& collision : collisions) {
      if (!fEMEventCut.IsSelected(collision))
        continue;

      fillEventInfo<1, EMBuilder>(collision);

      auto emSlice = emV0s.sliceBy(perCollision, collision.globalIndex());
      auto lfSlice = lfV0s.sliceBy(perCollisionMCDerived, collision.globalIndex());

      using EMIt = decltype(emSlice.begin());
      using LFIt = decltype(lfSlice.begin());
      struct Entry {
        std::optional<EMIt> emIt;
        std::optional<LFIt> lfIt;
      };
      std::unordered_map<int, Entry> table;

      for (EMIt it = emSlice.begin(); it != emSlice.end(); ++it) {
        auto posmc = it.posTrack_as<MyMCV0Legs>()
                       .emmcparticle_as<aod::EMMCParticles>();
        auto negmc = it.negTrack_as<MyMCV0Legs>()
                       .emmcparticle_as<aod::EMMCParticles>();
        int pid = FindCommonMotherFrom2Prongs(posmc, negmc,
                                              kPositron, kElectron, kGamma, mcparticles);
        if (pid >= 0)
          table[pid].emIt = it;
      }

      for (LFIt it = lfSlice.begin(); it != lfSlice.end(); ++it) {
        int posTrackIndex = it.posTrackId();
        auto negTrackIndex = it.negTrackId();

        if (!trackToMcLabel.count(posTrackIndex) || !trackToMcLabel.count(negTrackIndex))
          continue;
        auto posmc = mcparticles.iteratorAt(trackToMcLabel[posTrackIndex]);
        auto negmc = mcparticles.iteratorAt(trackToMcLabel[negTrackIndex]);
        int pid = FindCommonMotherFrom2Prongs(posmc, negmc,
                                              kPositron, kElectron, kGamma, mcparticles);
        if (pid >= 0)
          table[pid].lfIt = it;
      }

      for (auto const& [pid, entry] : table) {
        auto mcphoton = mcparticles.iteratorAt(pid);

        if (entry.emIt.has_value() && entry.lfIt.has_value()) {
          // --- Common V0 ---
          auto& lfV0 = *entry.lfIt.value();
          auto v0MC = lfV0.template v0MCCore_as<
            soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
          auto posTrack = lfV0.template posTrackExtra_as<DauTracks>();
          fillV0Info<Common>(lfV0, v0MC, posTrack);

        } else if (entry.emIt.has_value()) {
          // --- EM-only V0 ---
          auto& emV0 = *entry.emIt.value();
          auto negmc = emV0.negTrack_as<MyMCV0Legs>()
                         .emmcparticle_as<aod::EMMCParticles>();
          fillV0Info<EMOnly>(emV0, mcphoton, negmc);

        } else if (entry.lfIt.has_value()) {
          // --- LF-only V0 ---
          auto& lfV0 = *entry.lfIt.value();
          auto v0MC = lfV0.template v0MCCore_as<
            soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
          auto posTrack = lfV0.template posTrackExtra_as<DauTracks>();
          fillV0Info<LFOnly>(lfV0, v0MC, posTrack);
        }
      }
    }
  }

  PROCESS_SWITCH(Compconvbuilder, processMatchCategories, "Process V0s matched via MC photon", false);
  PROCESS_SWITCH(Compconvbuilder, processLFV0sMC, "Process LF Builder V0s", true);
  PROCESS_SWITCH(Compconvbuilder, processEMV0sMC, "Process EM Builder V0s", false);
  PROCESS_SWITCH(Compconvbuilder, processConvV0s, "Process generated converted V0s", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<Compconvbuilder>(cfg)};
}
