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
//
// ========================
//
// This code loops over v0 photons and makes pairs for neutral mesons analyses.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include <map>
#include "TString.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

struct PCMQC {
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
    },
  };
  void addhistograms()
  {
    const TString tracktype[2] = {"NonAmb", "Amb"};
    const TString pairtype[2] = {"NonAmb", "Amb"};

    // for single tracks
    for (int i = 0; i < 2; i++) {
      // registry.add(Form("%sTrack/hNbc", tracktype[i].Data()), "number of bcs", HistType::kTH1F, {{101, -0.5, 100.5}});
      // registry.add(Form("%sTrack/hNcoll", tracktype[i].Data()), "number of collisions", HistType::kTH1F, {{101, -0.5, 100.5}});
      registry.add(Form("%sTrack/hPt", tracktype[i].Data()), "pT", HistType::kTH1F, {{1000, 0.0f, 10}});
      registry.add(Form("%sTrack/hEtaPhi", tracktype[i].Data()), "#eta vs. #varphi", HistType::kTH2F, {{180, 0, TMath::TwoPi()}, {40, -2.0f, 2.0f}});
      registry.add(Form("%sTrack/hDCAxyz", tracktype[i].Data()), "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}});
      registry.add(Form("%sTrack/hNclsTPC", tracktype[i].Data()), "number of TPC clusters", HistType::kTH1F, {{161, -0.5, 160.5}});
      registry.add(Form("%sTrack/hNcrTPC", tracktype[i].Data()), "number of TPC crossed rows", HistType::kTH1F, {{161, -0.5, 160.5}});
      registry.add(Form("%sTrack/hChi2TPC", tracktype[i].Data()), "chi2/number of TPC clusters", HistType::kTH1F, {{100, 0, 10}});
      registry.add(Form("%sTrack/hTPCdEdx", tracktype[i].Data()), "TPC dE/dx", HistType::kTH2F, {{1000, 0, 10}, {200, 0, 200}});
      registry.add(Form("%sTrack/hTPCNsigmaEl", tracktype[i].Data()), "TPC n sigma el", HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}});
      registry.add(Form("%sTrack/hTPCNsigmaPi", tracktype[i].Data()), "TPC n sigma pi", HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}});
      registry.add(Form("%sTrack/hTPCNcr2Nf", tracktype[i].Data()), "TPC Ncr/Nfindable", HistType::kTH1F, {{200, 0, 2}});
      registry.add(Form("%sTrack/hNclsITS", tracktype[i].Data()), "number of ITS clusters", HistType::kTH1F, {{8, -0.5, 7.5}});
      registry.add(Form("%sTrack/hChi2ITS", tracktype[i].Data()), "chi2/number of ITS clusters", HistType::kTH1F, {{36, 0, 36}});
      registry.add(Form("%sTrack/hDCAxyPosToPV", pairtype[i].Data()), "hDCAxyPosToPV;DCA_{xy} (cm);", HistType::kTH1F, {{100, -5.0f, 5.0f}});
      registry.add(Form("%sTrack/hDCAxyNegToPV", pairtype[i].Data()), "hDCAxyNegToPV;DCA_{xy} (cm);", HistType::kTH1F, {{100, -5.0f, 5.0f}});
      registry.add(Form("%sTrack/hDCAzPosToPV", pairtype[i].Data()), "hDCAzPosToPV;DCA_{z} (cm);", HistType::kTH1F, {{100, -5.0f, 5.0f}});
      registry.add(Form("%sTrack/hDCAzNegToPV", pairtype[i].Data()), "hDCAzNegToPV;DCA_{z} (cm);", HistType::kTH1F, {{100, -5.0f, 5.0f}});
    }

    // for V0s
    for (int i = 0; i < 2; i++) {
      registry.add(Form("%sV0/hNgamma", pairtype[i].Data()), "Number of #photon candidates;Number of #gamma candidates", HistType::kTH1F, {{101, -0.5f, 100.5}});
      registry.add(Form("%sV0/hPt", pairtype[i].Data()), "pT;p_{T} (GeV/c)", HistType::kTH1F, {{1000, 0.0f, 10}});
      registry.add(Form("%sV0/hEtaPhi", pairtype[i].Data()), "#eta vs. #varphi;#varphi (rad.);#eta", HistType::kTH2F, {{180, 0, TMath::TwoPi()}, {40, -2.0f, 2.0f}});
      registry.add(Form("%sV0/hRadius", pairtype[i].Data()), "V0Radius; radius in Z (cm);radius in XY (cm)", HistType::kTH2F, {{500, -250, 250}, {2500, 0.0f, 250.0f}});
      registry.add(Form("%sV0/hCosPA", pairtype[i].Data()), "V0CosPA;cosine pointing angle", HistType::kTH1F, {{200, 0.8f, 1.0f}});
      registry.add(Form("%sV0/hPCA", pairtype[i].Data()), "distance between 2 legs; PCA (cm)", HistType::kTH1F, {{100, 0.0f, 10.0f}});
      registry.add(Form("%sV0/hAPplot", pairtype[i].Data()), "AP plot;#alpha;q_{T} (GeV/c)", HistType::kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}});
      registry.add(Form("%sV0/hGammaPsiPair", pairtype[i].Data()), "#psi_{pair} for photon conversion;#psi_{pair} (rad.);m_{ee} (GeV/c^{2})", HistType::kTH2F, {{160, 0.0, TMath::PiOver2()}, {100, 0.0f, 0.1f}});
      registry.add(Form("%sV0/hMassGamma", pairtype[i].Data()), "hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", HistType::kTH2F, {{900, 0.0f, 90.0f}, {100, 0.0f, 0.1f}});
      registry.add(Form("%sV0/hGammaRxy", pairtype[i].Data()), "conversion point in XY;V_{x} (cm);V_{y} (cm)", HistType::kTH2F, {{1800, -90.0f, 90.0f}, {1800, -90.0f, 90.0f}});
    }
  }

  void init(InitContext& context)
  {
    addhistograms();
  }

  template <int mode, typename T>
  void fillHistosLeg(const T& leg)
  {
    static constexpr std::string_view typenames[] = {"NonAmb", "Amb"};
    // registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hNbc"        ) , leg.numBC());
    // registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hNcoll"      ) , leg.numColl());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hPt"), leg.pt());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hEtaPhi"), leg.phi(), leg.eta());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hDCAxyz"), leg.dcaXY(), leg.dcaZ());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hNclsTPC"), leg.tpcNClsFound());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hNclsITS"), leg.itsNCls());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hNcrTPC"), leg.tpcNClsCrossedRows());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hTPCNcr2Nf"), leg.tpcCrossedRowsOverFindableCls());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hChi2TPC"), leg.tpcChi2NCl());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hChi2ITS"), leg.itsChi2NCl());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hTPCdEdx"), leg.tpcInnerParam(), leg.tpcSignal());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl());
    registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hTPCNsigmaPi"), leg.tpcInnerParam(), leg.tpcNSigmaPi());

    if (leg.sign() > 0) {
      registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hDCAxyPosToPV"), leg.dcaXY());
      registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hDCAzPosToPV"), leg.dcaZ());
    } else {
      registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hDCAxyNegToPV"), leg.dcaXY());
      registry.fill(HIST(typenames[mode]) + HIST("Track/") + HIST("hDCAzNegToPV"), leg.dcaZ());
    }
  }

  template <int mode, typename T>
  void fillHistosV0(const T& v0)
  {
    static constexpr std::string_view typenames[] = {"NonAmb", "Amb"};
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hPt"), v0.pt());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hEtaPhi"), v0.phi(), v0.eta());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hRadius"), v0.vz(), v0.v0radius());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hCosPA"), v0.cospa());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hPCA"), v0.pca());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hAPplot"), v0.alpha(), v0.qtarm());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hMassGamma"), v0.v0radius(), v0.mGamma());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hGammaPsiPair"), v0.psipair(), v0.mGamma());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hGammaRxy"), v0.vx(), v0.vy());
  }

  void process(aod::EMReducedEvents::iterator const& collision, aod::V0Photons const& V0Photons, aod::V0Legs const&)
  {
    registry.fill(HIST("hEventCounter"), 1.0); // all
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 2.0); // FT0VX i.e. FT0and

    if (collision.numContrib() < 0.5) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 3.0); // Ncontrib > 0

    if (abs(collision.posZ()) > 10.0) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 4.0); //|Zvtx| < 10 cm

    int ng_nonamb = 0;
    int ng_amb = 0;
    for (auto& g : V0Photons) {
      auto pos = g.posTrack_as<aod::V0Legs>();
      auto ele = g.negTrack_as<aod::V0Legs>();

      if (ele.isAmbTrack() || pos.isAmbTrack()) {
        fillHistosV0<1>(g);
        ng_amb++;
      } else {
        fillHistosV0<0>(g);
        ng_nonamb++;
      }
      for (auto& leg : {pos, ele}) {
        if (leg.isAmbTrack())
          fillHistosLeg<1>(leg);
        else
          fillHistosLeg<0>(leg);
      }

    } // end of v0 loop
    registry.fill(HIST("NonAmbV0/hNgamma"), ng_nonamb);
    registry.fill(HIST("AmbV0/hNgamma"), ng_amb);
  } // end of processSame

  void processDummy(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Photons const& v0photons)
  {
    // do nothing
  }

  PROCESS_SWITCH(PCMQC, process, "run PCM QC", true);
  PROCESS_SWITCH(PCMQC, processDummy, "Dummy function", false);
};

struct Pi0EtaToGammaGamma {
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
      {"hV0Pt", "pT", {HistType::kTH1F, {{1000, 0.0, 10}}}},
      {"hV0LegTPCNcl", "TPC Ncl", {HistType::kTH1F, {{161, -0.5, 160.5}}}},
      {"hV0LegTPCNcr", "TPC Ncr", {HistType::kTH1F, {{161, -0.5, 160.5}}}},
      {"hNgamma", "number of photon conversion per event", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"hV0EtaPhi", "#eta vs. #varphi", {HistType::kTH2F, {{180, 0, TMath::TwoPi()}, {40, -2.0f, 2.0f}}}},
      {"hV0AP", "AP", {HistType::kTH2F, {{200, -1, +1}, {300, 0, 0.3f}}}},
      {"h2MggPt_Same", "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", {HistType::kTH2F, {{400, 0.0, 0.8}, {100, 0.0, 10.}}}},
      {"h2MggPt_Mixed", "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", {HistType::kTH2F, {{400, 0.0, 0.8}, {100, 0.0, 10.}}}},
    },
  };

  void init(InitContext& context)
  {
    addhistograms();
  }

  void addhistograms()
  {
    const TString eventtype[2] = {"Same", "Mixed"};
    const TString pairtype[2] = {"NonAmb", "Amb"};
    for (int ie = 0; ie < 2; ie++) {
      for (int ip = 0; ip < 2; ip++) {
        registry.add(Form("h2MggPt_%s_%s", eventtype[ie].Data(), pairtype[ip].Data()), Form("M_{#gamma#gamma} vs. p_{T} %s %s;m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", eventtype[ie].Data(), pairtype[ip].Data()), HistType::kTH2F, {{400, 0.0, 0.8}, {100, 0.0, 10.}}, true);
      }
    }
  }

  void processSame(aod::EMReducedEvents::iterator const& collision, aod::V0Photons const& V0Photons, aod::V0Legs const&)
  {
    registry.fill(HIST("hEventCounter"), 1.0); // all
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 2.0); // FT0VX i.e. FT0and

    if (collision.numContrib() < 0.5) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 3.0); // Ncontrib > 0

    if (abs(collision.posZ()) > 10.0) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 4.0); //|Zvtx| < 10 cm
    registry.fill(HIST("hNgamma"), V0Photons.size());

    for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(V0Photons, V0Photons))) {
      ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), g1.mGamma());
      ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mGamma());
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      registry.fill(HIST("h2MggPt_Same"), v12.M(), v12.Pt());

      auto ele1 = g1.negTrack_as<aod::V0Legs>(); // Obtain the corresponding track
      auto pos1 = g1.posTrack_as<aod::V0Legs>(); // Obtain the corresponding track
      auto ele2 = g2.negTrack_as<aod::V0Legs>(); // Obtain the corresponding track
      auto pos2 = g2.posTrack_as<aod::V0Legs>(); // Obtain the corresponding track

      if ((!ele1.isAmbTrack() && !pos1.isAmbTrack()) && (!ele2.isAmbTrack() && !pos2.isAmbTrack())) {
        registry.fill(HIST("h2MggPt_Same_NonAmb"), v12.M(), v12.Pt());
      } else {
        registry.fill(HIST("h2MggPt_Same_Amb"), v12.M(), v12.Pt());
      }

    } // end of combination

  } // end of processSame

  Configurable<int> ndepth{"ndepth", 10, "depth of event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType colBinning{{ConfVtxBins}, true};

  Filter collisionFilter = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true;
  using MyFilteredCollisions = soa::Filtered<aod::EMReducedEvents>;
  SameKindPair<MyFilteredCollisions, aod::V0Photons, BinningType> pair{colBinning, ndepth, -1}; // indicates that ndepth events should be mixed and under/overflow (-1) to be ignored
  void processMixed(MyFilteredCollisions const& collisions, aod::V0Photons const& v0photons, aod::V0Legs const&)
  {
    for (auto& [coll1, gammas_coll1, coll2, gammas_coll2] : pair) {
      // LOGF(info, "Mixed event collisions: (%d, %d)", coll1.globalIndex(), coll2.globalIndex());

      for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(gammas_coll1, gammas_coll2))) {
        // LOGF(info, "Mixed event photon pair: (%d, %d) from events (%d, %d), photon event: (%d, %d)", g1.index(), g2.index(), coll1.index(), coll2.index(), g1.collision().index(), g2.collision().index());

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), g1.mGamma());
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mGamma());
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        registry.fill(HIST("h2MggPt_Mixed"), v12.M(), v12.Pt());

        auto ele1 = g1.negTrack_as<aod::V0Legs>(); // Obtain the corresponding track
        auto pos1 = g1.posTrack_as<aod::V0Legs>(); // Obtain the corresponding track
        auto ele2 = g2.negTrack_as<aod::V0Legs>(); // Obtain the corresponding track
        auto pos2 = g2.posTrack_as<aod::V0Legs>(); // Obtain the corresponding track

        if ((!ele1.isAmbTrack() && !pos1.isAmbTrack()) && (!ele2.isAmbTrack() && !pos2.isAmbTrack())) {
          registry.fill(HIST("h2MggPt_Mixed_NonAmb"), v12.M(), v12.Pt());
        } else {
          registry.fill(HIST("h2MggPt_Mixed_Amb"), v12.M(), v12.Pt());
        }
      }

    } // end of combination

  } // end of processMixed

  void processDummy(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Photons const& v0photons)
  {
    // do nothing
  }

  PROCESS_SWITCH(Pi0EtaToGammaGamma, processSame, "enable pairing in same event", true);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processMixed, "enable pairing in mixed event", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PCMQC>(cfgc, TaskName{"pcm-qc"}),
    adaptAnalysisTask<Pi0EtaToGammaGamma>(cfgc, TaskName{"pi0eta-to-gammagamma"}),
  };
}
