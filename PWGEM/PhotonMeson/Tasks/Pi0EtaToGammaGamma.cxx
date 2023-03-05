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
#include "TString.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzRotation.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"
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
//#include "DetectorsVertexing/DCAFitterN.h"
//#include "DetectorsBase/Propagator.h"
//#include "DetectorsBase/GeometryManager.h"
//#include "DataFormatsParameters/GRPObject.h"
//#include "DataFormatsParameters/GRPMagField.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

struct PCMQC {
  HistogramRegistry registry{
    "PCMQC",
    {
      {"hCollisionCounter", "hCollisionCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
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
      registry.add(Form("%sTrack/hEtaPhiPosDaughter", tracktype[i].Data()), "#eta vs. #varphi positive daughter;", HistType::kTH2F, {{40, -2.0f, 2.0f}, {180, 0, TMath::TwoPi()}});
      registry.add(Form("%sTrack/hEtaPhiNegDaughter", tracktype[i].Data()), "#eta vs. #varphi negative daughter;", HistType::kTH2F, {{40, -2.0f, 2.0f}, {180, 0, TMath::TwoPi()}});
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

  Preslice<aod::V0Datas> perCollision = aod::v0data::collisionId;
  Preslice<aod::V0DaughterTracks> perV0 = aod::v0data::v0Id;
  void processQC(aod::EMReducedEvents::iterator const& collision, aod::V0Datas const& v0photons, aod::V0DaughterTracks const& v0daughters)
  {
    registry.fill(HIST("hCollisionCounter"), 1.0); // all
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hCollisionCounter"), 2.0); // FT0VX i.e. FT0and

    if (collision.numContrib() < 0.5) {
      return;
    }
    registry.fill(HIST("hCollisionCounter"), 3.0); // Ncontrib > 0

    if (abs(collision.posZ()) > 10.0) {
      return;
    }
    registry.fill(HIST("hCollisionCounter"), 4.0); //|Zvtx| < 10 cm
    auto V0Photons_coll = v0photons.sliceBy(perCollision, collision.collisionId());

    // int ng_nonamb = 0;
    // int ng_amb = 0;
    for (auto& g : V0Photons_coll) {
      auto V0Daughters = v0daughters.sliceBy(perV0, g.globalIndex());
      for (auto& daughter : V0Daughters) {
        if (daughter.positivelyCharged() == true) {
          registry.fill(HIST("NonAmbTrack/hEtaPhiPosDaughter"), daughter.eta(), daughter.phi());
        } else {
          registry.fill(HIST("NonAmbTrack/hEtaPhiNegDaughter"), daughter.eta(), daughter.phi());
        }
      } // end of loop over V0 daughters
      // auto pos = g.posTrack_as<aod::V0DaughterTracks>();
      // auto ele = g.negTrack_as<aod::V0DaughterTracks>();

      // if (ele.isAmbTrack() || pos.isAmbTrack()) {
      //   fillHistosV0<1>(g);
      //   ng_amb++;
      // } else {
      //   fillHistosV0<0>(g);
      //   ng_nonamb++;
      // }
      // for (auto& leg : {pos, ele}) {
      //   if (leg.isAmbTrack())
      //     fillHistosLeg<1>(leg);
      //   else
      //     fillHistosLeg<0>(leg);
      // }

    } // end of v0 loop
    // registry.fill(HIST("NonAmbV0/hNgamma"), ng_nonamb);
    // registry.fill(HIST("AmbV0/hNgamma"), ng_amb);
  } // end of processSame

  void processDummy(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const& v0photons)
  {
    // do nothing
  }

  PROCESS_SWITCH(PCMQC, processQC, "run PCM QC", true);
  PROCESS_SWITCH(PCMQC, processDummy, "Dummy function", false);
};

struct Pi0EtaToGammaGamma {
  enum PairType {
    kPCMPCM = 0,
    kPHOSPHOS = 1,
    kEMCEMC = 2,
    kPCMPHOS = 3,
    kPCMEMC = 4,
    kPHOSEMC = 5,
  };

  Partition<aod::EMReducedEvents> goodEventsPCM = o2::aod::emreducedevent::ngpcm > 0; // && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true && o2::aod::emreducedevent::ngpcm > 0;

  HistogramRegistry registry{"Pi0EtaToGammaGamma"};

  Configurable<bool> useRotation{"useRotation", 0, "use rotation method for EMC-EMC background estimation"};
  Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};

  void init(InitContext& context)
  {
    addhistograms();
  }

  static constexpr std::string_view pairnames[6] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PHOSEMC"};
  void addhistograms()
  {
    for (int i = 0; i < 6; i++) {
      registry.add(Form("%s/hCollisionCounter", pairnames[i].data()), "Collision counter", HistType::kTH1F, {{5, 0.5f, 5.5f}});
      registry.add(Form("%s/hNgamma1", pairnames[i].data()), "Number of #gamma1 candidates per collision", HistType::kTH1F, {{101, -0.5f, 100.5f}});
      registry.add(Form("%s/hNgamma2", pairnames[i].data()), "Number of #gamma2 candidates per collision", HistType::kTH1F, {{101, -0.5f, 100.5f}});
      registry.add(Form("%s/h2MggPt_Same", pairnames[i].data()), "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", HistType::kTH2F, {{400, 0, 0.8}, {200, 0.0f, 40}});
      registry.add(Form("%s/h2MggPt_Mixed", pairnames[i].data()), "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", HistType::kTH2F, {{400, 0, 0.8}, {200, 0.0f, 40}});
    }
    registry.add("EMCEMC/h2MggPt_Rotated", "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/#it{c}^{2});p_{T,#gamma#gamma} (GeV/#it{c})", HistType::kTH2F, {{400, 0, 0.8}, {200, 0.0f, 40}});
  }

  Preslice<aod::V0Datas> perCollision = aod::v0data::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2>
  void SameEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2)
  {
    constexpr int itmp = pairtype;
    for (auto& collision : collisions) {
      registry.fill(HIST(pairnames[itmp]) + HIST("/hCollisionCounter"), 1.0); // all
      if (!collision.sel8()) {
        continue;
      }
      registry.fill(HIST(pairnames[itmp]) + HIST("/hCollisionCounter"), 2.0); // FT0VX i.e. FT0and

      if (collision.numContrib() < 0.5) {
        continue;
      }
      registry.fill(HIST(pairnames[itmp]) + HIST("/hCollisionCounter"), 3.0); // Ncontrib > 0

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      registry.fill(HIST(pairnames[itmp]) + HIST("/hCollisionCounter"), 4.0); //|Zvtx| < 10 cm

      auto photons1_coll = photons1.sliceBy(perCollision1, collision.collisionId());
      auto photons2_coll = photons2.sliceBy(perCollision2, collision.collisionId());

      registry.fill(HIST(pairnames[itmp]) + HIST("/hNgamma1"), photons1_coll.size());
      registry.fill(HIST(pairnames[itmp]) + HIST("/hNgamma2"), photons2_coll.size());
      // LOGF(info, "Number of photon candidates in a collision: %d", photons1_coll.size());

      if (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) {
        for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_coll, photons2_coll))) {

          ROOT::Math::PtEtaPhiMVector v1;
          ROOT::Math::PtEtaPhiMVector v2;

          if constexpr (pairtype == PairType::kPCMPCM) {
            v1.SetPt(g1.pt());
            v1.SetEta(g1.eta());
            v1.SetPhi(g1.phi());
            v1.SetM(g1.mGamma());
            v2.SetPt(g2.pt());
            v2.SetEta(g2.eta());
            v2.SetPhi(g2.phi());
            v2.SetM(g2.mGamma());
          } else if constexpr (pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) {
            v1.SetPt(g1.pt(0.));
            v1.SetEta(g1.eta());
            v1.SetPhi(g1.phi());
            v1.SetM(0.);
            v2.SetPt(g2.pt(0.));
            v2.SetEta(g2.eta());
            v2.SetPhi(g2.phi());
            v2.SetM(0.);
          }

          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          registry.fill(HIST(pairnames[itmp]) + HIST("/h2MggPt_Same"), v12.M(), v12.Pt());

        } // end of combination

      } else {
        for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_coll, photons2_coll))) {

          ROOT::Math::PtEtaPhiMVector v1;
          ROOT::Math::PtEtaPhiMVector v2;

          if constexpr (pairtype == PairType::kPCMPHOS || pairtype == PairType::kPCMEMC) {
            v1.SetPt(g1.pt());
            v1.SetEta(g1.eta());
            v1.SetPhi(g1.phi());
            v1.SetM(g1.mGamma());
            v2.SetPt(g2.pt(0.));
            v2.SetEta(g2.eta());
            v2.SetPhi(g2.phi());
            v2.SetM(0.);
          } else if constexpr (pairtype == PairType::kPHOSEMC) {
            v1.SetPt(g1.pt(0.));
            v1.SetEta(g1.eta());
            v1.SetPhi(g1.phi());
            v1.SetM(0.);
            v2.SetPt(g2.pt(0.));
            v2.SetEta(g2.eta());
            v2.SetPhi(g2.phi());
            v2.SetM(0.);
          }

          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          registry.fill(HIST(pairnames[itmp]) + HIST("/h2MggPt_Same"), v12.M(), v12.Pt());

        } // end of combination
      }

    } // end of collision loop
  }

  void SameEventPairingWithRotation(aod::EMReducedEvents const& collisions, aod::SkimEMCClusters const& photons)
  {
    for (auto& collision : collisions) {
      registry.fill(HIST("EMCEMC/hCollisionCounter"), 1.0); // all
      if (!collision.sel8()) {
        continue;
      }
      registry.fill(HIST("EMCEMC/hCollisionCounter"), 2.0); // FT0VX i.e. FT0and

      if (collision.numContrib() < 0.5) {
        continue;
      }
      registry.fill(HIST("EMCEMC/hCollisionCounter"), 3.0); // Ncontrib > 0

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      registry.fill(HIST("EMCEMC/hCollisionCounter"), 4.0); //|Zvtx| < 10 cm

      auto photons_coll = photons.sliceBy(perCollision_emc, collision.collisionId());

      registry.fill(HIST("EMCEMC/hNgamma1"), photons_coll.size());
      registry.fill(HIST("EMCEMC/hNgamma2"), photons_coll.size());
      // LOGF(info, "Number of photon candidates in a collision: %d", photons_coll.size());

      for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons_coll, photons_coll))) {

        ROOT::Math::PtEtaPhiMVector v1;
        ROOT::Math::PtEtaPhiMVector v2;
        v1.SetPt(g1.pt(0.));
        v1.SetEta(g1.eta());
        v1.SetPhi(g1.phi());
        v1.SetM(0.);
        v2.SetPt(g2.pt(0.));
        v2.SetEta(g2.eta());
        v2.SetPhi(g2.phi());
        v2.SetM(0.);

        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        registry.fill(HIST("EMCEMC/h2MggPt_Same"), v12.M(), v12.Pt());

        RotationBackground<aod::SkimEMCClusters, aod::EMReducedEvent>(v12, v1, v2, photons_coll, collision, g1.globalIndex(), g2.globalIndex());
      } // end of combination
    }   // end of collision loop
  }

  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType colBinning{{ConfVtxBins}, true};

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2>
  void MixedEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2)
  {
    constexpr int itmp = pairtype;
    // LOGF(info, "Number of collisions after filtering: %d", collisions.size());
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) {
      // LOGF(info, "Mixed event collisionId: (%d, %d)", collision1.collisionId(), collision2.collisionId());
      auto photons_coll1 = photons1.sliceBy(perCollision1, collision1.collisionId());
      auto photons_coll2 = photons2.sliceBy(perCollision2, collision2.collisionId());
      // LOGF(info, "collision1: posZ = %f, numContrib = %d , sel8 = %d, ngpcm = %d , ngphos = %d , ngemc = %d", collision1.posZ(), collision1.numContrib(), collision1.sel8(), collision1.ngpcm(), collision1.ngphos(), collision1.ngemc());

      for (auto& [g1, g2] : combinations(soa::CombinationsFullIndexPolicy(photons_coll1, photons_coll2))) {
        // LOGF(info, "Mixed event photon pair: (%d, %d) from events (%d, %d), photon event: (%d, %d)", g1.index(), g2.index(), collision1.index(), collision2.index(), g1.collisionId(), g2.collisionId());

        ROOT::Math::PtEtaPhiMVector v1;
        ROOT::Math::PtEtaPhiMVector v2;

        if constexpr (pairtype == PairType::kPCMPCM) {
          v1.SetPt(g1.pt());
          v1.SetEta(g1.eta());
          v1.SetPhi(g1.phi());
          v1.SetM(g1.mGamma());
          v2.SetPt(g2.pt());
          v2.SetEta(g2.eta());
          v2.SetPhi(g2.phi());
          v2.SetM(g2.mGamma());
        } else if constexpr (pairtype == PairType::kPHOSPHOS) {
          v1.SetPt(g1.pt(0.));
          v1.SetEta(g1.eta());
          v1.SetPhi(g1.phi());
          v1.SetM(0.);
          v2.SetPt(g2.pt(0.));
          v2.SetEta(g2.eta());
          v2.SetPhi(g2.phi());
          v2.SetM(0.);
        } else if constexpr (pairtype == PairType::kEMCEMC) {
          v1.SetPt(g1.pt(0.));
          v1.SetEta(g1.eta());
          v1.SetPhi(g1.phi());
          v1.SetM(0.);
          v2.SetPt(g2.pt(0.));
          v2.SetEta(g2.eta());
          v2.SetPhi(g2.phi());
          v2.SetM(0.);
        } else if constexpr (pairtype == PairType::kPCMPHOS) {
          v1.SetPt(g1.pt());
          v1.SetEta(g1.eta());
          v1.SetPhi(g1.phi());
          v1.SetM(g1.mGamma());
          v2.SetPt(g2.pt(0.));
          v2.SetEta(g2.eta());
          v2.SetPhi(g2.phi());
          v2.SetM(0.);
        } else if constexpr (pairtype == PairType::kPCMEMC) {
          v1.SetPt(g1.pt());
          v1.SetEta(g1.eta());
          v1.SetPhi(g1.phi());
          v1.SetM(g1.mGamma());
          v2.SetPt(g2.pt(0.));
          v2.SetEta(g2.eta());
          v2.SetPhi(g2.phi());
          v2.SetM(0.);
        } else if constexpr (pairtype == PairType::kPHOSEMC) {
          v1.SetPt(g1.pt(0.));
          v1.SetEta(g1.eta());
          v1.SetPhi(g1.phi());
          v1.SetM(0.);
          v2.SetPt(g2.pt(0.));
          v2.SetEta(g2.eta());
          v2.SetPhi(g2.phi());
          v2.SetM(0.);
        }

        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        registry.fill(HIST(pairnames[itmp]) + HIST("/h2MggPt_Mixed"), v12.M(), v12.Pt());

      } // end of different photon combinations

    } // end of different collision combinations
  }

  /// \brief Calculate background (using rotation background method only for EMCal!)
  template <typename TPhotons, typename TEvent>
  void RotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, TEvent const& collision, unsigned int ig1, unsigned int ig2)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photons_coll.size() < 3) {
      return;
    }
    const double rotationAngle = M_PI / 2.0; // rotaion angle 90°
    // ROOT::Math::XYZVector meson3D(meson.Px(), meson.Py(), meson.Pz());
    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), rotationAngle);
    // LOG(info) << "rotationAxis.Angle() = " << rotationAxis.Angle();
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    // ROOT::Math::XYZVector test(photon1.Px(), photon1.Py(), photon1.Pz());
    // ROOT::Math::XYZVector photon1_3D(photon1.Px(), photon1.Py(), photon1.Pz());
    // ROOT::Math::XYZVector photon2_3D(photon2.Px(), photon2.Py(), photon2.Pz());

    // float openingAngleTest = std::acos(photon1_3D.Dot(test) / (std::sqrt(photon1_3D.Mag2()) * std::sqrt(test.Mag2())));
    // LOG(info) << "openingAngleTest before rotation = " << openingAngleTest;
    photon1 = rotationMatrix * photon1;
    photon2 = rotationMatrix * photon2;

    // openingAngleTest = std::acos(photon1_3D.Dot(test) / (std::sqrt(photon1_3D.Mag2()) * std::sqrt(test.Mag2())));
    // LOG(info) << "openingAngleTest = " << openingAngleTest;

    for (auto& photon : photons_coll) {
      if (photon.globalIndex() == ig1 || photon.globalIndex() == ig2) {
        // only combine rotated photons with other photons
        continue;
      }
      ROOT::Math::PtEtaPhiMVector photon3;
      photon3.SetPt(photon.pt(0.));
      photon3.SetEta(photon.eta());
      photon3.SetPhi(photon.phi());
      photon3.SetM(0.);
      // ROOT::Math::XYZVector photon3_3D(photon3.Px(), photon3.Py(), photon3.Pz());
      ROOT::Math::PtEtaPhiMVector mother1 = photon1 + photon3;
      ROOT::Math::PtEtaPhiMVector mother2 = photon2 + photon3;

      float openingAngle1 = std::acos(photon1.Vect().Dot(photon3.Vect()) / (photon1.P() * photon3.P()));
      float openingAngle2 = std::acos(photon2.Vect().Dot(photon3.Vect()) / (photon2.P() * photon3.P()));
      // float openingAngle1_1 = std::acos(photon1_3D.Dot(photon3_3D) / (std::sqrt(photon1_3D.Mag2()) * std::sqrt(photon3_3D.Mag2())));
      // float openingAngle2_2 = std::acos(photon2_3D.Dot(photon3_3D) / (std::sqrt(photon2_3D.Mag2()) * std::sqrt(photon3_3D.Mag2())));
      // LOG(info) << "openingAngle1 = " << openingAngle1;
      // LOG(info) << "openingAngle2 = " << openingAngle2;
      // LOG(info) << "openingAngle1_1 = " << openingAngle1_1;
      // LOG(info) << "openingAngle2_2 = " << openingAngle2_2;

      // Fill histograms
      if (openingAngle1 > minOpenAngle) {
        registry.fill(HIST("EMCEMC/h2MggPt_Rotated"), mother1.M(), mother1.Pt());
      }
      if (openingAngle2 > minOpenAngle) {
        registry.fill(HIST("EMCEMC/h2MggPt_Rotated"), mother2.M(), mother2.Pt());
      }
    }
  }

  // Filter collisionFilter_PCM_mix = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true && o2::aod::emreducedevent::ngpcm > 0;
  // using MyFilteredCollisions_PCM_mix = soa::Filtered<aod::EMReducedEvents>;
  void processPCMPCM(aod::EMReducedEvents const& collisions, aod::V0Datas const& v0photons)
  {
    SameEventPairing<PairType::kPCMPCM>(collisions, v0photons, v0photons, perCollision, perCollision);
    MixedEventPairing<PairType::kPCMPCM>(collisions, v0photons, v0photons, perCollision, perCollision);
  }

  // Filter collisionFilter_phos_mix = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true && o2::aod::emreducedevent::ngphos > 0;
  // using MyFilteredCollisions_phos_mix = soa::Filtered<aod::EMReducedEvents>;
  void processPHOSPHOS(aod::EMReducedEvents const& collisions, aod::PHOSClusters const& phosclusters)
  {
    SameEventPairing<PairType::kPHOSPHOS>(collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos);
    MixedEventPairing<PairType::kPHOSPHOS>(collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos);
  }

  void processEMCEMC(aod::EMReducedEvents const& collisions, aod::SkimEMCClusters const& emcclusters)
  {
    if (useRotation) {
      SameEventPairingWithRotation(collisions, emcclusters);
    } else {
      SameEventPairing<PairType::kEMCEMC>(collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc);
    }
    MixedEventPairing<PairType::kEMCEMC>(collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc);
  }

  // Filter collisionFilter_pcm_phos_mix = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true && o2::aod::emreducedevent::ngpcm > 0 && o2::aod::emreducedevent::ngphos > 0;
  // using MyFilteredCollisions_pcm_phos_mix = soa::Filtered<aod::EMReducedEvents>;
  void processPCMPHOS(aod::EMReducedEvents const& collisions, aod::V0Datas const& v0photons, aod::PHOSClusters const& phosclusters)
  {
    SameEventPairing<PairType::kPCMPHOS>(collisions, v0photons, phosclusters, perCollision, perCollision_phos);
    MixedEventPairing<PairType::kPCMPHOS>(collisions, v0photons, phosclusters, perCollision, perCollision_phos);
  }

  void processPCMEMC(aod::EMReducedEvents const& collisions, aod::V0Datas const& v0photons, aod::SkimEMCClusters const& emcclusters)
  {
    SameEventPairing<PairType::kPCMEMC>(collisions, v0photons, emcclusters, perCollision, perCollision_emc);
    MixedEventPairing<PairType::kPCMEMC>(collisions, v0photons, emcclusters, perCollision, perCollision_emc);
  }

  void processPHOSEMC(aod::EMReducedEvents const& collisions, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters)
  {
    SameEventPairing<PairType::kPHOSEMC>(collisions, phosclusters, emcclusters, perCollision_phos, perCollision_emc);
    MixedEventPairing<PairType::kPHOSEMC>(collisions, phosclusters, emcclusters, perCollision_phos, perCollision_emc);
  }

  void processDummy(soa::Join<aod::Collisions, aod::EvSels> const& collision, aod::V0Datas const& v0photons)
  {
    // do nothing
  }

  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPCMPCM, "pairing PCM-PCM", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPHOSPHOS, "pairing PHOS-PHOS", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processEMCEMC, "pairing EMCal-EMCal", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPCMPHOS, "pairing PCM-PHOS", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPCMEMC, "pairing PCM-EMCal", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPHOSEMC, "pairing PHOS-EMCal", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PCMQC>(cfgc, TaskName{"pcm-qc"}),
    adaptAnalysisTask<Pi0EtaToGammaGamma>(cfgc, TaskName{"pi0eta-to-gammagamma"}),
  };
}
