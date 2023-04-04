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
// This code runs loop over v0 photons for PCM QC.
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
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

// namespace o2::aod
//{
// namespace v0photonflag // flag to distinguish 1 track belongs to 1 V0 or 2 (or more) V0s in a collision.
//{
// DECLARE_SOA_COLUMN(IsCloser, isCloser, bool); //! true if 2 legs of this v0 do not belong to other V0s in a collision or PCA between 2 legs is closer.
// } // namespace v0photonflag
// DECLARE_SOA_TABLE(V0PhotonFlags, "AOD", "V0PHOTONFLAG", v0photonflag::IsCloser);
// using V0PhotonFlag = V0PhotonFlags::iterator;
// } // namespace o2::aod
//
// struct LabelUniqueV0 {
//   Produces<aod::V0PhotonFlags> v0flags;
//
//   Preslice<aod::V0Photons> perCollision = aod::v0photon::collisionId;
//   void process(aod::EMReducedEvents::iterator const& collision, aod::V0Photons const& v0photons, aod::V0Legs const& v0legs)
//   {
//     auto v0photons_coll = v0photons.sliceBy(perCollision, collision.collisionId());
//
//     for (auto& g1 : v0photons_coll) {
//       auto pos1 = g1.posTrack_as<aod::V0Legs>();
//       auto ele1 = g1.negTrack_as<aod::V0Legs>();
//       bool flag = true;
//
//       int posid1 = pos1.trackId(); // unique index to point o2::track
//       int eleid1 = ele1.trackId(); // unique index to point o2::track
//       float pca1 = g1.pca();
//
//       for (auto& g2 : v0photons_coll) {
//         if (g2.index() == g1.index()) {
//           continue;
//         }
//
//         auto pos2 = g2.posTrack_as<aod::V0Legs>();
//         auto ele2 = g2.negTrack_as<aod::V0Legs>();
//         int posid2 = pos2.trackId(); // unique index to point o2::track
//         int eleid2 = ele2.trackId(); // unique index to point o2::track
//         float pca2 = g2.pca();
//
//         if ((posid2 == posid1 || eleid2 == eleid1) && pca1 > pca2) {
//           // LOGF(info, "g1 id = %d , g2 id = %d , posid1 = %d , eleid1 = %d , posid2 = %d , eleid2 = %d , pca1 = %f , pca2 = %f", g1.index(), g2.index(), posid1, eleid1, posid2, eleid2, pca1, pca2);
//           flag = false;
//           break;
//         }
//       }
//       v0flags(flag);
//     }
//   }
// };

using MyV0Photons = soa::Join<aod::V0Photons, aod::V0RecalculationAndKF, aod::V0PhotonFlags>;
using MyV0Photon = MyV0Photons::iterator;

struct PCMQCMC {
  HistogramRegistry registry{
    "PCMQCMC",
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
      registry.add(Form("%sV0/hRadius", pairtype[i].Data()), "V0Radius; radius in Z (cm);radius in XY (cm)", HistType::kTH2F, {{500, -250, 250}, {500, 0.0f, 250.0f}});
      registry.add(Form("%sV0/hRadius_recalc", pairtype[i].Data()), "V0Radius; radius in Z (cm);radius in XY (cm)", HistType::kTH2F, {{500, -250, 250}, {500, 0.0f, 250.0f}});
      registry.add(Form("%sV0/hCosPA", pairtype[i].Data()), "V0CosPA;cosine pointing angle", HistType::kTH1F, {{200, 0.8f, 1.0f}});
      registry.add(Form("%sV0/hPCA", pairtype[i].Data()), "distance between 2 legs; PCA (cm)", HistType::kTH1F, {{100, 0.0f, 10.0f}});
      registry.add(Form("%sV0/hAPplot", pairtype[i].Data()), "AP plot;#alpha;q_{T} (GeV/c)", HistType::kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}});
      registry.add(Form("%sV0/hGammaPsiPair", pairtype[i].Data()), "#psi_{pair} for photon conversion;#psi_{pair} (rad.);m_{ee} (GeV/c^{2})", HistType::kTH2F, {{160, 0.0, TMath::PiOver2()}, {100, 0.0f, 0.1f}});
      registry.add(Form("%sV0/hMassGamma", pairtype[i].Data()), "hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", HistType::kTH2F, {{2000, 0.0f, 200.0f}, {100, 0.0f, 0.1f}});
      registry.add(Form("%sV0/hMassGamma_recalc", pairtype[i].Data()), "recalc. KF hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", HistType::kTH2F, {{2000, 0.0f, 200.0f}, {100, 0.0f, 0.1f}});
      registry.add(Form("%sV0/hGammaRxy", pairtype[i].Data()), "conversion point in XY;V_{x} (cm);V_{y} (cm)", HistType::kTH2F, {{1000, -250.0f, 250.0f}, {1000, -250.0f, 250.0f}});
      registry.add(Form("%sV0/hGammaRxy_recalc", pairtype[i].Data()), "recalc. KF conversion point in XY;V_{x} (cm);V_{y} (cm)", HistType::kTH2F, {{1000, -250.0f, 250.0f}, {1000, -250.0f, 250.0f}});
      registry.add(Form("%sV0/hKFChi2vsR_recalc", pairtype[i].Data()), "recalc. KF conversion point in XY;R_{xy} (cm);KF chi2/NDF", HistType::kTH2F, {{250, 0.0f, 250.0f}, {5000, 0.f, 5000.0f}});
      registry.add(Form("%sV0/hKFChi2vsZ_recalc", pairtype[i].Data()), "recalc. KF conversion point in Z;Z (cm);KF chi2/NDF", HistType::kTH2F, {{500, -250.0f, 250.0f}, {5000, 0.f, 5000.0f}});
    }

    const float rxy[] = {0, 6, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    registry.add("Generated/hCollisionCounter", "Collision counter", HistType::kTH1F, {{5, 0.5f, 5.5f}});
    registry.add("Generated/hGammaRxy", "conversion point in XY MC;V_{x} (cm);V_{y} (cm)", HistType::kTH2F, {{2000, -100.0f, 100.0f}, {2000, -100.0f, 100.0f}});
    registry.add("Generated/hGammaRZ", "conversion point in RZ MC;V_{z} (cm);R_{xy} (cm)", HistType::kTH2F, {{5000, -250.0f, 250.0f}, {1000, 0.f, 100.0f}});

    const int n = sizeof(rxy) / sizeof(rxy[0]);
    for (int i = 0; i < n - 1; i++) {
      float rmin = rxy[i];
      float rmax = rxy[i + 1];
      registry.add(Form("Generated/hConvPhi_Rxy%d_%dcm", static_cast<int>(rmin), static_cast<int>(rmax)), Form("conversion point of #varphi MC in %d < R_{xy} < %d cm;#varphi (rad.);N_{e}", static_cast<int>(rmin), static_cast<int>(rmax)), HistType::kTH1F, {{360, 0.0f, TMath::TwoPi()}});
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
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hRadius_recalc"), v0.recalculatedVtxZ(), v0.recalculatedVtxR());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hCosPA"), v0.cospa());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hPCA"), v0.pca());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hAPplot"), v0.alpha(), v0.qtarm());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hMassGamma"), v0.v0radius(), v0.mGamma());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hMassGamma_recalc"), v0.recalculatedVtxR(), v0.mGamma());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hGammaPsiPair"), v0.psipair(), v0.mGamma());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hGammaRxy"), v0.vx(), v0.vy());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hGammaRxy_recalc"), v0.recalculatedVtxX(), v0.recalculatedVtxY());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hKFChi2vsR_recalc"), v0.recalculatedVtxR(), v0.chiSquareNDF());
    registry.fill(HIST(typenames[mode]) + HIST("V0/") + HIST("hKFChi2vsZ_recalc"), v0.recalculatedVtxZ(), v0.chiSquareNDF());
  }

  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int FindCommonMotherFrom2Prongs(TMCParticle1 const& p1, TMCParticle2 const& p2, const int expected_pdg1, const int expected_pdg2, const int expected_mother_pdg, TMCParticles const& mcparticles)
  {
    if (p1.globalIndex() == p2.globalIndex())
      return -1; // mc particle p1 and p2 is identical. reject.

    if (p1.pdgCode() != expected_pdg1)
      return -1;
    if (p2.pdgCode() != expected_pdg2)
      return -1;

    if (!p1.has_mothers())
      return -1;
    if (!p2.has_mothers())
      return -1;

    // LOGF(info,"original motherid1 = %d , motherid2 = %d", p1.mothersIds()[0], p2.mothersIds()[0]);

    int motherid1 = p1.mothersIds()[0];
    auto mother1 = mcparticles.iteratorAt(motherid1);
    int mother1_pdg = mother1.pdgCode();

    int motherid2 = p2.mothersIds()[0];
    auto mother2 = mcparticles.iteratorAt(motherid2);
    int mother2_pdg = mother2.pdgCode();

    // LOGF(info,"motherid1 = %d , motherid2 = %d", motherid1, motherid2);

    if (motherid1 != motherid2)
      return -1;
    if (mother1_pdg != mother2_pdg)
      return -1;
    if (mother1_pdg != expected_mother_pdg)
      return -1;
    return motherid1;
  }

  Filter v0filter = o2::aod::v0photonflag::isCloser == true;
  using MyFilteredV0Photons = soa::Filtered<MyV0Photons>;
  Preslice<MyV0Photons> perCollision = aod::v0photon::collisionId;

  using MyMCV0Legs = soa::Join<aod::V0Legs, aod::EMMCParticleLabels>;
  void processQCMC(soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels> const& collisions, MyFilteredV0Photons const& v0photons, MyMCV0Legs const& v0legs, aod::EMMCParticles const& mcparticles)
  {
    for (auto& collision : collisions) {
      registry.fill(HIST("hCollisionCounter"), 1.0); // all
      if (!collision.sel8()) {
        continue;
      }
      registry.fill(HIST("hCollisionCounter"), 2.0); // FT0VX i.e. FT0and

      if (collision.numContrib() < 0.5) {
        continue;
      }
      registry.fill(HIST("hCollisionCounter"), 3.0); // Ncontrib > 0

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      registry.fill(HIST("hCollisionCounter"), 4.0); //|Zvtx| < 10 cm
      auto V0Photons_coll = v0photons.sliceBy(perCollision, collision.collisionId());

      int ng_nonamb = 0;
      int ng_amb = 0;
      for (auto& g : V0Photons_coll) {
        auto pos = g.posTrack_as<MyMCV0Legs>();
        auto ele = g.negTrack_as<MyMCV0Legs>();
        auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
        auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

        int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
        if (photonid < 0) { // check swap, true electron is reconstructed as positron and vice versa.
          photonid = FindCommonMotherFrom2Prongs(posmc, elemc, 11, -11, 22, mcparticles);
        }

        if (photonid < 0) {
          continue;
        }

        registry.fill(HIST("NonAmbTrack/hEtaPhiPosDaughter"), pos.eta(), pos.phi());
        registry.fill(HIST("NonAmbTrack/hEtaPhiNegDaughter"), ele.eta(), ele.phi());

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
    }   // end of collision loop
  }     // end of process

  Preslice<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emreducedmceventId;
  void processGen(soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels> const& collisions, aod::EMReducedMCEvents const& mccollisions, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : collisions) {
      registry.fill(HIST("Generated/hCollisionCounter"), 1.0); // all
      if (!collision.sel8()) {
        continue;
      }
      registry.fill(HIST("Generated/hCollisionCounter"), 2.0); // FT0VX i.e. FT0and

      if (collision.numContrib() < 0.5) {
        continue;
      }
      registry.fill(HIST("Generated/hCollisionCounter"), 3.0); // Ncontrib > 0

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      registry.fill(HIST("Generated/hCollisionCounter"), 4.0); //|Zvtx| < 10 cm
      auto mccollision = collision.emreducedmcevent();

      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (auto& mctrack : mctracks_coll) {

        if (IsEleFromPC(mctrack, mcparticles) > 0) {
          float rxy = sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2));
          registry.fill(HIST("Generated/hGammaRZ"), mctrack.vz(), rxy);

          if (abs(mctrack.eta()) > 0.9)
            continue;
          registry.fill(HIST("Generated/hGammaRxy"), mctrack.vx(), mctrack.vy());
          if (rxy < 6) {
            registry.fill(HIST("Generated/hConvPhi_Rxy0_6cm"), mctrack.phi());
          } else if (rxy < 10) {
            registry.fill(HIST("Generated/hConvPhi_Rxy6_10cm"), mctrack.phi());

          } else if (rxy < 20) {
            registry.fill(HIST("Generated/hConvPhi_Rxy10_20cm"), mctrack.phi());

          } else if (rxy < 30) {
            registry.fill(HIST("Generated/hConvPhi_Rxy20_30cm"), mctrack.phi());

          } else if (rxy < 40) {
            registry.fill(HIST("Generated/hConvPhi_Rxy30_40cm"), mctrack.phi());

          } else if (rxy < 50) {
            registry.fill(HIST("Generated/hConvPhi_Rxy40_50cm"), mctrack.phi());

          } else if (rxy < 60) {
            registry.fill(HIST("Generated/hConvPhi_Rxy50_60cm"), mctrack.phi());

          } else if (rxy < 70) {
            registry.fill(HIST("Generated/hConvPhi_Rxy60_70cm"), mctrack.phi());

          } else if (rxy < 80) {
            registry.fill(HIST("Generated/hConvPhi_Rxy70_80cm"), mctrack.phi());

          } else if (rxy < 90) {
            registry.fill(HIST("Generated/hConvPhi_Rxy80_90cm"), mctrack.phi());
          }
        }
      }
    }
  }

  void processDummy(aod::EMReducedEvents::iterator const& collision)
  {
    // do nothing
  }

  PROCESS_SWITCH(PCMQCMC, processQCMC, "run PCM QC in MC", true);
  PROCESS_SWITCH(PCMQCMC, processGen, "run generated information", false);
  PROCESS_SWITCH(PCMQCMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    //    adaptAnalysisTask<LabelUniqueV0>(cfgc, TaskName{"label-unique-v0"}),
    adaptAnalysisTask<PCMQCMC>(cfgc, TaskName{"pcm-qc-mc"})};
}
