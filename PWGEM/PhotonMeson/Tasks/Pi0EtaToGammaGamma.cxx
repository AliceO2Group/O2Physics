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

#include <cstring>

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
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

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

struct Pi0EtaToGammaGamma {
  enum PairType {
    kPCMPCM = 0,
    kPHOSPHOS = 1,
    kEMCEMC = 2,
    kPCMPHOS = 3,
    kPCMEMC = 4,
    kPHOSEMC = 5,
  };

  Filter collisionFilter = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true;
  using MyFilteredCollisions = soa::Filtered<aod::EMReducedEvents>;

  HistogramRegistry registry{"Pi0EtaToGammaGamma"};

  Configurable<bool> useRotation{"useRotation", 0, "use rotation method for EMC-EMC background estimation"};
  Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};
  Configurable<std::string> fConfigEMCCuts{"fConfigEMCCuts", "custom,standard", "Comma separated list of EMCal photon cuts"};

  // Configurable for EMCal cuts
  Configurable<float> EMC_minTime{"EMC_minTime", -20., "Minimum cluster time for EMCal time cut"};
  Configurable<float> EMC_maxTime{"EMC_maxTime", +25., "Maximum cluster time for EMCal time cut"};
  Configurable<float> EMC_minM02{"EMC_minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
  Configurable<float> EMC_maxM02{"EMC_maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
  Configurable<float> EMC_minE{"EMC_minE", 0.7, "Minimum cluster energy for EMCal energy cut"};
  Configurable<int> EMC_minNCell{"EMC_minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
  Configurable<std::vector<float>> EMC_TM_Eta{"EMC_TM_Eta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
  Configurable<std::vector<float>> EMC_TM_Phi{"EMC_TM_Phi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
  Configurable<float> EMC_Eoverp{"EMC_Eoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
  Configurable<bool> EMC_UseExoticCut{"EMC_UseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};

  std::vector<EMCPhotonCut> fEMCCuts;

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
      registry.add(Form("%s/h2MggPt_Same", pairnames[i].data()), "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", HistType::kTH2F, {{400, 0, 0.8}, {400, 0.0f, 40}}, true);
      registry.add(Form("%s/h2MggPt_Mixed", pairnames[i].data()), "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", HistType::kTH2F, {{400, 0, 0.8}, {400, 0.0f, 40}}, true);
    }
    registry.add("EMCEMC/h2MggPt_Rotated", "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/#it{c}^{2});p_{T,#gamma#gamma} (GeV/#it{c})", HistType::kTH2F, {{400, 0, 0.8}, {400, 0.0f, 40}}, true);
  }

  void DefineCuts()
  {
    const float a = EMC_TM_Eta->at(0);
    const float b = EMC_TM_Eta->at(1);
    const float c = EMC_TM_Eta->at(2);

    const float d = EMC_TM_Phi->at(0);
    const float e = EMC_TM_Phi->at(1);
    const float f = EMC_TM_Phi->at(2);
    TString cutNamesStr = fConfigEMCCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        if (std::strcmp(cutname, "custom") == 0) {
          EMCPhotonCut custom_cut = EMCPhotonCut(cutname, cutname);
          custom_cut.SetMinE(EMC_minE);
          custom_cut.SetMinNCell(EMC_minNCell);
          custom_cut.SetM02Range(EMC_minM02, EMC_maxM02);
          custom_cut.SetTimeRange(EMC_minTime, EMC_maxTime);

          custom_cut.SetTrackMatchingEta([&a, &b, &c](float pT) {
            return a + pow(pT + b, c);
          });
          custom_cut.SetTrackMatchingPhi([&d, &e, &f](float pT) {
            return d + pow(pT + e, f);
          });
          custom_cut.SetMinEoverP(EMC_Eoverp);
          custom_cut.SetUseExoticCut(EMC_UseExoticCut);
          fEMCCuts.push_back(custom_cut);
        } else {
          fEMCCuts.push_back(*aod::emccuts::GetCut(cutname));
        }
      }
    }
    LOGF(info, "Number of EMCal cuts = %d", fEMCCuts.size());
  }

  Preslice<MyV0Photons> perCollision = aod::v0photon::collisionId;
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

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          registry.fill(HIST(pairnames[itmp]) + HIST("/h2MggPt_Same"), v12.M(), v12.Pt());

        } // end of combination

      } else {
        for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_coll, photons2_coll))) {

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
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

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
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

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
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
      photon3.SetPt(photon.pt());
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

  Filter v0filter = o2::aod::v0photonflag::isCloser == true;
  using MyFilteredV0Photons = soa::Filtered<MyV0Photons>;

  void processPCMPCM(aod::EMReducedEvents const& collisions, MyFilteredCollisions const& filtered_collisions, MyFilteredV0Photons const& v0photons)
  {
    SameEventPairing<PairType::kPCMPCM>(collisions, v0photons, v0photons, perCollision, perCollision);
    MixedEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision, perCollision);
  }

  void processPHOSPHOS(aod::EMReducedEvents const& collisions, MyFilteredCollisions const& filtered_collisions, aod::PHOSClusters const& phosclusters)
  {
    SameEventPairing<PairType::kPHOSPHOS>(collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos);
    MixedEventPairing<PairType::kPHOSPHOS>(filtered_collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos);
  }

  void processEMCEMC(aod::EMReducedEvents const& collisions, MyFilteredCollisions const& filtered_collisions, aod::SkimEMCClusters const& emcclusters)
  {
    if (useRotation) {
      SameEventPairingWithRotation(collisions, emcclusters);
    } else {
      SameEventPairing<PairType::kEMCEMC>(collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc);
    }
    MixedEventPairing<PairType::kEMCEMC>(filtered_collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc);
  }

  void processPCMPHOS(aod::EMReducedEvents const& collisions, MyFilteredCollisions const& filtered_collisions, MyFilteredV0Photons const& v0photons, aod::PHOSClusters const& phosclusters)
  {
    SameEventPairing<PairType::kPCMPHOS>(collisions, v0photons, phosclusters, perCollision, perCollision_phos);
    MixedEventPairing<PairType::kPCMPHOS>(filtered_collisions, v0photons, phosclusters, perCollision, perCollision_phos);
  }

  void processPCMEMC(aod::EMReducedEvents const& collisions, MyFilteredCollisions const& filtered_collisions, MyFilteredV0Photons const& v0photons, aod::SkimEMCClusters const& emcclusters)
  {
    SameEventPairing<PairType::kPCMEMC>(collisions, v0photons, emcclusters, perCollision, perCollision_emc);
    MixedEventPairing<PairType::kPCMEMC>(filtered_collisions, v0photons, emcclusters, perCollision, perCollision_emc);
  }

  void processPHOSEMC(aod::EMReducedEvents const& collisions, MyFilteredCollisions const& filtered_collisions, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters)
  {
    SameEventPairing<PairType::kPHOSEMC>(collisions, phosclusters, emcclusters, perCollision_phos, perCollision_emc);
    MixedEventPairing<PairType::kPHOSEMC>(filtered_collisions, phosclusters, emcclusters, perCollision_phos, perCollision_emc);
  }

  void processDummy(aod::EMReducedEvents::iterator const& collision)
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
    //    adaptAnalysisTask<LabelUniqueV0>(cfgc, TaskName{"label-unique-v0"}),
    adaptAnalysisTask<Pi0EtaToGammaGamma>(cfgc, TaskName{"pi0eta-to-gammagamma"})};
}
