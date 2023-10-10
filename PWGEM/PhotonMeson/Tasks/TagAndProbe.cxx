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
// This code is for data-driven efficiency for photon analyses. tag and probe method
//    Please write to: daiki.sekihata@cern.ch

#include <cstring>
#include <iterator>

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
#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/PairCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::photonpair;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0Recalculation>;
using MyV0Photon = MyV0Photons::iterator;

struct TagAndProbe {

  Configurable<float> CentMin{"CentMin", -1, "min. centrality"};
  Configurable<float> CentMax{"CentMax", 999, "max. centrality"};
  Configurable<std::string> CentEstimator{"CentEstimator", "FT0M", "centrality estimator"};

  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<std::string> fConfigTagPCMCut{"cfgTagPCMCuts", "tag_track", "Tag PCM photon cut. 1 per wagon"};
  Configurable<std::string> fConfigTagPHOSCut{"cfgTagPHOSCuts", "tag", "Tag PHOS photon cuts. 1 per wagon"};
  Configurable<std::string> fConfigTagEMCCut{"cfgTagEMCCuts", "tag", "Tag EMCal photon cuts. 1 per wagon"};
  Configurable<std::string> fConfigProbePCMCuts{"cfgProbePCMCuts", "analysis,analysis_wo_mee,qc,qc_w_mee", "Comma separated list of V0 photon cuts"};
  Configurable<std::string> fConfigProbePHOSCuts{"cfgProbePHOSCuts", "test02,test03", "Comma separated list of PHOS photon cuts"};
  Configurable<std::string> fConfigProbeEMCCuts{"cfgProbeEMCCuts", "standard", "Comma separated list of EMCal photon cuts"};
  Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  THashList* fMainList = new THashList();

  // std::vector<OutputObj<THashList>> fOutputListCent; //this does not work.

  V0PhotonCut fTagPCMCut;
  PHOSPhotonCut fTagPHOSCut;
  EMCPhotonCut fTagEMCCut;

  std::vector<V0PhotonCut> fProbePCMCuts;
  std::vector<PHOSPhotonCut> fProbePHOSCuts;
  std::vector<EMCPhotonCut> fProbeEMCCuts;
  std::vector<PairCut> fPairCuts;

  std::vector<std::string> fPairNames;
  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processPCMPCM")) {
      fPairNames.push_back("PCMPCM");
    }
    if (context.mOptions.get<bool>("processPHOSPHOS")) {
      fPairNames.push_back("PHOSPHOS");
    }
    if (context.mOptions.get<bool>("processEMCEMC")) {
      fPairNames.push_back("EMCEMC");
    }

    DefinePCMCuts();
    DefinePHOSCuts();
    DefineEMCCuts();
    DefinePairCuts();
    addhistograms();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputPair.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Pair")));
  }

  template <typename TTagCut, typename TProbeCuts, typename TPairCuts>
  void add_pair_histograms(THashList* list_pair, const std::string pairname, TTagCut const& tagcut, TProbeCuts const& probecuts, TPairCuts const& paircuts)
  {
    std::string cutname1 = tagcut.GetName();
    for (auto& cut2 : probecuts) {
      std::string cutname2 = cut2.GetName();
      std::string photon_cut_name = cutname1 + "_" + cutname2;
      THashList* list_pair_subsys_photoncut = o2::aod::emphotonhistograms::AddHistClass(list_pair, photon_cut_name.data());

      for (auto& cut3 : paircuts) {
        std::string pair_cut_name = cut3.GetName();
        o2::aod::emphotonhistograms::AddHistClass(list_pair_subsys_photoncut, pair_cut_name.data());
        THashList* list_pair_subsys_paircut = reinterpret_cast<THashList*>(list_pair_subsys_photoncut->FindObject(pair_cut_name.data()));
        o2::aod::emphotonhistograms::DefineHistograms(list_pair_subsys_paircut, "tag_and_probe", pairname.data());
      } // end of cut3 loop pair cut
    }   // end of cut2 loop
  }

  static constexpr std::string_view pairnames[6] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PHOSEMC"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Pair");
    THashList* list_pair = reinterpret_cast<THashList*>(fMainList->FindObject("Pair"));

    // create sub lists first.

    for (auto& pairname : fPairNames) {
      LOGF(info, "Enabled pairs = %s", pairname.data());

      THashList* list_ev_ss = o2::aod::emphotonhistograms::AddHistClass(list_ev, pairname.data());
      o2::aod::emphotonhistograms::DefineHistograms(list_ev_ss, "Event");

      THashList* list_pair_ss = o2::aod::emphotonhistograms::AddHistClass(list_pair, pairname.data());

      if (pairname == "PCMPCM") {
        add_pair_histograms(list_pair_ss, pairname, fTagPCMCut, fProbePCMCuts, fPairCuts);
      }
      if (pairname == "PHOSPHOS") {
        add_pair_histograms(list_pair_ss, pairname, fTagPHOSCut, fProbePHOSCuts, fPairCuts);
      }
      if (pairname == "EMCEMC") {
        add_pair_histograms(list_pair_ss, pairname, fTagEMCCut, fProbeEMCCuts, fPairCuts);
      }

    } // end of pair name loop
  }

  void DefinePCMCuts()
  {
    fTagPCMCut = *pcmcuts::GetCut(fConfigTagPCMCut.value.data());

    TString cutNamesStr = fConfigProbePCMCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fProbePCMCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of PCM cuts = %d", fProbePCMCuts.size());
  }
  void DefinePHOSCuts()
  {
    fTagPHOSCut = *phoscuts::GetCut(fConfigTagPHOSCut.value.data());

    TString cutNamesStr = fConfigProbePHOSCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fProbePHOSCuts.push_back(*phoscuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of PHOS cuts = %d", fProbePHOSCuts.size());
  }

  void DefineEMCCuts()
  {
    fTagEMCCut = *emccuts::GetCut(fConfigTagEMCCut.value.data());

    TString cutNamesStr = fConfigProbeEMCCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fProbeEMCCuts.push_back(*aod::emccuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of EMCal cuts = %d", fProbeEMCCuts.size());
  }

  void DefinePairCuts()
  {
    TString cutNamesStr = "nocut";
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPairCuts.push_back(*paircuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of Pair cuts = %d", fPairCuts.size());
  }

  Preslice<MyV0Photons> perCollision = aod::v0photon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TTagCut, typename TProbeCuts, typename TPairCuts, typename TLegs, typename TEMCMTs>
  void SameEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TTagCut const& tagcut, TProbeCuts const& probecuts, TPairCuts const& paircuts, TLegs const& legs, TEMCMTs const& emcmatchedtracks)
  {
    THashList* list_ev_pair = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data()));
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

    for (auto& collision : collisions) {
      if ((pairtype == PairType::kPHOSPHOS || pairtype == PairType::kPCMPHOS) && !collision.isPHOSCPVreadout()) {
        continue;
      }
      if ((pairtype == PairType::kEMCEMC || pairtype == PairType::kPCMEMC) && !collision.isEMCreadout()) {
        continue;
      }

      reinterpret_cast<TH1F*>(list_ev_pair->FindObject("hZvtx_before"))->Fill(collision.posZ());
      reinterpret_cast<TH1F*>(list_ev_pair->FindObject("hCollisionCounter"))->Fill(1.0); // all
      if (!collision.sel8()) {
        continue;
      }
      reinterpret_cast<TH1F*>(list_ev_pair->FindObject("hCollisionCounter"))->Fill(2.0); // FT0VX i.e. FT0and

      if (collision.numContrib() < 0.5) {
        continue;
      }
      reinterpret_cast<TH1F*>(list_ev_pair->FindObject("hCollisionCounter"))->Fill(3.0); // Ncontrib > 0

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      reinterpret_cast<TH1F*>(list_ev_pair->FindObject("hZvtx_after"))->Fill(collision.posZ());
      reinterpret_cast<TH1F*>(list_ev_pair->FindObject("hCollisionCounter"))->Fill(4.0); // |Zvtx| < 10 cm
      o2::aod::emphotonhistograms::FillHistClass<EMHistType::kEvent>(list_ev_pair, "", collision);

      auto photons1_coll = photons1.sliceBy(perCollision1, collision.collisionId());
      auto photons2_coll = photons2.sliceBy(perCollision2, collision.collisionId());

      for (auto& g1 : photons1_coll) {

        if constexpr (pairtype == PairType::kPCMPCM) {
          if (!tagcut.template IsSelected<TLegs>(g1)) {
            continue;
          }
        } else if constexpr (pairtype == PairType::kPHOSPHOS) {
          if (!tagcut.template IsSelected<int>(g1)) {
            continue;
          }
        } else if constexpr (pairtype == PairType::kEMCEMC) {
          if (!tagcut.template IsSelected<int>(g1)) {
            continue;
          }
        }

        for (auto& g2 : photons2_coll) {
          if ((pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) && (g1.globalIndex() == g2.globalIndex())) {
            continue;
          }

          for (auto& paircut : paircuts) {
            if (!paircut.IsSelected(g1, g2)) {
              continue;
            }

            for (auto& probecut : probecuts) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }
              reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", tagcut.GetName(), probecut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Probe_Same"))->Fill(v12.M(), v2.Pt());

              if constexpr (pairtype == PairType::kPCMPCM) {
                if (!probecut.template IsSelected<TLegs>(g2)) {
                  continue;
                }
              } else if constexpr (pairtype == PairType::kPHOSPHOS) {
                if (!probecut.template IsSelected<int>(g2)) {
                  continue;
                }
              } else if constexpr (pairtype == PairType::kEMCEMC) {
                if (!probecut.template IsSelected<int>(g2)) {
                  continue;
                }
              }

              reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", tagcut.GetName(), probecut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_PassingProbe_Same"))->Fill(v12.M(), v2.Pt());

              if constexpr (pairtype == PairType::kEMCEMC) {
                RotationBackground<aod::SkimEMCClusters>(v12, v1, v2, photons2_coll, g1.globalIndex(), g2.globalIndex(), probecut, paircut, emcmatchedtracks);
              }
            } // end of probe cut loop
          }   // end of pair cut loop
        }     // end of g2 loop
      }       // end of g1 loop
    }         // end of collision loop
  }

  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfMultBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.f, 20.0f, 30.f, 40.0f, 50.f, 60.0f, 70.f, 80.0f, 90.f, 100.0f, 999.0f}, "Mixing bins - centrality bins"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningType colBinning{{ConfVtxBins, ConfMultBins}, true};

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TTagCut, typename TProbeCuts, typename TPairCuts, typename TLegs, typename TEMCMTs>
  void MixedEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TTagCut const& tagcut, TProbeCuts const& probecuts, TPairCuts const& paircuts, TLegs const& legs, TEMCMTs const& emcmatchedtracks)
  {
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

    // LOGF(info, "Number of collisions after filtering: %d", collisions.size());
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.

      // LOGF(info, "Mixed event collisionId: (%d, %d) , ngpcm: (%d, %d), ngphos: (%d, %d), ngemc: (%d, %d)",
      //     collision1.collisionId(), collision2.collisionId(), collision1.ngpcm(), collision2.ngpcm(), collision1.ngphos(), collision2.ngphos(), collision1.ngemc(), collision2.ngemc());

      auto photons_coll1 = photons1.sliceBy(perCollision1, collision1.collisionId());
      auto photons_coll2 = photons2.sliceBy(perCollision2, collision2.collisionId());
      // LOGF(info, "collision1: posZ = %f, numContrib = %d , sel8 = %d | collision2: posZ = %f, numContrib = %d , sel8 = %d",
      //     collision1.posZ(), collision1.numContrib(), collision1.sel8(), collision2.posZ(), collision2.numContrib(), collision2.sel8());

      for (auto& g1 : photons_coll1) {
        if constexpr (pairtype == PairType::kPCMPCM) {
          if (!tagcut.template IsSelected<TLegs>(g1)) {
            continue;
          }
        } else if constexpr (pairtype == PairType::kPHOSPHOS) {
          if (!tagcut.template IsSelected<int>(g1)) {
            continue;
          }
        } else if constexpr (pairtype == PairType::kEMCEMC) {
          if (!tagcut.template IsSelected<int>(g1)) {
            continue;
          }
        }
        for (auto& g2 : photons_coll2) {
          // LOGF(info, "Mixed event photon pair: (%d, %d) from events (%d, %d), photon event: (%d, %d)", g1.index(), g2.index(), collision1.index(), collision2.index(), g1.collisionId(), g2.collisionId());

          for (auto& paircut : paircuts) {
            if (!paircut.IsSelected(g1, g2)) {
              continue;
            }
            for (auto& probecut : probecuts) {

              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }
              reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", tagcut.GetName(), probecut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Probe_Mixed"))->Fill(v12.M(), v2.Pt());

              if constexpr (pairtype == PairType::kPCMPCM) {
                if (!probecut.template IsSelected<TLegs>(g2)) {
                  continue;
                }
              } else if constexpr (pairtype == PairType::kPHOSPHOS) {
                if (!probecut.template IsSelected<int>(g2)) {
                  continue;
                }
              } else if constexpr (pairtype == PairType::kEMCEMC) {
                if (!probecut.template IsSelected<int>(g2)) {
                  continue;
                }
              }

              reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", tagcut.GetName(), probecut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_PassingProbe_Mixed"))->Fill(v12.M(), v2.Pt());

            } // end of probe cut loop
          }   // end of pair cut loop
        }     // end of g2 loop
      }       // end of g1 loop
    }         // end of different collision combinations
  }

  /// \brief Calculate background (using rotation background method only for EMCal!)
  template <typename TPhotons>
  void RotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, unsigned int ig1, unsigned int ig2, EMCPhotonCut const& cut, PairCut const& paircut, SkimEMCMTs const& emcmatchedtracks)
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
      if (!cut.template IsSelected<aod::SkimEMCMTs>(photon)) {
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

      if (openingAngle1 > minOpenAngle) {
        reinterpret_cast<TH2F*>(fMainList->FindObject("Pair")->FindObject("EMCEMC")->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Same_RotatedBkg"))->Fill(mother1.M(), mother1.Pt());
      }
      if (openingAngle2 > minOpenAngle) {
        reinterpret_cast<TH2F*>(fMainList->FindObject("Pair")->FindObject("EMCEMC")->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Same_RotatedBkg"))->Fill(mother2.M(), mother2.Pt());
      }
    }
  }

  Partition<aod::EMReducedEvents> grouped_collisions = CentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < CentMax; // this goes to same event.
  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true && CentMin < o2::aod::cent::centFT0M&& o2::aod::cent::centFT0M < CentMax;
  Filter collisionFilter_subsys = (o2::aod::emreducedevent::ngpcm >= 1) || (o2::aod::emreducedevent::ngphos >= 1) || (o2::aod::emreducedevent::ngemc >= 1);
  using MyFilteredCollisions = soa::Filtered<aod::EMReducedEvents>; // this goes to mixed event.

  void processPCMPCM(aod::EMReducedEvents const& collisions, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::V0Legs const& legs)
  {
    SameEventPairing<PairType::kPCMPCM>(grouped_collisions, v0photons, v0photons, perCollision, perCollision, fTagPCMCut, fProbePCMCuts, fPairCuts, legs, nullptr);
    MixedEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision, perCollision, fTagPCMCut, fProbePCMCuts, fPairCuts, legs, nullptr);
  }

  void processPHOSPHOS(aod::EMReducedEvents const& collisions, MyFilteredCollisions const& filtered_collisions, aod::PHOSClusters const& phosclusters)
  {
    SameEventPairing<PairType::kPHOSPHOS>(grouped_collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos, fTagPHOSCut, fProbePHOSCuts, fPairCuts, nullptr, nullptr);
    MixedEventPairing<PairType::kPHOSPHOS>(filtered_collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos, fTagPHOSCut, fProbePHOSCuts, fPairCuts, nullptr, nullptr);
  }

  void processEMCEMC(aod::EMReducedEvents const& collisions, MyFilteredCollisions const& filtered_collisions, aod::SkimEMCClusters const& emcclusters, aod::SkimEMCMTs const& emcmatchedtracks)
  {
    SameEventPairing<PairType::kEMCEMC>(grouped_collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc, fTagEMCCut, fProbeEMCCuts, fPairCuts, nullptr, emcmatchedtracks);
    MixedEventPairing<PairType::kEMCEMC>(filtered_collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc, fTagEMCCut, fProbeEMCCuts, fPairCuts, nullptr, emcmatchedtracks);
  }

  void processDummy(aod::EMReducedEvents::iterator const& collision) {}

  PROCESS_SWITCH(TagAndProbe, processPCMPCM, "tag and probe for PCM", false);
  PROCESS_SWITCH(TagAndProbe, processPHOSPHOS, "tag and probe for PHOS", false);
  PROCESS_SWITCH(TagAndProbe, processEMCEMC, "tag and probe for EMCal", false);
  PROCESS_SWITCH(TagAndProbe, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TagAndProbe>(cfgc, TaskName{"tag-and-probe"})};
}
