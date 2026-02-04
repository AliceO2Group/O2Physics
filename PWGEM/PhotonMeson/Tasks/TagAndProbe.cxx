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

#include "EMPhotonEventCut.h"

#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/PairCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/AxisAngle.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/Vector4Dfwd.h>
#include <THashList.h>
#include <TString.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::photonpair;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

struct TagAndProbe {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};

  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<std::string> fConfigTagPCMCut{"cfgTagPCMCuts", "tag_track", "Tag PCM photon cut. 1 per wagon"};
  Configurable<std::string> fConfigTagPHOSCut{"cfgTagPHOSCuts", "test03", "Tag PHOS photon cuts. 1 per wagon"};
  Configurable<std::string> fConfigTagEMCCut{"cfgTagEMCCuts", "standard", "Tag EMCal photon cuts. 1 per wagon"};
  Configurable<std::string> fConfigProbePCMCuts{"cfgProbePCMCuts", "analysis,qc", "Comma separated list of V0 photon cuts"};
  Configurable<std::string> fConfigProbePHOSCuts{"cfgProbePHOSCuts", "test02,test03", "Comma separated list of PHOS photon cuts"};
  Configurable<std::string> fConfigProbeEMCCuts{"cfgProbeEMCCuts", "standard", "Comma separated list of EMCal photon cuts"};
  Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMPhotonEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  THashList* fMainList = new THashList();

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

    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

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
      THashList* list_pair_subsys_photoncut = o2::aod::pwgem::photon::histogram::AddHistClass(list_pair, photon_cut_name.data());

      for (auto& cut3 : paircuts) {
        std::string pair_cut_name = cut3.GetName();
        o2::aod::pwgem::photon::histogram::AddHistClass(list_pair_subsys_photoncut, pair_cut_name.data());
        THashList* list_pair_subsys_paircut = reinterpret_cast<THashList*>(list_pair_subsys_photoncut->FindObject(pair_cut_name.data()));
        o2::aod::pwgem::photon::histogram::DefineHistograms(list_pair_subsys_paircut, "tag_and_probe", pairname.data());
      } // end of cut3 loop pair cut
    } // end of cut2 loop
  }

  static constexpr std::string_view pairnames[6] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PHOSEMC"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Pair");
    THashList* list_pair = reinterpret_cast<THashList*>(fMainList->FindObject("Pair"));

    // create sub lists first.

    for (auto& pairname : fPairNames) {
      LOGF(info, "Enabled pairs = %s", pairname.data());

      THashList* list_ev_pair = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev, pairname.data()));
      for (const auto& evtype : event_types) {
        THashList* list_ev_type = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev_pair, evtype.data()));
        o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_type, "Event", evtype.data());
      }

      THashList* list_pair_ss = o2::aod::pwgem::photon::histogram::AddHistClass(list_pair, pairname.data());

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
        fProbeEMCCuts.push_back(*emccuts::GetCut(cutname));
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

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emeventId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TTagCut, typename TProbeCuts, typename TPairCuts, typename TLegs>
  void SameEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TTagCut const& tagcut, TProbeCuts const& probecuts, TPairCuts const& paircuts, TLegs const& /*legs*/)
  {
    THashList* list_ev_pair_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject(event_types[0].data()));
    THashList* list_ev_pair_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject(event_types[1].data()));
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

    for (auto& collision : collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if ((pairtype == PairType::kPHOSPHOS || pairtype == PairType::kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }
      if ((pairtype == PairType::kEMCEMC || pairtype == PairType::kPCMEMC) && !collision.alias_bit(triggerAliases::kTVXinEMC)) {
        continue;
      }

      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_pair_before, "", collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_pair_after, "", collision);
      reinterpret_cast<TH1F*>(list_ev_pair_before->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);
      reinterpret_cast<TH1F*>(list_ev_pair_after->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);

      auto photons1_coll = photons1.sliceBy(perCollision1, collision.globalIndex());
      auto photons2_coll = photons2.sliceBy(perCollision2, collision.globalIndex());

      for (auto& g1 : photons1_coll) {

        if constexpr (pairtype == PairType::kPCMPCM) {
          if (!tagcut.template IsSelected<decltype(g1), TLegs>(g1)) {
            continue;
          }
        } else if constexpr (pairtype == PairType::kPHOSPHOS) {
          if (!tagcut.template IsSelected<decltype(g1)>(g1)) {
            continue;
          }
        } else if constexpr (pairtype == PairType::kEMCEMC) {
          if (!tagcut.template IsSelected<decltype(g1)>(g1)) {
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
                if (!probecut.template IsSelected<decltype(g2), TLegs>(g2)) {
                  continue;
                }
              } else if constexpr (pairtype == PairType::kPHOSPHOS) {
                if (!probecut.template IsSelected<decltype(g2)>(g2)) {
                  continue;
                }
              } else if constexpr (pairtype == PairType::kEMCEMC) {
                if (!probecut.template IsSelected<decltype(g2)>(g2)) {
                  continue;
                }
              }

              reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", tagcut.GetName(), probecut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_PassingProbe_Same"))->Fill(v12.M(), v2.Pt());

              if constexpr (pairtype == PairType::kEMCEMC) {
                RotationBackground<aod::SkimEMCClusters>(v12, v1, v2, photons2_coll, g1.globalIndex(), g2.globalIndex(), probecut, paircut);
              }
            } // end of probe cut loop
          } // end of pair cut loop
        } // end of g2 loop
      } // end of g1 loop
    } // end of collision loop
  }

  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 999.f}, "Mixing bins - centrality"};
  using BinningType_M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningType_A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningType_C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType_M colBinning_M{{ConfVtxBins, ConfCentBins}, true};
  BinningType_A colBinning_A{{ConfVtxBins, ConfCentBins}, true};
  BinningType_C colBinning_C{{ConfVtxBins, ConfCentBins}, true};

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TTagCut, typename TProbeCuts, typename TPairCuts, typename TLegs, typename TMixedBinning>
  void MixedEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TTagCut const& tagcut, TProbeCuts const& probecuts, TPairCuts const& paircuts, TLegs const& /*legs*/, TMixedBinning const& colBinning)
  {
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.

      const float centralities1[3] = {collision1.centFT0M(), collision1.centFT0A(), collision1.centFT0C()};
      const float centralities2[3] = {collision2.centFT0M(), collision2.centFT0A(), collision2.centFT0C()};

      if (centralities1[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities1[cfgCentEstimator]) {
        continue;
      }
      if (centralities2[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities2[cfgCentEstimator]) {
        continue;
      }
      if (!fEMEventCut.IsSelected(collision1) || !fEMEventCut.IsSelected(collision2)) {
        continue;
      }

      auto photons_coll1 = photons1.sliceBy(perCollision1, collision1.globalIndex());
      auto photons_coll2 = photons2.sliceBy(perCollision2, collision2.globalIndex());
      // LOGF(info, "collision1: posZ = %f, numContrib = %d , sel8 = %d | collision2: posZ = %f, numContrib = %d , sel8 = %d", collision1.posZ(), collision1.numContrib(), collision1.sel8(), collision2.posZ(), collision2.numContrib(), collision2.sel8());

      for (auto& g1 : photons_coll1) {
        if constexpr (pairtype == PairType::kPCMPCM) {
          if (!tagcut.template IsSelected<decltype(g1), TLegs>(g1)) {
            continue;
          }
        } else if constexpr (pairtype == PairType::kPHOSPHOS) {
          if (!tagcut.template IsSelected<decltype(g1)>(g1)) {
            continue;
          }
        } else if constexpr (pairtype == PairType::kEMCEMC) {
          if (!tagcut.template IsSelected<decltype(g1)>(g1)) {
            continue;
          }
        }
        for (auto& g2 : photons_coll2) {
          // LOGF(info, "Mixed event photon pair: (%d, %d) from events (%d, %d), photon event: (%d, %d)", g1.index(), g2.index(), collision1.index(), collision2.index(), g1.globalIndex(), g2.globalIndex());

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
                if (!probecut.template IsSelected<decltype(g2), TLegs>(g2)) {
                  continue;
                }
              } else if constexpr (pairtype == PairType::kPHOSPHOS) {
                if (!probecut.template IsSelected<decltype(g2)>(g2)) {
                  continue;
                }
              } else if constexpr (pairtype == PairType::kEMCEMC) {
                if (!probecut.template IsSelected<decltype(g2)>(g2)) {
                  continue;
                }
              }

              reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", tagcut.GetName(), probecut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_PassingProbe_Mixed"))->Fill(v12.M(), v2.Pt());

            } // end of probe cut loop
          } // end of pair cut loop
        } // end of g2 loop
      } // end of g1 loop
    } // end of different collision combinations
  }

  /// \brief Calculate background (using rotation background method only for EMCal!)
  template <typename TPhotons>
  void RotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, unsigned int ig1, unsigned int ig2, EMCPhotonCut const& cut, PairCut const& paircut)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photons_coll.size() < 3) {
      return;
    }
    const double rotationAngle = o2::constants::math::PIHalf; // rotaion angle 90°
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
      if (!cut.template IsSelected<decltype(photon)>(photon)) {
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

  Partition<MyCollisions> grouped_collisions = cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax; // this goes to same event.
  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > static_cast<uint16_t>(0) && o2::aod::evsel::sel8 == true;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using MyFilteredCollisions = soa::Filtered<MyCollisions>; // this goes to mixed event.

  void processPCMPCM(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::V0Legs const& legs)
  {
    SameEventPairing<PairType::kPCMPCM>(grouped_collisions, v0photons, v0photons, perCollision, perCollision, fTagPCMCut, fProbePCMCuts, fPairCuts, legs);
    if (cfgCentEstimator == 0) {
      MixedEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision, perCollision, fTagPCMCut, fProbePCMCuts, fPairCuts, legs, colBinning_M);
    } else if (cfgCentEstimator == 1) {
      MixedEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision, perCollision, fTagPCMCut, fProbePCMCuts, fPairCuts, legs, colBinning_A);
    } else if (cfgCentEstimator == 2) {
      MixedEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision, perCollision, fTagPCMCut, fProbePCMCuts, fPairCuts, legs, colBinning_C);
    }
  }

  void processPHOSPHOS(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, aod::PHOSClusters const& phosclusters)
  {
    SameEventPairing<PairType::kPHOSPHOS>(grouped_collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos, fTagPHOSCut, fProbePHOSCuts, fPairCuts, nullptr);
    MixedEventPairing<PairType::kPHOSPHOS>(filtered_collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos, fTagPHOSCut, fProbePHOSCuts, fPairCuts, nullptr, colBinning_C);
  }

  void processEMCEMC(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, aod::SkimEMCClusters const& emcclusters)
  {
    SameEventPairing<PairType::kEMCEMC>(grouped_collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc, fTagEMCCut, fProbeEMCCuts, fPairCuts, nullptr);
    MixedEventPairing<PairType::kEMCEMC>(filtered_collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc, fTagEMCCut, fProbeEMCCuts, fPairCuts, nullptr, colBinning_C);
  }

  void processDummy(MyCollisions const&) {}

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
