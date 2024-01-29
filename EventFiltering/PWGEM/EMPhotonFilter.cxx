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

// \brief software trigger for EM photon
// \author daiki.sekihata@cern.ch

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/CaloClusters.h"
#include "DataFormatsPHOS/TriggerRecord.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "EventFiltering/filterTables.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollision = MyCollisions::iterator;

// using MyPrimaryElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsPrefilterBit>;
// using MyPrimaryElectron = MyPrimaryElectrons::iterator;

struct EMPhotonFilter {

  enum EM_Filter_PhotonType {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
  };

  enum trigs {
    kPHOS_Photon = 0,
    kPHOS_El = 1,
    kPHOS_Pair = 2,
    kPHOS_Nbar = 3,
    kPCM_HighPtPhoton = 4,
    kPCM_MatCalib = 5,
    kPCM_EtaDalitz = 6,
    kPCM_EtaGG = 7,
    kPCM_EE = 8,
    kNtrg
  };

  Produces<aod::PhotonFilters> tags;

  Configurable<float> ePhot{"ePhot", 2.2, "Minimal photon energy (GeV)"};
  Configurable<float> eEl{"eEl", 1., "Minimal electron energy (GeV)"};
  Configurable<float> ePair{"ePair", 0.35, "Minimal photon pair mass (GeV)"};
  Configurable<int> nNbar{"nNbar", 2, "Minimal number of nbar clusters"};

  // for PCM
  Configurable<float> min_pt_tagging{"min_pt_tagging", 0.f, "min. pT for tagging"};
  Configurable<float> max_mee_pi0_dalitz{"max_mee_pi0_dalitz", 0.12, "max. mee for pi0 dalitz decay"};
  Configurable<float> min_meeg_pi0{"min_meeg_pi0", 0.04, "min. meeg for pi0"};
  Configurable<float> max_meeg_pi0{"max_meeg_pi0", 0.24, "max. meeg for pi0"};
  Configurable<float> max_mee_eta_dalitz{"max_mee_eta_dalitz", 0.5, "max. mee for eta dalitz decay"};
  Configurable<float> min_meeg_eta{"min_meeg_eta", 0.35, "min. meeg for eta"};
  Configurable<float> max_meeg_eta{"max_meeg_eta", 0.75, "max. meeg for eta"};
  Configurable<float> slope{"slope", 0.0185, "slope for m vs. phiv"};
  Configurable<float> intercept{"intercept", -0.0280, "intercept for m vs. phiv"};
  Configurable<float> min_pt_pcm_photon{"min_pt_pcm_photon", 4.f, "min. pT for PCM photon"};

  HistogramRegistry mHistManager{"events", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    auto scalers{std::get<std::shared_ptr<TH1>>(mHistManager.add("hEventCounter", "Number of filtered events", HistType::kTH1F, {{20, 0.5, 20.5}}))};
    scalers->GetXaxis()->SetBinLabel(1, "all events");
    scalers->GetXaxis()->SetBinLabel(2, "sel8");
    scalers->GetXaxis()->SetBinLabel(3, "|Z_{vtx}| < 10 cm");
    scalers->GetXaxis()->SetBinLabel(4, "sel8 && |Z_{vtx}| < 10 cm");
    scalers->GetXaxis()->SetBinLabel(5, "PHOS photon");
    scalers->GetXaxis()->SetBinLabel(6, "PHOS electron");
    scalers->GetXaxis()->SetBinLabel(7, "PHOS pair");
    scalers->GetXaxis()->SetBinLabel(8, "PHOS nbar");
    scalers->GetXaxis()->SetBinLabel(9, "PHOS photon & electron");
    scalers->GetXaxis()->SetBinLabel(10, "PHOS photon & pair");
    scalers->GetXaxis()->SetBinLabel(11, "events with PHOS");
    scalers->GetXaxis()->SetBinLabel(12, "PCM high p_{T} photon");
    scalers->GetXaxis()->SetBinLabel(13, "PCM Material budget calibration");
    scalers->GetXaxis()->SetBinLabel(14, "PCM #eta #rightarrow ee#gamma");
    scalers->GetXaxis()->SetBinLabel(15, "PCM #eta #rightarrow #gamma#gamma");
    scalers->GetXaxis()->SetBinLabel(16, "PCM DalitzEE #gamma-#gamma^{*} BEC");
  }

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  // Preslice<aod::DalitzEEs> perCollision_ee = aod::dalitzee::collisionId;
  Preslice<aod::CaloClusters> perCollision_phos = aod::calocluster::collisionId;
  // Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <uint8_t system, typename TCollisions, typename TPhotons1, typename TPhotons2, typename TPhotons3, typename TV0Legs, typename TDielectrons, typename TEMPrimaryElectrons>
  void runFilter(TCollisions const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPhotons3 const& photons3, TV0Legs const&, TDielectrons const& dielectrons, TEMPrimaryElectrons const& emprimaryelectrons)
  {
    for (auto& collision : collisions) {
      mHistManager.fill(HIST("hEventCounter"), 1.);
      bool keepEvent[kNtrg]{false};

      if (collision.sel8()) {
        mHistManager.fill(HIST("hEventCounter"), 2.);
      }
      if (abs(collision.posZ()) < 10.f) {
        mHistManager.fill(HIST("hEventCounter"), 3.);
      }
      if (collision.sel8() && abs(collision.posZ()) < 10.f) {
        mHistManager.fill(HIST("hEventCounter"), 4.);
      }

      if constexpr (static_cast<bool>(system & EM_Filter_PhotonType::kPCM)) {
        auto photons1_per_coll = photons1.sliceBy(perCollision_pcm, collision.globalIndex());
        // auto dielectrons_per_coll = dielectrons.sliceBy(perCollision_ee, collision.globalIndex());

        for (auto& v0photon : photons1_per_coll) {
          if (v0photon.pt() > min_pt_pcm_photon) {
            keepEvent[kPCM_HighPtPhoton] = true;
            mHistManager.fill(HIST("hEventCounter"), 12);
            break;
          }
        } // end of single v0 photon loop

        // for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_per_coll, dielectrons_per_coll))) {
        //   if (g2.pt() < min_pt_tagging) { // this is only to increase rejection factor
        //     continue;
        //   }
        //   if (g2.mass() > max_mee_pi0_dalitz) { // select only pi0 candidates
        //     continue;
        //   }
        //   if (g2.mass() < slope * g2.phiv() + intercept) {
        //     continue;
        //   }
        //   ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        //   ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
        //   ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        //   if (min_meeg_pi0 < v12.M() && v12.M() < max_meeg_pi0) {
        //     keepEvent[kPCM_MatCalib] = true;
        //     mHistManager.fill(HIST("hEventCounter"), 13);
        //     break;
        //   }

        // } // end of dielectron-photon pair loop

        // for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_per_coll, dielectrons_per_coll))) {
        //   if (g2.mass() > max_mee_eta_dalitz) { // select only eta candidates
        //     continue;
        //   }
        //   if (g2.mass() < slope * g2.phiv() + intercept) {
        //     continue;
        //   }

        //   ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        //   ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
        //   ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        //   if (min_meeg_eta < v12.M() && v12.M() < max_meeg_eta) { // eta -> eeg
        //     keepEvent[kPCM_EtaDalitz] = true;
        //     mHistManager.fill(HIST("hEventCounter"), 14);
        //     break;
        //   }
        // } // end of dielectron-photon pair loop

        // for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_per_coll, photons1_per_coll))) {
        //   ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        //   ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        //   ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        //   if (min_meeg_eta < v12.M() && v12.M() < max_meeg_eta) { // eta -> gg
        //     keepEvent[kPCM_EtaGG] = true;
        //     mHistManager.fill(HIST("hEventCounter"), 15);
        //     break;
        //   }
        // } // end of photon-photon pair loop

        // for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_per_coll, dielectrons_per_coll))) {
        //   if (g2.mass() < slope * g2.phiv() + intercept) {
        //     continue;
        //   }
        //   keepEvent[kPCM_EE] = true;
        //   mHistManager.fill(HIST("hEventCounter"), 16);
        //   break;
        // } // end of dielectron-photon pair loop

      } // end of PCM decision

      if constexpr (static_cast<bool>(system & EM_Filter_PhotonType::kPHOS)) {
        int nPHOSclu = 0;
        int nPHOSnbar = 0;
        auto photons2_coll = photons2.sliceBy(perCollision_phos, collision.globalIndex());
        for (const auto& clu : photons2_coll) {
          nPHOSclu++;

          // Scan current cluster
          //  photons
          keepEvent[kPHOS_Photon] |= (clu.e() > ePhot);
          // charged clusters above threshold
          keepEvent[kPHOS_El] |= (clu.trackdist() < 2. && clu.e() > eEl); // 2: Distance to CPV cluster in sigmas
          // antineutrons
          if ((clu.ncell() > 2 && clu.m02() > 0.2 && clu.e() > 0.7 && clu.trackdist() > 2.) &&
              ((clu.e() < 2. && clu.m02() > 4.5 - clu.m20()) ||
               (clu.e() > 2. && clu.m02() > 4. - clu.m20()))) {
            nPHOSnbar++;
          }

          // inv mass
          if (clu.trackdist() < 1.) {
            auto clu2 = clu;
            ++clu2;
            for (; !keepEvent[kPHOS_Pair] && clu2 != photons2.end(); clu2++) {
              // cluster selections
              if (clu2.trackdist() < 1.) { // select neutral clusters. Disp, Ncell cuts?
                continue;
              }
              double m = pow(clu.e() + clu2.e(), 2) - pow(clu.px() + clu2.px(), 2) -
                         pow(clu.py() + clu2.py(), 2) - pow(clu.pz() + clu2.pz(), 2);
              if (m > ePair * ePair) {
                keepEvent[kPHOS_Pair] |= true;
                break;
              }
            }
          }
        } // end of cluster loop
        keepEvent[kPHOS_Nbar] = (nPHOSnbar >= nNbar);

        // Collision processed, fill scalers here
        if (nPHOSclu) {
          mHistManager.fill(HIST("hEventCounter"), 11.);
        }
        // Can not fill with variable, have to fill manually
        if (keepEvent[kPHOS_Photon]) {
          mHistManager.fill(HIST("hEventCounter"), 5.);
          if (keepEvent[kPHOS_El]) {
            mHistManager.fill(HIST("hEventCounter"), 9.);
          }
          if (keepEvent[kPHOS_Pair]) {
            mHistManager.fill(HIST("hEventCounter"), 10.);
          }
        }
        if (keepEvent[kPHOS_El]) {
          mHistManager.fill(HIST("hEventCounter"), 6.);
        }
        if (keepEvent[kPHOS_Pair]) {
          mHistManager.fill(HIST("hEventCounter"), 7.);
        }
        if (keepEvent[kPHOS_Nbar]) {
          mHistManager.fill(HIST("hEventCounter"), 8.);
        }
      }

      // // EMC decision
      // if constexpr (static_cast<bool>(system & EM_Filter_PhotonType::kEMC)) {
      //   // so far, do nothing.
      // }
      tags(keepEvent[kPHOS_Photon], keepEvent[kPHOS_Nbar], keepEvent[kPCM_HighPtPhoton]);
    } // end of collision loop
  }

  // void process_PCM(MyCollisions const& collisions, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::DalitzEEs const& dielectrons, MyPrimaryElectrons const& emprimaryelectrons)
  void process_PCM(MyCollisions const& collisions, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs)
  {
    const uint8_t system = EM_Filter_PhotonType::kPCM;
    runFilter<system>(collisions, v0photons, nullptr, nullptr, v0legs, nullptr, nullptr);
  }

  Filter phosCluFilter = (o2::aod::calocluster::e > 0.3f);
  using CluCandidates = o2::soa::Filtered<o2::aod::CaloClusters>;
  void process_PHOS(MyCollisions const& collisions, CluCandidates const& clusters)
  {
    const uint8_t system = EM_Filter_PhotonType::kPHOS;
    runFilter<system>(collisions, nullptr, clusters, nullptr, nullptr, nullptr, nullptr);
  }

  void process_EMC(MyCollisions const& collisions, aod::SkimEMCClusters const& clusters)
  {
    const uint8_t system = EM_Filter_PhotonType::kEMC;
    runFilter<system>(collisions, nullptr, nullptr, clusters, nullptr, nullptr, nullptr);
  }

  // void process_PCM_PHOS(MyCollisions const& collisions, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::DalitzEEs const& dielectrons, MyPrimaryElectrons const& emprimaryelectrons, CluCandidates const& clusters)
  void process_PCM_PHOS(MyCollisions const& collisions, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, CluCandidates const& clusters)
  {
    const uint8_t system = EM_Filter_PhotonType::kPCM | EM_Filter_PhotonType::kPHOS;
    // runFilter<system>(collisions, v0photons, clusters, nullptr, v0legs, dielectrons, emprimaryelectrons);
    runFilter<system>(collisions, v0photons, clusters, nullptr, v0legs, nullptr, nullptr);
  }

  void processDummy(MyCollisions const& collisions)
  {
    for (int i = 0; i < collisions.size(); i++) {
      tags(false, false, false);
    }
  }

  PROCESS_SWITCH(EMPhotonFilter, process_PCM, "Process PCM software trigger decision", false);
  PROCESS_SWITCH(EMPhotonFilter, process_PHOS, "Process PHOS software trigger decision", false);
  PROCESS_SWITCH(EMPhotonFilter, process_EMC, "Process EMC software trigger decision", false);
  PROCESS_SWITCH(EMPhotonFilter, process_PCM_PHOS, "Process PCM and PHOS software trigger decision", false);
  PROCESS_SWITCH(EMPhotonFilter, processDummy, "Process dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<EMPhotonFilter>(cfg, TaskName{"em-photon-filter"})};
}
