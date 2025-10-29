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

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/DataModel/CaloClusters.h"
#include "EventFiltering/filterTables.h"

#include "DataFormatsPHOS/TriggerRecord.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollision = MyCollisions::iterator;

using MyPrimaryElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronsPrefilterBit>;
using MyPrimaryElectron = MyPrimaryElectrons::iterator;

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
    kPCM_EE = 5,
    kNtrg
  };

  Produces<aod::PhotonFilters> tags;
  Produces<aod::EMSwtInfosPCM> swtinfo_pcm;
  Produces<aod::EMSwtInfosPair> swtinfo_pair;

  Configurable<float> ePhot{"ePhot", 2.2, "Minimal photon energy (GeV)"};
  Configurable<float> eEl{"eEl", 1., "Minimal electron energy (GeV)"};
  Configurable<float> ePair{"ePair", 0.35, "Minimal photon pair mass (GeV)"};
  Configurable<int> nNbar{"nNbar", 2, "Minimal number of nbar clusters"};

  // for PCM
  Configurable<float> minpt_v0{"minpt_v0", 0.1, "min pt for v0"};
  Configurable<float> maxeta_v0{"maxeta_v0", 0.9, "eta acceptance for v0"};
  Configurable<float> min_pt_pcm_photon{"min_pt_pcm_photon", 0.f, "min. pT for PCM photon"};
  Configurable<float> minTPCNsigmaEl_v0{"minTPCNsigmaEl_v0", -3.5, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl_v0{"maxTPCNsigmaEl_v0", +3.5, "max. TPC n sigma for electron inclusion"};
  Configurable<float> max_dcatopv_xy_v0{"max_dcatopv_xy_v0", +1e+10, "max. DCAxy to PV for V0"};
  Configurable<float> max_dcatopv_z_v0{"max_dcatopv_z_v0", +1e+10, "max. DCAz to PV for V0"};

  // for prompt dielectron
  Configurable<float> minpt{"minpt", 0.2, "min pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> dca_3d_sigma_max{"dca_3d_sigma_max", 1.0f, "max DCA 3D in sigma"}; // for single track
  Configurable<int> mincrossedrows{"mincrossedrows", 80, "min crossed rows"};
  Configurable<float> minTPCNsigmaEl_primary{"minTPCNsigmaEl_primary", -3.0, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl_primary{"maxTPCNsigmaEl_primary", +4.0, "max. TPC n sigma for electron inclusion"};
  Configurable<float> minTOFNsigmaEl_primary{"minTOFNsigmaEl_primary", -4.0, "min. TOF n sigma for electron inclusion"}; // require TOF
  Configurable<float> maxTOFNsigmaEl_primary{"maxTOFNsigmaEl_primary", +4.0, "max. TOF n sigma for electron inclusion"}; // require TOF
  Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", -2.0, "min. TPC n sigma for pion exclusion"};                     // set to -2 for lowB, -999 for nominalB
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 2.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> min_pin_tof{"min_pin_tof", 0.4, "tof is required above this threshold in pin"};
  Configurable<float> slope{"slope", 0.0185, "slope for m vs. phiv"};
  Configurable<float> intercept{"intercept", -0.0280, "intercept for m vs. phiv"};
  Configurable<float> min_mee{"min_mee", 0.0, "min. mee"};
  Configurable<float> max_mee{"max_mee", 0.5, "max. mee"};
  Configurable<float> min_meeg{"min_meeg", 0.3, "min. meeg"};
  Configurable<float> max_meeg{"max_meeg", 0.8, "max. meeg"};
  Configurable<bool> applyPF{"applyPF", false, "apply pre-filter for primary electron"}; // i.e. reject electron from photon conversion with phiv

  HistogramRegistry mHistManager{"events", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    auto scalers{std::get<std::shared_ptr<TH1>>(mHistManager.add("hEventCounter", "Number of filtered events", HistType::kTH1F, {{20, 0.5, 20.5}}))};
    scalers->GetXaxis()->SetBinLabel(1, "all events");
    scalers->GetXaxis()->SetBinLabel(2, "sel8");
    scalers->GetXaxis()->SetBinLabel(3, "|Z_{vtx}| < 10 cm");
    scalers->GetXaxis()->SetBinLabel(4, "No TFB");
    scalers->GetXaxis()->SetBinLabel(5, "sel8 && |Z_{vtx}| < 10 cm");
    scalers->GetXaxis()->SetBinLabel(6, "sel8 && |Z_{vtx}| < 10 cm && No TFB");
    scalers->GetXaxis()->SetBinLabel(7, "PHOS photon");
    scalers->GetXaxis()->SetBinLabel(8, "PHOS electron");
    scalers->GetXaxis()->SetBinLabel(9, "PHOS pair");
    scalers->GetXaxis()->SetBinLabel(10, "PHOS nbar");
    scalers->GetXaxis()->SetBinLabel(11, "PHOS photon & electron");
    scalers->GetXaxis()->SetBinLabel(12, "PHOS photon & pair");
    scalers->GetXaxis()->SetBinLabel(13, "events with PHOS");
    scalers->GetXaxis()->SetBinLabel(14, "PCM high p_{T} photon");
    scalers->GetXaxis()->SetBinLabel(15, "PCM #gamma and dielectron");
  }

  template <typename TTrack>
  bool isSelectedSecondary(TTrack const& track)
  {
    if (track.hasTPC() && (track.tpcNSigmaEl() < minTPCNsigmaEl_v0 || maxTPCNsigmaEl_v0 < track.tpcNSigmaEl())) {
      return false;
    }
    return true;
  }
  template <typename TTrack>
  bool isSelectedPrimary(TTrack const& track)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.pt() < minpt || std::fabs(track.eta()) > maxeta) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }

    if (track.tpcNSigmaEl() < minTPCNsigmaEl_primary || maxTPCNsigmaEl_primary < track.tpcNSigmaEl()) {
      return false;
    }

    if ((track.tofNSigmaEl() < minTOFNsigmaEl_primary || maxTOFNsigmaEl_primary < track.tofNSigmaEl()) && track.tpcInnerParam() > min_pin_tof) {
      return false;
    }

    if (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi) {
      return false;
    }

    if (applyPF && track.pfb() > 0) {
      return false;
    }

    float dca_3d = 999.f;
    float det = track.cYY() * track.cZZ() - track.cZY() * track.cZY();
    if (det < 0) {
      dca_3d = 999.f;
    } else {
      float chi2 = (track.dcaXY() * track.dcaXY() * track.cZZ() + track.dcaZ() * track.dcaZ() * track.cYY() - 2. * track.dcaXY() * track.dcaZ() * track.cZY()) / det;
      dca_3d = std::sqrt(std::fabs(chi2) / 2.);
    }
    if (dca_3d > dca_3d_sigma_max) {
      return false;
    }
    return true;
  }

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::DalitzEEs> perCollision_ee = aod::dalitzee::collisionId;
  Preslice<aod::CaloClusters> perCollision_phos = aod::calocluster::collisionId;
  // Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <uint8_t system, typename TCollisions, typename TPhotons1, typename TPhotons2, typename TPhotons3, typename TV0Legs, typename TDielectrons, typename TEMPrimaryElectrons>
  void runFilter(TCollisions const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPhotons3 const& /*photons3*/, TV0Legs const&, TDielectrons const& dielectrons, TEMPrimaryElectrons const& /*emprimaryelectrons*/)
  {
    for (const auto& collision : collisions) {
      mHistManager.fill(HIST("hEventCounter"), 1.);
      bool keepEvent[kNtrg]{false};

      if (collision.sel8()) {
        mHistManager.fill(HIST("hEventCounter"), 2.);
      }
      if (std::fabs(collision.posZ()) < 10.f) {
        mHistManager.fill(HIST("hEventCounter"), 3.);
      }
      if (collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        mHistManager.fill(HIST("hEventCounter"), 4.);
      }
      if (collision.sel8() && std::fabs(collision.posZ()) < 10.f) {
        mHistManager.fill(HIST("hEventCounter"), 5.);
      }
      if (collision.sel8() && std::fabs(collision.posZ()) < 10.f && collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        mHistManager.fill(HIST("hEventCounter"), 6.);
      }

      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        tags(false, false, false, false);
        continue;
      }

      if constexpr (static_cast<bool>(system & EM_Filter_PhotonType::kPCM)) {
        auto photons1_per_coll = photons1.sliceBy(perCollision_pcm, collision.globalIndex());
        auto dielectrons_per_coll = dielectrons.sliceBy(perCollision_ee, collision.globalIndex());

        for (const auto& v0photon : photons1_per_coll) {
          auto pos_sv = v0photon.template posTrack_as<TV0Legs>();
          auto ele_sv = v0photon.template negTrack_as<TV0Legs>();
          if (!isSelectedSecondary(pos_sv) || !isSelectedSecondary(ele_sv)) {
            continue;
          }
          if (v0photon.pt() > min_pt_pcm_photon) {
            keepEvent[kPCM_HighPtPhoton] = true;
            swtinfo_pcm(collision.globalIndex(), v0photon.globalIndex());
          }
        } // end of single v0 photon loop

        for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_per_coll, dielectrons_per_coll))) {
          auto pos_sv = g1.template posTrack_as<TV0Legs>();
          auto ele_sv = g1.template negTrack_as<TV0Legs>();
          if (!isSelectedSecondary(pos_sv) || !isSelectedSecondary(ele_sv)) {
            continue;
          }

          auto pos_pv = g2.template posTrack_as<TEMPrimaryElectrons>();
          auto ele_pv = g2.template negTrack_as<TEMPrimaryElectrons>();
          if (!isSelectedPrimary(pos_pv) || !isSelectedPrimary(ele_pv)) {
            continue;
          }
          if (g2.mass() < min_mee || max_mee < g2.mass()) {
            continue;
          }
          if (g2.mass() < slope * g2.phiv() + intercept) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (min_meeg < v12.M() && v12.M() < max_meeg) {
            keepEvent[kPCM_EE] = true;
            swtinfo_pair(collision.globalIndex(), g1.globalIndex(), g2.globalIndex());
          }
        } // end of photon + dielectron pair loop

        if (keepEvent[kPCM_HighPtPhoton] == true) {
          mHistManager.fill(HIST("hEventCounter"), 14);
        }
        if (keepEvent[kPCM_EE] == true) {
          mHistManager.fill(HIST("hEventCounter"), 15);
        }
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
              double m = std::pow(clu.e() + clu2.e(), 2) - std::pow(clu.px() + clu2.px(), 2) -
                         std::pow(clu.py() + clu2.py(), 2) - std::pow(clu.pz() + clu2.pz(), 2);
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
          mHistManager.fill(HIST("hEventCounter"), 13.);
        }
        // Can not fill with variable, have to fill manually
        if (keepEvent[kPHOS_Photon]) {
          mHistManager.fill(HIST("hEventCounter"), 7.);
          if (keepEvent[kPHOS_El]) {
            mHistManager.fill(HIST("hEventCounter"), 11.);
          }
          if (keepEvent[kPHOS_Pair]) {
            mHistManager.fill(HIST("hEventCounter"), 12.);
          }
        }
        if (keepEvent[kPHOS_El]) {
          mHistManager.fill(HIST("hEventCounter"), 8.);
        }
        if (keepEvent[kPHOS_Pair]) {
          mHistManager.fill(HIST("hEventCounter"), 9.);
        }
        if (keepEvent[kPHOS_Nbar]) {
          mHistManager.fill(HIST("hEventCounter"), 10.);
        }
      }

      // // EMC decision
      // if constexpr (static_cast<bool>(system & EM_Filter_PhotonType::kEMC)) {
      //   // so far, do nothing.
      // }

      tags(keepEvent[kPHOS_Photon], keepEvent[kPHOS_Nbar], keepEvent[kPCM_HighPtPhoton], keepEvent[kPCM_EE]);

    } // end of collision loop
  }

  Filter DalitzEEFilter = o2::aod::dalitzee::sign == 0; // analyze only uls
  using filteredDalitzEEs = Filtered<aod::DalitzEEs>;

  void process_PCM(MyCollisions const& collisions, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, filteredDalitzEEs const& dielectrons, MyPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t system = EM_Filter_PhotonType::kPCM;
    runFilter<system>(collisions, v0photons, nullptr, nullptr, v0legs, dielectrons, emprimaryelectrons);
  }

  Filter phosCluFilter = (o2::aod::calocluster::e > 0.3f);
  using CluCandidates = Filtered<o2::aod::CaloClusters>;
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

  void process_PCM_PHOS(MyCollisions const& collisions, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, filteredDalitzEEs const& dielectrons, MyPrimaryElectrons const& emprimaryelectrons, CluCandidates const& clusters)
  {
    const uint8_t system = EM_Filter_PhotonType::kPCM | EM_Filter_PhotonType::kPHOS;
    runFilter<system>(collisions, v0photons, clusters, nullptr, v0legs, dielectrons, emprimaryelectrons);
  }

  void processDummy(MyCollisions const& collisions)
  {
    for (int i = 0; i < collisions.size(); i++) {
      tags(false, false, false, false);
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
