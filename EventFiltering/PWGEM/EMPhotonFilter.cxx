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
// O2 includes

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/CaloClusters.h"
#include "DataFormatsPHOS/TriggerRecord.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "../filterTables.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0Recalculation>;
using MyV0Photon = MyV0Photons::iterator;

struct EMPhotonFilter {

  enum EM_Filter_PhotonType {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
  };

  enum trigs {
    kPhot = 0,
    kEl = 1,
    kPair = 2,
    kNbar = 3,
    kPCM_Wwire = 4,
    kNtrg
  };

  Produces<aod::PhotonFilters> tags;

  Configurable<float> ePhot{"ePhot", 2.2, "Minimal photon energy (GeV)"};
  Configurable<float> eEl{"eEl", 1., "Minimal electron energy (GeV)"};
  Configurable<float> ePair{"ePair", 0.35, "Minimal photon pair mass (GeV)"};
  Configurable<int> nNbar{"nNbar", 2, "Minimal number of nbar clusters"};

  HistogramRegistry mHistManager{"events", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    auto scalers{std::get<std::shared_ptr<TH1>>(mHistManager.add("fProcessedEvents", "Number of filtered events", HistType::kTH1F, {{20, 0.5, 20.5}}))};
    scalers->GetXaxis()->SetBinLabel(1, "Processed events");
    scalers->GetXaxis()->SetBinLabel(2, "PHOS photon");
    scalers->GetXaxis()->SetBinLabel(3, "PHOS electron");
    scalers->GetXaxis()->SetBinLabel(4, "PHOS pair");
    scalers->GetXaxis()->SetBinLabel(5, "PHOS nbar");
    scalers->GetXaxis()->SetBinLabel(6, "PHOS photon & electron");
    scalers->GetXaxis()->SetBinLabel(7, "PHOS photon & pair");
    scalers->GetXaxis()->SetBinLabel(8, "events with PHOS");
    scalers->GetXaxis()->SetBinLabel(9, "PCM W wire"); // RF is expected to be 30k.
  }

  bool checkAP(float alpha, float qt)
  {
    const float alpha_max = 0.95;
    const float qt_max = 0.05;
    if (pow(alpha / alpha_max, 2) + pow(qt / qt_max, 2) < 1.0) {
      return true;
    } else {
      return false;
    }
  }

  V0PhotonCut cut_pcm_wwire = *pcmcuts::GetCut("wwire");

  Preslice<MyV0Photons> perCollision_pcm = aod::v0photon::collisionId;
  Preslice<aod::CaloClusters> perCollision_phos = aod::calocluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <uint8_t system, typename TCollisions, typename TPhotons1, typename TPhotons2, typename TPhotons3, typename TV0Legs>
  void runFilter(TCollisions const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPhotons3 const& photons3, TV0Legs const&)
  {

    for (auto& collision : collisions) {
      mHistManager.fill(HIST("fProcessedEvents"), 1.);
      bool keepEvent[kNtrg]{false};
      if constexpr (static_cast<bool>(system & EM_Filter_PhotonType::kPCM)) {
        auto photons1_coll = photons1.sliceBy(perCollision_pcm, collision.globalIndex());
        for (auto& v0 : photons1_coll) {
          if (!checkAP(v0.alpha(), v0.qtarm())) {
            continue;
          }
          if (cut_pcm_wwire.template IsSelected<TV0Legs>(v0)) {
            keepEvent[kPCM_Wwire] = true;
            mHistManager.fill(HIST("fProcessedEvents"), 9.);
            break; // if at least 1 photon conversion on W wire is found, break from v0 photon loop.
          }

        } // end of v0 photon loop
      } else if constexpr (static_cast<bool>(system & EM_Filter_PhotonType::kPHOS)) {
        int nPHOSclu = 0;
        int nPHOSnbar = 0;
        auto photons2_coll = photons2.sliceBy(perCollision_phos, collision.globalIndex());
        for (const auto& clu : photons2_coll) {
          nPHOSclu++;

          // Scan current cluster
          //  photons
          keepEvent[kPhot] |= (clu.e() > ePhot);
          // charged clusters above threshold
          keepEvent[kEl] |= (clu.trackdist() < 2. && clu.e() > eEl); // 2: Distance to CPV cluster in sigmas
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
            for (; !keepEvent[kPair] && clu2 != photons2.end(); clu2++) {
              // cluster selections
              if (clu2.trackdist() < 1.) { // select neutral clusters. Disp, Ncell cuts?
                continue;
              }
              double m = pow(clu.e() + clu2.e(), 2) - pow(clu.px() + clu2.px(), 2) -
                         pow(clu.py() + clu2.py(), 2) - pow(clu.pz() + clu2.pz(), 2);
              if (m > ePair * ePair) {
                keepEvent[kPair] |= true;
                break;
              }
            }
          }
        } // end of cluster loop
        keepEvent[kNbar] = (nPHOSnbar >= nNbar);

        // Collision processed, fill scalers here
        if (nPHOSclu) {
          mHistManager.fill(HIST("fProcessedEvents"), 8.);
        }
        // Can not fill with variable, have to fill manually
        if (keepEvent[kPhot]) {
          mHistManager.fill(HIST("fProcessedEvents"), 2.);
          if (keepEvent[kEl]) {
            mHistManager.fill(HIST("fProcessedEvents"), 6.);
          }
          if (keepEvent[kPair]) {
            mHistManager.fill(HIST("fProcessedEvents"), 7.);
          }
        }
        if (keepEvent[kEl]) {
          mHistManager.fill(HIST("fProcessedEvents"), 3.);
        }
        if (keepEvent[kPair]) {
          mHistManager.fill(HIST("fProcessedEvents"), 4.);
        }
        if (keepEvent[kNbar]) {
          mHistManager.fill(HIST("fProcessedEvents"), 5.);
        }
      } else if constexpr (static_cast<bool>(system & EM_Filter_PhotonType::kEMC)) {
        // so far, do nothing.
      }
      tags(keepEvent[kPhot], keepEvent[kEl], keepEvent[kPair], keepEvent[kNbar], keepEvent[kPCM_Wwire]);
    } // end of collision loop
  }

  void processPCM(aod::Collisions const& collisions, MyV0Photons const& v0photons, aod::V0Legs const& v0legs)
  {
    const uint8_t system = EM_Filter_PhotonType::kPCM;
    runFilter<system>(collisions, v0photons, nullptr, nullptr, v0legs);
  }

  Filter phosCluFilter = (o2::aod::calocluster::e > 0.3f);
  using CluCandidates = o2::soa::Filtered<o2::aod::CaloClusters>;
  void processPHOS(aod::Collisions const& collisions, CluCandidates const& clusters)
  {
    const uint8_t system = EM_Filter_PhotonType::kPHOS;
    runFilter<system>(collisions, nullptr, clusters, nullptr, nullptr);
  }

  void processEMC(aod::Collisions const& collisions, aod::SkimEMCClusters const& clusters)
  {
    const uint8_t system = EM_Filter_PhotonType::kEMC;
    runFilter<system>(collisions, nullptr, nullptr, clusters, nullptr);
  }

  PROCESS_SWITCH(EMPhotonFilter, processPCM, "Process PCM software trigger decision", true);
  PROCESS_SWITCH(EMPhotonFilter, processPHOS, "Process PHOS software trigger decision", true);
  PROCESS_SWITCH(EMPhotonFilter, processEMC, "Process EMC software trigger decision", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<EMPhotonFilter>(cfg, TaskName{"em-photon-filter"})};
}
