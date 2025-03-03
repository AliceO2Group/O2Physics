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
/// \brief Nuclei spectra analysis task
/// \author Bong-Hwi Lim
/// \since 2 December 2024

#include <iterator>
#include <algorithm>
#include "Framework/StaticFor.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::track;

struct NucleiSpectra {
  float collisionCentrality;
  constexpr static std::array<PID::ID, 4> nucleiTypes = {PID::Deuteron, PID::Triton, PID::Helium3, PID::Alpha};
  // inline static std::vector<PID::ID> nucleiTypes = {PID::Deuteron, PID::Triton, PID::Helium3, PID::Alpha};
  inline static std::vector<std::string> nucleiNames;
  inline static std::vector<std::string> hnsigmatof;
  inline static std::vector<std::string> hnsigmatpc;
  inline static std::vector<std::string> hnsigmatpctof;

  // Configuration parameters
  Configurable<float> cutEta{"cutEta", 0.8f, "Pseudorapidity cut"};
  Configurable<float> cutPtMin{"cutPtMin", 0.5f, "Minimum pT"};
  Configurable<float> cutPtMax{"cutPtMax", 10.0f, "Maximum pT"};
  Configurable<float> cutDCAz{"cutDCAz", 1.0f, "DCAz cut"};
  Configurable<float> cutDCAxy{"cutDCAxy", 0.1f, "DCAxy cut"};
  Configurable<float> cutNsigmaTPC{"cutNsigmaTPC", 3.0f, "TPC n-sigma PID cut"};
  Configurable<float> cutNsigmaTOF{"cutNsigmaTOF", 3.0f, "TOF n-sigma PID cut"};
  Configurable<bool> useTOF{"useTOF", true, "Use TOF for PID"};
  Configurable<bool> useTrackingPID{"useTrackingPID", false, "Use tracking PID"};

  // Collision candidate
  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels,
                                        aod::TPCMults, aod::PVMults, aod::MultZeqs,
                                        aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;

  // Track candidates
  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksDCA,
                                    aod::pidTPCDe, aod::pidTOFDe,
                                    aod::pidTPCTr, aod::pidTOFTr,
                                    aod::pidTPCHe, aod::pidTOFHe,
                                    aod::pidTPCAl, aod::pidTOFAl>;

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    collisionCentrality = -999.0f;
    for (const auto& pid : nucleiTypes) {
      nucleiNames.push_back(PID::getName(pid));
    }
    for (const auto& name : nucleiNames) {
      for (const auto& charge : {"pos", "neg"}) {
        hnsigmatof.push_back(name + "/" + charge + "/hnsigmatof");
        hnsigmatpc.push_back(name + "/" + charge + "/hnsigmatpc");
        hnsigmatpctof.push_back(name + "/" + charge + "/hnsigmatpctof");
      }
    }
    // Define axes
    const AxisSpec centAxis{106, 0.0, 105.0, "Centrality (%)"};
    const AxisSpec ptAxis{150, 0.0, 15.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec etaAxis{40, -1.0, 1.0, "#eta"};
    const AxisSpec phiAxis{72, 0.0, 2 * M_PI, "#varphi (rad)"};
    const AxisSpec dcaAxis{200, -1.0, 1.0, "DCA (cm)"};
    const AxisSpec nclsTPCaxis{160, 0.0, 160.0, "TPC clusters"};
    const AxisSpec nsigmaAxis{100, -5.0, 5.0, "n#sigma"};
    const AxisSpec mass2Axis{1000, -1.0, 10.0, "m^{2} (GeV^{2}/c^{4})"};

    // Create histograms for each particle species
    for (const auto& name : nucleiNames) {
      for (const auto& charge : {"pos", "neg"}) {
        histos.add((name + "/" + charge + "/hPt").c_str(), (name + " " + charge + " #it{p}_{T}").c_str(), kTH1F, {ptAxis});
        histos.add((name + "/" + charge + "/hEta").c_str(), (name + " " + charge + " #eta").c_str(), kTH1F, {etaAxis});
        histos.add((name + "/" + charge + "/hPhi").c_str(), (name + " " + charge + " #varphi").c_str(), kTH1F, {phiAxis});
        histos.add((name + "/" + charge + "/hDCAxy").c_str(), (name + " " + charge + " DCA_{xy}").c_str(), kTH2F, {ptAxis, dcaAxis});
        histos.add((name + "/" + charge + "/hDCAz").c_str(), (name + " " + charge + " DCA_{z}").c_str(), kTH2F, {ptAxis, dcaAxis});
      }
    }

    for (const auto& hname : hnsigmatof) {
      histos.add(hname.c_str(), hname.c_str(), kTH3F, {centAxis, ptAxis, nsigmaAxis});
    }
    for (const auto& hname : hnsigmatpc) {
      histos.add(hname.c_str(), hname.c_str(), kTH3F, {centAxis, ptAxis, nsigmaAxis});
    }
    for (const auto& hname : hnsigmatpctof) {
      histos.add(hname.c_str(), hname.c_str(), kTH3F, {centAxis, ptAxis, nsigmaAxis});
    }
  }

  // Template function to process each nuclei species
  template <PID::ID id, typename T>
  void processNucleiTrack(const T& track)
  {
    // Get TPC n-sigma
    const auto& nsigmaTPC = o2::aod::pidutils::tpcNSigma<id>(track);
    if (std::abs(nsigmaTPC) > cutNsigmaTPC) {
      return;
    }

    float nSigmaTOF = -999.0f;
    // If TOF is not used or track doesn't have TOF, proceed
    if (!useTOF || (useTOF && track.hasTOF())) {
      if (useTOF && track.hasTOF()) {
        const auto& nsigmaTOF = o2::aod::pidutils::tofNSigma<id>(track);
        if (std::abs(nsigmaTOF) > cutNsigmaTOF) {
          return;
        }
        nSigmaTOF = nsigmaTOF;
      }

      // QA
      std::string name = PID::getName(id);
      auto charge = (track.charge() > 0) ? "pos" : "neg";
      histos.fill((name + "/" + charge + "/hPt").c_str(), track.pt());
      histos.fill((name + "/" + charge + "/hEta").c_str(), track.eta());
      histos.fill((name + "/" + charge + "/hPhi").c_str(), track.phi());
      histos.fill((name + "/" + charge + "/hDCAxy").c_str(), track.dcaXY());
      histos.fill((name + "/" + charge + "/hDCAz").c_str(), track.dcaZ());
      // nSigma
      histos.fill((name + "/" + charge + "/hnsigmatof").c_str(), collisionCentrality, track.pt(), nSigmaTOF);
      histos.fill((name + "/" + charge + "/hnsigmatpc").c_str(), collisionCentrality, track.pt(), nsigmaTPC);
      histos.fill((name + "/" + charge + "/hnsigmatpctof").c_str(), collisionCentrality, track.pt(), nSigmaTOF);
    }
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  bool isTrackSelected(CollisionType const&, TrackType const& track)
  {
    if (std::abs(track.eta()) > cutEta) {
      return false;
    }
    if (track.pt() < cutPtMin || track.pt() > cutPtMax) {
      return false;
    }
    if (std::abs(track.dcaXY()) > cutDCAxy) {
      return false;
    }
    if (std::abs(track.dcaZ()) > cutDCAz) {
      return false;
    }
    if (track.tpcNClsFound() < 70) {
      return false;
    }
    return true;
  }

  void process(CollisionCandidates const& collisions, TrackCandidates const& tracks)
  {
    for (auto& collision : collisions) {
      // collisionCentrality = GetCentrality(collision);
      collisionCentrality = 999.0f;
      for (auto& track : tracks) {
        if (!isTrackSelected<false>(collision, track)) {
          continue;
        }

        for (const auto& iPID : nucleiTypes) {
          processNucleiTrack<iPID>(track);
        }
      } // End of track loop
    }
  }
}; // end of nuclei spectra task

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleiSpectra>(cfgc)};
}
