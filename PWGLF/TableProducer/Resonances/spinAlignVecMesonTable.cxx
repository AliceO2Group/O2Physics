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

/// \file spinAlignVecMesonTable.cxx
/// \brief Derived-data producer for K*0 and phi spin alignment

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFSpinAlignVecMesonTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>

#include <cmath>
#include <cstdint>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

struct SpinAlignVecMesonTable {

  Produces<aod::SpinAlignEvents> eventTable;
  Produces<aod::SAKstarCands> kstarTable;
  Produces<aod::SAPhiCands> phiTable;
  Produces<aod::SAVMTracks> trackTable;

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Maximum |Vz| (cm)"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.0f, "Minimum FT0C centrality (%)"};
  Configurable<float> cfgCentMax{"cfgCentMax", 80.0f, "Maximum FT0C centrality (%)"};
  Configurable<int> cfgCutOccupancy{"cfgCutOccupancy", 2000, "Maximum track occupancy"};
  Configurable<bool> cfgUseRCT{"cfgUseRCT", true, "Apply RCT flag checker"};
  Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "RCT flag checker label"};

  Configurable<float> cfgCutPt{"cfgCutPt", 0.15f, "Minimum track pT (GeV/c)"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Maximum track |eta|"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "Maximum |DCAxy| (cm)"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.1f, "Maximum |DCAz| (cm)"};
  Configurable<int> cfgMinTPCCrossedRows{"cfgMinTPCCrossedRows", 70, "Minimum TPC crossed rows"};
  Configurable<int> cfgMinITSCls{"cfgMinITSCls", 5, "Minimum ITS clusters"};

  Configurable<float> cfgNSigmaTPCKa{"cfgNSigmaTPCKa", 3.0f, "Kaon: maximum |n-sigma TPC|"};
  Configurable<float> cfgNSigmaTOFKa{"cfgNSigmaTOFKa", 3.0f, "Kaon: maximum |n-sigma TOF|"};
  Configurable<float> cfgNSigmaTPCPi{"cfgNSigmaTPCPi", 3.0f, "Pion: maximum |n-sigma TPC|"};
  Configurable<float> cfgNSigmaTOFPi{"cfgNSigmaTOFPi", 3.0f, "Pion: maximum |n-sigma TOF|"};
  Configurable<bool> cfgRequireTOF{"cfgRequireTOF", false, "Require a TOF hit above cfgPtTOFRequired"};
  Configurable<float> cfgPtTOFRequired{"cfgPtTOFRequired", 0.5f, "pT (GeV/c) above which TOF is required"};

  Configurable<float> cfgKstarMassMin{"cfgKstarMassMin", 0.70f, "Minimum K pi mass (GeV/c^2)"};
  Configurable<float> cfgKstarMassMax{"cfgKstarMassMax", 1.10f, "Maximum K pi mass (GeV/c^2)"};
  Configurable<float> cfgPhiMassMin{"cfgPhiMassMin", 0.99f, "Minimum K K mass (GeV/c^2)"};
  Configurable<float> cfgPhiMassMax{"cfgPhiMassMax", 1.07f, "Maximum K K mass (GeV/c^2)"};
  Configurable<float> cfgPairRapidity{"cfgPairRapidity", 0.5f, "Maximum pair |y|"};
  Configurable<bool> cfgStoreTracks{"cfgStoreTracks", true, "Store the selected single kaon and pion tracks (SAVMTracks), needed for the v2 like-sign and mixed-event pairs (not for rho00)"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  RCTFlagsChecker rctChecker;

  static constexpr double MassKa = o2::constants::physics::MassKaonCharged;
  static constexpr double MassPi = o2::constants::physics::MassPionCharged;

  void init(InitContext&)
  {
    rctChecker.init(cfgRCTLabel.value, false, true);

    histos.add("hEvents", "Number of events;;Events", kTH1F, {{2, 0.0f, 2.0f}});
    histos.add("hCent", "Selected events;FT0C centrality (%);Events", kTH1F, {{100, 0.0f, 100.0f}});
    histos.add("hKstarMass", "K* candidates;M_{K#pi} (GeV/c^{2});Counts", kTH1F, {{200, 0.6f, 1.2f}});
    histos.add("hPhiMass", "#phi candidates;M_{KK} (GeV/c^{2});Counts", kTH1F, {{200, 0.98f, 1.08f}});
    histos.add("hPsiFT0C", "Selected events;FT0C centrality (%);#Psi_{2}^{FT0C}", kTH2F, {{8, 0.0f, 80.0f}, {100, -o2::constants::math::PIHalf, o2::constants::math::PIHalf}});
    if (cfgStoreTracks) {
      histos.add("hNTracks", "Stored kaon/pion tracks per event;N_{tracks};Events", kTH1F, {{500, 0.0f, 2000.0f}});
    }
  }

  Filter vertexFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centFilter = aod::cent::centFT0C >= cfgCentMin && aod::cent::centFT0C < cfgCentMax;
  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta && aod::track::pt > cfgCutPt;
  Filter dcaFilter = nabs(aod::track::dcaXY) < cfgCutDCAxy && nabs(aod::track::dcaZ) < cfgCutDCAz;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::EPCalibrationTables>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullKa, aod::pidTPCFullPi,
                                                  aod::pidTOFFullKa, aod::pidTOFFullPi>>;

  template <typename T>
  bool isGoodEvent(T const& collision)
  {
    if (!collision.sel8()) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (collision.trackOccupancyInTimeRange() >= cfgCutOccupancy) {
      return false;
    }
    if (cfgUseRCT && !rctChecker(collision)) {
      return false;
    }
    if (!collision.triggereventep()) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isGoodTrack(T const& track)
  {
    if (!track.isGlobalTrack()) {
      return false;
    }
    if (!track.isPVContributor()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < cfgMinTPCCrossedRows) {
      return false;
    }
    if (track.itsNCls() < cfgMinITSCls) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isKaon(T const& track)
  {
    if (std::abs(track.tpcNSigmaKa()) > cfgNSigmaTPCKa) {
      return false;
    }
    if (cfgRequireTOF && track.pt() > cfgPtTOFRequired && !track.hasTOF()) {
      return false;
    }
    if (track.hasTOF() && std::abs(track.tofNSigmaKa()) > cfgNSigmaTOFKa) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isPion(T const& track)
  {
    if (std::abs(track.tpcNSigmaPi()) > cfgNSigmaTPCPi) {
      return false;
    }
    if (cfgRequireTOF && track.pt() > cfgPtTOFRequired && !track.hasTOF()) {
      return false;
    }
    if (track.hasTOF() && std::abs(track.tofNSigmaPi()) > cfgNSigmaTOFPi) {
      return false;
    }
    return true;
  }

  template <typename T>
  ROOT::Math::PxPyPzMVector makePair(T const& pos, double massPos, T const& neg, double massNeg)
  {
    ROOT::Math::PxPyPzMVector p4Pos(pos.px(), pos.py(), pos.pz(), massPos);
    ROOT::Math::PxPyPzMVector p4Neg(neg.px(), neg.py(), neg.pz(), massNeg);
    return p4Pos + p4Neg;
  }

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks,
                   aod::BCsWithTimestamps const&)
  {
    histos.fill(HIST("hEvents"), 0.5);
    if (!isGoodEvent(collision)) {
      return;
    }
    histos.fill(HIST("hEvents"), 1.5);
    histos.fill(HIST("hCent"), collision.centFT0C());
    histos.fill(HIST("hPsiFT0C"), collision.centFT0C(), collision.psiFT0C());

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    eventTable(bc.runNumber(), collision.posZ(), collision.centFT0M(), collision.centFT0C(),
               collision.psiFT0C(), collision.psiFT0A(), collision.psiTPC(),
               collision.qFT0C(), collision.qFT0A(), collision.qTPC());
    auto eventIndex = eventTable.lastIndex();

    if (cfgStoreTracks) {
      int nStored = 0;
      for (const auto& track : tracks) {
        if (track.sign() == 0 || !isGoodTrack(track)) {
          continue;
        }
        const bool kaon = isKaon(track);
        const bool pion = isPion(track);
        if (!kaon && !pion) {
          continue;
        }
        trackTable(eventIndex, track.px(), track.py(), track.pz(), static_cast<int8_t>(track.sign()),
                   track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPi(),
                   track.hasTOF(), kaon, pion);
        ++nStored;
      }
      histos.fill(HIST("hNTracks"), nStored);
    }

    for (const auto& pos : tracks) {
      if (pos.sign() <= 0 || !isGoodTrack(pos)) {
        continue;
      }

      for (const auto& neg : tracks) {
        if (neg.sign() >= 0 || !isGoodTrack(neg)) {
          continue;
        }

        if (isKaon(pos) && isKaon(neg)) {
          auto phi = makePair(pos, MassKa, neg, MassKa);
          float mass = phi.M();
          if (mass > cfgPhiMassMin && mass < cfgPhiMassMax && std::abs(phi.Rapidity()) < cfgPairRapidity) {
            phiTable(eventIndex, mass,
                     static_cast<float>(phi.Px()), static_cast<float>(phi.Py()), static_cast<float>(phi.Pz()),
                     pos.px(), pos.py(), pos.pz(),
                     neg.px(), neg.py(), neg.pz(),
                     pos.tpcNSigmaKa(), pos.tofNSigmaKa(), pos.hasTOF(),
                     neg.tpcNSigmaKa(), neg.tofNSigmaKa(), neg.hasTOF());
            histos.fill(HIST("hPhiMass"), mass);
          }
        }

        if (isKaon(pos) && isPion(neg)) {
          auto kstar = makePair(pos, MassKa, neg, MassPi);
          float mass = kstar.M();
          if (mass > cfgKstarMassMin && mass < cfgKstarMassMax && std::abs(kstar.Rapidity()) < cfgPairRapidity) {
            kstarTable(eventIndex, static_cast<int8_t>(1), mass,
                       static_cast<float>(kstar.Px()), static_cast<float>(kstar.Py()), static_cast<float>(kstar.Pz()),
                       pos.px(), pos.py(), pos.pz(),
                       neg.px(), neg.py(), neg.pz(),
                       pos.tpcNSigmaKa(), pos.tofNSigmaKa(), pos.hasTOF(),
                       neg.tpcNSigmaPi(), neg.tofNSigmaPi(), neg.hasTOF());
            histos.fill(HIST("hKstarMass"), mass);
          }
        }

        if (isPion(pos) && isKaon(neg)) {
          auto kstar = makePair(pos, MassPi, neg, MassKa);
          float mass = kstar.M();
          if (mass > cfgKstarMassMin && mass < cfgKstarMassMax && std::abs(kstar.Rapidity()) < cfgPairRapidity) {
            kstarTable(eventIndex, static_cast<int8_t>(-1), mass,
                       static_cast<float>(kstar.Px()), static_cast<float>(kstar.Py()), static_cast<float>(kstar.Pz()),
                       pos.px(), pos.py(), pos.pz(),
                       neg.px(), neg.py(), neg.pz(),
                       pos.tpcNSigmaPi(), pos.tofNSigmaPi(), pos.hasTOF(),
                       neg.tpcNSigmaKa(), neg.tofNSigmaKa(), neg.hasTOF());
            histos.fill(HIST("hKstarMass"), mass);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(SpinAlignVecMesonTable, processData, "Fill K* and phi derived tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<SpinAlignVecMesonTable>(cfgc)};
}
