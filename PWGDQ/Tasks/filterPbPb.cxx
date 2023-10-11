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
#include <iostream>
#include <fstream>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "CommonConstants/LHCConstants.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace std;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr>;

using MyMuons = aod::FwdTracks;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct DQFilterPbPbTask {
  Produces<aod::DQEventFilter> eventFilter;
  OutputObj<TH1I> fStats{"Statistics"};
  OutputObj<TH1F> fIsEventDGOutcome{TH1F("Filter outcome", "Filter outcome", 14, -1.5, 6.5)};

  Configurable<std::string> fConfigBarrelSelections{"cfgBarrelSels", "jpsiPID1:2:5", "<track-cut>:<nmin>:<nmax>,[<track-cut>:<nmin>:<nmax>],..."};
  Configurable<int> fConfigNDtColl{"cfgNDtColl", 4, "Number of standard deviations to consider in BC range"};
  Configurable<int> fConfigMinNBCs{"cfgMinNBCs", 7, "Minimum number of BCs to consider in BC range"};
  Configurable<int> fConfigMinNPVCs{"cfgMinNPVCs", 2, "Minimum number of PV contributors"};
  Configurable<int> fConfigMaxNPVCs{"cfgMaxNPVCs", 5, "Maximum number of PV contributors"};
  Configurable<float> fConfigMaxFITTime{"cfgMaxFITTime", 4, "Maximum time in FIT"};
  Configurable<float> fConfigFV0AmpLimit{"cfgFV0AmpLimit", 0, "FV0 amplitude limit for event selection"};
  Configurable<float> fConfigFT0AAmpLimit{"cfgFT0AAmpLimit", 0, "FT0A amplitude limit for event selection"};
  Configurable<float> fConfigFT0CAmpLimit{"cfgFT0CAmpLimit", 0, "FT0C amplitude limit for event selection"};
  Configurable<float> fConfigFDDAAmpLimit{"cfgFDDAAmpLimit", 0, "FDDA amplitude limit for event selection"};
  Configurable<float> fConfigFDDCAmpLimit{"cfgFDDCAmpLimit", 0, "FDDC amplitude limit for event selection"};
  Configurable<bool> fConfigUseFV0{"cfgUseFV0", true, "Whether to use FV0 for veto"};
  Configurable<bool> fConfigUseFT0{"cfgUseFT0", true, "Whether to use FT0 for veto"};
  Configurable<bool> fConfigUseFDD{"cfgUseFDD", true, "Whether to use FDD for veto"};
  Configurable<std::string> fConfigFITSides{"cfgFITSides", "both", "both, either, A, C, neither"};
  Configurable<bool> fConfigVetoForward{"cfgVetoForward", true, "Whether to veto on forward tracks"};
  Configurable<bool> fConfigVetoBarrel{"cfgVetoBarrel", false, "Whether to veto on barrel tracks"};

  // Filter filterBarrelTrackSelected = aod::dqppfilter::isDQBarrelSelected > uint32_t(0);

  int fNBarrelCuts;                   // number of barrel selections
  std::vector<int> fBarrelNminTracks; // minimal number of tracks in barrel
  std::vector<int> fBarrelNmaxTracks; // maximal number of tracks in barrel

  int FITVetoSides = -1; // Integer to encode which side(s) of the FIT to use for veto
  std::vector<std::string> FITVetoSidesOptions = {"both", "either", "A", "C", "neither"};

  // Helper function for selecting DG events
  template <typename TEvent, typename TBCs, typename TTracks, typename TMuons>
  int isEventDG(TEvent const& collision, TBCs const& bcs, TTracks const& tracks, TMuons const& muons,
                aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds,
                std::vector<float> FITAmpLimits, int nDtColl, int minNBCs, int minNPVCs, int maxNPVCs, float maxFITTime,
                bool useFV0, bool useFT0, bool useFDD, int FITSide, bool doVetoFwd, bool doVetoBarrel)
  {
    fIsEventDGOutcome->Fill(0., 1.);
    // Find BC associated with collision
    if (!collision.has_foundBC()) {
      return -1;
    }
    // foundBCId is stored in EvSels
    auto bc = collision.template foundBC_as<TBCs>();

    // Obtain slice of compatible BCs
    uint64_t mostProbableBC = bc.globalBC();
    uint64_t meanBC = mostProbableBC + std::lround(collision.collisionTime() / o2::constants::lhc::LHCBunchSpacingNS);
    int deltaBC = std::ceil(collision.collisionTimeRes() / o2::constants::lhc::LHCBunchSpacingNS * nDtColl);
    if (deltaBC < minNBCs) {
      deltaBC = minNBCs;
    }

    // Range of BCs to consider
    uint64_t minBC = (uint64_t)deltaBC < meanBC ? meanBC - (uint64_t)deltaBC : 0;
    uint64_t maxBC = meanBC + (uint64_t)deltaBC;

    int slicemin = 0;
    int slicemax = 0;
    // Check if there is overlap between acceptable and possible BCs
    if (maxBC > bcs.iteratorAt(0).globalBC() && minBC < bcs.iteratorAt(bcs.size() - 1).globalBC()) {
      // find slice of BCs table with BC in [minBC, maxBC]
      int moveCount = 0;
      int64_t minBCId = bc.globalIndex();
      int64_t maxBCId = bc.globalIndex();
      // lower limit
      if (bc.globalBC() < minBC) {
        while (bc != bcs.end() && bc.globalBC() < minBC) {
          ++bc;
          ++moveCount;
          minBCId = bc.globalIndex();
        }
      } else {
        while (bc.globalIndex() > 0 && bc.globalBC() >= minBC) {
          minBCId = bc.globalIndex();
          --bc;
          --moveCount;
        }
      }
      // upper limit
      if (bc.globalBC() < maxBC) {
        while (bc != bcs.end() && bc.globalBC() <= maxBC) {
          maxBCId = bc.globalIndex();
          ++bc;
          ++moveCount;
        }
      } else {
        while (bc.globalIndex() > 0 && bc.globalBC() > maxBC) {
          --bc;
          --moveCount;
          maxBCId = bc.globalIndex();
        }
      }
      // reset bc
      bc.moveByIndex(-moveCount);
      // Create BC slice
      slicemin = minBCId;
      slicemax = maxBCId - minBCId + 1;
    }
    MyBCs bcrange{{bcs.asArrowTable()->Slice(slicemin, slicemax)}, (uint16_t)slicemin};
    // DG condition: Check FIT activity in BC range
    bool isSideAClean = true;
    bool isSideCClean = true;

    for (auto const& bc : bcrange) {
      if (useFV0) {
        if (bc.has_foundFV0()) {
          auto fv0a = fv0as.iteratorAt(bc.foundFV0Id());
          float FV0Amplitude = 0;
          for (auto amp : fv0a.amplitude()) {
            FV0Amplitude += amp;
          }
          float FV0Time = std::abs(fv0a.time());
          if (FV0Amplitude > FITAmpLimits[0] || FV0Time > maxFITTime) {
            isSideAClean = false;
          }
        }
      }
      if (useFT0) {
        if (bc.has_foundFT0()) {
          auto ft0 = ft0s.iteratorAt(bc.foundFT0Id());
          float FT0AAmplitude = 0;
          float FT0CAmplitude = 0;
          for (auto amp : ft0.amplitudeA()) {
            FT0AAmplitude += amp;
          }
          for (auto amp : ft0.amplitudeC()) {
            FT0CAmplitude += amp;
          }
          float FT0ATime = std::abs(ft0.timeA());
          float FT0CTime = std::abs(ft0.timeC());
          if (FT0AAmplitude > FITAmpLimits[1] || FT0ATime > maxFITTime) {
            isSideAClean = false;
          }
          if (FT0CAmplitude > FITAmpLimits[2] || FT0CTime > maxFITTime) {
            isSideCClean = false;
          }
        }
      }
      if (useFDD) {
        if (bc.has_foundFDD()) {
          auto fdd = fdds.iteratorAt(bc.foundFDDId());
          float FDDAAmplitude = 0;
          float FDDCAmplitude = 0;
          for (auto amp : fdd.chargeA()) {
            FDDAAmplitude += amp;
          }
          for (auto amp : fdd.chargeC()) {
            FDDCAmplitude += amp;
          }
          float FDDATime = std::abs(fdd.timeA());
          float FDDCTime = std::abs(fdd.timeC());
          if (FDDAAmplitude > FITAmpLimits[3] || FDDATime > maxFITTime) {
            isSideAClean = false;
          }
          if (FDDCAmplitude > FITAmpLimits[4] || FDDCTime > maxFITTime) {
            isSideCClean = false;
          }
        }
      }
    }

    // Compute FIT decision
    bool FITDecision = 0;
    if (FITSide == 0) {
      FITDecision = isSideAClean && isSideCClean;
    } else if (FITSide == 1) {
      FITDecision = isSideAClean ^ isSideCClean;
    } else if (FITSide == 2) {
      FITDecision = isSideAClean;
    } else if (FITSide == 3) {
      FITDecision = isSideCClean;
    } else if (FITSide == 4) {
      FITDecision = 1;
    }
    if (!FITDecision) {
      return 1;
    }

    // Veto on activiy in either forward or barrel region
    if (doVetoFwd) {
      for (auto& muon : muons) {
        // Only care about muons with good timing (MID)
        if (muon.trackType() == 0 || muon.trackType() == 3) {
          return 2;
        }
      }
    }
    if (doVetoBarrel) {
      if (tracks.size() > 0) {
        return 3;
      }
    }

    // No global tracks which are not vtx tracks
    for (auto& track : tracks) {
      if (track.isGlobalTrack() && !track.isPVContributor()) {
        return 4;
      }
    }

    // Number of primary vertex contributors
    if (collision.numContrib() < minNPVCs || collision.numContrib() > maxNPVCs) {
      return 5;
    }

    // If we made it here, the event passed
    return 0;
  }

  void DefineCuts()
  {
    TString barrelSelsStr = fConfigBarrelSelections.value;
    std::unique_ptr<TObjArray> objArray(barrelSelsStr.Tokenize(","));
    fNBarrelCuts = objArray->GetEntries();
    if (fNBarrelCuts) {
      for (int icut = 0; icut < fNBarrelCuts; ++icut) {
        TString selStr = objArray->At(icut)->GetName();
        std::unique_ptr<TObjArray> sel(selStr.Tokenize(":"));
        if (sel->GetEntries() != 3) {
          continue;
        } else {
          fBarrelNminTracks.push_back(std::atoi(sel->At(1)->GetName()));
          fBarrelNmaxTracks.push_back(std::atoi(sel->At(2)->GetName()));
        }
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);

    // setup the Stats histogram
    fStats.setObject(new TH1I("Statistics", "Stats for DQ triggers", fNBarrelCuts + 2, -2.5, -0.5 + fNBarrelCuts));
    fStats->GetXaxis()->SetBinLabel(1, "Events inspected");
    fStats->GetXaxis()->SetBinLabel(2, "Events selected");
    if (fNBarrelCuts) {
      for (int ib = 3; ib < 3 + fNBarrelCuts; ib++) {
        fStats->GetXaxis()->SetBinLabel(ib, objArray->At(ib - 3)->GetName());
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    for (int i = 0; i < 5; i++) {
      if (fConfigFITSides.value == FITVetoSidesOptions[i]) {
        FITVetoSides = i;
      }
    }
    if (FITVetoSides == -1) {
      LOGF(fatal, "Invalid choice of FIT side(s) for veto: %s", fConfigFITSides.value);
    }
  }

  template <uint32_t TEventFillMap>
  void runFilterPbPb(MyEvents::iterator const& collision, MyBCs const& bcs, MyBarrelTracks const& tracks, MyMuons const& muons,
                     aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {
    fStats->Fill(-2.0);

    std::vector<float> FITAmpLimits = {fConfigFV0AmpLimit, fConfigFT0AAmpLimit, fConfigFT0CAmpLimit, fConfigFDDAAmpLimit, fConfigFDDCAmpLimit};
    int filterOutcome = isEventDG(collision, bcs, tracks, muons, ft0s, fv0as, fdds,
                                  FITAmpLimits, fConfigNDtColl, fConfigMinNBCs, fConfigMinNPVCs, fConfigMaxNPVCs, fConfigMaxFITTime,
                                  fConfigUseFV0, fConfigUseFT0, fConfigUseFDD, FITVetoSides, fConfigVetoForward, fConfigVetoBarrel);
    fIsEventDGOutcome->Fill(filterOutcome + 1, 1);

    // Don't need info on filter outcome anymore; reduce to a boolean
    bool isDG = !filterOutcome;
    fStats->Fill(-1.0, isDG);

    // Compute event filter
    uint64_t filter = 0;
    filter |= isDG;
    eventFilter(filter);
  }

  void processFilterPbPb(MyEvents::iterator const& collision, MyBCs const& bcs,
                         MyBarrelTracks const& tracks, MyMuons const& muons,
                         aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {
    runFilterPbPb<gkEventFillMap>(collision, bcs, tracks, muons, ft0s, fv0as, fdds);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQFilterPbPbTask, processFilterPbPb, "Run filter task", true);
  PROCESS_SWITCH(DQFilterPbPbTask, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQFilterPbPbTask>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    if (classStr.Contains("Track")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca");
    }
  }
}
