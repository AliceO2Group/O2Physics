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

/*
namespace o2::aod
{
namespace dppbpbfilter
{
DECLARE_SOA_COLUMN(IsDQBarrelSelected, isDQBarrelSelected, uint32_t);
} // namespace dqpbpbfilter

DECLARE_SOA_TABLE(DQBarrelTrackCuts, "AOD", "DQBARRELCUTS", dqpbpbfilter::IsDQBarrelSelected);
} // namespace o2::aod
*/
namespace o2::aod
{
namespace dqppfilter
{
DECLARE_SOA_COLUMN(IsDQBarrelSelected, isDQBarrelSelected, uint32_t);
} // namespace dqppfilter

DECLARE_SOA_TABLE(DQBarrelTrackCuts, "AOD", "DQBARRELCUTS", dqppfilter::IsDQBarrelSelected);
} // namespace o2::aod

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr>;

using MyBarrelTracksSelected = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                         aod::pidTPCFullEl, aod::pidTPCFullPi,
                                         aod::pidTPCFullKa, aod::pidTPCFullPr,
                                         aod::pidTOFFullEl, aod::pidTOFFullPi,
                                         aod::pidTOFFullKa, aod::pidTOFFullPr,
                                         aod::DQBarrelTrackCuts>;
using MyMuons = aod::FwdTracks;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct DQBarrelTrackSelection {
  Produces<aod::DQBarrelTrackCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigCuts{"cfgBarrelTrackCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/i/iarsene/Calib/TPCpostCalib", "base path to the ccdb object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<TString> fCutHistNames;

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  void init(o2::framework::InitContext&)
  {
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        AnalysisCompositeCut* cut = dqcuts::GetCompositeCut(objArray->At(icut)->GetName());
        if (cut) {
          fTrackCuts.push_back(*cut);
        } else {
          LOGF(fatal, "Invalid barrel track cut provided: %s", objArray->At(icut)->GetName());
        }
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      TString cutNames = "TrackBarrel_BeforeCuts;";
      for (auto& cut : fTrackCuts) {
        cutNames += Form("TrackBarrel_%s;", cut.GetName());
        fCutHistNames.push_back(Form("TrackBarrel_%s", cut.GetName()));
      }

      DefineHistograms(fHistMan, cutNames.Data());     // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());

      // CCDB configuration
      if (fConfigComputeTPCpostCalib) {
        fCCDB->setURL(fConfigCcdbUrl.value);
        fCCDB->setCaching(true);
        fCCDB->setLocalObjectValidityChecking();
        // Not later than now objects
        fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
      }
    }
  }

  // Templated function instantianed for all of the process functions
  template <uint32_t TTrackFillMap, typename TTracks>
  void runTrackSelection(aod::BCsWithTimestamps const& bcs, TTracks const& tracksBarrel)
  {
    auto bc = bcs.begin(); // check just the first bc to get the run number
    if (fConfigComputeTPCpostCalib && fCurrentRun != bc.runNumber()) {
      auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, bc.timestamp());
      VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
      VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
      VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
      VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
      VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
      VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
      fCurrentRun = bc.runNumber();
    }

    uint32_t filterMap = uint32_t(0);
    trackSel.reserve(tracksBarrel.size());

    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
    for (auto& track : tracksBarrel) {
      filterMap = uint32_t(0);
      if (!track.has_collision()) {
        trackSel(uint32_t(0));
      } else {
        VarManager::FillTrack<TTrackFillMap>(track);
        if (fConfigQA) {
          fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
        }
        int i = 0;
        for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            filterMap |= (uint32_t(1) << i);
            if (fConfigQA) {
              fHistMan->FillHistClass(fCutHistNames[i].Data(), VarManager::fgValues);
            }
          }
        }
        trackSel(filterMap);
      }
    } // end loop over tracks
  }

  void processSelection(aod::BCsWithTimestamps const& bcs, MyBarrelTracks const& tracks)
  {
    runTrackSelection<gkTrackFillMap>(bcs, tracks);
  }
  void processDummy(MyBarrelTracks&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQBarrelTrackSelection, processSelection, "Run barrel track selection", true);
  PROCESS_SWITCH(DQBarrelTrackSelection, processDummy, "Dummy function", false);
};

struct DQFilterPbPbTask {
  Produces<aod::DQEventFilter> eventFilter;
  OutputObj<TH1I> fStats{"Statistics"};

  Configurable<std::string> fConfigBarrelSelections{"cfgBarrelSels", "jpsiPID1:2:5", "<track-cut>:<nmin>:<nmax>,[<track-cut>:<nmin>:<nmax>],..."};
  Configurable<int> fConfigNDtColl{"cfgNDtColl", 4, "Number of standard deviations to consider in BC range"};
  Configurable<int> fConfigMinNBCs{"cfgMinNBCs", 7, "Minimum number of BCs to consider in BC range"};
  Configurable<float> fConfigFV0AmpLimit{"cfgFV0AmpLimit", 0, "FV0 amplitude limit for double gap event selection"};
  Configurable<float> fConfigFT0AAmpLimit{"cfgFT0AAmpLimit", 0, "FT0A amplitude limit for double gap event selection"};
  Configurable<float> fConfigFT0CAmpLimit{"cfgFT0CAmpLimit", 0, "FT0C amplitude limit for double gap event selection"};
  Configurable<float> fConfigFDDAAmpLimit{"cfgFDDAAmpLimit", 0, "FDDA amplitude limit for double gap event selection"};
  Configurable<float> fConfigFDDCAmpLimit{"cfgFDDCAmpLimit", 0, "FDDC amplitude limit for double gap event selection"};
  Configurable<bool> fConfigUseFV0{"cfgUseFV0", true, "Whether to use FV0 for DG veto"};
  Configurable<bool> fConfigUseFT0{"cfgUseFT0", true, "Whether to use FT0 for DG veto"};
  Configurable<bool> fConfigUseFDD{"cfgUseFDD", true, "Whether to use FDD for DG veto"};

  Filter filterBarrelTrackSelected = aod::dqppfilter::isDQBarrelSelected > uint32_t(0);

  int fNBarrelCuts;                   // number of barrel selections
  std::vector<int> fBarrelNminTracks; // minimal number of tracks in barrel
  std::vector<int> fBarrelNmaxTracks; // maximal number of tracks in barrel

  // Helper function for selecting DG events
  bool isEventDG(MyEvents::iterator const& collision, MyBCs const& bcs, MyBarrelTracksSelected const& tracks, MyMuons const& muons,
                 aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds,
                 std::vector<float> FITAmpLimits, int nDtColl, int minNBCs, bool useFV0, bool useFT0, bool useFDD)
  {
    // Find BC associated with collision
    if (!collision.has_foundBC()) {
      return 0;
    }
    // foundBCId is stored in EvSels
    auto bc = collision.foundBC_as<MyBCs>();

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

    // find slice of BCs table with BC in [minBC, maxBC]
    int64_t minBCId = bc.globalIndex();
    int64_t maxBCId = bc.globalIndex();

    // lower limit
    if (bc.globalBC() < minBC) {
      while (bc != bcs.end() && bc.globalBC() < minBC) {
        ++bc;
        minBCId = bc.globalIndex();
      }
    } else {
      while (bc.globalIndex() > 0 && bc.globalBC() >= minBC) {
        minBCId = bc.globalIndex();
        --bc;
      }
    }
    // upper limit
    if (bc.globalBC() < maxBC) {
      while (bc != bcs.end() && bc.globalBC() <= maxBC) {
        maxBCId = bc.globalIndex();
        ++bc;
      }
    } else {
      while (bc.globalIndex() > 0 && bc.globalBC() > maxBC) {
        --bc;
        maxBCId = bc.globalIndex();
      }
    }
    MyBCs bcrange{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, (uint64_t)minBCId};

    // DG condition: Check FIT activity in BC range
    bool cleanFIT = true;

    for (auto const& bc : bcrange) {
      if (useFV0) {
        if (bc.has_foundFV0()) {
          auto fv0a = fv0as.iteratorAt(bc.foundFV0Id());
          float FV0Amplitude = 0;
          for (auto amp : fv0a.amplitude()) {
            FV0Amplitude += amp;
          }
          if (FV0Amplitude > FITAmpLimits[0]) {
            cleanFIT = false;
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
          if (FT0AAmplitude > FITAmpLimits[1] || FT0CAmplitude > FITAmpLimits[2]) {
            cleanFIT = false;
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
          if (FDDAAmplitude > FITAmpLimits[3] || FDDCAmplitude > FITAmpLimits[4]) {
            cleanFIT = false;
          }
        }
      }
    }

    // DG condition: Veto on activiy in either forward or barrel region. Save info on both cases
    bool muonsEmpty = true;
    if (muons.size() > 0) {
      muonsEmpty = false;
    }
    bool barrelEmpty = true;
    if (tracks.size() > 0) {
      barrelEmpty = false;
    }

    // Compute decision. For now only check if FIT is empty and there are no forward tracks
    if (cleanFIT && muonsEmpty) {
      return 1;
    } else {
      return 0;
    }
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
  }

  template <uint32_t TEventFillMap>
  void runFilterPbPb(MyEvents::iterator const& collision, MyBCs const& bcs, MyBarrelTracksSelected const& tracks, MyMuons const& muons,
                     aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {
    fStats->Fill(-2.0);

    std::vector<float> FITAmpLimits = {fConfigFV0AmpLimit, fConfigFT0AAmpLimit, fConfigFT0CAmpLimit, fConfigFDDAAmpLimit, fConfigFDDCAmpLimit};
    bool isDG = isEventDG(collision, bcs, tracks, muons, ft0s, fv0as, fdds, FITAmpLimits, fConfigNDtColl, fConfigMinNBCs, fConfigUseFV0, fConfigUseFT0, fConfigUseFDD);
    fStats->Fill(-1.0, isDG);

    std::vector<int> objCountersBarrel(fNBarrelCuts, 0); // init all counters to zero
    // Count the number of barrel tracks fulfilling each cut
    for (auto track : tracks) {
      for (int i = 0; i < fNBarrelCuts; ++i) {
        if (track.isDQBarrelSelected() & (uint32_t(1) << i)) {
          objCountersBarrel[i] += 1;
        }
      }
    }

    // Compute event filter
    uint64_t filter = 0;
    filter |= isDG;
    for (int i = 0; i < fNBarrelCuts; i++) {
      if (objCountersBarrel[i] >= fBarrelNminTracks[i] && objCountersBarrel[i] <= fBarrelNmaxTracks[i]) {
        filter |= (uint64_t(1) << i + 1);
        fStats->Fill(static_cast<float>(i));
      }
    }
    eventFilter(filter);
  }

  void processFilterPbPb(MyEvents::iterator const& collision, MyBCs const& bcs,
                         soa::Filtered<MyBarrelTracksSelected> const& tracks, MyMuons const& muons,
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
    adaptAnalysisTask<DQBarrelTrackSelection>(cfgc),
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
