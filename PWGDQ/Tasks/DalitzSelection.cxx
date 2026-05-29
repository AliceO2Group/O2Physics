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
//
/// \file   DalitzSelection.cxx
/// \author Gauthier Legras, glegras@uni-muenster.de, gauthier.legras@cern.ch
/// \brief  Task to select electrons from dalitz decay
/// It can produce track and pair histograms for selected tracks
/// It creates a bitmap with selections to be stored during skimming
//
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <THashList.h>
#include <TList.h>
#include <TObjArray.h>
#include <TString.h>

#include <RtypesCore.h>

#include <chrono>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;

// using MyEvents = soa::Join<aod::Collisions, aod::EvSels>
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0As, aod::CentFT0Ms>;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                 aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullMu,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullMu,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

using MyBarrelTracksNoTOF = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                      aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullMu,
                                      aod::pidTPCFullKa, aod::pidTPCFullPr>;

// constexpr static uint32_t EventFillMap = VarManager::ObjTypes::Collision;
constexpr static uint32_t EventFillMapWithCent = VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent;
constexpr static uint32_t TrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t TrackFillMapNoTOF = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackTPCPID;

struct DalitzSelection {
  Produces<o2::aod::DalitzBits> dalitzbits;
  Preslice<MyBarrelTracks> perCollision = aod::track::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  // Configurables
  // cuts
  struct : ConfigurableGroup {
    Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandardNoINT7", "Event selection"};
    Configurable<std::string> fConfigDalitzTrackCuts{"cfgDalitzTrackCuts", "", "Dalitz track selection cuts for tag, separated by a comma"};
    Configurable<std::string> fConfigDalitzTrackCutsProbe{"cfgDalitzTrackCutsProbe", "", "Dalitz track selection cuts for probe, separated by a comma (if empty, use same selections as the tag)"};
    Configurable<std::string> fConfigDalitzPairCuts{"cfgDalitzPairCuts", "", "Dalitz pair selection cuts"};
    Configurable<std::string> fConfigTrackCutsJSON{"cfgTrackCutsJSON", "", "Additional list of barrel track cuts in JSON format"};
    Configurable<std::string> fConfigTrackCutsProbeJSON{"cfgTrackCutsProbeJSON", "", "Additional list of barrel track cuts for the probe in JSON format"};
    Configurable<std::string> fConfigPairCutsJSON{"cfgPairCutsJSON", "", "Additional list of barrel track cuts in JSON format"};
    Configurable<float> fConfigPtLow{"cfgLowPt", 0.1f, "Low pt cut for Dalitz tracks in the barrel"};
    Configurable<float> fConfigEtaCut{"cfgEtaCut", 0.9f, "Eta cut for Dalitz tracks in the barrel"};
    Configurable<float> fConfigTPCNSigLow{"cfgTPCNSigElLow", -3.f, "Low TPCNSigEl cut for Dalitz tracks in the barrel"};
    Configurable<float> fConfigTPCNSigHigh{"cfgTPCNSigElHigh", 3.f, "High TPCNsigEl cut for Dalitz tracks in the barrel"};
    Configurable<int> fConfigTrackSel{"cfgTrackSel", 0, "track selection requirement: 0: all, 1: reject ITS only, 2: reject TPC only, 3: reject ITS only and TPC only"};
  } fConfigCuts;

  // histograms
  struct : ConfigurableGroup {
    Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddPairHistogram{"cfgAddPairHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format for tracks"};
  } fConfigHistograms;

  // additional options
  struct : ConfigurableGroup {
    Configurable<bool> fConfigEnableLikeSign{"cfgEnableLikeSign", false, "Whether or not also add like-sign pairs (for studying combinatorial background which might contain misID or non-primary electrons)"};
    Configurable<bool> fQA{"cfgQA", true, "QA histograms"};
    Configurable<bool> fRemoveDoubleCounting{"cfgRemoveDoubleCounting", true, "If a track/pair is selected for several collisions, we still count it only once"};
    // Configurables for TPC post-calibration maps
    Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/i/iarsene/Calib/TPCpostCalib", "base path to the ccdb object"};
    Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
    Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas"};
    Configurable<bool> fUseRemoteField{"cfgUseRemoteField", true, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
    Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  } fConfigOptions;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::parameters::GRPMagField* grpmag = nullptr;
  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  Filter filterBarrelTrack = (o2::aod::track::pt >= fConfigCuts.fConfigPtLow) && (nabs(o2::aod::track::eta) <= fConfigCuts.fConfigEtaCut) && ((o2::aod::pidtpc::tpcNSigmaEl <= fConfigCuts.fConfigTPCNSigHigh && o2::aod::pidtpc::tpcNSigmaEl >= fConfigCuts.fConfigTPCNSigLow) || ((o2::aod::track::v001::detectorMap & 2) == 0)) && ((fConfigCuts.fConfigTrackSel != 1 && fConfigCuts.fConfigTrackSel != 3) || ((o2::aod::track::v001::detectorMap & 2) != 0)) && ((fConfigCuts.fConfigTrackSel != 2 && fConfigCuts.fConfigTrackSel != 3) || ((o2::aod::track::v001::detectorMap & 1) != 0));

  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  OutputObj<TList> fStatsList{"Statistics"};  //! skimming statistics

  std::map<int, uint8_t> fTrackmap;       // whether it is selected with symmetric or tag cut
  std::map<int, uint8_t> fTrackmapProbe;  // whether it is selected with probe cut
  std::map<int, uint8_t> fDalitzmap;      // whether it is selected as dalitz decay daughter with symmetric or tag cut
  std::map<int, uint8_t> fDalitzmapProbe; // whether it is selected as dalitz decay daughter with probe cut

  // maps to remove ambiguities
  std::map<int, uint8_t> fDalitzmapAmbiguity;
  std::map<int, uint8_t> fDalitzmapProbeAmbiguity;
  std::map<std::pair<int, int>, uint8_t> fAmbiguousPairs;

  AnalysisCompositeCut* fEventCut;
  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<AnalysisCompositeCut> fTrackCutsProbe;
  std::vector<AnalysisCompositeCut> fPairCuts;

  bool fIsTagAndProbe; // whether we are doing tag and probe, or just symmetric cuts
  bool fIsOutputRequested;
  bool fSkipEvent; // speed up by skipping next step of event if no track/pair is selected

  HistogramManager* fHistMan;

  void init(o2::framework::InitContext& context)
  {
    fIsTagAndProbe = false;
    VarManager::SetDefaultVarNames();

    // Event cuts
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigCuts.fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

    // Barrel track cuts
    // Tag cuts or symmetric cuts
    TString cutNamesStr = fConfigCuts.fConfigDalitzTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // extra cuts via JSON
    TString addTrackCutsStr = fConfigCuts.fConfigTrackCutsJSON.value;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (const auto& t : addTrackCuts) {
        fTrackCuts.push_back(*dynamic_cast<AnalysisCompositeCut*>(t));
      }
    }

    // Probe cuts
    TString cutNamesProbeStr = fConfigCuts.fConfigDalitzTrackCutsProbe.value;
    TString addTrackCutsProbeStr = fConfigCuts.fConfigTrackCutsProbeJSON.value;
    if (!cutNamesProbeStr.IsNull() && (cutNamesProbeStr.CompareTo(cutNamesStr) != 0 || addTrackCutsProbeStr.CompareTo(addTrackCutsStr) != 0)) {
      std::unique_ptr<TObjArray> objArray(cutNamesProbeStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCutsProbe.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
      fIsTagAndProbe = true;
    }
    // extra cuts via JSON
    if (addTrackCutsProbeStr != "" && (cutNamesProbeStr.CompareTo(cutNamesStr) != 0 || addTrackCutsProbeStr.CompareTo(addTrackCutsStr) != 0)) {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsProbeStr.Data());
      for (const auto& t : addTrackCuts) {
        fTrackCutsProbe.push_back(*dynamic_cast<AnalysisCompositeCut*>(t));
      }
      fIsTagAndProbe = true;
    }

    // Pair cuts
    TString cutNamesPairStr = fConfigCuts.fConfigDalitzPairCuts.value;
    if (!cutNamesPairStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesPairStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // extra cuts via JSON
    TString addPairCutsStr = fConfigCuts.fConfigPairCutsJSON.value;
    if (addPairCutsStr != "") {
      std::vector<AnalysisCut*> addPairCuts = dqcuts::GetCutsFromJSON(addPairCutsStr.Data());
      for (const auto& t : addPairCuts) {
        fPairCuts.push_back(*dynamic_cast<AnalysisCompositeCut*>(t));
      }
    }

    if (fTrackCuts.size() != fPairCuts.size() || (fTrackCuts.size() != fTrackCutsProbe.size() && fIsTagAndProbe)) {
      LOGF(fatal, "YOU SHOULD PROVIDE THE SAME NUMBER OF TRACK AND PAIR CUTS");
    }

    // CCDB configuration
    if (fConfigOptions.fConfigComputeTPCpostCalib) {
      fCCDB->setURL(fConfigOptions.fConfigCcdbUrl.value);
      fCCDB->setCaching(true);
      fCCDB->setLocalObjectValidityChecking();
      // Not later than now objects
      fCCDB->setCreatedNotAfter(fConfigOptions.fConfigNoLaterThan.value);
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Create the histogram class names to be added to the histogram manager
    TString histClasses = "";

    if (fConfigOptions.fQA) {
      auto trackCut = fTrackCuts.begin();
      int iCut = 0;
      for (auto pairCut = fPairCuts.begin(); pairCut != fPairCuts.end(); pairCut++, trackCut++, iCut++) {
        if (fIsTagAndProbe) {
          histClasses += Form("TrackBarrelTag_%s_%s_%s;", (*trackCut).GetName(), fTrackCutsProbe.at(iCut).GetName(), (*pairCut).GetName());
          histClasses += Form("TrackBarrelProbe_%s_%s_%s;", (*trackCut).GetName(), fTrackCutsProbe.at(iCut).GetName(), (*pairCut).GetName());
          histClasses += Form("Pair_%s_%s_%s;", (*trackCut).GetName(), fTrackCutsProbe.at(iCut).GetName(), (*pairCut).GetName());
          if (fConfigOptions.fConfigEnableLikeSign) {
            histClasses += Form("PairLS_%s_%s_%s;", (*trackCut).GetName(), fTrackCutsProbe.at(iCut).GetName(), (*pairCut).GetName());
          }
        } else {
          histClasses += Form("TrackBarrel_%s_%s;", (*trackCut).GetName(), (*pairCut).GetName());
          histClasses += Form("Pair_%s_%s;", (*trackCut).GetName(), (*pairCut).GetName());
          if (fConfigOptions.fConfigEnableLikeSign) {
            histClasses += Form("PairLS_%s_%s;", (*trackCut).GetName(), (*pairCut).GetName());
          }
        }
      }
    }

    // Define histograms

    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      fHistMan->AddHistClass(classStr.Data());

      if (classStr.Contains("Event")) {
        dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", "");
      }

      TString histTrackName = fConfigHistograms.fConfigAddTrackHistogram.value;
      if (classStr.Contains("Track")) {
        dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histTrackName);
        dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "software-trigger", histTrackName);
      }

      TString histPairName = fConfigHistograms.fConfigAddPairHistogram.value;
      if (classStr.Contains("Pair")) {
        dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "pair", histPairName);
      }
    }

    // Additional histogram via the JSON configurable
    TString addHistsStr = fConfigHistograms.fConfigAddJSONHistograms.value;
    if (addHistsStr != "") {
      dqhistograms::AddHistogramsFromJSON(fHistMan, addHistsStr.Data());
    }

    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);

    // Dalitz selection statistics: one bin for each (track,pair) selection
    TH1I* histTracks = new TH1I("TrackStats", "Dalitz selection statistics", fPairCuts.size() * (1 + fIsTagAndProbe) + 1, -0.5, fPairCuts.size() * (1 + fIsTagAndProbe) + 0.5);
    histTracks->GetXaxis()->SetBinLabel(1, "Events selected");
    auto trackCut = fTrackCuts.begin();
    int icut = 1;
    for (auto pairCut = fPairCuts.begin(); pairCut != fPairCuts.end(); pairCut++, trackCut++, icut++) {
      if (fIsTagAndProbe) {
        histTracks->GetXaxis()->SetBinLabel(2 * icut, Form("%s_%s_%s tag", (*trackCut).GetName(), fTrackCutsProbe.at(icut - 1).GetName(), (*pairCut).GetName()));
        histTracks->GetXaxis()->SetBinLabel(2 * icut + 1, Form("%s_%s_%s probe", (*trackCut).GetName(), fTrackCutsProbe.at(icut - 1).GetName(), (*pairCut).GetName()));
      } else {
        histTracks->GetXaxis()->SetBinLabel(icut + 1, Form("%s_%s", (*trackCut).GetName(), (*pairCut).GetName()));
      }
    }
    if (fConfigOptions.fQA) {
      fStatsList->Add(histTracks);
    }

    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    // autodetect whether the dalitz bits are requested
    fIsOutputRequested = false;
    auto& workflows = context.services().template get<o2::framework::RunningWorkflowInfo const>();
    // autodetect this table in other devices
    for (o2::framework::DeviceSpec const& device : workflows.devices) {
      // Check if this device subscribed to the dalitz table
      for (auto const& input : device.inputs) {
        if (o2::framework::DataSpecUtils::partialMatch(input.matcher, o2::header::DataOrigin("AOD"))) {
          std::string tableName = "DalitzBits";
          if (input.matcher.binding == tableName) {
            LOGF(info, "Device %s has subscribed to %s", device.name, "DalitzBits");
            fIsOutputRequested = true;
          }
        }
      }
    }
  }

  template <bool isReassoc, uint32_t TTrackFillMap, typename TFullTracks, typename TTracks, typename TEvent>
  void runTrackSelection(TTracks const& tracksBarrel, TEvent const& collision)
  {

    fSkipEvent = true;
    for (const auto& track1 : tracksBarrel) {
      uint8_t filterMap = uint8_t(0);
      uint8_t filterMapProbe = uint8_t(0);
      int fullTrackIdx = 0;
      if constexpr (isReassoc) {
        auto const& fullTrack = track1.template track_as<TFullTracks>();
        // Basic filter cuts which cannot be applied directly in case of reassoc
        if (fullTrack.pt() < fConfigCuts.fConfigPtLow) {
          continue;
        }
        if (std::fabs(fullTrack.eta()) > fConfigCuts.fConfigEtaCut) {
          continue;
        }
        if (fullTrack.hasTPC() && (fullTrack.tpcNSigmaEl() < fConfigCuts.fConfigTPCNSigLow || fullTrack.tpcNSigmaEl() > fConfigCuts.fConfigTPCNSigHigh)) {
          continue;
        }
        if ((fConfigCuts.fConfigTrackSel == 1 || fConfigCuts.fConfigTrackSel == 3) && !fullTrack.hasTPC()) {
          continue;
        }
        if ((fConfigCuts.fConfigTrackSel == 2 || fConfigCuts.fConfigTrackSel == 3) && !fullTrack.hasITS()) {
          continue;
        }
        VarManager::FillTrack<TTrackFillMap>(fullTrack);
        VarManager::FillTrackCollision<TTrackFillMap>(fullTrack, collision);
        fullTrackIdx = fullTrack.globalIndex();
      } else {
        VarManager::FillTrack<TTrackFillMap>(track1);
        fullTrackIdx = track1.globalIndex();
      }
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint8_t(1) << i);
          fSkipEvent = false;
        }
      }
      if (fIsTagAndProbe) {
        i = 0;
        for (auto cut = fTrackCutsProbe.begin(); cut != fTrackCutsProbe.end(); ++cut, ++i) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            filterMapProbe |= (uint8_t(1) << i);
          }
        }
      }
      if (filterMap) {
        fTrackmap[fullTrackIdx] = filterMap;
      }
      if (filterMapProbe) {
        fTrackmapProbe[fullTrackIdx] = filterMapProbe;
      }
    } // end loop over tracks
  }

  template <bool isReassoc, int TPairType, uint32_t TTrackFillMap, typename TFullTracks, typename TTracks, typename TEvent>
  void runDalitzPairing(TTracks const& tracks1, TTracks const& tracks2, TEvent collision)
  {
    if (isReassoc & fConfigOptions.fRemoveDoubleCounting) {
      fDalitzmapAmbiguity.clear();
      fDalitzmapProbeAmbiguity.clear();
    }
    fSkipEvent = true;

    for (const auto& track1 : tracks1) {
      int trackIdx1;
      bool isTrack1P;

      if constexpr (isReassoc) {
        auto const& fullTrack1 = track1.template track_as<TFullTracks>();
        trackIdx1 = fullTrack1.globalIndex();
        isTrack1P = fullTrack1.sign() > 0;
      } else {
        trackIdx1 = track1.globalIndex();
        isTrack1P = track1.sign() > 0;
      }

      if (!fTrackmap[trackIdx1]) {
        continue;
      }

      for (const auto& track2 : tracks2) {
        int trackIdx2;
        bool isLikeSign;
        if constexpr (isReassoc) {
          auto const& fullTrack2 = track2.template track_as<TFullTracks>();
          trackIdx2 = fullTrack2.globalIndex();
          isLikeSign = (isTrack1P && (fullTrack2.sign() > 0)) || (!isTrack1P && (fullTrack2.sign() < 0));
        } else {
          trackIdx2 = track2.globalIndex();
          isLikeSign = (isTrack1P && (track2.sign() > 0)) || (!isTrack1P && (track2.sign() < 0));
        }

        if (trackIdx1 == trackIdx2) {
          continue;
        }
        if (!fIsTagAndProbe && trackIdx1 >= trackIdx2) {
          continue;
        }
        if (!fConfigOptions.fConfigEnableLikeSign && isLikeSign) {
          continue;
        }

        uint8_t twoTracksFilterMap = fTrackmap[trackIdx1] & (fIsTagAndProbe ? fTrackmapProbe[trackIdx2] : fTrackmap[trackIdx2]);
        if (!twoTracksFilterMap) {
          continue;
        }

        // pairing
        if constexpr (isReassoc) {
          auto const& fullTrack1 = track1.template track_as<TFullTracks>();
          auto const& fullTrack2 = track2.template track_as<TFullTracks>();
          VarManager::FillPair<TPairType, TTrackFillMap>(fullTrack1, fullTrack2);
        } else {
          VarManager::FillPair<TPairType, TTrackFillMap>(track1, track2);
        }

        // Fill pair selection map and fill pair histogram
        int icut = 0;
        auto trackCut = fTrackCuts.begin();
        for (auto pairCut = fPairCuts.begin(); pairCut != fPairCuts.end(); pairCut++, trackCut++, icut++) {
          if (!(twoTracksFilterMap & (uint8_t(1) << icut))) {
            continue;
          }
          if ((*pairCut).IsSelected(VarManager::fgValues)) {
            fSkipEvent = false;
            bool isPairAlreadySelected = false;
            if (isReassoc & fConfigOptions.fRemoveDoubleCounting) {
              // if we remove double counting and the pair is already selected, we don't fill the histograms
              std::pair<int, int> iPair(trackIdx1, trackIdx2);
              if (fAmbiguousPairs.find(iPair) != fAmbiguousPairs.end()) {
                if (fAmbiguousPairs[iPair] & (static_cast<uint8_t>(1) << icut)) { // if this pair is already stored with this cut
                  isPairAlreadySelected = true;
                } else {
                  fAmbiguousPairs[iPair] |= static_cast<uint8_t>(1) << icut;
                }
              } else {
                fAmbiguousPairs[iPair] = static_cast<uint8_t>(1) << icut;
              }
            }
            if (!isLikeSign) {
              if (isReassoc & fConfigOptions.fRemoveDoubleCounting) {
                // if we want to remove double-counting, we fill separately the map for this event and for previous ones, in order to check if it was previously selected
                fDalitzmapAmbiguity[trackIdx1] |= (uint8_t(1) << icut);
              } else {
                fDalitzmap[trackIdx1] |= (uint8_t(1) << icut);
              }
              if (fIsTagAndProbe) {
                if (isReassoc & fConfigOptions.fRemoveDoubleCounting) {
                  fDalitzmapProbeAmbiguity[trackIdx1] |= (uint8_t(1) << icut);
                } else {
                  fDalitzmapProbe[trackIdx2] |= (uint8_t(1) << icut);
                }
                if (fConfigOptions.fQA && !isPairAlreadySelected) {
                  fHistMan->FillHistClass(Form("Pair_%s_%s_%s", (*trackCut).GetName(), fTrackCutsProbe.at(icut).GetName(), (*pairCut).GetName()), VarManager::fgValues);
                }
              } else {
                if (isReassoc & fConfigOptions.fRemoveDoubleCounting) {
                  fDalitzmapAmbiguity[trackIdx1] |= (uint8_t(1) << icut);
                } else {
                  fDalitzmap[trackIdx2] |= (uint8_t(1) << icut);
                }
                if (fConfigOptions.fQA && !isPairAlreadySelected) {
                  fHistMan->FillHistClass(Form("Pair_%s_%s", (*trackCut).GetName(), (*pairCut).GetName()), VarManager::fgValues);
                }
              }
            } else {
              if (fConfigOptions.fQA && !isPairAlreadySelected) {
                fHistMan->FillHistClass(fIsTagAndProbe ? Form("PairLS_%s_%s_%s", (*trackCut).GetName(), fTrackCutsProbe.at(icut).GetName(), (*pairCut).GetName()) : Form("PairLS_%s_%s", (*trackCut).GetName(), (*pairCut).GetName()), VarManager::fgValues);
              }
            } // end if like-sign
          } // end if isSelected
        } // end cut loop
      }
    } // end of tracksP,N loop

    // Fill Hists
    if (fConfigOptions.fQA && !fSkipEvent) {
      for (const auto& track : tracks1) {
        uint8_t filterMap;
        uint8_t filterMapProbe;
        if constexpr (isReassoc) {
          auto const& fullTrack = track.template track_as<TFullTracks>();
          filterMap = fDalitzmap[fullTrack.globalIndex()];
          filterMapProbe = fDalitzmapProbe[fullTrack.globalIndex()];
          if (fConfigOptions.fRemoveDoubleCounting) {
            // we remove track selections which were already selected before
            uint8_t currentFilterMap = fDalitzmapAmbiguity[fullTrack.globalIndex()];
            uint8_t currentFilterMapProbe = fDalitzmapProbeAmbiguity[fullTrack.globalIndex()];
            filterMap = currentFilterMap & ~filterMap;
            filterMapProbe = currentFilterMapProbe & ~filterMapProbe;
            // Then only we can fill the final map
            fDalitzmap[fullTrack.globalIndex()] |= fDalitzmapAmbiguity[fullTrack.globalIndex()];
            fDalitzmapProbe[fullTrack.globalIndex()] |= fDalitzmapProbeAmbiguity[fullTrack.globalIndex()];
          }

          if (!filterMap && !filterMapProbe) {
            continue;
          }

          VarManager::FillTrack<TTrackFillMap>(fullTrack);
          // The caveat here is that we only fill the DCA to the first selected collision
          VarManager::FillTrackCollision<TTrackFillMap>(fullTrack, collision);
        } else {
          filterMap = fDalitzmap[track.globalIndex()];
          filterMapProbe = fDalitzmapProbe[track.globalIndex()];

          if (!filterMap && !filterMapProbe) {
            continue;
          }

          VarManager::FillTrack<TTrackFillMap>(track);
        }

        int icut = 0;
        auto trackCut = fTrackCuts.begin();
        for (auto pairCut = fPairCuts.begin(); pairCut != fPairCuts.end(); pairCut++, trackCut++, icut++) {
          if (filterMap & (uint8_t(1) << icut)) {
            reinterpret_cast<TH1I*>(fStatsList->At(0))->Fill(fIsTagAndProbe ? 2 * icut + 1 : icut + 1);
            fHistMan->FillHistClass(fIsTagAndProbe ? Form("TrackBarrelTag_%s_%s_%s", (*trackCut).GetName(), fTrackCutsProbe.at(icut).GetName(), (*pairCut).GetName()) : Form("TrackBarrel_%s_%s", (*trackCut).GetName(), (*pairCut).GetName()), VarManager::fgValues);
          }
          if (filterMapProbe & (uint8_t(1) << icut)) {
            reinterpret_cast<TH1I*>(fStatsList->At(0))->Fill(2 * icut + 2);
            fHistMan->FillHistClass(Form("TrackBarrelProbe_%s_%s_%s", (*trackCut).GetName(), fTrackCutsProbe.at(icut).GetName(), (*pairCut).GetName()), VarManager::fgValues);
          }
        }
      }
    }
  }

  void initNewRun(int64_t timestamp)
  {

    VarManager::ResetValues(0, VarManager::kNRunWiseVariables);

    // We setup the magnetic field, because the conversion rejection cut might depend on it
    float magField = 0.;
    if (fConfigOptions.fUseRemoteField.value) {
      grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigOptions.grpmagPath, timestamp);
      if (grpmag != nullptr) {
        magField = grpmag->getNominalL3Field();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", timestamp);
      }
    } else {
      magField = fConfigOptions.fConfigMagField.value;
    }
    LOGF(info, "setting mag field to %f", magField);
    if (magField == 0.) {
      LOGF(fatal, "magnetic field not set correctly, please check");
    }
    VarManager::SetMagneticField(magField);

    if (fConfigOptions.fConfigComputeTPCpostCalib) {
      auto calibList = fCCDB->getForTimeStamp<TList>(fConfigOptions.fConfigCcdbPathTPC.value, timestamp);
      VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
      VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
      VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
      VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
      VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
      VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
    }
  }

  void processFullTracks(MyEventsWithCent const& collisions, aod::BCsWithTimestamps const&, soa::Filtered<MyBarrelTracks> const& filteredTracks, MyBarrelTracks const& tracks)
  {
    const int pairType = VarManager::kDecayToEE;
    fDalitzmap.clear();
    fDalitzmapProbe.clear();

    for (const auto& collision : collisions) {
      fTrackmap.clear();
      fTrackmapProbe.clear();
      VarManager::ResetValues(VarManager::kNRunWiseVariables, VarManager::kNBarrelTrackVariables);
      VarManager::FillEvent<EventFillMapWithCent>(collision);
      bool isEventSelected = fEventCut->IsSelected(VarManager::fgValues);

      if (isEventSelected) {
        reinterpret_cast<TH1I*>(fStatsList->At(0))->Fill(0);

        auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

        if (fCurrentRun != bc.runNumber()) {
          initNewRun(bc.timestamp());
          fCurrentRun = bc.runNumber();
        }

        auto groupedFilteredTracks = filteredTracks.sliceBy(perCollision, collision.globalIndex());
        runTrackSelection<false, TrackFillMap, MyBarrelTracks>(groupedFilteredTracks, nullptr);
        if (!fSkipEvent) {
          runDalitzPairing<false, pairType, TrackFillMap, MyBarrelTracks>(groupedFilteredTracks, groupedFilteredTracks, nullptr);
        }
      }
    }

    if (fIsOutputRequested) {
      for (const auto& track : tracks) { // Fill dalitz bits
        dalitzbits(fIsTagAndProbe ? fDalitzmapProbe[track.globalIndex()] : fDalitzmap[track.globalIndex()]);
      }
    }
  }

  void processFullTracksWithAssoc(MyEventsWithCent const& collisions, aod::BCsWithTimestamps const&, MyBarrelTracks const& tracks, TrackAssoc const& trackAssocs)
  {

    const int pairType = VarManager::kDecayToEE;
    fDalitzmap.clear();
    fDalitzmapProbe.clear();
    fDalitzmapAmbiguity.clear();
    fDalitzmapProbeAmbiguity.clear();
    fAmbiguousPairs.clear();

    for (const auto& collision : collisions) {
      fTrackmap.clear();
      fTrackmapProbe.clear();
      VarManager::ResetValues(VarManager::kNRunWiseVariables, VarManager::kNBarrelTrackVariables);
      VarManager::FillEvent<EventFillMapWithCent>(collision);
      bool isEventSelected = fEventCut->IsSelected(VarManager::fgValues);

      if (isEventSelected) {

        reinterpret_cast<TH1I*>(fStatsList->At(0))->Fill(0);

        auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

        if (fCurrentRun != bc.runNumber()) {
          initNewRun(bc.timestamp());
          fCurrentRun = bc.runNumber();
        }

        auto groupedTracksAssoc = trackAssocs.sliceBy(trackIndicesPerCollision, collision.globalIndex());
        runTrackSelection<true, TrackFillMap, MyBarrelTracks>(groupedTracksAssoc, collision);
        if (!fSkipEvent) {
          runDalitzPairing<true, pairType, TrackFillMap, MyBarrelTracks>(groupedTracksAssoc, groupedTracksAssoc, collision);
        }
      }
    }

    if (fIsOutputRequested) {
      for (const auto& track : tracks) { // Fill dalitz bits
        dalitzbits(fIsTagAndProbe ? fDalitzmapProbe[track.globalIndex()] : fDalitzmap[track.globalIndex()]);
      }
    }
  }

  void processFullTracksWithAssocNoTOF(MyEventsWithCent const& collisions, aod::BCsWithTimestamps const&, MyBarrelTracksNoTOF const& tracks, TrackAssoc const& trackAssocs)
  {

    const int pairType = VarManager::kDecayToEE;
    fDalitzmap.clear();
    fDalitzmapProbe.clear();
    fDalitzmapAmbiguity.clear();
    fDalitzmapProbeAmbiguity.clear();
    fAmbiguousPairs.clear();

    for (const auto& collision : collisions) {
      fTrackmap.clear();
      fTrackmapProbe.clear();
      VarManager::ResetValues(VarManager::kNRunWiseVariables, VarManager::kNBarrelTrackVariables);
      VarManager::FillEvent<EventFillMapWithCent>(collision);
      bool isEventSelected = fEventCut->IsSelected(VarManager::fgValues);

      if (isEventSelected) {

        reinterpret_cast<TH1I*>(fStatsList->At(0))->Fill(0);

        auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

        if (fCurrentRun != bc.runNumber()) {
          initNewRun(bc.timestamp());
          fCurrentRun = bc.runNumber();
        }

        auto groupedTracksAssoc = trackAssocs.sliceBy(trackIndicesPerCollision, collision.globalIndex());
        runTrackSelection<true, TrackFillMapNoTOF, MyBarrelTracksNoTOF>(groupedTracksAssoc, collision);
        if (!fSkipEvent) {
          runDalitzPairing<true, pairType, TrackFillMapNoTOF, MyBarrelTracksNoTOF>(groupedTracksAssoc, groupedTracksAssoc, collision);
        }
      }
    }

    if (fIsOutputRequested) {
      for (const auto& track : tracks) { // Fill dalitz bits
        dalitzbits(fIsTagAndProbe ? fDalitzmapProbe[track.globalIndex()] : fDalitzmap[track.globalIndex()]);
      }
    }
  }

  void processDummy(MyEventsWithCent const&)
  {
  }

  PROCESS_SWITCH(DalitzSelection, processFullTracksWithAssocNoTOF, "Run Dalitz selection on AO2D tables with reassociation and without TOF PID", false);
  PROCESS_SWITCH(DalitzSelection, processFullTracksWithAssoc, "Run Dalitz selection on AO2D tables with reassociation", false);
  PROCESS_SWITCH(DalitzSelection, processFullTracks, "Run Dalitz selection on AO2D tables", false);
  PROCESS_SWITCH(DalitzSelection, processDummy, "Do nothing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzSelection>(cfgc)};
}
