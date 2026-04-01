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
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//   Configurable workflow for running several DQ or other PWG analyses

#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/VarManager.h"

#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <THashList.h>
#include <TPDGCode.h>
#include <TString.h>

#include <RtypesCore.h>

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::McCollisionLabels>;
using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA,
                                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                        aod::McTrackLabels>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMapWithMults = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionMult | VarManager::ObjTypes::CollisionMultExtra;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;

// Forward declarations
void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups); // defines histograms for all tasks

struct AnalysisEnergyCorrelator {
  OutputObj<THashList> fOutputList{"output"};

  struct : ConfigurableGroup { // Event selection configurables
    Configurable<std::vector<double>> fConfigZBins{"cfgZBins", std::vector<double>{-10.0, -5.0, 0.0, 5.0, 10.0}, "Z vertex bins for mixing"};
    Configurable<std::vector<double>> fConfigMultBins{"cfgMultBins", std::vector<double>{0.0, 20.0, 40.0, 60.0, 100.0}, "Multiplicity bins for mixing"};
    Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
    Configurable<std::string> fConfigEventCutsJSON{"cfgEventCutsJSON", "", "Additional event cuts in JSON"};
    Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Event histograms"};
    Configurable<std::string> fConfigAddEventMCHistogram{"cfgAddEventMCHistogram", "generator", "MC Event histograms"};
    Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 5, "Event mixing pool depth"};
    Configurable<float> fConfigEventfilterVtz{"cfgEventfilterVtz", 10.0, "Event filter Vtz"};
    Configurable<bool> fConfigEventQA{"cfgEventQA", false, "If true, fill Event QA histograms"};
  } fConfigEventOptions;

  struct : ConfigurableGroup { // Track selection configurables
    Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "electronSelection1_ionut", "Comma separated list of barrel track cuts for electrons"};
    Configurable<std::string> fConfigTrackCutsJSON{"cfgTrackCutsJSON", "", "Additional track cuts in JSON"};
    Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Track histograms"};
    Configurable<bool> fConfigTrackQA{"cfgTrackQA", false, "If true, fill Track QA histograms"};
  } fConfigTrackOptions;

  struct : ConfigurableGroup { // Pair selection configurables
    Configurable<float> fConfigJpsiMassMin{"cfgJpsiMassMin", 2.8, "J/psi mass minimum"};
    Configurable<float> fConfigJpsiMassMax{"cfgJpsiMassMax", 3.3, "J/psi mass maximum"};
    Configurable<float> fConfigJpsiPtMin{"cfgJpsiPtMin", 0.0, "J/psi pt minimum"};
    Configurable<float> fConfigJpsiPtMax{"cfgJpsiPtMax", 100.0, "J/psi pt maximum"};
    Configurable<float> fConfigJpsiRapMax{"cfgJpsiRapMax", 0.9, "J/psi rapidity maximum"};
    Configurable<std::string> fConfigAddSEPHistogram{"cfgAddSEPHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> recSignals{"cfgMCRecSignals", "eeFromJpsi", "Comma separated list of MC signals (reconstructed)"};
    Configurable<std::string> recSignalsJSON{"cfgMCRecSignalsJSON", "", "Comma separated list of MC signals (reconstructed) via JSON"};
  } fConfigPairOptions;

  struct : ConfigurableGroup { // Dilepton selection configurables
    Configurable<std::string> fConfigHadronCuts{"cfgHadronCuts", "NoPID", "Comma separated list of hadron track cuts"};
    Configurable<std::string> fConfigHadronCutsJSON{"cfgHadronCutsJSON", "", "Additional hadron cuts in JSON"};
    Configurable<bool> fConfigApplyMassEC{"cfgApplyMassEC", false, "Apply fit mass for sideband for the energy correlator study"};
    Configurable<std::vector<int>> fConfigSavelessevents{"cfgSavelessevents", std::vector<int>{1, 0}, "Save less events for the energy correlator study"};
    Configurable<std::vector<float>> fConfigTransRange{"cfgTransRange", std::vector<float>{0.333333, 0.666667}, "Transverse region for the energy correlstor analysis"};
    Configurable<std::string> fConfigAddDileptonHadronHistogram{"cfgAddDileptonHadronHistogram", "", "Dilepton-hadron histograms"};
    Configurable<std::string> fConfigMCRecSignals{"cfgMCRecDileptonHadronSignals", "", "Comma separated list of MC signals (reconstructed)"};
    Configurable<std::string> fConfigMCGenSignals{"cfgMCGenDileptonHadronSignals", "", "Comma separated list of MC signals (generated)"};
    Configurable<std::string> fConfigMCRecSignalsJSON{"cfgMCRecDileptonHadronSignalsJSON", "", "Additional list of MC signals (reconstructed) via JSON"};
    Configurable<std::string> fConfigMCGenSignalsJSON{"cfgMCGenDileptonHadronSignalsJSON", "", "Comma separated list of MC signals (generated) via JSON"};
    Configurable<float> fConfigMCGenHadronEtaAbs{"cfgMCGenHadronEtaAbs", 0.9f, "eta abs range for the hadron"};
    Configurable<float> fConfigMCGenHadronPtMin{"cfgMCGenHadronPtMin", 0.1f, "minimum pt for the hadron"};
    Configurable<bool> fConfigContainlepton{"cfgContainlepton", false, "If true, require the hadron to contain the lepton in its decay tree for the energy correlator study"};
    Configurable<bool> fConfigUsePionMass{"cfgUsePionMass", false, "If true, use pion mass for the hadron in the energy correlator study"};
  } fConfigDileptonHadronOptions;

  // Histogram configurables
  Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON"};

  // CCDB configurables
  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB url"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "CCDB timestamp"};

  // Member variables
  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;

  AnalysisCompositeCut* fEventCut = nullptr;
  std::vector<AnalysisCompositeCut*> fTrackCuts;  // Electron cuts
  std::vector<AnalysisCompositeCut*> fHadronCuts; // Hadron cuts

  std::vector<TString> fTrackCutNames;
  std::vector<TString> fHadronCutNames;
  std::vector<TString> fHistNamesReco;

  std::map<int, std::vector<TString>> fTrackHistNames;
  std::map<int, std::vector<TString>> fBarrelHistNamesMCmatched;
  std::map<int64_t, bool> fSelMap;

  std::vector<MCSignal*> fRecMCSignals; // MC signals for reconstructed pairs
  std::vector<MCSignal*> fGenMCSignals;
  std::vector<MCSignal*> fRecMCTripleSignals; // MC signals for reconstructed triples

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  int fCurrentRun = -1;

  // Preslice for association table
  Preslice<aod::TrackAssoc> preslice = aod::track_association::collisionId;

  using MixingBinning = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultNTracksPV>;
  std::unique_ptr<MixingBinning> fMixingBinning;

  void init(o2::framework::InitContext& context)
  {
    std::vector<double> zBins = fConfigEventOptions.fConfigZBins.value;
    zBins.insert(zBins.begin(), VARIABLE_WIDTH);

    std::vector<double> multBins = fConfigEventOptions.fConfigMultBins.value;
    multBins.insert(multBins.begin(), VARIABLE_WIDTH);

    fMixingBinning = std::make_unique<MixingBinning>(std::array<std::vector<double>, 2>{zBins, multBins}, true);

    bool isBarrelME = context.mOptions.get<bool>("processBarrelMixedEvent");
    bool isMCGen_energycorrelators = context.mOptions.get<bool>("processMCGenEnergyCorrelators") || context.mOptions.get<bool>("processMCGenEnergyCorrelatorsPion");
    bool isMCGen_energycorrelatorsME = context.mOptions.get<bool>("processMCGenEnergyCorrelatorsME") || context.mOptions.get<bool>("processMCGenEnergyCorrelatorsPionME");

    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }
    VarManager::SetDefaultVarNames();

    // Setup Event Cuts
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventOptions.fConfigEventCuts.value;
    if (eventCutStr != "") {
      AnalysisCut* cut = dqcuts::GetAnalysisCut(eventCutStr.Data());
      if (cut != nullptr) {
        fEventCut->AddCut(cut);
      }
    }
    // Additional cuts via JSON
    TString eventCutJSONStr = fConfigEventOptions.fConfigEventCutsJSON.value;
    if (eventCutJSONStr != "") {
      std::vector<AnalysisCut*> jsonCuts = dqcuts::GetCutsFromJSON(eventCutJSONStr.Data());
      for (auto& cutIt : jsonCuts) {
        fEventCut->AddCut(cutIt);
      }
    }

    // Setup Electron Track Cuts
    TString trackCutStr = fConfigTrackOptions.fConfigTrackCuts.value;
    if (!trackCutStr.IsNull()) {
      std::unique_ptr<TObjArray> objArrayTrack(trackCutStr.Tokenize(","));
      for (int icut = 0; icut < objArrayTrack->GetEntries(); ++icut) {
        fTrackCuts.push_back(dqcuts::GetCompositeCut(objArrayTrack->At(icut)->GetName()));
        fTrackCutNames.push_back(objArrayTrack->At(icut)->GetName());
      }
    }

    TString trackCutsJSON = fConfigTrackOptions.fConfigTrackCutsJSON.value;
    if (trackCutsJSON != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(trackCutsJSON.Data());
      for (auto& t : addTrackCuts) {
        fTrackCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
        fTrackCutNames.push_back(t->GetName());
        trackCutStr += Form(",%s", t->GetName());
      }
    }

    // Setting the MC rec signal names
    TString sigNamesStr = fConfigPairOptions.recSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 2) { // NOTE: 2-prong signals required
          continue;
        }
        fRecMCSignals.push_back(sig);
      }
    }

    TString sigNamesHadronStr = fConfigDileptonHadronOptions.fConfigMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigTripleArray(sigNamesHadronStr.Tokenize(","));
    if (!sigNamesHadronStr.IsNull()) {
      for (int isig = 0; isig < objRecSigTripleArray->GetEntries(); ++isig) {
        MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigTripleArray->At(isig)->GetName());
        if (sig) {
          if (sig->GetNProngs() != 3) {
            LOG(fatal) << "Signal at reconstructed level requested (" << sig->GetName() << ") " << "does not have 3 prongs! Fix it";
          }
          fRecMCTripleSignals.push_back(sig);
        } else {
          LOG(fatal) << "Signal at reconstructed level requested (" << objRecSigTripleArray->At(isig)->GetName() << ") " << "could not be retrieved from the library! -> skipped";
        }
      }
    }

    // Add the MCSignals from the JSON config
    TString addMCSignalsStr = fConfigPairOptions.recSignalsJSON.value;
    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() != 2) { // NOTE: only 2 prong signals
          continue;
        }
        fRecMCSignals.push_back(mcIt);
      }
    }

    // Add the reco MCSignals from the JSON config
    TString addMCTripleSignalsStr = fConfigDileptonHadronOptions.fConfigMCRecSignalsJSON.value;
    if (addMCTripleSignalsStr != "") {
      std::vector<MCSignal*> addMCTripleSignals = dqmcsignals::GetMCSignalsFromJSON(addMCTripleSignalsStr.Data());
      for (auto& mcIt : addMCTripleSignals) {
        if (mcIt->GetNProngs() != 3) {
          LOG(fatal) << "Signal at reconstructed level requested (" << mcIt->GetName() << ") " << "does not have 3 prongs! Fix it";
        }
        fRecMCTripleSignals.push_back(mcIt);
      }
    }

    // Setup Hadron Cuts
    TString hadronCutStr = fConfigDileptonHadronOptions.fConfigHadronCuts.value;
    if (!hadronCutStr.IsNull()) {
      std::unique_ptr<TObjArray> objArrayHadron(hadronCutStr.Tokenize(","));
      for (int icut = 0; icut < objArrayHadron->GetEntries(); ++icut) {
        TString cutName = objArrayHadron->At(icut)->GetName();
        fHadronCuts.push_back(dqcuts::GetCompositeCut(cutName.Data()));
        fHadronCutNames.push_back(cutName);
      }
    }
    TString hadronCutsJSON = fConfigDileptonHadronOptions.fConfigHadronCutsJSON.value;
    if (hadronCutsJSON != "") {
      std::vector<AnalysisCut*> addHadronCuts = dqcuts::GetCutsFromJSON(hadronCutsJSON.Data());
      for (auto& t : addHadronCuts) {
        fHadronCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
        fHadronCutNames.push_back(t->GetName());
        hadronCutStr += Form(",%s", t->GetName());
      }
    }

    // Add histogram classes for each specified MCsignal at the generator level
    // TODO: create a std::vector of hist classes to be used at Fill time, to avoid using Form in the process function
    TString sigGenNamesStr = fConfigDileptonHadronOptions.fConfigMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        fGenMCSignals.push_back(sig);
      }
    }

    // Add the MCSignals from the JSON config
    TString addMCSignalsGenStr = fConfigDileptonHadronOptions.fConfigMCGenSignalsJSON.value;
    if (addMCSignalsGenStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsGenStr.Data());
      for (auto& mcIt : addMCSignals) {
        if (mcIt->GetNProngs() > 2) { // NOTE: only 2 prong signals
          continue;
        }
        fGenMCSignals.push_back(mcIt);
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars);

    fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(true);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Setup Histograms
    if (fConfigEventOptions.fConfigEventQA) {
      DefineHistograms(fHistMan, "TimeFrameStats;Event_BeforeCuts;Event_AfterCuts;", fConfigEventOptions.fConfigAddEventHistogram.value.data());
      DefineHistograms(fHistMan, "EventsMC", fConfigEventOptions.fConfigAddEventMCHistogram.value.data()); // mc
    }

    if (fConfigTrackOptions.fConfigTrackQA) {
      TString histClasses = "AssocsBarrel_BeforeCuts;";
      // Configure histogram classes for each track cut;
      // Add histogram classes for each track cut and for each requested MC signal (reconstructed tracks with MC truth)
      for (auto& cut : fTrackCuts) {
        TString nameStr = Form("AssocsBarrel_%s", cut->GetName());
        fHistNamesReco.push_back(nameStr);
        histClasses += Form("%s;", nameStr.Data());
      }
      DefineHistograms(fHistMan, histClasses.Data(), fConfigTrackOptions.fConfigAddTrackHistogram.value.data());
    }
    TString histNames = "";
    // check that the barrel track cuts array required in this task is not empty
    if (!trackCutStr.IsNull()) {
      // tokenize and loop over the barrel cuts produced by the barrel track selection task
      std::unique_ptr<TObjArray> objArray(trackCutStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        TString tempStr = objArray->At(icut)->GetName();
        // if the current barrel selection cut is required in this task, then switch on the corresponding bit in the mask
        // and assign histogram directories

        // assign the pair hist directories for the current cut
        std::vector<TString> names = {
          Form("PairsBarrelSEPM_%s", objArray->At(icut)->GetName()),
          Form("PairsBarrelSEPP_%s", objArray->At(icut)->GetName()),
          Form("PairsBarrelSEMM_%s", objArray->At(icut)->GetName())};
        for (auto& n : names) {
          histNames += Form("%s;", n.Data());
        }
        fTrackHistNames[icut] = names;
        // assign hist directories for the MC matched pairs for each (track cut,MCsignal) combination
        if (!sigNamesStr.IsNull()) {
          for (size_t isig = 0; isig < fRecMCSignals.size(); isig++) {
            auto sig = fRecMCSignals.at(isig);
            names = {
              Form("PairsBarrelSEPM_%s_%s", objArray->At(icut)->GetName(), sig->GetName()),
              Form("PairsBarrelSEPP_%s_%s", objArray->At(icut)->GetName(), sig->GetName()),
              Form("PairsBarrelSEMM_%s_%s", objArray->At(icut)->GetName(), sig->GetName())};
            for (auto& n : names) {
              histNames += Form("%s;", n.Data());
            }
            fBarrelHistNamesMCmatched.try_emplace(icut * fRecMCSignals.size() + isig, names);
          } // end loop over MC signals
        }
      }
      DefineHistograms(fHistMan, histNames.Data(), fConfigPairOptions.fConfigAddSEPHistogram.value.data());
    }

    for (size_t iCutTrack = 0; iCutTrack < fTrackCutNames.size(); iCutTrack++) {
      for (size_t iCutHadron = 0; iCutHadron < fHadronCutNames.size(); iCutHadron++) {
        DefineHistograms(fHistMan, Form("DileptonTrack_%s_%s", fTrackCutNames[iCutTrack].Data(), fHadronCutNames[iCutHadron].Data()), fConfigDileptonHadronOptions.fConfigAddDileptonHadronHistogram.value.data());
        for (auto& sig : fRecMCTripleSignals) {
          DefineHistograms(fHistMan, Form("DileptonTrackMCMatched_%s_%s_%s", fTrackCutNames[iCutTrack].Data(), fHadronCutNames[iCutHadron].Data(), sig->GetName()), fConfigDileptonHadronOptions.fConfigAddDileptonHadronHistogram.value.data());
          if (isBarrelME) {
            DefineHistograms(fHistMan, Form("DileptonTrackMCMatchedME_%s_%s_%s", fTrackCutNames[iCutTrack].Data(), fHadronCutNames[iCutHadron].Data(), sig->GetName()), fConfigDileptonHadronOptions.fConfigAddDileptonHadronHistogram.value.data());
          }
        }
      }
    }

    for (auto& sig : fGenMCSignals) {
      if (sig->GetNProngs() == 1) {
        if (isMCGen_energycorrelators) {
          DefineHistograms(fHistMan, Form("MCTruthGenSel_%s", sig->GetName()), "");
          DefineHistograms(fHistMan, Form("MCTruthEenergyCorrelators_%s", sig->GetName()), "");
          DefineHistograms(fHistMan, Form("MCTruthEenergyCorrelators_Pion_%s", sig->GetName()), "");
        }
        if (isMCGen_energycorrelatorsME) {
          DefineHistograms(fHistMan, Form("MCTruthEenergyCorrelatorsME_%s", sig->GetName()), "");
          DefineHistograms(fHistMan, Form("MCTruthEenergyCorrelatorsME_Pion_%s", sig->GetName()), "");
        }
      }
    }

    dqhistograms::AddHistogramsFromJSON(fHistMan, fConfigAddJSONHistograms.value.c_str()); // aditional histograms via JSON

    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    fCCDB->setURL(fConfigCcdbUrl.value);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
  }

  template <bool MixedEvent, uint32_t TTrackFillMap, typename TTrack1, typename TTrack2, typename THadron, typename TEvent>
  void runDileptonHadron(TTrack1 const& track1, TTrack2 const& track2, int iEleCut,
                         THadron const& hadron, TEvent const& event, aod::McParticles const& /*mcParticles*/)
  {
    VarManager::ResetValues(0, VarManager::kNVars); // reset variables before filling
    VarManager::FillEvent<gkEventFillMapWithMults>(event);
    VarManager::FillTrack<gkTrackFillMapWithCov>(hadron);
    VarManager::FillTrackCollision<gkTrackFillMapWithCov>(hadron, event);

    // Check that hadron is not one of the dilepton legs
    if (hadron.globalIndex() == track1.globalIndex() || hadron.globalIndex() == track2.globalIndex()) {
      return;
    }

    if (!track1.has_mcParticle() || !track2.has_mcParticle() || !hadron.has_mcParticle()) {
      return;
    }
    auto hadronMC = hadron.mcParticle();
    auto lepton1MC = track1.mcParticle();
    auto lepton2MC = track2.mcParticle();
    uint32_t mcDecision = 0;
    int isig = 0;
    for (auto sig = fRecMCTripleSignals.begin(); sig != fRecMCTripleSignals.end(); sig++, isig++) {
      if ((*sig)->CheckSignal(true, lepton1MC, lepton2MC, hadronMC)) {
        mcDecision |= (static_cast<uint32_t>(1) << isig);
      }
    }
    auto motherParticle = lepton1MC.template mothers_first_as<McParticles>();
    // Fill dilepton-hadron variables
    std::vector<float> fTransRange = fConfigDileptonHadronOptions.fConfigTransRange;
    VarManager::FillEnergyCorrelatorTriple(track1, track2, hadron, VarManager::fgValues, fTransRange[0], fTransRange[1], fConfigDileptonHadronOptions.fConfigApplyMassEC.value);
    if (fConfigDileptonHadronOptions.fConfigUsePionMass.value) {
      VarManager::FillEnergyCorrelatorsUnfoldingTriple<VarManager::kJpsiPionMass>(track1, track2, hadron, motherParticle, hadronMC, VarManager::fgValues, fConfigDileptonHadronOptions.fConfigApplyMassEC.value);
    } else {
      VarManager::FillEnergyCorrelatorsUnfoldingTriple<VarManager::kJpsiHadronMass>(track1, track2, hadron, motherParticle, hadronMC, VarManager::fgValues, fConfigDileptonHadronOptions.fConfigApplyMassEC.value);
    }

    int iHadronCut = 0;
    for (auto hCut = fHadronCuts.begin(); hCut != fHadronCuts.end(); hCut++, iHadronCut++) {
      if (!(*hCut)->IsSelected(VarManager::fgValues)) {
        continue;
      }
      // Fill the corresponding histogram
      if (!MixedEvent) {
        fHistMan->FillHistClass(
          Form("DileptonTrack_%s_%s", fTrackCutNames[iEleCut].Data(), fHadronCutNames[iHadronCut].Data()),
          VarManager::fgValues);
      }
      for (uint32_t isig = 0; isig < fRecMCTripleSignals.size(); isig++) {
        if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
          if (!MixedEvent) {
            fHistMan->FillHistClass(Form("DileptonTrackMCMatched_%s_%s_%s", fTrackCutNames[iEleCut].Data(), fHadronCutNames[iHadronCut].Data(), fRecMCTripleSignals[isig]->GetName()), VarManager::fgValues);
          }
          if (MixedEvent) {
            fHistMan->FillHistClass(Form("DileptonTrackMCMatchedME_%s_%s_%s", fTrackCutNames[iEleCut].Data(), fHadronCutNames[iHadronCut].Data(), fRecMCTripleSignals[isig]->GetName()), VarManager::fgValues);
          }
        } // end loop over MC signals
      }
    }
  }

  // Template function to run same event pairing (barrel-barrel, muon-muon, barrel-muon)
  template <uint32_t TTrackFillMap, int TPairType, typename TTrack1, typename TTrack2>
  void runSameEventPairing(TTrack1 const& t1, TTrack2 const& t2, int iEleCut, aod::McParticles const& /*mcParticles*/)
  {
    std::map<int, std::vector<TString>> histNames = fTrackHistNames;
    std::map<int, std::vector<TString>> histNamesMC = fBarrelHistNamesMCmatched;
    int sign1 = t1.sign();
    int sign2 = t2.sign();
    uint32_t mcDecision = static_cast<uint32_t>(0);
    // run MC matching for this pair
    int isig = 0;
    mcDecision = 0;
    for (auto sig = fRecMCSignals.begin(); sig != fRecMCSignals.end(); sig++, isig++) {
      if (t1.has_mcParticle() && t2.has_mcParticle()) {
        if ((*sig)->CheckSignal(true, t1.mcParticle(), t2.mcParticle())) {
          mcDecision |= (static_cast<uint32_t>(1) << isig);
        }
      }
    } // end loop over MC signals

    VarManager::FillPair<TPairType, TTrackFillMap>(t1, t2);

    if (sign1 * sign2 < 0) {                                                       // +- pairs
      fHistMan->FillHistClass(histNames[iEleCut][0].Data(), VarManager::fgValues); // reconstructed, unmatched
      for (size_t isig = 0; isig < fRecMCSignals.size(); isig++) {                 // loop over MC signals
        if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
          fHistMan->FillHistClass(histNamesMC[iEleCut * fRecMCSignals.size() + isig][0].Data(), VarManager::fgValues); // matched signal
        }
      }
    } else {
      if (sign1 > 0) { // ++ pairs
        fHistMan->FillHistClass(histNames[iEleCut][1].Data(), VarManager::fgValues);
        for (size_t isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
          if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
            fHistMan->FillHistClass(histNamesMC[iEleCut * fRecMCSignals.size() + isig][1].Data(), VarManager::fgValues);
          }
        }
      } else { // -- pairs
        fHistMan->FillHistClass(histNames[iEleCut][2].Data(), VarManager::fgValues);
        for (unsigned int isig = 0; isig < fRecMCSignals.size(); isig++) { // loop over MC signals
          if (mcDecision & (static_cast<uint32_t>(1) << isig)) {
            fHistMan->FillHistClass(histNamesMC[iEleCut * fRecMCSignals.size() + isig][2].Data(), VarManager::fgValues);
          }
        }
      }
    }
  }

  void processBarrel(MyEvents const& events, aod::TrackAssoc const& assocs, MyBarrelTracksWithCov const& /*tracks*/, soa::Join<aod::McCollisions, aod::McCollsExtra, aod::MultMCExtras> const& mcEvents, aod::McParticles const& mcParticles, BCsWithTimestamps const& bcs)
  {
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillTimeFrame(bcs);
    VarManager::FillTimeFrame(events);
    VarManager::FillTimeFrame(mcEvents);

    if (events.size() == 0)
      return;

    // CCDB initialization
    if (fCurrentRun != bcs.begin().runNumber()) {
      fCurrentRun = bcs.begin().runNumber();
    }

    if (fConfigEventOptions.fConfigEventQA) {
      fHistMan->FillHistClass("TimeFrameStats", VarManager::fgValues);
    }

    for (auto& event : mcEvents) {
      // Reset the fValues array and fill event observables
      VarManager::FillEvent<VarManager::ObjTypes::CollisionMC>(event);
      if (fConfigEventOptions.fConfigEventQA) {
        fHistMan->FillHistClass("EventsMC", VarManager::fgValues);
      }
    }

    // Event loop
    for (auto& event : events) {
      // Fill event variables first
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<gkEventFillMapWithMults>(event);

      if (fConfigEventOptions.fConfigEventQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }

      // Event selection
      if (!fEventCut->IsSelected(VarManager::fgValues))
        continue;

      if (fConfigEventOptions.fConfigEventQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }

      // saveless events for the energy correlator analysis
      std::vector<int> fSavelessevents = fConfigDileptonHadronOptions.fConfigSavelessevents.value;
      if (fSavelessevents[0] > 1 && event.globalIndex() % fSavelessevents[0] == fSavelessevents[1]) {
        continue;
      }

      // Get associated tracks for this event
      auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());

      // Triple loop: track1 (electron) x track2 (electron) x hadron
      for (auto& a1 : groupedAssocs) {
        auto t1 = a1.template track_as<MyBarrelTracksWithCov>();

        uint32_t filter1 = 0;
        // Fill track variables
        VarManager::FillTrack<gkTrackFillMapWithCov>(t1);
        VarManager::FillTrackCollision<gkTrackFillMapWithCov>(t1, event);
        if (t1.has_mcParticle()) {
          VarManager::FillTrackMC(mcParticles, t1.mcParticle());
        }

        if (fConfigTrackOptions.fConfigTrackQA) {
          fHistMan->FillHistClass("AssocsBarrel_BeforeCuts", VarManager::fgValues);
        }

        // Apply electron cuts and fill histograms
        int iCut1 = 0;
        for (auto cut1 = fTrackCuts.begin(); cut1 != fTrackCuts.end(); cut1++, iCut1++) {
          if ((*cut1)->IsSelected(VarManager::fgValues)) {
            filter1 |= (static_cast<uint32_t>(1) << iCut1);
            if (fConfigTrackOptions.fConfigTrackQA) {
              fHistMan->FillHistClass(fHistNamesReco[iCut1], VarManager::fgValues);
            }
          }
        }

        // Check opposite charge with t2
        for (auto& a2 : groupedAssocs) {
          auto t2 = a2.template track_as<MyBarrelTracksWithCov>();

          // Avoid double counting: use track globalIndex
          if (t2.globalIndex() <= t1.globalIndex()) {
            continue;
          }

          // Fill track variables for t2 (only once per t2)
          VarManager::FillTrack<gkTrackFillMapWithCov>(t2);
          VarManager::FillTrackCollision<gkTrackFillMapWithCov>(t2, event);

          // Compute filter2: which cuts t2 passes
          uint32_t filter2 = 0;
          int iCut2 = 0;
          for (auto cut2 = fTrackCuts.begin(); cut2 != fTrackCuts.end(); cut2++, iCut2++) {
            if ((*cut2)->IsSelected(VarManager::fgValues)) {
              filter2 |= (static_cast<uint32_t>(1) << iCut2);
            }
          }

          // Both tracks must pass at least one common cut
          uint32_t twoTrackFilter = filter1 & filter2;
          if (!twoTrackFilter) {
            continue;
          }

          // Fill pair histograms for all cuts that both tracks pass
          for (size_t iCut = 0; iCut < fTrackCuts.size(); iCut++) {
            if (twoTrackFilter & (static_cast<uint32_t>(1) << iCut)) {
              runSameEventPairing<gkTrackFillMapWithCov, VarManager::kDecayToEE>(t1, t2, iCut, mcParticles);
            }
          }

          float mass = VarManager::fgValues[VarManager::kMass];
          float pt = VarManager::fgValues[VarManager::kPt];
          float rap = VarManager::fgValues[VarManager::kRap];

          // Apply J/psi cuts
          if (mass < fConfigPairOptions.fConfigJpsiMassMin.value || mass > fConfigPairOptions.fConfigJpsiMassMax.value ||
              pt < fConfigPairOptions.fConfigJpsiPtMin.value || pt > fConfigPairOptions.fConfigJpsiPtMax.value ||
              std::abs(rap) > fConfigPairOptions.fConfigJpsiRapMax.value) {
            continue;
          }

          if (t1.sign() * t2.sign() >= 0) {
            continue; // Must be opposite charge
          }

          // correlate J/psi with hadrons
          for (auto& aHadron : groupedAssocs) {
            auto hadron = aHadron.template track_as<MyBarrelTracksWithCov>();
            // Process dilepton-hadron correlation for each common cut
            for (size_t iCut = 0; iCut < fTrackCuts.size(); iCut++) {
              if (twoTrackFilter & (static_cast<uint32_t>(1) << iCut)) {
                runDileptonHadron<false, gkTrackFillMapWithCov>(t1, t2, iCut, hadron, event, mcParticles);
              }
            }
          } // end hadron loop
        } // end track2 loop
      } // end track1 loop
    } // end event loop
  }

  Filter eventFilter = nabs(aod::collision::posZ) < fConfigEventOptions.fConfigEventfilterVtz && aod::evsel::sel8 == true;
  void processBarrelMixedEvent(soa::Filtered<MyEvents>& events, aod::TrackAssoc const& assocs, MyBarrelTracksWithCov const& /*tracks*/, aod::McCollisions const& /*mcCollisions*/, aod::McParticles const& mcParticles, BCsWithTimestamps const& bcs)
  {
    if (events.size() == 0) {
      return;
    }

    // CCDB initialization
    if (fCurrentRun != bcs.begin().runNumber()) {
      fCurrentRun = bcs.begin().runNumber();
    }

    fSelMap.clear();
    // Event loop
    for (auto& event : events) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMapWithMults>(event);
      if (event.has_mcCollision()) {
        VarManager::FillEvent<VarManager::ObjTypes::CollisionMC>(event.mcCollision());
      }
      bool decision = false;
      if (fEventCut->IsSelected(VarManager::fgValues)) {
        decision = true;
      }
      fSelMap[event.globalIndex()] = decision;
    }

    for (auto& [event1, event2] : selfCombinations(*fMixingBinning, fConfigEventOptions.fConfigMixingDepth.value, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      if (!fSelMap[event1.globalIndex()] || !fSelMap[event2.globalIndex()]) {
        continue;
      }

      // save less events if configured
      std::vector<int> fSavelessevents = fConfigDileptonHadronOptions.fConfigSavelessevents.value;
      if (fSavelessevents[0] > 1 && event1.globalIndex() % fSavelessevents[0] == fSavelessevents[1]) {
        continue;
      }

      // Get associated tracks for this event
      auto groupedAssocs1 = assocs.sliceBy(preslice, event1.globalIndex());
      if (groupedAssocs1.size() < 2) {
        continue; // Need at least 2 tracks for pairing
      }
      auto groupedAssocs2 = assocs.sliceBy(preslice, event2.globalIndex());

      // Triple loop: track1 (electron) x track2 (electron) x hadron
      for (auto& a1 : groupedAssocs1) {
        auto t1 = a1.template track_as<MyBarrelTracksWithCov>();

        uint32_t filter1 = 0;
        // Fill track variables
        VarManager::FillTrack<gkTrackFillMapWithCov>(t1);
        VarManager::FillTrackCollision<gkTrackFillMapWithCov>(t1, event1);

        // Apply electron cuts and fill histograms
        int iCut1 = 0;
        for (auto cut1 = fTrackCuts.begin(); cut1 != fTrackCuts.end(); cut1++, iCut1++) {
          if ((*cut1)->IsSelected(VarManager::fgValues)) {
            filter1 |= (static_cast<uint32_t>(1) << iCut1);
          }
        }

        // Check opposite charge with t2
        for (auto& a2 : groupedAssocs1) {
          auto t2 = a2.template track_as<MyBarrelTracksWithCov>();

          // Avoid double counting: use track globalIndex
          if (t2.globalIndex() <= t1.globalIndex())
            continue;

          // Fill track variables for t2 (only once per t2)
          VarManager::FillTrack<gkTrackFillMapWithCov>(t2);
          VarManager::FillTrackCollision<gkTrackFillMapWithCov>(t2, event1);

          // Compute filter2: which cuts t2 passes
          uint32_t filter2 = 0;
          int iCut2 = 0;
          for (auto cut2 = fTrackCuts.begin(); cut2 != fTrackCuts.end(); cut2++, iCut2++) {
            if ((*cut2)->IsSelected(VarManager::fgValues)) {
              filter2 |= (static_cast<uint32_t>(1) << iCut2);
            }
          }

          // Both tracks must pass at least one common cut
          uint32_t twoTrackFilter = filter1 & filter2;
          if (!twoTrackFilter) {
            continue;
          }
          // Fill pair variables for cut
          VarManager::FillPair<VarManager::kDecayToEE, gkTrackFillMapWithCov>(t1, t2, VarManager::fgValues);
          float mass = VarManager::fgValues[VarManager::kMass];
          float pt = VarManager::fgValues[VarManager::kPt];
          float rap = VarManager::fgValues[VarManager::kRap];
          // Apply J/psi cuts
          if (mass < fConfigPairOptions.fConfigJpsiMassMin.value || mass > fConfigPairOptions.fConfigJpsiMassMax.value ||
              pt < fConfigPairOptions.fConfigJpsiPtMin.value || pt > fConfigPairOptions.fConfigJpsiPtMax.value ||
              std::abs(rap) > fConfigPairOptions.fConfigJpsiRapMax.value) {
            continue;
          }
          if (t1.sign() * t2.sign() >= 0) {
            continue; // Must be opposite charge
          }
          // correlate J/psi with hadrons from different events
          for (auto& aHadron : groupedAssocs2) {
            auto hadron = aHadron.template track_as<MyBarrelTracksWithCov>();
            // Process dilepton-hadron correlation for each common cut
            for (size_t iCut = 0; iCut < fTrackCuts.size(); iCut++) {
              if (twoTrackFilter & (static_cast<uint32_t>(1) << iCut)) {
                runDileptonHadron<true, gkTrackFillMapWithCov>(t1, t2, iCut, hadron, event2, mcParticles);
              }
            }
          } // end hadron loop
        } // end track2 loop
      } // end track1 loop
    } // end event loop
  }

  PresliceUnsorted<aod::McParticles> perReducedMcEvent = aod::mcparticle::mcCollisionId;
  template <bool MixedEvent, bool PionMass, int THadronMassType, typename TEvent>
  void runEnergyCorrelators(TEvent const& event1, TEvent const& event2, McParticles const& mcTracks)
  {
    auto groupedMCTracks1 = mcTracks.sliceBy(perReducedMcEvent, event1.mcCollisionId());
    auto groupedMCTracks2 = mcTracks.sliceBy(perReducedMcEvent, event2.mcCollisionId());
    groupedMCTracks1.bindInternalIndicesTo(&mcTracks);
    groupedMCTracks2.bindInternalIndicesTo(&mcTracks);
    for (auto& t1 : groupedMCTracks1) {
      auto t1_raw = mcTracks.rawIteratorAt(t1.globalIndex());
      for (auto& sig : fGenMCSignals) {
        if (sig->CheckSignal(true, t1_raw)) {
          if (t1.mcCollisionId() != event1.mcCollisionId()) { // check that the mc track belongs to the same mc collision as the reconstructed event
            continue;
          }
          VarManager::FillTrackMC(mcTracks, t1_raw);
          if (!MixedEvent && !PionMass) {
            fHistMan->FillHistClass(Form("MCTruthGenSel_%s", sig->GetName()), VarManager::fgValues);
          }
        }
      }
      // apply kinematic cuts for signal
      if ((t1_raw.pt() < fConfigPairOptions.fConfigJpsiPtMin || t1_raw.pt() > fConfigPairOptions.fConfigJpsiPtMax))
        continue;
      if (std::abs(t1_raw.y()) > fConfigPairOptions.fConfigJpsiRapMax)
        continue;
      // for the energy correlators
      for (auto& t2 : groupedMCTracks2) {
        auto t2_raw = groupedMCTracks2.rawIteratorAt(t2.globalIndex());
        if (t2.mcCollisionId() != event2.mcCollisionId()) { // check that the mc track belongs to the same mc collision as the reconstructed event
          continue;
        }
        if (!t2_raw.isPhysicalPrimary()) {
          continue;
        }
        if (t2_raw.has_mothers()) {
          auto mother_raw = t2_raw.template mothers_first_as<McParticles>();
          if (mother_raw.globalIndex() == t1_raw.globalIndex()) {
            continue;
          }
        }
        if (fConfigDileptonHadronOptions.fConfigContainlepton && std::abs(t2_raw.pdgCode()) != PDG_t::kPiPlus && std::abs(t2_raw.pdgCode()) != PDG_t::kKPlus && std::abs(t2_raw.pdgCode()) != PDG_t::kProton && std::abs(t2_raw.pdgCode()) != PDG_t::kElectron && std::abs(t2_raw.pdgCode()) != PDG_t::kMuonMinus) {
          continue;
        }
        if (!fConfigDileptonHadronOptions.fConfigContainlepton && std::abs(t2_raw.pdgCode()) != PDG_t::kPiPlus && std::abs(t2_raw.pdgCode()) != PDG_t::kKPlus && std::abs(t2_raw.pdgCode()) != PDG_t::kProton) {
          continue;
        }
        if (t2_raw.pt() < fConfigDileptonHadronOptions.fConfigMCGenHadronPtMin.value || std::abs(t2_raw.eta()) > fConfigDileptonHadronOptions.fConfigMCGenHadronEtaAbs.value) {
          continue;
        }
        std::vector<float> fTransRange = fConfigDileptonHadronOptions.fConfigTransRange;
        VarManager::FillEnergyCorrelatorsMC<THadronMassType>(t1_raw, t2_raw, VarManager::fgValues, fTransRange[0], fTransRange[1]);
        for (auto& sig : fGenMCSignals) {
          if (sig->CheckSignal(true, t1_raw)) {
            if (!MixedEvent && !PionMass) {
              fHistMan->FillHistClass(Form("MCTruthEenergyCorrelators_%s", sig->GetName()), VarManager::fgValues);
            }
            if (MixedEvent && !PionMass) {
              fHistMan->FillHistClass(Form("MCTruthEenergyCorrelatorsME_%s", sig->GetName()), VarManager::fgValues);
            }
            if (!MixedEvent && PionMass) {
              fHistMan->FillHistClass(Form("MCTruthEenergyCorrelators_Pion_%s", sig->GetName()), VarManager::fgValues);
            }
            if (MixedEvent && PionMass) {
              fHistMan->FillHistClass(Form("MCTruthEenergyCorrelatorsME_Pion_%s", sig->GetName()), VarManager::fgValues);
            }
          }
        }
      }
    }
  }

  void processMCGenEnergyCorrelators(soa::Filtered<MyEvents>& events,
                                     McCollisions const& /*mcEvents*/, McParticles const& mcTracks)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    for (auto& event : events) {
      // Fill event variables first
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMapWithMults>(event);
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }
      if (!event.has_mcCollision()) {
        continue;
      }
      // saveless events for the energy correlator analysis
      std::vector<int> fSavelessevents = fConfigDileptonHadronOptions.fConfigSavelessevents.value;
      if (fSavelessevents[0] > 1 && event.globalIndex() % fSavelessevents[0] == fSavelessevents[1]) {
        continue;
      }
      runEnergyCorrelators<false, false, VarManager::kJpsiHadronMass>(event, event, mcTracks);
    }
  }

  void processMCGenEnergyCorrelatorsME(soa::Filtered<MyEvents>& events,
                                       McCollisions const& /*mcEvents*/, McParticles const& mcTracks)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    // loop over two event comibnations
    for (auto& [event1, event2] : selfCombinations(*fMixingBinning, fConfigEventOptions.fConfigMixingDepth.value, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMapWithMults>(event1);
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }
      VarManager::FillEvent<gkEventFillMapWithMults>(event2);
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }
      if (!event1.has_mcCollision() || !event2.has_mcCollision()) {
        continue;
      }
      // saveless events for the energy correlator analysis
      std::vector<int> fSavelessevents = fConfigDileptonHadronOptions.fConfigSavelessevents.value;
      if (fSavelessevents[0] > 1 && event1.globalIndex() % fSavelessevents[0] == fSavelessevents[1]) {
        continue;
      }
      runEnergyCorrelators<true, false, VarManager::kJpsiHadronMass>(event1, event2, mcTracks);
    }
  }

  void processMCGenEnergyCorrelatorsPion(soa::Filtered<MyEvents>& events,
                                         McCollisions const& /*mcEvents*/, McParticles const& mcTracks)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    for (auto& event : events) {
      // Fill event variables first
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMapWithMults>(event);
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }
      if (!event.has_mcCollision()) {
        continue;
      }
      // saveless events for the energy correlator analysis
      std::vector<int> fSavelessevents = fConfigDileptonHadronOptions.fConfigSavelessevents.value;
      if (fSavelessevents[0] > 1 && event.globalIndex() % fSavelessevents[0] == fSavelessevents[1]) {
        continue;
      }
      runEnergyCorrelators<false, true, VarManager::kJpsiPionMass>(event, event, mcTracks);
    }
  }

  void processMCGenEnergyCorrelatorsPionME(soa::Filtered<MyEvents>& events,
                                           McCollisions const& /*mcEvents*/, McParticles const& mcTracks)
  {
    if (events.size() == 0) {
      LOG(warning) << "No events in this TF, going to the next one ...";
      return;
    }
    // loop over two event comibnations
    for (auto& [event1, event2] : selfCombinations(*fMixingBinning, fConfigEventOptions.fConfigMixingDepth.value, -1, events, events)) {
      VarManager::ResetValues(0, VarManager::kNVars);
      VarManager::FillEvent<gkEventFillMapWithMults>(event1);
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }
      VarManager::FillEvent<gkEventFillMapWithMults>(event2);
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }
      if (!event1.has_mcCollision() || !event2.has_mcCollision()) {
        continue;
      }
      // saveless events for the energy correlator analysis
      std::vector<int> fSavelessevents = fConfigDileptonHadronOptions.fConfigSavelessevents.value;
      if (fSavelessevents[0] > 1 && event1.globalIndex() % fSavelessevents[0] == fSavelessevents[1]) {
        continue;
      }
      runEnergyCorrelators<true, true, VarManager::kJpsiPionMass>(event1, event2, mcTracks);
    }
  }

  void processDummy(aod::Collisions const&)
  {
    // Do nothing
  }

  PROCESS_SWITCH(AnalysisEnergyCorrelator, processBarrel, "Process barrel analysis", false);
  PROCESS_SWITCH(AnalysisEnergyCorrelator, processBarrelMixedEvent, "Run barrel dilepton-hadron mixed event pairing", false);
  PROCESS_SWITCH(AnalysisEnergyCorrelator, processMCGenEnergyCorrelators, "Loop over MC particle stack and fill generator level histograms(energy correlators)", false);
  PROCESS_SWITCH(AnalysisEnergyCorrelator, processMCGenEnergyCorrelatorsPion, "Loop over MC particle stack and fill generator level histograms(energy correlators)", false);
  PROCESS_SWITCH(AnalysisEnergyCorrelator, processMCGenEnergyCorrelatorsME, "Loop over MC particle stack and fill generator level histograms(energy correlators)", false);
  PROCESS_SWITCH(AnalysisEnergyCorrelator, processMCGenEnergyCorrelatorsPionME, "Loop over MC particle stack and fill generator level histograms(energy correlators)", false);
  PROCESS_SWITCH(AnalysisEnergyCorrelator, processDummy, "Dummy process function", true);
};

// Histogram definitions
void DefineHistograms(HistogramManager* histMan, TString histClasses, const char* histGroups)
{
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());
    TString histName = histGroups;
    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("TimeFrameStats")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "timeframe");
    }
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", histName);
    }
    if ((classStr.Contains("Track") || classStr.Contains("Assoc")) && !classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
      }
    }
    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }
    if (classStr.Contains("DileptonTrack")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-track", histName);
    }
    if (classStr.Contains("MCTruthEenergyCorrelators")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "energy-correlator-gen");
    }
    if (classStr.Contains("MCTruthGenSel")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth_track");
    }
  } // end loop over histogram classes
}

// Workflow definition
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisEnergyCorrelator>(cfgc)};
}
