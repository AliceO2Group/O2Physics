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
//
#include <iostream>
#include <vector>
#include <memory>
#include <cstring>
#include <TH1F.h>
#include <TH2I.h>
#include <THashList.h>
#include <TString.h>
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/filterTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "CommonConstants/LHCConstants.h"
#include "Common/Core/CollisionAssociation.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Some definitions
namespace
{
enum DQTriggers {
  kSingleE = 0,
  kLMeeIMR,
  kLMeeHMR,
  kDiElectron,
  kSingleMuLow,
  kSingleMuHigh,
  kDiMuon,
  kNTriggersDQ
};
} // namespace
namespace o2::aod
{
namespace dqppfilter
{
DECLARE_SOA_COLUMN(IsDQEventSelected, isDQEventSelected, int);
DECLARE_SOA_COLUMN(IsDQBarrelSelected, isDQBarrelSelected, uint32_t);
DECLARE_SOA_COLUMN(IsDQMuonSelected, isDQMuonSelected, uint32_t);
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //! Collision index
DECLARE_SOA_INDEX_COLUMN(Track, track);         //! Track index
DECLARE_SOA_INDEX_COLUMN(FwdTrack, fwdtrack);   //! FwdTrack index
} // namespace dqppfilter

DECLARE_SOA_TABLE(DQEventCuts, "AOD", "DQEVENTCUTS", dqppfilter::IsDQEventSelected);
DECLARE_SOA_TABLE(DQBarrelTrackCuts, "AOD", "DQBARRELCUTS", dqppfilter::IsDQBarrelSelected);
DECLARE_SOA_TABLE(DQMuonsCuts, "AOD", "DQMUONCUTS", dqppfilter::IsDQMuonSelected);
} // namespace o2::aod

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsSelected = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventCuts>;
// TODO: subscribe to the bare needed minimum, in particular for the CEFP task
// TODO: test working with TrackIU
// TODO: remove the TOF pid if not needed (requires changes in VarManager to separate TrackPID into TPC and TOF)
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                 aod::pidTPCFullEl, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr>;
using MyBarrelTracksSelected = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                         aod::pidTPCFullEl, aod::pidTPCFullPi,
                                         aod::pidTPCFullKa, aod::pidTPCFullPr,
                                         aod::pidTOFFullEl, aod::pidTOFFullPi,
                                         aod::pidTOFFullKa, aod::pidTOFFullPr,
                                         aod::DQBarrelTrackCuts>;

using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA>;
using MyMuonsAssocSelected = soa::Join<FwdTrackAssoc, aod::DQMuonsCuts>; // As the kinelatic values must be re-computed for the muons tracks everytime it is associated to a collision, the selection is done not on the muon, but on the muon-collision association

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct DQEventSelectionTask {
  Produces<aod::DQEventCuts> eventSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan = nullptr;
  AnalysisCompositeCut* fEventCut;

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Comma separated list of event cuts; multiple cuts are applied with a logical AND"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  // TODO: configure the histogram classes to be filled by QA

  void init(o2::framework::InitContext&)
  {
    // Construct the composite cut out of the user provided event selection cut(s)
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut = new AnalysisCompositeCut(true);
    if (!eventCutStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(eventCutStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fEventCut->AddCut(dqcuts::GetAnalysisCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSelection(TEvent const& collision, aod::BCs const& bcs)
  {
    // Reset the Values array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);

    VarManager::FillEvent<gkEventFillMap>(collision);
    if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
    }
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }
      eventSel(1);
    } else {
      eventSel(0);
    }
  }

  void processEventSelection(MyEvents::iterator const& collision, aod::BCs const& bcs)
  {
    runEventSelection<gkEventFillMap>(collision, bcs);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQEventSelectionTask, processEventSelection, "Run event selection", false);
  PROCESS_SWITCH(DQEventSelectionTask, processDummy, "Dummy function", false);
};

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
    if (fCurrentRun != bc.runNumber()) {
      fCurrentRun = bc.runNumber();
      if (fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, bc.timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
      }
    }

    uint32_t filterMap = uint32_t(0);
    trackSel.reserve(tracksBarrel.size());

    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
    for (auto& track : tracksBarrel) {
      filterMap = uint32_t(0);

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

  PROCESS_SWITCH(DQBarrelTrackSelection, processSelection, "Run barrel track selection", false);
  PROCESS_SWITCH(DQBarrelTrackSelection, processDummy, "Dummy function", false);
};

struct DQMuonsSelection {
  Produces<aod::DQMuonsCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigCuts{"cfgMuonsCuts", "muonQualityCuts", "Comma separated list of ADDITIONAL muon track cuts"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  Configurable<bool> fPropMuon{"cfgPropMuon", false, "Propgate muon tracks through absorber"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::parameters::GRPMagField* grpmag = nullptr; // for run 3, we access GRPMagField from GLO/Config/GRPMagField

  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  // TODO: configure the histogram classes to be filled by QA

  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<TString> fCutHistNames;

  void init(o2::framework::InitContext&)
  {
    fCCDB->setURL(fConfigCcdbUrl);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    if (fPropMuon) {
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        fCCDB->get<TGeoManager>(geoPath);
      }
    }

    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      TString cutNames = "Muon_BeforeCuts;";
      for (unsigned int i = 0; i < fTrackCuts.size(); i++) {
        cutNames += Form("Muon_%s;", fTrackCuts[i].GetName());
        fCutHistNames.push_back(Form("Muon_%s", fTrackCuts[i].GetName()));
      }

      DefineHistograms(fHistMan, cutNames.Data());     // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TMuonFillMap, typename TEvent, typename TMuons, typename AssocMuons>
  void runMuonSelection(TEvent const& collision, aod::BCsWithTimestamps const& bcs, TMuons const& muons, AssocMuons const& muonAssocs)
  {
    auto bc = bcs.begin(); // check just the first bc to get the run number
    if (fCurrentRun != bc.runNumber()) {
      fCurrentRun = bc.runNumber();
      if (fPropMuon) {
        VarManager::SetupMuonMagField();
      }
      grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
      if (grpmag != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
      }
    }

    uint32_t filterMap = uint32_t(0);
    trackSel.reserve(muons.size());

    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);

    for (auto& muonAssoc : muonAssocs) {
      filterMap = uint32_t(0);
      auto muon = muonAssoc.template fwdtrack_as<TMuons>();
      VarManager::FillTrack<TMuonFillMap>(muon);
      if (fPropMuon) {
        VarManager::FillPropagateMuon<TMuonFillMap>(muon, collision);
      }
      if (fConfigQA) {
        fHistMan->FillHistClass("Muon_BeforeCuts", VarManager::fgValues);
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
    } // end loop over muons
  }

  void processSelection(Collisions const& collisions, BCsWithTimestamps const& bcstimestamps, MyMuons const& muons, aod::FwdTrackAssoc const& muonAssocs)
  {
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = muonAssocs.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      runMuonSelection<gkMuonFillMap>(collision, bcstimestamps, muons, muonIdsThisCollision);
    }
  }
  void processDummy(MyMuons&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQMuonsSelection, processSelection, "Run muon selection", false);
  PROCESS_SWITCH(DQMuonsSelection, processDummy, "Dummy function", false);
};

/*
struct DQTrackToCollisionAssociation {

  Produces<TrackAssoc> association;
  Produces<TrackCompColls> reverseIndices;
  Produces<FwdTrackAssoc> fwdassociation;
  Produces<FwdTrkCompColls> fwdreverseIndices;

  // NOTE: the options for the collision associator are common for both the barrel and muon
  //       We should add separate ones if needed
  Configurable<float> nSigmaForTimeCompat{"nSigmaForTimeCompat", 4.f, "number of sigmas for time compatibility"};
  Configurable<float> timeMargin{"timeMargin", 0.f, "time margin in ns added to uncertainty because of uncalibrated TPC"};
  Configurable<bool> usePVAssociation{"usePVAssociation", true, "if the track is a PV contributor, use the collision time for it"};
  Configurable<bool> includeUnassigned{"includeUnassigned", false, "consider also tracks which are not assigned to any collision"};
  Configurable<bool> fillTableOfCollIdsPerTrack{"fillTableOfCollIdsPerTrack", false, "fill additional table with vector of collision ids per track"};
  Configurable<int> bcWindowForOneSigma{"bcWindowForOneSigma", 60, "BC window to be multiplied by the number of sigmas to define maximum window to be considered"};

  CollisionAssociation<true> collisionAssociatorBarrel;
  CollisionAssociation<false> collisionAssociatorMuon;

  Filter filterBarrelTrackSelected = aod::dqppfilter::isDQBarrelSelected > uint32_t(0);
  Filter filterMuonTrackSelected = aod::dqppfilter::isDQMuonSelected > uint32_t(0);

  void init(o2::framework::InitContext const&)
  {
    // set options in track-to-collision association
    collisionAssociatorBarrel.setNumSigmaForTimeCompat(nSigmaForTimeCompat);
    collisionAssociatorBarrel.setTimeMargin(timeMargin);
    collisionAssociatorBarrel.setTrackSelectionOptionForStdAssoc(track_association::TrackSelection::None);
    collisionAssociatorBarrel.setUsePvAssociation(usePVAssociation);
    collisionAssociatorBarrel.setIncludeUnassigned(includeUnassigned);
    collisionAssociatorBarrel.setFillTableOfCollIdsPerTrack(fillTableOfCollIdsPerTrack);
    collisionAssociatorBarrel.setBcWindow(bcWindowForOneSigma);
    // set options in muon-to-collision association
    collisionAssociatorMuon.setNumSigmaForTimeCompat(nSigmaForTimeCompat);
    collisionAssociatorMuon.setTimeMargin(timeMargin);
    collisionAssociatorMuon.setTrackSelectionOptionForStdAssoc(track_association::TrackSelection::None);
    collisionAssociatorMuon.setUsePvAssociation(false);
    collisionAssociatorMuon.setIncludeUnassigned(includeUnassigned);
    collisionAssociatorMuon.setFillTableOfCollIdsPerTrack(fillTableOfCollIdsPerTrack);
    collisionAssociatorMuon.setBcWindow(bcWindowForOneSigma);
  }

  void processAssocWithTime(Collisions const& collisions,
                            MyBarrelTracksSelected const& tracksUnfiltered, soa::Filtered<MyBarrelTracksSelected> const& tracks,
                            FwdTracks const& muons,
                            AmbiguousTracks const& ambiguousTracks, AmbiguousFwdTracks const& ambiguousFwdTracks, BCs const& bcs)
  {
    collisionAssociatorBarrel.runAssocWithTime(collisions, tracksUnfiltered, tracks, ambiguousTracks, bcs, association, reverseIndices);
    collisionAssociatorMuon.runAssocWithTime(collisions, muons, muons, ambiguousFwdTracks, bcs, fwdassociation, fwdreverseIndices);
  };
  void processDummy(Collisions&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQTrackToCollisionAssociation, processAssocWithTime, "Produce track-to-collision associations based on time", false);
  PROCESS_SWITCH(DQTrackToCollisionAssociation, processDummy, "Dummy function", false);
};
*/

struct DQFilterPPTask {
  Produces<aod::DQEventFilter> eventFilter;
  Produces<aod::DqFilters> dqtable;
  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TH1I> fStats{"Statistics"};
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigBarrelSelections{"cfgBarrelSels", "jpsiPID1:pairMassLow:1", "<track-cut>:[<pair-cut>]:<n>,[<track-cut>:[<pair-cut>]:<n>],..."};
  Configurable<std::string> fConfigMuonSelections{"cfgMuonSels", "muonQualityCuts:pairNoCut:1", "<muon-cut>:[<pair-cut>]:<n>"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigFilterLsBarrelTracksPairs{"cfgWithBarrelLS", "false", "Comma separated list of booleans for each trigger, If true, also select like sign (--/++) barrel track pairs"};
  Configurable<std::string> fConfigFilterLsMuonsPairs{"cfgWithMuonLS", "false", "Comma separated list of booleans for each trigger, If true, also select like sign (--/++) muon pairs"};

  Configurable<bool> fPropMuon{"cfgPropMuon", false, "Propgate muon tracks through absorber"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::parameters::GRPMagField* grpmag = nullptr; // for run 3, we access GRPMagField from GLO/Config/GRPMagField

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.

  int fNBarrelCuts;                                    // number of barrel selections
  int fNMuonCuts;                                      // number of muon selections
  std::vector<bool> fBarrelRunPairing;                 // bit map on whether the selections require pairing (barrel)
  std::vector<bool> fMuonRunPairing;                   // bit map on whether the selections require pairing (muon)
  std::vector<int> fBarrelNreqObjs;                    // minimal number of tracks/pairs required (barrel)
  std::vector<int> fMuonNreqObjs;                      // minimal number of tracks/pairs required (muon)
  std::map<int, AnalysisCompositeCut> fBarrelPairCuts; // map of barrel pair cuts
  std::map<int, AnalysisCompositeCut> fMuonPairCuts;   // map of muon pair cuts
  std::map<int, TString> fBarrelPairHistNames;         // map with names of the barrel pairing histogram directories
  std::map<int, TString> fMuonPairHistNames;           // map with names of the muon pairing histogram directories

  std::map<uint64_t, uint64_t> fFiltersMap;           // map of filters for events that passed at least one filter
  std::map<uint64_t, std::vector<bool>> fCEFPfilters; // map of CEFP filters for events that passed at least one filter

  void DefineCuts()
  {
    TString barrelSelsStr = fConfigBarrelSelections.value;
    std::unique_ptr<TObjArray> objArray(barrelSelsStr.Tokenize(","));
    fNBarrelCuts = objArray->GetEntries();
    if (fNBarrelCuts) {
      for (int icut = 0; icut < fNBarrelCuts; ++icut) {
        TString selStr = objArray->At(icut)->GetName();
        std::unique_ptr<TObjArray> sel(selStr.Tokenize(":"));
        if (sel->GetEntries() < 2 || sel->GetEntries() > 3) {
          continue;
        }
        // if "sel" contains 3 entries, it means the user provided both the track and pair cuts
        if (sel->GetEntries() == 3) {
          fBarrelPairCuts[icut] = (*dqcuts::GetCompositeCut(sel->At(1)->GetName()));
          fBarrelRunPairing.push_back(true);
          fBarrelNreqObjs.push_back(std::atoi(sel->At(2)->GetName()));
          fBarrelPairHistNames[icut] = Form("PairsBarrelSEPM_%s_%s", sel->At(0)->GetName(), sel->At(1)->GetName());
        } else {
          fBarrelNreqObjs.push_back(std::atoi(sel->At(1)->GetName()));
          fBarrelRunPairing.push_back(false);
        }
      }
    }
    TString muonSelsStr = fConfigMuonSelections.value;
    std::unique_ptr<TObjArray> objArray2(muonSelsStr.Tokenize(","));
    fNMuonCuts = objArray2->GetEntries();
    if (fNMuonCuts) {
      for (int icut = 0; icut < fNMuonCuts; ++icut) {
        TString selStr = objArray2->At(icut)->GetName();
        std::unique_ptr<TObjArray> sel(selStr.Tokenize(":"));
        if (sel->GetEntries() < 2 || sel->GetEntries() > 3) {
          continue;
        }
        // if "sel" contains 3 entries, it means the user provided both the track and pair cuts
        if (sel->GetEntries() == 3) {
          fMuonPairCuts[icut] = (*dqcuts::GetCompositeCut(sel->At(1)->GetName()));
          fMuonRunPairing.push_back(true);
          fMuonNreqObjs.push_back(std::atoi(sel->At(2)->GetName()));
          fMuonPairHistNames[icut] = Form("PairsForwardSEPM_%s_%s", sel->At(0)->GetName(), sel->At(1)->GetName());
        } else {
          fMuonNreqObjs.push_back(std::atoi(sel->At(1)->GetName()));
          fMuonRunPairing.push_back(false);
        }
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);

    // setup the Stats histogram
    fStats.setObject(new TH1I("Statistics", "Stats for DQ triggers", fNBarrelCuts + fNMuonCuts + 2, -2.5, -0.5 + fNBarrelCuts + fNMuonCuts));
    fStats->GetXaxis()->SetBinLabel(1, "Events inspected");
    fStats->GetXaxis()->SetBinLabel(2, "Events selected");
    if (fNBarrelCuts) {
      for (int ib = 3; ib < 3 + fNBarrelCuts; ib++) {
        fStats->GetXaxis()->SetBinLabel(ib, objArray->At(ib - 3)->GetName());
      }
    }
    if (fNMuonCuts) {
      for (int ib = 3 + fNBarrelCuts; ib < 3 + fNBarrelCuts + fNMuonCuts; ib++) {
        fStats->GetXaxis()->SetBinLabel(ib, objArray2->At(ib - 3 - fNBarrelCuts)->GetName());
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    fCCDB->setURL(fConfigCcdbUrl);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    if (fPropMuon) {
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        fCCDB->get<TGeoManager>(geoPath);
      }
    }
    DefineCuts();

    if (fConfigQA) {
      // initialize the variable manager
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      TString histNames = "";
      for (const auto& [key, value] : fBarrelPairHistNames) {
        histNames += value;
        histNames += ";";
      }
      for (const auto& [key, value] : fMuonPairHistNames) {
        histNames += value;
        histNames += ";";
      }
      DefineHistograms(fHistMan, histNames.Data());
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, typename TEvent, typename TTracks, typename TMuons, typename AssocMuons>
  uint64_t runFilterPP(TEvent const& collision,
                       aod::BCsWithTimestamps const& bcs,
                       TTracks const& tracksBarrel,
                       TMuons const& muons,
                       aod::TrackAssoc const& barrelAssocs, AssocMuons const& muonAssocs)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (fCurrentRun != bc.runNumber()) {
      fCurrentRun = bc.runNumber();
      if (fPropMuon) {
        VarManager::SetupMuonMagField();
      }
      grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
      if (grpmag != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
      }
    }
    // Reset the values array and compute event quantities
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<TEventFillMap>(collision); // event properties could be needed for cuts or histogramming

    std::vector<int> objCountersBarrel(fNBarrelCuts, 0); // init all counters to zero
    // count the number of barrel tracks fulfilling each cut
    for (auto trackAssoc : barrelAssocs) {
      auto track = trackAssoc.template track_as<TTracks>();
      for (int i = 0; i < fNBarrelCuts; ++i) {
        if (track.isDQBarrelSelected() & (uint32_t(1) << i)) {
          objCountersBarrel[i] += 1;
        }
      }
    }

    // check which selections require pairing
    uint32_t pairingMask = 0; // in order to know which of the selections actually require pairing
    for (int i = 0; i < fNBarrelCuts; i++) {
      if (fBarrelRunPairing[i]) {
        if (objCountersBarrel[i] > 1) { // pairing has to be enabled and at least two tracks are needed
          pairingMask |= (uint32_t(1) << i);
        }
        objCountersBarrel[i] = 0; // reset counters for selections where pairing is needed (count pairs instead)
      }
    }

    // check which selection should use like sign (LS) (--/++) barrel track pairs
    uint32_t pairingLS = 0; // used to set in which cut setting LS pairs will be analysed
    TString barrelLSstr = fConfigFilterLsBarrelTracksPairs.value;
    std::unique_ptr<TObjArray> objArrayLS(barrelLSstr.Tokenize(","));
    for (int icut = 0; icut < fNBarrelCuts; icut++) {
      TString objStr = objArrayLS->At(icut)->GetName();
      if (!objStr.CompareTo("true")) {
        pairingLS |= (uint32_t(1) << icut);
      }
    }

    // run pairing if there is at least one selection that requires it
    uint32_t pairFilter = 0;
    if (pairingMask > 0) {
      // run pairing on the collision grouped associations
      for (auto& [a1, a2] : combinations(barrelAssocs, barrelAssocs)) {

        // get the tracks from the index stored in the association
        auto t1 = a1.template track_as<TTracks>();
        auto t2 = a2.template track_as<TTracks>();

        // check the pairing mask and that the tracks share a cut bit
        pairFilter = pairingMask & t1.isDQBarrelSelected() & t2.isDQBarrelSelected();
        if (pairFilter == 0) {
          continue;
        }
        // construct the pair and apply pair cuts
        VarManager::FillPair<VarManager::kDecayToEE, TTrackFillMap>(t1, t2); // compute pair quantities
        for (int icut = 0; icut < fNBarrelCuts; icut++) {
          // select like-sign pairs if trigger has set boolean true within fConfigFilterLsBarrelTracksPairs
          if (!(pairingLS & (uint32_t(1) << icut))) {
            if (t1.sign() * t2.sign() > 0) {
              continue;
            }
          }

          if (!(pairFilter & (uint32_t(1) << icut))) {
            continue;
          }
          if (!fBarrelPairCuts[icut].IsSelected(VarManager::fgValues)) {
            continue;
          }
          objCountersBarrel[icut] += 1; // count the pair
          if (fConfigQA) {              // fill histograms if QA is enabled
            fHistMan->FillHistClass(fBarrelPairHistNames[icut].Data(), VarManager::fgValues);
          }
        }
      }
    }

    std::vector<int> objCountersMuon(fNMuonCuts, 0); // init all counters to zero
    // count the number of muon-collision associations fulfilling each selection
    for (auto muon : muonAssocs) {
      for (int i = 0; i < fNMuonCuts; ++i) {
        if (muon.isDQMuonSelected() & (uint32_t(1) << i)) {
          objCountersMuon[i] += 1;
        }
      }
    }

    // check which muon selections require pairing
    pairingMask = 0; // reset the mask for the muons
    for (int i = 0; i < fNMuonCuts; i++) {
      if (fMuonRunPairing[i]) { // pairing has to be enabled and at least two tracks are needed
        if (objCountersMuon[i] > 1) {
          pairingMask |= (uint32_t(1) << i);
        }
        objCountersMuon[i] = 0; // reset counters for selections where pairing is needed (count pairs instead)
      }
    }

    pairingLS = 0; // reset the decisions for muons
    TString musonLSstr = fConfigFilterLsMuonsPairs.value;
    std::unique_ptr<TObjArray> objArrayMuonLS(musonLSstr.Tokenize(","));
    for (int icut = 0; icut < fNMuonCuts; icut++) {
      TString objStr = objArrayMuonLS->At(icut)->GetName();
      if (!objStr.CompareTo("true")) {
        pairingLS |= (uint32_t(1) << icut);
      }
    }

    // run pairing if there is at least one selection that requires it
    pairFilter = 0;
    if (pairingMask > 0) {
      // pairing is done using the collision grouped muon associations
      for (auto& [a1, a2] : combinations(muonAssocs, muonAssocs)) {

        // check the pairing mask and that the tracks share a cut bit
        pairFilter = pairingMask & a1.isDQMuonSelected() & a2.isDQMuonSelected();
        if (pairFilter == 0) {
          continue;
        }

        // get the real muon tracks
        auto t1 = a1.template fwdtrack_as<TMuons>();
        auto t2 = a2.template fwdtrack_as<TMuons>();

        // construct the pair and apply cuts
        VarManager::FillPair<VarManager::kDecayToMuMu, TTrackFillMap>(t1, t2); // compute pair quantities
        if (fPropMuon) {
          VarManager::FillPairPropagateMuon<TTrackFillMap>(t1, t2, collision);
        }
        for (int icut = 0; icut < fNMuonCuts; icut++) {
          // select like-sign pairs if trigger has set boolean true within fConfigFilterLsMuonsPairs
          if (!(pairingLS & (uint32_t(1) << icut))) {
            if (t1.sign() * t2.sign() > 0) {
              continue;
            }
          }
          if (!(pairFilter & (uint32_t(1) << icut))) {
            continue;
          }
          if (!fMuonPairCuts[icut].IsSelected(VarManager::fgValues)) {
            continue;
          }
          objCountersMuon[icut] += 1;
          if (fConfigQA) {
            fHistMan->FillHistClass(fMuonPairHistNames[icut].Data(), VarManager::fgValues);
          }
        }
      }
    }

    // compute the decisions and publish
    uint64_t filter = 0;
    for (int i = 0; i < fNBarrelCuts; i++) {
      if (objCountersBarrel[i] >= fBarrelNreqObjs[i]) {
        filter |= (uint64_t(1) << i);
      }
    }
    for (int i = 0; i < fNMuonCuts; i++) {
      if (objCountersMuon[i] >= fMuonNreqObjs[i]) {
        filter |= (uint64_t(1) << (i + fNBarrelCuts));
      }
    }
    return filter;
  }

  Preslice<TrackAssoc> trackIndicesPerCollision = track_association::collisionId;
  Preslice<FwdTrackAssoc> muonIndicesPerCollision = track_association::collisionId;

  void processFilterPP(MyEventsSelected const& collisions,
                       aod::BCsWithTimestamps const& bcs,
                       MyBarrelTracksSelected const& tracks,
                       MyMuons const& muons,
                       TrackAssoc const& trackAssocs, MyMuonsAssocSelected const& muonAssocs)
  {
    fFiltersMap.clear();
    fCEFPfilters.clear();

    cout << "------------------- filterPP, n assocs barrel/muon :: " << trackAssocs.size() << " / " << muonAssocs.size() << endl;

    uint64_t barrelMask = 0;
    for (int i = 0; i < fNBarrelCuts; i++) {
      barrelMask |= (uint64_t(1) << i);
    }
    uint64_t muonMask = 0;
    for (int i = fNBarrelCuts; i < fNBarrelCuts + fNMuonCuts; i++) {
      muonMask |= (uint64_t(1) << i);
    }

    // Loop over collisions
    int event = 0;
    int eventsFired = 0;
    for (const auto& collision : collisions) {
      // skip those that do not pass our selection
      if (!collision.isDQEventSelected()) {
        event++;
        continue;
      }
      // group the tracks and muons for this collision
      auto groupedTrackIndices = trackAssocs.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      auto groupedMuonIndices = muonAssocs.sliceBy(muonIndicesPerCollision, collision.globalIndex());

      uint64_t filter = 0;
      // if there is at least one track or muon, run the filtering function and compute triggers
      if (groupedTrackIndices.size() > 0 || groupedMuonIndices.size() > 0) {
        filter = runFilterPP<gkEventFillMap, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracks, muons, groupedTrackIndices, groupedMuonIndices);
      }
      if (filter == 0) {
        event++;
        continue;
      }
      eventsFired++;
      // compute the CEPF decisions (this is done in a spacial setup with exactly kNTriggersDQ configured triggers)
      std::vector<bool> decisions(kNTriggersDQ, false); // event decisions to be transmitted to CEFP
      for (int i = 0; i < fNBarrelCuts; i++) {
        if (filter & (uint64_t(1) << i)) {
          if (i < kNTriggersDQ) {
            decisions[i] = true;
          }
        }
      }
      for (int i = fNBarrelCuts; i < fNBarrelCuts + fNMuonCuts; i++) {
        if (filter & (uint64_t(1) << i)) {
          if (i < kNTriggersDQ) {
            decisions[i] = true;
          }
        }
      }

      // if this collision fired at least one input, add it to the map, or if it is there already, update the decisions with a logical OR
      // This may happen in the case when some collisions beyond the iterator are added because they contain ambiguous tracks fired on by another collision
      if (fFiltersMap.find(collision.globalIndex()) == fFiltersMap.end()) {
        fFiltersMap[collision.globalIndex()] = filter;
        fCEFPfilters[collision.globalIndex()] = decisions;
      } else { // this collision was already fired, possible via collision - track association; add as an OR the new decisions
        fFiltersMap[collision.globalIndex()] |= filter;
        for (int i = 0; i < kNTriggersDQ; i++) {
          if (decisions[i]) {
            fCEFPfilters[collision.globalIndex()][i] = true;
          }
        }
      }

      // Now check through the associated tracks / fwdtracks and assign the same filter to their parent collisions
      // This is needed since if a collision was selected because of a track association from a neighbouring collision,
      //   then one needs to select also that collision in order to be able to redo the pairing at analysis time.
      if (filter & barrelMask) {
        for (auto& a : groupedTrackIndices) {
          auto t = a.template track_as<MyBarrelTracksSelected>();
          if (!t.has_collision()) {
            continue;
          }
          auto tColl = t.collisionId();
          if (tColl == collision.globalIndex()) { // track from this collision, nothing to do
            continue;
          } else {
            if (fFiltersMap.find(tColl) == fFiltersMap.end()) {
              fFiltersMap[tColl] = filter;
              fCEFPfilters[tColl] = decisions;
            } else { // this collision was already fired, possible via collision - track association; add as an OR the new decisions
              fFiltersMap[tColl] |= filter;
              for (int i = 0; i < kNTriggersDQ; i++) {
                if (decisions[i]) {
                  fCEFPfilters[tColl][i] = true;
                }
              }
            }
          }
        }
      }
      // Do the same for muons
      if (filter & muonMask) {
        for (auto& a : groupedMuonIndices) {
          auto t = a.template fwdtrack_as<MyMuons>();
          if (!t.has_collision()) {
            continue;
          }
          auto tColl = t.collisionId();
          if (tColl == collision.globalIndex()) { // track from this collision, nothing to do
            continue;
          } else {
            if (fFiltersMap.find(tColl) == fFiltersMap.end()) {
              fFiltersMap[tColl] = filter;
              fCEFPfilters[tColl] = decisions;
            } else { // this collision was already fired, possible via collision - track association; add as an OR the new decisions
              fFiltersMap[tColl] |= filter;
              for (int i = 0; i < kNTriggersDQ; i++) {
                if (decisions[i]) {
                  fCEFPfilters[tColl][i] = true;
                }
              }
            }
          }
        }
      }
      event++;
    }

    // At this point, we have all the non-null decisions for all collisions.
    // we loop again over collisions and create the decision tables
    // NOTE: For the CEFP decisions, decisions are placed in a vector of bool in an ordered way:
    //       start with all configured barrel selections and then continue with those from muons
    //       The configured order has to be in sync with that implemented in the cefp task and can be done
    //       by preparing a dedicated json configuration file
    int totalEventsTriggered = 0;
    for (const auto& collision : collisions) {
      fStats->Fill(-2.0);
      if (!collision.isDQEventSelected()) {
        eventFilter(0);
        dqtable(false, false, false, false, false, false, false);
        continue;
      }
      fStats->Fill(-1.0);

      if (fFiltersMap.find(collision.globalIndex()) == fFiltersMap.end()) {
        eventFilter(0);
        dqtable(false, false, false, false, false, false, false);
      } else {
        totalEventsTriggered++;
        for (int i = 0; i < fNBarrelCuts + fNMuonCuts; i++) {
          if (fFiltersMap[collision.globalIndex()] & (uint32_t(1) << i))
            fStats->Fill(static_cast<float>(i));
        }
        eventFilter(fFiltersMap[collision.globalIndex()]);
        auto dqDecisions = fCEFPfilters[collision.globalIndex()];
        dqtable(dqDecisions[0], dqDecisions[1], dqDecisions[2], dqDecisions[3], dqDecisions[4], dqDecisions[5], dqDecisions[6]);
      }
    }

    cout << "-------------------- In this TF, eventsFired / totalTriggered :: " << eventsFired << "/" << totalEventsTriggered << endl;
  }

  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQFilterPPTask, processFilterPP, "Run filter task", false);
  PROCESS_SWITCH(DQFilterPPTask, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQEventSelectionTask>(cfgc),
    adaptAnalysisTask<DQBarrelTrackSelection>(cfgc),
    adaptAnalysisTask<DQMuonsSelection>(cfgc),
    // adaptAnalysisTask<DQTrackToCollisionAssociation>(cfgc),
    adaptAnalysisTask<DQFilterPPTask>(cfgc)};
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

    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,vtxPbPb");
    }

    if (classStr.Contains("Track")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca");
    }

    if (classStr.Contains("Muon")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon");
    }

    if (classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "vertexing-barrel");
      }
      if (classStr.Contains("Forward")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "dimuon,vertexing-forward");
      }
    }
  }
}
