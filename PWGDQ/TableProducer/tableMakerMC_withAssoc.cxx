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
// Analysis task for skimming MC AODs
// Similar to tableMaker.cxx. The written skimmed data model includes, besides the reconstructed data tables, a skimmed MC stack.
//   The skimmed MC stack includes the MC truth particles corresponding to the list of user specified MC signals (see MCsignal.h)
//    and the MC truth particles corresponding to the reconstructed tracks selected by the specified track cuts on reconstructed data.

#include <iostream>
#include <map>
#include <string>
#include <memory>
#include <vector>
#include "TList.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CCDB/BasicCCDBManager.h"

using std::cout;
using std::endl;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Declare Joins used in the various process functions
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                 aod::McTrackLabels>;
using MyBarrelTracksWithDalitzBits = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                               aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                               aod::pidTPCFullKa, aod::pidTPCFullPr,
                                               aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                               aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                               aod::McTrackLabels, aod::DalitzBits>;
using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                        aod::McTrackLabels>;
using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA>;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MyEventsWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::McCollisionLabels>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::McCollisionLabels>;
using MyEventsWithCentAndMults = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::MultsExtra, aod::McCollisionLabels>;

// Declare bit maps containing information on the table joins content (used as argument in templated functions)
// constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkEventFillMapWithMults = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionMult | VarManager::ObjTypes::CollisionMultExtra;
// constexpr static uint32_t gkEventFillMapWithCent = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent;
constexpr static uint32_t gkEventFillMapWithCentAndMults = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent | VarManager::CollisionMult | VarManager::CollisionMultExtra;
// constexpr static uint32_t gkEventMCFillMap = VarManager::ObjTypes::CollisionMC;
// constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;
// constexpr static uint32_t gkTrackFillMapWithDalitzBits = gkTrackFillMap | VarManager::ObjTypes::DalitzBits;
// constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov;
// constexpr static uint32_t gkMuonFillMapWithAmbi = VarManager::ObjTypes::Muon | VarManager::ObjTypes::AmbiMuon;
// constexpr static uint32_t gkMuonFillMapWithCovAmbi = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov | VarManager::ObjTypes::AmbiMuon;
// constexpr static uint32_t gkTrackFillMapWithAmbi = VarManager::ObjTypes::Track | VarManager::ObjTypes::AmbiTrack;
constexpr static uint32_t gkMFTFillMap = VarManager::ObjTypes::TrackMFT;

struct TableMakerMC {

  Produces<ReducedMCEvents> eventMC;
  Produces<ReducedMCTracks> trackMC;

  Produces<ReducedEvents> event;
  Produces<ReducedEventsExtended> eventExtended;
  Produces<ReducedEventsVtxCov> eventVtxCov;
  Produces<ReducedEventsInfo> eventInfo;
  Produces<ReducedEventsMultPV> multPV;
  Produces<ReducedEventsMultAll> multAll;
  Produces<ReducedMCEventLabels> eventMClabels;

  Produces<ReducedTracksBarrelInfo> trackBarrelInfo;
  Produces<ReducedTracks> trackBasic;
  Produces<ReducedTracksBarrel> trackBarrel;
  Produces<ReducedTracksBarrelCov> trackBarrelCov;
  Produces<ReducedTracksBarrelPID> trackBarrelPID;
  Produces<ReducedTracksAssoc> trackBarrelAssoc;
  Produces<ReducedTracksBarrelLabels> trackBarrelLabels;

  Produces<ReducedMuons> muonBasic;
  Produces<ReducedMuonsExtra> muonExtra;
  Produces<ReducedMuonsCov> muonCov;
  Produces<ReducedMuonsAssoc> muonAssoc;
  Produces<ReducedMuonsLabels> muonLabels; // TODO: enable this once we have fwdtrack labels

  Produces<ReducedMFTs> mftTrack;
  Produces<ReducedMFTsExtra> mftTrackExtra;
  Produces<ReducedMFTAssoc> mftAssoc;

  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TList> fStatsList{"Statistics"}; //! skimming statistics
  HistogramManager* fHistMan;

  Configurable<bool> fIsRun2{"cfgIsRun2", false, "Whether we analyze Run-2 or Run-3 data"};

  // Event and track AnalysisCut configurables
  struct : ConfigurableGroup {
    Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
    Configurable<std::string> fConfigTrackCuts{"cfgBarrelTrackCuts", "jpsiPID1", "barrel track cut"};
    Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  } fConfigCuts;

  // MC signals to be skimmed
  Configurable<std::string> fConfigMCSignals{"cfgMCsignals", "", "Comma separated list of MC signals"};

  // Steer QA output
  struct : ConfigurableGroup {
    Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
    Configurable<bool> fConfigDetailedQA{"cfgDetailedQA", false, "If true, include more QA histograms (BeforeCuts classes)"};
    Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddMCTruthHistogram{"cfgAddMCTruthHistogram", "", "Comma separated list of histograms"};
  } fConfigHistOutput;

  // Selections to be applied as Filter on the Track and FwdTrack
  /*Configurable<float> fConfigBarrelTrackMaxAbsEta{"cfgBarrelMaxAbsEta", 0.9f, "Eta absolute value cut for tracks in the barrel"};
  Configurable<float> fConfigBarrelTrackMinPt{"cfgBarrelMinPt", 0.5f, "Minimum pt for tracks in the barrel"};
  Configurable<float> fConfigBarrelMinTPCncls{"cfgBarrelMinTPCncls", 50.0f, "Minimum TPC cls for tracks in the barrel"};
  Configurable<float> fConfigBarrelMaxTPCchi2{"cfgBarrelMaxTPCchi2", 10.0f, "Maximum TPC chi2/ndf for tracks in the barrel"};
  Configurable<float> fConfigBarrelMaxITSchi2{"cfgBarrelMaxITSchi2", 36.0f, "Maximum ITS chi2/ndf for tracks in the barrel"};
  Configurable<float> fConfigMuonPtLow{"cfgMuonLowPt", 1.0f, "Low pt cut for muons"};
  */

  // CCDB connection configurables
  struct : ConfigurableGroup {
    Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> fGeoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> fGrpMagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> fGrpMagPathRun2{"grpmagPathRun2", "GLO/GRP/GRP", "CCDB path of the GRPObject (Usage for Run 2)"};
  } fConfigCCDB;

  struct : ConfigurableGroup {
    // Track related options
    Configurable<bool> fPropTrack{"cfgPropTrack", true, "Propagate tracks to primary vertex"};
    // Muon related options
    Configurable<bool> fPropMuon{"cfgPropMuon", true, "Propagate muon tracks through absorber (do not use if applying pairing)"};
    Configurable<bool> fRefitGlobalMuon{"cfgRefitGlobalMuon", true, "Correct global muon parameters"};
    Configurable<float> fMuonMatchEtaMin{"cfgMuonMatchEtaMin", -4.0f, "Definition of the acceptance of muon tracks to be matched with MFT"};
    Configurable<float> fMuonMatchEtaMax{"cfgMuonMatchEtaMax", -2.5f, "Definition of the acceptance of muon tracks to be matched with MFT"};
  } fConfigVariousOptions;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  o2::parameters::GRPObject* fGrpMagRun2 = nullptr; // for run 2, we access the GRPObject from GLO/GRP/GRP
  o2::parameters::GRPMagField* fGrpMag = nullptr;   // for run 3, we access GRPMagField from GLO/Config/GRPMagField

  AnalysisCompositeCut* fEventCut;              //! Event selection cut
  std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut> fMuonCuts;  //! Muon track cuts

  bool fDoDetailedQA = false; // Bool to set detailed QA true, if QA is set true
  int fCurrentRun;            // needed to detect if the run changed and trigger update of calibrations etc.

  // list of MCsignal objects
  std::vector<MCSignal> fMCSignals;
  std::map<uint64_t, int> fLabelsMap;
  std::map<uint64_t, int> fLabelsMapReversed;
  std::map<uint64_t, uint16_t> fMCFlags;
  std::map<uint32_t, uint32_t> fCollIndexMap;             // key: old collision index, value: skimmed collision index
  std::map<uint32_t, uint32_t> fTrackIndexMap;            // key: old track global index, value: new track global index
  std::map<uint32_t, uint32_t> fFwdTrackIndexMap;         // key: fwd-track global index, value: new fwd-track global index
  std::map<uint32_t, uint32_t> fFwdTrackIndexMapReversed; // key: new fwd-track global index, value: fwd-track global index
  std::map<uint32_t, uint8_t> fFwdTrackFilterMap;         // key: fwd-track global index, value: fwd-track filter map
  std::map<uint32_t, uint32_t> fMftIndexMap;              // key: MFT tracklet global index, value: new MFT tracklet global index

  void init(o2::framework::InitContext& context)
  {
    // Check whether barrel or muon are enabled
    bool isProcessBCenabled = context.mOptions.get<bool>("processPP");
    bool isBarrelEnabled = (context.mOptions.get<bool>("processPP") || context.mOptions.get<bool>("processPPBarrelOnly") || context.mOptions.get<bool>("processPbPbBarrelOnly"));
    bool isMuonEnabled = (context.mOptions.get<bool>("processPP") || context.mOptions.get<bool>("processPPMuonOnly") || context.mOptions.get<bool>("processPbPbMuonOnly"));
    // Make sure at least one process function is enabled
    if (!(isProcessBCenabled || isBarrelEnabled || isMuonEnabled)) {
      LOG(fatal) << "No process function was enabled for TableMakerMC. Check it out!!!";
    }

    // Define user specified cut
    DefineCuts();

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Only use detailed QA when QA is set true
    if (fConfigHistOutput.fConfigQA && fConfigHistOutput.fConfigDetailedQA) {
      fDoDetailedQA = true;
    }

    // Create the histogram class names to be added to the histogram manager
    TString histClasses = "";
    // Event histograms before any cuts
    if (fDoDetailedQA) {
      histClasses += "Event_BeforeCuts;";
    }
    // Event histograms after cuts and for MC truth collisions
    if (fConfigHistOutput.fConfigQA) {
      histClasses += "Event_AfterCuts;";
      histClasses += "Event_MCTruth;";
    }

    if (isBarrelEnabled) {
      // Barrel track histograms before cuts
      if (fDoDetailedQA) {
        histClasses += "TrackBarrel_BeforeCuts;";
      }
      // Barrel track histograms after cuts; one directory per cut
      if (fConfigHistOutput.fConfigQA) {
        for (auto& cut : fTrackCuts) {
          histClasses += Form("TrackBarrel_%s;", cut.GetName());
        }
      }
    }

    if (isMuonEnabled) {
      // Muon track histograms before cuts
      if (fDoDetailedQA) {
        histClasses += "Muons_BeforeCuts;";
      }
      // Muon track histograms after cuts; one directory per cut
      if (fConfigHistOutput.fConfigQA) {
        for (auto& muonCut : fMuonCuts) {
          histClasses += Form("Muons_%s;", muonCut.GetName());
        }
      }
    }

    // Configure user specified MC signals and setup histogram classes
    TString configNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> objArray(configNamesStr.Tokenize(","));
    if (objArray->GetEntries() > 0) {
      // loop over MC signals
      for (int isig = 0; isig < objArray->GetEntries(); ++isig) {
        MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objArray->At(isig)->GetName());
        if (sig) {
          fMCSignals.push_back(*sig);
          // setup a histogram directory for this MC signal
          if (fConfigHistOutput.fConfigQA) {
            histClasses += Form("MCTruth_%s;", objArray->At(isig)->GetName());
          }
        } else {
          continue;
        }
        if (fDoDetailedQA) {
          if (isBarrelEnabled) {
            // in case of detailed QA, setup histogram directories for each combination of reconstructed track cuts and MC signals
            for (auto& cut : fTrackCuts) {
              histClasses += Form("TrackBarrel_%s_%s;", cut.GetName(), objArray->At(isig)->GetName());
            }
          }
          if (isMuonEnabled) {
            // in case of detailed QA, setup histogram directories for each combination of reconstructed muon cuts and MC signals
            for (auto& cut : fMuonCuts) {
              histClasses += Form("Muons_%s_%s;", cut.GetName(), objArray->At(isig)->GetName());
            }
          }
        }
      }
    }

    DefineHistograms(histClasses);                   // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    // Setup the CCDB
    fCCDB->setURL(fConfigCCDB.fConfigCcdbUrl);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      fCCDB->get<TGeoManager>(fConfigCCDB.fGeoPath);
    }
    fCCDBApi.init(fConfigCCDB.fConfigCcdbUrl.value);
  }

  void DefineCuts()
  {
    // Event cuts
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigCuts.fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

    // Barrel track cuts
    TString cutNamesStr = fConfigCuts.fConfigTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // Muon cuts
    cutNamesStr = fConfigCuts.fConfigMuonCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fMuonCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::MFTTrackAssoc> mfttrackIndicesPerCollision = aod::track_association::collisionId;

  void skimMCCollisions(aod::McCollisions const& mcCollisions)
  {
    // skim MC collisions
    // NOTE: So far, all MC collisions are skimmed. In case there will be filtering based on MC collisions,
    //       one has to do a mapping of the old vs new indices so that the skimmed labels are properly updated.
    VarManager::ResetValues(0, VarManager::kNVars);

    // Loop over MC collisions
    for (auto& mcCollision : mcCollisions) {
      // Get MC collision information into the VarManager
      VarManager::FillEvent<VarManager::ObjTypes::CollisionMC>(mcCollision);
      // Fill histograms
      fHistMan->FillHistClass("Event_MCTruth", VarManager::fgValues);
      // Create the skimmed table entry for this collision
      eventMC(mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
              mcCollision.t(), mcCollision.weight(), mcCollision.impactParameter());
    }
  }

  void skimMCParticles(aod::McParticles const& mcTracks, aod::McCollisions const&)
  {
    // Select MC particles which fulfill at least one of the user specified MC signals
    // In this function we just fill a map with the labels of selected particles, not creating the tables themselves.
    //  The reason is that in the skims we will additionally add any MC label connected to selected reconstructed tracks
    //      which were not selected already via the MC signals

    // Clear the label maps
    fLabelsMap.clear();
    fLabelsMapReversed.clear();
    fMCFlags.clear();

    uint16_t mcflags = static_cast<uint16_t>(0); // flags which will hold the decisions for each MC signal
    int trackCounter = 0;

    for (auto& mctrack : mcTracks) {
      // check all the requested MC signals and fill the decision bit map
      mcflags = 0;
      int i = 0;
      for (auto& sig : fMCSignals) {
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<aod::McParticles>) {
          auto mctrack_raw = mcTracks.rawIteratorAt(mctrack.globalIndex());
          checked = sig.CheckSignal(true, mctrack_raw);
        } else {
          checked = sig.CheckSignal(true, mctrack);
        }
        if (checked) {
          mcflags |= (static_cast<uint16_t>(1) << i);
        }
        i++;
      }
      // if no MC signals were matched, continue
      if (mcflags == 0) {
        continue;
      }

      // If this MC track was not already added to the map, add it now
      if (fLabelsMap.find(mctrack.globalIndex()) == fLabelsMap.end()) {
        fLabelsMap[mctrack.globalIndex()] = trackCounter;
        fLabelsMapReversed[trackCounter] = mctrack.globalIndex();
        fMCFlags[mctrack.globalIndex()] = mcflags;
        trackCounter++;

        // fill histograms for each of the signals, if found
        if (fConfigHistOutput.fConfigQA) {
          VarManager::FillTrackMC(mcTracks, mctrack);
          VarManager::FillEvent<VarManager::ObjTypes::CollisionMC>(mctrack.mcCollision());
          int j = 0;
          for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); signal++, j++) {
            if (mcflags & (static_cast<uint16_t>(1) << j)) {
              fHistMan->FillHistClass(Form("MCTruth_%s", (*signal).GetName()), VarManager::fgValues);
            }
          }
        }
      }
    } // end loop over mc stack
  }

  template <uint32_t TEventFillMap, typename TEvents>
  void skimCollisions(TEvents const& collisions, BCsWithTimestamps const& /*bcs*/)
  {
    // Skim reconstructed collisions which are selected by the user specified cuts

    // Create a collision index map to relate between the "old" AO2D indices and the skimmed ones
    fCollIndexMap.clear();
    int multTPC = -1.0;
    float multFV0A = -1.0;
    float multFV0C = -1.0;
    float multFT0A = -1.0;
    float multFT0C = -1.0;
    float multFDDA = -1.0;
    float multFDDC = -1.0;
    float multZNA = -1.0;
    float multZNC = -1.0;
    int multTracklets = -1.0;
    int multTracksPV = -1.0;
    float centFT0C = -1.0;

    // Loop over collisions
    for (const auto& collision : collisions) {

      // Fill the stats event histogram with the event selection bits
      for (int i = 0; i < o2::aod::evsel::kNsel; i++) {
        if (collision.selection_bit(i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(o2::aod::evsel::kNsel));

      auto bc = collision.template bc_as<BCsWithTimestamps>();
      // store the selection decisions
      uint64_t tag = static_cast<uint64_t>(0);
      // store some more information in the tag
      // if the BC found by event selection does not coincide with the collision.bc(), toggle the first bit
      auto bcEvSel = collision.template foundBC_as<BCsWithTimestamps>();
      if (bcEvSel.globalIndex() != bc.globalIndex()) {
        tag |= (static_cast<uint64_t>(1) << 0);
      }

      // Compute BC and event quantities and fill histograms
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillBC(bc);
      VarManager::FillEvent<TEventFillMap>(collision); // extract event information and place it in the fValues array
      if (collision.has_mcCollision()) {
        VarManager::FillEvent<VarManager::ObjTypes::CollisionMC>(collision.mcCollision());
      }
      if (fDoDetailedQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }

      // fill stats information, before selections
      for (int i = 0; i < o2::aod::evsel::kNsel; i++) {
        if (collision.selection_bit(i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(o2::aod::evsel::kNsel));

      // Apply the user specified event selection
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }

      // fill stats information, after selections
      for (int i = 0; i < o2::aod::evsel::kNsel; i++) {
        if (collision.selection_bit(i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(o2::aod::evsel::kNsel));

      // Fill historams after event cuts
      fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);

      // create the event tables
      event(tag, bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());
      if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0) {
        multTPC = collision.multTPC();
        multFV0A = collision.multFV0A();
        multFV0C = collision.multFV0C();
        multFT0A = collision.multFT0A();
        multFT0C = collision.multFT0C();
        multFDDA = collision.multFDDA();
        multFDDC = collision.multFDDC();
        multZNA = collision.multZNA();
        multZNC = collision.multZNC();
        multTracklets = collision.multTracklets();
        multTracksPV = collision.multNTracksPV();
      }
      if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
        centFT0C = collision.centFT0C();
      }
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                    multTPC, multFV0A, multFV0C, multFT0A, multFT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTracksPV, centFT0C);
      eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());
      eventMClabels(collision.mcCollisionId(), collision.mcMask());
      eventInfo(collision.globalIndex());
      if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMultExtra) > 0) {
        multPV(collision.multNTracksHasITS(), collision.multNTracksHasTPC(), collision.multNTracksHasTOF(), collision.multNTracksHasTRD(),
               collision.multNTracksITSOnly(), collision.multNTracksTPCOnly(), collision.multNTracksITSTPC(), collision.trackOccupancyInTimeRange());
        multAll(collision.multAllTracksTPCOnly(), collision.multAllTracksITSTPC(),
                0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0, 0.0, 0.0, 0.0);
      }

      // add an element for this collision into the map
      fCollIndexMap[collision.globalIndex()] = event.lastIndex();
    }
  }

  template <uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void skimTracks(TEvent const& collision, TTracks const& /*tracks*/, TrackAssoc const& assocs, aod::McParticles const& mcTracks)
  {
    // Skim the barrel track associations
    // Apply track cuts for each collision association and if it passes the cuts, we skim it.
    // NOTE: If selection cuts include conditions on quantities dependent on the associated collision (e.g. DCA),
    //         one track may pass for some association and fail for others.
    //       Tracks are written only once in the skims, even if they contribute to more than one association
    //         so in case of multiple associations, the variables depending on the collision association (e.g. DCA, secondary vertexing, etc)
    //         have to be recomputed at analysis time for each association.

    uint64_t trackFilteringTag = static_cast<uint64_t>(0);
    uint32_t trackTempFilterMap = static_cast<uint32_t>(0);
    uint16_t mcflags = static_cast<uint16_t>(0);
    int trackCounter = fLabelsMap.size();

    // Loop over associations
    for (const auto& assoc : assocs) {
      auto track = assoc.template track_as<TTracks>();

      // If the original collision of this track was not selected for skimming, then we skip this track.
      //  Normally, the filter-pp is selecting all collisions which contain the tracks which contributed to the triggering
      //    of an event, so this is rejecting possibly a few tracks unrelated to the trigger, originally associated with collisions distant in time.
      if (fCollIndexMap.find(track.collisionId()) == fCollIndexMap.end()) {
        continue;
      }

      trackFilteringTag = static_cast<uint64_t>(0);
      trackTempFilterMap = static_cast<uint32_t>(0);

      // Compute track quantities and fill histograms
      VarManager::FillTrack<TTrackFillMap>(track);
      if (fConfigVariousOptions.fPropTrack && (track.collisionId() != collision.globalIndex())) {
        VarManager::FillTrackCollision<TTrackFillMap>(track, collision);
      }
      if (fDoDetailedQA) {
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }

      // apply track cuts and fill histograms
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          trackTempFilterMap |= (static_cast<uint32_t>(1) << i);
          if (fConfigHistOutput.fConfigQA) {
            fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
          }
          (reinterpret_cast<TH1I*>(fStatsList->At(1)))->Fill(static_cast<float>(i));
        }
      }
      if (!trackTempFilterMap) {
        continue;
      }

      // If this track is already present in the index map, it means it was already skimmed,
      // so we just store the association and we skip the track
      if (fTrackIndexMap.find(track.globalIndex()) != fTrackIndexMap.end()) {
        trackBarrelAssoc(fCollIndexMap[collision.globalIndex()], fTrackIndexMap[track.globalIndex()]);
        continue;
      }

      // store V0 and Dalitz bits selection information in the track tag
      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackV0Bits)) { // BIT0-4: V0Bits
        trackFilteringTag |= static_cast<uint64_t>(track.pidbit());
        for (int iv0 = 0; iv0 < 5; iv0++) {
          if (track.pidbit() & (uint8_t(1) << iv0)) {
            (reinterpret_cast<TH1I*>(fStatsList->At(1)))->Fill(fTrackCuts.size() + static_cast<float>(iv0));
          }
        }
      } // end if V0Bits
      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::DalitzBits)) {
        trackFilteringTag |= (static_cast<uint64_t>(track.dalitzBits()) << VarManager::kDalitzBits); // BIT5-12: Dalitz
      }
      trackFilteringTag |= (static_cast<uint64_t>(trackTempFilterMap) << VarManager::kBarrelUserCutsBits); // BIT13-...:  user track filters

      // NOTE: The collision ID that is written in the table is the one originally assigned in the AOD.
      //       However, in data analysis one should loop over associations, so this one should not be used.
      //      In the case of Run2-like analysis, there will be no associations, so this ID will be the one originally assigned in the AO2Ds (updated for the skims)
      uint32_t reducedEventIdx = fCollIndexMap[track.collisionId()];

      // NOTE: trackBarrelInfo stores the index of the collision as in AO2D (for use in some cases where the analysis on skims is done
      //   in workflows where the original AO2Ds are also present)
      trackBarrelInfo(track.collisionId(), collision.posX(), collision.posY(), collision.posZ(), track.globalIndex());
      trackBasic(reducedEventIdx, trackFilteringTag, track.pt(), track.eta(), track.phi(), track.sign(), 0);
      trackBarrel(track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(),
                  track.tpcInnerParam(), track.flags(), track.itsClusterMap(), track.itsChi2NCl(),
                  track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                  track.tpcNClsShared(), track.tpcChi2NCl(),
                  track.trdChi2(), track.trdPattern(), track.tofChi2(),
                  track.length(), track.dcaXY(), track.dcaZ(),
                  track.trackTime(), track.trackTimeRes(), track.tofExpMom(),
                  track.detectorMap());
      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
        trackBarrelCov(track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                       track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                       track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2());
      }
      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackPID)) {
        trackBarrelPID(track.tpcSignal(),
                       track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                       track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                       track.trdSignal());
      }
      fTrackIndexMap[track.globalIndex()] = trackBasic.lastIndex();

      // Check whether the MCParticle corresponding to this reconstructed track was already selected for skimming
      // If not, add it to the skimming map
      if (!track.has_mcParticle()) {
        trackBarrelLabels(-1, 0, 0); // this is the case when there is no matched MCParticle
      } else {
        auto mctrack = track.template mcParticle_as<aod::McParticles>();
        VarManager::FillTrackMC(mcTracks, mctrack);

        mcflags = 0;
        int i = 0; // runs over the MC signals
        int j = 0; // runs over the track cuts
        // check all the specified signals and fill histograms for MC truth matched tracks
        for (auto& sig : fMCSignals) {
          if (sig.CheckSignal(true, mctrack)) {
            mcflags |= (static_cast<uint16_t>(1) << i);
            // If detailed QA is on, fill histograms for each MC signal and track cut combination
            if (fDoDetailedQA) {
              j = 0;
              for (auto& cut : fTrackCuts) {
                if (trackTempFilterMap & (uint8_t(1) << j)) {
                  fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", cut.GetName(), sig.GetName()), VarManager::fgValues); // fill the reconstructed truth
                }
                j++;
              }
            }
          }
          i++;
        }

        // if the MC truth particle corresponding to this reconstructed track is not already written,
        //   add it to the skimmed stack
        if (!(fLabelsMap.find(mctrack.globalIndex()) != fLabelsMap.end())) {
          fLabelsMap[mctrack.globalIndex()] = trackCounter;
          fLabelsMapReversed[trackCounter] = mctrack.globalIndex();
          fMCFlags[mctrack.globalIndex()] = mcflags;
          trackCounter++;
        }
        trackBarrelLabels(fLabelsMap.find(mctrack.globalIndex())->second, track.mcMask(), mcflags);
      }
      // write the skimmed collision - track association
      trackBarrelAssoc(fCollIndexMap[collision.globalIndex()], fTrackIndexMap[track.globalIndex()]);
    } // end loop over associations
  } // end skimTracks

  template <uint32_t TMFTFillMap, typename TEvent>
  void skimMFT(TEvent const& collision, MFTTracks const& /*mfts*/, MFTTrackAssoc const& mftAssocs)
  {
    // Skim MFT tracks
    // So far no cuts are applied here
    for (const auto& assoc : mftAssocs) {
      auto track = assoc.template mfttrack_as<MFTTracks>();

      if (fConfigHistOutput.fConfigQA) {
        VarManager::FillTrack<TMFTFillMap>(track);
        fHistMan->FillHistClass("MftTracks", VarManager::fgValues);
      }

      // write the MFT track global index in the map for skimming (to make sure we have it just once)
      if (fMftIndexMap.find(track.globalIndex()) == fMftIndexMap.end()) {
        uint32_t reducedEventIdx = fCollIndexMap[collision.globalIndex()];
        mftTrack(reducedEventIdx, static_cast<uint64_t>(0), track.pt(), track.eta(), track.phi());
        // TODO: We are not writing the DCA at the moment, because this depends on the collision association
        mftTrackExtra(track.mftClusterSizesAndTrackFlags(), track.sign(), 0.0, 0.0, track.nClusters());

        fMftIndexMap[track.globalIndex()] = mftTrack.lastIndex();
      }
      mftAssoc(fCollIndexMap[collision.globalIndex()], fMftIndexMap[track.globalIndex()]);
    }
  }

  template <uint32_t TMuonFillMap, uint32_t TMFTFillMap, typename TEvent, typename TMuons, typename TMFTTracks>
  void skimMuons(TEvent const& collision, TMuons const& muons, FwdTrackAssoc const& muonAssocs, aod::McParticles const& mcTracks, TMFTTracks const& /*mftTracks*/)
  {
    // Skim the fwd-tracks (muons)
    // Loop over the collision-track associations, recompute track properties depending on the collision assigned, and apply track cuts for selection
    //     Muons are written only once, even if they constribute to more than one association,
    //         which means that in the case of multiple associations, the track parameters are wrong and should be computed again at analysis time.
    uint8_t trackFilteringTag = static_cast<uint8_t>(0);
    uint8_t trackTempFilterMap = static_cast<uint8_t>(0);
    fFwdTrackIndexMapReversed.clear();
    uint16_t mcflags = static_cast<uint16_t>(0);
    int trackCounter = fLabelsMap.size();

    uint32_t offset = muonBasic.lastIndex();
    uint32_t counter = 0;
    for (const auto& assoc : muonAssocs) {
      // get the muon
      auto muon = assoc.template fwdtrack_as<TMuons>();

      trackFilteringTag = uint8_t(0);
      trackTempFilterMap = uint8_t(0);
      VarManager::FillTrack<TMuonFillMap>(muon);
      // NOTE: If a muon is associated to multiple collisions, depending on the selections,
      //       it may be accepted for some associations and rejected for other
      if (fConfigVariousOptions.fPropMuon) {
        VarManager::FillPropagateMuon<TMuonFillMap>(muon, collision);
      }
      // recalculte pDca and global muon kinematics
      if (static_cast<int>(muon.trackType()) < 2 && fConfigVariousOptions.fRefitGlobalMuon) {
        auto muontrack = muon.template matchMCHTrack_as<TMuons>();
        if (muontrack.eta() < fConfigVariousOptions.fMuonMatchEtaMin || muontrack.eta() > fConfigVariousOptions.fMuonMatchEtaMax) {
          continue;
        }
        auto mfttrack = muon.template matchMFTTrack_as<MFTTracks>();
        VarManager::FillTrackCollision<TMuonFillMap>(muontrack, collision);
        VarManager::FillGlobalMuonRefit<TMuonFillMap>(muontrack, mfttrack, collision);
      } else {
        VarManager::FillTrackCollision<TMuonFillMap>(muon, collision);
      }

      if (fDoDetailedQA) {
        fHistMan->FillHistClass("Muons_BeforeCuts", VarManager::fgValues);
      }
      // check the cuts and fill histograms for each fulfilled cut
      int i = 0;
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          trackTempFilterMap |= (uint8_t(1) << i);
          if (fConfigHistOutput.fConfigQA) {
            fHistMan->FillHistClass(Form("Muons_%s", (*cut).GetName()), VarManager::fgValues);
          }
          (reinterpret_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(i));
        }
      }

      // don't skim the muon if no cut has passed
      if (!trackTempFilterMap) {
        continue;
      }
      trackFilteringTag = trackTempFilterMap; // BIT0-7:  user selection cuts

      // update the index map if this is a new muon (it can already exist in the map from a different collision association)
      if (fFwdTrackIndexMap.find(muon.globalIndex()) == fFwdTrackIndexMap.end()) {
        counter++;
        fFwdTrackIndexMap[muon.globalIndex()] = offset + counter;
        fFwdTrackIndexMapReversed[offset + counter] = muon.globalIndex();
        fFwdTrackFilterMap[muon.globalIndex()] = trackFilteringTag;                                                    // store here the filtering tag so we don't repeat the cuts in the second iteration
        if (muon.has_matchMCHTrack() && (fFwdTrackIndexMap.find(muon.matchMCHTrackId()) == fFwdTrackIndexMap.end())) { // write also the matched MCH track
          counter++;
          fFwdTrackIndexMap[muon.matchMCHTrackId()] = offset + counter;
          fFwdTrackIndexMapReversed[offset + counter] = muon.matchMCHTrackId();
          fFwdTrackFilterMap[muon.matchMCHTrackId()] = trackFilteringTag; // store here the filtering tag so we don't repeat the cuts in the second iteration
        }

        if (muon.has_mcParticle()) {
          auto mctrack = muon.template mcParticle_as<aod::McParticles>();
          VarManager::FillTrackMC(mcTracks, mctrack);

          mcflags = 0;
          int i = 0; // runs over the MC signals
          int j = 0; // runs over the track cuts
          // check all the specified signals and fill histograms for MC truth matched tracks
          for (auto& sig : fMCSignals) {
            if (sig.CheckSignal(true, mctrack)) {
              mcflags |= (static_cast<uint16_t>(1) << i);
              if (fDoDetailedQA) {
                for (auto& cut : fMuonCuts) {
                  if (trackTempFilterMap & (uint8_t(1) << j)) {
                    fHistMan->FillHistClass(Form("Muons_%s_%s", cut.GetName(), sig.GetName()), VarManager::fgValues); // fill the reconstructed truth
                  }
                  j++;
                }
              } // end if do detailed QA
            }
            i++;
          } // end loop over MC signals

          // if the MC truth particle corresponding to this reconstructed muon is not already written,
          //   add it to the skimmed stack
          if (!(fLabelsMap.find(mctrack.globalIndex()) != fLabelsMap.end())) {
            fLabelsMap[mctrack.globalIndex()] = trackCounter;
            fLabelsMapReversed[trackCounter] = mctrack.globalIndex();
            fMCFlags[mctrack.globalIndex()] = mcflags;
            trackCounter++;
          }

        }      // end if (has_mcParticle)
      } else { // if muon already in the map, make a bitwise OR with previous existing cuts
        fFwdTrackFilterMap[muon.globalIndex()] |= trackFilteringTag;
      }
      // write the association table
      muonAssoc(fCollIndexMap[collision.globalIndex()], fFwdTrackIndexMap[muon.globalIndex()]);
    } // end loop over assocs

    // Now we have the full index map of selected muons so we can proceed with writing the muon tables
    // Special care needed for the MCH and MFT indices
    for (const auto& [skimIdx, origIdx] : fFwdTrackIndexMapReversed) {
      // get the muon
      auto muon = muons.rawIteratorAt(origIdx);
      uint32_t reducedEventIdx = -1;
      if (muon.has_collision() &&
          fCollIndexMap.find(muon.collisionId()) != fCollIndexMap.end()) { // if the collisionId of this muon was not skimmed, leave the skimmed event index to -1
        reducedEventIdx = fCollIndexMap[muon.collisionId()];
      }
      // NOTE: Currently, one writes the original AO2D momentum-vector (pt, eta and phi) in the tables because we write only one instance of the muon track,
      //       while multiple collision associations (and multiple mom vectors can exist)
      //       The momentum can be recomputed at the analysis time based on the associations written in the skims
      //   So all the information which pertains to collision association or MFT associations should not be taken from the skimmed data, but recomputed at analysis time.
      uint32_t mchIdx = -1;
      uint32_t mftIdx = -1;
      if (muon.trackType() == uint8_t(0) || muon.trackType() == uint8_t(2)) { // MCH-MID (2) or global (0)
        if (fFwdTrackIndexMap.find(muon.matchMCHTrackId()) != fFwdTrackIndexMap.end()) {
          mchIdx = fFwdTrackIndexMap[muon.matchMCHTrackId()];
        }
        if (fMftIndexMap.find(muon.matchMFTTrackId()) != fMftIndexMap.end()) {
          mftIdx = fMftIndexMap[muon.matchMFTTrackId()];
        }
      }
      VarManager::FillTrack<TMuonFillMap>(muon);
      if (fConfigVariousOptions.fPropMuon) {
        VarManager::FillPropagateMuon<TMuonFillMap>(muon, collision);
      }
      // recalculte pDca and global muon kinematics
      if (static_cast<int>(muon.trackType()) < 2 && fConfigVariousOptions.fRefitGlobalMuon) {
        auto muontrack = muon.template matchMCHTrack_as<TMuons>();
        auto mfttrack = muon.template matchMFTTrack_as<MFTTracks>();
        VarManager::FillTrackCollision<TMuonFillMap>(muontrack, collision);
        VarManager::FillGlobalMuonRefit<TMuonFillMap>(muontrack, mfttrack, collision);
      } else {
        VarManager::FillTrackCollision<TMuonFillMap>(muon, collision);
      }
      muonBasic(reducedEventIdx, mchIdx, mftIdx, fFwdTrackFilterMap[muon.globalIndex()], VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], muon.sign(), 0);
      muonExtra(muon.nClusters(), VarManager::fgValues[VarManager::kMuonPDca], VarManager::fgValues[VarManager::kMuonRAtAbsorberEnd],
                muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                muon.matchScoreMCHMFT(),
                muon.mchBitMap(), muon.midBitMap(),
                muon.midBoards(), muon.trackType(), VarManager::fgValues[VarManager::kMuonDCAx], VarManager::fgValues[VarManager::kMuonDCAy],
                muon.trackTime(), muon.trackTimeRes());
      if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {
        muonCov(VarManager::fgValues[VarManager::kX], VarManager::fgValues[VarManager::kY], VarManager::fgValues[VarManager::kZ], VarManager::fgValues[VarManager::kPhi], VarManager::fgValues[VarManager::kTgl], muon.sign() / VarManager::fgValues[VarManager::kPt],
                VarManager::fgValues[VarManager::kMuonCXX], VarManager::fgValues[VarManager::kMuonCXY], VarManager::fgValues[VarManager::kMuonCYY], VarManager::fgValues[VarManager::kMuonCPhiX], VarManager::fgValues[VarManager::kMuonCPhiY], VarManager::fgValues[VarManager::kMuonCPhiPhi],
                VarManager::fgValues[VarManager::kMuonCTglX], VarManager::fgValues[VarManager::kMuonCTglY], VarManager::fgValues[VarManager::kMuonCTglPhi], VarManager::fgValues[VarManager::kMuonCTglTgl], VarManager::fgValues[VarManager::kMuonC1Pt2X], VarManager::fgValues[VarManager::kMuonC1Pt2Y],
                VarManager::fgValues[VarManager::kMuonC1Pt2Phi], VarManager::fgValues[VarManager::kMuonC1Pt2Tgl], VarManager::fgValues[VarManager::kMuonC1Pt21Pt2]);
      }
      if (muon.has_mcParticle()) {
        auto mctrack = muon.template mcParticle_as<aod::McParticles>();
        muonLabels(fLabelsMap.find(mctrack.globalIndex())->second, muon.mcMask(), mcflags);
      } else {
        muonLabels(-1, 0, 0);
      }
    } // end loop over selected muons
  } // end skimMuons

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, uint32_t TMFTFillMap, typename TEvents, typename TTracks,
            typename TMuons, typename TMFTTracks, typename TTrackAssoc, typename TFwdTrackAssoc, typename TMFTTrackAssoc>
  void fullSkimming(TEvents const& collisions, BCsWithTimestamps const& bcs,
                    TTracks const& tracksBarrel, TMuons const& muons, TMFTTracks const& mftTracks,
                    TTrackAssoc const& trackAssocs, TFwdTrackAssoc const& fwdTrackAssocs, TMFTTrackAssoc const& mftAssocs,
                    aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    // Check whether the run changed and update CCDB if it did
    if (bcs.size() > 0 && fCurrentRun != bcs.begin().runNumber()) {
      if (fIsRun2 == true) {
        fGrpMagRun2 = fCCDB->getForTimeStamp<o2::parameters::GRPObject>(fConfigCCDB.fGrpMagPathRun2, bcs.begin().timestamp());
        if (fGrpMagRun2 != nullptr) {
          o2::base::Propagator::initFieldFromGRP(fGrpMagRun2);
        }
      } else {
        fGrpMag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fConfigCCDB.fGrpMagPath, bcs.begin().timestamp());
        if (fGrpMag != nullptr) {
          o2::base::Propagator::initFieldFromGRP(fGrpMag);
        }
      }
      std::map<string, string> metadataRCT, header;
      header = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", bcs.begin().runNumber()), metadataRCT, -1);
      uint64_t sor = std::atol(header["SOR"].c_str());
      uint64_t eor = std::atol(header["EOR"].c_str());
      VarManager::SetSORandEOR(sor, eor);

      fCurrentRun = bcs.begin().runNumber();
    }

    // skim MC Collisions
    eventMC.reserve(mcCollisions.size());
    skimMCCollisions(mcCollisions);

    // select MC particles to be written using the specified MC signals
    // NOTE: tables are not written at this point, only label maps are being created
    skimMCParticles(mcParticles, mcCollisions);

    // skim collisions
    event.reserve(collisions.size());
    eventExtended.reserve(collisions.size());
    eventVtxCov.reserve(collisions.size());
    eventMClabels.reserve(collisions.size());
    eventInfo.reserve(collisions.size());
    skimCollisions<TEventFillMap>(collisions, bcs);
    if (fCollIndexMap.size() == 0) {
      return;
    }

    // Clear index map and reserve memory for barrel tables
    if constexpr (static_cast<bool>(TTrackFillMap)) {
      fTrackIndexMap.clear();
      trackBarrelInfo.reserve(tracksBarrel.size());
      trackBasic.reserve(tracksBarrel.size());
      trackBarrel.reserve(tracksBarrel.size());
      trackBarrelCov.reserve(tracksBarrel.size());
      trackBarrelPID.reserve(tracksBarrel.size());
      trackBarrelAssoc.reserve(tracksBarrel.size());
      trackBarrelLabels.reserve(tracksBarrel.size());
    }

    // Clear index map and reserve memory for MFT tables
    if constexpr (static_cast<bool>(TMFTFillMap)) {
      fMftIndexMap.clear();
      mftTrack.reserve(mftTracks.size());
      mftTrackExtra.reserve(mftTracks.size());
      mftAssoc.reserve(mftTracks.size());
    }

    // Clear index map and reserve memory for muon tables
    if constexpr (static_cast<bool>(TMuonFillMap)) {
      fFwdTrackIndexMap.clear();
      fFwdTrackFilterMap.clear();
      muonBasic.reserve(muons.size());
      muonExtra.reserve(muons.size());
      muonCov.reserve(muons.size());
      muonAssoc.reserve(muons.size());
      muonLabels.reserve(muons.size());
    }

    // loop over selected collisions and select the tracks and fwd tracks to be skimmed
    if (fCollIndexMap.size() > 0) {
      for (auto const& [origIdx, skimIdx] : fCollIndexMap) {
        auto collision = collisions.rawIteratorAt(origIdx);
        // group the tracks and muons for this collision
        if constexpr (static_cast<bool>(TTrackFillMap)) {
          auto groupedTrackIndices = trackAssocs.sliceBy(trackIndicesPerCollision, origIdx);
          skimTracks<TTrackFillMap>(collision, tracksBarrel, groupedTrackIndices, mcParticles);
        }
        if constexpr (static_cast<bool>(TMFTFillMap)) {
          auto groupedMFTIndices = mftAssocs.sliceBy(mfttrackIndicesPerCollision, origIdx);
          skimMFT<TMFTFillMap>(collision, mftTracks, groupedMFTIndices);
        }
        if constexpr (static_cast<bool>(TMuonFillMap)) {
          if constexpr (static_cast<bool>(TMFTFillMap)) {
            auto groupedMuonIndices = fwdTrackAssocs.sliceBy(fwdtrackIndicesPerCollision, origIdx);
            skimMuons<TMuonFillMap, TMFTFillMap>(collision, muons, groupedMuonIndices, mcParticles, mftTracks);
          } else {
            auto groupedMuonIndices = fwdTrackAssocs.sliceBy(fwdtrackIndicesPerCollision, origIdx);
            skimMuons<TMuonFillMap, 0u>(collision, muons, groupedMuonIndices, mcParticles, nullptr);
          }
        }
      } // end loop over skimmed collisions
    }

    // Loop over the label map, create the mother/daughter relationships if these exist and write the skimmed MC stack
    for (const auto& [newLabel, oldLabel] : fLabelsMapReversed) {
      auto mctrack = mcParticles.iteratorAt(oldLabel);
      uint16_t mcflags = fMCFlags.find(oldLabel)->second;

      std::vector<int> mothers;
      if (mctrack.has_mothers()) {
        for (auto& m : mctrack.mothersIds()) {
          if (m < mcParticles.size()) { // protect against bad mother indices
            if (fLabelsMap.find(m) != fLabelsMap.end()) {
              mothers.push_back(fLabelsMap.find(m)->second);
            }
          } else {
            cout << "Mother label (" << m << ") exceeds the McParticles size (" << mcParticles.size() << ")" << endl;
            cout << " Check the MC generator" << endl;
          }
        }
      }

      // TODO: Check that the daughter slice in the skimmed table works as expected
      //       Note that not all daughters from the original table are preserved in the skimmed MC stack
      std::vector<int> daughters;
      if (mctrack.has_daughters()) {
        for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
          // TODO: remove this check as soon as issues with MC production are fixed
          if (d < mcParticles.size()) { // protect against bad daughter indices
            if (fLabelsMap.find(d) != fLabelsMap.end()) {
              daughters.push_back(fLabelsMap.find(d)->second);
            }
          } else {
            cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcParticles.size() << ")" << endl;
            cout << " Check the MC generator" << endl;
          }
        }
      }
      int daughterRange[2] = {-1, -1};
      if (daughters.size() > 0) {
        daughterRange[0] = daughters[0];
        daughterRange[1] = daughters[daughters.size() - 1];
      }

      // NOTE: Here we assume that MC collisions are not filtered, so there is no new vs old index map for translation
      trackMC(mctrack.mcCollision().globalIndex(), mctrack.pdgCode(), mctrack.statusCode(), mctrack.flags(),
              mothers, daughterRange,
              mctrack.weight(), mctrack.pt(), mctrack.eta(), mctrack.phi(), mctrack.e(),
              mctrack.vx(), mctrack.vy(), mctrack.vz(), mctrack.vt(), mcflags);
      for (unsigned int isig = 0; isig < fMCSignals.size(); isig++) {
        if (mcflags & (static_cast<uint16_t>(1) << isig)) {
          (reinterpret_cast<TH1I*>(fStatsList->At(3)))->Fill(static_cast<float>(isig));
        }
      }
      if (mcflags == 0) {
        (reinterpret_cast<TH1I*>(fStatsList->At(3)))->Fill(static_cast<float>(fMCSignals.size()));
      }
    } // end loop over labels
  }

  void DefineHistograms(TString histClasses)
  {
    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      if (fConfigHistOutput.fConfigQA) {
        fHistMan->AddHistClass(classStr.Data());
      }

      TString histEventName = fConfigHistOutput.fConfigAddEventHistogram.value;
      if (classStr.Contains("Event")) {
        if (fConfigHistOutput.fConfigQA && !classStr.Contains("MCTruth")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", histEventName);
        } else {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", "generator");
        }
      }

      TString histTrackName = fConfigHistOutput.fConfigAddTrackHistogram.value;
      if (classStr.Contains("Track")) {
        if (fConfigHistOutput.fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histTrackName);
        }
      }

      TString histMuonName = fConfigHistOutput.fConfigAddMuonHistogram.value;
      if (classStr.Contains("Muons")) {
        if (fConfigHistOutput.fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histMuonName);
        }
      }

      TString histMCTruthName = fConfigHistOutput.fConfigAddMCTruthHistogram.value;
      if (classStr.Contains("MCTruth") && !classStr.Contains("Event")) {
        if (fConfigHistOutput.fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "mctruth_track", histMCTruthName);
        }
      }
    }

    // create statistics histograms (event, tracks, muons, MCsignals)
    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);
    std::vector<TString> eventLabels{"BCs", "Collisions before filtering", "Before cuts", "After cuts"};
    TH2I* histEvents = new TH2I("EventStats", "Event statistics", eventLabels.size(), -0.5, eventLabels.size() - 0.5, o2::aod::evsel::kNsel + 1, -0.5, (float)o2::aod::evsel::kNsel + 0.5);
    int ib = 1;
    for (auto label = eventLabels.begin(); label != eventLabels.end(); label++, ib++) {
      histEvents->GetXaxis()->SetBinLabel(ib, (*label).Data());
    }
    for (int ib = 1; ib <= o2::aod::evsel::kNsel; ib++) {
      histEvents->GetYaxis()->SetBinLabel(ib, o2::aod::evsel::selectionLabels[ib - 1]);
    }
    histEvents->GetYaxis()->SetBinLabel(o2::aod::evsel::kNsel + 1, "Total");
    fStatsList->Add(histEvents);

    // Track statistics: one bin for each track selection and 5 bins for V0 tags (gamma, K0s, Lambda, anti-Lambda, Omega)
    TH1I* histTracks = new TH1I("TrackStats", "Track statistics", fTrackCuts.size() + 5.0, -0.5, fTrackCuts.size() - 0.5 + 5.0);
    ib = 1;
    for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, ib++) {
      histTracks->GetXaxis()->SetBinLabel(ib, (*cut).GetName());
    }
    const char* v0TagNames[5] = {"Photon conversion", "K^{0}_{s}", "#Lambda", "#bar{#Lambda}", "#Omega"};
    for (int ib = 0; ib < 5; ib++) {
      histTracks->GetXaxis()->SetBinLabel(fTrackCuts.size() + 1 + ib, v0TagNames[ib]);
    }
    fStatsList->Add(histTracks);
    TH1I* histMuons = new TH1I("MuonStats", "Muon statistics", fMuonCuts.size(), -0.5, fMuonCuts.size() - 0.5);
    ib = 1;
    for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, ib++) {
      histMuons->GetXaxis()->SetBinLabel(ib, (*cut).GetName());
    }
    fStatsList->Add(histMuons);
    TH1I* histMCsignals = new TH1I("MCsignals", "MC signals", fMCSignals.size() + 1, -0.5, fMCSignals.size() - 0.5 + 1.0);
    ib = 1;
    for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); signal++, ib++) {
      histMCsignals->GetXaxis()->SetBinLabel(ib, (*signal).GetName());
    }
    histMCsignals->GetXaxis()->SetBinLabel(fMCSignals.size() + 1, "Others (matched to reco tracks)");
    fStatsList->Add(histMCsignals);
  }

  void processPP(MyEventsWithMults const& collisions, aod::BCsWithTimestamps const& bcs,
                 MyBarrelTracksWithCov const& tracksBarrel, MyMuonsWithCov const& tracksMuon, aod::MFTTracks const& mftTracks,
                 aod::TrackAssoc const& trackAssocs, aod::FwdTrackAssoc const& fwdTrackAssocs, aod::MFTTrackAssoc const& mftAssocs,
                 aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    fullSkimming<gkEventFillMapWithMults, gkTrackFillMapWithCov, gkMuonFillMapWithCov, gkMFTFillMap>(collisions, bcs, tracksBarrel, tracksMuon, mftTracks, trackAssocs, fwdTrackAssocs, mftAssocs, mcCollisions, mcParticles);
  }

  void processPPBarrelOnly(MyEventsWithMults const& collisions, aod::BCsWithTimestamps const& bcs,
                           MyBarrelTracksWithCov const& tracksBarrel, aod::TrackAssoc const& trackAssocs,
                           aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    fullSkimming<gkEventFillMapWithMults, gkTrackFillMapWithCov, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, nullptr, trackAssocs, nullptr, nullptr, mcCollisions, mcParticles);
  }

  void processPPMuonOnly(MyEventsWithMults const& collisions, aod::BCsWithTimestamps const& bcs,
                         MyMuonsWithCov const& tracksMuon, aod::MFTTracks const& mftTracks,
                         aod::FwdTrackAssoc const& fwdTrackAssocs, aod::MFTTrackAssoc const& mftAssocs,
                         aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    fullSkimming<gkEventFillMapWithMults, 0u, gkMuonFillMapWithCov, gkMFTFillMap>(collisions, bcs, nullptr, tracksMuon, mftTracks, nullptr, fwdTrackAssocs, mftAssocs, mcCollisions, mcParticles);
  }

  void processPbPb(MyEventsWithCentAndMults const& collisions, aod::BCsWithTimestamps const& bcs,
                   MyBarrelTracksWithCov const& tracksBarrel, MyMuonsWithCov const& tracksMuon, aod::MFTTracks const& mftTracks,
                   aod::TrackAssoc const& trackAssocs, aod::FwdTrackAssoc const& fwdTrackAssocs, aod::MFTTrackAssoc const& mftAssocs,
                   aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMapWithCov, gkMuonFillMapWithCov, gkMFTFillMap>(collisions, bcs, tracksBarrel, tracksMuon, mftTracks, trackAssocs, fwdTrackAssocs, mftAssocs, mcCollisions, mcParticles);
  }

  void processPbPbBarrelOnly(MyEventsWithCentAndMults const& collisions, aod::BCsWithTimestamps const& bcs,
                             MyBarrelTracksWithCov const& tracksBarrel, aod::TrackAssoc const& trackAssocs,
                             aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMapWithCov, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, nullptr, trackAssocs, nullptr, nullptr, mcCollisions, mcParticles);
  }

  void processPbPbMuonOnly(MyEventsWithCentAndMults const& collisions, aod::BCsWithTimestamps const& bcs,
                           MyMuonsWithCov const& tracksMuon, aod::MFTTracks const& mftTracks,
                           aod::FwdTrackAssoc const& fwdTrackAssocs, aod::MFTTrackAssoc const& mftAssocs,
                           aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, 0u, gkMuonFillMapWithCov, gkMFTFillMap>(collisions, bcs, nullptr, tracksMuon, mftTracks, nullptr, fwdTrackAssocs, mftAssocs, mcCollisions, mcParticles);
  }

  // Process the BCs and store stats for luminosity retrieval -----------------------------------------------------------------------------------
  void processOnlyBCs(soa::Join<aod::BCs, aod::BcSels>::iterator const& bc)
  {
    for (int i = 0; i < o2::aod::evsel::kNsel; i++) {
      if (bc.alias_bit(i) > 0) {
        (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(0.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(0.0, static_cast<float>(o2::aod::evsel::kNsel));
  }

  PROCESS_SWITCH(TableMakerMC, processPP, "Produce both barrel and muon skims, pp settings", false);
  PROCESS_SWITCH(TableMakerMC, processPPBarrelOnly, "Produce only barrel skims, pp settings ", false);
  PROCESS_SWITCH(TableMakerMC, processPPMuonOnly, "Produce only muon skims, pp settings", false);
  PROCESS_SWITCH(TableMakerMC, processPbPb, "Produce both barrel and muon skims, PbPb settings", false);
  PROCESS_SWITCH(TableMakerMC, processPbPbBarrelOnly, "Produce only barrel skims, PbPb settings", false);
  PROCESS_SWITCH(TableMakerMC, processPbPbMuonOnly, "Produce only muon skims, PbPb settings", false);
  PROCESS_SWITCH(TableMakerMC, processOnlyBCs, "Analyze the BCs to store sampled lumi", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // TODO: For now TableMakerMC works just for PbPb (cent table is present)
  //      Implement workflow arguments for pp/PbPb and possibly merge the task with tableMaker.cxx
  return WorkflowSpec{
    adaptAnalysisTask<TableMakerMC>(cfgc)};
}
