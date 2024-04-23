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
///
/// \brief this is a code for the kstarqa resonance
/// \author prottay das, sawan
/// \since 13/03/2024

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct kstarqa {

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  SliceCache cache;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Confugrable for QA histograms
  Configurable<bool> QA{"QA", false, "QA"};
  Configurable<bool> QAbefore{"QAbefore", true, "QAbefore"};
  Configurable<bool> QAafter{"QAafter", true, "QAafter"};
  Configurable<bool> onlyTOF{"onlyTOF", false, "only TOF tracks"};
  Configurable<bool> onlyTOFHIT{"onlyTOFHIT", false, "accept only TOF hit tracks at high pt"};

  // Configurables for track selections
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPCPi{"nsigmacutTPCPi", 3.0, "Value of the TPC Nsigma cut for pions"};
  Configurable<float> nsigmaCutTPCKa{"nsigmacutTPCKa", 3.0, "Value of the TPC Nsigma cut for kaons"};
  Configurable<float> nsigmaCutTOFPi{"nsigmacutTOFPi", 3.0, "Value of the TOF Nsigma cut for pions"};
  Configurable<float> nsigmaCutTOFKa{"nsigmacutTOFKa", 3.0, "Value of the TOF Nsigma cut for kaons"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the Combined Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", false, "cfgMultFT0"};
  Configurable<bool> cfgCentFT0C{"cfgCentFT0C", true, "cfgCentFT0C"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> ismanulatrksel{"ismanulatrksel", true, "all manual track selection cuts"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<float> cfgRCRFC{"cfgRCRFC", 0.8f, "Crossed Rows to Findable Clusters"};
  Configurable<float> cfgITSChi2NCl{"cfgITSChi2NCl", 36.0, "ITS Chi2/NCl"};
  Configurable<float> cfgTPCChi2NCl{"cfgTPCChi2NCl", 4.0, "TPC Chi2/NCl"};
  Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
  Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
  ConfigurableAxis cMixMultBins{"cMixMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Mixing bins - multiplicity"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", false, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                       // kGoldenChi2 | kDCAxy | kDCAz

  // Event selection configurables
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> piluprejection{"piluprejection", false, "Pileup rejection"};
  Configurable<bool> goodzvertex{"goodzvertex", false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference."};
  Configurable<bool> itstpctracks{"itstpctracks", false, "selects collisions with at least one ITS-TPC track,"};
  Configurable<bool> MID{"MID", false, "Misidentification of tracks"};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  Configurable<int> nBinsinvMass{"nBinsinvMass", 360, "N bins in invMass histos"};
  Configurable<float> invMassbinlow{"invMassbinlow", 0.6, "invMass bin low"};
  Configurable<float> invMassbinhigh{"invMassbinhigh", 1.5, "invMass bin high"};
  Configurable<int> nBinspT{"nBinspT", 200, "N bins in pT histos"};
  Configurable<float> pTbinlow{"pTbinlow", 0.0, "pT bin low"};
  Configurable<float> pTbinhigh{"pTbinhigh", 20.0, "pT bin high"};
  ConfigurableAxis binsMultPlot{"binsCent", {200, 0.0f, 200.0f}, "Binning of the centrality axis for plots"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm] for plots"};
    AxisSpec ptAxis = {nBinspT, pTbinlow, pTbinhigh, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invmassAxis = {nBinsinvMass, invMassbinlow, invMassbinhigh, "Invariant mass (GeV/#it{c}^{2})"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmult", "Centrality distribution", kTH1F, {{binsMultPlot}});

    // for primary tracks
    if (QAbefore && QAafter) {
      histos.add("hNsigmaPionTPC_before", "NsigmaPion TPC distribution before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_before", "NsigmaPion TOF distribution before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      histos.add("hNsigmaKaonTPC_before", "NsigmaKaon TPC distribution before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      histos.add("hNsigmaKaonTOF_before", "NsigmaKaon TOF distribution before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});

      histos.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hCRFC_after", "CRFC after distribution", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("hCRFC_before", "CRFC before distribution", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("hNsigmaPionTPC_after", "NsigmaPion TPC distribution", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_after", "NsigmaPion TOF distribution", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      histos.add("hNsigmaKaonTPC_after", "NsigmaKaon TPC distribution", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      histos.add("hNsigmaKaonTOF_after", "NsigmaKaon TOF distribution", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
    }

    // KStar histograms
    histos.add("h3KstarInvMassUnlikeSign", "Invariant mass of kstar meson Unlike Sign", kTHnSparseF, {{binsMultPlot}, {ptAxis}, {invmassAxis}}, true);
    histos.add("h3KstarInvMasslikeSign", "Invariant mass of kstar meson like Sign", kTHnSparseF, {{binsMultPlot}, {ptAxis}, {invmassAxis}}, true);
    histos.add("h3KstarInvMassMixed", "Invariant mass of kstar meson Mixed", kTHnSparseF, {{binsMultPlot}, {ptAxis}, {invmassAxis}}, true);

    // MC generated histograms
    histos.add("k892Gen", "pT distribution of True MC K(892)0", kTH1D, {ptAxis});
    // histos.add("k892GenAnti", "pT distribution of True MC Anti-K(892)0", kTH1D, {ptAxis});
    // Reconstructed MC histogram
    histos.add("h3KstarRec", "pT of reconstructed kstar", kTH1D, {ptAxis});
    histos.add("h1KstarRecMass", "Invariant mass of kstar meson", kTH1D, {invmassAxis});
    // histos.add("h1KstarRecpt", "pT of kstar meson", kTH1D, {ptAxis});
    histos.add("h1genmass", "Invariant mass of generated kstar meson", kTH1D, {invmassAxis});
    histos.add("h1recpt", "pT of generated kstar meson", kTH1D, {ptAxis});
    histos.add("events_check", "No. of events in the reconstructed and generated MC", kTH1F, {{6, 0, 6}});
  }

  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass(); // FIXME: Get from the common header
  double massKa = o2::constants::physics::MassKPlus;
  ROOT::Math::PtEtaPhiMVector CKSVector;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    // if (iscustomDCAcut && (!candidate.isGlobalTrack() || !candidate.isPVContributor() || candidate.itsNCls() < cfgITScluster)) {
    //   return false;
    // }
    // if (ismanulatrksel &&
    //     !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster && candidate.tpcCrossedRowsOverFindableCls() > cfgRCRFC)) {
    //   return false;
    // }
    if (std::abs(candidate.pt()) < cfgCutPT)
      return false;
    if (std::abs(candidate.dcaXY()) > cfgCutDCAxy)
      return false;
    if (std::abs(candidate.dcaZ()) > cfgCutDCAz)
      return false;
    if (candidate.tpcCrossedRowsOverFindableCls() < cfgRCRFC)
      return false;
    if (candidate.itsNCls() < cfgITScluster)
      return false;
    if (candidate.tpcNClsFound() < cfgTPCcluster)
      return false;
    if (candidate.itsChi2NCl() >= cfgITSChi2NCl)
      return false;
    if (candidate.tpcChi2NCl() >= cfgTPCChi2NCl)
      return false;
    if (cfgUseITSRefit && !candidate.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !candidate.passedTPCRefit())
      return false;
    if (cfgPVContributor && !candidate.isPVContributor())
      return false;
    if (cfgPrimaryTrack && !candidate.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !candidate.isGlobalTrackWoDCA())
      return false;
    if (cfgGlobalTrack && !candidate.isGlobalTrack())
      return false;

    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (onlyTOF) {
        // LOG(info) << "************I am inside ONLYTOF****************";
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOFPi) {
          return true;
        }
      } else if (onlyTOFHIT) {
        // LOG(info) << "************I am inside ONLYTOFHIT****************";
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOFPi) {
          return true;
        }
        if (!candidate.hasTOF() &&
            std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      } else {
        // LOG(info) << "************I am neither in ONLYTOF or ONLYTOFHIT****************";
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      }
    } else if (PID == 1) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOFKa) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOFKa) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename T>
  bool MIDselectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < 3.0) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < 3.0) {
          return true;
        }
        if (!candidate.hasTOF() &&
            std::abs(candidate.tpcNSigmaPi()) < 3.0) {
          return true;
        }
      } else {
        // LOG(info) << "************I am neither in ONLYTOF or ONLYTOFHIT****************";
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (3.0 * 3.0)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < 3.0) {
          return true;
        }
      }
    } else if (PID == 1) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < 3.0) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < 3.0) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (3.0 * 3.0)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0) {
          return true;
        }
      }
    } else if (PID == 2) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < 3.0) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < 3.0) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0) {
          return true;
        }
      } else {
        // LOG(info) << "************I am neither in ONLYTOF or ONLYTOFHIT****************";
        if (candidate.hasTOF() && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (3.0 * 3.0)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0) {
          return true;
        }
      }
    }
    return false;
  }

  array<float, 3> pvec0;
  array<float, 3> pvec1;

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelectionExtension>>;
  using V0TrackCandidate = aod::V0Datas;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::McTrackLabels>>;

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for ME mixing"};
  // ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {10, 0, 100}, "multiplicity percentile for ME mixing"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity  for bin for ME mixing"};

  using BinningTypeTPCMultiplicity = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  // using BinningTypeVertexContributor =
  // ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  using BinningTypeCentralityM = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};

  SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};

  double totmomka = 0.0;
  double totmompi = 0.0;
  double totmomkamix = 0.0;
  double totmompimix = 0.0;
  // double openingangle = 0.0;
  // double openinganglemix = 0.0;

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)

  {

    if (!collision.sel8()) {
      return;
    }

    if (timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }

    if (piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    if (itstpctracks && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }

    std::vector<ROOT::Math::PtEtaPhiMVector> pions, kaons;
    std::vector<int64_t> PionIndex = {};
    std::vector<int64_t> KaonIndex = {};
    std::vector<int64_t> PioncollIndex = {};
    std::vector<int64_t> KaoncollIndex = {};
    std::vector<int64_t> PionSign = {};
    std::vector<int64_t> KaonSign = {};

    std::vector<double_t> PionPx = {};
    std::vector<double_t> KaonPx = {};
    std::vector<double_t> PionPy = {};
    std::vector<double_t> KaonPy = {};
    std::vector<double_t> PionPz = {};
    std::vector<double_t> KaonPz = {};
    std::vector<double_t> PionP = {};
    std::vector<double_t> KaonP = {};

    float multiplicity = 0.0f;
    /*  if (cfgMultFT0)
      multiplicity = collision.multZeqFT0A() + collision.multZeqFT0C();
    if (cfgMultFT0 == 0 && cfgCentFT0C == 1)
      multiplicity = collision.centFT0C();
      if (cfgMultFT0 == 0 && cfgCentFT0C == 0)*/
    multiplicity = collision.centFT0M();

    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmult"), multiplicity);

    for (auto track1 : tracks) {

      if (QAbefore) {
        histos.fill(HIST("hNsigmaKaonTPC_before"), track1.pt(), track1.tpcNSigmaKa());
        histos.fill(HIST("hNsigmaKaonTOF_before"), track1.pt(), track1.tofNSigmaKa());
        histos.fill(HIST("hCRFC_before"), track1.tpcCrossedRowsOverFindableCls());
      }

      if (!selectionTrack(track1)) {
        continue;
      }
      if (!selectionPID(track1, 1)) // kaon
        continue;

      if (MID) {
        if (MIDselectionPID(track1, 0)) // misidentified as pion
          continue;

        if (MIDselectionPID(track1, 2)) // misidentified as proton
          continue;
      }

      if (QAafter) {
        histos.fill(HIST("hEta_after"), track1.eta());
        histos.fill(HIST("hCRFC_after"), track1.tpcCrossedRowsOverFindableCls());
        histos.fill(HIST("hNsigmaKaonTPC_after"), track1.pt(), track1.tpcNSigmaKa());
        histos.fill(HIST("hNsigmaKaonTOF_after"), track1.pt(), track1.tofNSigmaKa());
      }

      totmomka = TMath::Sqrt(track1.px() * track1.px() + track1.py() * track1.py() + track1.pz() * track1.pz());

      ROOT::Math::PtEtaPhiMVector temp1(track1.pt(), track1.eta(), track1.phi(), massKa);
      kaons.push_back(temp1);
      KaonIndex.push_back(track1.globalIndex());
      KaoncollIndex.push_back(track1.collisionId());
      KaonSign.push_back(track1.sign());
      KaonPx.push_back(track1.px());
      KaonPy.push_back(track1.py());
      KaonPz.push_back(track1.pz());
      KaonP.push_back(totmomka);

    } // track loop ends

    for (auto track2 : tracks) {

      if (QAbefore) {
        histos.fill(HIST("hNsigmaPionTPC_before"), track2.pt(), track2.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_before"), track2.pt(), track2.tofNSigmaPi());
      }

      if (!selectionTrack(track2)) {
        continue;
      }

      if (!selectionPID(track2, 0)) // pion
        continue;

      if (MID) {
        if (MIDselectionPID(track2, 1)) // misidentified as kaon
          continue;
      }
      if (QAafter) {
        histos.fill(HIST("hNsigmaPionTPC_after"), track2.pt(), track2.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_after"), track2.pt(), track2.tofNSigmaPi());
      }

      totmompi = TMath::Sqrt(track2.px() * track2.px() + track2.py() * track2.py() + track2.pz() * track2.pz());

      ROOT::Math::PtEtaPhiMVector temp2(track2.pt(), track2.eta(), track2.phi(), massPi);
      pions.push_back(temp2);
      PionIndex.push_back(track2.globalIndex());
      PioncollIndex.push_back(track2.collisionId());
      PionSign.push_back(track2.sign());
      PionPx.push_back(track2.px());
      PionPy.push_back(track2.py());
      PionPz.push_back(track2.pz());
      PionP.push_back(totmompi);
    } // track loop ends

    if (!QA) {
      if (pions.size() != 0 && kaons.size() != 0) {
        for (auto ikaon = kaons.begin(); ikaon != kaons.end(); ++ikaon) {
          auto i1 = std::distance(kaons.begin(), ikaon);
          for (auto ipion = pions.begin(); ipion != pions.end(); ++ipion) {
            auto i3 = std::distance(pions.begin(), ipion);

            if (PionIndex.at(i3) <= KaonIndex.at(i1))
              continue;
            CKSVector = kaons.at(i1) + pions.at(i3);

            // openingangle = TMath::Abs((PionPx.at(i3) * KaonPx.at(i1) + PionPy.at(i3) * KaonPy.at(i1) + PionPz.at(i3) * KaonPz.at(i1)) / (PionP.at(i3) * KaonP.at(i1)));

            // openingangle = (PionPx.at(i3)*KaonPx.at(i1) + PionPy.at(i3)*KaonPy.at(i1) + PionPz.at(i3)*KaonPz.at(i1));
            // LOG(info) << "opening angle" << openingangle;

            if (TMath::Abs(CKSVector.Rapidity()) < 0.5) {
              if (PionSign.at(i3) * KaonSign.at(i1) < 0)
                histos.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, CKSVector.Pt(), CKSVector.M());
              else if (PionSign.at(i3) * KaonSign.at(i1) > 0)
                histos.fill(HIST("h3KstarInvMasslikeSign"), multiplicity, CKSVector.Pt(), CKSVector.M());
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(kstarqa, processSE, "Process Same event", true);

  void processME(EventCandidates const&, TrackCandidates const&)

  {

    for (auto& [c1, tracks1, c2, tracks2] : pair) {

      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }

      if (timFrameEvsel && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }

      if (piluprejection && (!c1.selection_bit(o2::aod::evsel::kNoSameBunchPileup) || !c2.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
        return;
      }
      if (goodzvertex && (!c1.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) || c2.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))) {
        return;
      }
      if (itstpctracks && (!c1.selection_bit(o2::aod::evsel::kIsVertexITSTPC) || c2.selection_bit(o2::aod::evsel::kIsVertexITSTPC))) {
        return;
      }

      // float multiplicity = 0.0f;
      /*      if (cfgMultFT0)
        multiplicity = c1.multZeqFT0A() + c1.multZeqFT0C();
        if (cfgMultFT0 == 0 && cfgCentFT0C == 1)
        multiplicity = c1.centFT0C();
        if (cfgMultFT0 == 0 && cfgCentFT0C == 0)*/
      auto multiplicity = c1.centFT0M();

      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!selectionTrack(t1))
          continue;
        if (!selectionTrack(t2))
          continue;
        if (!selectionPID(t1, 1))
          continue;
        if (!selectionPID(t2, 0))
          continue;
        if (MID) {
          if (MIDselectionPID(t1, 0)) // misidentified as pion
            continue;
          if (MIDselectionPID(t1, 2)) // misidentified as proton
            continue;
          if (MIDselectionPID(t2, 1)) // misidentified as kaon
            continue;
        }

        TLorentzVector KAON;
        KAON.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massKa);
        TLorentzVector PION;
        PION.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massPi);

        totmompimix = TMath::Sqrt(t2.px() * t2.px() + t2.py() * t2.py() + t2.pz() * t2.pz());
        totmomkamix = TMath::Sqrt(t1.px() * t1.px() + t1.py() * t1.py() + t1.pz() * t1.pz());

        // openinganglemix = TMath::Abs((t1.px() * t2.px() + t1.py() * t2.py() + t1.pz() * t2.pz()) / (totmomkamix * totmompimix));

        // LOG(info) << "mix angle" << openinganglemix;

        TLorentzVector CKSmix = KAON + PION;

        if (!QA) {
          if (TMath::Abs(CKSmix.Rapidity()) < 0.5) {
            if (t1.sign() * t2.sign() < 0)
              histos.fill(HIST("h3KstarInvMassMixed"), multiplicity, CKSmix.Pt(), CKSmix.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(kstarqa, processME, "Process Mixed event", true);

  void processGen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    histos.fill(HIST("events_check"), 0.5);
    if (std::abs(mcCollision.posZ()) < cutzvertex) {
      histos.fill(HIST("events_check"), 1.5);
    }
    // int Nchinel = 0;
    // for (auto& mcParticle : mcParticles) {
    //   auto pdgcode = std::abs(mcParticle.pdgCode());
    //   if (mcParticle.isPhysicalPrimary() && (pdgcode == 211 || pdgcode == 321 || pdgcode == 2212 || pdgcode == 11 || pdgcode == 13)) {
    //     if (std::abs(mcParticle.eta()) < 1.0) {
    //       Nchinel = Nchinel + 1;
    //     }
    //   }
    // }
    // if (Nchinel > 0 && std::abs(mcCollision.posZ()) < cutzvertex)
    //   histos.fill(HIST("events_check"), 2.5);
    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cutzvertex) {
        continue;
      }
      if (timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        continue;
      }

      if (piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        continue;
      }
      if (goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        continue;
      }
      if (itstpctracks && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        continue;
      }
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }

    // SelectedEvents.resize(nevts);
    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();
    histos.fill(HIST("events_check"), 3.5);
    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    histos.fill(HIST("events_check"), 4.5);
    for (auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) >= 0.5) {
        continue;
      }
      if (abs(mcParticle.pdgCode()) != 313) {
        continue;
      }
      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2) {
        continue;
      }
      auto passkaon = false;
      auto passpion = false;
      for (auto kCurrentDaughter : kDaughters) {
        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }
        if (abs(kCurrentDaughter.pdgCode()) == 321) {
          passkaon = true;
        } else if (abs(kCurrentDaughter.pdgCode()) == 211) {
          passpion = true;
        }
      }
      if (passkaon && passpion) {
        // if (mcParticle.pdgCode() > 0)
        histos.fill(HIST("k892Gen"), mcParticle.pt());
        // else
        //   histos.fill(HIST("k892GenAnti"), mcParticle.pt());
      }
    }
  }
  PROCESS_SWITCH(kstarqa, processGen, "Process Generated", false);

  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, aod::McCollisions const& /*mcCollisions*/)
  {

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

    if (!collision.has_mcCollision()) {
      return;
    }
    if (std::abs(collision.mcCollision().posZ()) > cutzvertex || !collision.sel8()) {
      return;
    }
    if (timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }

    if (piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    if (itstpctracks && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }

    // histos.fill(HIST("events_check"), 5.5);
    // auto oldindex = -999;
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!track1.has_mcParticle()) {
        continue;
      }
      auto track1ID = track1.index();
      for (auto track2 : tracks) {
        if (!track2.has_mcParticle()) {
          continue;
        }
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.index();
        if (track2ID == track1ID) {
          continue;
        }
        // if (!selectionPair(track1, track2)) {
        //   continue;
        // }
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());
        if (!mctrack1.isPhysicalPrimary()) {
          continue;
        }
        if (!mctrack2.isPhysicalPrimary()) {
          continue;
        }

        if (!(selectionPID(track1, 0) && selectionPID(track2, 1))) { // pion and kaon
          continue;
        }
        if (!(track1PDG == 211 && track2PDG == 321)) {
          continue;
        }
        for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
              continue;
            }
            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            if (std::abs(mothertrack1.y()) >= 0.5) {
              continue;
            }
            if (std::abs(mothertrack1.pdgCode()) != 313) {
              continue;
            }

            // if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
            //   histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
            //   continue;
            // }
            // oldindex = mothertrack1.globalIndex();
            pvec0 = array{track1.px(), track1.py(), track1.pz()};
            pvec1 = array{track2.px(), track2.py(), track2.pz()};
            auto arrMomrec = array{pvec0, pvec1};
            auto motherP = mothertrack1.p();
            auto motherE = mothertrack1.e();
            auto genMass = std::sqrt(motherE * motherE - motherP * motherP);
            auto recMass = RecoDecay::m(arrMomrec, array{massKa, massPi});
            // auto recpt = TMath::Sqrt((track1.px() + track2.px()) * (track1.px() + track2.px()) + (track1.py() + track2.py()) * (track1.py() + track2.py()));
            //// Resonance reconstruction
            lDecayDaughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), massPi);
            lDecayDaughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), massKa);
            lResonance = lDecayDaughter1 + lDecayDaughter2;
            histos.fill(HIST("h3KstarRec"), motherP);
            histos.fill(HIST("h1KstarRecMass"), recMass);
            // histos.fill(HIST("h1KstarRecpt"), recpt);
            histos.fill(HIST("h1genmass"), genMass);
            histos.fill(HIST("h1recpt"), lResonance.Pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(kstarqa, processRec, "Process Reconstructed", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstarqa>(cfgc)};
}
