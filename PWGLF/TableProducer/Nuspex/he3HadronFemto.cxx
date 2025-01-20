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
// Analysis task for he3-hadron femto analysis

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include <array>
#include <cstdlib>
#include <iterator> // std::prev

#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "ReconstructionDataFormats/Track.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFhe3HadronTables.h"
#include "PWGLF/Utils/svPoolCreator.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using CollBracket = o2::math_utils::Bracket<int>;

using McIter = aod::McParticles::iterator;
using CollBracket = o2::math_utils::Bracket<int>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::TOFSignal, aod::TOFEvTime>;
using TrackCandidatesMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::TOFSignal, aod::TOFEvTime, aod::McTrackLabels>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

// constexpr float he3Mass = o2::constants::physics::MassHelium3;
// constexpr float protonMass = o2::constants::physics::MassProton;
// constexpr float pionchargedMass = o2::constants::physics::MassPiPlus;
constexpr int li4PDG = 1000030040;
constexpr int prPDG = 2212;
constexpr int hePDG = 1000020030;
// constexpr int pichargedPDG = 211;

enum Selections {
  kNoCuts = 0,
  kTrackCuts,
  kPID,
  kAll
};

} // namespace

struct he3HadCandidate {

  float recoPtHe3() const { return signHe3 * std::hypot(momHe3[0], momHe3[1]); }
  float recoPhiHe3() const { return std::atan2(momHe3[1], momHe3[0]); }
  float recoEtaHe3() const { return std::asinh(momHe3[2] / std::abs(recoPtHe3())); }
  float recoPtHad() const { return signHad * std::hypot(momHad[0], momHad[1]); }
  float recoPhiHad() const { return std::atan2(momHad[1], momHad[0]); }
  float recoEtaHad() const { return std::asinh(momHad[2] / std::abs(recoPtHad())); }

  std::array<float, 3> momHe3 = {99.f, 99.f, 99.f};
  std::array<float, 3> momHad = {99.f, 99.f, 99.f};

  float signHe3 = 1.f;
  float signHad = 1.f;
  float invMass = -10.f;
  float DCAxyHe3 = -10.f;
  float DCAzHe3 = -10.f;
  float DCAxyHad = -10.f;
  float DCAzHad = -10.f;

  uint16_t tpcSignalHe3 = 0u;
  uint16_t tpcSignalHad = 0u;
  float momHe3TPC = -99.f;
  float momHadTPC = -99.f;
  uint8_t nTPCClustersHe3 = 0u;
  uint8_t sharedClustersHe3 = 0u;
  uint8_t sharedClustersHad = 0u;
  float chi2TPCHe3 = -10.f;
  float chi2TPCHad = -10.f;
  float nSigmaHe3 = -10.f;
  float nSigmaHad = -10.f;
  uint32_t PIDtrkHe3 = 0xFFFFF; // PID in tracking
  uint32_t PIDtrkHad = 0xFFFFF;
  float massTOFHe3 = -10;
  float massTOFHad = -10;
  uint32_t itsClSizeHe3 = 0u;
  uint32_t itsClSizeHad = 0u;

  bool isBkgUS = false; // unlike sign
  bool isBkgEM = false; // event mixing

  int trackIDHe3 = -1;
  int trackIDHad = -1;

  float l4MassMC = -10.f;
  float l4PtMC = -99.f;
  float momHe3MC = -99.f;
  float etaHe3MC = -99.f;
  float phiHe3MC = -99.f;
  float momHadMC = -99.f;
  float etaHadMC = -99.f;
  float phiHadMC = -99.f;

  // collision information
  int32_t collisionID = 0;
};

struct he3hadronfemto {

  Produces<aod::he3HadronTable> m_outputDataTable;
  Produces<aod::he3HadronTableMC> m_outputMCTable;
  Produces<aod::he3HadronMult> m_outputMultiplicityTable;

  // Selections
  Configurable<int> setting_HadPDGCode{"setting_HadPDGCode", 211, "Hadron - PDG code"};
  Configurable<float> setting_cutVertex{"setting_cutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> setting_cutRigidityMinHe3{"setting_cutRigidityMinHe3", 0.8f, "Minimum rigidity for He3"};
  Configurable<float> setting_cutEta{"setting_cutEta", 0.9f, "Eta cut on daughter track"};
  Configurable<float> setting_cutDCAxy{"setting_cutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> setting_cutDCAz{"setting_cutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> setting_cutChi2tpcLow{"setting_cutChi2tpcLow", 0.5f, "Low cut on TPC chi2"};
  Configurable<float> setting_cutInvMass{"setting_cutInvMass", 0.0f, "Invariant mass upper limit"};
  Configurable<float> setting_cutPtMinhe3Had{"setting_cutPtMinhe3Had", 0.0f, "Minimum PT cut on he3Had4"};
  Configurable<float> setting_cutClSizeItsHe3{"setting_cutClSizeItsHe3", 4.0f, "Minimum ITS cluster size for He3"};
  Configurable<float> setting_cutNsigmaTPC{"setting_cutNsigmaTPC", 3.0f, "Value of the TPC Nsigma cut"};
  Configurable<float> setting_cutNsigmaITS{"setting_cutNsigmaITS", -1.5f, "Value of the TPC Nsigma cut"};
  Configurable<float> setting_cutPtMinTOFHad{"setting_cutPtMinTOFHad", 0.4f, "Minimum pT to apply the TOF cut on hadrons"};
  Configurable<float> setting_cutNsigmaTOF{"setting_cutNsigmaTOF", 3.0f, "Value of the TOF Nsigma cut"};
  Configurable<int> setting_noMixedEvents{"setting_noMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> setting_enableBkgUS{"setting_enableBkgUS", false, "Enable US background"};
  Configurable<bool> setting_saveUSandLS{"setting_saveUSandLS", true, "Save All Pairs"};
  Configurable<bool> setting_isMC{"setting_isMC", false, "Run MC"};
  Configurable<bool> setting_fillMultiplicity{"setting_fillMultiplicity", false, "Fill multiplicity table"};

  // Zorro
  Configurable<bool> setting_skimmedProcessing{"setting_skimmedProcessing", false, "Skimmed dataset processing"};

  // svPool
  Configurable<bool> setting_skipAmbiTracks{"setting_skipAmbiTracks", false, "Skip ambiguous tracks"};
  Configurable<float> setting_customVertexerTimeMargin{"setting_customVertexerTimeMargin", 800, "Time margin for custom vertexer (ns)"};

  // CCDB options
  Configurable<double> setting_d_bz_input{"setting_d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> setting_ccdburl{"setting_ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> setting_grpPath{"setting_grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> setting_grpmagPath{"setting_grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> setting_lutPath{"setting_lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> setting_geoPath{"setting_geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> setting_pidPath{"setting_pidPath", "", "Path to the PID response object"};

  Configurable<LabeledArray<double>> setting_BetheBlochParams{"setting_BetheBlochParams", {betheBlochDefault[0], 1, 6, {"He3"}, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  Configurable<bool> setting_compensatePIDinTracking{"setting_compensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};
  Configurable<int> setting_materialCorrection{"setting_materialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Material correction type"};

  Preslice<TrackCandidates> m_perCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> m_perColMC = aod::track::collisionId;

  // binning for EM background
  ConfigurableAxis axisVertex{"axisVertex", {30, -10, 10}, "Binning for multiplicity"};
  ConfigurableAxis axisCentrality{"axisCentrality", {40, 0, 100}, "Binning for centrality"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType binningPolicy{{axisVertex, axisCentrality}, true};
  SliceCache cache;
  SameKindPair<CollisionsFull, TrackCandidates, BinningType> m_pair{binningPolicy, setting_noMixedEvents, -1, &cache};

  std::array<float, 6> m_BBparamsHe;

  std::vector<int> m_recoCollisionIDs;
  std::vector<bool> m_goodCollisions;
  std::vector<SVCand> m_trackPairs;
  o2::vertexing::DCAFitterN<2> m_fitter;
  svPoolCreator m_svPoolCreator{hePDG, prPDG};

  int m_runNumber;
  float m_d_bz;
  Service<o2::ccdb::BasicCCDBManager> m_ccdb;
  Zorro m_zorro;
  OutputObj<ZorroSummary> m_zorroSummary{"zorroSummary"};

  HistogramRegistry m_qaRegistry{
    "QA",
    {
      {"hVtxZ", "Vertex distribution in Z;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
      {"hNcontributor", "Number of primary vertex contributor", {HistType::kTH1F, {{2000, 0.0f, 2000.0f}}}},
      {"hTrackSel", "Accepted tracks", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
      {"hEvents", "; Events;", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
      {"hEmptyPool", "svPoolCreator did not find track pairs false/true", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"hDCAxyHe3", ";DCA_{xy} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
      {"hDCAzHe3", ";DCA_{z} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
      {"hhe3HadtInvMass", "; M(^{3}He + p) (GeV/#it{c}^{2})", {HistType::kTH1F, {{300, 3.74f, 4.34f}}}},
      {"hHe3Pt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hHadronPt", "Pt distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},
      {"h2dEdxHe3candidates", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
      {"h2NsigmaHe3TPC", "NsigmaHe3 TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(^{3}He)", {HistType::kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
      {"h2NsigmaHe3TPC_preselection", "NsigmaHe3 TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(^{3}He)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
      {"h2NSigmaHe3ITS_preselection", "NsigmaHe3 ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} ^{3}He", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {120, -3.0f, 3.0f}}}},
      {"h2NSigmaHe3ITS", "NsigmaHe3 ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} ^{3}He", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {120, -3.0f, 3.0f}}}},
      {"h2NsigmaHadronTPC", "NsigmaHadron TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
      {"h2NsigmaHadronTPC_preselection", "NsigmaHe3 TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(^{3}He)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
      {"h2NsigmaHadronTOF", "NsigmaHadron TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
      {"h2NsigmaHadronTOF_preselection", "NsigmaHadron TOF distribution; #iit{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  void init(o2::framework::InitContext&)
  {
    m_zorroSummary.setObject(m_zorro.getZorroSummary());
    m_runNumber = 0;

    m_ccdb->setURL(setting_ccdburl);
    m_ccdb->setCaching(true);
    m_ccdb->setLocalObjectValidityChecking();
    m_ccdb->setFatalWhenNull(false);

    m_fitter.setPropagateToPCA(true);
    m_fitter.setMaxR(200.);
    m_fitter.setMinParamChange(1e-3);
    m_fitter.setMinRelChi2Change(0.9);
    m_fitter.setMaxDZIni(1e9);
    m_fitter.setMaxChi2(1e9);
    m_fitter.setUseAbsDCA(true);
    int mat{static_cast<int>(setting_materialCorrection)};
    m_fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

    m_svPoolCreator.setTimeMargin(setting_customVertexerTimeMargin);
    if (setting_skipAmbiTracks) {
      m_svPoolCreator.setSkipAmbiTracks();
    }

    for (int i = 0; i < 5; i++) {
      m_BBparamsHe[i] = setting_BetheBlochParams->get("He3", Form("p%i", i));
    }
    m_BBparamsHe[5] = setting_BetheBlochParams->get("He3", "resolution");

    std::vector<std::string> selection_labels = {"All", "Track selection", "PID"};
    for (int i = 0; i < Selections::kAll; i++) {
      m_qaRegistry.get<TH1>(HIST("hTrackSel"))->GetXaxis()->SetBinLabel(i + 1, selection_labels[i].c_str());
    }

    std::vector<std::string> events_labels = {"All", "Selected", "Zorro He events"};
    for (int i = 0; i < Selections::kAll; i++) {
      m_qaRegistry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(i + 1, events_labels[i].c_str());
    }

    m_qaRegistry.get<TH1>(HIST("hEmptyPool"))->GetXaxis()->SetBinLabel(1, "False");
    m_qaRegistry.get<TH1>(HIST("hEmptyPool"))->GetXaxis()->SetBinLabel(2, "True");
  }

  void initCCDB(const aod::BCsWithTimestamps::iterator& bc)
  {
    if (m_runNumber == bc.runNumber()) {
      return;
    }
    if (setting_skimmedProcessing) {
      m_zorro.initCCDB(m_ccdb.service, bc.runNumber(), bc.timestamp(), "fHe");
      m_zorro.populateHistRegistry(m_qaRegistry, bc.runNumber());
    }
    m_runNumber = bc.runNumber();

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = m_ccdb->getForTimeStamp<o2::parameters::GRPObject>(setting_grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (setting_d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        m_d_bz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << m_d_bz << " kZG";
      } else {
        m_d_bz = setting_d_bz_input;
      }
    } else {
      grpmag = m_ccdb->getForTimeStamp<o2::parameters::GRPMagField>(setting_grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << setting_grpmagPath << " of object GRPMagField and " << setting_grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (setting_d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        m_d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << m_d_bz << " kZG";
      } else {
        m_d_bz = setting_d_bz_input;
      }
    }
  }

  // ==================================================================================================================

  template <bool isMC, typename Tcollision>
  bool selectCollision(const Tcollision& collision, const aod::BCsWithTimestamps&)
  {
    m_qaRegistry.fill(HIST("hEvents"), 0);

    if constexpr (isMC) {
      if (/*!collision.sel8() ||*/ std::abs(collision.posZ()) > setting_cutVertex) {
        return false;
      }
    } else {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8() || std::abs(collision.posZ()) > setting_cutVertex) {
        return false;
      }
      if (setting_skimmedProcessing) {
        bool zorroSelected = m_zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC());
        if (zorroSelected) {
          m_qaRegistry.fill(HIST("hEvents"), 2);
        }
      }
    }

    m_qaRegistry.fill(HIST("hEvents"), 1);
    m_qaRegistry.fill(HIST("hNcontributor"), collision.numContrib());
    m_qaRegistry.fill(HIST("hVtxZ"), collision.posZ());
    return true;
  }

  template <typename Ttrack>
  bool selectTrack(const Ttrack& candidate)
  {
    if (std::abs(candidate.eta()) > setting_cutEta) {
      return false;
    }
    if (candidate.itsNCls() < 5 ||
        candidate.tpcNClsFound() < 90 ||
        candidate.tpcNClsCrossedRows() < 70 ||
        candidate.tpcNClsCrossedRows() < 0.8 * candidate.tpcNClsFindable() ||
        candidate.tpcChi2NCl() > 4.f ||
        candidate.tpcChi2NCl() < setting_cutChi2tpcLow ||
        candidate.itsChi2NCl() > 36.f) {
      return false;
    }

    return true;
  }

  template <typename Ttrack>
  float computeTPCNSigmaHadron(const Ttrack& candidate)
  {
    float tpcNSigmaHad = 0;
    if (setting_HadPDGCode == 211) {
      tpcNSigmaHad = candidate.tpcNSigmaPi();
    } else if (setting_HadPDGCode == 2212) {
      tpcNSigmaHad = candidate.tpcNSigmaPr();
    } else {
      LOG(info) << "invalid PDG code for TPC";
    }
    return tpcNSigmaHad;
  }

  template <typename Ttrack>
  float computeTOFNSigmaHadron(const Ttrack& candidate)
  {
    float tofNSigmaHad = 0;
    if (setting_HadPDGCode == 211) {
      tofNSigmaHad = candidate.tofNSigmaPi();
    } else if (setting_HadPDGCode == 2212) {
      tofNSigmaHad = candidate.tofNSigmaPr();
    } else {
      LOG(info) << "invalid PDG code for TOF";
    }
    return tofNSigmaHad;
  }

  template <typename Ttrack>
  bool selectionPIDHadron(const Ttrack& candidate)
  {
    auto tpcNSigmaHad = computeTPCNSigmaHadron(candidate);
    m_qaRegistry.fill(HIST("h2NsigmaHadronTPC_preselection"), candidate.tpcInnerParam(), tpcNSigmaHad);
    if (candidate.hasTOF() && candidate.pt() > setting_cutPtMinTOFHad) {
      auto tofNSigmaHad = computeTOFNSigmaHadron(candidate);

      if (std::abs(tpcNSigmaHad) > setting_cutNsigmaTPC) {
        return false;
      }
      m_qaRegistry.fill(HIST("h2NsigmaHadronTOF_preselection"), candidate.pt(), tofNSigmaHad);
      if (std::abs(tofNSigmaHad) > setting_cutNsigmaTOF) {
        return false;
      }
      m_qaRegistry.fill(HIST("h2NsigmaHadronTPC"), candidate.pt(), tpcNSigmaHad);
      m_qaRegistry.fill(HIST("h2NsigmaHadronTOF"), candidate.pt(), tofNSigmaHad);
      return true;
    } else if (std::abs(tpcNSigmaHad) < setting_cutNsigmaTPC) {
      m_qaRegistry.fill(HIST("h2NsigmaHadronTPC"), candidate.pt(), tpcNSigmaHad);
      return true;
    }
    return false;
  }

  template <typename Ttrack>
  float computeNSigmaHe3(const Ttrack& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParam = (heliumPID && setting_compensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), m_BBparamsHe[0], m_BBparamsHe[1], m_BBparamsHe[2], m_BBparamsHe[3], m_BBparamsHe[4]);

    double resoTPC{expTPCSignal * m_BBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename Ttrack>
  bool selectionPIDHe3(const Ttrack& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParam = (heliumPID && setting_compensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();

    if (correctedTPCinnerParam < setting_cutRigidityMinHe3) {
      return false;
    }

    auto nSigmaHe3 = computeNSigmaHe3(candidate);
    m_qaRegistry.fill(HIST("h2NsigmaHe3TPC_preselection"), candidate.sign() * 2 * candidate.pt(), nSigmaHe3);
    if (std::abs(nSigmaHe3) > setting_cutNsigmaTPC) {
      return false;
    }
    //
    o2::aod::ITSResponse m_responseITS;
    auto ITSnSigmaHe3 = m_responseITS.nSigmaITS<o2::track::PID::Helium3>(candidate.itsClusterSizes(), 2 * candidate.p(), candidate.eta());
    //
    m_qaRegistry.fill(HIST("h2NSigmaHe3ITS_preselection"), candidate.sign() * 2 * candidate.pt(), ITSnSigmaHe3);
    if (ITSnSigmaHe3 < setting_cutNsigmaITS) {
      return false;
    }

    m_qaRegistry.fill(HIST("h2dEdxHe3candidates"), candidate.sign() * correctedTPCinnerParam, candidate.tpcSignal());
    m_qaRegistry.fill(HIST("h2NsigmaHe3TPC"), candidate.sign() * 2 * candidate.pt(), nSigmaHe3);
    m_qaRegistry.fill(HIST("h2NSigmaHe3ITS"), candidate.sign() * 2 * candidate.pt(), ITSnSigmaHe3);
    return true;
  }

  // ==================================================================================================================

  template <typename Ttrack, typename Tcollisions, typename Ttracks>
  bool fillCandidateInfo(const Ttrack& trackHe3, const Ttrack& trackHad, const CollBracket& collBracket, const Tcollisions& collisions, he3HadCandidate& he3Hadcand, const Ttracks& /*trackTable*/, bool isMixedEvent)
  {
    if (!isMixedEvent) {
      auto trackCovHe3 = getTrackParCov(trackHe3);
      auto trackCovHad = getTrackParCov(trackHad);
      int nCand = 0;
      try {
        nCand = m_fitter.process(trackCovHe3, trackCovHad);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        return false;
      }
      if (nCand == 0) {
        return false;
      }

      // associate collision id as the one that minimises the distance between the vertex and the PCAs of the daughters
      double distanceMin = -1;
      unsigned int collIdxMin = 0;

      for (int collIdx = collBracket.getMin(); collIdx <= collBracket.getMax(); collIdx++) {
        auto collision = collisions.rawIteratorAt(collIdx);
        std::array<float, 3> collVtx = {collision.posX(), collision.posY(), collision.posZ()};
        const auto& PCA = m_fitter.getPCACandidate();
        float distance = 0;
        for (int i = 0; i < 3; i++) {
          distance += (PCA[i] - collVtx[i]) * (PCA[i] - collVtx[i]);
        }
        if (distanceMin < 0 || distance < distanceMin) {
          distanceMin = distance;
          collIdxMin = collIdx;
        }
      }

      if (!m_goodCollisions[collIdxMin]) {
        return false;
      }
      he3Hadcand.collisionID = collIdxMin;
    } else {
      he3Hadcand.collisionID = collBracket.getMin();
    }

    he3Hadcand.momHe3 = std::array{trackHe3.px(), trackHe3.py(), trackHe3.pz()};
    for (int i = 0; i < 3; i++)
      he3Hadcand.momHe3[i] = he3Hadcand.momHe3[i] * 2;
    he3Hadcand.momHad = std::array{trackHad.px(), trackHad.py(), trackHad.pz()};
    float invMass = 0;
    if (setting_HadPDGCode == 211) {
      invMass = RecoDecay::m(std::array{he3Hadcand.momHe3, he3Hadcand.momHad}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassPiPlus});
    } else if (setting_HadPDGCode == 2212) {
      invMass = RecoDecay::m(std::array{he3Hadcand.momHe3, he3Hadcand.momHad}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassProton});
    } else {
      LOG(info) << "invalid PDG code for invMass";
    }
    // float invMass = RecoDecay::m(std::array{he3Hadcand.momHe3, he3Hadcand.momHad}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassPiPlus});
    if (setting_cutInvMass > 0 && invMass > setting_cutInvMass) {
      return false;
    }
    float pthe3Had = std::hypot(he3Hadcand.momHe3[0] + he3Hadcand.momHad[0], he3Hadcand.momHe3[1] + he3Hadcand.momHad[1]);
    if (pthe3Had < setting_cutPtMinhe3Had) {
      return false;
    }

    he3Hadcand.signHe3 = trackHe3.sign();
    he3Hadcand.signHad = trackHad.sign();

    he3Hadcand.DCAxyHe3 = trackHe3.dcaXY();
    he3Hadcand.DCAzHe3 = trackHe3.dcaZ();
    he3Hadcand.DCAxyHad = trackHad.dcaXY();
    he3Hadcand.DCAzHad = trackHad.dcaZ();

    he3Hadcand.tpcSignalHe3 = trackHe3.tpcSignal();
    bool heliumPID = trackHe3.pidForTracking() == o2::track::PID::Helium3 || trackHe3.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParamHe3 = (heliumPID && setting_compensatePIDinTracking) ? trackHe3.tpcInnerParam() / 2.f : trackHe3.tpcInnerParam();
    he3Hadcand.momHe3TPC = correctedTPCinnerParamHe3;
    he3Hadcand.tpcSignalHad = trackHad.tpcSignal();
    he3Hadcand.momHadTPC = trackHad.tpcInnerParam();

    he3Hadcand.nTPCClustersHe3 = trackHe3.tpcNClsFound();
    he3Hadcand.nSigmaHe3 = computeNSigmaHe3(trackHe3);
    he3Hadcand.nSigmaHad = computeTPCNSigmaHadron(trackHad);

    he3Hadcand.chi2TPCHe3 = trackHe3.tpcChi2NCl();
    he3Hadcand.chi2TPCHad = trackHad.tpcChi2NCl();

    he3Hadcand.PIDtrkHe3 = trackHe3.pidForTracking();
    he3Hadcand.PIDtrkHad = trackHad.pidForTracking();

    he3Hadcand.itsClSizeHe3 = trackHe3.itsClusterSizes();
    he3Hadcand.itsClSizeHad = trackHad.itsClusterSizes();

    he3Hadcand.sharedClustersHe3 = trackHe3.tpcNClsShared();
    he3Hadcand.sharedClustersHad = trackHad.tpcNClsShared();

    he3Hadcand.isBkgUS = trackHe3.sign() * trackHad.sign() < 0;
    he3Hadcand.isBkgEM = isMixedEvent;

    he3Hadcand.invMass = invMass;

    he3Hadcand.trackIDHe3 = trackHe3.globalIndex();
    he3Hadcand.trackIDHad = trackHad.globalIndex();

    if (trackHe3.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackHe3);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      bool heliumPID = trackHe3.pidForTracking() == o2::track::PID::Helium3 || trackHe3.pidForTracking() == o2::track::PID::Alpha;
      float correctedTPCinnerParamHe3 = (heliumPID && setting_compensatePIDinTracking) ? trackHe3.tpcInnerParam() / 2.f : trackHe3.tpcInnerParam();
      he3Hadcand.massTOFHe3 = correctedTPCinnerParamHe3 * 2.f * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    if (trackHad.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackHad);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      he3Hadcand.massTOFHad = trackHad.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
    }

    return true;
  }

  template <typename Mc>
  void fillCandidateInfoMC(const Mc& mctrackHe3, const Mc& mctrackHad, const Mc& mctrackMother, he3HadCandidate& he3Hadcand)
  {
    he3Hadcand.momHe3MC = mctrackHe3.pt() * (mctrackHe3.pdgCode() > 0 ? 1 : -1);
    he3Hadcand.etaHe3MC = mctrackHe3.eta();
    he3Hadcand.phiHe3MC = mctrackHe3.phi();
    he3Hadcand.momHadMC = mctrackHad.pt() * (mctrackHad.pdgCode() > 0 ? 1 : -1);
    he3Hadcand.etaHadMC = mctrackHad.eta();
    he3Hadcand.phiHadMC = mctrackHad.phi();
    he3Hadcand.l4PtMC = mctrackMother.pt() * (mctrackMother.pdgCode() > 0 ? 1 : -1);
    const double eLit = mctrackHe3.e() + mctrackHad.e();
    he3Hadcand.l4MassMC = std::sqrt(eLit * eLit - mctrackMother.p() * mctrackMother.p());
  }

  template <typename Ttrack>
  void pairTracksSameEvent(const Ttrack& tracks)
  {
    for (auto track0 : tracks) {

      m_qaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);

      if (!selectTrack(track0)) {
        continue;
      }
      m_qaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);

      if (!selectionPIDHe3(track0)) {
        continue;
      }
      m_qaRegistry.fill(HIST("hTrackSel"), Selections::kPID);

      for (auto track1 : tracks) {

        if (track0 == track1) {
          continue;
        }

        if (!setting_saveUSandLS) {
          if (!setting_enableBkgUS && (track0.sign() * track1.sign() < 0)) {
            continue;
          }
          if (setting_enableBkgUS && (track0.sign() * track1.sign() > 0)) {
            continue;
          }
        }

        if (!selectTrack(track1) || !selectionPIDHadron(track1)) {
          continue;
        }

        SVCand trackPair;
        trackPair.tr0Idx = track0.globalIndex();
        trackPair.tr1Idx = track1.globalIndex();
        const int collIdx = track0.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        trackPair.collBracket = collBracket;
        m_trackPairs.push_back(trackPair);
      }
    }
  }

  template <typename T>
  void pairTracksEventMixing(T& he3Cands, T& hadronCands)
  {
    for (auto& he3Cand : he3Cands) {
      if (!selectTrack(he3Cand) || !selectionPIDHe3(he3Cand)) {
        continue;
      }
      for (auto& hadronCand : hadronCands) {
        if (!selectTrack(hadronCand) || !selectionPIDHadron(hadronCand)) {
          continue;
        }

        SVCand trackPair;
        trackPair.tr0Idx = he3Cand.globalIndex();
        trackPair.tr1Idx = hadronCand.globalIndex();
        const int collIdx = he3Cand.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        trackPair.collBracket = collBracket;
        m_trackPairs.push_back(trackPair);
      }
    }
  }

  template <typename Tcoll>
  void fillTable(const he3HadCandidate& he3Hadcand, const Tcoll& collision, bool isMC = false)
  {
    m_outputDataTable(
      he3Hadcand.recoPtHe3(),
      he3Hadcand.recoEtaHe3(),
      he3Hadcand.recoPhiHe3(),
      he3Hadcand.recoPtHad(),
      he3Hadcand.recoEtaHad(),
      he3Hadcand.recoPhiHad(),
      he3Hadcand.DCAxyHe3,
      he3Hadcand.DCAzHe3,
      he3Hadcand.DCAxyHad,
      he3Hadcand.DCAzHad,
      he3Hadcand.tpcSignalHe3,
      he3Hadcand.momHe3TPC,
      he3Hadcand.tpcSignalHad,
      he3Hadcand.momHadTPC,
      he3Hadcand.nTPCClustersHe3,
      he3Hadcand.nSigmaHe3,
      he3Hadcand.nSigmaHad,
      he3Hadcand.chi2TPCHe3,
      he3Hadcand.chi2TPCHad,
      he3Hadcand.massTOFHe3,
      he3Hadcand.massTOFHad,
      he3Hadcand.PIDtrkHe3,
      he3Hadcand.PIDtrkHad,
      he3Hadcand.itsClSizeHe3,
      he3Hadcand.itsClSizeHad,
      he3Hadcand.sharedClustersHe3,
      he3Hadcand.sharedClustersHad,
      he3Hadcand.isBkgUS,
      he3Hadcand.isBkgEM);
    if (isMC) {
      m_outputMCTable(
        he3Hadcand.momHe3MC,
        he3Hadcand.etaHe3MC,
        he3Hadcand.phiHe3MC,
        he3Hadcand.momHadMC,
        he3Hadcand.etaHadMC,
        he3Hadcand.phiHadMC,
        he3Hadcand.l4PtMC,
        he3Hadcand.l4MassMC);
    }
    if (setting_fillMultiplicity) {
      m_outputMultiplicityTable(
        collision.globalIndex(),
        collision.posZ(),
        collision.numContrib(),
        collision.centFT0C(),
        collision.multFT0C());
    }
  }

  void fillHistograms(const he3HadCandidate& he3Hadcand)
  {
    m_qaRegistry.fill(HIST("hHe3Pt"), he3Hadcand.recoPtHe3());
    m_qaRegistry.fill(HIST("hHadronPt"), he3Hadcand.recoPtHad());
    m_qaRegistry.fill(HIST("hhe3HadtInvMass"), he3Hadcand.invMass);
    m_qaRegistry.fill(HIST("hDCAxyHe3"), he3Hadcand.DCAxyHe3);
    m_qaRegistry.fill(HIST("hDCAzHe3"), he3Hadcand.DCAzHe3);
  }

  // ==================================================================================================================

  template <typename Tcollisions, typename Ttracks>
  void fillPairs(const Tcollisions& collisions, const Ttracks& tracks, const bool isMixedEvent)
  {
    for (auto& trackPair : m_trackPairs) {

      auto heTrack = tracks.rawIteratorAt(trackPair.tr0Idx);
      auto hadTrack = tracks.rawIteratorAt(trackPair.tr1Idx);
      auto collBracket = trackPair.collBracket;

      he3HadCandidate he3Hadcand;
      if (!fillCandidateInfo(heTrack, hadTrack, collBracket, collisions, he3Hadcand, tracks, isMixedEvent)) {
        continue;
      }
      fillHistograms(he3Hadcand);
      auto collision = collisions.rawIteratorAt(he3Hadcand.collisionID);
      fillTable(he3Hadcand, collision, /*isMC*/ false);
    }
  }

  template <typename Tcollisions, typename TmcParticles>
  void fillMcParticles(const Tcollisions& collisions, const TmcParticles& mcParticles, std::vector<unsigned int>& filledMothers)
  {
    for (auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.pdgCode()) != li4PDG || std::abs(mcParticle.y()) > 1 || mcParticle.isPhysicalPrimary() == false) {
        continue;
      }

      if (std::find(filledMothers.begin(), filledMothers.end(), mcParticle.globalIndex()) != filledMothers.end()) {
        continue;
      }

      auto kDaughters = mcParticle.template daughters_as<aod::McParticles>();
      bool daughtHe3(false), daughtHad(false);
      McIter mcHe3, mcHad;
      for (auto kCurrentDaughter : kDaughters) {
        if (std::abs(kCurrentDaughter.pdgCode()) == hePDG) {
          daughtHe3 = true;
          mcHe3 = kCurrentDaughter;
        } else if (std::abs(kCurrentDaughter.pdgCode()) == prPDG) {
          daughtHad = true;
          mcHad = kCurrentDaughter;
        }
      }
      if (daughtHe3 && daughtHad) {
        he3HadCandidate he3Hadcand;
        fillCandidateInfoMC(mcHe3, mcHad, mcParticle, he3Hadcand);
        auto collision = collisions.rawIteratorAt(he3Hadcand.collisionID);
        fillTable(he3Hadcand, collision, /*isMC*/ true);
      }
    }
  }

  // ==================================================================================================================

  void processSameEvent(const CollisionsFull& collisions, const TrackCandidates& tracks, const aod::BCsWithTimestamps& bcs)
  {
    m_goodCollisions.clear();
    m_goodCollisions.resize(collisions.size(), false);

    for (auto& collision : collisions) {

      m_trackPairs.clear();

      if (!selectCollision</*isMC*/ false>(collision, bcs)) {
        continue;
      }

      m_goodCollisions[collision.globalIndex()] = true;
      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(m_perCol, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      pairTracksSameEvent(TrackTable_thisCollision);

      if (m_trackPairs.size() == 0) {
        continue;
      }

      fillPairs(collisions, tracks, /*isMixedEvent*/ false);
    }
  }
  PROCESS_SWITCH(he3hadronfemto, processSameEvent, "Process Same event", false);

  void processMixedEvent(const CollisionsFull& collisions, const TrackCandidates& tracks)
  {
    LOG(debug) << "Processing mixed event";
    m_trackPairs.clear();

    for (auto& [c1, tracks1, c2, tracks2] : m_pair) {
      if (!c1.sel8() || !c2.sel8()) {
        continue;
      }

      m_qaRegistry.fill(HIST("hNcontributor"), c1.numContrib());
      m_qaRegistry.fill(HIST("hVtxZ"), c1.posZ());

      pairTracksEventMixing(tracks1, tracks2);
      pairTracksEventMixing(tracks2, tracks1);
    }

    fillPairs(collisions, tracks, /*isMixedEvent*/ true);
  }
  PROCESS_SWITCH(he3hadronfemto, processMixedEvent, "Process Mixed event", false);

  void processMC(const CollisionsFullMC& collisions, const aod::BCsWithTimestamps& bcs, const TrackCandidatesMC& tracks, const aod::McParticles& mcParticles)
  {
    std::vector<unsigned int> filledMothers;

    m_goodCollisions.clear();
    m_goodCollisions.resize(collisions.size(), false);

    for (auto& collision : collisions) {

      m_trackPairs.clear();

      if (!selectCollision</*isMC*/ true>(collision, bcs)) {
        continue;
      }

      const uint64_t collIdx = collision.globalIndex();
      m_goodCollisions[collIdx] = true;
      auto TrackTable_thisCollision = tracks.sliceBy(m_perColMC, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      pairTracksSameEvent(TrackTable_thisCollision);

      for (auto& trackPair : m_trackPairs) {

        auto heTrack = tracks.rawIteratorAt(trackPair.tr0Idx);
        auto prTrack = tracks.rawIteratorAt(trackPair.tr1Idx);
        auto collBracket = trackPair.collBracket;

        if (!heTrack.has_mcParticle() || !prTrack.has_mcParticle()) {
          continue;
        }

        auto mctrackHe3 = heTrack.mcParticle();
        auto mctrackHad = prTrack.mcParticle();

        if (std::abs(mctrackHe3.pdgCode()) != hePDG || std::abs(mctrackHad.pdgCode()) != prPDG) {
          continue;
        }

        for (auto& mothertrack : mctrackHe3.mothers_as<aod::McParticles>()) {
          for (auto& mothertrackHad : mctrackHad.mothers_as<aod::McParticles>()) {

            if (mothertrack != mothertrackHad || std::abs(mothertrack.pdgCode()) != li4PDG || std::abs(mothertrack.y()) > 1) {
              continue;
            }

            he3HadCandidate he3Hadcand;
            if (!fillCandidateInfo(heTrack, prTrack, collBracket, collisions, he3Hadcand, tracks, /*mix*/ false)) {
              continue;
            }
            fillCandidateInfoMC(mctrackHe3, mctrackHad, mothertrack, he3Hadcand);
            fillHistograms(he3Hadcand);
            auto collision = collisions.rawIteratorAt(he3Hadcand.collisionID);
            fillTable(he3Hadcand, collision, /*isMC*/ true);
            filledMothers.push_back(mothertrack.globalIndex());
          }
        }
      }
    }

    fillMcParticles(collisions, mcParticles, filledMothers);
  }
  PROCESS_SWITCH(he3hadronfemto, processMC, "Process MC", false);

  void processSameEventPools(const CollisionsFull& collisions, const TrackCandidates& tracks, const aod::AmbiguousTracks& ambiguousTracks, const aod::BCsWithTimestamps& bcs)
  {
    m_goodCollisions.clear();
    m_goodCollisions.resize(collisions.size(), false);

    for (auto& collision : collisions) {
      if (selectCollision</*isMC*/ false>(collision, bcs)) {
        m_goodCollisions[collision.globalIndex()] = true;
      }
    }

    m_svPoolCreator.clearPools();
    m_svPoolCreator.fillBC2Coll(collisions, bcs);

    for (auto& track : tracks) {

      m_qaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);
      if (!selectTrack(track))
        continue;
      m_qaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);

      bool selHad = selectionPIDHadron(track);
      bool selHe = selectionPIDHe3(track);
      if ((!selHad && !selHe) || (selHad && selHe)) {
        continue;
      }
      m_qaRegistry.fill(HIST("hTrackSel"), Selections::kPID);

      int pdgHypo = selHe ? hePDG : prPDG;

      m_svPoolCreator.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }

    m_trackPairs = m_svPoolCreator.getSVCandPool(collisions, true);
    if (m_trackPairs.size() == 0) {
      m_qaRegistry.fill(HIST("hEmptyPool"), 1);
      return;
    }
    m_qaRegistry.fill(HIST("hEmptyPool"), 0);

    fillPairs(collisions, tracks, /*isMixedEvent*/ false);
  }
  PROCESS_SWITCH(he3hadronfemto, processSameEventPools, "Process Same event pools", false);

  void processMcPools(const CollisionsFullMC& collisions, const TrackCandidatesMC& tracks, const aod::AmbiguousTracks& ambiguousTracks, const aod::BCsWithTimestamps& bcs, const aod::McParticles& mcParticles, const aod::McTrackLabels& mcTrackLabels)
  {
    std::vector<unsigned int> filledMothers;

    m_goodCollisions.clear();
    m_goodCollisions.resize(collisions.size(), false);

    for (auto& collision : collisions) {
      if (selectCollision</*isMC*/ true>(collision, bcs)) {
        m_goodCollisions[collision.globalIndex()] = true;
      }
    }

    m_svPoolCreator.clearPools();
    m_svPoolCreator.fillBC2Coll(collisions, bcs);

    for (auto& track : tracks) {

      m_qaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);
      if (!selectTrack(track))
        continue;
      m_qaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);

      bool selHad = selectionPIDHadron(track);
      bool selHe = selectionPIDHe3(track);
      if ((!selHad && !selHe) || (selHad && selHe))
        continue;
      m_qaRegistry.fill(HIST("hTrackSel"), Selections::kPID);

      int pdgHypo = selHe ? hePDG : prPDG;

      m_svPoolCreator.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }

    auto& svPool = m_svPoolCreator.getSVCandPool(collisions, true);
    if (svPool.size() == 0) {
      m_qaRegistry.fill(HIST("hEmptyPool"), 1);
      return;
    }
    m_qaRegistry.fill(HIST("hEmptyPool"), 0);

    for (auto& svCand : svPool) {
      auto heTrack = tracks.rawIteratorAt(svCand.tr0Idx);
      auto prTrack = tracks.rawIteratorAt(svCand.tr1Idx);
      auto heTrackLabel = mcTrackLabels.rawIteratorAt(svCand.tr0Idx);
      auto prTrackLabel = mcTrackLabels.rawIteratorAt(svCand.tr1Idx);
      auto collBracket = svCand.collBracket;

      if (!heTrackLabel.has_mcParticle() || !prTrackLabel.has_mcParticle()) {
        continue;
      }

      auto mctrackHe3 = heTrackLabel.mcParticle_as<aod::McParticles>();
      auto mctrackHad = prTrackLabel.mcParticle_as<aod::McParticles>();

      if (std::abs(mctrackHe3.pdgCode()) != hePDG || std::abs(mctrackHad.pdgCode()) != prPDG || !mctrackHe3.has_mothers() || !mctrackHad.has_mothers()) {
        continue;
      }

      for (auto& mothertrackHe : mctrackHe3.mothers_as<aod::McParticles>()) {
        for (auto& mothertrackHad : mctrackHad.mothers_as<aod::McParticles>()) {

          if (mothertrackHe.globalIndex() != mothertrackHad.globalIndex() || std::abs(mothertrackHe.pdgCode()) != li4PDG || std::abs(mothertrackHe.y()) > 1) {
            continue;
          }

          he3HadCandidate he3Hadcand;
          if (!fillCandidateInfo(heTrack, prTrack, collBracket, collisions, he3Hadcand, tracks, /*mix*/ false)) {
            continue;
          }
          fillCandidateInfoMC(mctrackHe3, mctrackHad, mothertrackHe, he3Hadcand);
          fillHistograms(he3Hadcand);
          auto collision = collisions.rawIteratorAt(he3Hadcand.collisionID);
          fillTable(he3Hadcand, collision, /*isMC*/ true);
          filledMothers.push_back(mothertrackHe.globalIndex());
        }
      }
    }

    fillMcParticles(collisions, mcParticles, filledMothers);
  }
  PROCESS_SWITCH(he3hadronfemto, processMcPools, "Process MC pools", false);
};

WorkflowSpec defineDataProcessing(const ConfigContext& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<he3hadronfemto>(cfgc, TaskName{"he3hadronfemto"})};
}
