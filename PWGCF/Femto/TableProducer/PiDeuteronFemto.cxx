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

/// \file PiDeuteronFemto.cxx
/// \brief Analysis task for Deuteron-Pion femto analysis
/// \author CMY
/// \date 2025-04-10

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>

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
#include "PWGCF/Femto/DataModel/PionDeuteronTables.h"
#include "PWGLF/Utils/svPoolCreator.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using CollBracket = o2::math_utils::Bracket<int>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullDe, aod::pidTOFFullDe, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::TOFSignal, aod::TOFEvTime>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

enum Selections {
  kNoCuts = 0,
  kTrackCuts,
  kPID,
  kAll
};

} // namespace

struct PiDecandidate {

  float recoPtDe() const { return signDe * std::hypot(momDe[0], momDe[1]); }
  float recoPhiDe() const { return std::atan2(momDe[1], momDe[0]); }
  float recoEtaDe() const { return std::asinh(momDe[2] / std::abs(recoPtDe())); }
  float recoPtPi() const { return signPi * std::hypot(momPi[0], momPi[1]); }
  float recoPhiPi() const { return std::atan2(momPi[1], momPi[0]); }
  float recoEtaPi() const { return std::asinh(momPi[2] / std::abs(recoPtPi())); }

  std::array<float, 3> momDe = {99.f, 99.f, 99.f};
  std::array<float, 3> momPi = {99.f, 99.f, 99.f};

  float signDe = 1.f;
  float signPi = 1.f;
  float invMass = -10.f;
  float dcaxyDe = -10.f;
  float dcazDe = -10.f;
  float dcaxyPi = -10.f;
  float dcazPi = -10.f;

  uint16_t tpcSignalDe = 0u;
  uint16_t tpcSignalPi = 0u;
  float momDeTPC = -99.f;
  float momPiTPC = -99.f;
  uint8_t nTPCClustersDe = 0u;
  uint8_t sharedClustersDe = 0u;
  uint8_t sharedClustersPi = 0u;
  float chi2TPCDe = -10.f;
  float chi2TPCPi = -10.f;
  float nSigmaDe = -10.f;
  float nSigmaPi = -10.f;
  uint32_t pidTrkDe = 0xFFFFF; // PID in tracking
  uint32_t pidTrkPi = 0xFFFFF;
  float massTOFDe = -10;
  float massTOFPi = -10;
  uint32_t itsClSizeDe = 0u;
  uint32_t itsClSizePi = 0u;

  uint8_t nClsItsDe = 0u;
  uint8_t nClsItsPi = 0u;

  bool isBkgUS = false; // unlike sign
  bool isBkgEM = false; // event mixing

  int trackIDDe = -1;
  int trackIDPi = -1;

  // collision information
  int32_t collisionID = 0;
};

struct PiDeuteronFemto {

  Produces<aod::PionDeuteronTable> mOutputDataTable;
  Produces<aod::PionDeuteronMult> mOutputMultiplicityTable;

  // Selections
  Configurable<float> settingCutVertex{"settingCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> settingCutPinMinDe{"settingCutPinMinDe", 0.0f, "Minimum Pin for De"};
  Configurable<float> settingCutEta{"settingCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> settingCutChi2tpcLow{"settingCutChi2tpcLow", 0.0f, "Low cut on TPC chi2"};
  Configurable<float> settingCutChi2tpcHigh{"settingCutChi2tpcHigh", 999.f, "High cut on TPC chi2"};
  Configurable<float> settingCutInvMass{"settingCutInvMass", 0.0f, "Invariant mass upper limit"};
  Configurable<float> settingCutPtMinDePi{"settingCutPtMinDePi", 0.0f, "Minimum PT cut on DePi4"};
  Configurable<float> settingCutClSizeItsDe{"settingCutClSizeItsDe", 4.0f, "Minimum ITS cluster size for De"};
  Configurable<float> settingCutNCls{"settingCutNCls", 5.0f, "Minimum ITS Ncluster for tracks"};
  Configurable<float> settingCutChi2NClITS{"settingCutChi2NClITS", 999.f, "Maximum ITS Chi2 for tracks"};
  Configurable<float> settingCutNsigmaTPCPi{"settingCutNsigmaTPCPi", 3.0f, "Value of the TPC Nsigma cut on Pi"};
  Configurable<float> settingCutNsigmaTPCDe{"settingCutNsigmaTPCDe", 2.5f, "Value of the TPC Nsigma cut on De"};
  Configurable<float> settingCutNsigmaITSDe{"settingCutNsigmaITSDe", 2.5f, "Value of the ITD Nsigma cut on De"};
  Configurable<float> settingCutPinMinTOFPi{"settingCutPinMinTOFPi", 0.5f, "Minimum Pin to apply the TOF cut on Pions"};
  Configurable<float> settingCutPinMinTOFITSDe{"settingCutPinMinTOFITSDe", 1.2f, "Minimum pT to apply the TOF ITS cut on De"};
  Configurable<float> settingCutNsigmaTOFTPCDe{"settingCutNsigmaTOFTPCDe", 2.5f, "Value of the De TOF TPC Nsigma cut"};
  Configurable<float> settingCutNsigmaTOFTPCPi{"settingCutNsigmaTOFTPCPi", 3.0f, "Value of the Pion TOF TPC Nsigma cut"};
  Configurable<int> settingNoMixedEvents{"settingNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> settingEnableBkgUS{"settingEnableBkgUS", false, "Enable US background"};

  Configurable<bool> settingFillTable{"settingFillTable", false, "Enable table filling"};
  Configurable<float> settingCutPiptMin{"settingCutPiptMin", 0.14f, "Minimum PT cut on Pi"};
  Configurable<float> settingCutPiptMax{"settingCutPiptMax", 4.0f, "Maximum PT cut on Pi"};
  Configurable<float> settingCutDeptMin{"settingCutDeptMin", 0.6f, "Minimum PT cut on De"};
  Configurable<float> settingCutDeptMax{"settingCutDeptMax", 1.6f, "Maximum PT cut on De"};
  Configurable<float> settingCutPiDCAxyMin{"settingCutPiDCAxyMin", 0.3f, "DCAxy Min for Pi"};
  Configurable<float> settingCutPiDCAzMin{"settingCutPiDCAzMin", 0.3f, "DCAz Min for Pi"};
  Configurable<float> settingCutDeDCAzMin{"settingCutDeDCAzMin", 0.2f, "DCAxy Min for De"};
  Configurable<float> settingCutNsigTPCPrMin{"settingCutNsigTPCPrMin", 3.0f, "Minimum TPC Pr Nsigma cut on Pi"};
  Configurable<float> settingCutNsigTOFPrMin{"settingCutNsigTOFPrMin", 3.0f, "Minimum TOF Pr Nsigma cut on Pi"};

  Configurable<bool> settingSaveUSandLS{"settingSaveUSandLS", true, "Save All Pairs"};
  Configurable<bool> settingFillMultiplicity{"settingFillMultiplicity", false, "Fill multiplicity table"};

  // Zorro
  Configurable<bool> settingSkimmedProcessing{"settingSkimmedProcessing", false, "Skimmed dataset processing"};

  // CCDB options
  Configurable<double> settingDbz{"settingDbz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> settingCcdburl{"settingCcdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> settingGrpPath{"settingGrpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> settingGrpmagPath{"settingGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> settingLutPath{"settingLutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> settingGeoPath{"settingGeoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> settingPidPath{"settingPidPath", "", "Path to the PID response object"};

  Configurable<LabeledArray<double>> settingBetheBlochParams{"settingBetheBlochParams", {betheBlochDefault[0], 1, 6, {"De"}, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for De"};
  Configurable<bool> settingCompensatePIDinTracking{"settingCompensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};
  Configurable<int> settingMaterialCorrection{"settingMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Material correction type"};

  Preslice<TrackCandidates> mPerCol = aod::track::collisionId;

  // binning for EM background
  ConfigurableAxis axisVertex{"axisVertex", {30, -10, 10}, "Binning for multiplicity"};
  ConfigurableAxis axisCentrality{"axisCentrality", {40, 0, 100}, "Binning for centrality"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType binningPolicy{{axisVertex, axisCentrality}, true};
  SliceCache cache;
  SameKindPair<CollisionsFull, TrackCandidates, BinningType> mPair{binningPolicy, settingNoMixedEvents, -1, &cache};

  std::array<float, 6> mBBparamsDe;

  std::vector<int> mRecoCollisionIDs;
  std::vector<bool> mGoodCollisions;
  std::vector<SVCand> mTrackPairs;
  o2::vertexing::DCAFitterN<2> mFitter;

  int mRunNumber;
  float mDbz;
  Service<o2::ccdb::BasicCCDBManager> mCcdb;
  Zorro mZorro;
  OutputObj<ZorroSummary> mZorroSummary{"zorroSummary"};

  HistogramRegistry mQaRegistry{
    "QA",
    {{"hVtxZ", "Vertex distribution in Z;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
     {"hNcontributor", "Number of primary vertex contributor", {HistType::kTH1F, {{2000, 0.0f, 2000.0f}}}},
     {"hTrackSel", "Accepted tracks", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
     {"hEvents", "; Events;", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
     {"hEmptyPool", "svPoolCreator did not find track pairs false/true", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
     {"hdcaxyDe", ";DCA_{xy} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hdcazDe", ";DCA_{z} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hNClsDeITS", ";N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
     {"hNClsPiITS", ";N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
     {"hDePitInvMass", "; M(De + p) (GeV/#it{c}^{2})", {HistType::kTH1F, {{300, 3.74f, 4.34f}}}},
     {"hDePt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
     {"hPiPt", "Pt distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},
     {"h2dEdxDecandidates", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"h2NsigmaDeTPC", "NsigmaDe TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(De)", {HistType::kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaDeTPC_preselection", "NsigmaDe TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(De)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaDeTPC_preselecComp", "NsigmaDe TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(De)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NSigmaDeITS_preselection", "NsigmaDe ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} De", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {120, -3.0f, 3.0f}}}},
     {"h2NSigmaDeITS", "NsigmaDe ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} De", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {120, -3.0f, 3.0f}}}},
     {"h2NsigmaPiTPC", "NsigmaPi TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaPiTPC_preselection", "NsigmaDe TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(De)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaPiTOF", "NsigmaPi TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaPiTOF_preselection", "NsigmaPi TOF distribution; #iit{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
     {"hkStar_LS_M", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_LS_A", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_US_M", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_US_A", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hisBkgEM", "; isBkgEM;", {HistType::kTH1F, {{3, -1, 2}}}}},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  void init(o2::framework::InitContext&)
  {
    mZorroSummary.setObject(mZorro.getZorroSummary());
    mRunNumber = 0;

    mCcdb->setURL(settingCcdburl);
    mCcdb->setCaching(true);
    mCcdb->setLocalObjectValidityChecking();
    mCcdb->setFatalWhenNull(false);

    mFitter.setPropagateToPCA(true);
    mFitter.setMaxR(200.);
    mFitter.setMinParamChange(1e-3);
    mFitter.setMinRelChi2Change(0.9);
    mFitter.setMaxDZIni(1e9);
    mFitter.setMaxChi2(1e9);
    mFitter.setUseAbsDCA(true);
    int mat{static_cast<int>(settingMaterialCorrection)};
    mFitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

    const int numParticles = 5;
    for (int i = 0; i < numParticles; i++) {
      mBBparamsDe[i] = settingBetheBlochParams->get("De", Form("p%i", i));
    }
    mBBparamsDe[5] = settingBetheBlochParams->get("De", "resolution"); // wdf?

    std::vector<std::string> selectionLabels = {"All", "Track selection", "PID"};
    for (int i = 0; i < Selections::kAll; i++) {
      mQaRegistry.get<TH1>(HIST("hTrackSel"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i].c_str());
    }

    std::vector<std::string> eventsLabels = {"All", "Selected", "Zorro De events"};
    for (int i = 0; i < Selections::kAll; i++) {
      mQaRegistry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(i + 1, eventsLabels[i].c_str());
    }

    mQaRegistry.get<TH1>(HIST("hEmptyPool"))->GetXaxis()->SetBinLabel(1, "False");
    mQaRegistry.get<TH1>(HIST("hEmptyPool"))->GetXaxis()->SetBinLabel(2, "True");
  }

  void initCCDB(const aod::BCsWithTimestamps::iterator& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    if (settingSkimmedProcessing) {
      mZorro.initCCDB(mCcdb.service, bc.runNumber(), bc.timestamp(), "fDe");
      mZorro.populateHistRegistry(mQaRegistry, bc.runNumber());
    }
    mRunNumber = bc.runNumber();
    const float defaultBzValue = -999.0f;
    auto run3GrpTimestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = mCcdb->getForTimeStamp<o2::parameters::GRPObject>(settingGrpPath, run3GrpTimestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (settingDbz < defaultBzValue) {
        // Fetch magnetic field from ccdb for current collision
        mDbz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3GrpTimestamp << " with magnetic field of " << mDbz << " kZG";
      } else {
        mDbz = settingDbz;
      }
    } else {
      grpmag = mCcdb->getForTimeStamp<o2::parameters::GRPMagField>(settingGrpmagPath, run3GrpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << settingGrpmagPath << " of object GRPMagField and " << settingGrpPath << " of object GRPObject for timestamp " << run3GrpTimestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (settingDbz < defaultBzValue) {
        // Fetch magnetic field from ccdb for current collision
        mDbz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3GrpTimestamp << " with magnetic field of " << mDbz << " kZG";
      } else {
        mDbz = settingDbz;
      }
    }
  }

  // ==================================================================================================================

  template <bool isMC, typename Tcollision>
  bool selectCollision(const Tcollision& collision, const aod::BCsWithTimestamps&)
  {
    mQaRegistry.fill(HIST("hEvents"), 0);

    if constexpr (isMC) {
      if (/*!collision.sel8() ||*/ std::abs(collision.posZ()) > settingCutVertex) {
        return false;
      }
    } else {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8() || std::abs(collision.posZ()) > settingCutVertex) {
        return false;
      }
      if (settingSkimmedProcessing) {
        bool zorroSelected = mZorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC());
        if (zorroSelected) {
          mQaRegistry.fill(HIST("hEvents"), 2);
        }
      }
    }

    mQaRegistry.fill(HIST("hEvents"), 1);
    mQaRegistry.fill(HIST("hNcontributor"), collision.numContrib());
    mQaRegistry.fill(HIST("hVtxZ"), collision.posZ());
    return true;
  }

  template <typename Ttrack>
  bool selectTrack(const Ttrack& candidate)
  {
    if (std::abs(candidate.eta()) > settingCutEta) {
      return false;
    }
    const int minTPCNClsFound = 90;
    const int minTPCNClsCrossedRows = 100;
    const float crossedRowsToFindableRatio = 0.83f;
    if (candidate.itsNCls() < settingCutNCls ||
        candidate.tpcNClsFound() < minTPCNClsFound ||
        candidate.tpcNClsCrossedRows() < minTPCNClsCrossedRows ||
        candidate.tpcNClsCrossedRows() < crossedRowsToFindableRatio * candidate.tpcNClsFindable() ||
        candidate.tpcChi2NCl() > settingCutChi2tpcHigh ||
        candidate.tpcChi2NCl() < settingCutChi2tpcLow ||
        candidate.itsChi2NCl() > settingCutChi2NClITS) {
      return false;
    }

    return true;
  }

  template <typename Ttrack>
  bool selectionPIDPion(const Ttrack& candidate)
  {
    auto tpcNSigmaPi = candidate.tpcNSigmaPi();
    mQaRegistry.fill(HIST("h2NsigmaPiTPC_preselection"), candidate.tpcInnerParam(), tpcNSigmaPi);
    if (candidate.hasTOF() && candidate.tpcInnerParam() > settingCutPinMinTOFPi) {
      auto tofNSigmaPi = candidate.tofNSigmaPi();

      if (std::abs(tpcNSigmaPi) > settingCutNsigmaTPCPi) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaPiTOF_preselection"), candidate.pt(), tofNSigmaPi);
      if (std::sqrt(tofNSigmaPi * tofNSigmaPi + tpcNSigmaPi * tpcNSigmaPi) > settingCutNsigmaTOFTPCPi) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaPiTPC"), candidate.pt(), tpcNSigmaPi);
      mQaRegistry.fill(HIST("h2NsigmaPiTOF"), candidate.pt(), tofNSigmaPi);
      return true;
    } else if (std::abs(tpcNSigmaPi) < settingCutNsigmaTPCPi) {
      mQaRegistry.fill(HIST("h2NsigmaPiTPC"), candidate.pt(), tpcNSigmaPi);
      return true;
    }
    return false;
  }

  template <typename Ttrack>
  float computeNSigmaDe(const Ttrack& candidate)
  {
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(candidate.tpcInnerParam() / constants::physics::MassDeuteron), mBBparamsDe[0], mBBparamsDe[1], mBBparamsDe[2], mBBparamsDe[3], mBBparamsDe[4]);
    double resoTPC{expTPCSignal * mBBparamsDe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename Ttrack>
  bool selectionPIDDe(const Ttrack& candidate)
  {
    float tpcInnerParam = candidate.tpcInnerParam();

    if (tpcInnerParam < settingCutPinMinDe) {
      return false;
    }

    auto tpcNSigmaDe = computeNSigmaDe(candidate);
    mQaRegistry.fill(HIST("h2NsigmaDeTPC_preselection"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
    mQaRegistry.fill(HIST("h2NsigmaDeTPC_preselecComp"), candidate.sign() * candidate.pt(), candidate.tpcNSigmaDe());
    if (candidate.hasTOF() && candidate.tpcInnerParam() > settingCutPinMinTOFITSDe) {
      auto tofNSigmaDe = candidate.tofNSigmaDe();
      if (std::abs(tpcNSigmaDe) > settingCutNsigmaTPCDe) {
        return false;
      }
      if (std::sqrt(tofNSigmaDe * tofNSigmaDe + tpcNSigmaDe * tpcNSigmaDe) > settingCutNsigmaTOFTPCDe) {
        return false;
      }
      o2::aod::ITSResponse mResponseITS;
      auto itsnSigmaDe = mResponseITS.nSigmaITS<o2::track::PID::Deuteron>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());

      mQaRegistry.fill(HIST("h2NSigmaDeITS_preselection"), candidate.sign() * candidate.pt(), itsnSigmaDe);
      if (itsnSigmaDe < settingCutNsigmaITSDe) {
        return false;
      }
      mQaRegistry.fill(HIST("h2dEdxDecandidates"), candidate.sign() * tpcInnerParam, candidate.tpcSignal());
      mQaRegistry.fill(HIST("h2NsigmaDeTPC"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
      mQaRegistry.fill(HIST("h2NSigmaDeITS"), candidate.sign() * candidate.pt(), itsnSigmaDe);
      return true;
    } else if (std::abs(tpcNSigmaDe) < settingCutNsigmaTPCDe) {
      mQaRegistry.fill(HIST("h2NsigmaDeTPC"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
      return true;
    }
    return false;
  }

  // ==================================================================================================================

  template <typename Ttrack, typename Tcollisions, typename Ttracks>
  bool fillCandidateInfo(const Ttrack& trackDe, const Ttrack& trackPi, const CollBracket& collBracket, const Tcollisions& collisions, PiDecandidate& piDecand, const Ttracks& /*trackTable*/, bool isMixedEvent)
  {
    const int numCoordinates = 3;
    if (!isMixedEvent) {
      auto trackCovDe = getTrackParCov(trackDe);
      auto trackCovPi = getTrackParCov(trackPi);
      int nCand = 0;
      try {
        nCand = mFitter.process(trackCovDe, trackCovPi);
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
        const auto& pca = mFitter.getPCACandidate();
        float distance = 0;
        for (int i = 0; i < numCoordinates; i++) {
          distance += (pca[i] - collVtx[i]) * (pca[i] - collVtx[i]);
        }
        if (distanceMin < 0 || distance < distanceMin) {
          distanceMin = distance;
          collIdxMin = collIdx;
        }
      }

      if (!mGoodCollisions[collIdxMin]) {
        return false;
      }
      piDecand.collisionID = collIdxMin;
    } else {
      piDecand.collisionID = collBracket.getMin();
    }

    piDecand.momDe = std::array{trackDe.px(), trackDe.py(), trackDe.pz()};
    for (int i = 0; i < numCoordinates; i++)
      piDecand.momPi = std::array{trackPi.px(), trackPi.py(), trackPi.pz()};
    float invMass = 0;
    invMass = RecoDecay::m(std::array{piDecand.momDe, piDecand.momPi}, std::array{o2::constants::physics::MassDeuteron, o2::constants::physics::MassPiPlus});
    if (settingCutInvMass > 0 && invMass > settingCutInvMass) {
      return false;
    }
    float ptDePi = std::hypot(piDecand.momDe[0] + piDecand.momPi[0], piDecand.momDe[1] + piDecand.momPi[1]);
    if (ptDePi < settingCutPtMinDePi) {
      return false;
    }

    piDecand.signDe = trackDe.sign();
    piDecand.signPi = trackPi.sign();

    piDecand.dcaxyDe = trackDe.dcaXY();
    piDecand.dcaxyPi = trackPi.dcaXY();

    piDecand.dcazDe = trackDe.dcaZ();
    piDecand.dcazPi = trackPi.dcaZ();

    piDecand.tpcSignalDe = trackDe.tpcSignal();
    piDecand.momDeTPC = trackDe.tpcInnerParam();
    piDecand.tpcSignalPi = trackPi.tpcSignal();
    piDecand.momPiTPC = trackPi.tpcInnerParam();

    piDecand.nTPCClustersDe = trackDe.tpcNClsFound();
    piDecand.nSigmaDe = computeNSigmaDe(trackDe);
    piDecand.nSigmaPi = trackPi.tpcNSigmaPi();

    piDecand.chi2TPCDe = trackDe.tpcChi2NCl();
    piDecand.chi2TPCPi = trackPi.tpcChi2NCl();

    piDecand.pidTrkDe = trackDe.pidForTracking();
    piDecand.pidTrkPi = trackPi.pidForTracking();

    piDecand.itsClSizeDe = trackDe.itsClusterSizes();
    piDecand.itsClSizePi = trackPi.itsClusterSizes();

    piDecand.nClsItsDe = trackDe.itsNCls();
    piDecand.nClsItsPi = trackPi.itsNCls();

    piDecand.sharedClustersDe = trackDe.tpcNClsShared();
    piDecand.sharedClustersPi = trackPi.tpcNClsShared();

    piDecand.isBkgUS = trackDe.sign() * trackPi.sign() < 0;
    piDecand.isBkgEM = isMixedEvent;

    piDecand.invMass = invMass;

    piDecand.trackIDDe = trackDe.globalIndex();
    piDecand.trackIDPi = trackPi.globalIndex();

    if (trackDe.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackDe);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      float tpcInnerParamDe = trackDe.tpcInnerParam();
      piDecand.massTOFDe = tpcInnerParamDe * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    if (trackPi.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackPi);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      piDecand.massTOFPi = trackPi.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
    }

    return true;
  }

  template <typename Ttrack>
  void pairTracksSameEvent(const Ttrack& tracks)
  {
    // LOG(info) << "Number of tracks: " << tracks.size();
    for (const auto& track0 : tracks) {

      mQaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);

      if (!selectTrack(track0)) {
        continue;
      }
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);

      if (!selectionPIDDe(track0)) {
        continue;
      }
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kPID);

      for (const auto& track1 : tracks) {
        if (track0 == track1) {
          continue;
        }

        if (!settingSaveUSandLS) {
          if (!settingEnableBkgUS && (track0.sign() * track1.sign() < 0)) {
            continue;
          }
          if (settingEnableBkgUS && (track0.sign() * track1.sign() > 0)) {
            continue;
          }
        }

        if (!selectTrack(track1) || !selectionPIDPion(track1)) {
          continue;
        }

        SVCand trackPair;
        trackPair.tr0Idx = track0.globalIndex();
        trackPair.tr1Idx = track1.globalIndex();
        const int collIdx = track0.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        trackPair.collBracket = collBracket;
        mTrackPairs.push_back(trackPair);
      }
    }
  }

  template <typename T>
  void pairTracksEventMixing(T& DeCands, T& pionCands)
  {
    for (const auto& DeCand : DeCands) {
      if (!selectTrack(DeCand) || !selectionPIDDe(DeCand)) {
        continue;
      }
      for (const auto& pionCand : pionCands) {
        if (!selectTrack(pionCand) || !selectionPIDPion(pionCand)) {
          continue;
        }

        SVCand trackPair;
        trackPair.tr0Idx = DeCand.globalIndex();
        trackPair.tr1Idx = pionCand.globalIndex();
        const int collIdx = DeCand.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        trackPair.collBracket = collBracket;
        mTrackPairs.push_back(trackPair);
      }
    }
  }

  template <typename Tcoll>
  void fillTable(const PiDecandidate& piDecand, const Tcoll& collision)
  {
    mOutputDataTable(
      piDecand.recoPtDe(),
      piDecand.recoEtaDe(),
      piDecand.recoPhiDe(),
      piDecand.recoPtPi(),
      piDecand.recoEtaPi(),
      piDecand.recoPhiPi(),
      piDecand.dcaxyDe,
      piDecand.dcazDe,
      piDecand.dcaxyPi,
      piDecand.dcazPi,
      piDecand.tpcSignalDe,
      piDecand.momDeTPC,
      piDecand.tpcSignalPi,
      piDecand.momPiTPC,
      piDecand.nTPCClustersDe,
      piDecand.nSigmaDe,
      piDecand.nSigmaPi,
      piDecand.chi2TPCDe,
      piDecand.chi2TPCPi,
      piDecand.massTOFDe,
      piDecand.massTOFPi,
      piDecand.pidTrkDe,
      piDecand.pidTrkPi,
      piDecand.itsClSizeDe,
      piDecand.itsClSizePi,
      piDecand.sharedClustersDe,
      piDecand.sharedClustersPi,
      piDecand.isBkgUS,
      piDecand.isBkgEM);
    if (settingFillMultiplicity) {
      mOutputMultiplicityTable(
        collision.globalIndex(),
        collision.posZ(),
        collision.numContrib(),
        collision.centFT0C(),
        collision.multFT0C());
    }
  }

  void fillHistograms(const PiDecandidate& piDecand)
  {
    mQaRegistry.fill(HIST("hDePt"), piDecand.recoPtDe());
    mQaRegistry.fill(HIST("hPiPt"), piDecand.recoPtPi());
    mQaRegistry.fill(HIST("hDePitInvMass"), piDecand.invMass);
    mQaRegistry.fill(HIST("hdcaxyDe"), piDecand.dcaxyDe);
    mQaRegistry.fill(HIST("hdcazDe"), piDecand.dcazDe);
    mQaRegistry.fill(HIST("hNClsDeITS"), piDecand.nClsItsDe);
    mQaRegistry.fill(HIST("hNClsPiITS"), piDecand.nClsItsPi);
    mQaRegistry.fill(HIST("hisBkgEM"), piDecand.isBkgEM);
  }

  double computePrTPCnsig(double InnerParamTPCHad, double SignalTPCHad)
  {
    double m_BBparamsProton[6] = {-54.42066571222577, 0.2857381250239097, 1.247140602468868, 0.6297483918147729, 2.985438833884555, 0.09};

    float TPCinnerParam = InnerParamTPCHad;
    float expTPCSignal = o2::tpc::BetheBlochAleph((TPCinnerParam / 0.9382721), m_BBparamsProton[0], m_BBparamsProton[1], m_BBparamsProton[2], m_BBparamsProton[3], m_BBparamsProton[4]);
    double resoTPC{expTPCSignal * m_BBparamsProton[5]};
    return ((SignalTPCHad - expTPCSignal) / resoTPC);
  }

  double tofNSigmaCalculation(double MassTOFHad, double ptHad)
  {
    double fExpTOFMassHad = 0.9487; // Proton mass in TOF
    const float kp0 = 1.22204e-02;
    const float kp1 = 7.48467e-01;

    double fSigmaTOFMassHad = (kp0 * TMath::Exp(kp1 * TMath::Abs(ptHad))) * fExpTOFMassHad;
    double fNSigmaTOFHad = (MassTOFHad - fExpTOFMassHad) / fSigmaTOFMassHad;
    return fNSigmaTOFHad;
  }

  double computeKstar(const PiDecandidate& piDecand)
  {
    TLorentzVector he3, hadron;
    float massHe3 = 2.80839;
    float massHad = 0.1395704;
    he3.SetPtEtaPhiM(abs(piDecand.recoPtDe()), piDecand.recoEtaDe(), piDecand.recoPhiDe(), massHe3);
    hadron.SetPtEtaPhiM(abs(piDecand.recoPtPi()), piDecand.recoEtaPi(), piDecand.recoPhiPi(), massHad);

    TLorentzVector p_total_lab = he3 + hadron;
    TVector3 v_cm = p_total_lab.BoostVector();
    TLorentzVector p1_cm = he3;
    TLorentzVector p2_cm = hadron;
    p1_cm.Boost(-v_cm);
    p2_cm.Boost(-v_cm);
    TLorentzVector p_diff_cm = p1_cm - p2_cm;
    double kStar = sqrt(p_diff_cm.X() * p_diff_cm.X() + p_diff_cm.Y() * p_diff_cm.Y() + p_diff_cm.Z() * p_diff_cm.Z());
    return kStar / 2.0;
  }

  void fillKstar(const PiDecandidate& piDecand)
  {
    double PrTPCnsigma = computePrTPCnsig(piDecand.momPiTPC, piDecand.tpcSignalPi);
    double PrTOFnsigma = tofNSigmaCalculation(piDecand.massTOFPi, piDecand.recoPtPi());
    if (abs(PrTPCnsigma) < settingCutNsigTPCPrMin)
      return;
    if (abs(PrTOFnsigma) < settingCutNsigTOFPrMin)
      return;
    float DeDCAxyMin = 0.015 + 0.0305 / TMath::Power(piDecand.recoPtDe(), 1.1);
    if (abs(piDecand.dcaxyDe) > DeDCAxyMin || abs(piDecand.dcazDe) > settingCutDeDCAzMin || abs(piDecand.dcaxyPi) > settingCutPiDCAxyMin || abs(piDecand.dcazPi) > settingCutPiDCAzMin)
      return;
    if (std::abs(piDecand.recoPtPi()) < settingCutPiptMin || std::abs(piDecand.recoPtPi()) > settingCutPiptMax)
      return;
    if (std::abs(piDecand.recoPtDe()) < settingCutDeptMin || std::abs(piDecand.recoPtDe()) > settingCutDeptMax)
      return;

    fillHistograms(piDecand);

    double kstar = computeKstar(piDecand);
    if (piDecand.isBkgUS == 0) {
      if (piDecand.recoPtDe() > 0) {
        mQaRegistry.fill(HIST("hkStar_LS_M"), kstar);
      } else {
        mQaRegistry.fill(HIST("hkStar_LS_A"), kstar);
      }
    } else {
      if (piDecand.recoPtDe() > 0) {
        mQaRegistry.fill(HIST("hkStar_US_M"), kstar);
      } else {
        mQaRegistry.fill(HIST("hkStar_US_A"), kstar);
      }
    }
  }

  // ==================================================================================================================

  template <typename Tcollisions, typename Ttracks>
  void fillPairs(const Tcollisions& collisions, const Ttracks& tracks, const bool isMixedEvent)
  {
    for (const auto& trackPair : mTrackPairs) {

      auto deTrack = tracks.rawIteratorAt(trackPair.tr0Idx);
      auto piTrack = tracks.rawIteratorAt(trackPair.tr1Idx);
      auto collBracket = trackPair.collBracket;

      PiDecandidate piDecand;
      if (!fillCandidateInfo(deTrack, piTrack, collBracket, collisions, piDecand, tracks, isMixedEvent)) {
        continue;
      }

      fillKstar(piDecand);

      auto collision = collisions.rawIteratorAt(piDecand.collisionID);

      if (settingFillTable) {
        fillTable(piDecand, collision);
      }
    }
  }

  // ==================================================================================================================

  void processSameEvent(const CollisionsFull& collisions, const TrackCandidates& tracks, const aod::BCsWithTimestamps& bcs)
  {
    mGoodCollisions.clear();
    mGoodCollisions.resize(collisions.size(), false);

    for (const auto& collision : collisions) {

      mTrackPairs.clear();

      if (!selectCollision</*isMC*/ false>(collision, bcs)) {
        continue;
      }

      mGoodCollisions[collision.globalIndex()] = true;
      const uint64_t collIdx = collision.globalIndex();
      auto trackTableThisCollision = tracks.sliceBy(mPerCol, collIdx);
      trackTableThisCollision.bindExternalIndices(&tracks);

      pairTracksSameEvent(trackTableThisCollision);

      if (mTrackPairs.size() == 0) {
        continue;
      }

      fillPairs(collisions, tracks, /*isMixedEvent*/ false);
    }
  }
  PROCESS_SWITCH(PiDeuteronFemto, processSameEvent, "Process Same event", false);

  void processMixedEvent(const CollisionsFull& collisions, const TrackCandidates& tracks)
  {
    LOG(debug) << "Processing mixed event";
    mTrackPairs.clear();

    for (const auto& [c1, tracks1, c2, tracks2] : mPair) {
      if (!c1.sel8() || !c2.sel8()) {
        continue;
      }

      mQaRegistry.fill(HIST("hNcontributor"), c1.numContrib());
      mQaRegistry.fill(HIST("hVtxZ"), c1.posZ());

      pairTracksEventMixing(tracks1, tracks2);
      pairTracksEventMixing(tracks2, tracks1);
    }

    fillPairs(collisions, tracks, /*isMixedEvent*/ true);
  }
  PROCESS_SWITCH(PiDeuteronFemto, processMixedEvent, "Process Mixed event", false);
};

WorkflowSpec defineDataProcessing(const ConfigContext& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PiDeuteronFemto>(cfgc)};
}
