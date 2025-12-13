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

/// \file HadNucleiFemto.cxx
/// \brief Analysis task for Nuclei-Hadron femto analysis
/// \author CMY
/// \date 2025-04-10

#include "PWGCF/Femto/FemtoNuclei/DataModel/HadronNucleiTables.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldMath.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFHypernucleiTables.h"
#include "PWGLF/Utils/svPoolCreator.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Boost.h"
#include "Math/Vector4D.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TMath.h>
#include <TObjArray.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <iterator> // std::prev
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using CollBracket = o2::math_utils::Bracket<int>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullDe, aod::pidTOFFullDe, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::TOFSignal, aod::TOFEvTime>;

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

float MassHad = 0;

} // namespace

struct HadNucandidate {

  float recoPtNu() const { return signNu * std::hypot(momNu[0], momNu[1]); }
  float recoPhiNu() const { return std::atan2(momNu[1], momNu[0]); }
  float recoEtaNu() const { return std::asinh(momNu[2] / std::abs(recoPtNu())); }
  float recoPtHad() const { return signHad * std::hypot(momHad[0], momHad[1]); }
  float recoPhiHad() const { return std::atan2(momHad[1], momHad[0]); }
  float recoEtaHad() const { return std::asinh(momHad[2] / std::abs(recoPtHad())); }

  std::array<float, 3> momNu = {99.f, 99.f, 99.f};
  std::array<float, 3> momHad = {99.f, 99.f, 99.f};

  float ptHe3 = 1.f;
  float etaHe3 = 1.f;
  float signNu = 1.f;
  float signHad = 1.f;
  float invMass = -10.f;
  float dcaxyNu = -10.f;
  float dcazNu = -10.f;
  float dcaxyHad = -10.f;
  float dcazHad = -10.f;

  uint16_t tpcSignalNu = 0u;
  uint16_t tpcSignalHad = 0u;
  float momNuTPC = -99.f;
  float momHadTPC = -99.f;
  uint8_t nTPCClustersNu = 0u;
  uint8_t sharedClustersNu = 0u;
  uint8_t sharedClustersHad = 0u;
  float chi2TPCNu = -10.f;
  float chi2TPCHad = -10.f;
  float nSigmaNu = -10.f;
  float nSigmaHad = -10.f;
  float tpcPrnsigma = -10.f;
  float tofPrnsigma = -10.f;
  uint32_t pidTrkNu = 0xFFFFF; // PID in tracking
  uint32_t pidTrkHad = 0xFFFFF;
  float massTOFNu = -10;
  float massTOFHad = -10;
  uint32_t itsClSizeNu = 0u;
  uint32_t itsClSizeHad = 0u;

  uint8_t nClsItsNu = 0u;
  uint8_t nClsItsHad = 0u;

  bool isBkgUS = false; // unlike sign
  bool isBkgEM = false; // event mixing

  int trackIDNu = -1;
  int trackIDHad = -1;

  float kstar = 1.f;
  float mT = 1.f;

  // collision information
  int32_t collisionID = 0;
  float cent = 1.f;
};

struct HadNucleiFemto {

  Produces<aod::HadronNucleiTable> mOutputDataTable;
  Produces<aod::HadronHyperTable> mOutputHyperDataTable;
  Produces<aod::HadronNucleiMult> mOutputMultiplicityTable;

  // Selections
  Configurable<int> settingHadPDGCode{"settingHadPDGCode", 211, "Hadron - PDG code"};

  Configurable<float> settingCutVertex{"settingCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> settingCutPinMinDe{"settingCutPinMinDe", 0.0f, "Minimum Pin for De"};
  Configurable<float> settingCutEta{"settingCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> settingCutChi2tpcLow{"settingCutChi2tpcLow", 0.0f, "Low cut on TPC chi2"};
  Configurable<float> settingCutChi2tpcHigh{"settingCutChi2tpcHigh", 999.f, "High cut on TPC chi2"};
  Configurable<float> settingCutChi2tpcLowPion{"settingCutChi2tpcLowPion", 0.5f, "Low cut on TPC chi2 only for pion"};
  Configurable<float> settingCutChi2tpcHighPion{"settingCutChi2tpcHighPion", 4.f, "High cut on TPC chi2 only for pion"};
  Configurable<float> settingCutInvMass{"settingCutInvMass", 0.0f, "Invariant mass upper limit"};
  Configurable<float> settingCutPtMinDePi{"settingCutPtMinDePi", 0.0f, "Minimum PT cut on DePi4"};
  Configurable<float> settingCutClSizeItsDe{"settingCutClSizeItsDe", 4.0f, "Minimum ITS cluster size for De"};
  Configurable<float> settingCutNCls{"settingCutNCls", 5.0f, "Minimum ITS Ncluster for tracks"};
  Configurable<float> settingCutTPCChi2He{"settingCutTPCChi2He", 0.0f, "Minimum tpcChi2He for Hyper He3"};
  Configurable<float> settingCutAverClsSizeHe{"settingCutAverClsSizeHe", 0.0f, "Minimum averClusSizeHe for Hyper He3"};
  Configurable<float> settingCutChi2NClITS{"settingCutChi2NClITS", 999.f, "Maximum ITS Chi2 for tracks"};
  Configurable<float> settingCutChi2NClITSPion{"settingCutChi2NClITSPion", 36.f, "Maximum ITS Chi2 for tracks only for pion"};
  Configurable<float> settingCutNsigmaTPCHad{"settingCutNsigmaTPCHad", 3.0f, "Value of the TPC Nsigma cut on Had"};
  Configurable<float> settingCutNsigmaTOFHad{"settingCutNsigmaTOFHad", 3.0f, "Value of the hsdron TOF Nsigma cut"};
  Configurable<float> settingCutNsigmaTPCDe{"settingCutNsigmaTPCDe", 2.5f, "Value of the TPC Nsigma cut on De"};
  Configurable<float> settingCutNsigmaITSDe{"settingCutNsigmaITSDe", 2.5f, "Value of the ITD Nsigma cut on De"};
  Configurable<float> settingCutPinMinTOFHad{"settingCutPinMinTOFHad", 0.5f, "Minimum Pin to apply the TOF cut on hadrons"};
  Configurable<float> settingCutPinMinTOFITSDe{"settingCutPinMinTOFITSDe", 1.2f, "Minimum p to apply the TOF ITS cut on De"};
  Configurable<float> settingCutNsigmaTOFTPCDe{"settingCutNsigmaTOFTPCDe", 2.5f, "Value of the De TOF TPC combNsigma cut"};
  Configurable<float> settingCutNsigmaTOFTPCHad{"settingCutNsigmaTOFTPCHad", 3.0f, "Value of the hsdron TOF TPC combNsigma cut"};
  Configurable<int> settingNoMixedEvents{"settingNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> settingEnableBkgUS{"settingEnableBkgUS", false, "Enable US background"};
  Configurable<bool> settingSaferME{"settingSaferME", false, "For Safer ME"};

  Configurable<bool> settingFillTable{"settingFillTable", false, "Enable table filling"};
  Configurable<float> settingCutHadptMin{"settingCutHadptMin", 0.14f, "Minimum PT cut on Had"};
  Configurable<float> settingCutHadptMax{"settingCutHadptMax", 4.0f, "Maximum PT cut on Had"};
  Configurable<float> settingCutDeptMin{"settingCutDeptMin", 0.6f, "Minimum PT cut on De"};
  Configurable<float> settingCutDeptMax{"settingCutDeptMax", 1.6f, "Maximum PT cut on De"};
  Configurable<float> settingCutHadDCAxyMin{"settingCutHadDCAxyMin", 0.3f, "DCAxy Min for Had"};
  Configurable<float> settingCutHadDCAzMin{"settingCutHadDCAzMin", 0.3f, "DCAz Min for Had"};
  Configurable<float> settingCutDeDCAzMin{"settingCutDeDCAzMin", 0.2f, "DCAxy Min for De"};
  Configurable<float> settingCutNsigTPCPrMin{"settingCutNsigTPCPrMin", 3.0f, "Minimum TPC Pr Nsigma cut on Pi"};
  Configurable<float> settingCutNsigTOFPrMin{"settingCutNsigTOFPrMin", 3.0f, "Minimum TOF Pr Nsigma cut on Pi"};

  Configurable<bool> settingSaveUSandLS{"settingSaveUSandLS", true, "Save All Pairs"};
  Configurable<bool> settingFillMultiplicity{"settingFillMultiplicity", false, "Fill multiplicity table"};
  Configurable<bool> settingUseBBcomputeDeNsigma{"settingUseBBcomputeDeNsigma", false, "Use BB params to compute De TPC Nsigma"};

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
  PresliceUnsorted<o2::aod::DataHypCandsWColl> hypPerCol = o2::aod::hyperrec::collisionId;

  // binning for EM background
  ConfigurableAxis axisVertex{"axisVertex", {30, -10, 10}, "Binning for vtxz"};
  ConfigurableAxis axisCentrality{"axisCentrality", {40, 0, 100}, "Binning for centrality"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType binningPolicy{{axisVertex, axisCentrality}, true};
  SliceCache cache;
  SameKindPair<CollisionsFull, TrackCandidates, BinningType> mPair{binningPolicy, settingNoMixedEvents, -1, &cache};
  // Pair<CollisionsFull, TrackCandidates, o2::aod::DataHypCandsWColl, BinningType> hyperPair{binningPolicy, settingNoMixedEvents, -1, &cache};

  std::array<float, 6> mBBparamsDe;

  std::vector<int> mRecoCollisionIDs;
  std::vector<bool> mGoodCollisions;
  std::vector<SVCand> mTrackPairs;
  std::vector<SVCand> mTrackHypPairs;
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
     {"hCentrality", "Centrality", {HistType::kTH1F, {{100, 0.0f, 100.0f}}}},
     {"hTrackSel", "Accepted tracks", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
     {"hSkipReasons", "Why storedEvent skipped;Reason;Counts", {HistType::kTH1F, {{5, -0.5, 4.5}}}},
     {"hEvents", "; Events;", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
     {"hEmptyPool", "svPoolCreator did not find track pairs false/true", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
     {"hdcaxyNu", ";DCA_{xy} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hdcazNu", ";DCA_{z} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hdcazNu_min", ";DCA_{z}-min (cm)", {HistType::kTH1F, {{20, -1.0f, 1.0f}}}},
     {"hNClsNuITS", ";N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
     {"hNClsHadITS", ";N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
     {"hNuHadtInvMass", "; M(Nu + p) (GeV/#it{c}^{2})", {HistType::kTH1F, {{300, 3.74f, 4.34f}}}},
     {"hNuPt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
     {"hHadPt", "Pt distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},
     {"hSingleNuPt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
     {"hNuPin", "#it{p} distribution; #it{p} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
     {"hHadPin", "P distribution; #it{p} (GeV/#it{c})", {HistType::kTH1F, {{120, -4.0f, 4.0f}}}},
     {"hSingleNuPin", "#it{p} distribution; #it{p} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},

     {"hHe3TPCnsigma", "NsigmaHe3 TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(He3)", {HistType::kTH2F, {{100, -2.0f, 2.0f}, {200, -5.0f, 5.0f}}}},
     {"hHe3P", "Pin distribution; p (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},
     {"hHe3P_preselected", "Pin distribution_preselected; p (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},
     {"hNuEta", "eta distribution; #eta(Nu)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hHadEta", "eta distribution; #eta(had)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hNuPhi", "phi distribution; phi(Nu)", {HistType::kTH1F, {{600, -4.0f, 4.0f}}}},
     {"hHadPhi", "phi distribution; phi(had)", {HistType::kTH1F, {{600, -4.0f, 4.0f}}}},
     {"h2dEdxNucandidates", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"h2dEdx", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"h2NsigmaNuTPC", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{100, -2.0f, 2.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaNuComb", "NsigmaNu TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(Nu)", {HistType::kTH2F, {{100, -2.0f, 2.0f}, {100, 0.0f, 5.0f}}}},
     {"h2NsigmaHadComb", "NsigmaHad TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(had)", {HistType::kTH2F, {{100, -2.0f, 2.0f}, {100, 0.0f, 5.0f}}}},
     {"h2NsigmaNuTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaNuTPC_preselecComp", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NSigmaNuITS_preselection", "NsigmaNu ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} Nu", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {120, -3.0f, 3.0f}}}},
     {"h2NSigmaNuITS", "NsigmaNu ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} Nu", {HistType::kTH2F, {{100, -2.0f, 2.0f}, {120, -3.0f, 3.0f}}}},
     {"h2NsigmaHadTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaHadTOF", "NsigmaHad TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaNuTOF", "NsigmaNu TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(Nu)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadTOF_preselection", "NsigmaHad TOF distribution; #iit{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
     {"hkStaVsmTVsCent_LS_M", ";kStar (GeV/c);mT (GeV/#it{c}^{2});Centrality", {HistType::kTH3F, {{300, 0.0f, 3.0f}, {100, 0.2, 3.2}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsmTVsCent_LS_A", ";kStar (GeV/c);mT (GeV/#it{c}^{2});Centrality", {HistType::kTH3F, {{300, 0.0f, 3.0f}, {100, 0.2, 3.2}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsmTVsCent_US_M", ";kStar (GeV/c);mT (GeV/#it{c}^{2});Centrality", {HistType::kTH3F, {{300, 0.0f, 3.0f}, {100, 0.2, 3.2}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsmTVsCent_US_A", ";kStar (GeV/c);mT (GeV/#it{c}^{2});Centrality", {HistType::kTH3F, {{300, 0.0f, 3.0f}, {100, 0.2, 3.2}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsmT_LS_M", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsmT_LS_A", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsmT_US_M", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsmT_US_A", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},

     {"hNHypsPerPrevColl", "Number of V0Hypers in previous collision used for mixing;N_{V0Hypers};Entries", {HistType::kTH2F, {{4000, 0.0f, 4000.0f}, {50, -0.5, 49.5}}}},
     {"hkStar_LS_M", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_LS_A", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_US_M", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_US_A", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"h2NsigmaHadPrTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadPrTOF", "NsigmaHad TOF distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"hisBkgEM", "; isBkgEM;", {HistType::kTH1F, {{3, -1, 2}}}}},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  int numOfCentBins = 40;
  int numOfVertexZBins = 30;
  float Vz_low = -10.0f;
  float Vz_high = 10.0f;
  float Vz_step = (Vz_high - Vz_low) / numOfVertexZBins;

  struct EventRef {
    uint64_t collisionId;
  };

  struct PoolBin {
    std::deque<EventRef> events;
  };

  std::vector<PoolBin> All_Event_pool;
  bool isInitialized = false;

  int nPoolBins() const { return numOfVertexZBins * numOfCentBins; }

  void initializePools()
  {
    All_Event_pool.clear();
    All_Event_pool.resize(nPoolBins());
    isInitialized = true;
  }

  int where_pool(float vz, float v0Centr) const
  {
    float CentBinWidth = 100.0 / numOfCentBins; // = 2.5

    int iy = static_cast<int>(std::floor(v0Centr / CentBinWidth));
    if (iy < 0)
      iy = 0;
    if (iy >= numOfCentBins)
      iy = numOfCentBins - 1;

    int ix = static_cast<int>(std::floor((vz - Vz_low) / Vz_step));
    if (ix < 0)
      ix = 0;
    if (ix >= numOfVertexZBins)
      ix = numOfVertexZBins - 1;

    int bin = ix + numOfVertexZBins * iy;
    return bin;
  }

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
    mBBparamsDe[5] = settingBetheBlochParams->get("De", "resolution");

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
  bool selectionPIDKaon(const Ttrack& candidate)
  {
    if (abs(candidate.dcaXY()) > settingCutHadDCAxyMin || abs(candidate.dcaZ()) > settingCutHadDCAzMin)
      return false;

    auto tpcNSigmaKa = candidate.tpcNSigmaKa();
    mQaRegistry.fill(HIST("h2NsigmaHadTPC_preselection"), candidate.tpcInnerParam(), tpcNSigmaKa);
    if (std::abs(candidate.pt()) < settingCutHadptMin || std::abs(candidate.pt()) > settingCutHadptMax)
      return false;

    if (candidate.hasTOF() && candidate.tpcInnerParam() >= settingCutPinMinTOFHad) {
      auto tofNSigmaKa = candidate.tofNSigmaKa();

      mQaRegistry.fill(HIST("h2NsigmaHadTOF_preselection"), candidate.pt(), tofNSigmaKa);
      if (std::abs(tofNSigmaKa) > settingCutNsigmaTOFHad) {
        return false;
      }
      if (std::abs(tpcNSigmaKa) > settingCutNsigmaTPCHad) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaKa);
      mQaRegistry.fill(HIST("h2NsigmaHadTOF"), candidate.sign() * candidate.pt(), tofNSigmaKa);
      return true;
    } else if (candidate.tpcInnerParam() < settingCutPinMinTOFHad) {
      if (std::abs(tpcNSigmaKa) > settingCutNsigmaTPCHad) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaKa);
      return true;
    }
    return false;
  }

  template <typename Ttrack>
  bool selectionPIDPion(const Ttrack& candidate)
  {
    if (candidate.tpcChi2NCl() > settingCutChi2tpcHighPion || candidate.tpcChi2NCl() < settingCutChi2tpcLowPion || candidate.itsChi2NCl() > settingCutChi2NClITSPion)
      return false;
    if (abs(candidate.dcaXY()) > settingCutHadDCAxyMin || abs(candidate.dcaZ()) > settingCutHadDCAzMin)
      return false;

    auto tpcNSigmaPi = candidate.tpcNSigmaPi();
    mQaRegistry.fill(HIST("h2NsigmaHadTPC_preselection"), candidate.tpcInnerParam(), tpcNSigmaPi);
    if (std::abs(candidate.pt()) < settingCutHadptMin || std::abs(candidate.pt()) > settingCutHadptMax)
      return false;
    // reject protons
    if (std::abs(candidate.tpcNSigmaPr()) < settingCutNsigTPCPrMin)
      return false;
    mQaRegistry.fill(HIST("h2NsigmaHadPrTPC"), candidate.tpcNSigmaPr());
    if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < settingCutNsigTOFPrMin)
      return false;
    mQaRegistry.fill(HIST("h2NsigmaHadPrTOF"), candidate.tofNSigmaPr());

    if (candidate.hasTOF() && candidate.tpcInnerParam() >= settingCutPinMinTOFHad) {
      auto tofNSigmaPi = candidate.tofNSigmaPi();
      auto combNsigma = std::sqrt(tofNSigmaPi * tofNSigmaPi + tpcNSigmaPi * tpcNSigmaPi);

      mQaRegistry.fill(HIST("h2NsigmaHadTOF_preselection"), candidate.pt(), tofNSigmaPi);
      // if (combNsigma > settingCutNsigmaTOFTPCHad) {
      //   return false;
      // }
      if (std::abs(tofNSigmaPi) > settingCutNsigmaTOFHad) {
        return false;
      }
      if (std::abs(tpcNSigmaPi) > settingCutNsigmaTPCHad) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPi);
      mQaRegistry.fill(HIST("h2NsigmaHadTOF"), candidate.sign() * candidate.pt(), tofNSigmaPi);
      mQaRegistry.fill(HIST("h2NsigmaHadComb"), candidate.sign() * candidate.pt(), combNsigma);
      return true;
    } else if (candidate.tpcInnerParam() < settingCutPinMinTOFHad) {
      if (std::abs(tpcNSigmaPi) > settingCutNsigmaTPCHad) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPi);
      return true;
    }
    return false;
  }

  template <typename Ttrack>
  bool selectionPIDHadron(const Ttrack& candidate)
  {
    bool PID = false; 
    if (settingHadPDGCode == PDG_t::kPiPlus) {
      PID = selectionPIDPion(candidate);
      MassHad = o2::constants::physics::MassPiPlus;
    } else if (settingHadPDGCode == PDG_t::kKPlus) {
      PID = selectionPIDKaon(candidate);
      MassHad = o2::constants::physics::MassKPlus;
    } else {
      LOG(info) << "invalid PDG code";
    }
    return PID;
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
    mQaRegistry.fill(HIST("h2dEdx"), candidate.sign() * tpcInnerParam, candidate.tpcSignal());

    float DeDCAxyMin = 0.015 + 0.0305 / TMath::Power(candidate.pt(), 1.1);
    if (abs(candidate.dcaXY()) > DeDCAxyMin || abs(candidate.dcaXY()) > settingCutDeDCAzMin)
      return false;

    if (std::abs(tpcInnerParam) < settingCutPinMinDe) {
      return false;
    }
    float tpcNSigmaDe;
    if (settingUseBBcomputeDeNsigma) {
      tpcNSigmaDe = computeNSigmaDe(candidate);
    } else {
      tpcNSigmaDe = candidate.tpcNSigmaDe();
    }

    mQaRegistry.fill(HIST("h2NsigmaNuTPC_preselection"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
    mQaRegistry.fill(HIST("h2NsigmaNuTPC_preselecComp"), candidate.sign() * candidate.pt(), candidate.tpcNSigmaDe());
    if (std::abs(candidate.pt()) < settingCutDeptMin || std::abs(candidate.pt()) > settingCutDeptMax)
      return false;
    if (candidate.hasTOF() && candidate.tpcInnerParam() > settingCutPinMinTOFITSDe) {
      auto tofNSigmaDe = candidate.tofNSigmaDe();
      auto combNsigma = std::sqrt(tofNSigmaDe * tofNSigmaDe + tpcNSigmaDe * tpcNSigmaDe);
      if (combNsigma > settingCutNsigmaTOFTPCDe) {
        return false;
      }
      mQaRegistry.fill(HIST("h2dEdxNucandidates"), candidate.sign() * tpcInnerParam, candidate.tpcSignal());
      mQaRegistry.fill(HIST("h2NsigmaNuComb"), candidate.sign() * candidate.pt(), combNsigma);
      mQaRegistry.fill(HIST("h2NsigmaNuTPC"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
      mQaRegistry.fill(HIST("h2NsigmaNuTOF"), candidate.sign() * candidate.pt(), tofNSigmaDe);
      return true;
    } else if (candidate.tpcInnerParam() <= settingCutPinMinTOFITSDe) {
      if (std::abs(tpcNSigmaDe) > settingCutNsigmaTPCDe) {
        return false;
      }
      o2::aod::ITSResponse mResponseITS;
      auto itsnSigmaDe = mResponseITS.nSigmaITS<o2::track::PID::Deuteron>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
      mQaRegistry.fill(HIST("h2NSigmaNuITS_preselection"), candidate.sign() * candidate.pt(), itsnSigmaDe);
      if (std::abs(itsnSigmaDe) > settingCutNsigmaITSDe) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaNuTPC"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
      mQaRegistry.fill(HIST("h2NSigmaNuITS"), candidate.sign() * candidate.pt(), itsnSigmaDe);
      // mQaRegistry.fill(HIST("h2NsigmaNuComb"), candidate.sign() * candidate.pt(), combNsigma);
      mQaRegistry.fill(HIST("h2dEdxNucandidates"), candidate.sign() * tpcInnerParam, candidate.tpcSignal());
      return true;
    }
    return false;
  }

  float averageClusterSizeCosl(uint32_t itsClusterSizes, float eta)
  {
    float average = 0;
    int nclusters = 0;
    const float cosl = 1. / std::cosh(eta);
    const int nlayerITS = 7;

    for (int layer = 0; layer < nlayerITS; layer++) {
      if ((itsClusterSizes >> (layer * 4)) & 0xf) {
        nclusters++;
        average += (itsClusterSizes >> (layer * 4)) & 0xf;
      }
    }
    if (nclusters == 0) {
      return 0;
    }
    return average * cosl / nclusters;
  };

  bool selectionPIDHyper(const aod::DataHypCandsWColl::iterator& V0Hyper)
  {
    mQaRegistry.fill(HIST("hHe3P_preselected"), V0Hyper.tpcMomHe());
    float averClusSizeHe = averageClusterSizeCosl(V0Hyper.itsClusterSizesHe(), V0Hyper.etaHe3());
    if (averClusSizeHe <= settingCutAverClsSizeHe) {
      return false;
    }
    if (V0Hyper.tpcChi2He() <= settingCutTPCChi2He) {
      return false;
    }
    mQaRegistry.fill(HIST("hHe3P"), V0Hyper.tpcMomHe());
    mQaRegistry.fill(HIST("hHe3TPCnsigma"), V0Hyper.ptHe3(), V0Hyper.nSigmaHe());

    return true;
  }

  // ==================================================================================================================

  template <typename Ttrack, typename Tcollisions, typename Ttracks>
  bool fillCandidateInfo(const Ttrack& trackDe, const Ttrack& trackHad, const CollBracket& collBracket, const Tcollisions& collisions, HadNucandidate& hadNucand, const Ttracks& /*trackTable*/, bool isMixedEvent)
  {
    const int numCoordinates = 3;
    if (!isMixedEvent) {
      auto trackCovDe = getTrackParCov(trackDe);
      auto trackCovHad = getTrackParCov(trackHad);
      int nCand = 0;
      try {
        nCand = mFitter.process(trackCovDe, trackCovHad);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        mQaRegistry.fill(HIST("hSkipReasons"), 0);
        return false;
      }
      if (nCand == 0) {
        mQaRegistry.fill(HIST("hSkipReasons"), 1);
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
        mQaRegistry.fill(HIST("hSkipReasons"), 2);
        return false;
      }
      hadNucand.collisionID = collIdxMin;
    } else {
      hadNucand.collisionID = collBracket.getMin();
    }

    hadNucand.momNu = std::array{trackDe.px(), trackDe.py(), trackDe.pz()};
    hadNucand.momHad = std::array{trackHad.px(), trackHad.py(), trackHad.pz()};
    float invMass = 0;
    invMass = RecoDecay::m(std::array<std::array<float, 3>, 2>{hadNucand.momNu, hadNucand.momHad}, std::array<float, 2>{static_cast<float>(o2::constants::physics::MassDeuteron), MassHad});
    if (settingCutInvMass > 0 && invMass > settingCutInvMass) {
      mQaRegistry.fill(HIST("hSkipReasons"), 3);
      return false;
    }
    float ptDeHad = std::hypot(hadNucand.momNu[0] + hadNucand.momHad[0], hadNucand.momNu[1] + hadNucand.momHad[1]);
    if (ptDeHad < settingCutPtMinDePi) {
      mQaRegistry.fill(HIST("hSkipReasons"), 4);
      return false;
    }

    hadNucand.signNu = trackDe.sign();
    hadNucand.signHad = trackHad.sign();

    hadNucand.dcaxyNu = trackDe.dcaXY();
    hadNucand.dcaxyHad = trackHad.dcaXY();

    hadNucand.dcazNu = trackDe.dcaZ();
    hadNucand.dcazHad = trackHad.dcaZ();

    hadNucand.tpcSignalNu = trackDe.tpcSignal();
    hadNucand.momNuTPC = trackDe.tpcInnerParam();
    hadNucand.tpcSignalHad = trackHad.tpcSignal();
    hadNucand.momHadTPC = trackHad.tpcInnerParam();

    hadNucand.nTPCClustersNu = trackDe.tpcNClsFound();
    hadNucand.nSigmaNu = computeNSigmaDe(trackDe);
    hadNucand.nSigmaHad = trackHad.tpcNSigmaPi();

    hadNucand.chi2TPCNu = trackDe.tpcChi2NCl();
    hadNucand.chi2TPCHad = trackHad.tpcChi2NCl();

    hadNucand.pidTrkNu = trackDe.pidForTracking();
    hadNucand.pidTrkHad = trackHad.pidForTracking();

    hadNucand.itsClSizeNu = trackDe.itsClusterSizes();
    hadNucand.itsClSizeHad = trackHad.itsClusterSizes();

    hadNucand.nClsItsNu = trackDe.itsNCls();
    hadNucand.nClsItsHad = trackHad.itsNCls();

    hadNucand.sharedClustersNu = trackDe.tpcNClsShared();
    hadNucand.sharedClustersHad = trackHad.tpcNClsShared();

    hadNucand.isBkgUS = trackDe.sign() * trackHad.sign() < 0;
    hadNucand.isBkgEM = isMixedEvent;

    hadNucand.invMass = invMass;

    hadNucand.trackIDNu = trackDe.globalIndex();
    hadNucand.trackIDHad = trackHad.globalIndex();

    if (trackDe.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackDe);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      float tpcInnerParamDe = trackDe.tpcInnerParam();
      hadNucand.massTOFNu = tpcInnerParamDe * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    if (trackHad.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackHad);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      hadNucand.massTOFHad = trackHad.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
    }

    hadNucand.kstar = o2::analysis::femtoWorld::FemtoWorldMath::getkstar(trackHad, MassHad, trackDe, o2::constants::physics::MassDeuteron);
    hadNucand.mT = o2::analysis::femtoWorld::FemtoWorldMath::getmT(trackHad, MassHad, trackDe, o2::constants::physics::MassDeuteron);

    return true;
  }

  template <typename Ttrack>
  bool fillCandidateInfoHyper(const aod::DataHypCandsWColl::iterator& V0Hyper, const Ttrack& trackHad, HadNucandidate& hadHypercand, bool isMixedEvent)
  {
    hadHypercand.collisionID = V0Hyper.collisionId();
    // get hypertriton information
    // constexpr double mHe3 = o2::constants::physics::MassHelium3;
    // constexpr double mPi  = o2::constants::physics::MassPiPlus;
    //  --- He3
    float pxHe3 = V0Hyper.ptHe3() * std::cos(V0Hyper.phiHe3());
    float pyHe3 = V0Hyper.ptHe3() * std::sin(V0Hyper.phiHe3());
    float pzHe3 = V0Hyper.ptHe3() * std::sinh(V0Hyper.etaHe3());
    // float pHe3  = V0Hyper.ptHe3() * std::cosh(V0Hyper.etaHe3());
    // float enHe3 = std::sqrt(pHe3 * pHe3 + mHe3 * mHe3);
    //  --- pi
    float pxPi = V0Hyper.ptPi() * std::cos(V0Hyper.phiPi());
    float pyPi = V0Hyper.ptPi() * std::sin(V0Hyper.phiPi());
    float pzPi = V0Hyper.ptPi() * std::sinh(V0Hyper.etaPi());
    // float pPi  = V0Hyper.ptPi() * std::cosh(V0Hyper.etaPi());
    // float enPi = std::sqrt(pPi * pPi + mPi * mPi);
    //  --- hypertriton
    float px = pxHe3 + pxPi;
    float py = pyHe3 + pyPi;
    float pz = pzHe3 + pzPi;
    hadHypercand.momNu = std::array{px, py, pz};
    hadHypercand.momHad = std::array{trackHad.px(), trackHad.py(), trackHad.pz()};

    float invMass = 0;
    invMass = RecoDecay::m(std::array<std::array<float, 3>, 2>{hadHypercand.momNu, hadHypercand.momHad}, std::array<float, 2>{static_cast<float>(o2::constants::physics::MassHelium3), MassHad});
    if (settingCutInvMass > 0 && invMass > settingCutInvMass) {
      return false;
    }

    hadHypercand.signHad = trackHad.sign();
    if (V0Hyper.isMatter()) {
      hadHypercand.signNu = 1;
    } else {
      hadHypercand.signNu = -1;
    }
    hadHypercand.etaHe3 = V0Hyper.etaHe3();
    hadHypercand.ptHe3 = V0Hyper.ptHe3();
    hadHypercand.dcaxyHad = trackHad.dcaXY();
    hadHypercand.dcazHad = trackHad.dcaZ();
    hadHypercand.tpcSignalHad = trackHad.tpcSignal();
    hadHypercand.tpcSignalNu = V0Hyper.tpcSignalHe();
    hadHypercand.momHadTPC = trackHad.tpcInnerParam();
    hadHypercand.nSigmaHad = trackHad.tpcNSigmaPi();
    hadHypercand.nSigmaNu = V0Hyper.nSigmaHe();
    hadHypercand.chi2TPCHad = trackHad.tpcChi2NCl();
    hadHypercand.chi2TPCNu = V0Hyper.tpcChi2He();
    hadHypercand.pidTrkHad = trackHad.pidForTracking();
    hadHypercand.itsClSizeHad = trackHad.itsClusterSizes();
    hadHypercand.itsClSizeNu = V0Hyper.itsClusterSizesHe();
    hadHypercand.nClsItsHad = trackHad.itsNCls();
    hadHypercand.sharedClustersHad = trackHad.tpcNClsShared();

    hadHypercand.isBkgUS = hadHypercand.signNu * trackHad.sign() < 0;
    hadHypercand.isBkgEM = isMixedEvent;
    hadHypercand.invMass = invMass;

    hadHypercand.trackIDHad = trackHad.globalIndex();

    if (trackHad.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackHad);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      hadHypercand.massTOFHad = trackHad.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    return true;
  }

  template <typename Ttrack>
  void pairTracksSameEvent(const Ttrack& tracks, float /*cent*/)
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
      mQaRegistry.fill(HIST("hSingleNuPt"), track0.pt() * track0.sign());
      mQaRegistry.fill(HIST("hSingleNuPin"), track0.tpcInnerParam() * track0.sign());

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

        if (!selectTrack(track1) || !selectionPIDHadron(track1)) {
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

  template <typename Ttrack, typename Thypers>
  void pairTracksSameEventHyper(const Ttrack& hadTracks, const Thypers& V0Hypers)
  {
    for (const auto& V0Hyper : V0Hypers) {
      if (!selectionPIDHyper(V0Hyper)) {
        continue;
      }
      for (const auto& hadTrack : hadTracks) {

        mQaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);

        if (!selectTrack(hadTrack)) {
          continue;
        }
        mQaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);

        if (!selectionPIDHadron(hadTrack)) {
          continue;
        }
        mQaRegistry.fill(HIST("hTrackSel"), Selections::kPID);

        SVCand pair;
        pair.tr0Idx = V0Hyper.globalIndex();
        pair.tr1Idx = hadTrack.globalIndex();
        const int collIdx = V0Hyper.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        pair.collBracket = collBracket;
        mTrackHypPairs.push_back(pair);
      }
    }
  }

  template <typename T>
  void pairTracksEventMixing(T& DeCands, T& hadCands)
  {
    for (const auto& DeCand : DeCands) {
      if (!selectTrack(DeCand) || !selectionPIDDe(DeCand)) {
        continue;
      }
      for (const auto& hadCand : hadCands) {
        if (!selectTrack(hadCand) || !selectionPIDHadron(hadCand)) {
          continue;
        }

        SVCand trackPair;
        trackPair.tr0Idx = DeCand.globalIndex();
        trackPair.tr1Idx = hadCand.globalIndex();
        const int collIdx = DeCand.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        trackPair.collBracket = collBracket;
        mTrackPairs.push_back(trackPair);
      }
    }
  }

  template <typename T1, typename T2>
  void pairHyperEventMixing(T1& hadCands, T2& hypCands)
  {
    for (const auto& hypCand : hypCands) {
      if (!selectionPIDHyper(hypCand)) {
        continue;
      }
      for (const auto& hadCand : hadCands) {
        if (!selectTrack(hadCand) || !selectionPIDHadron(hadCand)) {
          continue;
        }

        SVCand pair;
        pair.tr0Idx = hypCand.globalIndex();
        pair.tr1Idx = hadCand.globalIndex();
        const int collIdx = hypCand.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        pair.collBracket = collBracket;
        mTrackHypPairs.push_back(pair);
      }
    }
  }

  template <typename Tcoll>
  void fillTable(const HadNucandidate& hadNucand, const Tcoll& collision)
  {
    mOutputDataTable(
      hadNucand.recoPtHad(),
      hadNucand.recoPtNu(),
      hadNucand.momHadTPC,
      hadNucand.momNuTPC,
      hadNucand.trackIDHad,
      hadNucand.trackIDNu);
    if (settingFillMultiplicity) {
      mOutputMultiplicityTable(
        collision.globalIndex(),
        collision.posZ(),
        collision.numContrib(),
        collision.centFT0C(),
        collision.multFT0C());
    }
  }

  template <typename Tcoll>
  void fillTableHyper(const HadNucandidate& hadNucand, const Tcoll& collision)
  {
    mOutputHyperDataTable(
      hadNucand.recoPtNu(),
      hadNucand.recoEtaNu(),
      hadNucand.ptHe3,
      hadNucand.etaHe3,
      hadNucand.recoPhiNu(),
      hadNucand.recoPtHad(),
      hadNucand.recoEtaHad(),
      hadNucand.recoPhiHad(),
      hadNucand.dcaxyHad,
      hadNucand.dcazHad,
      hadNucand.tpcSignalHad,
      hadNucand.tpcSignalNu,
      hadNucand.momHadTPC,
      hadNucand.nSigmaHad,
      hadNucand.nSigmaNu,
      hadNucand.chi2TPCHad,
      hadNucand.chi2TPCNu,
      hadNucand.massTOFHad,
      hadNucand.pidTrkHad,
      hadNucand.itsClSizeHad,
      hadNucand.itsClSizeNu,
      hadNucand.sharedClustersHad,
      hadNucand.trackIDHad,
      hadNucand.isBkgUS,
      hadNucand.isBkgEM);
    if (settingFillMultiplicity) {
      mOutputMultiplicityTable(
        collision.globalIndex(),
        collision.posZ(),
        collision.numContrib(),
        collision.centFT0C(),
        collision.multFT0C());
    }
  }

  void fillHistograms(const HadNucandidate& hadNucand)
  {
    mQaRegistry.fill(HIST("hNuPt"), hadNucand.recoPtNu());
    mQaRegistry.fill(HIST("hHadPt"), hadNucand.recoPtHad());
    mQaRegistry.fill(HIST("hNuPin"), hadNucand.momNuTPC * hadNucand.signNu);
    mQaRegistry.fill(HIST("hHadPin"), hadNucand.momHadTPC * hadNucand.signHad);
    mQaRegistry.fill(HIST("hNuEta"), hadNucand.recoEtaNu());
    mQaRegistry.fill(HIST("hHadEta"), hadNucand.recoEtaHad());
    mQaRegistry.fill(HIST("hNuPhi"), hadNucand.recoPhiNu());
    mQaRegistry.fill(HIST("hHadPhi"), hadNucand.recoPhiHad());
    mQaRegistry.fill(HIST("hNuHadtInvMass"), hadNucand.invMass);
    mQaRegistry.fill(HIST("hdcaxyNu"), hadNucand.dcaxyNu);
    mQaRegistry.fill(HIST("hdcazNu"), hadNucand.dcazNu);
    mQaRegistry.fill(HIST("hdcazNu_min"), (abs(hadNucand.dcazNu) - settingCutDeDCAzMin));
    mQaRegistry.fill(HIST("hNClsNuITS"), hadNucand.nClsItsNu);
    mQaRegistry.fill(HIST("hNClsHadITS"), hadNucand.nClsItsHad);
    mQaRegistry.fill(HIST("hisBkgEM"), hadNucand.isBkgEM);
  }

  template <typename Tcoll>
  void fillKstar(const HadNucandidate& hadNucand, const Tcoll& collision)
  {
    if (hadNucand.isBkgUS == 0) {
      if (hadNucand.recoPtNu() > 0) {
        mQaRegistry.fill(HIST("hkStar_LS_M"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmTVsCent_LS_M"), hadNucand.kstar, hadNucand.mT, collision.centFT0C());
        mQaRegistry.fill(HIST("hkStaVsmT_LS_M"), hadNucand.kstar, hadNucand.mT);
      } else {
        mQaRegistry.fill(HIST("hkStar_LS_A"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmTVsCent_LS_A"), hadNucand.kstar, hadNucand.mT, collision.centFT0C());
        mQaRegistry.fill(HIST("hkStaVsmT_LS_A"), hadNucand.kstar, hadNucand.mT);
      }
    } else {
      if (hadNucand.recoPtNu() > 0) {
        mQaRegistry.fill(HIST("hkStar_US_M"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmTVsCent_US_M"), hadNucand.kstar, hadNucand.mT, collision.centFT0C());
        mQaRegistry.fill(HIST("hkStaVsmT_US_M"), hadNucand.kstar, hadNucand.mT);
      } else {
        mQaRegistry.fill(HIST("hkStar_US_A"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmTVsCent_US_A"), hadNucand.kstar, hadNucand.mT, collision.centFT0C());
        mQaRegistry.fill(HIST("hkStaVsmT_US_A"), hadNucand.kstar, hadNucand.mT);
      }
    }
  }

  // ==================================================================================================================

  template <typename Tcollisions, typename Ttracks>
  void fillPairs(const Tcollisions& collisions, const Ttracks& tracks, const bool isMixedEvent)
  {
    for (const auto& trackPair : mTrackPairs) {

      auto deTrack = tracks.rawIteratorAt(trackPair.tr0Idx);
      auto hadTrack = tracks.rawIteratorAt(trackPair.tr1Idx);
      auto collBracket = trackPair.collBracket;

      HadNucandidate hadNucand;
      if (!fillCandidateInfo(deTrack, hadTrack, collBracket, collisions, hadNucand, tracks, isMixedEvent)) {
        continue;
      }

      auto collision = collisions.rawIteratorAt(hadNucand.collisionID);
      fillKstar(hadNucand, collision);
      fillHistograms(hadNucand);

      if (settingFillTable) {
        fillTable(hadNucand, collision);
      }
    }
  }

  template <typename Tcollisions, typename Ttracks>
  void fillPairsHyper(const Tcollisions& collisions, const Ttracks& hadTracks, const o2::aod::DataHypCandsWColl& V0Hypers, const bool isMixedEvent)
  {
    for (const auto& trackPair : mTrackHypPairs) {

      auto v0hyper = V0Hypers.rawIteratorAt(trackPair.tr0Idx);
      auto hadTrack = hadTracks.rawIteratorAt(trackPair.tr1Idx);
      // auto collBracket = trackPair.collBracket;

      HadNucandidate hadNucand;
      if (!fillCandidateInfoHyper(v0hyper, hadTrack, hadNucand, isMixedEvent)) {
        continue;
      }

      mQaRegistry.fill(HIST("hNuPt"), hadNucand.recoPtNu());
      mQaRegistry.fill(HIST("hHadPt"), hadNucand.recoPtHad());
      mQaRegistry.fill(HIST("hNuEta"), hadNucand.recoEtaNu());
      mQaRegistry.fill(HIST("hHadEta"), hadNucand.recoEtaHad());
      mQaRegistry.fill(HIST("hNuPhi"), hadNucand.recoPhiNu());
      mQaRegistry.fill(HIST("hHadPhi"), hadNucand.recoPhiHad());
      mQaRegistry.fill(HIST("hNuHadtInvMass"), hadNucand.invMass);
      mQaRegistry.fill(HIST("hNClsHadITS"), hadNucand.nClsItsHad);
      mQaRegistry.fill(HIST("hisBkgEM"), hadNucand.isBkgEM);

      auto collision = collisions.rawIteratorAt(hadNucand.collisionID);

      if (settingFillTable) {
        fillTableHyper(hadNucand, collision);
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

      pairTracksSameEvent(trackTableThisCollision, collision.centFT0C());

      if (mTrackPairs.size() == 0) {
        continue;
      }

      fillPairs(collisions, tracks, /*isMixedEvent*/ false);
    }
  }
  PROCESS_SWITCH(HadNucleiFemto, processSameEvent, "Process Same event", false);

  void processSameEventHyper(const CollisionsFull& collisions, const TrackCandidates& hadtracks, o2::aod::DataHypCandsWColl const& V0Hypers, const aod::BCsWithTimestamps& bcs)
  {
    mGoodCollisions.clear();
    mGoodCollisions.resize(collisions.size(), false);
    // LOG(info) << "Number of hyperCandidates read = " << V0Hypers.size();

    for (const auto& collision : collisions) {

      mTrackHypPairs.clear();

      if (!selectCollision</*isMC*/ false>(collision, bcs)) {
        continue;
      }

      mGoodCollisions[collision.globalIndex()] = true;
      const uint64_t collIdx = collision.globalIndex();
      auto trackTableThisCollision = hadtracks.sliceBy(mPerCol, collIdx);
      auto hypdTableThisCollision = V0Hypers.sliceBy(hypPerCol, collIdx);
      trackTableThisCollision.bindExternalIndices(&hadtracks);
      hypdTableThisCollision.bindExternalIndices(&V0Hypers);

      pairTracksSameEventHyper(trackTableThisCollision, hypdTableThisCollision);

      if (mTrackHypPairs.size() == 0) {
        continue;
      }

      fillPairsHyper(collisions, hadtracks, V0Hypers, /*isMixedEvent*/ false);
    }
  }
  PROCESS_SWITCH(HadNucleiFemto, processSameEventHyper, "Process Same event", false);

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
  PROCESS_SWITCH(HadNucleiFemto, processMixedEvent, "Process Mixed event", false);

  /*void processMixedEventHyper(const CollisionsFull& collisions, o2::aod::DataHypCandsWColl const& V0Hypers, const TrackCandidates& hadtracks)
  {
    LOG(debug) << "Processing mixed event for hypertriton";
    mTrackHypPairs.clear();

    for (const auto& [c1, tracks1, c2, V0Hypers2] : hyperPair) {
      if (!c1.sel8() || !c2.sel8()) {
        continue;
      }

      mQaRegistry.fill(HIST("hNcontributor"), c2.numContrib());
      //mQaRegistry.fill(HIST("hCentrality"), c2.centFT0C());
      mQaRegistry.fill(HIST("hVtxZ"), c2.posZ());

      pairHyperEventMixing(tracks1, V0Hypers2);
    }
}
PROCESS_SWITCH(HadNucleiFemto, processMixedEventHyper, "Process Mixed event", false);*/

  void processMixedEventHyperPool(const CollisionsFull& collisions, o2::aod::DataHypCandsWColl const& V0Hypers, const TrackCandidates& hadtracks)
  {
    mTrackHypPairs.clear();
    if (!isInitialized) {
      initializePools();
      LOG(info) << "Initialized event pool with size = " << All_Event_pool.size();
    }
    for (auto const& collision : collisions) {
      if (!collision.sel8()) {
        mQaRegistry.fill(HIST("hSkipReasons"), 0);
        continue;
      }
      mQaRegistry.fill(HIST("hNcontributor"), collision.numContrib());
      mQaRegistry.fill(HIST("hCentrality"), collision.centFT0C());
      mQaRegistry.fill(HIST("hVtxZ"), collision.posZ());

      int poolIndexHad = where_pool(collision.posZ(), collision.centFT0C());
      if (poolIndexHad < 0 || static_cast<size_t>(poolIndexHad) >= All_Event_pool.size()) {
        continue;
      }
      auto& pool = All_Event_pool[poolIndexHad];

      const uint64_t collIdxHad = collision.globalIndex();
      auto trackTableThisCollision = hadtracks.sliceBy(mPerCol, collIdxHad);
      trackTableThisCollision.bindExternalIndices(&hadtracks);

      for (auto const& storedEvent : pool.events) {
        const uint64_t collIdxHyp = storedEvent.collisionId;
        if (settingSaferME) {
          if (static_cast<int64_t>(collIdxHyp) > collisions.size()) {
            mQaRegistry.fill(HIST("hSkipReasons"), 4);
            continue;
          }
        }

        auto hypdTablepreviousCollision = V0Hypers.sliceBy(hypPerCol, collIdxHyp);
        hypdTablepreviousCollision.bindExternalIndices(&V0Hypers);
        if (hypdTablepreviousCollision.size() == 0) {
          mQaRegistry.fill(HIST("hSkipReasons"), 1);
          continue;
        }

        auto firstHyp = hypdTablepreviousCollision.iteratorAt(0);
        int poolIndexHyp = where_pool(firstHyp.zPrimVtx(), firstHyp.centralityFT0C());
        if (poolIndexHyp != poolIndexHad) {
          mQaRegistry.fill(HIST("hSkipReasons"), 2);
          continue;
        }
        mQaRegistry.fill(HIST("hNHypsPerPrevColl"), collIdxHyp, hypdTablepreviousCollision.size());

        pairHyperEventMixing(trackTableThisCollision, hypdTablepreviousCollision);
      }

      if (static_cast<int>(pool.events.size()) >= settingNoMixedEvents) {
        pool.events.pop_front();
      }
      pool.events.push_back({collIdxHad});
    }
    fillPairsHyper(collisions, hadtracks, V0Hypers, /*isMixedEvent*/ true);
  }
  PROCESS_SWITCH(HadNucleiFemto, processMixedEventHyperPool, "Process Mixed event", false);
};

WorkflowSpec defineDataProcessing(const ConfigContext& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HadNucleiFemto>(cfgc)};
}
