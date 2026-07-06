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
#include "PWGLF/DataModel/LFHypernucleiTables.h"
#include "PWGLF/Utils/svPoolCreator.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
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

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/BetheBlochAleph.h>
#include <MathUtils/Primitive2D.h>
#include <ReconstructionDataFormats/PID.h>

#include <THn.h>
#include <TPDGCode.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
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
constexpr double betheBlochDefault[2][6]{
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09}};
static const std::vector<std::string> betheBlochParticleNames{"De", "He3"};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static constexpr std::array<float, 9> tmpRadiiTPC{{85.f, 105.f, 125.f, 145.f, 165.f, 185.f, 205.f, 225.f, 245.f}};
constexpr int DeuteronPDG = o2::constants::physics::Pdg::kDeuteron;
constexpr int He3PDG = o2::constants::physics::Pdg::kHelium3;
constexpr float He3RigidityMin = 0.8f;
constexpr int He3TPCNClsFoundMin = 110;
constexpr float He3TPCChi2NClMin = 0.5f;
constexpr float He3TPCNSigmaMax = 3.0f;
constexpr float He3ITSNSigmaMin = -1.5f;

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
  float dcaPair = -10.f;

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
  float nSigmaTPCHadPi = -10.f;
  float nSigmaTPCHadKa = -10.f;
  float nSigmaTPCHadPr = -10.f;
  float nSigmaTOFHadPi = -10.f;
  float nSigmaTOFHadKa = -10.f;
  float nSigmaTOFHadPr = -10.f;
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

  // Particle species configuration
  Configurable<int> settingNuPDGCode{"settingNuPDGCode", static_cast<int>(DeuteronPDG), "Nucleus - PDG code"};
  Configurable<int> settingHadPDGCode{"settingHadPDGCode", 211, "Hadron - PDG code"};
  // Event selection and mixing configuration
  Configurable<float> settingCutVertex{"settingCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<int> settingNoMixedEvents{"settingNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> settingEnableBkgUS{"settingEnableBkgUS", false, "Enable US background"};
  Configurable<bool> settingSaveUSandLS{"settingSaveUSandLS", true, "Save All Pairs"};
  // Common track-quality cuts
  Configurable<float> settingCutEta{"settingCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> settingCutNCls{"settingCutNCls", 5.0f, "Minimum ITS Ncluster for tracks"};
  Configurable<float> settingCutChi2tpcLow{"settingCutChi2tpcLow", 0.5f, "Low cut on TPC chi2"};
  Configurable<float> settingCutChi2tpcHigh{"settingCutChi2tpcHigh", 4.f, "High cut on TPC chi2"};
  Configurable<float> settingCutChi2NClITS{"settingCutChi2NClITS", 36.f, "Maximum ITS Chi2 for tracks"};
  // Hadron purity and PID cuts
  Configurable<float> settingCutPinMinTOFHad{"settingCutPinMinTOFHad", 0.5f, "Minimum Pin to apply the TOF cut on hadrons"};
  Configurable<float> settingCutNsigmaTPCHad{"settingCutNsigmaTPCHad", 3.0f, "Value of the TPC Nsigma cut on Had"};
  Configurable<float> settingCutNsigmaTOFHad{"settingCutNsigmaTOFHad", 3.0f, "Value of the hsdron TOF Nsigma cut"};
  Configurable<float> settingCutNsigTPCPrMin{"settingCutNsigTPCPrMin", 3.0f, "Minimum TPC Pr Nsigma cut for rejection"};
  Configurable<float> settingCutNsigTPCPiMin{"settingCutNsigTPCPiMin", 3.0f, "Minimum TPC Pi Nsigma cut for rejection"};
  Configurable<float> settingCutNsigTOFPrMin{"settingCutNsigTOFPrMin", 3.0f, "Minimum TOF Pr Nsigma cut for rejection"};
  Configurable<float> settingCutNsigTOFPiMin{"settingCutNsigTOFPiMin", 3.0f, "Minimum TOF Pi Nsigma cut for rejection"};
  Configurable<float> settingHadptMin{"settingHadptMin", 0.14f, "Minimum pT for the reference pion track cuts"};
  Configurable<float> settingHadptMax{"settingHadptMax", 2.5f, "Maximum pT for the reference pion track cuts"};
  Configurable<int> settingPionITSInnerBarrelMin{"settingPionITSInnerBarrelMin", 3, "Minimum ITS inner barrel clusters for the reference pion track cuts"};
  Configurable<int> settingPionITSNClsMin{"settingPionITSNClsMin", 7, "Minimum ITS clusters for the reference pion track cuts"};
  Configurable<int> settingPionTPCNClsFoundMin{"settingPionTPCNClsFoundMin", 80, "Minimum found TPC clusters for the reference pion track cuts"};
  Configurable<int> settingPionTPCCrossedRowsMin{"settingPionTPCCrossedRowsMin", 90, "Minimum crossed TPC rows for the reference pion track cuts"};
  Configurable<float> settingPionDCAxyOffset{"settingPionDCAxyOffset", 0.004f, "DCAxy offset for the reference pion track cuts"};
  Configurable<float> settingPionDCAxyPtCoeff{"settingPionDCAxyPtCoeff", 0.013f, "DCAxy 1/pT coefficient for the reference pion track cuts"};
  Configurable<float> settingPionDCAzOffset{"settingPionDCAzOffset", 0.004f, "DCAz offset for the reference pion track cuts"};
  Configurable<float> settingPionDCAzPtCoeff{"settingPionDCAzPtCoeff", 0.013f, "DCAz 1/pT coefficient for the reference pion track cuts"};
  Configurable<float> settingPionMomCombMin{"settingPionMomCombMin", 0.5f, "Minimum momentum to use combined TPC+TOF PID for reference pions"};
  Configurable<float> settingPionTPCNsigMax{"settingPionTPCNsigMax", 3.0f, "Maximum TPC n-sigma for reference pions below the TOF threshold"};
  Configurable<float> settingPionCombNsigMax{"settingPionCombNsigMax", 3.0f, "Maximum combined TPC+TOF n-sigma for reference pions"};
  // Deuteron purity and PID cuts
  Configurable<float> settingCutPinMinDe{"settingCutPinMinDe", 0.0f, "Minimum Pin for De"};
  Configurable<float> settingCutClSizeItsDe{"settingCutClSizeItsDe", 4.0f, "Minimum ITS cluster size for De"};
  Configurable<float> settingCutDeptMin{"settingCutDeptMin", 0.6f, "Minimum PT cut on De"};
  Configurable<float> settingCutDeptMax{"settingCutDeptMax", 1.6f, "Maximum PT cut on De"};
  Configurable<float> settingCutPinMinTOFITSDe{"settingCutPinMinTOFITSDe", 1.2f, "Minimum p to apply the TOF ITS cut on De"};
  Configurable<float> settingCutNsigmaTPCDe{"settingCutNsigmaTPCDe", 2.5f, "Value of the TPC Nsigma cut on De"};
  Configurable<float> settingCutNsigmaITSDe{"settingCutNsigmaITSDe", 2.5f, "Value of the ITD Nsigma cut on De"};
  Configurable<float> settingCutNsigmaTOFTPCDe{"settingCutNsigmaTOFTPCDe", 2.5f, "Value of the De TOF TPC combNsigma cut"};
  Configurable<bool> settingReqSingleNsig{"settingReqSingleNsig", false, "If true, also require individual TPC and TOF n-sigma cuts in branches using combined TPC+TOF PID"};
  Configurable<bool> settingUseProtonMassForKstarMt{"settingUseProtonMassForKstarMt", false, "If true, use proton mass instead of deuteron mass for kstar and mT"};
  Configurable<bool> settingEnableClosePairRejection{"settingEnableClosePairRejection", false, "Enable close pair rejection for nucleus-hadron track pairs"};
  Configurable<float> settingClosePairDeltaPhiMax{"settingClosePairDeltaPhiMax", 0.01f, "Maximum delta phi star for close pair rejection"};
  Configurable<float> settingClosePairDeltaEtaMax{"settingClosePairDeltaEtaMax", 0.01f, "Maximum delta eta for close pair rejection"};
  Configurable<int> settingClosePairRadiusMode{"settingClosePairRadiusMode", 1, "Close pair rejection mode: 0 = PV, 1 = average phi star, 2 = specific TPC radius"};
  Configurable<float> settingClosePairSpecificRadius{"settingClosePairSpecificRadius", 85.f, "TPC radius in cm used when close pair rejection mode is 2"};
  // Hypertriton-specific cuts
  Configurable<float> settingCutTPCChi2He{"settingCutTPCChi2He", 0.0f, "Minimum tpcChi2He for Hyper He3"};
  Configurable<float> settingCutAverClsSizeHe{"settingCutAverClsSizeHe", 0.0f, "Minimum averClusSizeHe for Hyper He3"};
  // Output and QA controls
  Configurable<bool> settingFillTable{"settingFillTable", false, "Enable table filling"};
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

  Configurable<LabeledArray<double>> settingBetheBlochParams{"settingBetheBlochParams", {betheBlochDefault[0], 2, 6, betheBlochParticleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for the selected nucleus"};
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

  std::array<float, 6> mBBparamsNucleus;
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
    {// Event-level
     {"hVtxZ", "Vertex distribution in Z;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
     {"hNcontributor", "Number of primary vertex contributor", {HistType::kTH1F, {{2000, 0.0f, 2000.0f}}}},
     {"hCentrality", "Centrality", {HistType::kTH1F, {{100, 0.0f, 100.0f}}}},
     {"hSkipReasons", "Why storedEvent skipped;Reason;Counts", {HistType::kTH1F, {{5, -0.5, 4.5}}}},
     {"hEvents", "; Events;", {HistType::kTH1F, {{3, -0.5, 2.5}}}},

     // Candidate topology and kinematics
     {"hTrackSel", "Accepted hadron tracks", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
     {"hTrackSelNu", "Accepted nucleus tracks", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
     {"hNuPairFlow", "Nucleus pair-building flow;step;counts", {HistType::kTH1F, {{3, -0.5, 2.5}}}},

     {"hdcaxyNu", ";DCA_{xy} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hdcazNu", ";DCA_{z} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hNClsNuITS", ";N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
     {"hNuPt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{280, -7.0f, 7.0f}}}},
     {"hSingleNuPt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{280, -7.0f, 7.0f}}}},
     {"hNuPin", "#it{p} distribution; #it{p} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
     {"hSingleNuPin", "#it{p} distribution; #it{p} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
     {"hNuEta", "eta distribution; #eta(Nu)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hNuPhi", "phi distribution; phi(Nu)", {HistType::kTH1F, {{600, -4.0f, 4.0f}}}},

     {"hdcaxyHad", ";DCA_{xy} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hdcazHad", ";DCA_{z} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hNClsHadITS", ";N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
     {"hHadPt", "Pt distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{280, -7.0f, 7.0f}}}},
     {"hSingleHadPt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{280, -7.0f, 7.0f}}}},
     {"hHadPin", "P distribution; #it{p} (GeV/#it{c})", {HistType::kTH1F, {{120, -4.0f, 4.0f}}}},
     {"hHadEta", "eta distribution; #eta(had)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hHadPhi", "phi distribution; phi(had)", {HistType::kTH1F, {{600, -4.0f, 4.0f}}}},
     {"h2CPRBefore", "Close pair rejection before cut; #Delta#eta; #Delta#phi^{*}", {HistType::kTH2F, {{160, -2.0f, 2.0f}, {160, -3.2f, 3.2f}}}},
     {"h2CPRAfter", "Close pair rejection after cut; #Delta#eta; #Delta#phi^{*}", {HistType::kTH2F, {{160, -2.0f, 2.0f}, {160, -3.2f, 3.2f}}}},

     // dE/dx
     {"h2dEdxNucandidates", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"h2dEdxHadcandidates", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"h2dEdx", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},

     // Nucleus PID
     {"h2NsigmaNuTPC", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaNuComb", "NsigmaNu TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"h2NsigmaNuTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaNuTPC_preselecComp", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NSigmaNuITS_preselection", "NsigmaNu ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} Nu", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {120, -3.0f, 3.0f}}}},
     {"h2NSigmaNuITS", "NsigmaNu ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} Nu", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {120, -3.0f, 3.0f}}}},
     {"h2NsigmaNuTOF", "NsigmaNu TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaNuTOF_preselection", "NsigmaNu TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},

     // Hadron PID
     {"h2NsigmaHadComb", "NsigmaHad TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(had)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"h2NsigmaHadTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaHadTOF", "NsigmaHad TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadTOF_preselection", "NsigmaHad TOF distribution; #iit{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaHadComb_preselection", "NsigmaHad TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(had)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"h2NsigmaHadPrTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadPiTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(pi)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadKaTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(K)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadPrTOF", "NsigmaHad TOF distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadPiTOF", "NsigmaHad TOF distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(pi)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadKaTOF", "NsigmaHad TOF distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(K)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},

     // Purity
     {"purity/h2NsigmaNuTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NsigmaNuTPC_preselecComp", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NSigmaNuITS_preselection", "NsigmaNu ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} Nu", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {120, -3.0f, 3.0f}}}},
     {"purity/h2NsigmaNuTOF_preselection", "NsigmaNu TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NsigmaNuComb_preselection", "NsigmaNu TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"purity/h2NsigmaHadTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NsigmaHadTOF_preselection", "NsigmaHad TOF distribution; #iit{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NsigmaHadComb_preselection", "NsigmaHad TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(had)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},

     // Hypertriton
     {"hHe3TPCnsigma", "NsigmaHe3 TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(He3)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"hHe3P", "Pin distribution; p (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},
     {"hHe3P_preselected", "Pin distribution_preselected; p (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},

     // Correlation observables
     {"hkStar_LS_M", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_LS_A", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_US_M", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_US_A", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStaVsmT_LS_M", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsmT_LS_A", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsmT_US_M", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsmT_US_A", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsCent_LS_M", ";kStar (GeV/c);Centrality", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsCent_LS_A", ";kStar (GeV/c);Centrality", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsCent_US_M", ";kStar (GeV/c);Centrality", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsCent_US_A", ";kStar (GeV/c);Centrality", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {100, 0.0f, 100.0f}}}},
     {"hNuHadtInvMass", "; M(Nu + had) (GeV/#it{c}^{2})", {HistType::kTH1F, {{500, 2.5f, 4.5f}}}},

     // Mixed-event
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
    const char* betheBlochLabel = nucleusBetheBlochLabel();
    for (int i = 0; i < numParticles; i++) {
      mBBparamsNucleus[i] = settingBetheBlochParams->get(betheBlochLabel, Form("p%i", i));
    }
    mBBparamsNucleus[5] = settingBetheBlochParams->get(betheBlochLabel, "resolution");

    std::vector<std::string> selectionLabels = {"All", "Track selection", "PID"};
    for (int i = 0; i < Selections::kAll; i++) {
      mQaRegistry.get<TH1>(HIST("hTrackSel"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i].c_str());
      mQaRegistry.get<TH1>(HIST("hTrackSelNu"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i].c_str());
    }

    std::vector<std::string> eventsLabels = {"All", "Selected", "Zorro selected events"};
    for (int i = 0; i < Selections::kAll; i++) {
      mQaRegistry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(i + 1, eventsLabels[i].c_str());
    }
  }

  void initCCDB(const aod::BCsWithTimestamps::iterator& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    if (settingSkimmedProcessing) {
      mZorro.initCCDB(mCcdb.service, bc.runNumber(), bc.timestamp(), useHelium3Nucleus() ? "fHe" : "fDe");
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
    const int minTPCNClsCrossedRows = 70;
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
  bool selectTrackPion(const Ttrack& candidate)
  {
    if (std::abs(candidate.eta()) > settingCutEta) {
      return false;
    }

    const float absPt = std::abs(candidate.pt());
    if (absPt < settingHadptMin || absPt > settingHadptMax) {
      return false;
    }

    if (candidate.itsNClsInnerBarrel() < settingPionITSInnerBarrelMin ||
        candidate.itsNCls() < settingPionITSNClsMin ||
        candidate.tpcNClsFound() < settingPionTPCNClsFoundMin ||
        candidate.tpcNClsCrossedRows() < settingPionTPCCrossedRowsMin) {
      return false;
    }

    if (absPt <= 0.f) {
      return false;
    }

    const float pionDCAxyMax = settingPionDCAxyOffset + settingPionDCAxyPtCoeff / absPt;
    const float pionDCAzMax = settingPionDCAzOffset + settingPionDCAzPtCoeff / absPt;
    if (std::abs(candidate.dcaXY()) > pionDCAxyMax || std::abs(candidate.dcaZ()) > pionDCAzMax) {
      return false;
    }

    return true;
  }

  template <typename Ttrack>
  bool selectTrackHadron(const Ttrack& candidate)
  {
    if (settingHadPDGCode.value == static_cast<int>(PDG_t::kProton)) {
      return selectTrackProton(candidate);
    }
    if (settingHadPDGCode.value == static_cast<int>(PDG_t::kPiPlus)) {
      return selectTrackPion(candidate);
    }
    return selectTrack(candidate);
  }

  template <typename Ttrack>
  bool selectTrackProton(const Ttrack& candidate)
  {
    constexpr float protonEtaMax = 0.8f;
    constexpr int protonTPCNClsFoundMin = 90;
    constexpr int protonTPCCrossedRowsMin = 80;
    constexpr float protonDCAzMax = 0.2f;

    if (std::abs(candidate.eta()) >= protonEtaMax) {
      return false;
    }

    const float absPt = std::abs(candidate.pt());
    if (absPt <= 0.f) {
      return false;
    }

    if (candidate.tpcNClsFound() <= protonTPCNClsFoundMin ||
        candidate.tpcNClsCrossedRows() <= protonTPCCrossedRowsMin) {
      return false;
    }

    const float prDCAxyMax = 105.e-3f + 30.5e-3f / std::pow(absPt, 1.1f);
    if (std::abs(candidate.dcaXY()) >= prDCAxyMax || std::abs(candidate.dcaZ()) >= protonDCAzMax) {
      return false;
    }

    return true;
  }

  template <typename Ttrack>
  bool selectTrackDe(const Ttrack& candidate)
  {
    if (std::abs(candidate.eta()) > settingCutEta) {
      return false;
    }

    constexpr int minTPCNClsFound = 110;
    constexpr int minTPCNClsCrossedRows = 100;
    constexpr float minTPCCrossedRowsOverFound = 0.f;
    constexpr int maxTPCNClsShared = 160;
    constexpr float maxSharedTPCFraction = 1.f;
    constexpr int minITSNClsInnerBarrel = 1;
    const float tpcCrossedRowsOverFound = candidate.tpcNClsFound() > 0 ? static_cast<float>(candidate.tpcNClsCrossedRows()) / candidate.tpcNClsFound() : 0.f;

    if (candidate.tpcNClsFound() < minTPCNClsFound ||
        candidate.tpcNClsCrossedRows() < minTPCNClsCrossedRows ||
        tpcCrossedRowsOverFound < minTPCCrossedRowsOverFound ||
        candidate.tpcNClsShared() > maxTPCNClsShared ||
        candidate.tpcFractionSharedCls() > maxSharedTPCFraction ||
        candidate.itsNCls() < settingCutNCls ||
        candidate.itsNClsInnerBarrel() < minITSNClsInnerBarrel) {
      return false;
    }

    return true;
  }

  bool useDeuteronNucleus() const
  {
    return settingNuPDGCode.value == DeuteronPDG;
  }

  bool useHelium3Nucleus() const
  {
    return settingNuPDGCode.value == He3PDG;
  }

  const char* nucleusBetheBlochLabel() const
  {
    return useHelium3Nucleus() ? "He3" : "De";
  }

  float nucleusChargeFactor() const
  {
    return useHelium3Nucleus() ? 2.f : 1.f;
  }

  float nucleusMass() const
  {
    if (useHelium3Nucleus()) {
      return static_cast<float>(o2::constants::physics::MassHelium3);
    }
    return static_cast<float>(o2::constants::physics::MassDeuteron);
  }

  template <typename Ttrack>
  bool selectTrackHe3(const Ttrack& candidate)
  {
    if (std::abs(candidate.eta()) > settingCutEta) {
      return false;
    }

    constexpr int minTPCNClsCrossedRows = 70;
    constexpr float crossedRowsToFindableRatio = 0.8f;
    if (candidate.itsNCls() < settingCutNCls ||
        candidate.tpcNClsFound() < He3TPCNClsFoundMin ||
        candidate.tpcNClsCrossedRows() < minTPCNClsCrossedRows ||
        candidate.tpcNClsCrossedRows() < crossedRowsToFindableRatio * candidate.tpcNClsFindable() ||
        candidate.tpcChi2NCl() > settingCutChi2tpcHigh ||
        candidate.tpcChi2NCl() < He3TPCChi2NClMin ||
        candidate.itsChi2NCl() > settingCutChi2NClITS) {
      return false;
    }

    return true;
  }

  template <typename Ttrack>
  bool selectTrackNu(const Ttrack& candidate)
  {
    if (useHelium3Nucleus()) {
      return selectTrackHe3(candidate);
    }
    if (useDeuteronNucleus()) {
      return selectTrackDe(candidate);
    }
    LOG(info) << "invalid nucleus PDG code";
    return false;
  }

  void fillNucleusTrackSelection(const Selections selection)
  {
    mQaRegistry.fill(HIST("hTrackSelNu"), selection);
  }

  void fillNucleusPairFlow(const int step)
  {
    mQaRegistry.fill(HIST("hNuPairFlow"), step);
  }

  template <typename Ttrack>
  float phiAtSpecificRadiiTPC(const Ttrack& track, float radius) const
  {
    const float absPt = std::abs(track.pt());
    if (absPt <= 0.f) {
      return 999.f;
    }
    const float arg = 0.3f * static_cast<float>(track.sign()) * 0.1f * mDbz * radius * 0.01f / (2.f * absPt);
    if (std::fabs(arg) >= 1.f) {
      return 999.f;
    }
    return track.phi() - std::asin(arg);
  }

  float wrapDeltaPhi(float dphi) const
  {
    return std::atan2(std::sin(dphi), std::cos(dphi));
  }

  template <typename Ttrack1, typename Ttrack2>
  float averagePhiStar(const Ttrack1& track1, const Ttrack2& track2) const
  {
    constexpr float invalidPhiStar = 999.f;
    float dPhiAvg = 0.f;
    int meaningfulEntries = 0;
    for (const auto& radius : tmpRadiiTPC) {
      const float phi1 = phiAtSpecificRadiiTPC(track1, radius);
      const float phi2 = phiAtSpecificRadiiTPC(track2, radius);
      if (phi1 == invalidPhiStar || phi2 == invalidPhiStar) {
        continue;
      }
      dPhiAvg += wrapDeltaPhi(phi1 - phi2);
      meaningfulEntries++;
    }
    if (meaningfulEntries == 0) {
      return invalidPhiStar;
    }
    return dPhiAvg / static_cast<float>(meaningfulEntries);
  }

  template <typename Ttrack1, typename Ttrack2>
  bool isClosePair(const Ttrack1& track1, const Ttrack2& track2)
  {
    constexpr int closePairRadiusModePv = 0;
    constexpr int closePairRadiusModeSpecificTpc = 2;
    constexpr float invalidPhiStar = 999.f;
    if (!settingEnableClosePairRejection.value) {
      return false;
    }
    if (track1.sign() != track2.sign()) {
      return false;
    }

    const float deta = track1.eta() - track2.eta();
    const float dphiAtPV = wrapDeltaPhi(track1.phi() - track2.phi());
    const float dphiAtSpecificRadius = wrapDeltaPhi(phiAtSpecificRadiiTPC(track1, settingClosePairSpecificRadius.value) - phiAtSpecificRadiiTPC(track2, settingClosePairSpecificRadius.value));
    const float dphiAvg = averagePhiStar(track1, track2);

    float dphiToCut = dphiAvg;
    if (settingClosePairRadiusMode.value == closePairRadiusModePv) {
      dphiToCut = dphiAtPV;
    } else if (settingClosePairRadiusMode.value == closePairRadiusModeSpecificTpc) {
      dphiToCut = dphiAtSpecificRadius;
    }

    if (dphiToCut == invalidPhiStar) {
      return false;
    }

    mQaRegistry.fill(HIST("h2CPRBefore"), deta, dphiToCut);
    const bool isRejected = std::pow(dphiToCut, 2.f) / std::pow(settingClosePairDeltaPhiMax.value, 2.f) +
                              std::pow(deta, 2.f) / std::pow(settingClosePairDeltaEtaMax.value, 2.f) <
                            1.f;
    if (!isRejected) {
      mQaRegistry.fill(HIST("h2CPRAfter"), deta, dphiToCut);
    }
    return isRejected;
  }

  template <typename Ttrack>
  bool selectionPIDProton(const Ttrack& candidate)
  {
    constexpr float protonPtMin = 0.5f;
    constexpr float protonPtMax = 3.0f;
    constexpr float protonPCombMin = 0.75f;
    constexpr float protonTPCNsigmaMax = 3.0f;
    constexpr float protonCombNsigmaMax = 3.0f;

    const float tpcNSigmaPr = candidate.tpcNSigmaPr();
    mQaRegistry.fill(HIST("h2NsigmaHadTPC_preselection"), candidate.sign() * candidate.tpcInnerParam(), tpcNSigmaPr);

    if (std::abs(candidate.pt()) <= protonPtMin || std::abs(candidate.pt()) >= protonPtMax) {
      return false;
    }

    const float absPin = std::abs(candidate.tpcInnerParam());
    if (absPin < protonPCombMin) {
      if (std::abs(tpcNSigmaPr) > protonTPCNsigmaMax) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPr);
      mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
      return true;
    }

    if (!candidate.hasTOF()) {
      return false;
    }

    const float tofNSigmaPr = candidate.tofNSigmaPr();
    const float combNsigma = std::sqrt(tpcNSigmaPr * tpcNSigmaPr + tofNSigmaPr * tofNSigmaPr);
    mQaRegistry.fill(HIST("h2NsigmaHadTOF_preselection"), candidate.sign() * candidate.pt(), tofNSigmaPr);
    mQaRegistry.fill(HIST("h2NsigmaHadComb_preselection"), candidate.sign() * candidate.pt(), combNsigma);
    if (combNsigma > protonCombNsigmaMax) {
      return false;
    }
    if (settingReqSingleNsig.value &&
        (std::abs(tpcNSigmaPr) > protonCombNsigmaMax || std::abs(tofNSigmaPr) > protonCombNsigmaMax)) {
      return false;
    }

    mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPr);
    mQaRegistry.fill(HIST("h2NsigmaHadTOF"), candidate.sign() * candidate.pt(), tofNSigmaPr);
    mQaRegistry.fill(HIST("h2NsigmaHadComb"), candidate.sign() * candidate.pt(), combNsigma);
    mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
    return true;
  }

  template <typename Ttrack>
  bool selectionPIDKaon(const Ttrack& candidate)
  {
    auto tpcNSigmaKa = candidate.tpcNSigmaKa();
    float DeDCAxyMin = 0.004 + (0.013 / candidate.pt());
    float DeDCAzMin = 0.004 + (0.013 / candidate.pt());
    if (std::abs(candidate.dcaXY()) > DeDCAxyMin || std::abs(candidate.dcaZ()) > DeDCAzMin)
      return false;

    mQaRegistry.fill(HIST("h2NsigmaHadTPC_preselection"), candidate.tpcInnerParam(), tpcNSigmaKa);
    if (std::abs(candidate.pt()) < settingHadptMin || std::abs(candidate.pt()) > settingHadptMax)
      return false;

    // reject protons and pions
    if (std::abs(candidate.tpcNSigmaPr()) < settingCutNsigTPCPrMin || std::abs(candidate.tpcNSigmaPi()) < settingCutNsigTPCPiMin)
      return false;
    mQaRegistry.fill(HIST("h2NsigmaHadPrTPC"), candidate.tpcNSigmaPr());
    mQaRegistry.fill(HIST("h2NsigmaHadPiTPC"), candidate.tpcNSigmaPi());
    if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < settingCutNsigTOFPrMin)
      return false;
    if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < settingCutNsigTOFPiMin)
      return false;
    mQaRegistry.fill(HIST("h2NsigmaHadPrTOF"), candidate.tofNSigmaPr());
    mQaRegistry.fill(HIST("h2NsigmaHadPiTOF"), candidate.tofNSigmaPi());
    // rejection end

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
      mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
      return true;
    } else if (candidate.tpcInnerParam() < settingCutPinMinTOFHad) {
      if (std::abs(tpcNSigmaKa) > settingCutNsigmaTPCHad) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaKa);
      mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
      return true;
    }
    return false;
  }

  template <typename Ttrack>
  bool selectionPIDPion(const Ttrack& candidate)
  {
    const float tpcNSigmaPi = candidate.tpcNSigmaPi();
    const float absP = std::abs(candidate.p());
    mQaRegistry.fill(HIST("h2NsigmaHadTPC_preselection"), candidate.sign() * candidate.tpcInnerParam(), tpcNSigmaPi);

    if (absP <= settingPionMomCombMin) {
      if (std::abs(tpcNSigmaPi) > settingPionTPCNsigMax) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPi);
      mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
      return true;
    }

    if (!candidate.hasTOF()) {
      return false;
    }

    const float tofNSigmaPi = candidate.tofNSigmaPi();
    const float combNsigma = std::sqrt(tofNSigmaPi * tofNSigmaPi + tpcNSigmaPi * tpcNSigmaPi);
    mQaRegistry.fill(HIST("h2NsigmaHadTOF_preselection"), candidate.sign() * candidate.pt(), tofNSigmaPi);
    mQaRegistry.fill(HIST("h2NsigmaHadComb_preselection"), candidate.sign() * candidate.pt(), combNsigma);
    if (combNsigma > settingPionCombNsigMax) {
      return false;
    }
    if (settingReqSingleNsig.value &&
        (std::abs(tpcNSigmaPi) > settingPionCombNsigMax || std::abs(tofNSigmaPi) > settingPionCombNsigMax)) {
      return false;
    }

    mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPi);
    mQaRegistry.fill(HIST("h2NsigmaHadTOF"), candidate.sign() * candidate.pt(), tofNSigmaPi);
    mQaRegistry.fill(HIST("h2NsigmaHadComb"), candidate.sign() * candidate.pt(), combNsigma);
    mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
    return true;
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
    } else if (settingHadPDGCode == PDG_t::kProton) {
      PID = selectionPIDProton(candidate);
      MassHad = o2::constants::physics::MassProton;
    } else {
      LOG(info) << "invalid PDG code";
    }
    return PID;
  }

  template <typename Ttrack>
  float getHadronTPCNSigma(const Ttrack& candidate) const
  {
    if (settingHadPDGCode.value == static_cast<int>(PDG_t::kPiPlus)) {
      return candidate.tpcNSigmaPi();
    }
    if (settingHadPDGCode.value == static_cast<int>(PDG_t::kKPlus)) {
      return candidate.tpcNSigmaKa();
    }
    if (settingHadPDGCode.value == static_cast<int>(PDG_t::kProton)) {
      return candidate.tpcNSigmaPr();
    }
    return -10.f;
  }

  template <typename Ttrack>
  float computeNSigmaDe(const Ttrack& candidate)
  {
    float expTPCSignal = o2::common::BetheBlochAleph(static_cast<float>(candidate.tpcInnerParam() / constants::physics::MassDeuteron), mBBparamsNucleus[0], mBBparamsNucleus[1], mBBparamsNucleus[2], mBBparamsNucleus[3], mBBparamsNucleus[4]);
    double resoTPC{expTPCSignal * mBBparamsNucleus[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename Ttrack>
  float correctedTPCInnerParamHe3(const Ttrack& candidate) const
  {
    const bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    return (heliumPID && settingCompensatePIDinTracking.value) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
  }

  template <typename Ttrack>
  float computeNSigmaHe3(const Ttrack& candidate)
  {
    const float correctedTPCinnerParam = correctedTPCInnerParamHe3(candidate);
    float expTPCSignal = o2::common::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), mBBparamsNucleus[0], mBBparamsNucleus[1], mBBparamsNucleus[2], mBBparamsNucleus[3], mBBparamsNucleus[4]);
    double resoTPC{expTPCSignal * mBBparamsNucleus[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename Ttrack>
  bool selectionPIDDe(const Ttrack& candidate)
  {
    float tpcInnerParam = candidate.tpcInnerParam();
    mQaRegistry.fill(HIST("h2dEdx"), candidate.sign() * tpcInnerParam, candidate.tpcSignal());

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
    const float absPt = std::abs(candidate.pt());
    if (absPt <= 0.f) {
      return false;
    }
    const float deDCAxyMax = 0.004f + 0.013f / absPt;
    const float deDCAzMax = 0.004f + 0.013f / absPt;
    if (std::abs(candidate.dcaXY()) > deDCAxyMax || std::abs(candidate.dcaZ()) > deDCAzMax)
      return false;

    if (candidate.hasTOF() && candidate.tpcInnerParam() > settingCutPinMinTOFITSDe) {
      auto tofNSigmaDe = candidate.tofNSigmaDe();
      auto combNsigma = std::sqrt(tofNSigmaDe * tofNSigmaDe + tpcNSigmaDe * tpcNSigmaDe);
      mQaRegistry.fill(HIST("h2NsigmaNuTOF_preselection"), candidate.sign() * candidate.pt(), tofNSigmaDe);
      if (combNsigma > settingCutNsigmaTOFTPCDe) {
        return false;
      }
      if (settingReqSingleNsig.value &&
          (std::abs(tpcNSigmaDe) > settingCutNsigmaTOFTPCDe || std::abs(tofNSigmaDe) > settingCutNsigmaTOFTPCDe)) {
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
      o2::aod::ITSResponse itsResponse;
      auto itsnSigmaDe = itsResponse.nSigmaITS<o2::track::PID::Deuteron>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
      mQaRegistry.fill(HIST("h2NSigmaNuITS_preselection"), candidate.sign() * candidate.pt(), itsnSigmaDe);
      if (std::abs(itsnSigmaDe) > settingCutNsigmaITSDe) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaNuTPC"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
      mQaRegistry.fill(HIST("h2NSigmaNuITS"), candidate.sign() * candidate.pt(), itsnSigmaDe);
      mQaRegistry.fill(HIST("h2dEdxNucandidates"), candidate.sign() * tpcInnerParam, candidate.tpcSignal());
      return true;
    }
    return false;
  }

  template <typename Ttrack>
  bool selectionPIDHe3(const Ttrack& candidate)
  {
    const float correctedTPCinnerParam = correctedTPCInnerParamHe3(candidate);
    mQaRegistry.fill(HIST("h2dEdx"), candidate.sign() * correctedTPCinnerParam, candidate.tpcSignal());

    if (correctedTPCinnerParam < He3RigidityMin) {
      return false;
    }

    const float nSigmaHe3 = computeNSigmaHe3(candidate);
    mQaRegistry.fill(HIST("h2NsigmaNuTPC_preselection"), candidate.sign() * 2.f * candidate.pt(), nSigmaHe3);
    if (std::abs(nSigmaHe3) > He3TPCNSigmaMax) {
      return false;
    }

    o2::aod::ITSResponse itsResponse;
    const float itsNsigmaHe3 = itsResponse.nSigmaITS<o2::track::PID::Helium3>(candidate.itsClusterSizes(), 2.f * candidate.p(), candidate.eta());
    mQaRegistry.fill(HIST("h2NSigmaNuITS_preselection"), candidate.sign() * 2.f * candidate.pt(), itsNsigmaHe3);
    if (itsNsigmaHe3 < He3ITSNSigmaMin) {
      return false;
    }

    mQaRegistry.fill(HIST("h2dEdxNucandidates"), candidate.sign() * correctedTPCinnerParam, candidate.tpcSignal());
    mQaRegistry.fill(HIST("h2NsigmaNuTPC"), candidate.sign() * 2.f * candidate.pt(), nSigmaHe3);
    mQaRegistry.fill(HIST("h2NSigmaNuITS"), candidate.sign() * 2.f * candidate.pt(), itsNsigmaHe3);
    return true;
  }

  template <typename Ttrack>
  bool selectionPIDNu(const Ttrack& candidate)
  {
    if (useHelium3Nucleus()) {
      return selectionPIDHe3(candidate);
    }
    if (useDeuteronNucleus()) {
      return selectionPIDDe(candidate);
    }
    return false;
  }

  template <typename Ttrack>
  float getNucleusTPCNSigma(const Ttrack& candidate)
  {
    if (useHelium3Nucleus()) {
      return computeNSigmaHe3(candidate);
    }
    return computeNSigmaDe(candidate);
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

  float computePairKstar(const std::array<float, 3>& momHad, const float massHad, const std::array<float, 3>& momNu, const float massNu) const
  {
    const ROOT::Math::PxPyPzMVector vecHad(momHad[0], momHad[1], momHad[2], massHad);
    const ROOT::Math::PxPyPzMVector vecNu(momNu[0], momNu[1], momNu[2], massNu);
    const ROOT::Math::PxPyPzMVector trackSum = vecHad + vecNu;

    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());

    ROOT::Math::PxPyPzMVector partHadCMS(vecHad);
    ROOT::Math::PxPyPzMVector partNuCMS(vecNu);

    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    partHadCMS = boostPRF(partHadCMS);
    partNuCMS = boostPRF(partNuCMS);

    const ROOT::Math::PxPyPzMVector trackRelK = partHadCMS - partNuCMS;
    return 0.5f * trackRelK.P();
  }

  float computePairMT(const std::array<float, 3>& momHad, const float massHad, const std::array<float, 3>& momNu, const float massNu) const
  {
    const ROOT::Math::PxPyPzMVector vecHad(momHad[0], momHad[1], momHad[2], massHad);
    const ROOT::Math::PxPyPzMVector vecNu(momNu[0], momNu[1], momNu[2], massNu);
    const ROOT::Math::PxPyPzMVector trackSum = vecHad + vecNu;
    const float kT = 0.5f * trackSum.Pt();
    return std::sqrt(kT * kT + std::pow(0.5f * (massHad + massNu), 2.f));
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
      hadNucand.dcaPair = std::sqrt(std::abs(mFitter.getChi2AtPCACandidate()));

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

    const float nuChargeFactor = nucleusChargeFactor();
    hadNucand.momNu = std::array{trackDe.px(), trackDe.py(), trackDe.pz()};
    for (auto i = 0u; i < hadNucand.momNu.size(); ++i) {
      hadNucand.momNu[i] *= nuChargeFactor;
    }
    hadNucand.momHad = std::array{trackHad.px(), trackHad.py(), trackHad.pz()};
    float invMass = 0;
    invMass = RecoDecay::m(std::array<std::array<float, 3>, 2>{hadNucand.momNu, hadNucand.momHad}, std::array<float, 2>{nucleusMass(), MassHad});

    hadNucand.signNu = trackDe.sign();
    hadNucand.signHad = trackHad.sign();

    hadNucand.dcaxyNu = trackDe.dcaXY();
    hadNucand.dcaxyHad = trackHad.dcaXY();

    hadNucand.dcazNu = trackDe.dcaZ();
    hadNucand.dcazHad = trackHad.dcaZ();

    hadNucand.tpcSignalNu = trackDe.tpcSignal();
    hadNucand.momNuTPC = useHelium3Nucleus() ? correctedTPCInnerParamHe3(trackDe) : trackDe.tpcInnerParam();
    hadNucand.tpcSignalHad = trackHad.tpcSignal();
    hadNucand.momHadTPC = trackHad.tpcInnerParam();

    hadNucand.nTPCClustersNu = trackDe.tpcNClsFound();
    hadNucand.nSigmaNu = getNucleusTPCNSigma(trackDe);
    hadNucand.nSigmaHad = getHadronTPCNSigma(trackHad);
    hadNucand.nSigmaTPCHadPi = trackHad.tpcNSigmaPi();
    hadNucand.nSigmaTPCHadKa = trackHad.tpcNSigmaKa();
    hadNucand.nSigmaTPCHadPr = trackHad.tpcNSigmaPr();
    hadNucand.nSigmaTOFHadPi = trackHad.tofNSigmaPi();
    hadNucand.nSigmaTOFHadKa = trackHad.tofNSigmaKa();
    hadNucand.nSigmaTOFHadPr = trackHad.tofNSigmaPr();

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
      hadNucand.massTOFNu = hadNucand.momNuTPC * nuChargeFactor * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    if (trackHad.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackHad);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      hadNucand.massTOFHad = trackHad.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
    }

    const float massLightNucleusForKstarMt = useHelium3Nucleus() ? static_cast<float>(o2::constants::physics::MassHelium3) : (settingUseProtonMassForKstarMt ? static_cast<float>(o2::constants::physics::MassProton) : static_cast<float>(o2::constants::physics::MassDeuteron));
    hadNucand.kstar = computePairKstar(hadNucand.momHad, MassHad, hadNucand.momNu, massLightNucleusForKstarMt);
    hadNucand.mT = computePairMT(hadNucand.momHad, MassHad, hadNucand.momNu, massLightNucleusForKstarMt);

    return true;
  }

  template <typename Ttrack>
  void fillCandidateInfoHyper(const aod::DataHypCandsWColl::iterator& V0Hyper, const Ttrack& trackHad, HadNucandidate& hadHypercand, bool isMixedEvent)
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
    hadHypercand.nSigmaHad = getHadronTPCNSigma(trackHad);
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
  }

  template <typename Ttrack>
  void pairTracksSameEvent(const Ttrack& tracks, float /*cent*/)
  {
    // LOG(info) << "Number of tracks: " << tracks.size();
    for (const auto& track0 : tracks) {

      mQaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);
      fillNucleusTrackSelection(Selections::kNoCuts);

      if (!selectTrackNu(track0)) {
        continue;
      }
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);
      fillNucleusTrackSelection(Selections::kTrackCuts);

      if (!selectionPIDNu(track0)) {
        continue;
      }
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kPID);
      fillNucleusTrackSelection(Selections::kPID);
      mQaRegistry.fill(HIST("hSingleNuPt"), track0.pt() * track0.sign() * nucleusChargeFactor());
      mQaRegistry.fill(HIST("hSingleNuPin"), (useHelium3Nucleus() ? correctedTPCInnerParamHe3(track0) : track0.tpcInnerParam()) * track0.sign());
      fillNucleusPairFlow(0);

      bool hasHadronSelected = false;
      bool hasStoredPair = false;

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

        if (!selectTrackHadron(track1) || !selectionPIDHadron(track1)) {
          continue;
        }
        hasHadronSelected = true;
        if (isClosePair(track0, track1)) {
          continue;
        }

        SVCand trackPair;
        trackPair.tr0Idx = track0.globalIndex();
        trackPair.tr1Idx = track1.globalIndex();
        const int collIdx = track0.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        trackPair.collBracket = collBracket;
        mTrackPairs.push_back(trackPair);
        hasStoredPair = true;
      }

      if (hasHadronSelected) {
        fillNucleusPairFlow(1);
      }
      if (hasStoredPair) {
        fillNucleusPairFlow(2);
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

        if (!selectTrackHadron(hadTrack)) {
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
      if (!selectTrackNu(DeCand) || !selectionPIDNu(DeCand)) {
        continue;
      }
      for (const auto& hadCand : hadCands) {
        if (!selectTrackHadron(hadCand) || !selectionPIDHadron(hadCand)) {
          continue;
        }
        if (isClosePair(DeCand, hadCand)) {
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

  template <typename Tcoll>
  void fillTable(const HadNucandidate& hadNucand, const Tcoll& collision)
  {
    mOutputDataTable(
      hadNucand.recoPtNu(),
      hadNucand.recoEtaNu(),
      hadNucand.recoPhiNu(),
      hadNucand.recoPtHad(),
      hadNucand.recoEtaHad(),
      hadNucand.recoPhiHad(),
      hadNucand.dcaxyNu,
      hadNucand.dcazNu,
      hadNucand.dcaxyHad,
      hadNucand.dcazHad,
      hadNucand.dcaPair,
      hadNucand.tpcSignalNu,
      hadNucand.momNuTPC,
      hadNucand.tpcSignalHad,
      hadNucand.momHadTPC,
      hadNucand.nTPCClustersNu,
      hadNucand.nSigmaNu,
      hadNucand.nSigmaTPCHadPi,
      hadNucand.nSigmaTPCHadKa,
      hadNucand.nSigmaTPCHadPr,
      hadNucand.nSigmaTOFHadPi,
      hadNucand.nSigmaTOFHadKa,
      hadNucand.nSigmaTOFHadPr,
      hadNucand.chi2TPCNu,
      hadNucand.chi2TPCHad,
      hadNucand.massTOFNu,
      hadNucand.massTOFHad,
      hadNucand.pidTrkNu,
      hadNucand.pidTrkHad,
      hadNucand.itsClSizeNu,
      hadNucand.itsClSizeHad,
      hadNucand.sharedClustersNu,
      hadNucand.sharedClustersHad);
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
    mQaRegistry.fill(HIST("hdcaxyHad"), hadNucand.dcaxyHad);
    mQaRegistry.fill(HIST("hdcazHad"), hadNucand.dcazHad);
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
        mQaRegistry.fill(HIST("hkStaVsmT_LS_M"), hadNucand.kstar, hadNucand.mT);
        mQaRegistry.fill(HIST("hkStaVsCent_LS_M"), hadNucand.kstar, collision.centFT0C());
      } else {
        mQaRegistry.fill(HIST("hkStar_LS_A"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmT_LS_A"), hadNucand.kstar, hadNucand.mT);
        mQaRegistry.fill(HIST("hkStaVsCent_LS_A"), hadNucand.kstar, collision.centFT0C());
      }
    } else {
      if (hadNucand.recoPtNu() > 0) {
        mQaRegistry.fill(HIST("hkStar_US_M"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmT_US_M"), hadNucand.kstar, hadNucand.mT);
        mQaRegistry.fill(HIST("hkStaVsCent_US_M"), hadNucand.kstar, collision.centFT0C());
      } else {
        mQaRegistry.fill(HIST("hkStar_US_A"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmT_US_A"), hadNucand.kstar, hadNucand.mT);
        mQaRegistry.fill(HIST("hkStaVsCent_US_A"), hadNucand.kstar, collision.centFT0C());
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
      fillCandidateInfoHyper(v0hyper, hadTrack, hadNucand, isMixedEvent);

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

  void processPurity(const CollisionsFull& collisions, const TrackCandidates& tracks, const aod::BCsWithTimestamps& bcs)
  {
    for (const auto& collision : collisions) {
      if (!selectCollision</*isMC*/ false>(collision, bcs)) {
        continue;
      }

      const uint64_t collIdx = collision.globalIndex();
      auto trackTableThisCollision = tracks.sliceBy(mPerCol, collIdx);
      trackTableThisCollision.bindExternalIndices(&tracks);

      for (const auto& track : trackTableThisCollision) {
        const bool passTrackHad = selectTrackHadron(track);
        const bool passTrackNu = selectTrackNu(track);

        mQaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);
        if (passTrackHad) {
          mQaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);
        }

        fillNucleusTrackSelection(Selections::kNoCuts);
        if (passTrackNu) {
          fillNucleusTrackSelection(Selections::kTrackCuts);
        }

        if (passTrackHad && settingHadPDGCode == PDG_t::kPiPlus) {
          const float tpcNSigmaHad = track.tpcNSigmaPi();
          mQaRegistry.fill(HIST("purity/h2NsigmaHadTPC_preselection"), track.sign() * track.pt(), tpcNSigmaHad);
          if (track.hasTOF() && std::abs(track.p()) > settingPionMomCombMin) {
            const float tofNSigmaHad = track.tofNSigmaPi();
            const float combNsigmaHad = std::sqrt(tofNSigmaHad * tofNSigmaHad + tpcNSigmaHad * tpcNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadTOF_preselection"), track.sign() * track.pt(), tofNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadComb_preselection"), track.sign() * track.pt(), combNsigmaHad);
          }
        } else if (passTrackHad && settingHadPDGCode == PDG_t::kKPlus) {
          const float tpcNSigmaHad = track.tpcNSigmaKa();
          mQaRegistry.fill(HIST("purity/h2NsigmaHadTPC_preselection"), track.sign() * track.pt(), tpcNSigmaHad);
          if (track.hasTOF() && track.tpcInnerParam() >= settingCutPinMinTOFHad) {
            const float tofNSigmaHad = track.tofNSigmaKa();
            const float combNsigmaHad = std::sqrt(tofNSigmaHad * tofNSigmaHad + tpcNSigmaHad * tpcNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadTOF_preselection"), track.sign() * track.pt(), tofNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadComb_preselection"), track.sign() * track.pt(), combNsigmaHad);
          }
        } else if (passTrackHad && settingHadPDGCode == PDG_t::kProton) {
          constexpr float protonPCombMin = 0.75f;
          const float tpcNSigmaHad = track.tpcNSigmaPr();
          mQaRegistry.fill(HIST("purity/h2NsigmaHadTPC_preselection"), track.sign() * track.pt(), tpcNSigmaHad);
          if (track.hasTOF() && std::abs(track.tpcInnerParam()) >= protonPCombMin) {
            const float tofNSigmaHad = track.tofNSigmaPr();
            const float combNsigmaHad = std::sqrt(tofNSigmaHad * tofNSigmaHad + tpcNSigmaHad * tpcNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadTOF_preselection"), track.sign() * track.pt(), tofNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadComb_preselection"), track.sign() * track.pt(), combNsigmaHad);
          }
        }

        if (passTrackNu && useDeuteronNucleus()) {
          const float tpcNSigmaDe = settingUseBBcomputeDeNsigma ? computeNSigmaDe(track) : track.tpcNSigmaDe();
          mQaRegistry.fill(HIST("purity/h2NsigmaNuTPC_preselection"), track.sign() * track.pt(), tpcNSigmaDe);
          mQaRegistry.fill(HIST("purity/h2NsigmaNuTPC_preselecComp"), track.sign() * track.pt(), track.tpcNSigmaDe());
          if (track.hasTOF() && track.tpcInnerParam() > settingCutPinMinTOFITSDe) {
            const float tofNSigmaDe = track.tofNSigmaDe();
            const float combNsigmaDe = std::sqrt(tofNSigmaDe * tofNSigmaDe + tpcNSigmaDe * tpcNSigmaDe);
            mQaRegistry.fill(HIST("purity/h2NsigmaNuTOF_preselection"), track.sign() * track.pt(), tofNSigmaDe);
            mQaRegistry.fill(HIST("purity/h2NsigmaNuComb_preselection"), track.sign() * track.pt(), combNsigmaDe);
          } else if (track.tpcInnerParam() <= settingCutPinMinTOFITSDe) {
            o2::aod::ITSResponse itsResponse;
            const float itsNSigmaDe = itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track.itsClusterSizes(), track.p(), track.eta());
            mQaRegistry.fill(HIST("purity/h2NSigmaNuITS_preselection"), track.sign() * track.pt(), itsNSigmaDe);
          }
        } else if (passTrackNu && useHelium3Nucleus()) {
          const float tpcNSigmaHe3 = computeNSigmaHe3(track);
          const float signedPtHe3 = track.sign() * 2.f * track.pt();
          mQaRegistry.fill(HIST("purity/h2NsigmaNuTPC_preselection"), signedPtHe3, tpcNSigmaHe3);
          o2::aod::ITSResponse itsResponse;
          const float itsNSigmaHe3 = itsResponse.nSigmaITS<o2::track::PID::Helium3>(track.itsClusterSizes(), 2.f * track.p(), track.eta());
          mQaRegistry.fill(HIST("purity/h2NSigmaNuITS_preselection"), signedPtHe3, itsNSigmaHe3);
        }

        const bool isHadronSelected = passTrackHad && selectionPIDHadron(track);
        const bool isNucleusSelected = passTrackNu && selectionPIDNu(track);
        if (!isHadronSelected && !isNucleusSelected) {
          continue;
        }

        if (isHadronSelected) {
          mQaRegistry.fill(HIST("hTrackSel"), Selections::kPID);
          mQaRegistry.fill(HIST("hSingleHadPt"), track.pt() * track.sign());
        }

        if (isNucleusSelected) {
          fillNucleusTrackSelection(Selections::kPID);
          mQaRegistry.fill(HIST("hSingleNuPt"), track.pt() * track.sign() * nucleusChargeFactor());
          mQaRegistry.fill(HIST("hSingleNuPin"), (useHelium3Nucleus() ? correctedTPCInnerParamHe3(track) : track.tpcInnerParam()) * track.sign());
        }
      }
    }
  }
  PROCESS_SWITCH(HadNucleiFemto, processPurity, "Process for hadron and nucleus purity QA", false);
};

WorkflowSpec defineDataProcessing(const ConfigContext& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HadNucleiFemto>(cfgc)};
}
