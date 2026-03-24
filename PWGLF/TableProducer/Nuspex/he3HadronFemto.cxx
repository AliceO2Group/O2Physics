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

/// \file he3HadronFemto.cxx
/// \brief Femto analysis task for He3-hadron correlation
/// \author Your Name (your.email@cern.ch)
/// \since April 2025

#include "PWGCF/FemtoWorld/Core/FemtoWorldMath.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFhe3HadronTables.h"
#include "PWGLF/Utils/nucleiUtils.h"
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
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/BetheBlochAleph.h"
#include "ReconstructionDataFormats/Track.h"

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator> // std::prev
#include <string>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using CollBracket = o2::math_utils::Bracket<int>;

using McIter = aod::McParticles::iterator;
using CollBracket = o2::math_utils::Bracket<int>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::FT0Mults>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::FT0Mults>;
using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::TOFSignal, aod::TOFEvTime>;
using TrackCandidatesMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::TOFSignal, aod::TOFEvTime, aod::McTrackLabels>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

constexpr int Li4PDG = o2::constants::physics::Pdg::kLithium4;
constexpr int H3LPDG = o2::constants::physics::Pdg::kHyperTriton;
constexpr int ProtonPDG = PDG_t::kProton;
constexpr int He3PDG = o2::constants::physics::Pdg::kHelium3;
constexpr float CommonInite = 0.0f;

float getkstar(const float pt1, const float eta1, const float phi1, const float mass1, const float z1,
               const float pt2, const float eta2, const float phi2, const float mass2, const float z2)
{
  const ROOT::Math::PtEtaPhiMVector vecpart1(pt1 * z1, eta1, phi1, mass1);
  const ROOT::Math::PtEtaPhiMVector vecpart2(pt2 * z2, eta2, phi2, mass2);
  const ROOT::Math::PtEtaPhiMVector trackSum = vecpart1 + vecpart2;

  const float beta = trackSum.Beta();
  const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
  const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
  const float betaz = beta * std::cos(trackSum.Theta());

  ROOT::Math::PxPyPzMVector PartOneCMS(vecpart1);
  ROOT::Math::PxPyPzMVector PartTwoCMS(vecpart2);

  const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
  PartOneCMS = boostPRF(PartOneCMS);
  PartTwoCMS = boostPRF(PartTwoCMS);

  const ROOT::Math::PxPyPzMVector trackRelK = PartOneCMS - PartTwoCMS;
  return 0.5 * trackRelK.P();
}

enum Selections {
  kNoCuts = 0,
  kTrackCuts,
  kPID,
  kAll
};

enum Flags {
  kBothPrimaries = BIT(0),
  kBothFromLi4 = BIT(1),
  kBothFromHypertriton = BIT(2),
  kMixedPair = BIT(3), // a primary and one from Li4/hypertriton/material/other decays (or any other combination)
};

enum Species {
  kHe3 = 0,
  kHad,
  kAllSpecies
};

enum ParticleFlags {
  kPhysicalPrimary = BIT(0), // primary particle
  kFromLi4 = BIT(1),         // from Li4 decay
  kFromHypertriton = BIT(2), // from hypertriton decay
  kFromMaterial = BIT(3),    // from material
  kFromOtherDecays = BIT(4), // from other decays
};

std::array<float, 3> kDCAxyResolutionParams[static_cast<int>(Species::kAllSpecies)] = {
  {0.0118, 0.6889, 0.0017}, // He3
  {0.0032, 0.5206, 0.0012}  // Pr
};
std::array<float, 3> kDCAzResolutionParams[static_cast<int>(Species::kAllSpecies)] = {
  {0.1014, 1.7512, 0.0024}, // He3
  {0.0021, 1.1122, 0.0021}  // Pr
};

std::array<float, 2> kHePidTrkParams = {0.1593, -0.0445};

} // namespace

struct He3HadCandidate {

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
  float dcaxyHe3 = -10.f;
  float dcazHe3 = -10.f;
  float dcaxyHad = -10.f;
  float dcazHad = -10.f;
  float dcaPair = -10.f; // DCA between the two tracks

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
  float nSigmaTPCHadPi = -10.f;
  float nSigmaTPCHadKa = -10.f;
  float nSigmaTPCHadPr = -10.f;
  float nSigmaTOFHadPi = -10.f;
  float nSigmaTOFHadKa = -10.f;
  float nSigmaTOFHadPr = -10.f;
  uint32_t pidtrkHe3 = 0xFFFFF; // PID in tracking
  uint32_t pidtrkHad = 0xFFFFF;
  float massTOFHe3 = -10;
  float massTOFHad = -10;
  uint32_t itsClSizeHe3 = 0u;
  uint32_t itsClSizeHad = 0u;

  uint8_t nclsITSHe3 = 0u;
  uint8_t nclsITSHad = 0u;
  float chi2nclITSHe3 = -10.f;
  float chi2nclITSHad = -10.f;

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

  uint8_t flagsHe3 = 0; // flags for He3
  uint8_t flagsHad = 0; // flags for hadron
  uint8_t flags = 0;    // flags for the pair

  // collision information
  int32_t collisionID = 0;
};

struct he3HadronFemto {

  Produces<aod::he3HadronTable> outputDataTable;
  Produces<aod::he3HadronTableMC> outputMcTable;
  Produces<aod::he3HadronMult> outputMultiplicityTable;
  Produces<aod::he3HadronQa> outputQaTable;

  // Selections
  Configurable<int> settingHadPDGCode{"settingHadPDGCode", 211, "Hadron - PDG code"};

  Configurable<float> settingCutVertex{"settingCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> settingCutRigidityMinHe3{"settingCutRigidityMinHe3", 0.8f, "Minimum rigidity for He3"};
  Configurable<float> settingCutEta{"settingCutEta", 0.9f, "Eta cut on daughter track"};
  Configurable<float> settingCutDCAxy{"settingCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> settingCutDCAz{"settingCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> settingCutNClsTPC{"settingCutNClsTPC", 90, "number of TPC clusters for a generic track"};
  Configurable<float> settingCutNClsTPCHe3{"settingCutNClsTPCHe3", 110.0f, "number of TPC clusters for a He3 track"};
  Configurable<float> settingCutChi2tpcLow{"settingCutChi2tpcLow", 0.f, "Low cut on TPC chi2"};
  Configurable<float> settingCutChi2tpcLowHe3{"settingCutChi2tpcLowHe3", 0.5f, "Low cut on TPC chi2 for He3"};
  Configurable<float> settingCutInvMass{"settingCutInvMass", 0.0f, "Invariant mass upper limit"};
  Configurable<float> settingCutPtMinhe3Had{"settingCutPtMinhe3Had", 0.0f, "Minimum PT cut on he3Had4"};
  Configurable<float> settingCutClSizeItsHe3{"settingCutClSizeItsHe3", 4.0f, "Minimum ITS cluster size for He3"};
  Configurable<float> settingCutNCls{"settingCutNCls", 5.0f, "Minimum ITS Ncluster for tracks"};
  Configurable<float> settingCutChi2NClITS{"settingCutChi2NClITS", 36.f, "Maximum ITS Chi2 for tracks"};
  Configurable<float> settingCutNsigmaDcaXy{"settingCutNsigmaDcaXy", 3.0f, "Value of the DCA xy Nsigma cut"};
  Configurable<float> settingCutNsigmaDcaZ{"settingCutNsigmaDcaZ", 3.0f, "Value of the DCA z Nsigma cut"};
  Configurable<float> settingCutNsigmaTPC{"settingCutNsigmaTPC", 3.0f, "Value of the TPC Nsigma cut"};
  Configurable<float> settingCutNsigmaITSHad{"settingCutNsigmaITSHad", -2.f, "Value of the ITS Nsigma cutfor Had"};
  Configurable<float> settingCutNsigmaITSHe3{"settingCutNsigmaITSHe3", -1.5f, "Value of the ITS Nsigma cutfor He3"};
  Configurable<float> settingCutPtMinTOFHad{"settingCutPtMinTOFHad", 0.4f, "Minimum pT to apply the TOF cut on hadrons"};
  Configurable<float> settingCutNsigmaTOF{"settingCutNsigmaTOF", 3.0f, "Value of the TOF Nsigma cut"};

  Configurable<LabeledArray<int>> settingEventSelections{"settingEventSelections", {nuclei::EvSelDefault[0], 8, 1, nuclei::eventSelectionLabels, nuclei::eventSelectionTitle}, "Event selections"};
  Configurable<int> settingNoMixedEvents{"settingNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> settingEnableBkgUS{"settingEnableBkgUS", false, "Enable US background"};
  Configurable<bool> settingEnableDCAfitter{"settingEnableDCAfitter", false, "Enable DCA fitter"};
  Configurable<bool> settingSaveUSandLS{"settingSaveUSandLS", true, "Save All Pairs"};
  Configurable<bool> settingIsMC{"settingIsMC", false, "Run MC"};

  Configurable<bool> settingFillMultiplicity{"settingFillMultiplicity", false, "Fill multiplicity table"};
  Configurable<bool> settingFillQa{"settingFillQa", false, "Fill QA table"};
  Configurable<bool> settingFillPrimariesAndMixedMc{"settingFillPrimariesAndMixedMc", false, "Fill primary MC tracks and mixed tracks (e.g. a primary track and one from Li4)"};

  // Zorro
  Configurable<bool> settingSkimmedProcessing{"settingSkimmedProcessing", false, "Skimmed dataset processing"};

  // svPool
  Configurable<bool> settingSkipAmbiTracks{"settingSkipAmbiTracks", false, "Skip ambiguous tracks"};
  Configurable<float> settingCustomVertexerTimeMargin{"settingCustomVertexerTimeMargin", 800, "Time margin for custom vertexer (ns)"};

  // CCDB options
  Configurable<double> settingDbz{"settingDbz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> settingCcdburl{"settingCcdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> settingGrpPath{"settingGrpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> settingGrpmagPath{"settingGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> settingLutPath{"settingLutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> settingGeoPath{"settingGeoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> settingPidPath{"settingPidPath", "", "Path to the PID response object"};

  Configurable<LabeledArray<double>> settingBetheBlochParams{"settingBetheBlochParams", {betheBlochDefault[0], 1, 6, {"He3"}, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  Configurable<bool> settingCompensatePIDinTracking{"settingCompensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};
  Configurable<int> settingMaterialCorrection{"settingMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Material correction type"};

  Preslice<TrackCandidates> mPerCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> mPerColMC = aod::track::collisionId;

  // binning for EM background
  ConfigurableAxis axisVertex{"axisVertex", {30, -10, 10}, "Binning for multiplicity"};
  ConfigurableAxis axisCentrality{"axisCentrality", {40, 0, 100}, "Binning for centrality"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType binningPolicy{{axisVertex, axisCentrality}, true};
  SliceCache cache;
  SameKindPair<CollisionsFull, TrackCandidates, BinningType> mPair{binningPolicy, settingNoMixedEvents, -1, &cache};

  std::array<float, 6> mBBparamsHe;
  o2::aod::ITSResponse mResponseITS;

  std::vector<int> mRecoCollisionIDs;
  std::vector<bool> mGoodCollisions;
  std::vector<SVCand> mTrackPairs;
  o2::vertexing::DCAFitterN<2> mFitter;
  svPoolCreator mSvPoolCreator{He3PDG, ProtonPDG};
  int mRunNumber;
  float mDbz;
  Service<o2::ccdb::BasicCCDBManager> mCcdb;
  Zorro mZorro;
  OutputObj<ZorroSummary> mZorroSummary{"zorroSummary"};

  HistogramRegistry mQaRegistry{
    "QA",
    {
      {"hEventSelections", "Event selections; Selection step; Counts", {HistType::kTH1D, {{nuclei::evSel::kNevSels + 1, -0.5f, static_cast<float>(nuclei::evSel::kNevSels) + 0.5f}}}},
      {"hEvents", "Number of events processed;Counts", {HistType::kTH1F, {{1, -0.5f, 2.5f}}}},
      {"hVtxZBefore", "Vertex distribution in Z before selections;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
      {"hVtxZ", "Vertex distribution in Z;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
      {"hCentralityFT0A", ";Centrality FT0A (%)", {HistType::kTH1F, {{100, 0, 100.0}}}},
      {"hCentralityFT0C", ";Centrality FT0C (%)", {HistType::kTH1F, {{100, 0, 100.0}}}},
      {"hNcontributor", "Number of primary vertex contributor", {HistType::kTH1F, {{2000, 0.0f, 2000.0f}}}},
      {"hTrackSel", "Accepted tracks", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
      {"hEmptyPool", "svPoolCreator did not find track pairs false/true", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"hhe3HadtInvMass", "; M(^{3}He + p) (GeV/#it{c}^{2})", {HistType::kTH1F, {{300, 3.74f, 4.34f}}}},
      {"hKstarRecVsKstarGen", "; #it{k}*_{gen} (GeV/#it{c}); #it{k}*_{rec} (GeV/#it{c})", {HistType::kTH2F, {{400, 0.f, 0.8f}, {400, 0.f, 0.8f}}}},

      {"He3/hDCAxyHe3", "^{3}He;DCA_{xy} (cm)", {HistType::kTH1F, {{200, -0.5f, 0.5f}}}},
      {"He3/hDCAzHe3", "^{3}He;DCA_{z} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
      {"He3/hNClsHe3ITS", "^{3}He;N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
      {"He3/hChi2NClHe3ITS", "^{3}He;Chi2_{ITS} Ncluster", {HistType::kTH1F, {{100, 0, 100.0f}}}},
      {"He3/hHe3Pt", "^{3}He; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"He3/h2dEdxHe3candidates", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
      {"He3/h2NsigmaHe3ITS", "NsigmaHe3 ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} ^{3}He", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {120, -3.0f, 3.0f}}}},
      {"He3/h2NsigmaHe3ITS_preselection", "NsigmaHe3 ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} ^{3}He", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {120, -3.0f, 3.0f}}}},
      {"He3/h2NsigmaHe3TPC", "NsigmaHe3 TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(^{3}He)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
      {"He3/h2NsigmaHe3TPC_preselection", "NsigmaHe3 TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(^{3}He)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},

      {"Had/hNClsHadITS", "had;N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
      {"Had/hChi2NClHadITS", "had;Chi2_{ITS} Ncluster", {HistType::kTH1F, {{100, 0, 100.0f}}}},
      {"Had/hHadronPt", "had; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},
      {"Had/h2NsigmaHadronITS", "NsigmaHadron ITS distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{ITS}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
      {"Had/h2NsigmaHadronITS_preselection", "NsigmaHadron ITS distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
      {"Had/h2NsigmaHadronTPC", "NsigmaHadron TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
      {"Had/h2NsigmaHadronTPC_preselection", "NsigmaHadron TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
      {"Had/h2NsigmaHadronTPC_mcBackground", "NsigmaHadron TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
      {"Had/h2NsigmaHadronTPC_mcSignal", "NsigmaHadron TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
      {"Had/h2NsigmaHadronTOF", "NsigmaHadron TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
      {"Had/h2NsigmaHadronTOF_preselection", "NsigmaHadron TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
      {"Had/h2NsigmaHadronTOF_mcBackground", "NsigmaHadron TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
      {"Had/h2NsigmaHadronTOF_mcSignal", "NsigmaHadron TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(had)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
    },
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

    mSvPoolCreator.setTimeMargin(settingCustomVertexerTimeMargin);
    if (settingSkipAmbiTracks) {
      mSvPoolCreator.setSkipAmbiTracks();
    }

    mQaRegistry.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(1, "All");
    for (int iSel = 1; iSel < nuclei::evSel::kNevSels + 1; iSel++) {
      mQaRegistry.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(iSel + 1, nuclei::eventSelectionLabels[iSel].c_str());
    }

    const int nBetheBlochParameters = 5;
    for (int i = 0; i < nBetheBlochParameters; i++) {
      mBBparamsHe[i] = settingBetheBlochParams->get("He3", Form("p%i", i));
    }
    mBBparamsHe[5] = settingBetheBlochParams->get("He3", "resolution");

    std::vector<std::string> selectionLabels = {"All", "Track selection", "PID"};
    for (int i = 0; i < Selections::kAll; i++) {
      mQaRegistry.get<TH1>(HIST("hTrackSel"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i].c_str());
    }

    std::vector<std::string> eventsLabels = {"All", "Selected", "Zorro He events"};
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
      mZorro.initCCDB(mCcdb.service, bc.runNumber(), bc.timestamp(), "fHe");
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

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!nuclei::eventSelection(collision, mQaRegistry, settingEventSelections, settingCutVertex)) {
      return false;
    }
    if (settingSkimmedProcessing) {
      bool zorroSelected = mZorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC());
      if (zorroSelected) {
        mQaRegistry.fill(HIST("hEvents"), 2);
      }
    }

    mQaRegistry.fill(HIST("hEvents"), 1);
    mQaRegistry.fill(HIST("hNcontributor"), collision.numContrib());
    mQaRegistry.fill(HIST("hCentralityFT0A"), collision.centFT0A());
    mQaRegistry.fill(HIST("hCentralityFT0C"), collision.centFT0C());
    return true;
  }

  template <typename Ttrack>
  bool selectTrack(const Ttrack& candidate, const int ispecies)
  {
    if (std::abs(candidate.eta()) > settingCutEta) {
      return false;
    }
    const int minTPCNClsFound = ispecies == Species::kHe3 ? static_cast<int>(settingCutNClsTPCHe3) : static_cast<int>(settingCutNClsTPCHe3);
    const int minTPCNClsCrossedRows = 70;
    const float crossedRowsToFindableRatio = 0.8f;
    const float minChi2NCl = ispecies == Species::kHe3 ? static_cast<int>(settingCutChi2tpcLowHe3) : static_cast<int>(settingCutChi2tpcLow);
    const float maxChi2NCl = 4.f;
    if (candidate.itsNCls() < settingCutNCls ||
        candidate.tpcNClsFound() < minTPCNClsFound ||
        candidate.tpcNClsCrossedRows() < minTPCNClsCrossedRows ||
        candidate.tpcNClsCrossedRows() < crossedRowsToFindableRatio * candidate.tpcNClsFindable() ||
        candidate.tpcChi2NCl() > maxChi2NCl ||
        candidate.tpcChi2NCl() < minChi2NCl ||
        candidate.itsChi2NCl() > settingCutChi2NClITS) {
      return false;
    }

    return true;
  }

  template <typename Ttrack>
  float correctPtHe3TrackedAsTriton(const Ttrack& candidate)
  {
    if (candidate.pt() * 2. < 2.5 && candidate.pidForTracking() == o2::track::PID::Triton)
      return candidate.pt() * 2. * (1. - kHePidTrkParams[0] - kHePidTrkParams[1] * candidate.pt() * 2.);

    return candidate.pt() * 2.;
  }

  float correctPtHe3TrackedAsTriton(const float pt, const uint32_t pidForTracking)
  {
    if (pt < 2.5 && pidForTracking == o2::track::PID::Triton)
      return pt * (1. - kHePidTrkParams[0] - kHePidTrkParams[1] * pt);

    return pt;
  }

  float computeNsigmaDCA(const float pt, const float dca, const int iSpecies, const char* dcaType = "xy")
  {

    std::array<float, 3> parameters;
    if (std::strcmp(dcaType, "xy") == 0) {
      parameters = kDCAxyResolutionParams[iSpecies];
    } else if (std::strcmp(dcaType, "z") == 0) {
      parameters = kDCAzResolutionParams[iSpecies];
    } else {
      LOG(error) << "Invalid dcaType. Accepted types are 'xy' 'z'";
      parameters = {0., 0., 0.};
    }
    const float sigma = parameters[0] *
                          std::exp(-std::abs(pt) * parameters[1]) +
                        parameters[2];
    return dca / sigma;
  }

  template <typename Ttrack>
  bool selectDcaNsigmaCut(const Ttrack& candidate, const int ispecies)
  {
    const float pt = ispecies == Species::kHe3 ? 2. * candidate.pt() : candidate.pt();
    const float nsigmaDcaXy = computeNsigmaDCA(pt, candidate.dcaXY(), ispecies, "xy");
    const float nsigmaDcaZ = computeNsigmaDCA(pt, candidate.dcaZ(), ispecies, "z");

    if (std::abs(nsigmaDcaXy) > settingCutNsigmaDcaXy ||
        std::abs(nsigmaDcaZ) > settingCutNsigmaDcaZ)
      return false;

    return true;
  }

  template <typename Ttrack>
  float computeTPCNSigmaHadron(const Ttrack& candidate)
  {
    float tpcNSigmaHad = CommonInite;
    if (settingHadPDGCode == PDG_t::kPiPlus) {
      tpcNSigmaHad = candidate.tpcNSigmaPi();
    } else if (settingHadPDGCode == PDG_t::kProton) {
      tpcNSigmaHad = candidate.tpcNSigmaPr();
    } else {
      LOG(info) << "invalid PDG code for TPC";
    }
    return tpcNSigmaHad;
  }

  template <typename Ttrack>
  float computeTOFNSigmaHadron(const Ttrack& candidate)
  {
    float tofNSigmaHad = CommonInite;
    if (settingHadPDGCode == PDG_t::kPiPlus) {
      tofNSigmaHad = candidate.tofNSigmaPi();
    } else if (settingHadPDGCode == PDG_t::kProton) {
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
    mQaRegistry.fill(HIST("Had/h2NsigmaHadronTPC_preselection"), candidate.tpcInnerParam(), tpcNSigmaHad);
    if (candidate.hasTOF() && candidate.pt() > settingCutPtMinTOFHad) {
      auto tofNSigmaHad = computeTOFNSigmaHadron(candidate);

      if (std::abs(tpcNSigmaHad) > settingCutNsigmaTPC) {
        return false;
      }
      mQaRegistry.fill(HIST("Had/h2NsigmaHadronTOF_preselection"), candidate.pt(), tofNSigmaHad);
      if (std::abs(tofNSigmaHad) > settingCutNsigmaTOF) {
        return false;
      }
      mQaRegistry.fill(HIST("Had/h2NsigmaHadronTPC"), candidate.pt(), tpcNSigmaHad);
      mQaRegistry.fill(HIST("Had/h2NsigmaHadronTOF"), candidate.pt(), tofNSigmaHad);
      return true;
    } else if (std::abs(tpcNSigmaHad) < settingCutNsigmaTPC) {
      mQaRegistry.fill(HIST("Had/h2NsigmaHadronTPC"), candidate.pt(), tpcNSigmaHad);
      return true;
    }
    return false;
  }

  template <typename Ttrack>
  float computeNSigmaHe3(const Ttrack& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParam = (heliumPID && settingCompensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    float expTPCSignal = o2::common::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), mBBparamsHe[0], mBBparamsHe[1], mBBparamsHe[2], mBBparamsHe[3], mBBparamsHe[4]);

    double resoTPC{expTPCSignal * mBBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename Ttrack>
  bool selectionPIDHe3(const Ttrack& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParam = (heliumPID && settingCompensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();

    if (correctedTPCinnerParam < settingCutRigidityMinHe3) {
      return false;
    }

    auto nSigmaHe3 = computeNSigmaHe3(candidate);
    mQaRegistry.fill(HIST("He3/h2NsigmaHe3TPC_preselection"), candidate.sign() * 2 * candidate.pt(), nSigmaHe3);
    if (std::abs(nSigmaHe3) > settingCutNsigmaTPC) {
      return false;
    }

    auto itsNsigmaHe3 = mResponseITS.nSigmaITS<o2::track::PID::Helium3>(candidate.itsClusterSizes(), 2 * candidate.p(), candidate.eta());
    mQaRegistry.fill(HIST("He3/h2NsigmaHe3ITS_preselection"), candidate.sign() * 2 * candidate.pt(), itsNsigmaHe3);
    if (itsNsigmaHe3 < settingCutNsigmaITSHe3) {
      return false;
    }

    mQaRegistry.fill(HIST("He3/h2dEdxHe3candidates"), candidate.sign() * correctedTPCinnerParam, candidate.tpcSignal());
    mQaRegistry.fill(HIST("He3/h2NsigmaHe3TPC"), candidate.sign() * 2 * candidate.pt(), nSigmaHe3);
    mQaRegistry.fill(HIST("He3/h2NsigmaHe3ITS"), candidate.sign() * 2 * candidate.pt(), itsNsigmaHe3);
    return true;
  }

  // ==================================================================================================================

  template <typename Tcollisions>
  std::array<float, 3> getCollisionVertex(const Tcollisions& collisions, int32_t collisionID)
  {
    auto collision = collisions.rawIteratorAt(collisionID);
    std::array<float, 3> collisionVertex = {collision.posX(), collision.posY(), collision.posZ()};
    return collisionVertex;
  }

  template <typename Ttrack, typename Tcollisions>
  int32_t getCollisionID(const Ttrack& trackHe3, const Ttrack& trackHad, const CollBracket& collBracket, const Tcollisions& collisions, bool isMixedEvent)
  {
    if (isMixedEvent) {
      return collBracket.getMin();
    }

    auto trackCovHe3 = getTrackParCov(trackHe3);
    auto trackCovHad = getTrackParCov(trackHad);
    int nCand = CommonInite;
    try {
      nCand = mFitter.process(trackCovHe3, trackCovHad);
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
    const float defaultTodistance = 0.0f;
    for (int collIdx = collBracket.getMin(); collIdx <= collBracket.getMax(); collIdx++) {
      std::array<float, 3> collisionVertex = getCollisionVertex(collisions, collIdx);
      const auto& pca = mFitter.getPCACandidate();
      float distance = defaultTodistance;
      for (int i = 0; i < 3; i++) {
        distance += (pca[i] - collisionVertex[i]) * (pca[i] - collisionVertex[i]);
      }
      if (distanceMin < 0 || distance < distanceMin) {
        distanceMin = distance;
        collIdxMin = collIdx;
      }
    }

    if (!mGoodCollisions[collIdxMin]) {
      return false;
    }
    return collIdxMin;
  }

  template <typename Ttrack, typename Tcollisions, typename Ttracks>
  bool fillCandidateInfo(const Ttrack& trackHe3, const Ttrack& trackHad, const CollBracket& collBracket, const Tcollisions& collisions, He3HadCandidate& he3Hadcand, const Ttracks& /*trackTable*/, bool isMixedEvent)
  {
    he3Hadcand.collisionID = getCollisionID(trackHe3, trackHad, collBracket, collisions, isMixedEvent);
    std::array<float, 3> collisionVertex = getCollisionVertex(collisions, he3Hadcand.collisionID);

    he3Hadcand.momHe3 = std::array{trackHe3.px(), trackHe3.py(), trackHe3.pz()};
    for (int i = 0; i < 3; i++)
      he3Hadcand.momHe3[i] = he3Hadcand.momHe3[i] * 2;
    he3Hadcand.momHad = std::array{trackHad.px(), trackHad.py(), trackHad.pz()};
    float invMass = CommonInite;
    if (settingHadPDGCode == PDG_t::kPiPlus) {
      invMass = RecoDecay::m(std::array{he3Hadcand.momHe3, he3Hadcand.momHad}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassPiPlus});
    } else if (settingHadPDGCode == PDG_t::kProton) {
      invMass = RecoDecay::m(std::array{he3Hadcand.momHe3, he3Hadcand.momHad}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassProton});
    } else {
      LOG(info) << "invalid PDG code for invMass";
    }

    if (settingCutInvMass > 0 && invMass > settingCutInvMass) {
      return false;
    }
    float pthe3Had = std::hypot(he3Hadcand.momHe3[0] + he3Hadcand.momHad[0], he3Hadcand.momHe3[1] + he3Hadcand.momHad[1]);
    if (pthe3Had < settingCutPtMinhe3Had) {
      return false;
    }

    he3Hadcand.signHe3 = trackHe3.sign();
    he3Hadcand.signHad = trackHad.sign();
    if (!settingEnableDCAfitter) {
      he3Hadcand.dcaxyHe3 = trackHe3.dcaXY();
      he3Hadcand.dcaxyHad = trackHad.dcaXY();
      he3Hadcand.dcazHe3 = trackHe3.dcaZ();
      he3Hadcand.dcazHad = trackHad.dcaZ();
    } else {
      auto trackCovHe3 = getTrackParCov(trackHe3);
      auto trackCovHad = getTrackParCov(trackHad);
      std::array<float, 2> dcaInfo;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collisionVertex[0], collisionVertex[1], collisionVertex[2]}, trackCovHe3, 2.f, mFitter.getMatCorrType(), &dcaInfo);
      he3Hadcand.dcaxyHe3 = dcaInfo[0];
      he3Hadcand.dcazHe3 = dcaInfo[1];
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collisionVertex[0], collisionVertex[1], collisionVertex[2]}, trackCovHad, 2.f, mFitter.getMatCorrType(), &dcaInfo);
      he3Hadcand.dcaxyHad = dcaInfo[0];
      he3Hadcand.dcazHad = dcaInfo[1];
      he3Hadcand.dcaPair = std::sqrt(std::abs(mFitter.getChi2AtPCACandidate()));
    }

    he3Hadcand.tpcSignalHe3 = trackHe3.tpcSignal();
    bool heliumPID = trackHe3.pidForTracking() == o2::track::PID::Helium3 || trackHe3.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParamHe3 = (heliumPID && settingCompensatePIDinTracking) ? trackHe3.tpcInnerParam() / 2.f : trackHe3.tpcInnerParam();
    he3Hadcand.momHe3TPC = correctedTPCinnerParamHe3;
    he3Hadcand.tpcSignalHad = trackHad.tpcSignal();
    he3Hadcand.momHadTPC = trackHad.tpcInnerParam();

    he3Hadcand.nTPCClustersHe3 = trackHe3.tpcNClsFound();
    he3Hadcand.nSigmaHe3 = computeNSigmaHe3(trackHe3);
    he3Hadcand.nSigmaTPCHadPi = trackHad.tpcNSigmaPi();
    he3Hadcand.nSigmaTPCHadKa = trackHad.tpcNSigmaKa();
    he3Hadcand.nSigmaTPCHadPr = trackHad.tpcNSigmaPr();
    he3Hadcand.nSigmaTOFHadPi = trackHad.tofNSigmaPi();
    he3Hadcand.nSigmaTOFHadKa = trackHad.tofNSigmaKa();
    he3Hadcand.nSigmaTOFHadPr = trackHad.tofNSigmaPr();

    he3Hadcand.chi2TPCHe3 = trackHe3.tpcChi2NCl();
    he3Hadcand.chi2TPCHad = trackHad.tpcChi2NCl();

    he3Hadcand.pidtrkHe3 = trackHe3.pidForTracking();
    he3Hadcand.pidtrkHad = trackHad.pidForTracking();

    he3Hadcand.itsClSizeHe3 = trackHe3.itsClusterSizes();
    he3Hadcand.itsClSizeHad = trackHad.itsClusterSizes();

    he3Hadcand.nclsITSHe3 = trackHe3.itsNCls();
    he3Hadcand.nclsITSHad = trackHad.itsNCls();
    he3Hadcand.chi2nclITSHe3 = trackHe3.itsChi2NCl();
    he3Hadcand.chi2nclITSHad = trackHad.itsChi2NCl();

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
      float correctedTPCinnerParamHe3 = (heliumPID && settingCompensatePIDinTracking) ? trackHe3.tpcInnerParam() / 2.f : trackHe3.tpcInnerParam();
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
  void fillCandidateInfoMC(const Mc& mctrackHe3, const Mc& mctrackHad, He3HadCandidate& he3Hadcand)
  {
    he3Hadcand.momHe3MC = mctrackHe3.pt() * (mctrackHe3.pdgCode() > 0 ? 1 : -1);
    he3Hadcand.etaHe3MC = mctrackHe3.eta();
    he3Hadcand.phiHe3MC = mctrackHe3.phi();
    he3Hadcand.momHadMC = mctrackHad.pt() * (mctrackHad.pdgCode() > 0 ? 1 : -1);
    he3Hadcand.etaHadMC = mctrackHad.eta();
    he3Hadcand.phiHadMC = mctrackHad.phi();
  }

  template <typename Mc>
  void fillMotherInfoMC(const Mc& mctrackHe3, const Mc& mctrackHad, const Mc& mctrackMother, He3HadCandidate& he3Hadcand)
  {
    he3Hadcand.l4PtMC = mctrackMother.pt() * (mctrackMother.pdgCode() > 0 ? 1 : -1);
    const double eLit = mctrackHe3.e() + mctrackHad.e();
    he3Hadcand.l4MassMC = std::sqrt(eLit * eLit - mctrackMother.p() * mctrackMother.p());
  }

  template <typename Ttrack>
  void pairTracksSameEvent(const Ttrack& tracks)
  {
    for (const auto& track0 : tracks) {

      mQaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);

      if (!selectTrack(track0, Species::kHe3)) {
        continue;
      }
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);

      if (!selectionPIDHe3(track0)) {
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

        if (!selectTrack(track1, Species::kHad) || !selectionPIDHadron(track1)) {
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
  void pairTracksEventMixing(T& he3Cands, T& hadronCands)
  {
    for (const auto& he3Cand : he3Cands) {
      if (!selectTrack(he3Cand, Species::kHe3) || !selectionPIDHe3(he3Cand)) {
        continue;
      }
      for (const auto& hadronCand : hadronCands) {
        if (!selectTrack(hadronCand, Species::kHad) || !selectionPIDHadron(hadronCand)) {
          continue;
        }

        SVCand trackPair;
        trackPair.tr0Idx = he3Cand.globalIndex();
        trackPair.tr1Idx = hadronCand.globalIndex();
        const int collIdx = he3Cand.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        trackPair.collBracket = collBracket;
        mTrackPairs.push_back(trackPair);
      }
    }
  }

  template <typename Tcoll>
  void fillTable(const He3HadCandidate& he3Hadcand, const Tcoll& collision, bool isMC = false)
  {
    outputDataTable(
      he3Hadcand.recoPtHe3(),
      he3Hadcand.recoEtaHe3(),
      he3Hadcand.recoPhiHe3(),
      he3Hadcand.recoPtHad(),
      he3Hadcand.recoEtaHad(),
      he3Hadcand.recoPhiHad(),
      he3Hadcand.dcaxyHe3,
      he3Hadcand.dcazHe3,
      he3Hadcand.dcaxyHad,
      he3Hadcand.dcazHad,
      he3Hadcand.dcaPair,
      he3Hadcand.tpcSignalHe3,
      he3Hadcand.momHe3TPC,
      he3Hadcand.tpcSignalHad,
      he3Hadcand.momHadTPC,
      he3Hadcand.nTPCClustersHe3,
      he3Hadcand.nSigmaHe3,
      he3Hadcand.nSigmaTPCHadPi,
      he3Hadcand.nSigmaTPCHadKa,
      he3Hadcand.nSigmaTPCHadPr,
      he3Hadcand.nSigmaTOFHadPi,
      he3Hadcand.nSigmaTOFHadKa,
      he3Hadcand.nSigmaTOFHadPr,
      he3Hadcand.chi2TPCHe3,
      he3Hadcand.chi2TPCHad,
      he3Hadcand.massTOFHe3,
      he3Hadcand.massTOFHad,
      he3Hadcand.pidtrkHe3,
      he3Hadcand.pidtrkHad,
      he3Hadcand.itsClSizeHe3,
      he3Hadcand.itsClSizeHad,
      he3Hadcand.sharedClustersHe3,
      he3Hadcand.sharedClustersHad);
    if (isMC) {
      outputMcTable(
        he3Hadcand.momHe3MC,
        he3Hadcand.etaHe3MC,
        he3Hadcand.phiHe3MC,
        he3Hadcand.momHadMC,
        he3Hadcand.etaHadMC,
        he3Hadcand.phiHadMC,
        he3Hadcand.l4PtMC,
        he3Hadcand.l4MassMC,
        he3Hadcand.flags);
    }
    if (settingFillMultiplicity) {
      outputMultiplicityTable(
        collision.globalIndex(),
        collision.posZ(),
        collision.numContrib(),
        collision.centFT0C(),
        collision.multFT0C());
    }
    if (settingFillQa) {
      outputQaTable(
        he3Hadcand.trackIDHe3,
        he3Hadcand.trackIDHad);
    }
  }

  void fillHistograms(const He3HadCandidate& he3Hadcand, bool isMc = false)
  {
    mQaRegistry.fill(HIST("He3/hHe3Pt"), he3Hadcand.recoPtHe3());
    mQaRegistry.fill(HIST("Had/hHadronPt"), he3Hadcand.recoPtHad());
    mQaRegistry.fill(HIST("hhe3HadtInvMass"), he3Hadcand.invMass);
    mQaRegistry.fill(HIST("He3/hDCAxyHe3"), he3Hadcand.dcaxyHe3);
    mQaRegistry.fill(HIST("He3/hDCAzHe3"), he3Hadcand.dcazHe3);
    mQaRegistry.fill(HIST("He3/hNClsHe3ITS"), he3Hadcand.nclsITSHe3);
    mQaRegistry.fill(HIST("Had/hNClsHadITS"), he3Hadcand.nclsITSHad);
    mQaRegistry.fill(HIST("He3/hChi2NClHe3ITS"), he3Hadcand.chi2nclITSHe3);
    mQaRegistry.fill(HIST("Had/hChi2NClHadITS"), he3Hadcand.chi2nclITSHad);

    if (isMc) {
      const float correctedPtHe3 = correctPtHe3TrackedAsTriton(std::abs(he3Hadcand.recoPtHe3()), he3Hadcand.pidtrkHe3);
      const float kstarGen = getkstar(std::abs(he3Hadcand.momHe3MC), he3Hadcand.etaHe3MC, he3Hadcand.phiHe3MC, o2::constants::physics::MassHelium3, 1.,
                                      std::abs(he3Hadcand.momHadMC), he3Hadcand.etaHadMC, he3Hadcand.phiHadMC, settingHadPDGCode == PDG_t::kPiPlus ? o2::constants::physics::MassPiPlus : o2::constants::physics::MassProton, 1.);
      const float kstarRec = getkstar(std::abs(correctedPtHe3), he3Hadcand.recoEtaHe3(), he3Hadcand.recoPhiHe3(), o2::constants::physics::MassHelium3, 1.,
                                      std::abs(he3Hadcand.recoPtHad()), he3Hadcand.recoEtaHad(), he3Hadcand.recoPhiHad(), settingHadPDGCode == PDG_t::kPiPlus ? o2::constants::physics::MassPiPlus : o2::constants::physics::MassProton, 1.);
      mQaRegistry.fill(HIST("hKstarRecVsKstarGen"), kstarGen, kstarRec);
    }
  }

  // ==================================================================================================================

  template <typename Tcollisions, typename Ttracks>
  void fillPairs(const Tcollisions& collisions, const Ttracks& tracks, const bool isMixedEvent)
  {
    for (const auto& trackPair : mTrackPairs) {

      auto heTrack = tracks.rawIteratorAt(trackPair.tr0Idx);
      auto hadTrack = tracks.rawIteratorAt(trackPair.tr1Idx);
      auto collBracket = trackPair.collBracket;

      He3HadCandidate he3Hadcand;
      if (!fillCandidateInfo(heTrack, hadTrack, collBracket, collisions, he3Hadcand, tracks, isMixedEvent)) {
        continue;
      }
      fillHistograms(he3Hadcand);
      auto collision = collisions.rawIteratorAt(he3Hadcand.collisionID);
      fillTable(he3Hadcand, collision, /*isMC*/ false);
    }
  }

  template <typename TmcParticle>
  void setMcParticleFlag(const TmcParticle& mcParticle, std::vector<unsigned int>& mothers, uint8_t& flag)
  {
    if (mcParticle.isPhysicalPrimary()) {

      flag |= ParticleFlags::kPhysicalPrimary;
      if (!mcParticle.has_mothers()) {
        return;
      }

      for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
        mothers.push_back(mother.globalIndex());
        if (std::abs(mother.pdgCode()) == Li4PDG) {
          flag |= ParticleFlags::kFromLi4;
        } else if (std::abs(mother.pdgCode()) == H3LPDG) {
          flag |= ParticleFlags::kFromHypertriton;
        } else {
          flag |= ParticleFlags::kFromOtherDecays;
        }
      }

    } else {

      if (!mcParticle.has_mothers()) {
        flag |= ParticleFlags::kFromMaterial;
        return;
      }

      for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
        mothers.push_back(mother.globalIndex());
        if (std::abs(mother.pdgCode()) == Li4PDG) {
          flag |= ParticleFlags::kFromLi4;
        } else if (std::abs(mother.pdgCode()) == H3LPDG) {
          flag |= ParticleFlags::kFromHypertriton;
        } else {
          flag |= ParticleFlags::kFromOtherDecays;
        }
      }
    }
  }

  void searchForCommonMotherTrack(const std::vector<unsigned int>& motherHe3Idxs, const std::vector<unsigned int>& motherHadIdxs, const aod::McParticles& mcParticles, McIter& motherParticle, He3HadCandidate& he3Hadcand, bool& isMixedPair, const int motherPdgCode)
  {
    std::unordered_set<unsigned int> motherHe3SetIdxs(motherHe3Idxs.begin(), motherHe3Idxs.end());
    for (const auto& motherHadIdx : motherHadIdxs) {
      if (!motherHe3SetIdxs.contains(motherHadIdx)) {
        continue;
      }

      motherParticle = mcParticles.rawIteratorAt(motherHadIdx);
      if (std::abs(motherParticle.pdgCode()) != motherPdgCode || std::abs(motherParticle.y()) > 1) {
        continue;
      }
      isMixedPair = false;
      break;
    }
    if (!isMixedPair) {
      he3Hadcand.flags |= Flags::kBothFromLi4;
    }
  }

  template <typename Tcollisions, typename TmcParticles>
  void fillMcParticles(const Tcollisions& collisions, const TmcParticles& mcParticles, std::vector<unsigned int>& filledMothers)
  {
    for (const auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.pdgCode()) != Li4PDG || std::abs(mcParticle.y()) > 1 || mcParticle.isPhysicalPrimary() == false) {
        continue;
      }

      if (std::find(filledMothers.begin(), filledMothers.end(), mcParticle.globalIndex()) != filledMothers.end()) {
        continue;
      }

      auto kDaughters = mcParticle.template daughters_as<aod::McParticles>();
      bool daughtHe3(false), daughtHad(false);
      McIter mcHe3, mcHad;
      for (const auto& kCurrentDaughter : kDaughters) {
        if (std::abs(kCurrentDaughter.pdgCode()) == He3PDG) {
          daughtHe3 = true;
          mcHe3 = kCurrentDaughter;
        } else if (std::abs(kCurrentDaughter.pdgCode()) == ProtonPDG) {
          daughtHad = true;
          mcHad = kCurrentDaughter;
        }
      }
      if (daughtHe3 && daughtHad) {
        He3HadCandidate he3Hadcand;
        fillCandidateInfoMC(mcHe3, mcHad, he3Hadcand);
        fillMotherInfoMC(mcHe3, mcHad, mcParticle, he3Hadcand);
        auto collision = collisions.rawIteratorAt(he3Hadcand.collisionID);
        fillTable(he3Hadcand, collision, /*isMC*/ true);
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
  PROCESS_SWITCH(he3HadronFemto, processSameEvent, "Process Same event", false);

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
  PROCESS_SWITCH(he3HadronFemto, processMixedEvent, "Process Mixed event", false);

  void processMC(const CollisionsFullMC& collisions, const aod::BCsWithTimestamps& bcs, const TrackCandidatesMC& tracks, const aod::McParticles& mcParticles)
  {
    std::vector<unsigned int> filledMothers;

    mGoodCollisions.clear();
    mGoodCollisions.resize(collisions.size(), false);

    for (const auto& collision : collisions) {

      mTrackPairs.clear();

      if (!selectCollision</*isMC*/ true>(collision, bcs)) {
        continue;
      }

      const uint64_t collIdx = collision.globalIndex();
      mGoodCollisions[collIdx] = true;
      auto trackTableThisCollision = tracks.sliceBy(mPerColMC, collIdx);
      trackTableThisCollision.bindExternalIndices(&tracks);

      pairTracksSameEvent(trackTableThisCollision);

      for (const auto& trackPair : mTrackPairs) {

        auto heTrack = tracks.rawIteratorAt(trackPair.tr0Idx);
        auto prTrack = tracks.rawIteratorAt(trackPair.tr1Idx);
        auto collBracket = trackPair.collBracket;

        if (!heTrack.has_mcParticle() || !prTrack.has_mcParticle()) {
          continue;
        }

        auto mctrackHe3 = heTrack.mcParticle();
        auto mctrackHad = prTrack.mcParticle();

        if (std::abs(mctrackHe3.pdgCode()) != He3PDG || std::abs(mctrackHad.pdgCode()) != settingHadPDGCode) {
          continue;
        }

        He3HadCandidate he3Hadcand;
        McIter motherParticle;
        std::vector<unsigned int> motherHe3Idxs, motherHadIdxs;
        setMcParticleFlag(mctrackHe3, motherHe3Idxs, he3Hadcand.flagsHe3);
        setMcParticleFlag(mctrackHad, motherHadIdxs, he3Hadcand.flagsHad);

        bool isMixedPair = true;

        if ((he3Hadcand.flagsHe3 == ParticleFlags::kPhysicalPrimary && he3Hadcand.flagsHad == ParticleFlags::kPhysicalPrimary)) {
          he3Hadcand.flags |= Flags::kBothPrimaries;
          isMixedPair = false;

        } else if ((he3Hadcand.flagsHe3 & ParticleFlags::kFromLi4) && (he3Hadcand.flagsHad & ParticleFlags::kFromLi4)) {

          searchForCommonMotherTrack(motherHe3Idxs, motherHadIdxs, mcParticles, motherParticle, he3Hadcand, isMixedPair, Li4PDG);
          if (!isMixedPair) {
            he3Hadcand.flags |= Flags::kBothFromLi4;
          }

        } else if ((he3Hadcand.flagsHe3 & ParticleFlags::kFromHypertriton) && (he3Hadcand.flagsHad & ParticleFlags::kFromHypertriton)) {

          searchForCommonMotherTrack(motherHe3Idxs, motherHadIdxs, mcParticles, motherParticle, he3Hadcand, isMixedPair, H3LPDG);
          if (!isMixedPair) {
            he3Hadcand.flags |= Flags::kBothFromHypertriton;
          }
        }

        if (isMixedPair) {
          he3Hadcand.flags |= Flags::kMixedPair;
        }

        if (!settingFillPrimariesAndMixedMc && ((he3Hadcand.flags == Flags::kMixedPair) || he3Hadcand.flags == Flags::kBothPrimaries)) {
          continue;
        }

        if (!fillCandidateInfo(heTrack, prTrack, collBracket, collisions, he3Hadcand, tracks, /*mix*/ false)) {
          continue;
        }
        fillCandidateInfoMC(mctrackHe3, mctrackHad, he3Hadcand);

        if ((he3Hadcand.flags == Flags::kBothFromLi4) || (he3Hadcand.flags == Flags::kBothFromHypertriton)) {
          fillMotherInfoMC(mctrackHe3, mctrackHad, motherParticle, he3Hadcand);
          filledMothers.push_back(motherParticle.globalIndex());
        }

        fillHistograms(he3Hadcand, /*isMc*/ true);
        auto collision = collisions.rawIteratorAt(he3Hadcand.collisionID);
        fillTable(he3Hadcand, collision, /*isMC*/ true);
      }
    }

    fillMcParticles(collisions, mcParticles, filledMothers);
  }
  PROCESS_SWITCH(he3HadronFemto, processMC, "Process MC", false);

  void processSameEventPools(const CollisionsFull& collisions, const TrackCandidates& tracks, const aod::AmbiguousTracks& ambiguousTracks, const aod::BCsWithTimestamps& bcs)
  {
    mGoodCollisions.clear();
    mGoodCollisions.resize(collisions.size(), false);

    for (const auto& collision : collisions) {
      if (selectCollision</*isMC*/ false>(collision, bcs)) {
        mGoodCollisions[collision.globalIndex()] = true;
      }
    }

    mSvPoolCreator.clearPools();
    mSvPoolCreator.fillBC2Coll(collisions, bcs);

    for (const auto& track : tracks) {

      mQaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);
      if (!selectTrack(track, Species::kHad)) // specific he3 cuts skipped here, might need to refactor this
        continue;
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);

      bool selHad = selectionPIDHadron(track);
      bool selHe = selectionPIDHe3(track);
      if ((!selHad && !selHe) || (selHad && selHe)) {
        continue;
      }
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kPID);

      int pdgHypo = selHe ? He3PDG : ProtonPDG;

      mSvPoolCreator.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }

    mTrackPairs = mSvPoolCreator.getSVCandPool(collisions, true);
    if (mTrackPairs.size() == 0) {
      mQaRegistry.fill(HIST("hEmptyPool"), 1);
      return;
    }
    mQaRegistry.fill(HIST("hEmptyPool"), 0);

    fillPairs(collisions, tracks, /*isMixedEvent*/ false);
  }
  PROCESS_SWITCH(he3HadronFemto, processSameEventPools, "Process Same event pools", false);

  void processMcPools(const CollisionsFullMC& collisions, const TrackCandidatesMC& tracks, const aod::AmbiguousTracks& ambiguousTracks, const aod::BCsWithTimestamps& bcs, const aod::McParticles& mcParticles, const aod::McTrackLabels& mcTrackLabels)
  {
    std::vector<unsigned int> filledMothers;

    mGoodCollisions.clear();
    mGoodCollisions.resize(collisions.size(), false);

    for (const auto& collision : collisions) {
      if (selectCollision</*isMC*/ true>(collision, bcs)) {
        mGoodCollisions[collision.globalIndex()] = true;
      }
    }

    mSvPoolCreator.clearPools();
    mSvPoolCreator.fillBC2Coll(collisions, bcs);

    for (const auto& track : tracks) {

      mQaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);
      if (!selectTrack(track, Species::kHad)) // specific he3 cuts skipped here, might need to refactor this
        continue;
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);

      bool selHad = selectionPIDHadron(track);
      bool selHe = selectionPIDHe3(track);
      if ((!selHad && !selHe) || (selHad && selHe))
        continue;
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kPID);

      int pdgHypo = selHe ? He3PDG : ProtonPDG;

      mSvPoolCreator.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }

    auto& svPool = mSvPoolCreator.getSVCandPool(collisions, true);
    if (svPool.size() == 0) {
      mQaRegistry.fill(HIST("hEmptyPool"), 1);
      return;
    }
    mQaRegistry.fill(HIST("hEmptyPool"), 0);

    for (const auto& svCand : svPool) {
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

      if (std::abs(mctrackHe3.pdgCode()) != He3PDG || std::abs(mctrackHad.pdgCode()) != ProtonPDG || !mctrackHe3.has_mothers() || !mctrackHad.has_mothers()) {
        continue;
      }

      for (const auto& mothertrackHe : mctrackHe3.mothers_as<aod::McParticles>()) {
        for (const auto& mothertrackHad : mctrackHad.mothers_as<aod::McParticles>()) {

          if (mothertrackHe.globalIndex() != mothertrackHad.globalIndex() || std::abs(mothertrackHe.pdgCode()) != Li4PDG || std::abs(mothertrackHe.y()) > 1) {
            continue;
          }

          He3HadCandidate he3Hadcand;
          if (!fillCandidateInfo(heTrack, prTrack, collBracket, collisions, he3Hadcand, tracks, /*mix*/ false)) {
            continue;
          }
          fillCandidateInfoMC(mctrackHe3, mctrackHad, he3Hadcand);
          fillMotherInfoMC(mctrackHe3, mctrackHad, mothertrackHe, he3Hadcand);
          fillHistograms(he3Hadcand);
          auto collision = collisions.rawIteratorAt(he3Hadcand.collisionID);
          fillTable(he3Hadcand, collision, /*isMC*/ true);
          filledMothers.push_back(mothertrackHe.globalIndex());
        }
      }
    }

    fillMcParticles(collisions, mcParticles, filledMothers);
  }
  PROCESS_SWITCH(he3HadronFemto, processMcPools, "Process MC pools", false);

  void processPurity(const CollisionsFull::iterator& collision, const TrackCandidates& tracks, const aod::BCsWithTimestamps& bcs)
  {
    if (!selectCollision</*isMC*/ false>(collision, bcs))
      return;

    for (const auto& track : tracks) {

      if (!selectTrack(track, Species::kHad))
        continue;

      const float itsNSigmaHad = settingHadPDGCode == PDG_t::kProton ? mResponseITS.nSigmaITS<o2::track::PID::Proton>(track.itsClusterSizes(), track.p(), track.eta()) : mResponseITS.nSigmaITS<o2::track::PID::Pion>(track.itsClusterSizes(), track.p(), track.eta());
      mQaRegistry.fill(HIST("Had/h2NsigmaHadronITS_preselection"), track.sign() * track.pt(), itsNSigmaHad);

      if (selectDcaNsigmaCut(track, Species::kHad) && (itsNSigmaHad > settingCutNsigmaITSHad)) {

        mQaRegistry.fill(HIST("Had/hHadronPt"), track.sign() * track.pt());
        mQaRegistry.fill(HIST("Had/h2NsigmaHadronITS"), track.sign() * track.pt(), itsNSigmaHad);

        const float tpcNSigmaHad = computeTPCNSigmaHadron(track);
        mQaRegistry.fill(HIST("Had/h2NsigmaHadronTPC_preselection"), track.sign() * track.pt(), tpcNSigmaHad);

        if (track.hasTOF()) {
          const float tofNSigmaHad = computeTOFNSigmaHadron(track);
          mQaRegistry.fill(HIST("Had/h2NsigmaHadronTOF_preselection"), track.sign() * track.pt(), tofNSigmaHad);
        }
      }

      const float ptHe3Corrected = correctPtHe3TrackedAsTriton(track);
      const float itsNSigmaHe3 = mResponseITS.nSigmaITS<o2::track::PID::Helium3>(track.itsClusterSizes(), 2 * track.p(), track.eta());
      mQaRegistry.fill(HIST("He3/h2NsigmaHe3ITS_preselection"), ptHe3Corrected, itsNSigmaHe3);

      if (!selectTrack(track, Species::kHe3) || !selectDcaNsigmaCut(track, Species::kHe3) || (itsNSigmaHe3 < settingCutNsigmaITSHe3))
        continue;

      mQaRegistry.fill(HIST("He3/hHe3Pt"), track.sign() * ptHe3Corrected);
      mQaRegistry.fill(HIST("He3/hDCAxyHe3"), track.dcaXY());
      mQaRegistry.fill(HIST("He3/hDCAzHe3"), track.dcaZ());
      mQaRegistry.fill(HIST("He3/h2NsigmaHe3ITS"), track.sign() * ptHe3Corrected, itsNSigmaHe3);

      bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
      float correctedTPCinnerParam = (heliumPID && settingCompensatePIDinTracking) ? track.tpcInnerParam() / 2.f : track.tpcInnerParam();
      if (correctedTPCinnerParam < settingCutRigidityMinHe3) {
        continue;
      }

      const float tpcNSigmaHe3 = computeNSigmaHe3(track);
      mQaRegistry.fill(HIST("He3/h2NsigmaHe3TPC_preselection"), track.sign() * ptHe3Corrected, tpcNSigmaHe3);
    }
  }
  PROCESS_SWITCH(he3HadronFemto, processPurity, "Process for purity studies", false);

  void processPurityMc(const CollisionsFullMC::iterator& collision, const TrackCandidatesMC& tracks, const aod::BCsWithTimestamps& bcs, const aod::McParticles& /*mcParticles*/, const aod::McTrackLabels& /*mcTrackLabels*/)
  {
    if (!selectCollision</*isMC*/ false>(collision, bcs))
      return;

    for (const auto& track : tracks) {

      if (!track.has_mcParticle())
        continue;
      const auto& particle = track.mcParticle_as<aod::McParticles>();

      if (!selectTrack(track, Species::kHad))
        continue;

      if (selectDcaNsigmaCut(track, Species::kHad)) {
        mQaRegistry.fill(HIST("Had/hHadronPt"), track.pt());

        const float tpcNSigmaHad = computeTPCNSigmaHadron(track);
        mQaRegistry.fill(HIST("Had/h2NsigmaHadronTPC_preselection"), track.sign() * track.pt(), tpcNSigmaHad);

        if (std::abs(particle.pdgCode()) != settingHadPDGCode) {
          mQaRegistry.fill(HIST("Had/h2NsigmaHadronTPC_mcBackground"), track.sign() * track.pt(), tpcNSigmaHad);
        } else {
          mQaRegistry.fill(HIST("Had/h2NsigmaHadronTPC_mcSignal"), track.sign() * track.pt(), tpcNSigmaHad);
        }

        if (track.hasTOF()) {
          const float tofNSigmaHad = computeTOFNSigmaHadron(track);
          mQaRegistry.fill(HIST("Had/h2NsigmaHadronTOF_preselection"), track.sign() * track.pt(), tofNSigmaHad);

          if (std::abs(particle.pdgCode()) != settingHadPDGCode) {
            mQaRegistry.fill(HIST("Had/h2NsigmaHadronTOF_mcBackground"), track.sign() * track.pt(), tofNSigmaHad);
          } else {
            mQaRegistry.fill(HIST("Had/h2NsigmaHadronTOF_mcSignal"), track.sign() * track.pt(), tofNSigmaHad);
          }
        }
      }

      if (!selectTrack(track, Species::kHe3) || !selectDcaNsigmaCut(track, Species::kHe3))
        continue;

      mQaRegistry.fill(HIST("He3/hHe3Pt"), track.pt() * 2.f);
      mQaRegistry.fill(HIST("He3/hDCAxyHe3"), track.dcaXY());
      mQaRegistry.fill(HIST("He3/hDCAzHe3"), track.dcaZ());

      bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
      float correctedTPCinnerParam = (heliumPID && settingCompensatePIDinTracking) ? track.tpcInnerParam() / 2.f : track.tpcInnerParam();
      if (correctedTPCinnerParam < settingCutRigidityMinHe3) {
        continue;
      }

      const float nSigmaHe3 = computeNSigmaHe3(track);
      mQaRegistry.fill(HIST("He3/h2NsigmaHe3TPC_preselection"), track.sign() * 2 * track.pt(), nSigmaHe3);
    }
  }
  PROCESS_SWITCH(he3HadronFemto, processPurityMc, "Process for purity studies mc", false);
};

WorkflowSpec defineDataProcessing(const ConfigContext& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<he3HadronFemto>(cfgc)};
}
