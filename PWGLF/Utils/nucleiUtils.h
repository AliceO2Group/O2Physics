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

#ifndef PWGLF_UTILS_NUCLEIUTILS_H_
#define PWGLF_UTILS_NUCLEIUTILS_H_

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include "TMCProcess.h"

#include <algorithm>
#include <string>
#include <vector>

struct NucleusCandidate {
  int globalIndex;
  int collTrackIndex;
  float pt;
  float eta;
  float phi;
  float tpcInnerParam;
  float beta;
  float zVertex;
  int nContrib;
  float DCAxy;
  float DCAz;
  float TPCsignal;
  float ITSchi2;
  float TPCchi2;
  float TOFchi2;
  std::array<float, 5> nSigmaTPC;
  std::array<float, 5> tofMasses;
  bool fillTree;
  bool fillDCAHist;
  bool correctPV;
  bool isSecondary;
  bool fromWeakDecay;
  uint16_t flags;
  uint8_t TPCfindableCls;
  uint8_t TPCcrossedRows;
  uint8_t ITSclsMap;
  uint8_t TPCnCls;
  uint8_t TPCnClsShared;
  uint8_t ITSnCls;
  uint32_t clusterSizesITS;
};

struct NucleusCandidateFlow {
  float centFV0A;
  float centFT0M;
  float centFT0A;
  float centFT0C;
  float psiFT0A;
  float psiFT0C;
  float psiTPC;
  float psiTPCl;
  float psiTPCr;
  float qFT0A;
  float qFT0C;
  float qTPC;
  float qTPCl;
  float qTPCr;
};

namespace nuclei
{

struct SlimCandidate {
  float pt = -999.f;
  float eta = -999.f;
  float phi = -999.f;
  float tpcInnerParam = -999.f;
  uint32_t clusterSizesITS = 0.f;
  float TPCsignal = -999.f;
  float beta = -999.f;
  float DCAxy = -999.f;
  float DCAz = -999.f;
  uint16_t flags = 0;
  int pdgCode = 0;
  int motherPdgCode = 0;
  float ptGenerated = -999.f;
  float etaGenerated = -999.f;
  float yGenerated = -999.f;
  float phiGenerated = -999.f;
  float centrality = -1.f;
  uint64_t mcProcess = TMCProcess::kPNoProcess;
};

enum Species {
  kPr = 0,
  kDe = 1,
  kTr = 2,
  kHe = 3,
  kAl = 4,
  kNspecies = 5
};

enum Flags {
  kProton = BIT(0),
  kDeuteron = BIT(1),
  kTriton = BIT(2),
  kHe3 = BIT(3),
  kHe4 = BIT(4),
  kHasTOF = BIT(5),
  kHasTRD = BIT(6),
  kIsAmbiguous = BIT(7), /// just a placeholder now
  kITSrof = BIT(8),
  kIsPhysicalPrimary = BIT(9), /// MC flags starting from the second half of the short
  kIsSecondaryFromMaterial = BIT(10),
  kIsSecondaryFromWeakDecay = BIT(11) /// the last 4 bits are reserved for the PID in tracking
};

constexpr int getSpeciesFromPdg(int pdg)
{
  switch (std::abs(pdg)) {
    case PDG_t::kProton:
      return Species::kPr;
    case o2::constants::physics::Pdg::kDeuteron:
      return Species::kDe;
    case o2::constants::physics::Pdg::kTriton:
      return Species::kTr;
    case o2::constants::physics::Pdg::kHelium3:
      return Species::kHe;
    case o2::constants::physics::Pdg::kAlpha:
      return Species::kAl;
    default:
      return -1;
  }
}

bool checkSpeciesValidity(const int species)
{
  if (species < 0 || species > Species::kNspecies) {
    return false;
  }
  return true;
}

constexpr int speciesToProcessDefault[Species::kNspecies][1]{
  {0},
  {0},
  {0},
  {0},
  {0}};
constexpr int useCentralTpcCalibrationDefault[Species::kNspecies][1]{
  {1},
  {1},
  {1},
  {1},
  {1}};
constexpr int useTrackTuner[Species::kNspecies][1]{
  {1},
  {1},
  {1},
  {0},
  {0}};
constexpr double bbMomScalingDefault[Species::kNspecies][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr double betheBlochDefault[Species::kNspecies][6]{
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-239.99, 1.155, 1.099, 1.137, 1.006, 0.09},
  {-321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09},
  {-586.66, 1.859, 4.435, 0.282, 3.201, 0.09}};
constexpr double nSigmaTPCdefault[Species::kNspecies][2]{
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.}};
constexpr double nSigmaTOFdefault[Species::kNspecies][2]{
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.}};
constexpr double DCAcutDefault[Species::kNspecies][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr int TreeConfigDefault[Species::kNspecies][2]{
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}};
constexpr int FlowHistDefault[Species::kNspecies][1]{
  {0},
  {0},
  {0},
  {0},
  {0}};
constexpr int DCAHistDefault[Species::kNspecies][2]{
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}};
constexpr double DownscalingDefault[Species::kNspecies][1]{
  {1.},
  {1.},
  {1.},
  {1.},
  {1.}};
// constexpr bool storeTreesDefault[Species::kNspecies]{false, false, false, false, false};
constexpr float charges[Species::kNspecies]{1.f, 1.f, 1.f, 2.f, 2.f};
constexpr float masses[Species::kNspecies]{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron, o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3, o2::constants::physics::MassAlpha};
constexpr int pdgCodes[Species::kNspecies]{PDG_t::kProton, o2::constants::physics::Pdg::kDeuteron, o2::constants::physics::Pdg::kTriton, o2::constants::physics::Pdg::kHelium3, o2::constants::physics::Pdg::kAlpha};
static constexpr std::string_view cNames[] = {"proton", "deuteron", "triton", "He3", "alpha"};
static const std::vector<std::string> names{"proton", "deuteron", "triton", "He3", "alpha"};
static const std::vector<std::string> matter{"M", "A"};
static const std::vector<std::string> pidName{"TPC", "TOF"};
static const std::vector<std::string> treeConfigNames{"Filter trees", "Use TOF selection"};
static const std::vector<std::string> flowConfigNames{"Save flow hists"};
static const std::vector<std::string> DCAConfigNames{"Save DCA hist", "Matter/Antimatter"};
static const std::vector<std::string> nSigmaConfigName{"nsigma_min", "nsigma_max"};
static const std::vector<std::string> nDCAConfigName{"max DCAxy", "max DCAz"};
static const std::vector<std::string> DownscalingConfigName{"Fraction of kept candidates"};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> chargeLabelNames{"Positive", "Negative"};

float pidCutTPC[Species::kNspecies][2]; //[species][lower/upper limit]

o2::base::MatLayerCylSet* lut = nullptr;

std::vector<NucleusCandidate> candidates;
std::vector<NucleusCandidateFlow> candidates_flow;

enum centDetectors {
  kFV0A = 0,
  kFT0M = 1,
  kFT0A = 2,
  kFT0C = 3
};

static const std::vector<std::string> centDetectorNames{"FV0A", "FT0M", "FT0A", "FT0C"};

// Event selections

enum evSel {
  kTVX = 0,
  kZvtx,
  kTFborder,
  kITSROFborder,
  kNoSameBunchPileup,
  kIsGoodZvtxFT0vsPV,
  kIsGoodITSLayersAll,
  kIsEPtriggered,
  kNevSels
};

static const std::vector<std::string> eventSelectionTitle{"Event selections"};
static const std::vector<std::string> eventSelectionLabels{"TVX", "Z vtx", "TF border", "ITS ROF border", "No same-bunch pile-up", "kIsGoodZvtxFT0vsPV", "isGoodITSLayersAll", "isEPtriggered"};

constexpr int EvSelDefault[8][1]{
  {1},
  {1},
  {0},
  {0},
  {0},
  {0},
  {0},
  {0}};

template <typename Tcollision> // move to nucleiUtils
bool eventSelection(const Tcollision& collision, o2::framework::HistogramRegistry& registry, o2::framework::LabeledArray<int> eventSelections, const float cutVertex)
{
  if (!registry.contains(HIST("hVtxZBefore"))) {
    registry.add("hVtxZBefore", "Vertex distribution in Z before selections;Z (cm)", {o2::framework::HistType::kTH1F, {{400, -20.0, 20.0}}});
  }
  if (!registry.contains(HIST("hVtxZ"))) {
    registry.add("hVtxZ", "Vertex distribution in Z;Z (cm)", {o2::framework::HistType::kTH1F, {{400, -20.0, 20.0}}});
  }
  if (!registry.contains(HIST("hEventSelections"))) {
    registry.add("hEventSelections", "hEventSelections", {o2::framework::HistType::kTH1D, {{evSel::kNevSels + 1, -0.5f, static_cast<float>(evSel::kNevSels) + 0.5f}}});
  }
  registry.fill(HIST("hEventSelections"), 0);
  registry.fill(HIST("hVtxZBefore"), collision.posZ());

  if (eventSelections.get(evSel::kTVX) && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
    return false;
  }
  registry.fill(HIST("hEventSelections"), evSel::kTVX + 1);

  if (eventSelections.get(evSel::kZvtx) && std::abs(collision.posZ()) > cutVertex) {
    return false;
  }
  registry.fill(HIST("hEventSelections"), evSel::kZvtx + 1);

  if (eventSelections.get(evSel::kTFborder) && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    return false;
  }
  registry.fill(HIST("hEventSelections"), evSel::kTFborder + 1);

  if (eventSelections.get(evSel::kITSROFborder) && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
    return false;
  }
  registry.fill(HIST("hEventSelections"), evSel::kITSROFborder + 1);

  if (eventSelections.get(evSel::kNoSameBunchPileup) && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    return false;
  }
  registry.fill(HIST("hEventSelections"), evSel::kNoSameBunchPileup + 1);

  if (eventSelections.get(evSel::kIsGoodZvtxFT0vsPV) && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    return false;
  }
  registry.fill(HIST("hEventSelections"), evSel::kIsGoodZvtxFT0vsPV + 1);

  if (eventSelections.get(evSel::kIsGoodITSLayersAll) && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
    return false;
  }
  registry.fill(HIST("hEventSelections"), evSel::kIsGoodITSLayersAll + 1);

  if constexpr (
    requires {
      collision.triggereventep();
    }) {
    if (eventSelections.get(evSel::kIsEPtriggered) && !collision.triggereventep()) {
      return false;
    }
    registry.fill(HIST("hEventSelections"), evSel::kIsEPtriggered + 1);
  }
  registry.fill(HIST("hVtxZ"), collision.posZ());

  return true;
}

/**
 * hFailCentrality is an histogram which is filled with 0. each time this functio is called,
 * then fills 1. if the centrality filling fails (return = -1.)
 */
template <typename Tcollision>
float getCentrality(Tcollision const& collision, const int centralityEstimator, std::shared_ptr<TH1> hFailCentrality = nullptr)
{
  if (hFailCentrality) {
    hFailCentrality->Fill(0.);
  }
  if constexpr (!o2::aod::HasCentrality<Tcollision>) { // requires aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs
    return -1.f;
    if (hFailCentrality) {
      hFailCentrality->Fill(1.);
    }
  }
  if (centralityEstimator == centDetectors::kFV0A) {
    return collision.centFV0A();
  } else if (centralityEstimator == centDetectors::kFT0M) {
    return collision.centFT0M();
  } else if (centralityEstimator == centDetectors::kFT0A) {
    return collision.centFT0A();
  } else if (centralityEstimator == centDetectors::kFT0C) {
    return collision.centFT0C();
  } else {
    LOG(warning) << "Centrality estimator not valid. Possible values: (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3). Centrality set to -1.";
  }
  if (hFailCentrality) {
    hFailCentrality->Fill(1.);
  }
  return -1.f;
}

// Track selections

enum trackSelection {
  kNoCuts = 0,
  kTrackCuts = 1,
  kPidCuts = 2,
  kNtrackSelections = 3
};
static const std::array<std::string, static_cast<int>(trackSelection::kNtrackSelections)> trackSelectionLabels{"All", "Track cuts", "PID cuts"};

template <int iSpecies>
void createHistogramRegistryNucleus(o2::framework::HistogramRegistry& registry)
{

  constexpr int index = iSpecies;
  if (!checkSpeciesValidity(index)) {
    std::runtime_error("species contains invalid nucleus index");
  }

  registry.add(fmt::format("{}/hTrackSelections", cNames[index]).c_str(), (fmt::format("{} track selections;", cNames[index]) + std::string("Selection step; Counts")).c_str(), o2::framework::HistType::kTH1D, {{trackSelection::kNtrackSelections, -0.5f, static_cast<float>(trackSelection::kNtrackSelections) - 0.5f}});
  registry.add(fmt::format("{}/hPtReconstructed", cNames[index]).c_str(), (fmt::format("{} - reconstructed variables;", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c}); Counts")).c_str(), o2::framework::HistType::kTH1F, {{400, -10.0f, 10.0f}});
  registry.add(fmt::format("{}/h2PtVsCentralityReconstructed", cNames[index]).c_str(), (fmt::format("{} - reconstructed variables;", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c}); CentralityFT0C (%)")).c_str(), o2::framework::HistType::kTH2F, {{400, -10.0f, 10.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3PhiVsEtaVsCentralityReconstructed", cNames[index]).c_str(), (fmt::format("{} - reconstructed variables;", cNames[index]) + std::string("#phi (radians); #eta; CentralityFT0C (%)")).c_str(), o2::framework::HistType::kTH3F, {{40, 0, o2::constants::math::TwoPI}, {40, -1.0f, 1.f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3DCAxyVsPtVsCentrality", cNames[index]).c_str(), (fmt::format(";", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c}); DCA_{xy} (cm); CentralityFT0C (%)")).c_str(), o2::framework::HistType::kTH3F, {{400, -10.0f, 10.0f}, {200, -0.5f, 0.5f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3DCAzVsPtVsCentrality", cNames[index]).c_str(), (fmt::format("{};", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c}); DCA_{z} (cm); CentralityFT0C (%)")).c_str(), o2::framework::HistType::kTH3F, {{400, -10.0f, 10.0f}, {200, -0.5f, 0.5f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3NsigmaTPC_preselectionVsCentrality", cNames[index]).c_str(), (fmt::format("Nsigma{} TPC distribution;", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c});") + fmt::format("n#sigma_{{TPC}}({}); CentralityFT0C (%)", cNames[index])).c_str(), o2::framework::HistType::kTH3F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3NsigmaTPCVsCentrality", cNames[index]).c_str(), (fmt::format("Nsigma{} TPC distribution;", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c});") + fmt::format("n#sigma_{{TPC}}({}); Centrality FT0C (%)", cNames[index])).c_str(), o2::framework::HistType::kTH3F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3NsigmaITS_preselectionVsCentrality", cNames[index]).c_str(), (fmt::format("Nsigma{} ITS distribution;", cNames[index]) + std::string("signed #it{p}_{T} / |#it{Z}| (GeV/#it{c});") + fmt::format("n#sigma_{{ITS}}({}); Centrality FT0C (%)", cNames[index])).c_str(), o2::framework::HistType::kTH3F, {{50, -5.0f, 5.0f}, {120, -3.0f, 3.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3NsigmaITSVsCentrality", cNames[index]).c_str(), (fmt::format("Nsigma{} ITS distribution;", cNames[index]) + std::string("signed #it{p}_{T} / |#it{Z}| (GeV/#it{c});") + fmt::format("n#sigma_{{ITS}}({}); Centrality FT0C (%)", cNames[index])).c_str(), o2::framework::HistType::kTH3F, {{50, -5.0f, 5.0f}, {120, -3.0f, 3.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3NsigmaTOF_preselectionVsCentrality", cNames[index]).c_str(), (fmt::format("Nsigma{} TOF distribution;", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c});") + fmt::format("n#sigma_{{TOF}}({}); Centrality FT0C (%)", cNames[index])).c_str(), o2::framework::HistType::kTH3F, {{100, -5.0f, 5.0f}, {400, -10.0f, 10.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3NsigmaTOFVsCentrality", cNames[index]).c_str(), (fmt::format("Nsigma{} TOF distribution;", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c});") + fmt::format("n#sigma_{{TOF}}({}); Centrality FT0C (%)", cNames[index], cNames[index])).c_str(), o2::framework::HistType::kTH3F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3BetaVsPtVsCentrality", cNames[index]).c_str(), (fmt::format("{};", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c}); #beta; CentralityFT0C (%)")).c_str(), o2::framework::HistType::kTH3F, {{400, -10.0f, 10.0f}, {100, 0.0f, 1.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3dEdxVsPVsCentrality", cNames[index]).c_str(), (fmt::format("dEdx distribution for {};", cNames[index]) + std::string("#it{p} (GeV/#it{c}); d#it{E}/d#it{x} (a.u.); Centrality FT0C (%)")).c_str(), o2::framework::HistType::kTH3F, {{200, -6.0f, 6.0f}, {100, 0.0f, 2000.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3ClusterSizeVsPtVsCentrality", cNames[index]).c_str(), (fmt::format("{};", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c}); Cluster size ITS; CentralityFT0C (%)")).c_str(), o2::framework::HistType::kTH3F, {{400, -10.0f, 10.0f}, {90, 0.f, 15.f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/hPtGenerated", cNames[index]).c_str(), (fmt::format("{} - generated variables;", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c}); Counts")).c_str(), o2::framework::HistType::kTH1F, {{400, -10.0f, 10.0f}});
  registry.add(fmt::format("{}/h2PtVsCentralityGenerated", cNames[index]).c_str(), (fmt::format("{} - generated variables;", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c}); CentralityFT0C (%)")).c_str(), o2::framework::HistType::kTH2F, {{400, -10.0f, 10.0f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3PtVsRapidityVsCentralityGenerated", cNames[index]).c_str(), (fmt::format("{} - generated variables;", cNames[index]) + std::string("#it{p}_{T} / |#it{Z}| (GeV/#it{c}); y; CentralityFT0C (%)")).c_str(), o2::framework::HistType::kTH3F, {{400, -10.0f, 10.0f}, {40, -1.0f, 1.f}, {20, 0.0f, 100.0f}});
  registry.add(fmt::format("{}/h3PhiVsEtaVsCentralityGenerated", cNames[index]).c_str(), (fmt::format("{} - generated variables;", cNames[index]) + std::string("#phi (radians); #eta; CentralityFT0C (%)")).c_str(), o2::framework::HistType::kTH3F, {{40, 0, o2::constants::math::TwoPI}, {40, -1.0f, 1.f}, {20, 0.0f, 100.0f}});

  for (size_t iSel = 0; iSel < trackSelection::kNtrackSelections; iSel++) {
    registry.get<TH1>(HIST(cNames[index]) + HIST("/hTrackSelections"))->GetXaxis()->SetBinLabel(iSel + 1, trackSelectionLabels[iSel].c_str());
  }
};

// PID manager class

class PidManager
{

 public:
  explicit PidManager(const int species, const float* tpcBetheBlochParams = nullptr)
    : mSpecies(species)
  {
    if (!checkSpeciesValidity(species)) {
      std::runtime_error("species contains invalid nucleus index");
    }

    if (!tpcBetheBlochParams) {
      mUseTpcCentralCalibration = true;
      return;
    }

    for (int i = 0; i < 6; i++) {
      mTpcBetheBlochParams[i] = tpcBetheBlochParams[i];
    }
  }
  PidManager() = default;
  ~PidManager() = default;

  // TOF
  template <typename Ttrack>
  float getBetaTOF(const Ttrack& track)
  {
    if (!track.hasTOF())
      return -999.f;
    float beta = o2::pid::tof::Beta::GetBeta(track);
    return std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
  }

  template <typename Ttrack>
  float getMassTOF(const Ttrack& track)
  {
    if (!track.hasTOF())
      return -999.f;
    const float charge{1.f + static_cast<float>(mSpecies == Species::kHe || mSpecies == Species::kAl)};
    const float beta = getBetaTOF(track);
    return track.tpcInnerParam() * charge * std::sqrt(1.f / (beta * beta) - 1.f);
  }

  template <typename Ttrack>
  float getNSigmaTOF(const Ttrack& track)
  {
    if (!track.hasTOF())
      return -999.f;

    switch (mSpecies) {
      case Species::kPr:
        return track.tofNSigmaPr();
      case Species::kDe:
        return track.tofNSigmaDe();
      case Species::kTr:
        return track.tofNSigmaTr();
      case Species::kHe:
        return track.tofNSigmaHe();
      case Species::kAl:
        return track.tofNSigmaAl();
      default:
        return -999.f;
    }
  }

  // ITS
  template <typename Ttrack>
  float getClusterSizeCosLambdaITS(const Ttrack& track)
  {
    return mResponseITS.averageClusterSize(track.itsClusterSizes()) / std::cosh(track.eta());
  }

  float getClusterSizeCosLambdaITS(const u_int32_t clusterSizesITS, const float eta)
  {
    return mResponseITS.averageClusterSize(clusterSizesITS) / std::cosh(eta);
  }

  template <typename Ttrack>
  float getNSigmaITS(const Ttrack& track)
  {
    switch (mSpecies) {
      case Species::kPr:
        return mResponseITS.nSigmaITS<o2::track::PID::Proton>(track.itsClusterSizes(), track.p(), track.eta());
      case Species::kDe:
        return mResponseITS.nSigmaITS<o2::track::PID::Deuteron>(track.itsClusterSizes(), track.p(), track.eta());
      case Species::kTr:
        return mResponseITS.nSigmaITS<o2::track::PID::Triton>(track.itsClusterSizes(), track.p(), track.eta());
      case Species::kHe:
        return mResponseITS.nSigmaITS<o2::track::PID::Helium3>(track.itsClusterSizes(), 2. * track.p(), track.eta());
      case Species::kAl:
        return mResponseITS.nSigmaITS<o2::track::PID::Alpha>(track.itsClusterSizes(), 2. * track.p(), track.eta());
      default:
        return -999.f;
    }
  }

  // TPC
  float getExpectedTPCsignal(const float p)
  {
    if (mUseTpcCentralCalibration)
      return -999.f;

    float pScaled = p * mMomScaling[0] + mMomScaling[1];
    float betaGamma = pScaled / masses[mSpecies];
    return o2::tpc::BetheBlochAleph(betaGamma,
                                    mTpcBetheBlochParams[0],
                                    mTpcBetheBlochParams[1],
                                    mTpcBetheBlochParams[2],
                                    mTpcBetheBlochParams[3],
                                    mTpcBetheBlochParams[4]);
  }

  template <typename Ttrack>
  float getNSigmaTPC(const Ttrack& track)
  {
    if (mUseTpcCentralCalibration) {
      return getNSigmaTPCcentral(track);
    }
    float expectedSignal = getExpectedTPCsignal(track.tpcInnerParam());
    float resolution = mTpcBetheBlochParams[5];
    return (track.tpcSignal() - expectedSignal) / (expectedSignal * resolution);
  }

 protected:
  // TPC
  template <typename Ttrack>
  float getNSigmaTPCcentral(const Ttrack& track)
  {
    switch (mSpecies) {
      case Species::kPr:
        return track.tpcNSigmaPr();
      case Species::kDe:
        return track.tpcNSigmaDe();
      case Species::kTr:
        return track.tpcNSigmaTr();
      case Species::kHe:
        return track.tpcNSigmaHe();
      case Species::kAl:
        return track.tpcNSigmaAl();
      default:
        return -999.f;
    }
  }

 private:
  float mTpcBetheBlochParams[6];
  bool mUseTpcCentralCalibration = true; // this just becomes a check for the null pointer in the parameters
  o2::aod::ITSResponse mResponseITS;
  float mMomScaling[2]{1., 0.};
  int mSpecies;
};

} // namespace nuclei

#endif // PWGLF_UTILS_NUCLEIUTILS_H_
