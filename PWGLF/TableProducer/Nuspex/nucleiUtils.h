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

#ifndef PWGLF_TABLEPRODUCER_NUSPEX_NUCLEIUTILS_H_
#define PWGLF_TABLEPRODUCER_NUSPEX_NUCLEIUTILS_H_

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct NucleusCandidate {
  int globalIndex;
  int collTrackIndex;
  float pt;
  float eta;
  float phi;
  float tpcInnerParam;
  float beta;
  float zVertex;
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
constexpr double bbMomScalingDefault[5][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr double betheBlochDefault[5][6]{
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-239.99, 1.155, 1.099, 1.137, 1.006, 0.09},
  {-321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09},
  {-586.66, 1.859, 4.435, 0.282, 3.201, 0.09}};
constexpr double nSigmaTPCdefault[5][2]{
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.}};
// constexpr double nSigmaTOFdefault[5][2]{
//   {-5., 5.},
//   {-5., 5.},
//   {-5., 5.},
//   {-5., 5.},
//   {-5., 5.}};
constexpr double DCAcutDefault[5][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr int TreeConfigDefault[5][2]{
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}};
constexpr int FlowHistDefault[5][1]{
  {0},
  {0},
  {0},
  {0},
  {0}};
constexpr int DCAHistDefault[5][2]{
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}};
constexpr double DownscalingDefault[5][1]{
  {1.},
  {1.},
  {1.},
  {1.},
  {1.}};
// constexpr bool storeTreesDefault[5]{false, false, false, false, false};
constexpr int species{5};
constexpr float charges[5]{1.f, 1.f, 1.f, 2.f, 2.f};
constexpr float masses[5]{MassProton, MassDeuteron, MassTriton, MassHelium3, MassAlpha};
static const std::vector<std::string> matter{"M", "A"};
static const std::vector<std::string> pidName{"TPC", "TOF"};
static const std::vector<std::string> names{"proton", "deuteron", "triton", "He3", "alpha"};
static const std::vector<std::string> treeConfigNames{"Filter trees", "Use TOF selection"};
static const std::vector<std::string> flowConfigNames{"Save flow hists"};
static const std::vector<std::string> DCAConfigNames{"Save DCA hist", "Matter/Antimatter"};
static const std::vector<std::string> nSigmaConfigName{"nsigma_min", "nsigma_max"};
static const std::vector<std::string> nDCAConfigName{"max DCAxy", "max DCAz"};
static const std::vector<std::string> DownscalingConfigName{"Fraction of kept candidates"};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> chargeLabelNames{"Positive", "Negative"};

float pidCutTPC[5][2]; //[species][lower/upper limit]

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
} // namespace nuclei

#endif // PWGLF_TABLEPRODUCER_NUSPEX_NUCLEIUTILS_H_
