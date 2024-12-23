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
/// \file hypKfRecoTask.cxx
/// \brief Hypernuclei rconstruction using KFParticle package
/// \author Janik Ditzel <jditzel@cern.ch> and Michael Hartung <mhartung@cern.ch>

#include <limits>
#include <vector>
#include <string>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DCAFitter/DCAFitterN.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFHypernucleiKfTables.h"
#include "TRandom.h"
#include "Common/DataModel/CollisionAssociationTables.h"

// KFParticle
#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=[name/macro]
#endif
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, o2::aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCDe, aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl, aod::pidTOFmass>;

//----------------------------------------------------------------------------------------------------------------

namespace
{
static const int nDaughterParticles = 6;
enum DAUGHTERS { kPion,
                 kProton,
                 kDeuteron,
                 kTriton,
                 kHe3,
                 kAlpha };

static const std::vector<std::string> particleNames{"pion", "proton", "deuteron", "triton", "helion", "alpha"};
static const std::vector<int> particlePdgCodes{211, 2212, o2::constants::physics::kDeuteron, o2::constants::physics::kTriton, o2::constants::physics::kHelium3, o2::constants::physics::kAlpha};
static const std::vector<double> particleMasses{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron, o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3, o2::constants::physics::MassAlpha};
static const std::vector<int> particleCharge{1, 1, 1, 1, 2, 2};

const int nBetheParams = 6;
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
constexpr double betheBlochDefault[nDaughterParticles][nBetheParams]{
  {13.611469, 3.598765, -0.021138, 2.039562, 0.651040, 0.09},    // pion
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // proton
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // deuteron
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // triton
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09},  // helion
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09}}; // alpha

const int nTrkSettings = 15;
static const std::vector<std::string> trackPIDsettingsNames{"useBBparams", "minITSnCls", "minTPCnCls", "maxTPCchi2", "maxITSchi2", "minRigidity", "maxRigidity", "maxTPCnSigma", "TOFrequiredabove", "minTOFmass", "maxTOFmass", "minDcaToPvXY", "minDcaToPvZ", "minITSclsSize", "maxITSclsSize"};
constexpr double trackPIDsettings[nDaughterParticles][nTrkSettings]{
  {0, 0, 60, 3.0, 5000, 0.15, 1.2, 2.5, -1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 70, 2.5, 5000, 0.20, 4.0, 3.0, -1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 70, 5.0, 5000, 0.50, 5.0, 3.0, -1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 70, 5.0, 5000, 0.50, 5.0, 3.0, -1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 75, 1.5, 5000, 0.50, 5.0, 3.0, -1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 70, 1.5, 5000, 0.50, 5.0, 3.0, -1, 0, 100, 0., 0., 0., 1000}};

static const int nHyperNuclei = 10;
static const std::vector<std::string> hyperNucNames{"L->p+pi", "3LH->3He+pi", "3LH->d+p+pi", "4LH->4He+pi", "4LH->t+p+pi", "4LHe->3He+p+pi", "5LHe->4He+p+pi", "5LHe->3He+d+pi", "custom1", "custom2"};
static const std::vector<std::string> hyperNucEnabledLb{"enabled"};
static const std::vector<std::string> reduceLb{"reduce factor"};
constexpr int hyperNucEnabled[nHyperNuclei][1]{{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};
constexpr float reduceFactor[nHyperNuclei][1]{{1.}, {1.}, {1.}, {1.}, {1.}, {1.}, {1.}, {1.}, {1.}, {1.}};
static const std::vector<std::string> hyperNucPdgLb{"PDG code"};
static constexpr int hyperNucPdgCodes[nHyperNuclei][1]{
  {3122},
  {o2::constants::physics::kHyperTriton},
  {o2::constants::physics::kHyperTriton},
  {o2::constants::physics::kHyperHydrogen4},
  {o2::constants::physics::kHyperHydrogen4},
  {o2::constants::physics::kHyperHelium4},
  {o2::constants::physics::kHyperHelium5},
  {o2::constants::physics::kHyperHelium5},
  {0},
  {0}};
static const std::vector<std::string> hyperNucDaughtersLb{"daughter1", "daughter2", "daughter3", "daughter4"};
static const std::string hyperNucDaughters[nHyperNuclei][4]{{"proton", "pion", "none", "none"}, {"helion", "pion", "none", "none"}, {"deuteron", "proton", "pion", "none"}, {"alpha", "pion", "none", "none"}, {"triton", "proton", "pion", "none"}, {"helion", "proton", "pion", "none"}, {"alpha", "proton", "pion", "none"}, {"helion", "deuteron", "pion", "none"}, {"none", "none", "none", "none"}, {"none", "none", "none", "none"}}; // NOLINT: runtime/string
static const std::string hyperNucSigns[nHyperNuclei][4]{{"+", "-", "", ""}, {"+", "-", "", ""}, {"+", "+", "-", ""}, {"+", "-", "", ""}, {"+", "+", "-", ""}, {"+", "+", "-", ""}, {"+", "+", "-", ""}, {"+", "+", "-", ""}, {"", "", "", ""}, {"", "", "", ""}};                                                                                                                                                                            // NOLINT: runtime/string
const int nSelPrim = 8;
static const std::vector<std::string> preSelectionPrimNames{"minMass", "maxMass", "minCt", "maxCt", "minCosPa", "maxDcaTracks", "maxDcaMotherToPvXY", "maxDcaMotherToPvZ"};
constexpr double preSelectionsPrimaries[nHyperNuclei][nSelPrim]{
  {1.00, 1.30, 0, 50, 0.90, 100., 2.0, 5.0},
  {2.96, 3.04, 0, 30, 0.99, 100., 1.5, 4.0},
  {2.96, 3.04, 0, 30, 0.99, 100., 1.5, 4.0},
  {3.87, 3.97, 0, 30, 0.95, 100., 2.0, 5.0},
  {3.87, 3.97, 0, 30, 0.95, 100., 2.0, 5.0},
  {3.85, 3.99, 0, 30, 0.98, 100., 1.5, 4.0},
  {4.60, 5.20, 0, 100, -1., 100., 10., 10.},
  {4.60, 5.20, 0, 100, -1., 100., 10., 10.},
  {0.00, 9.90, 0, 100, -1., 100., 10., 10.},
  {0.00, 9.90, 0, 100, -1., 100., 10., 10.}};
const int nSelSec = 8;
static const std::vector<std::string> preSelectionSecNames{"minMass", "maxMass", "minCt", "maxCt", "minCosPaSv", "maxDcaTracks", "maxDcaMotherToSvXY", "maxDcaMotherToSvZ"};
constexpr double preSelectionsSecondaries[nHyperNuclei][nSelSec]{
  {1.00, 1.30, 0, 50, 0.90, 100., 2.0, 5.0},
  {2.96, 3.04, 0, 30, 0.99, 100., 1.5, 4.0},
  {2.96, 3.04, 0, 30, 0.99, 100., 1.5, 4.0},
  {3.87, 3.97, 0, 30, 0.95, 100., 2.0, 5.0},
  {3.87, 3.97, 0, 30, 0.95, 100., 2.0, 5.0},
  {3.85, 3.99, 0, 30, 0.98, 100., 1.5, 4.0},
  {4.60, 5.20, 0, 100, -1., 100., 10., 10.},
  {4.60, 5.20, 0, 100, -1., 100., 10., 10.},
  {0.00, 9.90, 0, 100, -1., 100., 10., 10.},
  {0.00, 9.90, 0, 100, -1., 100., 10., 10.}};

static const int nCascades = 6;
static const std::vector<std::string> cascadeNames{"4LLH->4LHe+pi", "4XHe->4LHe+pi", "custom1", "custom2", "custom3", "custom4"};
constexpr int cascadeEnabled[nCascades][1]{{0}, {0}, {0}, {0}, {0}, {0}};
constexpr int cascadePdgCodes[nCascades][1]{
  {1020010040},
  {1120010040},
  {0},
  {0},
  {0},
  {0}};
static const std::vector<std::string> cascadeHypDaughterLb{"hypernucleus"};
static const std::string cascadeHypDaughter[nCascades][1]{{"4LHe->3He+p+pi"}, {"4LHe->3He+p+pi"}, {"none"}, {"none"}, {"none"}, {"none"}}; // NOLINT: runtime/string
static const std::vector<std::string> cascadeDaughtersLb{"daughter2", "daughter3", "daughter4"};
static const std::string cascadeDaughters[nCascades][3]{{"pion", "none", "none"}, {"pion", "none", "none"}, {"none", "none", "none"}, {"none", "none", "none"}, {"none", "none", "none"}, {"none", "none", "none"}}; // NOLINT: runtime/string
static const std::string cascadeSigns[nCascades][4]{{"+", "-", "", ""}, {"+", "-", "", ""}, {"", "", "", ""}, {"", "", "", ""}, {"", "", "", ""}, {"", "", "", ""}};                                                 // NOLINT: runtime/string
const int nSelCas = 8;
static const std::vector<std::string> preSelectionCascadeNames{"minMass", "maxMass", "minCt", "maxCt", "minCosPa", "maxDcaTracks", "maxDcaMotherToPvXY", "maxDcaMotherToPvZ"};
constexpr double preSelectionsCascades[nCascades][nSelCas]{
  {4.00, 4.20, 0, 30, 0.95, 100., 2.0, 5.},
  {4.00, 4.20, 0, 30, 0.95, 100., 2.0, 5.},
  {0.00, 9.90, 0, 100, -1., 100., 10., 10.},
  {0.00, 9.90, 0, 100, -1., 100., 10., 10.},
  {0.00, 9.90, 0, 100, -1., 100., 10., 10.},
  {0.00, 9.90, 0, 100, -1., 100., 10., 10.}};

//----------------------------------------------------------------------------------------------------------------
struct DaughterParticle {
  TString name;
  int pdgCode;
  double mass;
  int charge;
  double resolution;
  std::vector<double> betheParams;

  DaughterParticle(std::string name_, int pdgCode_, double mass_, int charge_, LabeledArray<double> bethe)
  {
    name = TString(name_);
    pdgCode = pdgCode_;
    mass = mass_;
    charge = charge_;
    resolution = bethe.get(name, "resolution");
    betheParams.clear();
    for (unsigned int i = 0; i < 5; i++)
      betheParams.push_back(bethe.get(name, i));
  }
}; // class DaughterParticle

struct HyperNucleus {
  TString name;
  int pdgCode;
  double massMax;
  double massMin;
  bool active;
  std::vector<int> daughters;
  std::vector<int> daughterTrackSigns;

  HyperNucleus(std::string name_, int pdgCode_, bool active_, std::vector<int> daughters_, std::vector<int> daughterTrackSigns_)
  {
    name = TString(name_);
    pdgCode = pdgCode_;
    active = active_;
    for (const int& d : daughters_)
      daughters.push_back(d);
    for (const int& dc : daughterTrackSigns_)
      daughterTrackSigns.push_back(dc);
  }
  HyperNucleus(std::string name_, int pdgCode_, bool active_, int hypDaughter, std::vector<int> daughters_, std::vector<int> daughterTrackSigns_)
  {
    daughters.push_back(hypDaughter);
    name = TString(name_);
    pdgCode = pdgCode_;
    active = active_;
    for (const int& d : daughters_)
      daughters.push_back(d);
    for (const int& dc : daughterTrackSigns_)
      daughterTrackSigns.push_back(dc);
  }
  int getNdaughters() { return static_cast<int>(daughters.size()); }
  const char* motherName() { return name.Contains("->") ? ((TString)name(0, name.First("-"))).Data() : name.Data(); }
  const char* daughterNames() { return name.Contains("->") ? ((TString)name(name.First("-") + 2, name.Length())).Data() : ""; }
}; // class HyperNucleus

struct HyperNucCandidate {
  int species;
  KFParticle kfp;
  std::vector<KFParticle> kfpDaughters;
  HyperNucCandidate* hypNucDaughter;
  std::vector<int64_t> daughterTrackIds;
  std::vector<float> recoSV;
  float devToVtx;
  float dcaToVtxXY;
  float dcaToVtxZ;
  float chi2;
  bool mcTrue;
  bool isPhysPrimary;
  bool isPrimaryCandidate, isSecondaryCandidate, isUsedSecondary;
  int64_t mcParticleId;
  int tableId;

  HyperNucCandidate(int species_, std::vector<KFParticle> kfpDaughters_, std::vector<int64_t> daughterTrackIds_) : species(species_), hypNucDaughter(0), devToVtx(999), dcaToVtxXY(999), dcaToVtxZ(999), chi2(999), mcTrue(false), isPhysPrimary(false), isPrimaryCandidate(false), isSecondaryCandidate(false), isUsedSecondary(false), mcParticleId(-1), tableId(-1)
  {
    for (const auto& kfd : kfpDaughters_)
      kfpDaughters.push_back(kfd);
    for (const auto& dt : daughterTrackIds_)
      daughterTrackIds.push_back(dt);
    init();
  }
  HyperNucCandidate(int species_, HyperNucCandidate* hypNucDaughter_, std::vector<KFParticle> kfpDaughters_, std::vector<int64_t> daughterTrackIds_) : species(species_), hypNucDaughter(hypNucDaughter_), devToVtx(999), dcaToVtxXY(999), dcaToVtxZ(999), chi2(999), mcTrue(false), isPhysPrimary(false), isPrimaryCandidate(false), isSecondaryCandidate(false), isUsedSecondary(false), mcParticleId(-1), tableId(-1)
  {
    for (const auto& kfd : kfpDaughters_)
      kfpDaughters.push_back(kfd);
    for (const auto& dt : daughterTrackIds_)
      daughterTrackIds.push_back(dt);
    init();
  }

  void init()
  {
    kfp.SetConstructMethod(2);
    for (size_t i = 0; i < kfpDaughters.size(); i++)
      kfp.AddDaughter(kfpDaughters.at(i));
    kfp.TransportToDecayVertex();
    chi2 = kfp.GetChi2() / kfp.GetNDF();
    recoSV.clear();
    recoSV.push_back(kfp.GetX());
    recoSV.push_back(kfp.GetY());
    recoSV.push_back(kfp.GetZ());
  }
  bool checkKfp()
  {
    if (kfp.GetMass() == 0)
      return false;
    if (std::isnan(kfp.GetMass()))
      return false;
    return true;
  }
  int getDaughterTableId()
  {
    if (hypNucDaughter)
      return hypNucDaughter->tableId;
    return -1;
  }
  bool isCascade() { return hypNucDaughter != 0; }
  int getSign()
  {
    if (kfp.GetQ() == 0)
      return kfpDaughters.front().GetQ() / std::abs(kfpDaughters.front().GetQ());
    return kfp.GetQ() / std::abs(kfp.GetQ());
  }
  int getNdaughters() { return static_cast<int>(kfpDaughters.size()); }
  float getDcaTracks() { return getNdaughters() == 2 ? getDcaTracks2() : getMaxDcaToSv(); }
  float getDcaTracks2() { return kfpDaughters.at(0).GetDistanceFromParticle(kfpDaughters.at(1)); }
  float getMaxDcaToSv()
  {
    float maxDca = std::numeric_limits<float>::lowest();
    for (const auto& daughter : kfpDaughters) {
      float dca = daughter.GetDistanceFromVertex(&recoSV[0]);
      if (dca > maxDca)
        maxDca = dca;
    }
    return maxDca;
  }
  float getDcaMotherToVertex(std::vector<float> vtx) { return kfp.GetDistanceFromVertex(&vtx[0]); }
  double getCpa(std::vector<float> vtx)
  {
    kfp.TransportToDecayVertex();
    return RecoDecay::cpa(std::array{vtx[0], vtx[1], vtx[2]}, std::array{recoSV[0], recoSV[1], recoSV[2]}, std::array{kfp.GetPx(), kfp.GetPy(), kfp.GetPz()});
    ;
  }
  float getCt(std::vector<float> vtx)
  {
    float dl = 0;
    for (size_t i = 0; i < vtx.size(); i++) {
      float tmp = recoSV.at(i) - vtx.at(i);
      dl += (tmp * tmp);
    }
    return std::sqrt(dl) * kfp.GetMass() / kfp.GetP();
  }
  float getDcaMotherToVtxXY(std::vector<float> vtx) { return kfp.GetDistanceFromVertexXY(&vtx[0]); }
  float getDcaMotherToVtxZ(std::vector<float> vtx)
  {
    kfp.TransportToPoint(&vtx[0]);
    return kfp.GetZ() - vtx[2];
  }
  void getDaughterPosMom(int daughter, std::vector<float>& posMom)
  {
    kfpDaughters.at(daughter).TransportToPoint(&recoSV[0]);
    posMom.assign({kfpDaughters.at(daughter).GetX(), kfpDaughters.at(daughter).GetY(), kfpDaughters.at(daughter).GetZ(), kfpDaughters.at(daughter).GetPx(), kfpDaughters.at(daughter).GetPy(), kfpDaughters.at(daughter).GetPz()});
  }
  void calcDevToVtx(KFPVertex& vtx) { devToVtx = kfp.GetDeviationFromVertexXY(vtx); }
  void calcDevToVtx(HyperNucCandidate& cand)
  {
    devToVtx = kfp.GetDeviationFromParticleXY(cand.kfp);
    dcaToVtxXY = getDcaMotherToVtxXY(cand.recoSV);
    dcaToVtxZ = getDcaMotherToVtxZ(cand.recoSV);
  }
  float getSubDaughterMass(int d1, int d2)
  {
    KFParticle subDaughter;
    subDaughter.SetConstructMethod(2);
    subDaughter.AddDaughter(kfpDaughters.at(d1));
    subDaughter.AddDaughter(kfpDaughters.at(d2));
    subDaughter.TransportToDecayVertex();
    return subDaughter.GetMass();
  }
  KFParticle getDaughterTrackKfp(int track)
  {
    return kfpDaughters.at(track + (isCascade() ? 1 : 0));
  }
  float getDcaTrackToVtxXY(int track, std::vector<float> vtx)
  {
    return getDaughterTrackKfp(track).GetDistanceFromVertexXY(&vtx[0]);
  }
  float getDcaTrackToVtxZ(int track, std::vector<float> vtx)
  {
    auto dKfp = getDaughterTrackKfp(track);
    dKfp.TransportToPoint(&vtx[0]);
    return dKfp.GetZ() - vtx[2];
  }
}; // class HyperNucCandidate

struct IndexPairs {
  std::vector<std::pair<int64_t, int>> pairs;

  void add(int64_t a, int b) { pairs.push_back({a, b}); }
  void clear() { pairs.clear(); }
  bool getIndex(int64_t a, int& b)
  {
    for (const auto& pair : pairs) {
      if (pair.first == a) {
        b = pair.second;
        return true;
      }
    }
    return false;
  }
}; // class IndexPairs

struct McCollInfo {
  bool hasRecoColl;
  bool passedEvSel;
  bool hasRecoParticle;
  int tableIndex;
  McCollInfo() : hasRecoColl(false), passedEvSel(false), hasRecoParticle(false), tableIndex(-1) {}
}; // class McCollInfo

//----------------------------------------------------------------------------------------------------------------
std::vector<std::shared_ptr<TH2>> hDeDx;
std::vector<std::shared_ptr<TH1>> hInvMass;
} // namespace

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
struct hypKfRecoTask { // o2-linter: disable=[name/workflow-file][name/struct]

  Produces<aod::HypKfMcColls> outputMcCollisionTable;
  Produces<aod::HypKfMcParts> outputMcParticleTable;
  Produces<aod::HypKfColls> outputCollisionTable;
  Produces<aod::HypKfTracks> outputTrackTable;
  Produces<aod::HypKfDaughtAdds> outputDaughterAddonTable;
  Produces<aod::HypKfSubDs> outputSubDaughterTable;
  Produces<aod::HypKfHypNucs> outputHypNucTable;

  Preslice<aod::TrackAssoc> perCollision = aod::track_association::collisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> cfgSaveOnlyMcTrue{"cfgSaveOnlyMcTrue", true, "save only MCtrue candidates"};
  Configurable<int> cfgDebug{"cfgDebug", 1, "debug level"};

  Configurable<bool> cfgRigidityCorrection{"cfgRigidityCorrection", false, "apply rigidity correction"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
  Configurable<bool> cfgUsePVcontributors{"cfgUsePVcontributors", true, "use tracks that are PV contibutors"};

  Configurable<LabeledArray<int>> cfgHyperNucsActive{"cfgHyperNucsActive", {hyperNucEnabled[0], nHyperNuclei, 1, hyperNucNames, hyperNucEnabledLb}, "enable or disable reconstruction"};
  Configurable<LabeledArray<float>> cfgReduce{"cfgReduce", {reduceFactor[0], nHyperNuclei, 1, hyperNucNames, reduceLb}, "reconstruct only a percentage of all possible hypernuclei"};
  Configurable<LabeledArray<int>> cfgHyperNucPdg{"cfgHyperNucPdg", {hyperNucPdgCodes[0], nHyperNuclei, 1, hyperNucNames, hyperNucPdgLb}, "PDG codes"};
  Configurable<LabeledArray<std::string>> cfgHyperNucDaughters{"cfgHyperNucDaughters", {hyperNucDaughters[0], nHyperNuclei, 4, hyperNucNames, hyperNucDaughtersLb}, "Daughter particles"};
  Configurable<LabeledArray<std::string>> cfgHyperNucSigns{"cfgHyperNucSigns", {hyperNucSigns[0], nHyperNuclei, 4, hyperNucNames, hyperNucDaughtersLb}, "Daughter signs"};

  Configurable<LabeledArray<int>> cfgCascadesActive{"cfgCascadesActive", {cascadeEnabled[0], nCascades, 1, cascadeNames, hyperNucEnabledLb}, "enable or disable reconstruction"};
  Configurable<LabeledArray<int>> cfgCascadesPdg{"cfgCascadesPdg", {cascadePdgCodes[0], nCascades, 1, cascadeNames, hyperNucPdgLb}, "PDG codes"};
  Configurable<LabeledArray<std::string>> cfgCascadeHypDaughter{"cfgCascadeHypDaughter", {cascadeHypDaughter[0], nCascades, 1, cascadeNames, cascadeHypDaughterLb}, "Hyernuclei daugther"};
  Configurable<LabeledArray<std::string>> cfgCascadeDaughters{"cfgCascadeDaughters", {cascadeDaughters[0], nCascades, 3, cascadeNames, cascadeDaughtersLb}, "Daughter particles"};
  Configurable<LabeledArray<std::string>> cfgCascadeSigns{"cfgCascadeSigns", {cascadeSigns[0], nCascades, 4, cascadeNames, hyperNucDaughtersLb}, "Daughter signs"};

  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], nDaughterParticles, nBetheParams, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings{"cfgTrackPIDsettings", {trackPIDsettings[0], nDaughterParticles, nTrkSettings, particleNames, trackPIDsettingsNames}, "track selection and PID criteria"};
  Configurable<LabeledArray<double>> cfgPreSelectionsPrimaries{"cfgPreSelectionsPrimaries", {preSelectionsPrimaries[0], nHyperNuclei, nSelPrim, hyperNucNames, preSelectionPrimNames}, "selection criteria for primary hypernuclei"};
  Configurable<LabeledArray<double>> cfgPreSelectionsSecondaries{"cfgPreSelectionsSecondaries", {preSelectionsSecondaries[0], nHyperNuclei, nSelSec, hyperNucNames, preSelectionSecNames}, "selection criteria for secondary hypernuclei"};
  Configurable<LabeledArray<double>> cfgPreSelectionsCascades{"cfgPreSelectionsCascades", {preSelectionsCascades[0], nCascades, nSelCas, cascadeNames, preSelectionCascadeNames}, "selection criteria for cascade hypernuclei"};

  // CCDB
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<double> bField{"bField", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  std::vector<DaughterParticle> daughterParticles;
  std::vector<std::vector<int64_t>> foundDaughters;
  std::vector<std::vector<HyperNucCandidate>> singleHyperNucCandidates;  // hypernuclei candidates
  std::vector<std::vector<HyperNucCandidate>> cascadeHyperNucCandidates; // cascade candidates
  std::vector<HyperNucleus> singleHyperNuclei;
  std::vector<HyperNucleus> cascadeHyperNuclei;
  std::vector<float> primVtx;
  std::vector<float> cents;
  std::vector<McCollInfo> mcCollInfos;
  IndexPairs trackIndices;
  IndexPairs mcPartIndices;
  KFPVertex kfPrimVtx;
  bool collHasCandidate, collHasMcTrueCandidate;
  bool collPassedEvSel;
  int64_t mcCollTableIndex;
  int mRunNumber, occupancy;
  float dBz;
  TRandom rand;
  //----------------------------------------------------------------------------------------------------------------

  void init(InitContext const&)
  {
    mRunNumber = 0;
    dBz = 0;
    rand.SetSeed(0);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    for (int i = 0; i < nDaughterParticles; i++) { // create daughterparticles
      daughterParticles.push_back(DaughterParticle(particleNames.at(i), particlePdgCodes.at(i), particleMasses.at(i), particleCharge.at(i), cfgBetheBlochParams));
    }
    for (unsigned int i = 0; i < nHyperNuclei; i++) { // create hypernuclei
      singleHyperNuclei.push_back(HyperNucleus(hyperNucNames.at(i), cfgHyperNucPdg->get(i, 0u), cfgHyperNucsActive->get(i, 0u), getDaughterVec(i, cfgHyperNucDaughters), getDaughterSignVec(i, cfgHyperNucSigns)));
    }
    for (unsigned int i = 0; i < nCascades; i++) { // create cascades
      cascadeHyperNuclei.push_back(HyperNucleus(cascadeNames.at(i), cfgCascadesPdg->get(i, 0u), cfgCascadesActive->get(i, 0u), getHypDaughterVec(i, cfgCascadeHypDaughter), getDaughterVec(i, cfgCascadeDaughters), getDaughterSignVec(i, cfgCascadeSigns)));
    }
    // define histogram axes
    const AxisSpec axisMagField{10, -10., 10., "magnetic field"};
    const AxisSpec axisNev{3, 0., 3., "Number of events"};
    const AxisSpec axisRigidity{4000, -10., 10., "#it{p}^{TPC}/#it{z}"};
    const AxisSpec axisdEdx{2000, 0, 2000, "d#it{E}/d#it{x}"};
    const AxisSpec axisInvMass{1000, 1, 6, "inv mass"};
    const AxisSpec axisCent{100, 0, 100, "centrality"};
    const AxisSpec axisVtxZ{100, -10, 10, "z"};
    // create histograms
    histos.add("histMagField", "histMagField", kTH1F, {axisMagField});
    histos.add("histNev", "histNev", kTH1F, {axisNev});
    histos.add("histVtxZ", "histVtxZ", kTH1F, {axisVtxZ});
    histos.add("histCentFT0A", "histCentFT0A", kTH1F, {axisCent});
    histos.add("histCentFT0C", "histCentFT0C", kTH1F, {axisCent});
    histos.add("histCentFT0M", "histCentFT0M", kTH1F, {axisCent});
    hDeDx.resize(2 * nDaughterParticles + 2);
    for (int i = 0; i < nDaughterParticles + 1; i++) {
      TString histName = i < nDaughterParticles ? daughterParticles[i].name : "all";
      hDeDx[2 * i] = histos.add<TH2>(Form("histdEdx_%s", histName.Data()), ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x}", HistType::kTH2F, {axisRigidity, axisdEdx});
      hDeDx[2 * i + 1] = histos.add<TH2>(Form("histdEdx_%s_Cuts", histName.Data()), ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x}", HistType::kTH2F, {axisRigidity, axisdEdx});
    }
    // create invariant mass histograms
    hInvMass.resize(nHyperNuclei + nCascades);
    int histCount = 0;
    std::vector<std::vector<HyperNucleus>> hypNucVectors = {singleHyperNuclei, cascadeHyperNuclei};
    for (size_t i = 0; i < hypNucVectors.size(); i++) {
      for (size_t j = 0; j < hypNucVectors.at(i).size(); j++) {
        if (hypNucVectors.at(i).at(j).active) {
          hInvMass[histCount] = histos.add<TH1>(Form("h%d_%s", histCount, hypNucVectors.at(i).at(j).motherName()), ";;Counts", HistType::kTH1F, {axisInvMass});
        }
        histCount++;
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------

  void findDaughterParticles(aod::TrackAssoc const& tracksByColl, TracksFull const& tracks)
  {
    // track loop, store daughter candidates in std::vector
    for (const auto& trackId : tracksByColl) {
      const auto& track = tracks.rawIteratorAt(trackId.trackId());
      filldedx(track, nDaughterParticles);
      if (std::abs(track.eta()) > cfgCutEta)
        continue;
      if (!cfgUsePVcontributors && track.isPVContributor())
        continue;
      for (size_t i = 0; i < daughterParticles.size(); i++) {
        if (track.tpcNClsFound() < cfgTrackPIDsettings->get(i, "minTPCnCls"))
          continue;
        if (track.tpcChi2NCl() > cfgTrackPIDsettings->get(i, "maxTPCchi2"))
          continue;
        if (track.itsNCls() < cfgTrackPIDsettings->get(i, "minITSnCls"))
          continue;
        if (track.itsChi2NCl() > cfgTrackPIDsettings->get(i, "maxITSchi2"))
          continue;
        if (std::abs(getTPCnSigma(track, daughterParticles.at(i))) > cfgTrackPIDsettings->get(i, "maxTPCnSigma"))
          continue;
        filldedx(track, i);
        if (std::abs(track.dcaXY()) < cfgTrackPIDsettings->get(i, "minDcaToPvXY"))
          continue;
        if (std::abs(track.dcaZ()) < cfgTrackPIDsettings->get(i, "minDcaToPvZ"))
          continue;
        if (getMeanItsClsSize(track) < cfgTrackPIDsettings->get(i, "minITSclsSize"))
          continue;
        if (getMeanItsClsSize(track) > cfgTrackPIDsettings->get(i, "maxITSclsSize"))
          continue;
        if (getRigidity(track) < cfgTrackPIDsettings->get(i, "minRigidity") || getRigidity(track) > cfgTrackPIDsettings->get(i, "maxRigidity"))
          continue;
        if (cfgTrackPIDsettings->get(i, "TOFrequiredabove") >= 0 && getRigidity(track) > cfgTrackPIDsettings->get(i, "TOFrequiredabove") && (track.mass() < cfgTrackPIDsettings->get(i, "minTOFmass") || track.mass() > cfgTrackPIDsettings->get(i, "maxTOFmass")))
          continue;
        foundDaughters.at(i).push_back(track.globalIndex());
      }
    } // track loop
  }

  //----------------------------------------------------------------------------------------------------------------
  void checkMCTrueTracks(aod::McTrackLabels const& trackLabels, aod::McParticles const&)
  {
    for (int i = 0; i < nDaughterParticles; i++) {
      auto& daughterVec = foundDaughters.at(i);
      for (auto it = daughterVec.begin(); it < daughterVec.end(); it++) {
        bool mcTrue = true;
        const auto& mcLab = trackLabels.rawIteratorAt(*it);
        if (!mcLab.has_mcParticle()) {
          mcTrue = false;
        } else {
          const auto& mcPart = mcLab.mcParticle_as<aod::McParticles>();
          if (std::abs(mcPart.pdgCode()) != daughterParticles.at(i).pdgCode) {
            mcTrue = false;
          }
          if (!mcPart.has_mothers()) {
            mcTrue = false;
          }
        }
        if (!mcTrue) {
          daughterVec.erase(it);
        }
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  void createKFHypernuclei(TracksFull const& tracks)
  {
    // loop over all hypernuclei that are to be reconstructed
    for (size_t hyperNucIter = 0; hyperNucIter < singleHyperNuclei.size(); hyperNucIter++) {
      HyperNucleus* hyperNuc = &(singleHyperNuclei.at(hyperNucIter));
      if (!hyperNuc->active)
        continue;
      int nDaughters = hyperNuc->getNdaughters();

      std::vector<std::vector<int64_t>::iterator> it;
      int nCombinations = 1;
      for (int i = 0; i < nDaughters; i++) {
        nCombinations *= foundDaughters.at(hyperNuc->daughters.at(i)).size();
        it.push_back(foundDaughters.at(hyperNuc->daughters.at(i)).begin());
      }
      if (!nCombinations)
        continue;
      while (it[0] != foundDaughters.at(hyperNuc->daughters.at(0)).end()) {

        // check for correct signs, avoid double usage of tracks
        bool passedChecks = true;
        int checkSign = 0;
        std::vector<int64_t> vec;
        for (int i = 0; i < nDaughters; i++) {
          const auto& daughterTrack = tracks.rawIteratorAt(*(it[i]));
          if (!i)
            checkSign = daughterTrack.sign();
          if (daughterTrack.sign() != checkSign * hyperNuc->daughterTrackSigns.at(i) || std::find(vec.begin(), vec.end(), *it[i]) != vec.end()) {
            passedChecks = false;
            break;
          }
          vec.push_back(*it[i]);
        }
        if (passedChecks && rand.Rndm() <= cfgReduce->get((unsigned int)hyperNucIter, 0u)) {
          // create daugther KFParticles
          std::vector<int64_t> daughterIds;
          std::vector<KFParticle> daughterKfps;
          for (int i = 0; i < nDaughters; i++) {
            const auto& daughterTrack = tracks.rawIteratorAt(*(it[i]));
            daughterIds.push_back(*(it[i]));
            auto daughterMass = daughterParticles.at(hyperNuc->daughters.at(i)).mass;
            auto daughterCharge = daughterParticles.at(hyperNuc->daughters.at(i)).charge;
            daughterKfps.push_back(createKFParticle(daughterTrack, daughterMass, daughterCharge));
          }

          HyperNucCandidate candidate(hyperNucIter, daughterKfps, daughterIds);
          bool isPrimCandidate = true, isSecCandidate = true;
          if (candidate.checkKfp()) {
            // apply pre selections
            candidate.calcDevToVtx(kfPrimVtx);
            if (candidate.kfp.GetMass() < cfgPreSelectionsPrimaries->get(hyperNucIter, "minMass") || candidate.kfp.GetMass() > cfgPreSelectionsPrimaries->get(hyperNucIter, "maxMass"))
              isPrimCandidate = false;
            if (candidate.getDcaTracks() > cfgPreSelectionsPrimaries->get(hyperNucIter, "maxDcaTracks"))
              isPrimCandidate = false;
            if (candidate.getCt(primVtx) < cfgPreSelectionsPrimaries->get(hyperNucIter, "minCt") || candidate.getCt(primVtx) > cfgPreSelectionsPrimaries->get(hyperNucIter, "maxCt"))
              isPrimCandidate = false;
            if (candidate.getCpa(primVtx) < cfgPreSelectionsPrimaries->get(hyperNucIter, "minCosPa"))
              isPrimCandidate = false;
            if (std::abs(candidate.getDcaMotherToVtxXY(primVtx)) > cfgPreSelectionsPrimaries->get(hyperNucIter, "maxDcaMotherToPvXY"))
              isPrimCandidate = false;
            if (std::abs(candidate.getDcaMotherToVtxZ(primVtx)) > cfgPreSelectionsPrimaries->get(hyperNucIter, "maxDcaMotherToPvZ"))
              isPrimCandidate = false;
            if (isPrimCandidate) {
              candidate.isPrimaryCandidate = true;
              collHasCandidate = true;
            }
            if (candidate.kfp.GetMass() < cfgPreSelectionsSecondaries->get(hyperNucIter, "minMass") || candidate.kfp.GetMass() > cfgPreSelectionsSecondaries->get(hyperNucIter, "maxMass"))
              isSecCandidate = false;
            if (candidate.getDcaTracks() > cfgPreSelectionsSecondaries->get(hyperNucIter, "maxDcaTracks"))
              isSecCandidate = false;
            if (candidate.getCt(primVtx) < cfgPreSelectionsSecondaries->get(hyperNucIter, "minCt") || candidate.getCt(primVtx) > cfgPreSelectionsSecondaries->get(hyperNucIter, "maxCt"))
              isSecCandidate = false;
            if (isSecCandidate) {
              candidate.isSecondaryCandidate = true;
            }
            if (isPrimCandidate || isSecCandidate)
              singleHyperNucCandidates.at(hyperNucIter).push_back(candidate);
          }
        }
        it[nDaughters - 1]++;
        for (int i = nDaughters - 1; i && it[i] == foundDaughters.at(hyperNuc->daughters.at(i)).end(); i--) {
          it[i] = foundDaughters.at(hyperNuc->daughters.at(i)).begin();
          it[i - 1]++;
        }
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  void createKFCascades(TracksFull const& tracks)
  {

    // loop over all cascade hypernuclei that are to be reconstructed
    for (size_t hyperNucIter = 0; hyperNucIter < cascadeHyperNuclei.size(); hyperNucIter++) {
      HyperNucleus* hyperNuc = &(cascadeHyperNuclei.at(hyperNucIter));
      if (!hyperNuc->active)
        continue;
      int nDaughters = hyperNuc->getNdaughters();

      int nHypNucDaughters = singleHyperNucCandidates.at(hyperNuc->daughters.at(0)).size();
      std::vector<int64_t> vecHypNucDaughers;
      for (int64_t i = 0; i < static_cast<int64_t>(nHypNucDaughters); i++) {
        vecHypNucDaughers.push_back(i);
      }

      std::vector<std::vector<int64_t>::iterator> it;
      int nCombinations = 1;
      nCombinations *= nHypNucDaughters;
      it.push_back(vecHypNucDaughers.begin());
      for (int i = 1; i < nDaughters; i++) {
        nCombinations *= foundDaughters.at(hyperNuc->daughters.at(i)).size();
        it.push_back(foundDaughters.at(hyperNuc->daughters.at(i)).begin());
      }
      if (!nCombinations)
        continue;
      while (it[0] != vecHypNucDaughers.end()) {
        std::vector<int64_t> daughterIds;
        std::vector<KFParticle> daughterKfps;

        // select hypernuclei daughter KFParticle
        auto hypNucDaughter = &(singleHyperNucCandidates.at(hyperNuc->daughters.at(0)).at(*it[0]));
        // check for correct signs
        int checkSign = hypNucDaughter->getSign();
        bool passedChecks = true;
        std::vector<int64_t> vec = hypNucDaughter->daughterTrackIds;
        for (int i = 1; i < nDaughters; i++) {
          const auto& daughterTrack = tracks.rawIteratorAt(*(it[i]));
          if (!i)
            checkSign = daughterTrack.sign();
          if (daughterTrack.sign() != checkSign * hyperNuc->daughterTrackSigns.at(i) || std::find(vec.begin(), vec.end(), *it[i]) != vec.end()) {
            passedChecks = false;
            break;
          }
          vec.push_back(*it[i]);
        }
        if (passedChecks && hypNucDaughter->isSecondaryCandidate) {
          daughterKfps.push_back(hypNucDaughter->kfp);
          for (int i = 1; i < nDaughters; i++) {
            daughterIds.push_back(*(it[i]));
            const auto& daughterTrack = tracks.rawIteratorAt(*(it[i]));
            auto daughterMass = daughterParticles.at(hyperNuc->daughters.at(i)).mass;
            auto daughterCharge = daughterParticles.at(hyperNuc->daughters.at(i)).charge;
            daughterKfps.push_back(createKFParticle(daughterTrack, daughterMass, daughterCharge));
          }

          HyperNucCandidate candidate(hyperNucIter, hypNucDaughter, daughterKfps, daughterIds);
          if (candidate.checkKfp()) {
            hypNucDaughter->calcDevToVtx(candidate);
            bool isCandidate = true;
            // apply pre selections for hypernucleus daughter
            if (hypNucDaughter->getCpa(candidate.recoSV) < cfgPreSelectionsSecondaries->get(hyperNuc->daughters.at(0), "minCosPaSv"))
              isCandidate = false;
            if (hypNucDaughter->getDcaMotherToVtxXY(candidate.recoSV) > cfgPreSelectionsSecondaries->get(hyperNuc->daughters.at(0), "maxDcaMotherToSvXY"))
              isCandidate = false;
            if (hypNucDaughter->getDcaMotherToVtxZ(candidate.recoSV) > cfgPreSelectionsSecondaries->get(hyperNuc->daughters.at(0), "maxDcaMotherToSvZ"))
              isCandidate = false;
            // apply pre selections for cascade
            if (candidate.kfp.GetMass() < cfgPreSelectionsCascades->get(hyperNucIter, "minMass") || candidate.kfp.GetMass() > cfgPreSelectionsCascades->get(hyperNucIter, "maxMass"))
              isCandidate = false;
            if (candidate.getDcaTracks() > cfgPreSelectionsCascades->get(hyperNucIter, "maxDcaTracks"))
              isCandidate = false;
            if (candidate.getCt(primVtx) < cfgPreSelectionsCascades->get(hyperNucIter, "minCt") || candidate.getCt(primVtx) > cfgPreSelectionsCascades->get(hyperNucIter, "maxCt"))
              isCandidate = false;
            if (candidate.getCpa(primVtx) < cfgPreSelectionsCascades->get(hyperNucIter, "minCosPa"))
              isCandidate = false;
            if (std::abs(candidate.getDcaMotherToVtxXY(primVtx)) > cfgPreSelectionsCascades->get(hyperNucIter, "maxDcaMotherToPvXY"))
              isCandidate = false;
            if (std::abs(candidate.getDcaMotherToVtxZ(primVtx)) > cfgPreSelectionsCascades->get(hyperNucIter, "maxDcaMotherToPvZ"))
              isCandidate = false;

            if (isCandidate) {
              collHasCandidate = true;
              hypNucDaughter->isUsedSecondary = true;
              cascadeHyperNucCandidates.at(hyperNucIter).push_back(candidate);
            }
          }
        }
        it[nDaughters - 1]++;
        for (int i = nDaughters - 1; i && it[i] == foundDaughters.at(hyperNuc->daughters.at(i)).end(); i--) {
          it[i] = foundDaughters.at(hyperNuc->daughters.at(i)).begin();
          it[i - 1]++;
        }
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  void createMCinfo(aod::McTrackLabels const& trackLabels, aod::McCollisionLabels const&, aod::McParticles const& particlesMC, aod::McCollisions const&, bool cascadesOnly = false)
  {
    // check for mcTrue: single (primary & cascade daughter) and cascade hypernuclei

    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei, &cascadeHyperNuclei};
    std::vector<std::vector<std::vector<HyperNucCandidate>>*> candidateVectors = {&singleHyperNucCandidates, &cascadeHyperNucCandidates};
    const int nVecs = candidateVectors.size();
    const int startVec = cascadesOnly ? 1 : 0;
    for (int vec = startVec; vec < nVecs; vec++) {
      auto candidateVector = candidateVectors.at(vec);
      for (size_t hyperNucIter = 0; hyperNucIter < hypNucVectors.at(vec)->size(); hyperNucIter++) {
        HyperNucleus* hyperNuc = &(hypNucVectors.at(vec)->at(hyperNucIter));
        if (!hyperNuc->active)
          continue;
        for (auto& hypCand : candidateVector->at(hyperNucIter)) { // o2-linter: disable=[const-ref-in-for-loop]
          std::vector<int64_t> motherIds;
          int daughterCount = 0;
          if (hypCand.isCascade()) {
            if (!hypCand.hypNucDaughter->mcTrue)
              continue;
            const auto& mcPart = particlesMC.rawIteratorAt(hypCand.hypNucDaughter->mcParticleId);
            if (!mcPart.has_mothers())
              continue;
            daughterCount++;
            for (const auto& mother : mcPart.mothers_as<aod::McParticles>()) {
              if (mother.pdgCode() == hyperNuc->pdgCode * hypCand.getSign()) {
                motherIds.push_back(mother.globalIndex());
                break;
              }
            }
          }
          for (const auto& daughter : hypCand.daughterTrackIds) {
            const auto& mcLab = trackLabels.rawIteratorAt(daughter);
            if (!mcLab.has_mcParticle())
              continue;
            const auto& mcPart = mcLab.mcParticle_as<aod::McParticles>();
            if (std::abs(mcPart.pdgCode()) != daughterParticles.at(hyperNuc->daughters.at(daughterCount++)).pdgCode)
              continue;
            if (!mcPart.has_mothers())
              continue;
            for (const auto& mother : mcPart.mothers_as<aod::McParticles>()) {
              if (mother.pdgCode() == hyperNuc->pdgCode * hypCand.getSign()) {
                motherIds.push_back(mother.globalIndex());
                break;
              }
            }
          }
          if (motherIds.size() != hyperNuc->daughters.size()) {
            if (cfgSaveOnlyMcTrue)
              hypCand.isSecondaryCandidate = false;
            continue;
          }
          hypCand.mcTrue = true;
          for (auto iter = motherIds.begin(); iter != motherIds.end() - 1; iter++)
            if (*iter != *(iter + 1))
              hypCand.mcTrue = false;
          if (hypCand.mcTrue) {
            hypCand.mcParticleId = motherIds.front();
            collHasMcTrueCandidate = true;
          }
          if (!hypCand.mcTrue && cfgSaveOnlyMcTrue)
            hypCand.isSecondaryCandidate = false;
        }
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------

  void fillTree(TracksFull const& tracks, bool saveOnlyMcTrue = false)
  {

    outputCollisionTable(
      collPassedEvSel, mcCollTableIndex,
      primVtx.at(0), primVtx.at(1), primVtx.at(2),
      cents.at(0), cents.at(1), cents.at(2), occupancy);

    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei, &cascadeHyperNuclei};
    std::vector<std::vector<std::vector<HyperNucCandidate>>*> candidateVectors = {&singleHyperNucCandidates, &cascadeHyperNucCandidates};

    for (int vec = 0; vec < 2; vec++) {
      auto candidateVector = candidateVectors.at(vec);
      for (size_t hyperNucIter = 0; hyperNucIter < hypNucVectors.at(vec)->size(); hyperNucIter++) {
        HyperNucleus* hyperNuc = &(hypNucVectors.at(vec)->at(hyperNucIter));
        if (!hyperNuc->active)
          continue;
        for (auto& hypCand : candidateVector->at(hyperNucIter)) { // o2-linter: disable=[const-ref-in-for-loop]
          if (!hypCand.isPrimaryCandidate && !hypCand.isUsedSecondary && !hypCand.isCascade())
            continue;
          if (saveOnlyMcTrue && !hypCand.mcTrue)
            continue;
          hInvMass[vec * nHyperNuclei + hyperNucIter]->Fill(hypCand.kfp.GetMass());
          std::vector<int> vecDaugtherTracks, vecAddons, vecSubDaughters;
          int daughterCount = 0;
          for (const auto& daughterTrackId : hypCand.daughterTrackIds) {
            int trackTableId;
            if (!trackIndices.getIndex(daughterTrackId, trackTableId)) {
              auto daught = hyperNuc->daughters.at(daughterCount);
              const auto& track = tracks.rawIteratorAt(daughterTrackId);
              outputTrackTable(
                hyperNuc->daughters.at(daughterCount) * track.sign(),
                track.pt(), track.eta(), track.phi(),
                hypCand.getDcaTrackToVtxXY(daughterCount, primVtx), hypCand.getDcaTrackToVtxZ(daughterCount, primVtx),
                track.tpcNClsFound(), track.tpcChi2NCl(),
                track.itsClusterSizes(), track.itsChi2NCl(),
                getRigidity(track), track.tpcSignal(), getTPCnSigma(track, daughterParticles.at(daught)),
                daught == kAlpha ? -999 : getTPCnSigma(track, daughterParticles.at(daught + 1)),
                daught == kPion ? 999 : getTPCnSigma(track, daughterParticles.at(daught - 1)),
                track.mass(),
                track.isPVContributor());
              trackTableId = outputTrackTable.lastIndex();
              trackIndices.add(daughterTrackId, trackTableId);
            }
            vecDaugtherTracks.push_back(trackTableId);
            daughterCount++;
          }
          for (int i = 0; i < hypCand.getNdaughters(); i++) {
            std::vector<float> posMom;
            hypCand.getDaughterPosMom(i, posMom);
            outputDaughterAddonTable(
              posMom.at(0), posMom.at(1), posMom.at(2), posMom.at(3), posMom.at(4), posMom.at(5));
            vecAddons.push_back(outputDaughterAddonTable.lastIndex());
          }
          if (hypCand.getNdaughters() > 2) {
            for (int i = 0; i < hypCand.getNdaughters(); i++) {
              for (int j = i + 1; j < hypCand.getNdaughters(); j++) {
                outputSubDaughterTable(hypCand.getSubDaughterMass(i, j));
                vecSubDaughters.push_back(outputSubDaughterTable.lastIndex());
              }
            }
          }

          hypCand.kfp.TransportToDecayVertex();
          int mcPartTableId;
          outputHypNucTable(
            mcPartIndices.getIndex(hypCand.mcParticleId, mcPartTableId) ? mcPartTableId : -1,
            outputCollisionTable.lastIndex(), vecDaugtherTracks, vecAddons, hypCand.getDaughterTableId(), vecSubDaughters,
            (vec * nHyperNuclei + hyperNucIter + 1) * hypCand.getSign(),
            hypCand.isPrimaryCandidate, hypCand.kfp.GetMass(),
            hypCand.kfp.GetPx(), hypCand.kfp.GetPy(), hypCand.kfp.GetPz(),
            hypCand.getDcaMotherToVtxXY(primVtx), hypCand.getDcaMotherToVtxZ(primVtx),
            hypCand.devToVtx, hypCand.dcaToVtxXY, hypCand.dcaToVtxZ, hypCand.chi2,
            hypCand.recoSV.at(0), hypCand.recoSV.at(1), hypCand.recoSV.at(2));
          hypCand.tableId = outputHypNucTable.lastIndex();
        }
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------

  void processMC(CollisionsFullMC const& collisions, aod::McCollisions const& mcColls, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const& particlesMC, aod::McTrackLabels const& trackLabelsMC, aod::McCollisionLabels const& collLabels, aod::TrackAssoc const& tracksColl)
  {

    mcCollInfos.clear();
    mcCollInfos.resize(mcColls.size());
    mcPartIndices.clear();
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision())
        continue;
      if (collision.sel8() && std::abs(collision.posZ()) < 10)
        mcCollInfos.at(collision.mcCollisionId()).passedEvSel = true;
    }
    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei, &cascadeHyperNuclei};
    for (const auto& mcPart : particlesMC) {
      for (int vec = 0; vec < 2; vec++) {
        for (size_t hyperNucIter = 0; hyperNucIter < hypNucVectors.at(vec)->size(); hyperNucIter++) {
          HyperNucleus* hyperNuc = &(hypNucVectors.at(vec)->at(hyperNucIter));
          if (!hyperNuc->active)
            continue;
          if (std::abs(mcPart.pdgCode()) != hyperNuc->pdgCode)
            continue;
          bool isDecayMode = false;
          float svx, svy, svz;
          int daughterPdg;
          if (vec == 0)
            daughterPdg = daughterParticles.at(hyperNuc->daughters.at(0)).pdgCode;
          else
            daughterPdg = singleHyperNuclei.at(hyperNuc->daughters.at(0)).pdgCode;
          for (const auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
            if (std::abs(mcDaught.pdgCode()) == daughterPdg) {
              isDecayMode = true;
              svx = mcDaught.vx();
              svy = mcDaught.vy();
              svz = mcDaught.vz();
              break;
            }
          }
          if (!isDecayMode)
            continue;

          if (mcCollInfos.at(mcPart.mcCollisionId()).tableIndex < 0) {
            outputMcCollisionTable(
              mcCollInfos.at(mcPart.mcCollisionId()).passedEvSel,
              mcPart.mcCollision().posX(), mcPart.mcCollision().posY(), mcPart.mcCollision().posZ());
          }
          mcCollInfos.at(mcPart.mcCollisionId()).tableIndex = outputMcCollisionTable.lastIndex();

          outputMcParticleTable(
            mcCollInfos.at(mcPart.mcCollisionId()).tableIndex,
            (vec * nHyperNuclei + hyperNucIter + 1) * (mcPart.pdgCode() > 0 ? +1 : -1),
            mcPart.pdgCode(),
            mcPart.isPhysicalPrimary(),
            mcPart.px(), mcPart.py(), mcPart.pz(),
            mcPart.e(),
            svx, svy, svz);
          mcPartIndices.add(mcPart.globalIndex(), outputMcParticleTable.lastIndex());
        }
      }
    }

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      initCollision(collision);
      const uint64_t collIdx = collision.globalIndex();
      auto tracksByColl = tracksColl.sliceBy(perCollision, collIdx);
      findDaughterParticles(tracksByColl, tracks);
      if (cfgSaveOnlyMcTrue)
        checkMCTrueTracks(trackLabelsMC, particlesMC);
      createKFHypernuclei(tracks);
      createMCinfo(trackLabelsMC, collLabels, particlesMC, mcColls);
      createKFCascades(tracks);
      createMCinfo(trackLabelsMC, collLabels, particlesMC, mcColls, true);

      if (!collHasCandidate)
        continue;
      if (cfgSaveOnlyMcTrue && !collHasMcTrueCandidate)
        continue;

      mcCollTableIndex = -1;
      if (collision.has_mcCollision()) {
        mcCollTableIndex = mcCollInfos.at(collision.mcCollisionId()).tableIndex;
        if (mcCollTableIndex < 0) {
          outputMcCollisionTable(
            mcCollInfos.at(collision.mcCollisionId()).passedEvSel,
            collision.mcCollision().posX(), collision.mcCollision().posY(), collision.mcCollision().posZ());
          mcCollTableIndex = outputMcCollisionTable.lastIndex();
          mcCollInfos.at(collision.mcCollisionId()).tableIndex = mcCollTableIndex;
        }
      }
      fillTree(tracks, cfgSaveOnlyMcTrue);
    }
  }
  PROCESS_SWITCH(hypKfRecoTask, processMC, "MC analysis", false);
  //----------------------------------------------------------------------------------------------------------------
  void processData(CollisionsFull const& collisions, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::TrackAssoc const& tracksColl)
  {

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      initCollision(collision);
      if (!collPassedEvSel)
        continue;
      const uint64_t collIdx = collision.globalIndex();
      auto tracksByColl = tracksColl.sliceBy(perCollision, collIdx);
      findDaughterParticles(tracksByColl, tracks);
      createKFHypernuclei(tracks);
      createKFCascades(tracks);
      if (!collHasCandidate)
        continue;
      mcCollTableIndex = -1;
      fillTree(tracks);
    }
  }
  PROCESS_SWITCH(hypKfRecoTask, processData, "data analysis", true);
  //----------------------------------------------------------------------------------------------------------------
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grpTimestamp = bc.timestamp();
    dBz = 0;
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grpTimestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (bField < -990) {
        // Fetch magnetic field from ccdb for current collision
        dBz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
      } else {
        dBz = bField;
      }
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grpTimestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (bField < -990) {
        // Fetch magnetic field from ccdb for current collision
        dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
      } else {
        dBz = bField;
      }
    }
    mRunNumber = bc.runNumber();
    KFParticle::SetField(dBz);
  }
  //----------------------------------------------------------------------------------------------------------------
  template <typename T>
  void initCollision(const T& collision)
  {
    foundDaughters.clear();
    foundDaughters.resize(nDaughterParticles);
    singleHyperNucCandidates.clear();
    singleHyperNucCandidates.resize(nHyperNuclei);
    cascadeHyperNucCandidates.clear();
    cascadeHyperNucCandidates.resize(nCascades);
    trackIndices.clear();
    collHasCandidate = false;
    collHasMcTrueCandidate = false;
    histos.fill(HIST("histMagField"), dBz);
    histos.fill(HIST("histNev"), 0.5);
    collPassedEvSel = collision.sel8() && std::abs(collision.posZ()) < 10;
    if (collPassedEvSel) {
      histos.fill(HIST("histNev"), 1.5);
      histos.fill(HIST("histVtxZ"), collision.posZ());
      histos.fill(HIST("histCentFT0A"), collision.centFT0A());
      histos.fill(HIST("histCentFT0C"), collision.centFT0C());
      histos.fill(HIST("histCentFT0M"), collision.centFT0M());
    }
    occupancy = collision.trackOccupancyInTimeRange();
    kfPrimVtx = createKFPVertexFromCollision(collision);
    primVtx.assign({collision.posX(), collision.posY(), collision.posZ()});
    cents.assign({collision.centFT0A(), collision.centFT0C(), collision.centFT0M()});
  }

  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void filldedx(T const& track, int species)
  {
    const float rigidity = getRigidity(track);
    hDeDx[2 * species]->Fill(track.sign() * rigidity, track.tpcSignal());
    if (track.tpcNClsFound() < 100 || track.itsNCls() < 2)
      return;
    hDeDx[2 * species + 1]->Fill(track.sign() * rigidity, track.tpcSignal());
  }
  //----------------------------------------------------------------------------------------------------------------

  template <class T>
  float getTPCnSigma(T const& track, DaughterParticle const& particle)
  {
    const float rigidity = getRigidity(track);
    if (!track.hasTPC())
      return -999;

    if (particle.name == "pion" && cfgTrackPIDsettings->get("pion", "useBBparams") == 0)
      return track.tpcNSigmaPi();
    if (particle.name == "proton" && cfgTrackPIDsettings->get("proton", "useBBparams") == 0)
      return track.tpcNSigmaPr();
    if (particle.name == "deuteron" && cfgTrackPIDsettings->get("deuteron", "useBBparams") == 0)
      return track.tpcNSigmaDe();
    if (particle.name == "triton" && cfgTrackPIDsettings->get("triton", "useBBparams") == 0)
      return track.tpcNSigmaTr();
    if (particle.name == "helion" && cfgTrackPIDsettings->get("helion", "useBBparams") == 0)
      return track.tpcNSigmaHe();
    if (particle.name == "alpha" && cfgTrackPIDsettings->get("alpha", "useBBparams") == 0)
      return track.tpcNSigmaAl();

    double expBethe{tpc::BetheBlochAleph(static_cast<double>(particle.charge * rigidity / particle.mass), particle.betheParams[0], particle.betheParams[1], particle.betheParams[2], particle.betheParams[3], particle.betheParams[4])};
    double expSigma{expBethe * particle.resolution};
    float sigmaTPC = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
    return sigmaTPC;
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getMeanItsClsSize(T const& track)
  {
    int sum = 0, n = 0;
    for (int i = 0; i < 8; i++) {
      sum += (track.itsClusterSizes() >> (4 * i) & 15);
      if (track.itsClusterSizes() >> (4 * i) & 15)
        n++;
    }
    return n > 0 ? static_cast<float>(sum) / n : 0.f;
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getRigidity(T const& track)
  {
    if (!cfgRigidityCorrection)
      return track.tpcInnerParam();
    bool hePID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    return hePID ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
  }
  //----------------------------------------------------------------------------------------------------------------

  template <typename T>
  KFParticle createKFParticle(const T& track, float mass, int charge)
  {
    auto trackparCov = getTrackParCov(track);
    std::array<float, 3> fP;
    std::array<float, 3> fM;
    trackparCov.getXYZGlo(fP);
    trackparCov.getPxPyPzGlo(fM);
    float fPM[6];
    for (int i = 0; i < 3; i++) {
      fPM[i] = fP[i];
      fPM[i + 3] = fM[i] * std::abs(charge);
    }
    std::array<float, 21> fC;
    trackparCov.getCovXYZPxPyPzGlo(fC);
    KFParticle part;
    part.Create(fPM, fC.data(), std::abs(charge) * track.sign(), mass);
    return part;
  }
  //----------------------------------------------------------------------------------------------------------------

  template <typename T>
  KFPVertex createKFPVertexFromCollision(const T& collision)
  {
    KFPVertex kfpVertex;
    kfpVertex.SetXYZ(collision.posX(), collision.posY(), collision.posZ());
    kfpVertex.SetCovarianceMatrix(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    kfpVertex.SetChi2(collision.chi2());
    kfpVertex.SetNDF(2 * collision.numContrib() - 3);
    kfpVertex.SetNContributors(collision.numContrib());
    return kfpVertex;
  }
  //----------------------------------------------------------------------------------------------------------------

  int getHypDaughterVec(unsigned int cascade, LabeledArray<std::string> cfg)
  {
    std::string daughter = cfg.get(cascade, 0u);
    if (std::find(hyperNucNames.begin(), hyperNucNames.end(), daughter) == hyperNucNames.end())
      return -1;
    return std::find(hyperNucNames.begin(), hyperNucNames.end(), daughter) - hyperNucNames.begin();
  }
  //----------------------------------------------------------------------------------------------------------------
  std::vector<int> getDaughterVec(unsigned int hypNuc, LabeledArray<std::string> cfg)
  {
    std::vector<int> vec;
    for (unsigned int i = 0; i < 4; i++) {
      std::string daughter = cfg.get(hypNuc, i);
      if (std::find(particleNames.begin(), particleNames.end(), daughter) == particleNames.end())
        break;
      vec.push_back(std::find(particleNames.begin(), particleNames.end(), daughter) - particleNames.begin());
    }
    return vec;
  }
  //----------------------------------------------------------------------------------------------------------------

  std::vector<int> getDaughterSignVec(unsigned int hypNuc, LabeledArray<std::string> cfg)
  {
    std::vector<int> vec;
    for (unsigned int i = 0; i < 4; i++) {
      std::string sign = cfg.get(hypNuc, i);
      if (sign != "+" && sign != "-")
        break;
      vec.push_back(sign == "+" ? +1 : -1);
    }
    return vec;
  }
  //----------------------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------------------
};
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypKfRecoTask>(cfgc)};
}
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
