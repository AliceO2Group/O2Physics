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

#include "MetadataHelper.h"

#include "PWGLF/DataModel/LFHypernucleiKfTables.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTPCBase.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "MathUtils/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include "TRandom3.h"

#include <limits>
#include <map>
#include <string>
#include <vector>

// KFParticle
#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=name/macro (Name is defined in KFParticle package)
#endif
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsFull = soa::Join<aod::Collisions, aod::PIDMults, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, o2::aod::TracksDCA, aod::pidTOFmass>;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper
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
  {1120020040},
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
  int pdgCode, charge;
  double mass, resolution;
  std::array<float, 5> betheParams;
  bool active;
  DaughterParticle(std::string name_, int pdgCode_, double mass_, int charge_, LabeledArray<double> bethe) : name(name_), pdgCode(pdgCode_), charge(charge_), mass(mass_), active(false)
  {
    resolution = bethe.get(name, "resolution");
    for (unsigned int i = 0; i < betheParams.size(); i++)
      betheParams[i] = bethe.get(name, i);
  }
  int getCentralPIDIndex() { return getPIDIndex(pdgCode); }
}; // struct DaughterParticle

struct HyperNucleus {
  TString name;
  int pdgCode;
  bool active, savePrimary;
  std::vector<int> daughters, daughterTrackSigns;
  HyperNucleus(std::string name_, int pdgCode_, bool active_, std::vector<int> daughters_, std::vector<int> daughterTrackSigns_) : pdgCode(pdgCode_), active(active_), savePrimary(active_)
  {
    init(name_, daughters_, daughterTrackSigns_);
  }
  HyperNucleus(std::string name_, int pdgCode_, bool active_, int hypDaughter, std::vector<int> daughters_, std::vector<int> daughterTrackSigns_) : pdgCode(pdgCode_), active(active_), savePrimary(active_)
  {
    daughters.push_back(hypDaughter);
    init(name_, daughters_, daughterTrackSigns_);
  }
  void init(std::string name_, std::vector<int> daughters_, std::vector<int> daughterTrackSigns_)
  {
    name = TString(name_);
    for (const int& d : daughters_)
      daughters.push_back(d);
    for (const int& dc : daughterTrackSigns_)
      daughterTrackSigns.push_back(dc);
  }
  int getNdaughters() { return static_cast<int>(daughters.size()); }
  const char* motherName() { return name.Contains("->") ? ((TString)name(0, name.First("-"))).Data() : name.Data(); }
  const char* daughterNames() { return name.Contains("->") ? ((TString)name(name.First("-") + 2, name.Length())).Data() : ""; }
}; // struct HyperNucleus

struct DaughterKf {
  int64_t daughterTrackId;
  int species, hypNucId;
  KFParticle daughterKfp;
  float dcaToPv, dcaToPvXY, dcaToPvZ, tpcNsigma, tpcNsigmaNLP, tpcNsigmaNHP;
  bool active;
  std::vector<float> vtx;
  DaughterKf(int species_, int64_t daughterTrackId_, std::vector<float> vtx_, float tpcNsigma_, float tpcNsigmaNLP_, float tpcNsigmaNHP_) : daughterTrackId(daughterTrackId_), species(species_), hypNucId(-1), tpcNsigma(tpcNsigma_), tpcNsigmaNLP(tpcNsigmaNLP_), tpcNsigmaNHP(tpcNsigmaNHP_), vtx(vtx_) {}
  DaughterKf(int species_, KFParticle daughterKfp_, int hypNucId_) : daughterTrackId(-999), species(species_), hypNucId(hypNucId_), daughterKfp(daughterKfp_), dcaToPv(-999), dcaToPvXY(-999), dcaToPvZ(-999) {}
  void addKfp(KFParticle daughterKfp_)
  {
    daughterKfp = daughterKfp_;
    dcaToPvXY = daughterKfp.GetDistanceFromVertexXY(&vtx[0]);
    dcaToPv = daughterKfp.GetDistanceFromVertex(&vtx[0]);
    dcaToPvZ = std::sqrt(dcaToPv * dcaToPv - dcaToPvXY * dcaToPvXY);
  }

  bool isTrack() { return daughterTrackId >= 0; }
}; // struct DaughterKf

struct HyperNucCandidate {
  int species;
  KFParticle kfp;
  HyperNucCandidate* hypNucDaughter;
  std::vector<DaughterKf*> daughters;
  std::vector<float> recoSV;
  std::vector<std::vector<float>> daughterPosMoms;
  float mass, px, py, pz;
  float devToPvXY, dcaToPvXY, dcaToPvZ, dcaToVtxXY, dcaToVtxZ, chi2;
  bool mcTrue, isPhysPrimary, isPrimaryCandidate, isSecondaryCandidate, isUsedSecondary;
  int64_t mcParticleId;
  int tableId;
  HyperNucCandidate(int species_, HyperNucCandidate* hypNucDaughter_, std::vector<DaughterKf*> daughters_) : species(species_), hypNucDaughter(hypNucDaughter_), devToPvXY(999), dcaToPvXY(999), dcaToPvZ(999), dcaToVtxXY(999), dcaToVtxZ(999), chi2(999), mcTrue(false), isPhysPrimary(false), isPrimaryCandidate(false), isSecondaryCandidate(false), isUsedSecondary(false), mcParticleId(-1), tableId(-1)
  {
    for (const auto& d : daughters_)
      daughters.push_back(d);
    kfp.SetConstructMethod(2);
    for (size_t i = 0; i < daughters.size(); i++)
      kfp.AddDaughter(daughters.at(i)->daughterKfp);
    kfp.TransportToDecayVertex();
    chi2 = kfp.GetChi2() / kfp.GetNDF();
    recoSV.clear();
    recoSV.push_back(kfp.GetX());
    recoSV.push_back(kfp.GetY());
    recoSV.push_back(kfp.GetZ());
    mass = kfp.GetMass();
    px = kfp.GetPx();
    py = kfp.GetPy();
    pz = kfp.GetPz();
  }
  std::vector<int64_t> daughterTrackIds()
  {
    std::vector<int64_t> trackIds;
    for (const auto& daughter : daughters) {
      const auto& id = daughter->daughterTrackId;
      if (id >= 0)
        trackIds.push_back(id);
    }
    return trackIds;
  }
  bool checkKfp() { return mass != 0 && !std::isnan(mass); }
  int getDaughterTableId() { return hypNucDaughter ? hypNucDaughter->tableId : -1; }
  bool isCascade() { return hypNucDaughter != 0; }
  int getSign()
  {
    if (kfp.GetQ() == 0)
      return daughters.front()->daughterKfp.GetQ() / std::abs(daughters.front()->daughterKfp.GetQ());
    return kfp.GetQ() / std::abs(kfp.GetQ());
  }
  int getNdaughters() { return static_cast<int>(daughters.size()); }
  float getDcaTracks()
  {
    if (!daughterPosMoms.size())
      setDaughterPosMoms();
    float maxDca = std::numeric_limits<float>::lowest();
    for (size_t i = 0; i < daughters.size(); i++) {
      float dx = daughterPosMoms.at(i).at(0) - recoSV[0];
      float dy = daughterPosMoms.at(i).at(1) - recoSV[1];
      const float dca = std::sqrt(dx * dx + dy * dy);
      if (dca > maxDca)
        maxDca = dca;
    }
    return maxDca;
  }
  float getDcaMotherToVertex(std::vector<float> vtx) { return kfp.GetDistanceFromVertex(&vtx[0]); }
  double getCpa(std::vector<float> vtx)
  {
    return RecoDecay::cpa(std::array{vtx[0], vtx[1], vtx[2]}, std::array{recoSV[0], recoSV[1], recoSV[2]}, std::array{px, py, pz});
  }
  float getCt(std::vector<float> vtx)
  {
    float dl = 0;
    for (size_t i = 0; i < vtx.size(); i++) {
      float tmp = recoSV.at(i) - vtx.at(i);
      dl += (tmp * tmp);
    }
    return std::sqrt(dl) * mass / std::sqrt(px * px + py * py + pz * pz);
  }
  void setDaughterPosMoms()
  {
    for (size_t i = 0; i < daughters.size(); i++) {
      daughterPosMoms.push_back(getDaughterPosMom(i));
    }
  }
  std::vector<float> getDaughterPosMom(int daughter)
  {
    std::vector<float> posMom;
    auto kfpDaughter = daughters.at(daughter)->daughterKfp;
    kfpDaughter.TransportToPoint(&recoSV[0]);
    posMom.assign({kfpDaughter.GetX(), kfpDaughter.GetY(), kfpDaughter.GetZ(), kfpDaughter.GetPx(), kfpDaughter.GetPy(), kfpDaughter.GetPz()});
    return posMom;
  }
  float getDcaMotherToVtxXY(std::vector<float> vtx) { return kfp.GetDistanceFromVertexXY(&vtx[0]); }
  float getDcaMotherToVtxZ(std::vector<float> vtx)
  {
    kfp.TransportToPoint(&vtx[0]);
    return kfp.GetZ() - vtx[2];
  }
  void calcDcaToVtx(KFPVertex& vtx)
  {
    if (devToPvXY != 999) // o2-linter: disable=magic-number (To be checked)
      return;
    devToPvXY = kfp.GetDeviationFromVertexXY(vtx);
    dcaToPvXY = kfp.GetDistanceFromVertexXY(vtx);
    kfp.TransportToVertex(vtx);
    dcaToPvZ = kfp.GetZ() - vtx.GetZ();
  }
  void calcDcaToVtx(HyperNucCandidate& cand)
  {
    dcaToVtxXY = getDcaMotherToVtxXY(cand.recoSV);
    dcaToVtxZ = getDcaMotherToVtxZ(cand.recoSV);
  }
  float getSubDaughterMass(int d1, int d2)
  {
    return calcSubDaughterMass(daughters.at(d1)->daughterKfp, daughters.at(d2)->daughterKfp);
  }
  float getSubDaughterMassCascade(int d1, int d2)
  {
    if (!isCascade()) {
      LOGF(warning, "Primary hypernucleus has no hypernucleus daughter!");
      return -999;
    }
    return calcSubDaughterMass(daughters.at(d1)->daughterKfp, hypNucDaughter->daughters.at(d2)->daughterKfp);
  }
  float calcSubDaughterMass(KFParticle d1, KFParticle d2)
  {
    KFParticle subDaughter;
    subDaughter.SetConstructMethod(2);
    subDaughter.AddDaughter(d1);
    subDaughter.AddDaughter(d2);
    subDaughter.TransportToDecayVertex();
    return subDaughter.GetMass();
  }
}; // struct HyperNucCandidate

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
  bool hasIndex(int64_t a)
  {
    for (const auto& pair : pairs) {
      if (pair.first == a)
        return true;
    }
    return false;
  }
}; // struct IndexPairs

struct McCollInfo {
  bool hasRecoColl;
  bool passedEvSel;
  bool hasRecoParticle;
  int tableIndex;
  McCollInfo() : hasRecoColl(false), passedEvSel(false), hasRecoParticle(false), tableIndex(-1) {}
}; // struct McCollInfo

//----------------------------------------------------------------------------------------------------------------
std::vector<std::shared_ptr<TH2>> hDeDx;
std::vector<std::shared_ptr<TH1>> hInvMass;
} // namespace

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
struct HypKfRecoTask {

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

  // TPC PID Response
  bool usePidResponse;
  o2::pid::tpc::Response* response;
  std::map<std::string, std::string> metadata;
  std::array<float, 5> betheParams;

  // CCDB
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<double> bField{"bField", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "Analysis/PID/TPC/Response", "Path to the PID response object"};

  std::vector<int> activePdgs;
  std::vector<DaughterParticle> daughterParticles;
  std::vector<std::vector<DaughterKf>> foundDaughterKfs, hypNucDaughterKfs;
  std::vector<std::vector<HyperNucCandidate>> singleHyperNucCandidates, cascadeHyperNucCandidates;
  std::vector<HyperNucleus> singleHyperNuclei, cascadeHyperNuclei;
  std::vector<float> primVtx, cents;
  std::vector<McCollInfo> mcCollInfos;
  IndexPairs trackIndices, mcPartIndices;
  KFPVertex kfPrimVtx;
  bool collHasCandidate, collHasMcTrueCandidate, collPassedEvSel, activeCascade, isMC;
  int64_t mcCollTableIndex;
  int mRunNumber, occupancy;
  float dBz;
  TRandom3 rand;
  //----------------------------------------------------------------------------------------------------------------

  void init(InitContext const&)
  {
    isMC = false;
    mRunNumber = 0;
    dBz = 0;
    rand.SetSeed(0);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    usePidResponse = false;
    for (unsigned int i = 0; i < nDaughterParticles; i++) { // create daughterparticles
      daughterParticles.push_back(DaughterParticle(particleNames.at(i), particlePdgCodes.at(i), particleMasses.at(i), particleCharge.at(i), cfgBetheBlochParams));
      if (cfgTrackPIDsettings->get(i, "useBBparams") == 2 || cfgTrackPIDsettings->get(i, "useBBparams") == 0) // o2-linter: disable=magic-number (To be checked)
        usePidResponse = true;
    }

    for (unsigned int i = 0; i < nHyperNuclei; i++) { // create hypernuclei
      auto active = cfgHyperNucsActive->get(i, 0u);
      auto pdg = cfgHyperNucPdg->get(i, 0u);
      singleHyperNuclei.push_back(HyperNucleus(hyperNucNames.at(i), pdg, active, getDaughterVec(i, cfgHyperNucDaughters), getDaughterSignVec(i, cfgHyperNucSigns)));
      if (active)
        activePdgs.push_back(pdg);
    }

    activeCascade = false;
    for (unsigned int i = 0; i < nCascades; i++) { // create cascades
      auto active = cfgCascadesActive->get(i, 0u);
      auto pdg = cfgCascadesPdg->get(i, 0u);
      auto hypDaughter = getHypDaughterVec(i, cfgCascadeHypDaughter);
      cascadeHyperNuclei.push_back(HyperNucleus(cascadeNames.at(i), pdg, active, hypDaughter, getDaughterVec(i, cfgCascadeDaughters), getDaughterSignVec(i, cfgCascadeSigns)));
      if (active) {
        activePdgs.push_back(pdg);
        if (!singleHyperNuclei.at(hypDaughter).active) {
          singleHyperNuclei.at(hypDaughter).active = true;
          activePdgs.push_back(singleHyperNuclei.at(hypDaughter).pdgCode);
        }
        activeCascade = true;
      }
    }

    // define histogram axes
    const AxisSpec axisMagField{10, -10., 10., "magnetic field"};
    const AxisSpec axisNev{3, 0., 3., "Number of events"};
    const AxisSpec axisRigidity{4000, -10., 10., "#it{p}^{TPC}/#it{z}"};
    const AxisSpec axisdEdx{2000, 0, 2000, "d#it{E}/d#it{x}"};
    const AxisSpec axisInvMass{1000, 1, 6, "inv mass"};
    const AxisSpec axisCent{100, 0, 100, "centrality"};
    const AxisSpec axisOccupancy{5000, 0, 50000, "occupancy"};
    const AxisSpec axisVtxZ{100, -10, 10, "z"};
    // create histograms
    histos.add("histMagField", "histMagField", kTH1F, {axisMagField});
    histos.add("histNev", "histNev", kTH1F, {axisNev});
    histos.add("histVtxZ", "histVtxZ", kTH1F, {axisVtxZ});
    histos.add("histCentFT0A", "histCentFT0A", kTH1F, {axisCent});
    histos.add("histCentFT0C", "histCentFT0C", kTH1F, {axisCent});
    histos.add("histCentFT0M", "histCentFT0M", kTH1F, {axisCent});
    histos.add("histEvents", "histEvents", kTH2F, {axisCent, axisOccupancy});
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

  void findDaughterParticles(aod::TrackAssoc const& tracksByColl, TracksFull const& tracks, CollisionsFull::iterator const& coll)
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
        if (getRigidity(track) < cfgTrackPIDsettings->get(i, "minRigidity") || getRigidity(track) > cfgTrackPIDsettings->get(i, "maxRigidity"))
          continue;
        float tpcNsigma = getTPCnSigma(track, coll, daughterParticles.at(i));
        if (std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma"))
          continue;
        filldedx(track, i);
        if (getMeanItsClsSize(track) < cfgTrackPIDsettings->get(i, "minITSclsSize"))
          continue;
        if (getMeanItsClsSize(track) > cfgTrackPIDsettings->get(i, "maxITSclsSize"))
          continue;
        if (cfgTrackPIDsettings->get(i, "TOFrequiredabove") >= 0 && getRigidity(track) > cfgTrackPIDsettings->get(i, "TOFrequiredabove") && (track.mass() < cfgTrackPIDsettings->get(i, "minTOFmass") || track.mass() > cfgTrackPIDsettings->get(i, "maxTOFmass")))
          continue;
        float tpcNsigmaNHP = (i == kAlpha ? -999 : getTPCnSigma(track, coll, daughterParticles.at(i + 1)));
        float tpcNsigmaNLP = (i == kPion ? 999 : getTPCnSigma(track, coll, daughterParticles.at(i - 1)));
        foundDaughterKfs.at(i).push_back(DaughterKf(i, track.globalIndex(), primVtx, tpcNsigma, tpcNsigmaNLP, tpcNsigmaNHP));
      }
    } // track loop
  }
  //----------------------------------------------------------------------------------------------------------------

  void checkMCTrueTracks(aod::McTrackLabels const& trackLabels, aod::McParticles const&, TracksFull const& tracks, CollisionsFull::iterator const& coll)
  {
    for (int i = 0; i < nDaughterParticles; i++) {
      auto& daughterVec = foundDaughterKfs.at(i);
      if (!daughterVec.size())
        continue;
      for (auto it = daughterVec.end() - 1; it >= daughterVec.begin(); it--) {
        const auto& mcLab = trackLabels.rawIteratorAt(it->daughterTrackId);
        if (!mcLab.has_mcParticle()) {
          daughterVec.erase(it);
          continue;
        }
        const auto& mcPart = mcLab.mcParticle_as<aod::McParticles>();
        if (cfgSaveOnlyMcTrue) {
          if (std::abs(mcPart.pdgCode()) != daughterParticles.at(i).pdgCode) {
            daughterVec.erase(it);
            continue;
          }
          if (!mcPart.has_mothers()) {
            daughterVec.erase(it);
            continue;
          }
          bool isDaughter = false;
          for (const auto& mother : mcPart.mothers_as<aod::McParticles>()) {
            if (std::find(activePdgs.begin(), activePdgs.end(), std::abs(mother.pdgCode())) != activePdgs.end()) {
              isDaughter = true;
            }
          }
          if (!isDaughter) {
            daughterVec.erase(it);
            continue;
          }
        }
        if (cfgTrackPIDsettings->get(i, "useBBparams") == 0) {
          const auto& trk = tracks.rawIteratorAt(it->daughterTrackId);
          const auto tpcNsigmaMC = getTPCnSigmaMC(trk, coll, daughterParticles.at(i), daughterParticles.at(i));
          if (std::abs(tpcNsigmaMC) <= cfgTrackPIDsettings->get(i, "maxTPCnSigma")) {
            it->tpcNsigma = tpcNsigmaMC;
            it->tpcNsigmaNHP = (i == kAlpha ? -999 : getTPCnSigmaMC(trk, coll, daughterParticles.at(i), daughterParticles.at(i + 1)));
            it->tpcNsigmaNLP = (i == kPion ? 999 : getTPCnSigmaMC(trk, coll, daughterParticles.at(i), daughterParticles.at(i - 1)));
          } else {
            daughterVec.erase(it);
          }
        }
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  void createKFDaughters(TracksFull const& tracks)
  {
    for (size_t daughterCount = 0; daughterCount < daughterParticles.size(); daughterCount++) {
      daughterParticles.at(daughterCount).active = false;
    }
    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei, &cascadeHyperNuclei};
    bool singleHypNucActive = false;
    for (size_t vec = 0; vec < hypNucVectors.size(); vec++) {
      for (const auto& hyperNuc : *(hypNucVectors.at(vec))) {
        if (!hyperNuc.active)
          continue;
        for (size_t i = vec; i < hyperNuc.daughters.size(); i++) {
          if (foundDaughterKfs.at(hyperNuc.daughters.at(i)).size() > 0)
            daughterParticles.at(hyperNuc.daughters.at(i)).active = true;
          else
            break;
          if (i == hyperNuc.daughters.size() - 1)
            singleHypNucActive = true;
        }
      }
      if (!singleHypNucActive)
        break;
    }

    for (size_t daughterCount = 0; daughterCount < daughterParticles.size(); daughterCount++) {
      const auto& daughterParticle = daughterParticles.at(daughterCount);
      if (!daughterParticle.active) {
        foundDaughterKfs.at(daughterCount).clear();
        continue;
      }
      const auto& daughterMass = daughterParticle.mass;
      const auto& daughterCharge = daughterParticle.charge;
      auto& daughterVec = foundDaughterKfs.at(daughterCount);
      for (auto it = daughterVec.end() - 1; it >= daughterVec.begin(); it--) {
        const auto& daughterTrack = tracks.rawIteratorAt(it->daughterTrackId);
        it->addKfp(createKFParticle(daughterTrack, daughterMass, daughterCharge));
        if (std::abs(it->dcaToPvXY) < cfgTrackPIDsettings->get(daughterCount, "minDcaToPvXY") || std::abs(it->dcaToPvZ) < cfgTrackPIDsettings->get(daughterCount, "minDcaToPvZ"))
          daughterVec.erase(it);
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
      std::vector<std::vector<DaughterKf>::iterator> it;
      int nCombinations = 1;
      for (int i = 0; i < nDaughters; i++) {
        nCombinations *= foundDaughterKfs.at(hyperNuc->daughters.at(i)).size();
        it.push_back(foundDaughterKfs.at(hyperNuc->daughters.at(i)).begin());
      }
      if (!nCombinations)
        continue;
      const float reduceFactor = cfgReduce->get(hyperNucIter, 0u);
      const float minMassPrim = cfgPreSelectionsPrimaries->get(hyperNucIter, "minMass");
      const float maxMassPrim = cfgPreSelectionsPrimaries->get(hyperNucIter, "maxMass");
      const float minCtPrim = cfgPreSelectionsPrimaries->get(hyperNucIter, "minCt");
      const float maxCtPrim = cfgPreSelectionsPrimaries->get(hyperNucIter, "maxCt");
      const float minCosPaPrim = cfgPreSelectionsPrimaries->get(hyperNucIter, "minCosPa");
      const float maxDcaTracksPrim = cfgPreSelectionsPrimaries->get(hyperNucIter, "maxDcaTracks");
      const float maxDcaMotherToPvXYPrim = cfgPreSelectionsPrimaries->get(hyperNucIter, "maxDcaMotherToPvXY");
      const float maxDcaMotherToPvZPrim = cfgPreSelectionsPrimaries->get(hyperNucIter, "maxDcaMotherToPvZ");
      const float minMassSec = cfgPreSelectionsSecondaries->get(hyperNucIter, "minMass");
      const float maxMassSec = cfgPreSelectionsSecondaries->get(hyperNucIter, "maxMass");
      const float minCtSec = cfgPreSelectionsSecondaries->get(hyperNucIter, "minCt");
      const float maxCtSec = cfgPreSelectionsSecondaries->get(hyperNucIter, "maxCt");
      const float maxDcaTracksSec = cfgPreSelectionsSecondaries->get(hyperNucIter, "maxDcaTracks");
      while (it[0] != foundDaughterKfs.at(hyperNuc->daughters.at(0)).end()) {
        // check for correct signs, avoid double usage of tracks
        bool passedChecks = true;
        int checkSign = 0;
        std::vector<int64_t> vec;
        for (int i = 0; i < nDaughters; i++) {
          const auto& daughterTrack = tracks.rawIteratorAt(it[i]->daughterTrackId);
          if (!i)
            checkSign = daughterTrack.sign();
          if (daughterTrack.sign() != checkSign * hyperNuc->daughterTrackSigns.at(i) || std::find(vec.begin(), vec.end(), it[i]->daughterTrackId) != vec.end()) {
            passedChecks = false;
            break;
          }
          vec.push_back(it[i]->daughterTrackId);
        }
        if (passedChecks && rand.Rndm() <= reduceFactor) {
          std::vector<DaughterKf*> daughters;
          for (int i = 0; i < nDaughters; i++) {
            daughters.push_back(&(*it[i]));
          }
          HyperNucCandidate candidate(hyperNucIter, static_cast<HyperNucCandidate*>(0), daughters);
          // check preselections
          if (candidate.checkKfp()) {
            if (candidate.mass <= maxMassPrim && candidate.mass >= minMassPrim && candidate.getDcaTracks() <= maxDcaTracksPrim && candidate.getCt(primVtx) <= maxCtPrim && candidate.getCt(primVtx) >= minCtPrim && candidate.getCpa(primVtx) >= minCosPaPrim) {
              candidate.calcDcaToVtx(kfPrimVtx);
              if (std::abs(candidate.dcaToPvXY) <= maxDcaMotherToPvXYPrim && std::abs(candidate.dcaToPvZ) <= maxDcaMotherToPvZPrim) {
                candidate.isPrimaryCandidate = true;
                collHasCandidate = true;
              }
            }
            if (activeCascade && candidate.mass <= maxMassSec && candidate.mass >= minMassSec && candidate.getDcaTracks() <= maxDcaTracksSec && candidate.getCt(primVtx) <= maxCtSec && candidate.getCt(primVtx) >= minCtSec) {
              candidate.calcDcaToVtx(kfPrimVtx);
              candidate.isSecondaryCandidate = true;
            }
            if ((candidate.isPrimaryCandidate && hyperNuc->savePrimary) || (candidate.isSecondaryCandidate && activeCascade))
              singleHyperNucCandidates.at(hyperNucIter).push_back(candidate);
          }
        }
        it[nDaughters - 1]++;
        for (int i = nDaughters - 1; i && it[i] == foundDaughterKfs.at(hyperNuc->daughters.at(i)).end(); i--) {
          it[i] = foundDaughterKfs.at(hyperNuc->daughters.at(i)).begin();
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
      for (int64_t i = 0; i < static_cast<int64_t>(nHypNucDaughters); i++) {
        if (singleHyperNucCandidates.at(hyperNuc->daughters.at(0)).at(i).isSecondaryCandidate) {
          auto hypNucDaughter = &(singleHyperNucCandidates.at(hyperNuc->daughters.at(0)).at(i));
          hypNucDaughterKfs.at(hyperNucIter).push_back(DaughterKf(hyperNuc->daughters.at(0), hypNucDaughter->kfp, i));
        }
      }
      int nCombinations = hypNucDaughterKfs.at(hyperNucIter).size();
      std::vector<std::vector<DaughterKf>::iterator> it;
      it.push_back(hypNucDaughterKfs.at(hyperNucIter).begin());
      for (int i = 1; i < nDaughters; i++) {
        nCombinations *= foundDaughterKfs.at(hyperNuc->daughters.at(i)).size();
        it.push_back(foundDaughterKfs.at(hyperNuc->daughters.at(i)).begin());
      }
      if (!nCombinations)
        continue;
      const float minMassCas = cfgPreSelectionsCascades->get(hyperNucIter, "minMass");
      const float maxMassCas = cfgPreSelectionsCascades->get(hyperNucIter, "maxMass");
      const float minCtCas = cfgPreSelectionsCascades->get(hyperNucIter, "minCt");
      const float maxCtCas = cfgPreSelectionsCascades->get(hyperNucIter, "maxCt");
      const float minCosPaCas = cfgPreSelectionsCascades->get(hyperNucIter, "minCosPa");
      const float maxDcaTracksCas = cfgPreSelectionsCascades->get(hyperNucIter, "maxDcaTracks");
      const float maxDcaMotherToPvXYCas = cfgPreSelectionsCascades->get(hyperNucIter, "maxDcaMotherToPvXY");
      const float maxDcaMotherToPvZCas = cfgPreSelectionsCascades->get(hyperNucIter, "maxDcaMotherToPvZ");
      const float minCtSec = cfgPreSelectionsSecondaries->get(hyperNuc->daughters.at(0), "minCt");
      const float maxCtSec = cfgPreSelectionsSecondaries->get(hyperNuc->daughters.at(0), "maxCt");
      const float minCosPaSvSec = cfgPreSelectionsSecondaries->get(hyperNuc->daughters.at(0), "minCosPaSv");
      const float maxDcaMotherToSvXYSec = cfgPreSelectionsSecondaries->get(hyperNuc->daughters.at(0), "maxDcaMotherToSvXY");
      const float maxDcaMotherToSvZSec = cfgPreSelectionsSecondaries->get(hyperNuc->daughters.at(0), "maxDcaMotherToSvZ");

      while (it[0] != hypNucDaughterKfs.at(hyperNucIter).end()) {
        // select hypernuclei daughter KFParticle
        auto hypNucDaughter = &(singleHyperNucCandidates.at(hyperNuc->daughters.at(0)).at(it[0]->hypNucId));
        // check for correct signs
        int checkSign = hypNucDaughter->getSign();
        bool passedChecks = true;
        std::vector<int64_t> vec = hypNucDaughter->daughterTrackIds();
        for (int i = 1; i < nDaughters; i++) {
          const auto& daughterTrack = tracks.rawIteratorAt(it[i]->daughterTrackId);
          if (daughterTrack.sign() != checkSign * hyperNuc->daughterTrackSigns.at(i) || std::find(vec.begin(), vec.end(), it[i]->daughterTrackId) != vec.end()) {
            passedChecks = false;
            break;
          }
          vec.push_back(it[i]->daughterTrackId);
        }
        if (passedChecks) {
          std::vector<DaughterKf*> daughters;
          daughters.push_back(&(*it[0]));
          for (int i = 1; i < nDaughters; i++) {
            daughters.push_back(&(*it[i]));
          }
          HyperNucCandidate candidate(hyperNucIter, hypNucDaughter, daughters);
          if (candidate.checkKfp()) {
            // preselections for cascade and hypernucleus daughter
            if (candidate.mass <= maxMassCas && candidate.mass >= minMassCas && candidate.getDcaTracks() <= maxDcaTracksCas && candidate.getCt(primVtx) >= minCtCas && candidate.getCt(primVtx) <= maxCtCas && hypNucDaughter->getCt(candidate.recoSV) >= minCtSec && hypNucDaughter->getCt(candidate.recoSV) <= maxCtSec && candidate.getCpa(primVtx) >= minCosPaCas && hypNucDaughter->getCpa(candidate.recoSV) >= minCosPaSvSec) {
              candidate.calcDcaToVtx(kfPrimVtx);
              if (std::abs(candidate.dcaToPvXY) <= maxDcaMotherToPvXYCas && std::abs(candidate.dcaToPvZ) <= maxDcaMotherToPvZCas) {
                hypNucDaughter->calcDcaToVtx(candidate);
                if (hypNucDaughter->dcaToVtxXY <= maxDcaMotherToSvXYSec && hypNucDaughter->dcaToVtxZ <= maxDcaMotherToSvZSec) {
                  collHasCandidate = true;
                  hypNucDaughter->isUsedSecondary = true;
                  cascadeHyperNucCandidates.at(hyperNucIter).push_back(candidate);
                }
              }
            }
          }
        }
        it[nDaughters - 1]++;
        for (int i = nDaughters - 1; i && it[i] == foundDaughterKfs.at(hyperNuc->daughters.at(i)).end(); i--) {
          it[i] = foundDaughterKfs.at(hyperNuc->daughters.at(i)).begin();
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
        for (auto& hypCand : candidateVector->at(hyperNucIter)) { // o2-linter: disable=[const-ref-in-for-loop] (Object is non const and modified in loop)
          std::vector<int64_t> motherIds;
          if (hypCand.isCascade()) {
            if (!hypCand.hypNucDaughter->mcTrue)
              continue;
            const auto& mcPart = particlesMC.rawIteratorAt(hypCand.hypNucDaughter->mcParticleId);
            if (!mcPart.has_mothers())
              continue;
            for (const auto& mother : mcPart.mothers_as<aod::McParticles>()) {
              if (mother.pdgCode() == hyperNuc->pdgCode * hypCand.getSign()) {
                motherIds.push_back(mother.globalIndex());
                break;
              }
            }
          }
          for (const auto& daughter : hypCand.daughters) {
            if (!daughter->isTrack())
              continue;
            const auto& mcLab = trackLabels.rawIteratorAt(daughter->daughterTrackId);
            if (!mcLab.has_mcParticle())
              continue;
            const auto& mcPart = mcLab.mcParticle_as<aod::McParticles>();
            if (std::abs(mcPart.pdgCode()) != daughterParticles.at(daughter->species).pdgCode)
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
          if (!mcPartIndices.hasIndex(motherIds.front()))
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
      cents.at(0), cents.at(1), cents.at(2), occupancy, mRunNumber);

    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei, &cascadeHyperNuclei};
    std::vector<std::vector<std::vector<HyperNucCandidate>>*> candidateVectors = {&singleHyperNucCandidates, &cascadeHyperNucCandidates};

    for (unsigned int vec = 0; vec < candidateVectors.size(); vec++) {
      auto candidateVector = candidateVectors.at(vec);
      for (size_t hyperNucIter = 0; hyperNucIter < hypNucVectors.at(vec)->size(); hyperNucIter++) {
        HyperNucleus* hyperNuc = &(hypNucVectors.at(vec)->at(hyperNucIter));
        if (!hyperNuc->active)
          continue;
        for (auto& hypCand : candidateVector->at(hyperNucIter)) { // o2-linter: disable=const-ref-in-for-loop (Object is non const and modified in loop)
          if (!hypCand.isPrimaryCandidate && !hypCand.isUsedSecondary && !hypCand.isCascade())
            continue;
          if (saveOnlyMcTrue && !hypCand.mcTrue && !hypCand.isCascade())
            continue;
          hInvMass[vec * nHyperNuclei + hyperNucIter]->Fill(hypCand.mass);
          std::vector<int> vecDaugtherTracks, vecAddons, vecSubDaughters;
          for (const auto& daughter : hypCand.daughters) {
            if (!daughter->isTrack())
              continue;
            const auto& daughterTrackId = daughter->daughterTrackId;
            int trackTableId;
            if (!trackIndices.getIndex(daughterTrackId, trackTableId)) {
              const auto& track = tracks.rawIteratorAt(daughterTrackId);
              outputTrackTable(
                daughter->species * track.sign(), track.pt(), track.eta(), track.phi(), daughter->dcaToPvXY, daughter->dcaToPvZ, track.tpcNClsFound(), track.tpcChi2NCl(),
                track.itsClusterSizes(), track.itsChi2NCl(), getRigidity(track), track.tpcSignal(), daughter->tpcNsigma, daughter->tpcNsigmaNHP, daughter->tpcNsigmaNLP,
                track.mass(), track.isPVContributor());
              trackTableId = outputTrackTable.lastIndex();
              trackIndices.add(daughterTrackId, trackTableId);
            }
            vecDaugtherTracks.push_back(trackTableId);
          }
          for (int i = 0; i < hypCand.getNdaughters(); i++) {
            std::vector<float>& posMom = hypCand.daughterPosMoms.at(i);
            outputDaughterAddonTable(
              posMom.at(0), posMom.at(1), posMom.at(2), posMom.at(3), posMom.at(4), posMom.at(5));
            vecAddons.push_back(outputDaughterAddonTable.lastIndex());
          }
          if (hypCand.getNdaughters() > 2) { // o2-linter: disable=magic-number (To be checked)
            for (int i = 0; i < hypCand.getNdaughters(); i++) {
              for (int j = i + 1; j < hypCand.getNdaughters(); j++) {
                outputSubDaughterTable(hypCand.getSubDaughterMass(i, j));
                vecSubDaughters.push_back(outputSubDaughterTable.lastIndex());
              }
            }
          }
          if (hypCand.isCascade()) {
            for (int i = 1; i < hypCand.getNdaughters(); i++) {
              for (int j = 0; j < hypCand.hypNucDaughter->getNdaughters(); j++) {
                outputSubDaughterTable(hypCand.getSubDaughterMassCascade(i, j));
                vecSubDaughters.push_back(outputSubDaughterTable.lastIndex());
              }
            }
          }

          hypCand.kfp.TransportToDecayVertex();
          int mcPartTableId;
          outputHypNucTable(
            mcPartIndices.getIndex(hypCand.mcParticleId, mcPartTableId) ? mcPartTableId : -1,
            outputCollisionTable.lastIndex(), vecDaugtherTracks, vecAddons, hypCand.getDaughterTableId(), vecSubDaughters,
            (vec * nHyperNuclei + hyperNucIter + 1) * hypCand.getSign(), hypCand.isPrimaryCandidate, hypCand.mass,
            hypCand.px, hypCand.py, hypCand.pz, hypCand.dcaToPvXY, hypCand.dcaToPvZ, hypCand.devToPvXY,
            hypCand.dcaToVtxXY, hypCand.dcaToVtxZ, hypCand.chi2, hypCand.recoSV.at(0), hypCand.recoSV.at(1), hypCand.recoSV.at(2));
          hypCand.tableId = outputHypNucTable.lastIndex();
        }
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------

  void processMC(CollisionsFullMC const& collisions, aod::McCollisions const& mcColls, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const& particlesMC, aod::McTrackLabels const& trackLabelsMC, aod::McCollisionLabels const& collLabels, aod::TrackAssoc const& tracksColl, CollisionsFull const& colls)
  {
    isMC = true;
    mcCollInfos.clear();
    mcCollInfos.resize(mcColls.size());
    mcPartIndices.clear();
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision())
        continue;
      if (collision.sel8() && std::abs(collision.posZ()) < 10) // o2-linter: disable=magic-number (To be checked)
        mcCollInfos.at(collision.mcCollisionId()).passedEvSel = true;
    }
    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei, &cascadeHyperNuclei};
    for (const auto& mcPart : particlesMC) {
      if (!mcCollInfos.at(mcPart.mcCollisionId()).passedEvSel)
        continue;
      for (unsigned int vec = 0; vec < hypNucVectors.size(); vec++) {
        for (size_t hyperNucIter = 0; hyperNucIter < hypNucVectors.at(vec)->size(); hyperNucIter++) {
          HyperNucleus* hyperNuc = &(hypNucVectors.at(vec)->at(hyperNucIter));
          if (!hyperNuc->active)
            continue;
          if (std::abs(mcPart.pdgCode()) != hyperNuc->pdgCode)
            continue;
          bool isDecayMode = false;
          float svx, svy, svz;
          int daughterPdg = vec ? singleHyperNuclei.at(hyperNuc->daughters.at(0)).pdgCode : daughterParticles.at(hyperNuc->daughters.at(0)).pdgCode;
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
      if (!collision.has_mcCollision() || !mcCollInfos.at(collision.mcCollisionId()).passedEvSel)
        continue;

      const uint64_t collIdx = collision.globalIndex();
      auto tracksByColl = tracksColl.sliceBy(perCollision, collIdx);
      findDaughterParticles(tracksByColl, tracks, colls.rawIteratorAt(collision.globalIndex()));
      if (cfgSaveOnlyMcTrue || usePidResponse)
        checkMCTrueTracks(trackLabelsMC, particlesMC, tracks, colls.rawIteratorAt(collision.globalIndex()));
      createKFDaughters(tracks);
      createKFHypernuclei(tracks);
      createMCinfo(trackLabelsMC, collLabels, particlesMC, mcColls);
      createKFCascades(tracks);
      createMCinfo(trackLabelsMC, collLabels, particlesMC, mcColls, true);

      if (!collHasCandidate)
        continue;
      if (cfgSaveOnlyMcTrue && !collHasMcTrueCandidate)
        continue;

      mcCollTableIndex = mcCollInfos.at(collision.mcCollisionId()).tableIndex;
      if (mcCollTableIndex < 0) {
        outputMcCollisionTable(
          mcCollInfos.at(collision.mcCollisionId()).passedEvSel,
          collision.mcCollision().posX(), collision.mcCollision().posY(), collision.mcCollision().posZ());
        mcCollTableIndex = outputMcCollisionTable.lastIndex();
        mcCollInfos.at(collision.mcCollisionId()).tableIndex = mcCollTableIndex;
      }
      fillTree(tracks, cfgSaveOnlyMcTrue);
    }
  }
  PROCESS_SWITCH(HypKfRecoTask, processMC, "MC analysis", false);
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
      findDaughterParticles(tracksByColl, tracks, collision);
      createKFDaughters(tracks);
      createKFHypernuclei(tracks);
      createKFCascades(tracks);
      if (!collHasCandidate)
        continue;
      mcCollTableIndex = -1;
      fillTree(tracks);
    }
  }
  PROCESS_SWITCH(HypKfRecoTask, processData, "data analysis", true);
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
      if (bField < -990) { // o2-linter: disable=magic-number (To be checked)
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
      if (bField < -990) { // o2-linter: disable=magic-number (To be checked)
        // Fetch magnetic field from ccdb for current collision
        dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
      } else {
        dBz = bField;
      }
    }
    mRunNumber = bc.runNumber();
    KFParticle::SetField(dBz);

    // PID response
    if (!usePidResponse)
      return;
    if (metadataInfo.isFullyDefined()) {
      metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
      LOGP(info, "Automatically setting reco pass for TPC Response to {} from AO2D", metadata["RecoPassName"]);
    } else {
      LOG(info) << "Setting reco pass for TPC response to default name";
      metadata["RecoPassName"] = "apass5";
    }
    const std::string path = pidPath.value;
    ccdb->setTimestamp(run3grpTimestamp);
    response = ccdb->getSpecific<o2::pid::tpc::Response>(path, run3grpTimestamp, metadata);
    if (!response) {
      LOGF(warning, "Unable to find TPC parametrisation for specified pass name - falling back to latest object");
      response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(path, run3grpTimestamp);
      if (!response) {
        LOGF(fatal, "Unable to find any TPC object corresponding to timestamp {}!", run3grpTimestamp);
      }
    }
    LOG(info) << "Successfully retrieved TPC PID object from CCDB for timestamp " << run3grpTimestamp << ", recoPass " << metadata["RecoPassName"];
    response->PrintAll();
    betheParams = response->GetBetheBlochParams();
  }
  //----------------------------------------------------------------------------------------------------------------
  template <typename T>
  void initCollision(const T& collision)
  {
    foundDaughterKfs.clear();
    foundDaughterKfs.resize(nDaughterParticles);
    hypNucDaughterKfs.clear();
    hypNucDaughterKfs.resize(nCascades);
    singleHyperNucCandidates.clear();
    singleHyperNucCandidates.resize(nHyperNuclei);
    cascadeHyperNucCandidates.clear();
    cascadeHyperNucCandidates.resize(nCascades);
    trackIndices.clear();
    collHasCandidate = false;
    collHasMcTrueCandidate = false;
    histos.fill(HIST("histMagField"), dBz);
    histos.fill(HIST("histNev"), 0.5);
    collPassedEvSel = collision.sel8() && std::abs(collision.posZ()) < 10; // o2-linter: disable=magic-number (To be checked)
    occupancy = collision.trackOccupancyInTimeRange();
    if (collPassedEvSel) {
      histos.fill(HIST("histNev"), 1.5);
      histos.fill(HIST("histVtxZ"), collision.posZ());
      histos.fill(HIST("histCentFT0A"), collision.centFT0A());
      histos.fill(HIST("histCentFT0C"), collision.centFT0C());
      histos.fill(HIST("histCentFT0M"), collision.centFT0M());
      histos.fill(HIST("histEvents"), collision.centFT0C(), occupancy);
    }
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
    if (track.tpcNClsFound() < 100 || track.itsNCls() < 2) // o2-linter: disable=magic-number (To be checked)
      return;
    hDeDx[2 * species + 1]->Fill(track.sign() * rigidity, track.tpcSignal());
  }
  //----------------------------------------------------------------------------------------------------------------

  template <class T, class C>
  float getTPCnSigma(T const& track, C const& coll, DaughterParticle& particle)
  {
    const float rigidity = getRigidity(track);
    if (!track.hasTPC())
      return -999;
    float mMip = 1, chargeFactor = 1;
    float* parBB;

    switch (static_cast<int>(cfgTrackPIDsettings->get(particle.name, "useBBparams"))) {
      case -1:
        return 0;
      case 0:
        return isMC ? 0 : response->GetNumberOfSigma(coll, track, particle.getCentralPIDIndex());
      case 1:
        parBB = &particle.betheParams[0];
        break;
      case 2:
        mMip = response->GetMIP();
        chargeFactor = std::pow(particle.charge, response->GetChargeFactor());
        parBB = &betheParams[0];
        break;
      default:
        return -999;
    }
    double expBethe{mMip * chargeFactor * o2::common::BetheBlochAleph(static_cast<float>(particle.charge * rigidity / particle.mass), parBB[0], parBB[1], parBB[2], parBB[3], parBB[4])};
    double expSigma{expBethe * particle.resolution};
    float sigmaTPC = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
    return sigmaTPC;
  }
  //----------------------------------------------------------------------------------------------------------------

  template <class T, class C>
  float getTPCnSigmaMC(T const& trk, C const& coll, DaughterParticle& particle1, DaughterParticle& particle2)
  {
    const float pidval1 = particle1.getCentralPIDIndex();
    const float pidval2 = particle2.getCentralPIDIndex();
    const auto expSignal = response->GetExpectedSignal(trk, pidval2);
    const auto expSigma = response->GetExpectedSigma(coll, trk, pidval2);
    const auto mcTunedTPCSignal = gRandom->Gaus(expSignal, expSigma);
    return response->GetNumberOfSigmaMCTuned(coll, trk, pidval1, mcTunedTPCSignal);
  }
  //----------------------------------------------------------------------------------------------------------------

  template <class T>
  float getMeanItsClsSize(T const& track)
  {
    int sum = 0, n = 0;
    for (int i = 0; i < 0x08; i++) {
      sum += (track.itsClusterSizes() >> (0x04 * i) & 0x0f);
      if (track.itsClusterSizes() >> (0x04 * i) & 0x0f)
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
    for (int i = 0; i < 0x03; i++) {
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
    for (unsigned int i = 0; i < 0x04; i++) {
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
    for (unsigned int i = 0; i < 0x04; i++) {
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
  metadataInfo.initMetadata(cfgc); // Parse AO2D metadata
  return WorkflowSpec{
    adaptAnalysisTask<HypKfRecoTask>(cfgc)};
}
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
