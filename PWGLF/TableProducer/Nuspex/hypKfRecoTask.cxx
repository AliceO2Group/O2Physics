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

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTPCBase.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/BetheBlochAleph.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

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

using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTOFmass, aod::pidTPCFullPr, aod::pidTPCFullPi>;

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
enum Decays { kTwoBody = 2,
              kThreeBody = 3 };
constexpr float NoVal = -999.f;

static const std::vector<std::string> particleNames{"pion", "proton", "deuteron", "triton", "helion", "alpha"};
static const std::vector<int> particlePdgCodes{211, 2212, o2::constants::physics::kDeuteron, o2::constants::physics::kTriton, o2::constants::physics::kHelium3, o2::constants::physics::kAlpha};
static const std::vector<double> particleMasses{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron, o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3, o2::constants::physics::MassAlpha};
static const std::vector<int> particleCharge{1, 1, 1, 1, 2, 2};

const int nBetheParams = 8;
enum BBPAR { kP0,
             kP1,
             kP2,
             kP3,
             kP4,
             kResolution,
             kMip,
             kExp };
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution", "mip", "exp"};
constexpr double BetheBlochDefault[nDaughterParticles][nBetheParams]{
  {13.611469, 3.598765, -0.021138, 2.039562, 0.651040, 0.09, 1., 0.},    // pion
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09, 1., 0.},      // proton
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09, 1., 0.},      // deuteron
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09, 1., 0.},      // triton
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09, 1., 0.},  // helion
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09, 1., 0.}}; // alpha

const int nTrkSettings = 18;
enum TRACKPIDSETTINGS { kPIDmethodTPC,
                        kMinITSnCls,
                        kMinTPCnCls,
                        kMaxTPCchi2,
                        kMaxITSchi2,
                        kMinRigidity,
                        kMaxRigidity,
                        kMaxTPCnSigma,
                        kMaxITSnSigma,
                        kTOFrequiredabove,
                        kMinTOFmass,
                        kMaxTOFmass,
                        kMinDcaToPvXY,
                        kMinDcaToPvZ,
                        kMinITSmeanClsSize,
                        kMaxITSmeanClsSize,
                        kTrackCharge,
                        kUsePVcontributors };
static const std::vector<std::string> trackPIDsettingNames{"PIDmethodTPC", "minITSnCls", "minTPCnCls", "maxTPCchi2", "maxITSchi2", "minRigidity", "maxRigidity", "maxTPCnSigma", "maxITSnSigma", "TOFrequiredabove", "minTOFmass", "maxTOFmass", "minDcaToPvXY", "minDcaToPvZ", "minITSclsSize", "maxITSclsSize", "trackCharge", "usePVcontributors"};
constexpr double TrackPIDsettings[nDaughterParticles][nTrkSettings]{
  {0, 0, 60, 3.0, 5000, 0.15, 1.2, 2.5, 3.0, -1, 0, 100, 0., 0., 0., 1000, 0, 0},
  {1, 0, 70, 2.5, 5000, 0.20, 4.0, 3.0, 3.0, -1, 0, 100, 0., 0., 0., 1000, 0, 0},
  {1, 0, 70, 5.0, 5000, 0.50, 5.0, 3.0, 3.0, -1, 0, 100, 0., 0., 0., 1000, 0, 0},
  {1, 0, 70, 5.0, 5000, 0.50, 5.0, 3.0, 3.0, -1, 0, 100, 0., 0., 0., 1000, 0, 0},
  {1, 0, 75, 1.5, 5000, 0.50, 5.0, 3.0, 3.0, -1, 0, 100, 0., 0., 0., 1000, 0, 0},
  {1, 0, 70, 1.5, 5000, 0.50, 5.0, 3.0, 3.0, -1, 0, 100, 0., 0., 0., 1000, 0, 0}};

static const int nHyperNuclei = 10;
static const std::vector<std::string> hyperNucNames{"L->p+pi", "3LH->3He+pi", "3LH->d+p+pi", "4LH->4He+pi", "4LH->t+p+pi", "4LHe->3He+p+pi", "5LHe->4He+p+pi", "5LHe->3He+d+pi", "custom1", "custom2"};
static const int nHypNucDefs = 8;
enum HYPNUCDEFS { kEnabled,
                  kPdgCode,
                  kD1,
                  kD2,
                  kD3,
                  kD4,
                  kDsigns,
                  kUseV0for };
static const std::vector<std::string> hypNucDefsLb{"Enabled", "PDGCode", "d1", "d2", "d3", "d4", "daughterSigns", "useV0for"};
static const std::string hypNucDefs[nHyperNuclei][nHypNucDefs]{
  {"0", "3122", "proton", "pion", "none", "none", "+-", ""},
  {"0", "1010010030", "helion", "pion", "none", "none", "+-", ""},
  {"0", "1010010030", "deuteron", "proton", "pion", "none", "++-", ""},
  {"0", "1010010040", "alpha", "pion", "none", "none", "+-", ""},
  {"0", "1010010040", "triton", "proton", "pion", "none", "++-", ""},
  {"0", "1010020040", "helion", "proton", "pion", "none", "++-", ""},
  {"0", "1010020050", "alpha", "proton", "pion", "none", "++-", ""},
  {"0", "1010020050", "helion", "deuteron", "pion", "none", "++-", ""},
  {"0", "0", "none", "none", "none", "none", "", ""},
  {"0", "0", "none", "none", "none", "none", "", ""}}; // NOLINT: runtime/string

const int nSelPrim = 8;
enum PRESELECTIONSPRIMARIES { kMinMass,
                              kMaxMass,
                              kMinCt,
                              kMaxCt,
                              kMinCosPa,
                              kMaxDcaTracks,
                              kMaxDcaMotherToPvXY,
                              kMaxDcaMotherToPvZ };
static const std::vector<std::string> preSelectionPrimNames{"minMass", "maxMass", "minCt", "maxCt", "minCosPa", "maxDcaTracks", "maxDcaMotherToPvXY", "maxDcaMotherToPvZ"};
constexpr double PreSelectionsPrimaries[nHyperNuclei][nSelPrim]{
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
//----------------------------------------------------------------------------------------------------------------

struct DaughterParticle {
  TString name;
  int pdgCode, charge;
  double mass;
  std::array<float, nBetheParams> betheParams;
  std::array<float, nTrkSettings> trkSettings;
  bool active;
  DaughterParticle(std::string name_, int pdgCode_, double mass_, int charge_, LabeledArray<double> bethe, LabeledArray<double> settings) : name(name_), pdgCode(pdgCode_), charge(charge_), mass(mass_), active(false)
  {
    for (unsigned int i = 0; i < betheParams.size(); i++) {
      betheParams[i] = bethe.get(name, i);
    }
    for (unsigned int i = 0; i < trkSettings.size(); i++) {
      trkSettings[i] = settings.get(name, i);
    }
  }
  float getTPCnSigmaBB(float rigidity, float tpcSignal)
  {
    float expBethe = betheParams[kMip] * std::pow(charge, betheParams[kExp]) * o2::common::BetheBlochAleph(static_cast<float>(charge * rigidity / mass), betheParams[kP0], betheParams[kP1], betheParams[kP2], betheParams[kP3], betheParams[kP4]);
    return (tpcSignal - expBethe) / (expBethe * betheParams[kResolution]);
  }
}; // struct DaughterParticle

struct HyperNucleus {
  TString name;
  int pdgCode;
  bool active, savePrimary;
  std::vector<int> daughters, daughterTrackSigns, v0DaughterVec;
  std::vector<float> primSettings;
  HyperNucleus(std::string name_, int pdgCode_, bool active_, std::vector<int> daughters_, std::vector<int> daughterTrackSigns_, std::vector<int> v0DaughterVec_, LabeledArray<double> primSettings_) : pdgCode(pdgCode_), active(active_), savePrimary(active_)
  {
    init(name_, daughters_, daughterTrackSigns_, v0DaughterVec_);
    for (unsigned int i = 0; i < nSelPrim; i++) {
      primSettings.push_back(primSettings_.get(name, i));
    }
  }
  HyperNucleus(std::string name_, int pdgCode_, bool active_, int hypDaughter, std::vector<int> daughters_, std::vector<int> daughterTrackSigns_) : pdgCode(pdgCode_), active(active_), savePrimary(active_)
  {
    daughters.push_back(hypDaughter);
    init(name_, daughters_, daughterTrackSigns_);
  }
  void init(std::string name_, std::vector<int> daughters_, std::vector<int> daughterTrackSigns_, std::vector<int> v0DaughterVec_ = {})
  {
    name = TString(name_);
    for (const int& d : daughters_)
      daughters.push_back(d);
    for (const int& dc : daughterTrackSigns_)
      daughterTrackSigns.push_back(dc);
    for (const int& dv0 : v0DaughterVec_)
      v0DaughterVec.push_back(dv0 - 1);
  }
  int getNdaughters() { return static_cast<int>(daughters.size()); }
  std::vector<int> getV0daughters() { return v0DaughterVec; };
  std::vector<int> getNonV0daughters()
  {
    std::vector<int> vec;
    for (int daughter = 0; daughter < getNdaughters(); daughter++) {
      if (std::find(v0DaughterVec.begin(), v0DaughterVec.end(), daughter) == v0DaughterVec.end())
        vec.push_back(daughter);
    }
    return vec;
  };
  const char* motherName() { return name.Contains("->") ? ((TString)name(0, name.First("-"))).Data() : name.Data(); }
  const char* daughterNames() { return name.Contains("->") ? ((TString)name(name.First("-") + 2, name.Length())).Data() : ""; }
}; // struct HyperNucleus

struct DaughterKf {
  static int uniqueId;
  int64_t daughterTrackId;
  int id, species, sign, hypNucId;
  KFParticle daughterKfp;
  float dcaToPv, dcaToPvXY, dcaToPvZ, tpcNsigma, tpcNsigmaNLP, tpcNsigmaNHP;
  bool active;
  std::vector<float> vtx;
  DaughterKf(int species_, int64_t daughterTrackId_, int sign_, std::vector<float> vtx_, float tpcNsigma_, float tpcNsigmaNLP_, float tpcNsigmaNHP_) : daughterTrackId(daughterTrackId_), id(uniqueId++), species(species_), sign(sign_), hypNucId(-1), tpcNsigma(tpcNsigma_), tpcNsigmaNLP(tpcNsigmaNLP_), tpcNsigmaNHP(tpcNsigmaNHP_), vtx(vtx_) {}
  void addKfp(KFParticle daughterKfp_)
  {
    daughterKfp = daughterKfp_;
    dcaToPvXY = daughterKfp.GetDistanceFromVertexXY(&vtx[0]);
    dcaToPv = daughterKfp.GetDistanceFromVertex(&vtx[0]);
    dcaToPvZ = std::sqrt(dcaToPv * dcaToPv - dcaToPvXY * dcaToPvXY);
  }

  bool isTrack() { return daughterTrackId >= 0; }
};
int DaughterKf::uniqueId = 0;
// struct DaughterKf

struct HyperNucCandidate {
  int species;
  KFParticle kfp;
  HyperNucCandidate* hypNucDaughter;
  std::vector<DaughterKf*> daughters;
  std::vector<float> recoSV;
  std::vector<std::vector<float>> daughterPosMoms;
  float mass, px, py, pz;
  float devToPvXY, dcaToPvXY, dcaToPvZ, dcaToVtxXY, dcaToVtxZ, chi2, itsMeanClsSize;
  bool mcTrue, isPhysPrimary, isPrimaryCandidate, isSecondaryCandidate, isUsedSecondary;
  int64_t mcParticleId;
  int tableId;
  HyperNucCandidate(int species_, HyperNucCandidate* hypNucDaughter_, std::vector<DaughterKf*> daughters_) : species(species_), hypNucDaughter(hypNucDaughter_), devToPvXY(-999), dcaToPvXY(-999), dcaToPvZ(-999), dcaToVtxXY(-999), dcaToVtxZ(-999), chi2(-999), itsMeanClsSize(-1), mcTrue(false), isPhysPrimary(false), isPrimaryCandidate(false), isSecondaryCandidate(false), isUsedSecondary(false), mcParticleId(-1), tableId(-1)
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
  int getSign() // K0s !!!!!
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
    if (devToPvXY != NoVal)
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
      return NoVal;
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

struct DaughterCombinations {
  std::vector<std::vector<DaughterKf>::iterator> it, itBegin, itEnd;
  int nVecs, nCombinations;
  bool end;
  std::vector<int> nonV0daughters;
  DaughterCombinations(std::vector<std::vector<DaughterKf>*>& vecs, std::vector<int> nonV0daughters_) : nVecs(0), nCombinations(1), end(false), nonV0daughters(nonV0daughters_)
  {
    for (const auto& vec : vecs) {
      nVecs++;
      nCombinations *= vec->size();
      it.push_back(vec->begin());
      itEnd.push_back(vec->end());
    }
    itBegin = it;
  }
  void getNextCombination(std::vector<DaughterKf*>& vec)
  {
    int counter = 0;
    for (const auto& i : it) {
      vec.at(nonV0daughters.at(counter++)) = &(*i);
    }
    it[nVecs - 1]++;
    for (int i = nVecs - 1; i && it[i] == itEnd[i]; i--) {
      it[i] = itBegin[i];
      it[i - 1]++;
    }
    if (it[0] == itEnd[0])
      end = true;
  }
  void init() { it = itBegin; }
  bool isEmpty() { return nCombinations == 0; }
  int getNcombinations() { return nCombinations; }
}; // struct DaughterCombinations

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
  Preslice<aod::V0s> perCollisionV0 = o2::aod::v0::collisionId;
  Preslice<aod::Decay3Bodys> perCollision3b = o2::aod::decay3body::collisionId;
  PresliceUnsorted<aod::TrackedV0s> perV0 = aod::strangenesstracking::v0Id;
  PresliceUnsorted<aod::Tracked3Bodys> perDec3 = aod::strangenesstracking::decay3BodyId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> cfgSaveOnlyMcTrue{"cfgSaveOnlyMcTrue", true, "save only MCtrue candidates"};
  Configurable<int> cfgDebug{"cfgDebug", 1, "debug level"};
  Configurable<bool> cfgRigidityCorrection{"cfgRigidityCorrection", false, "apply rigidity correction"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
  Configurable<LabeledArray<std::string>> cfgHypNucDefs{"cfgHypNucDefs", {hypNucDefs[0], nHyperNuclei, nHypNucDefs, hyperNucNames, hypNucDefsLb}, "Hypernuclei definition"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {BetheBlochDefault[0], nDaughterParticles, nBetheParams, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings{"cfgTrackPIDsettings", {TrackPIDsettings[0], nDaughterParticles, nTrkSettings, particleNames, trackPIDsettingNames}, "track selection and PID criteria"};
  Configurable<LabeledArray<double>> cfgPreSelectionsPrimaries{"cfgPreSelectionsPrimaries", {PreSelectionsPrimaries[0], nHyperNuclei, nSelPrim, hyperNucNames, preSelectionPrimNames}, "selection criteria for primary hypernuclei"};
  Configurable<float> cfgVtxCutZ{"cfgVtxCutZ", 10.0f, "Accepted z-vertex range"};

  std::map<std::string, std::string> metadata;

  o2::aod::ITSResponse itsResponse;
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
  std::vector<std::vector<DaughterKf>> foundDaughterKfs;
  std::vector<std::vector<HyperNucCandidate>> singleHyperNucCandidates;
  std::vector<HyperNucleus> singleHyperNuclei, cascadeHyperNuclei;
  std::vector<float> primVtx, cents;
  std::vector<McCollInfo> mcCollInfos;
  IndexPairs trackIndices, mcPartIndices;
  KFPVertex kfPrimVtx;
  bool collHasCandidate, collHasMcTrueCandidate, collPassedEvSel, activeCascade, isMC;
  int64_t mcCollTableIndex;
  int mRunNumber, occupancy;
  float dBz;
  //----------------------------------------------------------------------------------------------------------------

  void init(o2::framework::InitContext& context)
  {
    isMC = doprocessMC ? true : false;
    mRunNumber = 0;
    dBz = 0;
    o2::aod::ITSResponse::setParameters(context, isMC);
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    for (unsigned int i = 0; i < nDaughterParticles; i++) { // create daughterparticles
      daughterParticles.push_back(DaughterParticle(particleNames.at(i), particlePdgCodes.at(i), particleMasses.at(i), particleCharge.at(i), cfgBetheBlochParams, cfgTrackPIDsettings));
    }

    for (unsigned int i = 0; i < nHyperNuclei; i++) { // create hypernuclei
      auto active = cfgHypNucDefs->get(i, "Enabled") != "0";
      auto pdg = std::stoi(cfgHypNucDefs->get(i, "PDGCode"));
      singleHyperNuclei.push_back(HyperNucleus(hyperNucNames.at(i), pdg, active, getDaughterVec(i, cfgHypNucDefs), getDaughterSignVec(i, cfgHypNucDefs), getV0DaughterVec(i, cfgHypNucDefs), cfgPreSelectionsPrimaries));
      if (active)
        activePdgs.push_back(pdg);
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
    hInvMass.resize(nHyperNuclei);
    int histCount = 0;
    std::vector<std::vector<HyperNucleus>> hypNucVectors = {singleHyperNuclei};
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

  void findDaughterParticles(aod::TrackAssoc const& tracksByColl, TracksFull const& tracks, CollisionsFull::iterator const&)
  {
    // track loop, store daughter candidates in std::vector
    for (const auto& trackId : tracksByColl) {
      const auto& track = tracks.rawIteratorAt(trackId.trackId());
      filldedx(track, nDaughterParticles);
      for (size_t i = 0; i < daughterParticles.size(); i++) {
        if (!checkTrack(track, daughterParticles.at(i)))
          continue;
        const float tpcNsigma = getTPCnSigma(track, daughterParticles.at(i));
        if (daughterParticles.at(i).trkSettings[kMaxTPCnSigma] >= 0 && std::abs(tpcNsigma) > daughterParticles.at(i).trkSettings[kMaxTPCnSigma])
          continue;
        const float itsNsigma = getITSnSigma(track, daughterParticles.at(i));
        if (daughterParticles.at(i).trkSettings[kMaxITSnSigma] >= 0 && std::abs(itsNsigma) > daughterParticles.at(i).trkSettings[kMaxITSnSigma])
          continue;
        filldedx(track, i);
        foundDaughterKfs.at(i).push_back(DaughterKf(i, track.globalIndex(), track.sign(), primVtx, 0, 0, 0));
      }
    } // track loop
  }
  //----------------------------------------------------------------------------------------------------------------

  void checkMCTrueTracks(aod::McTrackLabels const& trackLabels, aod::McParticles const&)
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
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  void createKFDaughters(TracksFull const& tracks)
  {
    for (size_t daughterCount = 0; daughterCount < daughterParticles.size(); daughterCount++) {
      daughterParticles.at(daughterCount).active = false;
    }
    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei};
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

  void createKFHypernuclei(aod::V0s const& V0s, aod::Decay3Bodys const& decay3Bodys, aod::TrackedV0s const& trackedV0s, aod::Tracked3Bodys const& tracked3Bodys)
  {
    // loop over all hypernuclei that are to be reconstructed
    for (size_t hyperNucIter = 0; hyperNucIter < singleHyperNuclei.size(); hyperNucIter++) {
      HyperNucleus* hyperNuc = &(singleHyperNuclei.at(hyperNucIter));
      if (!hyperNuc->active)
        continue;
      int nDaughters = hyperNuc->getNdaughters();
      auto nonV0daughters = hyperNuc->getNonV0daughters();
      std::vector<std::vector<DaughterKf>*> nonV0daughterKfs;
      for (const auto& d : nonV0daughters) {
        nonV0daughterKfs.push_back(&foundDaughterKfs.at(hyperNuc->daughters.at(d)));
      }
      bool hasNonV0daughters = nonV0daughters.size() > 0;
      DaughterCombinations nonV0daughterComb(nonV0daughterKfs, nonV0daughters);
      if (hasNonV0daughters && nonV0daughterComb.isEmpty())
        continue;

      auto v0daughters = hyperNuc->getV0daughters();
      bool useDecay3Body = (v0daughters.size() == Decays::kThreeBody), useV0 = (v0daughters.size() == Decays::kTwoBody);
      bool hasV0daughters = useDecay3Body || useV0;
      if ((useDecay3Body && !decay3Bodys.size()) || (useV0 && !V0s.size()))
        continue;
      auto v0Size = !hasV0daughters ? 0 : (useV0 ? V0s.size() : decay3Bodys.size());
      int v0Count = 0;

      do {
        std::vector<DaughterKf*> daughters;
        daughters.resize(nDaughters);
        std::vector<int> trackIds;
        if (useDecay3Body) {
          const auto& v0 = decay3Bodys.rawIteratorAt(v0Count);
          trackIds = std::vector{v0.track0Id(), v0.track1Id(), v0.track2Id()};
        }
        if (useV0) {
          const auto& v0 = V0s.rawIteratorAt(v0Count);
          trackIds = std::vector{v0.posTrackId(), v0.negTrackId()};
        }
        if (hasV0daughters && !findDaughterKfComb(daughters, hyperNuc, v0daughters, trackIds))
          continue;
        nonV0daughterComb.init();
        do {
          if (hasNonV0daughters)
            nonV0daughterComb.getNextCombination(daughters);

          // check for correct signs, avoid double usage of tracks
          bool passedChecks = true;
          int checkSign = daughters[0]->sign;
          ;
          std::vector<int> vec; ///!!! index + sign für daughterKFs???
          for (int i = 0; i < nDaughters; i++) {
            if (daughters[i]->sign != checkSign * hyperNuc->daughterTrackSigns.at(i) || std::find(vec.begin(), vec.end(), daughters[i]->id) != vec.end()) {
              passedChecks = false;
              break;
            }
            vec.push_back(daughters[i]->id);
          }
          if (!passedChecks)
            continue;
          HyperNucCandidate candidate(hyperNucIter, static_cast<HyperNucCandidate*>(0), daughters);
          // check preselections
          if (checkPrimaryHypNuc(candidate, hyperNuc->primSettings) && hyperNuc->savePrimary) {
            collHasCandidate = true;
            candidate.isPrimaryCandidate = true;
            if (useV0) {
              auto trackedByDecay = trackedV0s.sliceBy(perV0, v0Count);
              if (trackedByDecay.size())
                candidate.itsMeanClsSize = trackedByDecay.rawIteratorAt(0).itsClsSize();
            }
            if (useDecay3Body) {
              auto trackedByDecay = tracked3Bodys.sliceBy(perDec3, v0Count);
              if (trackedByDecay.size())
                candidate.itsMeanClsSize = trackedByDecay.rawIteratorAt(0).itsClsSize();
            }
            singleHyperNucCandidates.at(hyperNucIter).push_back(candidate);
          }
        } while (hasNonV0daughters && !nonV0daughterComb.end);
      } while (hasV0daughters && ++v0Count < v0Size);
    }
  }
  //----------------------------------------------------------------------------------------------------------------

  //----------------------------------------------------------------------------------------------------------------
  void createMCinfo(aod::McTrackLabels const& trackLabels, aod::McCollisionLabels const&, aod::McParticles const&, aod::McCollisions const&, bool cascadesOnly = false)
  {
    // check for mcTrue
    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei};
    std::vector<std::vector<std::vector<HyperNucCandidate>>*> candidateVectors = {&singleHyperNucCandidates};
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

    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei};
    std::vector<std::vector<std::vector<HyperNucCandidate>>*> candidateVectors = {&singleHyperNucCandidates};

    for (unsigned int vec = 0; vec < candidateVectors.size(); vec++) {
      auto candidateVector = candidateVectors.at(vec);
      for (size_t hyperNucIter = 0; hyperNucIter < hypNucVectors.at(vec)->size(); hyperNucIter++) {
        HyperNucleus* hyperNuc = &(hypNucVectors.at(vec)->at(hyperNucIter));
        if (!hyperNuc->active)
          continue;
        for (auto& hypCand : candidateVector->at(hyperNucIter)) { // o2-linter: disable=const-ref-in-for-loop (Object is non const and modified in loop)
          if (!hypCand.isPrimaryCandidate)
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
          if (hypCand.getNdaughters() > Decays::kTwoBody) {
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
            (vec * nHyperNuclei + hyperNucIter + 1) * hypCand.getSign(), hypCand.isPrimaryCandidate, hypCand.mass,
            hypCand.px, hypCand.py, hypCand.pz, hypCand.dcaToPvXY, hypCand.dcaToPvZ, hypCand.devToPvXY,
            hypCand.dcaToVtxXY, hypCand.dcaToVtxZ, hypCand.chi2, /* hypCand.itsMeanClsSize,*/ hypCand.recoSV.at(0), hypCand.recoSV.at(1), hypCand.recoSV.at(2));
          hypCand.tableId = outputHypNucTable.lastIndex();
        }
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------

  void processMC(CollisionsFullMC const& collisions, aod::McCollisions const& mcColls, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const& particlesMC, aod::McTrackLabels const& trackLabelsMC, aod::McCollisionLabels const& collLabels, aod::TrackAssoc const& tracksColl, CollisionsFull const& colls, aod::V0s const& V0s, aod::Decay3Bodys const& decay3bodys, aod::TrackedV0s const& trackedV0s, aod::Tracked3Bodys const& tracked3Bodys)
  {
    mcCollInfos.clear();
    mcCollInfos.resize(mcColls.size());
    mcPartIndices.clear();
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision())
        continue;
      if (collision.sel8() && std::abs(collision.posZ()) < cfgVtxCutZ)
        mcCollInfos.at(collision.mcCollisionId()).passedEvSel = true;
    }
    std::vector<std::vector<HyperNucleus>*> hypNucVectors = {&singleHyperNuclei};
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
      auto v0TableThisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      auto decay3bodyThisCollision = decay3bodys.sliceBy(perCollision3b, collIdx);

      findDaughterParticles(tracksByColl, tracks, colls.rawIteratorAt(collision.globalIndex()));
      if (cfgSaveOnlyMcTrue)
        checkMCTrueTracks(trackLabelsMC, particlesMC);
      createKFDaughters(tracks);
      createKFHypernuclei(v0TableThisCollision, decay3bodyThisCollision, trackedV0s, tracked3Bodys);
      createMCinfo(trackLabelsMC, collLabels, particlesMC, mcColls);

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
  void processData(CollisionsFull const& collisions, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::TrackAssoc const& tracksColl, aod::V0s const& V0s, aod::Decay3Bodys const& decay3bodys, aod::TrackedV0s const& trackedV0s, aod::Tracked3Bodys const& tracked3Bodys)
  {

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      initCollision(collision);
      if (!collPassedEvSel)
        continue;
      const uint64_t collIdx = collision.globalIndex();
      auto tracksByColl = tracksColl.sliceBy(perCollision, collIdx);
      auto v0TableThisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      auto decay3bodyThisCollision = decay3bodys.sliceBy(perCollision3b, collIdx);
      findDaughterParticles(tracksByColl, tracks, collision);
      createKFDaughters(tracks);
      createKFHypernuclei(v0TableThisCollision, decay3bodyThisCollision, trackedV0s, tracked3Bodys);
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
      if (bField <= NoVal) {
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
      if (bField <= NoVal) {
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
    foundDaughterKfs.clear();
    foundDaughterKfs.resize(nDaughterParticles);
    singleHyperNucCandidates.clear();
    singleHyperNucCandidates.resize(nHyperNuclei);
    trackIndices.clear();
    collHasCandidate = false;
    collHasMcTrueCandidate = false;
    histos.fill(HIST("histMagField"), dBz);
    histos.fill(HIST("histNev"), 0.5);
    collPassedEvSel = collision.sel8() && std::abs(collision.posZ()) < cfgVtxCutZ;
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
  bool checkPrimaryHypNuc(HyperNucCandidate& candidate, std::vector<float>& settings)
  {
    if (!candidate.checkKfp())
      return false;
    if (candidate.mass >= settings[kMaxMass] || candidate.mass <= settings[kMinMass])
      return false;
    if (candidate.getDcaTracks() >= settings[kMaxDcaTracks])
      return false;
    if (candidate.getCt(primVtx) >= settings[kMaxCt] || candidate.getCt(primVtx) <= settings[kMinCt])
      return false;
    if (candidate.getCpa(primVtx) <= settings[kMinCosPa])
      return false;
    candidate.calcDcaToVtx(kfPrimVtx);
    if (std::abs(candidate.dcaToPvXY) >= settings[kMaxDcaMotherToPvXY] || std::abs(candidate.dcaToPvZ) >= settings[kMaxDcaMotherToPvZ])
      return false;
    return true;
  }
  //----------------------------------------------------------------------------------------------------------------

  bool findDaughterKfComb(std::vector<DaughterKf*>& kfDaughters, HyperNucleus* hyperNuc, std::vector<int>& daughterPos, std::vector<int> trackIds)
  {
    if (daughterPos.size() != trackIds.size())
      return false;
    auto vecSize = daughterPos.size();
    std::set<int> foundMatches;
    for (size_t dpos = 0; dpos < vecSize; dpos++) {
      auto daughterParticle = hyperNuc->daughters.at(daughterPos.at(dpos));
      for (size_t trackId = 0; trackId < vecSize; trackId++) {
        auto daughterKf = findDaughterKfByTrackId(daughterParticle, trackIds.at(trackId));
        if (daughterKf >= 0) {
          kfDaughters.at(daughterPos.at(dpos)) = &foundDaughterKfs.at(daughterParticle).at(daughterKf);
          foundMatches.insert(dpos);
          break;
        }
      }
      if (foundMatches.size() != dpos + 1)
        return false;
    }
    return true;
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void filldedx(T const& track, int species)
  {
    constexpr int NTpcClsMin = 100;
    constexpr int NItsClsMin = 2;
    const float rigidity = getRigidity(track);
    hDeDx[2 * species]->Fill(track.sign() * rigidity, track.tpcSignal());
    if (track.tpcNClsFound() < NTpcClsMin || track.itsNCls() < NItsClsMin)
      return;
    hDeDx[2 * species + 1]->Fill(track.sign() * rigidity, track.tpcSignal());
  }

  template <class T>
  bool checkTrack(T const& track, DaughterParticle& particle)
  {
    if (std::abs(track.eta()) > cfgCutEta)
      return false;
    if (track.sign() * particle.trkSettings[kTrackCharge] < 0)
      return false;
    if (particle.trkSettings[kUsePVcontributors] >= 0 && particle.trkSettings[kUsePVcontributors] == track.isPVContributor())
      return false;
    if (track.tpcNClsFound() < particle.trkSettings[kMinTPCnCls])
      return false;
    if (track.tpcChi2NCl() > particle.trkSettings[kMaxTPCchi2])
      return false;
    if (track.itsNCls() < particle.trkSettings[kMinITSnCls])
      return false;
    if (track.itsChi2NCl() > particle.trkSettings[kMaxITSchi2])
      return false;
    if (getRigidity(track) < particle.trkSettings[kMinRigidity] || getRigidity(track) > particle.trkSettings[kMaxRigidity])
      return false;
    if (getMeanItsClsSize(track) < particle.trkSettings[kMinITSmeanClsSize])
      return false;
    if (getMeanItsClsSize(track) > particle.trkSettings[kMaxITSmeanClsSize])
      return false;
    if (particle.trkSettings[kTOFrequiredabove] >= 0 && getRigidity(track) > particle.trkSettings[kTOFrequiredabove] && (track.mass() < particle.trkSettings[kMinTOFmass] || track.mass() > particle.trkSettings[kMaxTOFmass]))
      return false;
    return true;
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getITSnSigma(T const& track, DaughterParticle& particle)
  {
    switch (std::abs(particle.pdgCode)) {
      case 211:
        return itsResponse.nSigmaITS<o2::track::PID::Pion>(track);
      case 2212:
        return itsResponse.nSigmaITS<o2::track::PID::Proton>(track);
      case o2::constants::physics::kDeuteron:
        return itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track);
      case o2::constants::physics::kTriton:
        return itsResponse.nSigmaITS<o2::track::PID::Triton>(track);
      case o2::constants::physics::kHelium3:
        return itsResponse.nSigmaITS<o2::track::PID::Helium3>(track);
      case o2::constants::physics::kAlpha:
        return itsResponse.nSigmaITS<o2::track::PID::Alpha>(track);
      default:
        return NoVal;
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getTPCnSigma(T const& track, DaughterParticle& particle)
  {
    const float rigidity = getRigidity(track);
    switch (static_cast<int>(cfgTrackPIDsettings->get(particle.name, "PIDmethodTPC"))) {
      case -1:
        return 0;
      case 0:
        if (particle.name == "proton")
          return track.tpcNSigmaPr();
        else if (particle.name == "pion")
          return track.tpcNSigmaPi();
        else
          return NoVal;
      case 1:
        return particle.getTPCnSigmaBB(rigidity, track.tpcSignal());
      default:
        return NoVal;
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getMeanItsClsSize(T const& track)
  {
    constexpr int NLayers = 8;
    constexpr int NBitsPerLayer = 4;
    constexpr int BitMask = (1 << NBitsPerLayer) - 1;
    int sum = 0, n = 0;
    for (int i = 0; i < NLayers; i++) {
      int clsSize = (track.itsClusterSizes() >> (NBitsPerLayer * i)) & BitMask;
      sum += clsSize;
      if (clsSize) {
        n++;
      }
    }
    const float lambda = 1. / std::cosh(track.eta());
    return n > 0 ? (static_cast<float>(sum) / n) * lambda : 0.f;
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
    constexpr int Ndim = 3;
    auto trackparCov = getTrackParCov(track);
    std::array<float, Ndim> fP;
    std::array<float, Ndim> fM;
    trackparCov.getXYZGlo(fP);
    trackparCov.getPxPyPzGlo(fM);
    float fPM[2 * Ndim];
    for (int i = 0; i < Ndim; i++) {
      fPM[i] = fP[i];
      fPM[i + Ndim] = fM[i] * std::abs(charge);
    }
    constexpr int NcovPars = 21;
    std::array<float, NcovPars> fC;
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
    for (unsigned int i = kD1; i <= kD4; i++) {
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
    std::string signs = cfg.get(hypNuc, "daughterSigns");
    for (size_t i = 0; i < signs.size(); i++) {
      std::string sign(1, signs.at(i));
      if (sign != "+" && sign != "-")
        break;
      vec.push_back(sign == "+" ? +1 : -1);
    }
    return vec;
  }
  //----------------------------------------------------------------------------------------------------------------
  std::vector<int> getV0DaughterVec(unsigned int hypNuc, LabeledArray<std::string> cfg)
  {
    std::vector<int> vec;
    std::string v0ds = cfg.get(hypNuc, "useV0for");
    for (size_t i = 0; i < v0ds.size(); i++) {
      std::string v0d(1, v0ds.at(i));
      vec.push_back(std::stoi(v0d));
    }
    return vec;
  }

  //----------------------------------------------------------------------------------------------------------------
  int findDaughterKfByTrackId(int daughter, int trackId)
  {
    const auto& daughterKfVector = foundDaughterKfs.at(daughter);
    int count = 0;
    for (const auto& d : daughterKfVector) {
      if (d.daughterTrackId == trackId)
        return count;
      count++;
    }
    return -1;
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
