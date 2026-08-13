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

#include "PWGLF/DataModel/LFDoubleOmegaTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsTable = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults>;
using Collisions = CollisionsTable::iterator;
using CollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults, aod::McCollisionLabels>;
using FullCascades = aod::CascDataExt;
using FullV0s = aod::V0Datas;
using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;

struct DoubleOmegaCandidate {
  float ptCasc1 = -999.f;
  float etaCasc1 = -999.f;
  float phiCasc1 = -999.f;
  float cascDecLength1 = -999.f;
  float omegaMassCasc1 = -999.f;
  float xiMassCasc1 = -999.f;
  float cosPACasc1 = -999.f;
  float dcaBachPVCasc1 = -999.f;
  float dcaV0BachCasc1 = -999.f;
  float nSigmaKBach1 = -999.f;

  float ptLambda = -999.f;
  float etaLambda = -999.f;
  float phiLambda = -999.f;
  float lambdaDecLength = -999.f;
  float lambdaMass = -999.f;
  float cosPALambda = -999.f;
  float dcaPosPVLambda = -999.f;
  float dcaNegPVLambda = -999.f;
  float dcaLambdaDaughters = -999.f;
  float nSigmaPrLambda = -999.f;
  float nSigmaPiLambda = -999.f;

  float ptKaon = -999.f;
  float etaKaon = -999.f;
  float phiKaon = -999.f;
  int8_t chargeKaon = 0;
  float dcaXYKaon = -999.f;
  float dcaZKaon = -999.f;
  float nSigmaKKaon = -999.f;

  float doubleOmegaMass = -999.f;
};

struct DoubleOmegaMCInfo {
  int64_t motherId = -1;
  float pt = -999.f;
  float eta = -999.f;
  float phi = -999.f;
  float decayLength = -999.f;
  int pdgCode = 0;
};

struct LambdaCandidate {
  float px = 0.f;
  float py = 0.f;
  float pz = 0.f;
  float decayLength = -999.f;
  float mass = -999.f;
  float cosPA = -999.f;
  float dcaPosToPV = -999.f;
  float dcaNegToPV = -999.f;
  float dcaDaughters = -999.f;
  float nSigmaPr = -999.f;
  float nSigmaPi = -999.f;
};

struct doubleOmegaTreeCreator {
  Produces<aod::DoubleOmegaTable> doubleOmegaTable;
  Produces<aod::DoubleOmegaTableMC> doubleOmegaTableMC;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  static constexpr int kDoubleOmegaPdg = 1060020020;
  enum FindableDaughter : uint8_t {
    kOmega,
    kLambda,
    kKaon,
    kNFindableDaughters
  };
  std::vector<int64_t> reconstructedDoubleOmegaIds;

  Preslice<TracksFull> tracksPerCollision = aod::track::collisionId;
  Preslice<FullCascades> cascadesPerCollision = aod::cascade::collisionId;
  Preslice<FullV0s> v0sPerCollision = aod::v0data::collisionId;

  int mRunNumber = 0;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis massOmegaAxis{"massOmegaAxis", {400, o2::constants::physics::MassOmegaMinus - 0.05f, o2::constants::physics::MassOmegaMinus + 0.05f}, "binning for the Omega invariant-mass"};
  ConfigurableAxis massXiAxis{"massXiAxis", {400, o2::constants::physics::MassXiMinus - 0.05f, o2::constants::physics::MassXiMinus + 0.05f}, "binning for the Xi invariant-mass"};
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, o2::constants::physics::MassLambda0 - 0.05f, o2::constants::physics::MassLambda0 + 0.05f}, "binning for the Lambda invariant-mass"};
  ConfigurableAxis momAxis{"momAxisFine", {5.e2, 0.f, 5.f}, "momentum axis binning"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.9f, "maximum eta"};
  Configurable<float> cascPtMin{"cascPtMin", 1.f, "minimum (anti)cascade pT (GeV/c)"};
  Configurable<float> cascPtMax{"cascPtMax", 5.f, "maximum (anti)cascade pT (GeV/c)"};

  Configurable<float> minNCrossedRows{"minNCrossedRows", 100, "Minimum number of crossed TPC rows"};
  Configurable<float> minNITSClus{"minNITSClus", 0., "Minimum number of ITS clusters"};
  Configurable<float> minNTPCClus{"minNTPCClus", 80, "Minimum number of TPC clusters"};
  Configurable<float> maxNSharedTPCClus{"maxNSharedTPCClus", 5, "Maximum number of shared TPC clusters"};

  Configurable<double> minCascCosPA{"minCascCosPA", 0.99f, "Minimum cosine of the pointing angle of the cascade"};
  Configurable<double> minLambdaCosPA{"minLambdaCosPA", 0.995f, "Minimum cosine of the pointing angle of the Lambda"};
  Configurable<float> nSigmaTPCCut{"nSigmaTPCCut", 3.f, "Number of sigmas for the TPC PID"};
  Configurable<float> dcaBachToPV{"dcaBachToPV", 0.05f, "Minimum DCA of a cascade bachelor to the primary vertex"};
  Configurable<float> dcaKaonToPV{"dcaKaonToPV", 0.05f, "Minimum transverse DCA of the kaon to the primary vertex"};
  Configurable<float> dcaV0DauToPV{"dcaV0DauToPV", 0.05f, "Minimum DCA of Lambda daughters to the primary vertex"};
  Configurable<float> dcaV0Bach{"dcaV0Bach", 1.f, "Maximum DCA between the V0 and cascade bachelor"};
  Configurable<float> dcaLambdaDaughters{"dcaLambdaDaughters", 1.f, "Maximum DCA between Lambda daughters"};
  Configurable<float> mXiWindow{"mXiWindow", 0.02f, "Xi mass window used by the cascade compatibility mode"};
  Configurable<float> mOmegaWindow{"mOmegaWindow", 0.01f, "Omega mass window"};
  Configurable<float> mLambdaWindow{"mLambdaWindow", 0.01f, "Lambda mass window"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  template <class T>
  bool selectTrack(T const& track)
  {
    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    if (track.itsNCls() < minNITSClus ||
        track.tpcNClsFound() < minNTPCClus ||
        track.tpcNClsCrossedRows() < minNCrossedRows ||
        track.tpcNClsCrossedRows() < 0.8f * track.tpcNClsFindable() ||
        track.tpcNClsShared() > maxNSharedTPCClus) {
      return false;
    }
    return true;
  }

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    LOG(info) << "Retrieved GRP for timestamp " << bc.timestamp();
    mRunNumber = bc.runNumber();
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fDoubleOmega,fOmegaXi");
      zorro.populateHistRegistry(histos, bc.runNumber());
    }
  }

  void init(InitContext const&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    zorroSummary.setObject(zorro.getZorroSummary());

    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    histos.add<TH2>("QA/massXi", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Lambda + #pi) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massXiAxis});
    histos.add<TH2>("QA/massOmega", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Lambda + K) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massOmegaAxis});
    histos.add<TH2>("QA/massLambda", ";#it{p}_{T} (GeV/#it{c});#it{M}(p + #pi) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massLambdaAxis});
    histos.add("MC/findableDaughters", ";Findable daughter;Entries", HistType::kTH1F, {{kNFindableDaughters, -0.5f, static_cast<float>(kNFindableDaughters) - 0.5f}});
    auto findableDaughters = histos.get<TH1>(HIST("MC/findableDaughters"));
    findableDaughters->GetXaxis()->SetBinLabel(kOmega + 1, "#Omega");
    findableDaughters->GetXaxis()->SetBinLabel(kLambda + 1, "#Lambda");
    findableDaughters->GetXaxis()->SetBinLabel(kKaon + 1, "K");
  }

  template <class C, class T, class Casc>
  bool isSelectedOmega(C const& collision, T const&, Casc const& casc)
  {
    auto bachelor = casc.template bachelor_as<T>();
    auto posDau = casc.template posTrack_as<T>();
    auto negDau = casc.template negTrack_as<T>();

    if (!selectTrack(bachelor) || !selectTrack(posDau) || !selectTrack(negDau)) {
      return false;
    }
    if (casc.sign() > 0) {
      if (std::abs(posDau.tpcNSigmaPi()) > nSigmaTPCCut || std::abs(negDau.tpcNSigmaPr()) > nSigmaTPCCut) {
        return false;
      }
    } else if (casc.sign() < 0) {
      if (std::abs(negDau.tpcNSigmaPi()) > nSigmaTPCCut || std::abs(posDau.tpcNSigmaPr()) > nSigmaTPCCut) {
        return false;
      }
    } else {
      return false;
    }
    if (std::abs(bachelor.tpcNSigmaKa()) > nSigmaTPCCut ||
        std::abs(casc.dcabachtopv()) < dcaBachToPV ||
        std::abs(casc.dcacascdaughters()) > dcaV0Bach ||
        casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < minCascCosPA ||
        std::abs(casc.eta()) > etaMax ||
        std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > mOmegaWindow) {
      return false;
    }
    return true;
  }

  template <class C, class T, class V0>
  bool isSelectedLambda(C const&, T const&, V0 const& v0, int8_t charge)
  {
    auto posDau = v0.template posTrack_as<T>();
    auto negDau = v0.template negTrack_as<T>();
    if (!selectTrack(posDau) || !selectTrack(negDau) ||
        std::abs(v0.eta()) > etaMax ||
        v0.v0cosPA() < minLambdaCosPA ||
        std::abs(v0.dcapostopv()) < dcaV0DauToPV ||
        std::abs(v0.dcanegtopv()) < dcaV0DauToPV ||
        std::abs(v0.dcaV0daughters()) > dcaLambdaDaughters) {
      return false;
    }

    const bool isMatter = charge < 0;
    const float mass = isMatter ? v0.mLambda() : v0.mAntiLambda();
    const float protonNSigma = isMatter ? posDau.tpcNSigmaPr() : negDau.tpcNSigmaPr();
    const float pionNSigma = isMatter ? negDau.tpcNSigmaPi() : posDau.tpcNSigmaPi();
    return std::abs(mass - o2::constants::physics::MassLambda0) <= mLambdaWindow &&
           std::abs(protonNSigma) <= nSigmaTPCCut &&
           std::abs(pionNSigma) <= nSigmaTPCCut;
  }

  template <class Track>
  bool isSelectedKaon(Track const& track, int8_t charge)
  {
    return track.sign() == charge &&
           selectTrack(track) &&
           std::abs(track.tpcNSigmaKa()) <= nSigmaTPCCut &&
           std::abs(track.dcaXY()) >= dcaKaonToPV;
  }

  template <class C, class T, class V0>
  LambdaCandidate getLambdaCandidateFromV0(C const& collision, T const&, V0 const& v0, int8_t charge)
  {
    auto posDau = v0.template posTrack_as<T>();
    auto negDau = v0.template negTrack_as<T>();
    const bool isMatter = charge < 0;
    return {v0.px(),
            v0.py(),
            v0.pz(),
            std::hypot(v0.x() - collision.posX(), v0.y() - collision.posY(), v0.z() - collision.posZ()),
            isMatter ? v0.mLambda() : v0.mAntiLambda(),
            v0.v0cosPA(),
            v0.dcapostopv(),
            v0.dcanegtopv(),
            v0.dcaV0daughters(),
            isMatter ? posDau.tpcNSigmaPr() : negDau.tpcNSigmaPr(),
            isMatter ? negDau.tpcNSigmaPi() : posDau.tpcNSigmaPi()};
  }

  template <class C, class T, class Casc>
  LambdaCandidate getLambdaCandidateFromCascade(C const& collision, T const&, Casc const& casc)
  {
    auto posDau = casc.template posTrack_as<T>();
    auto negDau = casc.template negTrack_as<T>();
    const bool isMatter = casc.sign() < 0;
    return {casc.pxlambda(),
            casc.pylambda(),
            casc.pzlambda(),
            std::hypot(casc.xlambda() - collision.posX(), casc.ylambda() - collision.posY(), casc.zlambda() - collision.posZ()),
            casc.mLambda(),
            casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
            casc.dcapostopv(),
            casc.dcanegtopv(),
            casc.dcaV0daughters(),
            isMatter ? posDau.tpcNSigmaPr() : negDau.tpcNSigmaPr(),
            isMatter ? negDau.tpcNSigmaPi() : posDau.tpcNSigmaPi()};
  }

  static float etaFromMomentum(float px, float py, float pz)
  {
    const float pt = std::hypot(px, py);
    return pt > 0.f ? std::asinh(pz / pt) : 0.f;
  }

  template <class T, class C, class Omega, class Kaon>
  DoubleOmegaCandidate makeCandidate(C const& collision, Omega const& omega, LambdaCandidate const& lambda, Kaon const& kaon)
  {
    DoubleOmegaCandidate cand;
    auto omegaBach = omega.template bachelor_as<T>();

    cand.ptCasc1 = omega.pt();
    cand.etaCasc1 = omega.eta();
    cand.phiCasc1 = omega.phi();
    cand.cascDecLength1 = std::hypot(omega.x() - collision.posX(), omega.y() - collision.posY(), omega.z() - collision.posZ());
    cand.omegaMassCasc1 = omega.mOmega();
    cand.xiMassCasc1 = omega.mXi();
    cand.cosPACasc1 = omega.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
    cand.dcaBachPVCasc1 = omega.dcabachtopv();
    cand.dcaV0BachCasc1 = omega.dcacascdaughters();
    cand.nSigmaKBach1 = omegaBach.tpcNSigmaKa();

    cand.ptLambda = std::hypot(lambda.px, lambda.py);
    cand.etaLambda = etaFromMomentum(lambda.px, lambda.py, lambda.pz);
    cand.phiLambda = std::atan2(lambda.py, lambda.px);
    cand.lambdaDecLength = lambda.decayLength;
    cand.lambdaMass = lambda.mass;
    cand.cosPALambda = lambda.cosPA;
    cand.dcaPosPVLambda = lambda.dcaPosToPV;
    cand.dcaNegPVLambda = lambda.dcaNegToPV;
    cand.dcaLambdaDaughters = lambda.dcaDaughters;
    cand.nSigmaPrLambda = lambda.nSigmaPr;
    cand.nSigmaPiLambda = lambda.nSigmaPi;

    cand.ptKaon = kaon.pt();
    cand.etaKaon = kaon.eta();
    cand.phiKaon = kaon.phi();
    cand.chargeKaon = kaon.sign();
    cand.dcaXYKaon = kaon.dcaXY();
    cand.dcaZKaon = kaon.dcaZ();
    cand.nSigmaKKaon = kaon.tpcNSigmaKa();

    const std::array<float, 3> momentum{
      omega.px() + lambda.px + kaon.px(),
      omega.py() + lambda.py + kaon.py(),
      omega.pz() + lambda.pz + kaon.pz()};
    const float energyOmega = std::sqrt(o2::constants::physics::MassOmegaMinus * o2::constants::physics::MassOmegaMinus + omega.p() * omega.p());
    const float momentumLambda2 = lambda.px * lambda.px + lambda.py * lambda.py + lambda.pz * lambda.pz;
    const float energyLambda = std::sqrt(o2::constants::physics::MassLambda0 * o2::constants::physics::MassLambda0 + momentumLambda2);
    const float energyKaon = std::sqrt(o2::constants::physics::MassKaonCharged * o2::constants::physics::MassKaonCharged + kaon.p() * kaon.p());
    const float energy = energyOmega + energyLambda + energyKaon;
    const float mass2 = energy * energy - momentum[0] * momentum[0] - momentum[1] * momentum[1] - momentum[2] * momentum[2];
    cand.doubleOmegaMass = std::sqrt(std::max(0.f, mass2));
    return cand;
  }

  template <class Omega, class V0, class Kaon, class CascLabels, class V0Labels, class TrackLabels>
  bool getMCInfo(DoubleOmegaMCInfo& mcInfo, Omega const& omega, V0 const& v0, Kaon const& kaon,
                 CascLabels const& cascLabels, V0Labels const& v0Labels, TrackLabels const& trackLabels)
  {
    if (omega.globalIndex() >= cascLabels.size() ||
        v0.globalIndex() >= v0Labels.size() ||
        kaon.globalIndex() >= trackLabels.size()) {
      LOG(info) << "Skipping MC info retrieval for double Omega candidate with invalid indices: "
                << "omega index = " << omega.globalIndex() << ", v0 index = " << v0.globalIndex()
                << ", kaon index = " << kaon.globalIndex();
      return false;
    }

    auto omegaLabel = cascLabels.rawIteratorAt(omega.globalIndex());
    auto v0Label = v0Labels.rawIteratorAt(v0.globalIndex());
    auto kaonLabel = trackLabels.rawIteratorAt(kaon.globalIndex());
    if (!omegaLabel.has_mcParticle() || !v0Label.has_mcParticle() || !kaonLabel.has_mcParticle()) {
      return false;
    }

    auto mcOmega = omegaLabel.template mcParticle_as<aod::McParticles>();
    auto mcLambda = v0Label.template mcParticle_as<aod::McParticles>();
    auto mcKaon = kaonLabel.template mcParticle_as<aod::McParticles>();
    const bool isMatter = omega.sign() < 0;
    if (mcOmega.pdgCode() != (isMatter ? 3334 : -3334) ||
        mcLambda.pdgCode() != (isMatter ? 3122 : -3122) ||
        mcKaon.pdgCode() != (isMatter ? -321 : 321)) {
      return false;
    }
    for (const auto& omegaMother : mcOmega.template mothers_as<aod::McParticles>()) {
      if (omegaMother.pdgCode() != (isMatter ? kDoubleOmegaPdg : -kDoubleOmegaPdg)) {
        continue;
      }
      for (const auto& lambdaMother : mcLambda.template mothers_as<aod::McParticles>()) {
        if (omegaMother.globalIndex() != lambdaMother.globalIndex()) {
          continue;
        }
        for (const auto& kaonMother : mcKaon.template mothers_as<aod::McParticles>()) {
          if (omegaMother.globalIndex() != kaonMother.globalIndex()) {
            continue;
          }
          mcInfo.motherId = omegaMother.globalIndex();
          mcInfo.pt = omegaMother.pt();
          mcInfo.eta = omegaMother.eta();
          mcInfo.phi = omegaMother.phi();
          mcInfo.decayLength = std::hypot(mcOmega.vx() - omegaMother.vx(),
                                          mcOmega.vy() - omegaMother.vy(),
                                          mcOmega.vz() - omegaMother.vz());
          mcInfo.pdgCode = omegaMother.pdgCode();
          return true;
        }
      }
    }
    return false;
  }

  template <class McParticle>
  bool getGeneratedMCInfo(DoubleOmegaMCInfo& info, McParticle const& particle)
  {
    if (std::abs(particle.pdgCode()) != kDoubleOmegaPdg) {
      return false;
    }

    const bool isMatter = particle.pdgCode() > 0;
    const int expectedOmegaPdg = isMatter ? 3334 : -3334;
    const int expectedLambdaPdg = isMatter ? 3122 : -3122;
    const int expectedKaonPdg = isMatter ? -321 : 321;
    bool foundOmega = false;
    bool foundLambda = false;
    bool foundKaon = false;
    std::array<float, 3> decayVertex{0.f, 0.f, 0.f};

    for (const auto& daughter : particle.template daughters_as<aod::McParticles>()) {
      if (daughter.pdgCode() == expectedOmegaPdg) {
        foundOmega = true;
        decayVertex = {daughter.vx(), daughter.vy(), daughter.vz()};
      } else if (daughter.pdgCode() == expectedLambdaPdg) {
        foundLambda = true;
      } else if (daughter.pdgCode() == expectedKaonPdg) {
        foundKaon = true;
      }
    }
    if (!foundOmega || !foundLambda || !foundKaon) {
      return false;
    }

    info.motherId = particle.globalIndex();
    info.pt = particle.pt();
    info.eta = particle.eta();
    info.phi = particle.phi();
    info.decayLength = std::hypot(decayVertex[0] - particle.vx(),
                                  decayVertex[1] - particle.vy(),
                                  decayVertex[2] - particle.vz());
    info.pdgCode = particle.pdgCode();
    return true;
  }

  template <class Candidates, class Labels>
  void fillFindableDaughters(Candidates const& candidates, Labels const& labels, int expectedAbsPdg, FindableDaughter daughterType)
  {
    for (const auto& candidate : candidates) {

      if (expectedAbsPdg == 3334) {
      LOG(info) << "Checking candidate with global index " << candidate.globalIndex() << " for findable daughter of type " << daughterType;
      LOG(info) << "Total number of labels: " << labels.size();
      }

      if (candidate.globalIndex() >= labels.size()) {
        continue;
      }
      auto label = labels.rawIteratorAt(candidate.globalIndex());
      if (!label.has_mcParticle()) {
        continue;
      }
      auto mcParticle = label.template mcParticle_as<aod::McParticles>();
      if (std::abs(mcParticle.pdgCode()) != expectedAbsPdg) {
        continue;
      }
      if(expectedAbsPdg == 3334){
          LOG(info) << "Found a cascade with pdg code " << mcParticle.pdgCode() << " and global index " << mcParticle.globalIndex();
          for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
              LOG(info) << "Mother pdg code: " << mother.pdgCode() << ", global index: " << mother.globalIndex();
          }
      }
      for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
        if (std::abs(mother.pdgCode()) == kDoubleOmegaPdg) {
          histos.fill(HIST("MC/findableDaughters"), daughterType);
          break;
        }
      }
    }
  }

  template <class Cascades, class V0s, class Tracks, class CascLabels, class V0Labels, class TrackLabels>
  void checkFindables(Cascades const& cascades, V0s const& v0s, Tracks const& tracks,
                      CascLabels const& cascLabels, V0Labels const& v0Labels, TrackLabels const& trackLabels)
  {
    fillFindableDaughters(cascades, cascLabels, 3334, kOmega);
    fillFindableDaughters(v0s, v0Labels, 3122, kLambda);
    fillFindableDaughters(tracks, trackLabels, 321, kKaon);
  }

  void writeDataCandidate(DoubleOmegaCandidate const& cand)
  {
    doubleOmegaTable(cand.ptCasc1,
                     cand.etaCasc1,
                     cand.phiCasc1,
                     cand.cascDecLength1,
                     cand.omegaMassCasc1,
                     cand.xiMassCasc1,
                     cand.cosPACasc1,
                     cand.dcaBachPVCasc1,
                     cand.dcaV0BachCasc1,
                     cand.nSigmaKBach1,
                     cand.ptLambda,
                     cand.etaLambda,
                     cand.phiLambda,
                     cand.lambdaDecLength,
                     cand.lambdaMass,
                     cand.cosPALambda,
                     cand.dcaPosPVLambda,
                     cand.dcaNegPVLambda,
                     cand.dcaLambdaDaughters,
                     cand.nSigmaPrLambda,
                     cand.nSigmaPiLambda,
                     cand.ptKaon,
                     cand.etaKaon,
                     cand.phiKaon,
                     cand.chargeKaon,
                     cand.dcaXYKaon,
                     cand.dcaZKaon,
                     cand.nSigmaKKaon,
                     cand.doubleOmegaMass);
  }

  void writeMCCandidate(DoubleOmegaCandidate const& cand, DoubleOmegaMCInfo const& mcInfo, bool isReco)
  {
    doubleOmegaTableMC(cand.ptCasc1,
                       cand.etaCasc1,
                       cand.phiCasc1,
                       cand.cascDecLength1,
                       cand.omegaMassCasc1,
                       cand.xiMassCasc1,
                       cand.cosPACasc1,
                       cand.dcaBachPVCasc1,
                       cand.dcaV0BachCasc1,
                       cand.nSigmaKBach1,
                       cand.ptLambda,
                       cand.etaLambda,
                       cand.phiLambda,
                       cand.lambdaDecLength,
                       cand.lambdaMass,
                       cand.cosPALambda,
                       cand.dcaPosPVLambda,
                       cand.dcaNegPVLambda,
                       cand.dcaLambdaDaughters,
                       cand.nSigmaPrLambda,
                       cand.nSigmaPiLambda,
                       cand.ptKaon,
                       cand.etaKaon,
                       cand.phiKaon,
                       cand.chargeKaon,
                       cand.dcaXYKaon,
                       cand.dcaZKaon,
                       cand.nSigmaKKaon,
                       cand.doubleOmegaMass,
                       mcInfo.pt,
                       mcInfo.eta,
                       mcInfo.phi,
                       mcInfo.decayLength,
                       mcInfo.pdgCode,
                       isReco);
  }

  template <bool isMC, class C, class T, class Cascades, class V0s, class... McLabels>
  void fillFromV0s(C const& collision, T const& tracks, Cascades const& cascades, V0s const& v0s, McLabels const&... mcLabels)
  {
    for (const auto& omega : cascades) {
      LOG(info) << "Processing omega with global index " << omega.globalIndex();
      if (!isSelectedOmega(collision, tracks, omega)) {
        continue;
      }
      histos.fill(HIST("QA/massXi"), omega.pt(), omega.mXi());
      histos.fill(HIST("QA/massOmega"), omega.pt(), omega.mOmega());
      const std::array<int64_t, 3> omegaTrackIds{omega.posTrackId(), omega.negTrackId(), omega.bachelorId()};
      for (const auto& v0 : v0s) {
        const int8_t charge = omega.sign();
        if (!isSelectedLambda(collision, tracks, v0, charge)) {
          continue;
        }
        if (std::find(omegaTrackIds.begin(), omegaTrackIds.end(), v0.posTrackId()) != omegaTrackIds.end() ||
            std::find(omegaTrackIds.begin(), omegaTrackIds.end(), v0.negTrackId()) != omegaTrackIds.end()) {
          continue;
        }
        const auto lambda = getLambdaCandidateFromV0(collision, tracks, v0, charge);
        histos.fill(HIST("QA/massLambda"), std::hypot(lambda.px, lambda.py), lambda.mass);
        for (const auto& kaon : tracks) {
          if (!isSelectedKaon(kaon, charge) ||
              kaon.globalIndex() == v0.posTrackId() ||
              kaon.globalIndex() == v0.negTrackId() ||
              std::find(omegaTrackIds.begin(), omegaTrackIds.end(), kaon.globalIndex()) != omegaTrackIds.end()) {
            continue;
          }
          auto cand = makeCandidate<T>(collision, omega, lambda, kaon);
          if constexpr (isMC) {
            static_assert(sizeof...(McLabels) == 3);
            DoubleOmegaMCInfo mcInfo;
            if (getMCInfo(mcInfo, omega, v0, kaon, mcLabels...)) {
              reconstructedDoubleOmegaIds.push_back(mcInfo.motherId);
              writeMCCandidate(cand, mcInfo, true);
            }
          } else {
            writeDataCandidate(cand);
          }
        }
      }
    }
  }

  template <class C, class T, class Cascades>
  void fillFromCascades(C const& collision, T const& tracks, Cascades const& cascades)
  {
    for (const auto& omega : cascades) {
      if (!isSelectedOmega(collision, tracks, omega)) {
        continue;
      }
      histos.fill(HIST("QA/massXi"), omega.pt(), omega.mXi());
      histos.fill(HIST("QA/massOmega"), omega.pt(), omega.mOmega());
      const std::array<int64_t, 3> omegaTrackIds{omega.posTrackId(), omega.negTrackId(), omega.bachelorId()};

      for (const auto& lambdaKaonSource : cascades) {
        if (lambdaKaonSource.globalIndex() == omega.globalIndex() ||
            lambdaKaonSource.sign() != omega.sign() ||
            !isSelectedOmega(collision, tracks, lambdaKaonSource)) {
          continue;
        }

        const std::array<int64_t, 3> sourceTrackIds{
          lambdaKaonSource.posTrackId(), lambdaKaonSource.negTrackId(), lambdaKaonSource.bachelorId()};
        bool sharesTrack = false;
        for (const auto omegaTrackId : omegaTrackIds) {
          if (std::find(sourceTrackIds.begin(), sourceTrackIds.end(), omegaTrackId) != sourceTrackIds.end()) {
            sharesTrack = true;
            break;
          }
        }
        if (sharesTrack) {
          continue;
        }

        auto kaon = lambdaKaonSource.template bachelor_as<T>();
        const auto lambda = getLambdaCandidateFromCascade(collision, tracks, lambdaKaonSource);
        histos.fill(HIST("QA/massLambda"), std::hypot(lambda.px, lambda.py), lambda.mass);
        writeDataCandidate(makeCandidate<T>(collision, omega, lambda, kaon));
      }
    }
  }

  template <class C>
  bool acceptCollision(C const& collision, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!collision.sel8() ||
        std::abs(collision.posZ()) > zVtxMax ||
        !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) ||
        !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (cfgSkimmedProcessing) {
      zorro.isSelected(bc.globalBC());
    }
    histos.fill(HIST("QA/zVtx"), collision.posZ());
    return true;
  }

  void processData(Collisions const& collision,
                   TracksFull const& tracks,
                   FullCascades const& cascades,
                   FullV0s const& v0s,
                   aod::BCsWithTimestamps const& bcs)
  {
    if (acceptCollision(collision, bcs)) {
      fillFromV0s<false>(collision, tracks, cascades, v0s);
    }
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processData, "Reconstruct Omega + Lambda(V0) + kaon(track)", true);

  void processDataFromCascades(Collisions const& collision,
                               TracksFull const& tracks,
                               FullCascades const& cascades,
                               aod::BCsWithTimestamps const& bcs)
  {
    if (acceptCollision(collision, bcs)) {
      fillFromCascades(collision, tracks, cascades);
    }
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processDataFromCascades, "Reconstruct Omega + Lambda + kaon from two cascade rows", false);

  void processMC(CollisionsMC const& collisions,
                 TracksFull const& tracks,
                 FullCascades const& cascades,
                 FullV0s const& v0s,
                 aod::McTrackLabels const& trackLabels,
                 aod::McCascLabels const& cascLabels,
                 aod::McV0Labels const& v0Labels,
                 aod::McParticles const& mcParticles,
                 aod::BCsWithTimestamps const&)
  {
    reconstructedDoubleOmegaIds.clear();

    for (const auto& collision : collisions) {
      // if (!acceptCollision(collision, bcs)) {
      //   continue;
      // }
      auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      auto cascadesThisCollision = cascades.sliceBy(cascadesPerCollision, collision.globalIndex());
      auto v0sThisCollision = v0s.sliceBy(v0sPerCollision, collision.globalIndex());
      cascadesThisCollision.bindExternalIndices(&tracks);
      v0sThisCollision.bindExternalIndices(&tracks);
      checkFindables(cascadesThisCollision, v0sThisCollision, tracksThisCollision, cascLabels, v0Labels, trackLabels);
      fillFromV0s<true>(collision, tracksThisCollision, cascadesThisCollision, v0sThisCollision, cascLabels, v0Labels, trackLabels);
    }

    for (const auto& mcParticle : mcParticles) {
      DoubleOmegaMCInfo mcInfo;
      if (!getGeneratedMCInfo(mcInfo, mcParticle)) {
        continue;
      }

      if (std::find(reconstructedDoubleOmegaIds.begin(), reconstructedDoubleOmegaIds.end(), mcInfo.motherId) != reconstructedDoubleOmegaIds.end()) {
        continue;
      }
      writeMCCandidate(DoubleOmegaCandidate{}, mcInfo, false);
    }
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processMC, "Reconstruct and truth-match Omega + Lambda + kaon in MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<doubleOmegaTreeCreator>(cfgc)};
}
