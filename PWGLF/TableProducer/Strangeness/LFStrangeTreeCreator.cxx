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

#include <vector>
#include <utility>
#include <random>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"
#include "DCAFitter/DCAFitterN.h"

#include "PWGLF/DataModel/LFSlimStrangeTables.h"

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFullIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr>;

namespace
{
void momTotXYZ(std::array<float, 3>& momA, std::array<float, 3> const& momB, std::array<float, 3> const& momC)
{
  for (int i = 0; i < 3; ++i) {
    momA[i] = momB[i] + momC[i];
  }
}
float invMass2Body(std::array<float, 3> const& momA, std::array<float, 3> const& momB, std::array<float, 3> const& momC, float const& massB, float const& massC)
{
  float p2B = momB[0] * momB[0] + momB[1] * momB[1] + momB[2] * momB[2];
  float p2C = momC[0] * momC[0] + momC[1] * momC[1] + momC[2] * momC[2];
  float eB = std::sqrt(p2B + massB * massB);
  float eC = std::sqrt(p2C + massC * massC);
  float eA = eB + eC;
  float massA = std::sqrt(eA * eA - momA[0] * momA[0] - momA[1] * momA[1] - momA[2] * momA[2]);
  return massA;
}
float alphaAP(std::array<float, 3> const& momA, std::array<float, 3> const& momB, std::array<float, 3> const& momC)
{
  float momTot = std::sqrt(std::pow(momA[0], 2.) + std::pow(momA[1], 2.) + std::pow(momA[2], 2.));
  float lQlPos = (momB[0] * momA[0] + momB[1] * momA[1] + momB[2] * momA[2]) / momTot;
  float lQlNeg = (momC[0] * momA[0] + momC[1] * momA[1] + momC[2] * momA[2]) / momTot;
  return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

float qtAP(std::array<float, 3> const& momA, std::array<float, 3> const& momB)
{
  float dp = momA[0] * momB[0] + momA[1] * momB[1] + momA[2] * momB[2];
  float p2A = momA[0] * momA[0] + momA[1] * momA[1] + momA[2] * momA[2];
  float p2B = momB[0] * momB[0] + momB[1] * momB[1] + momB[2] * momB[2];
  return std::sqrt(p2B - dp * dp / p2A);
}

float etaFromMom(std::array<float, 3> const& momA, std::array<float, 3> const& momB)
{
  if (std::sqrt((1.f * momA[0] + 1.f * momB[0]) * (1.f * momA[0] + 1.f * momB[0]) +
                (1.f * momA[1] + 1.f * momB[1]) * (1.f * momA[1] + 1.f * momB[1]) +
                (1.f * momA[2] + 1.f * momB[2]) * (1.f * momA[2] + 1.f * momB[2])) -
        (1.f * momA[2] + 1.f * momB[2]) <
      static_cast<float>(1e-7)) {
    if ((1.f * momA[2] + 1.f * momB[2]) < 0.f)
      return -100.f;
    return 100.f;
  }
  return 0.5f * std::log((std::sqrt((1.f * momA[0] + 1.f * momB[0]) * (1.f * momA[0] + 1.f * momB[0]) +
                                    (1.f * momA[1] + 1.f * momB[1]) * (1.f * momA[1] + 1.f * momB[1]) +
                                    (1.f * momA[2] + 1.f * momB[2]) * (1.f * momA[2] + 1.f * momB[2])) +
                          (1.f * momA[2] + 1.f * momB[2])) /
                         (std::sqrt((1.f * momA[0] + 1.f * momB[0]) * (1.f * momA[0] + 1.f * momB[0]) +
                                    (1.f * momA[1] + 1.f * momB[1]) * (1.f * momA[1] + 1.f * momB[1]) +
                                    (1.f * momA[2] + 1.f * momB[2]) * (1.f * momA[2] + 1.f * momB[2])) -
                          (1.f * momA[2] + 1.f * momB[2])));
}
float CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
{
  return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
}
} // namespace

struct CandidateV0 {
  float pt = -999.f;
  float eta = -999.f;
  float ct = -999.f;
  float len = -999.f;
  float mass = -999.f;
  float radius = -999.f;
  float cpa = -999.f;
  float alphaAP = -999.f;
  float qtAP = -999.f;
  uint8_t isFD = 0u;
  o2::track::TrackParCov trackv0;
  std::array<float, 3> mompos;
  std::array<float, 3> momneg;
  std::array<float, 3> momposMC;
  std::array<float, 3> momnegMC;
  float dcav0daugh = -999.f;
  float dcanegpv = -999.f;
  float dcapospv = -999.f;
  float dcav0pv = -999.f;
  float tpcnsigmaneg = -999.f;
  float tpcnsigmapos = -999.f;
  float genpt = -999.f;
  float geneta = -999.f;
  float genct = -999.f;
  float genlen = -999.f;
  int pdgcode = -999;
  int pdgcodemother = -999;
  int pdgposdau = -999;
  int pdgcodemotherdaupos = -999;
  int pdgnegdau = -999;
  int pdgcodemotherdauneg = -999;
  bool isreco = 0;
  int64_t pdgmatchmothersecondmother = -999;
  int64_t mcIndex = -999;
  int64_t globalIndex = -999;
  int64_t globalIndexPos = -999;
  int64_t globalIndexNeg = -999;
};

struct LFStrangeTreeCreator {
  Produces<o2::aod::LambdaTableML> lambdaTableML;
  Produces<o2::aod::V0TableAP> v0TableAP;
  Produces<o2::aod::McLambdaTableML> mcLambdaTableML;
  Produces<o2::aod::McV0TableAP> mcV0TableAP;
  std::mt19937 gen32;
  std::vector<CandidateV0> candidateV0s;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> fitter;

  int mRunNumber;
  float d_bz;
  // o2::base::MatLayerCylSet* lut = nullptr;
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};

  // binning of (anti)lambda mass QA histograms
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, o2::constants::physics::MassLambda0 - 0.03f, o2::constants::physics::MassLambda0 + 0.03f}, "binning for the lambda invariant-mass"};
  ConfigurableAxis massXiAxis{"massXiAxis", {400, o2::constants::physics::MassXiMinus - 0.05f, o2::constants::physics::MassXiMinus + 0.05f}, "binning for the Xi invariant-mass"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};
  Configurable<float> etaMaxV0dau{"etaMaxV0dau", 0.8f, "maximum eta V0 daughters"};
  ConfigurableAxis momAxis{"momAxisFine", {5.e2, 0.f, 5.f}, "momentum axis binning"};

  Configurable<float> downscaleFactor{"downscaleFactor", 1.f, "downscaling factor"};
  Configurable<bool> applyAdditionalEvSel{"applyAdditionalEvSel", false, "apply additional event selections"};

  Configurable<float> lambdaPtMin{"lambdaPtMin", 1.f, "minimum (anti)lambda pT (GeV/c)"};
  Configurable<float> lambdaPtMax{"lambdaPtMax", 4.f, "maximum (anti)lambda pT (GeV/c)"};

  Configurable<float> v0trackNcrossedRows{"v0trackNcrossedRows", 100, "Minimum number of crossed TPC rows for V0 daughter"};
  Configurable<float> v0trackNclusItsCut{"v0trackNclusITScut", 0, "Minimum number of ITS clusters for V0 daughter"};
  Configurable<float> v0trackNclusTpcCut{"v0trackNclusTPCcut", 100, "Minimum number of TPC clusters for V0 daughter"};
  Configurable<float> v0trackNsharedClusTpc{"v0trackNsharedClusTpc", 5, "Maximum number of shared TPC clusters for V0 daughter"};
  Configurable<float> vetoMassK0Short{"vetoMassK0Short", 0.01f, "veto for V0 compatible with K0s mass"};
  Configurable<float> v0radiusMax{"v0radiusMax", 100.f, "maximum V0 radius eccepted"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 0.5f, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcav0pv{"v0setting_dcav0pv", 1.f, "DCA V0 to Pv"};
  Configurable<float> v0setting_dcadaughtopv{"v0setting_dcadaughtopv", 0.1f, "DCA Pos To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.99f, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 5.f, "v0radius"};
  Configurable<float> v0setting_lifetime{"v0setting_lifetime", 40.f, "v0 lifetime cut"};
  Configurable<float> v0setting_nsigmatpc{"v0setting_nsigmatpc", 4.f, "nsigmatpc"};
  Configurable<float> cascsetting_dcabachpv{"cascsetting_dcabachpv", 0.1f, "cascdcabachpv"};
  Configurable<float> cascsetting_cospa{"cascsetting_cospa", 0.99f, "casc cospa cut"};
  Configurable<float> cascsetting_dcav0bach{"cascsetting_dcav0bach", 1.0f, "dcav0bach"};
  Configurable<float> cascsetting_vetoOm{"cascsetting_vetoOm", 0.01f, "vetoOm"};
  Configurable<float> cascsetting_mXi{"cascsetting_mXi", 0.02f, "mXi"};
  Configurable<float> lambdaMassCut{"lambdaMassCut", 0.02f, "maximum deviation from PDG mass (for QA histograms)"};
  Configurable<bool> k0short{"k0short", false, "process for k0short (true) or lambda (false)"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  uint32_t randomSeed = 0.;

  Preslice<TracksFullIU> perCollisionTracksFullIU = o2::aod::track::collisionId;
  Preslice<aod::V0s> perCollisionV0 = o2::aod::v0::collisionId;
  Preslice<aod::Cascades> perCollisionCasc = o2::aod::cascade::collisionId;
  Preslice<aod::McParticles> perCollisionMcParts = o2::aod::mcparticle::mcCollisionId;

  template <class T>
  bool selectV0Daughter(T const& track)
  {
    if (std::abs(track.eta()) > etaMaxV0dau) {
      return false;
    }
    if (track.itsNCls() < v0trackNclusItsCut ||
        track.tpcNClsFound() < v0trackNclusTpcCut ||
        track.tpcNClsCrossedRows() < v0trackNclusTpcCut ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcNClsShared() > v0trackNsharedClusTpc) {
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

    auto timestamp = bc.timestamp();
    o2::parameters::GRPMagField* grpmag = 0x0;

    auto grpmagPath{"GLO/Config/GRPMagField"};
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField for timestamp " << timestamp;
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);

    // Fetch magnetic field from ccdb for current collision
    d_bz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << d_bz << " kG";
    mRunNumber = bc.runNumber();
    fitter.setBz(d_bz);

    // o2::base::Propagator::Instance()->setMatLUT(lut);
  }

  void init(o2::framework::InitContext&)
  {

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    // lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(4);
    fitter.setMaxDXYIni(4);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);
    int mat{static_cast<int>(cfgMaterialCorrection)};
    fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

    // event QA
    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});

    // v0 QA
    histos.add<TH3>("QA/massLambda", ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, massLambdaAxis});
    histos.add<TH3>("QA/massXi", ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{M}(#Lambda + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, massXiAxis});
  }

  template <class C, class T>
  void fillRecoEvent(C const& collision, T const& /* tracks */, aod::V0s const& V0s, aod::V0s const& /* V0s_all */, aod::Cascades const& cascades, float const& centrality)
  {
    candidateV0s.clear();

    gpu::gpustd::array<float, 2> dcaInfo;

    for (const auto& v0 : V0s) {
      auto posTrack = v0.posTrack_as<T>();
      auto negTrack = v0.negTrack_as<T>();

      bool posSelect = selectV0Daughter(posTrack);
      bool negSelect = selectV0Daughter(negTrack);
      if (!posSelect || !negSelect)
        continue;

      auto posTrackCov = getTrackParCov(posTrack);
      auto negTrackCov = getTrackParCov(negTrack);

      int nCand = 0;
      try {
        nCand = fitter.process(posTrackCov, negTrackCov);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        continue;
      }
      if (nCand == 0) {
        continue;
      }

      auto& posPropTrack = fitter.getTrack(0);
      auto& negPropTrack = fitter.getTrack(1);

      std::array<float, 3> momPos = std::array{static_cast<float>(-999.), static_cast<float>(-999.), static_cast<float>(-999.)};
      std::array<float, 3> momNeg = std::array{static_cast<float>(-999.), static_cast<float>(-999.), static_cast<float>(-999.)};
      std::array<float, 3> momV0 = std::array{static_cast<float>(-999.), static_cast<float>(-999.), static_cast<float>(-999.)};
      posPropTrack.getPxPyPzGlo(momPos);
      negPropTrack.getPxPyPzGlo(momNeg);
      momTotXYZ(momV0, momPos, momNeg);

      auto ptV0 = std::hypot(momV0[0], momV0[1]);
      if (ptV0 < lambdaPtMin || ptV0 > lambdaPtMax) {
        continue;
      }

      auto etaV0 = etaFromMom(momPos, momNeg);
      if (std::abs(etaV0) > etaMax) {
        continue;
      }

      auto alpha = alphaAP(momV0, momPos, momNeg);
      bool matter = alpha > 0;
      auto massPos = matter ? o2::constants::physics::MassProton : o2::constants::physics::MassPionCharged;
      auto massNeg = matter ? o2::constants::physics::MassPionCharged : o2::constants::physics::MassProton;
      auto mLambda = invMass2Body(momV0, momPos, momNeg, massPos, massNeg);
      auto mK0Short = invMass2Body(momV0, momPos, momNeg, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged);

      // pid selections
      auto nSigmaTPCPos = matter ? posTrack.tpcNSigmaPr() : posTrack.tpcNSigmaPi();
      auto nSigmaTPCNeg = matter ? negTrack.tpcNSigmaPi() : negTrack.tpcNSigmaPr();
      // change for k0
      if (k0short) {
        nSigmaTPCPos = posTrack.tpcNSigmaPi();
        nSigmaTPCNeg = negTrack.tpcNSigmaPi();
      }

      if (std::abs(nSigmaTPCPos) > v0setting_nsigmatpc || std::abs(nSigmaTPCNeg) > v0setting_nsigmatpc) {
        continue;
      }

      // veto on K0s mass (only for lambda)
      if (!k0short && (std::abs(mK0Short - o2::constants::physics::MassK0Short) < vetoMassK0Short)) {
        continue;
      }

      float dcaV0dau = std::sqrt(fitter.getChi2AtPCACandidate());
      if (dcaV0dau > v0setting_dcav0dau) {
        continue;
      }

      std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
      const auto& vtx = fitter.getPCACandidate();

      float radiusV0 = std::hypot(vtx[0], vtx[1]);
      if (radiusV0 < v0setting_radius || radiusV0 > v0radiusMax) {
        continue;
      }

      float dcaV0Pv = CalculateDCAStraightToPV(
        vtx[0], vtx[1], vtx[2],
        momPos[0] + momNeg[0],
        momPos[1] + momNeg[1],
        momPos[2] + momNeg[2],
        collision.posX(), collision.posY(), collision.posZ());
      if (std::abs(dcaV0Pv) > v0setting_dcav0pv) {
        continue;
      }

      double cosPA = RecoDecay::cpa(primVtx, vtx, momV0);
      if (cosPA < v0setting_cospa) {
        continue;
      }

      auto ptotal = RecoDecay::sqrtSumOfSquares(momV0[0], momV0[1], momV0[2]);
      auto lengthTraveled = RecoDecay::sqrtSumOfSquares(vtx[0] - primVtx[0], vtx[1] - primVtx[1], vtx[2] - primVtx[2]);
      // change calculation of ML2P for k0 and lambda
      float particlemass;
      if (k0short) {
        particlemass = o2::constants::physics::MassK0;
      } else {
        particlemass = o2::constants::physics::MassLambda;
      }
      float ML2P = particlemass * lengthTraveled / ptotal;
      if (ML2P > v0setting_lifetime) {
        continue;
      }

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, posTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      auto posDcaToPv = std::hypot(dcaInfo[0], dcaInfo[1]);
      if (posDcaToPv < v0setting_dcadaughtopv && std::abs(dcaInfo[0]) < v0setting_dcadaughtopv) {
        continue;
      }

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, negTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      auto negDcaToPv = std::hypot(dcaInfo[0], dcaInfo[1]);
      if (negDcaToPv < v0setting_dcadaughtopv && std::abs(dcaInfo[0]) < v0setting_dcadaughtopv) {
        continue;
      }

      if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCut) { // for QA histograms
        continue;
      }
      histos.fill(HIST("QA/massLambda"), centrality, ptV0, mLambda);

      CandidateV0 candV0;
      candV0.pt = matter > 0. ? ptV0 : -ptV0;
      candV0.eta = etaV0;
      candV0.ct = ML2P;
      candV0.len = lengthTraveled;
      candV0.mass = mLambda;
      candV0.radius = radiusV0;
      candV0.cpa = cosPA;
      candV0.alphaAP = alpha;
      candV0.qtAP = qtAP(momV0, momPos);
      candV0.trackv0 = fitter.createParentTrackParCov();
      candV0.mompos = std::array{momPos[0], momPos[1], momPos[2]};
      candV0.momneg = std::array{momNeg[0], momNeg[1], momNeg[2]};
      candV0.dcav0daugh = dcaV0dau;
      candV0.dcav0pv = dcaV0Pv;
      candV0.dcanegpv = negDcaToPv;
      candV0.dcapospv = posDcaToPv;
      candV0.tpcnsigmaneg = nSigmaTPCNeg;
      candV0.tpcnsigmapos = nSigmaTPCPos;
      candV0.globalIndex = v0.globalIndex();
      candV0.globalIndexPos = posTrack.globalIndex();
      candV0.globalIndexNeg = negTrack.globalIndex();
      candidateV0s.push_back(candV0);
    }

    for (auto& casc : cascades) {
      auto v0 = casc.template v0_as<aod::V0s>();
      auto itv0 = find_if(candidateV0s.begin(), candidateV0s.end(), [&](CandidateV0 v0cand) { return v0cand.globalIndex == v0.globalIndex(); });
      if (itv0 == candidateV0s.end()) {
        continue;
      }
      auto bachTrack = casc.template bachelor_as<T>();
      auto bachTrackPar = getTrackPar(bachTrack);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, bachTrackPar, 2.f, fitter.getMatCorrType(), &dcaInfo);

      if (TMath::Abs(dcaInfo[0]) < cascsetting_dcabachpv)
        continue;

      auto bachelorTrack = getTrackParCov(bachTrack);
      int nCand = 0;
      try {
        nCand = fitter.process(itv0->trackv0, bachelorTrack);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        continue;
      }
      if (nCand == 0)
        continue;

      auto v0Track = fitter.getTrack(0);
      bachelorTrack = fitter.getTrack(1);

      std::array<float, 3> momV0;
      std::array<float, 3> momBach;
      std::array<float, 3> momCasc;
      v0Track.getPxPyPzGlo(momV0);
      bachelorTrack.getPxPyPzGlo(momBach);
      momTotXYZ(momCasc, momV0, momBach);

      auto dcacascv0bach = TMath::Sqrt(fitter.getChi2AtPCACandidate());
      if (dcacascv0bach > cascsetting_dcav0bach)
        continue;

      std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
      const auto& vtx = fitter.getPCACandidate();
      double cosPA = RecoDecay::cpa(primVtx, vtx, momCasc);
      if (cosPA < cascsetting_cospa)
        continue;

      float mXi = invMass2Body(momCasc, momV0, momBach, o2::constants::physics::MassLambda0, o2::constants::physics::MassPionCharged);
      float mOm = invMass2Body(momCasc, momV0, momBach, o2::constants::physics::MassLambda0, o2::constants::physics::MassKaonCharged);

      if (std::abs(mOm - o2::constants::physics::MassOmegaMinus) < cascsetting_vetoOm)
        continue;

      if (std::abs(mXi - o2::constants::physics::MassXiMinus) > cascsetting_mXi)
        continue;

      histos.fill(HIST("QA/massXi"), centrality, std::hypot(momCasc[0], momCasc[1]), mXi);
      itv0->isFD = 1u; // Xi
    }
  }

  template <class C, class T>
  void fillMcEvent(C const& collision, T const& tracks, aod::V0s const& V0s, aod::V0s const& V0s_all, aod::Cascades const& cascades, float const& centrality, aod::McParticles const&, aod::McTrackLabels const& mcLabels)
  {
    fillRecoEvent<C, T>(collision, tracks, V0s, V0s_all, cascades, centrality);

    for (auto& candidateV0 : candidateV0s) {
      candidateV0.isreco = true;
      auto mcLabPos = mcLabels.rawIteratorAt(candidateV0.globalIndexPos);
      auto mcLabNeg = mcLabels.rawIteratorAt(candidateV0.globalIndexNeg);

      if (mcLabPos.has_mcParticle() && mcLabNeg.has_mcParticle()) {
        auto mcTrackPos = mcLabPos.template mcParticle_as<aod::McParticles>();
        auto mcTrackNeg = mcLabNeg.template mcParticle_as<aod::McParticles>();
        candidateV0.pdgposdau = mcTrackPos.pdgCode();
        candidateV0.pdgnegdau = mcTrackNeg.pdgCode();
        auto pdgCodeMotherDauPos = -999;
        auto pdgCodeMotherDauNeg = -999;
        auto pdgMatchMotherSecondMother = -999;
        if (mcTrackPos.has_mothers() && mcTrackNeg.has_mothers()) {
          for (auto& negMother : mcTrackNeg.template mothers_as<aod::McParticles>()) {
            for (auto& posMother : mcTrackPos.template mothers_as<aod::McParticles>()) {
              if (posMother.globalIndex() != negMother.globalIndex()) {
                pdgCodeMotherDauPos = posMother.pdgCode();
                pdgCodeMotherDauNeg = negMother.pdgCode();
                if (negMother.pdgCode() == -211) {
                  if (negMother.has_mothers()) {
                    for (auto& negSecondMother : negMother.template mothers_as<aod::McParticles>()) {
                      if (negSecondMother.globalIndex() == posMother.globalIndex()) {
                        pdgMatchMotherSecondMother = negSecondMother.pdgCode();
                      }
                    }
                  }
                }
                if (posMother.pdgCode() == 211) {
                  if (posMother.has_mothers()) {
                    for (auto& posSecondMother : posMother.template mothers_as<aod::McParticles>()) {
                      if (posSecondMother.globalIndex() == negMother.globalIndex()) {
                        pdgMatchMotherSecondMother = posSecondMother.pdgCode();
                      }
                    }
                  }
                }
              } else {
                candidateV0.pdgcode = posMother.pdgCode();
                pdgCodeMotherDauPos = posMother.pdgCode();
                pdgCodeMotherDauNeg = negMother.pdgCode();
                // build  conditions for mother/daughter for k0short or lambda
                bool mother;
                bool daughter;
                if (k0short) {
                  // mother is k0short (310) and daughters are pions (211/-211)
                  mother = posMother.pdgCode() == 310;
                  daughter = (mcTrackPos.pdgCode() == 211 && mcTrackNeg.pdgCode() == -211);
                } else {
                  // mother is lambda (3122) and daughters are proton (2212) and pion(211)
                  mother = posMother.pdgCode() == 3122;
                  daughter = ((mcTrackPos.pdgCode() == 2212 && mcTrackNeg.pdgCode() == -211) || (mcTrackPos.pdgCode() == 211 && mcTrackNeg.pdgCode() == -2212));
                }
                // check conditions
                if (!mother || !daughter) {
                  continue;
                }

                if (!posMother.isPhysicalPrimary() && !posMother.has_mothers())
                  continue;

                auto pdgCodeMother = -999;
                if (posMother.isPhysicalPrimary()) {
                  pdgCodeMother = 0;
                } else if (posMother.has_mothers()) {
                  for (auto& mcMother : posMother.mothers_as<aod::McParticles>()) {
                    // feed-down: xi and omega decaying to lambda, ignore for k0
                    if (!k0short && (std::abs(mcMother.pdgCode()) == 3322 || std::abs(mcMother.pdgCode()) == 3312 || std::abs(mcMother.pdgCode()) == 3334)) {
                      pdgCodeMother = mcMother.pdgCode();
                      break;
                    }
                  }
                }
                auto genPt = std::hypot(posMother.px(), posMother.py());
                auto posPrimVtx = std::array{posMother.vx(), posMother.vy(), posMother.vz()};
                auto secVtx = std::array{mcTrackPos.vx(), mcTrackPos.vy(), mcTrackPos.vz()};
                auto mom = std::sqrt(std::pow(posMother.px(), 2) + std::pow(posMother.py(), 2) + std::pow(posMother.pz(), 2));
                auto len = std::sqrt(std::pow(secVtx[0] - posPrimVtx[0], 2) + std::pow(secVtx[1] - posPrimVtx[1], 2) + std::pow(secVtx[2] - posPrimVtx[2], 2));
                candidateV0.genpt = genPt;
                candidateV0.genlen = len;
                candidateV0.genct = len / (mom + 1e-10) * o2::constants::physics::MassLambda0;
                candidateV0.pdgcodemother = pdgCodeMother;
                candidateV0.geneta = posMother.eta();
                candidateV0.mcIndex = posMother.globalIndex();
              }
            }
          }
        }
        if ((!mcTrackPos.has_mothers()) && mcTrackNeg.has_mothers()) {
          pdgCodeMotherDauPos = -999;
          for (auto& negMother : mcTrackNeg.template mothers_as<aod::McParticles>()) {
            pdgCodeMotherDauNeg = negMother.pdgCode();
          }
        }
        if ((!mcTrackNeg.has_mothers()) && mcTrackPos.has_mothers()) {
          pdgCodeMotherDauNeg = -999;
          for (auto& posMother : mcTrackPos.template mothers_as<aod::McParticles>()) {
            pdgCodeMotherDauPos = posMother.pdgCode();
          }
        }
        if ((!mcTrackNeg.has_mothers()) && (!mcTrackPos.has_mothers())) {
          pdgCodeMotherDauNeg = -999;
          pdgCodeMotherDauPos = -999;
        }
        candidateV0.pdgcodemotherdauneg = pdgCodeMotherDauNeg;
        candidateV0.pdgcodemotherdaupos = pdgCodeMotherDauPos;
        candidateV0.pdgmatchmothersecondmother = pdgMatchMotherSecondMother;
        // momentum of daughters
        std::array<float, 3> momPosMC = std::array{static_cast<float>(-999.), static_cast<float>(-999.), static_cast<float>(-999.)};
        std::array<float, 3> momNegMC = std::array{static_cast<float>(-999.), static_cast<float>(-999.), static_cast<float>(-999.)};
        momPosMC[0] = mcTrackPos.px();
        momPosMC[1] = mcTrackPos.py();
        momPosMC[2] = mcTrackPos.pz();
        momNegMC[0] = mcTrackNeg.px();
        momNegMC[1] = mcTrackNeg.py();
        momNegMC[2] = mcTrackNeg.pz();
        candidateV0.momposMC = std::array{momPosMC[0], momPosMC[1], momPosMC[2]};
        candidateV0.momnegMC = std::array{momNegMC[0], momNegMC[1], momNegMC[2]};
      }
    }
  }

  void fillMcGen(aod::McParticles const& mcParticles, aod::McTrackLabels const& /*mcLab*/, uint64_t const& collisionId)
  {
    auto mcParticles_thisCollision = mcParticles.sliceBy(perCollisionMcParts, collisionId);
    for (auto& mcPart : mcParticles_thisCollision) {
      auto genEta = mcPart.eta();
      if (std::abs(genEta) > etaMax) {
        continue;
      }

      auto pdgCode = mcPart.pdgCode();
      std::array<float, 3> secVtx;
      std::array<float, 3> momPosMC = std::array{static_cast<float>(-999.), static_cast<float>(-999.), static_cast<float>(-999.)};
      std::array<float, 3> momNegMC = std::array{static_cast<float>(-999.), static_cast<float>(-999.), static_cast<float>(-999.)};

      // look for lambda (3122) or k0short (310)
      int pdg_test = 3122;
      if (k0short)
        pdg_test = 310;

      if (std::abs(pdgCode) == pdg_test) {
        if (!mcPart.isPhysicalPrimary() && !mcPart.has_mothers())
          continue;
        // check if its the right decay containing proton (2122) for lambda and charged pion (211) for k0short
        int pdg_particle;
        if (k0short) {
          pdg_particle = 211;
        } else {
          pdg_particle = 2212;
        }
        bool foundParticle = false;
        for (auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
          if (std::abs(mcDaught.pdgCode()) == pdg_particle) {
            foundParticle = true;
            secVtx = std::array{mcDaught.vx(), mcDaught.vy(), mcDaught.vz()};
            break;
          }
        }
        // momentum of daughters
        for (auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
          if (mcDaught.pdgCode() < 0) {
            momNegMC[0] = mcDaught.px();
            momNegMC[1] = mcDaught.py();
            momNegMC[2] = mcDaught.pz();
          } else {
            momPosMC[0] = mcDaught.px();
            momPosMC[1] = mcDaught.py();
            momPosMC[2] = mcDaught.pz();
          }
        }
        if (!foundParticle) {
          continue;
        }
        auto pdgCodeMother = -999;
        if (mcPart.isPhysicalPrimary()) {
          pdgCodeMother = 0;
        } else if (mcPart.has_mothers()) {
          for (auto& mcMother : mcPart.mothers_as<aod::McParticles>()) {
            // feed-down: xi and omega decaying to lambda, ignore for k0
            if (!k0short && (std::abs(mcMother.pdgCode()) == 3322 || std::abs(mcMother.pdgCode()) == 3312 || std::abs(mcMother.pdgCode()) == 3334)) {
              pdgCodeMother = mcMother.pdgCode();
              break;
            }
          }
        }
        auto genPt = std::hypot(mcPart.px(), mcPart.py());
        auto posPrimVtx = std::array{mcPart.vx(), mcPart.vy(), mcPart.vz()};
        auto mom = std::sqrt(std::pow(mcPart.px(), 2) + std::pow(mcPart.py(), 2) + std::pow(mcPart.pz(), 2));
        auto len = std::sqrt(std::pow(secVtx[0] - posPrimVtx[0], 2) + std::pow(secVtx[1] - posPrimVtx[1], 2) + std::pow(secVtx[2] - posPrimVtx[2], 2));

        CandidateV0 candV0;
        candV0.genpt = genPt;
        candV0.genlen = len;
        candV0.genct = len / (mom + 1e-10) * o2::constants::physics::MassLambda0;
        candV0.geneta = mcPart.eta();
        candV0.pdgcode = pdgCode;
        candV0.pdgcodemother = pdgCodeMother;
        candV0.momposMC = std::array{momPosMC[0], momPosMC[1], momPosMC[2]};
        candV0.momnegMC = std::array{momNegMC[0], momNegMC[1], momNegMC[2]};
        auto it = find_if(candidateV0s.begin(), candidateV0s.end(), [&](CandidateV0 v0) { return v0.mcIndex == mcPart.globalIndex(); });
        if (it == candidateV0s.end()) {
          candidateV0s.emplace_back(candV0);
        } else {
          continue;
        }
      }
    }
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults> const& collisions, TracksFullIU const& tracks, aod::V0s const& V0s, aod::Cascades const& cascades, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8())
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;

      if ((!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) && applyAdditionalEvSel)
        continue;

      if (downscaleFactor < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > downscaleFactor) {
        return;
      }

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      auto CascTable_thisCollision = cascades.sliceBy(perCollisionCasc, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);
      CascTable_thisCollision.bindExternalIndices(&tracks);
      CascTable_thisCollision.bindExternalIndices(&V0s);

      auto centrality = collision.centFT0C();
      fillRecoEvent(collision, tracks, V0Table_thisCollision, V0s, CascTable_thisCollision, centrality);

      for (auto& candidateV0 : candidateV0s) {
        lambdaTableML(
          candidateV0.pt,
          candidateV0.eta,
          candidateV0.mass,
          candidateV0.ct,
          candidateV0.radius,
          candidateV0.dcav0pv,
          candidateV0.dcapospv,
          candidateV0.dcanegpv,
          candidateV0.dcav0daugh,
          candidateV0.cpa,
          candidateV0.alphaAP,
          candidateV0.qtAP,
          candidateV0.tpcnsigmapos,
          candidateV0.tpcnsigmaneg,
          candidateV0.isFD);

        v0TableAP(
          candidateV0.eta,
          candidateV0.len,
          candidateV0.mompos[0],
          candidateV0.mompos[1],
          candidateV0.mompos[2],
          candidateV0.momneg[0],
          candidateV0.momneg[1],
          candidateV0.momneg[2],
          candidateV0.radius,
          candidateV0.dcav0pv,
          candidateV0.dcapospv,
          candidateV0.dcanegpv,
          candidateV0.dcav0daugh,
          candidateV0.cpa);
      }
    }
  }
  PROCESS_SWITCH(LFStrangeTreeCreator, processRun3, "process (Run 3)", false);

  void processMcRun3(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::McCollisions const& /*mcCollisions*/, TracksFullIU const& tracks, aod::V0s const& V0s, aod::Cascades const& cascades, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLab, aod::BCsWithTimestamps const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8())
        continue;

      if ((!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) && applyAdditionalEvSel)
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      auto centrality = collision.centFT0C();

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      auto CascTable_thisCollision = cascades.sliceBy(perCollisionCasc, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);
      CascTable_thisCollision.bindExternalIndices(&tracks);
      CascTable_thisCollision.bindExternalIndices(&V0s);

      fillMcEvent(collision, tracks, V0Table_thisCollision, V0s, CascTable_thisCollision, centrality, mcParticles, mcLab);
      fillMcGen(mcParticles, mcLab, collision.mcCollisionId());

      for (auto& candidateV0 : candidateV0s) {
        mcLambdaTableML(
          candidateV0.pt,
          candidateV0.eta,
          candidateV0.mass,
          candidateV0.ct,
          candidateV0.radius,
          candidateV0.dcav0pv,
          candidateV0.dcapospv,
          candidateV0.dcanegpv,
          candidateV0.dcav0daugh,
          candidateV0.cpa,
          candidateV0.alphaAP,
          candidateV0.qtAP,
          candidateV0.tpcnsigmapos,
          candidateV0.tpcnsigmaneg,
          candidateV0.genpt,
          candidateV0.geneta,
          candidateV0.genct,
          candidateV0.pdgposdau,
          candidateV0.pdgcodemotherdaupos,
          candidateV0.pdgnegdau,
          candidateV0.pdgcodemotherdauneg,
          candidateV0.pdgcode,
          candidateV0.pdgcodemother,
          candidateV0.isreco,
          candidateV0.pdgmatchmothersecondmother);

        mcV0TableAP(
          candidateV0.eta,
          candidateV0.len,
          candidateV0.mompos[0],
          candidateV0.mompos[1],
          candidateV0.mompos[2],
          candidateV0.momneg[0],
          candidateV0.momneg[1],
          candidateV0.momneg[2],
          candidateV0.momposMC[0],
          candidateV0.momposMC[1],
          candidateV0.momposMC[2],
          candidateV0.momnegMC[0],
          candidateV0.momnegMC[1],
          candidateV0.momnegMC[2],
          candidateV0.radius,
          candidateV0.dcav0pv,
          candidateV0.dcapospv,
          candidateV0.dcanegpv,
          candidateV0.dcav0daugh,
          candidateV0.cpa,
          candidateV0.geneta,
          candidateV0.genlen,
          candidateV0.pdgcode,
          candidateV0.isreco);
      }
    }
  }
  PROCESS_SWITCH(LFStrangeTreeCreator, processMcRun3, "process MC (Run 3)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LFStrangeTreeCreator>(cfgc)};
}
