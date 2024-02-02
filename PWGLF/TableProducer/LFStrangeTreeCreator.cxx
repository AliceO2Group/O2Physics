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
#include <random>
#include <array>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGLF/DataModel/LFSlimStrangeTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPr, aod::pidTPCPi>;
using TracksFullMC = soa::Join<TracksFull, aod::McTrackLabels>;

struct CandidateV0 {
  float pt;
  float eta;
  float centFT0C;
  bool isMatter;
  float mass;
  float ct;
  float radius;
  float dcaV0PV;
  float dcaPiPV;
  float dcaPrPV;
  float dcaV0Tracks;
  float cosPA;
  float tpcNsigmaPi;
  float tpcNsigmaPr;
  float genPt;
  float genEta;
  float genCt;
  int pdgCode;
  bool isReco;
  int64_t globalIndexPos = -999;
  int64_t globalIndexNeg = -999;
};

struct LFStrangeTreeCreator {

  std::mt19937 gen32;

  Produces<aod::LambdaTableML> lambdaTableML;
  Produces<aod::McLambdaTableML> mcLambdaTableML;
  std::vector<CandidateV0> candidateV0s;
  std::vector<int64_t> checkedV0s;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<float> downscaleFactor{"downscaleFactor", 1.f, "downscaling factor"};
  Configurable<bool> fillOnlySignal{"fillOnlySignal", true, "toggle table filling of signal-only candidates (for mc)"};
  Configurable<bool> fillNonRecSignal{"fillNonRecSignal", true, "toggle table filling of non-reco signal candidates (for mc)"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, o2::constants::physics::MassLambda0 - 0.05f, o2::constants::physics::MassLambda0 + 0.05f}, "binning for the lambda invariant-mass"};
  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};

  ConfigurableAxis cosPaAxis{"cosPaAxis", {1e3, 0.95f, 1.00f}, "binning for the cosPa axis"};
  ConfigurableAxis radiusAxis{"radiusAxis", {1e3, 0.f, 100.f}, "binning for the radius axis"};
  ConfigurableAxis dcaV0daughAxis{"dcaV0daughAxis", {2e2, 0.f, 2.f}, "binning for the dca of V0 daughters"};
  ConfigurableAxis dcaDaughPvAxis{"dcaDaughPvAxis", {1e3, 0.f, 10.f}, "binning for the dca of positive daughter to PV"};

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};

  Configurable<float> lambdaPtMin{"lambdaPtMin", 0.5f, "minimum (anti)lambda pT (GeV/c)"};
  Configurable<float> lambdaPtMax{"lambdaPtMax", 3.0f, "maximum (anti)lambda pT (GeV/c)"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1f, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1f, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5f, "v0radius"};
  Configurable<float> lambdaMassCut{"lambdaMassCut", 0.05f, "maximum deviation from PDG mass"};

  int mRunNumber;
  float d_bz;

  uint32_t randomSeed = 0.;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> perCollision = aod::v0data::collisionId;
  std::vector<unsigned int> filledMothers;

  template <class Collision, class RecV0s, class T>
  void fillRecoLambda(Collision const& collision, RecV0s const& v0s, T const& tracks)
  {
    candidateV0s.clear();

    CandidateV0 candV0;
    for (const auto& v0 : v0s) {
      if (v0.pt() < lambdaPtMin || v0.pt() > lambdaPtMax) {
        continue;
      }

      if (std::abs(v0.dcapostopv()) < v0setting_dcapostopv &&
          std::abs(v0.dcanegtopv()) < v0setting_dcanegtopv &&
          std::abs(v0.dcaV0daughters()) > v0setting_dcav0dau) {
        continue;
      }

      if (std::abs(v0.eta()) > etaMax ||
          v0.v0cosPA() < v0setting_cospa ||
          v0.v0radius() < v0setting_radius) {
        return;
      }
      auto mLambda = v0.alpha() > 0 ? v0.mLambda() : v0.mAntiLambda();
      if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCut) {
        return;
      }

      auto pos = v0.template posTrack_as<T>();
      auto neg = v0.template negTrack_as<T>();
      if (std::abs(pos.eta()) > etaMax || std::abs(pos.eta()) > etaMax) {
        return;
      }

      bool matter = v0.alpha() > 0;
      histos.fill(HIST("massLambda"), matter ? v0.mLambda() : v0.mAntiLambda());
      histos.fill(HIST("cosPa"), v0.v0cosPA());
      histos.fill(HIST("radius"), v0.v0radius());
      histos.fill(HIST("dcaV0daugh"), v0.dcaV0daughters());
      histos.fill(HIST("dcaPosPv"), v0.dcapostopv());
      histos.fill(HIST("dcaNegPv"), v0.dcanegtopv());

      float ct = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      float tpcNsigmaPi = matter ? neg.tpcNSigmaPi() : pos.tpcNSigmaPi();
      float tpcNsigmaPr = matter ? pos.tpcNSigmaPr() : neg.tpcNSigmaPr();

      candV0.pt = v0.pt();
      candV0.eta = v0.eta();
      candV0.centFT0C = collision.centFT0C();
      candV0.isMatter = matter;
      candV0.mass = matter ? v0.mLambda() : v0.mAntiLambda();
      candV0.ct = ct;
      candV0.radius = v0.v0radius();
      candV0.dcaV0PV = v0.dcav0topv();
      candV0.dcaPiPV = matter ? v0.dcanegtopv() : v0.dcapostopv();
      candV0.dcaPrPV = matter ? v0.dcapostopv() : v0.dcanegtopv();
      candV0.dcaV0Tracks = v0.dcaV0daughters();
      candV0.cosPA = v0.v0cosPA();
      candV0.tpcNsigmaPi = tpcNsigmaPi;
      candV0.tpcNsigmaPr = tpcNsigmaPr;
      candV0.globalIndexPos = pos.globalIndex();
      candV0.globalIndexNeg = neg.globalIndex();

      candidateV0s.emplace_back(candV0);
      // LOGF(info, "v0.pt = %f", candV0.pt);
    }
  }

  template <class Collision, class RecV0s, class T>
  void fillMcLambda(Collision const& collision, RecV0s const& v0s, T const& tracks, aod::McParticles const&, aod::McTrackLabels const& mcLabels)
  {
    fillRecoLambda(collision, v0s, tracks);
    for (auto& candidateV0 : candidateV0s) {
      candidateV0.isReco = true;
      auto mcLabPos = mcLabels.rawIteratorAt(candidateV0.globalIndexPos);
      auto mcLabNeg = mcLabels.rawIteratorAt(candidateV0.globalIndexNeg);

      if (mcLabPos.has_mcParticle() && mcLabNeg.has_mcParticle()) {
        auto mcTrackPos = mcLabPos.template mcParticle_as<aod::McParticles>();
        auto mcTrackNeg = mcLabNeg.template mcParticle_as<aod::McParticles>();
        if (mcTrackPos.has_mothers() && mcTrackNeg.has_mothers()) {
          for (auto& negMother : mcTrackNeg.template mothers_as<aod::McParticles>()) {
            for (auto& posMother : mcTrackPos.template mothers_as<aod::McParticles>()) {
              if (posMother.globalIndex() != negMother.globalIndex())
                continue;
              if (!((mcTrackPos.pdgCode() == 2212 && mcTrackNeg.pdgCode() == -211) || (mcTrackPos.pdgCode() == 211 && mcTrackNeg.pdgCode() == -1 * 2212))) // TODO: check warning for constant comparison
                continue;
              if (std::abs(posMother.pdgCode()) != 3122)
                continue;

              auto posPrimVtx = std::array{posMother.vx(), posMother.vy(), posMother.vz()};
              auto secVtx = std::array{mcTrackPos.vx(), mcTrackPos.vy(), mcTrackPos.vz()};
              candidateV0.genPt = std::hypot(posMother.px(), posMother.py());
              auto mom = std::sqrt(std::pow(posMother.px(), 2) + std::pow(posMother.py(), 2) + std::pow(posMother.pz(), 2));
              auto len = std::sqrt(std::pow(secVtx[0] - posPrimVtx[0], 2) + std::pow(secVtx[1] - posPrimVtx[1], 2) + std::pow(secVtx[2] - posPrimVtx[2], 2));
              candidateV0.pdgCode = posMother.pdgCode();
              candidateV0.genEta = posMother.eta();
              candidateV0.genCt = len / (mom + 1e-10) * o2::constants::physics::MassLambda0;
              checkedV0s.push_back(posMother.globalIndex());
            }
          }
        }
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    mRunNumber = 0;
    d_bz = 0;

    randomSeed = static_cast<unsigned int>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});

    // v0 QA
    histos.add<TH1>("massLambda", ";#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", {HistType::kTH1F, {massLambdaAxis}});
    histos.add<TH1>("cosPa", ";cosPa;Entries", {HistType::kTH1F}, {cosPaAxis});
    histos.add<TH1>("radius", ";radius;Entries", {HistType::kTH1F}, {radiusAxis});
    histos.add<TH1>("dcaV0daugh", ";dcaV0daugh;Entries", {HistType::kTH1F}, {dcaV0daughAxis});
    histos.add<TH1>("dcaPosPv", ";dcaPosPv;Entries", {HistType::kTH1F}, {dcaDaughPvAxis});
    histos.add<TH1>("dcaNegPv", ";dcaNegPv;Entries", {HistType::kTH1F}, {dcaDaughPvAxis});
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();

    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    }
    mRunNumber = bc.runNumber();
  }

  void processData(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>::iterator const& collision, TracksFull const& tracks, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!collision.sel8())
      return;

    if (std::abs(collision.posZ()) > zVtxMax)
      return;

    histos.fill(HIST("zVtx"), collision.posZ());

    if (downscaleFactor < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > downscaleFactor) {
      return;
    }

    fillRecoLambda(collision, V0s, tracks);
    for (auto& candidateV0 : candidateV0s) {
      // LOG(info) << candidateV0.pt;
      lambdaTableML(
        candidateV0.pt,
        candidateV0.eta,
        candidateV0.centFT0C,
        candidateV0.isMatter,
        candidateV0.mass,
        candidateV0.ct,
        candidateV0.radius,
        candidateV0.dcaV0PV,
        candidateV0.dcaPiPV,
        candidateV0.dcaPrPV,
        candidateV0.dcaV0Tracks,
        candidateV0.cosPA,
        candidateV0.tpcNsigmaPi,
        candidateV0.tpcNsigmaPr);
    }
  }
  PROCESS_SWITCH(LFStrangeTreeCreator, processData, "process data", false);

  void processMC(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions, TracksFull const& tracks, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLabels)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8())
        return;

      if (std::abs(collision.posZ()) > zVtxMax)
        return;

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollision, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      histos.fill(HIST("zVtx"), collision.posZ());

      if (downscaleFactor < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > downscaleFactor) {
        return;
      }

      fillMcLambda(collision, V0Table_thisCollision, tracks, mcParticles, mcLabels);

      for (auto& candidateV0 : candidateV0s) {
        if ((fillOnlySignal && std::abs(candidateV0.pdgCode) == 3122) || !fillOnlySignal) {
          mcLambdaTableML(
            candidateV0.pt,
            candidateV0.eta,
            candidateV0.centFT0C,
            candidateV0.isMatter,
            candidateV0.mass,
            candidateV0.ct,
            candidateV0.radius,
            candidateV0.dcaV0PV,
            candidateV0.dcaPiPV,
            candidateV0.dcaPrPV,
            candidateV0.dcaV0Tracks,
            candidateV0.cosPA,
            candidateV0.tpcNsigmaPi,
            candidateV0.tpcNsigmaPr,
            candidateV0.genPt,
            candidateV0.genEta,
            candidateV0.genCt,
            candidateV0.pdgCode,
            candidateV0.isReco);
        }
      }
    }

    if (fillNonRecSignal) {
      for (auto& mcPart : mcParticles) {

        if (std::abs(mcPart.pdgCode()) != 3122)
          continue;
        std::array<float, 3> secVtx;
        std::array<float, 3> primVtx = {mcPart.vx(), mcPart.vy(), mcPart.vz()};
        std::array<float, 3> momMother = {mcPart.px(), mcPart.py(), mcPart.pz()};
        for (auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
          if (std::abs(mcDaught.pdgCode()) == 2212) {
            secVtx = {mcDaught.vx(), mcDaught.vy(), mcDaught.vz()};
            break;
          }
        }
        if (std::find(checkedV0s.begin(), checkedV0s.end(), mcPart.globalIndex()) != std::end(checkedV0s)) {
          continue;
        }
        CandidateV0 candidateV0;
        auto momV0 = std::sqrt(std::pow(momMother[0], 2) + std::pow(momMother[1], 2) + std::pow(momMother[2], 2));
        candidateV0.genCt = std::sqrt(std::pow(secVtx[0] - primVtx[0], 2) + std::pow(secVtx[1] - primVtx[1], 2) + std::pow(secVtx[2] - primVtx[2], 2)) / (momV0 + 1e-10) * o2::constants::physics::MassLambda0;
        candidateV0.isReco = false;
        candidateV0.genPt = std::sqrt(std::pow(momMother[0], 2) + std::pow(momMother[1], 2));
        candidateV0.genEta = mcPart.eta();
        candidateV0.pdgCode = mcPart.pdgCode();
        mcLambdaTableML(
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          candidateV0.genPt,
          candidateV0.genEta,
          candidateV0.genCt,
          candidateV0.pdgCode,
          candidateV0.isReco);
      }
    }
  }
  PROCESS_SWITCH(LFStrangeTreeCreator, processMC, "process mc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LFStrangeTreeCreator>(cfgc)};
}
