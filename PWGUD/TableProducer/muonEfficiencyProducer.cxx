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

/// \file muonEfficiencyProducer.cxx
/// \brief Compact forward-muon skim for offline MCH and MCH-MID efficiency studies
/// \author Nazar Burmasov (JINR)

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHTracking/TrackExtrap.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>

#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>
#include <TGeoManager.h>

#include <array>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace o2::aod::mcheffskim
{
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(ReferenceBC, referenceBC, int64_t);
DECLARE_SOA_COLUMN(FT0AAmplitude, ft0aAmplitude, float);
DECLARE_SOA_COLUMN(FT0CAmplitude, ft0cAmplitude, float);
DECLARE_SOA_COLUMN(FV0AAmplitude, fv0aAmplitude, float);
DECLARE_SOA_COLUMN(CollisionId, collisionId, int32_t); // -1 for ambiguous tracks
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(RAbs, rAbs, float);
DECLARE_SOA_COLUMN(PDca, pDca, float);
DECLARE_SOA_COLUMN(Chi2Ndf, chi2Ndf, float); // -1 if propagation failed
DECLARE_SOA_COLUMN(NClusters, nClusters, int8_t);
DECLARE_SOA_COLUMN(MchBitMap, mchBitMap, uint16_t);
DECLARE_SOA_COLUMN(MidBitMap, midBitMap, uint8_t);
DECLARE_SOA_COLUMN(MidBoards, midBoards, uint32_t);
DECLARE_SOA_COLUMN(Chi2MatchMCHMID, chi2MatchMCHMID, float); // -1 for MCH-only tracks
DECLARE_SOA_COLUMN(TrackTime, trackTime, float);             // ns relative to referenceBC
DECLARE_SOA_COLUMN(TrackTimeRes, trackTimeRes, float);
} // namespace o2::aod::mcheffskim

namespace o2::aod
{
DECLARE_SOA_TABLE(MchEffBCs, "AOD", "MCHEFFBC",
                  o2::soa::Index<>,
                  mcheffskim::RunNumber, mcheffskim::ReferenceBC,
                  mcheffskim::FT0AAmplitude, mcheffskim::FT0CAmplitude,
                  mcheffskim::FV0AAmplitude);
} // namespace o2::aod

namespace o2::aod::mcheffskim
{
DECLARE_SOA_INDEX_COLUMN(MchEffBC, mchEffBC); // links a muon to the compact BC table
} // namespace o2::aod::mcheffskim

namespace o2::aod
{
DECLARE_SOA_TABLE(MchEffSkims, "AOD", "MCHEFFSKIM",
                  o2::soa::Index<>, mcheffskim::MchEffBCId,
                  mcheffskim::CollisionId,
                  mcheffskim::Pt, mcheffskim::Eta, mcheffskim::Phi, mcheffskim::Sign,
                  mcheffskim::RAbs, mcheffskim::PDca, mcheffskim::Chi2Ndf,
                  mcheffskim::NClusters, mcheffskim::MchBitMap,
                  mcheffskim::MidBitMap, mcheffskim::MidBoards,
                  mcheffskim::Chi2MatchMCHMID,
                  mcheffskim::TrackTime, mcheffskim::TrackTimeRes);
} // namespace o2::aod

struct MuonEfficiencyProducer {
  static constexpr float MinEta = -4.f;
  static constexpr float MaxEta = -2.5f;
  static constexpr float NoAmplitude = -1.f;
  static constexpr float MaxFITTime = 30.f;

  using Tracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
  using Matrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  using Vector5 = ROOT::Math::SVector<double, 5>;

  Produces<aod::MchEffBCs> outputBCs;
  Produces<aod::MchEffSkims> output;
  Service<ccdb::BasicCCDBManager> ccdb{};
  globaltracking::MatchGlobalFwd matcher;
  int cachedFieldRun = -1;

  Configurable<float> minPt{"minPt", 0.4f, "Minimum track pT (Gev/c)"};
  Configurable<float> minRAbs{"minRAbs", 17.6f, "Minimum absorber-end radius (cm)"};
  Configurable<float> maxRAbs{"maxRAbs", 89.5f, "Maximum absorber-end radius (cm)"};
  Configurable<float> maxPDca{"maxPDca", 400.f, "Maximum stored pDCA (GeV cm/c)"};
  Configurable<float> maxChi2{"maxChi2", 100.f, "Maximum MCH chi2/ndf"};

  void init(InitContext&)
  {
    ccdb->setLocalObjectValidityChecking();
  }

  void setField(int run)
  {
    if (cachedFieldRun == run) {
      return;
    }
    const auto runInfo = parameters::AggregatedRunInfo::buildAggregatedRunInfo(ccdb->instance(), run);
    const auto ts = (runInfo.sor + runInfo.eor) / 2;
    auto* grp = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", ts);
    if (!grp) {
      LOGF(fatal, "GRPMagField not found for run %d", run);
    }
    base::Propagator::initFieldFromGRP(grp);
    if (!base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");
    }
    mch::TrackExtrap::setField();
    cachedFieldRun = run;
  }

  struct Kinematics {
    float pt, eta, phi;
  };

  template <typename Track>
  bool propagateMCH(const Track& t, const std::array<float, 3>& vtx, const std::array<float, 2>& vtxCov, Kinematics& k)
  {
    const Vector5 params(t.x(), t.y(), t.phi(), t.tgl(), t.signed1Pt());
    const std::vector<double> covElements{t.cXX(), t.cXY(), t.cYY(), t.cPhiX(), t.cPhiY(), t.cPhiPhi(),
                                          t.cTglX(), t.cTglY(), t.cTglPhi(), t.cTglTgl(),
                                          t.c1PtX(), t.c1PtY(), t.c1PtPhi(), t.c1PtTgl(), t.c1Pt21Pt2()};
    const Matrix55 cov(covElements.begin(), covElements.end());
    dataformats::GlobalFwdTrack source;
    source.setParameters(params);
    source.setZ(t.z());
    source.setCovariances(cov);
    auto mchTrack = matcher.FwdtoMCH(source);
    if (!mch::TrackExtrap::extrapToVertex(mchTrack, vtx[0], vtx[1], vtx[2], vtxCov[0], vtxCov[1])) {
      return false;
    }
    const auto propagated = matcher.MCHtoFwd(mchTrack);
    k = {.pt = static_cast<float>(propagated.getPt()), .eta = static_cast<float>(propagated.getEta()), .phi = static_cast<float>(propagated.getPhi())};
    return true;
  }

  void process(aod::BCs const& bcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as,
               Tracks const& tracks, aod::Collisions const& collisions,
               aod::AmbiguousFwdTracks const& ambiguous)
  {
    if (bcs.size() == 0) {
      return;
    }
    const int run = bcs.begin().runNumber();

    std::vector<float> ft0aByBC(bcs.size(), NoAmplitude);
    std::vector<float> ft0cByBC(bcs.size(), NoAmplitude);
    std::vector<float> fv0aByBC(bcs.size(), NoAmplitude);
    std::vector<int32_t> skimBCIds(bcs.size(), -1);
    for (const auto& ft0 : ft0s) {
      if (ft0.timeA() < MaxFITTime) {
        ft0aByBC[ft0.bcId()] = ft0.sumAmpA();
      }
      if (ft0.timeC() < MaxFITTime) {
        ft0cByBC[ft0.bcId()] = ft0.sumAmpC();
      }
    }
    for (const auto& fv0 : fv0as) {
      if (fv0.time() >= MaxFITTime) {
        continue;
      }
      float amplitude = 0.f;
      for (const float a : fv0.amplitude()) {
        if (a > 0.f) {
          amplitude += a;
        }
      }
      fv0aByBC[fv0.bcId()] = amplitude;
    }

    std::vector<int32_t> ambiguousBCIds(tracks.size(), -1);
    for (const auto& entry : ambiguous) {
      const auto& ids = entry.bcIds();
      if (ids.empty()) {
        continue;
      }
      ambiguousBCIds[entry.fwdtrackId()] = *ids.begin();
    }

    for (const auto& t : tracks) {
      const auto type = t.trackType();
      if (type != aod::fwdtrack::MuonStandaloneTrack && type != aod::fwdtrack::MCHStandaloneTrack) {
        continue;
      }
      const float rAbs = t.rAtAbsorberEnd();
      if (rAbs < minRAbs || rAbs > maxRAbs || t.pDca() > maxPDca || t.chi2() > maxChi2) {
        continue;
      }
      Kinematics k{.pt = t.pt(), .eta = t.eta(), .phi = t.phi()};
      std::array<float, 3> vtx{0.f, 0.f, 0.f};
      std::array<float, 2> vtxCov{0.f, 0.f};
      int32_t referenceBCId = ambiguousBCIds[t.globalIndex()];
      if (t.has_collision()) {
        const auto col = collisions.iteratorAt(t.collisionId());
        vtx = {col.posX(), col.posY(), col.posZ()};
        vtxCov = {col.covXX(), col.covYY()};
        referenceBCId = col.bcId();
      }
      setField(run);
      if (!propagateMCH(t, vtx, vtxCov, k)) {
        continue;
      }
      if (k.eta < MinEta || k.eta > MaxEta) {
        continue;
      }
      if (k.pt < minPt) {
        continue;
      }
      int32_t skimBCId = -1;
      if (referenceBCId >= 0) {
        skimBCId = skimBCIds[referenceBCId];
        if (skimBCId < 0) {
          outputBCs(run, bcs.iteratorAt(referenceBCId).globalBC(),
                    ft0aByBC[referenceBCId], ft0cByBC[referenceBCId], fv0aByBC[referenceBCId]);
          skimBCId = outputBCs.lastIndex();
          skimBCIds[referenceBCId] = skimBCId;
        }
      }
      output(skimBCId, t.collisionId(), k.pt, k.eta, k.phi, static_cast<int8_t>(t.signed1Pt() > 0 ? 1 : -1),
             rAbs, t.pDca(), t.chi2(), t.nClusters(), t.mchBitMap(), t.midBitMap(), t.midBoards(),
             type == aod::fwdtrack::MuonStandaloneTrack ? t.chi2MatchMCHMID() : -1.f, t.trackTime(), t.trackTimeRes());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<MuonEfficiencyProducer>(context)};
}
