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

/// \file CheckMCV0.cxx
/// \brief check the v0 phase-space
/// \author daiki.sekihata@cern.ch felix.schlepper@cern.ch
/// \dependencies: o2-analysis-lf-lambdakzeromcfinder

#include "PWGEM/PhotonMeson/DataModel/mcV0Tables.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"
//
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <TDatabasePDG.h>
#include <TMath.h>

#include <cmath>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::pwgem::photonmeson;

struct CheckMCV0 {
  // Output Objects
  Produces<aod::MCV0> mcV0Table;
  HistogramRegistry registry{"output", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Track selection
  Configurable<bool> ptLogAxis{"ptLogAxis", false, "Flag to use a log momentum axis"};
  Configurable<float> minpt{"minpt", 0.001, "min pt for track in GeV/c"};
  Configurable<float> maxpt{"maxpt", 20.0, "max pt for track in GeV/c"};
  Configurable<float> maxeta{"maxeta", 999.0, "eta acceptance"};
  Configurable<float> dcamin{"dcamin", 0.1, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};
  Configurable<float> maxZ{"maxZ", 200.0, "max z for track"};
  Configurable<float> maxX{"maxX", 200.0, "maximum X (starting point of track X)"};
  Configurable<float> maxY{"maxY", 200.0, "maximum Y (starting point of track Y)"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 100.0, "max chi2/NclsTPC"};
  Configurable<int> minCrossedRowsTPC{"minCrossedRowsTPC", 10, "min crossed rows tpc"};
  Configurable<bool> cutSign{"cutSign", false, "wrong sign cut"};
  Configurable<bool> cutSameSign{"cutSameSign", false, "reject same sign"};

  // Filters
  Filter trackPt = aod::track::pt > minpt&& aod::track::pt < maxpt;
  Filter trackZ = nabs(aod::track::z) < maxZ;
  Filter trackX = aod::track::x < maxX;
  Filter trackY = nabs(aod::track::y) < maxY;
  Filter trackEta = nabs(aod::track::eta) < maxeta;
  Filter trackDCA = nabs(aod::track::dcaXY) > dcamin&& nabs(aod::track::dcaXY) < dcamax;
  using TracksMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TracksCovIU, aod::McTrackLabels>;
  using FilteredTracksMC = soa::Filtered<TracksMC>;
  using CollisionsMC = soa::Join<aod::McCollisionLabels, aod::Collisions>;
  using V0s = aod::V0Datas;

  // Histogram Parameters
  Configurable<int> tglNBins{"tglNBins", 500, "nBins for tgl"};
  Configurable<int> zNBins{"zNBins", 2 * static_cast<int>(maxZ), "nBins for z"};
  Configurable<int> xNBins{"xNBins", static_cast<int>(maxX), "nBins for x"};
  Configurable<int> yNBins{"yNBins", static_cast<int>(maxY), "nBins for y"};
  Configurable<int> ptNBins{"ptNBins", 200, "nBins for pt"};
  Configurable<float> ptLowCut{"ptLowCut", 0.1, "low pt cut"};
  // Axes
  AxisSpec axisTgl{tglNBins, -3.f, +3.f, "tan(#lambda)"};
  AxisSpec axisTglEle{tglNBins, -3.f, +3.f, "tan(#lambda) of e^{#minus}"};
  AxisSpec axisTglPos{tglNBins, -3.f, +3.f, "tan(#lambda) of e^{#plus}"};
  AxisSpec axisEta{300, -1.6, +1.6, "#eta"};
  AxisSpec axisX{xNBins, 0.0, maxX, "reco. x"};
  AxisSpec axisXMC{xNBins, 0.0, maxX, "MC prop. x"};
  AxisSpec axisXEle{xNBins, 0.0, maxX, "reco. x of e^{#minus}"};
  AxisSpec axisXPos{xNBins, 0.0, maxX, "reco. x of e^{#plus}"};
  AxisSpec axisXMCEle{xNBins, 0.0, maxX, "MC prop. x of e^{#minus}"};
  AxisSpec axisXMCPos{xNBins, 0.0, maxX, "MC prop. x of e^{#plus}"};
  AxisSpec axisY{yNBins, -maxY, maxY, "reco. y"};
  AxisSpec axisYMC{yNBins, -maxY, maxY, "MC prop. y"};
  AxisSpec axisYEle{yNBins, -maxY, maxY, "reco. y of e^{#minus}"};
  AxisSpec axisYPos{yNBins, -maxY, maxY, "reco. y of e^{#plus}"};
  AxisSpec axisYMCEle{yNBins, -maxY, maxY, "MC prop. y of e^{#minus}"};
  AxisSpec axisYMCPos{yNBins, -maxY, maxY, "MC prop. y of e^{#plus}"};
  AxisSpec axisZ{zNBins, -maxZ, maxZ, "reco. z"};
  AxisSpec axisZMC{zNBins, -maxZ, maxZ, "MC prop.z"};
  AxisSpec axisZEle{zNBins, -maxZ, maxZ, "reco. z of e^{#minus}"};
  AxisSpec axisZPos{zNBins, -maxZ, maxZ, "reco. z of e^{#plus}"};
  AxisSpec axisZMCEle{zNBins, -maxZ, maxZ, "MC prop. z of e^{#minus}"};
  AxisSpec axisZMCPos{zNBins, -maxZ, maxZ, "MC prop. z of e^{#plus}"};
  AxisSpec axisXDiff{xNBins, -maxX, maxX, "reco. X of e^{#plus} - reco. X of e^{#minus}"};
  AxisSpec axisYDiff{yNBins, -maxY, maxY, "reco. y of e^{#plus} - reco. y of e^{#minus}"};
  AxisSpec axisZDiff{zNBins, -maxZ, maxZ, "reco. z of e^{#plus} - reco. z of e^{#minus}"};
  AxisSpec axisRMCDiff{xNBins, -30, 30, "#delta(R_{reco} - R_{MC})"};
  AxisSpec axisXMCDiff{xNBins, -30, 30, "#delta(X_{reco} - X_{MC})"};
  AxisSpec axisYMCDiff{yNBins, -30, 30, "#delta(Y_{reco} - Y_{MC})"};
  AxisSpec axisZMCDiff{zNBins, -30, 30, "#delta(Z_{reco} - Z_{MC})"};
  AxisSpec axisPhi{tglNBins, 0, TMath::Pi() * 2, "#phi"};
  AxisSpec axisR{tglNBins, 0, 100, "R_{gen}"};
  AxisSpec axisPt{ptNBins, minpt, maxpt, "#it{p}_{T} GeV/#it{c}"};
  AxisSpec axisPtPos{ptNBins, minpt, maxpt, "#it{p}_{T} GeV/#it{c} of e^{#plus}"};
  AxisSpec axisPtEle{ptNBins, minpt, maxpt, "#it{p}_{T} GeV/#it{c} of e^{#minus}"};
  AxisSpec axisTPos{1000, -50'000, +50'000, "timeStamp ns of e^{#plus}"};
  AxisSpec axisTEle{1000, -50'000, +50'000, "timeStamp ns of e^{#minus}"};

  static constexpr std::array<std::string_view, 9> cutsBinLabels{"Pre", "checkV0leg", "sign", "sameSign", "photon", "propagationFailed", "lowPt", "highPt", "survived"};
  enum cutsBinEnum : uint8_t {
    PRE = 1,
    CHECKV0LEG,
    SIGN,
    SAMESIGN,
    PHOTON,
    PROPFAIL,
    LOWPT,
    HIGHPT,
    SURVIVED,
  };
  static_assert(cutsBinLabels.size() == cutsBinEnum::SURVIVED);

  static constexpr std::array<std::string_view, 8> checkV0legLabels{"Pt<minPt", "Pt>maxPt", "dca<dcamin", "dca>dcamax", "tpcChi2NCl>maxchi2tpc", "no ITS||TPC", "eta>maxEta", "tpcCrossedRows<minCrossedRowsTPC"};
  enum checkV0legEnum : uint8_t {
    PTMIN = 1,
    PTMAX,
    DCAMIN,
    DCAMAX,
    TPCCHI2NCL,
    NOITSTPC,
    MAXETA,
    MINCROSSEDROWSTPC,
  };
  static_assert(checkV0legLabels.size() == checkV0legEnum::MINCROSSEDROWSTPC);

  // CCDB
  Configurable<std::string> mCCDBPath{"ccdb-path", "GLO/GRP/GRP", "path to the ccdb object"};
  Configurable<std::string> mGRPMagPath{"grpmagPath", "GLO/Config/GRPMagField", "path to the GRPMagField object"};
  Configurable<std::string> mLUTPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<std::string> mCCDBUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> mCCDB;
  int mRunNumber{-1};
  o2::base::MatLayerCylSet* mLUT{nullptr};
  o2::parameters::GRPMagField* mGRPMagField{nullptr};

  // params
  std::array<float, 6> mcPosXYZEtaTglPtProp{};
  std::array<float, 6> mcEleXYZEtaTglPtProp{};
  std::array<float, 6> mcMotherXYZEtaTglPtProp{};

  // Track Types
  static constexpr std::array<std::string_view, 6> v0Types{"ITSTPC_ITSTPC/", "TPConly_TPConly/", "ITSonly_ITSonly/", "ITSTPC_TPConly/", "ITSTPC_ITSonly/", "TPConly_ITSonly/"};
  std::array<bool, v0Types.size()> v0TypesPassed{};
  static constexpr std::array<std::string_view, 4> cuts{"lowPt/", "highPt/", "zReco43/", "all/"};

  void init(InitContext const& /*unused*/)
  {
    // setup CCDB
    mCCDB->setURL(mCCDBUrl);
    mCCDB->setCaching(true);
    mCCDB->setLocalObjectValidityChecking();
    mLUT = o2::base::MatLayerCylSet::rectifyPtrFromFile(mCCDB->get<o2::base::MatLayerCylSet>(mLUTPath));

    // maybe logarithmic
    if (ptLogAxis) {
      axisPt.makeLogarithmic();
      axisPtPos.makeLogarithmic();
      axisPtEle.makeLogarithmic();
    }

    // Create histograms
    for (const auto& v0type : v0Types) {
      for (const auto& cut : cuts) {
        const auto path = Form("%s%s", v0type.data(), cut.data());
        registry.add(Form("%sTglTgl", path), "tan(#lambda) vs. tan(#lambda)", HistType::kTH2F, {axisTglPos, axisTglEle});
        registry.add(Form("%sPtPt", path), "p_{T} vs. p_{T}", HistType::kTH2F, {axisPtPos, axisPtEle});
        registry.add(Form("%sXX", path), "x vs. x", HistType::kTH2F, {axisXPos, axisXEle});
        registry.add(Form("%sYY", path), "y vs. y", HistType::kTH2F, {axisYPos, axisYEle});
        registry.add(Form("%sZZ", path), "z vs. z", HistType::kTH2F, {axisZPos, axisZEle});
        registry.add(Form("%sZPt", path), "z vs. p_{T}", HistType::kTH2F, {axisZ, axisPt});
        registry.add(Form("%sZTgl", path), "z vs. tan(#lambda)", HistType::kTH2F, {axisZ, axisTgl});
        registry.add(Form("%sPtTgl", path), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {axisPt, axisTgl});
        registry.add(Form("%sTT", path), "est. Time vs. Time", HistType::kTH2F, {axisTPos, axisTEle});
        registry.add(Form("%sMC/XX", path), Form("MC x vs. reconstructed x (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {axisXMC, axisX});
        registry.add(Form("%sMC/YY", path), Form("MC y vs. rconstructed y (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {axisYMC, axisY});
        registry.add(Form("%sMC/ZZ", path), Form("MC z vs. reconstructed z (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {axisZMC, axisZ});
        registry.add(Form("%sMC/VertexPropagationX", path), Form("MC vertex X propagated to track (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {{axisXMCPos, axisXMCEle}});
        registry.add(Form("%sMC/VertexPropagationY", path), Form("MC vertex Y propagated to track (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {{axisYMCPos, axisYMCEle}});
        registry.add(Form("%sMC/VertexPropagationZ", path), Form("MC vertex Z propagated to track (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {{axisZMCPos, axisZMCEle}});
      }
    }

    registry.add("check/DeltaRPhi", "#deltaR", HistType::kTH2F, {axisPhi, axisRMCDiff});
    registry.add("check/DeltaRR", "#deltaR", HistType::kTH2F, {axisR, axisRMCDiff});
    for (const auto& v0type : v0Types) {
      registry.add(Form("check/%sDeltaRPhi", v0type.data()), "#deltaR", HistType::kTH2F, {axisPhi, axisRMCDiff});
      registry.add(Form("check/%sDeltaXPhi", v0type.data()), "#deltaX", HistType::kTH2F, {axisPhi, axisXMCDiff});
      registry.add(Form("check/%sDeltaYPhi", v0type.data()), "#deltaY", HistType::kTH2F, {axisPhi, axisYMCDiff});
      registry.add(Form("check/%sDeltaZPhi", v0type.data()), "#deltaZ", HistType::kTH2F, {axisPhi, axisZMCDiff});
    }

    registry.add("V0Counter", "V0 counter", HistType::kTH1F, {{cutsBinLabels.size(), 0.5, 0.5 + cutsBinLabels.size()}});
    for (size_t iBin = 0; iBin < cutsBinLabels.size(); ++iBin) {
      registry.get<TH1>(HIST("V0Counter"))->GetXaxis()->SetBinLabel(iBin + 1, cutsBinLabels[iBin].data());
    }

    registry.add("V0TypeCounter", "V0 Type counter", HistType::kTH1F, {{v0Types.size(), 0.5, 0.5 + v0Types.size()}});
    for (size_t iBin = 0; iBin < v0Types.size(); ++iBin) {
      registry.get<TH1>(HIST("V0TypeCounter"))->GetXaxis()->SetBinLabel(iBin + 1, v0Types[iBin].data());
    }

    registry.add("CheckV0Leg", "CheckV0Leg", HistType::kTH1F, {{checkV0legLabels.size(), 0.5, 0.5 + checkV0legLabels.size()}});
    for (size_t iBin = 0; iBin < checkV0legLabels.size(); ++iBin) {
      registry.get<TH1>(HIST("CheckV0Leg"))->GetXaxis()->SetBinLabel(iBin + 1, checkV0legLabels[iBin].data());
    }
  }

  Preslice<aod::V0Datas> perCollision = aod::v0data::collisionId;
  void processMCV0(CollisionsMC const& collisions, V0s const& v0s, FilteredTracksMC const& /*unused*/, aod::McParticles const& /*unused*/, aod::McCollisions const& /*unused*/, aod::BCsWithTimestamps const& /*unused*/)
  {
    // Check for new ccdb parameters
    for (auto& collision : collisions) {
      const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      // Get the V0 candidates belonging to the current collision
      for (const auto& v0 : v0s.sliceBy(perCollision, static_cast<int>(collision.globalIndex()))) {
        registry.fill(HIST("V0Counter"), cutsBinEnum::PRE);

        // tracks
        const auto& pos = v0.template posTrack_as<TracksMC>(); // positive daughter
        const auto& ele = v0.template negTrack_as<TracksMC>(); // negative daughter
        if (!checkV0leg(pos) || !checkV0leg(ele)) {
          registry.fill(HIST("V0Counter"), cutsBinEnum::CHECKV0LEG);
          continue;
        }

        // set correct track types
        checkPassed(pos, ele);

        if (cutSign && (std::signbit(pos.sign()) || !std::signbit(ele.sign()))) { // wrong sign and same sign reject (e.g. loopers)
          registry.fill(HIST("V0Counter"), cutsBinEnum::SIGN);
          continue;
        }
        if (cutSameSign && (pos.sign() * ele.sign() > 0)) {
          registry.fill(HIST("V0Counter"), cutsBinEnum::SAMESIGN);
          continue;
        }
        const auto& posMC = pos.template mcParticle_as<aod::McParticles>();
        const auto& eleMC = ele.template mcParticle_as<aod::McParticles>();
        if (!checkMCParticles<kGamma>(posMC, eleMC)) {
          registry.fill(HIST("V0Counter"), cutsBinEnum::PHOTON);
          continue;
        }
        const auto& mother = posMC.mothers_first_as<aod::McParticles>();

        // Propagate the vertex position of the MC track to the reconstructed vertex position of the track
        if (!recacluateMCVertex(ele, eleMC, mcEleXYZEtaTglPtProp) || !recacluateMCVertex(pos, posMC, mcPosXYZEtaTglPtProp)) {
          registry.fill(HIST("V0Counter"), cutsBinEnum::PROPFAIL);
          continue;
        }
        registry.fill(HIST("V0Counter"), SURVIVED);

        float mcR = std::sqrt(posMC.vx() * posMC.vx() + posMC.vy() * posMC.vy());
        float dr = v0.v0radius() - mcR;
        float dx = v0.x() - posMC.vx();
        float dy = v0.y() - posMC.vy();
        float dz = v0.z() - posMC.vz();

        registry.fill(HIST("check/DeltaRPhi"), v0.phi(), dr);
        registry.fill(HIST("check/DeltaRR"), mcR, dr);

        // fill histograms
        static_for<0, v0Types.size() - 1>([&](auto i_idx) {
          constexpr unsigned int i = i_idx.value;
          if (!v0TypesPassed[i]) {
            return;
          }
          registry.fill(HIST("check/") + HIST(v0Types[i]) + HIST("DeltaRPhi"), v0.phi(), dr);
          registry.fill(HIST("check/") + HIST(v0Types[i]) + HIST("DeltaXPhi"), v0.phi(), dx);
          registry.fill(HIST("check/") + HIST(v0Types[i]) + HIST("DeltaYPhi"), v0.phi(), dy);
          registry.fill(HIST("check/") + HIST(v0Types[i]) + HIST("DeltaZPhi"), v0.phi(), dz);

          static_for<0, cuts.size() - 1>([&](auto j_idx) {
            constexpr unsigned int j = j_idx.value;
            fillHistograms<i, j>(pos, ele);
          });
        });

        mcV0Table(
          // V0
          v0.x(), v0.y(), v0.z(),
          v0.px(), v0.py(), v0.pz(),
          // Track Pos
          pos.pt(),
          pos.eta(), pos.tgl(), pos.x(), pos.y(), pos.z(), pos.trackTime(), pos.sign(),
          pos.hasITS(), pos.hasTPC(), pos.hasTRD(), pos.hasTOF(),
          // Track Ele
          ele.pt(), ele.eta(), ele.tgl(), ele.x(), ele.y(), ele.z(), ele.trackTime(), ele.sign(),
          ele.hasITS(), ele.hasTPC(), ele.hasTRD(), ele.hasTOF(),
          // MC particle Pos
          posMC.pt(), posMC.eta(), posMC.vx(), posMC.vy(), posMC.vz(),
          // MC particle Ele
          eleMC.pt(), eleMC.eta(), eleMC.vx(), eleMC.vy(), eleMC.vz(),
          // Propagated MC pos
          mcPosXYZEtaTglPtProp[0], mcPosXYZEtaTglPtProp[1], mcPosXYZEtaTglPtProp[2],
          mcPosXYZEtaTglPtProp[3], mcPosXYZEtaTglPtProp[4], mcPosXYZEtaTglPtProp[5],
          // Propagated MC ele
          mcEleXYZEtaTglPtProp[0], mcEleXYZEtaTglPtProp[1], mcEleXYZEtaTglPtProp[2],
          mcEleXYZEtaTglPtProp[3], mcEleXYZEtaTglPtProp[4], mcEleXYZEtaTglPtProp[5],
          // MC mother particle
          mother.isPhysicalPrimary(), mother.producedByGenerator(), mother.pt(), mother.eta(),
          mother.vx(), mother.vy(), mother.vz(),
          // MC Primary Vertex
          mother.vx(), mother.vy(), mother.vz(),
          // TrackType
          std::distance(v0TypesPassed.cbegin(), std::find(v0TypesPassed.cbegin(), v0TypesPassed.cend(), true)));
      }
    }
  }
  PROCESS_SWITCH(CheckMCV0, processMCV0, "process reconstructed MC V0 info", true);

  template <unsigned int idxV0Type, unsigned int idxCut, typename TTrack>
  inline void fillHistograms(TTrack const& pos, TTrack const& ele)
  {
    constexpr auto hist = HIST(v0Types[idxV0Type]) + HIST(cuts[idxCut]);

    if constexpr (idxCut == 0) { // lowPt
      if (pos.pt() > ptLowCut && ele.pt() > ptLowCut) {
        return;
      }
    } else if constexpr (idxCut == 1) { // highPt
      if (pos.pt() < ptLowCut && ele.pt() < ptLowCut) {
        return;
      }
    } else if constexpr (idxCut == 2) {
      if (!(std::abs(pos.z()) > 38.f && std::abs(pos.z()) < 46.f) && !(std::abs(ele.z()) > 38.f && std::abs(ele.z()) < 46.f)) {
        return;
      }
    }

    registry.fill(hist + HIST("TglTgl"), pos.tgl(), ele.tgl());
    registry.fill(hist + HIST("PtPt"), pos.pt(), ele.pt());
    registry.fill(hist + HIST("XX"), pos.x(), ele.x());
    registry.fill(hist + HIST("YY"), pos.y(), ele.y());
    registry.fill(hist + HIST("ZZ"), pos.z(), ele.z());
    registry.fill(hist + HIST("ZPt"), pos.z(), pos.pt());
    registry.fill(hist + HIST("ZPt"), ele.z(), ele.pt());
    registry.fill(hist + HIST("ZTgl"), pos.z(), pos.tgl());
    registry.fill(hist + HIST("ZTgl"), ele.z(), ele.tgl());
    registry.fill(hist + HIST("PtTgl"), pos.pt(), pos.tgl());
    registry.fill(hist + HIST("PtTgl"), ele.pt(), ele.tgl());
    registry.fill(hist + HIST("TT"), pos.trackTime(), ele.trackTime());
    if constexpr (idxCut == 0) {
      if (pos.pt() < ptLowCut || ele.pt() < ptLowCut) {
        registry.fill(HIST("V0Counter"), cutsBinEnum::LOWPT);
      }
      if (pos.pt() < ptLowCut) {
        registry.fill(hist + HIST("MC/") + HIST("XX"), mcPosXYZEtaTglPtProp[0], pos.x());
        registry.fill(hist + HIST("MC/") + HIST("YY"), mcPosXYZEtaTglPtProp[1], pos.y());
        registry.fill(hist + HIST("MC/") + HIST("ZZ"), mcPosXYZEtaTglPtProp[2], pos.z());
      }
      if (ele.pt() < ptLowCut) {
        registry.fill(hist + HIST("MC/") + HIST("XX"), mcEleXYZEtaTglPtProp[0], ele.x());
        registry.fill(hist + HIST("MC/") + HIST("YY"), mcEleXYZEtaTglPtProp[1], ele.y());
        registry.fill(hist + HIST("MC/") + HIST("ZZ"), mcEleXYZEtaTglPtProp[2], ele.z());
      }
      registry.fill(hist + HIST("MC/") + HIST("VertexPropagationX"), mcPosXYZEtaTglPtProp[0], mcEleXYZEtaTglPtProp[0]);
      registry.fill(hist + HIST("MC/") + HIST("VertexPropagationY"), mcPosXYZEtaTglPtProp[1], mcEleXYZEtaTglPtProp[1]);
      registry.fill(hist + HIST("MC/") + HIST("VertexPropagationZ"), mcPosXYZEtaTglPtProp[2], mcEleXYZEtaTglPtProp[2]);
    } else if constexpr (idxCut == 1) {
      if (pos.pt() > ptLowCut || ele.pt() > ptLowCut) {
        registry.fill(HIST("V0Counter"), cutsBinEnum::HIGHPT);
      }
      if (pos.pt() > ptLowCut) {
        registry.fill(hist + HIST("MC/") + HIST("XX"), mcPosXYZEtaTglPtProp[0], pos.x());
        registry.fill(hist + HIST("MC/") + HIST("YY"), mcPosXYZEtaTglPtProp[1], pos.y());
        registry.fill(hist + HIST("MC/") + HIST("ZZ"), mcPosXYZEtaTglPtProp[2], pos.z());
      }
      if (ele.pt() > ptLowCut) {
        registry.fill(hist + HIST("MC/") + HIST("XX"), mcEleXYZEtaTglPtProp[0], ele.x());
        registry.fill(hist + HIST("MC/") + HIST("YY"), mcEleXYZEtaTglPtProp[1], ele.y());
        registry.fill(hist + HIST("MC/") + HIST("ZZ"), mcEleXYZEtaTglPtProp[2], ele.z());
      }
      registry.fill(hist + HIST("MC/") + HIST("VertexPropagationX"), mcPosXYZEtaTglPtProp[0], mcEleXYZEtaTglPtProp[0]);
      registry.fill(hist + HIST("MC/") + HIST("VertexPropagationY"), mcPosXYZEtaTglPtProp[1], mcEleXYZEtaTglPtProp[1]);
      registry.fill(hist + HIST("MC/") + HIST("VertexPropagationZ"), mcPosXYZEtaTglPtProp[2], mcEleXYZEtaTglPtProp[2]);
    } else {
      registry.fill(hist + HIST("MC/") + HIST("XX"), mcPosXYZEtaTglPtProp[0], pos.x());
      registry.fill(hist + HIST("MC/") + HIST("YY"), mcPosXYZEtaTglPtProp[1], pos.y());
      registry.fill(hist + HIST("MC/") + HIST("ZZ"), mcPosXYZEtaTglPtProp[2], pos.z());
      registry.fill(hist + HIST("MC/") + HIST("XX"), mcEleXYZEtaTglPtProp[0], ele.x());
      registry.fill(hist + HIST("MC/") + HIST("YY"), mcEleXYZEtaTglPtProp[1], ele.y());
      registry.fill(hist + HIST("MC/") + HIST("ZZ"), mcEleXYZEtaTglPtProp[2], ele.z());
      registry.fill(hist + HIST("MC/") + HIST("VertexPropagationX"), mcPosXYZEtaTglPtProp[0], mcEleXYZEtaTglPtProp[0]);
      registry.fill(hist + HIST("MC/") + HIST("VertexPropagationY"), mcPosXYZEtaTglPtProp[1], mcEleXYZEtaTglPtProp[1]);
      registry.fill(hist + HIST("MC/") + HIST("VertexPropagationZ"), mcPosXYZEtaTglPtProp[2], mcEleXYZEtaTglPtProp[2]);
    }
  }

  // Templates
  template <typename TTrack>
  inline bool checkV0leg(TTrack const& track)
  {
    if (track.pt() < minpt) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::PTMIN);
      return false;
    }
    if (track.pt() > maxpt) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::PTMAX);
      return false;
    }
    if (std::abs(track.eta()) > maxeta) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::MAXETA);
      return false;
    }
    if (std::abs(track.dcaXY()) < dcamin) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::DCAMIN);
      return false;
    }
    if (std::abs(track.dcaXY()) > dcamax) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::DCAMAX);
      return false;
    }
    if (track.tpcChi2NCl() > maxChi2TPC) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::TPCCHI2NCL);
      return false;
    }
    if (track.tpcNClsCrossedRows() < minCrossedRowsTPC) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::MINCROSSEDROWSTPC);
      return false;
    }
    if (track.trackTimeRes() <= o2::constants::lhc::LHCBunchSpacingNS) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  inline void checkPassed(TTrack const& track0, TTrack const& track1)
  {
    v0TypesPassed.fill(false); // reset
    v0TypesPassed[0] = isITSTPC_ITSTPC(track0, track1);
    v0TypesPassed[1] = isTPConly_TPConly(track0, track1);
    v0TypesPassed[2] = isITSonly_ITSonly(track0, track1);
    v0TypesPassed[3] = isITSTPC_TPConly(track0, track1);
    v0TypesPassed[4] = isITSTPC_ITSonly(track0, track1);
    v0TypesPassed[5] = isTPConly_ITSonly(track0, track1);
    for (size_t i = 0; i < v0TypesPassed.size(); ++i) {
      if (v0TypesPassed[i]) {
        registry.fill(HIST("V0TypeCounter"), i + 1);
      }
    }
  }

  template <typename TTrack, typename MCTrack>
  inline bool recacluateMCVertex(TTrack const& track, MCTrack const& mcTrack, std::array<float, 6>& xyzEtaTglPt)
  {
    std::array<float, 3> xyzMC{mcTrack.vx(), mcTrack.vy(), mcTrack.vz()};
    std::array<float, 3> pxyzMC{mcTrack.px(), mcTrack.py(), mcTrack.pz()};
    auto pPDG = TDatabasePDG::Instance()->GetParticle(mcTrack.pdgCode());
    if (!pPDG) {
      return false;
    }
    o2::track::TrackPar mctrO2(xyzMC, pxyzMC, TMath::Nint(pPDG->Charge() / 3), false);
    // rotate track into padplane of detector and then propagate to corresponding x of track
    if (!mctrO2.rotate(track.alpha()) || !o2::base::Propagator::Instance()->PropagateToXBxByBz(mctrO2, track.x())) {
      return false;
    }
    xyzEtaTglPt = {mctrO2.getX(), mctrO2.getY(), mctrO2.getZ(), mctrO2.getEta(), mctrO2.getTgl(), mctrO2.getPt()};
    return true;
  }

  template <typename BC>
  inline void initCCDB(BC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mGRPMagField = mCCDB->getForTimeStamp<o2::parameters::GRPMagField>(mGRPMagPath, bc.timestamp());
    o2::base::Propagator::initFieldFromGRP(mGRPMagField);
    o2::base::Propagator::Instance()->setMatLUT(mLUT);
    mRunNumber = bc.runNumber();
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CheckMCV0>(cfgc, TaskName{"check-mc-v0"})};
}
