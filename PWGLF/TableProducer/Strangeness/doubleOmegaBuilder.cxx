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

/// \file doubleOmegaBuilder.cxx
/// \brief Builds the Omega, Lambda and direct-kaon candidates of the double-Omega analysis and stores
///        them in derived tables, consumed on the fly (or from saved derived data) by doubleOmegaTreeCreator
///        for the same-event candidates and the event mixing.

#include "PWGLF/DataModel/LFDoubleOmegaTables.h"
#include "PWGLF/Utils/doubleOmegaBuilderHelper.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TH2.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>

using namespace o2;
using namespace o2::framework;
using namespace o2::pwglf::doubleomega;

using CollisionsTable = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTOFKa>;

struct doubleOmegaBuilder {
  Produces<aod::DOmegaColls> dOmegaColls;
  Produces<aod::DOmegaTracks> dOmegaTracks;
  Produces<aod::DOmegaOmegas> dOmegaOmegas;
  Produces<aod::DOmegaLambdas> dOmegaLambdas;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  enum CollisionStep : uint8_t {
    kAllCollisions,
    kSelectedCollisions,
    kWithOmega,
    kStoredCollisions,
    kNCollisionSteps
  };

  Preslice<TracksFull> tracksPerCollision = aod::track::collisionId;
  Preslice<aod::Cascades> cascadesPerCollision = aod::cascade::collisionId;
  Preslice<aod::V0s> v0sPerCollision = aod::v0::collisionId;

  DoubleOmegaBuilderHelper helper;
  int mRunNumber = 0;
  float mBz = 0.f;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
  Configurable<bool> cfgApplyEventSelection{"cfgApplyEventSelection", true, "Apply the standard collision event selection"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Material correction for the V0/cascade fits"};

  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis massOmegaAxis{"massOmegaAxis", {400, o2::constants::physics::MassOmegaMinus - 0.05f, o2::constants::physics::MassOmegaMinus + 0.05f}, "binning for the Omega invariant-mass"};
  ConfigurableAxis massXiAxis{"massXiAxis", {400, o2::constants::physics::MassXiMinus - 0.05f, o2::constants::physics::MassXiMinus + 0.05f}, "binning for the Xi invariant-mass"};
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, o2::constants::physics::MassLambda0 - 0.05f, o2::constants::physics::MassLambda0 + 0.05f}, "binning for the Lambda invariant-mass"};
  ConfigurableAxis momAxis{"momAxisFine", {5.e2, 0.f, 5.f}, "momentum axis binning"};

  // Production selections: loose, they can be tightened when the derived tables are read.
  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.9f, "maximum eta"};
  Configurable<float> minNCrossedRows{"minNCrossedRows", 100, "Minimum number of crossed TPC rows"};
  Configurable<float> minNITSClus{"minNITSClus", 0., "Minimum number of ITS clusters"};
  Configurable<float> minNTPCClus{"minNTPCClus", 80, "Minimum number of TPC clusters"};
  Configurable<float> maxNSharedTPCClus{"maxNSharedTPCClus", 5, "Maximum number of shared TPC clusters"};
  Configurable<float> nSigmaTPCCut{"nSigmaTPCCut", 3.f, "Number of sigmas for the TPC PID"};
  Configurable<float> dcaBachToPV{"dcaBachToPV", 0.05f, "Minimum DCA of a cascade bachelor to the primary vertex"};
  Configurable<float> dcaV0DauToPV{"dcaV0DauToPV", 0.05f, "Minimum DCA of Lambda daughters to the primary vertex"};
  Configurable<float> dcaV0Bach{"dcaV0Bach", 1.f, "Maximum DCA between the V0 and cascade bachelor"};
  Configurable<float> dcaLambdaDaughters{"dcaLambdaDaughters", 1.f, "Maximum DCA between Lambda daughters"};
  Configurable<float> mOmegaWindow{"mOmegaWindow", 0.01f, "Omega mass window"};
  Configurable<float> mLambdaWindow{"mLambdaWindow", 0.01f, "Lambda mass window"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    auto& selections = helper.selections;
    selections.etaMax = etaMax;
    selections.minNITSClus = minNITSClus;
    selections.minNTPCClus = minNTPCClus;
    selections.minNCrossedRows = minNCrossedRows;
    selections.maxNSharedTPCClus = maxNSharedTPCClus;
    selections.nSigmaTPCCut = nSigmaTPCCut;
    selections.dcaBachToPV = dcaBachToPV;
    selections.dcaV0DauToPV = dcaV0DauToPV;
    selections.dcaV0Bach = dcaV0Bach;
    selections.dcaLambdaDaughters = dcaLambdaDaughters;
    selections.mOmegaWindow = mOmegaWindow;
    selections.mLambdaWindow = mLambdaWindow;
    helper.initFitters();
    helper.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value));

    zorroSummary.setObject(zorro.getZorroSummary());

    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    histos.add<TH2>("QA/massXi", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Lambda + #pi) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massXiAxis});
    histos.add<TH2>("QA/massOmega", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Lambda + K) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massOmegaAxis});
    histos.add<TH2>("QA/massLambda", ";#it{p}_{T} (GeV/#it{c});#it{M}(p + #pi) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massLambdaAxis});
    histos.add("QA/collisions", ";;Collisions", HistType::kTH1F, {{kNCollisionSteps, -0.5f, static_cast<float>(kNCollisionSteps) - 0.5f}});
    auto collisions = histos.get<TH1>(HIST("QA/collisions"));
    collisions->GetXaxis()->SetBinLabel(kAllCollisions + 1, "All");
    collisions->GetXaxis()->SetBinLabel(kSelectedCollisions + 1, "Event selection");
    collisions->GetXaxis()->SetBinLabel(kWithOmega + 1, "With #Omega");
    collisions->GetXaxis()->SetBinLabel(kStoredCollisions + 1, "With #Omega, #Lambda and K (stored)");
  }

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
    if (!grpmag) {
      LOG(fatal) << "Could not retrieve GRPMagField for timestamp " << bc.timestamp();
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);
    mBz = o2::base::Propagator::Instance()->getNominalBz();
    helper.setBz(mBz);

    if (static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>("GLO/Param/MatLUT", bc.timestamp()));
      if (!lut) {
        LOG(fatal) << "Could not retrieve material LUT for timestamp " << bc.timestamp();
      }
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    mRunNumber = bc.runNumber();
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fDoubleOmega,fOmegaXi");
      zorro.populateHistRegistry(histos, bc.runNumber());
    }
  }

  template <class C>
  bool acceptCollision(C const& collision)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (cfgApplyEventSelection &&
        (!collision.sel8() ||
         std::abs(collision.posZ()) > zVtxMax ||
         !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) ||
         !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
      return false;
    }
    if (cfgSkimmedProcessing) {
      zorro.isSelected(bc.globalBC());
    }
    histos.fill(HIST("QA/zVtx"), collision.posZ());
    return true;
  }

  void fillCandidateQA(CollisionCandidates const& candidates)
  {
    for (const auto& omega : candidates.omegas) {
      const float omegaPt = std::hypot(omega.px, omega.py);
      histos.fill(HIST("QA/massXi"), omegaPt, omega.massXi);
      histos.fill(HIST("QA/massOmega"), omegaPt, omega.massOmega);
    }
    for (const auto& lambdas : candidates.lambdas) {
      for (const auto& lambda : lambdas) {
        histos.fill(HIST("QA/massLambda"), std::hypot(lambda.candidate.px, lambda.candidate.py), lambda.candidate.mass);
      }
    }
  }

  // Writes the candidates of a collision. Tracks are written once per collision: the direct-kaon
  // candidates first, then the Omega and Lambda daughters not written yet.
  template <class C, class T>
  void writeCollision(C const& collision, CollisionCandidates const& candidates, T const& tracks)
  {
    dOmegaColls(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), mRunNumber, mBz);
    const auto collisionIndex = static_cast<int>(dOmegaColls.lastIndex());

    std::unordered_map<int64_t, int> trackRows; // Track index -> DOmegaTracks row.
    const auto writeTrack = [&](int64_t trackId, bool isKaon) {
      const auto [row, inserted] = trackRows.try_emplace(trackId, -1);
      if (!inserted) {
        return row->second;
      }
      const auto track = tracks.rawIteratorAt(trackId);
      float dcaXY = -999.f;
      helper.dcaXYToPV(collision, track, dcaXY);
      const auto trackParCov = getTrackParCov(track);
      const auto& cov = trackParCov.getCov();
      dOmegaTracks(collisionIndex,
                   trackParCov.getX(), trackParCov.getAlpha(), trackParCov.getY(), trackParCov.getZ(), trackParCov.getSnp(), trackParCov.getTgl(), trackParCov.getQ2Pt(),
                   cov[0], cov[1], cov[2], cov[3], cov[4], cov[5], cov[6], cov[7], cov[8], cov[9], cov[10], cov[11], cov[12], cov[13], cov[14],
                   track.hasTPC() ? track.tpcNSigmaKa() : -999.f,
                   track.hasTOF() ? track.tofNSigmaKa() : -999.f,
                   static_cast<uint8_t>(DoubleOmegaBuilderHelper::detectorMap(track)),
                   dcaXY,
                   isKaon);
      row->second = static_cast<int>(dOmegaTracks.lastIndex());
      return row->second;
    };

    for (const auto& kaons : candidates.kaons) {
      for (const auto& kaon : kaons) {
        writeTrack(kaon.trackId, true);
      }
    }
    for (const auto& omega : candidates.omegas) {
      const auto posTrackRow = writeTrack(omega.posTrackId, false);
      const auto negTrackRow = writeTrack(omega.negTrackId, false);
      const auto bachelorRow = writeTrack(omega.bachelorId, false);
      const auto& parent = omega.parentTrack;
      const auto& cov = parent.getCov();
      dOmegaOmegas(collisionIndex, posTrackRow, negTrackRow, bachelorRow,
                   parent.getX(), parent.getAlpha(), parent.getY(), parent.getZ(), parent.getSnp(), parent.getTgl(), parent.getQ2Pt(),
                   cov[0], cov[1], cov[2], cov[3], cov[4], cov[5], cov[6], cov[7], cov[8], cov[9], cov[10], cov[11], cov[12], cov[13], cov[14],
                   omega.x, omega.y, omega.z,
                   omega.px, omega.py, omega.pz,
                   omega.massOmega, omega.massXi,
                   omega.dcaDaughters,
                   omega.sign);
    }
    for (const int8_t charge : {int8_t{-1}, int8_t{1}}) {
      for (const auto& lambda : candidates.lambdas[DoubleOmegaBuilderHelper::chargeIndex(charge)]) {
        const auto posTrackRow = writeTrack(lambda.posTrackId, false);
        const auto negTrackRow = writeTrack(lambda.negTrackId, false);
        const auto& parent = lambda.parentTrack;
        const auto& cov = parent.getCov();
        dOmegaLambdas(collisionIndex, posTrackRow, negTrackRow,
                      parent.getX(), parent.getAlpha(), parent.getY(), parent.getZ(), parent.getSnp(), parent.getTgl(), parent.getQ2Pt(),
                      cov[0], cov[1], cov[2], cov[3], cov[4], cov[5], cov[6], cov[7], cov[8], cov[9], cov[10], cov[11], cov[12], cov[13], cov[14],
                      lambda.decayVertex[0], lambda.decayVertex[1], lambda.decayVertex[2],
                      lambda.candidate.px, lambda.candidate.py, lambda.candidate.pz,
                      lambda.candidate.mass,
                      lambda.dcaDaughters,
                      charge);
      }
    }
  }

  void processData(CollisionsTable const& collisions,
                   TracksFull const& tracks,
                   aod::V0s const& v0s,
                   aod::Cascades const& cascades,
                   aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      histos.fill(HIST("QA/collisions"), kAllCollisions);
      if (!acceptCollision(collision)) {
        continue;
      }
      histos.fill(HIST("QA/collisions"), kSelectedCollisions);
      auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      auto cascadesThisCollision = cascades.sliceBy(cascadesPerCollision, collision.globalIndex());
      auto v0sThisCollision = v0s.sliceBy(v0sPerCollision, collision.globalIndex());
      v0sThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&v0s);
      const auto candidates = helper.buildCollisionCandidates(collision, tracksThisCollision, cascadesThisCollision, v0sThisCollision, true);
      fillCandidateQA(candidates);
      if (candidates.omegas.empty()) {
        continue;
      }
      histos.fill(HIST("QA/collisions"), kWithOmega);
      if (!DoubleOmegaBuilderHelper::hasOmegaLambdaKaon(candidates)) {
        continue;
      }
      histos.fill(HIST("QA/collisions"), kStoredCollisions);
      writeCollision(collision, candidates, tracks);
    }
  }
  PROCESS_SWITCH(doubleOmegaBuilder, processData, "Build the Omega, Lambda and kaon candidates and store them in derived tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<doubleOmegaBuilder>(cfgc)};
}
