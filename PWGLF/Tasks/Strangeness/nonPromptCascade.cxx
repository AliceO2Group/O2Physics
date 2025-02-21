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

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "Math/Vector4D.h"

#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
// #include "PWGHF/Core/PDG.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "PWGLF/DataModel/LFNonPromptCascadeTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
struct NPCascCandidate {
  int64_t mcParticleId;
  int64_t trackGlobID;
  int64_t trackITSID;
  int64_t collisionID;
  float matchingChi2;
  float deltaPt;
  float itsClusSize;
  bool hasReassociatedCluster;
  bool isGoodMatch;
  bool isGoodCascade;
  int pdgCodeMom;
  int pdgCodeITStrack;
  bool isFromBeauty;
  bool isFromCharm;
  uint16_t pvContributors;
  float pvTimeResolution;
  float pvX;
  float pvY;
  float pvZ;
  float cascPt;
  float cascEta;
  float cascPhi;
  float protonPt;
  float protonEta;
  float pionPt;
  float pionEta;
  float bachPt;
  float bachEta;
  float cascDCAxy;
  float cascDCAz;
  float protonDCAxy;
  float protonDCAz;
  float pionDCAxy;
  float pionDCAz;
  float bachDCAxy;
  float bachDCAz;
  float casccosPA;
  float v0cosPA;
  double massXi;
  double massOmega;
  double massV0;
  float cascRadius;
  float v0radius;
  float cascLength;
  float v0length;
  int cascNClusITS;
  int protonNClusITS;
  int pionNClusITS;
  int bachNClusITS;
  int protonNClusTPC;
  int pionNClusTPC;
  int bachNClusTPC;
  float protonTPCNSigma;
  float pionTPCNSigma;
  float bachKaonTPCNSigma;
  float bachPionTPCNSigma;
  bool protonHasTOF;
  bool pionHasTOF;
  bool bachHasTOF;
  float protonTOFNSigma;
  float pionTOFNSigma;
  float bachKaonTOFNSigma;
  float bachPionTOFNSigma;
};

std::array<bool, 2> isFromHF(auto& particle)
{
  bool fromBeauty = false;
  bool fromCharm = false;
  if (particle.has_mothers()) {
    auto mom = particle.template mothers_as<aod::McParticles>()[0];
    int pdgCodeMom = mom.pdgCode();
    fromBeauty = std::abs(pdgCodeMom) / 5000 == 1 || std::abs(pdgCodeMom) / 500 == 1 || std::abs(pdgCodeMom) == 5;
    fromCharm = std::abs(pdgCodeMom) / 4000 == 1 || std::abs(pdgCodeMom) / 400 == 1 || std::abs(pdgCodeMom) == 4;
    while (mom.has_mothers()) {
      const auto grandma = mom.template mothers_as<aod::McParticles>()[0];
      int pdgCodeGrandma = std::abs(grandma.pdgCode());
      fromBeauty = fromBeauty || (pdgCodeGrandma / 5000 == 1 || pdgCodeGrandma / 500 == 1 || pdgCodeGrandma == 5);
      fromCharm = fromCharm || (pdgCodeGrandma / 4000 == 1 || pdgCodeGrandma / 400 == 1 || pdgCodeGrandma == 4);
      mom = grandma;
    }
  }
  return {fromBeauty, fromCharm};
}

static constexpr int nParticles{4};
static constexpr int nCutsPID{2};
static const std::vector<std::string> matterOrNot{"Matter", "Antimatter"};
static const std::vector<std::string> particlesNames{"K-bachelor", "Pi-bachelor", "Pr", "Pi"};
static const std::vector<std::string> cutsNames{"TPCnSigmaMin", "TPCnSigmaMax"};
static constexpr float cutsPID[nParticles][nCutsPID]{
  {-4.f, +4.f}, /*K bachelor*/
  {-4.f, +4.f}, /*Pi bachelor*/
  {-4.f, +4.f}, /*Pr*/
  {-4.f, +4.f}, /*Pi*/
};

std::vector<NPCascCandidate> gCandidates;
std::vector<NPCascCandidate> gCandidatesNT;

} // namespace

struct NonPromptCascadeTask {

  Produces<o2::aod::NPCascTable> NPCTable;
  Produces<o2::aod::NPCascTableMC> NPCTableMC;
  Produces<o2::aod::NPCascTableNT> NPCTableNT;
  Produces<o2::aod::NPCascTableMCNT> NPCTableMCNT;
  Produces<o2::aod::NPCascTableGen> NPCTableGen;

  using TracksExtData = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullKa, aod::pidTOFFullPi, aod::pidTOFFullPr>;
  using TracksExtMC = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullKa, aod::pidTOFFullPi, aod::pidTOFFullPr>;
  using CollisionCandidatesRun3 = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionCandidatesRun3MC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> cfgPropToPCA{"cfgPropToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> cfgUseAbsDCA{"cfgUseAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<double> cfgMaxR{"cfgMaxR", 200., "reject PCA's above this radius"};
  Configurable<double> cfgMaxDZIni{"cfgMaxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> cfgMinParamChange{"cfgMinParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> cfgMinRelChi2Change{"cfgMinRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> cfgSelectOnlyOmegas{"cfgSelectOnlyOmegas", false, "Toggle to select only Omegas"};

  Configurable<int> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};
  Configurable<float> cfgMinCosPA{"cfgMinCosPA", -1.f, "Minimum cosine of pointing angle"};
  Configurable<LabeledArray<float>> cfgCutsPID{"particlesCutsPID", {cutsPID[0], nParticles, nCutsPID, particlesNames, cutsNames}, "Nuclei PID selections"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", true, "Skimmed dataset processing"};

  Zorro mZorro;
  OutputObj<ZorroSummary> mZorroSummary{"ZorroSummary"};

  Service<o2::ccdb::BasicCCDBManager> mCCDB;
  int mRunNumber = 0;
  float mBz = 0.f;
  o2::vertexing::DCAFitterN<2> mDCAFitter;

  HistogramRegistry mRegistry;

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();

    if (o2::parameters::GRPMagField* grpmag = mCCDB->getForRun<o2::parameters::GRPMagField>(cfgGRPmagPath, mRunNumber)) {
      o2::base::Propagator::initFieldFromGRP(grpmag);
      mBz = static_cast<float>(grpmag->getNominalL3Field());
    }
    mDCAFitter.setBz(mBz);

    if (static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(mCCDB->getForRun<o2::base::MatLayerCylSet>("GLO/Param/MatLUT", mRunNumber));
      o2::base::Propagator::Instance(true)->setMatLUT(lut);
    }
  }

  void init(InitContext const&)
  {
    mZorroSummary.setObject(mZorro.getZorroSummary());
    mCCDB->setURL(ccdbUrl);
    mCCDB->setFatalWhenNull(true);
    mCCDB->setCaching(true);
    mCCDB->setLocalObjectValidityChecking();

    mDCAFitter.setPropagateToPCA(cfgPropToPCA);
    mDCAFitter.setMaxR(cfgMaxR);
    mDCAFitter.setMaxDZIni(cfgMaxDZIni);
    mDCAFitter.setMinParamChange(cfgMinParamChange);
    mDCAFitter.setMinRelChi2Change(cfgMinRelChi2Change);
    mDCAFitter.setUseAbsDCA(cfgUseAbsDCA);

    std::vector<double> ptBinning = {0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    std::array<std::string, 7> cutsNames{"# candidates", "hasTOF", "nClusTPC", "nSigmaTPCbach", "nSigmaTPCprotontrack", "nSigmaTPCpiontrack", "cosPA"};
    auto cutsOmega{std::get<std::shared_ptr<TH2>>(mRegistry.add("h_PIDcutsOmega", ";;Invariant mass (GeV/#it{c}^{2})", HistType::kTH2D, {{cutsNames.size(), -0.5, -0.5 + cutsNames.size()}, {125, 1.650, 1.700}}))};
    auto cutsXi{std::get<std::shared_ptr<TH2>>(mRegistry.add("h_PIDcutsXi", ";;Invariant mass (GeV/#it{c}^{2})", HistType::kTH2D, {{6, -0.5, 5.5}, {125, 1.296, 1.346}}))};

    for (size_t iBin{0}; iBin < cutsNames.size(); ++iBin) {
      cutsOmega->GetYaxis()->SetBinLabel(iBin + 1, cutsNames[iBin].c_str());
      cutsXi->GetYaxis()->SetBinLabel(iBin + 1, cutsNames[iBin].c_str());
    }
  }

  void zorroAccounting(const auto& collisions)
  {
    if (cfgSkimmedProcessing) {
      int runNumber{-1};
      for (const auto& coll : collisions) {
        auto bc = coll.template bc_as<aod::BCsWithTimestamps>();
        if (runNumber != bc.runNumber()) {
          mZorro.initCCDB(mCCDB.service, bc.runNumber(), bc.timestamp(), "fTrackedOmega");
          mZorro.populateHistRegistry(mRegistry, bc.runNumber());
          runNumber = bc.runNumber();
        }
        mZorro.isSelected(bc.globalBC()); /// Just let Zorro do the accounting
      }
    }
  }

  template <typename TrackType, typename CollisionType>
  void fillCandidatesVector(CollisionType const&, auto const& cascades, auto& candidates)
  {

    const auto& getCascade = [](auto const& candidate) {
      if constexpr (requires { candidate.cascade(); }) {
        return candidate.cascade();
      } else {
        return candidate;
      }
    };

    candidates.clear();
    for (const auto& candidate : cascades) {

      auto collision = candidate.template collision_as<CollisionType>();
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      const auto primaryVertex = getPrimaryVertex(collision);

      const auto& casc = getCascade(candidate);
      const auto& bachelor = casc.template bachelor_as<TrackType>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.template posTrack_as<TrackType>();
      const auto& ntrack = v0.template negTrack_as<TrackType>();
      const auto& protonTrack = bachelor.sign() > 0 ? ntrack : ptrack;
      const auto& pionTrack = bachelor.sign() > 0 ? ptrack : ntrack;

      mRegistry.fill(HIST("h_PIDcutsXi"), 0, 1.322);
      mRegistry.fill(HIST("h_PIDcutsOmega"), 0, 1.675);

      mRegistry.fill(HIST("h_PIDcutsXi"), 1, 1.322);
      mRegistry.fill(HIST("h_PIDcutsOmega"), 1, 1.675);

      if (protonTrack.tpcNClsFound() < cfgCutNclusTPC || pionTrack.tpcNClsFound() < cfgCutNclusTPC || bachelor.tpcNClsFound() < cfgCutNclusTPC) {
        continue;
      }
      mRegistry.fill(HIST("h_PIDcutsXi"), 2, 1.322);
      mRegistry.fill(HIST("h_PIDcutsOmega"), 2, 1.675);

      // QA PID
      float nSigmaTPC[nParticles]{bachelor.tpcNSigmaKa(), bachelor.tpcNSigmaPi(), protonTrack.tpcNSigmaPr(), pionTrack.tpcNSigmaPi()};

      bool isBachelorSurvived = false;
      if (nSigmaTPC[0] > cfgCutsPID->get(0u, 0u) && nSigmaTPC[0] < cfgCutsPID->get(0u, 1u)) {
        mRegistry.fill(HIST("h_PIDcutsOmega"), 3, 1.675);
        isBachelorSurvived = true;
      }

      if (!cfgSelectOnlyOmegas && nSigmaTPC[1] > cfgCutsPID->get(1u, 0u) && nSigmaTPC[1] < cfgCutsPID->get(1u, 1u)) {
        mRegistry.fill(HIST("h_PIDcutsXi"), 3, 1.322);
        isBachelorSurvived = true;
      }

      if (!isBachelorSurvived) {
        continue;
      }

      if (nSigmaTPC[2] < cfgCutsPID->get(2u, 0u) || nSigmaTPC[2] > cfgCutsPID->get(2u, 1u)) {
        continue;
      }

      mRegistry.fill(HIST("h_PIDcutsOmega"), 4, 1.675);
      mRegistry.fill(HIST("h_PIDcutsXi"), 4, 1.322);

      if (nSigmaTPC[3] < cfgCutsPID->get(3u, 0u) || nSigmaTPC[3] > cfgCutsPID->get(3u, 1u)) {
        continue;
      }

      mRegistry.fill(HIST("h_PIDcutsOmega"), 5, 1.675);
      mRegistry.fill(HIST("h_PIDcutsXi"), 5, 1.322);

      auto protonTrkParCov = getTrackParCov(protonTrack);
      auto pionTrkParCov = getTrackParCov(pionTrack);
      auto bachTrkParCov = getTrackParCov(bachelor);

      std::array<std::array<float, 3>, 2> momenta;
      std::array<float, 3> cascadeMomentum;
      o2::math_utils::SVector<double, 3> cascadePos, v0Pos;

      float cascCpa = -1, v0Cpa = -1;
      if (mDCAFitter.process(pionTrkParCov, protonTrkParCov)) {
        auto trackParCovV0 = mDCAFitter.createParentTrackParCov(0); // V0 track retrieved from p and pi daughters
        v0Pos = mDCAFitter.getPCACandidate();
        if (mDCAFitter.process(trackParCovV0, bachTrkParCov)) {
          mDCAFitter.getTrackParamAtPCA(0).getPxPyPzGlo(momenta[0]);
          mDCAFitter.getTrackParamAtPCA(1).getPxPyPzGlo(momenta[1]);
          mDCAFitter.createParentTrackParCov().getPxPyPzGlo(cascadeMomentum);
          std::array<float, 3> pvPos = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};
          cascadePos = mDCAFitter.getPCACandidate();
          cascCpa = RecoDecay::cpa(pvPos, mDCAFitter.getPCACandidate(), cascadeMomentum);
          v0Cpa = RecoDecay::cpa(pvPos, v0Pos, momenta[0]);
        } else {
          continue;
        }
      } else {
        continue;
      }
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> cascadeLvector;
      cascadeLvector.SetPxPyPzE(cascadeMomentum[0], cascadeMomentum[1], cascadeMomentum[2], std::hypot(cascadeMomentum[0], cascadeMomentum[1], cascadeMomentum[2])); /// 0 mass, used only for the momentum

      // Omega
      std::array<double, 2> masses{o2::constants::physics::MassLambda0, o2::constants::physics::MassKPlus};
      const auto massOmega = RecoDecay::m(momenta, masses);

      // Xi
      masses = {o2::constants::physics::MassLambda0, o2::constants::physics::MassPiPlus};
      const auto massXi = RecoDecay::m(momenta, masses);

      // Lambda
      masses = {o2::constants::physics::MassProton, o2::constants::physics::MassPiMinus};
      momenta[0] = {protonTrack.px(), protonTrack.py(), protonTrack.pz()};
      momenta[1] = {pionTrack.px(), pionTrack.py(), pionTrack.pz()};
      const auto v0mass = RecoDecay::m(momenta, masses);

      //// Omega hypohesis -> rejecting Xi, we don't do it in the MC as we can identify the particle with the MC truth
      bool isOmega{std::abs(massXi - constants::physics::MassXiMinus) > 0.005};
      if (cfgSelectOnlyOmegas && !isOmega) {
        continue;
      }

      std::array<bool, 2> fromHF{false, false};
      bool isGoodMatch{false}, isGoodCascade{false};
      int itsTrackPDG{0}, pdgCodeMom{0};
      int64_t mcParticleID{-1};

      if constexpr (TrackType::template contains<aod::McTrackLabels>()) {
        if (protonTrack.mcParticle().has_mothers() && pionTrack.mcParticle().has_mothers() && bachelor.mcParticle().has_mothers()) {
          if (protonTrack.mcParticle().mothersIds()[0] == pionTrack.mcParticle().mothersIds()[0]) {
            const auto v0part = protonTrack.mcParticle().template mothers_first_as<aod::McParticles>();
            if (std::abs(v0part.pdgCode()) == 3122 && v0part.has_mothers()) {
              const auto motherV0 = v0part.template mothers_as<aod::McParticles>()[0];
              if (v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
                if (std::abs(motherV0.pdgCode()) == 3312 || std::abs(motherV0.pdgCode()) == 3334) {
                  isGoodCascade = true;
                  isOmega = (std::abs(motherV0.pdgCode()) == 3334);
                  fromHF = isFromHF(motherV0);
                  mcParticleID = v0part.mothersIds()[0];
                }
              }
            }
          }
        }
      }

      if (cascCpa < cfgMinCosPA) {
        continue;
      }
      if (isOmega) {
        mRegistry.fill(HIST("h_PIDcutsOmega"), 6, massOmega);
      }
      mRegistry.fill(HIST("h_PIDcutsXi"), 6, massXi);

      const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);
      o2::dataformats::DCA motherDCA{-999.f, -999.f}, protonDCA{-999.f, -999.f}, pionDCA{-999.f, -999.f}, bachDCA{-999.f, -999.f};
      o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, protonTrkParCov, mBz, 2.f, matCorr, &protonDCA);
      o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, pionTrkParCov, mBz, 2.f, matCorr, &pionDCA);
      o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, bachTrkParCov, mBz, 2.f, matCorr, &bachDCA);

      float deltaPtITSCascade{-1.e10f}, cascITSclsSize{-1.e10f}, matchingChi2{-1.e10f};
      bool hasReassociatedClusters{false};
      int trackedCascGlobalIndex{-1}, itsTrackGlobalIndex{-1}, cascITSclusters{-1};
      if constexpr (requires { candidate.track(); }) {
        const auto& track = candidate.template track_as<TrackType>();
        const auto& ITStrack = candidate.template itsTrack_as<TrackType>();
        auto trackTrkParCov = getTrackParCov(track);
        o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackTrkParCov, mBz, 2.f, matCorr, &motherDCA);
        hasReassociatedClusters = (track.itsNCls() != ITStrack.itsNCls());
        cascadeLvector.SetCoordinates(track.pt(), track.eta(), track.phi(), 0);
        deltaPtITSCascade = std::hypot(cascadeMomentum[0], cascadeMomentum[1]) - ITStrack.pt();
        trackedCascGlobalIndex = track.globalIndex();
        itsTrackGlobalIndex = ITStrack.globalIndex();
        cascITSclusters = track.itsNCls();
        cascITSclsSize = candidate.itsClsSize();
        matchingChi2 = candidate.matchingChi2();
        cascadePos = {candidate.decayX(), candidate.decayY(), candidate.decayZ()};
        if constexpr (TrackType::template contains<aod::McTrackLabels>()) {
          isGoodMatch = ((mcParticleID == ITStrack.mcParticleId())) ? true : false;

          if (isGoodMatch) {
            pdgCodeMom = track.mcParticle().has_mothers() ? track.mcParticle().template mothers_as<aod::McParticles>()[0].pdgCode() : 0;
          }
          itsTrackPDG = ITStrack.has_mcParticle() ? ITStrack.mcParticle().pdgCode() : 0;
        }
      }
      candidates.emplace_back(NPCascCandidate{mcParticleID, trackedCascGlobalIndex, itsTrackGlobalIndex, candidate.collisionId(), matchingChi2, deltaPtITSCascade, cascITSclsSize, hasReassociatedClusters, isGoodMatch, isGoodCascade, pdgCodeMom, itsTrackPDG, fromHF[0], fromHF[1],
                                              collision.numContrib(), collision.collisionTimeRes(), primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                                              cascadeLvector.pt(), cascadeLvector.eta(), cascadeLvector.phi(),
                                              protonTrack.pt(), protonTrack.eta(), pionTrack.pt(), pionTrack.eta(), bachelor.pt(), bachelor.eta(),
                                              motherDCA.getY(), motherDCA.getZ(), protonDCA.getY(), protonDCA.getZ(), pionDCA.getY(), pionDCA.getZ(), bachDCA.getY(), bachDCA.getZ(),
                                              cascCpa, v0Cpa, massXi, massOmega, v0mass,
                                              static_cast<float>(std::hypot(cascadePos[0], cascadePos[1])), static_cast<float>(std::hypot(v0Pos[0], v0Pos[1])), static_cast<float>(std::hypot(cascadePos[0], cascadePos[1], cascadePos[2])), static_cast<float>(std::hypot(v0Pos[0], v0Pos[1], v0Pos[2])),
                                              cascITSclusters, protonTrack.itsNCls(), pionTrack.itsNCls(), bachelor.itsNCls(), protonTrack.tpcNClsFound(), pionTrack.tpcNClsFound(), bachelor.tpcNClsFound(),
                                              protonTrack.tpcNSigmaPr(), pionTrack.tpcNSigmaPi(), bachelor.tpcNSigmaKa(), bachelor.tpcNSigmaPi(),
                                              protonTrack.hasTOF(), pionTrack.hasTOF(), bachelor.hasTOF(),
                                              protonTrack.tofNSigmaPr(), pionTrack.tofNSigmaPi(), bachelor.tofNSigmaKa(), bachelor.tofNSigmaPi()});
    }
  }

  template <typename CascadeType>
  void fillDataTable(auto const& candidates)
  {
    for (const auto& c : candidates) {
      getDataTable<CascadeType>()(c.matchingChi2, c.deltaPt, c.itsClusSize, c.hasReassociatedCluster,
                                  c.pvContributors, c.pvTimeResolution, c.pvX, c.pvY, c.pvZ,
                                  c.cascPt, c.cascEta, c.cascPhi,
                                  c.protonPt, c.protonEta, c.pionPt, c.pionEta, c.bachPt, c.bachEta,
                                  c.cascDCAxy, c.cascDCAz, c.protonDCAxy, c.protonDCAz, c.pionDCAxy, c.pionDCAz, c.bachDCAxy, c.bachDCAz,
                                  c.casccosPA, c.v0cosPA,
                                  c.massXi, c.massOmega, c.massV0,
                                  c.cascRadius, c.v0radius, c.cascLength, c.v0length,
                                  c.cascNClusITS, c.protonNClusITS, c.pionNClusITS, c.bachNClusITS, c.protonNClusTPC, c.pionNClusTPC, c.bachNClusTPC,
                                  c.protonTPCNSigma, c.pionTPCNSigma, c.bachKaonTPCNSigma, c.bachPionTPCNSigma,
                                  c.protonHasTOF, c.pionHasTOF, c.bachHasTOF,
                                  c.protonTOFNSigma, c.pionTOFNSigma, c.bachKaonTOFNSigma, c.bachPionTOFNSigma);
    }
  }

  template <typename CascadeType>
  void fillMCtable(auto const& mcParticles, auto const& collisions, auto const& candidates)
  {
    for (size_t i = 0; i < candidates.size(); ++i) {
      auto& c = candidates[i];
      if (c.mcParticleId < 0) {
        continue;
      }
      auto particle = mcParticles.iteratorAt(c.mcParticleId);
      auto mcCollision = particle.template mcCollision_as<aod::McCollisions>();
      auto recCollision = collisions.iteratorAt(c.collisionID);

      getMCtable<CascadeType>()(c.matchingChi2, c.deltaPt, c.itsClusSize, c.hasReassociatedCluster, c.isGoodMatch, c.isGoodCascade, c.pdgCodeMom, c.pdgCodeITStrack, c.isFromBeauty, c.isFromCharm,
                                c.pvContributors, c.pvTimeResolution, c.pvX, c.pvY, c.pvZ, c.cascPt, c.cascEta, c.cascPhi,
                                c.protonPt, c.protonEta, c.pionPt, c.pionEta, c.bachPt, c.bachEta,
                                c.cascDCAxy, c.cascDCAz, c.protonDCAxy, c.protonDCAz, c.pionDCAxy, c.pionDCAz, c.bachDCAxy, c.bachDCAz,
                                c.casccosPA, c.v0cosPA, c.massXi, c.massOmega, c.massV0, c.cascRadius, c.v0radius, c.cascLength, c.v0length,
                                c.cascNClusITS, c.protonNClusITS, c.pionNClusITS, c.bachNClusITS, c.protonNClusTPC, c.pionNClusTPC, c.bachNClusTPC, c.protonTPCNSigma,
                                c.pionTPCNSigma, c.bachKaonTPCNSigma, c.bachPionTPCNSigma, c.protonHasTOF, c.pionHasTOF, c.bachHasTOF,
                                c.protonTOFNSigma, c.pionTOFNSigma, c.bachKaonTOFNSigma, c.bachPionTOFNSigma,
                                particle.pt(), particle.eta(), particle.phi(), particle.pdgCode(), mcCollision.posX() - particle.vx(), mcCollision.posY() - particle.vy(),
                                mcCollision.posZ() - particle.vz(), mcCollision.globalIndex() == recCollision.mcCollisionId());
    }
  }

  template <typename CascadeType>
  auto& getMCtable()
  {
    if constexpr (std::is_same_v<CascadeType, aod::Cascades>) {
      return NPCTableMCNT;
    } else {
      return NPCTableMC;
    }
  }

  template <typename CascadeType>
  auto& getDataTable()
  {
    if constexpr (std::is_same_v<CascadeType, aod::Cascades>) {
      return NPCTableNT;
    } else {
      return NPCTable;
    }
  }

  void processTrackedCascadesMC(CollisionCandidatesRun3MC const& collisions,
                                aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& /*cascades*/,
                                aod::V0s const& /*v0s*/, TracksExtMC const& /*tracks*/,
                                aod::McParticles const& mcParticles, aod::McCollisions const&, aod::BCsWithTimestamps const&)
  {
    fillCandidatesVector<TracksExtMC>(collisions, trackedCascades, gCandidates);
    fillMCtable<aod::AssignedTrackedCascades>(mcParticles, collisions, gCandidates);
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processTrackedCascadesMC, "process cascades from strangeness tracking: MC analysis", true);

  void processCascadesMC(CollisionCandidatesRun3MC const& collisions, aod::Cascades const& cascades,
                         aod::V0s const& /*v0s*/, TracksExtMC const& /*tracks*/,
                         aod::McParticles const& mcParticles, aod::McCollisions const&, aod::BCsWithTimestamps const&)
  {
    fillCandidatesVector<TracksExtMC>(collisions, cascades, gCandidatesNT);
    fillMCtable<aod::Cascades>(mcParticles, collisions, gCandidatesNT);
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processCascadesMC, "process cascades: MC analysis", false);

  void processGenParticles(aod::McParticles const& mcParticles, aod::McCollisions const&)
  {
    for (const auto& p : mcParticles) {
      auto absCode = std::abs(p.pdgCode());
      if (absCode != 3312 && absCode != 3334) {
        continue;
      }
      auto fromHF = isFromHF(p);
      int pdgCodeMom = p.has_mothers() ? p.template mothers_as<aod::McParticles>()[0].pdgCode() : 0;
      auto mcCollision = p.template mcCollision_as<aod::McCollisions>();

      NPCTableGen(p.pt(), p.eta(), p.phi(), p.pdgCode(), pdgCodeMom, mcCollision.posX() - p.vx(), mcCollision.posY() - p.vy(), mcCollision.posZ() - p.vz(), fromHF[0], fromHF[1]);
    }
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processGenParticles, "process gen cascades: MC analysis", false);

  void processTrackedCascadesData(CollisionCandidatesRun3 const& collisions,
                                  aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& /*cascades*/,
                                  aod::V0s const& /*v0s*/, TracksExtData const& /*tracks*/,
                                  aod::BCsWithTimestamps const&)
  {
    zorroAccounting(collisions);
    fillCandidatesVector<TracksExtData>(collisions, trackedCascades, gCandidates);
    fillDataTable<aod::AssignedTrackedCascades>(gCandidates);
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processTrackedCascadesData, "process cascades from strangeness tracking: Data analysis", false);

  void processCascadesData(CollisionCandidatesRun3 const& collisions, aod::Cascades const& cascades,
                           aod::V0s const& /*v0s*/, TracksExtData const& /*tracks*/,
                           aod::BCsWithTimestamps const&)
  {
    zorroAccounting(collisions);
    fillCandidatesVector<TracksExtData>(collisions, cascades, gCandidatesNT);
    fillDataTable<aod::Cascades>(gCandidatesNT);
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processCascadesData, "process cascades: Data analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NonPromptCascadeTask>(cfgc)};
}
