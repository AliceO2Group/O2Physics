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
/// \brief kinkAnalysis: Analysis task for kink topology, to be converted later in a table producer for common usage by other tasks
/// \author everyone

#include <cmath>

#include "Math/Vector4D.h"
#include "Common/Core/CollisionAssociation.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"
#include "Framework/DataTypes.h"
#include "DCAFitter/DCAFitterN.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/LHCConstants.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "DetectorsVertexing/SVertexer.h"
#include "DetectorsVertexing/SVertexerParams.h"
#include "TPCReconstruction/TPCFastTransformHelperO2.h"
#include "DataFormatsTPC/VDriftCorrFact.h"
#include "DataFormatsCalibration/MeanVertexObject.h"

#include "ReconstructionDataFormats/Track.h"

#include <iostream>
using namespace std;

using namespace o2;
using namespace vertexing;

using namespace o2::dataformats;
using namespace o2::aod;

using namespace o2::framework;
using namespace o2::framework::expressions;
using VBracket = o2::math_utils::Bracket<int>;

static constexpr int POS = 0, NEG = 1;

namespace kink
{
constexpr std::array<float, 7> LayerRadii{2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
}

struct KinkCandidates {
  int chargeMother;
  int globalIndexMother;
  int globalIndexDaughter;
  int mcParticleIdxMother;
  int mcParticleIdxDaughter;
  float prx, pry, prz;
  float decayvtxX, decayvtxY, decayvtxZ;
  float mthinnerpx, mthinnerpy, mthinnerpz;
  float mthouterpx, mthouterpy, mthouterpz;
  float dghtpx, dghtpy, dghtpz;
  unsigned int flags = 0u;
};

struct kinkAnalysis {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  enum strParticleAdDecay { SigmaMinus,
                            SigmaPlusToPi,
                            SigmaPlusToProton,
                            Kaon,
                            Pion,
                            Xi,
                            OmegaToL,
                            OmegaToXi };

  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<double> cfgBz{"cfgBz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> cfgGRPpath{"cfgGRPpath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> cfgCCDBurl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> cfgLUTpath{"cfgLUTpath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> dcaKinkDtopv{"dcaKinkDtopv", .2, "DCA kink daughter To PV"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<int> cfgMotherCharge{"cfgMotherCharge", -1, "mother charge"};
  Configurable<int> cfgDaughterCharge{"cfgDaughterCharge", -1, "mother charge"};
  Configurable<int> cfgParticleSpec{"cfgParticleSpec", SigmaPlusToProton, "particle species"};
  Configurable<float> cfgLowerHistLimit{"cfgLowerHistLimit", 1.1, "lower limit for inv mass histo"};
  Configurable<float> cfgUpperHistLimit{"cfgUpperHistLimit", 1.4, "upper limit for inv mass histo"};
  Configurable<float> cfgNsigmaTPCdaughter{"cfgNsigmaTPCdaughter", 3, "nSigma TPC for the daughter track"};
  Configurable<float> qTupper{"qTupper", 0.195, "upper qT cut"};
  Configurable<float> qTlower{"qTLower", 0.150, "upper qT cut"};
  Configurable<float> kinkAngle{"kinkAngle", 0.5, "mother-daughter angle cutoff"};
  Configurable<float> zDiff{"zDiff", 20., "mother-daughter z diff"};
  Configurable<float> phiDiff{"phiDiff", 100., "mother-daughter phi diff"};
  Configurable<bool> cfgIsMC{"cfgIsMC", true, "data or MC"};

  o2::base::MatLayerCylSet* lut = nullptr;
  o2::dataformats::MeanVertexObject* mVtx = nullptr;
  o2::vertexing::DCAFitterN<2> ft2;

  std::vector<KinkCandidates> mKinkCandidates;

  int mRunNumber;
  float mBz;

  float etaS = 0;
  float phiS = 0;
  float etaP = 0;
  float phiP = 0;

  std::array<float, 3> sigmaPDC = {0, 0, 0};
  std::array<float, 3> pionPDC = {0, 0, 0};
  std::array<float, 3> sigmaPin = {0, 0, 0};

  float pionPabsDC = 0;
  float sigmaPabsDC = 0;
  float sigmaPt = 0.;

  float costheta;
  float theta;
  float qT;

  float pionE = 0;
  double neutronE = 0;
  float sigmaE = 0;
  float neutronPabs = 0;
  float neutronM = 0;
  float mass;
  float radToDeg = 180. / M_PI;
  int particlePdgCode;

  o2::track::TrackParCov SigmaTr;
  o2::track::TrackParCov SigmaTr2;
  o2::track::TrackParCov PionTr;

  float mMother, mChargedDaughter, mNeutralDaughter;

  using CompleteTracks = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  using CompleteCollisions = soa::Join<aod::Collisions, aod::EvSels>;

  struct TrackCand {
    int Idxtr;
    int mcParticleIdx = -1;
    VBracket vBracket{};
  };

  float angleCutFunction(int particleName, float x)
  {

    float mSigmaMinus = 1.197449, mSigmaPlus = 1.189449, mKaon = 0.493677, mPion = 0.13957018, mXi = 1.321, mOmega = 1.67245;
    float par1 = mSigmaMinus; /// Default value is for SigmaMinus
    float par2 = 0.8;
    float par3 = M_PI;
    switch (particleName) {
      case SigmaPlusToPi:
        par1 = mSigmaPlus;
        par2 = 0.8;
        break;

      case SigmaPlusToProton:
        par1 = mSigmaPlus;
        par2 = 0.2;
        break;

      case Xi:
        par1 = mXi;
        par2 = 0.68;
        break;

      case OmegaToL:
        par1 = mOmega;
        par2 = 0.68;
        break;

      case OmegaToXi:
        par1 = mOmega;
        par2 = 0.68;
        break;

      case Kaon:
        par1 = mKaon;
        par2 = 0.9127037;
        break;

      case Pion:
        par1 = mPion;
        par2 = 0.2731374;
        break;
    }

    if ((particleName == SigmaMinus) || (particleName == SigmaPlusToPi) || (particleName == SigmaPlusToProton) || (particleName == Xi) || (particleName == OmegaToXi) || (particleName == OmegaToL))
      return par1 * (par2 / (sqrt((x * x) * (1 - (par2 * par2)) - ((par1 * par1) * (par2 * par2))))) * (180. / par3) + 1;

    return ((atan(par1 * par2 * (1.0 / (sqrt((x * x) * (1.0 - (par1 * par1)) - (par1 * par1) * (par2 * par2)))))) * 180.) / par3;
  }

  void init(InitContext const&)
  {
    ccdb->setURL(cfgCCDBurl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    int anParticleName = cfgParticleSpec;
    if (cfgIsMC) {
      switch (anParticleName) {
        case SigmaPlusToPi:
          if (cfgMotherCharge == 1)
            particlePdgCode = 3222;
          else
            particlePdgCode = -3222;
          break;

        case SigmaPlusToProton:
          if (cfgMotherCharge == 1)
            particlePdgCode = 3222;
          else
            particlePdgCode = -3222;
          break;

        case Xi:
          if (cfgMotherCharge == -1)
            particlePdgCode = 3312;
          else
            particlePdgCode = -3312;
          break;

        case OmegaToL:
          if (cfgMotherCharge == -1)
            particlePdgCode = 3334;
          else
            particlePdgCode = -3334;
          break;

        case OmegaToXi:
          if (cfgMotherCharge == -1)
            particlePdgCode = 3334;
          else
            particlePdgCode = -3334;
          break;

        case SigmaMinus:
          if (cfgMotherCharge == -1)
            particlePdgCode = 3112;
          else
            particlePdgCode = -3112;
          break;

        default:
          particlePdgCode = 3112;
      }
    }

    mRunNumber = 0;
    mBz = 0.f;

    // define axes you want to use
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisqT{1000, 0.0, 1.0, "q_{T}"};
    const AxisSpec axisRmother{500, 0.0, 50.0, "#it{R}_{mother} (cm)"};
    const AxisSpec axisDCAdaugh{1000, 0.0, 10.0, "#it{DCA}_{daughter} (cm)"};
    const AxisSpec axisSigmaMass{300, cfgLowerHistLimit, cfgUpperHistLimit, "#it{M}_{inv} (Gev/#it{c^{2}})"};
    const AxisSpec axisdtime{10000, 0, 50, "#delta t"};
    const AxisSpec axisPt{500, 0., 50., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisdphi{360, 0, 360, "#delta phi"};
    const AxisSpec axisdz{400, -20., 20, "#delta z"};
    const AxisSpec axisPdgCodes{8001, -4000.5, 4000.5, "mother pdg codes"};

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("hqT", "hqT", kTH1F, {axisqT});
    histos.add("hRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    histos.add("hVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    histos.add("hCollId", "collision id", HistType::kTH1F, {{1000, 0., 100000., "collision id"}});
    histos.add("hRadiusMth", "hRadiusMth", kTH1F, {axisRmother});
    histos.add("hMassMinusPt", "hMassMinusPt", kTH2F, {axisSigmaMass, axisPt});
    histos.add("hMassPlusPt", "hMassPlusPt", kTH2F, {axisSigmaMass, axisPt});
    histos.add("hBackgroundPosNegPt", "hBackgroundPosNegPt", kTH2F, {axisSigmaMass, axisPt});
    histos.add("hBackgroundNegPosPt", "hBackgroundNegPosPt", kTH2F, {axisSigmaMass, axisPt});
    histos.add("hDCAdaughterToPV", "hDCAdaughterToPV", kTH1F, {axisDCAdaugh});
    histos.add("hDCAMotherToPV", "hDCAMotherToPV", kTH1F, {axisDCAdaugh});
    histos.add("hdeltatime", "hdeltatime", kTH1F, {axisdtime});
    histos.add("hdelt_tthresh", "hdelt_tthresh", kTH2F, {axisdtime, axisdtime});
    histos.add("hdeltaphi", "hdeltaphi", kTH1F, {axisdphi});
    histos.add("hdeltaz", "hdeltaz", kTH1F, {axisdz});
    histos.add("generatedPt", "generatedPt", kTH1F, {axisPt});
    histos.add("hPtMinusRecMcTrth", "hPtMinusRecMcTrth", kTH2F, {axisSigmaMass, axisPt});
    histos.add("hPtPlusRecMcTrth", "hPtPlusRecMcTrth", kTH2F, {axisSigmaMass, axisPt});
    histos.add("hcodes", "hcodes", kTH2F, {axisPdgCodes, axisPdgCodes});
    histos.add("hptMtrue", "hptMtrue", kTH2F, {axisPt, axisPt});
    histos.add("hptMDtrue", "hptMDtrue", kTH2F, {axisPt, axisPt});
    histos.add("hptMDelse", "hptMDelse", kTH2F, {axisPt, axisPt});
    histos.add("hPtMinusRecMcTrthM", "hPtMinusRecMcTrthM", kTH2F, {axisSigmaMass, axisPt});
    histos.add("hPtMinusRecMcTrthelse", "hPtMinusRecMcTrthelse", kTH2F, {axisSigmaMass, axisPt});

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(cfgLUTpath));
    ft2.setMaxChi2(5);
    ft2.setUseAbsDCA(true);

    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    ft2.setMatCorrType(matCorr);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (cfgBz > -990) {
      mBz = cfgBz;
      ft2.setBz(mBz);
      o2::parameters::GRPMagField grpmag;
      if (fabs(mBz) > 1e-5) {
        grpmag.setL3Current(30000.f / (mBz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(cfgGRPpath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      mBz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgGRPmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << cfgGRPmagPath << " of object GRPMagField and " << cfgGRPpath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      mBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    }
    mVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    mRunNumber = bc.runNumber();

    ft2.setBz(mBz);

    o2::base::Propagator::Instance()->setMatLUT(lut);
  }

  std::array<std::vector<TrackCand>, 4> makeTracksPool(CompleteCollisions const& collisions, CompleteTracks const& tracks, o2::aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const&)
  {
    mKinkCandidates.clear();
    std::unordered_map<int, std::pair<int, int>> tmap;
    TrackCand trForpool;

    std::array<std::vector<TrackCand>, 4> pools; // pools of positive and negative seeds sorted in min VtxID

    std::vector<uint8_t> selected(tracks.size(), 0u);
    std::vector<uint64_t> globalBCvector;

    int index{0};
    for (const auto& track : tracks) {
      if (track.has_collision()) {
        if (track.collision_as<CompleteCollisions>().has_bc()) {
          globalBCvector.push_back(track.collision_as<CompleteCollisions>().bc_as<aod::BCsWithTimestamps>().globalBC());
        }
      } else {
        for (const auto& ambTrack : ambiTracks) {
          if (ambTrack.trackId() == track.globalIndex()) {
            globalBCvector.push_back(ambTrack.bc().begin().globalBC());
            break;
          }
        }
      }
      if (std::abs(track.eta()) < 0.8) {
        if (track.hasITS() && !track.hasTPC() && !track.hasTOF() && track.itsNCls() < 6 && track.itsNClsInnerBarrel() == 3 && track.itsChi2NCl() < 4) {
          selected[index] = 1;
        } else if (track.hasITS() && track.hasTPC() && track.itsNClsInnerBarrel() == 0 && track.itsNCls() < 4 && track.tpcNClsCrossedRows() >= 70 && track.itsChi2NCl() < 36.f && track.tpcChi2NCl() < 4.f && track.tpcNClsCrossedRows() > 0.8 * track.tpcNClsFindable()) {
          selected[index] = 2;
        }
      }
      index++;
    }

    constexpr auto bOffsetMax = 241; // 6 mus (ITS)

    for (const auto& collision : collisions) {
      if (!collision.sel8())
        continue;
      if (std::abs(collision.posZ()) > 10.)
        continue;

      const float collTime = collision.collisionTime();
      const float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.bc_as<aod::BCsWithTimestamps>().globalBC();
      const auto collIdx = collision.globalIndex();

      histos.fill(HIST("hVtxZData"), collision.posZ());
      index = -1;
      for (const auto& track : tracks) {
        index++;
        if (!selected[index] || !track.has_collision())
          continue;
        const int64_t bcOffset = (int64_t)globalBCvector[track.filteredIndex()] - (int64_t)collBC;
        if (std::abs(bcOffset) > bOffsetMax) {
          continue;
        }

        float trackTime{0.};
        float trackTimeRes{0.};
        if (track.isPVContributor()) {
          trackTime = track.collision_as<CompleteCollisions>().collisionTime(); // if PV contributor, we assume the time to be the one of the collision
          trackTimeRes = constants::lhc::LHCBunchSpacingNS;                     // 1 BC
        } else {
          trackTime = track.trackTime();
          trackTimeRes = track.trackTimeRes();
        }

        const float deltaTime = trackTime - collTime + bcOffset * constants::lhc::LHCBunchSpacingNS;
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;

        float thresholdTime = 0.;
        if (track.isPVContributor()) {
          thresholdTime = trackTimeRes;
        } else if (TESTBIT(track.flags(), o2::aod::track::TrackTimeResIsRange)) {
          thresholdTime = std::sqrt(sigmaTimeRes2); // + timeMargin;
        } else {
          // thresholdTime = nSigmaForTimeCompat * std::sqrt(sigmaTimeRes2); // + timeMargin;
          thresholdTime = 4. * std::sqrt(sigmaTimeRes2); // + timeMargin;
        }

        histos.fill(HIST("hdeltatime"), deltaTime);
        histos.fill(HIST("hdelt_tthresh"), deltaTime, thresholdTime);

        if (std::abs(deltaTime) < thresholdTime) {

          const auto& tref = tmap.find(track.globalIndex());
          if (tref != tmap.end()) {
            pools[tref->second.second][tref->second.first].vBracket.setMax(static_cast<int>(collIdx)); // this track was already processed with other vertex, account the latter
            continue;
          }
        }

        int poolIndex = (selected[index] - 1) * 2 + (track.sign() < 0); /// first the two mothers then the two daughters (mom pos 0, mom neg 1, dau pos 2, dau neg 3)
        trForpool.Idxtr = track.globalIndex();
        trForpool.vBracket = {static_cast<int>(collIdx), static_cast<int>(collIdx)};
        pools[poolIndex].emplace_back(trForpool);
        if (std::abs(deltaTime) < thresholdTime) { // track attached to >1 vertex, remember that it was already processed
          tmap[track.globalIndex()] = {pools[poolIndex].size() - 1, poolIndex};
        }

      } // track Mother loop
    }   // collision loop

    return pools;
  }

  void calculateInvMass(CompleteCollisions const& collisions, CompleteTracks const& tracks, o2::aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcWtmp, gsl::span<std::vector<TrackCand>> trackPoolM, gsl::span<std::vector<TrackCand>> trackPoolD, int chargeM, int chargeD, int particleName, const aod::McParticles* partTable = nullptr)
  {

    int ntrInner = chargeM < 0 ? trackPoolM[NEG].size() : trackPoolM[POS].size();
    int ntrOuter = chargeD < 0 ? trackPoolD[NEG].size() : trackPoolD[POS].size();

    const int poolCh1 = chargeM < 0 ? 1 : 0;
    const int poolCh2 = chargeD < 0 ? 1 : 0;
    int motherPdg = 999;
    int daughterPdg = 777;

    switch (particleName) {
      case SigmaMinus:
        mMother = 1.197449;
        mChargedDaughter = 0.13957039;
        mNeutralDaughter = 0.9395654205;
        break;

      case SigmaPlusToPi:
        mMother = 1.18937;
        mChargedDaughter = 0.13957039;
        mNeutralDaughter = 0.9395654205;
        break;

      case SigmaPlusToProton:
        mMother = 1.18937;
        mChargedDaughter = 0.93827208816;
        mNeutralDaughter = 0.1349768;
        break;

      case Kaon:
        mMother = 0.493677;
        mChargedDaughter = 0.1056583755;
        mNeutralDaughter = 0.;
        break;

      case Xi:
        mMother = 1.32171;
        mChargedDaughter = 0.13957039;
        mNeutralDaughter = 1.115683;
        break;

      case OmegaToL:
        mMother = 1.67245;
        mChargedDaughter = 0.493677;
        mNeutralDaughter = 1.115683;
        break;

      case OmegaToXi:
        mMother = 1.67245;
        mChargedDaughter = 0.13957039;
        mNeutralDaughter = 1.31486;
        break;
    }

    o2::dataformats::VertexBase primaryVertex;

    std::array<std::vector<int>, 2> mVtxSecondTrack{}; // 1st pos. and neg. track of the pools for each vertex

    for (int i = 0; i < 2; i++) {
      mVtxSecondTrack[i].clear();
      mVtxSecondTrack[i].resize(collisions.size(), -1);
    }

    for (int pn = 0; pn < 2; pn++) {
      auto& vtxFirstT = mVtxSecondTrack[pn];
      const auto& tracksPool = trackPoolD[pn];
      for (unsigned i = 0; i < tracksPool.size(); i++) {
        const auto& t = tracksPool[i];
        for (int j{t.vBracket.getMin()}; j <= t.vBracket.getMax(); ++j) {
          if (vtxFirstT[j] == -1) {
            vtxFirstT[j] = i;
          }
        }
      }
    }

    for (int itp = 0; itp < ntrInner; itp++) { // HERE change neg->pos to get the other charge!!!
      auto& seedM = trackPoolM[poolCh1][itp];
      int firstD = mVtxSecondTrack[poolCh2][seedM.vBracket.getMin()];
      if (firstD < 0) {
        LOG(debug) << "No partner is found for pos.track " << itp << " out of " << ntrInner;
        continue;
      }

      const auto& trackM = tracks.iteratorAt(seedM.Idxtr);
      if ((seedM.mcParticleIdx != -1) && partTable) {
        auto mcParticle = partTable->rawIteratorAt(seedM.mcParticleIdx);
        motherPdg = mcParticle.pdgCode();
      }

      if (trackM.has_collision()) {
        auto const& collision = trackM.collision_as<CompleteCollisions>();
        primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
        primaryVertex.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
        histos.fill(HIST("hRecVtxZData"), primaryVertex.getZ());
      } else {
        primaryVertex.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
        histos.fill(HIST("hRecVtxZData"), primaryVertex.getZ());
      }

      SigmaTr = getTrackParCov(trackM);
      o2::base::Propagator::Instance()->PropagateToXBxByBz(SigmaTr, kink::LayerRadii[trackM.itsNCls() - 1]);

      SigmaTr2 = getTrackParCov(trackM);

      gpu::gpustd::array<float, 2> dcaInfoM;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, SigmaTr2, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfoM);
      auto motherdcaXY = abs(dcaInfoM[0]);
      histos.fill(HIST("hDCAMotherToPV"), motherdcaXY);

      if (motherdcaXY > 0.03)
        continue;

      for (int itn = firstD; itn < ntrOuter; itn++) {
        auto& seedD = trackPoolD[poolCh2][itn];
        if (seedD.vBracket.isOutside(seedM.vBracket)) {
          LOG(debug) << "Brackets do not match";
          continue;
        }

        const auto& trackDgh = tracks.iteratorAt(static_cast<uint64_t>(seedD.Idxtr));

        if ((seedD.mcParticleIdx != -1) && partTable) {
          auto mcParticle = partTable->rawIteratorAt(seedD.mcParticleIdx);
          daughterPdg = mcParticle.pdgCode();
        }
        PionTr = getTrackParCov(trackDgh);

        SigmaTr2.getPxPyPzGlo(sigmaPin);

        if ((particleName == SigmaMinus) || (particleName == SigmaPlusToPi) || (particleName == Xi) || (particleName == OmegaToXi)) {
          if (trackDgh.tpcNSigmaPi() > cfgNsigmaTPCdaughter)
            continue;
        }
        if (particleName == SigmaPlusToProton) {
          if (trackDgh.tpcNSigmaPr() > cfgNsigmaTPCdaughter)
            continue;
        }
        if (particleName == OmegaToL) {
          if (trackDgh.tpcNSigmaKa() > cfgNsigmaTPCdaughter)
            continue;
        }

        gpu::gpustd::array<float, 2> dcaInfo;
        o2::base::Propagator::Instance()->propagateToDCABxByBz({primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, PionTr, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);
        auto kinkdcaXY = abs(dcaInfo[0]);
        histos.fill(HIST("hDCAdaughterToPV"), kinkdcaXY);

        if (kinkdcaXY < dcaKinkDtopv)
          continue;

        // check how much the tracks are close in space
        if ((cfgIsMC) && (motherPdg == particlePdgCode)) {
          histos.fill(HIST("hdeltaphi"), (SigmaTr.getPhi() - PionTr.getPhi()) * radToDeg);
          histos.fill(HIST("hdeltaz"), SigmaTr.getZ() - PionTr.getZ());
        }

        if (std::abs(SigmaTr.getZ() - PionTr.getZ()) > zDiff)
          continue;
        if ((std::abs(SigmaTr.getPhi() - PionTr.getPhi()) * radToDeg) > phiDiff)
          continue;

        {
          try {
            ft2.process(PionTr, SigmaTr);
            ft2.propagateTracksToVertex();
            if (ft2.isPropagateTracksToVertexDone() == true) {
              auto SigmaTrDCA = ft2.getTrack(1);
              auto PionTrDCA = ft2.getTrack(0);
              SigmaTrDCA.getPxPyPzGlo(sigmaPDC);
              sigmaPabsDC = SigmaTrDCA.getP();
              sigmaPt = SigmaTr.getPt();
              etaS = SigmaTr.getEta();
              phiS = SigmaTr.getPhi();

              PionTrDCA.getPxPyPzGlo(pionPDC);
              pionPabsDC = PionTrDCA.getP();
              etaP = PionTr.getEta();
              phiP = PionTr.getPhi();

              pionE = 0;
              neutronE = 0;

              if (ft2.getChi2AtPCACandidate() < 0)
                continue;
              std::array<float, 3> R = ft2.getPCACandidatePos();

              pionE = sqrt(mChargedDaughter * mChargedDaughter + pionPabsDC * pionPabsDC);

              costheta = (sigmaPDC[0] * pionPDC[0] + sigmaPDC[1] * pionPDC[1] + sigmaPDC[2] * pionPDC[2]) / (sigmaPabsDC * pionPabsDC);
              theta = std::acos(costheta);

              qT = pionPabsDC * std::sin(theta);

              histos.fill(HIST("hqT"), qT);

              sigmaE = sqrt(mMother * mMother + sigmaPabsDC * sigmaPabsDC);
              neutronPabs = sqrt(pow((sigmaPDC[2] - pionPDC[2]), 2) + pow((sigmaPDC[1] - pionPDC[1]), 2) + pow((sigmaPDC[0] - pionPDC[0]), 2));
              neutronM = sqrt((sigmaE - pionE) * (sigmaE - pionE) - neutronPabs * neutronPabs);

              if ((particleName == SigmaMinus) || (particleName == SigmaPlusToPi)) {
                if ((theta * radToDeg > (angleCutFunction(particleName, sigmaPabsDC))) && (sigmaPabsDC > 1.6))
                  continue;
              }
              if (particleName == SigmaPlusToProton) {
                if ((theta * radToDeg > (angleCutFunction(particleName, sigmaPabsDC))) && (sigmaPabsDC > 0.3))
                  continue;
              }
              if (particleName == Kaon) {
                Double_t maxDecAngKmu = angleCutFunction(particleName, sigmaPabsDC);
                Double_t maxDecAngpimu = angleCutFunction(Pion, sigmaPabsDC);

                if ((theta * radToDeg < maxDecAngpimu * 1.2))
                  continue;
                if ((theta * radToDeg > maxDecAngKmu * .98) && (sigmaPabsDC > 1.2))
                  continue;
              }
              if (particleName == Xi) {
                if (sigmaPabsDC < 2.)
                  if (theta * radToDeg < 2.5)
                    continue;
                if (sigmaPabsDC > 1.2)
                  if (theta * radToDeg > angleCutFunction(particleName, sigmaPabsDC))
                    continue;
              }
              if ((particleName == OmegaToL) || (particleName == OmegaToXi)) {
                if (sigmaPabsDC < 5.5)
                  if (theta * radToDeg < 4.)
                    continue;
                if (sigmaPabsDC > 0.6)
                  if (theta * radToDeg > angleCutFunction(particleName, sigmaPabsDC))
                    continue;
              }

              if (!cfgIsMC)
                if (sigmaPt < 1.6)
                  continue;

              if (theta * radToDeg < 0.5)
                continue;

              if ((qT < qTlower) || (qT > qTupper))
                continue;

              if (sqrt(R[0] * R[0] + R[1] * R[1]) < 18.)
                continue;

              if (chargeM == chargeD)
                mKinkCandidates.emplace_back(KinkCandidates{chargeM, seedM.Idxtr, seedD.Idxtr, seedM.mcParticleIdx, seedD.mcParticleIdx, primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(), R[0], R[1], R[2], sigmaPin[0], sigmaPin[1], sigmaPin[2], sigmaPDC[0], sigmaPDC[1], sigmaPDC[2], pionPDC[0], pionPDC[1], pionPDC[2]});

              neutronE = sqrt(mNeutralDaughter * mNeutralDaughter + pow((sigmaPDC[2] - pionPDC[2]), 2) + pow((sigmaPDC[1] - pionPDC[1]), 2) + pow((sigmaPDC[0] - pionPDC[0]), 2));

              mass = sqrt((neutronE + pionE) * (neutronE + pionE) - sigmaPabsDC * sigmaPabsDC);

              if ((chargeM == -1) && (chargeD == -1)) {
                if (cfgIsMC) {
                  histos.fill(HIST("hcodes"), motherPdg, daughterPdg);
                  if ((motherPdg == particlePdgCode || motherPdg == -3222) && (daughterPdg == -211)) {
                    histos.fill(HIST("hPtMinusRecMcTrth"), mass, sigmaPt);
                    histos.fill(HIST("hptMDtrue"), sigmaPt, PionTr.getPt());
                  } else if ((motherPdg == particlePdgCode || motherPdg == -3222) && (daughterPdg != -211)) {
                    histos.fill(HIST("hptMtrue"), sigmaPt, PionTr.getPt());
                    histos.fill(HIST("hPtMinusRecMcTrthM"), mass, sigmaPt);
                  } else { // if ((motherPdg != particlePdgCode)&&(daughterPdg!=-211)) {
                    histos.fill(HIST("hptMDelse"), sigmaPt, PionTr.getPt());
                    histos.fill(HIST("hPtMinusRecMcTrthelse"), mass, sigmaPt);
                  }
                }
                histos.fill(HIST("hMassMinusPt"), mass, sigmaPt);
              }
              if ((chargeM == 1) && (chargeD == 1)) {
                if (cfgIsMC) {
                  if (motherPdg == particlePdgCode)
                    histos.fill(HIST("hPtPlusRecMcTrth"), mass, sigmaPt);
                }
                histos.fill(HIST("hMassPlusPt"), mass, sigmaPt);
              }

              if ((chargeM == -1) && (chargeD == 1)) {
                histos.fill(HIST("hBackgroundNegPosPt"), mass, sigmaPt);
              }
              if ((chargeM == 1) && (chargeD == -1)) {
                histos.fill(HIST("hBackgroundPosNegPt"), mass, sigmaPt);
              }
            } // true
          } catch (std::runtime_error& e) {
            continue;
          }
          //} //dca fitter option

        } // try

      } // inner track loop
    }   // outer track loop
  }

  void processReco(CompleteCollisions const& collisions, CompleteTracks const& tracks, o2::aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcWtmp)
  {
    auto firstcollision = collisions.begin();
    auto bc1 = firstcollision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc1);

    auto trackPools = makeTracksPool(collisions, tracks, ambiTracks, bcWtmp);

    LOG(info) << "Collected for mother " << trackPools[POS].size() << " positive and " << trackPools[NEG].size() << " negative seeds";

    LOG(info) << "Collected for daughter " << trackPools[2 + POS].size() << " positive and " << trackPools[2 + NEG].size() << " negative seeds";

    gsl::span<std::vector<TrackCand>> trackPoolsMth{trackPools.data(), 2};
    gsl::span<std::vector<TrackCand>> trackPoolsDhgt{trackPools.data() + 2, 2};

    calculateInvMass(collisions, tracks, ambiTracks, bcWtmp, trackPoolsMth, trackPoolsDhgt, cfgMotherCharge, cfgDaughterCharge, cfgParticleSpec);
    calculateInvMass(collisions, tracks, ambiTracks, bcWtmp, trackPoolsMth, trackPoolsDhgt, -1, +1, cfgParticleSpec);
    calculateInvMass(collisions, tracks, ambiTracks, bcWtmp, trackPoolsMth, trackPoolsDhgt, +1, -1, cfgParticleSpec);

  } // process

  PROCESS_SWITCH(kinkAnalysis, processReco, "process reconstructed information", true);

  void processSim(CompleteCollisions const& collisions, CompleteTracks const& tracks, o2::aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcWtmp, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
  {
    auto firstcollision = collisions.begin();
    auto bc1 = firstcollision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc1);

    for (const auto& mcParticle : particlesMC) {
      if (mcParticle.pdgCode() != particlePdgCode)
        continue;
      histos.fill(HIST("generatedPt"), mcParticle.pt());
    }

    auto trackPools = makeTracksPool(collisions, tracks, ambiTracks, bcWtmp);

    for (auto& pool : trackPools) {
      for (auto& track : pool) {
        auto label = trackLabelsMC.iteratorAt(track.Idxtr);
        if (label.mcParticleId() < -1 || label.mcParticleId() >= particlesMC.size()) {
          continue;
        }
        track.mcParticleIdx = label.mcParticleId();
      }
    }

    LOG(info) << "Collected for mother " << trackPools[POS].size() << " positive and " << trackPools[NEG].size() << " negative seeds";

    LOG(info) << "Collected for daughter " << trackPools[2 + POS].size() << " positive and " << trackPools[2 + NEG].size() << " negative seeds";

    gsl::span<std::vector<TrackCand>> trackPoolsMth{trackPools.data(), 2};
    gsl::span<std::vector<TrackCand>> trackPoolsDhgt{trackPools.data() + 2, 2};

    calculateInvMass(collisions, tracks, ambiTracks, bcWtmp, trackPoolsMth, trackPoolsDhgt, cfgMotherCharge, cfgDaughterCharge, cfgParticleSpec, &particlesMC);

    // for (auto& candidates : mKinkCandidates) {

    // }

    calculateInvMass(collisions, tracks, ambiTracks, bcWtmp, trackPoolsMth, trackPoolsDhgt, -1, +1, cfgParticleSpec, &particlesMC);
    calculateInvMass(collisions, tracks, ambiTracks, bcWtmp, trackPoolsMth, trackPoolsDhgt, +1, -1, cfgParticleSpec, &particlesMC);
  } // process

  PROCESS_SWITCH(kinkAnalysis, processSim, "process sim information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kinkAnalysis>(cfgc)};
}
