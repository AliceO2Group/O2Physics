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

/// \file   kinkBuilder.cxx
/// \brief Builder task for kink decay topologies using ITS standalone tracks for the mother
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>

#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/Utils/svPoolCreator.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using VBracket = o2::math_utils::Bracket<int>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;

namespace
{
constexpr std::array<float, 7> LayerRadii{2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"Daughter"};

std::shared_ptr<TH2> h2ClsMapPtMoth;
std::shared_ptr<TH2> h2ClsMapPtDaug;
std::shared_ptr<TH2> h2DeDxDaugSel;
std::shared_ptr<TH2> h2KinkAnglePt;
std::shared_ptr<TH2> h2MothMassPt;
} // namespace

struct kinkCandidate {

  float recoPtMoth() const { return std::hypot(momMoth[0], momMoth[1]); }
  float recoPhiMoth() const { return std::atan2(momMoth[1], momMoth[0]); }
  float recoEtaMoth() const { return std::asinh(momMoth[2] / recoPtMoth()); }

  float recoPtDaug() const { return std::hypot(momDaug[0], momDaug[1]); }
  float recoPhiDaug() const { return std::atan2(momDaug[1], momDaug[0]); }
  float recoEtaDaug() const { return std::asinh(momDaug[2] / recoPtDaug()); }

  int mothTrackID;
  int daugTrackID;
  int collisionID;

  int mothSign;
  std::array<float, 3> momMoth = {-999, -999, -999};
  std::array<float, 3> momDaug = {-999, -999, -999};
  std::array<float, 3> primVtx = {-999, -999, -999};
  std::array<float, 3> decVtx = {-999, -999, -999};

  float dcaKinkTopo = -999;
  float nSigmaTPCDaug = -999;
  float nSigmaTOFDaug = -999;
  float dcaXYdaug = -999;
  float dcaXYmoth = -999;
  float kinkAngle = -999;
};

struct kinkBuilder {

  enum PartType { kSigmaMinus = 0,
                  kHypertriton,
                  kHyperhelium4sigma };

  Produces<aod::KinkCands> outputDataTable;
  Produces<aod::KinkCandsUnbound> outputDataTableUB;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> hypoMoth{"hypoMoth", kSigmaMinus, "Mother particle hypothesis"};
  Configurable<bool> fillDebugTable{"fillDebugTable", false, "If true, fill the debug table with all candidates unbound"};
  // Selection criteria
  Configurable<float> maxDCAMothToPV{"maxDCAMothToPV", 0.1, "Max DCA of the mother to the PV"};
  Configurable<float> minDCADaugToPV{"minDCADaugToPV", 0., "Min DCA of the daughter to the PV"};
  Configurable<float> minPtMoth{"minPtMoth", 0.5, "Minimum pT of the hypercandidate"};
  Configurable<float> maxZDiff{"maxZDiff", 20., "Max z difference between the kink daughter and the mother"};
  Configurable<float> maxPhiDiff{"maxPhiDiff", 100, "Max phi difference between the kink daughter and the mother"};
  Configurable<float> timeMarginNS{"timeMarginNS", 600, "Additional time res tolerance in ns"};
  Configurable<float> etaMax{"etaMax", 1., "eta daughter"};
  Configurable<float> nTPCClusMinDaug{"nTPCClusMinDaug", 80, "daug NTPC clusters cut"};
  Configurable<bool> askTOFforDaug{"askTOFforDaug", false, "If true, ask for TOF signal"};
  Configurable<bool> doSVRadiusCut{"doSVRadiusCut", true, "If true, apply the cut on the radius of the secondary vertex and tracksIU"};
  Configurable<bool> updateMothTrackUsePV{"updateMothTrackUsePV", false, "If true, update the mother track parameters using the primary vertex"};

  o2::vertexing::DCAFitterN<2> fitter;
  o2::base::MatLayerCylSet* lut = nullptr;

  // constants
  float radToDeg = o2::constants::math::Rad2Deg;
  svPoolCreator svCreator;

  // bethe bloch parameters
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for charged daughter"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};
  Configurable<float> customVertexerTimeMargin{"customVertexerTimeMargin", 800, "Time margin for custom vertexer (ns)"};
  Configurable<bool> skipAmbiTracks{"skipAmbiTracks", false, "Skip ambiguous tracks"};
  Configurable<bool> unlikeSignBkg{"unlikeSignBkg", false, "Use unlike sign background"};

  // CCDB options
  Configurable<std::string> ccdbPath{"ccdbPath", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};

  // PDG codes

  // histogram axes
  ConfigurableAxis rigidityBins{"rigidityBins", {200, -10.f, 10.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis dedxBins{"dedxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {200, -5.f, 5.f}, "Binning for n sigma"};

  // std vector of candidates
  std::vector<kinkCandidate> kinkCandidates;

  HistogramRegistry qaRegistry{"QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float mBz;
  std::array<float, 6> mBBparamsDaug;

  // mother and daughter tracks' properties (absolute charge and mass)
  int charge = 1;
  float mothMass = o2::constants::physics::MassSigmaMinus;
  float chargedDauMass = o2::constants::physics::MassPiMinus;
  float neutDauMass = o2::constants::physics::MassNeutron;

  void init(InitContext const&)
  {
    if (hypoMoth == kSigmaMinus) {
      charge = 1;
      mothMass = o2::constants::physics::MassSigmaMinus;
      chargedDauMass = o2::constants::physics::MassPiMinus;
      neutDauMass = o2::constants::physics::MassNeutron;
    } else if (hypoMoth == kHypertriton) {
      charge = 1;
      mothMass = o2::constants::physics::MassHyperTriton;
      chargedDauMass = o2::constants::physics::MassTriton;
      neutDauMass = o2::constants::physics::MassPi0;
    } else if (hypoMoth == kHyperhelium4sigma) {
      charge = 2;
      mothMass = o2::constants::physics::MassHyperHelium4;
      chargedDauMass = o2::constants::physics::MassAlpha;
      neutDauMass = o2::constants::physics::MassPi0;
    }

    // dummy values, 1 for mother, 0 for daughter
    svCreator.setPDGs(1, 0);

    mRunNumber = 0;
    mBz = 0;

    ccdb->setURL(ccdbPath);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);

    svCreator.setTimeMargin(customVertexerTimeMargin);
    if (skipAmbiTracks) {
      svCreator.setSkipAmbiTracks();
    }

    const AxisSpec itsClusterMapAxis(128, 0, 127, "ITS cluster map");
    const AxisSpec rigidityAxis{rigidityBins, "#it{p}^{TPC}/#it{z}"};
    const AxisSpec ptAxis{rigidityBins, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec kinkAngleAxis{100, 0, 180, "#theta_{kink} (deg)"};
    const AxisSpec dedxAxis{dedxBins, "d#it{E}/d#it{x}"};

    AxisSpec massAxis(100, 1.1, 1.4, "m (GeV/#it{c}^{2})");
    if (hypoMoth == kSigmaMinus) {
      massAxis = AxisSpec{100, 1.1, 1.4, "m (GeV/#it{c}^{2})"};
    } else if (hypoMoth == kHypertriton) {
      massAxis = AxisSpec{100, 2.94, 3.2, "m (GeV/#it{c}^{2})"};
    } else if (hypoMoth == kHyperhelium4sigma) {
      massAxis = AxisSpec{100, 3.85, 4.25, "m (GeV/#it{c}^{2})"};
    }

    h2DeDxDaugSel = qaRegistry.add<TH2>("h2DeDxDaugSel", ";p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dedxAxis});
    h2KinkAnglePt = qaRegistry.add<TH2>("h2KinkAnglePt", "; p_{T} (GeV/#it{c}); #theta_{kink} (deg)", HistType::kTH2F, {ptAxis, kinkAngleAxis});
    h2MothMassPt = qaRegistry.add<TH2>("h2MothMassPt", "; p_{T} (GeV/#it{c}); m (GeV/#it{c}^{2})", HistType::kTH2F, {ptAxis, massAxis});
    h2ClsMapPtMoth = qaRegistry.add<TH2>("h2ClsMapPtMoth", "; p_{T} (GeV/#it{c}); ITS cluster map", HistType::kTH2F, {ptAxis, itsClusterMapAxis});
    h2ClsMapPtDaug = qaRegistry.add<TH2>("h2ClsMapPtDaug", "; p_{T} (GeV/#it{c}); ITS cluster map", HistType::kTH2F, {ptAxis, itsClusterMapAxis});

    for (int i = 0; i < 5; i++) {
      mBBparamsDaug[i] = cfgBetheBlochParams->get("Daughter", Form("p%i", i));
    }
    mBBparamsDaug[5] = cfgBetheBlochParams->get("Daughter", "resolution");
  }

  template <typename T>
  bool selectMothTrack(const T& candidate)
  {
    if (candidate.has_collision() && candidate.hasITS() && !candidate.hasTPC() && !candidate.hasTOF() && candidate.itsNCls() < 6 &&
        candidate.itsNClsInnerBarrel() == 3 && candidate.itsChi2NCl() < 36 && candidate.pt() > minPtMoth) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectDaugTrack(const T& candidate)
  {
    if (!candidate.hasTPC() || !candidate.hasITS()) {
      return false;
    }

    if (askTOFforDaug && !candidate.hasTOF()) {
      return false;
    }

    bool isGoodTPCCand = false;
    if (candidate.itsNClsInnerBarrel() == 0 && candidate.itsNCls() < 4 &&
        candidate.tpcNClsCrossedRows() > 0.8 * candidate.tpcNClsFindable() && candidate.tpcNClsFound() > nTPCClusMinDaug) {
      isGoodTPCCand = true;
    }

    if (!isGoodTPCCand) {
      return false;
    }

    return true;
  }

  template <class Tcolls, class Ttracks>
  void fillCandidateData(const Tcolls& collisions, const Ttracks& tracks, aod::AmbiguousTracks const& ambiguousTracks, aod::BCs const& bcs)
  {
    svCreator.clearPools();
    svCreator.fillBC2Coll(collisions, bcs);

    for (const auto& track : tracks) {
      if (std::abs(track.eta()) > etaMax)
        continue;

      bool isDaug = selectDaugTrack(track);
      bool isMoth = selectMothTrack(track);

      if (!isDaug && !isMoth)
        continue;

      int pdgHypo = isMoth ? 1 : 0;
      svCreator.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }
    auto& kinkPool = svCreator.getSVCandPool(collisions, !unlikeSignBkg);

    for (const auto& svCand : kinkPool) {
      kinkCandidate kinkCand;

      auto trackMoth = tracks.rawIteratorAt(svCand.tr0Idx);
      auto trackDaug = tracks.rawIteratorAt(svCand.tr1Idx);

      auto const& collision = trackMoth.template collision_as<Tcolls>();
      auto const& bc = collision.template bc_as<aod::BCs>();
      initCCDB(bc);

      o2::dataformats::VertexBase primaryVertex;
      primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
      primaryVertex.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      kinkCand.primVtx = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};

      o2::track::TrackParCov trackParCovMoth = getTrackParCov(trackMoth);
      o2::track::TrackParCov trackParCovMothPV{trackParCovMoth};
      o2::base::Propagator::Instance()->PropagateToXBxByBz(trackParCovMoth, LayerRadii[trackMoth.itsNCls() - 1]);

      std::array<float, 2> dcaInfoMoth;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, trackParCovMothPV, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfoMoth);

      if (std::abs(dcaInfoMoth[0]) > maxDCAMothToPV) {
        continue;
      }

      o2::track::TrackParCov trackParCovDaug = getTrackParCov(trackDaug);

      // check if the kink daughter is close to the mother
      if (std::abs(trackParCovMoth.getZ() - trackParCovDaug.getZ()) > maxZDiff) {
        continue;
      }

      if ((std::abs(trackParCovMoth.getPhi() - trackParCovDaug.getPhi()) * radToDeg) > maxPhiDiff) {
        continue;
      }

      // propagate to PV
      std::array<float, 2> dcaInfoDaug;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, trackParCovDaug, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfoDaug);
      if (std::abs(dcaInfoDaug[0]) < minDCADaugToPV) {
        continue;
      }

      if (updateMothTrackUsePV) {
        // update the mother track parameters using the primary vertex
        trackParCovMoth = trackParCovMothPV;
        if (!trackParCovMoth.update(primaryVertex)) {
          continue;
        }
      }

      int nCand = 0;
      try {
        nCand = fitter.process(trackParCovMoth, trackParCovDaug);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        continue;
      }
      if (nCand == 0) {
        continue;
      }

      if (!fitter.propagateTracksToVertex()) {
        continue;
      }

      auto propMothTrack = fitter.getTrack(0);
      auto propDaugTrack = fitter.getTrack(1);
      kinkCand.decVtx = fitter.getPCACandidatePos();

      // cut on decay radius to 17 cm
      float decRad2 = kinkCand.decVtx[0] * kinkCand.decVtx[0] + kinkCand.decVtx[1] * kinkCand.decVtx[1];
      if (doSVRadiusCut && decRad2 < LayerRadii[3] * LayerRadii[3]) {
        continue;
      }

      // get last layer hitted by the mother and the first layer hitted by the daughter
      int lastLayerMoth = 0, firstLayerDaug = 0;
      for (int i = 0; i < 7; i++) {
        if (trackMoth.itsClusterMap() & (1 << i)) {
          lastLayerMoth = i;
        }
      }

      for (int i = 0; i < 7; i++) {
        if (trackDaug.itsClusterMap() & (1 << i)) {
          firstLayerDaug = i;
          break;
        }
      }

      if (doSVRadiusCut && lastLayerMoth >= firstLayerDaug) {
        continue;
      }

      if (doSVRadiusCut && decRad2 < LayerRadii[lastLayerMoth] * LayerRadii[lastLayerMoth]) {
        continue;
      }

      for (int i = 0; i < 3; i++) {
        kinkCand.decVtx[i] -= kinkCand.primVtx[i];
      }

      propMothTrack.getPxPyPzGlo(kinkCand.momMoth);
      propDaugTrack.getPxPyPzGlo(kinkCand.momDaug);
      for (int i = 0; i < 3; i++) {
        kinkCand.momMoth[i] *= charge;
        kinkCand.momDaug[i] *= charge;
      }
      float pMoth = propMothTrack.getP() * charge;
      float pDaug = propDaugTrack.getP() * charge;
      float spKink = kinkCand.momMoth[0] * kinkCand.momDaug[0] + kinkCand.momMoth[1] * kinkCand.momDaug[1] + kinkCand.momMoth[2] * kinkCand.momDaug[2];
      kinkCand.kinkAngle = std::acos(spKink / (pMoth * pDaug));

      std::array<float, 3> neutDauMom{0.f, 0.f, 0.f};
      for (int i = 0; i < 3; i++) {
        neutDauMom[i] = kinkCand.momMoth[i] - kinkCand.momDaug[i];
      }

      float chargedDauE = std::sqrt(pDaug * pDaug + chargedDauMass * chargedDauMass);
      float neutE = std::sqrt(neutDauMom[0] * neutDauMom[0] + neutDauMom[1] * neutDauMom[1] + neutDauMom[2] * neutDauMom[2] + neutDauMass * neutDauMass);
      float invMass = std::sqrt((chargedDauE + neutE) * (chargedDauE + neutE) - (pMoth * pMoth));

      h2DeDxDaugSel->Fill(trackDaug.tpcInnerParam() * trackDaug.sign(), trackDaug.tpcSignal());
      h2KinkAnglePt->Fill(trackMoth.pt() * charge * trackMoth.sign(), kinkCand.kinkAngle * radToDeg);
      h2MothMassPt->Fill(trackMoth.pt() * charge * trackMoth.sign(), invMass);
      h2ClsMapPtMoth->Fill(trackMoth.pt() * charge * trackMoth.sign(), trackMoth.itsClusterMap());
      h2ClsMapPtDaug->Fill(trackDaug.pt() * charge * trackDaug.sign(), trackDaug.itsClusterMap());

      kinkCand.collisionID = collision.globalIndex();
      kinkCand.mothTrackID = trackMoth.globalIndex();
      kinkCand.daugTrackID = trackDaug.globalIndex();

      kinkCand.dcaXYmoth = dcaInfoMoth[0];
      kinkCand.mothSign = trackMoth.sign();
      kinkCand.dcaXYdaug = dcaInfoDaug[0];
      kinkCand.dcaKinkTopo = std::sqrt(fitter.getChi2AtPCACandidate());
      kinkCandidates.push_back(kinkCand);
    }
  }

  void initCCDB(aod::BCs::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOG(info) << "Initializing CCDB for run " << mRunNumber;
    o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(grpmagPath, mRunNumber);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    mBz = grpmag->getNominalL3Field();
    fitter.setBz(mBz);

    if (!lut) {
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
      int mat{static_cast<int>(cfgMaterialCorrection)};
      fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));
    }
    o2::base::Propagator::Instance()->setMatLUT(lut);
    LOG(info) << "Task initialized for run " << mRunNumber << " with magnetic field " << mBz << " kZG";
  }

  void process(aod::Collisions const& collisions, TracksFull const& tracks, aod::AmbiguousTracks const& ambiTracks, aod::BCs const& bcs)
  {

    kinkCandidates.clear();
    fillCandidateData(collisions, tracks, ambiTracks, bcs);

    // sort kinkCandidates by collisionID to allow joining with collision table
    std::sort(kinkCandidates.begin(), kinkCandidates.end(), [](const kinkCandidate& a, const kinkCandidate& b) { return a.collisionID < b.collisionID; });

    for (const auto& kinkCand : kinkCandidates) {
      if (fillDebugTable) {
        outputDataTableUB(kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
                          kinkCand.mothSign, kinkCand.momMoth[0], kinkCand.momMoth[1], kinkCand.momMoth[2],
                          kinkCand.momDaug[0], kinkCand.momDaug[1], kinkCand.momDaug[2],
                          kinkCand.dcaXYmoth, kinkCand.dcaXYdaug, kinkCand.dcaKinkTopo);
      } else {
        outputDataTable(kinkCand.collisionID, kinkCand.mothTrackID, kinkCand.daugTrackID,
                        kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
                        kinkCand.mothSign, kinkCand.momMoth[0], kinkCand.momMoth[1], kinkCand.momMoth[2],
                        kinkCand.momDaug[0], kinkCand.momDaug[1], kinkCand.momDaug[2],
                        kinkCand.dcaXYmoth, kinkCand.dcaXYdaug, kinkCand.dcaKinkTopo);
      }
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kinkBuilder>(cfgc)};
}
