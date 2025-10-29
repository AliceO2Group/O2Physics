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

/// \file   spectraKinkPiKa.cxx
/// \brief Example of a simple task for the analysis of the muon from Kaon pion using kink topology
/// \author sandeep dudi sandeep.dudi@cern.ch

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

//////////////
#include "Common/DataModel/PIDResponse.h"

#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/PID.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TPDGCode.h"
#include "TVector3.h"
#include <TMath.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using VBracket = o2::math_utils::Bracket<int>;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCMu>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel, aod::FT0Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSel, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;
namespace
{
constexpr std::array<float, 7> LayerRadii{2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"Daughter"};

} // namespace

struct kinkCandidate {
  int mothTrackID;
  int daugTrackID;
  int collisionID;

  int mothSign;
  std::array<float, 3> momMoth = {-999, -999, -999};
  std::array<float, 3> momDaug = {-999, -999, -999};
  std::array<float, 3> primVtx = {-999, -999, -999};
  std::array<float, 3> decVtx = {-999, -999, -999};

  float dcaKinkTopo = -999;
  float dcaXYdaug = -999;
  float dcaXYmoth = -999;
  float kinkAngle = -999;
};
struct kinkBuilder {
  // kink analysis
  Produces<aod::KinkCands> outputDataTable;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // Selection criteria
  Configurable<float> maxDCAMothToPV{"maxDCAMothToPV", 0.2, "Max DCA of the mother to the PV"};
  Configurable<float> minDCADaugToPV{"minDCADaugToPV", 0., "Min DCA of the daughter to the PV"};
  Configurable<float> minPtMoth{"minPtMoth", 0.15, "Minimum pT of the hypercandidate"};
  Configurable<float> maxZDiff{"maxZDiff", 20., "Max z difference between the kink daughter and the mother"};
  Configurable<float> maxPhiDiff{"maxPhiDiff", 100, "Max phi difference between the kink daughter and the mother"};
  Configurable<float> timeMarginNS{"timeMarginNS", 600, "Additional time res tolerance in ns"};
  Configurable<float> etaMaxDaug{"etaMaxDaug", 1., "eta max daughter"};
  Configurable<float> etaMaxMoth{"etaMaxMoth", 1., "eta max Mother"};
  Configurable<float> nTPCClusMinDaug{"nTPCClusMinDaug", 30, "mother NTPC clusters cut"};
  Configurable<float> itsChi2cut{"itsChi2cut", 36, "mother itsChi2 cut"};
  Configurable<bool> askTOFforDaug{"askTOFforDaug", false, "If true, ask for TOF signal"};
  Configurable<bool> kaontopologhy{"kaontopologhy", true, "If true, selected mother have both ITS+TPC "};

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

  // std vector of candidates
  std::vector<kinkCandidate> kinkCandidates;
  int mRunNumber;
  float mBz;
  std::array<float, 6> mBBparamsDaug;

  // mother and daughter tracks' properties (absolute charge and mass)
  int charge = 1;
  void init(InitContext const&)
  {
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
    for (int i = 0; i < 5; i++) {
      mBBparamsDaug[i] = cfgBetheBlochParams->get("Daughter", Form("p%i", i));
    }
    mBBparamsDaug[5] = cfgBetheBlochParams->get("Daughter", "resolution");
  }

  template <typename T>
  bool selectMothTrack(const T& candidate)
  {
    // ITS-standalone (no TPC, no TOF)
    if (!kaontopologhy) {
      if (candidate.has_collision() && candidate.hasITS() && !candidate.hasTPC() && !candidate.hasTOF() &&
          candidate.itsNCls() < 6 &&
          candidate.itsNClsInnerBarrel() == 3 &&
          candidate.itsChi2NCl() < 36 &&
          candidate.pt() > minPtMoth) {
        return true;
      }
      return false;
    }
    // Kaon topology: ITS+TPC, no TOF
    if (kaontopologhy) {
      if (candidate.has_collision() && candidate.hasITS() && candidate.hasTPC() && !candidate.hasTOF() &&
          candidate.pt() > minPtMoth &&
          candidate.tpcNClsCrossedRows() >= nTPCClusMinDaug &&
          candidate.itsChi2NCl() <= itsChi2cut) {
        return true;
      }
      return false;
    }

    return false; // fallback
  }

  template <typename T>
  bool selectDaugTrack(const T& candidate)
  {
    if (!kaontopologhy && (!candidate.hasTPC() || !candidate.hasITS())) {
      return false;
    }

    if (kaontopologhy && (!candidate.hasTPC() || candidate.hasITS())) {
      return false;
    }

    if (askTOFforDaug && !candidate.hasTOF()) {
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

      bool isDaug = selectDaugTrack(track);
      bool isMoth = selectMothTrack(track);

      if (!isDaug && !isMoth)
        continue;
      if (isDaug && std::abs(track.eta()) > etaMaxDaug)
        continue;

      if (isMoth && std::abs(track.eta()) > etaMaxMoth)
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
      outputDataTable(kinkCand.collisionID, kinkCand.mothTrackID, kinkCand.daugTrackID,
                      kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
                      kinkCand.mothSign, kinkCand.momMoth[0], kinkCand.momMoth[1], kinkCand.momMoth[2],
                      kinkCand.momDaug[0], kinkCand.momDaug[1], kinkCand.momDaug[2],
                      kinkCand.dcaXYmoth, kinkCand.dcaXYdaug, kinkCand.dcaKinkTopo);
    }
  }
  PROCESS_SWITCH(kinkBuilder, process, "Produce kink tables", false);
};

struct spectraKinkPiKa {
  Service<o2::framework::O2DatabasePDG> pdg;
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rpiKkink{"rpiKkink", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutNSigmaPi{"cutNSigmaPi", 4, "NSigmaTPCPion"};
  Configurable<float> cutNSigmaKa{"cutNSigmaKa", 4, "NSigmaTPCKaon"};
  Configurable<float> cutNSigmaMu{"cutNSigmaMu", 4, "cutNSigmaMu"};
  Configurable<float> etaCut{"etaCut", 0.8, "etaCut"};
  Configurable<float> rapCut{"rapCut", 0.8, "rapCut"};
  Configurable<float> kinkanglecut{"kinkanglecut", 2.0, "kinkanglecut"};
  Configurable<float> minradius{"minradius", 130.0, "minradiuscut"};
  Configurable<float> maxradius{"maxradius", 200.0, "maxradiuscut"};
  Configurable<float> dcaXYcut{"dcaXYcut", 0.2, "dcaXYcut"};
  Configurable<float> dcaXYcutkink{"dcaXYcutkink", 0.2, "dcaXYcutkink"};
  Configurable<float> dcaZcut{"dcaZcut", 0.2, "dcaZcut"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4.0, "tpcChi2Cut"};
  Configurable<float> minqt{"minqt", 0.12, "min qt for kaon"};
  Configurable<float> maxqt{"maxqt", 0.3, "max qt for kaon"};

  Configurable<int> centestimator{"centestimator", 0, "Select multiplicity estimator: 0 - FT0C, 1 - FT0A, 2 - FT0M, 3 - FV0A, 4 - PVTracks"};
  Configurable<int> pid{"pidMother", 321, ""};
  Configurable<int> dpid{"pidDaughter", 13, ""};
  Configurable<bool> dptCut{"dptCut", 0, ""};
  Configurable<bool> qa{"qa", 0, ""};
  Configurable<int> maxtpcncle{"maxtpcncle", 0, "max tpc find ncle"};
  Configurable<int> mintpcncle{"mintpcncle", 0, "min tpc find ncle"};
  Configurable<bool> onlykaon{"onlykaon", 0, "kaon"};
  Configurable<bool> onlypion{"onlypion", 0, "pion"};
  Configurable<bool> additionalhist{"additionalhist", 1, "additional histogram"};

  ConfigurableAxis ptAxis{"ptAxis", {150, 0, 15}, ""};
  ConfigurableAxis qtAxis{"qtAxis", {2000, 0.0, 2.0}, ""};
  ConfigurableAxis kinkAxis{"kinkAxis", {200, 0.0, 4.0}, ""};
  ConfigurableAxis etaAxis{"etaAxis", {200, -5.0, 5.0}, ""};
  ConfigurableAxis vertexAxis{"vertexAxis", {1200, -300.0, 300.0}, ""};
  ConfigurableAxis radiusAxis{"radiusAxis", {600, 0.0, 300.0}, ""};
  ConfigurableAxis massAxis{"massAxis", {600, 0.1, 0.7}, ""};
  ConfigurableAxis multAxis{"multAxis", {120, 0.0, 120.0}, ""};

  Preslice<aod::KinkCands> mPerCol = aod::track::collisionId;
  Preslice<aod::Tracks> mtPerCol = aod::track::collisionId;

  void init(InitContext const&)
  {

    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {{200, -20.0, 20.0}}});
    rEventSelection.add("hMultiplicity", "hMultiplicity", {HistType::kTH1F, {multAxis}});

    rpiKkink.add("h2_dau_pt_vs_eta_rec", "pt_vs_eta_dau", {HistType::kTH3F, {ptAxis, etaAxis, multAxis}});
    rpiKkink.add("h2_moth_pt_vs_eta_rec", "pt_vs_eta_moth", {HistType::kTH3F, {ptAxis, etaAxis, multAxis}});
    rpiKkink.add("h2_pt_moth_vs_dau_rec", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});

    rpiKkink.add("h2_qt", "qt", {HistType::kTH1F, {qtAxis}});
    rpiKkink.add("h2_qt_vs_pt", "qt_pt", {HistType::kTH2F, {qtAxis, ptAxis}});
    rpiKkink.add("h2_kink_angle", "kink angle", {HistType::kTH2F, {kinkAxis, multAxis}});

    // pion
    rpiKkink.add("h2_dau_pt_vs_eta_rec_pion", "pt_vs_eta_dau", {HistType::kTH3F, {ptAxis, etaAxis, multAxis}});
    rpiKkink.add("h2_moth_pt_vs_eta_rec_pion", "pt_vs_eta_moth", {HistType::kTH3F, {ptAxis, etaAxis, multAxis}});
    rpiKkink.add("h2_pt_moth_vs_dau_rec_pion", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});

    rpiKkink.add("h2_qt_pion", "qt", {HistType::kTH1F, {qtAxis}});
    rpiKkink.add("h2_qt_vs_ptpion", "qt_pt", {HistType::kTH2F, {qtAxis, ptAxis}});
    rpiKkink.add("h2_kink_angle_pion", "kink angle", {HistType::kTH2F, {kinkAxis, multAxis}});

    // inv mass
    rpiKkink.add("h2_invmass_kaon", "Inv mass vs Pt", {HistType::kTH3F, {massAxis, ptAxis, ptAxis}});
    rpiKkink.add("h2_invmass_pion", "Inv mass vs Pt", {HistType::kTH3F, {massAxis, ptAxis, ptAxis}});

    rpiKkink.add("h2_kaon_data", "h2_kaon_data", HistType::kTHnSparseF, {massAxis, ptAxis, etaAxis, qtAxis, multAxis}, true);
    rpiKkink.add("h2_pion_data", "h2_pion_data", HistType::kTHnSparseF, {massAxis, ptAxis, etaAxis, qtAxis, multAxis}, true);

    // track qa
    if (qa) {
      rpiKkink.add("h2_kinkradius_vs_pt", "kink radius_vs pt", {HistType::kTH2F, {radiusAxis, ptAxis}});
      rpiKkink.add("h2_kinkradius_vs_ncl", "kink radius_vs ncl", {HistType::kTH2F, {radiusAxis, {300, 0.0, 300.0}}});

      rpiKkink.add("h2_kinkradius_vs_vz", "kink radius_vz", {HistType::kTH2F, {vertexAxis, radiusAxis}});
      rpiKkink.add("h2_kink_vx_vs_vy", "kink vx vs vz ", {HistType::kTH2F, {vertexAxis, vertexAxis}});

      rpiKkink.add("tpc_dedx", "p vs dE/dx", {HistType::kTH2F, {{500, 0.0, 10.0}, {5000, 0.0, 5000.0}}});
      rpiKkink.add("tpc_nsigma_kaon", "p#k n#sigma", {HistType::kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}}});
      rpiKkink.add("tpc_nsigma_pion", "p#pi n#sigma", {HistType::kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}}});

      rpiKkink.add("tr_dcaxyM", "dcaxym", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
      rpiKkink.add("tr_dcaxyD", "dcaxyd", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
      rpiKkink.add("tr_dcaxykink_topo", "tr_dcaxykink_topo", {HistType::kTH1F, {{1000, -5.0, 5.0}}});

      rpiKkink.add("tr_chi2nclM", "chi2nclm", {HistType::kTH1F, {{100, 0.0, 100.0}}});
      rpiKkink.add("tr_chi2nclD", "chi2ncld", {HistType::kTH1F, {{100, 0.0, 100.0}}});
      rpiKkink.add("tr_tpcnclfindM", "tpcnclfindm", {HistType::kTH1F, {{300, 0.0, 300.0}}});
      rpiKkink.add("tr_tpcnclfindD", "tpcnclfindd", {HistType::kTH1F, {{300, 0.0, 300.0}}});
      rpiKkink.add("tr_itsChi2NClM", "itsChi2NClm", {HistType::kTH1F, {{200, 0.0, 200.0}}});

      rpiKkink.add("h2_kinkradius_vs_vz_m", "kink radius_vz", {HistType::kTH2F, {vertexAxis, radiusAxis}});
      rpiKkink.add("h2_kink_vx_vs_vy_m", "kink vx vs vz ", {HistType::kTH2F, {vertexAxis, vertexAxis}});

      rpiKkink.add("tpc_dedx_m", "p vs dE/dx", {HistType::kTH2F, {{500, 0.0, 10.0}, {5000, 0.0, 5000.0}}});
      rpiKkink.add("tpc_nsigma_kaon_m", "p#k n#sigma", {HistType::kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}}});
      rpiKkink.add("tpc_nsigma_pion_m", "p#pi n#sigma", {HistType::kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}}});

      rpiKkink.add("tr_dcaxyM_m", "dcaxym", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
      rpiKkink.add("tr_dcaxyD_m", "dcaxyd", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
      rpiKkink.add("tr_dcaxykink_topo_m", "tr_dcaxykink_topo", {HistType::kTH1F, {{1000, -5.0, 5.0}}});

      rpiKkink.add("tr_chi2nclM_m", "chi2nclm", {HistType::kTH1F, {{100, 0.0, 100.0}}});
      rpiKkink.add("tr_chi2nclD_m", "chi2ncld", {HistType::kTH1F, {{100, 0.0, 100.0}}});
      rpiKkink.add("tr_tpcnclfindM_m", "tpcnclfindm", {HistType::kTH1F, {{300, 0.0, 300.0}}});
      rpiKkink.add("tr_tpcnclfindD_m", "tpcnclfindd", {HistType::kTH1F, {{300, 0.0, 300.0}}});
      rpiKkink.add("tr_itsChi2NClM_m", "itsChi2NClm", {HistType::kTH1F, {{200, 0.0, 200.0}}});

      rpiKkink.add("h2_kinkradius_vs_pt_m", "kinkradius_vs_pt", {HistType::kTH2F, {{250, 0.0, 250.0}, ptAxis}});
      rpiKkink.add("h2_kinkradius_vs_ncl_m", "kinkradius_vs_ncl", {HistType::kTH2F, {{250, 0.0, 250.0}, {300, 0.0, 300.0}}});
    }
    if (doprocessMC) {
      rpiKkink.add("h2_dau_pt_vs_eta_gen", "pt_vs_eta_dau", {HistType::kTH2F, {ptAxis, etaAxis}});
      rpiKkink.add("h2_moth_pt_vs_eta_gen", "pt_vs_eta_moth", {HistType::kTH2F, {ptAxis, etaAxis}});
      rpiKkink.add("h2_pt_moth_vs_dau_gen", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});

      rpiKkink.add("h2_qt_gen", "qt", {HistType::kTH1F, {qtAxis}});
      rpiKkink.add("h2_qt_rec", "qt", {HistType::kTH1F, {qtAxis}});
      rpiKkink.add("h2_kink_angle_gen", "kink angle", {HistType::kTH1F, {kinkAxis}});

      rpiKkink.add("h2_kaon_mc_gen", "h2_kaon_mc_gen", HistType::kTHnSparseF, {massAxis, ptAxis, etaAxis, qtAxis}, true);
      rpiKkink.add("h2_kaon_mc_rec", "h2_kaon_mc_rec", HistType::kTHnSparseF, {massAxis, ptAxis, etaAxis, qtAxis, multAxis}, true);

      rpiKkink.add("h2_dau_pt_vs_eta_rec_m", "pt_vs_eta_dau", {HistType::kTH3F, {ptAxis, etaAxis, multAxis}});
      rpiKkink.add("h2_moth_pt_vs_eta_rec_m", "pt_vs_eta_moth", {HistType::kTH3F, {ptAxis, etaAxis, multAxis}});
      rpiKkink.add("h2_pt_moth_vs_dau_rec_m", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});
      rpiKkink.add("h2_kink_angle_m", "kink angle", {HistType::kTH2F, {kinkAxis, multAxis}});
      rpiKkink.add("h2_invmass_kaon_m", "Inv mass vs Pt", {HistType::kTH3F, {massAxis, ptAxis, ptAxis}});
      rpiKkink.add("h2_kaon_mc_rec_m", "h2_kaon_mc_rec_m", HistType::kTHnSparseF, {massAxis, ptAxis, etaAxis, qtAxis, multAxis}, true);
    }
  }

  double computeMotherMass(ROOT::Math::PxPyPzMVector p_moth, ROOT::Math::PxPyPzMVector p_daug)
  {
    // Infer neutrino momentum from conservation
    ROOT::Math::XYZVector p_nu_vec = p_moth.Vect() - p_daug.Vect();

    // Neutrino energy (massless): E_nu = |p_nu|
    double E_nu = p_nu_vec.R();

    // Total energy of the system
    double E_total = p_daug.E() + E_nu;

    // Total momentum = p_nu + p_daug
    ROOT::Math::XYZVector p_total_vec = p_nu_vec + p_daug.Vect();
    double p_total_sq = p_total_vec.Mag2();

    // Invariant mass from E² - |p|²
    double m2 = E_total * E_total - p_total_sq;
    return (m2 > 0) ? std::sqrt(m2) : -1.0;
  }

  void processData(CollisionsFull::iterator const& collision, aod::KinkCands const& KinkCands, TracksFull const&)
  {
    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector v1;

    if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8() || !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }

    float multiplicity{-1};
    const int kCentFT0C = 0;
    const int kCentFT0A = 1;
    const int kCentFT0M = 2;
    const int kCentFV0A = 3;

    if (centestimator == kCentFT0C) {
      multiplicity = collision.centFT0C();
    } else if (centestimator == kCentFT0A) {
      multiplicity = collision.centFT0A();
    } else if (centestimator == kCentFT0M) {
      multiplicity = collision.centFT0M();
    } else if (centestimator == kCentFV0A) {
      multiplicity = collision.centFV0A();
    }

    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hMultiplicity"), multiplicity);
    for (const auto& kinkCand : KinkCands) {
      auto dauTrack = kinkCand.trackDaug_as<TracksFull>();
      auto mothTrack = kinkCand.trackMoth_as<TracksFull>();
      if (mothTrack.collisionId() != collision.globalIndex()) {
        continue; // not from this event
      }
      if (!mothTrack.has_collision() || !dauTrack.has_collision()) {
        continue;
      }
      if (mothTrack.collisionId() != dauTrack.collisionId()) {
        continue; // skip mismatched collision tracks
      }
      if (dauTrack.sign() != mothTrack.sign()) {
        LOG(info) << "Skipping kink candidate with opposite sign daughter and mother: " << kinkCand.globalIndex();
        continue; // Skip if the daughter has the opposite sign as the mother
      }
      bool kaon = false;
      bool pion = false;

      v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassKaonCharged);
      v1.SetCoordinates(dauTrack.px(), dauTrack.py(), dauTrack.pz(), o2::constants::physics::MassMuon);

      if (dptCut && v1.Pt() > v0.Pt())
        continue;
      if (qa) {

        rpiKkink.fill(HIST("tpc_dedx"), v0.P(), mothTrack.tpcSignal());
        rpiKkink.fill(HIST("tpc_nsigma_kaon"), v0.Pt(), mothTrack.tpcNSigmaKa());
        rpiKkink.fill(HIST("tpc_nsigma_pion"), v0.Pt(), mothTrack.tpcNSigmaPi());

        rpiKkink.fill(HIST("tr_chi2nclM"), mothTrack.tpcChi2NCl());
        rpiKkink.fill(HIST("tr_chi2nclD"), dauTrack.tpcChi2NCl());
        rpiKkink.fill(HIST("tr_tpcnclfindM"), mothTrack.tpcNClsFound());
        rpiKkink.fill(HIST("tr_tpcnclfindD"), dauTrack.tpcNClsFound());
        rpiKkink.fill(HIST("tr_itsChi2NClM"), mothTrack.itsChi2NCl());
      }
      if (mothTrack.tpcChi2NCl() > tpcChi2Cut)
        continue;

      if (mothTrack.tpcNClsFound() > maxtpcncle || mothTrack.tpcNClsFound() < mintpcncle)
        continue;

      if (std::abs(mothTrack.tpcNSigmaKa()) < cutNSigmaKa) {
        kaon = true;
      }
      if (std::abs(mothTrack.tpcNSigmaPi()) < cutNSigmaPi) {
        pion = true;
      }
      if (!kaon && !pion) {
        continue;
      }
      if (cutNSigmaMu != -1 && std::abs(dauTrack.tpcNSigmaMu()) > cutNSigmaMu) {
        continue;
      }
      double radiusxy = std::sqrt(kinkCand.xDecVtx() * kinkCand.xDecVtx() + kinkCand.yDecVtx() * kinkCand.yDecVtx());
      if (radiusxy < minradius || radiusxy > maxradius)
        continue;
      // dcaXYmoth        kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
      if (qa) {
        rpiKkink.fill(HIST("tr_dcaxyM"), kinkCand.dcaMothPv());
        rpiKkink.fill(HIST("tr_dcaxyD"), kinkCand.dcaDaugPv());
        rpiKkink.fill(HIST("tr_dcaxykink_topo"), kinkCand.dcaKinkTopo());

        rpiKkink.fill(HIST("h2_kinkradius_vs_vz"), kinkCand.zDecVtx(), radiusxy);
        rpiKkink.fill(HIST("h2_kink_vx_vs_vy"), kinkCand.xDecVtx(), kinkCand.yDecVtx());
        rpiKkink.fill(HIST("h2_kinkradius_vs_pt"), radiusxy, v0.Pt());
        rpiKkink.fill(HIST("h2_kinkradius_vs_ncl"), radiusxy, mothTrack.tpcNClsFound());
      }
      if (std::abs(kinkCand.dcaMothPv()) > dcaXYcut)
        continue;

      if (kinkCand.dcaKinkTopo() > dcaXYcutkink)
        continue;

      float pMoth = v0.P();
      float pDaug = v1.P();
      float spKink = mothTrack.px() * dauTrack.px() + mothTrack.py() * dauTrack.py() + mothTrack.pz() * dauTrack.pz();
      float kinkangle = std::acos(spKink / (pMoth * pDaug));
      float radToDeg = o2::constants::math::Rad2Deg;
      if (kinkangle * radToDeg < kinkanglecut)
        continue;
      TVector3 pdlab(v1.Px(), v1.Py(), v1.Pz());
      // Compute transverse component
      TVector3 motherDir(v0.Px(), v0.Py(), v0.Pz());
      double ptd = pdlab.Perp(motherDir); // or p_d_lab.Mag() * sin(theta)

      if (kaon && onlykaon && std::abs(v0.Rapidity()) < rapCut) {
        v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassKaonCharged);
        if (additionalhist) {
          rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec"), v0.Pt(), v0.Eta(), multiplicity);
          rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec"), v1.Pt(), v1.Eta(), multiplicity);
          rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec"), v0.Pt(), v1.Pt());
          rpiKkink.fill(HIST("h2_kink_angle"), kinkangle, multiplicity);
          rpiKkink.fill(HIST("h2_qt_vs_pt"), ptd, v1.Pt());
        }
        double mass = computeMotherMass(v0, v1);
        rpiKkink.fill(HIST("h2_kaon_data"), mass, v0.Pt(), v0.Rapidity(), ptd, multiplicity);

        if (qa && ptd > minqt && ptd < maxqt) {
          rpiKkink.fill(HIST("h2_kinkradius_vs_vz_m"), kinkCand.zDecVtx(), radiusxy);
          rpiKkink.fill(HIST("h2_kink_vx_vs_vy_m"), kinkCand.xDecVtx(), kinkCand.yDecVtx());

          rpiKkink.fill(HIST("tpc_dedx_m"), v0.P(), mothTrack.tpcSignal());
          rpiKkink.fill(HIST("tpc_nsigma_kaon_m"), v0.Pt(), mothTrack.tpcNSigmaKa());
          rpiKkink.fill(HIST("tpc_nsigma_pion_m"), v0.Pt(), mothTrack.tpcNSigmaPi());

          rpiKkink.fill(HIST("tr_chi2nclM_m"), mothTrack.tpcChi2NCl());
          rpiKkink.fill(HIST("tr_chi2nclD_m"), dauTrack.tpcChi2NCl());
          rpiKkink.fill(HIST("tr_tpcnclfindM_m"), mothTrack.tpcNClsFound());
          rpiKkink.fill(HIST("tr_tpcnclfindD_m"), dauTrack.tpcNClsFound());
          rpiKkink.fill(HIST("tr_itsChi2NClM_m"), mothTrack.itsChi2NCl());

          rpiKkink.fill(HIST("tr_dcaxyM_m"), kinkCand.dcaMothPv());
          rpiKkink.fill(HIST("tr_dcaxyD_m"), kinkCand.dcaDaugPv());
          rpiKkink.fill(HIST("tr_dcaxykink_topo_m"), kinkCand.dcaKinkTopo());

          rpiKkink.fill(HIST("h2_kinkradius_vs_pt"), radiusxy, v0.Pt());
          rpiKkink.fill(HIST("h2_kinkradius_vs_ncl"), radiusxy, mothTrack.tpcNClsFound());
        }
      }
      if (pion && onlypion && std::abs(v0.Rapidity()) < rapCut) {
        v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassPionCharged);
        if (additionalhist) {
          rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec_pion"), v0.Pt(), v0.Eta(), multiplicity);
          rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec_pion"), v1.Pt(), v1.Eta(), multiplicity);
          rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec_pion"), v0.Pt(), v1.Pt());
          rpiKkink.fill(HIST("h2_kink_angle_pion"), kinkangle, multiplicity);
          rpiKkink.fill(HIST("h2_qt_vs_ptpion"), ptd, v1.Pt());
        }
        double mass = computeMotherMass(v0, v1);
        rpiKkink.fill(HIST("h2_pion_data"), mass, v0.Pt(), v0.Rapidity(), ptd, multiplicity);
      }
    }
  }
  PROCESS_SWITCH(spectraKinkPiKa, processData, "Data processing", true);

  void processMC(CollisionsFullMC const& collisions, aod::KinkCands const& KinkCands, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC, TracksFull const&)
  {
    for (const auto& collision : collisions) {
      ROOT::Math::PxPyPzMVector v0;
      ROOT::Math::PxPyPzMVector v1;
      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8() || !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
        continue;
      }
      if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        continue;
      }

      float multiplicity{-1};
      const int kCentFT0C = 0;
      const int kCentFT0A = 1;
      const int kCentFT0M = 2;
      const int kCentFV0A = 3;

      if (centestimator == kCentFT0C) {
        multiplicity = collision.centFT0C();
      } else if (centestimator == kCentFT0A) {
        multiplicity = collision.centFT0A();
      } else if (centestimator == kCentFT0M) {
        multiplicity = collision.centFT0M();
      } else if (centestimator == kCentFV0A) {
        multiplicity = collision.centFV0A();
      }
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hMultiplicity"), multiplicity);

      auto kinkCandPerColl = KinkCands.sliceBy(mPerCol, collision.globalIndex());
      for (const auto& kinkCand : kinkCandPerColl) {
        auto dauTrack = kinkCand.trackDaug_as<TracksFull>();
        auto mothTrack = kinkCand.trackMoth_as<TracksFull>();
        if (dauTrack.sign() != mothTrack.sign()) {
          LOG(info) << "Skipping kink candidate with opposite sign daughter and mother: " << kinkCand.globalIndex();
          continue; // Skip if the daughter has the opposite sign as the mother
        }

        v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetCoordinates(dauTrack.px(), dauTrack.py(), dauTrack.pz(), o2::constants::physics::MassMuon);
        if (qa) {
          rpiKkink.fill(HIST("tpc_dedx"), v0.P(), mothTrack.tpcSignal());
          rpiKkink.fill(HIST("tpc_nsigma_kaon"), v0.Pt(), mothTrack.tpcNSigmaKa());
          rpiKkink.fill(HIST("tpc_nsigma_pion"), v0.Pt(), mothTrack.tpcNSigmaPi());

          rpiKkink.fill(HIST("tr_chi2nclM"), mothTrack.tpcChi2NCl());
          rpiKkink.fill(HIST("tr_chi2nclD"), dauTrack.tpcChi2NCl());
          rpiKkink.fill(HIST("tr_tpcnclfindM"), mothTrack.tpcNClsFound());
          rpiKkink.fill(HIST("tr_tpcnclfindD"), dauTrack.tpcNClsFound());
          rpiKkink.fill(HIST("tr_itsChi2NClM"), mothTrack.itsChi2NCl());
        }
        if (mothTrack.tpcChi2NCl() > tpcChi2Cut)
          continue;

        if (mothTrack.tpcNClsFound() > maxtpcncle || mothTrack.tpcNClsFound() < mintpcncle)
          continue;

        bool kaon = false;
        if (std::abs(mothTrack.tpcNSigmaKa()) < cutNSigmaKa) {
          kaon = true;
        }

        if (dptCut && v1.Pt() > v0.Pt())
          continue;

        if (!kaon)
          continue;
        double radiusxy = std::sqrt(kinkCand.xDecVtx() * kinkCand.xDecVtx() + kinkCand.yDecVtx() * kinkCand.yDecVtx());
        if (radiusxy < minradius || radiusxy > maxradius)
          continue;
        if (qa) {
          rpiKkink.fill(HIST("tr_dcaxyM"), kinkCand.dcaMothPv());
          rpiKkink.fill(HIST("tr_dcaxyD"), kinkCand.dcaDaugPv());
          rpiKkink.fill(HIST("tr_dcaxykink_topo"), kinkCand.dcaKinkTopo());

          rpiKkink.fill(HIST("h2_kinkradius_vs_vz"), kinkCand.zDecVtx(), radiusxy);
          rpiKkink.fill(HIST("h2_kink_vx_vs_vy"), kinkCand.xDecVtx(), kinkCand.yDecVtx());

          rpiKkink.fill(HIST("h2_kinkradius_vs_pt"), radiusxy, v0.Pt());
          rpiKkink.fill(HIST("h2_kinkradius_vs_ncl"), radiusxy, mothTrack.tpcNClsFound());
        }
        if (std::abs(kinkCand.dcaMothPv()) > dcaXYcut)
          continue;

        if (kinkCand.dcaKinkTopo() > dcaXYcutkink)
          continue;

        float pMoth = v0.P();
        float pDaug = v1.P();
        float spKink = mothTrack.px() * dauTrack.px() + mothTrack.py() * dauTrack.py() + mothTrack.pz() * dauTrack.pz();
        float kinkangle = std::acos(spKink / (pMoth * pDaug));
        float radToDeg = o2::constants::math::Rad2Deg;
        if (kinkangle * radToDeg < kinkanglecut)
          continue;
        if (additionalhist) {
          rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec"), v0.Pt(), v0.Eta(), multiplicity);
          rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec"), v1.Pt(), v1.Eta(), multiplicity);
          rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec"), v0.Pt(), v1.Pt());
          rpiKkink.fill(HIST("h2_kink_angle"), kinkangle, multiplicity);
        }
        TVector3 pdlab(v1.Px(), v1.Py(), v1.Pz());
        // Compute transverse component
        TVector3 motherDir(v0.Px(), v0.Py(), v0.Pz());
        double ptd = pdlab.Perp(motherDir); // or p_d_lab.Mag() * sin(theta)

        rpiKkink.fill(HIST("h2_qt"), ptd);

        double mass = computeMotherMass(v0, v1);
        rpiKkink.fill(HIST("h2_invmass_kaon"), mass, v0.Pt(), ptd);

        rpiKkink.fill(HIST("h2_kaon_mc_rec"), mass, v0.Pt(), v0.Eta(), ptd, multiplicity);
        // do MC association
        auto mcLabMoth = trackLabelsMC.rawIteratorAt(mothTrack.globalIndex());
        auto mcLabDau = trackLabelsMC.rawIteratorAt(dauTrack.globalIndex());

        if (mcLabMoth.has_mcParticle() && mcLabDau.has_mcParticle()) {

          auto mcTrackMoth = mcLabMoth.mcParticle_as<aod::McParticles>();
          auto mcTrackDau = mcLabDau.mcParticle_as<aod::McParticles>();

          if (!mcTrackDau.has_mothers()) {
            continue;
          }
          for (const auto& piMother : mcTrackDau.mothers_as<aod::McParticles>()) {
            if (piMother.globalIndex() != mcTrackMoth.globalIndex()) {
              continue;
            }
            if (std::abs(mcTrackMoth.pdgCode()) != pid || std::abs(mcTrackDau.pdgCode()) != dpid) {
              continue;
            }
            if (additionalhist) {
              rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec_m"), v0.Pt(), v0.Eta(), multiplicity);
              rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec_m"), v1.Pt(), v1.Eta(), multiplicity);
              rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec_m"), v0.Pt(), v1.Pt());
              rpiKkink.fill(HIST("h2_kink_angle_m"), kinkangle, multiplicity);
            }
            rpiKkink.fill(HIST("h2_kaon_mc_rec_m"), mass, v0.Pt(), v0.Rapidity(), ptd, multiplicity);

            if (qa && ptd > minqt && ptd < maxqt) {
              rpiKkink.fill(HIST("h2_kinkradius_vs_vz_m"), kinkCand.zDecVtx(), radiusxy);
              rpiKkink.fill(HIST("h2_kink_vx_vs_vy_m"), kinkCand.xDecVtx(), kinkCand.yDecVtx());

              rpiKkink.fill(HIST("tpc_dedx_m"), v0.P(), mothTrack.tpcSignal());
              rpiKkink.fill(HIST("tpc_nsigma_kaon_m"), v0.Pt(), mothTrack.tpcNSigmaKa());
              rpiKkink.fill(HIST("tpc_nsigma_pion_m"), v0.Pt(), mothTrack.tpcNSigmaPi());

              rpiKkink.fill(HIST("tr_chi2nclM_m"), mothTrack.tpcChi2NCl());
              rpiKkink.fill(HIST("tr_chi2nclD_m"), dauTrack.tpcChi2NCl());
              rpiKkink.fill(HIST("tr_tpcnclfindM_m"), mothTrack.tpcNClsFound());
              rpiKkink.fill(HIST("tr_tpcnclfindD_m"), dauTrack.tpcNClsFound());
              rpiKkink.fill(HIST("tr_itsChi2NClM_m"), mothTrack.itsChi2NCl());

              rpiKkink.fill(HIST("tr_dcaxyM_m"), kinkCand.dcaMothPv());
              rpiKkink.fill(HIST("tr_dcaxyD_m"), kinkCand.dcaDaugPv());
              rpiKkink.fill(HIST("tr_dcaxykink_topo_m"), kinkCand.dcaKinkTopo());

              rpiKkink.fill(HIST("h2_kinkradius_vs_pt_m"), radiusxy, v0.Pt());
              rpiKkink.fill(HIST("h2_kinkradius_vs_ncl_m"), radiusxy, mothTrack.tpcNClsFound());
            }
          }
        }
      }
    }

    for (const auto& mcPart : particlesMC) {
      ROOT::Math::PxPyPzMVector v0;
      ROOT::Math::PxPyPzMVector v1;

      if (std::abs(mcPart.pdgCode()) != pid) {
        continue;
      }

      if (!mcPart.has_daughters()) {
        continue; // Skip if no daughters
      }
      bool hasKaonpionDaughter = false;
      for (const auto& daughter : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(daughter.pdgCode()) == dpid) { // muon PDG code
          hasKaonpionDaughter = true;
          v1.SetCoordinates(daughter.px(), daughter.py(), daughter.pz(), o2::constants::physics::MassMuon);
          break; // Found a muon daughter, exit loop
        }
      }
      if (!hasKaonpionDaughter) {
        continue; // Skip if no muon daughter found
      }
      if (pid == kKPlus) {
        v0.SetCoordinates(mcPart.px(), mcPart.py(), mcPart.pz(), o2::constants::physics::MassKaonCharged);
      }

      if (pid == kPiPlus) {
        v0.SetCoordinates(mcPart.px(), mcPart.py(), mcPart.pz(), o2::constants::physics::MassPionCharged);
      }

      float pMoth = v0.P();
      float pDaug = v1.P();
      float spKink = v0.Px() * v1.Px() + v0.Py() * v1.Py() + v0.Pz() * v1.Pz();
      float kinkangle = std::acos(spKink / (pMoth * pDaug));
      float radToDeg = o2::constants::math::Rad2Deg;
      if (kinkangle * radToDeg < kinkanglecut)
        continue;
      if (additionalhist) {
        rpiKkink.fill(HIST("h2_moth_pt_vs_eta_gen"), v0.Pt(), v0.Eta());
        rpiKkink.fill(HIST("h2_dau_pt_vs_eta_gen"), v1.Pt(), v1.Eta());
        rpiKkink.fill(HIST("h2_pt_moth_vs_dau_gen"), v0.Pt(), v1.Pt());
        rpiKkink.fill(HIST("h2_kink_angle_gen"), kinkangle);
      }
      TVector3 pdlab(v1.Px(), v1.Py(), v1.Pz());
      // Compute transverse component
      TVector3 motherDir(v0.Px(), v0.Py(), v0.Pz());
      double ptd = pdlab.Perp(motherDir); // or p_d_lab.Mag() * sin(theta)
      rpiKkink.fill(HIST("h2_qt_gen"), ptd);
      double mass = computeMotherMass(v0, v1);
      rpiKkink.fill(HIST("h2_kaon_mc_gen"), mass, v0.Pt(), v0.Rapidity(), ptd);
    }
  }
  PROCESS_SWITCH(spectraKinkPiKa, processMC, "MC processing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto builderTask = adaptAnalysisTask<kinkBuilder>(cfgc);
  auto spectraTask = adaptAnalysisTask<spectraKinkPiKa>(cfgc);

  return {builderTask, spectraTask}; // Just return both tasks
}
