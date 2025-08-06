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
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSel>;
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
  Configurable<float> maxDCAMothToPV{"maxDCAMothToPV", 0.1, "Max DCA of the mother to the PV"};
  Configurable<float> minDCADaugToPV{"minDCADaugToPV", 0., "Min DCA of the daughter to the PV"};
  Configurable<float> minPtMoth{"minPtMoth", 0.5, "Minimum pT of the hypercandidate"};
  Configurable<float> maxZDiff{"maxZDiff", 20., "Max z difference between the kink daughter and the mother"};
  Configurable<float> maxPhiDiff{"maxPhiDiff", 100, "Max phi difference between the kink daughter and the mother"};
  Configurable<float> timeMarginNS{"timeMarginNS", 600, "Additional time res tolerance in ns"};
  Configurable<float> etaMax{"etaMax", 1., "eta daughter"};
  Configurable<float> nTPCClusMinDaug{"nTPCClusMinDaug", 30, "mother NTPC clusters cut"};
  Configurable<float> itsChi2cut{"itsChi2cut", 30, "mother itsChi2 cut"};
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
  Configurable<float> rapCut{"rapCut", 0.8, "rapCut"};
  Configurable<float> kinkanglecut{"kinkanglecut", 2.0, "kinkanglecut"};
  Configurable<float> minradius{"minradius", 130.0, "minradiuscut"};
  Configurable<float> maxradius{"maxradius", 200.0, "maxradiuscut"};
  Configurable<float> dcaXYcut{"dcaXYcut", 0.2, "dcaXYcut"};
  Configurable<float> dcaZcut{"dcaZcut", 0.2, "dcaZcut"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4.0, "tpcChi2Cut"};

  Configurable<int> pid{"pidMother", 321, ""};
  Configurable<int> dpid{"pidDaughter", 13, ""};
  Configurable<bool> d0pid{"dopid", 0, ""};

  Preslice<aod::KinkCands> mPerCol = aod::track::collisionId;
  Preslice<aod::Tracks> mtPerCol = aod::track::collisionId;

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec ptAxis{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec qtAxis{2000, 0, 2, "#it{q}_{T} (GeV/#it{c})"};
    const AxisSpec kinkAxis{200, 0, 4, "#theta"};
    const AxisSpec etaAxis{200, -5.0, 5.0, "#eta"};
    const AxisSpec vertexAxis{1200, -300., 300., "vrtx [cm]"};
    const AxisSpec radiusAxis{600, 0., 300., "vrtx [cm]"};
    const AxisSpec massAxis{600, 0.1, 0.7, "Inv mass (GeV/#it{c}^{2})"};

    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexAxis}});

    rpiKkink.add("h2_dau_pt_vs_eta_rec", "pt_vs_eta_dau", {HistType::kTH2F, {ptAxis, etaAxis}});
    rpiKkink.add("h2_moth_pt_vs_eta_rec", "pt_vs_eta_moth", {HistType::kTH2F, {ptAxis, etaAxis}});
    rpiKkink.add("h2_pt_moth_vs_dau_rec", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});

    rpiKkink.add("h2_qt", "qt", {HistType::kTH1F, {qtAxis}});
    rpiKkink.add("h2_qt_vs_pt", "qt_pt", {HistType::kTH2F, {qtAxis, ptAxis}});

    rpiKkink.add("h2_kink_angle", "kink angle", {HistType::kTH1F, {kinkAxis}});

    // pion
    rpiKkink.add("h2_dau_pt_vs_eta_rec_pion", "pt_vs_eta_dau", {HistType::kTH2F, {ptAxis, etaAxis}});
    rpiKkink.add("h2_moth_pt_vs_eta_rec_pion", "pt_vs_eta_moth", {HistType::kTH2F, {ptAxis, etaAxis}});
    rpiKkink.add("h2_pt_moth_vs_dau_rec_pion", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});

    // inv mass
    rpiKkink.add("h2_invmass_kaon", "Inv mass vs Pt", {HistType::kTH3F, {massAxis, ptAxis, ptAxis}});
    rpiKkink.add("h2_invmass_pion", "Inv mass vs Pt", {HistType::kTH3F, {massAxis, ptAxis, ptAxis}});

    rpiKkink.add("h2_qt_pion", "qt", {HistType::kTH1F, {qtAxis}});
    rpiKkink.add("h2_qt_vs_ptpion", "qt_pt", {HistType::kTH2F, {qtAxis, ptAxis}});
    rpiKkink.add("h2_kink_angle_pion", "kink angle", {HistType::kTH1F, {kinkAxis}});

    // track qa

    rpiKkink.add("h2_kinkradius_vs_vz", "kink radius_vz", {HistType::kTH2F, {vertexAxis, radiusAxis}});
    rpiKkink.add("h2_kink_vx_vs_vy", "kink vx vs vz ", {HistType::kTH2F, {vertexAxis, vertexAxis}});

    if (doprocessMC) {
      rpiKkink.add("h2_dau_pt_vs_eta_gen", "pt_vs_eta_dau", {HistType::kTH2F, {ptAxis, etaAxis}});
      rpiKkink.add("h2_moth_pt_vs_eta_gen", "pt_vs_eta_moth", {HistType::kTH2F, {ptAxis, etaAxis}});
      rpiKkink.add("h2_pt_moth_vs_dau_gen", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});

      rpiKkink.add("h2_qt_gen", "qt", {HistType::kTH1F, {qtAxis}});
      rpiKkink.add("h2_qt_rec", "qt", {HistType::kTH1F, {qtAxis}});
      rpiKkink.add("h2_kink_angle_gen", "kink angle", {HistType::kTH1F, {kinkAxis}});
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
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
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
      bool kaon = false;
      bool pion = false;
      /*
      if (mothTrack.dcaXY() > dcaXYcut)
  continue;
      if (mothTrack.dcaZ() > dcaZcut)
      continue;
      */
      if (mothTrack.tpcChi2NCl() > tpcChi2Cut)
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
      rpiKkink.fill(HIST("h2_kinkradius_vs_vz"), kinkCand.zDecVtx(), radiusxy);
      rpiKkink.fill(HIST("h2_kink_vx_vs_vy"), kinkCand.xDecVtx(), kinkCand.yDecVtx());

      v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassPionCharged);
      v1.SetCoordinates(dauTrack.px(), dauTrack.py(), dauTrack.pz(), o2::constants::physics::MassMuon);

      float pMoth = v0.P();
      float pDaug = v1.P();
      float spKink = mothTrack.px() * dauTrack.px() + mothTrack.py() * dauTrack.py() + mothTrack.pz() * dauTrack.pz();
      float kinkangle = std::acos(spKink / (pMoth * pDaug));
      float radToDeg = o2::constants::math::Rad2Deg;
      if (kinkangle * radToDeg < kinkanglecut)
        continue;

      if (kaon) {
        rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec"), v0.Pt(), v0.Eta());
        rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec"), v1.Pt(), v1.Eta());
        rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec"), v0.Pt(), v1.Pt());
        rpiKkink.fill(HIST("h2_kink_angle"), kinkangle);
      }
      if (pion) {
        rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec_pion"), v0.Pt(), v0.Eta());
        rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec_pion"), v1.Pt(), v1.Eta());
        rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec_pion"), v0.Pt(), v1.Pt());
        rpiKkink.fill(HIST("h2_kink_angle_pion"), kinkangle);
      }
      TVector3 pdlab(v1.Px(), v1.Py(), v1.Pz());
      // Compute transverse component
      TVector3 motherDir(v0.Px(), v0.Py(), v0.Pz());
      double ptd = pdlab.Perp(motherDir); // or p_d_lab.Mag() * sin(theta)

      if (kaon) {
        v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassKaonCharged);
        double mass = computeMotherMass(v0, v1);
        rpiKkink.fill(HIST("h2_qt"), ptd);
        rpiKkink.fill(HIST("h2_qt_vs_pt"), ptd, v1.Pt());
        rpiKkink.fill(HIST("h2_invmass_kaon"), mass, v0.Pt(), ptd);
      }
      if (pion) {
        v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassPionCharged);
        double mass = computeMotherMass(v0, v1);
        rpiKkink.fill(HIST("h2_qt_pion"), ptd);
        rpiKkink.fill(HIST("h2_qt_vs_ptpion"), ptd, v1.Pt());
        rpiKkink.fill(HIST("h2_invmass_pion"), mass, v0.Pt(), ptd);
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

      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      auto kinkCandPerColl = KinkCands.sliceBy(mPerCol, collision.globalIndex());
      for (const auto& kinkCand : kinkCandPerColl) {
        auto dauTrack = kinkCand.trackDaug_as<TracksFull>();
        auto mothTrack = kinkCand.trackMoth_as<TracksFull>();
        if (dauTrack.sign() != mothTrack.sign()) {
          LOG(info) << "Skipping kink candidate with opposite sign daughter and mother: " << kinkCand.globalIndex();
          continue; // Skip if the daughter has the opposite sign as the mother
        }
        bool kaon = false;
        bool pion = false;
        if (std::abs(mothTrack.tpcNSigmaKa()) < cutNSigmaKa) {
          kaon = true;
        }
        if (std::abs(mothTrack.tpcNSigmaPi()) < cutNSigmaPi) {
          pion = true;
        }
        if (!kaon && !pion)
          continue;
        double radiusxy = std::sqrt(kinkCand.xDecVtx() * kinkCand.xDecVtx() + kinkCand.yDecVtx() * kinkCand.yDecVtx());
        if (radiusxy < minradius || radiusxy > maxradius)
          continue;
        rpiKkink.fill(HIST("h2_kinkradius_vs_vz"), kinkCand.zDecVtx(), radiusxy);
        rpiKkink.fill(HIST("h2_kink_vx_vs_vy"), kinkCand.xDecVtx(), kinkCand.yDecVtx());
        v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassPionCharged);
        v1.SetCoordinates(dauTrack.px(), dauTrack.py(), dauTrack.pz(), o2::constants::physics::MassMuon);

        float pMoth = v0.P();
        float pDaug = v1.P();
        float spKink = mothTrack.px() * dauTrack.px() + mothTrack.py() * dauTrack.py() + mothTrack.pz() * dauTrack.pz();
        float kinkangle = std::acos(spKink / (pMoth * pDaug));
        float radToDeg = o2::constants::math::Rad2Deg;
        if (kinkangle * radToDeg < kinkanglecut)
          continue;

        rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec"), v0.Pt(), v0.Eta());
        rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec"), v1.Pt(), v1.Eta());
        rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec"), v0.Pt(), v1.Pt());
        rpiKkink.fill(HIST("h2_kink_angle"), kinkangle);

        TVector3 pdlab(v1.Px(), v1.Py(), v1.Pz());
        // Compute transverse component
        TVector3 motherDir(v0.Px(), v0.Py(), v0.Pz());
        double ptd = pdlab.Perp(motherDir); // or p_d_lab.Mag() * sin(theta)

        rpiKkink.fill(HIST("h2_qt"), ptd);

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
            rpiKkink.fill(HIST("h2_qt_rec"), ptd);
          }
        }
      }
    }
    for (const auto& mcPart : particlesMC) {
      ROOT::Math::PxPyPzMVector v0;
      ROOT::Math::PxPyPzMVector v1;
      if (!d0pid && (std::abs(mcPart.pdgCode()) != pid || std::abs(mcPart.eta()) > rapCut)) {
        continue;
      }
      bool isDmeson = std::abs(mcPart.pdgCode()) == kD0 || std::abs(mcPart.pdgCode()) == kDPlus || std::abs(mcPart.pdgCode()) == kDStar;
      if (d0pid && (!isDmeson || std::abs(mcPart.eta()) > rapCut)) {
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
      if (!d0pid && pid == kKPlus) {
        v0.SetCoordinates(mcPart.px(), mcPart.py(), mcPart.pz(), o2::constants::physics::MassKaonCharged);
      }

      if (!d0pid && pid == kPiPlus) {
        v0.SetCoordinates(mcPart.px(), mcPart.py(), mcPart.pz(), o2::constants::physics::MassPionCharged);
      }
      if (d0pid) {
        v0.SetCoordinates(mcPart.px(), mcPart.py(), mcPart.pz(), o2::constants::physics::MassD0);
      }
      float pMoth = v0.P();
      float pDaug = v1.P();
      float spKink = v0.Px() * v1.Px() + v0.Py() * v1.Py() + v0.Pz() * v1.Pz();
      float kinkangle = std::acos(spKink / (pMoth * pDaug));
      float radToDeg = o2::constants::math::Rad2Deg;
      if (kinkangle * radToDeg < kinkanglecut)
        continue;
      rpiKkink.fill(HIST("h2_moth_pt_vs_eta_gen"), v0.Pt(), v0.Eta());
      rpiKkink.fill(HIST("h2_dau_pt_vs_eta_gen"), v1.Pt(), v1.Eta());
      rpiKkink.fill(HIST("h2_pt_moth_vs_dau_gen"), v0.Pt(), v1.Pt());
      rpiKkink.fill(HIST("h2_kink_angle_gen"), kinkangle);
      TVector3 pdlab(v1.Px(), v1.Py(), v1.Pz());
      // Compute transverse component
      TVector3 motherDir(v0.Px(), v0.Py(), v0.Pz());
      double ptd = pdlab.Perp(motherDir); // or p_d_lab.Mag() * sin(theta)
      rpiKkink.fill(HIST("h2_qt_gen"), ptd);
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
