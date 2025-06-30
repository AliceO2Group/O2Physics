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
//
/// \file hyperhelium4sigmaRecoTask.cxx
/// \brief QA and analysis task for hyper-helium4sigma (He4S)
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#include "PWGLF/DataModel/LFHyperhelium4sigmaTables.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels>;
using MCLabeledCollisionsFull = soa::Join<CollisionsFull, aod::McCollisionLabels>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTPCFullAl, aod::pidTPCFullTr, aod::pidTPCFullPi>;
using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

enum Channel {
  k2body = 0, // helium4, pion0
  k3body_p,   // triton, proton, pion0
  k3body_n,   // triton, neutron, pion+
  kNDecayChannel
};

enum DaughterType {
  kDauAlpha = 0,
  kDauTriton,
  kDauProton,
  kDauChargedPion,
  kDauNeutron,
  kDauPion0,
  kNDaughterType,
  kNChargedDaughterType = kDauNeutron
};

namespace
{
constexpr std::array<float, 7> LayerRadii{2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
constexpr int kITSLayers = 7;
constexpr int kITSInnerBarrelLayers = 3;
// constexpr int kITSOuterBarrelLayers = 4;
constexpr std::array<int, kNDaughterType> kDaughterPDG = {
  o2::constants::physics::Pdg::kAlpha,
  o2::constants::physics::Pdg::kTriton,
  PDG_t::kProton,
  PDG_t::kPiPlus,
  PDG_t::kNeutron,
  PDG_t::kPi0};

std::shared_ptr<TH1> hMotherCounter;
std::shared_ptr<TH1> hMother2BCounter;
std::shared_ptr<TH1> hDauCounter[kNChargedDaughterType];
std::shared_ptr<TH2> hDauTPCNSigma[kNChargedDaughterType];
std::shared_ptr<TH1> hRecoMotherCounter;
std::shared_ptr<TH1> hRecoDauAlphaCounter;
} // namespace

//--------------------------------------------------------------
// Check the decay channel of hyperhelium4sigma
template <class TMCTrackTo, typename TMCParticle>
Channel getDecayChannelHe4S(TMCParticle const& particle, std::vector<int>& list)
{
  if (std::abs(particle.pdgCode()) != o2::constants::physics::Pdg::kHyperHelium4Sigma) {
    return kNDecayChannel;
  }

  // list: charged (alpha or triton), charged or empty, neutral
  list.clear();
  list.resize(3, -1);

  bool haveAlpha = false, haveTriton = false, haveProton = false, haveNeuteron = false;
  bool haveAntiAlpha = false, haveAntiTriton = false, haveAntiProton = false, haveAntiNeuteron = false;
  bool havePionPlus = false, havePionMinus = false, havePion0 = false;
  for (const auto& mcDaughter : particle.template daughters_as<TMCTrackTo>()) {
    if (mcDaughter.pdgCode() == o2::constants::physics::Pdg::kAlpha) {
      haveAlpha = true;
      list[0] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kAlpha) {
      haveAntiAlpha = true;
      list[0] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == o2::constants::physics::Pdg::kTriton) {
      haveTriton = true;
      list[0] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kTriton) {
      haveAntiTriton = true;
      list[0] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == PDG_t::kProton) {
      haveProton = true;
      list[1] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -PDG_t::kProton) {
      haveAntiProton = true;
      list[1] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == PDG_t::kNeutron) {
      haveNeuteron = true;
      list[2] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -PDG_t::kNeutron) {
      haveAntiNeuteron = true;
      list[2] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == PDG_t::kPiPlus) {
      havePionPlus = true;
      list[1] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -PDG_t::kPiPlus) {
      havePionMinus = true;
      list[1] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == PDG_t::kPi0) {
      havePion0 = true;
      list[2] = mcDaughter.globalIndex();
    }
  }

  if ((haveAlpha && havePion0) || (haveAntiAlpha && havePion0)) {
    return k2body;
  } else if ((haveTriton && haveProton && havePion0) || (haveAntiTriton && haveAntiProton && havePion0)) {
    return k3body_p;
  } else if ((haveTriton && haveNeuteron && havePionPlus) || (haveAntiTriton && haveAntiNeuteron && havePionMinus)) {
    return k3body_n;
  }

  return kNDecayChannel;
}

//--------------------------------------------------------------
// Extract track parameters from a mcparticle, use global coordinates as the local one
template <typename TrackPrecision = float, typename T>
o2::track::TrackParametrization<TrackPrecision> getTrackParFromMC(const T& mcparticle)
{
  int sign = mcparticle.pdgCode() > 0 ? 1 : -1; // ok for hyperhelium4sigma
  TrackPrecision snp = mcparticle.py() / (mcparticle.pt() + 1.e-10f);
  TrackPrecision tgl = mcparticle.pz() / (mcparticle.pt() + 1.e-10f);
  std::array<TrackPrecision, o2::track::kNParams> arraypar = {mcparticle.vy(), mcparticle.vz(), snp,
                                                              tgl, 2 * sign / (mcparticle.pt() + 1.e-10f)};
  return o2::track::TrackParametrization<TrackPrecision>(mcparticle.vx(), 0, std::move(arraypar));
}

//--------------------------------------------------------------
// construct index array from mcParticle to track
template <typename TTrackTable>
void setTrackIDForMC(std::vector<int64_t>& mcPartIndices, aod::McParticles const& particlesMC, TTrackTable const& tracks)
{
  mcPartIndices.clear();
  mcPartIndices.resize(particlesMC.size());
  std::fill(mcPartIndices.begin(), mcPartIndices.end(), -1);
  for (const auto& track : tracks) {
    if (track.has_mcParticle()) {
      auto mcparticle = track.template mcParticle_as<aod::McParticles>();
      if (mcPartIndices[mcparticle.globalIndex()] == -1) {
        mcPartIndices[mcparticle.globalIndex()] = track.globalIndex();
      } else {
        auto candTrack = tracks.rawIteratorAt(mcPartIndices[mcparticle.globalIndex()]);
        // Use the track which has innest information (also best quality?
        if (track.x() < candTrack.x()) {
          mcPartIndices[mcparticle.globalIndex()] = track.globalIndex();
        }
      }
    }
  }
}

//--------------------------------------------------------------
// get TPCNSigma for daughter track
template <typename TTrack>
float getTPCNSigma(const TTrack& track, const int daughterType)
{
  float nSigma = -999.f;
  switch (daughterType) {
    case kDauAlpha:
      nSigma = track.tpcNSigmaAl();
      break;
    case kDauTriton:
      nSigma = track.tpcNSigmaTr();
      break;
    case kDauProton:
      nSigma = track.tpcNSigmaPr();
      break;
    case kDauChargedPion:
      nSigma = track.tpcNSigmaPi();
      break;
    default:
      break;
  }
  return nSigma;
}

//--------------------------------------------------------------
struct Hyphe4sCandidate {

  bool isMatter = false;

  std::array<float, 3> primVtx = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> decVtx = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> lastPosMoth = {0.0f, 0.0f, 0.0f}; // last position of mother track at the radii of ITS layer which has the outermost update
  std::array<float, 3> momMoth = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> momDaug = {0.0f, 0.0f, 0.0f};

  float dcaXYMothPv = -999.f;
  float dcaXYDauPv = -999.f;
  float dcaKinkTopo = -999.f;

  float chi2ITSMoth = 0.0f;
  uint32_t itsClusterSizeMoth = 0u;
  uint32_t itsClusterSizeDau = 0u;
  float nSigmaTPCDau = -999.f;
  float nSigmaITSDau = -999.f;

  // mc information
  bool isSignal = false;
  bool isSignalReco = false;
  bool isCollReco = false;
  bool isSurvEvSelection = false;

  std::array<float, 3> trueDecVtx = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> gMomMoth = {0.0f, 0.0f, 0.0f};    // generated mother momentum
  std::array<float, 3> trueMomMoth = {0.0f, 0.0f, 0.0f}; // true mother momentum at decay vertex
  std::array<float, 3> gMomDau = {0.0f, 0.0f, 0.0f};

  bool isMothReco = false;
  float ptMoth = -999.f;
  float pzMoth = -999.f;
};

//--------------------------------------------------------------
// analysis task for hyperhelium4sigma 2-body decay
struct Hyperhelium4sigmaRecoTask {

  Produces<aod::He4S2BCands> outputDataTable;
  Produces<aod::MCHe4S2BCands> outputMCTable;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  std::vector<int> mcHe4sIndices;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry registry{"registry", {}};

  // Configurable for event selection
  Configurable<bool> doEventCut{"doEventCut", true, "Apply event selection"};
  Configurable<float> maxZVertex{"maxZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutNSigmaAl{"cutNSigmaAl", 5, "NSigmaTPCAlpha"};

  // CCDB options
  Configurable<double> inputBz{"inputBz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdbPath{"ccdbPath", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  int mRunNumber;
  float mBz;
  o2::base::MatLayerCylSet* lut = nullptr;

  o2::aod::ITSResponse itsResponse;

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec vertexZAxis{100, -15., 15., "vtx_{Z} [cm]"};
    const AxisSpec ptAxis{50, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nSigmaAxis{120, -6.f, 6.f, "n#sigma_{#alpha}"};
    const AxisSpec massAxis{100, 3.85, 4.25, "m (GeV/#it{c}^{2})"};
    const AxisSpec diffPxAxis{200, -10.f, 10.f, "#Delta #it{p}_{x} (GeV/#it{c})"};
    const AxisSpec diffPyAxis{200, -10.f, 10.f, "#Delta #it{p}_{y} (GeV/#it{c})"};
    const AxisSpec diffPtAxis{200, -10.f, 10.f, "#Delta #it{p}_{T} (GeV/#it{c})"};
    const AxisSpec diffPzAxis{200, -10.f, 10.f, "#Delta #it{p}_{z} (GeV/#it{c})"};
    const AxisSpec radiusAxis{40, 0.f, 40.f, "R (cm)"};

    registry.add<TH1>("hEventCounter", "hEventCounter", HistType::kTH1F, {{2, 0, 2}});
    registry.add<TH1>("hVertexZCollision", "hVertexZCollision", HistType::kTH1F, {vertexZAxis});
    registry.add<TH1>("hCandidateCounter", "hCandidateCounter", HistType::kTH1F, {{3, 0, 3}});

    if (doprocessMC == true) {
      itsResponse.setMCDefaultParameters();

      registry.add<TH1>("hTrueCandidateCounter", "hTrueCandidateCounter", HistType::kTH1F, {{3, 0, 3}});
      registry.add<TH1>("hDiffSVx", ";#Delta x (cm);", HistType::kTH1F, {{200, -10, 10}});
      registry.add<TH1>("hDiffSVy", ";#Delta y (cm);", HistType::kTH1F, {{200, -10, 10}});
      registry.add<TH1>("hDiffSVz", ";#Delta z (cm);", HistType::kTH1F, {{200, -10, 10}});
      registry.add<TH2>("h2RecSVRVsTrueSVR", ";Reconstruced SV R (cm);True SV R (cm);", HistType::kTH2F, {radiusAxis, radiusAxis});
      registry.add<TH2>("h2TrueMotherDiffPxVsRecSVR", ";Reconstruced SV R (cm);#Delta #it{p}_{T} (GeV/#it{c});", HistType::kTH2F, {radiusAxis, diffPxAxis});
      registry.add<TH2>("h2TrueMotherDiffPyVsRecSVR", ";Reconstruced SV R (cm);#Delta #it{p}_{T} (GeV/#it{c});", HistType::kTH2F, {radiusAxis, diffPyAxis});
      registry.add<TH2>("h2TrueMotherDiffPtVsRecSVR", ";Reconstruced SV R (cm);#Delta #it{p}_{T} (GeV/#it{c});", HistType::kTH2F, {radiusAxis, diffPtAxis});
      registry.add<TH2>("h2TrueMotherDiffPzVsRecSVR", ";Reconstruced SV R (cm);#Delta #it{p}_{z} (GeV/#it{c});", HistType::kTH2F, {radiusAxis, diffPzAxis});
      registry.add<TH2>("h2TrueMotherDiffTglVsRecSVR", ";Reconstruced SV R (cm);#Delta tan#lambda;", HistType::kTH2F, {radiusAxis, {200, -0.1, 0.1}});
      registry.add<TH2>("h2TrueMotherDiffEtaVsRecSVR", ";Reconstruced SV R (cm);#Delta #eta;", HistType::kTH2F, {radiusAxis, {200, -0.1, 0.1}});
      registry.add<TH1>("hDiffDauPx", ";#Delta p_{x} (GeV/#it{c}); ", HistType::kTH1D, {{200, -10, 10}});
      registry.add<TH1>("hDiffDauPy", ";#Delta p_{y} (GeV/#it{c}); ", HistType::kTH1D, {{200, -10, 10}});
      registry.add<TH1>("hDiffDauPz", ";#Delta p_{z} (GeV/#it{c}); ", HistType::kTH1D, {{200, -10, 10}});
      registry.add<TH2>("h2TrueSignalMassPt", "h2TrueSignalMassPt", HistType::kTH2F, {{ptAxis, massAxis}});
      registry.add<TH2>("h2TrueSignalNSigmaAlPt", "h2TrueSignalNSigmaAlPt", HistType::kTH2F, {{ptAxis, nSigmaAxis}});

      registry.add<TH1>("hDCAXYMothToRecSV", "hDCAXYMothToRecSV", HistType::kTH1F, {{200, -10, 10}});
      registry.add<TH1>("hDCAZMothToRecSV", "hDCAZMothToRecSV", HistType::kTH1F, {{200, -10, 10}});

      registry.add<TH2>("h2TrueMotherDiffPtVsKinkAngle", ";cos(#theta);#Delta p_{T} / p_{T};", HistType::kTH2F, {{100, 0.8f, 1.f}, {100, -5.f, 5.f}});
      registry.add<TH2>("h2TrueMotherDiffPtVsKinkAngleSV", ";cos(#theta);#Delta p_{T} / p_{T};", HistType::kTH2F, {{100, 0.8f, 1.f}, {100, -5.f, 5.f}});
      registry.add<TH2>("h2TrueMotherDiffPtVsKinkAngleXY", ";cos(#theta);#Delta p_{T} / p_{T};", HistType::kTH2F, {{100, 0.8f, 1.f}, {100, -5.f, 5.f}});
      registry.add<TH2>("h2TrueMotherDiffPtVsKinkAngleXYSV", ";cos(#theta);#Delta p_{T} / p_{T};", HistType::kTH2F, {{100, 0.8f, 1.f}, {100, -5.f, 5.f}});
    }

    registry.add<TH2>("h2MassHyperhelium4sigmaPt", "h2MassHyperhelium4sigmaPt", HistType::kTH2F, {{ptAxis, massAxis}});
    registry.add<TH2>("h2NSigmaAlPt", "h2NSigmaAlPt", HistType::kTH2F, {{ptAxis, nSigmaAxis}});

    ccdb->setURL(ccdbPath);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
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

    if (!lut) {
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }
    o2::base::Propagator::Instance()->setMatLUT(lut);
    LOG(info) << "Task initialized for run " << mRunNumber << " with magnetic field " << mBz << " kZG";
  }

  template <typename TCollision, typename TKindCandidate, typename TTrack>
  void fillCandidate(Hyphe4sCandidate& hyphe4sCand, TCollision const& collision, TKindCandidate const& kinkCand, TTrack const& trackMoth, TTrack const& trackDau)
  {
    hyphe4sCand.isMatter = kinkCand.mothSign() > 0;
    hyphe4sCand.primVtx[0] = collision.posX();
    hyphe4sCand.primVtx[1] = collision.posY();
    hyphe4sCand.primVtx[2] = collision.posZ();
    hyphe4sCand.decVtx[0] = kinkCand.xDecVtx() + collision.posX();
    hyphe4sCand.decVtx[1] = kinkCand.yDecVtx() + collision.posY();
    hyphe4sCand.decVtx[2] = kinkCand.zDecVtx() + collision.posZ();

    hyphe4sCand.momMoth[0] = kinkCand.pxMoth();
    hyphe4sCand.momMoth[1] = kinkCand.pyMoth();
    hyphe4sCand.momMoth[2] = kinkCand.pzMoth();
    hyphe4sCand.momDaug[0] = kinkCand.pxDaug();
    hyphe4sCand.momDaug[1] = kinkCand.pyDaug();
    hyphe4sCand.momDaug[2] = kinkCand.pzDaug();

    hyphe4sCand.dcaXYMothPv = kinkCand.dcaMothPv();
    hyphe4sCand.dcaXYDauPv = kinkCand.dcaDaugPv();
    hyphe4sCand.dcaKinkTopo = kinkCand.dcaKinkTopo();

    fillCandidateRecoMoth(hyphe4sCand, trackMoth);

    hyphe4sCand.itsClusterSizeDau = trackDau.itsClusterSizes();
    hyphe4sCand.nSigmaTPCDau = trackDau.tpcNSigmaAl();
    hyphe4sCand.nSigmaITSDau = itsResponse.nSigmaITS<o2::track::PID::Alpha>(trackDau);

    int lastLayerMoth = 0;
    for (int i = 6; i >= 0; i--) {
      if (trackMoth.itsClusterMap() & (1 << i)) {
        lastLayerMoth = i;
        break;
      }
    }
    auto trackparMother = getTrackParCov(trackMoth);
    o2::base::Propagator::Instance()->PropagateToXBxByBz(trackparMother, LayerRadii[lastLayerMoth]);
    std::array<float, 9> vecLab{0.f};
    if (trackparMother.getPosDirGlo(vecLab)) {
      hyphe4sCand.lastPosMoth[0] = vecLab[0];
      hyphe4sCand.lastPosMoth[1] = vecLab[1];
      hyphe4sCand.lastPosMoth[2] = vecLab[2];
    }
  }

  template <typename TTrack>
  void fillCandidateRecoMoth(Hyphe4sCandidate& hyphe4sCand, TTrack const& trackMoth)
  {
    hyphe4sCand.isMothReco = true;
    hyphe4sCand.chi2ITSMoth = trackMoth.itsChi2NCl();
    hyphe4sCand.itsClusterSizeMoth = trackMoth.itsClusterSizes();
    hyphe4sCand.ptMoth = trackMoth.pt();
    hyphe4sCand.pzMoth = trackMoth.pz();
  }

  template <typename TMCParticle>
  void fillCandidateMCInfo(Hyphe4sCandidate& hyphe4sCand, TMCParticle const& mcMothTrack, TMCParticle const& mcDauTrack, TMCParticle const& mcNeutDauTrack)
  {
    hyphe4sCand.trueDecVtx[0] = mcDauTrack.vx();
    hyphe4sCand.trueDecVtx[1] = mcDauTrack.vy();
    hyphe4sCand.trueDecVtx[2] = mcDauTrack.vz();
    hyphe4sCand.gMomMoth[0] = mcMothTrack.px();
    hyphe4sCand.gMomMoth[1] = mcMothTrack.py();
    hyphe4sCand.gMomMoth[2] = mcMothTrack.pz();
    hyphe4sCand.trueMomMoth[0] = mcDauTrack.px() + mcNeutDauTrack.px();
    hyphe4sCand.trueMomMoth[1] = mcDauTrack.py() + mcNeutDauTrack.py();
    hyphe4sCand.trueMomMoth[2] = mcDauTrack.pz() + mcNeutDauTrack.pz();
    hyphe4sCand.gMomDau[0] = mcDauTrack.px();
    hyphe4sCand.gMomDau[1] = mcDauTrack.py();
    hyphe4sCand.gMomDau[2] = mcDauTrack.pz();
  }

  void processData(CollisionsFull const& collisions, aod::KinkCands const& KinkCands, FullTracksExtIU const&, aod::BCs const&)
  {
    for (const auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 0);
      if (doEventCut && (std::abs(collision.posZ()) > maxZVertex || !collision.sel8())) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1);
      registry.fill(HIST("hVertexZCollision"), collision.posZ());
    }

    for (const auto& kinkCand : KinkCands) {
      registry.fill(HIST("hCandidateCounter"), 0);
      auto collision = kinkCand.collision_as<CollisionsFull>();
      if (doEventCut && (std::abs(collision.posZ()) > maxZVertex || !collision.sel8())) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 1);

      auto dauTrack = kinkCand.trackDaug_as<FullTracksExtIU>();
      if (std::abs(dauTrack.tpcNSigmaAl()) > cutNSigmaAl) {
        continue;
      }
      float invMass = RecoDecay::m(std::array{std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()}, std::array{kinkCand.pxDaugNeut(), kinkCand.pyDaugNeut(), kinkCand.pzDaugNeut()}}, std::array{o2::constants::physics::MassAlpha, o2::constants::physics::MassPi0});
      registry.fill(HIST("hCandidateCounter"), 2);
      registry.fill(HIST("h2MassHyperhelium4sigmaPt"), kinkCand.mothSign() * kinkCand.ptMoth(), invMass);
      registry.fill(HIST("h2NSigmaAlPt"), kinkCand.mothSign() * kinkCand.ptDaug(), dauTrack.tpcNSigmaAl());

      auto bc = collision.bc_as<aod::BCs>();
      initCCDB(bc);
      auto motherTrack = kinkCand.trackMoth_as<FullTracksExtIU>();
      Hyphe4sCandidate hyphe4sCand;
      fillCandidate(hyphe4sCand, collision, kinkCand, motherTrack, dauTrack);

      outputDataTable(
        hyphe4sCand.primVtx[0], hyphe4sCand.primVtx[1], hyphe4sCand.primVtx[2],
        hyphe4sCand.decVtx[0], hyphe4sCand.decVtx[1], hyphe4sCand.decVtx[2],
        hyphe4sCand.isMatter,
        hyphe4sCand.lastPosMoth[0], hyphe4sCand.lastPosMoth[1], hyphe4sCand.lastPosMoth[2],
        hyphe4sCand.momMoth[0], hyphe4sCand.momMoth[1], hyphe4sCand.momMoth[2],
        hyphe4sCand.momDaug[0], hyphe4sCand.momDaug[1], hyphe4sCand.momDaug[2],
        hyphe4sCand.dcaXYMothPv, hyphe4sCand.dcaXYDauPv, hyphe4sCand.dcaKinkTopo,
        hyphe4sCand.chi2ITSMoth, hyphe4sCand.itsClusterSizeMoth, hyphe4sCand.itsClusterSizeDau,
        hyphe4sCand.nSigmaTPCDau, hyphe4sCand.nSigmaITSDau);
    }
  }
  PROCESS_SWITCH(Hyperhelium4sigmaRecoTask, processData, "process data", true);

  void processMC(MCLabeledCollisionsFull const& collisions, aod::KinkCands const& KinkCands, MCLabeledTracksIU const& tracks, aod::McParticles const& particlesMC, aod::McCollisions const& mcCollisions, aod::BCs const&)
  {
    mcHe4sIndices.clear();
    std::vector<int64_t> mcPartIndices;
    setTrackIDForMC(mcPartIndices, particlesMC, tracks);
    std::vector<bool> isReconstructedMCCollisions(mcCollisions.size(), false);
    std::vector<bool> isSelectedMCCollisions(mcCollisions.size(), false);
    std::vector<bool> isGoodCollisions(collisions.size(), false);
    std::vector<int> dauIDList(3, -1);

    for (const auto& collision : collisions) {
      isReconstructedMCCollisions[collision.mcCollisionId()] = true;
      registry.fill(HIST("hEventCounter"), 0);
      if (doEventCut && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || std::abs(collision.posZ()) > maxZVertex)) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1);
      registry.fill(HIST("hVertexZCollision"), collision.posZ());
      isSelectedMCCollisions[collision.mcCollisionId()] = true;
      isGoodCollisions[collision.globalIndex()] = true;
    }

    for (const auto& kinkCand : KinkCands) {
      auto motherTrack = kinkCand.trackMoth_as<MCLabeledTracksIU>();
      auto dauTrack = kinkCand.trackDaug_as<MCLabeledTracksIU>();

      bool isTrueSignal = false;
      if (motherTrack.has_mcParticle() && dauTrack.has_mcParticle()) {
        auto mcMotherTrack = motherTrack.mcParticle_as<aod::McParticles>();
        auto mcDauTrack = dauTrack.mcParticle_as<aod::McParticles>();
        auto dChannel = getDecayChannelHe4S<aod::McParticles>(mcMotherTrack, dauIDList);
        if (dChannel == k2body && dauIDList[0] == mcDauTrack.globalIndex()) {
          isTrueSignal = true;
        }
      }

      registry.fill(HIST("hCandidateCounter"), 0);
      if (isTrueSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 0);
      }
      auto collision = kinkCand.collision_as<MCLabeledCollisionsFull>();
      if (!isGoodCollisions[collision.globalIndex()]) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 1);
      if (isTrueSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 1);
      }

      if (std::abs(dauTrack.tpcNSigmaAl()) > cutNSigmaAl) {
        continue;
      }
      float invMass = RecoDecay::m(std::array{std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()}, std::array{kinkCand.pxDaugNeut(), kinkCand.pyDaugNeut(), kinkCand.pzDaugNeut()}}, std::array{o2::constants::physics::MassAlpha, o2::constants::physics::MassPi0});
      registry.fill(HIST("hCandidateCounter"), 2);
      if (isTrueSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 2);
      }
      registry.fill(HIST("h2MassHyperhelium4sigmaPt"), kinkCand.mothSign() * kinkCand.ptMoth(), invMass);
      registry.fill(HIST("h2NSigmaAlPt"), kinkCand.mothSign() * kinkCand.ptDaug(), dauTrack.tpcNSigmaAl());

      auto bc = collision.bc_as<aod::BCs>();
      initCCDB(bc);
      Hyphe4sCandidate hyphe4sCand;
      fillCandidate(hyphe4sCand, collision, kinkCand, motherTrack, dauTrack);

      // qa for true signal
      if (isTrueSignal) {
        auto mcMotherTrack = motherTrack.mcParticle_as<aod::McParticles>();
        auto mcDauTrack = dauTrack.mcParticle_as<aod::McParticles>();
        auto mcNeutTrack = particlesMC.rawIteratorAt(dauIDList[2]);
        float posDecVtx[3] = {kinkCand.xDecVtx() + collision.posX(), kinkCand.yDecVtx() + collision.posY(), kinkCand.zDecVtx() + collision.posZ()};
        float recSVR = std::sqrt(posDecVtx[0] * posDecVtx[0] + posDecVtx[1] * posDecVtx[1]);
        registry.fill(HIST("hDiffSVx"), posDecVtx[0] - mcDauTrack.vx());
        registry.fill(HIST("hDiffSVy"), posDecVtx[1] - mcDauTrack.vy());
        registry.fill(HIST("hDiffSVz"), posDecVtx[2] - mcDauTrack.vz());
        registry.fill(HIST("h2RecSVRVsTrueSVR"), recSVR, std::hypot(mcDauTrack.vx(), mcDauTrack.vy()));
        registry.fill(HIST("h2TrueMotherDiffPtVsRecSVR"), recSVR, mcMotherTrack.pt() - kinkCand.ptMoth());
        registry.fill(HIST("h2TrueMotherDiffPzVsRecSVR"), recSVR, mcMotherTrack.pz() - kinkCand.pzMoth());
        registry.fill(HIST("h2TrueMotherDiffTglVsRecSVR"), recSVR, mcMotherTrack.pz() / mcMotherTrack.pt() - motherTrack.tgl());
        registry.fill(HIST("h2TrueMotherDiffEtaVsRecSVR"), recSVR, mcMotherTrack.eta() - motherTrack.eta());
        registry.fill(HIST("hDiffDauPx"), kinkCand.pxDaug() - mcDauTrack.px());
        registry.fill(HIST("hDiffDauPy"), kinkCand.pyDaug() - mcDauTrack.py());
        registry.fill(HIST("hDiffDauPz"), kinkCand.pzDaug() - mcDauTrack.pz());
        registry.fill(HIST("h2TrueSignalMassPt"), kinkCand.mothSign() * kinkCand.ptMoth(), invMass);
        registry.fill(HIST("h2TrueSignalNSigmaAlPt"), kinkCand.mothSign() * kinkCand.ptDaug(), dauTrack.tpcNSigmaAl());

        hyphe4sCand.isSignal = true;
        hyphe4sCand.isSignalReco = true;
        hyphe4sCand.isCollReco = true;
        hyphe4sCand.isSurvEvSelection = true;
        fillCandidateMCInfo(hyphe4sCand, mcMotherTrack, mcDauTrack, mcNeutTrack);
        mcHe4sIndices.push_back(mcMotherTrack.globalIndex());

        std::array<float, 2> dcaInfo;
        auto mcMotherTrackPar = getTrackParFromMC(mcMotherTrack);
        o2::base::Propagator::Instance()->propagateToDCABxByBz({posDecVtx[0], posDecVtx[1], posDecVtx[2]}, mcMotherTrackPar, 2.f, matCorr, &dcaInfo);
        registry.fill(HIST("hDCAXYMothToRecSV"), dcaInfo[0]);
        registry.fill(HIST("hDCAZMothToRecSV"), dcaInfo[1]);
        std::array<float, 3> pMotherAtSV = {-999.f, -999.f, -999.f};
        mcMotherTrackPar.getPxPyPzGlo(pMotherAtSV);
        registry.fill(HIST("h2TrueMotherDiffPxVsRecSVR"), recSVR, pMotherAtSV[0] - kinkCand.pxMoth());
        registry.fill(HIST("h2TrueMotherDiffPyVsRecSVR"), recSVR, pMotherAtSV[1] - kinkCand.pyMoth());

        float spKinkXY = kinkCand.pxMoth() * kinkCand.pxDaug() + kinkCand.pyMoth() * kinkCand.pyDaug();
        float spKink = spKinkXY + kinkCand.pzMoth() * kinkCand.pzDaug();
        float pMoth = std::hypot(kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth());
        float pDaug = std::hypot(kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug());

        float mothPDir[3] = {kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx()};
        float magMothPDirXY = std::hypot(mothPDir[0], mothPDir[1]);
        float magMothPDir = std::hypot(mothPDir[0], mothPDir[1], mothPDir[2]);
        float spKinkSV = mothPDir[0] * kinkCand.pxDaug() + mothPDir[1] * kinkCand.pyDaug() + mothPDir[2] * kinkCand.pzDaug();
        float sptKinkSV = mothPDir[0] * kinkCand.pxDaug() + mothPDir[1] * kinkCand.pyDaug();

        float kinkAngle = spKink / (pMoth * pDaug);
        float kinkAngleSV = spKinkSV / (magMothPDir * pDaug);
        float kinkAngleXY = spKinkXY / (kinkCand.ptMoth() * kinkCand.ptDaug());
        float kinkAngleXYSV = sptKinkSV / (magMothPDirXY * kinkCand.ptDaug());

        registry.fill(HIST("h2TrueMotherDiffPtVsKinkAngle"), kinkAngle, (mcMotherTrack.pt() - kinkCand.ptMoth()) / kinkCand.ptMoth());
        registry.fill(HIST("h2TrueMotherDiffPtVsKinkAngleSV"), kinkAngleSV, (mcMotherTrack.pt() - kinkCand.ptMoth()) / kinkCand.ptMoth());
        registry.fill(HIST("h2TrueMotherDiffPtVsKinkAngleXY"), kinkAngleXY, (mcMotherTrack.pt() - kinkCand.ptMoth()) / kinkCand.ptMoth());
        registry.fill(HIST("h2TrueMotherDiffPtVsKinkAngleXYSV"), kinkAngleXYSV, (mcMotherTrack.pt() - kinkCand.ptMoth()) / kinkCand.ptMoth());
      }

      outputMCTable(
        hyphe4sCand.primVtx[0], hyphe4sCand.primVtx[1], hyphe4sCand.primVtx[2],
        hyphe4sCand.decVtx[0], hyphe4sCand.decVtx[1], hyphe4sCand.decVtx[2],
        hyphe4sCand.isMatter,
        hyphe4sCand.lastPosMoth[0], hyphe4sCand.lastPosMoth[1], hyphe4sCand.lastPosMoth[2],
        hyphe4sCand.momMoth[0], hyphe4sCand.momMoth[1], hyphe4sCand.momMoth[2],
        hyphe4sCand.momDaug[0], hyphe4sCand.momDaug[1], hyphe4sCand.momDaug[2],
        hyphe4sCand.dcaXYMothPv, hyphe4sCand.dcaXYDauPv, hyphe4sCand.dcaKinkTopo,
        hyphe4sCand.chi2ITSMoth, hyphe4sCand.itsClusterSizeMoth, hyphe4sCand.itsClusterSizeDau,
        hyphe4sCand.nSigmaTPCDau, hyphe4sCand.nSigmaITSDau,
        hyphe4sCand.isSignal, hyphe4sCand.isSignalReco, hyphe4sCand.isCollReco, hyphe4sCand.isSurvEvSelection,
        hyphe4sCand.trueDecVtx[0], hyphe4sCand.trueDecVtx[1], hyphe4sCand.trueDecVtx[2],
        hyphe4sCand.gMomMoth[0], hyphe4sCand.gMomMoth[1], hyphe4sCand.gMomMoth[2],
        hyphe4sCand.trueMomMoth[0], hyphe4sCand.trueMomMoth[1], hyphe4sCand.trueMomMoth[2],
        hyphe4sCand.gMomDau[0], hyphe4sCand.gMomDau[1], hyphe4sCand.gMomDau[2],
        hyphe4sCand.isMothReco, hyphe4sCand.ptMoth, hyphe4sCand.pzMoth);
    }

    // fill hyperhelium4sigma signals which are not reconstructed
    for (auto const& mcparticle : particlesMC) {
      auto dChannel = getDecayChannelHe4S<aod::McParticles>(mcparticle, dauIDList);
      if (dChannel != k2body) {
        continue;
      }
      if (std::find(mcHe4sIndices.begin(), mcHe4sIndices.end(), mcparticle.globalIndex()) != mcHe4sIndices.end()) {
        continue;
      }

      Hyphe4sCandidate hyphe4sCand;
      auto mcDauTrack = particlesMC.rawIteratorAt(dauIDList[0]);
      auto mcNeutTrack = particlesMC.rawIteratorAt(dauIDList[2]);
      fillCandidateMCInfo(hyphe4sCand, mcparticle, mcDauTrack, mcNeutTrack);

      if (mcPartIndices[mcparticle.globalIndex()] != -1) {
        auto mothTrack = tracks.rawIteratorAt(mcPartIndices[mcparticle.globalIndex()]);
        fillCandidateRecoMoth(hyphe4sCand, mothTrack);
      }

      outputMCTable(
        -1, -1, -1,
        -1, -1, -1,
        -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1,
        true, false, isReconstructedMCCollisions[mcparticle.mcCollisionId()], isSelectedMCCollisions[mcparticle.mcCollisionId()],
        hyphe4sCand.trueDecVtx[0], hyphe4sCand.trueDecVtx[1], hyphe4sCand.trueDecVtx[2],
        hyphe4sCand.gMomMoth[0], hyphe4sCand.gMomMoth[1], hyphe4sCand.gMomMoth[2],
        hyphe4sCand.trueMomMoth[0], hyphe4sCand.trueMomMoth[1], hyphe4sCand.trueMomMoth[2],
        hyphe4sCand.gMomDau[0], hyphe4sCand.gMomDau[1], hyphe4sCand.gMomDau[2],
        hyphe4sCand.isMothReco, hyphe4sCand.ptMoth, hyphe4sCand.pzMoth);
    }
  }
  PROCESS_SWITCH(Hyperhelium4sigmaRecoTask, processMC, "process MC", false);
};

//--------------------------------------------------------------
// check the performance of mcparticle
struct Hyperhelium4sigmaQa {

  HistogramRegistry genQAHist{"genQAHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoQAHist{"recoQAHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 10.f}, "Binning for #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis ctBins{"ctBins", {100, 0.f, 25.f}, "Binning for c#it{t} (cm)"};
  ConfigurableAxis rigidityBins{"rigidityBins", {200, -10.f, 10.f}, "Binning for #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis nsigmaBins{"nsigmaBins", {120, -6.f, 6.f}, "Binning for n sigma"};
  ConfigurableAxis invMassBins{"invMassBins", {100, 3.85f, 4.15f}, "Binning for invariant mass (GeV/#it{c}^{2})"};
  ConfigurableAxis radiusBins{"radiusBins", {40, 0.f, 40.f}, "Binning for radius in xy plane (cm)"};

  o2::aod::ITSResponse itsResponse;

  void init(InitContext&)
  {
    if (doprocessMC == true) {
      itsResponse.setMCDefaultParameters();

      const AxisSpec pAxis{ptBins, "#it{p} (GeV/#it{c})"};
      const AxisSpec ptAxis{ptBins, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec ctAxis{ctBins, "c#it{t} (cm)"};
      const AxisSpec rigidityAxis{rigidityBins, "p/z (GeV/#it{c})"};
      const AxisSpec nsigmaAxis{nsigmaBins, "TPC n#sigma"};
      const AxisSpec itsnsigmaAxis{nsigmaBins, "ITS n#sigma"};
      const AxisSpec invMassAxis{invMassBins, "Inv Mass (GeV/#it{c}^{2})"};
      const AxisSpec diffPtAxis{200, -10.f, 10.f, "#Delta p_{T} (GeV/#it{c})"};
      const AxisSpec diffPzAxis{200, -10.f, 10.f, "#Delta p_{z} (GeV/#it{c})"};
      const AxisSpec itsRadiusAxis{radiusBins, "ITS R (cm)"};
      const AxisSpec svRadiuAxis{radiusBins, "Decay Vertex R (cm)"};

      auto hCollCounter = genQAHist.add<TH1>("hCollCounter", "hCollCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      hCollCounter->GetXaxis()->SetBinLabel(1, "Reconstructed Collisions");
      hCollCounter->GetXaxis()->SetBinLabel(2, "Selected");
      auto hMcCollCounter = genQAHist.add<TH1>("hMcCollCounter", "hMcCollCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      hMcCollCounter->GetXaxis()->SetBinLabel(1, "MC Collisions");
      hMcCollCounter->GetXaxis()->SetBinLabel(2, "Reconstructed");

      auto hGenHyperHelium4SigmaCounter = genQAHist.add<TH1>("hGenHyperHelium4SigmaCounter", "", HistType::kTH1F, {{11, 0.f, 11.f}});
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(1, "He4S All");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(2, "Matter");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(3, "AntiMatter");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(4, "#alpha + #pi^{0}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(5, "#bar{#alpha} + #pi^{0}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(6, "t + p + #pi^{0}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(7, "#bar{t} + #bar{p} + #pi^{0}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(8, "t + n + #pi^{+}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(9, "#bar{t} + #bar{n} + #pi^{+}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(10, "Tracks found");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(11, "Unexpected");

      auto hEvtSelectedHyperHelium4SigmaCounter = genQAHist.add<TH1>("hEvtSelectedHyperHelium4SigmaCounter", "", HistType::kTH1F, {{2, 0.f, 2.f}});
      hEvtSelectedHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(1, "Generated");
      hEvtSelectedHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(2, "Survived");

      genQAHist.add<TH1>("hGenHyperHelium4SigmaP", "", HistType::kTH1F, {pAxis});
      genQAHist.add<TH1>("hGenHyperHelium4SigmaPt", "", HistType::kTH1F, {ptAxis});
      genQAHist.add<TH1>("hGenHyperHelium4SigmaCt", "", HistType::kTH1F, {ctAxis});
      genQAHist.add<TH1>("hMcRecoInvMass", "", HistType::kTH1F, {invMassAxis});

      // efficiency/criteria studies for tracks which are true candidates
      hMotherCounter = recoQAHist.add<TH1>("hMotherCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hMother2BCounter = recoQAHist.add<TH1>("hMother2BCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      for (const auto& hist : {hMotherCounter, hMother2BCounter}) {
        hist->GetXaxis()->SetBinLabel(1, "Generated");
        hist->GetXaxis()->SetBinLabel(2, "Reconstructed");
        hist->GetXaxis()->SetBinLabel(3, "eta");
        hist->GetXaxis()->SetBinLabel(4, "has collision");
        hist->GetXaxis()->SetBinLabel(5, "ITSonly");
        hist->GetXaxis()->SetBinLabel(6, "ITS hits");
        hist->GetXaxis()->SetBinLabel(7, "ITS IR");
        hist->GetXaxis()->SetBinLabel(8, "ITS chi2");
        hist->GetXaxis()->SetBinLabel(9, "pt");
      }
      recoQAHist.add<TH2>("h2TrueMotherDiffPtVsTrueSVR", ";Decay Vertex R (cm);#Delta p_{T} (GeV/#it{c});", HistType::kTH2F, {svRadiuAxis, diffPtAxis});
      recoQAHist.add<TH2>("h2TrueMotherDiffPzVsTrueSVR", ";Decay Vertex R (cm);#Delta p_{z} (GeV/#it{c});", HistType::kTH2F, {svRadiuAxis, diffPzAxis});
      recoQAHist.add<TH2>("h2TrueMotherDiffTglVsTrueSVR", ";Decay Vertex R (cm);#Delta tgl;", HistType::kTH2F, {svRadiuAxis, {200, -1.f, 1.f}});
      recoQAHist.add<TH2>("h2TrueMotherDiffEtaVsTrueSVR", ";Decay Vertex R (cm);#Delta #eta;", HistType::kTH2F, {svRadiuAxis, {200, -1.f, 1.f}});
      recoQAHist.add<TH2>("h2GoodMotherDiffPtVsTrueSVR", ";Decay Vertex R (cm);#Delta p_{T} (GeV/#it{c});", HistType::kTH2F, {svRadiuAxis, diffPtAxis});
      recoQAHist.add<TH2>("h2GoodMotherDiffPzVsTrueSVR", ";Decay Vertex R (cm);#Delta p_{z} (GeV/#it{c});", HistType::kTH2F, {svRadiuAxis, diffPzAxis});

      hDauCounter[kDauAlpha] = recoQAHist.add<TH1>("hDauAlphaCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hDauCounter[kDauTriton] = recoQAHist.add<TH1>("hDauTritonCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hDauCounter[kDauProton] = recoQAHist.add<TH1>("hDauProtonCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hDauCounter[kDauChargedPion] = recoQAHist.add<TH1>("hDauPionCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hDauTPCNSigma[kDauAlpha] = recoQAHist.add<TH2>("hDauAlphaTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
      hDauTPCNSigma[kDauTriton] = recoQAHist.add<TH2>("hDauTritonTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
      hDauTPCNSigma[kDauProton] = recoQAHist.add<TH2>("hDauProtonTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
      hDauTPCNSigma[kDauChargedPion] = recoQAHist.add<TH2>("hDauPionTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
      recoQAHist.add<TH2>("hDauAlphaITSNSigmaCheck", "", HistType::kTH2F, {rigidityAxis, itsnsigmaAxis});

      hRecoMotherCounter = recoQAHist.add<TH1>("hRecoMotherCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hRecoDauAlphaCounter = recoQAHist.add<TH1>("hRecoDauAlphaCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      for (const auto& hist : {hDauCounter[kDauAlpha], hDauCounter[kDauTriton], hDauCounter[kDauProton], hDauCounter[kDauChargedPion], hRecoMotherCounter, hRecoDauAlphaCounter}) {
        hist->GetXaxis()->SetBinLabel(1, "Generated");
        hist->GetXaxis()->SetBinLabel(2, "Reconstructed");
        hist->GetXaxis()->SetBinLabel(3, "eta");
        hist->GetXaxis()->SetBinLabel(4, "has ITS&TPC");
        hist->GetXaxis()->SetBinLabel(5, "TPC crossed rows");
        hist->GetXaxis()->SetBinLabel(6, "TPC Ncls");
        hist->GetXaxis()->SetBinLabel(7, "TPC n#sigma");
        hist->GetXaxis()->SetBinLabel(8, "ITS hits");
        hist->GetXaxis()->SetBinLabel(9, "has TOF)");
      }

      recoQAHist.add<TH1>("hMotherIsPVContributer", "", HistType::kTH1F, {{2, 0.f, 2.f}});
      recoQAHist.add<TH1>("hMotherITSCls", "", HistType::kTH1F, {{8, 0.f, 8.f}});
      recoQAHist.add<TH1>("hDauAlphaIsPVContributer", "", HistType::kTH1F, {{2, 0.f, 2.f}});
      recoQAHist.add<TH1>("hDauAlphaITSCls", "", HistType::kTH1F, {{8, 0.f, 8.f}});
      recoQAHist.add<TH2>("hDauAlphaITSNSigma", "", HistType::kTH2F, {rigidityAxis, itsnsigmaAxis});
      recoQAHist.add<TH2>("hReco2BDauAlphaPVsITSNSigma", "", HistType::kTH2F, {rigidityAxis, itsnsigmaAxis});
      recoQAHist.add<TH1>("hReco2BCandidateCount", "", HistType::kTH1F, {{4, 0.f, 4.f}});
    }
  }

  Configurable<bool> skipRejectedEvents{"skipRejectedEvents", false, "Flag to skip events that fail event selection cuts"};
  Configurable<bool> doEventCut{"doEventCut", true, "Apply event selection"};
  Configurable<float> maxZVertex{"maxZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> only2BodyDecay{"only2BodyDecay", true, "Only consider 2-body decays for hyperhelium4sigma"};

  Configurable<float> etaMax{"etaMax", 1., "eta cut for tracks"};
  Configurable<float> minPtMoth{"minPtMoth", 0.5, "Minimum pT/z of the hyperhelium4sigma candidate"};
  Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
  Configurable<float> nTPCClusMinDaug{"nTPCClusMinDaug", 80, "daug NTPC clusters cut"};
  Configurable<float> itsMaxChi2{"itsMaxChi2", 36, "max chi2 for ITS"};
  Configurable<float> minRatioTPCNCls{"minRatioTPCNCls", 0.8, "min ratio of TPC crossed rows to findable clusters"};

  Preslice<aod::McParticles> permcCollision = o2::aod::mcparticle::mcCollisionId;

  // qa for mother track selection
  template <typename TTrack>
  bool motherTrackCheck(const TTrack& track, const std::shared_ptr<TH1> hist)
  {
    hist->Fill(1);

    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    hist->Fill(2);

    if (!track.has_collision()) {
      return false;
    }
    hist->Fill(3);

    if (!track.hasITS() || track.hasTPC() || track.hasTOF()) {
      return false;
    }
    hist->Fill(4);

    if (track.itsNCls() >= kITSLayers - 1) {
      return false;
    }
    hist->Fill(5);

    if (track.itsNClsInnerBarrel() != kITSInnerBarrelLayers) {
      return false;
    }
    hist->Fill(6);

    if (track.itsChi2NCl() >= itsMaxChi2) {
      return false;
    }
    hist->Fill(7);

    if (track.pt() <= minPtMoth) {
      return false;
    }
    hist->Fill(8);

    return true;
  }

  // qa for daughter track selection
  template <typename TTrack>
  bool daughterTrackCheck(const TTrack& track, const std::shared_ptr<TH1> hist, float tpcNSigma)
  {
    hist->Fill(1);

    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    hist->Fill(2);

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }
    hist->Fill(3);

    if (track.tpcNClsCrossedRows() <= minRatioTPCNCls * track.tpcNClsFindable()) {
      return false;
    }
    hist->Fill(4);

    if (track.tpcNClsFound() <= nTPCClusMinDaug) {
      return false;
    }
    hist->Fill(5);

    if (std::abs(tpcNSigma) > tpcPidNsigmaCut) {
      return false;
    }
    hist->Fill(6);

    if (track.itsNClsInnerBarrel() != 0 || track.itsNCls() > kITSInnerBarrelLayers) {
      return false;
    }
    hist->Fill(7);

    if (track.hasTOF()) {
      return false;
    }
    hist->Fill(8);

    return true;
  }

  void processData(o2::aod::Collisions const&)
  {
    // dummy process function;
  }
  PROCESS_SWITCH(Hyperhelium4sigmaQa, processData, "process data", true);

  void processMC(aod::McCollisions const& mcCollisions, aod::McParticles const& particlesMC, MCLabeledCollisionsFull const& collisions, MCLabeledTracksIU const& tracks, aod::BCs const&)
  {
    std::vector<int64_t> mcPartIndices;
    setTrackIDForMC(mcPartIndices, particlesMC, tracks);
    std::vector<bool> isSelectedMCCollisions(mcCollisions.size(), false);
    std::vector<int> dauIDList(3, -1);
    for (const auto& collision : collisions) {
      genQAHist.fill(HIST("hCollCounter"), 0.5);
      if (doEventCut && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || std::abs(collision.posZ()) > maxZVertex)) {
        continue;
      }
      genQAHist.fill(HIST("hCollCounter"), 1.5);
      isSelectedMCCollisions[collision.mcCollisionId()] = true;
    }

    for (const auto& mcCollision : mcCollisions) {
      genQAHist.fill(HIST("hMcCollCounter"), 0.5);
      if (isSelectedMCCollisions[mcCollision.globalIndex()]) { // Check that the event is reconstructed and that the reconstructed events pass the selection
        genQAHist.fill(HIST("hMcCollCounter"), 1.5);
      } else {
        if (skipRejectedEvents) {
          continue;
        }
      }

      const auto& dparticlesMC = particlesMC.sliceBy(permcCollision, mcCollision.globalIndex());

      for (const auto& mcparticle : dparticlesMC) {

        bool isMatter;
        if (mcparticle.pdgCode() == o2::constants::physics::Pdg::kHyperHelium4Sigma) {
          isMatter = true;
        } else if (mcparticle.pdgCode() == -o2::constants::physics::Pdg::kHyperHelium4Sigma) {
          isMatter = false;
        } else {
          continue;
        }

        auto dChannel = getDecayChannelHe4S<aod::McParticles>(mcparticle, dauIDList);
        if (dChannel == kNDecayChannel) {
          genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 10.5);
          continue;
        }

        // qa for mother tracks
        if (dChannel == k2body) {
          recoQAHist.fill(HIST("hMother2BCounter"), 0);
        } else {
          if (only2BodyDecay) {
            continue; // skip 3-body decays
          }
        }
        recoQAHist.fill(HIST("hMotherCounter"), 0);

        genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 0.5);
        genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), isMatter ? 1.5 : 2.5);
        genQAHist.fill(HIST("hEvtSelectedHyperHelium4SigmaCounter"), 0.5);
        if (isSelectedMCCollisions[mcCollision.globalIndex()]) {
          genQAHist.fill(HIST("hEvtSelectedHyperHelium4SigmaCounter"), 1.5);
        }

        float svPos[3] = {-999, -999, -999};
        std::vector<std::vector<float>> dauMom(kNDaughterType, std::vector<float>(3, -999.0f));
        for (const auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
          for (int type = 0; type < kNDaughterType; type++) {
            if (std::abs(mcparticleDaughter.pdgCode()) == kDaughterPDG[type]) {
              dauMom[type][0] = mcparticleDaughter.px();
              dauMom[type][1] = mcparticleDaughter.py();
              dauMom[type][2] = mcparticleDaughter.pz();

              if (type <= kDauTriton) {
                svPos[0] = mcparticleDaughter.vx();
                svPos[1] = mcparticleDaughter.vy();
                svPos[2] = mcparticleDaughter.vz();
              }

              if (type < kNChargedDaughterType) {
                hDauCounter[type]->Fill(0.f);
                // if daughter track is reconstructed
                if (type <= kNChargedDaughterType && mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
                  auto track = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
                  float tpcNSigma = getTPCNSigma(track, type);
                  daughterTrackCheck(track, hDauCounter[type], tpcNSigma);
                  if (track.hasTPC()) {
                    hDauTPCNSigma[type]->Fill(track.p() * track.sign(), tpcNSigma);
                  }
                  if (type == kDauAlpha && track.itsNCls() > kITSLayers - 2) {
                    recoQAHist.fill(HIST("hDauAlphaITSNSigmaCheck"), track.p() * track.sign(), itsResponse.nSigmaITS<o2::track::PID::Alpha>(track));
                  }
                }
              }
            }
          }
        }

        genQAHist.fill(HIST("hGenHyperHelium4SigmaP"), mcparticle.p());
        genQAHist.fill(HIST("hGenHyperHelium4SigmaPt"), mcparticle.pt());
        float ct = RecoDecay::sqrtSumOfSquares(svPos[0] - mcparticle.vx(), svPos[1] - mcparticle.vy(), svPos[2] - mcparticle.vz()) * o2::constants::physics::MassHyperHelium4Sigma / mcparticle.p();
        genQAHist.fill(HIST("hGenHyperHelium4SigmaCt"), ct);

        // if mother track is reconstructed
        if (mcPartIndices[mcparticle.globalIndex()] != -1) {
          auto motherTrack = tracks.rawIteratorAt(mcPartIndices[mcparticle.globalIndex()]);
          bool isGoodMother = motherTrackCheck(motherTrack, hMotherCounter);
          if (dChannel == k2body) {
            motherTrackCheck(motherTrack, hMother2BCounter);
          }
          float svR = RecoDecay::sqrtSumOfSquares(svPos[0], svPos[1]);
          float diffpt = mcparticle.pt() - 2 * motherTrack.pt();
          float diffpz = mcparticle.pz() - 2 * motherTrack.pz();

          recoQAHist.fill(HIST("h2TrueMotherDiffPtVsTrueSVR"), svR, diffpt);
          recoQAHist.fill(HIST("h2TrueMotherDiffPzVsTrueSVR"), svR, diffpz);
          recoQAHist.fill(HIST("h2TrueMotherDiffTglVsTrueSVR"), svR, mcparticle.pz() / mcparticle.pt() - motherTrack.tgl());
          recoQAHist.fill(HIST("h2TrueMotherDiffEtaVsTrueSVR"), svR, mcparticle.eta() - motherTrack.eta());
          if (isGoodMother) {
            recoQAHist.fill(HIST("h2GoodMotherDiffPtVsTrueSVR"), svR, diffpt);
            recoQAHist.fill(HIST("h2GoodMotherDiffPzVsTrueSVR"), svR, diffpz);
          }
          // if mother track and charged daughters are all reconstructed
          bool isDauReconstructed = mcPartIndices[dauIDList[0]] != -1 && (dChannel == k2body ? true : mcPartIndices[dauIDList[1]] != -1);
          if (isDauReconstructed) {
            genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 9.5);

            // qa for bc matching for reconstructed tracks
            if (dChannel == k2body) {
              auto daughterTrack = tracks.rawIteratorAt(mcPartIndices[dauIDList[0]]);
              bool isMoth = motherTrackCheck(motherTrack, hRecoMotherCounter);
              bool isDaug = daughterTrackCheck(daughterTrack, hRecoDauAlphaCounter, daughterTrack.tpcNSigmaAl());

              recoQAHist.fill(HIST("hReco2BCandidateCount"), 0.5);
              recoQAHist.fill(HIST("hRecoMotherCounter"), 0.5);
              recoQAHist.fill(HIST("hMotherITSCls"), motherTrack.itsNCls());
              recoQAHist.fill(HIST("hRecoDauAlphaCounter"), 0.5);
              recoQAHist.fill(HIST("hMotherIsPVContributer"), motherTrack.isPVContributor() ? 1.5 : 0.5);
              recoQAHist.fill(HIST("hDauAlphaIsPVContributer"), daughterTrack.isPVContributor() ? 1.5 : 0.5);

              float itsNSigma = itsResponse.nSigmaITS<o2::track::PID::Alpha>(daughterTrack);
              if (daughterTrack.hasITS()) {
                recoQAHist.fill(HIST("hDauAlphaITSNSigma"), daughterTrack.sign() * daughterTrack.p(), itsNSigma);
                recoQAHist.fill(HIST("hDauAlphaITSCls"), daughterTrack.itsNCls());
              }

              if (motherTrack.has_collision() && daughterTrack.has_collision()) {
                recoQAHist.fill(HIST("hReco2BCandidateCount"), 1.5);
                if (motherTrack.collisionId() == daughterTrack.collisionId()) {
                  recoQAHist.fill(HIST("hReco2BCandidateCount"), 2.5);
                }
              }

              if (isMoth && isDaug) {
                recoQAHist.fill(HIST("hReco2BCandidateCount"), 3.5);
                recoQAHist.fill(HIST("hReco2BDauAlphaPVsITSNSigma"), daughterTrack.sign() * daughterTrack.p(), itsNSigma);
              }
            }
          }
        }

        // qa for branching ratios and invariant mass
        if (dChannel == k2body) {
          genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), isMatter ? 3.5 : 4.5);
          float hyperHelium4SigmaMCMass = RecoDecay::m(std::array{std::array{dauMom[kDauAlpha][0], dauMom[kDauAlpha][1], dauMom[kDauAlpha][2]}, std::array{dauMom[kDauPion0][0], dauMom[kDauPion0][1], dauMom[kDauPion0][2]}}, std::array{o2::constants::physics::MassAlpha, o2::constants::physics::MassPi0});
          genQAHist.fill(HIST("hMcRecoInvMass"), hyperHelium4SigmaMCMass);
        } else if (dChannel == k3body_p) {
          genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), isMatter ? 5.5 : 6.5);
          float hyperHelium4SigmaMCMass = RecoDecay::m(std::array{std::array{dauMom[kDauTriton][0], dauMom[kDauTriton][1], dauMom[kDauTriton][2]}, std::array{dauMom[kDauProton][0], dauMom[kDauProton][1], dauMom[kDauProton][2]}, std::array{dauMom[kDauPion0][0], dauMom[kDauPion0][1], dauMom[kDauPion0][2]}}, std::array{o2::constants::physics::MassTriton, o2::constants::physics::MassProton, o2::constants::physics::MassPi0});
          genQAHist.fill(HIST("hMcRecoInvMass"), hyperHelium4SigmaMCMass);
        } else if (dChannel == k3body_n) {
          genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), isMatter ? 7.5 : 8.5);
          float hyperHelium4SigmaMCMass = RecoDecay::m(std::array{std::array{dauMom[kDauTriton][0], dauMom[kDauTriton][1], dauMom[kDauTriton][2]}, std::array{dauMom[kDauNeutron][0], dauMom[kDauNeutron][1], dauMom[kDauNeutron][2]}, std::array{dauMom[kDauChargedPion][0], dauMom[kDauChargedPion][1], dauMom[kDauChargedPion][2]}}, std::array{o2::constants::physics::MassTriton, o2::constants::physics::MassNeutron, o2::constants::physics::MassPiPlus});
          genQAHist.fill(HIST("hMcRecoInvMass"), hyperHelium4SigmaMCMass);
        }
      }
    }
  }
  PROCESS_SWITCH(Hyperhelium4sigmaQa, processMC, "do QA for MC prodcutions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Hyperhelium4sigmaRecoTask>(cfgc),
    adaptAnalysisTask<Hyperhelium4sigmaQa>(cfgc),
  };
}
