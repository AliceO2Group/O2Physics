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
/// \file hyperkinkRecoTask.cxx
/// \brief QA and analysis task for kink decay of hypernuclei
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#include "PWGLF/DataModel/LFHyperNucleiKinkTables.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFPIDTOFGenericTables.h"
#include "PWGLF/Utils/pidTOFGeneric.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

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

o2::common::core::MetadataHelper metadataInfo;

using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::EvTimeTOFFT0>;
using MCLabeledCollisionsFull = soa::Join<CollisionsFull, aod::McCollisionLabels>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullAl, aod::pidTPCFullTr, aod::EvTimeTOFFT0ForTrack>;
using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels, aod::pidTOFFullTr, aod::pidTOFFullAl>;

enum PartType {
  kHypertriton = 0,
  kHyperhelium4sigma
};

enum DaughterType {
  kDaugCharged = 0,
  kDaugNeutral,
  kNDaughterType
};

namespace
{
constexpr std::array<float, 7> LayerRadii{2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
constexpr int kITSLayers = 7;
constexpr int kITSInnerBarrelLayers = 3;
// constexpr int kITSOuterBarrelLayers = 4;

std::shared_ptr<TH1> hMothCounter;
std::shared_ptr<TH1> hDaugCounter;
std::shared_ptr<TH2> hDaugTPCNSigma;
std::shared_ptr<TH1> hRecoMothCounter;
std::shared_ptr<TH1> hRecoDaugCounter;
} // namespace

//--------------------------------------------------------------
struct H3LDecay {
  enum Channel {
    k2bodyNeutral = 0, // triton, pion0
    k2bodyCharged,
    k3bodyCharged,
    kNChannel
  };

  template <class TMCTrackTo, typename TMCParticle>
  static Channel getDecayChannel(TMCParticle const& particle, std::vector<int>& list)
  {
    if (std::abs(particle.pdgCode()) != o2::constants::physics::Pdg::kHyperTriton) {
      return kNChannel;
    }

    list.clear();
    list.resize(2, -1);

    bool haveHelium3 = false, haveAntiHelium3 = false, haveDeuteron = false, haveAntiDeuteron = false;
    bool haveProton = false, haveAntiProton = false, havePionPlus = false, havePionMinus = false;
    bool haveTriton = false, haveAntiTriton = false, havePion0 = false;
    for (const auto& mcDaughter : particle.template daughters_as<TMCTrackTo>()) {
      if (mcDaughter.pdgCode() == o2::constants::physics::Pdg::kTriton) {
        haveTriton = true;
        list[0] = mcDaughter.globalIndex();
      }
      if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kTriton) {
        haveAntiTriton = true;
        list[0] = mcDaughter.globalIndex();
      }
      if (mcDaughter.pdgCode() == PDG_t::kPi0) {
        havePion0 = true;
        list[1] = mcDaughter.globalIndex();
      }
      if (mcDaughter.pdgCode() == o2::constants::physics::Pdg::kHelium3) {
        haveHelium3 = true;
      }
      if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kHelium3) {
        haveAntiHelium3 = true;
      }
      if (mcDaughter.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
        haveDeuteron = true;
      }
      if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
        haveAntiDeuteron = true;
      }
      if (mcDaughter.pdgCode() == PDG_t::kProton) {
        haveProton = true;
      }
      if (mcDaughter.pdgCode() == -PDG_t::kProton) {
        haveAntiProton = true;
      }
      if (mcDaughter.pdgCode() == PDG_t::kPiPlus) {
        havePionPlus = true;
      }
      if (mcDaughter.pdgCode() == -PDG_t::kPiPlus) {
        havePionMinus = true;
      }
    }

    if ((haveTriton && havePion0) || (haveAntiTriton && havePion0)) {
      return H3LDecay::k2bodyNeutral;
    } else if ((haveHelium3 && havePionMinus) || (haveAntiHelium3 && havePionPlus)) {
      return H3LDecay::k2bodyCharged;
    } else if ((haveDeuteron && haveProton && havePionMinus) || (haveAntiDeuteron && haveAntiProton && havePionPlus)) {
      return H3LDecay::k3bodyCharged;
    } else {
      return kNChannel;
    }
  }
};

//--------------------------------------------------------------
// Check the decay channel of hyperhelium4sigma
struct He4SDecay {
  enum Channel {
    k2body = 0, // helium4, pion0
    k3body_p,   // triton, proton, pion0
    k3body_n,   // triton, neutron, pion+
    kNChannel
  };

  template <class TMCTrackTo, typename TMCParticle>
  static Channel getDecayChannel(TMCParticle const& particle, std::vector<int>& list)
  {
    if (std::abs(particle.pdgCode()) != o2::constants::physics::Pdg::kHyperHelium4Sigma) {
      return kNChannel;
    }

    list.clear();
    list.resize(2, -1);

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
      }
      if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kTriton) {
        haveAntiTriton = true;
      }
      if (mcDaughter.pdgCode() == PDG_t::kProton) {
        haveProton = true;
      }
      if (mcDaughter.pdgCode() == -PDG_t::kProton) {
        haveAntiProton = true;
      }
      if (mcDaughter.pdgCode() == PDG_t::kNeutron) {
        haveNeuteron = true;
      }
      if (mcDaughter.pdgCode() == -PDG_t::kNeutron) {
        haveAntiNeuteron = true;
      }
      if (mcDaughter.pdgCode() == PDG_t::kPiPlus) {
        havePionPlus = true;
      }
      if (mcDaughter.pdgCode() == -PDG_t::kPiPlus) {
        havePionMinus = true;
      }
      if (mcDaughter.pdgCode() == PDG_t::kPi0) {
        havePion0 = true;
        list[1] = mcDaughter.globalIndex();
      }
    }

    if ((haveAlpha && havePion0) || (haveAntiAlpha && havePion0)) {
      return He4SDecay::k2body;
    } else if ((haveTriton && haveProton && havePion0) || (haveAntiTriton && haveAntiProton && havePion0)) {
      return He4SDecay::k3body_p;
    } else if ((haveTriton && haveNeuteron && havePionPlus) || (haveAntiTriton && haveAntiNeuteron && havePionMinus)) {
      return He4SDecay::k3body_n;
    }

    return kNChannel;
  }
};

//--------------------------------------------------------------
// Extract track parameters from a mcparticle, use global coordinates as the local one
template <typename TrackPrecision = float, typename T>
o2::track::TrackParametrization<TrackPrecision> getTrackParFromMC(const T& mcparticle, int charge = 1)
{
  int sign = mcparticle.pdgCode() > 0 ? 1 : -1; // ok for hypernuclei
  TrackPrecision snp = mcparticle.py() / (mcparticle.pt() + 1.e-10f);
  TrackPrecision tgl = mcparticle.pz() / (mcparticle.pt() + 1.e-10f);
  std::array<TrackPrecision, o2::track::kNParams> arraypar = {mcparticle.vy(), mcparticle.vz(), snp,
                                                              tgl, charge * sign / (mcparticle.pt() + 1.e-10f)};
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
// get ITSNSigma for daughter track
template <typename TTrack>
float getITSNSigma(const TTrack& track, o2::aod::ITSResponse& itsResponse, o2::track::PID partType)
{
  float nSigma = -999.f;
  switch (partType) {
    case o2::track::PID::Alpha:
      nSigma = itsResponse.nSigmaITS<o2::track::PID::Alpha>(track);
      break;
    case o2::track::PID::Triton:
      nSigma = itsResponse.nSigmaITS<o2::track::PID::Triton>(track);
      break;
    default:
      break;
  }
  return nSigma;
}

//--------------------------------------------------------------
// get default TOFNSigma for daughter track
template <typename TTrack>
float getDefaultTOFNSigma(const TTrack& track, o2::track::PID partType)
{
  float nSigma = -999.f;
  switch (partType) {
    case o2::track::PID::Alpha:
      nSigma = track.tofNSigmaAl();
      break;
    case o2::track::PID::Triton:
      nSigma = track.tofNSigmaTr();
      break;
    default:
      break;
  }
  return nSigma;
}

//--------------------------------------------------------------
// get TPCNSigma for daughter track
template <typename TTrack>
float getTPCNSigma(const TTrack& track, o2::track::PID partType)
{
  float nSigma = -999.f;
  switch (partType) {
    case o2::track::PID::Alpha:
      nSigma = track.tpcNSigmaAl();
      break;
    case o2::track::PID::Triton:
      nSigma = track.tpcNSigmaTr();
      break;
    default:
      break;
  }
  return nSigma;
}

//--------------------------------------------------------------
struct HypKinkCandidate {

  bool isMatter = false;

  std::array<float, 3> posPV = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> posSV = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> lastPosMoth = {0.0f, 0.0f, 0.0f}; // last position of mother track at the radii of ITS layer which has the outermost update
  std::array<float, 3> momMothSV = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> momDaugSV = {0.0f, 0.0f, 0.0f};

  float dcaXYMothPv = -999.f;
  float dcaXYDaugPv = -999.f;
  float dcaKinkTopo = -999.f;

  float chi2ITSMoth = 0.0f;
  uint32_t itsClusterSizeMoth = 0u;
  uint32_t itsClusterSizeDaug = 0u;
  float tpcMomDaug = -999.f;
  float tpcSignalDaug = -999.f;
  int16_t tpcNClsPIDDaug = 0u;
  float nSigmaTPCDaug = -999.f;
  float nSigmaITSDaug = -999.f;
  float nSigmaTOFDaug = -999.f; // recalculated TOF NSigma

  // mc information
  bool isSignal = false;
  bool isSignalReco = false;
  bool isCollReco = false;
  bool isSurvEvSelection = false;

  std::array<float, 3> truePosSV = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> trueMomMothPV = {0.0f, 0.0f, 0.0f}; // generated mother momentum at primary vertex
  std::array<float, 3> trueMomMothSV = {0.0f, 0.0f, 0.0f}; // true mother momentum at decay vertex
  std::array<float, 3> trueMomDaugSV = {0.0f, 0.0f, 0.0f}; // true daughter momentum at decay vertex

  bool isMothReco = false;
  std::array<float, 3> momMothPV = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> updateMomMothPV = {0.0f, 0.0f, 0.0f}; // mother momentum at primary vertex after update using PV
};

//--------------------------------------------------------------
// analysis task for 2-body kink decay
struct HyperkinkRecoTask {

  Produces<aod::HypKinkCand> outputDataTable;
  Produces<aod::MCHypKinkCand> outputMCTable;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry registry{"registry", {}};

  Configurable<int> hypoMoth{"hypoMoth", kHypertriton, "Mother particle hypothesis"};
  // Configurable for event selection
  Configurable<bool> doEventCut{"doEventCut", true, "Apply event selection"};
  Configurable<float> maxZVertex{"maxZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutTPCNSigmaDaug{"cutTPCNSigmaDaug", 5, "TPC NSigma cut for daughter tracks"};
  Configurable<float> cutTOFNSigmaDaug{"cutTOFNSigmaDaug", 1000, "TOF NSigma cut for daughter tracks"};
  Configurable<bool> askTOFForDaug{"askTOFForDaug", false, "If true, ask for TOF signal"};
  Configurable<float> minDaugPt{"minDaugPt", 1.0f, "Minimum pT of daughter track (GeV/c)"};
  Configurable<float> minDCADaugToPV{"minDCADaugToPV", 0.f, "Minimum DCA of daughter track to primary vertex (cm)"};
  Configurable<float> maxDCADaugToPV{"maxDCADaugToPV", 999.0f, "Maximum DCA of daughter track to primary vertex (cm)"};
  Configurable<float> maxQtAP{"maxQtAP", 1.f, "Maximum qT of Armenteros-Podolanski Plot"};
  Configurable<float> maxDCAKinkTopo{"maxDCAKinkTopo", 1.f, "Maximum DCA of kink topology (cm)"};

  // CCDB options
  Configurable<double> inputBz{"inputBz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdbPath{"ccdbPath", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};

  int mRunNumber;
  float mBz;
  o2::base::MatLayerCylSet* lut = nullptr;

  o2::aod::ITSResponse itsResponse;

  float massMoth = 999.f;
  float massChargedDaug = 999.f;
  float massNeutralDaug = 999.f;
  int pdgMoth = 0;
  std::array<int, kNDaughterType> pdgDaug = {o2::constants::physics::Pdg::kTriton, PDG_t::kPi0}; // pdgcode of charged (0) and neutral (1) daughter particles
  o2::track::PID pidTypeDaug = o2::track::PID::Triton;

  // secondary TOF PID
  o2::aod::pidtofgeneric::TofPidNewCollision<FullTracksExtIU::iterator> secondaryTOFPID;          // to be updated in Init based on the hypothesis
  o2::aod::pidtofgeneric::TofPidNewCollision<MCLabeledTracksIU::iterator> secondaryTOFPIDLabeled; // to be updated in Init based on the hypothesis
  // TOF response and input parameters
  o2::pid::tof::TOFResoParamsV3 mRespParamsV3;
  o2::aod::pidtofgeneric::TOFCalibConfig mTOFCalibConfig; // TOF Calib configuration

  void init(InitContext& initContext)
  {
    if (hypoMoth == kHypertriton) {
      massMoth = o2::constants::physics::MassHyperTriton;
      massChargedDaug = o2::constants::physics::MassTriton;
      massNeutralDaug = o2::constants::physics::MassPi0;
      pdgMoth = o2::constants::physics::Pdg::kHyperTriton;
      pdgDaug[kDaugCharged] = o2::constants::physics::Pdg::kTriton;
      pdgDaug[kDaugNeutral] = PDG_t::kPi0;
      pidTypeDaug = o2::track::PID::Triton;
    } else if (hypoMoth == kHyperhelium4sigma) {
      massMoth = o2::constants::physics::MassHyperHelium4Sigma;
      massChargedDaug = o2::constants::physics::MassAlpha;
      massNeutralDaug = o2::constants::physics::MassPi0;
      pdgMoth = o2::constants::physics::Pdg::kHyperHelium4Sigma;
      pdgDaug[kDaugCharged] = o2::constants::physics::Pdg::kAlpha;
      pdgDaug[kDaugNeutral] = PDG_t::kPi0;
      pidTypeDaug = o2::track::PID::Alpha;
    } else {
      LOG(fatal) << "Unknown mother particle hypothesis";
    }

    // Axes
    const AxisSpec vertexZAxis{100, -15., 15., "vtx_{Z} [cm]"};
    const AxisSpec ptAxis{50, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nSigmaAxis{120, -6.f, 6.f, "n#sigma"};
    AxisSpec massAxis(100, 2.94, 3.2, "m (GeV/#it{c}^{2})");
    if (hypoMoth == kHyperhelium4sigma) {
      massAxis = AxisSpec{100, 3.85, 4.25, "m (GeV/#it{c}^{2})"};
    }
    const AxisSpec deltaPtAxis{200, -10.f, 10.f, "#Delta #it{p}_{T} (GeV/#it{c})"};
    const AxisSpec deltaPzAxis{200, -10.f, 10.f, "#Delta #it{p}_{z} (GeV/#it{c})"};
    const AxisSpec recRadiusAxis{40, 0.f, 40.f, "Rec SV R (cm)"};

    registry.add<TH1>("hEventCounter", "hEventCounter", HistType::kTH1F, {{2, 0, 2}});
    registry.add<TH1>("hVertexZCollision", "hVertexZCollision", HistType::kTH1F, {vertexZAxis});
    registry.add<TH1>("hCandidateCounter", "hCandidateCounter", HistType::kTH1F, {{8, 0, 8}});

    if (doprocessMC == true) {
      itsResponse.setMCDefaultParameters();

      registry.add<TH1>("hTrueCandidateCounter", "hTrueCandidateCounter", HistType::kTH1F, {{8, 0, 8}});
      registry.add<TH1>("hDeltaSVx", ";#Delta x (cm);", HistType::kTH1F, {{200, -2, 2}});
      registry.add<TH1>("hDeltaSVy", ";#Delta y (cm);", HistType::kTH1F, {{200, -2, 2}});
      registry.add<TH1>("hDeltaSVz", ";#Delta z (cm);", HistType::kTH1F, {{200, -2, 2}});
      registry.add<TH2>("hDeltaSVRVsTrueSVR", ";True SVR (cm);#Delta R (cm)", HistType::kTH2F, {{200, 0, 40}, {200, -2, 2}});
      registry.add<TH2>("h2RecSVRVsTrueSVR", ";Rec SV R (cm);True SV R (cm);", HistType::kTH2F, {recRadiusAxis, {40, 0, 40}});
      registry.add<TH2>("h2TrueMotherDeltaPtVsRecSVR", ";Rec SV R (cm);#Delta #it{p}_{T} (GeV/#it{c});", HistType::kTH2F, {recRadiusAxis, deltaPtAxis});
      registry.add<TH2>("h2TrueMotherDeltaEtaVsRecSVR", ";Rec SV R (cm);#Delta #eta;", HistType::kTH2F, {recRadiusAxis, {200, -0.1, 0.1}});
      registry.add<TH1>("hDeltaDauPx", ";#Delta p_{x} (GeV/#it{c}); ", HistType::kTH1F, {{200, -2, 2}});
      registry.add<TH1>("hDeltaDauPy", ";#Delta p_{y} (GeV/#it{c}); ", HistType::kTH1F, {{200, -2, 2}});
      registry.add<TH1>("hDeltaDauPz", ";#Delta p_{z} (GeV/#it{c}); ", HistType::kTH1F, {{200, -2, 2}});
      registry.add<TH2>("h2TrueSignalMassPt", "h2TrueSignalMassPt", HistType::kTH2F, {{ptAxis, massAxis}});
      registry.add<TH2>("h2TrueDaugTPCNSigmaPt", "h2TrueDaugTPCNSigmaPt", HistType::kTH2F, {{ptAxis, nSigmaAxis}});

      registry.add<TH1>("hDCAXYMothToRecSV", "hDCAXYMothToRecSV", HistType::kTH1F, {{200, -10, 10}});
      registry.add<TH1>("hDCAZMothToRecSV", "hDCAZMothToRecSV", HistType::kTH1F, {{200, -10, 10}});

      registry.add<TH1>("hDaugOldTOFNSigma_CorrectCol", "hDaugOldTOFNSigma_CorrectCol", HistType::kTH1F, {{600, -300, 300}});
      registry.add<TH1>("hDaugOldTOFNSigma_WrongCol", "hDaugOldTOFNSigma_WrongCol", HistType::kTH1F, {{600, -300, 300}});
      registry.add<TH1>("hDaugNewTOFNSigma_CorrectCol", "hDaugNewTOFNSigma_CorrectCol", HistType::kTH1F, {{600, -300, 300}});
      registry.add<TH1>("hDaugNewTOFNSigma_WrongCol", "hDaugNewTOFNSigma_WrongCol", HistType::kTH1F, {{600, -300, 300}});

      registry.add<TH1>("hTrueSignalInvMassNegNeutDaugE", "hTrueSignalInvMassNegNeutDaugE", HistType::kTH1F, {{1000, 2.8, 4, "m (GeV/#it{c}^{2})"}});
    }

    registry.add<TH2>("h2MothMassPt", "h2MothMassPt", HistType::kTH2F, {{ptAxis, massAxis}});
    registry.add<TH2>("h2DaugTPCNSigmaPt", "h2DaugTPCNSigmaPt", HistType::kTH2F, {{ptAxis, nSigmaAxis}});

    registry.add<TH1>("hInvMassNegNeutDaugE", "hInvMassNegNeutDaugE", HistType::kTH1F, {{1000, 2.8, 4, "m (GeV/#it{c}^{2})"}});

    ccdb->setURL(ccdbPath);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    mTOFCalibConfig.metadataInfo = metadataInfo;
    mTOFCalibConfig.inheritFromBaseTask(initContext);
    mTOFCalibConfig.initSetup(mRespParamsV3, ccdb);
    secondaryTOFPID.SetPidType(pidTypeDaug);
    secondaryTOFPIDLabeled.SetPidType(pidTypeDaug);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    ccdb->clearCache(grpmagPath.value.data());
    LOG(info) << "Initializing CCDB for run " << mRunNumber;
    o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(grpmagPath, mRunNumber);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    mBz = grpmag->getNominalL3Field();

    if (!lut) {
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }
    o2::base::Propagator::Instance()->setMatLUT(lut);
    LOG(info) << "Task initialized for run " << mRunNumber << " with magnetic field " << mBz << " kZG";

    mTOFCalibConfig.processSetup(mRespParamsV3, ccdb, bc);
  }

  // ______________________________________________________________
  // get recalulate TOFNSigma for daughter track
  template <bool isMC, typename TCollision, typename TTrack>
  double getTOFNSigma(o2::pid::tof::TOFResoParamsV3 const& parameters, TTrack const& track, TCollision const& collision, TCollision const& originalcol)
  {
    // TOF PID of deuteron
    if (track.has_collision() && track.hasTOF()) {
      if constexpr (isMC) {
        return secondaryTOFPIDLabeled.GetTOFNSigma(parameters, track, originalcol, collision);
      } else {
        return secondaryTOFPID.GetTOFNSigma(parameters, track, originalcol, collision);
      }
    }
    return -999;
  }

  template <typename TCollision, typename TKindCandidate, typename TTrack>
  void fillCandidate(HypKinkCandidate& hypkinkCand, TCollision const& collision, TKindCandidate const& kinkCand, TTrack const& trackMoth, TTrack const& trackDaug)
  {
    hypkinkCand.isMatter = kinkCand.mothSign() > 0;
    hypkinkCand.posPV[0] = collision.posX();
    hypkinkCand.posPV[1] = collision.posY();
    hypkinkCand.posPV[2] = collision.posZ();
    hypkinkCand.posSV[0] = kinkCand.xDecVtx() + collision.posX();
    hypkinkCand.posSV[1] = kinkCand.yDecVtx() + collision.posY();
    hypkinkCand.posSV[2] = kinkCand.zDecVtx() + collision.posZ();

    hypkinkCand.momMothSV[0] = kinkCand.pxMoth();
    hypkinkCand.momMothSV[1] = kinkCand.pyMoth();
    hypkinkCand.momMothSV[2] = kinkCand.pzMoth();
    hypkinkCand.momDaugSV[0] = kinkCand.pxDaug();
    hypkinkCand.momDaugSV[1] = kinkCand.pyDaug();
    hypkinkCand.momDaugSV[2] = kinkCand.pzDaug();

    hypkinkCand.dcaXYMothPv = kinkCand.dcaMothPv();
    hypkinkCand.dcaXYDaugPv = kinkCand.dcaDaugPv();
    hypkinkCand.dcaKinkTopo = kinkCand.dcaKinkTopo();

    fillCandidateRecoMoth(hypkinkCand, collision, trackMoth);

    hypkinkCand.itsClusterSizeDaug = trackDaug.itsClusterSizes();
    hypkinkCand.tpcMomDaug = trackDaug.tpcInnerParam() * trackDaug.sign();
    hypkinkCand.tpcSignalDaug = trackDaug.tpcSignal();
    hypkinkCand.tpcNClsPIDDaug = trackDaug.tpcNClsPID();
    hypkinkCand.nSigmaTPCDaug = getTPCNSigma(trackDaug, pidTypeDaug);
    hypkinkCand.nSigmaITSDaug = getITSNSigma(trackDaug, itsResponse, pidTypeDaug);

    int lastLayerMoth = 0;
    for (int i = 6; i >= 0; i--) {
      if (trackMoth.itsClusterMap() & (1 << i)) {
        lastLayerMoth = i;
        break;
      }
    }
    auto trackparMother = getTrackParCov(trackMoth);
    o2::base::Propagator::Instance()->PropagateToXBxByBz(trackparMother, LayerRadii[lastLayerMoth]);
    std::array<float, 3> posLastHit{-999.f};
    trackparMother.getXYZGlo(posLastHit);
    hypkinkCand.lastPosMoth[0] = posLastHit[0];
    hypkinkCand.lastPosMoth[1] = posLastHit[1];
    hypkinkCand.lastPosMoth[2] = posLastHit[2];
  }

  template <typename TCollision, typename TTrack>
  void fillCandidateRecoMoth(HypKinkCandidate& hypkinkCand, TCollision const& collision, TTrack const& trackMoth)
  {
    hypkinkCand.isMothReco = true;
    hypkinkCand.chi2ITSMoth = trackMoth.itsChi2NCl();
    hypkinkCand.itsClusterSizeMoth = trackMoth.itsClusterSizes();

    auto motherTrackPar = getTrackParCov(trackMoth);
    o2::dataformats::VertexBase primaryVtx = {{collision.posX(), collision.posY(), collision.posZ()}, {collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()}};
    std::array<float, 3> pMotherPv = {-999.f};
    std::array<float, 3> updatePMotherPv = {-999.f};
    if (o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVtx, motherTrackPar, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrLUT)) {
      motherTrackPar.getPxPyPzGlo(pMotherPv);
      if (motherTrackPar.update(primaryVtx)) {
        motherTrackPar.getPxPyPzGlo(updatePMotherPv);
      }
    }

    hypkinkCand.momMothPV[0] = pMotherPv[0];
    hypkinkCand.momMothPV[1] = pMotherPv[1];
    hypkinkCand.momMothPV[2] = pMotherPv[2];
    hypkinkCand.updateMomMothPV[0] = updatePMotherPv[0];
    hypkinkCand.updateMomMothPV[1] = updatePMotherPv[1];
    hypkinkCand.updateMomMothPV[2] = updatePMotherPv[2];
  }

  template <typename TMCParticle>
  void fillCandidateMCInfo(HypKinkCandidate& hypkinkCand, TMCParticle const& mcMothTrack, TMCParticle const& mcDaugTrack, TMCParticle const& mcNeutDaugTrack)
  {
    hypkinkCand.truePosSV[0] = mcDaugTrack.vx();
    hypkinkCand.truePosSV[1] = mcDaugTrack.vy();
    hypkinkCand.truePosSV[2] = mcDaugTrack.vz();
    hypkinkCand.trueMomMothPV[0] = mcMothTrack.px();
    hypkinkCand.trueMomMothPV[1] = mcMothTrack.py();
    hypkinkCand.trueMomMothPV[2] = mcMothTrack.pz();
    hypkinkCand.trueMomMothSV[0] = mcDaugTrack.px() + mcNeutDaugTrack.px();
    hypkinkCand.trueMomMothSV[1] = mcDaugTrack.py() + mcNeutDaugTrack.py();
    hypkinkCand.trueMomMothSV[2] = mcDaugTrack.pz() + mcNeutDaugTrack.pz();
    hypkinkCand.trueMomDaugSV[0] = mcDaugTrack.px();
    hypkinkCand.trueMomDaugSV[1] = mcDaugTrack.py();
    hypkinkCand.trueMomDaugSV[2] = mcDaugTrack.pz();
  }

  void processData(CollisionsFull const& collisions, aod::KinkCands const& KinkCands, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
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

      auto daugTrack = kinkCand.trackDaug_as<FullTracksExtIU>();
      float tpcNSigmaDaug = getTPCNSigma(daugTrack, pidTypeDaug);
      if (std::abs(tpcNSigmaDaug) > cutTPCNSigmaDaug) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 2);

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto motherTrack = kinkCand.trackMoth_as<FullTracksExtIU>();
      float nSigmaTOF = -999.f;
      if (daugTrack.hasTOF() && daugTrack.has_collision()) {
        auto originalDaugCol = daugTrack.collision_as<CollisionsFull>();
        nSigmaTOF = getTOFNSigma<false>(mRespParamsV3, daugTrack, collision, originalDaugCol);
      }
      if ((daugTrack.hasTOF() || askTOFForDaug) && std::abs(nSigmaTOF) > cutTOFNSigmaDaug) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 3);

      if (kinkCand.ptDaug() < minDaugPt) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 4);

      if (std::abs(kinkCand.dcaDaugPv()) < minDCADaugToPV || std::abs(kinkCand.dcaDaugPv()) > maxDCADaugToPV) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 5);

      if (kinkCand.dcaKinkTopo() > maxDCAKinkTopo) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 6);

      float p2Moth = kinkCand.pxMoth() * kinkCand.pxMoth() + kinkCand.pyMoth() * kinkCand.pyMoth() + kinkCand.pzMoth() * kinkCand.pzMoth();
      float p2Daug = kinkCand.pxDaug() * kinkCand.pxDaug() + kinkCand.pyDaug() * kinkCand.pyDaug() + kinkCand.pzDaug() * kinkCand.pzDaug();
      float sqKink = kinkCand.pxMoth() * kinkCand.pxDaug() + kinkCand.pyMoth() * kinkCand.pyDaug() + kinkCand.pzMoth() * kinkCand.pzDaug();
      float qt = std::sqrt(p2Daug - sqKink * sqKink / p2Moth);
      if (qt > maxQtAP) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 7);

      float invMass = RecoDecay::m(std::array{std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()}, std::array{kinkCand.pxDaugNeut(), kinkCand.pyDaugNeut(), kinkCand.pzDaugNeut()}}, std::array{massChargedDaug, massNeutralDaug});
      registry.fill(HIST("h2MothMassPt"), kinkCand.mothSign() * kinkCand.ptMoth(), invMass);
      registry.fill(HIST("h2DaugTPCNSigmaPt"), kinkCand.mothSign() * kinkCand.ptDaug(), tpcNSigmaDaug);

      HypKinkCandidate hypkinkCand;
      fillCandidate(hypkinkCand, collision, kinkCand, motherTrack, daugTrack);
      hypkinkCand.nSigmaTOFDaug = nSigmaTOF;

      outputDataTable(
        mBz > 0 ? 1 : -1,
        hypkinkCand.posPV[0], hypkinkCand.posPV[1], hypkinkCand.posPV[2],
        hypkinkCand.posSV[0], hypkinkCand.posSV[1], hypkinkCand.posSV[2],
        hypkinkCand.isMatter,
        hypkinkCand.lastPosMoth[0], hypkinkCand.lastPosMoth[1], hypkinkCand.lastPosMoth[2],
        hypkinkCand.momMothSV[0], hypkinkCand.momMothSV[1], hypkinkCand.momMothSV[2],
        hypkinkCand.momDaugSV[0], hypkinkCand.momDaugSV[1], hypkinkCand.momDaugSV[2],
        hypkinkCand.dcaXYMothPv, hypkinkCand.dcaXYDaugPv, hypkinkCand.dcaKinkTopo,
        hypkinkCand.chi2ITSMoth, hypkinkCand.itsClusterSizeMoth, hypkinkCand.itsClusterSizeDaug,
        hypkinkCand.tpcMomDaug, hypkinkCand.tpcSignalDaug, hypkinkCand.tpcNClsPIDDaug,
        hypkinkCand.nSigmaTPCDaug, hypkinkCand.nSigmaITSDaug, hypkinkCand.nSigmaTOFDaug,
        hypkinkCand.momMothPV[0], hypkinkCand.momMothPV[1], hypkinkCand.momMothPV[2],
        hypkinkCand.updateMomMothPV[0], hypkinkCand.updateMomMothPV[1], hypkinkCand.updateMomMothPV[2]);
    }
  }
  PROCESS_SWITCH(HyperkinkRecoTask, processData, "process data", true);

  void processMC(MCLabeledCollisionsFull const& collisions, aod::KinkCands const& KinkCands, MCLabeledTracksIU const& tracks, aod::McParticles const& particlesMC, aod::McCollisions const& mcCollisions, aod::BCsWithTimestamps const&)
  {
    std::vector<int64_t> mcPartIndices;
    setTrackIDForMC(mcPartIndices, particlesMC, tracks);
    std::vector<int64_t> signalIndicesPool;
    std::vector<bool> isReconstructedMCCollisions(mcCollisions.size(), false);
    std::vector<bool> isSelectedMCCollisions(mcCollisions.size(), false);
    std::vector<bool> isGoodCollisions(collisions.size(), false);
    std::vector<int> dauIDList(2, -1);

    for (const auto& collision : collisions) {
      if (collision.has_mcCollision()) {
        isReconstructedMCCollisions[collision.mcCollisionId()] = true;
      }
      registry.fill(HIST("hEventCounter"), 0);
      if (doEventCut && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || std::abs(collision.posZ()) > maxZVertex)) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1);
      registry.fill(HIST("hVertexZCollision"), collision.posZ());
      if (collision.has_mcCollision()) {
        isSelectedMCCollisions[collision.mcCollisionId()] = true;
      }
      isGoodCollisions[collision.globalIndex()] = true;
    }

    for (const auto& kinkCand : KinkCands) {
      auto motherTrack = kinkCand.trackMoth_as<MCLabeledTracksIU>();
      auto daugTrack = kinkCand.trackDaug_as<MCLabeledTracksIU>();

      bool isKinkSignal = false;
      if (motherTrack.has_mcParticle() && daugTrack.has_mcParticle()) {
        auto mcMothTrack = motherTrack.mcParticle_as<aod::McParticles>();
        auto mcDaugTrack = daugTrack.mcParticle_as<aod::McParticles>();
        if (hypoMoth == kHypertriton) {
          auto dChannel = H3LDecay::getDecayChannel<aod::McParticles>(mcMothTrack, dauIDList);
          if (dChannel == H3LDecay::k2bodyNeutral && dauIDList[0] == mcDaugTrack.globalIndex()) {
            isKinkSignal = true;
          }
        } else if (hypoMoth == kHyperhelium4sigma) {
          auto dChannel = He4SDecay::getDecayChannel<aod::McParticles>(mcMothTrack, dauIDList);
          if (dChannel == He4SDecay::k2body && dauIDList[0] == mcDaugTrack.globalIndex()) {
            isKinkSignal = true;
          }
        }
      }

      registry.fill(HIST("hCandidateCounter"), 0);
      if (isKinkSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 0);
      }
      auto collision = kinkCand.collision_as<MCLabeledCollisionsFull>();
      if (!isGoodCollisions[collision.globalIndex()]) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 1);
      if (isKinkSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 1);
      }

      float tpcNSigmaDaug = getTPCNSigma(daugTrack, pidTypeDaug);
      if (std::abs(tpcNSigmaDaug) > cutTPCNSigmaDaug) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 2);
      if (isKinkSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 2);
      }

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      float nSigmaTOF = -999.f;
      if (daugTrack.hasTOF() && daugTrack.has_collision()) {
        auto originalDaugCol = daugTrack.collision_as<MCLabeledCollisionsFull>();
        float defaultNSigmaTOF = getDefaultTOFNSigma(daugTrack, pidTypeDaug);
        nSigmaTOF = getTOFNSigma<true>(mRespParamsV3, daugTrack, collision, originalDaugCol);
        if (originalDaugCol.globalIndex() == collision.globalIndex()) {
          registry.fill(HIST("hDaugOldTOFNSigma_CorrectCol"), defaultNSigmaTOF);
          registry.fill(HIST("hDaugNewTOFNSigma_CorrectCol"), nSigmaTOF);
        } else {
          registry.fill(HIST("hDaugOldTOFNSigma_WrongCol"), defaultNSigmaTOF);
          registry.fill(HIST("hDaugNewTOFNSigma_WrongCol"), nSigmaTOF);
        }
      }
      if ((daugTrack.hasTOF() || askTOFForDaug) && std::abs(nSigmaTOF) > cutTOFNSigmaDaug) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 3);
      if (isKinkSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 3);
      }

      if (kinkCand.ptDaug() < minDaugPt) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 4);
      if (isKinkSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 4);
      }

      if (std::abs(kinkCand.dcaDaugPv()) < minDCADaugToPV || std::abs(kinkCand.dcaDaugPv()) > maxDCADaugToPV) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 5);
      if (isKinkSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 5);
      }

      if (kinkCand.dcaKinkTopo() > maxDCAKinkTopo) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 6);
      if (isKinkSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 6);
      }

      float p2Moth = kinkCand.pxMoth() * kinkCand.pxMoth() + kinkCand.pyMoth() * kinkCand.pyMoth() + kinkCand.pzMoth() * kinkCand.pzMoth();
      float p2Daug = kinkCand.pxDaug() * kinkCand.pxDaug() + kinkCand.pyDaug() * kinkCand.pyDaug() + kinkCand.pzDaug() * kinkCand.pzDaug();
      float sqKink = kinkCand.pxMoth() * kinkCand.pxDaug() + kinkCand.pyMoth() * kinkCand.pyDaug() + kinkCand.pzMoth() * kinkCand.pzDaug();
      float qt = std::sqrt(p2Daug - sqKink * sqKink / p2Moth);
      if (qt > maxQtAP) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 7);
      if (isKinkSignal) {
        registry.fill(HIST("hTrueCandidateCounter"), 7);
      }

      float invMass = RecoDecay::m(std::array{std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()}, std::array{kinkCand.pxDaugNeut(), kinkCand.pyDaugNeut(), kinkCand.pzDaugNeut()}}, std::array{massChargedDaug, massNeutralDaug});
      registry.fill(HIST("h2MothMassPt"), kinkCand.mothSign() * kinkCand.ptMoth(), invMass);
      registry.fill(HIST("h2DaugTPCNSigmaPt"), kinkCand.mothSign() * kinkCand.ptDaug(), tpcNSigmaDaug);

      // qa for energy of charged daughter greater than the mother track with mass hypothesis
      float pCharDaug = std::hypot(kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug());
      float pMoth = std::hypot(kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth());
      float chargedDauE = std::hypot(pCharDaug, massChargedDaug);
      float motherE = std::hypot(pMoth, massMoth);
      if (chargedDauE < motherE) {
        registry.fill(HIST("hInvMassNegNeutDaugE"), invMass);
        if (isKinkSignal) {
          registry.fill(HIST("hTrueSignalInvMassNegNeutDaugE"), invMass);
        }
      }

      HypKinkCandidate hypkinkCand;
      fillCandidate(hypkinkCand, collision, kinkCand, motherTrack, daugTrack);
      hypkinkCand.nSigmaTOFDaug = nSigmaTOF;

      std::array<float, 3> posDecVtx = {kinkCand.xDecVtx() + collision.posX(), kinkCand.yDecVtx() + collision.posY(), kinkCand.zDecVtx() + collision.posZ()};
      float recSVR = std::hypot(posDecVtx[0], posDecVtx[1]);

      // QA, store mcInfo for true signals
      if (isKinkSignal) {
        auto mcMothTrack = motherTrack.mcParticle_as<aod::McParticles>();
        auto mcDaugTrack = daugTrack.mcParticle_as<aod::McParticles>();
        auto mcNeutTrack = particlesMC.rawIteratorAt(dauIDList[1]);
        float trueSVR = std::hypot(mcDaugTrack.vx(), mcDaugTrack.vy());
        registry.fill(HIST("hDeltaSVx"), posDecVtx[0] - mcDaugTrack.vx());
        registry.fill(HIST("hDeltaSVy"), posDecVtx[1] - mcDaugTrack.vy());
        registry.fill(HIST("hDeltaSVz"), posDecVtx[2] - mcDaugTrack.vz());
        registry.fill(HIST("hDeltaSVRVsTrueSVR"), trueSVR, recSVR - trueSVR);
        registry.fill(HIST("h2RecSVRVsTrueSVR"), recSVR, trueSVR);
        registry.fill(HIST("h2TrueMotherDeltaPtVsRecSVR"), recSVR, mcMothTrack.pt() - kinkCand.ptMoth());
        registry.fill(HIST("h2TrueMotherDeltaEtaVsRecSVR"), recSVR, mcMothTrack.eta() - motherTrack.eta());
        registry.fill(HIST("hDeltaDauPx"), kinkCand.pxDaug() - mcDaugTrack.px());
        registry.fill(HIST("hDeltaDauPy"), kinkCand.pyDaug() - mcDaugTrack.py());
        registry.fill(HIST("hDeltaDauPz"), kinkCand.pzDaug() - mcDaugTrack.pz());
        registry.fill(HIST("h2TrueSignalMassPt"), kinkCand.mothSign() * kinkCand.ptMoth(), invMass);
        registry.fill(HIST("h2TrueDaugTPCNSigmaPt"), kinkCand.mothSign() * kinkCand.ptDaug(), tpcNSigmaDaug);

        hypkinkCand.isSignal = true;
        hypkinkCand.isSignalReco = true;
        fillCandidateMCInfo(hypkinkCand, mcMothTrack, mcDaugTrack, mcNeutTrack);
        signalIndicesPool.push_back(mcMothTrack.globalIndex());

        std::array<float, 2> dcaInfo;
        auto mcMothTrackPar = getTrackParFromMC(mcMothTrack, 2);
        o2::base::Propagator::Instance()->propagateToDCABxByBz({posDecVtx[0], posDecVtx[1], posDecVtx[2]}, mcMothTrackPar, 2.f, matCorr, &dcaInfo);
        registry.fill(HIST("hDCAXYMothToRecSV"), dcaInfo[0]);
        registry.fill(HIST("hDCAZMothToRecSV"), dcaInfo[1]);
      }

      hypkinkCand.isCollReco = true;
      hypkinkCand.isSurvEvSelection = true;

      outputMCTable(
        mBz > 0 ? 1 : -1,
        hypkinkCand.posPV[0], hypkinkCand.posPV[1], hypkinkCand.posPV[2],
        hypkinkCand.posSV[0], hypkinkCand.posSV[1], hypkinkCand.posSV[2],
        hypkinkCand.isMatter,
        hypkinkCand.lastPosMoth[0], hypkinkCand.lastPosMoth[1], hypkinkCand.lastPosMoth[2],
        hypkinkCand.momMothSV[0], hypkinkCand.momMothSV[1], hypkinkCand.momMothSV[2],
        hypkinkCand.momDaugSV[0], hypkinkCand.momDaugSV[1], hypkinkCand.momDaugSV[2],
        hypkinkCand.dcaXYMothPv, hypkinkCand.dcaXYDaugPv, hypkinkCand.dcaKinkTopo,
        hypkinkCand.chi2ITSMoth, hypkinkCand.itsClusterSizeMoth, hypkinkCand.itsClusterSizeDaug,
        hypkinkCand.tpcMomDaug, hypkinkCand.tpcSignalDaug, hypkinkCand.tpcNClsPIDDaug,
        hypkinkCand.nSigmaTPCDaug, hypkinkCand.nSigmaITSDaug, hypkinkCand.nSigmaTOFDaug,
        hypkinkCand.isSignal, hypkinkCand.isSignalReco, hypkinkCand.isCollReco, hypkinkCand.isSurvEvSelection,
        hypkinkCand.truePosSV[0], hypkinkCand.truePosSV[1], hypkinkCand.truePosSV[2],
        hypkinkCand.trueMomMothPV[0], hypkinkCand.trueMomMothPV[1], hypkinkCand.trueMomMothPV[2],
        hypkinkCand.trueMomMothSV[0], hypkinkCand.trueMomMothSV[1], hypkinkCand.trueMomMothSV[2],
        hypkinkCand.trueMomDaugSV[0], hypkinkCand.trueMomDaugSV[1], hypkinkCand.trueMomDaugSV[2],
        hypkinkCand.isMothReco, hypkinkCand.momMothPV[0], hypkinkCand.momMothPV[1], hypkinkCand.momMothPV[2],
        hypkinkCand.updateMomMothPV[0], hypkinkCand.updateMomMothPV[1], hypkinkCand.updateMomMothPV[2]);
    }

    // fill kink signals which are not reconstructed
    for (auto const& mcparticle : particlesMC) {
      bool isKinkSignal = false;
      if (hypoMoth == kHypertriton) {
        auto dChannel = H3LDecay::getDecayChannel<aod::McParticles>(mcparticle, dauIDList);
        if (dChannel == H3LDecay::k2bodyNeutral) {
          isKinkSignal = true;
        }
      } else if (hypoMoth == kHyperhelium4sigma) {
        auto dChannel = He4SDecay::getDecayChannel<aod::McParticles>(mcparticle, dauIDList);
        if (dChannel == He4SDecay::k2body) {
          isKinkSignal = true;
        }
      }
      if (!isKinkSignal) {
        continue;
      }

      if (std::find(signalIndicesPool.begin(), signalIndicesPool.end(), mcparticle.globalIndex()) != signalIndicesPool.end()) {
        continue;
      }

      HypKinkCandidate hypkinkCand;
      auto mcDaugTrack = particlesMC.rawIteratorAt(dauIDList[0]);
      auto mcNeutTrack = particlesMC.rawIteratorAt(dauIDList[1]);
      fillCandidateMCInfo(hypkinkCand, mcparticle, mcDaugTrack, mcNeutTrack);

      if (mcPartIndices[mcparticle.globalIndex()] != -1) {
        auto mothTrack = tracks.rawIteratorAt(mcPartIndices[mcparticle.globalIndex()]);
        if (mothTrack.has_collision()) {
          auto collision = mothTrack.collision_as<MCLabeledCollisionsFull>();
          auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
          initCCDB(bc);
          fillCandidateRecoMoth(hypkinkCand, collision, mothTrack);
        }
      }

      outputMCTable(
        mBz > 0 ? 1 : -1,
        -1, -1, -1,
        -1, -1, -1,
        -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        true, false, isReconstructedMCCollisions[mcparticle.mcCollisionId()], isSelectedMCCollisions[mcparticle.mcCollisionId()],
        hypkinkCand.truePosSV[0], hypkinkCand.truePosSV[1], hypkinkCand.truePosSV[2],
        hypkinkCand.trueMomMothPV[0], hypkinkCand.trueMomMothPV[1], hypkinkCand.trueMomMothPV[2],
        hypkinkCand.trueMomMothSV[0], hypkinkCand.trueMomMothSV[1], hypkinkCand.trueMomMothSV[2],
        hypkinkCand.trueMomDaugSV[0], hypkinkCand.trueMomDaugSV[1], hypkinkCand.trueMomDaugSV[2],
        hypkinkCand.isMothReco, hypkinkCand.momMothPV[0], hypkinkCand.momMothPV[1], hypkinkCand.momMothPV[2],
        hypkinkCand.updateMomMothPV[0], hypkinkCand.updateMomMothPV[1], hypkinkCand.updateMomMothPV[2]);
    }
  }
  PROCESS_SWITCH(HyperkinkRecoTask, processMC, "process MC", false);
};

//--------------------------------------------------------------
// check the performance of mcparticle
struct HyperkinkQa {

  Configurable<int> hypoMoth{"hypoMoth", kHypertriton, "Mother particle hypothesis"};

  HistogramRegistry genQAHist{"genQAHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoQAHist{"recoQAHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 10.f}, "Binning for #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis ctBins{"ctBins", {100, 0.f, 25.f}, "Binning for c#it{t} (cm)"};
  ConfigurableAxis rigidityBins{"rigidityBins", {200, -10.f, 10.f}, "Binning for #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis nsigmaBins{"nsigmaBins", {120, -6.f, 6.f}, "Binning for n sigma"};
  ConfigurableAxis invMassBins{"invMassBins", {100, 3.85f, 4.15f}, "Binning for invariant mass (GeV/#it{c}^{2})"};
  ConfigurableAxis radiusBins{"radiusBins", {40, 0.f, 40.f}, "Binning for radius in xy plane (cm)"};

  o2::aod::ITSResponse itsResponse;

  int charge = 1;
  float massMoth = 999.f;
  float massChargedDaug = 999.f;
  float massNeutralDaug = 999.f;
  int pdgMoth = 0;
  std::array<int, kNDaughterType> pdgDaug = {o2::constants::physics::Pdg::kTriton, PDG_t::kPi0}; // pdgcode of charged (0) and neutral (1) daughter particles
  o2::track::PID pidTypeDaug = o2::track::PID::Triton;

  void init(InitContext&)
  {
    if (doprocessMC == true) {
      itsResponse.setMCDefaultParameters();

      if (hypoMoth == kHypertriton) {
        charge = 1;
        massMoth = o2::constants::physics::MassHyperTriton;
        massChargedDaug = o2::constants::physics::MassTriton;
        massNeutralDaug = o2::constants::physics::MassPi0;
        pdgMoth = o2::constants::physics::Pdg::kHyperTriton;
        pdgDaug[kDaugCharged] = o2::constants::physics::Pdg::kTriton;
        pdgDaug[kDaugNeutral] = PDG_t::kPi0;
        pidTypeDaug = o2::track::PID::Triton;
      } else if (hypoMoth == kHyperhelium4sigma) {
        charge = 2;
        massMoth = o2::constants::physics::MassHyperHelium4Sigma;
        massChargedDaug = o2::constants::physics::MassAlpha;
        massNeutralDaug = o2::constants::physics::MassPi0;
        pdgMoth = o2::constants::physics::Pdg::kHyperHelium4Sigma;
        pdgDaug[kDaugCharged] = o2::constants::physics::Pdg::kAlpha;
        pdgDaug[kDaugNeutral] = PDG_t::kPi0;
        pidTypeDaug = o2::track::PID::Alpha;
      } else {
        LOG(fatal) << "Unknown mother particle hypothesis";
      }

      const AxisSpec pAxis{ptBins, "#it{p} (GeV/#it{c})"};
      const AxisSpec ptAxis{ptBins, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec ctAxis{ctBins, "c#it{t} (cm)"};
      const AxisSpec rigidityAxis{rigidityBins, "p/z (GeV/#it{c})"};
      const AxisSpec nsigmaAxis{nsigmaBins, "TPC n#sigma"};
      const AxisSpec itsnsigmaAxis{nsigmaBins, "ITS n#sigma"};
      const AxisSpec invMassAxis{invMassBins, "Inv Mass (GeV/#it{c}^{2})"};
      const AxisSpec deltaPtAxis{200, -10.f, 10.f, "#Delta p_{T} (GeV/#it{c})"};
      const AxisSpec deltaPzAxis{200, -10.f, 10.f, "#Delta p_{z} (GeV/#it{c})"};
      const AxisSpec itsRadiusAxis{radiusBins, "ITS R (cm)"};
      const AxisSpec svRadiuAxis{radiusBins, "Decay Vertex R (cm)"};

      auto hCollCounter = genQAHist.add<TH1>("hCollCounter", "hCollCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      hCollCounter->GetXaxis()->SetBinLabel(1, "Reconstructed Collisions");
      hCollCounter->GetXaxis()->SetBinLabel(2, "Selected");
      auto hMcCollCounter = genQAHist.add<TH1>("hMcCollCounter", "hMcCollCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      hMcCollCounter->GetXaxis()->SetBinLabel(1, "MC Collisions");
      hMcCollCounter->GetXaxis()->SetBinLabel(2, "Reconstructed");

      if (hypoMoth == kHypertriton) {
        auto hGenHyperMothCounter = genQAHist.add<TH1>("hGenHyperMothCounter", "", HistType::kTH1F, {{10, 0.f, 10.f}});
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(1, "H3L All");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(2, "Matter");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(3, "AntiMatter");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(4, "t + #pi^{0}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(5, "#bar{t} + #pi^{0}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(6, "{}^{3}He + #pi^{-}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(7, "{}^{3}#bar{He} + #pi^{+}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(8, "d + p + #pi^{-}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(9, "#bar{d} + #bar{p} + #pi^{+}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(10, "Others");
      } else if (hypoMoth == kHyperhelium4sigma) {
        auto hGenHyperMothCounter = genQAHist.add<TH1>("hGenHyperMothCounter", "", HistType::kTH1F, {{10, 0.f, 10.f}});
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(1, "He4S All");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(2, "Matter");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(3, "AntiMatter");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(4, "#alpha + #pi^{0}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(5, "#bar{#alpha} + #pi^{0}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(6, "t + p + #pi^{0}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(7, "#bar{t} + #bar{p} + #pi^{0}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(8, "t + n + #pi^{+}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(9, "#bar{t} + #bar{n} + #pi^{+}");
        hGenHyperMothCounter->GetXaxis()->SetBinLabel(10, "Others");
      }

      auto hEvtSelectedHyperMothCounter = genQAHist.add<TH1>("hEvtSelectedHyperMothCounter", "", HistType::kTH1F, {{2, 0.f, 2.f}});
      hEvtSelectedHyperMothCounter->GetXaxis()->SetBinLabel(1, "Generated");
      hEvtSelectedHyperMothCounter->GetXaxis()->SetBinLabel(2, "Survived");

      genQAHist.add<TH1>("hGenHyperMothP", "", HistType::kTH1F, {pAxis});
      genQAHist.add<TH1>("hGenHyperMothPt", "", HistType::kTH1F, {ptAxis});
      genQAHist.add<TH1>("hGenHyperMothCt", "", HistType::kTH1F, {ctAxis});
      genQAHist.add<TH1>("hMcRecoInvMass", "", HistType::kTH1F, {invMassAxis});

      // efficiency/criteria studies for tracks which are true candidates
      hMothCounter = recoQAHist.add<TH1>("hMothCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hMothCounter->GetXaxis()->SetBinLabel(1, "Generated");
      hMothCounter->GetXaxis()->SetBinLabel(2, "Reconstructed");
      hMothCounter->GetXaxis()->SetBinLabel(3, "eta");
      hMothCounter->GetXaxis()->SetBinLabel(4, "has collision");
      hMothCounter->GetXaxis()->SetBinLabel(5, "ITSonly");
      hMothCounter->GetXaxis()->SetBinLabel(6, "ITS hits");
      hMothCounter->GetXaxis()->SetBinLabel(7, "ITS IR");
      hMothCounter->GetXaxis()->SetBinLabel(8, "ITS chi2");
      hMothCounter->GetXaxis()->SetBinLabel(9, "pt");
      recoQAHist.add<TH2>("h2TrueMotherDeltaPtVsTrueSVR", ";Decay Vertex R (cm);#Delta p_{T} (GeV/#it{c});", HistType::kTH2F, {svRadiuAxis, deltaPtAxis});
      recoQAHist.add<TH2>("h2TrueMotherDeltaEtaVsTrueSVR", ";Decay Vertex R (cm);#Delta #eta;", HistType::kTH2F, {svRadiuAxis, {200, -1.f, 1.f}});
      recoQAHist.add<TH2>("h2GoodMotherDeltaPtVsTrueSVR", ";Decay Vertex R (cm);#Delta p_{T} (GeV/#it{c});", HistType::kTH2F, {svRadiuAxis, deltaPtAxis});

      hDaugCounter = recoQAHist.add<TH1>("hDaugCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hDaugTPCNSigma = recoQAHist.add<TH2>("hDaugTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});

      hRecoMothCounter = recoQAHist.add<TH1>("hRecoMothCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hRecoDaugCounter = recoQAHist.add<TH1>("hRecoDaugCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      for (const auto& hist : {hDaugCounter, hRecoMothCounter, hRecoDaugCounter}) {
        hist->GetXaxis()->SetBinLabel(1, "Generated");
        hist->GetXaxis()->SetBinLabel(2, "Reconstructed");
        hist->GetXaxis()->SetBinLabel(3, "eta");
        hist->GetXaxis()->SetBinLabel(4, "has ITS&TPC");
        hist->GetXaxis()->SetBinLabel(5, "TPC crossed rows");
        hist->GetXaxis()->SetBinLabel(6, "TPC Ncls");
        hist->GetXaxis()->SetBinLabel(7, "TPC n#sigma");
        hist->GetXaxis()->SetBinLabel(8, "ITS hits");
        hist->GetXaxis()->SetBinLabel(9, "has TOF");
      }

      recoQAHist.add<TH1>("hMothIsPVContributer", "", HistType::kTH1F, {{2, 0.f, 2.f}});
      recoQAHist.add<TH1>("hMothITSCls", "", HistType::kTH1F, {{8, 0.f, 8.f}});
      recoQAHist.add<TH1>("hDaugIsPVContributer", "", HistType::kTH1F, {{2, 0.f, 2.f}});
      recoQAHist.add<TH1>("hDaugITSCls", "", HistType::kTH1F, {{8, 0.f, 8.f}});
      recoQAHist.add<TH2>("hDaugITSNSigma", "", HistType::kTH2F, {rigidityAxis, itsnsigmaAxis});
      recoQAHist.add<TH2>("hRecoDaugPVsITSNSigma", "", HistType::kTH2F, {rigidityAxis, itsnsigmaAxis});
      recoQAHist.add<TH1>("hRecoCandidateCount", "", HistType::kTH1F, {{4, 0.f, 4.f}});

      recoQAHist.add<TH1>("hDiffZTracks", "", HistType::kTH1F, {{200, -100.f, 100.f}});
      recoQAHist.add<TH1>("hDiffAbsZTracks", "", HistType::kTH1F, {{200, -100.f, 100.f}});
      recoQAHist.add<TH1>("hDiffXTracks", "", HistType::kTH1F, {{200, -100.f, 100.f}});
    }
  }

  Configurable<bool> skipRejectedEvents{"skipRejectedEvents", false, "Flag to skip events that fail event selection cuts"};
  Configurable<bool> doEventCut{"doEventCut", true, "Apply event selection"};
  Configurable<float> maxZVertex{"maxZVertex", 10.0f, "Accepted z-vertex range (cm)"};

  Configurable<float> etaMax{"etaMax", 1., "eta cut for tracks"};
  Configurable<float> minPtMoth{"minPtMoth", 0.5, "Minimum pT/z of the mother track"};
  Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
  Configurable<float> nTPCClusMinDaug{"nTPCClusMinDaug", 80, "daug NTPC clusters cut"};
  Configurable<float> itsMaxChi2{"itsMaxChi2", 36, "max chi2 for ITS"};
  Configurable<float> minRatioTPCNCls{"minRatioTPCNCls", 0.8, "min ratio of TPC crossed rows to findable clusters"};

  Preslice<aod::McParticles> permcCollision = o2::aod::mcparticle::mcCollisionId;

  // QA for mother track selection
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
  PROCESS_SWITCH(HyperkinkQa, processData, "process data", true);

  void processMC(aod::McCollisions const& mcCollisions, aod::McParticles const& particlesMC, MCLabeledCollisionsFull const& collisions, MCLabeledTracksIU const& tracks, aod::BCsWithTimestamps const&)
  {
    std::vector<int64_t> mcPartIndices;
    setTrackIDForMC(mcPartIndices, particlesMC, tracks);
    std::vector<bool> isSelectedMCCollisions(mcCollisions.size(), false);
    std::vector<int> dauIDList(2, -1);

    for (const auto& collision : collisions) {
      genQAHist.fill(HIST("hCollCounter"), 0.5);
      if (doEventCut && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || std::abs(collision.posZ()) > maxZVertex)) {
        continue;
      }
      genQAHist.fill(HIST("hCollCounter"), 1.5);
      if (collision.has_mcCollision()) {
        isSelectedMCCollisions[collision.mcCollisionId()] = true;
      }
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
        if (std::abs(mcparticle.pdgCode()) != pdgMoth) {
          continue;
        }
        bool isMatter = mcparticle.pdgCode() > 0;
        genQAHist.fill(HIST("hGenHyperMothCounter"), 0.5);
        genQAHist.fill(HIST("hGenHyperMothCounter"), isMatter ? 1.5 : 2.5);

        // QA for decay channels
        bool isKinkSignal = false;
        if (hypoMoth == kHypertriton) {
          auto dChannel = H3LDecay::getDecayChannel<aod::McParticles>(mcparticle, dauIDList);
          if (dChannel == H3LDecay::k2bodyNeutral) {
            genQAHist.fill(HIST("hGenHyperMothCounter"), isMatter ? 3.5 : 4.5);
            isKinkSignal = true;
          } else if (dChannel == H3LDecay::k2bodyCharged) {
            genQAHist.fill(HIST("hGenHyperMothCounter"), isMatter ? 5.5 : 6.5);
          } else if (dChannel == H3LDecay::k3bodyCharged) {
            genQAHist.fill(HIST("hGenHyperMothCounter"), isMatter ? 7.5 : 8.5);
          } else if (dChannel == H3LDecay::kNChannel) {
            genQAHist.fill(HIST("hGenHyperMothCounter"), 9.5);
            continue;
          }
        } else if (hypoMoth == kHyperhelium4sigma) {
          auto dChannel = He4SDecay::getDecayChannel<aod::McParticles>(mcparticle, dauIDList);
          if (dChannel == He4SDecay::k2body) {
            genQAHist.fill(HIST("hGenHyperMothCounter"), isMatter ? 3.5 : 4.5);
            isKinkSignal = true;
          } else if (dChannel == He4SDecay::k3body_p) {
            genQAHist.fill(HIST("hGenHyperMothCounter"), isMatter ? 5.5 : 6.5);
          } else if (dChannel == He4SDecay::k3body_n) {
            genQAHist.fill(HIST("hGenHyperMothCounter"), isMatter ? 7.5 : 8.5);
          } else if (dChannel == He4SDecay::kNChannel) {
            genQAHist.fill(HIST("hGenHyperMothCounter"), 9.5);
            continue;
          }
        }

        if (!isKinkSignal) {
          continue;
        }
        recoQAHist.fill(HIST("hMothCounter"), 0);

        genQAHist.fill(HIST("hEvtSelectedHyperMothCounter"), 0.5);
        if (isSelectedMCCollisions[mcCollision.globalIndex()]) {
          genQAHist.fill(HIST("hEvtSelectedHyperMothCounter"), 1.5);
        }

        float svPos[3] = {-999, -999, -999};
        std::vector<std::vector<float>> dauMom(kNDaughterType, std::vector<float>(3, -999.0f));
        for (const auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
          for (int type = 0; type < kNDaughterType; type++) {
            if (std::abs(mcparticleDaughter.pdgCode()) == pdgDaug[type]) {
              dauMom[type][0] = mcparticleDaughter.px();
              dauMom[type][1] = mcparticleDaughter.py();
              dauMom[type][2] = mcparticleDaughter.pz();

              if (type == kDaugCharged) {
                svPos[0] = mcparticleDaughter.vx();
                svPos[1] = mcparticleDaughter.vy();
                svPos[2] = mcparticleDaughter.vz();

                // if daughter track is reconstructed
                if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
                  hDaugCounter->Fill(0.f);
                  auto track = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
                  float tpcNSigma = getTPCNSigma(track, pidTypeDaug);
                  daughterTrackCheck(track, hDaugCounter, tpcNSigma);
                  if (track.hasTPC()) {
                    hDaugTPCNSigma->Fill(track.p() * track.sign(), tpcNSigma);
                  }
                }
              }
            }
          }
        }

        genQAHist.fill(HIST("hGenHyperMothP"), mcparticle.p());
        genQAHist.fill(HIST("hGenHyperMothPt"), mcparticle.pt());
        float ct = RecoDecay::sqrtSumOfSquares(svPos[0] - mcparticle.vx(), svPos[1] - mcparticle.vy(), svPos[2] - mcparticle.vz()) * massMoth / mcparticle.p();
        genQAHist.fill(HIST("hGenHyperMothCt"), ct);
        float hypermothMCMass = RecoDecay::m(std::array{std::array{dauMom[kDaugCharged][0], dauMom[kDaugCharged][1], dauMom[kDaugCharged][2]}, std::array{dauMom[kDaugNeutral][0], dauMom[kDaugNeutral][1], dauMom[kDaugNeutral][2]}}, std::array{massChargedDaug, massNeutralDaug});
        genQAHist.fill(HIST("hMcRecoInvMass"), hypermothMCMass);

        // if mother track is reconstructed
        if (mcPartIndices[mcparticle.globalIndex()] != -1) {
          auto motherTrack = tracks.rawIteratorAt(mcPartIndices[mcparticle.globalIndex()]);
          bool isGoodMother = motherTrackCheck(motherTrack, hMothCounter);
          float svR = RecoDecay::sqrtSumOfSquares(svPos[0], svPos[1]);
          float deltapt = mcparticle.pt() - charge * motherTrack.pt();

          recoQAHist.fill(HIST("h2TrueMotherDeltaPtVsTrueSVR"), svR, deltapt);
          recoQAHist.fill(HIST("h2TrueMotherDeltaEtaVsTrueSVR"), svR, mcparticle.eta() - motherTrack.eta());
          if (isGoodMother) {
            recoQAHist.fill(HIST("h2GoodMotherDeltaPtVsTrueSVR"), svR, deltapt);
          }

          // if mother track and charged daughters are all reconstructed
          bool isDauReconstructed = mcPartIndices[dauIDList[0]] != -1;
          if (isDauReconstructed) {
            auto daughterTrack = tracks.rawIteratorAt(mcPartIndices[dauIDList[0]]);
            bool isMoth = motherTrackCheck(motherTrack, hRecoMothCounter);
            bool isDaug = daughterTrackCheck(daughterTrack, hRecoDaugCounter, getTPCNSigma(daughterTrack, pidTypeDaug));

            recoQAHist.fill(HIST("hRecoCandidateCount"), 0.5);
            recoQAHist.fill(HIST("hRecoMothCounter"), 0.5);
            recoQAHist.fill(HIST("hMothITSCls"), motherTrack.itsNCls());
            recoQAHist.fill(HIST("hRecoDaugCounter"), 0.5);
            recoQAHist.fill(HIST("hMothIsPVContributer"), motherTrack.isPVContributor() ? 1.5 : 0.5);
            recoQAHist.fill(HIST("hDaugIsPVContributer"), daughterTrack.isPVContributor() ? 1.5 : 0.5);

            if (svR > LayerRadii[3]) {
              recoQAHist.fill(HIST("hDiffZTracks"), daughterTrack.z() - motherTrack.z());
              recoQAHist.fill(HIST("hDiffAbsZTracks"), std::abs(daughterTrack.z()) - std::abs(motherTrack.z()));
              recoQAHist.fill(HIST("hDiffXTracks"), daughterTrack.x() - motherTrack.x());
            }

            float itsNSigma = getITSNSigma(daughterTrack, itsResponse, pidTypeDaug);
            if (daughterTrack.hasITS()) {
              recoQAHist.fill(HIST("hDaugITSNSigma"), daughterTrack.sign() * daughterTrack.p(), itsNSigma);
              recoQAHist.fill(HIST("hDaugITSCls"), daughterTrack.itsNCls());
            }

            if (motherTrack.has_collision() && daughterTrack.has_collision()) {
              recoQAHist.fill(HIST("hRecoCandidateCount"), 1.5);
              if (motherTrack.collisionId() == daughterTrack.collisionId()) {
                recoQAHist.fill(HIST("hRecoCandidateCount"), 2.5);
              }
            }

            if (isMoth && isDaug) {
              recoQAHist.fill(HIST("hRecoCandidateCount"), 3.5);
              recoQAHist.fill(HIST("hRecoDaugPVsITSNSigma"), daughterTrack.sign() * daughterTrack.p(), itsNSigma);
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HyperkinkQa, processMC, "do QA for MC prodcutions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{
    // Parse the metadata
    adaptAnalysisTask<HyperkinkRecoTask>(cfgc),
    adaptAnalysisTask<HyperkinkQa>(cfgc),
  };
}
